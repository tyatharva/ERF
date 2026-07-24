#!/usr/bin/env python3
"""Batched, resumable CDS downloader for the ChannelIslands hindcast.

WHY: erftools downloads ONE CDS request PER TIMESTEP. Fine for a day
(9 frames); fatal for a year of 3-hourly frames (5,840 queued requests,
~days of queue wall time). This script downloads NDAY-sized batches (one
request per stream per chunk), then splits them into the exact
era5_3d_YYYYMMDD_HHMM.grib / era5_surf_YYYYMMDD_HHMM.grib files the
erftools pipeline expects -- whose downloader then skips every existing
file and goes straight to processing. No erftools changes needed.

MEASURED (2026-07-24): a full-month 3-hourly 12-var 37-level subset
request (110,112 fields) is REJECTED by CDS ("cost limits exceeded");
16 days (56,832 fields) is the validated default. On a fresh "too
large" rejection the chunk is halved automatically and retried.

Resumability (run it again after any failure; it continues):
  1. a chunk whose per-timestep files ALL exist is skipped outright;
  2. a batch .grib that exists and passes the message-count check is
     not re-downloaded;
  3. splitting writes per-timestep files only when complete-for-that-
     timestep, so a killed split just re-splits from the batch file.
A batch file that fails validation (short download) is deleted so the
next run re-fetches it.

The variable lists MIRROR the patched erftools request lists (12
pressure-level vars x 37 levels; 7 surface vars INCLUDING
forecast_albedo -- keep in sync with erftools_fal_patch.py).

Usage (inside the erf-hindcast container, from era5_run/):
  python3 era5_batch_download.py --start 2023-01-01 --end 2024-01-01
  # then run WriteICFromERA5Data.py exactly as in RUNBOOK section 2.
"""
import argparse
import os
import sys
from datetime import datetime, timedelta

import cdsapi
import pygrib

PL_VARS = ["geopotential", "relative_humidity", "specific_cloud_ice_water_content",
           "specific_cloud_liquid_water_content", "specific_humidity",
           "specific_rain_water_content", "specific_snow_water_content",
           "temperature", "u_component_of_wind", "v_component_of_wind",
           "vertical_velocity", "vorticity"]
PL_LEVELS = ["1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100", "125",
             "150", "175", "200", "225", "250", "300", "350", "400", "450", "500",
             "550", "600", "650", "700", "750", "775", "800", "825", "850", "875",
             "900", "925", "950", "975", "1000"]
SFC_VARS = ["sea_surface_temperature", "surface_sensible_heat_flux",
            "surface_latent_heat_flux", "friction_velocity",
            "land_sea_mask", "forecast_albedo", "surface_pressure"]

STREAMS = {
    "3d": {
        "dataset": "reanalysis-era5-pressure-levels",
        "fields_per_time": len(PL_VARS) * len(PL_LEVELS),
        "extra": {"variable": PL_VARS, "pressure_level": PL_LEVELS},
    },
    "surf": {
        "dataset": "reanalysis-era5-single-levels",
        "fields_per_time": len(SFC_VARS),
        "extra": {"variable": SFC_VARS},
    },
}


def daterange_chunks(start, end, chunk_days):
    d = start
    while d < end:
        e = min(d + timedelta(days=chunk_days), end)
        yield d, e
        d = e


def chunk_timesteps(cs, ce, interval_hours):
    ts = []
    t = cs
    while t < ce:
        ts.append(t)
        t += timedelta(hours=interval_hours)
    return ts


def per_timestep_name(stream, t):
    return f"era5_{stream}_{t.strftime('%Y%m%d')}_{t.strftime('%H%M')}.grib"


def batch_name(stream, cs, ce):
    return f"era5_{stream}_batch_{cs.strftime('%Y%m%d')}_{(ce - timedelta(days=1)).strftime('%Y%m%d')}.grib"


def validate_batch(path, expected_msgs):
    try:
        with pygrib.open(path) as g:
            n = g.messages
    except Exception as e:
        print(f"  batch {path}: unreadable ({e})")
        return False
    if n != expected_msgs:
        print(f"  batch {path}: {n} messages, expected {expected_msgs}")
        return False
    return True


def download_batch(client, stream, spec, cs, ce, tsteps, area, out):
    days = sorted({t.strftime("%d") for t in tsteps})
    months = sorted({t.strftime("%m") for t in tsteps})
    years = sorted({t.strftime("%Y") for t in tsteps})
    if len(months) > 1 or len(years) > 1:
        # CDS day/month lists are cross-products; a chunk must not span months.
        raise RuntimeError(f"chunk {cs}..{ce} spans months -- reduce --chunk-days")
    times = sorted({t.strftime("%H:%M") for t in tsteps})
    req = {"product_type": ["reanalysis"],
           "year": years[0], "month": months[0], "day": days, "time": times,
           "area": area, "data_format": "grib", "download_format": "unarchived"}
    req.update(spec["extra"])
    client.retrieve(spec["dataset"], req, out)


def split_batch(stream, path, tsteps):
    """Split a batch grib into the per-timestep files erftools expects.

    Bucket by VALIDITY time, not dataDate/dataTime: forecast-stream
    surface variables (fluxes, friction velocity) are stamped with the
    forecast INIT time (06/18 UTC) in dataDate/dataTime; only
    validityDate/validityTime give the instant the field applies to
    (measured: init-time bucketing left counts 4..16 across timesteps).
    """
    want = {(int(t.strftime("%Y%m%d")), t.hour * 100): t for t in tsteps}
    buckets = {k: bytearray() for k in want}
    counts = {k: 0 for k in want}
    with pygrib.open(path) as g:
        for m in g:
            key = (m.validityDate, m.validityTime)
            if key in buckets:
                buckets[key] += m.tostring()
                counts[key] += 1
    per_time = min(counts.values())
    if per_time != max(counts.values()):
        raise RuntimeError(f"{path}: uneven message counts across timesteps ({min(counts.values())}..{max(counts.values())})")
    for key, t in want.items():
        fn = per_timestep_name(stream, t)
        if os.path.exists(fn):
            continue
        tmp = fn + ".part"
        with open(tmp, "wb") as f:
            f.write(bytes(buckets[key]))
        os.replace(tmp, fn)
    return len(want), per_time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="YYYY-MM-DD (inclusive)")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD (exclusive; 00Z frame of this day IS included as the final BC frame)")
    ap.add_argument("--interval-hours", type=int, default=3)
    ap.add_argument("--chunk-days", type=int, default=16,
                    help="days per CDS request (16 measured OK; a full month is rejected)")
    ap.add_argument("--area", default="36.0,-123.25,31.25,-115.25")
    ap.add_argument("--clean-batch", action="store_true",
                    help="delete batch gribs after a successful split")
    args = ap.parse_args()

    start = datetime.strptime(args.start, "%Y-%m-%d")
    # erftools' frame list is inclusive of the end instant (the final BC frame):
    end = datetime.strptime(args.end, "%Y-%m-%d") + timedelta(hours=args.interval_hours)
    area = [float(v) for v in args.area.split(",")]
    client = None

    for stream, spec in STREAMS.items():
        # Never let a chunk span a month boundary (CDS day-lists are per-month)
        month_edges = []
        d = start.replace(day=1)
        while d < end:
            month_edges.append(max(d, start))
            d = (d + timedelta(days=32)).replace(day=1)
        month_edges.append(end)

        for mi in range(len(month_edges) - 1):
            for cs, ce in daterange_chunks(month_edges[mi], month_edges[mi + 1], args.chunk_days):
                tsteps = chunk_timesteps(cs, ce, args.interval_hours)
                if not tsteps:
                    continue
                missing = [t for t in tsteps if not os.path.exists(per_timestep_name(stream, t))]
                tag = f"[{stream} {cs.date()}..{(ce - timedelta(seconds=1)).date()}]"
                if not missing:
                    print(f"{tag} complete ({len(tsteps)} frames) -- skip")
                    continue

                bf = batch_name(stream, cs, ce)
                expected = len(tsteps) * spec["fields_per_time"]
                if os.path.exists(bf) and not validate_batch(bf, expected):
                    print(f"{tag} deleting invalid batch file for re-download")
                    os.remove(bf)
                if not os.path.exists(bf):
                    if client is None:
                        client = cdsapi.Client()
                    days = args.chunk_days
                    while True:
                        try:
                            print(f"{tag} requesting {expected} fields ...")
                            download_batch(client, stream, spec, cs, ce, tsteps, area, bf)
                            break
                        except Exception as e:
                            msg = str(e)
                            if ("too large" in msg or "cost limits" in msg) and days > 1:
                                # Shrink THIS chunk in place: split it into two runs
                                # by just failing with instructions (simplest safe path).
                                sys.exit(f"{tag} FATAL: request rejected as too large.\n"
                                         f"Re-run with --chunk-days {days // 2}.")
                            raise
                    if not validate_batch(bf, expected):
                        os.remove(bf)
                        sys.exit(f"{tag} FATAL: downloaded batch failed validation "
                                 f"(deleted); re-run to retry.")

                nfiles, per_time = split_batch(stream, bf, tsteps)
                print(f"{tag} split -> {nfiles} frames ({per_time} fields each)")
                if args.clean_batch:
                    os.remove(bf)

    print("All requested frames present.")


if __name__ == "__main__":
    main()
