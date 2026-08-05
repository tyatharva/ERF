#!/usr/bin/env python3
"""d02 and MRMS accumulation references over an ARBITRARY hourly window.

  make_window_refs.py <h0> <h1> [<tag> <prob_lo_x> <prob_lo_y> <nx> <ny>]

With the optional grid arguments the references are built on THAT grid and named
refs/<tag>_{d02,mrms}_h<h0>_h<h1>.npy. Without them the Domain A grid is used and
the names are domA_*, exactly as before -- so existing references are unaffected.
Added because the ocean diagnostic domain is 96x96 at a different prob_lo, and
silently scoring it against Domain A references would compare different places.

Writes refs/domA_d02_h<h0>_h<h1>.npy and refs/domA_mrms_h<h0>_h<h1>.npy.

WHY. refs/domA_{d02,mrms}.npy are 23-h totals (00Z-23Z) and were the campaign's
only references. A shorter scored window needs its own: a 6-h model accumulation
cannot be compared against a 23-h total. Both raw sources support any window --
the wrfout has 24 hourly frames and MRMS one file per hour.

CONVENTIONS ARE COPIED, NOT REINVENTED, so these are directly comparable with
the existing 23-h references:
  * d02 total is RAINNC + BUCKET_MM*I_RAINNC, the bucket convention used by
    pod/make_domain_refs.py. (NOT RAINNC+RAINC -- that is a different quantity
    and would not line up with the 23-h reference.)
  * Binning is by PROJECTION into the pinned LCC and floor((x-prob_lo)/dx), as
    in make_domain_refs.py / mrms_to_domA.py. The target grid is Lambert, so
    binning on raw lat/lon would misplace every cell.
  * Sources are finer than the 3 km target, so bin-AVERAGE. Target cells with no
    source stay NaN: a missing observation is not a dry observation.
  * MRMS files are labelled with the hour they END and accumulate the preceding
    hour, so window [h0,h1] is the files h0+1 .. h1.
"""
import glob
import gzip
import os
import shutil
import sys
import tempfile

import numpy as np
import netCDF4 as nc
import pygrib
import pyproj

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
PROB_LO = (-389768.39, -89294.60)
NX, NY, DX = 192, 96, 3000.0
TAG = 'domA'
OUT = '/app/ERF/refs'
WRFOUT = '/app/ERF/wrfout_d02_2020-12-28_00_00_00'


def _bin(tr, lon, lat, vals):
    """Bin-average scattered source values onto the target grid."""
    x, y = tr.transform(lon, lat)
    i = np.floor((x - PROB_LO[0]) / DX).astype(int)
    j = np.floor((y - PROB_LO[1]) / DX).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY) & np.isfinite(vals)
    s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
    np.add.at(s, (i[ok], j[ok]), vals[ok])
    np.add.at(c, (i[ok], j[ok]), 1.0)
    return s, c


def d02_window(h0, h1, tr):
    f = nc.Dataset(WRFOUT)
    B = float(f.BUCKET_MM)
    R = lambda k: (np.asarray(f.variables['RAINNC'][k])
                   + B * np.asarray(f.variables['I_RAINNC'][k]))
    tot = R(h1) - R(h0)
    la = np.asarray(f.variables['XLAT'][0]); lo = np.asarray(f.variables['XLONG'][0])
    s, c = _bin(tr, lo.ravel(), la.ravel(), tot.ravel())
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def mrms_window(h0, h1, tr):
    tot = np.zeros((NX, NY)); cnt = np.zeros((NX, NY)); nhr = 0
    for h in range(h0 + 1, h1 + 1):
        fs = sorted(glob.glob(f'/app/ERF/mrms/*20201228-{h:02d}0000.grib2*'))
        if not fs:
            sys.exit(f'no MRMS file ending {h:02d}Z')
        with tempfile.NamedTemporaryFile(suffix='.grib2', delete=False) as tf:
            with gzip.open(fs[0], 'rb') as gz:
                shutil.copyfileobj(gz, tf)
            tmp = tf.name
        try:
            g = pygrib.open(tmp); m = g[1]
            vals = np.asarray(m.values, dtype=float)
            lat, lon = m.latlons(); g.close()
        finally:
            os.unlink(tmp)
        lon = np.where(lon > 180, lon - 360, lon)
        box = (np.isfinite(vals) & (vals >= 0.0) & (lat > 31.0) & (lat < 38.0)
               & (lon > -126.0) & (lon < -114.0))
        s, c = _bin(tr, lon[box], lat[box], vals[box])
        tot += s; cnt += c; nhr += 1
    print(f'  MRMS: {nhr} hourly files ({h0+1:02d}Z..{h1:02d}Z)')
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(cnt > 0, tot / np.maximum(cnt, 1) * nhr, np.nan)


def main():
    global PROB_LO, NX, NY, TAG
    if len(sys.argv) == 8:
        TAG = sys.argv[3]
        PROB_LO = (float(sys.argv[4]), float(sys.argv[5]))
        NX, NY = int(sys.argv[6]), int(sys.argv[7])
        print(f'grid: {TAG}  prob_lo={PROB_LO}  {NX}x{NY} at {DX:.0f} m')
    elif len(sys.argv) != 3:
        sys.exit(__doc__)
    h0, h1 = int(sys.argv[1]), int(sys.argv[2])

    tr = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(OURS), always_xy=True)
    os.makedirs(OUT, exist_ok=True)

    d = d02_window(h0, h1, tr)
    np.save(f'{OUT}/{TAG}_d02_h{h0}_h{h1}.npy', d)
    print(f'd02  {h0:02d}Z-{h1:02d}Z: coverage {100*np.isfinite(d).mean():5.1f}%  '
          f'mean {np.nanmean(d):6.3f} mm  max {np.nanmax(d):6.1f}')

    m = mrms_window(h0, h1, tr)
    np.save(f'{OUT}/{TAG}_mrms_h{h0}_h{h1}.npy', m)
    print(f'mrms {h0:02d}Z-{h1:02d}Z: coverage {100*np.isfinite(m).mean():5.1f}%  '
          f'mean {np.nanmean(m):6.3f} mm  max {np.nanmax(m):6.1f}')


if __name__ == '__main__':
    main()
