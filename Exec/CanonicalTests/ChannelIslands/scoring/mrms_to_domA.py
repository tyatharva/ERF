#!/usr/bin/env python3
"""Regrid MRMS 1-km hourly QPE onto the Domain A 3-km grid, 2020-12-28 00Z-23Z.

  mrms_to_domA.py <outdir>

REGRIDDING CONVENTION -- area-averaging (bin mean), NOT nearest neighbour.
MRMS is ~1 km and the target is 3 km, so each target cell receives ~9 source
cells and the mean over them is the correct conservative-ish estimate. This is
the safe direction: bin-averaging a FINER source onto a COARSER grid is well
posed. (Bin-averaging a coarser source onto a finer grid is what left 44% of
cells empty earlier in this campaign; that failure mode cannot occur here.)

Target cells receiving ZERO source cells are left NaN rather than filled, and
the count is reported. MRMS has no over-ocean coverage, so those NaNs are the
land/coastal mask and must be applied to every source before comparison.

The 23 files are hourly QPE valid at 01Z..23Z, each accumulating the PRECEDING
hour, so their sum is the 00Z-23Z total -- the same 23 h as ERF h48->h71.
"""
import glob
import gzip
import os
import shutil
import sys
import tempfile

import numpy as np
import pygrib
import pyproj

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
PROB_LO = (-389768.39, -89294.60)
NX, NY, DX = 192, 96, 3000.0


def main():
    out = sys.argv[1].rstrip('/') if len(sys.argv) > 1 else '/app/ERF/refs'
    files = sorted(glob.glob('/app/ERF/mrms/*.grib2.gz'))
    print(f'{len(files)} MRMS files')

    tr = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(OURS),
                                     always_xy=True)
    tot = np.zeros((NX, NY))
    cnt = np.zeros((NX, NY))
    nhr = 0

    for f in files:
        with tempfile.NamedTemporaryFile(suffix='.grib2', delete=False) as tf:
            with gzip.open(f, 'rb') as gz:
                shutil.copyfileobj(gz, tf)
            tmp = tf.name
        try:
            g = pygrib.open(tmp)
            m = g[1]
            vals = np.asarray(m.values, dtype=float)
            lat, lon = m.latlons()
            g.close()
        finally:
            os.unlink(tmp)

        lon = np.where(lon > 180, lon - 360, lon)
        # MRMS missing is -3 (or -999); anything negative is no-coverage
        good = np.isfinite(vals) & (vals >= 0.0)
        # crop to a generous lat/lon box first -- the full CONUS grid is huge
        box = good & (lat > 31.0) & (lat < 38.0) & (lon > -126.0) & (lon < -114.0)
        if not box.any():
            print(f'  {os.path.basename(f)}: no points in box'); continue
        x, y = tr.transform(lon[box], lat[box])
        v = vals[box]
        i = np.floor((x - PROB_LO[0]) / DX).astype(int)
        j = np.floor((y - PROB_LO[1]) / DX).astype(int)
        ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
        i, j, v = i[ok], j[ok], v[ok]
        np.add.at(tot, (i, j), v)
        np.add.at(cnt, (i, j), 1.0)
        nhr += 1
        print(f'  {os.path.basename(f)[-24:-11]}  {ok.sum():7d} src cells  '
              f'mean {v.mean():6.3f} mm')

    empty = int((cnt == 0).sum())
    # cnt accumulates across all hours; the per-hour mean is tot/cnt, and the
    # 23-h total is that mean times the number of hours.
    with np.errstate(invalid='ignore', divide='ignore'):
        mrms = np.where(cnt > 0, tot / np.maximum(cnt, 1) * nhr, np.nan)

    print(f'\n  hours used         : {nhr}')
    print(f'  target cells filled: {NX*NY - empty} / {NX*NY}  '
          f'({100.0*(NX*NY-empty)/(NX*NY):.1f}%)')
    print(f'  cells with NO MRMS : {empty}  (left NaN = the mask)')
    print(f'  mean over covered  : {np.nanmean(mrms):.3f} mm   '
          f'max {np.nanmax(mrms):.1f} mm')
    print(f'  median src cells per target per hour: '
          f'{np.median(cnt[cnt>0]/nhr):.1f}')

    np.save(f'{out}/domA_mrms.npy', mrms)
    print(f'\n  wrote {out}/domA_mrms.npy')


if __name__ == '__main__':
    main()
