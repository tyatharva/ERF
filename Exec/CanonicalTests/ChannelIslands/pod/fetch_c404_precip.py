#!/usr/bin/env python3
"""CONUS404 23-h precipitation on an ERF domain grid.

  fetch_c404_precip.py <x_lo> <y_lo> <x_hi> <y_hi> <nx> <ny> <out.npy>

Sums PREC_ACC_NC over 2020-12-28 01:00..23:00 (each file holds the preceding
hour) and area-averages it onto the ERF cell grid, matching how
wrf_d02_on_grid.npy was built -- bin-mean, which conserves the mean, so a
ratio of bin means is not biased by the 4 km -> 3 km change of support.

Needs C404_LAT / C404_LON / C404_BBOX, as conus404_to_bin.py does.

WHY THIS EXISTS: the .bin frames carry STATE only (TK/P/QVAPOR/...), never
precipitation, so a CONUS404-referenced score cannot be built from them.
"""
import os
import subprocess
import sys

import numpy as np
import pyproj

BASE = "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d559000"
OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")


def dods(url, shape, retries=3):
    enc = url.replace('[', '%5B').replace(']', '%5D')
    for a in range(retries):
        r = subprocess.run(['curl', '-g', '-sS', '--max-time', '900', enc],
                           capture_output=True)
        raw = r.stdout
        i = raw.find(b'\nData:\n')
        if i > 0:
            n = int(np.prod(shape))
            buf = raw[i + 7 + 8: i + 7 + 8 + 4 * n]
            if len(buf) == 4 * n:
                return np.frombuffer(buf, dtype='>f4').reshape(shape).astype(np.float64)
        print(f'    retry {a+1}', flush=True)
    raise RuntimeError(f'fetch failed: {url[:120]}')


def main():
    xl, yl, xh, yh = [float(v) for v in sys.argv[1:5]]
    nx, ny = int(sys.argv[5]), int(sys.argv[6])
    out = sys.argv[7]
    # Optional hour window [h0,h1] on 2020-12-28; default 0..23 keeps every
    # existing caller byte-identical. Files hold the PRECEDING hour, so the
    # window is files h0+1 .. h1 -- the same convention make_window_refs.py
    # uses for MRMS, so a 6-h driver total lines up with the 6-h references.
    H0 = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    H1 = int(sys.argv[9]) if len(sys.argv) > 9 else 23
    j0, j1, i0, i1 = [int(v) for v in np.load(os.environ['C404_BBOX'])]
    la = np.load(os.environ['C404_LAT'])[j0:j1+1, i0:i1+1]
    lo = np.load(os.environ['C404_LON'])[j0:j1+1, i0:i1+1]
    sub = f'[{j0}:1:{j1}][{i0}:1:{i1}]'
    shape = (1, j1 - j0 + 1, i1 - i0 + 1)

    tot = np.zeros(shape[1:])
    for h in range(H0 + 1, H1 + 1):
        stamp = f'2020-12-28_{h:02d}:00:00'.replace(':', '%3A')
        u = f'{BASE}/wy2021/202012/wrf2d_d01_{stamp}.nc.dods?PREC_ACC_NC[0:1:0]{sub}'
        tot += dods(u, shape)[0]
        print(f'  h{h:02d} cumulative mean {tot.mean():7.3f} mm', flush=True)

    # BILINEAR INTERPOLATION, not bin-averaging.
    #
    # wrf_d02_on_grid.npy is bin-averaged because d02 is 1.5 km -- FINER than
    # the 3 km target, so every target cell receives ~4 samples and the bin mean
    # both fills the grid and conserves the mean. CONUS404 is 4 km, COARSER than
    # the target: one source cell covers ~1.8 target cells, so binning leaves
    # 44% of the grid empty (measured). Coarse -> fine is an interpolation, and
    # this is the same bilinear map conus404_to_bin.py uses for the frames.
    C404 = ("+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.100006103515625 "
            "+lon_0=-97.9000015258789 +a=6370000 +b=6370000 +units=m +no_defs")
    to_c404 = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(C404), always_xy=True)
    dx = (xh - xl) / nx
    tx = xl + (np.arange(nx) + 0.5) * dx
    ty = yl + (np.arange(ny) + 0.5) * dx
    TX, TY = np.meshgrid(tx, ty, indexing='ij')
    tlon, tlat = pyproj.Transformer.from_crs(
        pyproj.CRS.from_proj4(OURS), 4326, always_xy=True).transform(TX, TY)
    cx, cy = to_c404.transform(tlon, tlat)
    sx, sy = to_c404.transform(lo[0, 0], la[0, 0])
    fi = (cx - sx) / 4000.0
    fj = (cy - sy) / 4000.0
    nys, nxs = tot.shape
    if fi.min() < -0.5 or fi.max() > nxs - 0.5 or fj.min() < -0.5 or fj.max() > nys - 0.5:
        raise SystemExit(f'target escapes the source subset: i {fi.min():.1f}..{fi.max():.1f} '
                         f'of {nxs}, j {fj.min():.1f}..{fj.max():.1f} of {nys}')
    from scipy.ndimage import map_coordinates
    grid = map_coordinates(tot, [fj, fi], order=1, mode='nearest')
    print(f'  interpolated: source mean {tot.mean():.3f} mm -> target mean '
          f'{grid.mean():.3f} mm (bilinear preserves this for a smooth field)')
    np.save(out, grid)
    print(f'wrote {out}: mean {np.nanmean(grid):.3f} mm, max {np.nanmax(grid):.3f} mm, '
          f'{int(np.isnan(grid).sum())} NaN of {grid.size}')


if __name__ == '__main__':
    main()
