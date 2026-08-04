#!/usr/bin/env python3
"""Compute the CONUS404 source-subset bbox for a target ERF domain.

  make_c404_bbox.py <x_lo> <y_lo> <x_hi> <y_hi> <out.npy>

Writes [j0, j1, i0, i1] for conus404_to_bin.py's C404_BBOX. The subset must
contain the whole 6 km frame grid; conus404_to_bin.build_mapping() hard-raises
if it does not, so this adds SLACK cells and then dry-runs that same check.

Reads the full CONUS404 XLAT/XLONG from C404_LAT / C404_LON.
"""
import os
import sys

import numpy as np
import pyproj

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
C404 = ("+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.100006103515625 "
        "+lon_0=-97.9000015258789 +a=6370000 +b=6370000 +units=m +no_defs")
SLACK = 12          # source cells of margin beyond the frame grid


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    xl, yl, xh, yh = [float(v) for v in sys.argv[1:5]]
    out = sys.argv[5]
    margin = float(os.environ.get('C404_MARGIN', '66000.'))
    # Must match conus404_to_bin.py's C404_DX, same default, same reason.
    # The bbox itself is spacing-independent (the frame EXTENT is set by margin,
    # not by dx) -- this is here so the dry-run below tests the grid that will
    # actually be built rather than a 6 km stand-in for it.
    dx = float(os.environ.get('C404_DX', '6000.'))
    snap = lambda v, up: (np.ceil(v / dx) if up else np.floor(v / dx)) * dx
    XS = np.arange(snap(xl - margin, False), snap(xh + margin, True) + 1., dx)
    YS = np.arange(snap(yl - margin, False), snap(yh + margin, True) + 1., dx)

    la = np.load(os.environ['C404_LAT'])
    lo = np.load(os.environ['C404_LON'])
    to_c404 = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(C404), always_xy=True)

    X, Y = np.meshgrid(XS, YS, indexing='ij')
    tlon, tlat = pyproj.Transformer.from_crs(
        pyproj.CRS.from_proj4(OURS), 4326, always_xy=True).transform(X, Y)
    tx, ty = to_c404.transform(tlon, tlat)
    sx, sy = to_c404.transform(lo, la)

    # source index of the target bounding box, via the source grid's own x/y
    ix = np.abs(sx[0, :] - tx.min()).argmin(), np.abs(sx[0, :] - tx.max()).argmin()
    iy = np.abs(sy[:, 0] - ty.min()).argmin(), np.abs(sy[:, 0] - ty.max()).argmin()
    i0 = max(min(ix) - SLACK, 0); i1 = min(max(ix) + SLACK, la.shape[1] - 1)
    j0 = max(min(iy) - SLACK, 0); j1 = min(max(iy) + SLACK, la.shape[0] - 1)

    # dry-run the exact check build_mapping() performs
    ssx, ssy = to_c404.transform(lo[j0, i0], la[j0, i0])
    fi = (tx - ssx) / 4000.0
    fj = (ty - ssy) / 4000.0
    nys, nxs = j1 - j0 + 1, i1 - i0 + 1
    ok = (fi.min() >= -0.5 and fi.max() <= nxs - 0.5
          and fj.min() >= -0.5 and fj.max() <= nys - 0.5)
    print(f'frame grid {len(XS)}x{len(YS)}  source subset j {j0}..{j1} ({nys}) '
          f'i {i0}..{i1} ({nxs})')
    print(f'  fi {fi.min():.1f}..{fi.max():.1f} of {nxs}   '
          f'fj {fj.min():.1f}..{fj.max():.1f} of {nys}   -> {"OK" if ok else "ESCAPES"}')
    if not ok:
        sys.exit('FATAL: target grid escapes the source subset; raise SLACK')
    np.save(out, np.array([j0, j1, i0, i1]))
    print(f'  wrote {out}')


if __name__ == '__main__':
    main()
