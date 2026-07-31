#!/usr/bin/env python3
"""Per-domain coordinate, terrain and d02 reference arrays on an ERF grid.

  make_domain_refs.py <tag> <x_lo> <y_lo> <x_hi> <y_hi> <nx> <ny> <terrain_file>

Writes <tag>_{lat,lon,terrain,d02}.npy next to /app/ERF/refs/.

wrf_d02_on_grid.npy is on the PARENT 192x96 grid and is NOT reusable for a
shifted domain -- the whole point of the item-56 domains is that the box moved.
d02 is 1.5 km, FINER than the 3 km target, so it is bin-averaged (conserves the
mean); that is the same operation the checked-in parent reference used. Cells
with no d02 sample are left NaN, which is how "outside d02" must be represented
-- Domain B is only ~10% covered and silently zero-filling it would read as a
dry reference rather than as absent data.

Terrain is read from the ERF terrain file that the run actually used, not
re-derived from the DEM, so the figure shows the model's own orography.
"""
import os
import sys

import numpy as np
import netCDF4 as nc
import pyproj

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
OUT = '/app/ERF/refs'
WRFOUT = '/app/ERF/wrfout_d02_2020-12-28_00_00_00'


def main():
    tag = sys.argv[1]
    xl, yl, xh, yh = [float(v) for v in sys.argv[2:6]]
    nx, ny = int(sys.argv[6]), int(sys.argv[7])
    terrfile = sys.argv[8]
    os.makedirs(OUT, exist_ok=True)
    dx = (xh - xl) / nx

    inv = pyproj.Transformer.from_crs(pyproj.CRS.from_proj4(OURS), 4326, always_xy=True)
    xc = xl + (np.arange(nx) + 0.5) * dx
    yc = yl + (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xc, yc, indexing='ij')
    lon, lat = inv.transform(X, Y)
    np.save(f'{OUT}/{tag}_lat.npy', lat)
    np.save(f'{OUT}/{tag}_lon.npy', lon)

    # terrain: the ERF file is  nx \n ny \n <nx x-coords> <ny y-coords> <nx*ny z>,
    # i.e. 2 + nx + ny + nx*ny lines (dem_to_erf_terrain.py's own line count).
    # The nx+ny coordinate lines MUST be skipped -- reading straight through
    # returns the x coordinates as elevations, which shows up as a terrain
    # minimum equal to prob_lo x.
    with open(terrfile) as f:
        vals = [float(v) for v in f.read().split()]
    nxt, nyt = int(vals[0]), int(vals[1])
    expect = 2 + nxt + nyt + nxt * nyt
    if len(vals) != expect:
        raise SystemExit(f'{terrfile}: expected {expect} numbers, got {len(vals)}')
    z = np.array(vals[2 + nxt + nyt:]).reshape(nxt, nyt)
    ti = np.clip(np.round(np.arange(nx) * (nxt - 1) / (nx - 1)).astype(int), 0, nxt - 1)
    tj = np.clip(np.round(np.arange(ny) * (nyt - 1) / (ny - 1)).astype(int), 0, nyt - 1)
    terr = z[np.ix_(ti, tj)]
    np.save(f'{OUT}/{tag}_terrain.npy', terr)

    # d02, bin-averaged (it is finer than the target)
    f = nc.Dataset(WRFOUT)
    B = float(f.BUCKET_MM)
    R = lambda k: (np.asarray(f.variables['RAINNC'][k])
                   + B * np.asarray(f.variables['I_RAINNC'][k]))
    tot = R(23) - R(0)
    la = np.asarray(f.variables['XLAT'][0]); lo = np.asarray(f.variables['XLONG'][0])
    x, y = pyproj.Transformer.from_crs(
        4326, pyproj.CRS.from_proj4(OURS), always_xy=True).transform(lo, la)
    i = np.floor((x - xl) / dx).astype(int)
    j = np.floor((y - yl) / dx).astype(int)
    ok = (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
    s = np.zeros((nx, ny)); c = np.zeros((nx, ny))
    np.add.at(s, (i[ok], j[ok]), tot[ok])
    np.add.at(c, (i[ok], j[ok]), 1.0)
    d02 = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    np.save(f'{OUT}/{tag}_d02.npy', d02)

    cov = 100.0 * np.isfinite(d02).mean()
    print(f'{tag}: lat {lat.min():.2f}..{lat.max():.2f} lon {lon.min():.2f}..{lon.max():.2f}')
    print(f'  terrain {terr.min():.0f}..{terr.max():.0f} m   land(>20 m) {(terr>20).sum()} cells')
    print(f'  d02 coverage {cov:.2f}%   mean over covered {np.nanmean(d02):.3f} mm')


if __name__ == '__main__':
    main()
