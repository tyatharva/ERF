#!/usr/bin/env python3
"""Build ERF's MOST roughness map from CONUS404 land use.

  make_roughness_map.py <x_lo> <y_lo> <x_hi> <y_hi> <nx> <ny> <out.txt>

Writes the file erf.most.roughness_file_name expects. Requires C404_LU
(LU_INDEX.npy from pod/fetch_c404_invariant.py) and C404_LAT / C404_LON.

WHY THIS EXISTS
  erf.most.z0 is a SINGLE constant for all land. Until the item-60 land-mask
  fix that constant was never even reached -- every column took the Charnock
  SEA branch, so terrain ran at z0 ~ 1e-4 m instead of 0.1 m. The mask fix
  recovered most of that, but 0.1 m is still 2-10x too smooth for the SoCal
  ranges: Noah-MP's own table gives 0.20 m for closed shrubland, 0.50-0.60 m
  for savanna and 0.80-1.10 m for forest.

  This needs NO land surface model. z0 is a geometric surface property read
  from a land-use category; NOAH-MP is a prognostic soil/vegetation model and
  is a completely separate (and much larger) question.

WHAT IT DOES
  LU_INDEX (MODIFIED_IGBP_MODIS_NOAH, 21 categories) -> Z0MVT from
  Submodules/Noah-MP/parameters/NoahmpTable.TBL, &noahmp_modis_parameters.
  The table is PARSED, not transcribed, so it cannot drift from the submodule.

  Resampling is NEAREST NEIGHBOUR because LU_INDEX is categorical -- averaging
  category 7 (open shrubland) and 13 (urban) into 10 (grassland) would be
  meaningless.

FILE FORMAT (read_custom_roughness, ERF_SurfaceLayer.cpp:1111)
  Whitespace-separated `x y z0` triples on NODES, not cell centres:
      x = ProbLo[0] + i*dx,  i = 0..nx      (nx+1 values)
      y = ProbLo[1] + j*dy,  j = 0..ny      (ny+1 values)
  ordered with x varying FASTEST (the reader indexes ii + jj*(nx+1)). ERF
  matches coordinates to a 1e-4 m tolerance and falls back to a brute-force
  search per cell if the order is wrong, so the ordering is worth getting
  right: it is the difference between O(N) and O(N^2) at startup.
"""
import os
import re
import sys

import numpy as np
import pyproj

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
C404 = ("+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.100006103515625 "
        "+lon_0=-97.9000015258789 +a=6370000 +b=6370000 +units=m +no_defs")

TBL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', '..', '..', '..',
                   'Submodules', 'Noah-MP', 'parameters', 'NoahmpTable.TBL')

# Categories whose Z0MVT is 0.00 in the table but which are LAND for MOST.
# A zero roughness is not merely inaccurate, it is undefined: MOST evaluates
# log(z/z0). Water (17) and lake (21) are excluded because the Charnock sea
# branch overwrites them -- their value here is never read.
LAND_FLOOR = {15: 0.001,   # snow and ice
              16: 0.010}   # barren / sparsely vegetated
WATER_CATS = (17, 21)

# What to write at water categories. NOT an ocean roughness, deliberately.
# Cells that ERF's land mask calls LAND are read from this file, and that mask
# comes from LANDMASK in the surface frames while these categories come from
# LU_INDEX -- two different CONUS404 fields that disagree on ~220 coastal cells
# of this domain. Writing an ocean value here would hand those cells z0 ~ 1e-4
# over land, which is precisely the defect this map exists to remove. True sea
# cells never read it: the Charnock branch overwrites them. So the safe value
# is the deck's own land constant.
WATER_FILL = 0.1


def modis_z0_table():
    """Parse Z0MVT out of the &noahmp_modis_parameters block."""
    txt = open(os.path.normpath(TBL)).read()
    i = txt.index('&noahmp_modis_parameters')
    blk = txt[i:]
    m = re.search(r'^\s*Z0MVT\s*=\s*([^\n!]+)', blk, re.M)
    if not m:
        sys.exit('FATAL: Z0MVT not found in the MODIS block of NoahmpTable.TBL')
    vals = [float(v) for v in m.group(1).replace(',', ' ').split()]
    # 1-indexed by category
    z0 = {k + 1: v for k, v in enumerate(vals)}
    for c, f in LAND_FLOOR.items():
        if z0.get(c, 0.0) <= 0.0:
            z0[c] = f
    for c in WATER_CATS:
        z0[c] = WATER_FILL      # see WATER_FILL: a LAND fallback, not an ocean value
    return z0


def main():
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    xl, yl, xh, yh = [float(v) for v in sys.argv[1:5]]
    nx, ny = int(sys.argv[5]), int(sys.argv[6])
    out = sys.argv[7]

    z0map = modis_z0_table()
    print('Z0MVT (m) by MODIS category, parsed from NoahmpTable.TBL:')
    print('  ' + '  '.join(f'{c}:{z0map[c]:.2f}' for c in sorted(z0map)))

    lu = np.load(os.environ['C404_LU'])
    la = np.load(os.environ['C404_LAT'])
    lo = np.load(os.environ['C404_LON'])

    dx, dy = (xh - xl) / nx, (yh - yl) / ny
    X, Y = np.meshgrid(xl + dx * np.arange(nx + 1),
                       yl + dy * np.arange(ny + 1), indexing='ij')

    to_wgs = pyproj.Transformer.from_crs(pyproj.CRS.from_proj4(OURS), 4326, always_xy=True)
    to_c404 = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(C404), always_xy=True)
    tlon, tlat = to_wgs.transform(X, Y)
    tx, ty = to_c404.transform(tlon, tlat)
    sx, sy = to_c404.transform(lo[0, 0], la[0, 0])

    # NEAREST neighbour -- LU_INDEX is categorical
    fi = np.rint((tx - sx) / 4000.0).astype(int)
    fj = np.rint((ty - sy) / 4000.0).astype(int)
    if fi.min() < 0 or fj.min() < 0 or fi.max() >= lu.shape[1] or fj.max() >= lu.shape[0]:
        sys.exit(f'FATAL: target grid escapes the CONUS404 domain: '
                 f'i {fi.min()}..{fi.max()} of {lu.shape[1]}, '
                 f'j {fj.min()}..{fj.max()} of {lu.shape[0]}')

    cat = lu[fj, fi].astype(int)
    z0 = np.vectorize(lambda c: z0map.get(c, 0.1))(cat)

    uniq, cnt = np.unique(cat, return_counts=True)
    print(f'\nland-use histogram over {cat.size} nodes:')
    for c, n in sorted(zip(uniq, cnt), key=lambda t: -t[1]):
        print(f'  cat {c:2d}  {n:6d} nodes ({100.0*n/cat.size:5.2f}%)  z0={z0map.get(c,0.1):.3f} m')
    land = ~np.isin(cat, WATER_CATS)
    print(f'\nland nodes {land.sum()} of {cat.size} ({100.0*land.mean():.1f}%)')
    if land.any():
        print(f'land z0: min {z0[land].min():.3f}  mean {z0[land].mean():.3f}  '
              f'max {z0[land].max():.3f} m   (the deck constant was 0.100)')
    else:
        # A wholly maritime domain is a legitimate configuration (Domain B is
        # one), not an error -- the map is then inert because the Charnock sea
        # branch overwrites every cell. Say so rather than reducing over an
        # empty array.
        print('no land nodes: this domain is entirely water, so the map is '
              'inert (Charnock overwrites every cell)')

    with open(out, 'w') as f:
        for j in range(ny + 1):            # x varies FASTEST
            for i in range(nx + 1):
                f.write(f'{X[i,j]:.6f} {Y[i,j]:.6f} {z0[i,j]:.6f}\n')
    print(f'\nwrote {out}  ({(nx+1)*(ny+1)} nodes)')
    print('Use with:  erf.most.roughness_file_name = "' + os.path.basename(out) + '"')


if __name__ == '__main__':
    main()
