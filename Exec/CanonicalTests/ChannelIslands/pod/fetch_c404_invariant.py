#!/usr/bin/env python3
"""Fetch the three CONUS404 invariant arrays conus404_to_bin.py needs.

  fetch_c404_invariant.py <outdir>

Writes XLAT.npy, XLONG.npy, LANDMASK.npy -- the full 1015x1367 source grid --
then prints the exports to feed the converter.

HANDOFF.md recorded these as the one campaign input that existed only on the
rented pod: "No committed script downloads XLAT/XLONG/LANDMASK." All three are
in a single file, INVARIANT/USGS404_geo_em_d01.nc, as XLAT_M / XLONG_M /
LANDMASK -- so this is one pass of ~17 MB, not the hour the handoff estimated.

LANDMASK is here and NOT in wrf2d/wrf3d, which is what item 55 was about.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from conus404_to_bin import BASE, dods           # noqa: E402  (path set above)

GEO = f'{BASE}/INVARIANT/USGS404_geo_em_d01.nc.dods'
NY, NX = 1015, 1367

# name in file -> output file, and the range it must land in for the fetch to
# be believed. A truncated .dods read returns short, not wrong, but a silently
# transposed or wrong-variable read returns the right SHAPE with absurd values.
# The latitude bound is deliberately loose: the grid is 1015x4 km = 4060 km
# tall about 39.1 N, and its Lambert corners reach 17.6 N -- well below the
# centre-edge latitude. A tighter bound rejects a correct fetch.
WANT = [('XLAT_M',   'XLAT.npy',     (10.0, 65.0)),
        ('XLONG_M',  'XLONG.npy',    (-145.0, -50.0)),
        ('LANDMASK', 'LANDMASK.npy', (0.0, 1.0)),
        # Dominant land-use category, MODIFIED_IGBP_MODIS_NOAH, 21 categories
        # (ISWATER=17, ISLAKE=21). Feeds make_roughness_map.py -- this is what
        # turns the uniform erf.most.z0 into a real roughness field, and it
        # needs no land surface model to use.
        ('LU_INDEX', 'LU_INDEX.npy', (1.0, 21.0))]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)

    for var, fn, (lo, hi) in WANT:
        print(f'fetching {var} ...', flush=True)
        a = dods(f'{GEO}?{var}[0:1:0][0:1:{NY-1}][0:1:{NX-1}]', (1, NY, NX))[0]
        if not (a.min() >= lo and a.max() <= hi):
            sys.exit(f'FATAL: {var} spans {a.min():.3f}..{a.max():.3f}, '
                     f'outside the expected {lo}..{hi}. Wrong variable or a bad read.')
        np.save(os.path.join(out, fn), a)
        print(f'  {fn}  {a.shape}  {a.min():.3f}..{a.max():.3f}', flush=True)

    lm = np.load(os.path.join(out, 'LANDMASK.npy'))
    print(f'\nland fraction over the full CONUS404 grid: {(lm >= 0.5).mean():.3f}')
    print(f"""
export C404_LAT={out}/XLAT.npy
export C404_LON={out}/XLONG.npy
export C404_LANDMASK={out}/LANDMASK.npy""")


if __name__ == '__main__':
    main()
