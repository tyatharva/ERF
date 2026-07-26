"""Write the HindCast IC hydrostatic anchor: ERA5 sp, surface orography and t2m
on the ERF cell-centre grid.

WHY THIS EXISTS.  init_thermo_from_hindcast integrates the hydrostatic equation
upward from a specified surface pressure.  Neither of the two anchors already in
the code is usable:

  * erf_enforce_hse assumes p = p_0 = 101325 Pa at z = 0 in EVERY column.  That
    is a fixed standard-atmosphere sea-level pressure -- it cannot represent the
    980 hPa low that is the whole point of a storm hindcast.
  * the frame binaries' bottom level is a byte-identical copy of the level above
    it (erftools height bug), and the lowest GOOD frame level inherits erftools'
    ~300 m downward displacement, so deriving p there via the EOS gives the
    pressure of ~455 m at the 155 m label: ~35 hPa too low.

ERA5's sp has neither problem.  Its companion `geopotential` (single-level, the
time-invariant surface geopotential) gives the height that sp is valid at, so p
can be carried hydrostatically to ERF's 3-km terrain height -- ERA5's 0.25 deg
orography and ERF's 3-km terrain differ by several hundred metres over the
San Gabriels, which is ~70 hPa if ignored.

OUTPUT FORMAT (little-endian, matches read_sfc_anchor in
Source/Utils/ERF_WeatherDataInterpolation.cpp):

    int32   magic = 0x45524653  ('ERFS')
    int32   nx, ny                       (must equal amr.n_cell[0], [1])
    float64 sp   [ny][nx]   Pa   surface pressure at the ERA5 surface
    float64 zorog[ny][nx]   m    height of that surface (geopotential / g)
    float64 t2m  [ny][nx]   K    2-m air temperature

Index order is j*nx + i, i fastest.
"""
import numpy as np, pygrib, struct, os, sys
from pyproj import CRS, Transformer

# ERF grid -- must match geometry.prob_lo/prob_hi and amr.n_cell in inputs_hindcast
PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59); NX, NY = 128, 64
AREA = [36.0, -123.25, 31.25, -115.25]
lat1, lat2, lon1, lon2 = AREA[2], AREA[0], AREA[1], AREA[3]
delta = lat2 - lat1; lon0 = (lon1 + lon2) / 2; lat0 = (lat1 + lat2) / 2
LCC = (f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
       f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")
to_ll = Transformer.from_crs(CRS.from_proj4(LCC), CRS.from_epsg(4326), always_xy=True)

GRIB = os.environ.get('GRIB', '/app/ERF/precip_check/era5_sfc_jan.grib')
DATE = int(os.environ.get('DATE', '20230109'))
TIME = int(os.environ.get('TIME', '0'))
OUT  = os.environ.get('OUT', '/app/ERF/precip_check/sfc_anchor_jan09.bin')
G0   = 9.80665                                    # WMO standard gravity: ERA5's z is in m2/s2

# ---- read the three fields at the requested valid time --------------------------
want = {'sp': None, 'z': None, '2t': None}
alat = alon = None
g = pygrib.open(GRIB)
for m in g:
    if m.shortName not in want:
        continue
    if m.validityDate != DATE or m.validityTime != TIME:
        continue
    want[m.shortName] = np.array(m.values, dtype=float)
    if alat is None:
        alat, alon = m.latlons()
g.close()
missing = [k for k, v in want.items() if v is None]
if missing:
    sys.exit(f'ERROR: {missing} not found in {GRIB} at {DATE} {TIME:04d}Z')

# ERA5 subsets are regular in lat/lon; collapse the 2-D coordinate arrays to vectors
lat_v = alat[:, 0]
lon_v = np.where(alon[0, :] > 180, alon[0, :] - 360, alon[0, :])
if not (np.allclose(np.diff(alat[0, :]), 0) and np.allclose(np.diff(alon[:, 0]), 0)):
    sys.exit('ERROR: ERA5 grid is not a regular lat/lon mesh')
print(f'ERA5 grid {alat.shape}  lat {lat_v[0]:.2f}..{lat_v[-1]:.2f}  '
      f'lon {lon_v[0]:.2f}..{lon_v[-1]:.2f}')


def bilinear(field, qlon, qlat):
    """Bilinear sample of `field` (on lat_v x lon_v) at the query points, with
    edge clamping -- the ERF footprint sits strictly inside the ERA5 area, so
    clamping only ever engages on round-off."""
    # lat_v may be descending; work in index space that increases with the vector
    la = np.interp(qlat, lat_v[::-1], np.arange(len(lat_v))[::-1]) if lat_v[0] > lat_v[-1] \
        else np.interp(qlat, lat_v, np.arange(len(lat_v)))
    lo = np.interp(qlon, lon_v, np.arange(len(lon_v)))
    i0 = np.clip(np.floor(la).astype(int), 0, len(lat_v) - 2)
    j0 = np.clip(np.floor(lo).astype(int), 0, len(lon_v) - 2)
    fa = np.clip(la - i0, 0.0, 1.0)
    fo = np.clip(lo - j0, 0.0, 1.0)
    return ((1 - fa) * (1 - fo) * field[i0,     j0    ] +
            (1 - fa) * fo       * field[i0,     j0 + 1] +
            fa       * (1 - fo) * field[i0 + 1, j0    ] +
            fa       * fo       * field[i0 + 1, j0 + 1])


# ---- ERF cell centres in lon/lat -----------------------------------------------
dxe = (PHI[0] - PLO[0]) / NX; dye = (PHI[1] - PLO[1]) / NY
ei, ej = np.meshgrid(np.arange(NX), np.arange(NY), indexing='xy')   # shape (NY, NX)
elon, elat = to_ll.transform(PLO[0] + (ei + 0.5) * dxe, PLO[1] + (ej + 0.5) * dye)
print(f'ERF footprint  lat {elat.min():.3f}..{elat.max():.3f}  '
      f'lon {elon.min():.3f}..{elon.max():.3f}')

sp    = bilinear(want['sp'], elon, elat)
zorog = bilinear(want['z'],  elon, elat) / G0
t2m   = bilinear(want['2t'], elon, elat)

for nm, a, unit in [('sp', sp, 'Pa'), ('zorog', zorog, 'm'), ('t2m', t2m, 'K')]:
    print(f'  {nm:6s} min {a.min():10.2f}  mean {a.mean():10.2f}  max {a.max():10.2f}  {unit}')
if not (np.all(sp > 5.0e4) and np.all(sp < 1.1e5)):
    sys.exit('ERROR: surface pressure outside 500-1100 hPa')
if not (np.all(t2m > 230.0) and np.all(t2m < 330.0)):
    sys.exit('ERROR: t2m outside 230-330 K')

with open(OUT, 'wb') as f:
    f.write(struct.pack('<iii', 0x45524653, NX, NY))
    for a in (sp, zorog, t2m):
        f.write(np.ascontiguousarray(a, dtype='<f8').tobytes())
print(f'wrote {OUT}  ({os.path.getsize(OUT)} bytes, expect {12 + 3*NX*NY*8})')
