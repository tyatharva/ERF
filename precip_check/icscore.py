"""Score 24-h precipitation at land d>=3 against MRMS, plus the island points.

Deliberately re-scores the OLD runs with the SAME code as the new ones. The
historical reference numbers (Davies 4.05x bias / 0.351 correlation) came from a
different script; re-deriving them here both validates this implementation and
makes old-vs-new a like-for-like comparison even where it does not reproduce them
exactly.

Remap: MRMS is 0.01 deg, the ERF grid is 3 km, so MRMS is ~10x finer and binning
MRMS cell CENTRES into ERF cells and taking the unweighted mean is area-weighting
to within the partial cells on the edges (~100 MRMS cells per ERF cell).

Radar QPE is not trustworthy over water and MRMS Pass 2 is gap-filled, so it
reports numbers offshore that are not observations. Scoring is restricted to land.

env: RUNS (comma list of directories under /app/ERF), LAND_M
"""
import numpy as np, pygrib, yt, glob, os, sys
from pyproj import CRS, Transformer
yt.set_log_level(50)

PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59); NX, NY = 128, 64
AREA = [36.0, -123.25, 31.25, -115.25]
lat1, lat2, lon1, lon2 = AREA[2], AREA[0], AREA[1], AREA[3]
delta = lat2 - lat1; lon0 = (lon1 + lon2) / 2; lat0 = (lat1 + lat2) / 2
LCC = (f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
       f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")
to_xy = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_proj4(LCC), always_xy=True)

RUNS   = os.environ.get('RUNS', 'bdyfix/ic_ctl,bdyfix/ic_nsc,bdyfix/sst_ctl,bdyfix/sst_nsc').split(',')
LAND_M = float(os.environ.get('LAND_M', '30'))
DMIN   = int(os.environ.get('DMIN', '3'))
MRMS   = os.environ.get('MRMS', '/app/ERF/mrms_pass2.grib2')

ISLANDS = [('Catalina',    33.39, -118.42),
           ('San Clemente', 32.90, -118.49),
           ('San Nicolas',  33.24, -119.50)]

dx = (PHI[0]-PLO[0])/NX; dy = (PHI[1]-PLO[1])/NY
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
dring = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])

# ---- MRMS on the ERF grid --------------------------------------------------------
g = pygrib.open(MRMS); m = list(g)[0]
mr = np.array(m.values, dtype=float)
mr = np.where(mr >= 9999, np.nan, mr)
mrlat = np.linspace(m['latitudeOfFirstGridPointInDegrees'],
                    m['latitudeOfLastGridPointInDegrees'], m['Nj'])
mrlon = np.linspace(m['longitudeOfFirstGridPointInDegrees'],
                    m['longitudeOfLastGridPointInDegrees'], m['Ni']) - 360.0
g.close()
MLON, MLAT = np.meshgrid(mrlon, mrlat)
mx, my = to_xy.transform(MLON, MLAT)
mi = np.floor((mx - PLO[0])/dx).astype(int)
mj = np.floor((my - PLO[1])/dy).astype(int)
inside = (mi >= 0) & (mi < NX) & (mj >= 0) & (mj < NY) & np.isfinite(mr)
flat = mj[inside]*NX + mi[inside]
tot = np.bincount(flat, weights=mr[inside], minlength=NX*NY)
cnt = np.bincount(flat, minlength=NX*NY)
mrms = np.where(cnt > 0, tot/np.maximum(cnt, 1), np.nan).reshape(NY, NX).T
print(f'MRMS -> ERF grid: {int((cnt>0).sum())} of {NX*NY} cells covered, '
      f'median {int(np.median(cnt[cnt>0]))} MRMS cells per ERF cell')


def erf_field(run):
    pl = sorted([p for p in glob.glob(f'/app/ERF/{run}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    if not pl:
        return None, None, None
    # ERF writes a final plotfile after a vanishingly small clipped step whose
    # rain_accum is NaN while every prognostic field is clean; take the last
    # output with finite rain_accum.
    for c in reversed(pl):
        ds = yt.load(c)
        gg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(gg[('boxlib', 'rain_accum')])[:, :, 0]
        if not np.isnan(ra).any():
            return ra, np.asarray(gg[('boxlib', 'z_phys')])[:, :, 0], (c, float(ds.current_time))
    return None, None, None


rows = []
for run in RUNS:
    ra, ter, meta = erf_field(run)
    if ra is None:
        print(f'\n=== {run}: no plotfile with finite rain_accum ===')
        continue
    land = (ter > LAND_M)
    m3 = land & (dring >= DMIN) & np.isfinite(mrms)
    mo, ob = ra[m3], mrms[m3]
    bias = mo.sum()/ob.sum() if ob.sum() > 0 else np.nan
    corr = np.corrcoef(mo, ob)[0, 1]
    rmse = float(np.sqrt(((mo-ob)**2).mean()))
    print(f'\n=== {run}  ({meta[0].split("/")[-1]}, t = {meta[1]/3600:.2f} h) ===')
    print(f'  land d>={DMIN}, valid MRMS: n = {int(m3.sum())}')
    print(f'    bias  model/MRMS = {bias:6.2f}x   (model {mo.mean():7.2f} mm, MRMS {ob.mean():6.2f} mm)')
    print(f'    corr               {corr:+.3f}')
    print(f'    RMSE               {rmse:7.2f} mm')
    print(f'    model p50/p90/max  {np.percentile(mo,50):.1f} / {np.percentile(mo,90):.1f} / {mo.max():.1f} mm')
    print(f'    MRMS  p50/p90/max  {np.percentile(ob,50):.1f} / {np.percentile(ob,90):.1f} / {ob.max():.1f} mm')
    print('  islands (MRMS observes ~0 over these):')
    isl = []
    for nm, la, lo in ISLANDS:
        x, y = to_xy.transform(lo, la)
        i = int((x-PLO[0])/dx); j = int((y-PLO[1])/dy)
        if not (0 <= i < NX and 0 <= j < NY):
            print(f'    {nm:13s} OUTSIDE the domain'); continue
        ob_i = mrms[i, j]
        print(f'    {nm:13s} (i,j)=({i:3d},{j:2d}) ter {ter[i,j]:6.1f} m  d={dring[i,j]:2d}  '
              f'model {ra[i,j]:8.2f} mm   MRMS {ob_i if np.isfinite(ob_i) else float("nan"):6.2f} mm')
        isl.append(ra[i, j])
    print(f'    island mean model {np.mean(isl):.2f} mm')
    rows.append((run, bias, corr, rmse, float(np.mean(isl))))

if len(rows) > 1:
    print('\n=== summary (land d>=%d vs MRMS) ===' % DMIN)
    print('  run                    bias      corr     RMSE    island mean')
    for r in rows:
        print(f'  {r[0]:<22s} {r[1]:6.2f}x  {r[2]:+7.3f} {r[3]:8.2f} {r[4]:11.2f} mm')
    print('\n  reference from the BROKEN IC: Davies 4.05x / 0.351')
