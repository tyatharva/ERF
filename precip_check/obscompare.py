"""Three-way 24-h precipitation comparison against gauge/radar QPE.

ERF (NSCBC), ERF (Davies) and ERA5 are all scored against NCEP Stage IV and MRMS
Pass-2 on ONE common grid over ONE masked footprint.

Design decisions that matter:

* Common grid = the ERA5 0.25 deg boxes. It is the coarsest of the four, so
  aggregating onto it avoids crediting/penalising ERF for point maxima that ERA5
  cannot represent, and it is the only grid on which ERA5 needs no interpolation.
* Aggregation is area-weighted (conservative), not bilinear, which is the correct
  operator for an accumulated flux. Source cells are equal-area within each
  product (ERF 3 km LCC, MRMS 0.01 deg, Stage IV 4 km HRAP), so binning source
  cell CENTRES into target boxes and taking the unweighted mean is equivalent to
  area weighting up to the partial cells on box edges -- with ~70 ERF cells and
  ~600 MRMS cells per box that edge effect is sub-percent.
* Radar QPE is not trustworthy over water. The domain is ~75% ocean, and both
  products degrade offshore (no radar, no gauges); MRMS Pass 2 is gap-filled, so
  it reports numbers over the ocean that are not observations. Boxes are kept
  only where the ERF terrain says land dominates.
"""
import numpy as np, netCDF4 as nc, pygrib, yt, glob, os
from pyproj import CRS, Transformer
yt.set_log_level(50)

PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59); NX, NY = 128, 64
AREA = [36.0, -123.25, 31.25, -115.25]
lat1, lat2, lon1, lon2 = AREA[2], AREA[0], AREA[1], AREA[3]
delta = lat2 - lat1; lon0 = (lon1 + lon2) / 2; lat0 = (lat1 + lat2) / 2
LCC = (f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
       f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")
to_ll = Transformer.from_crs(CRS.from_proj4(LCC), CRS.from_epsg(4326), always_xy=True)
LAND_M = float(os.environ.get('LAND_M', '30'))     # ERF terrain > this = land (ocean is 12 m)
LANDFRAC = float(os.environ.get('LANDFRAC', '0.5'))


def erf_field(cfg):
    pl = sorted([p for p in glob.glob(f'/app/ERF/bdyfix/{cfg}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    for c in reversed(pl):
        ds = yt.load(c); g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not np.isnan(ra).any():
            return ra, np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
    raise RuntimeError(cfg)


nsc, ter = erf_field('jan_nsc')
dav, _   = erf_field('jan_ctl')
dxe = (PHI[0]-PLO[0])/NX; dye = (PHI[1]-PLO[1])/NY
ei, ej = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
elon, elat = to_ll.transform(PLO[0]+(ei+0.5)*dxe, PLO[1]+(ej+0.5)*dye)

# ---- ERA5 target grid + field -------------------------------------------------
g = pygrib.open('/app/ERF/precip_check/era5_precip_jan.grib'); tp = {}
for m in g:
    if m.shortName != 'tp': continue
    tp[(m.validityDate, m.validityTime)] = np.array(m.values)*1000.0
    alat, alon = m.latlons()
g.close()
sel = [k for k in sorted(tp) if (k[0] == 20230109 and k[1] > 0) or (k[0] == 20230110 and k[1] == 0)]
era5 = sum(tp[k] for k in sel)
alon180 = np.where(alon > 180, alon-360, alon)
res = 0.25
print(f'ERA5: {len(sel)} hourly fields, grid {era5.shape}')

# ---- Stage IV ------------------------------------------------------------------
d = nc.Dataset('/app/ERF/ncep_stage4.nc')
s4 = np.array(d.variables['observation'][:], dtype=float)
s4 = np.where((s4 < -1e3) | (s4 > 1e30), np.nan, s4) * 25.4        # inches -> mm
s4x = np.array(d.variables['x'][:]); s4y = np.array(d.variables['y'][:])
hrap = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-105 +x_0=0 +y_0=0 '
                      '+a=6371200 +b=6371200 +units=m +no_defs')
s4_to_ll = Transformer.from_crs(hrap, CRS.from_epsg(4326), always_xy=True)
S4X, S4Y = np.meshgrid(s4x, s4y)
s4lon, s4lat = s4_to_ll.transform(S4X, S4Y)

# ---- MRMS ----------------------------------------------------------------------
g = pygrib.open('/app/ERF/mrms_pass2.grib2'); mm = list(g)[0]
mr = np.array(mm.values, dtype=float)
mr = np.where(mr >= 9999, np.nan, mr)
mrlat = np.linspace(mm['latitudeOfFirstGridPointInDegrees'],
                    mm['latitudeOfLastGridPointInDegrees'], mm['Nj'])
mrlon = np.linspace(mm['longitudeOfFirstGridPointInDegrees'],
                    mm['longitudeOfLastGridPointInDegrees'], mm['Ni']) - 360.0
g.close()

# ---- aggregate everything onto the ERA5 boxes inside the ERF footprint ----------
rows = []
for a in range(era5.shape[0]):
    for b in range(era5.shape[1]):
        clat, clon = alat[a, b], alon180[a, b]
        lo_la, hi_la = clat-res/2, clat+res/2
        lo_lo, hi_lo = clon-res/2, clon+res/2
        me = (elat >= lo_la) & (elat < hi_la) & (elon >= lo_lo) & (elon < hi_lo)
        if me.sum() < 20:            # box not substantially inside the ERF domain
            continue
        lf = float((ter[me] > LAND_M).mean())
        ms = (s4lat >= lo_la) & (s4lat < hi_la) & (s4lon >= lo_lo) & (s4lon < hi_lo)
        mm_ = (np.abs(mrlat[:, None] - clat) < res/2) & (np.abs(mrlon[None, :] - clon) < res/2)
        s4v = np.nanmean(s4[ms]) if ms.sum() and np.isfinite(s4[ms]).any() else np.nan
        mrv = np.nanmean(mr[mm_]) if mm_.sum() and np.isfinite(mr[mm_]).any() else np.nan
        rows.append(dict(lat=clat, lon=clon, lf=lf, n_erf=int(me.sum()),
                         nsc=float(nsc[me].mean()), dav=float(dav[me].mean()),
                         era5=float(era5[a, b]), s4=s4v, mrms=mrv,
                         ter=float(ter[me].mean())))

R = rows
print(f'\nboxes inside the ERF footprint: {len(R)}')


def stats(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if not v.size: return '   --      --      --      --      --'
    return ('%7.1f %7.1f %7.1f %7.1f %7.1f'
            % (v.mean(), np.percentile(v, 50), np.percentile(v, 90),
               np.percentile(v, 99), v.max()))


for label, sub in [('ALL boxes', R),
                   (f'LAND boxes (landfrac >= {LANDFRAC:g}, valid QPE)',
                    [r for r in R if r['lf'] >= LANDFRAC and np.isfinite(r['s4']) and np.isfinite(r['mrms'])])]:
    print(f'\n=== {label}: n = {len(sub)} ===')
    print('  field        mean     p50     p90     p99     max')
    for k, nm in [('mrms', 'MRMS'), ('s4', 'StageIV'), ('era5', 'ERA5'),
                  ('nsc', 'ERF NSCBC'), ('dav', 'ERF Davies')]:
        print(f'  {nm:<11s}' + stats([r[k] for r in sub]))
    if sub:
        for ref in ('mrms', 's4'):
            rr = [r for r in sub if np.isfinite(r[ref]) and r[ref] > 0.5]
            if not rr: continue
            tot = sum(r[ref] for r in rr)
            print(f'  ratio vs {ref.upper():7s}' + ''.join(
                ' %s %.2fx' % (nm, sum(r[k] for r in rr)/tot)
                for k, nm in [('era5', 'ERA5'), ('nsc', 'NSCBC'), ('dav', 'Davies')]))

# terrain stratification of the ERA5 ratio on the ERF grid
print('\n=== ERF vs ERA5 by terrain (ERF grid, d >= 20) ===')
dring = np.minimum.reduce([ei, ej, NX-1-ei, NY-1-ej])
for lab, lo, hi in [('flat <100 m', 0, 100), ('100-400 m', 100, 400)]:
    m = (dring >= 20) & (ter >= lo) & (ter < hi)
    if not m.sum(): continue
    box = [r for r in R if lo <= r['ter'] < hi]
    e5 = np.mean([r['era5'] for r in box]) if box else np.nan
    print('  %-12s n=%4d  NSCBC %6.2f  Davies %6.2f  ERA5(boxes,n=%d) %5.2f  -> NSCBC/ERA5 %.2fx'
          % (lab, m.sum(), nsc[m].mean(), dav[m].mean(), len(box), e5,
             nsc[m].mean()/e5 if e5 and np.isfinite(e5) else float('nan')))
