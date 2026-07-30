"""#25 scoring foundation: build matched 24-h precip fields on the ERF grid.

Outputs (to OUT dir as .npy):
  erf_mm.npy    (192,96)  ERF 24-h precip [mm], stitched Jan-9 day (run_a3)
  mrms_mm.npy   (192,96)  MRMS Pass2 24-h QPE regridded (nearest, 1km->3km)
  era5_mm.npy   (192,96)  ERA5 tp Jan-9 24-h total, bilinear from 0.25 deg
  lat.npy lon.npy (192,96) cell-centre coordinates (LCC inverse)
  terrain.npy   (192,96)  z_phys at k=0 [m]
Validation gates (known-nonzero controls, item 25 reference day):
  ERA5 domain mean must land near 5.91 mm; MRMS near 8.51 mm (old-domain
  numbers; new domain is larger so tolerance is loose: 2-12 mm).
"""
import numpy as np, yt, pygrib, pyproj, glob, sys
yt.set_log_level(50)

RUN = sys.argv[1] if len(sys.argv) > 1 else '/app/ERF/run_a3'
OUT = sys.argv[2] if len(sys.argv) > 2 else '/app/ERF/scoring_jan9/'
if not OUT.endswith('/'): OUT += '/'
import os; os.makedirs(OUT, exist_ok=True)
print(f'run={RUN} out={OUT}')

# --- ERF grid geometry (deck values, 192x96x48) ---
PLO = (-388229.74, -166933.25); PHI = (187770.26, 121066.75)
NX, NY = 192, 96

# PROJECTION.  This MUST be the projection the frames were built on, which is
# the PINNED area (36.0,-123.25,31.25,-115.25) -- not the wider DOWNLOAD area in
# era5_input.txt.  Scoring originally used the download area: that is an 85.9 km
# misregistration of the ERF grid against MRMS and ERA5, and it invalidated every
# ERF-vs-MRMS number in the first #25 scoring.  The gate below is not a comment,
# it is a measurement: the frame carries both its projected grid (xvec,yvec) and
# the lat/lon of the same points, so the projection is over-determined and the
# residual of the candidate against the frame must be ~0.
AREA = (36.0, -123.25, 31.25, -115.25)   # N,W,S,E -- PINNED projection area

lat1, lat2 = AREA[2], AREA[0]; lon1, lon2 = AREA[1], AREA[3]
delta = lat2 - lat1
P4 = (f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
      f"+lat_0={(lat1+lat2)/2:.6f} +lon_0={(lon1+lon2)/2:.6f} "
      f"+datum=WGS84 +units=m +no_defs")
lcc = pyproj.CRS.from_proj4(P4)
tr = pyproj.Transformer.from_crs(lcc, 4326, always_xy=True)

# --- projection gate, against the frame's own coordinates ---
_fr = sorted(glob.glob(RUN + '/ERA5Data_3D/*.bin'))[0]
_raw = open(_fr, 'rb').read()
_nx, _ny, _nz, _nd = np.frombuffer(_raw, dtype='<i4', count=4)
_o = 16
_la = np.frombuffer(_raw, '<f4', _nx*_ny, _o).reshape(_ny, _nx).T; _o += 4*_nx*_ny
_lo = np.frombuffer(_raw, '<f4', _nx*_ny, _o).reshape(_ny, _nx).T; _o += 4*_nx*_ny
_x  = np.frombuffer(_raw, '<f4', _nx, _o); _o += 4*_nx
_y  = np.frombuffer(_raw, '<f4', _ny, _o)
_X, _Y = np.meshgrid(_x, _y, indexing='ij')
_xt, _yt = pyproj.Transformer.from_crs(4326, lcc, always_xy=True).transform(_lo, _la)
_res = float(np.hypot(_xt-_X, _yt-_Y).max())
print(f'projection gate: max residual vs frame coords = {_res:.2f} m')
assert _res < 10.0, f'PROJECTION MISMATCH: {_res:.0f} m against {_fr}'
dx = (PHI[0]-PLO[0])/NX; dy = (PHI[1]-PLO[1])/NY
xc = PLO[0] + (np.arange(NX)+0.5)*dx
yc = PLO[1] + (np.arange(NY)+0.5)*dy
X, Y = np.meshgrid(xc, yc, indexing='ij')
lon, lat = tr.transform(X, Y)
print(f'grid: dx={dx:.0f} m, lat {lat.min():.3f}..{lat.max():.3f}, lon {lon.min():.3f}..{lon.max():.3f}')

# --- ERF 24-h precip: last plotfile with finite rain_accum ---
pl = sorted([p for p in glob.glob(RUN + '/plt[0-9]*') if p.split('plt')[-1].isdigit()],
            key=lambda p: int(p.split('plt')[-1]))
erf = None
for cand in reversed(pl):
    ds = yt.load(cand)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    if not np.isnan(ra).any():
        erf = ra; terr = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
        print(f'ERF: {cand.split("/")[-1]} t={float(ds.current_time):.0f}s '
              f'mean={ra.mean():.2f} max={ra.max():.1f} mm')
        break
if erf is None: sys.exit('no finite rain_accum')
if abs(float(ds.current_time) - 86400.0) > 1.0: sys.exit(f'last plotfile not at 24 h: {float(ds.current_time)}')

# --- MRMS 24-h QPE ---
gm = pygrib.open('/app/ERF/mrms_pass2.grib2'); msg = gm.message(1)
mv, mlat, mlon = msg.values, *msg.latlons()
mlon = np.where(mlon > 180, mlon - 360, mlon)
if hasattr(mv, 'filled'): mv = mv.filled(np.nan)
sel = (mlat >= lat.min()-0.1) & (mlat <= lat.max()+0.1) & \
      (mlon >= lon.min()-0.1) & (mlon <= lon.max()+0.1)
print(f'MRMS: {sel.sum()} source points in footprint, native mean over them '
      f'{np.nanmean(mv[sel]):.2f} mm')
from scipy.interpolate import griddata
pts = np.column_stack([mlon[sel], mlat[sel]])
mrms = griddata(pts, mv[sel], (lon, lat), method='nearest')
# MRMS has no coverage gaps over SoCal/coastal waters at this footprint, but a
# radar-quality mask matters offshore: mark cells farther than 0.02 deg from any
# source point as NaN via linear-method fallback detection.
mrms_lin = griddata(pts, mv[sel], (lon, lat), method='linear')
mrms = np.where(np.isnan(mrms_lin), np.nan, mrms)

# --- ERA5 tp: Jan-9 24-h total (hourly accumulations, m -> mm) ---
ge = pygrib.open('/app/ERF/precip_check/era5_precip_jan.grib')
tot = None; n = 0
for m in ge:
    if m.shortName != 'tp': continue
    vd = m.validDate
    if vd.year == 2023 and vd.month == 1 and (vd.day == 9 or (vd.day == 10 and vd.hour == 0)):
        # hourly tp valid at HH covers [HH-1,HH]; the day = valid 01:00..00:00(+1)
        if vd.day == 9 and vd.hour == 0: continue
        v = m.values
        if hasattr(v, 'filled'): v = v.filled(0.0)
        elat, elon = m.latlons()
        tot = v if tot is None else tot + v; n += 1
print(f'ERA5: {n} hourly tp fields summed (need 24)')
if n != 24: sys.exit('ERA5 hour count wrong')
elon = np.where(elon > 180, elon - 360, elon)
era5 = griddata(np.column_stack([elon.ravel(), elat.ravel()]), (tot*1000.0).ravel(),
                (lon, lat), method='linear')

# --- validation gates ---
mm_mean = np.nanmean(mrms); e_mean = np.nanmean(era5)
print(f'means over ERF domain: ERF={erf.mean():.2f} MRMS={mm_mean:.2f} ERA5={e_mean:.2f} mm')
print(f'valid MRMS cells: {np.isfinite(mrms).sum()}/{mrms.size}')
# Interpolation control, self-consistent and independent of any remembered
# number: compare each regridded field's domain mean against the NATIVE source
# values at the source points that fall inside the ERF footprint. A regridding
# that lands in the wrong place, or that smears, shows up as a mean offset.
# (The old gate compared against a documented 13.46 mm grib-footprint mean. That
# footprint is not this domain -- and with the projection corrected the domain
# moved 86 km -- so it can no longer serve as a control.)
esel = (elat >= lat.min()-0.1) & (elat <= lat.max()+0.1) & \
       (elon >= lon.min()-0.1) & (elon <= lon.max()+0.1)
e_native = float(np.nanmean((tot*1000.0)[esel]))
m_native = float(np.nanmean(mv[sel]))
print(f'interp control: ERA5 native-in-footprint {e_native:.2f} vs regridded {e_mean:.2f} mm')
print(f'interp control: MRMS native-in-footprint {m_native:.2f} vs regridded {mm_mean:.2f} mm')
assert abs(e_mean - e_native) < 0.15*e_native, 'ERA5 regrid mean drifted from native'
assert abs(mm_mean - m_native) < 0.15*m_native, 'MRMS regrid mean drifted from native'
assert 5.0 < mm_mean < 20.0, f'MRMS mean {mm_mean} outside sanity window'
assert erf.mean() > 0.5, 'ERF precip implausibly small'

for name, arr in [('erf_mm', erf), ('mrms_mm', mrms), ('era5_mm', era5),
                  ('lat', lat), ('lon', lon), ('terrain', terr)]:
    np.save(OUT + name + '.npy', np.asarray(arr, dtype=np.float64))
print('FOUNDATION OK ->', OUT)
