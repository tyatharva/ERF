"""Verification package for the nested ChannelIslands run.

Usage: python3 nestverify.py <run_dir> <era5_precip_tag> <d0> <d1>
e.g.   python3 nestverify.py /app/ERF/run_nested_jan jan 20230109 20230110

Reports, all on the NEST (level 1):
 1. per-shell w>1 fraction / w_rms at 6/12/24 h (contamination-growth check)
 2. 24-h rain_accum stats (nest interior d>=12) vs ERA5 tp over the same
    lat/lon footprint
 3. mass drift from the run log
 4. PNG maps: nest 24-h precip, k=0 temp, k=3 wind speed; ERA5 tp
"""
import sys, glob, re
import numpy as np
import yt
yt.set_log_level(50)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer

run, tag, d0, d1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

# --- run-frame LCC (erftools convention, area 36.0,-123.25,31.25,-115.25) ---
AREA = [36.0, -123.25, 31.25, -115.25]
lat1, lat2, lon1, lon2 = AREA[2], AREA[0], AREA[1], AREA[3]
delta = lat2 - lat1
lcc = (f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
       f"+lat_0={(lat1+lat2)/2:.6f} +lon_0={(lon1+lon2)/2:.6f} "
       f"+datum=WGS84 +units=m +no_defs")
to_ll = Transformer.from_crs(CRS.from_proj4(lcc), CRS.from_epsg(4326), always_xy=True)

PLO = (-219131.04, -144000.0)
DXF = 1000.0   # fine

def mosaic(ds, field, lev=1):
    gs = [g for g in ds.index.grids if g.Level == lev]
    lo = np.min([g.get_global_startindex() for g in gs], axis=0)
    hi = np.max([g.get_global_startindex() + g.ActiveDimensions for g in gs], axis=0)
    a = np.full(tuple(hi - lo), np.nan)
    for g in gs:
        s = g.get_global_startindex() - lo
        d = g.ActiveDimensions
        a[s[0]:s[0]+d[0], s[1]:s[1]+d[1], s[2]:s[2]+d[2]] = np.asarray(g['boxlib', field])
    return a, lo

pls = sorted([p for p in glob.glob(f'{run}/plt0*') if '.old' not in p],
             key=lambda q: int(q.split('plt')[-1]))
print(f'{len(pls)} plotfiles; last {pls[-1]}')

# --- 1. shell tables over time -------------------------------------------
picks = [p for frac in (0.25, 0.5, 1.0) for p in [pls[min(len(pls)-1, int(frac*(len(pls)-1)))]]]
seen = []
for p in picks:
    if p in seen: continue
    seen.append(p)
    ds = yt.load(p)
    w, lo = mosaic(ds, 'z_velocity')
    nx, ny, nz = w.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    d3 = np.repeat(np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])[:, :, None], nz, axis=2)
    t_hr = float(ds.current_time) / 3600.0
    row = [(100*(w[d3 == k] > 1).mean()) for k in (0, 2, 5, 10)]
    interior = 100*(w[d3 >= 20] > 1).mean()
    wr = [np.sqrt((w[d3 == k]**2).mean()) for k in (0, 2, 5, 10)]
    print(f't={t_hr:5.1f} h  w>1%% d0/d2/d5/d10: {row[0]:5.2f} {row[1]:5.2f} {row[2]:5.2f} {row[3]:5.2f}'
          f'  interior(d>=20): {interior:5.2f}   w_rms d0/d10: {wr[0]:5.3f}/{wr[3]:5.3f}')

# --- 2. precipitation -----------------------------------------------------
ds = yt.load(pls[-1])
rain, lo = mosaic(ds, 'rain_accum')
r2 = rain[:, :, 0]
nx, ny = r2.shape
ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
d2 = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])
core = d2 >= 12
v = r2[core]
print(f'\nNEST 24-h rain_accum (interior d>=12, n={core.sum()}): '
      f'mean {v.mean():6.2f}  p95 {np.percentile(v,95):6.2f}  p99 {np.percentile(v,99):6.2f}  MAX {v.max():6.2f} mm')

# nest footprint in lat/lon
Xe = PLO[0] + (lo[0] + np.arange(nx) + 0.5) * DXF
Ye = PLO[1] + (lo[1] + np.arange(ny) + 0.5) * DXF
XX, YY = np.meshgrid(Xe, Ye, indexing='ij')
LON, LAT = to_ll.transform(XX, YY)
la0, la1 = LAT[core].min(), LAT[core].max()
lo0, lo1 = LON[core].min(), LON[core].max()
print(f'nest-core footprint: lat {la0:.3f}..{la1:.3f} lon {lo0:.3f}..{lo1:.3f}')

import pygrib
tp = {}
lats = lons = None
for g in pygrib.open(f'/app/ERF/precip_check/era5_precip_{tag}.grib'):
    if g.shortName != 'tp': continue
    tp[(g.validityDate, g.validityTime)] = np.array(g.values) * 1000.0
    if lats is None: lats, lons = g.latlons()
sel = [k for k in sorted(tp) if (k[0] == d0 and k[1] > 0) or (k[0] == d1 and k[1] == 0)]
tot = sum(tp[k] for k in sel)
lons = np.where(lons > 180, lons - 360, lons)
m = (lats >= la0) & (lats <= la1) & (lons >= lo0) & (lons <= lo1)
ev = tot[m]
print(f'ERA5 24-h tp over same footprint (n={m.sum()}): '
      f'mean {ev.mean():6.2f}  p95 {np.percentile(ev,95):6.2f}  MAX {ev.max():6.2f} mm')
print(f'RATIOS model/ERA5: mean {v.mean()/max(ev.mean(),1e-9):5.2f}x   max {v.max()/max(ev.max(),1e-9):5.2f}x')

# --- 3. mass drift --------------------------------------------------------
masses = []
for ln in open(f'{run}/run24.log'):
    mm = re.match(r'\s*MASS\s+SL/ML = (\S+) (\S+)', ln)
    if mm: masses.append((float(mm.group(1)), float(mm.group(2))))
if masses:
    print(f'\nMASS drift: SL {100*(masses[-1][0]/masses[0][0]-1):+7.3f}%   '
          f'ML {100*(masses[-1][1]/masses[0][1]-1):+7.3f}%  over {len(masses)} samples')

# --- 4. maps --------------------------------------------------------------
temp, _ = mosaic(ds, 'temp')
u, _ = mosaic(ds, 'x_velocity'); vv, _ = mosaic(ds, 'y_velocity')
spd = np.sqrt(u[:, :, 3]**2 + vv[:, :, 3]**2)
fig, axs = plt.subplots(2, 2, figsize=(15, 10))
for ax, fld, ttl, cm in ((axs[0,0], r2.T, 'nest 24-h precip (mm)', 'viridis'),
                          (axs[0,1], temp[:, :, 0].T - 273.15, 'nest T k=0 (C)', 'RdYlBu_r'),
                          (axs[1,0], spd.T, 'nest wind speed k=3 (m/s)', 'plasma')):
    im = ax.pcolormesh(Xe/1e3, Ye/1e3, fld, cmap=cm, shading='auto')
    plt.colorbar(im, ax=ax); ax.set_title(ttl); ax.set_aspect('equal')
im = axs[1,1].pcolormesh(lons, lats, tot, cmap='viridis', shading='auto')
axs[1,1].set_xlim(lo0-0.3, lo1+0.3); axs[1,1].set_ylim(la0-0.3, la1+0.3)
plt.colorbar(im, ax=axs[1,1]); axs[1,1].set_title('ERA5 24-h tp (mm)')
plt.tight_layout()
plt.savefig(f'{run}/verify_maps.png', dpi=110)
print(f'wrote {run}/verify_maps.png')
