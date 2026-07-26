"""Validate the verification metric before trusting it.

Pearson on 24-h precipitation is a bad metric here: MRMS p50 is 0.2 mm against the
model's 2.1 mm, so the coefficient is set by a few extremes and a small spatial
displacement is punished as hard as no skill at all. Three checks, all
post-processing:

  1. spatial lag correlation over +-60 km -- a peak away from (0,0) is a
     georeferencing error, not a physics failure
  2. FSS at 3/9/15/30/60 km and POD/FAR/CSI at 1 and 5 mm -- scale-aware skill
  3. MRMS reliability, using Stage IV as an independent second opinion: where two
     products disagree, neither is a trustworthy target
"""
import numpy as np, pygrib, yt, glob, os, netCDF4 as nc
from pyproj import CRS, Transformer
yt.set_log_level(50)

PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59); NX, NY = 128, 64
A = [36.0, -123.25, 31.25, -115.25]; la1, la2, lo1, lo2 = A[2], A[0], A[1], A[3]
dl = la2-la1
LCC = (f'+proj=lcc +lat_1={la1+dl/6:.6f} +lat_2={la2-dl/6:.6f} '
       f'+lat_0={(la1+la2)/2:.6f} +lon_0={(lo1+lo2)/2:.6f} +datum=WGS84 +units=m +no_defs')
to_xy = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_proj4(LCC), always_xy=True)
dx = (PHI[0]-PLO[0])/NX; dy = (PHI[1]-PLO[1])/NY


def to_erf(vals, lat2d, lon2d):
    x, y = to_xy.transform(lon2d, lat2d)
    i = np.floor((x-PLO[0])/dx).astype(int); j = np.floor((y-PLO[1])/dy).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY) & np.isfinite(vals)
    tot = np.bincount(j[ok]*NX+i[ok], weights=vals[ok], minlength=NX*NY)
    cnt = np.bincount(j[ok]*NX+i[ok], minlength=NX*NY)
    return (np.where(cnt > 0, tot/np.maximum(cnt, 1), np.nan).reshape(NY, NX).T,
            cnt.reshape(NY, NX).T)


g = pygrib.open('/app/ERF/mrms_pass2.grib2'); m = list(g)[0]
mrv = np.where(np.array(m.values, dtype=float) >= 9999, np.nan, np.array(m.values, dtype=float))
mlat = np.linspace(m['latitudeOfFirstGridPointInDegrees'], m['latitudeOfLastGridPointInDegrees'], m['Nj'])
mlon = np.linspace(m['longitudeOfFirstGridPointInDegrees'], m['longitudeOfLastGridPointInDegrees'], m['Ni'])-360.0
nmsg = len(list(pygrib.open('/app/ERF/mrms_pass2.grib2')))
g.close()
MLON, MLAT = np.meshgrid(mlon, mlat)
mrms, mcnt = to_erf(mrv, MLAT, MLON)

d4 = nc.Dataset('/app/ERF/ncep_stage4.nc')
s4v = np.array(d4.variables['observation'][:], dtype=float)
s4v = np.where((s4v < -1e3) | (s4v > 1e30), np.nan, s4v)*25.4
s4x = np.array(d4.variables['x'][:]); s4y = np.array(d4.variables['y'][:])
hrap = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-105 +x_0=0 +y_0=0 '
                      '+a=6371200 +b=6371200 +units=m +no_defs')
s4ll = Transformer.from_crs(hrap, CRS.from_epsg(4326), always_xy=True)
S4X, S4Y = np.meshgrid(s4x, s4y)
s4lon, s4lat = s4ll.transform(S4X, S4Y)
stage4, _ = to_erf(s4v, s4lat, s4lon)


def last(run):
    pl = sorted([p for p in glob.glob(f'/app/ERF/{run}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    for c in reversed(pl):
        ds = yt.load(c); gg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(gg[('boxlib', 'rain_accum')])[:, :, 0]
        if not np.isnan(ra).any():
            return ra, np.asarray(gg[('boxlib', 'z_phys')])[:, :, 0]


MOD, ter = last('bdyfix/ic_hyd')
BAS, _ = last('bdyfix/sst_ctl')
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
dring = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])
BASE = (ter > 30) & (dring >= 3) & np.isfinite(mrms)
print(f'MRMS grib messages: {nmsg}   scoring cells (land d>=3): {BASE.sum()}')
print(f'MRMS source cells per ERF cell: median {int(np.median(mcnt[mcnt>0]))}, '
      f'min {int(mcnt[BASE].min())}')

# ---- 1. spatial lag correlation -------------------------------------------------
print('\n1. SPATIAL LAG CORRELATION  (model shifted; peak away from 0,0 = georeferencing error)')
best = (None, -9)
grid = {}
for li in range(-20, 21, 2):
    for lj in range(-20, 21, 2):
        sm = np.roll(np.roll(MOD, li, axis=0), lj, axis=1)
        msk = BASE.copy()
        if li > 0: msk[:li, :] = False
        elif li < 0: msk[li:, :] = False
        if lj > 0: msk[:, :lj] = False
        elif lj < 0: msk[:, lj:] = False
        if msk.sum() < 100: continue
        c = np.corrcoef(sm[msk], mrms[msk])[0, 1]
        grid[(li, lj)] = c
        if c > best[1]: best = ((li, lj), c)
print('     lag_y (km) ->      -60    -36    -12    +12    +36    +60')
for li in (-20, -12, -4, 4, 12, 20):
    row = [grid.get((li, lj), np.nan) for lj in (-20, -12, -4, 4, 12, 20)]
    print(f'   lag_x {li*3:+4d} km   ' + ' '.join(f'{v:+6.3f}' for v in row))
print(f'   zero lag correlation    {grid.get((0,0), np.corrcoef(MOD[BASE],mrms[BASE])[0,1]):+.3f}')
print(f'   PEAK at lag ({best[0][0]*3:+d}, {best[0][1]*3:+d}) km  corr {best[1]:+.3f}')

# ---- 2. FSS and contingency -----------------------------------------------------
def boxmean(f, n):
    if n <= 1: return f.copy()
    k = np.ones((n, n))/(n*n)
    pad = n//2
    fp = np.pad(f, pad, mode='edge')
    out = np.zeros_like(f, dtype=float)
    cs = np.cumsum(np.cumsum(fp, axis=0), axis=1)
    cs = np.pad(cs, ((1, 0), (1, 0)))
    for a in range(f.shape[0]):
        for b in range(f.shape[1]):
            out[a, b] = (cs[a+n, b+n]-cs[a, b+n]-cs[a+n, b]+cs[a, b])/(n*n)
    return out


print('\n2. FRACTIONS SKILL SCORE   (FSS > 0.5 = useful at that scale)')
print('   thresh  scale   ic_hyd    baseline   |  FSS_uniform (no-skill reference)')
for thr in (1.0, 5.0):
    for n, km in [(1, 3), (3, 9), (5, 15), (10, 30), (20, 60)]:
        bo = (mrms >= thr).astype(float); bo[~np.isfinite(mrms)] = 0
        res = []
        for fld in (MOD, BAS):
            bf = (fld >= thr).astype(float)
            Pf, Po = boxmean(bf, n), boxmean(bo, n)
            num = np.mean((Pf[BASE]-Po[BASE])**2)
            den = np.mean(Pf[BASE]**2)+np.mean(Po[BASE]**2)
            res.append(1-num/den if den > 0 else np.nan)
        f0 = bo[BASE].mean()
        print(f'   {thr:4.0f}mm {km:4d}km   {res[0]:+.3f}    {res[1]:+.3f}    |  {f0:.3f} base rate')

print('\n   POD / FAR / CSI')
print('   thresh          ic_hyd POD  FAR   CSI   |  baseline POD  FAR   CSI')
for thr in (1.0, 5.0):
    out = []
    for fld in (MOD, BAS):
        f = fld[BASE] >= thr; o = mrms[BASE] >= thr
        h = (f & o).sum(); fa = (f & ~o).sum(); ms = (~f & o).sum()
        pod = h/max(h+ms, 1); far = fa/max(h+fa, 1); csi = h/max(h+fa+ms, 1)
        out += [pod, far, csi]
    print(f'   {thr:4.0f}mm          {out[0]:.3f}  {out[1]:.3f} {out[2]:.3f}   |'
          f'           {out[3]:.3f}  {out[4]:.3f} {out[5]:.3f}')

# ---- 3. MRMS reliability --------------------------------------------------------
print('\n3. MRMS RELIABILITY   Stage IV as an independent second opinion')
both = BASE & np.isfinite(stage4)
r = np.corrcoef(mrms[both], stage4[both])[0, 1]
print(f'   cells with both products: {both.sum()}   corr(MRMS, StageIV) = {r:+.3f}')
print(f'   means: MRMS {mrms[both].mean():.2f} mm   StageIV {stage4[both].mean():.2f} mm')
rel = np.abs(mrms-stage4)/np.maximum(0.5*(mrms+stage4), 0.5)
for lab, sub in [('all scoring cells', both),
                 ('products agree within 50%', both & (rel < 0.5)),
                 ('products agree within 25%', both & (rel < 0.25))]:
    if sub.sum() < 20: continue
    cm = np.corrcoef(MOD[sub], mrms[sub])[0, 1]
    cs = np.corrcoef(MOD[sub], stage4[sub])[0, 1]
    cb = np.corrcoef(BAS[sub], mrms[sub])[0, 1]
    print(f'   {lab:28s} n={sub.sum():4d}  corr(ic_hyd,MRMS) {cm:+.3f}  '
          f'corr(ic_hyd,StageIV) {cs:+.3f}  corr(base,MRMS) {cb:+.3f}')
print(f'\n   corr(MRMS, terrain) {np.corrcoef(mrms[BASE], ter[BASE])[0,1]:+.3f}   '
      f'corr(StageIV, terrain) {np.corrcoef(stage4[both], ter[both])[0,1]:+.3f}')
