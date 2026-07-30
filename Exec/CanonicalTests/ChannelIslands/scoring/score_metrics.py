"""#25 metrics battery on the scoring_jan9 foundation fields.

FSS windows are ODD cell counts; nominal 3/9/15/30/60 km map to 1/3/5/11/21
cells = actual 3/9/15/33/63 km (footnoted). Self-FSS==1 is asserted as the
instrument control. Believable threshold: FSS_useful = 0.5 + f/2 (f = obs base
rate). All comparisons on valid-MRMS cells; d>=20 variant excludes the
forcing-dominated 10-cell band with margin (prior convention).
"""
import numpy as np
from scipy.ndimage import uniform_filter, label

import sys
OUT = (sys.argv[1] if len(sys.argv)>1 else '/app/ERF/scoring_jan9') .rstrip('/') + '/'
erf  = np.load(OUT+'erf_mm.npy'); mrms = np.load(OUT+'mrms_mm.npy')
era5 = np.load(OUT+'era5_mm.npy'); terr = np.load(OUT+'terrain.npy')
lat = np.load(OUT+'lat.npy'); lon = np.load(OUT+'lon.npy')
NX, NY = erf.shape
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
d = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])
valid = np.isfinite(mrms)

def fss(fcst, obs, thr, w, mask):
    Pf = uniform_filter((fcst >= thr).astype(float), size=w, mode='constant')
    Po = uniform_filter((obs  >= thr).astype(float), size=w, mode='constant')
    m = mask
    num = np.nanmean((Pf[m]-Po[m])**2)
    den = np.nanmean(Pf[m]**2) + np.nanmean(Po[m]**2)
    return 1.0 - num/den if den > 0 else np.nan

WINS = [(3,1),(9,3),(15,5),(30,11),(60,21)]
assert abs(fss(mrms, mrms, 1.0, 3, valid) - 1.0) < 1e-12, 'self-FSS control failed'

for mtag, mask in [('d>=0', valid), ('d>=20', valid & (d >= 20))]:
    print(f'\n=== FSS ({mtag}, n={mask.sum()}) ===')
    print('thr  km   ERF     ERA5    useful')
    for thr in (1.0, 5.0):
        f_obs = float((mrms[mask] >= thr).mean()); useful = 0.5 + f_obs/2
        for km, w in WINS:
            print(f'{thr:3.0f}  {km:3d}  {fss(erf,mrms,thr,w,mask):+.3f}  {fss(era5,mrms,thr,w,mask):+.3f}   {useful:.3f}')

def cat(fcst, obs, thr, mask):
    F = fcst[mask] >= thr; O = obs[mask] >= thr
    h = (F&O).sum(); fa = (F&~O).sum(); ms = (~F&O).sum()
    pod = h/(h+ms) if h+ms else np.nan
    far = fa/(h+fa) if h+fa else np.nan
    csi = h/(h+ms+fa) if h+ms+fa else np.nan
    return pod, far, csi

print('\n=== categorical (3 km cells) & bias ===')
for mtag, mask in [('d>=0', valid), ('d>=20', valid & (d >= 20))]:
    for thr in (1.0, 5.0):
        pe = cat(erf, mrms, thr, mask); pa = cat(era5, mrms, thr, mask)
        print(f'{mtag} {thr:.0f}mm  ERF POD/FAR/CSI {pe[0]:.3f}/{pe[1]:.3f}/{pe[2]:.3f}'
              f'   ERA5 {pa[0]:.3f}/{pa[1]:.3f}/{pa[2]:.3f}')
    print(f'{mtag} bias: ERF {np.nanmean(erf[mask])/np.nanmean(mrms[mask]):.2f}x '
          f'({np.nanmean(erf[mask]):.2f} vs {np.nanmean(mrms[mask]):.2f} mm), '
          f'ERA5 {np.nanmean(era5[mask])/np.nanmean(mrms[mask]):.2f}x; '
          f'p90 ERF/MRMS {np.nanpercentile(erf[mask],90):.1f}/{np.nanpercentile(mrms[mask],90):.1f}, '
          f'max {np.nanmax(erf[mask]):.1f}/{np.nanmax(mrms[mask]):.1f}')

print('\n=== power spectrum (radial, precip anomaly, interior d>=20) ===')
def spec(a, mask):
    x = np.where(mask, np.nan_to_num(a), 0.0)
    core = x[20:NX-20, 20:NY-20]
    core = core - core.mean()
    F = np.abs(np.fft.rfft2(core * np.hanning(core.shape[0])[:,None] * np.hanning(core.shape[1])[None,:]))**2
    kx = np.fft.fftfreq(core.shape[0], 3.0)[:,None]
    ky = np.fft.rfftfreq(core.shape[1], 3.0)[None,:]
    k = np.sqrt(kx**2 + ky**2)
    bins = np.linspace(0, k.max(), 40)
    idx = np.digitize(k.ravel(), bins)
    P = np.array([F.ravel()[idx==i].mean() if (idx==i).any() else np.nan for i in range(1,len(bins))])
    lam = 1.0/(0.5*(bins[1:]+bins[:-1]))
    return lam, P
lam, Pe = spec(erf, valid); _, Pm = spec(mrms, valid); _, Pa = spec(era5, valid)
band = (lam >= 8) & (lam <= 19)
print(f'spectral ratio ERF/MRMS 8-19 km: {np.nanmean(Pe[band]/Pm[band]):.3f}')
print(f'spectral ratio ERA5/MRMS 8-19 km: {np.nanmean(Pa[band]/Pm[band]):.3f}')
np.save(OUT+'spec_lambda.npy', lam); np.save(OUT+'spec_erf.npy', Pe)
np.save(OUT+'spec_mrms.npy', Pm); np.save(OUT+'spec_era5.npy', Pa)

print('\n=== islands (land components not touching domain edge land mass) ===')
land = terr > 1.0
lab, nl = label(land)
edge_labels = set(lab[0,:]) | set(lab[-1,:]) | set(lab[:,0]) | set(lab[:,-1])
for L in range(1, nl+1):
    m = lab == L
    if L in edge_labels or m.sum() < 2: continue
    print(f'island comp {L}: n={m.sum()} at lat {lat[m].mean():.2f} lon {lon[m].mean():.2f} '
          f'zmax {terr[m].max():.0f}m | ERF {np.nanmean(erf[m]):.1f} MRMS {np.nanmean(mrms[m]):.1f} '
          f'ERA5 {np.nanmean(era5[m]):.1f} mm')

print('\n=== north-south gradient (row means over valid, d>=20 interior) ===')
mask = valid & (d >= 20)
prof_e = [np.nanmean(np.where(mask[:,j], erf[:,j], np.nan)) for j in range(NY)]
prof_m = [np.nanmean(np.where(mask[:,j], mrms[:,j], np.nan)) for j in range(NY)]
prof_a = [np.nanmean(np.where(mask[:,j], era5[:,j], np.nan)) for j in range(NY)]
np.save(OUT+'ns_erf.npy', prof_e); np.save(OUT+'ns_mrms.npy', prof_m); np.save(OUT+'ns_era5.npy', prof_a)
third = (NY-40)//3
s = slice(20, 20+third); n = slice(NY-20-third, NY-20)
for nm, p in [('ERF',prof_e),('MRMS',prof_m),('ERA5',prof_a)]:
    p = np.array(p)
    print(f'{nm}: south-third {np.nanmean(p[s]):.2f} mm, north-third {np.nanmean(p[n]):.2f} mm, '
          f'N/S ratio {np.nanmean(p[n])/np.nanmean(p[s]):.2f}')

print('\n=== max location / band check ===')
i0, j0 = np.unravel_index(np.nanargmax(erf), erf.shape)
print(f'ERF max {np.nanmax(erf):.0f} mm at d={d[i0,j0]} terrain={terr[i0,j0]:.0f}m '
      f'lat={lat[i0,j0]:.2f} lon={lon[i0,j0]:.2f}; MRMS there {mrms[i0,j0]:.1f} mm')
for dd in (0,1,2,5,10,15,25):
    m = valid & (d == dd)
    print(f'd={dd:2d}: ERF mean {np.nanmean(erf[m]):6.1f} max {np.nanmax(erf[m]):7.1f} | MRMS mean {np.nanmean(mrms[m]):5.1f}')
print('METRICS OK')
