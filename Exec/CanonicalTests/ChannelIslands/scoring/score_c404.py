#!/usr/bin/env python3
"""Score the CONUS404-driven arms at a fixed window against a chosen reference.

  score_c404.py <ref.npy> <reflabel> <arm_label>=<plotdir/plotfile> ...

The reference is wrfout d02 for the full-domain table and MRMS for the
observational table; both are pre-accumulated over the SAME window as the arms.
Reported: FSS ladder, percentile-matched FSS, PCC/bias/RMSE, interior mean,
spectrum, wall-band profile, fetch-binned bias, islands.

Scoring against d02 measures AGREEMENT WITH ANOTHER MODEL, not skill: d02 carries
a documented 1.09x wet bias concentrated on orographic peaks and +68% small-scale
variance against MRMS over land. MRMS-over-land is the observational check.
"""
import sys, numpy as np, yt
from scipy.ndimage import uniform_filter
yt.set_log_level(50)
NX, NY = 192, 96
WINS = [(3, 1), (9, 3), (15, 5), (30, 11), (60, 21)]
THRS = [1., 5., 15., 30.]

ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
D = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
FETCH = np.minimum(ii, jj)                 # distance from the xlo/ylo inflow walls
OUTF = np.minimum(NX - 1 - ii, NY - 1 - jj)
terr = np.load('/app/ERF/scoring_ab_dav/terrain.npy')
LAND = terr > terr.min() + 20
lat = np.load('/app/ERF/scoring_ab_dav/lat.npy')
lon = np.load('/app/ERF/scoring_ab_dav/lon.npy')


def fss(f, o, thr, w, m):
    Pf = uniform_filter((f >= thr).astype(float), size=w, mode='constant')
    Po = uniform_filter((o >= thr).astype(float), size=w, mode='constant')
    n = np.nanmean((Pf[m] - Po[m]) ** 2); d = np.nanmean(Pf[m] ** 2) + np.nanmean(Po[m] ** 2)
    return 1 - n / d if d > 0 else np.nan


def binz(f, m, r):
    return f >= float(np.nanquantile(f[m], 1 - r))


def spec(x, m):
    x = np.where(m, np.nan_to_num(x), 0.); c = x[20:NX-20, 20:NY-20]; c = c - c.mean()
    F = np.abs(np.fft.rfft2(c * np.hanning(c.shape[0])[:, None] * np.hanning(c.shape[1])[None, :])) ** 2
    kx = np.fft.fftfreq(c.shape[0], 3.)[:, None]; ky = np.fft.rfftfreq(c.shape[1], 3.)[None, :]
    k = np.sqrt(kx**2 + ky**2); b = np.linspace(0, k.max(), 40); idx = np.digitize(k.ravel(), b)
    P = np.array([F.ravel()[idx == i].mean() if (idx == i).any() else np.nan for i in range(1, len(b))])
    return 1. / (0.5 * (b[1:] + b[:-1])), P


def load_arm(spec_):
    lab, p = spec_.split('=', 1)
    ds = yt.load(p); g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    return lab, ra, float(ds.current_time)


def main():
    ref = np.load(sys.argv[1]); reflab = sys.argv[2]
    arms = [load_arm(s) for s in sys.argv[3:]]
    for lab, a, t in arms:
        print(f'  arm {lab:14s} t={t/3600:6.2f} h  mean {np.nanmean(a):7.3f} mm')
    print(f'  reference {reflab}: mean {np.nanmean(ref):.3f} mm\n')
    masks = [('full domain', np.isfinite(ref)),
             ('interior d>=20', (D >= 20) & np.isfinite(ref)),
             ('LAND', LAND & np.isfinite(ref))]
    for mname, m in masks:
        print(f'=== {mname} (n={int(m.sum())}) vs {reflab} ===')
        for lab, a, _ in arms:
            v = np.isfinite(a) & m
            x, y = a[v], ref[v]
            print(f'  {lab:14s} bias {x.mean()/y.mean():5.2f}x  PCC {np.corrcoef(x,y)[0,1]:6.3f}  '
                  f'RMSE {np.sqrt(((x-y)**2).mean()):6.2f}  mean {x.mean():6.2f} vs {y.mean():6.2f}')
        print(f'  {"FSS":>14} ' + ' '.join(f'{k[0]:>6}km' for k in WINS))
        for thr in THRS:
            base = float((ref[m] >= thr).mean())
            if base < 0.02 or base > 0.98:
                print(f'  {"thr %g"%thr:>14}  degenerate base rate {base:.3f} -- not scoreable')
                continue
            for lab, a, _ in arms:
                print(f'  {lab[:9]:>9} {thr:4.0f} ' +
                      ' '.join(f'{fss(a,ref,thr,w,m):8.3f}' for _, w in WINS) +
                      f'   base {base:.3f} useful {0.5+base/2:.3f}')
        print(f'  {"PM-FSS":>14} ' + ' '.join(f'{k[0]:>6}km' for k in WINS))
        for thr in (1., 15.):
            base = float((ref[m] >= thr).mean())
            if base < 0.02 or base > 0.98: continue
            ob = binz(ref, m, base)
            for lab, a, _ in arms:
                fb = binz(a, m, base)
                print(f'  {lab[:9]:>9} r={base:.2f} ' +
                      ' '.join(f'{fss(fb.astype(float),ob.astype(float),0.5,w,m):8.3f}' for _, w in WINS))
        print()
    lam, Pr = spec(ref, np.isfinite(ref)); b = (lam >= 8) & (lam <= 19)
    print('spectrum 8-19 km (ratio to reference):')
    for lab, a, _ in arms:
        _, Pa = spec(a, np.isfinite(ref))
        print(f'  {lab:14s} {np.nanmean(Pa[b]/Pr[b]):6.3f}')
    print('\nwall-band profile (mean mm by distance from ANY wall):')
    print(f'  {"d":>4} ' + ' '.join(f'{l[:9]:>10}' for l, _, _ in arms) + f' {reflab[:10]:>10}')
    for d in (0, 1, 2, 5, 10, 15, 20, 25):
        print(f'  {d:>4} ' + ' '.join(f'{np.nanmean(a[D==d]):10.2f}' for _, a, _ in arms) +
              f' {np.nanmean(ref[D==d]):10.2f}')
    print('\nfetch-binned ratio to reference (distance from INFLOW walls xlo/ylo):')
    for lo_, hi in ((0,0),(1,2),(6,9),(15,19),(20,29),(45,95)):
        mm = (FETCH >= lo_) & (FETCH <= hi) & np.isfinite(ref)
        print(f'  {lo_:>3}-{hi:<3} ' + ' '.join(
            f'{np.nanmean(a[mm])/np.nanmean(ref[mm]):8.2f}' for _, a, _ in arms))
    print('\nislands (three western):')
    from scipy.ndimage import label
    lab_, n = label(LAND, structure=np.ones((3,3)))
    edge = set(lab_[0,:]) | set(lab_[-1,:]) | set(lab_[:,0]) | set(lab_[:,-1])
    for i in range(1, n+1):
        if i in edge: continue
        s = lab_ == i
        if s.sum() < 4: continue
        print(f'  lat {lat[s].mean():6.2f} lon {lon[s].mean():8.2f} n={int(s.sum()):3d}  ' +
              ' '.join(f'{np.nanmean(a[s]):7.2f}' for _, a, _ in arms) +
              f' | {reflab[:8]} {np.nanmean(ref[s]):7.2f}')


if __name__ == '__main__':
    main()
