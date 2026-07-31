#!/usr/bin/env python3
"""Terrain-stratified bias, Santa Ynez / domain max, and a windward-lee split.

  terrain_windward.py <arm_label>=<plotfile> ...

Covers the parts of the battery that score_c404.py does not: elevation-binned
bias, the two named point/max checks in inches, and a windward/lee partition of
the ranges. Also reports column snow+graupel, because the classic WSM6
cool-season failure is frozen hydrometeor falling out too fast on the windward
slope instead of advecting over the crest.

The flow direction is taken from the arm's OWN low-level wind (mass-weighted
below 1500 m over land), not assumed, so the windward/lee split is derived from
each run rather than imposed. Windward = terrain rising along the flow.
"""
import sys
import numpy as np
import yt
from plt_guard import reject_if_poisoned

yt.set_log_level(50)
NX, NY, DX = 192, 96, 3000.0
MM2IN = 1.0 / 25.4
SY_LAT, SY_LON = 34.486, -119.802

S = '/app/ERF/scoring_ab_dav/'
lat = np.load(S + 'lat.npy')
lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy')
LAND = terr > terr.min() + 20.0

d02 = np.load('/app/ERF/wrf_d02_on_grid.npy')
mrms = np.load('/app/ERF/mrms_20201228_on_grid.npy')

_d = np.hypot(lat - SY_LAT, lon - SY_LON)
SY = np.unravel_index(np.argmin(_d), _d.shape)

# terrain gradient in m per m, on the model grid
gx, gy = np.gradient(terr, DX, DX)

BINS = [('flat   <100 m', (terr < 100) & LAND),
        ('100-500 m', (terr >= 100) & (terr < 500) & LAND),
        ('>500 m', (terr >= 500) & LAND)]



def load(spec):
    lab, p = spec.split('=', 1)
    ds = yt.load(p)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    reject_if_poisoned(lab, p, ra)
    u = np.asarray(g[('boxlib', 'x_velocity')])
    v = np.asarray(g[('boxlib', 'y_velocity')])
    rho = np.asarray(g[('boxlib', 'density')])
    z = np.asarray(g[('boxlib', 'z_phys')])
    try:
        frozen = (np.asarray(g[('boxlib', 'qsnow')])
                  + np.asarray(g[('boxlib', 'qgraup')]))
        col = (frozen * rho).sum(axis=2) * (z[:, :, 1] - z[:, :, 0])[..., None].mean()
    except Exception:
        col = None
    # mass-weighted mean wind below 1500 m AGL, land columns only
    agl = z - z[:, :, :1]
    m = (agl < 1500.0) & LAND[:, :, None]
    w = np.where(m, rho, 0.0)
    ubar = (u * w).sum() / max(w.sum(), 1e-9)
    vbar = (v * w).sum() / max(w.sum(), 1e-9)
    return lab, ra, float(ds.current_time), ubar, vbar, col


def main():
    arms = [load(s) for s in sys.argv[1:]]

    print('=== arm summary (23-h accumulation) ===')
    for lab, a, t, ub, vb, _ in arms:
        sp = np.hypot(ub, vb)
        drc = (np.degrees(np.arctan2(-ub, -vb)) + 360) % 360
        print(f'  {lab:14s} t={t/3600:6.2f} h  domain-mean {np.nanmean(a)*MM2IN:5.3f} in'
              f'   low-level flow {sp:4.1f} m/s from {drc:5.1f} deg')

    print('\n=== Santa Ynez point and domain max (inches) ===')
    print(f'  {"":14s} {"SantaYnez":>10} {"domain max":>11}')
    for lab, a, _, _, _, _ in arms:
        print(f'  {lab:14s} {a[SY]*MM2IN:10.2f} {np.nanmax(a)*MM2IN:11.2f}')
    print(f'  {"d02 (ref)":14s} {d02[SY]*MM2IN:10.2f} {np.nanmax(d02)*MM2IN:11.2f}')
    print(f'  {"MRMS (obs)":14s} {mrms[SY]*MM2IN:10.2f} {np.nanmax(mrms)*MM2IN:11.2f}')

    for rname, ref in (('d02', d02), ('MRMS', mrms)):
        print(f'\n=== terrain-stratified bias vs {rname} (ratio of means) ===')
        print(f'  {"bin":14s} {"n":>6} ' + ' '.join(f'{l[:11]:>12}' for l, _, _, _, _, _ in arms)
              + f' {"ref mean in":>12}')
        for bname, m in BINS:
            mm = m & np.isfinite(ref)
            row = ' '.join(f'{np.nanmean(a[mm])/np.nanmean(ref[mm]):12.2f}'
                           for _, a, _, _, _, _ in arms)
            print(f'  {bname:14s} {int(mm.sum()):6d} {row} {np.nanmean(ref[mm])*MM2IN:12.3f}')

    print('\n=== windward / lee over the ranges (terrain > 100 m) ===')
    for lab, a, _, ub, vb, _ in arms:
        sp = max(np.hypot(ub, vb), 1e-9)
        uh, vh = ub / sp, vb / sp
        upslope = gx * uh + gy * vh          # >0 : terrain rises along the flow
        rng = (terr > 100) & LAND
        wind = rng & (upslope > 0.002)
        lee = rng & (upslope < -0.002)
        for rname, ref in (('d02', d02), ('MRMS', mrms)):
            bw = np.nanmean(a[wind]) / np.nanmean(ref[wind])
            bl = np.nanmean(a[lee]) / np.nanmean(ref[lee])
            print(f'  {lab:14s} vs {rname:5s} windward {bw:5.2f}x (n={int(wind.sum())})'
                  f'   lee {bl:5.2f}x (n={int(lee.sum())})   W/L ratio {bw/bl:5.2f}')
        print(f'  {"":14s} {"":8s} mean mm: windward {np.nanmean(a[wind]):6.2f}'
              f'  lee {np.nanmean(a[lee]):6.2f}')

    if any(c is not None for _, _, _, _, _, c in arms):
        print('\n=== column snow+graupel proxy (arb. units, windward vs lee) ===')
        for lab, _, _, ub, vb, col in arms:
            if col is None:
                continue
            sp = max(np.hypot(ub, vb), 1e-9)
            upslope = gx * (ub / sp) + gy * (vb / sp)
            rng = (terr > 100) & LAND
            w_, l_ = rng & (upslope > 0.002), rng & (upslope < -0.002)
            print(f'  {lab:14s} windward {np.nanmean(col[w_]):10.4f}'
                  f'   lee {np.nanmean(col[l_]):10.4f}'
                  f'   W/L {np.nanmean(col[w_])/max(np.nanmean(col[l_]),1e-12):6.2f}')


if __name__ == '__main__':
    main()
