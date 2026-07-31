#!/usr/bin/env python3
"""Score a 72-h arm over the h48-h71 window against d02 and CONUS404.

  score_71h.py <run_dir> [label]

Model time is 72 h: start 2020-12-26 00Z, stop 2020-12-29 00Z. The SCORED
window is h48-h71 = 2020-12-28 00Z-23Z, the identical 23 h that every earlier
arm ran as its entire life -- so domA_d02.npy / domA_c404.npy apply directly.

Window accumulation is rain_accum(h71) - rain_accum(h48), which is what the
23-h arms measured as rain_accum(end) - rain_accum(0).

Tests the predictions recorded in PREDICTIONS_71h.md before launch:
  P5  near-wall enhancement in dN 0-6  -> predicted 1.0 +/- 0.7
  P7  Santa Ynez                       -> predicted 0.7-1.2x
"""
import sys

import numpy as np
import yt

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
REFS = '/app/ERF/refs'
NX, NY, DXM = 192, 96, 3000.0
OCEAN_Z = 12.5
SY_LAT, SY_LON = 34.486, -119.802
H0, H1 = 48, 71


def grab(run, hour):
    for p in plotfiles_by_time(f'/app/ERF/{run}'):
        ds = yt.load(p)
        if abs(float(ds.current_time) / 3600.0 - hour) < 0.02:
            g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
            ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
            if is_poisoned(ra):
                raise SystemExit(f'{run} h{hour}: plotfile poisoned, refusing')
            z = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
            return ra, z, p.split('/')[-1]
    raise SystemExit(f'{run}: no plotfile at h{hour}')


def main():
    run = sys.argv[1].rstrip('/')
    label = sys.argv[2] if len(sys.argv) > 2 else run

    ra0, _, p0 = grab(run, H0)
    ra1, zs, p1 = grab(run, H1)
    acc = ra1 - ra0
    if acc.min() < -1e-6:
        print(f'  WARNING: {int((acc < -1e-6).sum())} cells with NEGATIVE window '
              f'accumulation (restart discontinuity?) -- clipped')
        acc = np.maximum(acc, 0.0)

    lat = np.load(f'{REFS}/domA_lat.npy')
    lon = np.load(f'{REFS}/domA_lon.npy')
    d02 = np.load(f'{REFS}/domA_d02.npy')
    c404 = np.load(f'{REFS}/domA_c404.npy')
    terr = zs - OCEAN_Z
    LAND = terr > 20.0
    ok = LAND & np.isfinite(d02) & (d02 > 0.1)

    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dN = NY - 1 - jj

    print(f'=== {label}: h{H0}-h{H1} (2020-12-28 00Z-23Z) ===')
    print(f'    from {p0} -> {p1}')
    print(f'    terrain {terr.min():.0f}..{terr.max():.0f} m, {LAND.sum()} land cells')

    print(f'\n--- A. DOMAIN MEANS (mm over the 23 h window) ---')
    fin = np.isfinite(d02)
    print(f'  ERF      {acc.mean():8.3f}   (land only {acc[LAND].mean():7.3f})')
    print(f'  d02      {np.nanmean(d02):8.3f}   ERF/d02 on d02 cells '
          f'{acc[fin].mean()/np.nanmean(d02):6.2f}x')
    print(f'  CONUS404 {np.nanmean(c404):8.3f}   ERF/c404 '
          f'{acc[np.isfinite(c404)].mean()/np.nanmean(c404):6.2f}x')

    print(f'\n--- B. TERRAIN-STRATIFIED ERF/d02 (never aggregated) ---')
    DB = [(0, 6), (7, 12), (13, 20), (21, 29), (30, 47), (48, 95)]
    print(f'  {"terrain":13s}' + ''.join(f'{f"dN {a}-{b}":>15}' for a, b in DB))
    for hlo, hhi in [(20, 100), (100, 300), (300, 600), (600, 1000), (1000, 3000)]:
        row = ''
        for lo, hi in DB:
            m = ok & (terr >= hlo) & (terr < hhi) & (dN >= lo) & (dN <= hi)
            row += (f'{acc[m].mean()/d02[m].mean():7.2f}x (n{m.sum():<4d})'
                    if m.sum() >= 12 else f'{"--":>15}')
        print(f'  {f"{hlo}-{hhi} m":13s}{row}')

    print(f'\n--- C. P5: near-wall enhancement, dN 0-6 vs far field ---')
    print(f'    PREDICTED 1.0 +/- 0.7 (decay fit asymptote; h48 >> tau = 10.6 h)')
    CTRL = ok & (terr >= 300) & (terr < 1000)
    vals = {}
    for lo, hi in [(0, 6), (7, 12), (13, 20), (21, 29), (30, 47)]:
        m = CTRL & (dN >= lo) & (dN <= hi)
        if m.sum() >= 12:
            vals[(lo, hi)] = acc[m].mean() / d02[m].mean()
            print(f'    dN {lo:>2}-{hi:<2}  n={m.sum():<5d} d02 {d02[m].mean():7.2f}  '
                  f'ERF {acc[m].mean():7.2f}  {vals[(lo,hi)]:6.2f}x')
    near = [v for k, v in vals.items() if k[1] <= 12]
    far = [v for k, v in vals.items() if k[0] >= 21]
    if near and far:
        r = np.mean(near) / np.mean(far)
        print(f'    near/far {r:.2f}   (23-h arms measured ERF 4.49, d02 1.83, MRMS 1.76)')
        print(f'    P5 {"HELD" if 0.3 <= np.mean(near) <= 1.7 else "FALSIFIED"}: '
              f'near-wall {np.mean(near):.2f}x vs predicted 1.0 +/- 0.7')

    print(f'\n--- D. P7: Santa Ynez ---')
    print(f'    PREDICTED 0.7-1.2x (rotation fix moved it 2.13 -> 0.83 at 23 h)')
    d = np.hypot(lat - SY_LAT, lon - SY_LON)
    si, sj = np.unravel_index(d.argmin(), d.shape)
    nb = (np.abs(ii - si) <= 2) & (np.abs(jj - sj) <= 2) & LAND
    r1 = acc[si, sj] / max(d02[si, sj], 1e-9)
    rn = acc[nb].mean() / d02[nb].mean()
    print(f'    cell ({si},{sj}) lat {lat[si,sj]:.3f} lon {lon[si,sj]:.3f} '
          f'terrain {terr[si,sj]:.0f} m  dN {dN[si,sj]}')
    print(f'    single cell : ERF {acc[si,sj]:7.2f}  d02 {d02[si,sj]:7.2f}  {r1:5.2f}x')
    print(f'    5x5 (n={nb.sum()}) : ERF {acc[nb].mean():7.2f}  d02 {d02[nb].mean():7.2f}  '
          f'{rn:5.2f}x')
    print(f'    P7 {"HELD" if 0.7 <= rn <= 1.2 else "FALSIFIED"} on the 5x5')


if __name__ == '__main__':
    main()
