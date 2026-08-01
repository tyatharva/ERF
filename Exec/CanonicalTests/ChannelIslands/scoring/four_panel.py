#!/usr/bin/env python3
"""One 2x2 comparison: Davies / NSCBC over MRMS / d02, whole domain, 0-100 mm.

  four_panel.py <outpath.png>

Whole domain, every cell, common 0-100 mm scale so the four panels are directly
comparable by eye. Window h48-h71 = 2020-12-28 00Z-23Z.

Also prints a bias decomposition: how much of each arm's excess over d02 sits
in the four 10-cell lateral bands versus the interior.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
NX, NY, BAND = 192, 96, 10
VMAX = 100.0
H0, H1 = 48, 71


def window(run):
    acc, zs = {}, None
    for p in plotfiles_by_time(f'/app/ERF/{run}'):
        ds = yt.load(p)
        h = float(ds.current_time) / 3600.0
        for hr in (H0, H1):
            if abs(h - hr) < 0.02 and hr not in acc:
                g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
                ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
                if is_poisoned(ra):
                    raise SystemExit(f'{run} h{hr} poisoned')
                acc[hr] = ra
                if hr == H1:
                    zs = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
    return np.maximum(acc[H1] - acc[H0], 0.0), zs - 12.5


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else '/app/ERF/figs_compare/four_panel.png'
    dav, terr = window('run_A71_davies')
    nsc, _ = window('run_A71_nscbc')
    d02 = np.load('/app/ERF/refs/domA_d02.npy')
    mrms = np.load('/app/ERF/refs/domA_mrms.npy')

    LAND = terr > 20.0
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dmin = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    M = LAND & np.isfinite(d02) & np.isfinite(mrms)
    INT = M & (dmin >= BAND)
    BND = M & (dmin < BAND)

    panels = [(dav, 'ERF  NSCBC=off  (DAVIES control)'),
              (nsc, 'ERF  NSCBC  sigma=0.03'),
              (mrms, 'MRMS  (observations)'),
              (d02, 'WRF d02  (1.5 km reference)')]

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))
    for a, (arr, title) in zip(ax.flat, panels):
        p = a.imshow(arr.T, origin='lower', vmin=0, vmax=VMAX,
                     cmap='turbo', aspect='auto')
        a.set_title(f'{title}\nwhole domain {np.nanmean(arr):.1f} mm   |   '
                    f'land {np.nanmean(arr[M]):.1f} mm   |   '
                    f'max {np.nanmax(arr):.0f} mm', fontsize=10)
        plt.colorbar(p, ax=a, label='mm / 23 h', extend='max')
        a.contour(LAND.T.astype(float), levels=[0.5], colors='w',
                  linewidths=0.6, alpha=0.9)
        for e in (BAND, NX - 1 - BAND):
            a.axvline(e, color='k', ls=':', lw=1.1)
        for e in (BAND, NY - 1 - BAND):
            a.axhline(e, color='k', ls=':', lw=1.1)
        a.set_xlabel('i'); a.set_ylabel('j')
    fig.suptitle('h48-h71 = 2020-12-28 00Z-23Z   |   whole domain   |   turbo, '
                 'common 0-100 mm (saturating; see per-panel max)\n'
                 'white = coastline (20 m)   |   black dotted = the four '
                 '10-cell lateral bands', fontsize=9.5)
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig.savefig(outp, dpi=135)
    plt.close(fig)
    print(f'wrote {outp}')

    print(f'\n=== WHERE THE EXCESS LIVES (land cells, vs d02) ===')
    print(f'    band = within {BAND} cells of ANY lateral wall')
    print(f'    n: band {BND.sum()}, interior {INT.sum()}, total {M.sum()}')
    print(f'\n  {"arm":8s}{"region":10s}{"n":>6}{"ERF mm":>9}{"d02 mm":>9}'
          f'{"ratio":>8}{"excess mm":>11}{"% of excess":>13}')
    for lab, e in (('DAVIES', dav), ('NSCBC', nsc)):
        tot_ex = (e[M] - d02[M]).sum()
        for nm, q in (('band', BND), ('interior', INT), ('ALL', M)):
            ex = (e[q] - d02[q]).sum()
            print(f'  {lab:8s}{nm:10s}{q.sum():>6}{e[q].mean():>9.1f}'
                  f'{d02[q].mean():>9.1f}{e[q].mean()/d02[q].mean():>8.2f}'
                  f'{ex/M.sum():>11.1f}{100*ex/tot_ex:>12.1f}%')
        print()

    print('=== INTERIOR ONLY: is the residual bias orographic? ===')
    print(f'  {"terrain":14s}{"n":>6}{"d02":>8}{"MRMS":>8}{"DAVIES":>9}{"NSCBC":>9}'
          f'{"Dav/d02":>9}{"NSC/d02":>9}')
    for lo, hi in [(20, 100), (100, 300), (300, 600), (600, 1000), (1000, 3000)]:
        q = INT & (terr >= lo) & (terr < hi)
        if q.sum() < 12:
            continue
        print(f'  {f"{lo}-{hi} m":14s}{q.sum():>6}{d02[q].mean():>8.1f}'
              f'{mrms[q].mean():>8.1f}{dav[q].mean():>9.1f}{nsc[q].mean():>9.1f}'
              f'{dav[q].mean()/d02[q].mean():>9.2f}{nsc[q].mean()/d02[q].mean():>9.2f}')

    print('\n=== INTENSITY DISTRIBUTION, interior land cells (mm/23h) ===')
    print(f'  {"source":10s}' + ''.join(f'{f"p{p}":>8}' for p in (50, 75, 90, 99))
          + f'{"max":>9}{"frac>25mm":>11}')
    for nm, a in (('d02', d02), ('MRMS', mrms), ('DAVIES', dav), ('NSCBC', nsc)):
        v = a[INT]
        print(f'  {nm:10s}' + ''.join(f'{np.percentile(v,p):>8.1f}'
                                      for p in (50, 75, 90, 99))
              + f'{v.max():>9.1f}{100*(v>25).mean():>10.1f}%')


if __name__ == '__main__':
    main()
