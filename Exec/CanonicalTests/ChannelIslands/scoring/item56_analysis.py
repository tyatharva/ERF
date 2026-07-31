#!/usr/bin/env python3
"""Item-56 discriminator: did the excess follow the WALL or follow the TERRAIN?

  item56_analysis.py <tag>          e.g. domA

The parent domain could not answer this: 99.6% of its >500 m cells sat outside
the scored interior and 74.8% inside the d<10 band, so "near the north wall" and
"on the Transverse Range front" were the same cells. Domain A moves the pin
0.7 deg north, which moves Santa Ynez from dN 0-6 to dN ~21-29 and brings new,
lower terrain into the near-wall band.

  excess stays at dN 0-12  -> boundary artifact
  excess follows the range -> real orographic response

STRATIFIED BY TERRAIN HEIGHT, never aggregated. Domain A's dN 0-12 band tops out
far below the parent's 2368 m, so an aggregate near-wall number would confound a
weaker wall effect with lower terrain. Item 56 measured ~4x at EVERY height from
100-1000 m, so matched-height bins are the only valid comparison.
"""
import os
import sys

import numpy as np
import yt

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
REFS = '/app/ERF/refs'
NX, NY, DXM = 192, 96, 3000.0
OCEAN_Z = 12.5
SY_LAT, SY_LON = 34.486, -119.802
FLOW_DEG = 242.0

RUNS = {'domA': [('sigma=0.1', 'run_domA_sig01'), ('sigma=0.03', 'run_domA_sig003')]}


def load(run):
    for p in reversed(plotfiles_by_time(f'/app/ERF/{run}')):
        ds = yt.load(p)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not is_poisoned(ra):
            z = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
            return ra, z, os.path.basename(p)
    raise SystemExit(f'no clean plotfile in {run}')


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'domA'
    lat = np.load(f'{REFS}/{tag}_lat.npy'); lon = np.load(f'{REFS}/{tag}_lon.npy')
    d02 = np.load(f'{REFS}/{tag}_d02.npy')
    c404 = np.load(f'{REFS}/{tag}_c404.npy')
    arms = [(lab,) + load(run) for lab, run in RUNS[tag]]
    terr = arms[0][2] - OCEAN_Z
    LAND = terr > 20.0

    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dN, dE = NY - 1 - jj, NX - 1 - ii
    dW, dS = ii, jj
    gx, gy = np.gradient(terr, DXM, DXM)
    th = np.deg2rad(FLOW_DEG)
    ups = gx * (-np.sin(th)) + gy * (-np.cos(th))
    ERF = np.mean([a[1] for a in arms], axis=0)
    ok = LAND & np.isfinite(d02) & (d02 > 0.1)

    print(f'=== {tag}: terrain {terr.min():.0f}..{terr.max():.0f} m, '
          f'{LAND.sum()} land cells, >500 m {(LAND & (terr>=500)).sum()} ===')
    for lab, ra, _, pf in arms:
        print(f'  {lab:11s} {pf:12s} mean {ra.mean():7.3f} mm   vs d02 '
              f'{ra[np.isfinite(d02)].mean()/np.nanmean(d02):.2f}x')
    print(f'  d02 mean {np.nanmean(d02):.3f} mm   CONUS404 {np.nanmean(c404):.3f} mm')

    print('\n=== A. NEAR-WALL RATIO STRATIFIED BY TERRAIN HEIGHT (not aggregated) ===')
    print('    parent domain measured ~4x at EVERY height 100-1000 m in dN 0-12')
    HB = [(100, 300), (300, 600), (600, 1000), (1000, 3000)]
    DB = [(0, 6), (7, 12), (13, 20), (21, 29), (30, 47), (48, 95)]
    print(f'  {"terrain":13s}' + ''.join(f'{f"dN {a}-{b}":>16}' for a, b in DB))
    for hlo, hhi in HB:
        row = ''
        for lo, hi in DB:
            m = ok & (terr >= hlo) & (terr < hhi) & (dN >= lo) & (dN <= hi)
            row += f'{ERF[m].mean()/d02[m].mean():8.2f}x (n{m.sum():<4d})' \
                if m.sum() >= 12 else f'{"--":>16}'
        print(f'  {f"{hlo}-{hhi} m":13s}{row}')

    print('\n=== B. SANTA YNEZ -- the same mountain, moved ===')
    d = np.hypot(lat - SY_LAT, lon - SY_LON)
    si, sj = np.unravel_index(d.argmin(), d.shape)
    print(f'  cell ({si},{sj})  lat {lat[si,sj]:.3f} lon {lon[si,sj]:.3f}  '
          f'terrain {terr[si,sj]:.0f} m')
    print(f'  dN = {dN[si,sj]}   (parent domain: dN 0-6)')
    print(f'  {"":14s}{"ERF mm":>9}{"d02":>9}{"MRMS-":>9}{"ratio":>8}')
    for lab, ra, _, _ in arms:
        print(f'  {lab:14s}{ra[si,sj]:9.2f}{d02[si,sj]:9.2f}{"n/a":>9}'
              f'{ra[si,sj]/max(d02[si,sj],1e-9):7.2f}x')
    # a small neighbourhood, since a single cell is noisy
    nb = (np.abs(ii - si) <= 2) & (np.abs(jj - sj) <= 2) & LAND
    print(f'  5x5 neighbourhood (n={nb.sum()}, mean terrain {terr[nb].mean():.0f} m, '
          f'dN {dN[nb].min()}-{dN[nb].max()}):')
    for lab, ra, _, _ in arms:
        print(f'    {lab:12s} ERF {ra[nb].mean():7.2f} mm  d02 {d02[nb].mean():7.2f}  '
              f'ratio {ra[nb].mean()/d02[nb].mean():5.2f}x')

    print('\n=== C. PER-WALL ratio, terrain >= 100 m ===')
    for wn, dd in [('north yhi', dN), ('east xhi', dE), ('west xlo', dW), ('south ylo', dS)]:
        row = ''
        for lo, hi in [(0, 4), (5, 9), (10, 14), (15, 24), (25, 47), (48, 95)]:
            m = ok & (terr >= 100) & (dd >= lo) & (dd <= hi)
            row += f'{ERF[m].mean()/d02[m].mean():11.2f}x' if m.sum() >= 15 else f'{"--":>12}'
        print(f'  {wn:11s}{row}')
    print(f'  {"":11s}' + ''.join(f'{f"{a}-{b}":>12}' for a, b in
                                  [(0, 4), (5, 9), (10, 14), (15, 24), (25, 47), (48, 95)]))

    print('\n=== D. CONTROLLED near/far, matched terrain + upslope + inland ===')
    CTRL = ok & (terr >= 300) & (terr < 1000) & (ups > 0.01) & (ups < 0.06)
    print(f'  {"dN":>8}{"n":>6}{"terr":>7}{"ups":>7}{"d02":>8}{"ERF":>8}{"ratio":>8}')
    vals = {}
    for lo, hi in [(0, 6), (7, 12), (13, 20), (21, 29), (30, 47)]:
        m = CTRL & (dN >= lo) & (dN <= hi)
        if m.sum() < 12:
            print(f'  {f"{lo}-{hi}":>8}{m.sum():6d}   (too few)'); continue
        r = ERF[m].mean() / d02[m].mean(); vals[(lo, hi)] = r
        print(f'  {f"{lo}-{hi}":>8}{m.sum():6d}{terr[m].mean():7.0f}{ups[m].mean():7.3f}'
              f'{d02[m].mean():8.2f}{ERF[m].mean():8.2f}{r:7.2f}x')
    near = [v for k, v in vals.items() if k[1] <= 12]
    far = [v for k, v in vals.items() if k[0] >= 21]
    if near and far:
        print(f'\n  near (dN<=12) {np.mean(near):.2f}x   far (dN>=21) {np.mean(far):.2f}x   '
              f'near/far {np.mean(near)/np.mean(far):.2f}')
        print('  parent domain: ERF near/far 4.49x, d02 1.83x, MRMS 1.76x')


if __name__ == '__main__':
    main()
