#!/usr/bin/env python3
"""Item-56 figure set: Domain A (shifted north) and Domain B (all-ocean control).

  item56_figs.py [domA|domB] ...

Same treatment as the parent-domain figures: linear inches, 0-3 and 0-12 in,
model land mask as coastline, d=20 interior boundary solid white, d=10
relaxation-band edge dotted white. Adds terrain contours so the position of the
ranges relative to the wall is readable, which is the whole point of Domain A.

Terrain comes from the plotfile's own z_phys at k=0, not from the terrain file:
that is what the model actually ran on, and it carries the 12.5 m first-cell
offset the checked-in parent array uses (ocean reads 12.50 exactly).

Writes into /app/ERF/figs/.
"""
import os
import sys

import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
OUT = '/app/ERF/figs'
REFS = '/app/ERF/refs'
os.makedirs(OUT, exist_ok=True)
MM2IN = 1.0 / 25.4
NX, NY, DX = 192, 96, 3.0
OCEAN_Z = 12.5

CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']

DOMAINS = {
    'domA': dict(title='Domain A -- NE pin 35.4 / -117.2 (shifted north)',
                 runs=[('sigma=0.1', 'run_domA_sig01'), ('sigma=0.03', 'run_domA_sig003')],
                 primary='d02'),
    'domB': dict(title='Domain B -- NE pin 32.3 / -119.3 (all-ocean control)',
                 runs=[('sigma=0.1', 'run_domB_sig01'), ('sigma=0.03', 'run_domB_sig003')],
                 primary='c404'),
}


def load_arm(run):
    for p in reversed(plotfiles_by_time(f'/app/ERF/{run}')):
        ds = yt.load(p)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not is_poisoned(ra):
            z = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
            return ra * MM2IN, z, float(ds.current_time), os.path.basename(p)
    raise SystemExit(f'no clean plotfile in {run}')


def geom():
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    return dict(ii=ii, jj=jj, dW=ii, dS=jj, dE=NX - 1 - ii, dN=NY - 1 - jj,
                dmin=np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj]))


def decorate(ax, lon, lat, land, dmin, terr, first):
    ax.contour(lon, lat, land.astype(float), levels=[0.5], colors='k', linewidths=0.7, zorder=3)
    ax.contour(lon, lat, dmin.astype(float), levels=[19.5], colors='w', linewidths=1.2, zorder=4)
    ax.contour(lon, lat, dmin.astype(float), levels=[9.5], colors='w', linewidths=0.9,
               linestyles=':', zorder=4)
    if terr.max() > 100:
        ax.contour(lon, lat, terr, levels=[500, 1000, 1500, 2000], colors='#00000060',
                   linewidths=0.5, zorder=2)
    ax.set_xlim(lon.min(), lon.max()); ax.set_ylim(lat.min(), lat.max())
    ax.set_xlabel('lon', fontsize=8); ax.tick_params(labelsize=7)
    if first:
        ax.set_ylabel('lat', fontsize=8)
    else:
        ax.set_yticklabels([])


def build(tag):
    D = DOMAINS[tag]
    lat = np.load(f'{REFS}/{tag}_lat.npy'); lon = np.load(f'{REFS}/{tag}_lon.npy')
    d02 = np.load(f'{REFS}/{tag}_d02.npy') * MM2IN
    c404 = np.load(f'{REFS}/{tag}_c404.npy') * MM2IN
    G = geom()
    arms = [(lab,) + load_arm(run) for lab, run in D['runs']]
    terr = arms[0][2] - OCEAN_Z
    land = terr > 20.0
    panels = [('d02 (1.5 km)', d02), ('CONUS404 (4 km, driver)', c404)] + \
             [(f'{lab}   t={t:.0f}s', ra) for lab, ra, _, t, _ in
              [(a[0], a[1], a[2], a[3], a[4]) for a in arms]]

    # 0-3 and 0-12 in are the parent-domain treatment and are kept for
    # comparability, but Domain B tops out near 0.4 in for the references and
    # renders as a black rectangle at those limits. Add a scale matched to the
    # data so the spatial structure is actually readable.
    scales = (3.0, 12.0) if tag == 'domA' else (0.5, 3.0, 12.0)
    for vmax in scales:
        fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0), constrained_layout=True)
        cmap = plt.get_cmap('viridis').copy(); cmap.set_over(cmap(1.0))
        norm = Normalize(0, vmax); m = None
        for k, (name, f) in enumerate(panels):
            ax = axes[k // 2, k % 2]
            m = ax.pcolormesh(lon, lat, np.where(np.isfinite(f), f, np.nan), cmap=cmap,
                              norm=norm, shading='nearest', zorder=1)
            decorate(ax, lon, lat, land, G['dmin'], terr, k % 2 == 0)
            cov = 100.0 * np.isfinite(f).mean()
            sub = '' if cov > 99.9 else f'  [{cov:.0f}% coverage]'
            ax.set_title(f'{name}{sub}\nmean {np.nanmean(f):.3f} in   max {np.nanmax(f):.2f} in',
                         fontsize=9)
        cb = fig.colorbar(m, ax=axes, extend='max', shrink=0.7)
        cb.set_label(f'23-h precipitation (in), 0-{vmax:g}, above saturates')
        fig.suptitle(f'{D["title"]}\nblack = model land mask   |   white solid = d=20   |   '
                     f'white dotted = d=10 band edge   |   grey = terrain 500/1000/1500/2000 m',
                     fontsize=10)
        p = f'{OUT}/item56_{tag}_fields_0-{vmax:g}in.png'
        fig.savefig(p, dpi=150); plt.close(fig); print('wrote', p)

    # difference vs the domain's primary reference
    ref = d02 if D['primary'] == 'd02' else c404
    refname = 'd02' if D['primary'] == 'd02' else 'CONUS404'
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0), constrained_layout=True)
    vlim = 3.0 if tag == 'domA' else 0.6
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0.0, vmax=vlim)
    for k, (lab, ra, _, t, _) in enumerate(arms):
        ax = axes[k]
        dd = ra - ref
        m = ax.pcolormesh(lon, lat, dd, cmap='RdBu_r', norm=norm, shading='nearest', zorder=1)
        decorate(ax, lon, lat, land, G['dmin'], terr, k == 0)
        ax.set_title(f'{lab} - {refname}\nmean {np.nanmean(dd):+.3f} in   '
                     f'RMSE {np.sqrt(np.nanmean(dd**2)):.3f} in', fontsize=9)
    cb = fig.colorbar(m, ax=axes, extend='both', shrink=0.9)
    cb.set_label(f'arm minus {refname} (in);  red = ERF wetter')
    fig.suptitle(f'{D["title"]}  --  difference from {refname}', fontsize=10)
    p = f'{OUT}/item56_{tag}_diff.png'
    fig.savefig(p, dpi=150); plt.close(fig); print('wrote', p)

    # wall-normal profiles, all four walls separately
    fig, axes = plt.subplots(1, 4, figsize=(19.0, 4.6), constrained_layout=True)
    WALLS = [('xlo W (inflow)', G['dW']), ('ylo S (inflow)', G['dS']),
             ('xhi E (outflow)', G['dE']), ('yhi N (outflow)', G['dN'])]
    for k, (wn, d) in enumerate(WALLS):
        ax = axes[k]
        series = [('d02', d02, CAT[0], '--'), ('CONUS404', c404, CAT[1], '--')] + \
                 [(lab, ra, CAT[2 + i], '-') for i, (lab, ra, _, _, _) in enumerate(arms)]
        for nm, f, c, ls in series:
            ds_, ys = [], []
            for dd in range(0, int(d.max()) + 1):
                mm = d == dd
                if mm.sum() and np.isfinite(f[mm]).any():
                    ds_.append(dd); ys.append(np.nanmean(f[mm]))
            ax.plot(ds_, ys, ls, color=c, lw=2.0, label=nm)
        ax.axvline(9.5, color='#888888', lw=0.8, ls=':')
        ax.axvline(19.5, color='#888888', lw=1.0)
        ax.set_xlabel(f'distance from {wn.split()[0]} (cells)')
        if k == 0: ax.set_ylabel('23-h mean precip (in)')
        ax.set_title(wn, fontsize=10)
        ax.grid(alpha=0.25, lw=0.6)
        if k == 0: ax.legend(fontsize=8, frameon=False)
    fig.suptitle(f'{D["title"]}  --  wall-normal profiles, each wall separately '
                 f'(dotted grey d=10, solid grey d=20)', fontsize=10)
    p = f'{OUT}/item56_{tag}_wall_profiles.png'
    fig.savefig(p, dpi=150); plt.close(fig); print('wrote', p)

    # spectra
    def spec(x):
        x = np.where(np.isfinite(x), x, 0.0)
        c = x[20:NX-20, 20:NY-20]; c = c - c.mean()
        F = np.abs(np.fft.rfft2(c * np.hanning(c.shape[0])[:, None]
                                * np.hanning(c.shape[1])[None, :])) ** 2
        kx = np.fft.fftfreq(c.shape[0], DX)[:, None]; ky = np.fft.rfftfreq(c.shape[1], DX)[None, :]
        kk = np.sqrt(kx**2 + ky**2); b = np.linspace(0, kk.max(), 40)
        idx = np.digitize(kk.ravel(), b)
        P = np.array([F.ravel()[idx == i].mean() if (idx == i).any() else np.nan
                      for i in range(1, len(b))])
        return 1.0 / (0.5 * (b[1:] + b[:-1])), P
    fig, ax = plt.subplots(figsize=(7.4, 5.4), constrained_layout=True)
    for nm, f, c, ls in ([('d02', d02, CAT[0], '--'), ('CONUS404', c404, CAT[1], '--')] +
                         [(lab, ra, CAT[2 + i], '-') for i, (lab, ra, _, _, _) in enumerate(arms)]):
        if not np.isfinite(f).all() and np.isfinite(f).mean() < 0.9:
            continue
        lam, P = spec(f)
        ax.loglog(lam, P, ls, color=c, lw=2.0, label=nm)
    ax.axvspan(8, 19, color='#888888', alpha=0.13, lw=0)
    ax.set_xlabel('wavelength (km)'); ax.set_ylabel('radially averaged power')
    ax.set_title(f'{D["title"]}\nradial power spectrum, interior only', fontsize=10)
    ax.grid(alpha=0.25, lw=0.6, which='both'); ax.legend(fontsize=8, frameon=False)
    p = f'{OUT}/item56_{tag}_spectra.png'
    fig.savefig(p, dpi=150); plt.close(fig); print('wrote', p)


for tag in (sys.argv[1:] or ['domA', 'domB']):
    build(tag)
