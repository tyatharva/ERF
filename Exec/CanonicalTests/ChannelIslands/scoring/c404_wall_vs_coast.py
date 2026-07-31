#!/usr/bin/env python3
"""Is the d=5-9 precipitation spike a boundary artifact or a windward response?

  c404_wall_vs_coast.py

Post-processing only. The spike: all four 2020-12-28 arms reach 1.24-1.52 in at
d = 5-9 cells from the domain edge while d02 is flat at 0.35-0.46 there.

Three explanations, confounded because the Transverse Range coast runs nearly
parallel to the northern wall:
  (a) boundary artifact  -> tracks distance from the DOMAIN EDGE
  (b) windward response  -> tracks UPSLOPE EXPOSURE / distance inland
  (c) orographic         -> tracks TERRAIN HEIGHT

Everything is binned on the ERF/reference RATIO, which divides out the storm's
own spatial structure (Santa Barbara was hit far harder than the southern Bight,
and that is real, not model error).

CAUTION, learned the hard way here: do NOT use d = min(distance to any wall).
The north (yhi) and east (xhi) walls behave OPPOSITELY -- 3.49x wet against
1.15x near-neutral -- so the min collapses two different populations into one
number and manufactures a "wall effect" out of their average. Every wall test
below is per-wall.

Writes figs/c404_wall_vs_coast.png.
"""
import os

import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.ndimage import distance_transform_edt as edt

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
OUT = '/app/ERF/figs'
os.makedirs(OUT, exist_ok=True)
S = '/app/ERF/scoring_ab_dav/'
MM2IN = 1.0 / 25.4
DXM = 3000.0
FLOW_DEG = 242.0          # d02 mean sub-1500 m direction over the ranges, h13-h17

lat = np.load(S + 'lat.npy'); lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy')
NX, NY = terr.shape
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
DN, DE = NY - 1 - jj, NX - 1 - ii          # distance from the north / east walls
DWALL = np.minimum.reduce([ii, jj, DE, DN])
LAND = terr > terr.min() + 20.0
DCOAST = edt(LAND)

gx, gy = np.gradient(terr, DXM, DXM)
_t = np.deg2rad(FLOW_DEG)
UPS = gx * (-np.sin(_t)) + gy * (-np.cos(_t))

d02 = np.load('/app/ERF/wrf_d02_on_grid.npy') * MM2IN
mrms = np.load('/app/ERF/mrms_20201228_on_grid.npy') * MM2IN


def last_clean(run):
    for p in reversed(plotfiles_by_time(f'/app/ERF/{run}')):
        ds = yt.load(p)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not is_poisoned(ra):
            return ra * MM2IN
    raise SystemExit(f'no clean plotfile in {run}')


NAMES = ['sigma=1', 'sigma=0.1', 'MYNN25', 'w_damp']
ARMS = [last_clean(r) for r in ('run_c404_sig1', 'run_c404_sig01',
                                'run_c404_mynn25', 'run_c404_wdamp')]
ERF = np.mean(ARMS, axis=0)
BASE = LAND & np.isfinite(d02) & (d02 > 0.1)

CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
fig, axes = plt.subplots(2, 3, figsize=(19.2, 10.2), constrained_layout=True)

# (a) the requested plot: references vs distance from the edge, terrain fixed
ax = axes[0, 0]
WB = [(0, 4), (5, 9), (10, 14), (15, 24), (25, 47)]
WX = [2, 7, 12, 19.5, 36]
for c, (hlo, hhi) in zip(CAT, [(100, 300), (300, 600), (600, 1000), (1000, 3000)]):
    for fld, ls, ms in ((d02, '-', 'o'), (mrms, ':', 's')):
        xs, ys = [], []
        for (wlo, whi), x in zip(WB, WX):
            m = LAND & (terr >= hlo) & (terr < hhi) & (DWALL >= wlo) & (DWALL <= whi) & np.isfinite(d02)
            if m.sum() >= 15:
                xs.append(x); ys.append(np.nanmean(fld[m]))
        ax.plot(xs, ys, ms + ls, color=c, lw=2.0 if ls == '-' else 1.3, ms=6,
                label=f'{hlo}-{hhi} m' if ls == '-' else None, alpha=1.0 if ls == '-' else 0.7)
ax.axvspan(5, 9, color='#888888', alpha=0.16, lw=0)
ax.set_xlabel('distance from nearest domain edge (cells)')
ax.set_ylabel('23-h mean precipitation (in)')
ax.set_ylim(0, 1.8)
ax.set_title('(a) THE REFERENCES ARE FLAT across the spike\n'
             'solid d02, dotted MRMS -- no near-wall enhancement at any terrain height',
             fontsize=10)
ax.grid(alpha=0.25, lw=0.6); ax.legend(fontsize=7.5, frameon=False, ncol=2)

# (b) a real windward maximum exists -- inland, in both references
ax = axes[0, 1]
for nm, f, c in (('d02', d02, CAT[0]), ('MRMS', mrms, CAT[1]), ('ERF (4-arm mean)', ERF, CAT[2])):
    xs, ys = [], []
    for c_ in range(1, 16):
        m = LAND & (DCOAST == c_) & np.isfinite(d02)
        if m.sum() >= 10:
            xs.append(c_ * 3); ys.append(np.nanmean(f[m]))
    ax.plot(xs, ys, 'o-', color=c, lw=2.2, ms=6, label=nm)
ax.axvline(33, color='#444444', lw=1.0, ls='--')
ax.text(34, 3.6, 'reference windward max\n33 km inland', fontsize=8, color='#444444')
ax.set_xlabel('distance inland from the coast (km)')
ax.set_ylabel('23-h mean precipitation (in)')
ax.set_title('(b) A REAL windward maximum DOES exist\n'
             'both references peak 24-33 km inland at ~2.5x the immediate coast', fontsize=10)
ax.grid(alpha=0.25, lw=0.6); ax.legend(fontsize=8, frameon=False)

# (c) per-wall -- the north and east walls behave oppositely
ax = axes[0, 2]
S1 = BASE & (terr >= 100)
for c, (nm, dd) in zip(CAT, [('north (yhi)', DN), ('east (xhi)', DE)]):
    xs, ys = [], []
    for (lo, hi), x in zip([(0, 4), (5, 9), (10, 14), (15, 24), (25, 47)], WX):
        m = S1 & (dd >= lo) & (dd <= hi)
        if m.sum() >= 15:
            xs.append(x); ys.append(np.nanmean(ERF[m]) / np.nanmean(d02[m]))
    ax.plot(xs, ys, 'o-', color=c, lw=2.4, ms=8, label=nm)
ax.axhline(1.0, color='#444444', lw=1.0, ls='--')
ax.axvspan(5, 9, color='#888888', alpha=0.16, lw=0)
ax.set_xlabel('distance from that wall (cells)')
ax.set_ylabel('ERF / d02')
ax.set_title('(c) IT IS ONE WALL, NOT "walls"\n'
             'north peaks at 4.1x; east is flat near 1.2x -- internal control', fontsize=10)
ax.grid(alpha=0.25, lw=0.6); ax.legend(fontsize=9, frameon=False)

# (d) THE CONTROLLED TEST
ax = axes[1, 0]
CTRL = (BASE & (terr >= 300) & (terr < 1000) & (UPS > 0.01) & (UPS < 0.06)
        & (DCOAST >= 6) & (DCOAST <= 15))
DB = [(0, 6), (7, 12), (13, 19), (20, 34)]
DBX = [3, 9.5, 16, 27]
for nm, f, c, lw in (('ERF (4-arm mean)', ERF, CAT[2], 2.6),
                     ('d02', d02, CAT[0], 2.0), ('MRMS', mrms, CAT[1], 2.0)):
    xs, ys = [], []
    for (lo, hi), x in zip(DB, DBX):
        m = CTRL & (DN >= lo) & (DN <= hi)
        if m.sum() >= 12:
            xs.append(x); ys.append(np.nanmean(f[m]))
    ax.plot(xs, ys, 'o-', color=c, lw=lw, ms=8, label=nm)
for f, c in zip(ARMS, ['#999999'] * 4):
    xs, ys = [], []
    for (lo, hi), x in zip(DB, DBX):
        m = CTRL & (DN >= lo) & (DN <= hi)
        if m.sum() >= 12:
            xs.append(x); ys.append(np.nanmean(f[m]))
    ax.plot(xs, ys, '-', color=c, lw=0.9, alpha=0.7, zorder=0)
ax.set_xlabel('distance from the NORTH wall (cells)')
ax.set_ylabel('23-h mean precipitation (in)')
ax.set_title('(d) THE TEST -- terrain 300-1000 m, upslope 0.01-0.06,\n'
             '18-45 km inland ALL HELD FIXED; only wall distance varies\n'
             'ERF falls 4.5x, references only 1.8x', fontsize=10)
ax.grid(alpha=0.25, lw=0.6); ax.legend(fontsize=8, frameon=False)

# (e) the same, as a ratio
ax = axes[1, 1]
for nm, ref, c in (('ERF / d02', d02, CAT[0]), ('ERF / MRMS', mrms, CAT[1])):
    xs, ys = [], []
    for (lo, hi), x in zip(DB, DBX):
        m = CTRL & (DN >= lo) & (DN <= hi)
        if m.sum() >= 12:
            xs.append(x); ys.append(np.nanmean(ERF[m]) / np.nanmean(ref[m]))
    ax.plot(xs, ys, 'o-', color=c, lw=2.4, ms=8, label=nm)
ax.axhline(1.0, color='#444444', lw=1.0, ls='--')
ax.axvline(12.5, color='#444444', lw=1.0, ls=':')
ax.text(13, 3.4, 'break at\ndN ~ 12-13\n(36-39 km)', fontsize=8, color='#444444')
ax.set_xlabel('distance from the NORTH wall (cells)')
ax.set_ylabel('ratio to reference')
ax.set_title('(e) Both references agree, and the cutoff is sharp\n'
             '3.3-4.0x inside, 1.1-1.3x outside', fontsize=10)
ax.grid(alpha=0.25, lw=0.6); ax.legend(fontsize=9, frameon=False)

# (f) the map
ax = axes[1, 2]
ratio = np.where(BASE, ERF / np.maximum(d02, 1e-9), np.nan)
m_ = ax.pcolormesh(lon, lat, np.log2(ratio), cmap='RdBu_r',
                   norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2),
                   shading='nearest', zorder=1)
ax.contour(lon, lat, LAND.astype(float), levels=[0.5], colors='k', linewidths=1.0, zorder=3)
ax.contour(lon, lat, DN.astype(float), levels=[12.5], colors='#00a000', linewidths=2.0, zorder=4)
ax.set_xlim(lon.min(), lon.max()); ax.set_ylim(lat.min(), lat.max())
ax.set_xlabel('lon'); ax.set_ylabel('lat')
ax.set_title('(f) ERF/d02 over land, log2\n'
             'green = 12 cells from the NORTH wall; black = coast', fontsize=10)
cb = fig.colorbar(m_, ax=ax, extend='both')
cb.set_ticks([-2, -1, 0, 1, 2]); cb.set_ticklabels(['0.25x', '0.5x', '1x', '2x', '4x'])

fig.suptitle('The d=5-9 spike is the NORTH WALL, not a windward response.  It survives holding terrain, '
             'upslope exposure and distance-inland fixed;\nit is absent at the east wall; and both references '
             'show only a 1.8x real enhancement where ERF shows 4.5x.', fontsize=11)
p = f'{OUT}/c404_wall_vs_coast.png'
fig.savefig(p, dpi=130); plt.close(fig)
print('wrote', p)
