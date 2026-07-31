#!/usr/bin/env python3
"""23-h precipitation field comparison for the 2020-12-28 CONUS404 arms. INCHES.

  c404_fields.py

Same treatment as the Jan-9 `make_maps.py`: one shared 192x96 LCC grid so panels
compare cell for cell, coastline and Channel Islands taken from the MODEL's own
terrain (no network data -- it shows exactly what the model treats as land),
d=20 interior boundary solid white, d=10 relaxation-band edge dotted white.

Three figures:
  c404_fields_0-3in.png / 0-12in.png  2x3 field panels
  c404_diff_vs_d02.png                each arm minus d02, diverging
  c404_profiles_spectra.png           wall-normal profile + radial spectra

MATCHED TIME. Every arm is read from its last NaN-free plotfile (plt_guard), so
arms that took a degenerate zero-length final step are read one plotfile earlier
(item 55b). That leaves a 15.89 s spread across the four arms -- 82800.11 to
82816.00 s, or 0.019% of the 23 h accumulation. Panel titles carry each arm's
own time so the mismatch is visible rather than assumed away.

The sigma=0.03 BASELINE is not reproducible here: `run_c404_nsc/` is not in this
checkout (the pod transfer bundle carried driving frames and reference arrays,
not run output). Its panel is drawn as an explicit placeholder rather than
dropped, so the layout still reads against HANDOFF section 2.
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
os.makedirs(OUT, exist_ok=True)
S = '/app/ERF/scoring_ab_dav/'
NX, NY, DX = 192, 96, 3.0
MM2IN = 1.0 / 25.4

lat = np.load(S + 'lat.npy')
lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy')
LAND = terr > terr.min() + 20.0
d02 = np.load('/app/ERF/wrf_d02_on_grid.npy') * MM2IN

ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
DWALL = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
INTERIOR = DWALL >= 20
HI = (terr >= 500) & LAND

# Categorical slots 1-6 of the validated default order, used in fixed order.
CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']


def last_clean(run):
    """Last NaN-free plotfile in `run` -- never a degenerate-final-step one."""
    for p in reversed(plotfiles_by_time(f'/app/ERF/{run}')):
        ds = yt.load(p)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not is_poisoned(ra):
            return ra * MM2IN, float(ds.current_time), os.path.basename(p)
    sys.exit(f'FATAL: no clean plotfile in {run}')


ARMS = [('d02 (reference, 1.5 km)', d02, None, None),
        ('sigma=0.03 baseline', None, None, None)]
for label, run in [('sigma=1', 'run_c404_sig1'),
                   ('sigma=0.1', 'run_c404_sig01'),
                   ('MYNN25 (sigma=0.03)', 'run_c404_mynn25'),
                   ('w_damp cfl=0.15 (sigma=0.03)', 'run_c404_wdamp')]:
    f, t, pf = last_clean(run)
    ARMS.append((label, f, t, pf))
    print(f'{label:30s} {pf:12s} t={t:9.2f} s  max {np.nanmax(f):5.2f} in')

times = [t for _, f, t, _ in ARMS if t is not None]
SPREAD = max(times) - min(times)
print(f'matched-time spread across arms: {SPREAD:.2f} s '
      f'({100 * SPREAD / 82800:.3f}% of the accumulation)')


def decorate(ax, first_col):
    ax.contour(lon, lat, LAND.astype(float), levels=[0.5], colors='k',
               linewidths=0.6, zorder=3)
    ax.contour(lon, lat, DWALL.astype(float), levels=[19.5], colors='w',
               linewidths=1.1, zorder=4)
    ax.contour(lon, lat, DWALL.astype(float), levels=[9.5], colors='w',
               linewidths=0.9, linestyles=':', zorder=4)
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())
    ax.set_xlabel('lon', fontsize=8)
    ax.tick_params(labelsize=7)
    if first_col:
        ax.set_ylabel('lat', fontsize=8)
    else:
        ax.set_yticklabels([])


def missing_panel(ax, note):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linestyle('--'); s.set_color('#888888')
    ax.text(0.5, 0.5, note, ha='center', va='center', fontsize=9,
            color='#555555', wrap=True)


# ------------------------------------------------------------ fields --------
def field_figure(vmax):
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 8.4), constrained_layout=True)
    cmap = plt.get_cmap('viridis').copy(); cmap.set_over(cmap(1.0))
    norm = Normalize(0, vmax)
    m = None
    for k, (name, f, t, pf) in enumerate(ARMS):
        ax = axes[k // 3, k % 3]
        if f is None:
            missing_panel(ax, 'sigma = 0.03 baseline\n\nNOT AVAILABLE on this pod\n'
                              'run_c404_nsc/ is not in this checkout\n\n'
                              'published values (HANDOFF S2):\n'
                              'interior 1.03x d02,  >500 m 2.21x\n'
                              'Santa Ynez 6.07 in,  domain max 17.47 in')
            ax.set_title('sigma=0.03 baseline', fontsize=9)
            continue
        m = ax.pcolormesh(lon, lat, np.where(np.isfinite(f), f, np.nan),
                          cmap=cmap, norm=norm, shading='nearest', zorder=1)
        decorate(ax, k % 3 == 0)
        sub = '' if t is None else f'   t={t:.1f} s'
        ax.set_title(f'{name}{sub}\ninterior mean {np.nanmean(f[INTERIOR]):.2f} in   '
                     f'>500 m mean {np.nanmean(f[HI]):.2f} in\n'
                     f'domain max {np.nanmax(f):.2f} in', fontsize=9)
    cb = fig.colorbar(m, ax=axes, extend='max', shrink=0.7)
    cb.set_label(f'23-h precipitation (in), 0-{vmax:g}, above saturates')
    fig.suptitle('black = model land mask (coast + Channel Islands)   |   solid white = d=20 '
                 'interior boundary   |   dotted white = d=10 relaxation-band edge\n'
                 f'all arms from their last NaN-free plotfile; matched-time spread '
                 f'{SPREAD:.2f} s ({100*SPREAD/82800:.3f}% of the accumulation)',
                 fontsize=9)
    p = f'{OUT}/c404_fields_0-{vmax:g}in.png'
    fig.savefig(p, dpi=140); plt.close(fig)
    print('wrote', p)


# -------------------------------------------------------- differences -------
def diff_figure():
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 8.4), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-3.0, vcenter=0.0, vmax=3.0)
    m = None
    for k, (name, f, t, pf) in enumerate(ARMS):
        ax = axes[k // 3, k % 3]
        if k == 0:
            missing_panel(ax, 'd02 is the REFERENCE\n\nevery other panel is\n'
                              'arm minus d02\n\n'
                              f'd02 interior mean {np.nanmean(d02[INTERIOR]):.2f} in\n'
                              f'd02 >500 m mean {np.nanmean(d02[HI]):.2f} in')
            ax.set_title('d02 (reference)', fontsize=9)
            continue
        if f is None:
            missing_panel(ax, 'sigma = 0.03 baseline\n\nNOT AVAILABLE on this pod\n'
                              'run_c404_nsc/ is not in this checkout')
            ax.set_title('sigma=0.03 baseline - d02', fontsize=9)
            continue
        dd = f - d02
        m = ax.pcolormesh(lon, lat, dd, cmap='RdBu_r', norm=norm,
                          shading='nearest', zorder=1)
        decorate(ax, k % 3 == 0)
        ax.set_title(f'{name} - d02\ninterior mean diff '
                     f'{np.nanmean(dd[INTERIOR]):+.2f} in   '
                     f'>500 m {np.nanmean(dd[HI]):+.2f} in\n'
                     f'RMSE {np.sqrt(np.nanmean(dd[INTERIOR]**2)):.2f} in (interior)',
                     fontsize=9)
    cb = fig.colorbar(m, ax=axes, extend='both', shrink=0.7)
    cb.set_label('arm minus d02 (in);  red = ERF wetter,  blue = ERF drier')
    fig.suptitle('Difference from d02, diverging about zero, same panel slots as the field figure\n'
                 'black = model land mask   |   solid white = d=20 interior boundary   |   '
                 'dotted white = d=10 relaxation-band edge', fontsize=9)
    p = f'{OUT}/c404_diff_vs_d02.png'
    fig.savefig(p, dpi=140); plt.close(fig)
    print('wrote', p)


# ------------------------------------------- profiles and spectra -----------
def radial_spectrum(x, m):
    x = np.where(m, np.nan_to_num(x), 0.0)
    c = x[20:NX - 20, 20:NY - 20]
    c = c - c.mean()
    win = np.hanning(c.shape[0])[:, None] * np.hanning(c.shape[1])[None, :]
    F = np.abs(np.fft.rfft2(c * win)) ** 2
    kx = np.fft.fftfreq(c.shape[0], DX)[:, None]
    ky = np.fft.rfftfreq(c.shape[1], DX)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    b = np.linspace(0, k.max(), 40)
    idx = np.digitize(k.ravel(), b)
    P = np.array([F.ravel()[idx == i].mean() if (idx == i).any() else np.nan
                  for i in range(1, len(b))])
    return 1.0 / (0.5 * (b[1:] + b[:-1])), P


def profile_figure():
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), constrained_layout=True)
    present = [(n, f) for n, f, _, _ in ARMS if f is not None]

    ax = axes[0]
    ds = np.arange(0, 48)
    for c, (name, f) in zip(CAT, present):
        prof = [np.nanmean(f[DWALL == d]) for d in ds]
        ax.plot(ds, prof, lw=2.0, color=c, label=name,
                ls='--' if name.startswith('d02') else '-')
    ax.axvline(9.5, color='#888888', lw=0.9, ls=':')
    ax.axvline(19.5, color='#888888', lw=1.1)
    ax.text(9.7, ax.get_ylim()[1] * 0.97, 'd=10 band edge', fontsize=7,
            rotation=90, va='top', color='#555555')
    ax.text(19.7, ax.get_ylim()[1] * 0.97, 'd=20 interior', fontsize=7,
            rotation=90, va='top', color='#555555')
    ax.set_xlabel('d = distance from nearest wall (cells; 1 cell = 3 km)')
    ax.set_ylabel('23-h mean precipitation (in)')
    ax.set_title('Wall-normal profile\nmean over all cells at each distance from the boundary',
                 fontsize=10)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    lam_r, Pr = radial_spectrum(d02, np.isfinite(d02))
    for c, (name, f) in zip(CAT, present):
        lam, P = radial_spectrum(f, np.isfinite(d02))
        ax.loglog(lam, P, lw=2.0, color=c, label=name,
                  ls='--' if name.startswith('d02') else '-')
    ax.axvspan(8, 19, color='#888888', alpha=0.13, lw=0)
    ax.text(12, ax.get_ylim()[1] * 0.3, '8-19 km\nscored band', fontsize=7,
            ha='center', color='#555555')
    ax.set_xlabel('wavelength (km)')
    ax.set_ylabel('radially averaged power')
    ax.set_title('Radial power spectrum, interior only\nratio over 8-19 km is the '
                 'scored spectrum metric', fontsize=10)
    ax.grid(alpha=0.25, lw=0.6, which='both')
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle('sigma=0.03 baseline absent from both panels -- run_c404_nsc/ is not in this checkout',
                 fontsize=9)
    p = f'{OUT}/c404_profiles_spectra.png'
    fig.savefig(p, dpi=140); plt.close(fig)
    print('wrote', p)


field_figure(3.0)
field_figure(12.0)
diff_figure()
profile_figure()
