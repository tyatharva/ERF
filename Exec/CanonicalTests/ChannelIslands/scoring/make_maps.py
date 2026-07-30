#!/usr/bin/env python3
"""24-h precipitation maps and diagnostics, Jan-9 2023 arms. Units: INCHES.

All fields share one 192x96 LCC grid, so panels compare cell for cell. Coastline
and islands come from the MODEL's own terrain (z_phys at k=0 exceeds the 12.5 m
ocean value over land) -- no network data, and it shows exactly what the model
treats as land.

Two copies of the field figure: 0-12 in spans the full range including the
offshore MRMS Pass-2 gap-filled maximum; 0-3 in resolves the scored interior,
whose observed maximum is 70 mm = 2.8 in.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.ndimage import gaussian_filter, label, center_of_mass

OUT = '/app/ERF/figs'
os.makedirs(OUT, exist_ok=True)
S = '/app/ERF/scoring_ab_'
NX, NY, DX = 192, 96, 3.0
MM2IN = 1.0 / 25.4

lat = np.load(S + 'dav/lat.npy'); lon = np.load(S + 'dav/lon.npy')
terr = np.load(S + 'dav/terrain.npy')
mrms = np.load(S + 'dav/mrms_mm.npy') * MM2IN
era5 = np.load(S + 'dav/era5_mm.npy') * MM2IN
land = terr > terr.min() + 20.0

ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
dwall = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
interior = dwall >= 20

arms = [('MRMS (observed)', mrms), ('ERA5 (driver, 25 km)', era5),
        ('Davies', np.load(S + 'dav/erf_mm.npy') * MM2IN),
        ('NSCBC no band (sig=0.03)', np.load(S + 'sig003/erf_mm.npy') * MM2IN)]
if os.path.exists(S + 'rw2/erf_mm.npy'):
    arms.append(('NSCBC + 2-cell band', np.load(S + 'rw2/erf_mm.npy') * MM2IN))


def decorate(ax, first=False):
    ax.contour(lon, lat, land.astype(float), levels=[0.5], colors='k',
               linewidths=0.6, zorder=3)
    ax.contour(lon, lat, dwall.astype(float), levels=[19.5], colors='w',
               linewidths=1.1, zorder=4)
    ax.contour(lon, lat, dwall.astype(float), levels=[9.5], colors='w',
               linewidths=0.9, linestyles=':', zorder=4)
    ax.set_xlim(lon.min(), lon.max()); ax.set_ylim(lat.min(), lat.max())
    ax.set_xlabel('lon')
    if first:
        ax.set_ylabel('lat')
    else:
        ax.set_yticklabels([])


def field_figure(vmax, tag):
    n = len(arms)
    fig, axes = plt.subplots(1, n, figsize=(4.1 * n, 4.7), constrained_layout=True)
    cmap = plt.get_cmap('viridis').copy(); cmap.set_over(cmap(1.0))
    norm = Normalize(0, vmax)
    for k, (name, f) in enumerate(arms):
        ax = axes[k]
        m = ax.pcolormesh(lon, lat, np.where(np.isfinite(f), f, np.nan), cmap=cmap,
                          norm=norm, shading='nearest', zorder=1)
        decorate(ax, k == 0)
        ax.set_title(f'{name}\ninterior mean {np.nanmean(f[interior]):.2f} in   '
                     f'interior max {np.nanmax(np.where(interior, f, np.nan)):.2f} in\n'
                     f'domain max {np.nanmax(f):.2f} in', fontsize=9)
    cb = fig.colorbar(m, ax=axes, extend='max', shrink=0.85)
    cb.set_label(f'24-h precip (in), 0-{vmax:g}, above saturates')
    fig.suptitle('black = model land mask (coast + Channel Islands);   solid white = d=20 '
                 'interior boundary;   dotted white = d=10 relaxation-band edge', fontsize=10)
    p = f'{OUT}/precip_fields_0-{vmax:g}in.png'
    fig.savefig(p, dpi=150); plt.close(fig)
    print('wrote', p)


field_figure(12.0, 'full')
field_figure(3.0, 'interior')

# ---------------- diagnostics ----------------
models = arms[1:]
fig = plt.figure(figsize=(4.1 * len(models), 12.8), constrained_layout=True)
gs = fig.add_gridspec(3, len(models), height_ratios=[1.0, 1.0, 0.85])
cmap = plt.get_cmap('viridis').copy(); cmap.set_over(cmap(1.0))
norm3 = Normalize(0, 3.0)

dnorm = TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
for k, (name, f) in enumerate(models):
    ax = fig.add_subplot(gs[0, k])
    m = ax.pcolormesh(lon, lat, f - mrms, cmap='RdBu_r', norm=dnorm,
                      shading='nearest', zorder=1)
    decorate(ax, k == 0)
    ax.set_title(f'{name} - MRMS', fontsize=10)
    if k == len(models) - 1:
        fig.colorbar(m, ax=ax, extend='both', label='inches')

THR_IN = 5.0 * MM2IN


def objects(field, thr=THR_IN):
    sm = gaussian_filter(np.nan_to_num(field), 1.0)
    binm = (sm >= thr) & interior
    lab, nl = label(binm, structure=np.ones((3, 3)))
    keep = np.zeros_like(binm, dtype=float); info = []
    for i in range(1, nl + 1):
        sel = lab == i
        if sel.sum() < 10:
            continue
        keep[sel] = 1.0
        cy, cx = center_of_mass(sel)
        info.append((cy, cx, int(sel.sum())))
    return keep, info


om, oinfo = objects(mrms)
aobs = sum(a for _, _, a in oinfo) or 1
for k, (name, f) in enumerate(models):
    ax = fig.add_subplot(gs[1, k])
    ax.pcolormesh(lon, lat, np.where(np.isfinite(f), f, np.nan), cmap=cmap,
                  norm=norm3, shading='nearest', zorder=1)
    decorate(ax, k == 0)
    ax.contour(lon, lat, om, levels=[0.5], colors='r', linewidths=2.0, zorder=5)
    fm, finfo = objects(f)
    ax.contour(lon, lat, fm, levels=[0.5], colors='cyan', linewidths=1.6, zorder=5)
    for cy, cx, a in oinfo:
        ax.plot(lon[int(cy), int(cx)], lat[int(cy), int(cx)], 'r*', ms=15, zorder=6)
    for cy, cx, a in finfo:
        ax.plot(lon[int(cy), int(cx)], lat[int(cy), int(cx)], 'c*', ms=13, zorder=6)
    dsp = ''
    if oinfo and finfo:
        oy, ox, _ = oinfo[0]
        best = min(finfo, key=lambda t: np.hypot((t[0]-oy)*DX, (t[1]-ox)*DX))
        dsp = f'   displ {np.hypot((best[0]-oy)*DX,(best[1]-ox)*DX):.0f} km'
    ax.set_title(f'MODE objects, 5 mm (0.20 in)\n'
                 f'red MRMS ({len(oinfo)})  cyan model ({len(finfo)})\n'
                 f'area {sum(a for _,_,a in finfo)/aobs:.2f}x{dsp}', fontsize=9)

ax = fig.add_subplot(gs[2, :])
sel = (jj >= 20) & (jj < NY - 20)
xkm = (np.arange(NX) + 0.5) * DX
for name, f in [arms[0]] + models:
    prof = np.array([np.nanmean(f[i, :][sel[i, :]]) for i in range(NX)])
    ax.semilogy(xkm, np.maximum(prof, 1e-4), lw=2.4 if 'MRMS' in name else 1.5,
                ls='-' if 'MRMS' in name else '--', label=name)
for v, s in ((10 * DX, ':'), (20 * DX, '-'), ((NX - 20) * DX, '-'), ((NX - 10) * DX, ':')):
    ax.axvline(v, color='0.4', ls=s, lw=1)
ax.set_xlabel('distance along the inflow (xlo) -> outflow (xhi) axis, km    '
              '[dotted = d=10 band edge, solid = d=20 interior boundary]')
ax.set_ylabel('24-h precipitation, inches\n(mean over interior y)')
ax.set_title('Cross-section: wall band and interior decay (log scale, inches)', fontsize=10)
ax.legend(fontsize=9, ncol=5); ax.grid(alpha=0.3, which='both')
ax.set_xlim(0, NX * DX)
p = f'{OUT}/precip_diagnostics.png'
fig.savefig(p, dpi=150)
print('wrote', p)
