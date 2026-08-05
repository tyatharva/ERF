#!/usr/bin/env python3
"""Accumulation comparison for the 29-h Domain A pair.

  plot_pair_29h.py <out.png> <label>=<rundir> [<label>=<rundir> ...]

Panels: every arm given, then MRMS and WRF d02. Works with one arm (Davies
alone, after the first run) or two (the final Davies-vs-NSCBC figure).

WINDOW. The run starts 12/27 18Z, so the scored window 12/28 00Z-23Z is
h6 -> h29 of model time. Both references are 23-h accumulations ending 23Z,
which is why the run stops there and why nothing later is plotted.

The accumulation is rain_accum(h29) - rain_accum(h6): a DIFFERENCE of two
plotfiles, not the final field, because rain_accum runs from t=0 and the first
6 h are spin-up that neither reference covers.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt

from plt_guard import plotfiles_by_time, reject_if_poisoned

yt.set_log_level(50)

H0, H1 = 6.0, 12.0          # model hours bounding the scored window (00Z-06Z)
TOL = 0.02                  # hours
NX, NY, BAND = 192, 96, 10
REFS = '/app/ERF/refs'


def window(label, rundir):
    """rain_accum(H1) - rain_accum(H0) for one arm."""
    got = {}
    for p in plotfiles_by_time(rundir):
        ds = yt.load(p)
        h = float(ds.current_time) / 3600.0
        for hr in (H0, H1):
            if abs(h - hr) < TOL and hr not in got:
                g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
                ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
                reject_if_poisoned(f'{label} h{hr:g}', p, ra)
                got[hr] = ra
    missing = [h for h in (H0, H1) if h not in got]
    if missing:
        raise SystemExit(f'{label}: no plotfile at h{missing} in {rundir} -- '
                         f'the run did not reach the scored window')
    return np.maximum(got[H1] - got[H0], 0.0)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    outp = sys.argv[1]
    arms = [a.split('=', 1) for a in sys.argv[2:]]

    d02 = np.load(f'{REFS}/domA_d02_h0_h6.npy')
    mrms = np.load(f'{REFS}/domA_mrms_h0_h6.npy')
    terr = np.load(f'{REFS}/domA_terrain.npy')
    LAND = terr > 20.0

    panels = [(lab, window(lab, rd)) for lab, rd in arms]
    panels += [('MRMS  (observed)', mrms), ('WRF d02  (1.5 km)', d02)]

    # Interior mask: the outer real_width cells are forcing-dominated and are
    # excluded from every scored number in this campaign.
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dmin = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    M = LAND & np.isfinite(d02) & np.isfinite(mrms) & (dmin >= BAND)

    n = len(panels)
    ncol = 2
    nrow = (n + 1) // 2
    fig, ax = plt.subplots(nrow, ncol, figsize=(13.5, 4.6 * nrow), squeeze=False)
    vmax = 20.0
    for a, (lab, arr) in zip(ax.flat, panels):
        p = a.imshow(arr.T, origin='lower', vmin=0, vmax=vmax, cmap='turbo',
                     aspect='auto')
        # Report the INTERIOR max, not the domain max. The two differ by more
        # than an order of magnitude for the ERF arms -- Davies peaks at 3565 mm
        # ON the xhi wall (i=191, zero cells from the boundary) while its
        # interior max is 250 mm -- so a domain max next to an interior mean
        # invites reading a boundary artifact as a physical result. Every cell
        # above 500 mm in this run is inside the excluded band.
        a.set_title(f'{lab}\ninterior land mean {np.nanmean(arr[M]):.2f} mm   |   '
                    f'interior max {np.nanmax(arr[M]):.0f} mm   |   '
                    f'band max {np.nanmax(arr[~M & LAND]) if (~M & LAND).any() else float("nan"):.0f} mm',
                    fontsize=9)
        plt.colorbar(p, ax=a, label='mm / 6 h', extend='max')
        a.contour(LAND.T.astype(float), levels=[0.5], colors='w', linewidths=0.6)
        for e in (BAND, NX - 1 - BAND):
            a.axvline(e, color='k', ls=':', lw=1.0)
        for e in (BAND, NY - 1 - BAND):
            a.axhline(e, color='k', ls=':', lw=1.0)
        a.set_xlabel('i'); a.set_ylabel('j')
    for a in ax.flat[n:]:
        a.axis('off')

    fig.suptitle('Domain A, COMPRESSIBLE, 2020-12-28 00Z-06Z (h6-h12 of a 12 h run)\n'
                 'interior land only in the statistics; dotted lines = the '
                 f'{BAND}-cell forcing-dominated band, excluded', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    fig.savefig(outp, dpi=130)
    print(f'wrote {outp}')

    print(f'\n{"panel":26s} {"int-land mean":>14s} {"int max":>9s} '
          f'{"band max":>9s} {"vs d02":>8s} {"vs MRMS":>9s}')
    dm, mm = np.nanmean(d02[M]), np.nanmean(mrms[M])
    BND = (~M) & LAND
    for lab, arr in panels:
        m = np.nanmean(arr[M])
        bmax = np.nanmax(arr[BND]) if BND.any() else float('nan')
        print(f'{lab:26s} {m:14.2f} {np.nanmax(arr[M]):9.0f} {bmax:9.0f} '
              f'{m/dm:7.2f}x {m/mm:8.2f}x')
    # A matching mean is not a matching field -- the anelastic Davies arm scored
    # 1.09x on the 23 h window while its pattern correlation against MRMS was
    # only +0.18, because a too-dry median and a too-wet tail cancelled.
    print(f'\n{"panel":26s} {"corr vs d02":>12s} {"corr vs MRMS":>13s} '
          f'{"median":>8s} {"p99":>8s}')
    for lab, arr in panels:
        cd = np.corrcoef(arr[M], d02[M])[0, 1]
        cm = np.corrcoef(arr[M], mrms[M])[0, 1]
        print(f'{lab:26s} {cd:+12.3f} {cm:+13.3f} '
              f'{np.median(arr[M]):8.2f} {np.percentile(arr[M], 99):8.2f}')


if __name__ == '__main__':
    main()
