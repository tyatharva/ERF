#!/usr/bin/env python3
"""Score one or more arms over an ARBITRARY model-hour window.

  score_window.py <out.png> <h0> <h1> <ref_tag> [--lead H] <label>=<rundir> [...]

<h0>/<h1> are WALL-CLOCK hours from 00Z and name the reference files. --lead H
is the spin-up length, added to those to reach MODEL hours when reading
plotfiles: a run starting 12/27 22Z scores wall 00Z-06Z as model hours 2-8, so
    score_window.py out.png 0 6 ocean --lead 2 'ERF'=run_x

<ref_tag> selects refs/<tag>_{d02,mrms}_h<h0>_h<h1>.npy ("domA", "ocean", ...).

Same conventions as plot_pair_29h.py, which this is a parameterized copy of --
they must stay comparable:

  * Accumulation is rain_accum(h1) - rain_accum(h0), a DIFFERENCE of two
    plotfiles. rain_accum runs from t=0, so the final field alone would include
    the spin-up hours the window deliberately excludes.
  * Scored over interior LAND only: terrain > 20 m, and >= BAND cells from every
    lateral boundary. The outer real_width cells are forcing-dominated and are
    excluded from every scored number in this campaign.
  * References are refs/domA_{d02,mrms}_h<h0>_h<h1>.npy, built by
    pod/make_window_refs.py with the same bucket and binning conventions as the
    23-h references.
  * Plotfiles are ordered by MODEL TIME and guarded (plt_guard), so a degraded
    file stops the scoring rather than averaging away.

The headline number is PATTERN CORRELATION against MRMS. The mean is reported
but is not the criterion: it has read 1.09x while the correlation was +0.18.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt

from plt_guard import plotfiles_by_time, reject_if_poisoned

yt.set_log_level(50)

TOL = 0.02
BAND = 10
REFS = '/app/ERF/refs'


def window(label, rundir, h0, h1):
    got = {}
    for p in plotfiles_by_time(rundir):
        ds = yt.load(p)
        h = float(ds.current_time) / 3600.0
        for hr in (h0, h1):
            if abs(h - hr) < TOL and hr not in got:
                g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
                ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
                reject_if_poisoned(f'{label} h{hr:g}', p, ra)
                got[hr] = ra
    missing = [h for h in (h0, h1) if h not in got]
    if missing:
        raise SystemExit(f'{label}: no plotfile at h{missing} in {rundir} -- '
                         f'the run did not reach the scored window')
    return np.maximum(got[h1] - got[h0], 0.0)


def main():
    argv = list(sys.argv[1:])
    # --lead H separates WALL-CLOCK hours from MODEL hours. With a spin-up the
    # two stop being the same number: a run starting 12/27 22Z reaches 00Z at
    # model hour 2, so the 00Z-06Z window is model hours 2-8. h0/h1 stay
    # WALL-CLOCK (they name the reference files); the lead is added only when
    # looking up plotfiles. Getting this backwards silently scores the wrong
    # hours against the right references, which is item 83's failure mode.
    lead = 0.0
    for i, a in enumerate(argv):
        if a == '--lead':
            lead = float(argv[i+1]); del argv[i:i+2]; break
        if a.startswith('--lead='):
            lead = float(a.split('=', 1)[1]); del argv[i]; break
    if len(argv) < 5:
        sys.exit(__doc__)
    outp = argv[0]
    h0, h1 = float(argv[1]), float(argv[2])
    ref_tag = argv[3]
    arms = [a.split('=', 1) for a in argv[4:]]
    m0, m1 = h0 + lead, h1 + lead
    if lead:
        print(f'lead {lead:g} h: wall {h0:g}Z-{h1:g}Z = MODEL hours {m0:g}-{m1:g}')

    tag = f'h{h0:g}_h{h1:g}'
    d02 = np.load(f'{REFS}/{ref_tag}_d02_{tag}.npy')
    mrms = np.load(f'{REFS}/{ref_tag}_mrms_{tag}.npy')
    # Land mask only where a terrain reference exists for this grid. The ocean
    # diagnostic domain has NO land, so masking to terrain > 20 m would leave an
    # empty set and score nothing -- there the whole interior is the scored set.
    tp = f'{REFS}/{ref_tag}_terrain.npy'
    if os.path.exists(tp):
        LAND = np.load(tp) > 20.0
        scored_set = 'interior land (terrain > 20 m)'
        print(f'scoring over LAND (terrain > 20 m) from {ref_tag}_terrain.npy')
    else:
        LAND = np.ones_like(d02, dtype=bool)
        scored_set = 'whole interior (no land in this domain)'
        print(f'no {ref_tag}_terrain.npy -- scoring the WHOLE interior '
              f'(all-ocean domain)')

    panels = [(lab, window(lab, rd, m0, m1)) for lab, rd in arms]
    panels += [('MRMS  (observed)', mrms), ('WRF d02  (1.5 km)', d02)]

    NX, NY = d02.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dmin = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    M = LAND & np.isfinite(d02) & np.isfinite(mrms) & (dmin >= BAND)

    n = len(panels)
    nrow = (n + 1) // 2
    fig, ax = plt.subplots(nrow, 2, figsize=(13.5, 4.6 * nrow), squeeze=False)
    vmax = float(np.nanpercentile(np.concatenate(
        [p[M] for _, p in panels]), 99.5))
    vmax = max(vmax, 1.0)
    for a, (lab, arr) in zip(ax.flat, panels):
        im = a.pcolormesh(np.ma.masked_invalid(arr).T, cmap='turbo',
                          vmin=0.0, vmax=vmax)
        a.set_title(f'{lab}\n{h0:g}Z-{h1:g}Z accumulation (mm)', fontsize=10)
        a.set_aspect('equal')
        for e in (BAND, NX - 1 - BAND):
            a.axvline(e, color='w', lw=0.7, ls='--', alpha=0.7)
        for e in (BAND, NY - 1 - BAND):
            a.axhline(e, color='w', lw=0.7, ls='--', alpha=0.7)
        fig.colorbar(im, ax=a, shrink=0.85)
    for a in ax.flat[n:]:
        a.axis('off')
    fig.suptitle(f'{ref_tag}, 2020-12-28 {h0:g}Z-{h1:g}Z; dashed = '
                 f'{BAND}-cell forcing-dominated band, excluded', fontsize=11)
    fig.tight_layout()
    fig.savefig(outp, dpi=130)
    print(f'wrote {outp}')

    print(f'\nScored over {scored_set}, {int(M.sum())} cells\n')
    hdr = f'{"field":26s} {"mean mm":>9s} {"vs MRMS":>9s} {"vs d02":>9s} ' \
          f'{"corr MRMS":>10s} {"corr d02":>9s}'
    print(hdr)
    print('-' * len(hdr))
    mm, md = np.nanmean(mrms[M]), np.nanmean(d02[M])
    for lab, arr in panels:
        a = arr[M]
        print(f'{lab:26s} {np.nanmean(a):9.3f} {np.nanmean(a)/mm:9.2f} '
              f'{np.nanmean(a)/md:9.2f} '
              f'{np.corrcoef(a, mrms[M])[0,1]:+10.3f} '
              f'{np.corrcoef(a, d02[M])[0,1]:+9.3f}')
    print('\nReference points, 23-h Domain A window: CONUS404 driver +0.709 vs '
          'MRMS, WRF d02 +0.910.\nThe ERF numbers from that window (Davies '
          '+0.061, NSCBC +0.487) are VOID -- they were\nmeasured before item 83, '
          'when frame indexing ignored start_datetime. Do not compare\nagainst '
          'them. Ocean h1-h2 references: MRMS vs d02 +0.708.')


if __name__ == '__main__':
    main()
