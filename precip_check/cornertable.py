"""Corner pass/fail table: do the relaxation-band corners still hold the domain
maximum |w|?

The test for UPSTREAM_ISSUES #27. Before the corner blend, the four corners held
the largest |w| in the domain on BOTH domains -- 7.75 against an interior 2.09 on
128x64 (on runs scored as clean), and 18.00 on 192x96, which killed the 24-h run.
The corner term is a grad(F).(A-B) source with no physical counterpart, so a
corner that dominates the domain is a defect however stable the run looks.

PASS = corners no longer hold the domain maximum, i.e. corner |w| falls to the
same order as the interior maximum.

usage:  RUNS='run_a run_b' python3 cornertable.py
"""
import numpy as np, yt, glob, os
yt.set_log_level(50)

RUNS = os.environ.get('RUNS', 'run_192x96_24h run_192x96_fix').split()
BAND = int(os.environ.get('BAND', '10'))

print('=' * 96)
print('CORNER PASS/FAIL   max|w| at each band corner vs the interior maximum')
print(f'  band = outer {BAND} cells; corner = the {BAND}x{BAND} square where two bands overlap')
print('  PASS = corners no longer hold the domain max')
print('=' * 96)
print()
print(f'  {"run":26s} {"t(h)":>6s} {"interior":>9s} {"xlo/ylo":>9s} {"xlo/yhi":>9s} '
      f'{"xhi/yhi":>9s} {"xhi/ylo":>9s} {"max/int":>8s}  verdict')

for run in RUNS:
    ps = sorted([q for q in glob.glob(f'/app/ERF/{run}/plt[0-9]*') if q.split('plt')[-1].isdigit()],
                key=lambda q: int(q.split('plt')[-1]))
    if not ps:
        print(f'  {run:26s}  no plotfiles')
        continue
    for p in [ps[-1]]:                       # last good output
        ds = yt.load(p)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        w = np.abs(np.asarray(g[('boxlib', 'z_velocity')]))
        t = float(ds.current_time) / 3600.0
        if not np.isfinite(w).all():
            print(f'  {run:26s} {t:6.2f}   NON-FINITE')
            continue
        NX, NY = w.shape[0], w.shape[1]
        wm = w.max(axis=2)
        ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
        dring = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])
        interior = wm[dring >= BAND].max()
        # the corner SQUARES, not just the corner cell
        sq = {'xlo/ylo': (slice(0, BAND), slice(0, BAND)),
              'xlo/yhi': (slice(0, BAND), slice(NY-BAND, NY)),
              'xhi/yhi': (slice(NX-BAND, NX), slice(NY-BAND, NY)),
              'xhi/ylo': (slice(NX-BAND, NX), slice(0, BAND))}
        cv = {k: wm[s].max() for k, s in sq.items()}
        ratio = max(cv.values()) / max(interior, 1e-9)
        verdict = 'PASS' if ratio <= 1.0 else ('marginal' if ratio <= 1.5 else 'FAIL')
        print(f'  {run:26s} {t:6.2f} {interior:9.2f} ' +
              ' '.join(f'{cv[k]:9.2f}' for k in ['xlo/ylo', 'xlo/yhi', 'xhi/yhi', 'xhi/ylo']) +
              f' {ratio:8.2f}  {verdict}')

print()
print('  max/int = largest corner-square max|w| divided by the interior max|w|.')
print('  <= 1.00 means the corners no longer hold the domain maximum.')
