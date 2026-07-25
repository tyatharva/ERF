"""Shell profile of band contamination vs distance from the lateral boundary.

Scores the w>1 m/s fraction and w_rms in shells d = 0..DMAX cells inward, for
the control and mass-consistent boundary runs, plus the domain mass drift.
"""
import yt, numpy as np, glob, sys
yt.set_log_level(50)

DMAX = 26
import os
CFGS = os.environ.get('CFGS','ctl2,wfc').split(',')
root = '/app/ERF/bdyfix'

res, mass = {}, {}
for cfg in CFGS:
    pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*')
                 if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    ds = yt.load(pl[-1])
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w = np.asarray(g[('boxlib', 'z_velocity')])
    nx, ny, nz = w.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    d2 = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj])
    d3 = np.repeat(d2[:, :, None], nz, axis=2)
    res[cfg] = [(100 * (np.abs(w[d3 == k]) > 1).mean(),
                 float(np.sqrt((w[d3 == k] ** 2).mean()))) for k in range(DMAX)]
    # interior background: everything at least 30 cells in
    bg = np.abs(w[d3 >= 30])
    res[cfg].append((100 * (bg > 1).mean(), float(np.sqrt((w[d3 >= 30] ** 2).mean()))))
    print(f'{cfg}: {pl[-1].split("/")[-1]}  grid {nx}x{ny}x{nz}', file=sys.stderr)

print(f'\n  shell |   {CFGS[0]:>10s} w>1   w_rms | {CFGS[1]:>16s} w>1   w_rms |  ratio')
print('  ------+-------------------------+-----------------------------+-------')
for k in range(DMAX):
    c, m = res[CFGS[0]][k], res[CFGS[1]][k]
    r = (m[0] / c[0]) if c[0] > 1e-9 else float('nan')
    print(f'  d={k:2d}  | {c[0]:15.2f}% {c[1]:7.3f} | {m[0]:19.2f}% {m[1]:7.3f} | {r:5.2f}')
c, m = res[CFGS[0]][DMAX], res[CFGS[1]][DMAX]
print('  ------+-------------------------+-----------------------------+-------')
print(f'  d>=30 | {c[0]:15.2f}% {c[1]:7.3f} | {m[0]:19.2f}% {m[1]:7.3f} |   (interior background)')
