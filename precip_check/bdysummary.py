"""Consolidated boundary-artifact summary.

Scores the |w|>1 fraction by distance restricted to the two OCEAN walls
(xlo, ylo -- flat 12 m terrain), which isolates the relaxation-zone artifact
from orographic vertical motion at the two land walls.
"""
import yt, numpy as np, glob, os, sys
yt.set_log_level(50)
root = '/app/ERF/bdyfix'
CFGS = os.environ.get('CFGS', 'ctl2,wfc,w15,w20,nf50').split(',')
DMAX = 26

prof, bg = {}, {}
for cfg in CFGS:
    pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*')
                 if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    ds = yt.load(pl[-1])
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w = np.asarray(g[('boxlib', 'z_velocity')])
    nx, ny, nz = w.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    dring = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj])
    ocean = ((ii == dring) | (jj == dring))          # nearest wall is xlo or ylo
    row = []
    for d in range(DMAX):
        m = (dring == d) & ocean
        m3 = np.repeat(m[:, :, None], nz, axis=2)
        row.append(100 * (np.abs(w[m3]) > 1).mean() if m3.sum() else float('nan'))
    prof[cfg] = row
    m3 = np.repeat((dring >= 25)[:, :, None], nz, axis=2)
    bg[cfg] = 100 * (np.abs(w[m3]) > 1).mean()
    print(f'{cfg}: {pl[-1].split("/")[-1]}', file=sys.stderr)

hdr = '   d  |' + ''.join(f' {c:>8s}' for c in CFGS)
print('\n  ocean-wall (xlo/ylo, flat 12 m terrain) |w|>1 fraction, by distance')
print(hdr)
print('  -----+' + '-' * (9 * len(CFGS)))
for d in range(DMAX):
    print(f'  {d:3d}  |' + ''.join(f' {prof[c][d]:7.2f}%' for c in CFGS))
print('  -----+' + '-' * (9 * len(CFGS)))
print('  bg   |' + ''.join(f' {bg[c]:7.2f}%' for c in CFGS) + '   (d>=25, all walls)')
