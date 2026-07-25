"""Decompose the band w signal by which wall it is adjacent to.

The relaxation zone is a ring, but the four walls sit over very different
surfaces (ocean vs the Transverse Ranges), so a ring-averaged shell profile
cannot distinguish a boundary artifact from orography.
"""
import yt, numpy as np, glob, os, sys
yt.set_log_level(50)
cfg = os.environ.get('CFG', 'ctl2')
root = '/app/ERF/bdyfix'

pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*')
             if p.split('plt')[-1].isdigit()],
            key=lambda p: int(p.split('plt')[-1]))
ds = yt.load(pl[-1])
g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
w = np.asarray(g[('boxlib', 'z_velocity')])
try:
    ter = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
except Exception:
    ter = None
nx, ny, nz = w.shape
ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
dring = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj])

walls = {
    'xlo': ii, 'xhi': nx - 1 - ii,
    'ylo': jj, 'yhi': ny - 1 - jj,
}
print(f'{cfg}: {pl[-1].split("/")[-1]}  grid {nx}x{ny}x{nz}', file=sys.stderr)
print(f'\n  {cfg}: |w|>1 fraction by wall and distance '
      f'(cells scored only where that wall is the nearest)')
print('    d |     xlo      xhi      ylo      yhi  |   ring')
print('  ----+------------------------------------+--------')
for d in range(0, 13):
    row = []
    for name, dist in walls.items():
        m = (dist == d) & (dring == d)
        m3 = np.repeat(m[:, :, None], nz, axis=2)
        row.append(100 * (np.abs(w[m3]) > 1).mean() if m3.sum() else float('nan'))
    m3 = np.repeat((dring == d)[:, :, None], nz, axis=2)
    ring = 100 * (np.abs(w[m3]) > 1).mean()
    print('  %3d | %7.2f%% %7.2f%% %7.2f%% %7.2f%% | %6.2f%%'
          % (d, row[0], row[1], row[2], row[3], ring))
m3 = np.repeat((dring >= 25)[:, :, None], nz, axis=2)
print('  interior background (d>=25): %.2f%%' % (100 * (np.abs(w[m3]) > 1).mean()))

if ter is not None:
    print('\n  mean terrain height (m) at each wall, d=0..8:')
    for name, dist in walls.items():
        hs = [float(ter[(dist == d) & (dring == d)].mean()) for d in range(9)]
        print('    %-4s ' % name + ' '.join('%6.0f' % h for h in hs))
