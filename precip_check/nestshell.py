"""Per-shell w diagnostic relative to the NEST boundary (level 1).

Reads level-1 FABs directly (yt covering_grid mis-places level-1 z for
anisotropic ref_ratio, but raw grid data + integer indices are exact).
Usage: python3 nestshell.py <run_dir> [plt_prefix]
"""
import sys, glob
import yt, numpy as np
yt.set_log_level(50)

run = sys.argv[1] if len(sys.argv) > 1 else '/app/ERF/run_stress_sep/nest_smoke'
pl = sorted(glob.glob(f'{run}/plt*'), key=lambda p: int(p.split('plt')[-1]))
ds = yt.load(pl[-1])
grids1 = [g for g in ds.index.grids if g.Level == 1]
if not grids1:
    sys.exit('no level-1 grids in ' + pl[-1])

# Assemble level-1 box (single box expected)
lo = np.min([g.get_global_startindex() for g in grids1], axis=0)
hi = np.max([g.get_global_startindex() + g.ActiveDimensions for g in grids1], axis=0)
nx, ny, nz = hi - lo
w = np.full((nx, ny, nz), np.nan, dtype=np.float64)
qv = np.full((nx, ny, nz), np.nan, dtype=np.float64)
for g in grids1:
    s = g.get_global_startindex() - lo
    d = g.ActiveDimensions
    w[s[0]:s[0]+d[0], s[1]:s[1]+d[1], s[2]:s[2]+d[2]] = np.asarray(g['boxlib','z_velocity'])
    try:
        qv[s[0]:s[0]+d[0], s[1]:s[1]+d[1], s[2]:s[2]+d[2]] = np.asarray(g['boxlib','qv'])
    except Exception:
        pass
assert np.isfinite(w).all(), 'level-1 mosaic has holes'

ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
d2 = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])
d3 = np.repeat(d2[:, :, None], nz, axis=2)
print(f'level-1 box {nx}x{ny}x{nz}, file {pl[-1]}')
print(' shell |  w>1 frac   w_rms   |w|max')
for k in range(12):
    s = d3 == k
    print(f'  d={k:2d} | {100*(w[s]>1).mean():8.2f}%  {np.sqrt((w[s]**2).mean()):6.3f}  {np.abs(w[s]).max():6.2f}')
s = d3 >= 20
print(f'  d>=20| {100*(w[s]>1).mean():8.2f}%  {np.sqrt((w[s]**2).mean()):6.3f}  {np.abs(w[s]).max():6.2f}')

# Level-0 reference: same diagnostic relative to the DOMAIN boundary
g0 = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
w0 = np.asarray(g0[('boxlib','z_velocity')])
n0x, n0y, n0z = w0.shape
i0, j0 = np.meshgrid(np.arange(n0x), np.arange(n0y), indexing='ij')
d0 = np.repeat(np.minimum.reduce([i0, j0, n0x-1-i0, n0y-1-j0])[:, :, None], n0z, axis=2)
b = d0 < 10; i2 = d0 >= 20
print(f'\nlevel-0: band d<10 w>1 {100*(w0[b]>1).mean():5.2f}% w_rms {np.sqrt((w0[b]**2).mean()):6.3f} | '
      f'interior d>=20 w>1 {100*(w0[i2]>1).mean():5.2f}% w_rms {np.sqrt((w0[i2]**2).mean()):6.3f}')
