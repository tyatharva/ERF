"""Vertical structure of the near-wall w signal, control vs mass-consistent."""
import yt, numpy as np, glob, sys
yt.set_log_level(50)
root = '/app/ERF/bdyfix'

def load(cfg):
    pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*')
                 if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    ds = yt.load(pl[-1])
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    return (np.asarray(g[('boxlib', 'z_velocity')]),
            np.asarray(g[('boxlib', 'density')]))

w_c, r_c = load('ctl2')
w_m, r_m = load('mc')
nx, ny, nz = w_c.shape
ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
d2 = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj])

print('\n  w_rms(k) in the first three shells, and rho there')
print('    k |   ctl d<3   mc d<3 |  ctl d>=25  mc d>=25 |  rho ctl   rho mc')
print('  ----+--------------------+----------------------+------------------')
m_band = d2 < 3
m_int = d2 >= 25
for k in range(nz):
    a = float(np.sqrt((w_c[m_band, k] ** 2).mean()))
    b = float(np.sqrt((w_m[m_band, k] ** 2).mean()))
    c = float(np.sqrt((w_c[m_int, k] ** 2).mean()))
    d = float(np.sqrt((w_m[m_int, k] ** 2).mean()))
    rc = float(r_c[m_band, k].mean())
    rm = float(r_m[m_band, k].mean())
    print(f'  {k:3d} | {a:8.3f} {b:8.3f} | {c:9.3f} {d:8.3f} | {rc:8.4f} {rm:8.4f}')

print('\n  band-mean density (d<3): ctl %.5f  mc %.5f  ratio %.4f'
      % (r_c[m_band].mean(), r_m[m_band].mean(), r_m[m_band].mean() / r_c[m_band].mean()))
print('  domain-mean density    : ctl %.5f  mc %.5f  ratio %.4f'
      % (r_c.mean(), r_m.mean(), r_m.mean() / r_c.mean()))
