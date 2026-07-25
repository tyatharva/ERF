"""Step-2 kill condition: Helmholtz decomposition of the relaxation forcing.

C_raw = F*(A-B) is the momentum forcing the Davies relaxation adds (dumped by
erf.realbdy_dump_relax_step).  The proposed fix replaces it with its
divergence-free part

    C_proj = C_raw - grad(phi),   lap(phi) = div(C_raw)

so the spurious mass source grad(F).(A-B) disappears by construction.  The
question this script answers is whether anything USEFUL survives that: if
C_raw is mostly a gradient, projecting it away also projects away the nudging
and the band stops constraining the solution.

BCs on phi, per the design: phi = 0 outside the band (where C_raw = 0, so the
correction must vanish -- this is exact, not an approximation), and homogeneous
Neumann at the domain wall so the wall-normal flux is left to the wall flux
correction.  Mixed BCs make the operator non-singular, so there is no
compatibility constant to dispose of.

Discrete operator matches the continuity equation's face-flux divergence.
"""
import yt, numpy as np, glob, os, sys
yt.set_log_level(50)

root  = '/app/ERF/bdyfix'
CFG   = os.environ.get('CFG', 'dump10')
WIDTH = int(os.environ.get('WIDTH', '10'))
TOL   = 1e-10


def cg(apply, b, mask, tol=TOL, maxit=5000):
    """CG on the SPD operator -apply, restricted to mask."""
    x = np.zeros_like(b)
    r = b - (-apply(x))
    r = np.where(mask, r, 0.0)
    p = r.copy()
    rs = float((r * r).sum())
    rs0 = rs
    if rs0 == 0.0:
        return x, 0.0, 0
    for it in range(maxit):
        Ap = np.where(mask, -apply(p), 0.0)
        den = float((p * Ap).sum())
        if den == 0.0:
            break
        a = rs / den
        x += a * p
        r -= a * Ap
        rs_new = float((r * r).sum())
        if rs_new / rs0 < tol * tol:
            rs = rs_new
            break
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x, np.sqrt(rs / rs0), it + 1


ds = yt.load(f'{root}/{CFG}/relaxdump')
g  = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
cxlo = np.asarray(g[('boxlib', 'Cx_lo')], dtype=np.float64)
cxhi = np.asarray(g[('boxlib', 'Cx_hi')], dtype=np.float64)
cylo = np.asarray(g[('boxlib', 'Cy_lo')], dtype=np.float64)
cyhi = np.asarray(g[('boxlib', 'Cy_hi')], dtype=np.float64)
nx, ny, nz = cxlo.shape
dx = float(ds.domain_width[0]) / nx
dy = float(ds.domain_width[1]) / ny

# Reconstruct the staggered face arrays exactly.
Cx = np.zeros((nx + 1, ny, nz)); Cx[:nx] = cxlo; Cx[nx] = cxhi[nx - 1]
Cy = np.zeros((nx, ny + 1, nz)); Cy[:, :ny] = cylo; Cy[:, ny] = cyhi[:, ny - 1]
assert np.allclose(Cx[1:nx], cxhi[:nx - 1]), 'face reconstruction mismatch (x)'
assert np.allclose(Cy[:, 1:ny], cyhi[:, :ny - 1]), 'face reconstruction mismatch (y)'

ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
band = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj]) < WIDTH


def grad(p):
    gx = np.zeros((nx + 1, ny)); gx[1:nx] = (p[1:] - p[:-1]) / dx
    gy = np.zeros((nx, ny + 1)); gy[:, 1:ny] = (p[:, 1:] - p[:, :-1]) / dy
    return gx, gy                     # gx[0]=gx[nx]=0 : Neumann at the wall


def lap(p):
    gx, gy = grad(np.where(band, p, 0.0))
    return (gx[1:] - gx[:-1]) / dx + (gy[:, 1:] - gy[:, :-1]) / dy


print(f'\n  {CFG}: relaxation forcing Helmholtz decomposition '
      f'(grid {nx}x{ny}x{nz}, dx {dx/1000:.1f} km, band width {WIDTH})')
print('   k   ||C_raw||    ||C_proj||/||C_raw||   retained   div reduction   CG')
print('  ------------------------------------------------------------------------')

tot_raw = tot_proj = tot_dot = 0.0
rows = []
for k in range(nz):
    cx, cy = Cx[:, :, k], Cy[:, :, k]
    div = ((cx[1:] - cx[:-1]) / dx + (cy[:, 1:] - cy[:, :-1]) / dy)
    b = np.where(band, div, 0.0)
    phi, res, nit = cg(lap, -b, band)   # cg solves -lap(phi)=rhs, we want lap(phi)=div
    gx, gy = grad(np.where(band, phi, 0.0))
    px, py = cx - gx, cy - gy

    n_raw = np.sqrt((cx ** 2).sum() + (cy ** 2).sum())
    n_prj = np.sqrt((px ** 2).sum() + (py ** 2).sum())
    dot = (cx * px).sum() + (cy * py).sum()

    dnew = ((px[1:] - px[:-1]) / dx + (py[:, 1:] - py[:, :-1]) / dy)
    d0 = np.sqrt((np.where(band, div, 0.0) ** 2).sum())
    d1 = np.sqrt((np.where(band, dnew, 0.0) ** 2).sum())

    tot_raw += n_raw ** 2; tot_proj += n_prj ** 2; tot_dot += dot
    rows.append((k, n_raw, n_prj / n_raw if n_raw else 0.0,
                 dot / n_raw ** 2 if n_raw else 0.0, d1 / d0 if d0 else 0.0, res, nit))

for k, n_raw, ratio, ret, dred, res, nit in rows:
    if k % 4 == 0 or k == nz - 1:
        print('  %3d  %10.3e      %8.4f          %8.4f     %9.2e   %4d (%.0e)'
              % (k, n_raw, ratio, ret, dred, nit, res))

print('  ------------------------------------------------------------------------')
print('  ALL  %10.3e      %8.4f          %8.4f'
      % (np.sqrt(tot_raw), np.sqrt(tot_proj / tot_raw), tot_dot / tot_raw))
print('\n  ||C_proj||/||C_raw||  = fraction of the correction that SURVIVES projection')
print('  retained              = <C_proj,C_raw>/||C_raw||^2, the fraction of the')
print('                          original nudging kept IN ITS ORIGINAL DIRECTION')
print('  div reduction         = ||div C_proj|| / ||div C_raw|| over the band')
