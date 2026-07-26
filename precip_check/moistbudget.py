"""Moisture budget across the d=20 contour, and the spin-up-overshoot test.

Turns "Davies loses moisture to spurious band precipitation that NSCBC delivers
inward instead" into a measurement.  For the region d >= 20:

    IN  - P  - dStorage  =  E   (surface evaporation, not diagnosed here)

so the budget cannot be closed absolutely, but E is nearly common to both runs
(same SST, same surface scheme, similar winds), so it cancels in the DIFFERENCE.
The hypothesis predicts:  d(IN) ~ d(P) + d(dStorage) between the two runs.

Second test, spin-up overshoot vs redistribution.  With 10 hydrometeors zeroed at
inflow, air enters condensate-free and the model must spin condensate up from
scratch, which can OVERSHOOT downstream rather than only deficit at the edge.
That excess would depend on distance from the INFLOW walls (xlo, ylo here).
Redistribution would not.  So bin the NSCBC-minus-Davies precipitation excess by
distance from the inflow corner, stratified by terrain.
"""
import yt, numpy as np, glob, os

yt.set_log_level(50)
root = '/app/ERF/bdyfix'
CFGS = os.environ.get('CFGS', 'jan_nsc,jan_ctl').split(',')
CUT  = int(os.environ.get('CUT', '20'))
STRIDE = int(os.environ.get('STRIDE', '2'))          # sample every N hours
QS = ['qv', 'qc', 'qrain', 'qgraup', 'qsnow']
RHO_W = 1000.0                                        # kg/m3, for kg -> mm


def plist(cfg):
    return sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                  key=lambda p: int(p.split('plt')[-1]))


def fields(ds):
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    rho = np.asarray(g[('boxlib', 'density')])
    q = np.zeros_like(rho)
    for f in QS:
        try:
            q += np.asarray(g[('boxlib', f)])
        except Exception:
            pass
    u = np.asarray(g[('boxlib', 'x_velocity')])
    v = np.asarray(g[('boxlib', 'y_velocity')])
    z = np.asarray(g[('boxlib', 'z_phys')])
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    return rho, q, u, v, z, ra


res = {}
for cfg in CFGS:
    pl = plist(cfg)
    # drop a trailing plotfile whose rain_accum is NaN (clipped final step)
    ds = yt.load(pl[-1])
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    if np.isnan(np.asarray(g[('boxlib', 'rain_accum')])).any():
        pl = pl[:-1]

    ds0 = yt.load(pl[0])
    dx = float(ds0.domain_width[0]) / ds0.domain_dimensions[0]
    dy = float(ds0.domain_width[1]) / ds0.domain_dimensions[1]
    rho, q, u, v, z, ra0 = fields(ds0)
    nx, ny, nz = rho.shape
    # layer thickness from cell-centre heights (identical grid in both runs, so any
    # systematic error in dz cancels in the difference between them)
    dz = np.empty_like(z)
    dz[:, :, 1:-1] = 0.5 * (z[:, :, 2:] - z[:, :, :-2])
    dz[:, :, 0] = z[:, :, 1] - z[:, :, 0]
    dz[:, :, -1] = z[:, :, -1] - z[:, :, -2]

    lo, hix, hiy = CUT, nx - 1 - CUT, ny - 1 - CUT
    area_int = (hix - lo + 1) * (hiy - lo + 1) * dx * dy

    def storage(rho, q):
        return float((rho[lo:hix+1, lo:hiy+1, :] * q[lo:hix+1, lo:hiy+1, :]
                      * dz[lo:hix+1, lo:hiy+1, :]).sum() * dx * dy)

    S0 = storage(rho, q)

    # inward total-water flux across the four sides of the d=CUT contour
    times, flux = [], []
    for p in pl[::STRIDE]:
        ds = yt.load(p)
        rho, q, u, v, _, _ = fields(ds)
        rq = rho * q
        f  = float(( rq[lo,  lo:hiy+1, :] * u[lo,  lo:hiy+1, :] * dz[lo,  lo:hiy+1, :]).sum() * dy)
        f += float((-rq[hix, lo:hiy+1, :] * u[hix, lo:hiy+1, :] * dz[hix, lo:hiy+1, :]).sum() * dy)
        f += float(( rq[lo:hix+1, lo,  :] * v[lo:hix+1, lo,  :] * dz[lo:hix+1, lo,  :]).sum() * dx)
        f += float((-rq[lo:hix+1, hiy, :] * v[lo:hix+1, hiy, :] * dz[lo:hix+1, hiy, :]).sum() * dx)
        times.append(float(ds.current_time)); flux.append(f)

    IN = float(np.trapezoid(flux, times)) if hasattr(np, 'trapezoid') else float(np.trapz(flux, times))
    rho, q, u, v, _, raN = fields(yt.load(pl[-1]))
    S1 = storage(rho, q)
    P_mm = float(raN[lo:hix+1, lo:hiy+1].mean())
    P_kg = P_mm / 1000.0 * RHO_W * area_int

    res[cfg] = dict(IN=IN, dS=S1 - S0, P=P_kg, P_mm=P_mm, area=area_int, ra=raN,
                    nx=nx, ny=ny)
    tomm = lambda x: x / RHO_W / area_int * 1000.0
    print(f'\n=== {cfg}: total-water budget for d >= {CUT} '
          f'({hix-lo+1}x{hiy-lo+1} cells, {area_int/1e9:.1f} x 10^9 m2) ===')
    print(f'  inward flux across contour   IN = {tomm(IN):9.2f} mm   ({IN:.3e} kg)')
    print(f'  storage change               dS = {tomm(S1-S0):9.2f} mm')
    print(f'  precipitation                P  = {P_mm:9.2f} mm')
    print(f'  residual (IN - dS - P) = E     = {tomm(IN-(S1-S0)-P_kg):9.2f} mm '
          f'  [surface evaporation, not diagnosed]')

if len(CFGS) == 2:
    a, b = CFGS
    A, B = res[a], res[b]
    tomm = lambda x: x / RHO_W / A['area'] * 1000.0
    print(f'\n=== difference, {a} - {b} (evaporation largely cancels) ===')
    print(f'  d(IN)       = {tomm(A["IN"]-B["IN"]):8.2f} mm')
    print(f'  d(P) + d(dS)= {A["P_mm"]-B["P_mm"] + tomm(A["dS"]-B["dS"]):8.2f} mm'
          f'   [ d(P) = {A["P_mm"]-B["P_mm"]:.2f}, d(dS) = {tomm(A["dS"]-B["dS"]):.2f} ]')
    print('  -> hypothesis predicts these two match: extra moisture delivered inward'
          '\n     accounts for the extra interior precipitation.')

    # spin-up overshoot vs redistribution
    ra_a, ra_b = A['ra'], B['ra']
    nx, ny = A['nx'], A['ny']
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    d = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])
    d_in = np.minimum(ii, jj)                      # distance from the INFLOW walls
    try:
        ter = np.asarray(yt.load(plist(a)[0]).covering_grid(
            0, yt.load(plist(a)[0]).domain_left_edge,
            yt.load(plist(a)[0]).domain_dimensions)[('boxlib', 'z_phys')])[:, :, 0]
    except Exception:
        ter = np.zeros_like(ra_a)

    print(f'\n=== excess ({a} - {b}) vs distance from the INFLOW walls, d >= {CUT} ===')
    print('   d_inflow    n   flat<100m   100-400m    all')
    for b0, b1 in [(20, 30), (30, 40), (40, 55), (55, 75), (75, 200)]:
        m = (d >= CUT) & (d_in >= b0) & (d_in < b1)
        if not m.sum():
            continue
        mf = m & (ter < 100); mt = m & (ter >= 100)
        ef = (ra_a-ra_b)[mf].mean() if mf.sum() else float('nan')
        et = (ra_a-ra_b)[mt].mean() if mt.sum() else float('nan')
        print(f'   {b0:3d}-{b1:<4d} {m.sum():5d} {ef:10.2f} {et:10.2f} {(ra_a-ra_b)[m].mean():8.2f}')
    print('  -> falling with d_inflow = spin-up overshoot; flat = redistribution')
