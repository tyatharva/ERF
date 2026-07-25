"""Shell-resolved 24-h precipitation, and the spin-up / contamination split.

Two different things both look like "wrong precipitation near the inflow edge":

  spin-up       -- a DEFICIT, only on INFLOW walls, recovering inward over some
                   distance, caused by driving with hydrometeors zeroed
                   (Roberge et al. 2024, GMD 17, 1497-1510: up to 300 km).
  contamination -- an EXCESS or noise at the wall, on ALL walls regardless of
                   flow direction, caused by the boundary treatment.

Scoring per wall separates them: a signal present only where the flow enters is
spin-up, a signal present on every wall is the boundary.

env: CFGS (comma list), ROOT, ERA5 (grib tag, e.g. "jan"), D0/D1 (grib dates)
"""
import yt, numpy as np, glob, os, sys
yt.set_log_level(50)

root = os.environ.get('ROOT', '/app/ERF/bdyfix')
CFGS = os.environ.get('CFGS', 'jan_nsc,jan_ctl').split(',')
DMAX = 26

PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59); NX, NY = 128, 64


def load(cfg):
    pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    # ERF writes an extra plotfile after a vanishingly small final step (dt clipped
    # onto stop_datetime), and rain_accum comes out NaN in it while every prognostic
    # field is clean.  Take the last output whose rain_accum is finite.
    for cand in reversed(pl):
        ds = yt.load(cand)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        if not np.isnan(ra).any():
            return ra, np.asarray(g[('boxlib', 'z_velocity')]), cand, float(ds.current_time)
    raise RuntimeError('no plotfile with finite rain_accum')


for cfg in CFGS:
    try:
        ra, w, last, tcur = load(cfg)
    except Exception as e:
        print(f'\n{cfg}: LOAD FAILED ({e})'); continue
    nx, ny = ra.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    d = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])

    print(f'\n=== {cfg}  ({last.split("/")[-1]}, t = {tcur/3600:.1f} h) ===')
    print(f'  domain-mean 24-h precip: {ra.mean():.2f} mm   max {ra.max():.1f} mm')
    print(f'  {"cutoff":>8s} {"mean":>8s} {"max":>8s}   (cells at d >= cutoff)')
    for c in (0, 1, 2, 3, 5, 10, 15, 20):
        m = d >= c
        print(f'  {c:>8d} {ra[m].mean():8.2f} {ra[m].max():8.1f}')

    # Per-wall profile: precip and w, scored only where that wall is nearest.
    walls = {'xlo': ii, 'xhi': nx-1-ii, 'ylo': jj, 'yhi': ny-1-jj}
    print(f'\n  precip (mm) by wall and distance  [deficit only on inflow walls = spin-up;'
          f'\n                                     signal on all walls = boundary]')
    print('    d |     xlo      xhi      ylo      yhi |    ring')
    print('  ----+------------------------------------+--------')
    for dd in range(0, DMAX):
        row = []
        for name, dist in walls.items():
            m = (dist == dd) & (d == dd)
            row.append(ra[m].mean() if m.sum() else float('nan'))
        m = (d == dd)
        if dd < 20 or dd % 5 == 0:
            print('  %3d | %7.2f %7.2f %7.2f %7.2f | %7.2f'
                  % (dd, row[0], row[1], row[2], row[3], ra[m].mean()))
    deep = ra[d >= 20].mean()
    print(f'  deep interior (d>=20) mean: {deep:.2f} mm  '
          f'-- per-wall recovery is measured against this')

    # band w on a wet day, same metric as the dry control
    nz = w.shape[2]
    ocean = ((ii == d) | (jj == d))
    print('\n  |w|>1 fraction (ocean walls) : ' + ' '.join(
        '%d:%.1f%%' % (dd, 100*(np.abs(w[np.repeat(((d == dd) & ocean)[:, :, None], nz, axis=2)]) > 1).mean())
        for dd in (0, 1, 2, 3, 6, 10, 13, 16)))
    m3 = np.repeat((d >= 25)[:, :, None], nz, axis=2)
    print('  interior background (d>=25) : %.2f%%' % (100*(np.abs(w[m3]) > 1).mean()))
