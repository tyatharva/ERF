"""Inspect level-1 near-surface fields before the blowup."""
import sys, glob
import yt, numpy as np
yt.set_log_level(50)

run = sys.argv[1] if len(sys.argv) > 1 else '/app/ERF/run_nested'

def mosaic(ds, field, lev=1):
    gs = [g for g in ds.index.grids if g.Level == lev]
    lo = np.min([g.get_global_startindex() for g in gs], axis=0)
    hi = np.max([g.get_global_startindex() + g.ActiveDimensions for g in gs], axis=0)
    a = np.full(tuple(hi - lo), np.nan)
    for g in gs:
        s = g.get_global_startindex() - lo
        d = g.ActiveDimensions
        a[s[0]:s[0]+d[0], s[1]:s[1]+d[1], s[2]:s[2]+d[2]] = np.asarray(g['boxlib', field])
    return a, lo

for p in sorted(glob.glob(f'{run}/plt0*'), key=lambda q: int(q.split('plt')[-1].split('.')[0])):
    if '.old' in p: continue
    ds = yt.load(p)
    try:
        th, lo = mosaic(ds, 'rhotheta'); rho, _ = mosaic(ds, 'density')
        th = th / rho
        w, _ = mosaic(ds, 'z_velocity')
    except Exception as e:
        print(p, 'skip:', e); continue
    k = 1
    ths = th[:, :, k]
    i = np.unravel_index(np.nanargmin(ths), ths.shape)
    j = np.unravel_index(np.nanargmax(np.abs(w[:, :, k])), ths.shape)
    print(f"{p.split('/')[-1]}: k={k} theta min {np.nanmin(ths):7.2f} at fine{tuple(int(v)+int(l) for v,l in zip(i,lo[:2]))} "
          f"max {np.nanmax(ths):7.2f}; |w| max {np.nanmax(np.abs(w[:,:,k])):6.2f} at fine{tuple(int(v)+int(l) for v,l in zip(j,lo[:2]))}")
    # column profile at the theta-min point
    prof = th[i[0], i[1], :8]
    print('   theta(k=0..7) at min point:', ' '.join(f'{v:.1f}' for v in prof))
    rhs = rho[i[0], i[1], :4]
    print('   rho(k=0..3) at min point:  ', ' '.join(f'{v:.4f}' for v in rhs))
