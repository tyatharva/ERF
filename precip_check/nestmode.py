"""Watch the growing mode at the nest rim across snapshots."""
import sys, glob
import yt, numpy as np
yt.set_log_level(50)

run = sys.argv[1] if len(sys.argv) > 1 else '/app/ERF/run_stress_sep/nest_smoke'

def mosaic(ds, field):
    g1 = [g for g in ds.index.grids if g.Level == 1]
    lo = np.min([g.get_global_startindex() for g in g1], axis=0)
    hi = np.max([g.get_global_startindex() + g.ActiveDimensions for g in g1], axis=0)
    a = np.full(tuple(hi - lo), np.nan)
    for g in g1:
        s = g.get_global_startindex() - lo
        d = g.ActiveDimensions
        a[s[0]:s[0]+d[0], s[1]:s[1]+d[1], s[2]:s[2]+d[2]] = np.asarray(g['boxlib', field])
    return a

for p in sorted(glob.glob(f'{run}/plt0*'), key=lambda q: int(q.split('plt')[-1])):
    ds = yt.load(p)
    w = mosaic(ds, 'z_velocity')
    th = mosaic(ds, 'rhotheta') / mosaic(ds, 'density')
    nx, ny, nz = w.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    d3 = np.repeat(np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])[:, :, None], nz, axis=2)
    thbar = th.mean(axis=(0, 1))          # per-level mean
    thp = th - thbar[None, None, :]
    aw = np.abs(w); iw = np.unravel_index(np.nanargmax(aw), aw.shape)
    at = np.abs(thp); it = np.unravel_index(np.nanargmax(at), at.shape)
    print(f'{p.split("/")[-1]}: |w|max {aw[iw]:7.3f} at {iw} d={d3[iw]}   '
          f"|th'|max {at[it]:7.2f} at {it} d={d3[it]}")
    for dd, lab in ((0, 'd=0'), (1, 'd=1'), (2, 'd=2'), (5, 'd=5')):
        s = d3 == dd
        print(f'    {lab}: |w|max {np.nanmax(aw[s]):8.3f}  w_rms {np.sqrt(np.nanmean(w[s]**2)):7.4f}  '
              f"|th'|max {np.nanmax(at[s]):7.2f}")
