"""Interior background |w|>1 vs time, for the acoustic-vs-gravity attribution.

d>=25 is 75 km from the nearest wall at 3 km.  Sound covers that in ~220 s;
the dominant gravity-wave modes (c* ~ 25-30 m/s) need ~2500-3000 s.
"""
import yt, numpy as np, glob, os, sys
yt.set_log_level(50)
root = '/app/ERF/bdyfix'
CFGS = os.environ.get('CFGS', 'lgctl,lgnorel').split(',')
DEEP = int(os.environ.get('DEEP', '25'))

print('\n  interior background |w|>1 (d>=%d) and mid-band peak vs time' % DEEP)
print('   cfg      step      time(s)   bg(d>=%d)   band peak (d)' % DEEP)
print('  ---------------------------------------------------------')
for cfg in CFGS:
    pl = sorted([p for p in glob.glob(f'{root}/{cfg}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
                key=lambda p: int(p.split('plt')[-1]))
    for p in pl:
        ds = yt.load(p)
        g  = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        w  = np.asarray(g[('boxlib', 'z_velocity')])
        nx, ny, nz = w.shape
        ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        dring  = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])
        ocean  = ((ii == dring) | (jj == dring))
        m3  = np.repeat((dring >= DEEP)[:, :, None], nz, axis=2)
        bg  = 100*(np.abs(w[m3]) > 1).mean()
        prof = []
        for d in range(20):
            m = np.repeat(((dring == d) & ocean)[:, :, None], nz, axis=2)
            prof.append(100*(np.abs(w[m]) > 1).mean() if m.sum() else 0.0)
        dpk = int(np.argmax(prof))
        print('  %-8s %6d %11.1f %10.2f%% %9.2f%% (d=%d)'
              % (cfg, int(p.split('plt')[-1]), float(ds.current_time), bg, prof[dpk], dpk))
