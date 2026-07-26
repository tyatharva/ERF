"""Stability probe scoring for the 192x96 domain: is the lateral relaxation band
generating breaking mountain waves over the terrain it now contains?

WHY. Moving the NE pin 34.2 -> 34.7 took Mt. Baldy out of the north relaxation
band but pulled far more terrain in behind it: cells above 1000 m inside the
10-cell band went 7 -> 322 on the north wall and 4 -> 90 on the east wall, with
the maximum rising 1483 -> 2095 m. The 34.5 pin was rejected for exactly this
configuration, with a measured failure at 4.3 h.

The 34.5 failure signature was a grid-scale w dipole in the band over terrain,
so this scores |w| and |v| by LOCATION -- band vs interior, and by wall -- and
asks specifically whether the extremes sit on the Alamo/Liebre block at
34.676/-118.952, the highest point now inside the band.

env: RUN (default run_192x96)
"""
import numpy as np, yt, glob, os, re, sys
from pyproj import CRS, Transformer
yt.set_log_level(50)

RUN = os.environ.get('RUN', '/app/ERF/run_192x96')
LOG = os.environ.get('LOG', '')
PLO = (-388229.74, -166933.25); PHI = (187770.26, 121066.75); NX, NY = 192, 96
BAND = 10                                   # erf.real_width
LCC = ('+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 '
       '+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs')
to_ll = Transformer.from_crs(CRS.from_proj4(LCC), CRS.from_epsg(4326), always_xy=True)
ALAMO = (34.676, -118.952)                  # highest terrain inside the band

dx = (PHI[0]-PLO[0])/NX; dy = (PHI[1]-PLO[1])/NY
xc = PLO[0]+(np.arange(NX)+0.5)*dx; yc = PLO[1]+(np.arange(NY)+0.5)*dy
LON, LAT = to_ll.transform(*np.meshgrid(xc, yc, indexing='ij'))
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
dring = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])
INBAND = dring < BAND
WALL = np.choose(np.argmin(np.stack([ii, jj, NX-1-ii, NY-1-jj]), axis=0),
                 np.array(['xlo', 'ylo', 'xhi', 'yhi'])[:, None, None])

plots = sorted([p for p in glob.glob(f'{RUN}/plt[0-9]*') if p.split('plt')[-1].isdigit()],
               key=lambda p: int(p.split('plt')[-1]))
if not plots:
    sys.exit(f'no plotfiles in {RUN}')

print('='*100)
print(f'BAND STABILITY PROBE   {RUN}   {len(plots)} plotfiles   band = outer {BAND} cells')
print(f'  watching the Alamo/Liebre block at {ALAMO[0]}, {ALAMO[1]} (2095 m, 1 cell from the yhi wall)')
print('='*100)
print()
print('  t(h)   max|w|  where            max|v|  where            band max|w|  wall  ter(m)  d_alamo')

for p in plots:
    ds = yt.load(p)
    t_h = float(ds.current_time)/3600.0
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w = np.asarray(g[('boxlib', 'z_velocity')])
    v = np.asarray(g[('boxlib', 'y_velocity')])
    ter = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
    if not np.isfinite(w).all() or not np.isfinite(v).all():
        nb = np.argwhere(~np.isfinite(w))
        print(f'  {t_h:5.2f}   NON-FINITE w at {len(nb)} cells; first at '
              f'(i,j,k)={tuple(nb[0])} lat {LAT[nb[0][0],nb[0][1]]:.3f} '
              f'lon {LON[nb[0][0],nb[0][1]]:.3f} inband={INBAND[nb[0][0],nb[0][1]]}')
        break

    aw = np.abs(w).max(axis=2); av = np.abs(v).max(axis=2)
    kw = np.unravel_index(aw.argmax(), aw.shape); kv = np.unravel_index(av.argmax(), av.shape)
    bw = np.where(INBAND, aw, -1.0); kb = np.unravel_index(bw.argmax(), bw.shape)
    # great-circle-ish distance in km from the Alamo/Liebre point
    dl = np.hypot((LAT[kb]-ALAMO[0])*111.0, (LON[kb]-ALAMO[1])*111.0*np.cos(np.radians(ALAMO[0])))
    print(f'  {t_h:5.2f}  {aw[kw]:7.2f}  {LAT[kw]:6.3f},{LON[kw]:8.3f}'
          f'{"B" if INBAND[kw] else "I"}  {av[kv]:7.2f}  {LAT[kv]:6.3f},{LON[kv]:8.3f}'
          f'{"B" if INBAND[kv] else "I"}   {bw[kb]:8.2f}  {str(WALL[kb]):4s} {ter[kb]:6.0f}  {dl:6.1f} km')

print()
print('  B = in the relaxation band, I = interior. A band max|w| that grows while the')
print('  interior does not, localised on high terrain, is the 34.5 failure signature.')

# interior vs band ratio over time, and the Alamo neighbourhood specifically
print()
print('  t(h)   interior max|w|   band max|w|   ratio   |   max|w| within 30 km of Alamo/Liebre')
near = (np.hypot((LAT-ALAMO[0])*111.0, (LON-ALAMO[1])*111.0*np.cos(np.radians(ALAMO[0]))) < 30.0)
for p in plots:
    ds = yt.load(p); t_h = float(ds.current_time)/3600.0
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w = np.asarray(g[('boxlib', 'z_velocity')])
    if not np.isfinite(w).all():
        break
    aw = np.abs(w).max(axis=2)
    wi = aw[~INBAND].max(); wb = aw[INBAND].max()
    print(f'  {t_h:5.2f}   {wi:13.2f}   {wb:11.2f}   {wb/max(wi,1e-9):5.2f}   |   {aw[near].max():8.2f}')

if LOG and os.path.exists(LOG):
    dts = []
    for line in open(LOG, errors='ignore'):
        m = re.match(r'Coarse STEP (\d+) ends. TIME = ([\d.]+) DT = ([\d.]+)', line)
        if m:
            dts.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    if dts:
        d = np.array([x[2] for x in dts]); t = np.array([x[1] for x in dts])/3600.0
        print()
        print('  dt trajectory (a collapsing dt is the model fighting an instability)')
        for lo in np.arange(0, t.max()+1e-9, 0.5):
            m = (t >= lo) & (t < lo+0.5)
            if m.sum():
                print(f'    t {lo:4.1f}-{lo+0.5:4.1f} h   dt min {d[m].min():.4f}  '
                      f'median {np.median(d[m]):.4f}  max {d[m].max():.4f}   n={m.sum()}')
        print(f'    dt at the end: {d[-1]:.4f} s (started {d[0]:.4f} s)')
    nw = sum(1 for l in open(LOG, errors='ignore') if 'damping' in l.lower())
    nt = sum(1 for l in open(LOG, errors='ignore') if 'low temp' in l.lower() or 'Tmin' in l)
    print()
    print(f'  w-damping mentions in log: {nw}    low-temperature warnings: {nt}')
    print('  (UPSTREAM_ISSUES #13: at marginal cfl this fork is not bit-reproducible in SP.')
    print('   A clean run is necessary but NOT sufficient -- require zero damping events.)')
