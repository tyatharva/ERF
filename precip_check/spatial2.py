import yt, numpy as np, glob, sys
yt.set_log_level(50)
run=sys.argv[1]
ds=yt.load(sorted(glob.glob(f'/app/ERF/{run}/plt*'))[-1])
g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
ra=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
z =np.asarray(g[('boxlib','z_phys')])[:,:,0]
nx,ny=ra.shape
ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
print(f'=== {run} ===')
for cut in (10,15,20,25):
    m=d>=cut
    print(f'  interior d>={cut:2d}: n={m.sum():5d} mean {ra[m].mean():7.1f}  p95 {np.percentile(ra[m],95):7.1f}  p99 {np.percentile(ra[m],99):7.1f}  max {ra[m].max():8.1f} mm')
m=d>=20
print(f'  --- clean interior (d>=20) terrain bins ---')
for lo,hi in [(0,50),(50,200),(200,500),(500,3000)]:
    s=m&(z>=lo)&(z<hi)
    if s.sum(): print(f'    {lo:4d}-{hi:4d} m: n={s.sum():5d} mean {ra[s].mean():7.1f} p99 {np.percentile(ra[s],99):7.1f} max {ra[s].max():8.1f}')
print(f'  corr(rain,terrain) d>=20: {np.corrcoef(ra[m],z[m])[0,1]:+.3f}')
# coarse-grain the CLEAN interior to ERA5 scale (8 cells = 24 km)
s=ra[24:104,24:40]
cg=s.reshape(s.shape[0]//8,8,s.shape[1]//8,8).mean(axis=(1,3))
print(f'  clean interior coarse-grained to 24 km: mean {cg.mean():.1f} p99 {np.percentile(cg,99):.1f} MAX {cg.max():.1f} mm')
