import yt, numpy as np, glob, sys
yt.set_log_level(50)
run=sys.argv[1]
ds=yt.load(sorted(glob.glob(f'/app/ERF/{run}/plt*'))[-1])
g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
ra=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
z =np.asarray(g[('boxlib','z_phys')])[:,:,0]      # terrain height at k=0
nx,ny=ra.shape
# distance (in cells) from nearest lateral boundary
ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
print(f'--- {run} ---   domain {nx}x{ny}, terrain {z.min():.0f}-{z.max():.0f} m')
print('  dist-from-boundary shells (cells): mean / p99 / max rain (mm), n cells')
for lo,hi in [(0,4),(5,9),(10,14),(15,19),(20,29),(30,63)]:
    m=(d>=lo)&(d<=hi)
    if m.sum(): print(f'    {lo:2d}-{hi:2d}: {ra[m].mean():8.1f} {np.percentile(ra[m],99):8.1f} {ra[m].max():9.1f}   n={m.sum()}')
I=(d>=10)   # interior
print('  interior terrain bins: mean / p99 / max rain (mm), n')
for lo,hi in [(0,50),(50,200),(200,500),(500,1000),(1000,3000)]:
    m=I&(z>=lo)&(z<hi)
    if m.sum(): print(f'    {lo:4d}-{hi:4d} m: {ra[m].mean():8.1f} {np.percentile(ra[m],99):8.1f} {ra[m].max():9.1f}   n={m.sum()}')
print(f'  interior corr(rain, terrain) = {np.corrcoef(ra[I], z[I])[0,1]:+.3f}')
# coarse-grain interior to ~25 km (ERA5 scale) = 8x8 cells of 3 km
sub=ra[8:120,8:56]
cg=sub.reshape(sub.shape[0]//8,8,sub.shape[1]//8,8).mean(axis=(1,3))
print(f'  coarse-grained to 24 km: mean {cg.mean():.1f} p99 {np.percentile(cg,99):.1f} MAX {cg.max():.1f} mm')
print(f'  native 3 km interior:    mean {ra[I].mean():.1f} p99 {np.percentile(ra[I],99):.1f} MAX {ra[I].max():.1f} mm')
