import yt, numpy as np, glob
yt.set_log_level(50)
ds=yt.load(sorted(glob.glob('/app/ERF/run_jan09/plt*'))[-1])
g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
z=np.asarray(g[('boxlib','z_phys')])[:,:,0]
nx,ny=z.shape
ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
print(f'terrain {z.min():.0f}-{z.max():.0f} m over {nx}x{ny}')
for cut in (0,5,10,15,20,25):
    m=d>=cut
    hi=(z>200)&m
    print(f'  d>={cut:2d}: cells {m.sum():5d}, of which terrain>200 m: {hi.sum():4d} ({100*hi.sum()/m.sum():.1f}%), max terrain {z[m].max():6.0f} m')
k=np.argmax(z); i,j=np.unravel_index(k,z.shape)
print(f'  highest terrain {z[i,j]:.0f} m at (i={i}, j={j}) -> distance from boundary = {d[i,j]} cells')
