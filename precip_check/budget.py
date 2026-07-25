import yt, numpy as np, glob, sys
yt.set_log_level(50)
for run in ['run_stress_aug','run_stress_sep','run_jan09']:
    ds=yt.load(sorted(glob.glob(f'/app/ERF/{run}/plt*'))[-1])
    g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
    nx,ny=ra.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
    tot=ra.mean(); band=ra[d<10].sum(); allsum=ra.sum()
    print(f'{run:16s} domain-mean {tot:7.1f} mm | outer-10-cell band holds {100*band/allsum:5.1f}% of all precip | d>=20 mean {ra[d>=20].mean():6.1f}')
