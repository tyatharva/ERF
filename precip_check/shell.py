import yt, numpy as np, glob
yt.set_log_level(50)
res={}
for cfg in ['p1_base','p1_rho','p1_rhow']:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    res[cfg]=[(100*(w[d3==k]>1).mean(), np.sqrt((w[d3==k]**2).mean())) for k in range(12)]
print(' shell |   base w>1  w_rms |    rho w>1  w_rms |   rho+w w>1  w_rms')
for k in range(12):
    b,r,rw=res['p1_base'][k],res['p1_rho'][k],res['p1_rhow'][k]
    print(f'  d={k:2d} | {b[0]:8.2f}% {b[1]:6.3f} | {r[0]:8.2f}% {r[1]:6.3f} | {rw[0]:9.2f}% {rw[1]:6.3f}')
