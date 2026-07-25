import yt, numpy as np, glob
yt.set_log_level(50)
print(' shell | blend-off base w>1 w_rms | blend-off rho+w w>1 w_rms')
res={}
for cfg in ['nb_base','nb_rhow']:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')]); r=np.asarray(g[('boxlib','density')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    res[cfg]=([(100*(w[d3==k]>1).mean(), np.sqrt((w[d3==k]**2).mean())) for k in range(6)], r)
for k in range(6):
    b,r=res['nb_base'][0][k],res['nb_rhow'][0][k]
    print(f'  d={k:2d} | {b[0]:14.2f}% {b[1]:6.3f} | {r[0]:15.2f}% {r[1]:6.3f}')
print('max|drho| between them:', np.abs(res['nb_rhow'][1]-res['nb_base'][1]).max())
