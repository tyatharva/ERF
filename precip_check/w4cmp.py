import yt, numpy as np, glob
yt.set_log_level(50)
print('OUR dry control (Sep-9), effect of relaxation width:')
print(' shell |  width 10 w>1  w_rms |   width 4 w>1  w_rms')
res={}
for cfg,lab in [('p1_base','w10'),('w4','w4')]:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    res[lab]=[(100*(w[d3==k]>1).mean(), np.sqrt((w[d3==k]**2).mean())) for k in range(10)]
for k in range(10):
    a,b=res['w10'][k],res['w4'][k]
    print(f'  d={k:2d} | {a[0]:12.2f}% {a[1]:6.3f} | {b[0]:12.2f}% {b[1]:6.3f}')
