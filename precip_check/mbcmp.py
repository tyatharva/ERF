import yt, numpy as np, glob
yt.set_log_level(50)
res={}
for cfg in ['p1_base','nb_base','mb_mom']:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    res[cfg]=[(100*(w[d3==k]>1).mean(), np.sqrt((w[d3==k]**2).mean())) for k in range(12)]
print(' shell | blend-on(old) w>1 w_rms | blend-off w>1 w_rms | blend-mom w>1  w_rms')
for k in range(12):
    b,n,m=res['p1_base'][k],res['nb_base'][k],res['mb_mom'][k]
    print(f'  d={k:2d} | {b[0]:10.2f}% {b[1]:6.3f} | {n[0]:9.2f}% {n[1]:6.3f} | {m[0]:9.2f}% {m[1]:6.3f}')
