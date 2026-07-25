import yt, numpy as np, glob
yt.set_log_level(50)
for run in ['run_stress_aug','run_jan09','run_stress_sep']:
    pl=sorted(glob.glob(f'/app/ERF/{run}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')]); qc=np.asarray(g[('boxlib','qc')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    b=d3<10; m=(d3>=10)&(d3<20); i2=d3>=20
    f=lambda x,s: (np.sqrt((x[s]**2).mean()), np.abs(x[s]).max())
    print(f'--- {run} (t=24h) ---')
    for nm,s in (('band  d<10',b),('mid 10-19',m),('interior>=20',i2)):
        wr,wm=f(w,s)
        print(f'  {nm:12s} w_rms {wr:6.3f}  |w|max {wm:6.2f} m/s   qc mean {qc[s].mean()*1000:7.4f} g/kg  updraft frac(w>1) {100*(w[s]>1).mean():5.2f}%')
