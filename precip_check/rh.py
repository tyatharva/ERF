import yt, numpy as np, glob, sys
yt.set_log_level(50)
def esat(T):  # Bolton, Pa
    Tc=T-273.15
    return 611.2*np.exp(17.67*Tc/(Tc+243.5))
for run in ['run_stress_aug','run_jan09','run_stress_sep']:
    pl=sorted(glob.glob(f'/app/ERF/{run}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    print(f'--- {run} ---')
    for p in [pl[1], pl[len(pl)//2], pl[-1]]:
        ds=yt.load(p); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        qv=np.asarray(g[('boxlib','qv')]); T=np.asarray(g[('boxlib','temp')]); P=np.asarray(g[('boxlib','pressure')])
        qs=0.622*esat(T)/np.maximum(P-esat(T),1.0)
        rh=qv/np.maximum(qs,1e-12)*100.0
        nx,ny,nz=rh.shape
        ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
        d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
        d3=np.repeat(d[:,:,None],nz,axis=2)
        lo=slice(0,12)   # lowest 12 levels ~ below 1 km
        b=rh[:,:,lo][d3[:,:,lo]<10]; i2=rh[:,:,lo][d3[:,:,lo]>=20]
        print(f'  t={float(ds.current_time)/3600:5.1f} h  band RH mean {b.mean():5.1f}% p99 {np.percentile(b,99):6.1f}% frac>100% {100*(b>100).mean():5.2f}%'
              f' | interior RH mean {i2.mean():5.1f}% frac>100% {100*(i2>100).mean():5.2f}%')
