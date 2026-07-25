import yt, numpy as np, glob, sys
yt.set_log_level(50)
print(f'{"config":8s} {"band w>1":>9s} {"mid w>1":>9s} {"int w>1":>9s} {"band w_rms":>11s} {"int w_rms":>10s} {"band |w|max":>12s}')
for cfg in ['p1_base','p1_rho','p1_rhow']:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    w=np.asarray(g[('boxlib','z_velocity')])
    nx,ny,nz=w.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d3=np.repeat(np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])[:,:,None],nz,axis=2)
    b=d3<10; m=(d3>=10)&(d3<20); i2=d3>=20
    print(f'{cfg.replace("p1_",""):8s} {100*(w[b]>1).mean():8.2f}% {100*(w[m]>1).mean():8.2f}% {100*(w[i2]>1).mean():8.2f}% '
          f'{np.sqrt((w[b]**2).mean()):11.4f} {np.sqrt((w[i2]**2).mean()):10.4f} {np.abs(w[b]).max():12.2f}')
