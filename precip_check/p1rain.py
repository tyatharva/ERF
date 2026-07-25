import yt, numpy as np, glob
yt.set_log_level(50)
print(f'{"config":10s} {"blend":6s} {"bdy rho/w":10s} {"dom-mean":>9s} {"band d<10":>10s} {"d>=20":>8s} {"band share":>11s}')
for cfg,bl,rw in [('p1_base','ON','no'),('p1_rho','ON','rho'),('p1_rhow','ON','rho+w'),
                  ('nb_base','OFF','no'),('nb_rhow','OFF','rho+w')]:
    pl=sorted(glob.glob(f'/app/ERF/run_stress_sep/{cfg}/plt*'), key=lambda p:int(p.split('plt')[-1]))
    ds=yt.load(pl[-1]); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
    nx,ny=ra.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
    tot=ra.sum()
    print(f'{cfg:10s} {bl:6s} {rw:10s} {ra.mean():9.2f} {ra[d<10].mean():10.2f} {ra[d>=20].mean():8.2f} {100*ra[d<10].sum()/max(tot,1e-9):10.1f}%')
