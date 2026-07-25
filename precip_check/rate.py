import yt, numpy as np, glob, sys
yt.set_log_level(50)
run=sys.argv[1]
pl=sorted(glob.glob(f'/app/ERF/{run}/plt*'), key=lambda p:int(p.split('plt')[-1]))
prev_band=prev_int=0.0; prev_t=0.0
print(f'--- {run}: precip RATE (mm/h) by region, per output hour ---')
print('   t(h)   band(d<10)   mid(10-19)   interior(d>=20)')
for p in pl[::3]:
    ds=yt.load(p); g=ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
    nx,ny=ra.shape
    ii,jj=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij')
    d=np.minimum.reduce([ii,jj,nx-1-ii,ny-1-jj])
    t=float(ds.current_time)/3600.0
    b=ra[d<10].mean(); m=ra[(d>=10)&(d<20)].mean(); i2=ra[d>=20].mean()
    if t>prev_t:
        dt=t-prev_t
        print(f'  {t:5.1f}   {(b-prev_band)/dt:10.2f}   {(m-prev_mid)/dt:10.2f}   {(i2-prev_int)/dt:12.2f}')
    prev_band, prev_mid, prev_int, prev_t = b, m, i2, t
