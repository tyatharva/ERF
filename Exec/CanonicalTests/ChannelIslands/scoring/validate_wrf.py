import numpy as np, pyproj
from scipy.ndimage import uniform_filter
PLO=(-388229.74,-166933.25); DX=3000.0; NX,NY=192,96
P4=("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
    "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
acc=np.load('/app/ERF/wrf_d02_accum_mm.npy')
la,lo=np.load('/app/ERF/wrf_d02_latlon.npy')
tr=pyproj.Transformer.from_crs(4326,pyproj.CRS.from_proj4(P4),always_xy=True)
x,y=tr.transform(lo,la)
i=np.floor((x-PLO[0])/DX).astype(int); j=np.floor((y-PLO[1])/DX).astype(int)
ok=(i>=0)&(i<NX)&(j>=0)&(j<NY)
s=np.zeros((NX,NY)); c=np.zeros((NX,NY))
np.add.at(s,(i[ok],j[ok]),acc[ok]); np.add.at(c,(i[ok],j[ok]),1.0)
wrf=np.where(c>0,s/np.maximum(c,1),np.nan)
print(f'regrid: {int(ok.sum())} d02 cells binned; target cells filled {int((c>0).sum())}/{NX*NY}; '
      f'source per target min/med/max {int(c[c>0].min())}/{int(np.median(c[c>0]))}/{int(c.max())}')
print(f'conservation check: d02 mean over our box {acc[ok].mean():.3f} mm vs regridded mean {np.nanmean(wrf):.3f} mm')
np.save('/app/ERF/wrf_d02_on_grid.npy',wrf)

mrms=np.load('/app/ERF/mrms_20201228_on_grid.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy'); land=terr>terr.min()+20
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
d=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])

def fss(f,o,thr,w,m):
    Pf=uniform_filter((f>=thr).astype(float),size=w,mode='constant')
    Po=uniform_filter((o>=thr).astype(float),size=w,mode='constant')
    num=np.nanmean((Pf[m]-Po[m])**2); den=np.nanmean(Pf[m]**2)+np.nanmean(Po[m]**2)
    return 1-num/den if den>0 else np.nan

for nm,m in (('LAND all (n=2831)',land),('LAND interior d>=20 (n=509)',land&(d>=20))):
    m=m&np.isfinite(mrms)&np.isfinite(wrf)
    a,b=wrf[m],mrms[m]
    print(f'\n=== {nm} ===')
    print(f'  wrfout mean {a.mean():6.2f}  MRMS mean {b.mean():6.2f}  bias {a.mean()/b.mean():5.2f}x   '
          f'PCC {np.corrcoef(a,b)[0,1]:5.3f}   RMSE {np.sqrt(((a-b)**2).mean()):5.2f} mm')
    print(f'  p90 {np.percentile(a,90):6.2f} vs {np.percentile(b,90):6.2f}   max {a.max():6.2f} vs {b.max():6.2f}')
    print(f'  {"thr":>5} ' + ' '.join(f'{k:>7}' for k in ('3km','9km','15km','30km','60km')))
    for thr in (1.,5.,15.,30.):
        r=[fss(wrf,mrms,thr,w,m) for w in (1,3,5,11,21)]
        base=float((mrms[m]>=thr).mean()); useful=0.5+base/2
        print(f'  {thr:5.0f} ' + ' '.join(f'{v:7.3f}' for v in r) + f'   base {base:.3f} useful {useful:.3f}')
