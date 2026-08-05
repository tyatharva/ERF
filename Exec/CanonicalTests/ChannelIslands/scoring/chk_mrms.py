import glob, gzip, os, shutil, numpy as np, pygrib
from scipy.interpolate import griddata
lat=np.load('/app/ERF/scoring_ab_dav/lat.npy'); lon=np.load('/app/ERF/scoring_ab_dav/lon.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy')
NX,NY=lat.shape
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
d=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])
files=sorted(glob.glob('/app/ERF/mrms/*.grib2.gz'))
print(f'{len(files)} hourly files: {os.path.basename(files[0])[-21:-9]} .. {os.path.basename(files[-1])[-21:-9]}')
tot=None; nneg=0
for k,fz in enumerate(files):
    tmp='/tmp/m.grib2'
    with gzip.open(fz,'rb') as a, open(tmp,'wb') as b: shutil.copyfileobj(a,b)
    g=pygrib.open(tmp); m=g.message(1); v=m.values
    if hasattr(v,'filled'): v=v.filled(np.nan)
    v=np.asarray(v,dtype=float); v[v<0]=np.nan          # MRMS uses -1/-3 for no-coverage
    if tot is None:
        mlat,mlon=m.latlons(); mlon=np.where(mlon>180,mlon-360,mlon)
        tot=np.zeros_like(v)
    tot=np.where(np.isnan(v),tot,tot+np.nan_to_num(v))
    g.close()
print(f'native MRMS grid {tot.shape}, summed total mean over CONUS {np.nanmean(tot):.2f} mm')
sel=(mlat>=lat.min()-0.1)&(mlat<=lat.max()+0.1)&(mlon>=lon.min()-0.1)&(mlon<=lon.max()+0.1)
print(f'source points in our footprint: {int(sel.sum())}, native mean there {np.nanmean(tot[sel]):.2f} mm')
pts=np.column_stack([mlon[sel],mlat[sel]])
near=griddata(pts,tot[sel],(lon,lat),method='nearest')
lin =griddata(pts,tot[sel],(lon,lat),method='linear')
mrms=np.where(np.isnan(lin),np.nan,near)
print(f'regridded: valid {int(np.isfinite(mrms).sum())}/{mrms.size} cells; '
      f'native-vs-regridded control {np.nanmean(tot[sel]):.2f} vs {np.nanmean(mrms):.2f} mm')
land=terr>terr.min()+20
for nm,msk in (('whole domain',np.ones_like(land,bool)),('interior d>=20',d>=20),
               ('band d<10',d<10),('interior LAND',(d>=20)&land),('LAND all',land)):
    a=mrms[msk&np.isfinite(mrms)]
    if a.size: print(f'  {nm:16s} n={a.size:6d} mean {a.mean():6.2f}  p90 {np.percentile(a,90):6.2f}  max {a.max():7.2f} mm')
k=np.unravel_index(np.nanargmax(mrms),mrms.shape)
print(f'\ndomain max {np.nanmax(mrms):.1f} mm at lat {lat[k]:.3f} lon {lon[k]:.3f}, '
      f'd={int(d[k])} -> {"INSIDE relax band (d<10)" if d[k]<10 else ("BAND MARGIN 10<=d<20" if d[k]<20 else "scored interior")}')
ki=np.unravel_index(np.nanargmax(np.where(d>=20,mrms,np.nan)),mrms.shape)
print(f'interior(d>=20) max {mrms[ki]:.1f} mm at lat {lat[ki]:.3f} lon {lon[ki]:.3f}')
np.save('/app/ERF/mrms_20201228_on_grid.npy',mrms)
print('saved /app/ERF/mrms_20201228_on_grid.npy')
