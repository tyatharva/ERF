#!/usr/bin/env python3
"""Hourly domain-mean precipitation RATE for all four fields, 2020-12-28.

Answers whether the Davies/NSCBC divergence is present from hour 1 or develops.
Rates are per-hour increments, so each series is differenced in its own native
accumulation: ERF rain_accum between consecutive hourly plotfiles, d02
RAINNC+bucket between consecutive stamps, MRMS the per-file hourly total.
"""
import glob, gzip, shutil, os, numpy as np, netCDF4 as nc, pyproj, pygrib, yt, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
yt.set_log_level(50)
OUT='/app/ERF/figs'; NX,NY,DX=192,96,3000.
PLO=(-388229.74,-166933.25)
PIN=("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
     "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
lat=np.load('/app/ERF/scoring_ab_dav/lat.npy'); lon=np.load('/app/ERF/scoring_ab_dav/lon.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy'); LAND=terr>terr.min()+20
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
D=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])

def erf_series(run):
    got={}
    for p in sorted(glob.glob(f'/app/ERF/{run}/plt[0-9]*')):
        try: ds=yt.load(p)
        except Exception: continue
        t=float(ds.current_time); h=int(round(t/3600.))
        if abs(t-h*3600.)>400. or h<1 or h>23: continue
        if h in got: continue
        g=ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
        got[h]=np.asarray(g[('boxlib','rain_accum')])[:,:,0]
    hs=sorted(got); out={}
    for k,h in enumerate(hs):
        prev=got[hs[k-1]] if k>0 else 0.0
        out[h]=got[h]-prev
    return out

def d02_series():
    f=nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
    ts=[b''.join(r).decode() for r in f.variables['Times'][:]]; B=float(f.BUCKET_MM)
    R=lambda k: np.asarray(f.variables['RAINNC'][k])+B*np.asarray(f.variables['I_RAINNC'][k])
    la=np.asarray(f.variables['XLAT'][0]); lo=np.asarray(f.variables['XLONG'][0])
    x,y=pyproj.Transformer.from_crs(4326,pyproj.CRS.from_proj4(PIN),always_xy=True).transform(lo,la)
    i=np.floor((x-PLO[0])/DX).astype(int); j=np.floor((y-PLO[1])/DX).astype(int)
    ok=(i>=0)&(i<NX)&(j>=0)&(j<NY)
    out={}
    for h in range(1,24):
        d=R(h)-R(h-1)
        s=np.zeros((NX,NY)); c=np.zeros((NX,NY))
        np.add.at(s,(i[ok],j[ok]),d[ok]); np.add.at(c,(i[ok],j[ok]),1.)
        out[h]=np.where(c>0,s/np.maximum(c,1),np.nan)
    return out

def mrms_series():
    out={}
    pts=None
    for fz in sorted(glob.glob('/app/ERF/mrms/*.grib2.gz')):
        h=int(os.path.basename(fz).split('-')[1][0:2])
        with gzip.open(fz,'rb') as a, open('/tmp/m.g2','wb') as b: shutil.copyfileobj(a,b)
        g=pygrib.open('/tmp/m.g2'); m=g.message(1)
        v=m.values
        v=v.filled(np.nan) if hasattr(v,'filled') else np.asarray(v,float)
        v=np.asarray(v,float); v[v<0]=np.nan
        if pts is None:
            mla,mlo=m.latlons(); mlo=np.where(mlo>180,mlo-360,mlo)
            sel=(mla>=lat.min()-.1)&(mla<=lat.max()+.1)&(mlo>=lon.min()-.1)&(mlo<=lon.max()+.1)
            pts=np.column_stack([mlo[sel],mla[sel]])
        near=griddata(pts,np.nan_to_num(v[sel]),(lon,lat),method='nearest')
        out[h]=near; g.close()
    return out

S={'Davies':erf_series('run_c404'),'NSCBC':erf_series('run_c404_nsc'),
   'd02':d02_series(),'MRMS':mrms_series()}
cols={'d02':'k','MRMS':'tab:blue','Davies':'tab:green','NSCBC':'tab:red'}
fig,ax=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True)
for c,(mn,mk) in enumerate((('full domain',np.ones_like(LAND,bool)),
                            ('interior d>=20',D>=20),('LAND',LAND))):
    for name,ser in S.items():
        hs=sorted(ser); v=[np.nanmean(ser[h][mk]) for h in hs]
        ax[c].plot(hs,v,label=name,color=cols[name],
                   lw=2.2 if name in ('d02','MRMS') else 1.6,
                   ls='-' if name in ('d02','MRMS') else '--',marker='o',ms=3)
    ax[c].axvline(18,color='0.5',ls=':'); ax[c].set_xlabel('hour (UTC) of 2020-12-28')
    ax[c].set_title(f'{mn}  [dotted = Davies FPE at 18 h]',fontsize=10)
    ax[c].grid(alpha=.3); ax[c].set_ylabel('domain-mean rate, mm/h') if c==0 else None
    ax[c].legend(fontsize=8)
p=f'{OUT}/c404_hourly_rate.png'; fig.savefig(p,dpi=150); print('wrote',p)
with open(f'{OUT}/tables_18h.txt','a') as T:
    T.write('\n\n### HOURLY DOMAIN-MEAN RATE (mm/h) ###\n')
    for mn,mk in (('full domain',np.ones_like(LAND,bool)),('interior d>=20',D>=20),('LAND',LAND)):
        T.write(f'\n-- {mn} --\n{"h":>3} ' + ' '.join(f'{k:>9}' for k in S) + '\n')
        for h in range(1,24):
            T.write(f'{h:>3} ' + ' '.join(
                f'{np.nanmean(S[k][h][mk]):9.3f}' if h in S[k] else f'{"--":>9}' for k in S) + '\n')
print('appended hourly table')
for mn,mk in (('full domain',np.ones_like(LAND,bool)),('LAND',LAND)):
    print(f'\n{mn}: hourly mean rate mm/h')
    print('  h  ' + ' '.join(f'{k:>8}' for k in S))
    for h in (1,2,3,6,9,12,15,18):
        print(f' {h:>2}  ' + ' '.join(
            f'{np.nanmean(S[k][h][mk]):8.3f}' if h in S[k] else f'{"--":>8}' for k in S))
