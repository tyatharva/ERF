import os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
OUT='/app/ERF/figs'; os.makedirs(OUT,exist_ok=True); MM2IN=1/25.4
NX,NY=192,96
lat=np.load('/app/ERF/scoring_ab_dav/lat.npy'); lon=np.load('/app/ERF/scoring_ab_dav/lon.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy'); land=terr>terr.min()+20
wrf=np.load('/app/ERF/wrf_d02_on_grid.npy')*MM2IN
mrms=np.load('/app/ERF/mrms_20201228_on_grid.npy')*MM2IN
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
d=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])
def dec(ax,first):
    ax.contour(lon,lat,land.astype(float),levels=[0.5],colors='k',linewidths=0.6,zorder=3)
    ax.contour(lon,lat,d.astype(float),levels=[19.5],colors='w',linewidths=1.1,zorder=4)
    ax.contour(lon,lat,d.astype(float),levels=[9.5],colors='w',linewidths=0.9,linestyles=':',zorder=4)
    ax.set_xlim(lon.min(),lon.max()); ax.set_ylim(lat.min(),lat.max()); ax.set_xlabel('lon')
    ax.set_ylabel('lat') if first else ax.set_yticklabels([])
for vmax in (12.0,3.0):
    fig,ax=plt.subplots(1,3,figsize=(15.5,4.8),constrained_layout=True)
    cm=plt.get_cmap('viridis').copy(); cm.set_over(cm(1.0)); nrm=Normalize(0,vmax)
    for k,(nm,f) in enumerate((('MRMS Pass-2 (observed)',mrms),('wrfout d02 1.5 km -> 3 km',wrf))):
        m=ax[k].pcolormesh(lon,lat,f,cmap=cm,norm=nrm,shading='nearest',zorder=1); dec(ax[k],k==0)
        ax[k].set_title(f'{nm}\nland mean {np.nanmean(f[land]):.2f} in   domain mean {np.nanmean(f):.2f} in\n'
                        f'land max {np.nanmax(f[land]):.2f} in',fontsize=9)
    fig.colorbar(m,ax=ax[:2],extend='max',shrink=0.85,label=f'24-h precip (in), 0-{vmax:g}')
    dn=TwoSlopeNorm(vmin=-1.5,vcenter=0.,vmax=1.5)
    md=ax[2].pcolormesh(lon,lat,wrf-mrms,cmap='RdBu_r',norm=dn,shading='nearest',zorder=1); dec(ax[2],False)
    ax[2].set_title('wrfout - MRMS\nland bias 1.09x, PCC 0.756\n(red = wrfout wetter)',fontsize=9)
    fig.colorbar(md,ax=ax[2],extend='both',label='inches')
    fig.suptitle('2020-12-28 00Z-23Z.  black = model land mask (coast + Channel Islands);  '
                 'solid white = d=20;  dotted white = d=10',fontsize=10)
    p=f'{OUT}/wrf_vs_mrms_0-{vmax:g}in.png'; fig.savefig(p,dpi=150); plt.close(fig); print('wrote',p)
