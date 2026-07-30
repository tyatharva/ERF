#!/usr/bin/env python3
"""Full 0-18 h analysis and visualisation of the CONUS404-driven arms.

Arms: Davies (18 h), NSCBC sigma=0.03 (18 h), wrfout d02 (00Z-18Z), MRMS (00Z-18Z).
All on the 192x96 ERF grid, all over the identical window.

Scoring against d02 measures AGREEMENT WITH ANOTHER MODEL, not skill: d02 carries a
1.09x wet bias concentrated on orographic peaks and +68% small-scale variance
against MRMS over land (39a). MRMS-over-land is the observational check, and MRMS
is gap-filled offshore, which is why the observational tables are land-only.
"""
import os, numpy as np, yt, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.ndimage import uniform_filter, label
yt.set_log_level(50)

OUT = '/app/ERF/figs'; os.makedirs(OUT, exist_ok=True)
TBL = open(f'{OUT}/tables_18h.txt', 'w')
MM2IN = 1/25.4
NX, NY, DX = 192, 96, 3.0
WINS = [(3,1),(9,3),(15,5),(30,11),(60,21)]
THRS = [1.,5.,15.,30.]

lat = np.load('/app/ERF/scoring_ab_dav/lat.npy'); lon = np.load('/app/ERF/scoring_ab_dav/lon.npy')
terr = np.load('/app/ERF/scoring_ab_dav/terrain.npy'); LAND = terr > terr.min()+20
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
D = np.minimum.reduce([ii, jj, NX-1-ii, NY-1-jj])
FETCH = np.minimum(ii, jj); OUTF = np.minimum(NX-1-ii, NY-1-jj)

def plt_field(p):
    ds = yt.load(p); g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    return np.asarray(g[('boxlib','rain_accum')])[:,:,0]

F = {'d02': np.load('/app/ERF/ref18_d02.npy'),
     'MRMS': np.load('/app/ERF/ref18_mrms.npy'),
     'Davies': plt_field('/app/ERF/run_c404/plt69876'),
     'NSCBC': plt_field('/app/ERF/run_c404_nsc/plt57508')}
ARMS = ['Davies','NSCBC']; REFS = ['d02','MRMS']
def W(s): TBL.write(s+'\n')

MASKS = [('full domain', np.ones_like(LAND,bool)),
         ('interior d>=20', D>=20), ('LAND', LAND)]

def fss(f,o,thr,w,m):
    Pf=uniform_filter((f>=thr).astype(float),size=w,mode='constant')
    Po=uniform_filter((o>=thr).astype(float),size=w,mode='constant')
    n=np.nanmean((Pf[m]-Po[m])**2); d=np.nanmean(Pf[m]**2)+np.nanmean(Po[m]**2)
    return 1-n/d if d>0 else np.nan
def binz(f,m,r): return f>=float(np.nanquantile(f[m],1-r))

W('='*78); W('0-18 h CONUS404-driven comparison. d02 = model reference (agreement,')
W('not skill; 1.09x wet on peaks, +68% small-scale variance). MRMS land-only.'); W('='*78)

# ---- bulk stats ----
W('\n### BULK STATISTICS (mm) ###')
for ref in REFS:
    for mn, mk in MASKS:
        if ref=='MRMS' and mn!='LAND': continue
        m = mk & np.isfinite(F[ref])
        W(f'\n-- vs {ref}, {mn} (n={int(m.sum())}) --')
        W(f'{"arm":10s} {"mean":>8} {"bias":>7} {"PCC":>7} {"RMSE":>8} {"MAE":>8}')
        W(f'{ref:10s} {F[ref][m].mean():8.2f} {"--":>7} {"--":>7} {"--":>8} {"--":>8}')
        for a in ARMS:
            x,y = F[a][m], F[ref][m]
            W(f'{a:10s} {x.mean():8.2f} {x.mean()/y.mean():7.2f} '
              f'{np.corrcoef(x,y)[0,1]:7.3f} {np.sqrt(((x-y)**2).mean()):8.2f} '
              f'{np.abs(x-y).mean():8.2f}')

# ---- quantiles ----
W('\n\n### QUANTILES (mm) ###')
for mn, mk in MASKS:
    m = mk & np.isfinite(F['d02']) & np.isfinite(F['MRMS'])
    W(f'\n-- {mn} (n={int(m.sum())}) --')
    W(f'{"field":10s} ' + ' '.join(f'{q:>8}' for q in ('q10','q25','q50','q75','q90','q99','max')))
    for k in ['d02','MRMS']+ARMS:
        v=F[k][m]
        W(f'{k:10s} ' + ' '.join(f'{np.percentile(v,q):8.2f}' for q in (10,25,50,75,90,99))
          + f' {v.max():8.2f}')

# ---- wet-area fraction ----
W('\n\n### WET-AREA FRACTION ###')
for mn, mk in MASKS:
    m = mk & np.isfinite(F['d02']) & np.isfinite(F['MRMS'])
    W(f'\n-- {mn} --');  W(f'{"field":10s} ' + ' '.join(f'{">=%gmm"%t:>9}' for t in THRS))
    for k in ['d02','MRMS']+ARMS:
        W(f'{k:10s} ' + ' '.join(f'{(F[k][m]>=t).mean():9.3f}' for t in THRS))

# ---- FSS tables ----
for ref in REFS:
    for mn, mk in MASKS:
        if ref=='MRMS' and mn!='LAND': continue
        m = mk & np.isfinite(F[ref])
        W(f'\n\n### FSS vs {ref}, {mn} (n={int(m.sum())}) ###')
        for thr in THRS:
            base=float((F[ref][m]>=thr).mean())
            if base<0.02 or base>0.98:
                W(f'  thr {thr:g}: degenerate base rate {base:.3f} -- not scoreable'); continue
            W(f'  thr {thr:g}  base {base:.3f}  useful {0.5+base/2:.3f}')
            for a in ARMS:
                W(f'    {a:8s} fixed ' + ' '.join(f'{fss(F[a],F[ref],thr,w,m):7.3f}' for _,w in WINS))
            ob=binz(F[ref],m,base)
            for a in ARMS:
                fb=binz(F[a],m,base)
                W(f'    {a:8s} PM    ' + ' '.join(
                    f'{fss(fb.astype(float),ob.astype(float),0.5,w,m):7.3f}' for _,w in WINS))

# ---- fetch-binned ----
W('\n\n### FETCH-BINNED BIAS (ratio to reference) ###')
for ref in REFS:
    for cn, C in (('INFLOW walls (xlo,ylo)',FETCH),('OUTFLOW walls (xhi,yhi)',OUTF)):
        W(f'\n-- vs {ref}, from {cn} --')
        W(f'{"cells":>9} {"km":>9} ' + ' '.join(f'{a:>9}' for a in ARMS))
        for lo,hi in ((0,0),(1,2),(3,5),(6,9),(10,14),(15,19),(20,29),(30,44),(45,95)):
            m=(C>=lo)&(C<=hi)&np.isfinite(F[ref])
            if ref=='MRMS': m=m&LAND
            if m.sum()<20: continue
            r=F[ref][m].mean()
            W(f'{lo:>4}-{hi:<4} {3*lo:>4}-{3*hi:<4} ' +
              ' '.join(f'{F[a][m].mean()/r:9.2f}' for a in ARMS))

# ---- islands ----
W('\n\n### ISLANDS (all interior land components) ###')
lb,n = label(LAND, structure=np.ones((3,3)))
edge = set(lb[0,:])|set(lb[-1,:])|set(lb[:,0])|set(lb[:,-1])
W(f'{"lat":>7} {"lon":>9} {"n":>4} ' + ' '.join(f'{k:>9}' for k in ['d02','MRMS']+ARMS))
for i in range(1,n+1):
    if i in edge: continue
    s = lb==i
    if s.sum()<4: continue
    W(f'{lat[s].mean():7.2f} {lon[s].mean():9.2f} {int(s.sum()):4d} ' +
      ' '.join(f'{np.nanmean(F[k][s]):9.2f}' for k in ['d02','MRMS']+ARMS))
TBL.close()

# ================= FIGURES =================
def dec(ax, first):
    ax.contour(lon,lat,LAND.astype(float),levels=[0.5],colors='k',linewidths=0.6,zorder=3)
    ax.contour(lon,lat,D.astype(float),levels=[19.5],colors='w',linewidths=1.1,zorder=4)
    ax.contour(lon,lat,D.astype(float),levels=[9.5],colors='w',linewidths=0.9,linestyles=':',zorder=4)
    ax.set_xlim(lon.min(),lon.max()); ax.set_ylim(lat.min(),lat.max()); ax.set_xlabel('lon')
    ax.set_ylabel('lat') if first else ax.set_yticklabels([])

ORDER=['d02','MRMS','Davies','NSCBC']
for vmax in (12.,3.):
    fig,ax=plt.subplots(1,4,figsize=(17.5,4.8),constrained_layout=True)
    cm=plt.get_cmap('viridis').copy(); cm.set_over(cm(1.0)); nm=Normalize(0,vmax)
    for k,name in enumerate(ORDER):
        f=F[name]*MM2IN
        m=ax[k].pcolormesh(lon,lat,f,cmap=cm,norm=nm,shading='nearest',zorder=1); dec(ax[k],k==0)
        ax[k].set_title(f'{name}\ndomain {np.nanmean(f):.2f} in  interior {np.nanmean(f[D>=20]):.2f} in\n'
                        f'land {np.nanmean(f[LAND]):.2f} in  max {np.nanmax(f):.2f} in',fontsize=9)
    fig.colorbar(m,ax=ax,extend='max',shrink=0.85,label=f'0-18 h precip (in), 0-{vmax:g}')
    fig.suptitle('2020-12-28 00Z-18Z, CONUS404-driven.  black = model land mask; '
                 'solid white = d=20; dotted white = d=10 band edge',fontsize=10)
    p=f'{OUT}/c404_18h_fields_0-{vmax:g}in.png'; fig.savefig(p,dpi=150); plt.close(fig); print('wrote',p)

fig,ax=plt.subplots(2,2,figsize=(11,9.4),constrained_layout=True)
dn=TwoSlopeNorm(vmin=-2.,vcenter=0.,vmax=2.)
for r,ref in enumerate(REFS):
    for c,a in enumerate(ARMS):
        d=(F[a]-F[ref])*MM2IN
        if ref=='MRMS': d=np.where(LAND,d,np.nan)
        m=ax[r,c].pcolormesh(lon,lat,d,cmap='RdBu_r',norm=dn,shading='nearest',zorder=1)
        dec(ax[r,c],c==0)
        ax[r,c].set_title(f'{a} - {ref}' + ('  (land only)' if ref=='MRMS' else ''),fontsize=10)
fig.colorbar(m,ax=ax,extend='both',shrink=0.8,label='inches (red = model wetter)')
p=f'{OUT}/c404_18h_differences.png'; fig.savefig(p,dpi=150); plt.close(fig); print('wrote',p)

fig,ax=plt.subplots(1,2,figsize=(15,5.2),constrained_layout=True)
cols={'d02':'k','MRMS':'tab:blue','Davies':'tab:green','NSCBC':'tab:red'}
for name in ORDER:
    prof=[np.nanmean(F[name][D==d])*MM2IN for d in range(48)]
    ax[0].semilogy(np.arange(48)*DX,np.maximum(prof,1e-4),label=name,color=cols[name],
                   lw=2.2 if name in REFS else 1.6, ls='-' if name in REFS else '--')
ax[0].axvline(10*DX,color='0.5',ls=':'); ax[0].axvline(20*DX,color='0.5',ls='-')
ax[0].set_xlabel('distance inward from the NEAREST wall, km  [dotted d=10, solid d=20]')
ax[0].set_ylabel('0-18 h precip, inches'); ax[0].legend(fontsize=9); ax[0].grid(alpha=.3,which='both')
ax[0].set_title('Wall-normal profile: the band, all arms',fontsize=10)
sel=(jj>=20)&(jj<NY-20)
for name in ORDER:
    pr=np.array([np.nanmean(F[name][i,:][sel[i,:]]) for i in range(NX)])*MM2IN
    ax[1].semilogy((np.arange(NX)+.5)*DX,np.maximum(pr,1e-4),label=name,color=cols[name],
                   lw=2.2 if name in REFS else 1.6, ls='-' if name in REFS else '--')
for v,s in ((10*DX,':'),(20*DX,'-'),((NX-20)*DX,'-'),((NX-10)*DX,':')):
    ax[1].axvline(v,color='0.5',ls=s)
ax[1].set_xlabel('inflow (xlo) -> outflow (xhi), km'); ax[1].legend(fontsize=9)
ax[1].grid(alpha=.3,which='both'); ax[1].set_title('Cross-section along the flow axis',fontsize=10)
p=f'{OUT}/c404_18h_boundary.png'; fig.savefig(p,dpi=150); plt.close(fig); print('wrote',p)

def psd(x):
    c=np.nan_to_num(x)[20:NX-20,20:NY-20]; c=c-c.mean()
    Fq=np.abs(np.fft.rfft2(c*np.hanning(c.shape[0])[:,None]*np.hanning(c.shape[1])[None,:]))**2
    kx=np.fft.fftfreq(c.shape[0],DX)[:,None]; ky=np.fft.rfftfreq(c.shape[1],DX)[None,:]
    k=np.sqrt(kx**2+ky**2); b=np.linspace(0,k.max(),40); idx=np.digitize(k.ravel(),b)
    P=np.array([Fq.ravel()[idx==i].mean() if (idx==i).any() else np.nan for i in range(1,len(b))])
    return 1./(0.5*(b[1:]+b[:-1])), P
fig,ax=plt.subplots(figsize=(8.4,5.6),constrained_layout=True)
for name in ORDER:
    lam,P=psd(F[name]); ax.loglog(lam,P,label=name,color=cols[name],
                                  lw=2.2 if name in REFS else 1.6, ls='-' if name in REFS else '--')
ax.axvspan(8,19,color='0.85',zorder=0,label='8-19 km band')
ax.set_xlabel('wavelength, km'); ax.set_ylabel('radial power spectral density')
ax.set_title('0-18 h precipitation PSD, interior core',fontsize=10)
ax.legend(fontsize=9); ax.grid(alpha=.3,which='both'); ax.invert_xaxis()
p=f'{OUT}/c404_18h_spectra.png'; fig.savefig(p,dpi=150); plt.close(fig); print('wrote',p)
print(f'tables -> {OUT}/tables_18h.txt')
