import numpy as np
NX,NY=192,96
wrf=np.load('/app/ERF/wrf_d02_on_grid.npy'); mrms=np.load('/app/ERF/mrms_20201228_on_grid.npy')
terr=np.load('/app/ERF/scoring_ab_dav/terrain.npy'); land=terr>terr.min()+20
ii,jj=np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
d=np.minimum.reduce([ii,jj,NX-1-ii,NY-1-jj])
def spec(x,mask):
    x=np.where(mask,np.nan_to_num(x),0.0); core=x[20:NX-20,20:NY-20]; core=core-core.mean()
    F=np.abs(np.fft.rfft2(core*np.hanning(core.shape[0])[:,None]*np.hanning(core.shape[1])[None,:]))**2
    kx=np.fft.fftfreq(core.shape[0],3.0)[:,None]; ky=np.fft.rfftfreq(core.shape[1],3.0)[None,:]
    k=np.sqrt(kx**2+ky**2); bins=np.linspace(0,k.max(),40); idx=np.digitize(k.ravel(),bins)
    P=np.array([F.ravel()[idx==i].mean() if (idx==i).any() else np.nan for i in range(1,len(bins))])
    return 1.0/(0.5*(bins[1:]+bins[:-1])),P
v=np.isfinite(mrms)
lam,Pw=spec(wrf,v); _,Pm=spec(mrms,v)
b=(lam>=8)&(lam<=19)
print(f'interior-core spectral ratio wrfout/MRMS 8-19 km: {np.nanmean(Pw[b]/Pm[b]):.3f}')
print('  (caveat: this core includes ocean, where MRMS is gap-filled -- the very thing under test)')
m=land&np.isfinite(mrms)
print('\nland quantiles (mm):  q ' + ' '.join(f'{q:>7}' for q in (10,25,50,75,90,99)))
print('  wrfout              ' + ' '.join(f'{np.percentile(wrf[m],q):7.2f}' for q in (10,25,50,75,90,99)))
print('  MRMS                ' + ' '.join(f'{np.percentile(mrms[m],q):7.2f}' for q in (10,25,50,75,90,99)))
dif=wrf[m]-mrms[m]
print(f'\nland difference: mean {dif.mean():+.2f}  median {np.median(dif):+.2f}  '
      f'|d|<5mm in {100*np.mean(np.abs(dif)<5):.0f}% of cells, |d|<10mm in {100*np.mean(np.abs(dif)<10):.0f}%')
# where do the big disagreements sit?
big=m&(np.abs(wrf-mrms)>15)
print(f'cells disagreeing by >15 mm: {int(big.sum())} of {int(m.sum())} ({100*big.sum()/m.sum():.1f}%)')
if big.sum():
    print(f'  of those, wrfout wetter in {int((wrf[big]>mrms[big]).sum())}, drier in {int((wrf[big]<mrms[big]).sum())}')
    print(f'  their wall-distance d: min {int(d[big].min())} median {int(np.median(d[big]))} max {int(d[big].max())}')
