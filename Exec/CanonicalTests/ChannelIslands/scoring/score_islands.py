import numpy as np
from scipy.ndimage import label
import sys
OUT=(sys.argv[1] if len(sys.argv)>1 else '/app/ERF/scoring_jan9').rstrip('/')+'/'
erf=np.load(OUT+'erf_mm.npy'); mrms=np.load(OUT+'mrms_mm.npy'); era5=np.load(OUT+'era5_mm.npy')
terr=np.load(OUT+'terrain.npy'); lat=np.load(OUT+'lat.npy'); lon=np.load(OUT+'lon.npy')
print('z_phys k=0 over ocean (min):', terr.min())
land = terr > terr.min() + 20.0
lab, nl = label(land)
edge = set(lab[0,:]) | set(lab[-1,:]) | set(lab[:,0]) | set(lab[:,-1])
for L in range(1, nl+1):
    m = lab == L
    if L in edge or m.sum() < 2: continue
    print(f'island n={m.sum():3d} lat {lat[m].mean():.2f} lon {lon[m].mean():.2f} zmax {terr[m].max():4.0f}m'
          f' | ERF {np.nanmean(erf[m]):6.1f}  MRMS {np.nanmean(mrms[m]):5.1f}  ERA5 {np.nanmean(era5[m]):5.1f} mm')
