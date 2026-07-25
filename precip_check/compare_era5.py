import pygrib, numpy as np, sys, glob
tag=sys.argv[1]
f=f'era5_precip_{tag}.grib'
grbs=pygrib.open(f)
tp={}
for g in grbs:
    if g.shortName!='tp': continue
    key=(g.validityDate,g.validityTime)
    tp[key]=np.array(g.values)*1000.0   # m -> mm
grbs.close()
keys=sorted(tp)
# 24-h window: sum the 24 hourly accumulations ending 01Z..00Z next day
d0=int(sys.argv[2]); d1=int(sys.argv[3])
sel=[k for k in keys if (k[0]==d0 and k[1]>0) or (k[0]==d1 and k[1]==0)]
tot=sum(tp[k] for k in sel)
print(f'{tag}: {len(sel)} hourly fields summed ({sel[0]} .. {sel[-1]}), grid {tot.shape}')
print(f'  ERA5 24-h precip over domain: mean {tot.mean():.1f}  p95 {np.percentile(tot,95):.1f}  p99 {np.percentile(tot,99):.1f}  MAX {tot.max():.1f} mm')
ny,nx=tot.shape
c=tot[3:ny-3,3:nx-3]
print(f'  ERA5 inner subregion:         mean {c.mean():.1f}  p95 {np.percentile(c,95):.1f}  MAX {c.max():.1f} mm')
