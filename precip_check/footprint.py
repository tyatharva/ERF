import pygrib, numpy as np, sys
from pyproj import CRS, Transformer

AREA=[36.0,-123.25,31.25,-115.25]           # N,W,S,E  (as fed to CreateLCCMapping)
lat1,lat2,lon1,lon2=AREA[2],AREA[0],AREA[1],AREA[3]
delta=lat2-lat1; lon0=(lon1+lon2)/2; lat0=(lat1+lat2)/2
lcc=(f"+proj=lcc +lat_1={lat1+delta/6:.6f} +lat_2={lat2-delta/6:.6f} "
     f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")
to_ll=Transformer.from_crs(CRS.from_proj4(lcc), CRS.from_epsg(4326), always_xy=True)

PLO=(-195131.04,-126372.41); PHI=(188868.96,65627.59); NX,NY=128,64
dx=(PHI[0]-PLO[0])/NX; dy=(PHI[1]-PLO[1])/NY

def cell_ll(i,j):
    x=PLO[0]+(i+0.5)*dx; y=PLO[1]+(j+0.5)*dy
    lon,lat=to_ll.transform(x,y); return lat,lon

for cut in (0,20):
    pts=[cell_ll(i,j) for i in range(cut,NX-cut) for j in range(cut,NY-cut)]
    la=[p[0] for p in pts]; lo=[p[1] for p in pts]
    print(f"d>={cut:2d} footprint: lat {min(la):.3f}..{max(la):.3f}  lon {min(lo):.3f}..{max(lo):.3f}")

tag,d0,d1=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
grbs=pygrib.open(f'era5_precip_{tag}.grib'); tp={}; lats=lons=None
for g in grbs:
    if g.shortName!='tp': continue
    tp[(g.validityDate,g.validityTime)]=np.array(g.values)*1000.0
    if lats is None: lats,lons=g.latlons()
grbs.close()
sel=[k for k in sorted(tp) if (k[0]==d0 and k[1]>0) or (k[0]==d1 and k[1]==0)]
tot=sum(tp[k] for k in sel)
lons=np.where(lons>180,lons-360,lons)

for cut,label in ((0,'full ERF footprint'),(20,'clean interior d>=20')):
    pts=[cell_ll(i,j) for i in range(cut,NX-cut) for j in range(cut,NY-cut)]
    la=[p[0] for p in pts]; lo=[p[1] for p in pts]
    m=(lats>=min(la))&(lats<=max(la))&(lons>=min(lo))&(lons<=max(lo))
    v=tot[m]
    print(f"  ERA5 over {label}: n={m.sum():4d} cells  mean {v.mean():6.1f}  p95 {np.percentile(v,95):6.1f}  MAX {v.max():6.1f} mm")
