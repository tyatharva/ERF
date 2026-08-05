import numpy as np, netCDF4 as nc
f=nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
print('dims:', {k:len(v) for k,v in f.dimensions.items() if k in
      ('Time','bottom_top','south_north','west_east','bottom_top_stag','soil_layers_stag')})
for a in ('DX','DY','MAP_PROJ','TRUELAT1','TRUELAT2','STAND_LON','CEN_LAT','CEN_LON','MOAD_CEN_LAT'):
    if hasattr(f,a): print(f'  {a} = {getattr(f,a)}')
t=f.variables['Times'][:]
ts=[b''.join(r).decode() for r in t]
print(f'times: {len(ts)}  {ts[0]} .. {ts[-1]}')
la=f.variables['XLAT'][0]; lo=f.variables['XLONG'][0]
print(f'd02 lat {la.min():.3f}..{la.max():.3f}   lon {lo.min():.3f}..{lo.max():.3f}')
print('corners (SW,SE,NW,NE):',
      [(round(float(la[j,i]),3),round(float(lo[j,i]),3)) for j,i in
       ((0,0),(0,-1),(-1,0),(-1,-1))])
OUR=dict(lat=(32.064,34.703), lon=(-123.470,-117.217))
inside=((la>=OUR['lat'][0])&(la<=OUR['lat'][1])&(lo>=OUR['lon'][0])&(lo<=OUR['lon'][1]))
print(f'\nour box lat {OUR["lat"]} lon {OUR["lon"]}')
print(f'  d02 cells inside our box: {int(inside.sum())}')
print(f'  our box fully inside d02? '
      f'lat {la.min()<=OUR["lat"][0] and la.max()>=OUR["lat"][1]}, '
      f'lon {lo.min()<=OUR["lon"][0] and lo.max()>=OUR["lon"][1]}')
v3=[k for k,v in f.variables.items() if v.ndim==4]
print('\n4-D (3-D+time) variables:', sorted(v3)[:24])
