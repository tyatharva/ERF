"""Recover the LCC the run_a3 frames were actually built on, from the frame itself.

The .bin carries BOTH the projected grid (xvec[nx], yvec[ny]) and the geographic
coordinates (lat[nx*ny], lon[nx*ny]) of the same points, so the projection is
over-determined: transform (lon,lat) under a candidate proj4 and compare to the
stored (x,y). The true one has ~0 residual; anything else is metres of error.
"""
import numpy as np, pyproj

f = '/app/ERF/run_a3/ERA5Data_3D/ERF_IC_2023_01_09_00_00.bin'
raw = open(f, 'rb').read()
nx, ny, nz, nd = np.frombuffer(raw, dtype='<i4', count=4)
print(f'nx={nx} ny={ny} nz={nz} ndata={nd}')
o = 16
g = lambda n: np.frombuffer(raw, dtype='<f4', count=n, offset=o)
lat = g(nx*ny); o += 4*nx*ny
lon = g(nx*ny); o += 4*nx*ny
x   = g(nx);    o += 4*nx
y   = g(ny);    o += 4*ny
z   = g(nz)
print(f'x  {x.min():.0f} .. {x.max():.0f}   y {y.min():.0f} .. {y.max():.0f}')
print(f'lat {lat.min():.3f} .. {lat.max():.3f}  lon {lon.min():.3f} .. {lon.max():.3f}')
print(f'z  {z[:4]}')

def lcc(area):
    lat1, lat2 = area[2], area[0]; lon1, lon2 = area[1], area[3]
    d = lat2 - lat1
    return (f"+proj=lcc +lat_1={lat1+d/6:.6f} +lat_2={lat2-d/6:.6f} "
            f"+lat_0={(lat1+lat2)/2:.6f} +lon_0={(lon1+lon2)/2:.6f} "
            f"+datum=WGS84 +units=m +no_defs")

cands = {
    'PINNED  (36.0,-123.25,31.25,-115.25)': lcc((36.0, -123.25, 31.25, -115.25)),
    'SCORING (36.0,-125.0 ,30.75,-115.25)': lcc((36.0, -125.0,  30.75, -115.25)),
}

for order in ('lat[j*nx+i]', 'lat[i*ny+j]'):
    if order == 'lat[j*nx+i]':
        LA = lat.reshape(ny, nx).T; LO = lon.reshape(ny, nx).T   # -> (nx,ny)
    else:
        LA = lat.reshape(nx, ny);   LO = lon.reshape(nx, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    print(f'\n--- ordering {order} ---')
    for name, p4 in cands.items():
        tr = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(p4), always_xy=True)
        xt, yt = tr.transform(LO, LA)
        r = np.hypot(xt - X, yt - Y)
        print(f'  {name}: residual mean {r.mean():10.1f} m   max {r.max():10.1f} m')
