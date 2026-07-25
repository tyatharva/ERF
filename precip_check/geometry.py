"""Trade-curve geometry: island positions in LCC meters from dem.tif.

Classifies land into mainland (components touching DEM N/E edges) and
islands, prints each island's LCC bounding box, then for each candidate
resolution reports whether a domain within the 256x128-column budget can
hold island sets with >=20-cell clearance to every boundary while keeping
the mainland outside.
"""
import numpy as np, rasterio
from pyproj import Transformer
from scipy import ndimage

# erftools-convention LCC for the current campaign box (dem_to_erf_terrain.py)
PROJ = ("+proj=lcc +lat_1=33.055453 +lat_2=34.211453 +lat_0=33.633453 "
        "+lon_0=-119.291839 +datum=WGS84 +units=m +no_defs")
# Recompute exactly as the generator does
import sys
sys.path.insert(0, '/app/ERF/Exec/CanonicalTests/ChannelIslands')
from dem_to_erf_terrain import geodesic_box, lcc_proj_string
c = geodesic_box(34.5, -117.2, 384.0, 192.0)
lat_n = max(c["NE"][0], c["NW"][0]); lat_s = min(c["SE"][0], c["SW"][0])
lon_w = min(c["NW"][1], c["SW"][1]); lon_e = max(c["NE"][1], c["SE"][1])
PROJ = lcc_proj_string(lat_s, lat_n, lon_w, lon_e)
print("LCC:", PROJ)

src = rasterio.open('/app/ERF/Exec/CanonicalTests/ChannelIslands/dem.tif')
band = src.read(1)
band = np.where(np.isfinite(band), band, 0.0)
land = band > 0.5   # meters; ocean/nodata are 0

# Downsample for component analysis (DEM is 1-arcsec-ish; keep it manageable)
step = max(1, land.shape[0] // 2000)
land_s = land[::step, ::step]
lab, nlab = ndimage.label(land_s)
sizes = ndimage.sum(land_s, lab, range(1, nlab + 1))

# Mainland: any component that touches the N or E edge of the DEM,
# or is huge (>25% of land pixels)
edge_ids = set(np.unique(lab[0, :])) | set(np.unique(lab[:, -1])) | set(np.unique(lab[-1, :]))
edge_ids.discard(0)
big_ids = {i + 1 for i, s in enumerate(sizes) if s > 0.25 * land_s.sum()}
mainland_ids = edge_ids | big_ids

tf = Transformer.from_crs("EPSG:4326", PROJ, always_xy=True)

def rc_to_lcc(rows, cols):
    xs, ys = rasterio.transform.xy(src.transform, rows * step, cols * step)
    return tf.transform(np.asarray(xs), np.asarray(ys))

print(f"\ncomponents: {nlab}, mainland ids: {len(mainland_ids)}")
islands = []
for cid in range(1, nlab + 1):
    if cid in mainland_ids: continue
    if sizes[cid - 1] < 20: continue  # skip rocks < ~20 coarse pixels
    rows, cols = np.where(lab == cid)
    x, y = rc_to_lcc(rows, cols)
    hmax = band[::step, ::step][lab == cid].max()
    islands.append((sizes[cid - 1], x.min(), x.max(), y.min(), y.max(), hmax))
islands.sort(reverse=True)
print(f"{'px':>6} {'xlo(km)':>9} {'xhi(km)':>9} {'ylo(km)':>9} {'yhi(km)':>9} {'hmax(m)':>8}")
for s, xlo, xhi, ylo, yhi, h in islands:
    print(f"{int(s):6d} {xlo/1e3:9.1f} {xhi/1e3:9.1f} {ylo/1e3:9.1f} {yhi/1e3:9.1f} {h:8.0f}")

# Mainland pixels in LCC (for min-distance checks)
mrows, mcols = np.where(np.isin(lab, list(mainland_ids)))
mx, my = rc_to_lcc(mrows, mcols)
print(f"\nmainland extent: x {mx.min()/1e3:.1f}..{mx.max()/1e3:.1f} km, "
      f"y {my.min()/1e3:.1f}..{my.max()/1e3:.1f} km, n={len(mx)}")
np.savez('/app/ERF/precip_check/geometry_lcc.npz',
         islands=np.array([(s, xlo, xhi, ylo, yhi, h) for s, xlo, xhi, ylo, yhi, h in islands]),
         mainland_x=mx, mainland_y=my)
print("saved geometry_lcc.npz")
