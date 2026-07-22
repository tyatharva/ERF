#!/usr/bin/env python3
"""
dem_to_erf_terrain.py -- convert a GeoTIFF DEM to ERF's plain terrain format
for the Channel Islands hindcast, using the SAME Lambert Conformal Conic
projection convention as erftools' ERA5 preprocessing, so the terrain and the
ERA5 IC/BC data live in one consistent coordinate system.

The domain is anchored at its NORTHEAST corner (lat/lon) with fixed E-W / N-S
extents in km. Corner geometry uses the WGS84 geodesic (pyproj.Geod).

Modes
-----
1) Print the lat/lon corners, the DEM download box, and the erftools ERA5
   'area' line (run this BEFORE downloading anything):

     python3 dem_to_erf_terrain.py --print-bounds

2) Convert a downloaded GeoTIFF (EPSG:4326 lat/lon grid, e.g. USGS 3DEP,
   SRTM 1-arcsec, or Copernicus GLO-30) to ERF's terrain_file_name format:

     python3 dem_to_erf_terrain.py --dem dem.tif --out channel_islands_terrain.txt

   If you already ran erftools and have Output/domain_extents.txt, pass
   --xlo/--xhi/--ylo/--yhi with its prob_lo/prob_hi values so the terrain grid
   matches the run domain exactly; otherwise a symmetric box about the LCC
   origin is used.

ERF file format written (Source/ERF_ProbCommon.H, read_custom_terrain,
is_usgs = false):   nx, ny, x[0..nx-1], y[0..ny-1], then nx*ny elevations
with index = i*ny + j (x-major, contiguous in y). Coordinates in meters in
the LCC plane; ocean / nodata elevations are clamped to 0.

Projection convention (replicating erftools WriteICFromERA5Data.py):
  +proj=lcc +lat_1=<latS + range/6> +lat_2=<latN - range/6>
            +lat_0=<center lat> +lon_0=<center lon> +datum=WGS84 +units=m
VERIFY against the erftools version you run (its CreateLCCMapping) -- if it
differs, pass the exact string via --proj-str.
"""
import argparse
import sys

# ---- Domain anchor (northeast corner) and extents --------------------------
NE_LAT_DEFAULT = 34.5
NE_LON_DEFAULT = -117.2
EW_KM_DEFAULT = 384.0
NS_KM_DEFAULT = 192.0


def geodesic_box(ne_lat, ne_lon, ew_km, ns_km):
    """Corners of the domain box via WGS84 geodesic from the NE anchor."""
    from pyproj import Geod
    g = Geod(ellps="WGS84")
    # west along the northern edge; south along the eastern edge
    nw_lon, nw_lat, _ = g.fwd(ne_lon, ne_lat, 270.0, ew_km * 1e3)
    se_lon, se_lat, _ = g.fwd(ne_lon, ne_lat, 180.0, ns_km * 1e3)
    # southwest: go west along the southern edge (wider in degrees than north)
    sw_lon, sw_lat, _ = g.fwd(se_lon, se_lat, 270.0, ew_km * 1e3)
    return {
        "NE": (ne_lat, ne_lon), "NW": (nw_lat, nw_lon),
        "SE": (se_lat, se_lon), "SW": (sw_lat, sw_lon),
    }


def lcc_proj_string(lat_s, lat_n, lon_w, lon_e):
    """erftools-convention LCC for the lat/lon box."""
    rng = lat_n - lat_s
    lat_1 = lat_s + rng / 6.0
    lat_2 = lat_n - rng / 6.0
    lat0 = 0.5 * (lat_s + lat_n)
    lon0 = 0.5 * (lon_w + lon_e)
    return (f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
            f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ne-lat", type=float, default=NE_LAT_DEFAULT)
    ap.add_argument("--ne-lon", type=float, default=NE_LON_DEFAULT)
    ap.add_argument("--ew-km", type=float, default=EW_KM_DEFAULT)
    ap.add_argument("--ns-km", type=float, default=NS_KM_DEFAULT)
    ap.add_argument("--print-bounds", action="store_true",
                    help="print corners, DEM download box, ERA5 area; exit")
    ap.add_argument("--dem", help="input GeoTIFF (EPSG:4326)")
    ap.add_argument("--out", default="channel_islands_terrain.txt")
    ap.add_argument("--target-dx", type=float, default=333.0,
                    help="terrain grid spacing in meters (default 333, ~dx/3)")
    ap.add_argument("--proj-str", default=None,
                    help="override LCC proj string (use erftools' exact one)")
    ap.add_argument("--xlo", type=float, default=None, help="prob_lo x from domain_extents.txt")
    ap.add_argument("--xhi", type=float, default=None, help="prob_hi x from domain_extents.txt")
    ap.add_argument("--ylo", type=float, default=None, help="prob_lo y from domain_extents.txt")
    ap.add_argument("--yhi", type=float, default=None, help="prob_hi y from domain_extents.txt")
    args = ap.parse_args()

    corners = geodesic_box(args.ne_lat, args.ne_lon, args.ew_km, args.ns_km)
    lat_n = max(corners["NE"][0], corners["NW"][0])
    lat_s = min(corners["SE"][0], corners["SW"][0])
    lon_w = min(corners["NW"][1], corners["SW"][1])
    lon_e = max(corners["NE"][1], corners["SE"][1])

    proj = args.proj_str or lcc_proj_string(lat_s, lat_n, lon_w, lon_e)

    if args.print_bounds:
        print("Domain corners (WGS84 geodesic from NE anchor):")
        for k in ("NE", "NW", "SE", "SW"):
            print(f"  {k}: lat {corners[k][0]:9.4f}   lon {corners[k][1]:10.4f}")
        margin = 0.15  # deg, covers projection curvature + interpolation halo
        print("\nDEM download box (with 0.15 deg margin):")
        print(f"  lat: {lat_s - margin:.3f} .. {lat_n + margin:.3f}")
        print(f"  lon: {lon_w - margin:.3f} .. {lon_e + margin:.3f}")
        print("\nerftools ERA5 config 'area' line (latN, lonW, latS, lonE):")
        print(f"  area: {lat_n:.2f},{lon_w - margin:.2f},{lat_s - margin:.2f},{lon_e:.2f}")
        print(f"\nLCC projection to be used for terrain (verify vs erftools):\n  {proj}")
        return 0

    if not args.dem:
        ap.error("--dem is required unless --print-bounds")

    import numpy as np
    import rasterio
    from pyproj import Transformer

    # Terrain grid in LCC meters. Default: symmetric box about the LCC origin
    # (erftools centers its projection on the domain); override with the exact
    # prob_lo/prob_hi from Output/domain_extents.txt when available.
    xlo = args.xlo if args.xlo is not None else -args.ew_km * 1e3 / 2
    xhi = args.xhi if args.xhi is not None else args.ew_km * 1e3 / 2
    ylo = args.ylo if args.ylo is not None else -args.ns_km * 1e3 / 2
    yhi = args.yhi if args.yhi is not None else args.ns_km * 1e3 / 2

    nx = int(round((xhi - xlo) / args.target_dx)) + 1
    ny = int(round((yhi - ylo) / args.target_dx)) + 1
    xs = np.linspace(xlo, xhi, nx)
    ys = np.linspace(ylo, yhi, ny)

    # LCC (x,y) -> lat/lon, then sample the DEM
    tf = Transformer.from_crs(proj, "EPSG:4326", always_xy=True)
    X, Y = np.meshgrid(xs, ys, indexing="ij")           # shape (nx, ny)
    lons, lats = tf.transform(X.ravel(), Y.ravel())

    with rasterio.open(args.dem) as src:
        if src.crs is not None and src.crs.to_epsg() != 4326:
            sys.exit(f"DEM CRS is {src.crs}; expected EPSG:4326 lat/lon. "
                     "Reproject first: gdalwarp -t_srs EPSG:4326 in.tif out.tif")
        band = src.read(1).astype(np.float64)
        nodata = src.nodata
        inv = ~src.transform
        cols, rows = inv * (np.asarray(lons), np.asarray(lats))
        rows = np.round(rows).astype(int)
        cols = np.round(cols).astype(int)
        # Points outside the DEM are zero-filled (ocean), NOT edge-clamped --
        # clamping would smear the border column/row across the uncovered area.
        # Verify any uncovered region really is ocean before relying on this.
        oob = ((rows < 0) | (rows >= band.shape[0]) |
               (cols < 0) | (cols >= band.shape[1]))
        rows_c = np.clip(rows, 0, band.shape[0] - 1)
        cols_c = np.clip(cols, 0, band.shape[1] - 1)
        elev = band[rows_c, cols_c]
        elev[oob] = 0.0
        n_oob = int(oob.sum())
        if n_oob:
            print(f"NOTE: {n_oob}/{oob.size} grid points "
                  f"({100.0 * n_oob / oob.size:.2f}%) fall outside the DEM "
                  "and were zero-filled (assumed ocean).")

    if nodata is not None:
        elev = np.where(np.isclose(elev, nodata), 0.0, elev)
    elev = np.where(np.isfinite(elev), elev, 0.0)
    elev = np.maximum(elev, 0.0)                        # ocean / voids -> 0
    Z = elev.reshape(nx, ny)                            # index = (i, j)

    with open(args.out, "w") as f:
        f.write(f"{nx}\n{ny}\n")
        for v in xs:
            f.write(f"{v:.6f}\n")
        for v in ys:
            f.write(f"{v:.6f}\n")
        for i in range(nx):                             # ERF plain format:
            for j in range(ny):                         # index = i*ny + j
                f.write(f"{Z[i, j]:.3f}\n")

    print(f"Wrote {args.out}: nx={nx} ny={ny} "
          f"({2 + nx + ny + nx * ny} lines), "
          f"elev min/max = {Z.min():.1f}/{Z.max():.1f} m")
    print(f"LCC used: {proj}")
    if Z.max() < 100.0:
        print("WARNING: max elevation < 100 m -- wrong box or wrong DEM?")
    return 0


if __name__ == "__main__":
    sys.exit(main())
