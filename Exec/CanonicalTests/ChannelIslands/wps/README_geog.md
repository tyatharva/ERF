# WPS_GEOG subset for the ChannelIslands NOAHMP wrfinput

The ERF Noah-MP driver reads `wrfinput_d01` (NOT `geo_em.d01.nc` — it needs
initial soil moisture/temperature, snow, TSK, CANWAT, TMN, which only
`real.exe` produces). Pipeline: geogrid -> ungrib(ERA5 GRIB) -> metgrid ->
real.exe. All scripted in `build_and_run_wps.sh`.

## Datasets to download (subset — NOT the full ~60 GB archive)

From https://www2.mmm.ucar.edu/wrf/src/wps_files/ download and extract ALL of
these into one directory, staged at `~/ERF/wps_geog/` (the extracted geog_high_res_mandatory bundle covers all of these; islope is a WPS-V3-only field and is NOT needed for V4):

| tarball | provides |
|---|---|
| `topo_gmted2010_30s.tar.bz2`                  | HGT (terrain for geogrid; ERF uses its own DEM) |
| `modis_landuse_20class_30s_with_lakes.tar.bz2`| IVGTYP / LU_INDEX, land mask (XLAND) |
| `soiltype_top_30s.tar.bz2`                    | ISLTYP (top) |
| `soiltype_bot_30s.tar.bz2`                    | soil bottom category |
| `greenfrac_fpar_modis.tar.bz2`                | VEGFRA + SHDMIN/SHDMAX (monthly min/max) |
| `lai_modis_30s.tar.bz2`                       | LAI |
| `soiltemp_1deg.tar.bz2`                       | TMN (deep soil temperature) |
| `albedo_modis.tar.bz2`                        | ALBBCK |
| `maxsnowalb_modis.tar.bz2`                    | SNOALB |

Extraction: each tarball unpacks to a directory (e.g. `topo_gmted2010_30s/`);
`geog_data_path` in `namelist.wps` must point at their common parent
(`/app/ERF/wps_geog` inside the container = `~/ERF/wps_geog` on the host).

## ERA5 GRIB reuse

`ungrib` consumes the SAME `era5_3d_*.grib` / `era5_surf_*.grib` files the
erftools step downloads — stage (or symlink) them at `~/ERF/era5_grib/`.
Only the IC time is needed for `wrfinput_d01`, so a single-time run of
real.exe suffices (namelists are set up that way).

## Order of operations

    docker run --rm -v ~/ERF:/app/ERF -w /app/ERF erf-hindcast \
        Exec/CanonicalTests/ChannelIslands/wps/build_and_run_wps.sh all

(or stage-by-stage: build, geogrid, ungrib, metgrid, real). Fill the dates in
`namelist.wps` + `namelist.input.real` first; after metgrid, verify
`num_metgrid_levels` per the script's note. Output: `wps/wrfinput_d01`, which
the ERF run directory needs alongside `NoahmpTable.TBL` and `namelist.erf`
(edit its setup-file entry to point at `wrfinput_d01`).

## Known approximation (flagged)

WPS assumes a spherical earth (R = 6370 km); erftools/ERF use WGS84. The
resulting position offset is O(100–300 m) across this 384 km domain —
negligible for 30-arcsec (~900 m) soil/vegetation categories, but do not use
this geogrid output as a source of truth for terrain height (ERF's terrain
comes from the GLO-30 DEM via dem_to_erf_terrain.py, which is
WGS84-consistent).
