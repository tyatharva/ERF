# NOAHMP / WPS assets — DEFERRED, not active

The hindcast now uses `erf.land_surface_model = MM5` (soil-column model, no
external data), so nothing in this directory is used by the current pipeline.
It is kept intact for a future NOAHMP configuration (e.g. the multi-year run,
where land-surface fidelity matters more than it does for the 75%-ocean
24-hour benchmark).

Contents (all validated to the point noted):
- `namelist.wps`          — geogrid/ungrib/metgrid config matching the exact
                            ERF LCC domain (center 33.651056/-119.290876).
- `namelist.input.real`   — minimal real.exe config (Noah-MP 4-layer soil).
- `build_and_run_wps.sh`  — builds WRF (real.exe) + WPS and runs the
                            geogrid -> ungrib -> metgrid -> real chain.
                            WRF v4.6.1 real.exe BUILT SUCCESSFULLY with this
                            script (HDF5 env must stay unset; netcdff .so
                            configure patch included). WPS geogrid was the
                            next step when the path was deferred.
- `namelist.erf`          — ERF-side Noah-MP driver namelist (Jan 9 2023,
                            crop/irrigation off => no optional geog bundle).
- `README_geog.md`        — the ~29 GB WPS_GEOG subset actually required.

To resume: re-add to the Dockerfile the pieces removed when this was parked
(gfortran, m4, csh, NetCDF-Fortran build, and the WRF/WPS build stages — see
git history of ../Dockerfile), rebuild ERF with -DERF_ENABLE_NOAHMP=ON, and
follow ../RUNBOOK.md's git history for the step-4 WPS instructions.

Note: `~/ERF/wps_geog/` (29 GB, gitignored) is still on disk and unused by
the MM5 path — kept to avoid re-downloading if NOAHMP returns.
