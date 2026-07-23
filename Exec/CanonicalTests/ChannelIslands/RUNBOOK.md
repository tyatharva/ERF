# ChannelIslands 24-h hindcast — RUNBOOK

Single-precision ERF, 384×192×48 @ 1 km, NE-anchored (34.5, −117.2).
Morrison microphysics, MYNN-EDMF PBL, RRTMGP radiation, MM5 land surface
(soil-column model -- compiled in, no external data; NOAHMP assets parked in
wps_noahmp_deferred/).
2023-01-09 00:00 UTC → 2023-01-10 00:00 UTC, hourly plotfiles.

All commands run in the `erf-hindcast` container with the repo bind-mounted:

    alias erfrun='docker run --rm --gpus all -v ~/ERF:/app/ERF \
        -v ~/.cdsapirc:/root/.cdsapirc:ro -w /app/ERF erf-hindcast'

The CDS key is mounted read-only at runtime; it is never in an image layer.

---

## 0. Build the image (once)

    cd ~/ERF/Exec/CanonicalTests/ChannelIslands
    docker build -t erf-hindcast -f Dockerfile .

* **What it does:** extends `erf-env` with the Python geo stack and erftools
  (source-tree via PYTHONPATH, with the eccodes flux-name compat shim applied
  and verified). No Fortran toolchain: the MM5 land model needs none.
* **How long:** ~5–10 min first time (pip layer dominates); tweaks ~1–2 min.
* **Success:** `Successfully tagged erf-hindcast:latest`; the build hard-fails
  if erftools doesn't import.
* **Failure modes:**
  - erftools import error → upstream erftools changed; see
    `erftools_eccodes_compat.py` (it hard-fails with instructions).
  - pip resolution errors → the dep list in the Dockerfile mirrors erftools'
    declared + undeclared imports; a new upstream import shows up as a
    ModuleNotFoundError at step 2 instead.

## 1. Build single-precision ERF (once per source change)

    erfrun bash -c 'git config --global --add safe.directory "*"; \
        ERF_HOME=/app/ERF Build/cmake_single_precision_cuda.sh'

* **What it does:** applies the RRTMGP single-precision patch to the submodule
  (byte-exact sha256 check — the build hard-fails on an unknown submodule
  state), configures with `ERF_PRECISION=SINGLE`, RRTMGP/NetCDF/FFT on
  (no NOAHMP flag — MM5 is always compiled in), and builds into `~/ERF/build`.
* **How long:** ~35 min clean; seconds–minutes incremental.
* **Success:** `[rrtmgp-sp patch] already applied (hash verified)` (or
  `applied and hash-verified`), then `[100%] Built target erf_exec`.
* **Failure modes:**
  - `FATAL [rrtmgp-sp patch] ... matches NEITHER ...` → the RRTMGP submodule
    moved; do NOT force — re-derive the patch.
  - Verify precision: `grep '#define AMREX_USE_FLOAT'
    build/Submodules/AMReX/AMReX_Config_3D.H` must match (double builds show
    `/* #undef ... */`).

## 2. ERA5 download + IC/BC processing (erftools)

    mkdir -p ~/ERF/era5_run && cd ~/ERF/era5_run
    cat > era5_input.txt <<'EOT'
    year: 2023
    month: 01
    day: 09
    time: 00:00
    area: 36.0,-123.25,31.25,-115.25
    EOT
    erfrun bash -c 'cd era5_run && \
        cp /opt/erftools/notebooks/era5/WriteICFromERA5Data.py . && \
        cp -r /opt/erftools/notebooks/gfs/TypicalAtmosphereData . && \
        mpirun -n 8 python3 WriteICFromERA5Data.py era5_input.txt \
            --do_forecast=true --forecast_time_hours=24 --interval_hours=3'

* **What it does:** downloads ERA5 surface + pressure-level GRIBs for 9
  three-hourly frames (00Z Jan 9 → 00Z Jan 10) over the analysis box, projects
  them into the LCC frame, and writes ERF binary IC/BC frames.
* **How long:** download 5–40 min (CDS queue-dependent; GRIBs are cached in
  `era5_run/` so reruns skip finished files); processing ~3–5 min on 8 ranks.
* **Output:** `era5_run/era5_{surf,3d}_*.grib` (9+9),
  `era5_run/Output/ERA5Data_3D/*.bin` (9), `Output/ERA5Data_Surface/*.bin` (9),
  `Output/domain_extents.txt`, VTK previews under `Output/VTK/`.
* **Success:** exit 0; 9 `.bin` files in each Output subdir;
  `domain_extents.txt` present. Sanity: our fixed domain
  (x −195790..188210, y −93110..98890) must sit inside the printed
  prob_lo/prob_hi box — it does, with ~15–27 km to spare; ERF re-checks at
  startup and aborts with "The xlo value of the domain has to be greater
  than ..." if not.
* **Failure modes:**
  - `403 ... required licences not accepted` → accept the licence on BOTH
    ERA5 dataset pages (single-levels and pressure-levels) in the CDS web UI.
  - `need at least one array to stack` → eccodes renamed the flux variables;
    the image applies `erftools_eccodes_compat.py` (rerun in a fresh image if
    you see this).
  - `TypicalAtmosphereData/... not found` → the `cp -r` of that directory was
    skipped (the era5 notebook doesn't ship it; it's taken from the gfs one).
  - `ModuleNotFoundError` → erftools has undeclared imports; the image
    preinstalls the known set (cartopy, herbie-data, pygrib, tqdm, dateutil,
    yt, ussa1976, ...). A new one means upstream added a dependency.

## 3. Terrain (once; already generated and validated)

    erfrun bash -c 'cd Exec/CanonicalTests/ChannelIslands && \
        python3 dem_to_erf_terrain.py --dem dem.tif \
            --out channel_islands_terrain.txt \
            --proj-str "+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 +lon_0=-119.250000 +datum=WGS84 +units=m +no_defs" \
            --xlo -195790.23 --xhi 188209.77 --ylo -93110.31 --yhi 98889.69'

* **Output:** `channel_islands_terrain.txt` (1154×578 @ 333 m, LCC meters).
* **Success:** "elev min/max = 0.0/3040.7 m"; ~8.75% zero-filled (verified
  open-ocean strip the DEM doesn't cover). ERF prints
  `Reading terrain file: ...` + `Expecting 1154 values of x, 578 ... 667012 ...`
  at startup — if the file were missing ERF aborts (no silent flat fallback).

## 4. Land surface: nothing to prepare

MM5 needs no static or init files. (The former step 4 — the WPS/real.exe
chain producing `wrfinput_d01` for NOAHMP — is parked in
`wps_noahmp_deferred/`; see its README to resurrect it.)

## 5. Assemble the run directory

    erfrun Exec/CanonicalTests/ChannelIslands/stage_run.sh

* Collects into `~/ERF/run_hindcast/`: `inputs_hindcast`, terrain, ERA5 3D +
  surface frames (as `ERA5Data_3D/`, `ERA5Data_Surface/`), and the 4 RRTMGP
  tables (from the repo — no download). Hard-fails naming the missing
  ingredient otherwise.
* **Success:** `staged: 9 3D frames, 9 surface frames` + `Run directory ready`.

## 6. The 24-hour run

    erfrun bash -c 'cd run_hindcast && \
        /app/ERF/build/Exec/erf_exec inputs_hindcast |& tee run.log'

See "What a healthy run looks like" below.

---

## What a healthy run looks like

*(numbers below are from the validation runs on the RTX 4080; see the
handoff message for the measured values)*
