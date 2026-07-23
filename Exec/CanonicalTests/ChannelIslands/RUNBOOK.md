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

## What a healthy run looks like (measured, RTX 4080, single precision)

- **Startup**: ~2 min before step 1 (init + first-call RRTMGP/CUDA JIT).
  Watch for: `Reading terrain file` + `Expecting 1154 values...` (terrain in),
  `Reading weather data/surface data 0 0 1 9` (ERA5 frames in), MM5/Morrison/
  RRTMGP/MYNNEDMF banners, then `Coarse STEP 1 starts`.
- **Timestep**: starts at ~0.2 s (`init_shrink` x CFL), grows 10%/step
  (`change_max`) toward the CFL ceiling. During the first hour of violent
  spin-up it plateaus around ~2 s, limited by vertical CFL over the San
  Gabriels (w ~ 5-7 m/s in 18.6 m terrain-compressed cells); as the interior
  equilibrates it should rise toward the horizontal advective ceiling
  (~0.7*dx/u ~ 15-20 s at jet speed). Step count for 24 h: roughly
  5,000-15,000 depending on where dt settles -- treat the first benchmark as
  the measurement.
- **Per-step wall time**: ~15-17 s/step during spin-up with
  `amrex.max_gpu_streams=1` (the stream-race workaround costs roughly 3x vs
  multi-stream; revisit when the race is fixed upstream).
- **Per-step log**: `Coarse STEP N starts/ends` + `DT =`; with
  `erf.sum_interval=1` a `TIME=/MASS/RHO THETA/...` block each step -- MASS
  ~7.2e14 and slowly drifting is healthy; `nan` anywhere = stop. Radiation
  steps add `Radiation advancing level 0 ... DONE` (every
  `rad_freq_in_steps`).
- **Plotfiles**: `plt00000` at start, then hourly (`plot_per_1 = 3600`) --
  25 plotfiles for 24 h; `chk*` every 21600 s.
- **Boundary spin-up signature** (first model hour): sponge-band winds reach
  the ERA5 values (|u| up to ~35-44 m/s) while the interior fills in from
  rest -- see Known issues #3.
- **GPU memory**: arena ~7 GB + Kokkos ~3-4 GB peaks; `nvidia-smi` ~12-14 GB
  of 16 GB. OOM in a radiation step => lower `erf.rad_ncol_chunk` (512 set)
  or the arena cap.
- **Radiation cost**: at rad_freq=30 (shipped), radiation contributed
  ~2-10%% of a clean-GPU 60-step segment (2 calls in 1045 s). rad_freq=30 at
  spin-up dt (~1.4-2 s) gives a ~1 min model-time radiation cadence,
  matching the WRF radt ~ dx(km) rule at 1 km; once dt grows post-spin-up
  the cadence coarsens (~5-8 min) — acceptable for benchmarking, revisit
  for science runs (a time-based trigger would be better).
  CAUTION: measure timings on an idle GPU — a forgotten compute-sanitizer
  container skewed an earlier timing set by >30%%; `docker ps` first.

## Known issues (all documented in the session/commits)

1. **Cross-stream GPU race** (open): with default AMReX streams the step path
   NaNs non-deterministically; `amrex.max_gpu_streams=1` (pinned in the deck)
   is bit-reproducibly clean at 1 and 30 steps. compute-sanitizer initcheck
   is clean (0 errors) after the five uninitialized-memory fixes, isolating
   a true execution race. Costs ~3x per-step; upstream-report material.
2. **Spin-up intensity — RESOLVED with cfl=0.5**: the at-rest IC + strong
   boundary jets make the first model hour violent; cfl=0.7 blew up
   deterministically near step ~45 (rad-schedule-independent, one-stream).
   cfl=0.5 is finite through the window (60-step validation, ~17.4 s/step,
   dt ~1.4 s while vertical-CFL-limited during spin-up). Double precision
   was finite as far as it fit (aborted ~step 30 on memory: full-domain
   double + RRTMGP exceeds 16 GB — double is NOT a viable fallback on this
   card, and precision is not implicated in any remaining issue). The
   flagged bulk-coeff velmag defect (Known issue #5 / UPSTREAM_ISSUES #4)
   is now fixed in this fork.
3. **HindCast pathway couples momenta only** (upstream gap, confirmed against
   upstream `development`): no ERA5 field enters the initial state (u=v=0,
   theta=300 K isentrope, qv=0 verified via max_step=0 plotfile), and the
   lateral sponge relaxes rho_u/rho_v/rho_w only -- theta/moisture are never
   anchored to ERA5. Fine for timing benchmarks; for science runs the theta/qv
   boundary-sponge extension is scoped separately. The docs (HindCast.rst)
   describe an `erf.hindcast_IC_filename` IC-from-file step that no code
   parses -- the docs-vs-code contradiction is part of the upstream report.
4. **ERA5 SST land fill values** (9999 K) are sanitized at the wiring into
   the MOST surface layer (land->288 K placeholder, sea clamped to
   [271,305] K); erftools should mask these upstream.
5. **rain_accum / bulk-coeff quirks**: the hindcast bulk-coefficient surface
   treatment read rhotheta as a velocity component -- FIXED in this fork
   (face velocities plumbed into the CC treatment; UPSTREAM_ISSUES #4).
   Upstream still has the defect; only active with hindcast_surface_bcs.

---

## Measured performance levers (2 km, 192x96x32, RTX 4080, one stream)

All numbers below are measured on 400-step runs of the shipped deck unless
noted. Baseline AFTER the check_for_low_temp printf fix: 0.0348 s/step
median (dycore+physics), 0.38 s per radiation call, 28.8 s total for 400
cold-start steps at rad_freq=10; dt plateau 2.60 s at cfl=0.5. (Before the
printf fix the same deck measured 0.159 s/step — the per-cell device
printf flood from the isentropic-IC cold top was 4.6x of runtime. Historic
numbers in earlier notes/commits carry that inflation.)

- **check_for_low_temp printf flood (FIXED)**: with the 300 K-isentrope IC,
  every cell above ~11.4 km is below the 188 K microphysics floor; the
  per-cell device printf (2 checks/step, ~30k cells) serialized the GPU.
  Now a reduction + one summary line per check (UPSTREAM_ISSUES #7). State
  evolution verified identical (MASS to 9 digits over 400 steps). After the
  theta/qv coupling task the real atmosphere clears 188 K everywhere and
  the check goes quiet on its own.
- **Stream pin re-measured**: with the flood gone, `amrex.max_gpu_streams=1`
  costs ~1.2x (0.0348 vs 0.0287 s/step at 4 streams), NOT the historic 3x
  (that ratio was flood-inflated). The race is still real: one NaN event in
  ~390 cold-start steps at 4 streams at 2 km (immediate at 1 km). Keep the
  pin.
- **rad_freq recalibrated to 30**: rad is 0.38 s/call vs 0.035 s dycore
  steps; at freq 10 it was half of wall. 30 steps = 80-135 s model-time
  cadence, matching WRF radt practice (1 min per km of dx).

- **dt is horizontal-ACOUSTIC bound**: dt = cfl*dx/(c+|u|). The per-step
  log proves it (compressible dt 2.60 vs anelastic/advective estimate ~21 s
  at step 399). Consequences:
  - `initial_dz` does NOT move dt. Measured: dz0=50 m -> dt 2.607,
    0.180 s/step (+13%); dz0=100 m -> dt 2.608, 0.275 s/step (+73%).
    Coarsening the surface layer is strictly a loss: no dt gain, higher
    per-step cost, and worse MYNN-EDMF surface-layer resolution. Keep 25 m.
  - `erf.cfl` is the lever that works. The historical 0.7 blowup is ONLY
    the IC-adjustment transient (~50-65 s model time at both 1 km and 2 km).
    Recipe: cold-start cfl=0.5 with `erf.check_int=100`, then restart
    chk00100 with `erf.cfl=0.7` (dt 3.62 s) or `erf.cfl=0.9` (dt 4.46 s);
    both validated 300 steps clean at identical per-step cost. Ship 0.7;
    0.9 has less margin.
  - Exceeding the acoustic-proportional dt fails: `fixed_dt=8` NaN'd the
    state within 20 steps (substep ratio stays 4; not an acoustic-substep
    issue).
- **Anelastic**: 0.140 s/step (~12% cheaper) and ~8x dt headroom in
  principle, but NOT usable with the HindCast pathway today -- three
  independent failures (UPSTREAM_ISSUES #6). Revisit after the theta/qv
  boundary-coupling task.
- **Stale run-dir deck**: `stage_run.sh` recopies `inputs_hindcast` for a
  reason. A stale copy silently reran the old surface_bcs=true config and
  confounded a whole measurement round. Always re-stage (or diff the run-dir
  deck against the repo deck) before benchmarking.
