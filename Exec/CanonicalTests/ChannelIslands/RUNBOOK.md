# ChannelIslands hindcast — RUNBOOK

**PRODUCTION CONFIG (locked 2026-07-24): single-precision compressible,
384×192 km @ 3 km (128×64×32), NE-anchored (34.2, −117.2), hourly
plotfiles, 6-hourly checkpoints — see the PRODUCTION section at the end.**
Morrison microphysics, MYNN-EDMF PBL, RRTMGP radiation, and MM5 nominally
selected as the land surface model -- but **MM5 is INERT: land surface
temperature is exactly the constant `erf.mm5.soil_theta`, which must be set
per segment to that month's mean ERA5 skin temperature.** Read "MM5 IS
INERT" below before interpreting any land field. (NOAHMP permanently
dropped 2026-07-24, assets removed.)
Sections 0–6 below document the pipeline (written in the 2-km/34.5 era;
where coordinates differ, the deck and the PRODUCTION section govern).

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
* **PRODUCTION (34.2 pin, 3 km):** the deck uses
  `channel_islands_terrain_3km.txt` (385×193 @ 1 km, max 1972 m), generated
  with the same command but `--target-dx 1000` and the 34.2 box
  `--xlo -195131.04 --xhi 188868.96 --ylo -126372.41 --yhi 65627.59`.
  (`channel_islands_terrain_3km_345.txt` is the measured-unstable 34.5
  variant — see the pin-move section.)

## 4. Land surface: nothing to prepare

MM5 needs no static or init files -- and does no work either: it is inert
(see "MM5 IS INERT" below). The only land-surface input that matters is the
per-segment `erf.mm5.soil_theta` constant. (NOAHMP and its WPS/real.exe
input chain were dropped permanently on 2026-07-24 and the assets removed;
recover from git history before that date if ever needed.)

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

---

## Theta/qv boundary coupling (implemented; 3-km status)

The HindCast pathway now carries thermodynamics, not just momenta:

1. **ERA5 initialization** (`ERF::init_thermo_from_hindcast`): the base state
   and the state (rho, rho theta, rho qv) are rebuilt from the interpolated
   first ERA5 frame at init (ERA5 rho -> enforce_hse -> consistent p/pi/theta;
   momenta stay zero). This removes the dry-300 K-isentrope start entirely:
   no sub-188 K cold top (zero check_for_low_temp warnings), realistic
   stratification (theta 286->443 K), and it is REQUIRED for the sponge --
   relaxing band theta toward ERA5 against a 300 K interior was measured
   driving +-57 m/s w within 15 model minutes (baroclinic wall at the bands).
2. **Lateral theta/qv sponge** (`ApplyBndryForcingCC_Forecast`): relaxes
   (rho theta) and (rho qv) toward the time-interpolated frames in the same
   quadratic bands as the momentum sponge (model-rho-weighted targets; the
   momentum sponge now uses max-weight instead of additive band stacking).
3. **MYNN TKE bound** (`erf.bound_pbl_tke`, default true): prognostic RhoKE
   clipped to the scheme's internal 150 m2/s2 qke design limit.

**Validated**: 3.2-h run fully clean; hourly-plotfile physics healthy at 4 h
(band and interior theta agree to 1-2 K at every level, w bounded by mountain
waves, qv realistic). **Known limit**: the run ends at ~9.1 model hours when
boundary mass influx breaks the EOS (global mass 2.35x initial by then) --
that is UPSTREAM_ISSUES item 8 (no inflow treatment on the hindcast faces),
the remaining blocker for 24-h and multi-year runs. Deck ships the
best-measured mitigation (sponge strength 0.05 + 6th-order num_diff 0.12).

---

## Anelastic (WORKING; measured 2026-07-24)

The anelastic dycore now completes the full-physics 24-h hindcast (Morrison
+ MYNN-EDMF + RRTMGP + MM5 + real BCs). Chain of fixes, in order:

1. GMRES false convergence (dJ-deflation in TerrainPoisson::apply) --
   see UPSTREAM_ISSUES 6d/6e/6f and erf-model/ERF#3487.
2. Prescribed-inflow projection, minimal mask (32e4d49f): zero correction
   fluxes ONLY on the domain-boundary normal faces (the only faces the
   set_width=1 real-BC fill overwrites) plus terrain/lid. Do NOT mask the
   second face layer -- that traps each boundary column's integrated
   convergence and pumps rho/theta/qv (dJ-weighted, worst in the deep top
   cells). Post-solve physical divergence: max ~1.5e-9.
3. THE key bug (32e4d49f): the anelastic branch of SlowRhsPost copied raw
   momenta into avg_{x,y,z}mom, but scalar advection consumes those as
   area-weighted metric fluxes with OMEGA in the z slot. Scalars were
   advected with rho*w while rho used Omega -> exponential qv/theta pumping
   over terrain slopes (e-fold ~65 s; invisible on flat idealized cases).
4. Radiation p0 (d82162d9): with rho frozen at rho0, EOS pressure inherits
   the flat ERA5 extrapolation below its lowest level (dp ~ 0 over the first
   5 levels; exactly 0 in SP -> non-finite RRTMGP taus). Anelastic radiation
   now uses the hydrostatic reference p0(z), T = theta*Exner(p0).

**Measured (double precision build, deck + erf.anelastic=1 erf.cfl=0.3
erf.rad_ncol_chunk=1024 amrex.the_arena_init_size=3000000000):**
- 24 h complete: 8,941 steps, mean dt 9.66 s (vs 2.24 s compressible =
  4.3x fewer steps), wall 17.4 min, qv <= 0.026, T 191-304 K at h24,
  rho sum bit-identical over 24 h, zero NaNs.
- Wall-clock parity with SP compressible (17.3 min): the dt gain is eaten
  by double precision (16-GB card forces DP chunk 1024) + per-stage GMRES.
- Production options at 3 km on one RTX 4080: SP compressible or DP
  anelastic, both ~17.4 min/day -> ~4.4 GPU-days/year.

**Known limit (SP anelastic)**: single-precision anelastic runs ~3.5 h
then hits a fast (~12-step) local thermal collapse (nondeterministic
onset; transient 185-K cold pools recover earlier in the run). DP is
unaffected. Needs its own investigation before SP anelastic production.

## Pin-move measurement: 34.5 NE pin (2026-07-24)

Terrain + box for the ORIGINAL 34.5/-117.2 pin are ready:
`channel_islands_terrain_3km_345.txt` (385x193 @ 1 km, max 2906.7 m;
generated from dem.tif in the ERA5-frame LCC), box
prob_lo = -195790.23 -93110.31, prob_hi = 188209.77 98889.69
(NE pin projected into the ERA5 LCC minus 384x192 km; the 34.2 recipe
reproduces the current deck exactly, verifying the construction).

**Measured (SP compressible, production deck, cfl 0.3): UNSTABLE.**
Dies at t ~ 4.3 h (step 12,781, rad aborts on an already-NaN state).
Last healthy snapshot shows w = -32/+21 m/s breaking mountain waves
directly over the restored San Gabriels (2,900 m) at (114-118, 54-57,
k~27) -- INSIDE the northern relaxation band (10 cells from yhi). The
34.2 deck's |w| stays ~2x smaller. The 34.2 pin remains the validated
production configuration. If anyone revisits 34.5: the breaking waves
sit over peaks INSIDE the northern relaxation band, so the boundary
nudging and the wave physics are likely fighting each other -- the
first thing to try is moving or widening the band (erf.real_width, or
shifting the domain so the peaks clear the 30-km zone), not lowering
cfl. Lower cfl / stronger num_diff / vert_implicit_fac are second-line
knobs.

---

## PRODUCTION: the 1-year run (locked config, measured 2026-07-24)

**Configuration decision (measured basis):** 34.2 NE pin + SP compressible.
Both dycores cost the same wall (17.3 vs 17.4 min/day); SP compressible has
the longest validated record. 34.5 rejected (see pin-move section). The deck
in this directory IS the production deck: hourly plotfiles, 6-hourly
checkpoints.

**Cost and storage per simulated year (RTX 4080, one GPU):**
- Wall: ~4.4 days (17.3 min/day x 365).
- Plotfiles: hourly = 8,760 files x ~23 MB ~ 200 GB.
- Checkpoints: 6-hourly, 36 MB each = 1,460 x 36 MB ~ 53 GB if all kept.
  Keep a rolling tail (e.g. last 8) plus monthly keepers; a cleanup cron of
  `ls -d chk* | head -n -8 | xargs rm -rf` between segments is sufficient.

**ERA5 inputs for a year:** rerun the erftools step (RUNBOOK section 2) with
the year's dates; frames are 3-hourly. Mind CDS request limits -- chunk the
download by month. The .bin frames are interpolated to the grid at runtime,
so grid changes never require re-downloading GRIBs (only re-running the
erftools conversion if the AREA changes; the area covers both pins).

**Launch (production):**

    docker run --rm --gpus all -v ~/ERF:/app/ERF -w /app/ERF/run_hindcast \
        erf-hindcast bash -c '/app/ERF/build/Exec/erf_exec inputs_hindcast \
        amrex.max_gpu_streams=1 > run_prod.log 2>&1'

start_datetime / stop_datetime in the deck govern the segment; run in
month-long segments (stop_datetime at month end, restart into the next
month) so failures cost at most a month's tail and logs stay manageable.

**Checkpoint/restart (VERIFIED end to end, not assumed):**

    # restart from the latest checkpoint after a crash:
    .../erf_exec inputs_hindcast erf.restart=chk<NNNNN> amrex.max_gpu_streams=1

Measured verification (400-step A/B with chk at step 200, production deck):
restart loads in 0.04 s, rebuilds the boundary planes from the .bin frames,
and runs to completion with no NaNs and physical fields. Fidelity: the
restarted run differs from the uninterrupted one by up to ~5 K theta /
4 m/s u locally after 200 steps, which is ~16x the SP run-to-run
nondeterminism floor (theta 0.31 K max over 400 steps, measured A-vs-A').
So a restart is a RECOVERY, not a bitwise continuation: it splices a
slightly different realization (radiation is recomputed on the restart
step's cadence; same class of change as switching GPUs). For OSSE-truth
statistics this is benign; do not expect trajectory-level reproducibility
across a crash boundary.

**Dockerfile:** the image recipe in this directory is the production
environment (erf-env + geo stack + erftools + eccodes shim). Rebuild with
`docker build -t erf-hindcast -f Dockerfile .`; the RRTMGP SP patch is
applied via apply_rrtmgp_patch.sh (hash-gated; the runtime guard aborts
with instructions if it is missing). Build both trees inside the image
bind-mount as in section 1 (`git config --global --add safe.directory "*"`
first -- the kokkos AlwaysCheckGit probe fails on dubious ownership
otherwise).

## Radiative surface state (audited 2026-07-24; CORRECTED same day)

With zlo=surface_layer, radiation takes t_sfc per column from the surface
layer (verified in Radiation::mf_to_kokkos_buffers), and erf.rad_t_sfc=288
is an inert fallback. What that per-column value actually is:

- **over ocean:** the live, time-varying ERA5 SST. Correct.
- **over land:** the CONSTANT `erf.mm5.soil_theta`. There is no diurnal
  cycle. MM5 is inert -- see "MM5 IS INERT" below.

An earlier revision of this section claimed a "live MM5 diurnal cycle over
land". That was wrong and is retracted: land skin temperature is a single
number for the whole segment.

**Surface albedo (FIXED 2026-07-24): per-column, time-varying from ERA5
forecast albedo.** The surface frames carry fal as field 6
(erftools_fal_patch.py adds it to BOTH CDS download functions -- the
--do_forecast pipeline calls Download_ERA5_ForecastSurfaceData, which
has its own variable list; regenerate frames after patching, since
cached surface GRIBs predate the field). ERF reads it as comp 2 of the
surface state, registers alb_lev, and radiation fills the four SW
albedo slots per column (startup line: "Hindcast surface frames provide
albedo ..."). Content is the MODIS-derived climatological annual cycle
+ model snow: it captures the land-sea step (dominant), spatial spread
(in-domain land 0.08-0.16), snow events, and the (modest, ~0.02)
seasonal cycle -- NOT actual-year vegetation anomalies. Validated:
seasonal samples land mean 0.178-0.194 / ocean 0.061 on the raw grid;
24-h production run clean (38,557 steps, 16.8 min, no NaNs). Fallback
for 5-field frames: erf.rad_alb_land (default 0.17) over land via the
ERA5 mask. Remaining simplifications: uniform emissivity 0.98; one
broadband albedo fills all four SW (dir/dif x vis/nir) slots.

## Year-scale ERA5 acquisition: batched downloader (validated 2026-07-24)

erftools issues ONE CDS request PER TIMESTEP -- fine for a day, fatal for
a year (5,840 queued requests). `era5_batch_download.py` (this dir) fixes
the request pattern while keeping the whole validated processing stack:

    cd ~/ERF/era5_run
    cp ../Exec/CanonicalTests/ChannelIslands/era5_batch_download.py .
    python3 era5_batch_download.py --start 2023-01-01 --end 2024-01-01
    # then run WriteICFromERA5Data.py exactly as in section 2 -- its own
    # downloader sees every per-timestep file already present and skips
    # straight to processing.

It downloads 16-day batch requests (one per stream per chunk) and splits
them into the exact era5_{3d,surf}_YYYYMMDD_HHMM.grib files erftools
expects. MEASURED: a full-month 3-hourly request (110,112 fields) is
REJECTED by CDS ("cost limits exceeded"); 16 days (56,832 fields) is
accepted -- 75 MB, 47 min including queue. Year extrapolation: ~24 3D +
~24 surface requests, overnight hands-off. Split output verified
BIT-IDENTICAL to per-timestep downloads (444/444 pl fields and 7/7
surface fields, max abs diff 0.0). Forecast-stream surface variables
are bucketed by GRIB validity time (init-time bucketing scrambles
fluxes/zust -- the built-in uniformity check catches it).

Resumable by design: re-run after any failure -- complete chunks are
skipped via their per-timestep files, valid batch gribs are not
re-downloaded, short/corrupt batches are deleted and re-fetched. The
surface variable list includes forecast_albedo and must stay in sync
with erftools_fal_patch.py.

## Month-segment launch recipe (FRESH-INIT PARALLEL SEGMENTS -- the only supported mode)

**Each month is an INDEPENDENT fresh-init run**: its own run directory,
its own start/stop datetimes, its own ERA5 frames, its own
erf.mm5.soil_theta anchor. **Do NOT restart-chain segments**: the LSM
state is checkpointed (with the anchor ghost), so a restarted segment
silently keeps the PREVIOUS month's soil and the soil_theta knob is
inert; the boundary-plane startup cost also grows with chained frame
series. Restarts are for CRASH RECOVERY WITHIN a segment only (same
directory, erf.restart=chk<last> -- soil then correctly restores that
segment's own state).

Per segment N (month M):

    1. Batch-download + process month M (+ spin-up lead days before it,
       + the first frame of month M+1) in its own work dir; the frames
       for the segment MUST start exactly at the segment's
       start_datetime (positional indexing; the preflight enforces it).
    2. Segment deck: start_datetime = month start MINUS the spin-up
       lead (see the measured spin-up section); stop_datetime = first
       instant of month M+1. max_step = -1.
    3. Launch with erf.mm5.soil_theta=<month M value from the SKT table
       below> plus the stable-segment-config flags.
    4. In post, DISCARD the spin-up lead; analysis uses month M proper.

Segments are independent, so they parallelize freely ACROSS MACHINES/GPUs
-- but **one segment per GPU**. Two ERF processes do not fit on a 16 GB
card: RRTMGP's per-call working set makes the second process fail a 136.7
MiB Kokkos allocation mid-run (measured on an RTX 4080; both stress tests
died this way on the first attempt). Sequence segments on a single GPU.
Startup cost per segment is a month's frames (~250), not the year's.

## MM5 IS INERT. LAND TEMPERATURE IS A CONSTANT YOU SET. (audit item #3)

State this plainly, because every earlier version of this section
understated it:

**The MM5 land surface model does nothing in this configuration.** It is
not "damped", not "simplified", not "partially coupled". It is switched
off in effect, and the land surface temperature is EXACTLY the constant
`erf.mm5.soil_theta` for the entire segment.

Mechanism (measured, not inferred). SurfaceLayer only wires itself to an
LSM when that LSM exports at least four flux fields:

    AMREX_ALWAYS_ASSERT(n_valid_lsm_flux == 0 || n_valid_lsm_flux >= 4);

MM5 exports one. The count is therefore 0-or-1, the coupling pointers stay
null, and NOTHING passes in either direction -- no flux down into the soil,
no soil temperature back up. Verified directly with the MM5SOIL diagnostic
(`erf.mm5.diag_interval`): the soil column is bit-identical at all four
diagnostic depths (0.1 / 0.5 / 1.5 / 3 m) after three simulated days.

What that means in practice:

- Land surface temperature = `erf.mm5.soil_theta`, constant in time and
  uniform across all land points.
- That constant is what MOST uses to compute land surface-layer fluxes.
- That constant is what RRTMGP uses as the land radiating temperature.
- There is NO land diurnal cycle at all: not damped, absent.
- There is NO soil spin-up, because nothing evolves. Do not budget time
  for one; the 266-day relaxation timescale estimated earlier is moot.
- Absorbed shortwave (now with correct ERA5 albedo) heats the ATMOSPHERE
  via MOST only; the ground never receives it.

Ocean -- 75% of this domain -- is unaffected: SST comes from ERA5 and is
live and time-varying. Land-driven flows (sea breeze strength, nocturnal
drainage, inland convective triggering) carry a no-diurnal-cycle bias of
unquantified size. **Do not use land skin or soil temperature from the
output as a physical field.** It is a boundary condition you chose.

The only way to get real land radiative driving in this fork is Noah-MP,
dropped permanently 2026-07-24 with its geog-data dependency. ERF's other
LSM (SLM) is the same bare column and exports the same single flux, so
there is no intermediate option.

## Per-segment soil anchor -- USE MONTHLY-MEAN SKIN TEMPERATURE (skt)

Because land temperature IS `erf.mm5.soil_theta`, this is the single most
consequential land setting in the deck, and it must be right per segment.

Set it to that month's domain land-mean ERA5 **skin temperature (skt)**.
Not the hardcoded 300 K default. Not deep-soil stl4.

Why skt and not stl4: an inert MM5 means the knob is not initializing a
soil column that will evolve toward radiative equilibrium -- it is directly
prescribing the radiating/flux-exchanging surface. The physical quantity
that plays that role in ERA5 is skt. stl4 (100-289 cm) is the right depth
to anchor a *working* soil column and the wrong quantity to impose as a
surface temperature; it lags the surface by months (in this domain stl4
peaks in September while skt peaks in July) and would systematically
mis-time the seasonal cycle.

MEASURED 2023 domain land-mean skt (K) -- use these:

    Jan 283.2  Feb 283.2  Mar 284.4  Apr 289.2  May 291.5  Jun 293.9
    Jul 300.1  Aug 298.9  Sep 296.1  Oct 293.9  Nov 288.9  Dec 286.5

The hardcoded 300 K default happens to be right for July and is ~17 K too
warm in January. Fetch with CDS reanalysis-era5-single-levels-monthly-means,
variables skin_temperature + land_sea_mask, area 34.2,-121.2,32.8,-117.2,
masking to land points. For other years fetch that year's own means; the
table above is an acceptable climatology fallback.

    (Superseded stl4 table, kept only so nobody re-derives it and thinks it
    is the recommendation: Jan 289.6 Feb 287.5 Mar 286.2 Apr 286.2 May 287.9
    Jun 289.6 Jul 292.0 Aug 294.9 Sep 296.1 Oct 295.7 Nov 294.5 Dec 292.3.
    DO NOT USE.)

The LSM state IS checkpointed (with ghosts), so on a RESTART-chained
segment the restored soil wins and this knob is silently inert. The
production workflow is fresh-init parallel segments precisely so the knob
takes effect -- see the launch recipe above.

## FINAL SEGMENT CONFIG (locked 2026-07-24)

Every segment carries all of these. They are each traceable to a measured
failure; do not drop one because a particular month runs without it.

    erf.cfl                       = 0.2       # set by Hurricane Hilary
    erf.max_dt                    = 2.5       # set by a calm January start
    erf.moistscal_horiz_adv_type  = Upwind_3rd
    erf.moistscal_vert_adv_type   = Upwind_3rd
    erf.advect_tke                = false
    erf.hindcast_blend_band_density = true
    erf.mm5.soil_theta            = <month skt from the table below>
    erf.check_for_nans            = 1
    erf.check_for_nans_int        = 100
    max_step                      = -1

The two dt controls do different jobs and BOTH are needed: `max_dt` binds on
quiet days (where the CFL estimate would let dt run to ~15 s), `cfl` binds on
active days. Because quiet days are capped by `max_dt` anyway, lowering cfl
from 0.3 to 0.2 costs wall time only on convectively active days.

## Stress tests: the year's extreme regimes (2026-07-24)

Jan 9, the only day validated during bring-up, is mid-pack for this domain.
A survey of all 2023 ERA5 (4x daily, whole domain) picked out the real
extremes, and the config above was validated against the worst two:

| Case | Date | Why | Result |
|---|---|---|---|
| Hurricane Hilary | 2023-08-20 | TCWV 66 mm, 2.5x the Jan-9 storm | see below |
| Max CAPE | 2023-09-09 | 1404 J/kg, TCWV 44 mm | see below |
| Windstorms | 2023-02-15, 02-22 | surface winds 18 m/s | not run |

**What Hilary broke.** At the bring-up cfl of 0.3 the Hilary case dies ~9
model minutes in: a grid-scale vertical dipole (w = +48.9 / -29.5 m/s) in the
lowest three cells of the ylo relaxation band at (54,1,1..3), w-CFL 2.87 in
the 18.5-m surface-compressed cells, driving RhoTheta to NaN via a -220 K
temperature. The dt estimate never anticipates it. cfl 0.25 still dies at the
same point; cfl 0.20 clears the whole transient with zero w-damping events
and zero low-temperature warnings.

**It is marginal, not cleanly reproducible.** At cfl 0.3 with the density
blend off, two byte-identical invocations gave 250 clean steps and a step-240
blowup respectively (UPSTREAM_ISSUES #13). When probing stability at the
margin in this SP GPU build, the acceptance bar is **zero w-damping and zero
low-temperature warnings**, not merely "no NaN in one run".

**One segment per GPU.** Both stress tests initially died from a failed
136.7 MiB Kokkos allocation because two ERF processes shared a 16 GB RTX
4080 (one process resident set: ~7.4 GB, plus RRTMGP's per-call working
set). Parallelize segments across GPUs, never within one.

## Calm-start stability findings (segment probes, 2026-01-01 start)

A Jan-1 (calm) start exposed three stacked failure modes that every
Jan-9 (storm) validation missed; all are resolved and the STABLE
SEGMENT CONFIG is:

    erf.max_dt=2.5  (dt ceiling: calm-day advective bound otherwise
        lets dt grow to ~15 s and cross the explicit vertical-diffusion
        stability limit)
    erf.moistscal_horiz_adv_type=Upwind_3rd
    erf.moistscal_vert_adv_type=Upwind_3rd
        (WENOZ5 on Morrison scalars NaN'd a number-concentration
        component at step 16 in boundary-band rows; upwind is clean)
    erf.advect_tke=false
        (RhoKE went Inf at the xhi relax band ~380 s after the first
        3-h frame transition with TKE advected under WENOZ5 dryscal;
        not advecting TKE is standard for MYNN-class schemes --
        production/dissipation dominate at 3 km)

DO NOT set dryscal_*_adv_type=Upwind_3rd: that path is broken in this
fork (density NaN at step 1 in corner cells).

Verified: 24 h Jan-1 with erf.check_for_nans=1 armed: 37,798 steps,
8 frame transitions, zero non-finite trips. erf.check_for_nans also
detects Inf and prints locations (contains_nan alone misses Inf -- an
SP overflow becomes NaN only one step later via Inf-Inf).
