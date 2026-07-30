# Upstream-reportable findings (erf-model/ERF and erftools)

Compiled during the ChannelIslands single-precision hindcast bring-up
(this fork, branch `ERF`; each item verified against upstream
`erf-model/ERF@development` where noted).

## 1. HindCast pathway: state is never initialized from data; theta/moisture never forced (docs contradict code)

**Documentation says** (`Docs/sphinx_doc/theory/HindCast.rst`, step 4):

> "Copy the initial condition file to the ERF run directory. This can be the
> first file in the ``hindcast_boundary_data_dir``."

and its inputs example lists:

> ```
> // Initial condition filename -
> // obtained from running the python script
> erf.hindcast_IC_filename = "ERF_IC_2025_08_18_00_00_000.bin"
> ```

**Code does**: no source file parses `hindcast_IC_filename` (grep of
`Source/` — zero hits; confirmed identical in upstream `development`'s
`ERF_WeatherDataInterpolation.cpp`). The interpolation machinery reads rho,
u, v, w, theta, qv, qc from the `.bin` frames but writes **only `Rho_comp`**
(plus lat/lon) into the forecast arrays; theta and moisture are interpolated
into temporaries and discarded. Nothing fills `vars_new` from data at
initialization, and the lateral sponge (`ERF_ApplyBndryForcing_Forecast.cpp`)
relaxes **momenta only**.

**Evidence** (this fork, `max_step=0` plotfile of the full hindcast config):
`x_velocity = y_velocity = 0` identically, `qv = 0` identically,
`theta = 300.000 +/- 0.01 K` (the dry-isentropic HSE base state), i.e. the
"hindcast" initial state contains no ERA5 information whatsoever, and during
the run only momentum is anchored to ERA5. Verified live boundary forcing:
after 30 s of model time the sponge bands carry |u| up to 35.7 m/s (matching
the ERA5 data) while the interior spins up from rest — so the momenta path
works; theta/qv have no large-scale anchor at all.

**Suggested fix**: either implement the documented IC read + extend the
sponge to `RhoTheta`/`RhoQ1` (the fill kernel already interpolates them), or
correct the documentation. WRF-style practice (wrfbdy carries theta and
moisture) argues for the former.

## 2. GPU cross-stream race in the hindcast step path

With default AMReX stream counts, the full hindcast configuration produces
NaN in memory-dependent, run-to-run-varying regions within one step. With
`amrex.max_gpu_streams=1` the identical binary/config is clean and
bit-reproducible (1-step and 30-step verified). `compute-sanitizer
--tool initcheck` reports zero errors on the failing configuration (after
the fixes in item 3), isolating a true execution-order race rather than
uninitialized memory. Racing kernel pair not yet identified.

Additional bracketing (2-km config, all item-3/4 fixes in): the race still
exists but is RARE at 6 boxes — one NaN event in ~390 cold-start steps at
4 streams (corruption in a plain dycore step, 8 steps after the last
radiation call, with the surface path disabled — so it is in the ordinary
hindcast step path, not radiation-interop or surface code), while two
200-step multi-stream runs were clean and mass-identical to one-stream.
At 1 km (24 boxes) it fires within one step. Frequency scaling with box
count fits a cross-stream timing-window race. Reproduction: this fork @
`ERF` branch, `Exec/CanonicalTests/ChannelIslands/`, remove
`amrex.max_gpu_streams=1`; use the 1-km config (git history, commit
0bd051d5) for an immediate repro.

## 3. Uninitialized-memory defects found via compute-sanitizer initcheck

All fixed in this fork (commits noted in git log); upstream should verify
which apply to `development`:

- `ERF_SurfaceDataInterpolation.cpp`: `ProbLo()`/`CellSize()` host pointers
  captured in a device lambda (CUDA error 700 on first surface-data read).
- `FillSurfaceStateMultiFabs`: writes only the k==0 plane of a 2D slab that
  inherits the 3D state's ghost vector (z-ghost planes left unwritten and
  then copied/blended).
- `SurfaceDataInterpolation`: the `surface_state_interp` time-blend LinComb
  is commented out upstream while the slow-RHS consumes the MultiFab
  whenever `hindcast_surface_bcs` is on (garbage land-sea mask).
- `WeatherDataInterpolation`: forecast velocity fill kernels cover valid
  face boxes only; the interp LinComb blended their uninitialized ghosts.
- `ERF_ProbCommon.H read_custom_terrain` (plain format): query points at or
  epsilon past the last terrain grid line read past the coordinate arrays
  (OOB). Terrain files that end exactly at the domain edge trigger it.
- `ERF_TerrainMetrics.cpp` STF branch: `h_mf`/`h_mf_old` FillBoundary'd
  before all planes are written.

## 4. Hindcast bulk-coefficient surface treatment misreads state components

`ERF_ApplySurfaceTreatment_BulkCoeff.cpp` (CC variant) computes
`uvel = cons_state(i,j,k,1)/cons_state(i,j,k,0)` — component 1 is
`RhoTheta`, not momentum, so "velmag" evaluates to ~theta ≈ 300 (m/s scale)
over sea. Latent/sensible bulk fluxes are then ~an order of magnitude too
strong. Only active with `erf.hindcast_surface_bcs = true`; all-ocean
hurricane cases may have tolerated it. **Fixed in this fork**: the CC variant
now takes the face-velocity arrays (`u_arr`/`v_arr`, averaged to cell
centers) plumbed through `make_sources`; see
`ERF_ApplySurfaceTreatment_BulkCoeff.cpp` / `ERF_MakeSources.cpp` on branch
`ERF`. Upstream still has the defect.

## 5. erftools

- Modern eccodes renamed accumulated flux fields
  ("Surface sensible heat flux" -> "Time-integrated surface sensible heat
  net flux"); `ReadERA5_SurfaceData`'s substring matches collect nothing and
  crash with "need at least one array to stack"
  (`erftools_eccodes_compat.py` in this directory patches it).
- ERA5 SST carries ~9999 K fill values over land; erftools writes them into
  the surface binaries unmasked, and bilinear interpolation smears them into
  coastal sea cells (values of several thousand K observed). Should be
  masked/filled at write time.
- `pyproject.toml` packaging omits subpackages (a pip install succeeds but
  `erftools.preprocessing` is missing); several imports (cartopy,
  herbie-data, pygrib, tqdm, dateutil, yt) are undeclared.
- The era5 notebook requires `TypicalAtmosphereData/` but only the gfs
  notebooks ship it.

## 6. Anelastic mode is incompatible with the HindCast pathway (three independent failures)

Measured on the 2-km ChannelIslands config (this fork, branch `ERF`):

- **Adaptive dt is unbounded from the at-rest IC.** The HindCast IC starts
  the state at rest (item 1), so the anelastic dt estimator (advective CFL
  only — no acoustic constraint) returns "undefined"/unbounded; the run
  reached 12 model-hours in 4 steps, NaN'd, and aborted in the hindcast
  frame-interpolation weights (`alpha1 = -0.084`).
- **Pressure blow-up at the lateral momentum sponge.** With `fixed_dt=2.6`
  the run survives ~50 steps, then the state fed to radiation shows
  `p_lay` up to 1718 hPa and temperatures at the gas-optics table ceiling
  (355 K), localized at the north-edge sponge — the anelastic pressure
  projection and the hindcast momentum sponge are inconsistent at the
  boundary. (Compressible with the identical deck runs clean.)
- **No compressible→anelastic restart.** Anelastic restart requires
  `Level_0/PP_Inc_H`, which compressible checkpoints do not write, so the
  obvious workaround (spin up compressible, switch to anelastic) aborts at
  I/O.

Worth fixing upstream: measured anelastic per-step cost on this config is
~12% BELOW compressible (no acoustic substepping outweighs the FFT-
preconditioned GMRES Poisson solve), and the advective dt estimate is ~8x
the compressible (acoustic-bound) dt — a large win if the boundary
treatment supported it.

## 7. check_for_low_temp floods stdout with the HindCast IC

The dry-isentropic 300 K IC (item 1) puts every cell above ~11.4 km below
the 188 K microphysics validity floor, so `ERF::check_for_low_temp`
(called twice per step with any moisture model) emits ~30,000 device
printfs per step — ~12M warning lines in a 400-step run. The `Abort()` in
the device lambda does not fire in release GPU builds, so the run
continues; the flood is log bloat plus SEVERE printf overhead: fixing it
(reduction + one summary line per check, done in this fork's `ERF.cpp`)
took the full 2-km hindcast step from 0.159 to 0.0348 s/step — the
per-cell device printf was 4.6x of total runtime, and it also scales
pathologically with GPU stream count (4 streams: 0.68 s/step WITH the
flood vs 0.0287 WITHOUT). State evolution verified identical (MASS to 9
digits over 400 steps). Upstream should adopt the summary form.

## 8. HindCast lateral boundary has no inflow treatment: domain mass grows without bound

Measured at 3 km (all numbers from `run_hindcast` logs on branch `ERF`,
with the theta/qv coupling of item 1 implemented so runs survive past
spin-up): with `Outflow` (foextrap) lateral faces and the momentum sponge
forcing ERA5 winds, mass flows in through inflow faces with no
characteristic/specified treatment and never leaves. **Global mass reaches
2.35x its initial value after 9 model hours**; density in the sponge bands
grows ~0.1-0.3 kg/m3 per hour (worst at band-overlap corners) until the
surface pressure breaks the EOS/radiation (~1.9-2.4 kg/m3 at blowup).

Sensitivity (24-h attempts, 3-km config, blowup time):
- sponge strength 0.3, no diffusion: 4.6 h
- + 6th-order numerical diffusion (0.12): 5.2 h
- + max-weight (non-additive) band overlap: 5.5 h
- strength 0.05 + diffusion: 9.1 h (rate halves, does not vanish)

Approaches that made things WORSE (documented so nobody retries them):
- RHS mass-equation relaxation toward forecast rho: destabilizes the
  acoustic substepping within ~60 steps.
- Forecast-rho weighting of the theta/qv sponge targets (rho_f*theta_f):
  NaN within ~50 steps.
- Direct post-step density blending in the bands: NaN at ~800 steps
  (alpha~0.9/step) and ~420 steps (alpha~0.04/step).

**Fix needed**: a real specified/characteristic inflow boundary (WRF-style
specified zone at the face, or characteristic-based inflow BC) for the
HindCast pathway. The momentum sponge alone over-determines the velocity
while leaving the mass/pressure at the face unconstrained.

Related port defect found on the way: the MYNN-EDMF column solver clips its
internal qke at 150 m2/s2, but the clip never reaches the dycore-advected
prognostic RhoKE (in WRF the clipped qke IS the prognostic state). Measured
consequence: TKE ran past 200 m2/s2 at 8-15 km (2.5-km cells, strong jet
shear) on this config. Fixed in this fork by bounding the prognostic state
each step (`bound_mynn_tke` in ERF_Advance.cpp, `erf.bound_pbl_tke`).

## 9. RESOLVED IN FORK: real-BC machinery repointed at hindcast frames (fixes item 8)

Rather than writing a new inflow BC, this fork drives the existing
specified/relaxation-zone machinery (`fill_from_realbdy`,
`realbdy_compute_interior_ghost_rhs`) from the interpolated hindcast frames:
`ERF::fill_bdy_data_from_hindcast` sizes and fills
`bdy_data_{xlo,xhi,ylo,yhi}[time][U,V,T,QV]` exactly as `init_from_metgrid`
does (plain u/v/theta/qv strips; ParallelCopy gather + broadcast), gated by
`erf.use_real_bcs` with `init_type = HindCast` (sanity-check assert relaxed
accordingly). **Measured result: global mass drift +0.25% over a full 24-h
run** (item 8 measured +135% in 9 h before), first complete 24-h hindcast of
this campaign. The machinery is init-type-agnostic by construction —
upstream may want this adapter, since it gives the ERA5/GFS route the same
boundary treatment as the WRF-file route with ~150 lines of marshaling.

Note for upstream: the interior relaxation (`realbdy_compute_interior_ghost_rhs`)
nudges U/V/T but NOT QV, although the boundary planes carry QV and the
specified zone sets it. WRF relaxes moisture in the zone. Not yet extended
here (a sharp qv step stands at the specified/relax interface; it was NOT
the cause of our storm-time NaN, but is worth closing for consistency).

## 10. estTimeStep omits the vertical ADVECTIVE constraint under implicit substepping

`ERF::estTimeStep`'s compressible branch drops the z-direction entirely when
`substepping_type = Implicit` ("the z-direction does not contribute") --
correct for vertical ACOUSTICS (integrated implicitly) but wrong for
vertical ADVECTION, which remains explicit in the slow RHS. Measured
consequence (3-km ChannelIslands, 18.5-m terrain-compressed surface cells):
a resolved convective cell reached w ~ 10 m/s while the estimator held
dt = 3.9 s (advective limit ~0.9 s) -> NaN. This also explains previously
observed cfl-0.9 and fixed_dt=8 failures on this config. Fixed in this fork
by adding a |w|/dz_local advective term (detJ-based local cell thickness so
deep cells aloft are not over-constrained) to the substepping branch.
With the fix the adaptive dt visibly tracks convective pulses (drops to
~1.4 s during storms). Explicit moist convection additionally needs
cfl <= 0.3 on this grid (w can grow faster than one step's margin).

## 6b. Anelastic update (with the theta/qv + real-BC atmosphere)

Re-tested after items 8-10 were resolved: the at-rest adaptive-dt failure
and the sponge pressure blow-up of item 6 are CLEARED (anelastic dt is
well-defined at 9.3-9.7 s vs compressible 2.1-2.3 s -- a measured ~4.5x dt
headroom -- and T/p stay physical). A third blocker remains: with
`use_real_bcs` specified-inflow boundaries, the anelastic pressure
projection mis-evolves the state from the start (theta gains +130 K and u
reaches 97 m/s within 30 steps, silent NaN by ~step 100 with radiation
off). The Poisson solve's lateral BCs do not account for specified
(non-divergence-free-consistent) inflow. Needs upstream-grade work on the
projection boundary conditions before anelastic + real BCs is usable.

## 6c. Anelastic scoping result: structural, precision-independent top-level runaway

Systematic elimination (3-km config, real BCs, theta/qv init; all runs
reproduce from cold start):

- theta runs away in the TOP TWO LEVELS only (healthy 387 K at k=29 ->
  452 -> 567 K at the lid within 30 steps, ~5 K/step; NaN by ~step 100).
  Lateral bands are NOT the locus (9% of hot cells).
- Radiation off: unchanged. zhi w-sponge off: unchanged. cfl 0.3: unchanged.
- DOUBLE PRECISION reproduces theta=567.33 at the same cell to 5 digits:
  precision-independent.
- The Poisson correction is structurally inconsistent with the recomputed
  divergence: GMRES converges its own residual to tolerance (7 iters,
  2.8e-4) yet the post-correction divergence is only 4x (SP) / 8-25x
  (double) smaller than pre-solve -- the operator and the
  getFluxes/compute_divergence path disagree. Also the post-solve
  volume-weighted divergence sum jumps to -1.37e6 IDENTICALLY in both
  precisions (structural, possibly a diagnostic artifact -- unverified).
- Suspect space that remains: Omega/w handling in the deep stretched top
  cells (dz ~2.5 km), buoyancy type 3 against a 3-D data-derived base
  state, theta advection with the residual divergence. The anelastic x
  terrain x real-data combination appears never to have been exercised
  upstream (anelastic is used for idealized LES).

Scope estimate for a fix: week-class with tail risk (structural dycore
excavation, not BC glue). The measured prize if fixed: anelastic dt
9.3-9.7 s vs compressible 2.1-2.3 s at identical per-step cost -- ~4x.

## 6d. Anelastic root cause narrowed to the GMRES solve itself (operator/flux exonerated)

Synthetic-phi consistency test (erf.poisson_consistency_test=1, this fork):
TerrainPoisson::apply and the applied-correction path
(getFluxes -> momenta increment -> compute_divergence) agree to round-off
(norms match to 10 digits; div(flux) = -A phi exactly, the intended sign).
The item-6c "operator vs correction inconsistency" is EXONERATED.

The true defect (erf.poisson_consistency_test=2): ||rhs - A phi_solved||
equals the leftover post-projection divergence to every digit, i.e. the
GMRES recurrence residual is fictitious: the solver reports ~1e-8 relative
convergence while the TRUE operator residual stalls at 2-44 percent per
solve (double precision; SP similar). Iterative refinement on the true
residual DIVERGES deterministically (0.236 -> 0.45 L2 after 10 passes),
so each solve's output can be anti-correlated with its own residual --
the Krylov mechanics are broken for this operator/preconditioner pair,
not merely slow. Eliminated: rhs in-place mutation by the FFT
preconditioner (solving on a copy is bit-identical); apply_bcs
inhomogeneity (it is linear/homogeneous); precision (double reproduces).
Remaining suspects: Krylov-basis corruption via the apply() const_cast
ghost mutation, catastrophic orthogonality loss in the Gram-Schmidt for
this A*Minv, or an inconsistent singular system (all-Neumann null space
vs the dJ-weighted mean subtraction). Consequence in production: the
un-projected divergence remainder pumps (rho theta) at ~5 K/step near
the lid until NaN (~step 100).

## 6e. Probes (c) and (b): the defect is inside AMReX GMRES's recurrence

Probe (c) -- null-space compatibility (erf.poisson_consistency_test=3):
- The discrete left-null vector of the terrain operator is dJ, NOT the
  constant vector: <A x, dJ> normalized ~ 2e-6..4e-9 for arbitrary x while
  <A x, 1> ~ 7e-5.
- The production mean subtraction removes exactly the right component:
  <rhs, dJ> ~ 1e-17 after subtraction. The system is CONSISTENT.
- The stalled residual is orthogonal to BOTH candidate null vectors
  (cosines ~1e-5): the stall is not a compatibility floor. CLEAN.

Probe (b) -- unpreconditioned A/B (erf.poisson_consistency_test=4):
- WITHOUT the FFT preconditioner, GMRES is HONEST: the recurrence stagnates
  at ~16 percent relative after 2000 iterations and reports so (true
  residual matches, 0.045 vs reported 0.038).
- WITH the preconditioner it reports 1e-8 while the true residual is
  2-44 percent.

Mode 5 -- preconditioner properties:
- Deterministic to the bit (|M v - M v| = 0 across calls) and linear to
  double round-off (relative 3e-16).
- Amplification is large and legitimate: |M^-1 v| / |v| ~ 1e8 (Poisson
  inverse of the lowest modes).
- Rescaling M^-1 by a frozen constant does NOT restore honesty (tried,
  reverted).

Conclusion: every GMRES precondition (operator linearity/determinism,
preconditioner linearity/determinism, flux/operator consistency, system
compatibility, precision) is verified; honesty appears/disappears with the
preconditioner alone. The false-convergence defect is inside
amrex::GMRES's recurrence/orthogonalization bookkeeping when paired with
this legitimate M^-1 (possibly its Gram-Schmidt under the ~1e8 spectral
spread of A*M^-1). Recommend filing against AMReX with the two one-flag
reproducers in this fork (poisson_consistency_test=2 vs =4).

## 6f. FIXED IN FORK: singular-pair GMRES false convergence (root cause of 6d/6e)

Root cause found and fixed ERF-side (no AMReX patch needed): the terrain
Poisson operator is singular (constants in null(A), dJ spans null(A^T))
and the FFT preconditioner is singular too (pinned mean). GMRES on the
singular pair lets the null component re-enter through round-off, the
Hessenberg goes near-singular, and amrex::GMRES's unprotected backsolve
(`m_grs[it] /= m_hh(it,it)`, exact-zero check only) produces corrupted
updates while the Givens recurrence reports convergence -- explaining the
false 1e-8 reports, the 2-44 percent true-residual stalls, and the
minimizer-property violation (updates worse than zero).

Fix: deflate the singular mode in TerrainPoisson::apply -- project the
operator output orthogonal to dJ (solves P A phi = rhs; equivalent for
the compatible rhs whose dJ-component is ~1e-17 after the production
mean subtraction). MEASURED: true residual after solve drops from 2-44
percent to 5e-9 (meets tolerance, zero refinement passes); post-projection
divergence drops from 7.5e-4 to 1.9e-9 max -- seven orders. Upstream
recommendations: (1) this deflation in ERF's TerrainPoisson; (2) a
relative near-zero pivot guard in amrex::GMRES build_solution; (3) a true-
residual verification before declaring convergence (defense in depth).

REMAINING anelastic blocker (next layer, precisely characterized): the
real-BC specified zone overwrites boundary-strip velocities AFTER each
projection, re-injecting divergence there every stage. Compressible
absorbs this acoustically; anelastic advects scalars with the divergent
velocity and they accumulate without bound (measured with the FIXED
solver: qv reaches 0.8 kg/kg at a boundary-adjacent surface cell by step
80; the top-level theta runaway persists at ~1/3 its former rate).
Proper fix: projection with prescribed inflow -- treat specified-zone
faces as fixed-flux (inhomogeneous Neumann) data in the Poisson solve and
exclude them from correction. Standard construction; est. 2-3 days.
Also note: double + RRTMGP + anelastic needs rad_ncol_chunk ~1024 and
arena ~3e9 on 16 GB (Kokkos OOM otherwise).

## 6g. RESOLVED IN FORK: the two defects behind the anelastic real-BC blowups

Fixed in 32e4d49f + d82162d9; anelastic full-physics 24-h hindcast now
completes in double precision (8,941 steps, wall 17.4 min, mass exact).

1. **Anelastic scalar advection used raw momenta as metric fluxes**
   (upstream-relevant, `ERF_SlowRhsPost.cpp`): the `l_anelastic` branch
   copied cur_{x,y,z}mom into avg_{x,y,z}mom, but AdvectionSrcForScalars
   consumes those as area-weighted, map-factor-scaled fluxes with OMEGA in
   the z slot (as AdvectionSrcForRho builds them in the compressible path).
   Feeding rho*w instead of Omega means scalars see a different effective
   divergence than rho over terrain slopes: measured exponential qv/theta
   pumping (e-fold ~65 s at 3 km over the San Gabriels / Channel Islands),
   independent of dt, PBL scheme, and microphysics phase; absent in the
   compressible A/B; invisible in flat idealized anelastic cases where
   Omega == rho*w. Fix: rebuild avg_* exactly as AdvectionSrcForRho
   (OmegaFromW; Omega(klo)=0).

2. **Anelastic radiation pressure from EOS of the frozen-rho state**
   (`ERF_Radiation.cpp`): with rho pinned to rho0, getPgivenRTh inherits
   base-state artifacts. The ERA5-interpolated rho0 is flat-extrapolated
   below the lowest ERA5 level, so EOS p_lay is flat/non-monotone over the
   first ~5 model levels (dp ~ 0.4 Pa in double -- RRTMGP survived by
   luck; exactly 0 in single precision -> col_gas = 0 -> non-finite
   optical depths on the first rad call). Fix: when anelastic, radiation
   uses the hydrostatic reference p0(z) and T = theta*Exner(p0) -- the
   thermodynamically correct anelastic pressure.

Also fixed en route: the prescribed-inflow mask (6f follow-on) must zero
correction fluxes ONLY on the domain-boundary normal faces (the only ones
the set_width=1 realbdy fill overwrites -- verified from realbdy_bc_bxs_xy
box arithmetic) plus terrain/lid; masking the second face layer traps
each boundary column's integrated convergence (no lateral relief) and
pumps scalars dJ-weighted into the deep top cells.

**Remaining (fork known-limit, not yet root-caused): single-precision
anelastic + RRTMGP** dies at ~3.5 model hours via a fast (~12-step) local
thermal collapse (transient 185-K cold pools at the stratospheric cold
point recover earlier in the run; onset is nondeterministic). Double
precision is unaffected (full 24 h clean); single precision with
radiation_model=None is also clean past the crash point (1300 steps).
So the defect lives in the SP radiation-anelastic interplay at the
clamped cold top.

## 11. `Stop datetime` startup line prints `start_datetime` (cosmetic, but misleading)

`Source/ERF.cpp` (in the datetime parse block) does:

```cpp
stop_time = static_cast<amrex::Real>(getEpochTime(stop_datetime, datetime_format));
Print() << "Stop  datetime : " << start_datetime << std::endl;   // <-- start_
```

`stop_time` itself is correct, so runs stop at the right instant; only the
log line is wrong. It matters because the startup banner is the natural
place to verify a segment's window before committing GPU-days to it, and
it currently reports every segment as zero-length. Fixed in this fork
(one identifier). Found 2026-07-24 while staging the Aug/Sep stress tests,
where the banner claimed `Stop datetime : 2023-08-20 00:00:00` for a run
that correctly integrated to 2023-08-21.

## 12. Hindcast band-density blending is not momentum-consistent

`hindcast_blend_band_density` (this fork, `ERF_Advance.cpp`) nudges `Rho`
in the lateral relaxation band toward the interpolated forecast and
rescales `RhoTheta` and `RhoQ1` so that theta and qv are preserved. It does
NOT rescale the momenta, which live in separate MultiFabs (`rU_old` etc.)
and have already been filled by the time the blend runs. Velocity is
therefore perturbed implicitly by `rho_old/rho_new` each step. The
perturbation is small (alpha <= 0.021/step) and the blend is required --
without it a calm segment inflates domain mass +33%/day -- but the
inconsistency is real and should be closed by scaling the momenta with the
same factor if this is ever upstreamed.

Not the cause of the August blowup investigated on 2026-07-24: that case
fails at cfl 0.3 with the blend OFF as well (nondeterministically -- see
item 13). Lowering cfl to 0.2 makes it stable with the blend on.

## 13. Run-to-run nondeterminism decides stability at marginal cfl (SP + GPU)

Two byte-identical invocations of the same single-precision GPU build, same
inputs, same `amrex.max_gpu_streams=1`, on the 2023-08-20 (Hurricane
Hilary) case at `erf.cfl=0.3`: one completed 250 steps clean, the other
produced whole-array non-finite RRTMGP optical depths at step 240. GPU
reduction/atomic ordering is not bit-reproducible, and in single precision
that is enough to decide whether a marginally-stable configuration
survives. Practical consequence for anyone tuning this fork: **a single
clean short run is not evidence of stability at the margin.** Require a
regime to run with zero w-damping events and zero low-temperature warnings,
not merely without NaNs. (The ChannelIslands production config now uses
cfl 0.2, which meets that bar on the year's most violent day.)

## 14. Lateral relaxation zone generates spurious updrafts, which precipitate out any moisture advected through it

**This is the most consequential defect found in the hindcast pathway so
far, and it invalidates precipitation output.** Measured 2026-07-24 on
three full 24-h runs (2023-08-20 Hilary, 2023-09-09, 2023-01-09) with the
locked production config.

### The artifact

At t=24 h, fraction of cells with w > 1 m/s, by distance from the lateral
boundary (cells), and cloud water in the same regions:

| run | band d<10 | mid 10-19 | interior d>=20 | band qc | interior qc |
|---|---|---|---|---|---|
| Hilary 08-20 | 11.62% | 1.38% | 0.82%  | 0.105 | 0.077 g/kg |
| storm 01-09  | 11.55% | 4.26% | 2.63%  | 0.311 | 0.253 g/kg |
| dry   09-09  |  9.18% | 0.56% | 0.28%  | 0.011 | 0.047 g/kg |

The band carries ~9-12% of cells in updraft in EVERY regime, including the
dry one -- 4x to 33x the interior fraction. w_rms in the band is 1.9-3.7x
interior. This is a property of the relaxation zone, not of the weather.

### The consequence

Precipitation is the product of that numerical updraft field and whatever
moisture is being advected through the boundary, so it appears only when
both are present:

| run | domain-mean 24-h precip | ERA5 same footprint | ratio | share of all precip falling in outer 10 cells |
|---|---|---|---|---|
| Hilary 08-20 | 181.9 mm | 16.4 mm | **11.1x** | **94.5%** |
| storm 01-09  | 153.8 mm |  5.9 mm | **26.1x** | **87.9%** |
| dry   09-09  |   7.5 mm |  0.0 mm | n/a       | 21.4% |

The dry case is the control that exonerates the microphysics: same
Morrison configuration, same spurious updrafts, but no moisture to condense
-- band rain rate DECAYS to 0.05 mm/h and ERA5 also reports zero. Morrison
is correctly raining out condensate that a numerical updraft produced.

### It grows with integration time

Band-mean precipitation rate (mm/h) through the Hilary run:

    t (h)     3     6     9    12    15    18    21    24
    band   1.87  2.29  2.85 10.27 25.32 24.20 33.40 36.25
    interior 0.25 0.25 0.21  0.19  0.22  0.37  0.67  1.08

Still accelerating at 24 h, and the interior is now rising too. **A
month-long segment will be far worse than these 24-h tests show.**

### Contamination radius

The excess decays inward but not to zero within the relaxation width.
Hilary 24-h rain by shell: d 0-4: 668 mm; 5-9: 118; 10-14: 28.1; 15-19:
17.6; 20-29: 9.6; 30+: 10.6. So `real_width = 10` is NOT a sufficient
analysis exclusion -- contamination is still 3x at d=10-19. Use d >= 20
(60 km at 3 km) as the minimum discard, and note that even there the
January interior runs 7.6x wetter than ERA5 (24.4 vs 3.2 mm).

Corollary: an apparent orographic precipitation signal is largely this
artifact. Interior corr(rain, terrain) collapses from +0.577 at d>=10 to
+0.105 at d>=20 for Hilary (January retains a real +0.380).

### Same root cause as the stability failures

The Hilary cfl-0.3 blowup was a w dipole reaching +48.9 m/s at (54,1,1..3)
-- in this same band. Lowering cfl to 0.2 keeps the artifact numerically
stable; it does not remove it. The unbounded band mass accumulation of
item 8 is very likely the same defect seen through a different variable.
`erf.hindcast_blend_band_density` constrains band DENSITY only -- it
explicitly rescales RhoQ1 to preserve qv, so it does nothing for this.

### Domain-geometry aggravation (this deck specifically)

The highest terrain in the ChannelIslands domain (1,377 m) sits at j=63,
i.e. **distance 0 cells from the boundary**. Terrain above 200 m is 8.6% of
the domain but only 0.4% of the d>=20 interior, where the tallest remaining
peak is 359 m. So every mountain is inside the contaminated zone and
orographic precipitation cannot be validated in this configuration at all.

### Mechanism (traced 2026-07-24, scoping only -- no fix attempted)

The relaxation specifies an incomplete state. Both the specified zone
(`ERF_BoundaryConditionsRealbdy.cpp`) and the RHS relaxation
(`realbdy_compute_interior_ghost_rhs`, `ERF_InteriorGhostCells.cpp`) act on
exactly four fields:

    cons_read = {0, 1, 0, 0, 1, 0, ...}   // Rho NOT read; RhoTheta, RhoQ1 read
    is_read: xvel 1, yvel 1, zvel 0       // w NOT read
    var_map  = {xvel, yvel, cons, cons}
    comp_map = {0, 0, RhoTheta_comp, RhoQ1_comp}

So u, v, theta and qv are forced toward ERA5 while **rho and w are never
constrained**. The imposed (rho, u) pair does not satisfy ERF's discrete
continuity equation -- the horizontal mass-flux divergence carried by the
imposed winds was never in balance with the local density -- and w is the
only remaining degree of freedom that can absorb the column imbalance.
Hence a persistent band updraft field that is present in every regime,
scales with wind speed rather than with moisture, and at cfl 0.3 ran away
to the +48.9 m/s dipole that NaN'd the Hilary case.

Difference from the WRF pathway, which shares this relaxation code:
`ERF_ReadFromWRFBdy.cpp` additionally reads PH (geopotential) and MU (dry
column mass) and performs an explicit mass coupling at read time
(`mu = mu_arr + mub_arr`, with PH/PHB) to re-derive boundary column heights
and re-interpolate U/V/T/QV onto them. The hindcast pathway has no
analogue: `RealBdyVars` is only {U, V, T, QV}. Whether wrfbdy-driven runs
are actually clean has NOT been tested here (no wrfinput/wrfbdy available
since the WPS chain was dropped) -- but the WRF path has a consistency step
this one entirely lacks.

**The needed data is already present and discarded.** The ERA5 .bin frames
carry 8 fields [rho, u, v, w, theta, qv, qc, qr];
`FillForecastStateMultiFabs` interpolates all of them onto the ERF grid;
`fill_bdy_data_from_hindcast` then copies only u, v, theta, qv into the
boundary planes. rho and w are computed and thrown away.

## 15. SP + AMR: t_old = time - Real(1.e200) overflows float to -inf and
## ErrorEst's FillPatchCrseLevel then wipes the just-initialized level 0

In a single-precision build, `Real(1.e200)` is +inf, so the "never been
advanced" sentinel `t_old[lev] = time - Real(1.e200)` (ERF_MakeNewLevel.cpp,
ERF.cpp init_only, WRFInput/Metgrid init) sets t_old = -inf.
`amrex::almostEqual(time, -inf)` is TRUE for any finite time (|x-y| = inf
<= eps*|x+y| = inf), so `FillPatchCrseLevel` takes the `time == t_old`
branch and sources `vars_old` -- which at init is never-written define
memory. Single-level runs never hit this because the only pre-first-step
caller of FillPatchCrseLevel on the full state is `ErrorEst`, which runs
only when `max_level > 0`. Net effect: with amr.max_level=1 (even with no
refinement box at all), the freshly HindCast-initialized level-0 state is
replaced by garbage (rho ~ 0) during InitFromScratch; radiation then
aborts on all-NaN gas-optics inputs while the dycore itself limps along.
Diagnosed with MASS SL/ML = 0 / (nest-fraction) at t=0.
FIX (this fork): use a finite sentinel `Real(1.e30)` at all five sites.
DP builds are unaffected (1.e200 is finite in double).

## 16. SP + AMR: ERFFillPatcher::Fill time-window assert uses an absolute
## float-epsilon; subcycled fine times legitimately drift a few ulp of t

`ERF_FillPatcher.H` asserted `time` within `[m_crse_times[0]-eps,
m_crse_times[1]+eps]` with `eps = numeric_limits<float>::epsilon()`
(absolute, ~1.19e-7). In SP, one ulp of elapsed time t is 1.19e-7*t --
already 4e-7 at t = 3.35 s, 8e-3 s at t = 86400 s. Fine-level and
coarse-level clocks accumulate a few ulp apart under subcycling, so the
assert fires within seconds of model time (observed: level-1 step 36,
t = 3.2032070 vs window edge 3.2032068). FIX (this fork): tolerance made
relative (8 ulp of the window magnitude) and the time-interpolation
factors clamped to [0,1]. DP hits the same issue only after ~1e9 s.

## 17. Nested domains have NO active boundary relaxation -- the nest is a
## sealed cavity and rings itself to NaN (likely the mechanism behind #3135)

Upstream ERF's live nesting support fills fine-level ghosts by interpolation
and hard-sets interface-normal momenta (ERFFillPatcher FillSet), with no
blending zone: `cf_set_width != 0` is an unconditional Abort
(Source/ERF.cpp:2970), `cf_width > 0` only triggers RegisterCoarseData whose
cons data is never consumed, and `fine_compute_interior_ghost_rhs`
(Source/Utils/ERF_InteriorGhostCells.cpp) -- a WRF-style relax zone over all
IntVars INCLUDING rho and zmom -- has zero call sites. Measured consequence
(ChannelIslands hindcast, 3 km parent + 1 km ocean nest): a coherent
whole-nest w oscillation (w_rms identical at every shell d=0,1,2,5), growing
and collapsing to NaN at t~216 s under BOTH compressible and anelastic nest
solvers -- i.e. not an acoustic-scheme problem but wave energy trapped by a
reflective interface. FIX (this fork): wired fine_compute_interior_ghost_rhs
into the slow RHS for level>0 under cf_width>0 (one-way coupling), after
correcting bit-rot: RegisterCoarseData registers MOMENTA, but the dead
routine multiplied the FillRelax result by rho again (rho^2*u). With
cf_width=10 the same nest runs 3500 steps with zero warnings and the nest
rim shows w>1 fractions of 1.3-2.7% (= interior background), vs 12-24% for
the lev-0 ERA5 band. Box-aligned stripe artifacts at nest boundaries over
terrain (#3135, closed without diagnosis) are consistent with this
reflective-interface mechanism.

### #17 addendum: measured behavior of five c/f interface configurations
### under a strong winter jet (Jan-9 2023, 3 km parent / 1 km nest, SP GPU)

The convective regime (Sep-9 max-CAPE) is STABLE and rim-clean with
cf_width=10 relaxation (3,500 steps, zero warnings). The winter-jet regime
(40+ m/s crossing the nest) fails under every interface tried, each in a
different, diagnosable way:

1. relax incl. rho (cf_width=5):    NaN ~330 fine steps  (acoustic
   destabilization by mass-equation RHS relaxation -- same failure class
   the lev-0 band hit, which is why lev 0 uses a post-step rho blend)
2. relax incl. rho (cf_width=10):   NaN ~600 (wider band delays, same mode)
3. relax excl. rho:                 band rho ratchets 1.2 -> 2.5 kg/m3 in
   ~2.5 min (continuity violation from momentum forcing with no mass
   closure), EOS/radiation NaN ~450
4. + post-step band rho blend toward parent (tau ~ 20 s), weaker momentum
   relax (tau 50*dt), cf_width=20:  slows pile-up 4x, NaN ~1140
5. + WRF-style specified zone (cf_set_width=3, fork-enabled past the
   upstream abort; FillSet of full cons incl rho): mass budget contained
   (rho max 1.43 vs 2.5), cold-top warnings vanish, but a near-surface
   seam instability forms at the relax-band INNER edge (fine cell 20 of a
   20-cell band, upwind side) and NaNs ~720.

Interpretation: ERF currently has no complete nest-boundary closure; the
pieces (FillSet, masks, relax RHS) exist but were never finished. The
specified-zone + relax combination is closest -- the remaining defect is
the band-to-interior transition (linear Factor taper; WRF uses exponential
decay and applies horizontal diffusion in the relax zone).

### #17 second addendum: the nested instability is STOCHASTIC, not
### configuration-bound (2026-07-25)

Re-running the exact configuration that survived 3,500 coarse steps on the
Sep-9 convective case (same binary lineage, same knobs: cf_width=10 relax
incl. rho, cfl 0.15, original terrain) produced deaths at coarse steps 151,
181 (with amrex.max_gpu_streams=1 -- stream race excluded), and 211 across
repeats; the Jan-9 jet case dies at 330-1140 depending on interface tuning.
Failure signature is always the same: the level-1 state goes wholesale NaN
within a few fine steps between interval checks, first observed by the
level-1 radiation input marshaling. Per-step NaN checking (which inserts
device syncs every step) extends survival (~600-1800 fine steps) but does
not prevent death. The one 3,500-step clean run was a lucky draw of the
same SP-GPU run-to-run nondeterminism documented in #13.

Status: the two-level HindCast path (this fork's wiring of upstream's
incomplete nest machinery) is not production-viable in single precision on
GPU. Next diagnostics, in order of information value: (1) compute-sanitizer
initcheck/racecheck over ~20 steps of the 2-level case; (2) a
double-precision control build (DP stable => SP conditioning of the c/f
interpolation/relaxation; DP unstable => algorithmic defect in the nest
path); (3) upstream escalation with this ladder -- upstream CI has no SP
GPU multilevel real-case coverage that would have caught any of #15-#17.

---

## 18. Lateral relaxation manufactures spurious vertical motion: the ramp-gradient term

**Severity:** invalidates precipitation and vertical-motion statistics in and near
the relaxation zone for every `use_real_bcs` run. Affects wrfbdy/metgrid equally;
found on `init_type=HindCast`.

### Symptom
In the level-0 relaxation band, 9-50 % of cells carry |w| > 1 m/s against an
interior background near 11 %, in every regime tested. It is not orographic and
not inherited from the driving data (both shown below).

### Derivation (operator level)
`Omega` is hard-zeroed at `k=klo` (SlowRhsPre, SlowRhsPost, Substep_T,
PoissonSolve) and the lid is a SlipWall, so the vertical mass flux telescopes
out of a column sum of the discrete continuity equation in
`AdvectionSrcForRho`:

    d/dt sum_k (detJ/m^2) rho = - sum_k [ Dx(ax rho_u/mf_u)/dx + Dy(ay rho_v/mf_v)/dy ]

The column mass tendency is determined ENTIRELY by the horizontal fluxes -- this
is ERF's implicit analogue of WRF's mu equation. (ERF does not reproduce WRF's
mu-coupled lateral BC: `WRFBdyVars::MU` is used only to de-couple at read time in
`convert_wrfbdy_data` and is never consulted again; `PC` is never used at all.)

Prescribing rho_u and rho_v in the band therefore already prescribes the column
mass tendency, while `Rho_comp` is relaxed toward nothing at all -- `comp_map` in
`realbdy_compute_interior_ghost_rhs` covers U,V,T,QV only. The level-by-level
residual has only `Dz(rho Omega)` left to absorb it, so **Omega, and hence w, is
the residual variable.**

The concrete source is the ramp. The relaxation adds `S = F*F1*(A-B)` to the
MOMENTUM rhs, and

    div(F*A + (1-F)*B) = F div A + (1-F) div B + grad(F).(A-B)

The last term has no physical counterpart. It peaks where |grad F| * |momentum
error| is largest -- mid-band, since the error grows as the forcing weakens.

### Evidence
1. **The driving data is clean.** A discrete continuity residual of the
   interpolated ERA5 target, measured with ERF's own operator, gives an implied
   |w| of 6.4 mm/s in the band and 5.2 mm/s in the interior -- 150x below the
   1 m/s criterion, and the same inside the band as outside. The target also
   equals the initialized state exactly (max|drho| = 0).
2. **The wall is clean; the mid-band is not.** Decomposing by wall, over flat
   water (12 m terrain) the wall itself is 0.00 % at d=1-2 while d=5-7 reaches
   44-50 %. The only loud cells at d=0 are on the two land walls (250 m terrain).
3. **The scaling law holds.** Peak location tracks `real_width` and peak
   amplitude tracks 1/width, to ~10 % across three widths:

   | real_width | peak location | peak \|w\|>1 | predicted | mass drift |
   |---|---|---|---|---|
   | 10 | d = 6  | 49.9 % | --     | +0.53 %/day |
   | 15 | d = 10 | 36.1 % | 33.3 % | +0.92 %/day |
   | 20 | d = 15 | 21.6 % | 25.0 % | +3.82 %/day |

4. **Weakening the nudge makes it worse.** `bdy_nudge_factor` 10 -> 50 (larger
   momentum error) raises it to 38.9 % right at the wall.

### Three fixes attempted, all reported honestly
- **(a) Complete, mass-consistent relaxation target** (relax `Rho_comp` toward
  rho*, momentum target rho* u* instead of rho_model u*). Mass drift bounded,
  but near-wall w_rms went 0.71 -> 4.33 m/s with the excess growing monotonically
  with height (0.04 m/s at k=0 to 8.2 m/s at k=28). A relaxation source in the
  MASS equation forces the acoustic solver. Knob `erf.hindcast_mass_consistent_bdy`,
  default off. Prescribing rho* u* WITHOUT the rho constraint NaNs at step 1 --
  the momentum target must remain velocity-stabilising (rho_model u*).
- **(b) Barotropic wall mass-flux correction** delivered through the boundary
  FLUX, never a volumetric source: `dvel = dx*(M_tgt/M - 1)/tau` on the wall
  face. Dynamically invisible (shell profile identical to control to two
  decimals at every d) and halves mass drift, 0.53 -> 0.28 %/day, with no blend.
  But it cannot touch the band w, because the defect is not at the wall.
  Knob `erf.hindcast_wall_flux_correction` (+`_tau`), default off. Worth having
  for mass; tau=900 s destabilises, tau=3600 s is stable.
- **(c) Cancelling grad(F).(A-B) with a companion MOMENTUM increment** (no mass
  source), 1-D wall-normal since grad F is wall-normal:
  `C(d+1) = C(d) + 2*F1*(w-d)/w^2 * err(d)`. **Both disposals of the integration
  residual fail, for one reason:** the net is weighted by (w-d)/w^2 times an
  error that GROWS inward, so it peaks mid-band exactly where the artifact does
  and is comparable to it.
  - taper C to 0 at both ends -> peaks 47.0/32.0/15.8 % at width 10/15/20, still
    tracking 1/width; near-wall got worse.
  - route the residual out through the wall -> globally destabilising, interior
    background 11 % -> 32-56 %.
  Knob `erf.hindcast_ramp_div_correction` (0 off / 1 taper / 2 wall), default 0.

### Conclusion
No LOCAL cancellation exists: the wall-normal integral of grad(F).(A-B) has a
nonzero net that is as large as the artifact, and it cannot be disposed of
either inside the band or through the boundary. The only untried route is an
exact 2-D Helmholtz projection of the correction restricted to the band, which
would leave only the global compatibility constant. Absent that, the practical
mitigation is a WIDER, gentler ramp -- which works (49.9 -> 21.6 %) but costs
domain area and worsens mass drift 7x, partially recoverable with (b).

### Two independent defects found along the way
- **Hi-wall product-rule RHS reads a permanently-zero ghost.** In the normal-face
  block (PR #3209), `ihi = domainx.bigEnd(0) = nx`, so `rhs_cons(ihi,j,k)` is a
  ghost of `F_slow[cons]`: never written (every RHS writer uses `mfi.tilebox()`)
  and `setVal(0)` once at construction. The `u*d(rho)/dt` term is therefore live
  on the lo walls and silently dropped on the hi walls. Fixed by reading the
  adjacent valid cell.
- **The HindCast frame `w` slot is ERA5 omega in Pa/s, not m/s.** Measured from
  the frame binary: range -1.8..+2.5, mean|.| 0.094 -- consistent with Pa/s, not
  with a 25-km geometric vertical velocity. Any consumer of that slot needs the
  omega -> w conversion first. Also `forecast_state_interp[lev][Vars::zvel]` is
  never written (the LinComb is commented out), so no w target exists on the ERF
  grid at all. This is the trap the abandoned lateral WfromOmega experiment
  (PR #2872) walked into.

### Ruled out: the relaxation's placement in the time splitting

ERF is the outlier among split-explicit LAMs in *where* it applies the Davies
relaxation. ERF adds it to `S_rhs` in `slow_rhs_pre`; `S_rhs` is then held fixed
while every acoustic substep integrates it. COSMO instead excludes the
relaxation from the slow-mode forcing and applies it Marchuk-split once per big
step, explicitly "for stability reasons" (Leps, Brauch & Ahrens 2019, JAMES 11,
2694-2707); TRAM (Romero et al. 2024, QJRMS, doi:10.1002/qj.4639) likewise
applies its Newtonian relaxation after each completed time step.

Tested with `erf.realbdy_relax_split` (default 0): the identical operator, same
ramp, same coefficient, same target, moved out of `S_rhs` and applied once per
timestep as a sequential state update after the substep loop. Net per-step
amplitude is unchanged -- the RHS path uses `F1 = 1/(nudge*slow_dt)` with the
*stage* dt and RK3 restarts each stage from `S_old`, so both paths deliver
`Factor/nudge * (A - B)` per step.

Sep-9 dry control, 1200 steps, blend off, wall flux correction on. Ocean-wall
|w| > 1 fraction (the two flat-water walls, 12 m terrain), and mass drift:

| real_width | ordering | peak (d) | peak value | interior bg | mass drift |
|---|---|---|---|---|---|
| 10 | in `S_rhs` | 6 | 49.85% | 10.88% | +0.453 %/day |
| 10 | Marchuk split | 6 | 50.17% | 10.86% | +0.382 %/day |
| 15 | in `S_rhs` | 10 | 36.07% | 11.72% | +1.469 %/day |
| 15 | Marchuk split | 10 | 36.28% | 11.68% | +1.476 %/day |

Runs are deterministic here (repeats agree to 0.02 points / 3 decimals), so the
differences are real and they are nil: the peak does not move, its amplitude
rises 0.6% relative, and the 1/width scaling is untouched (49.85/36.07 = 1.38
either way). The mass-drift gain at width 10 (-16%) does not reproduce at width
15 (+0.5%). **The ordering is a genuine structural difference from COSMO/TRAM
and is worth fixing on its own terms, but it is not the source of the band
artifact.**

### Control: the artifact is entirely relaxation-sourced, and it is not confined to the band

Same deck with `erf.bdy_nudge_factor = 1e7` (relaxation effectively off; the
specified-zone Dirichlet wall and the wall flux correction still active):

| d | relaxation on | relaxation off |
|---|---|---|
| 0 | 0.20% | **55.48%** |
| 1-4 | 0.00 / 4.99 / 19.13 / 33.62% | 7.60 / 7.68 / 3.84 / 5.34% |
| 5-9 | 44.10 / **49.85** / 46.41 / 21.44 / 8.71% | 0.24 / 0.69 / 0.02 / 0.06 / 0.06% |
| interior bg (d>=25) | **10.88%** | **2.86%** |
| mass drift | +0.453 %/day | +5.346 %/day |

Two things follow. (1) The mid-band peak is 100% attributable to the relaxation
-- with the ramp gone it does not merely shrink, it vanishes into the noise
floor. (2) The relaxation also raises the *whole-domain interior* background by
3.8x, from 2.86% to 10.88%, 75 km from the nearest wall. At 1200 steps
(31 model minutes) advection cannot carry the band signal that far (19 km at
10 m/s), but acoustic (220 s) and gravity-wave (~2500 s) propagation can
[inferred]. The `grad(F).(A-B)` source is therefore contaminating the free
interior, not just the sponge.

The cost of switching it off is the reason Davies relaxation exists: the bare
Dirichlet wall goes to 55% at d=0 and mass drift worsens 12x. This is the
over-specification tradeoff, measured -- not an argument for keeping the ramp.

### The interior contamination is acoustic, and it is source-limited not reflection-limited

Domain 128x64x32 at 3 km = 384 x 192 km. `d >= 25` is >= 75 km from the nearest
wall. At 1200 steps, t = 1869 s. Carrier arrival times over 75 km:

| carrier | speed | time to 75 km | arrived by t=1869 s? |
|---|---|---|---|
| acoustic | c = sqrt(gamma R T) = 335 m/s at 280 K | 224 s | yes, 8.3x over |
| gravity wave | c* ~ 25-30 m/s (the speed WRF/ERF radiation BCs use) | 2500-3000 s | **no** |
| advection | ~10 m/s | 7500 s | no |

So the 10.88% interior background at 1869 s can only have been carried by sound.
Extending the control to 3600 steps (t = 5649 s), by which time gravity waves
*have* arrived, separates the two:

| t (s) | bg (d>=25), relax on | bg, relax off |
|---|---|---|
| 604 | 5.60% | 2.57% |
| 1233 | **36.76%** | 1.99% |
| 1869 | 10.88% | 2.86% |
| 3793 | 11.27% | -- |
| 5649 | 11.43% | -- |

The level is set by 1869 s and moves 0.55 points over the next 3780 s, so the
gravity-wave share of the deep-interior contamination is <= 5%: **it is ~95%
acoustic, measured rather than inferred.** The 36.76% transient at t = 1233 s
is within 8% of the 1146 s acoustic transit of the long axis, consistent with a
single coherent pulse sweep [inferred].

That raises the question of whether the acoustically rigid lateral wall (`open`
== `foextrap`, ghost rho_theta = first interior value, dp/dn = 0, plus `gpx`
forced to zero on the real-BC path) is trapping the radiated energy. Two
independent lines say no:

1. **Saturation.** The interior level is flat from 1.6 to 4.9 long-axis acoustic
   transits. A high-Q cavity would keep filling. It does not.
2. **Direct test.** The existing lateral Rayleigh layer (`erf.hindcast_lateral_
   sponge_strength`, 30 km deep) at 0.003 and 0.01 /s. At 0.01 /s a wave takes
   90 s to traverse the layer, so it is damped by exp(-0.9) ~ 0.41 per pass,
   ~60% energy absorption. Interior background: 10.88% (off) -> 10.86% (0.003)
   -> 10.77% (0.01). **Nothing.** Caveat: that sponge damps rho_u toward the
   *initial* state with its own xi^2 ramp, so it adds a source while absorbing,
   and being velocity-only it is impedance-mismatched and partially reflecting.
   A null result alone would not be conclusive; combined with (1) it is.

**Conclusion: the band radiates directly inward and that first pass sets the
interior level. A non-reflecting lateral treatment cannot reach the 2.86% floor
-- only removing the ramp source can.** The NSCBC outflow radiation condition
remains a well-posedness fix worth making, but it is not the lever for this
metric.

### Kill condition for the Helmholtz-projection fix: FIRES

The proposed fix replaces the relaxation forcing C_raw = F*(A-B) with its
divergence-free part, C_proj = C_raw - grad(phi), lap(phi) = div(C_raw), so that
grad(F).(A-B) vanishes by construction. C_raw was dumped from the control
trajectory at step 1200 (`erf.realbdy_dump_relax_step`) and decomposed offline
per level on the band annulus, phi = 0 outside the band (exact, since C_raw = 0
there -- verified: energy outside the band is identically 0) and homogeneous
Neumann at the wall so the wall flux is left to the wall flux correction.
Divergence is removed to 1e-10 and `retained` reproduces `ratio^2` to 0.3%,
confirming a true orthogonal projection.

| real_width | \|\|C_proj\|\|/\|\|C_raw\|\| | energy retained | curl-free share |
|---|---|---|---|
| 10 | **0.376** | 14.2% | **85.8%** |
| 15 | **0.250** | 6.3% | **93.8%** |

**86-94% of the relaxation forcing is a pure gradient.** The reason is
structural, not incidental: C_raw = F(n)*(A-B) with F a monotone function of
wall-normal distance alone. For a locally uniform error (A-B) = (a,0) the
forcing is C = (F(x)a, 0) = grad(a * integral F dx) -- *exactly* a gradient. The
ramp IS the gradient content. And it gets worse as the ramp gets smoother
(0.376 -> 0.250 from width 10 to 15), the same 1/width signature seen from the
other side.

Projecting therefore destroys 86-94% of the nudging, and specifically destroys
the divergent part -- which is precisely what mass drift is made of. The
surviving rotational remnant carries no mass constraint at all, so mass drift
would be at or worse than the relax-off +5.346 %/day. The streamfunction
fallback (solve lap(psi) = curl C_raw, build C = (-d psi/dy, d psi/dx)) dies on
the same number: it reconstructs the rotational Helmholtz component, which *is*
the 14%.

**The projection approach is not viable. Nothing was built.**

### ERF's real-BC default already gets the one thing that matters right

The characteristic analysis says exactly one variable must be left FREE at a
subsonic lateral inflow: rho (equivalently p), because it rides the OUTGOING
u_n - c wave and must be computed from the interior. Everything else -- u_n,
both tangential velocities, theta, and every scalar -- is admissible.

`fill_from_realbdy` already does this. `cons_read[Rho_comp] = 0`, so rho is
zero-gradient (computed from the interior), and the driver's theta and q are
multiplied by that *local* rho (`dest_arr(...,comp) *= dest_arr(...,Rho_comp)`).
So theta is specified and rho is free -- precisely the admissible split.

That reframes two proposals. Attempt 2 here (`erf.hindcast_mass_consistent_bdy`)
and upstream PR #3483 both set rho at the boundary from the driving deck. Both
therefore *break a behaviour ERF already gets right*, and the characteristic
analysis says so independently of the empirical failure we measured (near-wall
w_rms 0.71 -> 4.33 m/s, excess growing monotonically with height). The
recommendation upstream is not "add rho at the boundary" -- it is "do not".

The accurate accounting of what ERF specifies today, for a Morrison deck
(N = 13 scalars, 18 variables, admissible 17 in / 1 out):

| variable | treatment | specified? |
|---|---|---|
| rho | zero-gradient | **free (correct)** |
| theta | driver | yes |
| u, v | driver | yes x2 |
| w | zero-gradient | free |
| q_v | driver | yes |
| q_c .. q_11 | **zeroed** | yes x10, to a wrong value |
| KE, scalar | zero-gradient | free x2 |

14 specified, 4 free. So ERF is **under**-specified by 3 at inflow (w, KE,
scalar float when they must be given) and **over**-specified by 13 at outflow.
Not the "over-specified everywhere" picture; the inflow side is nearly right.

### Step 3 result: NSCBC eliminates the artifact, and breaks the mass budget

Built behind `erf.nscbc_lateral` (default 0, so every non-real-BC path is
byte-identical). Regime from the sign of the just-prescribed driver normal
velocity, blended over a +-0.5 m/s window; inflow specifies 4+N with rho free;
outflow computes from the interior with either plain extrapolation
(`erf.nscbc_outflow=0`) or one incoming-acoustic Riemann condition
(`=1`); the ramp is removed. Sep-9 dry control, 1200 steps, blend off, wall
flux correction on.

| config | d=0 | band peak (d>=3) | interior bg | mass drift |
|---|---|---|---|---|
| control (Davies) | 0.20% | **49.85%** (d=6) | 10.88% | +0.453 %/day |
| relax-off floor | 55.48% | 0.69% | 2.86% | +5.346 %/day |
| NSCBC, extrap outflow | 49.13% | **0.00%** (d=3-11) | **2.46%** | +174.6 %/day |
| NSCBC, Riemann outflow | 61.42% | 0.97% | **1.50%** | +96.0 %/day |
| NSCBC + ramp kept, w=10 | 8.43% | 50.51% (d=6) | 15.22% | +637 %/day |
| NSCBC + ramp kept, w=15 | 4.83% | 34.83% (d=10) | 7.66% | +559 %/day |

Two results and one failure.

1. **The band artifact is eliminated, not reduced.** With the ramp gone the
   profile is 0.00% from d=3 to d=11. And the keep-ramp rows show the 1/width
   scaling returning the moment the ramp is restored (50.51% at width 10,
   34.83% at width 15, ratio 1.45 against 1.5 predicted). The mechanism is
   confirmed by construction and by ablation.
2. **The interior background beats the relax-off floor**: 1.50% against 2.86%,
   7x better than the Davies control's 10.88%.
3. **The mass budget is destroyed** -- +96 to +175 %/day, 200-400x the control.

The mass failure is isolated, not speculative. Bisecting the pass with
`erf.nscbc_parts` (bit 2 = the velocity pass):

| parts | wall velocities | mass drift |
|---|---|---|
| 13 (no velocity pass) | driver, as today | **-18.3 %/day** |
| 7 / 15 (velocity pass on) | characteristic | **+94.7 %/day** |

and retuning the wall flux correction does almost nothing (tau 3600 -> 100, a
36x faster correction, moves drift only 96.0 -> 86.4 %/day).

**Cause: ERA5's boundary winds are approximately mass-balanced around the
domain, and the characteristic outflow condition replaces half of them with
model-derived values, destroying that balance.** This is the classical
limited-area solvability problem, and it is also why d=0 stays loud -- the same
defect seen locally rather than globally. The characteristic count tells you how
many conditions are admissible; it says nothing about the discrete global mass
budget, and in a limited-area domain that has to be imposed separately.

The remedy is a global outflow rescaling of exactly the kind ERF already
implements for the anelastic path: `enforceInOutSolvability`
(`ERF_PoissonSolve.cpp:465`), which scales outflow to match inflow. Applied to
the NSCBC wall fluxes with the target net flux set from the ERA5 column-mass
tendency, it is the missing third piece. Not yet built.

### Global mass constraint: fixes the budget, does NOT fix the wall

Added `erf.nscbc_mass_tau`. The total wall mass flux is corrected to the value the
ERA5 column-mass tendency asks for -- standard practice for limited-area models
that use characteristic/radiation OBCs (Flather 1976 is derived from mass
conservation; ROMS/NEMO adjust barotropic inflow-outflow to preserve volume), and
the same thing ERF already does for the anelastic path in `enforceInOutSolvability`.
The target is trustworthy because the ERA5 field is mass-consistent under ERF's own
discrete operator to an equivalent |w| of 6.4 mm/s.

    Phi_desired = (M - M_tgt)/tau ;  du = (Phi_desired - Phi_now)/sum_outflow(rho*A)
    u_n -> u_n + du   on outflow faces only, du a single uniform scalar

**The correction must be additive, not multiplicative.** Rescaling outflow by
lambda = (Phi_desired - Phi_in)/Phi_out was tried first: with Phi_desired ~ 0 it
demands lambda = -Phi_in/Phi_out, which in a net-convergent synoptic regime is
many-fold. It hit the guard rail, collapsed the CFL and killed every run inside
30-180 steps. The additive form is bounded by the imbalance divided by the wall's
mass-flux capacity, ~0.6 m/s here.

| config | d=0 | band d=3-10 | interior bg | mass drift |
|---|---|---|---|---|
| control (Davies) | 0.20% | **49.85%** | 10.88% | +0.453 %/day |
| relax-off floor | 55.48% | 0.69% | 2.86% | +5.346 %/day |
| NSCBC extrap, no constraint | 49.13% | 0.00% | 2.46% | +174.6 %/day |
| **NSCBC extrap + constraint, tau=600** | **46.12%** | 0.34% | 3.55% | **-7.26 %/day** |
| NSCBC extrap + constraint, tau=1800 | 46.14% | 0.37% | 3.60% | -13.34 %/day |
| NSCBC Riemann + constraint, tau=600 | 61.75% | 2.07% | 2.52% | +65.75 %/day |

Three things.

1. **The budget improves 24x** for the extrapolation outflow (+174.6 -> -7.26 %/day)
   and is now bounded rather than runaway. It is still 16x worse than the Davies
   control's +0.453 %/day. The residual is most likely metric error in the flux sum
   (rho taken at the wall CELL not the face; uniform reference dz against a stretched
   grid with dz_min = 18.5 m) -- the constraint drives the flux *it computes* to the
   target exactly, so whatever it mismeasures shows up as drift.
2. **No band structure reappeared.** d=3-10 stays at 0.14-0.79%, so the uniform
   additive correction did not rebuild a wall-normal gradient. There is a mild
   shoulder at d=11-16 (2.5-4.7% against nsc0's 0.06-3.4%) which is diffuse and does
   not have the monotone ramp signature; it is not a returning grad(F).
3. **The d=0 prediction is FALSIFIED, and this is a separate mechanism.**
   The prediction on record was that the loud wall cell is the local face of the
   global imbalance, so fixing the budget should relieve it. The budget improved 24x
   and d=0 moved 49.13% -> 46.12% -- three points, i.e. not at all. For the Riemann
   variant it went the wrong way, 61.42% -> 61.75%.

**The loud wall cell is therefore an unexplained defect independent of the global
mass budget.** It is not the ramp (removing the ramp is what exposed it), not the
global budget (fixed, no effect), and not the w specification (measured: the
ERA5-consistent w is 6.4 mm/s, and the MPAS w=0 A/B moved d=0 by 0.25 points).
The remaining untested candidate is that rho at the boundary is only ZEROTH-order
extrapolated: the derivation says it must be computed from the interior via the
outgoing u_n - c characteristic, and zero-gradient is the crudest possible stand-in
for that. That is the next thing to try, and it is stated here as a hypothesis, not
a conclusion.

### Characteristic inflow density: the coupling was real, and mass drift now lands

Two fixes, deliberately separated with `erf.nscbc_parts` bit 16.

**Flux metric (rho at the face, not the wall cell): no effect** -- -7.265 -> -7.297
%/day. And the dz half of that concern was unfounded: `ax` already carries the
vertical stretching (`ax = 0.5*(z_nd(k+1)-z_nd(k))/dz_ref`), so `ax*dy*dz_ref` is
already the true face area. Only the rho location was ever wrong, and it was worth
0.03 %/day.

**Characteristic inflow density: real.** At a subsonic lateral boundary, in the
outward-normal frame, `u_n + c` is OUTGOING in BOTH regimes (at inflow u_n < 0 but
|u_n| < c). So `J+ = u_n + 2c/(gamma-1)` carries the interior state out and fixes
the boundary sound speed once u_n is specified. With theta also specified, rho
follows as a ratio against the interior with no constants needed:

    c_b = (gamma-1)/2 * (J+_interior - u_n,driver)
    rho_b/rho_i = (th_i/th_b) * [ (c_b/c_i)^2 (th_i/th_b) ]^(1/(gamma-1))

and every rho-weighted quantity is rescaled onto it so theta and each q keep their
specified values. Zero-gradient rho -- what ERF does today -- is the crudest
possible stand-in, and it puts a density discontinuity on the wall face.

Drift -7.30 -> -5.48 %/day at tau=600 from that change alone. Drift then scales
with tau, which is the signature of a static bias offset rather than ongoing loss:
the constraint settles where `(M - M_tgt)/tau` balances the flux-measurement bias,
so the offset is proportional to tau and the measured "drift" is the approach to it.

| tau (s) | drift, whole run | drift, second half |
|---|---|---|
| 300 | -3.153 %/day | +0.552 %/day |
| 100 | -1.363 %/day | +0.132 %/day |
| **60** | **-0.917 %/day** | **-0.079 %/day** |
| Davies control | +0.453 %/day | -6.019 %/day |
| relax-off floor | +5.346 %/day | -4.229 %/day |

Note the Davies control's own mass trace is non-monotonic (+0.453 whole-run,
-6.019 over the second half), so which window you score on matters. On the settled
rate the NSCBC path at tau=60 is -0.079 %/day, two orders better than the control;
on the whole-run number it is 2x the control's magnitude with the opposite sign,
because it includes the constraint's startup transient.

Full state at tau=60, against the two reference points:

| metric | Davies control | relax-off floor | NSCBC + constraint |
|---|---|---|---|
| d=0 | 0.20% | 55.48% | 42.75% |
| band d=3-10 | **49.85%** | 0.69% | **0.09-0.27%** |
| interior bg (d>=25) | 10.88% | 2.86% | 3.61% |
| drift (settled) | -6.019 %/day | -4.229 %/day | **-0.079 %/day** |

The characteristic rho also moved the near-wall cells that the global constraint
could not: d=0 49.13% -> 42.75%, d=1 7.57% -> 1.25%, d=2 6.40% -> 3.74%. So the
coupling was real -- the same zeroth-order rho was hurting both the wall and the
budget. d=0 at ~43% remains unexplained and is NOT claimed as fixed.

The d=11-16 shoulder did **not** sharpen (1.32/2.78/4.43/4.17/4.79/4.33% against
the previous 2.49/3.82/4.74/4.14/4.72/4.04%) -- it is unchanged, diffuse, and still
carries no monotone ramp signature.

**Unexplained and not chased, per instruction:** the Riemann outflow variant
(`erf.nscbc_outflow=1`) reversed under the global constraint. Without the
constraint it halved the drift versus plain extrapolation (+96.0 vs +174.6 %/day);
with it, it is far worse (+65.8 vs -7.3 %/day). The extrapolation variant is used
throughout.

### Jan-9 wet test: the band precipitation artifact is gone; the interior is wetter

24 h, Jan-9 2023, both on the same binary. NSCBC = `nscbc_lateral=1 nscbc_outflow=0
nscbc_parts=31 nscbc_mass_tau=60`. Both ran clean: 24 h complete, exit 0, zero
w-damping / low-temperature / negative-theta warnings.

**Reference correction.** ERA5 over the *whole grib footprint* is 13.46 mm, but that
footprint is much larger than the ERF domain and extends into wetter terrain. Over
the matching ERF footprint it is **5.91 mm** (d>=0), **3.17 mm** (d>=20), peak
33.1 mm. Every ratio below uses the matched footprint. Note also that ERA5 at 25 km
peaks at 33.1 mm against ~130 mm observed, so it under-resolves the maximum by ~4x
and is a poor denominator for peak comparison.

**The historical 26x does not reproduce.** Davies on the current binary gives 5.97x
(matched footprint) or 2.61x (grib footprint) -- not 26x. Whatever that number
measured, it is not what this deck and this code now do, and the like-for-like
control below is the number to use.

| metric | ERA5 (matched) / obs | Davies | NSCBC |
|---|---|---|---|
| domain mean | 5.91 mm | 35.30 mm (5.97x) | 30.36 mm (5.14x) |
| d>=20 mean | 3.17 mm | 13.61 mm (4.29x) | 22.77 mm (7.18x) |
| d>=20 max | ~130 mm obs | 222.1 mm (1.71x) | 244.9 mm (1.88x) |
| mass drift, settled | -- | +0.775 %/day | **-0.098 %/day** |
| band w, d=6 (ocean walls) | -- | **15.5%** | **0.5%** |

**The decisive result is the per-wall precipitation.** 24-h wall-normal winds:
xlo +10.53, ylo +12.65 m/s inward (inflow, ocean, 12 m terrain); xhi -10.55,
yhi -15.03 m/s outward (outflow, land, 264/279 m).

| d | Davies xlo | NSCBC xlo | Davies ylo | NSCBC ylo |
|---|---|---|---|---|
| 0 | **141.09** | 2.35 | 33.99 | 0.07 |
| 1 | **175.66** | 5.13 | 60.11 | 0.35 |
| 5 | 82.26 | 8.96 | 22.61 | 2.18 |
| 25 | 9.78 | 25.12 | 9.68 | 20.05 |

Davies dumps **141-176 mm of spurious 24-h rain in the first two cells of the flat
ocean inflow wall**, decaying inward over ~10 cells -- the band w artifact expressed
in precipitation. NSCBC puts 2.35 mm there. That artifact is eliminated.

**Spin-up and contamination separate cleanly.** Under NSCBC the deficit appears only
on the two INFLOW walls (xlo 2.35, ylo 0.07 mm at d=0, recovering over ~10-12 and
~14-16 cells = 30-48 km); the excess appears only on the two OUTFLOW walls, which
are also the two LAND walls, and stratifying by terrain at d>=20 shows it is
orographic: flat cells (<100 m) mean 22.11 / max 101.1 mm, 100-400 m cells mean
72.34 / max 244.9 mm. The interior maximum sits at d=31 on 283 m terrain -- deep
interior, not boundary. Over flat water at d>=20 the max is 101 mm, *below* the
~130 mm observed. So the 30-48 km recovery is Roberge spin-up from the 10 zeroed
hydrometeors (an order of magnitude short of their 300 km worst case), and the
boundary scheme's own contaminated region is d=0-2.

**The interior is wetter under NSCBC (22.77 vs 13.61 mm) while the total is lower
(30.36 vs 35.30 mm).** Hypothesis, not established: the Davies band's spurious
ascent wrings moisture out at the inflow edge before it enters, so removing the
artifact lets that moisture reach the interior and rain there. Davies' near-exact
d>=20 match to the grib-footprint ERA5 mean (13.61 vs 13.46) was a coincidence of
denominators, and against the matched footprint it is 4.29x, not 1.01x.

**The d=11-16 shoulder does not appear in precipitation**, and on the wet day it is
not distinct in w either (d=13 3.5%, *below* the 5.37% interior background; d=16
6.7%, marginally above).

### Minor ERF bug: a clipped final step corrupts rain_accum

When the last step is clipped onto `stop_datetime` (here dt = 0.0117 s), ERF writes
an extra plotfile in which `rain_accum` is NaN over 5852 of 8192 columns while every
prognostic field is clean and the run exits 0. The preceding output at the same
wall-clock time is correct. Any analysis that takes the last plotfile silently gets
a NaN precipitation field.

### Moisture budget: redistribution is measured, and the interior is still spinning up

**Budget across the d = 20 contour** (88x24 cells, 19.0 x 10^9 m2), total water
(qv + qc + qr + qgraup + qsnow), 2-hourly sampling over 24 h:

| term | NSCBC | Davies |
|---|---|---|
| inward flux across contour, IN | -1.29 mm | -14.86 mm |
| storage change, dS | +6.42 mm | +4.60 mm |
| precipitation, P | 22.77 mm | 13.61 mm |
| residual IN - dS - P | -30.48 mm | -33.08 mm |

The residual should be surface evaporation and should be POSITIVE, so the budget
does **not** close in absolute terms -- the flux uses plotfile (cell-centred)
velocities on the contour, a differenced-height dz, and 2-hourly sampling. Only the
DIFFERENCE between the runs is trustworthy, since both use an identical grid, the
same code path and near-identical evaporation.

    d(IN)         = 13.57 mm
    d(P) + d(dS)  = 10.98 mm      [ d(P) = 9.16, d(dS) = 1.82 ]

**These match to 19%**, with the 2.6 mm gap being the difference in the (undiagnosed)
evaporation terms. So the redistribution hypothesis is confirmed as far as this can
confirm it: NSCBC delivers ~13.6 mm more total water across d = 20 than Davies, and
~11 mm of that appears as extra interior precipitation plus storage. The extra
interior rain is not manufactured -- it is moisture Davies was destroying at the
inflow edge.

**Spin-up: there is no overshoot, and recovery is far from complete.** Binning by
distance from the two INFLOW walls (xlo, ylo), flat cells (<100 m) only, d >= 20:

| d_inflow | n | NSCBC | Davies | excess |
|---|---|---|---|---|
| 20-30 | 1016 | 20.68 | 10.16 | **+10.52** |
| 30-40 | 796 | 22.47 | 14.63 | **+7.84** |
| 40-55 | 272 | 26.37 | 19.88 | **+6.49** |
| 20-30 (terrain >=100 m) | 4 | -- | -- | +31.87 |
| 30-40 (terrain >=100 m) | 24 | -- | -- | +21.64 |

Two readings, and the second matters more.

1. The excess falls monotonically with distance from the inflow walls (10.52 ->
   6.49), so it is inflow-related rather than spatially uniform. Orographic cells
   amplify it 2-3x, though on a small sample (4 and 24 cells).
2. But absolute precipitation **rises** monotonically with distance from the inflow
   walls in BOTH runs, out to d_inflow = 40-55 (120-165 km). Nothing overshoots and
   comes back down. **This is not spin-up overshoot -- it is a spin-up DEFICIT that
   has not finished recovering at 165 km.** That is much closer to Roberge et al.'s
   300 km than the 30-48 km estimated earlier from the ring-distance profile, which
   mixed inflow and outflow walls together.

**Consequence, and it inverts the expected argument for supplying hydrometeors.**
If the interior is still under-precipitating from the zeroed hydrometeors, then
supplying them would push interior precipitation UP -- making the interior wet bias
(currently ~6-8x the matched-footprint ERA5) **worse**, not better. Supplying
hydrometeors fixes the near-edge deficit and the usable-interior question; it should
not be expected to fix the wet bias, and may aggravate it.

### Supplying hydrometeors is far cheaper than estimated: the frames already carry qc and qr

`FillForecastStateMultiFabs` already reads and interpolates cloud and rain water --
the frame variable list is `{"rho","uvel","vvel","wvel","theta","qv","qc","qr"}` and
`tmp_qc` / `tmp_qr` are computed on the ERF grid -- and then **discards them**: only
`fine_cons_arr(RhoQ1_comp) = tmp_qv` is stored
(`ERF_WeatherDataInterpolation.cpp:239`). So qc and qr require no new data, no WPS
chain and no re-download:

1. store `tmp_qc` / `tmp_qr` into their Morrison component slots;
2. add QC/QR to the boundary planes (enum, sizing, fill) alongside `HindcastBdyVars::RHO`;
3. set `cons_read` / `ind_map` for those components in `fill_from_realbdy`;
4. drop them from the blanket `zero_here = (comp_idx > RhoQ1_comp)`.

Roughly 30-40 lines across three files.

What is NOT cheap is **ice**, which is the species Roberge found dominates the winter
benefit. The `.bin` frame format carries only qv/qc/qr, so ciwc/cswc would require
regenerating every frame from ERA5 (which does have them) plus the number
concentrations Morrison wants. That is a data-regeneration job, not a plumbing job.

### Filed upstream

The clipped-final-step `rain_accum` NaN is filed as erf-model/ERF#3491.

### Observational verification: the interior wet bias is real

**The Stage IV file does not appear to be the 24-h field.** Three independent
signs: (a) its CF `time` variable says 2023-01-10 00:00 UTC but its `data_time`
global attribute says 12:00; (b) against MRMS on 0.25 deg SoCal boxes it correlates
at **r = 0.879** -- same storm, same spatial pattern -- but runs **2.87x larger in
the median**; (c) its CONUS max is 651 mm (25.63 in) against MRMS's 240.6 mm. Good
spatial correlation with a ~3x magnitude offset is the signature of a LONGER
accumulation window over the same event, and the file carries PRISM `normal` /
`departure_from_normal` / `percent_of_normal`, i.e. it is the AHPS observed-precip
product which is served for 1/3/7/14/30-day windows. January 2023 was a sequence of
atmospheric rivers, so a multi-day total would be ~3x the daily one.

MRMS's metadata is self-consistent (validityDate 20230110, validityTime 0, 24 h
ending 2023-01-10 00Z = exactly the run window), so **MRMS is used as the reference
and Stage IV is set aside** -- not by preference but because its window cannot be
confirmed and the evidence says it differs.

Common grid = the ERA5 0.25 deg boxes (coarsest of the four, and the only grid on
which ERA5 needs no interpolation). Aggregation is area-weighted/conservative, not
bilinear. Land mask = ERF terrain > 30 m in >= 50% of the box (ocean is exactly
12 m in the terrain file); radar QPE is not scored over water, and MRMS Pass 2 is
gap-filled so its ocean values are not observations.

Land boxes, n = 20:

| field | mean | p50 | p90 | p99 | max |
|---|---|---|---|---|---|
| MRMS | 9.1 | 1.8 | 26.3 | 48.1 | 50.8 |
| ERA5 | 12.9 | 7.5 | 35.2 | 46.0 | 47.5 |
| ERF NSCBC | 108.5 | 72.3 | 265.9 | 291.1 | 295.9 |
| ERF Davies | 56.0 | 50.4 | 90.2 | 99.5 | 100.5 |

    vs MRMS:  ERA5 1.35x   NSCBC 8.06x   Davies 3.65x

**ERA5 is not 4x low against observations -- it is 1.35x HIGH.** So the interior
ratio against ERA5 is not an artefact of ERA5 under-resolution: the model has a
genuine wet bias, and NSCBC's is larger than Davies'.

Caveats, all material: only 20 land boxes survive the mask, and they sit on the
mainland fringe at the north/east edges -- spatially biased toward both orography
and the outflow boundary where NSCBC carries its orographic signal. MRMS p50 over
them is only 1.8 mm, and coastal SoCal radar coverage is imperfect.

**Terrain stratification of the ERA5 ratio (ERF grid, d >= 20):**

| terrain | n | NSCBC | Davies | ERA5 | NSCBC/ERA5 |
|---|---|---|---|---|---|
| flat < 100 m | 2084 | 22.11 | 13.14 | 8.50 | **2.60x** |
| 100-400 m | 28 | 72.34 | 49.24 | 20.59 | **3.51x** |

The excess is **not terrain-only**: flat-water precipitation is 2.60x ERA5, which
by the stated criterion is a genuine bias rather than a resolution artefact.
Orography amplifies it to 3.51x but does not create it. (These are lower than the
earlier 7.18x because that used the d >= 20 ERA5 subset of only 22 boxes, which is
drier than the ERA5 cells actually overlying the flat interior.)

### Root cause of the wet bias: hardcoded tropical surface fluxes

`ApplySurfaceTreatment_BulkCoeff_CC` (ERF_ApplySurfaceTreatment_BulkCoeff.cpp:88-93),
active whenever `init_type = HindCast && hindcast_surface_bcs` (MakeSources.cpp:438):

    Real Ch = 0.0015*(1-ls_mask);  Real Ce = 0.0015*(1-ls_mask);
    Real dT = max(0, 301.0 - temp);        // 301 K = 28 C
    Real dq = max(0, 0.024  - qv);         // 0.024 kg/kg
    cell_rhs(RhoQ1_comp) += rho*Ce*velmag*dq/dz;

`surface_state_arr` is consulted ONLY as a land/sea mask. The ERA5 SST is never
used in the flux -- both driving values are hardcoded tropical constants.

Measured at hour 6 over ocean cells (k=0): T = 294.7 K (21.6 C), qv = 0.01698
kg/kg, |U| = 9.1 m/s, rho = 1.12. q_sat at that temperature is **0.01621 kg/kg**,
so the near-surface air is ALREADY SUPERSATURATED with respect to its own
temperature and physical evaporation should be zero or negative.

| target | dq | E |
|---|---|---|
| hardcoded 0.024 | +0.00702 | **9.2 mm/day** |
| q_sat at model T | -0.00077 | 0 (clamped) |

So the term injects ~9 mm/day of moisture into already-saturated air. The
companion heat term drives temperature toward 301 K against an actual 294.7 K,
adding 6.3 K of spurious sensible heating, which raises q_sat and lets the column
hold the spurious moisture until it rains out -- a closed positive feedback.

Of the four candidate terms: **q_sat carries it, and it is a hardcoded constant
rather than a closure**. C_q = 1.5e-3 against a typical ocean 1.2e-3 is only 25%
high; |U| = 9.1 m/s is reasonable; and q_air is not too dry (the MYNN
under-mixing hypothesis has the wrong sign -- the surface air is too MOIST).

This is independent of all the boundary-condition work and affects every HindCast
run with `hindcast_surface_bcs = true`.

### Jan-9 is the wrong validation case for this domain

Domain corners (lon, lat): SW (-121.326, 32.468), SE (-117.240, 32.469),
NW (-121.368, 34.199), NE (-117.200, 34.200). **North edge is 34.20 N.**

**The Santa Ynez Mountains (34.50 N) are OUTSIDE the domain**, 10 cells (~30 km)
north of the north boundary. MRMS there is 165.9 mm -- the storm's maximum, and
the source of the "~130 mm observed" reference used earlier, fell outside the
domain entirely. That comparison was never valid here.

Point verification against MRMS (24-h, mm):

| point | MRMS | NSCBC | Davies |
|---|---|---|---|
| Santa Catalina Is | **0.0** | **115.7** | 123.2 |
| San Clemente Is | **0.0** | **283.6** | 69.6 |
| San Nicolas Is | **0.0** | 51.3 | 12.5 |
| Santa Cruz Is | 63.9 | 102.5 | 252.3 |
| LA / Santa Monica | 15.1 | 66.9 | 125.9 |
| Long Beach | 0.3 | 31.9 | 56.8 |
| Oxnard coast | 34.4 | 44.0 | 126.4 |

**Direct falsification, no regridding or ratio required: the model produces
116-284 mm of 24-h rain at three island sites where essentially nothing fell.**
(Offshore MRMS coverage is imperfect, but Catalina is within LA-area radar range
and is independently reported as near-zero.)

So Jan-9 samples the EDGE of the event, not its core: the maximum is north of the
domain, and in-domain observations run 0-64 mm. That makes it an excellent
falsification case and a poor skill case. The earlier "20 land boxes" were
sampling the LA/Santa Ana fringe, which is why MRMS maxed at 50.8 mm there.

### Two surface flux closures run simultaneously, and NEITHER uses the ERA5 SST

`ApplySurfaceTreatment_BulkCoeff_CC` is gated on
`init_type == InitType::HindCast && solverChoice.hindcast_surface_bcs`
(ERF_MakeSources.cpp:438) with **no check on `zlo.type`**. So any HindCast deck that
also selects a surface-layer lower boundary gets TWO surface flux closures at once:

* MOST / moeng, as the `zlo.type = "surface_layer"` boundary condition, with
  stability-dependent coefficients; and
* the bulk source term, added to `cell_rhs` at k = 0 with fixed Ch = Ce = 0.0015.

They are structurally additive -- one is a boundary flux, the other an interior
source on the same cell -- so the surface flux is applied twice. Removing
`erf.most.Cd/Ch/Cq` from the deck (done earlier in this campaign to avoid the
bulk_coeff hard abort) disabled MOST's *bulk_coeff flux option*; it did nothing to
this source term, which has been running the whole time.

**And neither closure uses the ERA5 SST.**

* The bulk source consults `surface_state` only for the land/sea mask (comp 0) and
  drove toward hardcoded 301 K / 0.024 kg/kg. The SST is sitting unread in comp 1.
* MOST's `t_surf` is initialised to `default_land_surf_temp`
  (ERF_SurfaceLayer.H:445) = `erf.most.surf_temp`, and `set_t_surf` is called ONLY
  for `InitType::Input_Sounding` (ERF.cpp:1615). Over ocean the LSM never updates
  it, so `t_surf` stays at the deck's **constant 288 K everywhere**, and
  `fill_qsurf_with_qsat` builds q_surf from that constant.

The deck comment `erf.most.surf_temp = 288.0  # fallback only; MM5/SST provide the
live surface state` is therefore wrong -- nothing overrides it over water.

MOST is NOT inert: it computes q_sat over sea and applies stability-dependent
transfer. Its 288 K (q_sat ~ 0.0106 kg/kg) is a physically plausible January SoCal
value, unlike the bulk path's 301 K / 0.024. So the right move is to disable the
bulk source and keep MOST -- but MOST still needs `t_surf` wired to the ERA5 SST to
be spatially and temporally correct, which is a separate change from removing the
double count.

### This bug was found once before and misdiagnosed

The earlier surface-BC instability was correctly traced to "fixed-coefficient ocean
moistening" but was attributed to `qv = 0` in the inert initial condition rather than
to the moistening TARGET being tropical. Fixing the IC removed the visible symptom
and left the constant live for weeks, during which it was silently adding ~9 mm/day
of moisture to already-saturated marine air.

### Correction for the record

The "~130 mm observed interior maximum" used as a reference throughout this campaign
is Santa Ynez at 34.50 N. The domain's north edge is 34.20 N, so that station is
~30 km OUTSIDE the domain and its value was never comparable. Every ratio derived
from it -- including reported passes -- is void.

### THE LARGER FINDING: no HindCast run has ever used the ERA5 SST

Both surface flux closures ran on constants for the entire campaign:

* the bulk source read `surface_state` comp 0 (the land-sea mask) and ignored comp 1
  (the SST), driving toward hardcoded 301 K / 0.024 kg/kg;
* MOST's `t_surf` was initialised to `default_land_surf_temp` = `erf.most.surf_temp`
  (ERF_SurfaceLayer.H:445) and `set_t_surf` is reachable only from
  `InitType::Input_Sounding` (ERF.cpp:1615), so HindCast never updated it. Over ocean
  the LSM never touches it either, so it stayed at the deck's 288 K everywhere, and
  `fill_qsurf_with_qsat` built q_surf from that constant.

So **every HindCast run to date has had a spatially and temporally uniform sea
surface**, on a domain that is ~75% ocean with a real SST gradient across it. That is
a substantial deficiency independent of the tropical constants, and it means the deck
comment `erf.most.surf_temp = 288.0  # fallback only; MM5/SST provide the live surface
state` was wrong for the whole campaign -- nothing ever provided the live state.

Fixed by pushing the ERA5 SST into `t_surf` over ocean points each time
`surface_state_interp` refreshes (ERF_SurfaceDataInterpolation.cpp). Land points are
left to the LSM (here the MM5 `soil_theta` constant, since MM5 is inert), selected on
the land-sea mask; the SST guard is two-sided because the field carries fill values.

### Measured double-count factor

`erf.hindcast_bulk_surface_flux` (new, default 0 = off) gates the bulk source so the
two closures can be scored separately. Jan-9, first 1.11 h, ocean cells at k = 0,
**with the SST fix already applied to both paths**:

| config | d(qv) kg/kg | d(T) K | precip mm |
|---|---|---|---|
| MOST only | 0.003234 | -7.837 | 0.023 |
| MOST + bulk | 0.003497 | -7.695 | 0.040 |

Bulk-path contribution: d(qv) +0.000263, d(T) +0.142 K. As a factor:
**qv 1.08x, T 0.98x, precipitation 1.71x.**

The moisture tendency is only 8% larger but precipitation is 71% larger -- rain is a
threshold response, so a small extra moistening of near-saturated marine air produces
a disproportionate precipitation increase. That is with the SST fix in place; with the
original 301 K / 0.024 constants the bulk contribution was far larger.

The d(T) = -7.8 K cooling in both configurations is the SST fix working: the initial
294.7 K near-surface air is now relaxing toward the real ~288 K sea surface instead of
being driven toward 301 K.

Default is now MOST only. The `hindcast_surface_bcs` DATA path is untouched -- it is
the flux closure being removed, not the SST/land-mask supply.

### CORRECTION: the SST guard DOES have false accepts, on the grid that matters

Instrumenting the fill to score filled cells against originals produced:

    SST filled   cells: n=1204  min 13.62  mean 25.09  max 62.41 C
    SST original cells: n=5832  min 13.42  mean 14.50  max 62.41 C

against a true clean frame range of 13.1-15.7 C, mean 14.3 C.

**A 62.41 C (335.5 K) value is passing the 200-340 K guard, and it is in the
ORIGINAL (guard-passing) set, not the filled set.** The earlier claim of "zero false
accepts" was measured on the 33x20 FRAME grid; the guard actually operates on the
128x64 ERF grid, after `FillSurfaceStateMultiFabs` bilinearly interpolates the frame.
On the coarse frame grid the smear jumped straight from ~287 K to >340 K with nothing
in the plausible band. Interpolating onto the finer ERF mesh fills that gap
continuously, so intermediate values now land at every temperature between clean SST
and the 9999 fill -- including inside the guard.

So the contamination is worse on the ERF grid than in the source data, and the fill
is not the culprit: it faithfully propagated already-bad neighbours, which is why the
filled mean (25.09 C) sits so far above the original mean. The fill IS mildly
flattening as well (it is distance-ordered, so cells filled on later sweeps average
means-of-means), but that is second-order next to the guard admitting 60 C water.

Tightening the guard to a physical SST range (271-305 K) helps but cannot be
sufficient: an interpolation weight of ~0.002 against a 9999 fill still lands inside
any plausible band.

**The correct fix is ERF-side and cheap** -- and unlike the erftools generator fix, it
is in this tree. `FillSurfaceStateMultiFabs` calls `bilinear_interpolation_2d` on the
raw SST array. It should instead do a MASKED interpolation: accumulate only source
points that are water and in range, renormalising the weights, and flag cells whose
stencil contained no valid source. Then no contaminated value is ever created, the
guard becomes a backstop rather than the primary defence, and q_star/t_star would be
fixed by the same change if they are ever wired.

24-h runs were NOT launched on this binary -- with 60 C water reaching the surface
closure over ~20% of the ocean they would not have been interpretable.

### Masked SST interpolation: fabricated values eliminated

`FillSurfaceStateMultiFabs` now interpolates SST with a masked bilinear stencil --
accumulating only source points that are BOTH water (ls_mask < 0.5) AND in physical
range (271-305 K), renormalising the weights. Both criteria are needed: a land-masked
point can still carry a fill, and a water point can still be smeared in the source.
Stencils with no valid source are FLAGGED (-1) rather than given a fabricated value,
and the nearest-valid fill handles them. The downstream guard is unified to 271-305 K
and is now a backstop that should never fire.

Source audit (logged once at startup, so generator defects are separable from ours):

    HindCast SST source frame: 423 water points, 44 of them outside 271-305 K
    (generator smear); 237 land points carrying fills

Result:

| | before (raw bilinear) | after (masked) |
|---|---|---|
| filled cells | n=1204, mean **25.09** C, max **62.41** C | n=478, mean **14.74** C, max **15.27** C |
| original cells | n=5832, mean 14.50 C, max **62.41** C | n=6565, mean 14.41 C, max **15.38** C |
| cells needing fill | 1211 (~20% of ocean) | 478 (~6.8%) |
| cells left unfilled | 7 | **0** |
| anything outside 271-305 K | yes | **none** |

Whole-field max is now 15.38 C = 288.5 K, against a true clean frame range of
13.1-15.7 C. No value anywhere in the field is outside the physical band, so the
guard never fires.

Note the fill count DROPPED below the source contamination rate (6.8% of ERF water
cells vs 10.4% of frame water points): masked interpolation recovers any cell whose
stencil retains at least one valid neighbour, so it both stops fabricating values and
shrinks the genuine gaps. This is a two-stage design -- masked interpolation creates
no bad values, nearest-valid fill covers real gaps, guard is a backstop.

### Tangential-only band relaxation: FAILS, and explains why no local fix exists

Dropping the wall-normal momentum from the band relaxation (keeping tangential/w/
theta/qv, normal momentum still specified at the wall face) does eliminate
grad(F).(A-B) at straight faces by construction. It also destroys the run.

Inflow-wall 24-h precipitation (mm), corner cells excluded:

| d | TANGENTIAL xlo/ylo | Davies xlo/ylo | NSCBC xlo/ylo |
|---|---|---|---|
| 0 | **3137.7 / 864.4** | 153.6 / 33.6 | 1.5 / 0.1 |
| 1 | **3859.1 / 1314.0** | 183.0 / 56.4 | 3.8 / 0.2 |
| 2 | **3199.9 / 1175.2** | 143.0 / 49.7 | 5.1 / 0.5 |
| 8 | 684.2 / 227.1 | 8.8 / 3.1 | 7.4 / 1.8 |
| 20 | 135.0 / 18.2 | 0.5 / 0.7 | 10.8 / 10.4 |

Skill vs MRMS, land d>=3 (obs 8.36 mm): bias **+546.88**, RMSE **974.82**,
ratio **66.40x** (Davies 4.05x). Corners 1687 mm against Davies' 54.7. Mass drift
-3.222 %/day whole / -2.237 settled, against Davies +2.315 / +0.697.

The pass condition asked for the 140/172/138 mm inflow-wall signal to COLLAPSE. It
went up 20x. Answering the second watch item -- "is the domain still adequately
forced?" -- emphatically no.

**Why: the artifact and the forcing are the same term.** Normal momentum relaxation
is not merely the thing that generates grad(F).(A-B); over a 10-cell band it is what
holds the inflow mass flux consistent. Specifying it at the single wall face leaves
it unconstrained across the remaining nine cells, so with 10-12 m/s inflow the normal
velocity drifts freely, converges, and rains out.

This is the SAME result the Helmholtz kill condition gave from the other direction:
||C_proj||/||C_raw|| = 0.376, i.e. 86% of the relaxation forcing is the curl-free
(normal/divergent) part. Removing the divergent part removes the nudging -- measured
there as a norm, measured here as a 66x wet bias. Two independent methods, same
conclusion: **grad(F).(A-B) cannot be separated from the forcing by any local
operation on the relaxed variable set.**

Note the fallbacks do NOT share this defect: an exponential ramp profile (Marbaix)
and diffusive relaxation (TRAM/Tatsumi) both RETAIN normal-momentum relaxation and
only change the ramp shape or the operator. They remain untested and cheap.

### ROOT CAUSE of the interior wet bias: the IC theta was the base state's theta

`ERF_WeatherDataInterpolation.cpp:597` read

    cons_arr(i,j,k,RhoTheta_comp) = r_arr(i,j,k) * th_arr(i,j,k);   // th_arr = th_hse

so theta for the initial state came from the hydrostatic base state that
`erf_enforce_hse` derives while integrating dp/dz from the interpolated rho -- NOT
from the interpolated ERA5 frame. The frame's theta was read, horizontally
interpolated, vertically interpolated onto the ERF levels, and then discarded; only
rho and qv survived. `f_arr` (the frame) was already in scope and used for qv on the
very next line. Fixed to `f_arr(i,j,k,RhoTheta_comp)`.

**Only one instance exists** -- the boundary planes take theta from `fcons`
(line ~523), so the lateral forcing was never affected. The interior was initialised
~4 K warmer than the theta the boundaries were relaxing it toward.

Measured before/after at t=0 over ocean:

| check | before | after | target |
|---|---|---|---|
| theta monotonic with height | **NO** (falls 3.35 K over lowest 150 m) | **YES** | stable |
| T at 12 m | 24.61 C (**+10.14 K** vs SST) | 14.38 C (**-0.73 K** vs SST+0.64) | ~15.1 C |
| worst \|theta - ERA5\| | **+5.30 K** | **1.41 K** | -- |
| \|p - p_hse\| | 17 hPa | 4.86 hPa | 0 |
| RhoTheta/Rho vs frame theta | +3.77 K | **-0.15 to -0.21 K** | ~0 |

The `RhoTheta/Rho` column now equals the stored theta exactly and sits within 0.21 K
of the frame value interpolated to the same height (residual is
interpolation-scheme difference, not substitution). The domain-wide absolutely
unstable surface layer is gone.

**The 1.41 K residual is the erftools offset, as predicted**: ~300 m x ERA5's
~5 K/km theta gradient = ~1.5 K.

### erftools ~300 m vertical offset (unfixable here, source not in this tree)

Signature, for whoever has the source:

* **Constant to +-7 m over 1100 m**: offsets of +309.4 / +302.6 / +297.1 / +295.6 /
  +301.6 / +304.9 m at frame levels 1-6, derived by taking the frame's (T, p) and
  locating that pair in the ERA5 column.
* **The z array itself is correct**: frame z = 155.1 / 369.2 / 588.2 / 812.2 against
  ERA5 geopotential heights 161.6 / 375.5 / 594.3 / 818.1, i.e. ~6 m low and
  shrinking -- consistent with a small surface-geopotential term, not the bug. It is
  the DATA assigned to each z that comes from ~300 m higher, not the z values.
* **~1.3-1.4 level spacings** (levels are 214-235 m apart), so NOT a clean index slip.
* Constant rather than growing rules out hypsometric integration from a wrong surface
  pressure; it points at a wrong z<->p relation (e.g. a standard atmosphere) used when
  placing pressure-level data onto the z grid.
* Frame theta is +1.51 to +1.56 K uniform vs ERA5 while frame T is -1.4 K, at the same
  level, with p 35 hPa low -- self-consistent (theta = T(p0/p)^kappa checks out at
  964.8 hPa), i.e. the right values for the wrong height, not corruption.

**This residual survives the ERF fix and is worth ~1.5 K.** Any scoring after the ERF
fix should expect ~1.4 K of warm bias to remain, not zero.

### Band relaxation: closed

Third and final band formulation tested. Band-only rho relaxation with the wall left
free (never previously separable -- `hindcast_mass_consistent_bdy` was OR'd into
`l_bdy_rho`, so attempt 2 was always a compound experiment; now decoupled) degrades
the interior 4x: 17.11x bias at land d>=3 against Davies' 4.05x, RMSE 239.66 vs
51.20, and 8.2/2.8 mm at d=20 against Davies' 0.5/0.7. Correlation rose (0.443 vs
0.351) but that does not survive a 4x bias increase. Same domain-wide failure mode as
tangential-only. **Davies is the base; the band is closed.**

## 19. HindCast initialization solved for the wrong variable, and anchored the
## pressure field on a fixed 101325 Pa sea-level value in every column

`init_thermo_from_hindcast` ran the thermodynamic initialization backwards. It took
the interpolated ERA5 **density** as the base state, called `erf_enforce_hse` to
integrate `dp/dz` from it, and let that integration **derive** theta. The frame's
theta was read, horizontally interpolated, vertically interpolated onto the ERF
levels, and then discarded (`ERF_WeatherDataInterpolation.cpp:597` used `th_arr` =
`th_hse`, not `f_arr`). Whatever error the pressure integration accumulated became
the state's theta error.

Three defects compounded:

1. **Solving for the wrong variable.** Measured Jan-9 2023: model theta 293.07 K at
   155 m against the frame's 289.30 K, +3.77 K, growing with height. Worse than the
   magnitude, the SIGN of the near-surface gradient was wrong -- theta *fell* 3.35 K
   through the lowest 150 m over the entire ocean, i.e. a domain-wide absolutely
   unstable marine layer at t = 0 where ERA5 has a stable one.

2. **`erf_enforce_hse` anchors p = p_0 = 101325 Pa at z = 0 in every column**
   (`ERF_Init1D.cpp:276`). That is a fixed standard-atmosphere sea-level pressure. It
   cannot represent a synoptic pressure field, which is the entire content of a storm
   hindcast. Measured ERA5 sp over this domain on Jan-9: 903.7 to 1022.7 hPa.

3. **The frame's bottom level is fabricated.** erftools writes a surface level at
   z = 0 that is BYTE-IDENTICAL to the level above it (z = 155.07 m): maxdiff exactly
   0.000e+00 in all eight fields (rho, uvel, vvel, wvel, theta, qv, qc, qr). So the
   lowest 155 m carries no data at all. Substituting the frame theta without fixing
   this gives constant rho AND constant theta below 155 m, hence constant p over the
   lowest five ERF levels, hence zero layer thickness -- 1,818,000 non-finite RRTMGP
   optical depths and an abort at step 1. The old path masked this only because
   `erf_enforce_hse` made theta vary hydrostatically, at the cost of being +3.77 K wrong.

### Fix (in fork, gated on `erf.hindcast_sfc_anchor_file`)

Specify the thermodynamic profile and solve for the mass field, the direction
real.exe works in:

* theta from the frame; p integrated hydrostatically upward from ERA5's `sp`;
  rho from the equation of state. Base state and state are then the same field, so
  buoyancy is exactly zero at init.
* The integration uses the SAME trapezoidal discretization as `erf_enforce_hse`
  (`p(k) = p(k-1) - dz*g*(rho(k)+rho(k-1))/2` with `rho = rho(p,theta)`, solved by a
  4-sweep fixed point that contracts ~4e-3 per sweep), so the base state sits in the
  model's own DISCRETE hydrostatic balance, not merely a continuous one. p is then
  strictly decreasing by construction.
* The anchor is ERA5 `sp` carried from ERA5's own orography height to ERF's 3-km
  terrain height. Deriving p from the lowest good frame level instead would inherit
  erftools' ~300 m downward displacement -- ~35 hPa of surface-pressure error.
* The fabricated bottom frame level is detected by exact equality across all eight
  fields and dropped. This is numerically a no-op for the interpolation (the
  interpolator clamps below `zvec[0]` either way); what it changes is that nothing
  downstream can mistake the duplicate for a surface observation.
* The layer below the lowest frame level with real data is a linear-in-z blend from
  the ERA5 2-m air temperature. The blend base is the FIRST CELL CENTRE, not ERA5's
  orography height: ERA5's 0.25 deg orography smears land elevation over coastal
  water, and blending from it drove the weight negative in a band of ocean columns
  where it clamped to zero and produced constant-theta layers -- 847 ocean layers
  with d(theta) = 0, the same fabricated-uniform-layer defect being removed.

### Measured at t = 0 (plt00000, Jan-9 2023 00Z)

| condition | result |
|---|---|
| p monotonically decreasing | PASS -- 0 non-decreasing layers, largest dp -206.1 Pa |
| theta increasing over ocean | PASS -- 0 of 215,264 ocean layers non-increasing; min d(theta) +0.058 K |
| non-zero layer thickness | PASS -- min dz 19.54 m, min dp 206.1 Pa |
| T(12 m) vs SST + 0.64 K | PASS -- mean -0.34 K, p5/p95 -1.30/+0.28 K, max 1.90 K |
| theta vs frame, z >= 155 m | mean 0.019 K, p99 0.42 K; **ocean mean 0.0017 K, max 0.42 K** |

Ocean surface-layer lapse is now **+17.8 K/km** (stable) against the old path's
-22 K/km (absolutely unstable). `|p - p_hse|` max 0.16 Pa, i.e. base and state agree.
RRTMGP runs clean at step 1.

### A separate, pre-existing interpolator defect surfaced by the check

The theta-vs-frame comparison has a tail on steep terrain: mean |dtheta| rises
monotonically with terrain slope, 0.0027 K on flat ground to 0.172 K where the
terrain changes >75 m/cell, max 2.378 K. **This is not the initialization.** A
control run with `erf.hindcast_ic_frame_theta = 1` -- raw interpolated frame theta,
no hydrostatic integration, no blend, no anchor -- reproduces the same worst cell
(i,j,k = 90,15,2), the same model value 291.345 K, and the same slope distribution to
four decimals (0.1718 vs 0.1717).

Cause: `FillForecastStateMultiFabs` samples the frame at
`z = (z_phys_nd(i,j,k) + z_phys_nd(i,j,k+1))/2` -- the height of the (i,j) NODE
column, averaged in k only -- while x and y are the CELL CENTRE. On flat ground the
node and cell-centre heights coincide; on a slope they differ by O(terrain
gradient x dx/2), which against a strong inversion is worth whole kelvins. The
cell-centred height is already available as `z_phys_cc`.

## 20. `erf_enforce_hse` hardcodes p = 101325 Pa at z = 0 in every column

Independent of item 19 and of anything hindcast-specific. `ERF_Init1D.cpp:276`:

```cpp
pres_arr(i,j,klo) = p_0 - hz * rho_arr(i,j,klo) * l_gravity;
```

`p_0` is the reference pressure from `ERF_Constants.H`, 1.0e5 Pa. So the base
state's pressure is anchored to a single fixed sea-level value everywhere, with no
way to supply another. There is no input to override it and no diagnostic that it
happened.

**Why it matters beyond the hindcast.** For an idealized problem with a uniform
background this is the right convention and costs nothing. For any case driven by
real data it silently discards the synoptic pressure field:

* ERA5 `sp` over this 384 x 192 km domain on 2023-01-09 00Z spans **903.7 to
  1022.7 hPa** (min over the 1020 m orography in the NE, max offshore). A single
  101325 Pa anchor is wrong by up to 20 hPa at sea level in this one small domain,
  and the error is a smooth horizontal field -- i.e. exactly a spurious synoptic
  pressure gradient.
* A deepening cyclone is a surface-pressure anomaly. Initializing every hindcast
  from the same sea-level pressure removes the feature being hindcast.
* The error is invisible in the usual checks: the base state is still in perfect
  hydrostatic balance, `p - p_hse` is still ~0, and nothing is non-finite. It only
  shows up against an independent surface-pressure observation.

The anchor also interacts with terrain. `hz = z_cc(i,j,klo)` is the height of the
first cell centre above **sea level**, not above ground, so the routine integrates
from z = 0 up to the first cell centre using the local first-cell density -- a
reduce-to-sea-level extrapolation through terrain that may be 1700 m thick, using
the density of the air at the top of it.

**Fix used in the fork** (item 19): anchor on an observed surface pressure at an
observed height and integrate from there. A general fix upstream would be an
optional 2-D surface-pressure field, defaulting to the current behaviour when
absent.

### 19a. The t2m blend was anchored to the frame level's LABELLED height, not its true one

The fix in item 19 blends from the ERA5 2-m air temperature up to the lowest frame
level carrying data, at z = 155.07 m. But that level's data is not from 155 m. Item 5
records erftools displacing the levels **~300 m downward**, so the values labelled
155 m are really from ~455 m. Blending across the labelled depth therefore compresses
a 455 m temperature difference into 153 m:

    imposed:  2.94 K over 153 m  = 19.2 K/km
    true:     2.94 K over 453 m  =  6.5 K/km       ratio 2.96x too steep

The PBL scheme mixes that spurious inversion out immediately. Measured, ocean mean,
corrected-IC Davies run:

| t (h) | th(k=0) | th(k=4) | d(theta) 0->4 |
|---|---|---|---|
| 0 | 286.50 | 288.94 | +2.44 K |
| 1 | 287.68 | 288.92 | +1.24 K |
| 2 | 287.71 | 288.94 | +1.23 K |
| 4 | 287.83 | 289.06 | +1.23 K |

theta at the first cell centre rises **+1.18 K within the first hour** and the
stratification settles at half what was imposed. Predicted surface warming from
mixing the over-steep profile rather than the true one is ~1.0 K; observed +1.18 K.

**Consequence: precipitation switches off.** At 24 h the run is +0.65 K against ERA5
t2m (Jan-10 00Z, ocean), which raises q_sat enough to drop mean low-level RH from
87% to 80%. Cloud coverage collapses from 17.8% of the domain to 0.62%, and
domain-mean 24-h precipitation falls from 22.4 mm to 0.18 mm -- of which 0.16 mm
falls in hour 1 and the rest of the run produces 0.02 mm. Winds and vertical motion
are unaffected (|v| 10.25 vs 9.80 m/s, w p99 3.90 vs 3.98 m/s over land), and column
water vapour is HIGHER than the broken run (29.74 vs 29.18 mm). The model is not
short of moisture or lift; it is ~1 K too warm to saturate.

This is the same class of error item 19 was written to remove -- accommodating the
erftools displacement instead of correcting it, in a new place. Two ways out:

1. Blend to `z_low + offset` with the ~300 m displacement as an input. Local, cheap,
   still built on a known-wrong level table.
2. Correct `zvec` on read by the measured offset. Fixes the blend AND the ~1.4 K
   theta residual AND the ~35 hPa pressure error at every level, not just near the
   ground. The offset is documented as constant (item 5), which is what makes this
   tractable.

The displacement is not the cosmetic ~1.4 K item 5 described. It breaks the
initialization.

### 19b. RETRACTION: the blend defect is NOT what suppressed precipitation

Item 19a diagnosed the Davies precipitation collapse as a consequence of the
over-steep t2m blend. **That causal claim is wrong.** The control: NSCBC and Davies
were run from the byte-identical initial condition (same anchor file, same blend,
same dropped level, same `min layer pressure thickness = 206.0625 Pa`), differing
only in lateral boundary treatment.

| t (h) | NSCBC rain | NSCBC th(k=0) | NSCBC qc>0 | Davies rain | Davies th(k=0) | Davies qc>0 |
|---|---|---|---|---|---|---|
| 0 | 0.0000 | 286.50 | 0.00% | 0.0000 | 286.50 | 0.00% |
| 1 | 1.8399 | 288.60 | 13.95% | 0.1615 | 287.68 | 2.24% |
| 4 | 2.8077 | 287.79 | 19.07% | 0.1650 | 287.83 | 2.37% |
| 7 | 9.2849 | 288.32 | 18.46% | 0.1778 | 287.97 | 0.71% |

NSCBC carries the identical blend defect and rains ~50x more. It also warms MORE in
the first hour (+2.10 K against Davies' +1.18 K) and stays warmer throughout, which
is the opposite of what the proposed "warming raises q_sat and shuts off
condensation" mechanism requires.

**What survives from 19a:** the blend IS ~3x too steep, the geometry argument is
unchanged, and the inversion demonstrably mixes out within the first hour in both
runs. Correcting `zvec` on read remains justified on its own terms -- the ~1.4 K
theta residual and ~35 hPa pressure error are real at every level. It is simply not
the explanation for the precipitation collapse.

**What is now open:** Davies + corrected IC produces 0.18 mm in 24 h while
Davies + broken IC produced 22.4 mm and NSCBC + corrected IC is on pace for ~30 mm.
The collapse is specific to the COMBINATION of Davies relaxation with the corrected
IC, which points at the interaction between the relaxation target and the interior
state rather than at the interior state alone. Note the boundary planes still take
theta from the CLAMPED frame profile below 155 m, i.e. the un-blended one, so the
band is relaxed toward a different near-surface profile than the interior was
initialized with. That is a candidate, not a conclusion.

### 19c. NSCBC NaN forensics: a localized boundary spike, visible before the crash

Reported earlier that the 7.00 h plotfile "looks like the healthy 5 h state". That
was wrong -- it was read on domain means and max|w| only. The signal is in v, at the
wall:

| t (h) | max abs(v) | location | d | n(abs(v)>30) |
|---|---|---|---|---|
| 5.00 | 34.77 | (52,63,22) | 0 | 14 |
| 6.00 | 32.77 | (86,31,3) | 31 | 242 |
| 7.00 | **295.30** | **(73,63,1)** | **0** | **3397** |

max abs(u) confirms the same site: at 5 h its maximum sits at (127,37,31), the model
top, which is the jet and physical; by 7 h it has moved to **(72,63,1)** -- the same
near-surface yhi wall cell, 79.23 m/s. max abs(w) at 7 h is 13.60 at (73,58,2), the
same i, five cells inboard of the same wall.

So the failure is a localized instability on the **yhi boundary near i = 72-73 at
k = 1-2**, already at 295 m/s four minutes before the NaN -- not a domain-wide
convective blowup. It is visible in the plotfile; the earlier read looked at the
wrong variable.

Consistent with the timestep collapse: dz_min is 18.5 m, so w ~ 23 m/s gives
dt = 0.2 * 18.5 / 23 = 0.16 s, against the observed floor of 0.158 s.

### 19d. NSCBC fails when a lateral face reverses from outflow to inflow

19c localized the NSCBC NaN to (72,63,1), the yhi wall. The trigger is not terrain
and not marginal CFL. It is a **regime reversal of the face**.

Wall-mean v at the yhi boundary, k = 1 (v > 0 = outflow at yhi), and the mean over
the six cells that blow up (i = 68..73):

| run | 0 h | 1 h | 2 h | 3 h | 4 h | 5 h | 6 h | 7 h | outcome |
|---|---|---|---|---|---|---|---|---|---|
| `sst_nsc` (broken IC) | +0.5 | +7.2 | +19.2 | +30.9 | | | | | survived 24 h |
| `ic_nsc` (corrected IC) | +1.0 | +18.9 | +14.6 | +16.8 | +22.1 | **+4.7** | **-3.4** | **-15.5** | NaN at 7h04m |
| `ic_ctl` (Davies, same IC) | +2.6 | +3.9 | +4.2 | +4.4 | +4.8 | +5.1 | +5.5 | +5.7 | survived 24 h |

The NSCBC run that survived never reversed -- v at that wall is strictly positive and
growing to +37 m/s. The NSCBC run that died reversed between 5 h and 6 h and blew up
one hour later. Cell-level detail at i = 68..78, k = 1:

    5 h:   +9.5  +8.6  +6.2  +3.7  +4.9  +23.1 +23.8 +26.6 +21.7 +20.2 +17.3
    6 h:   -9.6 -10.8 -18.3 -16.6 -11.9   -7.1  -4.1  -6.7  -7.9  -4.0  -4.8
    7 h: -174.6 -217.7 -237.4 -209.6 -215.2 -295.3  -8.1  -9.2 -15.5 -19.5 -18.3

The six cells that reverse most strongly at 6 h are exactly the six that blow up at
7 h; i = 74 onward reversed only weakly and stayed bounded. The 6 h reversal
magnitude predicts the 7 h failure location.

**Why this is the expected failure mode.** A subsonic face imposes 4+N conditions on
inflow and exactly 1 on outflow. A face that reverses must switch between those two
counts mid-run. That switch is the hard case in any characteristic treatment, and it
is the one this configuration never exercised before: the earlier NSCBC runs all had
persistently outflowing lateral faces.

**Why the corrected IC exposed it.** Davies relaxes v toward the frame, and the
frame's v at yhi stays positive, so a Davies run *cannot* reverse there. NSCBC leaves
the face free, so interior dynamics can reverse it. The corrected IC changed the
interior enough to produce that reversal; the broken IC's convecting interior did
not. So the IC is the proximate cause of exposure but not of the defect.

Note `erf.nscbc_outflow=0` in these runs, i.e. the outflow branch was already
disabled -- so the face was being treated as inflow-or-nothing while physically
switching between the two.

**Falsified along the way:** terrain steepness. The steepest along-wall terrain step
on yhi is at i = 111->112 (271.8 m per 3-km cell); the failure is at i = 72, forty
cells away, where the step is 108 m and does not rank in the top six. The xlo and ylo
walls are flat ocean throughout (12.5 m), so only xhi and yhi carry terrain at all,
but within yhi the failure site is not distinguished by terrain.

### 19e. The boundary planes were the source of Davies' near-surface warm drift

Discriminator, Davies, 7 h, same binary, `erf.hindcast_blend_bdy_theta` off vs on.
Off = interior initialized with the blended profile while the band is relaxed toward
the raw frame (clamped and vertically uniform below the lowest frame level). On =
the blend is applied to the frame itself in `FillForecastStateMultiFabs`, so the
interior and the relaxation target are one field.

| t (h) | off: rain | off: th(k=0) | off: RH(k=4) | on: rain | on: th(k=0) | on: RH(k=4) |
|---|---|---|---|---|---|---|
| 0 | 0.0000 | 286.50 | 68.4 | 0.0000 | 286.50 | 68.2 |
| 1 | 0.1615 | 287.68 | 70.0 | 0.2523 | 286.79 | 70.3 |
| 4 | 0.1650 | 287.83 | 78.5 | 0.3096 | 286.70 | 79.6 |
| 7 | 0.1777 | **287.97** | 81.5 | 0.3396 | **286.75** | 83.3 |

**Confirmed:** the +1.47 K near-surface warm drift over 7 h is caused entirely by the
relaxation target disagreeing with the interior initialization. With the planes
blended the drift is +0.25 K. Measured mismatch at t = 0 was +1.90 K at the first
cell centre in the band, decaying to zero above ~250 m.

**Not confirmed:** that this explains the precipitation collapse. Rain doubles
(0.178 -> 0.340 mm at 7 h) and RH gains 1.8 points, which leaves the run ~15x too dry
against MRMS instead of ~50x. NSCBC on the same IC had 9.28 mm at 7 h, 27x more.

So the 19a mechanism is real and operative but an order of magnitude too small to
account for a 124x collapse. Something else dominates. The remaining candidate with
a measured magnitude is the erftools level displacement: air taken from ~300 m higher
than labelled is drier as well as warmer in theta, and it enters through both the
boundary forcing and the interior initialization.

Caveat on the A/B: the two arms' ICs are not bit-identical. The frame-side blend uses
`(z_nd(i,j,0)+z_nd(i,j,1))/2`, the node convention that routine already samples the
frame with, while the initializer used cell-centred `z_phys_cc` -- item 19c's
inconsistency again. They coincide on flat ocean (85% of the domain) and differ by
~0.1 K in the worst terrain cell, which cannot produce a 2x rain change.

### 19f. The erftools level displacement was the dominant cause of the dry collapse

`erf.hindcast_frame_z_offset` relabels the frame levels upward by the measured
erftools displacement (item 5), correcting every consumer at once. Single-variable
step from 19e's `bt_on` arm -- same IC construction, same plane blend, offset the only
difference. Davies, 7 h, same binary.

| t (h) | bt_off | bt_on (planes blended) | **bt_zoff (+305 m)** | bt_zoff RH(k=4) | bt_zoff qc>0 |
|---|---|---|---|---|---|
| 0 | 0.0000 | 0.0000 | 0.0000 | 76.0 | 0.00% |
| 1 | 0.1615 | 0.2523 | **3.7267** | 77.5 | 7.15% |
| 3 | 0.1619 | 0.2970 | **7.3774** | 85.9 | 16.88% |
| 7 | 0.1777 | 0.3396 | **10.9182** | 89.5 | 9.53% |

Domain-mean 24-h-style accumulation at 7 h rises **32x** over `bt_on` and 61x over
`bt_off`. Land-only mean at 7 h: 0.1188 -> 0.1253 -> **8.5990 mm**. Mean low-level RH
goes 81.5% -> 83.3% -> **89.5%**, and cloud coverage 0.71% -> 1.20% -> 9.53%.

**Mechanism.** Placing pressure-level data ~305 m below where it belongs means both
the interior initialization and the lateral forcing supply air from ~305 m higher
than labelled. That air is drier as well as warmer in theta. The whole domain was
being held ~8 RH points below saturation by a systematic vertical mislabelling, so
condensation never triggered. Nothing was wrong with the moisture budget, the winds,
or the vertical motion -- all three matched the broken-IC run throughout.

**Item 5's assessment was a large understatement.** It recorded the displacement as
"worth ~1.5 K" and advised expecting "~1.4 K of warm bias to remain" after the ERF
fix. The 1.5 K is real but it is a symptom. The same displacement was suppressing
precipitation by a factor of ~30 through the humidity field, which no theta-only
accounting would surface.

**Geometric confirmation at t = 0.** With the offset the near-surface blend spans the
levels' true depth, and the spurious inversion 19a identified disappears:

    ocean-mean theta, k = 0..4
      blend, no offset:   286.50 286.99 287.56 288.22 288.99   d = +2.48 K
      blend + 305 m:      286.50 286.66 286.84 287.05 287.30   d = +0.79 K

2.48/0.79 = **3.14x**, against 2.96x predicted from the geometry alone.

**Offset value.** 305 m, the midpoint of two independent estimates that agree to 1%:
310 m from the frame's 35 hPa pressure error against dp/dz = -11.3 Pa/m, and 306 m
from its +1.53 K theta error against a ~5 K/km gradient. It is an input, defaulting
to 0, because it is a property of the erftools build that wrote the frames.

**Not yet scored.** This is 7 h; MRMS is a 24-h accumulation, so no bias/correlation
number is available yet. Land-mean 8.60 mm at 7 h against MRMS's 8.51 mm over 24 h
means the run may now overshoot -- accumulation is decelerating (increments 3.73,
2.12, 1.53, 1.42, 0.80, 0.71, 0.61 mm/h) but a 24-h run is needed to say.

## 21. PIPELINE AUDIT: item 5's height-displacement diagnosis is wrong, and two
## compensating errors have been masking each other

`precip_check/pipeline_audit.py` compares every variable across ERA5 -> erftools
frame -> ERF interpolation -> ERF state at t = 0. Its first run overturns three
things at once.

### 21a. erftools does NOT displace data in height. It gets the PRESSURE wrong.

ERA5 -> frame delta, ocean columns, at matched geometric height:

| z (m) | T | qv | u | v | theta | rho | implied p |
|---|---|---|---|---|---|---|---|
| 398 | **-0.03 K** | -5e-5 | +0.03 | -0.02 | **+1.60 K** | **-0.023** | **-1898 Pa** |
| 846 | **-0.02 K** | -6e-5 | +0.04 | -0.03 | **+1.64 K** | **-0.022** | **-1815 Pa** |
| 3130 | **-0.03 K** | -4e-5 | +0.13 | +0.03 | **+1.97 K** | **-0.020** | **-1598 Pa** |

**T, qv, u and v are correct at the labelled height.** A height displacement would
move all of them; T alone rules it out (a 300 m error is ~2 K in T). What is wrong is
the pressure implied by the (rho, theta) pair the frame stores: ~18 hPa low,
uniformly. Everything else follows arithmetically:

    theta = T (p0/p)^kappa    p 1.914% low -> theta +287.5*0.2857*0.01914 = +1.57 K   (obs +1.57)
    rho   = p/(R_d T (1+..))  p 1.914% low -> rho -1.914% = -0.0229                   (obs -0.0224)

Item 5 read the +1.5 K theta and the 35 hPa pressure error as "~1.3-1.4 level
spacings of downward displacement". Both symptoms are real; the cause is not
displacement. The frame's T is the trustworthy field and its theta and rho are not.

### 21b. The item-19 restructure inherits that error and is 1.7-2.2 K too warm

The restructure takes theta from the frame and anchors pressure on ERA5 `sp`. Since
the frame's theta is +1.6 K too high *because* it was built with a pressure 19 hPa
too low, pairing it with a CORRECT pressure converts the pressure error directly into
a temperature error. ERF state minus frame, no-offset run:

| z (m) | dT | dp |
|---|---|---|
| 150 | **+1.58 K** | +2028 Pa |
| 398 | **+1.66 K** | +1948 Pa |
| 846 | **+1.72 K** | +1923 Pa |
| 3130 | **+2.16 K** | +1900 Pa |

Against ERA5 the state is +1.71 K at 846 m and +2.13 K at 3130 m. That is
**the cause of the precipitation collapse**: ~2 K raises q_sat ~13% and drops RH ~8
points, which is exactly the measured 87% -> 80% and exactly enough to stop
condensation. Not the blend (19a, retracted), not the boundary planes (19e, worth 2x).

The OLD path was accidentally self-consistent: it used the frame's rho (also -2%
wrong) with a p_0-at-sea-level anchor (item 20, also wrong), and the two errors
partially cancelled to give roughly the right T. Replacing one of them alone exposed
the other.

### 21c. `erf.hindcast_frame_z_offset` is a compensating error, not a fix

Shifting the levels up 305 m makes the model sample air from 305 m LOWER, which is
warmer and moister, offsetting 21b's warm-dry bias and generating rain. It does not
correct anything. With the offset applied, the fields that WERE right become wrong:

| z (m) | ERA5 T | frame T | with +305 m | ERA5 qv | frame qv | with +305 m |
|---|---|---|---|---|---|---|
| 3130 | 272.04 | 272.01 | **273.86** | 0.00268 | 0.00264 | **0.00341** |

+1.85 K and +29% moisture aloft, both spurious. The 32x precipitation gain in item
19f is real as a measurement and wrong as an attribution: it came from wetting the
column, not from fixing a displacement. **19f's conclusion is retracted; the knob
should not be used.**

### The actual fix

Use the frame's **T**, which is correct, rather than its theta, which is not.
Derive theta from frame T and the hydrostatically-integrated pressure that item 19
already computes from ERA5 `sp`. That leaves every trustworthy field untouched and
drops the erftools pressure error entirely, instead of cancelling it against
another error.

## 22. 24-h Jan-9 Davies on the corrected IC: amplitude fixed, placement worse

First scoring taken against an initial condition verified correct against ERA5
(T 0.00 K, theta -0.05 K, rho +0.05%, p +28 Pa at every level; statically stable
everywhere). Land d>=3, 941 cells, against MRMS.

| metric | baseline (broken IC) | corrected IC | MRMS |
|---|---|---|---|
| bias | 3.99x | **0.55x** | 1.00 |
| correlation | +0.357 | **-0.107** | |
| RMSE | 51.03 mm | **19.39 mm** | |
| mean | 33.90 mm | 4.69 mm | 8.51 mm |
| p50 / p90 / max | 20.1 / 66.1 / 470.7 | 2.1 / 12.4 / 33.1 | 0.2 / 36.1 / 69.1 |
| island mean (obs 0.00) | 5.82 mm | **0.20 mm** | 0.00 |

**Amplitude is fixed.** Bias 7x closer to unity, RMSE down 62%, and the island
points -- where MRMS observes exactly zero -- go from 5.82 mm to 0.20 mm.

**Placement is worse, and worst where it matters.** Correlation degrades as the
subset is restricted toward where it actually rained:

| subset | n | corrected | baseline | corr(MRMS, terrain) |
|---|---|---|---|---|
| all land d>=3 | 941 | -0.107 | +0.357 | -0.135 |
| north band j>=48 | 594 | -0.287 | +0.327 | -0.130 |
| MRMS > 1 mm | 357 | -0.620 | +0.183 | +0.058 |
| MRMS > 5 mm | 221 | **-0.826** | -0.067 | -0.032 |

The hypothesis that the baseline's +0.357 was merely the north-south gradient is
**falsified**: it holds +0.327 within the north band alone. The baseline had modest
real pattern skill; the corrected run has anti-skill where the rain is.

The large-scale gradient is right, though. Land d>=3 means by j band (S->N):

    MRMS      0.00   0.00   0.05  13.45
    corrected 0.10   0.03   2.45   6.28
    baseline  2.95   4.78  22.39  42.79

Correct shape, about half the amount, against a baseline that was wet everywhere.

### Leading hypothesis: hydrometeor spin-up, not the IC

Flow is northward (v at ylo +7.2 m/s inflow, yhi outflow). MRMS puts the event at
j = 48-64, i.e. **144-192 km from the inflow edge** -- squarely inside the measured
120-165 km spin-up recovery distance for a model driven with hydrometeors zeroed at
inflow. The model cannot have mature precipitation there; it is still growing
condensate from scratch. That predicts exactly what is observed: roughly the right
domain total, the right large-scale gradient, and the worst pattern error
concentrated where the observed rain is.

qc and qr are already interpolated onto the ERF mesh and then discarded (item 21,
section 4 of the audit). Wiring them into the cons components and the boundary
planes is the next test, and it is now the leading candidate rather than a
nice-to-have.

The IC work is not implicated in the placement error: the IC is verified correct
against ERA5 to 0.05 K and the amplitude metrics all moved the right way.

## 23. Displacement hypothesis falsified; FSS shows the model does have skill

### 23a. There is no georeferencing error

Island terrain maxima against known positions -- unambiguous point features:

| island | expected (i,j) | terrain there | local max | offset |
|---|---|---|---|---|
| Santa Cruz | (49,56) | 358.3 m | (48,56) | -3, +0 km |
| Catalina | (90,33) | 346.3 m | (91,32) | +3, -3 km |
| San Clemente | (88,15) | 307.6 m | (89,14) | +3, -3 km |
| San Nicolas | (57,27) | 117.8 m | (56,28) | -3, +3 km |

All four within ONE 3-km cell. The LCC parameters, the DEM converter and the
domain anchor are all correct. A 60-150 km georeferencing error is ruled out.

### 23b. The lag-correlation peak was a subset artifact

Extending the scan to +-150 km does not bring the peak interior -- it moves it to
the search edge (lag_x -150 km, corr +0.800). The giveaway is that the whole
lag_x = -150 km ROW is uniformly high (+0.61 to +0.72 across every lag_y). A real
rigid displacement produces a localized peak; a uniformly elevated row is what
happens when the retained overlap shrinks to a narrow strip whose internal
gradient dominates the correlation. Pure-y best is only +0.130 at +45 km, i.e. no
coherent advective displacement either.

**The displacement hypothesis is retracted.** It was reported as a strong candidate
on the basis of a +-60 km scan whose peak sat at the boundary; bracketing it
properly, as instructed, killed it.

### 23c. What survives: FSS

| thresh | scale | corrected | baseline |
|---|---|---|---|
| 1 mm | 3 km | **0.576** | 0.563 |
| 1 mm | 30 km | **0.666** | 0.589 |
| 1 mm | 60 km | **0.700** | 0.611 |
| 5 mm | 15 km | **0.501** | 0.440 |
| 5 mm | 60 km | **0.600** | 0.486 |

Better at every scale and threshold, above 0.5 at all scales for 1 mm. CSI 0.405
vs 0.391. The baseline's POD of exactly 1.000 with FAR 0.609 is what a 4x wet bias
looks like: it detects everything because it rains everywhere.

Pearson's -0.107 is not a measurement of placement skill. It is a metric dominated
by extremes (MRMS p50 is 0.2 mm) applied to a field with a different intensity
distribution.

### 23d. Stage IV comparison is NOT established

Units are confirmed `Inches` from the variable metadata, so the x25.4 conversion is
right. But the file is a single 2-D field with no time dimension and a maximum of
25.63 in = 651 mm, and it reads 3.4x MRMS over the scoring footprint. Whether that
is a genuine inter-product disagreement or a mismatched accumulation window cannot
be determined without the file's provenance. **The earlier claim that "two
observational products disagree by 3.6x" is withdrawn pending that check**, and
with it the +0.554 inter-product correlation, which is equally window-dependent.

## 24. Moist-scalar diffusion: mechanism confirmed, but it does not buy skill

`moistscal_*_adv_type` is ALL-OR-NOTHING. `ERF_SlowRhsPost.cpp:405` sets
`num_comp = n_qstate` at `ivar == RhoQ1_comp`, advecting the whole moisture block --
six Morrison mass mixing ratios AND every number concentration -- in one call with
one scheme. Splitting mass from number would need two calls plus two more AdvType
inputs and their stencil-width validation.

**WENOZ5 fails deterministically.** Both attempts NaN'd at step 100 in component 11
(RhoQ8, a number concentration) at k = 26-27. Single-precision overflow in the
smoothness indicators on number-concentration gradients, exactly as predicted. The
corrected IC gentles the thermodynamic gradients but not these -- they are set by the
microphysics.

**Upwind_5th runs clean for 24 h and confirms the diagnosis:**

| metric | Upwind_3rd (ic_hyd) | Upwind_5th (wz) | MRMS |
|---|---|---|---|
| spectral ratio 8-19 km | 0.496 | **0.661** | 1.000 |
| ratio at 10 km | 0.405 | **0.598** | |
| ratio at 8 km | 0.486 | **0.942** | |
| p90 | 12.95 | **16.62** | 36.07 |
| max | 34.09 | **39.53** | 69.05 |
| bias | 0.566x | **0.727x** | 1.00 |
| CV | 1.29 | 1.28 | 2.04 |

**But every placement metric moves the wrong way:**

| metric | Upwind_3rd | Upwind_5th |
|---|---|---|
| FSS 1 mm @ 3 km | 0.576 | **0.566** |
| FSS 5 mm @ 3 km | 0.423 | **0.395** |
| FSS 5 mm @ 60 km | 0.600 | **0.576** |
| CSI @ 1 mm | 0.405 | **0.394** |
| FAR @ 1 mm | 0.534 | **0.554** |
| correlation | -0.106 | -0.124 |

FSS at 3 km was the stated check on the diagnosis and it did NOT improve.

**Reading.** Implicit diffusion was genuinely suppressing 8-19 km variance -- removing
it restored the spectrum and raised peak intensity and bias exactly as predicted. But
the restored variance lands in the wrong places, so hits and false alarms rise
together (FAR 0.534 -> 0.554) and the categorical scores fall. CV is unchanged at
1.28 against MRMS's 2.04: the field got more intense everywhere rather than more
concentrated.

Diffusion was a real defect and is now half-fixed. It was not what limits skill.
The residual is placement at 3 km, which no advection scheme will supply.

Note the 6-7 km overshoot (ratio 1.13 and 1.48) -- Upwind_5th may be adding
grid-scale noise at the shortest resolved scales, worth watching before adopting it.

## 25. The 3-km run scores BELOW its own 25-km driver at every scale

FSS of ERA5's own 24-h precipitation against MRMS, on the identical footprint,
scales and thresholds used to score the model. ERA5 is the information the run is
forced with, so its FSS is the skill available from the driver.

| thresh | scale | ERA5 (driver) | Upwind_5th | Upwind_3rd | useful thr |
|---|---|---|---|---|---|
| 1 mm | 3 km | **+0.633** | +0.566 | +0.576 | 0.690 |
| 1 mm | 30 km | **+0.679** | +0.649 | +0.666 | 0.690 |
| 1 mm | 60 km | **+0.725** | +0.685 | +0.700 | 0.690 |
| 5 mm | 3 km | **+0.625** | +0.395 | +0.423 | 0.617 |
| 5 mm | 30 km | **+0.722** | +0.522 | +0.549 | 0.617 |
| 5 mm | 60 km | **+0.838** | +0.576 | +0.600 | 0.617 |

Believable scale (smallest neighbourhood with FSS > 0.5 + f/2):

| thresh | ERA5 | Upwind_5th | Upwind_3rd |
|---|---|---|---|
| 1 mm | 60 km | never | 60 km |
| 5 mm | **3 km** | never | never |

ERA5 also has the better amplitude: mean 9.91 mm over the scoring cells against
MRMS's 8.51 (bias 1.16x) versus the model's 6.19 mm (0.73x).

**The downscaling is subtracting skill, not adding it.** A 0.25 deg driver has
useful FSS at 3 km for the 5 mm threshold; the 3-km run built from it does not, at
any scale. This is the "above us" branch: the residual is NOT the configuration's
floor, and something in the model is actively degrading information the driver
already contains.

Consequences for the open decisions:

* Choosing between Upwind_5th and Upwind_3rd is premature -- both sit below the
  driver at every scale, and the gap (0.20-0.26 FSS at 5 mm) dwarfs the difference
  between them (0.03).
* CONUS404 cannot be justified on skill grounds while the run scores below its
  present driver. Better forcing does not help a pipeline that degrades what it is
  given. That decision should wait.
* The next question is not which physics option to tune but WHERE the driver's
  information is being lost -- boundary relaxation, the 10-cell band, vertical
  interpolation onto 32 stretched levels, or the microphysics.

### 24a. advect_tke: unpinned by the IC fix, but worth nothing to skill

`erf.advect_tke=true` was set false during stability work on the broken IC. On the
corrected IC with Upwind_5th it runs the full 24 h clean. Scored against the
otherwise identical `wz` arm:

| metric | tke=T | tke=F |
|---|---|---|
| spectral ratio 8-19 km | 0.658 | 0.661 |
| bias | 0.723x | 0.727x |
| p90 / max | 16.59 / 39.28 | 16.62 / 39.53 |
| CV | 1.25 | 1.28 |
| FSS 1 mm @ 3 km | +0.571 | +0.566 |
| FSS 5 mm @ 3 km | +0.399 | +0.395 |
| FSS 5 mm @ 60 km | +0.580 | +0.576 |

Every metric moves under 1%. No skill basis to prefer either; advecting a
prognostic TKE field is the physically consistent choice, so `true` is defensible,
but it is a configuration decision and not a result.

Same pattern as the moist-scalar fallback (item 24): a knob pinned on a broken IC,
freed once the IC was correct, and worth ~nothing to skill. Both were real defects
in provenance; neither was what limits the model.

## 26. `rad_freq_in_time` is not honoured in single precision: absolute float32 clock where a difference was intended

Measured 2026-07-26 on the 192x96 benchmark. The deck asks for radiation every
180 s of model time. **The effective cadence is 256 s**, and every scored run in
this campaign has had it.

`ERF_Radiation.cpp:210`:

```cpp
m_update_rad = (m_step == 0) ||
               ((m_time - m_last_rad_time) >= m_rad_freq_in_time);
```

`m_time` is ABSOLUTE epoch seconds -- `gmtime(time_t(time))` a few lines above
resolves it to 2023-01-09, so it is ~1.673e9. In a single-precision build
`amrex::Real` is float32, whose ULP at 1.673e9 is **128 s**. Model times are
therefore stored quantised to multiples of 128 s, the difference of two of them
is a multiple of 128, and the first value that satisfies `>= 180.0` is 256.

Any request in (128, 256] gives 256 s. Any request below 128 s gives 128 s. The
knob is not merely imprecise -- across a factor-of-two range of inputs it
returns the same answer.

Predicted vs measured, 600-step benchmark, dt ~ 1.63 s:

| | predicted | measured |
|---|---|---|
| first radiation call after t=0 | crossing into the 256 s bucket | model t = 194.1 s |
| second | crossing into 512 | 450.0 s |
| third | crossing into 768 | 707.0 s |
| gaps | 256 s | 255.9, 257.0 s |

**This is the campaign's recurring defect class, not a new one.** Every instance
has the same shape -- *an absolute single-precision value used where only a
difference was ever meaningful*:

* item 19/21: the frame's `(rho, theta)` pair carries an absolute pressure that
  is ~18 hPa wrong; only the T recovered from the pair is meaningful, because
  the two errors share the wrong p and cancel in the ratio.
* item 20: `erf_enforce_hse` anchors an absolute p = 101325 Pa at z = 0 in every
  column, when only the vertical *increment* is physically constrained.
* item 21 (retracted 19f): profile-matching read a pressure error as an absolute
  305 m height offset, and a model knob was added for it.
* this item: an absolute epoch clock differenced to get an interval.

The general lesson for an SP build: **carry differences as differences.** An
absolute time, pressure, or height in float32 has ULP proportional to its
magnitude, and every one of these bugs was invisible because the absolute value
looked plausible.

**Deliberately NOT fixed yet.** Every scored run used `180.0` and therefore ran
at 256 s, so the runs are mutually consistent and #25 stays comparable. Fixing
it changes the radiative forcing of every future run and requires rescoring the
baseline. Fix after the #25 measurement lands.

**Fix when addressed:** test against MODEL-RELATIVE time. ERF's own `cur_time`
is small and precise (the step log prints `TIME = 951.4483779`); only the
radiation interface adds the epoch offset. Either pass `cur_time` for the
cadence test and keep the epoch value solely for the solar-geometry call, or
accumulate `m_time_since_rad += dt` and compare that. Note the run-termination
test at `ERF.cpp:626` (`start_time + cur_time < stop_time`) has the same shape
and should be checked for the same quantisation.

## 27. Relaxation-band CORNER instability kills the 192x96 24-h run at 18 h

Measured 2026-07-26. The 192x96 Jan-9 run reached **18.002 h** (46,270 steps)
and then went non-finite in all 15 conserved components at once. The tripwire
(`check_for_nans_int = 10`) caught it; there was no crash.

**Where.** A single cell: `(i, j) = (0, 95)` -- the corner where the xlo and yhi
relaxation bands intersect, at 34.632, -123.470, over open ocean with 12 m of
terrain, at mid-levels (k = 32-34, z = 4585-5530 m).

At t = 18 h, along the entire xlo inflow wall:

| j | 75 | 79 | 83 | 87 | 91 | 93 | **95** |
|---|---|---|---|---|---|---|---|
| max\|w\| | 0.27 | 0.26 | 0.77 | 0.31 | 0.20 | 0.80 | **18.00** |

A factor of 22 over its own neighbours two cells away. It decays inward from the
corner (18.00 at i=0, 13.05 at i=12, 10.62 at i=24), which is the corner
spreading, not a wall-wide problem.

**How.** Monotonic growth at a FIXED location over three hours, while the argmax
before that was in the interior/east:

| t (h) | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|
| argmax \|w\| | (182,82,31) | (182,83,32) | **(0,95,34)** | **(0,95,33)** | **(0,95,32)** |
| max \|w\| | 11.31 | 11.38 | 13.05 | 15.41 | **18.00** |
| v there | -3.39 | -3.38 | +15.43 | +16.29 | +16.56 |

The corner has westerly INFLOW on x and +16.6 m/s northward flow into the yhi
face simultaneously -- the relaxation has to reconcile two specified states in
one cell. dt had been declining for twelve hours before the failure (median
1.599 at 0-3 h, 1.556 at 6-9 h, 1.400 at 9-12 h, 1.204 at 12-15 h, 1.191 at
15-18 h): a slow squeeze, not a sudden event.

**This is NOT the terrain-in-band risk flagged before the run, and NOT a
timestep problem.** The terrain-loaded band (xhi, 1204-1322 m) oscillated in
7-11 m/s for the whole run and never led. The Alamo/Liebre block at
34.676/-118.952, the specific worry, peaked at 4.63 m/s and decayed; the failure
is 412 km away from it over water.

**The mechanism is chronic and pre-existing.** Corner |w| on the completed
128x64 runs at t = 24 h:

| run | interior max\|w\| | xlo/ylo | xlo/yhi | xhi/yhi | xhi/ylo |
|---|---|---|---|---|---|
| bdyfix/wz | 2.09 | 0.17 | 4.84 | 5.60 | 7.75 |
| bdyfix/ic_hyd | 2.08 | 0.17 | 4.84 | 5.59 | 7.75 |
| bdyfix/sst_ctl | 12.47 | 0.26 | 4.77 | 5.26 | 7.92 |

The four corners hold the largest |w| in the baseline domain -- 7.75 against an
interior 2.09, a 3.7x ratio -- on runs we scored and called clean. The 192x96
domain did not create this; it moved the same corner from 4.84 to 18.00 and past
the margin. **Every scored result in this campaign was produced with corner
amplification already the dominant |w| signal in the domain.**

Related to item 14 (the relaxation zone manufactures spurious updrafts) and item
18 (the ramp-gradient term), but distinct: those are wall-wide, this is
specifically where two bands overlap and each wants a different state.

**Not yet diagnosed:** why this corner and not the other three, and why worse on
the larger domain. Candidates, untested -- the NW corner now sits 2 deg further
west in open ocean under stronger cross-corner flow; the corner is inside the
6.92% zero-filled terrain strip (lon < -123.0), though terrain there is
genuinely 0; the larger domain places the corner in a different synoptic
position relative to the storm.

**#25 IS NOT ANSWERED.** The run needed 24 h of accumulation to score against
the MRMS 24-h total. It has 18 h. No FSS comparison was made.

### 27a. How the corner weights actually combine: already `max`, not compounding

Read before changing anything. `Source/Utils/ERF_Utils.H`,
`realbdy_compute_relaxation`, four kernels:

```cpp
// bx_xlo / bx_xhi  -- "Corners with x boxes"
Real eta    = std::max(eta_lo, eta_hi);
Real Factor = std::max(xi*xi, eta*eta);      // <-- combination
rhs_arr(i,j,k,n+icomp) += Factor*F1*delta;

// bx_ylo / bx_yhi  -- "No corners for y boxes"
Real Factor = eta*eta;                        // x weight not consulted
rhs_arr(i,j,k,n+icomp) += Factor*F1*delta;
```

**It is not a product, not a sum, and not compounding sequential application.**
The x boxes OWN the corners and the y boxes explicitly EXCLUDE them, so a corner
cell is written exactly once, by the x kernel, with `max(xi^2, eta^2)`. There is
no double `+=`. The nearest-wall-weight fix is already implemented.

Note the sign convention before proposing `min`: the weight is 1 AT the wall and
0 at the band's inner edge, so the NEAREST wall is the LARGER weight. `max` IS
the nearest-wall rule; `min` would be the weaker/farther wall, which is the
opposite of what is there. These are not the same change.

**What `max` does geometrically -- this is the actual defect.** Level sets of a
max are L-shaped, so the corner square carries nested L-contours. Factor with
width = 10, NW corner:

| | j=95 | j=93 | j=91 | j=89 | j=87 | j=85 |
|---|---|---|---|---|---|---|
| i=0 | 0.902 | 0.902 | 0.902 | 0.902 | 0.902 | 0.902 |
| i=2 | 0.902 | 0.562 | 0.562 | 0.562 | 0.562 | 0.562 |
| i=4 | 0.902 | 0.562 | 0.303 | 0.303 | 0.303 | 0.303 |
| i=6 | 0.902 | 0.562 | 0.303 | 0.122 | 0.122 | 0.122 |
| i=8 | 0.902 | 0.562 | 0.303 | 0.122 | 0.022 | 0.022 |

Along a straight face F depends on ONE coordinate, grad(F) is wall-normal, and
the tangential-only fix kills grad(F).(A-B) identically. Inside the corner
square grad(F) has components in BOTH directions and flips direction
discontinuously across the diagonal xi = eta, where the max switches branch.

**This is the corner case of grad(F).(A-B) identified during the tangential-only
work: corners survive BY CONSTRUCTION, because grad(F) varies in both
directions and no wall-normal projection can remove it.** What is new here is
not that the term survives -- that was predicted -- but that it is now the
DOMINANT |w| signal in the domain rather than a residual.

### 27b. Every scored result carries this

The four corners hold the largest |w| in the 128x64 baseline (7.75 against an
interior 2.09) on runs scored as clean, and the corner term is a grad(F) source
with no physical counterpart. **If a corner fix changes the baselines
materially, every comparison in items 22-25 needs re-reading**, including the
#25 conclusion that the 3 km run scores below its own driver. The fix must
therefore be scored on the 128x64 domain too, not only on 192x96.

### 27c. Ocean elevation verified

Not a fill artifact. In `channel_islands_terrain_3km_192x96.txt`: 141,122 of
166,753 points (84.63%) are EXACTLY 0.0, with no negative values, no -0.0, and
nothing in the 1e-9 range. The 8 points in (0, 1) m and the 0.402 m minimum
nonzero are real coastal DEM values.

Correcting one number in the report above: the failure cell's "terrain 12 m" is
`z_phys` at k = 0, which is the FIRST CELL CENTRE -- terrain plus half of the
25 m initial layer. Terrain there is 0.0. The corner sits over open ocean at
exactly zero elevation.

### 27d. The corner blend is FALSIFIED, and the null result exonerates the weighting

`F = xi^2 + eta^2 - xi^2*eta^2` replacing `max(xi^2, eta^2)` (commit f54ec3a2).
Straight faces verified bit-identical (0 of 760 cells changed); only the
200-cell corner square differs, by up to 0.2461.

Result: **no effect.** The run died at 18.001 h; the unfixed run died at
18.002 h. Same corner, same curve, ~0.5% higher throughout:

| t (h) | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|
| max\|w\| at xlo/yhi, `max()` | 8.69 | 11.34 | 13.05 | 15.41 | 18.00 |
| max\|w\| at xlo/yhi, blend | 8.73 | 11.42 | 13.14 | 15.50 | 18.09 |

**The null result is the useful part.** F in the corner square changed by up to
27% at individual cells (0.562 -> 0.809 at (2,93)) and the response moved 0.6%.
That is near-zero sensitivity of the corner |w| to the relaxation weight. If the
runaway were driven by the weighting -- by its magnitude OR by the grad(F)
branch discontinuity -- a 27% change in the forcing coefficient could not
produce a 0.6% change in the response.

**This exonerates the corner weighting as the mechanism, and predicts that the
other two weighting options will also fail:** `min(F_x,F_y)` and a corner taper
both act on the same coefficient the response is insensitive to. Neither is
worth a run. Two experiments saved by a null result.

Note also that corners exceed the interior from t = 1 h in BOTH runs (ratio 1.20
at t=1, 1.98 at t=3), so "corners hold the domain max" is not a threshold the
fix has to restore -- it is a property this scheme has never had.

**Correction to 27:** the twelve-hour dt decline was reported there as the
failure signature. It is not. The blend run reproduces the decline to within 1%
(6-9 h: 1.5561 -> 1.5596; 9-12 h: 1.3998 -> 1.4108) while the corner weighting
changed, so the decline is independent of the corner and most likely tracks the
storm's intensification. The corner runaway begins at 14 h and is a separate,
later event. The dt trajectory is NOT a diagnostic for this failure.

**Where to look next.** Not the relaxation ramp. The cell (0,95) is in the
width-1 SPECIFIED region of both the xlo and yhi faces, not merely in both
relaxation bands. In `ERF_BoundaryConditionsRealbdy.cpp` the x-face and y-face
fills are separate `ParallelFor` launches writing `dest_arr`, and where their
boxes meet at a corner the same cell is written by both, so the later launch
wins and the corner silently takes one face's data while its neighbours take the
other's. That is a target inconsistency independent of F, and it is consistent
with a response insensitive to F. Verify the box overlap before assuming it.

## 28. Boundary time interpolation is quantised to 128 s: the float32-epoch class, fifth appearance

Read-only audit, 2026-07-27, prompted by three 24-h runs dying within one step of
each other at 18.001-18.002 h.

**The frames are clean.** All 9 3-hourly frames, both streams: zero non-finite
values in any of rho/u/v/w/theta/qv/qc/qr or sst/q_star/t_star/u_star/ls_mask/alb.
The 18:00 frame's frame-to-frame deltas sit in line with its seven predecessors
(theta 109 K against 51-127 across the series; u 15.8, v 20.5, w 1.97, all
mid-range). Sampled at the ERF corner cell (0,95) the forcing ramps smoothly
through 18:00 with no discontinuity:

| frame | 09Z | 12Z | 15Z | **18Z** | 21Z | 00Z |
|---|---|---|---|---|---|---|
| v at (0,95), k=20 | +11.30 | +15.06 | +22.92 | **+29.56** | +38.54 | +45.59 |

So the data hypothesis is dead: nothing is wrong with the 18:00 frame.

**The interpolation is NOT sound.** `ERF_BoundaryConditionsRealbdy.cpp:31-46`:

```cpp
Real time_tot = time + start_time;                      // ABSOLUTE epoch ~1.673e9
Real time_since_start_bdy = time_tot - start_bdy_time;  // difference of float32 absolutes
int  n_time = static_cast<int>( time_since_start_bdy / dT );
Real alpha  = (time_since_start_bdy - n_time * dT) / dT;
```

float32 ULP at 1.673e9 is 128 s, so `time_tot` is quantised to 128 s and every
quantity derived from it inherits that. With dT = 10800 s:

* alpha resolution = 128/10800 = **0.011852**
* **676** distinct (n_time, alpha) states over 24 h, against ~66,000 for a
  continuous alpha
* the boundary target therefore advances as a STAIRCASE: it holds constant for
  ~128 s, then jumps by 1.19% of the whole frame-to-frame difference in a single
  step. At the corner between 15Z and 18Z that is an instantaneous +0.079 m/s in
  v every 128 s; between 18Z and 21Z, +0.107 m/s.

This is the same defect as items 19/21 (frame rho,theta absolute pressure), 20
(erf_enforce_hse absolute 101325 Pa), retracted 19f (absolute 305 m), and 26
(rad_freq_in_time absolute epoch clock): **an absolute single-precision value
used where only a difference was ever meaningful.** Fifth appearance.

**But it does NOT explain the 18 h timing.** Every frame transition shows the
identical pattern -- alpha climbing to ~0.99, n_time incrementing, alpha
restarting near 0.01:

| transition | 3 h | 6 h | 9 h | 12 h | 15 h | **18 h** | 21 h |
|---|---|---|---|---|---|---|---|
| alpha before | 0.9956 | 0.9911 | 0.9985 | 0.9941 | 0.9896 | **0.9970** | 0.9926 |
| alpha after | 0.0074 | 0.0030 | 0.0104 | 0.0059 | 0.0015 | **0.0089** | 0.0044 |

18 h is unremarkable among them. The quantisation is a real defect worth fixing
on its own terms -- it applies impulsive boundary forcing 675 times per day, and
the impulse size scales with the frame-to-frame difference, which grows through
this storm -- but it is uniform in time and does not by itself select 18 h.

**Fix:** compute `time_since_start_bdy` from MODEL-RELATIVE time. `time` is
already small and precise; only the epoch offset destroys it. Same fix as 26,
and `ERF.cpp:626` (`start_time + cur_time < stop_time`) has the same shape and
matters for month-long segments.

**Status: the corner line is closed and the frame-data hypothesis is closed.**
What selects 18 h is not yet explained.

### 28a. The quantisation fix is correct and changes nothing: fourth death at 18.001 h

Fix applied (frame index and alpha from model-relative time; three absolute-epoch
stop comparisons in ERF.cpp/ERF_Coupling.cpp converted to elapsed spans).

**The fix is live** -- clamp-only and clamp+fix diverge at step 25 in the 8th
significant digit (TIME 17.15820269 vs 17.15820448), which is exactly the size of
perturbation a continuous-instead-of-staircase alpha should produce.

**And it changes nothing.** Fourth consecutive death at 18.0017 h, and the dt
trajectory is indistinguishable from the unfixed run:

| window | 0-3 | 3-6 | 6-9 | 9-12 | 12-15 | 15-18 |
|---|---|---|---|---|---|---|
| original | 1.5990 | 1.6201 | 1.5561 | 1.3998 | 1.2037 | 1.1907 |
| +corner clamp +28 fix | 1.5990 | 1.6200 | 1.5571 | 1.4013 | 1.2029 | 1.1898 |

Under 0.1% in every window across 18 hours. The staircase impulses were
dynamically negligible: the model integrates through them and the response is
set by the smooth trend, not by the 675 kicks.

**The fix is still worth keeping** -- it is correct, it removes a real defect
from the forcing path, and the elapsed-span stop conditions are REQUIRED for
month-long segments where the epoch magnitude is larger. But it is not the
cause, and the "threshold reading" it was meant to test is falsified with it.

### 28b. Four runs, one death time, six falsified hypotheses

| # | hypothesis | how tested | result |
|---|---|---|---|
| 1 | corner weight magnitude | max -> C1 blend, F changed up to 27% | 0.6% response, died 18.001 |
| 2 | grad(F) branch discontinuity at the diagonal | same change (C1 by construction) | same |
| 3 | double-write / later-launch-wins at corners | source read | boxes explicitly trimmed; written once |
| 4 | disagreeing per-face targets | source read | same MultiFab, bit-identical |
| 5 | bad 18:00 frame data | all 9 frames, both streams | clean; corner forcing ramps smoothly |
| 6 | float32 epoch quantisation of the forcing | fixed and rerun | fix live, no effect |

Death times: 18.002, 18.001, 18.001, 18.0017 h. dt trajectories agree to <0.1%
across all four. The failure is robust and deterministic and is insensitive to
corner weighting, corner |w|, boundary data content, and boundary time
interpolation.

Whatever selects 18 h is in the bulk solution, not the boundary machinery.
NOT INVESTIGATED FURTHER -- reporting per instruction before any next step.

## 29. The 18 h failure is an ILLEGAL DEVICE MEMORY ACCESS, not a NaN. Six hypotheses were tested in the wrong neighbourhood.

Restart from `chk46225` (step 46225, t = 17.9999 h, five steps before the death)
with `check_for_nans_int = 1`. Reproduces in ~30 s and two steps.

**Step 46226** completes normally. **Step 46227**:

```
Reading weather data 64801.21094 6 7 9          <- frame pair advances (5,6) -> (6,7)
The values of alpha1 and alpha2 are 0.9998878837 0.0001121163368
amrex::Abort::0::CUDA error 700 ... an illegal memory access was encountered !!!
```

**There is no first non-finite value.** The NaNs reported in all four 24-h runs
were downstream of memory corruption: with `check_for_nans_int = 10` the model
ran up to ten steps on a corrupted device heap before anything looked, which is
why all 15 conserved components went non-finite simultaneously across ~87% of the
domain. We were reading the debris, not the event.

**compute-sanitizer names the faulting kernel:**

```
Invalid __global__ read of size 4 bytes
  at AMReX_FBI.H:553  (FabArray<FArrayBox>::FB_local_copy_gpu)
  by thread (128,0,0) in block (0,0,0)
  Address 0xac7d163cb27cde8 is out of bounds
  and is 776685958333975529 bytes after the nearest allocation
Host: FB_local_copy_gpu <- FBEP_nowait <- FillBoundary
      <- FillPatchSingleLevel <- ERF::FillPatchCrseLevel <- ERF::timeStep
```

An address 7.8e17 bytes past the nearest allocation is not an indexing slip; it
is a garbage pointer. A MultiFab reaching `FillPatchCrseLevel` carries a
dangling or corrupt FAB in its FillBoundary copy tags, and the trigger is
deterministically the frame-pair advance to (6,7).

**This explains every observation that defeated six hypotheses:**

| observation | explanation |
|---|---|
| four runs died at 18.002/18.001/18.001/18.0017 h | triggered by a frame INDEX transition, not by dynamics |
| insensitive to corner weighting, corner \|w\|, boundary data, time interpolation | none of them touch the faulting path |
| all 15 components non-finite at once | memory corruption, not a physical instability |
| dt trajectories identical to <0.1% across all four runs | the solution was healthy right up to the fault |

The corner |w| growth from 14 h is real but is a separate, benign feature of the
storm's intensification. It was never the cause, and #27's containment, #27d's
blend, and the #28 fix were all treating a solution that was not failing.

**Ruled out already:** frame indices are in range (6, 7 of 9, guarded);
all 9 frames have the level-0 duplicate and drop uniformly to nz = 37, so no
frame-to-frame size mismatch.

**NOT diagnosed:** which MultiFab, and why the (6,7) advance specifically. This
is a memory-lifetime bug, a different class from everything investigated so far.
Reported without further investigation, per instruction.

**Prediction for the two-segment test.** The trigger is index-driven, not
elapsed-time-driven and not amplitude-driven. A fresh-init 12-h segment gets its
own frame list and never reaches idx2 = 7, so it should survive -- and if it
does, that is evidence for an index/count-driven resource bug rather than
anything physical. That makes the segment run diagnostic as well as productive.

**Keep the #28 fix regardless:** the elapsed-span stop conditions are required
for month segments at production epoch magnitudes.

### 29a. Fix attempt 1 (lifetime barrier): FAILED. What it eliminated.

`FillForecastStateMultiFabs` and `FillSurfaceStateMultiFabs` both declare
FUNCTION-LOCAL `Gpu::DeviceVector`s and capture their raw `.data()` pointers in
ParallelFors. The only `streamSynchronize` in either function sits after the H2D
copies, not after the consuming kernels, so both returned -- destroying the
buffers, which AMReX frees to the arena WITHOUT stream ordering -- while those
kernels could still be in flight.

That is a genuine lifetime hazard and the barrier is kept. **It is not this bug:**
the reproducer still faults identically.

ELIMINATED: use-after-free of those functions' own local device buffers by their
own kernels.

### 29b. The fault tracks arena layout

| arena setting | outcome |
|---|---|
| default | CUDA 700 at step 46227 |
| `the_arena_is_managed=1` | CUDA 700, still faults |
| `the_arena_init_size=1000000` (1 MB) | CUDA 700 **at step 46226** -- one step EARLIER |

Moving the fault by changing allocator layout confirms a stale/dangling pointer
whose visibility depends on what the arena hands out where. It does not identify
the owner.

Sanitizer report is stable across attempts -- same kernel, same call stack, a
different garbage address each time:

```
Invalid __global__ read of size 4 bytes at AMReX_FBI.H:553
  FB_local_copy_gpu <- FBEP_nowait <- FillBoundary
  <- FillPatchSingleLevel <- ERF::FillPatchCrseLevel <- ERF::timeStep
```

`FillPatchCrseLevel` FillBoundaries `vars_old/vars_new[0]`, which the frame-read
path does not touch -- so the frame advance is corrupting something those
MultiFabs depend on, not the forecast state itself.

**Leading hypothesis, NOT yet tested:** AMReX caches FillBoundary communication
metadata (`m_TheFBCache`) keyed by BDKey (BoxArray + DistributionMap ids + ngrow
+ periodicity). Those ids are RECYCLED when the objects are destroyed. A cached
FB whose key collides with a recycled id returns copy tags built for a different
BoxArray, and the local box indices in those tags then resolve to fabs that do
not exist -- producing exactly a garbage source pointer in `FB_local_copy_gpu`.
This is count-driven, which is why it first bites on the seventh frame read and
why it is deterministic. Next step is to find what BoxArray/DistributionMap the
frame-advance path creates and destroys per read.

**Status: NOT FIXED. The 24-h gate is not reached.**

### 29c. It is a GPU MEMORY fault, not physics. Still not fixed.

**Found: per-call BoxArray/DistributionMapping/MultiFab churn.** `strip_to_global_fab`
(ERF_WeatherDataInterpolation.cpp:725) constructs and destroys a fresh `BoxArray`,
`DistributionMapping` and `MultiFab` on EVERY call, and `ParallelCopy` registers
cache metadata keyed on those ids:

```cpp
BoxArray sba(strip);
DistributionMapping sdm(Vector<int>({0}));
MultiFab smf(sba, sdm, 1, 0);
smf.ParallelCopy(src, scomp, 0, 1);
```

It is called 4 planes x BdyEnd vars x 9 times ~= 250 times, but only at
INITIALIZATION, not per frame advance. So the BDKey-recycling shape the search
predicted is present -- just not on the per-frame path.

**Precedent already in this codebase.** ERF_SurfaceDataInterpolation.cpp:462
carries this comment, from an earlier fix in this campaign:

> "Allocating these per call -- and worse, allocating `nxt` inside the 12-sweep
> loop -- meant 12 GPU MultiFab allocations every timestep, which segfaulted a
> 24-h run at ~27.5k steps (SIGSEGV, no NaN, no assert)."

Same class, same signature, fixed there by hoisting scratch to statics.

**Fix attempt 2 (hoist the 13 per-frame DeviceVectors to reused statics): FAILED,
and made it WORSE.** The reproducer then failed EARLIER, at step 46226, with
`Kokkos ERROR: Cuda memory space failed to allocate 205.1 MiB`. Reverted.
ELIMINATED: allocation churn of the frame scratch buffers is not the cause, and
removing it costs headroom that something else needs.

**The failure is memory, not physics.** Symptoms across probes:

| probe | result |
|---|---|
| default | CUDA 700 at 46227 |
| `the_arena_init_size=1000000` | CUDA 700 one step EARLIER |
| `the_arena_is_managed=1` | still faults |
| static frame scratch | Kokkos OOM, earlier |
| `rad_ncol_chunk` 2048 / 1024 | CUDA 700, unchanged |
| restart from chk28351 (12 h) | Kokkos OOM immediately, in `Radiation::run_impl()` |

A 205 MiB allocation failing with 8.5 GB free is not a true OOM -- it is a CUDA
context already poisoned by an earlier illegal access, after which every
allocation fails. compute-sanitizer on that restart reports
`cudaErrorMemoryAllocation` inside `Radiation::run_impl()`.

**There is no leak.** GPU memory over the full 18-h run is flat: mean 7.7-8.2 GB
per octile, peak 10617 MiB of 16376, growth first-to-last octile **+206 MiB**.
5.8 GB free at peak. So exhaustion-by-accumulation is ruled out.

**Symptom is not run-to-run stable** (CUDA 700 vs Kokkos OOM for the same build
and inputs), consistent with item 13's finding that this SP/GPU fork is not
bit-reproducible.

**Kept:** the lifetime barrier from attempt 1 (real use-after-free window).
**Reverted:** the static frame scratch.

**STATUS: NOT FIXED. The 24-h gate is not reached.** Two attempts, both
falsified, both recorded. The remaining strong lead is the init-time
BoxArray/DistributionMapping churn in `strip_to_global_fab` -- ~250 create/destroy
cycles whose recycled BDKeys can collide with cached FillBoundary metadata for
the long-lived `vars_old`/`vars_new`. That is consistent with the sanitizer stack
(`FB_local_copy_gpu` under `FillPatchCrseLevel`) and with a garbage source
pointer. Untested fix: give `strip_to_global_fab` a single persistent
BoxArray/DistributionMapping/MultiFab reused across all ~250 calls.

### 29d. CORRECTION: 29c's memory conclusions were contaminated. Attempt 3 falsified.

**A leftover Docker container held 7072 MiB of GPU memory for 11 hours.** It was
a survivor of a `timeout`-killed compute-sanitizer run: `timeout` kills the
client, not the container. With 7.2 GB unavailable against a ~10.6 GB peak, the
machine could not run the model at all.

**Every OOM result in 29c is therefore invalid**, including:

* attempt 2 (static frame scratch) "made it worse, Kokkos OOM" -- untested;
* the chk28351 restart failing immediately in `Radiation::run_impl()`;
* the `rad_ncol_chunk=5000` Kokkos OOM;
* the inference that "the failure is memory, not physics" and that a poisoned
  context explained a 205 MiB allocation failing with 8.5 GB free. **There was
  not 8.5 GB free.** That reading was wrong, and 29c is retracted on that point.

The one 29c finding that survives is from an uncontaminated source: GPU memory
over the real 18-h run is flat, +206 MiB first-to-last octile, peak 10617 MiB.

**On a clean GPU the original signature is exactly reproducible again:**
step 46226 completes, `CUDA error 700` on the (6,7) frame advance at step 46227,
`Compressible dt = 1.224641323` throughout -- a healthy solution up to the fault.

**Attempt 3 (persistent box-keyed scratch in `strip_to_global_fab`): FAILED, and
it is incorrect.** One persistent MultiFab per distinct box means the SAME
MultiFab serves T, QV, RHO, QC and QR -- all cell-centred, all the same strip
box -- across all 9 frame times. The state is corrupted immediately on restart:
`Compressible dt` collapses from 1.2246 s to **0.0093 s**, a 130x drop, before
any fault. Reverted.

ELIMINATED: not the mechanism, and box-keyed persistence is the wrong shape --
any persistent-scratch fix here must key on (box, variable, time) or explicitly
clear between uses, or it cross-contaminates the boundary planes.

**Process lesson, and it is the expensive one this session.** Three consecutive
"results" were environmental. The tell was available and ignored: the symptom
changed across builds on identical inputs (CUDA 700 -> Kokkos OOM), and item 13
already establishes that this fork is not bit-reproducible, so a CHANGED FAILURE
MODE should have prompted an environment check before a code conclusion. Check
`nvidia-smi --query-compute-apps` and `docker ps` before interpreting any GPU
memory symptom.

**STATUS: NOT FIXED.** Three attempts, all falsified, all recorded. The
strip_to_global_fab churn is real and remains the leading suspect, but the fix
must not reuse one buffer across variables.

### 29e. BDKey collision is DEAD. Flushing AMReX's caches changes nothing.

AMReX's `flushFBCache()` / `flushCPCache()` are private, so a hash-gated local
patch exposes them (`apply_amrex_fbcache_patch.sh`, same machinery as the RRTMGP
fix, wired into the build). `erf.realbdy_flush_fb_cache=1` drops the entire
FillBoundary and ParallelCopy metadata cache immediately after boundary-plane
setup, forcing every later FillBoundary to rebuild its copy tags from the
CURRENT BoxArray.

Reproducer from `chk46225`, verified-clean GPU, flush confirmed firing in the log:

| | last step | fault |
|---|---|---|
| `realbdy_flush_fb_cache=0` | 46226 ends | CUDA 700 at 46227 |
| `realbdy_flush_fb_cache=1` | 46226 ends | CUDA 700 at 46227 |

Identical. **A stale cache entry seeded by a recycled BoxArray/DistributionMapping
id is NOT the mechanism**, and with it dies the `strip_to_global_fab` churn as a
suspect -- the churn is real, but it is not what produces the garbage pointer.

That also disposes of the timing objection raised against it: the churn is
init-only and never explained why the fault waits for frame 7. It does not have
to, because it is not the cause.

**Eliminated so far, on this bug:** corner weight magnitude; grad(F) branch
discontinuity; corner double-write; per-face target disagreement; bad frame data;
float32 epoch quantisation of the forcing; frame-scratch lifetime; frame-scratch
allocation churn; BDKey cache collision. Nine.

**What survives:** an invalid 4-byte global read in `FB_local_copy_gpu` under
`ERF::FillPatchCrseLevel`, deterministically on the (6,7) frame advance, with a
healthy solution (`dt = 1.2246 s`) up to the fault and a garbage source pointer
(offsets of order 1e17-1e19 bytes). Since the tags are rebuilt fresh after a
flush and are still wrong, the corruption is in what the tags are built FROM --
the FabArray's own fab pointers -- not in cached metadata describing them.

### 29f. Structural guard against environment contamination

`gpu_preflight.sh` refuses to start any experiment if a stray erf-hindcast
container is running or if the GPU holds more than a desktop baseline, printing
the holders. `--clean` kills strays first. Wired into every experimental
invocation. This is the 29d lesson made structural: `timeout` kills the docker
client, not the container, and eleven hours of results were invalidated by one
survivor.

### 29g. Path trace: nothing found that satisfies the timing constraint

Read-only trace of what touches `vars_old`/`vars_new[0]` fab storage around a
frame advance. The filter throughout: it must fire on the (6,7) advance and NOT
on (0,1) through (5,6). Six clean advances, then a fault.

| checked | finding | verdict |
|---|---|---|
| growth/resize/push_back in both interpolation TUs | every resize is at INIT (`bdy_data_*` to `ntimes`, `next/last_read_forecast_time` to `nlevs`); nothing grows per advance | no count-triggered reallocation |
| `bdy_data_xlo[n_time][ivar]` indexing | outer dim is `ntimes` = 9, inner is `BdyEnd` = 7. `n_time_p1` = 7 indexes the OUTER (9) dim -- in range | the `BdyEnd == 7` coincidence is only a coincidence |
| `ivar` range against the 7-element inner dim | `ivar` comes from `cons_map` at the read components only: RHO(4), T(2), QV(3), QC(5), QR(6). Max 6 | in range, no OOB write to corrupt neighbours |
| surface per-read `resize(1)` blocks (`sst_lev`, `tsk_lev`, `lmask_lev`, `alb_lev`) | all guarded `if (empty())`; allocate once, refresh contents after. Comment confirms the surface layer holds raw pointers and the allocation must persist | fires on read 1, not read 7 |
| `define`/`clear`/`resize`/`swap`/`move` of `vars_new[0]` on the advance branch | none | -- |
| alias/non-owning views over `vars_new` | `ERF_MakeNewArrays.cpp:880` builds `MultiFab(cons_mf, make_alias, 0, ncomp_cons)` into the MRI integrator state, so the integrator holds an alias sharing `vars_new[lev][Vars::cons]`'s fab pointers. Real latent hazard if the parent is ever redefined | called at init/regrid, NOT on the advance -- does not satisfy the timing filter |

**Nothing surfaced that fires on the seventh advance and not the first six.**

The two-deep-cache idea checked out negative specifically: the frame containers
are sized once to `ntimes` at init and only indexed thereafter; there is no
small container growing per read whose reallocation could invalidate a held
pointer.

What remains true and unexplained: an invalid 4-byte global read in
`FB_local_copy_gpu` under `FillPatchCrseLevel`, deterministic on the (6,7)
advance, with `dt = 1.2246 s` healthy up to the fault, a garbage source pointer,
and no improvement from flushing the metadata caches -- so the FabArray's own fab
pointers are corrupt, from a writer this trace has not found.

**Ten hypotheses eliminated. Recommend deciding between the segment-length
workaround and handing this upstream rather than an eleventh guess.**

### 29h. Hunt for the WRITE: prime suspect ruled out, and the tool cannot see it

Reframed per instruction: look for an out-of-bounds WRITE into neighbouring
allocation, not the read we already know about.

**`strip_to_global_fab` is definitively ruled out -- on frequency, not on
plausibility.** Its only caller is `fill_bdy_data_from_hindcast()`, which has
exactly two call sites: `ERF.cpp:2240` (restart) and `ERF.cpp:2398` (init). Both
are one-time, before any stepping. It does NOT run per advance and therefore
cannot corrupt anything at the seventh one.

**compute-sanitizer STRUCTURALLY CANNOT see this write.** `AMReX_Arena.cpp:427`:

```cpp
the_arena_init_size = Gpu::Device::totalGlobalMem() / numDevicePartners() / 4L * 3L;
```

AMReX pre-allocates 3/4 of device memory as ONE `cudaMalloc` and sub-allocates
from it. An overrun from one sub-allocation into the next is INSIDE a valid
allocation, so memcheck reports nothing. This is why eleven sanitizer runs have
only ever shown the far-out-of-bounds READ (the eventual dereference of the
already-corrupted pointer) and never a write. The instrument is blind to the
hypothesis.

**Loop-bound audit of the per-advance path** (`FillForecastStateMultiFabs`), all
negative:

* `ParallelFor(gbx, ...)` writes `fine_cons_arr` over `mfi.growntilebox()` of the
  cons MultiFab itself -- in bounds by construction.
* `fine_latlon_arr` is written over the same `gbx` but is a DIFFERENT MultiFab,
  defined with `ngrow = src0.nGrow()` (isotropic int) while cons uses
  `ng = nGrowVect()` (IntVect). They are sized by different code paths -- a real
  smell, and `ngrow` is a DEAD STORE in the cons loop at
  `ERF_MakeNewArrays.cpp:397` -- but isotropic `ngrow` is >= cons's per-direction
  ghosts, so latlon is over-provisioned, not under.
* Face writes over `mfi.tilebox(IntVect(1,0,0))` reach i = hi+1, which lands in
  the ghost region of an array with ng >= 1. In bounds.
* nz after the level-0 drop: the drop erases `nxy` from every field AND one entry
  from `zvec_h`, so `nz = 37` and the buffers sized `nx*ny*nz` agree. No 37/38
  mismatch.
* `strip_to_global_fab`'s ParallelCopy: `strip` IS `dest.box()`, 1 component both
  sides, so extents and counts are equal by construction.
* Host->device writes: `Gpu::copy` (synchronous) in `strip_to_global_fab`;
  `copyAsync` + `streamSynchronize` in the frame fill.

**Layout perturbations move the fault; serialization does not fix it.**

| perturbation | result |
|---|---|
| default arena (12 GB slab) | CUDA 700 at 46227 |
| `the_arena_init_size` 1 MB | CUDA 700 one step EARLIER |
| 8 MB / 128 MB, no sanitizer | CUDA 700, no steps |
| 8 MB + compute-sanitizer | **passes 46227**, zero invalid accesses |
| `CUDA_LAUNCH_BLOCKING=1` | CUDA 700 at 46227 -- NOT a launch race |

The sanitizer "pass" is memcheck's allocation padding changing adjacency, the
same class of effect as the arena knobs -- not evidence of a race, since
serializing launches does not help. Adjacency-dependence is confirmed; the
writer is not localized.

**Eleven hypotheses eliminated. The blocker is now tooling, not ideas:** the
write cannot be caught until each buffer is individually bounded. The concrete
next step is a diagnostic AMReX build whose arena performs a direct `cudaMalloc`
per allocation instead of slab sub-allocation, using the same hash-gated patch
machinery -- then memcheck bounds every buffer and the first bad write surfaces.

## 29i. MECHANISM: it is a USE-AFTER-FREE on a recycled arena block, not an overrun

Established by making every allocation its own exact-sized `cudaMalloc`.

**The patch.** `CArena` sub-allocates from 8 MB hunks
(`AMReX_CArena.H:161`), and `AMReX_CArena.cpp:81` sizes each new block
`N = max(m_hunk, nbytes)`. A hash-gated diagnostic patch adds
`amrex.the_arena_hunk_size` (default 0 = stock, verified bit-identical): with a
tiny hunk, every allocation gets its own exact-sized `cudaMalloc` and every free
a real `cudaFree`, so there is neither sub-allocation nor block REUSE.

**Two measurements settle it:**

| run | result |
|---|---|
| hunk=64, `init_size=0`, **with** compute-sanitizer, to 46227 | **zero invalid accesses**, step 46227 completes |
| hunk=64, `init_size=0`, **no** sanitizer, to 46400 | **reaches 46400**, 173 steps past the fault, zero CUDA 700 |

The second is the important one: it is not the sanitizer's padding and it is not
bounding. **Eliminating arena block reuse removes the fault**, and with every
buffer individually bounded there is **NO out-of-bounds write to find**.

**Mechanism.** Something retains a device pointer into an arena block after that
block is freed. CArena returns the block to its free list and hands it to a later
allocation, whose writes overwrite the retained data. The stale holder is then
used -- specifically the fab pointers that `FB_local_copy_gpu` dereferences under
`FillPatchCrseLevel` -- producing a garbage address. With reuse eliminated the
block is never handed out and never overwritten, so nothing breaks.

**This explains every observation, including the ones that killed earlier models:**

| observation | use-after-free explanation |
|---|---|
| deterministic step number | fixed alloc/free sequence; collision occurs at a fixed point |
| garbage address varies run to run | whatever the new owner of the block wrote |
| `the_arena_init_size` moves the fault | changes which block is recycled when |
| cache flushing changes nothing | the pointers are overwritten, not the tags |
| `CUDA_LAUNCH_BLOCKING=1` changes nothing | host-side free/realloc ordering, not a launch race |
| per-allocation bounding finds NO write | correct -- there is no overrun |
| no reuse => no fault | the defining test |

It also explains why 29a's lifetime barrier did not help: that fixed
kernels-still-in-flight for two functions, whereas this is a STORED pointer
outliving its allocation.

**NOT YET IDENTIFIED: which pointer.** The mechanism is established; the specific
retained pointer is not. `hunk=64` is a DIAGNOSTIC, not a configuration -- it
masks a live use-after-free and must not ship.

**Next step, concrete:** binary-search the holder by keeping candidate
allocations alive. AMReX's `Arena::free` is the choke point; instrumenting it to
log (ptr, size) alongside the FabArray fab pointers at each frame advance would
name the block whose reuse coincides with the corruption. The reproducer is 30
seconds, so this is bisectable rather than speculative.

## 29j. THE POINTER IS NAMED. Culprit: `FabArray::clear()` frees fab data ~25 lines before its only stream sync

The search 29i called for was run: `CArena::alloc_protected` and `CArena::free`
instrumented to log `(seq, op, ptr, size, 4 resolved callers)` for a size window
(`apply_amrex_carena_trace_patch.sh`), plus a validator on the cached FillBoundary
tag vector (`apply_amrex_fbi_tagcheck_patch.sh`). Both are diagnostics; neither
ships.

**First, two dead ends killed by measurement, not argument.**

The "512-byte size class" that framed the search was an artifact. The sanitizer
reports `Address 0xfc87c88ebce80728 ... is 18196602103903290665 bytes after the
nearest allocation at 0x775a92a00000 of size 512 bytes`. That address is wild
(high bits set); the "nearest allocation" is simply the closest thing memcheck
knew about, 1.8e19 bytes away, and has no relationship to the fault. There are no
arena records for it at all -- it is not even a CArena sub-allocation. **The real
object is 4352 bytes.** Filtering to 512 would have missed it entirely.

The victim is not stale, it is SCRIBBLED. On every cache hit the validator
rebuilt the tags from the current fabs and compared against the device buffer:

```
ERFTAGCHECK fb_id=0 ntags=22 tag=0 sizeof_tag=176 first_diff_byte=0
  d_tags=0x7faa83133900 fabarray=0x577f959ac640 nbytes=3872
  device dfab.p=0x43905337439051b1 sfab.p=0x4390036043900360 dindex=1133511520
  fresh  dfab.p=0x7faa473c3f00     sfab.p=0x7faa48460700     dindex=0
```

`0x43905337`, `0x439051b1`, `0x43900360` are float32 **288.65, 288.64, 288.03** --
field data written over the fab-pointer array. It fires immediately after
`Reading surface data 64801.21094 6 7 9`, the frame pair advancing to (6,7): the
seventh advance, exactly the documented timing.

(Caveat for whoever reads the raw output: `sizeof(Array4CopyTag<float,float>)` is
176 with a **4-byte padding hole at offset 68**. `TagVector::define` memcpy's the
struct including padding, so 26 of the 27 reported mismatches say
`first_diff_byte=68` and are noise. Only `first_diff_byte=0` is real.)

**The block's history**, from the arena trace, symbols via `addr2line`:

| seq | op | size | call site |
|-----|----|------|-----------|
| 3287 | ALLOC | 768 B | `FabArray::setVal` <- `compute_wall_flux_correction` <- `ERF::fill_from_realbdy` <- `FillIntermediatePatch` |
| 3291 | FREE | 768 B | **`FabArray::~FabArray()`** <- `ERF::fill_from_realbdy` <- `FillIntermediatePatch` <- `advance_dycore` |
| 3335 | ALLOC | 4352 B | `FB_get_local_copy_tag_vector` <- `FB_local_copy_gpu` <- `FBEP_nowait` <- `ERF::advance_microphysics` |
| -- | never freed | | live at the fault |

**The holder** is `FabArray::m_fb_local_copy_handler` (`AMReX_FabArray.H:1594`), a
`map<uint64_t, unique_ptr<TagVector>>` populated by
`FB_get_local_copy_tag_vector` (`AMReX_FBI.H:491`). It caches raw `Array4` fab
pointers for the **lifetime of the FabArray** and is cleared only in
`FabArray::clear()`. That is why the block never reappears in the trace, and why
one scribble is permanent rather than transient: every later `FillBoundary` on
that MultiFab re-reads the same corrupted tags.

**The defect**, `AMReX_FabArray.H:1942-1982`:

```cpp
template <class FAB> void FabArray<FAB>::clear () {
    ...
    for (auto *x : m_fabs_v) {
        if (x) { nbytes += amrex::nBytesOwned(*x);
                 m_factory->destroy(x); }   // <-- line 1952: fab data -> arena free list.
    }                                       //     NO stream synchronization.
    ...
    m_fb_local_copy_handler.clear();        // <-- line 1977: ~TagVector -> undefine()
                                            //     -> Gpu::streamSynchronize().
```

The only synchronization in the destructor path runs **after** the fab memory has
already been returned to the arena. `FillBoundary`'s local-copy path is
fire-and-forget on the stream -- `FB_local_copy_gpu` launches and returns, and it
is not inside an `MFIter`, so it gets none of `~MFIter`'s all-stream sync
(`AMReX_MFIter.cpp:249-255`). So a short-lived MultiFab that ends with a
`FillBoundary` -- exactly what `compute_wall_flux_correction` does at
`ERF_InteriorGhostCells.cpp:219`, on a temporary owned by `fill_from_realbdy` --
releases memory a live kernel is still writing to.

**This accounts for every prior observation:**

- deterministic step number -- allocation order is deterministic, so which block
  lands under which cache is deterministic. Not a race.
- `hunk=64` + `init_size=0` makes it vanish (29i) -- no reuse, so the freed block
  is never handed to the tag vector.
- eleven compute-sanitizer runs found no out-of-bounds write (29h) -- there isn't
  one. The write is perfectly in bounds of a block that was legitimately freed.
- the victim is a fab-pointer array -- it is the cached TagVector.
- the lifetime barriers added to `FillForecastStateMultiFabs` and
  `FillSurfaceStateMultiFabs` did not help. They are in ERF's interpolation
  functions; this is in AMReX's FabArray destructor. Both barriers are still
  correct and stay.

**Not yet proven:** that the specific kernel writing those 288 K bytes is the
`FB_local_copy_gpu` launched on the temporary. What is proven is the release site,
the ordering defect, and that the buffer is overwritten with float field data
rather than going stale. The candidate writers all sit in the same window.

**Standing constraint reminder:** no upstream activity. This is an AMReX defect
and it stays in this file. Both diagnostic patches revert before the 24-h gate.

### 29j-bis. VERIFICATION FAILED: the `clear()` ordering defect is REAL but is NOT the cause

Per the step-1 plan, `FabArray::clear()` was patched to hoist
`Gpu::Device::synchronize()` above the fab-release loop
(`apply_amrex_clear_sync_patch.sh`, hash-gated). If the mechanism in 29j were the
cause, the fault had to disappear.

**It did not.** Same death at step 46227, `CUDA error 700`, and ERFTAGCHECK still
reports exactly one real mismatch (`first_diff_byte=0`) alongside the 26 padding
artifacts.

The instrument was verified live before accepting the negative -- this campaign
has been burned by unverified results six times on this item:

- `~FabArray()` routes through `clear()` (`AMReX_FabArray.H:2187`).
- The build recompiled all 249 objects, so the header change was picked up.
- `Device::synchronize()` is a real `cudaDeviceSynchronize()` +
  `streamSynchronizeAll()` (`AMReX_GpuDevice.cpp:844-852`), not a no-op.
- **Decisive:** the abort now surfaces at `AMReX_GpuDevice.cpp:848`, which IS the
  `cudaDeviceSynchronize()` in the patched path. The sync is executing.

**What this rules out, and what it leaves.** If no kernel can be in flight when
the block is released and the tag buffer is still overwritten, then the
corrupting write happens *after* the free. It is a stale pointer captured and
launched later, NOT a live kernel racing a destructor. Every "in-flight kernel"
framing of item 29 -- including the two lifetime barriers already in the tree --
is aimed at the wrong half of the problem. (Those barriers are independently
correct and stay; they are just not this.)

The block-provenance chain in 29j stands as measured. What does not stand is the
inference from "previous owner of the block" to "writer of the bytes". The
previous owner was identified correctly; it is simply not who wrote.

**Next instrument** (not yet run, needs sign-off per the standing no-new-hunt
rule): a device-side canary written into the tag buffer immediately after
`TagVector::define`'s `htod_memcpy_async`, checked at each cache hit. That
brackets the write to a single call and can be bisected within a step, naming the
writer directly instead of inferring it from block provenance.

### 29j-ter. Upstream note: two real AMReX defects, not ERF-specific

Not filed -- standing constraint is no upstream activity, and this stays local.
Recording so it is not lost, because both are genuine and neither depends on
anything ERF does:

1. **`FabArray::clear()` frees before it syncs.** Fab storage goes back to the
   arena at `AMReX_FabArray.H:1952` with no stream ordering; the only
   synchronization in the destructor path is `m_fb_local_copy_handler.clear()`
   ~25 lines later (`~TagVector` -> `undefine()` -> `Gpu::streamSynchronize()`).
   Any GPU work outstanding on those fabs is racing the free.

2. **`FB_local_copy_gpu` gets no `MFIter` sync.** `FillBoundary`'s local-copy path
   launches and returns; it is not inside an `MFIter`, so it never sees
   `~MFIter`'s all-stream synchronization (`AMReX_MFIter.cpp:249-255`). A
   short-lived MultiFab whose last operation is a `FillBoundary` therefore has
   outstanding device work at destruction by construction.

Together these are a latent use-after-free for any temporary MultiFab ending in a
FillBoundary. That it is not what kills THIS run does not make it safe; on a
year-long hindcast it is exactly the class of thing that surfaces at month seven.
Worth reporting upstream when the campaign's no-upstream freeze lifts.

## 29m. ROOT CAUSE, FIXED: `MultiFab::Copy` writes 82944 bytes outside the SST fab on restart

Item 29 is **not** a use-after-free. It is a plain out-of-bounds write, and every
lifetime hypothesis (29a-29l) was chasing the wrong class of bug.

**The defect**, `ERF_SurfaceDataInterpolation.cpp`:

```cpp
MultiFab::Copy(*sst_lev[lev][0], surf_mf, 1, 0, 1, surf_mf.nGrowVect());
```

The ghost vector comes from the SOURCE. A fresh start allocates `sst_lev[lev][0]`
a few lines above with `surf_mf.nGrowVect()` = (4,4,4), so source and destination
agree. **The restart path allocates it first** -- `ERF_Checkpoint.cpp:851`,
`ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;` -- and the
`if (sst_lev[lev].empty())` branch is then skipped. Destination ghosts are
(4,4,0); the copy still addresses k = -4..4.

Measured geometry across the three calls in one reproducer run:

```
2 x  surf_ng=(4,4,4)  sst_ng=(4,4,0)   <-- MISMATCH
1 x  surf_ng=(4,4,4)  sst_ng=(4,4,4)
```

The destination fab is 72x72x1 = 20736 B. `kstride` = 72*72*4 = 20736 B, so k=-4
indexes **82944 bytes BELOW** `p`, and k>0 the same distance above. The arithmetic
closes exactly: the second-set fabs start at `0x775a47138800`; minus 82944 is
`0x775a47124400`, and all seven corrupted blocks (`0x775a4712b900` ...
`0x775a47133900`) lie inside `[0x47124400, 0x47138800)`.

**The victim** is the CACHED FillBoundary tag vector of `vars_new[0][cons]`
(`m_fb_local_copy_handler`, `AMReX_FabArray.H:1594`), which holds raw fab pointers
for the FabArray's lifetime. Its `Array4CopyTag` array is overwritten with SST
floats, so `dfab.p` becomes `0x439051b1`/`0x43905337` = 288.64 / 288.65 K. Because
the tag vector is cached and never rebuilt, the damage is permanent: the next
`FillPatchCrseLevel` -> `FillBoundary` -> `FB_local_copy_gpu` dereferences a
temperature as a pointer -> CUDA 700.

**Why nothing caught it for eleven diagnostics.**

- AMReX guards exactly this: `BL_ASSERT(dst.nGrowVect().allGE(nghost) && ...)`,
  `AMReX_MultiFab.cpp:205`. `BL_ASSERT` is **compiled out in Release**. A Debug
  build would have aborted on the first frame with the right message.
- compute-sanitizer is blind to it: with 8 MB `CArena` hunks the write stays
  inside a valid `cudaMalloc`. Eleven memcheck runs found nothing (29h).
- `hunk=64` "fixing" it (29i) was a red herring -- it changed which allocation sat
  in the blast radius, not whether the overrun happened. That single false signal
  is what created the use-after-free model and cost 29i-29l.
- The fault is deterministic because arena layout is deterministic.

**The fix** (3 lines): copy only what the destination owns.

```cpp
IntVect ng_sst_copy = surf_mf.nGrowVect();
ng_sst_copy.min(sst_lev[lev][0]->nGrowVect());
MultiFab::Copy(*sst_lev[lev][0], surf_mf, 1, 0, 1, ng_sst_copy);
```

Not chosen: reallocating `sst_lev` to match. The comment above that block warns
that `m_SurfaceLayer` "stores raw pointers at construction, so the allocation must
persist and only its contents change" -- reallocating would dangle those.

**Verified.** Reproducer restart from `chk46225` to `max_step=46260`, i.e. 33
steps past the fault that killed five 24-h runs at 46227: `cuda700=0`, and the
tag-vector canary reports **zero** modifications (previously one real scribble at
the (6,7) advance).

**How it was found**, since inference failed repeatedly and measurement did not:
a canary comparing each cached tag buffer against its own `h_buffer` on every
cache hit, bracketing the write between the last clean and first dirty call site,
then bisecting that window with probes -- `C:after-surface` -> `S2:after-sst-lmask`
-> `T1:after-copy-sst`. Note the trap this avoided: the SST *sanitize* ParallelFor
20 lines below contains the literal `Real(288.0)`, matching the corrupting bytes
perfectly. It is not the writer. Stopping at the matching literal would have named
the wrong statement for the third time on this bug.

**Lesson for the campaign.** `hunk=64` was treated as a mechanism-confirming
measurement in 29i. It was a mechanism-*perturbing* one. A change that makes a
symptom disappear constrains the mechanism far less than it appears to, because it
also relocates every allocation. Prefer instruments that OBSERVE (canary, tag
check) over instruments that PERTURB (hunk size, added syncs).

## 30. GATE BLOCKER: +20%/day global mass gain. Supersedes the corner framing entirely

The 24-h gate run (fresh, post-29m) died at step 46181, **18.003 h**, first NaN in
density at (0,22,7) -- the xlo boundary -- spreading to 300 cells across i=0..39
and all 15 conserved components. `cuda700=0`, compute-sanitizer `0 errors`,
clean exit. **The item-29 memory fault is genuinely gone. This is a different,
larger problem that it was hiding.**

### The drift

| | t=0 | t=18 h | change |
|-----|-----|--------|--------|
| interior mean rho | 0.8685 | 1.0166 | **+17.1%** |
| global mass (yt)  | 7.637e5 | 8.948e5 | **+17.2%** |
| solver `MASS`     | 1.5957e15 | 1.9239e15 | **+20.6%** |

Monotonic from the first hour. Two independent measurements agree. Interior
column gain is **17.06% +/- 0.24%** -- uniform to 1.4% relative across the whole
interior, so it is NOT a boundary influx advecting inward (that leaves gradients)
and NOT terrain-slope pumping (that concentrates over slopes). Gain is present at
every level, 15.6% through the low/mid column rising to 36.6% at the lid.

**A 20%/day mass gain corrupts every scored result from this deck, including runs
that do not crash.** It is upstream of #25, of the corner work, and of item 29.

### Term attribution (all by measured null, not by argument)

Protocol: 800 steps, `MASS` vs `TIME` from the solver's own diagnostic.

| config | drift | verdict |
|--------|-------|---------|
| baseline | **+2.374 %/h** | -- |
| `hindcast_wall_flux_correction=false` | +2.401 | excluded |
| `hindcast_zhi_sponge_damping=false` | +2.375 | excluded |
| `w_damping=false` | +2.374 | excluded |
| `moisture_model="None"` (Morrison off) | +2.373 | excluded |
| `hindcast_blend_bdy_theta=0` | +2.374 | excluded |
| `cfl` 0.2 -> 0.1 (2x the steps) | +2.545 | **not per-step; not SP rounding** |

The cfl test matters: a per-step accumulation bias would roughly DOUBLE the
drift per unit time when dt halves. It moved 7%. So this is a continuous-time
source, single-precision rounding is excluded, and **no DP build is needed** to
establish that.

(Early rate +2.374 %/h vs the 18-h average ~0.86 %/h: there is a fast spin-up
transient on top of a steady drift. Toggle comparisons are all same-window.)

### Leading mechanism: net lateral boundary mass flux is UNENFORCED in the compressible path

The anelastic control could not run -- it aborts in the initial projection:

```
Projecting initial velocity field at level 0
 TOTAL INFLUX / OUTFLOW 238790832 234497616
Erroneous arithmetic operation
```

influx exceeds outflux by 4.29e6, a **1.8% imbalance in the boundary data**, and
the projection dies trying to remove it. `enforceInOutSolvability`
(`ERF_ConvertForProjection.cpp:325`) exists precisely to rescale outflow by
`alpha_fcf = influx/outflux` so net flux is zero -- **and it is called only from
the projection path, i.e. only when anelastic. The compressible path has no
equivalent.** A net influx therefore just accumulates.

This is consistent with every measurement: uniform in space (acoustic waves cross
576 km in ~28 min, fast vs 18 h, so an imbalance equilibrates into a uniform
compression), continuous-time, immune to every physics toggle. It also explains
why the ONLY clean mass result on record (RUNBOOK ~line 329, "rho sum
bit-identical over 24 h") was measured on a **double-precision, `anelastic=1`**
build -- a configuration that ENFORCES solvability, and in which rho is fixed to
rho0 anyway. That note is vacuous as evidence for the compressible path.

**NOT YET CONFIRMED.** Direct flux integration from the gate plotfiles did not
reproduce the steady drift: net flux came out -0.228 / +0.781 / -0.128 %/h at
4/9/18 h, fluctuating in sign against a steady +0.86 %/h. That estimate uses
cell-centred velocities from the outermost cell instead of true face fluxes and
ignores terrain-following dz, so it is too crude to settle magnitude (the mass
integral itself checks out to 0.7% against the solver). **The decisive next
measurement is to call `compute_influx_outflux` every step in the compressible
path and log it** -- the model's own definition, no reconstruction.

### Why the wall flux correction does not save this

It is enabled and firing (`cons_only=false` at `ERF.cpp:1383`), but measurably
does nothing (+2.401 vs +2.374 with it off), and it has a real metric defect:

```cpp
const Real dz = geom.CellSize(2);      // UNIFORM
M  += rm(i,j,k,Rho_comp) * dz;
Mt += rt(i,j,k,Rho_comp) * dz;
```

This domain is stretched (ratio 1.09311734, dz_k spans ~55x) and
terrain-following. The constant cancels in `Mt/M`, so amplitude is fine, but every
level is weighted EQUALLY where true column mass weights by the actual `dz_k`. It
constrains an unweighted density sum to the ERA5 target, not column mass -- and
drives it exactly, which is how a metric error becomes persistent drift. Worth
fixing on its own merits; it is not the remedy for 20%/day, and it only touches
the outermost ring in any case.

### This supersedes the corner framing

Band-maximum rho sits at **(0,95,0)** -- the exact corner the original 24-h runs
blamed -- climbing 1.3536 -> 1.4143 over seven hours while theta there falls
288.4 -> 286.6. **The corner clamp (27e) constrained |w| and did nothing about rho
accumulation at that cell.** It moved the failure rather than removing it, exactly
as suspected. The corner was never the mechanism; it is where a global mass drift
concentrates first.

### 30a. FALSIFIED: the boundary flux imbalance does NOT carry the drift (and a correction)

Instrumented the boundary flux budget every step with the model's own convention
(`erf_bdy_mass_flux_diag`, gated by `ERF_MASSFLUX_DIAG`), logging both volume flux
and rho-weighted mass flux with the true metric face areas `ax`/`ay`, then
integrated net mass flux and compared against the solver's own `MASS`:

| quantity | value |
|----------|-------|
| integrated net boundary mass flux | 3.370e10 kg |
| measured mass gain, same window | 1.303e13 kg |
| **ratio** | **0.0026** |

**The boundary carries 0.26% of the gain -- off by a factor of ~390.** The
imbalance is real (mass in/out 2.391e8 / 2.342e8 = **+2.07%**) but far too small
to matter: 4.8e6 kg/s against a 1.6e15 kg domain is ~1e-5 %/h, four orders below
the observed 2.374 %/h. Enforcing solvability in the compressible path would not
have fixed anything, and shipping it would have been a plausible-looking no-op.

**CORRECTION to 30.** I claimed `compute_influx_outflux` carries no density and so
enforces only a VOLUME balance, making `enforceInOutSolvability` untransferable to
the compressible path. That was wrong. The function is called on `vels_vec` AFTER
`ConvertForProjection` has converted velocity to rho-weighted momentum, so it is
already a mass-like flux. Confirmed numerically: this diagnostic's mass flux
(2.391e8 / 2.342e8) matches the projection's printed `238790832 / 234497616`
almost exactly, while its volume flux is 2.3x larger. **The function is more
transferable than 30 said** -- it simply is not the remedy for this defect.

### 30b. Therefore: mass is created in the INTERIOR

With the boundary excluded quantitatively and all seven physics terms excluded by
measured null, the only surviving reading is that **the compressible dycore's
density update is not conservative in this configuration**. Consistent with every
measurement: uniform to +/-0.24% in space, continuous-time (cfl halving moved it
7%, not 2x), immune to every option toggle.

Note the ERA5 net-flux check at R ~ 6.4 mm/s is CONSISTENT with this: the driving
data being mass-consistent under ERF's own operator is exactly what you expect
when the defect lives downstream of the boundary entirely.

### 30c. Scope: this is upstream of #25 and invalidates prior scoring

**Every scored result from this compressible deck was measured on a run inflating
at ~20%/day**, including runs that completed without crashing. That includes the
FSS comparison against ERA5 in #25: both the model field and any derived
precipitation diagnostics were computed from a state whose density was drifting
by ~1%/h. #25 cannot be answered until item 30 is fixed, and prior FSS numbers
from this deck should be treated as unscored, not as a baseline to re-use.

### 30d. Metric excluded; band-density overwrite is a PARTIAL carrier (~58%)

Grid discriminator (800 steps, same window, baseline +2.374 %/h):

| config | drift | verdict |
|--------|-------|---------|
| `grid_stretching_ratio=1.0` (uniform dz) | +2.371 %/h | stretching excluded |
| `terrain_type=None` | +2.332 %/h | terrain excluded |

The drift survives an unstretched grid AND terrain removal, so `detJ` / `ax,ay,az`
inconsistency is NOT the mechanism. It is the update, not the metric. This also
retires the analogy to the estTimeStep and Omega-vs-rho*w bugs -- both were metric
classes, this is not.

**Blind spot in 30a's instrument.** `erf_bdy_mass_flux_diag` integrates advective
flux `rho*u*A` through the domain faces. The relaxation zone does not add mass by
flux -- it **overwrites rho directly in the band cells**. A direct overwrite is
invisible to a flux integral by construction. So "boundary flux carries 0.26%" and
"the boundary is a source" were never in conflict; the wrong channel was measured.

**Measured:** `hindcast_blend_band_density = true` -> **+1.003 %/h**, vs +2.374
baseline. **The band density overwrite carries ~58% of the drift.** Band is ~29% of
cells; acoustic equilibration crosses 576 km in ~28 min, so a band overwrite
spreads uniformly within the hour -- which is why the signature looked volumetric.

**NOT THE WHOLE STORY: +1.003 %/h survives with the blend on.** A second carrier
of comparable size remains unidentified. The cell-local mass budget is still the
next instrument, now with the band source separable as a known term.

**Do not simply ship `blend_band_density=true`.** Item #12 holds it false
deliberately: it perturbs momenta implicitly via `rho_old/rho_new`, and EVERY
scored run including the baseline was made with it false. Turning it on trades a
mass bug for a momentum bug unless #12's objection is addressed first, and it
would invalidate the existing baseline comparison on top of 30c.

### 30e. LOCALISED: mass is created inside `advance_dycore`. Fill and microphysics are exactly zero

Mass ledger (`erf_mass_ledger`, gated by `ERF_MASS_LEDGER`) sums rho*dV at tagged
points in a step; deltas attribute mass change to a code region. Run on the FLAT
UNSTRETCHED grid (`grid_stretching_ratio=1.0`, `terrain_type=None`) where the
metric is already excluded and dV is uniform, 60 steps, `sum_interval=1`.

| transition | region | delta mass |
|------------|--------|-----------|
| `Z->D` | FillPatchCrseLevel / boundary + relaxation FILL | **0.000e0** |
| `D->A` | old/new swap | -2.68e8 (buffer relabel) |
| `A->B` | **`advance_dycore`** | **+2.68e8**, then **+9.40e8** |
| `B->C` | `advance_microphysics` | **0.000e0** |

**Instrument validated before use:** `B->C` is exactly zero, independently
matching the Morrison toggle null; `D->A` is exactly minus the prior `A->B`,
which is `std::swap(vars_old,vars_new)` relabelling buffers; per-step relative
gain 1.68e-7 -> 5.9e-7 as dt ramps from 0.2 s, extrapolating to ~1e-6 at
dt~1.6 s as predicted.

**The direct-write channel at FillPatch level is CLEAN** -- zero to twelve digits.
That does not contradict 30d: `realbdy_compute_interior_ghost_rhs` runs as part of
the SLOW RHS *inside* `advance_dycore`, so the band relaxation source and the
advection/acoustic update are both inside the single `A->B` interval. The band's
58% and the remaining 42% are both in there.

**Next split (not yet done):** instrument inside `advance_dycore` -- slow RHS
(with the relaxation source separable) vs the acoustic substep loop. The substeps
advance only rho and rho*theta, so a residual there indicts the acoustic solve;
one in the slow path indicts advection or the relaxation source.

### 30f. RETRACTION of 30e, and the ledger's precision floor

**30e is WITHDRAWN as a firm result.** Both ledgers accumulated in `amrex::Real`
= float32. A float32 sum of ~1.6e15 has a **1.34e8 ULP**; the per-step mass change
is ~2.7e8, i.e. **two ULP**. The MRI stage sums were worse: ~4.5e5 with a 0.0625
ULP against a ~0.075 change, so every stage delta printed was 0.5-1 ULP of
quantisation noise -- and the two ledgers disagreed in SIGN, which is what two
independent noise floors look like.

So 30e's "boundary/relaxation fill and microphysics are exactly zero" was
RESOLUTION, not measurement: zero to twelve printed digits was zero to ~2 ULP, and
any contribution below ~1e-8 relative was invisible. The localisation of the gain
to `advance_dycore` rested on the same 2-ULP delta and must also be re-measured.

**What survives:** every AGGREGATE comparison. The toggle matrix (wall-flux,
sponges, w_damping, Morrison, theta blend, blend_band_density, stretching,
terrain), the cfl dt-scaling test, and the +0.73%/800-step and +20%/18h drifts are
differences of ~1e-2 relative -- seven orders above the noise floor. Those stand.

Both ledgers now force `ReduceData<double>` regardless of build precision.

**Caution recorded:** a `cmake` exit of 0 confirms nothing about whether an
intended edit landed. A patch script threw `substring not found`, no edit was
applied, and the rebuild still reported `exit=0 errors=0` because it compiled
unchanged source. Verify the EDIT (grep the new symbol), not the build status.

### 30g. Stage ledger still does not close -- buffer mismatch, attribution BLOCKED

With double accumulation:

| quantity | value | relative |
|----------|-------|----------|
| `A->B` (whole `advance_dycore`) | +3.523e8 | +2.2e-7 |
| MRI stages `R0(nrk=0) -> R3(nrk=2)` | -0.01155 | -2.6e-8 |
| MRI sum x dV at `R0` | 1.60805e15 | -- |
| `MASSLEDGER` at `A` | 1.59634e15 | -- |

**The absolute values differ by 0.73%, so the two ledgers are measuring different
arrays.** `MASSLEDGER` sums `vars_new[lev][Vars::cons]`; the MRI ledger sums
`S_new[0]`, the integrator's **IntVars** state (`state_new`, `ERF_Advance.cpp:466`)
-- separate storage that `advance_dycore` fills from vars, integrates, and copies
back. Their deltas are not comparable and the sub-tags do not sum to `A->B`.

**Attribution is blocked until this closes.** The suggestive reading -- MRI
integration LOSES mass while `advance_dycore` GAINS it, putting the source in the
vars<->IntVars conversion or the pre/post code rather than in the RK stages or
acoustic substeps -- is NOT claimed. With a 0.73% buffer mismatch it could be two
unrelated quantities being differenced.

**Next:** tag the SAME array on both sides. Either ledger `vars_new[lev][cons]`
inside `advance_dycore` around the `integrator.advance()` call, or ledger the
IntVars `state_new[IntVars::cons]` at the `A`/`B` points too. Then confirm the
sub-tags sum to `A->B` before reading anything from the split.

### 30h. VALIDATION CLOSED: mass is created inside `mri_integrator.advance()`

Validation chain completed in order, double accumulation, flat unstretched grid.

**1. Step stamps agree.** `erf_ledger_step` is stamped by `Evolve`; every tag
reports the same step. The 0.73% gap is NOT a one-step offset.

**2. Positive control PASSES EXACTLY.** Step 1 `D->A` = -352313405.0, which is
bit-exactly -(step 0 `A->B` = +352313405.0). The swap is a relabel and the ledger
reproduces it to full double precision -- so it tracks real state changes, not
merely failing to invent them. (`B->C = 0` alone could not have shown this.)

**3. Closure holds, and is sharper than expected:**

| interval | step 0 | step 1 |
|----------|--------|--------|
| `Z->D` boundary/relaxation FILL | **0** | **0** |
| `D->A` swap | 0 | -3.523e8 (relabel) |
| `A->M0` pre-integrator | **0** | **0** |
| **`M0->M1` `mri_integrator.advance()`** | **+3.523e8** | **+7.744e8** |
| `M1->B` post-integrator | **0** | -- |
| `B->C` microphysics | **0** | -- |

`A->B` equals `M0->M1` exactly; every other interval is identically zero.

**4. `Z->D` re-tested at full resolution: EXACTLY ZERO**, bit-identical at 17
significant digits. 30e's retracted claim is now CONFIRMED properly rather than by
artifact. The direct-write channel is clean; the fill adds no mass.

**5. The unclaimed 30g finding is REVERSED.** The gain is entirely INSIDE
`mri_integrator.advance()`, so it is in the RK stages / acoustic substeps, NOT in
the vars<->IntVars conversion or the pre/post code. The earlier "MRI loses mass"
reading came from the offset ledger and was correctly not asserted.

### 30i. STILL BLOCKED: the substep-vs-slow-path split

The MRI stage ledger sums `S_new[0]` (integrator IntVars state); the MASSLEDGER
sums `vars_new[level][Vars::cons]`. They remain 0.73% apart in absolute value
(MRI 451309.6998 x dV 3.56306e9 = 1.60803e15 vs MASSLEDGER 1.59634e15), and their
deltas disagree in sign: MRI net over step 0 is **-0.01155** against `M0->M1` =
**+0.09888** in the same units. **They are not the same array**, so the stage split
(`R1->R2` acoustic substeps vs `R2->R3` slow path) cannot be interpreted.

**Next:** print a double-precision sum of `state_new[IntVars::cons]` at the M0/M1
points as well. That establishes the relationship between the two arrays directly
-- if it matches MRI `R0`/`R3`, the stage deltas are real and simply describe a
different (IntVars) representation, and the conversion between them is where the
sign flips. Only then split substeps vs slow path, and difference blend-on against
blend-off to separate the known 58% band term.

### 30j. 0.73% offset EXPLAINED; sign lead is DEAD; and a question that threatens #30 itself

**The two ledgers were always the same array.** `ERF_Advance.cpp:397`:
`state_new.push_back(MultiFab(S_new, amrex::make_alias, 0, nvars))` -- `state_new
[IntVars::cons]` is an ALIAS of `vars_new[lev][Vars::cons]`. No second array, no
ghost mismatch, no different BoxArray. The gap is only that `erf_mass_ledger`
weights by `detJ` and the MRI ledger does not.

**Measured:** `meanJ = 0.99272` (= rho-weighted mean of detJ) = the whole 0.73%.

**The sign-disagreement lead (30g/30i) is DEAD.** Sigma(rho) and Sigma(rho*detJ)
are not two representations of one conserved quantity, so their disagreeing is
expected. Nothing indicts the vars<->IntVars map. Correctly never asserted.

**But the decomposition is alarming:**

```
M0: mass= 1596340559782521.8  rawsum= 451309.69977148622  meanJ= 0.99272150294
M1: mass= 1596340912095926.8  rawsum= 451309.68822192401  meanJ= 0.99272174744
```

Self-consistent: dmass/mass = drawsum/rawsum + dmeanJ/meanJ
= -2.56e-8 + 2.46e-7 = **+2.21e-7**.

**Raw Sigma(rho) DECREASES (-2.6e-8/step, essentially conserved). The entire
"+2.2e-7 mass gain" is meanJ rising** -- i.e. rho redistributing toward
larger-detJ cells, NOT rho being created.

**OPEN AND CRITICAL: which quantity is the mass?** If `Rho_comp` is physical
density and cell volume is `detJ*dV`, mass = Sigma(rho*detJ*dV) and it grows. If
`Rho_comp` is already a computational density (rho*detJ), mass = Sigma(rho)*dV --
the rawsum -- and it is NEARLY CONSERVED. I assumed the former when writing the
ledger and did NOT verify it. **If the latter is true, the ~20%/day "mass gain"
of item 30 may be a diagnostic artifact of detJ weighting rather than a defect.**
ERF's own `sum_integrated_quantities` MASS is the authority -- read what IT
weights by before anything else in #30 is trusted.

**ALSO: `detJ_null = 0` and `meanJ = 0.9927` with `terrain_type = None`.** detJ is
neither absent nor unity on runs treated as FLAT. The "metric excluded"
discriminator (30d) was measured on runs that were not unit-Jacobian and must be
revisited.

**Next, in order:** (1) read `sum_integrated_quantities` to establish ERF's mass
definition; (2) if it is detJ-weighted, confirm the conserved variable's identity
in the continuity update; (3) only then resume the substep/slow split.

### 30k. #30 SURVIVES: ERF's MASS is detJ-weighted and IS the conserved quantity

**(1) ERF's own definition.** `ERF::volWgtSumMF` (`Source/Utils/ERF_VolWgtSum.cpp:20`):

```cpp
if (SolverChoice::mesh_type == MeshType::ConstantDz) {
    dst_arr(i,j,k,0) = src_arr(i,j,k,comp) / (mfx_arr*mfy_arr);
} else {
    dst_arr(i,j,k,0) = src_arr(i,j,k,comp) * dJ_arr(i,j,k) / (mfx_arr*mfy_arr);
}
```

MASS **is** detJ-weighted, and also divided by m^2 -- the source comment states
"The quantity that is conserved is not (rho S), but rather (rho S / m^2)".

**(2) Rho_comp is PHYSICAL density.** From the continuity form
`rho_t = -(m^2/detJ)*[Dx(ax*rho u/m_u)/dx + ...]`, multiply by `detJ/m^2` (both
static): `d/dt[rho*detJ/m^2] = -div(flux)`. Summed over cells the divergence
telescopes to boundary fluxes, so **Sigma(rho*detJ/m^2) is the conserved
quantity** -- precisely what volWgtSumMF computes. detJ is the volume factor; it
is NOT already absorbed into Rho_comp.

**Therefore #30 DOES NOT DISSOLVE.** The rawsum being nearly conserved
(-2.6e-8/step) is incidental -- it is not the conserved quantity. The conserved
quantity grows at +2.2e-7/step, and pure flux form should hold it to round-off.
The defect is real and remains localised to `mri_integrator.advance()` (30h).

My ledger omitted the 1/m^2 factor. Map factors are static, so this cannot create
drift, but the ledger should carry it for exact agreement with ERF's MASS.

**(3) The detJ anomaly is explained, and 30d is INVALID.** The `ConstantDz` branch
is selected by **`mesh_type`**, not by `terrain_type`. detJ is populated from the
VERTICAL MESH independently of terrain, so `erf.terrain_type = None` never made
the grid unit-Jacobian (measured: `detJ_null=0`, `meanJ=0.9927`). **30d's "metric
excluded" toggled the wrong knob** and must be re-run against `mesh_type` /
`MeshType::ConstantDz`, not `terrain_type`.

**Next:** (a) re-run the metric discriminator on a genuine `MeshType::ConstantDz`
configuration; (b) resume the substep-vs-slow-path split inside
`mri_integrator.advance()` using a detJ- and m^2-weighted stage ledger, with
blend-on/blend-off differencing for the known 58% band term.

### 30l. The metric discriminator was never valid, and a genuine ConstantDz run is blocked

**Why 30d was wrong.** `ERF_DataStruct.H:370-376`:

```cpp
if (grid_stretching_ratio >= 1) {
    if (terrain_type == TerrainType::None) {
        terrain_type = TerrainType::StaticFittedMesh;   // silently re-enables terrain
    }
    if (mesh_type == MeshType::ConstantDz) {
        mesh_type = MeshType::StretchedDz;
    }
}
```

Setting `grid_stretching_ratio = 1.0` to mean "no stretching" satisfies `>= 1`, so
it **re-enabled terrain and forced StretchedDz**, overriding the
`terrain_type = None` set in the same deck. The unstretched value is **0**, not 1.
Both the "nostretch" and "flat" runs in 30d were therefore neither, which is why
`detJ_null=0` and `meanJ=0.9927`. **30d's "metric excluded" is withdrawn.**

**A correct ConstantDz run is blocked.** With `grid_stretching_ratio = 0` and
`terrain_type = None`, uniform `dz = 395.9 m` puts the first cell centre at
197.95 m and MOST aborts:

```
Assertion `zref_tmp >= m_zlo + myhalf * m_dz' failed
  ERF_MOSTAverage.cpp:430  Msg: Query point must be past first z-cell!
```

The stretched mesh exists to give 25 m near the surface, which is what MOST
requires; removing stretching removes it. The two are incompatible at 48 levels
over 19 km. Options for whoever resumes: raise `erf.most.zref` above ~198 m (
changes surface physics -- acceptable for a conservation discriminator if noted),
or raise nz enough to keep dz small (expensive), or disable the surface layer.

**The metric is now the LEADING candidate, not an excluded one.** Sigma(rho*detJ)
growing while Sigma(rho) is flat to 2.6e-8 is the signature of an inconsistent
Jacobian between where flux is constructed and where divergence is taken -- the
same class as the estTimeStep bug and the Omega-vs-rho*w scalar advection bug,
both of which were exactly zero on uniform grids and nonzero on this one.

**Targeted read available without any run:** compare the detJ and a-factors used
in `AdvectionSrcForRho`'s FLUX construction against the detJ used in the
DIVERGENCE and in the update. A factor applied in one place and not the other, or
applied at cell centres where the flux needs faces, produces exactly this
signature.

### 30m. METRIC EXCLUDED (for real this time): drift survives a genuine unit-Jacobian grid, the interior algebra telescopes in M, and no volumetric rho source is active

Three parallel probes (2026-07-28), all on the production deck in `run_a3`
(192x96x48, cfl 0.2). Production dispatch is `StaticFittedMesh` ->
**`VariableDz`** (`ERF_DataStruct.H:352-355` forces it; the deck's explicit
`terrain_type = StaticFittedMesh` means the StretchedDz variants never run) ->
`erf_substep_T`.

**A. Jacobian consistency audit (6-lens code trace, adversarial verify pass
lost to session limits -- key claims re-verified by hand below).** For the
VariableDz path, BOTH rho updates telescope EXACTLY in ERF's own
M = Sigma rho*detJ/(mfx*mfy):

- slow: `advectionSrc = -mfsq/detJ * div(ax*rho_u/mf_uy, ay*rho_v/mf_vx,
  az*Omega/mfsq)` (`ERF_AdvectionSrcForState.cpp:85-88`), applied unweighted;
  weighting a cell by `detJ/mfsq` cancels the prefactor and leaves shared face
  fluxes.
- fast: `Substep_T` horizontal perturbation fluxes carry the identical
  h_zeta/mf weighting (`ERF_Substep_T.cpp:467-476`), the update divides by detJ
  (`:708`); same cancellation. Vertical telescopes within columns (Omega
  hard-zeroed at k=klo, zero at the SlipWall lid).

So the interior advection CANNOT change M: every unit of M change is either a
domain-wall face flux or a volumetric source. Raw Sigma(rho) is NOT an
invariant of the scheme at all -- its flatness on the production grid (30k) was
flow coincidence, not algebra.

**Volumetric rho sources: none active.** `cc_src` zeroes Rho
(`ERF_MakeSources.cpp:69`), nothing writes it back; band rho relaxation is
gated by `erf.hindcast_mass_consistent_bdy` (default FALSE, not in the deck,
and rejected in item 25's attempt (a) on w_rms + characteristic grounds);
`fill_from_realbdy` leaves rho zero-gradient (ghost only, not valid cells).

**C. Unit-Jacobian control RAN (50 steps, clean).** The MOST block in 30l was
cleared with `erf.most.zref=250` (first cell centre 197.95 m < 250 m < lid).
`meanW = 1` to 17 digits: detJ AND map-factor weights are exactly unity, so M
and raw Sigma(rho) are the same quantity. **The drift SURVIVES: ~+9.9e-5 over
50 steps, ~2.9e-6/step late (dt ramping toward 2.5 s). The metric is OUT.**
(30l's "leading candidate" framing is withdrawn in turn. The signature
argument was wrong because raw-Sigma flatness was never protected.)

**B. detJ/m^2-weighted stage ledger (shared reduction with MASSLEDGER,
`erf_weighted_mass_sum`).** Closure is bit-exact: R3@nrk2 == M1 == B == C ==
Z(next) == D(next) == R0(next), every step. The step>=1 M0/A anomalies are the
known buffer-cycle artifact (vars_new holds stale data between the swap and the
MRI's internal init) -- budget on R0->R3 and the Z/D/B/C chain only. 100% of
the M gain sits in the R1->R2 windows of the three RK stages. Crucially the
substep loop applies `dtau*(slow_rhs - fast_div)` per substep
(`ERF_Substep_T.cpp:710`), so the R1->R2 window contains the SLOW RHS
application too -- substep-locality never distinguished acoustic fluxes from
slow sources. Blend on/off (`hindcast_blend_band_density`): per-stage deltas
differ <0.3% at 6 steps -- not the carrier of the substep-window gain.

**Consequence.** M gain == net discretized domain-wall flux, integrated over
the substeps. The 30-era "boundary flux excluded (0.26%)" null came from
`erf_bdy_mass_flux_diag` -- a ONCE-PER-STEP SNAPSHOT of post-step face fluxes,
not the substep-integrated fluxes the update actually applied. That null is
withdrawn as instrument-weak. The 58%-band/42%-interior split from the per-cell
budget instrument is likewise suspect (its flux reconstruction is not the
model's own).

**Next split (30i tap, zero derivation risk):** the slow rho RHS is constant
within a stage, so its M contribution is exactly `stagedt * W(F_slow[cons])`,
logged at R1. fast-wall = measured (R1->R2) - slow part. This names which path
carries the wall imbalance: the WENOZ5 slow fluxes or the acoustic substep
fluxes.

### 30n. MECHANISM NAMED: net wall inflow in the SLOW rho advection -- a 3.6% in/out asymmetry of the x-throughflow

Two taps, both validated against a known-nonzero control before reading.

**30i (slow-vs-fast split).** The slow rho RHS is constant across a stage's
substeps, so its M contribution is exactly `stagedt * W(F_slow[cons])`. Logged
at R1 with the shared weighted reduction. Per stage, per step (production deck,
6 steps): `stagedt*W_rhs` reproduces the measured stage M gain (referenced to
S_old -- each RK stage restarts from S_old) to 91-99%; the acoustic remainder
is +1-9% and shrinking. **The slow WENOZ5 rho advection carries the drift; the
acoustic substeps are a minor positive correction.**

**30j (per-wall decomposition).** The slow wall fluxes summed per wall
(`ax*xmom/mf_uy` at the x-faces, `ay*ymom/mf_vx` at the y-faces, in double):
`net = xlo - xhi + ylo - yhi` equals `W_rhs` to ~1e-7 RELATIVE at every stage
of every step. The four lateral walls are the ENTIRE slow rho RHS integral --
no lid leak, no interior residual. Numbers (step 0, M-units/s):

| wall | flux | meaning |
|---|---|---|
| xlo | +6.05e10 | west wall, inflow (AR westerlies) |
| xhi | +5.83e10 | east wall, outflow |
| ylo | -3.4e8  | south wall, net outflow |
| yhi | +0.8e8  | north wall, net outflow |
| net | +1.7e9  | = W_rhs exactly |

The drift is a **3.6% imbalance between gross x-inflow and gross x-outflow**
(net = 1/37th of gross). Real synoptic mass convergence is an order of
magnitude smaller than this. The rate grew 1.7e9 -> 3.0e9/s over the first six
steps (cold-start adjustment transient); settled rate measured separately.

**Attribution boundary (open, non-blocking).** Which side of the interface owns
the 3.6% -- the interpolated driver data as discretized on these faces
(erftools ~300 m offset, terrain-following areas, map factors) or the model
state's deviation from the driver at the walls -- is NOT settled. It does not
block the fix: both validated remedies (below) drive M toward the ERA5 target
through the wall flux regardless of which side the imbalance enters from.

**Fix candidates, both already implemented and validated in this fork, neither
in the production deck:**
- `erf.hindcast_wall_flux_correction=true` (+`_tau=3600`): per-wall-column
  barotropic dvel, flux-form, dynamically invisible, stable at tau=3600.
  Historically halved drift (item 25 attempt b; float32-era numbers).
- NSCBC combo `nscbc_lateral=1 nscbc_outflow=0 nscbc_parts=31 nscbc_mass_tau=60`:
  validated on a 24-h Jan-9 wet run -- settled drift -0.098%/day, band-w
  artifact 15.5% -> 0.5%, spurious wall precipitation eliminated, clean exit.
  Caveat: interior d>=20 precip 7.18x vs Davies 4.29x on that comparison.

Settled-drift A/B on the current production deck (2000 steps each: base / wfc /
nscbc, double ledger) decides which ships for the gate.

### 30o. Gate take 2 (NSCBC): mass BOUNDED for 4.7 h, then a dt-collapse FPE -- the NSCBC 24-h validation does not transfer to this deck

Gate take 1 was invalid (silent deck-copy failure: run_a3 is root-owned from
Docker, the host-side `cp` hit Permission denied, and the script did not check
it -- the run executed plain Davies and climbed +2.1% by 2 h. Lesson repeated:
verify the ARTIFACT, not the exit code. Take 2 verified the deck by grep
before launch.)

**Take 2 (deck-verified NSCBC combo):**
- Mass: bounded the entire run. +0.011% at t=1241 s (bit-matching the A/B
  signature), +0.0033% at 2 h, -0.0054% at 4 h, -0.007% at death. The
  constraint works.
- Death: dt collapsed 1.57 -> 0.77 (t=11.7 ks) -> 0.79 (16.0 ks) -> 0.19 s
  (17.07 ks) and the run died with "Erroneous arithmetic operation" (FPE) in
  the fast integration at model t=17073 s (4.74 h), step 16630. NO warning
  storm preceded it (zero w-damping / low-T / negative-theta lines). Evidence:
  run_a3/gate_nscbc_fail/ (log, Backtrace.0, plt at 0/3534/6143/8824/12459 s).
- The 24-h NSCBC validation (item 25 era) was on the 128x64x32 deck. On the
  production 192x96x48 deck the scheme is unstable by ~4.7 h into Jan-9.
  Davies on THIS deck ran 18 h before dying of the (now-fixed) 29m UAF -- the
  4.7 h collapse is NSCBC-specific or NSCBC-x-storm, not the deck itself.

**Decision: decouple the global mass constraint from the BC scheme.** The
nscbc_mass_tau block is scheme-agnostic (measures M vs the ERA5 target, adds a
clamped uniform additive du on outflow faces, runs last). New knob
`erf.hindcast_global_mass_tau` enables the same block under Davies; the NSCBC
knobs and semantics are untouched. Deck switched to Davies +
`hindcast_global_mass_tau = 60`. Validation: 2000-step A/B expecting the
bounded signature with Davies' dt; then gate take 3.

Open (parked): why NSCBC collapses on this deck; the 3.6% attribution
(driver-data-side vs model-side) from 30n.

### 30p. Gate take 3 (Davies + global controller): MASS BOUNDED FOR 18 H; the run still dies at 18.0003 h -- and the death is NOT mass, NOT data, NOT deterministic from a checkpoint

**The mass fix holds under storm conditions.** Davies + `hindcast_global_mass_tau=60`:
M stayed within [-0.018%, +0.012%] of start for the full 18 h (unfixed: +9.1%/day),
tracking the ERA5 target through the Jan-9 AR. dt healthy throughout (1.2-1.6 s).
2000-step validation had already shown peak transient +0.0124% (better than NSCBC's
+0.78% -- the controller acts from step 1). The MASS GATE CONDITION IS MET.

**The completion gate is not.** Death at TIME = 64801.05 s = 18.0003 h, step 46633,
"Erroneous arithmetic operation" (SIGFPE; deck has amrex.fpe_trap_invalid=1),
immediately after "Reading weather data 64801 6 7 9" and "Reading surface data
64801 6 7 9" -- the FIRST substep after the boundary window rotated to frames
[6,7] (18:00 -> 21:00 UTC). This is the same model time as the ORIGINAL 18h death
(pre-mass-fix), which the mass hypothesis was supposed to explain. **Mass is now
bounded and the run still dies there: the mass explanation for the 18 h death is
FALSIFIED. The frame rotation itself is the trigger.**

Facts assembled (all reads/cheap runs):
- Backtrace (addr2line vs the live binary): SIGFPE handler ->
  `ERF::fill_from_realbdy` <- `FillIntermediatePatch` <- substep apply_bcs
  lambda <- advance_dycore. Host-side invalid FP op (device kernels cannot
  raise host SIGFPE).
- All 18 frame files (3-D + surface) scanned: zero NaN/inf, identical sizes;
  3-D frames have identical raw min/max. Frame 7 is not visibly poisoned.
- **Restart from chk46632 (the death step) does NOT reproduce**: 68 steps
  through the crossing, clean, with the controller ON and OFF (both). Same
  time, same frames, same binary -- different outcome. The trap is therefore
  HISTORY-DEPENDENT (fresh run carries 46k steps of allocation/rotation
  history; restart carries none) -- the #29 memory-class signature, not an
  arithmetic property of the data.
- gate take 2 (NSCBC) also died by FPE but at 17073 s mid-fast-integration
  (not a frame boundary); its backtrace addresses are unresolvable (binary
  rebuilt since). Same underlying class or not: UNKNOWN.
- Reproduction attempt from chk29146 (12 h, crosses the 15 h and 18 h
  rotations with 17k steps of history): launched, result pending.

Segment-relevant property: a restart that brackets the 18 h rotation clears
it. Ugly but load-bearing for production sequencing if the root cause holds
out.

STOPPED here per standing rule: this is a new mechanism (memory-history class,
#29 family) and the hunt needs sign-off. Evidence: run_a3/gate_nscbc_fail/
(take 2), scratchpad gate24c.full / repro18_*.full, chk13785/29146/46632, plt
hourly through 17 h.

30p addendum: the chk29146 (12 h) restart ALSO survived -- 17,554 steps through
BOTH the 15 h and 18 h rotations, clean to t=64983 s. 17k steps of history + two
rotations is insufficient. The trap needs the full cold-init history, or is a
rare race (consistent with item 13 non-reproducibility and with take 2's FPE at
an unrelated time under NSCBC). Gate relaunched as-is: a second death at exactly
18.0003 h would argue cold-init x rotation-7 interaction over a race; a pass is
a legitimate continuous 24 h. Take-3 evidence moved to run_a3/gate3_evidence/.

30p final: gate take 4 (identical config, independent fresh run) died at
TIME = 64801.07 s, step 46617 -- within 0.02 model-seconds of take 3's death
(64801.05, step 46632), with MASS values matching take 3 digit-for-digit at
every 4 h checkpoint. Same resolved FPE site (fill_from_realbdy <-
FillIntermediatePatch <- substep apply_bcs). The surviving 12 h restart
executed the IDENTICAL "Reading weather/surface data 6 7 9" rotation
(ERF_WeatherDataInterpolation.cpp:1492 / ERF_SurfaceDataInterpolation.cpp:293)
at t=64800.18 and ran on cleanly.

**VERDICT: deterministic cold-init x rotation-7 interaction. Not a race, not
mass, not data.** Two fresh runs die at the same model event; two restarts
(12 h and death-step) cross it cleanly. The discriminating variable is the
INIT PATH: some object touched by the [6,7] rotation is initialized/sized
differently at cold init than on restart -- the mirror image of 29m (which was
restart-path-only). Hunt lead for whoever picks this up: enumerate frame-
indexed allocations in ERF_WeatherDataInterpolation.cpp /
ERF_SurfaceDataInterpolation.cpp whose init-path sizing differs from the
restart path, and what the 7th entry first touches. HELD pending sign-off (new
mechanism rule).

Production impact tonight: 0-18 h runs continuously with bounded mass; a
restart bracketing the 18 h rotation completes the remainder. The mass fix
(hindcast_global_mass_tau=60) is independent of this bug and stands.

30p coda: the 18->24 h completion from chk46617 ran to TIME = 86400 exactly,
exit 0, zero FPE, crossing the 21 h rotation cleanly on the restart path as
predicted. Full Jan-9 day mass excursion, 24 storm hours: -0.032% (unfixed:
~+9%/day). Day artifacts in run_a3: hourly plots 0-24 h, chk at 6/12/18 h +
final. The restart-bracketed day is production-viable pending the 30p hunt.

### 25 SCORED 2026-07-28: first mass-bounded, verified-IC day. The interior is precipitation-STARVED; the band eats the storm at the walls

**Every prior FSS/precip number in this file (item 25 tables, 24a, the NSCBC
wet-test comparison) was measured on runs inflating ~9-20%/day and is hereby
UNSCORED -- historical context only.**

Scored: run_a3 Jan-9 2023 day (0-18 h continuous + 18-24 h restart leg,
rain_accum continuous through the checkpoint), Davies + global_mass_tau=60,
mass -0.032%/24 h. MRMS Pass2 24-h QPE regridded nearest onto the ERF grid
(18432/18432 valid cells); ERA5 tp summed from 24 hourly fields, bilinear,
scored IDENTICALLY. Controls: self-FSS == 1 asserted; ERA5 domain mean
13.45 mm vs the documented 13.46 grib-footprint value. FSS windows are odd
cell counts: nominal 3/9/15/30/60 km = actual 3/9/15/33/63 km.

**FSS vs MRMS (d>=20 interior, the science cells):**

| thr | scale | ERF | ERA5 (driver) | useful |
|---|---|---|---|---|
| 1 mm | 3 km  | +0.608 | +0.737 | 0.814 |
| 1 mm | 63 km | +0.685 | +0.796 | 0.814 |
| 5 mm | 3 km  | +0.403 | +0.777 | 0.737 |
| 5 mm | 33 km | +0.469 | +0.837 | 0.737 |
| 5 mm | 63 km | +0.527 | +0.874 | 0.737 |

ERA5 exceeds the useful threshold at 5 mm from 9 km up; ERF reaches it at NO
threshold and NO scale. Full-domain (d>=0) numbers flatter the model (1 mm ERF
+0.738 vs ERA5 +0.729) only because the band's spurious rain inflates 1-mm
coverage -- band-contaminated, do not score.

**The spatial story is one artifact + its shadow:**
- d=0 mean 210.4 mm/day (MRMS 22.6), max 2404.7 mm at d=0 on the 1471 m east
  wall (MRMS there: 1.8 mm). Decays inward through d=15.
- Interior d>=20: ERF 3.24 mm vs MRMS 7.24 (bias 0.45x DRY); p90 9.7 vs 19.7.
- South interior third: ERF 0.14 mm vs MRMS 6.35 -- precipitation-dead. N/S
  ratio ERF 66.6 vs MRMS 1.72 (ERA5 4.70).
- Channel Islands: ERF 0.9-5.1 mm where MRMS has 10-28 mm (3-8x dry).
- Interior power spectrum, 8-19 km: ERF/MRMS = 0.024 (ERA5/MRMS = 0.604). The
  interior field is near-featureless at storm scales.
- Categorical, interior 5 mm: ERF POD/FAR/CSI 0.274/0.239/0.253 vs ERA5
  0.821/0.261/0.636.

**Reading:** the Davies band wrings the inflowing moisture out at the walls
(210 mm/day at d=0 vs 22.6 observed) and the interior downstream starves --
dry, southern half dead, no fine-scale variance. This is the record's earlier
hypothesis ("the band's spurious ascent wrings moisture out at the inflow
edge") now QUANTIFIED on the production deck with bounded mass and a verified
IC. Mass control did not fix the hydrologic artifact: the drift and the band
rain were two symptoms of the same wall treatment, and the controller bounds
the first only. The skill blocker is the band, not the physics options --
consistent with every physics toggle scoring null.

Artifacts: /home/atyagi/ERF/scoring_jan9/*.npy (fields, spectra, N-S
profiles), scripts in the session scratchpad (score_foundation.py,
score_metrics.py, score_islands.py -- foundation gates: self-FSS, ERA5-mean
control, last-finite-rain_accum plotfile rule).

30p segment gate (2026-07-28): 48-h bracketed run PASSES. Leg 1 fresh 0-18 h
clean (stop before the fatal rotation, chk46614); leg 2 restart 18-48 h ran to
TIME=172800 exactly, exit 0, zero FPE, TWELVE weather-frame rotations crossed
on the restart path, mass +0.077% over the 30 h. (Data prep exposed a separate
pipeline defect, fixed+committed: process_pool.py derived the LCC mapping from
the DOWNLOAD area; the widened area shifted lon_0 0.875 deg and put frames
55 km off baseline -> ERF init abort. Now pinned like the stock driver;
Jan-9 control bit-identical.) The cold-init reading holds at scale: ONE
bracketing restart per segment is a legitimate production pattern pending the
root-cause fix.

### Performance pins re-tested (2026-07-28, post-#25 scoring, user-directed)

- **rad_freq_in_time: CONFIRMED still broken.** Measured from the 48-h leg-2
  log: 423 "Radiation advancing" calls over 108031 model-s = **255.4 s
  effective** against the 180 s request -- item 26's float32-epoch
  quantisation is still live at ERF_Radiation.cpp:209-211 (m_time is absolute
  epoch). Fix is the item-28 pattern (model-relative clock, run-constant
  offset formed once in double). NOT yet applied.
- **max_gpu_streams: the pin CANNOT be lifted, and the reason rewrites the
  story.** amrex.max_gpu_streams=4 on the current binary SEGFAULTS AT STEP 0,
  resolved stack identical to the 30p site: fill_from_realbdy <-
  FillIntermediatePatch <- substep apply_bcs. The pin was never a 1.21x
  luxury; it suppresses a live unsynchronised access in the boundary-fill
  path. UNIFYING HYPOTHESIS for 30p: at 1 stream the only remaining
  host/device overlap window is the frame ROTATION (host rebuilds the pinned
  bdy/forecast buffers while asynchronously-launched substep kernels still
  read them) -- which would explain the deterministic 18.0003 h cold-init
  death, the restart survival (different launch timing), the deck's "one
  event in ~390 cold-start steps at 4 streams", and item 13's
  non-reproducibility. 4-stream step-0 death = a 2-second reproducer;
  compute-sanitizer run launched.
- **cfl probes (4000 steps each, cold start through t~2.5-2.8 model-h):**
  cfl=0.3 clean (dt -> 2.43 s, +61% model-time/step vs 0.2); cfl=0.5 clean
  (dt rides the max_dt=2.5 cap, +77%). Wall/step within 5% of baseline. The
  old cfl pin premises are gone with the fixed IC, but the storm-peak window
  (12-21 h) is untested: full bracketed-day validation required before
  production adoption.

### 30q. The 4-stream reproducer is a NaN BLOW-UP, not a memory fault -- and the `fill_from_realbdy` site is a DETECTOR, not a source

compute-sanitizer memcheck on the 2-second 4-stream reproducer
(`amrex.max_gpu_streams=4 max_step=2 amrex.fpe_trap_invalid=0`):

    ========= ERROR SUMMARY: 0 errors

Zero invalid accesses, zero out-of-bounds, zero misaligned. What the run
produced instead, immediately after `Done with advance_dycore at level 0` on
**step 1**:

    Component 0 ... contains NaNs/Infs   total NaN count comp 0: 776496
    Component 1 ... 776496
    ... every component through 14, all 776496

All fifteen conserved components, the whole field, on the first step. So the
4-stream failure is a **data race producing NaN** (or a systematic missing
stream sync), not the illegal-access class of #29. The perf-pin section's
"unsynchronised access in the boundary-fill path" is retracted in its
memory-fault form; the pin still cannot be lifted, but the reason is different.

**The more useful half of this result.** With the trap ON (the deck default),
that same NaN field kills the run inside `ERF::fill_from_realbdy` -- the exact
resolved stack of the 30p 18.0003 h death. The reason is now explicit:
`fill_from_realbdy` performs HOST ordered comparisons on device-reduced values
(the wall-flux correction's `if (f > 0)` accumulation and the mass controller's
`if (cap > 0)` / `min(max(du,-2),2)` clamp, ERF_BoundaryConditionsRealbdy.cpp
~947-973). An ordered compare against a quiet NaN raises FE_INVALID, and
`amrex.fpe_trap_invalid=1` turns that into SIGFPE. **The function is the first
host code to look at a reduced field value, so it is where any NaN anywhere in
the state announces itself.** It is a smoke alarm, not a fire.

This reframes 30p: the question is no longer "what does the [6,7] rotation
corrupt in `fill_from_realbdy`" but "what goes non-finite at 18.0003 h, or --
if nothing does -- which host operation traps on a value that is locally
invalid but harmless". Instrument: rerun the cold-init gate with
`amrex.fpe_trap_invalid=0` and `erf.check_for_nans_int=1`, which either names
the first non-finite component and cells (ERF.cpp:3294 reports both) or runs
the day through. Both outcomes are decisive; the 4-stream run above is the
known-nonzero control proving the instrument reports real content.

## 19d FIXED (2026-07-28): the regime switch was driven by the DRIVER's velocity, not the model's

**Reproduced on the production deck.** 19d was diagnosed on the 128x64 deck.
Gate take 2 (NSCBC, 192x96, corrected IC, died 4.7 h by dt collapse) carries
the identical signature, measured from its hourly plotfiles as the
wall-adjacent normal velocity made OUTWARD-POSITIVE over k=0..3 -- with the
Davies day (same deck, same hours, survived) as the control:

| t (h) | face | NSCBC mean u_n | NSCBC n_inflow | NSCBC max abs | Davies mean u_n | Davies n_inflow |
|---|---|---|---|---|---|---|
| 1 | xhi | +8.14 | 129/384 | 27.7 | +2.83 | **0/384** |
| 2 | xhi | +0.05 | 275/384 | 17.2 | +2.70 | **0/384** |
| 3 | xhi | **-3.85** | 355/384 | 23.3 | +2.74 | **0/384** |
| 4 | xhi | **-6.02** | 258/384 | 28.9 | +3.99 | **0/384** |
| 4 | yhi | +2.44 | 276/768 | **28.6** | +7.09 | 8/768 |

Davies literally cannot reverse -- it relaxes toward the frame, and the frame's
u_n at xhi is outflowing all day, so `n_inflow` is exactly zero at every hour
after the first. NSCBC leaves the face free, xhi reverses wholesale between
2 h and 3 h, and the run dies at 4.7 h. Same sequence as the 19c/19d 128x64
run (reversal at 6 h, 295 m/s spike at 7 h, NaN at 7h04m), one deck up.

**Root cause.** `sigma` -- the inflow/outflow weight -- was computed from
`un_w`, the velocity the DRIVER had just prescribed on the wall face. That was
chosen because it is smooth in time and will not chatter. But it is not the
model's regime: the driver's u_n at xhi/yhi stays positive all day, so sigma
stayed 0 while the model's own flow ran INWARD at 6 m/s. Sigma 0 selects the
outflow branch, and with the shipped default `nscbc_outflow=0` that branch is a
zeroth-order extrapolation from the interior. A face carrying inflow needs 4+N
conditions; extrapolation supplies **none**. The boundary was unconstrained,
and the resulting mode grew until dt collapsed.

**Fix, both halves:**
1. The regime is taken from `un_i`, the model's own normal velocity one face
   inboard, through a Schmitt trigger with a +-`nscbc_eps` deadband, latched
   per column and per face in a persistent `nscbc_regime` MultiFab (ERF member,
   4 comps = xlo/xhi/ylo/yhi, defined lazily, -1 = unset). Outside the deadband
   the regime follows the flow; inside it the previous regime is held, so it
   cannot chatter on acoustic noise. sigma is now exactly 0 or 1 -- never
   between -- so the condition COUNT is always 4+N or 1, which the old linear
   blend could not promise (its own APPROXIMATION 2 note said so).
2. `erf.nscbc_outflow` now defaults to **1**, the characteristic (Riemann)
   outflow solve -- the actual one-incoming-condition treatment -- instead of
   the raw extrapolation that 19d showed is unconstrained on a reversed face.

Files: ERF_BoundaryConditionsRealbdy.cpp (regime latch + default), ERF.H
(`nscbc_regime`). Validation pending: NSCBC must first clear 4.7 h, then a full
day scored against Davies on the same binary.

## 26 FIXED (2026-07-28): the radiation cadence now runs on a model-relative clock

Confirmed broken first, by measurement on the 48-h leg-2 log: 423 "Radiation
advancing" calls over 108031 model-s = **255.4 s effective** against the 180 s
request. The trigger differenced `m_time`, which `ERF_AdvanceRadiation.cpp` had
formed as `t_old[lev] + start_time` -- a float32 sum of a ~1.673e9 epoch value,
ULP 128 s. The difference could therefore only be a multiple of 128 s, and 180
first cleared at 256. Fifth appearance of the float32-epoch class (#26, #28,
and the two the item-28 fix already covered).

Fix is item 28's pattern: keep the absolute instant for what genuinely needs it
(the orbital date/time, and the plotfile stamp, where 128 s is immaterial), and
give the cadence a MODEL-relative clock. `IRadiation::set_model_time()` (no-op
default), called from `advance_radiation` with `t_old[lev]`; the trigger
differences that. ULP at 86400 s is 0.008 s.

Consequence for scoring: every result in this file, #25 included, was produced
with radiation firing at 255 s. The Davies baseline must be re-run on the fixed
binary before it can be compared against anything built on it.

## 25 RE-SCORED (2026-07-28, same day): the first scoring was geolocated 86 km off. Verdict survives, numbers do not

`score_foundation.py` built the ERF grid's lat/lon from the LCC of the
**download** area in era5_input.txt, `(36.0,-125.0,30.75,-115.25)` -> lon_0
-120.125. The frames -- and therefore the model grid -- are on the **pinned**
projection area `(36.0,-123.25,31.25,-115.25)` -> lon_0 -119.25. This is the
same defect as the process_pool.py projection drift found earlier the same day,
in the scoring script instead of the pipeline, and it went unnoticed because the
only geolocation gate was a domain-mean sanity band, which a rigid shift barely
moves.

**Measured, not argued.** The frame .bin carries both its projected grid
(xvec[nx], yvec[ny]) and the lat/lon of the same points, so the projection is
over-determined:

| candidate LCC | max residual vs frame coordinates |
|---|---|
| pinned  (lon_0 -119.25)  | **0.42 m** |
| scoring (lon_0 -120.125) | **88,211 m** |

Independent check on the corrected mapping: the model's own terrain islands now
land on the real Channel Islands -- 34.04/-120.38 San Miguel, 33.97/-120.11
Santa Rosa, 34.02/-119.74 Santa Cruz, 34.01/-119.43 Anacapa, 33.24/-119.50 San
Nicolas, 33.47/-119.04 Santa Barbara, 33.39/-118.44 Catalina, 32.90/-118.48 San
Clemente -- every one within 0.02 deg of its true position.

The old sampling pulled MRMS and ERA5 from 86 km west, i.e. largely offshore,
where MRMS QPE is smooth: that flattened the observed field (old MRMS interior
N/S ratio 1.72, i.e. no north-south storm gradient at all, against 9.70 truth)
and inflated every ratio taken against its variance.

**Corrected numbers.** Gates: projection residual 0.42 m; regrid-vs-native
interpolation control ERA5 15.06 vs 15.04 mm and MRMS 14.84 vs 15.02 mm;
self-FSS == 1.

| thr | scale | ERF | ERA5 | useful | (old ERF / old ERA5) |
|---|---|---|---|---|---|
| 1 mm | 3 km  | +0.646 | +0.773 | 0.792 | 0.608 / 0.737 |
| 1 mm | 63 km | +0.742 | +0.822 | 0.792 | 0.685 / 0.796 |
| 5 mm | 3 km  | +0.466 | +0.829 | 0.713 | 0.403 / 0.777 |
| 5 mm | 33 km | +0.546 | +0.890 | 0.713 | 0.469 / 0.837 |
| 5 mm | 63 km | +0.617 | +0.926 | 0.713 | 0.527 / 0.874 |

(d >= 20 interior.) ERF still clears the believable threshold at NO threshold
and NO scale; ERA5 clears it at 5 mm everywhere and at 1 mm from 30 km up.

- Interior d>=20 bias: ERF 3.24 mm vs MRMS **8.48** -> **0.38x dry** (was 0.45x).
- Whole domain: ERF 34.20 mm vs MRMS 14.84 -> 2.31x WET, all of it the band.
- Band unchanged (ERF-internal, so the projection never touched it): d=0 mean
  210.4 mm/day, max 2404.7 mm, decaying inward through d=15.
- South interior third: ERF 0.14 mm vs MRMS 2.20; N/S ratio ERF 66.6, MRMS 9.7,
  ERA5 6.6. The south is still dead.
- Islands: ERF 0.0-5.1 mm where MRMS has up to 59.7 on the three western
  islands -- **10-20x dry**, worse than the old 3-8x reading.
- Interior 8-19 km spectral power: ERF/MRMS 0.003, ERA5/MRMS 0.008. The old
  0.024 / 0.604 pair was the offshore-flattened MRMS variance in the
  denominator. Read it as ERF carrying **0.4x the small-scale variance of its
  own 25 km driver**.
- Categorical, interior 5 mm: ERF POD/FAR/CSI 0.326/0.189/0.303 vs ERA5
  0.998/0.291/0.708.

**The verdict is unchanged and sharper:** the band wrings the inflowing moisture
out at the walls and the interior starves -- drier (0.38x), southern half dead,
essentially no storm-scale variance. Re-opening NSCBC is the right move against
exactly this failure.

Guards added so this cannot recur: `score_foundation.py` now pins the projection
area, asserts the LCC against the frame's own coordinates (<10 m), and replaces
the remembered-mean gate with a regrid-vs-native interpolation control on both
observation fields. Fields in /home/atyagi/ERF/scoring_jan9_fixed/.

### 19d premise check: NSCBC does remove the band, but the interior gain at 4 h is small

Before spending a day of GPU on the A/B, the reason for re-opening NSCBC was
tested on the pair already on disk -- gate take 2 (NSCBC, died 4.7 h) against
the Davies day, same deck, same hours. Interior = d >= 20 (the #25 science
cells), band = d = 0.

| t (h) | rain d=0 NSCBC | rain d=0 Davies | rain d>=20 NSCBC | rain d>=20 Davies | CWV d>=20 NSCBC | CWV d>=20 Davies |
|---|---|---|---|---|---|---|
| 1 | 1.09 | 16.66 | 0.097 | 0.026 | 23.15 | 25.43 |
| 2 | 1.50 | 27.65 | 0.130 | 0.107 | 21.53 | 25.20 |
| 3 | 1.67 | 37.99 | 0.181 | 0.265 | 23.47 | 25.44 |
| 4 | 1.98 | **48.17** | 0.379 | 0.334 | 25.25 | 26.08 |

(rain_accum mm; CWV = column-integrated qv, kg/m2.)

**The band is a Davies artifact and NSCBC does not have it.** By 4 h Davies has
already put 48.2 mm on the d=0 ring -- more than twice MRMS's entire 24-h total
there (18.6 mm) -- while NSCBC has 2.0 mm. A factor of **24**. This is the #25
skill blocker measured directly at its source, four hours in.

**The interior claim is weaker than the record implied.** Interior rain is 3.8x
Davies at 1 h but only +13% by 4 h, and interior column water vapour is
consistently *lower* under NSCBC in this window (21.5-25.2 vs 25.2-26.1), not
higher. So "NSCBC delivers more moisture to the interior" is not established at
4 h; what is established is that it does not wring the inflow out at the wall.
Whether that converts into interior skill over a full day is exactly what the
A/B has to answer -- it is not assumed here.

### 30r. The 18 h death is a REAL NaN, and a two-step bit-identical bracket proves it is memory, not arithmetic

Cold-init gate, identical to takes 3/4 except `amrex.fpe_trap_invalid=0` and the
NaN reporter armed (`check_for_nans_int=10`):

    Coarse STEP 46650 ends. TIME = 64807.13494
    Component 0..14 of conserved variables contains NaNs/Infs
      total NaN count, every component: 884736      (= 192*96*48, the whole field)

So the trap was not spurious. The state genuinely goes non-finite at the [6,7]
rotation; with the trap on it dies at t=64801.05 because a host compare in
`fill_from_realbdy` sees the NaN first (30q), and with the trap off it runs four
more steps and the whole field is gone.

**A checkpoint landed one step before the event** (`chk46645`, t = 64800.42936,
written by the 6-hourly cadence). Restarting from it re-executes the rotation
identically -- same read `64800.42969 6 7 9`, same weights alpha1 = 0.9999602437,
alpha2 = 3.975629807e-05 -- and runs on cleanly. Comparing the two dt sequences:

| step | cold init | restart from chk46645 |
|---|---|---|
| 46646 | 1.188644409 | 1.188644409 |
| 46647 | 1.188738823 | 1.188738823 |
| 46648 | **1.307612777** | 1.188817263 |
| 46649 | 1.438374043 | 1.188907385 |
| 46650 | 1.582211494 | 1.188995838 |

**Bit-identical for two full steps after the rotation, then divergence.** The
restart's dt continues its smooth trend; the cold run's dt *rises* 10% per step,
which is the NaN itself -- `fmax` ignores NaN operands, so as cells go non-finite
the max-wave-speed reduction runs over a shrinking finite population and the
timestep grows. (The first NaN therefore appears in step 46646 or 46647, in cells
that do not hold the domain maximum; it reaches every cell by 46650.)

**What this rules out.** Two runs executing identical arithmetic on identical
data cannot diverge. The checkpointed state is identical, the frame data is
identical, the interpolation weights are identical, the binary is identical.
Therefore something read during those two steps is **not part of the checkpointed
state and differs between the runs** -- i.e. memory content, which is the #29
family, not a formulation error. Combined with 30q (memcheck: 0 errors) the class
is now pinned: a legal read of memory whose *contents* depend on allocation
history, not an out-of-bounds and not an arithmetic defect.

**Also ruled out, this round:** uninitialized FArrayBox data in the hot path.
`fab.init_snan=1` (every newly allocated FAB starts as signalling NaN) plus a
per-step NaN check runs 30 cold-start steps clean. That does not cover
`Gpu::DeviceVector`/raw arena allocations, which is where
`FillForecastStateMultiFabs` churns ~13 device vectors on every rotation.

**Cheap reproducer for whoever continues:** run_nan/chk46645 plus the cold log.
The failing window is two steps wide and the clean counterfactual is one command.
Next discriminator on the list: vary `amrex.the_arena_init_size` on a cold run --
#29b already established that this fault class tracks arena layout, and a shift
in the failure step would confirm it here.

Production is unaffected: one bracketing restart per segment, validated to 48 h.

### 19d validated, and a second blocker found and fixed: the mass controller SATURATES under NSCBC

**Survival: the fix works.** NSCBC with the regime latch and the Riemann outflow
variant ran 16000 steps to TIME = 24373.9 s = **6.77 h**, exit 0, no NaN, dt
1.571 at the end. Gate take 2, same deck and knobs but with the driver-driven
sigma and `nscbc_outflow=0`, died at TIME = 17073 with dt collapsing to 0.158.
Cleared by 2.0 h of model time.

**Mechanism confirmed at the wall.** Face-mean outward-positive u_n at xhi:

| t (h) | broken NSCBC | fixed NSCBC | Davies |
|---|---|---|---|
| 1 | +8.14 (129/384 in) | +3.42 (60/384 in) | +2.83 (0/384 in) |
| 3 | **-3.85** (355/384 in) | +2.71 (83/384 in) | +2.74 (0/384 in) |
| 5 | — (dead) | +2.72 (69/384 in) | — |
| 6.77 | — | +2.82 (61/384 in) | — |

The fixed run allows *local* reversal -- 60-83 of 384 columns inflow at any
time, which is physical -- while the face mean stays outflowing and bounded, and
it now tracks the Davies control. The broken run reversed wholesale and died.
This is the condition-count property the latch was built to give.

**Second blocker: the global mass constraint saturates under NSCBC.** The 6.77 h
probe ran mass +2.217%, i.e. ~8%/day -- unscoreable by the #30 rule. Gain sweep
over 3000 steps, matched model time:

| `nscbc_mass_tau` | mass @ 3000 steps | peak abs(dM) |
|---|---|---|
| 60 | +0.834% | 1.671% |
| 10 | +0.832% | 1.673% |
| 5  | +0.830% | 1.672% |

A 12x change in gain with no effect at three decimal places. That is a pinned
clamp, not a slow controller -- confirmed against the null: with the controller
OFF, +2.215% / peak 2.714%, so it is acting, just railed. The clamp was the
hardcoded `du = min(max(du,-2),2)`. Under Davies the demanded du is ~0.6 m/s and
never binds; under a characteristic boundary the correction is partly re-derived
away by the next call's Riemann solve, so more headroom is needed.

Exposed as `erf.hindcast_mass_du_max` (default 2.0 -- the long-validated Davies
behaviour is bit-unchanged). Sweep:

| `hindcast_mass_du_max` | mass @ 3000 steps | peak abs(dM) | min dt |
|---|---|---|---|
| off (no controller) | +2.215% | 2.714% | — |
| 2 (default) | +0.834% | 1.671% | 1.142 |
| **8** | **+0.012%** | **0.048%** | 1.116 |
| 25 | +0.012% | 0.020% | 1.113 |

8 and 25 agree to three decimals, so the true demanded du is under 8 and the
clamp no longer binds there; dt is unaffected. NSCBC now holds mass to Davies
quality (Davies: +-0.02% over 18 h). Production NSCBC setting:
`erf.nscbc_mass_tau=60 erf.hindcast_mass_du_max=8`.

## NSCBC vs Davies, scored on the corrected deck (2026-07-28)

Both arms: same binary (post-#26 rad clock, post-#19d regime latch), same
bracketing recipe (cold 0-18 h, restart 18-24 h), Jan-9 2023, cfl 0.2.
Davies = `hindcast_global_mass_tau=60`. NSCBC = `nscbc_lateral=1
nscbc_outflow=1 nscbc_parts=31 nscbc_mass_tau=60 hindcast_mass_du_max=8`.
Both reached TIME = 86400 exactly, exit 0, zero FPE, zero NaN.
**Mass: Davies -0.032%, NSCBC -0.030%** -- both scoreable, for the first time.
Control: the Davies arm reproduces the re-scored #25 baseline to within 0.001
FSS on a different binary, so the rad-clock fix is precipitation-neutral and the
pipeline is reproducible.

**Interior FSS (d >= 20), the science cells:**

| thr | scale | Davies | NSCBC | ERA5 | useful |
|---|---|---|---|---|---|
| 1 mm | 3 km  | +0.645 | **+0.696** | +0.773 | 0.792 |
| 1 mm | 63 km | +0.741 | **+0.751** | +0.822 | 0.792 |
| 5 mm | 3 km  | +0.465 | **+0.595** | +0.829 | 0.713 |
| 5 mm | 33 km | +0.546 | **+0.650** | +0.890 | 0.713 |
| 5 mm | 63 km | +0.617 | **+0.677** | +0.926 | 0.713 |

NSCBC improves interior FSS at **every threshold and every scale**, most at
5 mm / 3 km (+0.130). Neither arm reaches the believable threshold, and ERA5
still beats both everywhere.

**What changed, and it is not subtle:**

| quantity | Davies | NSCBC | MRMS |
|---|---|---|---|
| band d=0 mean (mm/day) | 210.3 | 91.9 | 18.6 |
| band d=0 max | 2402.0 | 774.3 | 239.9 |
| interior d>=20 mean | 3.25 (0.38x DRY) | 22.34 (2.63x WET) | 8.48 |
| south interior third | 0.14 | 2.01 | 2.20 |
| N/S ratio | 67.1 | 29.0 | 9.70 |
| interior 8-19 km spectral power vs MRMS | 0.003 | 2.297 | 1.0 |
| islands, three western (mm) | 5.0 / 3.0 / 3.0 | 91.3 / 93.6 / 111.2 | 57.6 / 56.2 / 59.7 |
| interior 5 mm POD/FAR/CSI | 0.326/0.189/0.303 | 0.670/0.464/0.424 | -- |
| steps for the day | 62602 | 56095 | -- |

The wall-normal decay profile inverts. Davies falls monotonically inward
(210 -> 125 -> 92 -> 63 -> 40 -> 21 -> 4.4 at d = 0,1,2,5,10,15,25): the band
wrings the inflow out and the interior starves. NSCBC is 92 -> 57 -> 54 -> 61
-> 72 -> 62 -> 33: the precipitation is *in the domain* instead of on its edge.

**The hypothesis behind re-opening NSCBC is confirmed -- and it overshoots.**
The dead southern third comes back to life (0.14 -> 2.01 against 2.20 observed,
essentially exact). The featureless interior gains storm-scale structure
(spectral ratio 0.003 -> 2.297). The islands go from 10-20x DRY to ~1.7x wet.
Interior POD doubles. But the interior bias swings from 0.38x dry to **2.63x
wet**, FAR rises 0.189 -> 0.464, and the band is reduced 2.3x rather than
removed -- 91.9 mm/day at d=0 is still 5x observed.

**Read the FSS gain with that bias in mind.** A 2.6x wet field is rewarded by
FSS at a fixed low threshold simply through coverage. Three things argue the
gain is not only that: it is larger at 5 mm than at 1 mm, the spectrum and the
north-south gradient both move *toward* observations rather than merely up, and
the island pattern goes from absent to present. A percentile-matched-threshold
FSS would separate placement skill from bias cleanly and is the obvious next
measurement; it has not been made.

Verdict: NSCBC is the better lateral treatment for this configuration on every
structural measure, and it is now mass-bounded and stable for a full day. The
remaining error changes character -- from "the boundary eats the storm" to "the
interior rains too much" -- which is a tractable physics problem rather than a
boundary artifact.

## 31. A checkpoint written on the stop_datetime truncation step POISONS every restart from it

Found while validating cfl 0.3. Leg 1 stopped at `stop_datetime = 18:00:00`:

    Coarse STEP 31413 ends. TIME = 64767.19909 DT = 1.816368222
    Coarse STEP 31414 ends. TIME = 64768        DT = 0.8009066582
    Coarse STEP 31415 ends. TIME = 64768        DT = 2.980232239e-08   <- extra step
    Writing native checkpoint chk31415

ERF lands on the stop instant with a truncated step (31414), then takes ONE MORE
step of ~3e-08 s to close a float32 residual, and the end-of-run checkpoint is
written from that step. Restarting from it:

    31416 DT = 3.278255534e-08
    31417 DT = 3.606081123e-08     each exactly 1.1x the previous --
    31418 DT = 3.966689377e-08     the dt growth limiter
    31419 DT = 4.363358386e-08
    31420 DT = 4.799694153e-08  -> all 15 conserved components NaN

dt recovery is rate-limited to 1.1x per step, so from 3e-08 it needs ~180 steps
to return to O(1) -- and the run goes non-finite after five. Whatever breaks at
sub-microsecond dt (not chased here), the operational consequence is what
matters: **a segment that ends on stop_datetime can write a final checkpoint
that no restart can use.** Whether it happens depends on whether the truncated
step lands exactly on the stop instant: the cfl 0.2 Davies arm's final step hit
64768 exactly with dt = 0.917 and its checkpoint restarts fine, the cfl 0.3 arm's
did not and its checkpoint is poison.

This is directly in the path of the year-scale segment recipe, which chains
stop_datetime-bounded segments through restarts. **Rule: never restart from the
checkpoint written on a run's final step without checking its dt; use the
previous periodic checkpoint.** A cheap structural guard would be to skip the
end-of-run checkpoint when the final dt is a small fraction of the previous one.

### The NSCBC FSS gain was bias, not placement. Percentile-matched scoring REVERSES it

Fixed-threshold FSS credited NSCBC with +0.05 to +0.13 over Davies at every
threshold and scale. With each field thresholded at its OWN quantile so all
binary fields carry the SAME base rate as the observations, the sign flips:

Base rate 0.5834 (= MRMS exceedance of 1 mm), believable >= 0.792:

| scale | Davies | NSCBC | ERA5 | NSCBC - Davies |
|---|---|---|---|---|
| 3 km  | **0.822** | 0.649 | 0.938 | **-0.173** |
| 15 km | **0.849** | 0.676 | 0.962 | -0.174 |
| 60 km | **0.895** | 0.712 | 0.986 | -0.182 |

Base rate 0.4250 (= MRMS exceedance of 5 mm), believable >= 0.713:

| scale | Davies | NSCBC | ERA5 | NSCBC - Davies |
|---|---|---|---|---|
| 3 km  | **0.743** | 0.604 | 0.890 | -0.139 |
| 15 km | **0.778** | 0.638 | 0.928 | -0.140 |
| 60 km | **0.841** | 0.684 | 0.973 | -0.157 |

Matched thresholds (mm) at the 1 mm base rate: MRMS 0.96, Davies 0.23, NSCBC
3.56, ERA5 5.40 -- Davies must be cut at 0.23 mm to reach the observed coverage,
which is the dry bias stated as a threshold.

**Davies CLEARS the believable threshold on placement at every scale and both
base rates. NSCBC clears it nowhere.** The fixed-threshold result was the wet
bias earning coverage, exactly the confound it was suspected of being.

Independent confirmation on a statistic sharing no machinery with the
neighbourhood score -- with base rates matched, n_forecast == n_obs so FAR =
1 - POD and POD carries the comparison:

| base rate | Davies POD/CSI | NSCBC POD/CSI | ERA5 POD/CSI |
|---|---|---|---|
| 0.5834 | 0.822 / 0.698 | 0.649 / 0.481 | 0.938 / 0.882 |
| 0.4250 | 0.743 / 0.591 | 0.604 / 0.433 | 0.890 / 0.802 |

Instrument controls, all asserted in the script: self-FSS == 1; bias invariance
verified directly (multiplying a forecast by 2.63 leaves the score
bit-identical); displacement response on a known-nonzero case (MRMS shifted
5 cells scores 0.908 at 3 km, recovering to 0.992 at 63 km).

**Corrected verdict.** Davies places precipitation better and gets the amount
badly wrong (0.38x dry). NSCBC gets closer on amount and structure -- the dead
southern third, the missing islands, the flat spectrum all improve -- but places
it worse. Neither is preferable outright, and the earlier "NSCBC wins on every
structural measure" reading must not be carried forward as a skill claim.

### Where the NSCBC excess sits: the OUTFLOW walls, not a residual inflow band

Interior precipitation binned by distance from the inflow walls (xlo, ylo --
face-mean u_n negative there all day) and separately from the outflow walls
(xhi, yhi), as ERF/MRMS ratio so the observed field's own gradient divides out:

| distance (cells) | Davies, from INFLOW | NSCBC, from INFLOW | NSCBC, from OUTFLOW |
|---|---|---|---|
| 0 | 7.21 | **0.00** | **6.35** |
| 1-2 | 9.84 | 0.00 | 3.55 |
| 6-9 | 4.66 | 0.08 | 3.59 |
| 15-19 | 1.10 | 0.34 | 3.78 |
| 30-44 | 2.56 | 2.86 | 2.07 |
| 45-95 | 1.91 | 3.68 | 0.34 |

The Davies control behaves exactly as its known band should: 7-10x within 6 km
of the inflow walls, decaying to 1.1 by 45-57 km.

**Neither of the two hypotheses is right.** NSCBC has no residual inflow band --
0.03 mm/day against 8.25 observed at the inflow wall, ratio 0.00, the band is
genuinely gone. Nor is the over-delivery uniform: the ratio spans 0.00 to 6.35.
It is a NEW structure concentrated toward the OUTFLOW walls.

**The mass controller amplifies it but does not create it.** The controller adds
its uniform du on outflow faces ONLY, and raising `hindcast_mass_du_max` from 2
to 8 is what made NSCBC mass-bounded -- so it is the first suspect. Comparison at
matched model time (t = 6.00 h), mean rain_accum by distance from the outflow
walls, cells 0 / 1-2 / 6-9 / 10-14 / 30-95:

    du_max = 2 (saturated, mass +2.2%):  3.0  2.0  6.1   8.9  0.7   domain 2.71
    du_max = 8 (bounded,   mass -0.03%): 7.5  5.3  8.6  10.6  1.0   domain 3.67

Same shape -- both peak 10-14 cells inboard of the outflow walls -- at 1.35x the
amplitude domain-wide and 2.5x at the wall itself. So the structure belongs to
the outflow condition (the Riemann variant), and the controller's outflow-only du
inflates its wall-adjacent part. That is a cost of the du_max=8 setting that the
mass numbers alone do not show.

## 32. CONSOLIDATED: the float32-absolute-time class -- four instances, one root

Every one of these is the same mistake: an ABSOLUTE epoch instant (~1.673e9 s for
a 2023 case) held or differenced in single precision, where the float32 ULP is
**128 s**. Recorded here in one place because it has now been rediscovered four
times from four different symptoms, each time costing a hunt.

| # | where | symptom | status |
|---|---|---|---|
| 28 | `fill_from_realbdy` boundary interpolation | `time - start_bdy_time` quantised alpha to 128/10800; 675 impulsive boundary kicks per day | FIXED -- run-constant offset formed once in double |
| 26 | `Radiation::set_grids` cadence | `t_old + start_time` differenced; 180 s request fired at a measured **255.4 s** | FIXED -- `set_model_time()`, model-relative clock |
| -- | `ERF.cpp:626` evolve stop condition | `start_time + cur_time` vs `stop_time` did not change for ~100 consecutive steps | FIXED -- `stop_elapsed` formed once in double |
| **31** | **`ERF::stop_time` itself** | **the requested stop instant is silently displaced by up to +-64 s** | **LIVE** |

**The fourth is new and is the root of item 31.** `stop_time` is declared
`amrex::Real` (ERF.cpp:45) and holds `getEpochTime(stop_datetime)` -- an absolute
epoch, in float32. So the stop instant a deck asks for is not the stop instant
the run uses. Predicted and then measured:

| stop_datetime | exact elapsed | what ERF uses | error |
|---|---|---|---|
| 2023-01-09 18:00:00 | 64800 s | **64768 s** | **-32 s** |
| 2023-01-09 21:00:00 | 75600 s | 75648 s | +48 s |
| 2023-01-10 00:00:00 | 86400 s | 86400 s | 0 |

Both observed cases confirm it. Every leg-1 of every bracketed run in this
campaign stopped at TIME = 64768, which was read as "before the fatal rotation,
by design" -- it is actually 18:00:00 rounded down by float32. The full-day runs
end at exactly 86400 because that instant happens to be representable. The same
rounding is visible in the log timestamps: at elapsed 64799.24 s the run prints
`2023-01-09 17:59:28`, 31 s early, because the printer forms `start_time +
cur_time` in float32 too.

`ERF_ComputeTimestep.cpp:46-47` then still forms `stop_time - start_time` as a
float32 difference before casting to double, and the final truncated dt is
rounded to float32, which leaves `cur_time` a hair short of the target and
spawns the extra ~3e-08 s step that item 31 documents.

**Fix for the class:** `start_time` and `stop_time` must be `double`, and every
elapsed span formed from them in double. 98 references to `start_time` alone, so
it is a contained but non-trivial change -- and it is the only change that
retires the class instead of patching its next symptom. Until then, a deck's
stop_datetime is accurate only to +-64 s, which matters for segment boundaries
that must align with a frame instant.

### 31a. Guard shipped

Two layers, because the failure is silent and a remembered rule fails at
segment 7:

1. **In ERF** (`ERF_ComputeTimestep.cpp`): on the first `ComputeDt` of a
   restarted run, if the dt carried in from the checkpoint is more than 1000x
   below the CFL timestep for the restored state, abort with the measured ratio
   and the reason. The check sits where the CFL estimate is already in hand, so
   the comparison needs no notion of a "typical" dt, and it fires before any
   physics runs rather than five steps into a NaN.
2. **In the harness** (`pick_restart_chk.sh`): reads the fixed-layout checkpoint
   Header (line 7 istep, line 8 dt, line 9 time) and returns the LATEST-IN-TIME
   checkpoint whose dt is within 1e-3 of the largest dt across the run's
   checkpoints.

Validated against both controls before use. On the known-poisoned directory it
rejects exactly the two bad checkpoints and nothing else:

    chk31396  t=64800.8  dt=1.804          OK
    chk31415  t=64768    dt=2.98e-08       REJECT: 83886080x below the run maximum
    chk31416  t=64768    dt=3.28e-08       REJECT: 76260071x below the run maximum
    chk42597  t=86400    dt=1.034          OK

and on three known-healthy directories (run_ab_dav, run_ab_nsc, run_48h) it
rejects nothing. A first version of the script read dt from Header line 7 and
passed the poisoned checkpoints; the control caught it. Selection is by TIME,
not step number -- a run directory accumulates checkpoints from several legs and
after a restart the step count for a given model time differs between them
(run_ab_cfl03 holds chk31396 at t=64800.8 from leg 2 and chk31415 at t=64768
from leg 1, so step order is not time order).

Not applied, but the one-line source fix that would retire the hazard rather
than guard it: stop the evolve loop when the remaining span is a negligible
fraction of the CFL dt, so the degenerate step is never taken and no checkpoint
can be written from it.

## 33. NSCBC outflow variant unconfounded: the LATCH buys survival, the RIEMANN solve buys the outflow excess, and NEITHER variant fixes placement

Third full day, same binary, same bracketing recipe, identical knobs to the
Riemann arm except `nscbc_outflow=0`. Completed TIME = 86400, exit 0, zero FPE,
zero NaN, mass **-0.032%**.

**Survival is the latch, not the Riemann solve.** Gate take 2 was extrapolation
WITHOUT the regime latch and died at 4.7 h by dt collapse. This run is
extrapolation WITH the latch and completes the day. The two were confounded in
the first NSCBC result; they are not any more.

| | Davies | NSCBC Riemann | NSCBC extrapolation | MRMS |
|---|---|---|---|---|
| **PM-FSS 1 mm, 3 km** | **0.822** | 0.649 | **0.586** | useful 0.792 |
| **PM-FSS 5 mm, 3 km** | **0.743** | 0.604 | **0.478** | useful 0.713 |
| PM-FSS 5 mm, 60 km | 0.841 | 0.684 | 0.589 | useful 0.713 |
| matched-rate POD (1 mm) | 0.822 | 0.649 | 0.586 | -- |
| fixed-thr FSS 1 mm 3 km | 0.645 | 0.696 | 0.374 | 0.792 |
| interior bias | 0.38x dry | 2.63x wet | **0.19x dry** | 8.48 mm |
| domain mean (mm) | 34.14 | 45.48 | 12.83 | 14.84 |
| band d=0 (mm/day) | 210.3 | 91.9 | 12.5 | 18.6 |
| inflow-wall ratio | 7.21 | 0.00 | 0.08 | 1.0 |
| outflow-wall peak ratio | -- | **6.35** (at the wall) | **1.85** (10-14 cells in) | 1.0 |
| interior 8-19 km spectrum | 0.003 | 2.297 | 0.002 | 1.0 |
| N/S ratio | 67.1 | 29.0 | 4.72 | 9.70 |
| islands, three western (mm) | 5.0/3.0/3.0 | 91.3/93.6/111.2 | 3.5/2.3/3.2 | 57.6/56.2/59.7 |
| steps for the day | 62602 | 56095 | **79560** | -- |
| dt at day end (s) | 0.99 | 1.30 | **0.167** | -- |

**FAILS the pass condition.** Placement had to reach Davies' 0.822 / 0.743 with
the inflow band gone. The band IS gone -- 0.64 mm/day at the inflow wall against
8.25 observed, ratio 0.08, so band elimination belongs to NSCBC generally and not
to either outflow variant. But placement is 0.586 / 0.478: **worse than the
Riemann arm and far below Davies.** Percentile-matched FSS ranks the three
Davies > Riemann > extrapolation at every base rate and every scale, and
matched-rate POD ranks them identically.

**What the contrast does establish, cleanly:**
- The **outflow-wall excess is a Riemann artifact**. Its peak ratio drops 6.35 ->
  1.85 and moves off the wall to 10-14 cells inboard when the Riemann solve is
  replaced by extrapolation. Task 2's localisation is confirmed by substitution.
- The **outflow variant sets the interior amount**, and it swings hugely: 2.63x
  WET under Riemann, 0.19x DRY under extrapolation -- drier even than Davies.
  Both are wrong, from opposite sides, with the same latch and the same mass
  controller.
- Extrapolation costs **27% more steps than Davies and 42% more than Riemann**,
  with dt degrading to 0.167 s by day end against ~1.0-1.3 s for the others. It
  completes, but it is a marginal boundary, not a comfortable one.

**Verdict: neither NSCBC outflow variant is usable as-is.** Davies remains the
best placement by a wide margin and the only arm clearing the believable
threshold; its defect is amplitude (0.38x dry) and a 210 mm/day wall band. The
untried lever is the Riemann sigma relaxation coefficient, which has never been
tuned -- the Riemann arm is the better of the two NSCBC starting points on both
placement and cost.

## 34. NSCBC CLOSED: sigma tuning does not reach Davies placement, and intermediate sigma is unstable

`erf.nscbc_sigma` added: the LODI relaxation coefficient on the incoming
characteristic, blending J- between the far-field value (sigma=1, the shipped
Riemann behaviour) and the interior-extrapolated value (sigma=0, which reduces
ALGEBRAICALLY to variant 0 -- un_s = un_i, c_s = c_i, hence rho_s = rho_i and
p_s = p_i). Both endpoints were already measured over a full day, so the sweep
interpolates between known points. Default 1.0, so nothing else changes.

Endpoint controls before use: sigma=1 reproduced the Riemann arm's dt
BIT-IDENTICALLY at steps 1/2/5/20/100, drifting only by step 400 (round-off under
fast-math reassociation of the blend); sigma=0 matched the extrapolation arm to
the EOS round-trip the sigma path takes and variant 0 does not.

**Stability: both endpoints are stable, every intermediate value is not.**

| sigma | outcome | dt at death / day end |
|---|---|---|
| 0 (extrapolation) | completes 24 h | 0.167 |
| 0.03 | completes 24 h | 1.111 |
| 0.1 | **FPE at 7.9 h** | 1.656 (healthy to the last step) |
| 0.3 | **FPE at 1.9 h** | 1.643 (healthy to the last step) |
| 1 (Riemann) | completes 24 h | 1.303 |

Neither failure is a dt collapse -- both trap with a perfectly healthy timestep,
and time-to-death lengthens as sigma falls. That is the signature of the blend
itself being ill-posed: a mixture of the far-field and interior incoming
invariants is neither properly non-reflecting nor properly clamped.

**Placement is monotone in sigma and bounded by the endpoints, as predicted.**
Percentile-matched FSS, interior d >= 20:

| arm | 1 mm, 3 km | 5 mm, 3 km | 5 mm, 60 km |
|---|---|---|---|
| **Davies** | **0.822** | **0.743** | **0.841** |
| sigma=1 Riemann | 0.649 | 0.604 | 0.684 |
| sigma=0.03 | 0.628 | 0.547 | 0.634 |
| sigma=0 extrapolation | 0.586 | 0.478 | 0.589 |
| useful | 0.792 | 0.713 | 0.713 |

No NSCBC arm comes within 0.17 of Davies at any base rate or scale, and the
sigma arm sits BETWEEN its two endpoints -- the bounding argument holds. Since
percentile matching divides out exactly the bias that sigma tunes, this is the
expected result and it is now measured rather than argued.

**sigma = 0.03 is nonetheless the best NSCBC arm ever run on amplitude and
structure**, and by a wide margin:

| | Davies | sigma=0 | **sigma=0.03** | sigma=1 | MRMS |
|---|---|---|---|---|---|
| interior bias | 0.38x | 0.19x | **1.42x** | 2.63x | 1.0 |
| interior 8-19 km spectrum | 0.003 | 0.002 | **0.493** | 2.297 | 1.0 |
| islands, three western (mm) | 5.0/3.0/3.0 | 3.5/2.3/3.2 | **36.2/46.3/55.5** | 91.3/93.6/111.2 | 57.6/56.2/59.7 |
| inflow-wall ratio | 7.21 | 0.08 | **0.01** | 0.00 | 1.0 |
| N/S ratio | 67.1 | 4.7 | 60.5 | 29.0 | 9.70 |
| steps for the day | 62602 | 79560 | **55965** | 56095 | -- |

It is the only arm with a near-unity interior amount, the only one with
storm-scale variance within a factor of two of observed, and the only one that
puts realistic rain on the Channel Islands (36-56 mm against 58/56/60 observed,
where Davies has 3-5 and Riemann 91-111). It is also the cheapest day run. And
its placement is still 0.19 below Davies.

**Verdict, and the reason to stop.** Every structural deficiency #25 identified
in the Davies configuration -- starved interior, dead islands, featureless
spectrum, 210 mm/day wall band -- is fixed by NSCBC at sigma = 0.03. None of that
buys placement skill. Davies puts its (far too little) rain in the right places;
NSCBC puts a realistic amount in the wrong ones. These are independent failures
and the boundary condition addresses only the second.

Retained results, independent of how sigma landed:
1. The regime latch buys survival on its own -- extrapolation dies at 4.7 h
   without it, completes the day with it.
2. The outflow-wall excess is specifically a Riemann artifact -- peak ratio
   6.35 at the wall becomes 1.85 at 10-14 cells inboard on substitution.
3. The outflow condition sets the interior AMOUNT, swinging it 0.19x to 2.63x
   with the latch and the mass controller held fixed.
4. The LODI blend is unstable at intermediate sigma while both endpoints are
   stable, with time-to-death lengthening as sigma falls.

## 35. Davies' placement advantage is INHERITED ADVECTIVELY: absent for 4 h, saturated at one domain transit

Two discriminators on the same pair of scored days (Davies vs NSCBC sigma=0.03),
both percentile-matched so bias is divided out, both interior d >= 20.

**Distance-binned: FLAT.** Davies minus NSCBC, matched within each bin:

| bin (cells) | km from wall | n | 1 mm 3 km | 5 mm 3 km |
|---|---|---|---|---|
| 20-29 | 60-87 | 3760 | +0.225 | +0.160 |
| 30-39 | 90-117 | 2960 | +0.214 | +0.086 |
| 40-47 | 120-141 | 1792 | +0.205 | +0.506 |

A 9% decline over 80 km at 1 mm, and at 5 mm the deepest bin carries the LARGEST
advantage. No decay toward the interior.

**Hour-binned: the gap is BUILT, not present from the start.** Reference is the
MRMS 24 h total held fixed across hours, so absolute scores are depressed early;
both arms share it, so the gap is the meaningful quantity.

| hour | 1 mm 3 km gap | 5 mm 3 km gap |
|---|---|---|
| 1 | +0.005 | **-0.049** |
| 2 | **-0.010** | **-0.068** |
| 3 | **-0.024** | -0.010 |
| 4 | **-0.032** | +0.006 |
| 5 | +0.054 | +0.108 |
| 6 | +0.206 | +0.280 |
| **7** | **+0.258** | **+0.339** |
| 8 | +0.243 | +0.280 |
| 12 | +0.209 | +0.216 |
| 18-24 | +0.186 to +0.194 | +0.196 to +0.198 |

For the first four hours the two schemes are indistinguishable, and at 5 mm
**NSCBC is ahead**. The gap appears between hours 5 and 7, peaks at hour 7, then
relaxes to a plateau it holds for the rest of the day. Domain transit is
576 km / ~20 m/s = **~8 h**; the gap saturates at hour 7.

**Reading.** Both arms start from the same IC and are indistinguishable while
the interior is still IC-dominated. The advantage appears exactly as
boundary-influenced air fills the domain, and stops growing once it has. Davies'
interior placement is **inherited advectively from the air it admits**, not
generated by the interior dynamics and not imposed locally by the band.

**This reconciles the two tests, and shows why the distance test was the wrong
instrument.** Advected structure fills the interior roughly uniformly -- air
crosses the domain in ~8 h regardless of where it ends up -- so an advective
mechanism produces a FLAT distance profile. The distance test could not see the
mechanism it was designed to test; only the time axis resolves it.

**Consequence.** The lever is what arrives at the inflow wall: the driving data
(ERA5, 25 km, 3-hourly) and the width and strength over which it is imposed.
Honest caveat for the record: this partially rehabilitates the narrow-band idea
(item 2a of the anchoring plan), because a relaxation band conditions the air
entering the domain and that is precisely the mechanism now implicated -- band
WIDTH was never the operative variable, air conditioning was. The open question
that makes it worth or not worth doing: as the band recovers Davies' placement it
also recovers Davies, so the target configuration is the narrowest band that
buys placement while keeping NSCBC's amplitude, spectrum and islands.

**Geometry constraint, for the record.** On the 192x96 domain,
d = min(i, j, NX-1-i, NY-1-j) maxes at **47**. Distance-binned diagnostics
therefore reach only 141 km from the wall, a "60+ cells" bin is empty, and
"45-60" holds three rows. Any distance-resolved boundary diagnostic on this
domain is limited to that range; a wider domain would be needed to extend it.

## 36. THE RESULT: the relaxation ramp's grad(F) is simultaneously the anchoring mechanism and the artifact source

NSCBC (sigma=0.03) with `nscbc_keep_ramp=1` and `real_width=3` -- a 3-cell Davies
band behind the characteristic wall. Completed the day, mass -0.032%, dt healthy
(min 0.83), no FPE/NaN.

| | Davies | NSCBC no band | **NSCBC + 3-cell band** | MRMS |
|---|---|---|---|---|
| PM-FSS 1 mm 3 km | 0.822 | 0.628 | **0.819** | useful 0.792 |
| PM-FSS 5 mm 3 km | 0.743 | 0.547 | **0.747** | useful 0.713 |
| interior bias | 0.38x | 1.42x | **0.31x** | 1.0 |
| interior spectrum | 0.003 | 0.493 | **0.002** | 1.0 |
| islands, three western (mm) | 5.0/3.0/3.0 | 36/46/56 | **4.3/3.3/3.3** | 58/56/60 |
| inflow-wall ratio | 7.21 | 0.01 | **5.50** | 1.0 |
| d=0 band (mm/day) | 210.3 | 28.4 | **87.2** | 18.6 |

**A three-cell band buys 100% of the placement and 76% of the pathology.**
Placement matches Davies at both base rates and all five scales (5 mm is +0.004
to +0.014 ABOVE Davies). Every property NSCBC bought is simultaneously lost.

The hourly curve confirms the mechanism is inflow conditioning: against Davies
the no-band arm's gap grows to +0.258 by hour 7 and plateaus at +0.19, while the
3-cell arm's gap never exceeds +0.06 and closes to +0.003 by hour 17 -- it opens
and shuts inside the same hours 4-7 window in which the no-band gap opened.

**The mechanism, stated once.** The relaxation ramp's grad(F) is BOTH the
anchoring and the artifact. The spatial gradient of the relaxation coefficient is
what transmits the driver's placement into the interior; that same gradient is
what drives the spurious wall ascent that wrings the inflowing moisture out and
starves the interior downstream. They are not two effects to be traded against
each other -- they are one mechanism observed twice, which is why no band WIDTH
separates them and why width 3, already near the minimum, gives essentially all
of both.

This is also why the field's convention is to tolerate the relaxation zone's
noise and DISCARD the zone in analysis rather than try to remove it: the noise is
not a defect of the implementation, it is the anchoring being paid for. And it
makes spectral nudging the indicated construction rather than one option among
several -- it is the only formulation that anchors with NO spatial gradient.

Full-day arms scored to date, all same binary, same bracketing, all mass-bounded:

| arm | PM-FSS 1mm/3km | interior bias | spectrum | islands | inflow ratio |
|---|---|---|---|---|---|
| Davies | 0.822 | 0.38x | 0.003 | 5/3/3 | 7.21 |
| NSCBC sigma=1 (Riemann) | 0.649 | 2.63x | 2.297 | 91/94/111 | 0.00 |
| NSCBC sigma=0.03 | 0.628 | 1.42x | 0.493 | 36/46/56 | 0.01 |
| NSCBC sigma=0 (extrap) | 0.586 | 0.19x | 0.002 | 3.5/2.3/3.2 | 0.08 |
| NSCBC sigma=0.03 + width 3 | 0.819 | 0.31x | 0.002 | 4.3/3.3/3.3 | 5.50 |
| ERA5 (driver) | 0.938 | 1.50x | 0.008 | 24/22/32 | -- |

The driver outscores every model arm on placement at every scale, which is the
same statement from the other side: the boundary formulation controls how much of
the driver's placement the interior inherits, and the driver sets the ceiling.
Better driving data, not better BC formulation, is the lever on the ceiling.

### 36a. Spectral nudging: scope (NOT built)

**Integration point -- already exists.** ERF has a post-step Marchuk-split
relaxation path, `erf.realbdy_relax_split` (ERF_AdvanceDycore.cpp:371-436): it
assembles a tendency, applies it to the completed state with `MultiFab::Saxpy`,
then calls `MomentumToVelocity` to push the momentum increment back into
xvel/yvel -- the prognostic carrier between steps is the VELOCITY, so an
increment left only in momentum is silently discarded. That subtlety is already
solved there. Spectral nudging should reuse this hook verbatim, which also puts
it in the COSMO/TRAM ordering (relaxation excluded from the slow-mode forcing,
applied once per big step) rather than in the slow RHS where it would be held
fixed across every acoustic substep.

**Target field -- already exists.** `forecast_state_interp[lev][Vars::cons]` is
the ERA5 state already interpolated onto the ERF grid and already time-blended
between frames; the Davies relaxation consumes it today. No new I/O.

**FFT -- reusable, and non-periodic is supported.** `ERF_USE_FFT` is on in the
current build (FFT::PoissonHybrid symbols are in the shipped binary).
`amrex::FFT::R2X` (Submodules/AMReX/Src/FFT/AMReX_FFT_R2X.H) provides real
transforms with per-direction `FFT::Boundary::{periodic,even,odd}`, i.e.
cosine/sine series for a non-periodic limited-area domain -- which is what this
case needs and what WRF's own spectral nudging uses. ERF already has
`get_fft_bc()` (ERF_SolverUtils.H:56) mapping domain BC types onto those
boundary kinds. `R2X::forwardThenBackward(in, out, post_forward)` is exactly the
required primitive: transform, zero the modes above cutoff inside the callback,
transform back. Only constraint found: `domain.smallEnd() == 0`, which ERF
satisfies.

**Separability removes the per-level plumbing.** A horizontal-only filter does
not need per-slab 2-D transforms: the transform is separable, so a single 3-D
R2X that zeroes only (kx, ky) above cutoff and retains every kz is mathematically
identical to filtering each level in 2-D. The "above the boundary layer only"
requirement is then applied as a vertical taper on the resulting tendency in
PHYSICAL space, after the inverse transform, where it costs nothing and
introduces no horizontal gradient.

**Wavenumber cutoff -- and this is the finding that matters.** Lx = 576 km,
Ly = 288 km:

| k | lambda_x | lambda_y |
|---|---|---|
| 1 | 576 km | 288 km |
| 2 | 288 km | 144 km |
| 3 | 192 km | 96 km |

Conventional practice retains only scales the driver credibly resolves --
typically lambda > 600-1000 km, i.e. 3-5 modes on domains of several thousand km.
**On this domain the LARGEST available wavelength is 576 km, already below that
cutoff: by the standard criterion there are ZERO modes to nudge.** To nudge
anything the cutoff must be relaxed to lambda >= ~150-200 km -- defensible
against ERA5's ~110-170 km effective resolution (0.25 deg, ~4-6 dx) -- which
retains kx <= 3 and ky <= 2, six modes plus the mean.

So at 576x288 km spectral nudging degenerates to nudging the domain mean plus a
few tilts. The scale separation the technique exists to exploit -- driver-credible
large scales cleanly distinguished from model-generated small scales -- barely
exists at this domain size. That is a science risk, not an implementation risk,
and it should be weighed before any code is written.

**Variables:** u, v, theta. NOT moisture -- nudging qv would reimpose direct
moisture control, which is the lever this whole line of work is trying to avoid,
and standard practice excludes it.

**Vertical:** taper to zero below ~1 km and full above ~2 km. Prefer a fixed
height over a diagnosed PBL depth: a time-varying diagnosed depth introduces its
own moving gradient, which is the failure mode this construction exists to avoid.

**Effort:** ~300-400 lines (tendency assembly, cached R2X plan, spectral mask,
vertical taper, wiring into the split path) plus validation, call it 2-3 days.
Validation must include: round-trip identity with no modes filtered (to
round-off); a known-nonzero control nudging toward a synthetic field of known
spectrum, verifying only retained modes move; and tau -> infinity reproducing the
un-nudged run bit-for-bit. Main implementation risk is the R2X redistribution
cost per step on GPU -- the Poisson preconditioner already pays this, so it is
measurable rather than unknown.

### 36b. Width 2 confirms it: the transition is a STEP, not a trade-off curve

`real_width=2` -- one nudged cell behind the specified cell, the narrowest ramp
the scheme admits. Completed the day, mass -0.032%, no FPE/NaN.

| width | PM-FSS 1mm/3km | PM-FSS 5mm/3km | inflow ratio | interior bias | spectrum | islands (3 W, mm) |
|---|---|---|---|---|---|---|
| 0 (no band) | 0.628 | 0.547 | **0.01** | **1.42x** | **0.493** | **36/46/56** |
| **2** | **0.803** | **0.741** | 4.95 | 0.31x | 0.002 | 4.3/3.4/3.8 |
| 3 | 0.819 | 0.747 | 5.50 | 0.31x | 0.002 | 4.3/3.3/3.3 |
| 10 (Davies) | 0.822 | 0.743 | 7.21 | 0.38x | 0.003 | 5.0/3.0/3.0 |
| MRMS | -- | -- | 1.0 | 1.0 | 1.0 | 58/56/60 |

**Going from no band to the narrowest possible band recovers 89% of the placement
deficit and destroys essentially all of the structure gain, in one step.**
Spectrum collapses 0.493 -> 0.002, islands 36-56 mm -> 3-4 mm, interior bias
1.42x -> 0.31x, inflow-wall ratio 0.01 -> 4.95 (69% of Davies' 7.21). Widening
from 2 to 3 to 10 then changes little on either axis: +0.019 of placement and a
further 1.5x of wall band.

So this is not a trade-off curve with an optimum somewhere in the middle -- it is
a switch. Any spatial gradient of the relaxation coefficient, however narrow,
delivers the driver's placement AND the wall artifact together. Width is not the
operative variable; the existence of a gradient is. That is the sharpest possible
statement of 36's mechanism, and it closes the narrow-band line: there is no
width to find.

Note width 2 is still USABLE on its own terms -- 0.803 and 0.741 both clear the
believability thresholds (0.792, 0.713) -- it is simply Davies with extra steps.

Operational note: leg 1 of both the width-2 and width-3 runs ended on a
degenerate stop-truncation step (DT = 1.49e-08), and `pick_restart_chk.sh`
rejected that checkpoint and fell back to the previous periodic one. The item-31
guard fired in production, twice, as designed.

### 36c. Scale-filtered placement: DO NOT BUILD 2b. The advantage lives below the reachable scales, and inverts above them

Both model and observation fields low-passed with a 2-D DCT-II/III (a cosine
series on a non-periodic domain -- the same transform `amrex::FFT::R2X` with
`Boundary::even` performs, so the retained modes are exactly what a spectral
nudging implementation here could act on), then scored percentile-matched on the
interior. Controls: DCT round-trip exact to 1e-9, identity filter is the
identity, self-FSS == 1 at every cutoff.

Percentile-matched FSS at 3 km, interior d >= 20:

| cutoff | modes kx,ky | Davies | NSCBC no band | NSCBC width 2 | ERA5 |
|---|---|---|---|---|---|
| unfiltered | 191,95 | **0.822** | 0.628 | 0.803 | 0.938 |
| lambda >= 200 km | 5,2 | 0.720 | **0.891** | 0.824 | 0.952 |
| lambda >= 576 km | 2,1 | 0.680 | **0.802** | 0.787 | 0.902 |

At the 5 mm base rate the same reversal, more strongly: unfiltered Davies 0.743
vs NSCBC 0.547; at lambda >= 576 km, Davies 0.519 vs NSCBC 0.726.

**Davies' placement advantage does not merely vanish under filtering -- it
INVERTS.** Filtered to the scales spectral nudging could touch, the band-free
NSCBC arm is the BETTER placed field by +0.17 (1 mm) and +0.21 (5 mm), and it
sits within 0.06 of ERA5's own filtered score (0.891 vs 0.952). At those scales
NSCBC is already 94% of the way to the driver.

**So spectral nudging would act precisely where NSCBC needs no help.** The whole
of Davies' unfiltered advantage (0.822 vs 0.628) lives at wavelengths BELOW
200 km, which a nudging cutoff excludes by construction -- that is what the
cutoff is for. Building it would spend 2-3 days to improve a quantity already
near its ceiling, and would leave the deficit that matters untouched. **Do not
build.** This is the decision rule the test was set up to answer, answered
against the build.

**A caveat this raises about the campaign's central number.** Davies wins
unfiltered and loses filtered, and Davies has almost no small-scale variance
(interior spectral ratio 0.003 against NSCBC's 0.493). Binarising a nearly smooth
field at its own quantile produces large contiguous blobs, which overlap a
smoothed observation field more readily than a structured field does even when
bias is matched. So part of Davies' sub-200 km "placement advantage" may be the
metric rewarding smoothness rather than skill. Percentile matching removes
amplitude bias; it does not remove variance bias. Any future use of these
numbers should carry that qualification, and a variance-matched or
object-oriented score would be the way to settle it.

### 36d. The smoothness caveat is FALSIFIED. Davies' placement advantage is real -- and it is an advantage over an envelope, not over a storm

Three tests, all post-processing, against the three scored arms plus ERA5.

**1. Variance-matched scoring -- caveat withdrawn.** Each field was Gaussian-
smoothed until its 8-19 km spectral ratio matched Davies' 0.0031 (NSCBC no-band
needed sigma = 1.94 cells, ERA5 sigma = 0.71, the two already-smooth arms needed
none), then re-scored. Controls: sigma = 0 reproduces the native score exactly.

| arm | native | variance-matched | change |
|---|---|---|---|
| Davies | 0.822 | 0.822 | +0.000 |
| NSCBC no band | 0.628 | 0.624 | **-0.003** |
| NSCBC width 2 | 0.803 | 0.803 | +0.000 |
| ERA5 | 0.938 | 0.937 | -0.000 |

Smoothing NSCBC to Davies' variance moves its score by three thousandths. **The
sub-200 km advantage is not a smoothness artifact.** Percentile matching at a
fixed base rate turns out to be robust to this: the quantile adapts, so the
retained set moves very little under smoothing. The caveat raised in 36c is
withdrawn -- it was a real risk, it was tested, and it does not hold.

**2. MODE-style object verification -- confirms the placement ranking, and adds
the finding that matters.** Convolve (sigma 1 cell), threshold, label, discard
specks under 10 cells, match on centroid distance + area similarity. Control: a
field matched against itself gives 100% matched at 0 km.

| threshold | arm | n_obj | matched | median displacement | area bias | intensity bias |
|---|---|---|---|---|---|---|
| 5 mm | ERA5 | 1 | 1/1 | **25.6 km** | 1.40 | 0.78 |
| 5 mm | NSCBC width 2 | 3 | 1/1 | 42.2 km | 0.20 | 0.36 |
| 5 mm | Davies | 2 | 1/1 | 42.7 km | 0.25 | 0.41 |
| 5 mm | NSCBC no band | 1 | 1/1 | **139.9 km** | **0.80** | **1.57** |
| 15 mm | Davies | 2 | **0/2** | -- | -- | -- |
| 15 mm | NSCBC width 2 | 2 | **0/2** | -- | -- | -- |
| 15 mm | NSCBC no band | 1 | 1/2 | 70.8 km | 1.17 | 1.84 |
| 30 mm | Davies | 1 | **0/1** | -- | -- | -- |
| 30 mm | NSCBC no band | 1 | 1/1 | 72.1 km | 1.60 | 2.02 |

An independent method, sharing no machinery with FSS, reproduces the placement
ranking at 5 mm: Davies and width-2 put the storm envelope ~42 km from observed,
NSCBC no-band ~140 km. So the placement advantage is real.

**But the same table shows what that advantage is over.** At 5 mm Davies'
matched object is **25% of the observed area and 41% of its intensity**; at 15 mm
and 30 mm Davies and width-2 produce **no matchable objects at all** -- they have
no heavy-rain systems where MRMS has them. NSCBC no-band is the only model arm
that produces heavy-rain objects, at 0.80-1.60 area bias and ~70 km displacement.

**3. ERA5 control.** Variance matching moves ERA5 by -0.000, and ERA5 has both
the smallest displacement (25.6 km) and the most realistic object area (1.40).
The metric is not implicated: a field that is smooth relative to MRMS can still
score well or badly on its merits. ERA5 outscoring every model arm on placement
AND object attributes is unchanged and remains the campaign's least comfortable
result.

**Limitation, stated plainly:** MRMS forms only 1-2 objects in the interior at
these thresholds, so each displacement figure is a single realisation, not a
distribution. The agreement with the independent FSS ranking is what gives it
weight, not its own sample size.

## 37. Which configuration is the better forecast, and the headline restated

**Restated headline.** The earlier form -- "Davies places precipitation better;
NSCBC gets the amount and structure right" -- survives, with the smoothness
caveat withdrawn (36d) and one clause added that changes what it means:

> Davies places the light-rain ENVELOPE better, and produces no storm inside it.
> NSCBC produces a storm of nearly the right size and intensity, ~100 km from
> where it belongs. Neither beats the 25 km driver at either job.

The added clause is the object result: at 5 mm Davies' storm is a quarter of the
observed area and 41% of its intensity, and at 15 mm and 30 mm it has no
matchable objects at all. Its placement advantage is an advantage over an
envelope, not over a system. That was invisible to FSS, which scores a binary
coverage pattern and therefore cannot distinguish "right place, no storm" from
"right place, right storm".

**Which is the better precipitation forecast -- it depends on the use, and this
should not be collapsed into one number:**

- For **where it rained at all** (light-rain footprint, flood-watch framing):
  Davies, or equivalently NSCBC + a 2-cell band. 42 km envelope displacement,
  percentile-matched placement 0.822/0.743, clearing the believability threshold.
- For **how hard it rained and over what area** (QPF, hydrology, anything driven
  by intensity): NSCBC no-band, sigma = 0.03. Interior bias 1.42x against
  Davies' 0.38x, spectral variance 0.493 against 0.003, islands 36-56 mm against
  Davies' 3-5 (observed 58/56/60), and the only arm producing heavy-rain objects
  at all.
- For **either purpose on this domain and day: ERA5**, the 25 km driver, beats
  every model arm on placement (0.938/0.890), on envelope displacement (25.6 km)
  and on object area (1.40). The 3 km hindcast does not add value over its own
  driver at any scored quantity.

**What the campaign established, independent of which arm is preferred:**
1. The relaxation ramp's grad(F) is simultaneously the anchoring and the artifact;
   width is not the operative variable, the existence of a gradient is (36, 36b).
2. Interior placement is inherited advectively from the inflow, saturating at one
   domain transit (35).
3. Above ~200 km the placement ranking inverts and NSCBC is within 0.06 of the
   driver, so spectral nudging would act where nothing needs fixing (36c).
4. The NSCBC regime latch buys survival; the Riemann solve buys the outflow-wall
   excess; the outflow condition sets interior amount (33, 34).

**The lever these four point at is not the boundary condition.** Placement is
inherited from what arrives at the wall and is capped by the driver; amount and
structure are set by the interior scheme and are the one thing the 3 km run does
better than ERA5. A configuration that had NSCBC's amplitude and structure with
Davies' envelope placement would beat both -- and 36/36b say no boundary
formulation with a spatial gradient can deliver it. Driving data at higher
resolution is the lever on placement; the interior physics is the lever on
amount.

### 36e. Width sweep completed: pathology is monotone in width, placement is a STEP then flat with scatter

`real_width=5` completed the day, mass -0.032%, no FPE/NaN. Full series, all same
binary, same bracketing, same scoring:

| width | PM-FSS 1mm/3km | PM-FSS 5mm/3km | inflow-wall ratio | interior bias | spectrum | islands (3 W, mm) |
|---|---|---|---|---|---|---|
| 0 (no band) | 0.628 | 0.547 | **0.01** | **1.42x** | **0.493** | **36/46/56** |
| 2 | 0.803 | 0.741 | 4.95 | 0.31x | 0.002 | 4.3/3.4/3.8 |
| 3 | **0.819** | **0.747** | 5.50 | 0.31x | 0.002 | 4.3/3.3/3.3 |
| 5 | 0.784 | 0.713 | 6.53 | 0.33x | 0.003 | 4.2/3.0/3.0 |
| 10 (Davies) | 0.822 | 0.743 | 7.21 | 0.38x | 0.003 | 5.0/3.0/3.0 |

**The pathology is monotone in width** -- the inflow-wall ratio climbs 0.01 ->
4.95 -> 5.50 -> 6.53 -> 7.21 without exception, so a wider band always costs more
wall artifact.

**Placement is not monotone.** It steps from 0.628 to 0.803 between width 0 and
width 2, then wanders inside 0.784-0.822 across widths 2, 3, 5 and 10 with no
ordering: width 5 (0.784) scores BELOW width 3 (0.819) and below Davies (0.822),
and width 3 scores above width 5 and width 2. The spread across the plateau
(0.038 at 1 mm, 0.034 at 5 mm) is the same size as the differences between
adjacent widths, so those differences are not systematic.

This strengthens 36b rather than complicating it: once any relaxation gradient
exists, placement is set and width does not tune it, while width continues to buy
wall artifact monotonically. **Width is strictly a cost knob with no placement
benefit beyond its existence.** There is no optimum to find, and width 5 -- the
one point that could have shown a mid-range peak -- shows the opposite.

Interior structure is unchanged across every non-zero width: bias 0.31-0.33x,
spectrum 0.002-0.003, islands 3.0-4.3 mm against 58/56/60 observed. All of it is
lost at width 2 and none of it returns.

## 38. CONUS404 driving data: SCOPE (not built). One verified blocker

### Verified from store metadata

Store: `s3://hytest/conus404/conus404_hourly.zarr` via
`https://usgs.osn.mghpcc.org/` (anonymous, no egress fees), HyTEST catalog
`conus404-catalog.yml`. Read `.zmetadata` directly; lat/lon chunks decoded
locally (zstd level 9, 175x175 chunks).

- **Native hybrid model levels: YES for the simulation.** The store carries the
  full WRF hybrid vertical machinery -- `ZNU(50)`, `ZNW(51)`, hybrid coefficients
  `C1H..C4H(50)` and `C1F..C4F(51)`, plus `MUB(y,x)`, `PB(50,y,x)`,
  `PHB(51,y,x)`. Those arrays exist only for a run on hybrid sigma-pressure
  model levels, not pressure levels. Verified, not inferred.
- **Time: 394488 hourly steps, 1979-10-01 00Z -> 2024-09-30 23Z.** 2023-01-09 00Z
  is index 379344, in range. Our scored day IS covered.
- **Grid: Lambert Conformal, DX=DY=4000 m, TRUELAT1/2 = 30/50,
  STAND_LON = -97.9, 1367 x 1015.** `TITLE = OUTPUT FROM WRF V3.9.1.1 MODEL`.
- **Domain coverage including the offshore western edge: CONFIRMED.** Decoded
  chunk (2,0) spans lon -129.00..-119.57, lat 29.36..37.23 and contains **5429
  cells inside our footprint** (-123.47..-117.22, 32.06..34.70). ~5.5 deg of
  margin west of our -123.5 edge.

### The blocker: the accessible store has NO 3-D atmospheric state

Every atmospheric field in it is single-level despite carrying WRF names:
`U['time','y','x_stag']`, `V['time','y_stag','x']`,
`QVAPOR/TK/T2/Q2/PSFC/TSK['time','y','x']`. **Absent entirely: W, T (theta), P,
PH, QCLOUD, QRAIN, QICE, QSNOW, QGRAUP.** The only time-varying 3-D variables are
soil (`SMOIS`, `SH2O`, `TSLB`, 4 layers). `bottom_top` exists solely to carry the
time-invariant coordinate arrays. The catalog says so in words -- "CONUS404 Hydro
Variable subset" -- and the metadata confirms it.

So the premise is true of the SIMULATION and false of the PRODUCT: the
hydrometeors that motivate the move are exactly what this store omits. Obtaining
full 3-D wrfout (NCAR RDA/GDEX, the ScienceBase release, or a request to the
producers) is UNRESOLVED and is the gate on everything else.

### What our pipeline needs -- cheaper than expected

ERF reads only the custom `.bin` frame format (ERF_ReadCustomBinaryIC.H): 4 int32
`(nx, ny, nz, ndata)`, then `lat[nx*ny]`, `lon[nx*ny]`, `x[nx]`, `y[ny]`,
`z[nz]`, then `ndata` blocks of `nx*ny*nz`. **So the work is a new CONVERTER, not
a new ERF reader path -- ERF needs zero changes** (with one exception below).
That turns a C++ project into a Python one.

Which fixes survive:

| fix | fate under CONUS404 |
|---|---|
| LCC pinning + frame-coordinate residual gate (32) | SURVIVES, and matters more: this becomes Lambert->Lambert reprojection, not lat/lon->Lambert. Reuse the residual check as the acceptance test. |
| `hindcast_frame_from_T` | UNNECESSARY. It exists because ERA5 gives T and we need theta/rho; WRF gives perturbation theta plus P/PB directly. Emit theta/rho in the converter, leave the knob off. |
| `hindcast_sfc_anchor_file` surface anchor | UNNECESSARY. It exists because ERA5's lowest data level is 150.19 m (measured from our own frames: z = [0, 150.19, 364.42, 583.47]); WRF model level 1 is ~10-30 m. |
| masked SST interpolation | REWORK, not delete. `TSK` and `LANDMASK` are present; `SST` is not. TSK-over-water substitutes for SST and the fallback-to-nearest-valid logic still applies. |
| erftools ~300 m level displacement (19f/21) and the fabricated duplicate bottom level | BYPASSED. Both are erftools artifacts; our own converter never creates them. |

**The one place ERF source MUST change.** `ndata` is hard-capped at 8
(ERF_ReadCustomBinaryIC.H:70) and our frames already use all 8 -- verified,
header reads `40 22 38 8` for rho,u,v,w,theta,qv,qc,qr. So qc/qr come free, but
**adding qi/qs requires raising the cap and extending the reader and its
consumers**. "Resolved hydrometeors" beyond cloud+rain is not a data-only change.

### What we gain, and what caps it

- Hydrometeors: our ERA5 hydrometeor test failed only because ERA5 carries
  essentially none (measured 4.0e-4 against 4.8e-3 the model generates).
  Contingent on 3-D access.
- No sub-150 m extrapolation: removes the whole surface-anchor class.
- Effective driver resolution goes from ~110-170 km (ERA5, 4-6 dx at 0.25 deg) to
  ~16-24 km. **This is the right lever by our own measurement**: the scale-filtered
  test (36c) put the model within 0.06 of ERA5 above 200 km and located the entire
  placement deficit BELOW 200 km -- where a 4 km driver has information and ERA5
  has none.
- Cap on the upside: CONUS404 is itself ERA5-driven, so at scales ERA5 resolves it
  inherits ERA5's placement. The gain is confined to roughly 4-200 km.

### Effort

| phase | days | note |
|---|---|---|
| 0. Resolve 3-D access | 0.5-1 | GATING. Deliverable: one hour of 3-D wrfout over our footprint on disk + variable inventory. |
| 1. Converter -> .bin | 3-5 | 3-D read, theta/rho, Lambert->Lambert reprojection, native levels, exact frame layout, residual gate |
| 2. Surface frames | 1-2 | 6-field frame from TSK/LANDMASK/T2/Q2/PSFC/albedo; masked fill reworked |
| 3. Validation + the two scored days | 2-3 | init sanity, 6-h run, then Davies and NSCBC no-band, identical scoring |
| **total** | **8-13** | if phase 0 resolves; **unbounded if it does not** |

Volume estimate: our footprint is 144 x 72 CONUS404 cells, so 50 levels x ~10
fields is ~207 MB/hour, 1.9 GB for 9 three-hourly frames or 5 GB hourly -- IF
server-side subsetting exists. Without it, full-CONUS fields are ~27 GB/hour.

**Recommendation: do not start phase 1 until phase 0 produces a real 3-D file.**
Every other cost here is contingent on it, and the accessible store has been
verified not to contain what the campaign needs.

### 38a. Phase 0 RESOLVED: the 3-D distribution exists and is subsettable -- but it ends 2022-09-30

**1. NCAR GDEX d559000 -- the 3-D archive, confirmed.** DOI 10.5065/ZYY0-Y036,
815.67 TB total, 43 water-year products at ~18.9 TB each, netCDF4. Both `wrf2d`
and `wrf3d` files, **744 of each per month = hourly**. OPeNDAP DDS of
`wrf3d_d01_2021-10-01_00:00:00.nc` gives the full 3-D atmospheric state:

    U[50][1015][1368]   V[50][1016][1367]   W[51][1015][1367]
    TK[50]  P[50]  P_HYD[50]  Z[51]
    QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP  (all [50])
    QNICE QNRAIN CLDFRA REFL_10CM O3RAD
    XLAT/XLONG (+ _U/_V staggered)  lev[50]  ilev[51]

Better than the scope assumed on two counts: **Z is supplied directly** as
geopotential height on staggered levels, so no PH/PHB reconstruction and exact
per-column level heights; and **TK is diagnosed temperature**, so theta comes from
TK and P offline in the converter -- where it can be verified -- instead of inside
ERF at runtime.

**2. Subsetting: YES.** The THREDDS catalog for d559000 advertises
`OpenDAP`, `NetcdfSubset` (ncssGrid), `CdmRemote` and `cdmrFeatureGrid`.

**3. Volume: comfortably viable, and my scope figure was pessimistic.** Our
footprint is 144 x 72 CONUS404 cells; 50 levels x 4 bytes = 2.07 MB per 3-D
field, so the ~10 fields we need are **20.7 MB/hour**, i.e. **187 MB for a
9-frame 3-hourly day**. (Correcting 38: full-CONUS is ~4.4 GB/hour for 16 3-D
fields uncompressed, not 27 GB/hour -- I over-counted the field list.)

**4. THE BLOCKER, verified from the dataset page:**

    Temporal Range: 1979-10-01 00:00 +0000 to 2022-09-30 23:00 +0000 (Entire dataset)

Last product is `202110-202209` (WY2022). **Jan 9 2023 -- the day every arm in
this campaign is scored on -- is outside the 3-D archive.** Meanwhile the HyTEST
2-D hydro subset runs to 2024-09-30, which implies the simulation WAS extended
past WY2022 and the 3-D output exists but is not published on GDEX.

**So all three prior avenues are settled:** GDEX has 3-D and subsetting but stops
at Sep 2022; the USGS/HyTEST products past that date are the 2-D hydro subset
only (38); ScienceBase lists the same release lineage. No accessible 3-D
distribution covers Jan 2023.

**Two ways forward, both needing a decision:**

- **(a) Move the scored day inside the archive.** MRMS Pass-2 QPE starts ~2014
  and the 3-D archive ends Sep 2022, giving an **8-year window (2014 - Sep 2022)**
  of MRMS-verifiable SoCal AR events. Cost: the comparison must be internally
  consistent, so it needs ERA5-driven Davies + NSCBC no-band on the NEW day as
  the control, plus the two CONUS404-driven arms -- 4 full days of GPU (~6 h) plus
  MRMS and frame prep. The Jan-9-2023 absolute numbers would not transfer, but
  every within-day conclusion in this file would still stand as its own result.
- **(b) Ask the producers** (Rasmussen / Gutmann at NCAR, or the USGS stewards)
  whether post-WY2022 3-D wrfout is retrievable. The 2-D extension to 2024 is
  direct evidence it was produced.

Not proceeding to the converter: the precondition -- 3-D data for the scored day
-- is not met, and building a converter before knowing which day it targets would
be work at risk.

## 39. New case 2020-12-28, CONUS404-driven. Step 0: all three checks PASS, with two scoring caveats

### 1. CONUS404 3-D access: UNBLOCKED for this date

The Jan-2023 blocker (38a: GDEX d559000 ends 2022-09-30) does not apply --
**2020-12-28 sits in WY2021**. Verified:

- All **24 hourly `wrf3d_d01_2020-12-28_HH:00:00.nc`** files present in
  `d559000/wy2021/202012` (744 wrf2d + 744 wrf3d for the month).
- 3-D state confirmed by OPeNDAP DDS: `U/V/W`, `TK`, `P`, `P_HYD`, `Z`,
  `QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP`, `CLDFRA`, `REFL_10CM`.
- **Server-side subsetting PROVEN, not just advertised**: anonymous OPeNDAP
  hyperslab requests against this exact file returned HTTP 200 for both
  `XLAT/XLONG[500:502][300:302]` and a real 3-D `QCLOUD[0][0:1][500:501][300:301]`
  slab. No authentication.
- Practical note for the converter: THREDDS/Tomcat rejects raw `[`/`]`/`:` in the
  request target -- percent-encode as `%5B %5D %3A`, and pass `curl -g` to stop
  glob expansion. Both cost a debugging cycle here.
- Volume unchanged from 38a: ~20.7 MB/hour for our 144x72 footprint, ~500 MB for
  24 hourly frames.

The zarr Hydro subset remains 2-D only (38) and is not usable; native subsettable
wrfout via OPeNDAP is the route.

### 2. Domain overlap: our box is fully inside wrfout d02

`wrfout_d02_2020-12-28_00_00_00` (27.7 GB): 745 x 765 at **DX = 1500 m**, 54
levels, Lambert TRUELAT1/2 = 30/60, STAND_LON = -120, **24 hourly times
2020-12-28_00 through _23** -- exactly our window. Full 3-D state present
(T, THM, P, PB, PH, PHB, QVAPOR, QCLOUD, QRAIN, QICE, QSNOW, QGRAUP).

    d02 extent   lat 31.827 .. 42.570   lon -126.963 .. -112.905
    our box      lat 32.064 .. 34.703   lon -123.470 .. -117.217
    fully inside: lat YES, lon YES;  73,729 d02 cells inside our box

**Caveat: the southern margin is thin.** Our south edge is 0.24 deg = **26 km**
from d02's south boundary. d02's own lateral relaxation zone at 1.5 km is
typically 5-10 cells (7.5-15 km), so our southern boundary samples d02 within
~11-18 km of d02's own sponge. Clear, but not generously. Every other margin is
large (west 3.5 deg, north 7.9 deg, east 4.3 deg).

### 3. MRMS: full coverage, but the event is a LAND event and the peak is unscoreable

23 hourly Pass-2 files (01Z..23Z) covering 00Z->23Z exactly. Regridded to our grid:
**18432/18432 valid cells**; native-in-footprint vs regridded control **7.59 vs
7.63 mm** (0.5%).

| mask | n | mean | p90 | max |
|---|---|---|---|---|
| whole domain | 18432 | 7.63 | 21.0 | 66.3 |
| interior d>=20 | 8512 | 6.03 | 15.1 | 36.1 |
| band d<10 | 5360 | 9.03 | 25.3 | **66.3** |
| **interior d>=20 LAND** | **509** | **22.00** | 30.3 | 36.1 |
| all LAND | 2831 | 22.29 | 34.0 | 66.3 |

**Confirmed: the domain maximum (66.3 mm, 34.486N -119.802W) sits at d = 8 --
inside the relaxation band and outside d >= 20.** The event peak will not be
scored, as anticipated.

**Caveat A -- the interior land sample is 509 cells.** "MRMS primary over land"
restricted to the scored interior leaves 6% of the 8512 cells used on Jan 9 2023.
That is enough for bias, PCC and object work, but FSS on 509 cells at 60 km
neighbourhoods is nearly a single-neighbourhood statistic and should be read with
that in mind.

**Caveat B -- the ladder does not survive intact on interior land.** Exceedance
counts:

| mask | n | >=1mm | >=5mm | >=15mm | >=30mm |
|---|---|---|---|---|---|
| interior d>=20 (all) | 8512 | 84.3% | 40.6% | 10.2% | **0.7% (63)** |
| interior d>=20 LAND | 509 | **100%** | **99.2%** | 79.8% | 12.4% (63) |
| all LAND | 2831 | 100% | 99.5% | 73.2% | 20.9% |
| whole domain | 18432 | 86.2% | 46.1% | 16.0% | 3.3% |

On interior LAND the 1 mm and 5 mm rungs are **degenerate** -- 100% and 99.2% of
cells exceed, so those rungs carry no information there. The 15 and 30 mm rungs
are the informative ones over land, which is the opposite of Jan 9 2023 and
exactly the reason for extending the ladder. Conversely on the full interior the
30 mm rung has only 63 cells (0.7%) and is not scoreable.

**Recommended scoring split, to be confirmed:** report 1/5 mm on the full
interior (d>=20, all surfaces) where they are informative, and 15/30 mm on land
(all-land n=2831 preferred over interior-land n=509 for sample size, with the
band caveat stated). Percentile-matched FSS sidesteps the degeneracy by
construction and should carry the like-for-like comparison.

**Verdict: proceed to Step 1.** No blocker. The two caveats are scoring-design
issues, not data problems.

### 39a. wrfout d02 validated against MRMS over land: USABLE as a reference, with two limits

**Time convention, both sides, stated before comparing.**
MRMS Pass-2 `..._01H_..._HHMMSS` is a PAST-HOUR accumulation valid at its stamp,
so file `-010000` covers 00:00-01:00 and `-230000` covers 22:00-23:00. The 23
staged files therefore sum to exactly **00:00Z -> 23:00Z**. The wrfout carries 24
hourly stamps `2020-12-28_00:00:00 .. _23:00:00` and RAINNC is a RUN-TOTAL
(non-zero at t0: mean 15.7 mm, so the run started before our window), so the
window accumulation is the difference between the last and first stamps -- also
exactly **00:00Z -> 23:00Z**. Identical periods, no offset.

**The bucket had to be handled.** `BUCKET_MM = 100` and `I_RAINNC` is live:
RAINNC wraps at 100 mm. Total = `RAINNC + 100*I_RAINNC`. A naive `RAINNC(23) -
RAINNC(0)` gives **517 cells at -99.8 mm**; bucket-corrected the field is
min 0.000, max 78.75, mean 4.114 mm. `RAINC = RAINSH = 0` everywhere, as expected
for a 1.5 km convection-permitting nest, so grid-scale precipitation is the whole
signal.

**Regridding** is an area-average of the 1.5 km field into our 3 km cells (median
4 source cells per target, all 18432 targets filled), not nearest-neighbour
sampling, which would have thrown away half the information. Conservation check:
d02 mean over our box 9.526 mm vs regridded 9.554 mm (0.3%).

**Agreement over land:**

| mask | bias | PCC | RMSE | FSS 15 mm (3->60 km) | FSS 30 mm (3->60 km) |
|---|---|---|---|---|---|
| all land, n=2831 | **1.09x** | **0.756** | 9.03 mm | 0.912 -> 0.989 (useful 0.866) | 0.535 -> 0.857 (useful 0.604) |
| interior land, n=509 | **0.93x** | **0.773** | 6.66 mm | 0.902 -> 0.989 (useful 0.899) | 0.286 -> 0.829 (useful 0.562) |

FSS clears the believability threshold at 1, 5 and 15 mm at **every** scale on
both masks. The 30 mm rung clears from 9 km (all land) / 15-30 km (interior land).

**The disagreement is AMPLITUDE IN THE TAIL, not placement.** Quantiles over land
(mm): wrfout 5.20/15.68/23.75/32.51/40.22/60.62 against MRMS
8.90/14.50/23.10/28.90/34.00/45.91 at q10/25/50/75/90/99 -- the middle matches
closely and the **upper tail overshoots** (q99 60.6 vs 45.9, land max 77.0 vs
66.3) while the low end undershoots. 76% of land cells agree within 10 mm; of the
9.1% that disagree by more than 15 mm, **244 of 258 are wrfout-wetter**, at median
wall-distance d = 8, i.e. concentrated on the northern coastal ranges. PCC 0.76
and FSS >= 0.9 at 1/5/15 mm say the pattern is right; the 30 mm weakness at 3 km
that recovers by 30-60 km says the heaviest cells are slightly displaced and too
intense -- the familiar convection-permitting over-intensification of orographic
peaks. Interior-core spectral ratio wrfout/MRMS 8-19 km = 1.677 (68% more
small-scale variance), though that core includes ocean and so partly measures the
artifact under test.

Figures: `figs/wrf_vs_mrms_0-12in.png`, `figs/wrf_vs_mrms_0-3in.png`.

**Verdict: usable as a full-domain reference.** Over land it is within 9% in the
mean, correlates at 0.76, and is believable at 1/5/15 mm at all scales.

**Limit 1, as required: offshore trust is extrapolation, not measurement.** The
validation could only be done over land (2831 of 18432 cells, 15% of the box).
Nothing here measures wrfout's skill over the ocean, which is 85% of the domain
and the region the reference is wanted for. The land result makes offshore
plausible; it does not make it verified.

**Limit 2, structural: d02 cannot be both the reference and an arm.** Step 3
lists wrfout d02 as the nesting counterfactual. If it is also the reference it
scores 1.0 against itself by construction and the arm is degenerate. These must
be separate tables: MRMS-over-land as truth with d02 as an arm, OR d02 as a
full-domain reference with d02 excluded from the arms. Scoring ERF against d02
also measures agreement with a different model -- ERA5-driven WRF at 1.5 km,
carrying its own wet tail and +68% small-scale variance -- not skill against
observations.

## 40. 30p history bisect and the 17 -> 90 min performance attribution (no runs)

### 30p: the failure did NOT always exist, and the bisect range is 3 days wide

Twelve mechanism hypotheses were burned. Searching history instead:

**Earliest continuous 24 h with no restart-bracketing: commit `7f917ebd`
(2026-07-23), "Real-BC inflow for HindCast + vertical-advective dt fix: first
complete 24-h run", consolidated in `83228eec` the same day -- 38,550 steps,
17.3 min, mass -0.26%/24 h, on the 128x64x32 deck.** Further clean 24-h
completions follow through items 22-25 (all scored on 24-h accumulations, so all
completed) and at line 1164 ("Both ran clean: 24 h complete, exit 0").

**First recorded failure: item 27, measured 2026-07-26 -- "the 192x96 Jan-9 run
reached 18.002 h (46,270 steps) and then went non-finite".** The domain-change
commit is `1ab7739a` (2026-07-26), "Domain chain steps 1-7: 192x96 domain on a
pinned projection".

So the bisect range is **`7f917ebd` .. `1ab7739a`**, ~10 commits over three days,
with four candidate boundaries:

| commit | date | change |
|---|---|---|
| `87c25a24` | 07-24 | production launch config: **cfl 0.3 -> 0.2**, NaN tripwire, max_dt |
| `617017b0` | 07-26 | five audit SEV-1/SEV-2 corrections wired into the deck |
| `aca969fd` | 07-26 | validated-configuration manifest |
| `1ab7739a` | 07-26 | **192x96 domain + pinned projection** |

**The failure post-dates the domain change and does not predate it.** It also
post-dates the IC fixes (items 19-21, scored on completing runs), the mass
controller (item 30, later still), and the frame-count change (frames stayed
9 three-hourly across the boundary -- so rotation 7 existed in the runs that
COMPLETED).

That last point is the sharp one: **the 128x64 runs also crossed rotation 7 and
survived.** So "rotation 7" is not sufficient by itself; it is rotation 7 *on the
192x96 domain*. Item 27b already noted the corner |w| precursor existed at 128x64
(7.75 vs interior 2.09) on runs scored clean.

**A git bisect is the wrong instrument here.** Deck and code changed together:
checking out a pre-`1ab7739a` commit gets the 128x64 deck back, so the bisect
would only re-confirm that 128x64 does not fail. The clean experiment is a 2x2 of
{old code, new code} x {128x64 deck, 192x96 deck}, and the single most
informative cell is **current code + 128x64 deck, 24 h continuous** (~20 min GPU,
one run). If it completes, the failure is domain-dependent and no commit bisect
is needed at all; if it dies at 18 h, the failure is code-dependent and the
four-boundary bisect above is worth its ~4 runs.

**Does CONUS404 driving sidestep it? Expectation: NO, and the premise needs
correcting.** Under the Step-1 design the converter emits the existing `.bin`
frame format and ERF is unchanged -- so `ERF_WeatherDataInterpolation` and the
rotation machinery are **bit-identical code on an identical format**. Only frame
CONTENT changes, and CADENCE is our choice. Emitting 3-hourly frames reproduces
the rotation structure exactly and should reproduce the failure; emitting hourly
frames TRIPLES the rotations, and since the suspect mechanism (30r) is
allocation-history dependence in the ~13 `Gpu::DeviceVector`s that
`FillForecastStateMultiFabs` churns per rotation, more rotations should make it
more likely, not less. The one thing that would genuinely sidestep it is a
different domain size, per the bisect finding above.

### Performance: 17.3 -> 64.6 min is FULLY attributed, with no residual

Baseline `83228eec` (128x64x32): 38,550 steps, 17.3 min, 0.0269 s/step.
Current (192x96x48, Davies A/B day): 62,602 steps, median **0.0619 s/step**,
implied 64.6 min. (Observed wall ~75 min including I/O, checkpointing and the
leg-2 restart.)

| factor | ratio | source |
|---|---|---|
| cells 262,144 -> 884,736 | 3.375x | `1ab7739a` domain change |
| ...but measured s/step | **2.301x** | per-cell-per-step cost **IMPROVED 0.68x** -- the 262K baseline was launch-latency bound |
| step count | **1.624x** | cfl 0.3 -> 0.2 (`87c25a24`) is 1.50x of this; max_dt=2.5 and domain winds the remaining 1.08x |
| **product** | **3.74x** | 17.3 x 3.74 = **64.6 min** -- matches |

**There is no unexplained 1.5x.** The 1.5x is the cfl change, and the domain
change is *sub*-linear in cells because GPU efficiency rises with problem size.

Levers, measured not guessed:
- **cfl 0.2 -> 0.3: the one large win, already validated and never adopted.**
  A full 18 h leg ran 31,415 steps against 46,586 (**33% fewer**) and the full day
  completed with mass -0.032%, identical to cfl 0.2, no NaN/FPE. Would take
  64.6 -> ~43 min. Caveat: cfl 0.2 was set by Hurricane Hilary (2023-08-20) where
  0.3 died in ~9 model-minutes, so 0.3 is validated for AR regimes and NOT for the
  convective summer case. Per-regime cfl, not a global change.
- **`max_gpu_streams` 1 -> 4: unavailable.** Retested this session on the CURRENT
  binary in the production config: exit 127, **zero steps**, FPE at launch. The
  1.2x is not a tunable cost; it is a hard failure.
- **`check_for_nans_int` 10 -> 100: ~9%.** Measured this session: int=1 costs 3x
  (4.7 vs 14.9 steps/s), so int=10 is ~10% and int=100 ~1%. The launch config
  specified 100; the deck runs 10.
- **Radiation: not a regression source, and #26 is confirmed live.** Measured
  from the current 24-h logs: 359 calls, intervals min 180.0 / median 180.7 /
  max 181.6 s. The fix did raise the call count ~1.4x versus the broken 255 s
  cadence, worth ~1% of wall.

## 41. 30p RECURS on CONUS404 driving -- driving data, date and event are all eliminated

Davies, 2020-12-28, CONUS404-driven, cfl 0.3, `real_width=10`. Died at

    Reading weather data 64800.10547 6 7 9
    Reading surface data 64800.10547 6 7 9
    Coarse STEP 69876 ends. TIME = 64800.10671 DT = 0.8675203919
    Erroneous arithmetic operation

**Same rotation, same model time, everything else different.** Against the ERA5
Jan-9-2023 case:

| | ERA5 Jan-9 2023 | CONUS404 Dec-28 2020 |
|---|---|---|
| death time | 64801.05 / 64801.07 s | **64800.11 s** |
| death step | 46617 / 46632 | **69876** |
| rotation | `6 7 9` | `6 7 9` |
| cfl | 0.2 | 0.3 |
| driving data | ERA5 0.25 deg via erftools | CONUS404 4 km via our converter |
| frame grid | 40 x 22 x 38 | 118 x 71 x 46 |
| event | atmospheric river | winter soaker |
| dt at death | ~1.19 | 0.868 |
| mass over the leg | -- | +0.014% (1.588663e15 -> 1.588892e15) |

The prediction stated before the run held exactly. What this ELIMINATES, none of
which was previously ruled out: the driving dataset, the frame producer
(erftools vs our converter), the frame grid dimensions (nx,ny,nz all changed),
the frame content, the date, the meteorological regime, the timestep size, and
the step count (69,876 vs 46,617 -- so it is not an accumulated-step-count
threshold either).

What it does NOT eliminate, and now stands almost alone: the rotation machinery
in `ERF_WeatherDataInterpolation` and the allocation churn it performs. Note that
with 3-hourly frames "rotation 7" and "18 h" are the SAME event by construction
(`idx1 = t / 10800`, so idx1 = 6 first occurs at t = 64800), so this run cannot
separate them. **Emitting hourly frames would separate them for the first time**
-- rotation 7 would fall at t = 6 h -- and the converter can do that with a
one-line change. That is the cheapest remaining discriminator and it needs no
model change at all.

Combined with 40 (the failure post-dates the 192x96 domain change and 128x64 runs
crossed rotation 7 and survived) the surviving hypothesis is narrow: something in
the rotation path that is sensitive to the 192x96 allocation layout.

Per instruction the run was NOT restart-bracketed. It reached 18 h of the 23 h
window, so hours 0-18 are available for scoring; plotfiles through plt69876.

## 42. NSCBC SURVIVES rotation 7 where Davies dies -- and completes the first unbracketed full day of this campaign

CONUS404-driven, 2020-12-28, cfl 0.3, same binary, same frames from the same files.

| | Davies (real_width=10) | NSCBC (sigma=0.03, keep_ramp=0) |
|---|---|---|
| outcome | **FPE at TIME 64800.11, step 69876** | **completed TIME 82816, step 73172** |
| rotation `6 7 9` | died on it | **crossed it at t=64800.82 and ran on** |
| rotations crossed | 6 | **10** |
| NaN / FPE | 0 / 1 | **0 / 0** |
| mass | +0.014% to 18 h | **+0.044% over 23 h** |

**This is the first continuous 23 h run of the campaign with no restart
bracketing.**

**What it eliminates, and it is the sharpest cut yet.** Both arms call the same
`WeatherDataInterpolation` -> `FillForecastStateMultiFabs`, so in the SURVIVING
run the frame read, the ~13 `Gpu::DeviceVector` allocations churned per rotation,
and the `LinComb` into `forecast_state_interp` all executed identically. **The
allocation-churn hypothesis (30r/41) is therefore weakened, not confirmed** --
that code ran and did not kill NSCBC. Further, per 30q `fill_from_realbdy` is a
NaN DETECTOR, and NSCBC's mass controller performs the same host ordered compares
on reduced values yet never trapped: so **no NaN was produced under NSCBC at all**.

What differs is downstream and is Davies-specific: the specified-zone fill and the
relaxation ramp (`keep_ramp=0`, so NSCBC has no ramp and replaces the wall with
the characteristic solution). The non-finite value is generated by the Davies
fill/relax path consuming the rotated frame, not by the rotation machinery.
Item 40 still stands -- the 128x64 runs that crossed rotation 7 were Davies runs,
so domain size matters too -- but the mechanism now has a much narrower home.

### Preliminary like-for-like at 18 h, full domain vs wrfout d02

**Scoring against d02 measures agreement with another model, not skill.** d02
carries a measured 1.09x wet bias concentrated on orographic peaks and +68%
small-scale variance against MRMS over land (39a). MRMS-over-land is the
observational check.

| mask | arm | bias | PCC | RMSE |
|---|---|---|---|---|
| full domain | Davies | **5.26x** | **-0.016** | 131.74 |
| full domain | NSCBC | 1.46x | **0.590** | 17.44 |
| interior d>=20 | Davies | 0.88x | 0.434 | 8.02 |
| interior d>=20 | NSCBC | 1.16x | 0.290 | 7.98 |
| LAND | Davies | 3.38x | 0.031 | 157.98 |
| LAND | NSCBC | 1.79x | **0.591** | 39.33 |

Davies' full-domain correlation with the reference is **zero** -- the wall band
(252 mm/day at d=0 against d02's 7.7) swamps the field entirely. Its interior is
fine (0.88x); the domain-mean 41.6 mm against d02's 7.9 is all band.

FSS, land, the informative rungs: at 15 mm Davies 0.634->0.876 vs NSCBC
0.687->0.859; at **30 mm Davies 0.267->0.442 vs NSCBC 0.420->0.537** -- NSCBC
clearly better on the heavy tail. Percentile-matched splits: at the low base rate
Davies is marginally ahead (interior r=0.88: 0.877 vs 0.876), at the high rate
NSCBC is far ahead full-domain (r=0.18: 0.452 vs 0.281) and behind in the
interior (r=0.10: 0.168 vs 0.235).

Spectrum vs d02: Davies 0.681, **NSCBC 1.053** -- NSCBC matches the 4 km driver's
small-scale variance almost exactly. Wall band d=0: Davies 252.4, NSCBC 12.8,
d02 7.7. Fetch ratio at the inflow wall: Davies **211.6x**, NSCBC **1.72x**.
Islands (three western): Davies 30.9/33.9/34.7, NSCBC 27.0/22.1/12.6,
d02 8.6/9.7/12.3.

MRMS-over-land at 18 h (secondary, observational): Davies bias 3.94x PCC 0.004;
NSCBC bias 2.09x PCC 0.558. FSS 15 mm: Davies 0.589->0.905, NSCBC 0.612->0.885.
The 1 and 5 mm rungs are degenerate over land as anticipated (base rates 1.000
and 0.993).

### NSCBC full 23 h -- separate table, do not mix windows

vs d02(23h), full domain: bias 1.44x, **PCC 0.616**, RMSE 23.86. FSS 1/5/15/30 mm
at 3 km: 0.848 / 0.740 / 0.505 / 0.477, rising to 0.926 / 0.880 / 0.675 / 0.719
at 60 km. Interior d>=20: bias **1.03x**, PCC 0.269.
vs MRMS over land: bias 1.92x, PCC 0.527; FSS 15 mm 0.728->0.847, 30 mm
0.478->0.789.

**Reading.** On CONUS404 driving NSCBC is the better configuration on nearly every
axis: it runs continuously, holds mass, agrees with both the 4 km driver and MRMS
where Davies does not correlate at all, matches the driver's small-scale variance,
and has effectively no wall band (1.72x vs 211.6x at the inflow wall). Davies
retains a small edge only in the percentile-matched interior at low base rate --
the same "places the envelope, produces no storm" signature as item 37, and the
same caveat applies.

## 43. ROTATION 7 IS EXONERATED. The failure is model time / accumulated history, not a frame index

The confound present in every run of this campaign: with 3-hourly frames
`idx1 = t/10800`, so idx1 first reaches 6 at exactly t = 64800 s. "Rotation 7"
and "18 h" have therefore been the SAME EVENT in every prior run, and the record
has consistently framed it as a count ("the failure is at rotation 7 -- a count,
not a time").

**Separated for the first time.** The converter regenerated the same day at
HOURLY cadence (10 frames, 00-09Z), putting the [6,7] transition at t = 6 h
instead of 18 h. Davies, same binary, same deck, cfl 0.3:

    Reading weather data 21600.56055 6 7 10      <- the exact transition that kills it at 18 h
    ...
    Coarse STEP 25714 ends. TIME = 28800         <- ran to the 8 h stop
    NaN = 0   FPE = 0   rotations crossed = 10   mass -0.163%

**The [6,7] rotation is harmless.** Crossed at t = 6 h and the run continued for
another 2 h and three further rotations ([7,8], [8,9] and the initial read),
completing its window cleanly. Every one of the ten rotations was clean.

**So the failure is NOT a frame-index count.** It is model time near 18 h, or the
accumulated history that goes with it. Every hypothesis built on "what does the
7th frame first touch" -- including the enumeration of frame-indexed allocations
in 30p, and the allocation-churn reading in 30r/41 -- was aimed at the wrong
variable. Item 41's elimination table stands (driving data, date, regime, frame
grid, frame producer, step count are all still excluded); what changes is the
surviving hypothesis, which is no longer "rotation 7 does something special".

Combined with 42 (NSCBC crosses the 18 h event and completes 23 h, with the same
FillForecastStateMultiFabs churn executing in the surviving arm), the two live
readings are now:
  * something in the DAVIES specified-fill / relaxation-ramp path degrades with
    accumulated model time and becomes fatal near 18 h; and
  * it is domain-dependent (item 40: 128x64 completed 24 h continuously).

Note the step counts at failure are NOT constant either: 46,617 (ERA5, cfl 0.2)
vs 69,876 (CONUS404, cfl 0.3), a 1.5x spread at the same model time. So it tracks
MODEL TIME, not step count, not rotation index, not wall-clock.

**The one experiment that would close this**, and it is cheap: run the hourly-frame
Davies configuration out to 18 h (needs 20 hourly frames, ~20 min of conversion,
~75 min of GPU). If it dies at 18 h with idx1 = 17, model time is confirmed as the
variable and the frame index is fully eliminated. If it survives, the difference is
the frame cadence itself and the interval becomes a workaround for the Davies path.

### 42a. Full 0-18 h analysis: the Davies excess is PRESENT FROM HOUR 1, not developed

d02 accumulation contract re-verified: stamps index 0 = `2020-12-28_00:00:00`,
index 18 = `2020-12-28_18:00:00`, last-minus-first, `BUCKET_MM = 100` applied to
517 wrapped cells (naive differencing gives min -99.835 mm, corrected min 0.000).
`RAINC` = 0 over the window. MRMS = 18 past-hour files valid 01Z..18Z. Identical
periods on all four fields.

**Bulk (mm), all four on our grid:**

| mask | field | mean | bias | PCC | RMSE | MAE |
|---|---|---|---|---|---|---|
| full | d02 | 7.90 | -- | -- | -- | -- |
| full | Davies | 41.59 | **5.26x** | **-0.016** | 131.74 | 39.16 |
| full | NSCBC | 11.51 | 1.46x | **0.590** | 17.44 | 8.66 |
| interior | d02 | 6.82 | -- | -- | -- | -- |
| interior | Davies | 5.98 | 0.88x | 0.434 | 8.02 | 5.28 |
| interior | NSCBC | 7.92 | 1.16x | 0.290 | 7.98 | 6.03 |
| LAND | d02 | 17.79 | -- | -- | -- | -- |
| LAND | Davies | 60.15 | 3.38x | 0.031 | 157.98 | 54.71 |
| LAND | NSCBC | 31.86 | 1.79x | **0.591** | 39.33 | 23.76 |
| LAND | MRMS | 15.26 | -- | -- | -- | -- |
| LAND | Davies vs MRMS | 60.15 | 3.94x | **0.004** | 158.74 | 55.31 |
| LAND | NSCBC vs MRMS | 31.86 | 2.09x | 0.558 | 42.22 | 25.41 |

**d02 and MRMS agree closely on this window** -- land means 17.79 vs 15.26, and
their wall-normal profiles nearly overlay -- which is the strongest support yet for
using d02 as the full-domain reference. Davies correlates with NEITHER (0.031 and
0.004 over land): the band destroys the field.

**Wall-normal profile (figure 3, left) is the boundary story in one line.**
Davies 10.0 in at d=0 falling to 0.35 by 45 km and then BELOW both references in
the deep interior (0.08 in at 140 km against d02 0.18, MRMS 0.17) -- band plus
starved interior. NSCBC 0.44 in at the wall, a mild bump to 0.85 at 15-20 km, then
tracking the references at 0.32-0.35. d02 and MRMS overlay each other throughout.

**Hourly rate settles the question the record could not: the divergence is there
from the first hour.** Domain-mean mm/h:

| h | Davies | NSCBC | d02 | MRMS |
|---|---|---|---|---|
| 1 | **5.843** | 0.260 | 0.365 | 0.319 |
| 2 | 6.843 | 0.392 | 0.434 | 0.525 |
| 3 | 7.099 | 0.580 | 0.421 | 0.506 |
| 6 | 2.719 | 1.118 | 0.422 | 0.377 |
| 9 | 1.119 | 0.605 | 0.480 | 0.492 |
| 12 | 0.666 | 0.501 | 0.421 | 0.343 |
| 18 | 0.783 | 0.468 | 0.437 | 0.211 |

Davies is **16-18x the reference in hour 1** and decays through hour 9; NSCBC
tracks d02 and MRMS from the start. Over LAND the point is sharper still: at hours
1-3 d02 and MRMS are both **exactly 0.000** mm/h -- the event had not reached land
-- while Davies produces 3.5, 3.5 and 5.5 mm/h and NSCBC 0.54, 0.55, 0.24. So the
Davies band is not a slow drift that accumulates; it is an initialisation transient
that fires immediately, dwarfs the event for the first three hours, and then decays
into a persistent ~2x land wet bias.

Spectra (figure 4): interior-core PSD ratio in the 8-19 km band, Davies 0.681,
**NSCBC 1.053** against d02 -- NSCBC reproduces the 4 km driver's small-scale
variance almost exactly.

Figures: `c404_18h_fields_0-12in.png`, `c404_18h_fields_0-3in.png`,
`c404_18h_differences.png`, `c404_18h_boundary.png`, `c404_18h_spectra.png`,
`c404_hourly_rate.png`. Full tables (FSS fixed + percentile-matched on three masks,
quantiles, wet-area fractions, fetch bins from both wall pairs, all eight islands,
hourly series): `figs/tables_18h.txt`.

## 44. du_max is NOT the lever: the CONUS404 excess is OROGRAPHIC, not boundary-adjacent

### Pre-run diagnostic (free, on existing output): the excess separates cleanly

du_max=8 NSCBC 23 h vs d02, stratified:

| terrain bin | n | NSCBC | d02 | bias |
|---|---|---|---|---|
| flat <100 m, LAND ONLY | 408 | 14.56 | 19.41 | **0.75 (DRY)** |
| 100-500 m | 1101 | 38.72 | 25.49 | 1.52 |
| >500 m | 1322 | 55.04 | 24.82 | **2.22** |

By distance from the inflow walls: 1.82 (0-6 km), 1.54 (9-27), 0.88 (30-57),
1.10 (60-117), 1.56 (120-285 km). **The cross-tab settles which is which: every
cell above 100 m lies in the 40-95 band, so the far-field 1.56 IS the orography.**
Within flat terrain alone the inflow profile is 1.82 / 1.54 / 0.88 / 1.10 / 1.26.

So the domain-mean 1.44x wet bias is not uniform over-production at all: it is
**mountains too wet (2.22x) and flat land too dry (0.75x)**, averaging to
something that looks uniform and is not. The boundary-adjacent part is 855 cells
(4.6% of the domain); the orographic part is 2423 cells sitting 120-285 km from
any inflow wall -- which a wall-normal velocity clamp cannot plausibly reach.

Prediction registered before the run: du_max=2 trims the near-wall excess, leaves
the orographic excess intact, and costs mass.

### The run: du_max=2, 23 h, everything else identical

Completed TIME=82816, 72236 steps, zero NaN/FPE, 10 rotations.

**Mass: +0.0440% final, peak +0.044%, min -0.199% -- IDENTICAL to du_max=8's
+0.044%.** The Jan-9 railing (+2.2% at du_max=2) does NOT transfer: on this case
the demanded du rarely reaches 2 m/s, so the clamp costs nothing. That half of the
prediction was wrong.

**And it buys nothing.** Fields are 98.4% correlated (r=0.9842); domain mean 14.01
vs 13.78 mm -- du_max=2 is marginally WETTER.

| metric vs d02 | du_max=8 | du_max=2 |
|---|---|---|
| full-domain bias / PCC | 1.44 / 0.616 | **1.47** / 0.626 |
| interior bias / PCC | 1.03 / 0.269 | 1.01 / 0.325 |
| LAND bias / PCC | 1.76 / 0.551 | **1.83** / 0.544 |
| LAND vs MRMS bias | 1.92 | **2.00** |
| bias >500 m | 2.22 | **2.34** |
| bias 0-6 km from inflow | 1.82 | **1.93** |
| spectrum 8-19 km | **0.957** | **0.641** |
| FSS 15 mm 3 km (fixed) | 0.505 | 0.534 |
| PM-FSS 30 mm 3 km | 0.511 | 0.531 |

**The near-wall excess got WORSE, not better** (1.82 -> 1.93 at 0-6 km), so the
prediction failed in direction as well. The only clear degradation is the
spectrum, which falls from 0.957 to 0.641 of d02 -- du_max=2 loses a third of the
small-scale variance that made du_max=8 match the 4 km driver. Against that,
heavy-threshold FSS improves slightly (15 mm 3 km 0.505 -> 0.534).

**Verdict: du_max=8 stays.** It is not carrying the excess, it costs nothing in
mass on this case, and it holds the small-scale variance that du_max=2 loses.

### What the excess actually is

Wall-normal profile (inches): NSCBC peaks at **1.22 at d=5** against d02's 0.41 --
a 3x bump 15 km inboard, not at the wall -- then falls to 0.34 by d=20 while d02
holds 0.42. Santa Ynez point (MRMS max, 34.486N -119.802W): du8 **6.07 in**,
du2 6.28, d02 1.89, MRMS **2.61**. Domain max: du8 **17.47 in**, du2 17.33,
d02 3.03, MRMS 2.61 -- a **6x over-production at the peak**.

So the outstanding NSCBC defect is orographic over-intensification, worst above
500 m and at the domain maximum, with a secondary bump ~15 km inboard of the
walls. Neither is a mass-controller artifact. Figures:
`c404_dumax_fields_0-3in.png`, `c404_dumax_fields_0-12in.png` (matched to
`c404_18h_fields_0-3in.png`).

## 45. BLOCKED: ERF's Smagorinsky2D is unstable in this configuration at every Cs tested

The manifest flagged this as open -- "[OPEN, flagged in audit item 4] no LES
closure at 3 km ... 3 km is grey-zone ... Revisit after step 7" -- and the
remaining defect (44) is orographic, so horizontal subgrid mixing is the right
next lever. It does not run.

**The configuration is correct and ERF confirms it.** `erf.les_type=Smagorinsky2D`
with MYNNEDMF prints exactly the intended pairing:

    Selected a PBL model and an LES model: Using PBL for vertical transport, LES for horizontal
    Turning off mix_isotropic for 2-D Smagorinsky
    Using 2D Smagorinsky LES model at level 0

ERF *requires* the 2-D form with a PBL model (it errors otherwise), which is the
correct grey-zone choice at dx/dz = 150.

**Cs has no default and 0.25 is the reference-matched value.** `Cs = 0` in
ERF_TurbStruct.H and ERF errors "Need to specify Cs for Smagorsinky LES", so the
coefficient must be chosen. Read from the references rather than guessed: **both**
wrfout d02 (1.5 km) and CONUS404 d01 (4 km) run `DIFF_OPT=2, KM_OPT=4` -- WRF's
2-D horizontal Smagorinsky, default `c_s = 0.25`. And ERF's formulation is
identical in form: with `mix_isotropic` forced off, `DeltaH = sqrt(dx*dy)` and
`nu_h = Cs^2 * DeltaH^2 * |S|`, the same as WRF's `KM_OPT=4`. So 0.25 transfers
exactly. At 3 km that is nu_h ~ 5.6e5*|S|, about 560 m2/s at typical strain.

**It dies at every Cs, with time-to-death monotone in Cs:**

| Cs | died at step | model time | first NaN cell |
|---|---|---|---|
| 0.25 | 7 | 2.4 s | (0,4,8) |
| 0.20 | 10 | 4.0 s | (0,95,1) |
| 0.15 | 30 | 14.4 s | (0,20,27) |
| 0.10 | 140 | 145 s | (32,10,7) |
| 0.05 | 150 | 259 s | (0,20,27) |

**Not a diffusion CFL violation.** The explicit horizontal-diffusion limit here is
dt < dx^2/(4*nu_h) = 4/|S| seconds; at dt = 0.44 s that needs |S| > 9 s^-1, a
27,000 m/s velocity difference across a cell. Nothing near that exists.

**Not boundary-scheme specific.** Davies fails too, at step 10 with the first NaN
at (118,95,1) -- the yhi wall -- against NSCBC's step 7 at (0,4,8), the xlo wall.
Four of the five first-NaN locations sit ON a lateral wall (i=0 or j=95).

**Not a tunable stability margin.** Lowering Cs only postpones it: 7 -> 10 -> 30
-> 140 -> 150 steps. Every arm still goes fully non-finite (884,736 cells at
Cs = 0.15 and 0.05). This is the same signature as the intermediate-sigma NSCBC
failures (34): monotone in the coefficient, fatal at every value, which reads as a
defect the coefficient paces rather than a stability limit it crosses.

**Consequence: the horizontal-mixing lever cannot be exercised until this is
fixed.** Shrinking Cs further is not a workaround -- by the trend it would still
fail, and a coefficient far below what both reference models use would not be a
meaningful test of the hypothesis anyway. The orographic over-intensification
(44: 2.22x above 500 m, 6x at the domain max) therefore stands unaddressed, and
the natural next step is a diagnosis of the Smagorinsky2D lateral-boundary path
rather than another parameter sweep.

Not committed to the deck or manifest: `erf.les_type` remains `"None"` and the
manifest MUST line is unchanged, since nothing here validates a new value.

## 46. FALSIFIED: the terrain ghost IS filled. The mechanism claimed here was wrong; the test-coverage gap is real

**This item previously claimed that `z_nd(-1)` is never filled at a non-periodic,
non-symmetry wall, and that horizontal LES is its first consumer. That is wrong.
The instrumented probe run killed it before the "fix" was applied.**

**The measurement.** A temporary reduction was added at the end of
`make_terrain_fitted_coords`, comparing the exterior ghost plane against the first
interior plane of `z_phys_nd`, over the FULL column, on all four lateral walls,
on the production CONUS404 deck (`xlo/xhi/ylo/yhi.type = "Outflow"`,
`use_real_bcs = true`, i.e. exactly the configuration the claim was about):

    [ghost46] FULL COLUMN (k=0..48) terrain z_nd, ghost vs first interior, lev 0:
      xlo  GHOST[0,19003.00977]  INTERIOR[0,19003.00977]  max|g-v| 0
      xhi  GHOST[0,19003.00977]  INTERIOR[0,19003.00977]  max|g-v| 0
      ylo  GHOST[0,19003.00977]  INTERIOR[0,19003.00977]  max|g-v| 0
      yhi  GHOST[0,19003.00977]  INTERIOR[0,19003.00977]  max|g-v| 0

`max|ghost - interior| = 0` at every k on every wall. The ghost already holds
exactly `z_nd(-1) = z_nd(0)` -- bit-for-bit the value the proposed fix would have
written. **The fix was a no-op.**

**What actually fills it.** `ERF_TerrainMetrics.cpp:367-370`, inside
`init_which_terrain_grid`, in the same kernel that computes the max terrain
height:

    int ii = amrex::max(amrex::min(i,imax),imin);
    int jj = amrex::max(amrex::min(j,jmax),jmin);
    // Fill the lateral boundaries
    z_arr(i,j,k0) = z_arr(ii,jj,k0);

An unconditional clamp-to-domain index -- zero-gradient surface terrain into every
lateral ghost, for every BC type, no `phys_bc_type` guard. The 3-D coordinate is
then built from it on grown boxes, which is why the equality holds at all k, not
just at k0.

**Why the original diagnosis was wrong.** It surveyed the two places that look
like ghost fills -- `FillBoundary(periodicity())` and the symmetry block at
`ERF_TerrainMetrics.cpp:98-134` -- and never looked inside
`init_which_terrain_grid`. The `setDomainBndry(Real(1.234e20),...)` poison at
line 58 made the wrong conclusion look confirmed: it reads exactly like a
deliberate tripwire that only the symmetry block ever clears. It is cleared, just
somewhere the survey did not cover.

**The trap this avoided.** With `les_type = "None"` the planned validation would
have returned a near-identity result, and that identity would have been read as
"the fix is correct and inert as predicted". It would have been identity because
the fix changes nothing at all. Instrumenting the ghost before fixing it, rather
than after, is the only reason that did not happen.

**45 still stands -- re-confirmed on the current binary.** Smagorinsky2D,
Cs = 0.25, production CONUS404 deck, clean HEAD build: SIGFPE ("Erroneous
arithmetic operation") at coarse step 16. Per the FPE-trap signature the SIGFPE is
the NaN alarm, not the fault. A 12-step run completes cleanly, which is only
because 12 < 16 -- short runs do not clear this.

**Ruled out by this pass (three candidates, all clean):**
1. *Terrain ghost* -- measured above, exactly zero-gradient.
2. *Strain-rate boundary handling* -- `ERF_ComputeStrain_T.cpp:61-130` handles
   lateral Dirichlet explicitly and includes `ext_dir_ingested`, which is the
   real-BC type this deck uses. It is not BC-blind.
3. *Eddy-viscosity exterior ghost* -- `mu_turb` is computed only on
   `growntilebox(1) & domain`, so it is genuinely never computed outside the
   domain; but `ComputeTurbulentViscosity` sets `impose_phys_bcs = true` whenever
   `les_type != None` and then does `FillBoundary` plus a masked physical-boundary
   extrapolation into those planes. Handled.

**The durable finding: no test in the tree exercises this combination.** Strict
census of `Exec/**/inputs*` that SET a real `erf.les_type`:

| | count |
|---|---|
| inputs setting a real les_type | 52 |
| of those, doubly periodic in x and y | 46 |
| of those, NOT doubly periodic | 6 |
| using `Smagorinsky2D` | **0** |
| using `use_real_bcs` | **0** |

So `Smagorinsky2D` has no test at all, and no LES test anywhere in the tree runs
against ingested real lateral boundaries. In the 46 doubly-periodic cases the
lateral ghost is filled by wrap, so nothing there can expose a lateral-wall LES
defect even in principle. That gap is real regardless of this item's dead
mechanism, and it is the reason a lateral-wall LES failure can sit in the tree
unnoticed.

**Status: the Smagorinsky2D mechanism is NOT localized.** Three candidates
eliminated, no remaining hypothesis from source reading alone. Per the standing
rule for this diagnosis, we stop here rather than continue guessing. The
horizontal-mixing lever stays blocked and the orographic over-intensification
(44) passes to the microphysics path (WSM6) instead.

## 47. TERRAIN EXONERATED: ours is comparable in height and 8-10% SMOOTHER than d02 over the over-producing cells

No runs. All three terrains compared on our 192x96 grid, each as its model sees it.

**What each model's smoothing actually is:**
- **ERF, `terrain_smoothing = 1` (STF).** This is a VERTICAL coordinate method --
  Smoothed Terrain Following -- that decays the terrain-following deformation with
  height (ERF_TerrainMetrics.cpp:196+, `case 1`). It does **not** horizontally
  filter the surface: the k0 plane is the raw terrain file. So our surface
  elevation is exactly what `dem_to_erf_terrain.py` wrote.
- **`dem_to_erf_terrain.py` point-samples.** `rows = np.round(rows).astype(int)`
  -- nearest-neighbour sampling of the 30 m GLO-30 DEM at 3 km spacing, with no
  area averaging.
- **d02 (1.5 km):** `SMOOTH_OPTION = 0`, i.e. no extra WPS topo smoothing, but
  geogrid's topo interpolation area-averages the source DEM onto each 1.5 km cell.
- **CONUS404 d01 (4 km):** same geogrid approach at 4 km, and it measures as the
  smoothest of the three.

**Elevation on our grid (LAND cells, m):**

| field | mean | p50 | p90 | p99 | max |
|---|---|---|---|---|---|
| ERF 3 km (GLO-30 point-sampled) | 592.5 | 452.8 | 1277.9 | 1915.0 | 2355.2 |
| d02 1.5 km -> 3 km (area-avg) | 594.2 | 448.5 | 1291.7 | 1974.8 | 2647.1 |
| CONUS404 4 km -> 3 km | 592.9 | 451.1 | 1276.5 | 1908.1 | 2290.0 |

Essentially identical -- domain-mean land elevation agrees to 0.3%.

**Slope |grad h| (m per km), LAND:**

| field | mean | p50 | p90 | p99 | max |
|---|---|---|---|---|---|
| ERF 3 km | 37.25 | 27.83 | 81.02 | 136.87 | 200.33 |
| d02 -> 3 km | **40.69** | 30.11 | **91.48** | 155.33 | **214.25** |
| CONUS404 -> 3 km | 31.00 | 23.44 | 67.36 | 113.08 | 156.46 |

**Over the over-producing cells (our terrain >= 500 m, n = 1322):**

| field | elev mean | slope mean | slope p90 | slope max |
|---|---|---|---|---|
| ERF 3 km | 1019.5 | 49.85 | 102.81 | 200.33 |
| d02 -> 3 km | 1025.7 | **55.12** | **111.73** | **214.25** |
| CONUS404 -> 3 km | 1015.7 | 41.72 | 83.15 | 156.46 |

**Ratios ERF/d02 over those cells: elevation 0.99, mean slope 0.90, p90 slope
0.92.** Santa Ynez point: ERF 402 m at 94.1 m/km, d02 331 m at 109.5 m/km.

**Verdict: terrain is exonerated.** Over exactly the cells carrying the 2.22x wet
bias, our terrain is the same height as d02's and 8-10% LESS steep -- and d02
produces 2.6x LESS precipitation there. Steeper forced ascent cannot explain the
over-production; if terrain slope were the driver, d02 should over-produce more
than we do, and it does the opposite. No regridding of d02's terrain and no
smoothing change is indicated.

**A prediction of mine that the data falsified.** Having found the point-sampling
in `dem_to_erf_terrain.py` I expected aliasing to make our 3 km terrain steeper
than an area-averaged field. It does not: point-sampling a 30 m DEM at 3 km
spacing yields a field 10% smoother in slope than d02's 1.5 km area-average
re-averaged to 3 km, because the 1.5 km field retains resolved structure that a
single 3 km point sample never captures. The point-sampling remains a real
methodological wart worth fixing for its own sake, but it is not making our
terrain rough, and it is not the source of the orographic bias.

**Consequence: the LES / microphysics path stands as the remaining explanation**
for the orographic over-intensification, and the horizontal-mixing lever is still
blocked behind item 46.

## 48. Smagorinsky closed on evidence: 3-D is REFUSED with a PBL model, both formulations fail above a common horizontal-viscosity threshold, and the NaN is real

Follow-on to 45/46. All runs on the production CONUS404 deck, clean HEAD build,
`max_step = 500`, plot/check off.

**1. ERF REFUSES 3-D Smagorinsky with any PBL model. It does not warn, and it
cannot silently double-count.** `ERF_TurbStruct.H:115-118`, inside
`if (pbl_type != PBLType::None)`:

    if (les_type == LESType::Smagorinsky) {
      if (!smag2d)
        amrex::Error("If using Smagorinsky with a PBL model, the 2-D "
                     "formulation should be used");
    }

`Smagorinsky2D` maps to `les_type = Smagorinsky, smag2d = true` at line 92-95, so
the guard fires on exactly 3-D + PBL. Confirmed at runtime: the arm aborts at
step 0 with that message. **So the double-counting concern is answered by
construction -- the combination is unreachable, not silently permitted.**

**2. The MetGrid / WPS_Test precedent does not transfer.** Neither
`RegTests/MetGrid/inputs_metgrid` nor `RegTests/WPS_Test/inputs_wps` sets
`erf.pbl_type` at all, so both run with `pbl_type = None`, which is the only
reason they clear the guard above. They are also at Cs = 0.1 and **Cs = 0.005**
respectively -- the second is 50x below the reference-model value, which is itself
a hint at the threshold below.

**3. 3-D is NOT immune -- it fails too, at higher Cs. The defect is not specific
to the 2-D path.**

| LES | PBL | Cs | last step | model time | outcome |
|---|---|---|---|---|---|
| 2-D | None | 0.002 | 500 | 563.6 s | clean |
| 2-D | None | 0.01 | 500 | 563.4 s | clean |
| 2-D | None | 0.04 | 151 | 214.7 s | SIGFPE |
| 2-D | None | 0.10 | 19 | 8.8 s | SIGFPE |
| 3-D | None | 0.10 | 500 | 552.6 s | clean |
| 3-D | None | 0.25 | 24 | 15.3 s | SIGFPE |
| 2-D | MYNNEDMF | 0.10 | 246 | 192.0 s | SIGFPE |
| 2-D | MYNNEDMF | 0.25 | 16 | -- | SIGFPE |
| 3-D | MYNNEDMF | any | 0 | -- | REFUSED |

**4. It is a threshold in horizontal eddy viscosity, common to both
formulations.** `nu_h = Cs^2 * DeltaH^2 * |S|`, and DeltaH differs by
construction: 2-D uses `sqrt(dx*dy) = 3000 m`, 3-D uses `cbrt(dx*dy*dz) ~ 500 m`
near the surface on this stretched grid. That is a ~36x ratio in `DeltaH^2` at
equal Cs. Normalizing, the 2-D stable/unstable window (0.01 lives, 0.04 dies) and
the 3-D window (0.10 lives, 0.25 dies) **overlap**. The apparent "2-D is broken,
3-D is fine" is a Cs-normalization artifact, not a code-path difference.

This revises 45's "not a tunable stability margin". It IS a threshold; 45 could
not see it because it swept Cs within one formulation only, where DeltaH is
constant and the threshold looks like "fails at every value".

**5. The PBL only paces it.** 2-D at Cs = 0.10 dies at step 19 without a PBL model
and step 246 with MYNNEDMF. MYNNEDMF's vertical mixing delays the blowup by ~13x
in steps but does not prevent it, and its presence is not required to produce it.

**6. The NaN is real, not just the FPE alarm.** Re-run of the 2-D Cs = 0.10
no-PBL arm with `amrex.fpe_trap_invalid = 0` and `check_for_nans_int = 1`: same
death step (19), and the NaN checker reports **262,704 NaN cells in every state
component (comps 0-14)**. The field is already globally non-finite when detected.
The FPE trap is reporting a genuine blowup, not tripping on a benign operation.

**Consequence: there is no usable horizontal-mixing lever, and the reason is now
evidential rather than a blocked mechanism.**
- The reference-matched value (Cs = 0.25, 2-D, as d02 and CONUS404 both use for
  `KM_OPT=4`) sits ~25x above the 2-D stability threshold in Cs, i.e. ~600x in
  viscosity.
- The largest stable 2-D value is between Cs 0.01 and 0.04. At 0.01 that is
  ~1/600 of the reference viscosity -- physically negligible against a 2.22x
  orographic bias, so running it would not test the hypothesis.
- 3-D at a stable Cs = 0.1 is reachable ONLY by deleting MYNNEDMF. That replaces
  the vertical closure rather than adding a horizontal one, so it is a different
  model, not a lever, and it would confound the orographic comparison outright.
- The "clean" runs are 500 steps ~ 9.4 model minutes. A 23-h arm is ~80,000
  steps. These results bound instability, they do not establish stability.

**Not committed to the deck or manifest:** `erf.les_type` stays `"None"`.

**The orographic over-intensification (44) therefore passes to microphysics
(WSM6), with terrain exonerated (47), du_max exonerated (44), and horizontal
mixing closed here.**
