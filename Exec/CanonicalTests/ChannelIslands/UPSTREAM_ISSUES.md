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
