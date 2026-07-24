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
