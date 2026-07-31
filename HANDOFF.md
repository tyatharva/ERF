# ChannelIslands ERF Hindcast — Handoff

**Written 2026-07-30.** Assumes you have none of the prior session's context.
Everything here is reconstructable from `Exec/CanonicalTests/ChannelIslands/UPSTREAM_ISSUES.md`,
which is the campaign's numbered evidence log (items 1–51). Where this document
says "item N", read that entry for the full evidence.

---

## 1. What this is

A convection-permitting hindcast of the Channel Islands / Santa Barbara region
using [ERF](https://github.com/erf-model/ERF), driven by reanalysis lateral
boundary conditions, scored against WRF `wrfout` d02 and MRMS observations.

Two cases exist:

| case | driver | window | status |
|---|---|---|---|
| 2023-01-09 | ERA5 (25 km) | 24 h | scored, closed (item 37) |
| **2020-12-28** | **CONUS404 (4 km)** | **23 h (00–23Z)** | **active** |

All execution is inside the `erf-hindcast` Docker image with the repo bind-mounted:

```bash
docker run --rm --gpus all -v ~/ERF:/app/ERF \
    -v ~/.cdsapirc:/root/.cdsapirc:ro -w /app/ERF erf-hindcast <cmd>
```

Build (single precision, CUDA):

```bash
docker run --rm --gpus all -v ~/ERF:/app/ERF -w /app/ERF erf-hindcast bash -c \
  'git config --global --add safe.directory "*"; ERF_HOME=/app/ERF Build/cmake_single_precision_cuda.sh'
```

Hardware assumption throughout: **one 16 GB GPU (RTX 4080)**. Single precision
everywhere. Do not switch to double precision without asking — memory will not fit.

---

## 2. Current best configuration

**NSCBC σ = 0.03, `du_max` = 8, Morrison microphysics, `les_type = "None"`.**

Run directory `run_c404_nsc/`, final plotfile `plt73172` (t = 82816 s = 23.00 h).
Launched as:

```
erf.cfl=0.3 erf.moisture_model=Morrison erf.les_type=None
erf.hindcast_mass_du_max=8 erf.nscbc_lateral=1 erf.nscbc_outflow=1
erf.nscbc_parts=31 erf.nscbc_sigma=0.03 erf.nscbc_mass_tau=60 erf.nscbc_keep_ramp=0
```

Note `cfl=0.3` is a **run-time override**; the deck (`inputs_c404`) ships `cfl=0.2`.
That was deliberate — the manifest's bar for changing the deck value (reproducing
on the Hurricane Hilary case, more than once) was never met.

### What it scores (23 h accumulation)

| mask | ref | bias | PCC | RMSE |
|---|---|---|---|---|
| full domain | d02 | 1.44× | 0.616 | 23.86 |
| interior d≥20 | d02 | 1.03× | 0.269 | 8.95 |
| LAND | d02 | 1.76× | 0.551 | 56.92 |
| full domain | MRMS | 1.81× | 0.576 | 25.13 |
| interior d≥20 | MRMS | 1.42× | 0.070 | 10.36 |
| LAND | MRMS | 1.92× | 0.527 | 59.25 |

| | value | d02 | MRMS |
|---|---|---|---|
| Santa Ynez point (34.486 N, −119.802 W) | **6.07 in** | 1.89 | 2.61 |
| domain max | **17.47 in** | 3.03 | 2.61 |
| spectrum ratio 8–19 km | — | **0.957** | 1.351 |

Terrain-stratified bias vs d02: flat (<100 m) **0.74×**, 100–500 m **1.48×**,
>500 m **2.21×**. This orographic excess is the campaign's remaining defect.

PM-FSS (full domain, vs d02), 3/9/15/30/60 km at r = 0.87:
0.867 / 0.889 / 0.900 / 0.919 / 0.941.

### It is not the best *placement* arm

Davies scores better on placement (item 34: "no NSCBC arm comes within 0.17 of
Davies at any base rate or scale"). NSCBC is in use because **Davies dies at 18 h**
(see §5) and NSCBC survives the full day. That is a survival choice, not a skill
choice.

---

## 3. THE OPEN QUESTION — item 51

**The orographic excess is a boundary spin-down failure, and it is NSCBC-specific.**

The 2.21× bias is not a steady multiplier. It is accumulated entirely between
**h13 and h17** (item 50). During the storm's main period (h7–h11) the model is
essentially unbiased over high terrain (0.82×, 0.91×, then 1.5–1.65×). Both
references spin down after h12; ERF does not.

Hourly water-vapour flux (integrated 0–6000 m ASL) through the lateral faces,
ERF ÷ the CONUS404 driving frames sampled on the **same faces at the same heights**:

| h | NSCBC σ=0.03 | Davies | driver kg/s |
|---|---|---|---|
| 3 | 1.021 | 1.018 | 6.90e7 |
| 6 | 1.083 | 0.944 | 5.99e7 |
| 9 | 1.486 | **1.029** | 3.72e7 |
| 12 | 1.853 | **1.040** | 2.46e7 |
| 15 | **2.125** | **1.109** | 1.62e7 |
| 18 | 3.885 | 2.613¹ | 3.32e6 |

¹ Davies was dying at h18; through h15 it is clean.

The driver's inflow decays **21×** from h3 to h18. NSCBC decays only **5.5×**.

**The divergence is entirely the ylo (south) face.** xlo — the main inflow — is
faithful in NSCBC at every frame hour (ratio 0.92–0.96). But at h15 the driver
exports −4.24e7 southward while NSCBC manages only −2.06e7, and at h9 the driver
exports 2.02e7 while NSCBC *imports* 8.2e5. **NSCBC will not let moisture leave
through the southern boundary.**

Column water vapour over the >500 m cells is ~1.0× d02 for both arms, so the
excess is neither arriving as extra moisture nor being retained. It is being
*processed faster*, driven by a sub-1500 m wind over the ranges that runs
**1.5–2.0× d02** through exactly the h13–h17 window.

### Neither scheme is a free win — keep this framing

**Davies has the faithful boundary flux (1.02–1.11) and a 0.3–0.6× too-weak
interior response over the ranges. NSCBC has the unfaithful flux and roughly the
right orographic response. They fail at different points in the same chain.**
Do not read "Davies tracks the driver" as "use Davies" — its low-level wind over
the ranges is 3.3–7.4 m/s against d02's 11–13, and too westerly (279–291° vs 242°).

### The next test — NOT YET RUN

**NSCBC σ = 1, 23 h, Morrison, du_max = 8.** σ is the coupling to the far field
(σ=1 fully driver-driven, σ=0 fully interior-extrapolated); 0.03 is nearly
free-running, which is why the boundary does not feel the driver's decay. The
flux says we need *more* coupling. From the item-34 sweep:

| σ | stability | PM-FSS |
|---|---|---|
| 0 | stable 24 h | 0.586 / 0.478 / 0.589 |
| 0.03 | stable 24 h | 0.628 / 0.547 / 0.634 |
| 0.1 | **FPE at 7.9 h** | — |
| 0.3 | **FPE at 1.9 h** | — |
| 1 | stable 24 h | 0.649 / 0.604 / 0.684 |

σ = 1 is the only stable strong-coupling value; everything that would interpolate
is unstable.

**A σ=1 run was started on 2026-07-30 and deliberately stopped by the operator
before completion.** Its partial output is **preserved on purpose** in
`run_c404_sig1/` — do not clean that directory. It has not been scored and
nothing was concluded from it, but its latest checkpoint is a valid restart point
for the σ=1 test. Select the checkpoint with `pick_restart_chk.sh` (see §5.2 —
never restart from a checkpoint written on a truncation step).

**The primary result to report from that run is the flux comparison, not the
accumulation**: hourly ERF/driver inflow ratio at both faces, and specifically
whether ylo now exports when the driver exports.

**Watch the tradeoff:** σ = 1 was 2.63× wet on the Jan 9 ERA5 case. Stronger
coupling should fix the spin-down failure but may reintroduce over-retention. If
both appear, σ cannot satisfy both and the outflow condition needs reformulating
rather than tuning.

---

## 4. What is exhausted, and why each was closed

Do not re-open these without new evidence. Each has an item number with data.

### Boundary formulations
| approach | item | why closed |
|---|---|---|
| Davies relaxation band | 36 | The ramp's ∇F is *simultaneously* the anchoring mechanism and the artifact source — they cannot be separated |
| Band width sweep (`real_width` 1–7) | 36b/36e | Width is a cost knob only; pathology monotone, placement a step then flat |
| NSCBC σ sweep | 34 | Placement never reaches Davies; intermediate σ unstable at both 0.1 and 0.3 |
| NSCBC outflow variants (latch vs Riemann) | 33 | Latch buys survival, Riemann buys the outflow excess, neither fixes placement |
| Helmholtz projection of the band correction | 36 (§"kill condition") | Kill condition FIRES — it reconstructs the rotational component, which *is* the artifact |
| Tangential-only band relaxation | 36 (§"tangential-only") | Fails; explains why no *local* fix exists |
| Spectral nudging | 36c | Declined on evidence — Davies' advantage inverts above 200 km, so the premise fails |

**Considered and declined without implementation — NOT tested:**

- **Mesinger / Eta boundary scheme.** Scoped early in the campaign and declined
  on three grounds, none of which have changed:
  1. The published specification is **five sentences of prose** in Leps et al.
     (2019), with no equations and no stencils.
  2. The primary sources — Mesinger (1977) and Black (1988) — **are not
     digitized**, so the specification cannot be recovered from them.
  3. Its variable count was inherited from **Sundström's hydrostatic
     derivation**. The characteristic analysis for our configuration showed it
     would leave ERF **over-specified by 16 conditions at outflow with Morrison
     active** — i.e. the scheme as published does not carry over to a
     non-hydrostatic compressible solver with this many prognostic moisture
     variables.

  Record this as a scoping decision, not a negative result. Nothing was built and
  nothing was measured. If someone recovers the primary sources or redoes the
  characteristic count for a reduced moisture set, the decision is revisitable.

### Physics levers for the orographic excess
| lever | item | why closed |
|---|---|---|
| Terrain | **47** | Ours is *comparable in height and 8–10% smoother* than d02 over the over-producing cells (mean slope 49.85 vs 55.12 m/km at equal elevation). d02 produces 2.6× less rain there on steeper terrain. Exonerated. |
| `du_max` (mass-constraint clamp) | **44** | Prediction failed both ways: mass unaffected (+0.044%, identical) and the near-wall excess got *worse* (1.82 → 1.93). The excess is orographic, not boundary-adjacent. |
| Horizontal mixing (Smagorinsky) | **48** | A **viscosity threshold**, not a configuration guard. Both 2-D and 3-D fail above a common `Cs²·DeltaH²`. Reference Cs=0.25 is ~600× above threshold in viscosity; the largest stable value is physically negligible. 3-D + any PBL scheme is hard-errored by ERF (filter-width reason, not double-counting — the PBL *overwrites* `Mom_v` over the full column). |
| Microphysics swap (WSM6) | **49** | Worse on nearly every measure. Doubles the domain max (17.47 → 34.67 in), worsens every terrain bin, nearly doubles LAND RMSE. PM-FSS flat, so it changes amount not placement. |
| Spin-up transient | **50** | Excluding h1–4 moves the >500 m bias only 2.21× → 2.13×. Genuine and steady. |

### Ruled out by the item-51 discriminator — do NOT run
The vertical branch (`pbl_type = MYNN25`, a `w_damping` sweep, a `pbl_type=None`
+ 3-D Smagorinsky diagnostic) was scoped but is **not indicated**. The
discriminator returned "boundary", not "internal". Running them would target a
mechanism the evidence excludes.

---

## 5. Known unfixed defects

### 5.1 The 18 h Davies death (item 30p / 41 / 42 / 43)
Davies dies at ~18 h. Established:
- **It tracks model time / accumulated history, not frame index.** Rotation 7 is
  exonerated (item 43): with hourly frames the [6,7] transition happens at t=6 h
  and the run continues past it to the 8 h stop cleanly.
- It **recurs on CONUS404** driving (item 41), eliminating driving data, date,
  regime, frame grid and step count.
- **NSCBC survives it** (item 42) — the first unbracketed full day of the campaign.
- Not a race: `amrex.max_gpu_streams=4` sanitizer run showed memcheck 0 errors
  but full-field NaN.

This is the single blocker preventing use of the better-placing scheme.

### 5.2 The float32-absolute-time class (item 32) — ONE LIVE INSTANCE
An absolute epoch instant (~1.673e9 s) held or differenced in single precision,
where the float32 ULP is **128 s**. Rediscovered four times from four symptoms.

| # | where | status |
|---|---|---|
| 28 | `fill_from_realbdy` interpolation | FIXED |
| 26 | `Radiation::set_grids` cadence | FIXED (`set_model_time()`) |
| — | `ERF.cpp:626` evolve stop condition | FIXED |
| **31** | **`ERF::stop_time` itself (`ERF.cpp:45`)** | **LIVE** |

`stop_time` is `amrex::Real` holding `getEpochTime(stop_datetime)`. The stop
instant a deck asks for is not the one the run uses: 18:00:00 → 64768 s (−32 s),
21:00:00 → 75648 s (+48 s), 00:00:00 → exact. Every leg-1 of every bracketed run
stopped at 64768 and that was misread as "before the fatal rotation, by design."

Related: **item 31** — a checkpoint written on the `stop_datetime` truncation step
POISONS every restart from it. `pick_restart_chk.sh` guards against this by
reading the checkpoint Header (line 7 istep, **line 8 dt**, line 9 time) and
rejecting checkpoints with an insane dt. A guard in `ERF_ComputeTimestep.cpp`
aborts on a restart dt >1000× below the CFL estimate.

### 5.3 Point-sampled terrain
`dem_to_erf_terrain.py` (~line 148) does `rows = np.round(rows).astype(int)` —
nearest-neighbour point sampling of the 30 m GLO-30 DEM at 3 km, with no area
averaging (unlike WRF geogrid). **This is a real methodological wart but it is NOT
the source of the orographic bias** — item 47 measured the resulting field as
*smoother* than d02's, falsifying the aliasing prediction. Worth fixing on its own
merits; do not expect it to change scores.

Also note `erf.terrain_smoothing = 1` (STF) is a **vertical** coordinate method
(`ERF_TerrainMetrics.cpp:196`, case 1). It does not horizontally filter the
surface. The k=0 plane is the raw terrain file.

---

## 6. Verification lessons — these were each learned the hard way

1. **Validate every instrument against a known-nonzero control before trusting
   it.** `pick_restart_chk.sh` read dt from Header line 7 (istep) instead of line 8
   and silently passed poisoned checkpoints; it was caught only by running it
   against a known-poisoned control.
2. **Verify edits by `grep`, not by build exit code.** A build once reported
   `exit=0` on failure because the output was piped to `tail`, which swallowed the
   status — and the subsequent knob-string check passed against the *stale*
   binary. Check artifact mtime after every build.
3. **Confirm knobs do what their names say.** Item 46 is the cautionary tale: a
   complete, plausible source-level mechanism ("the terrain ghost is never
   filled") was **falsified by a 20-line instrumented probe** that measured
   `max|ghost − interior| = 0`. The fix would have been a bit-exact no-op, and the
   planned validation would have returned near-identity — which would have been
   read as confirmation. **Instrument before fixing, not after.**
4. **`gpu_preflight.sh` before every run.** Checks for stray containers and GPU
   memory. Two 24-h jobs on one card corrupts timing — score timing only from
   solo runs.
5. **Filter bash output at the source.** Logs contain binary backtrace bytes;
   plain `grep` silently treats such files as binary and prints nothing. Use
   `grep -a`. A step count was mis-read as 0 this way.
   Also: a filter of `/^==/d` will eat your own `=== ` section headers.
6. **The AMReX submodule must be CLEAN for any scored run.** Five diagnostic
   patches exist, all re-appliable from hash-gated scripts. See the header of
   `validated_config.txt`.
7. **A SIGFPE ("Erroneous arithmetic operation") is the NaN alarm, not the
   fault.** Confirmed: with `amrex.fpe_trap_invalid=0` the same run dies at the
   same step with 262,704 NaN cells in every state component.

---

## 7. Where things live

```
Exec/CanonicalTests/ChannelIslands/
├── UPSTREAM_ISSUES.md      # THE evidence log, items 1–51. Read this first.
├── RUNBOOK.md              # build + run procedures, failure modes
├── validated_config.txt    # canonical config; MUST/EXCEPT/DOMAIN lines
├── inputs_hindcast         # ERA5 (Jan 9) deck
├── stage_run.sh            # preflight: asserts validated_config.txt vs staged deck
├── gpu_preflight.sh        # run before EVERY run
├── pick_restart_chk.sh     # safe checkpoint selection (dt sanity)
├── conus404_to_bin.py      # CONUS404 → ERF .bin frames
├── dem_to_erf_terrain.py   # GLO-30 → terrain file
├── scoring/                # 27 scripts, see scoring/README.md
└── harness/                # 9 run-driver scripts
```

Run directories live in the repo root (`run_c404_nsc/`, `run_c404_wsm6/`, …) and
are gitignored via the `plt*` / `chk*` patterns.

### The manifest guard
`validated_config.txt` is the single source of truth for deck knobs. `stage_run.sh`'s
preflight asserts every `MUST` line against the staged deck and **refuses to stage
a mismatch**. `EXCEPT` marks a deliberate, documented divergence. Adding a knob
makes it mandatory; there is no third option. Every knob carries a `prov:`
provenance line — what it was set for, when, and against which configuration.
`[BROKEN-IC]` marks knobs set before the `frame_from_T` thermodynamic correction
(commit `617017b0`) and never re-validated — treat those as open questions, not
settled configuration.

### Key scoring scripts
| script | what it does |
|---|---|
| `score_c404.py` | main battery: bias/PCC/RMSE per mask, FSS ladder (1/5/15/30 mm × 3/9/15/30/60 km), PM-FSS, spectrum 8–19 km, wall-band profile, fetch-binned bias, islands |
| `terrain_windward.py` | terrain-stratified bias, Santa Ynez / domain max in inches, windward-lee split (flow derived from each arm's own sub-1500 m wind), column snow+graupel |
| `spinup_split.py` | re-accumulate excluding spin-up hours; hour-by-hour and cumulative >500 m bias |
| `inflow_flux.py` | **the item-51 discriminator**: hourly lateral vapour flux per face for ERF arms *and* the driving frames on the same faces; low-level wind over the ranges vs d02; column vapour |
| `hourly_series_arms.py` | hourly domain-mean rate for arbitrary arms + d02 + MRMS |

Reference arrays (23 h, on the 192×96 grid):
`wrf_d02_on_grid.npy` (d02) and `mrms_20201228_on_grid.npy` (MRMS). Their
hour-by-hour reconstruction reproduces these to 0.00% — `spinup_split.py` asserts
this as a control.

### Item numbering
`UPSTREAM_ISSUES.md` items are append-only and referenced throughout. Sub-letters
(29i, 29j, 29m, 30p, 36b) mark follow-ups within one investigation. When an item
is later falsified, it is **rewritten in place with the falsification** (see item
46) rather than deleted — the wrong conclusion and its refutation both stay on the
record.

---

## 7b. Running on a rented pod (RunPod)

**Image:** `ghcr.io/tyatharva/erf-hindcast:cuda126-sm89` -- environment only, no
ERF source (ERF is cloned by the bootstrap and bind-mounted/checked out at
`/app/ERF`). The pod has no Docker; RunPod starts the container FROM this image
and the scripts run natively inside it.

**Hardware: sm_89 only** -- RTX 4090 / 4080 / L40S. See the Blackwell note below.

```
pod/pod_transfer.sh pack            # on the workstation -> 79 MB tarball
scp .../erf_pod_transfer.tar.gz pod:/app/ERF/
pod/pod_transfer.sh verify /app/ERF # on the pod, checksums
pod/pod_bootstrap.sh                # clone, patch, build, verify, VALIDATE
```

`pod_bootstrap.sh` refuses to continue past any failure and ends in a **hard
validation gate** (`pod/pod_validate.py`): a 20-step deterministic run compared
against a reference measured on sm_89, with physical (not bitwise) tolerances.
Cross-architecture bitwise equality is not expected; the gate exists to catch
gross corruption -- wrong arch, bad codegen. If it fails, no production arm runs.

It also verifies more than the CMake cache. On the Blackwell attempt the cache
read back `CMAKE_CUDA_ARCHITECTURES=120` while AMReX actually emitted
`compute_60..86`, so the script additionally asserts `AMREX_CUDA_ARCHS` contains
89 and that `AMREX_USE_FLOAT` is really defined.

**Four GPUs = four independent arms, not one job across four.** The standing
constraint holds: one 23-h job per card, and score timing only from solo runs.

### THE BLACKWELL (RTX 5090 / sm_120) PORT IS BLOCKED -- do not retry casually

Attempted 2026-07-30 and dropped. `Dockerfile.blackwell`,
`apply_kokkos_blackwell_patch.sh` and `Build/cmake_single_precision_cuda_blackwell.sh`
are kept **as documentation, clearly marked non-working**. The image builds and
is a valid environment (CUDA 12.8.1, nvcc 12.8.93 does emit `compute_120`);
**ERF does not compile for sm_120** against the pinned submodules. Three
independent pre-Blackwell components block it:

1. **Kokkos 4.5 `cmake/kokkos_arch.cmake`** -- CUDA arch list ends at
   `HOPPER90`/sm_90; no Blackwell option. With none set, Kokkos auto-detects
   from the BUILD machine's GPU, so a wrong-arch Kokkos would be produced
   silently. (`apply_kokkos_blackwell_patch.sh` fixes this one -- necessary, not
   sufficient.)
2. **Kokkos 4.5 `core/src/impl/Kokkos_NvidiaGpuArchitectures.hpp`** -- the
   `KOKKOS_ARCH_* -> KOKKOS_IMPL_ARCH_NVIDIA_GPU` chain also ends at
   `HOPPER90`, then `#error NVIDIA GPU arch not recognized`. Configure succeeds;
   every Kokkos translation unit then fails to compile.
3. **AMReX routes `AMReX_CUDA_ARCH` through `cuda_select_nvcc_arch_flags`** from
   CMake's **deprecated** `FindCUDA/select_compute_arch`. CMake 4.4.1's copy
   knows 9.0 and 10.0 but **not 12.0**, so it does not error -- it silently falls
   back to a "common" list. Measured: cache said 120/12.0, nvcc got
   `compute_60,61,70,75,80,86`, `AMREX_CUDA_ARCHS=60;61;70;75;80;86`.

**All three are the same silent-failure class this campaign keeps hitting: the
knob reads back correct and the wrong thing happens underneath.** Only #2
announced itself. #3 was caught only by reading the emitted nvcc flags instead
of trusting the cache -- the same discipline as item 46.

Patching through would require patching **AMReX**, which `validated_config.txt`
requires to be clean for any scored run. **The honest route is bumping AMReX and
Kokkos to versions with native Blackwell support, and that is a re-validation
project, not a config change.**

## 7c. BLAST RADIUS — read this before quoting any number in this document

**Item 60: the surface was never ingested.** `t_surf` is pinned at 271.00 K with
zero variance, land and ocean identical, for the whole run, in every
CONUS404-driven result. `Qstar ≈ 0` — essentially no surface moisture flux
anywhere. The reference d02 has a real SST field (mean 286.55 K, std 1.632).
**ERF ran 10.7 K cold with a dead surface; the references did not.**

Also **item 58**: the CONUS404 wind rotation was never applied, so every frame
carried vectors mis-oriented by 8.5–22°.

This partitions every result in this campaign:

### SURVIVES — ERF-vs-ERF comparisons (same broken surface in every arm)

Rankings hold; **magnitudes do not**.

- Davies vs NSCBC (items 30–36, 42) — the 18 h Davies death, NSCBC's survival
- the σ sweep (item 34) and σ=0.1 stability (item 53)
- the band-width sweep (36b/36e), Helmholtz projection, tangential-only (36)
- `du_max` (44), WSM6 vs Morrison (49), MYNN25 vs MYNNEDMF (54)
- item 57's **relative** finding that the excess follows the wall rather than
  the terrain — both bands share the same surface

### SUSPECT — absolute bias against d02/CONUS404

Every one of these compares ERF to a reference that had a working surface.

- the orographic excess (2.21×, 2.12–2.61×) — items 44, 47–57
- Domain B's marine over-production (6.9–9.2×) — item 57
- the interior ascent excess (2.41× d02) — item 59
- domain max 17.47 in / 18–21 in, Santa Ynez 6.07 in — HANDOFF §2, items 49, 57
- the full §2 scoring table (bias, RMSE; PCC and PM-FSS are placement measures
  and degrade more gracefully, but are not clean either)
- **item 59's velocity table (ERF vs driver) — SUSPECT.** An earlier assignment
  of it to the surviving bucket was wrong: it compares ERF to the driver, not
  ERF to ERF.

**Nothing absolute should be quoted until a fixed-surface run exists.**

## 8. Standing constraints

- **No upstream activity.** Do not file issues, open PRs, or comment on
  `erf-model/ERF`. Everything goes in the local `UPSTREAM_ISSUES.md`. (#3487 and
  #3491 are already out; leave them, nothing further.)
- **CDS credentials are never copied into the image.** Always mount read-only:
  `-v ~/.cdsapirc:/root/.cdsapirc:ro`. Never ask the user to paste the key.
- **16 GB GPU. Single precision. No double-precision builds without asking.**
- **Never run two 24-h jobs on one card**, and score timing only from solo runs.
- **No new mechanism hunt without checking with the operator first.**
- Diagnostics stay stripped from committed source. Add, measure, revert, rebuild.
