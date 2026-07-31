# Predictions for the 71 h NSCBC-vs-Davies comparison

Written **before** either arm was launched, at commit `dc4a29fe` (the item-60
surface fix). Recorded so the result can falsify something rather than be
explained after the fact. Every number below is from the campaign record; where
the record is a single measurement rather than a fit, that is said.

**Configuration.** Domain A, pin 35.4, 192x96x48 at 3 km. Start 2020-12-26 00Z,
stop 2020-12-29 00Z (71 h). Scored window 2020-12-28 00Z-23Z = **h48-h71**, the
same 23 h that every earlier arm ran as its entire life. Rotated CONUS404
frames (item 58 fix), item-60 surface fix in both arms. ARM 1 NSCBC sigma=0.03,
keep_ramp=0, GPUs 0-1. ARM 2 Davies, identical in every other respect, GPUs 2-3.

---

## 1. What the surface fix should do (the strongest prediction)

Item 63 measured the interior running **2.5-3.6x too fast at 100 m** against
both d02 and CONUS404, with a drag coefficient **9-17x below physical**. The
30-step proof run now measures **Cd = 1.217e-3** (median 1.217e-3, range
5.4e-4..1.9e-3), inside the 1.0-1.5e-3 target.

> **P1.** The interior 100 m speed excess falls from 2.5-3.6x to **below 1.5x**
> in both arms. If it does not, item 63's causal reading is wrong and the
> excess has a second, independent source.

This is the prediction I most expect to be wrong in magnitude and right in
sign. Drag acting for 71 h on a spun-up boundary layer is a much stronger
constraint than drag acting for 23 h from rest, so a partial recovery would
still leave the sign intact while falsifying the "surface was the whole story"
reading.

> **P2.** Qstar over water stays nonzero for the whole run (proof: nonzero on
> 12899/12899 water cells, mean -7.62e-05). t_surf tracks SST over water
> (285.18-288.92 K) and sits at T0 = 288.00 K over land. Zero cells at the
> 271.0 K clamp at any output time.

## 2. NSCBC vs Davies at the inflow wall

The record has NSCBC at **210 mm/day** and Davies at **0.03 mm/day** in the
inflow wall cells. That is a factor of 7000 and it is the single largest
scheme-dependent difference the campaign has measured.

> **P3.** The difference survives at 71 h and stays above 100x. Both numbers
> move -- Davies up, NSCBC down -- because 48 h of lead removes the cold-start
> transient that both were sitting in, but they do not converge.

> **P4.** Davies remains the better-behaved scheme at the wall and the worse
> one in the interior. If Davies wins on both, NSCBC has no case in this
> configuration and the paper's conclusion inverts.

## 3. The wall enhancement

Item 56 measured ~4x at every terrain height 100-1000 m in the dN 0-12 band,
and the excess followed the **wall**, not the terrain. The decay fit gives
A = 0.99 +/- 0.70, B = 2.73 +/- 0.65, tau = 10.6 +/- 5.1 h.

> **P5.** With 48 h of lead, the h48-h71 window sits far outside the fitted
> transient (t >> tau), so the fit's asymptote A applies: near-wall enhancement
> in dN 0-6 lands at **1.0 +/- 0.7**, i.e. plausibly gone.

The uncertainty on A spans 0.3-1.7 and does not discriminate. This is the
weakest prediction here and I am recording it as a range, not a value. If the
enhancement is still ~4x at h48-h71, the decay fit is refuted outright and the
wall effect is a standing feature, not a transient.

> **P6.** The north/east asymmetry (3.79x vs 1.27x) persists in ARM 1. It is a
> property of which faces admit air, and 48 h of lead does not change the
> synoptic flow direction.

## 4. Santa Ynez

Item 62: the rotation fix moved Santa Ynez from **2.13x to 0.83x** against d02,
meeting a falsification criterion stated in advance. At pin 35.4 it sits at
dN 34, well outside the near-wall band.

> **P7.** Santa Ynez stays in **0.7-1.2x** in both arms. The rotation fix is
> upstream of the boundary scheme, so the two arms should agree here to within
> ~15%. Arms disagreeing by more than that would mean the boundary scheme
> reaches 34 cells into the interior, which nothing in the record supports.

## 5. Davies stability

Davies is the incumbent and is expected to survive; NSCBC is the scheme with
the measured 4.7 h failure at 192x96x48 under a global mass constraint
(UPSTREAM_ISSUES 30o), which `nscbc_mass_tau=60` and `du_max=8` were tuned to
suppress.

> **P8.** If either arm crashes it is **NSCBC**, not Davies, and it happens in
> the first 15 h of model time or not at all. Davies checkpoints hourly
> throughout precisely so this prediction can be wrong cheaply.

Restarts are part of the result, not an embarrassment to be hidden: a scheme
that needs three restarts to cover 71 h has said something about itself.

## 6. What would make this whole comparison uninteresting

> **P9.** If the two arms agree to within 10% on every interior metric, then
> the boundary scheme does not matter at this resolution and domain size, and
> the honest paper is a null result. That outcome is on the record here as a
> real possibility, not a fallback.

---

## Known-open items that could contaminate this run

- **Item 64** (init at rest, `ERF.cpp:2391`): unfixed by design. The 48 h lead
  is the mitigation, not a fix. If h48 still shows a rest signature, the lead
  is too short.
- **Item 32** (float32 absolute time): at t = 255,600 s the float32 ulp is
  2^-6 = 0.015625 s against dt ~ 0.43 s at cfl 0.2, so ~27 ulps per step. No
  degenerate steps expected mid-run. The truncation step at `stop_datetime` is
  the known hazard, and `plt_guard` covers it in all seven scoring scripts.
- **Item 60's second defect** (`lmask` never populated) is fixed here but was
  only visible after the first; there may be further consumers of the surface
  state that were never exercised while it was degenerate.
