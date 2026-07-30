# Scoring instruments and run harnesses

Every number in UPSTREAM_ISSUES items 25-39 came from these. They lived in a
session scratchpad and are checked in here so the results are reproducible.

Paths are container-absolute (`/app/ERF/...`); run them inside the `erf-hindcast`
image with the repo bind-mounted at `/app/ERF`.

## Pipeline
1. `score_foundation.py <rundir> <outdir>` -- builds matched 24-h fields on the
   ERF grid: `erf_mm`, `mrms_mm`, `era5_mm`, `lat`, `lon`, `terrain`. Gates on a
   PROJECTION residual against the frame's own coordinates (<10 m) and on
   regrid-vs-native interpolation controls. **The projection gate is not
   optional** -- item 25 was first scored 85.9 km off because the LCC came from
   the ERA5 download area instead of the pinned projection area.
2. `score_metrics.py <scoredir>` -- fixed-threshold FSS, categorical, bias,
   spectrum, N/S gradient, wall-band scan.
3. `score_islands.py <scoredir>` -- per-island precipitation.

## Bias-free placement (item 36d: fixed-threshold FSS rewards coverage)
- `pm_fss.py <dav> <arm>` -- percentile-matched FSS. Asserts self-FSS==1, bias
  invariance, and a known-nonzero displacement response.
- `pm_fss_bydist.py`, `pm_fss_hourly.py` -- the distance- and time-resolved
  variants that separated advective inheritance (item 35) from local anchoring.
- `variance_matched.py` -- smooths to a matched spectral ratio; falsified the
  smoothness-artifact caveat (36d).
- `scale_filtered_fss.py` -- DCT low-pass then score; decided against building
  spectral nudging (36c).
- `mode_objects.py` -- MODE-style object verification.

## Diagnostics
`fetch_bias.py` (inflow/outflow binning), `nscbc_walls.py` (wall regime),
`outflow_probe.py`, `interior_moisture.py`, `sig_probe.py`, `proj_truth.py`
(recovers the true frame LCC from a frame's own coordinates).

## CONUS404 / wrfout case (item 39)
`chk_wrfout.py`, `chk_mrms.py`, `chk_bucket.py` (RAINNC wraps at BUCKET_MM=100 --
naive differencing gives -99.8 mm cells), `chk_ladder.py`, `validate_wrf.py`,
`wrf_struct.py`, `wrf_maps.py`.

## Harnesses (`../harness/`)
`ab_day.sh <arm> [knobs...]` -- the bracketed full-day runner used for every
scored arm; selects its leg-2 restart via `pick_restart_chk.sh` (item 31).
`sig_sweep.sh`, `nsc_probe.sh`, `nsc_tau.sh`, `radcheck.sh`, `build_verify.sh`
(checks the binary mtime, not the exit code -- a piped build once reported
success on a failure), `nanloc.sh`, `run48.sh`, `perf_probes.sh`.
