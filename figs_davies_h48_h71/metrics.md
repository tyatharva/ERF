# Davies control, h48-h71 (2020-12-28 00Z-23Z) -- metrics

Mask: land-only, terrain > 20 m AND finite in all three sources.
n = 5536 of 18432 cells (30.0%).
Regridding: both references BIN-AVERAGED (area-averaged) from finer
grids -- d02 1.5 km, MRMS ~1 km, median 9.0 source cells per target
per hour. No nearest-neighbour sampling. NN would retain 1 km extrema
and inflate reference maxima against a 3 km model field.
Restart discontinuity: the arm rolled back to the h8 checkpoint once;
the scored window lies entirely inside the surviving leg and window
accumulation had zero negative cells.

| metric | Davies vs d02 | Davies vs MRMS |
|---|---|---|
| ERF domain mean (mm) | 47.226 | 47.226 |
| reference mean (mm) | 18.583 | 18.849 |
| ratio | 2.54x | 2.51x |
| bias (mm) | +28.643 | +28.376 |
| RMSE (mm) | 184.707 | 184.992 |
| correlation | 0.048 | 0.012 |
| near/far (dN<=12 / dN>=21) | ERF 1.25, d02 0.46 | MRMS 0.56 |

## Terrain-stratified ERF/reference, dN 0-6

| terrain | n | ERF mm | d02 mm | ERF/d02 | MRMS mm | ERF/MRMS |
|---|---|---|---|---|---|---|
| 20-100 m | 64 | 20.9 | 9.4 | 2.24x | 10.2 | 2.05x |
| 100-300 m | 156 | 41.3 | 8.2 | 5.06x | 8.9 | 4.64x |
| 300-600 m | 140 | 114.3 | 18.3 | 6.25x | 21.1 | 5.43x |
| 600-1000 m | 261 | 100.7 | 9.5 | 10.65x | 13.2 | 7.63x |
| 1000-3000 m | 159 | 120.6 | 6.3 | 19.01x | 5.8 | 20.69x |

## Instrument control (run BEFORE trusting the correlations)

d02 vs MRMS on the identical grid and mask: **corr 0.840**, ratio 0.99,
RMSE 8.06 mm. The two independent references agree strongly, so the grids are
aligned and the near-zero ERF correlations below are REAL, not a regridding or
orientation bug.

## The aggregate numbers are dominated by the lateral bands

| mask | n | ERF mm | d02 mm | ratio | corr(d02) | corr(MRMS) |
|---|---|---|---|---|---|---|
| all land | 5536 | 47.2 | 18.6 | 2.54x | 0.048 | 0.012 |
| dN > 10 (NORTH band only excluded) | 4324 | 41.7 | 20.8 | 2.00x | 0.045 | -0.003 |
| **all four bands excluded, d >= 10** | 3692 | 28.8 | 22.2 | **1.29x** | **0.541** | 0.405 |
| all four bands, d >= 20 | 2221 | 32.5 | 26.8 | 1.22x | 0.542 | 0.351 |

Excluding the NORTH band alone changes nothing (0.048 -> 0.045). Excluding ALL
FOUR lateral bands moves the correlation from 0.05 to **0.54** and the ratio
from 2.54x to **1.29x**.

**This corrects the campaign's framing.** Items 56/57 stratified by dN and
concluded the excess follows the NORTH wall. On this arm it is present at every
lateral boundary, and the worst single cell is at the SOUTH wall: ERF's maximum
is **5761 mm / 23 h at dN = 94** (j ~ 1, the ylo band) on 42 m terrain -- the
same band where the dead legs' w dipole fired. RMSE 184.7 mm and the p99 of
579 mm are driven by these band cells, not by the interior.

Read plainly: in the interior this arm is a **1.2-1.3x wet bias with real
spatial skill (r ~ 0.54)**. The 2.5x and r ~ 0 headline is a boundary-band
artifact contaminating a domain that is only 192x96, where the four 10-cell
bands are 33% of the land cells.
