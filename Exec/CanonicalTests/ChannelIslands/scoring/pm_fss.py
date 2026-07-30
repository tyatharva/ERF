#!/usr/bin/env python3
"""Percentile-matched FSS: placement skill with the bias divided out.

Fixed-threshold FSS rewards a wet field through coverage alone, which is exactly
the confound in the NSCBC-vs-Davies comparison (NSCBC is 2.63x wet interior,
Davies 0.38x dry). Here each field is thresholded at its OWN quantile, chosen so
every binary field has the SAME base rate as the observations do at 1 mm and
5 mm. Bias is then removed by construction and what is left is placement.

Instrument controls, all asserted, not assumed:
  1. self-FSS == 1 for MRMS against itself.
  2. BIAS INVARIANCE: multiplying a forecast by any positive constant must leave
     percentile-matched FSS bit-identical. This is the property the metric is
     being used for, so it is tested directly.
  3. DISPLACEMENT RESPONSE (known-nonzero): MRMS shifted by 5 cells must score
     BELOW 1 at small scales and recover with scale. A metric that cannot see a
     15 km displacement cannot be trusted to see a real one.
"""
import sys
import numpy as np
from scipy.ndimage import uniform_filter

WINS = [(3, 1), (9, 3), (15, 5), (30, 11), (60, 21)]


def fss_binary(fb, ob, w, mask):
    Pf = uniform_filter(fb.astype(float), size=w, mode='constant')
    Po = uniform_filter(ob.astype(float), size=w, mode='constant')
    num = np.nanmean((Pf[mask] - Po[mask]) ** 2)
    den = np.nanmean(Pf[mask] ** 2) + np.nanmean(Po[mask] ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def binarize_at_rate(field, mask, rate):
    """Threshold `field` so exactly `rate` of the masked cells exceed it."""
    if rate <= 0:
        return np.zeros_like(field, dtype=bool), np.inf
    thr = float(np.nanquantile(field[mask], 1.0 - rate))
    return (field >= thr), thr


def load(d):
    d = d.rstrip('/') + '/'
    return (np.load(d + 'erf_mm.npy'), np.load(d + 'mrms_mm.npy'),
            np.load(d + 'era5_mm.npy'))


def main():
    dav_dir, nsc_dir = sys.argv[1], sys.argv[2]
    erf_d, mrms, era5 = load(dav_dir)
    erf_n, mrms_n, _ = load(nsc_dir)
    assert np.allclose(mrms, mrms_n, equal_nan=True), 'the two arms must share one observation field'

    NX, NY = mrms.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    mask = (d >= 20) & np.isfinite(mrms)
    print(f'interior mask: {mask.sum()} cells')

    # --- controls ---
    for rate in (0.5,):
        ob, _ = binarize_at_rate(mrms, mask, rate)
        assert abs(fss_binary(ob, ob, 3, mask) - 1.0) < 1e-12, 'self-FSS control failed'
    fb_a, _ = binarize_at_rate(erf_n, mask, 0.3)
    fb_b, _ = binarize_at_rate(erf_n * 2.63, mask, 0.3)
    assert np.array_equal(fb_a, fb_b), 'bias-invariance control failed (binarization)'
    ob, _ = binarize_at_rate(mrms, mask, 0.3)
    assert abs(fss_binary(fb_a, ob, 5, mask) - fss_binary(fb_b, ob, 5, mask)) < 1e-15, \
        'bias-invariance control failed (score)'
    shifted = np.roll(mrms, 5, axis=0)
    sb, _ = binarize_at_rate(shifted, mask, 0.3)
    s3, s63 = fss_binary(sb, ob, 1, mask), fss_binary(sb, ob, 21, mask)
    assert s3 < 0.99 and s63 > s3, f'displacement control failed ({s3:.3f}, {s63:.3f})'
    print(f'controls OK: self-FSS=1; bias-invariant; 15 km displacement scores '
          f'{s3:.3f} at 3 km rising to {s63:.3f} at 63 km')

    # --- base rates taken from the OBSERVATIONS at the two physical thresholds ---
    for phys in (1.0, 5.0):
        rate = float((mrms[mask] >= phys).mean())
        ob, thr_o = binarize_at_rate(mrms, mask, rate)
        useful = 0.5 + rate / 2.0
        print(f'\n=== percentile-matched FSS, base rate {rate:.4f} '
              f'(= MRMS exceedance of {phys:.0f} mm), useful >= {useful:.3f} ===')
        cols = {}
        for name, fld in (('Davies', erf_d), ('NSCBC', erf_n), ('ERA5', era5)):
            fb, thr = binarize_at_rate(fld, mask, rate)
            cols[name] = (fb, thr)
        print(f'  matched thresholds (mm): MRMS {thr_o:6.2f} | ' +
              ' | '.join(f'{n} {cols[n][1]:6.2f}' for n in cols))
        print(f'  {"scale":>6} {"Davies":>8} {"NSCBC":>8} {"ERA5":>8}   {"NSCBC-Davies":>12}')
        for km, w in WINS:
            v = {n: fss_binary(cols[n][0], ob, w, mask) for n in cols}
            print(f'  {km:>4} km {v["Davies"]:>8.3f} {v["NSCBC"]:>8.3f} {v["ERA5"]:>8.3f}   '
                  f'{v["NSCBC"] - v["Davies"]:>+12.3f}')

        # Independent bias-free placement check, on a different statistic.
        # With the base rates matched, n_forecast == n_obs, so FAR == 1 - POD and
        # CSI = H/(2n-H): POD alone carries the whole comparison, and it shares no
        # machinery with the neighbourhood score above.
        line = []
        for n in cols:
            H = int((cols[n][0] & ob & mask).sum())
            no = int((ob & mask).sum())
            nf = int((cols[n][0] & mask).sum())
            pod = H / no
            line.append(f'{n} POD {pod:.3f} CSI {H/(no+nf-H):.3f}')
        print('  grid-scale, matched rate: ' + ' | '.join(line))


if __name__ == '__main__':
    main()
