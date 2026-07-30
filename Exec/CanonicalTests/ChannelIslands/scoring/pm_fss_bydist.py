#!/usr/bin/env python3
"""Percentile-matched FSS as a function of distance from the lateral boundary.

Hypothesis under test: Davies places better because its relaxation band anchors
the interior toward ERA5's synoptic pattern, while NSCBC constrains only the wall
face. If so, Davies' advantage should be LARGEST near the band and DECAY inward.
If the advantage is uniform with distance, anchoring is not the mechanism.

Base rates are matched WITHIN EACH BIN (each field thresholded at its own
quantile over that bin's cells), so bias is divided out separately at every
distance -- the fetch-binned bias measurement already showed bias varies with
distance, so a single global threshold would leak that variation into the score.

Second read-out: the same bins scored against ERA5's precipitation field instead
of MRMS. Davies tracking ERA5 more closely than NSCBC does, with the gap decaying
inward, is direct evidence of inherited driver placement.

Geometry note: d = min(i, j, NX-1-i, NY-1-j) on a 192x96 grid maxes at 47, so a
"60+" bin is empty and "45-60" holds only d = 45..47.
"""
import sys
import numpy as np
from scipy.ndimage import uniform_filter

WINS = [(3, 1), (9, 3), (15, 5), (30, 11), (60, 21)]
BINS = [(20, 29), (30, 39), (40, 47)]


def fss_binary(fb, ob, w, mask):
    Pf = uniform_filter(fb.astype(float), size=w, mode='constant')
    Po = uniform_filter(ob.astype(float), size=w, mode='constant')
    num = np.nanmean((Pf[mask] - Po[mask]) ** 2)
    den = np.nanmean(Pf[mask] ** 2) + np.nanmean(Po[mask] ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def binarize(field, mask, rate):
    thr = float(np.nanquantile(field[mask], 1.0 - rate))
    return field >= thr


def load(d):
    d = d.rstrip('/') + '/'
    return (np.load(d + 'erf_mm.npy'), np.load(d + 'mrms_mm.npy'),
            np.load(d + 'era5_mm.npy'))


def main():
    dav_dir, nsc_dir = sys.argv[1], sys.argv[2]
    erf_d, mrms, era5 = load(dav_dir)
    erf_n, _, _ = load(nsc_dir)
    NX, NY = mrms.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    print(f'domain {NX}x{NY}; d ranges 0..{int(d.max())} (a 60+ bin is empty by geometry)')

    for refname, ref in (('MRMS', mrms), ('ERA5', era5)):
        for phys in (1.0, 5.0):
            base = np.isfinite(ref) & (d >= 20)
            rate = float((ref[base] >= phys).mean())
            print(f'\n=== vs {refname}, base rate {rate:.4f} (= {refname} exceedance of '
                  f'{phys:.0f} mm over d>=20), matched within each bin ===')
            print(f'{"bin (cells)":>12} {"km":>9} {"n":>6} '
                  f'{"3km Dav":>8} {"3km NSC":>8} {"delta":>7}   '
                  f'{"60km Dav":>9} {"60km NSC":>9} {"delta":>7}')
            for lo, hi in BINS:
                m = (d >= lo) & (d <= hi) & np.isfinite(ref)
                if m.sum() < 50:
                    continue
                r = float((ref[m] >= phys).mean())
                if r <= 0 or r >= 1:
                    print(f'{lo}-{hi}: degenerate base rate {r}')
                    continue
                ob = binarize(ref, m, r)
                assert abs(fss_binary(ob, ob, 3, m) - 1.0) < 1e-12, 'self-FSS control failed'
                row = {}
                for name, fld in (('Dav', erf_d), ('NSC', erf_n)):
                    fb = binarize(fld, m, r)
                    fb2 = binarize(fld * 3.7, m, r)
                    assert np.array_equal(fb, fb2), 'bias-invariance control failed'
                    row[name] = {km: fss_binary(fb, ob, w, m) for km, w in WINS}
                print(f'{lo:>5}-{hi:<6} {3*lo:>4}-{3*hi:<4} {int(m.sum()):>6} '
                      f'{row["Dav"][3]:>8.3f} {row["NSC"][3]:>8.3f} '
                      f'{row["Dav"][3]-row["NSC"][3]:>+7.3f}   '
                      f'{row["Dav"][60]:>9.3f} {row["NSC"][60]:>9.3f} '
                      f'{row["Dav"][60]-row["NSC"][60]:>+7.3f}')


if __name__ == '__main__':
    main()
