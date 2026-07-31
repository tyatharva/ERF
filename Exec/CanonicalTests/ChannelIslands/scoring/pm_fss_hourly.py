#!/usr/bin/env python3
"""Percentile-matched FSS hour by hour, interior only -- the time-based
discriminator for whether interior placement is inherited advectively.

Transit time across this domain is ~576 km / ~20 m/s = ~8 h, so:
  gap GROWS over the first ~8 h  -> placement is inherited from what arrives at
                                    the inflow wall; the lever is the driving
                                    data, not the BC scheme
  gap present early and FLAT     -> placement is intrinsic to the schemes'
                                    interior behaviour; no boundary work helps

Reference is the MRMS 24 h total, held FIXED across hours. That depresses the
absolute score early (a 3 h accumulation cannot match a 24 h pattern), but both
arms are scored against the SAME reference at every hour, so the GAP between them
is the quantity that carries meaning -- and the gap is what the test reads.

Degeneracy guard: in the first hours a model may have fewer wet cells than the
target base rate, so its quantile threshold collapses to 0 and the binarization
cannot hit the rate. Any hour where either arm's achieved rate misses the target
by more than 10% is reported as degenerate and excluded, not scored.
"""
import glob
import sys
import numpy as np
import yt
from scipy.ndimage import uniform_filter
from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)


def fss_binary(fb, ob, w, mask):
    Pf = uniform_filter(fb.astype(float), size=w, mode='constant')
    Po = uniform_filter(ob.astype(float), size=w, mode='constant')
    num = np.nanmean((Pf[mask] - Po[mask]) ** 2)
    den = np.nanmean(Pf[mask] ** 2) + np.nanmean(Po[mask] ** 2)
    return 1.0 - num / den if den > 0 else np.nan


def binarize(field, mask, rate):
    thr = float(np.nanquantile(field[mask], 1.0 - rate))
    fb = field >= thr
    return fb, float(fb[mask].mean())


def series(run):
    out = {}
    for p in plotfiles_by_time(run):
        try:
            ds = yt.load(p)
        except Exception:
            continue
        t = float(ds.current_time)
        h = int(round(t / 3600.0))
        # `h in out` matches inflow_flux.py:51. Without it BOTH the h23 plotfile
        # and the end-of-run one land on h=23 (t=82816 is only 16 s off), and
        # the later one wins -- so a run that took a degenerate zero-length
        # final step silently OVERWROTE hour 23 with its poisoned rain_accum.
        if abs(t - h * 3600.0) > 400.0 or h == 0 or h in out:
            continue
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
        # Ordering alone is too fragile to rely on; reject on the data.
        if is_poisoned(ra):
            continue
        out[h] = ra
    return out


def main():
    dav_run, nsc_run, scoredir = sys.argv[1], sys.argv[2], sys.argv[3]
    mrms = np.load(scoredir.rstrip('/') + '/mrms_mm.npy')
    NX, NY = mrms.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    mask = (d >= 20) & np.isfinite(mrms)
    print(f'interior mask {mask.sum()} cells; d maxes at {int(d.max())} on this {NX}x{NY} grid')

    A, B = series(dav_run), series(nsc_run)
    hours = sorted(set(A) & set(B))
    print(f'hours with both arms: {hours[0]}..{hours[-1]} ({len(hours)})')

    for phys in (1.0, 5.0):
        rate = float((mrms[mask] >= phys).mean())
        ob, _ = binarize(mrms, mask, rate)
        assert abs(fss_binary(ob, ob, 3, mask) - 1.0) < 1e-12, 'self-FSS control failed'
        print(f'\n=== matched base rate {rate:.4f} (MRMS >= {phys:.0f} mm), interior d>=20 ===')
        print(f'{"hour":>5} {"Dav 3km":>8} {"NSC 3km":>8} {"gap":>7}   '
              f'{"Dav 60km":>9} {"NSC 60km":>9} {"gap":>7}   note')
        for h in hours:
            fa, ra = binarize(A[h], mask, rate)
            fb, rb = binarize(B[h], mask, rate)
            bad = (abs(ra - rate) > 0.1 * rate) or (abs(rb - rate) > 0.1 * rate)
            if bad:
                print(f'{h:>5} {"":>8} {"":>8} {"":>7}   {"":>9} {"":>9} {"":>7}   '
                      f'DEGENERATE (achieved {ra:.3f}/{rb:.3f} vs {rate:.3f}) -- too few wet cells')
                continue
            a3, b3 = fss_binary(fa, ob, 1, mask), fss_binary(fb, ob, 1, mask)
            a6, b6 = fss_binary(fa, ob, 21, mask), fss_binary(fb, ob, 21, mask)
            print(f'{h:>5} {a3:>8.3f} {b3:>8.3f} {a3-b3:>+7.3f}   '
                  f'{a6:>9.3f} {b6:>9.3f} {a6-b6:>+7.3f}')


if __name__ == '__main__':
    main()
