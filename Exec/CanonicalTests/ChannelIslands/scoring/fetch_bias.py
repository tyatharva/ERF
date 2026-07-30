#!/usr/bin/env python3
"""Fetch-binned bias: is the NSCBC wet interior a residual band, or uniform?

Fetch = distance in cells from the nearest INFLOW wall. Measured on both arms,
the inflow walls for this day are xlo and ylo (face-mean outward-positive u_n is
negative there all day: xlo -1.2 to -8.3, ylo -5.7 to -9.1); xhi and yhi carry
outflow. So fetch = min(i, j).

Read-out is the ratio ERF/MRMS per fetch bin, which divides out the observed
field's own spatial gradient (the day has a strong north-south gradient, so a
raw ERF profile would confound the two).

  ratio DECAYS with fetch -> the residual band at the inflow wall is still the
                             source, and the fix is local to the wall treatment
  ratio FLAT             -> uniform over-delivery, implicating the outflow
                             condition or the mass controller instead

Control: the Davies arm, whose band is known and steep (210 mm/day at d=0
decaying monotonically to 4.4 by d=25). If the instrument cannot resolve THAT,
it cannot resolve a subtler one.
"""
import sys
import numpy as np

BINS = [(0, 0), (1, 2), (3, 5), (6, 9), (10, 14), (15, 19),
        (20, 29), (30, 44), (45, 95)]


def load(d):
    d = d.rstrip('/') + '/'
    return np.load(d + 'erf_mm.npy'), np.load(d + 'mrms_mm.npy')


def profile(erf, mrms, fetch, label):
    print(f'\n=== {label} ===')
    print(f'{"fetch (cells)":>14} {"km":>7} {"n":>6} {"ERF":>9} {"MRMS":>8} {"ratio":>8}')
    for lo, hi in BINS:
        m = (fetch >= lo) & (fetch <= hi) & np.isfinite(mrms)
        if not m.any():
            continue
        e, o = float(erf[m].mean()), float(mrms[m].mean())
        print(f'{lo:>6}-{hi:<7} {3*lo:>4}-{3*hi:<3} {int(m.sum()):>6} '
              f'{e:>9.2f} {o:>8.2f} {e/o if o > 0 else float("nan"):>8.2f}')


def main():
    dav_dir, nsc_dir = sys.argv[1], sys.argv[2]
    erf_d, mrms = load(dav_dir)
    erf_n, _ = load(nsc_dir)
    NX, NY = mrms.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    fetch = np.minimum(ii, jj)          # distance from the xlo / ylo inflow walls
    outf = np.minimum(NX - 1 - ii, NY - 1 - jj)   # distance from the outflow walls

    profile(erf_d, mrms, fetch, 'DAVIES (control: known steep band)')
    profile(erf_n, mrms, fetch, 'NSCBC by INFLOW fetch')
    profile(erf_n, mrms, outf, 'NSCBC by OUTFLOW distance (does the far wall matter too?)')

    # Interior-only slopes, so the band cells cannot dominate the trend.
    m = (fetch >= 20) & np.isfinite(mrms)
    for name, erf in (('Davies', erf_d), ('NSCBC', erf_n)):
        r = np.where(mrms > 0, erf / np.maximum(mrms, 1e-9), np.nan)
        f = fetch[m].astype(float)
        rr = r[m]
        ok = np.isfinite(rr)
        slope = np.polyfit(f[ok], rr[ok], 1)[0]
        print(f'\n{name}: interior (fetch>=20) ratio vs fetch, least-squares slope '
              f'{slope:+.4f} per cell ({slope*3:+.4f} per 3 km); '
              f'mean ratio {np.nanmean(rr):.2f}')


if __name__ == '__main__':
    main()
