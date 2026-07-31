#!/usr/bin/env python3
"""Read the sigma sweep at a matched model time.

Reports, per arm, on the 6 h plotfile:
  interior d>=20 mean rain_accum   -- the quantity whose 24 h value is the bias
  band d=0 mean                    -- the wall artifact
  mean by distance from the OUTFLOW walls -- the Riemann structure task 2 found

The two endpoints (sigma=1 = run_ab_nsc, sigma=0 = run_ab_nscx) anchor the scale:
their 24 h interior biases are 2.63x and 0.19x, so the sigma whose 6 h interior
mean falls between theirs in the right proportion is the candidate for unity.
"""
import sys
import numpy as np
import yt

from plt_guard import reject_if_poisoned

yt.set_log_level(50)
BINS = [(0, 0), (1, 2), (3, 5), (6, 9), (10, 14), (15, 19), (20, 29), (30, 95)]


def one(pltdir, label):
    ds = yt.load(pltdir)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    reject_if_poisoned(label, pltdir, ra)
    t = float(ds.current_time) / 3600.0
    NX, NY = ra.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    d = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    outf = np.minimum(NX - 1 - ii, NY - 1 - jj)
    prof = '  '.join(f'{lo}-{hi}:{np.nanmean(ra[(outf >= lo) & (outf <= hi)]):.1f}'
                     for lo, hi in BINS)
    print(f'{label:<24} t={t:5.2f}h  interior(d>=20)={np.nanmean(ra[d >= 20]):8.3f}  '
          f'band(d=0)={np.nanmean(ra[d == 0]):8.2f}  domain={np.nanmean(ra):7.2f}')
    print(f'{"":24}   by outflow distance: {prof}')


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        p, lab = spec.split('=', 1)
        try:
            one(p, lab)
        except Exception as e:
            print(f'{lab}: FAILED {e}')
