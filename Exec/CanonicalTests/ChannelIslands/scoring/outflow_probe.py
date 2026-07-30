#!/usr/bin/env python3
"""Is the NSCBC outflow-wall excess driven by the mass controller?

The controller adds a uniform du on OUTFLOW FACES ONLY. Raising its clamp from
2 to 8 m/s is what made NSCBC mass-bounded -- and the excess precipitation sits
against the outflow walls. So: does the excess depend on the clamp?

Two runs already on disk, no GPU needed:
  run_nsc    -- 6.77 h probe, nscbc_mass_tau=60 at the DEFAULT du_max=2, where
                the controller is SATURATED (mass ran +2.2%)
  run_ab_nsc -- the scored day, du_max=8, controller unsaturated (mass -0.030%)

Same physics otherwise. If the outflow-wall profile is the same in both, the
controller is exonerated and the Riemann outflow condition owns it.
"""
import sys
import numpy as np
import yt

yt.set_log_level(50)
BINS = [(0, 0), (1, 2), (3, 5), (6, 9), (10, 14), (15, 19), (20, 29), (30, 95)]


def profile(pltdir, label):
    ds = yt.load(pltdir)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    t = float(ds.current_time) / 3600.0
    NX, NY = ra.shape
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    outf = np.minimum(NX - 1 - ii, NY - 1 - jj)   # distance from xhi / yhi
    print(f'\n=== {label}  t={t:.2f} h  domain mean {np.nanmean(ra):.2f} mm ===')
    row = []
    for lo, hi in BINS:
        m = (outf >= lo) & (outf <= hi)
        row.append(f'{lo}-{hi}:{np.nanmean(ra[m]):.1f}')
    print('  mean rain_accum by distance from the OUTFLOW walls (cells:mm)')
    print('   ' + '  '.join(row))
    return t, np.nanmean(ra)


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        d, lab = spec.split('=', 1)
        profile(d, lab)
