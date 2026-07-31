#!/usr/bin/env python3
"""Is the premise for re-opening NSCBC true on the PRODUCTION deck?

The claim: NSCBC eliminates the Davies band and delivers more moisture to the
interior, which is exactly what #25 shows the interior lacks. Measure it on the
one pair we already have -- gate take 2 (NSCBC, died 4.7 h) against the Davies
day, same deck, same hours, hours 1-4 only.

Reported per hour, on the interior mask d >= 20 (the #25 science cells) and on
the wall band d = 0:
  rain_accum mean  [mm]          -- what actually reached the ground
  column qv        [kg/m2]       -- integrated water vapour
"""
import sys
import numpy as np
import yt

from plt_guard import reject_if_poisoned

yt.set_log_level(50)
NX, NY = 192, 96
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
D = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
INT = D >= 20
BAND = D == 0


def one(pltdir):
    ds = yt.load(pltdir)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    t = float(ds.current_time) / 3600.0
    ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    reject_if_poisoned(pltdir, pltdir, ra)
    qv = np.asarray(g[('boxlib', 'qv')])
    rho = np.asarray(g[('boxlib', 'density')])
    z = np.asarray(g[('boxlib', 'z_phys')])
    dz = np.diff(z, axis=2, append=z[:, :, -1:] + (z[:, :, -1:] - z[:, :, -2:-1]))
    cwv = (qv * rho * dz).sum(axis=2)
    return t, ra, cwv


def main():
    label, dirs = sys.argv[1], sys.argv[2:]
    print(f"=== {label} ===")
    print(f"{'t(h)':>5} {'rain_int':>9} {'rain_band':>10} {'cwv_int':>9} {'cwv_band':>9}")
    for p in dirs:
        t, ra, cwv = one(p)
        print(f"{t:5.2f} {np.nanmean(ra[INT]):9.4f} {np.nanmean(ra[BAND]):10.3f} "
              f"{np.nanmean(cwv[INT]):9.3f} {np.nanmean(cwv[BAND]):9.3f}")


if __name__ == "__main__":
    main()
