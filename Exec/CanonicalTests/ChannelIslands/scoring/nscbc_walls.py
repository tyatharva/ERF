#!/usr/bin/env python3
"""19d verification on the PRODUCTION deck: does the NSCBC gate-take-2 run show
the wall-reversal signature before its 4.7 h dt-collapse death?

For each lateral face and each hourly plotfile, report the wall-adjacent-cell
normal velocity made OUTWARD-POSITIVE, at the lowest few levels:
  mean, min, and the count of INFLOW cells (outward-negative).

The 19d claim is that the failure follows a face reversing from outflow to
inflow.  Control: the same measurement on the Davies day (run_a3), which
survived -- a face that never reverses there is the known-nonzero control for
this instrument.
"""
import sys
import numpy as np
import yt

yt.set_log_level(50)


def walls(pltdir, nk=4):
    ds = yt.load(pltdir)
    dims = ds.domain_dimensions
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)
    u = np.array(cg["x_velocity"])   # cell-centred in the plotfile
    v = np.array(cg["y_velocity"])
    t = float(ds.current_time)
    out = {}
    # outward-positive wall-adjacent normal velocity, lowest nk levels
    out["xlo"] = -u[0, :, :nk]
    out["xhi"] = +u[-1, :, :nk]
    out["ylo"] = -v[:, 0, :nk]
    out["yhi"] = +v[:, -1, :nk]
    return t, out


def main():
    label = sys.argv[1]
    print(f"=== {label} ===")
    print(f"{'t(h)':>6} {'face':>5} {'mean_un':>9} {'min_un':>9} {'max|un|':>9} {'n_inflow':>9} {'n_tot':>6}")
    for pltdir in sys.argv[2:]:
        t, out = walls(pltdir)
        for face in ("xlo", "xhi", "ylo", "yhi"):
            a = out[face]
            print(f"{t/3600:6.2f} {face:>5} {a.mean():9.3f} {a.min():9.3f} "
                  f"{np.abs(a).max():9.3f} {int((a < 0).sum()):9d} {a.size:6d}")


if __name__ == "__main__":
    main()
