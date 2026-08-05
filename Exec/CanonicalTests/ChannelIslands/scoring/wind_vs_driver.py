#!/usr/bin/env python3
"""ERF wind against its CONUS404 driver, on a common ASL axis. Domain-agnostic.

  wind_vs_driver.py <rundir> <plotfile_prefix> [<hours_into_run>]

Replaces the hardcoded-192x96 checks used on Domain A so the same metric works
on the 96x96x96 ocean domain. Everything is read from the plotfile and the frame
headers; nothing about the grid is assumed.

WHY THIS AXIS. The frame's z array and ERF's z_phys are BOTH heights above sea
level (ERF_WeatherDataInterpolation.cpp:426 interpolates the frame at the cell's
z_phys), so they are directly comparable. Comparing ERF's AGL profile against the
frame's ASL levels is NOT valid and produced a spurious "the IC is wrong" reading
earlier in this campaign.

THE NUMBER THAT MATTERS is the ratio |V|_ERF / |V|_driver at the same ASL height
and the same time. Reference points measured on Domain A, 1 h, full physics:
    buoyancy_type 1 : 3.12 at 26 m       water-class ratio 4.43
    buoyancy_type 2 : 1.76 at 26 m       water-class ratio 2.30
    driver target   : ~1.24-1.27
The ocean domain's job is to say whether the residual ~2.3 survives once terrain
is absent. If it does, the remaining error is terrain-independent.
"""
import glob
import os
import sys

import numpy as np
import yt

yt.set_log_level(50)

FRAMES = '/app/ERF/pod_data_A_run30h/CONUS404Data_3D'
BAND = 10                     # relaxation cells excluded at every lateral edge
_KEEP = []


def read_frame(path):
    with open(path, 'rb') as fh:
        nx, ny, nz, nd = np.fromfile(fh, dtype=np.int32, count=4)
        np.fromfile(fh, dtype=np.float32, count=nx * ny)          # lat
        np.fromfile(fh, dtype=np.float32, count=nx * ny)          # lon
        x = np.fromfile(fh, dtype=np.float32, count=nx)
        y = np.fromfile(fh, dtype=np.float32, count=ny)
        z = np.fromfile(fh, dtype=np.float32, count=nz)
        b = [np.fromfile(fh, dtype=np.float32, count=nx * ny * nz).reshape(nz, ny, nx)
             for _ in range(nd)]
    return x, y, z, np.hypot(b[1], b[2])


def frame_times():
    """Frame files sorted by the timestamp in their name."""
    out = []
    for p in sorted(glob.glob(os.path.join(FRAMES, 'ERF_IC_*.bin'))):
        s = os.path.basename(p)[7:-4]                             # YYYY_MM_DD_HH_MM
        out.append((s, p))
    return out


def plotfiles(rundir, prefix):
    out = []
    for p in glob.glob(os.path.join(rundir, prefix + '[0-9]*')):
        if not os.path.isdir(p) or '.old' in p:
            continue
        with open(os.path.join(p, 'Header')) as fh:
            L = fh.read().split('\n')
        out.append((float(L[3 + int(L[1])]), p))
    return sorted(out)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rundir = sys.argv[1]
    prefix = sys.argv[2]
    want_h = float(sys.argv[3]) if len(sys.argv) > 3 else None

    pfs = plotfiles(rundir, prefix)
    if not pfs:
        sys.exit(f'no plotfiles in {rundir}')
    t0, p0 = pfs[0]
    t1, p1 = (pfs[-1] if want_h is None else
              min(pfs, key=lambda a: abs(a[0] / 3600.0 - want_h)))

    dom_edges = {}

    def erf_prof(path, edges):
        ds = yt.load(path)
        _KEEP.append(ds)
        dom_edges['lo'] = np.array(ds.domain_left_edge)
        dom_edges['hi'] = np.array(ds.domain_right_edge)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        A = lambda n: np.asarray(g[('boxlib', n)], dtype=np.float64)
        u, v, zp = A('x_velocity'), A('y_velocity'), A('z_phys')
        nx, ny, _ = u.shape
        ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        M = np.minimum.reduce([ii, jj, nx - 1 - ii, ny - 1 - jj]) >= BAND
        sp = np.hypot(u, v)
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (zp >= lo) & (zp < hi) & M[:, :, None]
            out.append(np.median(sp[sel]) if sel.sum() > 20 else np.nan)
        return np.array(out), (nx, ny)

    # ERF selects frames by MODEL TIME FROM FRAME 0, not by start_datetime:
    #     idx1 = int(time / hindcast_data_interval)   (ERF_WeatherDataInterpolation.cpp:1839)
    # so a run whose deck says 2020-12-28 00Z but whose pool begins 2020-12-27
    # 18Z is actually initialised from 18Z. Comparing against the start_datetime
    # frame therefore mismatches by the pool offset -- which made a faithful IC
    # read as 2.4x too slow. Match what the model DOES, not what the deck says.
    ft = frame_times()
    i0 = 0
    print(f'driver frames (as ERF indexes them, from frame 0): '
          f'{ft[i0][0]} -> {ft[i0+1][0]}')
    x, y, z, s_a = read_frame(ft[i0][1])
    _, _, _, s_b = read_frame(ft[i0 + 1][1])
    dT = 3.0 * 3600.0                       # frames are 3-hourly in this pool
    w = min(max((t1 - t0) / dT, 0.0), 1.0)
    lv = [k for k in range(len(z)) if z[k] <= 3000.0]
    edges = [0.0] + [0.5 * (z[k] + z[k + 1]) for k in lv[:-1]] + [z[lv[-1]] * 1.05]
    e0, dims = erf_prof(p0, edges)
    e1, _ = erf_prof(p1, edges)

    # Restrict the driver to the SAME GEOGRAPHIC BOX as the ERF domain, minus the
    # same relaxation band. Taking the driver median over the whole frame would
    # compare an ocean-box ERF profile against a driver average that includes
    # land and mountains -- not a ratio of anything.
    lo, hi = dom_edges['lo'], dom_edges['hi']
    bx = BAND * 3000.0
    ix = np.where((x >= lo[0] + bx) & (x <= hi[0] - bx))[0]
    jy = np.where((y >= lo[1] + bx) & (y <= hi[1] - bx))[0]
    if ix.size < 5 or jy.size < 5:
        sys.exit('ERF domain does not overlap the frame as expected')
    print(f'driver restricted to frame i={ix[0]}..{ix[-1]}, j={jy[0]}..{jy[-1]} '
          f'({ix.size} x {jy.size} cells) to match the ERF domain')
    dr = lambda s: np.array([np.median(s[np.ix_([k], jy, ix)]) for k in range(len(z))])
    d0, d1 = dr(s_a), dr(s_a) + w * (dr(s_b) - dr(s_a))

    print(f'run {rundir}   grid {dims[0]}x{dims[1]}   '
          f't0={t0:.0f}s  t1={t1:.0f}s  (driver interpolated w={w:.2f})\n')
    hdr = (f'{"z ASL":>8s} {"drv t0":>8s} {"ERF t0":>8s} {"r0":>6s}  '
           f'{"drv t1":>8s} {"ERF t1":>8s} {"r1":>6s}')
    print(hdr); print('-' * len(hdr))
    r1s = []
    for n, k in enumerate(lv):
        a = e0[n] / d0[k] if d0[k] > 0 else np.nan
        b = e1[n] / d1[k] if d1[k] > 0 else np.nan
        r1s.append(b)
        print(f'{z[k]:8.1f} {d0[k]:8.2f} {e0[n]:8.2f} {a:6.2f}  '
              f'{d1[k]:8.2f} {e1[n]:8.2f} {b:6.2f}')
    r1s = np.array(r1s, dtype=float)
    print(f'\nmean |ratio-1| over 0-3 km at t1: {np.nanmean(np.abs(r1s - 1)):.3f}'
          f'   spread {np.nanstd(r1s):.3f}')
    print('Domain A reference, same metric: buoyancy_type 1 -> 1.263 / 0.889, '
          'type 2 -> 0.439 / 0.280')


if __name__ == '__main__':
    main()
