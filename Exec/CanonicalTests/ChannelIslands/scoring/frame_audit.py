#!/usr/bin/env python3
"""Audit a directory of hindcast frames before committing GPU time to them.

  frame_audit.py <frame_dir> [3d|sfc]

Three things, none of which the model checks for you:

  1. NaN / inf anywhere in any field.
  2. Physical range per variable. A frame can parse, interpolate and run while
     carrying a field that is quietly wrong -- item 60's surface frames were
     structurally perfect and read as longitudes.
  3. Timestamp parse vs the filename, and even 3-hourly spacing with no gaps.
     ERF indexes frames POSITIONALLY from t=0, so one missing frame silently
     shifts every later frame three hours early.

Layout (write_frame): 4 int32 header, lat[nx*ny], lon[nx*ny], x[nx], y[ny],
z[nz], then nfields blocks of nx*ny*nz, each in k,j,i order.
"""
import datetime as dt
import glob
import os
import sys

import numpy as np

RANGE_3D = [('rho', 0.05, 1.6), ('u', -120., 120.), ('v', -120., 120.),
            ('w', -30., 30.), ('theta', 220., 800.), ('qv', 0., 0.05),
            ('qc', 0., 0.02), ('qr', 0., 0.02)]
RANGE_SFC = [('sst', 250., 320.), ('zero1', -1e-6, 1e-6), ('zero2', -1e-6, 1e-6),
             ('zero3', -1e-6, 1e-6), ('lmask', 0., 1.), ('alb', 0., 1.)]


def read(path):
    with open(path, 'rb') as f:
        nx, ny, nz, nd = np.fromfile(f, dtype='<i4', count=4)
        lat = np.fromfile(f, dtype='<f4', count=nx * ny)
        lon = np.fromfile(f, dtype='<f4', count=nx * ny)
        xs = np.fromfile(f, dtype='<f4', count=nx)
        ys = np.fromfile(f, dtype='<f4', count=ny)
        zs = np.fromfile(f, dtype='<f4', count=nz)
        flds = [np.fromfile(f, dtype='<f4', count=nx * ny * nz) for _ in range(nd)]
        trailing = len(f.read())
    return dict(nx=nx, ny=ny, nz=nz, nd=nd, lat=lat, lon=lon, xs=xs, ys=ys,
                zs=zs, flds=flds, trailing=trailing)


def main():
    d = sys.argv[1].rstrip('/')
    kind = sys.argv[2] if len(sys.argv) > 2 else ('sfc' if 'Surface' in d else '3d')
    names = RANGE_SFC if kind == 'sfc' else RANGE_3D
    paths = sorted(glob.glob(f'{d}/ERF_IC_*.bin'))
    print(f'=== {d}  ({kind}, {len(paths)} frames) ===')
    if not paths:
        raise SystemExit('no frames')

    bad = 0
    times, geom0 = [], None
    for p in paths:
        b = os.path.basename(p)
        stamp = b[len('ERF_IC_'):-len('.bin')]
        t = dt.datetime.strptime(stamp, '%Y_%m_%d_%H_%M')
        times.append(t)
        r = read(p)
        geom = (r['nx'], r['ny'], r['nz'], r['nd'],
                round(float(r['xs'][0]), 1), round(float(r['xs'][-1]), 1),
                round(float(r['ys'][0]), 1), round(float(r['ys'][-1]), 1))
        msgs = []
        if geom0 is None:
            geom0 = geom
            print(f'  geometry: nx={geom[0]} ny={geom[1]} nz={geom[2]} nfields={geom[3]}')
            print(f'            xs {geom[4]} .. {geom[5]}   ys {geom[6]} .. {geom[7]}')
            print(f'            lat {r["lat"].min():.3f}..{r["lat"].max():.3f}  '
                  f'lon {r["lon"].min():.3f}..{r["lon"].max():.3f}')
        elif geom != geom0:
            msgs.append(f'GEOMETRY CHANGED {geom} vs {geom0}')
        if r['trailing'] != 0:
            msgs.append(f'{r["trailing"]} trailing bytes (field count wrong?)')
        if len(r['flds']) != len(names):
            msgs.append(f'{len(r["flds"])} fields, audit knows {len(names)}')
        for (nm, lo, hi), a in zip(names, r['flds']):
            if a.size == 0:
                msgs.append(f'{nm}: EMPTY (short read)'); continue
            nbad = int(np.sum(~np.isfinite(a)))
            if nbad:
                msgs.append(f'{nm}: {nbad} non-finite')
            if a.min() < lo or a.max() > hi:
                msgs.append(f'{nm}: {a.min():.4g}..{a.max():.4g} outside [{lo},{hi}]')
        if msgs:
            bad += 1
            print(f'  FAIL {b}')
            for m in msgs:
                print(f'         {m}')

    print(f'\n  content: {len(paths)-bad} clean, {bad} failing')

    print('\n  cadence (ERF indexes frames POSITIONALLY -- a gap shifts everything):')
    gaps = 0
    for a, b in zip(times, times[1:]):
        h = (b - a).total_seconds() / 3600.
        if abs(h - 3.0) > 1e-6:
            gaps += 1
            print(f'    GAP {a:%m-%d %HZ} -> {b:%m-%d %HZ} = {h:.2f} h')
    span = (times[-1] - times[0]).total_seconds() / 3600.
    print(f'    {times[0]:%Y-%m-%d %HZ} -> {times[-1]:%Y-%m-%d %HZ}  '
          f'span {span:.0f} h, {gaps} gaps')

    print(f'\n  VERDICT: {"PASS" if (bad == 0 and gaps == 0) else "FAIL"}')


if __name__ == '__main__':
    main()
