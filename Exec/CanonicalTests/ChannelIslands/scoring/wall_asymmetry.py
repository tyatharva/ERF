#!/usr/bin/env python3
"""Why is the north wall 3.79x and the east wall 1.27x? Same code, same BC.

  wall_asymmetry.py <run> <x_lo> <y_lo> <x_hi> <y_hi> <frame_dir>

Post-processing only. Four questions, in order of how cheaply they could kill the
"scheme" explanation:

  1. INFLOW vs OUTFLOW.   Hourly wall-normal MASS flux at all four faces. If the
     north face is inflowing while the east is outflowing through the h13-h17
     enhancement window, the asymmetry is just which faces admit air.
  2. FRAME FIDELITY.      Hourly ERF vs CONUS404 vapour flux, per face. Does the
     north face over-supply where the east does not?
  3. TERRAIN OUTSIDE.     What sits immediately beyond each wall. If the frames
     arriving at the north face already carry orographic precipitation generated
     outside the domain, that is a DATA problem, not a scheme problem.
  4. WALL vs COAST.       At the north face, where the two diverge.

Everything is inward-positive, integrated 0-6000 m ASL, matching inflow_flux.py.
"""
import glob
import os
import sys

import numpy as np
import yt

from plt_guard import plotfiles_by_time

yt.set_log_level(50)
NX, NY, DX = 192, 96, 3000.0
ZTOP = 6000.0


def erf_hourly(run, faces_out=None):
    out = {}
    for p in plotfiles_by_time(f'/app/ERF/{run}'):
        ds = yt.load(p)
        t = float(ds.current_time); h = int(round(t / 3600.))
        if abs(t - h * 3600.) > 400. or h < 1 or h > 23 or h in out:
            continue
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        rho = np.asarray(g[('boxlib', 'density')])
        qv = np.asarray(g[('boxlib', 'qv')])
        u = np.asarray(g[('boxlib', 'x_velocity')])
        v = np.asarray(g[('boxlib', 'y_velocity')])
        z = np.asarray(g[('boxlib', 'z_phys')])
        dz = np.gradient(z, axis=2)
        below = z <= ZTOP
        wt = np.where(below, rho * dz, 0.0)          # mass weight
        wq = wt * qv                                  # vapour weight
        M = {'xlo': (wt * u)[0].sum() * DX, 'xhi': -(wt * u)[-1].sum() * DX,
             'ylo': (wt * v)[:, 0].sum() * DX, 'yhi': -(wt * v)[:, -1].sum() * DX}
        Q = {'xlo': (wq * u)[0].sum() * DX, 'xhi': -(wq * u)[-1].sum() * DX,
             'ylo': (wq * v)[:, 0].sum() * DX, 'yhi': -(wq * v)[:, -1].sum() * DX}
        out[h] = dict(M=M, Q=Q)
    return out


def frame_flux(path, xl, yl, xh, yh):
    with open(path, 'rb') as f:
        nx, ny, nz, nd = np.fromfile(f, dtype='<i4', count=4)
        f.read(4 * 2 * nx * ny)
        xs = np.fromfile(f, dtype='<f4', count=nx).astype(float)
        ys = np.fromfile(f, dtype='<f4', count=ny).astype(float)
        zs = np.fromfile(f, dtype='<f4', count=nz).astype(float)
        blk = lambda: np.fromfile(f, dtype='<f4',
                                  count=nx * ny * nz).reshape(nz, ny, nx).astype(float)
        rho = blk(); u = blk(); v = blk(); _w = blk(); _th = blk(); qv = blk()
    dzs = np.gradient(zs); keep = zs <= ZTOP

    def samp(a, xq, yq):
        fi = np.clip((np.asarray(xq) - xs[0]) / (xs[1] - xs[0]), 0, nx - 1.001)
        fj = np.clip((np.asarray(yq) - ys[0]) / (ys[1] - ys[0]), 0, ny - 1.001)
        i0 = fi.astype(int); j0 = fj.astype(int); ti = fi - i0; tj = fj - j0
        return (a[:, j0, i0] * (1 - ti) * (1 - tj) + a[:, j0, i0 + 1] * ti * (1 - tj)
                + a[:, j0 + 1, i0] * (1 - ti) * tj + a[:, j0 + 1, i0 + 1] * ti * tj)

    yc = yl + (np.arange(NY) + 0.5) * DX
    xc = xl + (np.arange(NX) + 0.5) * DX

    def face(xq, yq, comp, moist):
        r = samp(rho, xq, yq); c = samp(u if comp == 'u' else v, xq, yq)
        wt = (dzs * keep)[:, None] * r
        if moist:
            wt = wt * samp(qv, xq, yq)
        return (wt * c).sum() * DX

    mk = lambda moist: {
        'xlo': face(np.full(NY, xl), yc, 'u', moist),
        'xhi': -face(np.full(NY, xh), yc, 'u', moist),
        'ylo': face(xc, np.full(NX, yl), 'v', moist),
        'yhi': -face(xc, np.full(NX, yh), 'v', moist)}
    return dict(M=mk(False), Q=mk(True))


def main():
    run = sys.argv[1]
    xl, yl, xh, yh = [float(v) for v in sys.argv[2:6]]
    fdir = sys.argv[6]
    E = erf_hourly(run)
    FR = {}
    for p in sorted(glob.glob(f'{fdir}/ERF_IC_*.bin')):
        parts = os.path.basename(p).split('_')
        if parts[4] == '28':
            FR[int(parts[5])] = frame_flux(p, xl, yl, xh, yh)
    W = ['xlo', 'ylo', 'xhi', 'yhi']

    print(f'=== 1. WALL-NORMAL MASS FLUX, inward positive, 0-6000 m (kg/s) -- {run} ===')
    print('    which faces admit air, and when')
    print(f'  {"h":>3}' + ''.join(f'{w:>12}' for w in W) + f'{"  sign pattern":>18}')
    for h in sorted(E):
        M = E[h]['M']
        sign = ''.join(('IN ' if M[w] > 0 else 'out') for w in W)
        print(f'  {h:>3}' + ''.join(f'{M[w]:12.3e}' for w in W) + f'{sign:>18}')

    print(f'\n=== 2. VAPOUR FLUX, ERF vs CONUS404 frame, per face (kg/s) ===')
    print('    ratio > 1 = ERF admits more than the driver through that face')
    print(f'  {"h":>3} {"src":>9}' + ''.join(f'{w:>12}' for w in W))
    for h in sorted(FR):
        if h not in E or h < 1 or h > 23:
            continue
        print(f'  {h:>3} {"ERF":>9}' + ''.join(f'{E[h]["Q"][w]:12.3e}' for w in W))
        print(f'  {"":>3} {"CONUS404":>9}' + ''.join(f'{FR[h]["Q"][w]:12.3e}' for w in W))
        print(f'  {"":>3} {"ratio":>9}' + ''.join(
            f'{E[h]["Q"][w]/FR[h]["Q"][w]:12.2f}' if abs(FR[h]['Q'][w]) > 1e3
            else f'{"--":>12}' for w in W))

    print(f'\n=== 2b. MEAN |ratio| over frame hours, per face ===')
    for w in W:
        rs = [E[h]['Q'][w] / FR[h]['Q'][w] for h in sorted(FR)
              if h in E and 1 <= h <= 23 and abs(FR[h]['Q'][w]) > 1e3]
        if rs:
            print(f'  {w:>5}: n={len(rs)}  mean {np.mean(rs):6.2f}  '
                  f'median {np.median(rs):6.2f}  range {min(rs):6.2f}..{max(rs):6.2f}')


if __name__ == '__main__':
    main()
