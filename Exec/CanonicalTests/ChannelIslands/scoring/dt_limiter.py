#!/usr/bin/env python3
"""Where is the timestep set -- boundary cells or interior?

  dt_limiter.py <run_dir> [nplot]

The compressible dt is set by the max over cells of the advective CFL sum
  s = |u|/dx + |v|/dy + |w|/dz
and dz at the surface is ~18.7 m against dx = 3000 m, so a 1 m/s vertical
velocity in a surface cell outweighs a 150 m/s horizontal one. Finding the
argmax of s and asking whether it sits inside the 10-cell relaxation/NSCBC
band is what separates "NSCBC costs timestep at its own boundary" (a real
cost of the scheme) from "the interior is the limiter" (which is not).

real_width = 10 cells, so band = the outer 10 cells on each lateral face.
"""
import sys

import numpy as np
import yt

from plt_guard import plotfiles_by_time

yt.set_log_level(50)
NX, NY, DX = 192, 96, 3000.0
BAND = 10


def main():
    run = sys.argv[1].rstrip('/')
    npl = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    paths = plotfiles_by_time(f'/app/ERF/{run}')[-npl:]
    print(f'=== {run}: dt limiter location, last {len(paths)} plotfiles ===')
    print(f'    band = outer {BAND} cells (real_width); dz_min ~18.7 m vs dx 3000 m')

    for p in paths:
        ds = yt.load(p)
        t = float(ds.current_time)
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        u = np.abs(np.asarray(g[('boxlib', 'x_velocity')]))
        v = np.abs(np.asarray(g[('boxlib', 'y_velocity')]))
        w = np.abs(np.asarray(g[('boxlib', 'z_velocity')]))
        z = np.asarray(g[('boxlib', 'z_phys')])
        dz = np.gradient(z, axis=2)
        dz = np.maximum(np.abs(dz), 1e-3)

        su, sv, sw = u / DX, v / DX, w / dz
        s = su + sv + sw
        i, j, k = np.unravel_index(np.nanargmax(s), s.shape)

        dN, dE, dW, dS = NY - 1 - j, NX - 1 - i, i, j
        dmin = min(dN, dE, dW, dS)
        where = 'BAND' if dmin < BAND else 'INTERIOR'
        face = ['yhi', 'xhi', 'xlo', 'ylo'][int(np.argmin([dN, dE, dW, dS]))]

        # what fraction of the CFL sum is vertical, at the limiting cell
        tot = s[i, j, k]
        print(f'\n  t = {t:9.1f} s (h{t/3600:5.2f})   dt_cfl = {0.2/tot:.4f} s')
        print(f'    argmax at (i={i}, j={j}, k={k})  -> {where}, '
              f'{dmin} cells from {face}')
        print(f'    |u|={u[i,j,k]:7.2f} |v|={v[i,j,k]:7.2f} |w|={w[i,j,k]:7.3f} m/s   '
              f'dz={dz[i,j,k]:6.2f} m')
        print(f'    contributions: horiz {100*(su[i,j,k]+sv[i,j,k])/tot:5.1f}%   '
              f'VERT {100*sw[i,j,k]/tot:5.1f}%')

        # how much of the top-100 limiting cells are in the band
        flat = s.ravel()
        top = np.argpartition(flat, -100)[-100:]
        ii, jj, kk = np.unravel_index(top, s.shape)
        d = np.minimum.reduce([NY - 1 - jj, NX - 1 - ii, ii, jj])
        nb = int(np.sum(d < BAND))
        print(f'    of the 100 most CFL-limiting cells: {nb} in band, '
              f'{100-nb} interior;  k range {kk.min()}-{kk.max()}')


if __name__ == '__main__':
    main()
