#!/usr/bin/env python3
"""Interior wind and moisture: ERF arms vs WRF d02, on the ERF grid.

  interior_wind_check.py <d02_hour> <label>=<rundir> [...]

Item 63 measured the ERF interior running 2.5-3.6x too fast with Cd 9-17x below
physical. Today's land-mask and roughness fixes changed surface drag over
terrain by 2-3 orders of magnitude, and that has never been re-measured. The
precipitation is downstream of this: if the interior flow is wrong, no
microphysics setting can put rain in the right place.

Compares at two heights AGL (500 m, inside the drag-affected layer, and 3000 m,
above it) so a surface-drag error and a synoptic error can be told apart:
a drag error should shrink with height, a synoptic error should not.

d02 is destaggered and bin-averaged onto the ERF grid, the same convention as
every other reference in this campaign.
"""
import sys

import numpy as np
import netCDF4 as nc
import pyproj
import yt

from plt_guard import plotfiles_by_time

yt.set_log_level(50)

OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
PROB_LO = (-389768.39, -89294.60)
NX, NY, DX, BAND = 192, 96, 3000.0, 10
WRFOUT = '/app/ERF/wrfout_d02_2020-12-28_00_00_00'
LEVELS = (500.0, 3000.0)


def erf_at(rundir, hour):
    """u, v, qv interpolated to LEVELS above local terrain, on the ERF grid."""
    for p in plotfiles_by_time(rundir):
        ds = yt.load(p)
        if abs(float(ds.current_time) / 3600.0 - hour) < 0.02:
            g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
            u = np.asarray(g[('boxlib', 'x_velocity')])
            v = np.asarray(g[('boxlib', 'y_velocity')])
            qv = np.asarray(g[('boxlib', 'qv')])
            z = np.asarray(g[('boxlib', 'z_phys')])
            zs = z[:, :, 0:1]                       # terrain height
            agl = z - zs
            out = {}
            for L in LEVELS:
                k = np.abs(agl - L).argmin(axis=2)
                I, J = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
                out[L] = (u[I, J, k], v[I, J, k])
            # column vapour, mass-weighted proxy: sum qv*dz
            dz = np.diff(z, axis=2, prepend=z[:, :, 0:1])
            out['cwv'] = (qv * dz).sum(axis=2)
            return out
    raise SystemExit(f'{rundir}: no plotfile at h{hour}')


def d02_at(hour):
    f = nc.Dataset(WRFOUT)
    U = np.asarray(f['U'][hour]); V = np.asarray(f['V'][hour])
    u = 0.5 * (U[:, :, :-1] + U[:, :, 1:])          # destagger x
    v = 0.5 * (V[:, :-1, :] + V[:, 1:, :])          # destagger y
    ph = np.asarray(f['PH'][hour]); phb = np.asarray(f['PHB'][hour])
    zf = (ph + phb) / 9.81
    zc = 0.5 * (zf[:-1] + zf[1:])
    agl = zc - np.asarray(f['HGT'][hour])[None, :, :]
    qv = np.asarray(f['QVAPOR'][hour])
    la = np.asarray(f['XLAT'][0]); lo = np.asarray(f['XLONG'][0])
    tr = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(OURS), always_xy=True)
    x, y = tr.transform(lo, la)
    i = np.floor((x - PROB_LO[0]) / DX).astype(int)
    j = np.floor((y - PROB_LO[1]) / DX).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)

    def binmean(field2d):
        s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
        np.add.at(s, (i[ok], j[ok]), field2d[ok])
        np.add.at(c, (i[ok], j[ok]), 1.0)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)

    out = {}
    for L in LEVELS:
        k = np.abs(agl - L).argmin(axis=0)
        J, I = np.meshgrid(np.arange(agl.shape[1]), np.arange(agl.shape[2]), indexing='ij')
        out[L] = (binmean(u[k, J, I]), binmean(v[k, J, I]))
    dzf = np.diff(zf, axis=0)
    out['cwv'] = binmean((qv * dzf).sum(axis=0))
    return out


def main():
    hour = int(sys.argv[1])
    arms = [a.split('=', 1) for a in sys.argv[2:]]
    d = d02_at(hour)
    terr = np.load('/app/ERF/refs/domA_terrain.npy')
    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dmin = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])
    M = (dmin >= BAND) & np.isfinite(d[LEVELS[0]][0])
    ML = M & (terr > 20.0)

    for L in LEVELS:
        du, dv = d[L]
        dsp = np.hypot(du, dv)
        print(f'\n=== {L:.0f} m AGL ===   (interior, {M.sum()} cells; land subset {ML.sum()})')
        print(f'{"arm":26s} {"|V| mean":>9s} {"ratio":>7s} {"corr u":>8s} {"corr v":>8s} '
              f'{"land |V|":>9s} {"land ratio":>11s}')
        print(f'{"WRF d02":26s} {np.nanmean(dsp[M]):9.2f} {1.0:7.2f} '
              f'{1.0:+8.3f} {1.0:+8.3f} {np.nanmean(dsp[ML]):9.2f} {1.0:11.2f}')
        for lab, rd in arms:
            e = erf_at(rd, 12.0)
            eu, ev = e[L]
            esp = np.hypot(eu, ev)
            print(f'{lab:26s} {np.nanmean(esp[M]):9.2f} '
                  f'{np.nanmean(esp[M])/np.nanmean(dsp[M]):7.2f} '
                  f'{np.corrcoef(eu[M], du[M])[0,1]:+8.3f} '
                  f'{np.corrcoef(ev[M], dv[M])[0,1]:+8.3f} '
                  f'{np.nanmean(esp[ML]):9.2f} '
                  f'{np.nanmean(esp[ML])/np.nanmean(dsp[ML]):11.2f}')

    print(f'\n=== column water vapour (kg/m2 proxy) ===')
    print(f'{"WRF d02":26s} {np.nanmean(d["cwv"][M]):9.2f}')
    for lab, rd in arms:
        e = erf_at(rd, 12.0)
        print(f'{lab:26s} {np.nanmean(e["cwv"][M]):9.2f} '
              f'ratio {np.nanmean(e["cwv"][M])/np.nanmean(d["cwv"][M]):5.2f}  '
              f'corr {np.corrcoef(e["cwv"][M], d["cwv"][M])[0,1]:+.3f}')


if __name__ == '__main__':
    main()
