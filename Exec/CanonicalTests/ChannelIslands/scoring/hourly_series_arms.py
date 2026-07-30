#!/usr/bin/env python3
"""Hourly domain-mean precipitation RATE for arbitrary arms, 2020-12-28.

  hourly_series_arms.py <label>=<rundir> ...

Generalises hourly_series_c404.py, which hardwired the Davies/NSCBC pair. Each
series is differenced in its own native accumulation: ERF rain_accum between
consecutive hourly plotfiles, d02 RAINNC+bucket between consecutive stamps, MRMS
the per-file hourly total. Answers whether a difference between arms is present
from hour 1 or develops with time.
"""
import glob, gzip, shutil, os, sys
import numpy as np
import netCDF4 as nc
import pyproj
import pygrib
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

yt.set_log_level(50)
OUT = '/app/ERF/figs'
os.makedirs(OUT, exist_ok=True)
NX, NY, DX = 192, 96, 3000.
PLO = (-388229.74, -166933.25)
PIN = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
       "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
S = '/app/ERF/scoring_ab_dav/'
lat = np.load(S + 'lat.npy'); lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy'); LAND = terr > terr.min() + 20
ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
D = np.minimum.reduce([ii, jj, NX - 1 - ii, NY - 1 - jj])


def erf_series(run):
    got = {}
    for p in sorted(glob.glob(f'/app/ERF/{run}/plt[0-9]*')):
        try:
            ds = yt.load(p)
        except Exception:
            continue
        t = float(ds.current_time); h = int(round(t / 3600.))
        if abs(t - h * 3600.) > 400. or h < 1 or h > 23 or h in got:
            continue
        g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        got[h] = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
    hs = sorted(got); out = {}
    for k, h in enumerate(hs):
        out[h] = got[h] - (got[hs[k - 1]] if k > 0 else 0.0)
    return out


def d02_series():
    f = nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
    B = float(f.BUCKET_MM)
    R = lambda k: (np.asarray(f.variables['RAINNC'][k])
                   + B * np.asarray(f.variables['I_RAINNC'][k]))
    la = np.asarray(f.variables['XLAT'][0]); lo = np.asarray(f.variables['XLONG'][0])
    x, y = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(PIN),
                                       always_xy=True).transform(lo, la)
    i = np.floor((x - PLO[0]) / DX).astype(int)
    j = np.floor((y - PLO[1]) / DX).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
    out = {}
    for h in range(1, 24):
        d = R(h) - R(h - 1)
        s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
        np.add.at(s, (i[ok], j[ok]), d[ok]); np.add.at(c, (i[ok], j[ok]), 1.)
        out[h] = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    return out


def mrms_series():
    out = {}; pts = None; sel = None
    for fz in sorted(glob.glob('/app/ERF/mrms/*.grib2.gz')):
        h = int(os.path.basename(fz).split('-')[1][0:2])
        with gzip.open(fz, 'rb') as a, open('/tmp/m.g2', 'wb') as b:
            shutil.copyfileobj(a, b)
        g = pygrib.open('/tmp/m.g2'); m = g.message(1)
        v = m.values
        v = v.filled(np.nan) if hasattr(v, 'filled') else np.asarray(v, float)
        v = np.asarray(v, float); v[v < 0] = np.nan
        if pts is None:
            mla, mlo = m.latlons(); mlo = np.where(mlo > 180, mlo - 360, mlo)
            sel = ((mla >= lat.min() - .1) & (mla <= lat.max() + .1)
                   & (mlo >= lon.min() - .1) & (mlo <= lon.max() + .1))
            pts = np.column_stack([mlo[sel], mla[sel]])
        out[h] = griddata(pts, np.nan_to_num(v[sel]), (lon, lat), method='nearest')
        g.close()
    return out


def main():
    S_ = {}
    for spec in sys.argv[1:]:
        lab, run = spec.split('=', 1)
        S_[lab] = erf_series(run)
    S_['d02'] = d02_series()
    S_['MRMS'] = mrms_series()

    masks = (('full domain', np.ones_like(LAND, bool)),
             ('interior d>=20', D >= 20), ('LAND', LAND))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    for c, (mn, mk) in enumerate(masks):
        for name, ser in S_.items():
            hs = sorted(ser)
            ref = name in ('d02', 'MRMS')
            ax[c].plot(hs, [np.nanmean(ser[h][mk]) for h in hs], label=name,
                       lw=2.2 if ref else 1.6, ls='-' if ref else '--',
                       marker='o', ms=3)
        ax[c].set_xlabel('hour (UTC) of 2020-12-28'); ax[c].set_title(mn, fontsize=10)
        ax[c].grid(alpha=.3); ax[c].legend(fontsize=8)
        if c == 0:
            ax[c].set_ylabel('domain-mean rate, mm/h')
    p = f'{OUT}/wsm6_hourly_rate.png'
    fig.savefig(p, dpi=150)
    print('wrote', p)

    for mn, mk in masks:
        print(f'\n-- {mn}: hourly mean rate mm/h --')
        print('  h  ' + ' '.join(f'{k:>12}' for k in S_))
        for h in range(1, 24):
            print(f' {h:>2}  ' + ' '.join(
                f'{np.nanmean(S_[k][h][mk]):12.3f}' if h in S_[k] else f'{"--":>12}'
                for k in S_))


if __name__ == '__main__':
    main()
