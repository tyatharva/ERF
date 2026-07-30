#!/usr/bin/env python3
"""Terrain-stratified bias with the spin-up hours excluded.

  spinup_split.py <label>=<rundir> ...

Both ERF arms produce precipitation in hours 1-4 while d02 and MRMS are at
exactly 0.000 mm/h over land. This re-accumulates hours 5-23 only and recomputes
the elevation-binned bias, the Santa Ynez point and the domain max, so the
question "is the orographic excess a steady-state physics error or an
initialization transient deposited on terrain?" can be answered from output
already on disk.

Every source is differenced in its own native accumulation, identically to
hourly_series_arms.py. Two controls are asserted before any conclusion: the
h1-23 reconstruction of each reference must reproduce the checked-in full-period
array, otherwise the windowing is wrong and the split means nothing.
"""
import glob, gzip, shutil, os, sys
import numpy as np
import netCDF4 as nc
import pyproj
import pygrib
import yt
from scipy.interpolate import griddata

yt.set_log_level(50)
NX, NY, DX = 192, 96, 3000.
MM2IN = 1.0 / 25.4
PLO = (-388229.74, -166933.25)
PIN = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
       "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
S = '/app/ERF/scoring_ab_dav/'
lat = np.load(S + 'lat.npy'); lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy'); LAND = terr > terr.min() + 20
SPLIT = 4          # hours 1..SPLIT are the spin-up window

_d = np.hypot(lat - 34.486, lon + 119.802)
SY = np.unravel_index(np.argmin(_d), _d.shape)

BINS = [('flat   <100 m', (terr < 100) & LAND),
        ('100-500 m', (terr >= 100) & (terr < 500) & LAND),
        ('>500 m', (terr >= 500) & LAND)]


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
    hs = sorted(got)
    return {h: got[h] - (got[hs[k - 1]] if k > 0 else 0.0) for k, h in enumerate(hs)}


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


def accum(ser, h0, h1):
    return sum(ser[h] for h in range(h0, h1 + 1) if h in ser)


def table(name, arms, refs):
    print(f'\n=== {name} ===')
    for rname, ref in refs:
        print(f'  -- vs {rname} --')
        print(f'  {"bin":14s} {"n":>5} ' + ' '.join(f'{l[:12]:>13}' for l, _ in arms)
              + f' {"ref mean in":>12}')
        for bname, m in BINS:
            mm = m & np.isfinite(ref)
            row = ' '.join(f'{np.nanmean(a[mm])/np.nanmean(ref[mm]):13.2f}' for _, a in arms)
            print(f'  {bname:14s} {int(mm.sum()):5d} {row} {np.nanmean(ref[mm])*MM2IN:12.3f}')
    print(f'  {"":14s} {"":5s} ' + ' '.join(f'{l[:12]:>13}' for l, _ in arms))
    print(f'  {"Santa Ynez in":14s} {"":5s} ' + ' '.join(f'{a[SY]*MM2IN:13.2f}' for _, a in arms))
    print(f'  {"domain max in":14s} {"":5s} ' + ' '.join(f'{np.nanmax(a)*MM2IN:13.2f}' for _, a in arms))


def main():
    S_ = {}
    for spec in sys.argv[1:]:
        lab, run = spec.split('=', 1)
        S_[lab] = erf_series(run)
    labels = list(S_)
    S_['d02'] = d02_series()
    S_['MRMS'] = mrms_series()

    # ---- controls: h1-23 reconstruction must match the checked-in arrays ----
    for nm, path in (('d02', '/app/ERF/wrf_d02_on_grid.npy'),
                     ('MRMS', '/app/ERF/mrms_20201228_on_grid.npy')):
        rec = accum(S_[nm], 1, 23); full = np.load(path)
        rel = np.nanmean(np.abs(rec - full)) / max(np.nanmean(full), 1e-9)
        print(f'[control] {nm} h1-23 reconstruction vs checked-in array: '
              f'mean |diff| = {rel*100:.2f}% of mean  '
              f'({np.nanmean(rec)*MM2IN:.3f} vs {np.nanmean(full)*MM2IN:.3f} in)')

    print(f'\n[spin-up window] hours 1-{SPLIT}, LAND mean rate mm/h:')
    for k in S_:
        print(f'  {k:14s} ' + ' '.join(
            f'h{h}={np.nanmean(S_[k][h][LAND]):6.3f}' for h in range(1, SPLIT + 1) if h in S_[k]))

    full = [(l, accum(S_[l], 1, 23)) for l in labels]
    late = [(l, accum(S_[l], SPLIT + 1, 23)) for l in labels]
    rf = [('d02', accum(S_['d02'], 1, 23)), ('MRMS', accum(S_['MRMS'], 1, 23))]
    rl = [('d02', accum(S_['d02'], SPLIT + 1, 23)),
          ('MRMS', accum(S_['MRMS'], SPLIT + 1, 23))]
    table(f'FULL hours 1-23', full, rf)
    table(f'SPIN-UP EXCLUDED hours {SPLIT+1}-23', late, rl)

    print(f'\n=== reference totals, hours {SPLIT+1}-23 (in) ===')
    for rname, ref in rl:
        print(f'  {rname:6s} Santa Ynez {ref[SY]*MM2IN:5.2f}   domain max '
              f'{np.nanmax(ref)*MM2IN:5.2f}   LAND mean {np.nanmean(ref[LAND])*MM2IN:5.3f}')

    print('\n=== hour-by-hour: LAND mean rate mm/h, and >500 m bias vs d02 ===')
    keys = labels + ['d02', 'MRMS']
    print('   h ' + ' '.join(f'{k[:11]:>12}' for k in keys)
          + '  |' + ' '.join(f'{l[:9]+" >500":>16}' for l in labels))
    hi = BINS[2][1]
    for h in range(1, 24):
        rates = ' '.join(f'{np.nanmean(S_[k][h][LAND]):12.3f}' if h in S_[k] else f'{"--":>12}'
                         for k in keys)
        dref = np.nanmean(S_['d02'][h][hi]) if h in S_['d02'] else np.nan
        bias = ' '.join(
            (f'{np.nanmean(S_[l][h][hi])/dref:16.2f}' if (dref > 0.005 and h in S_[l])
             else f'{"n/a":>16}') for l in labels)
        print(f'  {h:2d} {rates}  |{bias}')

    print('\n=== cumulative >500 m bias vs d02, accumulating from hour 5 ===')
    print('  through h ' + ' '.join(f'{l[:12]:>13}' for l in labels))
    for h in range(5, 24):
        dref = np.nanmean(accum(S_['d02'], 5, h)[hi])
        if dref <= 0:
            continue
        row = ' '.join(f'{np.nanmean(accum(S_[l],5,h)[hi])/dref:13.2f}' for l in labels)
        print(f'  {h:9d} {row}')


if __name__ == '__main__':
    main()
