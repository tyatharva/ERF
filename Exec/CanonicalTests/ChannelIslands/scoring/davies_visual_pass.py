#!/usr/bin/env python3
"""Quick visual + metrics pass on the Davies control, h48-h71.

  davies_visual_pass.py <run_dir> <outdir>

Not the end-analysis suite: three figures and one metrics table.

MASK. MRMS has no reliable over-ocean QPE, so every comparison here is
LAND-ONLY: terrain > 20 m AND finite in all three sources. The same mask is
applied to Davies, d02 and MRMS so the ratios are over identical cells; the
cell count is printed and captioned.

REGRIDDING. Both references are bin-averaged (area-averaged) onto the 3 km
Domain A grid from a FINER source -- d02 at 1.5 km, MRMS at ~1 km, median 9.0
source cells per target per hour. That is the well-posed direction. No
nearest-neighbour sampling is used anywhere, which matters because NN would
preserve 1 km extremes and inflate the reference maxima relative to a 3 km
model field.

WINDOW. h48-h71 = 2020-12-28 00Z-23Z, 23 h. ERF window accumulation is
rain_accum(h71) - rain_accum(h48); MRMS is the sum of the 23 hourly QPE files
valid 01Z..23Z; d02 is RAINNC+I_RAINNC differenced over the same span.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyproj
import yt

from plt_guard import is_poisoned, plotfiles_by_time

yt.set_log_level(50)
REFS = '/app/ERF/refs'
NX, NY, DX = 192, 96, 3000.0
OCEAN_Z, BAND = 12.5, 10
OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
PROB_LO = (-389768.39, -89294.60)
H0, H1 = 48, 71


def erf_hourly(run):
    """rain_accum at each hour h48..h71, plus the terrain."""
    want = {h: None for h in range(H0, H1 + 1)}
    zs = None
    for p in plotfiles_by_time(f'/app/ERF/{run}'):
        ds = yt.load(p)
        h = float(ds.current_time) / 3600.0
        hr = int(round(h))
        if hr in want and abs(h - hr) < 0.02 and want[hr] is None:
            g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
            ra = np.asarray(g[('boxlib', 'rain_accum')])[:, :, 0]
            if is_poisoned(ra):
                raise SystemExit(f'{run} h{hr}: poisoned plotfile')
            want[hr] = ra
            if zs is None:
                zs = np.asarray(g[('boxlib', 'z_phys')])[:, :, 0]
    missing = [h for h, v in want.items() if v is None]
    if missing:
        raise SystemExit(f'missing hours: {missing}')
    return want, zs


def d02_hourly():
    f = nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
    la = np.asarray(f.variables['XLAT'][0])
    lo = np.asarray(f.variables['XLONG'][0])
    B = 100.0 if 'I_RAINNC' in f.variables else 0.0
    def R(k):
        v = np.asarray(f.variables['RAINNC'][k], dtype=float)
        if B:
            v = v + B * np.asarray(f.variables['I_RAINNC'][k], dtype=float)
        return v
    x, y = pyproj.Transformer.from_crs(
        4326, pyproj.CRS.from_proj4(OURS), always_xy=True).transform(lo, la)
    i = np.floor((x - PROB_LO[0]) / DX).astype(int)
    j = np.floor((y - PROB_LO[1]) / DX).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
    ii, jj = i[ok], j[ok]

    def grid(v):
        s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
        np.add.at(s, (ii, jj), v[ok]); np.add.at(c, (ii, jj), 1.0)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)

    out = [grid(R(k) - R(k - 1)) for k in range(1, 24)]
    f.close()
    return out


def mrms_hourly():
    import glob, gzip, shutil, tempfile
    import pygrib
    tr = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(OURS),
                                     always_xy=True)
    out = []
    for fp in sorted(glob.glob('/app/ERF/mrms/*.grib2.gz')):
        with tempfile.NamedTemporaryFile(suffix='.grib2', delete=False) as tf:
            with gzip.open(fp, 'rb') as gz:
                shutil.copyfileobj(gz, tf)
            tmp = tf.name
        try:
            g = pygrib.open(tmp); m = g[1]
            vals = np.asarray(m.values, dtype=float)
            la, lo = m.latlons(); g.close()
        finally:
            os.unlink(tmp)
        lo = np.where(lo > 180, lo - 360, lo)
        box = (np.isfinite(vals) & (vals >= 0) & (la > 31) & (la < 38)
               & (lo > -126) & (lo < -114))
        x, y = tr.transform(lo[box], la[box])
        i = np.floor((x - PROB_LO[0]) / DX).astype(int)
        j = np.floor((y - PROB_LO[1]) / DX).astype(int)
        ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
        s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
        np.add.at(s, (i[ok], j[ok]), vals[box][ok])
        np.add.at(c, (i[ok], j[ok]), 1.0)
        out.append(np.where(c > 0, s / np.maximum(c, 1), np.nan))
    return out


def main():
    run = sys.argv[1].rstrip('/')
    out = sys.argv[2].rstrip('/')
    os.makedirs(out, exist_ok=True)

    acc, zs = erf_hourly(run)
    erf_tot = np.maximum(acc[H1] - acc[H0], 0.0)
    erf_hr = [np.maximum(acc[h + 1] - acc[h], 0.0) for h in range(H0, H1)]
    d02_hr = d02_hourly()
    mrms_hr = mrms_hourly()
    d02_tot = np.nansum(np.dstack(d02_hr), axis=2)
    mrms_tot = np.nansum(np.dstack(mrms_hr), axis=2)

    terr = zs - OCEAN_Z
    M = ((terr > 20.0) & np.isfinite(d02_tot) & np.isfinite(mrms_tot)
         & np.isfinite(erf_tot))
    n = int(M.sum())
    CAP = (f'Land-only mask: terrain > 20 m AND finite in all three sources. '
           f'n = {n} of {NX*NY} cells ({100.0*n/(NX*NY):.1f}%). '
           f'References bin-averaged (area-averaged) from finer grids: '
           f'd02 1.5 km, MRMS ~1 km (median 9.0 src cells/target/hour).')

    ii, jj = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    dN = NY - 1 - jj

    # ---------- FIG 1: maps ----------
    vmax = float(np.nanpercentile(np.concatenate(
        [erf_tot[M], d02_tot[M], mrms_tot[M]]), 99))
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    show = lambda a: np.where(M, a, np.nan).T
    for k, (a, t) in enumerate([(erf_tot, 'ERF Davies'), (d02_tot, 'WRF d02'),
                                (mrms_tot, 'MRMS')]):
        p = ax.flat[k].imshow(show(a), origin='lower', vmin=0, vmax=vmax,
                              cmap='viridis', aspect='auto')
        ax.flat[k].set_title(f'{t}  (mean {np.nanmean(a[M]):.1f} mm)')
        plt.colorbar(p, ax=ax.flat[k], label='mm / 23 h')
    d = erf_tot - d02_tot
    lim = float(np.nanpercentile(np.abs(d[M]), 99))
    p = ax.flat[3].imshow(show(d), origin='lower', vmin=-lim, vmax=lim,
                          cmap='RdBu_r', aspect='auto')
    ax.flat[3].set_title(f'Davies - d02  (mean {np.nanmean(d[M]):+.1f} mm)')
    plt.colorbar(p, ax=ax.flat[3], label='mm / 23 h')
    for a in ax.flat:
        a.axhline(NY - 1 - BAND, color='w', ls='--', lw=1)
        a.set_xlabel('i'); a.set_ylabel('j')
    fig.suptitle('Davies control, h48-h71 = 2020-12-28 00Z-23Z. '
                 'Dashed = 10-cell band edge at the north wall.', fontsize=10)
    fig.text(0.5, 0.005, CAP, ha='center', fontsize=7.5, wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(f'{out}/01_maps_h48_h71.png', dpi=130)
    plt.close(fig)

    # ---------- FIG 2: wall-distance profile ----------
    fig, a2 = plt.subplots(figsize=(9, 5))
    ds_ = np.arange(0, 60)
    for arr, lab, c in ((erf_tot, 'ERF Davies', 'C3'), (d02_tot, 'WRF d02', 'C0'),
                        (mrms_tot, 'MRMS', 'k')):
        prof = [np.nanmean(arr[M & (dN == q)]) if (M & (dN == q)).any() else np.nan
                for q in ds_]
        a2.plot(ds_, prof, label=lab, color=c, lw=1.8)
    a2.axvline(BAND, color='grey', ls='--', label=f'band edge (dN={BAND})')
    a2.set_xlabel('dN  (cells from the NORTH wall)')
    a2.set_ylabel('mean accumulation, mm / 23 h')
    a2.set_title('Wall-distance profile, land-only, h48-h71')
    a2.legend(); a2.grid(alpha=0.3)
    fig.text(0.5, 0.005, CAP, ha='center', fontsize=7.5, wrap=True)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(f'{out}/02_wall_profile.png', dpi=130)
    plt.close(fig)

    # ---------- FIG 3: hourly series ----------
    hrs = np.arange(1, 24)
    se = [np.nanmean(a[M]) for a in erf_hr]
    sd = [np.nanmean(a[M]) for a in d02_hr]
    sm = [np.nanmean(a[M]) for a in mrms_hr]
    fig, a3 = plt.subplots(figsize=(9, 5))
    a3.plot(hrs, se, 'C3-o', ms=3, label=f'ERF Davies (sum {np.sum(se):.1f})')
    a3.plot(hrs, sd, 'C0-o', ms=3, label=f'WRF d02 (sum {np.sum(sd):.1f})')
    a3.plot(hrs, sm, 'k-o', ms=3, label=f'MRMS (sum {np.sum(sm):.1f})')
    a3.set_xlabel('hour ending, 2020-12-28 UTC  (= ERF h48+n)')
    a3.set_ylabel('domain-mean hourly precip, mm')
    a3.set_title('Hourly domain-mean, land-only mask, h48-h71')
    a3.legend(); a3.grid(alpha=0.3)
    fig.text(0.5, 0.005, CAP, ha='center', fontsize=7.5, wrap=True)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(f'{out}/03_hourly_series.png', dpi=130)
    plt.close(fig)

    # ---------- metrics ----------
    L = []
    P = L.append
    P('# Davies control, h48-h71 (2020-12-28 00Z-23Z) -- metrics\n')
    P(f'Mask: land-only, terrain > 20 m AND finite in all three sources.')
    P(f'n = {n} of {NX*NY} cells ({100.0*n/(NX*NY):.1f}%).')
    P('Regridding: both references BIN-AVERAGED (area-averaged) from finer')
    P('grids -- d02 1.5 km, MRMS ~1 km, median 9.0 source cells per target')
    P('per hour. No nearest-neighbour sampling. NN would retain 1 km extrema')
    P('and inflate reference maxima against a 3 km model field.')
    P('Restart discontinuity: the arm rolled back to the h8 checkpoint once;')
    P('the scored window lies entirely inside the surviving leg and window')
    P('accumulation had zero negative cells.\n')

    e, d0, m0 = erf_tot[M], d02_tot[M], mrms_tot[M]
    P('| metric | Davies vs d02 | Davies vs MRMS |')
    P('|---|---|---|')
    P(f'| ERF domain mean (mm) | {e.mean():.3f} | {e.mean():.3f} |')
    P(f'| reference mean (mm) | {d0.mean():.3f} | {m0.mean():.3f} |')
    P(f'| ratio | {e.mean()/d0.mean():.2f}x | {e.mean()/m0.mean():.2f}x |')
    P(f'| bias (mm) | {(e-d0).mean():+.3f} | {(e-m0).mean():+.3f} |')
    P(f'| RMSE (mm) | {np.sqrt(((e-d0)**2).mean()):.3f} | '
      f'{np.sqrt(((e-m0)**2).mean()):.3f} |')
    P(f'| correlation | {np.corrcoef(e,d0)[0,1]:.3f} | '
      f'{np.corrcoef(e,m0)[0,1]:.3f} |')

    def nearfar(a):
        nr = M & (dN <= 12); fr = M & (dN >= 21)
        return a[nr].mean() / a[fr].mean()
    P(f'| near/far (dN<=12 / dN>=21) | ERF {nearfar(erf_tot):.2f}, '
      f'd02 {nearfar(d02_tot):.2f} | MRMS {nearfar(mrms_tot):.2f} |')
    P('')

    P('## Terrain-stratified ERF/reference, dN 0-6\n')
    P('| terrain | n | ERF mm | d02 mm | ERF/d02 | MRMS mm | ERF/MRMS |')
    P('|---|---|---|---|---|---|---|')
    for lo, hi in [(20, 100), (100, 300), (300, 600), (600, 1000), (1000, 3000)]:
        q = M & (terr >= lo) & (terr < hi) & (dN <= 6)
        if q.sum() < 12:
            continue
        P(f'| {lo}-{hi} m | {q.sum()} | {erf_tot[q].mean():.1f} | '
          f'{d02_tot[q].mean():.1f} | {erf_tot[q].mean()/d02_tot[q].mean():.2f}x | '
          f'{mrms_tot[q].mean():.1f} | {erf_tot[q].mean()/mrms_tot[q].mean():.2f}x |')

    txt = '\n'.join(L)
    open(f'{out}/metrics.md', 'w').write(txt + '\n')
    print(txt)
    print(f'\nwrote {out}/')


if __name__ == '__main__':
    main()
