#!/usr/bin/env python3
"""Hourly lateral moisture flux, low-level flow and column vapour over the ranges.

  inflow_flux.py <label>=<rundir> ...

Discriminator for item 50: the >500 m bias is built entirely in h13-h17, where
both references spin down and ERF does not. Either the lateral boundary keeps
supplying moisture after the driver stops, or the excess is generated internally.

Computes, hourly:
  A. lateral water-vapour flux through each face, integrated 0-6000 m ASL,
     for each ERF arm AND for the CONUS404 driving frames sampled on the SAME
     faces at the SAME heights (frames are 3-hourly, so 8 of the 24 hours).
  B. mass-weighted sub-1500 m AGL wind over the >500 m cells, ERF vs d02.
  C. column water vapour over the >500 m cells, ERF arms vs d02.

The ERF-vs-frame comparison in A is the discriminator. B and C say whether any
excess is arriving or being retained.
"""
import glob, os, sys
import numpy as np
import netCDF4 as nc
import pyproj
import yt

yt.set_log_level(50)
NX, NY, DX = 192, 96, 3000.0
PLOX, PLOY = -388229.75, -166933.25
G = 9.81
RD, CP, P0 = 287.0, 1004.5, 1.0e5
ZTOP_FLUX = 6000.0          # integrate the flux over 0-6000 m ASL
PIN = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
       "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")

S = '/app/ERF/scoring_ab_dav/'
lat = np.load(S + 'lat.npy'); lon = np.load(S + 'lon.npy')
terr = np.load(S + 'terrain.npy')
LAND = terr > terr.min() + 20
HI = (terr >= 500) & LAND


# ----------------------------------------------------------------- ERF -----
def erf_hourly(run):
    out = {}
    for p in sorted(glob.glob(f'/app/ERF/{run}/plt[0-9]*')):
        try:
            ds = yt.load(p)
        except Exception:
            continue
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
        w = rho * qv * dz
        below = z <= ZTOP_FLUX
        fx = np.where(below, w * u, 0.0)
        fy = np.where(below, w * v, 0.0)
        # inward-positive through each lateral face
        F = {'xlo': fx[0, :, :].sum() * DX, 'xhi': -fx[-1, :, :].sum() * DX,
             'ylo': fy[:, 0, :].sum() * DX, 'yhi': -fy[:, -1, :].sum() * DX}
        agl = z - z[:, :, :1]
        m = (agl < 1500.0) & HI[:, :, None]
        ww = np.where(m, rho, 0.0); sw = max(ww.sum(), 1e-9)
        pw = (rho * qv * dz).sum(axis=2)
        out[h] = dict(F=F, ub=(u * ww).sum() / sw, vb=(v * ww).sum() / sw,
                      pw=float(np.nanmean(pw[HI])))
    return out


# -------------------------------------------------------------- frames -----
def frame_flux(path):
    """Lateral vapour flux from a CONUS404 .bin frame, on the ERF faces."""
    with open(path, 'rb') as f:
        nx, ny, nz, nd = np.fromfile(f, dtype='<i4', count=4)
        f.read(4 * 2 * nx * ny)                       # lat, lon
        xs = np.fromfile(f, dtype='<f4', count=nx).astype(float)
        ys = np.fromfile(f, dtype='<f4', count=ny).astype(float)
        zs = np.fromfile(f, dtype='<f4', count=nz).astype(float)
        blk = lambda: np.fromfile(f, dtype='<f4',
                                  count=nx * ny * nz).reshape(nz, ny, nx).astype(float)
        rho = blk(); u = blk(); v = blk(); _w = blk(); _th = blk(); qv = blk()

    dzs = np.gradient(zs)
    keep = zs <= ZTOP_FLUX

    def samp(a, xq, yq):
        """Bilinear sample of a (nz,ny,nx) block at scalar/array (xq,yq)."""
        fi = np.clip((np.asarray(xq) - xs[0]) / (xs[1] - xs[0]), 0, nx - 1.001)
        fj = np.clip((np.asarray(yq) - ys[0]) / (ys[1] - ys[0]), 0, ny - 1.001)
        i0 = fi.astype(int); j0 = fj.astype(int); ti = fi - i0; tj = fj - j0
        return (a[:, j0, i0] * (1 - ti) * (1 - tj) + a[:, j0, i0 + 1] * ti * (1 - tj)
                + a[:, j0 + 1, i0] * (1 - ti) * tj + a[:, j0 + 1, i0 + 1] * ti * tj)

    yc = PLOY + (np.arange(NY) + 0.5) * DX
    xc = PLOX + (np.arange(NX) + 0.5) * DX
    xlo_x, xhi_x = PLOX, PLOX + NX * DX
    ylo_y, yhi_y = PLOY, PLOY + NY * DX

    def face(xq, yq, comp):
        r = samp(rho, xq, yq); q = samp(qv, xq, yq)
        c = samp(u if comp == 'u' else v, xq, yq)
        return ((r * q * c) * (dzs * keep)[:, None]).sum() * DX

    return {'xlo': face(np.full(NY, xlo_x), yc, 'u'),
            'xhi': -face(np.full(NY, xhi_x), yc, 'u'),
            'ylo': face(xc, np.full(NX, ylo_y), 'v'),
            'yhi': -face(xc, np.full(NX, yhi_y), 'v')}


# ----------------------------------------------------------------- d02 -----
def d02_hourly():
    f = nc.Dataset('/app/ERF/wrfout_d02_2020-12-28_00_00_00')
    la = np.asarray(f.variables['XLAT'][0]); lo = np.asarray(f.variables['XLONG'][0])
    x, y = pyproj.Transformer.from_crs(4326, pyproj.CRS.from_proj4(PIN),
                                       always_xy=True).transform(lo, la)
    ii = np.floor((x - PLOX) / DX).astype(int)
    jj = np.floor((y - PLOY) / DX).astype(int)
    ok = (ii >= 0) & (ii < NX) & (jj >= 0) & (jj < NY)
    hgt = np.asarray(f.variables['HGT'][0])

    def to_grid(a2):
        s = np.zeros((NX, NY)); c = np.zeros((NX, NY))
        np.add.at(s, (ii[ok], jj[ok]), a2[ok]); np.add.at(c, (ii[ok], jj[ok]), 1.0)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)

    out = {}
    for h in range(1, 24):
        U = np.asarray(f.variables['U'][h]); V = np.asarray(f.variables['V'][h])
        u = 0.5 * (U[:, :, :-1] + U[:, :, 1:]); v = 0.5 * (V[:, :-1, :] + V[:, 1:, :])
        q = np.asarray(f.variables['QVAPOR'][h])
        p = np.asarray(f.variables['P'][h]) + np.asarray(f.variables['PB'][h])
        th = np.asarray(f.variables['T'][h]) + 300.0
        tk = th * (p / P0) ** (RD / CP)
        rho = p / (RD * tk * (1.0 + 0.608 * q))
        zst = (np.asarray(f.variables['PH'][h])
               + np.asarray(f.variables['PHB'][h])) / G
        dz = np.diff(zst, axis=0)
        zc = 0.5 * (zst[:-1] + zst[1:]) - hgt[None]
        low = zc < 1500.0
        wgt = np.where(low, rho * dz, 0.0)
        sw = np.maximum(wgt.sum(axis=0), 1e-9)
        ub = to_grid((u * wgt).sum(axis=0) / sw)
        vb = to_grid((v * wgt).sum(axis=0) / sw)
        pw = to_grid((rho * q * dz).sum(axis=0))
        out[h] = dict(ub=float(np.nanmean(ub[HI])), vb=float(np.nanmean(vb[HI])),
                      pw=float(np.nanmean(pw[HI])))
        del U, V, u, v, q, p, th, tk, rho, zst, dz, zc, wgt
    return out


def drc(ub, vb):
    return (np.degrees(np.arctan2(-ub, -vb)) + 360) % 360


def main():
    arms = {}
    runs = []
    for spec in sys.argv[1:]:
        lab, run = spec.split('=', 1)
        arms[lab] = erf_hourly(run)
        runs.append(run)
    print('[loaded ERF arms]', ', '.join(f'{k}:{len(v)}h' for k, v in arms.items()), flush=True)

    # The driver frames are the ones the FIRST arm was actually driven by, read
    # through its own run directory. This was hardcoded to run_c404_nsc, which
    # does not exist on a fresh pod -- and pointing it at some other arm's copy
    # would silently compare against the wrong frames.
    fdir = f'/app/ERF/{runs[0]}/CONUS404Data_3D'
    if not os.path.isdir(fdir):
        sys.exit(f'FATAL: no driving frames at {fdir}')
    fr = {}
    for p in sorted(glob.glob(f'{fdir}/ERF_IC_*.bin')):
        # ERF_IC_2020_12_28_00_00.bin -> [ERF, IC, yyyy, mm, dd, HH, MM.bin]
        parts = os.path.basename(p).split('_')
        if parts[4] == '28':
            fr[int(parts[5])] = frame_flux(p)
    print('[loaded frames]', sorted(fr), flush=True)

    d02 = d02_hourly()
    print('[loaded d02]', len(d02), 'hours\n', flush=True)

    lab0 = list(arms)[0]
    print('=== A. lateral water-vapour flux, 0-6000 m ASL, kg/s (inward positive) ===')
    print('  frames are the DRIVER; ERF should track it if the boundary is faithful')
    print(f'  {"h":>3} {"src":>14} {"xlo":>11} {"ylo":>11} {"INFLOW":>11}'
          f' {"xhi":>11} {"yhi":>11} {"NET":>11}')
    for h in range(1, 24):
        rows = [(k, arms[k][h]['F']) for k in arms if h in arms[k]]
        if h in fr:
            rows.append(('CONUS404', fr[h]))
        for nm, F in rows:
            inflow = F['xlo'] + F['ylo']
            net = sum(F.values())
            print(f'  {h:>3} {nm:>14} {F["xlo"]:11.3e} {F["ylo"]:11.3e} {inflow:11.3e}'
                  f' {F["xhi"]:11.3e} {F["yhi"]:11.3e} {net:11.3e}')

    print('\n=== A2. ERF inflow / driver inflow, at frame hours ===')
    print(f'  {"h":>3} ' + ' '.join(f'{k[:12]:>13}' for k in arms) + f' {"driver kg/s":>13}')
    for h in sorted(fr):
        if h < 1 or h > 23:
            continue
        dv = fr[h]['xlo'] + fr[h]['ylo']
        row = ' '.join(f'{(arms[k][h]["F"]["xlo"]+arms[k][h]["F"]["ylo"])/dv:13.3f}'
                       if h in arms[k] else f'{"--":>13}' for k in arms)
        print(f'  {h:>3} {row} {dv:13.3e}')

    print('\n=== B. sub-1500 m AGL wind over >500 m cells (speed m/s, dir deg) ===')
    print(f'  {"h":>3} ' + ' '.join(f'{k[:11]+" spd":>12}{k[:11]+" dir":>12}' for k in arms)
          + f' {"d02 spd":>12}{"d02 dir":>12}')
    for h in range(1, 24):
        cells = ''
        for k in arms:
            if h in arms[k]:
                a = arms[k][h]
                cells += f'{np.hypot(a["ub"],a["vb"]):12.2f}{drc(a["ub"],a["vb"]):12.1f}'
            else:
                cells += f'{"--":>12}{"--":>12}'
        d = d02[h]
        print(f'  {h:>3} {cells}{np.hypot(d["ub"],d["vb"]):12.2f}{drc(d["ub"],d["vb"]):12.1f}')

    print('\n=== C. column water vapour over >500 m cells (kg/m2) ===')
    print(f'  {"h":>3} ' + ' '.join(f'{k[:12]:>13}' for k in arms) + f' {"d02":>13}'
          + ''.join(f'{k[:8]+"/d02":>13}' for k in arms))
    for h in range(1, 24):
        vals = [arms[k][h]['pw'] if h in arms[k] else np.nan for k in arms]
        r = d02[h]['pw']
        print(f'  {h:>3} ' + ' '.join(f'{x:13.2f}' for x in vals) + f' {r:13.2f}'
              + ''.join(f'{x/r:13.2f}' for x in vals))


if __name__ == '__main__':
    main()
