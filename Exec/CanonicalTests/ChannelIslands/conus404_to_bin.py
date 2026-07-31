#!/usr/bin/env python3
"""CONUS404 (GDEX d559000 wrf3d/wrf2d) -> ERF hindcast .bin frames.

ERF needs NO source changes: this emits the existing custom binary frame format
(ERF_ReadCustomBinaryIC.H), so ERF_WeatherDataInterpolation reads CONUS404 frames
with the identical code path it uses for ERA5.

  3-D frame: int32 nx,ny,nz,ndata=8
             lat[nx*ny] lon[nx*ny] x[nx] y[ny] z[nz]
             then 8 blocks of nx*ny*nz: rho,u,v,w,theta,qv,qc,qr
  surface  : int32 nx,ny,1,ndata=6
             lat lon x y z, then sst,q_star,t_star,u_star,ls_mask,alb
             (q_star/t_star/u_star are read but never uploaded by ERF -- zeros)

What this REPLACES from the ERA5 path, and why each is no longer needed:
  * erf.hindcast_frame_from_T  -- WRF carries TK and P, so theta and rho are
    computed here, offline, where they can be checked. Leave the knob OFF.
  * the surface anchor          -- WRF's level 1 is ~10-30 m, so there is no
    150 m data floor to extrapolate beneath. Leave the anchor unset.
  * erftools level displacement and the fabricated duplicate bottom level --
    erftools is not in this path at all.
Kept and asserted: the PINNED LCC. It matters more here, because this is a
Lambert -> Lambert reprojection, and CONUS404's Lambert is on a SPHERE of radius
6,370,000 m (verified: that choice makes the source grid regular to 4000.00 +- 1 m;
WGS84 gives 4009.5 +- 2.8 m, a 0.24% scale error that would accumulate to ~2 km
across the subset).
"""
import os, sys, subprocess, datetime as dt
import numpy as np, pyproj
from scipy.ndimage import map_coordinates

BASE = "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d559000"
OURS = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 +lat_0=33.625000 "
        "+lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")
C404 = ("+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.100006103515625 "
        "+lon_0=-97.9000015258789 +a=6370000 +b=6370000 +units=m +no_defs")
NLEV_SRC = 50
RD, CP, P0 = 287.0, 1004.5, 1.0e5

# target frame grid: ERF prob_lo/hi with > 4 ERF cells of margin on every side.
#
# The ERF box comes from C404_PROB (x_lo y_lo x_hi y_hi, metres in the PINNED
# LCC). Defaults reproduce the channelislands-3km-192x96 grid this file shipped
# with, bit for bit -- the 2020-12-28 frames in pod_data were built from those
# literals, so the default path must not move. A second domain only needs
# C404_PROB set; nothing else here is domain-specific.
#
# MARGIN is generous on purpose: erftools insets the frame and ERF refuses a
# domain within 4 cells (12 km) of the frame edge
# (ERF_WeatherDataInterpolation.cpp:312-323). build_mapping() hard-raises if the
# target grid escapes the source subset, so an under-sized margin fails loudly.
if 'C404_PROB' in os.environ:
    _pl = [float(v) for v in os.environ['C404_PROB'].split()]
    if len(_pl) != 4:
        raise SystemExit('C404_PROB must be "x_lo y_lo x_hi y_hi" in metres')
    MARGIN = float(os.environ.get('C404_MARGIN', '66000.'))
    _snap = lambda v, up: (np.ceil(v / 6000.) if up else np.floor(v / 6000.)) * 6000.
    XS = np.arange(_snap(_pl[0] - MARGIN, False), _snap(_pl[2] + MARGIN, True) + 1., 6000.)
    YS = np.arange(_snap(_pl[1] - MARGIN, False), _snap(_pl[3] + MARGIN, True) + 1., 6000.)
else:
    # The 192x96 literals this file shipped with. The pod_data frames were built
    # from exactly these, so the unset path stays byte-for-byte identical rather
    # than being re-derived from a margin rule that would round differently.
    XS = np.arange(-450000., 252001., 6000.)
    YS = np.arange(-230000., 190001., 6000.)
# stretched heights ASL; must exceed the ERF domain top (19003 m) with margin
_s = np.arange(46) / 45.0
ZS = (25000.0 * (np.exp(3.6 * _s) - 1.0) / (np.exp(3.6) - 1.0)).astype(np.float64)
ZS[0] = 0.0


def dods(url_var, shape, retries=3):
    """Fetch one variable via OPeNDAP .dods and return it as float64.

    THREDDS/Tomcat rejects raw [ ] : in the request target, and curl glob-expands
    brackets, so both must be handled: percent-encode and pass -g.
    """
    enc = url_var.replace('[', '%5B').replace(']', '%5D')
    for a in range(retries):
        r = subprocess.run(['curl', '-g', '-sS', '--max-time', '900', enc],
                           capture_output=True)
        raw = r.stdout
        i = raw.find(b'\nData:\n')
        if i > 0:
            n = int(np.prod(shape))
            buf = raw[i + 7 + 8: i + 7 + 8 + 4 * n]
            if len(buf) == 4 * n:
                return np.frombuffer(buf, dtype='>f4').reshape(shape).astype(np.float64)
        print(f'    retry {a+1} ({len(raw)} bytes)', flush=True)
    raise RuntimeError(f'OPeNDAP fetch failed: {url_var[:120]}')


def build_mapping(j0, j1, i0, i1):
    """Fractional source indices for every target point, plus target lat/lon."""
    X, Y = np.meshgrid(XS, YS, indexing='ij')
    tlon, tlat = pyproj.Transformer.from_crs(
        pyproj.CRS.from_proj4(OURS), 4326, always_xy=True).transform(X, Y)
    tx, ty = pyproj.Transformer.from_crs(
        4326, pyproj.CRS.from_proj4(C404), always_xy=True).transform(tlon, tlat)
    la = np.load(os.environ['C404_LAT']); lo = np.load(os.environ['C404_LON'])
    sx, sy = pyproj.Transformer.from_crs(
        4326, pyproj.CRS.from_proj4(C404), always_xy=True).transform(
        lo[j0, i0], la[j0, i0])
    fi = (tx - sx) / 4000.0
    fj = (ty - sy) / 4000.0
    ny_s, nx_s = j1 - j0 + 1, i1 - i0 + 1
    if fi.min() < -0.5 or fi.max() > nx_s - 0.5 or fj.min() < -0.5 or fj.max() > ny_s - 0.5:
        raise RuntimeError(f'target grid escapes the source subset: '
                           f'i {fi.min():.1f}..{fi.max():.1f} of {nx_s}, '
                           f'j {fj.min():.1f}..{fj.max():.1f} of {ny_s}')
    return fi, fj, tlat, tlon


def hinterp(a3, fj, fi):
    """Bilinear horizontal interpolation of a (nlev, ny, nx) source block."""
    out = np.empty((a3.shape[0], fi.shape[0], fi.shape[1]))
    for k in range(a3.shape[0]):
        out[k] = map_coordinates(a3[k], [fj, fi], order=1, mode='nearest')
    return out


def write_frame(path, lat, lon, x, y, z, fields):
    nx, ny, nz = len(x), len(y), len(z)
    with open(path, 'wb') as f:
        np.array([nx, ny, nz, len(fields)], dtype='<i4').tofile(f)
        # lat/lon are stored row-major as [j*nx + i]
        np.asarray(lat.T, dtype='<f4').ravel().tofile(f)
        np.asarray(lon.T, dtype='<f4').ravel().tofile(f)
        np.asarray(x, dtype='<f4').tofile(f)
        np.asarray(y, dtype='<f4').tofile(f)
        np.asarray(z, dtype='<f4').tofile(f)
        for a in fields:                       # (nx,ny,nz) -> k,j,i order
            np.asarray(np.transpose(a, (2, 1, 0)), dtype='<f4').ravel().tofile(f)


def main():
    # Fail on a missing input NOW, not with a bare KeyError partway through the
    # surface step of frame 1 -- by which point several minutes of OPeNDAP
    # fetching have already been spent.
    missing = [k for k in ('C404_BBOX', 'C404_LAT', 'C404_LON', 'C404_LANDMASK')
               if k not in os.environ]
    if missing:
        sys.exit(f'FATAL: unset environment variable(s): {", ".join(missing)}\n'
                 f'  C404_LAT / C404_LON / C404_LANDMASK are the full 1015x1367\n'
                 f'  CONUS404 XLAT, XLONG and LANDMASK saved as .npy. LANDMASK is in\n'
                 f'  INVARIANT/USGS404_geo_em_d01.nc (it is NOT in wrf2d/wrf3d).\n'
                 f'  C404_BBOX comes from pod/make_c404_bbox.py.')

    outdir3 = sys.argv[1]; outdirS = sys.argv[2]
    os.makedirs(outdir3, exist_ok=True); os.makedirs(outdirS, exist_ok=True)
    j0, j1, i0, i1 = [int(v) for v in np.load(os.environ['C404_BBOX'])]
    fi, fj, tlat, tlon = build_mapping(j0, j1, i0, i1)
    print(f'target {len(XS)}x{len(YS)}x{len(ZS)}  z {ZS[0]:.0f}..{ZS[-1]:.0f} m '
          f'(ERF top 19003) ; source subset y{j0}..{j1} x{i0}..{i1}', flush=True)

    NH = int(os.environ.get("C404_FRAME_HOURS", "3"))
    NF = int(os.environ.get("C404_NFRAMES", "9"))
    times = [dt.datetime(2020, 12, 28) + dt.timedelta(hours=NH * k) for k in range(NF)]
    S3 = f'[{j0}:1:{j1}][{i0}:1:{i1}]'
    for t in times:
        wy = 'wy2021'; mo = t.strftime('%Y%m')
        stamp = t.strftime('%Y-%m-%d_%H:00:00').replace(':', '%3A')
        u3 = f'{BASE}/{wy}/{mo}/wrf3d_d01_{stamp}.nc.dods?'
        g = lambda v, extra='': dods(u3 + f'{v}[0:1:0][0:1:{NLEV_SRC-1}]' + (extra or S3),
                                     (1, NLEV_SRC, j1-j0+1, i1-i0+1))[0]
        TK = g('TK'); P = g('P'); QV = g('QVAPOR'); QC = g('QCLOUD'); QR = g('QRAIN')
        U = dods(u3 + f'U[0:1:0][0:1:{NLEV_SRC-1}][{j0}:1:{j1}][{i0}:1:{i1+1}]',
                 (1, NLEV_SRC, j1-j0+1, i1-i0+2))[0]
        V = dods(u3 + f'V[0:1:0][0:1:{NLEV_SRC-1}][{j0}:1:{j1+1}][{i0}:1:{i1}]',
                 (1, NLEV_SRC, j1-j0+2, i1-i0+1))[0]
        W = dods(u3 + f'W[0:1:0][0:1:{NLEV_SRC}][{j0}:1:{j1}][{i0}:1:{i1}]',
                 (1, NLEV_SRC+1, j1-j0+1, i1-i0+1))[0]
        Z = dods(u3 + f'Z[0:1:0][0:1:{NLEV_SRC}][{j0}:1:{j1}][{i0}:1:{i1}]',
                 (1, NLEV_SRC+1, j1-j0+1, i1-i0+1))[0]
        U = 0.5 * (U[:, :, :-1] + U[:, :, 1:])
        V = 0.5 * (V[:, :-1, :] + V[:, 1:, :])
        W = 0.5 * (W[:-1] + W[1:]); ZC = 0.5 * (Z[:-1] + Z[1:])
        TH = TK * (P0 / P) ** (RD / CP)
        RHO = P / (RD * TK * (1.0 + 0.608 * QV))

        H = {k: hinterp(v, fj, fi) for k, v in
             (('rho', RHO), ('u', U), ('v', V), ('w', W),
              ('th', TH), ('qv', QV), ('qc', QC), ('qr', QR), ('z', ZC))}
        nx, ny = len(XS), len(YS)
        out = {k: np.empty((nx, ny, len(ZS))) for k in H if k != 'z'}
        for a in range(nx):
            for b in range(ny):
                zc = H['z'][:, a, b]
                for k in out:
                    out[k][a, b] = np.interp(ZS, zc, H[k][:, a, b])
        fn = os.path.join(outdir3, 'ERF_IC_' + t.strftime('%Y_%m_%d_%H_%M') + '.bin')
        write_frame(fn, tlat, tlon, XS, YS, ZS,
                    [out['rho'], out['u'], out['v'], out['w'],
                     out['th'], out['qv'], out['qc'], out['qr']])
        print(f'  wrote {os.path.basename(fn)}  rho {out["rho"].mean():.4f} '
              f'th {out["th"].mean():.1f} qv {out["qv"].mean():.5f}', flush=True)

        # ---- surface frame ----
        u2 = f'{BASE}/{wy}/{mo}/wrf2d_d01_{stamp}.nc.dods?'
        s2 = lambda v: dods(u2 + f'{v}[0:1:0]{S3}', (1, j1-j0+1, i1-i0+1))[0]
        SST = s2('SST'); TSK = s2('TSK'); ALB = s2('ALBEDO')
        LM = np.load(os.environ['C404_LANDMASK'])[j0:j1+1, i0:i1+1]
        # masked SST, reworked onto TSK/LANDMASK: use SST over water where it is
        # physical, else TSK; land cells carry TSK and are ignored by the ocean path.
        water = LM < 0.5
        bad = ~((SST > 271.0) & (SST < 305.0))
        sst = np.where(water & ~bad, SST, TSK)
        SH = {k: hinterp(v[None], fj, fi)[0] for k, v in
              (('sst', sst), ('lm', LM.astype(float)), ('alb', ALB))}
        z1 = np.array([0.0]); zero = np.zeros((nx, ny, 1))
        write_frame(os.path.join(outdirS, 'ERF_IC_' + t.strftime('%Y_%m_%d_%H_%M') + '.bin'),
                    tlat, tlon, XS, YS, z1,
                    [SH['sst'][..., None], zero, zero, zero,
                     SH['lm'][..., None], SH['alb'][..., None]])
    print('done')


if __name__ == '__main__':
    main()
