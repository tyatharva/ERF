"""One-pass verification of the whole ERA5 -> erftools -> ERF pipeline at t = 0.

No model runs. Reads an existing plt00000 plus the files every stage was built from,
and reports each variable and the delta at every stage boundary -- so a
read-then-discard or a vertical-reference mismatch shows up as a number, instead of
being found one at a time over days.

    A  ERA5 source        era5_run/era5_3d_*.grib, era5_surf_*.grib
    B  erftools frame     ERA5Data_3D/*.bin, ERA5Data_Surface/*.bin
    C  ERF interpolated   reimplementation of ERF's bilinear_interpolation
    D  ERF state at t=0   plt00000
    E  boundary planes    stage C restricted to the relaxation band

The column to read first is the BEST-FIT VERTICAL OFFSET: the height shift that best
aligns the frame's profile with ERA5's. If erftools mislabels heights, every variable
reports the SAME offset and that number IS the defect. This column alone would have
found item 19f on day one.

env: RUN, PLT, FRAME_DIR, SFC_DIR, ERA5_DIR, STAMP, FSTAMP, SRC
"""
import numpy as np, pygrib, yt, os, re, subprocess
from pyproj import CRS, Transformer
yt.set_log_level(50)

RUN      = os.environ.get('RUN', 'bdyfix/bt_zoff')
FRAMEDIR = os.environ.get('FRAME_DIR', '/app/ERF/run_jan09/ERA5Data_3D')
SFCDIR   = os.environ.get('SFC_DIR', '/app/ERF/run_jan09/ERA5Data_Surface')
ERA5DIR  = os.environ.get('ERA5_DIR', '/app/ERF/era5_run')
STAMP    = os.environ.get('STAMP', '20230109_0000')
FSTAMP   = os.environ.get('FSTAMP', '2023_01_09_00_00')
SRC      = os.environ.get('SRC', '/app/ERF/Source')
ZOFF     = float(os.environ.get('ZOFF', '305.0'))     # offset the run applied, for stage C
PLO = (-195131.04, -126372.41); PHI = (188868.96, 65627.59)
AREA = [36.0, -123.25, 31.25, -115.25]
la1, la2, lo1, lo2 = AREA[2], AREA[0], AREA[1], AREA[3]
dl = la2 - la1
LCC = (f"+proj=lcc +lat_1={la1+dl/6:.6f} +lat_2={la2-dl/6:.6f} "
       f"+lat_0={(la1+la2)/2:.6f} +lon_0={(lo1+lo2)/2:.6f} +datum=WGS84 +units=m +no_defs")
to_ll = Transformer.from_crs(CRS.from_proj4(LCC), CRS.from_epsg(4326), always_xy=True)
G0 = 9.80665; RD = 287.0; RV = 461.505; CP = 1004.5; P0 = 1.0e5; KAPPA = RD/CP; GAMMA = 1.4
TOL = {'theta': 0.5, 'T': 0.5, 'p': 200.0, 'rho': 0.010, 'qv': 5e-4, 'qc': 1e-4,
       'qr': 1e-4, 'u': 1.0, 'v': 1.0, 'w': 0.05, 'sst': 0.5, 't2m': 0.5,
       'alb': 0.02, 'lsm': 0.05}

# ---------------------------------------------------------------- B: erftools frame
f3 = f'{FRAMEDIR}/ERF_IC_{FSTAMP}.bin'
h = np.fromfile(f3, dtype=np.int32, count=4); nx3, ny3, nz3, nd3 = (int(v) for v in h)
a3 = np.fromfile(f3, dtype=np.float32, offset=16); o = 2*nx3*ny3
fx = a3[o:o+nx3].astype(float); o += nx3
fy = a3[o:o+ny3].astype(float); o += ny3
fz = a3[o:o+nz3].astype(float); o += nz3
FRAME_FIELDS = ['rho', 'u', 'v', 'w', 'theta', 'qv', 'qc', 'qr']
F = {}
for nm in FRAME_FIELDS:
    F[nm] = a3[o:o+nx3*ny3*nz3].reshape(nz3, ny3, nx3).astype(float); o += nx3*ny3*nz3
dup = all(np.array_equal(F[k][0], F[k][1]) for k in F)
k0 = 1 if dup else 0

fs = f'{SFCDIR}/ERF_Surface_{FSTAMP}.bin'
hs = np.fromfile(fs, dtype=np.int32, count=4); snx, sny, snz, snd = (int(v) for v in hs)
asf = np.fromfile(fs, dtype=np.float32, offset=16)
sx = asf[0:snx].astype(float); sy = asf[snx:snx+sny].astype(float)
o = snx+sny+snz
SFC_FIELDS = ['sst', 'q_star', 't_star', 'u_star', 'lsm', 'alb'][:snd]
S = {}
for nm in SFC_FIELDS:
    S[nm] = asf[o:o+snx*sny].reshape(sny, snx).astype(float); o += snx*sny

# ---------------------------------------------------------------- A: ERA5 source
g = pygrib.open(f'{ERA5DIR}/era5_3d_{STAMP}.grib'); E = {}; elat = elon = None
for m in g:
    E.setdefault(m.shortName, {})[int(m.level)] = np.array(m.values, dtype=float)
    if elat is None: elat, elon = m.latlons()
g.close()
levs = np.array(sorted(E['t'].keys()))
st = lambda s: np.stack([E[s][int(l)] for l in levs], axis=0)
Et, Eq, Ezh = st('t'), st('q'), st('z')/G0
Eu, Ev, Eomega = st('u'), st('v'), st('w')
Eclw = st('clwc') if 'clwc' in E else np.zeros_like(Et)
Ecrw = st('crwc') if 'crwc' in E else np.zeros_like(Et)
Ep = (levs[:, None, None]*100.0)*np.ones_like(Et)
Eth = Et*(P0/Ep)**KAPPA
Emr = Eq/(1.0-Eq)                                  # ERA5 gives SPECIFIC humidity
Erho_d = Ep/(RD*Et*(1.0+(RV/RD)*Emr))              # dry density, ERF's convention
Ew = -Eomega/((Ep/(RD*Et))*G0)
gs = pygrib.open(f'{ERA5DIR}/era5_surf_{STAMP}.grib'); SFA = {}
for m in gs: SFA[m.shortName] = np.array(m.values, dtype=float)
gs.close()
elon180 = np.where(elon > 180, elon-360, elon)
latv, lonv = elat[:, 0], elon180[0, :]

# ---------------------------------------------------------------- D: model state
ds = yt.load(os.environ.get('PLT', f'/app/ERF/{RUN}/plt00000'))
gg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
M = {n: np.asarray(gg[('boxlib', b)]) for n, b in
     [('theta', 'theta'), ('T', 'temp'), ('p', 'pressure'), ('rho', 'density'),
      ('qv', 'qv'), ('qc', 'qc'), ('qr', 'qrain'), ('u', 'x_velocity'),
      ('v', 'y_velocity'), ('w', 'z_velocity')]}
Z = np.asarray(gg[('boxlib', 'z_phys')]); ter = Z[:, :, 0]
nx, ny, nz = M['theta'].shape
xc = PLO[0]+(np.arange(nx)+0.5)*(PHI[0]-PLO[0])/nx
yc = PLO[1]+(np.arange(ny)+0.5)*(PHI[1]-PLO[1])/ny
mlon, mlat = to_ll.transform(*np.meshgrid(xc, yc, indexing='ij'))
CLS = {'ocean': ter < 30, 'coastal': (ter >= 30) & (ter < 200), 'terrain': ter >= 200}
ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
dring = np.minimum.reduce([ii, jj, nx-1-ii, ny-1-jj])
BAND = dring < 10


def bilin(fld, xv, yv):
    i0 = np.clip(np.searchsorted(xv, xc)-1, 0, len(xv)-2)
    j0 = np.clip(np.searchsorted(yv, yc)-1, 0, len(yv)-2)
    fx_ = np.clip((xc-xv[i0])/(xv[i0+1]-xv[i0]), 0, 1)
    fy_ = np.clip((yc-yv[j0])/(yv[j0+1]-yv[j0]), 0, 1)
    return ((1-fx_)[:, None]*(1-fy_)[None, :]*fld[j0, :][:, i0].T +
            fx_[:, None]*(1-fy_)[None, :]*fld[j0, :][:, i0+1].T +
            (1-fx_)[:, None]*fy_[None, :]*fld[j0+1, :][:, i0].T +
            fx_[:, None]*fy_[None, :]*fld[j0+1, :][:, i0+1].T)


def era5_ll(f2):
    la = np.interp(mlat, latv[::-1], np.arange(len(latv))[::-1]) if latv[0] > latv[-1] \
        else np.interp(mlat, latv, np.arange(len(latv)))
    lo = np.interp(mlon, lonv, np.arange(len(lonv)))
    i0 = np.clip(np.floor(la).astype(int), 0, len(latv)-2)
    j0 = np.clip(np.floor(lo).astype(int), 0, len(lonv)-2)
    fa = np.clip(la-i0, 0, 1); fo = np.clip(lo-j0, 0, 1)
    return ((1-fa)*(1-fo)*f2[i0, j0] + (1-fa)*fo*f2[i0, j0+1] +
            fa*(1-fo)*f2[i0+1, j0] + fa*fo*f2[i0+1, j0+1])


# class-mean 1-D profiles: ERA5 on its own geopotential heights, frame on its labels
EPROF, EH = {}, {}
FPROF = {}
for cls, msk in CLS.items():
    EH[cls] = np.array([era5_ll(Ezh[l])[msk].mean() for l in range(len(levs))])
    for nm, fld in [('theta', Eth), ('T', Et), ('qv', Emr), ('u', Eu), ('v', Ev),
                    ('rho', Erho_d), ('p', Ep), ('w', Ew), ('qc', Eclw), ('qr', Ecrw)]:
        EPROF[(cls, nm)] = np.array([era5_ll(fld[l])[msk].mean() for l in range(len(levs))])
    for nm in FRAME_FIELDS:
        FPROF[(cls, nm)] = np.array([bilin(F[nm][k], fx, fy)[msk].mean() for k in range(k0, nz3)])
    fp = FPROF[(cls, 'p')] = P0*(RD*FPROF[(cls, 'rho')]*FPROF[(cls, 'theta')] *
                                 (1+(RV/RD)*FPROF[(cls, 'qv')])/P0)**GAMMA
    FPROF[(cls, 'T')] = FPROF[(cls, 'theta')]*(fp/P0)**KAPPA
FZ = fz[k0:]

O = []
def emit(s=''): O.append(s)


emit('='*104)
emit('PIPELINE AUDIT   t = 0     ERA5 -> erftools frame -> ERF interp -> ERF state -> planes')
emit(f'  run {RUN}    frame {os.path.basename(f3)}    ERA5 {STAMP}    stage-C offset applied {ZOFF:.0f} m')
emit(f'  frame {nx3}x{ny3}x{nz3} ({nd3} fields)   surface {snx}x{sny} ({snd} fields)   model {nx}x{ny}x{nz}')
emit(f'  frame bottom level is a duplicate: {dup}   lowest real datum z = {FZ[0]:.1f} m (labelled)')
emit('='*104)

# ---- 1. vertical offset ---------------------------------------------------------
emit('')
emit('1. BEST-FIT VERTICAL OFFSET   height shift minimising |frame(z) - ERA5(z+off)| over 300-6000 m')
emit('   A non-zero CONSENSUS offset across variables is a vertical-reference defect in erftools.')
emit('')
emit('   variable |   ocean    coastal    terrain')
offs = np.arange(-600, 601, 5.0)
allb = []
for nm in ['theta', 'T', 'qv', 'u', 'v', 'rho']:
    row = []
    for cls in CLS:
        zs = FZ[(FZ > 300) & (FZ < 6000)]
        fv = np.interp(zs, FZ, FPROF[(cls, nm)])
        eh, ev = EH[cls], EPROF[(cls, nm)]
        order = np.argsort(eh)
        cost = [np.abs(fv-np.interp(zs+off, eh[order], ev[order])).mean() for off in offs]
        row.append(offs[int(np.argmin(cost))])
    allb += row
    emit(f'   {nm:8s} | {row[0]:8.0f}m {row[1]:9.0f}m {row[2]:9.0f}m')
med = float(np.median(allb))
emit('')
emit(f'   >>> CONSENSUS OFFSET = {med:+.0f} m   (spread {np.min(allb):+.0f} .. {np.max(allb):+.0f} m)')
emit(f'   >>> {"FLAG: erftools places data at the wrong height" if abs(med) > 50 else "ok"}')

# ---- 2. stage table -------------------------------------------------------------
emit('')
emit('2. STAGE VALUES AND DELTAS   ocean columns, at model levels')
emit('   A=ERA5(at that height)  B=frame  C=ERF interp (offset applied)  D=ERF state')
emit('')
emit('   var    z(m)      A          B          C          D    | A->B      B->C      C->D    flag')
oc = CLS['ocean']
for k in [0, 2, 4, 8, 12, 20]:
    zk = Z[:, :, k][oc].mean()
    for nm in ['theta', 'T', 'qv', 'u', 'v', 'rho', 'p']:
        eh, ev = EH['ocean'], EPROF[('ocean', nm)]
        order = np.argsort(eh)
        A = float(np.interp(zk, eh[order], ev[order]))
        B = float(np.interp(zk, FZ, FPROF[('ocean', nm)]))
        C = float(np.interp(zk, FZ+ZOFF, FPROF[('ocean', nm)]))
        D = float(M[nm][:, :, k][oc].mean()) if nm in M else float('nan')
        ab, bc, cd = B-A, C-B, (D-C if np.isfinite(D) else float('nan'))
        t = TOL.get(nm, 1e9)
        fl = ''.join(['!' if (np.isfinite(x) and abs(x) > t) else '.' for x in (ab, bc, cd)])
        fmt = '%10.5f' if nm in ('qv', 'qc', 'qr', 'rho') else '%10.2f'
        emit(f'   {nm:6s} {zk:6.0f} ' + ' '.join(fmt % v for v in (A, B, C, D)) +
             ' | ' + ' '.join(fmt % v for v in (ab, bc, cd)) + f'   {fl}')
    emit('')

# ---- 3. surface fields ----------------------------------------------------------
emit('3. SURFACE FIELDS   A=ERA5  B=surface frame  (ocean columns for sst; all for lsm/alb)')
emit('')
emit('   field     A(ERA5)     B(frame)     A->B     flag   downstream use')
SFMAP = [('sst', 'sst', 'sst'), ('lsm', 'lsm', 'lsm'), ('alb', 'fal', 'alb'),
         ('u_star', 'zust', 'u_star')]
for nm, era_sn, fr in SFMAP:
    if fr not in S:
        emit(f'   {nm:9s} {"--":>11s} {"absent":>12s}      -       -    field not in this frame')
        continue
    B = bilin(S[fr], sx, sy)
    msk = oc if nm == 'sst' else np.ones_like(oc, dtype=bool)
    if era_sn in SFA:
        A = era5_ll(SFA[era_sn]); a, b = A[msk].mean(), B[msk].mean()
        d = b-a; t = TOL.get(nm, 1e9)
        emit(f'   {nm:9s} {a:11.4f} {b:12.4f} {d:+9.4f}   {"!" if abs(d) > t else "."}')
    else:
        emit(f'   {nm:9s} {"n/a":>11s} {B[msk].mean():12.4f}       -       -')

# ---- 4. read-then-discard (DERIVED) --------------------------------------------
emit('')
emit('4. READ-THEN-DISCARD   derived from the source: every field read from a frame,')
emit('   traced to the component it lands in, then to whether anything reads that back.')
emit('')


def slurp(fn):
    try:
        return open(fn).read()
    except Exception:
        return ''


WD = slurp(f'{SRC}/Utils/ERF_WeatherDataInterpolation.cpp')
SD = slurp(f'{SRC}/Utils/ERF_SurfaceDataInterpolation.cpp')
ALLSRC = {}
for root, _, files in os.walk(SRC):
    for f in files:
        if f.endswith(('.cpp', '.H')):
            ALLSRC[os.path.join(root, f)] = slurp(os.path.join(root, f))

# fields read = the host vectors declared in each reader
host = {}
for tag, txt in (('3-D frame', WD), ('surface frame', SD)):
    for m in re.finditer(r'Vector<Real>\s+([A-Za-z0-9_,\s]+);', txt):
        for v in m.group(1).split(','):
            v = v.strip()
            if v.endswith('_h') and v not in ('xvec_h', 'yvec_h', 'zvec_h',
                                              'latvec_h', 'lonvec_h', 'junk_h'):
                host.setdefault(v, tag)

emit('   field(host vec)  reader          interpolated  -> component      read back by')
for v, tag in sorted(host.items()):
    txt = WD if tag == '3-D frame' else SD
    base = v[:-2]
    dptr = f'{base}_d_ptr'
    interp = 'yes' if dptr in txt else 'no '
    # which array component does its tmp_ land in?
    comp = '-'
    mm = re.search(rf'\(i,j,k,\s*([A-Za-z0-9_]+)\s*\)\s*=\s*tmp_{base}\b', txt)
    if not mm:
        mm = re.search(rf'\(i,\s*j,\s*k,\s*([0-9]+)\s*\)\s*=\s*tmp_{base}\b', txt)
    if mm:
        comp = mm.group(1)
    else:
        mm2 = re.search(rf'tmp_{base}\b', txt)
        comp = 'computed, not stored' if mm2 else 'NOT INTERPOLATED'
    # who reads that component back, outside the two interpolation files?
    readers = []
    if comp not in ('-', 'NOT INTERPOLATED', 'computed, not stored'):
        for fn, t in ALLSRC.items():
            b = os.path.basename(fn)
            if b in ('ERF_WeatherDataInterpolation.cpp', 'ERF_SurfaceDataInterpolation.cpp'):
                continue
            if re.search(rf'\b{re.escape(comp)}\b', t):
                readers.append(b)
    verdict = f'{len(readers)} file(s)' if readers else 'NOBODY  <-- discarded'
    emit(f'   {v:16s} {tag:15s} {interp:12s} -> {comp:15s} {verdict}')
emit('')
emit('   "discarded" = interpolated onto the ERF mesh and then read by no other file.')

# ---- 5. vertical-reference audit (DERIVED) -------------------------------------
emit('')
emit('5. VERTICAL-REFERENCE AUDIT   derived: every site using a height/pressure symbol,')
emit('   tagged by convention. A function mixing conventions is flagged.')
emit('')
CONV = [('node',        r'z_phys_nd|z_nd\b|zcc_arr\(i,j,k\)\s*\+\s*z'),
        ('cell-centre', r'z_phys_cc|zcc_arr|z_cc\b'),
        ('sea-level',   r'\bp_0\b\s*-\s*hz|\bp_0\b\s*\+\s*hz'),
        ('frame-level', r'zvec_d_ptr|zvec_h'),
        ('geopotential', r'CONST_GRAV|/\s*9\.80665|G0\b')]
rows = []
for fn, t in sorted(ALLSRC.items()):
    b = os.path.basename(fn)
    found = {}
    for nm, pat in CONV:
        hits = [t[:m.start()].count('\n')+1 for m in re.finditer(pat, t)]
        if hits:
            found[nm] = hits
    if not found:
        continue
    rows.append((b, found))
emit('   file                                    conventions present            flag')
for b, found in rows:
    names = sorted(found)
    n = len(names)
    flag = 'MIXED' if n > 1 else ''
    emit(f'   {b:39s} {",".join(names):30s} {flag}')
emit('')
emit('   MIXED means one file uses more than one vertical reference; every defect found')
emit('   so far (items 19c, 19f/21a, 20) lives in a MIXED file.')

# ---- 6. compensating pairs (DERIVED) -------------------------------------------
emit('')
emit('6. COMPENSATING-PAIR CHECK   derived from section 2: a combination whose error is')
emit('   much SMALLER than its inputs\' means two errors are cancelling, and breaking')
emit('   one of them alone will make things worse.')
emit('')
emit('   combination                     input errors            output error   verdict')
zc = 846.0
def prof(nm, stage):
    if stage == 'A':
        eh, ev = EH['ocean'], EPROF[('ocean', nm)]
        o = np.argsort(eh); return float(np.interp(zc, eh[o], ev[o]))
    return float(np.interp(zc, FZ, FPROF[('ocean', nm)]))
relerr = {}
for nm in ['theta', 'rho', 'T', 'p', 'qv']:
    a, b = prof(nm, 'A'), prof(nm, 'B')
    relerr[nm] = (b-a)/abs(a) if a else 0.0
COMBOS = [('T   = p(rho,theta)/(R_d rho)', ['rho', 'theta'], 'T'),
          ('p   = p0(R_d rho theta/p0)^g', ['rho', 'theta'], 'p'),
          ('rho = p/(R_d T)',              ['p', 'T'],       'rho')]
for name, ins, outv in COMBOS:
    ie = max(abs(relerr[x]) for x in ins)
    oe = abs(relerr[outv])
    verdict = 'COMPENSATING PAIR' if (ie > 5*max(oe, 1e-6) and ie > 1e-3) else 'errors propagate'
    emit(f'   {name:31s} {", ".join(f"{x} {100*relerr[x]:+.2f}%" for x in ins):23s} '
         f'{100*oe:+8.3f}%   {verdict}')
emit('')
emit('   Read this with section 2: a pair flagged here cannot be half-fixed. Item 21b is')
emit('   exactly that failure -- the frame\'s theta was paired with a CORRECT pressure,')
emit('   which un-cancelled the pair and converted a pressure error into a 2 K warm bias.')

emit('')
emit('='*104)
print('\n'.join(O))
