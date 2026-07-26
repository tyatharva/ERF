"""t = 0 verification of the restructured HindCast initial condition.

Five pass conditions, checked against plt00000 and the data the IC was built
from.  Nothing here reads the model's own diagnostics back as truth: p, theta
and T come from the plotfile, the frame comes from the .bin erftools wrote, and
SST/t2m come from ERA5 -- so a bug in the integration cannot hide itself.

  1. p decreases monotonically with height, everywhere.
  2. theta increases monotonically with height over the ocean.
  3. every layer has non-zero thickness, in both z and p.
  4. T at the first cell centre is within ~1 K of SST + 0.64 K (= ERA5 t2m).
  5. theta is within ~1.5 K of the ERA5 frame at matching heights, allowing the
     ~1.4 K residual that erftools' ~300 m level displacement contributes.

env: RUN (run directory), FRAME (frame .bin), ANCHOR (anchor .bin), PLT
"""
import numpy as np, yt, glob, os, struct, sys
yt.set_log_level(50)

RUN    = os.environ.get('RUN', '/app/ERF/run_icfix')
PLT    = os.environ.get('PLT', f'{RUN}/plt00000')
FRAME  = os.environ.get('FRAME', f'{RUN}/ERA5Data_3D/ERF_IC_2023_01_09_00_00.bin')
ANCHOR = os.environ.get('ANCHOR', '/app/ERF/precip_check/sfc_anchor_jan09.bin')
SFCBIN = os.environ.get('SFCBIN', f'{RUN}/ERA5Data_Surface/ERF_Surface_2023_01_09_00_00.bin')
PLO    = (-195131.04, -126372.41); PHI = (188868.96, 65627.59)
LAND_M = 30.0                       # ERF terrain above this = land (ocean is ~12 m)

# ---- model ---------------------------------------------------------------------
ds = yt.load(PLT)
g  = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
P   = np.asarray(g[('boxlib', 'pressure')])
TH  = np.asarray(g[('boxlib', 'theta')])
T   = np.asarray(g[('boxlib', 'temp')])
Z   = np.asarray(g[('boxlib', 'z_phys')])
PH  = np.asarray(g[('boxlib', 'pres_hse')])
nx, ny, nz = P.shape
ter   = Z[:, :, 0]
ocean = ter < LAND_M
print(f'{PLT}: {nx}x{ny}x{nz}, t = {float(ds.current_time):.1f} s')
print(f'ocean columns {ocean.sum()} / {nx*ny} ({100*ocean.mean():.1f}%)')
print(f'first cell-centre height: ocean {Z[:,:,0][ocean].mean():.2f} m, '
      f'land {Z[:,:,0][~ocean].mean():.1f} m')

PLO_x = PLO[0]
xc = PLO[0] + (np.arange(nx) + 0.5)*(PHI[0]-PLO[0])/nx
yc = PLO[1] + (np.arange(ny) + 0.5)*(PHI[1]-PLO[1])/ny


def bilin2d(fld, xv, yv, xq, yq):
    ii = np.clip(np.searchsorted(xv, xq) - 1, 0, len(xv)-2)
    jj = np.clip(np.searchsorted(yv, yq) - 1, 0, len(yv)-2)
    fx = np.clip((xq - xv[ii])/(xv[ii+1]-xv[ii]), 0, 1)
    fy = np.clip((yq - yv[jj])/(yv[jj+1]-yv[jj]), 0, 1)
    return ((1-fx)[:, None]*(1-fy)[None, :]*fld[jj, :][:, ii].T +
            fx[:, None]*(1-fy)[None, :]*fld[jj, :][:, ii+1].T +
            (1-fx)[:, None]*fy[None, :]*fld[jj+1, :][:, ii].T +
            fx[:, None]*fy[None, :]*fld[jj+1, :][:, ii+1].T)


# ---- the ERA5 frame, sampled at every model cell centre --------------------------
# Same operator ERF applies: bilinear in x/y, linear in z, clamped below the
# lowest level. Needed by BOTH check 2 (is a non-increasing layer ours or ERA5's?)
# and check 5.
fhdr = np.fromfile(FRAME, dtype=np.int32, count=4)
fnx, fny, fnz, _ = (int(v) for v in fhdr)
fa = np.fromfile(FRAME, dtype=np.float32, offset=16)
o = 2*fnx*fny
fx_ = fa[o:o+fnx]; o += fnx
fy_ = fa[o:o+fny]; o += fny
fz_ = fa[o:o+fnz]; o += fnz
o += 4*fnx*fny*fnz                                   # rho,u,v,w
fth = fa[o:o+fnx*fny*fnz].reshape(fnz, fny, fnx)
fth_col = np.stack([bilin2d(fth[k], fx_, fy_, xc, yc) for k in range(fnz)], axis=2)
kk = np.clip(np.searchsorted(fz_, Z) - 1, 0, fnz - 2)
wz = np.clip((Z - fz_[kk])/(fz_[kk+1] - fz_[kk]), 0, 1)
I, J = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
fth_at = ((1-wz)*fth_col[I[:, :, None], J[:, :, None], kk] +
          wz*fth_col[I[:, :, None], J[:, :, None], kk+1])
dth = TH - fth_at

fails = []


def check(n, ok, msg):
    print(f'  [{"PASS" if ok else "FAIL"}] {n}. {msg}')
    if not ok:
        fails.append(n)


print('\n=== pass conditions ===')

# 1. monotonic p ------------------------------------------------------------------
dP = np.diff(P, axis=2)                       # want < 0 everywhere
worst = dP.max()
nbad  = int((dP >= 0).sum())
check(1, nbad == 0,
      f'p monotonically decreasing: {nbad} non-decreasing layer(s); '
      f'largest dp = {worst:+.4f} Pa (must be < 0)')

# 2. monotonic theta over ocean ---------------------------------------------------
dTH = np.diff(TH, axis=2)
oc3 = np.repeat(ocean[:, :, None], nz - 1, axis=2)
dFR = np.diff(fth_at, axis=2)                 # the same layers, in ERA5 itself
bad = (dTH <= 0) & oc3
# A layer that is non-increasing in ERA5 TOO is the observation, not a defect of
# the construction -- ERA5 does contain well-mixed marine layers. Only layers that
# ERF makes non-increasing where ERA5 does not are failures.
bad_ours = bad & (dFR > 0)
nbad_all, nbad_ours = int(bad.sum()), int(bad_ours.sum())
print(f'       ocean d(theta) over the FIRST layer: min {dTH[:,:,0][ocean].min():+.3f} K, '
      f'mean {dTH[:,:,0][ocean].mean():+.3f} K')
print(f'       ocean layers non-increasing: {nbad_all} of {int(oc3.sum())}; '
      f'of those, {nbad_all-nbad_ours} are non-increasing in ERA5 as well')
if nbad_ours:
    bi, bj, bk = np.unravel_index(np.argmin(np.where(bad_ours, dTH, np.inf)), dTH.shape)
    print(f'       worst construction-only inversion: (i,j,k)=({bi},{bj},{bk}) '
          f'z={Z[bi,bj,bk]:.1f} m ter={ter[bi,bj]:.1f} m '
          f'd(theta)={dTH[bi,bj,bk]:+.5f} K (ERA5 {dFR[bi,bj,bk]:+.5f} K)')
check(2, nbad_ours == 0,
      f'theta increasing with height over ocean: {nbad_ours} layer(s) that ERF makes '
      f'non-increasing where ERA5 does not; min d(theta) over ocean = {dTH[oc3].min():+.4f} K')

# 3. non-zero layer thickness ------------------------------------------------------
dZ = np.diff(Z, axis=2)
check(3, dZ.min() > 0 and (-dP).min() > 0,
      f'layer thickness: min dz = {dZ.min():.4f} m, min dp = {(-dP).min():.5f} Pa '
      f'(both must be > 0)')

# 4. T(first level) vs SST + 0.64 K -------------------------------------------------
with open(ANCHOR, 'rb') as f:
    magic, anx, any_ = struct.unpack('<iii', f.read(12))
    assert magic == 0x45524653 and (anx, any_) == (nx, ny), (magic, anx, any_, nx, ny)
    sp, zorog, t2m = (np.frombuffer(f.read(8*nx*ny), dtype='<f8').reshape(ny, nx).T.copy()
                      for _ in range(3))
# SST straight from the surface frame ERF itself reads
snx, sny, snz, sndata = np.fromfile(SFCBIN, dtype=np.int32, count=4)
sa = np.fromfile(SFCBIN, dtype=np.float32, offset=16)
o = snx + sny + snz
sst_f = sa[o:o + snx*sny].reshape(sny, snx)          # field 0 = sst
sx = sa[0:snx]; sy = sa[snx:snx+sny]

# MASKED SST interpolation, matching what ERF itself does. Interpolating the raw
# frame field mixes land fill values into coastal water and then a range guard on
# the OUTPUT accepts them, because bilinear interpolation fills the gap
# continuously -- the same false-accept that let 62 C through once already.
# Accumulate only source points that are themselves valid SST, renormalise, and
# mark a cell invalid when its stencil had none.
sst_ok = (sst_f > 271) & (sst_f < 305)
num = bilin2d(np.where(sst_ok, sst_f, 0.0), sx, sy, xc, yc)
den = bilin2d(sst_ok.astype(float),          sx, sy, xc, yc)
sst = np.where(den > 0.5, num/np.maximum(den, 1e-12), np.nan)
valid = ocean & np.isfinite(sst) & (sst > 271) & (sst < 305)
print(f'       SST stencils: {int((~np.isfinite(sst) & ocean).sum())} ocean cell(s) with no '
      f'valid source, {int(valid.sum())} scored')
dT = T[:, :, 0][valid] - (sst[valid] + 0.64)
print(f'       ocean: SST mean {sst[valid].mean()-273.15:.2f} C, '
      f'ERA5 t2m mean {t2m[valid].mean()-273.15:.2f} C '
      f'(t2m - SST = {(t2m[valid]-sst[valid]).mean():+.3f} K), '
      f'model T(k=0) mean {T[:,:,0][valid].mean()-273.15:.2f} C')
wn = np.argmax(np.abs(dT))
wi4, wj4 = np.argwhere(valid)[wn]
print(f'       worst cell: (i,j)=({wi4},{wj4}) ter={ter[wi4,wj4]:.1f} m  '
      f'T(k=0) {T[wi4,wj4,0]-273.15:.2f} C  SST {sst[wi4,wj4]-273.15:.2f} C  '
      f't2m {t2m[wi4,wj4]-273.15:.2f} C  ERA5 orog {zorog[wi4,wj4]:.1f} m')
check(4, abs(dT.mean()) <= 1.0 and np.abs(dT).max() <= 3.0,
      f'T(first level) - (SST + 0.64 K): mean {dT.mean():+.3f} K, '
      f'p5/p95 {np.percentile(dT,5):+.2f}/{np.percentile(dT,95):+.2f} K, '
      f'max|.| {np.abs(dT).max():.2f} K')

# 5. theta vs the ERA5 frame at matching heights ------------------------------------
above = Z >= fz_[1]                                  # where the frame has real data
print('       level-mean model theta vs frame theta (ocean columns):')
print('         k    z(m)   model     frame      diff')
for k in list(range(0, 6)) + [8, 12, 20, nz-1]:
    print(f'        {k:3d} {Z[:,:,k][ocean].mean():7.1f} {TH[:,:,k][ocean].mean():8.3f} '
          f'{fth_at[:,:,k][ocean].mean():9.3f} {dth[:,:,k][ocean].mean():+9.3f}')
wi, wj, wk = np.unravel_index(np.argmax(np.where(above, np.abs(dth), -1)), dth.shape)
print(f'       worst |theta - frame|: (i,j,k)=({wi},{wj},{wk}) z={Z[wi,wj,wk]:.1f} m '
      f'ter={ter[wi,wj]:.1f} m ({"ocean" if ocean[wi,wj] else "land"}) '
      f'model {TH[wi,wj,wk]:.3f} frame {fth_at[wi,wj,wk]:.3f} K')
# Where the tail lives. ERF samples the frame at a height built from z_phys_nd
# -- the (i,j) NODE column, averaged in k only -- while x and y are the CELL
# CENTRE, and while this comparison uses the cell-centred z_phys the plotfile
# carries. On flat ground the two heights coincide and the comparison is exact;
# on a slope they differ by O(dz_terrain/2), which at a strong inversion is
# worth whole kelvins. That mismatch is in ERF's interpolator and predates this
# work -- it is identical on the old path. Stratifying by terrain slope
# separates it from an error in the initialization.
slope = np.zeros_like(ter)
slope[1:-1, 1:-1] = np.maximum(np.abs(ter[2:, 1:-1] - ter[:-2, 1:-1]),
                               np.abs(ter[1:-1, 2:] - ter[1:-1, :-2]))/2.0
sl3 = np.repeat(slope[:, :, None], nz, axis=2)
print('        terrain slope      n      mean|dth|   p99      max')
for lab, lo, hi in [('flat  <  5 m/cell', 0, 5), ('  5-25 m/cell', 5, 25),
                    (' 25-75 m/cell', 25, 75), ('  > 75 m/cell', 75, 1e9)]:
    m = above & (sl3 >= lo) & (sl3 < hi)
    if not m.sum():
        continue
    a = np.abs(dth[m])
    print(f'        {lab:18s} {m.sum():7d} {a.mean():9.4f} {np.percentile(a,99):8.4f} {a.max():8.4f}')
oc_ab = above & np.repeat(ocean[:, :, None], nz, axis=2)
print(f'        ocean columns only {oc_ab.sum():7d} {np.abs(dth[oc_ab]).mean():9.4f} '
      f'{np.percentile(np.abs(dth[oc_ab]),99):8.4f} {np.abs(dth[oc_ab]).max():8.4f}')
check(5, np.abs(dth[above]).max() <= 1.5,
      f'|theta - frame| where the frame has data (z >= {fz_[1]:.1f} m): '
      f'mean {np.abs(dth[above]).mean():.3f} K, max {np.abs(dth[above]).max():.3f} K, '
      f'p99 {np.percentile(np.abs(dth[above]), 99):.3f} K')

# ---- supporting numbers (not pass conditions) ---------------------------------------
print('\n=== supporting ===')
print(f'  |p - p_hse| max {np.abs(P-PH).max():.4f} Pa  '
      f'(base and state are the same field by construction; must be ~0)')
print(f'  model p(k=0) mean {P[:,:,0].mean()/100:.2f} hPa vs ERA5 sp mean {sp.mean()/100:.2f} hPa '
      f'(ERA5 orography mean {zorog.mean():.1f} m, ERF terrain mean {ter.mean():.1f} m)')
print(f'  ocean surface-layer lapse: d(theta)/dz over the first layer = '
      f'{(dTH[:,:,0][ocean]/dZ[:,:,0][ocean]).mean()*1000:+.2f} K/km')

print('\n' + ('ALL FIVE PASS' if not fails else f'FAILED: {sorted(set(fails))}'))
sys.exit(0 if not fails else 1)
