#!/usr/bin/env python3
"""Cross-architecture smoke test: does the sm_120 build reproduce sm_89 physics?

  pod_validate.py <plotfile>            check against the embedded reference
  pod_validate.py --emit <plotfile>     print a fresh reference block

WHY THIS EXISTS
  Kokkos 4.5 has no Blackwell architecture; the sm_120 arch entry is OUR patch
  (apply_kokkos_blackwell_patch.sh), and Kokkos 4.5 has never been validated on
  Blackwell by anyone. If the numbers below disagree, nothing computed on this
  machine is trustworthy and no production arm should be started.

WHAT THIS IS AND IS NOT
  This is a SMOKE TEST for gross corruption -- wrong architecture, bad codegen,
  a miscompiled kernel. It is NOT a bitwise reproducibility check. Single
  precision on a different GPU architecture reorders FMAs and reductions, so
  bitwise equality across sm_89 and sm_120 is NOT expected and its absence is
  not a defect.

  Tolerances are therefore physical, not bitwise, and are set well inside the
  margin where a real codegen fault would show. A wrong-architecture build does
  not drift by 1e-4 -- it produces NaNs, zeros, or values wrong in the first
  significant figure.

REFERENCE RUN (regenerate with --emit if the deck or driver data ever change)
  Deck      : Exec/CanonicalTests/ChannelIslands (CONUS404 2020-12-28 case)
  Overrides : max_step=20 amr.check_int=-1 amr.plot_int=20 erf.cfl=0.3
              erf.moisture_model=Morrison erf.les_type=None
              erf.hindcast_mass_du_max=8 erf.nscbc_lateral=1 erf.nscbc_outflow=1
              erf.nscbc_parts=31 erf.nscbc_sigma=0.03 erf.nscbc_mass_tau=60
  Produced  : 2026-07-30 on RTX 4080 (sm_89), single precision, CUDA 12.6.2
  End state : Coarse STEP 20, TIME = 14.31875262 s
"""
import sys
import numpy as np
import yt

yt.set_log_level(50)

# --- expected end-of-step-20 state, measured on sm_89 -----------------------
REFERENCE = {
    #  field          (mean, min, max) at end of step 20, sm_89
    'density':    (0.8643798863924901, 0.11145045608282089, 1.2886161804199219),
    'theta':      (301.2162477212303, 280.9743957519531, 542.8499145507812),
    'qv':         (0.0029642279277506123, 2.9445748617185075e-10, 0.009189499542117119),
    'x_velocity': (0.3237632131482443, -1.824497938156128, 48.07758331298828),
    'y_velocity': (0.14055079262225412, -9.77464485168457, 16.13129425048828),
}
REFERENCE_TIME = 14.31875228881836

# Relative tolerance on each statistic. Generous against SP cross-architecture
# reordering, tight against a codegen fault.
RTOL_MEAN = 1.0e-4
RTOL_EXTREME = 1.0e-3      # min/max are single cells, so noisier than a mean
RTOL_TIME = 1.0e-6         # dt is CFL-derived; a real arch fault moves it far more


def stats(path):
    ds = yt.load(path)
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    out = {}
    for f in REFERENCE:
        a = np.asarray(g[('boxlib', f)], dtype=np.float64)
        out[f] = (float(a.mean()), float(a.min()), float(a.max()))
    return float(ds.current_time), out


def main():
    if sys.argv[1] == '--emit':
        t, s = stats(sys.argv[2])
        print(f'REFERENCE_TIME = {t!r}')
        print('REFERENCE = {')
        for f, (m, lo, hi) in s.items():
            print(f'    {f!r:16s}: ({m!r}, {lo!r}, {hi!r}),')
        print('}')
        return 0

    if any(v[0] is None for v in REFERENCE.values()):
        print('FATAL: reference block is unpopulated. Run with --emit on a '
              'known-good sm_89 plotfile and paste the result in.', file=sys.stderr)
        return 2

    t, s = stats(sys.argv[1])
    bad = []

    dt_rel = abs(t - REFERENCE_TIME) / max(abs(REFERENCE_TIME), 1e-30)
    ok_t = dt_rel <= RTOL_TIME
    print(f'{"model time":14s} {t:>18.9g} vs {REFERENCE_TIME:>18.9g}'
          f'  rel {dt_rel:9.2e}  {"OK" if ok_t else "*** FAIL ***"}')
    if not ok_t:
        bad.append('model time')

    for f, (rm, rlo, rhi) in REFERENCE.items():
        m, lo, hi = s[f]
        for label, got, exp, rtol in (('mean', m, rm, RTOL_MEAN),
                                      ('min', lo, rlo, RTOL_EXTREME),
                                      ('max', hi, rhi, RTOL_EXTREME)):
            denom = max(abs(exp), 1e-30)
            rel = abs(got - exp) / denom
            ok = np.isfinite(got) and rel <= rtol
            print(f'{f:12s} {label:4s} {got:>18.9g} vs {exp:>18.9g}'
                  f'  rel {rel:9.2e}  {"OK" if ok else "*** FAIL ***"}')
            if not ok:
                bad.append(f'{f}.{label}')

    print()
    if bad:
        print('=' * 72)
        print('VALIDATION FAILED:', ', '.join(bad))
        print('The sm_120 build does NOT reproduce the sm_89 reference.')
        print('DO NOT run production arms. Nothing computed here is trustworthy.')
        print('Most likely causes, in order:')
        print('  1. Kokkos built for the wrong architecture (check the CMake cache)')
        print('  2. The Kokkos Blackwell patch did not apply as intended')
        print('  3. A genuine Kokkos 4.5 / sm_120 codegen problem -- this is the')
        print('     scenario the patch header warns about and it is unexplored')
        print('=' * 72)
        return 1

    print('VALIDATION PASSED -- the build under test reproduces the sm_89 reference within tolerance.')
    print('NOTE: this is a 20-step smoke test, not a guarantee over 23 h.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
