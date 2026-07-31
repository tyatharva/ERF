#!/bin/bash
# Run the full reporting battery for ONE landed arm.
#
#   pod/analyze_arm.sh <label> <run_dir>
#
# Emits, in the order the campaign asks for them:
#   1. inflow_flux.py     hourly ERF/driver flux ratio at both faces -- and the
#                         ylo sign, which is the item-51 question; plus sub-1500 m
#                         wind over the ranges vs d02, and column vapour
#   2. spinup_split.py    hour-by-hour and cumulative >500 m bias
#   3. terrain_windward.py  terrain-stratified bias, Santa Ynez, domain max
#   4. score_c404.py      accumulation battery vs d02 AND vs MRMS
#   5. hourly_series_arms.py  hourly domain-mean rate vs d02 and MRMS
#
# NO WALL-CLOCK CONCLUSIONS. These arms ran four-to-a-box; fields, flux ratios
# and scores are valid, timing is not.
set -uo pipefail

LAB=${1:?usage: analyze_arm.sh <label> <run_dir>}
RUN=${2:?usage: analyze_arm.sh <label> <run_dir>}
ERF_ROOT=${ERF_ROOT:-/app/ERF}
SC=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/scoring
DEST=$ERF_ROOT/$RUN
OUT=$DEST/ANALYSIS

die () { printf '\nFATAL: %s\n' "$*" >&2; exit 1; }
[ -d "$DEST" ] || die "no run dir $DEST"
mkdir -p "$OUT"

# The final plotfile, chosen by TIME not by name -- plt numbering is by step, so
# lexical order is not time order across runs with different dt histories.
#
# AND IT MUST NOT BE THE stop_datetime TRUNCATION-STEP PLOTFILE. Measured on
# arm A (sigma=1): the run lands on stop_time with a final step of
# dt = 2.98e-08 s, and the plotfile written on that step carries NaN in 51.5%
# of the surface rain_accum cells while every prognostic field is clean. It is
# a diagnostic artifact of the near-zero dt, not state corruption -- but
# rain_accum is the field the whole accumulation battery scores, so using it
# silently produces a wrong bias, a wrong domain max and a NaN at Santa Ynez.
#
# This is item 31 (a checkpoint written on the truncation step poisons every
# restart from it) reaching PLOTFILES. pick_restart_chk.sh guards checkpoints;
# nothing guarded plotfiles. Reject on the data itself rather than on dt, so
# the guard cannot be fooled by a different truncation size.
FINAL=$(python3 - "$DEST" <<'PY'
import glob, os, sys, numpy as np, yt
yt.set_log_level(50)
cand = []
for p in sorted(glob.glob(os.path.join(sys.argv[1], 'plt[0-9]*'))):
    try:
        cand.append((float(yt.load(p).current_time), p))
    except Exception:
        continue
for t, p in sorted(cand, reverse=True):
    ds = yt.load(p)
    r = np.array(ds.covering_grid(0, ds.domain_left_edge,
                                  ds.domain_dimensions)['rain_accum'])
    n = int(np.isnan(r).sum())
    if n == 0:
        print(p)
        break
    print(f'REJECT {p} t={t:.2f}s -- {n} NaN in rain_accum', file=sys.stderr)
PY
)
[ -n "$FINAL" ] || die "no readable plotfile in $DEST"
echo "=== $LAB : final plotfile $FINAL ==="
python3 -c "import yt,sys; yt.set_log_level(50); print('  t =', float(yt.load('$FINAL').current_time), 's =', float(yt.load('$FINAL').current_time)/3600, 'h')"

run_step () {
    local name=$1; shift
    echo; echo "########## $name ##########"
    timeout 7200 python3 "$@" 2>&1 | tee "$OUT/$name.txt" | tail -200
    echo "[wrote $OUT/$name.txt]"
}

run_step inflow_flux    "$SC/inflow_flux.py"    "$LAB=$RUN"
run_step spinup_split   "$SC/spinup_split.py"   "$LAB=$RUN"
run_step terrain_windward "$SC/terrain_windward.py" "$LAB=$FINAL"
run_step score_d02      "$SC/score_c404.py" "$ERF_ROOT/wrf_d02_on_grid.npy"       d02  "$LAB=$FINAL"
run_step score_mrms     "$SC/score_c404.py" "$ERF_ROOT/mrms_20201228_on_grid.npy" MRMS "$LAB=$FINAL"
run_step hourly_series  "$SC/hourly_series_arms.py" "$LAB=$RUN"

echo; echo "=== $LAB analysis complete -> $OUT ==="
