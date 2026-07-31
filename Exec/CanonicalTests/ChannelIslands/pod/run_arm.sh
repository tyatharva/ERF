#!/bin/bash
# Launch ONE production arm on ONE GPU, preflighted.
#
#   pod/run_arm.sh <gpu> <run_dir> [extra erf overrides ...]
#
# The baseline every 2020-12-28 arm shares is applied here so the four arms
# cannot drift apart by a typo; per-arm knobs come in as trailing arguments and
# are appended AFTER the baseline, so a repeated key wins.
#
# One 23-h job per card -- the preflight enforces it for the TARGET card only,
# which is what makes four concurrent arms legal. Score timing only from solo
# runs; these are not solo.
set -uo pipefail

GPU=${1:?usage: run_arm.sh <gpu> <run_dir> [overrides...]}
RUN=${2:?usage: run_arm.sh <gpu> <run_dir> [overrides...]}
shift 2

ERF_ROOT=${ERF_ROOT:-/app/ERF}
CI=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands
BIN=$ERF_ROOT/build/Exec/erf_exec
DEST=$ERF_ROOT/$RUN

die () { printf '\nFATAL: %s\n' "$*" >&2; exit 1; }

[ -x "$BIN" ]  || die "no binary at $BIN"
[ -d "$DEST" ] || die "$DEST not staged -- run pod/stage_arm.sh first"
[ -e "$DEST/run.log" ] && die "$DEST/run.log exists -- refusing to clobber a result"

bash "$CI/gpu_preflight.sh" "$GPU" || die "preflight refused GPU $GPU"

BASELINE=(
    erf.cfl=0.3
    erf.moisture_model=Morrison
    erf.les_type=None
    erf.hindcast_mass_du_max=8
    erf.nscbc_lateral=1
    erf.nscbc_outflow=1
    erf.nscbc_parts=31
    erf.nscbc_mass_tau=60
    erf.nscbc_keep_ramp=0
)

cd "$DEST" || die "cannot cd $DEST"
{
    echo "=== arm $RUN on GPU $GPU ==="
    echo "binary : $BIN ($(stat -c %y "$BIN"))"
    echo "commit : $(git -C "$ERF_ROOT" rev-parse --short HEAD)"
    echo "cmd    : CUDA_VISIBLE_DEVICES=$GPU $BIN inputs_c404 ${BASELINE[*]} $*"
} > ARM_PROVENANCE.txt

CUDA_VISIBLE_DEVICES=$GPU "$BIN" inputs_c404 "${BASELINE[@]}" "$@" > run.log 2>&1
rc=$?
echo "EXIT=$rc" >> run.log
echo "[run_arm] $RUN finished on GPU $GPU with EXIT=$rc"
exit $rc
