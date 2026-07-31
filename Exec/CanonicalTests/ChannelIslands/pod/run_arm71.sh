#!/bin/bash
# Launch ONE 71-h arm of the NSCBC-vs-Davies comparison on TWO GPUs.
#
#   pod/run_arm71.sh nscbc "0,1" run_A71_nscbc
#   pod/run_arm71.sh davies "2,3" run_A71_davies
#
# Sibling of run_arm.sh rather than an edit to it: run_arm.sh's baseline hard-
# codes NSCBC and one card, and the existing 23-h record was produced by it.
#
# The two arms MUST differ in exactly one thing -- the lateral boundary scheme.
# Everything else lives in SHARED below so a typo cannot make them diverge
# silently. Davies is erf.nscbc_lateral=0: ERF_Utils.H:50 says nscbc_lateral
# "replaces the Davies specified+relaxation zone", so the deck's use_real_bcs +
# real_width=10 path IS Davies when NSCBC is off.
set -uo pipefail

SCHEME=${1:?usage: run_arm71.sh <nscbc|davies> <gpus> <run_dir>}
GPUS=${2:?usage: run_arm71.sh <nscbc|davies> <gpus> <run_dir>}
RUN=${3:?usage: run_arm71.sh <nscbc|davies> <gpus> <run_dir>}

ERF_ROOT=${ERF_ROOT:-/app/ERF}
CI=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands
BIN=$ERF_ROOT/build/Exec/erf_exec
DEST=$ERF_ROOT/$RUN

die () { printf '\nFATAL: %s\n' "$*" >&2; exit 1; }

[ -x "$BIN" ]  || die "no binary at $BIN"
[ -d "$DEST" ] || die "$DEST not staged -- run pod/stage_arm.sh first"

# cfl 0.3 matches every scored arm in the campaign and is demonstrated on this
# exact domain and date by run_domA_rot. The deck's own default is 0.2 and its
# comment records blow-ups at 0.25/0.3 on a different (convectively active)
# day. Holding 0.3 keeps continuity with the record; a crash is recoverable
# from an hourly checkpoint, a silent config drift is not.
SHARED=(
    erf.cfl=0.3
    erf.moisture_model=Morrison
    erf.les_type=None
    erf.hindcast_mass_du_max=8
)

case "$SCHEME" in
  nscbc)  ARM=(erf.nscbc_lateral=1 erf.nscbc_outflow=1 erf.nscbc_parts=31
               erf.nscbc_mass_tau=60 erf.nscbc_keep_ramp=0 erf.nscbc_sigma=0.03
               amr.check_per=10800.0) ;;
  # Hourly throughout rather than only from 15/33/51 h: a superset of the
  # requested cadence, ~117 MB per checkpoint, and it makes every restart
  # cheap for the scheme that is expected to need them.
  davies) ARM=(erf.nscbc_lateral=0 erf.nscbc_outflow=0
               amr.check_per=3600.0) ;;
  *) die "scheme must be nscbc or davies, got '$SCHEME'" ;;
esac

NG=$(awk -F, '{print NF}' <<<"$GPUS")
for g in ${GPUS//,/ }; do
    bash "$CI/gpu_preflight.sh" "$g" || die "preflight refused GPU $g"
done

cd "$DEST" || die "cannot cd $DEST"
{
    echo "=== arm $RUN ($SCHEME) on GPUs $GPUS ==="
    echo "binary : $BIN ($(stat -c %y "$BIN"))"
    echo "commit : $(git -C "$ERF_ROOT" rev-parse --short HEAD)"
    echo "amrex  : $(git -C "$ERF_ROOT/Submodules/AMReX" rev-parse --short HEAD)$(
        [ -n "$(git -C "$ERF_ROOT/Submodules/AMReX" status --porcelain)" ] && echo ' DIRTY')"
    echo "cmd    : CUDA_VISIBLE_DEVICES=$GPUS mpiexec -n $NG $BIN inputs_c404 ${SHARED[*]} ${ARM[*]}"
} > ARM_PROVENANCE.txt
cat ARM_PROVENANCE.txt

CUDA_VISIBLE_DEVICES=$GPUS nohup mpiexec -n "$NG" "$BIN" inputs_c404 \
    "${SHARED[@]}" "${ARM[@]}" >> run.log 2>&1 &
echo "launched $RUN pid $! on GPUs $GPUS -> $DEST/run.log"
