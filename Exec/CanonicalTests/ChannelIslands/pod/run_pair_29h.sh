#!/bin/bash
# Run ONE arm of the 29-h anelastic Davies-vs-NSCBC pair, in two legs.
#
#   pod/run_pair_29h.sh davies run_A29_davies
#   pod/run_pair_29h.sh nscbc  run_A29_nscbc
#
#   leg 1  12/27 18Z -> 12/28 00Z    6 h spin-up, ends on a frame instant
#   leg 2  restart    -> 12/28 23Z   23 h scored (matches BOTH references,
#                                    which are 23-h accumulations ending 23Z)
#
# THE TWO ARMS MUST DIFFER IN EXACTLY ONE THING: the lateral boundary scheme.
# That is structural here, not a promise -- every shared setting lives in the
# DECK, both arms read the SAME deck file, and the only command-line arguments
# are the ARM[] flags below. There is no SHARED[] array to drift.
#
# In particular BOTH arms are anelastic at cfl 0.1. Davies tolerates 0.2 and
# NSCBC does not (3 of 6 runs NaN'd), and anelastic is unstable above ~0.10 --
# so the pair runs at the CFL the stricter member needs. Handicapping the
# tolerant arm is the price of a controlled comparison.
set -uo pipefail

SCHEME=${1:?usage: run_pair_29h.sh <davies|nscbc> <run_dir>}
RUN=${2:?usage: run_pair_29h.sh <davies|nscbc> <run_dir>}

ERF_ROOT=${ERF_ROOT:-/app/ERF}
CI=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands
BIN=$ERF_ROOT/build/Exec/erf_exec
DATA=${DATA:-$ERF_ROOT/pod_data_A_run30h}
DECK=${DECK:-inputs_c404_domA}
TERRAIN=${TERRAIN:-terrain_3km_192x96_domA.txt}
DEST=$ERF_ROOT/$RUN

LEG1_STOP="2020-12-28 00:00:00"
LEG2_STOP="2020-12-28 23:00:00"

die () { printf '\nFATAL: %s\n' "$*" >&2; exit 1; }
say () { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }

case "$SCHEME" in
  davies) ARM=(erf.nscbc_lateral=0 erf.nscbc_outflow=0) ;;
  nscbc)  ARM=(erf.nscbc_lateral=1 erf.nscbc_outflow=1 erf.nscbc_parts=31
               erf.nscbc_sigma=0.03 erf.nscbc_mass_tau=60) ;;
  *) die "scheme must be davies or nscbc" ;;
esac

[ -x "$BIN" ] || die "no binary at $BIN"
[ -e "$DEST" ] && die "$DEST exists -- run dirs are results, not scratch"
for f in CONUS404Data_3D CONUS404Data_Surface; do
    [ -d "$DATA/$f" ] || die "missing $DATA/$f"
done

# ---------------------------------------------------------------- stage -----
say "stage $RUN ($SCHEME)"
mkdir -p "$DEST" && cd "$DEST" || die "cannot create $DEST"
ln -sf "$DATA/CONUS404Data_3D" .
ln -sf "$DATA/CONUS404Data_Surface" .
cp "$CI/$TERRAIN" .        || die "no terrain $CI/$TERRAIN"
cp "$CI/$DECK" inputs_c404 || die "no deck $CI/$DECK"
grep -q "erf.terrain_file_name = \"$TERRAIN\"" inputs_c404 \
    || die "deck does not name the staged terrain $TERRAIN"
ROUGH=$(grep -oP '(?<=^erf.most.roughness_file_name = ")[^"]+' inputs_c404 || true)
[ -n "$ROUGH" ] && { cp "$CI/$ROUGH" . || die "missing roughness map $CI/$ROUGH"; }
for t in rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/rrtmgp/data/$t" . || die "missing $t"
done
for t in rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/extensions/cloud_optics/$t" . || die "missing $t"
done

# Refuse to run a pair that is not actually a pair.
# The pair is COMPRESSIBLE: NSCBC's LODI relations are built on an acoustic
# sound speed that anelastic does not have (see the deck). Refuse the invalid
# pairing rather than let it run and NaN at a wall 23 model minutes in.
grep -q '^erf.anelastic = 0' inputs_c404 \
    || die "deck is anelastic: NSCBC is a compressible-equations BC and fails at a wall"
echo "deck sha256: $(sha256sum inputs_c404 | cut -c1-16)  <- must match the other arm"
echo "arm flags  : ${ARM[*]}"

# ------------------------------------------------------------ run_leg -------
# Run one leg, restarting from a checkpoint if it dies.
#
# Two distinct failure modes are known and they need different responses:
#   * NONDETERMINISTIC (measured today at max_dt 0.7): the same state can run
#     clean on a second attempt, so retrying the SAME checkpoint is worthwhile.
#   * PATH-DEPENDENT (UPSTREAM_ISSUES item 66): restarting from the checkpoint
#     AT the failure reproduced it exactly, and only an EARLIER one cleared it.
# So: retry the latest checkpoint once, then walk backwards. That is also why
# amr.check_per is hourly -- each step back costs at most one model hour.
run_leg () {                       # $1 = leg label, $2 = log, $3 = optional restart
    local label=$1 log=$2 restart=${3:-}
    local attempt=0 back=0
    while [ $attempt -lt 5 ]; do
        attempt=$((attempt+1))
        local args=(inputs_c404 "${ARM[@]}")
        [ -n "$restart" ] && args+=(amr.restart="$restart")
        echo "  attempt $attempt${restart:+ (restart from $restart)}"
        "$BIN" "${args[@]}" > "$log" 2>&1
        local rc=$?
        if [ $rc -eq 0 ] && ! grep -q 'contains NaNs' "$log"; then
            return 0
        fi
        grep -q 'contains NaNs' "$log" && echo "  $label: NaN" || echo "  $label: exit $rc"
        cp "$log" "$log.attempt$attempt"
        # Choose the checkpoint to resume from, walking back one further each time.
        mapfile -t CHKS < <(ls -d chk* 2>/dev/null | grep -E 'chk[0-9]+$' | sort)
        local n=${#CHKS[@]}
        [ $n -eq 0 ] && { tail -20 "$log"; die "$label failed with no checkpoint to resume from"; }
        back=$((back+1))
        local idx=$((n-back))
        [ $idx -lt 0 ] && { tail -20 "$log"; die "$label: exhausted checkpoints"; }
        restart=${CHKS[$idx]}
        echo "  resuming from $restart"
    done
    die "$label: 5 attempts exhausted"
}

# ---------------------------------------------------------------- leg 1 -----
say "leg 1 -> $LEG1_STOP"
sed -i "s|^stop_datetime  = .*|stop_datetime  = \"$LEG1_STOP\"|" inputs_c404
run_leg "leg 1" leg1.log

CHK=$(ls -d chk* 2>/dev/null | grep -E 'chk[0-9]+$' | sort | tail -1)
[ -n "$CHK" ] || die "leg 1 wrote no checkpoint"
echo "leg 1 done, restarting from $CHK"

# ---------------------------------------------------------------- leg 2 -----
# start_datetime is deliberately NOT touched -- model time continues across the
# restart, so changing it would reinterpret every frame index. ERF now aborts if
# it disagrees with the checkpoint, so this is enforced rather than trusted.
say "leg 2 -> $LEG2_STOP (restart from $CHK)"
sed -i "s|^stop_datetime  = .*|stop_datetime  = \"$LEG2_STOP\"|" inputs_c404
run_leg "leg 2" leg2.log "$CHK"

say "ARM COMPLETE: $RUN ($SCHEME)"
grep -a 'Coarse STEP' leg2.log | tail -1
echo "plotfiles: $(ls -d plt* 2>/dev/null | wc -l)"
