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
grep -q '^erf.anelastic = 1' inputs_c404 || die "deck is not anelastic"
grep -q '^amrex.fpe_trap_invalid = 0' inputs_c404 \
    || die "anelastic needs amrex.fpe_trap_invalid = 0 (PTX JIT raises FP exceptions)"
echo "deck sha256: $(sha256sum inputs_c404 | cut -c1-16)  <- must match the other arm"
echo "arm flags  : ${ARM[*]}"

# ---------------------------------------------------------------- leg 1 -----
say "leg 1 -> $LEG1_STOP"
sed -i "s|^stop_datetime  = .*|stop_datetime  = \"$LEG1_STOP\"|" inputs_c404
"$BIN" inputs_c404 "${ARM[@]}" > leg1.log 2>&1
rc=$?
grep -q 'contains NaNs' leg1.log && die "leg 1 produced NaNs -- see leg1.log"
[ $rc -eq 0 ] || { tail -20 leg1.log; die "leg 1 exited $rc"; }

CHK=$(ls -d chk* 2>/dev/null | grep -E 'chk[0-9]+$' | sort | tail -1)
[ -n "$CHK" ] || die "leg 1 wrote no checkpoint"
echo "leg 1 done, restarting from $CHK"

# ---------------------------------------------------------------- leg 2 -----
# start_datetime is deliberately NOT touched -- model time continues across the
# restart, so changing it would reinterpret every frame index. ERF now aborts if
# it disagrees with the checkpoint, so this is enforced rather than trusted.
say "leg 2 -> $LEG2_STOP (restart from $CHK)"
sed -i "s|^stop_datetime  = .*|stop_datetime  = \"$LEG2_STOP\"|" inputs_c404
"$BIN" inputs_c404 "${ARM[@]}" amr.restart="$CHK" > leg2.log 2>&1
rc=$?
grep -q 'contains NaNs' leg2.log && die "leg 2 produced NaNs -- see leg2.log"
[ $rc -eq 0 ] || { tail -20 leg2.log; die "leg 2 exited $rc"; }

say "ARM COMPLETE: $RUN ($SCHEME)"
grep -a 'Coarse STEP' leg2.log | tail -1
echo "plotfiles: $(ls -d plt* 2>/dev/null | wc -l)"
