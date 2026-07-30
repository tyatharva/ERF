#!/bin/bash
# Full-day NSCBC-vs-Davies A/B, both arms on the SAME binary (the #26 rad-clock
# fix changes the trajectory, so the #25 Davies baseline cannot be reused).
#
# Each arm uses the validated bracketing recipe: leg 1 cold 0-18 h stopping
# before the fatal [6,7] rotation, leg 2 restart 18-24 h. rain_accum is
# continuous through the checkpoint, which is what the scoring reads.
#
#   arm dav : Davies + hindcast_global_mass_tau=60           (deck default)
#   arm nsc : NSCBC + 19d regime fix + nscbc_mass_tau=60
#             + hindcast_mass_du_max=8 (the ±2 default saturates -- measured
#             +0.834% vs +0.012% over 3000 steps)
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
ARM=$1; shift
EXTRA="$*"
RUN=/home/atyagi/ERF/run_ab_$ARM
CRUN=/app/ERF/run_ab_$ARM

docker run --rm -v /home/atyagi/ERF:/app/ERF erf-hindcast bash -c "
  set -e
  rm -rf $CRUN && mkdir -p $CRUN && cd $CRUN
  cp /app/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast .
  cp /app/ERF/run_a3/channel_islands_terrain_3km_192x96.txt .
  cp /app/ERF/run_a3/rrtmgp-*.nc /app/ERF/run_a3/sfc_anchor_jan09_192x96.bin .
  cp -r /app/ERF/run_a3/ERA5Data_3D /app/ERF/run_a3/ERA5Data_Surface .
  sed 's|^stop_datetime.*|stop_datetime = \"2023-01-09 18:00:00\"|' inputs_hindcast > inputs_leg1
  cp inputs_hindcast inputs_leg2"
for f in inputs_leg1 inputs_leg2 sfc_anchor_jan09_192x96.bin; do
  [ -s "$RUN/$f" ] || { echo "PRECOND: missing $f"; exit 4; }
done
grep -q '2023-01-09 18:00:00' "$RUN/inputs_leg1" || { echo "PRECOND: leg1 stop not set"; exit 4; }
grep -q '2023-01-10 00:00:00' "$RUN/inputs_leg2" || { echo "PRECOND: leg2 stop wrong"; exit 4; }
echo "PRECOND OK ($ARM): extra=[$EXTRA]"

"$PF" || { echo "PREFLIGHT FAIL"; exit 3; }
echo "=== $ARM LEG1 $(date -u +%H:%M:%S) ==="
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w $CRUN erf-hindcast \
  /app/ERF/build/Exec/erf_exec inputs_leg1 $EXTRA > "$SP/ab_${ARM}_leg1.full" 2>&1
E1=$?
# UPSTREAM_ISSUES 31: never take the newest checkpoint blindly -- a run ending on
# stop_datetime can write a final one with a degenerate dt that poisons the restart.
CHK=$RUN/$(/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/pick_restart_chk.sh $RUN)
echo "$ARM LEG1 exit=$E1 last=[$(grep -a 'Coarse STEP .* ends' "$SP/ab_${ARM}_leg1.full" | tail -1)] chk=$(basename ${CHK:-none})"
[ "$E1" -eq 0 ] && [ -n "$CHK" ] || { echo "$ARM LEG1 FAILED"; exit 4; }

"$PF" || { echo "PREFLIGHT FAIL"; exit 3; }
echo "=== $ARM LEG2 $(date -u +%H:%M:%S) restart=$(basename $CHK) ==="
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w $CRUN erf-hindcast \
  /app/ERF/build/Exec/erf_exec inputs_leg2 amr.restart=$(basename $CHK) $EXTRA \
  > "$SP/ab_${ARM}_leg2.full" 2>&1
E2=$?
M0=$(grep -a ' MASS       =' "$SP/ab_${ARM}_leg1.full" | head -1 | awk '{print $3}')
MN=$(grep -a ' MASS       =' "$SP/ab_${ARM}_leg2.full" | tail -1 | awk '{print $3}')
echo "$ARM LEG2 exit=$E2 last=[$(grep -a 'Coarse STEP .* ends' "$SP/ab_${ARM}_leg2.full" | tail -1)]"
echo "$ARM MASS over the day: $M0 -> $MN  $(awk -v a=$M0 -v b=$MN 'BEGIN{printf "%+.3f%%", 100*(b-a)/a}')"
echo "$ARM FPE=$(grep -ac 'Erroneous' "$SP/ab_${ARM}_leg1.full" "$SP/ab_${ARM}_leg2.full" | paste -sd+ | bc)"
echo "AB_${ARM}_DONE"
