#!/bin/bash
# 48-h segment-production gate: leg 1 fresh 0-18h (clean stop before the fatal
# rotation), leg 2 restart 18-48h crossing ten further rotations. Mass and
# stability are the read-outs.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh

"$PF" || { echo "PREFLIGHT FAIL leg1"; exit 3; }
echo "=== LEG1 $(date -u +%H:%M:%S) ==="
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_48h erf-hindcast bash -c '
  /app/ERF/build/Exec/erf_exec inputs_leg1' > "$SP/run48_leg1.full" 2>&1
E1=$?
L1=$(grep -E "Coarse STEP .* ends" "$SP/run48_leg1.full" | tail -1 | sed 's/Coarse //')
CHK=$(ls -d /home/atyagi/ERF/run_48h/chk* 2>/dev/null | sort -V | tail -1)
FPE1=$(grep -c "Erroneous arithmetic" "$SP/run48_leg1.full")
echo "LEG1 exit=$E1 FPE=$FPE1 last=[$L1] chk=$CHK"
[ "$E1" -eq 0 ] && [ -n "$CHK" ] || { echo "LEG1 FAILED - no leg2"; exit 4; }

"$PF" || { echo "PREFLIGHT FAIL leg2"; exit 3; }
echo "=== LEG2 $(date -u +%H:%M:%S) restart=$(basename $CHK) ==="
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_48h erf-hindcast bash -c "
  /app/ERF/build/Exec/erf_exec inputs_leg2 amr.restart=$(basename $CHK)" > "$SP/run48_leg2.full" 2>&1
E2=$?
L2=$(grep -E "Coarse STEP .* ends" "$SP/run48_leg2.full" | tail -1 | sed 's/Coarse //')
FPE2=$(grep -c "Erroneous arithmetic" "$SP/run48_leg2.full")
ROT=$(grep -cE "Reading weather data" "$SP/run48_leg2.full")
M0=$(grep -m1 " MASS       =" "$SP/run48_leg2.full" | awk '{print $3}')
MN=$(grep " MASS       =" "$SP/run48_leg2.full" | tail -1 | awk '{print $3}')
echo "LEG2 exit=$E2 FPE=$FPE2 rotations=$ROT last=[$L2] M:$M0->$MN"
echo RUN48_DONE
