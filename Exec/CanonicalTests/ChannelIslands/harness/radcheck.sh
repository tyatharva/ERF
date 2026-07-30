#!/bin/bash
# Verify the #26 rad-clock fix on the running model, against the known-nonzero
# control: the SAME measurement on the pre-fix 48-h leg-2 log gave 255.4 s
# effective for a 180 s request. Anything near 180 here is the fix; anything
# near 256 means it did not take.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
"$PF" || { echo "PREFLIGHT FAIL"; exit 3; }

docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_nan erf-hindcast \
  /app/ERF/build/Exec/erf_exec inputs_hindcast max_step=2400 \
    amr.plot_int=100000 erf.plot_per_1=-1 amr.check_per=-1 erf.sum_interval=1200 \
    amr.check_file=radchk > "$SP/radcheck.full" 2>&1
EX=$?

NCALL=$(grep -ac "Radiation advancing" "$SP/radcheck.full")
TEND=$(grep -a "Coarse STEP .* ends" "$SP/radcheck.full" | tail -1 | awk '{print $7}')
echo "RADCHECK exit=$EX calls=$NCALL model_s=$TEND"
awk -v n="$NCALL" -v t="$TEND" 'BEGIN{ if (n>1) printf "  effective cadence = %.1f s (request 180, pre-fix control 255.4)\n", t/(n-1) }'
echo RADCHECK_DONE
