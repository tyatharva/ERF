#!/bin/bash
# Performance-pin probes (post-#25, user-directed):
#  A. amrex.max_gpu_streams unpinned (default) vs pinned baseline -- wall/step + health
#  B. cfl 0.3 and 0.5 probes (4000 steps, covers the cold-start transient + storm ramp)
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh

run () {
  local NAME=$1 MS=$2; shift 2
  "$PF" || { echo "PREFLIGHT FAIL $NAME"; exit 3; }
  local T0=$(date +%s)
  docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_a3 erf-hindcast bash -c "
    /app/ERF/build/Exec/erf_exec inputs_hindcast max_step=$MS $* \
      amr.plot_int=100000 erf.plot_per_1=-1 amr.check_per=-1 erf.sum_interval=200" \
    > "$SP/perf_$NAME.full" 2>&1
  local EX=$?
  local T1=$(date +%s)
  local STEPS=$(grep -cE "Coarse STEP .* ends" "$SP/perf_$NAME.full")
  local LAST=$(grep -E "Coarse STEP .* ends" "$SP/perf_$NAME.full" | tail -1 | awk '{print $7, $10}')
  local BAD=$(grep -cE "Erroneous|SIGABRT|illegal|= -?nan" "$SP/perf_$NAME.full")
  local MEANSTEP=$(grep -oE "Timestep time = [0-9.]+" "$SP/perf_$NAME.full" | awk '{s+=$4; n++} END {if(n)printf "%.4f", s/n}')
  echo "RESULT $NAME: exit=$EX steps=$STEPS wall=$((T1-T0))s meanstep=${MEANSTEP}s bad=$BAD last=[$LAST]"
}

run pin1   2000
run pin4   2000 amrex.max_gpu_streams=4
run cfl03  4000 erf.cfl=0.3
run cfl05  4000 erf.cfl=0.5
echo PERF_DONE
