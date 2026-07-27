#!/bin/bash
# Refuse to start a GPU experiment if anything else holds device memory.
#
# WHY THIS EXISTS. A container that survived a `timeout`-killed compute-sanitizer
# run held 7072 MiB for ELEVEN HOURS. `timeout` kills the docker client, not the
# container. Against this model's ~10.6 GB peak that left the machine unable to
# run at all, and three consecutive experiment "results" were environmental
# rather than real -- including a retracted root-cause conclusion
# (UPSTREAM_ISSUES 29c/29d). The tell was there: the failure MODE changed across
# builds on identical inputs, and item 13 already says this fork is not
# bit-reproducible, so a changed failure mode means check the environment, not
# the source. This makes that check structural instead of remembered.
#
#   gpu_preflight.sh          assert clean, else fail loudly
#   gpu_preflight.sh --clean  kill stray erf-hindcast containers first
#
# Wrap EVERY experimental run:
#   Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh && docker run ... || exit 1
set -u

BASELINE_MIB=${GPU_BASELINE_MIB:-1200}   # desktop/compositor headroom
fail () { echo "FATAL [gpu_preflight]: $*" >&2; exit 1; }

if [ "${1:-}" = "--clean" ]; then
    stray=$(docker ps -q --filter ancestor=erf-hindcast)
    if [ -n "$stray" ]; then
        echo "[gpu_preflight] killing stray erf-hindcast containers: $stray"
        docker kill $stray >/dev/null 2>&1
        sleep 3
    fi
fi

stray=$(docker ps -q --filter ancestor=erf-hindcast)
[ -z "$stray" ] || fail "erf-hindcast container(s) STILL RUNNING: $stray
  A previous run was killed by \`timeout\` or Ctrl-C, which kills the client and
  leaves the container holding GPU memory. Re-run with --clean, or:
      docker kill $stray"

used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null) \
    || fail "nvidia-smi unavailable"

if [ "$used" -gt "$BASELINE_MIB" ]; then
    echo "--- processes holding GPU memory ---" >&2
    nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv >&2
    fail "GPU already holds ${used} MiB (baseline ${BASELINE_MIB} MiB).
  Anything above the desktop baseline will silently distort or invalidate this
  experiment. Identify and clear the holder before running."
fi

echo "[gpu_preflight] OK -- ${used} MiB in use, no stray erf-hindcast containers"
