#!/bin/bash
# Refuse to start a GPU experiment if anything else holds the TARGET device.
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
#   gpu_preflight.sh          check every visible GPU (single-card workstation)
#   gpu_preflight.sh 2        check ONLY GPU 2 (one arm per card, four in flight)
#   gpu_preflight.sh --clean  kill stray erf-hindcast containers first
#
# CUDA_VISIBLE_DEVICES, if set to a single index, is honoured as the target, so
#   CUDA_VISIBLE_DEVICES=2 gpu_preflight.sh && CUDA_VISIBLE_DEVICES=2 erf_exec ...
# checks the card the run will actually land on.
#
# Wrap EVERY experimental run:
#   Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh && docker run ... || exit 1
#
# 2026-07-31 -- REWRITTEN FOR MULTI-GPU. The previous version queried
# memory.used with no --id and compared the result with `[ "$used" -gt ... ]`.
# On a multi-GPU host nvidia-smi returns one line PER CARD, so that test got a
# newline-separated string, errored with "integer expression expected", took the
# else branch, and printed "OK" -- unconditionally, whatever the cards held. It
# FAILED OPEN, which is worse than refusing: it is the same silent-failure class
# as the Blackwell arch fallback (HANDOFF 7b) and item 46. It also called
# `docker ps` unguarded; on a pod with no Docker that is `command not found`,
# an empty stray list, and another silent pass.
set -u

BASELINE_MIB=${GPU_BASELINE_MIB:-1200}   # desktop/compositor headroom
fail () { echo "FATAL [gpu_preflight]: $*" >&2; exit 1; }

# ------------------------------------------------------------- target GPU ---
TARGET=""
if [ "${1:-}" = "--clean" ]; then
    shift
    if command -v docker >/dev/null 2>&1; then
        stray=$(docker ps -q --filter ancestor=erf-hindcast)
        if [ -n "$stray" ]; then
            echo "[gpu_preflight] killing stray erf-hindcast containers: $stray"
            docker kill $stray >/dev/null 2>&1
            sleep 3
        fi
    fi
fi
if [ -n "${1:-}" ]; then
    TARGET=$1
elif [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
    case "$CUDA_VISIBLE_DEVICES" in
        *,*) : ;;                       # multi-card job: check all of them
        *)   TARGET=$CUDA_VISIBLE_DEVICES ;;
    esac
fi
case "$TARGET" in
    ''|*[!0-9]*) [ -z "$TARGET" ] || fail "target GPU must be an index, got '$TARGET'" ;;
esac

# ------------------------------------------------------ stray containers ----
# Only meaningful where Docker exists. On a RunPod pod the container IS the
# environment, so say that out loud instead of passing a check that never ran.
if command -v docker >/dev/null 2>&1; then
    stray=$(docker ps -q --filter ancestor=erf-hindcast 2>/dev/null)
    [ -z "$stray" ] || fail "erf-hindcast container(s) STILL RUNNING: $stray
  A previous run was killed by \`timeout\` or Ctrl-C, which kills the client and
  leaves the container holding GPU memory. Re-run with --clean, or:
      docker kill $stray"
    docker_note="no stray erf-hindcast containers"
else
    docker_note="docker absent (native/pod) -- container check N/A, process check below is authoritative"
fi

# ------------------------------------------------------------ device state --
command -v nvidia-smi >/dev/null 2>&1 || fail "nvidia-smi unavailable"

if [ -n "$TARGET" ]; then
    mem=$(nvidia-smi --id="$TARGET" --query-gpu=index,memory.used --format=csv,noheader,nounits 2>&1) \
        || fail "nvidia-smi failed for GPU $TARGET: $mem"
    apps=$(nvidia-smi --id="$TARGET" --query-compute-apps=pid,process_name,used_memory \
                      --format=csv,noheader 2>/dev/null)
else
    mem=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits 2>&1) \
        || fail "nvidia-smi failed: $mem"
    apps=$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
                      --format=csv,noheader 2>/dev/null)
fi
[ -n "$mem" ] || fail "nvidia-smi returned no GPU rows${TARGET:+ for GPU $TARGET}"

# One card per line. NEVER compare the whole blob -- that is the bug this
# rewrite exists to kill.
bad=0
while IFS=, read -r idx used; do
    idx=${idx// /}; used=${used// /}
    case "$used" in
        ''|*[!0-9]*) fail "unparseable memory.used for GPU '$idx': '$used'" ;;
    esac
    if [ "$used" -gt "$BASELINE_MIB" ]; then
        echo "[gpu_preflight] GPU $idx holds ${used} MiB (baseline ${BASELINE_MIB} MiB)" >&2
        bad=1
    else
        echo "[gpu_preflight] GPU $idx: ${used} MiB"
    fi
done <<EOF
$mem
EOF

# A stray holder is a PROCESS. On a headless pod the memory baseline is ~1 MiB,
# so thresholding alone would wave through a 1 GB leftover. Any compute process
# on the target card fails the check regardless of how much it holds.
if [ -n "$apps" ]; then
    if [ -n "$TARGET" ]; then where="GPU $TARGET"; else where="all GPUs"; fi
    echo "--- compute processes on $where ---" >&2
    echo "$apps" >&2
    bad=1
fi

[ "$bad" -eq 0 ] || fail "GPU${TARGET:+ $TARGET} is not clean.
  Anything above the baseline, or any compute process at all, will silently
  distort or invalidate this experiment. Identify and clear the holder first.
  (One 23-h job per card. Never two.)"

echo "[gpu_preflight] OK -- ${TARGET:+GPU $TARGET }clean, $docker_note"
