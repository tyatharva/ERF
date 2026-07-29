#!/bin/bash
# Pick a SAFE restart checkpoint from a run directory (UPSTREAM_ISSUES 31).
#
# A run that ends on stop_datetime can take one final degenerate step of ~1e-08 s
# (the stop-time residual), and the end-of-run checkpoint is written from it,
# storing that dt. Restarting from such a checkpoint is fatal: dt is rate-limited
# to change_max per step, so it needs ~180 steps to climb back to O(1) and the
# run goes non-finite after about five. Whether it happens depends on where the
# truncated step lands, so it cannot be predicted -- it has to be checked.
#
# ERF's checkpoint Header is text and fixed-layout:
#   line 7 = istep, line 8 = dt, line 9 = time.
# This picks the NEWEST checkpoint whose dt is within 1e-3 of the largest dt seen
# across the run's checkpoints -- the fallback rule stated in item 31.
#
# Usage:  pick_restart_chk.sh <rundir>       -> prints the checkpoint name
#         pick_restart_chk.sh <rundir> -v    -> plus a per-checkpoint listing
set -u
RUN=${1:?usage: pick_restart_chk.sh <rundir> [-v]}
VERBOSE=${2:-}

hdr () { sed -n "${2}p" "$1/Header" | tr -d ' '; }

mapfile -t CHKS < <(ls -d "$RUN"/chk[0-9]* 2>/dev/null | sort -V)
[ ${#CHKS[@]} -gt 0 ] || { echo "no checkpoints in $RUN" >&2; exit 2; }

REF=0
for c in "${CHKS[@]}"; do
  [ -r "$c/Header" ] || continue
  REF=$(awk -v a="$REF" -v b="$(hdr "$c" 8)" 'BEGIN{print (b>a)?b:a}')
done
[ "$(awk -v r="$REF" 'BEGIN{print (r>0)?1:0}')" = 1 ] || { echo "no usable dt in $RUN" >&2; exit 2; }

# Select by TIME, not by step number. A run directory accumulates checkpoints
# from several legs, and after a restart the step count for a given model time
# differs between legs -- run_ab_cfl03 holds chk31396 at t=64800.8 (leg 2) and
# chk31415 at t=64768 (leg 1), so step order is NOT time order.
PICK=""; PICKT=-1
for c in "${CHKS[@]}"; do
  [ -r "$c/Header" ] || continue
  d=$(hdr "$c" 8); t=$(hdr "$c" 9)
  ok=$(awk -v d="$d" -v r="$REF" 'BEGIN{print (d >= 1e-3*r)?1:0}')
  if [ "$ok" = 1 ] && [ "$(awk -v a="$t" -v b="$PICKT" 'BEGIN{print (a>b)?1:0}')" = 1 ]; then
    PICK=$c; PICKT=$t
  fi
  [ -n "$VERBOSE" ] && printf '%-12s t=%-14s dt=%-26s %s\n' \
      "$(basename "$c")" "$t" "$d" \
      "$([ "$ok" = 1 ] && echo OK || echo "REJECT: dt is $(awk -v d="$d" -v r="$REF" 'BEGIN{printf "%.0f", r/d}')x below the run maximum")"
done

[ -n "$PICK" ] || { echo "every checkpoint in $RUN has a degenerate dt" >&2; exit 3; }
basename "$PICK"
