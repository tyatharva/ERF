#!/bin/bash
# Wait for frame generation, GATE on the frame audit, then stage and launch
# both 71-h arms. Autonomous: the audit is a hard gate, not a report.
set -uo pipefail
cd /app/ERF
CI=Exec/CanonicalTests/ChannelIslands
D=pod_data_A_lead

log () { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }

log "waiting for frame generation"
while pgrep -f conus404_to_bin >/dev/null || \
      [ "$(ls $D/CONUS404Data_Surface/*.bin 2>/dev/null | wc -l)" -lt 25 ]; do
    sleep 30
done
n3=$(ls $D/CONUS404Data_3D/*.bin | wc -l)
ns=$(ls $D/CONUS404Data_Surface/*.bin | wc -l)
log "generation finished: 3D=$n3 SFC=$ns"
[ "$n3" -eq 25 ] && [ "$ns" -eq 25 ] || { log "ABORT: expected 25 frames each"; exit 1; }

for k in 3D:3d Surface:sfc; do
    dir=${k%%:*}; kind=${k##*:}
    out=$(python3 $CI/scoring/frame_audit.py $D/CONUS404Data_$dir $kind 2>&1)
    echo "$out"
    grep -q "VERDICT: PASS" <<<"$out" || { log "ABORT: $dir frames failed audit"; exit 1; }
done
log "frame audit PASSED both sets"

export DATA=/app/ERF/$D
export DECK=inputs_c404_domA_71h
export TERRAIN=terrain_3km_192x96_domA.txt
for r in run_A71_nscbc run_A71_davies; do
    bash $CI/pod/stage_arm.sh $r || { log "ABORT: staging $r failed"; exit 1; }
done

bash $CI/pod/run_arm71.sh nscbc  "0,1" run_A71_nscbc  || { log "ABORT: nscbc launch"; exit 1; }
sleep 20
bash $CI/pod/run_arm71.sh davies "2,3" run_A71_davies || { log "ABORT: davies launch"; exit 1; }
log "both arms launched"
