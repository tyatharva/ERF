#!/bin/bash
# Stage and run ONE month segment of the year hindcast.
#
#   Exec/CanonicalTests/ChannelIslands/launch_segment.sh 2023 08
#
# Each segment is an INDEPENDENT FRESH-INIT run -- never a restart of the
# previous month. The LSM state is checkpointed, so a restart-chained
# segment silently keeps the previous month's soil and erf.mm5.soil_theta
# becomes inert. Restart is for crash recovery WITHIN a segment only:
#   cd run_<YYYY>-<MM> && erf_exec inputs_hindcast erf.restart=chk<last>
#
# ONE SEGMENT PER GPU: two ERF processes do not fit on a 16 GB card.
set -euo pipefail

YEAR=${1:?usage: launch_segment.sh YEAR MM}
MM=${2:?usage: launch_segment.sh YEAR MM}
SPINUP_DAYS=2
ERF=/home/atyagi/ERF
SEG="$YEAR-$MM"

# Domain land-mean ERA5 skin temperature by month (K). Because MM5 is inert
# this IS the land surface temperature for the whole segment -- see RUNBOOK
# "MM5 IS INERT". Refetch for a year other than 2023 (this is then climatology).
SKT=(x 283.2 283.2 284.4 289.2 291.5 293.9 300.1 298.9 296.1 293.9 288.9 286.5)
SOIL=${SKT[$((10#$MM))]}

eval "$(python3 -c "
from datetime import datetime, timedelta
y, m = $YEAR, $((10#$MM))
s = datetime(y, m, 1) - timedelta(days=$SPINUP_DAYS)
e = datetime(y + (m == 12), (m % 12) + 1, 1)
print(f'START=\"{s:%Y-%m-%d %H:%M:%S}\"; STOP=\"{e:%Y-%m-%d %H:%M:%S}\"')")"

echo "=== segment $SEG : $START -> $STOP (lead ${SPINUP_DAYS}d), soil_theta=$SOIL K ==="
[ -d "$ERF/era5_year_$YEAR/$SEG/Output/ERA5Data_3D" ] || {
    echo "FATAL: frames not processed -- run pull_year.sh process $YEAR" >&2; exit 1; }

docker run --rm --gpus all -v "$ERF:/app/ERF" -w /app/ERF erf-hindcast \
    Exec/CanonicalTests/ChannelIslands/stage_run.sh \
    "/app/ERF/era5_year_$YEAR/$SEG/Output" "/app/ERF/run_$SEG" "$START" "$STOP" "$SOIL"

echo "launching $SEG"
docker run --rm --gpus all -v "$ERF:/app/ERF" -w "/app/ERF/run_$SEG" erf-hindcast \
    bash -lc "/app/ERF/build/Exec/erf_exec inputs_hindcast amrex.max_gpu_streams=1 > run.log 2>&1; echo EXIT=\$? >> run.log"
tail -1 "$ERF/run_$SEG/run.log"
