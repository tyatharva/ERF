#!/bin/bash
# Acquire a full year of ERA5 frames as 12 independent monthly segments.
#
#   Exec/CanonicalTests/ChannelIslands/pull_year.sh download 2023
#   Exec/CanonicalTests/ChannelIslands/pull_year.sh process  2023
#
# Two phases on purpose: `download` is network/CDS-queue bound and can run
# alongside GPU work; `process` is an 8-rank CPU job and should not.
#
# Each segment is a FRESH-INIT run (see RUNBOOK "Month-segment launch
# recipe"), so each one gets its own directory, its own frames starting
# exactly at its own start_datetime, and its own soil anchor. Segment M
# starts SPINUP_DAYS before month M and ends at the first instant of month
# M+1; the lead is discarded in post.
#
# CDS credentials are mounted read-only at run time and are never copied
# into the image.
set -euo pipefail

PHASE=${1:?usage: pull_year.sh download-or-process YEAR}
YEAR=${2:?usage: pull_year.sh download-or-process YEAR}
SPINUP_DAYS=2
AREA="36.0,-123.25,31.25,-115.25"
ROOT=/home/atyagi/ERF/era5_year_$YEAR
CI=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands
mkdir -p "$ROOT"

# start/stop/hours for each month, computed once
mapfile -t SEGMENTS < <(python3 - "$YEAR" "$SPINUP_DAYS" <<'PY'
import sys
from datetime import datetime, timedelta
year, lead = int(sys.argv[1]), int(sys.argv[2])
for m in range(1, 13):
    start = datetime(year, m, 1) - timedelta(days=lead)
    stop  = datetime(year + (m == 12), (m % 12) + 1, 1)
    print(f"{m:02d} {start:%Y-%m-%d} {stop:%Y-%m-%d} {int((stop-start).total_seconds()//3600)}")
PY
)

for seg in "${SEGMENTS[@]}"; do
    read -r MM START STOP HOURS <<<"$seg"
    DIR="$ROOT/$YEAR-$MM"
    mkdir -p "$DIR"

    if [ "$PHASE" = download ]; then
        echo "=== $YEAR-$MM: download $START -> $STOP ($HOURS h) ==="
        cp "$CI/era5_batch_download.py" "$DIR/"
        docker run --rm -v /home/atyagi/ERF:/app/ERF -v ~/.cdsapirc:/root/.cdsapirc:ro \
            -w "/app/ERF/era5_year_$YEAR/$YEAR-$MM" erf-hindcast \
            python3 era5_batch_download.py --start "$START" --end "$STOP" --area "$AREA" \
            2>&1 | tail -5
        echo "    frames on disk: $(ls "$DIR"/era5_3d_*.grib 2>/dev/null | wc -l) 3D / $(ls "$DIR"/era5_surf_*.grib 2>/dev/null | wc -l) surf"

    elif [ "$PHASE" = process ]; then
        echo "=== $YEAR-$MM: process ($HOURS h from $START) ==="
        printf 'year: %s\nmonth: %s\nday: %s\ntime: 00:00\narea: %s\n' \
            "${START:0:4}" "${START:5:2}" "${START:8:2}" "$AREA" > "$DIR/era5_input.txt"
        docker run --rm -v /home/atyagi/ERF:/app/ERF -v ~/.cdsapirc:/root/.cdsapirc:ro \
            -w "/app/ERF/era5_year_$YEAR/$YEAR-$MM" erf-hindcast bash -lc "
                cp /opt/erftools/notebooks/era5/WriteICFromERA5Data.py . &&
                cp -r /opt/erftools/notebooks/gfs/TypicalAtmosphereData . &&
                mpirun -n 8 python3 WriteICFromERA5Data.py era5_input.txt \
                    --do_forecast=true --forecast_time_hours=$HOURS --interval_hours=3" 2>&1 | tail -5
        echo "    bins: $(ls "$DIR"/Output/ERA5Data_3D/*.bin 2>/dev/null | wc -l) 3D / $(ls "$DIR"/Output/ERA5Data_Surface/*.bin 2>/dev/null | wc -l) surf"
    else
        echo "unknown phase: $PHASE" >&2; exit 2
    fi
done
echo "PHASE $PHASE COMPLETE for $YEAR"
