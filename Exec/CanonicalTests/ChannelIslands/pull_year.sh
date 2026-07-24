#!/bin/bash
# Acquire and preprocess a full year of ERA5 BC/IC frames into ONE FLAT POOL.
#
#   Exec/CanonicalTests/ChannelIslands/pull_year.sh download [START END]
#   Exec/CanonicalTests/ChannelIslands/pull_year.sh process  [START END]
#
# Defaults cover calendar 2023 with a ONE-MONTH spin-up cushion in front
# (so the year can also be run as a single continuous job if we ever want
# to) and one spare day at the end:  2022-12-01 -> 2024-01-02.
#
# Two phases on purpose: `download` is CDS-queue bound and safe to run
# alongside GPU work; `process` is an 8-rank CPU job and is not.
#
# WHY A FLAT POOL: erftools is a pure PREPROCESSOR -- WriteICFromERA5Data.py
# globs era5_{3d,surf}_*.grib in its cwd, shards the list across MPI ranks,
# and writes .bin frames to Output/. ERF only ever reads those .bin files at
# run time; nothing is generated on the fly. So the whole year is processed
# ONCE, here, and every later run just points at the result.
#
# NOTE for consumers: ERF indexes boundary frames POSITIONALLY from t=0, so
# a run's frame[0] must BE its start_datetime. A run that does not start at
# the pool's first frame therefore needs a directory holding just its own
# slice -- launch_segment.sh builds that as symlinks into this pool (no
# copying). The preflight in stage_run.sh enforces the invariant.
#
# CDS credentials are mounted read-only at run time, never baked into the image.
set -euo pipefail

PHASE=${1:?usage: pull_year.sh download-or-process [START END]}
START=${2:-2022-12-01}
END=${3:-2024-01-02}
AREA="36.0,-123.25,31.25,-115.25"
POOL=/home/atyagi/ERF/era5_year
CI=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands
mkdir -p "$POOL"

HOURS=$(python3 -c "
from datetime import datetime
a=datetime.strptime('$START','%Y-%m-%d'); b=datetime.strptime('$END','%Y-%m-%d')
print(int((b-a).total_seconds()//3600))")

case "$PHASE" in
download)
    echo "=== download $START -> $END ($HOURS h, $((HOURS/3+1)) frames/stream) ==="
    cp "$CI/era5_batch_download.py" "$POOL/"
    docker run --rm -v /home/atyagi/ERF:/app/ERF -v ~/.cdsapirc:/root/.cdsapirc:ro \
        -w /app/ERF/era5_year erf-hindcast \
        python3 era5_batch_download.py --start "$START" --end "$END" --area "$AREA"
    echo "pool now: $(ls "$POOL"/era5_3d_*.grib 2>/dev/null | wc -l) 3D / $(ls "$POOL"/era5_surf_*.grib 2>/dev/null | wc -l) surf gribs"
    ;;
process)
    # erftools' own downloader skips every GRIB already present, so this is
    # pure processing. Re-runnable: rerunning reprocesses, it does not refetch.
    echo "=== process $START -> $END ($HOURS h) ==="
    printf 'year: %s\nmonth: %s\nday: %s\ntime: 00:00\narea: %s\n' \
        "${START:0:4}" "${START:5:2}" "${START:8:2}" "$AREA" > "$POOL/era5_input.txt"
    docker run --rm -v /home/atyagi/ERF:/app/ERF -v ~/.cdsapirc:/root/.cdsapirc:ro \
        -w /app/ERF/era5_year erf-hindcast bash -lc "
            cp /opt/erftools/notebooks/era5/WriteICFromERA5Data.py . &&
            cp -r /opt/erftools/notebooks/gfs/TypicalAtmosphereData . &&
            mpirun -n 8 python3 WriteICFromERA5Data.py era5_input.txt \
                --do_forecast=true --forecast_time_hours=$HOURS --interval_hours=3"
    echo "pool now: $(ls "$POOL"/Output/ERA5Data_3D/*.bin 2>/dev/null | wc -l) 3D / $(ls "$POOL"/Output/ERA5Data_Surface/*.bin 2>/dev/null | wc -l) surf frames"
    ;;
*)
    echo "unknown phase: $PHASE" >&2; exit 2 ;;
esac
echo "PHASE $PHASE COMPLETE"
