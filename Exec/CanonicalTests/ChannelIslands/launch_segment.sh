#!/bin/bash
# Stage and run ONE month segment of the year hindcast, from the flat frame
# pool built by pull_year.sh.
#
#   Exec/CanonicalTests/ChannelIslands/launch_segment.sh 2023 08
#   Exec/CanonicalTests/ChannelIslands/launch_segment.sh 2023 08 --stage-only
#
# Each segment is an INDEPENDENT FRESH-INIT run -- never a restart of the
# previous month. The LSM state is checkpointed, so a restart-chained
# segment silently keeps the previous month's soil and erf.mm5.soil_theta
# becomes inert. Restart is for crash recovery WITHIN a segment only:
#   cd run_<YYYY>-<MM> && erf_exec inputs_hindcast erf.restart=chk<last>
#
# ONE SEGMENT PER GPU: two ERF processes do not fit on a 16 GB card.
set -euo pipefail

YEAR=${1:?usage: launch_segment.sh YEAR MM [--stage-only]}
MM=${2:?usage: launch_segment.sh YEAR MM [--stage-only]}
STAGE_ONLY=${3:-}
SPINUP_DAYS=2
ERF=/home/atyagi/ERF
POOL=$ERF/era5_year/Output
SEG="$YEAR-$MM"
RUN="$ERF/run_$SEG"

# Domain land-mean ERA5 skin temperature by month (K). Because MM5 is inert
# this IS the land surface temperature for the whole segment -- see RUNBOOK
# "MM5 IS INERT". Refetch for a year other than 2023 (else it is climatology).
SKT=(x 283.2 283.2 284.4 289.2 291.5 293.9 300.1 298.9 296.1 293.9 288.9 286.5)
SOIL=${SKT[$((10#$MM))]}

eval "$(python3 -c "
from datetime import datetime, timedelta
y, m = $YEAR, $((10#$MM))
s = datetime(y, m, 1) - timedelta(days=$SPINUP_DAYS)
e = datetime(y + (m == 12), (m % 12) + 1, 1)
print(f'START=\"{s:%Y-%m-%d %H:%M:%S}\"; STOP=\"{e:%Y-%m-%d %H:%M:%S}\"')")"

echo "=== segment $SEG : $START -> $STOP (lead ${SPINUP_DAYS}d), soil_theta=$SOIL K ==="
[ -d "$POOL/ERA5Data_3D" ] || { echo "FATAL: frame pool missing -- run pull_year.sh process" >&2; exit 1; }

# ERF indexes frames positionally from t=0, so this segment needs a directory
# whose FIRST frame is its own start_datetime. Symlink the slice; never copy.
mkdir -p "$RUN/ERA5Data_3D" "$RUN/ERA5Data_Surface"
rm -f "$RUN"/ERA5Data_3D/*.bin "$RUN"/ERA5Data_Surface/*.bin
python3 - "$POOL" "$RUN" "$START" "$STOP" <<'PY'
import os, sys, glob, re
from datetime import datetime
pool, run, t0, t1 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
t0 = datetime.strptime(t0, "%Y-%m-%d %H:%M:%S"); t1 = datetime.strptime(t1, "%Y-%m-%d %H:%M:%S")
def stamp(p):
    m = re.search(r"(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})\.bin$", p)
    return datetime(*(int(g) for g in m.groups()))
for sub in ("ERA5Data_3D", "ERA5Data_Surface"):
    n = 0
    for src in sorted(glob.glob(os.path.join(pool, sub, "*.bin"))):
        if t0 <= stamp(src) <= t1:
            os.symlink(src, os.path.join(run, sub, os.path.basename(src))); n += 1
    if n == 0:
        sys.exit(f"FATAL: no {sub} frames in pool for {t0}..{t1}")
    print(f"  linked {n} {sub} frames")
PY

# stage_run.sh copies the deck, applies the segment window + soil anchor, and
# runs the preflight (max_step, frame coverage, frame[0]==start, 6-field
# surface frames). It skips the frame copy when the run dir already holds them.
docker run --rm --gpus all -v "$ERF:/app/ERF" -w /app/ERF erf-hindcast \
    Exec/CanonicalTests/ChannelIslands/stage_run.sh \
    "USE_EXISTING" "/app/ERF/run_$SEG" "$START" "$STOP" "$SOIL"

[ "$STAGE_ONLY" = "--stage-only" ] && { echo "staged only: $RUN"; exit 0; }

echo "launching $SEG"
docker run --rm --gpus all -v "$ERF:/app/ERF" -w "/app/ERF/run_$SEG" erf-hindcast \
    bash -lc "/app/ERF/build/Exec/erf_exec inputs_hindcast amrex.max_gpu_streams=1 > run.log 2>&1; echo EXIT=\$? >> run.log"
tail -1 "$RUN/run.log"
