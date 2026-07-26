#!/bin/bash
# Assemble the ChannelIslands hindcast run directory from all pipeline outputs.
# Run INSIDE the container:
#   docker run --rm --gpus all -v ~/ERF:/app/ERF -w /app/ERF erf-hindcast \
#       Exec/CanonicalTests/ChannelIslands/stage_run.sh
# Hard-fails (with a named reason) on any missing ingredient.
set -euo pipefail
# Defaults reproduce the original single-run behaviour; a month segment
# passes its own frame source and run directory:
#   stage_run.sh /app/ERF/era5_year_2023/2023-08/Output /app/ERF/run_2023-08 \
#       "2023-07-30 00:00:00" "2023-09-01 00:00:00" 298.9
# The datetimes and soil anchor are applied to the staged deck BEFORE the
# preflight runs -- the frame-coverage assertions are only meaningful
# against the segment's real window.
CI=/app/ERF/Exec/CanonicalTests/ChannelIslands
ERA5_OUT=${1:-/app/ERF/era5_run/Output}
RUN=${2:-/app/ERF/run_hindcast}
SEG_START=${3:-}
SEG_STOP=${4:-}
SEG_SOIL=${5:-}

fail () { echo "FATAL [stage_run]: $*" >&2; exit 1; }

# USE_EXISTING: the caller already populated $RUN/ERA5Data_* (launch_segment.sh
# symlinks a slice of the flat year pool, because ERF indexes frames
# positionally and each segment must start at its own frame[0]).
if [ "$ERA5_OUT" = USE_EXISTING ]; then
    [ -d "$RUN/ERA5Data_3D" ] && [ -d "$RUN/ERA5Data_Surface" ] \
        || fail "USE_EXISTING but $RUN/ERA5Data_* not populated"
else
    [ -d "$ERA5_OUT/ERA5Data_3D" ]      || fail "$ERA5_OUT/ERA5Data_3D missing -- run the erftools ERA5 step first"
    [ -d "$ERA5_OUT/ERA5Data_Surface" ] || fail "$ERA5_OUT/ERA5Data_Surface missing -- run the erftools ERA5 step first"
fi
[ -f "$CI/channel_islands_terrain.txt" ] || fail "terrain file missing -- run dem_to_erf_terrain.py"
[ -x /app/ERF/build/Exec/erf_exec ] || fail "erf_exec not built -- run Build/cmake_single_precision_cuda.sh"

mkdir -p $RUN
cp    $CI/inputs_hindcast                          $RUN/
if [ -n "$SEG_START" ]; then
    sed -i "s|^start_datetime = .*|start_datetime = \"$SEG_START\"|;    \
            s|^stop_datetime  = .*|stop_datetime  = \"$SEG_STOP\"|;     \
            s|^erf.mm5.soil_theta     = .*|erf.mm5.soil_theta     = $SEG_SOIL|" \
        $RUN/inputs_hindcast
    echo "segment deck: $SEG_START -> $SEG_STOP, soil_theta=$SEG_SOIL K"
fi
# Stage every terrain variant present (the deck names which one it uses;
# previously only the 2-km file was staged while the deck referenced the
# 3-km one).
cp    $CI/channel_islands_terrain*.txt             $RUN/
# ERA5 surface anchor for the IC hydrostatic integration (audit item 4 SEV-1).
# Case-specific: regenerate with precip_check/make_sfc_anchor.py per date/domain.
ANCHOR=$(sed -n 's/^erf.hindcast_sfc_anchor_file[ =]*"\([^"]*\)".*/\1/p' $CI/inputs_hindcast)
[ -n "$ANCHOR" ] || fail "erf.hindcast_sfc_anchor_file not set in the deck"
[ -f "/app/ERF/precip_check/$ANCHOR" ] \
    || fail "anchor $ANCHOR missing -- run precip_check/make_sfc_anchor.py"
cp    /app/ERF/precip_check/$ANCHOR                 $RUN/
cp    /app/ERF/Submodules/RRTMGP/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc  $RUN/
cp    /app/ERF/Submodules/RRTMGP/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc  $RUN/
cp    /app/ERF/Submodules/RRTMGP/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc $RUN/
cp    /app/ERF/Submodules/RRTMGP/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc $RUN/
if [ "$ERA5_OUT" != USE_EXISTING ]; then
    rm -rf $RUN/ERA5Data_3D $RUN/ERA5Data_Surface
    cp -r $ERA5_OUT/ERA5Data_3D      $RUN/ERA5Data_3D
    cp -r $ERA5_OUT/ERA5Data_Surface $RUN/ERA5Data_Surface
fi

n3=$(ls $RUN/ERA5Data_3D/*.bin 2>/dev/null | wc -l)
ns=$(ls $RUN/ERA5Data_Surface/*.bin 2>/dev/null | wc -l)
echo "staged: $n3 3D frames, $ns surface frames"

# ---------------------------------------------------------------------------
# PRE-LAUNCH PREFLIGHT (audit 2026-07-24). Everything here fails LOUDLY at
# stage time; each check guards a measured silent-failure mode:
#   1. max_step must be -1: any positive cap ends the run early with a
#      clean exit (ERF stops at min(max_step, stop_datetime); istep
#      continues across restarts).
#   2. ERA5 frame coverage must span start_datetime..stop_datetime
#      INCLUSIVE (N*8+1 frames; the stop-instant 00Z frame is consumed by
#      the final interpolation window), and frame[0] must equal
#      start_datetime (ERF indexes frames as sorted-list position from t=0).
#   3. Every surface frame must carry 6 fields (field 6 = forecast
#      albedo); a 5-field frame would silently leave radiation on a stale
#      or constant albedo.
# ---------------------------------------------------------------------------
python3 - "$RUN" <<'PYEOF' || fail "preflight failed (see messages above)"
import glob, os, re, struct, sys
from datetime import datetime, timedelta
run = sys.argv[1]

deck = open(os.path.join(run, "inputs_hindcast")).read()
def deck_val(key):
    m = re.search(rf"^\s*{key}\s*=\s*\"?([^\"\n#]+)", deck, re.M)
    return m.group(1).strip() if m else None

# 1. max_step
ms = deck_val("max_step")
if ms != "-1":
    sys.exit(f"PREFLIGHT: max_step = {ms}; must be -1 for production "
             "(positive caps truncate the run with a clean exit)")

# 1b. the five IC/boundary corrections must be present and enabled.
#     These were validated in bdyfix/ but lived only on run-script command
#     lines for weeks, so the shipped deck silently produced the original
#     broken IC (UPSTREAM_ISSUES 19-24, audit item 4 SEV-1). Fail loudly
#     rather than let a deck without them look normal.
#     Asserted against validated_config.txt, the single source of truth.
manifest = os.path.join(os.path.dirname(os.path.abspath(run.rstrip("/"))),
                        "Exec/CanonicalTests/ChannelIslands/validated_config.txt")
if not os.path.exists(manifest):
    manifest = "/app/ERF/Exec/CanonicalTests/ChannelIslands/validated_config.txt"
if not os.path.exists(manifest):
    sys.exit(f"PREFLIGHT: validated_config.txt not found; cannot verify the deck")
bad = []
nmust = 0
for line in open(manifest):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    parts = line.split(None, 2)
    if parts[0] == "MUST":
        nmust += 1
        key, want = parts[1], parts[2].strip().strip('"')
        got = deck_val(key)
        got = got.strip().strip('"') if got is not None else None
        if got is None:
            bad.append(f"  {key}: ABSENT, must be {want}")
        elif got != want:
            bad.append(f"  {key}: deck has {got}, must be {want}")
if bad:
    sys.exit("PREFLIGHT: the staged deck diverges from the validated "
             "configuration (validated_config.txt):\n" + "\n".join(bad) +
             "\n  Either fix the deck or move the knob to an EXCEPT line "
             "with a reason. Eight of these had drifted silently before the "
             "2026-07-26 audit.")
print(f"PREFLIGHT: {nmust} validated-config knobs match")
anchor = deck_val("erf.hindcast_sfc_anchor_file").strip('"')
if not os.path.exists(os.path.join(run, anchor)):
    sys.exit(f"PREFLIGHT: anchor file {anchor} not staged into {run}. "
             "Generate it with precip_check/make_sfc_anchor.py for this "
             "date and domain.")
if deck_val("erf.rad_freq_in_steps") is not None:
    sys.exit("PREFLIGHT: erf.rad_freq_in_steps is still in the deck; it is "
             "superseded by erf.rad_freq_in_time and the two together make the "
             "radiative timescale dt-dependent. Remove it.")

# 2. frame coverage vs deck datetimes
t0 = datetime.strptime(deck_val("start_datetime"), "%Y-%m-%d %H:%M:%S")
t1 = datetime.strptime(deck_val("stop_datetime"), "%Y-%m-%d %H:%M:%S")
iv = timedelta(hours=float(deck_val("erf.hindcast_data_interval_in_hrs") or 3.0))
need = int((t1 - t0) / iv) + 1
def stamp(p):
    m = re.search(r"(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})\.bin$", p)
    return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)),
                    int(m.group(4)), int(m.group(5)))
for sub, tagname in (("ERA5Data_3D", "3D"), ("ERA5Data_Surface", "surface")):
    files = sorted(glob.glob(os.path.join(run, sub, "*.bin")))
    if len(files) < need:
        sys.exit(f"PREFLIGHT: {tagname} frames: {len(files)} < {need} needed for "
                 f"{t0}..{t1} at {iv} (the stop-instant frame is required)")
    if stamp(files[0]) != t0:
        sys.exit(f"PREFLIGHT: first {tagname} frame is {stamp(files[0])}, but "
                 f"start_datetime is {t0}: ERF indexes frames positionally "
                 f"from t=0 -- frame[0] MUST be the start instant")
    if stamp(files[need - 1]) != t1:
        sys.exit(f"PREFLIGHT: {tagname} frame[{need-1}] is {stamp(files[need-1])}, "
                 f"expected {t1}: gap or misordering inside the frame series")

# 3. surface frames all 6-field (albedo present)
bad = []
for p in sorted(glob.glob(os.path.join(run, "ERA5Data_Surface", "*.bin"))):
    with open(p, "rb") as f:
        nx, ny, nz, nd = struct.unpack("iiii", f.read(16))
    if nd != 6:
        bad.append((os.path.basename(p), nd))
if bad:
    sys.exit(f"PREFLIGHT: {len(bad)} surface frame(s) do not carry 6 fields "
             f"(albedo would silently fall back/stale): {bad[:3]} ...")

print(f"PREFLIGHT OK: max_step=-1; {need} frames span {t0}..{t1} in both "
      f"streams; all surface frames carry 6 fields (albedo present)")
PYEOF

echo "Run directory ready: $RUN"
echo "Launch:  cd $RUN && /app/ERF/build/Exec/erf_exec inputs_hindcast"
