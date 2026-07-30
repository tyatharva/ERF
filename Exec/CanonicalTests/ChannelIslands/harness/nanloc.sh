#!/bin/bash
# 30p NaN-localization gate.
#
# Sanitizer verdict on the 4-stream reproducer: ERROR SUMMARY 0 errors, but ALL
# 15 conserved components fully NaN after advance_dycore at step 1.  So the
# fill_from_realbdy SIGFPE is not a memory fault -- it is the first HOST ordered
# compare (wall-flux / mass-controller clamp) touching an already-NaN reduction.
#
# This run removes the trap and turns the per-step NaN reporter on, so the
# 18.0003 h event either (a) names the first non-finite component + cells, or
# (b) does not happen at all -- in which case the trap was spurious and this IS
# the continuous 24 h.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
DECK=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast

"$PF" || { echo "PREFLIGHT FAIL"; exit 3; }

# All run-dir mutation in-container (dir is root-owned); then VERIFY by content.
docker run --rm -v /home/atyagi/ERF:/app/ERF erf-hindcast bash -c '
  set -e
  rm -rf /app/ERF/run_nan && mkdir -p /app/ERF/run_nan
  cd /app/ERF/run_nan
  cp /app/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast .
  cp /app/ERF/run_a3/channel_islands_terrain_3km_192x96.txt .
  cp /app/ERF/run_a3/rrtmgp-*.nc /app/ERF/run_a3/sfc_anchor_jan09_192x96.bin .
  cp -r /app/ERF/run_a3/ERA5Data_3D /app/ERF/run_a3/ERA5Data_Surface .
'
# Preconditions, checked on CONTENT not on exit codes.
grep -q 'erf.hindcast_global_mass_tau *= *60'   /home/atyagi/ERF/run_nan/inputs_hindcast || { echo "PRECOND: mass tau missing"; exit 4; }
grep -q 'stop_datetime *= *"2023-01-10 00:00:00"' /home/atyagi/ERF/run_nan/inputs_hindcast || { echo "PRECOND: stop_datetime wrong"; exit 4; }
grep -q 'max_step *= *-1'                       /home/atyagi/ERF/run_nan/inputs_hindcast || { echo "PRECOND: max_step wrong"; exit 4; }
diff -q "$DECK" /home/atyagi/ERF/run_nan/inputs_hindcast || { echo "PRECOND: deck copy differs"; exit 4; }
[ "$(ls /home/atyagi/ERF/run_nan/ERA5Data_3D/*.bin | wc -l)" = 9 ] || { echo "PRECOND: wrong frame count"; exit 4; }
# every file the deck names must exist in the run dir (the first take died on a
# missing RRTMGP coefficient set -- check the artifact, not the exit code)
for f in rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc \
         rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc \
         sfc_anchor_jan09_192x96.bin channel_islands_terrain_3km_192x96.txt; do
  [ -s "/home/atyagi/ERF/run_nan/$f" ] || { echo "PRECOND: missing $f"; exit 4; }
done
echo "PRECOND OK"

T0=$(date +%s)
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_nan erf-hindcast \
  /app/ERF/build/Exec/erf_exec inputs_hindcast \
    amrex.fpe_trap_invalid=0 \
    erf.check_for_nans=1 erf.check_for_nans_int=10 \
    amr.plot_int=100000 erf.plot_per_1=-1 \
  > "$SP/nanloc.full" 2>&1
EX=$?
T1=$(date +%s)

echo "NANLOC exit=$EX wall=$(( (T1-T0)/60 ))min"
echo "last: $(grep -a 'Coarse STEP .* ends' "$SP/nanloc.full" | tail -1)"
echo "nan reports: $(grep -ac 'contains NaNs/Infs' "$SP/nanloc.full")"
grep -a -m3 -A3 'contains NaNs/Infs' "$SP/nanloc.full" | head -20
grep -a 'Erroneous\|SIGABRT\|SIGFPE\|illegal' "$SP/nanloc.full" | head -5
echo NANLOC_DONE
