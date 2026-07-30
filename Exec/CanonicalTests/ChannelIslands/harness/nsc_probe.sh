#!/bin/bash
# NSCBC with the 19d regime fix: does it clear the 4.7 h dt-collapse death?
# Gate take 2 (same deck, same knobs but nscbc_outflow=0 and the driver-driven
# sigma) died at TIME = 17073 s. Run past it with margin.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
MS=${1:-16000}
"$PF" || { echo "PREFLIGHT FAIL"; exit 3; }

docker run --rm -v /home/atyagi/ERF:/app/ERF erf-hindcast bash -c '
  set -e
  rm -rf /app/ERF/run_nsc && mkdir -p /app/ERF/run_nsc && cd /app/ERF/run_nsc
  cp /app/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast .
  cp /app/ERF/run_a3/channel_islands_terrain_3km_192x96.txt .
  cp /app/ERF/run_a3/rrtmgp-*.nc /app/ERF/run_a3/sfc_anchor_jan09_192x96.bin .
  cp -r /app/ERF/run_a3/ERA5Data_3D /app/ERF/run_a3/ERA5Data_Surface .
'
for f in rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc \
         rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc \
         sfc_anchor_jan09_192x96.bin channel_islands_terrain_3km_192x96.txt inputs_hindcast; do
  [ -s "/home/atyagi/ERF/run_nsc/$f" ] || { echo "PRECOND: missing $f"; exit 4; }
done
echo "PRECOND OK"

T0=$(date +%s)
docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_nsc erf-hindcast \
  /app/ERF/build/Exec/erf_exec inputs_hindcast max_step="$MS" \
    erf.nscbc_lateral=1 erf.nscbc_outflow=1 erf.nscbc_parts=31 erf.nscbc_mass_tau=60 \
    amr.check_per=-1 erf.sum_interval=200 \
  > "$SP/nsc_probe.full" 2>&1
EX=$?
T1=$(date +%s)

echo "NSCPROBE exit=$EX wall=$(( (T1-T0)/60 ))min"
echo "last: $(grep -a 'Coarse STEP .* ends' "$SP/nsc_probe.full" | tail -1)"
echo "control: gate take 2 (nscbc_outflow=0, driver-driven sigma) died at TIME = 17073"
echo "mass first: $(grep -a ' MASS       =' "$SP/nsc_probe.full" | head -1)"
echo "mass last : $(grep -a ' MASS       =' "$SP/nsc_probe.full" | tail -1)"
echo "min dt    : $(grep -a 'Coarse STEP .* ends' "$SP/nsc_probe.full" | awk '{print $10}' | sort -g | head -1)"
grep -a 'Erroneous\|SIGABRT\|contains NaNs\|illegal' "$SP/nsc_probe.full" | head -5
echo NSCPROBE_DONE
