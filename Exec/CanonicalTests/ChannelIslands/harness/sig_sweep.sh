#!/bin/bash
# Sweep the LODI relaxation coefficient on the Riemann outflow branch.
#
# sigma = 1 is the current Riemann arm (interior 2.63x WET over the day);
# sigma = 0 reduces algebraically to the extrapolation arm (0.19x DRY). Both
# endpoints are already measured on a full day, so this sweep is an interpolation
# between two known points, not an extrapolation into the unknown.
# Lower sigma = more non-reflecting = drains more freely.
#
# Every arm stops at the same MODEL time (06:00), not the same step count, so the
# interior amounts are comparable.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh

for SIG in "$@"; do
  TAG=$(echo "$SIG" | tr -d '.')
  CRUN=/app/ERF/run_sig$TAG; HRUN=/home/atyagi/ERF/run_sig$TAG
  docker run --rm -v /home/atyagi/ERF:/app/ERF erf-hindcast bash -c "
    set -e
    rm -rf $CRUN && mkdir -p $CRUN && cd $CRUN
    cp /app/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast .
    cp /app/ERF/run_a3/channel_islands_terrain_3km_192x96.txt .
    cp /app/ERF/run_a3/rrtmgp-*.nc /app/ERF/run_a3/sfc_anchor_jan09_192x96.bin .
    cp -r /app/ERF/run_a3/ERA5Data_3D /app/ERF/run_a3/ERA5Data_Surface .
    sed 's|^stop_datetime.*|stop_datetime = \"2023-01-09 06:00:00\"|' inputs_hindcast > inputs_6h"
  grep -q '2023-01-09 06:00:00' "$HRUN/inputs_6h" || { echo "PRECOND sigma=$SIG: stop not set"; exit 4; }

  "$PF" || { echo "PREFLIGHT FAIL"; exit 3; }
  T0=$(date +%s)
  docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w $CRUN erf-hindcast \
    /app/ERF/build/Exec/erf_exec inputs_6h \
      erf.nscbc_lateral=1 erf.nscbc_outflow=1 erf.nscbc_parts=31 erf.nscbc_mass_tau=60 \
      erf.hindcast_mass_du_max=8 erf.nscbc_sigma=$SIG \
      amr.check_per=-1 amr.check_int=-1 \
    > "$SP/sig$TAG.full" 2>&1
  EX=$?; T1=$(date +%s)
  M0=$(grep -a ' MASS       =' "$SP/sig$TAG.full" | head -1 | awk '{print $3}')
  MN=$(grep -a ' MASS       =' "$SP/sig$TAG.full" | tail -1 | awk '{print $3}')
  echo "sigma=$SIG exit=$EX wall=$(( (T1-T0)/60 ))min last=[$(grep -a 'Coarse STEP .* ends' "$SP/sig$TAG.full" | tail -1)]"
  echo "   mass $(awk -v a=$M0 -v b=$MN 'BEGIN{printf "%+.3f%%", 100*(b-a)/a}')  NaN=$(grep -ac 'contains NaNs' "$SP/sig$TAG.full")  FPE=$(grep -ac Erroneous "$SP/sig$TAG.full")"
done
echo SIGSWEEP_DONE
