#!/bin/bash
# NSCBC global-mass-controller gain sweep.
#
# With the 19d regime fix, NSCBC clears the 4.7 h death but mass runs +2.04% in
# 5.3 h at nscbc_mass_tau=60 -- ~9%/day, the same order as the unfixed Davies
# drift, and by the #30 rule that makes the run unscoreable. Simplest hypothesis
# first: the gain is too low, because a non-reflecting boundary partially undoes
# an imposed net wall flux (the characteristic solve re-derives u_n from the
# interior on the next call), so it needs a shorter timescale than Davies does.
# Same step count for every arm so the comparison is at matched model time.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
PF=/home/atyagi/ERF/Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
MS=3000

for TAU in 10 5; do
  "$PF" || { echo "PREFLIGHT FAIL"; exit 3; }
  docker run --rm -v /home/atyagi/ERF:/app/ERF erf-hindcast bash -c "
    rm -rf /app/ERF/run_tau && mkdir -p /app/ERF/run_tau && cd /app/ERF/run_tau
    cp /app/ERF/Exec/CanonicalTests/ChannelIslands/inputs_hindcast .
    cp /app/ERF/run_a3/channel_islands_terrain_3km_192x96.txt .
    cp /app/ERF/run_a3/rrtmgp-*.nc /app/ERF/run_a3/sfc_anchor_jan09_192x96.bin .
    cp -r /app/ERF/run_a3/ERA5Data_3D /app/ERF/run_a3/ERA5Data_Surface ."
  docker run --rm --gpus all -v /home/atyagi/ERF:/app/ERF -w /app/ERF/run_tau erf-hindcast \
    /app/ERF/build/Exec/erf_exec inputs_hindcast max_step=$MS \
      erf.nscbc_lateral=1 erf.nscbc_outflow=1 erf.nscbc_parts=31 \
      erf.nscbc_mass_tau=$TAU \
      amr.plot_int=100000 erf.plot_per_1=-1 amr.check_per=-1 erf.sum_interval=200 \
    > "$SP/nsc_tau$TAU.full" 2>&1
  EX=$?
  M0=$(grep -a ' MASS       =' "$SP/nsc_tau$TAU.full" | head -1 | awk '{print $3}')
  MN=$(grep -a ' MASS       =' "$SP/nsc_tau$TAU.full" | tail -1 | awk '{print $3}')
  MX=$(grep -a ' MASS       =' "$SP/nsc_tau$TAU.full" | awk -v m="$M0" '{d=100*($3-m)/m; if(d<0)d=-d; if(d>x)x=d} END{printf "%.3f", x}')
  LAST=$(grep -a 'Coarse STEP .* ends' "$SP/nsc_tau$TAU.full" | tail -1)
  MINDT=$(grep -a 'Coarse STEP .* ends' "$SP/nsc_tau$TAU.full" | tail -n +50 | awk '{print $10}' | sort -g | head -1)
  echo "TAU=$TAU exit=$EX  M: $M0 -> $MN  $(awk -v a="$M0" -v b="$MN" 'BEGIN{printf "%+.3f%%", 100*(b-a)/a}')  peak|dM| ${MX}%  mindt(post-transient)=$MINDT"
  echo "   last: $LAST"
done
echo NSCTAU_DONE
