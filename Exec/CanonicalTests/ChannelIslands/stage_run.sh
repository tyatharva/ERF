#!/bin/bash
# Assemble the ChannelIslands hindcast run directory from all pipeline outputs.
# Run INSIDE the container:
#   docker run --rm --gpus all -v ~/ERF:/app/ERF -w /app/ERF erf-hindcast \
#       Exec/CanonicalTests/ChannelIslands/stage_run.sh
# Hard-fails (with a named reason) on any missing ingredient.
set -euo pipefail
CI=/app/ERF/Exec/CanonicalTests/ChannelIslands
RUN=/app/ERF/run_hindcast
ERA5_OUT=/app/ERF/era5_run/Output

fail () { echo "FATAL [stage_run]: $*" >&2; exit 1; }

[ -d "$ERA5_OUT/ERA5Data_3D" ]      || fail "$ERA5_OUT/ERA5Data_3D missing -- run the erftools ERA5 step first"
[ -d "$ERA5_OUT/ERA5Data_Surface" ] || fail "$ERA5_OUT/ERA5Data_Surface missing -- run the erftools ERA5 step first"
[ -f "$CI/channel_islands_terrain.txt" ] || fail "terrain file missing -- run dem_to_erf_terrain.py"
[ -x /app/ERF/build/Exec/erf_exec ] || fail "erf_exec not built -- run Build/cmake_single_precision_cuda.sh"

mkdir -p $RUN
cp    $CI/inputs_hindcast                          $RUN/
cp    $CI/channel_islands_terrain.txt              $RUN/
cp    /app/ERF/Submodules/RRTMGP/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc  $RUN/
cp    /app/ERF/Submodules/RRTMGP/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc  $RUN/
cp    /app/ERF/Submodules/RRTMGP/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc $RUN/
cp    /app/ERF/Submodules/RRTMGP/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc $RUN/
rm -rf $RUN/ERA5Data_3D $RUN/ERA5Data_Surface
cp -r $ERA5_OUT/ERA5Data_3D      $RUN/ERA5Data_3D
cp -r $ERA5_OUT/ERA5Data_Surface $RUN/ERA5Data_Surface

n3=$(ls $RUN/ERA5Data_3D/*.bin 2>/dev/null | wc -l)
ns=$(ls $RUN/ERA5Data_Surface/*.bin 2>/dev/null | wc -l)
echo "staged: $n3 3D frames, $ns surface frames"
[ "$n3" -ge 9 ] || fail "expected >= 9 three-hourly 3D frames for 24 h, found $n3"
echo "Run directory ready: $RUN"
echo "Launch:  cd $RUN && /app/ERF/build/Exec/erf_exec inputs_hindcast"
