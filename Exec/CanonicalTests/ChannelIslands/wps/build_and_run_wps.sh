#!/bin/bash
# Build WRF (real.exe only) + WPS inside the erf-hindcast container and produce
# wrfinput_d01 for the ERF Noah-MP driver.
#
#   docker run --rm -v ~/ERF:/app/ERF -w /app/ERF erf-hindcast \
#       Exec/CanonicalTests/ChannelIslands/wps/build_and_run_wps.sh [stage]
#
# Stages: build | geogrid | ungrib | metgrid | real | all (default: build)
#
# Prerequisites (user-staged, on the host):
#   ~/ERF/wps_geog/          <- extracted WPS_GEOG subset (see README_geog.md)
#   ~/ERF/era5_grib/         <- the era5_3d_*.grib + era5_surf_*.grib files
#                               already downloaded by the erftools step
#   namelist.wps + namelist.input.real in this directory, dates filled in.
set -euo pipefail

WORK=/app/ERF/wps_work
SRC_WRF=$WORK/WRF
SRC_WPS=$WORK/WPS
HERE=/app/ERF/Exec/CanonicalTests/ChannelIslands/wps
STAGE="${1:-build}"

export NETCDF=/usr/local
export HDF5=/usr/local
export J="-j $(nproc)"

build_all () {
    mkdir -p $WORK && cd $WORK
    if [ ! -d $SRC_WRF ]; then
        curl -fsSL -o wrf.tar.gz \
          https://github.com/wrf-model/WRF/releases/download/v4.6.1/v4.6.1.tar.gz
        mkdir -p $SRC_WRF && tar xf wrf.tar.gz -C $SRC_WRF --strip-components=1
    fi
    if [ ! -d $SRC_WPS ]; then
        curl -fsSL -o wps.tar.gz \
          https://github.com/wrf-model/WPS/archive/refs/tags/v4.6.0.tar.gz
        mkdir -p $SRC_WPS && tar xf wps.tar.gz -C $SRC_WPS --strip-components=1
    fi
    # WRF: serial GNU (configure option 32), em_real -- we only need real.exe
    cd $SRC_WRF
    if [ ! -x main/real.exe ]; then
        printf '32\n1\n' | ./configure
        ./compile em_real 2>&1 | tail -5
        [ -x main/real.exe ] || { echo "FATAL: real.exe did not build"; exit 1; }
    fi
    # WPS: serial GNU (configure option 1)
    cd $SRC_WPS
    if [ ! -x geogrid/src/geogrid.exe ]; then
        printf '1\n' | ./configure
        # WPS links against the WRF build above
        sed -i "s|^WRF_DIR.*|WRF_DIR = $SRC_WRF|" configure.wps
        ./compile 2>&1 | tail -5
        [ -x geogrid/src/geogrid.exe ] || { echo "FATAL: geogrid.exe did not build"; exit 1; }
        [ -x ungrib/src/ungrib.exe ]  || { echo "FATAL: ungrib.exe did not build"; exit 1; }
        [ -x metgrid/src/metgrid.exe ] || { echo "FATAL: metgrid.exe did not build"; exit 1; }
    fi
    echo "WRF + WPS built."
}

run_geogrid () {
    cd $SRC_WPS
    cp $HERE/namelist.wps namelist.wps
    [ -d /app/ERF/wps_geog ] || { echo "FATAL: /app/ERF/wps_geog not staged"; exit 1; }
    ./geogrid/src/geogrid.exe
    ls -la geo_em.d01.nc
}

run_ungrib () {
    cd $SRC_WPS
    [ -d /app/ERF/era5_grib ] || { echo "FATAL: /app/ERF/era5_grib not staged"; exit 1; }
    # ERA5 pressure-level + surface GRIB: the ERA-interim pressure-level Vtable
    # matches ERA5 field codes
    ln -sf ungrib/Variable_Tables/Vtable.ERA-interim.pl Vtable
    ./link_grib.csh /app/ERF/era5_grib/era5_*
    ./ungrib/src/ungrib.exe
}

run_metgrid () {
    cd $SRC_WPS
    ./metgrid/src/metgrid.exe
    ls -la met_em.d01.*.nc
    echo "NOTE: set num_metgrid_levels in namelist.input.real to the value"
    echo "      shown by: ncdump -h met_em.d01.*.nc | grep num_metgrid_levels"
}

run_real () {
    cd $SRC_WRF/test/em_real
    cp $HERE/namelist.input.real namelist.input
    ln -sf $SRC_WPS/met_em.d01.*.nc .
    ./real.exe || { tail -20 rsl.error.0000 2>/dev/null; exit 1; }
    [ -f wrfinput_d01 ] || { echo "FATAL: wrfinput_d01 not produced"; exit 1; }
    cp wrfinput_d01 /app/ERF/Exec/CanonicalTests/ChannelIslands/wps/wrfinput_d01
    echo "wrfinput_d01 written to Exec/CanonicalTests/ChannelIslands/wps/"
}

case "$STAGE" in
    build)   build_all ;;
    geogrid) run_geogrid ;;
    ungrib)  run_ungrib ;;
    metgrid) run_metgrid ;;
    real)    run_real ;;
    all)     build_all; run_geogrid; run_ungrib; run_metgrid; run_real ;;
    *) echo "unknown stage: $STAGE"; exit 1 ;;
esac
