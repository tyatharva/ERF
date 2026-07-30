#!/bin/bash
# *** DOES NOT PRODUCE A WORKING BUILD -- KEPT AS DOCUMENTATION ONLY. ***
# The sm_120 port is blocked by three pre-Blackwell components; see the header of
# Exec/CanonicalTests/ChannelIslands/Dockerfile.blackwell. In particular AMReX
# silently ignores AMReX_CUDA_ARCH=12.0 and builds for compute_60..86, so this
# script's flags read back correct in the CMake cache and are not what nvcc gets.
# Port dropped 2026-07-30 in favour of sm_89 hardware. Use
# cmake_single_precision_cuda.sh.
# Single-precision CUDA build for ERF.
# Based on Build/cmake_with_kokkos_many_cuda.sh with:
#   - ERF_PRECISION=SINGLE               (build in single precision -> AMReX defines AMREX_USE_FLOAT)
#   - Blackwell (RTX 5090, compute 12.0) CUDA arch flags for AMReX + Kokkos
#   - FCOMPARE / TESTS disabled           (not needed for this workflow)
#   - parallel build with -j$(nproc)
#
# Defaults (customize here or set ERF_BUILD_DIR/ERF_SOURCE_DIR/ERF_INSTALL_DIR in environment)
# If ERF_HOME is set, use it as base for absolute paths
if [ -n "$ERF_HOME" ]; then
  : ${ERF_BUILD_DIR:="$ERF_HOME/build"}
  : ${ERF_SOURCE_DIR:="$ERF_HOME"}
  : ${ERF_INSTALL_DIR:="$ERF_HOME/install"}
else
  : ${ERF_BUILD_DIR:="."}
  : ${ERF_SOURCE_DIR:=".."}
  : ${ERF_INSTALL_DIR:="install"}
fi

# Single-precision RRTMGP requires the minor-gas-scaling overflow fix in the
# RRTMGP submodule. Applied with byte-exact verification; the build HARD-FAILS
# if the file is in an unknown state. See Exec/CanonicalTests/ChannelIslands/.
"$ERF_SOURCE_DIR/Exec/CanonicalTests/ChannelIslands/apply_rrtmgp_patch.sh" || exit 1

# Kokkos 4.5 has no Blackwell architecture (its CUDA list ends at HOPPER90/sm_90).
# Without an explicit arch it auto-detects from the BUILD machine's GPU, which
# would silently compile Kokkos for the wrong architecture while AMReX and ERF
# are built for sm_120. This adds the arch option; hash-gated, hard-fails if the
# Kokkos submodule has moved. See the script header for what is NOT verified.
"$ERF_SOURCE_DIR/Exec/CanonicalTests/ChannelIslands/apply_kokkos_blackwell_patch.sh" || exit 1

echo "Source: $ERF_SOURCE_DIR | Build: $ERF_BUILD_DIR | Install: $ERF_INSTALL_DIR | PWD: $(pwd)"
echo "Customize: export ERF_BUILD_DIR=... ERF_SOURCE_DIR=... ERF_INSTALL_DIR=... or ERF_HOME=..."

cmake -DCMAKE_INSTALL_PREFIX:PATH=$ERF_INSTALL_DIR \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_PRECISION:STRING=SINGLE \
      -DCMAKE_CUDA_ARCHITECTURES=120 \
      -DAMReX_CUDA_ARCH=12.0 \
      -DKokkos_ARCH_BLACKWELL120:BOOL=ON \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=OFF \
      -DERF_ENABLE_FCOMPARE:BOOL=OFF \
      -DERF_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_HDF5:BOOL=ON \
      -DERF_ENABLE_FFT:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      --log-context --log-level STATUS \
      -B $ERF_BUILD_DIR -S $ERF_SOURCE_DIR && \
cmake --build $ERF_BUILD_DIR -j$(nproc) && \
cmake --install $ERF_BUILD_DIR --prefix=$ERF_INSTALL_DIR
