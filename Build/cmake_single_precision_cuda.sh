#!/bin/bash
# Single-precision CUDA build for ERF.
# Based on Build/cmake_with_kokkos_many_cuda.sh with:
#   - ERF_PRECISION=SINGLE               (build in single precision -> AMReX defines AMREX_USE_FLOAT)
#   - Ada (RTX 4080, compute 8.9) CUDA arch flags for AMReX + Kokkos
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
# DIAGNOSTIC (UPSTREAM_ISSUES 29e): exposes AMReX flushFBCache/flushCPCache so
# erf.realbdy_flush_fb_cache=1 can test the BDKey-collision hypothesis. Inert
# with the knob off (default), but the source references the symbols, so the
# patch must be applied for the build to compile.
"$ERF_SOURCE_DIR/Exec/CanonicalTests/ChannelIslands/apply_amrex_fbcache_patch.sh" || exit 1

echo "Source: $ERF_SOURCE_DIR | Build: $ERF_BUILD_DIR | Install: $ERF_INSTALL_DIR | PWD: $(pwd)"
echo "Customize: export ERF_BUILD_DIR=... ERF_SOURCE_DIR=... ERF_INSTALL_DIR=... or ERF_HOME=..."

cmake -DCMAKE_INSTALL_PREFIX:PATH=$ERF_INSTALL_DIR \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_PRECISION:STRING=SINGLE \
      -DCMAKE_CUDA_ARCHITECTURES=89 \
      -DAMReX_CUDA_ARCH=8.9 \
      -DKokkos_ARCH_ADA89:BOOL=ON \
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
