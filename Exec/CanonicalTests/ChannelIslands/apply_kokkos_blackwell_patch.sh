#!/bin/bash
# *** NECESSARY BUT NOT SUFFICIENT -- THE BLACKWELL PORT DOES NOT WORK. ***
# Kept as documentation. This patch fixes only blocker #1 of three; see the
# header of Dockerfile.blackwell for the full diagnosis. Applying this alone
# gets you a clean CONFIGURE and then a compile failure in
# Kokkos_NvidiaGpuArchitectures.hpp. Port dropped 2026-07-30 in favour of sm_89.
# Teach the bundled Kokkos 4.5 about NVIDIA Blackwell CC 12.0 (sm_120).
#
# WHY THIS EXISTS
#   Kokkos 4.5's CUDA architecture list ends at HOPPER90/sm_90
#   (cmake/kokkos_arch.cmake:836). There is NO Blackwell option. With no
#   Kokkos_ARCH_* set, Kokkos auto-detects from the BUILD MACHINE's GPU
#   (kokkos_arch.cmake:1071) -- so building on an sm_89 box for an sm_120
#   target silently yields a Kokkos compiled for the wrong architecture while
#   AMReX and ERF are compiled for sm_120. That mixed binary is the failure
#   mode this patch prevents.
#
#   Upstream Kokkos added Blackwell in a later release. Bumping the submodule
#   is the clean fix, but ekat pins Kokkos and the ekat/RRTMGP combination at a
#   newer Kokkos is untested here. This patch is the minimal alternative: it
#   adds the arch OPTION so the correct -arch=sm_120 is emitted. It changes no
#   kernel code and no codegen logic.
#
# WHAT IS AND IS NOT VERIFIED
#   Verified: the patched tree CONFIGURES and COMPILES for sm_120.
#   NOT verified: correct EXECUTION on Blackwell hardware. No sm_120 device was
#   available when this was written. Kokkos 4.5 has never been validated on
#   Blackwell by anyone. Treat the first real run on a 5090 as unvalidated and
#   check results against a known-good sm_89 run before trusting them.
#
# Byte-exact hash gating, same contract as apply_rrtmgp_patch.sh:
#   already patched -> OK (no-op)
#   pristine        -> patch, OK
#   anything else   -> HARD FAIL (submodule moved; re-derive deliberately)
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
K="$ERF_ROOT/Submodules/ekat/extern/kokkos/cmake"
ARCH="$K/kokkos_arch.cmake"
CONF="$K/KokkosCore_config.h.in"

PRISTINE_ARCH=c0900932a49bddf8809df80dddb51d5634904837e19549f374ba1172254334c8
PRISTINE_CONF=0e3e26cba1144b59ff10559eaf7687b0c317af28402c1a6f97db146ffd64e0d1

fail () { echo "FATAL [kokkos-blackwell patch]: $*" >&2; exit 1; }

[ -f "$ARCH" ] || fail "not found: $ARCH (submodule not initialized?)"
[ -f "$CONF" ] || fail "not found: $CONF"

a_done=0; c_done=0
grep -q 'BLACKWELL120' "$ARCH" && a_done=1
grep -q 'KOKKOS_ARCH_BLACKWELL120' "$CONF" && c_done=1

if [ "$a_done" = 1 ] && [ "$c_done" = 1 ]; then
    echo "[kokkos-blackwell patch] already applied."
    exit 0
fi
if [ "$a_done" != "$c_done" ]; then
    fail "half-applied state (arch=$a_done config=$c_done) -- revert the kokkos submodule and re-run"
fi

cur_a=$(sha256sum "$ARCH" | cut -d' ' -f1)
cur_c=$(sha256sum "$CONF" | cut -d' ' -f1)
[ "$cur_a" = "$PRISTINE_ARCH" ] || fail "kokkos_arch.cmake hash $cur_a matches neither pristine nor patched -- Kokkos moved; re-derive this patch"
[ "$cur_c" = "$PRISTINE_CONF" ] || fail "KokkosCore_config.h.in hash $cur_c matches neither pristine nor patched -- Kokkos moved; re-derive this patch"

# 1. declare the arch option, immediately after HOPPER90
sed -i 's|^\(kokkos_arch_option(HOPPER90 GPU .*\)$|\1\nkokkos_arch_option(BLACKWELL120 GPU "NVIDIA Blackwell generation CC 12.0" "KOKKOS_SHOW_CUDA_ARCHS")|' "$ARCH" \
    || fail "sed failed on the arch option"

# 2. map it to sm_120, immediately after the HOPPER90 mapping
sed -i 's|^\(check_cuda_arch(HOPPER90 sm_90)\)$|\1\ncheck_cuda_arch(BLACKWELL120 sm_120)|' "$ARCH" \
    || fail "sed failed on check_cuda_arch"

# 3. expose the define so KokkosCore_config.h carries it
sed -i 's|^\(#cmakedefine KOKKOS_ARCH_HOPPER90\)$|\1\n#cmakedefine KOKKOS_ARCH_BLACKWELL120|' "$CONF" \
    || fail "sed failed on the config define"

grep -q 'kokkos_arch_option(BLACKWELL120' "$ARCH" || fail "verification failed: arch option not present after patch"
grep -q 'check_cuda_arch(BLACKWELL120 sm_120)' "$ARCH" || fail "verification failed: sm_120 mapping not present after patch"
grep -q 'KOKKOS_ARCH_BLACKWELL120' "$CONF" || fail "verification failed: config define not present after patch"

echo "[kokkos-blackwell patch] applied and verified (3 insertions)."
