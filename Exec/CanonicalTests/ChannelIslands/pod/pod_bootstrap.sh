#!/bin/bash
# Bootstrap an ERF hindcast pod (sm_89: RTX 4090 / 4080 / L40S).
#
# The pod has NO Docker -- RunPod starts the container FROM our image, so this
# runs natively inside it. It expects the image
# ghcr.io/tyatharva/erf-hindcast:cuda126-sm89, which carries the environment
# only (no ERF source; ERF is cloned here).
#
#   bash pod_bootstrap.sh              clone + build + verify + validate
#   bash pod_bootstrap.sh --no-clone   skip the clone (repo already present)
#
# Exits non-zero and LOUDLY at the first failure. The validation gate at the end
# is a HARD GATE: if it fails, do not start production arms.
#
# NOTE ON ARCHITECTURE: sm_89 needs no patches. The Blackwell (sm_120) port is
# blocked by three independent pre-Blackwell components -- see the header of
# ../Dockerfile.blackwell. Do not try to retarget this script by changing the
# arch flags; it will appear to work and silently build for the wrong device.
set -uo pipefail

ERF_ROOT=${ERF_ROOT:-/app/ERF}
REPO=${REPO:-https://github.com/tyatharva/ERF.git}
BRANCH=${BRANCH:-ERF}
CI=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands

say  () { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
die  () { printf '\n\033[1;31mFATAL: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- 0. clone --
if [ "${1:-}" != "--no-clone" ]; then
    say "0. clone $REPO ($BRANCH)"
    [ -e "$ERF_ROOT/.git" ] && die "$ERF_ROOT already a git repo -- pass --no-clone"
    mkdir -p "$(dirname "$ERF_ROOT")"
    git clone --branch "$BRANCH" --recurse-submodules "$REPO" "$ERF_ROOT" \
        || die "clone failed"
else
    say "0. clone skipped (--no-clone)"
    [ -d "$ERF_ROOT/.git" ] || die "$ERF_ROOT is not a git checkout"
fi

git config --global --add safe.directory '*'
cd "$ERF_ROOT" || die "cannot cd $ERF_ROOT"
echo "HEAD: $(git rev-parse --short HEAD)  branch: $(git rev-parse --abbrev-ref HEAD)"

# ---------------------------------------------------- 1. submodule contract --
say "1. submodule state (validated_config.txt contract)"
git submodule update --init --recursive || die "submodule init failed"
# AMReX MUST be clean for any scored run; RRTMGP MUST carry the SP patch.
if [ -n "$(git -C Submodules/AMReX status --porcelain)" ]; then
    git -C Submodules/AMReX status --porcelain | head
    die "AMReX submodule is DIRTY. validated_config.txt requires it clean for scored runs."
fi
echo "AMReX clean: OK"

# --------------------------------------------------- 2. hash-gated patches --
say "2. RRTMGP single-precision patch (hash-gated)"
# Hard-fails if the submodule has moved; never applies blind.
"$CI/apply_rrtmgp_patch.sh" || die "RRTMGP patch failed -- do NOT force, re-derive it"

# ------------------------------------------------------------- 3. build ERF --
say "3. build ERF (single precision, CUDA, sm_89)"
t0=$(date +%s)
ERF_HOME=$ERF_ROOT "$ERF_ROOT/Build/cmake_single_precision_cuda.sh"
rc=$?
t1=$(date +%s)
BUILD_SEC=$((t1 - t0))
[ $rc -eq 0 ] || die "build failed after ${BUILD_SEC}s (exit $rc)"
BIN=$ERF_ROOT/build/Exec/erf_exec
[ -x "$BIN" ] || die "build reported success but $BIN is missing"
# Do not trust the exit code alone -- this campaign has had a build report
# exit=0 while leaving a stale binary. Check the artifact is newer than the run.
[ "$(stat -c %Y "$BIN")" -ge "$t0" ] || die "$BIN is older than this build -- STALE BINARY"
printf 'build OK in %dm %ds\n' $((BUILD_SEC / 60)) $((BUILD_SEC % 60))

# -------------------------------------------------- 4. verify configuration --
say "4. verify the build is actually sm_89 single precision"
CACHE=$ERF_ROOT/build/CMakeCache.txt
[ -f "$CACHE" ] || die "no CMakeCache.txt"

grep -q '^CMAKE_CUDA_ARCHITECTURES:[A-Z]*=89$' "$CACHE" \
    || die "CMAKE_CUDA_ARCHITECTURES is not 89: $(grep '^CMAKE_CUDA_ARCHITECTURES' "$CACHE")"
grep -q '^AMReX_CUDA_ARCH:STRING=8\.9$' "$CACHE" \
    || die "AMReX_CUDA_ARCH is not 8.9: $(grep '^AMReX_CUDA_ARCH' "$CACHE")"
grep -q '^Kokkos_ARCH_ADA89:BOOL=ON$' "$CACHE" \
    || die "Kokkos_ARCH_ADA89 is not ON: $(grep '^Kokkos_ARCH_ADA89' "$CACHE")"
echo "cache: CUDA_ARCHITECTURES=89, AMReX_CUDA_ARCH=8.9, Kokkos_ARCH_ADA89=ON"

# The cache is NOT sufficient. On the Blackwell attempt the cache read back
# correct while AMReX silently emitted compute_60..86. Check the flags AMReX
# actually recorded, which is the thing nvcc was really given.
ARCHS=$(grep '^AMREX_CUDA_ARCHS:INTERNAL=' "$CACHE" | cut -d= -f2)
echo "AMREX_CUDA_ARCHS=$ARCHS"
case ";$ARCHS;" in
    *";89;"*) echo "AMReX really built for sm_89: OK" ;;
    *) die "AMReX did NOT build for sm_89 (got '$ARCHS'). This is the silent
   fallback in CMake's deprecated FindCUDA/select_compute_arch. The binary
   would run on the wrong architecture or JIT. Do not proceed." ;;
esac

# Single precision must be real, not assumed.
grep -q '^#define AMREX_USE_FLOAT' "$ERF_ROOT/build/Submodules/AMReX/AMReX_Config_3D.H" \
    || die "AMREX_USE_FLOAT not defined -- this is a DOUBLE precision build; it will not fit in GPU memory"
echo "single precision: OK"

# --------------------------------------------------------- 5. HARD GATE -----
say "5. VALIDATION GATE -- reproduce the known-good sm_89 reference"
DATA=${DATA:-$ERF_ROOT/pod_data}
for f in CONUS404Data_3D CONUS404Data_Surface; do
    [ -d "$DATA/$f" ] || die "missing $DATA/$f -- run pod_fetch_manifest.sh / copy the transfer set first"
done
RUNV=$ERF_ROOT/run_validate
rm -rf "$RUNV"; mkdir -p "$RUNV"; cd "$RUNV" || die "cannot cd $RUNV"
ln -sf "$DATA/CONUS404Data_3D" .
ln -sf "$DATA/CONUS404Data_Surface" .
cp "$CI/channel_islands_terrain_3km_192x96.txt" .
cp "$CI/inputs_c404" . 2>/dev/null || cp "$DATA/inputs_c404" . \
    || die "no inputs_c404 available"
for t in rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/rrtmgp/data/$t" . || die "missing RRTMGP table $t"
done
for t in rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/extensions/cloud_optics/$t" . || die "missing $t"
done

# Exactly the overrides the reference was produced with. Changing ANY of these
# invalidates the comparison -- regenerate the reference instead.
"$BIN" inputs_c404 max_step=20 amr.check_int=-1 amr.plot_int=20 \
    erf.cfl=0.3 erf.moisture_model=Morrison erf.les_type=None \
    erf.hindcast_mass_du_max=8 erf.nscbc_lateral=1 erf.nscbc_outflow=1 \
    erf.nscbc_parts=31 erf.nscbc_sigma=0.03 erf.nscbc_mass_tau=60 \
    > "$RUNV/validate.log" 2>&1 \
    || { tail -20 "$RUNV/validate.log"; die "validation run did not complete"; }

[ -d "$RUNV/plt00020" ] || die "validation run produced no plt00020"
python3 "$CI/pod/pod_validate.py" "$RUNV/plt00020"
vrc=$?

echo
if [ $vrc -ne 0 ]; then
    die "VALIDATION GATE FAILED -- this build does not reproduce the sm_89 reference.
   Do NOT start production arms. Nothing computed here is trustworthy."
fi

say "BOOTSTRAP COMPLETE"
printf 'build time      : %dm %ds\n' $((BUILD_SEC / 60)) $((BUILD_SEC % 60))
printf 'binary          : %s\n' "$BIN"
printf 'validation      : PASSED\n'
cat <<'EOF'

Before any production run:
  * bash Exec/CanonicalTests/ChannelIslands/gpu_preflight.sh
  * one 23-h job per GPU -- never two on one card, and score timing only from
    solo runs (four 4090s = four independent arms, not one job across four)
  * read HANDOFF.md section 3 for the open question (item 51, sigma=1 next)
EOF
