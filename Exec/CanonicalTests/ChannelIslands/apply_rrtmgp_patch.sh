#!/bin/bash
# Apply the RRTMGP single-precision minor-gas-scaling fix to the submodule,
# with byte-exact verification -- stricter than rejecting patch fuzz/offsets:
# the target file must hash to one of exactly two known states.
#
#   already patched  -> OK (no-op)
#   pristine         -> apply patch, verify resulting hash, OK
#   anything else    -> HARD FAIL (submodule moved / diverged; re-derive the
#                       patch deliberately rather than guessing)
#
# Called from Build/cmake_single_precision_cuda.sh; safe to run standalone.
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/RRTMGP/cpp/rrtmgp/kernels/mo_gas_optics_kernels.h"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/rrtmgp-sp-minor-scaling.patch"

# Known byte-exact states of the target file:
UNPATCHED_SHA=07cd73825f110285aea50a03612ed25af845148112e76ab6c90881043e06ed55
PATCHED_SHA=078a3ea7e6001c35eac5ce289e1bb85f80ecf9ac3100971fa4bc71f1effeaddd

fail () { echo "FATAL [rrtmgp-sp patch]: $*" >&2; exit 1; }

[ -f "$TARGET" ] || fail "target not found: $TARGET (submodule not initialized?)"
[ -f "$PATCH" ]  || fail "patch file not found: $PATCH"

cur=$(sha256sum "$TARGET" | cut -d' ' -f1)

if [ "$cur" = "$PATCHED_SHA" ]; then
    echo "[rrtmgp-sp patch] already applied (hash verified)."
    exit 0
fi

if [ "$cur" != "$UNPATCHED_SHA" ]; then
    fail "mo_gas_optics_kernels.h matches NEITHER the known pristine nor the \
known patched state (sha256=$cur). The RRTMGP submodule has changed; re-derive \
the patch against the new revision instead of force-applying. Refusing to build \
single precision without the overflow fix."
fi

git -C "$ERF_ROOT/Submodules/RRTMGP" apply --whitespace=error \
    "$PATCH" || fail "git apply failed on a pristine file (should not happen)"

cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$cur" = "$PATCHED_SHA" ] || fail "post-apply hash mismatch (got $cur)"
echo "[rrtmgp-sp patch] applied and hash-verified."
