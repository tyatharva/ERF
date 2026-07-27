#!/bin/bash
# DIAGNOSTIC ONLY (UPSTREAM_ISSUES 29e). Expose AMReX's private
# flushFBCache()/flushCPCache() so an instrumented build can drop the
# FillBoundary/ParallelCopy metadata caches after hindcast boundary-plane setup.
#
# Tests one thing: whether a recycled BoxArray/DistributionMapping id from the
# ~250 create/destroy cycles in strip_to_global_fab has seeded a stale cache
# entry that a later FillBoundary matches by key. Clears the item-29 fault ->
# BDKey collision confirmed. Does not clear -> BDKey is dead as a hypothesis.
#
# NOT part of any fix. Revert before shipping:
#   git -C Submodules/AMReX checkout Src/Base/AMReX_FabArrayBase.H
#
# Same byte-exact gating as apply_rrtmgp_patch.sh: the target must hash to one
# of exactly two known states, or this hard-fails rather than guessing.
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/AMReX/Src/Base/AMReX_FabArrayBase.H"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/amrex-expose-fbcache-flush.patch"

UNPATCHED_SHA=c44f5e9c74b18e21ccfc3efcd32b8ddc08238053cd6680444f9c9618cd087abd
PATCHED_SHA=3ad22f21e4299621b50b921250dc6f300ed370d302f20f9b7823448781e2d287

fail () { echo "FATAL [amrex-fbcache patch]: $*" >&2; exit 1; }
[ -f "$TARGET" ] || fail "target not found: $TARGET"
[ -f "$PATCH" ]  || fail "patch not found: $PATCH"

cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
if [ "$cur" = "$PATCHED_SHA" ]; then echo "[amrex-fbcache patch] already applied (hash verified)."; exit 0; fi
[ "$cur" = "$UNPATCHED_SHA" ] || fail "target hashes to $cur, which is neither the
  known pristine nor the known patched state. The AMReX submodule has moved;
  re-derive this patch deliberately rather than applying it blind."

patch -p1 -d "$ERF_ROOT/Submodules/AMReX" < "$PATCH" || fail "patch application failed"
new=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$new" = "$PATCHED_SHA" ] || fail "post-patch hash $new != expected $PATCHED_SHA"
echo "[amrex-fbcache patch] applied and hash-verified."
