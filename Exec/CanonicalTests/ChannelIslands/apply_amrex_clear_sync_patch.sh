#!/bin/bash
# DIAGNOSTIC / VERIFICATION ONLY (UPSTREAM_ISSUES 29j). DO NOT SHIP.
#
# Hoists a device-wide sync to the TOP of FabArray::clear(), before the loop that
# returns fab storage to the arena.
#
# The defect: clear() frees fab data at AMReX_FabArray.H:1952 with no stream
# ordering, and the only synchronization in the whole destructor path is
# m_fb_local_copy_handler.clear() ~25 lines later (~TagVector -> undefine() ->
# Gpu::streamSynchronize()). FillBoundary's local-copy path is fire-and-forget --
# FB_local_copy_gpu is not inside an MFIter, so it gets none of ~MFIter's
# all-stream sync -- so a short-lived MultiFab whose last operation is a
# FillBoundary hands memory back to the arena while a kernel is still writing it.
# The arena then gives that block to the CACHED FillBoundary tag vector of an
# unrelated MultiFab, and the in-flight kernel scribbles float field data over the
# fab-pointer array.
#
# This patch exists to confirm the mechanism INDEPENDENTLY of any ERF change: with
# it applied the item-29 fault must disappear. It is not the fix -- it stalls the
# device on every MultiFab destruction, and we are not maintaining a submodule
# patch. The shipping fix is ERF-side buffer persistence.
#
# Revert after the verification run:
#   git -C Submodules/AMReX checkout Src/Base/AMReX_FabArray.H
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/AMReX/Src/Base/AMReX_FabArray.H"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/amrex-fabarray-clear-sync.patch"
UNPATCHED_SHA=a29eca97ff8165824b1158112b971ad1d0d2ddd6fe2e8bdb06953c2f6fb85383
PATCHED_SHA=9c560fc3edc9dff915b412f1c55bc6e7e92cdbceb17b176632f8b23e383a1509
fail () { echo "FATAL [amrex-clear-sync patch]: $*" >&2; exit 1; }
[ -f "$TARGET" ] || fail "target not found"
cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$cur" = "$PATCHED_SHA" ] && { echo "[amrex-clear-sync patch] already applied (hash verified)."; exit 0; }
[ "$cur" = "$UNPATCHED_SHA" ] || fail "target hashes to $cur, neither known state; re-derive deliberately"
patch -p1 -d "$ERF_ROOT/Submodules/AMReX" < "$PATCH" || fail "patch failed"
new=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$new" = "$PATCHED_SHA" ] || fail "post-patch hash mismatch"
echo "[amrex-clear-sync patch] applied and hash-verified."
