#!/bin/bash
# DIAGNOSTIC ONLY (UPSTREAM_ISSUES 29j). Traces CArena alloc/free in a size window.
#
# Item 29 is a use-after-free on a RECYCLED arena block: with
# amrex.the_arena_hunk_size=64 (every allocation its own cudaMalloc, no reuse)
# the fault disappears; with stock 8 MB hunks it kills every 24-h run at the
# seventh frame advance. To name the retained block we need the (ptr, size,
# caller) history of every allocation in the relevant size class, which AMReX
# does not expose -- CArena keeps m_busylist/m_freelist private and the profiler
# hooks report totals, not per-block identity.
#
# Enabled only by AMREX_ARENA_TRACE_LO / AMREX_ARENA_TRACE_HI (bytes) in the
# environment. Unset => the branch is never taken and the build is bit-identical
# to stock, but the STDERR VOLUME when set is large and the fprintf serialises
# allocation, so this must not be present in a production build.
#
# Records go to stderr as:  ARENATRACE {ALLOC|FREE} <ptr> <bytes> [<ret0> <ret1> <ret2>]
# Resolve the return addresses against the executable with addr2line -Cfie.
#
# CAVEAT (cost me a false signal once already): alloc logs the REQUESTED nbytes,
# free logs busy_it->size(), which is the ALIGNED block size. A request just
# below the filter floor is therefore invisible on alloc but visible on free, so
# ALLOC and FREE counts do NOT balance and unmatched frees are expected. Do not
# read that asymmetry as a double free.
#
# NOT part of any fix. Revert before shipping:
#   git -C Submodules/AMReX checkout Src/Base/AMReX_CArena.cpp
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/AMReX/Src/Base/AMReX_CArena.cpp"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/amrex-carena-alloc-trace.patch"
UNPATCHED_SHA=7b715c3f6ccfc5a9309afeca8880b32bcf805a74d480b96b6cd302f7e96e91c4
PATCHED_SHA=352ef2dd2248e3c99e34685208a0ebaa8eee507d6208120ef3ceb26978878fa0
fail () { echo "FATAL [amrex-carena-trace patch]: $*" >&2; exit 1; }
[ -f "$TARGET" ] || fail "target not found"
cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$cur" = "$PATCHED_SHA" ] && { echo "[amrex-carena-trace patch] already applied (hash verified)."; exit 0; }
[ "$cur" = "$UNPATCHED_SHA" ] || fail "target hashes to $cur, neither known state; re-derive deliberately"
patch -p1 -d "$ERF_ROOT/Submodules/AMReX" < "$PATCH" || fail "patch failed"
new=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$new" = "$PATCHED_SHA" ] || fail "post-patch hash mismatch"
echo "[amrex-carena-trace patch] applied and hash-verified."
