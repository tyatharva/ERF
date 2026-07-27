#!/bin/bash
# DIAGNOSTIC ONLY (UPSTREAM_ISSUES 29i). Adds amrex.the_arena_hunk_size.
#
# CArena sub-allocates from 8 MB hunks, so an overrun from one sub-allocation
# into the next lies INSIDE a valid cudaMalloc and compute-sanitizer is blind to
# it -- which is why eleven sanitizer runs showed only the eventual bad READ.
# With a tiny hunk, CArena's N = max(m_hunk, nbytes) gives every allocation its
# own exact-sized cudaMalloc, so memcheck bounds each buffer individually.
#
# Default 0 = stock behaviour, so the patched build is bit-identical unless the
# knob is set. NOT part of any fix. Revert before shipping:
#   git -C Submodules/AMReX checkout Src/Base/AMReX_Arena.cpp
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/AMReX/Src/Base/AMReX_Arena.cpp"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/amrex-arena-hunk-size.patch"
UNPATCHED_SHA=a5e86f4e17535856e3afb63c3dc80f4272ae12d4e37807d6e08218a19d532ee2
PATCHED_SHA=5eb61d69a83cf0ecd4cef6e23fbb96e771c7def8567fdd0c595e502ac7da50a2
fail () { echo "FATAL [amrex-arena patch]: $*" >&2; exit 1; }
[ -f "$TARGET" ] || fail "target not found"
cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$cur" = "$PATCHED_SHA" ] && { echo "[amrex-arena patch] already applied (hash verified)."; exit 0; }
[ "$cur" = "$UNPATCHED_SHA" ] || fail "target hashes to $cur, neither known state; re-derive deliberately"
patch -p1 -d "$ERF_ROOT/Submodules/AMReX" < "$PATCH" || fail "patch failed"
new=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$new" = "$PATCHED_SHA" ] || fail "post-patch hash mismatch"
echo "[amrex-arena patch] applied and hash-verified."
