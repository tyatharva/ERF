#!/bin/bash
# DIAGNOSTIC ONLY (UPSTREAM_ISSUES 29k). Validates the CACHED FillBoundary tag
# vector on every cache hit.
#
# FB_get_local_copy_tag_vector caches a TagVector per FB id in the FabArray's
# m_fb_local_copy_handler, holding raw Array4 fab pointers for the lifetime of the
# FabArray. Item 29 needed to distinguish two failure modes that look identical
# from the crash site (a wild address inside FB_local_copy_gpu):
#
#   device == cached but != fresh  -> cached tags are STALE (the fabs moved)
#   device != cached               -> the device buffer was SCRIBBLED by whoever
#                                     was handed its arena block
#
# On a hit this rebuilds the tags from the CURRENT fabs, copies the device buffer
# back, and prints on mismatch. VERDICT: scribbled. The device tag array held
# 0x43905337 / 0x439051b1 / 0x43900360 -- float32 288.65 / 288.64 / 288.03 K --
# where fab pointers belong.
#
# READING THE OUTPUT: sizeof(Array4CopyTag<float,float>) is 176 with a 4-byte
# PADDING HOLE at offset 68 (between dindex at 64..67 and sfab at 72). define()
# memcpy's the struct including padding, so a mismatch reported at
# first_diff_byte=68 is padding noise, NOT corruption. Only diffs at meaningful
# offsets count -- the real one reports first_diff_byte=0.
#
# COST: a dtoh_memcpy plus a full host-side tag rebuild on EVERY FillBoundary
# cache hit. This is enormously slow and MUST NOT SHIP.
#
# Revert before any scored or production run:
#   git -C Submodules/AMReX checkout Src/Base/AMReX_FBI.H
set -u
ERF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARGET="$ERF_ROOT/Submodules/AMReX/Src/Base/AMReX_FBI.H"
PATCH="$ERF_ROOT/Exec/CanonicalTests/ChannelIslands/patches/amrex-fbi-tagcheck.patch"
UNPATCHED_SHA=c0ebf690c972bb696103db407931977106826b0fcf2094fc3424147b3ab4a826
PATCHED_SHA=ff3b5ca225552f9b3ce37408f299f87692102963796c4e2f49fe6e4f0999e58d
fail () { echo "FATAL [amrex-fbi-tagcheck patch]: $*" >&2; exit 1; }
[ -f "$TARGET" ] || fail "target not found"
cur=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$cur" = "$PATCHED_SHA" ] && { echo "[amrex-fbi-tagcheck patch] already applied (hash verified)."; exit 0; }
[ "$cur" = "$UNPATCHED_SHA" ] || fail "target hashes to $cur, neither known state; re-derive deliberately"
patch -p1 -d "$ERF_ROOT/Submodules/AMReX" < "$PATCH" || fail "patch failed"
new=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ "$new" = "$PATCHED_SHA" ] || fail "post-patch hash mismatch"
echo "[amrex-fbi-tagcheck patch] applied and hash-verified."
