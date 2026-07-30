#!/bin/bash
# Pack / verify the data that must reach the pod beyond the git clone.
#
#   bash pod_transfer.sh pack   [outdir]   on the WORKSTATION -- build the tarball
#   bash pod_transfer.sh verify [dir]      on the POD -- check contents + checksums
#
# WHAT SHIPS WITH THE REPO (do NOT put these in the tarball):
#   channel_islands_terrain_3km_192x96.txt   terrain, tracked
#   inputs_c404                              the CONUS404 deck, tracked
#   RRTMGP coefficient tables (46 MB)        in the RRTMGP submodule
#   all scoring/ and harness/ scripts        tracked
#
# WHAT IS NOT INCLUDED AND WHY:
#   wrfout_d02_2020-12-28_00_00_00  26 GB. Transfer separately, and ONLY if you
#       want the d02-dependent parts of the battery: the hourly d02 series and
#       inflow_flux.py's low-level wind / column-vapour comparisons. All
#       accumulation scoring works from wrf_d02_on_grid.npy, which IS included.
#   dem.tif  897 MB. Only regenerates the terrain file, which is already tracked.
set -uo pipefail

SRC=${SRC:-$HOME/ERF}
MODE=${1:-}
DIR=${2:-}

die () { printf '\033[1;31mFATAL: %s\033[0m\n' "$*" >&2; exit 1; }

# path (relative to $SRC) : regenerable-on-pod? : note
ITEMS=(
  "run_c404_nsc/CONUS404Data_3D|no-but-slow|CONUS404 3-D frames, 3-hourly, 2020-12-28 (conus404_to_bin.py can rebuild via GDEX OPeNDAP)"
  "run_c404_nsc/CONUS404Data_Surface|no-but-slow|CONUS404 surface frames (same source)"
  "mrms|yes|MRMS hourly QPE grib2 -- re-downloadable"
  "mrms_20201228_on_grid.npy|yes|MRMS 23-h reference on our grid (chk_mrms.py rebuilds from mrms/)"
  "wrf_d02_on_grid.npy|no|d02 23-h reference on our grid -- ONLY regenerable from the 26 GB wrfout"
  "scoring_ab_dav/lat.npy|no|grid definition -- every scoring script loads it"
  "scoring_ab_dav/lon.npy|no|grid definition"
  "scoring_ab_dav/terrain.npy|no|grid definition + LAND / >500 m masks"
)

case "$MODE" in
# ---------------------------------------------------------------- pack -----
pack)
    OUT=${DIR:-$PWD/erf_pod_transfer}
    mkdir -p "$OUT" || die "cannot create $OUT"
    STAGE=$OUT/payload
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    echo "source: $SRC"
    for entry in "${ITEMS[@]}"; do
        p=${entry%%|*}
        [ -e "$SRC/$p" ] || die "missing on workstation: $SRC/$p"
        mkdir -p "$STAGE/$(dirname "$p")"
        cp -r "$SRC/$p" "$STAGE/$(dirname "$p")/" || die "copy failed: $p"
        printf '  + %-42s %s\n' "$p" "$(du -sh "$SRC/$p" | cut -f1)"
    done

    # flatten the frame dirs to where pod_bootstrap.sh expects them ($DATA)
    mkdir -p "$STAGE/pod_data"
    mv "$STAGE/run_c404_nsc/CONUS404Data_3D" "$STAGE/pod_data/"
    mv "$STAGE/run_c404_nsc/CONUS404Data_Surface" "$STAGE/pod_data/"
    rmdir "$STAGE/run_c404_nsc" 2>/dev/null

    ( cd "$STAGE" && find . -type f -exec sha256sum {} + | sort -k2 > "$OUT/MANIFEST.sha256" )
    n=$(wc -l < "$OUT/MANIFEST.sha256")
    cp "$OUT/MANIFEST.sha256" "$STAGE/MANIFEST.sha256"

    TAR=$OUT/erf_pod_transfer.tar.gz
    tar -czf "$TAR" -C "$STAGE" . || die "tar failed"
    echo
    echo "wrote $TAR  ($(du -h "$TAR" | cut -f1), $n files)"
    echo "manifest: $OUT/MANIFEST.sha256"
    cat <<EOF

Send it:
  scp $TAR pod:/app/ERF/

On the pod:
  mkdir -p /app/ERF/pod_data && tar -xzf /app/ERF/erf_pod_transfer.tar.gz -C /app/ERF
  bash Exec/CanonicalTests/ChannelIslands/pod/pod_transfer.sh verify /app/ERF
EOF
    ;;

# -------------------------------------------------------------- verify -----
verify)
    ROOT=${DIR:-/app/ERF}
    [ -f "$ROOT/MANIFEST.sha256" ] || die "no MANIFEST.sha256 in $ROOT -- untar the transfer set there first"
    cd "$ROOT" || die "cannot cd $ROOT"
    echo "verifying $(wc -l < MANIFEST.sha256) files against MANIFEST.sha256 ..."
    if sha256sum --quiet -c MANIFEST.sha256; then
        echo "all checksums OK"
    else
        die "CHECKSUM MISMATCH -- the transfer is corrupt or incomplete. Re-send."
    fi
    for d in pod_data/CONUS404Data_3D pod_data/CONUS404Data_Surface; do
        n=$(ls "$ROOT/$d" 2>/dev/null | wc -l)
        [ "$n" -gt 0 ] || die "$d is empty"
        printf '  %-34s %d files\n' "$d" "$n"
    done
    echo
    echo "transfer verified. Next: bash Exec/CanonicalTests/ChannelIslands/pod/pod_bootstrap.sh"
    ;;

# ---------------------------------------------------------------- list -----
*)
    echo "usage: pod_transfer.sh {pack [outdir] | verify [dir]}"
    echo
    printf '%-44s %-13s %s\n' "ITEM" "REGENERABLE" "NOTE"
    for entry in "${ITEMS[@]}"; do
        p=${entry%%|*}; rest=${entry#*|}; r=${rest%%|*}; note=${rest#*|}
        sz=$( [ -e "$SRC/$p" ] && du -sh "$SRC/$p" 2>/dev/null | cut -f1 || echo '-' )
        printf '%-44s %-13s %s\n' "$p ($sz)" "$r" "$note"
    done
    exit 1
    ;;
esac
