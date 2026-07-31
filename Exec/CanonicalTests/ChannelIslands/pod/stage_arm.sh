#!/bin/bash
# Stage ONE production c404 arm directory on the pod (no Docker).
#
#   pod/stage_arm.sh run_c404_sig1
#
# Mirrors exactly what pod_bootstrap.sh section 5 builds for run_validate:
# symlinked frame directories, the 192x96 terrain, the deck, the four RRTMGP
# tables. Refuses to overwrite an existing directory -- run dirs are results.
set -uo pipefail

RUN=${1:?usage: stage_arm.sh <run_dir_name>}
ERF_ROOT=${ERF_ROOT:-/app/ERF}
CI=$ERF_ROOT/Exec/CanonicalTests/ChannelIslands
# DATA / DECK / TERRAIN select the domain. Defaults are the parent
# channelislands-3km-192x96 box, so existing callers are unaffected.
DATA=${DATA:-$ERF_ROOT/pod_data}
DECK=${DECK:-inputs_c404}
TERRAIN=${TERRAIN:-channel_islands_terrain_3km_192x96.txt}
DEST=$ERF_ROOT/$RUN

die () { printf '\nFATAL: %s\n' "$*" >&2; exit 1; }

[ -e "$DEST" ] && die "$DEST already exists -- run dirs are results, not scratch"
for f in CONUS404Data_3D CONUS404Data_Surface; do
    [ -d "$DATA/$f" ] || die "missing $DATA/$f"
done

mkdir -p "$DEST" || die "cannot create $DEST"
cd "$DEST" || die "cannot cd $DEST"

ln -sf "$DATA/CONUS404Data_3D" .
ln -sf "$DATA/CONUS404Data_Surface" .
cp "$CI/$TERRAIN" .   || die "no terrain file $CI/$TERRAIN"
# ERF reads the deck by the name the run command passes, always inputs_c404.
cp "$CI/$DECK" inputs_c404 || die "no deck $CI/$DECK"
# The deck must name the terrain file that was actually staged, or ERF silently
# falls back to the problem's custom terrain (see inputs_c404 line 203).
grep -q "erf.terrain_file_name = \"$TERRAIN\"" inputs_c404 \
    || die "deck $DECK does not name the staged terrain $TERRAIN"
for t in rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/rrtmgp/data/$t" . || die "missing RRTMGP table $t"
done
for t in rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc; do
    cp "$ERF_ROOT/Submodules/RRTMGP/extensions/cloud_optics/$t" . || die "missing $t"
done

# The frame count ERF will index positionally from t=0 must cover the window.
n3=$(ls "$DATA/CONUS404Data_3D"/*.bin | wc -l)
ns=$(ls "$DATA/CONUS404Data_Surface"/*.bin | wc -l)
[ "$n3" -eq "$ns" ] || die "frame count mismatch: 3D=$n3 Surface=$ns"
echo "staged $DEST  ($n3 3D + $ns surface frames)"
