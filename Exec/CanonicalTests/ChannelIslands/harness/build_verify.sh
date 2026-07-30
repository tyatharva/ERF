#!/bin/bash
# Build the rad-clock (#26) + NSCBC-regime (#19d) changes and verify them on the
# BINARY, not on the build exit code.
#
# Two traps this script has already fallen into, both fixed here:
#  * `cmake --build . | tail` returns TAIL's status, so a failed build reported
#    exit 0. No pipe.
#  * kokkos' AlwaysCheckGit target runs git inside the container, where the
#    bind-mounted repo is owned by the host user -> "dubious ownership" -> the
#    build dies at 47%. Mark it safe for the (throwaway) container.
set -u
SP=/tmp/claude-1000/-home-atyagi-ERF/bde1210f-5368-4296-897d-ecbd7c312e37/scratchpad
BIN=/home/atyagi/ERF/build/Exec/erf_exec
BEFORE=$(stat -c %Y "$BIN" 2>/dev/null || echo 0)

docker run --rm -v /home/atyagi/ERF:/app/ERF -w /app/ERF/build erf-hindcast \
  bash -c 'git config --global --add safe.directory "*"; cmake --build . -j 12' \
  > "$SP/build.log" 2>&1
EX=$?
echo "BUILD exit=$EX"
if [ "$EX" != 0 ]; then
  echo "BUILD FAILED -- last errors:"; grep -a -i -B2 -A4 "error" "$SP/build.log" | tail -30; exit 2
fi

AFTER=$(stat -c %Y "$BIN")
echo "binary mtime: before=$BEFORE after=$AFTER"
[ "$AFTER" != "$BEFORE" ] || { echo "ARTIFACT UNCHANGED -- the build did not relink"; exit 2; }
ls -la --time-style=+%H:%M "$BIN"
tail -3 "$SP/build.log"
echo BUILD_DONE
