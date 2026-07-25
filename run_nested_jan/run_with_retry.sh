#!/bin/bash
# 24-h nested verification with restart-on-failure (SP nondeterministic NaN,
# UPSTREAM_ISSUES #13: a restart is a RECOVERY splice, documented in RUNBOOK).
cd /app/ERF/run_nested_jan
EXE=/app/ERF/build/Exec/erf_exec
tries=0
while [ $tries -lt 5 ]; do
    latest=$(ls -d chk????? 2>/dev/null | sort | tail -1)
    if [ -n "$latest" ]; then
        RESTART="erf.restart=$latest"
        echo "=== attempt $tries restarting from $latest ==="
    else
        RESTART=""
        echo "=== attempt $tries cold start ==="
    fi
    $EXE inputs_nested $RESTART >> run24.log 2>&1
    ec=$?
    echo "attempt $tries EXIT=$ec" >> run24.log
    if [ $ec -eq 0 ]; then
        echo "COMPLETE after $tries retries" >> run24.log
        exit 0
    fi
    tries=$((tries+1))
done
echo "FAILED after 5 attempts" >> run24.log
exit 1
