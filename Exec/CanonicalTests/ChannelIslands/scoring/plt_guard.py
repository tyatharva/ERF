#!/usr/bin/env python3
"""Shared plotfile guard for the scoring scripts.

RECONSTRUCTED 2026-08-04. The original was committed in d17e6dfb ("Close the
plotfile guard gap: one shared module, all seven scoring scripts") but git never
actually tracked it: .gitignore carried a bare `plt*`, meant for AMReX plotfile
DIRECTORIES, which also matched `plt_guard.py`. The commit message describes a
file that was never in the tree, so the module survived only on the rented pod
and every one of its fourteen importers was broken on any fresh clone. The
ignore patterns are now anchored on digits (plt[0-9]* / chk[0-9]*).

This reconstruction implements the two rules the original commit message states.
It is behaviourally equivalent for the cases those scripts exercise, but it is
NOT byte-recovered -- if the original rejected on some signature beyond
non-finite and negative accumulation, that extra strictness is not reproduced.
Both rules below are quoted from the commit that introduced it:

  "reject on the DATA, never on dt. The degenerate dt is not a fixed value
   (2^-25 on items 52-55, 2^-26 on the item-29-era width runs) and an AMReX
   plotfile Header carries no dt field at all"

  "order by TIME, never by name, so first-wins dedup cannot be inverted by
   lexicographic step numbering"

The second rule is not hypothetical: sorted(glob('plt*')) puts 'plt10000'
before 'plt9999', so any first-wins hour dedup silently inverts once a run
crosses a digit-width boundary. A 29 h run at dt 0.5 reaches step ~209000, so
it crosses that boundary four times.
"""
import glob
import os

import numpy as np

# CORRUPTION sentinels, NOT physical plausibility bounds.
#
# The first version of this reconstruction used a 1e4 mm ceiling as a "no run
# produces 400 inches of rain" check. That was WRONG and it fired on a healthy
# run: the forcing-dominated lateral band routinely reaches such values --
# measured 12219 mm of cumulative rain_accum at (191,2), i.e. ON the xhi wall,
# in a 12 h compressible Davies run whose INTERIOR max was 54 mm and whose field
# was entirely finite. Every cell above 1000 mm was inside the 10-cell band that
# every scored number in this campaign already excludes.
#
# So magnitude does not discriminate corruption from the known band artifact,
# and a guard that rejects healthy files is worse than no guard -- it trains you
# to loosen it. The reliable signature of the degenerate-dt runs this module
# exists to catch is NON-FINITE data. The bounds below stay only to catch gross
# nonsense (uninitialised memory, a mis-parsed field), several orders of
# magnitude beyond anything the band produces.
_ACCUM_FLOOR = -1.0e-3      # mm; round-off below zero only
_ACCUM_CEIL = 1.0e7         # mm; ~1000x the largest observed band value


def plotfile_time(path):
    """Model time from an AMReX plotfile Header, without loading the data.

    Layout (verified against this campaign's plotfiles):
        0        version string
        1        ncomp
        2..2+n-1 variable names
        2+ncomp  spatial dimension
        3+ncomp  time
    """
    with open(os.path.join(path, 'Header')) as fh:
        lines = fh.read().split('\n')
    ncomp = int(lines[1])
    return float(lines[3 + ncomp])


def plotfiles_by_time(rundir):
    """Every plotfile in rundir, ordered by MODEL TIME.

    Never by name -- see the module docstring. Files whose Header cannot be
    read are dropped rather than sorted to an arbitrary position: a partially
    written plotfile (the run was killed mid-write) has no defensible time.
    """
    out = []
    for p in glob.glob(os.path.join(rundir, 'plt[0-9]*')):
        if not os.path.isdir(p):
            continue
        try:
            out.append((plotfile_time(p), p))
        except (OSError, ValueError, IndexError):
            continue
    return [p for _, p in sorted(out)]


def is_poisoned(arr):
    """True if this field cannot be scored.

    Rejects on the DATA. A run that has gone degenerate produces non-finite
    values; nanmean() would hide that and return a plausible number from the
    surviving cells, which is exactly how the defect this module exists to
    prevent went unnoticed.
    """
    a = np.asarray(arr, dtype=np.float64)
    if not np.isfinite(a).all():
        return True
    if a.min() < _ACCUM_FLOOR or a.max() > _ACCUM_CEIL:
        return True
    return False


def reject_if_poisoned(label, path, arr):
    """Abort loudly. Callers must NOT be able to continue on a poisoned file."""
    a = np.asarray(arr, dtype=np.float64)
    if is_poisoned(a):
        nbad = int((~np.isfinite(a)).sum())
        raise SystemExit(
            f'POISONED plotfile: {label} at {path}\n'
            f'  non-finite cells {nbad} of {a.size}, '
            f'range {np.nanmin(a):.4g} .. {np.nanmax(a):.4g}\n'
            f'  Refusing to score. This is the guard from UPSTREAM_ISSUES item '
            f'55b -- a degraded file must stop the run, not average away.')
