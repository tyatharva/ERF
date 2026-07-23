#!/usr/bin/env python3
"""Compatibility fix for erftools vs modern eccodes GRIB tables.

Newer eccodes renamed the accumulated surface fluxes:
    "Surface sensible heat flux" -> "Time-integrated surface sensible heat net flux"
    "Surface latent heat flux"   -> "Time-integrated surface latent heat net flux"
erftools' ReadERA5_SurfaceData matches the old names by substring, collects
nothing, and crashes with "need at least one array to stack".

This script rewrites the two matches to a case-insensitive core that accepts
BOTH namings, and HARD-FAILS if the expected source lines are not found
(so an erftools update that changes this code is caught, not papered over).

Usage: python3 erftools_eccodes_compat.py /opt/erftools
"""
import sys

root = sys.argv[1] if len(sys.argv) > 1 else "/opt/erftools"
target = f"{root}/erftools/preprocessing/era5/ReadERA5DataAndWriteERF_SurfBC.py"

with open(target) as f:
    src = f.read()

replacements = [
    ('if "Surface sensible heat flux" in grb.name:',
     'if "surface sensible heat" in grb.name.lower():'),
    ('if "Surface latent heat flux" in grb.name:',
     'if "surface latent heat" in grb.name.lower():'),
]

for old, new in replacements:
    if new in src:
        continue  # already patched
    if old not in src:
        sys.exit(f"FATAL [erftools eccodes compat]: expected line not found in "
                 f"{target}:\n  {old}\nerftools has changed; re-derive this fix.")
    src = src.replace(old, new)

with open(target, "w") as f:
    f.write(src)
print("[erftools eccodes compat] flux-name matching patched (or already OK).")
