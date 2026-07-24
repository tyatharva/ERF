#!/usr/bin/env python3
"""Add ERA5 forecast albedo (fal) to the erftools surface-BC pipeline.

The production ERF config feeds RRTMGP a per-column, time-varying surface
albedo from the ERA5 forecast-albedo field (MODIS-derived climatological
annual cycle + model snow). Validated over this domain (2026-07-24): land
mean 0.178-0.194 through the year (expected SoCal 0.15-0.20), spatial
p10-p90 = 0.11-0.26, ocean 0.061, winter snow max 0.52.

This script patches /opt/erftools' ReadERA5DataAndWriteERF_SurfBC.py to:
  1. request "forecast_albedo" in the CDS surface download,
  2. extract it from the GRIBs (name match: "Forecast albedo"),
  3. append it as the SIXTH field of the surface .bin frames
     (after sea_surf_temp, surf_sen_hf, surf_lat_hf, fric_vel, ls_mask --
     positions of the existing five fields are unchanged, and the .bin
     header's ndata lets the ERF reader detect old 5-field files).

Idempotent; HARD-FAILS if the expected source lines are not found (an
erftools update that moves this code is caught, not papered over).
NOTE: cached era5 surface GRIBs downloaded before this patch do NOT
contain fal -- delete them so the request re-downloads with the field.

Usage: python3 erftools_fal_patch.py /opt/erftools
"""
import sys

root = sys.argv[1] if len(sys.argv) > 1 else "/opt/erftools"
target = f"{root}/erftools/preprocessing/era5/ReadERA5DataAndWriteERF_SurfBC.py"

with open(target) as f:
    src = f.read()

replacements = [
    # 1a. CDS request list, Download_ERA5_SurfaceData
    ('            "land_sea_mask",\n            "surface_pressure"  # <- added',
     '            "land_sea_mask",\n            "forecast_albedo",  # <- fal patch\n            "surface_pressure"  # <- added'),
    # 1b. CDS request list, Download_ERA5_ForecastSurfaceData (the one the
    #     production --do_forecast pipeline actually calls; different
    #     formatting, so it needs its own literal)
    ('        "land_sea_mask",\n        "surface_pressure"\n    ]',
     '        "land_sea_mask",\n        "forecast_albedo",  # <- fal patch\n        "surface_pressure"\n    ]'),
    # 2. GRIB accumulator declaration
    ('    ls_mask_era5 = []',
     '    ls_mask_era5 = []\n    forecast_albedo_era5 = []'),
    # 3. GRIB message match
    ('            if "Land-sea mask" in grb.name:\n                ls_mask_era5.append(grb.values)',
     '            if "Land-sea mask" in grb.name:\n                ls_mask_era5.append(grb.values)\n\n'
     '            if "albedo" in grb.name.lower():\n                forecast_albedo_era5.append(grb.values)'),
    # 4. Stack alongside the others
    ('    ls_mask_era5 = np.stack(ls_mask_era5, axis=0)',
     '    ls_mask_era5 = np.stack(ls_mask_era5, axis=0)\n'
     '    forecast_albedo_era5 = np.stack(forecast_albedo_era5, axis=0)'),
    # 5. Per-level extraction + array
    ('    ls_mask = np.zeros((nx, ny, nz))',
     '    ls_mask = np.zeros((nx, ny, nz))\n    forecast_albedo = np.zeros((nx, ny, nz))'),
    ('        ls_mask_at_lev = ls_mask_era5[k]',
     '        ls_mask_at_lev = ls_mask_era5[k]\n        forecast_albedo_at_lev = forecast_albedo_era5[k]'),
    ('        ls_mask[:, :, k] = ls_mask_at_lev',
     '        ls_mask[:, :, k] = ls_mask_at_lev\n        forecast_albedo[:, :, k] = forecast_albedo_at_lev'),
]

for old, new in replacements:
    if new in src:
        continue  # already patched
    if old not in src:
        sys.exit(f"FATAL [erftools fal patch]: expected line not found in "
                 f"{target}:\n---\n{old}\n---\nerftools has changed; re-derive this fix.")
    src = src.replace(old, new)

# 6. Append to BOTH scalars dicts (ERA5-grid writer and ERF-grid writer).
#    Same literal appears twice; str.replace patches both.
old_dicts = src.count('         "ls_mask": ls_mask,\n    }')
if '"forecast_albedo": forecast_albedo,' not in src:
    if old_dicts < 1:
        sys.exit("FATAL [erftools fal patch]: scalars dict pattern not found; re-derive.")
    src = src.replace('         "ls_mask": ls_mask,\n    }',
                      '         "ls_mask": ls_mask,\n         "forecast_albedo": forecast_albedo,\n    }')
    # ERF-grid dict uses the *_erf arrays; add the target array + dict entry
    if '    ls_mask_erf = np.zeros((nx_erf, ny_erf, nz_erf))' in src:
        src = src.replace('    ls_mask_erf = np.zeros((nx_erf, ny_erf, nz_erf))',
                          '    ls_mask_erf = np.zeros((nx_erf, ny_erf, nz_erf))\n'
                          '    forecast_albedo_erf = np.zeros((nx_erf, ny_erf, nz_erf))')
    else:
        sys.exit("FATAL [erftools fal patch]: ls_mask_erf allocation not found; re-derive.")
    if '         "ls_mask": ls_mask_erf,\n    }' in src:
        src = src.replace('         "ls_mask": ls_mask_erf,\n    }',
                          '         "ls_mask": ls_mask_erf,\n'
                          '         "forecast_albedo": forecast_albedo_erf,\n    }')
    else:
        sys.exit("FATAL [erftools fal patch]: ERF-grid scalars dict not found; re-derive.")

with open(target, "w") as f:
    f.write(src)

# Verify
with open(target) as f:
    out = f.read()
need = ["forecast_albedo_era5.append(grb.values)",
        '"forecast_albedo",  # <- fal patch',
        '"forecast_albedo": forecast_albedo,',
        '"forecast_albedo": forecast_albedo_erf,']
for n in need:
    if n not in out:
        sys.exit(f"FATAL [erftools fal patch]: verification failed for: {n}")
print("[erftools fal patch] applied and verified.")
