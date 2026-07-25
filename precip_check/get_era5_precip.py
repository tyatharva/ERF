#!/usr/bin/env python3
"""Fetch ERA5 total_precipitation for the stress-test windows."""
import cdsapi

AREA = [36.0, -123.25, 31.25, -115.25]
c = cdsapi.Client()

for tag, y, m, days in [("aug", "2023", "08", ["20", "21"]),
                        ("jan", "2023", "01", ["09", "10"])]:
    c.retrieve("reanalysis-era5-single-levels", {
        "product_type": "reanalysis",
        "variable": ["total_precipitation", "convective_precipitation",
                     "large_scale_precipitation"],
        "year": y, "month": m, "day": days,
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": AREA,
        "format": "grib",
    }, f"era5_precip_{tag}.grib")
    print(f"{tag}: done")

c.retrieve("reanalysis-era5-single-levels", {
    "product_type": "reanalysis",
    "variable": ["total_precipitation"],
    "year": "2023", "month": "09", "day": ["09","10"],
    "time": [f"{h:02d}:00" for h in range(24)],
    "area": AREA, "format": "grib"}, "era5_precip_sep.grib")
print("sep: done")
