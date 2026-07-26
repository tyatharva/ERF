import cdsapi, os
# AREA/OUT are overridable because the 192x96 domain's footprint reaches
# lon -123.49, outside the 128x64 box below; make_sfc_anchor.py's bilinear()
# would clamp the westernmost ~7 columns to the edge value without complaining.
# Defaults reproduce the 128x64 file byte-for-byte.
AREA = [float(v) for v in os.environ.get('AREA', '36.0,-123.25,31.25,-115.25').split(',')]
OUT = os.environ.get('OUT', '/app/ERF/precip_check/era5_sfc_jan.grib')
c = cdsapi.Client()
c.retrieve('reanalysis-era5-single-levels', {
    'product_type': 'reanalysis',
    # surface_pressure and geopotential are the IC hydrostatic anchor: sp gives
    # p at ERA5's own surface, geopotential/g gives the height that surface sits
    # at, so p can be moved to ERF's 3-km terrain height.  The frame binaries
    # cannot supply this -- their bottom level is a byte-identical copy of the
    # 1000 hPa level (erftools height bug), and deriving p from the lowest GOOD
    # frame level instead inherits erftools' ~300 m downward displacement, i.e.
    # ~35 hPa of surface-pressure error.
    'variable': ['surface_latent_heat_flux','surface_sensible_heat_flux',
                 '2m_temperature','sea_surface_temperature',
                 'surface_pressure','geopotential'],
    'year': '2023', 'month': '01', 'day': ['09','10'],
    'time': [f'{h:02d}:00' for h in range(24)],
    'area': AREA,
    'format': 'grib',
}, OUT)
print('DOWNLOAD_OK')
