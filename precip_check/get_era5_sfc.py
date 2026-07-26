import cdsapi
c = cdsapi.Client()
c.retrieve('reanalysis-era5-single-levels', {
    'product_type': 'reanalysis',
    'variable': ['surface_latent_heat_flux','surface_sensible_heat_flux',
                 '2m_temperature','sea_surface_temperature'],
    'year': '2023', 'month': '01', 'day': ['09','10'],
    'time': [f'{h:02d}:00' for h in range(24)],
    'area': [36.0, -123.25, 31.25, -115.25],
    'format': 'grib',
}, '/app/ERF/precip_check/era5_sfc_jan.grib')
print('DOWNLOAD_OK')
