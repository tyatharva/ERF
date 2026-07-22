# ============================================================================
# ERF Channel Islands hindcast environment
#   erf-env (CUDA 12.6.2 + MPICH + HDF5 + NetCDF-C + CMake)
#   + gfortran & NetCDF-Fortran   (required by ERF_ENABLE_NOAHMP)
#   + erftools + Python deps     (ERA5 -> ERF IC/BC pipeline)
#   + RRTMGP coefficient tables  (fetched at build time)
#
# Build:  docker build -t erf-hindcast -f Dockerfile .
# Run:    docker run --rm --gpus all -v ~/ERF:/app/ERF \
#             -v ~/.cdsapirc:/root/.cdsapirc:ro -w /app/ERF erf-hindcast <cmd>
# NOTE: the CDS key is mounted read-only at runtime -- never baked into a layer.
# ============================================================================
FROM erf-env

# ---- Fortran toolchain + NetCDF-Fortran (built against /usr/local netcdf-c) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran m4 curl ca-certificates python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN cd /tmp \
    && curl -fsSL -o netcdf-fortran.tar.gz \
       https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz \
    && tar xf netcdf-fortran.tar.gz && cd netcdf-fortran-4.6.1 \
    && FC=gfortran CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --disable-static \
    && make -j"$(nproc)" && make install && ldconfig \
    && cd / && rm -rf /tmp/netcdf-fortran*

# ---- Python deps for erftools (ERA5 pipeline) and the terrain converter ----
RUN python3 -m pip install --no-cache-dir \
        cdsapi xarray cfgrib eccodes netCDF4 numpy scipy pandas \
        pyproj rasterio mpi4py matplotlib

# ---- erftools (ERA5 -> ERF IC/BC preprocessing) ----
RUN git clone --depth 1 https://github.com/erf-model/erftools.git /opt/erftools \
    && python3 -m pip install --no-cache-dir -e /opt/erftools || true
ENV ERFTOOLS_HOME=/opt/erftools

# ---- RRTMGP coefficient tables (gas optics + cloud optics) ----
# E3SM inputdata mirror; the four defaults expected by ERF_Radiation.H.
RUN mkdir -p /opt/rrtmgp-data && cd /opt/rrtmgp-data \
    && for f in rrtmgp-data-sw-g224-2018-12-04.nc \
                rrtmgp-data-lw-g256-2018-12-04.nc \
                rrtmgp-cloud-optics-coeffs-sw.nc \
                rrtmgp-cloud-optics-coeffs-lw.nc ; do \
         curl -fsSL -o $f \
           "https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/cam/rad/$f" \
         || echo "WARNING: could not fetch $f -- supply manually"; \
       done
ENV RRTMGP_DATA=/opt/rrtmgp-data

WORKDIR /app/ERF
