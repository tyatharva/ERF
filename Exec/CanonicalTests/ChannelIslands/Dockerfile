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
        gfortran m4 csh curl ca-certificates python3-pip file \
    && rm -rf /var/lib/apt/lists/*

RUN cd /tmp \
    && curl -fsSL -o netcdf-fortran.tar.gz \
       https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz \
    && tar xf netcdf-fortran.tar.gz && cd netcdf-fortran-4.6.1 \
    && FC=gfortran CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --disable-static \
    && make -j"$(nproc)" && make install && ldconfig \
    && cd / && rm -rf /tmp/netcdf-fortran*

# ---- WRF (real.exe) + WPS, prebuilt for the NOAHMP wrfinput pipeline ----
# Serial GNU builds; we only need geogrid/ungrib/metgrid/real, never wrf.exe.
ENV NETCDF=/usr/local
RUN cd /opt \
    && curl -fsSL -o wrf.tar.gz \
         https://github.com/wrf-model/WRF/releases/download/v4.6.1/v4.6.1.tar.gz \
    && mkdir WRF && tar xf wrf.tar.gz -C WRF --strip-components=1 && rm wrf.tar.gz \
    && cd WRF && printf '32\n1\n' | ./configure \
    && { ./compile em_real > /tmp/wrf_compile.log 2>&1 || true; } \
    && { test -x main/real.exe || { echo "FATAL: real.exe did not build; last errors:" >&2; \
         grep -nE "error|Error" /tmp/wrf_compile.log | tail -20 >&2; exit 1; }; } \
    && rm -f /tmp/wrf_compile.log && find . -name "*.o" -delete
RUN cd /opt \
    && curl -fsSL -o wps.tar.gz \
         https://github.com/wrf-model/WPS/archive/refs/tags/v4.6.0.tar.gz \
    && mkdir WPS && tar xf wps.tar.gz -C WPS --strip-components=1 && rm wps.tar.gz \
    && cd WPS && printf '1\n' | ./configure \
    && sed -i "s|^WRF_DIR.*|WRF_DIR = /opt/WRF|" configure.wps \
    && ./compile 2>&1 | tail -4 \
    && test -x geogrid/src/geogrid.exe || { echo "FATAL: geogrid.exe did not build" >&2; exit 1; } \
    && test -x ungrib/src/ungrib.exe  || { echo "FATAL: ungrib.exe did not build"  >&2; exit 1; } \
    && test -x metgrid/src/metgrid.exe || { echo "FATAL: metgrid.exe did not build" >&2; exit 1; }
ENV WRF_DIR=/opt/WRF WPS_DIR=/opt/WPS

# ---- Python deps for erftools (ERA5 pipeline) and the terrain converter ----
# Upgrade pip/setuptools first: the base image's setuptools cannot parse
# erftools' PEP-621 pyproject (it built an empty "UNKNOWN-0.0.0" wheel).
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && python3 -m pip install --no-cache-dir \
        cdsapi xarray cfgrib eccodes netCDF4 numpy scipy pandas \
        pyproj rasterio mpi4py matplotlib \
        h5netcdf bottleneck pydantic f90nml click ussa1976 cartopy \
        herbie-data pygrib tqdm python-dateutil yt

# ---- erftools (ERA5 -> ERF IC/BC preprocessing) ----
# Used from the source tree via PYTHONPATH: erftools' wheel packaging omits
# its subpackages (pip install produces a package without .preprocessing),
# so a pip install would "succeed" and then fail at import time. The import
# self-test fails the image build if the source-tree import breaks.
RUN git clone --depth 1 https://github.com/erf-model/erftools.git /opt/erftools
COPY erftools_eccodes_compat.py /opt/erftools_eccodes_compat.py
RUN python3 /opt/erftools_eccodes_compat.py /opt/erftools
ENV ERFTOOLS_HOME=/opt/erftools PYTHONPATH=/opt/erftools
RUN python3 -c "import erftools.preprocessing; print('erftools import OK')"

# ---- RRTMGP coefficient tables ----
# The four tables ERF needs ship inside the ERF repo itself (bind-mounted at
# runtime), so nothing is fetched here:
#   Submodules/RRTMGP/rrtmgp/data/rrtmgp-data-{sw-g224,lw-g256}-2018-12-04.nc
#   Submodules/RRTMGP/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-{sw,lw}.nc
# Run scripts copy them into the run directory (see benchmark instructions).



WORKDIR /app/ERF
