# ============================================================================
# ERF Channel Islands hindcast environment
#   erf-env (CUDA 12.6.2 + MPICH + HDF5 + NetCDF-C + CMake)
#   + erftools + Python deps     (ERA5 -> ERF IC/BC pipeline)
#   + RRTMGP coefficient tables  (fetched at build time)
#
# Build:  docker build -t erf-hindcast -f Dockerfile .
# Run:    docker run --rm --gpus all -v ~/ERF:/app/ERF \
#             -v ~/.cdsapirc:/root/.cdsapirc:ro -w /app/ERF erf-hindcast <cmd>
# NOTE: the CDS key is mounted read-only at runtime -- never baked into a layer.
# ============================================================================
FROM erf-env

# ---- Minimal additions over erf-env (MM5 land model needs no Fortran) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates python3-pip \
    && rm -rf /var/lib/apt/lists/*

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
