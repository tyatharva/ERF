import sys
import os
import argparse
from mpi4py import MPI
import glob

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from erftools.preprocessing import Download_ERA5_Data
from erftools.preprocessing import Download_ERA5_ForecastData
from erftools.preprocessing import Download_ERA5_SurfaceData
from erftools.preprocessing import Download_ERA5_ForecastSurfaceData
from erftools.preprocessing import ReadERA5_3DData
from erftools.preprocessing import ReadERA5_SurfaceData

from pyproj import CRS, Transformer
from numpy import *
import time

# ---------------------------------------------------------------------------
#  PINNED PROJECTION
#
#  The DOWNLOAD area and the PROJECTION area are two different things and this
#  domain change is the first time they diverge.
#
#  erftools derives the LCC from whatever box it downloaded, so enlarging the
#  ERA5 box would silently move lat_1/lat_2/lat_0/lon_0 and put every 192x96
#  result on a different map from the 128x64 baseline and from UPSTREAM_ISSUES
#  #25. Nothing would error; the numbers would just stop being comparable.
#
#  So: download the larger box, but build the projection from the ORIGINAL box
#  and refuse to proceed if the emitted string is not character-for-character
#  the one the baseline was run on.
#
#  Why the download box is this much larger than the domain: erftools'
#  find_erf_domain_extents() inscribes an axis-aligned rectangle in the
#  projected lat/lon box and then insets it 100 km on every side, and ERF
#  additionally refuses a domain within 4 cells (12 km) of the frame edge
#  (ERF_WeatherDataInterpolation.cpp:312-323). The 192x96 domain needs frame
#  coverage to x=-388229.74 and y=-166933.25; the original box only reaches
#  x=-265027.4, y=-159204.5.  MEASURED margins with the box below: 37.9 km W,
#  77.3 km E, 42.5 km S, 42.4 km N -- all > 12 km.
# ---------------------------------------------------------------------------
PROJECTION_AREA = [36.0, -123.25, 31.25, -115.25]
PINNED_LCC = ("+proj=lcc +lat_1=32.041667 +lat_2=35.208333 "
              "+lat_0=33.625000 +lon_0=-119.250000 +datum=WGS84 +units=m +no_defs")


def CreateLCCMapping(area):

    lat1 = area[2]
    lat2 = area[0]
    lon1 = area[1]
    lon2 = area[3]

    # Build CRS
    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2

    lat_1 = lat1 + delta/6
    lat_2 = lat2 - delta/6

    lambert_conformal = (
        f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
        f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"
    )

    return lambert_conformal
    
if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    start_time = time.time()

    if len(sys.argv) == 1:
        print("Usage: python3 WriteICFromERA5Data.py <input_filename> [--do_forecast=true]")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Download and process ERA5 data.")
    parser.add_argument("input_filename", help="Input filename, e.g. inputs_for_Laura")

    parser.add_argument(
        "--do_forecast", type=lambda x: x.lower() == "true", default=False,
        help="Set to true to download forecast data"
    )   

    # Forecast args (only required if --do_forecast is set)
    parser.add_argument("--forecast_time_hours", type=int,
                        help="Forecast length in hours, e.g. 72")
    parser.add_argument("--interval_hours", type=int,
                        help="Forecast interval in hours, e.g. 6")

    args = parser.parse_args()

    # Enforce requirement only if do_forecast is true
    if args.do_forecast:
        if args.forecast_time_hours is None or args.interval_hours is None:
            parser.error("--do_forecast requires --forecast_time_hours and --interval_hours")

    input_filename = args.input_filename
    do_forecast = args.do_forecast
    forecast_time_hours = args.forecast_time_hours
    interval_hours = args.interval_hours

    os.path.exists("Output/domain_extents.txt") and os.remove("Output/domain_extents.txt")

    os.makedirs("Output/ERA5Data_3D", exist_ok=True)
    os.makedirs("Output/ERA5Data_Surface", exist_ok=True)
    os.makedirs("Output/VTK/3D/ERA5Domain", exist_ok=True)
    os.makedirs("Output/VTK/3D/ERFDomain", exist_ok=True)
    os.makedirs("Output/VTK/Surface/ERA5Domain", exist_ok=True)
    os.makedirs("Output/VTK/Surface/ERFDomain", exist_ok=True)

   # Download surface data
    if do_forecast:
        if(rank==0):
            print("Running forecast (pressure-level) download + processing...")
        filenames, area = Download_ERA5_ForecastSurfaceData(input_filename, forecast_time_hours, interval_hours)
    else:
        print("Running surface download + processing...")
        filenames, area = Download_ERA5_SurfaceData(input_filename)
    comm.Barrier();

    lambert_conformal = CreateLCCMapping(PROJECTION_AREA)
    if lambert_conformal != PINNED_LCC:
        sys.exit("FATAL: projection drift.\n"
                 f"  emitted: {lambert_conformal}\n"
                 f"  pinned : {PINNED_LCC}\n"
                 "Refusing to write frames on a projection the baseline was not run on.")
    if rank == 0:
        print(f"Projection PINNED to {PROJECTION_AREA} (download area {area}):\n"
              f"  {lambert_conformal}", flush=True)

    # Each rank scans the local directory for matching files
    filenames = sorted(glob.glob("era5_surf_*.grib"))

    my_files = filenames[rank::size]

    for filename in my_files:
        print(f"[Rank {rank}] Processing file: {filename}")
        ReadERA5_SurfaceData(filename, lambert_conformal)

    comm.Barrier();
    if rank == 0:
        print("All ranks finished successfully. Exiting.", flush=True)
 
    # Download 3d data over pressure levels
    if do_forecast:
        filenames, area = Download_ERA5_ForecastData(input_filename, forecast_time_hours, interval_hours)
        comm.Barrier();

        filenames = sorted(glob.glob("era5_3d_*.grib"))

        my_files = filenames[rank::size]
        # Create the directory if it doesn't exist
        for filename in my_files:
            print(f"Processing file: {filename}")
            ReadERA5_3DData(filename, lambert_conformal)

    else:
        filename, area = Download_ERA5_Data(input_filename)
        print("Filename is ", filename)
        print(f"Processing file: {filename}")
        ReadERA5_3DData(filename, lambert_conformal) 

    comm.Barrier()
    
    end_time = time.time()
    elapsed = end_time - start_time

    # Gather or print only on rank 0
    max_elapsed = comm.reduce(elapsed, op=MPI.MAX, root=0)  # max over all ranks

    if rank == 0:
        print(f"Total runtime (wall-clock, across ranks): {max_elapsed:.2f} seconds", flush=True)
