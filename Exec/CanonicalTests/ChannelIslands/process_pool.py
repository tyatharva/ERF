#!/usr/bin/env python3
"""Incremental, idempotent erftools preprocessing of the flat ERA5 GRIB pool.

Run inside the container, from the pool directory:

    mpirun -n 8 python3 process_pool.py era5_input.txt

WHY NOT WriteICFromERA5Data.py: that driver (a) always calls the CDS
download functions first, so pointing it at a partially-downloaded pool
makes it issue per-timestep CDS requests -- exactly what
era5_batch_download.py exists to avoid -- and (b) reprocesses every GRIB
in the glob on every invocation, with no check on whether the .bin
already exists. Both make it unusable for processing-while-downloading.

This calls the same erftools readers with the same LCC mapping, but:
  * never downloads anything;
  * skips any GRIB whose .bin output already exists;
  * takes whatever is on disk right now, per stream independently.

So it can be run repeatedly while the download proceeds, and each frame
is converted exactly once. Output layout is identical to the driver's:
Output/ERA5Data_3D/ERF_IC_<stamp>.bin and Output/ERA5Data_Surface/...
"""
import glob
import os
import re
import sys
import time

from mpi4py import MPI
from erftools.preprocessing import ReadERA5_3DData, ReadERA5_SurfaceData

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()


def lcc_from_area(area):
    """Verbatim reproduction of CreateLCCMapping in WriteICFromERA5Data.py.

    area is the era5_input.txt list [north, west, south, east]; keep the
    index order and the .6f formatting identical so frames produced here
    are bit-compatible with frames produced by the stock driver.
    """
    lat1, lat2 = area[2], area[0]
    lon1, lon2 = area[1], area[3]
    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2
    lat_1 = lat1 + delta / 6
    lat_2 = lat2 - delta / 6
    return (f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
            f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs")


def expected_bin(grib, outdir):
    """era5_3d_20230820_0300.grib -> Output/<outdir>/ERF_IC_2023_08_20_03_00.bin"""
    m = re.search(r"_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})\.grib$", grib)
    if not m:
        return None
    return os.path.join("Output", outdir, "ERF_IC_{}_{}_{}_{}_{}.bin".format(*m.groups()))


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "era5_input.txt"
    area = None
    for line in open(input_file):
        if line.strip().startswith("area:"):
            area = [float(v) for v in line.split(":", 1)[1].split(",")]
    if area is None:
        sys.exit("no 'area:' line in " + input_file)
    lcc = lcc_from_area(area)

    streams = [("era5_3d_*.grib", "ERA5Data_3D", ReadERA5_3DData),
               ("era5_surf_*.grib", "ERA5Data_Surface", ReadERA5_SurfaceData)]

    # The readers also emit VTK next to the .bin and do NOT create these
    # directories themselves -- the stock driver does it for them.
    for d in ("Output/VTK/3D/ERA5Domain", "Output/VTK/3D/ERFDomain",
              "Output/VTK/Surface/ERA5Domain", "Output/VTK/Surface/ERFDomain"):
        os.makedirs(d, exist_ok=True)

    for pattern, outdir, reader in streams:
        os.makedirs(os.path.join("Output", outdir), exist_ok=True)
        todo = [g for g in sorted(glob.glob(pattern))
                if not os.path.exists(expected_bin(g, outdir))]
        if rank == 0:
            have = len(glob.glob(pattern))
            print(f"[{outdir}] {have} gribs on disk, {len(todo)} to convert",
                  flush=True)
        comm.Barrier()
        t0 = time.time()
        for g in todo[rank::size]:
            reader(g, lcc)
        comm.Barrier()
        if rank == 0 and todo:
            dt = time.time() - t0
            print(f"[{outdir}] converted {len(todo)} frames in {dt:.1f} s "
                  f"({dt / len(todo) * size:.2f} s/frame/rank)", flush=True)

    if rank == 0:
        for _, outdir, _ in streams:
            n = len(glob.glob(os.path.join("Output", outdir, "*.bin")))
            print(f"[{outdir}] total frames now {n}", flush=True)


if __name__ == "__main__":
    main()
