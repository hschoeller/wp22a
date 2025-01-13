# Script to calculate eofs from control ensemble member of ERA5
# eda.
import logging
from dask.diagnostics import ResourceProfiler, CacheProfiler, Profiler, visualize
from dask.distributed import Client, performance_report
import sys
import os
import data_utils as du
sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))


d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

lon_min, lon_max = -100, 20  # Longitude range
lat_min, lat_max = 20, 90    # Latitude range
year_min = 2000
year_max = 2020
chunkdict = {'time': -1, 'lat': 25, 'lon': 30}

if __name__ == "__main__":
    logging.basicConfig(
        filename="dask_resource_usage.log",
        level=logging.INFO,
        format="%(asctime)s - %(message)s"
    )

    logging.info(os.cpu_count())
    client = Client(n_workers=40, threads_per_worker=1, nanny=True)

    logging.info(client)
    logging.info(f"Dask dashboard available at: {client.dashboard_link}")

    logging.info("Starting computations...")
    with performance_report(filename="dask_performance_report.html"):
        with ResourceProfiler() as rp, Profiler() as p, CacheProfiler() as cp:
            results = du.retrieve_ens(year_min, year_max, d_path, lat_min,
                                      lat_max, lon_min, lon_max)

            logging.info("Saving profiler results...")
            visualize([p, rp, cp], file_path="dask_profile.html")

    # Close the Dask client
    logging.info("Job completed. Closing Dask client.")
    client.close()
