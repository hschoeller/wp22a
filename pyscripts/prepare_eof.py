from dask.distributed import Client, performance_report
from dask.diagnostics import ResourceProfiler, CacheProfiler, Profiler, visualize
import logging
import os
import sys

d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

# Define the lat-lon box bounds
lon_min, lon_max = -100, 20  # Longitude range
lat_min, lat_max = 20, 90    # Latitude range
# bc we have 427 lon and 249 lat
chunkdict = {'time': -1, 'lat': 25, 'lon': 30}


if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))
    import data_utils as du
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

            zg = du.compute_clim(d_path + "data.grib", chunkdict, lat_min,
                                 lat_max, lon_min, lon_max, d_path, mf=False)

            du.compute_anom(zg, chunkdict, d_path)
            du.compute_detrend(d_path, chunkdict)
            du.compute_norm(d_path, chunkdict)
            logging.info("Saving profiler results...")
            visualize([p, rp, cp], file_path="dask_profile.html")

    # Close the Dask client
    logging.info("Job completed. Closing Dask client.")
    client.close()
