from dask.distributed import Client, performance_report
from dask.diagnostics import (ResourceProfiler, CacheProfiler, Profiler,
                              visualize)
import logging
import os
import sys

d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

# Define the lat-lon box bounds
lon_min, lon_max = -100, 20  # Longitude range
lat_min, lat_max = 20, 90    # Latitude range
parallel = False
if parallel:
    chunkdict = {'time': -1, 'lat': 50, 'lon': 30}
else:
    chunkdict = {'time': -1, 'lat': 'auto', 'lon': 'auto'}

if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))
    import data_utils as du
    logging.basicConfig(
        filename="dask_resource_usage.log",
        level=logging.INFO,
        format="%(asctime)s - %(message)s"
    )

    logging.info(os.cpu_count())
    if parallel:
        client = Client(n_workers=12, threads_per_worker=1, nanny=True)

        logging.info(client)
        logging.info(f"Dask dashboard available at: {client.dashboard_link}")

    logging.info("Member %s", sys.argv[1])
    logging.info("Starting computations...")
    if parallel:
        with performance_report(filename="dask_performance_report.html"):
            with (ResourceProfiler() as rp,
                  Profiler() as p, CacheProfiler() as cp):

                zg = du.compute_clim(d_path + "data.grib", chunkdict, lat_min,
                                     lat_max, lon_min, lon_max,
                                     d_path, mf=False)

                du.compute_anom(zg, chunkdict, d_path)
                du.compute_detrend(d_path, chunkdict)
                du.compute_norm(d_path, chunkdict)
                logging.info("Saving profiler results...")
                visualize([p, rp, cp], file_path="dask_profile.html")

        # Close the Dask client
        logging.info("Job completed. Closing Dask client.")
        client.close()
    else:
        # zg = du.compute_clim(d_path + "data.grib", chunkdict, lat_min,
        #  lat_max, lon_min, lon_max, d_path, mf=False,
        #  member=int(sys.argv[1]))

        # du.compute_anom(zg, chunkdict, d_path, member=int(sys.argv[1]))
        du.compute_detrend(d_path, chunkdict, member=int(sys.argv[1]))
        du.compute_norm(d_path, chunkdict, member=int(sys.argv[1]))
        logging.info("Job completed.")
