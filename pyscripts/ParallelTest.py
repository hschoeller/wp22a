import xarray as xr
from dask.distributed import Client, performance_report
from dask.diagnostics import ResourceProfiler, CacheProfiler, Profiler, visualize
import logging
import os
from scipy.signal import detrend
import numpy as np

d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

# Define the lat-lon box bounds
lon_min, lon_max = -100, 20  # Longitude range
lat_min, lat_max = 20, 90    # Latitude range
ERA5_DIR = ("/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/" +
            "ERA5/day/atmos/zg/r1i1p1-050000Pa/")
files = [ERA5_DIR + f"zg_day_reanalysis_era5_r1i1p1-050000Pa_{y}0101-{y}1231.nc"
         for y in range(1950, 1954)]
# bc we have 427 lon and 249 lat
chunkdict = {'time': -1, 'lat': 25, 'lon': 30}


def detrend_1d(data):
    # Apply detrend
    return detrend(data, axis=-1)


def normalize_data(data, variable, window):
    std_day = (
        data[variable].groupby("time.dayofyear")
        .std(dim="time")
    )

    extended = xr.concat(
        [std_day.isel(dayofyear=slice(-int(np.ceil(window/2)), None)), std_day,
         std_day.isel(dayofyear=slice(None, int(np.ceil(window/2))))], dim="dayofyear"
    )

    roll_std = extended.rolling(dayofyear=window, center=True).mean()

    cos_lat = np.cos(np.deg2rad(roll_std.lat))

    mean_std = ((
        roll_std.isel(dayofyear=slice(15, 366+15)) * cos_lat)
        .mean(dim=['lat', 'lon']) / cos_lat.mean()
    )

    return data.groupby("time.dayofyear") / mean_std


def compute_clim():
    ds = xr.open_mfdataset(files, parallel=True,
                           combine="by_coords")
    ds = ds.assign_coords(
        lon=(((ds.lon + 180) % 360) - 180)
    )

    # Sort the dataset by longitude for proper slicing
    ds = ds.sortby('lon')
    zg = ds['zg'].sel(lat=slice(lat_min, lat_max),
                      lon=slice(lon_min, lon_max)).chunk(chunkdict)
    logging.info(zg.chunks)

    zg_climatology = zg.groupby("time.dayofyear").mean(dim="time")
    zg_climatology = zg_climatology.persist()
    logging.info(zg_climatology.chunks)
    logging.info("saving to disk")
    zg_climatology.to_netcdf(d_path + "zg_climatology.nc")
    logging.info("Climatology computation completed.")
    return zg


def compute_anom(zg):
    zg_climatology = xr.open_dataset(
        d_path + "zg_climatology.nc", chunks=chunkdict)
    logging.info("Calculate anomalies")
    zg_anom = zg.groupby("time.dayofyear") - zg_climatology
    logging.info("Write anomalies")
    zg_anom.to_netcdf(d_path + "zg_anom.nc", compute=True)


def compute_detrend():
    zg_anom = xr.open_dataset(d_path + "zg_anom.nc",
                              chunks=chunkdict)
    logging.info("Detrending")

    # Apply the function using apply_ufunc
    zg_detrend = xr.apply_ufunc(
        detrend_1d,
        zg_anom,
        input_core_dims=[["time"]],
        output_core_dims=[["time"]],
        dask="parallelized",
        vectorize=True,
        output_dtypes=[float],  # Ensure output dtype is float
        keep_attrs=True
    )

    logging.info(type(zg_detrend))
    zg_detrend = zg_detrend.astype("float32")
    logging.info("Writing detrended")
    zg_detrend.to_netcdf(d_path + "zg_detrend.nc")


def compute_norm():
    zg_detrend = xr.open_dataset(
        d_path + "zg_detrend.nc", chunks=chunkdict)
    logging.info("Normalizing")
    zg_norm = normalize_data(zg_detrend, "zg", window=30)
    logging.info("Writing normalized")
    zg_norm.to_netcdf(d_path + "zg_norm.nc")


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

            zg = compute_clim()

            compute_anom(zg)
            compute_detrend()
            compute_norm()
            logging.info("Saving profiler results...")
            visualize([p, rp, cp], file_path="dask_profile.html")

    # Close the Dask client
    logging.info("Job completed. Closing Dask client.")
    client.close()
