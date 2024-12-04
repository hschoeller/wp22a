import numpy as np
import xarray as xr
from scipy.signal import detrend
import logging
import cdsapi
import zipfile
import dask


def retrieve_ens(years, lat_min, lat_max,
                 lon_min, lon_max, d_path, f_name):
    # download data

    c = cdsapi.Client()

    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": ["ensemble_members"],
        "variable": ["geopotential"],
        "year": years,
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "day": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"
        ],
        "time": ["12:00"],
        "pressure_level": ["500"],
        "data_format": "grib",
        "download_format": "zip",
        "area": [lat_max, lon_min, lat_min, lon_max]
    }

    c.retrieve(dataset, request, d_path + f_name)
    # Unzip the downloaded file
    with zipfile.ZipFile(f"{d_path}/{f_name}", 'r') as zip_ref:
        zip_ref.extractall(d_path)


def detrend_1d(data):
    # Apply detrend
    return detrend(data, axis=-1)


def normalize_data(data, variable, window):
    std_day = (
        data[variable].groupby("time.dayofyear")
        .std(dim="time")
    )

    extended = xr.concat(
        [std_day.isel(dayofyear=slice(-int(np.ceil(window/2)), None)),
         std_day,
         std_day.isel(dayofyear=slice(None, int(np.ceil(window/2))))],
        dim="dayofyear"
    )

    roll_std = extended.rolling(dayofyear=window, center=True).mean()

    cos_lat = np.cos(np.deg2rad(roll_std.lat))

    mean_std = ((
        roll_std.isel(dayofyear=slice(15, 366+15)) * cos_lat)
        .mean(dim=['lat', 'lon']) / cos_lat.mean()
    )

    return data.groupby("time.dayofyear") / mean_std


def compute_clim(files, chunkdict, lat_min, lat_max, lon_min, lon_max,
                 d_path, mf=True):
    if mf:
        ds = xr.open_mfdataset(files, parallel=True,
                               combine="by_coords",
                               chunks=chunkdict)
    else:
        ds = xr.open_dataset(files, engine="cfgrib",
                             chunks=chunkdict).isel(number=0)
        ds = ds.rename({'z': 'zg', 'longitude': 'lon',
                        'latitude': 'lat'})
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


def compute_anom(zg, chunkdict, d_path):
    zg_climatology = xr.open_dataset(
        d_path + "zg_climatology.nc", chunks=chunkdict)
    logging.info("Calculate anomalies")
    zg_anom = zg.groupby("time.dayofyear") - zg_climatology
    logging.info("Write anomalies")
    zg_anom.to_netcdf(d_path + "zg_anom.nc", compute=True)


def compute_detrend(d_path, chunkdict):
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


def compute_norm(d_path, chunkdict):
    zg_detrend = xr.open_dataset(
        d_path + "zg_detrend.nc", chunks=chunkdict)
    logging.info("Normalizing")
    zg_norm = normalize_data(zg_detrend, "zg", window=30)
    logging.info("Writing normalized")
    zg_norm.to_netcdf(d_path + "zg_norm.nc")
