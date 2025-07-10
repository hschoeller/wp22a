import numpy as np
import xarray as xr
from scipy.signal import detrend
import os
import cfgrib
import logging
import cdsapi
import zipfile
import dask


def retrieve_era5(
    years, lat_min, lat_max, lon_min, lon_max, d_path, f_name, var,
    dataset="reanalysis-era5-single-levels", product_type="reanalysis", pressure_level=None
):
    """
    Retrieve ERA5 data from CDS API.

    Parameters:
    -----------
    years : list or str
        Years to retrieve.
    lat_min, lat_max, lon_min, lon_max : float
        Bounding box for area.
    d_path : str
        Directory to save file.
    f_name : str
        Name of the zip file to download.
    var : str
        Variable to retrieve.
    dataset : str
        CDS dataset name.
    product_type : str
        Product type for request.
    pressure_level : str or None
        Pressure level (e.g., "500") if required.

    Returns:
    --------
    str
        Name of the extracted GRIB file.
    """
    c = cdsapi.Client()

    request = {
        "product_type": [product_type],
        "variable": [var],
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
        "data_format": "grib",
        "download_format": "zip",
        "area": [lat_max, lon_min, lat_min, lon_max]
    }

    if pressure_level is not None:
        request["pressure_level"] = [str(pressure_level)]

    c.retrieve(dataset, request, d_path + f_name)

    # Unzip the downloaded file
    with zipfile.ZipFile(f"{d_path}/{f_name}", 'r') as zip_ref:
        zip_ref.extract("data.grib", path=d_path, pwd=None)
        extracted_file = os.path.join(d_path, "data.grib")
        target_file = os.path.join(d_path, f"{var}.grib")
        if os.path.exists(target_file):
            os.remove(target_file)
        os.rename(extracted_file, target_file)
    os.remove(f"{d_path}/{f_name}")
    return f"{var}.grib"


def convert_grib_to_nc(d_path, grib_name, cleanup=True):
    """
    Convert a GRIB file to NetCDF format, subsetting to 0.5° resolution grid points.

    Parameters:
    -----------
    d_path : str
        Directory path of the input GRIB file
    grib_name : str
        Name of the input GRIB file
    cleanup : bool
        Whether to delete the original GRIB file after conversion

    Returns:
    --------
    bool
        True if conversion was successful, False otherwise
    """
    try:
        print(f"opening file {os.path.join(d_path, grib_name)}")
        # Open the GRIB file using cfgrib and xarray

        ds = xr.open_dataset(os.path.join(d_path, grib_name), engine='cfgrib')
        print(ds)

        lat_name = [dim for dim in ds.dims if 'lat' in dim][0]
        lon_name = [dim for dim in ds.dims if 'lon' in dim][0]
        print(f"Latitude name: {lat_name}, Longitude name: {lon_name}")

        def is_half_degree(arr):
            return np.isclose(np.mod(arr, 1), 0.0) | np.isclose(np.mod(arr, 1), 0.5)

        # Create boolean masks for each coordinate
        lat_mask = is_half_degree(ds[lat_name])
        lon_mask = is_half_degree(ds[lon_name])

        # Method 1: Using isel with boolean indexing
        ds = ds.isel({
            lat_name: lat_mask,
            lon_name: lon_mask
        })

        # Save as NetCDF
        base_name = os.path.splitext(grib_name)[0]
        nc_filename = f"{base_name}.nc"
        nc_filepath = os.path.join(d_path, nc_filename)
        print(f"Saving to NetCDF file: {nc_filepath}")
        ds.to_netcdf(nc_filepath)

        # Close and optionally clean up
        ds.close()
        print(f"Successfully converted {grib_name} to {nc_filepath}")

        if cleanup:
            os.remove(os.path.join(d_path, grib_name))

        return True

    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        return False


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
                 d_path, mf=True, member=None):
    if mf:
        ds = xr.open_mfdataset(files, parallel=True,
                               combine="by_coords",
                               chunks=chunkdict)
    else:
        ds = xr.open_dataset(files, engine="cfgrib",
                             chunks=chunkdict).isel(number=member)
        ds = ds.rename({'z': 'zg', 'longitude': 'lon',
                        'latitude': 'lat'})
    ds = ds.assign_coords(
        lon=(((ds.lon + 180) % 360) - 180)
    )

    # Sort the dataset by longitude for proper slicing
    ds = ds.sortby('lon')
    ds = ds.sortby('lat')
    zg = ds['zg'].sel(lat=slice(lat_min, lat_max),
                      lon=slice(lon_min, lon_max)).chunk(chunkdict)
    logging.info(zg.chunks)
    if member == 0:
        zg_climatology = zg.groupby("time.dayofyear").mean(dim="time")
        zg_climatology = zg_climatology.persist()
        logging.info(zg_climatology.chunks)
        logging.info("saving to disk")
        zg_climatology.to_netcdf(d_path + "zg_climatology.nc")
        logging.info("Climatology computation completed.")
    return zg


def compute_anom(zg, chunkdict, d_path, member=None):
    zg_climatology = xr.open_dataset(
        d_path + "zg_climatology.nc", chunks=chunkdict)
    logging.info("Calculate anomalies")
    zg_anom = zg.groupby("time.dayofyear") - zg_climatology
    logging.info("Write anomalies")
    zg_anom.to_netcdf(d_path + f"zg_anom{member}.nc", compute=True)


def compute_detrend(d_path, chunkdict, member=None):
    zg_anom = xr.open_dataset(d_path + f"zg_anom{member}.nc",
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
    zg_detrend.to_netcdf(d_path + f"zg_detrend{member}.nc")


def compute_norm(d_path, chunkdict, member=None):
    zg_detrend = xr.open_dataset(
        d_path + f"zg_detrend{member}.nc", chunks=chunkdict)
    logging.info("Normalizing")
    zg_norm = normalize_data(zg_detrend, "zg", window=30)
    logging.info("Writing normalized")
    zg_norm.to_netcdf(d_path + f"zg_norm{member}.nc")


def calculate_ensemble_spread(input_file, output_file, year, chunkdict):
    """
    Calculate the ensemble var for the specified year and save as .nc file.

    Args:
        input_file (str): Path to the input GRIB file.
        output_file (str): Path to save the resulting NetCDF file.
        year (int): Year for which to calculate the ensemble spread.
    """
    # Load data using xarray and cfgrib
    ds = xr.open_dataset(input_file, engine="cfgrib", chunks=chunkdict)

    # Filter data for the specified year
    yearly_data = ds.sel(time=ds.time.dt.year == year).where(
        ds["number"] != 0, drop=True)

    # Compute ensemble spread (var across ensemble dimension)
    spread = yearly_data.var(dim="number")

    # Save the result as a NetCDF file
    spread.to_netcdf(output_file)
    print(f"Saved ensemble spread for {year} to {output_file}")
