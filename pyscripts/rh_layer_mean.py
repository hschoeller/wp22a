#!/usr/bin/env python3
import sys
import xarray as xr

# usage:
# python rh_layer_mean.py rh700.nc rh750.nc rh800.nc rh850.nc

files = sys.argv[1:]
datasets = [xr.open_dataset(f, chunks="auto") for f in files]

# assumes each file has one pressure level and variable name "r"
arrs = [ds["r"] for ds in datasets]

rh_mean = xr.concat(arrs, dim="level").mean("level").rename("r")
rh_mean.to_dataset().to_netcdf("relative_humidity_700-850.nc")