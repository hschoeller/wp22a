#!/usr/bin/env python3
import sys
import xarray as xr

# usage:
# python hydrosum.py tcsw.nc tcrw.nc tclw.nc tciw.nc

tcsw, tcrw, tclw, tciw = sys.argv[1:5]

ds_tcsw = xr.open_dataset(tcsw, chunks="auto")
ds_tcrw = xr.open_dataset(tcrw, chunks="auto")
ds_tclw = xr.open_dataset(tclw, chunks="auto")
ds_tciw = xr.open_dataset(tciw, chunks="auto")

hydrosum = (
    ds_tcsw["tcsw"] +
    ds_tcrw["tcrw"] +
    ds_tclw["tclw"] +
    ds_tciw["tciw"]
).rename("hydrosum")

hydrosum.to_dataset().to_netcdf("hydrosum.nc")