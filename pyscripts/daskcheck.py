#!/usr/bin/env python3

from pathlib import Path

import xarray as xr
import dask.array as da

path = Path("u_component_of_wind500.nc")

print(f"Testing file: {path.resolve()}")
print()

# Open lazily with dask chunks
ds = xr.open_dataset(path, chunks={"time": 100})

print("Dataset summary:")
print(ds)
print()

print("Coordinates:")
for name, coord in ds.coords.items():
    print(f"  {name}: shape={coord.shape}, dtype={coord.dtype}")
print()

print("Data variables:")
for name, var in ds.data_vars.items():
    backend = type(var.data)
    is_dask = isinstance(var.data, da.Array)
    print(f"  {name}: shape={var.shape}, dims={var.dims}, dtype={var.dtype}")
    print(f"     backend={backend}")
    print(f"     dask-backed={is_dask}")
    if is_dask:
        print(f"     chunks={var.chunks}")
print()

# Pick main variable
if len(ds.data_vars) == 1:
    vname = list(ds.data_vars)[0]
else:
    candidates = []
    for v in ds.data_vars:
        dims_lower = [d.lower() for d in ds[v].dims]
        if any("time" in d for d in dims_lower) and any("lat" in d for d in dims_lower) and any("lon" in d for d in dims_lower):
            candidates.append(v)
    if len(candidates) == 1:
        vname = candidates[0]
    else:
        raise ValueError(f"Could not uniquely detect main variable. Found: {list(ds.data_vars)}")

da_u = ds[vname]
print(f"Selected variable: {vname}")
print()

# Check time coordinate
time_name = next((c for c in ds.coords if "time" in c.lower()), None)
if time_name is None:
    raise ValueError("No time coordinate found")

time = ds[time_name]
print(f"Time coordinate name: {time_name}")
print(f"First 3 time values: {time.values[:3]}")
print(f"Last 3 time values:  {time.values[-3:]}")
print()

# Check for lat/lon coords
lat_name = next((c for c in ds.coords if "lat" in c.lower()), None)
lon_name = next((c for c in ds.coords if "lon" in c.lower()), None)
print(f"Latitude coord:  {lat_name}")
print(f"Longitude coord: {lon_name}")
print()

# Tiny dask computation
print("Running a tiny lazy computation...")
sample = da_u.isel(time=slice(0, 10)).mean()
print("Before compute:", sample)
result = sample.compute()
print("After compute:", result.values)
print()

print("Variable attributes:")
for k, v in da_u.attrs.items():
    print(f"  {k}: {v}")