#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np
import xarray as xr


RKAPPA = 0.286
OMEG2 = 1.4584e-4
DAYSEC = 86400.0
PI180 = np.pi / 180.0
G = 9.80665


def build_file_path(base_dir: str, var: str, level_pa: int, year: int) -> str:
    level_str = f"{level_pa:07d}Pa"
    return os.path.join(
        base_dir,
        var,
        f"r1i1p1-{level_str}",
        f"{var}_day_reanalysis_era5_r1i1p1-{level_str}_{year}0101-{year}1231.nc",
    )


def detect_data_var(ds: xr.Dataset, expected: str) -> str:
    if expected in ds.data_vars:
        return expected
    data_vars = list(ds.data_vars)
    if len(data_vars) == 1:
        return data_vars[0]
    raise ValueError(f"Could not identify data variable '{expected}' in dataset. Found: {data_vars}")


def subset_box(ds: xr.Dataset, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> xr.Dataset:
    lat_name = next((n for n in ds.coords if "lat" in n.lower()), None)
    lon_name = next((n for n in ds.coords if "lon" in n.lower()), None)
    if lat_name is None or lon_name is None:
        raise ValueError("Could not find latitude/longitude coordinates")

    lat = ds[lat_name]
    lon = ds[lon_name]

    # Handle longitude convention
    if float(lon.max()) > 180 and lon_min < 0:
        lon_min_use = lon_min % 360
        lon_max_use = lon_max % 360
    elif float(lon.max()) <= 180 and lon_min >= 0:
        lon_min_use = ((lon_min + 180) % 360) - 180
        lon_max_use = ((lon_max + 180) % 360) - 180
    else:
        lon_min_use = lon_min
        lon_max_use = lon_max

    # Latitude may be descending
    if float(lat[0]) > float(lat[-1]):
        ds = ds.sel({lat_name: slice(lat_max, lat_min)})
    else:
        ds = ds.sel({lat_name: slice(lat_min, lat_max)})

    # Longitude selection, including dateline-crossing case
    if lon_min_use <= lon_max_use:
        ds = ds.sel({lon_name: slice(lon_min_use, lon_max_use)})
    else:
        ds1 = ds.sel({lon_name: slice(lon_min_use, float(lon.max()))})
        ds2 = ds.sel({lon_name: slice(float(lon.min()), lon_max_use)})
        ds = xr.concat([ds1, ds2], dim=lon_name)

    return ds


def open_field(base_dir: str, var: str, level_pa: int, year: int,
               lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> xr.DataArray:
    path = build_file_path(base_dir, var, level_pa, year)
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ds = xr.open_dataset(path)
    ds = subset_box(ds, lon_min, lon_max, lat_min, lat_max)

    data_var = detect_data_var(ds, var)
    da = ds[data_var]

    # Standardize coordinate names
    rename_map = {}
    for c in da.coords:
        cl = c.lower()
        if "lat" in cl and c != "latitude":
            rename_map[c] = "latitude"
        if "lon" in cl and c != "longitude":
            rename_map[c] = "longitude"
        if c == "time":
            continue
    if rename_map:
        da = da.rename(rename_map)

    return da


def to_geopotential_height(z_or_zg: xr.DataArray, varname: str) -> xr.DataArray:
    v = varname.lower()
    units = str(z_or_zg.attrs.get("units", "")).lower()

    if v == "zg":
        return z_or_zg
    if v == "z":
        return z_or_zg / G
    if "m**2 s**-2" in units or "m2 s-2" in units or "m^2 s^-2" in units:
        return z_or_zg / G
    return z_or_zg


def compute_eady(ta_down, ta_up, ua_down, ua_up, va_down, va_up, zg_down, zg_up):
    pkappa_down = ta_down.attrs.get("_level_pa", None)
    pkappa_up = ta_up.attrs.get("_level_pa", None)
    if pkappa_down is None or pkappa_up is None:
        raise ValueError("Pressure level metadata missing on temperature fields")

    p_down = float(pkappa_down)
    p_up = float(pkappa_up)

    pkappa1 = p_down ** RKAPPA
    pkappa2 = p_up ** RKAPPA

    fac1 = (0.31 / np.sqrt(2.0)) * DAYSEC * OMEG2

    lat_rad = ta_down["latitude"] * PI180
    fac = np.abs(fac1 * np.sin(lat_rad))
    fac = xr.DataArray(fac, coords={"latitude": ta_down["latitude"]}, dims=("latitude",))

    qtheta = ta_down * pkappa1 + ta_up * pkappa2
    dtheta = ta_down * pkappa1 - ta_up * pkappa2
    qtheta = qtheta / dtheta

    dz = np.abs(zg_down - zg_up)
    shear2 = (ua_down - ua_up) ** 2 + (va_down - va_up) ** 2

    val2 = qtheta * shear2 / dz
    val2 = xr.where(val2 >= 0, val2, np.nan)

    eady = fac * np.sqrt(val2)
    eady.name = "eady_growth_rate"
    eady.attrs["long_name"] = "Eady growth rate"
    eady.attrs["units"] = "day-1"
    eady.attrs["description"] = "Two-level approximation of Eady growth rate"
    eady.attrs["level_down_pa"] = int(p_down)
    eady.attrs["level_up_pa"] = int(p_up)

    return eady


def main():
    parser = argparse.ArgumentParser(description="Compute two-level Eady growth rate from ERA5 daily fields")
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--base-dir", type=str, required=True,
                        help="Base directory ending at .../ERA5/day/atmos")
    parser.add_argument("--out-dir", type=str, required=True)
    parser.add_argument("--lon-min", type=float, required=True)
    parser.add_argument("--lon-max", type=float, required=True)
    parser.add_argument("--lat-min", type=float, required=True)
    parser.add_argument("--lat-max", type=float, required=True)
    parser.add_argument("--level-down", type=int, default=85000)
    parser.add_argument("--level-up", type=int, default=50000)
    parser.add_argument("--z-var", type=str, default="zg", choices=["zg", "z"],
                        help="ERA5 height variable to use: zg if available, otherwise z")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    year = args.year
    p_down = args.level_down
    p_up = args.level_up

    print(f"Processing year {year}")
    print(f"Levels: {p_down} Pa and {p_up} Pa")

    ta_down = open_field(args.base_dir, "ta", p_down, year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    ta_up   = open_field(args.base_dir, "ta", p_up,   year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    ua_down = open_field(args.base_dir, "ua", p_down, year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    ua_up   = open_field(args.base_dir, "ua", p_up,   year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    va_down = open_field(args.base_dir, "va", p_down, year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    va_up   = open_field(args.base_dir, "va", p_up,   year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    z_down_raw = open_field(args.base_dir, args.z_var, p_down, year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)
    z_up_raw   = open_field(args.base_dir, args.z_var, p_up,   year, args.lon_min, args.lon_max, args.lat_min, args.lat_max)

    ta_down.attrs["_level_pa"] = p_down
    ta_up.attrs["_level_pa"] = p_up

    zg_down = to_geopotential_height(z_down_raw, args.z_var)
    zg_up   = to_geopotential_height(z_up_raw, args.z_var)

    # Align all arrays
    ta_down, ta_up, ua_down, ua_up, va_down, va_up, zg_down, zg_up = xr.align(
        ta_down, ta_up, ua_down, ua_up, va_down, va_up, zg_down, zg_up,
        join="exact"
    )

    eady = compute_eady(ta_down, ta_up, ua_down, ua_up, va_down, va_up, zg_down, zg_up)

    out_name = f"eady_day_era5_{p_down//100:03d}_{p_up//100:03d}_{year}.nc"
    out_path = os.path.join(args.out_dir, out_name)

    encoding = {
        "eady_growth_rate": {
            "zlib": True,
            "complevel": 4,
            "_FillValue": np.float32(np.nan),
            "dtype": "float32",
        }
    }

    ds_out = eady.to_dataset()
    ds_out.attrs["source"] = "ERA5 daily pressure-level fields"
    ds_out.attrs["formula_note"] = "Two-level Eady approximation from ta, ua, va, zg/z"
    ds_out.to_netcdf(out_path, encoding=encoding)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()