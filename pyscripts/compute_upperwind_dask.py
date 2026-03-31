#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np
import xarray as xr
from dask.distributed import Client, LocalCluster


def build_file_path(in_dir: str, var_name: str, level_hpa: int) -> str:
    return str(Path(in_dir) / f"{var_name}{level_hpa}.nc")


def detect_data_var(ds: xr.Dataset, preferred_name: str | None = None) -> str:
    if preferred_name and preferred_name in ds.data_vars:
        return preferred_name

    data_vars = list(ds.data_vars)
    if len(data_vars) == 1:
        return data_vars[0]

    for var in data_vars:
        dims = [d.lower() for d in ds[var].dims]
        has_time = any("time" in d for d in dims)
        has_lat = any("lat" in d for d in dims)
        has_lon = any("lon" in d for d in dims)
        if has_time and has_lat and has_lon:
            return var

    raise ValueError(f"Could not uniquely identify data variable. Found: {data_vars}")


def standardize_da(da: xr.DataArray) -> xr.DataArray:
    """
    Standardize to coordinates/dims named:
      - time
      - latitude
      - longitude

    Handles ERA5-style cases like:
      dims: valid_time
      coords: time, valid_time, latitude, longitude

    without causing rename conflicts.
    """

    # ---- latitude / longitude ----
    lat_candidates = [c for c in da.coords if "lat" in c.lower()]
    lon_candidates = [c for c in da.coords if "lon" in c.lower()]

    if "latitude" not in da.coords and lat_candidates:
        src = lat_candidates[0]
        if src != "latitude":
            da = da.rename({src: "latitude"})

    if "longitude" not in da.coords and lon_candidates:
        src = lon_candidates[0]
        if src != "longitude":
            da = da.rename({src: "longitude"})

    # ---- time ----
    # Case 1: already has a time dimension -> nothing to do
    if "time" in da.dims:
        return da

    # Case 2: has valid_time dimension
    if "valid_time" in da.dims:
        # If there is already a separate coord called "time", drop it first to avoid conflict
        if "time" in da.coords and "time" not in da.dims:
            da = da.drop_vars("time")

        da = da.rename({"valid_time": "time"})
        return da

    # Case 3: some other time-like dimension
    time_dim_candidates = [d for d in da.dims if "time" in d.lower()]
    if time_dim_candidates:
        src = time_dim_candidates[0]
        if src != "time":
            if "time" in da.coords and "time" not in da.dims:
                da = da.drop_vars("time")
            da = da.rename({src: "time"})
        return da

    # Case 4: no time dimension, but maybe time-like coord attached to a dimension
    time_coord_candidates = [c for c in da.coords if "time" in c.lower()]
    for c in time_coord_candidates:
        if c in da.dims:
            if c != "time":
                if "time" in da.coords and "time" not in da.dims:
                    da = da.drop_vars("time")
                da = da.rename({c: "time"})
            return da

    raise ValueError(
        f"Could not identify a usable time dimension. "
        f"dims={da.dims}, coords={list(da.coords)}"
    )


def maybe_drop_singleton_vertical_coord(da: xr.DataArray) -> xr.DataArray:
    for coord in list(da.coords):
        cl = coord.lower()
        if ("level" in cl or "pressure" in cl or cl in {"plev", "lev"}) and da[coord].size == 1:
            da = da.squeeze(coord, drop=True)
    return da


def open_field(in_dir: str, file_var_name: str, level_hpa: int, time_chunk: int) -> xr.DataArray:
    path = build_file_path(in_dir, file_var_name, level_hpa)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    ds = xr.open_dataset(path, chunks={"time": time_chunk})
    data_var = detect_data_var(ds)
    da = ds[data_var]
    da = standardize_da(da)
    da = maybe_drop_singleton_vertical_coord(da)

    if "time" not in da.coords:
        raise ValueError(f"'time' coordinate not found in {path}")

    return da


def compute_upper_wind(ua: xr.DataArray, va: xr.DataArray, level_hpa: int) -> xr.DataArray:
    ua, va = xr.align(ua, va, join="exact")
    wind = np.sqrt(ua**2 + va**2)
    wind.name = "upper_wind_speed"
    wind.attrs["long_name"] = "Upper-level wind speed"
    wind.attrs["units"] = "m s-1"
    wind.attrs["pressure_level_hpa"] = int(level_hpa)
    wind.attrs["description"] = "Wind speed computed as sqrt(u^2 + v^2)"
    return wind


def get_year_range(da: xr.DataArray) -> tuple[int, int]:
    years = np.unique(da["time"].dt.year.values)
    return int(years.min()), int(years.max())


def setup_dask(n_workers: int, threads_per_worker: int) -> Client:
    cluster = LocalCluster(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        processes=True,
        dashboard_address=None,
    )
    return Client(cluster)


def write_single_file(ds_out: xr.Dataset, out_path: str) -> None:
    encoding = {
        "upper_wind_speed": {
            "zlib": True,
            "complevel": 4,
            "_FillValue": np.float32(np.nan),
            "dtype": "float32",
        }
    }
    delayed = ds_out.to_netcdf(out_path, encoding=encoding, compute=False)
    delayed.compute()


def main():
    parser = argparse.ArgumentParser(
        description="Compute all-years upper-level wind speed from single NetCDF files per variable and pressure level using dask"
    )
    parser.add_argument("--in-dir", type=str, required=True)
    parser.add_argument("--out-dir", type=str, required=True)

    parser.add_argument("--level", type=int, default=300)
    parser.add_argument("--u-name", type=str, default="u_component_of_wind")
    parser.add_argument("--v-name", type=str, default="v_component_of_wind")

    parser.add_argument("--n-workers", type=int, default=4)
    parser.add_argument("--threads-per-worker", type=int, default=1)
    parser.add_argument("--time-chunk", type=int, default=180)
    parser.add_argument(
        "--write-mode",
        type=str,
        default="single",
        choices=["single", "yearly"],
        help="Use 'single' for one all-years file. 'yearly' kept only for compatibility.",
    )

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    client = setup_dask(args.n_workers, args.threads_per_worker)
    print(client)

    print("Opening input fields...")
    ua = open_field(args.in_dir, args.u_name, args.level, args.time_chunk)
    va = open_field(args.in_dir, args.v_name, args.level, args.time_chunk)

    y0, y1 = get_year_range(ua)
    print(f"Time range detected: {y0} - {y1}")

    print("Computing upper wind lazily...")
    upper_wind = compute_upper_wind(ua, va, args.level)

    ds_out = upper_wind.to_dataset()
    ds_out.attrs["source"] = "ERA5 single-file pressure-level fields"
    ds_out.attrs["formula_note"] = "Upper-level wind speed from u_component_of_wind and v_component_of_wind"
    ds_out.attrs["time_coverage_start_year"] = y0
    ds_out.attrs["time_coverage_end_year"] = y1
    ds_out.attrs["dask_time_chunk"] = args.time_chunk

    if args.write_mode != "single":
        raise NotImplementedError(
            "This version is intended for all-years output only. "
            "Use --write-mode single."
        )

    out_name = f"upper_wind_day_era5_{args.level:03d}_{y0}_{y1}.nc"
    out_path = os.path.join(args.out_dir, out_name)

    print(f"Writing {out_path}")
    write_single_file(ds_out, out_path)
    print(f"Wrote {out_path}")

    client.close()


if __name__ == "__main__":
    main()