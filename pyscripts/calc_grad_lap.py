#!/usr/bin/env python3
"""
Compute horizontal gradient magnitude and spherical Laplacian
for a scalar field on a regular latitude-longitude grid.

Designed for z500 / geopotential height on daily data, but can be used
for any scalar variable on a regular lat-lon grid.

Example:
    python compute_z500_diagnostics.py \
        --input z500_daily_1940_2024.nc \
        --output z500_diagnostics.nc \
        --var z500 \
        --lat latitude \
        --lon longitude \
        --time time

Notes
-----
- Assumes a regular lat-lon grid.
- Uses spherical geometry:
      d/dx = 1 / (R cos(phi)) * d/dlambda
      d/dy = 1 / R * d/dphi
- Computes the spherical Laplacian:
      ∇²f = 1/(R² cosφ) ∂/∂φ (cosφ ∂f/∂φ) + 1/(R² cos²φ) ∂²f/∂λ²
- Longitude is not treated as explicitly periodic here; interior points are
  generally fine, but the seam can be slightly less accurate. If you want a
  cyclic-longitude version later, that can be done too.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr


EARTH_RADIUS_M = 6_371_000.0
GRAVITY = 9.80665


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute gradient magnitude and Laplacian on a regular lat-lon grid."
    )

    parser.add_argument("--input", required=True, help="Input NetCDF file")
    parser.add_argument("--grad-output", required=True, help="Output NetCDF file for gradient magnitude")
    parser.add_argument("--lap-output", required=True, help="Output NetCDF file for Laplacian")
    parser.add_argument("--var", required=True, help="Variable name, e.g. z500 or z")
    parser.add_argument("--lat", default="latitude", help="Latitude dimension/coord name")
    parser.add_argument("--lon", default="longitude", help="Longitude dimension/coord name")
    parser.add_argument("--time", default="time", help="Time dimension name")
    parser.add_argument(
        "--convert-geopotential-to-height",
        action="store_true",
        help="Convert geopotential (m^2 s^-2) to geopotential height (m) by dividing by g",
    )
    parser.add_argument(
        "--chunks",
        nargs="*",
        default=None,
        help=(
            "Optional dask chunking as key=value pairs, e.g. "
            "--chunks time=100 latitude=181 longitude=360"
        ),
    )
    parser.add_argument(
        "--engine",
        default=None,
        help="Optional xarray engine, e.g. netcdf4 or h5netcdf"
    )

    return parser.parse_args()


def parse_chunks(chunk_args: list[str] | None) -> dict[str, int] | None:
    if not chunk_args:
        return None

    chunks: dict[str, int] = {}
    for item in chunk_args:
        if "=" not in item:
            raise ValueError(f"Invalid chunk specification: {item!r}. Use key=value.")
        key, value = item.split("=", 1)
        chunks[key] = int(value)
    return chunks


def validate_grid(da: xr.DataArray, lat_name: str, lon_name: str) -> None:
    if lat_name not in da.coords:
        raise ValueError(f"Latitude coordinate {lat_name!r} not found in data array.")
    if lon_name not in da.coords:
        raise ValueError(f"Longitude coordinate {lon_name!r} not found in data array.")

    if da[lat_name].ndim != 1 or da[lon_name].ndim != 1:
        raise ValueError("This script only supports 1D latitude and longitude coordinates.")

    lat_vals = da[lat_name].values
    lon_vals = da[lon_name].values

    if lat_vals.size < 3 or lon_vals.size < 3:
        raise ValueError("Latitude and longitude must each have at least 3 points.")

    dlat = np.diff(lat_vals)
    dlon = np.diff(lon_vals)

    if not np.allclose(dlat, dlat[0], rtol=0, atol=1e-8):
        raise ValueError("Latitude grid is not regular.")
    if not np.allclose(dlon, dlon[0], rtol=0, atol=1e-8):
        raise ValueError("Longitude grid is not regular.")


def compute_diagnostics(
    da: xr.DataArray,
    lat_name: str,
    lon_name: str,
    radius: float = EARTH_RADIUS_M,
) -> xr.Dataset:
    """
    Compute first derivatives, gradient magnitude, and spherical Laplacian.
    """
    da = da.astype(np.float64)

    validate_grid(da, lat_name, lon_name)

    original_units = da.attrs.get("units", "")

    # Convert coordinates to radians for differentiation.
    lat_rad = np.deg2rad(da[lat_name])
    lon_rad = np.deg2rad(da[lon_name])

    f = da.assign_coords({
        lat_name: lat_rad,
        lon_name: lon_rad,
    })

    coslat = np.cos(f[lat_name])
    eps = 1e-12
    coslat_safe = xr.where(np.abs(coslat) < eps, np.nan, coslat)

    # First derivatives in spherical coordinates
    dfdphi = f.differentiate(coord=lat_name)       # ∂f/∂φ
    dfdlambda = f.differentiate(coord=lon_name)    # ∂f/∂λ

    # Convert to physical eastward/northward derivatives
    dfdx = dfdlambda / (radius * coslat_safe)
    dfdy = dfdphi / radius

    grad_mag = np.hypot(dfdx, dfdy)

    # Spherical Laplacian
    term_phi = (coslat * dfdphi).differentiate(coord=lat_name) / (radius**2 * coslat_safe)
    d2fdlambda2 = dfdlambda.differentiate(coord=lon_name)
    term_lambda = d2fdlambda2 / (radius**2 * coslat_safe**2)
    laplacian = term_phi + term_lambda

    var_name = da.name if da.name is not None else "field"

    dfdx = dfdx.rename(f"{var_name}_dx")
    dfdy = dfdy.rename(f"{var_name}_dy")
    grad_mag = grad_mag.rename(f"{var_name}_grad_mag")
    laplacian = laplacian.rename(f"{var_name}_laplacian")
    abs_laplacian = np.abs(laplacian).rename(f"{var_name}_abs_laplacian")

    dfdx.attrs.update({
        "long_name": f"Eastward derivative of {var_name}",
        "units": f"{original_units} m^-1".strip(),
    })
    dfdy.attrs.update({
        "long_name": f"Northward derivative of {var_name}",
        "units": f"{original_units} m^-1".strip(),
    })
    grad_mag.attrs.update({
        "long_name": f"Horizontal gradient magnitude of {var_name}",
        "units": f"{original_units} m^-1".strip(),
    })
    laplacian.attrs.update({
        "long_name": f"Spherical Laplacian of {var_name}",
        "units": f"{original_units} m^-2".strip(),
    })
    abs_laplacian.attrs.update({
        "long_name": f"Absolute spherical Laplacian of {var_name}",
        "units": f"{original_units} m^-2".strip(),
    })

    out = xr.Dataset(
        data_vars={
            grad_mag.name: grad_mag,
            laplacian.name: laplacian,
        }
    )

    # Restore original degree coordinates
    out = out.assign_coords({
        lat_name: da[lat_name],
        lon_name: da[lon_name],
    })

    out.attrs["description"] = (
        "Diagnostics computed on a regular lat-lon grid using spherical geometry"
    )
    out.attrs["earth_radius_m"] = radius

    return out


def main() -> int:
    args = parse_args()

    input_path = Path(args.input)
    grad_output_path = Path(args.grad_output)
    lap_output_path = Path(args.lap_output)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    chunks = parse_chunks(args.chunks)

    open_kwargs = {}
    if args.engine is not None:
        open_kwargs["engine"] = args.engine
    if chunks is not None:
        open_kwargs["chunks"] = chunks

    print(f"Opening {input_path}")
    ds = xr.open_dataset(input_path, **open_kwargs)

    if args.var not in ds:
        raise KeyError(f"Variable {args.var!r} not found in dataset. Variables: {list(ds.data_vars)}")

    da = ds[args.var]

    if args.convert_geopotential_to_height:
        print("Converting geopotential to geopotential height by dividing by g.")
        da = da / GRAVITY
        da.attrs = dict(da.attrs)
        da.attrs["units"] = "m"
        if "long_name" in da.attrs:
            da.attrs["long_name"] = f"{da.attrs['long_name']} converted to geopotential height"

    print("Computing diagnostics...")
    out = compute_diagnostics(
        da=da,
        lat_name=args.lat,
        lon_name=args.lon,
        radius=EARTH_RADIUS_M,
    )

    # Compression settings for NetCDF output
    encoding = {}
    for var in out.data_vars:
        encoding[var] = {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
        }

    grad_var = f"{da.name if da.name is not None else 'field'}_grad_mag"
    lap_var = f"{da.name if da.name is not None else 'field'}_laplacian"

    print(f"Writing {grad_output_path}")
    out[[grad_var]].to_netcdf(grad_output_path, encoding={grad_var: encoding[grad_var]})

    print(f"Writing {lap_output_path}")
    out[[lap_var]].to_netcdf(lap_output_path, encoding={lap_var: encoding[lap_var]})

    print("Done.")

    return 0


if __name__ == "__main__":
    sys.exit(main())