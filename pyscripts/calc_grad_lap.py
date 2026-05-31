#!/usr/bin/env python3
"""
Compute horizontal gradient magnitude, spherical Laplacian,
and absolute spherical Laplacian for a scalar field on a regular
latitude-longitude grid.

FIXED VERSION:
- Output grid is guaranteed to lie exactly on target spacing (e.g. 0.5°)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import dask.array as dsa
import numpy as np
import xarray as xr


EARTH_RADIUS_M = 6_371_000.0
GRAVITY = 9.80665

DEFAULT_SOURCE_DEG = 0.25
DEFAULT_TARGET_DEG = 0.5
DEFAULT_RES_TOL_DEG = 0.02
DEFAULT_MAX_TIME_CHUNK = 50
DEFAULT_TARGET_MB_PER_WORKER = 256


def log(msg: str) -> None:
    print(f"[INFO] {msg}", flush=True)


# =========================
# ARGUMENTS
# =========================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--grad-output", required=True)
    parser.add_argument("--lap-output", required=True)
    parser.add_argument("--abs-lap-output", required=True)
    parser.add_argument("--var", required=True)

    parser.add_argument("--lat", default="latitude")
    parser.add_argument("--lon", default="longitude")
    parser.add_argument("--time", default="time")

    parser.add_argument("--convert-geopotential-to-height", action="store_true")

    parser.add_argument("--engine", default=None)

    parser.add_argument("--auto-coarsen-from-deg", type=float, default=DEFAULT_SOURCE_DEG)
    parser.add_argument("--coarsen-to-deg", type=float, default=DEFAULT_TARGET_DEG)
    parser.add_argument("--resolution-tol-deg", type=float, default=DEFAULT_RES_TOL_DEG)

    parser.add_argument("--edge-order", type=int, default=2, choices=(1, 2))

    parser.add_argument("--max-time-chunk", type=int, default=DEFAULT_MAX_TIME_CHUNK)
    parser.add_argument("--target-mb-per-worker", type=int, default=DEFAULT_TARGET_MB_PER_WORKER)

    return parser.parse_args()


# =========================
# GRID HELPERS
# =========================
def _regular_spacing_deg(coord: xr.DataArray) -> float:
    vals = np.asarray(coord.values, dtype=np.float64)
    diffs = np.abs(np.diff(vals))
    spacing = float(np.median(diffs))
    if not np.allclose(diffs, spacing, atol=1e-8):
        raise ValueError("Coordinate not regular")
    return spacing


def detect_grid_spacing_deg(da: xr.DataArray, lat: str, lon: str):
    return _regular_spacing_deg(da[lat]), _regular_spacing_deg(da[lon])


def regular_target_coord_1d(source_vals, target_deg, n_out):
    """
    Construct EXACT regular grid (no inherited offsets).
    """
    source_vals = np.asarray(source_vals, dtype=np.float64)

    ascending = source_vals[-1] > source_vals[0]
    start = np.round(source_vals[0] / target_deg) * target_deg

    if ascending:
        out = start + target_deg * np.arange(n_out)
    else:
        out = start - target_deg * np.arange(n_out)

    return np.round(out, 10)


# =========================
# COARSENING
# =========================
def cosine_weighted_coarsen_numpy(field, lat_vals, lon_vals, target_deg):
    nt, ny, nx = field.shape

    ny2 = ny // 2
    nx2 = nx // 2

    field = field[:, : ny2 * 2, : nx2 * 2]
    lat_vals = lat_vals[: ny2 * 2]
    lon_vals = lon_vals[: nx2 * 2]

    block = field.reshape(nt, ny2, 2, nx2, 2)

    lat_w = np.cos(np.deg2rad(lat_vals)).reshape(ny2, 2)

    lon_mean = block.mean(axis=-1)
    numerator = (lon_mean * lat_w[None, :, :, None]).sum(axis=2)
    denominator = lat_w.sum(axis=1)[None, :, None]

    out = numerator / denominator

    lat_out = regular_target_coord_1d(lat_vals, target_deg, ny2)
    lon_out = regular_target_coord_1d(lon_vals, target_deg, nx2)

    return out, lat_out, lon_out


# =========================
# CORE COMPUTATION
# =========================
def compute_block_diagnostics(
    block,
    lat_name,
    lon_name,
    time_name,
    source_deg,
    target_deg,
    tol_deg,
    radius,
    edge_order,
):
    field = np.asarray(block.data, dtype=np.float64)

    lat_vals = np.asarray(block[lat_name].values)
    lon_vals = np.asarray(block[lon_name].values)

    lat_step, lon_step = detect_grid_spacing_deg(block, lat_name, lon_name)

    if abs(lat_step - source_deg) <= tol_deg:
        field, lat_vals, lon_vals = cosine_weighted_coarsen_numpy(
            field, lat_vals, lon_vals, target_deg
        )
    elif abs(lat_step - target_deg) <= tol_deg:
        lat_vals = regular_target_coord_1d(lat_vals, target_deg, len(lat_vals))
        lon_vals = regular_target_coord_1d(lon_vals, target_deg, len(lon_vals))
    else:
        raise ValueError("Unsupported grid resolution")

    lat_rad = np.deg2rad(lat_vals)
    lon_rad = np.deg2rad(lon_vals)

    coslat = np.cos(lat_rad)[None, :, None]

    dfdphi = np.gradient(field, lat_rad, axis=1, edge_order=edge_order)
    dfdlambda = np.gradient(field, lon_rad, axis=2, edge_order=edge_order)

    dfdx = dfdlambda / (radius * coslat)
    dfdy = dfdphi / radius

    grad = np.hypot(dfdx, dfdy)

    term_phi = np.gradient(coslat * dfdphi, lat_rad, axis=1) / (radius**2 * coslat)
    term_lambda = np.gradient(dfdlambda, lon_rad, axis=2) / (radius**2 * coslat**2)

    lap = term_phi + term_lambda
    abs_lap = np.abs(lap)

    coords = {
        time_name: block[time_name],
        lat_name: lat_vals,
        lon_name: lon_vals,
    }

    return xr.Dataset(
        {
            f"{block.name}_grad_mag": ((time_name, lat_name, lon_name), grad),
            f"{block.name}_laplacian": ((time_name, lat_name, lon_name), lap),
            f"{block.name}_abs_laplacian": ((time_name, lat_name, lon_name), abs_lap),
        },
        coords=coords,
    )


# =========================
# TEMPLATE
# =========================
def build_template(da, lat, lon, time, source_deg, target_deg, tol_deg):
    lat_step, lon_step = detect_grid_spacing_deg(da, lat, lon)

    if abs(lat_step - source_deg) <= tol_deg:
        nlat = da.sizes[lat] // 2
        nlon = da.sizes[lon] // 2
    else:
        nlat = da.sizes[lat]
        nlon = da.sizes[lon]

    lat_out = regular_target_coord_1d(da[lat].values, target_deg, nlat)
    lon_out = regular_target_coord_1d(da[lon].values, target_deg, nlon)

    # ✅ CRITICAL: use SAME chunking as input
    time_chunks = da.chunksizes[time]

    chunks = (time_chunks, nlat, nlon)
    shape = (da.sizes[time], nlat, nlon)

    grad_name = f"{da.name}_grad_mag"
    lap_name = f"{da.name}_laplacian"
    abs_lap_name = f"{da.name}_abs_laplacian"

    template = xr.Dataset(
        {
            grad_name: (
                (time, lat, lon),
                dsa.empty(shape, chunks=chunks, dtype=np.float64),
            ),
            lap_name: (
                (time, lat, lon),
                dsa.empty(shape, chunks=chunks, dtype=np.float64),
            ),
            abs_lap_name: (
                (time, lat, lon),
                dsa.empty(shape, chunks=chunks, dtype=np.float64),
            ),
        },
        coords={
            time: da[time],
            lat: xr.DataArray(lat_out, dims=(lat,)),
            lon: xr.DataArray(lon_out, dims=(lon,)),
        },
    )

    return template


# =========================
# MAIN
# =========================
def main():
    args = parse_args()

    ds = xr.open_dataset(args.input, chunks={args.time: 1})

    da = ds[args.var]

    if args.convert_geopotential_to_height:
        da = da / GRAVITY

    da = da.chunk({args.time: -1, args.lat: -1, args.lon: -1})

    template = build_template(
        da,
        args.lat,
        args.lon,
        args.time,
        args.auto_coarsen_from_deg,
        args.coarsen_to_deg,
        args.resolution_tol_deg,
    )

    out = xr.map_blocks(
        compute_block_diagnostics,
        da,
        kwargs=dict(
            lat_name=args.lat,
            lon_name=args.lon,
            time_name=args.time,
            source_deg=args.auto_coarsen_from_deg,
            target_deg=args.coarsen_to_deg,
            tol_deg=args.resolution_tol_deg,
            radius=EARTH_RADIUS_M,
            edge_order=args.edge_order,
        ),
        template=template,
    )

    out.to_netcdf(args.grad_output)


if __name__ == "__main__":
    sys.exit(main())