#!/usr/bin/env python3
"""
detect_checkerboard_loess.py

- compute sd(log(z)) over time per grid point
- compute zonal (longitude) spectrum of sd(log(z)) per latitude
- LOESS detrending along longitude (like R span=0.3)
- output:
    - CSV summary (power at requested wavelengths per latitude)
    - heatmap latitude x wavelength (LINEAR colour scale)
    - heatmap lon x lat of temporal sd(log(z))
    - optional per-lat PNGs (variance/spectrum)
"""
import os
import argparse
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.fft import fft
from statsmodels.nonparametric.smoothers_lowess import lowess
import pandas as pd

def ensure_lat_dim(ds):
    """Ensure dataset has 'latitude' and 'longitude' dims."""
    coord_names = list(ds.coords)
    lat_name = next((n for n in coord_names if 'lat' in n.lower()), None)
    lon_name = next((n for n in coord_names if 'lon' in n.lower()), None)
    if lat_name is None or lon_name is None:
        raise ValueError("Could not find latitude/longitude coordinates in dataset")

    # latitude
    if lat_name not in ds.dims:
        lat_vals = np.atleast_1d(ds[lat_name].values)
        ds = ds.expand_dims({'latitude': lat_vals})
        if lat_name != 'latitude':
            ds = ds.rename({lat_name: 'latitude'})
    else:
        if lat_name != 'latitude':
            ds = ds.rename({lat_name: 'latitude'})

    # longitude
    if lon_name not in ds.dims:
        lon_vals = np.atleast_1d(ds[lon_name].values)
        ds = ds.expand_dims({'longitude': lon_vals})
        if lon_name != 'longitude':
            ds = ds.rename({lon_name: 'longitude'})
    else:
        if lon_name != 'longitude':
            ds = ds.rename({lon_name: 'longitude'})
    return ds

def spectrum_along_lon(var_lon, dx=0.5, pad_factor=8):
    """
    Compute FFT spectrum of a 1D zonal array.
    Returns wavelengths (deg) ascending and linear power.
    """
    n = len(var_lon)
    if n < 3 or np.all(np.isnan(var_lon)):
        return np.array([]), np.array([])
    # replace NaNs by linear interpolation
    y = np.array(var_lon, dtype=float)
    nans = np.isnan(y)
    if np.any(nans):
        x = np.arange(n)
        good = ~nans
        y[nans] = np.interp(x[nans], x[good], y[good])
    # zero-pad
    nfft = int(max(4, pad_factor * n))
    ypad = np.zeros(nfft)
    ypad[:n] = y
    Y = fft(ypad)
    P = np.abs(Y)**2 / n
    half = nfft // 2
    pos_idx = np.arange(1, half)
    if pos_idx.size == 0:
        return np.array([]), np.array([])
    freqs = pos_idx / (nfft * dx)  # cycles per degree
    with np.errstate(divide='ignore', invalid='ignore'):
        wavelengths = 1.0 / freqs
    power = P[pos_idx]
    # sort ascending wavelengths
    sort_idx = np.argsort(wavelengths)
    wavelengths_sorted = wavelengths[sort_idx]
    power_sorted = power[sort_idx]
    finite_mask = np.isfinite(wavelengths_sorted)
    return wavelengths_sorted[finite_mask], power_sorted[finite_mask]

def find_nearest_idx(arr, value):
    arr_nonan = np.array(arr)
    if arr_nonan.size == 0 or np.all(np.isnan(arr_nonan)):
        return None
    idx = np.nanargmin(np.abs(arr_nonan - value))
    return int(idx)

def main(infile, outdir, varname, lon_spacing, target_wavelengths,
         plot_per_lat, pad_factor, max_wavelength, wgrid_res, loess_span):
    os.makedirs(outdir, exist_ok=True)
    ds = xr.open_dataset(infile)
    ds = ensure_lat_dim(ds)
    if varname is None:
        varname = list(ds.data_vars)[0]
    da = ds[varname]
    if 'time' not in da.dims:
        raise ValueError("Variable must have a 'time' dimension")

    # sort longitudes
    lons = da['longitude'].values
    lon_sort_idx = np.argsort(lons)
    lons = lons[lon_sort_idx]
    da = da.sel(longitude=lons)

    lats = da['latitude'].values
    nlat, nlon = len(lats), len(lons)

    # --- compute temporal sd of log(z) ---
    eps = max(np.median(da.where(da > 0).values) * 1e-8, 1e-12)
    da_pos = da.where(da > 0, other=eps)
    sd_vals = da_pos.std(dim='time', skipna=True).values  # shape: (lat, lon)

    # --- plot temporal sd heatmap (lon x lat) ---
    plt.figure(figsize=(10,4.5))
    extent = [lons.min(), lons.max(), lats.min(), lats.max()]
    im = plt.imshow(sd_vals, origin='lower', aspect='auto', extent=extent)
    plt.colorbar(im, label='sd(variable)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Temporal sd of variable')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'sd_lon_lat_heatmap.png'), dpi=200)
    plt.close()

    # --- compute FFT per latitude with LOESS detrending ---
    wgrid = np.arange(wgrid_res, max_wavelength + wgrid_res/2, wgrid_res)
    all_power_interp = np.full((nlat, wgrid.size), np.nan)
    summary_rows = []

    for ia, lat in enumerate(lats):
        var_lon = sd_vals[ia, :]
        if np.all(np.isnan(var_lon)):
            continue
        # LOESS detrending
        x = np.arange(nlon)
        y = np.array(var_lon, dtype=float)
        nans = np.isnan(y)
        if np.any(nans):
            good = ~nans
            y[nans] = np.interp(x[nans], x[good], y[good])
        trend = lowess(y, x, frac=loess_span, return_sorted=False)
        y_detr = y - trend

        wavelengths, power = spectrum_along_lon(y_detr, dx=lon_spacing, pad_factor=pad_factor)
        if wavelengths.size == 0:
            continue
        mask_keep = (wavelengths <= max_wavelength) & np.isfinite(wavelengths)
        if not np.any(mask_keep):
            continue
        wl_sorted = wavelengths[mask_keep]
        p_sorted = power[mask_keep]
        # remove duplicates
        _, unique_idx = np.unique(wl_sorted, return_index=True)
        wl_unique = wl_sorted[unique_idx]
        p_unique = p_sorted[unique_idx]
        interp_power = np.interp(wgrid, wl_unique, p_unique, left=np.nan, right=np.nan)
        all_power_interp[ia, :] = interp_power

        # save summary for target wavelengths
        for tw in target_wavelengths:
            idx = find_nearest_idx(wl_unique, tw)
            pval = float(p_unique[idx]) if idx is not None else np.nan
            summary_rows.append({'latitude': float(lat), 'wavelength': float(tw), 'power': pval})

        # optional per-lat PNGs
        if plot_per_lat:
            # raw sd_log plot
            plt.figure(figsize=(9,3))
            plt.plot(lons, var_lon, '-o', ms=2)
            plt.plot(lons, trend, 'r-', lw=1)
            plt.xlabel('Longitude')
            plt.ylabel('sd(log(z))')
            plt.title(f'Latitude {lat:.3f}')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'sd_lat_{ia:03d}_{lat:.3f}.png'), dpi=150)
            plt.close()
            # spectral plot
            plt.figure(figsize=(6,3))
            plt.plot(wl_unique, p_unique, lw=1)
            plt.xlim(0, min(40, np.nanmax(wl_unique)))
            plt.gca().invert_xaxis()
            for tw in target_wavelengths:
                plt.axvline(tw, color='red', alpha=0.5, lw=0.8)
            plt.xlabel('wavelength (deg)')
            plt.ylabel('power')
            plt.title(f'Spectrum (lat={lat:.3f})')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'spec_lat_{ia:03d}_{lat:.3f}.png'), dpi=150)
            plt.close()

    # --- save summary CSV ---
    df = pd.DataFrame(summary_rows)
    if not df.empty:
        df_pivot = df.pivot_table(index='latitude', columns='wavelength', values='power')
        df_pivot.to_csv(os.path.join(outdir, 'power_by_lat_wavelength.csv'))
    else:
        pd.DataFrame().to_csv(os.path.join(outdir, 'power_by_lat_wavelength.csv'))

    # --- heatmap latitude × wavelength (linear scale) ---
    if not np.all(np.isnan(all_power_interp)):
        finite_mask_cols = np.any(np.isfinite(all_power_interp), axis=0)
        Wm = wgrid[finite_mask_cols]
        Pm = all_power_interp[:, finite_mask_cols]
        plt.figure(figsize=(9,6))
        vmin = np.nanpercentile(Pm, 1)
        vmax = np.nanpercentile(Pm, 99)
        im = plt.imshow(Pm, origin='lower', aspect='auto',
                        extent=[Wm.min(), Wm.max(), lats.min(), lats.max()],
                        vmin=vmin, vmax=vmax, cmap='viridis')
        plt.colorbar(im, label='power (linear)')
        plt.xlabel('wavelength (deg longitude)')
        plt.ylabel('latitude')
        plt.title('Latitude × Wavelength spectral power (linear scale)')
        plt.gca().invert_xaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'heatmap_lat_wavelength_linear.png'), dpi=200)
        plt.close()

    print("All outputs written to", outdir)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute zonal spectra of sd(log(z))")
    parser.add_argument('infile', help='input NetCDF file (geopotential.nc)')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--var', help='variable name (default: first data var)', default=None)
    parser.add_argument('--lon_spacing', type=float, default=0.5, help='longitude grid spacing (deg)')
    parser.add_argument('--wavelengths', type=float, nargs='*',
                        default=[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],
                        help='wavelengths to report (deg)')
    parser.add_argument('--plot_per_lat', action='store_true', help='save per-latitude PNGs')
    parser.add_argument('--pad_factor', type=int, default=8, help='zero-padding factor for FFT')
    parser.add_argument('--max_wavelength', type=float, default=10.0, help='max wavelength to consider (deg)')
    parser.add_argument('--wgrid_res', type=float, default=0.02, help='wavelength grid resolution (deg)')
    parser.add_argument('--loess_span', type=float, default=0.3, help='LOESS span fraction')
    args = parser.parse_args()
    main(args.infile, args.outdir, args.var, args.lon_spacing, args.wavelengths,
         args.plot_per_lat, args.pad_factor, args.max_wavelength, args.wgrid_res,
         args.loess_span)
