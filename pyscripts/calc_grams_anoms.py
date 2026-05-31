import numpy as np
import xarray as xr
from scipy.ndimage import uniform_filter1d

from dask.distributed import Client, LocalCluster

cluster = LocalCluster(
    n_workers=1,
    threads_per_worker=16,
    memory_limit="56GB",
    processes=False,       # threads only, no subprocess spawning
)
client = Client(cluster)

print(client.dashboard_link)

# ── 1. Load data ──────────────────────────────────────────────────────────────

ds = xr.open_dataset("/scratch/schoelleh96/wp22a/data/geopotential500.nc", 
                     chunks={"time": 1000})
zg500 = ds["z"]  # adjust variable name as needed


# ── Helper functions ──────────────────────────────────────────────────────────

def compute_daily_climatology(
    da: xr.DataArray,
    start_year: int,
    end_year: int,
) -> xr.DataArray:
    """Mean over calendar days for a given reference period."""
    ref_period = da.sel(time=slice(str(start_year), str(end_year)))
    return ref_period.groupby("time.dayofyear").mean("time")


def smooth_climatology_circular(
    da_doy: xr.DataArray,
    window: int,
) -> xr.DataArray:
    """
    Apply a running mean along the dayofyear axis with wrap-around
    padding so that Jan and Dec are treated as neighbours.
    """
    smoothed_values = uniform_filter1d(
        da_doy.values,
        size=window,
        axis=0,
        mode="wrap",
    )
    return da_doy.copy(data=smoothed_values)


def subtract_smoothed_climatology(
    da: xr.DataArray,
    smoothed_clim: xr.DataArray,
) -> xr.DataArray:
    """Compute anomalies by subtracting the smoothed climatology."""
    return da.groupby("time.dayofyear") - smoothed_clim


def build_lanczos_weights(
    cutoff_freq: float,
    n_weights: int,
) -> np.ndarray:
    """
    Lanczos low-pass filter weights.

    Uses numpy's normalised sinc (sin(πx)/πx), so the ideal
    low-pass kernel  sin(2π f₀ k)/(π k) = 2f₀ · sinc(2f₀ k).
    The Lanczos window is  sinc(k / half_width).
    """
    half = (n_weights - 1) // 2
    k = np.arange(-half, half + 1)
    ideal_lowpass = 2.0 * cutoff_freq * np.sinc(2.0 * cutoff_freq * k)
    lanczos_window = np.sinc(k / half)
    weights = ideal_lowpass * lanczos_window
    return weights / weights.sum()


def apply_lanczos_lowpass(
    da: xr.DataArray,
    cutoff_days: float,
    n_weights: int,
) -> xr.DataArray:
    """
    Apply a centred Lanczos low-pass filter along the time axis.
    Edge values within half the window width are returned as NaN.
    """
    weights = build_lanczos_weights(1.0 / cutoff_days, n_weights)
    weight_da = xr.DataArray(weights, dims=["window"])
    return (
        da.rolling(time=n_weights, center=True)
        .construct("window")
        .dot(weight_da)
    )


# ── 2. Processing pipeline ────────────────────────────────────────────────────

# Step 1 – calendar-day climatology (1979-2019)
zg500_m = compute_daily_climatology(zg500, 1979, 2019)

# Step 2 – 90-day circular running mean of the climatology
zg500_m90 = smooth_climatology_circular(zg500_m, window=90)

# Step 3 – anomalies relative to the smoothed climatology
zg500_prime = subtract_smoothed_climatology(zg500, zg500_m90)

# Step 4 – Lanczos low-pass filter (10-day cutoff, N=21)
zg500_prime_lp = apply_lanczos_lowpass(
    zg500_prime,
    cutoff_days=10,
    n_weights=21,
)


# ── 3. Persist results ────────────────────────────────────────────────────────

out = xr.Dataset(
    {
        "zg500_m":        zg500_m,
        "zg500_m90":      zg500_m90,
        "zg500_prime":    zg500_prime,
        "zg500_prime_lp": zg500_prime_lp,
    }
)
out.to_netcdf("zg500_processed.nc")