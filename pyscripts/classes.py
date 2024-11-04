#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:44:41 2024

@author: schoelleh96
"""

# class and function definitions
from scipy import stats
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import cmcrameri.cm as cmc  # Import Fabio Crameri's colormap package
import pandas as pd


class MidpointNormalize(mcolors.Normalize):
    """Custom Norm to allow for midpoint at 0 with constant gradient."""

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax,   clip)

    def __call__(self, value, clip=None):
        """Return the value in the given norm."""
        if -self.vmin > self.vmax:
            x, y = [self.vmin, self.midpoint, -self.vmin], [0, 0.5, 1]
        else:
            x, y = [-self.vmax, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


class DataHandler:
    """Class to load a NetCDF file and plot."""

    def __init__(self, filepath):
        """
        Initialize the DataHandler with the path to the NetCDF file.

        Parameters
        ----------
        filepath (str): Path to the NetCDF file containing data.
        """
        self.filepath = filepath
        self.dataset = None
        self.trends = None
        self.p_values = None

    def load_data(self):
        """Load the data from file. Initialize trend and p_value datasets."""
        self.dataset = xr.open_dataset(self.filepath)
        self.trends = xr.Dataset()
        self.p_values = xr.Dataset()

    def calc_trend(self, var):
        """Calculate gridpoint-wise trend in the data.

        Parameters
        ----------
        var (str): The variable name.
        """
        # Define a function to apply linear regression for each time series
        def linregress_func(time_series):
            # Assume the time dimension is the first dimension (axis=0)
            n_time = time_series.shape[0]
            # Create a time array (assuming evenly spaced time points)
            time = np.arange(n_time)
            # Apply linregress to each time series (ignoring NaNs)
            mask = ~np.isnan(time_series)
            if np.sum(mask) < 2:
                return np.nan, np.nan  # Insufficient data points
            result = stats.linregress(time[mask], time_series[mask])
            # Return slope and p-value
            return result.slope, result.pvalue
        # Use xr.apply_ufunc to apply the function over the spatial dimensions
        self.trends[var], self.p_values[var] = xr.apply_ufunc(
            linregress_func,
            self.dataset[var],  # Input data
            input_core_dims=[['time']],  # Apply over the 'time' dimension
            output_core_dims=[[], []],  # Output two scalars: slope, p-value
            vectorize=True,  # Vectorize the function for all grid points
            dask="parallelized",  # Enable parallelization
            output_dtypes=[float, float]  # Output types
        )

    def plot(self, data, clabel, title=None, stipple_mask=None):
        """
        Plot data.

        Parameters
        ----------
        data (xarray.DataArray): The data.
        clabel (str): Colorbar label.
        title (str, optional): Title of the plot.
        stipple_mask (xarray.DataArray, optional): Mask for stippling.
        """
        # Extract latitude and longitude
        lons = self.dataset['lon']
        lats = self.dataset['lat']

        # Set up the plot with Cartopy
        fig, ax = plt.subplots(subplot_kw={
            'projection': ccrs.Orthographic(
                central_longitude=-45, central_latitude=45)})
        ax.coastlines()  # Add coastlines
        ax.add_feature(cfeature.BORDERS)  # Add country borders
        ax.gridlines(draw_labels=True)  # Add lat/lon grid lines

        # Plot the data
        norm = MidpointNormalize(vmin=np.min(data), vmax=np.max(data),
                                 midpoint=0)
        plot = ax.contourf(lons, lats, data, 60, transform=ccrs.PlateCarree(),
                           cmap=cmc.vik, norm=norm)
        # Add stippling if stipple_mask is provided
        if stipple_mask is not None:
            stipple_mask = np.asarray(stipple_mask, dtype=bool)
            # Create a masked array for the hatching regions
            hatching_data = np.ma.masked_where(~stipple_mask,
                                               np.ones_like(data))
            # Overlay hatching
            ax.contourf(lons, lats, hatching_data, levels=[0, 1],
                        hatches=['.'], transform=ccrs.PlateCarree(), alpha=0)

        # Add a color bar
        cbar = plt.colorbar(plot, ax=ax, orientation='horizontal', pad=0.05,
                            fraction=.046)
        cbar.set_label(clabel)
        # Set fewer ticks using MaxNLocator
        cbar.locator = MaxNLocator(nbins=5)  # Reduce to 5 ticks
        cbar.update_ticks()

        # Set tick labels to scientific notation
        cbar.formatter = ScalarFormatter(useMathText=True)
        cbar.formatter.set_powerlimits((-2, 2))
        cbar.update_ticks()

        # Set the title
        if title is not None:
            ax.set_title(title, fontsize=14)

        plt.show()
        return fig, ax

    def plot_data(self, var, t_i, clabel, title=None):
        """
        Plot data from self.dataset at a specific time index.

        Parameters
        ----------
        var (str): Variable name.
        t_i (int): Time index.
        clabel (str): Colorbar label.
        title (str, optional): Plot title.
        """
        if self.dataset is None:
            raise ValueError("Dataset not loaded. Call load_data() first.")

        # Select the data for the variable and time index
        data = self.dataset[var].isel(time=t_i)

        # Call the abstract plot function
        self.plot(data, clabel, title)

    def plot_trend(self, var, clabel, title=None):
        """
        Plot the trend data with stippling for significance.

        Parameters
        ----------
        var (str): Variable name.
        clabel (str): Colorbar label.
        title (str, optional): Plot title.
        """
        if self.trends is None or self.p_values is None:
            raise ValueError("Trend data or p-values not loaded.")

        # Select the trend data and p-values for the variable
        trend_data = self.trends[var]
        p_values = self.p_values[var]

        # Create a stippling mask where p-values indicate significance
        stipple_mask = p_values < 0.05

        # Call the abstract plot function with the stippling mask
        self.plot(trend_data, clabel, title, stipple_mask=stipple_mask)

    def plot_1d(self, var, xlabel, ylabel, title=None, n=None):
        """
        Plot some 1D data.

        Parameters
        ----------
        var (str): Variable name
        n (int): How many data points should be plotted?
        xlabel (str): x axis label
        ylabel (str): y axis label
        title (str, optional): Plot title
        """
        if self.dataset is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        # Extract eigenvalues from the 'zg' variable
        data = self.dataset[var].values.flatten()

        if n is not None:
            # Ensure we do not exceed available eigenvalues
            n = min(n, len(data))
        else:
            n = len(data)

        # Generate x values (e.g., 1 to n)
        x_values = range(1, n + 1)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(x_values, data[:n], 'o', color='b', markersize=8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title is not None:
            plt.title(title)
        # Set x-ticks to be integers
        ax.set_xticks(x_values)
        # Adjust grid lines to match tick positions
        ax.grid(True, which='both', linestyle='--', linewidth=0.7)

        # Ensure tick labels are properly formatted
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        plt.show()
        return fig, ax

    def plot_stacked_lines(self, var, ylabel="Plot {}", n=None):
        """
        Plot the time series of the given data as stacked subplots.

        Parameters
        ----------
        var (str): Variable name
        ylabel (str): y axis label
        n (int): Plot first n time series only.
        """
        if self.dataset is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        data = self.dataset[var].values
        time = self.dataset['time'].values
        time = pd.to_datetime(time)
        if n is not None:
            n = min(n, data.shape[1])
        else:
            n = data.shape[1]
        fig, axs = plt.subplots(n, 1,
                                figsize=(10, 2 * n),
                                sharex=True, sharey=True)
        for i in range(n):
            ax = axs[i]
            # Plot all time series in grey, with the i-th in black
            for j in range(20):
                color = 'k' if j == i else 'grey'
                alpha = 1 if j == i else .4
                zorder = 10 if j == i else 1
                ax.plot(time, data[:, j].flatten(), color=color, alpha=alpha,
                        zorder=zorder, linewidth=.1)
            ax.set_ylabel(ylabel.format(i))
            ax.grid(True, linestyle='--', linewidth=0.5)
            # ax.get_yaxis().set_visible(False)
            # Remove spines
            for spine in ax.spines.values():
                spine.set_visible(False)
            # Remove horizontal grid lines
            ax.yaxis.grid(False)
        # Set x-axis date formatting
        axs[-1].xaxis.set_major_locator(mdates.YearLocator())
        axs[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

        # Rotate and format x-tick labels
        plt.setp(axs[-1].xaxis.get_majorticklabels(), rotation=45, ha='right')

        # Set common x-label
        axs[-1].set_xlabel('Date')

        plt.tight_layout()
        plt.show()
        return fig, axs
