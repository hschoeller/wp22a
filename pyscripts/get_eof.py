import xarray as xr
import dask.array as da
import numpy as np
from scipy import stats
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from dask_ml.decomposition import IncrementalPCA


class EOFAnalysis:
    """Class to handle EOF analysis on ERA5 data."""

    def __init__(self, data_path, latitude_bounds, longitude_bounds, n_components, n_chunks):
        """
        Initialize the EOFAnalysis instance with dataset path, spatial bounds, 
        and EOF configuration.

        Parameters:
            data_path (str): Path to ERA5 dataset files.
            latitude_bounds (tuple): Latitude bounds for data selection.
            longitude_bounds (tuple): Longitude bounds for data selection.
            n_components (int): Number of EOF components to retain.
        """
        self.data_path = data_path
        self.latitude_bounds = latitude_bounds
        self.longitude_bounds = longitude_bounds
        self.n_components = n_components
        self.n_chunks = n_chunks
        self.ds = None
        self.z500_anomaly = None
        self.z500_anomaly_detrend = None
        self.z500_anomaly_norm = None

    def setup_dask_cluster(self, cores, memory, walltime, queue, nodelist):
        """
        Initialize and connect to a Dask cluster with SLURM.

        Parameters:
            cores (int): Number of CPU cores.
            memory (str): Memory limit for the job.
            walltime (str): Walltime for SLURM job.
            queue (str): SLURM queue/partition.
            nodelist (str): Specific node to use.

        Returns:
            Client: Dask distributed client instance.
        """
        cluster = SLURMCluster(
            job_name='anos',
            cores=cores,
            memory=memory,
            walltime=walltime,
            queue=queue,
            nodelist=nodelist
        )
        cluster.scale(jobs=1)
        client = Client(cluster)
        return client

    def load_and_preprocess_data(self):
        """Load ERA5 data, select spatial domain, and calculate anomalies."""
        self.ds = xr.open_mfdataset(self.data_path, parallel=True,
                                    combine="by_coords")
        z500 = self.ds['z500'].sel(lat=slice(*self.latitude_bounds),
                                   lon=slice(*self.longitude_bounds))
        climatology = self._calculate_climatology(z500)
        self.z500_anomaly = self._calculate_anomaly(z500, climatology)
        self.z500_anomaly_detrend = self._detrend_data(self.z500_anomaly)
        self.z500_anomaly_norm = self._normalize_data(
            self.z500_anomaly_detrend)

    def _calculate_climatology(self, data):
        """
        Calculate the climatology for each grid point and each day of the year
        (across all years), i.e., the mean of each calendar day's value over all years.

        Parameters:
            data (xarray.DataArray): The data array containing time, lat, and lon dimensions.

        Returns:
            xarray.DataArray: The climatology for each day of the year at each grid point.
        """
        # Group by the day of the year and calculate the mean over all years
        # """Calculate the 90-day rolling mean climatology."""
        # return (
        #     data.groupby("time.dayofyear")
        #         .rolling(time=90, center=True, min_periods=1)
        #         .mean()
        #         .groupby("time.dayofyear")
        #         .mean(dim='time')
        # )
        climatology = data.groupby("time.dayofyear").mean(dim="time")

        return climatology

    def _calculate_anomaly(self, data, climatology):
        """Calculate anomalies by subtracting the climatology."""
        return data.groupby("time.dayofyear") - climatology

    def _detrend_data(self, data):
        """
        Detrend the data at each grid point by removing a linear trend.

        Parameters:
            data (xarray.DataArray): The data array containing time, lat, and lon dimensions.

        Returns:
            xarray.DataArray: The detrended data with long-term trends removed.
        """
        # Convert time to an integer index (for linear regression)
        time_idx = np.arange(len(data.time))

        # Define a function to apply linear detrending to each grid point
        def detrend_grid_point(ts):
            # Perform linear regression to find the slope and intercept
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                time_idx, ts)
            # Subtract the trend (slope * time + intercept) to obtain anomalies
            detrended_ts = ts - (slope * time_idx + intercept)
            return detrended_ts

        # Apply the detrending function to each grid point
        detrended_data = data.map(detrend_grid_point)

        return detrended_data

    def _normalize_data(self, data):
        """
        Apply a 10-day low-pass filter, calculate 30-day running standard deviation, 
        and normalize the data by removing seasonal amplitude.
        """
        # 10-day low-pass filter
        # z500_filtered = self.z500_anomaly.rolling(time=10, center=True).mean()

        # 30-day running standard deviation
        std_30day = (
            data.groupby("time.dayofyear")
            .rolling(time=30, center=True, min_periods=1)
            .std()
            .groupby("time.dayofyear")
            .mean(dim='time')
        )

        # Calculate the cosine of the latitude for each latitude point
        cos_lat = np.cos(np.deg2rad(std_30day.lat))

        # Weight the data by the cosine of the latitude and then calculate the mean over lat and lon
        weighted_std = std_30day * cos_lat

        # Now calculate the spatial mean
        mean_std = weighted_std.mean(dim=['lat', 'lon']) / cos_lat.mean()
        return data / mean_std

    def reshape_for_eof(self):
        """Flatten data for EOF analysis, returning time and space dimensions."""
        return self.z500_anomaly_norm.stack(space=('lat', 'lon')).transpose('time', 'space')

    def perform_eof_analysis(self, data):
        """
        Perform EOF analysis using IncrementalPCA.

        Parameters:
            data (dask.array): Data array reshaped for EOF analysis.

        Returns:
            tuple: PCs (time x n_components) and EOFs (n_components x latitude x longitude).
        """
        pca = IncrementalPCA(n_components=self.n_components)
        pcs = pca.fit_transform(data)
        eofs = pca.components_.reshape(
            self.n_components,
            self.z500_anomaly_norm.sizes['lat'],
            self.z500_anomaly_norm.sizes['lon']
        )
        return pcs, eofs

    def save_results(self, pcs, eofs):
        """Save the EOF and PC results as an xarray Dataset to a NetCDF file."""
        eof_ds = xr.Dataset(
            {
                'EOFs': (['mode', 'latitude', 'longitude'], eofs),
                'PCs': (['time', 'mode'], pcs)
            },
            coords={
                'mode': np.arange(self.n_components),
                'lat': self.z500_anomaly_norm['lat'],
                'lon': self.z500_anomaly_norm['lon'],
                'time': self.z500_anomaly_norm['time']
            }
        )
        eof_ds.to_netcdf("EOF_results.nc")


def main():
    # Initialize EOF Analysis
    ERA5_dir = ("/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/" +
                "ERA5/day/atmos/zg/r1i1p1-050000Pa/")

    files = [ERA5_dir + f"zg_day_reanalysis_era5_r1i1p1-050000Pa_{y}0101-{y}1231.nc"
             for y in range(1950, 1953)]
    eof_analysis = EOFAnalysis(
        data_path=files,
        latitude_bounds=(20., 90.),  # 249 =3*83 grid points
        longitude_bounds=(250., 360.),
        n_components=50,
        chunks={
            'lon': 49,
            'lat': 83}
    )

    # Setup Dask Cluster
    client = eof_analysis.setup_dask_cluster(
        cores=48,
        memory="180GB",
        walltime="02:00:00",
        queue='main',
        nodelist='calc04'
    )

    print("Dask cluster and client initialized:", client)

    # Load, preprocess, filter, and normalize data
    eof_analysis.load_and_preprocess_data()
    eof_analysis.filter_and_normalize_data()

    # Reshape for EOF Analysis
    reshaped_data = eof_analysis.reshape_for_eof()

    # Perform EOF Analysis
    pcs, eofs = eof_analysis.perform_eof_analysis(reshaped_data)

    # Save Results
    eof_analysis.save_results(pcs, eofs)

    # Close Dask client and cluster
    client.close()


if __name__ == "__main__":
    main()
