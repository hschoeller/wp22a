import xarray as xr
import numpy as np
import os
import xesmf as xe


def split_nc_file(input_file_path, output_dir, num_chunks):
    """
    Splits a large NetCDF file into smaller chunks along latitude and longitude.

    Args:
        input_file_path (str): Path to the input NetCDF file.
        output_dir (str): Directory where the smaller chunks will be saved.
        num_chunks (int): Number of chunks to split the file into.

    Returns:
        None
    """
    # Load the NetCDF file
    dataset = xr.open_dataset(input_file_path)
    print(dataset)
    lat_start = np.ceil(dataset.latitude.min())
    lon_start = np.ceil(dataset.longitude.min())

    target_lats = np.arange(float(lat_start),
                            float(dataset.latitude.max()) + 1, 1)
    target_lons = np.arange(float(lon_start),
                            float(dataset.longitude.max()) + 1, 1)

    target_grid = xr.Dataset(
        {"latitude": (["latitude"], target_lats),
         "longitude": (["longitude"], target_lons)}
    )

    regridder = xe.Regridder(dataset, target_grid, method="conservative")
    dataset = regridder(dataset)
    # Get latitude and longitude coordinates
    lats = dataset['latitude'].values
    lons = dataset['longitude'].values

    # Calculate number of chunks along each axis (latitude and longitude)
    lat_chunks = num_chunks[0]
    lon_chunks = num_chunks[1]

    # Split the latitude and longitude ranges
    lat_indices = np.linspace(0, len(lats), lat_chunks + 1, dtype=int)
    lon_indices = np.linspace(0, len(lons), lon_chunks + 1, dtype=int)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize chunk ID counter
    chunk_id = 1

    # Loop through latitude and longitude ranges to create subsets
    for i in range(lat_chunks):
        for j in range(lon_chunks):
            print("chunk_id", chunk_id)
            # Define lat/lon ranges for the current chunk
            lat_start, lat_end = lat_indices[i], lat_indices[i + 1]
            lon_start, lon_end = lon_indices[j], lon_indices[j + 1]
            print(f"Latitude: {lat_start} to {lat_end}")
            print(f"Longitude: {lon_start} to {lon_end}")
            # Subset the data
            subset = dataset.isel(
                latitude=slice(lat_start, lat_end),
                longitude=slice(lon_start, lon_end)
            )

            for var in ['latitude', 'longitude']:
                if var in dataset.variables:
                    subset[var].encoding = dataset[var].encoding

            # Generate a meaningful filename for the chunk
            output_file_path = os.path.join(
                output_dir,
                f"chunk_{chunk_id:02d}.nc"  # Sequential numbering
            )
            if subset.sizes['latitude'] == 0 or subset.sizes['longitude'] == 0:
                print(f"Skipping empty subset for chunk {chunk_id}")
                continue
            print(f"Saving chunk to {output_file_path}")
            # Save the subset to a new NetCDF file
            subset.to_netcdf(output_file_path)

            chunk_id += 1

    # Close the dataset
    dataset.close()

    print(f"Successfully split the file into {chunk_id - 1} chunks.")


if __name__ == "__main__":
    # Define input parameters
    # Path to the input NetCDF file
    input_file_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/data.nc"
    # Directory for output chunks
    output_dir = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/chunks"
    num_chunks = [3, 13]  # Number of desired chunks

    # Call the function to split the file
    split_nc_file(input_file_path, output_dir, num_chunks)
