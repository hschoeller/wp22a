import xarray as xr
import numpy as np
import os
import xesmf as xe
import argparse
import sys


def split_nc_file(input_file_path, output_dir, num_chunks, do_regrid=True):
    """
    Splits a large NetCDF file into smaller chunks along latitude and longitude.
    Args:
        input_file_path (str): Path to the input NetCDF file.
        output_dir (str): Directory where the smaller chunks will be saved.
        num_chunks (list): Number of chunks to split the file into [lat_chunks, lon_chunks].
        do_regrid (bool): Whether to perform regridding (default: True for backwards compatibility).
    Returns:
        None
    """
    os.makedirs(output_dir, exist_ok=True)
    # Load the NetCDF file
    dataset = xr.open_dataset(input_file_path)
    print(dataset)

    if do_regrid:
        print("Performing regridding...")
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
        print("Regridding completed.")
    else:
        print("Skipping regridding.")

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


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Split NetCDF files into smaller chunks along latitude and longitude",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using command line arguments
  python script.py input.nc output_folder --chunks 1 48
  python script.py input.nc output_folder --chunks 4 4 --no-regrid
  python script.py input.nc output_folder --chunks 2 3 --regrid
  
  # The script will also work with the original hardcoded values if no arguments provided
        """
    )

    parser.add_argument(
        "input_file",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Path to input NetCDF file"
    )

    parser.add_argument(
        "output_folder",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Output folder for chunked files"
    )

    parser.add_argument(
        "--chunks",
        nargs=2,
        type=int,
        metavar=("LAT_CHUNKS", "LON_CHUNKS"),
        help="Number of chunks along latitude and longitude dimensions"
    )

    parser.add_argument(
        "--regrid",
        action="store_true",
        help="Enable regridding (default behavior when using original hardcoded values)"
    )

    parser.add_argument(
        "--no-regrid",
        action="store_true",
        help="Disable regridding"
    )

    return parser.parse_args()


def validate_inputs(args):
    """Validate command line arguments."""
    if args.input_file and not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    if args.input_file and not args.input_file.lower().endswith(('.nc', '.netcdf')):
        raise ValueError("Input file must be a NetCDF file (.nc or .netcdf)")

    if args.regrid and args.no_regrid:
        raise ValueError("Cannot specify both --regrid and --no-regrid")

    if args.chunks and (args.chunks[0] <= 0 or args.chunks[1] <= 0):
        raise ValueError("Number of chunks must be positive integers")

    # If using command line args, all required arguments must be present
    if args.input_file or args.output_folder or args.chunks:
        if not all([args.input_file, args.output_folder, args.chunks]):
            raise ValueError(
                "When using command line arguments, input_file, output_folder, and --chunks are all required")


def main():
    """Main function that handles both command line and backwards compatibility."""
    args = parse_arguments()

    try:
        validate_inputs(args)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    input_file_path = args.input_file
    output_dir = args.output_folder
    num_chunks = args.chunks

    # Determine regridding setting
    if args.no_regrid:
        do_regrid = False
    elif args.regrid:
        do_regrid = True
    else:
        # Default to True for backwards compatibility
        do_regrid = True

    print(f"Input file: {input_file_path}")
    print(f"Output directory: {output_dir}")
    print(f"Chunks: {num_chunks}")
    print(f"Regridding: {'Enabled' if do_regrid else 'Disabled'}")
    print("=" * 40)

    # Call the function to split the file
    try:
        split_nc_file(input_file_path, output_dir, num_chunks, do_regrid)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
