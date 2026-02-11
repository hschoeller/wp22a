import xarray as xr
import numpy as np
import os
import argparse
import sys

def split_nc_file(input_file_path, output_dir, num_chunks):
    """
    Split a NetCDF file into contiguous chunks along latitude and longitude.
    - If latitude is not a dimension, a singleton latitude dimension is added.
    - Ensures chunk counts do not exceed available points.
    - Preserves each variable's original dimension ordering and encodings.

    Parameters
    ----------
    input_file_path : str
        Path to input NetCDF.
    output_dir : str
        Directory to write chunk files.
    num_chunks : sequence of two ints
        Desired number of chunks along (latitude, longitude).

    Returns
    -------
    list
        List of saved chunk file paths.
    """
    os.makedirs(output_dir, exist_ok=True)

    ds = xr.open_dataset(input_file_path)

    # capture original variable dims and encodings BEFORE we rename/expand anything
    orig_var_dims = {name: tuple(ds[name].dims) for name in ds.variables}
    orig_encodings = {name: dict(ds[name].encoding) if hasattr(ds[name], "encoding") else {} for name in ds.variables}

    # find coordinate names that look like lat/lon
    coord_names = list(ds.coords)
    lat_name = next((n for n in coord_names if 'lat' in n.lower()), None)
    lon_name = next((n for n in coord_names if 'lon' in n.lower()), None)

    if lat_name is None or lon_name is None:
        ds.close()
        raise ValueError("Latitude or longitude coordinate not found in dataset.")

    # We'll map original coordinate names to canonical names used internally
    name_map = {lat_name: 'latitude', lon_name: 'longitude'}

    # Ensure latitude exists as a dimension; if not, add a singleton latitude dimension
    lat_vals = ds[lat_name].values
    if lat_name not in ds.dims:
        lat_arr = np.atleast_1d(lat_vals)
        if lat_arr.size == 1:
            ds = ds.expand_dims({'latitude': lat_arr})
            if lat_name != 'latitude':
                ds = ds.rename({lat_name: 'latitude'})
            lat_coord = 'latitude'
        else:
            # length > 1 but not a dim: add singleton latitude equal to mean
            ds = ds.expand_dims({'latitude': np.array([float(np.mean(lat_arr))])})
            lat_coord = 'latitude'
            if lat_name != 'latitude':
                # if original coord still exists as variable, rename it to avoid collision
                if lat_name in ds.variables:
                    ds = ds.rename({lat_name: 'latitude_coord_orig'})
    else:
        # lat_name is already a dimension
        if lat_name != 'latitude':
            ds = ds.rename({lat_name: 'latitude'})
        lat_coord = 'latitude'

    # Ensure longitude is a dimension; if not, convert it to one
    lon_vals = ds[lon_name].values
    if lon_name not in ds.dims:
        lon_arr = np.atleast_1d(lon_vals)
        if lon_arr.size == 1:
            ds = ds.expand_dims({'longitude': lon_arr})
            if lon_name != 'longitude':
                ds = ds.rename({lon_name: 'longitude'})
            lon_coord = 'longitude'
        else:
            # create an integer dimension and assign the coordinate values
            ds = ds.expand_dims({'longitude': np.arange(lon_arr.size)})
            ds = ds.assign_coords(longitude = ('longitude', lon_arr))
            lon_coord = 'longitude'
    else:
        if lon_name != 'longitude':
            ds = ds.rename({lon_name: 'longitude'})
        lon_coord = 'longitude'

    # canonical coordinate arrays
    lat_arr = np.atleast_1d(ds[lat_coord].values)
    lon_arr = np.atleast_1d(ds[lon_coord].values)
    n_lat = len(lat_arr)
    n_lon = len(lon_arr)

    # clamp chunk counts
    lat_chunks_req = int(num_chunks[0])
    lon_chunks_req = int(num_chunks[1])
    lat_chunks = max(1, min(lat_chunks_req, n_lat))
    lon_chunks = max(1, min(lon_chunks_req, n_lon))

    # compute integer edges for contiguous slices
    lat_edges = np.linspace(0, n_lat, lat_chunks + 1, dtype=int)
    lon_edges = np.linspace(0, n_lon, lon_chunks + 1, dtype=int)

    # create a mapping from original variable names to their new names (only possible renames are lat/lon)
    # This is used to remap encodings captured earlier
    encoding_map = {}
    for orig_name, enc in orig_encodings.items():
        new_name = name_map.get(orig_name, orig_name)
        encoding_map[new_name] = enc

    out_files = []
    chunk_id = 1

    for i in range(lat_chunks):
        lat_start = int(lat_edges[i])
        lat_end   = int(lat_edges[i + 1])   # exclusive
        if lat_end <= lat_start:
            continue
        for j in range(lon_chunks):
            lon_start = int(lon_edges[j])
            lon_end   = int(lon_edges[j + 1])  # exclusive
            if lon_end <= lon_start:
                continue
            
            subset = ds.isel(latitude = slice(lat_start, lat_end),
                             longitude = slice(lon_start, lon_end))

            # skip empty subsets defensively
            if subset.sizes.get('latitude', 0) == 0 or subset.sizes.get('longitude', 0) == 0:
                continue

            # restore encodings for coordinates/variables when available
            for var in list(subset.variables):
                if var in encoding_map:
                    # copy encoding dict to avoid accidental shared references
                    try:
                        subset[var].encoding = dict(encoding_map[var])
                    except Exception:
                        # some encodings may be non-assignable; ignore in that case
                        pass

            # Ensure each variable's dim order matches the original variable dim order (mapped to new names)
            for var in list(subset.variables):
                if var in orig_var_dims:
                    orig_order = list(orig_var_dims[var])  # original dims (may include original lat/lon names)
                    # map lat/lon original names to their new canonical names if present
                    mapped_order = [name_map.get(d, d) for d in orig_order]
                    # keep only dims that actually exist in the subset variable (preserve relative order)
                    target = [d for d in mapped_order if d in subset[var].dims]
                    # append any remaining dims from the subset's current dims (they were not in original order)
                    remaining = [d for d in subset[var].dims if d not in target]
                    final_order = target + remaining

                    # if final_order differs from current order, transpose
                    if tuple(final_order) != tuple(subset[var].dims):
                        try:
                            subset[var] = subset[var].transpose(*final_order)
                        except Exception:
                            # If transpose fails for some reason, leave as-is (defensive)
                            pass

            output_file_path = os.path.join(output_dir, f"chunk_{chunk_id:02d}.nc")
            subset.to_netcdf(output_file_path)
            out_files.append(output_file_path)
            chunk_id += 1

    ds.close()
    return out_files





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

    return parser.parse_args()


def validate_inputs(args):
    """Validate command line arguments."""
    if args.input_file and not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    if args.input_file and not args.input_file.lower().endswith(('.nc', '.netcdf')):
        raise ValueError("Input file must be a NetCDF file (.nc or .netcdf)")

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

    print(f"Input file: {input_file_path}")
    print(f"Output directory: {output_dir}")
    print(f"Chunks: {num_chunks}")
    print("=" * 40)

    # Call the function to split the file
    try:
        split_nc_file(input_file_path, output_dir, num_chunks)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
