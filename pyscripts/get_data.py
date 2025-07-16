# Script to download control ensemble member of ERA5
# eda.
import data_utils as du
import sys
import os
import argparse

sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))


# d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/"

lon_min, lon_max = -80, 40  # Longitude range
lat_min, lat_max = 30, 90    # Latitude range
year_min = 1940
year_max = 2024


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download specified ERA5 variable",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        """
    )

    parser.add_argument(
        "variable_name",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Variable name in ERA5"
    )

    parser.add_argument(
        "dataset",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Dataset name in ERA5"
    )

    parser.add_argument(
        "product_type",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Product type in ERA5"
    )
    parser.add_argument(
        "ens",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Ensemble Spread of Analysis?"
    )

    return parser.parse_args()


def main():

    args = parse_arguments()

    var_name = args.variable_name
    d_path_full = f"{d_path}/ens_data/" if args.ens.lower() == 'true' else f"{d_path}/data/"

    years = [f'{y}' for y in range(year_min, year_max + 1)]
    grib_file = du.retrieve_era5(
        years=years,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        d_path=d_path_full,
        f_name=f"{var_name}.zip",
        var=var_name,
        dataset=args.dataset,
        product_type=args.product_type,
        pressure_level=None if args.dataset == "reanalysis-era5-single-levels" else "500"
    )
    print(f"grib_file: {grib_file}")
    du.convert_grib_to_nc(d_path_full, grib_file, cleanup=True)


if __name__ == "__main__":
    main()
