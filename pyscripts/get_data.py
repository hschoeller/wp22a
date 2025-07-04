# Script to download control ensemble member of ERA5
# eda.
import data_utils as du
import sys
import os
import argparse

sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))


d_path_ens = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"
d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/data/"

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
        "ens",
        nargs='?',  # Make optional for backwards compatibility
        type=str,
        help="Ensemble Spread of Analysis?"
    )

    return parser.parse_args()


def main():

    args = parse_arguments()

    var_name = args.variable_name
    ens = args.ens.lower() == 'true' if args.ens else False

    years = [f'{y}' for y in range(year_min, year_max + 1)]
    if ens:
        print("ens")
        grib_file = du.retrieve_ens(years, lat_min, lat_max, lon_min,
                                    lon_max, d_path_ens,
                                    f'{year_min}_{year_max}_ens.zip', var_name)
        print(f"grib_file: {grib_file}")
        du.convert_grib_to_nc(d_path_ens, grib_file, cleanup=True)

    else:
        print("no ens")
        grib_file = du.retrieve_single(
            years, lat_min, lat_max, lon_min, lon_max, d_path,
            f'{year_min}_{year_max}_{var_name}.zip', var_name)
        print(f"grib_file: {grib_file}")
        du.convert_grib_to_nc(d_path, grib_file, cleanup=True)


if __name__ == "__main__":
    main()
