# Script to download control ensemble member of ERA5
# eda.
import sys
import os
import data_utils as du
sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))


d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

lon_min, lon_max = -100, 20  # Longitude range
lat_min, lat_max = 20, 90    # Latitude range
year_min = 1940
year_max = 2024

if __name__ == "__main__":

    years = [f'{y}' for y in range(year_min, year_max + 1)]
    du.retrieve_ens(years, lat_min, lat_max, lon_min,
                    lon_max, d_path,
                    f'{year_min}_{year_max}_ens.zip')
