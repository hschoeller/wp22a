# -*- coding: utf-8 -*-
import sys
import os

d_path = "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"
chunkdict = 'auto'

# Input and output file paths
input_file = "data.grib"
not_years = [2018, 2017, 2016, 2007, 1997, 1987, 1977, 1967, 1957, 1947]
years = [year for year in range(1980, 2024) if year not in not_years]
year = years[int(sys.argv[1])]

output_file = f"spread_{year}.nc"

if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join('..', 'pyscripts')))
    import data_utils as du
    du.calculate_ensemble_spread(d_path + input_file, d_path + output_file,
                                 year, chunkdict)
