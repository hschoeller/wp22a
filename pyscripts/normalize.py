import xarray as xr
import sys

def normalize_data(anosdetrend_file, fldmean_file, output_file):
    # Load the datasets
    anosdetrend_var = xr.open_dataset(anosdetrend_file)
    fldmean_var = xr.open_dataset(fldmean_file)

    fldmean_var = fldmean_var.squeeze()

    # Compute the mean field (if not already a time series, we assume it is)
    # Broadcasting the mean over the time dimension
    normalized_ds = anosdetrend_var / fldmean_var

    # Save the output to a new NetCDF file
    normalized_ds.to_netcdf(output_file)

    # Close datasets
    anosdetrend_var.close()
    fldmean_var.close()
    normalized_ds.close()

if __name__ == '__main__':
    # Get command line arguments
    anosdetrend_file = sys.argv[1]
    fldmean_file = sys.argv[2]
    output_file = sys.argv[3]

    # Call the normalization function
    normalize_data(anosdetrend_file, fldmean_file, output_file)
