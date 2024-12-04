import xarray as xr

xr.open_mfdataset(["/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/" +
                   "ERA5/day/atmos/zg/r1i1p1-050000Pa/zg_day_reanalysis_era5_r1i1p1-050000Pa_19500101-19501231.nc", "/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/" +
                   "ERA5/day/atmos/zg/r1i1p1-050000Pa/zg_day_reanalysis_era5_r1i1p1-050000Pa_19510101-19511231.nc"], combine='by_coords')
