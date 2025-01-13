import xarray as xr
import xeofs as xe
import sys

if __name__ == "__main__":
    d_path = "../data_ens/"
    chunkdict = 'auto'
    if len(sys.argv) > 1:
        # Open the dataset (already preprocessed)
        zg_norm = xr.open_dataset(d_path + f'zg_norm{sys.argv[1]}.nc')['zg']
    else:
        zg_norm = xr.open_dataset(d_path + 'zg_norm.nc')['zg']
    # Ensure data is chunked along time (Dask-friendly)
    zg_norm = zg_norm.chunk(chunkdict)

    if len(sys.argv) > 1:
        solver = xe.single.EOF.load(d_path + 'solver')
        pcs = solver.transform(zg_norm)
        pcs.to_netcdf(d_path + 'projected_pcs.nc')
    else:
        # Instantiate the EOF solver with 50 EOFs
        solver = xe.single.EOF(n_modes=50, use_coslat=True)

        # Fit the solver (computes EOFs)
        solver.fit(zg_norm, "time")

        solver.save(d_path + 'solver', overwrite=True)

        # Access the EOFs and PCs
        eofs = solver.components()  # EOF spatial patterns
        pcs = solver.scores()    # Principal components
        ev = solver.explained_variance_ratio()

        eofs.attrs.pop("solver_kwargs", None)
        pcs.attrs.pop("solver_kwargs", None)
        ev.attrs.pop("solver_kwargs", None)

        eofs.to_netcdf(d_path + 'eofs.nc')
        pcs.to_netcdf(d_path + 'pcs.nc')
        ev.to_netcdf(d_path + 'evals.nc')
    print("EOF computation completed.")
