library(ncdf4)

# Function to split the data.nc into 39 smaller files (including time dimension)
split_data_into_chunks <- function(input_file, num_chunks = 39) {
    print("Open file")
    # Open the original NetCDF file
    nc_file <- nc_open(input_file)

    # Get dimensions of the data
    lon_data <- ncvar_get(nc_file, "longitude")
    lon_atts <- ncatt_get(nc_file, "longitude")
    lat_data <- ncvar_get(nc_file, "latitude")
    lat_atts <- ncatt_get(nc_file, "latitude")
    time_data <- ncvar_get(nc_file, "valid_time")
    time_atts <- ncatt_get(nc_file, "valid_time")
    number_data <- ncvar_get(nc_file, "number")
    number_atts <- ncatt_get(nc_file, "number")
    print("get z data")
    z_data <- ncvar_get(nc_file, "z")
    global_atts <- ncatt_get(nc_file, 0)
    z_atts <- ncatt_get(nc_file, "z")

    # Close the original file
    nc_close(nc_file)
    # Create grid points (expand.grid of longitude and latitude)
    grid_points <- expand.grid(
        lon = 1:length(lon_data),
        lat = 1:length(lat_data)
    )

    # Calculate the chunk size based on grid points
    chunk_size <- ceiling(nrow(grid_points) / num_chunks)
    print("get chunks")
    # Split the grid points into chunks
    grid_chunks <- split(
        grid_points,
        ceiling(seq_along(1:nrow(grid_points)) / chunk_size)
    )

    # Loop through each chunk and create a new NetCDF file
    for (i in seq_along(grid_chunks)) {
        print(paste0("chunk", i))
        chunk <- grid_chunks[[i]]

        # Create new dimension definitions for this chunk
        lon_dim <- ncdim_def("longitude", "degrees_east", lon_data[chunk$lon])
        lat_dim <- ncdim_def("latitude", "degrees_north", lat_data[chunk$lat])
        time_dim <- ncdim_def("time", "seconds since 1970-01-01", time_data)
        number_dim <- ncdim_def("number", "number", number_data)

        # Create variables for the chunk file
        z_var <- ncvar_def("z", "units_of_z",
            list(lon_dim, lat_dim, time_dim, number_dim),
            missval = NA
        )
        print("create file")
        # Create a new NetCDF file for the chunk (skip the existing longitude,
        # latitude, and time)
        output_file <- paste0("../ens_data/data_", i, ".nc")
        nc_out <- nc_create(output_file, z_data[chunk$lon, chunk$lat, , ])
        for (att_name in names(global_atts)) {
            ncatt_put(nc_out, 0, att_name, global_atts[[att_name]])
        }
        for (att_name in names(z_atts)) {
            ncatt_put(nc_out, "z", att_name, z_atts[[att_name]])
        }
        for (att_name in names(lon_atts)) {
            ncatt_put(nc_out, "longitude", att_name, lon_atts[[att_name]])
        }
        for (att_name in names(lat_atts)) {
            ncatt_put(nc_out, "latitude", att_name, lat_atts[[att_name]])
        }
        for (att_name in names(time_atts)) {
            ncatt_put(nc_out, "valid_time", att_name, time_atts[[att_name]])
        }
        for (att_name in names(number_atts)) {
            ncatt_put(nc_out, "number", att_name, number_atts[[att_name]])
        }
        print("put z data")
        # Write the variables to the new chunk file
        ncvar_put(nc_out, "z", z_data)

        # Write the longitude, latitude, and time variables to the chunk file (no need to redefine them)
        ncvar_put(nc_out, "longitude", lon_data[chunk$lon])
        ncvar_put(nc_out, "latitude", lat_data[chunk$lat])
        ncvar_put(nc_out, "valid_time", time_data)
        ncvar_put(nc_out, "number", number_data)

        # Close the chunk file
        nc_close(nc_out)
        print("close file")
    }
}

# Run the function to split the original data.nc into 39 smaller files
split_data_into_chunks("../ens_data/data.nc")
