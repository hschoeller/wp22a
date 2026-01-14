source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")


args <- commandArgs(trailingOnly = TRUE)
CHUNK_DIR <- paste0(args[1], "zg_chunks/")
CHUNK_NO <- as.integer(args[2]) +1
OUT_DIR <- paste0(args[1], "mods/")

# Change points
# Convert to Date objects for easier comparison
change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

# grid-point-wise function
process_grid_point <- function(lat, lon, z, time, cps, smooth_mod) {
    print("pull data together")

    z_df <- data.frame(time = time, z = z)
    print("Data frame created")
    z_df$year <- as.integer(format(z_df$time, "%Y"))
    z_df$doy <- as.integer(format(z_df$time, "%j"))
    print("now combine")

    z_df$sin_doy <- sin(2 * pi * z_df$doy / 365)
    z_df$cos_doy <- cos(2 * pi * z_df$doy / 365)
    z_df$log_variance <- log(z_df$z)


    z_df$segment <- cut(as.Date(z_df$time),
        breaks = c(as.Date(c(
            z_df$time[1],
            cps,
            z_df$time[nrow(z_df)]
        ))),
        labels = 0:(length(cps)) + 1,
        include.lowest = TRUE,
        right = FALSE
    )
    z_df <- z_df %>%
        group_by(segment) %>%
        mutate(day_no = row_number()) %>%
        ungroup()

    if (smooth_mod) {
        print("Estimate smooth model")
        # Estimate smooth model
        lmod_seas_smooth <- gamm(
            log_variance ~
                s(doy, by = segment, bs = "cc") +
                s(year, by = segment, bs = "cr") +
                ti(year, doy, by = segment, bs = c("cr", "cc")),
            data = z_df,
            method = "REML",
            correlation = corAR1(form = ~ day_no | segment),
            weights = varIdent(form = ~ 1 | segment),
            keepData = FALSE
        )
        print("Saving smooth model")
        saveRDS(lmod_seas_smooth, file = paste0(
            OUT_DIR, "slm", lat, "_", lon,
            ".Rds"
        ))
        print("Smooth model saved")
    } else {
        print("Estimate linear model")
        lmod_seas <- nlme::gls(
            log_variance ~ segment +
                segment:year + segment:sin_doy + segment:cos_doy - 1,
            data = z_df,
            correlation = corAR1(form = ~ day_no | segment),
            weights = varIdent(form = ~ 1 | segment)
        )
        print("Saving linear model")
        saveRDS(lmod_seas, file = paste0(
            OUT_DIR, "lm", lat, "_", lon,
            ".Rds"
        ))

        print("Linear model saved")
    }
    print("Done")
}

# Define the wrapper for parallel execution
process_wrapper <- function(chunk_idx) # , grid_points, lat_data, lon_data, z_data,
# time_data, change_points)
{
    # redirect output into separate files
    output_file <- paste0(OUT_DIR, "logs/output_process_", chunk_idx, ".log")
    error_file <- paste0(OUT_DIR, "logs/error_process_", chunk_idx, ".log")
    output_con <- file(output_file, open = "wt")
    error_con <- file(error_file, open = "wt")
    # Redirect standard output and standard error
    sink(output_con) # Redirect standard output
    sink(error_con, type = "message") # Redirect standard error

    on.exit({
        sink() # Reset standard output
        sink(type = "message") # Reset standard error
        close(output_con) # Close output connection
        close(error_con) # Close error connection
    })
    # Read in the data
    chunk_file <- paste0(CHUNK_DIR, sprintf("chunk_%02d.nc", chunk_idx))

    print("Reading data")
    nc_file <- nc_open(chunk_file)
    print("Data opened; reading z")
    z_data <- ncvar_get(nc_file, "z", collapse_degen = FALSE)
    print("reading latlon")
    # Get coordinates
    lon_data <- ncvar_get(nc_file, "longitude")
    lat_data <- ncvar_get(nc_file, "latitude")
    print(lat_data)
    print(lon_data)
    grid_points <- expand.grid(1:length(lon_data), 1:length(lat_data))
    print("reading time")
    time_data <- ncvar_get(nc_file, "valid_time")

    # Convert time
    time_origin <- sub(
        "seconds since ", "",
        ncatt_get(nc_file, "valid_time", "units")$value
    )
    time <- as.Date(as.POSIXct(time_data, origin = time_origin, tz = "UTC"))

    nc_close(nc_file)
    print("Data read")
    for (smooth_mod in c(FALSE)) { # }, TRUE)) {
        for (i in seq_len(nrow(grid_points))) {
            lon_idx <- grid_points[i, 1]
            lat_idx <- grid_points[i, 2]
            print(paste0("Processing grid point ", lat_idx, ", ", lon_idx))
            print(paste0(
                "Lat: ", lat_data[lat_idx],
                ", Lon: ", lon_data[lon_idx]
            ))
            file_path <- paste0(
                OUT_DIR, if (smooth_mod) "slm" else "lm",
                lat_data[lat_idx], "_", lon_data[lon_idx], ".Rds"
            )

            if (file.exists(file_path)) {
                print("File exists")
                next
            }
            process_grid_point(
                lat_data[lat_idx], lon_data[lon_idx],
                z_data[lon_idx, lat_idx, ], time, change_points, smooth_mod
            )
        }
    }
}


# Sys.setenv(OMP_NUM_THREADS = "1")
# Sys.setenv(MKL_NUM_THREADS = "1")
# Sys.setenv(OPENBLAS_NUM_THREADS = "1")
# Sys.setenv(BLIS_NUM_THREADS = "1")
# Sys.setenv(NUMEXPR_NUM_THREADS = "1")
# Sys.setenv(R_THREADS = "1")

process_wrapper(CHUNK_NO)

# num_chunks <- 48
# num_cores <- 48
# mclapply(1:num_chunks, process_wrapper, mc.cores = num_cores)

# parallelization didnt work bc of oom errors

# for (i in seq_len(nrow(grid_points))) {
#     process_wrapper(i)
# }
