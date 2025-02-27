# Set library path
# List of required packages
required_packages <- c(
    "ncdf4", "Matrix", "mgcv", "tidyr",
    "lubridate", "dplyr", "parallel"
)

# Install missing packages and load libraries
missing_packages <- required_packages[
    !(required_packages %in% installed.packages()[, "Package"])
]
if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

# Load the required libraries
lapply(required_packages, require, character.only = TRUE)

print("All packages loaded")

# Constants
FILE <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/data.nc"
CP <- c("1948-05", "1958-05", "1979-02", "1998-08", "2009-07")
OUT_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/mods/"
CHUNK_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/chunks/"

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
    z_df <- z_df %>%
        pivot_longer(
            cols = starts_with("z"), # Select all columns that are part of 'z'
            names_to = "member",
            values_to = "value"
        ) %>%
        group_by(time) %>% # Group by time to calculate differences per row
        mutate(
            value = if_else(
                member != "z.1", # If not 'z1',
                value - value[member == "z.1"], # Subtract 'z1' value
                value # Keep 'z1' values unchanged
            )
        ) %>%
        ungroup()

    z_df$sin_doy <- sin(2 * pi * z_df$doy / 365)
    z_df$cos_doy <- cos(2 * pi * z_df$doy / 365)

    print("calculate observed variance")
    z_df <- z_df %>% filter(member != "z.1")
    observed_variance <- z_df %>%
        group_by(time) %>%
        summarize(observed_variance = var(value, na.rm = TRUE)) %>%
        ungroup()

    observed_variance <- observed_variance %>%
        mutate(
            year = year(time), # Extract year
            doy = yday(time), # Extract day of the year (1-366)
            log_variance = log(observed_variance), # Log-transform
            sin_doy = sin(2 * pi * doy / 365), # Sinusoidal component
            cos_doy = cos(2 * pi * doy / 365), # Cosine component
        )
    observed_variance$segment <- cut(as.Date(observed_variance$time),
        breaks = c(as.Date(c(
            observed_variance$time[1],
            cps,
            observed_variance$time[nrow(observed_variance)]
        ))),
        labels = 0:(length(cps)) + 1,
        include.lowest = TRUE,
        right = FALSE
    )
    if (smooth_mod) {
        print("Estimate smooth model")
        # Estimate smooth model
        lmod_seas_smooth <- gam(
            log_variance ~ s(doy, by = segment, bs = "cc") +
                s(year, by = segment, bs = "cr") +
                ti(year, doy, by = segment, bs = c("cr", "cc")),
            data = observed_variance,
            method = "REML",
            keepData = TRUE
        )
        print("Saving smooth model")
        saveRDS(lmod_seas_smooth, file = paste0(
            OUT_DIR, "slm", lat, "_", lon,
            ".Rds"
        ))
        print("Smooth model saved")
    } else {
        print("Estimate linear model")
        lmod_seas <- lm(
            log_variance ~ segment +
                segment:year + segment:sin_doy + segment:cos_doy +
                segment:year:sin_doy + segment:year:cos_doy - 1,
            data = observed_variance
        )
        # lmod_seas$residuals <- NULL # Remove residuals
        # lmod_seas$fitted.values <- NULL # Remove fitted values
        # lmod_seas$model <- NULL # Remove the model frame
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
    for (smooth_mod in c(FALSE, TRUE)) {
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
                z_data[lon_idx, lat_idx, , ], time, change_points, smooth_mod
            )
        }
    }
}


Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(R_THREADS = "1")

num_chunks <- 48
num_cores <- 48
mclapply(1:num_chunks, process_wrapper, mc.cores = num_cores)

# parallelization didnt work bc of oom errors

# for (i in seq_len(nrow(grid_points))) {
#     process_wrapper(i)
# }
