# Set library path
.libPaths("/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib")

# List of required packages
required_packages <- c(
    "ncdf4", "tidync", "Matrix", "mgcv", "tidyr",
    "lubridate", "dplyr"
)

options(repos = c(CRAN = "https://cran.r-project.org"))

# Install missing packages and load libraries
missing_packages <- required_packages[
    !(required_packages %in% installed.packages()[, "Package"])
]
if (length(missing_packages) > 0) {
    install.packages(missing_packages,
        lib = "/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib"
    )
}

# Load the required libraries
lapply(required_packages, require, character.only = TRUE)


# Constants
FILE <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/data.nc"
CP <- c("1947-12", "1955-01", "1979-03", "1998-05")
OUT_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/mods/"

# Read in the data
dims <- tidync(FILE)$dimension
time_steps <- filter(dims, name == "time") %>% pull(length)
lat_steps <- filter(dims, name == "latitude") %>% pull(length)
lon_steps <- filter(dims, name == "longitude") %>% pull(length)

nc_file <- nc_open(FILE)

z_data <- ncvar_get(nc_file, "z")

# Get coordinates
lon_data <- ncvar_get(nc_file, "longitude", count = lon_steps)
lat_data <- ncvar_get(nc_file, "latitude", count = lat_steps)
grid_points <- expand.grid(lon_data, lat_data)

time_data <- ncvar_get(nc_file, "valid_time", count = time_steps)

# Convert time
time_origin <- sub(
    "seconds since ", "",
    ncatt_get(nc_file, "valid_time", "units")$value
)
time <- as.Date(as.POSIXct(time_data, origin = time_origin, tz = "UTC"))

nc_close(nc_file)

# Change points
# Convert to Date objects for easier comparison
change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

# grid-point-wise function
process_grid_point <- function(lat, lon, z, time, cps) {
    # pull data together

    z_df <- data.frame(time = time, z = z)
    z_df$year <- as.integer(format(z_df$time, "%Y"))
    z_df$doy <- as.integer(format(z_df$time, "%j"))

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
        right = TRUE
    )

    # Estimate linear model
    lmod_seas <- lm(log_variance ~ segment * year * (sin_doy + cos_doy) - 1,
        data = observed_variance
    )

    saveRDS(lmod_seas, file = paste0(
        OUT_DIR, "lm", lat, "_", lon,
        ".Rds"
    ))

    print("Linear model saved")
    # Estimate smooth model
    lmod_seas_smooth <- gam(
        log_variance ~ s(doy, by = segment, bs = "cc") +
            s(year, by = segment, bs = "cr") +
            ti(year, doy, by = segment, bs = c("cr", "cc")),
        data = observed_variance,
        method = "REML"
    )
    saveRDS(lmod_seas_smooth, file = paste0(
        OUT_DIR, "slm", lat, "_", lon,
        ".Rds"
    ))
    print("Smooth model saved")
}

# Define the wrapper for parallel execution
process_wrapper <- function(grid_idx) {
    lon_idx <- grid_points[grid_idx, 1]
    lat_idx <- grid_points[grid_idx, 2]

    process_grid_point(
        lat_data[lat_idx], lon_data[lon_idx],
        z_data[lon_idx, lat_idx, , ], time, change_points
    )
}

num_cores <- detectCores() - 1
mclapply(seq_len(nrow(grid_points)), processWrapper, mc.cores = numCores)
