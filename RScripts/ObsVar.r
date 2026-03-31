source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")


args <- commandArgs(trailingOnly = TRUE)
CHUNK_DIR <- args[1]
CHUNK_NO <- as.integer(args[2]) +1
SEAS <- as.logical(args[3])
WRS  <- as.logical(args[4])
CP_NO_SEAS <- c("1948-09", "1957-08", "1979-03", "1998-08", "2009-01")

OUT_DIR <- paste0(CHUNK_DIR, "mods")
if (SEAS) {
    OUT_DIR <- paste0(OUT_DIR, "Seas")
    change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

} else {
  change_points <- as.Date(paste0(CP_NO_SEAS, "-01"), format = "%Y-%m-%d")

}
if (WRS) {
    OUT_DIR <- paste0(OUT_DIR, "WR")
    wr_min <- readRDS("/home/schoelleh96/wp22a/data/wrnames.rds")
    wr_min$date <- as.Date(wr_min$date)
} else {
    wr_min = NULL
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(OUT_DIR, "/logs/"), showWarnings = FALSE, recursive = TRUE)
# Change points
# Convert to Date objects for easier comparison

# grid-point-wise function
process_grid_point <- function(lat, lon, z, time, cps, wrs) {
    print("pull data together")

    z_df <- data.frame(time = time, z = z)
    print("Data frame created")
    z_df$year <- as.integer(format(z_df$time, "%Y"))
    z_df$doy <- as.integer(format(z_df$time, "%j"))
    print("now combine")

    z_df$sin_doy <- sin(2 * pi * z_df$doy / 365)
    z_df$cos_doy <- cos(2 * pi * z_df$doy / 365)
    if (SEAS){
        z_df$sin2_doy <- sin(4 * pi * z_df$doy / 365)
        z_df$cos2_doy <- cos(4 * pi * z_df$doy / 365)
        z_df$sin3_doy <- sin(6 * pi * z_df$doy / 365)
        z_df$cos3_doy <- cos(6 * pi * z_df$doy / 365)
    }
    if (WRS) {
        z_df <- dplyr::inner_join(
            z_df,
            wrs,
            by = c("time" = "date")
        )
    }
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

    print("Estimate linear model")
    rhs <- c(
    "segment",
    "segment:year",
    "segment:sin_doy + segment:cos_doy"
    )

    # add seasonal harmonics 2 and 3 only if SEAS flag is TRUE
    if (SEAS) {
    rhs <- c(rhs,
            "segment:sin2_doy + segment:cos2_doy",
            "segment:sin3_doy + segment:cos3_doy")
    }

    # add WR predictor only if WRS flag is TRUE
    if (WRS) {
        rhs <- c(rhs, "wrname")   # or "segment:wrname" if you want regime effects by segment
    }

    fml <- as.formula(
        paste0("log_variance ~ ", paste(rhs, collapse = " + ")) # , " - 1")
    )

    lmod_seas <- nlme::gls(
        fml,
        data = z_df,
        correlation = corAR1(form = ~ day_no | segment),
        weights = varIdent(form = ~ 1 | segment),
          control = glsControl(
    maxIter = 500,
    msMaxIter = 500,
    tolerance = 1e-5,
    msTol = 1e-5,
    msVerbose = TRUE
  )
    )
    lmod_seas$data <- z_df
    print("Saving linear model")
    saveRDS(lmod_seas, file = paste0(
        OUT_DIR, "/lm", lat, "_", lon,
        ".Rds"
    ))

    print("Linear model saved")

}

# Define the wrapper for parallel execution
process_wrapper <- function(chunk_idx) # , grid_points, lat_data, lon_data, z_data,
# time_data, change_points)
{
    # redirect output into separate files
    output_file <- paste0(OUT_DIR, "/logs/output_process_", chunk_idx, ".log")
    error_file <- paste0(OUT_DIR, "/logs/error_process_", chunk_idx, ".log")
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
    if (WRS) {
    years <- as.integer(format(time, "%Y"))
    year_idx <- which(years >= YEAR_BOUND_WR[1] & years <= YEAR_BOUND_WR[2])

    time <- time[year_idx]
    z_data <- z_data[, , year_idx, drop = FALSE]
    }
    nc_close(nc_file)
    print("Data read")
    
        for (i in seq_len(nrow(grid_points))) {
            lon_idx <- grid_points[i, 1]
            lat_idx <- grid_points[i, 2]
            print(paste0("Processing grid point ", lat_idx, ", ", lon_idx))
            print(paste0(
                "Lat: ", lat_data[lat_idx],
                ", Lon: ", lon_data[lon_idx]
            ))
            file_path <- paste0(
                OUT_DIR, "/lm",
                lat_data[lat_idx], "_", lon_data[lon_idx], ".Rds"
            )

            # if (file.exists(file_path)) {
            #     print("File exists")
            #     next
            # }
            process_grid_point(
                lat_data[lat_idx], lon_data[lon_idx],
                z_data[lon_idx, lat_idx, ], time, change_points, wr_min
            )
        }
    
}

process_wrapper(CHUNK_NO)
