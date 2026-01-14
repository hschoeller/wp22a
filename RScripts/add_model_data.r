# add_model_data.R
source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")

args <- commandArgs(trailingOnly = TRUE)
CHUNK_DIR <- paste0(args[1], "zg_chunks/")
CHUNK_NO <- as.integer(args[2]) +1
OUT_DIR <- paste0(args[1], "mods/")

# convert CP (vector of years or strings) to Date objects like original
change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

# function to create the z_df exactly like the original script
make_z_df <- function(z_vec, time_vec, cps) {
    z_df <- data.frame(time = time_vec, z = z_vec)
    z_df$year <- as.integer(format(z_df$time, "%Y"))
    z_df$doy <- as.integer(format(z_df$time, "%j"))
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

    # uses dplyr verbs in your environment (assumed available from sourced scripts)
    z_df <- z_df %>%
        group_by(segment) %>%
        mutate(day_no = row_number()) %>%
        ungroup()

    return(z_df)
}

# Attach data to model object in-place
attach_data_to_model <- function(mod, z_df) {
    # put z_df on top-level object
    mod$data <- z_df

    # if object contains an lme component (gamm returns list with $lme)
    if (!is.null(mod$lme)) {
        mod$lme$data <- z_df
    }

    # if object contains a gam component (gamm returns $gam)
    if (!is.null(mod$gam)) {
        # many mgcv::gam objects store model frame in $model
        mod$gam$model <- z_df
        # also, if $gam$call$data exists, replace by the actual data (not necessary but helpful)
        if (!is.null(mod$gam$call) && !is.null(mod$gam$call$data)) {
            mod$gam$call$data <- z_df
        }
    }

    # For lme/gls objects directly, ensure the internal call$data is left as-is but data slot exists.
    # Return the modified object
    return(mod)
}

process_grid_point <- function(lat, lon, z, time, cps) {
    # Build file names for both possible prefixes (lm/slm)
    lm_file <- paste0(OUT_DIR, "lm", lat, "_", lon, ".Rds")

    # Create z_df once (used for either model type if present)
    z_df <- make_z_df(z, time, cps)

    # If linear model exists and lacks embedded data, attach and save
    if (file.exists(lm_file)) {
        print(paste("Processing file:", lm_file))
        mod <- readRDS(lm_file)
        has_data <- !is.null(mod$data)
        if (!has_data) {
            mod <- attach_data_to_model(mod, z_df)
            saveRDS(mod, file = lm_file)
            message("Attached data and saved: ", lm_file)
        } else {
            message("Data already present in: ", lm_file)
        }
    }
}

process_wrapper <- function(chunk_idx) {
    output_file <- paste0(OUT_DIR, "logs/output_process_", chunk_idx, ".log")
    error_file <- paste0(OUT_DIR, "logs/error_process_", chunk_idx, ".log")
    output_con <- file(output_file, open = "wt")
    error_con <- file(error_file, open = "wt")
    sink(output_con)
    sink(error_con, type = "message")
    on.exit({
        sink()
        sink(type = "message")
        close(output_con)
        close(error_con)
    })

    chunk_file <- paste0(CHUNK_DIR, sprintf("chunk_%02d.nc", chunk_idx))

    message("Reading chunk: ", chunk_file)
    nc_file <- nc_open(chunk_file)
    z_data <- ncvar_get(nc_file, "z", collapse_degen = FALSE)
    lon_data <- ncvar_get(nc_file, "longitude")
    lat_data <- ncvar_get(nc_file, "latitude")
    grid_points <- expand.grid(1:length(lon_data), 1:length(lat_data))
    time_data <- ncvar_get(nc_file, "valid_time")

    time_origin <- sub("seconds since ", "", ncatt_get(nc_file, "valid_time", "units")$value)
    time <- as.Date(as.POSIXct(time_data, origin = time_origin, tz = "UTC"))

    nc_close(nc_file)

    # iterate grid points (same loop structure)
    for (i in seq_len(nrow(grid_points))) {
        lon_idx <- grid_points[i, 1]
        lat_idx <- grid_points[i, 2]

        message(
            "Processing grid point ", lat_idx, ", ", lon_idx,
            "  Lat: ", lat_data[lat_idx], ", Lon: ", lon_data[lon_idx]
        )

        # only proceed if at least one of the model files exists
        file_lm <- paste0(OUT_DIR, "lm", lat_data[lat_idx], "_", lon_data[lon_idx], ".Rds")
        if (!file.exists(file_lm)) {
            next
        }

        # Extract the z time series for this grid point (note ordering used in original)
        z_vec <- z_data[lon_idx, lat_idx, ]

        process_grid_point(lat_data[lat_idx], lon_data[lon_idx], z_vec, time, change_points)
    }
}

process_wrapper(CHUNK_NO)


# Sys.setenv(OMP_NUM_THREADS = "1")
# Sys.setenv(MKL_NUM_THREADS = "1")
# Sys.setenv(OPENBLAS_NUM_THREADS = "1")
# Sys.setenv(BLIS_NUM_THREADS = "1")
# Sys.setenv(NUMEXPR_NUM_THREADS = "1")
# Sys.setenv(R_THREADS = "1")

# num_chunks <- 48
# num_cores <- 48
# mclapply(1:num_chunks, process_wrapper, mc.cores = num_cores)
