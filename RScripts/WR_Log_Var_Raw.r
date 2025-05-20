source("RScripts/config.r")
source("RScripts/data_functions.r")

# Constants
FILE <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/data.nc"
OUT_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"
CHUNK_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/chunks/"

# Change points
# Convert to Date objects for easier comparison
change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

# WRs
result <- wrera(
    start = "19500111_00",
    end = "20250113_21",
    hours = c("12"),
    tformat = "string",
    setup = "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019",
    dataset = "era5",
    basepath = paste0(
        OUT_DIR,
        "../WR_read_example_package/wr_era5_update_1950_latwgt/"
    )
)

wr_df <- result$data$LC
wr_df$date <- as.Date(wr_df$time)

print("WRs loaded")

# grid-point-wise function
process_grid_point <- function(lat, lon, z, time, cps, wr_df) {
    print("pull data together")

    z_df <- data.frame(time = time, z = z)
    print("Data frame created")
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

    print("calculate observed variance")
    z_df <- z_df %>% filter(member != "z.1")
    observed_variance <- z_df %>%
        group_by(time) %>%
        summarize(observed_variance = var(value, na.rm = TRUE)) %>%
        ungroup()

    observed_variance <- observed_variance %>%
        mutate(
            log_variance = log(observed_variance), # Log-transform
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
    observed_variance <- observed_variance %>%
        inner_join(wr_df %>% select(date, wrindex, wrname),
            by = c("time" = "date")
        )

    print("Calculating means and pvals")

    overall_mean <- mean(observed_variance$log_variance, na.rm = TRUE)
    results <- data.frame()
    # 1. Overall analysis (across whole time series)
    for (wr in unique(observed_variance$wrindex)) {
        wr_dates <- observed_variance$time[
            observed_variance$wrindex == wr
        ]

        if (length(wr_dates) == 0) next

        result <- permute_test(observed_variance, wr_dates, N_PERM)

        results <- rbind(results, data.frame(
            segment = -999,
            wrindex = wr,
            mean = result$mean,
            p_val = result$p_val
        ))
    }
    # 2. Segment-wise analysis
    for (seg in unique(observed_variance$segment)) {
        seg_data <- observed_variance[observed_variance$segment == seg, ]
        seg_mean <- mean(seg_data$log_variance, na.rm = TRUE)

        results <- rbind(results, data.frame(
            segment = seg,
            wrindex = -999,
            mean = seg_mean,
            p_val = NA
        ))

        for (wr in unique(seg_data$wrindex)) {
            category_data_seg <- seg_data[seg_data$wrindex == wr, ]
            category_dates_seg <- category_data_seg$time

            if (length(category_dates_seg) == 0) next

            result <- permute_test(seg_data, category_dates_seg, N_PERM)

            results <- rbind(results, data.frame(
                segment = seg,
                wrindex = wr,
                mean = result$mean,
                p_val = result$p_val
            ))
        }
    }

    results <- rbind(data.frame(
        segment = -999,
        wrindex = -999,
        mean = overall_mean,
        p_val = NA
    ), results)
    results <- results %>% mutate(
        lat = lat, lon = lon,
        segment = factor(segment),
        wrindex = factor(wrindex)
    )
    results$segment <- factor(results$segment)
    results$wrindex <- factor(results$wrindex)
    print("Done")
    return(results)
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
    results <- data.frame()

    for (i in seq_len(nrow(grid_points))) {
        lon_idx <- grid_points[i, 1]
        lat_idx <- grid_points[i, 2]
        print(paste0("Processing grid point ", lat_idx, ", ", lon_idx))
        print(paste0(
            "Lat: ", lat_data[lat_idx],
            ", Lon: ", lon_data[lon_idx]
        ))
        results <- bind_rows(
            results,
            process_grid_point(
                lat_data[lat_idx], lon_data[lon_idx],
                z_data[lon_idx, lat_idx, , ], time, change_points, wr_df
            )
        )
    }
    return(results)
}


Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(R_THREADS = "1")

num_chunks <- 48
num_cores <- 48
grid_result <- mclapply(1:num_chunks, process_wrapper, mc.cores = num_cores)
results <- bind_rows(grid_result)
saveRDS(results, paste0(OUT_DIR, "Raw_WR_Comp.RDS"))
