library(dplyr)
library(purrr)
library(tidyr)
library(broom)
library(stringr)
library(ncdf4)
library(lubridate)

load_lm_objects <- function(lm_dir = LM_DIR) {
    lm_files <- list.files(lm_dir, pattern = "^lm.*\\.Rds", full.names = TRUE)
    lm_files <- lm_files

    # Extract lat/lon from filenames
    coords <- stringr::str_extract_all(lm_files, "-?\\d+\\.?\\d*",
        simplify = TRUE
    )
    lats <- as.numeric(coords[, 1])
    lons <- as.numeric(coords[, 2])

    # Load models and extract diagnostics
    lm_data <- purrr::map_dfr(seq_along(lm_files), function(i) {
        lm_obj <- readRDS(lm_files[i])
        segment_range <- range(as.numeric(unique(lm_obj$model$segment)))
        # Model-level diagnostics (broom::glance)
        model_summary <- broom::glance(lm_obj)

        # Coefficients and p-values (broom::tidy)
        coef_table <- broom::tidy(lm_obj) %>%
            dplyr::select(term, estimate, p.value) %>%
            tidyr::pivot_wider(
                names_from = term,
                values_from = c(estimate, p.value),
                names_glue = "{term}_{.value}"
            )

        new_data <- tibble::tibble(
            year = YEAR_BOUND,
            segment = as.factor(segment_range),
            sin_doy = 0,
            cos_doy = 0
        )
        fitted_values <- predict(lm_obj, newdata = new_data)

        # Compute the difference
        fitted_diff <- -diff(fitted_values)

        # Compute segment jumps at change points (CPs)
        cp_years <- as.integer(sub("-.*", "", CP))
        segment_jumps <- tibble::tibble()
        j <- 1
        for (cp_year in cp_years) {
            # Generate new data for adjacent segments at the given CP year
            new_data_CP <- tibble::tibble(
                year = cp_year,
                segment = as.factor(c(j, j + 1)),
                sin_doy = 0,
                cos_doy = 0
            )
            fitted_CP <- predict(lm_obj, newdata = new_data_CP)
            jump_name <- paste0("jump_", cp_year, "_", j, "_to_", j + 1)
            # convention: positive jump means improvement
            jump_value <- unname(fitted_CP[1] - fitted_CP[2])
            # Use add_column() to add the new jump value
            if (ncol(segment_jumps) == 0) {
                segment_jumps <- tibble::tibble(!!jump_name := jump_value)
            } else {
                segment_jumps <- segment_jumps %>%
                    tibble::add_column(!!jump_name := jump_value)
            }
            j <- j + 1
        }

        # Combine all extracted data
        tibble::tibble(
            lat = lats[i],
            lon = lons[i],
            obs_var = var(lm_obj$model$log_variance),
            sigma = model_summary$sigma, # Residual standard error
            aic = model_summary$AIC,
            bic = model_summary$BIC,
            f_statistic = unname(model_summary$statistic[1]),
            f_p_value = unname(model_summary$p.value[1]),
            fit1940 = unname(fitted_values[1]),
            fit2024 = unname(fitted_values[2]),
            fitted_diff = unname(fitted_diff)
        ) %>%
            mutate(
                RSS = model_summary$sigma^2 * (TIME_STEPS - N_COEFFS - 1),
                TSS = obs_var * (TIME_STEPS - 1),
                r_squared = 1 - (RSS / TSS),
                adj_r_squared = 1 - ((1 - r_squared) * (TIME_STEPS - 1) /
                    (TIME_STEPS - N_COEFFS - 1))
            ) %>%
            dplyr::bind_cols(coef_table) %>%
            dplyr::bind_cols(segment_jumps)
    })

    return(lm_data)
}

fdr_adjust_pvalues <- function(df, method = "fdr") {
    # Identify p.value columns
    p_cols <- names(df)[grepl("p\\.value$", names(df))]

    # Adjust p-values using the specified method
    df <- df %>%
        mutate(across(all_of(p_cols),
            ~ p.adjust(., method = method),
            .names = "{.col}_adj"
        ))

    # Identify pairs of columns for seasonal magnitude calculation
    estimate_cols <- names(df)[grepl("estimate$", names(df))]
    sin_cols <- estimate_cols[grepl("sin_doy", estimate_cols)]
    cos_cols <- estimate_cols[grepl("cos_doy", estimate_cols)]

    # Create seasonal magnitude columns
    for (sin_col in sin_cols) {
        base_name <- sub(":sin_doy_estimate", "", sin_col)
        cos_col <- paste0(base_name, ":cos_doy_estimate")

        if (cos_col %in% estimate_cols) {
            new_col <- paste0(base_name, "_seas")
            df[[new_col]] <- sqrt(df[[sin_col]]^2 + df[[cos_col]]^2)
        }
    }

    # Identify pairs of adjusted p-value columns
    p_adj_cols <- names(df)[grepl("p\\.value_adj$", names(df))]
    sin_p_cols <- p_adj_cols[grepl("sin_doy", p_adj_cols)]
    cos_p_cols <- p_adj_cols[grepl("cos_doy", p_adj_cols)]

    # Create columns with minimum adjusted p-value
    for (sin_p_col in sin_p_cols) {
        base_name <- sub(":sin_doy_p\\.value_adj", "", sin_p_col)
        cos_p_col <- paste0(base_name, ":cos_doy_p.value_adj")

        if (cos_p_col %in% p_adj_cols) {
            new_col <- paste0(base_name, "_p.value_min")
            df[[new_col]] <- pmin(df[[sin_p_col]], df[[cos_p_col]],
                na.rm = TRUE
            )
        }
    }

    return(df)
}

add_obs_var <- function(data, time_steps, coeffs) {
    calculate_variance <- function(r_squared, sigma, n, p) {
        # Residual sum of squares (RSS)
        rss <- sigma^2 * (n - p)

        # Total sum of squares (TSS)
        tss <- rss / (1 - r_squared)

        # Total variance is TSS / (n - 1)
        variance_original_data <- tss / (n - 1)

        return(variance_original_data)
    }

    # Add a new column with the total variance for each model
    data <- data %>%
        mutate(variance_original_data = mapply(calculate_variance, r_squared,
            sigma,
            MoreArgs = list(time_steps, coeffs)
        ))
    return(data)
}

preprocess_pval <- function(data,
                            pattern = "^segment[1-5].*_p\\.value_adj$",
                            alpha = ALPHA) {
    # Step 1: Extract relevant p-value columns
    pval_cols <- grep(pattern, colnames(data), value = TRUE)

    # Step 2: Count significant values for each column
    count_data <- set_names(pval_cols) |>
        map_dbl(~ sum(data[[.x]] < alpha, na.rm = TRUE))

    # Step 3: Extract metadata from column names
    col_metadata <- tibble(original_name = pval_cols) |>
        mutate(
            Segment = str_extract(original_name, "(?<=^segment)[1-5]"),
            Metric = str_replace(original_name, "^segment[1-5]", "seg") |>
                str_remove("_p\\.value_adj$")
        )

    # Step 4: Create summary data
    summary_df <- col_metadata |>
        mutate(Count = count_data[original_name]) |>
        select(-original_name) |>
        mutate(Segment = as.numeric(Segment))

    # Step 5: Reshape data for heatmap
    heatmap_data <- summary_df |>
        pivot_wider(names_from = Metric, values_from = Count) |>
        pivot_longer(-Segment, names_to = "Metric", values_to = "Count")

    return(heatmap_data)
}

#---------------------- Dataset loading functions ----------------------------

calculate_monthly_averages <- function(file_name) {
    #' Calculate monthly averages of variable "z" after spatial averaging.
    #'
    #' The function opens a netCDF file, reads variable "z", averages over
    #' the spatial dimensions ("latitude" and "longitude"), extracts the month
    #' from the time coordinate, and then computes the monthly mean.
    #'
    #' @param file_name Character. NetCDF file name.
    #' @return A data frame with columns: month and avg_z.
    nc <- nc_open(file_name)
    on.exit(nc_close(nc))

    z_data <- ncvar_get(nc, "z")

    lat <- ncvar_get(nc, "latitude")
    lon <- ncvar_get(nc, "longitude")
    time_data <- ncvar_get(nc, "valid_time")
    time_origin <- sub(
        "seconds since ", "",
        ncatt_get(nc, "valid_time", "units")$value
    )
    time <- as.Date(as.POSIXct(time_data,
        origin = time_origin,
        tz = "UTC"
    ))
    months <- as.numeric(format(time, "%m"))

    z_info <- nc$var[["z"]]
    z_dim_names <- sapply(z_info$dim, function(x) x$name)
    time_dim_index <- which(z_dim_names %in% c("time", "valid_time"))
    spatial_dim_indices <- setdiff(seq_along(z_dim_names), time_dim_index)

    spatial_avg <- apply(z_data, time_dim_index, mean, na.rm = TRUE)
    df <- data.frame(month = months, avg_z = spatial_avg)
    monthly_df <- df %>%
        group_by(month) %>%
        summarise(avg_z = mean(avg_z, na.rm = TRUE)) %>%
        ungroup() %>%
        arrange(month)

    return(monthly_df)
}

combine_datasets <- function(datasets, years) {
    #' Combine a list of monthly datasets into one data frame.
    #'
    #' Each dataset (a data frame with month and avg_z columns) is tagged with
    #' its corresponding year.
    #'
    #' @param datasets List of data frames.
    #' @param years Integer vector corresponding to each data frame.
    #' @return A combined data frame with columns: year, month, avg_z.

    combined <- mapply(function(df, yr) {
        df$year <- yr
        df
    }, datasets, years, SIMPLIFY = FALSE)
    combined_df <- do.call(rbind, combined)
    combined_df <- combined_df[order(
        combined_df$year,
        combined_df$month
    ), ]
    rownames(combined_df) <- NULL
    return(combined_df)
}

load_model_data <- function(lon, lat, lm_dir, conf_level = 0.95) {
    # Find indices for which the model file exists
    valid_indices <- which(purrr::map2_lgl(lon, lat, ~ {
        file_name <- file.path(lm_dir, paste0("lm", .y, "_", .x, ".Rds"))
        file.exists(file_name)
    }))

    # Use the first valid grid point to extract the common predictor data.
    first_idx <- valid_indices[1]
    file_first <- file.path(
        lm_dir,
        paste0("lm", lat[first_idx], "_", lon[first_idx], ".Rds")
    )
    master_mod <- readRDS(file_first)

    # Determine the start date as January 1 of the earliest year in the data.
    min_year <- min(master_mod$model$year)
    start_date <- as.Date(paste0(min_year, "-01-01"))

    # Arrange the model data in chronological order across all years.
    # (Here we assume that sorting by 'year' gives the correct order.)
    df_master_unique <- master_mod$model %>%
        arrange(year) %>%
        mutate(
            global_day = row_number(), # Sequential day index over the full dataset
            date = start_date + (global_day - 1), # Create dates starting from the earliest Jan 1
            doy = lubridate::yday(date) # Day-of-year from the constructed date
        ) %>%
        select(year, segment, sin_doy, cos_doy, doy, date)

    # For each grid point, extract observed and fitted values and create new columns.
    additional_cols <- purrr::map2_dfc(lon, lat, function(lon_i, lat_i) {
        file_name <- file.path(
            lm_dir,
            paste0("lm", lat_i, "_", lon_i, ".Rds")
        )

        mod <- readRDS(file_name)
        df_mod <- mod$model

        # Extract the observed response (log_variance) and fitted values.
        observed_vals <- df_mod$log_variance
        fitted_vals <- as.vector(fitted(mod))
        residuals <- observed_vals - fitted_vals

        # Calculate confidence intervals for fitted values.
        ci <- predict(mod, interval = "confidence", level = conf_level)
        lower_ci <- ci[, "lwr"]
        upper_ci <- ci[, "upr"]

        # Create column names that include the grid point information.
        obs_col_name <- paste0("observed_", lat_i, "_", lon_i)
        fit_col_name <- paste0("fitted_", lat_i, "_", lon_i)
        lwr_col_name <- paste0("lwr_", lat_i, "_", lon_i)
        upr_col_name <- paste0("upr_", lat_i, "_", lon_i)
        res_col_name <- paste0("res_", lat_i, "_", lon_i)

        tibble(
            !!obs_col_name := observed_vals,
            !!fit_col_name := fitted_vals,
            !!lwr_col_name := lower_ci,
            !!upr_col_name := upper_ci,
            !!res_col_name := residuals
        )
    })

    # Combine the master predictor data with the additional observed/fitted columns.
    final_df <- dplyr::bind_cols(df_master_unique, additional_cols)

    return(final_df)
}


aggregate_monthly <- function(df) {
    # Create a year_month column (e.g. "2020-03")
    df <- df %>%
        mutate(date = format(date, "%Y-%m-01"))

    # Identify grid points from observed columns (e.g. "observed_40_-60")
    obs_cols <- grep("^observed_", names(df), value = TRUE)
    grid_points <- sub("^observed_", "", obs_cols)

    # For each grid point, compute the monthly aggregates
    aggregated_list <- lapply(grid_points, function(grid) {
        # Define the column names for this grid
        obs_col <- paste0("observed_", grid)
        fitted_col <- paste0("fitted_", grid)
        lwr_col <- paste0("lwr_", grid)
        upr_col <- paste0("upr_", grid)
        res_col <- paste0("res_", grid)

        df %>%
            group_by(year, segment, date) %>%
            summarise(
                !!paste0("obs_median_", grid) := median(.data[[obs_col]],
                    na.rm = TRUE
                ),
                !!paste0("obs_q1_", grid) := quantile(.data[[obs_col]],
                    probs = 0.25, na.rm = TRUE
                ),
                !!paste0("obs_q3_", grid) := quantile(.data[[obs_col]],
                    probs = 0.75, na.rm = TRUE
                ),
                !!paste0("fitted_", grid) := median(.data[[fitted_col]],
                    na.rm = TRUE
                ),
                !!paste0("lwr_", grid) := median(.data[[lwr_col]],
                    na.rm = TRUE
                ),
                !!paste0("upr_", grid) := median(.data[[upr_col]],
                    na.rm = TRUE
                ),
                !!paste0("res_rmse_", grid) := sqrt(mean((.data[[res_col]])^2,
                    na.rm = TRUE
                ))
            ) %>%
            ungroup()
    })

    # Merge the aggregated data for each grid point by year, segment, and year_month
    aggregated_df <- reduce(aggregated_list,
        full_join,
        by = c("year", "segment", "date")
    )

    aggregated_df <- aggregated_df %>% arrange(year, segment, date)

    return(aggregated_df)
}

tibble_to_long <- function(df) {
    # Pivot the grid-specific columns into long format.
    # The regex below captures:
    #   - the type (which might be "observed", "obs_median", "obs_q1", "obs_q3",
    #     "fitted", "lwr", "upr", "res", or "res_rmse")
    #   - the latitude and longitude (assumed to be integers, possibly negative)
    df_long <- df %>%
        pivot_longer(
            cols = matches("^(observed|obs_median|obs_q1|obs_q3|fitted|lwr|upr|res|res_rmse)_-?\\d+_-?\\d+$"),
            names_to = c("type", "lat", "lon"),
            names_pattern = "^(observed|obs_median|obs_q1|obs_q3|fitted|lwr|upr|res|res_rmse)_(-?\\d+)_(-?\\d+)$",
            values_to = "value"
        ) %>%
        # Standardize the type names: for daily data, map "observed" -> "obs_median" and "res" -> "res_rmse"
        mutate(
            type = case_when(
                type == "observed" ~ "obs_median",
                type == "res" ~ "res_rmse",
                TRUE ~ type
            ),
            lat = as.numeric(lat),
            lon = as.numeric(lon)
        ) %>%
        # Select only the required columns: date, lon, lat, value, and type.
        select(date, lon, lat, value, type)

    return(df_long)
}
