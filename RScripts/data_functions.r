library(dplyr)
library(purrr)
library(tidyr)
library(broom)
library(stringr)
library(ncdf4)
library(lubridate)
library(data.table)
library(parallel)

#---------------------- Dataset loading functions ----------------------------

read_nc_variable <- function(file_path, var_name) {
    # Open the NetCDF file
    nc <- ncdf4::nc_open(file_path)
    on.exit(ncdf4::nc_close(nc))

    # Get longitude and latitude (optional, depending on variable grid)
    lon <- ncdf4::ncvar_get(nc, "longitude")
    lat <- ncdf4::ncvar_get(nc, "latitude")

    # Read the raw time data and convert to Date
    time_data <- ncdf4::ncvar_get(nc, "valid_time")
    time_units <- ncdf4::ncatt_get(nc, "valid_time", "units")$value
    # Extract origin string, e.g., "seconds since YYYY-MM-DD hh:mm:ss"
    origin <- sub("seconds since ", "", time_units)
    time_dates <- as.Date(
        as.POSIXct(time_data,
            origin = origin,
            tz = "UTC"
        )
    )

    # Read the requested variable data
    data_values <- ncdf4::ncvar_get(nc, var_name)

    # Return a list with data and time
    list(
        data = data_values,
        time = time_dates,
        lon = lon,
        lat = lat
    )
}

#' Load Linear Model Objects and Extract Diagnostics
#'
#' This function loads linear model objects from a specified directory, extracts
#' relevant diagnostics, and computes additional metrics such as fitted differences
#' and segment jumps at change points.
#'
#' @param lm_dir Character. The directory containing the linear model `.Rds` files.
#'   Defaults to the global variable `LM_DIR`.
#'
#' @return A tibble containing the following information for each linear model:
#'   - `lat`: Latitude extracted from the file name.
#'   - `lon`: Longitude extracted from the file name.
#'   - `obs_var`: Observed variance of the `log_variance` in the model.
#'   - `sigma`: Residual standard error of the model.
#'   - `aic`: Akaike Information Criterion (AIC) of the model.
#'   - `bic`: Bayesian Information Criterion (BIC) of the model.
#'   - `f_statistic`: F-statistic of the model.
#'   - `f_p_value`: P-value associated with the F-statistic.
#'   - `fit1940`: Fitted value for the year 1940.
#'   - `fit2024`: Fitted value for the year 2024.
#'   - `fitted_diff`: Difference between fitted values for 1940 and 2024.
#'   - `RSS`: Residual Sum of Squares.
#'   - `TSS`: Total Sum of Squares.
#'   - `r_squared`: R-squared value of the model.
#'   - `adj_r_squared`: Adjusted R-squared value of the model.
#'   - Coefficients and p-values for each term in the model.
#'   - Segment jumps at specified change points (e.g., `jump_<year>_<segment>_to_<segment>`).
#'
#' @details
#' The function performs the following steps:
#'   1. Reads `.Rds` files matching the pattern `^lm.*\\.Rds` from the specified directory.
#'   2. Extracts latitude and longitude from the file names.
#'   3. Loads each linear model object and computes diagnostics using the `broom` package.
#'   4. Predicts fitted values for specific years and computes their differences.
#'   5. Computes segment jumps at specified change points.
#'   6. Combines all extracted data into a single tibble.
#'
#' @note
#' This function relies on several global variables:
#'   - `YEAR_BOUND`: A vector specifying the years for fitted value predictions (e.g., 1940 and 2024).
#'   - `CP`: A vector of change points (e.g., years where segment jumps are computed).
#'   - `TIME_STEPS`: Total number of time steps in the model.
#'   - `N_COEFFS`: Number of coefficients in the model.
#'
#' @importFrom stringr str_extract_all
#' @importFrom purrr map_dfr
#' @importFrom broom glance tidy
#' @importFrom dplyr select mutate bind_cols
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble add_column
#'
#' @examples
#' # Example usage:
#' lm_data <- load_lm_objects("/path/to/lm/files")
#'
#' @export
load_lm_objects <- function(lm_dir = LM_DIR) {
    lm_files <- list.files(lm_dir, pattern = "^lm.*\\.Rds", full.names = FALSE)
    lm_files <- lm_files

    # Extract lat/lon from filenames
    coords <- stringr::str_extract_all(lm_files, "-?\\d+\\.?\\d*",
        simplify = TRUE
    )
    lats <- as.numeric(coords[, 1])
    lons <- as.numeric(coords[, 2])

    # Load models and extract diagnostics
    lm_data <- purrr::map_dfr(seq_along(lm_files), function(i) {
        lm_obj <- readRDS(paste0(LM_DIR, lm_files[i]))
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



# Some helper functions for data loaded as above

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

#' Load Model Data for Specified Grid Points
#'
#' This function loads and processes model data for multiple grid points,
#' combining predictor data and observed/fitted values into a single data frame.
#'
#' @param lon A numeric vector of longitudes for the grid points.
#' @param lat A numeric vector of latitudes for the grid points.
#' @param lm_dir A string specifying the directory containing the model files.
#' @param conf_level A numeric value specifying the confidence level for
#'   prediction intervals (default is 0.95).
#'
#' @return A data frame containing:
#'   - Predictor data (e.g., year, segment, sin_doy, cos_doy, doy, date).
#'   - Observed response values (`log_variance`) for each grid point.
#'   - Fitted values, residuals, and confidence intervals for each grid point.
#'
#' @details
#' The function performs the following steps:
#'   1. Identifies valid grid points for which model files exist.
#'   2. Extracts common predictor data from the first valid grid point.
#'   3. Constructs a master data frame with predictor variables arranged
#'      chronologically.
#'   4. For each grid point, extracts observed response values, fitted values,
#'      residuals, and confidence intervals.
#'   5. Combines the predictor data with the observed/fitted data into a
#'      single data frame.
#'
#' @examples
#' \dontrun{
#' lon <- c(10, 20, 30)
#' lat <- c(50, 60, 70)
#' lm_dir <- "/path/to/model/files"
#' conf_level <- 0.95
#' result <- load_model_data(lon, lat, lm_dir, conf_level)
#' head(result)
#' }
#'
#' @importFrom purrr map2_lgl map2_dfc
#' @importFrom dplyr arrange mutate select bind_cols
#' @importFrom lubridate yday
#' @importFrom tibble tibble
#' @export
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

calculate_monthly_averages <- function(file_name) {
    #' Calculate monthly averages of variable "z" after area-weighted spatial averaging.
    #'
    #' @param file_name Character. NetCDF file name.
    #' @return A data frame with columns: month and avg_z.

    nc <- nc_open(file_name)
    on.exit(nc_close(nc))

    z_data <- ncvar_get(nc, "z")
    lat <- ncvar_get(nc, "latitude")
    time_data <- ncvar_get(nc, "valid_time")

    time_origin <- sub("seconds since ", "", ncatt_get(nc, "valid_time", "units")$value)
    time <- as.Date(as.POSIXct(time_data, origin = time_origin, tz = "UTC"))
    months <- as.numeric(format(time, "%m"))

    # Get dimension info and create area weights
    z_info <- nc$var[["z"]]
    z_dim_names <- sapply(z_info$dim, function(x) x$name)
    time_dim_index <- which(z_dim_names %in% c("time", "valid_time"))

    # Create area weight array matching z_data dimensions
    lat_weights <- cos(lat * pi / 180)
    weight_dims <- dim(z_data)
    weight_array <- array(1, dim = weight_dims)

    # Apply latitude weights to appropriate dimension
    lat_dim_pos <- which(z_dim_names %in% c("latitude", "lat"))
    if (lat_dim_pos == 1) {
        weight_array <- sweep(weight_array, 1, lat_weights, "*")
    } else if (lat_dim_pos == 2) {
        weight_array <- sweep(weight_array, 2, lat_weights, "*")
    }

    # Vectorized area-weighted spatial averaging
    weighted_z <- z_data * weight_array
    valid_mask <- !is.na(z_data)

    spatial_avg <- apply(weighted_z, time_dim_index, sum, na.rm = TRUE) /
        apply(weight_array * valid_mask, time_dim_index, sum, na.rm = TRUE)

    # Calculate monthly averages
    data.frame(month = months, avg_z = spatial_avg) %>%
        group_by(month) %>%
        summarise(avg_z = mean(avg_z, na.rm = TRUE), .groups = "drop") %>%
        arrange(month)
}

# --- For gls data

add_fdr_pvalues_to_coefs <- function(models_df, method = "fdr") {
    stopifnot("coefs" %in% names(models_df))

    library(dplyr)
    library(tidyr)

    models_df <- models_df %>%
        mutate(.id = row_number())

    coef_long <- models_df %>%
        select(.id, coefs) %>%
        unnest(coefs)

    coef_long <- coef_long %>%
        group_by(term) %>%
        mutate(p.value_adj = p.adjust(p.value, method = method)) %>%
        ungroup()

    coefs_adj <- coef_long %>%
        group_by(.id) %>%
        summarise(
            coefs = list(
                data.frame(
                    term,
                    estimate,
                    std.error,
                    t.value,
                    p.value,
                    p.value_adj
                )
            ),
            .groups = "drop"
        )

    models_df %>%
        select(-coefs) %>% # <-- do NOT drop .id yet
        left_join(coefs_adj, by = ".id") %>%
        select(-.id) # <-- drop it only at the end
}


wrera <- function(start, end, hours, tformat, setup, dataset, basepath) {
    #' Load weather regime indices and life cycles from era-interim or era5
    #' for a specific period (for details about the weather regimes
    #' contact Christian Grams, christian.grams@gmx.de)
    #'
    #' @param start Character string with start date in format "YYYYMMDD_HH"
    #' @param end Character string with end date in format "YYYYMMDD_HH"
    #' @param hours Hour or list of hours (e.g., "00","03","06") to include
    #' @param tformat "string" for string format dates, "dtime" for Date objects
    #'        (Note: time columns are now always POSIXct objects regardless of this parameter)
    #' @param setup Setup string defining the weather regime data source
    #' @param dataset Either "erainterim" or "era5"
    #' @param basepath Directory path to the raw data
    #'
    #' @return List with components:
    #'   - IWR: weather regime index with fields tsince, time (POSIXct), cci and
    #'          the different weather regime projections
    #'   - MAXIWR: maximum weather regime index with fields tsince, time (POSIXct),
    #'             wrindex (factor), and wrname (factor)
    #'   - LC: full life cycle with fields tsince, time (POSIXct), wrindex (factor),
    #'         and wrname (factor)
    #'
    #' @author R Translation by Henry Schoeller, Original by Dominik Bueeler
    #' @date April 2025

    # Helper function to extract indices for desired hours (inside main function)
    extract_hour_indices <- function(dtimes, hours) {
        if (!is.list(hours) && !is.vector(hours)) {
            hours <- c(hours)
        }

        # Vectorized extraction of hours from datetime vector
        dt_hours <- as.integer(format(dtimes, "%H"))
        numeric_hours <- as.integer(hours)

        # Vectorized matching of hours
        which(dt_hours %in% numeric_hours)
    }

    ##################
    # Initialization #
    ##################

    # Set basic filenames based on dataset and setup
    if (dataset == "erainterim") {
        if (setup == "z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015") {
            # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr <- "Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local"
            basefname_maxiwr <-
                "Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local"
            basefname_lc <-
                "Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local"
        }
    } else if (dataset == "era5") {
        if (setup == "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019") {
            # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr <- "Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local"
            basefname_maxiwr <-
                "Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local"
            basefname_lc <-
                "Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local"
        }
    }

    # Set input files
    infile_iwr <- base::file.path(basepath, paste0(basefname_iwr, ".txt"))
    infile_maxiwr <- base::file.path(basepath, paste0(basefname_maxiwr, ".txt"))
    infile_lc <- base::file.path(basepath, paste0(basefname_lc, ".txt"))

    # Get order of weather regimes from data (do not change!)
    lines <- base::readLines(infile_iwr, n = 5)
    wrsorder <- base::strsplit(lines[5], "\\s+")[[1]]
    wrsorder <- wrsorder[17:length(wrsorder)] # Adjust for R 1-indexing

    # Replace specific regime names if needed (vectorized)
    replace_indices <- base::match(c("ZOEA", "ZOWE", "BL"), wrsorder)
    if (!all(base::is.na(replace_indices))) {
        wrsorder[replace_indices[!base::is.na(replace_indices)]] <-
            c("ScTr", "EuBL", "ScBL")[!base::is.na(replace_indices)]
    }

    # Get indices of weather regimes from data (do not change!)
    lines <- base::readLines(infile_iwr, n = 3)
    wrsind <- base::strsplit(lines[3], "\\s+")[[1]]
    wrsind <- wrsind[5:length(wrsind)] # Adjust for R 1-indexing

    # Replace specific regime indices (vectorized)
    replace_indices <- base::match(c("ZOEA", "ZOWE", "BL"), wrsind)
    if (!all(base::is.na(replace_indices))) {
        wrsind[replace_indices[!base::is.na(replace_indices)]] <-
            c("ScTr", "EuBL", "ScBL")[!base::is.na(replace_indices)]
    }

    # Define weather regime indices and corresponding names
    regime_names <- c("AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "no")
    wr_indices <- base::c(
        which(wrsind == "AT"),
        which(wrsind == "ZO"),
        which(wrsind == "ScTr"),
        which(wrsind == "AR"),
        which(wrsind == "EuBL"),
        which(wrsind == "ScBL"),
        which(wrsind == "GL"),
        0
    )

    # Define factor levels for consistency
    wrname_levels <- regime_names
    wrindex_levels <- wr_indices

    wrmeta <- base::data.frame(
        wrname = regime_names,
        wrindex = wr_indices,
        stringsAsFactors = FALSE
    )

    #############
    # Load data #
    #############

    data_iwr <- data.table::fread(
        infile_iwr,
        skip = 7,
        select = c(1, 2, 3, 5:11), # Skip column 4
        colClasses = list(integer = 1, character = 2, integer = 3, numeric = 5:11),
        col.names = c("tsince", "time", "cci", wrsorder)
    )

    # For data_maxiwr - selecting specific columns
    data_maxiwr <- data.table::fread(
        infile_maxiwr,
        skip = 7,
        select = c(1, 2, 4), # Skip column 3
        colClasses = list(integer = 1, character = 2, integer = 4),
        col.names = c("tsince", "time", "wrindex")
    )

    # For data_lc - selecting specific columns
    data_lc <- data.table::fread(
        infile_lc,
        skip = 7,
        select = c(1, 2, 5), # Skip columns 3 and 4
        colClasses = list(integer = 1, character = 2, integer = 5),
        col.names = c("tsince", "time", "wrindex")
    )

    ####################
    # Postprocess data #
    ####################

    # Convert datestring to datetime objects
    dtimes <- as.POSIXct(
        strptime(data_iwr$time, format = "%Y%m%d_%H"),
        tz = "UTC"
    )

    # Create output data structures
    data <- list()

    # Optimized regime name lookup using match (vectorized)
    lookup_regime_name <- function(indices) {
        regime_names[base::match(indices, wr_indices)]
    }

    # Create output data structures more efficiently
    data$IWR <- data_iwr
    # Convert time column to POSIXct for IWR data
    data$IWR$time <- dtimes

    # Using data.table for more efficient operations (avoid copies)
    if (requireNamespace("data.table", quietly = TRUE)) {
        # Convert to data.table for in-place operations
        data$MAXIWR <- data.table::as.data.table(data_maxiwr)
        data$MAXIWR[, wrname := lookup_regime_name(wrindex)]

        # Convert to factors with defined levels
        data$MAXIWR[, wrindex := factor(wrindex, levels = wrindex_levels)]
        data$MAXIWR[, wrname := factor(wrname, levels = wrname_levels)]

        # Convert time column to POSIXct
        data$MAXIWR[, time := dtimes]

        data$LC <- data.table::as.data.table(data_lc)
        data$LC[, wrname := lookup_regime_name(wrindex)]

        # Convert to factors with defined levels
        data$LC[, wrindex := factor(wrindex, levels = wrindex_levels)]
        data$LC[, wrname := factor(wrname, levels = wrname_levels)]

        # Convert time column to POSIXct
        data$LC[, time := dtimes]
    } else {
        # Fallback if data.table is not available
        data$MAXIWR <- data_maxiwr
        data$MAXIWR$wrname <- lookup_regime_name(data$MAXIWR$wrindex)

        # Convert to factors with defined levels
        data$MAXIWR$wrindex <- factor(data$MAXIWR$wrindex, levels = wrindex_levels)
        data$MAXIWR$wrname <- factor(data$MAXIWR$wrname, levels = wrname_levels)

        # Convert time column to POSIXct
        data$MAXIWR$time <- dtimes

        data$LC <- data_lc
        data$LC$wrname <- lookup_regime_name(data$LC$wrindex)

        # Convert to factors with defined levels
        data$LC$wrindex <- factor(data$LC$wrindex, levels = wrindex_levels)
        data$LC$wrname <- factor(data$LC$wrname, levels = wrname_levels)

        # Convert time column to POSIXct
        data$LC$time <- dtimes
    }

    # Note: time column is now always in POSIXct format (date objects)
    # The tformat parameter is preserved for backward compatibility

    # Improved date handling for exact matching
    if (tformat == "string") {
        # Convert search dates to POSIXct for comparison
        start_dt <- base::as.POSIXct(
            base::strptime(start, format = "%Y%m%d_%H"),
            tz = "UTC"
        )
        end_dt <- base::as.POSIXct(
            base::strptime(end, format = "%Y%m%d_%H"),
            tz = "UTC"
        )

        # Find exact matches
        ind_start <- base::which(dtimes == start_dt)[1]
        ind_end <- base::which(dtimes == end_dt)[1]

        # Fallback to closest match if exact match isn't found
        if (base::is.na(ind_start)) {
            ind_start <- base::which.min(base::abs(dtimes - start_dt))
        }
        if (base::is.na(ind_end)) {
            ind_end <- base::which.min(base::abs(dtimes - end_dt))
        }
    } else if (tformat == "dtime") {
        # Convert search dates to POSIXct for exact matching
        start_dt <- base::as.POSIXct(
            base::strptime(start, format = "%Y%m%d_%H"),
            tz = "UTC"
        )
        end_dt <- base::as.POSIXct(
            base::strptime(end, format = "%Y%m%d_%H"),
            tz = "UTC"
        )

        # Exact matching instead of using difftime
        ind_start <- base::which(dtimes == start_dt)[1]
        ind_end <- base::which(dtimes == end_dt)[1]

        # Fallback to closest match if exact match isn't found
        if (base::is.na(ind_start)) {
            ind_start <- base::which.min(base::abs(dtimes - start_dt))
        }
        if (base::is.na(ind_end)) {
            ind_end <- base::which.min(base::abs(dtimes - end_dt))
        }
    }

    # Trim data for time period
    range <- ind_start:ind_end

    # More efficient subsetting
    if (requireNamespace("data.table", quietly = TRUE)) {
        # For data.table objects
        data$IWR <- data$IWR[range]
        data$MAXIWR <- data$MAXIWR[range]
        data$LC <- data$LC[range]
    } else {
        # For data.frames
        data$IWR <- data$IWR[range, ]
        data$MAXIWR <- data$MAXIWR[range, ]
        data$LC <- data$LC[range, ]
    }
    dtimes <- dtimes[range]

    # More efficient hour extraction using internal helper function
    hour_indices <- extract_hour_indices(dtimes, hours)

    # Extract data with desired hours if not all hours
    if (dataset == "erainterim" &&
        !identical(as.character(hours), c("00", "06", "12", "18"))) {
        if (requireNamespace("data.table", quietly = TRUE)) {
            data$LC <- data$LC[hour_indices]
            data$MAXIWR <- data$MAXIWR[hour_indices]
            data$IWR <- data$IWR[hour_indices]
        } else {
            data$LC <- data$LC[hour_indices, ]
            data$MAXIWR <- data$MAXIWR[hour_indices, ]
            data$IWR <- data$IWR[hour_indices, ]
        }
        dtimes <- dtimes[hour_indices]
    } else if (dataset == "era5" &&
        !identical(
            as.character(hours),
            c("00", "03", "06", "09", "12", "15", "18", "21")
        )) {
        if (requireNamespace("data.table", quietly = TRUE)) {
            data$LC <- data$LC[hour_indices]
            data$MAXIWR <- data$MAXIWR[hour_indices]
            data$IWR <- data$IWR[hour_indices]
        } else {
            data$LC <- data$LC[hour_indices, ]
            data$MAXIWR <- data$MAXIWR[hour_indices, ]
            data$IWR <- data$IWR[hour_indices, ]
        }
        dtimes <- dtimes[hour_indices]
    }

    # Return result
    return(list(dtimes = dtimes, data = data))
}

#---------------------- Dataset manipulation functions ------------------------

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



aggregate_monthly <- function(df) {
    # Create a year_month column (e.g. "2020-03-01") and ensure it's in date format
    df <- df %>%
        mutate(date = as.Date(format(date, "%Y-%m-01")))

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

#---------------------- WR Composite Functions --------------------------------

permute_test <- function(
    data, category_dates, var_name = "log_variance",
    n_perm = 1000) {
    # Extract data for the category
    category_data <- data[data$time %in% category_dates, ]
    observed_mean <- mean(category_data[[var_name]], na.rm = TRUE)
    # Get baseline mean
    baseline_mean <- mean(data[[var_name]], na.rm = TRUE)
    # Get run blocks for this category
    full_dates <- sort(unique(data$time))
    blocks <- get_run_blocks(full_dates, category_dates)

    # Generate surrogate means
    surrogate_means <- replicate(n_perm, {
        surrogate_dates <- permute_blocks(blocks, full_dates)
        surrogate_data <- data[data$time %in% surrogate_dates, ]
        mean(surrogate_data[[var_name]], na.rm = TRUE)
    })
    # Calculate p-value (two-tailed test)
    # Test if observed anomaly is significantly different from surrogate anomalies
    observed_anomaly <- observed_mean - baseline_mean
    surrogate_anomalies <- surrogate_means - baseline_mean

    p_value <- mean(abs(surrogate_anomalies) >= abs(observed_anomaly))
    return(list(mean = observed_mean, p_val = p_value))
}


# Function to process a single weather regime for a single grid point.
# This function calculates the observed composite value and the permutation
# p-value.
process_wr_composite <- function(
    wr, model_data, residuals_ordered, df,
    wr_lookup, surrogate_dates_list) {
    # Get weather regime name from lookup
    wrname <- wr_lookup$wrname[wr_lookup$wrindex == wr][1]

    # Observed composite: get dates for the wr and calculate the mean residual
    wr_dates_obs <- df$date[df$wrindex == wr]
    date_indices <- which(model_data$date %in% wr_dates_obs)
    obs_resid_subset <- residuals_ordered[date_indices]
    obs_composite <- mean(obs_resid_subset, na.rm = TRUE)

    # Permutation composites
    surrogate_list <- surrogate_dates_list[[as.character(wr)]]
    permuted_composites <- sapply(surrogate_list, function(surr_dates) {
        indices <- which(model_data$date %in% surr_dates)
        res <- residuals_ordered[indices]
        mean(res, na.rm = TRUE)
    })

    # Two-tailed p-value based on the absolute observed composite value.
    # p = probability that a composite mean as least as large (in absolute
    # values) as the observed composite mean occurs by chance.
    p_value <- mean(abs(permuted_composites) >= abs(obs_composite),
        na.rm = TRUE
    )

    tibble(
        wr = wr,
        wrname = wrname,
        mean = obs_composite,
        p_value = p_value
    )
}

# A few helpers for permutation tested composites

get_run_blocks <- function(full_dates, category_dates) {
    # Identify the blocks of consecutive dates
    indicator <- full_dates %in% category_dates
    r <- rle(indicator)
    # Extract only run lengths for the "true" segments (i.e. when the category is present)
    blocks <- r$lengths[r$values]
    return(blocks)
}

# Helper function that takes a vector of block lengths and the full sorted date vector
# Returns a surrogate date set that preserves roughly
# the block length distribution. This simple implementation randomly partitions the gaps.
# Deviations from the original block length distribution occurs only if the
# final block overflows.
permute_blocks <- function(blocks, full_dates) {
    N <- length(full_dates)
    total_ones <- sum(blocks) # total days in this category
    k <- length(blocks)
    # Total gap (days not in category)
    total_gap <- N - total_ones
    # Randomly partition the total gap into (k+1) parts
    # Here we use the "stars and bars" idea.
    if (total_gap > 0) {
        # Draw k random numbers between 0 and total_gap and sort
        dividers <- sort(sample(0:total_gap, k, replace = FALSE))
        gaps <- c(dividers[1], diff(dividers), total_gap - dividers[k])
    } else {
        gaps <- rep(0, k + 1)
    }

    # To add some additional randomness, you can also randomly permute the order of the blocks.
    # (This does not change the overall total days but shuffles run orders.)
    blocks_shuffled <- sample(blocks, k, replace = FALSE)
    # Initialize surrogate indicator vector
    surrogate_indicator <- rep(0, N)

    # Track current position in the vector
    current_pos <- 1 + gaps[1] # Start after the first gap

    # Place each block sequentially, accounting for its length and the following gap
    for (i in seq_along(blocks_shuffled)) {
        block_length <- blocks_shuffled[i]
        # Place the current block
        end_pos <- min(current_pos + block_length - 1, N)
        surrogate_indicator[current_pos:end_pos] <- 1

        # Move position to after this block and its following gap
        current_pos <- end_pos + 1 + ifelse(i < length(gaps), gaps[i + 1], 0)

        # Break if we've reached the end of the vector
        if (current_pos > N) break
    }
    surrogate_dates <- full_dates[which(surrogate_indicator == 1)]
    return(surrogate_dates)
}

# Function to process one grid point (one lm file)
process_grid_point <- function(
    lat, lon, file_path, df, wr_lookup, wr_indices,
    surrogate_dates_list) {
    # Load the linear model
    current_lm <- readRDS(file_path)

    # Reconstruct dates from the model data
    model_data <- current_lm$model
    min_year <- min(model_data$year)
    start_date <- as.Date(paste0(min_year, "-01-01"))
    model_data <- model_data %>%
        arrange(year) %>%
        mutate(
            global_day = row_number(),
            date = start_date + (global_day - 1)
        )
    # Order residuals consistently with the dates
    residuals_ordered <- current_lm$residuals[order(model_data$year)]

    # Process each weather regime for this grid point
    category_results <- lapply(wr_indices, function(wr) {
        process_wr_composite(
            wr, model_data, residuals_ordered, df,
            wr_lookup, surrogate_dates_list
        )
    })

    # Combine the results and add grid point coordinates
    grid_result <- bind_rows(category_results) %>%
        mutate(lon = lon, lat = lat)

    return(grid_result)
}

# Top-level function to calculate composite values.
calculate_composite_values <- function(df, lm_dir, n_perm = 10000, n_cores = 4) {
    # Ensure df$date is a Date object
    df$date <- as.Date(df$date)

    # Precompute the full sorted date vector from df (assumed common for all time series)
    full_dates <- sort(unique(df$date))

    # Create lookup for weather regime indices to names
    wr_lookup <- unique(df[, c("wrindex", "wrname")])
    wr_indices <- sort(unique(df$wrindex))

    # Precompute surrogate date series (permutations) for each weather regime once.
    surrogate_dates_list <- list()
    for (wr in wr_indices) {
        actual_dates <- df$date[df$wrindex == wr]
        blocks <- get_run_blocks(full_dates, actual_dates)

        surrogate_list <- vector("list", n_perm)
        for (p in 1:n_perm) {
            surrogate_list[[p]] <- permute_blocks(blocks, full_dates)
        }
        surrogate_dates_list[[as.character(wr)]] <- surrogate_list
    }

    print("Surrogate dates precomputed.")

    # Create grid point info from LM files
    lm_files <- list.files(lm_dir, pattern = "^lm.*\\.Rds", full.names = FALSE)
    coords <- str_extract_all(lm_files, "-?\\d+\\.?\\d*", simplify = TRUE)
    print(coords)
    print(str(coords))
    lm_info <- tibble(
        file = file.path(lm_dir, lm_files),
        lat = as.numeric(coords[, 1]),
        lon = as.numeric(coords[, 2])
    )

    # Process grid points in parallel using mclapply
    grid_results <- parallel::mclapply(1:nrow(lm_info), function(i) {
        file_path <- lm_info$file[i]
        lat <- lm_info$lat[i]
        lon <- lm_info$lon[i]
        print(paste("Processing grid point:", lat, lon))
        process_grid_point(
            lat, lon, file_path, df, wr_lookup, wr_indices,
            surrogate_dates_list
        )
    }, mc.cores = n_cores)

    # Combine results from all grid points
    results <- bind_rows(grid_results)
    results <- results %>% mutate(p_value_adj = p.adjust(p_value, method = "fdr"))
    return(results)
}

calculate_wr_composites <- function(
    z_data, z_dates, lon, lat, wr_data, var_name = "z",
    calculate_variance = FALSE, n_perm = 1000,
    n_cores = parallel::detectCores() - 1) {
    # prepare dates & regimes
    date_lookup <- data.table(
        date = as.Date(z_dates),
        time_idx = seq_along(z_dates)
    )
    merged <- merge(wr_data[, .(date = as.Date(date), wrindex, wrname)],
        date_lookup,
        by = "date"
    )
    regimes <- unique(merged[, .(wrindex, wrname)])
    # precompute variance & surrogates
    if (calculate_variance) {
        overall_var <- apply(z_data, 1:2, var, na.rm = TRUE)
    }

    # pivot helper using pivot_longer, ensure data.table
    proc <- function(mat, name) {
        dt <- as.data.table(mat)
        dt[, lon_idx := seq_len(.N)]
        dt2 <- pivot_longer(dt,
            cols = -lon_idx,
            names_to = "lat_idx", values_to = name
        )
        dt2 <- as.data.table(dt2)
        dt2[, lat_idx := as.integer(gsub("V", "", lat_idx))]
        dt2[, lon_idx := as.integer(lon_idx)]
        dt2
    }

    out <- list()
    for (i in seq_len(nrow(regimes))) {
        wr <- regimes$wrindex[i]
        wname <- regimes$wrname[i]
        idx <- merged[wrindex == wr, time_idx]
        if (length(idx) == 0) next

        # mean composite
        mean_m <- apply(z_data[, , idx, drop = FALSE], 1:2, mean, na.rm = TRUE)
        dt <- proc(mean_m, var_name)[, `:=`(
            lon = lon[lon_idx], lat = lat[lat_idx],
            wrindex = wr, wrname = wname
        )]
        if (calculate_variance) {
            # observed anomaly
            obs <- apply(z_data[, , idx, drop = FALSE], 1:2,
                var,
                na.rm = TRUE
            ) - overall_var

            blocks <- get_run_blocks(z_dates, merged[wrindex == wr, date])
            surrogates <- replicate(
                n_perm,
                permute_blocks(blocks, z_dates), FALSE
            )

            # Parallelized permutation loop
            perm_results <- mclapply(seq_len(n_perm), function(p) {
                tix <- which(as.Date(z_dates) %in% surrogates[[p]])
                apply(z_data[, , tix, drop = FALSE], 1:2, var, na.rm = TRUE) - overall_var
            }, mc.cores = n_cores)

            # Combine results into array
            perm_arr <- array(unlist(perm_results),
                dim = c(length(lon), length(lat), n_perm)
            )
            pmat <- rowMeans(
                sweep(abs(perm_arr), c(1, 2), abs(obs), FUN = "<"),
                dims = 2
            )

            dt <- merge(dt, proc(obs, paste0(var_name, "_var"))[, .(lon_idx, lat_idx, z_var)],
                by = c("lon_idx", "lat_idx")
            )
            dt <- merge(dt, proc(pmat, "pval")[, .(lon_idx, lat_idx, pval)],
                by = c("lon_idx", "lat_idx")
            )
        }

        dt[, c("lon_idx", "lat_idx") := NULL]
        out[[i]] <- dt
    }

    rbindlist(out)
}
