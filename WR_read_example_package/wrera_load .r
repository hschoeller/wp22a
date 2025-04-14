#' Load Weather Regime Indices and Life Cycles
#'
#' Loads weather regime indices and life cycles from ERA-Interim or ERA5
#' for a specific period. For more information, contact Christian Grams
#' (christian.grams@gmx.de).
#'
#' @param start Character, start date in format "YYYYMMDD_HH"
#' @param end Character, end date in format "YYYYMMDD_HH"
#' @param hours Character vector of hours to include (e.g., c("00", "06"))
#' @param tformat Either "string" or "dtime" for time formatting
#' @param setup String defining the weather regime data source
#' @param dataset Either "erainterim" or "era5"
#' @param basepath Directory path to the raw data
#'
#' @return List with:
#'   - dtimes: POSIXct vector of filtered time points
#'   - data: List with components IWR, MAXIWR, and LC
#' @author R translation by Henry Schoeller, Original code by Dominik Büeler,
#'  Christian Grams
#' @date April 2025
wrera <- function(start, end, hours, tformat, setup, dataset, basepath) {
    install.packages(c("data.table"), dependencies = TRUE)
    # --- Helper function for extracting hour indices ---
    extract_hour_indices <- function(dtimes, hours) {
        hours <- as.character(hours)
        dt_hours <- as.integer(format(dtimes, "%H"))
        which(dt_hours %in% as.integer(hours))
    }

    # --- File name logic based on dataset/setup ---
    if (dataset == "erainterim" && setup == "z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015") {
        basefname <- "Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all"
    } else if (dataset == "era5" && setup == "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019") {
        basefname <- "Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all"
    }

    infile_iwr <- file.path(basepath, paste0(basefname, "_proj_local.txt"))
    infile_maxiwr <- file.path(basepath, paste0(basefname, "_LCO_local.txt"))
    infile_lc <- file.path(basepath, paste0(basefname, "_LCO_local.txt"))

    # --- Extract weather regime order and index ---
    lines_iwr <- readLines(infile_iwr, n = 5)
    wrsorder <- strsplit(lines_iwr[5], "\\s+")[[1]][17:length(strsplit(lines_iwr[5], "\\s+")[[1]])]

    wrsorder[match(c("ZOEA", "ZOWE", "BL"), wrsorder, nomatch = 0)] <- c("ScTr", "EuBL", "ScBL")

    lines_iwr_short <- readLines(infile_iwr, n = 3)
    wrsind <- strsplit(lines_iwr_short[3], "\\s+")[[1]][5:length(strsplit(lines_iwr_short[3], "\\s+")[[1]])]
    wrsind[match(c("ZOEA", "ZOWE", "BL"), wrsind, nomatch = 0)] <- c("ScTr", "EuBL", "ScBL")

    regime_names <- c("AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "no")
    wr_indices <- c(
        which(wrsind == "AT"),
        which(wrsind == "ZO"),
        which(wrsind == "ScTr"),
        which(wrsind == "AR"),
        which(wrsind == "EuBL"),
        which(wrsind == "ScBL"),
        which(wrsind == "GL"),
        0
    )

    wrmeta <- data.frame(
        wrname = regime_names,
        wrindex = wr_indices,
        stringsAsFactors = FALSE
    )

    # --- Load data files using data.table::fread ---
    data_iwr <- data.table::fread(
        infile_iwr,
        skip = 7, select = c(1, 2, 3, 5:11),
        colClasses = list(integer = 1, character = 2, integer = 3, numeric = 5:11),
        col.names = c("tsince", "time", "cci", wrsorder)
    )

    data_maxiwr <- data.table::fread(
        infile_maxiwr,
        skip = 7, select = c(1, 2, 4),
        colClasses = list(integer = 1, character = 2, integer = 4),
        col.names = c("tsince", "time", "wrindex")
    )

    data_lc <- data.table::fread(
        infile_lc,
        skip = 7, select = c(1, 2, 5),
        colClasses = list(integer = 1, character = 2, integer = 5),
        col.names = c("tsince", "time", "wrindex")
    )

    # --- Create unified time vector ---
    dtimes <- as.POSIXct(strptime(data_iwr$time, format = "%Y%m%d_%H"), tz = "UTC")

    # --- Regime name lookup function ---
    lookup_regime_name <- function(indices) {
        regime_names[match(indices, wr_indices)]
    }

    # --- Attach wrname to data sets ---
    data <- list(
        IWR = data_iwr,
        MAXIWR = data_maxiwr[, wrname := lookup_regime_name(wrindex)],
        LC = data_lc[, wrname := lookup_regime_name(wrindex)]
    )

    # --- Time format handling ---
    if (tformat == "dtime") {
        data$IWR$time_obj <- dtimes
        data$MAXIWR$time_obj <- dtimes
        data$LC$time_obj <- dtimes
    }

    # --- Determine time range indices ---
    if (tformat == "string") {
        ind_start <- which(data$IWR$time == start)[1]
        ind_end <- which(data$IWR$time == end)[1]
    } else {
        start_dt <- as.POSIXct(strptime(start, "%Y%m%d_%H"), tz = "UTC")
        end_dt <- as.POSIXct(strptime(end, "%Y%m%d_%H"), tz = "UTC")
        ind_start <- which(dtimes == start_dt)[1]
        ind_end <- which(dtimes == end_dt)[1]

        if (is.na(ind_start)) ind_start <- which.min(abs(dtimes - start_dt))
        if (is.na(ind_end)) ind_end <- which.min(abs(dtimes - end_dt))
    }

    range <- ind_start:ind_end
    dtimes <- dtimes[range]

    # --- Subset data to selected time range ---
    for (name in names(data)) {
        data[[name]] <- data[[name]][range]
    }

    # --- Filter by selected hours ---
    hour_indices <- extract_hour_indices(dtimes, hours)

    if (!identical(
        as.character(hours),
        if (dataset == "erainterim") c("00", "06", "12", "18") else c("00", "03", "06", "09", "12", "15", "18", "21")
    )) {
        dtimes <- dtimes[hour_indices]
        for (name in names(data)) {
            data[[name]] <- data[[name]][hour_indices]
        }
    }

    return(list(dtimes = dtimes, data = data))
}
