#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ncdf4)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
    out <- list()
    i <- 1L
    while (i <= length(x)) {
        key <- x[i]
        if (!startsWith(key, "--")) stop("Expected --key, got: ", key)
        key <- sub("^--", "", key)
        if (i == length(x) || startsWith(x[i + 1L], "--")) {
            out[[key]] <- TRUE
            i <- i + 1L
        } else {
            out[[key]] <- x[i + 1L]
            i <- i + 2L
        }
    }
    out
}

opt <- parse_args(args)

need <- function(name) {
    if (is.null(opt[[name]])) stop("Missing required argument --", name)
    opt[[name]]
}

ROOT <- need("root")
OUT_CSV <- need("out")
START_YEAR <- if (!is.null(opt[["start-year"]])) as.integer(opt[["start-year"]]) else NA_integer_
END_YEAR <- if (!is.null(opt[["end-year"]])) as.integer(opt[["end-year"]]) else NA_integer_

# -----------------------------
# helpers
# -----------------------------
list_hit_files <- function(root, start_year = NA_integer_, end_year = NA_integer_) {
    years <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    years <- years[grepl("/[0-9]{4}$", years)]

    if (!is.na(start_year)) {
        yy <- as.integer(basename(years))
        keep <- yy >= start_year
        years <- years[keep]
    }
    if (!is.na(end_year)) {
        yy <- as.integer(basename(years))
        keep <- yy <= end_year
        years <- years[keep]
    }

    files <- character()

    for (yd in sort(years)) {
        months <- list.dirs(yd, recursive = FALSE, full.names = TRUE)
        months <- months[grepl("/[0-9]{2}$", months)]
        if (length(months) == 0L) next

        for (md in sort(months)) {
            ff <- list.files(md, full.names = TRUE, recursive = FALSE)
            if (length(ff) == 0L) next

            # only _00 or _12, optionally with .nc suffix
            keep <- grepl("hit_[0-9]{8}_(00|12)(\\.nc)?$", basename(ff))
            ff <- ff[keep]
            if (length(ff) > 0L) files <- c(files, ff)
        }
    }

    sort(unique(files))
}

pick_dim_name <- function(nc, candidates) {
    dn <- names(nc$dim)
    hit <- candidates[candidates %in% dn]
    if (length(hit) > 0L) {
        return(hit[1])
    }

    vn <- names(nc$var)
    hit <- candidates[candidates %in% vn]
    if (length(hit) > 0L) {
        return(hit[1])
    }

    NA_character_
}

read_coords <- function(path) {
    nc <- nc_open(path)
    on.exit(nc_close(nc), add = TRUE)

    lon_name <- pick_dim_name(nc, c("lon", "longitude", "x"))
    lat_name <- pick_dim_name(nc, c("lat", "latitude", "y"))

    if (is.na(lon_name) || is.na(lat_name)) {
        stop("Could not find lon/lat names in ", path)
    }

    lon <- ncvar_get(nc, lon_name)
    lat <- ncvar_get(nc, lat_name)

    list(
        lon = as.numeric(lon),
        lat = as.numeric(lat),
        lon_name = lon_name,
        lat_name = lat_name
    )
}

same_vec <- function(a, b, tol = 0) {
    if (length(a) != length(b)) {
        return(FALSE)
    }
    if (tol <= 0) {
        return(identical(a, b))
    }
    all(abs(a - b) <= tol)
}

is_constant_or_all_na <- function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0L) {
        return(TRUE)
    }
    length(unique(x2)) == 1L
}

summarize_problem <- function(lon, lat, ref_lon, ref_lat) {
    issues <- character()

    if (length(lon) != length(ref_lon)) issues <- c(issues, "lon_length")
    if (length(lat) != length(ref_lat)) issues <- c(issues, "lat_length")

    if (is_constant_or_all_na(lon)) issues <- c(issues, "lon_constant")
    if (is_constant_or_all_na(lat)) issues <- c(issues, "lat_constant")

    if (length(lon) == length(ref_lon) && !identical(lon, ref_lon)) {
        if (all(sort(lon) == sort(ref_lon))) {
            issues <- c(issues, "lon_permuted")
        } else {
            issues <- c(issues, "lon_values")
        }
    }

    if (length(lat) == length(ref_lat) && !identical(lat, ref_lat)) {
        if (all(sort(lat) == sort(ref_lat))) {
            issues <- c(issues, "lat_permuted")
        } else {
            issues <- c(issues, "lat_values")
        }
    }

    if (length(issues) == 0L) issues <- "unknown_mismatch"
    paste(issues, collapse = ";")
}

# -----------------------------
# main
# -----------------------------
files <- list_hit_files(ROOT, START_YEAR, END_YEAR)

if (length(files) == 0L) {
    stop("No files found under ", ROOT)
}

cat("Found", length(files), "candidate files\n")

ref_file <- files[1]
cat("Reference file:", ref_file, "\n")
ref <- read_coords(ref_file)

results <- vector("list", length(files))
bad_n <- 0L
err_n <- 0L

for (i in seq_along(files)) {
    f <- files[i]

    rec <- tryCatch(
        {
            cc <- read_coords(f)

            lon_ok <- same_vec(cc$lon, ref$lon)
            lat_ok <- same_vec(cc$lat, ref$lat)
            ok <- lon_ok && lat_ok

            problem <- if (ok) "" else summarize_problem(cc$lon, cc$lat, ref$lon, ref$lat)

            data.frame(
                file = f,
                ok = ok,
                problem = problem,
                lon_len = length(cc$lon),
                lat_len = length(cc$lat),
                lon_min = suppressWarnings(min(cc$lon, na.rm = TRUE)),
                lon_max = suppressWarnings(max(cc$lon, na.rm = TRUE)),
                lat_min = suppressWarnings(min(cc$lat, na.rm = TRUE)),
                lat_max = suppressWarnings(max(cc$lat, na.rm = TRUE)),
                stringsAsFactors = FALSE
            )
        },
        error = function(e) {
            err_n <<- err_n + 1L
            data.frame(
                file = f,
                ok = FALSE,
                problem = paste0("read_error: ", conditionMessage(e)),
                lon_len = NA_integer_,
                lat_len = NA_integer_,
                lon_min = NA_real_,
                lon_max = NA_real_,
                lat_min = NA_real_,
                lat_max = NA_real_,
                stringsAsFactors = FALSE
            )
        }
    )

    if (!rec$ok) bad_n <- bad_n + 1L
    results[[i]] <- rec

    if (i %% 500L == 0L || i == length(files)) {
        cat(sprintf(
            "%d / %d checked | bad=%d | read_errors=%d\n",
            i, length(files), bad_n, err_n
        ))
    }
}

tab <- do.call(rbind, results)
write.csv(tab, OUT_CSV, row.names = FALSE)

cat("\nDone.\n")
cat("Reference file: ", ref_file, "\n", sep = "")
cat("Scanned files:  ", nrow(tab), "\n", sep = "")
cat("Bad files:      ", sum(!tab$ok), "\n", sep = "")
cat("Output CSV:     ", OUT_CSV, "\n", sep = "")

if (sum(!tab$ok) > 0L) {
    cat("\nProblem summary:\n")
    print(sort(table(tab$problem[!tab$ok]), decreasing = TRUE))
}
