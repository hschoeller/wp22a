#!/usr/bin/env Rscript

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

args <- parse_args(commandArgs(trailingOnly = TRUE))

need <- function(name) {
  if (is.null(args[[name]])) stop("Missing required argument --", name)
  args[[name]]
}

PARTIAL_DIR <- need("partial-dir")
OUT_FILE    <- need("out-file")

files <- sort(list.files(
  PARTIAL_DIR,
  pattern = "^partial_[0-9]{4}\\.rds$",
  full.names = TRUE
))

if (length(files) == 0L) {
  stop("No partial files found in ", PARTIAL_DIR)
}

message("Found ", length(files), " partial files")

acc <- readRDS(files[1])

sum_hit   <- acc$hit_counts
sum_valid <- acc$valid_counts
sum_event <- acc$event_counts
sum_files <- acc$file_counts

lon <- acc$lon
lat <- acc$lat
lag_names <- acc$lag_names
lag_hours <- acc$lag_hours
regimes <- acc$regimes
variables <- acc$variables
nc_variables <- acc$nc_variables

if (length(files) > 1L) {
  for (ff in files[-1]) {
    message("Adding ", basename(ff))
    x <- readRDS(ff)

    if (!identical(dim(x$hit_counts), dim(sum_hit))) {
      stop("Dimension mismatch in ", ff)
    }
    if (!identical(x$regimes, regimes)) {
      stop("Regime mismatch in ", ff)
    }
    if (!identical(x$variables, variables)) {
      stop("Variable mismatch in ", ff)
    }
    if (!identical(x$lag_names, lag_names)) {
      stop("Lag mismatch in ", ff)
    }

    sum_hit   <- sum_hit + x$hit_counts
    sum_valid <- sum_valid + x$valid_counts
    sum_event <- sum_event + x$event_counts
    sum_files <- sum_files + x$file_counts
  }
}

climatology <- sum_hit / sum_valid
climatology[sum_valid == 0] <- NA_real_

out <- list(
  lon = lon,
  lat = lat,
  lag_names = lag_names,
  lag_hours = lag_hours,
  regimes = regimes,
  variables = variables,
  nc_variables = nc_variables,
  climatology = climatology,   # dims: lon x lat x lag x regime x variable
  hit_counts = sum_hit,
  valid_counts = sum_valid,
  event_counts = sum_event,
  file_counts = sum_files,
  partial_files = basename(files)
)

message("Writing final result: ", OUT_FILE)
saveRDS(out, OUT_FILE, compress = FALSE)
message("Done.")