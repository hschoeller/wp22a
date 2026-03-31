#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(dplyr)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript build_residual_cube.R <model_dir> <output_rds> [weight_mode] [ncores]")
}

model_dir   <- args[[1]]
output_rds  <- args[[2]]
weight_mode <- if (length(args) >= 3) args[[3]] else "none"       # none | stabilize
ncores      <- if (length(args) >= 4) as.integer(args[[4]]) else detectCores()

if (!weight_mode %in% c("none", "stabilize")) {
  stop("weight_mode must be one of: none, stabilize")
}

extract_lon_lat <- function(path) {
  nm <- basename(path)
  nm <- sub("\\.Rds$", "", nm)
  nm <- sub("^lm", "", nm)
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  c(lat = as.numeric(parts[1]), lon = as.numeric(parts[2]))
}

make_time_index <- function(dat) {
  min_year   <- min(dat$year)
  start_date <- as.Date(sprintf("%d-01-01", min_year))

  dat %>%
    arrange(year) %>%
    mutate(
      time_index = row_number(),
      date = start_date + (time_index - 1)
    )
}

segment_sd_mult <- function(mod) {
  vs <- mod$modelStruct$varStruct
  seg_levels <- attr(vs, "groupNames")
  sd_mult_partial <- coef(vs, unconstrained = FALSE)

  sd_mult <- setNames(rep(1, length(seg_levels)), seg_levels)
  sd_mult[names(sd_mult_partial)] <- as.numeric(sd_mult_partial)
  sd_mult
}

adjust_residuals <- function(mod, dat, res, mode = "none") {
  if (mode == "none") return(res)

  sd_mult <- segment_sd_mult(mod)
  seg_sd  <- unname(sd_mult[as.character(dat$segment)])

  res_adj <- res / seg_sd

  # preserve total variance
  s0 <- sd(res)
  s1 <- sd(res_adj)
  res_adj * (s0 / s1)
}

read_one_model <- function(path, weight_mode) {
  mod <- readRDS(path)
  dat <- make_time_index(mod$data)
  res <- as.numeric(residuals(mod, type = "response"))
  res <- adjust_residuals(mod, dat, res, mode = weight_mode)

  crd <- extract_lon_lat(path)

  list(
    lon  = unname(crd["lon"]),
    lat  = unname(crd["lat"]),
    time = dat$date,
    res  = res
  )
}

files <- list.files(model_dir, pattern = "\\.Rds$", full.names = TRUE)
if (length(files) == 0) stop("No .Rds files found in ", model_dir)

message("Found ", length(files), " model files")
message("Using ", ncores, " cores")

mods <- mclapply(files, read_one_model, weight_mode = weight_mode, mc.cores = ncores)

lon_vals <- sort(unique(vapply(mods, `[[`, numeric(1), "lon")))
lat_vals <- sort(unique(vapply(mods, `[[`, numeric(1), "lat")))
time_ref <- mods[[1]]$time
ntime    <- length(time_ref)

same_time <- vapply(mods, function(x) identical(x$time, time_ref), logical(1))
if (!all(same_time)) stop("Not all models share the same reconstructed time axis")

residuals_array <- array(
  NA_real_,
  dim = c(length(lon_vals), length(lat_vals), ntime),
  dimnames = list(
    lon  = as.character(lon_vals),
    lat  = as.character(lat_vals),
    time = as.character(time_ref)
  )
)

for (x in mods) {
  i <- match(x$lon, lon_vals)
  j <- match(x$lat, lat_vals)
  residuals_array[i, j, ] <- x$res
}

dimnames(residuals_array) <- NULL
out <- list(
  residuals = residuals_array,   # [lon, lat, time]
  lon = lon_vals,
  lat = lat_vals,
  time = time_ref,
  weight_mode = weight_mode
)

saveRDS(out, output_rds)
message("Saved: ", output_rds)
message("Dims [lon, lat, time] = ", paste(dim(residuals_array), collapse = " x "))