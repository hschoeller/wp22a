#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript yearly_grid_corrs.R <RESID_RDS> <NC_FILE> <VAR_NAME> <OUT_RDS>\n",
    "Example:\n",
    "  Rscript yearly_grid_corrs.R residual_cube_weighted.rds relative_humidity_700-850.nc r yearly_corr_rh.rds"
  )
}

RESID_RDS <- args[1]
NC_FILE   <- args[2]
VAR_NAME  <- args[3]
OUT_RDS   <- args[4]

# -----------------------------
# helpers
# -----------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

to_posix_utc <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC") + 12 * 3600)
  y <- as.POSIXct(x, tz = "UTC")
  if (all(is.na(y))) stop("Could not parse response time vector")
  y
}

find_nc_coord_var <- function(nc, patterns) {
  cand <- names(nc$dim)
  hit <- cand[grepl(paste(patterns, collapse = "|"), cand, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])

  candv <- names(nc$var)
  hitv <- candv[grepl(paste(patterns, collapse = "|"), candv, ignore.case = TRUE)]
  if (length(hitv) > 0) return(hitv[1])

  NULL
}

get_nc_coord <- function(nc, patterns) {
  nm <- find_nc_coord_var(nc, patterns)
  if (is.null(nm)) return(NULL)
  if (nm %in% names(nc$dim)) return(nc$dim[[nm]]$vals)
  ncvar_get(nc, nm)
}

parse_nc_time <- function(nc) {
  time_name <- find_nc_coord_var(nc, c("^time$", "time"))
  if (is.null(time_name)) stop("Could not find time variable in NetCDF.")

  vals <- if (time_name %in% names(nc$dim)) nc$dim[[time_name]]$vals else ncvar_get(nc, time_name)
  units_att <- ncatt_get(nc, time_name, "units")$value
  if (is.null(units_att) || is.na(units_att)) stop("Time variable has no units attribute.")

  parts <- strsplit(trimws(units_att), " since ", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("Unsupported NetCDF time units: ", units_att)

  unit <- tolower(trimws(parts[1]))
  origin_txt <- trimws(parts[2])
  origin <- as.POSIXct(origin_txt, tz = "UTC")
  if (is.na(origin)) origin <- as.POSIXct(paste0(origin_txt, " 00:00:00"), tz = "UTC")
  if (is.na(origin)) stop("Could not parse NetCDF time origin: ", origin_txt)

  mult <- switch(unit,
    "seconds" = 1,
    "second"  = 1,
    "secs"    = 1,
    "sec"     = 1,
    "hours"   = 3600,
    "hour"    = 3600,
    "days"    = 86400,
    "day"     = 86400,
    stop("Unsupported NetCDF time unit: ", unit)
  )

  origin + vals * mult
}

convert_lon_to_source <- function(target_lon, source_lon) {
  src_0360 <- all(source_lon >= 0, na.rm = TRUE)
  trg_0360 <- all(target_lon >= 0, na.rm = TRUE)

  if (src_0360 && !trg_0360) return(target_lon %% 360)
  if (!src_0360 && trg_0360) return(ifelse(target_lon > 180, target_lon - 360, target_lon))
  target_lon
}

match_coord_indices <- function(source, target, tol = 1e-4) {
  out <- integer(length(target))
  out[] <- NA_integer_
  for (i in seq_along(target)) {
    hit <- which(abs(source - target[i]) <= tol)
    if (length(hit) > 0) out[i] <- hit[1]
  }
  out
}

match_time_indices <- function(source_time, target_time) {
  idx <- match(target_time, source_time)
  if (all(!is.na(idx))) return(idx)

  out <- integer(length(target_time))
  for (i in seq_along(target_time)) {
    same_day <- as.Date(source_time) == as.Date(target_time[i])
    if (!any(same_day)) {
      out[i] <- NA_integer_
      next
    }
    cand <- which(same_day)
    out[i] <- cand[which.min(abs(as.numeric(difftime(source_time[cand], target_time[i], units = "secs"))))]
  }
  out
}

row_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  n  <- rowSums(ok)

  a[!ok] <- NA_real_
  b[!ok] <- NA_real_

  ma <- rowMeans(a, na.rm = TRUE)
  mb <- rowMeans(b, na.rm = TRUE)

  ac <- a - ma
  bc <- b - mb

  num <- rowSums(ac * bc, na.rm = TRUE)
  den <- sqrt(rowSums(ac^2, na.rm = TRUE) * rowSums(bc^2, na.rm = TRUE))

  r <- num / den
  r[n < 3] <- NA_real_
  r[!is.finite(r)] <- NA_real_
  r
}

read_nc_cube_aligned <- function(nc_path, var_name, target_lon, target_lat, target_time) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  src_lon  <- get_nc_coord(nc, c("^lon$", "longitude"))
  src_lat  <- get_nc_coord(nc, c("^lat$", "latitude"))
  src_time <- parse_nc_time(nc)

  if (is.null(src_lon) || is.null(src_lat)) {
    stop("Could not detect lon/lat in ", nc_path)
  }

  trg_lon_src <- convert_lon_to_source(target_lon, src_lon)
  lon_idx  <- match_coord_indices(src_lon, trg_lon_src)
  lat_idx  <- match_coord_indices(src_lat, target_lat)
  time_idx <- match_time_indices(src_time, target_time)

  if (any(is.na(lon_idx)))  stop("Some target longitudes could not be matched to source file.")
  if (any(is.na(lat_idx)))  stop("Some target latitudes could not be matched to source file.")
  if (any(is.na(time_idx))) stop("Some target times could not be matched to source file.")

  var <- nc$var[[var_name]]
  if (is.null(var)) stop("Variable '", var_name, "' not found in ", nc_path)

  dim_names <- vapply(var$dim, function(d) d$name, character(1))
  dim_lens  <- vapply(var$dim, function(d) d$len, integer(1))

  lon_pos  <- which(grepl("lon|longitude", dim_names, ignore.case = TRUE))[1]
  lat_pos  <- which(grepl("lat|latitude", dim_names, ignore.case = TRUE))[1]
  time_pos <- which(grepl("time", dim_names, ignore.case = TRUE))[1]

  if (is.na(lon_pos))  lon_pos  <- which(dim_lens == length(src_lon))[1]
  if (is.na(lat_pos))  lat_pos  <- which(dim_lens == length(src_lat))[1]
  if (is.na(time_pos)) time_pos <- which(dim_lens == length(src_time))[1]

  if (is.na(lon_pos) || is.na(lat_pos) || is.na(time_pos)) {
    stop("Could not identify lon/lat/time dimensions for variable ", var_name)
  }

  raw <- ncvar_get(nc, var_name, collapse_degen = FALSE)
  raw3 <- aperm(raw, c(lon_pos, lat_pos, time_pos))

  out <- raw3[lon_idx, lat_idx, time_idx, drop = FALSE]
  out
}

# -----------------------------
# load residuals
# -----------------------------

message("Reading residual cube: ", RESID_RDS)
resp <- readRDS(RESID_RDS)

if (!is.list(resp) || is.null(resp$residuals) || is.null(resp$lon) || is.null(resp$lat) || is.null(resp$time)) {
  stop("Residual RDS must contain: residuals, lon, lat, time")
}

res <- resp$residuals   # [lon, lat, time]
lon <- resp$lon
lat <- resp$lat
time_vals <- to_posix_utc(resp$time)

if (length(dim(res)) != 3) {
  stop("Expected residuals array with dims [lon, lat, time]")
}

nlon  <- dim(res)[1]
nlat  <- dim(res)[2]
ntime <- dim(res)[3]
ngrid <- nlon * nlat

if (length(lon) != nlon || length(lat) != nlat || length(time_vals) != ntime) {
  stop("Residual dimensions do not match lon/lat/time lengths")
}

# -----------------------------
# load aligned nc variable
# -----------------------------

message("Reading and aligning NetCDF variable: ", VAR_NAME)
cov_arr <- read_nc_cube_aligned(NC_FILE, VAR_NAME, lon, lat, time_vals)

if (!identical(dim(cov_arr), dim(res))) {
  stop("Aligned covariate cube dims do not match residual cube dims")
}

# -----------------------------
# yearly correlations
# -----------------------------

yrs <- as.integer(format(time_vals, "%Y", tz = "UTC"))
year_levels <- sort(unique(yrs))
nyears <- length(year_levels)

message("Computing yearly correlations for ", nyears, " years")

res2 <- matrix(res, nrow = ngrid, ncol = ntime)
cov2 <- matrix(cov_arr, nrow = ngrid, ncol = ntime)

r_year <- array(NA_real_, dim = c(nlon, nlat, nyears),
                dimnames = list(NULL, NULL, as.character(year_levels)))
n_year <- array(0L, dim = c(nlon, nlat, nyears),
                dimnames = list(NULL, NULL, as.character(year_levels)))

for (k in seq_along(year_levels)) {
  yy <- year_levels[k]
  tt <- which(yrs == yy)

  a <- res2[, tt, drop = FALSE]
  b <- cov2[, tt, drop = FALSE]

  ok <- is.finite(a) & is.finite(b)
  n_ok <- rowSums(ok)

  r <- row_cor(a, b)

  r_year[, , k] <- array(r, dim = c(nlon, nlat))
  n_year[, , k] <- array(n_ok, dim = c(nlon, nlat))

  message("  year ", yy, ": ", length(tt), " timesteps")
}

# -----------------------------
# output
# -----------------------------

out <- list(
  residual_file = RESID_RDS,
  nc_file = NC_FILE,
  variable = VAR_NAME,
  lon = lon,
  lat = lat,
  years = year_levels,
  correlation = r_year,      # [lon, lat, year]
  n_pairs = n_year,          # [lon, lat, year]
  metadata = list(
    dims = c("lon", "lat", "year"),
    note = "Yearly gridpoint-wise Pearson correlations between residuals and aligned covariate"
  )
)

saveRDS(out, OUT_RDS, compress = "xz")
message("Wrote: ", OUT_RDS)