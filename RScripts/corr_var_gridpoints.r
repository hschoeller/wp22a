#!/usr/bin/env Rscript

library(ncdf4)

args       <- commandArgs(trailingOnly = TRUE)
resid_file <- args[1]
cov_file   <- args[2]
cov_var    <- args[3]
out_file   <- args[4]
start_year <- as.integer(args[5])
end_year   <- as.integer(args[6])

x    <- readRDS(resid_file)
res  <- x$residuals
lon  <- x$lon
lat  <- x$lat
time <- as.POSIXct(x$time, origin = "1970-01-01", tz = "UTC")
year <- as.integer(format(time, "%Y"))

nc        <- nc_open(cov_file)
cov_raw   <- ncvar_get(nc, cov_var)
dim_names <- sapply(nc$var[[cov_var]]$dim, function(d) d$name)
cov_lon   <- as.numeric(ncvar_get(nc, "longitude"))
cov_lat   <- as.numeric(ncvar_get(nc, "latitude"))
nc_close(nc)

lon_pos  <- which(dim_names %in% c("lon", "longitude"))
lat_pos  <- which(dim_names %in% c("lat", "latitude"))
time_pos <- which(dim_names %in% c("time", "t"))

if (length(lon_pos) != 1 || length(lat_pos) != 1 || length(time_pos) != 1) {
  stop(sprintf(
    "Could not unambiguously identify lon/lat/time dims. Found: %s",
    paste(dim_names, collapse = ", ")
  ))
}

cov <- aperm(cov_raw, c(lon_pos, lat_pos, time_pos))

lon_idx <- match(round(lon, 6), round(cov_lon, 6))
lat_idx <- match(round(lat, 6), round(cov_lat, 6))

if (anyNA(lon_idx) || anyNA(lat_idx)) {
  stop("Lon/lat mismatch between residuals and covariate file.")
}

cov <- cov[lon_idx, lat_idx, ]

period_idx <- which(year >= start_year & year <= end_year)
if (length(period_idx) < 3) {
  stop(sprintf("Fewer than 3 time steps in %d-%d.", start_year, end_year))
}

res <- res[, , period_idx]
cov <- cov[, , period_idx]

nlon  <- dim(res)[1]
nlat  <- dim(res)[2]
ntime <- dim(res)[3]
ngrid <- nlon * nlat

res2 <- matrix(res, nrow = ngrid, ncol = ntime)
cov2 <- matrix(cov, nrow = ngrid, ncol = ntime)

gridpoint_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

r_vec <- mapply(
  gridpoint_cor,
  asplit(res2, 1),
  asplit(cov2, 1)
)

correlation_map <- array(r_vec, dim = c(nlon, nlat))

saveRDS(
  list(
    variable    = cov_var,
    correlation = correlation_map,
    lon         = lon,
    lat         = lat
  ),
  out_file,
  compress = "xz"
)