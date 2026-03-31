#!/usr/bin/env Rscript

library(ncdf4)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
resid_file <- args[1]
hydro_file   <- args[2]
out_file   <- args[3]

lsm_file <- "/scratch/schoelleh96/wp22a/data/land_sea_mask.nc"

x <- readRDS(resid_file)
res <- x$residuals
lon <- x$lon
lat <- x$lat
time <- as.POSIXct(x$time, origin = "1970-01-01", tz = "UTC")
year <- as.integer(format(time, "%Y"))

nc <- nc_open(hydro_file)
hydro <- ncvar_get(nc, "hydrosum")
nc_close(nc)

nc <- nc_open(lsm_file)
lsm <- ncvar_get(nc, "lsm")   # [lat, lon]
nc_close(nc)

dim_res <- dim(res)
nlon <- dim_res[1]
nlat <- dim_res[2]
ntime <- dim_res[3]
ngrid <- nlon * nlat

res2 <- matrix(res, nrow = ngrid, ncol = ntime)
hydro2 <- matrix(hydro, nrow = ngrid, ncol = ntime)

row_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  n  <- rowSums(ok)
  a[!ok] <- NA
  b[!ok] <- NA
  ma <- rowMeans(a, na.rm = TRUE)
  mb <- rowMeans(b, na.rm = TRUE)
  ac <- a - ma
  bc <- b - mb
  num <- rowSums(ac * bc, na.rm = TRUE)
  den <- sqrt(rowSums(ac^2, na.rm = TRUE) * rowSums(bc^2, na.rm = TRUE))
  r <- num / den
  r[n < 3] <- NA
  r[!is.finite(r)] <- NA
  r
}

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

# land mask [lon, lat]
land <- t(lsm) >= 0.5
w2d <- outer(cos(lat * pi / 180), rep(1, nlon))
w2d <- t(w2d)  # [lon, lat]
w <- as.vector(w2d)
land_vec <- as.vector(land)

years <- sort(unique(year))
ncores <- detectCores()

calc_year <- function(y) {
  tt <- which(year == y)
  a <- res2[, tt, drop = FALSE]
  b <- hydro2[, tt, drop = FALSE]

  idx <- split(seq_len(ngrid), cut(seq_len(ngrid), ncores, labels = FALSE))
  r_list <- mclapply(
    idx,
    function(ii) row_cor(a[ii, , drop = FALSE], b[ii, , drop = FALSE]),
    mc.cores = ncores
  )
  r <- unlist(r_list)

  make_row <- function(mask, label) {
    qs <- weighted_quantile(r[mask], w[mask], c(0.25, 0.5, 0.75))
    data.frame(
      year = y,
      median = qs[2],
      q1 = qs[1],
      q3 = qs[3],
      spatial_aggregate = label
    )
  }

  rbind(
    make_row(rep(TRUE, ngrid), "all"),
    make_row(land_vec, "land"),
    make_row(!land_vec, "sea")
  )
}

summary_df <- do.call(rbind, lapply(years, calc_year))

saveRDS(
  list(
    summary = summary_df
  ),
  out_file
)