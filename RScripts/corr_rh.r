#!/usr/bin/env Rscript

library(ncdf4)

args <- commandArgs(trailingOnly = TRUE)
resid_file <- args[1]
cov_file   <- args[2]
out_file   <- args[3]

x <- readRDS(resid_file)
res  <- x$residuals
lon  <- x$lon
lat  <- x$lat
time <- as.POSIXct(x$time, origin = "1970-01-01", tz = "UTC")
year <- as.integer(format(time, "%Y"))

nc <- nc_open(cov_file)
cov <- ncvar_get(nc, "r")
nc_close(nc)

dim_res <- dim(res)
nlon  <- dim_res[1]
nlat  <- dim_res[2]
ntime <- dim_res[3]

# flatten space
ngrid <- nlon * nlat
res2 <- matrix(res, nrow = ngrid, ncol = ntime)
cov2 <- matrix(cov, nrow = ngrid, ncol = ntime)

# latitude per grid point and cosine weights
lat_vec <- rep(lat, each = nlon)
w_space <- cos(lat_vec * pi / 180)

# bands <- list(
#   c(59.5, 59.5),
#   c(59, 60),
#   c(58.5, 60.5),
#   c(58, 61),
#   c(57.5, 61.5),
#   c(57, 62),
#   c(56.5, 62.5),
#   c(56, 63),
#   c(55.5, 63.5),
#   c(55, 64),
#   c(54.5, 64.5)
# )

bands <- list(
  c(30, 90)
)

band_label <- function(b) {
  if (b[1] == b[2]) as.character(b[1]) else paste0(b[1], "-", b[2])
}

weighted_cor_pool <- function(xmat, ymat, w) {
  x <- as.vector(xmat)
  y <- as.vector(ymat)
  w <- rep(w, times = ncol(xmat))

  ok <- is.finite(x) & is.finite(y) & is.finite(w) & (w > 0)
  x <- x[ok]; y <- y[ok]; w <- w[ok]

  if (length(x) < 3) return(NA_real_)

  wx <- sum(w * x) / sum(w)
  wy <- sum(w * y) / sum(w)

  xc <- x - wx
  yc <- y - wy

  num <- sum(w * xc * yc)
  den <- sqrt(sum(w * xc^2) * sum(w * yc^2))

  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den
}

years <- sort(unique(year))

out <- do.call(rbind, lapply(bands, function(b) {
  band_rows <- which(lat_vec >= b[1] & lat_vec <= b[2])
  if (!length(band_rows)) return(NULL)

  do.call(rbind, lapply(years, function(y) {
    tt <- which(year == y)
    r <- weighted_cor_pool(
      res2[band_rows, tt, drop = FALSE],
      cov2[band_rows, tt, drop = FALSE],
      w_space[band_rows]
    )

    data.frame(
      band = band_label(b),
      year = y,
      r = r
    )
  }))
}))

saveRDS(out, out_file)

### Old script:

# #!/usr/bin/env Rscript

# library(ncdf4)
# library(parallel)

# args <- commandArgs(trailingOnly = TRUE)
# resid_file <- args[1]
# cov_file   <- args[2]
# out_file   <- args[3]

# lsm_file <- "/scratch/schoelleh96/wp22a/data/land_sea_mask.nc"

# x <- readRDS(resid_file)
# res <- x$residuals
# lon <- x$lon
# lat <- x$lat

# nc <- nc_open(cov_file)
# cov <- ncvar_get(nc, "r")
# nc_close(nc)

# nc <- nc_open(lsm_file)
# lsm <- ncvar_get(nc, "lsm")   # [lat, lon]
# nc_close(nc)

# dim_res <- dim(res)
# nlon <- dim_res[1]
# nlat <- dim_res[2]
# ntime <- dim_res[3]
# ngrid <- nlon * nlat

# res2 <- matrix(res, nrow = ngrid, ncol = ntime)
# cov2 <- matrix(cov, nrow = ngrid, ncol = ntime)

# row_cor <- function(a, b) {
#   ok <- is.finite(a) & is.finite(b)
#   n  <- rowSums(ok)
#   a[!ok] <- NA
#   b[!ok] <- NA
#   ma <- rowMeans(a, na.rm = TRUE)
#   mb <- rowMeans(b, na.rm = TRUE)
#   ac <- a - ma
#   bc <- b - mb
#   num <- rowSums(ac * bc, na.rm = TRUE)
#   den <- sqrt(rowSums(ac^2, na.rm = TRUE) * rowSums(bc^2, na.rm = TRUE))
#   r <- num / den
#   r[n < 3] <- NA
#   r[!is.finite(r)] <- NA
#   r
# }

# weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
#   ok <- is.finite(x) & is.finite(w) & (w > 0)
#   x <- x[ok]; w <- w[ok]
#   if (!length(x)) return(rep(NA_real_, length(probs)))
#   o <- order(x)
#   x <- x[o]; w <- w[o]
#   cw <- cumsum(w) / sum(w)
#   sapply(probs, function(p) x[which(cw >= p)[1]])
# }

# ncores <- detectCores()
# idx <- split(seq_len(ngrid), cut(seq_len(ngrid), ncores, labels = FALSE))

# r_list <- mclapply(
#   idx,
#   function(ii) row_cor(res2[ii, , drop = FALSE], cov2[ii, , drop = FALSE]),
#   mc.cores = ncores
# )

# r <- unlist(r_list)
# r_mat <- array(r, dim = c(nlon, nlat))

# # land mask to [lon, lat]
# land <- t(lsm) >= 0.5
# w2d <- outer(cos(lat * pi / 180), rep(1, nlon))
# w2d <- t(w2d)  # [lon, lat]

# get_summary <- function(mask, label) {
#   qs <- weighted_quantile(r_mat[mask], w2d[mask], c(0.25, 0.5, 0.75))
#   data.frame(
#     spatial_aggregate = label,
#     q1 = qs[1],
#     median = qs[2],
#     q3 = qs[3]
#   )
# }

# summary_df <- rbind(
#   get_summary(matrix(TRUE, nlon, nlat), "all"),
#   get_summary(land, "land"),
#   get_summary(!land, "sea")
# )

# saveRDS(
#   list(
#     correlation = r_mat,
#     lon = lon,
#     lat = lat,
#     summary = summary_df
#   ),
#   out_file
# )