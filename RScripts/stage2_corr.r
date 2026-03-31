#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript stage2Corrs.r <CHUNK_DIR> <CHUNK_INDEX> <OUT_DIR> [OVERWRITE]")
}

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2]) + 1
OUT_DIR   <- if (length(args) >= 3) args[3] else file.path(CHUNK_DIR, "corr")
OVERWRITE <- if (length(args) >= 4) tolower(args[4]) %in% c("true", "1", "yes", "y") else FALSE

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR, sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) stop("Chunk file not found: ", chunk_file)

outfile <- file.path(OUT_DIR, sprintf("stage2_corr_chunk_%02d.rds", CHUNK_NO))

if (file.exists(outfile) && !OVERWRITE) {
  message("Output exists, skipping: ", outfile)
  quit(save = "no", status = 0)
}

chunk <- readRDS(chunk_file)

res <- chunk$residuals   # [lon, lat, time]
lon <- chunk$lon
lat <- chunk$lat
time <- chunk$time

if (length(dim(res)) != 3) stop("chunk$residuals must be 3D [lon, lat, time]")

nlon  <- dim(res)[1]
nlat  <- dim(res)[2]
ntime <- dim(res)[3]
ngrid <- nlon * nlat

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

predictor_names <- setdiff(
  names(chunk),
  c("metadata", "lon", "lat", "time", "residuals")
)

corrs <- list()

res2 <- matrix(res, nrow = ngrid, ncol = ntime)

for (pred in predictor_names) {
  x <- chunk[[pred]]
  if (is.null(x)) next
  if (!identical(dim(x), dim(res))) {
    warning("Skipping ", pred, " because dimensions do not match residuals")
    next
  }

  x2 <- matrix(x, nrow = ngrid, ncol = ntime)
  r  <- row_cor(res2, x2)
  corrs[[pred]] <- array(r, dim = c(nlon, nlat))
}

out <- list(
  metadata = chunk$metadata,
  lon = lon,
  lat = lat,
  time = time,
  correlations = corrs
)

saveRDS(out, outfile, compress = "xz")
message("Wrote: ", outfile)