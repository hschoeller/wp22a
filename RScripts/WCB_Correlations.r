#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Simplified assumptions:
# - Linux (use mclapply only)
# - residuals is regular lon x lat x time
# - time axes differ in length → we align by intersection via match()
#   (only overlapping dates are used automatically)

compute_one <- function(dt, residuals, time_res, lon, lat, wr, var) {
  nx <- length(lon)
  ny <- length(lat)

  dsub <- dt[wrname == wr]
  if (nrow(dsub) == 0L) return(NULL)

  # Align times: only keep dates that exist in residuals
  idx <- match(dsub$date, time_res)
  keep <- !is.na(idx) & !is.na(dsub[[var]])
  if (!any(keep)) return(NULL)

  x <- dsub[[var]][keep]
  idx <- idx[keep]

  n   <- matrix(0L, nx, ny)
  sx  <- matrix(0, nx, ny)
  sx2 <- matrix(0, nx, ny)
  sy  <- matrix(0, nx, ny)
  sy2 <- matrix(0, nx, ny)
  sxy <- matrix(0, nx, ny)

  for (k in seq_along(idx)) {
    y <- residuals[,,idx[k]]
    xk <- x[k]

    valid <- !is.na(y)
    if (!any(valid)) next

    y0 <- y
    y0[!valid] <- 0

    n   <- n   + valid
    sx  <- sx  + xk * valid
    sx2 <- sx2 + (xk^2) * valid
    sy  <- sy  + y0
    sy2 <- sy2 + y0^2
    sxy <- sxy + xk * y0
  }

  num <- n * sxy - sx * sy
  den <- sqrt((n * sx2 - sx^2) * (n * sy2 - sy^2))
  r <- num / den

  r[den == 0 | n < 3] <- NA

  data.table(
    wrname = wr,
    variable = var,
    lon = rep(lon, times = ny),
    lat = rep(lat, each = nx),
    r = as.vector(r)
  )
}

compute_all <- function(dt, res_obj, ncores = parallel::detectCores() - 1) {
  dt <- as.data.table(dt)

  residuals <- res_obj$residuals
  lon <- res_obj$lon
  lat <- res_obj$lat
  time_res <- as.Date(res_obj$time)

  # Optional micro-optimization: drop dt rows before residuals start
  dt <- dt[date >= min(time_res)]

  vars <- setdiff(names(dt), c("date", "wrname"))
  wrs <- unique(dt$wrname)

  tasks <- CJ(wr = wrs, var = vars)

  res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    compute_one(dt, residuals, time_res,
                lon, lat,
                tasks$wr[i], tasks$var[i])
  }, mc.cores = ncores)

  rbindlist(res, use.names = TRUE, fill = TRUE)
}

# CLI
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script.R dt.rds residuals.rds out.rds [ncores]")
}

ncores <- if (length(args) >= 4) as.integer(args[4]) else parallel::detectCores() - 1

dt <- readRDS(args[1])
res_obj <- readRDS(args[2])

out <- compute_all(dt, res_obj, ncores)

saveRDS(out, args[3])
