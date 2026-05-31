#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

source("RScripts/config.r")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript stage2_hac_ultralight.r ",
       "<RESID_FILE> <WR_RDS> <OUT_FILE> ",
       "[N_CORES] [OVERWRITE] [BANDWIDTH]")
}

RESID_FILE <- args[1]
WR_RDS     <- args[2]
OUT_FILE   <- args[3]
N_CORES <- if (length(args) >= 4) {
  as.integer(args[4])
} else {
  max(1L, parallel::detectCores() - 1L)
}
OVERWRITE <- if (length(args) >= 5) {
  tolower(args[5]) %in% c("true", "1", "yes", "y")
} else {
  FALSE
}
BANDWIDTH_OVERRIDE <- if (length(args) >= 6) {
  as.integer(args[6])
} else {
  NA_integer_
}

if (file.exists(OUT_FILE) && !OVERWRITE) {
  stop("Output exists (set OVERWRITE=TRUE to replace): ",
       OUT_FILE)
}

change_points <- as.Date(paste0(CP, "-01"),
                          format = "%Y-%m-%d")
dir.create(dirname(OUT_FILE), recursive = TRUE,
           showWarnings = FALSE)

message("Reading residuals: ", RESID_FILE)
resid_data <- readRDS(RESID_FILE)
message("Reading WR file: ", WR_RDS)
wr_min <- readRDS(WR_RDS)
wr_min$date <- as.Date(wr_min$date)

newey_west_bandwidth <- function(n) {
  max(1L, as.integer(floor(4 * (n / 100)^(2 / 9))))
}

pattern_size_rms <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(NA_real_)
  sqrt(mean((v - mean(v))^2))
}

base_df <- data.frame(
  time     = as.POSIXct(resid_data$time, tz = "UTC"),
  date     = as.Date(resid_data$time),
  time_idx = seq_along(resid_data$time)
) |>
  left_join(wr_min, by = c("date" = "date")) |>
  filter(!is.na(wrname)) |>
  mutate(
    segment = cut(
      date,
      breaks = c(date[1], change_points,
                 date[length(date)] + 1),
      labels = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )
  )

base_df$wrname <- factor(base_df$wrname)
if (!"no" %in% levels(base_df$wrname)) {
  stop("Reference level 'no' not found in wrname.")
}
base_df$wrname <- relevel(base_df$wrname, ref = "no")

X <- model.matrix(~ wrname, data = base_df)
wr_levels <- levels(base_df$wrname)
K <- length(wr_levels)
p <- ncol(X)
wr_cols <- grep("^wrname", colnames(X))
n_t <- nrow(base_df)
time_idx_used <- base_df$time_idx
seg_indices <- split(seq_len(n_t), base_df$segment)

bandwidth <- if (is.na(BANDWIDTH_OVERRIDE)) {
  newey_west_bandwidth(n_t)
} else {
  BANDWIDTH_OVERRIDE
}

XtX     <- crossprod(X)
XtX_inv <- solve(XtX)

nx <- length(resid_data$lon)
ny <- length(resid_data$lat)

message("HAC bandwidth: ", bandwidth)
message("Time points used: ", n_t)
message("Grid points: ", nx * ny)
message("Cores: ", N_CORES)

process_gridpoint <- function(y_full) {
  y <- y_full[time_idx_used]
  if (anyNA(y) || any(!is.finite(y))) return(NULL)
  if (sd(y) == 0) return(NULL)

  beta <- drop(XtX_inv %*% crossprod(X, y))
  u    <- y - drop(X %*% beta)

  W <- X * u
  S <- crossprod(W)

  for (h in seq_len(bandwidth)) {
    w_h   <- 1 - h / (bandwidth + 1L)
    Gamma <- matrix(0, p, p)
    for (idx in seg_indices) {
      ns <- length(idx)
      if (ns <= h) next
      W_seg <- W[idx, , drop = FALSE]
      Gamma <- Gamma +
        crossprod(W_seg[(h + 1L):ns, , drop = FALSE],
                  W_seg[1L:(ns - h), , drop = FALSE])
    }
    S <- S + w_h * (Gamma + t(Gamma))
  }

  V_beta <- XtX_inv %*% S %*% XtX_inv

  emm <- rep(beta[1], K)
  for (k_idx in seq_along(wr_levels)[-1L]) {
    dummy_name <- paste0("wrname", wr_levels[k_idx])
    emm[k_idx] <- beta[1] + beta[dummy_name]
  }
  S_raw <- pattern_size_rms(emm)

  b_wr  <- beta[wr_cols]
  V_wr  <- V_beta[wr_cols, wr_cols, drop = FALSE]
  V_inv <- tryCatch(solve(V_wr), error = function(e) NULL)
  if (is.null(V_inv)) {
    wald_stat <- NA_real_
    wald_p    <- NA_real_
  } else {
    wald_stat <- drop(crossprod(b_wr, V_inv %*% b_wr))
    wald_p    <- 1 - pchisq(wald_stat, length(wr_cols))
  }

  c(S_raw = S_raw, wr_wald_stat = wald_stat,
    wr_wald_p = wald_p)
}

jobs <- expand.grid(ii = seq_len(nx), jj = seq_len(ny),
                    KEEP.OUT.ATTRS = FALSE)

t0 <- Sys.time()
results <- mclapply(seq_len(nrow(jobs)), function(idx) {
  ii <- jobs$ii[idx]
  jj <- jobs$jj[idx]
  y_full <- resid_data$residuals[ii, jj, ]
  res <- process_gridpoint(y_full)
  if (is.null(res)) return(NULL)
  list(
    lon          = resid_data$lon[ii],
    lat          = resid_data$lat[jj],
    i_local      = ii,
    j_local      = jj,
    S_raw        = unname(res["S_raw"]),
    wr_wald_stat = unname(res["wr_wald_stat"]),
    wr_wald_p    = unname(res["wr_wald_p"])
  )
}, mc.cores = N_CORES)
elapsed <- difftime(Sys.time(), t0, units = "mins")
message("Computation finished in ",
        round(as.numeric(elapsed), 2), " min")

results <- Filter(Negate(is.null), results)
summary_df <- do.call(rbind, lapply(results, as.data.frame))
summary_df$n             <- n_t
summary_df$n_wr          <- K
summary_df$hac_bandwidth <- bandwidth
summary_df$wr_wald_df    <- length(wr_cols)

summary_df <- summary_df[, c("lon", "lat", "i_local",
                              "j_local", "n", "n_wr",
                              "hac_bandwidth", "S_raw",
                              "wr_wald_stat", "wr_wald_df",
                              "wr_wald_p")]
rownames(summary_df) <- NULL

out <- list(
  metadata = list(
    resid_file       = RESID_FILE,
    wr_rds           = WR_RDS,
    n_cores          = N_CORES,
    hac_bandwidth    = bandwidth,
    n_time_points    = n_t,
    n_grid_points    = nx * ny,
    inference_method = "ols_hac_raw_only"
  ),
  summary = summary_df
)

saveRDS(out, OUT_FILE, compress = "xz")
message("Wrote: ", OUT_FILE)