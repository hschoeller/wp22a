#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(mvtnorm)
  library(ncdf4)
  library(parallel)
})

source("RScripts/config.r")

log_mem <- function(tag = "") {
  x <- readLines("/proc/self/status")
  vals <- grep("^(VmRSS|VmHWM|VmSize):", x, value = TRUE)
  to_gb <- function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    name <- sub(":", "", parts[1])
    kb <- as.numeric(parts[2])
    sprintf("%s: %.2f GB", name, kb / 1024^2)
  }
  cat(tag, "\n", sep = "")
  cat(vapply(vals, to_gb, character(1)), sep = "\n")
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop(
        "Usage: Rscript GramsAtt.r ",
        "<RESPONSE_RDS> <COVARIATE_NC> <WR_RDS> <OUTFILE> ",
        "[COVARIATE_VAR] [N_SIM] [N_CORES] [OVERWRITE] ",
        "[BANDWIDTH]"
    )
}

RESPONSE_RDS <- args[1]
COVARIATE_NC <- args[2]
WR_RDS <- args[3]
OUTFILE <- args[4]
COVARIATE_VAR <- if (length(args) >= 5) {
    args[5]
} else {
    "zg500_prime_lp"
}
N_SIM <- if (length(args) >= 6) as.integer(args[6]) else 5000L
BANDWIDTH_OVERRIDE <- if (length(args) >= 9) {
    as.integer(args[9])
} else {
    NA_integer_
}
OVERWRITE <- if (length(args) >= 8) {
    tolower(args[8]) %in% c("true", "1", "yes", "y")
} else {
    FALSE
}
N_CORES <- if (length(args) >= 7) {
    as.integer(args[7])
} else {
    max(1L, detectCores() - 1L)
}

change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")
dir.create(dirname(OUTFILE), recursive = TRUE,
           showWarnings = FALSE)

newey_west_bandwidth <- function(n) {
  max(1L, as.integer(floor(4 * (n / 100)^(2 / 9))))
}

compute_hac_meat <- function(X_a, X_b, u, segment, bandwidth) {
  W_a <- X_a * u
  W_b <- X_b * u
  S <- crossprod(W_a, W_b)
  if (bandwidth < 1L) return(S)

  seg_indices <- split(seq_along(u), segment)
  pa <- ncol(X_a)
  pb <- ncol(X_b)

  for (h in seq_len(bandwidth)) {
    w_h <- 1 - h / (bandwidth + 1L)
    Gamma_pos <- matrix(0, pa, pb)
    Gamma_neg <- matrix(0, pa, pb)

    for (idx in seg_indices) {
      ns <- length(idx)
      if (ns <= h) next
      Wa_seg <- W_a[idx, , drop = FALSE]
      Wb_seg <- W_b[idx, , drop = FALSE]

      Gamma_pos <- Gamma_pos +
        crossprod(Wa_seg[(h + 1L):ns, , drop = FALSE],
                  Wb_seg[1L:(ns - h), , drop = FALSE])
      Gamma_neg <- Gamma_neg +
        crossprod(Wa_seg[1L:(ns - h), , drop = FALSE],
                  Wb_seg[(h + 1L):ns, , drop = FALSE])
    }

    S <- S + w_h * (Gamma_pos + Gamma_neg)
  }
  S
}

make_linfct <- function(X, wr_levels, predictor_cols, data) {
  K <- length(wr_levels)
  p <- ncol(X)
  cn <- colnames(X)
  L <- matrix(0, K, p)
  rownames(L) <- wr_levels
  colnames(L) <- cn
  if ("(Intercept)" %in% cn) L[, "(Intercept)"] <- 1
  for (k in wr_levels[-1L]) {
    dummy_name <- paste0("wrname", k)
    if (dummy_name %in% cn) L[k, dummy_name] <- 1
  }
  for (pred in predictor_cols) {
    if (pred %in% cn) L[, pred] <- mean(data[[pred]])
  }
  L
}

pattern_size_rms <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(NA_real_)
  sqrt(mean((v - mean(v))^2))
}

two_sided_p <- function(draws) {
  draws <- draws[is.finite(draws)]
  if (length(draws) == 0L) return(NA_real_)
  min(1, 2 * min(mean(draws <= 0), mean(draws >= 0)))
}

wald_p_vec <- function(est, se) {
  out <- rep(NA_real_, length(est))
  ok <- is.finite(est) & is.finite(se) & se > 0
  out[ok] <- 2 * pnorm(-abs(est[ok] / se[ok]))
  out
}

bind_prefix <- function(prefix, df) {
  if (nrow(df) == 0L) {
    return(cbind(prefix[0, , drop = FALSE], df))
  }
  out <- cbind(prefix[rep(1L, nrow(df)), , drop = FALSE], df)
  rownames(out) <- NULL
  out
}

compute_joint_inference <- function(d, predictor, n_sim,
                                     bandwidth = NULL) {
  if (nrow(d) < 80) return(NULL)
  if (sd(d$resid) == 0) return(NULL)
  if (length(unique(d$wrname)) < 2) return(NULL)

  d$wrname <- factor(d$wrname)
  if (!"no" %in% levels(d$wrname)) return(NULL)
  d$wrname <- relevel(d$wrname, ref = "no")

  if (is.null(bandwidth) || is.na(bandwidth)) {
    bandwidth <- newey_west_bandwidth(nrow(d))
  }
  d[[predictor]] <- scale(d[[predictor]])[, 1]

  formulas <- list(
    raw  = resid ~ wrname,
    full = as.formula(paste0("resid ~ wrname + ", predictor))
  )

  designs <- lapply(formulas,
                    function(f) model.matrix(f, data = d))
  y <- d$resid

  fits <- vector("list", length(formulas))
  names(fits) <- names(formulas)
  for (nm in names(formulas)) {
    X <- designs[[nm]]
    XtX <- crossprod(X)
    XtX_inv <- tryCatch(solve(XtX),
                        error = function(e) NULL)
    if (is.null(XtX_inv)) return(NULL)
    fits[[nm]] <- list(
      beta    = drop(XtX_inv %*% crossprod(X, y)),
      XtX_inv = XtX_inv
    )
  }

  u_full <- y - drop(designs[["full"]] %*% fits[["full"]]$beta)

  model_names <- names(fits)
  n_models <- length(model_names)
  p_per <- vapply(fits, function(f) length(f$beta),
                  integer(1L))
  offsets <- c(0L, cumsum(p_per))
  p_total <- sum(p_per)

  joint_cov_beta <- matrix(0, p_total, p_total)
  for (i in seq_len(n_models)) {
    rows <- (offsets[i] + 1L):offsets[i + 1L]
    for (j in seq_len(n_models)) {
      cols <- (offsets[j] + 1L):offsets[j + 1L]
      S_ij <- compute_hac_meat(designs[[i]], designs[[j]],
                                u_full, d$segment, bandwidth)
      joint_cov_beta[rows, cols] <-
        fits[[i]]$XtX_inv %*% S_ij %*% fits[[j]]$XtX_inv
    }
  }

  wr_levels <- levels(d$wrname)
  K <- length(wr_levels)

  linfcts <- lapply(model_names, function(nm) {
    make_linfct(designs[[nm]], wr_levels, predictor, d)
  })
  names(linfcts) <- model_names

  emmeans_list <- mapply(function(L, f) drop(L %*% f$beta),
                          linfcts, fits, SIMPLIFY = FALSE)

  big_L <- matrix(0, K * n_models, p_total)
  for (i in seq_len(n_models)) {
    big_L[((i - 1L) * K + 1L):(i * K),
          (offsets[i] + 1L):offsets[i + 1L]] <- linfcts[[i]]
  }
  emm_joint_cov <- big_L %*% joint_cov_beta %*% t(big_L)

  get_block <- function(i, j) {
    rows <- ((i - 1L) * K + 1L):(i * K)
    cols <- ((j - 1L) * K + 1L):(j * K)
    emm_joint_cov[rows, cols, drop = FALSE]
  }

  RAW <- 1L
  FULL <- 2L

  cov_rr <- get_block(RAW, RAW)
  cov_ff <- get_block(FULL, FULL)
  cov_rf <- get_block(RAW, FULL)

  emm_raw <- emmeans_list[[RAW]]
  emm_full <- emmeans_list[[FULL]]
  se_raw <- sqrt(pmax(diag(cov_rr), 0))
  se_full <- sqrt(pmax(diag(cov_ff), 0))

  delta_emm <- emm_raw - emm_full
  var_delta <- diag(cov_rr) + diag(cov_ff) - 2 * diag(cov_rf)
  se_delta <- sqrt(pmax(var_delta, 0))

  z_crit <- qnorm(0.975)

  wr_emmeans_df <- data.frame(
    wrname           = wr_levels,
    emmean_raw       = emm_raw,
    se_raw           = se_raw,
    p_raw            = wald_p_vec(emm_raw, se_raw),
    lowerCL_raw      = emm_raw - z_crit * se_raw,
    upperCL_raw      = emm_raw + z_crit * se_raw,
    emmean_full      = emm_full,
    se_full          = se_full,
    p_full           = wald_p_vec(emm_full, se_full),
    lowerCL_full     = emm_full - z_crit * se_full,
    upperCL_full     = emm_full + z_crit * se_full,
    delta_emmean     = delta_emm,
    se_delta_emmean  = se_delta,
    p_delta_emmean   = wald_p_vec(delta_emm, se_delta),
    lowerCL_delta    = delta_emm - z_crit * se_delta,
    upperCL_delta    = delta_emm + z_crit * se_delta,
    row.names = NULL
  )

  emm_mean_vec <- unlist(emmeans_list)
  emm_cov_sym <- (emm_joint_cov + t(emm_joint_cov)) / 2
  ridge_scale <- max(abs(diag(emm_cov_sym)))
  ridge <- if (is.finite(ridge_scale) && ridge_scale > 0) {
    1e-10 * ridge_scale
  } else {
    1e-12
  }
  emm_cov_safe <- emm_cov_sym + diag(ridge, nrow(emm_cov_sym))

  emm_sim <- tryCatch(
    mvtnorm::rmvnorm(n_sim, emm_mean_vec, emm_cov_safe),
    error = function(e) NULL
  )
  if (is.null(emm_sim)) return(NULL)

  pattern_size_sim <- matrix(NA_real_, n_sim, n_models)
  for (i in seq_len(n_models)) {
    cols <- ((i - 1L) * K + 1L):(i * K)
    block <- emm_sim[, cols, drop = FALSE]
    rowm <- rowMeans(block)
    pattern_size_sim[, i] <- sqrt(rowMeans((block - rowm)^2))
  }

  s_obs <- vapply(emmeans_list, pattern_size_rms,
                  numeric(1L))

  deltaS_obs <- s_obs[RAW] - s_obs[FULL]
  deltaS_draws <- pattern_size_sim[, RAW] -
                  pattern_size_sim[, FULL]
  deltaS_ci <- quantile(deltaS_draws, c(0.025, 0.975),
                        na.rm = TRUE)
  deltaS_p <- two_sided_p(deltaS_draws)

  pred_col <- which(colnames(designs[["full"]]) == predictor)
  pred_global_idx <- offsets[FULL] + pred_col
  beta_pred <- fits[["full"]]$beta[pred_col]
  se_pred <- sqrt(pmax(
    joint_cov_beta[pred_global_idx, pred_global_idx], 0
  ))
  pred_p <- wald_p_vec(beta_pred, se_pred)

  rss_per <- vapply(model_names, function(nm) {
    r <- y - drop(designs[[nm]] %*% fits[[nm]]$beta)
    sum(r^2)
  }, numeric(1L))
  sigma2_per <- rss_per / nrow(d)
  var_y <- var(y)
  delta_r2 <- if (is.finite(var_y) && var_y > 0) {
    unname((sigma2_per[RAW] - sigma2_per[FULL]) / var_y)
  } else {
    NA_real_
  }

  predictor_effect_df <- data.frame(
    predictor = predictor,
    beta      = unname(beta_pred),
    se        = unname(se_pred),
    z         = unname(beta_pred / se_pred),
    p_value   = unname(pred_p),
    lowerCL   = unname(beta_pred - z_crit * se_pred),
    upperCL   = unname(beta_pred + z_crit * se_pred),
    row.names = NULL
  )

  summary_row <- data.frame(
    n                = nrow(d),
    n_wr             = K,
    hac_bandwidth    = bandwidth,
    sigma2_full      = unname(sigma2_per[FULL]),
    S_raw            = unname(s_obs[RAW]),
    S_full           = unname(s_obs[FULL]),
    deltaS           = unname(deltaS_obs),
    lowerCL_deltaS   = unname(deltaS_ci[1]),
    upperCL_deltaS   = unname(deltaS_ci[2]),
    deltaS_p         = deltaS_p,
    delta_r2         = delta_r2,
    row.names = NULL
  )

  list(
    summary       = summary_row,
    wr_emmeans    = wr_emmeans_df,
    predictor_eff = predictor_effect_df
  )
}

log_mem("start")

message("Reading response RDS: ", RESPONSE_RDS)
response <- readRDS(RESPONSE_RDS)

required_response_fields <- c("time", "lon", "lat", "residuals")
missing_fields <- setdiff(required_response_fields,
                          names(response))
if (length(missing_fields) > 0L) {
  stop("Response RDS missing fields: ",
       paste(missing_fields, collapse = ", "))
}

response_time <- as.POSIXct(response$time, tz = "UTC")
response_dates <- as.Date(response_time)
lon_vals <- response$lon
lat_vals <- response$lat
nx <- length(lon_vals)
ny <- length(lat_vals)
nt <- length(response_time)

response_array_dims <- dim(response$residuals)
if (length(response_array_dims) != 3L ||
    !all(response_array_dims == c(nx, ny, nt))) {
  stop("Expected residuals array [nx, ny, nt] = [",
       nx, ", ", ny, ", ", nt, "], got [",
       paste(response_array_dims, collapse = ", "), "]")
}

log_mem("response loaded")

message("Reading WR RDS: ", WR_RDS)
wr_min <- readRDS(WR_RDS)
wr_min$date <- as.Date(wr_min$date)

message("Inspecting netCDF metadata: ", COVARIATE_NC)
nc_meta <- nc_open(COVARIATE_NC)

if (!COVARIATE_VAR %in% names(nc_meta$var)) {
  stop("Variable '", COVARIATE_VAR,
       "' not found in netCDF. Available: ",
       paste(names(nc_meta$var), collapse = ", "))
}

nc_var <- nc_meta$var[[COVARIATE_VAR]]
nc_dim_names <- vapply(nc_var$dim, function(d) d$name,
                       character(1L))
message("netCDF dims for ", COVARIATE_VAR, ": ",
        paste(nc_dim_names, collapse = ", "))

find_dim <- function(candidates) {
  hit <- which(tolower(nc_dim_names) %in% tolower(candidates))
  if (length(hit) == 0L) {
    stop("Could not find dimension among: ",
         paste(candidates, collapse = ", "))
  }
  hit[1]
}

LON_AXIS  <- find_dim(c("lon", "longitude", "x"))
LAT_AXIS  <- find_dim(c("lat", "latitude", "y"))
TIME_AXIS <- find_dim(c("time", "t"))

nc_lon_raw <- nc_var$dim[[LON_AXIS]]$vals
nc_lat_raw <- nc_var$dim[[LAT_AXIS]]$vals

message("---- coordinate axis diagnostics ----")
message("Response lon: n=", length(lon_vals),
        " first=", lon_vals[1],
        " last=", lon_vals[length(lon_vals)],
        " order=",
        if (lon_vals[1] < lon_vals[length(lon_vals)])
          "increasing" else "decreasing")
message("netCDF   lon: n=", length(nc_lon_raw),
        " first=", nc_lon_raw[1],
        " last=", nc_lon_raw[length(nc_lon_raw)],
        " order=",
        if (nc_lon_raw[1] < nc_lon_raw[length(nc_lon_raw)])
          "increasing" else "decreasing")
message("Response lat: n=", length(lat_vals),
        " first=", lat_vals[1],
        " last=", lat_vals[length(lat_vals)],
        " order=",
        if (lat_vals[1] < lat_vals[length(lat_vals)])
          "increasing" else "decreasing")
message("netCDF   lat: n=", length(nc_lat_raw),
        " first=", nc_lat_raw[1],
        " last=", nc_lat_raw[length(nc_lat_raw)],
        " order=",
        if (nc_lat_raw[1] < nc_lat_raw[length(nc_lat_raw)])
          "increasing" else "decreasing")

if (length(nc_lon_raw) != nx) {
  stop("Longitude length mismatch: response=", nx,
       " netCDF=", length(nc_lon_raw))
}
if (length(nc_lat_raw) != ny) {
  stop("Latitude length mismatch: response=", ny,
       " netCDF=", length(nc_lat_raw))
}

normalize_lon <- function(x) {
  ((x + 180) %% 360) - 180
}

build_axis_map <- function(target, source, tol = 1e-3,
                            label = "axis",
                            wrap_lon = FALSE) {
  if (wrap_lon) {
    target_n <- normalize_lon(target)
    source_n <- normalize_lon(source)
  } else {
    target_n <- target
    source_n <- source
  }
  idx <- vapply(target_n, function(v) {
    j <- which.min(abs(source_n - v))
    if (abs(source_n[j] - v) > tol) NA_integer_ else j
  }, integer(1L))
  n_bad <- sum(is.na(idx))
  if (n_bad > 0L) {
    bad_i <- which(is.na(idx))[1]
    stop(label, ": ", n_bad,
         " response coords have no netCDF match (tol=",
         tol, "). First bad target value: ",
         target[bad_i])
  }
  if (anyDuplicated(idx)) {
    stop(label, ": some netCDF indices matched multiple ",
         "response coords; axes may differ in resolution")
  }
  idx
}

lon_map <- build_axis_map(lon_vals, nc_lon_raw,
                           label = "lon",
                           wrap_lon = TRUE)
lat_map <- build_axis_map(lat_vals, nc_lat_raw,
                           label = "lat",
                           wrap_lon = FALSE)

LON_REORDERED <- !all(lon_map == seq_along(lon_map))
LAT_REORDERED <- !all(lat_map == seq_along(lat_map))

message("lon index map: ",
        if (LON_REORDERED) "reordered" else "identity",
        " (e.g. response[1] -> netCDF[", lon_map[1],
        "], response[", nx, "] -> netCDF[",
        lon_map[nx], "])")
message("lat index map: ",
        if (LAT_REORDERED) "reordered" else "identity",
        " (e.g. response[1] -> netCDF[", lat_map[1],
        "], response[", ny, "] -> netCDF[",
        lat_map[ny], "])")
message("--------------------------------------")

nc_time_raw <- nc_var$dim[[TIME_AXIS]]$vals
nc_time_units <- nc_var$dim[[TIME_AXIS]]$units
NC_N_TIME <- length(nc_time_raw)
NC_N_DIM  <- length(nc_dim_names)
message("netCDF time units: ", nc_time_units)

parse_nc_time <- function(vals, units) {
  m <- regmatches(units,
    regexec(
      "^(\\w+)\\s+since\\s+(\\d{4}-\\d{2}-\\d{2})",
      units
    )
  )[[1]]
  if (length(m) < 3L) {
    stop("Could not parse time units: ", units)
  }
  unit <- tolower(m[2])
  origin <- as.POSIXct(m[3], tz = "UTC")
  mult <- switch(unit,
    seconds = 1,
    minutes = 60,
    hours   = 3600,
    days    = 86400,
    stop("Unsupported time unit: ", unit)
  )
  origin + vals * mult
}

nc_time <- parse_nc_time(nc_time_raw, nc_time_units)
nc_dates <- as.Date(nc_time)
response_to_nc <- match(response_dates, nc_dates)
n_unmatched <- sum(is.na(response_to_nc))
if (n_unmatched > 0L) {
  message("Warning: ", n_unmatched,
          " response timestamps have no netCDF match")
}

nc_close(nc_meta)

base_df <- data.frame(
  time     = response_time,
  date     = response_dates,
  time_idx = seq_along(response_time),
  nc_idx   = response_to_nc
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
  ) |>
  group_by(segment) |>
  mutate(day_no = row_number()) |>
  ungroup()

existing_out <- if (file.exists(OUTFILE) && !OVERWRITE) {
  message("Existing output found: ", OUTFILE)
  readRDS(OUTFILE)
} else {
  list(
    summary       = data.frame(),
    wr_emmeans    = data.frame(),
    predictor_eff = data.frame()
  )
}
for (slot in c("summary", "wr_emmeans", "predictor_eff")) {
  if (is.null(existing_out[[slot]])) {
    existing_out[[slot]] <- data.frame()
  }
}

done_keys <- character(0)
if (nrow(existing_out$summary) > 0L) {
  done_keys <- unique(paste(existing_out$summary$i_local,
                            existing_out$summary$j_local,
                            sep = "::"))
}

all_tasks <- expand.grid(ii = seq_len(nx), jj = seq_len(ny),
                         KEEP.OUT.ATTRS = FALSE)
task_keys <- paste(all_tasks$ii, all_tasks$jj, sep = "::")
if (!OVERWRITE) {
  todo <- !task_keys %in% done_keys
  all_tasks <- all_tasks[todo, , drop = FALSE]
}
n_tasks <- nrow(all_tasks)
message("Tasks to run: ", n_tasks, " on ", N_CORES, " cores")
log_mem("before fork")

worker_env <- new.env()
worker_env$nc_handle <- NULL

open_worker_nc <- function() {
  if (is.null(worker_env$nc_handle)) {
    worker_env$nc_handle <- nc_open(COVARIATE_NC)
  }
  worker_env$nc_handle
}

read_covariate_series <- function(ii, jj) {
  nc <- open_worker_nc()
  start_vec <- integer(NC_N_DIM)
  count_vec <- integer(NC_N_DIM)
  start_vec[LON_AXIS]  <- lon_map[ii]
  start_vec[LAT_AXIS]  <- lat_map[jj]
  start_vec[TIME_AXIS] <- 1L
  count_vec[LON_AXIS]  <- 1L
  count_vec[LAT_AXIS]  <- 1L
  count_vec[TIME_AXIS] <- NC_N_TIME
  as.numeric(ncvar_get(nc, COVARIATE_VAR,
                       start = start_vec, count = count_vec))
}

fit_gridpoint <- function(ii, jj) {
  resid_vec <- response$residuals[ii, jj, ]
  if (all(!is.finite(resid_vec))) return(NULL)

  cov_full <- tryCatch(read_covariate_series(ii, jj),
                       error = function(e) {
                         message("nc read err (", ii, ",",
                                 jj, "): ",
                                 conditionMessage(e))
                         NULL
                       })
  if (is.null(cov_full)) return(NULL)

  df <- base_df
  df$resid <- resid_vec[df$time_idx]
  df[[COVARIATE_VAR]] <- ifelse(
    is.na(df$nc_idx), NA_real_, cov_full[df$nc_idx]
  )

  keep <- is.finite(df$resid) &
    is.finite(df[[COVARIATE_VAR]]) &
    is.finite(df$day_no) &
    !is.na(df$wrname)
  d <- df[keep, , drop = FALSE]

  result <- tryCatch(
    compute_joint_inference(d, COVARIATE_VAR, N_SIM,
                            BANDWIDTH_OVERRIDE),
    error = function(e) {
      message("fit err (", ii, ",", jj, "): ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(result)) return(NULL)

  prefix <- data.frame(
    lon = lon_vals[ii],
    lat = lat_vals[jj],
    i_local = ii,
    j_local = jj
  )

  list(
    summary       = bind_prefix(prefix, result$summary),
    wr_emmeans    = bind_prefix(prefix, result$wr_emmeans),
    predictor_eff = bind_prefix(prefix,
                                result$predictor_eff)
  )
}

worker_task <- function(row_idx) {
  ii <- all_tasks$ii[row_idx]
  jj <- all_tasks$jj[row_idx]
  fit_gridpoint(ii, jj)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(20260518L)

t0 <- Sys.time()
results <- if (N_CORES > 1L) {
  mclapply(seq_len(n_tasks), worker_task,
           mc.cores = N_CORES,
           mc.preschedule = FALSE,
           mc.cleanup = TRUE)
} else {
  lapply(seq_len(n_tasks), worker_task)
}
t1 <- Sys.time()
message("Parallel section done in ",
        round(as.numeric(difftime(t1, t0, units = "mins")),
              2), " min")
log_mem("after fork join")

errored <- vapply(results, inherits, logical(1L),
                  "try-error")
if (any(errored)) {
  message("Workers returning try-error: ", sum(errored))
  results <- results[!errored]
}

results <- Filter(Negate(is.null), results)
message("Successful gridpoint fits: ", length(results))

summary_out  <- lapply(results, `[[`, "summary")
wr_emm_out   <- lapply(results, `[[`, "wr_emmeans")
pred_eff_out <- lapply(results, `[[`, "predictor_eff")

out <- list(
  metadata = list(
    response_rds = RESPONSE_RDS,
    covariate_nc = COVARIATE_NC,
    covariate_var = COVARIATE_VAR,
    wr_rds = WR_RDS,
    n_sim = N_SIM,
    hac_bandwidth_override = BANDWIDTH_OVERRIDE,
    n_cores = N_CORES,
    inference_method = "ols_hac_single_covariate_parallel",
    incremental_update = file.exists(OUTFILE) && !OVERWRITE
  ),
  summary = bind_rows(existing_out$summary,
                      bind_rows(summary_out)) |>
    distinct(i_local, j_local, .keep_all = TRUE),
  wr_emmeans = bind_rows(existing_out$wr_emmeans,
                         bind_rows(wr_emm_out)) |>
    distinct(i_local, j_local, wrname, .keep_all = TRUE),
  predictor_eff = bind_rows(existing_out$predictor_eff,
                            bind_rows(pred_eff_out)) |>
    distinct(i_local, j_local, predictor, .keep_all = TRUE)
)

saveRDS(out, OUTFILE, compress = "xz")
message("Wrote: ", OUTFILE)