#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(mvtnorm)
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
if (length(args) < 3) {
  stop("Usage: Rscript stage2_hac.r <CHUNK_DIR> ",
       "<CHUNK_INDEX> <OUT_DIR> [WR_RDS] [OVERWRITE] ",
       "[N_SIM] [BANDWIDTH]")
}

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2])
OUT_DIR   <- args[3]
WR_RDS    <- if (length(args) >= 4) args[4] else
  "/home/schoelleh96/wp22a/data/wrnames.rds"
OVERWRITE <- if (length(args) >= 5) {
  tolower(args[5]) %in% c("true", "1", "yes", "y")
} else {
  FALSE
}
N_SIM <- if (length(args) >= 6) as.integer(args[6]) else 5000L
BANDWIDTH_OVERRIDE <- if (length(args) >= 7) {
  as.integer(args[7])
} else {
  NA_integer_
}

change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR,
                       sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) {
  stop("Chunk file not found: ", chunk_file)
}
outfile <- file.path(OUT_DIR,
                     sprintf("stage2_estimates_chunk_%02d.rds",
                             CHUNK_NO))

log_mem("start")
message("Reading chunk: ", chunk_file)
chunk <- readRDS(chunk_file)
message("Reading WR file: ", WR_RDS)
wr_min <- readRDS(WR_RDS)
wr_min$date <- as.Date(wr_min$date)

PREDICTOR_NAMES <- c("rh_500_850", "z_laplacian", "z_grad_mag")
predictors <- PREDICTOR_NAMES[PREDICTOR_NAMES %in% names(chunk)]
if (length(predictors) == 0) stop("No predictors found in chunk.")
message("Predictors used: ", paste(predictors, collapse = ", "))
message("MVN simulation draws: ", N_SIM)
if (!is.na(BANDWIDTH_OVERRIDE)) {
  message("HAC bandwidth (override): ", BANDWIDTH_OVERRIDE)
}

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

compute_joint_inference <- function(d, predictors, n_sim,
                                     bandwidth = NULL) {
  if (nrow(d) < 80) {
    message("rejected: nrow(d) < 80"); return(NULL)
  }
  if (sd(d$resid) == 0) {
    message("rejected: residual sd == 0"); return(NULL)
  }
  if (length(unique(d$wrname)) < 2) {
    message("rejected: <2 wrname levels"); return(NULL)
  }

  d$wrname <- factor(d$wrname)
  if (!"no" %in% levels(d$wrname)) {
    message("rejected: missing 'no' reference"); return(NULL)
  }
  d$wrname <- relevel(d$wrname, ref = "no")

  if (is.null(bandwidth) || is.na(bandwidth)) {
    bandwidth <- newey_west_bandwidth(nrow(d))
  }
  for (pred in predictors) {
    d[[pred]] <- scale(d[[pred]])[, 1]
  }
  full_formula <- as.formula(paste0(
    "resid ~ wrname + ", paste(predictors, collapse = " + ")
  ))

  formulas <- list(
    raw  = resid ~ wrname,
    full = full_formula
  )
  for (pname in predictors) {
    remaining <- setdiff(predictors, pname)
    formulas[[paste0("loo_", pname)]] <- if (
      length(remaining) == 0L
    ) {
      resid ~ wrname
    } else {
      as.formula(paste0(
        "resid ~ wrname + ",
        paste(remaining, collapse = " + ")
      ))
    }
  }

  designs <- lapply(formulas,
                    function(f) model.matrix(f, data = d))
  y <- d$resid

  fits <- vector("list", length(formulas))
  names(fits) <- names(formulas)
  for (nm in names(formulas)) {
    X <- designs[[nm]]
    XtX <- crossprod(X)

XtX_inv <- tryCatch(
  solve(XtX),
  error = function(e) {
    message("Singular XtX in model: ", nm)

    message("Dimensions of X:")
    print(dim(X))

    message("Column names:")
    print(colnames(X))

    message("Column SDs:")
    print(apply(X, 2, sd))

    message("Rank of X:")
    print(qr(X)$rank)

    message("ncol(X): ", ncol(X))

    message("Condition number:")
    print(kappa(XtX))

    message("Correlation matrix:")
    print(cor(X))

    NULL
  }
)
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
    make_linfct(designs[[nm]], wr_levels, predictors, d)
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
  LOO_INDS <- 2L + seq_along(predictors)

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

  loo_attr_blocks <- lapply(seq_along(predictors), function(j) {
    loo_idx <- LOO_INDS[j]
    emm_loo <- emmeans_list[[loo_idx]]
    cov_ll <- get_block(loo_idx, loo_idx)
    cov_lf <- get_block(loo_idx, FULL)

    attr_diff <- emm_loo - emm_full
    var_attr <- diag(cov_ll) + diag(cov_ff) - 2 * diag(cov_lf)
    se_attr <- sqrt(pmax(var_attr, 0))

    data.frame(
      wrname       = wr_levels,
      predictor    = predictors[j],
      attr_emmean  = attr_diff,
      se_attr      = se_attr,
      attr_p       = wald_p_vec(attr_diff, se_attr),
      lowerCL_attr = attr_diff - z_crit * se_attr,
      upperCL_attr = attr_diff + z_crit * se_attr,
      row.names = NULL
    )
  })
  loo_attr_df <- do.call(rbind, loo_attr_blocks)

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
  if (is.null(emm_sim)) {
    message("mvtnorm::rmvnorm failed; HAC cov non-PSD?")
    return(NULL)
  }

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

  loo_rms_blocks <- lapply(seq_along(predictors), function(j) {
    loo_idx <- LOO_INDS[j]
    rms_obs <- s_obs[loo_idx] - s_obs[FULL]
    rms_draws <- pattern_size_sim[, loo_idx] -
                 pattern_size_sim[, FULL]
    ci <- quantile(rms_draws, c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      predictor   = predictors[j],
      rms_attr    = unname(rms_obs),
      lowerCL_rms = unname(ci[1]),
      upperCL_rms = unname(ci[2]),
      rms_attr_p  = two_sided_p(rms_draws),
      row.names = NULL
    )
  })
  loo_rms_df <- do.call(rbind, loo_rms_blocks)

  pred_cols <- which(colnames(designs[["full"]]) %in% predictors)
  if (length(pred_cols) > 0L) {
    pred_global_idx <- offsets[FULL] + pred_cols
    beta_pred <- fits[["full"]]$beta[pred_cols]
    V_pred <- joint_cov_beta[pred_global_idx,
                              pred_global_idx, drop = FALSE]
    V_inv <- tryCatch(solve(V_pred), error = function(e) NULL)
    if (is.null(V_inv)) {
      pred_wald_stat <- NA_real_
      pred_wald_p <- NA_real_
    } else {
      pred_wald_stat <- drop(
        crossprod(beta_pred, V_inv %*% beta_pred)
      )
      pred_wald_p <- 1 - pchisq(pred_wald_stat,
                                 length(pred_cols))
    }
    pred_wald_df <- length(pred_cols)
  } else {
    pred_wald_stat <- NA_real_
    pred_wald_p <- NA_real_
    pred_wald_df <- 0L
  }

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

# After computing pred_wald_*, before building summary_row:
full_pred_cols <- which(colnames(designs[["full"]]) %in% predictors)
full_pred_global <- offsets[FULL] + full_pred_cols
full_V_pred <- joint_cov_beta[full_pred_global,
                               full_pred_global, drop = FALSE]
full_beta <- fits[["full"]]$beta[full_pred_cols]
full_se   <- sqrt(pmax(diag(full_V_pred), 0))

predictor_effects_df <- data.frame(
  predictor = colnames(designs[["full"]])[full_pred_cols],
  beta      = full_beta,
  se        = full_se,
  z         = full_beta / full_se,
  p_value   = wald_p_vec(full_beta, full_se),
  lowerCL   = full_beta - z_crit * full_se,
  upperCL   = full_beta + z_crit * full_se,
  row.names = NULL
)

  summary_row <- data.frame(
    n                = nrow(d),
    n_wr             = K,
    hac_bandwidth    = bandwidth,
    sigma2_full      = unname(sigma2_per[FULL]),
    pred_wald_stat   = pred_wald_stat,
    pred_wald_df     = pred_wald_df,
    pred_wald_p      = pred_wald_p,
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
    summary    = summary_row,
    wr_emmeans = wr_emmeans_df,
    loo_attr   = loo_attr_df,
    loo_rms    = loo_rms_df,
    predictor_effects = predictor_effects_df
  )
}

base_df <- data.frame(
  time     = as.POSIXct(chunk$time, tz = "UTC"),
  date     = as.Date(chunk$time),
  time_idx = seq_along(chunk$time)
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

fit_gridpoint <- function(df) {
  lon_val     <- df$lon[1]
  lat_val     <- df$lat[1]
  i_local_val <- df$i_local[1]
  j_local_val <- df$j_local[1]

  d <- df[, c("resid", "segment", "day_no", "wrname",
              predictors), drop = FALSE]

  keep <- is.finite(d$resid) & is.finite(d$day_no)
  for (pred in predictors) {
    keep <- keep & is.finite(d[[pred]])
  }
  keep <- keep & !is.na(d$wrname)
  d <- d[keep, , drop = FALSE]

  log_mem(paste0("gridpoint (i,j) = ",
                 i_local_val, ", ", j_local_val))

  result <- compute_joint_inference(d, predictors, N_SIM,
                                     BANDWIDTH_OVERRIDE)
  if (is.null(result)) return(NULL)

  prefix <- data.frame(
    lon = lon_val,
    lat = lat_val,
    i_local = i_local_val,
    j_local = j_local_val
  )

  list(
    summary    = bind_prefix(prefix, result$summary),
    wr_emmeans = bind_prefix(prefix, result$wr_emmeans),
    loo_attr   = bind_prefix(prefix, result$loo_attr),
    loo_rms    = bind_prefix(prefix, result$loo_rms),
    pred_eff    = bind_prefix(prefix, result$predictor_effects)
  )
}

existing_out <- if (file.exists(outfile) && !OVERWRITE) {
  message("Existing output found: ", outfile)
  readRDS(outfile)
} else {
  list(
    summary    = data.frame(),
    wr_emmeans = data.frame(),
    loo_attr   = data.frame(),
    loo_rms    = data.frame(),
    pred_eff = data.frame()
    )
}
for (slot in c("summary", "wr_emmeans", "loo_attr", "loo_rms", "pred_eff")) {
  if (is.null(existing_out[[slot]])) {
    existing_out[[slot]] <- data.frame()
  }
}

nx <- length(chunk$lon)
ny <- length(chunk$lat)

summary_out  <- list()
wr_emm_out   <- list()
loo_attr_out <- list()
loo_rms_out  <- list()
pred_eff_out <- list()
k1 <- k2 <- k3 <- k4 <- k5 <- 1L

done_keys <- character(0)
if (nrow(existing_out$summary) > 0) {
  done_keys <- unique(paste(existing_out$summary$i_local,
                            existing_out$summary$j_local,
                            sep = "::"))
}

for (ii in seq_len(nx)) {
  for (jj in seq_len(ny)) {
    key <- paste(ii, jj, sep = "::")
    if (!OVERWRITE && key %in% done_keys) next

    df <- base_df
    df$lon <- chunk$lon[ii]
    df$lat <- chunk$lat[jj]
    df$i_local <- ii
    df$j_local <- jj
    df$resid <- chunk$residuals[ii, jj, df$time_idx]

    for (pred in predictors) {
      df[[pred]] <- chunk[[pred]][ii, jj, df$time_idx]
    }

    ans <- fit_gridpoint(df)
    if (is.null(ans)) next

    summary_out[[k1]]  <- ans$summary;    k1 <- k1 + 1L
    wr_emm_out[[k2]]   <- ans$wr_emmeans; k2 <- k2 + 1L
    loo_attr_out[[k3]] <- ans$loo_attr;   k3 <- k3 + 1L
    loo_rms_out[[k4]]  <- ans$loo_rms;    k4 <- k4 + 1L
    pred_eff_out[[k5]]  <- ans$pred_eff;    k5 <- k5 + 1L
  }
}

out <- list(
  metadata = list(
    chunk_no = CHUNK_NO,
    chunk_file = chunk_file,
    source_chunk_metadata = chunk$metadata,
    wr_rds = WR_RDS,
    predictors = predictors,
    n_sim = N_SIM,
    hac_bandwidth_override = BANDWIDTH_OVERRIDE,
    inference_method = "ols_hac_joint",
    incremental_update = file.exists(outfile) && !OVERWRITE
  ),
  summary = bind_rows(existing_out$summary,
                      bind_rows(summary_out)) |>
    distinct(i_local, j_local, .keep_all = TRUE),
  wr_emmeans = bind_rows(existing_out$wr_emmeans,
                          bind_rows(wr_emm_out)) |>
    distinct(i_local, j_local, wrname, .keep_all = TRUE),
  loo_attr = bind_rows(existing_out$loo_attr,
                       bind_rows(loo_attr_out)) |>
    distinct(i_local, j_local, wrname, predictor,
             .keep_all = TRUE),
  loo_rms = bind_rows(existing_out$loo_rms,
                      bind_rows(loo_rms_out)) |>
    distinct(i_local, j_local, predictor, .keep_all = TRUE),
  pred_eff = bind_rows(existing_out$pred_eff,
                      bind_rows(pred_eff_out)) |>
    distinct(i_local, j_local, predictor, .keep_all = TRUE)
)

saveRDS(out, outfile, compress = "xz")
message("Wrote: ", outfile)