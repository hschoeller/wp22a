#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(emmeans)
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
  stop("Usage: Rscript second_stage_chunk.R <CHUNK_DIR> ",
       "<CHUNK_INDEX> <OUT_DIR> [WR_RDS] [OVERWRITE]")
}

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2])
OUT_DIR   <- args[3]
WR_RDS    <- if (length(args) >= 4) args[4] else 
  "/home/schoelleh96/wp22a/data/wrnames.rds"
OVERWRITE <- if (length(args) >= 5) {
  tolower(args[5]) %in% c("true","1","yes","y")
} else {
  FALSE
}

change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR, sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) {
  stop("Chunk file not found: ", chunk_file)
}
outfile <- file.path(OUT_DIR, 
                     sprintf("stage2_estimates_chunk_%02d.rds", CHUNK_NO))

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

# =========================
# NEW: simulation helpers
# =========================

simulate_sd_diff <- function(emm_raw_obj, emm_full_obj, B = 500) {
  
  raw_df  <- as.data.frame(emm_raw_obj)
  full_df <- as.data.frame(emm_full_obj)
  
  mu_raw  <- raw_df$emmean
  mu_full <- full_df$emmean
  
  Sigma_raw  <- vcov(emm_raw_obj)
  Sigma_full <- vcov(emm_full_obj)
  
  raw_sim  <- mvtnorm::rmvnorm(B, mu_raw,  Sigma_raw)
  full_sim <- mvtnorm::rmvnorm(B, mu_full, Sigma_full)
  
  sd_raw  <- apply(raw_sim, 1, sd)
  sd_full <- apply(full_sim, 1, sd)
  
  delta_sim <- sd_raw - sd_full
  delta_obs <- sd(mu_raw) - sd(mu_full)
  
  p_val <- mean(abs(delta_sim) >= abs(delta_obs))
  
  list(delta = delta_obs, p = p_val)
}

simulate_loo_diff <- function(emm_loo_obj, emm_full_obj, B = 300) {
  
  loo_df  <- as.data.frame(emm_loo_obj)
  full_df <- as.data.frame(emm_full_obj)
  
  mu_loo  <- loo_df$emmean
  mu_full <- full_df$emmean
  
  Sigma_loo  <- vcov(emm_loo_obj)
  Sigma_full <- vcov(emm_full_obj)
  
  loo_sim  <- mvtnorm::rmvnorm(B, mu_loo,  Sigma_loo)
  full_sim <- mvtnorm::rmvnorm(B, mu_full, Sigma_full)
  
  delta_sim <- loo_sim - full_sim
  delta_obs <- mu_loo - mu_full
  
  p_vals <- sapply(seq_along(delta_obs), function(k) {
    mean(abs(delta_sim[, k]) >= abs(delta_obs[k]))
  })
  
  list(delta = delta_obs, p = p_vals)
}

fit_gls <- function(formula, data) {
  gls(
    formula,
    data = data,
    method = "ML",
    correlation = corAR1(form = ~ day_no | segment),
    na.action = na.omit,
    control = glsControl(msMaxIter = 100, msVerbose = FALSE)
  )
}

pattern_size <- function(emm_df) {
  m <- emm_df$emmean
  m <- m[is.finite(m)]
  if (length(m) < 2L) return(NA_real_)
  centered <- m - mean(m, na.rm = TRUE)
  sqrt(mean(centered^2, na.rm = TRUE))
}

extract_lrt_info <- function(anova_tbl) {
  if (is.null(anova_tbl) || nrow(anova_tbl) < 2L) {
    return(list(df = NA_real_, chisq = NA_real_, p = NA_real_))
  }

  nm <- names(anova_tbl)
  p_col  <- grep("^p", nm, ignore.case = TRUE, value = TRUE)
  lr_col <- grep("ratio", nm, ignore.case = TRUE, value = TRUE)
  df_col  <- grep("^df$", nm, ignore.case = TRUE, value = TRUE)

  p_val <- if (length(p_col) >= 1) {
    as.numeric(anova_tbl[[p_col[1]]][nrow(anova_tbl)])
  } else {
    NA_real_
  }
  lr <- if (length(lr_col) >= 1) {
    as.numeric(anova_tbl[[lr_col[1]]][nrow(anova_tbl)])
  } else {
    NA_real_
  }
  df <- if (length(df_col) >= 1) {
    as.numeric(anova_tbl[[df_col[1]]][nrow(anova_tbl)])
  } else {
    NA_real_
  }

  list(df = df, chisq = lr, p = p_val)
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

  if (nrow(d) < 80 || sd(d$resid) == 0 || 
      length(unique(d$wrname)) < 2) {
    return(NULL)
  }

  d$wrname <- factor(d$wrname)
  if (!"no" %in% levels(d$wrname)) return(NULL)
  d$wrname <- relevel(d$wrname, ref = "no")

  log_mem("fit raw model")
  fit_raw <- fit_gls(resid ~ wrname, d)

  log_mem("fit full model")
  full_formula <- as.formula(
    paste0("resid ~ wrname + ", paste(predictors, collapse = " + "))
  )
  fit_full <- fit_gls(full_formula, d)

emm_raw_obj  <- emmeans(fit_raw,  ~ wrname, mode = "asymptotic", data = d)
  emm_full_obj <- emmeans(fit_full, ~ wrname, mode = "asymptotic", data = d)

  emm_raw  <- as.data.frame(summary(emm_raw_obj,  infer = c(TRUE, TRUE)))
  emm_full <- as.data.frame(summary(emm_full_obj, infer = c(TRUE, TRUE)))

  emm_raw$wrname  <- as.character(emm_raw$wrname)
  emm_full$wrname <- as.character(emm_full$wrname)

  merged_emm <- merge(
    emm_raw, emm_full,
    by = "wrname",
    suffixes = c(".raw", ".full"),
    all = FALSE
  )

  if (nrow(merged_emm) < 1) return(NULL)

  wr_emmeans <- data.frame(
    lon = lon_val,
    lat = lat_val,
    i_local = i_local_val,
    j_local = j_local_val,
    wrname = merged_emm$wrname,
    emmean_raw = merged_emm$emmean.raw,
    se_raw = merged_emm$SE.raw,
    p_raw = merged_emm$p.value.raw,
    lowerCL_raw = merged_emm$asymp.LCL.raw,
    upperCL_raw = merged_emm$asymp.UCL.raw,
    emmean_full = merged_emm$emmean.full,
    se_full = merged_emm$SE.full,
    p_full = merged_emm$p.value.full,
    lowerCL_full = merged_emm$asymp.LCL.full,
    upperCL_full = merged_emm$asymp.UCL.full,
    delta_emmean = merged_emm$emmean.raw - merged_emm$emmean.full
  )

  sd_inf <- simulate_sd_diff(emm_raw_obj, emm_full_obj)
  deltaS_p <- sd_inf$p


  log_mem("fit leave-one-out models")
  loo_blocks <- lapply(predictors, function(left_out) {
    remaining <- setdiff(predictors, left_out)
    loo_formula <- if (length(remaining) == 0) {
      resid ~ wrname
    } else {
      as.formula(
        paste0("resid ~ wrname + ", 
               paste(remaining, collapse = " + "))
      )
    }
    fit_loo <- fit_gls(loo_formula, d)
    emm_loo_obj <- emmeans(fit_loo, ~ wrname, mode = "asymptotic", data = d)
    emm_loo <- as.data.frame(summary(emm_loo_obj, infer = c(TRUE, TRUE)))

    emm_loo$wrname <- as.character(emm_loo$wrname)
    lrt_tbl <- tryCatch(
  anova(fit_loo, fit_full),
  error = function(e) NULL
)
lrt_info <- extract_lrt_info(lrt_tbl)

list(
  obj = emm_loo_obj,
  df = emm_loo,
  left_out = left_out,
  lrt_df = lrt_info$df,
  lrt_chisq = lrt_info$chisq,
  lrt_p = lrt_info$p
)
  })

  # ---- LOO per wr inference ----
  loo_attr <- do.call(rbind, lapply(loo_blocks, function(block) {

    sim <- simulate_loo_diff(block$obj, emm_full_obj)

    data.frame(
      lon = lon_val,
      lat = lat_val,
      i_local = i_local_val,
      j_local = j_local_val,
      wrname = block$df$wrname,
      predictor = block$left_out,
      attr_emmean = sim$delta,
      attr_p = sim$p,
      lrt_p = block$lrt_p
    )
  }))

  # ---- LOO RMS inference ----
  loo_rms <- do.call(rbind, lapply(loo_blocks, function(block) {

    sd_inf <- simulate_sd_diff(block$obj, emm_full_obj)

    data.frame(
      lon = lon_val,
      lat = lat_val,
      i_local = i_local_val,
      j_local = j_local_val,
      predictor = block$left_out,
      rms_attr = sd_inf$delta,
      rms_attr_p = sd_inf$p
    )
  }))

  s_raw  <- pattern_size(emm_raw)
  s_full <- pattern_size(emm_full)
  delta_s <- if (is.finite(s_raw) && is.finite(s_full)) {
    s_raw - s_full
  } else {
    NA_real_
  }

  var_y <- var(d$resid)
  sigma_raw  <- fit_raw$sigma^2
  sigma_full <- fit_full$sigma^2

  delta_r2 <- if (is.finite(var_y) && var_y > 0 &&
                  is.finite(sigma_raw) && is.finite(sigma_full)) {
    (sigma_raw - sigma_full) / var_y
  } else {
    NA_real_
  }

  log_mem("lrt")
  lrt_tbl <- tryCatch(
    anova(fit_raw, fit_full),
    error = function(e) NULL
  )
  lrt_info <- extract_lrt_info(lrt_tbl)

  wr_summary <- data.frame(
    n = nrow(d),
    n_wr = nlevels(d$wrname),
    lrt_df = lrt_info$df,
    lrt_chisq = lrt_info$chisq,
    lrt_pvalue = lrt_info$p,
    S_raw = s_raw,
    S_full = s_full,
    deltaS = delta_s,
    deltaS_p = deltaS_p, 
    delta_r2 = delta_r2,
    lon = lon_val,
    lat = lat_val,
    i_local = i_local_val,
    j_local = j_local_val,
    row.names = NULL
  )

  list(
    summary = wr_summary,
    wr_emmeans = wr_emmeans,
    loo_attr = loo_attr,
    loo_rms = loo_rms
  )
}

existing_out <- if (file.exists(outfile) && !OVERWRITE) {
  message("Existing output found, loading for incremental update: ",
          outfile)
  readRDS(outfile)
} else {
  list(
    summary = data.frame(),
    wr_emmeans = data.frame(),
    loo_attr = data.frame(),
    loo_rms = data.frame()
  )
}

if (is.null(existing_out$summary)) {
  existing_out$summary <- data.frame()
}
if (is.null(existing_out$wr_emmeans)) {
  existing_out$wr_emmeans <- data.frame()
}
if (is.null(existing_out$loo_attr)) {
  existing_out$loo_attr <- data.frame()
}
if (is.null(existing_out$loo_rms)) {
  existing_out$loo_rms <- data.frame()
}

nx <- length(chunk$lon)
ny <- length(chunk$lat)

summary_out <- list()
wr_emm_out <- list()
loo_attr_out <- list()
loo_rms_out <- list()

k1 <- k2 <- k3 <- k4 <- 1L

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

    log_mem(paste0("now fitting gridpoint (i,j) = ", ii, ", ", jj))
    ans <- fit_gridpoint(df)
    if (is.null(ans)) next

    summary_out[[k1]] <- ans$summary
    k1 <- k1 + 1
    wr_emm_out[[k2]] <- ans$wr_emmeans
    k2 <- k2 + 1
    loo_attr_out[[k3]] <- ans$loo_attr
    k3 <- k3 + 1
    loo_rms_out[[k4]] <- ans$loo_rms
    k4 <- k4 + 1
  }
}

out <- list(
  metadata = list(
    chunk_no = CHUNK_NO,
    chunk_file = chunk_file,
    source_chunk_metadata = chunk$metadata,
    wr_rds = WR_RDS,
    predictors = predictors,
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
    distinct(i_local, j_local, predictor, .keep_all = TRUE)
)

saveRDS(out, outfile, compress = "xz")
message("Wrote: ", outfile)