#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

WR_RDS <- get_arg("--wr-rds", Sys.getenv("WR_RDS"))
TASK_ID <- as.integer(get_arg("--task-id",
              Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

IN_DIR  <- get_arg("--in-dir",  Sys.getenv("IN_DIR", "./stage2_chunks_res"))
OUT_DIR <- get_arg("--out-dir", Sys.getenv("OUT_DIR", "./stage3_results"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CHUNK_RDS <- file.path(IN_DIR, sprintf("chunk_%02d.rds", TASK_ID))
OUT_RDS   <- file.path(OUT_DIR, sprintf("chunk_%02d.rds", TASK_ID))

message("Task ", TASK_ID, " | reading chunk")

chunk <- readRDS(CHUNK_RDS)
wr_obj <- readRDS(WR_RDS)

# --- WR mapping ---
wr_map_df <- data.frame(
  date = as.Date(wr_obj$date),
  wr = factor(wr_obj$wrname)
)
wr_map_df <- wr_map_df[!duplicated(wr_map_df$date), ]
wr_levels <- levels(wr_map_df$wr)

wr_map <- setNames(as.character(wr_map_df$wr), as.character(wr_map_df$date))

time_vals <- as.POSIXct(chunk$time, tz = "UTC")
date_vals <- as.Date(time_vals)
wr_per_time <- factor(wr_map[as.character(date_vals)], levels = wr_levels)

# --- GLS ---
fit_gls <- function(formula, data) {
  tryCatch(
    gls(
      formula,
      data = data,
      correlation = corAR1(form = ~ day_no | segment),
      method = "ML"
    ),
    error = function(e) NULL
  )
}

extract_cell <- function(i, j) {

  res <- chunk$residuals[i, j, ]
  inf <- chunk$wcb_in_12utc[i, j, ]
  asc <- chunk$wcb_asc_12utc[i, j, ]
  out <- chunk$wcb_out_12utc[i, j, ]

  keep <- is.finite(res) & is.finite(inf) &
          is.finite(asc) & is.finite(out) &
          !is.na(wr_per_time)

  if (!any(keep)) return(NULL)

  dat <- data.frame(
    residual = res[keep],
    inf = inf[keep],
    asc = asc[keep],
    out = out[keep],
    wr = wr_per_time[keep]
  )

  if (nrow(dat) < 20) return(NULL)

  dat$day_no <- seq_len(nrow(dat))
  dat$segment <- factor(1L)

  dat
}

pattern_size <- function(emm_df) {
  m <- emm_df$emmean
  m <- m[is.finite(m)]
  if (length(m) < 2L) return(NA_real_)
  centered <- m - mean(m, na.rm = TRUE)
  sqrt(mean(centered^2, na.rm = TRUE))
}
pattern_size_contrast <- function(con_df) {
  c <- con_df$estimate
  c <- c[is.finite(c)]
  if (length(c) < 1L) return(NA_real_)
  sqrt(mean(c^2, na.rm = TRUE))
}

extract_lrt_info <- function(anova_tbl) {
  if (is.null(anova_tbl) || nrow(anova_tbl) < 2L) {
    return(list(df = NA_real_, chisq = NA_real_, p = NA_real_))
  }

  nm <- names(anova_tbl)
  p_col  <- grep("^p", nm, ignore.case = TRUE, value = TRUE)
  lr_col <- grep("ratio", nm, ignore.case = TRUE, value = TRUE)
  df_col  <- grep("^df$", nm, ignore.case = TRUE, value = TRUE)

  p_val <- if (length(p_col) >= 1) as.numeric(anova_tbl[[p_col[1]]][nrow(anova_tbl)]) else NA_real_
  lr    <- if (length(lr_col) >= 1) as.numeric(anova_tbl[[lr_col[1]]][nrow(anova_tbl)]) else NA_real_
  df    <- if (length(df_col) >= 1) as.numeric(anova_tbl[[df_col[1]]][nrow(anova_tbl)]) else NA_real_

  list(df = df, chisq = lr, p = p_val)
}

process_cell <- function(i, j) {

  dat <- extract_cell(i, j)
  if (is.null(dat)) return(NULL)

  lon <- chunk$lon[i]
  lat <- chunk$lat[j]

  # =========================
  # 1) WCB model
  # =========================
  m_wcb <- fit_gls(residual ~ inf + asc + out, dat)
  if (is.null(m_wcb)) return(NULL)

  eff_inf <- as.data.frame(summary(
    contrast(emmeans(m_wcb, ~ inf, data = dat, mode = "df.error"), "revpairwise"),
    infer = c(TRUE, TRUE)
  ))
  eff_asc <- as.data.frame(summary(
    contrast(emmeans(m_wcb, ~ asc, data = dat, mode = "df.error"), "revpairwise"),
    infer = c(TRUE, TRUE)
  ))
  eff_out <- as.data.frame(summary(
    contrast(emmeans(m_wcb, ~ out, data = dat, mode = "df.error"), "revpairwise"),
    infer = c(TRUE, TRUE)
  ))

  wcb_df <- data.frame(
    lon = lon,
    lat = lat,
    inf = eff_inf$estimate,
    inf_se = eff_inf$SE,
    inf_df = eff_inf$df,
    inf_t = eff_inf$t.ratio,
    inf_p = eff_inf$p.value,
    asc = eff_asc$estimate,
    asc_se = eff_asc$SE,
    asc_df = eff_asc$df,
    asc_t = eff_asc$t.ratio,
    asc_p = eff_asc$p.value,
    out = eff_out$estimate,
    out_se = eff_out$SE,
    out_df = eff_out$df,
    out_t = eff_out$t.ratio,
    out_p = eff_out$p.value
  )

  # =========================
  # 2) WR reduction (vs "no")
  # =========================
  m_wr  <- fit_gls(residual ~ wr, dat)
  m_all <- fit_gls(residual ~ wr + inf + asc + out, dat)

  if (is.null(m_wr) || is.null(m_all)) return(NULL)

  # EMMs for WR levels in each model
  emm_wr1 <- as.data.frame(summary(
    emmeans(m_wr, ~ wr, data = dat, mode = "df.error"),
    infer = c(TRUE, TRUE)
  ))
  emm_wr2 <- as.data.frame(summary(
    emmeans(m_all, ~ wr, data = dat, mode = "df.error"),
    infer = c(TRUE, TRUE)
  ))

  # Align by WR level
  emm_wr1$wr <- as.character(emm_wr1$wr)
  emm_wr2$wr <- as.character(emm_wr2$wr)

  merged_emm <- merge(
    emm_wr1,
    emm_wr2,
    by = "wr",
    suffixes = c(".raw", ".full"),
    all = FALSE,
    sort = FALSE
  )

  if (nrow(merged_emm) < 1L) return(NULL)
  wr_emm_df <- data.frame(
    lon = rep(lon, nrow(merged_emm)),
    lat = rep(lat, nrow(merged_emm)),
    wr = merged_emm$wr,

    emmean_raw = merged_emm$emmean.raw,
    se_raw = merged_emm$SE.raw,
    df_raw = merged_emm$df.raw,
  lowerCL_raw = merged_emm$lower.CL.raw,
  upperCL_raw = merged_emm$upper.CL.raw,
    t_raw = merged_emm$t.ratio.raw,
    p_raw = merged_emm$p.value.raw,

    emmean_full = merged_emm$emmean.full,
    se_full = merged_emm$SE.full,
    df_full = merged_emm$df.full,
  lowerCL_full = merged_emm$lower.CL.full,
  upperCL_full = merged_emm$upper.CL.full,
    t_full = merged_emm$t.ratio.full,
    p_full = merged_emm$p.value.full,

    delta_emmean = merged_emm$emmean.raw - merged_emm$emmean.full,
    row.names = NULL
  )

  # Summary of overall WR-pattern variability change
  s_raw  <- pattern_size(emm_wr1)
  s_full <- pattern_size(emm_wr2)

  delta_s <- if (is.finite(s_raw) && is.finite(s_full)) s_raw - s_full else NA_real_

  var_y <- var(dat$residual)
  sigma_raw  <- m_wr$sigma^2
  sigma_full <- m_all$sigma^2

  delta_r2 <- if (is.finite(var_y) && var_y > 0 &&
                  is.finite(sigma_raw) && is.finite(sigma_full)) {
    (sigma_raw - sigma_full) / var_y
  } else {
    NA_real_
  }

  # Contrasts vs "no" in each model
  con_wr1 <- as.data.frame(summary(
    contrast(
      emmeans(m_wr, ~ wr, data = dat, mode = "df.error"),
      method = "trt.vs.ctrl",
      ref = "no"
    ),
    infer = c(TRUE, TRUE)
  ))

  con_wr2 <- as.data.frame(summary(
    contrast(
      emmeans(m_all, ~ wr, data = dat, mode = "df.error"),
      method = "trt.vs.ctrl",
      ref = "no"
    ),
    infer = c(TRUE, TRUE)
  ))

  con_wr1$contrast <- as.character(con_wr1$contrast)
  con_wr2$contrast <- as.character(con_wr2$contrast)

  merged_con <- merge(
    con_wr1,
    con_wr2,
    by = "contrast",
    suffixes = c(".raw", ".full"),
    all = FALSE,
    sort = FALSE
  )

  if (nrow(merged_con) < 1L) return(NULL)

  wr_contrast_df <- data.frame(
    lon = rep(lon, nrow(merged_con)),
    lat = rep(lat, nrow(merged_con)),
    contrast = merged_con$contrast,

    estimate_raw = merged_con$estimate.raw,
    se_raw = merged_con$SE.raw,
    df_raw = merged_con$df.raw,
    lowerCL_raw = merged_con$lower.CL.raw,
    upperCL_raw = merged_con$upper.CL.raw,
    t_raw = merged_con$t.ratio.raw,
    p_raw = merged_con$p.value.raw,

    estimate_full = merged_con$estimate.full,
    se_full = merged_con$SE.full,
    df_full = merged_con$df.full,
    lowerCL_full = merged_con$lower.CL.full,
    upperCL_full = merged_con$upper.CL.full,
    t_full = merged_con$t.ratio.full,
    p_full = merged_con$p.value.full,

    delta_contrast = merged_con$estimate.raw - merged_con$estimate.full,
    row.names = NULL
  )
s_con_raw  <- pattern_size_contrast(con_wr1)
s_con_full <- pattern_size_contrast(con_wr2)

delta_s_con <- if (is.finite(s_con_raw) && is.finite(s_con_full)) {
  s_con_raw - s_con_full
} else {
  NA_real_
}
  lrt_tbl <- tryCatch(
    anova(fit_raw, fit_full),
    error = function(e) NULL
  )
  lrt_info <- extract_lrt_info(lrt_tbl)
wr_summary_df <- data.frame(
  lon = lon,
  lat = lat,
  lrt_df = lrt_info$df,
  lrt_chisq = lrt_info$chisq,
  lrt_pvalue = lrt_info$p,

  S_raw = s_raw,
  S_full = s_full,
  deltaS = delta_s,

  S_con_raw = s_con_raw,
  S_con_full = s_con_full,
  deltaS_con = delta_s_con,

  delta_r2 = delta_r2,
  row.names = NULL
)

  rm(dat, m_wcb, m_wr, m_all)
  invisible(gc(FALSE))

  list(
    wcb = wcb_df,
    wr_emm = wr_emm_df,
    wr_contrast = wr_contrast_df,
    sum = wr_summary_df
  )
}

# --- loop ---
lon_n <- length(chunk$lon)
lat_n <- length(chunk$lat)

res_wcb <- list()
res_emm <- list()
res_con <- list()
res_sum <- list()

k <- 0
for (i in seq_len(lon_n)) {
  for (j in seq_len(lat_n)) {
    message("Processing cell ", i, ",", j)

    out <- process_cell(i, j)
    if (is.null(out)) next

    k <- k + 1
    res_wcb[[k]] <- out$wcb
    res_emm[[k]] <- out$wr_emm
    res_con[[k]] <- out$wr_contrast
    res_sum[[k]] <- out$sum

    if (k %% 5 == 0) gc(FALSE)
  }
}

safe_rbind <- function(x) {
  if (length(x) == 0L) return(NULL)
  do.call(rbind, x)
}

res_wcb <- safe_rbind(res_wcb)
res_emm <- safe_rbind(res_emm)
res_con <- safe_rbind(res_con)
res_sum <- safe_rbind(res_sum)

saveRDS(
  list(
    wcb_effects = res_wcb,
    wr_emmeans = res_emm,
    wr_contrasts = res_con,
    wr_summary = res_sum
  ),
  OUT_RDS,
  compress = "xz"
)