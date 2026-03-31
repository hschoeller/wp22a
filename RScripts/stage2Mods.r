#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(emmeans)
  library(dplyr)
})

source("RScripts/config.r")

log_mem <- function(tag = "") {
  x <- readLines("/proc/self/status")
  vals <- grep("^(VmRSS|VmHWM|VmSize):", x, value = TRUE)

  to_gb <- function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    name <- sub(":", "", parts[1])
    kb <- as.numeric(parts[2])
    gb <- kb / 1024^2
    sprintf("%s: %.2f GB", name, gb)
  }

  cat(tag, "\n", sep = "")
  cat(vapply(vals, to_gb, character(1)), sep = "\n")
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript second_stage_chunk.R <CHUNK_DIR> <CHUNK_INDEX> <OUT_DIR> [WR_RDS] [OVERWRITE]")
}

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2]) + 1
OUT_DIR   <- args[3]
WR_RDS    <- if (length(args) >= 4) args[4] else "/home/schoelleh96/wp22a/data/wrnames.rds"
OVERWRITE <- if (length(args) >= 5) tolower(args[5]) %in% c("true", "1", "yes", "y") else FALSE

change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR, sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) stop("Chunk file not found: ", chunk_file)

outfile <- file.path(OUT_DIR, sprintf("stage2_estimates_chunk_%02d.rds", CHUNK_NO))

log_mem("start")

message("Reading chunk: ", chunk_file)
chunk <- readRDS(chunk_file)

message("Reading WR file: ", WR_RDS)
wr_min <- readRDS(WR_RDS)
wr_min$date <- as.Date(wr_min$date)

available_chunk_vars <- names(chunk)

candidate_predictors <- c("rh_700_850", "upper_wind", "eady", "grad")
predictors <- candidate_predictors[candidate_predictors %in% available_chunk_vars]

if (length(predictors) == 0) {
  stop("None of the requested predictors were found in the chunk. Available vars: ",
       paste(available_chunk_vars, collapse = ", "))
}

message("Predictors used: ", paste(predictors, collapse = ", "))

is_binary01 <- function(x) {
  ux <- sort(unique(x[!is.na(x)]))
  length(ux) > 0 && all(ux %in% c(0, 1))
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

var_reduction <- function(fit0, fit1) {
  v0 <- var(residuals(fit0, type = "response"), na.rm = TRUE)
  v1 <- var(residuals(fit1, type = "response"), na.rm = TRUE)
  if (is.na(v0) || v0 <= 0 || is.na(v1)) return(NA_real_)
  1 - v1 / v0
}

base_df <- data.frame(
  time = as.POSIXct(chunk$time, tz = "UTC"),
  date = as.Date(chunk$time)
) |>
  left_join(wr_min, by = c("date" = "date")) |>
  mutate(
    segment = cut(
      date,
      breaks = c(date[1], change_points, date[length(date)] + 1),
      labels = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )
  ) |>
  group_by(segment) |>
  mutate(day_no = row_number()) |>
  ungroup()

fit_gridpoint <- function(df) {
  d <- df[, c("resid", "segment", "day_no", "wrname", predictors), drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]

  if (nrow(d) < 80) return(NULL)
  if (sd(d$resid) == 0) return(NULL)
  if (length(unique(d$wrname)) < 2) return(NULL)

  d$wrname <- factor(d$wrname)

  fit_wr <- fit_gls(resid ~ wrname, d)
  fit_full <- fit_gls(
    as.formula(paste0("resid ~ wrname + ", paste(predictors, collapse = " + "))),
    d
  )

  emm_wr  <- emmeans(fit_wr,  ~ wrname, data = d, mode = "asymptotic")
  emm_full <- emmeans(fit_full, ~ wrname, data = d, mode = "asymptotic")

  ctr_wr <- emmeans::contrast(emm_wr, "pairwise", data=d)
  ctr_full <- emmeans::contrast(emm_full, "pairwise", data=d)

  # simple scalar summary of how much the wr-pattern shrinks
  ctr_wr_df <- as.data.frame(ctr_wr)
  ctr_full_df <- as.data.frame(ctr_full)

  mean_abs_ctr_wr <- mean(abs(ctr_wr_df$estimate), na.rm = TRUE)
  mean_abs_ctr_full <- mean(abs(ctr_full_df$estimate), na.rm = TRUE)
  wr_contrast_attenuation <- if (is.finite(mean_abs_ctr_wr) && mean_abs_ctr_wr > 0) {
    1 - mean_abs_ctr_full / mean_abs_ctr_wr
  } else {
    NA_real_
  }

  coef_tab <- summary(fit_full)$tTable
  coef_df <- data.frame(
    predictor = rownames(coef_tab),
    estimate = coef_tab[, "Value"],
    se = coef_tab[, "Std.Error"],
    t_value = coef_tab[, "t-value"],
    p_value = coef_tab[, "p-value"],
    row.names = NULL
  )

  list(
    summary = data.frame(
      n = nrow(d),
      var_reduction = var_reduction(fit_wr, fit_full),
      wr_contrast_attenuation = wr_contrast_attenuation,
      n_wr = nlevels(d$wrname)
    ),
    wr_only_emmeans = as.data.frame(emm_wr),
    full_emmeans = as.data.frame(emm_full),
    wr_only_contrasts = ctr_wr_df,
    full_contrasts = ctr_full_df,
    predictor_coefficients = coef_df
  )
}

nx <- length(chunk$lon)
ny <- length(chunk$lat)

summary_out <- list()
wr_only_emm_out <- list()
full_emm_out <- list()
wr_only_con_out <- list()
full_con_out <- list()
coef_out <- list()

k1 <- k2 <- k3 <- k4 <- k5 <- 1L

log_mem("before iteration")

for (ii in seq_len(nx)) {
  for (jj in seq_len(ny)) {

    df <- base_df
    df$lon <- chunk$lon[ii]
    df$lat <- chunk$lat[jj]
    df$i_local <- ii
    df$j_local <- jj
    df$resid <- chunk$residuals[ii, jj, ]

    for (pred in predictors) {
      df[[pred]] <- chunk[[pred]][ii, jj, ]
    }

    ans <- fit_gridpoint(df)
    if (is.null(ans)) next

    ans$summary$lon <- df$lon[1]
    ans$summary$lat <- df$lat[1]
    ans$summary$i_local <- ii
    ans$summary$j_local <- jj
    summary_out[[k1]] <- ans$summary
    k1 <- k1 + 1L

    ans$wr_only_emmeans$lon <- df$lon[1]
    ans$wr_only_emmeans$lat <- df$lat[1]
    ans$wr_only_emmeans$i_local <- ii
    ans$wr_only_emmeans$j_local <- jj
    wr_only_emm_out[[k2]] <- ans$wr_only_emmeans
    k2 <- k2 + 1L

    ans$full_emmeans$lon <- df$lon[1]
    ans$full_emmeans$lat <- df$lat[1]
    ans$full_emmeans$i_local <- ii
    ans$full_emmeans$j_local <- jj
    full_emm_out[[k3]] <- ans$full_emmeans
    k3 <- k3 + 1L

    ans$wr_only_contrasts$lon <- df$lon[1]
    ans$wr_only_contrasts$lat <- df$lat[1]
    ans$wr_only_contrasts$i_local <- ii
    ans$wr_only_contrasts$j_local <- jj
    wr_only_con_out[[k4]] <- ans$wr_only_contrasts
    k4 <- k4 + 1L

    ans$full_contrasts$lon <- df$lon[1]
    ans$full_contrasts$lat <- df$lat[1]
    ans$full_contrasts$i_local <- ii
    ans$full_contrasts$j_local <- jj
    full_con_out[[k5]] <- ans$full_contrasts
    k5 <- k5 + 1L

    coef_df <- ans$predictor_coefficients
    coef_df$lon <- df$lon[1]
    coef_df$lat <- df$lat[1]
    coef_df$i_local <- ii
    coef_df$j_local <- jj
    coef_out[[k5]] <- coef_df
    k5 <- k5 + 1L
  }
}

out <- list(
  metadata = list(
    chunk_no = CHUNK_NO,
    chunk_file = chunk_file,
    source_chunk_metadata = chunk$metadata,
    wr_rds = WR_RDS,
    predictors = predictors
  ),
  summary = bind_rows(summary_out),
  wr_only_emmeans = bind_rows(wr_only_emm_out),
  full_emmeans = bind_rows(full_emm_out),
  wr_only_contrasts = bind_rows(wr_only_con_out),
  full_contrasts = bind_rows(full_con_out),
  predictor_coefficients = bind_rows(coef_out)
)

if (file.exists(outfile) && !OVERWRITE) {
  message("Output exists and OVERWRITE=FALSE, skipping write: ", outfile)
} else {
  saveRDS(out, outfile, compress = "xz")
  message("Wrote: ", outfile)
}