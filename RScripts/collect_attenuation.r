#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript collect_stage2_attenuation.R <CHUNK_DIR> <OUT_FILE> <N_CORES>")
}

CHUNK_DIR <- args[1]
OUT_FILE  <- args[2]
N_CORES   <- as.integer(args[3])
if(is.na(N_CORES) || N_CORES < 1) N_CORES <- 1

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

chunk_files <- list.files(
  CHUNK_DIR,
  pattern = "^stage2_estimates_chunk_[0-9]+\\.rds$",
  full.names = TRUE
)
if(length(chunk_files) == 0) stop("No chunk files found in: ", CHUNK_DIR)

chunk_ids <- as.integer(sub("^.*chunk_([0-9]+)\\.rds$", "\\1", basename(chunk_files)))
chunk_files <- chunk_files[order(chunk_ids)]

mc_cores <- min(N_CORES, length(chunk_files))
message("Using ", mc_cores, " cores for parallel collection.")

setDTthreads(1)  # prevent data.table from multi-threading

sig_class <- function(effect, p, alpha = 0.05) {
  fifelse(
    !is.na(p) & p < alpha & effect > 0, "pos_sig",
    fifelse(!is.na(p) & p < alpha & effect < 0, "neg_sig", "ns")
  )
}

extract_summary <- function(chunk, chunk_file) {
  if (is.null(chunk$summary) || !nrow(chunk$summary)) return(NULL)

  x <- as.data.table(chunk$summary)

  req <- c(
    "n", "joint_wr_pattern_var_reduction", "joint_wr_pattern_mad_reduction", 
    "lon", "lat", "i_local", "j_local"
  )
  if (!all(req %in% names(x))) return(NULL)
  x[, .(
    n, lon, lat, i_local, j_local,
    joint_wr_pattern_var_reduction, joint_wr_pattern_mad_reduction
  )]
}

extract_wr_attenuation <- function(chunk, chunk_file) {
  if (is.null(chunk$wr_only_emmeans) || !nrow(chunk$wr_only_emmeans)) return(NULL)
  if (is.null(chunk$full_emmeans) || !nrow(chunk$full_emmeans)) return(NULL)

  wr0 <- as.data.table(chunk$wr_only_emmeans)
  wr1 <- as.data.table(chunk$full_emmeans)

  req <- c("wrname", "emmean", "SE", "asymp.LCL", "asymp.UCL", "p.value",
           "lon", "lat", "i_local", "j_local")
  if (!all(req %in% names(wr0))) return(NULL)
  if (!all(req %in% names(wr1))) return(NULL)

  wr0 <- wr0[, .(
    wrname, lon, lat, i_local, j_local,
    emmean_wr = emmean,
    SE_wr = SE,
    lcl_wr = asymp.LCL,
    ucl_wr = asymp.UCL,
    p_value_wr = p.value
  )]

  wr1 <- wr1[, .(
    wrname, lon, lat, i_local, j_local,
    emmean_full = emmean,
    SE_full = SE,
    lcl_full = asymp.LCL,
    ucl_full = asymp.UCL,
    p_value_full = p.value
  )]

  m <- merge(
    wr0, wr1,
    by = c("lon", "lat", "i_local", "j_local", "wrname"),
    all = FALSE
  )

  m[, `:=`(
    delta_abs = emmean_full - emmean_wr,
    delta_rel = ifelse(abs(emmean_wr) > 0, 1 - abs(emmean_full) / abs(emmean_wr), NA_real_),
    se_delta = sqrt(SE_wr^2 + SE_full^2)
  )]

  m[, z_delta := delta_abs / se_delta]
  m[, p_delta := 2 * pnorm(-abs(z_delta))]
  m[, p_delta_fdr := p.adjust(p_delta, method = "fdr"), by = wrname]
  m[, significant_delta := !is.na(p_delta_fdr) & p_delta_fdr < 0.05]
  m[, sign_class_delta := sig_class(delta_abs, p_delta_fdr)]

  m[, source_file := chunk_file]

  m[, .(
    source_file, wrname, lon, lat, i_local, j_local,
    emmean_wr, SE_wr, lcl_wr, ucl_wr, p_value_wr,
    emmean_full, SE_full, lcl_full, ucl_full, p_value_full,
    delta_abs, delta_rel, se_delta, z_delta, p_delta, p_delta_fdr,
    significant_delta, sign_class_delta
  )]
}

extract_predictor_attenuation <- function(chunk, chunk_file) {
  if (is.null(chunk$predictor_attenuation) || !nrow(chunk$predictor_attenuation)) return(NULL)

  x <- as.data.table(chunk$predictor_attenuation)

  req <- c(
    "predictor", "n", "wr_pattern_var_reduction", "wr_pattern_mad_reduction",
    "wr_only_mean_abs_dev", "pred_mean_abs_dev", "lon", "lat", "i_local", "j_local"
  )
  if (!all(req %in% names(x))) return(NULL)

  x[, `:=`(
    mad_abs_reduction = wr_only_mean_abs_dev - pred_mean_abs_dev,
    mad_rel_reduction = ifelse(
      abs(wr_only_mean_abs_dev) > 0,
      1 - pred_mean_abs_dev / wr_only_mean_abs_dev,
      NA_real_
    ),
    source_file = chunk_file
  )]

  x[, .(
    source_file, predictor, n, lon, lat, i_local, j_local,
    wr_only_mean_abs_dev, pred_mean_abs_dev,
    mad_abs_reduction, mad_rel_reduction,
    wr_pattern_var_reduction, wr_pattern_mad_reduction
  )]
}

process_one <- function(f) {
  chunk <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(chunk)) return(NULL)

  list(
    summary = tryCatch(extract_summary(chunk, f), error = function(e) NULL),
    wr_attenuation = tryCatch(extract_wr_attenuation(chunk, f), error = function(e) NULL),
    predictor_attenuation = tryCatch(extract_predictor_attenuation(chunk, f), error = function(e) NULL)
  )
}

message("Processing ", length(chunk_files), " chunk files with ", mc_cores, " cores.")


res <- mclapply(chunk_files, process_one, mc.cores = mc_cores, mc.preschedule = FALSE)

summary <- rbindlist(lapply(res, `[[`, "summary"), fill = TRUE)
wr_attenuation <- rbindlist(lapply(res, `[[`, "wr_attenuation"), fill = TRUE)
predictor_attenuation <- rbindlist(lapply(res, `[[`, "predictor_attenuation"), fill = TRUE)

out <- list(
  metadata = list(
    chunk_dir = CHUNK_DIR,
    n_chunks = length(chunk_files),
    n_cores = mc_cores
  ),
  summary = summary,
  wr_attenuation = wr_attenuation,
  predictor_attenuation = predictor_attenuation
)

saveRDS(out, OUT_FILE, compress = "xz")
message("Wrote: ", OUT_FILE)