#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript collect_wcb_plot_data_parallel.R <CHUNK_DIR> <OUT_FILE>")
}

CHUNK_DIR <- args[1]
OUT_FILE  <- args[2]

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

chunk_files <- list.files(
  CHUNK_DIR,
  pattern = "^wcb_stage2_chunk_[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(chunk_files) == 0) {
  stop("No chunk files found in: ", CHUNK_DIR)
}

# keep files in numeric order
chunk_ids <- as.integer(sub("^.*chunk_([0-9]+)\\.rds$", "\\1", basename(chunk_files)))
chunk_files <- chunk_files[order(chunk_ids)]

mc_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (!is.finite(mc_cores) || mc_cores < 1) mc_cores <- 1
mc_cores <- min(mc_cores, length(chunk_files))

sig_class <- function(effect, p, alpha = 0.05) {
  fifelse(
    !is.na(p) & p < alpha & effect > 0, "pos_sig",
    fifelse(!is.na(p) & p < alpha & effect < 0, "neg_sig", "ns")
  )
}

extract_wr_full <- function(chunk, chunk_file) {
  if (is.null(chunk$wr_means) || !nrow(chunk$wr_means)) return(NULL)

  x <- as.data.table(chunk$wr_means)

needed <- c(
  "wrname", "emmean_wr", "SE_wr", "asymp.LCL_wr", "asymp.UCL_wr", "p.value_wr",
  "emmean_full", "SE_full", "asymp.LCL_full", "asymp.UCL_full", "p.value_full",
  "delta_emmean", "regime_attenuation",
  "lon", "lat", "i_local", "j_local"
)
  if (!all(needed %in% names(x))) return(NULL)

  se_delta <- sqrt(x$SE_wr^2 + x$SE_full^2)
  delta_emmean <- x$emmean_full - x$emmean_wr
  z_delta <- delta_emmean / se_delta
  p_delta <- 2 * pnorm(-abs(z_delta))

  x[, .(
    source_file = chunk_file,
    wrname,
    lon, lat, i_local, j_local,

    effect_wr = emmean_wr,
    se_wr = SE_wr,
    lcl_wr = asymp.LCL_wr,
    ucl_wr = asymp.UCL_wr,
    p_value_wr = p.value_wr,
    sign_class_wr = sig_class(emmean_wr, p.value_wr),

    effect_full = emmean_full,
    se_full = SE_full,
    lcl_full = asymp.LCL_full,
    ucl_full = asymp.UCL_full,
    p_value_full = p.value_full,
    sign_class_full = sig_class(emmean_full, p.value_full),

    delta_emmean = delta_emmean,
    se_delta = se_delta,
    z_delta = z_delta,
    p_delta = p_delta,
    significant_delta = !is.na(p_delta) & p_delta < 0.05,

    regime_attenuation = regime_attenuation
  )]
}

extract_wcb_effects <- function(chunk, chunk_file) {
  if (is.null(chunk$wcb_effects) || !nrow(chunk$wcb_effects)) return(NULL)

  x <- as.data.table(chunk$wcb_effects)
  if (!("predictor" %in% names(x))) return(NULL)

  # x <- x[predictor %in% c("in_24", "asc_12", "out_00")]
  x <- x[predictor %in% c("inf_12", "asc_00", "out_12")]

  if (!nrow(x)) return(NULL)

  if (!"p.value" %in% names(x)) x[, p.value := NA_real_]
  if (!"sd_predictor" %in% names(x)) x[, sd_predictor := NA_real_]

  # Normalize the effect column to "effect_1sd"
  x[, effect_unit := NA_real_]

  if ("in_24_c.trend" %in% names(x))  x[predictor == "in_24",  effect_unit := `in_24_c.trend`]
  if ("asc_12_c.trend" %in% names(x)) x[predictor == "asc_12", effect_unit := `asc_12_c.trend`]
  if ("out_00_c.trend" %in% names(x)) x[predictor == "out_00", effect_unit := `out_00_c.trend`]

  # fall back to scaled_effect if present
  if ("scaled_effect" %in% names(x)) {
    x[, effect_1sd := scaled_effect]
  } else {
    x[, effect_1sd := effect_unit * sd_predictor]
  }

  if ("se_1sd" %in% names(x)) {
    # keep as is
  } else if ("SE" %in% names(x) && "sd_predictor" %in% names(x)) {
    x[, se_1sd := SE * sd_predictor]
  } else {
    x[, se_1sd := NA_real_]
  }

  if (!"lcl_1sd" %in% names(x)) {
    if ("asymp.LCL" %in% names(x) && "sd_predictor" %in% names(x)) {
      x[, lcl_1sd := asymp.LCL * sd_predictor]
    } else {
      x[, lcl_1sd := NA_real_]
    }
  }

  if (!"ucl_1sd" %in% names(x)) {
    if ("asymp.UCL" %in% names(x) && "sd_predictor" %in% names(x)) {
      x[, ucl_1sd := asymp.UCL * sd_predictor]
    } else {
      x[, ucl_1sd := NA_real_]
    }
  }

  x[, `:=`(
    source_file = chunk_file,
    significant = !is.na(p.value) & p.value < 0.05,
    sign_class = sig_class(effect_1sd, p.value)
  )]

  keep <- c(
    "source_file", "lon", "lat", "i_local", "j_local",
    "wrname", "predictor",
    "effect_unit", "effect_1sd", "se_1sd", "lcl_1sd", "ucl_1sd",
    "p.value", "significant", "sign_class", "sd_predictor"
  )
  keep <- keep[keep %in% names(x)]
  x[, ..keep]
}
combine_three_predictors <- function(wcb_long) {
  if (is.null(wcb_long) || !nrow(wcb_long)) return(NULL)

  wcb_long <- as.data.table(wcb_long)
  # wcb_long <- wcb_long[predictor %in% c("in_24", "asc_12", "out_00")]
  wcb_long <- wcb_long[predictor %in% c("inf_12", "asc_00", "out_12")
]

  if (!nrow(wcb_long)) return(NULL)

  # compute combined effect and propagated standard error
  wcb_comb <- wcb_long[, .(
    combined_effect_1sd = sum(effect_1sd, na.rm = TRUE),
    se_combined = sqrt(sum(se_1sd^2, na.rm = TRUE)),
    n_predictors = sum(!is.na(effect_1sd)),
    all_component_sig = all(significant, na.rm = TRUE),
    same_sign = {
      s <- sign(effect_1sd[!is.na(effect_1sd)])
      length(unique(s[s != 0])) <= 1
    }
  ), by = .(source_file, lon, lat, i_local, j_local, wrname)]

  # calculate z-score and p-value for the combined effect
  wcb_comb[, z_combined := combined_effect_1sd / se_combined]
  wcb_comb[, p_combined := 2 * pnorm(-abs(z_combined))]
  wcb_comb[, significant_combined := p_combined < 0.05]

  # retain the previous combined screening flag
  wcb_comb[, combined_screen_sig := (n_predictors == 3L) & all_component_sig & same_sign]

  wcb_comb[]
}

process_one <- function(f) {
  chunk <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(chunk)) return(NULL)

  wr_part <- tryCatch(extract_wr_full(chunk, f), error = function(e) NULL)
  wcb_part <- tryCatch(extract_wcb_effects(chunk, f), error = function(e) NULL)

  wcb_comb <- tryCatch(combine_three_predictors(wcb_part), error = function(e) NULL)

  list(
    wr_full_map = wr_part,
    wcb_effect_map = wcb_part,
    wcb_combined_map = wcb_comb
  )
}

message("Processing ", length(chunk_files), " chunk files with ", mc_cores, " cores.")
res <- mclapply(chunk_files, process_one, mc.cores = mc_cores, mc.preschedule = FALSE)

wr_full_map <- rbindlist(lapply(res, `[[`, "wr_full_map"), fill = TRUE)
wcb_effect_map <- rbindlist(lapply(res, `[[`, "wcb_effect_map"), fill = TRUE)
wcb_combined_map <- rbindlist(lapply(res, `[[`, "wcb_combined_map"), fill = TRUE)

out <- list(
  metadata = list(
    chunk_dir = CHUNK_DIR,
    n_chunks = length(chunk_files),
    n_cores = mc_cores
  ),
  wr_full_map = wr_full_map,
  wcb_effect_map = wcb_effect_map,
  wcb_combined_map = wcb_combined_map
)

saveRDS(out, OUT_FILE, compress = "xz")
message("Wrote: ", OUT_FILE)