#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(emmeans)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript extract_wr_effects_emmeans.R <model_dir> <output_rds> [ncores]"
  )
}

model_dir  <- args[[1]]
output_rds <- args[[2]]
ncores     <- if (length(args) >= 3) as.integer(args[[3]]) else 24L

extract_lon_lat <- function(path) {
  nm <- basename(path)
  nm <- sub("\\.Rds$", "", nm)
  nm <- sub("^lm", "", nm)
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  c(lat = as.numeric(parts[1]), lon = as.numeric(parts[2]))
}

process_one_model <- function(fpath) {
  mod <- readRDS(fpath)
  crd <- extract_lon_lat(fpath)
  mod$call$model <- formula(mod)
  emm <- emmeans(mod, ~ wrname, data=mod$data, mode = "df.error")
  emm_df <- as.data.frame(summary(emm))

  # Effect of each WR level relative to the average over WR levels
  eff <- contrast(emm, method = "eff")
  eff_df <- as.data.frame(summary(eff, infer = c(TRUE, TRUE)))

  # Try to recover wrname cleanly from contrast labels
  # For "eff" contrasts, labels are usually just the level names
  wr_levels <- as.character(emm_df$wrname)
  if (nrow(eff_df) == length(wr_levels)) {
    eff_df$wrname <- wr_levels
  } else {
    eff_df$wrname <- eff_df$contrast
  }

  data.frame(
    lon = unname(crd["lon"]),
    lat = unname(crd["lat"]),
    wrname = eff_df$wrname,
    composite_mean = eff_df$estimate,   # analogue to composite anomaly
    SE = eff_df$SE,
    df = eff_df$df,
    t_ratio = eff_df$t.ratio,
    p_value = eff_df$p.value,
    emmean = emm_df$emmean[match(eff_df$wrname, as.character(emm_df$wrname))]
  )
}

files <- list.files(model_dir, pattern = "\\.Rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("No .Rds files found in ", model_dir)
}

message("Found ", length(files), " model files")
message("Using ", ncores, " cores")

# res_list <- lapply(files[1:10], process_one_model)

res_list <- mclapply(
  files,
  process_one_model,
  mc.cores = ncores
)

results <- do.call(rbind, res_list)
results <- results %>% group_by(wrname) %>% mutate(p_value_adj = p.adjust(p_value, method = "fdr")) %>% ungroup()

saveRDS(results, output_rds)
message("Saved: ", output_rds)