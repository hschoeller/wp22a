library(data.table)

chunk_dir <- "/scratch/schoelleh96/wp22a/stage1_models_perm/"
output_file <- "/home/schoelleh96/wp22a/ens_data/perm_models_collect.rds"

summary_cols <- c(
    "lon", "lat", "n", "n_wr", "hac_bandwidth",
    "pred_wald_p",
    "S_raw", "S_full", "deltaS",
    "lowerCL_deltaS", "upperCL_deltaS", "deltaS_p"
)

emmeans_cols <- c(
    "lon", "lat", "wrname",
    "emmean_raw", "se_raw","p_raw",
    "lowerCL_raw", "upperCL_raw",
    "emmean_full", "se_full","p_full",
    "lowerCL_full", "upperCL_full",
    "delta_emmean", "se_delta_emmean", "p_delta_emmean", 
    "lowerCL_delta", "lowerCL_delta"
)

loo_attr_cols <- c(
    "lon", "lat", "wrname", "predictor", "attr_emmean", "attr_p"
)

loo_rms_cols <- c(
    "lon", "lat", "predictor", "rms_attr", "rms_attr_p"
)

pred_eff_cols <- c(
    "lon", "lat", "predictor", "p_value"
)

chunk_files <- list.files(
    chunk_dir,
    pattern = "^stage2_estimates_chunk_\\d+\\.rds$",
    full.names = TRUE
)

# Allocate containers
summaries     <- vector("list", length(chunk_files))
wr_emmeans    <- vector("list", length(chunk_files))
loo_attrs     <- vector("list", length(chunk_files))
loo_rms_list  <- vector("list", length(chunk_files))
pred_effs  <- vector("list", length(chunk_files))

for (i in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[[i]])

    summaries[[i]]    <- chunk$summary[, summary_cols, drop = FALSE]
    wr_emmeans[[i]]   <- chunk$wr_emmeans[, emmeans_cols, drop = FALSE]

    loo_attrs[[i]]    <- chunk$loo_attr[, loo_attr_cols, drop = FALSE]
    loo_rms_list[[i]] <- chunk$loo_rms[, loo_rms_cols, drop = FALSE]
    pred_effs[[i]] <- chunk$pred_eff[, pred_eff_cols, drop = FALSE]

}

# Bind safely
result <- list(
    summary      = rbindlist(summaries, use.names = TRUE, fill = TRUE),
    wr_emmeans   = rbindlist(wr_emmeans, use.names = TRUE, fill = TRUE),
    loo_attr     = rbindlist(loo_attrs, use.names = TRUE, fill = TRUE),
    loo_rms      = rbindlist(loo_rms_list, use.names = TRUE, fill = TRUE),
    pred_effs   = rbindlist(pred_effs, use.names = TRUE, fill = TRUE)
)

saveRDS(result, output_file)