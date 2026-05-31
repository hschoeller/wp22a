#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript extract_wr_emmeans.R <STAGE2_DIR> <OUTFILE> [N_CORES]")
}

stage2_dir <- args[1]
outfile    <- args[2]
n_cores    <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores() - 1L

chunk_files <- list.files(
    stage2_dir,
    pattern    = "^stage2_estimates_chunk_\\d+\\.rds$",
    full.names = TRUE
)

if (length(chunk_files) == 0) {
    stop("No stage2 chunk files found in: ", stage2_dir)
}

message(
    "Found ", length(chunk_files), " chunk files. ",
    "Using ", n_cores, " cores."
)

extract_chunk_emmeans <- function(chunk_file) {
    chunk_out <- tryCatch(
        readRDS(chunk_file),
        error = function(e) {
            message("Failed to read: ", chunk_file, " — ", e$message)
            return(NULL)
        }
    )

    if (is.null(chunk_out)) return(NULL)

    emm <- chunk_out$wr_only_emmeans
    if (is.null(emm) || nrow(emm) == 0) return(NULL)

    dplyr::transmute(
        emm,
        lon, lat, wrname,
        emmean,
        p_value = 2 * pnorm(-abs(emmean / SE))
    )
}

emmeans_list <- parallel::mclapply(
    chunk_files,
    extract_chunk_emmeans,
    mc.cores     = n_cores,
    mc.preschedule = TRUE
)

emmeans_df <- dplyr::bind_rows(emmeans_list) |>
    dplyr::distinct(lon, lat, wrname, .keep_all = TRUE) |>
    dplyr::group_by(wrname) |>
    dplyr::mutate(p_value_adj = p.adjust(p_value, method = "fdr")) |>
    dplyr::ungroup()

message(
    "Extracted ", nrow(emmeans_df), " rows across ",
    dplyr::n_distinct(emmeans_df[c("lon", "lat")]), " grid points and ",
    dplyr::n_distinct(emmeans_df$wrname), " WR levels."
)

saveRDS(emmeans_df, outfile, compress = "xz")
message("Saved to: ", outfile)