#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

IN_DIR  <- get_arg("--in-dir",  Sys.getenv("IN_DIR",  "./wcb_wr_gls_results"))
OUT_RDS <- get_arg("--out-rds", Sys.getenv("OUT_RDS", "./combined_wcb_wr_effects.rds"))
PATTERN <- get_arg("--pattern", Sys.getenv("PATTERN", "^chunk_[0-9]+\\.rds$"))
JOBS    <- as.integer(get_arg("--jobs", Sys.getenv("JOBS", max(1L, detectCores(logical = FALSE) - 1L))))

if (!dir.exists(IN_DIR)) {
  stop("Input directory does not exist: ", IN_DIR)
}

files <- list.files(IN_DIR, pattern = PATTERN, full.names = TRUE)
files <- sort(files)

if (length(files) == 0L) {
  stop("No chunk result files found in: ", IN_DIR)
}

message("Found ", length(files), " result files.")
message("Using ", JOBS, " worker(s).")

extract_one <- function(f) {
  obj <- readRDS(f)

  meta <- obj$metadata %||% list()
  chunk_file <- meta$chunk_file %||% f
  task_id <- meta$task_id %||% NA_integer_

  eff <- obj$effect_table
  if (!is.null(eff) && nrow(eff) > 0L) {
    eff$source_file <- basename(f)
    eff$chunk_file  <- chunk_file
    eff$task_id     <- task_id
  }

  sk <- obj$skipped_cells
  if (!is.null(sk) && nrow(sk) > 0L) {
    sk$source_file <- basename(f)
    sk$chunk_file  <- chunk_file
    sk$task_id     <- task_id
  }

  list(
    effect_table = eff,
    skipped_cells = sk,
    metadata = meta
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

res_list <- if (.Platform$OS.type == "windows" || JOBS <= 1L) {
  lapply(files, extract_one)
} else {
  mclapply(files, extract_one, mc.cores = JOBS, mc.preschedule = FALSE)
}

combine_dfs <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0L) return(NULL)

  cols <- unique(unlist(lapply(lst, names), use.names = FALSE))

  align <- function(df) {
    miss <- setdiff(cols, names(df))
    for (nm in miss) df[[nm]] <- NA
    df[cols]
  }

  lst <- lapply(lst, align)
  do.call(rbind, lst)
}

effect_table <- combine_dfs(lapply(res_list, `[[`, "effect_table"))
skipped_cells <- combine_dfs(lapply(res_list, `[[`, "skipped_cells"))

out <- list(
  metadata = list(
    input_dir = IN_DIR,
    n_files = length(files),
    files = basename(files),
    created = Sys.time()
  ),
  effect_table = effect_table,
  skipped_cells = skipped_cells
)

saveRDS(out, OUT_RDS, compress = "xz")
message("Wrote combined object: ", OUT_RDS)