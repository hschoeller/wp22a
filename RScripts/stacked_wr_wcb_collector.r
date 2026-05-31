#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

IN_DIR  <- get_arg("--in-dir",  "./stage3_results")
OUT_RDS <- get_arg("--out-rds", file.path(IN_DIR, "combined_results.rds"))

files <- list.files(IN_DIR, pattern = "^chunk_.*\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No chunk files found in: ", IN_DIR)
}

message("Found ", length(files), " chunk files")

# Collect any data.frame-like component if present and non-empty
collect_field <- function(obj, field) {
  x <- obj[[field]]
  if (is.null(x)) return(NULL)
  if (!is.data.frame(x)) return(NULL)
  if (nrow(x) == 0L) return(NULL)
  x
}

bind_or_null <- function(x) {
  if (length(x) > 0L) do.call(rbind, x) else NULL
}

wcb_list    <- list()
wr_list     <- list()
summary_list <- list()

k_wcb <- 0L
k_wr  <- 0L
k_sum <- 0L

for (f in files) {
  message("Reading: ", basename(f))

  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) next

  # New chunk layout:
  # $wcb_effects
  # $wr_reduction
  # $wr_summary

  x <- collect_field(obj, "wcb_effects")
  if (!is.null(x)) {
    k_wcb <- k_wcb + 1L
    wcb_list[[k_wcb]] <- x
  }

  x <- collect_field(obj, "wr_reduction")
  if (!is.null(x)) {
    k_wr <- k_wr + 1L
    wr_list[[k_wr]] <- x
  }

  x <- collect_field(obj, "wr_summary")
  if (!is.null(x)) {
    k_sum <- k_sum + 1L
    summary_list[[k_sum]] <- x
  }
}

message("Combining results...")

wcb_all     <- bind_or_null(wcb_list)
wr_all      <- bind_or_null(wr_list)
summary_all <- bind_or_null(summary_list)

out <- list(
  wcb_effects = wcb_all,
  wr_reduction = wr_all,
  wr_summary   = summary_all
)

saveRDS(out, OUT_RDS, compress = "xz")

message("Saved combined results to: ", OUT_RDS)
message("WCB rows: ", if (!is.null(wcb_all)) nrow(wcb_all) else 0)
message("WR rows: ", if (!is.null(wr_all)) nrow(wr_all) else 0)
message("Summary rows: ", if (!is.null(summary_all)) nrow(summary_all) else 0)