#!/usr/bin/env Rscript

source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript compute_composites_from_cube_parallel.R <cube_rds> <wr_rds> <output_rds> <n_perm> [ncores] [perm_chunk]"
  )
}

cube_file  <- args[[1]]
wr_file    <- args[[2]]
out_file   <- args[[3]]
n_perm     <- as.integer(args[[4]])
ncores     <- if (length(args) >= 5) as.integer(args[[5]]) else 24
perm_chunk <- if (length(args) >= 6) as.integer(args[[6]]) else 50

split_perm_indices <- function(n_perm, chunk_size) {
  split(seq_len(n_perm), ceiling(seq_len(n_perm) / chunk_size))
}

cube <- readRDS(cube_file)
wrdf <- readRDS(wr_file)

print("Files read:")
print(str(cube))
print(str(wrdf))

lon <- cube$lon
lat <- cube$lat

arr <- cube$residuals
dims <- dim(arr)
nlon <- dims[1]
nlat <- dims[2]
ntime <- dims[3]
ngrid <- nlon * nlat

# reshape full cube to [time, grid]
Rmat_full <- matrix(aperm(arr, c(3, 1, 2)), nrow = ntime, ncol = ngrid)
grid <- expand.grid(lon = lon, lat = lat)

wr_levels <- levels(factor(wrdf$wrname))
print(wr_levels)

cube_dates <- as.Date(cube$time)
wrdf <- wrdf %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

common_dates <- sort(cube_dates[cube_dates %in% wrdf$date])

if (length(common_dates) == 0) {
  stop("No overlapping dates between residual cube and WR time series")
}

cube_keep <- cube_dates %in% common_dates
Rmat <- Rmat_full[cube_keep, , drop = FALSE]
full_dates <- cube_dates[cube_keep]

wrdf <- wrdf %>%
  filter(date %in% common_dates)

process_one_category <- function(wr) {
  actual_dates <- sort(unique(wrdf$date[wrdf$wrname == wr]))
  obs_idx <- which(full_dates %in% actual_dates)
  blocks <- get_run_blocks(full_dates, actual_dates)

  obs_mean <- colMeans(Rmat[obs_idx, , drop = FALSE], na.rm = TRUE)

  perm_groups <- split_perm_indices(n_perm, perm_chunk)
  print("Now in parallel")
#   chunk_counts <- lapply(
#   seq_along(perm_groups),
#   function(i) {
#     pg <- perm_groups[[i]]
#     message("Processing permutation chunk ", i, " / ", length(perm_groups))

#     cnt <- integer(ngrid)

#     for (p in pg) {
#       surr_dates <- permute_blocks(blocks, full_dates)
#       surr_idx <- which(full_dates %in% surr_dates)
#       print(head(surr_idx))
# print(range(surr_idx))
# print(is.integer(surr_idx))
#       surr_mean <- colMeans(Rmat[surr_idx, , drop = FALSE], na.rm = TRUE)
#       cnt <- cnt + (abs(surr_mean) >= abs(obs_mean))
#     }

#     cnt
#   }
# )
chunk_counts <- mclapply(
  seq_along(perm_groups),
  function(i) {
    pg <- perm_groups[[i]]

    cnt <- integer(ngrid)

    for (p in pg) {
      surr_dates <- permute_blocks(blocks, full_dates)
      surr_idx <- which(full_dates %in% surr_dates)
      surr_mean <- colMeans(Rmat[surr_idx, , drop = FALSE], na.rm = TRUE)
      cnt <- cnt + (abs(surr_mean) >= abs(obs_mean))
    }

    cnt
  },
  mc.cores = min(ncores, length(perm_groups))
)

  extreme_count <- Reduce(`+`, chunk_counts)
  p_value <- (extreme_count + 1) / (n_perm + 1)

  print(str(grid))
  print(str(wr))
  print(str(obs_mean))
  print(str(p_value))
  data.frame(
    lon = grid$lon,
    lat = grid$lat,
    wrname = wr,
    composite_mean = obs_mean,
    p_value = p_value
  )
}

print("Calculating")

results_list <- lapply(wr_levels, process_one_category)
results <- bind_rows(results_list)
results <- results %>% group_by(wrname) %>% mutate(p_value_adj = p.adjust(p_value, method = "fdr")) %>% ungroup()

saveRDS(results, out_file)
cat("Saved:", out_file, "\n")