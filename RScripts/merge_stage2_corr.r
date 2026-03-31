#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript merge_stage2Corrs.r <CHUNK_CORR_DIR> <RESIDUAL_RDS> <OUT_DIR>")
}

CHUNK_CORR_DIR <- args[1]
RESID_RDS      <- args[2]
OUT_DIR        <- args[3]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(RESID_RDS)
global_lon <- x$lon
global_lat <- x$lat

chunk_files <- sort(list.files(CHUNK_CORR_DIR, pattern = "^stage2_corr_chunk_[0-9]+\\.rds$", full.names = TRUE))
if (!length(chunk_files)) stop("No chunk correlation files found in ", CHUNK_CORR_DIR)

# discover variable names
first <- readRDS(chunk_files[1])
vars <- names(first$correlations)
print(vars)
# initialize full maps
full_maps <- lapply(vars, function(v) {
  array(NA_real_, dim = c(length(global_lon), length(global_lat)))
})
names(full_maps) <- vars

for (f in chunk_files) {
  z <- readRDS(f)

  ix <- z$metadata$global_lon_index
  iy <- z$metadata$global_lat_index

  for (v in names(z$correlations)) {
    full_maps[[v]][ix, iy] <- z$correlations[[v]]
  }
}

for (v in names(full_maps)) {
  saveRDS(
    list(
      variable = v,
      correlation = full_maps[[v]],
      lon = global_lon,
      lat = global_lat
    ),
    file.path(OUT_DIR, paste0(v, "_corr_map.rds")),
    compress = "xz"
  )
  message("Wrote: ", file.path(OUT_DIR, paste0(v, "_corr_map.rds")))
}