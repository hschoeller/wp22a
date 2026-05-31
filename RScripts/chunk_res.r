#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

infile  <- "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
outdir  <- "/scratch/schoelleh96/wp22a/stage2_chunks_res"
lat_chunks <- 11
lon_chunks <- 241

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(infile)

time <- as.Date(x$time)
keep <- time >= as.Date("1979-01-01")

resid <- x$residuals[, , keep, drop = FALSE]
time  <- time[keep]
lon   <- x$lon
lat   <- x$lat

split_indices <- function(n, k) {
  k <- max(1, min(k, n))
  edges <- floor(seq(0, n, length.out = k + 1))
  lapply(seq_len(k), function(i) seq.int(edges[i] + 1, edges[i + 1]))
}

lon_idx <- split_indices(length(lon), lon_chunks)
lat_idx <- split_indices(length(lat), lat_chunks)

chunk_id <- 1L
for (iy in seq_along(lat_idx)) {
  for (ix in seq_along(lon_idx)) {
    jj <- lat_idx[[iy]]
    ii <- lon_idx[[ix]]

    out <- list(
      residuals = resid[ii, jj, , drop = FALSE],
      lon = lon[ii],
      lat = lat[jj],
      time = time,
      metadata = list(
        source_file = infile,
        time_start = as.character(min(time)),
        time_end = as.character(max(time)),
        lon_chunk = ix,
        lat_chunk = iy,
        chunk_id = chunk_id
      )
    )

    saveRDS(out, file.path(outdir, sprintf("chunk_%02d.rds", chunk_id)), compress = "xz")
    message("Wrote chunk", chunk_id + 1L)
    chunk_id <- chunk_id + 1L
  }
}

message("Wrote ", chunk_id - 1L, " chunks to ", outdir)