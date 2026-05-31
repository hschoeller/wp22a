#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ncdf4))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

TASK_ID <- as.integer(get_arg("--task-id",
                    Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

IN_DIR  <- get_arg("--in-dir",
            Sys.getenv("IN_DIR", "./residual_chunks"))

OUT_DIR <- get_arg("--out-dir",
            Sys.getenv("OUT_DIR", "./stage2_chunks"))

WCB_ROOT <- get_arg("--wcb-root",
             Sys.getenv("WCB_ROOT",
               "/scratch/schoelleh96/wp22a/ELIAS_data"))

OVERWRITE <- as.logical(as.integer(
  get_arg("--overwrite", Sys.getenv("OVERWRITE", "0"))
))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(IN_DIR, sprintf("chunk_%02d.rds", TASK_ID))
out_file   <- file.path(OUT_DIR, sprintf("chunk_%02d.rds", TASK_ID))

if (!file.exists(chunk_file)) {
  stop("Residual chunk not found: ", chunk_file)
}

if (file.exists(out_file) && !OVERWRITE) {
  message("Chunk exists and OVERWRITE=FALSE, skipping: ", out_file)
  quit(save = "no", status = 0L)
}

message("Reading residual chunk: ", chunk_file)
chunk <- readRDS(chunk_file)

if (is.null(chunk$residuals) || is.null(chunk$lon) || is.null(chunk$lat) || is.null(chunk$time)) {
  stop("Input chunk must contain residuals, lon, lat, and time.")
}

resid_arr <- chunk$residuals
lon_vals  <- as.numeric(chunk$lon)
lat_vals  <- as.numeric(chunk$lat)
time_vals <- as.POSIXct(chunk$time, tz = "UTC") + 12L * 3600L

find_wcb_file <- function(root, dt) {
  ddir <- file.path(root, format(dt, "%Y", tz = "UTC"),
                          format(dt, "%m", tz = "UTC"))
  base <- sprintf("hit_%s_%s",
                  format(dt, "%Y%m%d", tz = "UTC"),
                  format(dt, "%H",     tz = "UTC"))
  cand <- file.path(ddir, c(base, paste0(base, ".nc")))
  hit  <- cand[file.exists(cand)]
  if (length(hit) > 0L) return(hit[1L])
  alt <- Sys.glob(file.path(ddir, paste0(base, "*")))
  if (length(alt) > 0L) return(alt[1L])
  NA_character_
}

get_dim_vals <- function(nc, candidates) {
  for (nm in candidates) {
    if (!is.null(nc$dim[[nm]]) && !is.null(nc$dim[[nm]]$vals)) {
      return(nc$dim[[nm]]$vals)
    }
  }
  NULL
}

read_wcb_snapshot <- function(nc_path, var_name, ix_wcb, iy_wcb) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)
  raw <- ncvar_get(nc, var_name)  # assumed [lon, lat]
  raw[ix_wcb, iy_wcb]
}

# Find WCB files for each residual timestamp, then keep only the intersection in time.
wcb_files_all <- lapply(time_vals, find_wcb_file, root = WCB_ROOT)
keep_time <- !is.na(unlist(wcb_files_all))

if (!any(keep_time)) {
  stop("No WCB files found for any timestamps in this residual chunk.")
}

time_vals <- time_vals[keep_time]
resid_arr <- resid_arr[, , keep_time, drop = FALSE]
wcb_files <- unlist(wcb_files_all[keep_time], use.names = FALSE)

message("Keeping ", length(time_vals), " time steps with WCB data.")

# Use the first available WCB file to read the exact WCB grid.
first_file <- wcb_files[[1L]]
nc_ref <- nc_open(first_file)
wcb_lon <- get_dim_vals(nc_ref, c("lon", "longitude", "x"))
wcb_lat <- get_dim_vals(nc_ref, c("lat", "latitude", "y"))
nc_close(nc_ref)

if (is.null(wcb_lon) || is.null(wcb_lat)) {
  stop("Could not identify WCB lon/lat dimensions in: ", first_file)
}

# Exact intersection only: keep coordinates present in both datasets.
common_lon <- intersect(wcb_lon, lon_vals)
common_lat <- intersect(wcb_lat, lat_vals)

if (length(common_lon) == 0L || length(common_lat) == 0L) {
  message(
    "Skipping chunk ", TASK_ID,
    ": no overlapping 1-degree grid points with WCB grid.\n",
    "  lon range: [", min(lon_vals), ", ", max(lon_vals), "]\n",
    "  lat range: [", min(lat_vals), ", ", max(lat_vals), "]"
  )
  quit(save = "no", status = 0L)
}

# Indices in the residual chunk.
ix_res <- match(common_lon, lon_vals)
iy_res <- match(common_lat, lat_vals)

# Indices in the WCB grid.
ix_wcb <- match(common_lon, wcb_lon)
iy_wcb <- match(common_lat, wcb_lat)

if (anyNA(ix_res) || anyNA(iy_res)) {
  stop("Internal error: residual intersection indices could not be resolved.")
}
if (anyNA(ix_wcb) || anyNA(iy_wcb)) {
  stop("Internal error: WCB intersection indices could not be resolved.")
}

# Subset residuals to the common lon/lat box and time intersection.
resid_sel <- resid_arr[ix_res, iy_res, , drop = FALSE]

wcb_specs <- list(
  wcb_in_12utc  = "GT800",
  wcb_asc_12utc = "MIDTROP",
  wcb_out_12utc = "LT400"
)

out_nlon <- length(common_lon)
out_nlat <- length(common_lat)
nt       <- length(time_vals)

message("Reading WCB cubes for ", nt, " time steps...")
wcb_cubes <- lapply(names(wcb_specs), function(cube_name) {
  var_name <- wcb_specs[[cube_name]]
  cube <- array(NA_integer_, dim = c(out_nlon, out_nlat, nt))

  for (tt in seq_len(nt)) {
    cube[, , tt] <- read_wcb_snapshot(
      wcb_files[[tt]], var_name, ix_wcb, iy_wcb
    )

    if (tt %% 200 == 0L || tt == nt) {
      message("  ", cube_name, ": ", tt, "/", nt)
    }
  }

  cube
})
names(wcb_cubes) <- names(wcb_specs)

out <- c(
  list(
    metadata = c(
      chunk$metadata %||% list(),
      list(
        task_id      = TASK_ID,
        input_file   = chunk_file,
        wcb_root     = WCB_ROOT,
        lon_count    = out_nlon,
        lat_count    = out_nlat,
        time_count   = nt,
        lon_min      = min(common_lon),
        lon_max      = max(common_lon),
        lat_min      = min(common_lat),
        lat_max      = max(common_lat)
      )
    ),
    lon       = common_lon,
    lat       = common_lat,
    time      = time_vals,
    residuals = resid_sel
  ),
  wcb_cubes
)

saveRDS(out, out_file, compress = "xz")
message("Wrote: ", out_file)