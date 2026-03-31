#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

TASK_ID   <- as.integer(get_arg("--task-id", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
OUT_DIR   <- get_arg("--out-dir", Sys.getenv("OUT_DIR", "./stage2_chunks"))
RESP_RDS  <- get_arg("--response-rds", Sys.getenv("RESP_RDS",
               "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"))
MCC_NC    <- get_arg("--mcc-nc", Sys.getenv("MCC_NC",
               "/scratch/schoelleh96/wp22a/data/medium_cloud_cover.nc"))
HCC_NC    <- get_arg("--hcc-nc", Sys.getenv("HCC_NC",
               "/scratch/schoelleh96/wp22a/data/high_cloud_cover.nc"))
EADY_NC   <- get_arg("--eady-nc", Sys.getenv("EADY_NC",
               "/scratch/schoelleh96/wp22a/data/eady_day_era5_all_years.nc"))
WIND_NC   <- get_arg("--wind-nc", Sys.getenv("WIND_NC",
               "/scratch/schoelleh96/wp22a/data/upper_wind_day_era5_all_years.nc"))
RH_NC     <- get_arg("--rh-nc", Sys.getenv("RH_NC",
               "/scratch/schoelleh96/wp22a/data/relative_humidity_700-850.nc"))
HYDRO_NC  <- get_arg("--hydro-nc", Sys.getenv("HYDRO_NC",
               "/scratch/schoelleh96/wp22a/data/hydrosum.nc"))
GRAD_NC  <- get_arg("--grad-nc", Sys.getenv("GRAD_NC",
               "/scratch/schoelleh96/wp22a/data/z500_grad_mag.nc"))
LAP_NC  <- get_arg("--lap-nc", Sys.getenv("LAP_NC",
               "/scratch/schoelleh96/wp22a/data/z500_laplacian.nc"))

EADY_VAR  <- get_arg("--eady-var", Sys.getenv("EADY_VAR", "eady_growth_rate"))
WIND_VAR  <- get_arg("--wind-var", Sys.getenv("WIND_VAR", "upper_wind_speed"))
MCC_VAR   <- get_arg("--mcc-var", Sys.getenv("MCC_VAR", "medium_cloud_cover"))
HCC_VAR   <- get_arg("--hcc-var", Sys.getenv("HCC_VAR", "high_cloud_cover"))
RH_VAR    <- get_arg("--rh-var", Sys.getenv("RH_VAR", "r"))
HYDRO_VAR <- get_arg("--hydro-var", Sys.getenv("HYDRO_VAR", "hydrosum"))
GRAD_VAR <- get_arg("--grad-var", Sys.getenv("GRAD_VAR", "z_grad_mag"))
LAP_VAR <- get_arg("--lap-var", Sys.getenv("LAP_VAR", "z_laplacian"))

WCB_ROOT  <- get_arg("--wcb-root", Sys.getenv("WCB_ROOT",
               "/scratch/schoelleh96/wp22a/ELIAS_data"))
OVERWRITE <- as.logical(as.integer(get_arg("--overwrite", Sys.getenv("OVERWRITE", "0"))))

chunks_arg <- get_arg("--chunks", Sys.getenv("CHUNKS", "1,1"))
chunks_vec <- as.integer(strsplit(chunks_arg, ",")[[1]])
if (length(chunks_vec) != 2 || any(is.na(chunks_vec)) || any(chunks_vec < 1)) {
  stop("Chunks must be provided as --chunks=LAT_CHUNKS,LON_CHUNKS")
}
LAT_CHUNKS <- chunks_vec[1]
LON_CHUNKS <- chunks_vec[2]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Reading response cube: ", RESP_RDS)
resp <- readRDS(RESP_RDS)

if (!is.list(resp) || is.null(resp$residuals) || is.null(resp$lon) || is.null(resp$lat) || is.null(resp$time)) {
  stop("Response RDS must contain: residuals, lon, lat, time")
}

resid_arr <- resp$residuals   # [lon, lat, time]
lon_vals  <- resp$lon
lat_vals  <- resp$lat
time_vals <- resp$time

if (length(dim(resid_arr)) != 3) {
  stop("Expected residuals array with dims [lon, lat, time]")
}

to_posix_utc <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC") + 12 * 3600)
  y <- as.POSIXct(x, tz = "UTC")
  if (all(is.na(y))) stop("Could not parse response time vector")
  y
}

time_vals <- to_posix_utc(time_vals)

nlon <- length(lon_vals)
nlat <- length(lat_vals)

split_indices <- function(n, k) {
  k <- max(1, min(k, n))
  edges <- floor(seq(0, n, length.out = k + 1))
  out <- vector("list", k)
  for (i in seq_len(k)) out[[i]] <- seq.int(edges[i] + 1, edges[i + 1])
  out
}

lat_chunks <- split_indices(nlat, LAT_CHUNKS)
lon_chunks <- split_indices(nlon, LON_CHUNKS)
n_chunks_total <- length(lat_chunks) * length(lon_chunks)

if (TASK_ID < 1 || TASK_ID > n_chunks_total) {
  stop("TASK_ID out of range: ", TASK_ID, " not in 1:", n_chunks_total)
}

chunk_id <- 1
sel_lat_chunk <- NA_integer_
sel_lon_chunk <- NA_integer_
for (i in seq_along(lat_chunks)) {
  for (j in seq_along(lon_chunks)) {
    if (chunk_id == TASK_ID) {
      sel_lat_chunk <- i
      sel_lon_chunk <- j
      break
    }
    chunk_id <- chunk_id + 1
  }
  if (!is.na(sel_lat_chunk)) break
}

iy <- lat_chunks[[sel_lat_chunk]]
ix <- lon_chunks[[sel_lon_chunk]]

chunk_lon <- lon_vals[ix]
chunk_lat <- lat_vals[iy]

out_file <- file.path(OUT_DIR, sprintf("chunk_%02d.rds", TASK_ID))

find_nc_coord_var <- function(nc, patterns) {
  cand <- names(nc$dim)
  hit <- cand[grepl(paste(patterns, collapse = "|"), cand, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[1])

  candv <- names(nc$var)
  hitv <- candv[grepl(paste(patterns, collapse = "|"), candv, ignore.case = TRUE)]
  if (length(hitv) > 0) return(hitv[1])

  NULL
}

get_nc_coord <- function(nc, patterns) {
  nm <- find_nc_coord_var(nc, patterns)
  if (is.null(nm)) return(NULL)
  if (nm %in% names(nc$dim)) return(nc$dim[[nm]]$vals)
  ncvar_get(nc, nm)
}

parse_nc_time <- function(nc) {
  time_name <- find_nc_coord_var(nc, c("^time$", "time"))
  if (is.null(time_name)) stop("Could not find time variable in NetCDF.")

  vals <- if (time_name %in% names(nc$dim)) nc$dim[[time_name]]$vals else ncvar_get(nc, time_name)
  units_att <- ncatt_get(nc, time_name, "units")$value
  if (is.null(units_att) || is.na(units_att)) stop("Time variable has no units attribute.")

  parts <- strsplit(trimws(units_att), " since ", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("Unsupported NetCDF time units: ", units_att)

  unit <- tolower(trimws(parts[1]))
  origin_txt <- trimws(parts[2])
  origin <- as.POSIXct(origin_txt, tz = "UTC")
  if (is.na(origin)) origin <- as.POSIXct(paste0(origin_txt, " 00:00:00"), tz = "UTC")
  if (is.na(origin)) stop("Could not parse NetCDF time origin: ", origin_txt)

  mult <- switch(unit,
    "seconds" = 1,
    "second"  = 1,
    "secs"    = 1,
    "sec"     = 1,
    "hours"   = 3600,
    "hour"    = 3600,
    "days"    = 86400,
    "day"     = 86400,
    stop("Unsupported NetCDF time unit: ", unit)
  )

  origin + vals * mult
}

convert_lon_to_source <- function(target_lon, source_lon) {
  src_0360 <- all(source_lon >= 0, na.rm = TRUE)
  trg_0360 <- all(target_lon >= 0, na.rm = TRUE)

  if (src_0360 && !trg_0360) return(target_lon %% 360)
  if (!src_0360 && trg_0360) return(ifelse(target_lon > 180, target_lon - 360, target_lon))
  target_lon
}

match_coord_indices <- function(source, target, tol = 1e-4) {
  out <- integer(length(target))
  out[] <- NA_integer_
  for (i in seq_along(target)) {
    hit <- which(abs(source - target[i]) <= tol)
    if (length(hit) > 0) out[i] <- hit[1]
  }
  out
}

match_time_indices <- function(source_time, target_time) {
  idx <- match(target_time, source_time)
  if (all(!is.na(idx))) return(idx)

  out <- integer(length(target_time))
  for (i in seq_along(target_time)) {
    same_day <- as.Date(source_time) == as.Date(target_time[i])
    if (!any(same_day)) {
      out[i] <- NA_integer_
      next
    }
    cand <- which(same_day)
    out[i] <- cand[which.min(abs(as.numeric(difftime(source_time[cand], target_time[i], units = "secs"))))]
  }
  out
}

read_nc_cube_subset <- function(nc_path, var_name, target_lon, target_lat, target_time) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  src_lon  <- get_nc_coord(nc, c("^lon$", "longitude"))
  src_lat  <- get_nc_coord(nc, c("^lat$", "latitude"))
  src_time <- parse_nc_time(nc)

  if (is.null(src_lon) || is.null(src_lat)) stop("Could not detect lon/lat in ", nc_path)

  trg_lon_src <- convert_lon_to_source(target_lon, src_lon)
  lon_idx  <- match_coord_indices(src_lon, trg_lon_src)
  lat_idx  <- match_coord_indices(src_lat, target_lat)
  time_idx <- match_time_indices(src_time, target_time)

  var <- nc$var[[var_name]]
  if (is.null(var)) stop("Variable '", var_name, "' not found in ", nc_path)

  dim_names <- vapply(var$dim, function(d) d$name, character(1))
  dim_lens  <- vapply(var$dim, function(d) d$len, integer(1))

  lon_pos  <- which(grepl("lon|longitude", dim_names, ignore.case = TRUE))[1]
  lat_pos  <- which(grepl("lat|latitude", dim_names, ignore.case = TRUE))[1]
  time_pos <- which(grepl("time", dim_names, ignore.case = TRUE))[1]

  if (is.na(lon_pos))  lon_pos  <- which(dim_lens == length(src_lon))[1]
  if (is.na(lat_pos))  lat_pos  <- which(dim_lens == length(src_lat))[1]
  if (is.na(time_pos)) time_pos <- which(dim_lens == length(src_time))[1]

  good_lon <- which(!is.na(lon_idx))
  good_lat <- which(!is.na(lat_idx))
  good_t   <- which(!is.na(time_idx))

  if (length(good_lon) == 0 || length(good_lat) == 0 || length(good_t) == 0) {
    return(array(NA_real_, dim = c(length(target_lon), length(target_lat), length(target_time))))
  }

  lon_rng  <- range(lon_idx[good_lon])
  lat_rng  <- range(lat_idx[good_lat])
  time_rng <- range(time_idx[good_t])

  start <- rep(1L, length(dim_names))
  count <- dim_lens
  start[lon_pos]  <- lon_rng[1]
  start[lat_pos]  <- lat_rng[1]
  start[time_pos] <- time_rng[1]

  count[lon_pos]  <- diff(lon_rng) + 1L
  count[lat_pos]  <- diff(lat_rng) + 1L
  count[time_pos] <- diff(time_rng) + 1L

  raw <- ncvar_get(nc, var_name, start = start, count = count, collapse_degen = FALSE)

  perm <- c(lon_pos, lat_pos, time_pos)
  raw3 <- aperm(raw, perm)

  lon_local  <- lon_idx - lon_rng[1] + 1L
  lat_local  <- lat_idx - lat_rng[1] + 1L
  time_local <- time_idx - time_rng[1] + 1L

  out <- array(NA_real_, dim = c(length(target_lon), length(target_lat), length(target_time)))

  good_lon_local <- which(!is.na(lon_local))
  good_lat_local <- which(!is.na(lat_local))

  for (tt in seq_along(target_time)) {
    if (is.na(time_local[tt])) next
    if (length(good_lon_local) == 0 || length(good_lat_local) == 0) next

    out[good_lon_local, good_lat_local, tt] <-
      raw3[lon_local[good_lon_local], lat_local[good_lat_local], time_local[tt], drop = FALSE][, , 1]
  }
  out
}

find_wcb_file <- function(root, dt) {
  y <- format(dt, "%Y", tz = "UTC")
  m <- format(dt, "%m", tz = "UTC")
  base <- sprintf("hit_%s_%s", format(dt, "%Y%m%d", tz = "UTC"), format(dt, "%H", tz = "UTC"))
  ddir <- file.path(root, y, m)

  cand <- c(
    file.path(ddir, base),
    file.path(ddir, paste0(base, ".nc"))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) > 0) return(hit[1])

  alt <- Sys.glob(file.path(ddir, paste0(base, "*")))
  if (length(alt) > 0) return(alt[1])

  NA_character_
}

read_wcb_snapshot <- function(nc_path, var_name, target_lon, target_lat) {
  if (is.na(nc_path) || !file.exists(nc_path)) {
    return(matrix(NA_integer_, nrow = length(target_lon), ncol = length(target_lat)))
  }

  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  src_lon <- get_nc_coord(nc, c("^lon$", "longitude"))
  src_lat <- get_nc_coord(nc, c("^lat$", "latitude"))
  if (is.null(src_lon) || is.null(src_lat)) stop("Could not detect lon/lat in ", nc_path)

  trg_lon_src <- convert_lon_to_source(target_lon, src_lon)
  lon_idx <- match_coord_indices(src_lon, trg_lon_src)
  lat_idx <- match_coord_indices(src_lat, target_lat)

  var <- nc$var[[var_name]]
  if (is.null(var)) stop("Variable '", var_name, "' not found in ", nc_path)

  dim_names <- vapply(var$dim, function(d) d$name, character(1))
  dim_lens  <- vapply(var$dim, function(d) d$len, integer(1))

  lon_pos <- which(grepl("lon|longitude", dim_names, ignore.case = TRUE))[1]
  lat_pos <- which(grepl("lat|latitude", dim_names, ignore.case = TRUE))[1]

  if (is.na(lon_pos)) lon_pos <- which(dim_lens == length(src_lon))[1]
  if (is.na(lat_pos)) lat_pos <- which(dim_lens == length(src_lat))[1]

  raw <- ncvar_get(nc, var_name)
  raw2 <- aperm(raw, c(lon_pos, lat_pos, setdiff(seq_along(dim(raw)), c(lon_pos, lat_pos))))
  if (length(dim(raw2)) > 2) raw2 <- raw2[, , 1, drop = FALSE][, , 1]

  out <- matrix(NA_integer_, nrow = length(target_lon), ncol = length(target_lat))

  good_lon <- which(!is.na(lon_idx))
  good_lat <- which(!is.na(lat_idx))

  if (length(good_lon) > 0 && length(good_lat) > 0) {
    out[good_lon, good_lat] <- raw2[lon_idx[good_lon], lat_idx[good_lat], drop = FALSE]
  }

  storage.mode(out) <- "integer"
  out
}

read_or_keep_cube <- function(existing, field_name, nc_path, var_name, target_lon, target_lat, target_time) {
  if (!is.null(existing[[field_name]]) && !OVERWRITE) {
    message("Keeping existing field: ", field_name)
    return(existing[[field_name]])
  }
  message("Reading ", field_name, " subset...")
  read_nc_cube_subset(nc_path, var_name, target_lon, target_lat, target_time)
}

read_or_keep_wcb <- function(existing, field_name, files, var_name, target_lon, target_lat, nt) {
  if (!is.null(existing[[field_name]]) && !OVERWRITE) {
    message("Keeping existing field: ", field_name)
    return(existing[[field_name]])
  }

  message("Reading WCB field: ", field_name)
  out <- array(NA_integer_, dim = c(length(target_lon), length(target_lat), nt))
  for (tt in seq_len(nt)) {
    out[, , tt] <- read_wcb_snapshot(files[[tt]], var_name, target_lon, target_lat)
    if (tt %% 200 == 0 || tt == nt) {
      message("  ", field_name, ": ", tt, "/", nt)
    }
  }
  out
}

existing <- NULL
if (file.exists(out_file) && !OVERWRITE) {
  message("Existing chunk found, loading for incremental update: ", out_file)
  existing <- readRDS(out_file)
}

if (file.exists(out_file) && OVERWRITE) {
  message("Existing chunk found, but OVERWRITE=TRUE, rebuilding: ", out_file)
}

resid_chunk <- if (!is.null(existing$residuals) && !OVERWRITE) {
  existing$residuals
} else {
  resid_arr[ix, iy, , drop = FALSE]
}

mcc_chunk <- read_or_keep_cube(existing %||% list(), "mcc", MCC_NC, MCC_VAR, chunk_lon, chunk_lat, time_vals)
hcc_chunk <- read_or_keep_cube(existing %||% list(), "hcc", HCC_NC, HCC_VAR, chunk_lon, chunk_lat, time_vals)
eady_chunk <- read_or_keep_cube(existing %||% list(), "eady", EADY_NC, EADY_VAR, chunk_lon, chunk_lat, time_vals)
upper_wind_chunk <- read_or_keep_cube(existing %||% list(), "upper_wind", WIND_NC, WIND_VAR, chunk_lon, chunk_lat, time_vals)

rh_chunk <- if (!is.null(RH_NC) && nzchar(RH_NC) && file.exists(RH_NC)) {
  read_or_keep_cube(existing %||% list(), "rh_700_850", RH_NC, RH_VAR, chunk_lon, chunk_lat, time_vals)
} else {
  existing$rh_700_850 %||% NULL
}

hydro_chunk <- if (!is.null(HYDRO_NC) && nzchar(HYDRO_NC) && file.exists(HYDRO_NC)) {
  read_or_keep_cube(existing %||% list(), "hydrosum", HYDRO_NC, HYDRO_VAR, chunk_lon, chunk_lat, time_vals)
} else {
  existing$hydrosum %||% NULL
}

grad_chunk <- if (!is.null(GRAD_NC) && nzchar(GRAD_NC) && file.exists(GRAD_NC)) {
  read_or_keep_cube(existing %||% list(), "z_grad_mag", GRAD_NC, GRAD_VAR, 
  chunk_lon, chunk_lat, time_vals)
} else {
  existing$z_grad_mag %||% NULL
}
lap_chunk <- if (!is.null(LAP_NC) && nzchar(LAP_NC) && file.exists(LAP_NC)) {
  read_or_keep_cube(existing %||% list(), "z_laplacian", LAP_NC, LAP_VAR, 
  chunk_lon, chunk_lat, time_vals)
} else {
  existing$z_laplacian %||% NULL
}

nt <- length(time_vals)
nx <- length(chunk_lon)
ny <- length(chunk_lat)

need_any_wcb <- OVERWRITE || any(vapply(
  c("wcb_in_12utc", "wcb_asc_12utc", "wcb_out_12utc",
    "wcb_in_lag12", "wcb_asc_lag12", "wcb_out_lag12",
    "wcb_in_lag24", "wcb_asc_lag24", "wcb_out_lag24"),
  function(v) is.null(existing[[v]]),
  logical(1)
))

if (need_any_wcb) {
  message("Preparing WCB file lookup...")
  f0  <- vector("list", nt)
  f12 <- vector("list", nt)
  f24 <- vector("list", nt)

  for (tt in seq_along(time_vals)) {
    t0  <- time_vals[tt]
    t12 <- t0 - 12 * 3600
    t24 <- t0 - 24 * 3600

    f0[[tt]]  <- find_wcb_file(WCB_ROOT, t0)
    f12[[tt]] <- find_wcb_file(WCB_ROOT, t12)
    f24[[tt]] <- find_wcb_file(WCB_ROOT, t24)
  }
} else {
  f0 <- f12 <- f24 <- NULL
}

wcb_in_12utc  <- read_or_keep_wcb(existing %||% list(), "wcb_in_12utc",  f0,  "GT800",   chunk_lon, chunk_lat, nt)
wcb_asc_12utc <- read_or_keep_wcb(existing %||% list(), "wcb_asc_12utc", f0,  "MIDTROP", chunk_lon, chunk_lat, nt)
wcb_out_12utc <- read_or_keep_wcb(existing %||% list(), "wcb_out_12utc", f0,  "LT400",   chunk_lon, chunk_lat, nt)

wcb_in_lag12  <- read_or_keep_wcb(existing %||% list(), "wcb_in_lag12",  f12, "GT800",   chunk_lon, chunk_lat, nt)
wcb_asc_lag12 <- read_or_keep_wcb(existing %||% list(), "wcb_asc_lag12", f12, "MIDTROP", chunk_lon, chunk_lat, nt)
wcb_out_lag12 <- read_or_keep_wcb(existing %||% list(), "wcb_out_lag12", f12, "LT400",   chunk_lon, chunk_lat, nt)

wcb_in_lag24  <- read_or_keep_wcb(existing %||% list(), "wcb_in_lag24",  f24, "GT800",   chunk_lon, chunk_lat, nt)
wcb_asc_lag24 <- read_or_keep_wcb(existing %||% list(), "wcb_asc_lag24", f24, "MIDTROP", chunk_lon, chunk_lat, nt)
wcb_out_lag24 <- read_or_keep_wcb(existing %||% list(), "wcb_out_lag24", f24, "LT400",   chunk_lon, chunk_lat, nt)

out <- list(
  metadata = list(
    task_id = TASK_ID,
    lat_chunks = LAT_CHUNKS,
    lon_chunks = LON_CHUNKS,
    lat_chunk_id = sel_lat_chunk,
    lon_chunk_id = sel_lon_chunk,
    global_lon_index = ix,
    global_lat_index = iy,
    weight_mode = resp$weight_mode %||% NA_character_,
    incremental_update = !OVERWRITE && file.exists(out_file)
  ),
  lon = chunk_lon,
  lat = chunk_lat,
  time = time_vals,
  residuals = resid_chunk,
  mcc = mcc_chunk,
  hcc = hcc_chunk,
  eady = eady_chunk,
  upper_wind = upper_wind_chunk,
  wcb_in_12utc = wcb_in_12utc,
  wcb_asc_12utc = wcb_asc_12utc,
  wcb_out_12utc = wcb_out_12utc,
  wcb_in_lag12 = wcb_in_lag12,
  wcb_asc_lag12 = wcb_asc_lag12,
  wcb_out_lag12 = wcb_out_lag12,
  wcb_in_lag24 = wcb_in_lag24,
  wcb_asc_lag24 = wcb_asc_lag24,
  wcb_out_lag24 = wcb_out_lag24
)

if (!is.null(rh_chunk)) out$rh_700_850 <- rh_chunk
if (!is.null(hydro_chunk)) out$hydrosum <- hydro_chunk
if (!is.null(grad_chunk)) out$grad <- grad_chunk
if (!is.null(lap_chunk)) out$lap <- lap_chunk

saveRDS(out, out_file, compress = "xz")
message("Wrote: ", out_file)