#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
  library(parallel)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

OUT_DIR  <- get_arg(
  "--out-dir", Sys.getenv("OUT_DIR", "./stage2_chunks"))
RESP_RDS <- get_arg("--response-rds", Sys.getenv(
  "RESP_RDS",
  "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
))
RH_NC    <- get_arg("--rh-nc", Sys.getenv(
  "RH_NC",
  "/scratch/schoelleh96/wp22a/data/relative_humidity_500-850.nc"
))
GRAD_NC  <- get_arg("--grad-nc", Sys.getenv(
  "GRAD_NC",
  "/scratch/schoelleh96/wp22a/data/z500_grad_mag.nc"
))
LAP_NC   <- get_arg("--lap-nc", Sys.getenv(
  "LAP_NC",
  "/scratch/schoelleh96/wp22a/data/z500_laplacian.nc"
))
RH_VAR   <- get_arg("--rh-var",   Sys.getenv("RH_VAR",   "r"))
GRAD_VAR <- get_arg("--grad-var", Sys.getenv("GRAD_VAR", "z_grad_mag"))
LAP_VAR  <- get_arg("--lap-var",  Sys.getenv("LAP_VAR",  "z_laplacian"))
OVERWRITE <- as.logical(as.integer(
  get_arg("--overwrite", Sys.getenv("OVERWRITE", "0"))
))
N_WORKERS <- as.integer(get_arg(
  "--workers",
  Sys.getenv("N_WORKERS", as.character(detectCores(logical = FALSE)))
))

chunks_arg <- get_arg("--chunks", Sys.getenv("CHUNKS", "1,1"))
chunks_vec <- as.integer(strsplit(chunks_arg, ",")[[1]])
if (length(chunks_vec) != 2 ||
    any(is.na(chunks_vec)) ||
    any(chunks_vec < 1)) {
  stop("--chunks must be LAT_CHUNKS,LON_CHUNKS")
}
LAT_CHUNKS <- chunks_vec[1]
LON_CHUNKS <- chunks_vec[2]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# NetCDF helpers  (identical logic to template script)
# ---------------------------------------------------------------------------

find_nc_coord_var <- function(nc, patterns) {
  combined  <- paste(patterns, collapse = "|")
  hit_dims  <- names(nc$dim)[
    grepl(combined, names(nc$dim), ignore.case = TRUE)]
  if (length(hit_dims) > 0) return(hit_dims[1])
  hit_vars  <- names(nc$var)[
    grepl(combined, names(nc$var), ignore.case = TRUE)]
  if (length(hit_vars) > 0) return(hit_vars[1])
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
  if (is.null(time_name)) stop("Could not find time variable.")

  vals <- if (time_name %in% names(nc$dim)) {
    nc$dim[[time_name]]$vals
  } else {
    ncvar_get(nc, time_name)
  }
  units_att <- ncatt_get(nc, time_name, "units")$value
  if (is.null(units_att) || is.na(units_att)) {
    stop("Time variable has no units attribute.")
  }

  parts <- strsplit(trimws(units_att), " since ", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Unsupported NetCDF time units: ", units_att)
  }
  unit       <- tolower(trimws(parts[1]))
  origin_txt <- trimws(parts[2])
  origin     <- as.POSIXct(origin_txt, tz = "UTC")
  if (is.na(origin)) {
    origin <- as.POSIXct(paste0(origin_txt, " 00:00:00"), tz = "UTC")
  }
  if (is.na(origin)) {
    stop("Could not parse NC time origin: ", origin_txt)
  }

  mult <- switch(unit,
    seconds = 1, second = 1, secs = 1, sec = 1,
    hours   = 3600, hour = 3600,
    days    = 86400, day = 86400,
    stop("Unsupported NC time unit: ", unit)
  )
  origin + vals * mult
}

convert_lon_to_source <- function(target_lon, source_lon) {
  src_0360 <- all(source_lon >= 0, na.rm = TRUE)
  trg_0360 <- all(target_lon >= 0, na.rm = TRUE)
  if (src_0360 && !trg_0360)  return(target_lon %% 360)
  if (!src_0360 && trg_0360) {
    return(ifelse(target_lon > 180, target_lon - 360, target_lon))
  }
  target_lon
}

match_coord_indices <- function(source, target, tol = 1e-4) {
  out    <- integer(length(target))
  out[]  <- NA_integer_
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
    if (!any(same_day)) { out[i] <- NA_integer_; next }
    cand  <- which(same_day)
    diffs <- abs(as.numeric(
      difftime(source_time[cand], target_time[i], units = "secs")
    ))
    out[i] <- cand[which.min(diffs)]
  }
  out
}

# ---------------------------------------------------------------------------
# Full-cube reader and subsetter
# ---------------------------------------------------------------------------

# Reads the entire variable into a list(lon, lat, time, data)
# with data permuted to [lon, lat, time].  Called once per NC file.
read_nc_full_cube <- function(nc_path, var_name) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  src_lon  <- get_nc_coord(nc, c("^lon$", "longitude"))
  src_lat  <- get_nc_coord(nc, c("^lat$", "latitude"))
  src_time <- parse_nc_time(nc)

  if (is.null(src_lon) || is.null(src_lat)) {
    stop("Could not detect lon/lat in ", nc_path)
  }

  var <- nc$var[[var_name]]
  if (is.null(var)) {
    stop("Variable '", var_name, "' not found in ", nc_path)
  }

  dim_names <- vapply(var$dim, function(d) d$name, character(1))
  dim_lens  <- vapply(var$dim, function(d) d$len,  integer(1))

  lon_pos  <- which(
    grepl("lon|longitude", dim_names, ignore.case = TRUE))[1]
  lat_pos  <- which(
    grepl("lat|latitude",  dim_names, ignore.case = TRUE))[1]
  time_pos <- which(
    grepl("time",          dim_names, ignore.case = TRUE))[1]

  if (is.na(lon_pos))  lon_pos  <- which(dim_lens == length(src_lon))[1]
  if (is.na(lat_pos))  lat_pos  <- which(dim_lens == length(src_lat))[1]
  if (is.na(time_pos)) time_pos <- which(dim_lens == length(src_time))[1]

  raw  <- ncvar_get(nc, var_name, collapse_degen = FALSE)
  raw3 <- aperm(raw, c(lon_pos, lat_pos, time_pos))

  list(lon = src_lon, lat = src_lat, time = src_time, data = raw3)
}

# Extracts a [lon, lat, time] sub-array from a full cube using
# vectorised R indexing — no per-time-step loop needed.
subset_nc_cube <- function(
    full_cube, target_lon, target_lat, target_time, tol = 1e-4) {

  trg_lon_src <- convert_lon_to_source(target_lon, full_cube$lon)
  lon_idx     <- match_coord_indices(full_cube$lon, trg_lon_src, tol)
  lat_idx     <- match_coord_indices(full_cube$lat, target_lat,  tol)
  time_idx    <- match_time_indices(full_cube$time, target_time)

  out <- array(
    NA_real_,
    dim = c(
      length(target_lon), length(target_lat), length(target_time)
    )
  )

  good_lon <- which(!is.na(lon_idx))
  good_lat <- which(!is.na(lat_idx))
  good_t   <- which(!is.na(time_idx))

  if (length(good_lon) == 0 ||
      length(good_lat) == 0 ||
      length(good_t)   == 0) {
    return(out)
  }

  # R array indexing here is a cross-product, which is exactly what we
  # want: extract the [good_lon x good_lat x good_t] sub-cube in one shot
  out[good_lon, good_lat, good_t] <- full_cube$data[
    lon_idx[good_lon], lat_idx[good_lat], time_idx[good_t]
  ]
  out
}

# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------

split_indices <- function(n, k) {
  k     <- max(1, min(k, n))
  edges <- floor(seq(0, n, length.out = k + 1))
  lapply(seq_len(k), function(i) seq.int(edges[i] + 1, edges[i + 1]))
}

to_posix_utc <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC") + 12 * 3600)
  y <- as.POSIXct(x, tz = "UTC")
  if (all(is.na(y))) stop("Could not parse response time vector")
  y
}

# ---------------------------------------------------------------------------
# Load response cube
# ---------------------------------------------------------------------------

message("Reading response cube: ", RESP_RDS)
resp <- readRDS(RESP_RDS)

required_fields <- c("residuals", "lon", "lat", "time")
if (!is.list(resp) || !all(required_fields %in% names(resp))) {
  stop("Response RDS must contain: ",
       paste(required_fields, collapse = ", "))
}

resid_arr <- resp$residuals   # [lon, lat, time]
lon_vals  <- resp$lon
lat_vals  <- resp$lat
time_vals <- to_posix_utc(resp$time)

if (length(dim(resid_arr)) != 3) {
  stop("Expected residuals array with dims [lon, lat, time]")
}

start_cutoff <- as.POSIXct("1950-01-01 00:00:00", tz = "UTC")
keep_time <- time_vals >= start_cutoff
time_vals <- time_vals[keep_time]
resid_arr <- resid_arr[, , keep_time, drop = FALSE]

# ---------------------------------------------------------------------------
# Read all NC files into memory once
# ---------------------------------------------------------------------------

nc_source_defs <- list(
  rh_500_850  = list(path = RH_NC,   var = RH_VAR),
  z_grad_mag  = list(path = GRAD_NC, var = GRAD_VAR),
  z_laplacian = list(path = LAP_NC,  var = LAP_VAR)
)

full_cubes <- list()
for (field in names(nc_source_defs)) {
  src <- nc_source_defs[[field]]
  if (is.null(src$path) ||
      !nzchar(src$path %||% "") ||
      !file.exists(src$path)) next
  message("Loading full cube: ", field, "  [", src$path, "]")
  full_cubes[[field]] <- read_nc_full_cube(src$path, src$var)
}

# ---------------------------------------------------------------------------
# Build chunk specification table
# ---------------------------------------------------------------------------

lat_idx_chunks <- split_indices(length(lat_vals), LAT_CHUNKS)
lon_idx_chunks <- split_indices(length(lon_vals), LON_CHUNKS)

chunk_specs <- list()
task_id     <- 1L
for (i in seq_along(lat_idx_chunks)) {
  for (j in seq_along(lon_idx_chunks)) {
    chunk_specs[[task_id]] <- list(
      task_id      = task_id,
      lat_chunk_id = i,
      lon_chunk_id = j,
      iy           = lat_idx_chunks[[i]],
      ix           = lon_idx_chunks[[j]]
    )
    task_id <- task_id + 1L
  }
}

# ---------------------------------------------------------------------------
# Per-chunk worker  (forked by mclapply — inherits all in-memory data
# via copy-on-write, so NC files are never re-read from disk)
# ---------------------------------------------------------------------------

process_chunk <- function(spec) {
  ix        <- spec$ix
  iy        <- spec$iy
  chunk_lon <- lon_vals[ix]
  chunk_lat <- lat_vals[iy]
  out_path  <- file.path(
    OUT_DIR, sprintf("chunk_%02d.rds", spec$task_id)
  )

  existing <- if (file.exists(out_path) && !OVERWRITE) {
    readRDS(out_path)
  } else {
    NULL
  }

  resid_chunk <- if (!is.null(existing$residuals) && !OVERWRITE) {
    existing$residuals
  } else {
    resid_arr[ix, iy, , drop = FALSE]
  }

  field_chunks <- lapply(names(full_cubes), function(field) {
    if (!is.null(existing[[field]]) && !OVERWRITE) {
      return(existing[[field]])
    }
    subset_nc_cube(
      full_cubes[[field]], chunk_lon, chunk_lat, time_vals
    )
  })
  names(field_chunks) <- names(full_cubes)

  out <- c(
    list(
      metadata = list(
        task_id            = spec$task_id,
        lat_chunks         = LAT_CHUNKS,
        lon_chunks         = LON_CHUNKS,
        lat_chunk_id       = spec$lat_chunk_id,
        lon_chunk_id       = spec$lon_chunk_id,
        global_lon_index   = ix,
        global_lat_index   = iy,
        weight_mode        = resp$weight_mode %||% NA_character_,
        incremental_update = !OVERWRITE && file.exists(out_path)
      ),
      lon       = chunk_lon,
      lat       = chunk_lat,
      time      = time_vals,
      residuals = resid_chunk
    ),
    field_chunks
  )

  saveRDS(out, out_path, compress = "xz")
  message(
    sprintf("Chunk %02d/%02d -> %s",
            spec$task_id, length(chunk_specs), out_path))
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

n_chunks <- length(chunk_specs)
message(sprintf(
  "Dispatching %d chunk(s) across %d worker(s)...",
  n_chunks, min(N_WORKERS, n_chunks)
))

results <- mclapply(
  chunk_specs,
  process_chunk,
  mc.cores       = N_WORKERS,
  mc.preschedule = FALSE   # each chunk may differ in size
)

failed <- vapply(results, inherits, logical(1), "error")
if (any(failed)) {
  first_err <- conditionMessage(results[[which(failed)[1]]])
  stop(sum(failed), " chunk(s) failed. First error: ", first_err)
}

message(sprintf("Done — %d chunk(s) written to %s.", n_chunks, OUT_DIR))