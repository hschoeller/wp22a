#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
})

# -----------------------------
# simple CLI parser: --key value
# -----------------------------
parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[i]
    if (!startsWith(key, "--")) stop("Expected --key, got: ", key)
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[i + 1L], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[i + 1L]
      i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

need <- function(name) {
  if (is.null(args[[name]])) stop("Missing required argument --", name)
  args[[name]]
}

YEAR    <- as.integer(need("year"))
WR_RDS  <- need("wr-rds")
HIT_ROOT <- need("hit-root")
OUT_DIR <- need("out-dir")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FILE <- file.path(OUT_DIR, sprintf("partial_%04d.rds", YEAR))

var_map <- c(
  inflow  = "GT800",
  ascent = "MIDTROP",
  outflow = "LT400"
)

lag_names <- c("lag0", "lag12", "lag24")
lag_hours <- c(0L, 12L, 24L)

# -----------------------------
# helpers
# -----------------------------
find_hit_file <- function(root, dt) {
  y <- format(dt, "%Y", tz = "UTC")
  m <- format(dt, "%m", tz = "UTC")
  base <- sprintf("hit_%s_%s",
                  format(dt, "%Y%m%d", tz = "UTC"),
                  format(dt, "%H", tz = "UTC"))
  ddir <- file.path(root, y, m)

  cand <- c(
    file.path(ddir, base),
    file.path(ddir, paste0(base, ".nc"))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) > 0L) return(hit[1])

  alt <- Sys.glob(file.path(ddir, paste0(base, "*")))
  if (length(alt) > 0L) return(alt[1])

  NA_character_
}

make_targets <- function(d) {
  d <- as.Date(d)
  list(
    lag0  = as.POSIXct(sprintf("%s 12:00:00", format(d, "%F")), tz = "UTC"),
    lag12 = as.POSIXct(sprintf("%s 00:00:00", format(d, "%F")), tz = "UTC"),
    lag24 = as.POSIXct(sprintf("%s 12:00:00", format(d - 1, "%F")), tz = "UTC")
  )
}

pick_dim_name <- function(avail, candidates, what) {
  hit <- candidates[candidates %in% avail]
  if (length(hit) == 0L) {
    stop("Could not find ", what, " dimension. Available dims: ",
         paste(avail, collapse = ", "))
  }
  hit[1]
}

read_2d_var_lonlat <- function(nc, var_name, lon_name, lat_name) {
  if (!var_name %in% names(nc$var)) {
    stop("Variable ", var_name, " not found in file")
  }

  v <- nc$var[[var_name]]
  x <- ncvar_get(nc, var_name)

  if (length(dim(x)) != 2L) {
    stop("Variable ", var_name, " is not 2D")
  }

  dim_names <- vapply(v$dim, function(dd) dd$name, character(1))
  lon_pos <- match(lon_name, dim_names)
  lat_pos <- match(lat_name, dim_names)

  if (is.na(lon_pos) || is.na(lat_pos)) {
    stop("Could not identify lon/lat dimension order for variable ", var_name)
  }

  if (!identical(c(lon_pos, lat_pos), c(1L, 2L))) {
    x <- aperm(x, c(lon_pos, lat_pos))
  }

  storage.mode(x) <- "double"
  x
}

read_snapshot <- function(nc_path, var_map) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)

  dim_names <- names(nc$dim)
  lon_name <- pick_dim_name(dim_names, c("lon", "longitude", "x"), "lon")
  lat_name <- pick_dim_name(dim_names, c("lat", "latitude", "y"), "lat")

  lon <- ncvar_get(nc, lon_name)
  lat <- ncvar_get(nc, lat_name)

  dat <- lapply(unname(var_map), function(vn) {
    read_2d_var_lonlat(nc, vn, lon_name, lat_name)
  })
  names(dat) <- names(var_map)

  list(
    lon = lon,
    lat = lat,
    data = dat
  )
}

message("Reading WR file: ", WR_RDS)
wr <- readRDS(WR_RDS)

if (!is.data.frame(wr)) stop("WR file must contain a data.frame")
if (!all(c("date", "wrname") %in% names(wr))) {
  stop("WR data.frame must contain columns: date, wrname")
}

wr$date <- as.Date(wr$date)

if (is.factor(wr$wrname)) {
  regime_levels <- levels(wr$wrname)
  wr$wrname <- as.character(wr$wrname)
  regimes <- regime_levels[regime_levels %in% wr$wrname]
} else {
  wr$wrname <- as.character(wr$wrname)
  regimes <- sort(unique(wr$wrname[!is.na(wr$wrname) & nzchar(wr$wrname)]))
}

wr <- wr[!is.na(wr$date) & !is.na(wr$wrname) & nzchar(wr$wrname), c("date", "wrname")]
wr_year <- wr[format(wr$date, "%Y") == sprintf("%04d", YEAR), , drop = FALSE]

if (nrow(wr_year) == 0L) {
  stop("No WR entries found for year ", YEAR)
}

message("WR rows for ", YEAR, ": ", nrow(wr_year))

# -----------------------------
# find one sample file to get grid
# -----------------------------
sample_file <- NA_character_
for (ii in seq_len(nrow(wr_year))) {
  targets <- make_targets(wr_year$date[ii])
  for (nm in lag_names) {
    ff <- find_hit_file(HIT_ROOT, targets[[nm]])
    if (!is.na(ff)) {
      sample_file <- ff
      break
    }
  }
  if (!is.na(sample_file)) break
}

if (is.na(sample_file)) {
  stop("No hit files found for year ", YEAR)
}

message("Using sample file for grid: ", sample_file)
sample <- read_snapshot(sample_file, var_map)

lon <- sample$lon
lat <- sample$lat

# -----------------------------
# Apply spatial cropping
# -----------------------------
LON_BOUND <- c(-80, 40)
LAT_BOUND <- c(30, 90)

# handle possible 0–360 longitude
lon_full <- sample$lon
lat_full <- sample$lat

# -----------------------------
# Apply spatial cropping
# -----------------------------
LON_BOUND <- c(-80, 40)
LAT_BOUND <- c(30, 90)

# handle possible 0–360 longitude
lon_adj <- ifelse(lon_full > 180, lon_full - 360, lon_full)

lon_idx <- which(lon_adj >= LON_BOUND[1] & lon_adj <= LON_BOUND[2])
lat_idx <- which(lat_full >= LAT_BOUND[1] & lat_full <= LAT_BOUND[2])

if (length(lon_idx) == 0 || length(lat_idx) == 0) {
  stop("Cropping bounds produced empty domain")
}

# keep both full and cropped coords
lon <- lon_full[lon_idx]
lat <- lat_full[lat_idx]

nlon_full <- length(lon_full)
nlat_full <- length(lat_full)

nlon <- length(lon)
nlat <- length(lat)
nlag <- length(lag_names)
nreg <- length(regimes)
nvar <- length(var_map)

message("Cropped domain: ",
        nlon, " x ", nlat,
        " (", round(100 * nlon * nlat / (nlon_full * nlat_full), 1), "% of grid)")
arr_dim <- c(nlon, nlat, nlag, nreg, nvar)
arr_dn  <- list(NULL, NULL, lag_names, regimes, names(var_map))

hit_counts   <- array(0,  dim = arr_dim, dimnames = arr_dn)
valid_counts <- array(0L, dim = arr_dim, dimnames = arr_dn)

event_counts <- matrix(
  0L, nrow = nlag, ncol = nreg,
  dimnames = list(lag = lag_names, regime = regimes)
)

file_counts <- matrix(
  0L, nrow = nlag, ncol = nreg,
  dimnames = list(lag = lag_names, regime = regimes)
)

# tiny 1-file cache:
# if we read lag0 of day d as the last file,
# lag24 of day d+1 is the same file
last_path <- NULL
last_snap <- NULL

get_snapshot_cached <- function(path) {
  if (!is.null(last_path) && identical(path, last_path)) {
    return(last_snap)
  }
  snap <- read_snapshot(path, var_map)
  last_path <<- path
  last_snap <<- snap
  snap
}

# read order chosen so that lag0 is read last;
# then next day's lag24 can reuse the cache
read_order <- c(3L, 2L, 1L)  # lag24, lag12, lag0

message("Processing year ", YEAR, " ...")

for (ii in seq_len(nrow(wr_year))) {
  d <- wr_year$date[ii]
  reg_name <- wr_year$wrname[ii]
  reg_i <- match(reg_name, regimes)

  if (is.na(reg_i)) next

  event_counts[, reg_i] <- event_counts[, reg_i] + 1L
  targets <- make_targets(d)

  for (lag_i in read_order) {
    lag_name <- lag_names[lag_i]
    f <- find_hit_file(HIT_ROOT, targets[[lag_name]])
    if (is.na(f)) next

    snap <- get_snapshot_cached(f)

if (!identical(as.numeric(snap$lon), as.numeric(lon_full)) ||
    !identical(as.numeric(snap$lat), as.numeric(lat_full))) {
  stop("Grid coordinate mismatch in file: ", f)
}

    file_counts[lag_i, reg_i] <- file_counts[lag_i, reg_i] + 1L

    for (v in seq_len(nvar)) {
      x <- snap$data[[v]][lon_idx, lat_idx, drop = FALSE]
      ok <- !is.na(x)
      x[!ok] <- 0

      hit_counts[, , lag_i, reg_i, v] <-
        hit_counts[, , lag_i, reg_i, v] + x

      valid_counts[, , lag_i, reg_i, v] <-
        valid_counts[, , lag_i, reg_i, v] + ok
    }
  }

  if (ii %% 50L == 0L || ii == nrow(wr_year)) {
    message(sprintf("[%04d] %d / %d dates done", YEAR, ii, nrow(wr_year)))
  }
}

res <- list(
  year = YEAR,
  lon = lon,
  lat = lat,
  lag_names = lag_names,
  lag_hours = lag_hours,
  regimes = regimes,
  variables = names(var_map),
  nc_variables = unname(var_map),
  hit_counts = hit_counts,
  valid_counts = valid_counts,
  event_counts = event_counts,
  file_counts = file_counts
)

message("Writing partial result: ", OUT_FILE)
saveRDS(res, OUT_FILE, compress = FALSE)
message("Done.")