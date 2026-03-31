#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

RAW_ROOT  <- get_arg("--raw-root", "/scratch/schoelleh96/wp22a/ELIAS_data")
WR_RDS    <- get_arg("--wr-rds", "/home/schoelleh96/wp22a/data/wrnames.rds")
MASKS_RDS <- get_arg("--masks-rds", "/home/schoelleh96/wp22a/data/wcb_masks_df.rds")
OUT_FILE  <- get_arg("--out-file", "/scratch/schoelleh96/wp22a/data/wcb_predictors_all_years.rds")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

wr_df <- as.data.table(readRDS(WR_RDS))
wr_df[, date := as.Date(date)]

mask_df <- as.data.table(readRDS(MASKS_RDS))
stopifnot(all(c("lon", "lat", "variable", "regime", "lag_name", "mask_value") %in% names(mask_df)))

find_wcb_file <- function(root, date, hour) {
  y <- format(date, "%Y")
  m <- format(date, "%m")
  base <- sprintf("hit_%s_%02d", format(date, "%Y%m%d"), hour)
  f <- file.path(root, y, m, base)
  if (file.exists(f)) return(f)
  if (file.exists(paste0(f, ".nc"))) return(paste0(f, ".nc"))
  NA_character_
}

read_wcb_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  list(
    GT800   = ncvar_get(nc, "GT800"),
    MIDTROP = ncvar_get(nc, "MIDTROP"),
    LT400   = ncvar_get(nc, "LT400"),
    lon     = ncvar_get(nc, "lon"),
    lat     = ncvar_get(nc, "lat")
  )
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

make_mask_lookup <- function(mask_df, raw_lon, raw_lat) {
  out <- list()

  for (reg in unique(mask_df$regime)) {
    out[[reg]] <- list()
    for (var in unique(mask_df$variable)) {
      out[[reg]][[var]] <- list()
      for (lag in unique(mask_df$lag_name)) {
        sub <- mask_df[regime == reg & variable == var & lag_name == lag & mask_value > 0]
        if (nrow(sub) == 0) next

        lon_idx <- match_coord_indices(raw_lon, convert_lon_to_source(sub$lon, raw_lon))
        lat_idx <- match_coord_indices(raw_lat, sub$lat)
        ok <- !is.na(lon_idx) & !is.na(lat_idx)
        if (!any(ok)) next

        # cosine weighting normalized to 45N = 1
        w <- cos(sub$lat[ok] * pi / 180) / cos(45 * pi / 180)

        out[[reg]][[var]][[lag]] <- list(
          i = lon_idx[ok],
          j = lat_idx[ok],
          w = w
        )
      }
    }
  }

  out
}

weighted_sum_45 <- function(field, idx) {
  if (is.null(idx) || is.null(field)) return(NA_real_)
  v <- field[cbind(idx$i, idx$j)]
  ok <- !is.na(v) & !is.na(idx$w)
  if (!any(ok)) return(NA_real_)
  sum(v[ok] * idx$w[ok])
}

# --- discover all available WCB files and dates from RAW_ROOT ---
all_files <- list.files(
  RAW_ROOT,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "^hit_[0-9]{8}_(00|12)(\\.nc)?$"
)

if (length(all_files) == 0) {
  stop("No WCB files found under RAW_ROOT: ", RAW_ROOT)
}

parse_file_date <- function(path) {
  bn <- basename(path)
  dstr <- sub("^hit_([0-9]{8})_[0-9]{2}(\\.nc)?$", "\\1", bn)
  as.Date(dstr, format = "%Y%m%d")
}

file_dates <- parse_file_date(all_files)
file_dates <- file_dates[!is.na(file_dates)]
if (length(file_dates) == 0) stop("Could not parse any dates from WCB file names.")

first_wcb_date <- min(file_dates)
last_wcb_date  <- max(file_dates)

# use the first available WCB file to infer raw lon/lat grid
probe_file <- all_files[order(all_files)][1]
probe_nc <- nc_open(probe_file)
raw_lon <- ncvar_get(probe_nc, "lon")
raw_lat <- ncvar_get(probe_nc, "lat")
nc_close(probe_nc)

mask_lookup <- make_mask_lookup(mask_df, raw_lon, raw_lat)

# daily date sequence over the available WCB archive period
dates <- seq.Date(first_wcb_date, last_wcb_date, by = "day")

# wrnames merged onto the WCB period
dt <- data.table(date = dates)
dt <- merge(dt, wr_df[, .(date, wrname)], by = "date", all.x = TRUE, sort = TRUE)

out <- data.table(
  date = dt$date,
  wrname = dt$wrname,
  in_24 = NA_real_,
  asc_12 = NA_real_,
  out_00 = NA_real_
)

for (k in seq_len(nrow(dt))) {
  reg <- dt$wrname[k]
  if (is.na(reg) || !(reg %in% names(mask_lookup))) next

  # predictor timing:
  # inflow  -> lag24 -> previous day 12 UTC
  # ascent  -> lag12 -> previous day 00 UTC
  # outflow -> lag00 -> current day 12 UTC
  date_in  <- dt$date[k] - 1
  date_asc <- dt$date[k] - 1
  date_out <- dt$date[k]

  f12_in  <- read_wcb_file(find_wcb_file(RAW_ROOT, date_in, 12))
  f00_asc <- read_wcb_file(find_wcb_file(RAW_ROOT, date_asc, 0))
  f12_out <- read_wcb_file(find_wcb_file(RAW_ROOT, date_out, 12))

  out$in_24[k]  <- weighted_sum_45(f12_in$GT800,     mask_lookup[[reg]][["inflow"]][["lag24"]])
  out$asc_12[k] <- weighted_sum_45(f00_asc$MIDTROP,  mask_lookup[[reg]][["ascent"]][["lag12"]])
  out$out_00[k] <- weighted_sum_45(f12_out$LT400,    mask_lookup[[reg]][["outflow"]][["lag0"]])
}

saveRDS(out, OUT_FILE, compress = "xz")
message("Saved combined WCB predictors to: ", OUT_FILE)