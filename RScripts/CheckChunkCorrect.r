# # #!/usr/bin/env Rscript

# # suppressPackageStartupMessages({
# #   library(ncdf4)
# # })

# # resid_file <- "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
# # rh_file    <- "/scratch/schoelleh96/wp22a/data/relative_humidity_700-850.nc"
# # chunk_file <- "/scratch/schoelleh96/wp22a/stage2_chunks/chunk_10.rds"

# # cat("=== Loading files ===\n")
# # full_resid <- readRDS(resid_file)
# # chunk      <- readRDS(chunk_file)

# # nc <- nc_open(rh_file)
# # rh <- ncvar_get(nc, "r")
# # rh_dim_names <- sapply(nc$var[["r"]]$dim, function(d) d$name)
# # rh_dim_lens  <- sapply(nc$var[["r"]]$dim, function(d) d$len)
# # nc_close(nc)

# # res <- full_resid$residuals

# # cat("\n=== Full residual object ===\n")
# # cat("names(full_resid):\n")
# # print(names(full_resid))
# # cat("dim(full_resid$residuals): ", paste(dim(res), collapse = " x "), "\n", sep = "")
# # cat("length(full_resid$lon): ", length(full_resid$lon), "\n", sep = "")
# # cat("length(full_resid$lat): ", length(full_resid$lat), "\n", sep = "")

# # if (!is.null(full_resid$time)) {
# #   cat("length(full_resid$time): ", length(full_resid$time), "\n", sep = "")
# # } else {
# #   cat("full_resid$time: not present\n")
# # }

# # cat("\n=== RH NetCDF variable ===\n")
# # cat("var name: r\n")
# # cat("dim(rh): ", paste(dim(rh), collapse = " x "), "\n", sep = "")
# # cat("nc dim names: ", paste(rh_dim_names, collapse = ", "), "\n", sep = "")
# # cat("nc dim lens : ", paste(rh_dim_lens, collapse = " x "), "\n", sep = "")

# # cat("\n=== Chunk object ===\n")
# # cat("names(chunk):\n")
# # print(names(chunk))
# # cat("dim(chunk$residuals): ", paste(dim(chunk$residuals), collapse = " x "), "\n", sep = "")
# # cat("length(chunk$lon): ", length(chunk$lon), "\n", sep = "")
# # cat("length(chunk$lat): ", length(chunk$lat), "\n", sep = "")
# # cat("length(chunk$time): ", length(chunk$time), "\n", sep = "")

# # if (!("rh_700_850" %in% names(chunk))) {
# #   stop("chunk does not contain a variable named 'rh_700_850'")
# # }

# # cat("dim(chunk$rh_700_850): ", paste(dim(chunk$rh_700_850), collapse = " x "), "\n", sep = "")

# # cat("\n=== Checks that correspond to your two scripts ===\n")

# # # Script 1 logic: residual cube vs NC variable should match as [lon, lat, time]
# # same_dim_full_vs_nc <- identical(dim(res), dim(rh))
# # cat("full residual dims == RH nc dims: ", same_dim_full_vs_nc, "\n", sep = "")

# # if (!same_dim_full_vs_nc) {
# #   cat("WARNING: full residuals and RH nc variable do not have identical dims.\n")
# #   cat("This is often a dim-order issue, e.g. [lon,lat,time] vs [time,lat,lon].\n")
# # }

# # # Script 2 logic: within a chunk, residuals vs predictor must match exactly
# # same_dim_chunk_res_vs_rh <- identical(dim(chunk$residuals), dim(chunk$rh_700_850))
# # cat("chunk residual dims == chunk rh dims: ", same_dim_chunk_res_vs_rh, "\n", sep = "")

# # if (!same_dim_chunk_res_vs_rh) {
# #   cat("WARNING: chunk residuals and chunk rh do not have identical dims.\n")
# # }

# # cat("\n=== Length consistency checks ===\n")
# # cat("full: dim1 == length(lon): ", dim(res)[1] == length(full_resid$lon), "\n", sep = "")
# # cat("full: dim2 == length(lat): ", dim(res)[2] == length(full_resid$lat), "\n", sep = "")

# # if (!is.null(full_resid$time)) {
# #   cat("full: dim3 == length(time): ", dim(res)[3] == length(full_resid$time), "\n", sep = "")
# # }

# # cat("chunk: dim1 == length(lon): ", dim(chunk$residuals)[1] == length(chunk$lon), "\n", sep = "")
# # cat("chunk: dim2 == length(lat): ", dim(chunk$residuals)[2] == length(chunk$lat), "\n", sep = "")
# # cat("chunk: dim3 == length(time): ", dim(chunk$residuals)[3] == length(chunk$time), "\n", sep = "")

# # cat("\n=== Compare chunk against corresponding full-cube subset ===\n")

# # # Try to locate chunk lon/lat inside full lon/lat
# # lon_idx <- match(chunk$lon, full_resid$lon)
# # lat_idx <- match(chunk$lat, full_resid$lat)

# # cat("all chunk lon found in full lon: ", all(!is.na(lon_idx)), "\n", sep = "")
# # cat("all chunk lat found in full lat: ", all(!is.na(lat_idx)), "\n", sep = "")

# # if (all(!is.na(lon_idx)) && all(!is.na(lat_idx))) {
# #   res_sub <- res[lon_idx, lat_idx, , drop = FALSE]

# #   cat("dim(full residual subset): ", paste(dim(res_sub), collapse = " x "), "\n", sep = "")
# #   cat("subset dims identical to chunk residual dims: ",
# #       identical(dim(res_sub), dim(chunk$residuals)), "\n", sep = "")

# #   # Compare actual values for residuals
# #   res_equal <- isTRUE(all.equal(res_sub, chunk$residuals, check.attributes = FALSE))
# #   cat("full residual subset equals chunk residuals: ", res_equal, "\n", sep = "")

# #   if (!res_equal) {
# #     cat("Residual subset and chunk residuals differ.\n")
# #     cat("all.equal says:\n")
# #     print(all.equal(res_sub, chunk$residuals, check.attributes = FALSE))
# #   }

# #   # Compare RH nc subset to chunk rh, assuming nc RH is already [lon,lat,time]
# #   if (length(dim(rh)) == 3 &&
# #       max(lon_idx, na.rm = TRUE) <= dim(rh)[1] &&
# #       max(lat_idx, na.rm = TRUE) <= dim(rh)[2]) {

# #     rh_sub_direct <- rh[lon_idx, lat_idx, , drop = FALSE]
# #     cat("\ndim(RH nc direct subset [lon,lat,time] assumption): ",
# #         paste(dim(rh_sub_direct), collapse = " x "), "\n", sep = "")
# #     cat("direct RH subset dims identical to chunk rh dims: ",
# #         identical(dim(rh_sub_direct), dim(chunk$rh)), "\n", sep = "")

# #     rh_equal_direct <- isTRUE(all.equal(rh_sub_direct, chunk$rh_700_850, check.attributes = FALSE))
# #     cat("direct RH subset equals chunk rh: ", rh_equal_direct, "\n", sep = "")

# #     if (!rh_equal_direct) {
# #       cat("direct comparison all.equal says:\n")
# #       print(all.equal(rh_sub_direct, chunk$rh, check.attributes = FALSE))
# #     }
# #   }

# #   # Also try common alternative permutations
# #   perms <- list(
# #     c(1, 2, 3),
# #     c(2, 1, 3),
# #     c(3, 2, 1),
# #     c(3, 1, 2),
# #     c(2, 3, 1),
# #     c(1, 3, 2)
# #   )

# #   cat("\n=== Testing common RH dim permutations against chunk$rh ===\n")
# #   for (p in perms) {
# #     rhp <- aperm(rh, p)
# #     lbl <- paste0("[", paste(p, collapse = ","), "]")

# #     cat("\nPermutation ", lbl, " -> dim: ",
# #         paste(dim(rhp), collapse = " x "), "\n", sep = "")

# #     ok_shape <- length(dim(rhp)) == 3 &&
# #       max(lon_idx, na.rm = TRUE) <= dim(rhp)[1] &&
# #       max(lat_idx, na.rm = TRUE) <= dim(rhp)[2]

# #     cat("subset possible with [lon_idx, lat_idx, ]: ", ok_shape, "\n", sep = "")

# #     if (ok_shape) {
# #       rh_sub <- rhp[lon_idx, lat_idx, , drop = FALSE]
# #       same_dim <- identical(dim(rh_sub), dim(chunk$rh_700_850))
# #       same_val <- isTRUE(all.equal(rh_sub, chunk$rh_700_850, check.attributes = FALSE))

# #       cat("subset dims identical to chunk rh_700_850: ", same_dim, "\n", sep = "")
# #       cat("subset values equal to chunk rh_700_850: ", same_val, "\n", sep = "")

# #       if (same_val) {
# #         cat(">>> MATCH FOUND with permutation ", lbl, "\n", sep = "")
# #       }
# #     }
# #   }

# # } else {
# #   cat("Could not map chunk lon/lat back to full residual lon/lat, so subset comparison was skipped.\n")
# # }

# # cat("\n=== Matrix-shape check used by both scripts ===\n")
# # ngrid_full  <- dim(res)[1] * dim(res)[2]
# # ntime_full  <- dim(res)[3]
# # res2_full   <- matrix(res, nrow = ngrid_full, ncol = ntime_full)
# # cat("dim(matrix(full residuals)): ", paste(dim(res2_full), collapse = " x "), "\n", sep = "")

# # ngrid_chunk <- dim(chunk$residuals)[1] * dim(chunk$residuals)[2]
# # ntime_chunk <- dim(chunk$residuals)[3]
# # res2_chunk  <- matrix(chunk$residuals, nrow = ngrid_chunk, ncol = ntime_chunk)
# # rh2_chunk   <- matrix(chunk$rh_700_850, nrow = ngrid_chunk, ncol = ntime_chunk)

# # cat("dim(matrix(chunk residuals)): ", paste(dim(res2_chunk), collapse = " x "), "\n", sep = "")
# # cat("dim(matrix(chunk rh)): ", paste(dim(rh2_chunk), collapse = " x "), "\n", sep = "")
# # cat("matrix dims identical inside chunk: ",
# #     identical(dim(res2_chunk), dim(rh2_chunk)), "\n", sep = "")

# library(ncdf4)

# resid_file <- "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
# rh_file    <- "/scratch/schoelleh96/wp22a/data/relative_humidity_700-850.nc"
# chunk_file <- "/scratch/schoelleh96/wp22a/stage2_chunks/chunk_10.rds"

# resp  <- readRDS(resid_file)
# chunk <- readRDS(chunk_file)

# to_posix_utc <- function(x) {
#   if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
#   if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC") + 12 * 3600)
#   y <- as.POSIXct(x, tz = "UTC")
#   y
# }

# find_nc_coord_var <- function(nc, patterns) {
#   cand <- names(nc$dim)
#   hit <- cand[grepl(paste(patterns, collapse = "|"), cand, ignore.case = TRUE)]
#   if (length(hit) > 0) return(hit[1])

#   candv <- names(nc$var)
#   hitv <- candv[grepl(paste(patterns, collapse = "|"), candv, ignore.case = TRUE)]
#   if (length(hitv) > 0) return(hitv[1])

#   NULL
# }

# get_nc_coord <- function(nc, patterns) {
#   nm <- find_nc_coord_var(nc, patterns)
#   if (is.null(nm)) return(NULL)
#   if (nm %in% names(nc$dim)) return(nc$dim[[nm]]$vals)
#   ncvar_get(nc, nm)
# }

# parse_nc_time <- function(nc) {
#   time_name <- find_nc_coord_var(nc, c("^time$", "time"))
#   vals <- if (time_name %in% names(nc$dim)) nc$dim[[time_name]]$vals else ncvar_get(nc, time_name)
#   units_att <- ncatt_get(nc, time_name, "units")$value

#   parts <- strsplit(trimws(units_att), " since ", fixed = TRUE)[[1]]
#   unit <- tolower(trimws(parts[1]))
#   origin_txt <- trimws(parts[2])

#   origin <- as.POSIXct(origin_txt, tz = "UTC")
#   if (is.na(origin)) origin <- as.POSIXct(paste0(origin_txt, " 00:00:00"), tz = "UTC")

#   mult <- switch(unit,
#     "seconds" = 1,
#     "second"  = 1,
#     "secs"    = 1,
#     "sec"     = 1,
#     "hours"   = 3600,
#     "hour"    = 3600,
#     "days"    = 86400,
#     "day"     = 86400,
#     stop("Unsupported unit: ", unit)
#   )

#   origin + vals * mult
# }

# match_time_indices <- function(source_time, target_time) {
#   idx <- match(target_time, source_time)
#   if (all(!is.na(idx))) return(idx)

#   out <- integer(length(target_time))
#   for (i in seq_along(target_time)) {
#     same_day <- as.Date(source_time) == as.Date(target_time[i])
#     if (!any(same_day)) {
#       out[i] <- NA_integer_
#       next
#     }
#     cand <- which(same_day)
#     out[i] <- cand[which.min(abs(as.numeric(difftime(source_time[cand], target_time[i], units = "secs"))))]
#   }
#   out
# }

# nc <- nc_open(rh_file)
# rh <- ncvar_get(nc, "r")
# src_time <- parse_nc_time(nc)
# nc_close(nc)

# target_time <- to_posix_utc(resp$time)
# time_idx <- match_time_indices(src_time, target_time)

# lon_idx <- match(chunk$lon, resp$lon)
# lat_idx <- match(chunk$lat, resp$lat)

# rh_matched <- rh[lon_idx, lat_idx, time_idx, drop = FALSE]

# cat("identical dims: ", identical(dim(rh_matched), dim(chunk$rh_700_850)), "\n", sep = "")
# cat("all.equal matched-vs-chunk:\n")
# print(all.equal(rh_matched, chunk$rh_700_850, check.attributes = FALSE))

# cat("\nHow many exact time matches?\n")
# cat(sum(!is.na(match(target_time, src_time))), " / ", length(target_time), "\n", sep = "")

# cat("\nHow many duplicated source timesteps used?\n")
# cat(sum(duplicated(time_idx[!is.na(time_idx)])), "\n")
library(ncdf4)

resid_file <- "/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
rh_file    <- "/scratch/schoelleh96/wp22a/data/relative_humidity_700-850.nc"
chunk_file <- "/scratch/schoelleh96/wp22a/stage2_chunks/chunk_10.rds"

resp  <- readRDS(resid_file)
chunk <- readRDS(chunk_file)

to_posix_utc <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC") + 12 * 3600)
  y <- as.POSIXct(x, tz = "UTC")
  if (all(is.na(y))) stop("Could not parse time")
  y
}

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
  vals <- if (time_name %in% names(nc$dim)) nc$dim[[time_name]]$vals else ncvar_get(nc, time_name)
  units_att <- ncatt_get(nc, time_name, "units")$value

  parts <- strsplit(trimws(units_att), " since ", fixed = TRUE)[[1]]
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
    stop("Unsupported time unit: ", unit)
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

nc <- nc_open(rh_file)
rh <- ncvar_get(nc, "r")

src_lon <- get_nc_coord(nc, c("^lon$", "longitude"))
src_lat <- get_nc_coord(nc, c("^lat$", "latitude"))
src_time <- parse_nc_time(nc)

var <- nc$var[["r"]]
dim_names <- vapply(var$dim, function(d) d$name, character(1))
dim_lens  <- vapply(var$dim, function(d) d$len, integer(1))
nc_close(nc)

lon_pos  <- which(grepl("lon|longitude", dim_names, ignore.case = TRUE))[1]
lat_pos  <- which(grepl("lat|latitude", dim_names, ignore.case = TRUE))[1]
time_pos <- which(grepl("time", dim_names, ignore.case = TRUE))[1]

if (is.na(lon_pos))  lon_pos  <- which(dim_lens == length(src_lon))[1]
if (is.na(lat_pos))  lat_pos  <- which(dim_lens == length(src_lat))[1]
if (is.na(time_pos)) time_pos <- which(dim_lens == length(src_time))[1]

rh3 <- aperm(rh, c(lon_pos, lat_pos, time_pos))

target_lon  <- chunk$lon
target_lat  <- chunk$lat
target_time <- to_posix_utc(chunk$time)

trg_lon_src <- convert_lon_to_source(target_lon, src_lon)
lon_idx_src <- match_coord_indices(src_lon, trg_lon_src)
lat_idx_src <- match_coord_indices(src_lat, target_lat)
time_idx    <- match_time_indices(src_time, target_time)

cat("Exact time matches: ", sum(!is.na(match(target_time, src_time))), " / ", length(target_time), "\n", sep = "")
cat("All lon matched: ", all(!is.na(lon_idx_src)), "\n", sep = "")
cat("All lat matched: ", all(!is.na(lat_idx_src)), "\n", sep = "")

rh_sub_correct <- rh3[lon_idx_src, lat_idx_src, time_idx, drop = FALSE]

cat("dim(rh_sub_correct): ", paste(dim(rh_sub_correct), collapse = " x "), "\n", sep = "")
cat("dim(chunk$rh_700_850): ", paste(dim(chunk$rh_700_850), collapse = " x "), "\n", sep = "")

cat("all.equal(correct subset vs chunk):\n")
print(all.equal(rh_sub_correct, chunk$rh_700_850, check.attributes = FALSE))

cat("\nCoordinate diagnostics:\n")
cat("resp lon head/tail: ", head(resp$lon, 3), " ... ", tail(resp$lon, 3), "\n")
cat("src  lon head/tail: ", head(src_lon, 3), " ... ", tail(src_lon, 3), "\n")
cat("resp lat head/tail: ", head(resp$lat, 3), " ... ", tail(resp$lat, 3), "\n")
cat("src  lat head/tail: ", head(src_lat, 3), " ... ", tail(src_lat, 3), "\n")

cat("\nMatched source indices:\n")
print(lon_idx_src)
print(lat_idx_src)