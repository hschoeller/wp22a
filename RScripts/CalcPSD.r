#!/usr/bin/env Rscript

source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript psd_summary_chunks_mc.R <chunk_dir> <nc_var> <out_dir> [n_cores]")
}

CHUNK_DIR <- args[1]
NC_VAR    <- args[2]
OUT_DIR   <- args[3]
N_CORES   <- if (length(args) >= 4) as.integer(args[4]) else 32L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH_RDS <- file.path(OUT_DIR, "psd_summary_1000bins.rds")
OUT_PATH_CSV <- file.path(OUT_DIR, "psd_summary_1000bins.csv")

# -------------------------
# Load coeff tibble -> beta_wide
# -------------------------
COEFF_RDS <- "/home/schoelleh96/wp22a/gls_model_summary.rds"
coeff_tbl <- readRDS(COEFF_RDS) %>%
  dplyr::select(lat, lon, coefs)

beta_long <- coeff_tbl %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::select(id, lat, lon, coefs) %>%
  tidyr::unnest(coefs) %>%                 # expects term, estimate
  dplyr::transmute(id, lat, lon, term, estimate)

all_terms <- sort(unique(beta_long$term))

beta_wide <- beta_long %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

missing_cols <- setdiff(all_terms, names(beta_wide))
if (length(missing_cols) > 0) beta_wide[missing_cols] <- 0

beta_wide <- beta_wide %>%
  dplyr::select(lat, lon, dplyr::all_of(all_terms))

# Faster lookup key
beta_wide <- beta_wide %>% mutate(key = paste(lat, lon))
beta_mat <- as.matrix(beta_wide[, all_terms, drop = FALSE])
rownames(beta_mat) <- beta_wide$key

# -------------------------
# CPs -> dates
# -------------------------
change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

# -------------------------
# Discover chunks
# -------------------------
chunk_files <- sort(Sys.glob(file.path(CHUNK_DIR, "chunk_*.nc")))
if (length(chunk_files) == 0) stop("No chunks found in: ", CHUNK_DIR)

# -------------------------
# Infer time axis & build X + bins from first chunk
# -------------------------
nc0 <- ncdf4::nc_open(chunk_files[1])
on.exit(ncdf4::nc_close(nc0), add = TRUE)

# time var name assumption consistent with earlier code; adjust if needed
time_data <- ncdf4::ncvar_get(nc0, "valid_time")
time_origin <- sub("seconds since ", "", ncdf4::ncatt_get(nc0, "valid_time", "units")$value)
dates <- as.Date(as.POSIXct(time_data, origin = time_origin, tz = "UTC"))
ntime <- length(dates)

ncdf4::nc_close(nc0)

# Build design matrix X once (uses your make_design_data)
# make_design_data expects time, z, cps and creates segment/year/sin/cos etc.
tmp_df <- make_design_data(time = dates, z = rep(1, ntime), cps = change_points)
form <- log_variance ~ segment + segment:year + segment:sin_doy + segment:cos_doy - 1
X <- model.matrix(form, data = tmp_df)

# align X columns to all_terms
miss_X <- setdiff(all_terms, colnames(X))
if (length(miss_X) > 0) {
  X <- cbind(X, matrix(0, nrow(X), length(miss_X), dimnames = list(NULL, miss_X)))
}
X <- X[, all_terms, drop = FALSE]

# -------------------------
# Frequency bins (log-binned)
# -------------------------
n_freq_target <- 1000L

freq_full <- (1:(floor(ntime / 2))) / ntime
freq_full <- freq_full[freq_full > 0]

freq_bins <- exp(seq(log(min(freq_full)), log(max(freq_full)), length.out = n_freq_target + 1))
freq_mid  <- sqrt(freq_bins[-1] * freq_bins[-length(freq_bins)])

bin_id <- findInterval(freq_full, freq_bins, rightmost.closed = FALSE, all.inside = TRUE)

compute_binned_psd <- function(x, bin_id, n_bins) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 512) return(rep(NA_real_, n_bins))

  x <- x - mean(x)

  Xf <- fft(x)
  P  <- (Mod(Xf)^2) / n
  P  <- P[2:(floor(n / 2) + 1)]   # positive freqs, drop DC

  s   <- rowsum(P, group = bin_id, reorder = FALSE)
  cts <- rowsum(rep(1, length(P)), group = bin_id, reorder = FALSE)

  out <- rep(NA_real_, n_bins)
  out[as.integer(rownames(s))] <- as.numeric(s) / as.numeric(cts)
  out
}

# -------------------------
# Chunk worker
# -------------------------
process_chunk <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  lon <- as.numeric(ncdf4::ncvar_get(nc, "longitude"))
  lat <- as.numeric(ncdf4::ncvar_get(nc, "latitude"))
  z   <- ncdf4::ncvar_get(nc, NC_VAR)

  d <- dim(z)
  if (length(d) == 3) {
    nlon <- d[1]; nlat <- d[2]; nt <- d[3]
    if (nt != ntime) stop("ntime mismatch in ", f, " got ", nt, " expected ", ntime)

    z_mat <- matrix(z, nrow = nlon * nlat, ncol = nt) # points x time
    lon_rep <- rep(lon, times = nlat)
    lat_rep <- rep(lat, each  = nlon)
  } else if (length(d) == 2) {
    npt <- d[1]; nt <- d[2]
    if (nt != ntime) stop("ntime mismatch in ", f, " got ", nt, " expected ", ntime)

    z_mat <- z
    lon_rep <- lon
    lat_rep <- lat
  } else {
    stop("Unsupported dim(z) in ", f, ": ", paste(d, collapse = "x"))
  }

  keys <- paste(lat_rep, lon_rep)
  have <- keys %in% rownames(beta_mat)
  if (!any(have)) return(NULL)

  keys_use <- keys[have]
  z_mat <- z_mat[have, , drop = FALSE]

  # obs/resid in time x points
  obs_mat <- log(z_mat)
  obs_mat[!is.finite(obs_mat)] <- NA_real_

  B <- beta_mat[keys_use, , drop = FALSE]        # points x p
  fit_tp <- X %*% t(B)                           # time x points
  obs_tp <- t(obs_mat)                           # time x points
  res_tp <- obs_tp - fit_tp

  npts <- ncol(obs_tp)
  psd_obs <- matrix(NA_real_, nrow = n_freq_target, ncol = npts)
  psd_res <- matrix(NA_real_, nrow = n_freq_target, ncol = npts)

  for (j in seq_len(npts)) {
    psd_obs[, j] <- compute_binned_psd(obs_tp[, j], bin_id, n_freq_target)
    psd_res[, j] <- compute_binned_psd(res_tp[, j], bin_id, n_freq_target)
  }

  list(psd_obs = psd_obs, psd_res = psd_res)
}

# -------------------------
# Parallel over chunks (mcapply)
# -------------------------
N_CHUNK_CORES <- min(N_CORES, 32L)  # avoid FS meltdown
message("Chunks: ", length(chunk_files), " | mc.cores=", N_CHUNK_CORES)

res_chunks <- parallel::mclapply(chunk_files, process_chunk, mc.cores = N_CHUNK_CORES)

# filter NULLs
res_chunks <- Filter(Negate(is.null), res_chunks)
if (length(res_chunks) == 0) stop("All chunks returned NULL (no matching betas?)")

PSD_obs_all <- do.call(cbind, lapply(res_chunks, `[[`, "psd_obs"))
PSD_res_all <- do.call(cbind, lapply(res_chunks, `[[`, "psd_res"))

mean_obs <- rowMeans(PSD_obs_all, na.rm = TRUE)
iqr_obs  <- apply(PSD_obs_all, 1, IQR, na.rm = TRUE)
mean_res <- rowMeans(PSD_res_all, na.rm = TRUE)
iqr_res  <- apply(PSD_res_all, 1, IQR, na.rm = TRUE)

out_df <- tibble::tibble(
  freq_mid = freq_mid,
  mean_obs = mean_obs,
  iqr_obs  = iqr_obs,
  mean_res = mean_res,
  iqr_res  = iqr_res
)

saveRDS(out_df, OUT_PATH_RDS)
write.csv(out_df, OUT_PATH_CSV, row.names = FALSE)

message("Saved:\n  ", OUT_PATH_RDS, "\n  ", OUT_PATH_CSV)
