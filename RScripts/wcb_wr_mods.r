#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

WR_RDS <- get_arg("--wr-rds", Sys.getenv("WR_RDS"))

TASK_ID <- as.integer(get_arg(
  "--task-id",
  Sys.getenv("SLURM_ARRAY_TASK_ID", "1")
))

IN_DIR <- get_arg(
  "--in-dir",
  Sys.getenv("IN_DIR", "./stage2_chunks_res")
)

OUT_DIR <- get_arg(
  "--out-dir",
  Sys.getenv("OUT_DIR", "./wcb_wr_gls_results")
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CHUNK_RDS <- file.path(IN_DIR, sprintf("chunk_%02d.rds", TASK_ID))
OUT_RDS   <- file.path(OUT_DIR, sprintf("chunk_%02d.rds", TASK_ID))

if (is.null(CHUNK_RDS) || CHUNK_RDS == "" || !file.exists(CHUNK_RDS)) {
  stop("Missing or unreadable chunk file: ", CHUNK_RDS)
}
if (is.null(WR_RDS) || WR_RDS == "" || !file.exists(WR_RDS)) {
  stop("Missing or unreadable --wr-rds")
}

message("Reading chunk: ", CHUNK_RDS)
chunk <- readRDS(CHUNK_RDS)

needed <- c("residuals", "lon", "lat", "time",
            "wcb_in_12utc", "wcb_asc_12utc", "wcb_out_12utc")
missing_needed <- setdiff(needed, names(chunk))
if (length(missing_needed) > 0L) {
  stop("Chunk is missing required fields: ", paste(missing_needed, collapse = ", "))
}

message("Reading WR series: ", WR_RDS)
wr_obj <- readRDS(WR_RDS)

load_wr_series <- function(x) {
  if (is.data.frame(x)) {
    if (!all(c("date", "wrname") %in% names(x))) {
      stop("WR data frame must contain columns named 'date' and 'wrname'.")
    }
    out <- data.frame(
      date = as.Date(x$date),
      wrname = as.factor(x$wrname),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  if (is.vector(x) && !is.list(x)) {
    if (is.null(names(x))) {
      stop("Named vector WR input must have date names.")
    }
    out <- data.frame(
      date = as.Date(names(x)),
      wrname = as.factor(as.character(x)),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  if (is.list(x) && all(c("date", "wrname") %in% names(x))) {
    out <- data.frame(
      date = as.Date(x$date),
      wrname = as.factor(x$wrname),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  stop("Unsupported WR object. Use a data frame with date/wrname, a named vector, or a list with date/wrname.")
}

wr_df <- load_wr_series(wr_obj)
wr_df <- wr_df[!is.na(wr_df$date) & !is.na(wr_df$wrname), , drop = FALSE]
wr_df <- wr_df[!duplicated(wr_df$date), , drop = FALSE]
wr_levels <- levels(factor(wr_df$wrname))

wr_map <- setNames(as.character(wr_df$wrname), as.character(wr_df$date))

fit_gls <- function(formula, data) {
  gls(
    formula,
    data = data,
    method = "ML",
    correlation = corAR1(form = ~ day_no | segment),
    na.action = na.omit,
    control = glsControl(msMaxIter = 100, msVerbose = FALSE)
  )
}

contrast_from_fit <- function(fit, wr_level, inf = 0, asc = 0, out = 0) {
  tt <- delete.response(terms(fit))
  beta <- coef(fit)
  V <- vcov(fit)

  nd0 <- data.frame(
    wr = factor(wr_level, levels = wr_levels),
    inf = 0, asc = 0, out = 0
  )
  nd1 <- data.frame(
    wr = factor(wr_level, levels = wr_levels),
    inf = inf, asc = asc, out = out
  )

  X0 <- model.matrix(tt, nd0, xlev = fit$xlevels)
  X1 <- model.matrix(tt, nd1, xlev = fit$xlevels)

  X0 <- X0[, names(beta), drop = FALSE]
  X1 <- X1[, names(beta), drop = FALSE]

  cvec <- drop(X1 - X0)
  est  <- as.numeric(cvec %*% beta)
  se   <- sqrt(as.numeric(cvec %*% V %*% cvec))

  if (!is.finite(se) || se <= 0) {
    return(data.frame(
      estimate  = est,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value   = NA_real_,
      conf_low  = NA_real_,
      conf_high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  tval <- est / se
  pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
  crit <- qnorm(0.975)

  data.frame(
    estimate  = est,
    std_error = se,
    statistic = tval,
    p_value   = pval,
    conf_low  = est - crit * se,
    conf_high = est + crit * se,
    stringsAsFactors = FALSE
  )
}

# Precompute time-dependent objects once.
time_vals <- as.POSIXct(chunk$time, tz = "UTC")
date_vals <- as.Date(time_vals, tz = "UTC")
wr_per_time <- unname(wr_map[as.character(date_vals)])

extract_cell_data <- function(lon_i, lat_i) {
  residual <- as.numeric(chunk$residuals[lon_i, lat_i, ])
  inf_wcb  <- as.numeric(chunk$wcb_in_12utc[lon_i, lat_i, ])
  asc_wcb  <- as.numeric(chunk$wcb_asc_12utc[lon_i, lat_i, ])
  out_wcb  <- as.numeric(chunk$wcb_out_12utc[lon_i, lat_i, ])

  keep <- is.finite(residual) &
    !is.na(wr_per_time) &
    is.finite(inf_wcb) &
    is.finite(asc_wcb) &
    is.finite(out_wcb)

  if (!any(keep)) return(NULL)

  dat <- data.frame(
    time = time_vals[keep],
    residual = residual[keep],
    wr = factor(wr_per_time[keep], levels = wr_levels),
    inf = inf_wcb[keep],
    asc = asc_wcb[keep],
    out = out_wcb[keep],
    stringsAsFactors = FALSE
  )

  if (nrow(dat) == 0L) return(NULL)

  dat <- dat[order(dat$time), , drop = FALSE]
  dat$day_no <- seq_len(nrow(dat))
  dat$segment <- factor(1L)

  dat
}

fit_cell <- function(lon_i, lat_i) {
  dat <- extract_cell_data(lon_i, lat_i)

  if (is.null(dat) || nrow(dat) < 10L) {
    return(list(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      nobs = 0L,
      skipped = TRUE,
      reason = "Too few complete observations after joining WR data."
    ))
  }

  if (nlevels(droplevels(dat$wr)) < 2L) {
    return(list(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      nobs = nrow(dat),
      skipped = TRUE,
      reason = "Fewer than 2 WR levels present in this cell."
    ))
  }

  fit <- tryCatch(
    fit_gls(residual ~ wr * (inf + asc + out), data = dat),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      nobs = nrow(dat),
      skipped = TRUE,
      reason = paste("GLS failed:", fit$message)
    ))
  }

  effect_rows <- vector("list", length(wr_levels) * 4L)
  k <- 0L

  for (wr_level in levels(dat$wr)) {
    k <- k + 1L
    effect_rows[[k]] <- cbind(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      wr = wr_level,
      effect = "inf",
      contrast_from_fit(fit, wr_level, inf = 1, asc = 0, out = 0)
    )

    k <- k + 1L
    effect_rows[[k]] <- cbind(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      wr = wr_level,
      effect = "asc",
      contrast_from_fit(fit, wr_level, inf = 0, asc = 1, out = 0)
    )

    k <- k + 1L
    effect_rows[[k]] <- cbind(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      wr = wr_level,
      effect = "out",
      contrast_from_fit(fit, wr_level, inf = 0, asc = 0, out = 1)
    )

    k <- k + 1L
    effect_rows[[k]] <- cbind(
      lon = as.numeric(chunk$lon[lon_i]),
      lat = as.numeric(chunk$lat[lat_i]),
      wr = wr_level,
      effect = "all",
      contrast_from_fit(fit, wr_level, inf = 1, asc = 1, out = 1)
    )
  }

  effect_table <- do.call(rbind, effect_rows)

  # Explicitly drop large objects as soon as possible.
  rm(dat, fit, effect_rows)
  invisible(gc(FALSE))

  list(
    effect_table = effect_table
  )
}

lon_n <- length(chunk$lon)
lat_n <- length(chunk$lat)
n_cells <- lon_n * lat_n

message("Fitting ", n_cells, " grid-point model(s)...")

effect_tables <- vector("list", n_cells)
skip_table <- vector("list", n_cells)

idx <- 0L
for (i in seq_len(lon_n)) {
  for (j in seq_len(lat_n)) {
    idx <- idx + 1L
    message(
      "  Cell ", idx, "/", n_cells,
      "  (lon=", chunk$lon[i], ", lat=", chunk$lat[j], ")"
    )

    res <- fit_cell(i, j)

    if (is.null(res$effect_table)) {
      skip_table[[idx]] <- data.frame(
        lon = as.numeric(chunk$lon[i]),
        lat = as.numeric(chunk$lat[j]),
        nobs = NA_integer_,
        skipped = TRUE,
        reason = "Skipped or failed before effect extraction.",
        stringsAsFactors = FALSE
      )
    } else {
      effect_tables[[idx]] <- res$effect_table
      skip_table[[idx]] <- NULL
    }

    rm(res)
    if (idx %% 2L == 0L) invisible(gc(FALSE))
  }
}

effect_tables <- Filter(Negate(is.null), effect_tables)
skip_table <- Filter(Negate(is.null), skip_table)

effect_table_all <- if (length(effect_tables) > 0L) do.call(rbind, effect_tables) else NULL
skip_table_all   <- if (length(skip_table) > 0L) do.call(rbind, skip_table) else NULL

out <- list(
  metadata = list(
    chunk_file = CHUNK_RDS,
    wr_file = WR_RDS,
    formula = "residual ~ wr * (inf + asc + out)",
    correlation = "AR1 over day_no within segment",
    model_scope = "one GLS model per grid point",
    wr_levels = wr_levels,
    n_grid_points = n_cells
  ),
  effect_table = effect_table_all,
  skipped_cells = skip_table_all
)

saveRDS(out, OUT_RDS, compress = "xz")
message("Wrote: ", OUT_RDS)