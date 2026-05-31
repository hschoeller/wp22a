#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(emmeans)
  library(dplyr)
})

source("RScripts/config.r")

log_mem <- function(tag = "") {
  x <- readLines("/proc/self/status")
  vals <- grep("^(VmRSS|VmHWM|VmSize):", x, value = TRUE)

  to_gb <- function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    name <- sub(":", "", parts[1])
    kb <- as.numeric(parts[2])
    gb <- kb / 1024^2
    sprintf("%s: %.2f GB", name, gb)
  }

  cat(tag, "\n", sep = "")
  cat(vapply(vals, to_gb, character(1)), sep = "\n")
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript wcb_stage2.r <CHUNK_DIR> <CHUNK_INDEX> <OUT_DIR> <WCB_RDS> [PREDICTORS] [OVERWRITE]"
  )
}

is_boolish <- function(x) {
  tolower(x) %in% c("true", "false", "1", "0", "yes", "no", "y", "n")
}

default_predictors <- c("in_24", "asc_12", "out_00")

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2]) + 1
OUT_DIR    <- args[3]
WCB_RDS    <- args[4]

if (length(args) >= 5 && !is_boolish(args[5])) {
  predictors <- trimws(strsplit(args[5], "[,[:space:]]+")[[1]])
  predictors <- predictors[nzchar(predictors)]
  overwrite_arg_pos <- 6
} else {
  predictors <- default_predictors
  overwrite_arg_pos <- 5
}

if (length(predictors) == 0) {
  stop("No predictors selected.")
}

OVERWRITE <- if (length(args) >= overwrite_arg_pos) {
  tolower(args[overwrite_arg_pos]) %in% c("true", "1", "yes", "y")
} else {
  FALSE
}

change_points <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR, sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) stop("Chunk file not found: ", chunk_file)

outfile <- file.path(OUT_DIR, sprintf("wcb_stage2_chunk_%02d.rds", CHUNK_NO))
if (file.exists(outfile) && !OVERWRITE) {
  stop("Output exists and OVERWRITE=FALSE: ", outfile)
}

log_mem("start")

message("Reading chunk: ", chunk_file)
chunk <- readRDS(chunk_file)

message("Reading WCB predictor file: ", WCB_RDS)
wcb_df <- readRDS(WCB_RDS)
wcb_df$date <- as.Date(wcb_df$date)

needed_cols <- c("date", "wrname", predictors)
if (!all(needed_cols %in% names(wcb_df))) {
  stop("WCB file must contain columns: ", paste(needed_cols, collapse = ", "))
}

fit_gls <- function(formula, data) {
  nlme::gls(
    formula,
    data = data,
    method = "ML",
    correlation = corAR1(form = ~ day_no | segment),
    na.action = na.omit,
    control = glsControl(msMaxIter = 100, msVerbose = FALSE)
  )
}

trend_col_name <- function(df) {
  cand <- grep("\\.trend$|^trend$", names(df), value = TRUE)
  if (length(cand) > 0) return(cand[1])
  stop("Could not identify trend column in emtrends output.")
}

common_start <- max(min(as.Date(chunk$time), na.rm = TRUE), min(wcb_df$date, na.rm = TRUE))
common_end   <- min(max(as.Date(chunk$time), na.rm = TRUE), max(wcb_df$date, na.rm = TRUE))
common_idx <- which(as.Date(chunk$time) >= common_start & as.Date(chunk$time) <= common_end)

make_base_df <- function(chunk, wcb_df, predictors, common_start, common_end, common_idx) {
  chunk_dates <- data.frame(date = as.Date(chunk$time)[common_idx])
  left_join(wcb_df[, c("date", "wrname", predictors)], chunk_dates, by = "date")
}

base_df <- make_base_df(chunk, wcb_df, predictors, common_start, common_end, common_idx)

fit_gridpoint <- function(ii, jj, chunk, base_df, predictors) {
  d <- base_df
  d$resid <- d$resid <- chunk$residuals[ii, jj, common_idx]
  d$lon <- chunk$lon[ii]
  d$lat <- chunk$lat[jj]
  d$i_local <- ii
  d$j_local <- jj

  keep_cols <- c("resid", "wrname", predictors)
  d <- d[complete.cases(d[, keep_cols]), , drop = FALSE]
  if (nrow(d) < 80) return(NULL)
  if (sd(d$resid, na.rm = TRUE) == 0) return(NULL)
  if (length(unique(d$wrname)) < 2) return(NULL)
  d <- d[order(d$date), , drop = FALSE]

  cp_use <- change_points[change_points > min(d$date) & change_points < max(d$date)]
  d$segment <- cut(
    d$date,
    breaks = c(min(d$date), cp_use, max(d$date) + 1),
    labels = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )

  d <- dplyr::group_by(d, segment) |>
    dplyr::mutate(day_no = dplyr::row_number()) |>
    dplyr::ungroup()

  d$wrname <- factor(d$wrname)

for (pred in predictors) {
  cname <- paste0(pred, "_c")
  d[[cname]] <- d[[pred]] - ave(d[[pred]], d$wrname, FUN = function(x) mean(x, na.rm = TRUE))
}
  message("fitting wr only")
  fit_wr <- tryCatch(fit_gls(resid ~ wrname, d), error = function(e) NULL)

  centered_terms <- paste0(predictors, "_c")
fit_formula <- as.formula(
    paste(
        "resid ~ wrname + wrname*(",
        paste(centered_terms, collapse = " + "),
        ")"
    )
)
  message("fitting full")
  fit_full <- tryCatch(fit_gls(fit_formula, d), error = function(e) NULL)
  if (is.null(fit_wr) || is.null(fit_full)) return(NULL)

  emm_wr <- tryCatch(
    as.data.frame(summary(emmeans::emmeans(fit_wr, ~ wrname, data = d, mode = "asymptotic"),
                          infer = c(TRUE, TRUE))),
    error = function(e) NULL
  )

  at_list <- as.list(setNames(rep(0, length(centered_terms)), centered_terms))

  emm_full <- tryCatch(
    as.data.frame(summary(
      emmeans::emmeans(
        fit_full,
        ~ wrname,
        at = at_list,
        data = d,
        mode = "asymptotic"
      ),
      infer = c(TRUE, TRUE)
    )),
    error = function(e) NULL
  )

  if (is.null(emm_wr) || is.null(emm_full)) return(NULL)

  wr_means <- merge(
    emm_wr,
    emm_full,
    by = "wrname",
    suffixes = c("_wr", "_full"),
    all = TRUE
  )

  dev_wr <- wr_means$emmean_wr - mean(wr_means$emmean_wr, na.rm = TRUE)
  dev_full <- wr_means$emmean_full - mean(wr_means$emmean_full, na.rm = TRUE)

  var_wr <- stats::var(dev_wr, na.rm = TRUE)
  var_full <- stats::var(dev_full, na.rm = TRUE)

  wr_pattern_var_reduction <- if (is.finite(var_wr) && var_wr > 0 && is.finite(var_full)) {
    1 - var_full / var_wr
  } else {
    NA_real_
  }

  wr_means$dev_wr <- dev_wr
  wr_means$dev_full <- dev_full
  wr_means$delta_emmean <- wr_means$emmean_full - wr_means$emmean_wr
  wr_means$regime_attenuation <- ifelse(
    abs(wr_means$dev_wr) > 0,
    1 - abs(wr_means$dev_full) / abs(wr_means$dev_wr),
    NA_real_
  )

  wcb_effects <- list()
  kk <- 1L

  for (pred in centered_terms) {
    tr <- tryCatch(
      as.data.frame(summary(
        emmeans::emtrends(fit_full, specs = ~ wrname, var = pred, data = d, mode = "asymptotic"),
        infer = c(TRUE, TRUE)
      )),
      error = function(e) NULL
    )
    if (is.null(tr)) next

    tcol <- trend_col_name(tr)
    base_pred <- sub("_c$", "", pred)
    wr_sd_predictor <- tapply(
    d[[base_pred]],
    d$wrname,
    sd,
    na.rm = TRUE
)
resid_sd <- sd(d$resid, na.rm = TRUE)

tr$predictor       <- base_pred
tr$slope           <- tr[[tcol]]
tr$sd_predictor    <- stats::sd(d[[base_pred]], na.rm = TRUE)
tr$within_wr_sd    <- wr_sd_predictor[as.character(tr$wrname)]
tr$resid_sd        <- resid_sd
tr$effect_1sd      <- tr$slope * tr$sd_predictor
tr$effect_within_wr_1sd <- tr$slope * tr$within_wr_sd
tr$effect_std_beta <- tr$effect_1sd / resid_sd
    tr$scaled_effect <- tr[[tcol]] * tr$sd_predictor
    tr$lon <- d$lon[1]
    tr$lat <- d$lat[1]
    tr$i_local <- ii
    tr$j_local <- jj
    wcb_effects[[kk]] <- tr
    kk <- kk + 1L
  }

  # Combined effect: all predictors +1 SD together, with inference
# baseline (all predictors at mean)
message(sprintf("[ii=%d jj=%d] starting combined effect block", ii, jj))

  at_0 <- as.list(setNames(rep(0, length(centered_terms)), centered_terms))
  message(sprintf("[ii=%d jj=%d] at_0 built, length=%d, names=%s",
                  ii, jj, length(at_0),
                  paste(names(at_0), collapse = ",")))

  wr_levels <- levels(d$wrname)
  message(sprintf("[ii=%d jj=%d] wr_levels: %s",
                  ii, jj, paste(wr_levels, collapse = ",")))

  sd_by_wr <- lapply(predictors, function(pred) {
    tapply(d[[pred]], d$wrname, sd, na.rm = TRUE)
  })
  names(sd_by_wr) <- predictors
  message(sprintf("[ii=%d jj=%d] sd_by_wr computed for %d predictors",
                  ii, jj, length(sd_by_wr)))
  for (pred in predictors) {
    message(sprintf("  %s: %s",
                    pred,
                    paste(round(sd_by_wr[[pred]], 4), collapse = ",")))
  }

  any_na_sd <- any(vapply(sd_by_wr,
                          function(v) any(!is.finite(v)),
                          logical(1)))
  if (any_na_sd) {
    message(sprintf("[ii=%d jj=%d] WARNING: non-finite SDs present", ii, jj))
  }

beta <- coef(fit_full)
V    <- vcov(fit_full)

combined_effects <- tryCatch({
  rows <- lapply(wr_levels, function(wr) {
    a <- setNames(numeric(length(beta)), names(beta))

    for (k in seq_along(centered_terms)) {
      ct  <- centered_terms[k]
      pr  <- predictors[k]
      sd_l <- sd_by_wr[[pr]][[wr]]

      main_term        <- ct
      interaction_term <- paste0("wrname", wr, ":", ct)

      if (main_term %in% names(a)) {
        a[main_term] <- a[main_term] + sd_l
      }
      if (interaction_term %in% names(a)) {
        a[interaction_term] <- a[interaction_term] + sd_l
      }
    }

    est  <- sum(a * beta)
    se   <- sqrt(as.numeric(t(a) %*% V %*% a))
    z    <- est / se
    p    <- 2 * pnorm(-abs(z))
    lcl  <- est - 1.96 * se
    ucl  <- est + 1.96 * se

    data.frame(
      wrname               = wr,
      combined_effect_1sd  = est,
      se                   = se,
      z                    = z,
      p.value              = p,
      asymp.LCL            = lcl,
      asymp.UCL            = ucl
    )
  })
  do.call(rbind, rows)
}, error = function(e) {
  message(sprintf("[ii=%d jj=%d] manual contrast FAILED: %s",
                  ii, jj, conditionMessage(e)))
  NULL
})

if (!is.null(combined_effects)) {
  combined_effects$lon     <- d$lon[1]
  combined_effects$lat     <- d$lat[1]
  combined_effects$i_local <- ii
  combined_effects$j_local <- jj
}

  wcb_effects <- bind_rows(wcb_effects)
    wcb_effects <- wcb_effects[, !grepl("_c\\.trend$", names(wcb_effects))]
  list(
    summary = data.frame(
      lon = d$lon[1],
      lat = d$lat[1],
      i_local = ii,
      j_local = jj,
      n = nrow(d),
      n_wr = nlevels(d$wrname),
      wr_pattern_var_reduction = wr_pattern_var_reduction,
      wr_only_mean_abs_dev = mean(abs(dev_wr), na.rm = TRUE),
      full_mean_abs_dev = mean(abs(dev_full), na.rm = TRUE)
    ),
    wr_means = transform(
      wr_means,
      lon = d$lon[1],
      lat = d$lat[1],
      i_local = ii,
      j_local = jj
    ),
    wcb_effects = wcb_effects,
        combined_effects = combined_effects
  )
}

nx <- length(chunk$lon)
ny <- length(chunk$lat)

bind_safe <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (length(x) == 0) return(data.frame())
  bind_rows(x)
}

res_list <- vector("list", nx * ny)
kk <- 1L


for (ii in seq_len(nx)) {
  for (jj in seq_len(ny)) {
    # res_list[[kk]] <- tryCatch(
    #   fit_gridpoint(ii, jj, chunk, base_df, predictors),
    #   error = function(e) NULL
    # )
    res_list[[kk]] <- fit_gridpoint(ii, jj, chunk, base_df, predictors)
    kk <- kk + 1L
  }
}

res_list <- Filter(Negate(is.null), res_list)

out <- list(
  metadata = list(
    chunk_no = CHUNK_NO,
    chunk_file = chunk_file,
    source_chunk_metadata = chunk$metadata,
    wcb_rds = WCB_RDS,
    predictors = predictors,
    centered_within_gridpoint = TRUE
  ),
  summary = bind_safe(lapply(res_list, `[[`, "summary")),
  wr_means = bind_safe(lapply(res_list, `[[`, "wr_means")),
  wcb_effects = bind_safe(lapply(res_list, `[[`, "wcb_effects")),
    combined_effects = bind_safe(lapply(res_list, `[[`, "combined_effects"))
)

saveRDS(out, outfile, compress = "xz")
message("Wrote: ", outfile)