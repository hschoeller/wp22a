#!/usr/bin/env Rscript

## -----------------------------
## configuration
## -----------------------------
source("RScripts/config.r")
source("RScripts/data_functions.r")
source("RScripts/algo_functions.r")

LM_DIR <- "/scratch/schoelleh96/wp22a/ens_data/mods/"
# OUT_FILE <- file.path(LM_DIR, "gls_model_summary.rds")
OUT_FILE <- "/home/schoelleh96/wp22a/gls_model_summary.rds"
# number of parallel workers (can be overridden by env var)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))

## required globals (must exist)
YEAR_BOUND <- c(1940, 2024)

## -----------------------------
## summarize_gls_folder
## -----------------------------

summarize_gls_folder <- function(folder, n_cores = 1L) {
  files <- list.files(folder, pattern = "^lm.*\\.Rds$", full.names = TRUE)
  if (length(files) == 0) stop("No model files found")

  extract_nums <- function(x) {
    m <- gregexpr("-?\\d+\\.?\\d*", x)
    regmatches(x, m)[[1]]
  }

  worker <- function(fpath) {
    fname <- basename(fpath)
    nums <- extract_nums(fname)
    lat <- as.numeric(nums[1])
    lon <- as.numeric(nums[2])

    mod <- readRDS(fpath)

    dat <- mod$data
    if (is.null(dat)) {
      return(tibble::tibble(
        lat = lat, lon = lon,
        obs_var = NA_real_,
        sigma = NA_real_,
        aic = AIC(mod),
        bic = BIC(mod),
        fit1940 = NA_real_,
        fit2024 = NA_real_,
        fitted_diff = NA_real_,
        RSS = NA_real_,
        TSS = NA_real_,
        r_squared = NA_real_,
        adj_r_squared = NA_real_,
        coefs = list(NULL),
        GRSS = NA_real_,
        GTSS = NA_real_,
        n = NA_real_,
        p = NA_real_
      ))
    }

    obs_var <- var(dat$log_variance, na.rm = TRUE)

    ## response residuals and RSS (existing)
    res <- as.numeric(residuals(mod, type = "response"))
    RSS <- sum(res^2, na.rm = TRUE)

    TSS <- sum((dat$log_variance - mean(dat$log_variance, na.rm = TRUE))^2)

    n <- nrow(dat)
    p <- if (!is.null(N_COEFFS)) N_COEFFS else length(coef(mod))

    sigma <- sqrt(RSS / (n - p))
    r_squared <- 1 - RSS / TSS
    adj_r_squared <- 1 - ((RSS / (n - p)) / (TSS / (n - 1)))

    fit_vals <- numeric(2)
    for (i in seq_along(YEAR_BOUND)) {
      nd <- dat
      nd$year <- YEAR_BOUND[i]
      fit_vals[i] <- mean(predict(mod, newdata = nd), na.rm = TRUE)
    }

    s <- summary(mod)
    coef_tbl <- as.data.frame(s$tTable)
    names(coef_tbl) <- c("estimate", "std.error", "t.value", "p.value")
    coef_tbl$term <- rownames(coef_tbl)
    coef_tbl <- coef_tbl[, c("term", "estimate", "std.error", "t.value", "p.value")]

    ## --- NEW: compute GRSS and GTSS (GLS-whitened sums of squares) ---
    ## Try fast path: use normalized residuals (+ normalized fitted) if available
    GRSS <- NA_real_
    GTSS <- NA_real_

    res_norm <- tryCatch(
      as.numeric(residuals(mod, type = "normalized")),
      error = function(e) NULL
    )

    fit_norm <- tryCatch(
      as.numeric(fitted(mod, type = "normalized")),
      error = function(e) NULL
    )

    if (!is.null(res_norm) && !is.null(fit_norm)) {
      ## normalized residuals exist and normalized fitted exists -> direct computation
      GRSS <- sum(res_norm^2, na.rm = TRUE)
      y_star <- res_norm + fit_norm
      GTSS <- sum((y_star - mean(y_star, na.rm = TRUE))^2, na.rm = TRUE)
    } else if (!is.null(res_norm)) {
      ## At minimum we can compute GRSS from normalized residuals
      GRSS <- sum(res_norm^2, na.rm = TRUE)
      ## Attempt fallback to compute GTSS by constructing Sigma and whitening y
      gtss_try <- tryCatch({
        Sigma <- nlme::getVarCov(mod, type = "conditional")
        # If Sigma is block-diagonal or more complex, getVarCov may fail or return
        # per-group matrices; wrap in tryCatch
        if (!is.matrix(Sigma)) stop("getVarCov did not return a matrix")
        W <- chol(solve(Sigma))
        y <- tryCatch({
          # try common locations for the response vector
          if (!is.null(mod$model) && "log_variance" %in% names(mod$model)) {
            as.numeric(mod$model$log_variance)
          } else if (!is.null(dat$log_variance)) {
            as.numeric(dat$log_variance)
          } else stop("cannot locate response vector for whitening")
        }, error = function(e) stop(e))
        y_star <- as.numeric(W %*% y)
        sum((y_star - mean(y_star, na.rm = TRUE))^2, na.rm = TRUE)
      }, error = function(e) NA_real_)
      GTSS <- gtss_try
    } else {
      ## Neither normalized residuals nor normalized fitted available -> attempt full whitening
      full_try <- tryCatch({
        Sigma <- nlme::getVarCov(mod, type = "conditional")
        if (!is.matrix(Sigma)) stop("getVarCov did not return a matrix")
        W <- chol(solve(Sigma))
        y <- if (!is.null(mod$model) && "log_variance" %in% names(mod$model)) {
          as.numeric(mod$model$log_variance)
        } else if (!is.null(dat$log_variance)) {
          as.numeric(dat$log_variance)
        } else stop("cannot locate response vector for whitening")
        yhat <- as.numeric(fitted(mod)) # fallback fitted in response space
        y_star <- as.numeric(W %*% y)
        yhat_star <- as.numeric(W %*% yhat)
        GRSS_calc <- sum((y_star - yhat_star)^2, na.rm = TRUE)
        GTSS_calc  <- sum((y_star - mean(y_star, na.rm = TRUE))^2, na.rm = TRUE)
        list(grss = GRSS_calc, gtss = GTSS_calc)
      }, error = function(e) NULL)

      if (!is.null(full_try)) {
        GRSS <- full_try$grss
        GTSS <- full_try$gtss
      } else {
        ## leave GRSS and GTSS as NA
        GRSS <- NA_real_
        GTSS <- NA_real_
      }
    }
    phi <- unname(coef(mod$modelStruct$corStruct, unconstrained = FALSE))
    n_eff <- n * (1 - phi) / (1 + phi)
    tibble::tibble(
      lat = lat,
      lon = lon,
      obs_var = obs_var,
      sigma = sigma,
      aic = AIC(mod),
      bic = BIC(mod),
      fit1940 = fit_vals[1],
      fit2024 = fit_vals[2],
      fitted_diff = fit_vals[1] - fit_vals[2],
      RSS = RSS,
      TSS = TSS,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      coefs = list(coef_tbl),
      GRSS = GRSS,
      GTSS = GTSS,
      n = n,
      p = p,
      n_eff = n_eff,
      phi = phi
    )
  }

  message("Summarizing ", length(files), " models using ", n_cores, " cores")

  rows <- parallel::mclapply(
    files,
    worker,
    mc.cores = n_cores
  )

  dplyr::bind_rows(rows)
}


## -----------------------------
## run + save
## -----------------------------

res <- summarize_gls_folder(LM_DIR, n_cores = N_CORES)
res <- add_fdr_pvalues_to_coefs(res, method = "fdr")

saveRDS(res, OUT_FILE)

message("Saved combined summary to: ", OUT_FILE)
