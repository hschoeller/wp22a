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
    sprintf("%s: %.2f GB", name, kb / 1024^2)
  }
  cat(tag, "\n", sep = "")
  cat(vapply(vals, to_gb, character(1)), sep = "\n")
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript second_stage_chunk.R <CHUNK_DIR> <CHUNK_INDEX> <OUT_DIR> [WR_RDS] [OVERWRITE]")

CHUNK_DIR <- args[1]
CHUNK_NO  <- as.integer(args[2]) + 1
OUT_DIR   <- args[3]
WR_RDS    <- if (length(args) >= 4) args[4] else "/home/schoelleh96/wp22a/data/wrnames.rds"
OVERWRITE <- if (length(args) >= 5) tolower(args[5]) %in% c("true","1","yes","y") else FALSE

change_points <- as.Date(paste0(CP,"-01"), format="%Y-%m-%d")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chunk_file <- file.path(CHUNK_DIR, sprintf("chunk_%02d.rds", CHUNK_NO))
if (!file.exists(chunk_file)) stop("Chunk file not found: ", chunk_file)
outfile <- file.path(OUT_DIR, sprintf("stage2_estimates_chunk_%02d.rds", CHUNK_NO))

log_mem("start")
message("Reading chunk: ", chunk_file)
chunk <- readRDS(chunk_file)
message("Reading WR file: ", WR_RDS)
wr_min <- readRDS(WR_RDS)
wr_min$date <- as.Date(wr_min$date)

candidate_predictors <- c("rh_500_850","upper_wind","eady","grad")
predictors <- candidate_predictors[candidate_predictors %in% names(chunk)]
if (length(predictors)==0) stop("No predictors found in chunk.")

message("Predictors used: ", paste(predictors, collapse=", "))

fit_gls <- function(formula, data) {
  gls(formula, data=data, method="ML", correlation=corAR1(form=~day_no|segment),
      na.action=na.omit, control=glsControl(msMaxIter=100, msVerbose=FALSE))
}

pattern_metrics <- function(emm_ref, emm_adj){
  m0 <- emm_ref$emmean
  m1 <- emm_adj$emmean
  dev0 <- m0 - mean(m0, na.rm=TRUE)
  dev1 <- m1 - mean(m1, na.rm=TRUE)
  v0 <- var(dev0, na.rm=TRUE)
  v1 <- var(dev1, na.rm=TRUE)
  mad0 <- mean(abs(dev0), na.rm=TRUE)
  mad1 <- mean(abs(dev1), na.rm=TRUE)
  list(
    var_reduction = if(is.finite(v0) && v0>0 && is.finite(v1)) 1-v1/v0 else NA_real_,
    mad_reduction = if(is.finite(mad0) && mad0>0 && is.finite(mad1)) 1-mad1/mad0 else NA_real_,
    mad0 = mad0, mad1 = mad1
  )
}

base_df <- data.frame(
  time     = as.POSIXct(chunk$time, tz = "UTC"),
  date     = as.Date(chunk$time),
  time_idx = seq_along(chunk$time)          # <-- track original position
) |>
  left_join(wr_min, by = c("date" = "date")) |>
  filter(!is.na(wrname)) |>
  mutate(
    segment = cut(
      date,
      breaks        = c(date[1], change_points, date[length(date)] + 1),
      labels        = FALSE,
      include.lowest = TRUE,
      right         = FALSE
    )
  ) |>
  group_by(segment) |>
  mutate(day_no = row_number()) |>
  ungroup()

fit_gridpoint <- function(df){
  lon_val     <- df$lon[1]
  lat_val     <- df$lat[1]
  i_local_val <- df$i_local[1]
  j_local_val <- df$j_local[1]

  d <- df[, c("resid","segment","day_no","wrname",predictors), drop=FALSE]
  d <- d |> filter(if_all(everything(), is.finite))
  if(nrow(d)<80 || sd(d$resid)==0 || length(unique(d$wrname))<2) return(NULL)
  d$wrname <- factor(d$wrname)
  log_mem("fit wr model")

  fit_wr <- fit_gls(resid ~ wrname, d)
  log_mem("fit full model")

  fit_full <- fit_gls(as.formula(paste0("resid ~ wrname + ", paste(predictors, collapse=" + "))), d)
  log_mem("emm means")

  emm_wr <- as.data.frame(summary(emmeans(fit_wr, ~wrname, data=d, mode="asymptotic"), infer=c(TRUE,TRUE)))
  emm_full <- as.data.frame(summary(emmeans(fit_full, ~wrname, data=d, mode="asymptotic"), infer=c(TRUE,TRUE)))
  log_mem("joint metrics")

  joint_metrics <- pattern_metrics(emm_wr, emm_full)

  pred_att_list <- lapply(predictors, function(pred){
    fit_pred <- fit_gls(as.formula(paste0("resid ~ wrname + ", pred)), d)
    emm_pred <- as.data.frame(summary(emmeans(fit_pred, ~wrname, data=d, mode="asymptotic"), infer=c(TRUE,TRUE)))
    met <- pattern_metrics(emm_wr, emm_pred)
    data.frame(
      predictor=pred,
      n=nrow(d),
      wr_pattern_var_reduction=met$var_reduction,
      wr_pattern_mad_reduction=met$mad_reduction,
      wr_only_mean_abs_dev=met$mad0,
      pred_mean_abs_dev=met$mad1,
      row.names=NULL
    )
  })
  pred_att <- bind_rows(pred_att_list)

  list(
    summary = data.frame(
      n=nrow(d), n_wr=nlevels(d$wrname),
      joint_wr_pattern_var_reduction=joint_metrics$var_reduction,
      joint_wr_pattern_mad_reduction=joint_metrics$mad_reduction,
      lon=lon_val, lat=lat_val,
      i_local=i_local_val, j_local=j_local_val,
      row.names=NULL
    ),
    wr_only_emmeans = transform(emm_wr,
      lon=lon_val, lat=lat_val,
      i_local=i_local_val, j_local=j_local_val),
    full_emmeans = transform(emm_full,
      lon=lon_val, lat=lat_val,
      i_local=i_local_val, j_local=j_local_val),
    predictor_attenuation = transform(pred_att,
      lon=lon_val, lat=lat_val,
      i_local=i_local_val, j_local=j_local_val)
  )
}

# --- incremental update for predictor_attenuation only ---
existing_out <- if(file.exists(outfile) && !OVERWRITE){
  message("Existing output found, loading for incremental update: ", outfile)
  readRDS(outfile)
} else {
  list(summary=data.frame(), wr_only_emmeans=data.frame(),
       full_emmeans=data.frame(), predictor_attenuation=data.frame())
}

nx <- length(chunk$lon)
ny <- length(chunk$lat)

summary_out <- list(); wr_only_emm_out <- list()
full_emm_out <- list(); pred_att_out <- list()
k1 <- k2 <- k3 <- k4 <- 1L

done_keys <- unique(paste(existing_out$predictor_attenuation$i_local,
                          existing_out$predictor_attenuation$j_local, sep="::"))

message("Done keys: ", paste(done_keys, collapse=", "))

for(ii in seq_len(nx)){
  for(jj in seq_len(ny)){
    key <- paste(ii,jj,sep="::")
    if(!OVERWRITE && key %in% done_keys) next

    df       <- base_df
    df$lon   <- chunk$lon[ii]
    df$lat   <- chunk$lat[jj]
    df$i_local <- ii
    df$j_local <- jj
    df$resid <- chunk$residuals[ii, jj, df$time_idx]   # <-- sliced
    for (pred in predictors) {
      df[[pred]] <- chunk[[pred]][ii, jj, df$time_idx] # <-- sliced
    }
    log_mem(paste0("now fitting gridpoint (i,j) = ",  ii, jj))

    ans <- tryCatch(fit_gridpoint(df), error=function(e) NULL)
    if(is.null(ans)) next

    summary_out[[k1]] <- ans$summary; k1 <- k1+1
    wr_only_emm_out[[k2]] <- ans$wr_only_emmeans; k2 <- k2+1
    full_emm_out[[k3]] <- ans$full_emmeans; k3 <- k3+1
    pred_att_out[[k4]] <- ans$predictor_attenuation; k4 <- k4+1
  }
}

out <- list(
  metadata=list(
    chunk_no=CHUNK_NO,
    chunk_file=chunk_file,
    source_chunk_metadata=chunk$metadata,
    wr_rds=WR_RDS,
    predictors=predictors,
    incremental_update=file.exists(outfile) && !OVERWRITE
  ),
  summary=bind_rows(existing_out$summary, bind_rows(summary_out)) |> distinct(i_local,j_local,.keep_all=TRUE),
  wr_only_emmeans=bind_rows(existing_out$wr_only_emmeans, bind_rows(wr_only_emm_out)) |> distinct(i_local,j_local,wrname,.keep_all=TRUE),
  full_emmeans=bind_rows(existing_out$full_emmeans, bind_rows(full_emm_out)) |> distinct(i_local,j_local,wrname,.keep_all=TRUE),
  predictor_attenuation=bind_rows(existing_out$predictor_attenuation, bind_rows(pred_att_out)) |> distinct(i_local,j_local,predictor,.keep_all=TRUE)
)

saveRDS(out, outfile, compress="xz")
message("Wrote: ", outfile)