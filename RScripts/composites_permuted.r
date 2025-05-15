source("RScripts/config.r")
source("RScripts/data_functions.r")

result <- wrera(
    start = "19500111_00",
    end = "20250113_21",
    hours = c("12"),
    tformat = "string",
    setup = "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019",
    dataset = "era5",
    basepath = "WR_read_example_package/wr_era5_update_1950_latwgt/"
)
wr_df <- result$data$LC
wr_df$date <- as.Date(wr_df$time)

print("WRs loaded")
cps <- as.Date(paste0(CP, "-01"), format = "%Y-%m-%d")

wr_df$segment <- cut(as.Date(wr_df$date),
    breaks = c(as.Date(c(
        wr_df$date[1],
        cps,
        wr_df$date[nrow(wr_df)]
    ))),
    labels = 0:(length(cps)) + 1,
    include.lowest = TRUE,
    right = FALSE
)

print("CPs loaded")

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(R_THREADS = "1")

for (seg in c(6)) {
    wr_df_i <- wr_df %>%
        filter(segment == seg)
    comp_df <- calculate_composite_values(
        df = wr_df_i,
        lm_dir = LM_DIR,
        n_cores = parallel::detectCores() - 1
    )
    saveRDS(comp_df, file = paste0(
        "ens_data/composites_permutedSeg",
        seg, ".RDS"
    ))
}

# composite_df <- calculate_composite_values(wr_df, LM_DIR,
#     n_perm = 10000, n_cores = parallel::detectCores()
# )

print("Composites calculated")

# saveRDS(composite_df, file = "../ens_data/composites_permuted.RDS")

# print("Composites saved")
