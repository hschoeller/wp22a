source("../RScripts/config.r")
# source("../RScripts/plot_functions.r")
source("../RScripts/data_functions.r")
# source("../RScripts/algo_functions.r")

result <- wrera(
    start = "19500111_00",
    end = "20250113_21",
    hours = c("00", "06", "12", "18"),
    tformat = "string",
    setup = "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019",
    dataset = "era5",
    basepath = "../WR_read_example_package/wr_era5_update_1950_latwgt/"
)

wr_df <- result$data$LC %>%
    filter(grepl("12$", time))
wr_df$date <- as.Date(wr_df$time, format = "%Y%m%d_%H")

print("WRs loaded")

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(R_THREADS = "1")

composite_df <- calculate_composite_values(wr_df, LM_DIR,
    n_perm = 10000, n_cores = parallel::detectCores()
)

print("Composites calculated")

saveRDS(composite_df, file = "../ens_data/composites_permuted.RDS")

print("Composites saved")
