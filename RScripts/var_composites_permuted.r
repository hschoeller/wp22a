source("RScripts/config.r")
source("RScripts/data_functions.r")

result <- wrera(
    start = "19500111_00",
    end = "20250113_21",
    hours = c("00", "06", "12", "18"),
    tformat = "string",
    setup = "z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019",
    dataset = "era5",
    basepath = "WR_read_example_package/wr_era5_update_1950_latwgt/"
)

wr_df <- result$data$LC %>%
    filter(grepl("12$", time))
wr_df$date <- as.Date(wr_df$time, format = "%Y%m%d_%H")

print("WRs loaded")

nc <- nc_open("ens_data/data.nc")

# Identify the index for number == 0
num <- ncvar_get(nc, "number")
number_idx <- which(num == 0)
# Use the first occurrence (modify if needed)
number_idx <- number_idx[1]
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
time_data <- ncvar_get(nc, "valid_time")
time_origin <- sub(
    "seconds since ", "",
    ncatt_get(nc, "valid_time", "units")$value
)
time_whole <- as.Date(as.POSIXct(time_data,
    origin = time_origin,
    tz = "UTC"
))

start_vec <- c(1, 1, 1, number_idx)
count_vec <- c(-1, -1, -1, 1)
z_control <- ncvar_get(nc, "z", start = start_vec, count = count_vec)

nc_close(nc)


Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(R_THREADS = "1")



z_comp <- calculate_wr_composites(z_control, time_whole, lon, lat, wr_df,
    calculate_variance = TRUE, n_perm = 10000,
    n_cores = parallel::detectCores() - 1
)
z_comp$wr <- z_comp$wrindex
z_comp$p_value_adj <- p.adjust(z_comp$pval, method = "fdr")

saveRDS(comp_df, file = "ens_data/VarianceComposite.RDS")


print("Composites calculated")
