source("RScripts/config.r")
source("RScripts/data_functions.r")

lm_files <- list.files(LM_DIR, pattern = "^lm.*\\.Rds$", full.names = TRUE)

lm_info <- tibble(file = lm_files) %>%
    mutate(
        coords = str_extract_all(basename(file),
            "-?\\d+\\.?\\d*",
            simplify = TRUE
        ),
        lat = as.numeric(coords[, 1]),
        lon = as.numeric(coords[, 2])
    ) %>%
    select(-coords)

extract_residuals <- function(info_row) {
    file <- info_row$file
    lat <- info_row$lat
    lon <- info_row$lon
    m <- readRDS(file)
    d <- m$model %>%
        arrange(year) %>%
        mutate(date = as.Date(paste0(min(year), "-01-01")) + row_number() - 1)
    tibble(
        lat = lat,
        lon = lon,
        date = d$date,
        residual = m$residuals[order(d$year)]
    )
}

print("files screened")
num_cores <- detectCores(logical = FALSE) - 1
print(paste0("using ", num_cores, "cores"))
residuals_list <- mclapply(
    split(lm_info, seq_len(nrow(lm_info))),
    extract_residuals,
    mc.cores = num_cores
)

print("residuals loaded")
all_residuals <- bind_rows(residuals_list)
print("saving residuals")
saveRDS(all_residuals, "ens_data/LM_res.Rds")
