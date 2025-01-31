library(dplyr)
library(purrr)
library(tidyr)
library(broom)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(stringr)

# Plot aesthetics
THEME_PUB <- theme_minimal() +
    theme(
        text = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
    )

COLOR_SCALE <- scale_fill_viridis_c(option = "plasma", na.value = "grey90")


load_lm_objects <- function(lm_dir = LM_DIR) {
    lm_files <- list.files(lm_dir, pattern = "^lm.*\\.Rds", full.names = TRUE)

    # Extract lat/lon from filenames
    coords <- stringr::str_extract_all(lm_files, "-?\\d+\\.?\\d*",
        simplify = TRUE
    )
    lats <- as.numeric(coords[, 1])
    lons <- as.numeric(coords[, 2])

    # Load models and extract diagnostics
    purrr::map_dfr(lm_files, ~ {
        lm_obj <- readRDS(.x)

        # Model-level diagnostics (broom::glance)
        model_summary <- broom::glance(lm_obj)

        # Coefficients and p-values (broom::tidy)
        coef_table <- broom::tidy(lm_obj)

        tibble::tibble(
            lat = lats[which(lm_files == .x)],
            lon = lons[which(lm_files == .x)],
            r_squared = model_summary$r.squared,
            adj_r_squared = model_summary$adj.r.squared,
            sigma = model_summary$sigma, # Residual standard error
            aic = model_summary$AIC,
            bic = model_summary$BIC,
            f_statistic = model_summary$statistic, # Overall F-statistic
            f_p_value = model_summary$p.value,
            coefficients = list(coef_table$estimate),
            p_values = list(coef_table$p.value)
        )
    })
}

plot_spatial <- function(data, var_name,
                         title = "",
                         geo_path = file.path(GEO_DIR, "coastlines.shp")) {
    # Load geographic data
    coastlines <- sf::st_read(geo_path, quiet = TRUE)

    # Create plot
    p <- ggplot2::ggplot(data) +
        ggplot2::geom_tile(aes(x = lon, y = lat, fill = .data[[var_name]])) +
        ggplot2::geom_sf(
            data = coastlines, fill = NA,
            color = "black", linewidth = 0.2
        ) +
        THEME_PUB +
        COLOR_SCALE +
        labs(title = title)

    return(p)
}

plot_diagnostics <- function(data, var_name, title,
                             geo_path = file.path(GEO_DIR, "coastlines.shp")) {
    coastlines <- sf::st_read(geo_path, quiet = TRUE)

    ggplot2::ggplot(data) +
        ggplot2::geom_tile(aes(x = lon, y = lat, fill = .data[[var_name]])) +
        ggplot2::geom_sf(
            data = coastlines, fill = NA, color = "black",
            linewidth = 0.2
        ) +
        THEME_PUB +
        scale_fill_viridis_c(
            option = ifelse(var_name == "r_squared", "viridis", "plasma"),
            na.value = "grey90"
        ) +
        labs(title = title, fill = "")
}

plot_model_significance <- function(data) {
    coastlines <- sf::st_read(file.path(GEO_DIR, "coastlines.shp"),
        quiet = TRUE
    )

    data %>%
        mutate(
            significant = ifelse(f_p_value < 0.05,
                "Significant", "Not significant"
            )
        ) %>%
        ggplot() +
        geom_tile(aes(x = lon, y = lat, fill = significant)) +
        geom_sf(
            data = coastlines, fill = NA, color = "black",
            linewidth = 0.2
        ) +
        THEME_PUB +
        scale_fill_manual(values = c(
            "Significant" = "#440154",
            "Not significant" = "grey90"
        )) +
        labs(title = "Overall Model Significance (F-test, p < 0.05)")
}

save_publication_plot <- function(plot_obj, filename, width = 6, height = 4) {
    cairo_pdf(
        file = file.path(OUT_DIR, filename),
        width = width,
        height = height,
        family = "Helvetica"
    )
    print(plot_obj)
    dev.off()
}
