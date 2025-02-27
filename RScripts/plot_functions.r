library(sf)
library(ggplot2)
library(mapproj)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(rnaturalearthdata)
library(latex2exp)
library(colorspace)

# Plot aesthetics
THEME_PUB <- theme_minimal() +
    theme(
        text = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank(),
        # axis.text = element_blank(),
        # axis.title = element_blank(),
        legend.position = "bottom",
        # legend.title = element_blank()
    )

get_continuous_scale <- function(clims = NULL) {
    if (is.null(clims)) {
        return(scale_fill_scico(palette = CONT_SEQ_SCALE))
    } else {
        return(scale_fill_scico(palette = CONT_SEQ_SCALE, limits = clims))
    }
}

get_categorical_scale <- function(clims = NULL) {
    return(scale_color_scico_d(
        aesthetics = c("color", "fill"),
        palette = CONT_SEQ_SCALE
    ))
}

plot_heatmap <- function(heatmap_data,
                         total_cases = 8591,
                         text_size = 5,
                         title = "Number of grid points with significant coefficients by segment") {
    p <- ggplot(heatmap_data, aes(
        x = Metric,
        y = as.factor(Segment),
        fill = Count
    )) +
        geom_tile() +
        geom_text(
            aes(label = round(Count / total_cases, 2)),
            color = "black",
            size = text_size
        ) +
        scale_fill_scico(palette = "nuuk") +
        labs(
            title = title,
            x = "Metric",
            y = "Segment",
            fill = "Count"
        ) +
        THEME_PUB
    return(p)
}

plot_spatial <- function(data, var_name, legend_name = "",
                         sig_name = NULL, title = "", alpha = ALPHA,
                         clims = NULL) {
    # Create WGS84 grid from data
    grid_wgs84 <- st_bbox(c(
        xmin = min(data$lon) - 10,
        ymin = min(data$lat),
        xmax = max(data$lon),
        ymax = max(data$lat)
    ), crs = 4326) %>%
        st_make_grid(
            n = c(length(unique(data$lon)), length(unique(data$lat))),
            what = "polygons"
        ) %>%
        st_sf() %>%
        st_join(
            st_as_sf(data, coords = c("lon", "lat"), crs = 4326),
            join = st_contains
        )

    # Apply significance filter
    if (!is.null(sig_name)) {
        grid_wgs84 <- grid_wgs84 %>%
            mutate(!!var_name := ifelse(.data[[sig_name]] < alpha,
                .data[[var_name]], NA
            ))
    }

    # Transform grid to target projection
    grid_proj <- grid_wgs84 %>%
        st_transform(CRS) %>%
        filter(!is.na(.data[[var_name]]))

    # Prepare map elements
    coastlines <- ne_countries(scale = "medium", returnclass = "sf") %>%
        st_transform(CRS)

    graticule <- st_graticule(
        lat = GRAT_LAT,
        lon = GRAT_LON,
        crs = 4326
    ) %>%
        st_transform(CRS)

    # Calculate plot limits
    bbox_proj <- st_bbox(grid_proj)

    # Create plot
    p <- ggplot() +
        geom_sf(
            data = grid_proj,
            aes(fill = .data[[var_name]]),
            color = NA
        ) +
        geom_sf(data = graticule, color = "gray60", linewidth = 0.3) +
        geom_sf(data = coastlines, fill = NA, color = "black", linewidth = 0.3) +
        get_continuous_scale(clims = clims) +
        THEME_PUB +
        coord_sf(
            xlim = c(bbox_proj["xmin"], bbox_proj["xmax"]),
            ylim = c(bbox_proj["ymin"], bbox_proj["ymax"]),
            expand = FALSE
        ) +
        scale_x_continuous(
            breaks = GRAT_LON
        ) +
        scale_y_continuous(
            breaks = GRAT_LAT
        ) +
        labs(
            fill = if (legend_name == "") TeX(var_name) else TeX(legend_name),
            title = title
        )

    return(p)
}

#--- Function to plot the data and change points ------------------------------
plot_change_points <- function(data,
                               cp_df,
                               theme_pub = NULL) {
    p <- ggplot(data = data, aes(x = date, y = avg_z)) +
        geom_line() +
        scale_y_log10() +
        labs(
            x = "Time",
            y = "Log(Var)",
            title = "Original Data with Detected Change Points",
            color = "Hierarchy"
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
        get_categorical_scale()

    if (!is.null(theme_pub)) {
        p <- p + theme_pub
    } else {
        p <- p + THEME_PUB
    }

    if (nrow(cp_df) > 0) {
        y_pos <- max(data$avg_z, na.rm = TRUE) * 0.7
        p <- p +
            geom_vline(
                data = cp_df,
                aes(
                    xintercept = cp_date,
                    color = factor(max_cp)
                ),
                linetype = "dashed",
                linewidth = 1
            ) +
            geom_text(
                data = cp_df,
                aes(
                    x = cp_date, y = y_pos,
                    label = format(cp_date, "%Y-%m"),
                    color = factor(max_cp)
                ),
                angle = 90,
                vjust = -0.5,
                hjust = 0,
                size = 5,
                show.legend = FALSE
            ) +
            guides(color = guide_legend(
                override.aes = list(
                    shape = NA,
                    linetype = "dashed"
                )
            )) +
            theme(
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 18),
                legend.key.size = unit(2, "lines")
            )
    }
    return(p)
}

save_plot <- function(plot_obj, filename, width = 6, height = 4) {
    cairo_pdf(
        file = file.path(OUT_DIR, filename),
        width = width,
        height = height,
        family = "Helvetica"
    )
    print(plot_obj)
    dev.off()
}
