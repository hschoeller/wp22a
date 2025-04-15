library(sf)
library(ggplot2)
library(mapproj)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(rnaturalearthdata)
library(latex2exp)
library(colorspace)
library(tikzDevice)
library(cowplot)

# Plot aesthetics
THEME_PUB <- theme_minimal() +
    theme(
        text = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank(),
        # axis.text = element_blank(),
        # axis.title = element_blank(),
        legend.position = "bottom",
        # legend.title = element_blank(),
        panel.background = element_rect(fill = "white", color = NA)
    )

THEME_PUB_LARGE <- theme_minimal(base_size = 16) +
    theme(
        text = element_text(family = "Helvetica", size = 16),
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "lines"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0)
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

get_diverging_scale <- function(clims = NULL) {
    if (is.null(clims)) {
        return(scale_fill_scico(palette = "vik", direction = 1, midpoint = 0))
    } else {
        # Ensure the color scale is properly centered at 0
        max_abs <- max(abs(clims))
        if (any(clims < 0) && any(clims > 0)) {
            # If clims spans negative and positive values
            return(scale_fill_scico(
                palette = "vik",
                direction = 1,
                midpoint = 0,
                limits = clims,
                na.value = "transparent"
            ))
        } else {
            # If clims is only positive or only negative
            return(scale_fill_scico(
                palette = "vik",
                direction = 1,
                limits = clims,
                na.value = "transparent"
            ))
        }
    }
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
                         clims = NULL, show_graticule_labels = TRUE,
                         colorbar_ticks = 2, use_diverging = FALSE) {
    # Create WGS84 grid from data
    grid_wgs84 <- st_bbox(c(
        xmin = min(data$lon),
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
        st_transform(CRS) # %>%
    # filter(!is.na(.data[[var_name]]))

    # Get the extent of the projected grid
    bbox_proj <- st_bbox(grid_proj)

    # Create raster from sf object
    # First ensure your data has the variable with proper name
    names(grid_proj)[names(grid_proj) == var_name] <- "value"

    # Convert to stars object (raster)
    resolution_factor <- 8
    grid_rast <- stars::st_rasterize(
        grid_proj["value"], # Use the column we need
        nx = length(unique(data$lon)) * resolution_factor,
        ny = length(unique(data$lat)) * resolution_factor,
        bounds = bbox_proj
    )

    grid_rast_df <- as.data.frame(grid_rast, xy = TRUE)

    # Now grid_rast_df will have columns x, y, and value.
    # You can then remove NA cells if you wish:
    grid_rast_df <- grid_rast_df %>% filter(!is.na(value))

    # Prepare map elements
    coastlines <- ne_countries(scale = "medium", returnclass = "sf") %>%
        st_transform(CRS)

    graticule <- st_graticule(
        lat = GRAT_LAT,
        lon = GRAT_LON,
        crs = 4326
    ) %>%
        st_transform(CRS)

    # Create plot
    p <- ggplot() +
        geom_raster(
            data = grid_rast_df, aes(x = x, y = y, fill = value),
            interpolate = TRUE
        ) +
        # ggplot() +
        # geom_sf(
        #     data = grid_proj,
        #     aes(fill = .data[[var_name]]),
        #     color = NA,
        #     linewidth = 0
        # )
        # stars::geom_stars(data = grid_rast, interpolate = TRUE) +
        # scale_fill_continuous(
        #     name = if (legend_name == "") TeX(var_name) else TeX(legend_name)
        # ) +
        geom_sf(data = graticule, color = "gray60", linewidth = 0.3) +
        geom_sf(data = coastlines, fill = NA, color = "black", linewidth = 0.3) +
        coord_sf(
            xlim = c(bbox_proj["xmin"], bbox_proj["xmax"]),
            ylim = c(bbox_proj["ymin"], bbox_proj["ymax"]),
            expand = FALSE
        )

    if (use_diverging) {
        p <- p + get_diverging_scale(clims = clims)
    } else {
        p <- p + get_continuous_scale(clims = clims)
    }
    # Option to show or hide graticule labels
    if (!show_graticule_labels) {
        p <- p +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank()
            )
    } else {
        p <- p +
            scale_x_continuous(breaks = GRAT_LON) +
            scale_y_continuous(breaks = GRAT_LAT) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }

    p <- p + labs(
        fill = if (legend_name == "") TeX(var_name) else TeX(legend_name),
        title = title,
        x = "",
        y = ""
    ) #+
    # theme(
    #     axis.title.x = element_blank(),
    #     axis.title.y = element_blank(),

    # legend.position = "bottom",
    # )
    if (!is.null(clims)) {
        p <- p +
            guides(fill = guide_colorbar(
                barwidth = 10,
                barheight = 1,
                ticks.linewidth = 1,
                # Set exactly colorbar_ticks number of breaks on the colorbar
                breaks = seq(clims[1], clims[2], length.out = colorbar_ticks)
            ))
    } else {
        # If no clims provided, use the range of the data
        value_range <- range(grid_rast_df$value, na.rm = TRUE)
        # print(seq(value_range[1], value_range[2],
        #         length.out = colorbar_ticks))
        p <- p +
            guides(fill = guide_colorbar(
                barwidth = 10,
                barheight = 1,
                ticks.linewidth = 1,
                breaks = seq(value_range[1], value_range[2],
                    length.out = colorbar_ticks
                )
            ))
    }

    return(p)
}

add_contour <- function(
    plot_obj, data, contour_var, alpha = ALPHA, CRS,
    resolution_factor = 8,
    contour_color = "red", contour_alpha = 0.8,
    contour_linewidth = 0.5) {
    # Extract coordinate limits from original plot
    plot_build <- ggplot_build(plot_obj)
    xlim <- plot_build$layout$panel_params[[1]]$x.range
    ylim <- plot_build$layout$panel_params[[1]]$y.range

    # Create WGS84 grid for contour data
    grid_wgs84_contour <- st_bbox(c(
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

    # Transform to target CRS using explicit CRS parameter
    grid_contour_proj <- grid_wgs84_contour %>%
        st_transform(crs = CRS) %>% # Explicit CRS specification
        filter(!is.na(.data[[contour_var]]))

    # Rasterize contour data
    nx <- length(unique(data$lon)) * resolution_factor
    ny <- length(unique(data$lat)) * resolution_factor

    grid_contour_rast <- stars::st_rasterize(
        grid_contour_proj[contour_var],
        nx = nx,
        ny = ny,
        bounds = c(xlim[1], ylim[1], xlim[2], ylim[2])
    )

    # Convert to data frame
    grid_contour_df <- as.data.frame(grid_contour_rast, xy = TRUE) %>%
        filter(!is.na(.data[[contour_var]]))

    # Add contour layer to existing plot
    plot_obj +
        geom_contour(
            data = grid_contour_df,
            aes(x = x, y = y, z = .data[[contour_var]]),
            color = contour_color,
            linewidth = contour_linewidth,
            alpha = contour_alpha,
            binwidth = 750
        )
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

save_plot <- function(
    plot_obj, filename, width = 6, height = 4,
    type = "pdf") {
    if (type == "tikz") {
        plot_obj <- plot_obj + theme_minimal(base_size = 33)
        tikz(
            file = file.path(OUT_DIR, filename),
            width = width,
            height = height,
            standAlone = TRUE
        )
    } else {
        cairo_pdf(
            file = file.path(OUT_DIR, filename),
            width = width,
            height = height,
            # family = "Helvetica",
            fallback_resolution = 300
        )
    }
    print(plot_obj)
    dev.off()
}
