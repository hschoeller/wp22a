library(sf)
library(ggplot2)
library(mapproj)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(rnaturalearthdata)
library(latex2exp)
library(colorspace)
# library(tikzDevice)
library(cowplot)
library(stars)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

THEME_PUB <- theme_minimal() +
    theme(
        text = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank(),
        # axis.text = element_blank(),
        # axis.title = element_blank(),
        legend.position = "bottom",
        # legend.title = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "white", color = NA),
        axis.title.y = element_text(size = 10),
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
        return(scale_fill_scico(
            palette = CONT_SEQ_SCALE,
            limits = clims,
            expand = expansion(mult = 0, add = 0)
        ))
    }
}

get_categorical_scale <- function(clims = NULL) {
    return(scale_color_brewer(
        aesthetics = c("color", "fill"),
        palette = "Dark2"
    ))
    # return(scale_color_scico_d(
    #     aesthetics = c("color", "fill"),
    #     palette = CONT_SEQ_SCALE
    # ))
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
    # First, create a complete grid
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
        mutate(id = row_number())

    # Convert data points to SF
    data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

    # Join data to grid
    grid_with_data <- st_join(grid_wgs84, data_sf)

    # Aggregate values per grid cell
    if (is.null(sig_name)) {
        grid_values <- grid_with_data %>%
            group_by(id) %>%
            summarize(
                !!var_name := mean(.data[[var_name]], na.rm = TRUE)
            )
    } else {
        grid_values <- grid_with_data %>%
            group_by(id) %>%
            summarize(
                !!var_name := mean(.data[[var_name]], na.rm = TRUE),
                !!sig_name := if (!is.null(sig_name)) {
                    mean(.data[[sig_name]],
                        na.rm = TRUE
                    )
                } else {
                    NULL
                }
            )
    }
    # Important fix: remove sf class from grid_values before joining
    grid_values <- grid_values %>% st_drop_geometry()

    # Merge back with the complete grid
    grid_complete <- grid_wgs84 %>%
        left_join(grid_values, by = "id")

    # Fill in NA values using nearest non-NA neighbor
    # First, identify cells with data
    has_data <- !is.na(grid_complete[[var_name]])
    grid_with_data <- grid_complete[has_data, ]

    # For cells without data, find nearest cell with data
    grid_no_data <- grid_complete[!has_data, ]
    if (nrow(grid_no_data) > 0) {
        nearest_idx <- st_nearest_feature(grid_no_data, grid_with_data)
        grid_no_data[[var_name]] <- grid_with_data[[var_name]][nearest_idx]

        # Combine both sets
        grid_complete <- rbind(grid_with_data, grid_no_data)
    }

    # Apply significance filter if needed
    if (!is.null(sig_name)) {
        grid_complete <- grid_complete %>%
            mutate(!!var_name := ifelse(.data[[sig_name]] < alpha,
                .data[[var_name]], NA
            ))
    }

    # Transform grid to target projection
    grid_proj <- grid_complete %>%
        st_transform(CRS) # %>%
    # filter(!is.na(.data[[var_name]]))

    # Get the extent of the projected grid
    bbox_proj <- st_bbox(grid_proj)

    # Create raster from sf object
    # First ensure your data has the variable with proper name
    names(grid_proj)[names(grid_proj) == var_name] <- "value"

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
        geom_sf(
            data = graticule,
            color = "gray60", linewidth = 0.2, alpha = .5
        ) +
        geom_sf(
            data = coastlines, fill = NA,
            color = "grey30", linewidth = 0.15, alpha = .5
        ) +
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
    )
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
        p <- p +
            guides(fill = guide_colorbar(
                barwidth = 8,
                barheight = .8,
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
    resolution_factor = 8, contour_binwidth = 750,
    contour_color = "red", contour_alpha = 0.8,
    contour_linewidth = 0.5, contour_linetype = "solid") {
    # Extract coordinate limits from original plot
    plot_build <- ggplot_build(plot_obj)
    xlim <- plot_build$layout$panel_params[[1]]$x.range
    ylim <- plot_build$layout$panel_params[[1]]$y.range

    # Create WGS84 grid for contour data
    grid_wgs84_contour <- st_bbox(c(
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
            binwidth = 750,
            linetype = contour_linetype
        )
}


#--- Function to plot the data and change points ------------------------------
plot_change_points <- function(data,
                               cp_df,
                               theme_pub = NULL,
                               show_ci = TRUE,
                               title = "North Atlantic monthly mean g500 ensemble variance") {
    p <- ggplot(data = data, aes(x = date, y = log_variance)) +
        geom_line(aes(size = "Assimilated"), linewidth = .2) +
        labs(
            x = "Year",
            y = TeX("$\\log(\\sigma_{EDA})$"),
            title = title,
            color = "Break Points"
        ) +
        scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
        get_categorical_scale()

    if (!is.null(theme_pub)) {
        p <- p + theme_pub +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else {
        p <- p + THEME_PUB +

            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }

    if (nrow(cp_df) > 0) {
        # Compute a y position for the labels
        cp_df$y_pos <- max(data$log_variance, na.rm = TRUE) * 0.8
        cp_df$y_pos[cp_df$cp_date == min(cp_df$cp_date)] <- max(
            data$log_variance,
            na.rm = TRUE
        ) * 0.5

        # If showing uncertainty, add a shaded rectangle using cp_date_lower and cp_date_upper
        if (show_ci &&
            all(c("cp_date_lower", "cp_date_upper") %in% names(cp_df))) {
            p <- p +
                geom_rect(
                    data = cp_df,
                    mapping = aes(
                        xmin = cp_date_lower, xmax = cp_date_upper,
                        fill = factor(cp_no)
                    ),
                    ymin = -Inf, ymax = Inf,
                    alpha = 0.25,
                    inherit.aes = FALSE,
                    show.legend = FALSE
                )
        }

        # Add vertical dashed lines for the change points
        p <- p +
            geom_vline(
                data = cp_df,
                mapping = aes(xintercept = cp_date, color = factor(cp_no)),
                linetype = "dashed",
                linewidth = .5
            ) +
            geom_text(
                data = cp_df,
                mapping = aes(
                    x = cp_date, y = y_pos,
                    label = format(cp_date, "%Y-%m"),
                    color = factor(cp_no)
                ),
                angle = 90,
                vjust = -0.5,
                hjust = 0,
                size = 3,
                show.legend = FALSE
            ) +
            guides(color = guide_legend(
                override.aes = list(
                    shape = NA,
                    linetype = "dashed"
                )
            )) +
            theme(
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.key.size = unit(1, "lines")
            ) +
            get_categorical_scale()
    }
    return(p)
}

add_fitted_line_ci <- function(model, data, line_color = "red",
                               fill_color = "red", fill_alpha = 0.25) {
    if (inherits(model, "gls")) {
        fit_vals <- predict(model, newdata = data)
        X_new <- model.matrix(formula(model), data = data)
        V_beta <- vcov(model)
        SE_fit <- sqrt(rowSums((X_new %*% V_beta) * X_new))
        df_resid <- model$dims$N - model$dims$p
        t_crit <- qt(0.975, df = df_resid)
        lower_CI <- fit_vals - t_crit * SE_fit
        upper_CI <- fit_vals + t_crit * SE_fit
        preds <- data.frame(
            date = data$date,
            fit  = fit_vals,
            lwr  = lower_CI,
            upr  = upper_CI
        )
    } else if (inherits(model, "lm")) {
        preds <- as.data.frame(predict(model,
            newdata = data,
            interval = "confidence"
        ))
        preds$date <- data$date
    }
    # Create the fitted line and confidence ribbon layers
    list(
        geom_line(
            data = preds, aes(x = date, y = fit, size = "Fitted"),
            color = line_color, linewidth = .2, alpha = .8
        ),
        geom_ribbon(
            data = preds, aes(x = date, ymin = lwr, ymax = upr),
            fill = fill_color, alpha = fill_alpha, inherit.aes = FALSE
        ),
        scale_size_manual(
            name = "",
            values = c("Assimilated" = 0.25, "Fitted" = 0.25),
            guide = guide_legend(override.aes = list(color = c(
                "black",
                line_color
            )))
        )
    )
}


#--- Plot Composites ----------------------------------------------------------

grid_and_legend <- function(
    comp_df, var = "mean", sig_name = "p_value_adj", alpha = .05,
    clims = c(-.357, 0.357), legend_name = "$\\log(\\sigma^{2}_{EDA})$") {
    plots <- list()
    i <- 1
    for (wri in c(1, 6, 7, 2, 4, 5, 3, 0)) {
        wrn <- comp_df$wrname[comp_df$wr == wri][1]
        wr_df_i <- comp_df %>%
            filter(wr == wri)
        # Plotting the mean
        p <- plot_spatial(wr_df_i, var,
            legend_name = legend_name,
            sig_name = sig_name,
            alpha = alpha,
            show_graticule_labels = FALSE,
            use_diverging = TRUE,
            clims = clims
        ) + THEME_PUB_LARGE +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                plot.margin = margin(0, 0, 0, 0),
                panel.spacing = unit(0, "lines")
            ) + labs(title = wrn)
        plots[[i]] <- add_contour(
            p, wr_df_i,
            contour_var = "z",
            contour_color = "darkgreen",
            contour_linewidth = 0.4,
            CRS = CRS,
            contour_linetype = "solid",
            resolution_factor = 16,
            contour_binwidth = 9807 # every 100 gpdm
        )
        i <- i + 1
    }

    # Get legend from one plot with legend
    p_legend <- plot_spatial(wr_df_i, var,
        legend_name = legend_name,
        sig_name = sig_name,
        alpha = alpha,
        show_graticule_labels = FALSE,
        use_diverging = TRUE,
        clims = clims
    ) + theme(legend.position = "bottom")
    # Extract just the legend grob
    legend_grob <- ggplotGrob(
        p_legend
    )$grobs[[which(sapply(
        ggplotGrob(p_legend)$grobs,
        function(x) x$name
    ) == "guide-box")]]
    legend_grob$grobs[[1]]$grobs[[2]]$children[[1]]$gp$fontsize <- 24 # Title size
    legend_grob$grobs[[1]]$grobs[[1]]$children[[1]]$grobs[[1]]$children[[1]]$gp$fontsize <- 20 # Text size

    # Manually create a new plot containing just the legend
    legend_plot <- ggplot() +
        annotation_custom(legend_grob) +
        theme_void()
    return(list(plots, legend_plot))
}
combine_plots <- function(plots, legend_plot, orientation = "horizontal") {
    # Combine plots
    if (orientation == "horizontal") {
        combined <- plot_grid(plotlist = plots, ncol = 4)
        final <- plot_grid(combined, legend_plot, ncol = 1, rel_heights = c(1, 0.1), greedy = TRUE)
    } else {
        combined <- plot_grid(plotlist = plots, ncol = 2)
        final <- plot_grid(combined, legend_plot, ncol = 1, rel_heights = c(1, 0.1), greedy = TRUE)
    }
    # Print to see if it works
    return(final)
}


# ─────────────────────────────────────────────────────────────────────────────
# 3. Plot helper functions (simplified)
# ─────────────────────────────────────────────────────────────────────────────
make_grid_data <- function(R_max) {
    levels <- pretty(c(0, R_max), 4)[2:3]

    circles <- purrr::map_dfr(levels, ~ {
        theta <- seq(0, 2 * pi, length.out = 200)
        tibble::tibble(x = .x * cos(theta), y = .x * sin(theta), level = .x)
    })

    month_angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = 13)[1:12]
    month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    lines <- tibble::tibble(
        angle = month_angles,
        xend = R_max * cos(angle),
        yend = R_max * sin(angle),
        label_x = R_max * 1.1 * cos(angle),
        label_y = R_max * 1.1 * sin(angle),
        label = month_names
    )

    radius_labels <- tibble::tibble(
        x = -levels,
        y = 0,
        label = as.character(round(levels, 2))
    )

    list(circles = circles, lines = lines, radius_labels = radius_labels)
}


make_arcs <- function(df_ci) {
    purrr::map_dfr(1:nrow(df_ci), ~ {
        row <- df_ci[.x, ]
        angle_seq <- if (row$phi_upr < row$phi_lwr) {
            seq(row$phi_lwr, row$phi_upr + 2 * pi, length.out = 80)
        } else {
            seq(row$phi_lwr, row$phi_upr, length.out = 80)
        }
        tibble::tibble(
            segment_clean = row$segment_clean,
            x = row$R * cos(angle_seq),
            y = row$R * sin(angle_seq),
            group = .x
        )
    })
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. Main plotting function (simplified)
# ─────────────────────────────────────────────────────────────────────────────
plot_seasonal_polar <- function(df_ci_obj, marker = NULL, p = NULL) {
    df_ci <- if (is.null(marker)) {
        df_ci_obj
    } else {
        df_ci_obj %>% dplyr::filter(marker_shape == marker)
    }

    if (nrow(df_ci) == 0) stop("No rows to plot (after filtering by marker).")

    # base plot only once
    if (is.null(p)) {
        R_max <- max(df_ci_obj$R_upr, na.rm = TRUE) # keep scaling stable across additions
        grid_data <- make_grid_data(R_max)

        p <- ggplot2::ggplot() +
            ggplot2::geom_path(
                data = grid_data$circles,
                ggplot2::aes(x, y, group = level),
                color = "grey10", linewidth = 0.5
            ) +
            ggplot2::geom_segment(
                data = grid_data$lines,
                ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
                color = "grey10", linewidth = 0.3
            ) +
            ggplot2::geom_text(
                data = grid_data$lines,
                ggplot2::aes(label_x, label_y, label = label),
                size = 3.5, color = "grey10"
            ) +
            ggplot2::geom_text(
                data = grid_data$radius_labels,
                ggplot2::aes(x, y, label = label),
                size = 3, color = "grey10", hjust = 1.2, vjust = -0.5
            ) +
            ggplot2::scale_color_brewer(type = "qual", palette = "Dark2", name = "Segment") +
            ggplot2::coord_equal() +
            ggplot2::theme_void(base_size = 14) +
            ggplot2::theme(legend.position = "right") +
            ggplot2::labs(color = "Segment")
    }

    arcs <- make_arcs(df_ci)

    p <- p +
        ggplot2::geom_path(
            data = arcs,
            ggplot2::aes(x, y, group = group, color = segment_clean),
            linewidth = 0.7, alpha = 0.9
        ) +
        ggplot2::geom_segment(
            data = df_ci,
            ggplot2::aes(
                x = R_lwr * cos(phi), y = R_lwr * sin(phi),
                xend = R_upr * cos(phi), yend = R_upr * sin(phi),
                color = segment_clean
            ),
            linewidth = 0.7, alpha = 0.9
        )

    # points: plain if marker is NULL, otherwise triangles (or whatever integer)
    if (is.null(marker)) {
        p <- p + ggplot2::geom_point(
            data = df_ci,
            ggplot2::aes(x, y, fill = segment_clean),
            shape = 16, size = 2.5
        )
    } else {
        p <- p + ggplot2::geom_point(
            data = df_ci,
            ggplot2::aes(x, y, fill = segment_clean, color = segment_clean),
            shape = marker, size = 3, stroke = 0.8
        )
    }

    p <- p +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Dark2", name = "Segment")

    p
}

plot_multiple_grid_points_daily_one_year <- function(df_obj, year) {
    # helper: darken hex by mixing toward black
    darken_hex <- function(hex, amount = 0.50) {
        rgb <- grDevices::col2rgb(hex) / 255
        rgb2 <- rgb * (1 - amount)
        grDevices::rgb(rgb2[1, ], rgb2[2, ], rgb2[3, ])
    }

    long <- df_obj %>%
        dplyr::mutate(grid = paste(lat, lon, sep = "_")) %>%
        tidyr::unnest(obj) %>%
        tidyr::pivot_longer(
            cols = c(obs_log, fit_log, lwr_log, upr_log),
            names_to = "type",
            values_to = "value"
        ) %>%
        dplyr::filter(
            date >= as.Date(sprintf("%d-01-01", year)),
            date <= as.Date(sprintf("%d-12-31", year))
        )

    grids <- unique(long$grid)
    n_grids <- length(grids)

    # Okabe-Ito palette (up to 8; extend if needed)
    base_palette <- c(
        "#E69F00", "#56B4E9", "#009E73", "#999999",
        "#CC79A7", "#F0E442", "#0072B2", "#D55E00"
    )
    if (n_grids > length(base_palette)) {
        base_palette <- grDevices::colorRampPalette(base_palette)(n_grids)
    } else {
        base_palette <- base_palette[seq_len(n_grids)]
    }

    pal_obs <- setNames(base_palette, grids)
    pal_ci <- vapply(pal_obs, darken_hex, character(1), amount = 0.50)

    ribbon_df <- long %>%
        dplyr::filter(type %in% c("lwr_log", "upr_log")) %>%
        dplyr::select(date, grid, type, value) %>%
        tidyr::pivot_wider(names_from = type, values_from = value)

    obs_pts <- long %>%
        dplyr::filter(type == "obs_log")

    ggplot() +
        # CI ribbon behind
        geom_ribbon(
            data = ribbon_df,
            aes(x = date, ymin = lwr_log, ymax = upr_log, fill = grid, group = grid),
            alpha = 0.95,
            show.legend = FALSE
        ) +
        # Observed POINTS (shape from marker_shape; colored by grid)
        geom_point(
            data = obs_pts,
            aes(x = date, y = value, color = grid, shape = marker_shape),
            size = 1.8,
            fill = "white",
            stroke = 0.8,
            alpha = 1,
            show.legend = FALSE
        ) +
        scale_color_manual(values = pal_obs, name = NULL) +
        scale_fill_manual(values = pal_ci) +
        scale_shape_identity() +
        labs(x = "Month", y = TeX("$\\log(\\sigma_{EDA})$")) +
        scale_x_date(
            limits = c(
                as.Date(sprintf("%d-01-01", year)),
                as.Date(sprintf("%d-12-31", year))
            ),
            breaks = seq(
                as.Date(sprintf("%d-01-15", year)),
                as.Date(sprintf("%d-12-15", year)),
                by = "1 month"
            ),
            labels = format(
                seq(
                    as.Date(sprintf("%d-01-15", year)),
                    as.Date(sprintf("%d-12-15", year)),
                    by = "1 month"
                ),
                "%b"
            ),
            expand = c(.01, 0)
        ) +
        THEME_PUB +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
}

plot_multiple_grid_lines_monthly_inset <- function(df_obj,
                                                   split_date = NULL,
                                                   late_ylim_quant = c(0.02, 0.98),
                                                   late_ylim_pad = 0.05,
                                                   inset_pos = c(left = 0.58, bottom = 0.55, right = 0.98, top = 0.98)) {
    # helper: darken hex by mixing toward black
    darken_hex <- function(hex, amount = 0.50) {
        rgb <- grDevices::col2rgb(hex) / 255
        rgb2 <- rgb * (1 - amount)
        grDevices::rgb(rgb2[1, ], rgb2[2, ], rgb2[3, ])
    }

    long <- df_obj %>%
        dplyr::mutate(grid = paste(lat, lon, sep = "_")) %>%
        tidyr::unnest(obj) %>%
        tidyr::pivot_longer(
            cols = c(obs_log, fit_log, lwr_log, upr_log),
            names_to = "type",
            values_to = "value"
        )

    # choose split_date if not provided
    if (is.null(split_date)) {
        split_date <- stats::median(long$date, na.rm = TRUE)
    } else {
        split_date <- as.Date(split_date)
    }
    start_pts <- long %>%
        dplyr::filter(type == "obs_log") %>%
        dplyr::group_by(grid) %>%
        dplyr::slice_min(date, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(date, value, marker_shape, grid)
    start_pts_inset <- long %>%
        dplyr::filter(type == "obs_log", date >= split_date + 0) %>%
        dplyr::group_by(grid) %>%
        dplyr::slice_min(date, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(date, value, marker_shape, grid)

    grids <- unique(long$grid)
    n_grids <- length(grids)

    # Okabe-Ito palette (up to 8; extend if needed)
    base_palette <- c(
        "#E69F00", "#56B4E9", "#009E73", "#999999",
        "#CC79A7", "#F0E442", "#0072B2", "#D55E00"
    )
    if (n_grids > length(base_palette)) {
        base_palette <- grDevices::colorRampPalette(base_palette)(n_grids)
    } else {
        base_palette <- base_palette[seq_len(n_grids)]
    }

    pal_obs <- setNames(base_palette, grids)
    pal_ci <- vapply(pal_obs, darken_hex, character(1), amount = 0.50)

    ribbon_df <- long %>%
        dplyr::filter(type %in% c("lwr_log", "upr_log")) %>%
        dplyr::select(date, grid, type, value) %>%
        tidyr::pivot_wider(names_from = type, values_from = value)

    # --- compute late-period ylim from obs + CI ---
    late_vals <- dplyr::bind_rows(
        long %>% dplyr::filter(type == "obs_log", date >= split_date) %>%
            dplyr::transmute(v = value),
        ribbon_df %>% dplyr::filter(date >= split_date) %>%
            dplyr::transmute(v = lwr_log),
        ribbon_df %>% dplyr::filter(date >= split_date) %>%
            dplyr::transmute(v = upr_log)
    ) %>%
        dplyr::filter(is.finite(v))

    q <- stats::quantile(late_vals$v, probs = late_ylim_quant, na.rm = TRUE, names = FALSE)
    span <- diff(q)
    if (!is.finite(span) || span == 0) span <- 1
    late_ylim <- c(q[1] - late_ylim_pad * span, q[2] + late_ylim_pad * span)

    x_min <- min(long$date, na.rm = TRUE)
    x_max <- max(long$date, na.rm = TRUE)

    zoom_box <- tibble::tibble(
        xmin = split_date, xmax = x_max,
        ymin = late_ylim[1], ymax = late_ylim[2]
    )

    # --- common layers (NO zoom box here) ---
    p_layers_main <- ggplot() +
        geom_ribbon(
            data = ribbon_df,
            aes(x = date, ymin = lwr_log, ymax = upr_log, fill = grid, group = grid),
            alpha = 1,
            show.legend = FALSE
        ) +
        geom_line(
            data = long %>% dplyr::filter(type == "obs_log"),
            aes(x = date, y = value, color = grid, group = grid),
            linewidth = 0.175,
            alpha = 0.75
        ) +
        geom_point(
            data = start_pts,
            aes(x = date, y = value, shape = marker_shape),
            inherit.aes = FALSE,
            size = 3,
            fill = "white",
            color = "black",
            stroke = 0.8,
            show.legend = FALSE
        ) +
        scale_color_manual(values = pal_obs, name = NULL) +
        scale_fill_manual(values = pal_ci) +
        scale_shape_identity() +
        labs(x = "Year", y = TeX("$\\log(\\sigma_{EDA})$")) +
        scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
        THEME_PUB +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none" # "bottom"
        )
    # --- inset layers (same as main, but no zoom box) ---
    p_layers_inset <- ggplot() +
        geom_ribbon(
            data = ribbon_df,
            aes(x = date, ymin = lwr_log, ymax = upr_log, fill = grid, group = grid),
            alpha = 1,
            show.legend = FALSE
        ) +
        geom_line(
            data = long %>% dplyr::filter(type == "obs_log"),
            aes(x = date, y = value, color = grid, group = grid),
            linewidth = .25,
            alpha = 0.95
        ) +
        geom_point(
            data = start_pts_inset,
            aes(x = date, y = value, shape = marker_shape),
            inherit.aes = FALSE,
            size = 3,
            fill = "white",
            color = "black",
            stroke = 0.8,
            show.legend = FALSE
        ) +
        scale_color_manual(values = pal_obs, name = NULL) +
        scale_fill_manual(values = pal_ci) +
        scale_shape_identity() +
        labs(x = "Year", y = TeX("$\\log(\\sigma_{EDA})$")) +
        scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
        THEME_PUB +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom"
        )
    # --- main plot WITH zoom box ---
    p_main <- p_layers_main +
        geom_rect(
            data = zoom_box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA,
            color = "black",
            linewidth = 0.35
        ) +
        coord_cartesian(xlim = c(x_min, x_max)) # explicit; optional

    # --- inset plot (zoomed), with ONLY a panel border around the data ---
    p_inset <- p_layers_inset +
        scale_x_date(
            date_breaks = "5 year", date_labels = "%Y",
            expand = ggplot2::expansion(mult = 0, add = 0)
        ) +
        scale_y_continuous(expand = ggplot2::expansion(mult = 0, add = 0)) +
        coord_cartesian(xlim = c(split_date, x_max), ylim = late_ylim, expand = FALSE) +
        labs(x = NULL, y = NULL) +
        theme(
            legend.position = "none",
            axis.text = element_text(size = rel(0.7)),
            axis.title = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.background = element_blank()
        )

    p_main +
        patchwork::inset_element(
            p_inset,
            left   = inset_pos["left"],
            bottom = inset_pos["bottom"],
            right  = inset_pos["right"],
            top    = inset_pos["top"]
        )
}

plot_spectral_power_obs_and_resid <- function(df_obj,
                                              detrend = TRUE,
                                              taper = 0.1,
                                              smooth_span = NULL,
                                              min_period_days = 2) {
    base_palette <- c(
        "#E69F00", "#56B4E9", "#009E73", "#999999",
        "#CC79A7", "#F0E442", "#0072B2", "#D55E00"
    )

    grids <- paste(df_obj$lat, df_obj$lon, sep = "_")
    n_grids <- length(grids)

    if (n_grids > length(base_palette)) {
        base_palette <- grDevices::colorRampPalette(base_palette)(n_grids)
    } else {
        base_palette <- base_palette[seq_len(n_grids)]
    }

    pal <- setNames(base_palette, grids)

    spec_list <- lapply(seq_len(nrow(df_obj)), function(i) {
        g <- grids[i]
        shape_val <- df_obj$marker_shape[i]

        dat <- df_obj$obj[[i]] %>%
            dplyr::transmute(
                date = date,
                obs = obs_log,
                resid = obs_log - fit_log
            ) %>%
            dplyr::arrange(date)

        full_dates <- seq(min(dat$date), max(dat$date), by = "1 day")
        dat_full <- tibble::tibble(date = full_dates) %>%
            dplyr::left_join(dat, by = "date")

        mk_spec <- function(x) {
            x <- x[is.finite(x)]
            sp <- stats::spec.pgram(
                x,
                taper = taper,
                detrend = detrend,
                demean = TRUE,
                log = "no",
                plot = FALSE,
                spans = smooth_span
            )
            tibble::tibble(
                freq = sp$freq,
                period_days = 1 / sp$freq,
                power = sp$spec
            )
        }

        s_obs <- mk_spec(dat_full$obs) %>% dplyr::mutate(series = "Observations")
        s_resid <- mk_spec(dat_full$resid) %>% dplyr::mutate(series = "Residuals")

        dplyr::bind_rows(s_obs, s_resid) %>%
            dplyr::mutate(
                grid = g,
                marker_shape = shape_val
            )
    })

    spec_df <- dplyr::bind_rows(spec_list) %>%
        dplyr::filter(
            is.finite(period_days),
            is.finite(power),
            period_days >= min_period_days
        )

    # ---- Right-end marker positions ----
    end_pts <- spec_df %>%
        dplyr::group_by(grid, series) %>%
        dplyr::slice_max(period_days, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    # ---- Natural time unit breaks ----
    x_breaks <- c(1, 7, 30, 365, 3650)
    x_labels <- c("Day", "Week", "Month", "Year", "Decade")
    month_lines <- (1:12) * 30.4375 # avg days per month
    year_lines <- (1:10) * 365.25 # incl leap years approx
    vlines <- c(month_lines, year_lines)
    ggplot2::ggplot(
        spec_df,
        ggplot2::aes(
            x = period_days,
            y = power,
            color = grid,
            linetype = series
        )
    ) +
        ggplot2::geom_vline(
            xintercept = vlines,
            linewidth = 0.25,
            alpha = 0.1
        ) +
        ggplot2::geom_line(linewidth = 0.35, alpha = 0.95) +

        # Right-end markers
        ggplot2::geom_point(
            data = end_pts,
            ggplot2::aes(x = period_days, y = power, shape = marker_shape),
            inherit.aes = FALSE,
            size = 3,
            fill = "white",
            color = "black",
            stroke = 0.8,
            show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = pal, name = NULL) +
        ggplot2::scale_linetype_manual(
            values = c(
                Residuals = "solid",
                Observations = "dashed"
            ),
            name = NULL
        ) +
        ggplot2::scale_shape_identity() +
        ggplot2::scale_x_log10(
            breaks = x_breaks,
            labels = x_labels
        ) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            x = "Time scale",
            y = TeX("Power Spectral Density")
        ) +
        THEME_PUB +
        ggplot2::theme(
            legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1), # anchor legend's top-left corner
            legend.background = element_rect(fill = "white", color = NA),
            legend.key = element_blank(), axis.text.x = ggplot2::element_text(vjust = 0.5)
        ) +
        ggplot2::guides(color = "none")
}

save_plot <- function(
    plot_obj, filename, width = 6, height = 4,
    type = "pdf") {
    tight_theme <- theme(
        plot.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0, "pt")
    )

    plot_obj <- plot_obj + tight_theme
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
            fallback_resolution = 600
        )
    }
    print(plot_obj)
    dev.off()
}
