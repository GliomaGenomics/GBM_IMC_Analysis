clean_ia_data <- function(ia_df,
                          filter_by = c("none", "surgery", "region_type", "region_type_new"),
                          filter_val = NA,
                          count_method = c("classic", "patch"),
                          label_order = levels(lab_spe$manual_gating)) {
  filter_on <- match.arg(filter_by, several.ok = FALSE)

  if (length(filter_val) > 1) stop("Only one filter value can be selected")
  filter_on_vals <- if (is.na(filter_val)) "none" else filter_val

  count_filter <- match.arg(count_method, several.ok = FALSE)

  if (filter_on != "none") {
    plot_data <- ia_df %>%
      dplyr::filter(!!sym(filter_on) %in% filter_val) %>%
      dplyr::filter(count_method == count_filter)

    if (nrow(plot_data) == 0) {
      stop("No data found for the given filter value")
    }
  } else {
    plot_data <- ia_df %>%
      dplyr::filter(count_method == count_filter)
  }

  plot_data <- plot_data %>%
    dplyr::group_by(from_label, to_label) %>%
    dplyr::summarize(
      sum_sigval = sum(sigval, na.rm = TRUE),
      pct = abs(sum_sigval) / n() * 100,
      type = ifelse(sum_sigval == 0, NA, ifelse(sum_sigval > 0, "Interacting", "Avoiding")),
      .groups = "keep"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      filter_on = filter_on,
      filter_vals = paste0(filter_on_vals, collapse = ","),
      count_method = count_filter,
    ) %>%
    dplyr::mutate(
      across(type, ~ factor(., levels = c("Interacting", "Avoiding"))),
      across(pct, ~ ifelse(. == 0, NA, .)),
      across(c("from_label", "to_label"), ~ factor(., levels = label_order))
    )

  return(plot_data)
}


plot_ia <- function(ia_clean_df, highlight_colors = lab_spe@metadata$v2_colours$cell_groups){
    
    plot_subtitle <- glue::glue(
        "filters: {ia_clean_df$filter_on}",
        "filter values: {ia_clean_df$filter_vals}",
        "graph: delaunay triangulation",
        "count method: {unique(ia_clean_df$count_method)}",
        "p threshold: < 0.01",
        .sep = "\n"
    )
    
    p <- ia_clean_df %>%
        ggplot(
            aes(x = from_label, y = to_label)
        ) +
        geom_tile(fill = "#fffdfa", color = "grey50", linewidth = 0.1, linetype = 2) +
        geom_point(
            aes(size = pct, fill = type),
            color = "black", shape = 21, stroke = 0.5, na.rm = TRUE
        ) +
        scale_fill_manual(
            values = c("Interacting" = "darkgreen", "Avoiding" = "darkred"),
            name = "", na.translate = FALSE
        ) +
        scale_size_continuous(name = "% ROI\nsignificant\ninteractions\n", range = c(5, 15)) +
        ggplot2::labs(
            title = "Cell-Cell Interactions",
            subtitle = plot_subtitle,
            caption = "*empty tiles denote non-significant interactions"
        ) +
        ggplot2::xlab("from cell type ...") +
        ggplot2::ylab("to cell type ...") +
        theme_minimal(base_size = 16) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
            axis.text.y = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 20, face = "bold", colour = "grey50"),
            plot.title = element_text(size = 20, face = "bold"),
            plot.subtitle = element_text(size = 16, face = "italic"),
            plot.caption = element_text(size = 12, face = "italic")
        ) +
        guides(
            fill = guide_legend(
                order = 1,
                override.aes = list(size = 15)
            ),
            size = guide_legend(
                order = 2,
                override.aes = list(fill = "black")
            )
        )
    
    p <- add_highlight_regions(
        baseplot = p,
        highlight_colors = highlight_colors
        )
    
    return(p)
    
}


plot_ia_tile_fill <- function(ia_clean_df,
                              tile_fill = c("#FFF2CC", "#FFE699", "#FFD966", "#BF9000", "#7F6000")) {
    plot_subtitle <- glue::glue(
        "filters: {ia_clean_df$filter_on}",
        "filter values: {ia_clean_df$filter_vals}",
        "graph: delaunay triangulation",
        "count method: {unique(ia_clean_df$count_method)}",
        "p threshold: < 0.01",
        .sep = "\n"
    )
    
    p <- ia_clean_df %>%
        ggplot(
            aes(x = from_label, y = to_label)
        ) +
        geom_tile(
            aes(fill = pct),
            color = "grey50",
            linewidth = 0.1, linetype = 2
        ) +
        ggplot2::scale_fill_gradientn(
            colours = tile_fill,
            na.value = "white",
            name = " ",
            guide = guide_colourbar(
                theme = theme(
                    legend.key.height = unit(5, "cm"),
                    legend.key.width = unit(1, "cm"),
                ),
                frame.colour = "black",
                frame.linewidth = 0.35,
                ticks.colour = "black",
                ticks.linewidth = 0.35
            )
        ) +
        geom_point(
            aes(size = pct, color = type),
            fill = "white", shape = 21, stroke = 0.5, na.rm = TRUE
        ) +
        scale_color_manual(
            values = c("Interacting" = "darkgreen", "Avoiding" = "darkred"),
            name = "",
            na.translate = FALSE
        ) +
        scale_size_continuous(name = "% ROI\nsignificant\ninteractions\n", range = c(5, 15)) +
        ggplot2::labs(
            title = "Cell-Cell Interactions",
            subtitle = plot_subtitle,
            caption = "*empty tiles denote non-significant interactions"
        ) +
        ggplot2::xlab("from cell type ...") +
        ggplot2::ylab("to cell type ...") +
        theme_minimal(base_size = 16) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
            axis.text.y = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 20, face = "bold", colour = "grey50"),
            plot.title = element_text(size = 20, face = "bold"),
            plot.subtitle = element_text(size = 16, face = "italic"),
            plot.caption = element_text(size = 12, face = "italic")
        )
    
    p$layers[[2]]$aes_params$fill <- ifelse(ia_clean_df$type == "Interacting", "darkgreen", "darkred")
    
    p <- p +
        guides(
            color = guide_legend(
                order = 1,
                override.aes = list(
                    fill = c("darkgreen", "darkred"),
                    size = 12
                )
            ),
            size = guide_legend(
                order = 2,
                override.aes = list(fill = "black")
            )
        )
    
    p <- add_highlight_regions(baseplot = p)
    
    return(p)
}

add_highlight_regions <- function(baseplot,
                                  groups = rep(
                                      names(lab_spe@metadata$labels$cell_types),
                                      lengths(lab_spe@metadata$labels$cell_types)
                                  ),
                                  highlight_colors = "black",
                                  highlight_width = 3,
                                  legend_key_size = 6) {
    unique_groups <- unique(groups)
    
    highlight_regions <- data.frame(
        group = factor(unique_groups, levels = names(highlight_colors)),
        xmin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        xmax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        ymin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        ymax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        color = highlight_colors[unique(groups)]
    )
    
    if (length(highlight_colors) == 1) {
        baseplot +
            geom_rect(
                data = highlight_regions,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                color = highlight_colors, fill = NA,
                linewidth = highlight_width,
                inherit.aes = FALSE
            )
    } else {
        baseplot +
            geom_rect(
                data = highlight_regions,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
                fill = NA, linewidth = highlight_width, inherit.aes = FALSE
            ) +
            scale_color_manual(
                name = "", 
                values = highlight_colors
            ) +
            guides(
                color = guide_legend(override.aes = list(size =  legend_key_size)),
            )
    }
}