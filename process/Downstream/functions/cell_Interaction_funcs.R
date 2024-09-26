.clean_ia_data <- function(ia_df,
                           filter_by = NA,
                           filter_val = NA,
                           filter_min_patients = TRUE,
                           min_patients = 3,
                           count_method = c("histocat", "classic", "patch"),
                           label_order = levels(lab_spe$manual_gating)) {
  if (length(filter_val) > 1) stop("Only one filter value can be selected")

  filter_on <- if (is.na(filter_by)) "none" else filter_by
  filter_on_vals <- if (is.na(filter_val)) "none" else filter_val
  count_filter <- match.arg(count_method, several.ok = FALSE)

  if (filter_on != "none") {
    plot_data <- ia_df %>%
      dplyr::filter(!!sym(filter_on) %in% filter_on_vals) %>%
      dplyr::filter(count_method == count_filter)

    if (nrow(plot_data) == 0) stop("No data found for the given filter value")
  } else {
    plot_data <- ia_df %>%
      dplyr::filter(count_method == count_filter)
  }

  plot_data <- plot_data %>%
    dplyr::group_by(from_label, to_label) %>%
    dplyr::summarize(
      sum_sigval = sum(sigval, na.rm = TRUE),
      pct_sigval = abs(sum_sigval) / n() * 100,
      type = ifelse(sum_sigval == 0, NA, ifelse(sum_sigval > 0, "Interacting", "Avoiding")),
      .groups = "keep"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      filter_on = filter_on,
      filter_vals = paste0(filter_on_vals, collapse = ","),
      count_method = count_filter,
      min_patient_filtered = filter_min_patients,
      min_patients = ifelse(filter_min_patients, min_patients, NA),
      cell_state = ifelse(filter_on %in% "state_level", filter_on_vals, NA)
    ) %>%
    dplyr::mutate(
      across(type, ~ factor(., levels = c("Interacting", "Avoiding"))),
      across(c("from_label", "to_label"), ~ factor(., levels = label_order))
    )

  if (filter_min_patients) plot_data$type[abs(plot_data$sum_sigval) < min_patients] <- NA

  # Create a list of labels for the plot:
  title_suffix <- ifelse(
    filter_on_vals == "Prim", "Primary Samples",
    ifelse(filter_on_vals == "Rec", "Recurrent Samples",
      ifelse(filter_on_vals %in% c("low", "high"),
        paste(tools::toTitleCase(filter_on_vals), unique(ia_df$cell_state), sep = " "),
        "All Samples"
      )
    )
  )

  label_list <- list(
    title = glue::glue("Cell-Cell Interactions - {title_suffix}"),
    subtitle = glue::glue(
      "graph: delaunay triangulation",
      "count method: {unique(ia_df$count_method)}",
      "minimum patient filter: {min_patients}",
      "p threshold: < 0.01",
      .sep = "\n"
    ),
    caption = "\n*only significant interactions are displayed"
  )

  return(
    list(
      plot_df = plot_data,
      plot_labs = label_list
    )
  )
}


delta_surgery_interactions <- function(ia_df,
                                       patient_min = 3,
                                       title_suffix = NA) {
  surgery_split <- split(ia_df, ia_df$surgery)

  surgery_split <- lapply(surgery_split, \(x) .clean_ia_data(ia_df = x, min_patients = patient_min))

  delta <- surgery_split$Rec$plot_df %>%
    select(from_label, to_label, rec_type = type)

  delta <- surgery_split$Prim$plot_df %>%
    select(from_label, to_label,
      prim_sum_sigval = sum_sigval,
      prim_pct_sigval = pct_sigval,
      min_patient_filtered,
      min_patients,
      prim_type = type
    ) %>%
    left_join(delta, by = c("from_label", "to_label")) %>%
    mutate(
      interact_in = NA,
      interact_type = NA
    )

  delta$interact_in[!is.na(delta$prim_type) & is.na(delta$rec_type)] <- "Primary"
  delta$interact_in[is.na(delta$prim_type) & !is.na(delta$rec_type)] <- "Recurrent"
  delta$interact_in[!is.na(delta$prim_type) & !is.na(delta$rec_type)] <- "Both"

  delta$interact_type[which(delta$interact_in == "Primary")] <- as.character(delta$prim_type[which(delta$interact_in == "Primary")])
  delta$interact_type[which(delta$interact_in == "Recurrent")] <- as.character(delta$rec_type[which(delta$interact_in == "Recurrent")])


  # index of both interacting
  both_interacting <- which(delta$interact_in == "Both")

  prim_both <- as.character(delta$prim_type[both_interacting])
  rec_both <- as.character(delta$rec_type[both_interacting])
  both_same <- both_interacting[prim_both == rec_both]

  delta$interact_type[both_same] <- as.character(delta$prim_type[both_same])

  delta$interact_type[both_interacting[prim_both == "Interacting" & rec_both == "Avoiding"]] <- "Prim(I) -> Rec(A)"
  delta$interact_type[both_interacting[prim_both == "Avoiding" & rec_both == "Interacting"]] <- "Prim(A) -> Rec(I)"

  delta$interact_in <- factor(delta$interact_in, levels = c("Primary", "Recurrent", "Both"))
  delta$interact_type <- factor(delta$interact_type, levels = c("Interacting", "Avoiding", "Prim(I) -> Rec(A)", "Prim(A) -> Rec(I)"))

  if (!is.na(title_suffix)) {
    surgery_split$Prim$plot_labs$title <-  glue::glue("Cell-Cell Interactions - {title_suffix}")
  }

  return(
    list(
      plot_df = delta,
      plot_labs = surgery_split$Prim$plot_labs
    )
  )
}



plot_ia <- function(ia_clean_list,
                    point_size = 15,
                    highlight_colors = lab_spe@metadata$v2_colours$cell_groups) {
  p <- ia_clean_list$plot_df %>%
    ggplot(
      aes(x = from_label, y = to_label)
    ) +
    geom_tile(fill = "#fffdfa", color = "grey50", linewidth = 0.1, linetype = 2) +
    geom_point(
      data = subset(ia_clean_list$plot_df, !is.na(interact_type)),
      aes(fill = interact_type, shape = interact_in),
      color = "black", size = point_size, stroke = 0.5, na.rm = TRUE
    ) +
    scale_shape_manual(
      name = "",
      values = c(
        "Primary" = 21,
        "Recurrent" = 22,
        "Both" = 23
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c(
        "Interacting" = "darkgreen",
        "Avoiding" = "darkred",
        "Prim(I) -> Rec(A)" = "blue",
        "Prim(A) -> Rec(I)" = "purple"
      ),
      na.value = "white"
    ) +
    ggplot2::labs(
      title = ia_clean_list$plot_labs$title,
      subtitle = ia_clean_list$plot_labs$subtitle,
      caption = ia_clean_list$plot_labs$caption
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
      fill = guide_legend(order = 1, override.aes = list(shape = 21, size = point_size, color = "black")),
      shape = guide_legend(order = 2, override.aes = list(size = point_size))
    ) 
  
  p <- .add_highlight_regions(baseplot = p, highlight_colors = highlight_colors)

  return(p)
}


.add_highlight_regions <- function(baseplot,
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
        color = guide_legend(override.aes = list(size = legend_key_size)),
      )
  }
}



# code to plot points with shapes
# p <- ia_df_clean$prim %>%
#   ggplot(
#     aes(x = from_label, y = to_label)
#   ) +
#   geom_tile(fill = "#fffdfa", color = "grey50", linewidth = 0.1, linetype = 2) +
#   geom_point(
#     data = subset(ia_df_clean$prim, !is.na(type)),
#     aes(fill = type, shape = filter_vals),
#     color = "black", size = 15, stroke = 0.5, na.rm = TRUE
#   ) +
#   scale_shape_manual(
#     name = "",
#     values = c("Primary" = 21, "Recurrent" = 22, "Both" = 23)
#   ) +
#   scale_fill_manual(
#     name = "",
#     values = c("Interacting" = "darkgreen", "Avoiding" = "darkred"),
#     na.value = "white"
#   ) +
#   ggplot2::xlab("from cell type ...") +
#   ggplot2::ylab("to cell type ...") +
#   theme_minimal(base_size = 16) +
#   theme(
#     panel.grid.major = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
#     axis.text.y = element_text(size = 16, face = "bold"),
#     axis.title = element_text(size = 20, face = "bold", colour = "grey50"),
#     plot.title = element_text(size = 20, face = "bold"),
#     plot.subtitle = element_text(size = 16, face = "italic"),
#     plot.caption = element_text(size = 12, face = "italic")
#   ) +
#   guides(
#     fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 15, color = "black")),
#     shape = guide_legend(order = 2, override.aes = list(size = 15))
#   )
