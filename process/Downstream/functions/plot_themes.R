.horizonal_facet_prop_theme <- function(text_size = 20,
                                       show_facet_strip_text = FALSE,
                                       show_x_gridlines = TRUE,
                                       show_y_axis_text = FALSE,
                                       show_x_axis_text = TRUE,
                                       show_panel_border = TRUE,
                                       panel_background = "#FFFFFA",
                                       panel_border_color = "black",
                                       panel_border_type = 1,
                                       panel_border_size = 1,
                                       panel_gridline_x_color = "grey",
                                       panel_gridline_major_x_line_width = 0.5,
                                       panel_gridline_minor_x_line_width = 0.5,
                                       facet_stip_text_size = 25,
                                       panel_spacing_mm = 0.5,
                                       legend_key_size = 15,
                                       x_axis_text_angle = 0,
                                       x_axis_tick_color = "grey",
                                       x_axis_tick_size = 0.5,
                                       ...
                                       ) {
  out <- ggplot2::theme_classic(base_size = text_size) +
    ggplot2::theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = "black", hjust = 0.5
      ),
      legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(legend_key_size, "mm"),
      legend.text = ggplot2::element_text(size = text_size),
      strip.background = ggplot2::element_rect(linetype = "blank"),
      panel.spacing = ggplot2::unit(panel_spacing_mm, "mm")
    )


  if (show_panel_border) {
    out <- out + ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = panel_border_color,
        linewidth = panel_border_size,
        linetype = panel_border_type
      ),
      panel.background = ggplot2::element_rect(fill = panel_background)
    )
  } else {
    out <- out + ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
  }

  if (show_facet_strip_text) {
    out <- out + ggplot2::theme(
      strip.text = ggplot2::element_text(
        size = facet_stip_text_size
      )
    )
  } else {
    out <- out + ggplot2::theme(
      strip.text = ggplot2::element_blank()
    )
  }

  if (show_x_gridlines) {
    out <- out + ggplot2::theme(
      panel.grid.major.x = element_line(
        colour = panel_gridline_x_color,
        linewidth = panel_gridline_major_x_line_width,
        linetype = 1
      ),
      panel.grid.minor.x = element_line(
        colour = panel_gridline_x_color,
        linewidth = panel_gridline_minor_x_line_width,
        linetype = 1
      )
    )
  } else {
    out <- out + ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  }


  if (show_x_axis_text) {
    out <- out + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, vjust = 0.5),
      axis.ticks.x = element_line(colour = x_axis_tick_color, size = x_axis_tick_size),
    )
  } else {
    out <- out + ggplot2::theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  if (show_y_axis_text) {
    out <- out + ggplot2::theme(
      axis.text.y = element_text(
        color = "black",
        face = "bold",
        size = 32,
        margin = margin(t = 0, r = 0, l = 5, b = 0, unit = "mm")
      )
    )
  } else {
    out <- out + ggplot2::theme(
      axis.text.y = element_blank()
    )
  }
  
  out <- out + theme(...)

  return(out)
}
