.custom_cols <- function(values, n_breaks = 10, 
                         low_color = "grey97", high_color = "darkred", 
                         low_color_pct = 0.5) {
    
    breaks = pretty(values, n = n_breaks)
    
    color_palette = colorRampPalette(c(low_color, high_color))
    
    colors = color_palette(length(breaks))
    
    n_low_colors <- round(n_breaks * low_color_pct)
    colors[1:n_low_colors] <- low_color
    
    return(list(colors = colors, breaks = breaks))
    
}


cluster_dotplot = function(spe, assay, markers, cluster_col, 
                           summary_fun = c("mean", "median"),
                           dot_size = 50,
                           base_font_size = 20,
                           high_color = "darkred", low_color = "grey97",
                           low_color_pct = 0.5){
    
    summary_function = switch(
        EXPR = match.arg(summary_fun),
        mean = function(x) {mean(x[x != 0])},
        median = function(x) {median(x[x != 0])}
    )
    
   summary_expression = dittoSeq::dittoDotPlot(
       object = spe,
       assay = "minmax",
       vars = markers,
       group.by = cluster_col,
       summary.fxn.color = summary_function,
       data.out = TRUE)
   
    plot_cols = .custom_cols(summary_expression$data$color,
                             low_color = low_color, high_color = high_color, 
                             low_color_pct = low_color_pct)
   
    ggplot(summary_expression$data, 
           aes(x = grouping, y = var, color = color, group = grouping)) +
        geom_point(stat = "identity", size = dot_size) +
        theme_classic(base_size = base_font_size) +
        scale_color_gradientn(colors = plot_cols$colors,
                              breaks = plot_cols$breaks,
                              guide = guide_colorbar(
                                  title = "Z-Score",
                                  title.position = "top",
                                  direction = "vertical",
                                  barwidth = 3.5,
                                  barheight = 25,
                                  label.hjust = 1,
                                  ticks.colour = "grey15",
                                  frame.colour = "black",
                                  ticks.linewidth = 1.5,
                                  ticks = TRUE)
        ) +
        theme(axis.title = element_blank()) 
    
    
}




