# Create the Heatmap Data
heatmap_data = anno_expression(object = initial_undefined, 
                               assay = "exprs", 
                               anno_label = "undefined_anno_1",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "minmax",
                               heatmap_scale_across = "features")


set_heatmap_global_opts(
    row_title_fill = heatmap_data$plot_opts$row_title_bkg_colours,
    row_title_text = heatmap_data$plot_opts$row_title_text_colours)

heatmap_data$ht = Heatmap(
    matrix = heatmap_data$heatmap_data, 
    col = heatmap_data$plot_opts$cols,
    rect_gp = gpar(col = "grey75", lwd = 0.70),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    heatmap_legend_param = list(
        title = heatmap_data$plot_opts$scale_title,
        at = heatmap_data$plot_opts$scale_breaks,
        legend_height = unit(25, "cm"),
        grid_width = unit(12.5, "mm")
    ),
    
    top_annotation = heatmap_data$column_anno,
    
    row_split = heatmap_data$plot_opts$rowsplits,
    
    row_names_side = "left",
    row_title_side = "left",
    row_names_max_width = unit(500, "mm"),
    row_gap = unit(5, "mm"),
    column_names_centered = TRUE
)


pdf(.new_file("Initial_undefined_anno", io$output$som_pheno),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "right")

dev.off()



