anno_expression = function(object, assay, anno_label,
                           subset_rows = NULL, subset_columns = NULL,
                           anno_summary_stat = c("mean", "median"),
                           heatmap_scale_fun = c("minmax","zscore","mean", "none", 
                                                 "robust_scalar", "unit_lengths"),
                           heatmap_scale_across = c("features", "samples"),
                           heatmap_scale_breaks = 10,
                           scale_colors = c("#A50026","#FFFFBF","#313695"),
                           colanno_name = "Number of cells",
                           col_anno_name_size = 35,
                           col_anno_name_face = "italic",
                           col_anno_side = "left"){
    
    stopifnot(exprs = {
        
        !missing(object)
        !missing(assay)
        !missing(anno_label)
        
        class(object) %in% c("SpatialExperiment","SingleCellExperiment")
        assay %in% SummarizedExperiment::assayNames(object)
        anno_label %in% names(object@colData)
    })
    
    anno_summary_stat = match.arg(anno_summary_stat, several.ok = FALSE)
    
    # Aggregate(using summary statistic) across the annotation labels
    summarised_exprs = scuttle::summarizeAssayByGroup(
        x = object, 
        ids = object@colData[[anno_label]],
        assay.type = assay,
        subset.row = subset_rows,
        subset.col = subset_columns,
        statistics = anno_summary_stat,
        store.number = "ncells"
        )
    
    # Generate the heatmap data
    heatmap_data = SummarizedExperiment::assay(summarised_exprs, anno_summary_stat)

    # Scale the heatmap data using a defined function
    heatmap_scale_fun = match.arg(heatmap_scale_fun, several.ok = FALSE)
    
    heatmap_data = switch( 
        EXPR = match.arg(heatmap_scale_across, several.ok = FALSE),
        
        "features" = t(apply(t(heatmap_data), 2, .norm_fun(heatmap_scale_fun))),
        
        "samples" =  apply(heatmap_data, 2, .norm_fun(heatmap_scale_fun))
        
        )
    
    # Add column annotation for ncells
    column_anno = ComplexHeatmap::HeatmapAnnotation(
        colanno_name = ComplexHeatmap::anno_barplot(colData(summarised_exprs)$ncells,
                             which = "column",
                             bar_width = 0.8,
                             axis = FALSE,
                             add_numbers = TRUE,
                             numbers_gp = grid::gpar(fontsize = 20,
                                                     fontface = "bold"),
                             numbers_rot = 0,
                             numbers_offset = grid::unit(5, "mm"),
                             height = grid::unit(65, "mm"),
                             gp = grid::gpar(fill = "slateblue"),
                             border = TRUE),

        annotation_label = colanno_name,
        annotation_name_gp = grid::gpar(fontsize = col_anno_name_size,
                                        fontface = col_anno_name_face),
        annotation_name_offset = grid::unit(3, "mm"),
        annotation_name_side = col_anno_side)
    
    
    # Heatmap option ----
    plot_opts = list()
    
    alt_colors = function(x){
        
        num_colors = length(x)
        
        alt_index = c(matrix(c(1:num_colors, num_colors:1), nrow = 2, byrow = TRUE))
        x[alt_index[1:num_colors]]
        
    }
    
    plot_opts$row_title_bkg_colours = alt_colors(viridis::turbo(length(unique(names(subset_rows)))))
    
    get_text_color = function(bg_color, luminance_threshold = 0.5) {
        # Convert background color to RGB format
        bg_rgb <- grDevices::col2rgb(bg_color)
        
        # Calculate the luminance of the background color
        luminance <- (bg_rgb[1] * 0.299 + bg_rgb[2] * 0.587 + bg_rgb[3] * 0.114) / 255
        
        # Return white text color for dark background, black text color for light background
        if (luminance > luminance_threshold) {
            return("black")
        } else {
            return("white")
        }
    }
    get_text_color = Vectorize(get_text_color)
    
    if(!is.null(subset_rows)){
        
        plot_opts$rowsplits = factor(names(subset_rows), 
                                     levels = unique(names(subset_rows)))
        
    }else plot_opts$rowsplits = subset_rows
    
    plot_opts$row_title_text_colours = get_text_color(plot_opts$row_title_bkg_colours)
    
    plot_opts$scale_title = glue::glue("scaled ({heatmap_scale_fun}) expression")
    
    plot_opts$scale_breaks = pretty(heatmap_data, n = heatmap_scale_breaks)

    if(heatmap_scale_fun == "minmax"){  
        
        plot_opts$cols =  circlize::colorRamp2(
            
            breaks = plot_opts$scale_breaks, 
            
            colors = c(rep("#000004FF", 6),
                       "#420A68FF","#781C6DFF","#ED6925FF","#FCA50AFF", "#FCFFA4FF")
            )
        
    }else{
        
        generate_color_palette = function(breaks, colors) {
            n = length(breaks) - 1
            color_func = colorRampPalette(colors)
            color_palette = color_func(n)
            return(color_palette)
        }
        
        plot_opts$cols = generate_color_palette(
            breaks = plot_opts$scale_breaks,
            colors = scale_colors
            ) 
        
        
    }
    
    
    
    

    # ----
    
    
    return(
        list(heatmap_data = heatmap_data, 
             column_anno = column_anno,
             plot_opts = plot_opts)
        )
    
}




















set_heatmap_global_opts = function(dim_fontsize=30, 
                                   row_title_size = 35,
                                   row_title_face = "bold",
                                   row_title_fill = NULL,
                                   row_title_text = NULL,
                                   row_title_border = "dashed",
                                   row_title_border_col = "black",
                                   title_padding_mm = 3.5,
                                   row_col_name_padding_mm = 3,
                                   column_anno_padding = 3,
                                   legend_title_size = 30,
                                   legend_title_face = "bold",
                                   legend_text_size = 25,
                                   legend_text_face = "bold",
                                   legend_title_position = "leftcenter-rot",
                                   legend_border_col = "grey85",
                                   heatmap_border = TRUE
                                   ){
    
    # Rest all previous options
    ComplexHeatmap::ht_opt(RESET = TRUE)
    
    # Row and column names
    ht_opt$heatmap_row_names_gp = gpar(fontsize = dim_fontsize)
    ht_opt$heatmap_column_names_gp = gpar(fontsize = dim_fontsize)
    
    # Row title
    ht_opt$heatmap_row_title_gp = gpar(
        fontsize = row_title_size, 
        fontface = row_title_face,
        fill = row_title_fill,
        col = row_title_text,
        lty = row_title_border,
        border = row_title_border_col)
    
    # space between row/column titles and heatmap body.
    ht_opt$TITLE_PADDING = unit(title_padding_mm, "mm")
    
    # space between row/column names and heatmap body
    ht_opt$DIMNAME_PADDING = unit(row_col_name_padding_mm, "mm")
    
    # space between column annotation and heatmap body
    ht_opt$COLUMN_ANNO_PADDING = unit(column_anno_padding, "mm")
    
    # legend title 
    ht_opt$legend_title_gp =  gpar(fontsize = legend_title_size, 
                                   fontface = legend_title_face)
    # legend text 
    ht_opt$legend_labels_gp = gpar(fontsize = legend_text_size,
                                   fontface = legend_text_face)
    # legend title position
    ht_opt$legend_title_position = legend_title_position
    
    # legend border
    ht_opt$legend_border = legend_border_col
    
    # heatmap border
    ht_opt$heatmap_border = heatmap_border
}
