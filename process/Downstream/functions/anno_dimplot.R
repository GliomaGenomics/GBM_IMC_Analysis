.anno_dimplot = function(object, anno_label, reduction,
                         show_labels = FALSE, dot_size = 0.75, opacity = 0.75,
                         base_font_size = 25, anno_colors = NULL){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(missing(anno_label)) cli::cli_abort("missing anno_label argument")
    if(!reduction %in% SingleCellExperiment::reducedDimNames(object)) cli::cli_abort("{reduction} not found in {object}")
    
    number_of_colors = length(unique(object[[anno_label]]))
    
    if(is.null(anno_colors)){
        
        anno_colors = viridis::turbo(n = number_of_colors)
        
    }else{
        
        if(length(anno_colors) != number_of_colors) cli::cli_abort("incorrect number of anno colours supplied")
    }
    
    if(!anno_label %in% names(object@colData)){
        
        cli::cli_abort("{anno_label} not found in colData(object)")
    }
    
    out_plot_theme = ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    
    dittoSeq::dittoDimPlot(object, 
                           var = anno_label, 
                           reduction.use = reduction, 
                           size = dot_size,
                           opacity = opacity,
                           do.label = show_labels, 
                           color.panel = anno_colors) +
        ggplot2::ggtitle(glue::glue("{reduction} cells ({anno_label})")) +
        out_plot_theme
}


.plot_batches = function(object, batches, reduction,
                         show_labels = FALSE, dot_size = 0.75, opacity = 0.75,
                         base_font_size = 25, anno_colors = NULL){
    
    out_plots = vector("list", length = length(batches))
    
    for (batch in seq_along(batches)) {
        
        out_plots[[batch]]  = .anno_dimplot(object = object, 
                                            anno_label = batches[batch],
                                            reduction = reduction,
                                            show_labels = show_labels, 
                                            dot_size = dot_size,
                                            opacity = opacity,
                                            base_font_size = base_font_size
        )
        
        names(out_plots)[batch] = batches[batch]
    }
    
    out_plots = list(out_plots)
    names(out_plots) = reduction
    
    return(out_plots)
    
}

.anno_counts = function(object, anno_label, anno_colors = NULL, 
                        base_font_size = 25){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(missing(anno_label)) cli::cli_abort("missing anno_label argument")
    
    if(!anno_label %in% names(object@colData)){
        cli::cli_abort("anno_label not found in colData({quote(object)})")
    }
    
    if(!is.null(anno_colors)){
        
        if(length(anno_colors) != length(unique(object[[anno_label]]))) cli::cli_abort("length(anno_colors) != length(unique anno_label)")
        
    }else{
        
        anno_colors = viridis::turbo(length(unique(object[[anno_label]])))

    }
    
    as.data.frame(SummarizedExperiment::colData(object)) %>%
        dplyr::count(!!dplyr::sym(anno_label), name = "cells") %>%
        ggplot2::ggplot(aes(x = cells, 
                   y = !!dplyr::sym(anno_label), 
                   fill = !!dplyr::sym(anno_label))) +
        ggplot2::geom_bar(stat = "identity") +
        scale_fill_manual(values = anno_colors) +
        theme_classic(base_size = base_font_size) +
        theme(legend.position = "none")
    
}



visualise_clusters = function(object, anno_label, reduction,
                              dot_size = 0.5, opacity = 0.7,
                              base_font_size = 25,
                              show_labels = FALSE, 
                              anno_colors = NULL,
                              return_plots = TRUE,
                              save_plots = FALSE, 
                              outfile = NULL,
                              width = 20, height = 15){
    
    plots = list()
    
    plots$dimplot = .anno_dimplot(object = object, 
                                 anno_label = anno_label,
                                 reduction = reduction,
                                 dot_size = dot_size,
                                 opacity =  opacity,
                                 base_font_size = base_font_size,
                                 show_labels =show_labels, 
                                 anno_colors = anno_colors)
    
    plots$anno_counts = .anno_counts(object = object, 
                                    anno_label = anno_label,
                                    base_font_size = base_font_size,
                                    anno_colors = anno_colors)
    
    if(save_plots){
        
        if(is.null(outfile)) cli::cli_abort("please supply an outfile path/name")
        
        pdf(outfile, width = width, height = height, onefile = TRUE)
        
        print(plots)
        
        dev.off()
        
    }
    
    if(return_plots) return(plots) else return(invisible(NULL))
    
} 

