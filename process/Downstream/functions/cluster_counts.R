# For a given cluster label return the counts of members in each cluster 
cluster_cell_counts <- function(spe, cluster_label, 
                                cluster_colors = NULL,
                                plot_title ="",
                                scale_cols = T){
    
    if(!class(spe) %in% c("SpatialExperiment", "SingleCellExperiment")) stop("spe is not the correct class")
    
    if(!cluster_label %in% names(colData(spe))) stop("cluster_label not found in spe")
    
    plot_table = table(colData(spe)[cluster_label]) %>%
        as.data.frame() %>%
        dplyr::rename(cluster = 1,
                      cells = Freq) %>%
        dplyr::mutate(dplyr::across(cluster, factor))
        
    if(scale_cols){
        scale_breaks <- c(500,
                          seq(2500, 5000, by = 2500),
                          seq(10000, max(plot_table$cells), by = 5000))
    }
    
 
    
    outplot <- ggplot(plot_table, aes(cells, cluster, fill = cluster)) +
        
        geom_col() +
        
        ggtitle(plot_title) +
        
        {if(is.null(cluster_colors)) scale_fill_viridis(discrete = T, option = "H")} +
        
        {if(!is.null(cluster_colors)) scale_fill_manual(values = cluster_colors)} +
        
        {if(scale_cols) scale_x_continuous(breaks = scale_breaks)} +
        
        theme_classic(base_size = 12) +
        
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", 
                                        colour = "black", 
                                        size = 16, 
                                        hjust = 0.5),
              panel.grid.major.x = element_line(colour = "lightgrey",
                                                size = 0.5,
                                                linetype = "dotted")
        )
    
    return(outplot)
    
}
