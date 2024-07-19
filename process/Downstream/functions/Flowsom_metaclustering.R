som_clusters = function(matrix, xdim, ydim){
    
    cli::cli_inform("SOM Clutering Using {xdim * ydim} Codes...")
    
    cli::cli_alert_info("Input matrix of {dim(matrix)[1]} rows and {dim(matrix)[2]} columns")
    
    cli::cli_progress_step("Creating a {xdim} * {ydim} grid")
    
    som_grid = kohonen::somgrid(xdim = xdim, ydim = ydim, topo = "hexagonal")
    
    set.seed(123)
    
    cli::cli_progress_step("Creating self-organising maps")
    
    som_model = kohonen::som(X = matrix, grid = som_grid)
    
    som_out = cbind(matrix, "SOM_cluster" = as.numeric(som_model$unit.classif))
    
    cli::cli_progress_step("Calculating median expression of SOM clusters")
    
    # Get the median expression for each som cluster
    cluster_center_expr <- data.frame(som_out) %>%
        dplyr::group_by_at("SOM_cluster") %>%
        dplyr::summarise_if(is.numeric, median)
    
    cli::cli_progress_step("Mapping SOM clusters back to cells")
    
    som_out = as.data.frame(som_out) %>%
        select(SOM_cluster) %>%
        tibble::rownames_to_column(var = "cell_id")

    cli::cli_progress_done()
    
    return(list(som_cluster_medians = cluster_center_expr,
                som_cluster_ids = som_out))
    
}

metaclustering = function(cluster_medians, 
                          cluster_col = "SOM_cluster", 
                          cell_cluster_ids, k = 30){
    
    cli::cli_progress_step("Metaclustering SOM codes using phenograph (k={k})")
    
    cluster_ids <- which(colnames(cluster_medians) == cluster_col)
    
    input_dims = dim(cluster_medians[,-cluster_ids])
    
    cli::cli_alert_info("Input data of {input_dims[1]} rows and {input_dims[2]} columns")
    
    # Open a sink to capture output
    sink(nullfile())
    
    phenograph_clusters = suppressMessages({
            Rphenograph::Rphenograph(data = cluster_medians[,-cluster_ids],  k = k)
            })
    
    # Close the sink to restore output
    sink()
    
    cli::cli_alert_info("Number of clusters: {length(unique(phenograph_clusters[[2]]$membership))}")
    
    cli::cli_progress_step("Mapping metaclusters back to cells")
    
    # Add meta clusters back to the input matrix
    cluster_medians$metacluster = phenograph_clusters[[2]]$membership
    
    mapped_clusters = cluster_medians %>%
        select(SOM_cluster, metacluster) %>%
        left_join(cell_cluster_ids, ., by = "SOM_cluster")
    
    cli::cli_progress_done()
    
    return(mapped_clusters)
    
}


som_metacluster = function(cluster_matrix, object, som_xdim = 10, som_ydim = 10,
                           metacluster_k = 30, outname = "SOM_metaclusters"){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(!class(object) %in% c("SpatialExperiment","SingleCellExperiment")) cli::cli_abort("object must be one of class: SingleCellExperiment, SpatialExperiment")
    
    if(missing(cluster_matrix)) cli::cli_abort("missing assay argument")
    
    som_codes = som_clusters(matrix = cluster_matrix,
                             xdim = som_xdim, ydim = som_ydim)
    
    som_metacodes = metaclustering(cluster_medians = som_codes$som_cluster_medians, 
                                   cell_cluster_ids = som_codes$som_cluster_ids, 
                                   k = metacluster_k)
    
    if(all(rownames(SummarizedExperiment::colData(object)) == som_metacodes$cell_id)){
        
        if(outname %in% names(object@metadata)){
            
            cli::cli_alert_warning("Overwritten {outname} within object metadata")
            
        }else cli::cli_alert_info("Added {outname} to object metadata")
        
        
        object@metadata[[outname]] = list(som_codes = som_codes,
                                          som_metacodes = som_metacodes)
        
    }
    
    
    beepr::beep(2)
    
    return(object)
    
}
