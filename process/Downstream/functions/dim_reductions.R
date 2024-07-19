runUMAP <- function(object, reduct_matrix, reducedDim_suffix, n_neighbors = 35){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(!class(object) %in% c("SpatialExperiment","SingleCellExperiment")) cli::cli_abort("object must be one of class: SingleCellExperiment, SpatialExperiment")
    
    if(missing(reduct_matrix)) cli::cli_abort("missing reduct_matrix argument")

    set.seed(123)
    
    reduced_dimensions <- uwot::umap(reduct_matrix, 
                                     n_neighbors = n_neighbors,
                                     min_dist = 0.01,
                                     metric = "cosine",
                                     n_epochs = 200,
                                     learning_rate = 1,
                                     verbose = TRUE)
    
    reduction_label = glue::glue("UMAP_{reducedDim_suffix}")
    
    if(reduction_label %in% reducedDimNames(object)){
        
        cli::cli_progress_step("Overwriting {reduction_label} reducedDim")
    
    }else{
        
        cli::cli_progress_step("Adding {reduction_label} to reducedDims")
    }
    
    reducedDim(object, reduction_label) = reduced_dimensions

    return(object)
    
}
