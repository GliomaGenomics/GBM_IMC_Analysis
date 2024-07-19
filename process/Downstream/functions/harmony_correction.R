harmony_correction = function(object, assay, features = NULL, correct_vars,
                              pca_before_correction = TRUE, npcs = 20,
                              dim_reduction = TRUE, n_neighbors = 35){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(!class(object) %in% c("SpatialExperiment","SingleCellExperiment")) cli::cli_abort("object must be one of class: SingleCellExperiment, SpatialExperiment")
    
    if(missing(assay)) cli::cli_abort("missing assay argument")
    if(!assay %in% SummarizedExperiment::assayNames(object)) cli::cli_abort("assay not in assayNames(object)")
    
    if(missing(correct_vars)) cli::cli_abort("correct_vars are missing")
    
    if(!all(correct_vars %in% names(SummarizedExperiment::colData(object)))) cli::cli_abort("all correct_vars not in colData(object)")
    
    harmony_matrix = t(SummarizedExperiment::assay(object, assay))
    
    if(!is.null(features)){
        
        if(any(duplicated(features))) cli::cli_abort("features are not all unique")  
        
        if(!all(features %in% colnames(harmony_matrix))) cli::cli_abort("features missing in assay")
        
        harmony_matrix = harmony_matrix[,features]
        
    }
    
    harmony_emb <- harmony::HarmonyMatrix(data_mat = harmony_matrix, 
                                          meta_data = colData(object),
                                          vars_use =  correct_vars,
                                          do_pca = pca_before_correction,
                                          npcs = npcs)
    
    
    corrected_embeddings_label = glue::glue("harmony_{assay}")
    
    if(corrected_embeddings_label %in% reducedDimNames(object)){
        
        cli::cli_progress_step("Overwriting {corrected_embeddings_label} reducedDim")
        
    }else{
        
        cli::cli_progress_step("Adding {corrected_embeddings_label} to reducedDims")
    }
    
    SingleCellExperiment::reducedDim(object, corrected_embeddings_label) <- harmony_emb
    
    if(dim_reduction){
        
        set.seed(1234)
        
        object = runUMAP(
            object = object,
            reduct_matrix = SingleCellExperiment::reducedDim(object, corrected_embeddings_label),
            reducedDim_suffix = corrected_embeddings_label,
            n_neighbors = n_neighbors)
    }
    
    return(object)
    
}