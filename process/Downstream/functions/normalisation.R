.norm_fun = function(fun = c("zscore","minmax",
                             "mean","none",
                             "robust_scalar", 
                             "unit_lengths")){
    
    norm_function = switch(match.arg(fun, several.ok = FALSE),
                           
                           zscore = function(x){ scale(x, center = T, scale = T) },
                           
                           robust_scalar = function(x){ (x- median(x)) /(quantile(x,probs = .75)-quantile(x,probs = .25)) },
                           
                           minmax = function(x){ (x- min(x)) /(max(x)-min(x)) },
                           
                           mean = function(x){ (x- mean(x)) /(max(x)-min(x)) },
                           
                           unit_lengths = function(x) { x / sqrt(sum(x^2)) },
                           
                           none = function(x) {x}
                           
    )
    
    return(norm_function)
    
}

.clip_expression_vals = function(exprs_matrix, quantile = 0.999){
    
    # Clip non-zero values at a desired quantile
    q = stats::quantile(exprs_matrix[exprs_matrix>0], quantile)
    
    cli::cli_progress_step("Clipping non-zero values at {quantile * 100}%, q={round(q, 3)}")
    
    exprs_matrix[exprs_matrix > q] = q
    
    return(exprs_matrix)
}


.scale_expression = function(exprs_matrix, scale_fun = "zscore"){
    
    cli::cli_progress_step("Scaling expression matrix using {scale_fun}")
    
    scaled_expr = apply(exprs_matrix, 2, .norm_fun(scale_fun))

    if(length(rownames(scaled_expr)) == 0) rownames(scaled_expr) = rownames(exprs_matrix)
    
    cli::cli_progress_done()
    
    t(scaled_expr)    
    
}


clean_expression = function(object, assay, scale_fun = "zscore", 
                            clip_exprs_vals = TRUE, clip_quantile = 0.999){
    
    if(missing(object)) cli::cli_abort("missing object argument")
    if(!class(object) %in% c("SpatialExperiment","SingleCellExperiment")) cli::cli_abort("object must be one of class: SingleCellExperiment, SpatialExperiment")
    
    if(missing(assay)) cli::cli_abort("missing assay argument")
    if(!assay %in% SummarizedExperiment::assayNames(object)) cli::cli_abort("assay not in assayNames(object)")
    
    expr = t(assay(object, assay))
    
    cli::cli_alert_info("Input matrix rows: {dim(expr)[1]}")
    cli::cli_alert_info("Input matrix columns: {dim(expr)[2]}")
    
    
    if(clip_exprs_vals){
        
        expr = .clip_expression_vals(exprs_matrix = expr, 
                                    quantile = clip_quantile) 
        
    }
    
    expr = .scale_expression(exprs_matrix = expr, 
                            scale_fun = scale_fun)
    
    
    if(scale_fun %in% SummarizedExperiment::assayNames(object)){
        
        cli::cli_progress_step("Overwriting {scale_fun} assay")
    
        }else cli::cli_progress_step("Adding {scale_fun} assay to {quote(object)}")
    
    
    assay(object, scale_fun) = expr
    
    cli::cli_progress_done()
    
    return(object)

}
