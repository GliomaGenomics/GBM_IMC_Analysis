plot_cluster_expression = function(spe, 
                                   assay = "exprs", 
                                   markers, 
                                   cluster_labels, 
                                   fun = c("mean", "median"),
                                   scale_fun = c("zscore", "minmax"),
                                   scale = c("none","row","column"),
                                   plot_title = "",
                                   outfile){
    
    if(!class(spe) %in% c("SpatialExperiment", "SingleCellExperiment")) stop("spe is not the correct class")
    if(!all(markers %in% rownames(spe))) stop("all marker not present in spe")
    if(!assay %in% assayNames(spe)) stop("assay not found in spe")
    if(!cluster_labels %in% names(colData(spe))) stop("cluster_labels not found in spe")
    
    expr = t(assay(spe, assay))[, markers, drop = F]
    rownames(expr) = spe@colData[[cluster_labels]]
    
    expr = switch (match.arg(fun, several.ok = FALSE),
        mean = t(sapply(by(expr, rownames(expr), colMeans), identity)),
        median = t(sapply(by(expr, rownames(expr), \(x) apply(x,2, median)), identity))
        )
    
    if(is.factor(spe@colData[[cluster_labels]])){
       
        expr = expr[match(levels(spe@colData[[cluster_labels]]), rownames(expr)),] 
        
    }
    
    should_scale_on = match.arg(scale, several.ok = F)
    
    if(should_scale_on == "none"){
    
        print("No scaling performed")
            
        pheatmap::pheatmap(mat = t(expr), 
                          cluster_rows = F,
                          cluster_cols = F,
                          scale = should_scale_on,
                          main = plot_title,
                          labels_row = colnames(expr),
                          color = viridis::inferno(1000),
                          silent = T,
                          filename = outfile)
        }
    
    what_scaling_fun = match.arg(scale_fun, several.ok = F)
    
    if(what_scaling_fun == "minmax"){
        
        print("Min/Max scaling performed")
        
       expr = switch(should_scale_on,
            row = t(apply(expr, 1, \(x) scales::rescale(x, to=c(0,1)))),  # min/max by row,
            column = apply(expr, 2, \(x) scales::rescale(x, to=c(0,1)))     # min/max by col
            )
       
      pheatmap::pheatmap(mat = t(expr), 
                          cluster_rows = F,
                          cluster_cols = F,
                          scale = "none",
                          main = plot_title,
                          labels_row = colnames(expr),
                          color = viridis::inferno(1000),
                          silent = T, 
                         filename = outfile)
    }

    if(what_scaling_fun == "zscore"){
        
        print("Zscore scaling performed")
        
        pheatmap::pheatmap(mat = t(expr), 
                                 cluster_rows = F,
                                 cluster_cols = F,
                                 scale = should_scale_on,
                                 main = plot_title,
                                 labels_row = colnames(expr),
                                 color = viridis::inferno(1000),
                                 silent = T,
                                 filename = outfile)
    }
    
    }


