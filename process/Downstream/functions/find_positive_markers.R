find_positive_cells <- function(marker_exprs, cell_markers, 
                                high_thresholds, low_thresholds, 
                                percentage_threshold = 0.4) {
    
    num_low_genes <- sum(marker_exprs[, !(colnames(marker_exprs) %in% cell_markers)] < low_thresholds[!(colnames(marker_exprs) %in% cell_markers)])
    required_low_genes <- ceiling((ncol(marker_exprs) - length(cell_markers)) * percentage_threshold)
    
    positive_cells <- rowSums(marker_exprs[, cell_markers, drop = F] > high_thresholds[cell_markers]) == length(cell_markers) &
        num_low_genes >= required_low_genes
    
    return(which(positive_cells))
}



find_high_expressing_cells <- function(matrix_data, 
                                       quantile_threshold = rep(0.5, ncol(matrix_data))) {
    # Check if quantile threshold length matches the number of columns in the data matrix
    if (length(quantile_threshold) != ncol(matrix_data)) {
        stop("Quantile threshold length must match the number of columns in the data matrix.")
    }
    
    # Calculate high expression thresholds based on quantile
    high_thresholds <- sapply(seq_along(quantile_threshold), function(i) {
        gene_exprs <- matrix_data[, i]
        quantile_value <- quantile(gene_exprs, quantile_threshold[i])
        quantile_value
    })
    
    # Find cells with high expression for each marker
    high_cells <- lapply(seq_along(quantile_threshold), function(i) {
        gene_exprs <- matrix_data[, i]
        threshold <- high_thresholds[i]
        gene_exprs > threshold
    })
    
    # Find cells with high expression for all markers
    positive_cells <- Reduce(`&`, high_cells)
    
    positive_cells = which(positive_cells)
    
    cli::cli_alert_info("{length(positive_cells)} cells found")
    
    return(positive_cells)
}
