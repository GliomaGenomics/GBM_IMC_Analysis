marker_exprs <- function(spe, assay = "exprs", reducedDim = "UMAP_harmony",
                         markers = rownames(spe), scale = FALSE,
                         col = c("darkblue", "yellow"), 
                         pointsize = 0.2, ncol = 6,
                         xlab = "UMAP1", ylab = "UMAP2"){
  
  th <- theme_bw(base_size = 24) + 
    theme(
      legend.background = element_rect(), 
      title = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", vjust = 1, hjust = 0.5), 
      plot.subtitle = element_text(size = 12, vjust = 1, hjust = 0.5),
      plot.caption = element_text(size = 12, vjust = 1), 
      axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5), 
      axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5), 
      axis.title = element_text(size = 16), 
      legend.title = element_blank(), 
      legend.position = "right", 
      legend.key = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "cm"), 
      legend.text = element_text(size = 12), 
      legend.key.height = unit(2.5, "cm"), 
      strip.text.x = element_text(size = 16, face = "bold", margin = margin(b = 5, t = 5))
      
    )
  
  plotobj <- reducedDim(spe, reducedDim)
  
  colnames(plotobj) <- paste0("Dim", 1:ncol(plotobj))
  
  exprs <- t(assay(spe, assay))[,markers]
  
  if(scale) exprs <- apply(exprs, 2, function(x)  x/max(x))
  
  plotobj <- as.data.frame(cbind(plotobj, exprs))
  
  plotobj <- pivot_longer(plotobj, -starts_with("Dim"), 
                          names_to = "marker", values_to = "exprs")
  
  plot <- ggplot(plotobj, aes(x = Dim1, y = Dim2, alpha = exprs)) + 
    th +
    guides(fill = guide_legend(), shape = guide_legend(), alpha = "none")
  
  plot <- plot + geom_point(aes(colour = exprs), size = pointsize)
  
  if (length(col) == 2) {
    
    plot <- plot + scale_colour_gradient(low = col[1], high = col[2], 
                                         name = "Expression")
  }
  
  if (length(col) == 3) {
    
    plot <- plot + scale_colour_gradient2(low = col[1], mid = col[2], high = col[3], 
                                          midpoint = colMidpoint, 
                                          limits = c(min(plotobj$Expression), 
                                                     max(plotobj$Expression)), 
                                          space = "Lab", 
                                          name = "Expression")
  }
  
  if (length(col) >3){
    
    print("Colour Length >3... Using first two cols")
    
    plot <- plot + scale_colour_gradient(low = col[1], high = col[2], 
                                         name = "Expression")
    
  }
  
  nrow <- ceiling(length(markers)/ncol)
  
  plot <- plot + facet_wrap(~marker, nrow = nrow, ncol = ncol)
  
  plot <- plot + xlab(xlab) + ylab(ylab)
  
  return(plot)
  
}



cluster_exprs <- function(spe, assay = "counts", 
                          cluster_var, 
                          markers = rownames(spe)[rowData(spe)$downstream == 1],
                          exprs_stat = c("median","mean"),
                          scale = TRUE,
                          plot_title = "Marker Expression") {
  
  if (!class(spe) %in% c("SingleCellExperiment","SpatialExperiment")){
    
    stop("spe not of class SingleCellExperiment or SpatialExperiment")
    
  }
  
  if (!assay %in% assayNames(spe)) stop("Assay not found in assayNames(spe)")
  
  if (!cluster_var %in% names(colData(spe))) stop("cluster_var not found in colData(spe)")
  
  
  stat_method <- match.arg(exprs_stat, 
                           choices = c("median", "mean"),
                           several.ok = FALSE)
  
  if(scale == TRUE) {
    
    expr <- t(assay(spe, assay))[,markers]
    
    # Take the upper-lower quantile
    col_quant <- colQuantiles(expr, probs = c(0.01, 0.99))
    
    # Scale to min max
    scaled_expression_mat <- t((t(expr) - col_quant[, 1])/(col_quant[,2] - col_quant[, 1]))
    
    # set upper lower bounds
    scaled_expression_mat[scaled_expression_mat < 0] <- 0
    scaled_expression_mat[scaled_expression_mat > 1] <- 1
    
    expr <- scaled_expression_mat
    
  }else expr <- t(assay(spe, assay))[,markers]
  
  # Add grouping variable to expression 
  expr <- as.data.frame(expr)
  
  expr <- cbind(expr, clusters = as.character(spe[[cluster_var]]))
  
  if(stat_method %in% "median"){
    
    cluster_stat <- expr %>% 
      
      group_by(clusters) %>% 
      
      summarise(across(.fns = median))
  }
  
  if(stat_method %in% "mean"){
    
    cluster_stat <- expr %>% 
      
      group_by(clusters) %>% 
      
      summarise(across(.fns = mean))
  }
  
  
  # Calculate cluster stat expression
  cluster_stat <- tidyr::pivot_longer(cluster_stat, 
                                      cols = -c(clusters), 
                                      names_to = "marker", 
                                      values_to = "exprs") 
  
  plot_data <- cluster_stat %>%
    
    mutate(across(clusters, ~factor(.x, levels = levels(spe[[cluster_var]])))) %>%
    
    mutate(across(marker, ~factor(.x, levels = markers))) %>%
    
    drop_na(clusters)
    
  
  legend_title <- ifelse(scale, 
                         glue::glue("{assay}\n(scaled)"), 
                         glue::glue("{assay}")
  ) %>%
    tools::toTitleCase()
  
  
  out <-  ggplot(plot_data, aes(clusters, marker, fill = exprs)) + 
    geom_tile(color = "grey",
              lwd = 0.1,
              linetype = 1) + 
    ggtitle(plot_title) +
    
    scale_fill_viridis(option = "B", limits = c(0.2, 1), oob = scales::squish) +
    
    # scale_fill_viridis(option = "B") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(face = "bold", colour = "black", size = 14),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(colour = "black",
                                      vjust = 3.5, hjust = 0.5),
          legend.key.height = unit(1.5, "cm"))
  
  out$labels$fill <- legend_title
  
  
  return(out)
  
}


cluster_cell_counts <- function(spe, clusters = c("som","pg")){
  
  cluster_method <- match.arg(clusters, 
                              choices = c("som", "pg"),
                              several.ok = FALSE)
  
  clusters <- grep(glue::glue("(?i)^{cluster_method}"), 
                   names(colData(spe)),value = T)
  
  cluster_cols <- spe@metadata[["color_vectors"]][[clusters]]
  
  if(length(clusters) != 1) stop("Clusters not found!", call. = FALSE)
  
  if(grepl("(?i)^pg", clusters)) plot_title <- "Phenograph Clustering"

  if(grepl("(?i)^som", clusters)) plot_title <- "SOM Clustering"
  
  
  plot_table <- table(colData(spe)[clusters]) %>% 
    
    as.data.frame() %>%
    
    rename(Cluster = Var1,
           Ncells = Freq) %>%
    
    mutate(across(Cluster, factor))
  
  scale_breaks <- c(500,
                    seq(2500, 5000, by = 2500,),
                    seq(10000, max(plot_table$Ncells), by = 5000)
                    )
  
  outplot <- ggplot(plot_table, aes(Ncells, Cluster, fill = Cluster)) +
    
    geom_col() +
    
    ggtitle(plot_title) +
    
    scale_fill_manual(values = cluster_cols) +
    
    scale_x_continuous(breaks = scale_breaks) +
    
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


create_cluster_labels <- function(clust_anno, labels, clusters){
  
  clust_split  <- str_split(clust_anno[[Clusters]], ",")
  
  clust_labs <- unlist(clust_split, use.names = FALSE)
  
  dup_clusters <- unique(clust_labs[duplicated(clust_labs)])
  
  clust_split <- lapply(clust_split,  function(x) x[!x %in% dup_clusters])
  
  clust_labs <- clust_labs[!clust_labs %in% dup_clusters]
  
  out <- clust_labs
  
  names(out) <- rep(clust_anno$Label, lengths(clust_split))
  
  # out <- list(old = clust_labs, 
  #             new = rep(clust_anno$Label, lengths(clust_split))
  #             )
  
  return(out)
}


combine_annotations <- function(x){
  
  if(!is.data.frame(neoplastic_clusters)) stop("x must be of class dataframe")
  
  unite(data = x, dplyr::everything(),
        col = "combined_anno",
        na.rm = TRUE) %>%
    
    mutate(across(combined_anno, ~ifelse(.x %in% "", NA, .x)))
  
}


scale01 <-function(spe, assay){
  
  if (!assay %in% assayNames(spe)) stop("Assay not found in assayNames(spe)", call. = FALSE)
  
  expr <- t(assay(spe, assay))
  
  # Take the upper-lower quantile
  col_quant <- colQuantiles(expr, probs = c(0.01, 0.99))
  
  # Scale to min max
  scaled_expression_mat <- t((t(expr) - col_quant[, 1])/(col_quant[,2] - col_quant[, 1]))
  
  # set upper lower bounds
  scaled_expression_mat[scaled_expression_mat < 0] <- 0
  scaled_expression_mat[scaled_expression_mat > 1] <- 1
  
  return(scaled_expression_mat)
  
}




collect_MEM_scores <- function(MEM_values, 
                               display_thresh = 3,
                               output_dir,
                               filename = "MEM_scores"
){
  
  if(missingArg(output_dir)) stop("output_dir is missing", call. = F)

  heatmap_data <- (MEM_values[[5]])[[1]]
  
  MEM_vals_scale <- as.matrix(round(heatmap_data, 0))
  
  new_rownames <- cytoMEM:::create.labels(MEM_vals_scale, 
                                          display_thresh, 
                                          heatmap_data)
  
  matrix.test <- as.matrix(new_rownames)
  
  rownames(matrix.test) <- 1:nrow(matrix.test)
  
  matrix.test <- matrix.test[match(rownames(matrix.test), 
                                   str_extract(matrix.test, "^\\d+")),]
  
  matrix.test <- as.matrix(matrix.test)

  savefilename <- paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),
                         "_", filename,".txt")
  
  
  apply(matrix.test, 1, function(x){ 
    
    cat("\n", x, sep = "\n", append = TRUE, 
        file = file.path(output_dir, savefilename))
    
    })

}




