# Subset cell-type markers -----------------------------------------------------
markers = getmarkers(c("immune","tumour","stroma"), unique = T)

# Remove CD56, CD8, P2Y12R 
markers = markers[!markers %in% c("CD56","CD8", "P2Y12R")]

# Subset cell-state markers ----------------------------------------------------

markers = getmarkers(names(spe@metadata$markers$cell_states))

# Subset marker expression -----------------------------------------------------

marker_expression = .scale_expression(
    exprs_matrix =  t(assay(spe, "zscore"))[, markers], scale_fun = "zscore"
    )

marker_expression = t(marker_expression)


# Plotting marker expression distribution --------------------------------------
# io$plots$expression_plots = purrr::map(markers, ~{
#     
#     plot_marker_exprs(spe, assay = "minmax", .x, 
#                       quantile_cuttoff = 0.5, 
#                       save_plot = F, return_plot = T)
# })
# names(io$plots$expression_plots) = markers
# 
# pdf(.new_file("marker_expression", path = io$output$data_out_dir), onefile = T)
# print(io$plots$expression_plots)
# dev.off()



# Define the global marker thresholds ------------------------------------------
.return_quantiles = function(marker, quantile, all_thresholds){
    
    all_thresholds[paste0(quantile*100,"%"), marker]
    
}

set_base_thresholds = function(expression, markers, default_quantile = 0.75, 
                               default_rank_quantile = 0.5){
    
    stopifnot(exprs = {
        
        all(markers %in% colnames(expression))
        
    })
    
    thresholds = list()
    
    thresholds$marker_expression = expression
    
    
    thresholds$maker_ranks = apply(thresholds$marker_expression, 
                                   MARGIN = 1, 
                                   FUN = rank, 
                                   ties.method = "min")
    
    thresholds$maker_ranks = as.data.frame(t(thresholds$maker_ranks))
    
    
    thresholds$quantiles = apply(X = thresholds$marker_expression,
                                 MARGIN = 2,
                                 FUN =  quantile, 
                                 probs = seq(0.01, 1, by = 0.01)
                                 )
    
    thresholds$markers = vector(mode = "list", length(markers)) |> setNames(markers)
    
    
    default_rank = floor(default_rank_quantile * ncol(thresholds$maker_ranks))
        
    thresholds$markers = purrr::map(thresholds$markers, ~{.x = list(quantile = default_quantile, 
                                                                    threshold = NA,
                                                                    rank = default_rank)})
    
    thresholds$markers = purrr::imap(thresholds$markers, ~{
        
        .x[["threshold"]] = .return_quantiles(marker = .y, 
                                             quantile = .x[["quantile"]], 
                                             all_thresholds = thresholds$quantiles)
        
        return(.x)
        
    })
    
    return(thresholds)
    
}

update_thresholds = function(thresholds, marker, update_quantile){
    
    stopifnot(exprs = {
        
        marker %in% names(thresholds[["markers"]])
        update_quantile > 0 
        update_quantile < 1 
        
    })
    
    current_threshold = thresholds[["markers"]][[marker]]
    
    updated_threshold = current_threshold
    
    updated_threshold[["quantile"]] = update_quantile
    
    updated_threshold[["threshold"]] = .return_quantiles(
        marker = marker, 
        quantile = update_quantile,
        all_thresholds = thresholds[["quantiles"]]
    )
    
    cli::cli_alert_info("Current {marker} threshold: {round(current_threshold$threshold, 3)} (quantile: {current_threshold$quantile}) ")
    
    cli::cli_alert("Updated {marker} threshold: {round(updated_threshold$threshold, 3)} (quantile: {updated_threshold$quantile}) ")
    
    thresholds[["markers"]][[marker]] = updated_threshold
    
    return(thresholds)
    
}

get_thresholds = function(gating_info, type = c("quantile", "rank")){
    
    out_thresholds = switch(EXPR = match.arg(type, several.ok = FALSE),
                            quantile = {vapply(gating_info[["markers"]], 
                                               \(x) x[["threshold"]],
                                               vector("numeric", 1))},
                            
                            rank = {vapply(gating_info[["markers"]], 
                                           \(x) x[["rank"]],
                                           vector("numeric", 1))}
    )
    
    return(out_thresholds)
    
} 

gating_info = set_base_thresholds(marker_expression, markers)

rm(marker_expression)


# Update thresholds ----
gating_info = update_thresholds(gating_info, "CD45", 0.5)
gating_info = update_thresholds(gating_info, "CD31", 0.8)
gating_info = update_thresholds(gating_info, "SMA", 0.8)

gating_info$pos_exprs = sweep(gating_info$marker_expression, 
                              MARGIN = 2, 
                              get_thresholds(gating_info),
                              FUN =  `>=`)

labels = list()
labels$undefined = rownames(gating_info$marker_expression)

# Positive thresholding --------------------------------------------------------
identify_cells = function(expression_matrix, markers, precentage_matching = NULL){
    
    total_markers = ncol(expression_matrix[,markers,drop = F])
    
    if(!is.null(precentage_matching)){
        
      cols_matching =  ceiling(total_markers * precentage_matching)
      
      cli::cli_alert_info("showing hits in {cols_matching}/{total_markers} marker")
    
      } else cols_matching = total_markers
    
    identified_cells = which(rowSums(expression_matrix[,markers, drop = F]) == cols_matching)
    
    rownames(expression_matrix)[identified_cells]
}

compare_vectors = function(label_list, vector1, vector2, return_common = TRUE) {
    
    if(length(attr(label_list, "names")) == 0) cli::cli_abort("label_list vectors are not named")
    if(!vector1 %in% names(label_list)) cli::cli_abort("vector1 not in label_list")
    if(!vector2 %in% names(label_list)) cli::cli_abort("vector2 not in label_list")
    
    intersect_name = unique(c(unlist(str_split(vector1, "_")),
                              unlist(str_split(vector2, "_"))))
    
    intersect_name = paste0(intersect_name, collapse = "_")
    
    intersect_markers = intersect(label_list[[vector1]], label_list[[vector2]])
    updated_vector1_markers = setdiff(label_list[[vector1]], intersect_markers)
    updated_vector2_markers = setdiff(label_list[[vector2]], intersect_markers)
    
    outlist = label_list
    
    if(return_common) outlist[[intersect_name]] = intersect_markers 
    outlist[[vector1]] = updated_vector1_markers
    outlist[[vector2]] = updated_vector2_markers
    
    return(outlist) 
    
}

update_undefined = function(label_list, 
                            all_cells = rownames(gating_info$marker_expression), 
                            undefined_name = "undefined", 
                            exclude_cells = NULL){
    
    if(length(attr(label_list, "names")) == 0) cli::cli_abort("label_list vectors are not named")
    if(!undefined_name %in% names(label_list)) cli::cli_abort("undefined_name not in label_list")
    if(!all(exclude_cells %in% names(label_list))) cli::cli_abort("all exclude cells not in label_list")
    
    
    outlist = label_list
    
    if(is.null(exclude_cells)) exclude_cells = undefined_name else exclude_cells = c(undefined_name, exclude_cells)
    
    updated_undefined_markers = setdiff(
            x = all_cells,
            y = unlist(outlist[which(!names(outlist) %in% exclude_cells)])
        )
        
    outlist[["undefined"]] = updated_undefined_markers
    
    return(outlist)
        
}

nest_labels = function(to_nest_list, nest_items, nest_name){
    
    out_list = to_nest_list
    
    nested_item = list(to_nest_list[nest_items])
    
    out_list = out_list[which(!names(out_list) %in% nest_items)]
    
    out_list[[nest_name]] = nested_item[[1]]
    
    return(out_list)    
}


combine_labels = function(label_list, combine_elements, combined_name){
    
    stopifnot(exprs = {
        
        !missing(label_list)
        !missing(combine_elements)
        !missing(combined_name)
        
        all(combine_elements %in% names(label_list))
        
        !combined_name %in% names(label_list)
        
    })
    
    combined = unlist(unlist(label_list[combine_elements], use.names = FALSE))
    
    outlist = label_list[which(!names(label_list) %in% combine_elements)]
    
    outlist[[combined_name]] = combined
    
    return(outlist)
    
}


check_distinct_vectors <- function(my_list) {
    all_distinct <- TRUE
    
    for (i in seq_along(my_list)) {
        for (j in seq_along(my_list)) {
            if (i != j && length(intersect(my_list[[i]], my_list[[j]])) > 0) {
                all_distinct <- FALSE
                break
            }
        }
    }
    
    return(all_distinct)
}


# MAP LABELS BACK TO CELLS  ----------------------------------------------------

spe$manual_gating = "undefined"

for (label in seq_along(labels)){
    
    colData(spe)[labels[[label]],"manual_gating"] = names(labels)[label]

}

spe$manual_gating[spe$manual_gating == "undefined"] = NA

spe$manual_gating = factor(spe$manual_gating, 
       levels = c("T cell", "NK cell", "Macrophage","Microglia",
                  "AC", "MES","NPC","OPC",
                  "Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial")
       )


table(spe$manual_gating)

# HEATMAP ----------------------------------------------------------------------

# Create the Heatmap Data
heatmap_data = anno_expression(object = spe, 
                               assay = "zscore", 
                               anno_label = "manual_gating",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "none",
                               # scale_colors = c("#A50026","white","#313695"),
                               scale_colors = c("grey99","grey99","#313695"),
                               # scale_colors = c("grey99","grey99","darkred"),
                               # scale_colors = c("black", "black", "#FCFFA4FF"),
                               heatmap_scale_across = "features"
                               )

set_heatmap_global_opts(
    row_title_fill = heatmap_data$plot_opts$row_title_bkg_colours,
    row_title_text = heatmap_data$plot_opts$row_title_text_colours
)



heatmap_data$ht = Heatmap(
    matrix = heatmap_data$heatmap_data, 
    col = heatmap_data$plot_opts$cols,
    rect_gp = gpar(col = "grey75", lwd = 0.70),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    heatmap_legend_param = list(
        title = heatmap_data$plot_opts$scale_title,
        at = heatmap_data$plot_opts$scale_breaks,
        legend_height = unit(25, "cm"),
        grid_width = unit(12.5, "mm")
    ),
    
    top_annotation = heatmap_data$column_anno,
    
    row_split = heatmap_data$plot_opts$rowsplits,
    
    row_names_side = "left",
    row_title_side = "left",
    row_names_max_width = unit(500, "mm"),
    row_gap = unit(5, "mm"),
    column_names_centered = TRUE
)



# svglite::svglite(
#   filename = .new_file("Manual_gating", io$output$data_out_dir, ext = "svg"),
#   width = 27.5, height = 27.5
# )
# 
# draw(heatmap_data$ht, 
#      padding = unit(c(20, 5, 5, 20), "mm"),    # bottom, left, top, right
#      heatmap_legend_side = "right")
# 
# dev.off()
# 



pdf(.new_file("Manual_gating", io$output$data_out_dir),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 10), "mm"),
     heatmap_legend_side = "right")

dev.off()


# END --------------------------------------------------------------------------