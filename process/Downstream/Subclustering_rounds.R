# PACKAGES ---------------------------------------------------------------------
library(kohonen)
library(Rphenograph)
library(CATALYST)
library(cytoMEM)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(viridis)
library(dittoSeq)

# I/O --------------------------------------------------------------------------
io <- list(
    
    inputs = list(
        functions = "process/Downstream/functions",
        comp_data = "data/downstream/compensated",
        dim_reductions = "data/downstream/dim_reductions"
    ),
    
    output= list(
        data_out_dir = "outputs/cell_phenotyping",
        batches = "outputs/QC/batch",
        som_pheno = "outputs/cell_phenotyping/SOM_Phenograph",
        mem = "outputs/cell_phenotyping/mem"
    ),
    plots = list(batches = list())
)

# LOAD CURRENT SINGLE CELL DATA LOAD -------------------------------------------
spe = readRDS(file.path(io$inputs$comp_data, "spe_comp.rds"))

# LOAD IN REQUIRED FUNCTIONS ---------------------------------------------------
source(file.path(io$inputs$functions, "utils.R"))
source(file.path(io$inputs$functions, "normalisation.R"))
source(file.path(io$inputs$functions, "dim_reductions.R"))
source(file.path(io$inputs$functions, "anno_dimplot.R"))
source(file.path(io$inputs$functions, "Flowsom_metaclustering.R"))
source(file.path(io$inputs$functions, "get_markers.R"))
source(file.path(io$inputs$functions, "anno_expression.R"))
source(file.path(io$inputs$functions, "marker_expression_histogram.R"))
source(file.path(io$inputs$functions, "gather_MEM_scores.R"))

# SET CLUSTER ROUND PARAMS -----------------------------------------------------
cluster_round = "round_1"

somxy = 100
pheno_K = 10

# Create new round dir (if it does not exist)
io$output$cluster_round_out = file.path(io$output$som_pheno, cluster_round)
if(!dir.exists(io$output$cluster_round_out)) dir.create(io$output$cluster_round_out)

# Set the Markers
markers = getmarkers(marker_list = spe@metadata$markers, anno_level = 2,
                     pheno = c("Immune", "Tumour", "Stroma"),
                     unique = TRUE)

# SUB-CLUSTER ROUND START ----

# Filter out the unknown cells
filt_spe = spe[markers, spe$cell_anno_v2 == "unknown"]
# re-level the markers 
filt_spe$cell_anno_v2 = factor(filt_spe$cell_anno_v2)

# Rescale the unknown cells
filt_spe = clean_expression(object = filt_spe, 
                            assay = "exprs", 
                            scale_fun = "minmax", 
                            clip_exprs_vals = FALSE)

# THIS IS OPTIONAL AND WILL BE DETERMINED BY CLUSTERS ----
# HAVING EXPRESSION ACROSSS ALL CELLS 

# identify cells which have high expression across multiple markers
loco_cells = function(expr_matrix, # rows = cells and columns = markers
                      quant_threshold = 0.75, 
                      hig_exprs_threshold = 0.50){
    
    thresholds = apply(expr_matrix, 2,\(marker) quantile(marker, quant_threshold))
    
    # Identify cells with high expression across multiple markers using varying thresholds
    high_exprs_cells = apply(expr_matrix, 1, \(cells) sum(cells >= thresholds)/ncol(expr_matrix))
    
    
    high_exprs_cells = high_exprs_cells[which(high_exprs_cells >= hig_exprs_threshold)]
    
    return(high_exprs_cells)
    
}

# high_exprs_cells = loco_cells(expr_matrix = t(assay(filt_spe, "minmax")),
#                               quant_threshold = 0.90, 
#                               hig_exprs_threshold = 0.70)
# 
# filt_spe = filt_spe[,!colnames(filt_spe) %in% names(high_exprs_cells)]
# 
# 
# # Rescale the unknown cells
# filt_spe = clean_expression(object = filt_spe, 
#                             assay = "exprs", 
#                             scale_fun = "minmax", 
#                             clip_exprs_vals = FALSE)

# # Look at CD45 marker expression in the unknown immune cells
# plot_marker_exprs(filt_spe, assay = "minmax", 
#                   marker = "CD45", 
#                   quantile_cuttoff = 0.5,
#                   save_plot = TRUE,
#                   outdir = io$output$cluster_round_out,
#                   filename = paste0(cluster_round, "_CD45_expression")
# )
# 
# # Filter the cells which have CD45 expression below CD45 threshold = 0.299
# filt_spe = filt_spe[,t(assay(filt_spe, "minmax"))[,"CD45"] >= 0.319]
# 
# # Rescale the unknown cells
# filt_spe = clean_expression(object = filt_spe, 
#                             assay = "exprs", 
#                             scale_fun = "minmax", 
#                             clip_exprs_vals = FALSE)
# 
# # Update the UMAP visualization using the batch-corrected embeddings
# filt_spe = runUMAP(object = filt_spe, 
#                    reduct_matrix = reducedDim(filt_spe, "harmony_minmax"),
#                    reducedDim_suffix = "harmony_minmax",
#                    n_neighbors = 30)
# 
# ----

pdf(.new_file("Immune_clustering_batches", path = io$output$cluster_round_out), 
    onefile = T, width = 15, height = 10)

.plot_batches(object = filt_spe,
              batches = filt_spe@metadata$batches,
              reduction = "UMAP_harmony_minmax")

dev.off()


# Metacluster the Immune filtered cells
filt_spe = som_metacluster(cluster_matrix = reducedDim(filt_spe, "harmony_minmax"),
                           object = filt_spe,
                           som_xdim = somxy, som_ydim = somxy, 
                           metacluster_k = pheno_K)

filt_spe$SOM_cluster = factor(filt_spe@metadata$SOM_metaclusters$som_metacodes$SOM_cluster)
filt_spe$metacluster = factor(filt_spe@metadata$SOM_metaclusters$som_metacodes$metacluster)


visualise_clusters(object = filt_spe, 
                   anno_label = "metacluster", 
                   reduction = "UMAP_harmony_minmax",
                   dot_size = 0.5,
                   opacity = 0.7,
                   show_labels = TRUE, 
                   return_plots = FALSE, 
                   save_plots = TRUE,
                   outfile = .new_file("Immune_clusters",io$output$cluster_round_out)
                   )


# Create the Heatmap Data
heatmap_data = anno_expression(object = filt_spe, 
                               assay = "exprs", 
                               anno_label = "metacluster",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "minmax",
                               heatmap_scale_across = "features")

set_heatmap_global_opts(
    row_title_fill = heatmap_data$plot_opts$row_title_bkg_colours,
    row_title_text = heatmap_data$plot_opts$row_title_text_colours)


if("Pan-Immune" %in% unique(heatmap_data$plot_opts$rowsplits)){
    
    heatmap_data$plot_opts$rowsplits= recode(heatmap_data$plot_opts$rowsplits,
                                             "Pan-Immune" = "Immune")  
}

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


pdf(.new_file("Immune_cluster_exprs", io$output$cluster_round_out),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "right")

dev.off()


heatmap_data$MEM_values <- MEM(
    exp_data = cbind(t(heatmap_data$heatmap_data), 
                     cluster = as.numeric(colnames(heatmap_data$heatmap_data))),
    transform = FALSE, 
    choose.markers = FALSE, # Change choose.markers to TRUE to see and select channels in console
    markers = "all",
    choose.ref = FALSE,
    zero.ref = FALSE,
    rename.markers = FALSE, # Change rename.markers to TRUE to see and choose new names for channels in console
    new.marker.names = "none",
    IQR.thresh = NULL)


heatmap_data$MEM_scores = gather_MEM_scores(
    MEM_score_matrix = heatmap_data$MEM_values$MEM_matrix[[1]],
    threshold = 1)

write.csv(x = heatmap_data$MEM_scores, 
          file = .new_file("Immune_cluster_MEM", io$output$cluster_round_out, ext = "csv"), 
          row.names = FALSE)

saveRDS(filt_spe, .new_file(cluster_round, io$inputs$comp_data, ext = "rds"))


# SUB-CLUSTER ROUND END ----

# Varying the K for Rphenograph: 
new_k_pheno = metaclustering(cluster_medians = filt_spe@metadata$SOM_metaclusters$som_codes$som_cluster_medians,
                             cell_cluster_ids = filt_spe@metadata$SOM_metaclusters$som_codes$som_cluster_ids,
                             cluster_col = "SOM_cluster",
                             k = 10)

filt_spe$metacluster = factor(new_k_pheno $metacluster)


markers = getmarkers(anno_level = 3, unique = T,
                     pheno =  c("Pan-Immune", "Lymphoid", "Myeloid", "Endothelial"))

markers = markers[markers != "CD45"]


# ASSIGN CLUSTER LABELS FOR IMMUNE CELLS ----
cell_labels = sort(list.files(io$output$cluster_round_out, 
                              pattern = "MEM", full.names = TRUE),
                   decreasing = TRUE)

cell_labels = read.csv(file = cell_labels[1], 
                       colClasses =  rep("character", 4))

cell_labels = cell_labels[c("cluster", "label")]

filt_spe$metacluster = factor(cell_labels$label[match(filt_spe$metacluster, cell_labels$cluster)])

# SUBCLUSTER UNDEFINED IMMUNE CELLS ----
filt_sub_spe = filt_spe[,filt_spe$metacluster == "SUB"]

# Metacluster the Immune filtered cells
filt_sub_spe = som_metacluster(
    cluster_matrix = reducedDim(filt_sub_spe, "harmony_minmax"),
    object = filt_sub_spe,
    som_xdim = 50, som_ydim = 50, 
    metacluster_k = 10)

filt_sub_spe$SOM_cluster = factor(filt_sub_spe@metadata$SOM_metaclusters$som_metacodes$SOM_cluster)
filt_sub_spe$metacluster = factor(filt_sub_spe@metadata$SOM_metaclusters$som_metacodes$metacluster)


visualise_clusters(object = filt_sub_spe, 
                   anno_label = "metacluster", 
                   reduction = "UMAP_harmony_minmax",
                   dot_size = 0.5,
                   opacity = 0.7,
                   show_labels = TRUE, 
                   return_plots = FALSE, 
                   save_plots = TRUE,
                   outfile = .new_file("Immune_sub_clusters",io$output$cluster_round_out)
                   )

# Create the Heatmap Data
heatmap_data = anno_expression(object = filt_sub_spe, 
                               assay = "exprs", 
                               anno_label = "metacluster",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "minmax",
                               heatmap_scale_across = "features")

set_heatmap_global_opts(
    row_title_fill = heatmap_data$plot_opts$row_title_bkg_colours,
    row_title_text = heatmap_data$plot_opts$row_title_text_colours)


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


pdf(.new_file("Immune_sub_cluster_exprs", io$output$cluster_round_out),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "right")

dev.off()


heatmap_data$MEM_values <- MEM(
    exp_data = cbind(t(heatmap_data$heatmap_data), 
                     cluster = as.numeric(colnames(heatmap_data$heatmap_data))),
    transform = FALSE, 
    choose.markers = FALSE, # Change choose.markers to TRUE to see and select channels in console
    markers = "all",
    choose.ref = FALSE,
    zero.ref = FALSE,
    rename.markers = FALSE, # Change rename.markers to TRUE to see and choose new names for channels in console
    new.marker.names = "none",
    IQR.thresh = NULL)


heatmap_data$MEM_scores = gather_MEM_scores(
    MEM_score_matrix = heatmap_data$MEM_values$MEM_matrix[[1]],
    threshold = 1)

write.csv(x = heatmap_data$MEM_scores, 
          file = .new_file("Immune_sub_cluster_MEM", io$output$cluster_round_out, ext = "csv"), 
          row.names = FALSE)

# ASSIGN SUB CLUSTER LABELS FOR IMMUNE CELLS ----
cell_labels = sort(list.files(io$output$cluster_round_out, 
                              pattern = "sub_cluster_MEM", full.names = TRUE),
                   decreasing = TRUE)

cell_labels = read.csv(file = cell_labels[1], 
                       colClasses =  rep("character", 4))

cell_labels = cell_labels[c("cluster", "label")]

filt_sub_spe$metacluster = factor(cell_labels$label[match(filt_sub_spe$metacluster, 
                                                          cell_labels$cluster)])


filt_spe$metacluster[rownames(colData(filt_spe)) %in% rownames(colData(filt_sub_spe))] = as.character(filt_sub_spe$metacluster)

filt_spe$cell_anno_v2 = factor(filt_spe$metacluster)


cell_labels = as.data.frame(colData(filt_spe)["cell_anno_v2"])
cell_labels = tibble::rownames_to_column(cell_labels, var = "cell_id")
cell_labels = dplyr::rename(cell_labels, filt_anno = cell_anno_v2)
cell_labels$filt_anno = as.character(cell_labels$filt_anno)


main_anno = as.data.frame(colData(spe)["cell_anno_v2"])
main_anno = tibble::rownames_to_column(main_anno, var = "cell_id")
main_anno$cell_anno_v2 = as.character(main_anno$cell_anno_v2)

combined_anno = left_join(main_anno, cell_labels, by = "cell_id")

filt_anno_index = !is.na(combined_anno$filt_anno)

combined_anno$cell_anno_v2[filt_anno_index] = combined_anno$filt_anno[filt_anno_index]

all(rownames(colData(spe)) == combined_anno$cell_id)

spe$cell_anno_v2 = factor(combined_anno$cell_anno_v2)


saveRDS(spe, file.path(io$inputs$comp_data, "spe_comp.rds"))


# END ----

