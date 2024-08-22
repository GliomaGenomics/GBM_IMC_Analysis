# PACKAGES ---------------------------------------------------------------------
# library(kohonen)
# library(Rphenograph)
# library(cytoMEM)
library(CATALYST)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(viridis)
library(dittoSeq)
library(ComplexHeatmap)


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
        current_anno = "outputs/cell_phenotyping/SOM_Phenograph/cell_by_cell",
        mem = "outputs/cell_phenotyping/mem"
    ),
    plots = list(batches = list())
)

# Load in utility functions
source(file.path(io$inputs$functions, "anno_expression.R"))


# SINGLE CELL DATA LOAD --------------------------------------------------------
spe = list.files(io$inputs$comp_data, pattern = "spe_comp_\\d{4}", full.names = T)
spe = readRDS(sort(spe,decreasing = T)[1])

# spe = readRDS(list.files(io$inputs$comp_data, pattern = "comp_raw", full.names = T))

# FILTER HIGH EXPRESSING TYPE MARKERS ------------------------------------------
markers = getmarkers(marker_list = spe@metadata$markers, 
                     anno_level = 2, 
                     pheno = c("Immune", "Tumour", "Stroma"), 
                     unique = TRUE)

# Identify cells which have high expression across multiple markers
loco_cells = function(expr_matrix, # rows = cells and columns = markers
                      quant_threshold = 0.75, 
                      hig_exprs_threshold = 0.50){
    
    thresholds = apply(expr_matrix, 2,\(marker) quantile(marker, quant_threshold))
    
    # Identify cells with high expression across multiple markers using varying thresholds
    high_exprs_cells = apply(expr_matrix, 1, \(cells) sum(cells >= thresholds)/ncol(expr_matrix))
    
    
    high_exprs_cells = high_exprs_cells[which(high_exprs_cells >= hig_exprs_threshold)]
    
    return(high_exprs_cells)
    
}


high_exprs_cells = loco_cells(expr_matrix = t(assay(spe, "exprs")),
                              quant_threshold = 0.9,
                              hig_exprs_threshold = 0.5)

spe = spe[,!colnames(spe) %in% names(high_exprs_cells)]

rm(high_exprs_cells, loco_cells)

# NORMALIZATION ----------------------------------------------------------------
spe = clean_expression(object = spe, 
                       assay = "exprs", 
                       scale_fun = "minmax", 
                       clip_exprs_vals = T,
                       clip_quantile = 0.99)


# DIMENSION REDUCTION (UMAP) ---------------------------------------------------
source(file.path(io$inputs$functions, "dim_reductions.R"))

spe = runUMAP(object = spe, 
              reduct_matrix = t(assay(spe, "minmax")),
              reducedDim_suffix = "minmax",
              n_neighbors = 35)

# HARMONY BATCH CORRECTION -----------------------------------------------------
source(file.path(io$inputs$functions, "harmony_correction.R"))

spe = harmony_correction(object = spe, 
                         assay = "minmax",
                         correct_vars = "patient", 
                         pca_before_correction = TRUE, 
                         npcs = 30, 
                         dim_reduction = TRUE, n_neighbors = 35)

# VISUALISING BATCH EFFECTS ----------------------------------------------------
source(file.path(io$inputs$functions, "anno_dimplot.R"))

spe@metadata$batches = spe@metadata$batches[spe@metadata$batches != "region_type_new"]

io$plots$batches$before = .plot_batches(object = spe, 
                                        batches = spe@metadata$batches,
                                        reduction = "UMAP_minmax")

io$plots$batches$after = .plot_batches(object = spe,
                                       batches = spe@metadata$batches,
                                       reduction = "UMAP_harmony_minmax")

io$plots$batches$combined <- purrr::map2(
    io$plots$batches$before$UMAP_minmax, 
    io$plots$batches$after$UMAP_harmony_minmax, ~{
        
        combined <- .x + .y + patchwork::plot_layout(guides = "collect") 
        
        })


pdf(.new_file("Harmony_batch_correction", path = io$output$som_pheno), 
    onefile = T, width = 35, height = 15)

io$plots$batches$combined

dev.off()

# NOTE TO FUTURE SELF ----------------------------------------------------------

# after visualising the batch effects and completing the batch correction, the below
# steps were initally completed but these did not yeild distinct enough cell type
# labels, therefore a more sytematic gating method was employed to label the various
# cell type labels


# At this point the next logical step in the pipline was that detailed in the script
# named "5b_Gating_based_Cell_Annotation.R" and the logic contained within:
# "5c_Expression_Ranks_Gating_Logic.R"



# METACLUSTERING CORRECTED CELLS -----------------------------------------------
source(file.path(io$inputs$functions, "Flowsom_metaclustering.R"))

spe = som_metacluster(cluster_matrix = reducedDim(spe, "harmony_minmax"),
                      object = spe,
                      som_xdim = 100, som_ydim = 100, 
                      metacluster_k = 15)

# Running just the phenograph clustering 
adjust_pheno_K = 6

adjust_pheno_K = metaclustering(
    cluster_medians = spe@metadata$SOM_metaclusters$som_codes$som_cluster_medians,
    cell_cluster_ids = spe@metadata$SOM_metaclusters$som_codes$som_cluster_ids,
    cluster_col = "SOM_cluster",
    k = adjust_pheno_K)

spe@metadata$SOM_metaclusters$inital_som_metaclusters = adjust_pheno_K
spe$initial_metacluster = factor(spe@metadata$SOM_metaclusters$inital_som_metaclusters$metacluster)

rm(adjust_pheno_K)

# VIZUALISE CLUSTER MEMBERSHIP -------------------------------------------------
visualise_clusters(object = spe, 
                   anno_label = "initial_metacluster", 
                   reduction = "UMAP_harmony_minmax",
                   dot_size = 0.5,
                   opacity = 0.7,
                   show_labels = TRUE, 
                   return_plots = FALSE, 
                   save_plots = TRUE,
                   outfile = .new_file("Initial_metaclusters",io$output$som_pheno)
                   )
# CLUSTER EXPRESSION HEATMAP ---------------------------------------------------

# Create the Heatmap Data
heatmap_data = anno_expression(object = spe, 
                               assay = "exprs", 
                               anno_label = "initial_metacluster",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "minmax",
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


pdf(.new_file("Initial_metaclusters", io$output$som_pheno),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "right")

dev.off()

# MARKER ENRICHMENT MODELLING --------------------------------------------------
source(file.path(io$inputs$functions, "gather_MEM_scores.R"))

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
    threshold = 3)

write.csv(x = heatmap_data$MEM_scores, row.names = FALSE,
          file = .new_file("Initial_metaclusters_MEM",
                           io$output$som_pheno, ext = "csv")
          )

# CLUSTER EXPRESSION DOTPLOT ---------------------------------------------------
source(file.path(io$inputs$functions, "cluster_dotplots.R"))

pdf(.new_file("Initial_metaclusters", io$output$som_pheno),
    width = 90, height = 50,
    onefile = T)

cluster_dotplot(spe = spe, 
                assay = "minmax",
                summary_fun = "median", 
                markers = rev(markers), 
                cluster_col = "initial_metacluster",
                low_color = "grey90",
                high_color = "darkred",
                low_color_pct = 0.5,
                base_font_size = 50)

dev.off()

# RE-CODE THE INITIAL CLUSTERS -------------------------------------------------
cluster_labels = sort(list.files(io$output$som_pheno, pattern = "MEM", 
                                 full.names = TRUE), decreasing = TRUE)

cluster_labels = read.csv(file = cluster_labels[1], 
                       colClasses =  rep("character", 4))

cluster_labels = cluster_labels[c("cluster", "label")]

cluster_labels$label[cluster_labels$label == ""] = "undefined"

spe$initial_anno = factor(cluster_labels$label[match(spe$initial_metacluster, 
                                                     cluster_labels$cluster)])

table(spe$initial_anno)


compare_anno = function(spe_object, compare_anno, compare_label, lookup_anno){
    
    if(!all(c(compare_anno,lookup_anno) %in% names(colData(spe_object)))) cli::cli_abort("one or more annotations not found in colData(spe_object)")
    
    if(!all(is.factor(spe[[compare_anno]]) & is.factor(spe[[lookup_anno]]))) cli::cli_abort("one or more annotations are not factors")
    
    if(!compare_label %in% levels(spe_object[[compare_anno]])) cli::cli_abort("anno label not in levels({compare_anno})")
    
    cli::cli_inform("Comparing {compare_label} labels:\n")
    
    table(spe_object[[lookup_anno]][spe_object[[compare_anno]] == compare_label])
    
}

compare_anno(spe, 
             compare_anno = "cell_anno_v2",
             compare_label = "Macrophage",
             lookup_anno =  "cell_anno")

# SAVE INITIAL CLUSTERS SPATIAL OBJECT -----------------------------------------
saveRDS(spe, .new_file("spe_comp", path = io$inputs$comp_data, ext = "rds"))
# END --------------------------------------------------------------------------