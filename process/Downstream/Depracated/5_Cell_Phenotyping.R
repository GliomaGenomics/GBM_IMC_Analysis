# The strategy for annotation involves first identifying three general groups:
#  CD45+ (Immune); CD31+, SMA+ (Endothelial); and all other remaining cells.
#  That way i should be able to resolve the immune, malignant/normal brain fractions
#  without the signal causing interference, especially when clustering.

# PACKAGES ---------------------------------------------------------------------
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(viridis)
library(dittoSeq)
library(CATALYST)

# I/O --------------------------------------------------------------------------
io <- list(
    
    inputs = list(
        functions = "process/Downstream/functions",
        comp_data = "data/downstream/compensated",
        dim_reductions = "data/downstream/dim_reductions"
    ),
    
    output= list(
        data_out_dir = "outputs/cell_phenotyping",
        manual_gating = "outputs/cell_phenotyping/manual",
        phenograph = "outputs/cell_phenotyping/phenograph",
        som = "outputs/cell_phenotyping/som",
        som_immune = "outputs/cell_phenotyping/som/immune",
        mem = "outputs/cell_phenotyping/mem"
    )
)

if(!dir.exists(io$output$data_out)) dir.create(io$output$data_out)

# Create QC out sub directories 
if(!dir.exists(io$output$manual_gating)) dir.create(io$output$manual_gating)
if(!dir.exists(io$output$phenograph)) dir.create(io$output$phenograph)
if(!dir.exists(io$output$som)) dir.create(io$output$som)
if(!dir.exists(io$output$som_immune)) dir.create(io$output$som_immune)
# if(!dir.exists(io$output$)) dir.create(io$output$)

# READ IN COMPENSATED SINGLE CELL DATA  ----------------------------------------
spe <- readRDS(file.path(io$inputs$comp_data, "spe_comp.rds"))

# VISUALISING MARKER EXPRESSION (DIMENSIONALITY REDUCTION)  --------------------

plot_list <- multi_dittoDimPlot(spe, 
                                var = rowData(spe)$name[rowData(spe)$downstream == 1],
                                reduction.use = "pacMap_harmony",
                                assay = "exprs", 
                                size = 0.2, 
                                list.out = TRUE)

plot_list <- lapply(plot_list, function(x) x + viridis::scale_color_viridis())

marker_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)

pdf(file = file.path(io$output$manual_gating, 
                     "PacMap_Harmony_Marker_expressions.pdf"),
    width = 12, height = 40, onefile = T)

marker_plot

dev.off()

# SETTING MARKER INFO METADATA FOR CLUSTERING   --------------------------------

spe@metadata$markers = list()

spe@metadata$markers$immune <- c("CD45",                      # Pan Immune
                                 "CD3","CD8",                 # T cells
                                 "NKP46","CD56",              # NK cells
                                 "IBA1",                      # Macrophages
                                 "P2Y12R","TMEM119")          # Microglia

spe@metadata$markers$neoplastic <- c("SLC1A3_EAAT1", "HOPX",  # AC
                                     "SOD2","CHI3L1",         # MES
                                     "ANEXIN_A2", "ANXA1",    # MES-hypoxia
                                     "DLL3", "BCAN",          # NPC
                                     "SCD5", "OLIG1")         # OPC

spe@metadata$markers$endothelial  <- c("CD31", "SMA")

spe@metadata$markers$stroma  <- c("NeuN_FOX3", "CD56",        # Neuron
                                  "GFAP",                     # Astrocyte
                                  "MOG")                      # Oligodendrocyte

spe@metadata$markers$all <- unlist(spe@metadata$markers)


# INITIAL PHENOGRAPH CLUSTERING (ALL MARKERS)  ---------------------------------

# The clustering will be run using the batch corrected harmony embeddings
phenograph_clusters <- Rphenograph::Rphenograph(
    reducedDim(spe, "harmony"), k = 45)

phenograph_clusters <- factor(igraph::membership(phenograph_clusters[[2]]))

spe$pg_clusters <- phenograph_clusters


spe@metadata$color_vectors$pg_clusters <- setNames(
    
    dittoSeq:dittoColors(1)[1:length(levels(spe$pg_clusters))], levels(spe$pg_clusters)
    
)


# Set the options for the ggrepl labels
options(ggrepel.max.overlaps = Inf,
        min.segment.length = 0.1)


pdf(file = file.path(io$output$phenograph, "Initial_clustering.pdf"), 
    onefile = T, width = 15, height = 10)

dittoDimPlot(spe, var = "pg_clusters", 
             reduction.use = "pacMap_harmony", 
             size = 0.2,
             color.panel = spe@metadata$color_vectors$pg_clusters,
             do.label = TRUE) +
    ggtitle("Phenograph Clusters (Harmony Embeddings)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", hjust = 0.5)
    )

dev.off()


# function to plot the counts for each cluster
source(file.path(io$inputs$functions, "cluster_counts.R"))


pdf(file = file.path(io$output$phenograph, "Initial_clustering_counts.pdf"), 
    onefile = T, width = 9, height = 7)

cluster_cell_counts(spe, cluster_label = "pg_clusters", 
                    cluster_colors = spe@metadata$color_vectors$pg_clusters)

dev.off()

source(file.path(io$inputs$functions, "cluster_marker_expression.R"))

imap(spe@metadata$markers, ~{
    plot_cluster_expression(
        spe, 
        markers = .x,
        cluster_labels = "pg_clusters",
        fun = "mean",
        scale_fun = "zscore",
        scale = "row",
        plot_title = .y,
        outfile = file.path(io$output$phenograph,
                            glue::glue("{.y}_marker_expression.pdf"))
        )
})


# INITIAL SOM CLUSTERING (ALL MARKERS)  ----------------------------------------
library(bluster)
# Perform SOM clustering
som_clusters <- clusterRows(reducedDim(spe, "harmony"), SomParam(100), full = TRUE)

# Cluster the 100 SOM codes into larger clusters
ccp <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = t(som_clusters$objects$som$codes[[1]]),
    maxK = 45,
    reps = 100,
    pItem = 1.0,
    pFeature = 1.0, 
    clusterAlg = "km", 
    distance = "euclidean", 
    seed = 1234, 
    plot = NULL)

final_som_clusters <- lapply(2:length(ccp), function(x){
    
    ccp[[x]][["consensusClass"]][som_clusters$clusters]
    
})  

names(final_som_clusters) <- glue::glue("cluster_{1:length(final_som_clusters)+1}")

spe$som_clusters = factor(as.numeric(final_som_clusters[[length(final_som_clusters)]]))


spe@metadata$color_vectors$som_clusters <- setNames(
    
    viridis::turbo(length(levels(spe$som_clusters))), levels(spe$som_clusters)
    
)


pdf(file = file.path(io$output$som, "Initial_clustering.pdf"), 
    onefile = T, width = 15, height = 10)

dittoDimPlot(spe, var = "som_clusters", 
             reduction.use = "pacMap_harmony", 
             size = 0.2,
             color.panel = spe@metadata$color_vectors$som_clusters,
             do.label = TRUE) +
    ggtitle("SOM Clusters (Harmony Embeddings)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", hjust = 0.5)
    )

dev.off()


pdf(file = file.path(io$output$som, "Initial_clustering_counts.pdf"), 
    onefile = T, width = 9, height = 7)

cluster_cell_counts(spe, cluster_label = "som_clusters", 
                    cluster_colors = spe@metadata$color_vectors$som_clusters)

dev.off()


imap(spe@metadata$markers, ~{
    plot_cluster_expression(
        spe, 
        markers = .x,
        cluster_labels = "som_clusters",
        fun = "mean",
        scale_fun = "zscore",
        scale = "row",
        plot_title = .y,
        outfile = file.path(io$output$som,
                            glue::glue("{.y}_marker_expression.pdf"))
    )
})


# MARKER ENRICHMENT MODELLING (MEM) --------------------------------------------

library(cytoMEM)

expr = t(spe@assays@data$exprs)
expr = expr[,rownames(spe)[rowData(spe)$downstream ==1]]
expr <- cbind(expr, cluster = as.numeric(spe@colData[["som_clusters"]]))

# expr <- cbind(expr, cluster = as.numeric(spe@colData[["pg_clusters"]]))

MEM_values <- MEM(exp_data = expr,
                  transform = FALSE, 
                  choose.markers = FALSE, # Change choose.markers to TRUE to see and select channels in console
                  markers = "all",
                  choose.ref = FALSE,
                  zero.ref = FALSE,
                  rename.markers = FALSE, # Change rename.markers to TRUE to see and choose new names for channels in console
                  new.marker.names = "none",
                  IQR.thresh = NULL)


build_heatmaps(MEM_values, 
               cluster.MEM = "none", 
               cluster.medians = "none",
               cluster.IQRs = "none", 
               display.thresh = 3, 
               output.files = TRUE,
               labels = FALSE, 
               only.MEMheatmap = FALSE)


# After running the MEM scoring on the arcsinh transformed counts along with The
# SOM cluster labels a number of cluster were showing high expression of 
# endothelial cell markers:
#     
# 10 : ▲ SMA+5 ANEXIN_A2+3 
# 8  : ▲ SMA+4 CD31+3
# 11 : ▲ CD31+3
# 7  : ▲ SMA+6
# 
# Because cluster 10 was also showing high expression for marker that are being 
# used to delineate the neoplastic cells, 7,8 and 11 were labelled as endothelial
# cells

spe@colData$current_anno = ifelse(spe@colData$som_clusters %in% c("7","8","11"), "endothelial", "undefined")

table(spe@colData$current_anno)

# LABELLING THE IMMUNE CELLS ---------------------------------------------------

source(file.path(io$inputs$functions, "marker_expression_histogram.R"))

io$output$exprs_histograms = file.path(io$output$manual_gating, "expression_histograms")
if(!dir.exists(io$output$exprs_histograms)) dir.create(io$output$exprs_histograms)

threshold = function(spe, marker, quantile = 0.5, assay = "exprs"){
    
    quantile(t(assay(spe, assay))[,marker], quantile)
    
}

walk(spe@metadata$markers$immune, ~{
    
    plot_marker_exprs(spe_object = spe, 
                      marker = .x, 
                      quantile_cuttoff = 0.5,  
                      outdir = io$output$exprs_histograms)
})

plot_marker_exprs(spe_object = spe, marker = "CD45", 
                  quantile_cuttoff = 0.25,  outdir = io$output$exprs_histograms)


immune <- spe[,which(assay(spe, "exprs")["CD45",] >= threshold(spe, "CD45", 0.25))]
dim(immune)[2]

anno <- list()

# NK cells
anno$nk_cells = assay(immune, "exprs")["NKP46",] >= threshold(immune, "NKP46") | assay(immune, "exprs")["CD56",] >= threshold(immune, "CD56")
table(anno$nk_cells)

# T cells
anno$t_cells = assay(immune, "exprs")["CD3",] >= threshold(immune, "CD3") | assay(immune, "exprs")["CD8",] >= threshold(immune, "CD8")
table(anno$t_cells)

# Macrophages
anno$macro = assay(immune, "exprs")["IBA1",] >= threshold(immune, "IBA1")
table(anno$macro)


immune$immune_anno <- "undefined"

anno = lapply(anno, which)

immune$immune_anno[setdiff(setdiff(anno$nk_cells, anno$t_cells), anno$macro)] <- "NK cells"

immune$immune_anno[setdiff(setdiff(anno$t_cells, anno$nk_cells), anno$macro)] <- "T cells"

immune$immune_anno[setdiff(setdiff(anno$macro, anno$nk_cells), anno$t_cells)] <- "Macrophage"

table(immune@colData$immune_anno)


# Microglia
anno$microglia = assay(immune, "exprs")["TMEM119",] >= threshold(immune, "TMEM119")| assay(immune, "exprs")["P2Y12R",] >= threshold(immune, "P2Y12R")
table(anno$microglia)

immune$immune_anno[intersect(which(immune$immune_anno == "Macrophage"), which(anno$microglia))] <- "Microglia"

table(immune@colData$immune_anno)


all_cells <- as.data.frame(colData(spe)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname)

all_cells <- as.data.frame(colData(immune)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname, immune_anno) %>%
    left_join(all_cells, ., by = "rowname")

all_cells$immune_anno = ifelse(is.na(all_cells$immune_anno), 
                               "undefined", all_cells$immune_anno)

all(rownames(colData(spe)) == all_cells$rowname)

spe$cell_anno = all_cells$immune_anno 

table(spe$cell_anno)


# LABELLING THE ENDOTHELIAL CELLS ----------------------------------------------

walk(spe@metadata$markers$endothelial, ~{
    
    plot_marker_exprs(spe_object = spe, 
                      marker = .x, 
                      quantile_cuttoff = 0.5,  
                      outdir = io$output$exprs_histograms)
})

endothelial <- spe[,colData(spe)$cell_anno == "undefined"]
dim(endothelial)[2]

anno <- list()

anno$endothelial = assay(endothelial, "exprs")["CD31",] >= threshold(endothelial, "CD31", 0.8) & assay(endothelial, "exprs")["SMA",] >= threshold(endothelial, "SMA", 0.8) 
table(anno$endothelial)

endothelial$endo_anno <- "undefined"

endothelial$endo_anno[which(anno$endothelial)] = "Endothelial"
table(endothelial$endo_anno)


all_cells <- as.data.frame(colData(endothelial)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname, endo_anno) %>%
    left_join(all_cells, ., by = "rowname")

all_cells$endo_anno = ifelse(is.na(all_cells$endo_anno), 
                               "undefined", all_cells$endo_anno)

all(rownames(colData(spe)) == all_cells$rowname)

unique(all_cells$immune_anno[all_cells$endo_anno == "Endothelial"])

all_cells$immune_anno[all_cells$endo_anno == "Endothelial"] = "Endothelial"

table(all_cells$immune_anno)

spe$cell_anno = all_cells$immune_anno 

table(spe$cell_anno)

# Cleaning up some unused variables

rm(immune, endothelial)

# LABELLING THE NEOPLASTIC CELLS -----------------------------------------------

walk(spe@metadata$markers$neoplastic, ~{
    
    plot_marker_exprs(spe_object = spe, 
                      marker = .x, 
                      quantile_cuttoff = 0.5,  
                      outdir = io$output$exprs_histograms)
})

neoplastic <- spe[,colData(spe)$cell_anno == "undefined"]
dim(neoplastic)[2]

anno <- list()

anno$AC = assay(neoplastic, "exprs")["SLC1A3_EAAT1",] >= threshold(spe, "SLC1A3_EAAT1", 0.5) & assay(neoplastic, "exprs")["HOPX",] >= threshold(spe, "HOPX", 0.5) 
table(anno$AC)

anno$MES = assay(neoplastic, "exprs")["SOD2",] >= threshold(spe, "SOD2" , 0.5) & assay(neoplastic, "exprs")["CHI3L1",] >= threshold(spe, "CHI3L1", 0.5) 
table(anno$MES)

anno$MES_hypoxia = assay(neoplastic, "exprs")["ANXA1",] >= threshold(spe, "ANXA1", 0.5) & assay(neoplastic, "exprs")["ANEXIN_A2",] >= threshold(spe, "ANEXIN_A2", 0.75) 
table(anno$MES_hypoxia)

anno$NPC = assay(neoplastic, "exprs")["DLL3",] >= threshold(spe, "DLL3", 0.5) & assay(neoplastic, "exprs")["BCAN",] >= threshold(spe, "BCAN", 0.5) 
table(anno$NPC)

anno$OPC = assay(neoplastic, "exprs")["OLIG1",] >= threshold(spe, "OLIG1", 0.5) & assay(neoplastic, "exprs")["SCD5",] >= threshold(spe, "SCD5", 0.5) 
table(anno$OPC)

neoplastic$neoplastic_anno <- "undefined"

anno = lapply(anno, which)

neoplastic$neoplastic_anno[setdiff(setdiff(setdiff(setdiff(anno$MES, anno$NPC), anno$OPC), anno$AC), anno$MES_hypoxia)] <- "MES"

neoplastic$neoplastic_anno[setdiff(setdiff(setdiff(setdiff(anno$NPC, anno$MES), anno$OPC), anno$AC), anno$MES_hypoxia)] <- "NPC"

neoplastic$neoplastic_anno[setdiff(setdiff(setdiff(setdiff(anno$OPC, anno$MES), anno$NPC), anno$AC), anno$MES_hypoxia)] <- "OPC"

neoplastic$neoplastic_anno[setdiff(setdiff(setdiff(setdiff(anno$AC, anno$MES), anno$NPC), anno$OPC), anno$MES_hypoxia)] <- "AC"

neoplastic$neoplastic_anno[setdiff(setdiff(setdiff(setdiff(anno$MES_hypoxia, anno$NPC), anno$OPC), anno$AC), anno$MES)] <- "MES Hypoxia"

table(neoplastic$neoplastic_anno)


all_cells <- as.data.frame(colData(neoplastic)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname, neoplastic_anno) %>%
    left_join(all_cells, ., by = "rowname")

all_cells$neoplastic_anno = ifelse(is.na(all_cells$neoplastic_anno), 
                               "undefined", all_cells$neoplastic_anno)

all_cells$neoplastic_anno = ifelse(all_cells$neoplastic_anno == "MES Hypoxia", 
                                   "MES", all_cells$neoplastic_anno)

all(rownames(colData(spe)) == all_cells$rowname)

unique(all_cells$immune_anno[all_cells$neoplastic_anno != "undefined"])

all_cells$immune_anno[all_cells$neoplastic_anno != "undefined"] = all_cells$neoplastic_anno[all_cells$neoplastic_anno != "undefined"]

table(all_cells$immune_anno)

spe$cell_anno = all_cells$immune_anno 

table(spe$cell_anno)


# LABELLING THE NORMAL BRAIN CELLS ---------------------------------------------

rm(neoplastic)

walk(spe@metadata$markers$stroma, ~{
    
    plot_marker_exprs(spe_object = spe, 
                      marker = .x, 
                      quantile_cuttoff = 0.5,  
                      outdir = io$output$exprs_histograms)
})

stroma <- spe[,colData(spe)$cell_anno == "undefined"]
dim(stroma)[2]

anno <- list()

anno$neuron = assay(stroma, "exprs")["NeuN_FOX3",] >= threshold(spe, "NeuN_FOX3", 0.6) & assay(stroma, "exprs")["CD56",] >= threshold(spe, "CD56", 0.6) 
table(anno$neuron)

anno$astrocyte = assay(stroma, "exprs")["GFAP",] >= threshold(spe, "GFAP", 0.5)
table(anno$astrocyte)

anno$oligo = assay(stroma, "exprs")["MOG",] >= threshold(spe, "MOG", 0.5)
table(anno$oligo)

stroma$stroma_anno <- "undefined"

anno = lapply(anno, which)

stroma$stroma_anno[setdiff(setdiff(anno$neuron, anno$astrocyte), anno$oligo)] <- "Neuron"

stroma$stroma_anno[setdiff(setdiff(anno$astrocyte, anno$neuron), anno$oligo)] <- "Astrocyte"

stroma$stroma_anno[setdiff(setdiff(anno$oligo, anno$astrocyte), anno$neuron)] <- "Oligodendrocyte"

table(stroma$stroma_anno)


all_cells <- as.data.frame(colData(stroma)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname, stroma_anno) %>%
    left_join(all_cells, ., by = "rowname")

all_cells$stroma_anno = ifelse(is.na(all_cells$stroma_anno), 
                                   "undefined", all_cells$stroma_anno)

all(rownames(colData(spe)) == all_cells$rowname)

unique(all_cells$immune_anno[all_cells$stroma_anno != "undefined"])

all_cells$immune_anno[all_cells$stroma_anno != "undefined"] = all_cells$stroma_anno[all_cells$stroma_anno != "undefined"]

table(all_cells$immune_anno)

spe$cell_anno = all_cells$immune_anno 

table(spe$cell_anno)

rm(stroma)


imap(spe@metadata$markers, ~{
    plot_cluster_expression(
        spe, 
        markers = .x,
        cluster_labels = "cell_anno",
        fun = "mean",
        scale_fun = "zscore",
        scale = "row",
        plot_title = .y,
        outfile = file.path(io$output$manual_gating,
                            glue::glue("{.y}_marker_expression.pdf"))
    )
})


# LABELLING THE CELL STATES ----------------------------------------------------

plot_marker_exprs(spe_object = spe, marker = "JARID2_C_Terminus", 
                  quantile_cuttoff = 0.5,  outdir = io$output$exprs_histograms)

spe$high_JARID2_C = assay(spe, "exprs")["JARID2_C_Terminus",] >= threshold(spe, "JARID2_C_Terminus", 0.6) 
table(spe$high_JARID2_C)


plot_marker_exprs(spe_object = spe, marker = "JARID2_N_Terminus", 
                  quantile_cuttoff = 0.75,  outdir = io$output$exprs_histograms)

spe$high_JARID2_N = assay(spe, "exprs")["JARID2_N_Terminus",] >= threshold(spe, "JARID2_N_Terminus", 0.75) 
table(spe$high_JARID2_N)


plot_marker_exprs(spe_object = spe, marker = "Ki67", 
                  quantile_cuttoff = 0.75,  outdir = io$output$exprs_histograms)

spe$high_proliferating = assay(spe, "exprs")["Ki67",] >= threshold(spe, "Ki67", 0.75) 
table(spe$high_proliferating)


plot_marker_exprs(spe_object = spe, marker = "HIF1A", 
                  quantile_cuttoff = 0.60,  outdir = io$output$exprs_histograms)

spe$high_hypoxia = assay(spe, "exprs")["HIF1A",] >= threshold(spe, "HIF1A", 0.6) 
table(spe$high_hypoxia)



plot_marker_exprs(spe_object = spe, marker = "EZH2", 
                  quantile_cuttoff = 0.50,  outdir = io$output$exprs_histograms)

spe$transcriptionally_repressive = assay(spe, "exprs")["EZH2",] >= threshold(spe, "EZH2", 0.5) 
table(spe$transcriptionally_repressive)


plot_marker_exprs(spe_object = spe, marker = "SNAI1", 
                  quantile_cuttoff = 0.50,  outdir = io$output$exprs_histograms)

spe$emt = assay(spe, "exprs")["SNAI1",] >= threshold(spe, "SNAI1", 0.5) 
table(spe$emt)


plot_marker_exprs(spe_object = spe, marker = "SOX2", 
                  quantile_cuttoff = 0.50,  outdir = io$output$exprs_histograms)

spe$proliferating_stem_cell = assay(spe, "exprs")["SOX2",] >= threshold(spe, "SOX2", 0.5) 
table(spe$proliferating_stem_cell)



plot_marker_exprs(spe_object = spe, marker = "TNC", 
                  quantile_cuttoff = 0.50,  outdir = io$output$exprs_histograms)

plot_marker_exprs(spe_object = spe, marker = "TGFBeta", 
                  quantile_cuttoff = 0.50,  outdir = io$output$exprs_histograms)


spe$quiescent_stem_cell = assay(spe, "exprs")["TNC",] >= threshold(spe, "TNC", 0.5) & assay(spe, "exprs")["TGFBeta",] >= threshold(spe, "TGFBeta", 0.5)
table(spe$quiescent_stem_cell)


colData(spe)


# CLUSTER THE UNDEFINED CELL LABELS  -------------------------------------------

undefined = spe[,colData(spe)$cell_anno == "undefined"]
dim(undefined)[2]

# Perform SOM clustering
som_clusters <- clusterRows(reducedDim(undefined, "harmony"), 
                            SomParam(100), 
                            full = TRUE)

# Cluster the 100 SOM codes into larger clusters
ccp <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = t(som_clusters$objects$som$codes[[1]]),
    maxK = 45,
    reps = 100,
    pItem = 1.0,
    pFeature = 1.0, 
    clusterAlg = "km", 
    distance = "euclidean", 
    seed = 1234, 
    plot = NULL)

final_som_clusters <- lapply(2:length(ccp), function(x){
    
    ccp[[x]][["consensusClass"]][som_clusters$clusters]
    
})  

names(final_som_clusters) <- glue::glue("cluster_{1:length(final_som_clusters)+1}")

undefined$som_clusters = factor(as.numeric(final_som_clusters[[length(final_som_clusters)]]))

undefined@metadata$color_vectors$som_clusters <- setNames(
    
    viridis::turbo(length(levels(undefined$som_clusters))), 
    levels(undefined$som_clusters)
    )



pdf(file = file.path(io$output$manual_gating, "undefined_label_clusters.pdf"), 
    onefile = T, width = 15, height = 10)

dittoDimPlot(undefined, 
             var = "som_clusters", 
             reduction.use = "pacMap_harmony", 
             size = 0.2,
             color.panel = undefined@metadata$color_vectors$som_clusters,
             do.label = TRUE) +
    ggtitle("SOM Clusters (Harmony Embeddings)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", hjust = 0.5)
    )

dev.off()


pdf(file = file.path(io$output$manual_gating, "immune_clustering_counts.pdf"), 
    onefile = T, width = 9, height = 7)

cluster_cell_counts(undefined, cluster_label = "som_clusters", 
                    cluster_colors = undefined@metadata$color_vectors$som_clusters,
                    scale_cols = F)

dev.off()


expr = t(undefined@assays@data$exprs)[,rownames(undefined)[rowData(undefined)$downstream == 1]]
expr <- cbind(expr, cluster = as.numeric(undefined@colData[["som_clusters"]]))

expr[1:5,]

MEM_values <- MEM(exp_data = expr,
                  transform = FALSE, 
                  choose.markers = FALSE, # Change choose.markers to TRUE to see and select channels in console
                  markers = "all",
                  choose.ref = FALSE,
                  zero.ref = FALSE,
                  rename.markers = FALSE, # Change rename.markers to TRUE to see and choose new names for channels in console
                  new.marker.names = "none",
                  IQR.thresh = NULL)

build_heatmaps(MEM_values, 
               cluster.MEM = "none", 
               cluster.medians = "none",
               cluster.IQRs = "none", 
               display.thresh = 5, 
               output.files = TRUE,
               labels = FALSE, 
               only.MEMheatmap = FALSE)




# After running the MEM scoring on the arcsinh transformed counts along with The
# SOM cluster labels for the undefined labels, the following clusters were assigned:

undefined$undef_labs = "undefined"

# AC
# 44 : ▲ SLC1A3_EAAT1+5
undefined$undef_labs[undefined$som_clusters == "44"] = "AC"


# AC_MES
# 43 : ▲ SLC1A3_EAAT1+7 ANEXIN_A2+7 ANXA1+7 SOD2+6 CHI3L1+5 JARID2_N_Terminus+5
# 32 : ▲ ANXA1+6 SLC1A3_EAAT1+5 ANEXIN_A2+5
# 36 : ▲ ANXA1+6 SLC1A3_EAAT1+5 ANEXIN_A2+5
# 39 : ▲ SLC1A3_EAAT1+6 ANXA1+6 TGFBeta+5
# 26 : ▲ SLC1A3_EAAT1+6 ANXA1+5
undefined$undef_labs[undefined$som_clusters %in% c("43","32","36","39","26")] = "AC_MES"


# Endothelial
# 5 : ▲ SMA+6
undefined$undef_labs[undefined$som_clusters %in% c("5")] = "Endothelial" 


# MES
# 4 : ▲ ANEXIN_A2+5
# 34 : ▲ ANXA1+7 SOD2+5 ANEXIN_A2+5
# 37 : ▲ ANXA1+6 SOD2+5 ANEXIN_A2+5
undefined$undef_labs[undefined$som_clusters %in% c("4","34","37")] = "MES" 


# Macrophage
# 13 : ▲ IBA1+5
# 1 : ▲ IBA1+6
# 2 : ▲ IBA1+5
undefined$undef_labs[undefined$som_clusters %in% c("1","2","13")] = "Macrophage" 

# Astrocyte
# 31 : ▲ GFAP+5
undefined$undef_labs[undefined$som_clusters %in% c("31")] = "Astrocyte" 

table(undefined$undef_labs)


all_cells <- as.data.frame(colData(undefined)) %>% 
    tibble::rownames_to_column() %>%
    select(rowname, undef_labs) %>%
    left_join(all_cells, ., by = "rowname")

all_cells$undef_labs = ifelse(is.na(all_cells$undef_labs), 
                               "undefined", all_cells$undef_labs)

all(rownames(colData(spe)) == all_cells$rowname)

unique(all_cells$immune_anno[all_cells$undef_labs != "undefined"])

all_cells$immune_anno[all_cells$undef_labs != "undefined"] = all_cells$undef_labs[all_cells$undef_labs != "undefined"]

table(all_cells$immune_anno)

spe$cell_anno = all_cells$immune_anno 

table(spe$cell_anno)

# Percent annotation
sum(spe$cell_anno != "undefined")/ncol(spe) * 100

# SAVE OUTPUT  -----------------------------------------------------------------

saveRDS(spe, file.path(io$inputs$comp_data, "spe_comp.rds"))

# END ----






