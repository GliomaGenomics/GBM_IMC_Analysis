# Author: Shoaib Ajaib
# Date: 02/09/2024

# USAGE ------------------------------------------------------------------------
# This script is designed to perform spatial analysis on the data generated in
# the previous steps. The following types of analysis will be performed:
#
# 1. Create spatial interaction graphs
# 2. Spatial communities
# 3. Cellular neighbourhoods/niches
# 4. Cellular interactions

# OPTIONS ----------------------------------------------------------------------
# options(scipen = 999)

# PACKAGES ---------------------------------------------------------------------
library(ggplot2)
library(viridis)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(SpatialExperiment)
library(IMCfuncs)
library(imcRtools)

# I/O --------------------------------------------------------------------------
io <- list(
  inputs = list(
    input_dir = "outputs/cell_phenotyping",
    prevelance_out = "outputs/cell_prevalences"
  ),
  outputs = list(
    out_dir = "outputs/spatial_analysis"
  ),
  plots = list()
)

# create out directory
ndirs(io$outputs)

# create time-stamped output directory
io$outputs$temp_out <- nd(path = io$outputs$out_dir)

# obtain the most-recent data from a time-stamped directory
find_file <- function(dir_path,
                      file_pattern,
                      dir_pattern = "^[T0-9-]+$",
                      dir_select = c("recent", "oldest")) {
  if (missing(dir_path)) stop("No directory path provided")
  if (missing(file_pattern)) stop("No file pattern provided")

  dirs_found <- list.files(dir_path, pattern = dir_pattern)

  if (length(dirs_found) == 0) stop("No directories found")

  dir_selected <- switch(match.arg(dir_select, several.ok = FALSE),
    recent = dirs_found[order(dirs_found, decreasing = TRUE)][[1]],
    oldest = dirs_found[order(dirs_found, decreasing = FALSE)][[1]]
  )

  file_found <- list.files(
    path = file.path(dir_path, dir_selected),
    pattern = file_pattern,
    ignore.case = TRUE,
    recursive = FALSE,
    full.names = TRUE
  )

  if (length(file_found) == 0) {
    stop("No files found")
  } else if (length(file_found) > 1) {
    stop("Multiple files found")
  } else {
    return(file_found)
  }
}

io$inputs$data <- find_file(io$inputs$input_dir, file_pattern = "spe_downstream")
io$inputs$prevelance_out <- find_file(io$inputs$prevelance_out, file_pattern = "prevelance_data")

rm(find_file)

# LOAD DATA --------------------------------------------------------------------
spe <- readRDS(io$inputs$data)
regions <- readRDS(io$inputs$prevelance_out)

# LABELLED CELL COUNTS ---------------------------------------------------------
pdf(
  file = nf("labelled_cell_counts.pdf", io$outputs$temp_out),
  width = 10,
  height = 10,
  onefile = TRUE
)

lab_spe@colData[c("patient", "manual_gating")] %>%
  as.data.frame() %>%
  group_by(patient) %>%
  summarise(cells = n()) %>%
  arrange(cells) %>%
  ungroup() %>%
  mutate(name = factor(patient, levels = patient)) %>%
  ggplot(aes(x = cells, y = name)) +
  xlab("ncells") +
  ylab("") +
  geom_col(fill = "slateblue") +
  theme_classic()

lab_spe@colData[c("patient", "surgery", "manual_gating")] %>%
  as.data.frame() %>%
  mutate(patient_surgery = paste(patient, surgery, sep = "_")) %>%
  group_by(patient_surgery) %>%
  summarise(cells = n()) %>%
  arrange(cells) %>%
  ungroup() %>%
  mutate(name = factor(patient_surgery, levels = patient_surgery)) %>%
  ggplot(aes(x = cells, y = name)) +
  xlab("ncells") +
  ylab("") +
  geom_col(fill = "slateblue") +
  theme_classic()

lab_spe@colData[c("sample_id", "manual_gating")] %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarise(cells = n()) %>%
  arrange(cells) %>%
  ungroup() %>%
  mutate(name = factor(sample_id, levels = sample_id)) %>%
  ggplot(aes(x = cells, y = name)) +
  xlab("ncells") +
  ylab("") +
  geom_col(fill = "slateblue") +
  theme_classic()

dev.off()

# LABELLED CELL SPATIAL COORDINATES --------------------------------------------
plot_cells <- function(spe_obj,
                       patient = c("64", "67", "71", "82", "84"),
                       anno = c("fine", "main"),
                       sample_regions = c("001", "002", "003")) {
  anno_info <- switch(
    EXPR = match.arg(anno, several.ok = FALSE),
    fine = list(
      label = "manual_gating",
      colors = spe_obj@metadata$v2_colours$cells
    ),
    main = list(
      label = "main_anno_v2",
      colors = spe_obj@metadata$v2_colours$cell_groups
    )
  )


  lab_spe <- spe[, spe$ROI %in% sample_regions & !is.na(spe[[anno_info$label]])]

  lab_spe <- as_tibble(spatialCoords(lab_spe)) %>%
    dplyr::rename(x = Pos_X, y = Pos_Y) %>%
    cbind(
      id = lab_spe$sample_id,
      patient_id = lab_spe$patient,
      surgery = lab_spe$surgery,
      label = lab_spe[[anno_info$label]]
    ) %>%
    dplyr::filter(patient_id %in% patient)

  out <- list()

  out$surgery <- lab_spe %>%
    ggplot(aes(x = x, y = y, color = label)) +
    geom_point(size = 2, alpha = 0.75) +
    facet_grid(surgery ~ patient_id) +
    scale_color_manual(values = anno_info$colors) +
    IMCfuncs::facetted_comp_bxp_theme() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 10)))

  out$sample <- lab_spe %>%
    ggplot(aes(x = x, y = y, color = label)) +
    geom_point(size = 2, alpha = 0.75) +
    facet_wrap(~id) +
    scale_color_manual(values = anno_info$colors) +
    IMCfuncs::facetted_comp_bxp_theme() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 10)))


  return(out)
}


cell_plots <- purrr::map(unique(spe$patient), ~ {
  out <- list(
    surgery = list(),
    sample = list()
  )

  out$surgery$main <- plot_cells(spe_obj = spe, patient = .x, anno = "main")$surgery
  out$surgery$fine <- plot_cells(spe_obj = spe, patient = .x, anno = "fine")$surgery

  out$sample$main <- plot_cells(spe_obj = spe, patient = .x, anno = "main")$sample
  out$sample$fine <- plot_cells(spe_obj = spe, patient = .x, anno = "fine")$sample

  return(out)
})
names(cell_plots) <- paste("patient", unique(spe$patient), sep = "_")

tosave <- map(cell_plots, ~ pluck(.x, "surgery")) %>% list_flatten()

pdf(
  file = nf("surgey_labelled_cell_spatial_coords.pdf", io$outputs$temp_out),
  width = 12,
  height = 14,
  onefile = TRUE
)
print(tosave)
dev.off()

tosave <- map(cell_plots, ~ pluck(.x, "sample")) %>% list_flatten()

pdf(
  file = nf("sample_labelled_cell_spatial_coords.pdf", io$outputs$temp_out),
  width = 18,
  height = 12,
  onefile = TRUE
)
print(tosave)
dev.off()

rm(cell_plots, plot_cells, tosave)

# KNN SPATIAL INTERACTION GRAPHS -----------------------------------------------
lab_spe <- spe[, spe$ROI %in% c("001", "002", "003") & !is.na(spe$manual_gating)]

for (i in seq(5, 30, 5)) {
  lab_spe <- buildSpatialGraph(
    object = lab_spe,
    img_id = "sample_id",
    type = "knn",
    k = i,
    name = paste0("k_", i, collapse = "")
  )
}

plot_interaction_graphs <- function(spe_obj,
                                    graph_name,
                                    patients = c("64", "67", "71", "82", "84"),
                                    node_label = "manual_gating",
                                    node_colours = spe@metadata$v2_colours$cells) {
  out <- vector(mode = "list", length = length(patients))
  names(out) <- paste("patient", patients, sep = "_")

  for (i in seq_along(patients)) {
    out[[i]] <- imcRtools::plotSpatial(
      object = spe_obj[, spe_obj$patient == patients[[i]] & spe_obj$surgery %in% c("Prim", "Rec")],
      node_color_by = "manual_gating",
      img_id = "sample_id",
      colPairName = graph_name,
      nodes_first = FALSE,
      ncols = 3,
      draw_edges = TRUE,
      edge_color_fix = "grey"
    ) +
      ggtitle(glue::glue("{graph_name} interaction graph")) +
      scale_color_manual(values = node_colours) +
      IMCfuncs::facetted_comp_bxp_theme() +
      theme(legend.position = "right") +
      guides(
        color = guide_legend(
          override.aes = list(size = 10),
          ncol = 1,
          bycol = TRUE
        )
      )
  }

  return(out)
}

outplots <- vector("list", length = length(colPairNames(lab_spe)))
names(outplots) <- colPairNames(lab_spe)

for (i in seq_along(colPairNames(lab_spe))) {
  outplots[[i]] <- plot_interaction_graphs(
    spe_obj = lab_spe,
    graph_name = colPairNames(lab_spe)[[i]],
    patients = c("64", "82") # only plotting the most and least dense patient samples
  )
}

outplots <- list_flatten(outplots)

pdf(
  file = nf("knn_sweep_graphs.pdf", io$outputs$temp_out),
  width = 20,
  height = 15,
  onefile = TRUE
)
print(outplots)
dev.off()

# Each image comprises on a varying number of cells and so a suitable K value
# need to be established to identify large and small cluster. After visually
# inspecting a number of different K values using the most dense (82) and
# least dense (64) samples, k = 15 seems to be a adequate number of
# neighbours to use:
lab_spe <- spe[, spe$ROI %in% c("001", "002", "003") & !is.na(spe$manual_gating)]

lab_spe <- buildSpatialGraph(
  object = lab_spe,
  img_id = "sample_id",
  type = "knn",
  k = 15,
  name = "k_15"
)

outplot <- plot_interaction_graphs(spe_obj = lab_spe, graph_name = "k_15")

pdf(
  file = nf("k_15_graphs.pdf", io$outputs$temp_out),
  width = 20,
  height = 15,
  onefile = TRUE
)
print(outplot)
dev.off()

rm(outplots, i)
# DELAUNAY SPATIAL INTERACTION GRAPHS ------------------------------------------
# The Delaunay triangulation connects points in such a way that no point is inside
# the circumcircle of any triangle. Neighbours are defined as points that share
# an edge in the triangulation. This make is a good method for capturing the
# inherent spatial structure of the especially when the points are irregularly spaced.
#
# Pros:
#     Captures the underlying spatial structure well.
#     No arbitrary distance thresholds or fixed neighbour counts.
# Cons:
#     Can be computationally intensive for large datasets.
#     May result in some connections that are not meaningful in a biological context.

# Although the method does not inherently involve any distance thresholds we can
# prune the graph to only include edges that are within a certain distance of each
# other. This can be useful for removing spurious connections that may be present
# and limit the number of connections to only those that are biologically relevant.

lab_spe <- buildSpatialGraph(
  object = lab_spe,
  img_id = "sample_id",
  type = "delaunay",
  max_dist = 50,
  name = "delaunay_50"
)

outplot <- plot_interaction_graphs(
  spe_obj = lab_spe,
  graph_name = "delaunay_50"
)

pdf(
  file = nf("delaunay_50_graphs.pdf", io$outputs$temp_out),
  width = 20,
  height = 15,
  onefile = TRUE
)
print(outplot)
dev.off()

rm(outplot, plot_interaction_graphs)

# SPATIAL COMMUNITY ANALYSIS ---------------------------------------------------
# This method was first described in:
# Jackson, H. W. et al. The single-cell pathology landscape of breast cancer. Nature 578, 615–620 (2020).
#
# Cells are clustered solely based on their interactions as defined by the
# spatial object graph. Communities with less than a set number of interactions
# are excluded.

set.seed(123)
lab_spe <- detectCommunity(
  object = lab_spe,
  colPairName = "delaunay_50",
  name = "delaunay_ca",
  size_threshold = 10
)

set.seed(123)
lab_spe <- detectCommunity(
  object = lab_spe,
  colPairName = "k_15",
  name = "knn_ca",
  size_threshold = 10
)

pdf(
  file = nf("spatial_communities.pdf", io$outputs$temp_out),
  width = 10,
  height = 35,
  onefile = TRUE
)

plotSpatial(lab_spe,
  node_color_by = "delaunay_ca",
  img_id = "sample_id",
  node_size_fix = 0.75,
  ncols = 3
) +
  ggtitle("Spatial Communities - Delaunay (50)") +
  scale_color_manual(values = rev(colors())) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(legend.position = "none")


plotSpatial(lab_spe,
  node_color_by = "knn_ca",
  img_id = "sample_id",
  node_size_fix = 0.75,
  ncols = 3
) +
  ggtitle("Spatial Communities - KNN (15)") +
  scale_color_manual(values = rev(colors())) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(legend.position = "none")

dev.off()


plot_ca_exprs <- function(spe_obj,
                          surgery,
                          anno,
                          community,
                          graph_name) {
  filt_spe <- lab_spe[, lab_spe$surgery == surgery]

  filt_spe <- prop.table(
    x = table(filt_spe[[community]], filt_spe[[anno]]),
    margin = 1
  )

  graph_lab <- stringr::str_split_i(community, "_", 1)
  graph_param_lab <- stringr::str_split_i(graph_name, "_", 2)

  out_heatmap <- pheatmap::pheatmap(
    mat = filt_spe,
    color = colorRampPalette(c("dark blue", "white", "dark red"))(15),
    show_rownames = FALSE,
    cluster_cols = FALSE,
    scale = "column",
    main = glue::glue("{surgery} communities - {graph_lab} ({graph_param_lab})"),
    silent = TRUE
  )

  return(out_heatmap)
}


plot_params <- expand.grid(
  surgery = c("Prim", "Rec"),
  anno = c("main_anno_v2"),
  community = c("knn_ca", "delaunay_ca"),
  stringsAsFactors = FALSE
) %>%
  mutate(graph_name = case_when(
    community == "knn_ca" ~ "k_15",
    community == "delaunay_ca" ~ "delaunay_50"
  ))


tosave <- pmap(plot_params, plot_ca_exprs, spe_obj = lab_spe)

pdf(
  file = nf("spatial_community_main_exprs.pdf", io$outputs$temp_out),
  width = 5,
  height = 10,
  onefile = TRUE
)
for (p in tosave) {
  grid::grid.newpage()
  print(p)
}
dev.off()

plot_params <- expand.grid(
    surgery = c("Prim", "Rec"),
    anno = c("manual_gating"),
    community = c("knn_ca", "delaunay_ca"),
    stringsAsFactors = FALSE
) %>%
    mutate(graph_name = case_when(
        community == "knn_ca" ~ "k_15",
        community == "delaunay_ca" ~ "delaunay_50"
    ))

tosave <- pmap(plot_params, plot_ca_exprs, spe_obj = lab_spe)

pdf(
    file = nf("spatial_community_fine_exprs.pdf", io$outputs$temp_out),
    width = 10,
    height = 10,
    onefile = TRUE
)
for (p in tosave) {
    grid::grid.newpage()
    print(p)
}
dev.off()

rm(p, plot_params, tosave, plot_ca_exprs)

# CELLULAR NEIGHBORHOOD ANALYSIS -----------------------------------------------
# Rather than clustering solely based on the interaction graph this method first 
# aggregates cells based on information contained in their direct neighbourhood and
# then performs clustering to define cellular neighbourhoods.
# 
# This method we will employ here has previously been used in:
# 
# Goltsev et al. 2018. Cell 174: 968–81.
# Schürch et al. 2020. Cell 182: 1341–59.
# 
# Aggregation is performed in 2 different ways:
# 
# - For each cell, compute the fraction of cells of a certain type among its neighbours.
# - For each cell it aggregates (mean) the expression counts across all neighbouring cells.


# SAVE DATA --------------------------------------------------------------------
# END --------------------------------------------------------------------------
