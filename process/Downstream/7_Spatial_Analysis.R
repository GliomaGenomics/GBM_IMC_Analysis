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
library(patchwork)
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

# Considering the first three regions across each sample and the labelled cells
# lab_spe <- spe[, spe$ROI %in% c("001", "002", "003") & !is.na(spe$manual_gating)]

# PREVIOUSLY CREATED DATA ------------------------------------------------------
lab_spe <- list.files(
  path = "outputs/spatial_analysis/2024-09-05T12-01-52",
  pattern = "lab_spe",
  full.names = TRUE
) %>%
  readRDS()

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

# SCORE CELL STATES ------------------------------------------------------------
state_markers <- list(
  proliferative = c(
    spe@metadata$markers$cell_states$Proliferating,
    spe@metadata$markers$cell_states$Proliferating_stem_cell
  ),
  hypoxia = spe@metadata$markers$cell_states$Hypoxia,
  queiescence = spe@metadata$markers$cell_states$Quiescent_stem_cell,
  EMT = spe@metadata$markers$cell_states$Epithelial_mesenchymal_transition
)

classify_cell_states <- function(marker_list,
                                 expression_matrix,
                                 high_threshold = 2, low_threshold = -2) {
  # Ensure expression_matrix has row names
  if (is.null(rownames(expression_matrix))) {
    stop("The expression_matrix must have row names corresponding to marker genes.")
  }

  # Ensure marker_list is a named list
  if (is.null(names(marker_list)) || any(names(marker_list) == "")) {
    stop("marker_list must be a named list with cell state names.")
  }


  cell_ids <- colnames(expression_matrix)
  if (is.null(cell_ids)) {
    stop("The expression_matrix must have column names corresponding to cell IDs.")
  }

  state_scores <- data.frame(cell = cell_ids, stringsAsFactors = FALSE)

  for (state in names(marker_list)) {
    markers <- marker_list[[state]]

    missing_markers <- setdiff(markers, rownames(expression_matrix))
    if (length(missing_markers) > 0) {
      warning(paste(
        "Markers not found for state", state, ":",
        paste(missing_markers, collapse = ", ")
      ))
      markers <- setdiff(markers, missing_markers)

      if (length(markers) == 0) {
        warning(paste(
          "No valid markers left for state", state,
          ". Skipping this state."
        ))
        next
      }
    }

    # Extract the expression data for the markers
    marker_expr <- expression_matrix[markers, , drop = FALSE]

    # Sum the z-scores if multiple markers, else take the z-score directly
    if (length(markers) > 1) {
      # Sum the z-scores across markers for each cell
      state_sum <- colSums(marker_expr, na.rm = TRUE)
    } else {
      # Single marker: take the z-score directly
      state_sum <- as.numeric(marker_expr)
      names(state_sum) <- colnames(marker_expr)
    }

    # Classify cells as High or Low based on thresholds
    high_class <- state_sum > high_threshold
    low_class <- state_sum < low_threshold

    # Add classification to the state_scores dataframe
    state_scores[[paste0(state, "_high")]] <- high_class
    state_scores[[paste0(state, "_low")]] <- low_class
  }

  return(state_scores)
}

summarize_states <- function(state_scores, marker_list) {
  for (state in names(marker_list)) {
    high_col <- paste0(state, "_high")
    low_col <- paste0(state, "_low")

    # Check if the columns exist (in case some states were skipped due to missing markers)
    if (!(high_col %in% colnames(state_scores)) || !(low_col %in% colnames(state_scores))) {
      next
    }

    high_count <- sum(state_scores[[high_col]], na.rm = TRUE)
    low_count <- sum(state_scores[[low_col]], na.rm = TRUE)

    cat(sprintf("State: %s\n", state))
    cat(sprintf("  High: %d cells\n", high_count))
    cat(sprintf("  Low: %d cells\n\n", low_count))
  }
}


plot_state_distribution <- function(expression_matrix = t(assay(spe, "zscore")),
                                    state_name,
                                    state_markers,
                                    high_threshold = 1.2,
                                    low_threshold = -1.2) {
  histogram_df <- expression_matrix[, state_markers, drop = FALSE] %>%
    as.data.frame()

  if (ncol(histogram_df) > 1) {
    histogram_df <- histogram_df %>%
      mutate(state_zscore = rowSums(across(everything()))) %>%
      select(state_zscore)
  } else {
    histogram_df <- histogram_df %>%
      dplyr::rename(state_zscore = 1)
  }


  histogram_df %>%
    ggplot(aes(x = state_zscore)) +
    geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
    geom_vline(xintercept = c(-1.2), linetype = "dashed", linewidth = 1, color = "darkred") +
    geom_vline(xintercept = c(1.2), linetype = "dashed", linewidth = 1, color = "darkblue") +
    labs(
      title = glue::glue("{tools::toTitleCase(state_name)} Score Distribution"),
      subtitle = glue::glue("low threshold (-1.2); high threshold (1.2)"),
      x = "z-score",
      y = "Number of cells"
    ) +
    scale_x_continuous(limits = c(-5, 5)) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(face = "italic")
    )
}

state_scores <- classify_cell_states(
  marker_list = state_markers,
  expression_matrix = assay(spe, "zscore"),
  high_threshold = 1.2,
  low_threshold = -1.2
)

state_histograms <- purrr::imap(
  state_markers,
  ~ plot_state_distribution(state_name = .y, state_markers = .x)
)

pdf(
  file = nf("cell_state_score_histograms.pdf", io$outputs$temp_out),
  width = 10,
  height = 10,
  onefile = TRUE
)
print(state_histograms)
dev.off()

summarize_states(state_scores = state_scores, marker_list = state_markers)

saveRDS(state_scores, file = nf("cell_state_scores.rds", io$outputs$temp_out))

rm(
  state_histograms, state_markers,
  classify_cell_states, plot_state_distribution, summarize_states
)

# MALIGNANT CELL STATE SCORES PROPORTIONS --------------------------------------
state_scores <- readRDS("outputs/spatial_analysis/2024-09-25T15-19-12/cell_state_scores_2024-09-25T15-26-51.rds")

all_cancer <- colData(lab_spe)[, c("main_anno_v2", "manual_gating", "patient", "surgery")] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(state_scores, by = "cell") %>%
  filter(main_anno_v2 == "Cancer") %>%
  select(-main_anno_v2)

all_cancer_long <- all_cancer %>%
  mutate(proliferative = case_when(
    proliferative_high == TRUE ~ "high",
    proliferative_low == TRUE ~ "low",
    TRUE ~ "medium"
  )) %>%
  mutate(hypoxia = case_when(
    hypoxia_high == TRUE ~ "high",
    hypoxia_low == TRUE ~ "low",
    TRUE ~ "medium"
  )) %>%
  mutate(quiescence = case_when(
    queiescence_high == TRUE ~ "high",
    queiescence_low == TRUE ~ "low",
    TRUE ~ "medium"
  )) %>%
  mutate(EMT = case_when(
    EMT_high == TRUE ~ "high",
    EMT_low == TRUE ~ "low",
    TRUE ~ "medium"
  )) %>%
  select(-ends_with(c("_high", "_low"))) %>%
  tidyr::pivot_longer(
    cols = ends_with(c("proliferative", "hypoxia", "quiescence", "EMT")),
    names_to = "state",
    values_to = "category"
  ) %>%
  mutate(across(manual_gating, as.character))

all_cancer_long_props <- all_cancer_long %>%
  group_by(manual_gating, surgery, state) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  left_join(all_cancer_long, ., by = c("manual_gating", "surgery", "state")) %>%
  group_by(manual_gating, surgery, state, category) %>%
  summarise(
    cells = n(),
    percentage = n() / total_cells[1],
    .groups = "drop"
  ) %>%
  mutate(across(category, ~ factor(., levels = c("low", "medium", "high"))))


svglite::svglite(
  file = nf("cancer_cell_state_props.svg", "outputs/cell_prevalences/2024-08-23T10-01-10"),
  width = 18,
  height = 15
)

all_cancer_long_props %>%
  mutate(across(state, ~ ifelse(. != "EMT", tools::toTitleCase(.), .))) %>%
  ggplot(
    aes(
      x = forcats::fct_rev(surgery),
      y = percentage,
      fill = forcats::fct_rev(category)
    )
  ) +
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  coord_flip() +
  facet_grid(manual_gating ~ state) +
  labs(y = "") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "low" = "#440154FF",
      "medium" = "white",
      "high" = "#FDE725FF"
    )
  ) +
  IMCfuncs::facetted_cell_prop_theme(text_size = 24) +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, ),
    axis.text.x = element_text(angle = 0)
  ) +
  guides(fill = guide_legend(reverse = TRUE))

dev.off()


cancer_state_stats <- all_cancer_long_props %>%
  filter(!category %in% c("medium")) %>%
  select(-percentage) %>%
  tidyr::unite(col = state, state, category, sep = "_") %>%
  tidyr::pivot_wider(
    names_from = surgery,
    values_from = cells,
    values_fill = list(cells = 0)
  ) %>%
  split(.$manual_gating) %>%
  lapply(\(x){
    y <- x[, c("state", "Prim", "Rec")]
    y <- as.data.frame(y)
    rownames(y) <- x$state
    return(y[, c("Prim", "Rec")])
  })

cancer_state_stats <- purrr::imap(cancer_state_stats, ~ {
  out <- rstatix::pairwise_fisher_test(as.matrix(.x), p.adjust.method = "fdr")
  out$cancer_cell_type <- .y
  out <- dplyr::relocate(out, cancer_cell_type)
  return(out)
})
cancer_state_stats <- bind_rows(cancer_state_stats)

openxlsx::write.xlsx(
  x = cancer_state_stats,
  file = nf(
    "cancer_cell_state_fisher_test_stats.xlsx",
    "outputs/cell_prevalences/2024-08-23T10-01-10"
  ),
  asTable = TRUE
)

rm(all_cancer, all_cancer_long, all_cancer_long_props, cancer_state_stats)

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
for (i in seq(5, 30, 5)) {
  lab_spe <- buildSpatialGraph(
    object = lab_spe,
    img_id = "sample_id",
    type = "knn",
    k = i,
    name = paste0("k_", i, collapse = "")
  )
}

# Set the geom point to match the pervious plots
update_geom_defaults("point", list(alpha = 0.75))

plot_interaction_graphs <- function(spe_obj,
                                    graph_name,
                                    patients = c("64", "67", "71", "82", "84"),
                                    node_label = "manual_gating",
                                    node_colours = spe_obj@metadata$v2_colours$cells) {
  out <- vector(mode = "list", length = length(patients))
  names(out) <- paste("patient", patients, sep = "_")

  for (i in seq_along(patients)) {
    out[[i]] <- imcRtools::plotSpatial(
      object = spe_obj[, spe_obj$patient == patients[[i]] & spe_obj$surgery %in% c("Prim", "Rec")],
      node_color_by = "manual_gating",
      img_id = "sample_id",
      colPairName = graph_name,
      nodes_first = FALSE,
      node_size_fix = 2,
      ncols = 3,
      draw_edges = TRUE,
      flip_y = FALSE,
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

rm(outplots, i)

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
io$outputs$temp_ca <- nd(
  directory_name = "spatial_communities",
  path = io$outputs$temp_out,
  add_timestamp = FALSE
)

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
  file = nf("spatial_communities.pdf", io$outputs$temp_ca),
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
  file = nf("spatial_community_main_exprs.pdf", io$outputs$temp_ca),
  width = 5,
  height = 10,
  onefile = TRUE
)
for (p in tosave) {
  grid::grid.newpage()
  print(p)
}
dev.off()

plot_params$anno <- "manual_gating"
tosave <- pmap(plot_params, plot_ca_exprs, spe_obj = lab_spe)

pdf(
  file = nf("spatial_community_fine_exprs.pdf", io$outputs$temp_ca),
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
# Rather than clustering solely based on the interaction graph, we can also
# aggregates cells based on information contained in their direct neighbourhood and
# then perform clustering to define cellular neighbourhoods.
#
# This method has previously been used in the following studies:
#
# Goltsev et al. 2018. Cell 174: 968–81.
# Schürch et al. 2020. Cell 182: 1341–59.
#
# There are two types of cell neighbour information we can use:
#
# 1. The fraction of each cell type in the defined neighbourhood.
# 2. The aggregated (mean, median, etc.) expression for each cell type in the defined neighbourhood.

# Neighbourhoods by cell type fractions:
lab_spe <- aggregateNeighbors(
  object = lab_spe,
  colPairName = "delaunay_50",
  aggregate_by = "metadata",
  count_by = "manual_gating",
  name = "delaunay_cn"
)

lab_spe <- aggregateNeighbors(
  object = lab_spe,
  colPairName = "k_15",
  aggregate_by = "metadata",
  count_by = "manual_gating",
  name = "knn_cn"
)

set.seed(1234)
cn_kmeans <- kmeans(lab_spe$delaunay_cn, centers = 12)
lab_spe$delaunay_cn_clusters <- as.factor(cn_kmeans$cluster)

set.seed(1234)
cn_kmeans <- kmeans(lab_spe$knn_cn, centers = 12)
lab_spe$knn_cn_clusters <- as.factor(cn_kmeans$cluster)


# Neighbourhoods by cell type marker expression:
lab_spe <- aggregateNeighbors(
  object = lab_spe,
  colPairName = "delaunay_50",
  aggregate_by = "expression",
  assay_type = "exprs",
  subset_row = rowData(lab_spe)$name %in% unlist(lab_spe@metadata$labels$markers),
  name = "delaunay_cn_exprs"
)

lab_spe$delaunay_cn_exprs <- as.matrix(lab_spe$delaunay_cn_exprs)
lab_spe$delaunay_cn_exprs[which(is.na(lab_spe$delaunay_cn_exprs))] <- 0

lab_spe <- aggregateNeighbors(
  object = lab_spe,
  colPairName = "k_15",
  aggregate_by = "expression",
  assay_type = "exprs",
  subset_row = rowData(lab_spe)$name %in% unlist(lab_spe@metadata$labels$markers),
  name = "knn_cn_exprs"
)

lab_spe$knn_cn_exprs <- as.matrix(lab_spe$knn_cn_exprs)
lab_spe$knn_cn_exprs[which(is.na(lab_spe$knn_cn_exprs))] <- 0

set.seed(1234)
cn_kmeans <- kmeans(lab_spe$delaunay_cn_exprs, centers = 12)
lab_spe$delaunay_cn_exprs_clusters <- as.factor(cn_kmeans$cluster)

set.seed(1234)
cn_kmeans <- kmeans(lab_spe$knn_cn_exprs, centers = 12)
lab_spe$knn_cn_exprs_clusters <- as.factor(cn_kmeans$cluster)

rm(cn_kmeans)

# VISUALISE CELLULAR NEIGHBORHOODS ---------------------------------------------
io$outputs$temp_cn <- nd(
  directory_name = "cellular_neighbourhoods",
  path = io$outputs$temp_out,
  add_timestamp = FALSE
)

lab_spe@metadata$v2_colours$cn_colours <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(lab_spe$delaunay_cn_clusters))),
  levels(lab_spe$delaunay_cn_clusters)
)

plot_cn <- function(spe_obj,
                    node_label = "manual_gating",
                    patients = c("64", "67", "71", "82", "84"),
                    image_id = "sample_id",
                    plot_title = "Cellular Neighborhoods",
                    plot_subtitle = node_label,
                    brewer_pal = "Paired") {
  out <- vector(mode = "list", length = length(patients))
  names(out) <- paste("patient", patients, sep = "_")

  for (i in seq_along(patients)) {
    out[[i]] <- imcRtools::plotSpatial(
      object = spe_obj[, spe_obj$patient == patients[[i]] & spe_obj$surgery %in% c("Prim", "Rec")],
      node_color_by = node_label,
      img_id = image_id,
      node_size_fix = 2,
      flip_y = FALSE,
      colPairName = NULL,
      ncols = 3
    ) +
      ggtitle(
        label = plot_title,
        subtitle = plot_subtitle
      ) +
      scale_color_brewer(palette = brewer_pal) +
      IMCfuncs::facetted_comp_bxp_theme() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0, size = 14, face = "italic"),
      ) +
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

# Delaunay graph neighbourhoods
pdf(
  file = nf("delaunay_cn.pdf", io$outputs$temp_cn),
  width = 20,
  height = 15,
  onefile = TRUE
)
plot_cn(
  spe_obj = lab_spe,
  node_label = "delaunay_cn_clusters",
  plot_subtitle = glue::glue(
    "graph: Delanuay (50)",
    "neighbour measure: cell fraction",
    "cell label: manual_gating",
    .sep = "\n"
  )
)

plot_cn(
  spe_obj = lab_spe,
  node_label = "delaunay_cn_exprs_clusters",
  plot_subtitle = glue::glue(
    "graph: Delanuay (50)",
    "neighbour measure: marker expression",
    "cell label: manual_gating",
    .sep = "\n"
  )
)
dev.off()

# KNN graph neighbourhoods
pdf(
  file = nf("knn_cn.pdf", io$outputs$temp_cn),
  width = 20,
  height = 15,
  onefile = TRUE
)
plot_cn(
  spe_obj = lab_spe,
  node_label = "knn_cn_clusters",
  plot_subtitle = glue::glue(
    "graph: KNN (15)",
    "neighbour measure: cell fraction",
    "cell label: manual_gating",
    .sep = "\n"
  )
)
plot_cn(
  spe_obj = lab_spe,
  node_label = "knn_cn_exprs_clusters",
  plot_subtitle = glue::glue(
    "graph: KNN (15)",
    "neighbour measure: marker expression",
    "cell label: manual_gating",
    .sep = "\n"
  )
)
dev.off()

# VISUALISE INTERACTION GRAPH AND CN SAMPLE PAIRS (PUBLICATION) ----------------

plot_sample <- function(spe_obj,
                        graph_name = "delaunay_50",
                        sample_ids = c("64Prim_001", "71Rec_001", "84Rec_002"),
                        node_label = "manual_gating",
                        node_size = 2,
                        plot_title = "",
                        node_colours = spe_obj@metadata$v2_colours$cells) {
  single_plot <- imcRtools::plotSpatial(
    object = spe_obj[, spe_obj$sample_id %in% sample_ids],
    node_color_by = "manual_gating",
    img_id = "sample_id",
    colPairName = graph_name,
    nodes_first = FALSE,
    node_size_fix = node_size,
    ncols = 3,
    draw_edges = TRUE,
    flip_y = FALSE,
    edge_color_fix = "grey"
  ) +
    ggtitle(plot_title) +
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

  return(single_plot)
}

plot_sample_cn <- function(spe_obj,
                           sample_ids = c("64Prim_001", "71Rec_001", "84Rec_002"),
                           node_label = "delaunay_cn_clusters",
                           node_size = 2,
                           plot_title = "",
                           plot_subtitle = "",
                           brewer_pal = "Paired") {
  imcRtools::plotSpatial(
    object = spe_obj[, spe_obj$sample_id %in% sample_ids],
    node_color_by = node_label,
    img_id = "sample_id",
    node_size_fix = node_size,
    flip_y = FALSE,
    colPairName = NULL,
    ncols = 3
  ) +
    ggtitle(
      label = plot_title,
      subtitle = plot_subtitle
    ) +
    scale_color_brewer(palette = brewer_pal) +
    IMCfuncs::facetted_comp_bxp_theme() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0, size = 14, face = "italic"),
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 10),
        ncol = 1,
        bycol = TRUE
      )
    )
}

x <- plot_sample(spe_obj = lab_spe)
y <- plot_sample_cn(spe_obj = lab_spe)

svglite::svglite(
  file = nf("interaction_graph_sample_pairs.svg", io$outputs$temp_cn),
  width = 20,
  height = 15
)

x/y + plot_layout(nrow = 2, guides = "collect")

dev.off()

# VISUALISE CELLULAR NEIGHBORHOOD CELL COUNTS ----------------------------------
plot_cn_counts <- function(spe_obj,
                           group_col = "surgery",
                           cn_col = "delaunay_cn_clusters",
                           cn_colours = lab_spe@metadata$v2_colours$cn_colours) {
  cns <- colData(spe_obj)[, c(group_col, cn_col)] %>%
    as.data.frame() %>%
    dplyr::rename(cn = !!sym(cn_col))

  cn_counts <- cns %>%
    group_by(cn) %>%
    summarise(total_cells = n(), .groups = "drop")

  cns %>%
    group_by(!!sym(group_col), cn) %>%
    summarise(cells = n(), .groups = "drop") %>%
    left_join(cn_counts, by = "cn") %>%
    mutate(cn_prop = cells / total_cells) %>%
    ggplot(aes(x = forcats::fct_rev(!!sym(group_col)), y = cn_prop, fill = cn)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_y_continuous(
      labels = scales::percent_format(),
      breaks = seq(0, 1, 0.1)
    ) +
    facet_wrap(~cn) +
    scale_fill_manual(values = cn_colours) +
    coord_flip() +
    IMCfuncs::facetted_cell_prop_theme() +
    ylab("Fraction of total cells in each CN") +
    theme(
      legend.position = "none",
      axis.title.y = element_blank()
    )
}

pdf(
  file = nf("delaunay_cn_counts.pdf", io$outputs$temp_cn),
  width = 20,
  height = 15,
  onefile = TRUE
)
plot_cn_counts(spe_obj = lab_spe, group_col = "surgery")
plot_cn_counts(spe_obj = lab_spe, group_col = "patient")
plot_cn_counts(spe_obj = lab_spe, group_col = "patient_surgery")
dev.off()

# VISUALISE CELLULAR NEIGHBORHOOD ENRICHMENT -----------------------------------
plot_cn_enrichment <- function(spe_obj,
                               filter_col = NA,
                               filter_val = NA,
                               cell_label = "manual_gating",
                               cn_label = "delaunay_cn_clusters",
                               scale_on = c("cell", "cn", "none"),
                               clip_min = NULL,
                               fill_pal = RColorBrewer::brewer.pal(9, "Reds"),
                               plot_type = c("bubble", "heatmap"),
                               min_bubble_size = 6,
                               max_bubble_size = 16,
                               show_heatmap_text = TRUE,
                               rev_cells = FALSE,
                               plot_title = "Cellular Neighborhood (CN) Enrichment",
                               plot_caption = ifelse(plot_type == "heatmap" & show_heatmap_text,
                                 "*Tile text denotes ncells",
                                 ""
                               ),
                               x_axis_title = "\n\nCN (total cells per CN)",
                               y_axis_title = "",
                               fill_lab = ifelse(scale_on == "none", "", "zscore")) {
  data_info <- list(
    samples = "samples: All",
    graph = str_split_i(cn_label, "_", 1),
    neighbor_measure = ifelse(grepl("_exprs_", cn_label), "marker expression", "cell fraction"),
    scale_on_label = switch(match.arg(scale_on, several.ok = FALSE, choices = c("cell", "cn", "none")),
      "cell" = "cell types",
      "cn" = "neighborhoods",
      "none" = "none"
    )
  )

  if (!is.na(filter_col) & !is.na(filter_val)) {
    if (length(filter_val) > 1) stop("only one filter_val can be used")
    if (!filter_col %in% names(colData(spe_obj))) stop("filter_col not found in colData(spe_obj)")

    df <- as.data.frame(colData(spe_obj))[, c(cell_label, cn_label, filter_col)]
    if (!filter_val %in% unique(df[, filter_col])) stop("filter_val not found in filter_col")

    df <- df %>%
      dplyr::filter(!!rlang::sym(filter_col) %in% filter_val)

    filt_vals <- tolower(paste0(unique(df[, filter_col]), collapse = ", "))
    data_info$samples <- paste0(filter_col, ": ", filt_vals)
  } else {
    df <- as.data.frame(colData(spe_obj))[, c(cell_label, cn_label)]
  }

  data_info <- glue::glue(
    "{data_info$samples}",
    "graph: {data_info$graph}",
    "neighbour measure: {data_info$neighbor_measure}",
    "exprs scaled across: {data_info$scale_on_label}",
    "minimum zscore: {clip_min}",
    .null = "n/a",
    .sep = "\n"
  )

  lab_df <- table(df[, cn_label], df[, cell_label]) %>%
    as.data.frame() %>%
    dplyr::rename(text_lab = Freq)

  tab <- switch(
    EXPR = match.arg(scale_on, several.ok = FALSE, choices = c("cell", "cn", "none")),
    "cell" = scale(table(df[, cn_label], df[, cell_label])),
    "cn" = t(scale(table(df[, cell_label], df[, cn_label]))),
    "none" = table(df[, cn_label], df[, cell_label])
  )

  tab <- as.data.frame(tab) %>%
    dplyr::left_join(lab_df, by = c("Var1", "Var2")) %>%
    dplyr::rename(
      cn_label = Var1,
      cell_label = Var2,
      value = Freq
    )

  if (rev_cells) tab$cell_label <- forcats::fct_rev(tab$cell_label)

  if (!is.null(clip_min)) {
    tab$value <- ifelse(tab$value < clip_min, NA, tab$value)
    scale_min_val <- clip_min
  } else {
    scale_min_val <- floor(min(tab$value, na.rm = TRUE)) * 2 / 2
  }

  scale_max_val <- ceiling(max(tab$value, na.rm = TRUE) * 2) / 2
  scale_breaks <- seq(scale_min_val, scale_max_val, 0.5)


  tab$fill_color <- fill_pal[as.numeric(cut(tab$value, breaks = length(fill_pal)))]
  tab$text_color <- IMCfuncs::get_text_color(tab$fill_color)

  tab <- tab %>%
    group_by(cn_label) %>%
    summarise(cn_total = sum(text_lab, na.rm = TRUE), .groups = "keep") %>%
    ungroup() %>%
    dplyr::left_join(tab, ., by = "cn_label") %>%
    dplyr::mutate(
      x_axis_lab = glue::glue("{cn_label}\n({cn_total})"),
      text_lab_percent = round(text_lab / cn_total * 100)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("text_"), ~ ifelse(is.na(value), NA, .)
      )
    )

  outplot_theme <- ggplot2::theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16, face = "italic"),
      plot.caption = element_text(size = 12, face = "italic")
    )

  if (match.arg(plot_type, several.ok = FALSE) == "bubble") {
    tab %>%
      ggplot(aes(x = cn_label, y = cell_label)) +
      geom_tile(fill = "#fffdfa", colour = "grey50", linewidth = 0.1, linetype = 2) +
      geom_point(
        aes(size = text_lab_percent, fill = value),
        shape = 21, colour = "black", stroke = 0.75, na.rm = TRUE
      ) +
      ggplot2::scale_size(
        range = c(min_bubble_size, max_bubble_size),
        breaks = seq(0, 100, 20)
      ) +
      ggplot2::scale_fill_gradientn(
        breaks = scale_breaks,
        limits = c(scale_min_val, scale_max_val),
        colours = fill_pal,
        guide = guide_colourbar(
          barwidth = 2,
          barheight = 15,
          frame.colour = "black",
          ticks.colour = "black",
          ticks.linewidth = 0.35,
          frame.linewidth = 0.35
        )
      ) +
      labs(
        title = plot_title,
        subtitle = data_info,
        caption = plot_caption,
        x = x_axis_title,
        y = y_axis_title,
        fill = fill_lab,
        size = "% cells in CN"
      ) +
      scale_x_discrete(labels = unique(tab$x_axis_lab)) +
      outplot_theme +
      guides(
        size = guide_legend(
          order = 1,
          override.aes = list(fill = "black")
        )
      )
  } else {
    tab <- tab %>%
      mutate(
        across(
          text_lab,
          ~ ifelse(!is.na(.),
            glue::glue("{text_lab}", "({text_lab_percent}%)", .na = "", .sep = "\n"),
            .
          )
        )
      )

    tab %>%
      ggplot(aes(x = cn_label, y = cell_label)) +
      geom_tile(
        aes(fill = value),
        colour = "grey50", linewidth = 0.1, linetype = 2,
        na.rm = TRUE
      ) +
      {
        if (show_heatmap_text) {
          geom_text(
            aes(label = text_lab, colour = text_color),
            size = 4,
            fontface = "bold",
            na.rm = TRUE,
            show.legend = FALSE
          )
        }
      } +
      ggplot2::scale_fill_gradientn(
        colours = fill_pal,
        breaks = scale_breaks,
        limits = c(scale_min_val, scale_max_val),
        na.value = "white",
        guide = guide_colourbar(
          barwidth = 2,
          barheight = 15,
          frame.colour = "black",
          ticks.colour = "black",
          ticks.linewidth = 0.35,
          frame.linewidth = 0.35
        )
      ) +
      ggplot2::scale_color_identity() +
      labs(
        title = plot_title,
        subtitle = data_info,
        caption = plot_caption,
        x = x_axis_title,
        y = y_axis_title,
        fill = fill_lab
      ) +
      scale_x_discrete(labels = unique(tab$x_axis_lab)) +
      outplot_theme
  }
}


cn_enrich_params <- expand.grid(
  filter_col = c(NA, "surgery"),
  filter_val = c(NA, "Prim", "Rec"),
  cell_label = c("main_anno_v2", "manual_gating"),
  stringsAsFactors = FALSE
) %>%
  tibble() %>%
  filter(
    !is.na(filter_col) & !is.na(filter_val) |
      is.na(filter_col) & is.na(filter_val)
  )

update_geom_defaults("point", list(alpha = 1))

cn_enrichment <- list()

cn_enrichment$bubble <- purrr::pmap(cn_enrich_params, ~ {
  plot_cn_enrichment(
    spe_obj = lab_spe,
    filter_col = ..1,
    filter_val = ..2,
    cell_label = ..3,
    scale_on = "cn",
    clip_min = 0.5
  )
})

cn_enrichment$heatmap <- purrr::pmap(cn_enrich_params, ~ {
  plot_cn_enrichment(
    spe_obj = lab_spe,
    plot_type = "heatmap",
    filter_col = ..1,
    filter_val = ..2,
    cell_label = ..3,
    scale_on = "cn",
    clip_min = 0.5
  )
})

pdf(
  file = nf("delaunay_cn_enrichment_bubble.pdf", io$outputs$temp_cn),
  width = 15,
  height = 10,
  onefile = TRUE
)
print(cn_enrichment$bubble)
dev.off()


pdf(
  file = nf("delaunay_cn_enrichment_heatmap.pdf", io$outputs$temp_cn),
  width = 15,
  height = 10,
  onefile = TRUE
)
print(cn_enrichment$heatmap)
dev.off()

rm(cn_enrich_params, cn_enrichment)

plot_cn_labels <- function(spe_obj,
                           anno = c("fine", "main"),
                           patient_col = "patient",
                           surgery_col = "surgery",
                           roi_col = "ROI",
                           cn_label = "delaunay_cn_clusters") {
  if (!cn_label %in% names(colData(spe_obj))) stop("cn_label not found in colData(spe_obj)")

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

  data_info <- list(
    graph = str_split_i(cn_label, "_", 1),
    neighbor_measure = ifelse(grepl("_exprs_", cn_label), "marker expression", "cell fraction")
  )

  data_info <- glue::glue(
    "graph: {data_info$graph}",
    "neighbour measure: {data_info$neighbor_measure}",
    .sep = "\n"
  )

  filt_spe <- lab_spe[, !is.na(lab_spe[[anno_info$label]])]

  filt_spe <- as_tibble(spatialCoords(filt_spe)) %>%
    dplyr::rename(x = Pos_X, y = Pos_Y) %>%
    cbind(
      patient = filt_spe[[patient_col]],
      surgery = filt_spe[[surgery_col]],
      roi = filt_spe[[roi_col]],
      cluster = filt_spe[[cn_label]],
      label = filt_spe[[anno_info$label]]
    ) %>%
    mutate(
      across(surgery, ~ ifelse(. == "Prim", "Primary", "Recurrent"))
    )

  out <- filt_spe %>%
    group_by(patient, surgery) %>%
    group_split(.keep = TRUE)

  names(out) <- sapply(out, \(x) paste(unique(x$patient), unique(x$surgery), sep = "_"))

  out <- purrr::map(out, ~ {
    .x %>%
      ggplot(aes(x = x, y = y, color = label)) +
      geom_point(size = 2, alpha = 0.75) +
      facet_wrap(~cluster, ncol = 3) +
      labs(
        title = glue::glue("Cellular Neighborhoods - {unique(.x$patient)}{ unique(.x$surgery)}"),
        subtitle = data_info
      ) +
      scale_color_manual(values = anno_info$colors) +
      IMCfuncs::facetted_comp_bxp_theme() +
      guides(color = guide_legend(override.aes = list(size = 10))) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0, size = 14, face = "italic"),
      )
  })

  return(out)
}

pdf(
  file = nf("delaunay_cn_surgery_main.pdf", io$outputs$temp_cn),
  width = 20,
  height = 20,
  onefile = TRUE
)
plot_cn_labels(
  spe_obj = lab_spe,
  anno = "main"
)
dev.off()

pdf(
  file = nf("delaunay_cn_surgery_fine.pdf", io$outputs$temp_cn),
  width = 20,
  height = 20,
  onefile = TRUE
)
plot_cn_labels(
  spe_obj = lab_spe,
  anno = "fine"
)
dev.off()

# VISUALISE CELLULAR NEIGHBORHOOD CENTROID FRACTIONS ---------------------------
source("process/Downstream/functions/plot_themes.R")

plot_centroid_fractions <- function(spe_obj,
                                    anno = c("fine", "main"),
                                    cn_label = "delaunay_cn_clusters",
                                    surgery_col = "surgery",
                                    patchwork_ylab = "Cellular Neighborhoods") {
  if (!cn_label %in% names(colData(spe_obj))) stop("cn_label not found in colData(spe_obj)")
  if (!surgery_col %in% names(colData(spe_obj))) stop("surgery col not found in colData(spe_obj)")

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

  cn_prop <- as.data.frame(colData(spe_obj))[, c(anno_info$label, cn_label, surgery_col)]

  cn_label_prop <- cn_prop %>%
    dplyr::select(
      label = !!sym(anno_info$label),
      cn = !!sym(cn_label)
    ) %>%
    dplyr::group_by(label, cn) %>%
    dplyr::summarise(freq = n(), .groups = "drop") %>%
    dplyr::group_by(cn) %>%
    dplyr::mutate(
      cn_total_freq = sum(freq),
      cn_prop = freq / cn_total_freq * 100
    ) %>%
    dplyr::ungroup()


  cn_label_prop <- cn_prop %>%
    dplyr::select(
      label = !!sym(anno_info$label),
      cn = !!sym(cn_label),
      pheno = !!sym(surgery_col)
    ) %>%
    dplyr::group_by(cn, pheno) %>%
    dplyr::summarise(pheno_freq = n(), .groups = "drop") %>%
    dplyr::left_join(cn_label_prop, ., by = "cn", relationship = "many-to-many")

  cn_pheno_counts <- cn_label_prop %>%
    dplyr::select(label, cn, pheno, pheno_freq) %>%
    tidyr::pivot_wider(names_from = pheno, values_from = pheno_freq) %>%
    dplyr::mutate(gap = Rec - Prim) %>%
    dplyr::group_by(cn) %>%
    dplyr::mutate(max = max(Prim, Rec)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(c(Prim, Rec), names_to = "pheno", values_to = "pheno_freq")

  # Left: facetted stacked barplot of cellular neighborhood proportions
  p1 <- cn_label_prop %>%
    ggplot(
      aes(
        x = cn, y = freq,
        fill = forcats::fct_rev(label),
        color = forcats::fct_rev(label)
      )
    ) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      x = patchwork_ylab
    ) +
    coord_flip() +
    facet_wrap(~cn, scales = "free_y", ncol = 1) +
    scale_y_continuous(
      labels = scales::percent,
      minor_breaks = seq(0, 1, 0.1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_fill_manual(
      values = anno_info$colors
    ) +
    scale_color_manual(
      values = anno_info$colors
    ) +
    .horizonal_facet_prop_theme(
      panel_spacing_mm = 10,
      show_y_axis_text = TRUE,
      show_panel_border = FALSE,
      axis.title.y = element_text(size = 25, face = "italic", angle = 90)
    ) +
    guides(
      fill = guide_legend(reverse = TRUE),
      color = guide_legend(reverse = TRUE)
    )

  # Right: facetted dumbell plot showing the cell counts per cellular neighborhood
  p2 <- cn_pheno_counts %>%
    ggplot(aes(x = pheno_freq, y = cn)) +
    geom_line(aes(group = cn), color = "#E7E7E7", linewidth = 10) +
    geom_point(aes(color = pheno), size = 10) +
    geom_text(aes(label = pheno_freq, color = pheno),
      size = 5,
      nudge_x = if_else(cn_pheno_counts$pheno_freq == cn_pheno_counts$max, 275, -275),
      hjust = if_else(cn_pheno_counts$pheno_freq == cn_pheno_counts$max, 0, 1)
    ) +
    facet_wrap(~cn, scales = "free_y", ncol = 1) +
    scale_color_manual(values = spe_obj@metadata$v2_colours$dataset_pheno) +
    .horizonal_facet_prop_theme(
      panel_spacing_mm = 10,
      show_y_axis_text = FALSE,
      show_panel_border = FALSE,
      show_x_axis_text = FALSE,
      show_x_gridlines = FALSE
    )

  return(
    p1 + p2 + plot_layout(guides = "collect")
  )
}

svglite::svglite(
  filename = nf("delaunay_cn_main_anno_props.svg", io$outputs$temp_cn),
  width = 25,
  height = 15
)
plot_centroid_fractions(
  spe_obj = lab_spe,
  anno = "main"
)
dev.off()

svglite::svglite(
  filename = nf("delaunay_cn_fine_anno_props.svg", io$outputs$temp_cn),
  width = 25,
  height = 15
)
plot_centroid_fractions(
  spe_obj = lab_spe,
  anno = "fine"
)
dev.off()

# SPATIAL CONTEXT ANALYSIS -----------------------------------------------------
# The spatial context analysis is a follow-on method to the cellular neighbourhood
# and identifies tissue regions in which distinct cellular neighbourhoods may be interacting.
#
# The method was first introduced and described in:
#
# Bhate, Salil S. et al.Cell Systems (2022), 13, 109-130.e6
#
# The previously defined CN labels will be used to assign cells to spatial
# context for a set of CNs, say CN1,.,CNn if more than 90% of the cells in a
# window of size 100 are assigned to one of those CNs.
#
# This was done using the kNN method in the above paper, however, in our case
# we will use the cellular neighbours obtained using the Delaunay triangulation
# graph, and where the neighbourhoods were defined by the cell fractions.

# drop_cols <- c("knn_cn_neighborhood_agg", "knn_sc", "knn_sc_filt")
# colData(lab_spe) <- colData(lab_spe)[, !colnames(colData(lab_spe)) %in% drop_cols]
# colPair(lab_spe, "k_sc_graph") <- NULL

lab_spe <- aggregateNeighbors(
  object = lab_spe,
  colPairName = "delaunay_50",
  aggregate_by = "metadata",
  count_by = "delaunay_cn_clusters",
  name = "delaunay_cn_neighborhood_agg"
)

lab_spe <- detectSpatialContext(
  object = lab_spe,
  entry = "delaunay_cn_neighborhood_agg",
  threshold = 0.90,
  name = "delaunay_sc"
)

# determine a minimum cell threshold to filter the spatial contexts
cell_threshold <- colData(lab_spe)[, c("patient", "surgery", "ROI")] %>%
  as.data.frame() %>%
  group_by_all() %>%
  dplyr::count() %>%
  dplyr::group_by(surgery) %>%
  dplyr::summarise(min_cells = round(mean(n) * 0.1))

cell_threshold <- setNames(cell_threshold$min_cells, cell_threshold$surgery)

min_sc_patients <- 3

# we will filter the spatial context to include only those contexts which are
# present in at least three separate patients and contain more than
# the minimum cell threshold determined for each surgery type.
prim_filt_sc <- filterSpatialContext(
  object = lab_spe[, lab_spe$surgery == "Prim"],
  entry = "delaunay_sc",
  group_by = "patient",
  group_threshold = min_sc_patients,
  cells_threshold = cell_threshold[["Prim"]],
  name = "prim_delaunay_sc_filt"
)

rec_filt_sc <- filterSpatialContext(
  object = lab_spe[, lab_spe$surgery == "Rec"],
  entry = "delaunay_sc",
  group_by = "patient",
  group_threshold = min_sc_patients,
  cells_threshold = cell_threshold[["Rec"]],
  name = "rec_delaunay_sc_filt"
)

order_unique_scs <- function(vec, index = FALSE) {
  unique_vec <- na.omit(unique(vec))
  split_vec <- strsplit(unique_vec, "_")

  max_len <- max(sapply(split_vec, length))

  padded_matrix <- t(sapply(split_vec, \(x) as.numeric(c(x, rep(NA, max_len - length(x))))))
  padded_matrix <- as.data.frame(padded_matrix)
  padded_matrix$length <- apply(padded_matrix, 1, \(x) length(which(!is.na(x))))

  order_cols_by <- c(
    grep("(?i)^(v1|length)$", colnames(padded_matrix), value = TRUE),
    grep("(?i)^(v1|length)$", colnames(padded_matrix), value = TRUE, invert = TRUE)
  )


  orderd_val <- padded_matrix %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(order_cols_by))) %>%
    dplyr::select(-length) %>%
    tidyr::unite(col = "ordered_val", dplyr::everything(), sep = "_", na.rm = TRUE) %>%
    pull(ordered_val)

  ordered_index <- match(orderd_val, unique_vec)

  if (index) {
    return(ordered_index)
  } else {
    return(orderd_val)
  }
}


# lab_spe$delaunay_sc_filt <- factor(
#     x = lab_spe$delaunay_sc_filt,
#     levels = order_unique_scs(lab_spe$delaunay_sc_filt)
# )
#
# col_sc <- setNames(
# colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(lab_spe$delaunay_sc_filt))),
# levels(lab_spe$delaunay_sc_filt)
# )

prim_filt_sc$prim_delaunay_sc_filt <- factor(
  x = prim_filt_sc$prim_delaunay_sc_filt,
  levels = order_unique_scs(prim_filt_sc$prim_delaunay_sc_filt)
)

rec_filt_sc$rec_delaunay_sc_filt <- factor(
  x = rec_filt_sc$rec_delaunay_sc_filt,
  levels = order_unique_scs(rec_filt_sc$rec_delaunay_sc_filt)
)

sc_colors <- list(
  prim = setNames(
    lab_spe@metadata$v2_colours$cn_colours[str_extract(levels(prim_filt_sc$prim_delaunay_sc_filt), "^\\d+")],
    levels(prim_filt_sc$prim_delaunay_sc_filt)
  ),
  rec = setNames(
    lab_spe@metadata$v2_colours$cn_colours[str_extract(levels(rec_filt_sc$rec_delaunay_sc_filt), "^\\d+")],
    levels(rec_filt_sc$rec_delaunay_sc_filt)
  )
)

# VISUALISE SPATIAL CONTEXTS ---------------------------------------------------
io$outputs$temp_sc <- nd(
  directory_name = "spatial_contexts",
  path = io$outputs$temp_out,
  add_timestamp = FALSE
)

update_geom_defaults("point", list(alpha = 0.75))

plot_sc <- function(spe_obj,
                    spe_filt_col = "surgery",
                    node_label = "delaunay_sc_filt",
                    node_colours = col_sc,
                    graph_node_size = "20",
                    image_id = "sample_id",
                    plot_title = "Filtered Spatial Contexts",
                    graph_method = "delaunay triangulation",
                    min_patient_context = min_sc_patients,
                    min_cells = cell_threshold,
                    plot_subtitle = glue::glue(
                      "graph method: {graph_method}",
                      "min. patients context detected in: {min_patient_context}",
                      "min. cells per context: {min_cells}",
                      .sep = "\n"
                    ),
                    plot_caption = "*NA points did not meet the filtering criteria") {
  if (!spe_filt_col %in% names(colData(spe_obj))) stop("Filter column not found in colData")

  filt_labs <- unique(spe_obj[[spe_filt_col]])

  if (spe_filt_col == "surgery") {
    sample_label <- ifelse(
      length(filt_labs) > 1, "All Samples",
      ifelse(filt_labs == "Prim", "Primary Samples", "Recurrent Samples")
    )
    plot_title <- paste(plot_title, sample_label, sep = " - ")
  } else if (spe_filt_col %in% c("region_type", "region_type_new")) {
    sample_label <- ifelse(
      test = length(filt_labs) > 1,
      yes = "All Regions",
      no = paste0(stringr::str_to_title(unlist(strsplit(filt_labs, "_"))), collapse = "_")
    )

    sample_label <- paste(sample_label, "Regions", sep = " ")

    plot_title <- paste(plot_title, sample_label, sep = " - ")
  } else {
    plot_title <- paste(plot_title, "All Samples", sep = " - ")
  }

  spatial_locs <- imcRtools::plotSpatial(
    object = spe_obj,
    node_color_by = node_label,
    node_size_fix = 2,
    img_id = image_id,
    ncols = 3,
    flip_y = FALSE
  ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = plot_caption
    ) +
    scale_color_manual(values = node_colours) +
    IMCfuncs::facetted_comp_bxp_theme() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0, size = 14, face = "italic"),
      plot.caption = element_text(size = 10, face = "italic")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 10)
      )
    )

  sc_graph_data <- plotSpatialContext(
    object = spe_obj,
    entry = node_label,
    group_by = image_id,
    return_data = TRUE
  )

  sc_graph <- plotSpatialContext(
    object = spe_obj,
    entry = node_label,
    group_by = image_id,
    edge_color_fix = "grey75",
    node_label_color_fix = "black",
    node_color_by = "name",
    node_size_fix = graph_node_size,
  ) +
    scale_color_manual(values = node_colours) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    IMCfuncs::facetted_comp_bxp_theme(legend_key_size = 10) +
    theme(
      # legend.position = "right",
      # legend.title = element_text(size = 20, face = "bold"),
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0, size = 14, face = "italic"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank()
    )

  return(
    list(
      sc_data = sc_graph_data,
      locs = spatial_locs,
      graph = sc_graph
    )
  )
}

sc_plots <- list(
  primary = plot_sc(
    spe_obj = prim_filt_sc,
    spe_filt_col = "surgery",
    node_label = "prim_delaunay_sc_filt",
    node_colours = sc_colors$prim,
    min_cells = cell_threshold[["Prim"]]
  ),
  recurrent = plot_sc(
    spe_obj = rec_filt_sc,
    spe_filt_col = "surgery",
    node_label = "rec_delaunay_sc_filt",
    node_colours = sc_colors$rec,
    min_cells = cell_threshold[["Rec"]]
  )
)

sc_plots <- list_flatten(sc_plots)
sc_graphs <- sc_plots[grep("_data$", names(sc_plots))]
names(sc_graphs) <- str_extract(names(sc_graphs), "(?i)^[a-z]+")

sc_plots <- sc_plots[!grepl("_data$", names(sc_plots))]

pdf(
  file = nf("spatial_context_locs.pdf", io$outputs$temp_sc),
  width = 20,
  height = 25,
  onefile = TRUE
)
print(sc_plots[grep("(primary|recurrent)_locs$", names(sc_plots))])
dev.off()


pdf(
  file = nf("spatial_context_graphs.pdf", io$outputs$temp_sc),
  width = 12,
  height = 9,
  onefile = TRUE
)
print(sc_plots[grep("_graph$", names(sc_plots))])
dev.off()

plot_sc_graph_measures <- function(graph_list, point_colors, title_prefix = "") {
  v <- graph_list$vertices

  g <- igraph::graph_from_data_frame(
    d = graph_list$edges,
    vertices = graph_list$vertices,
    directed = TRUE
  )

  # Calculate various centrality measures
  v$degree_centrality <- igraph::degree(g, mode = "all")
  v$closeness_centrality <- igraph::closeness(g, mode = "all", )
  v$betweenness_centrality <- igraph::betweenness(g, directed = TRUE)

  plot_data <- v %>%
    select(
      sc = 1,
      n_cells = n_cells,
      n_regions = n_group,
      degree_centrality,
      closeness_centrality,
      betweenness_centrality
    ) %>%
    tidyr::pivot_longer(
      cols = -sc,
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(
      measure = factor(measure, levels = c("n_cells", "n_regions", "degree_centrality", "closeness_centrality", "betweenness_centrality"))
    )

  ggplot(plot_data, aes(x = value, y = sc, fill = sc)) +
    geom_point(size = 8, shape = 21, color = "black") +
    facet_wrap(~measure, scales = "free_x", ncol = 1) +
    scale_fill_manual(values = point_colors) +
    IMCfuncs::facetted_comp_bxp_theme() +
    labs(
      title = glue::glue("{title_prefix} Spatial Context Measures"),
      x = "",
      y = ""
    ) +
    theme(legend.position = "none")
}

sc_plots <- imap(sc_graphs, ~ {
  colour_to_plot <- switch(
    EXPR = match.arg(.y, several.ok = FALSE),
    "primary" = sc_colors$prim,
    "recurrent" = sc_colors$rec
  )

  plot_sc_graph_measures(
    graph_list = .x,
    point_colors = colour_to_plot,
    title_prefix = tools::toTitleCase(.y)
  )
})

pdf(
  file = nf("spatial_context_measures.pdf", io$outputs$temp_sc),
  width = 10,
  height = 45,
  onefile = TRUE
)
print(sc_plots)
dev.off()

rm(sc_plots, sc_graphs)

# VISUALISE FILTERED SPATIAL CONTEXT CENTROID FRACTIONS ------------------------
svglite::svglite(
  filename = nf("delaunay_sc_main_anno_props.svg", io$outputs$temp_sc),
  width = 30,
  height = 50
)
plot_centroid_fractions(
  spe_obj = lab_spe[, !is.na(lab_spe$delaunay_sc_filt)],
  anno = "main",
  cn_label = "delaunay_sc_filt",
  patchwork_ylab = "Filtered Spatial Contexts"
)
dev.off()


svglite::svglite(
  filename = nf("delaunay_sc_fine_anno_props.svg", io$outputs$temp_sc),
  width = 30,
  height = 50
)
plot_centroid_fractions(
  spe_obj = lab_spe[, !is.na(lab_spe$delaunay_sc_filt)],
  anno = "fine",
  cn_label = "delaunay_sc_filt",
  patchwork_ylab = "Filtered Spatial Contexts"
)
dev.off()

# CELL INTERACTION ANALYSIS ----------------------------------------------------
# The cell interaction analysis method we will use was first described in:
#
# Schapiro, D. et al. Nat Methods 14, 873–876 (2017).
#
# This method computes the averaged cell type/cell type interaction count and
# compares this count against an empirical null distribution which is generated
# by permuting all cell labels (while maintaining the tissue structure), across
# each sample.
#
# As with the previous analysis we will use the Delaunay triangulation graph to
# define the cell interactions.

io$outputs$temp_ia <- nd(
  directory_name = "cellular_interactions",
  path = io$outputs$temp_out,
  add_timestamp = FALSE
)

library(scales)

source("process/Downstream/functions/cell_Interaction_funcs.R", local = TRUE)

# Test cell interactions across patient/surgery groups
interactions_histocat <- testInteractions(
  object = lab_spe,
  group_by = "patient_surgery",
  label = "manual_gating",
  colPairName = "delaunay_50",
  method = "histocat",
  iter = 1000,
  p_threshold = 0.01,
  BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
)

interactions_histocat_df <- interactions_histocat %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    count_method = "histocat",
    surgery = stringr::str_extract(group_by, "(?i)(prim|rec)"),
    patient = stringr::str_extract(group_by, "^\\d+") %>% stringr::str_c("Patient", ., sep = " ")
  )


# update the labelled spatial object to include the cell states
add_states <- function(spe_obj, state_df) {
  states <- str_split_i(colnames(state_df)[-1], "_", 1) %>% unique()

  for (i in states) {
    spe_obj[[i]] <- NA
    spe_obj[[i]][rownames(colData(spe_obj)) %in% state_df$cell[state_df[[paste0(i, "_high")]]]] <- "high"
    spe_obj[[i]][rownames(colData(spe_obj)) %in% state_df$cell[state_df[[paste0(i, "_low")]]]] <- "low"

    # add the non-NA state labels to each patient/surgery pair
    spe_obj[[i]] <- str_c(spe_obj[["patient_surgery"]], spe_obj[[i]], sep = "_")
  }

  return(spe_obj)
}
lab_spe <- add_states(spe_obj = lab_spe, state_df = state_scores)

state_interactions <- purrr::map(
  list("EMT", "proliferative", "hypoxia", "queiescence"),
  ~ get_state_interactions(spe_obj = lab_spe, state = .x)
)

save_interactions <- list_flatten(state_interactions)
save_interactions$patient_surgery <- interactions_histocat_df

saveRDS(save_interactions, nf("histocat_cell_interactions.rds", io$outputs$temp_ia))

rm(
  interactions_histocat,
  interactions_histocat_df,
  save_interactions,
  state_interactions,
  state_scores
)

# VISUALISE CELL INTERACTION ANALYSIS ------------------------------------------
# load previously saved cell interactions
interacts <- readRDS("outputs/spatial_analysis/2024-09-25T15-19-12/cellular_interactions/histocat_cell_interactions_2024-09-26T15-07-13.rds")

patients <- split(interacts$patient_surgery, interacts$patient_surgery$patient)

patients <- purrr::imap(patients, ~ {
  delta_surgery_interactions(
    ia_df = .x,
    patient_min = 1,
    title_suffix = .y
  )
})
patients <- purrr::imap(patients, ~ plot_ia(ia_clean_list = .x))

all <- delta_surgery_interactions(interacts$patient_surgery)
all <- plot_ia(ia_clean_list = all)

pdf(
  file = nf("patient_surgery.pdf", io$outputs$temp_ia),
  width = 15,
  height = 15,
  onefile = TRUE
)

print(all)
print(patients)
dev.off()

rm(all, patients)

# remove the patient_surgery interactions and just keep the state interactions
interacts <- interacts[-grep("patient_surgery", names(interacts))]

states_min_1 <- purrr::imap(interacts, ~ {
  delta_surgery_interactions(
    ia_df = .x,
    patient_min = 1,
    title_suffix = str_replace(.y, "_", " ")
  )
})
states_min_1 <- purrr::imap(states_min_1, ~ plot_ia(ia_clean_list = .x))

states_min_2 <- purrr::imap(interacts, ~ {
  delta_surgery_interactions(
    ia_df = .x,
    patient_min = 2,
    title_suffix = str_replace(.y, "_", " ")
  )
})
states_min_2 <- purrr::imap(states_min_2, ~ plot_ia(ia_clean_list = .x))

states_min_3 <- purrr::imap(interacts, ~ {
  delta_surgery_interactions(
    ia_df = .x,
    patient_min = 3,
    title_suffix = str_replace(.y, "_", " ")
  )
})
states_min_3 <- purrr::imap(states_min_3, ~ plot_ia(ia_clean_list = .x))

pdf(
  file = nf("min_patient=1_states.pdf", io$outputs$temp_ia),
  width = 15,
  height = 15,
  onefile = TRUE
)
print(states_min_1)
dev.off()

pdf(
  file = nf("min_patient=2_states.pdf", io$outputs$temp_ia),
  width = 15,
  height = 15,
  onefile = TRUE
)
print(states_min_2)
dev.off()

pdf(
  file = nf("min_patient=3_states.pdf", io$outputs$temp_ia),
  width = 15,
  height = 15,
  onefile = TRUE
)
print(states_min_3)
dev.off()

rm(states_min_1, states_min_2, states_min_3, interacts)

# SAVE DATA --------------------------------------------------------------------
saveRDS(lab_spe, nf("lab_spe.rds", io$outputs$temp_out))

# END --------------------------------------------------------------------------
