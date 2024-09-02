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

# VISUALISE LABELLED CELLS COORDINATES -----------------------------------------
lab_spe <- spe[, spe$ROI %in% c("001", "002", "003") & !is.na(spe$manual_gating)]

lab_spe <- as_tibble(spatialCoords(lab_spe)) %>%
  dplyr::rename(x = Pos_X, y = Pos_Y) %>%
  cbind(
    id = lab_spe$sample_id,
    patient = lab_spe$patient,
    surgery = lab_spe$surgery,
    label = lab_spe$manual_gating
  )

svglite::svglite(
  filename = nf("labelled_samples.svg", io$outputs$temp_out),
  width = 25,
  height = 20
)
lab_spe %>%
  ggplot(aes(x = x, y = y, color = label)) +
  geom_point(size = 1, alpha = 0.75) +
  facet_wrap(~id) +
  scale_color_manual(values = spe@metadata$v2_colours$cells) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 10)))
dev.off()


svglite::svglite(
  filename = nf("labelled_sugery.svg", io$outputs$temp_out),
  width = 27,
  height = 10
)
lab_spe %>%
  ggplot(aes(x = x, y = y, color = label)) +
  geom_point(size = 1, alpha = 0.75) +
  facet_grid(surgery ~ patient) +
  scale_color_manual(values = spe@metadata$v2_colours$cells) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 10)))
dev.off()

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

# CREATE SPATIAL INTERACTION GRAPHS --------------------------------------------
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

outplots <- vector("list", length = length(colPairNames(lab_spe)))
names(outplots) <- colPairNames(lab_spe)

for (i in seq_along(colPairNames(lab_spe))) {
  colpair <- colPairNames(lab_spe)[[i]]

  outplots[[i]] <- plotSpatial(lab_spe[, lab_spe$patient %in% c("64") & lab_spe$surgery %in% c("Prim", "Rec")],
    node_color_by = "manual_gating",
    img_id = "sample_id",
    draw_edges = TRUE,
    colPairName = colpair,
    nodes_first = FALSE,
    ncols = 3,
    edge_color_fix = "grey"
  ) +
    ggtitle(glue::glue("{colpair} interaction graph")) +
    scale_color_manual(values = spe@metadata$v2_colours$cells) +
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

pdf(
  file = nf("knn_sweep_graphs.pdf", io$outputs$temp_out),
  width = 20,
  height = 15,
  onefile = TRUE
)

print(outplots)

dev.off()

# SAVE DATA --------------------------------------------------------------------
# END --------------------------------------------------------------------------
