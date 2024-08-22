# This script needs to be re-factored from the cell-type prevalences onwards -
# The intention is that all of the downstream analysis will be split into distinct
# files and this script will serve as a downstream preparation and initial
# exploratory analysis script.

# PACKAGES ---------------------------------------------------------------------
library(SpatialExperiment)
# library(dittoSeq)
# library(imcRtools)
library(ggplot2)
library(viridis)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(IMCfuncs)

# I/O --------------------------------------------------------------------------
io <- list(
  inputs = list(
    comp_data = "data/downstream/compensated",
    functions = "process/Downstream/functions"
  ),
  output = list(
    cell_pheno_out = "outputs/cell_phenotyping"
  ),
  plots = list()
)

# create the output directories (if they do not exist)
ndirs(io$output$cell_pheno_out)
io$output$temp_out <- nd(path = io$output$cell_pheno_out)

# use the latest compensated, labelled data
io$inputs$comp_spe <- list.files(io$inputs$comp_data,
  pattern = "spe_comp_[T0-9-]+",
  full.names = T
)

io$inputs$comp_spe <- io$inputs$comp_spe[order(io$inputs$comp_spe, decreasing = TRUE)][[1]]

# LOAD LABELLED DATA  ----------------------------------------------------------
spe <- readRDS(io$inputs$comp_spe)

# ADD LABELS AND COLOURS -------------------------------------------------------
io$inputs$colors <- list(
  dataset_pheno = c(
    Prim = "#2166AC", Rec = "#B2182B", up = "#004529", down = "#FF7F00"
  ),
  cell_groups = c(
    Immune = "#1AE4B6FF", Cancer = "#30123BFF", Normal = "#FABA39FF", Vasculature = "#7A0403FF"
  ),
  cells = c(
    `T cell` = "#0a5746",
    `NK cell` = "#129e7e",
    Microglia = "#1AE4B6FF",
    Macrophage = "#a6f5e3",
    AC = "#30123BFF",
    MES = "#9e3bc2",
    NPC = "#d2a4e3",
    OPC = "#f0e0f5",
    Neuron = "#b97d05",
    Astrocyte = "#FABA39FF",
    Oligodendrocyte = "#fcd586",
    Endothelial = "#7A0403FF"
  ),
  diverging = c(
    "#446faa", "#FFFFFF", "#BB4444"
  ),
  positive = c(
    "#D1E6EF", "#ABC4DE", "#9099CA", "#8566B1", "#762D81", "#540046"
  ),
  imc_visual = c(
    "#0000FF", "#00FF00", "#FF0000", "#FF00FF", "#00FFFF", "#FFFF00", "#FFFFFF", "#FFA500"
  )
)

io$inputs$labels <- list()

io$inputs$labels$markers <- list(
  Immune = c("CD45", "CD3", "NKP46", "IBA1", "TMEM119"),
  Cancer = c(
    "SLC1A3_EAAT1", "HOPX", "SOD2", "CHI3L1", "ANEXIN_A2", "ANXA1",
    "DLL3", "BCAN", "SCD5", "OLIG1"
  ),
  Normal = c("NeuN_FOX3", "GFAP", "MOG"),
  Vasculature = c("CD31", "SMA")
)

io$inputs$labels$cell_types <- list(
  Immune = c("T cell", "NK cell", "Macrophage", "Microglia"),
  Cancer = c("AC", "MES", "NPC", "OPC"),
  Normal = c("Neuron", "Astrocyte", "Oligodendrocyte"),
  Vasculature = c("Endothelial")
)


spe$main_anno_v2 <- ifelse(
  spe$manual_gating %in% io$inputs$labels$cell_types$Immune, "Immune",
  ifelse(
    spe$manual_gating %in% io$inputs$labels$cell_types$Cancer, "Cancer",
    ifelse(
      spe$manual_gating %in% io$inputs$labels$cell_types$Normal, "Normal",
      ifelse(
        spe$manual_gating %in% io$inputs$labels$cell_types$Vasculature, "Vasculature", NA
      )
    )
  )
)


spe$main_anno_v2 <- factor(spe$main_anno_v2,
  levels = names(io$inputs$labels$cell_types)
)

# LABELLED CELL DIMPLOTS  ------------------------------------------------------
dimplot_data <- reducedDim(spe, "UMAP_harmony_minmax") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  mutate(
    main_anno = spe$main_anno_v2,
    fine_anno = spe$manual_gating
  ) %>%
  dplyr::rename(
    "UMAP1" = V1,
    "UMAP2" = V2
  ) %>%
  mutate(
    across(
      ends_with("_anno"),
      ~ forcats::fct_na_value_to_level(.x, "Unassigned")
    )
  )


main_anno <- dimplot_data %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = main_anno)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(
    values = io$inputs$colors$cell_groups,
    na.value = "grey80"
  ) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(
    legend.position = "right"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  )

fine_anno <- dimplot_data %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = fine_anno)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(
    values = io$inputs$colors$cells,
    na.value = "grey80"
  ) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(
    legend.position = "right"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  )

svglite::svglite(
  file = nf("main_anno_dimplot.svg", io$output$temp_out),
  width = 15,
  height = 12
)
print(main_anno)
dev.off()

svglite::svglite(
  file = nf("fine_anno_dimplot.svg", io$output$temp_out),
  width = 15,
  height = 12
)
print(fine_anno)
dev.off()

rm(main_anno, fine_anno, dimplot_data)

# LABELLED CELL HEATMAP  -------------------------------------------------------
label_exprs <- scuttle::summarizeAssayByGroup(
  x = spe,
  ids = spe$manual_gating,
  assay.type = "zscore",
  subset.row = unlist(io$inputs$labels$markers, use.names = FALSE),
  subset.col = NULL,
  statistics = "median",
  store.number = "ncells"
) %>%
  SummarizedExperiment::assay("median")


heatmap_data <- as.data.frame(label_exprs) %>%
  tibble::rownames_to_column(var = "marker") %>%
  tidyr::pivot_longer(
    cols = -marker,
    names_to = "cell_type",
    values_to = "score"
  ) %>%
  mutate(
    across(
      cell_type,
      ~ factor(.x, levels = colnames(label_exprs))
    ),
    across(
      marker,
      ~ factor(.x, levels = rev(rownames(label_exprs)))
    )
  )

# Label counts boxplot
count_data <- spe$manual_gating[spe$manual_gating %in% unlist(io$inputs$labels$cell_types)] %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(cell_type = 1, count = Freq)

p1 <- count_data %>%
  ggplot(aes(x = cell_type, y = count)) +
  geom_bar(stat = "identity", fill = "slateblue", width = 0.8) + # Use `stat = "identity"` to plot actual values
  geom_text(aes(label = count), vjust = -0.5, size = 5, fontface = "bold") + # Add labels on top of bars
  ylim(0, max(count_data$count) * 1.15) +
  IMCfuncs::facetted_comp_bxp_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()
  )


heatmap_theme <-
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(color = "grey50", linewidth = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = ggplot2::unit(10, "mm"),
    panel.grid.major = element_blank()
  )



p2 <- heatmap_data %>%
  ggplot(mapping = aes(x = cell_type, y = marker)) +
  geom_tile(aes(fill = score), color = "grey25", linetype = 3, na.rm = TRUE) +
  ggplot2::scale_fill_gradientn(
    colors = c("grey99", "grey99", "#313695"),
    limits = c(-2, 2.5),
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = "grey90",
      frame.colour = "grey90",
      barwidth = 2.5,
      barheight = 15
    ),
    name = "z-score",
    na.value = "white"
  ) +
  ggplot2::labs(
    title = "Marker Expression",
    subtitle = "Median expression of all labelled cells"
  ) +
  heatmap_theme


add_highlight_regions <- function(baseplot,
                                  row_groups,
                                  col_groups,
                                  highlight_colors,
                                  highlight_width = 3,
                                  legend_name = "") {
  unique_row_groups <- unique(row_groups)
  unique_col_groups <- unique(col_groups)

  highlight_regions <- data.frame(
    group = factor(unique_row_groups, levels = names(highlight_colors)),
    xmin = NA,
    xmax = NA,
    ymin = NA,
    ymax = NA,
    color = highlight_colors[unique_row_groups]
  )

  for (i in seq_along(unique_row_groups)) {
    # Calculate ymin and ymax for the row groups (in reverse order)
    ymin_position <- min(which(row_groups == unique_row_groups[i])) # top-most position
    ymax_position <- max(which(row_groups == unique_row_groups[i])) # bottom-most position

    highlight_regions$ymin[i] <- length(row_groups) - ymax_position + 0.5
    highlight_regions$ymax[i] <- length(row_groups) - ymin_position + 1.5

    highlight_regions$xmin[i] <- min(which(col_groups == unique_col_groups[i])) - 0.5
    highlight_regions$xmax[i] <- max(which(col_groups == unique_col_groups[i])) + 0.5
  }

  out_plot <- baseplot +
    geom_rect(
      data = highlight_regions,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
      fill = NA, linewidth = highlight_width, inherit.aes = FALSE
    ) +
    scale_color_manual(name = legend_name, values = highlight_colors)

  return(out_plot)
}

p2 <- add_highlight_regions(
  baseplot = p2,
  row_groups = rep(
    x = names(io$inputs$labels$cell_types),
    times = lengths(io$inputs$labels$markers)
  ),
  col_groups = rep(
    x = names(io$inputs$labels$markers),
    times = lengths(io$inputs$labels$cell_types)
  ),
  highlight_colors = io$inputs$colors$cell_groups,
  legend_name = ""
)

library(patchwork)

svglite::svglite(
  file = nf("label_expression_heatmap.svg", io$output$temp_out),
  width = 12,
  height = 15
)

p1 / p2 + plot_layout(heights = c(1, 7))

dev.off()

rm(
  label_exprs, count_data, heatmap_data,
  add_highlight_regions, heatmap_theme,
  p1, p2
)

# SAVE DOWNSTREAM DATA ---------------------------------------------------------
spe@metadata$v2_colours <- io$inputs$colors
spe@metadata$labels <- io$inputs$labels

saveRDS(spe, file = nf("spe_downstream.rds", io$output$temp_out))

# REGION-TYPE CELL LABELS   ----------------------------------------------------
count_anno_labels <- function(grouping_var, anno_label,
                              outcol = "cells",
                              remove_na = TRUE) {
  out_data <- as.data.frame(spe@colData) %>%
    dplyr::count(!!sym(anno_label),
      group = !!sym(grouping_var),
      name = outcol
    )

  if (remove_na) out_data <- out_data %>% tidyr::drop_na()

  return(out_data)
}

anno_labels <- c("main_anno", "manual_gating")

counts <- lapply(anno_labels, \(x) count_anno_labels(
  grouping_var = "sample_id",
  anno_label = x,
  remove_na = TRUE
))

names(counts) <- anno_labels

# Add SOX2 and HIF1A to the proliferating group

# cell_states = IMCfuncs::combine_vectors(
#   list_of_vectors = cell_states,
#   combined_vec = "high_prolif",
#   unique_vals = TRUE,
#   "high_proliferation",
#   "proliferating_stem_cell"
# )

cell_states <- spe@metadata$cell_states

spe$high_prolif <- FALSE
spe@colData[cell_states$high_proliferation, "high_prolif"] <- TRUE
table(spe$high_prolif)

spe$high_hypoxia <- FALSE
spe@colData[cell_states$high_hypoxia, "high_hypoxia"] <- TRUE
table(spe$high_hypoxia)

anno_labels <- c("high_prolif", "high_hypoxia")

counts$state_anno <- lapply(anno_labels, \(x) count_anno_labels(
  grouping_var = "sample_id",
  anno_label = x,
  remove_na = FALSE
))


counts$state_anno <- lapply(counts$state_anno, \(x){
  state_anno <- colnames(x)[1]

  out <- x %>%
    filter(.[[1]] == TRUE) %>%
    mutate(state_anno = state_anno) %>%
    select(state_anno, group, cells)

  return(out)
})

counts$state_anno <- bind_rows(counts$state_anno)

total_cells_per_sample <- as.data.frame(spe@colData) %>%
  count(sample_id, name = "total_cells") %>%
  rename(group = sample_id)

counts <- map(counts, ~ {
  left_join(.x, total_cells_per_sample, by = "group") %>%
    mutate(percent = round(cells / total_cells * 100)) %>%
    rename_with(~"anno_label", .cols = 1)
})


# REGION-TYPE PLOT DATA  -------------------------------------------------------
plot_data <- bind_rows(counts$main_anno, counts$state_anno)

plot_data <- plot_data %>%
  filter(!anno_label %in% c("Tumour", "Stroma")) %>%
  mutate(across(anno_label, tolower)) %>%
  arrange(group)

plot_data$anno_label <- gsub("(?i)^high_", "", plot_data$anno_label)

plot_data <- as.data.frame(spe@colData) %>%
  select(sample_id, region_type) %>%
  distinct() %>%
  rename(group = sample_id) %>%
  left_join(plot_data, ., by = "group")

# RE-LABEL THE REGION-TYPES   --------------------------------------------------

region_counts <- split(plot_data, as.factor(plot_data$group))

# Logic for assigning the region labels
#
# Get all labels with prevalences >= min_percent:
#   for any labels meet this criteria:
#       assign labels as any prevalences within 10% of the highest label
#   if none found:
#        assign "unknown"


region_annotate <- function(x,
                            immune_min = 20,
                            state_min = 25,
                            within = 10,
                            label_col = "anno_label",
                            percent_col = "percent") {
  stopifnot(exprs = {
    state_min >= 0 & immune_min <= 100
    immune_min >= 0 & immune_min <= 100
    within >= 1 & within <= 99
  })

  out <- x[order(x[[percent_col]], decreasing = T), ][c(label_col, percent_col)]
  out <- setNames(out[[percent_col]], out[[label_col]])

  region_labels <- vector("character", 0)

  if (out[["immune"]] >= immune_min) region_labels[[1]] <- "immune"

  out <- out[out >= state_min]
  if (length(out) == 0) {
    return("unknown")
  }

  region_out <- append(region_labels, names(out[out >= out[1] * (1 - within)]))

  return(unique(region_out))
}


min_state <- 25
min_immune <- 20

labels <- lapply(region_counts, \(x) region_annotate(
  x,
  immune_min = min_immune,
  state_min = min_state,
  within = 10
))

labels <- unlist(sapply(labels, \(x) paste0(x, collapse = "_")))

plot_data <- left_join(
  plot_data,
  data.frame(group = names(labels), new_label = labels),
  by = "group"
)

new_labels <- strsplit(plot_data$new_label, "_")

new_labels <- lapply(new_labels, \(x){
  if (length(x) > 1) out <- paste0(sort(x), collapse = "_") else out <- x
  return(out)
})

plot_data$new_label <- unlist(new_labels, use.names = FALSE)

plot_data <- plot_data %>%
  mutate(facet_label = glue::glue("{group}\n{region_type} --> {new_label}"))

spe$region_type_new <- NULL

spe$region_type_new <- plot_data %>%
  rename(sample_id = group) %>%
  group_by(sample_id) %>%
  summarise(region_type_new = unique(new_label)) %>%
  left_join(as.data.frame(spe@colData), ., by = "sample_id") %>%
  pull(region_type_new)


region_counts <- plot_data %>%
  group_by(group) %>%
  summarise(region_type_new = unique(new_label)) %>%
  ungroup() %>%
  mutate(surgery = str_extract(group, "Prim|Rec")) %>%
  mutate(across(surgery, ~ ifelse(.x %in% "Prim", "Primary", "Recurrent"))) %>%
  count(region_type_new, surgery)

# PLOT RE-LABELLED REGION TYPES   ----------------------------------------------

plot_theme <- theme_classic() +
  ggplot2::theme(
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(
      size = 32, hjust = 0.5,
      face = "bold", colour = "black",
      margin = margin(25, 0, 10, 0)
    ),
    plot.subtitle = element_text(
      size = 25, hjust = 0.5,
      face = "italic", colour = "black",
      margin = margin(10, 0, 50, 0)
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 20, colour = "black"),

    # axis.title.x = element_blank(),
    strip.text = element_text(size = 24, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "#EBFFFF"),
    panel.spacing.y = unit(2.5, "lines"),
    legend.text = element_text(size = 25, colour = "black", face = "bold"),
    legend.key.size = unit(3, "line")
  )


svglite::svglite(
  filename = nf("region_annotations.svg", io$output$cell_pheno),
  width = 12,
  height = 250
)

ggplot(
  plot_data,
  aes(x = anno_label, y = percent, group = anno_label, fill = anno_label)
) +
  geom_col() +
  facet_wrap("facet_label", ncol = 1) +
  ggtitle("Region Type Estimation",
    subtitle = glue::glue("Minimum prevelances (immune={min_immune}% | cell states={min_state}%)")
  ) +
  scale_fill_manual(
    values = spe@metadata$color_vectors$region_type,
    breaks = unique(plot_data$anno_label)
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  plot_theme

dev.off()

# Update colours and save the spatial object
spe@metadata$color_vectors$region_type

spe@metadata$color_vectors$region_type_new <- set_names(
  viridis::viridis(length(unique(spe$region_type_new))),
  unique(spe$region_type_new)
)


svglite::svglite(
  filename = nf("region_counts.svg", io$output$cell_pheno),
  width = 15,
  height = 10
)

region_counts %>%
  ggplot(aes(x = region_type_new, y = n, group = surgery, fill = surgery)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("") +
  ylab("Number of Regions") +
  scale_fill_manual(values = setNames(
    c("#4F94CD", "#CD4F39"),
    c("Primary", "Recurrent")
  )) +
  scale_y_continuous(breaks = 1:10) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank())

dev.off()

# SAVE OBJECT AND EXPRESSION MATRIX --------------------------------------------
outfile <- sort(
  list.files(io$inputs$comp_data,
    pattern = "spe_comp_\\d{4}",
    full.names = T
  ),
  decreasing = T
)[1]

saveRDS(spe, file = outfile)

# Create an abundance matrix from a spatial object
cell_exprs <- as.data.frame(t(assay(spe, "exprs")))

expression_markers <- getmarkers(c("Immune", "Tumour", "Stroma", "cell_states"))
names(expression_markers) <- NULL

cell_exprs <- cell_exprs %>% tibble::rownames_to_column(var = "cells_id")
cell_exprs <- cbind(as.data.frame(spe@colData), cell_exprs)


cell_exprs <- cell_exprs %>% select(
  sample_id,
  patient,
  surgery,
  ROI,
  region_type,
  region_type_new,
  responder_type,
  cell_id,
  cell_pheno = main_anno,
  cell_anno = manual_gating,
  high_prolif,
  high_hypoxia,
  dplyr::all_of(expression_markers)
)

cell_exprs$surgery <- ifelse(cell_exprs$surgery == "Prim", "primary", "recurrent")
cell_exprs$cell_pheno <- as.character(cell_exprs$cell_pheno)
cell_exprs$cell_anno <- as.character(cell_exprs$cell_anno)


readr::write_csv(
  x = cell_exprs,
  file = nf("IMC_cell_abundance.csv", io$output$cell_abundances)
)


# CELL TYPE PREVALENCES ----
rm(list = grep("^spe$|^io$", ls(), value = T, invert = T))

# For the first bit of exploratory analysis we can check the cell type prevalence
# across various stratifications, starting with the the specific tissue samples.
# In order to do this we will first obtain counts for each of the cells and then
# subset/aggregate them to show the cell prevalence across different groupings.

# In our previous analysis we clustered and annotated the imaging mass optometry
# single cells. During this process a number of cells remained unannotated as their
# identity could not be easily inferred using one of more of the cell type marker
# which were used. Also, in some cases the abundance of the metal isotopes corresponding
# to the cell type markers was indiscriminate in distinguishing one particular cell
# type or cell state and so was left "unknown".
# This initial cell type prevalence analysis will include these "unknown" markers.

# Creating a data.frame of prevalences
source("../Scripts/Functions/Spatial_Analysis/calculate_cell_prevalances.R")

prevelances <- calculate_cell_prevalences(spe,
  anno_field = "manual_gating",
  id_field = "sample_id",
  remove_undefined = F
)

# Adding metadata to the prevalences
prevelances <- prevelances %>%
  mutate(across(sample_id, ~ str_replace(.x, "^(\\d{2})", "\\1_"))) %>%
  tidyr::separate(
    col = sample_id, into = c("patient", "surgery", "region"),
    sep = "_", remove = FALSE
  ) %>%
  mutate(across(sample_id, ~ str_replace(.x, "_", ""))) %>%
  mutate(patient_id = paste0(patient, surgery)) %>%
  relocate(patient_id, .before = sample_id)

prevelances <- as.data.frame(spe@colData) %>%
  select(sample_id, responder_type,
    region_type = region_type_new
  ) %>%
  distinct() %>%
  left_join(prevelances, ., by = "sample_id")


# Adding high-level annotations information
prevelances$high_level_anno <-
  ifelse(prevelances$cell_type %in% unlist(spe@metadata$cell_types$immune), "Immune",
    ifelse(prevelances$cell_type %in% unlist(spe@metadata$cell_types$tumour), "Tumour",
      ifelse(prevelances$cell_type %in% unlist(spe@metadata$cell_types$stroma), "Stroma",
        "Undefined"
      )
    )
  )


prevelances <- prevelances %>% relocate(high_level_anno, .after = cell_type)

# Plotting prevalences
io$plots <- list()

source("../Scripts/Functions/Spatial_Analysis/plot_cell_prevalences.R")

io$plot_params <- data.frame(
  x_val = rep(
    c("surgery", "region_type", "responder_type"),
    c(10, 10, 5)
  ),
  fill = "cell_type",
  facet_val = rep(c(NA, "patient", NA, "surgery", NA), each = 5),
  show = rep(c("all", "labelled", "immune", "tumour", "stroma"),
    times = 5
  )
)

# Add the high-level annotation params
io$plot_params <- rbind(
  data.frame(
    x_val = rep(
      c("surgery", "region_type", "responder_type"),
      c(4, 4, 2)
    ),
    fill = "high_level_anno",
    facet_val = rep(c(NA, "patient", NA, "surgery", NA), each = 2),
    show = rep(c("all", "labelled"))
  ),
  io$plot_params
)

io$plot_params

names(spe@metadata$color_vectors$cell_anno)[names(spe@metadata$color_vectors$cell_anno) %in% c("T cells", "NK cells")] <- c("T cell", "NK cell")


spatial_anno_cols <- c(
  spe@metadata$color_vectors$main_anno,
  spe@metadata$color_vectors$cell_anno
)


# plot_cell_prevalences(prevelances,
#                       x_val = "surgery",
#                       y_val = "count",
#                       fill = "high_level_anno",
#                       colours = spatial_anno_cols,
#                       show = "all"
#                       )


io$plots <- pmap(io$plot_params,
  .f = plot_cell_prevalences,
  df = prevelances,
  colours = spatial_anno_cols
)

names(io$plots) <- tidyr::unite(io$plot_params,
  col = "plots",
  sep = "_",
  na.rm = T
) %>% pull(plots)

pdf(
  onefile = T, width = 15, height = 15,
  file = IMCfuncs:::new_file("cell_prevalences", io$output$cell_prevalences)
)

io$plots

dev.off()


source("../Scripts/Functions/Spatial_Analysis/aggregate_cell_prevalances.R")

aggregate_prevalences(
  prevalence_df = prevelances,
  main_group = "surgery",
  seconday_group = "patient",
  cell_counts_col = "count",
  cell_labels_col = "high_level_anno",
  show = "labelled"
)

# Obtain aggregate counts
io$cell_counts <- io$plot_params %>%
  rename(
    main_group = x_val,
    seconday_group = facet_val,
    cell_labels_col = fill
  ) %>%
  pmap(.f = aggregate_prevalences, prevalence_df = prevelances)


names(io$cell_counts) <- tidyr::unite(io$plot_params,
  col = "plots",
  sep = "_",
  na.rm = T
) %>%
  pull(plots) %>%
  str_replace_all("high_level_anno", "low") %>%
  str_replace_all("cell_type", "high") %>%
  str_replace_all("region_type", "regions") %>%
  str_replace_all("responder_type", "responders")


openxlsx::write.xlsx(io$cell_counts,
  file = file.path(io$out, "cell_counts.xlsx")
)


# QUANTIFYING CELL TYPE PREVALENCES ----

# we can quantify the cell type prevalences stratified across various groups in
# order to see if any particular cell types are significantly different across
# comparisons groups. The identities of the stratifications are as follows:
#
# surgery (primary and recurrent)
# region_type (hypoxia, proliferative and immune)
# responder_type (up and down)


# Sub-setting the prevalences to include only the columns which I intend to compare
prevelances <- prevelances %>%
  select(
    patient, surgery, region,
    responder_type, region_type,
    high_level_anno, cell_type, count
  )



split_prevalences <- function(data, group_by = NULL,
                              return_plot = F,
                              split_surgery = F) {
  if (is.null(group_by)) stop("group_by not supplied")

  cols <- unique(data[[group_by]])

  if (split_surgery) {
    message(glue::glue("grouping on:\t patient, region, surgery, cell_type and {group_by}"))

    split_data <- data %>%
      group_by(patient, region, surgery, .data[[group_by]], cell_type) %>%
      summarize(count = sum(count), .groups = "keep") %>%
      ungroup() %>%
      tidyr::pivot_wider(names_from = group_by, values_from = count) %>%
      mutate(across(all_of(cols), ~ ifelse(is.na(.), 0, .))) %>%
      tidyr::unite(region, patient, region) %>%
      split(.$cell_type)
  } else {
    message(glue::glue("grouping on:\t patient, region, cell_type and {group_by}"))

    split_data <- data %>%
      group_by(patient, region, .data[[group_by]], cell_type) %>%
      summarize(count = sum(count), .groups = "keep") %>%
      ungroup() %>%
      tidyr::pivot_wider(names_from = group_by, values_from = count) %>%
      mutate(across(all_of(cols), ~ ifelse(is.na(.), 0, .))) %>%
      tidyr::unite(region, patient:region) %>%
      split(.$cell_type)
  }

  if (return_plot) {
    plots <- purrr::imap(split_data, ~ {
      n_fill_colors <- ncol(.x) - 2

      outplot <- .x %>%
        tidyr::pivot_longer(
          cols = -c("region", "cell_type"),
          names_to = "group", values_to = "count"
        ) %>%
        ggplot(aes(x = group, y = count, fill = group)) +
        geom_boxplot(alpha = .3, outlier.shape = NA) +
        geom_jitter(width = .1, size = 3, alpha = .5) +
        theme_classic() +
        theme(
          plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 16)
        ) +
        ggtitle(glue::glue("{.y}")) +
        scale_fill_manual(values = viridis::turbo(n = n_fill_colors, begin = .2))

      return(outplot)
    })

    out <- list(data = split_data, plots = plots)

    return(out)
  } else {
    return(split_data)
  }
}


split_data <- setNames(
  c("surgery", "region_type", "responder_type"),
  c("surgery", "region_type", "responder_type")
) |>
  as.list() |>
  purrr::map(~ split_prevalences(
    data = prevelances,
    group_by = .x, return_plot = T
  ))


# region_type split by surgery
prim_region_type <- split_prevalences(prevelances, "region_type",
  return_plot = F, split_surgery = T
)


split_data$prim_region_type <- list(data = map(prim_region_type, ~ {
  .x %>%
    filter(surgery %in% "Prim") %>%
    select(-surgery)
}))
names(split_data$prim_region_type$data) <- paste0(names(split_data$prim_region_type$data), "_Prim")


split_data$rec_region_type <- list(data = map(prim_region_type, ~ {
  .x %>%
    filter(surgery %in% "Rec") %>%
    select(-surgery)
}))
names(split_data$rec_region_type$data) <- paste0(names(split_data$rec_region_type$data), "_Rec")


plot_surgery_split_data <- function(.x, .y) {
  n_fill_colors <- ncol(.x) - 2

  plot_title <- unlist(str_split(.y, "_"))

  plot_title <- glue::glue("{plot_title[1]} ({plot_title[2]})")

  outplot <- .x %>%
    tidyr::pivot_longer(
      cols = -c("region", "cell_type"),
      names_to = "group", values_to = "count"
    ) %>%
    ggplot(aes(x = group, y = count, fill = group)) +
    geom_boxplot(alpha = .3, outlier.shape = NA) +
    geom_jitter(width = .1, size = 3, alpha = .5) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_text(size = 16)
    ) +
    ggtitle(plot_title) +
    scale_fill_manual(values = viridis::turbo(n = n_fill_colors, begin = .2))

  return(outplot)
}

split_data$prim_region_type$plots <- imap(split_data$prim_region_type$data, ~ plot_surgery_split_data(.x, .y))

split_data$rec_region_type$plots <- imap(split_data$rec_region_type$data, ~ plot_surgery_split_data(.x, .y))


# responder_type split by surgery
prim_responder_type <- split_prevalences(prevelances, "responder_type",
  return_plot = F, split_surgery = T
)


split_data$prim_responder_type <- list(data = map(prim_responder_type, ~ {
  .x %>%
    filter(surgery %in% "Prim") %>%
    select(-surgery)
}))
names(split_data$prim_responder_type$data) <- paste0(names(split_data$prim_responder_type$data), "_Prim")


split_data$rec_responder_type <- list(data = map(prim_responder_type, ~ {
  .x %>%
    filter(surgery %in% "Rec") %>%
    select(-surgery)
}))
names(split_data$rec_responder_type$data) <- paste0(names(split_data$rec_responder_type$data), "_Rec")


split_data$prim_responder_type$plots <- imap(split_data$prim_responder_type$data, ~ plot_surgery_split_data(.x, .y))

split_data$rec_responder_type$plots <- imap(split_data$rec_responder_type$data, ~ plot_surgery_split_data(.x, .y))



compare_two_groups <- function(data) {
  cols <- grep("region|cell_type", colnames(data), invert = T, value = T)

  if (length(cols) != 2) stop(glue::glue("{length(cols)} groups in data - only 2 allowed"))

  # Paired Wilcoxon Signed Rank Tests
  test <- wilcox.test(
    x = data[[cols[1]]],
    y = data[[cols[2]]],
    paired = T,
    conf.int = T, # with CIs
    exact = F # gets rid of obnoxious warnings
  )

  effect_data <- data %>%
    tidyr::pivot_longer(all_of(cols)) %>%
    arrange(name, region)

  effect <- rcompanion::wilcoxonPairedR(
    x = effect_data$value,
    g = effect_data$region, ci = T
  )

  out <- data.frame(
    celltype = unique(data$cell_type),
    comparion = paste0(cols, collapse = "_vs_"),
    pval = test$p.value,
    effect_size = abs(effect$r)
  )

  return(out)
}

compare_two_plus_groups <- function(data) {
  cols <- grep("region|cell_type", colnames(data), invert = T, value = T)

  comps <- combn(cols, 2, simplify = F)

  comp_stats <- function(data, comps) {
    # Paired Wilcoxon Signed Rank Tests
    test <- wilcox.test(
      x = data[[comps[1]]],
      y = data[[comps[2]]],
      paired = T,
      conf.int = T, # with CIs
      exact = F # gets rid of obnoxious warnings
    )

    effect_data <- tidyr::pivot_longer(
      data = data[c("region", "cell_type", comps)],
      cols = all_of(comps)
    ) %>%
      arrange(name, region)


    effect <- rcompanion::wilcoxonPairedR(
      x = effect_data$value,
      g = effect_data$region
    )

    out <- data.frame(
      celltype = unique(data$cell_type),
      comparion = paste0(comps, collapse = "_vs_"),
      pval = test$p.value,
      effect_size = abs(effect)
    )

    return(out)
  }

  out <- purrr::map(comps, ~ comp_stats(data = data, comps = .x)) |> bind_rows()

  return(out)
}


## The Wilcoxon Signed Rank Test can be used to see if the population median of
## the difference scores is equal to 0 or not

# We can also check the effect size measure for this test, using the
# wilcoxonPairedR() function from the rcompanion package.
# The function requires the data to be ordered specifically: The top half of the
# data are the values for group1. The bottom half of the data are the values for group2.
# Then the id variable needs to be in the same order for both groups.
#
# An effect size above 0.7 is considered to be a large,
# Anything above 0.5 is a moderate effect size,
# and anything above 0.3 is a small effect size.

# compare_two_groups(split_data$prim_responder_type$data$AC_Prim)
# compare_two_plus_groups(split_data$prim_region_type$data$AC_Prim)

split_data$surgery$stats <- purrr::map(
  split_data$surgery$data,
  ~ compare_two_groups(.x)
) |> bind_rows()

split_data$responder_type$stats <- purrr::map(
  split_data$responder_type$data,
  ~ compare_two_groups(.x)
) |> bind_rows()


split_data$region_type$stats <- purrr::map(
  split_data$region_type$data,
  ~ compare_two_plus_groups(.x)
) |> bind_rows()


split_data$prim_region_type$stats <- purrr::map(
  split_data$prim_region_type$data,
  ~ compare_two_plus_groups(.x)
) |> bind_rows()


split_data$rec_region_type$stats <- purrr::map(
  split_data$rec_region_type$data,
  ~ compare_two_plus_groups(.x)
) |> bind_rows()


split_data$prim_responder_type$stats <- purrr::map(
  split_data$prim_responder_type$data,
  ~ compare_two_groups(.x)
) |> bind_rows()


split_data$rec_responder_type$stats <- purrr::map(
  split_data$rec_responder_type$data,
  ~ compare_two_groups(.x)
) |> bind_rows()


io$prevalances <- list()

io$prevalances$stats <- map(split_data, ~ {
  prevelance_stats <- pluck(.x, "stats")

  prevelance_stats$pval_adj <- p.adjust(prevelance_stats$pval, method = "fdr")

  return(prevelance_stats)
})

io$prevalances$data <- map(split_data, ~ {
  pluck(.x, "data")
})

io$prevalances$plots <- map(split_data, ~ {
  pluck(.x, "plots")
})

# Save the prevalence data
purrr::iwalk(io$prevalances$data, ~ {
  filename <- glue::glue("{.y}_prevalence_data.xlsx")

  openxlsx::write.xlsx(.x,
    file = file.path(
      "Outputs/Spatial_Analysis/Prevalences",
      filename
    )
  )
})

# Save the prevalence statistics
openxlsx::write.xlsx(io$prevalances$stats,
  file = file.path(
    io$output$cell_prevalences,
    "prevalence_stats.xlsx"
  )
)

# Save the prevalence pairwise comparison plots
pdf(file = "Outputs/Spatial_Analysis/Prevalences/prevalence_comps.pdf", onefile = T)

io$prevalances$plots

dev.off()

# clear memory
rm(
  split_data, compare_two_groups, compare_two_plus_groups, split_prevalences,
  prim_region_type, plot_surgery_split_data, prim_responder_type
)

io$prevalances <- NULL
io$plot_params <- NULL


# Test the normality assumption
# qqnorm(foo$diff)
# qqline(foo$diff, col = "steelblue", lwd = 2) # bit better
#
# # test for nomality
# shapiro.test(foo$diff)



ggsave(
  filename = file.path(io$out, "Surgery_cell_differences.pdf"),
  plot = gridExtra::marrangeGrob(io$diff_plots, nrow = 1, ncol = 1),
  width = 10, height = 10
)



# SPATIAL INTERACTION GRAPHS ----

# The spatial graphs are stored in colPair(spe, name) slots.
# These slots store SelfHits objects representing edge lists in which the first
# column indicates the index of the “from” cell and the second column the index
# of the “to” cell. Each edge list is newly constructed when sub-setting the object.
# colPair(spe, "neighborhood") stores the spatial graph constructed by steinbock
#
# During the image processing step a neighbourhood graph was generated based on
# object centroids using dmax = 15 pixels . Aside from this there are also a number
# of other metrics which we can use to estimate clustering or interaction frequencies
# between cell types.

# Determining a suitable value for K neighbors
# Each image will comprise on a varying number of cells and so a suitable K value
# need to be established to identify large and small cluster. As an initial check
# we can see how many cell we have in each image to see what k may be appropriate.


for (i in seq(5, 25, 5)) {
  spe <- buildSpatialGraph(spe,
    img_id = "ImageName", type = "knn", k = i,
    name = paste0("k_", i, collapse = "")
  )
}

# spe <- buildSpatialGraph(spe, img_id = "ImageName", name = "k_25", type = "knn", k = 25)
# spe <- buildSpatialGraph(spe, img_id = "ImageName", type = "expansion", threshold = 20)
# spe <- buildSpatialGraph(spe, img_id = "ImageName", type = "delaunay", max_dist = 50)

colPairNames(spe)


# SPATIAL INTERACTION VISUALISATION ----

pdf(file = "Outputs/Spatial_Analysis/Interactions/Interaction_graphs.pdf")

spe@colData[c("ImageName", "manual_anno")] %>%
  as.data.frame() %>%
  group_by(ImageName) %>%
  summarise(cells = n()) %>%
  arrange(cells) %>%
  mutate(name = factor(ImageName, levels = ImageName)) %>%
  ggplot(aes(x = cells, y = name)) +
  geom_col(fill = "slateblue") +
  ggtitle("Total Cells per Image") +
  ylab("Image") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


spe@colData[c("patient_id_extn", "manual_anno")] %>%
  as.data.frame() %>%
  group_by(patient_id_extn) %>%
  summarise(cells = n()) %>%
  arrange(cells) %>%
  mutate(name = factor(patient_id_extn, levels = patient_id_extn)) %>%
  ggplot(aes(x = cells, y = name)) +
  geom_col(fill = "slateblue") +
  ggtitle("Total Cells per Patient ID") +
  ylab("Patient id") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

dev.off()

io$plots$all_spatial <- plotSpatial(spe,
  node_color_by = "manual_anno",
  img_id = "sample_id",
  node_size_fix = 0.5
) +
  scale_color_manual(values = spe@metadata$color_vectors$v2_all_anno) +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave("Outputs/Spatial_Analysis/Interactions/all_spatial.pdf",
  plot = io$plots$all_spatial,
  device = "pdf",
  width = 35, height = 30, units = "cm", limitsize = FALSE
)


plot_interactions <- function(interaction, sample, color_by = "manual_anno",
                              colors = spe@metadata$color_vectors$v2_all_anno) {
  if (!sample %in% unique(spe@colData$ImageName)) stop("sample not found in object")

  if (!color_by %in% names(spe@colData)) stop("color_by not found in colData()")

  if (!interaction %in% colPairNames(spe)) stop("interaction not found in colPairNames()")

  if (interaction == "neighborhood") {
    plot_title <- glue::glue("Centriod dmax = 15 ({sample})")
  } else {
    plot_title <- glue::glue(str_replace(interaction, "_", " = "), " ({sample})")
  }

  plotSpatial(spe[, spe$ImageName == sample],
    node_color_by = color_by,
    img_id = "sample_id",
    draw_edges = TRUE,
    colPairName = interaction,
    nodes_first = FALSE,
    edge_color_fix = "grey"
  ) +
    scale_color_manual(values = colors) +
    ggtitle(plot_title) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.title = element_blank()
    )
}

pdf("Outputs/Spatial_Analysis/Interactions/test.pdf",
  width = 15, height = 12
)

purrr::map(colPairNames(spe), ~ plot_interactions(.x, sample = "64Rec_002"))

purrr::map(colPairNames(spe), ~ plot_interactions(.x, sample = "82Prim_003"))

dev.off()


# INTERACTION ANALYSIS ----

# sampling the data for testing purposes
# spe_sample <- spe[,sample(1:ncol(spe), ncol(spe)* 0.05)]

# spe$all_v2 <- ifelse(is.na(spe$manual_anno),
#                       "undefined", spe$manual_anno) |> as.factor()
#
# spe$responder_type <- ifelse(spe$patient_id == "71",
#                              "down", "up") |> as.factor()


# surgery_interactions <- test_interactions_2(
#   object = spe_sample,
#   group_by = "surgery",
#   label = "all_v2",
#   method = "histocat",
#   p_threshold = 0.05,
#   colPairName = "k_15",
#   BPPARAM = BiocParallel::SerialParam(RNGseed = 221029)
# ) |> as.data.frame()



interactions <- set_names(
  c("surgery", "responder_type", "region"),
  c("surgery", "responder_type", "region")
) %>%
  as.list() %>%
  purrr::map(~ testInteractions(
    object = spe,
    group_by = .x,
    label = "all_v2",
    method = "histocat",
    p_threshold = 0.05,
    colPairName = "k_15",
    BPPARAM = BiocParallel::SerialParam(RNGseed = 1234)
  ))

interactions_2 <- map(interactions, as.data.frame)

openxlsx::write.xlsx(
  interactions_2,
  "Outputs/Spatial_Analysis/Interactions/Interaction_stats.xlsx"
)

interaction_heatmap <- function(interactions_tbl) {
  split_groups <- interactions_tbl %>%
    group_split(group_by)

  out_plots <- purrr::map(split_groups, ~ {
    plot_title <- unique(.x$group_by)

    .x %>%
      group_by(from_label, to_label) %>%
      summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "keep") %>%
      ggplot() +
      geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
      scale_fill_gradient2(
        low = scales::muted("blue"),
        mid = "white",
        high = scales::muted("red")
      ) +
      ggtitle(plot_title) +
      xlab("Cell phenotype in neighborhood") +
      ylab("Cell phenotype of interest") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_text(size = 14, face = "italic")
      )
  })

  return(out_plots)
}

io$plots$heatmaps <- purrr::map(interactions_2, ~ interaction_heatmap(.x)) |> purrr::flatten()

pdf("Outputs/Spatial_Analysis/Interactions/interaction_heatmaps.pdf",
  width = 10, height = 10
)

io$plots$heatmaps

dev.off()


# SPATIAL VISUALISATION ----

# Steinbock interaction graph
plotSpatial(spe[, spe$sample_id == "64Prim_001"],
  node_color_by = "all_anno",
  img_id = "sample_id",
  draw_edges = TRUE,
  colPairName = "neighborhood",
  nodes_first = FALSE,
  node_size_fix = 2.5,
  edge_color_fix = "darkgrey"
) +
  scale_color_manual(values = metadata(spe)$color_vectors$all_anno) +
  ggtitle("steinbock interaction graph") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  )

# knn interaction graph
plotSpatial(spe[, spe$sample_id == "64Prim_001"],
  node_color_by = "all_anno",
  img_id = "sample_id",
  draw_edges = TRUE,
  colPairName = "knn_interaction_graph",
  nodes_first = FALSE,
  node_size_fix = 2.5,
  edge_color_fix = "darkgrey"
) +
  scale_color_manual(values = metadata(spe)$color_vectors$all_anno) +
  ggtitle("knn interaction graph") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  )

# expansion interaction graph
plotSpatial(spe[, spe$sample_id == "64Prim_001"],
  node_color_by = "all_anno",
  img_id = "sample_id",
  draw_edges = TRUE,
  colPairName = "expansion_interaction_graph",
  nodes_first = FALSE,
  node_size_fix = 2.5,
  edge_color_fix = "darkgrey"
) +
  scale_color_manual(values = metadata(spe)$color_vectors$all_anno) +
  ggtitle("expansion interaction graph") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  )

# delaunay interaction graph
plotSpatial(spe[, spe$sample_id == "64Prim_001"],
  node_color_by = "all_anno",
  img_id = "sample_id",
  draw_edges = TRUE,
  colPairName = "delaunay_interaction_graph",
  nodes_first = FALSE,
  node_size_fix = 2.5,
  edge_color_fix = "darkgrey"
) +
  scale_color_manual(values = metadata(spe)$color_vectors$all_anno) +
  ggtitle("delaunay interaction graph") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  )


p <- plotSpatial(spe,
  node_color_by = "all_anno",
  img_id = "sample_id",
  node_size_fix = 0.5
) +
  scale_color_manual(values = metadata(spe)$color_vectors$all_anno) +
  theme(legend.title = element_blank())

ggsave(
  plot = p,
  filename = "Outputs/Spatial_Analysis/All_samples.svg",
  width = 20, height = 18, units = "in"
)


# SPATIAL COMMINITY ANALYSIS ----
library(igraph)
set.seed(220819)

# Spatial community detection - Primary samples
prim_spe <- spe[, spe$surgery == "Prim"]

gr <- graph_from_data_frame(as.data.frame(colPair(prim_spe, "neighborhood")),
  directed = FALSE,
  vertices = data.frame(index = seq_len(ncol(prim_spe)))
)

cl_comm <- cluster_louvain(gr)
comm_tumor <- paste0("Primary_", membership(cl_comm))
comm_tumor[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_tumor) <- colnames(prim_spe)


plotSpatial(spe[, spe$celltype == "Tumor"],
  node_color_by = "spatial_community",
  img_id = "sample_id",
  node_size_fix = 0.5
) +
  theme(legend.position = "none") +
  ggtitle("Spatial tumor communities") +
  scale_color_manual(values = rev(colors()))


# Spatial community detection - Recurrent samples
rec_spe <- spe[, spe$surgery == "Rec"]

gr <- graph_from_data_frame(as.data.frame(colPair(rec_spe, "neighborhood")),
  directed = FALSE,
  vertices = data.frame(index = seq_len(ncol(rec_spe)))
)

cl_comm <- cluster_louvain(gr)
comm_stroma <- paste0("Recurrent_", membership(cl_comm))
comm_stroma[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_stroma) <- colnames(rec_spe)


comm <- c(comm_tumor, comm_stroma)
spe$spatial_community <- comm[colnames(spe)]


plotSpatial(prim_spe,
  node_color_by = "spatial_community",
  img_id = "sample_id",
  node_size_fix = 0.5
) +
  theme(legend.position = "none") +
  ggtitle("Spatial Primary communities") +
  scale_color_manual(values = rev(colors()))



# END----
