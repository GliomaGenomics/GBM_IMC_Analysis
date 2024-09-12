# This script needs to be re-factored from the cell-type prevalences onwards -
# The intention is that all of the downstream analysis will be split into distinct
# files and this script will serve as a downstream preparation and initial
# exploratory analysis script.

# PACKAGES ---------------------------------------------------------------------
library(SpatialExperiment)
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
    comp_data = "data/downstream/compensated"
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

# END----
