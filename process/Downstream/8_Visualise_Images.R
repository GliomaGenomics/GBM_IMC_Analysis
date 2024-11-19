# Author: Shoaib Ajaib
# Date: 15/10/2024

# USAGE ------------------------------------------------------------------------
# This script is designed to visualise the spatial analysis regions of interest

# OPTIONS ----------------------------------------------------------------------
# options(scipen = 999)

# PACKAGES ---------------------------------------------------------------------
library(viridis)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(IMCfuncs)
library(SpatialExperiment)
library(cytomapper)

# I/O --------------------------------------------------------------------------
io <- list(
  inputs = list(
    spe_dir = "outputs/cell_phenotyping",
    images = "data/downstream/compensated/images_comp.rds",
    masks = "data/downstream/raw/masks.rds",
    interactions = "outputs/spatial_analysis/2024-09-25T15-19-12/cellular_interactions/histocat_cell_interactions_2024-09-26T15-07-13.rds"
  ),
  outputs = list(
    out_dir = "outputs/ROI_images"
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

# get the most recent spe data
io$inputs$spe_data <- find_file(io$inputs$spe_dir, file_pattern = "spe_downstream")

# get the most recent labelled spe data
io$inputs$lab_spe_data <- list.files(
  path = "outputs/spatial_analysis/2024-09-05T12-01-52",
  pattern = "lab_spe",
  full.names = TRUE
)

rm(find_file)

# LOAD DATA --------------------------------------------------------------------
images <- readRDS(io$inputs$images)

masks <- readRDS(io$inputs$masks)
names(masks) <- str_extract(names(masks), "(?i)[a-z0-9]+_[0-9]+")
masks@elementMetadata$sample_id <- str_extract(masks@elementMetadata$sample_id, "(?i)[a-z0-9]+_[0-9]+")

spe <- readRDS(io$inputs$spe_data)

lab_spe <- readRDS(io$inputs$lab_spe_data)
lab_spe <- lab_spe[channelNames(images), ]

images <- images[unique(lab_spe$sample_id)]
masks <- masks[unique(lab_spe$sample_id)]

interactions <- readRDS(io$inputs$interactions)

interactions <- interactions$patient_surgery

# IDENTIFY THE PATIENTS WHERE CELLS ARE INTERACING -----------------------------
signif_patients <- function(interact_df = interactions,
                            surgery_filt = c("Prim", "Rec"),
                            from_cell,
                            to_cell) {
  if (missing(from_cell) | missing(to_cell)) {
    cli::cli_abort("Please provide the from_cell and to_cell arguments!")
  }

  vector_list <- vector("list", length(to_cell))
  names(vector_list) <- to_cell

  for (i in seq_along(to_cell)) {
    vector_list[[i]] <- interact_df %>%
      dplyr::filter(surgery %in% surgery_filt) %>%
      dplyr::filter(sigval >= 1) %>%
      dplyr::filter(from_label %in% from_cell) %>%
      dplyr::filter(to_label %in% to_cell[[i]]) %>%
      dplyr::arrange(p) %>%
      dplyr::pull(group_by)
  }

  common_elements <- Reduce(intersect, vector_list)

  if (length(common_elements) == 0) {
    cli::cli_alert_danger("No common patients!")
    cli::cat_line()
    return(vector_list)
  } else {
    cli::cli_alert_success("Common patients:\t{common_elements}")
    cli::cat_line()
    return(vector_list)
  }
}

# 1. ENDOTHELIAL CELLS INTERACTIONS --------------------------------------------
signif_patients(
  surgery_filt = "Prim",
  from_cell = "Endothelial",
  to_cell = c("Microglia", "MES")
)
cur_id <- grep("^71Prim", names(images), value = TRUE)
cell_types <- c("Endothelial", "Microglia", "MES")


signif_patients(
  from_cell = "Endothelial",
  to_cell = c("Macrophage")
)
cur_id <- grep("^71Rec", names(images), value = TRUE)
cell_types <- c("Endothelial", "Macrophage")

# 2. OLIGODENDROCYTES INTERACTIONS ---------------------------------------------
signif_patients(
  surgery_filt = "Prim",
  from_cell = "Oligodendrocyte",
  to_cell = c("MES")
)

cur_id <- grep("^64Prim", names(images), value = TRUE)
cell_types <- c("Oligodendrocyte", "MES")

signif_patients(
  surgery_filt = "Rec",
  from_cell = "Oligodendrocyte",
  to_cell = c("T cell", "Microglia", "Endothelial")
)

cur_id <- grep("^64Rec", names(images), value = TRUE)
cell_types <- c("Oligodendrocyte", "T cell", "Microglia", "Endothelial")


# 3. ASTROCYTE INTERACTIONS ----------------------------------------------------
signif_patients(
  surgery_filt = "Prim",
  from_cell = "Astrocyte",
  to_cell = c("Microglia", "Macrophage")
)

cur_id <- grep("^84Prim", names(images), value = TRUE)
cell_types <- c("Astrocyte", "Microglia", "Macrophage")



signif_patients(
  surgery_filt = "Rec",
  from_cell = "Astrocyte",
  to_cell = c("OPC", "NPC")
)

cur_id <- grep("^71Rec", names(images), value = TRUE)
cell_types <- c("Astrocyte", "OPC", "NPC")

# 4. NEURON INTERACTIONS -------------------------------------------------------
signif_patients(
  surgery_filt = "Prim",
  from_cell = "Neuron",
  to_cell = c("Macrophage")
)

cur_id <- grep("^82Prim", names(images), value = TRUE)
cell_types <- c("Neuron", "Macrophage")


signif_patients(
  surgery_filt = "Rec",
  from_cell = "Neuron",
  to_cell = c("Astrocyte", "MES")
)

cur_id <- grep("^82Rec", names(images), value = TRUE)
cell_types <- c("Neuron", "Astrocyte", "MES")
# CELL VISUALSATION ------------------------------------------------------------
# c(
#   "T cell", "NK cell", "Macrophage", "Microglia", "AC", "MES", "NPC", "OPC",
#   "Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"
# )

cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

cell_colours <- c("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", "#FFFF00")
names(cell_colours) <- cell_types

filt_spe <- lab_spe[, lab_spe$manual_gating %in% cell_types]
filt_spe$outline_col <- filt_spe$manual_gating

# plotting the cell masks coloured by specific cells
plotCells(cur_masks,
  object = filt_spe,
  cell_id = "cell_id",
  img_id = "sample_id",
  outline_by = "outline_col",
  colour_by = "manual_gating",
  background_colour = "black",
  missing_colour = "grey20",
  colour = list(
    manual_gating = cell_colours
  ),
  display = "single",
  image_title = NULL,
  scale_bar = NULL,
  thick = FALSE,
  save_plot = list(
    filename = nf(
      filename = paste0(unique(str_extract(cur_id, "(?i)^[a-z0-9]+")), ".png"),
      filepath = io$output$temp_out
    ),
    # filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)
# image_title  = list(
#     text = mcols(cur_images)$sample_id,
#     colour = "white",
#     margin = c(10, 10),
#     position = "topleft",
#     font = 2,
#     cex = 4
# ),
# scale_bar = list(
#     length = 100,
#     label = expression("100" ~ mu *"m"),
#     colour = "white",
#     lwidth = 5,
#     position = "bottomright",
#     margin = c(75, 20),
#     frame = length(cur_id),
#     cex = 3
# ),
# margin = 40

# INTERACTIVE IMAGE VISUALSATION -----------------------------------------------

library(cytoviewer)

app <- cytoviewer(
  image = images,
  mask = masks,
  object = lab_spe,
  cell_id = "cell_id",
  img_id = "sample_id"
)

if (interactive()) {
  shiny::runApp(app)
}

# 1. ENDOTHELIAL PIXELS --------------------------------------------------------
cur_id <- "71Prim_003"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "SMA",
    "SOD2"
  ),
  bcg = list(
    SMA = c(0, 5, 1),
    SOD2 = c(0, 5, 1)
  ),
  colour = list(
    SMA = c("black", "red2"),
    SOD2 = c("black", "#00FF00")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

cur_id <- "71Rec_002"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "SMA",
    "IBA1"
  ),
  bcg = list(
    SMA = c(0, 5, 1),
    IBA1 = c(0, 5, 1)
  ),
  colour = list(
    SMA = c("black", "red2"),
    IBA1 = c("black", "#0000FF")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

# 2. OLIGODENDROCYTE PIXELS ----------------------------------------------------
cur_id <- "64Prim_001"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "MOG",
    "SOD2",
  ),
  bcg = list(
    MOG = c(0, 5, 1),
    SOD2 = c(0, 5, 1)
  ),
  colour = list(
    MOG = c("black", "red2"),
    SOD2 = c("black", "#00FF00")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)


cur_id <- "64Rec_001"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "MOG",
    "CD3",
    "SMA"
  ),
  bcg = list(
    MOG = c(0, 6, 1),
    CD3 = c(0, 6, 1),
    SMA = c(0, 5, 1)
  ),
  colour = list(
    MOG = c("black", "red2"),
    CD3 = c("black", "#0000FF"),
    SMA = c("black", "#FF00FF")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

# 3. ASTROCYTE PIXELS ----------------------------------------------------------

cur_id <- "82Prim_002"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "GFAP",
    "IBA1"
  ),
  bcg = list(
    GFAP = c(0, 3, 1),
    IBA1 = c(0, 5, 1)
  ),
  colour = list(
    GFAP = c("black", "red2"),
    IBA1 = c("black", "#0000FF")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

cur_id <- "71Rec_001"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "GFAP",
    "SCD5",
    "BCAN"
  ),
  bcg = list(
    GFAP = c(0, 3, 1),
    SCD5 = c(0, 2, 1),
    BCAN = c(0, 8, 1)
  ),
  colour = list(
    GFAP = c("black", "red2"),
    SCD5 = c("black", "#00FFFF"),
    BCAN = c("black", "#FF00FF")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

# 4. NEURON PIXELS -------------------------------------------------------------
cur_id <- "82Prim_002"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]


plotPixels(
  cur_images,
  colour_by = c(
    "NeuN_FOX3",
    "IBA1"
  ),
  bcg = list(
    NeuN_FOX3 = c(0, 10, 1),
    IBA1 = c(0, 2, 1)
  ),
  colour = list(
    NeuN_FOX3 = c("black", "red2"),
    IBA1 = c("black", "#00FF00")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)



cur_id <- "82Rec_001"
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

plotPixels(
  cur_images,
  colour_by = c(
    "NeuN_FOX3",
    "GFAP",
    "SOD2"
  ),
  bcg = list(
    NeuN_FOX3 = c(0, 9, 1),
    GFAP = c(0, 2, 1),
    SOD2 = c(0, 5, 1)
  ),
  colour = list(
    NeuN_FOX3 = c("black", "red2"),
    GFAP = c("black", "#0000FF"),
    SOD2 = c("black", "#00FFFF")
  ),
  image_title = NULL,
  scale_bar = NULL,
  display = "single",
  interpolate = TRUE,
  save_plot = list(
    filename = nf(glue::glue("{cur_id}.png"), io$output$temp_out),
    scale = 3
  )
)

# VISUALISING ASTROCYTES/NEURONS -----------------------------------------------

cur_id <- c(
  # "82Prim_002"
  "82Rec_002"
)

cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

# normalise the channel intensities
# norm_images = cytomapper::normalize(cur_images,
#                                     separateChannels = TRUE,
#                                     separateImages = FALSE,
#                                     inputRange = NULL
#                                     )
# # clip channel intensities
# norm_images = normalize(norm_images, inputRange = c(0, 0.4))


png(
  filename = nf("82_rec.png", io$output$image_visualisation$paper_2),
  units = "in",
  width = 30,
  height = 15,
  res = 600
)

plotPixels(

  cur_images,
  colour_by = c(
    "NeuN_FOX3",
    "GFAP",
    "SMA"
  ),
  bcg = list(
    GFAP = c(0, 3, 1),
    NeuN_FOX3 = c(0, 2, 1),
    SMA = c(0, 10, 1)
  ),
  colour = list(
    SMA = c("black", "red2"),
    GFAP = c("black", "magenta4"),
    NeuN_FOX3 = c("black", "green")
  ),
  legend = NULL,
  image_title = list(
    text = mcols(cur_images)$`patient_id.V1`,
    colour = "white",
    margin = c(10, 10),
    position = "topleft",
    font = 2,
    cex = 5
  ),
  scale_bar = list(
    length = 100,
    label = expression("100 " ~ mu * "m"),
    colour = "white",
    lwidth = 10,
    position = "bottomright",
    margin = c(75, 20),
    frame = length(cur_id),
    cex = 5
  ),
  margin = 40
)

dev.off()



# VISUALISING NEOPLASTIC MARKERS ----

# SOD2:MES-LIKE, DLL3:NPC-LIKE, OLIG1:OPC-LIKE

cur_id <- c("67Prim_001")
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]


# normalise the channel intensities
norm_images <- normalize(cur_images)

# clip channel intensities
norm_images <- normalize(norm_images, inputRange = c(0, 0.175))

tiff("Outputs/neoplastic_markers.tiff", units = "in", width = 10, height = 10, res = 600)

plotPixels(norm_images,
  return_plot = FALSE,
  return_images = FALSE,
  scale_bar = list(
    length = 100,
    label = "",
    colour = "white",
    cex = 1
  ),
  legend = NULL,
  image_title = NULL,
  display = "single",
  colour_by = c("SOD2", "DLL3", "OLIG1"),
  colour = list(
    SOD2 = c("black", "green"),
    DLL3 = c("black", "blue"),
    OLIG1 = c("black", "red")
  ),
)
dev.off()

# VISUALISING VARIOUS OTHER MARKERS ----

# KI67:proliferating cells, IBA1: pan macrophage, TNC:Quiecent stem cell


cur_id <- c("64Prim_003")
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]


# normalise the channel intensities
norm_images <- normalize(cur_images)

# clip channel intensities
norm_images <- normalize(norm_images, inputRange = c(0, 0.15))


tiff("Outputs/various_markers.tiff", units = "in", width = 10, height = 10, res = 600)

plotPixels(norm_images,
  return_plot = FALSE,
  return_images = FALSE,
  scale_bar = list(
    length = 100,
    label = "",
    colour = "white",
    cex = 1
  ),
  legend = NULL,
  image_title = NULL,
  display = "single",
  colour_by = c("Ki67", "TNC", "IBA1"),
  colour = list(
    Ki67 = c("black", "red"),
    TNC = c("black", "green"),
    IBA1 = c("black", "blue")
  ),
)

dev.off()

# END --------------------------------------------------------------------------
