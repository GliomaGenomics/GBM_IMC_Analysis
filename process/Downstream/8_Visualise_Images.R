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

# spe <- readRDS(io$inputs$spe_data)

lab_spe <- readRDS(io$inputs$lab_spe_data)
lab_spe <- lab_spe[channelNames(images), ]

images <- images[unique(lab_spe$sample_id)]
masks <- masks[unique(lab_spe$sample_id)]

interactions <- readRDS(io$inputs$interactions)

interactions <- interactions$patient_surgery

# CELL VISUALSATION ------------------------------------------------------------
cur_id <- c("64Prim_001")

cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]


plotCells(cur_masks,
  object = lab_spe,
  cell_id = "cell_id",
  img_id = "sample_id",
  # outline_by = "manual_gating",
  # colour_by = "manual_gating",
  colour_by = "delaunay_cn_clusters",
  background_colour = "black",
  missing_colour = "black",
  # colour = list(
  #     manual_gating = lab_spe@metadata$v2_colours$cells
  # ),
  image_title = NULL,
  scale_bar = list(
    length = 100,
    label = expression("100 " ~ mu * "m"),
    cex = 2,
    lwidth = 10,
    colour = "white",
    position = "bottomright",
    margin = c(10, 10)
  ),
  legend = NULL,
  thick = TRUE
)

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

# VISUALISING VASCULAR/OLIGODENDROCYTES ----------------------------------------

# normalise the channel intensities
# norm_images = cytomapper::normalize(cur_images,
#                                     separateChannels = TRUE,
#                                     separateImages = FALSE,
#                                     inputRange = NULL
#                                     )
# clip channel intensities
# norm_images = normalize(norm_images, inputRange = c(0, 0.2))

png(
  filename = nf("64_Rec.png", io$output$image_visualisation$paper_2),
  units = "in",
  width = 30,
  height = 15,
  res = 600
)

plotPixels(

  cur_images,
  colour_by = c(
    "MOG",
    "SMA"
  ),
  bcg = list(
    MOG = c(0, 5, 1),
    SMA = c(0, 10, 1)
  ),
  colour = list(
    SMA = c("black", "red2"),
    MOG = c("black", "cyan")
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
