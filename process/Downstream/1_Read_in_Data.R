# This script takes the single-cell and image data and generates the single cell
# data as a SingeCellExperiment object, as well as reading in the images and segmentation masks.

# PACKAGES ---------------------------------------------------------------------

library(imcRtools)
library(cytomapper)
library(dplyr)
library(stringr)
library(magrittr)


# I/O --------------------------------------------------------------------------

io <- list(
  
  inputs = list(
    functions = "process/Downstream/functions",
    data = "data",
    fullstack_images = file.path("data/fullstacks"),
    segmentation_masks = file.path("data/masks/tiff")),
  
  output= list(raw_data_out = "data/downstream/raw")
)


if(!dir.exists(io$output$raw_data_out)) dir.create(io$output$raw_data_out)


# READ IN SINGLE CELL DATA -----------------------------------------------------
# this will exclude the small regions which were run as part of the 67Prim samples.
img_to_read <- "(?i)^(64|67|71|82|84)(Prim|Rec)_00[1-7]{1}"

spe <- imcRtools::read_steinbock(
  path = io$input$data,
  intensities_folder = "intensities",
  regionprops_folder = "regionprops",
  graphs_folder = NULL,
  pattern = img_to_read,
  extract_cellid_from = "Object",
  extract_coords_from = c("centroid-1", "centroid-0"),
  image_file = "panel/image_info.csv",
  extract_imagemetadata_from = c("region_type", "width_px", "height_px"),
  panel_file = "panel/panel.csv",
  extract_names_from = "name",
  return_as = c("spe", "sce")
)


# Summarized pixel intensities per channel and cell (here mean intensity) 
spe@assays@data$counts[1:5,1:5]

# Metadata associated to individual cells
head(spe@colData)

# Spatial locations of all cells
head(SpatialExperiment::spatialCoords(spe))

# Channel metadata
head(rowData(spe))


# ADD SINGLE CELL METADATA -----------------------------------------------------
# unique column names
names(colData(spe))[2] = "cell_id"

colnames(spe) <- paste0(spe$sample_id, "_", spe$cell_id)

# Patient
spe$patient = unlist(
  stringr::str_extract_all(spe@colData$sample_id, pattern = "^\\d{2}"))

# Surgery
spe$surgery = unlist(
  stringr::str_extract_all(spe@colData$sample_id, pattern = "(?i)prim|rec"))

# Region of interest 
spe$ROI = unlist(
  stringr::str_extract_all(spe@colData$sample_id, pattern = "\\d{3}$"))

# Responder type
spe$responder_type = "up"
spe$responder_type = ifelse(spe@colData$patient == 71, "down", "up")



# Changing the colData column order
coldata_order = c("sample_id","patient", "surgery","ROI", "cell_id", 
                  "region_type", "responder_type", "area", "eccentricity",
                  "axis_major_length","axis_minor_length", "width_px", "height_px") 

if(length(coldata_order) != ncol(colData(spe))) stop("coldata_order and colData ncol is not the same")

colData(spe) = colData(spe)[coldata_order]

colData(spe)

# Add colours
spe@metadata$colour_vectors = list()

spe@metadata$color_vectors$all_anno <- setNames(

  c("#E69F00","#009E73","#CC79A7","#BD7FF8",
    "#0072B2","#56B4E9","#F0E442","#D55E00","#AD7700",
    "#30123BFF","#1AE4B6FF","#A2FC3CFF","#7A0403FF",
    "gray75"),

  c("AC", "MES", "NPC", "OPC",
    "Monocytes","Macrophages","Microglia", "NK Cells", "T Cells",
    "Neuron", "Oligodendrocyte", "Astrocyte", "Vasculature",
    "Undefined")
  )


spe@metadata$color_vectors$high_level_anno <- setNames(
  c("#440154FF", "#21908CFF", "#FDE725FF","gray75"),
  c("Immune", "Neoplastic","Stroma", "Undefined")
  )



# READ IN IMAGES ---------------------------------------------------------------
# source(file.path(io$inputs$functions, "IMC_Imaging.R"))
# 
# img_to_read <- c(64, 67, 71, 82, 84)
# img_to_read <- paste0("(?i)", img_to_read, "(Prim|Rec)_00[1-7]{1}[_a-z]+\\.tiff$")
# 
# process_images <- function(x){
#   
#   # Read in the 32-bit multi-channel Tiff images
#   images <- cytomapper::loadImages(x = io$inputs$fullstack_images, pattern = x)
#   
#   # Set multichannel image names for visualization
#   channelNames(images) <- rownames(spe)
#   
#   images <- cytomapper::getChannels(images, 
#                                     which(channelNames(images) != "DNA1"))
#   
#   return(images)
#   
# }
# 
# images <- lapply(img_to_read, process_images) 
# 
# images <- do.call(c, images)
# 
# # Check if the images were read in correctly 
# check_images(images)
# 
# # Remove the DNA1 channel(channel 34)
# spe <- spe[rownames(spe)[rownames(spe) != "DNA1"],]
# spe

# READ IN SEGMENTATION MASKS ---------------------------------------------------

img_to_read <- c(64, 67, 71, 82, 84)
img_to_read <- paste0("(?i)", img_to_read, "(Prim|Rec)_00[1-7]{1}[_a-z]+\\.tiff$")


# Read in 16-bit unsigned integer segmentation masks
masks <- loadImages(io$inputs$segmentation_masks, 
                    pattern = img_to_read,
                    as.is = TRUE)

purrr::walk(masks@listData, ~print(.x@.Data[1:20,1:5]))


length(masks)

names(masks)
# 
# # Check if names of images and masks are the same
# all.equal(names(images), names(masks))
# 
# patient_id <- stringr::str_extract_all(names(images), 
#                                        "(?i)^\\d+(prim|rec)", 
#                                        simplify = TRUE)
# 
# indication <- ifelse(grepl("(?i)prim$", patient_id), "primary", "recurrent")
# 
# mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
#                                            patient_id = patient_id,
#                                            indication = indication)

# SAVE OBJECTS -----------------------------------------------------------------

saveRDS(spe, file.path(io$output$raw_data_out, "spe.rds")) 
# saveRDS(images, file.path(io$output$raw_data_out, "images.rds"))
# saveRDS(masks, file.path(io$output$raw_data_out, "masks.rds"))

# END ----
