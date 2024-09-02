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

    if(missing(dir_path)) stop("No directory path provided")
    if(missing(file_pattern)) stop("No file pattern provided")
        
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

io$inputs$data <-  find_file(io$inputs$input_dir, file_pattern = "spe_downstream")
io$inputs$prevelance_out <-  find_file(io$inputs$prevelance_out, file_pattern = "prevelance_data")

# LOAD DATA --------------------------------------------------------------------
spe <- readRDS(io$inputs$data)
regions <-  readRDS(io$inputs$prevelance_out)

# ADD LABELS AND COLOURS -------------------------------------------------------
io$plots$labels <- list()

# CLEAN DATA -------------------------------------------------------------------
# SAVE DATA --------------------------------------------------------------------
# END --------------------------------------------------------------------------
