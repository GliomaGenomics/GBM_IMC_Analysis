# This script carries out spillover correction on the metal spots and also 
# compensates the images by accounting for the above correction. This 
# Single cell experiment object that is loaded will have raw counts, denoted by the 
# assay label "counts" and asinh () transformed counts, denoted by the
# assay label "exprs".

# PACKAGES ---------------------------------------------------------------------
library(imcRtools)
library(CATALYST)
library(cytomapper)
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(dittoSeq)
library(patchwork)
library(BiocParallel)

# I/O --------------------------------------------------------------------------
io <- list(
  
  inputs = list(
    functions = "process/Downstream/functions",
    data = "data/downstream/raw",
    spillover_path = file.path("data/spillover_correction/spillover_matrix.csv"),
    compensation_files = file.path("data/spillover_correction/Compensation_files"),
    fullstack_images = file.path("data/fullstacks"),
    segmentation_masks = file.path("data/masks/tiff")),
  
  output= list(
    data_out = "data/downstream/compensated",
    plots_out = "outputs/QC/spillover_correction")
)

if(!dir.exists(io$output$data_out)) dir.create(io$output$data_out)
if(!dir.exists(io$output$plots_out)) dir.create(io$output$plots_out)

source(file.path(io$inputs$functions, "utils.R"))

# The spillover matrix was previously created for the steinbock process and so does 
# not need to be regenerated. Instead, the values can be used to correct the single
# cell expression data.

# A Hyperbolic arcsine (arcsinh) scaling factor (5) was applied to the data with 
# the intention of ensuring the high variance of a potentially 
# larger value (positive signal) is given the same significance and 
# weighting as a negative (low) signal that has low variance. 
# This scale normalisation is essential for dimensionality reduction and 
# clustering algorithms to work well.

# SINGLE-CELL COMPENSATION -----------------------------------------------------
spe = readRDS(file.path(io$inputs$data, "spe.rds"))

# From the previous cofactor tuning a asinh cofactor of 1 was determined for 
# transforming the raw counts
sm <- as.matrix(read.csv(io$inputs$spillover_path, row.names = 1))

rowData(spe)$channel_name <- paste0(rowData(spe)$channel, "Di")

# Compensate the raw counts
spe <-  CATALYST::compCytof(x = spe, 
                            sm = sm,
                            assay = "counts",
                            transform = TRUE,
                            cofactor = 1, 
                            isotope_list = CATALYST::isotope_list, 
                            overwrite = FALSE)

assay(spe, "exprs") <- asinh(counts(spe)/1)

assayNames(spe)

# PLOTTING COMPENSATION EXPRESSION ---------------------------------------------
panel <- read.csv("data/panel/panel.csv")

# Plot adjacent channels before and after compensation
panel <- panel[panel$keep == 1,]
panel <- panel[order(panel$channel_number),]

x =  panel$name[1:nrow(panel)-1]
y = panel$name[seq(2,nrow(panel))]


adjacent_channels = function(marker_x, marker_y, before_assay, after_assay){
  
  plot_theme <- theme_classic() + 
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  
  before <- dittoScatterPlot(spe, x.var = marker_x, y.var = marker_y,
                             assay.x = before_assay, assay.y = before_assay) +
    ggtitle(glue::glue("Before compensation ({before_assay})")) +
    plot_theme
  
  after <- dittoScatterPlot(spe, x.var = marker_x, y.var = marker_y,
                            assay.x = after_assay, assay.y = after_assay) +
    ggtitle(glue::glue("After compensation ({after_assay})")) +
    plot_theme
  
  return(before + after)
  
  
}

plot_adjacent_channels <- function(adj_plots, out_filename, out_dir){
  
  purrr::walk(adj_plots, ~{
    
    ggsave(filename = tempfile(tmpdir = out_dir,
                               fileext = ".pdf",
                               pattern = "adjcomp_"),
           plot = .x,
           width = 17, height = 8, dpi = 300, 
           limitsize = T, 
           device = "pdf")
    
  })
  
  pdftools::pdf_combine(input = list.files(out_dir,pattern = "^adjcomp_", full.names = T),
                        output = file.path(out_dir, out_filename))
  
  Sys.sleep(5)
  
  fs::file_delete(
    list.files(out_dir,pattern = "^adjcomp_",full.names = T)
  )
  
}

plots = list()

plots$counts <- purrr::map2(x, y, 
                            ~ adjacent_channels(marker_x = .x,
                                                marker_y = .y,
                                                before_assay = "counts",
                                                after_assay = "compcounts"))

plot_adjacent_channels(plots$counts,
                       out_dir = io$output$plots_out, 
                       out_filename = "counts_biscatter.pdf")


plots$exprs <- purrr::map2(x, y, 
                           ~adjacent_channels(marker_x = .x,
                                              marker_y = .y,
                                              before_assay = "exprs",
                                              after_assay = "compexprs"))

plot_adjacent_channels(plots$exprs,
                       out_dir = io$output$plots_out, 
                       out_filename = "asinh_exprs_biscatter.pdf")


# SAVE COMPENSATED SINGLE CELL EXPRESSION --------------------------------------

assay(spe, "counts") <- assay(spe, "compcounts") 
assay(spe, "exprs") <- assay(spe, "compexprs") 
assay(spe, "compcounts") <- assay(spe, "compexprs") <- NULL

spe

saveRDS(spe, .new_file("spe_comp",path = io$output$data_out, ext = "rds"))

# The images do not need to be compensated as this step has already been done
# previously and the images haven't changed, but just the single cells and the 
# segmentation masks.

# END --------------------------------------------------------------------------