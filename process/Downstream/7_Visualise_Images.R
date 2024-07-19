# PACKAGES ----

library(SpatialExperiment)
library(cytomapper)

# LOAD DATA ----

spe <- readRDS("data/Downstream_Analysis/spe_comp_subclustered.rds")
images <- readRDS("data/Downstream_Analysis/images_comp.rds")
masks <- readRDS("data/Downstream_Analysis/masks.rds")

# SAMPLE IMAGES 3 IMAGES ----

unique(spe$sample_id)

cur_id <- c("82Rec_001")
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

# VISUALISING NORMAL BRAIN STUCTURES ----

# SMA:vasculature, GFAP:astrocyte, MOG:oligodendrocyte

cur_id <- c("67Prim_003")
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]


# normalise the channel intensities
norm_images <- normalize(cur_images)

# clip channel intensities
norm_images <- normalize(norm_images, inputRange = c(0, 0.15))


tiff("Outputs/normal_brain_markers.tiff", units="in", width=10, height=10, res=600)

plotPixels(norm_images,
           return_plot = FALSE,
           
           return_images = FALSE,
           
           scale_bar = list(length = 100,
                            label = "",
                            colour = "white",
                            cex = 1),
           legend = NULL, 
           
           image_title = NULL,
           
           display = "single",
           
           colour_by = c("SMA","GFAP", "MOG"),
          
           colour = list(SMA = c("black", "red"),
                         GFAP = c("black", "magenta4"),
                         MOG = c("black", "chartreuse")
                         ),
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

tiff("Outputs/neoplastic_markers.tiff", units="in", width=10, height=10, res=600)

plotPixels(norm_images,
           return_plot = FALSE,
           
           return_images = FALSE,
           
           scale_bar = list(length = 100,
                            label = "",
                            colour = "white",
                            cex = 1),
           legend = NULL, 
           
           image_title = NULL,
           
           display = "single",
           
           colour_by = c("SOD2","DLL3", "OLIG1"),
           
           colour = list(SOD2 = c("black","green"),
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


tiff("Outputs/various_markers.tiff", units="in", width=10, height=10, res=600)

plotPixels(norm_images,
           return_plot = FALSE,
           
           return_images = FALSE,
           
           scale_bar = list(length = 100,
                            label = "",
                            colour = "white",
                            cex = 1),
           legend = NULL, 
           
           image_title = NULL,
           
           display = "single",
           
           colour_by = c("Ki67","TNC", "IBA1"),
           
           colour = list(Ki67 = c("black","red"),
                         TNC = c("black", "green"),
                         IBA1 = c("black", "blue")
           ),
)

dev.off()