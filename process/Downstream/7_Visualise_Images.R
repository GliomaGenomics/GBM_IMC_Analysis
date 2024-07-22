options(scipen = 999)

# PACKAGES ---------------------------------------------------------------------
library(SpatialExperiment)
library(cytomapper)
library(IMCfuncs)

# IO ---------------------------------------------------------------------------
io = list(
    inputs = list(
        spe = "data/downstream/compensated",
        images = "data/downstream/compensated/images_comp.rds",
        masks = "data/downstream/raw/masks.rds"
    ),
    
    output= list(
        
    ),
    
    plots = list(
        args = list()
    )
)

io$inputs$spe = list.files(io$inputs$spe, "^spe_comp_\\d+", full.names = TRUE)
io$inputs$spe = io$inputs$spe[length(io$inputs$spe)]

# Create the output directories (if they do not exist)
source("process/misc/Output_directory_trees.R")

io$output$image_visualisation = dir_tree(
    root_directory = roots$image_visualisation,
    directory_tree = directory_trees$image_visualisation
)

rm(roots, directory_trees)

# LOAD DATA --------------------------------------------------------------------
spe = readRDS(io$inputs$spe)
images = readRDS(io$inputs$images)
masks = readRDS(io$inputs$masks)

# VISUALISING VASCULAR/OLIGODENDROCYTES ----------------------------------------

cur_id = c(
    # "64Prim_001"
    "64Rec_001"
    )

# cur_id = grep("^82", unique(spe$sample_id),value = TRUE)

cur_images = images[names(images) %in% cur_id]
cur_masks = masks[names(masks) %in% cur_id]

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
    width= 30,
    height= 15,
    res=600
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
    
    image_title = list(text = mcols(cur_images)$`patient_id.V1`,
                       colour = "white",
                       margin = c(10,10),
                       position = "topleft",
                       font = 2,
                       cex = 5),

    scale_bar = list(length = 100,
                     label = expression("100 " ~ mu * "m"),
                     colour = "white",
                     lwidth = 10,
                     position = "bottomright",
                     margin = c(75,20),
                     frame = length(cur_id),
                     cex = 5),
    margin = 40
    
    )

dev.off()


# VISUALISING ASTROCYTES/NEURONS -----------------------------------------------

cur_id = c(
    # "82Prim_002"
    "82Rec_002"
           )

cur_images = images[names(images) %in% cur_id]
cur_masks = masks[names(masks) %in% cur_id]

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
    width= 30,
    height= 15,
    res=600
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
    
    image_title = list(text = mcols(cur_images)$`patient_id.V1`,
                       colour = "white",
                       margin = c(10,10),
                       position = "topleft",
                       font = 2,
                       cex = 5),
    
    scale_bar = list(length = 100,
                     label = expression("100 " ~ mu * "m"),
                     colour = "white",
                     lwidth = 10,
                     position = "bottomright",
                     margin = c(75,20),
                     frame = length(cur_id),
                     cex = 5),
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

# END --------------------------------------------------------------------------