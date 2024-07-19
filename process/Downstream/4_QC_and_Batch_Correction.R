# PACKAGES ---------------------------------------------------------------------
library(cytomapper)
library(tidyverse)
library(EBImage)
library(patchwork)
library(RColorBrewer)
library(scuttle)
library(dittoSeq)
library(ggrepel)

# I/O --------------------------------------------------------------------------
io <- list(
  
  inputs = list(
    functions = "process/Downstream/functions",
    raw_data = "data/downstream/raw",
    comp_data = "data/downstream/compensated",
    dim_reductions = "data/downstream/dim_reductions"
    ),
  
  output= list(
    data_out_dir = "outputs/QC",
    segmentation = "outputs/QC/segmentation",
    image_cell = "outputs/QC/image_cell",
    batch = "outputs/QC/batch"
    )
)

if(!dir.exists(io$output$data_out)) dir.create(io$output$data_out)

# Create QC out sub directories 
if(!dir.exists(io$output$segmentation)) dir.create(io$output$segmentation)
if(!dir.exists(io$output$image_cell)) dir.create(io$output$image_cell)
if(!dir.exists(io$output$batch)) dir.create(io$output$batch)

# READ IN IMAGES AND MASKS  ----------------------------------------------------
images <- readRDS(file.path(io$inputs$comp_data, "images_comp.rds"))
masks <- readRDS(file.path(io$inputs$raw_data, "masks.rds"))

names(masks) <- str_remove_all(names(masks), "_combined_cp_masks")

mcols(masks)$sample_id <- str_remove_all(mcols(masks)$sample_id, "_fullstack")

split_CytoImageList <- function(img_list){
  
  if(!(class(img_list) %in% "CytoImageList")){
    
    stop("image_list is not of class 'CytoImageList'",call. = F)
    
  }
  
  image_list_split <- list()
  
  for (image in seq_along(img_list)) {
    
    image_list_split[[image]] <- cytomapper::getImages(img_list, image)
    
  }
  
  names(image_list_split) <- names(img_list)
  
  return(image_list_split)
  
}

# SEGMENTATION MASK OVERLAYS VISUALISATION -------------------------------------

images <- split_CytoImageList(images)

plot_segement_masks <- function(img, mask, out_dir){

  sample_name <- names(img)
    
  print(glue::glue("Processing {sample_name}..."))
  
  norm_image <- normalize(img, separateImages = TRUE)
  norm_image <- normalize(norm_image, inputRange = c(0, 0.2))
  
  outfilename <- file.path(out_dir, paste0(sample_name, ".tiff")) 
  
  print(glue::glue("\tPlotting {sample_name}..."))
  
  plotPixels(norm_image,
             mask = masks[sample_name],
             missing_colour = "red",
             img_id = "sample_id", 
             colour_by = "DNA2", 
             colour = list(DNA2 = c("black", "white")), 
             margin = 3,
             image_title = list(cex = 3, 
                                colour = "yellow",
                                margin = c(10, 10)),
             scale_bar = NULL, 
             legend = NULL,
             return_images = FALSE,
             return_plot = FALSE, 
             display = "single",
             save_plot = list(filename = outfilename,
                              scale = 5)
  )
  
  print(glue::glue("\tDone!"))
  
  
}

all(names(images) == names(masks))

lapply(images, function(x) plot_segement_masks(img = x, masks, out_dir = io$output$segmentation))

# IMAGE LEVEL QC METRICS -------------------------------------------------------

rm(images)

images <- readRDS(file.path(io$inputs$comp_data, "images_comp.rds"))

qc_plots <- list()

# Image-level signal to noise ratios
# This will not change with segmentation as it's based soley on the images themselves
# there it only needs to be calculated once
# 
# cur_snr <- lapply(images, function(img){
#   mat <- apply(img, 3, function(ch){
#     # Otsu threshold
#     thres <- otsu(ch, range = c(min(ch), max(ch)))
#     # Signal-to-noise ratio
#     snr <- mean(ch[ch > thres]) / mean(ch[ch <= thres])
#     # Signal intensity
#     ps <- mean(ch[ch > thres])
#     
#     return(c(snr = snr, ps = ps))
#   })
#   t(mat) %>% as.data.frame() %>% 
#     mutate(marker = colnames(mat)) %>% 
#     pivot_longer(cols = c(snr, ps))
# })
# 
# cur_snr <- do.call(rbind, cur_snr)
# 
# image_snr_plot = cur_snr %>% 
#   group_by(marker, name) %>%
#   summarize(mean = mean(value),
#             ci = qnorm(0.975)*sd(value)/sqrt(n()),.groups = "keep") %>%
#   pivot_wider(names_from = name, values_from = c(mean, ci)) %>%
#   left_join(., as.data.frame(rowData(spe)), by =  dplyr::join_by(marker == name)) %>%
#   ggplot(aes(color = cell_type)) +
#   geom_point(aes(log2(mean_ps), log2(mean_snr), size = 5), show.legend = FALSE) +
#   geom_label_repel(aes(log2(mean_ps), log2(mean_snr), 
#                        label = marker), max.overlaps = 20) +
#   viridis::scale_color_viridis(option = "H",discrete = T) +
#   theme_classic(base_size = 15) +
#   theme(legend.title = element_blank()) +
#   ylab("Signal-to-noise ratio [log2]") +
#   xlab("Signal intensity [log2]")
# 
# 
# ggsave(file.path(io$output$image_cell,"Image_pixel_singal_to_noise_ratio.svg"),
#        plot = image_snr_plot,
#        device = "svg",dpi = 300,
#        width = 15, height = 10)

# SINGLE CELL QC METRICS -------------------------------------------------------

spe <- readRDS(file.path(io$inputs$comp_data, "spe_comp.rds"))

# Defining plotting colour schemes
spe@metadata$color_vectors$ROI <- 
  setNames(brewer.pal(length(unique(spe$ROI)), name = "BrBG"), unique(spe$ROI))

spe@metadata$color_vectors$patient <- 
  setNames(brewer.pal(length(unique(spe$patient)), name = "Set1"), 
           unique(spe$patient))

spe@metadata$color_vectors$surgery <- 
  setNames(viridis::viridis(length(unique(spe$surgery)), option = "H"), 
                    unique(spe$surgery))

spe@metadata$color_vectors$region_type <- 
  setNames(brewer.pal(length(unique(spe$region_type)), name = "BrBG"), 
         unique(spe$region_type))


# Total cells per region
qc_plots$total_cells = colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(cells = n()) %>%
  mutate(patient = str_extract(sample_id, "^\\d{2}")) %>%
  arrange(cells) 

qc_plots$total_cells = ggplot(qc_plots$total_cells, aes(group = patient, fill = patient)) +
  geom_col(aes(sample_id, cells)) + 
  theme_classic(base_size = 15) +
  coord_flip(ylim = c(0, 10000)) +
  viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text( hjust = 1, size = 8),
        legend.title = element_blank()) +
  ylab("total cells") + xlab("")

# Total area covered
qc_plots$cells_area_covered = colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(cell_area = sum(area),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(covered_area = cell_area / no_pixels) %>%
  mutate(patient = str_extract(sample_id, "^\\d{2}")) %>%
  arrange(cell_area) 


qc_plots$cells_area_covered = ggplot(qc_plots$cells_area_covered, aes(group = patient, fill = patient)) +
  geom_col(aes(sample_id, covered_area)) + 
  theme_classic(base_size = 15) +
  coord_flip(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::percent) +
  viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text( hjust = 1, size = 8),
        legend.title = element_blank()) +
  ylab("% area covered(cells)") + xlab("")


qc_plots$cell_areas <- colData(spe) %>%
  as.data.frame() %>%
  select(sample_id, area) %>%
  mutate(diameter = sqrt(area/pi)*2) %>%
  group_by(sample_id) %>%
  # summarise(min(diameter),
  #           mean(diameter),
  #           median(diameter),
  #           max(diameter),
  #           sd(diameter))
  ggplot() +
  geom_boxplot(aes(sample_id, diameter)) +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell Diameter") + xlab("")

ggsave(file.path(io$output$image_cell, "cell_QC.pdf"),
       plot = gridExtra::marrangeGrob(qc_plots, ncol = 1, nrow = 1, top = ""),
       device = "pdf", width = 10, height = 6, dpi = 100)


qc_plots$exprs <- dittoSeq::multi_dittoPlot(
  spe, 
  vars = rowData(spe)$name[as.logical(rowData(spe)$keep)],
  group.by = c("patient"), 
  plots = c("ridgeplot"), 
  assay = "exprs", 
  color.panel = metadata(spe)$color_vectors$patient)

ggsave(file.path(io$output$image_cell, "Staining_patterns.pdf"),
       plot = qc_plots$exprs,
       device = "pdf", width = 15, height = 40, dpi = 100)


image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$sample_id, 
                                   statistics= "mean",
                                   use.assay.type = "exprs")


qc_plots$mean_marker_exprs <- dittoHeatmap(object = image_mean, 
                                           genes = rownames(spe),
                                           assay = "exprs", cluster_cols = FALSE, scale = "none",
                                           heatmap.colors = viridis::viridis(100), 
                                           annot.by = c("surgery", "patient", "region_type"),
                                           order.by =  c("surgery", "patient", "region_type"),
                                           annotation_colors = list(surgery = metadata(spe)$color_vectors$surgery,
                                                                    patient = metadata(spe)$color_vectors$patient,
                                                                    region_type = metadata(spe)$color_vectors$region_type),
                                           show_colnames = TRUE)


ggsave(file.path(io$output$image_cell, "Mean_Expression.pdf"),
       plot = qc_plots$mean_marker_exprs,
       device = "pdf", width = 12, height = 8, dpi = 100)

# DIMENSIONALITY REDUCTION -----------------------------------------------------

optimise_neighbors <- function(input_mat, n_neighbors){
  
  set.seed(123)
  
  foo <- uwot::umap(input_mat, 
                    n_neighbors = n_neighbors,
                    min_dist = 0.01,
                    metric = "cosine",
                    n_epochs = 200,
                    learning_rate = 1,
                    verbose = TRUE)
}

dimplot_vars <- set_names(c("patient","surgery", "ROI", "region_type"),
                          c("patient","surgery", "ROI", "region_type"))

dimplot_vars <- as.list(dimplot_vars)

plot_batch_dims <- function(var, reduct = "UMAP"){
  
  if(grepl("_", var)){
    
    plot_title <- gsub("_"," ", dimplot_vars[1]) %>% 
      gsub("id", "ID", .) %>%
      tools::toTitleCase()
    
  } else plot_title <- tools::toTitleCase(var)
  
  if(!reduct %in% SingleCellExperiment::reducedDimNames(spe)){
    
    stop("Reduction not found!", call. = FALSE)
    
  }
  
  reduced_dim_plot <- dittoDimPlot(spe, 
                                   var = var, 
                                   reduction.use = reduct, 
                                   size = 0.2,
                                   color.panel = metadata(spe)[["color_vectors"]][[var]],
                                   legend.title = plot_title) + 
    ggtitle("Before Batch Correction") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", hjust = 0.5)
    )
  
  
  return(reduced_dim_plot)
  
}

batch_dimplots <- list()

# UMAP -------------------------------------------------------------------------

reducedDim(spe, "UMAP") <- optimise_neighbors(
  input_mat = t(assay(spe, "exprs")[rownames(spe)[rowData(spe)$downstream == 1],]),
  n_neighbors = 40)

batch_dimplots$UMAP<- lapply(dimplot_vars, plot_batch_dims, reduct = "UMAP")

# PACMAP -----------------------------------------------------------------------
dim_reduc_mat = t(assay(spe, "exprs")[rownames(spe)[rowData(spe)$downstream == 1],])

# calculate the number of neighbors for PacMAP
# round(10 + 15 * (log10(dim(dim_reduc_mat)[1]) - 4))     # 27 for theses data

write.csv(dim_reduc_mat, 
          file.path(io$inputs$dim_reductions,"dim_reduction_matrix.csv"), row.names = F)

pacMAP_embeddings = read.csv("data/downstream/dim_reductions/PacMAP_embeddings.csv", 
                             header = F)

pacMAP_embeddings = as.matrix(pacMAP_embeddings)
rownames(pacMAP_embeddings) = colnames(spe)
colnames(pacMAP_embeddings) = NULL

reducedDim(spe, "pacMap") <- pacMAP_embeddings

batch_dimplots$pacMAP<- lapply(dimplot_vars, plot_batch_dims, reduct = "pacMap")

# Z-SCORE NORMALISATION --------------------------------------------------------

mat = t(assay(spe, "exprs")[rownames(spe)[rowData(spe)$downstream == 1],])

reducedDim(spe, "arcsinh_Zscore") <- scale(mat, center = T, scale = T)

set.seed(1234)

UMAP_arcsinh_Zscore <- optimise_neighbors(input_mat = reducedDim(spe, "arcsinh_Zscore"), 
                                          n_neighbors = 40)

rownames(UMAP_arcsinh_Zscore) <- colnames(spe)

reducedDim(spe, "UMAP_arcsinh_Zscore") <- UMAP_arcsinh_Zscore


batch_dimplots$zscore <- lapply(dimplot_vars, plot_batch_dims, 
                                reduct = "UMAP_arcsinh_Zscore")

# HARMONY CORRECTION -----------------------------------------------------------
library(harmony)

mat <- t(assay(spe, "exprs")[rownames(spe)[rowData(spe)$downstream == 1],])

harmony_emb <- HarmonyMatrix(mat, 
                             meta_data = colData(spe),
                             vars_use =  "patient",
                             do_pca = TRUE)

reducedDim(spe, "harmony") <- harmony_emb

# Saving the harmony embeddings for PacMAP
write.csv(harmony_emb, file.path(io$inputs$dim_reductions,
                                 "harmony_embeddings_matrix.csv"), 
          row.names = F)

pacMAP_harmony_embeddings = 
  read.csv("data/downstream/dim_reductions/PacMAP_harmony_embeddings.csv", header = F)

pacMAP_harmony_embeddings = as.matrix(pacMAP_harmony_embeddings)
rownames(pacMAP_harmony_embeddings) = colnames(spe)
colnames(pacMAP_harmony_embeddings) = NULL

reducedDim(spe, "pacMap_harmony") <- pacMAP_harmony_embeddings

batch_dimplots$pacMAP_harmony <- lapply(dimplot_vars, plot_batch_dims, reduct = "pacMap_harmony")


set.seed(1234)
UMAP_harmony <- optimise_neighbors(input_mat = reducedDim(spe, "harmony"), 
                                   n_neighbors = 40)

rownames(UMAP_harmony) <- colnames(spe)

reducedDim(spe, "UMAP_harmony") <- UMAP_harmony


batch_dimplots$harmony <- lapply(dimplot_vars, plot_batch_dims, reduct = "UMAP_harmony")

# BATCH CORRECTION VISUALISATION -----------------------------------------------

batch_dimplots$zscore = purrr::map(batch_dimplots$zscore, ~ .x +   ggtitle("After z-score correction"))
batch_dimplots$harmony = purrr::map(batch_dimplots$harmony, ~ .x +   ggtitle("After Harmony correction"))
batch_dimplots$pacMAP_harmony = purrr::map(batch_dimplots$pacMAP_harmony, ~ .x +   ggtitle("After Harmony correction"))


zscore_plots <- purrr::map2(batch_dimplots$UMAP, batch_dimplots$zscore, ~{
  
  combined <- .x + .y + patchwork::plot_layout(guides = "collect") 
  
})


harmony_plots <- map2(batch_dimplots$UMAP, batch_dimplots$harmony, ~{
  
  combined <- .x + .y + patchwork::plot_layout(guides = "collect") 
  
})

harmony_pacmap_plots <- map2(batch_dimplots$pacMAP, batch_dimplots$pacMAP_harmony, ~{
  
  combined <- .x + .y + patchwork::plot_layout(guides = "collect") 
  
})

pdf(file = file.path(io$output$batch, "UMAP_zscore_correction.pdf"), 
    onefile = T, width = 14, height = 7)

zscore_plots

dev.off()


pdf(file = file.path(io$output$batch, "UMAP_Harmony_Correction.pdf"), 
    onefile = T, width = 14, height = 7)

harmony_plots

dev.off()


pdf(file = file.path(io$output$batch, "PacMAP_Harmony_Correction.pdf"), 
    onefile = T, width = 14, height = 7)

harmony_pacmap_plots

dev.off()

# SAVE PROCESSED OBJECTS -------------------------------------------------------

saveRDS(spe, file.path(io$inputs$comp_data,"spe_comp.rds"))

# END --------------------------------------------------------------------------


