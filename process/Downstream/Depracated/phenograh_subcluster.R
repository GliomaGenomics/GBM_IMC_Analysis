initial_undefined = spe[markers,spe$initial_anno == "undefined"]

# CLUSTERING BATCH CORRECTED EMBEDDINGS ----

set.seed(123)
phenograph_clusters = Rphenograph(
    data = reducedDim(initial_undefined, "harmony_minmax"), k = 15
    )

# Add meta clusters back to the input matrix
initial_undefined$undefined_anno_1 = phenograph_clusters[[2]]$membership

# CREATE HEATMAPS ----
# Create the Heatmap Data
heatmap_data = anno_expression(object = initial_undefined, 
                               assay = "exprs", 
                               anno_label = "undefined_anno_1",
                               subset_rows = markers,
                               anno_summary_stat = "median",
                               heatmap_scale_fun = "minmax",
                               heatmap_scale_across = "features")


set_heatmap_global_opts(
    row_title_fill = heatmap_data$plot_opts$row_title_bkg_colours,
    row_title_text = heatmap_data$plot_opts$row_title_text_colours)

heatmap_data$ht = Heatmap(
    matrix = heatmap_data$heatmap_data, 
    col = heatmap_data$plot_opts$cols,
    rect_gp = gpar(col = "grey75", lwd = 0.70),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    heatmap_legend_param = list(
        title = heatmap_data$plot_opts$scale_title,
        at = heatmap_data$plot_opts$scale_breaks,
        legend_height = unit(25, "cm"),
        grid_width = unit(12.5, "mm")
    ),
    
    top_annotation = heatmap_data$column_anno,
    
    row_split = heatmap_data$plot_opts$rowsplits,
    
    row_names_side = "left",
    row_title_side = "left",
    row_names_max_width = unit(500, "mm"),
    row_gap = unit(5, "mm"),
    column_names_centered = TRUE
)


pdf(.new_file("Initial_undefined_anno", io$output$som_pheno),
    width = 40, height = 27.5, onefile = T)

draw(heatmap_data$ht, 
     padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "right")

dev.off()

# MEM SCORING ----

heatmap_data$MEM_values <- MEM(
    exp_data = cbind(t(heatmap_data$heatmap_data), 
                     cluster = as.numeric(colnames(heatmap_data$heatmap_data))),
    transform = FALSE, 
    choose.markers = FALSE, # Change choose.markers to TRUE to see and select channels in console
    markers = "all",
    choose.ref = FALSE,
    zero.ref = FALSE,
    rename.markers = FALSE, # Change rename.markers to TRUE to see and choose new names for channels in console
    new.marker.names = "none",
    IQR.thresh = NULL)

heatmap_data$MEM_scores = gather_MEM_scores(
    MEM_score_matrix = heatmap_data$MEM_values$MEM_matrix[[1]],
    threshold = 3)

write.csv(x = heatmap_data$MEM_scores, row.names = FALSE,
          file = .new_file("Initial_undefined_MEM",
                           io$output$som_pheno, ext = "csv")
          )


# MAP ANNOTATED CELLS ----

cluster_labels = sort(list.files(io$output$som_pheno, pattern = "MEM", 
                                 full.names = TRUE), decreasing = TRUE)

cluster_labels = read.csv(file = cluster_labels[1], 
                          colClasses =  rep("character", 4))

cluster_labels = cluster_labels[c("cluster", "label")]

cluster_labels$label[cluster_labels$label == ""] = "undefined"

initial_undefined$undefined_anno = factor(cluster_labels$label[match(initial_undefined$undefined_anno_1, 
                                                                     cluster_labels$cluster)])

initial_undefined$undefined_anno = as.character(initial_undefined$undefined_anno)

spe$initial_anno = as.character(spe$initial_anno)

colData(spe)[rownames(initial_undefined@colData) , "initial_anno"] = initial_undefined$undefined_anno

table(spe$initial_anno)


# END ----