# For the first bit of exploratory analysis we can check the cell type prevalence 
# across various stratifications, starting with the the specific tissue samples. 
# In order to do this we will first obtain counts for each of the cells and then 
# subset/aggregate them to show the cell prevalence across different groupings.

# Previously, i clustered and annotated the imaging mass cytometry single cells. 
# During this process a number of cells remained unannotated as their 
# identity could not be easily inferred using one of more of the cell type markers 
# which were used. Also, in some cases the abundance of the metal isotopes corresponding 
# to the cell type markers was indiscriminate in distinguishing one particular cell 
# type or cell state and so was left "unknown". 
# This initial cell type prevalence analysis will include these "unknown" markers.


# PACKAGES ---------------------------------------------------------------------
library(SpatialExperiment)
library(dittoSeq)
library(imcRtools)
library(ggplot2)
library(viridis)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(IMCfuncs)

# I/O --------------------------------------------------------------------------
io = list(
    inputs = list(
        functions = "process/Downstream/functions",
        comp_data = "data/downstream/compensated",
        dim_reductions = "data/downstream/dim_reductions"
    ),
    output= list(
        data_out_dir = "outputs/spatial_analysis",
        cell_pheno = "outputs/cell_phenotyping",
        cell_abundances = "outputs/spatial_analysis/labelled_cell_abundances",
        cell_prevalences = "outputs/spatial_analysis/prevalences",
        cell_prevalence_comps = "outputs/spatial_analysis/prevalences/comp_stats",
        cell_interactions ="outputs/spatial_analysis/interactions"
    ),
    plots = list(
        
    )
)

# Create the directories if they do not exist
ndirs(io$output)


# LOAD LABELLED DATA  ----------------------------------------------------------
spe = list.files(
  io$inputs$comp_data,
  pattern = "spe_comp_\\d{4}",
  full.names = T
)
spe = readRDS(sort(spe,decreasing = T)[1])

# CREATE CELL TYPES PREVELENCE DATA  -------------------------------------------

# Get the region level metadata
regions = as.data.frame(colData(spe)) %>%
  distinct(
    image = sample_id,
    patient,
    surgery,
    ROI,
    responder_type,
    hist_region_label =  region_type,
    imc_region_label = region_type_new
  ) %>%
  # The proceeding line excludes the additional 67Prim images
  filter(!image %in% paste0("67Prim_00", 4:7)) 


regions = as.data.frame(spe@colData) %>%
  group_by(image = sample_id) %>%
  summarise(
    total_cells = n(),
    undefined = length(which(is.na(manual_gating))),
    labelled = total_cells - undefined,
    region_immune = length(which(main_anno == "Immune")),
    region_hypoxia = sum(high_hypoxia),
    region_prolif = sum(high_prolif),
    .groups = "keep"
  ) %>% 
  ungroup() %>%
  left_join(regions, ., by = "image")


# Add anno labels to region level metadata
regions$anno_labels = as.data.frame(colData(spe)) %>%
  filter(sample_id %in% regions$image) %>%
  mutate(across(c("main_anno", "manual_gating"), as.character)) %>%
  select(image = sample_id, cell_main = main_anno, cell_fine = manual_gating) %>%
  tidyr::drop_na() %>%
  split(.$image) %>%
  lapply(\(x) setNames(x$cell_fine, x$cell_main))

count_labels = function(annotation, labels){
  
  table(factor(annotation, levels = labels))
  
}

regions$main_anno = lapply(regions$anno_labels, \(x){
  
  count_labels(names(x),  tools::toTitleCase(names(spe@metadata$cell_types)))
  
})

regions$fine_anno = lapply(regions$anno_labels, \(x){
  
  # remove the hybrid anno labels 
  fine_anno_labels = unlist(spe@metadata$cell_types)[unlist(spe@metadata$cell_types) != "Hybrid"]
  
  count_labels(x, fine_anno_labels)
  
  })


# Set factors
regions$surgery = factor(regions$surgery, names(plot_colours$surgery))
regions$responder_type = factor(regions$responder_type, names(plot_colours$responder_type))
regions$patient = factor(regions$patient, names(plot_colours$patient))
regions$ROI = factor(regions$ROI, names(plot_colours$ROI))
regions$imc_region_label = factor(regions$imc_region_label, names(plot_colours$region_type))
regions$hist_region_label = factor(regions$hist_region_label, names(plot_colours$region_type))

rm(count_labels)

# PLOT LABELLED/UNDEFINED LABEL COUNTS -----------------------------------------

io$plots$undefined_cells = ggplot(regions, 
                                 aes(x = image, y = undefined, fill = patient)) +
  geom_col() +
  xlab("") +
  ylab("undefined cells") +
  theme_classic(base_size = 20) +
  scale_fill_manual(values = IMCfuncs::plot_colours$patient) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        legend.title = element_blank())


io$plots$labelled_cells = ggplot(regions, 
                                aes(x = image, y = labelled, fill = patient)) +
  geom_col() +
  xlab("") +
  ylab("labelled cells") +
  theme_classic(base_size = 20) +
  scale_fill_manual(values = IMCfuncs::plot_colours$patient) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        legend.title = element_blank())


# Save labelled and undefined count plots

# pdf(nf("cell_labelling_counts.pdf", io$output$cell_prevalences),
#     onefile = T, width = 15, height = 10)
# 
# io$plots$undefined_cells
# io$plots$labelled_cells
# 
# dev.off()


# RE-LABEL CELLS BASED ON CELL COUNTS ------------------------------------------
regions_labels = regions %>% 
  select(image, total_cells, starts_with("region_")) %>%
  tidyr::pivot_longer(
    cols = starts_with("region_"),
    names_to = "region",
    names_prefix = "region_",
    values_to = "cells",
  ) %>%
  mutate(percent = round(cells/total_cells * 100)) %>%
  mutate(across(region, ~factor(.x, levels = c("immune", "prolif", "hypoxia"))))


# Function to add new region labels based on the percentage of cell that are 
# labelled as either immune, hypoxia or proliferative.  
add_region_cols = function(regions, count_col, region_labels, 
                           immune_min = 20, 
                           state_min = 50){
  
  out = as.list(setNames(levels(regions[[region_labels]]),
                         levels(regions[[region_labels]]))
                )
  
  out = lapply(out, \(x){
    
    enrichment = ifelse(regions[[region_labels]] == x, regions[[count_col]], 0)
    if(grepl("(?i)immune", x)) enrichment >= immune_min else enrichment >= state_min
    
    
  })
  
  return(cbind(regions, dplyr::bind_cols(out)))
  
}

# Apply new region labels based on min defined thresholds
io$plots$args = list(immune_min = 20, state_min = 25)

regions_labels = add_region_cols(
  regions_labels,
  count_col = "percent",
  region_labels = "region",
  immune_min = io$plots$args$immune_min,
  state_min = io$plots$args$state_min
)

# 82Prim samples were a wee suspect because they had a large amount of 
# unlabelled cells, yet the number of labelled cells is not too dissimilar when 
# compared to the other samples which do also have lower labelled cells. 
# This is likely due to a very high cell density in the particular sample. 
# After looking at 82prim_002 this is going to be re-labelled as having a 
# high immune fraction

regions_labels$immune[regions_labels$image == "82Prim_002" & 
                      regions_labels$region == "immune"] = TRUE

# Add colour info for facets based on the region labels
regions_labels$facet_background = plot_colours$patient[str_extract(regions_labels$image, "^\\d{2}")]
regions_labels$facet_text = IMCfuncs::get_text_color(regions_labels$facet_background)

# PLOT REGION LABELS -----------------------------------------------------------

# Plot the region proportions along with the min thresholds used for each of the
# cut-off values.
io$plots$regions =
  ggplot(data = regions_labels, 
         aes(x = region, y = percent, fill = region, group = image)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_abline(aes(intercept = io$plots$args$immune_min, slope = 0, color = "Immune Threshold"),
              lty = 2) +
  geom_abline(aes(intercept = io$plots$args$state_min, slope = 0, color = "State Threshold"),
                lty = 2) +
  ggh4x::facet_wrap2(~image, ncol = 3, 
              strip =  ggh4x::strip_themed(
                
                background_x = lapply(split(regions_labels, regions_labels$image), \(x) {
                  element_rect(fill = unique(x[["facet_background"]]))
                  }),
                
                text_x = lapply(split(regions_labels, regions_labels$image), \(x) {
                  element_text(colour = unique(x[["facet_text"]]),
                               face = "bold")
                })
                )) +
  xlab("") +
  ylab("Percentage") +
  scale_fill_manual(values = plot_colours$region_type[levels(regions_labels$region)]) +
  scale_color_manual(values = c("Immune Threshold" = "lightpink2",
                                "State Threshold" = "skyblue1")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
        )


# Each image represents a tissue section that was selected based on 
# positive immunohistochemistry markers for either immune-rich, proliferative or 
# hypoxic regions. Each patient and surgery sample pairing was selected to have at least
# one of each region type. We can use the imaging mass cytometry labelled cell-type
# and state prevalences to corroborate if the histochemistry region labels were 
# indeed reflective of the cells initially labelled based on IHC single cell-type/state
# canonical markers. The region types and therefore the labels are not mutually 
# exclusive, e.g, any given region can be high (as determined by an appropriate threshold) 
# in immune cells and also comprise cells that are in a highly proliferate/hypoxic state.   

regions = regions_labels %>%
  group_by(image) %>%
  summarise(immune = any(immune),
            prolif = any(prolif),
            hypoxia = any(hypoxia), 
            .groups = "keep"
            ) %>%
  ungroup() %>%
  left_join(regions, ., by = "image")


region_present_data = regions %>%
  select(image, patient, surgery, ROI, immune, prolif, hypoxia)

# Add colour info for facets based on the region labels
region_present_data$facet_background = plot_colours$patient[region_present_data$patient]
region_present_data$facet_text = IMCfuncs::get_text_color(region_present_data$facet_background)


region_present_data = region_present_data %>% 
  tidyr::pivot_longer(c("immune", "prolif", "hypoxia"), 
                      names_to = "region", 
                      values_to = "present") %>%
  mutate(patient_surgery = paste0(patient, surgery))


io$plots$regions_present = region_present_data %>%
  ggplot(aes(x = region, y = ROI, fill = present)) +
  geom_point(shape = 21, size = 45, color = "grey75") +
  xlab("") +
  ylab("") +
  ggh4x::facet_wrap2(~ patient_surgery, ncol = 2,
    strip =  ggh4x::strip_themed(
      
      background_x = lapply(split(region_present_data, region_present_data$patient_surgery), \(x) {
                         element_rect(fill = unique(x[["facet_background"]]))
                       }),
      
      text_x = lapply(split(region_present_data, region_present_data$patient_surgery), \(x) {
                         element_text(colour = unique(x[["facet_text"]]),
                                      face = "bold")
                       })
    )
  ) +
  
  scale_fill_manual(values = c("FALSE" = "#000004FF","TRUE" = "#FCFDBFFF")) +
  theme_classic(base_size = 20) +
  theme(legend.title = element_blank())


# Save re-labelled region plots
# svglite::svglite(nf("regions_labels.svg",io$output$cell_prevalences),
#                  width = 12, height = 25)
# io$plots$regions
# dev.off()
# 
# 
# 
# svglite::svglite(nf("regions_per_sample.svg",io$output$cell_prevalences),
#                  width = 15, height = 30)
# io$plots$regions_present
# dev.off()

rm(region_present_data, regions_labels, add_region_cols)

# CALCULATE SHANNON ENTROPIES --------------------------------------------------
# Here we are seeking to measure heterogeneity by 
# using Shannon entropy (H) based on annotated cell subtypes as determined from
# the cell clustering/phenotyping process. To account for the different number 
# of cells per sample, we will sub-sample n (1000) cells from each patient sample (i), 
# over ten (three is default) rounds can then calculate the Shannon entropy of 
# each cell-type frequency occurrence.

calculate_props = function(x, labels, sample_size = 1000, rounds = 3){
  
  base = setNames(object = rep(0, length(labels)), nm= labels)

  out = vector("list", rounds)
  
  for (i in seq_len(rounds)){
    
    exp = prop.table(table(sample(x, sample_size)))
    exp = set_names(as.numeric(exp), names(exp))
    
    exp = c(base[setdiff(names(base), names(exp))], exp)
  
    out[[i]] = philentropy::H(exp, unit = "log2")
    
  }
  
  return(unlist(out, use.names = FALSE))
  
}

# Calculate the entropy for each cell type
set.seed(123)
regions$entropy = map(regions$anno_labels, ~calculate_props(
  .x,
  labels = levels(spe$manual_gating),
  sample_size = 1000,
  rounds = 10
))

rm(calculate_props)

# PLOT SHANNON ENTROPIES -------------------------------------------------------
# As we are mainly interested in how things change in response to treatment in each
# patient we will initially seek to make the following comparisons:

# convert the regions data.frame to a long format based on the region type
regions_long = regions %>% 
  tidyr::pivot_longer(
    immune:hypoxia,
    names_to = "region",
    values_to = "region_present"
  ) %>% 
  filter(region_present == TRUE)


# P vs R by patient
# up vs down by surgery
# P vs R by regions
io$plots$args$shannon_entropy_comps = tribble(
  
  ~df,              ~response_var,      ~comp_var,       ~group_var,         ~plot_names,
  "regions",          "entropy",         "surgery",       "patient",          "patient_surgery", 
  "regions",          "entropy",         "surgery",       "responder_type",   "responder_surgery",
  "regions_long",     "entropy",         "surgery",       "region",           "region_surgery"
  
)

io$plots$shannon_entropies = pmap(io$plots$args$shannon_entropy_comps, ~{
  
  df_to_use = get(..1, envir = globalenv())
  
  stats = compare_groups(
      df = df_to_use,
      response_var = ..2,
      nested_response_var = TRUE,
      comp_var = ..3,
      group_var = ..4,
      facetted_plot_dims = TRUE
      )
  
  out_boxplot = comp_boxplot(
      df = tidyr::unnest(df_to_use, ..2),
      x = ..3,
      y = ..2,
      facet_by = ..4,
      fill = ..3,
      color_palette = plot_colours$surgery
      )
  
  
  out = out_boxplot +
    # Add statistics
    ggpubr::stat_pvalue_manual(stats,  label = "p_signif", tip.length = 0, size = 7) +
    
    # Add 10% spaces between the p-value labels and the plot border
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
    
    # Add custom theme
    IMCfuncs::facetted_comp_bxp_theme(hide_x_labs = T) 
    
    # Expand the y axis scale to zoom into the box plots
    # # + ggplot2::coord_cartesian(ylim = c(0, stats$y.position[1])) 
  
  
  return(out)
  
})

names(io$plots$shannon_entropies) = io$plots$args$shannon_entropy_comps$plot_names


iwalk(io$plots$shannon_entropies, ~{
  
  # Save the Shannon entropy comparisons
  svglite::svglite(nf(paste0(.y, "_sh_entropy.svg"), io$output$cell_prevalences),
                   width = 15, height = 10)
  
  print(.x)
  
  dev.off()
  
  
})

# PLOT CELL PREVELENCES --------------------------------------------------------

# groups: patient, regions, responder types
io$plots$args$prevelence_comps = tribble(
  
  ~df,              ~x,                ~facet_by,                    ~fill,          ~color_palette,            ~plot_names,
  "regions",        "ROI",             "patient ~ surgery",          "main_anno",    plot_colours$main_anno,   "patient_surgery_main",
  "regions",        "ROI",             "patient ~ surgery",          "fine_anno",    plot_colours$cell_anno,   "patient_surgery_fine",
  "regions",        "responder_type",  "patient ~ surgery",          "main_anno",    plot_colours$main_anno,   "patient_responder_surgery_main",
  "regions",        "responder_type",  "patient ~ surgery",          "fine_anno",    plot_colours$cell_anno,   "patient_responder_surgery_fine",
  "regions_long",   "region",          "patient ~ surgery",          "main_anno",    plot_colours$main_anno,   "region_patient_surgery_main",
  "regions_long",   "region",          "patient ~ surgery",          "fine_anno",    plot_colours$cell_anno,   "region_patient_surgery_fine",
  "regions_long",   "region",          "responder_type ~ surgery",   "main_anno",    plot_colours$main_anno,   "region_responder_surgery_main",
  "regions_long",   "region",          "responder_type ~ surgery",   "fine_anno",    plot_colours$cell_anno,   "region_responder_surgery_fine",
)

# cell_props_facetted(
#   df = regions,
#   x = "responder_type",
#   facet_by = "patient ~ surgery",
#   fill = "main_anno",
#   unnest_labels = TRUE,
#   color_palette = plot_colours$main_anno
# ) +
#   facetted_cell_prop_theme(text_size = 20)

io$plots$cell_prevelences = pmap(io$plots$args$prevelence_comps, ~{
  
 out = cell_props_facetted(
   df = get(..1, envir = globalenv()),
   x = ..2,
   facet_by = ..3,
   fill = ..4,
   unnest_labels = TRUE,
   color_palette = ..5
 ) +
   facetted_cell_prop_theme(text_size = 20)
  
  return(out)
  
})

names(io$plots$cell_prevelences) = io$plots$args$prevelence_comps$plot_names

# Save plots as one single pdf
pdf(onefile = T, width = 25, height = 25,
    file = nf("cell_prevalences.pdf",io$output$cell_prevalences))

io$plots$cell_prevelences

dev.off()

# QUANTIFYING CELL TYPE PREVALENCES --------------------------------------------

# Convert anno counts to proportions
regions$main_anno = lapply(regions$main_anno, \(x) prop.table(x))
regions$fine_anno = lapply(regions$fine_anno, \(x) prop.table(x))

regions_long$main_anno = lapply(regions_long$main_anno, \(x) prop.table(x))
regions_long$fine_anno = lapply(regions_long$fine_anno, \(x) prop.table(x))

add_stats_to_plot = function(plot_stats, box_plot){

 out = box_plot +
   
    # Add stats to the boxplot
    ggpubr::stat_pvalue_manual(plot_stats,  label = "p_signif", tip.length = 0, size = 7) +
    
    # Add 10% spaces between the p-value labels and the plot border
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1)),
                                labels = scales::percent_format()) +
    # Add custom theme
    IMCfuncs::facetted_comp_bxp_theme(hide_x_labs = T) 

 return(out)
 
  }


# SURGERY PREVALENCES STATS ----------------------------------------------------

stats = compare_groups(
  df = regions,
  response_var = "fine_anno",
  nested_response_var = TRUE,
  comp_var = "surgery",
  group_var = "fine_anno",
  facetted_plot_dims = TRUE,
  stat_y_pos_multiplier = 3,
  p_signif_col = "p"
)




io$plots$args$comps_stats = tribble(
  
  ~df,         ~tbl,    ~response_var,   ~comp_var,   ~group_var,    ~unnest_y_levels,                ~plot_names,
  "regions",   NULL,    "main_anno",     "surgery",   "main_anno",    names(plot_colours$main_anno),  "all_regions_main",
  "regions",   NULL,    "fine_anno",     "surgery",   "fine_anno",    names(plot_colours$cell_anno),  "all_regions_fine",
)




# New version of the function
foo = purrr::pmap(plot_args, ~{
  
  if(is.null(..2)){
    
    df_to_use = get(..1, envir = globalenv())
    
  }else df_to_use = get(..1, envir = globalenv())[[..2]]
  
  return(
    list(
      comparison = ifelse(
        ..3 == ..5,
        yes = paste(..4, ..3, sep = "_"),
        no =  paste(..4, ..3, ..5, sep = "_")
      ) ,
      stats = compare_groups(
        df = df_to_use,
        response_var = ..3,
        nested_response_var = TRUE,
        unnest_levels = ..6,
        comp_var = ..4,
        group_var = ..5,
        facetted_plot_dims = TRUE,
        stat_y_pos_multiplier = 3,
        p_signif_col = "p"),
      
      boxplot = comp_boxplot(
        df = df_to_use,
        x = ..4,
        y = ..3,
        facet_by = ..3,
        unnest_y = TRUE,
        unnest_y_levels = ..6,
        fill = ..4,
        ylab = "Cell Proportions (%)",
        color_palette = plot_colours$surgery)
    )
  )
  
})

add_stats_to_plot(plot_stats = foo[[1]]$stats, box_plot = foo[[1]]$boxplot)

foo = pmap(io$plots$args$comps_stats, ~{
  
  if(is.null(..2)){
    
    df_to_use = get(..1, envir = globalenv())
    
  }else df_to_use = get(..1, envir = globalenv())[[..2]]  
  
  stats = compare_groups(
    df = df_to_use,
    response_var = ..3,
    nested_response_var = TRUE,
    comp_var = ..4,
    group_var = ..5,
    facetted_plot_dims = TRUE,
    stat_y_pos_multiplier = 3,
    p_signif_col = "p"
  )
  
  stats[[..3]] = factor(stats[[..3]], levels = ..6)
  
  
  bxp = comp_boxplot(
    df = df_to_use,
    x = ..4,
    y = ..3,
    facet_by = ..3,
    unnest_y = TRUE,
    unnest_y_levels = ..6,
    fill = ..4,
    ylab = "Cell Proportions (%)",
    color_palette = plot_colours$surgery) 
  
  out = add_stats_to_plot(stats, bxp) 
  
  return(out)
  
})

iwalk(io$plots$cell_comp_stats, ~{
  
  if(grepl("(?i)_fine$", x = .y)){
    
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
                     width = 20, height = 20)
    print(.x)
    
    dev.off()
    
  }else{
    
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
                     width = 15, height = 10)
    print(.x)
    
    dev.off()
                    
  }
  
  })
    
rm(stats, bxp)

# REGION/RESPONDER PREVALENCES STATS -------------------------------------------

# We can quantify the cell type prevalences stratified across various groups in
# order to see if any particular cell types are significantly different across 
# comparisons groups. The identities of the stratifications are as follows:
# 
# surgery (primary and recurrent)
# region_type (hypoxia, proliferative and immune)
# responder_type (up and down)

region_filt_props = as.list(setNames(
  object = c("immune", "prolif", "hypoxia"),
  nm = c("Immune Regions", "Prolif Regions", "Hypoxia Regions")
))

region_filt_props = lapply(region_filt_props, \(x) regions[regions[[x]] == T,])

region_filt_props$`Up Responders` = regions[regions[["responder_type"]] == "up",]
region_filt_props$`Down Responders` = regions[regions[["responder_type"]] == "down",]


io$plots$args$comps_stats = tribble(
  
  ~df,                     ~tbl,                ~response_var,      ~comp_var,       ~group_var,    ~unnest_y_levels,                ~plot_names,
  "region_filt_props",     "Immune Regions",    "main_anno",        "surgery",       "main_anno",    names(plot_colours$main_anno),  "immune_regions_main",
  "region_filt_props",     "Immune Regions",    "fine_anno",        "surgery",       "fine_anno",    names(plot_colours$cell_anno),  "immune_regions_fine",
  "region_filt_props",     "Prolif Regions",    "main_anno",        "surgery",       "main_anno",    names(plot_colours$main_anno),  "prolif_regions_main",
  "region_filt_props",     "Prolif Regions",    "fine_anno",        "surgery",       "fine_anno",    names(plot_colours$cell_anno),  "prolif_regions_fine",
  "region_filt_props",     "Hypoxia Regions",   "main_anno",        "surgery",       "main_anno",    names(plot_colours$main_anno),  "hypoxia_regions_main",
  "region_filt_props",     "Hypoxia Regions",   "fine_anno",        "surgery",       "fine_anno",    names(plot_colours$cell_anno),  "hypoxia_regions_fine",
  "region_filt_props",     "Up Responders",     "main_anno",        "surgery",       "main_anno",    names(plot_colours$main_anno),  "up_responders_main",
  "region_filt_props",     "Up Responders",     "fine_anno",        "surgery",       "fine_anno",    names(plot_colours$cell_anno),  "up_responders_fine",
  "region_filt_props",     "Down Responders",   "main_anno",        "surgery",       "main_anno",    names(plot_colours$main_anno),  "down_responders_main",
  "region_filt_props",     "Down Responders",   "fine_anno",        "surgery",       "fine_anno",    names(plot_colours$cell_anno),  "down_responders_fine",
  
  )
    

io$plots$cell_comp_region_stats = pmap(io$plots$args$comps_stats, ~{
  
  df_to_use = get(..1, envir = globalenv())[[..2]]  
  
  stats = compare_groups(
    df = df_to_use,
    response_var = ..3,
    nested_response_var = TRUE,
    comp_var = ..4,
    group_var = ..5,
    facetted_plot_dims = TRUE,
    stat_y_pos_multiplier = 3,
    p_signif_col = "p"
  )
  
  stats[[..3]] = factor(stats[[..3]], levels = ..6)
  
  
  bxp = comp_boxplot(
    df = df_to_use,
    x = ..4,
    y = ..3,
    facet_by = ..3,
    unnest_y = TRUE,
    unnest_y_levels = ..6,
    fill = ..4,
    ylab = "Cell Proportions (%)",
    color_palette = plot_colours$surgery) 
  
  out = add_stats_to_plot(stats, bxp, title = ..2) 
  
  return(out)

})

names(io$plots$cell_comp_region_stats) = io$plots$args$comps_stats$plot_names


iwalk(io$plots$cell_comp_region_stats, ~{
  
  if(grepl("(?i)_fine$", x = .y)){
    
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
                     width = 20, height = 20)
    print(.x)
    
    dev.off()
    
  }else{
    
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
                     width = 15, height = 10)
    print(.x)
    
    dev.off()
    
  }
  
})


## The Wilcoxon Signed Rank Test can be used to see if the population median of 
## the difference scores is equal to 0 or not

# We can also check the effect size measure for this test, using the 
# wilcoxonPairedR() function from the rcompanion package. 
# The function requires the data to be ordered specifically: The top half of the 
# data are the values for group1. The bottom half of the data are the values for group2. 
# Then the id variable needs to be in the same order for both groups.
# 
# An effect size above 0.7 is considered to be a large, 
# Anything above 0.5 is a moderate effect size,
# and anything above 0.3 is a small effect size.

# END --------------------------------------------------------------------------