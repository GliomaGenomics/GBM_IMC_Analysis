# For the first bit of exploratory analysis we can check the cell type prevalence
# across various stratifications, starting with the the specific tissue samples.
# In order to do this we will first obtain counts for each of the cells and then
# subset/aggregate them to show the cell prevalence across different groupings.

# Previously, I clustered and annotated the imaging mass cytometry single cells.
# During this process a number of cells remained un-annotated as their
# identity could not be easily inferred using one of more of the cell type markers
# which were used. Also, in some cases the abundance of the metal isotopes corresponding
# to the cell type markers was indiscriminate in distinguishing one particular cell
# type or cell state and so was left "unknown".
# This initial cell type prevalence analysis will include these "unknown" markers.

# options(scipen = 999)

# PACKAGES ---------------------------------------------------------------------
library(SpatialExperiment)
# library(dittoSeq)
# library(imcRtools)
library(ggplot2)
library(viridis)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(IMCfuncs)

# I/O --------------------------------------------------------------------------
io <- list(
  inputs = list(
    input_dir = "outputs/cell_phenotyping"
  ),
  outputs = list(
    out_dir = "outputs/cell_prevalences"
  ),
  plots = list()
)

# create out directory
ndirs(io$outputs)

# create time-stamped output directory
io$outputs$temp_out <- nd(path = io$outputs$out_dir)

# obtain the most-recent data from a time-stamped directory
io$inputs$data <- list.files(io$inputs$input_dir, pattern = "^[T0-9-]+$")
io$inputs$data <- io$inputs$data[order(io$inputs$input_dir, decreasing = TRUE)][[1]]

io$inputs$data <- list.files(
  path = file.path(io$inputs$input_dir, io$inputs$data),
  pattern = "(?i)^spe_downstream",
  full.names = TRUE
)

# LOAD DATA  -------------------------------------------------------------------
spe <- readRDS(io$inputs$data)

# CREATE PREVELENCE DATA  ------------------------------------------------------
regions <- as.data.frame(colData(spe)) %>%
  distinct(
    sample_id,
    patient,
    surgery,
    ROI,
    responder_type,
    hist_region_label = region_type,
    imc_region_label = region_type_new
  ) %>%
  filter(!sample_id %in% paste0("67Prim_00", 4:7)) # exclude additional 67Prim images


regions <- as.data.frame(spe@colData) %>%
  group_by(sample_id) %>%
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
  left_join(regions, ., by = "sample_id")


# Add labels to region level metadata
regions$anno_labels <- as.data.frame(colData(spe)) %>%
  dplyr::filter(sample_id %in% regions$sample_id) %>%
  dplyr::mutate(across(c("main_anno_v2", "manual_gating"), as.character)) %>%
  dplyr::select(sample_id, cell_main = main_anno_v2, cell_fine = manual_gating) %>%
  tidyr::drop_na() %>%
  split(.$sample_id) %>%
  lapply(\(x) setNames(x$cell_fine, x$cell_main))

count_labels <- function(annotation, labels) {
  count_tbl <- table(factor(annotation, levels = labels))
  count_tbl %>%
    as.data.frame() %>%
    dplyr::rename(label = 1, freq = 2)
}

regions$main_anno <- purrr::map(regions$anno_labels, ~ {
  count_labels(
    annotation = names(.x),
    labels = names(spe@metadata$labels$cell_types)
  )
})

regions$fine_anno <- purrr::map(regions$anno_labels, ~ {
  count_labels(
    annotation = .x,
    labels = unlist(spe@metadata$labels$cell_types, use.names = FALSE)
  )
})


factor(regions$surgery, names(spe@metadata$v2_colours$dataset_pheno)) %>%
  droplevels()


regions <- regions %>%
  mutate(
    across(
      surgery, ~ factor(.x, levels = names(plot_colours$surgery))
    ),
    across(
      responder_type, ~ factor(.x, levels = c("up", "down"))
    ),
    across(
      patient, ~ factor(.x, levels = names(plot_colours$patient))
    ),
    across(
      ROI, ~ factor(.x, levels = names(plot_colours$ROI))
    ),
    across(
      c(imc_region_label, hist_region_label),
      ~ factor(.x, levels = names(plot_colours$region_type))
    )
  )

rm(count_labels)

# PLOT LABELLED/UNDEFINED CELL COUNTS ------------------------------------------
svglite::svglite(
  nf("cell_labelling_counts.svg", io$output$temp_out),
  width = 15, height = 10
)

regions %>%
  select(patient, surgery, labelled, undefined) %>%
  tidyr::pivot_longer(
    cols = c("labelled", "undefined"),
    names_to = "anno",
    values_to = "count"
  ) %>%
  group_by(patient, surgery, anno) %>%
  summarise(count = sum(count), .groups = "keep") %>%
  ungroup() %>%
  ggplot(aes(x = anno, y = count, fill = anno)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(surgery ~ patient) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("cell count") +
  scale_y_continuous(breaks = seq(0, 14000, 2000)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.5) +
  IMCfuncs::facetted_cell_prop_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()

# RE-LABEL REGIONS BASED ON CELL COUNTS ----------------------------------------
# Each image represents a tissue section that was selected based on
# positive immunohistochemistry markers for either immune-rich, proliferative or
# hypoxic regions. Each patient and surgery sample pairing was selected to have at least
# one of each region type. We can use the imaging mass cytometry labelled cell-type
# and state prevalences to corroborate if the histochemistry region labels were
# indeed reflective of the cells initially labelled based on IHC single cell-type/state
# canonical markers. The region types and therefore the labels are not mutually
# exclusive, e.g, any given region can be high (as determined by an appropriate threshold)
# in immune cells and also comprise cells that are in a highly proliferate/hypoxic state.

regions_labels <- regions %>%
  select(sample_id, total_cells, starts_with("region_")) %>%
  tidyr::pivot_longer(
    cols = starts_with("region_"),
    names_to = "region",
    names_prefix = "region_",
    values_to = "cells",
  ) %>%
  mutate(percent = round(cells / total_cells * 100)) %>%
  mutate(across(region, ~ factor(.x, levels = c("immune", "prolif", "hypoxia"))))

add_region_cols <- function(regions,
                            count_col,
                            region_labels,
                            immune_min = 20,
                            state_min = 50) {
  out <- as.list(setNames(
    levels(regions[[region_labels]]),
    levels(regions[[region_labels]])
  ))

  out <- lapply(out, \(x){
    enrichment <- ifelse(regions[[region_labels]] == x, regions[[count_col]], 0)
    if (grepl("(?i)immune", x)) enrichment >= immune_min else enrichment >= state_min
  })

  return(cbind(regions, dplyr::bind_cols(out)))
}

# Apply new region labels based on min defined thresholds
io$plots$args <- list(immune_min = 20, state_min = 25)

regions_labels <- add_region_cols(
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

regions_labels$immune[regions_labels$sample_id == "82Prim_002" &
  regions_labels$region == "immune"] <- TRUE

# Add colour info for facets based on the region labels
regions_labels$facet_background <- plot_colours$patient[str_extract(regions_labels$sample_id, "^\\d{2}")]
regions_labels$facet_text <- IMCfuncs::get_text_color(regions_labels$facet_background)

rm(add_region_cols)
# PLOT RE-LABELLED REGIONS -----------------------------------------------------
svglite::svglite(
  filename = nf("regions_labels.svg", io$output$temp_out),
  width = 12,
  height = 25
)

regions_labels %>%
  ggplot(aes(x = region, y = percent, fill = region, group = sample_id)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_abline(aes(intercept = io$plots$args$immune_min, slope = 0, color = "Immune Threshold"),
    lty = 2
  ) +
  geom_abline(aes(intercept = io$plots$args$state_min, slope = 0, color = "State Threshold"),
    lty = 2
  ) +
  ggh4x::facet_wrap2(
    facets = " sample_id",
    ncol = 3,
    strip = ggh4x::strip_themed(
      background_x = lapply(
        split(regions_labels, regions_labels$sample_id), \(x) {
          element_rect(fill = unique(x[["facet_background"]]))
        }
      ),
      text_x = lapply(
        split(regions_labels, regions_labels$sample_id), \(x) {
          element_text(
            colour = unique(x[["facet_text"]]),
            face = "bold"
          )
        }
      )
    )
  ) +
  xlab("") +
  ylab("Percentage") +
  scale_fill_manual(values = plot_colours$region_type[levels(regions_labels$region)]) +
  scale_color_manual(values = c(
    "Immune Threshold" = "lightpink2",
    "State Threshold" = "skyblue1"
  )) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  )

dev.off()



regions <- regions_labels %>%
  group_by(sample_id) %>%
  summarise(
    immune = any(immune),
    prolif = any(prolif),
    hypoxia = any(hypoxia),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  left_join(regions, ., by = "sample_id")


plot_data <- regions %>%
  select(sample_id, patient, surgery, ROI, immune, prolif, hypoxia) %>%
  mutate(
    facet_background = plot_colours$patient[patient],
    facet_text = IMCfuncs::get_text_color(facet_background)
  ) %>%
  tidyr::pivot_longer(c("immune", "prolif", "hypoxia"),
    names_to = "region",
    values_to = "present"
  ) %>%
  mutate(
    patient_surgery = paste(patient, surgery, sep = "_"),
    across(present, ~ factor(ifelse(.x, "yes", "no"), levels = c("yes", "no")))
  )

svglite::svglite(
  filename = nf("regions_per_sample.svg", io$output$temp_out),
  width = 15,
  height = 30
)

plot_data %>%
  ggplot(aes(x = region, y = ROI, fill = present)) +
  geom_point(shape = 21, size = 45, color = "grey75") +
  xlab("") +
  ylab("") +
  ggh4x::facet_wrap2(
    facets = "patient_surgery",
    ncol = 2,
    strip = ggh4x::strip_themed(
      background_x = lapply(
        split(plot_data, plot_data$patient_surgery), \(x) {
          element_rect(fill = unique(x[["facet_background"]]))
        }
      ),
      text_x = lapply(
        split(plot_data, plot_data$patient_surgery), \(x) {
          element_text(
            colour = unique(x[["facet_text"]]),
            face = "bold"
          )
        }
      )
    )
  ) +
  scale_fill_manual(
    values = c("yes" = "#FCFDBFFF", "no" = "#000004FF")
  ) +
  IMCfuncs::grouped_comp_bxp_theme()

dev.off()

rm(plot_data, regions_labels)
# SHANNON ENTROPY --------------------------------------------------------------
# Here we are seeking to measure the intra-patient heterogeneity by
# using Shannon entropy (H) and the annotated cell type labels. In order to
# account for the different number of cells per sample, we will
# sub-sample n (1000) cells from each patient sample (i). This will be repeated
# over ten rounds, to obtain the Shannon entropy of the cell type frequencies
# in each patient sample region.

calculate_props <- function(x, labels, sample_size = 1000, rounds = 3) {
  base <- setNames(object = rep(0, length(labels)), nm = labels)

  out <- vector("list", rounds)

  for (i in seq_len(rounds)) {
    exp <- prop.table(table(sample(x, sample_size)))
    exp <- set_names(as.numeric(exp), names(exp))

    exp <- c(base[setdiff(names(base), names(exp))], exp)

    out[[i]] <- philentropy::H(exp, unit = "log2")
  }

  return(unlist(out, use.names = FALSE))
}

# Calculate the entropy for each cell type label (per sample ID)
set.seed(123)
regions$entropy <- map(regions$anno_labels, ~ calculate_props(
  .x,
  labels = levels(spe$manual_gating),
  sample_size = 1000,
  rounds = 10
))

# plot Shannon entropy results

# convert the regions data.frame to a long format based on the region type
# regions_long <- regions %>%
#   tidyr::pivot_longer(
#     immune:hypoxia,
#     names_to = "region",
#     values_to = "region_present"
#   ) %>%
#   filter(region_present == TRUE)

# we are mainly interested in how things change in response to treatment in each
# patient we will initially seek to make the following comparisons:
# P vs R
# P vs R (grouped by patient)

se_comps <- tribble(
  ~df, ~response_var, ~comp_var, ~group_var, ~label,
  "regions", "entropy", "surgery", NULL, "surgery",
  "regions", "entropy", "surgery", "patient", "surgery_patient"
  # "regions",          "entropy",         "surgery",       "responder_type",  "surgery_responder",
  # "regions_long",     "entropy",         "surgery",       "region"           "surgery_region"
)

sh_ent <- pmap(se_comps, ~ {
  df_to_use <- get(..1, envir = globalenv())

  stats <- compare_groups(
    df = df_to_use,
    df_name = ..1,
    response_var = ..2,
    nested_response_var = TRUE,
    comp_var = ..3,
    group_var = ..4,
    stat = "wilcox",
    p_adjust_method = "fdr",
    facetted_plot_dims = TRUE
  )

  out_boxplot <- comp_boxplot(
    df = tidyr::unnest(df_to_use, ..2),
    x = ..3,
    y = ..2,
    facet_by = ..4,
    fill = ..3,
    color_palette = plot_colours$surgery
  )

  out_boxplot <- out_boxplot +
    ggpubr::stat_pvalue_manual(stats, label = "p_signif", tip.length = 0, size = 7) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
    ylab("Shannon Entropy") +
    IMCfuncs::facetted_comp_bxp_theme(hide_x_labs = T)

  return(
    list(
      stats = stats,
      boxplot = out_boxplot,
      label = ..5
    )
  )
})

purrr::walk(sh_ent, ~ {
  svglite::svglite(
    filename = nf(paste0(.x$label, "_sh_entropy.svg"), io$output$temp_out),
    width = 15,
    height = 10
  )

  print(.x$boxplot)

  dev.off()
})

rm(calculate_props, se_comps)

# PLOT CELL PROPORTIONS --------------------------------------------------------
labelled_props <- function(df,
                           label_col,
                           colgroups = NULL,
                           rowgroups = NULL,
                           color_palette,
                           y_axis_title = "labelled cell proportions",
                           plot_title = "",
                           plot_subtitle = "",
                           x_axis_title = "") {
  plot_data <- df %>%
    tidyr::unnest_longer(col = label_col, values_to = "count") %>%
    tidyr::unnest_wider(count)

  plot_data %>%
    ggplot2::ggplot(aes(x = surgery, y = freq, fill = label, colour = label)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(
      rows = if (!is.null(rowgroups)) ggplot2::vars(!!ggplot2::sym(rowgroups)) else NULL,
      cols = if (!is.null(colgroups)) ggplot2::vars(!!ggplot2::sym(colgroups)) else NULL
    ) +
    ylab(y_axis_title) +
    xlab(x_axis_title) +
    labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::scale_color_manual(values = color_palette) +
    IMCfuncs::facetted_cell_prop_theme()
}

prop_comps <- tribble(
  ~df, ~label_col, ~colgroups, ~rowgroups, ~color_palette,
  "regions", "main_anno", NULL, NULL, spe@metadata$v2_colours$cell_groups,
  "regions", "main_anno", "patient", NULL, spe@metadata$v2_colours$cell_groups,
  "regions", "main_anno", "patient", "ROI", spe@metadata$v2_colours$cell_groups,
  "regions", "fine_anno", NULL, NULL, spe@metadata$v2_colours$cells,
  "regions", "fine_anno", "patient", NULL, spe@metadata$v2_colours$cells,
  "regions", "fine_anno", "patient", "ROI", spe@metadata$v2_colours$cells
)

prop_plots <- pmap(prop_comps, ~ {
  labelled_props(
    df = get(..1, envir = globalenv()),
    label_col = ..2,
    colgroups = ..3,
    rowgroups = ..4,
    color_palette = ..5
  )
})

pdf(
  onefile = T, width = 25, height = 25,
  file = nf("label_proportions.pdf", io$outputs$temp_out)
)

prop_plots

dev.off()

# QUANTIFYING CELL TYPE PREVALENCES --------------------------------------------
# Convert anno counts to proportions
regions$main_anno <- lapply(regions$main_anno, \(x) prop.table(x))
regions$fine_anno <- lapply(regions$fine_anno, \(x) prop.table(x))

regions_long$main_anno <- lapply(regions_long$main_anno, \(x) prop.table(x))
regions_long$fine_anno <- lapply(regions_long$fine_anno, \(x) prop.table(x))

add_stats_to_plot <- function(plot_stats, box_plot) {
  out <- box_plot +

    # Add stats to the boxplot
    ggpubr::stat_pvalue_manual(plot_stats, label = "p_signif", tip.length = 0, size = 7) +

    # Add 10% spaces between the p-value labels and the plot border
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.1, 0.1)),
      labels = scales::percent_format()
    ) +
    # Add custom theme
    IMCfuncs::facetted_comp_bxp_theme(hide_x_labs = T)

  return(out)
}

# SURGERY PREVALENCES STATS ----------------------------------------------------
stats <- compare_groups(
  df = regions,
  response_var = "fine_anno",
  nested_response_var = TRUE,
  comp_var = "surgery",
  group_var = "fine_anno",
  facetted_plot_dims = TRUE,
  stat_y_pos_multiplier = 3,
  p_signif_col = "p"
)

io$plots$args$comps_stats <- tribble(
  ~df, ~tbl, ~response_var, ~comp_var, ~group_var, ~unnest_y_levels, ~plot_names,
  "regions", NULL, "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "all_regions_main",
  "regions", NULL, "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "all_regions_fine",
)

# New version of the function
foo <- purrr::pmap(plot_args, ~ {
  if (is.null(..2)) {
    df_to_use <- get(..1, envir = globalenv())
  } else {
    df_to_use <- get(..1, envir = globalenv())[[..2]]
  }

  return(
    list(
      comparison = ifelse(
        ..3 == ..5,
        yes = paste(..4, ..3, sep = "_"),
        no = paste(..4, ..3, ..5, sep = "_")
      ),
      stats = compare_groups(
        df = df_to_use,
        response_var = ..3,
        nested_response_var = TRUE,
        unnest_levels = ..6,
        comp_var = ..4,
        group_var = ..5,
        facetted_plot_dims = TRUE,
        stat_y_pos_multiplier = 3,
        p_signif_col = "p"
      ),
      boxplot = comp_boxplot(
        df = df_to_use,
        x = ..4,
        y = ..3,
        facet_by = ..3,
        unnest_y = TRUE,
        unnest_y_levels = ..6,
        fill = ..4,
        ylab = "Cell Proportions (%)",
        color_palette = plot_colours$surgery
      )
    )
  )
})

add_stats_to_plot(plot_stats = foo[[1]]$stats, box_plot = foo[[1]]$boxplot)

foo <- pmap(io$plots$args$comps_stats, ~ {
  if (is.null(..2)) {
    df_to_use <- get(..1, envir = globalenv())
  } else {
    df_to_use <- get(..1, envir = globalenv())[[..2]]
  }

  stats <- compare_groups(
    df = df_to_use,
    response_var = ..3,
    nested_response_var = TRUE,
    comp_var = ..4,
    group_var = ..5,
    facetted_plot_dims = TRUE,
    stat_y_pos_multiplier = 3,
    p_signif_col = "p"
  )

  stats[[..3]] <- factor(stats[[..3]], levels = ..6)


  bxp <- comp_boxplot(
    df = df_to_use,
    x = ..4,
    y = ..3,
    facet_by = ..3,
    unnest_y = TRUE,
    unnest_y_levels = ..6,
    fill = ..4,
    ylab = "Cell Proportions (%)",
    color_palette = plot_colours$surgery
  )

  out <- add_stats_to_plot(stats, bxp)

  return(out)
})

iwalk(io$plots$cell_comp_stats, ~ {
  if (grepl("(?i)_fine$", x = .y)) {
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
      width = 20, height = 20
    )
    print(.x)

    dev.off()
  } else {
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
      width = 15, height = 10
    )
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

region_filt_props <- as.list(setNames(
  object = c("immune", "prolif", "hypoxia"),
  nm = c("Immune Regions", "Prolif Regions", "Hypoxia Regions")
))

region_filt_props <- lapply(region_filt_props, \(x) regions[regions[[x]] == T, ])

region_filt_props$`Up Responders` <- regions[regions[["responder_type"]] == "up", ]
region_filt_props$`Down Responders` <- regions[regions[["responder_type"]] == "down", ]


io$plots$args$comps_stats <- tribble(
  ~df, ~tbl, ~response_var, ~comp_var, ~group_var, ~unnest_y_levels, ~plot_names,
  "region_filt_props", "Immune Regions", "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "immune_regions_main",
  "region_filt_props", "Immune Regions", "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "immune_regions_fine",
  "region_filt_props", "Prolif Regions", "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "prolif_regions_main",
  "region_filt_props", "Prolif Regions", "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "prolif_regions_fine",
  "region_filt_props", "Hypoxia Regions", "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "hypoxia_regions_main",
  "region_filt_props", "Hypoxia Regions", "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "hypoxia_regions_fine",
  "region_filt_props", "Up Responders", "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "up_responders_main",
  "region_filt_props", "Up Responders", "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "up_responders_fine",
  "region_filt_props", "Down Responders", "main_anno", "surgery", "main_anno", names(plot_colours$main_anno), "down_responders_main",
  "region_filt_props", "Down Responders", "fine_anno", "surgery", "fine_anno", names(plot_colours$cell_anno), "down_responders_fine",
)


io$plots$cell_comp_region_stats <- pmap(io$plots$args$comps_stats, ~ {
  df_to_use <- get(..1, envir = globalenv())[[..2]]

  stats <- compare_groups(
    df = df_to_use,
    response_var = ..3,
    nested_response_var = TRUE,
    comp_var = ..4,
    group_var = ..5,
    facetted_plot_dims = TRUE,
    stat_y_pos_multiplier = 3,
    p_signif_col = "p"
  )

  stats[[..3]] <- factor(stats[[..3]], levels = ..6)


  bxp <- comp_boxplot(
    df = df_to_use,
    x = ..4,
    y = ..3,
    facet_by = ..3,
    unnest_y = TRUE,
    unnest_y_levels = ..6,
    fill = ..4,
    ylab = "Cell Proportions (%)",
    color_palette = plot_colours$surgery
  )

  out <- add_stats_to_plot(stats, bxp, title = ..2)

  return(out)
})

names(io$plots$cell_comp_region_stats) <- io$plots$args$comps_stats$plot_names


iwalk(io$plots$cell_comp_region_stats, ~ {
  if (grepl("(?i)_fine$", x = .y)) {
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
      width = 20, height = 20
    )
    print(.x)

    dev.off()
  } else {
    svglite::svglite(nf(paste0(.y, ".svg"), io$output$cell_prevalence_comps),
      width = 15, height = 10
    )
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
