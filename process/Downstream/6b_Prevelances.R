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

# we are mainly interested in how things change in response to treatment in each
# patient we will initially seek to make the following comparisons:
# P vs R
# P vs R (grouped by patient)

se_comps <- tribble(
  ~df, ~response_var, ~comp_var, ~group_var, ~label,
  "regions", "entropy", "surgery", NULL, "surgery",
  "regions", "entropy", "surgery", "patient", "surgery_patient"
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

sh_stats <- purrr::map(sh_ent, ~ .x$stats)
names(sh_stats) <- sapply(sh_ent, \(x) x$label)

openxlsx::write.xlsx(
  x = sh_stats,
  file = nf("shannon_entropy_stats.xlsx", io$output$temp_out)
)

rm(calculate_props, se_comps, sh_stats, sh_ent)

# VISUALISE LABELLED CELL PROPORTIONS ------------------------------------------
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

rm(prop_comps, prop_plots, labelled_props)
# QUANTIFY LABELLED CELL PROPORTIONS -------------------------------------------
add_proportions <- function(x) {
  x %>%
    dplyr::mutate(
      total = sum(freq),
      prop = freq / total * 100
    ) %>%
    dplyr::select(-total)
}

regions <- regions %>%
  mutate(
    across(c("main_anno", "fine_anno"), ~ map(.x, add_proportions))
  )

compare_props <- function(df,
                          response_var,
                          comp_var,
                          stat = c("wilcox", "t"),
                          paired_test = FALSE,
                          p_adjust_method = "fdr",
                          signif_on = c("p_adj", "p"),
                          stat_y_pos_multiplier = 1) {
  df_name <- deparse(substitute(df))
  stat_used <- match.arg(stat, several.ok = FALSE)
  paired_out <- ifelse(paired_test, "(Paired)", "")
  comp_formula <- reformulate(response = "prop", termlabels = comp_var)

  long_data <- df %>%
    dplyr::select(patient, surgery, ROI, !!sym(response_var)) %>%
    tidyr::unnest_longer(col = response_var, values_to = "response") %>%
    tidyr::unnest_wider(response) %>%
    dplyr::group_by(label)

  df_grouped_by <- as.character(dplyr::groups(long_data))

  # print message start --
  cli::cli_div(theme = list(span.emph = list(color = "orange")))
  cli::cli_div(theme = list(span.strong = list(color = "darkred")))
  cli::cli_h1("Comparison Summary")
  cli::cli_text("")
  cli::cli_alert_info("{.strong DATA:} {.emph {df_name}}")
  cli::cli_alert_info("{.strong STAT:} {.emph {stat_used}{paired_out}, p.adjust = '{p_adjust_method}'}")
  cli::cli_alert_info("{.strong FORMULA:} {.emph {comp_formula[[2]]} ~ {comp_formula[[3]]}}")
  cli::cli_alert_info("{.strong GROUP(S):} {.emph '{response_var}' ({dplyr::n_groups(long_data)})}")
  cli::cli_h1("")
  # print message end --

  stat_test <- switch(stat_used,
    wilcox = rstatix::wilcox_test(
      data = long_data,
      formula = comp_formula,
      paired = paired_test,
      p.adjust.method = "none"
    ),
    t = rstatix::t_test(
      data = long_data,
      formula = comp_formula,
      paired = paired_test,
      p.adjust.method = "none"
    )
  )

  stat_test <- stat_test %>%
    rstatix::adjust_pvalue(
      output.col = "p_adj",
      method = p_adjust_method
    ) %>%
    rstatix::add_significance(
      p.col = signif_on,
      output.col = "p_signif"
    ) %>%
    rstatix::add_xy_position(x = comp_var) %>%
    dplyr::mutate(
      method = stat_used,
      paired_test = paired_test,
      p_adj_method = p_adjust_method,
      `.y.` = "prop",
      dataset = df_name,
      across(
        p_signif,
        ~ ifelse(. == "ns", glue::glue("{signif_on} = {round(stat_test[[signif_on]], 3)}"), .)
      )
    ) %>%
    dplyr::relocate(dataset)

  outliers <- sapply(df[[response_var]], \(x) x$prop, simplify = FALSE) %>%
    unlist(use.names = FALSE) %>%
    boxplot.stats()

  if (length(outliers$out) == 0) {
    stat_test %>% dplyr::mutate(y.position = max(stat_test$y.position))
  } else {
    stat_test %>% dplyr::mutate(y.position = min(outliers$out) * stat_y_pos_multiplier)
  }

  stat_test[[comp_var]] <- factor(stat_test$group1, levels = levels(long_data[[comp_var]]))


  return(
    list(
      data = long_data,
      stats = stat_test
    )
  )
}


create_plot <- function(label_data, label_stats, label_name) {
  ggplot(label_data, aes(x = surgery, y = prop, fill = surgery)) +
    geom_boxplot() +
    scale_fill_manual(values = plot_colours$surgery) +
    ylab("Cell Proportions (%)") +
    ggtitle(label_name) +
    IMCfuncs::facetted_comp_bxp_theme() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5)
    ) +
    stat_pvalue_manual(
      data = label_stats,
      label = "p_signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0,
      size = 7,
      fontface = "italic"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.1, 0.1))
    )
}

comp_patchwork <- function(comp_list, ncols = 2, nrows = NULL) {
  facets <- unique(comp_list$data$label)

  plots <- purrr::map(facets, ~ {
    label_data <- comp_list$data[comp_list$data$label == .x, ]
    label_stats <- comp_list$stats[comp_list$stats$label == .x, ]

    create_plot(label_data, label_stats, .x)
  })

  patchwork::wrap_plots(plots, ncol = ncols, nrow = nrows) +
    patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

comps <- list()

comps$main <- compare_props(
  df = regions,
  response_var = "main_anno",
  comp_var = "surgery",
  signif_on = "p"
)

comps$fine <- compare_props(
  df = regions,
  response_var = "fine_anno",
  comp_var = "surgery",
  signif_on = "p"
)

stats <- purrr::map(comps, ~ .x$stats)

openxlsx::write.xlsx(
  x = stats,
  file = nf("label_prop_stats.xlsx", io$output$temp_out)
)

svglite::svglite(
  filename = nf("surgery_main_anno.svg", io$output$temp_out),
  width = 15,
  height = 15
)
comp_patchwork(comps$main)
dev.off()


svglite::svglite(
  filename = nf("surgery_fine_anno.svg", io$output$temp_out),
  width = 20,
  height = 18
)
comp_patchwork(comps$fine, ncols = 4)
dev.off()

# SAVE OUTPUT DATA -------------------------------------------------------------
saveRDS(regions, nf("prevelance_data.rds", io$output$temp_out))

# END --------------------------------------------------------------------------
