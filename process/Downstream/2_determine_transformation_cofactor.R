# PACKAGES ---------------------------------------------------------------------
library(imcRtools)
library(dittoSeq)
library(cytomapper)
library(dplyr)
library(stringr)
library(magrittr)

# I/O --------------------------------------------------------------------------
io <- list(
  
  inputs = list(
    functions = "process/Downstream/functions",
    data = file.path("data/downstream/raw")
    ),
  
  output= list(outdir = "outputs/QC/cofactor_tuning")
)

if(!dir.exists(io$output$outdir)) dir.create(io$output$outdir)

# LOAD DATA --------------------------------------------------------------------
spe = readRDS(file.path(io$inputs$data, "spe.rds"))

# ITERATE OVER COFACTORS -------------------------------------------------------
cofactors = c(0.1,0.5,seq(1:10))

for (i in seq_along(cofactors)) {
  
  assay(spe, glue::glue("exprs_cofactor_{cofactors[[i]]}")) = asinh(counts(spe)/cofactors[[i]])
  
}

cofactors = assayNames(spe)[grepl("cofactor", assayNames(spe))]

# PLOT OUTPUTS -----------------------------------------------------------------
plot_marker_exprs = function(spe, assay){
  
  exprs = as.data.frame(t(spe@assays@data[[assay]]))
  
  out = purrr::imap(exprs, ~{
    
    data.frame(.x) %>%
    ggplot(aes(x = .x)) +
      geom_histogram(fill="lightblue", bins = 1000) +
      xlab(glue::glue("{.y}")) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text = element_text(face = "bold"),
        legend.position = "none"
      )})
  
  out =  purrr::map(out, ~ .x + ggtitle(assay))
  
  return(out)
  
  }

cofactor_plots = purrr::map(cofactors, ~ plot_marker_exprs(spe = spe, assay = .x))

names(cofactor_plots) = cofactors


purrr::iwalk(cofactor_plots, ~{
  
  
  pdf(file = file.path(io$output$outdir,glue::glue("{.y}.pdf")), 
      onefile = T)
  
  print(.x)
  
  dev.off()
  
  
})

# END --------------------------------------------------------------------------
