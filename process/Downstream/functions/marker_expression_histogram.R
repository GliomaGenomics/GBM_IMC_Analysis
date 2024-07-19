# This function plot the distribution of expression values as a histogram 
# for a given marker and also allows a cuttof to be placed, which can be useful 
# for manual gating apporches. 
plot_marker_exprs <- function(spe_object, assay="exprs", marker,
                              
                              quantile_cuttoff = 0.75, 
                              return_plot = TRUE,
                              save_plot = FALSE,
                              outdir = NULL, filename = NULL){
    
    marker_exprs <- data.frame(exprs = assay(spe_object, assay)[marker,]) 
    
    cuttoff_value <-  signif(quantile(marker_exprs$exprs, quantile_cuttoff),
                             digits = 3)
    
    above_threshold <- length(which(marker_exprs$exprs >= cuttoff_value))
    
    percent_total <- round(above_threshold/nrow(marker_exprs) * 100, 2)
    
    x_axis_scale_breaks = pretty(marker_exprs$exprs)
    
    outplot <- ggplot(marker_exprs, aes(x = exprs)) +
        geom_histogram(bins = 100, binwidth =  0.05,
                       color = "grey20", fill = "skyblue") +
        geom_vline(xintercept = cuttoff_value, color = "maroon", size = 1, 
                   linetype = "dashed") +
        xlab(glue::glue("assay({assay}) Expression")) +
        ylab("Cells") +
        ggtitle(glue::glue("{marker} Expression"),
                subtitle = glue::glue("cuttoff (quantile) = {cuttoff_value} ({quantile_cuttoff})\nncells = {above_threshold} ({percent_total}%)")) +
        scale_x_continuous(breaks = x_axis_scale_breaks) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
              plot.subtitle = element_text(size = 14,colour = "grey30", 
                                           face = "italic", hjust = 0.5)) 
    
    cli::cli_alert_info("{marker} threshold = {cuttoff_value}")
    
    if(save_plot){
       
        if(is.null(filename)){
            
            ggsave(filename = file.path(outdir, glue::glue("{marker}.png")))
            
        }else{
            
            ggsave(filename = file.path(outdir, glue::glue("{filename}.png")))
            
        } 
        
        
    }
    
    if(return_plot) return(outplot)

}

