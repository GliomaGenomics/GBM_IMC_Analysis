getmarkers = function(
    pheno = NULL,
    unique = TRUE,
    markers = spe@metadata$markers
) {
    
    all_markers = unlist(markers)
    
    marker_pheno = gsub("\\d+$","", names(all_markers)) |> 
        strsplit(split = "\\.")
    
    if(!is.null(pheno)){
        
        if(length(pheno) > 1) pheno = paste0("(", paste0(pheno, collapse = "|"), ")")
        
        lookup = paste0("(?i)^", pheno)
        
        outnames =  unlist(lapply(marker_pheno, \(x){ 
            
          name =  grep(lookup, x, value = T)

          if(length(name) == 1) return(name)
          if(length(name) > 1) return(paste0(name, collapse = "_"))
          
        }))
        
        out = all_markers[vapply(marker_pheno, \(x) any(grepl(lookup, x)), logical(1))]
        
        names(out) = outnames
        
    } else out = all_markers
    
    if(unique) return(out[!duplicated(out)]) else return(out)
    
}
