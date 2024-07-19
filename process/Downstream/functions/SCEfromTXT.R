SCEfromTXT <- function (x, pattern = ".txt$", 
          metadata_cols = c("Start_push", "End_push", 
                            "Pushes_duration", "X", "Y", "Z"), 
          verbose = TRUE, 
          read_metal_from_filename = TRUE,
          filter_channels = NULL) {
  
  if (all(is.character(x)) & length(x) == 1) {
    
    if (!dir.exists(x)) {
      stop("Path does not exist.")
    }
    
    cur_names <- list.files(x, pattern = pattern, full.names = FALSE)
    
    if (length(cur_names) == 0) {
      stop("Files could not be read in.")
    }
    
    if (read_metal_from_filename) {
      cur_names <- str_extract(cur_names, "[A-Z]{1}[a-z]{0,1}[0-9]{2,3}")
    }else {
      cur_names <- sub("\\.[^.]*$", "", basename(cur_names))
    }
    
    txt_list <- list.files(x, pattern = pattern, full.names = TRUE)
    
    txt_list <- lapply(txt_list, readr::read_delim, delim = "\t", 
                       col_types = cols())
    
    txt_list <- lapply(txt_list, as.data.frame)
    
    names(txt_list) <- cur_names
    
  }else if (is.list(x)) {
    
    if (is.null(names(x))) {
      stop("If 'x' is a list, it needs to be named.")
    }
    cur_names <- names(x)
    txt_list <- lapply(x, as.data.frame)
  
    }else stop("Input 'x' is not of the correct format.")

  imcRtools:::.valid.readSCEfromTXT.input(txt_list, cur_names, metadata_cols, 
                                          verbose, read_metal_from_filename)
  cur_out <- do.call(rbind, txt_list)
  
  cell_meta <- DataFrame(cur_out[metadata_cols])
  
  if (read_metal_from_filename) {
    cell_meta$sample_id <- str_extract(rownames(cell_meta), 
                                       "^[A-Z]{1}[a-z]{0,1}[0-9]{2,3}")
    cell_meta$sample_metal <- str_extract(cell_meta$sample_id, 
                                          "^[A-Z]{1}[a-z]{0,1}")
    cell_meta$sample_mass <- str_extract(cell_meta$sample_id, 
                                         "[0-9]{2,3}$")
  }
  else {
    cell_meta$sample_id <- str_split(rownames(cell_meta), 
                                     "\\.", simplify = TRUE)[, 1]
  }
  
  cur_counts <- cur_out[grepl("[A-Z]{1}[a-z]{0,1}[0-9]{2,3}", 
                              colnames(cur_out))]
  cur_counts <- t(cur_counts)
  channel_name <- str_extract(rownames(cur_counts), "[A-Z]{1}[a-z]{0,1}[0-9]{2,3}Di")
  channel_meta <- DataFrame(channel_name = channel_name, marker_name = sub("Di", 
                                                                           "", channel_name))
  sce <- SingleCellExperiment(assays = list(counts = cur_counts))
  colData(sce) <- cell_meta
  rowData(sce) <- channel_meta
  return(sce)
}
