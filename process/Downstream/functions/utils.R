.new_file = function(name, path = NULL, 
                     timestamp_format = "%Y-%m-%dT%H-%M-%S",
                     ext = c("pdf","svg","png","csv","xlsx", "rds")) {
    
    if(missing(name)) cli::cli_abort("please supply a name")
    if(!is.null(path)) if(!dir.exists(path)) cli::cli_abort("path does not exist!")
    
    file_extension = paste0(".", match.arg(ext, several.ok = FALSE))

    time_stamp = format(Sys.time(), timestamp_format)
    
    outfilename = paste0(paste(name, time_stamp, sep = "_"), file_extension)
    
    if(is.null(path)) return(outfilename) else return(file.path(path, outfilename))
    
}
