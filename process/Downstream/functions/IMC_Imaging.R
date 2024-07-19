check_images <- function(image_list){
  
  if(!(class(image_list) %in% "CytoImageList")){
    
    stop("image_list is not of class 'CytoImageList'",call. = F)
  
    }
  single_image_list <- image_list@listData
  
  purrr::iwalk(single_image_list,~{
    
    total_image_pixels <- dim(.x)[1] * dim(.x)[2] 
    
    detected_pixels <- apply(.x, MARGIN = 3, sum, na.rm = T)
    
    if(!all(detected_pixels > 0)){
      
      channels_not_read = paste0(which(detected_pixels < 1), 
                                 collapse = ", ")
      
      print(glue::glue("{.y} channels not read in correctly:\t",
                 "{channels_not_read}")
            )
    } else print(glue::glue("{.y} read in correctly"))
  
  })
  
}


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


visualise_channel_pixels <- function(split_images, channels){
  
  outfilename <- glue::glue("{names(split_images)}.pdf")
  
  foo <- lapply(channel_vars, function(x){
    
    plotPixels(split_images, 
               colour_by = x,
               scale_bar = list(length = 100,
                                label = expression("100 " ~ mu * "m"),
                                colour = "white",
                                cex = 1),
               image_title = list(colour = "white",
                                  cex = 1,
                                  text = x),
               legend = NULL,
               # legend = list(colour_by.legend.cex = 2,
               #               colour_by.labels.cex = 2,
               #               colour_by.title.cex = 1.5,
               #               margin = 50),
               display = "single",
               return_plot = TRUE)
  })
  
  foo <- foo %>% purrr::flatten(.) %>% purrr::flatten(.)
  
  pdf(file= file.path("Outputs/Channel_pixel_intesity_images/",outfilename))
  
  purrr::walk(foo, ~print(.x))
  
  dev.off()
  
}


visualise_masks <- function(split_mask, output_dir = file.path("Outputs/Segmentation_Masks")){
  
  plot_out_name <- glue::glue(file.path(output_dir, "{names(split_mask)}.png"))
  
  plotCells(split_mask, 
            scale_bar = list(length = 100,
                             label = expression("100 " ~ mu * "m"),
                             colour = "red",
                             cex = 3),
            image_title = list(colour = "red",
                               cex = 3),
            legend = NULL,
            display = "single", 
            save_plot = list(filename = plot_out_name)
            )
}


clean_mask_names <- function(clean_dir = file.path("Outputs/Segmentation_Masks")){
  
  replace_names <- list.files(clean_dir, full.names = T,
                              pattern = "_\\d+\\.png$")
  
  new_name <- stringr::str_replace_all(replace_names, 
                                       "_\\d+\\.png$", 
                                       "\\.png")
  
  file.rename(from = replace_names, to = new_name)
  
  
}


# read_images <- function(sample_pattern){
#   
#   images <- loadImages("Data/img/", pattern = sample_pattern)
#   
#   # Check if the images were read in correctly 
#   check_images(images)
#   
#   masks <- loadImages("data/masks/", pattern = sample_pattern,
#                       as.is = TRUE)
#   
#   purrr::walk(masks@listData, ~print(.x@.Data[1:10,1:5]))
#   
#   
#   split_masks <- split_CytoImageList(masks)
#   
#   walk(split_masks, ~ visualise_masks(split_mask = .x))
#   
#   clean_mask_names()
#   
#   return(list(images = images, masks = split_masks))
#   
# }
