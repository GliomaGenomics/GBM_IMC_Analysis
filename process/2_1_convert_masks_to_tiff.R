# CODE to convert PNG files to Tiff
# 
to_rename = list.files("data/combined_nuclear_cytoplasm/Masks/tiff", full.names = T)
new_names = str_replace_all(to_rename, "\\.png$", "\\.tiff")

library(magick)

lapply(new_names, function(x){

  foo = magick::image_read(x)

  image_write(foo, path = x, format = "tiff")

})