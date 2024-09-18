library(factoextra)
library(FactoMineR)
library(googledrive)
library(raster)
library(RStoolbox)
library(tidyverse)
library(viridis)

setwd("/home/alan/Dropbox/Alan/INMA/Dados/MNE_BHRD")

###
repo <- file.path("/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/South_America")
list <- list.files(repo, pattern = ".tif$", full.names = TRUE, recursive = TRUE)

var_names <- c("bio1", paste0("bio", 10:19), paste0("bio", 2:9),
               "gdd0", "gdd10", "gdd5", "gls", "gsp", "gst", "npp")

pres <- grep("1981-2010", list, value = TRUE)
r_pres <- stack(pres)
names(r_pres) <- var_names

### run via modified function  ####
raster_pres = r_pres
#raster_proj = r_future
cum_sum = .95
graphics = FALSE
graphic_names = ""
graphic_var_pres_names = names(r_pres)
prefix_pca_raster_pres = "pres"
prefix_pca_raster_proj = "future"
path_output_pca_pres = file.path(getwd(), "pres")
path_output_pca_proj = file.path(getwd(), "proj")

# file paths
f_in <- pres
f_out <- tempfile(fileext = ".tif")

# input and output rasters
r_in <- stack(f_in); names(r_in) <- var_names
r_out <- raster(r_in)

# blocks
b <- blockSize(r_in)
print(b)

# open files
r_in <- readStart(r_in)
r_out <- writeStart(r_out, filename = f_out)

for (i in seq_along(b$row)) {
  # read values for block
  # format is a matrix with rows the cells values and columns the layers
  if(i == 1){
    x <- raster::xFromCol(r_in)
    y <- raster::yFromRow(r_in, row = seq(from = b$row[i], to = b$row[i]+b$nrows[i]-1))
    v <- getValues(r_in, row = b$row[i], nrows = b$nrows[i])
    v <- cbind(x, y = rep(y, each = length(x)), v)
  } else {
    y <- raster::yFromRow(r_in, row = seq(from = b$row[i], to = b$row[i]+b$nrows[i]-1))
    v2 <- getValues(r_in, row = b$row[i], nrows = b$nrows[i])
    v2 <- cbind(x, y = rep(y, each = length(x)), v2)
    v <- rbind(v, v2)
  }
  
  v <- v[complete.cases(v),]
  
  # write to output file
  #r_out <- writeValues(r_out, v, b$row[i])
}



# close files
r_out <- writeStop(r_out)
r_in <- readStop(r_in)



