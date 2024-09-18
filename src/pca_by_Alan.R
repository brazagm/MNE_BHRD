###
###  Title: ENVIRONMENTAL VARIABLES PCA
###
###  Description: PCA calculation and projection from environmental variables
###     using 10,000 random points to sample the environmental conditions in
###     South America. This was a way to avoid the highly consumption using the
###     entire raster cells for PCA calculation
###
###  Created & Edited by: Alan Braz (@brazagm) // 17 Nov 2022
###
###  Observations:
###
###  Next tasks:
###

library(dismo)
library(raster)
library(tidyverse)


# 0. INPUTS & OUTPUTS  #####----------------------------------------------------
## functions
source("./src/functions/function_pca.r")

## repositories
repo <- "/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/South_America"

## path to the environmental variables 
list <- list.files(repo, pattern = ".tif$", full.names = TRUE, recursive = TRUE)


# 1. ENVIRONMENTAL VARIABLES PREPARATION  ####----------------------------------
## list Chelsa variables in 1km in global extent
list <- list.files("/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/Global", pattern = ".tif$", full.names = TRUE, recursive = TRUE)

## subset by a pattern
#list <- grep("/MPI_ssp370/", list, value= TRUE)

## import variables
r <- stack(list)

## set extent
rr <- extent(raster("/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/South_America/1981-2010/_CHELSA_bio1_1981.2010_V.2.1.tif"))

## crop and mask
r <- crop(r, extent(rr))
r <- mask(r, rr)

## get variable names
names <- names(r)

## write variables as .tif raster
for(n in names){
  writeRaster(r[[n]], filename = paste0("/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/South_America/_", n, ".tif"), format = "GTiff", overwrite = TRUE)
}


# 2. IMPORT ENVIRONMENTAL VARIABLES  ####---------------------------------------
## variable names
var_names <- c("bio1", paste0("bio", 10:19), paste0("bio", 2:9),
               "gdd0", "gdd10", "gdd5", "gls", "gsp", "gst", "npp")

## scenarios
current <- "1981-2010"
time <- c("2011-2040", "2041-2070", "2071-2100")
gcm <- c("MPI_ssp126", "MPI_ssp370", "MPI_ssp585")

## present rasters
pres <- grep("1981-2010", list, value = TRUE)
r_pres <- stack(pres)
names(r_pres) <- var_names
r_pres <- aggregate(r_pres, fact = 3, fun = mean, na.rm = TRUE)

## future rasters
for(i in time){
  for(j in gcm){
    
    future <- grep(file.path(i, j), list, value = TRUE)
    r_future <- stack(future)
    names(r_future) <- var_names
    r_future <- aggregate(r_future, fact = 3, fun = mean, na.rm = TRUE)
    
    
    ### run via var_pca_proj function
    var_pca_proj(raster_pres = r_pres,
                 raster_proj = r_future,
                 cum_sum = .95,
                 graphics = TRUE,
                 graphic_names = paste(i, j, sep = "/"),
                 graphic_var_pres_names = names(r_pres),
                 prefix_pca_raster_pres = "",
                 prefix_pca_raster_proj = "",
                 path_output_pca_pres = file.path(getwd(), "processed_data/pca", current),
                 path_output_pca_proj = file.path(getwd(), "processed_data/pca", i, j))
    
  }
}

