###-----------------------------------------------------------------------------
###  Title: ENVIRONMENTAL VARIABLES
###
###  Description: Import environmental variables from WorldClim, crop them by
###    the Latin America extent, and summarize them in PCA axes.
###
###  Obs:
###-----------------------------------------------------------------------------

## load required libraries
library(dismo)
library(terra)
library(tidyverse)


# 0. INPUTS & OUTPUTS  #####----------------------------------------------------
## import required functions
source("./src/functions/function_pca.r")

## repositories
repo <- "/home/alan/Documentos/Repositório/GIS/layers/WorldClim_v2.1/2.5min_Global"

## inputs
in_records <- "./processed_data/02-3_Gbif_records_clean_human-revised.csv"

## outputs
out_pca <- "./processed_data/pca"


# 1. IMPORT AND PREPARE ENVIRONMENTAL VARIABLES  ####---------------------------
## list WorldClim variables in 5km in global extent
list <- list.files(repo, pattern = ".tif$", full.names = TRUE, recursive = TRUE)

## import occurrence records
records <- read.csv(in_records) %>%
  subset(., .valid == TRUE) %>%
  vect(., geom = c("decimalLongitude", "decimalLatitude"))

## set Latin America extent
extent <- ext(records)

## variable names
var_names <- c(paste0("bio_", c(1:19)))

## scenarios
current <- "1970-2000"
time <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
gcm <- c(paste0("MPI-ESM1-2-HR_", c("ssp126", "ssp370", "ssp585")))

## present rasters
pres <- grep(current, list, value = TRUE)
r_pres <- rast(pres) %>%
  crop(., ext(extent)) # crop by extent
names(r_pres) <- sort(var_names)
#r_pres <- aggregate(r_pres, fact = 3, fun = mean, na.rm = TRUE)

ifelse(dir.exists(file.path(out_pca, current)), "Results directory already exists!",
       dir.create(file.path(out_pca, current), recursive = TRUE))

## import variables
#r <- rast(list) %>%
#  crop(., ext(extent)) # crop by extent

## get variable names
#names <- names(r)

## write variables as .tif raster
#for(n in names){
#  writeRaster(r[[n]], filename = paste0("/home/alan/Documentos/Repositório/GIS/layers/Chelsa_1km_v2/South_America/_", n, ".tif"), format = "GTiff", overwrite = TRUE)
#}


# 2. PCA OF ENVIRONMENTAL VARIABLES  ####---------------------------------------
## run PCA analyses to extract environmental variables summary
## using 'function_pca.R'
## PCA runs for current rasters and projected for each future scenarios
for(i in time){
  for(j in gcm){
    
    future <- grep(file.path(i, j), list, value = TRUE)
    r_future <- rast(future) %>%
      crop(., ext(extent)) # crop by extent
    names(r_future) <- var_names
    r_future <- r_future[[sort(var_names)]]
    #r_future <- aggregate(r_future, fact = 3, fun = mean, na.rm = TRUE)
    
    ifelse(dir.exists(file.path(out_pca, i, j)), "Results directory already exists!",
           dir.create(file.path(out_pca, i, j), recursive = TRUE))
    
    ### run via var_pca_proj function
    var_pca_proj(raster_pres = r_pres,
                 raster_proj = r_future,
                 cum_sum = .95,
                 graphics = TRUE,
                 graphic_names = file.path(out_pca, paste(i, j, sep = "/")),
                 graphic_var_pres_names = names(r_pres),
                 prefix_pca_raster_pres = "",
                 prefix_pca_raster_proj = "",
                 path_output_pca_pres = file.path(getwd(), "processed_data/pca", current),
                 path_output_pca_proj = file.path(getwd(), "processed_data/pca", i, j))
    
  }
}
