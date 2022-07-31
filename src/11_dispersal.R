###
###  DISPERSAL SIMULATION
###
###  Description: simulate species dispersal from their current geographic 
###     distribution to future scenarios using MigClim model. See MigClim userguide
###     in https://github.com/robinengler/MigClim/blob/master/inst/doc/MigClim_userGuide.pdf
###
###  Created & Edited by: Alan Braz (@brazagm) // 31 May 2022
###
###  Observations: MigClim.migrate() fails when using .asc files exported with
###     writeRaster(). To solve that, I converted raster objects in data.frames.
###     Please, see MigClim userguide to more details on input formats.
###
###  Next tasks:
###

library(MigClim)
library(remotes)
library(tidyverse)

remotes::install_version("SDMTools", "1.1-221")
install_github("robinengler/MigClim")

#_____________________
# 1. IMPORT DATA  ####
#_____________________
## 1.1. Inputs and outputs  ####
input <- "./processed_data/Migclim_test"
input_bhrd <- "/home/alan/Documentos/azure_MNE_BHRD/azure/shapes/bhrd_wgs84_dissolv.shp"
output <- "./results"

# import Doce Riverbasin
doce <- shapefile("/home/alan/Documentos/azure_MNE_BHRD/azure/shapes/bhrd_wgs84_dissolv.shp")
proj4string(doce) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # set projection

# species name
spp_names <- list.files("/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo")


## 1.2. Data importation  ####
list <- list.files(input, full.names = TRUE, recursive = TRUE)

# import initial distribution
init <- stack(list[grep("present/ensemble/", list)])
init <- raster("/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo/Acnistus_arborescens/present/ensemble/Acnistus_arborescens_ensemble_0.5_consensus.tif")
# init[is.na(init[])] <- -9999 # change NA values by -9999

# import habitat suitability distribution in future scenarios
hsuit <- raster("/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo/Acnistus_arborescens/mpi_585/ensemble/Acnistus_arborescens_TSSmax_ensemble_weighted_average.tif")
# hsuit[is.na(hsuit[])] <- -9999 # change NA values by -9999


#__________________________
# 2. DATA PREPARATION  ####
#__________________________

azure_dir <- "/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo"

area <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(area) <- c("species", "decolonized", "never colonized", "stable", "colonized")

for(sp in spp_names){
  
  ## 2.1. Initial distribution  ####
  # import
  pattern <- file.path(azure_dir, sp, "present/ensemble", paste0(sp, "_ensemble_0.5_consensus.tif"))
  init <- raster(pattern)
  # init[is.na(init[])] <- -9999 # change NA values by -9999
  
  # crop distributions by the riverbasin
  init_crop <- crop(mask(init, doce), extent(doce))
  
  # convert raster to data.frame
  init_dt <- raster::as.data.frame(init_crop, xy = TRUE)
  colnames(init_dt) <- c("Xcoordinate", "Ycoordinate", "InitDist")
  init_dt <- na.omit(init_dt) # remove NA rows
  
  
  ## 2.2. Habitat suitability in future scenarios
  # import
  pattern <- file.path(azure_dir, sp, "mpi_585/ensemble/", paste0(sp, "_TSSmax_ensemble_weighted_average.tif"))
  hsuit <- raster(pattern)
  # hsuit[is.na(hsuit[])] <- -9999 # change NA values by -9999
  
  # crop distribution by riverbasin
  hsuit_crop <- crop(mask(hsuit, doce), extent(doce))
  
  # convert raster to data.frame
  hsuit_dt <- raster::as.data.frame(hsuit_crop, xy = TRUE)
  colnames(hsuit_dt) <- c("Xcoordinate", "Ycoordinate", "HSmap1")
  hsuit_dt <- na.omit(hsuit_dt) # remove NA rows
  hsuit_dt <- as.data.frame(as.integer(hsuit_dt$HSmap1*1000))
  
  
  ## 2.3. Test if current and future data have the same length
  if(nrow(init_dt) != nrow(hsuit_dt)){
    
    print(paste0(sp, ": current and future distributions does not match"))
    
  } else {
    
    #______________________________
    # 3. DISPERSAL SIMULATION  ####
    #______________________________
    # simulate migration
    MigClim.migrate(iniDist = init_dt,
                    hsMap = hsuit_dt, 
                    rcThreshold = 300, # continuous hsuit values assumed if = 0
                    envChgSteps = 1, # number of environmental changes in suitability
                    dispSteps = 1, # number of dispersal steps within envChgSteps
                    dispKernel = c(1.0, 0.5), # prob function over distance
                    barrier = "", barrierType = "strong", 
                    iniMatAge = 1, # initial maturity age
                    propaguleProd = c(1.0), # propagule production potential
                    lddFreq = 0.0, # prob of long distance dispersal events
                    lddMinDist = NULL, lddMaxDist = NULL,
                    simulName = paste0("MigClimTest_", sp), # base name
                    replicateNb = 1, # number of simulation replicates
                    overWrite = TRUE, testMode = FALSE,
                    fullOutput = FALSE, keepTempFiles = FALSE)
    
    # create comparison rasters
    final <- raster(paste0("./MigClimTest_", sp, "/MigClimTest_", sp, "_raster.asc"))
    
    # see MigClim tutorial for details on raster output values
    # basically, decolonized cells < 0; colonized > 1; stable cells = 1; and never colonized = 0
    values(final)[values(final) < 0] <- -1
    values(final)[values(final) > 1] <- 2
    writeRaster(final, paste0("./results/Migclim_results_", sp, ".asc"), format = "ascii", overwrite = TRUE)
    
    
   # MigClim.plot(asciiFile = paste0("./MigClimTest_", sp, "/MigClimTest_", sp, "_raster.asc"),
  #               outDir = "", fileFormat = "inR", fullOutput = FALSE)

    decol <- final ; decol[decol != -1] <- NA
    cell_size <- raster::area(decol, na.rm = TRUE, weights = FALSE)
    decol <- sum(values(cell_size), na.rm = TRUE)
    
    
    never_col <- final ; never_col[never_col != 0] <- NA
    cell_size <- raster::area(never_col, na.rm = TRUE, weights = FALSE)
    never_col <- sum(values(cell_size), na.rm = TRUE)
    
    stable <- final ; stable[stable != 1] <- NA
    cell_size <- raster::area(stable, na.rm = TRUE, weights = FALSE)
    stable <- sum(values(cell_size), na.rm = TRUE)
    
    col <- final ; col[col != 2] <- NA
    cell_size <- raster::area(col, na.rm = TRUE, weights = FALSE)
    col <- sum(values(cell_size), na.rm = TRUE)
    
    area <- rbind(area, c(sp, decol, never_col, stable, col))
    
      }

  
  
}

write.csv(area, "./results/Migclim_area_results.csv")






stack(init_crop, hsuit_crop)

writeRaster(init_crop, "./processed_data/initial_distr.asc", format = "ascii", overwrite = TRUE)
writeRaster(hsuit_crop*1000, "./processed_data/suitability1.asc", format = "ascii", overwrite = TRUE)



