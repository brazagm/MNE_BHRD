###
###  Title: POST-PROCESSING
###
###  Description: code for calculate area gain and loss of species geographic
###    ranges under different future scenarios, export png maps and raster 
###    files.
###
###  Created by: Alan Braz (@brazagm) & Danielle Moreira (@daniomoreira)
###
###  Observations:
###   1. this code is based and only works with ensemble_0.5_consensus models 
###      from ModleR
###   2. this code version includes only present, mpi_126 and mpi_585 as 
###      projections
###
###  Next steps:
###   1. include legend in exported png (maps)
###   2. find a better color palette for maps...
###

library(adehabitatHR)
library(raster)
library(rgeos)
library(rgdal)
library(tidyverse)

# 0. INPUTS AND OUTPUTS  ####
## 0.1. Input  ####
# get the output directory from ModleR
path <- "/home/alan/Documentos/Repositório/MNE_BHRD/niche_models" 

# get species name and projections from local folders
species <- list.files(path) # './modleR' must contain only species folders!
time <- c("2011-2040", "2041-2070", "2071-2100")
scenarios <- c("MPI_ssp126", "MPI_ssp585")
futures <- apply(expand.grid(time, scenarios), 1, paste, collapse = "/")
proj <- c("current", futures)

# get occurrence records
#records <- read.csv("./species/6_Gualaxo_ModleR.csv")

# get shapefiles for the Doce riverbasin and world borders
riverbasin <- "./data/shape/bhrd_wgs84_dissolv.shp"
#border <- "./shapes/ne_10m_admin_0_countries.shp"

# get file path for final models and ensemble of each species
input <- list()
for(i in species){
  list <- list()
  
  for(j in proj){
    list <- append(list, list(file.path(path, i, j, paste0(i, "_maxent_bin.asc"))))
  }
  
  names(list) <- proj
  input <- append(input, list(list))
}

# set list names by species
names(input) <- species


## 0.2. Output  ####
output <- "./results/distribution_area.csv"


# 1. ESTIMATE AREA GAIN AND LOSSES  ####
results <- data.frame(matrix(ncol = 4, nrow = 0))

# calculate range area for each species, in current and future scenarios
# note that this code version includes only c("present", "mpi_126", "mpi_585") as projections
for(i in species){
  
  # raster manipulation
#  raster_present <- list.files(input[[i]]$current, pattern = "0.5_consensus.tif", full.names = TRUE) %>%
#    stack(.) # raster files for present
#  raster_mpi_126 <- list.files(input[[i]]$mpi_126, pattern = "0.5_consensus.tif", full.names = TRUE) %>%
#    stack(.) # raster files for mpi_126
#  raster_mpi_585 <- list.files(input[[i]]$mpi_585, pattern = "0.5_consensus.tif", full.names = TRUE) %>%
#    stack(.) # raster files for mpi_585
 
  raster_present <- raster(input[[i]]$current)
  raster_mpi_126 <- raster(input[[i]]$`2071-2100/MPI_ssp126`)
  raster_mpi_585 <- raster(input[[i]]$`2071-2100/MPI_ssp585`)
  raster_present <- mask(crop(raster_present, extent(raster_mpi_126)), raster_mpi_126)
  
  ## 1.1. Calculate area gain/loss at the South America extent  ####
  # get sizes of all cells in raster (km2)
  present <- raster_present; present[present == 0] <- NA
  cell_size <- raster::area(present, na.rm = TRUE, weights = FALSE)
  present_area <- sum(values(cell_size), na.rm = TRUE)
  
  future_126 <- raster_mpi_126 ; future_126[future_126 == 0] <- NA
  cell_size <- raster::area(future_126, na.rm = TRUE, weights = FALSE)
  mpi_126_area <- sum(values(cell_size), na.rm = TRUE)
  
  future_585 <- raster_mpi_585 ; future_585[future_585 == 0] <- NA
  cell_size <- raster::area(future_585, na.rm = TRUE, weights = FALSE)
  mpi_585_area <- sum(values(cell_size), na.rm = TRUE)
  
  raster_sa <- stack(present, future_126, future_585)
  names(raster_sa) <- c("present", "mpi_126", "mpi_585")
  
  # compute the results
  results <- rbind(results, c(i, "ssp126", present_area, mpi_126_area))
  results <- rbind(results, c(i, "ssp585", present_area, mpi_585_area))
  
  ## 1.2. Calculate area gain/loss assuming no dispersal  ####
#  occ <- subset(records, sp == i, select = c("lon", "lat"))
  
  # delimit geographic range as the Minimum Convex Polygon (MCP) with a 100 km buffer
  # create a MCP including all records
#  msk <- adehabitatHR::mcp(SpatialPoints(occ, CRS("+proj=longlat +datum=WGS84 +no_defs")), percent = 100)
#  msk <- spTransform(msk, CRS("+init=epsg:32724")) # project MCP shapefile
#  msk <- raster::buffer(msk, width = 100*1000, dissolve = TRUE) # create a 100 km buffer
#  msk <- spTransform(msk, CRS("+proj=longlat +datum=WGS84 +no_defs")) # back to WGS84 unprojected
  
  # get sizes of all cells in raster (km2)
#  present <- crop(mask(raster_present, msk), extent(msk))
#  present[present == 0] <- NA
#  cell_size <- raster::area(present, na.rm = TRUE, weights = FALSE)
#  present_area <- sum(values(cell_size), na.rm = TRUE)
  
#  future_126 <- crop(mask(raster_mpi_126, msk), extent(msk))
#  future_126[future_126 == 0] <- NA
#  cell_size <- raster::area(future_126, na.rm = TRUE, weights = FALSE)
#  mpi_126_area <- sum(values(cell_size), na.rm = TRUE)
  
#  future_585 <- crop(mask(raster_mpi_585, msk), extent(msk))
#  future_585[future_585 == 0] <- NA
#  cell_size <- raster::area(future_585, na.rm = TRUE, weights = FALSE)
#  mpi_585_area <- sum(values(cell_size), na.rm = TRUE)
  
#  raster_no <- stack(present, future_126, crop(future_585, extent(present)))
#  names(raster_no) <- c("present", "mpi_126", "mpi_585")
  
  # compute the results
#  results <- rbind(results, c(i, "no_dispersal", present_area, mpi_126_area, mpi_585_area))
  
  
  ## 1.3. Calculate area gain/loss for the Doce Riverbasin  ####
  # import riverbasin shapefile
  doce <- shapefile(riverbasin)
#  proj4string(doce) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # get sizes of all cells in raster (km2)
#  present <- crop(mask(raster_present, doce), extent(doce))
#  present[present == 0] <- NA
#  cell_size <- raster::area(present, na.rm = TRUE, weights = FALSE)
#  present_area <- sum(values(cell_size), na.rm = TRUE)
  
#  future_126 <- crop(mask(raster_mpi_126, doce), extent(doce))
#  future_126[future_126 == 0] <- NA
#  cell_size <- raster::area(future_126, na.rm = TRUE, weights = FALSE)
#  mpi_126_area <- sum(values(cell_size), na.rm = TRUE)
  
#  future_585 <- crop(mask(raster_mpi_585, doce), extent(doce))
#  future_585[future_585 == 0] <- NA
#  cell_size <- raster::area(future_585, na.rm = TRUE, weights = FALSE)
#  mpi_585_area <- sum(values(cell_size), na.rm = TRUE)
  
#  raster_bhrd <- stack(present, future_126, crop(future_585, extent(present)))
#  names(raster_bhrd) <- c("present", "mpi_126", "mpi_585")
  
  # compute the results
#  results <- rbind(results, c(i, "BHRD", present_area, mpi_126_area, mpi_585_area))
  
  
  ## 1.4. Map gain/losses in geographic distributions  ####
  # set layout for the image exportation
  png(paste("./results/", i, ".png", sep = ""), units = "in", width = 10, height = 4, res = 300)
  
  # set figure layout
  par(mfrow = c(1, 2),
      oma = c(0, 1.5, 1.5, 0), mar = c(0.5, 0.5, 0.5, 0), mgp = c(1.8, 0.5, 0),
      cex.lab = 1.4, family = "sans", font.lab = 2)
  
  # color palette
  # blue #0568BF = area gain // yellow #ED553B = area stability // red #F5EA54 = area loss
  # suggestion: c("#0568BF", "#ED553B", "#F5EA54")
  palette <- c("purple", "red", "yellow")
  
  # import world borders shapefiles
#  world <- shapefile(border)
  
  # plot map for each extent (South America, no dispersal and BHRD)
  for(j in "raster_sa"){
    
    x <- get(j)
    x[["present"]][x[["present"]] == 1] <- 2
    
    # map area gain, loss and stability for both scenarios
    mpi_126 <- sum(x[["present"]], x[["mpi_126"]], na.rm = TRUE)
    mpi_126[mpi_126 == 0] <- NA
    
    mpi_585 <- sum(x[["present"]], x[["mpi_585"]], na.rm = TRUE)
    mpi_585[mpi_585 == 0] <- NA
    
    # export raster files
    writeRaster(mpi_126, paste("./results/", i, "_mpi126_", j, ".tif", sep = ""), 
                format = "GTiff", overwrite = TRUE)
    writeRaster(mpi_585, paste("./results/", i, "_mpi585_", j, ".tif", sep = ""), 
                format = "GTiff", overwrite = TRUE)
    
    # plot maps
    if(j == "raster_sa"){
      plot(mpi_126, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
      plot(doce, add = TRUE)
      
      mtext("MPI 126", side = 3, line = 0.5, cex = 1, font = 2)
#      mtext("South Am", side = 2, line = 0.5, cex = 1, font = 2)
      
      plot(mpi_585, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
      plot(doce, add = TRUE)
      mtext("MPI 585", side = 3, line = 0.5, cex = 1, font = 2)
      
#    } else if(j == "raster_no"){
#      plot(mpi_126, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
#      plot(world, add = TRUE)
      
#      mtext("no dispersal", side = 2, line = 0.5, cex = 1, font = 2)
      
#      plot(mpi_585, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
#      plot(world, add = TRUE)
      
#    } else if(j == "raster_bhrd"){
#      plot(mpi_126, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
#      plot(doce, add = TRUE)
      
#      mtext("BHRD", side = 2, line = 0.5, cex = 1, font = 2)
      
#      plot(mpi_585, yaxt = "n", xaxt = "n", col = palette, ylab = NA, legend = FALSE)
#      plot(doce, add = TRUE)
      
    }
  }
  dev.off()
}

# give column names for the output dataframe
colnames(results) <- c("species", "scenario", "present_area", "future_area")

# correct column class
results$species <- as.factor(results$species)
results$scenario <- as.factor(results$scenario)
results$present_area <- as.numeric(results$present_area)
results$future_area <- as.numeric(results$future_area)

# calculate ratio of area gain/loss
results$area_ratio <- results$future_area/results$present_area
results$area_percent <- ((results$future_area - results$present_area)/results$present_area)*100

# export results
write.csv(results, output)
