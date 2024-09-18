
library(dplyr)
library(raster)
library(rgeos)
library(tidyverse)


#_____________________
# 1. IMPORT DATA  ####
#_____________________
## 1.1. Inputs and outputs  ####
input_niche <- "./results/niche_models"
input_migclim <- "./results/migclim"
output <- "./results/migclim/suitable_colonized"

# species name
spp_names <- list.files("./results/niche_models/")


#__________________________
# 2. DATA PREPARATION  ####
#__________________________

scenarios <- c("MPI_ssp126", "MPI_ssp585")

area <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(area) <- c("species", "scenario", "colonized", "colonized_unsuitable", "suitable_inaccesible", "absence")

for(sp in spp_names){
  
  for(s in scenarios){

clim <- paste0("./results/niche_models/", sp, "/2071-2100/", s, "/", sp, "_maxent_bin.asc") %>%
  raster()

migclim <- paste0("./results/migclim/Migclim_results_", sp, "_", s, ".asc") %>%
  raster() %>%
  crop(., extent(clim))
migclim <- migclim > 0

clim <- mask(clim, migclim)

# 0 = absence; 1 = climatic suitability; 10 = colonized but without suitability; 11 = colonized
result <- clim + (migclim*10)

col <- result ; col[col != 11] <- NA
cell_size <- raster::area(col, na.rm = TRUE, weights = FALSE)
col <- sum(values(cell_size), na.rm = TRUE)

col_unsuit <- result ; col_unsuit[col_unsuit != 10] <- NA
cell_size <- raster::area(col_unsuit, na.rm = TRUE, weights = FALSE)
col_unsuit <- sum(values(cell_size), na.rm = TRUE)

suit <- result ; suit[suit != 1] <- NA
cell_size <- raster::area(suit, na.rm = TRUE, weights = FALSE)
suit <- sum(values(cell_size), na.rm = TRUE)

absence <- result ; absence[absence != 0] <- NA
cell_size <- raster::area(absence, na.rm = TRUE, weights = FALSE)
absence <- sum(values(absence), na.rm = TRUE)

area <- rbind(area, c(sp, s, col, col_unsuit, suit, absence))
colnames(area) <- c("species", "scenario", "colonized", "colonized_unsuitable", "suitable_inaccesible", "absence")

  }
}



