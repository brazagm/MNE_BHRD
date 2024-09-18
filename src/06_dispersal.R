###
###  DISPERSAL SIMULATION
###
###  DESCRIPTION: estimate the dispersal parameters for simulations based on the
###     species' functional traits and simulate species dispersal from their
###     current geographic distribution to future scenarios using MigClim model. 
###     See MigClim userguide at
###     in https://github.com/robinengler/MigClim/blob/master/inst/doc/MigClim_userGuide.pdf
###
###  Created & Edited by: Alan Braz (@brazagm)
###
###  OBSERVATIONS: MigClim.migrate() fails when using .asc files exported with
###     writeRaster(). To solve that, I converted raster objects in data.frames.
###     Please, see MigClim userguide to more details on input formats.
###
###  ----------------------------- IMPORTANT!! --------------------------------- 
###     Importing MigClim functions is not working anymore.
###     Running R v3.6.3 in Ubuntu 20.04 and install MigClim from GitHub repo
###     was the only solution that worked :(
###
###  NEXT TASKS:
###     1. Replace strings with abbreviated species names by their complete names
###     (1) The datasheet "./results/Dispersal_dist.csv" stills need the
###         addition of species name manually following the code names/abbrev.
###

## older package versions must be acquired through "remotes::install_version("<NAME>", "<VERSION>")
## consult the package version for R 3.6.3 in https://cran-archive.r-project.org/bin/macosx/el-capitan/contrib/3.6/

# install SDMTools and MigClim from GitHub repositories
devtools::install_github("statsbomb/SDMTools") # install SDMTools
#devtools::install_github("robinengler/MigClim") # install MigClim
#install.packages("MigClim", "contriburl = https://cran.r-project.org/src/contrib/Archive/MigClim")
remotes::install_version("MigClim", "1.6.1")

# MigClim compilation from robinengler/ does not work in > R 4.3
# Alternatively, import MigClim functions from my forked repo
#source("https://raw.githubusercontent.com/brazagm/MigClim/master/R/MigClim.migrate.R")

# load packages

library(MigClim)
library(scales)
library(SDMTools)
library(terra)
library(tidyverse)

# set my working directory
#setwd("/home/alan/Dropbox/Alan/INMA/Dados/MNE_BHRD")


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories for model results
project_repo <- "/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD"

## inputs
#in_records <- "./processed_data/02-2_Gbif_records_clean_Revised.csv"
in_riverbasin <- "./data/shape/limite_BHRD.shp"
in_traits <- "./results/Dispersal_dist.csv"
in_models <- file.path(project_repo, "results", "niche_models")

## outputs
out_migclim <- "./results/migclim"
out_migclim_area <- file.path(out_migclim, "Migclim_area_results.csv")


# 1. DATA PREPARATION FOR MIGCLIM  ####-----------------------------------------
## 1.1. Import data  ####
## import Doce Riverbasin
doce <- vect(in_riverbasin) %>%
  project(., "+init=EPSG:4326") # reproject from SIRGAS2000 to WGS84

## import dispersal distance values
dist <- read.csv(in_traits) %>%
  mutate(Species = str_replace(Species, " ", "_"))

## species' name
species <- list.files(in_models) %>%
  grep(".csv|.png", ., value = TRUE, invert = TRUE)

## future projection names
time <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
scenarios <- c(paste0("MPI-ESM1-2-HR_", c("ssp126", "ssp370", "ssp585")))
futures <- apply(expand.grid(time, scenarios), 1, paste, collapse = "/")

## results dataframe
area <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(area) <- c("species", "scenario", "decolonized", "never_colonized", "stable", "colonized", "not_colonizable")

## 2.2. Import initial distribution and habitat suitability in future scenarios  ####
# Lafoensia_glyptocarpa distribution is outside the BHRD!!
for(sp in species){

  no_traits <- dist[which(is.na(dist$GF) | is.na(dist$log10MDD)), "Species"]
  no_traits <- sp %in% no_traits
  
  no_models <- list.files(file.path(project_repo, "results", "niche_models", sp))
  no_models <- length(no_models) == 1
  
  
  if(no_models == TRUE){
    print(paste(sp, "does not have niche models!!"))
  }else if(no_traits == TRUE){
    print(paste(sp, "does not have all neccessary traits!!"))
  }else{
  
  for(s in scenarios){
    
  ### 1.1.1. Initial distribution  ####
  # import raster from current distribution
  pattern <- list.files(file.path(project_repo, "results", "niche_models", sp, "current"), full.names = TRUE, recursive = TRUE) %>%
    grep("ensemble", ., value = TRUE)
  init <- rast(pattern)
  
  # initial distribution is given only by presence/absence
  init <- init[[2]]
  
  # crop distributions by the riverbasin
  init_crop <- crop(mask(init, doce), ext(doce))
  
  # convert raster to data.frame
  init_dt <- raster::as.data.frame(init_crop, xy = TRUE)
  colnames(init_dt) <- c("Xcoordinate", "Ycoordinate", "InitDist")
  init_dt <- na.omit(init_dt) # remove NA rows
  
  
  ## 1.1.2. Habitat suitability in future scenarios  ####
  # import
hsuit_dt2 <- map(.x = time,
      .f = function(x){
  
  pattern_future <- list.files(file.path(project_repo, "results", "niche_models", sp, x, s), full.names = TRUE) %>%
    grep("ensemble", ., value = TRUE)
  hsuit <- rast(pattern_future)
  
  # create filtered model
#  hsuit_filtered <- mask(hsuit[[1]], hsuit[[2]], maskvalue = 0, updatevalue = 0)
  hsuit_filtered <- mask(hsuit[[1]], hsuit[[2]], maskvalue = 0, updatevalue = NA)
  
  # crop distributions by the riverbasin
  hsuit_crop <- crop(mask(hsuit_filtered, doce), ext(doce))
  
  # rescale to max and 0
  range <- range(values(hsuit_crop), na.rm = TRUE)
  values(hsuit_crop) <- scales::rescale(values(hsuit_crop), to = c(0, max(values(hsuit_crop), na.rm = TRUE)), from = range)
  
  # replace NAs by zero
  hsuit_crop[is.na(hsuit_crop[])] <- 0 
  hsuit_crop <- mask(hsuit_crop, hsuit[[2]])
  
  # convert raster to data.frame
  hsuit_dt <- raster::as.data.frame(hsuit_crop, xy = FALSE)
  colnames(hsuit_dt) <- "HSmap"
  hsuit_dt <- na.omit(hsuit_dt) # remove NA rows
  hsuit_dt$HSmap <- as.data.frame(as.integer(hsuit_dt$HSmap*10)) 
  
  return(hsuit_dt)
  
      })
  
  hsuit_dt3 <- cbind(HSmap1 = hsuit_dt2[[1]]$HSmap,
        HSmap2 = hsuit_dt2[[2]]$HSmap,
        HSmap3 = hsuit_dt2[[3]]$HSmap,
        HSmap4 = hsuit_dt2[[4]]$HSmap)
  colnames(hsuit_dt3) <- c("HSmap1", "HSmap2", "HSmap3", "HSmap4")
  
  # check for any values < 0 or > 1000 and correct then
  hsuit_dt3 <- hsuit_dt3 %>%
    mutate(across(everything(), ~ifelse(. < 0, 0, ifelse(. > 1000, 1000, .))))
  

  ## 1.1.3 Import dispersal distance values
  dist_value <- dist[which(dist$Species == sp), "log10MDD"]
  dist_value <- 10^dist_value
  if(!sp %in% dist$Species){
  
    dist_kernel <- c(1, 0.5)
  
  }else if(dist_value <= 5*1000){
    
    dist_kernel <- dist_value/c(5*1000, 10*1000, 15*1000, 20*1000)
    dist_kernel <- dist_kernel[dist_kernel > 0.25 | dist_kernel == dist_kernel[1]]
    
  }else if(dist_value > 5*1000 & dist_value <= 10*1000){
    
    dist_kernel <- dist_value/c(5*1000, 10*1000, 15*1000, 20*1000)
    dist_kernel[dist_kernel > 1] <- 1
    
  }else if(dist_value > 10*1000){
    
    dist_kernel <- dist_value/c(5*1000, 10*1000, 15*1000, 20*1000, 25*1000)
    dist_kernel[dist_kernel > 1] <- 1
    
  }
  
  ## 1.1.3. Propagule potential
  # Define the logistic function
  logistic_function <- function(x, a, b, c, d) {
    return(a + b / (1 + exp(-c * (x - d))))
  }
  
  # Parameters to control the sigmoid curve
  #a <- 1  # Minimum value (lower asymptote)
  #b <- 9  # Maximum value (upper asymptote)
  #c <- 1   # Controls the steepness of the curve
  #d <- 5   # Controls the midpoint of the curve

  
  growth_form <- dist[which(dist$Species == sp), "GF"]
  
  if(growth_form == "shrub"){
    
    init_maturity <- 1
    
    x_values <- seq(1, 5)
    
    propag_prod <- logistic_function(x_values, a = 0, b = 1, c = 1, d = 5/2)
    
  }else if(growth_form == "tree"){
    
    init_maturity <- 11
    
    x_values <- seq(1, 30)
    
    propag_prod <- logistic_function(x_values, a = 0, b = 1, c = 1, d = 30/2)
    
  }
  
  
  ## 2.3. Test if current and future data have the same length
  if(nrow(init_dt) != nrow(hsuit_dt3)){
    
    print(paste0(sp, ": current and future distributions does not match"))
    
  } else {
    
    #______________________________
    # 2. DISPERSAL SIMULATION  ####
    #______________________________
    path_results <- file.path(out_migclim, sp, s)
    
    ifelse(dir.exists(path_results), paste0(sp, " directory already exists!"),
           dir.create(path_results, recursive = TRUE))
    
    # simulate migration
    MigClim.migrate(iniDist = init_dt,
                    hsMap = hsuit_dt3, 
                    rcThreshold = 0, # continuous hsuit values assumed if = 0
                    envChgSteps = 4, # number of environmental changes in suitability
                    dispSteps = 20, # number of dispersal steps within envChgSteps
                    dispKernel = dist_kernel, # prob function over distance
                    barrier = "", barrierType = "strong", 
                    iniMatAge = init_maturity, # initial maturity age
                    propaguleProd = propag_prod, # propagule production potential
                    lddFreq = 0.0, # prob of long distance dispersal events
                    lddMinDist = NULL, lddMaxDist = NULL,
                    simulName = paste0("MigClim_", sp), # base name
                    replicateNb = 1, # number of simulation replicates
                    overWrite = TRUE, testMode = FALSE,
                    fullOutput = TRUE, keepTempFiles = TRUE)
    
    # create comparison rasters
    final <- rast(paste0("./MigClim_", sp, "/MigClim_", sp, "_raster.asc"))
    
    # see MigClim tutorial for details on raster output values
    # basically, decolonized cells < 0; colonized > 1 and < 30.000; stable cells = 1; and never colonized = 0; potentially suitable but not colonized = 30.000
    values(final)[values(final) < 0] <- -1
    values(final)[values(final) > 1 & values(final) < 30000] <- 2
    writeRaster(final, file.path(path_results, paste0(sp, "_", s, ".tif")), overwrite = TRUE)
    
    
   # MigClim.plot(asciiFile = paste0("./MigClimTest_", sp, "/MigClimTest_", sp, "_raster.asc"),
  #               outDir = "", fileFormat = "inR", fullOutput = FALSE)

    decol <- final ; decol[decol != -1] <- NA
    cell_size <- terra::cellSize(decol, mask = TRUE)
    decol <- sum(values(cell_size), na.rm = TRUE)
    
    never_col <- final ; never_col[never_col != 0] <- NA
    cell_size <- terra::cellSize(never_col, mask = TRUE)
    never_col <- sum(values(cell_size), na.rm = TRUE)
    
    stable <- final ; stable[stable != 1] <- NA
    cell_size <- terra::cellSize(stable, mask = TRUE)
    stable <- sum(values(cell_size), na.rm = TRUE)
    
    col <- final ; col[col != 2] <- NA
    cell_size <- terra::cellSize(col, mask = TRUE)
    col <- sum(values(cell_size), na.rm = TRUE)
    
    no_col <- final ; no_col[no_col != 30000] <- NA
    cell_size <- terra::cellSize(no_col, mask = TRUE)
    no_col <- sum(values(cell_size), na.rm = TRUE)
    
    sp_area <- data.frame(species = sp,
                          scenario = s,
                          decolonized = decol,
                          never_colonized = never_col,
                          stable = stable,
                          colonized = col,
                          not_colonizable = no_col)
    
    area <- rbind(area, sp_area)
    
    #colnames(area) <- c("species", "scenario", "decolonized", "never colonized", "stable", "colonized")
    
    migclim_files <- list.files(file.path(paste0("MigClim_", sp)))
    migclim_files <- c(migclim_files, list.files("./", pattern = ".asc"))
    
    steps_files <- grep("_step_", migclim_files, value = TRUE)
    ifelse(dir.exists(file.path(path_results, "steps")), paste0(sp, " directory already exists!"),
           dir.create(file.path(path_results, "steps"), recursive = TRUE))
    file.copy(from = file.path(paste0("MigClim_", sp), steps_files),
    to = file.path(path_results, "steps", steps_files))
    
    initial_files <- grep("HSmap|InitialDist", migclim_files, value = TRUE)
    ifelse(dir.exists(file.path(path_results, "initial")), paste0(sp, " directory already exists!"),
           dir.create(file.path(path_results, "initial"), recursive = TRUE))
    file.copy(from = file.path(initial_files),
    to = file.path(path_results, "initial", initial_files))
 
    results_files <- grep("raster|.txt", migclim_files, value = TRUE)
    file.copy(from = file.path(paste0("MigClim_", sp), results_files),
    to = file.path(path_results, results_files))
    
    unlink(file.path(paste0("MigClim_", sp)), recursive = TRUE)
    unlink(initial_files)
      }

  }
}
}

## get initial distribution area from potential climatic range
enm_area <- read.csv("/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD/results/niche_models/niche_models_area_results.csv") %>%
  filter(time == "2081-2100") %>%
  select(species, scenario, initial)

area <- left_join(area, enm_area, by = c("species", "scenario"))
  
data <- area %>% 
  mutate_at(vars(3:6), as.numeric) %>%
  mutate(final = stable + colonized) %>%
  mutate(gain_loss = ((final - initial)/initial)*100)

#filtered <- data %>% 
#  mutate(initial = stable + decolonized, final = stable + colonized) %>%
#  mutate(gain_loss = ((final - initial)/initial)*100) %>%
  # Rank enumerates values by their rank (the lower the number , the lower the rank)
#  filter(
#    ( 10 >= rank(gain_loss)) |  # We want the low ranks  OR ("|")
#      (rank(gain_loss) > nrow(filtered)-10) # the high ranks
#  )

write_csv(data, out_migclim_area)


