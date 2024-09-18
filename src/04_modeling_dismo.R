###
###  ECOLOGICAL NICHE MODELLING USING DISMO/BIOMOD
###
###  Description: calibrate ecological niche models using dismo and biomod2
###     package.
###
###  Created & Edited by:
###
###  Observations:
###     (1) 'in_variables' and 'out_models' are located in a local repository
###         because to avoid saving these files in my Dropbox directory. Please,
###         remove the 'proj_repo' to inputs and outputs within the working
###         directory.
###     (2) Niche models are projected to future scenarios only in the Doce
###         riverbasin. Add a comment (#) in line 88 to project to the entire
###         South America.
###     (3) Calibration/training area is defined as a 500-km buffer from the
###         presence records, while pseudo-absences were defined with a minimum
###         distance of 100 km from presences. That is, pseudo-absences were
###         created with a min and max distance of 100 and 500 km from presences,
###         respectively, and background points within a 500 km buffer.
###     (4) Model evaluating is using pseudo-absences as absences to calculate
###         performance metrics.
###     *** BIOMOD is working using 'terra' but dismo still needs 'raster'
###
###  Next tasks:
###

library(adehabitatHR)
library(dismo)
library(biomod2)
library(raster)
library(pROC)
#library(rgeos)
library(doParallel)
library(furrr)
library(scales)
library(spThin)
library(terra)
library(tidyverse)

setwd("/home/alan/Dropbox/Alan/INMA/Dados/MNE_BHRD")

# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories
# obs: 'in_variables' and 'out_models' are located in a local repository because to avoid saving these files in my Dropbox directory. Please, remove the 'proj_repo' to inputs and outputs within the working directory.
project_repo <- "/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD"

## inputs
in_records <- "./processed_data/02-3_Gbif_records_clean_human-revised.csv"
in_variables <- file.path(project_repo, "processed_data/pca")
riverbasin <- "./data/shape/limite_BHRD.shp"

## outputs
out_models <- file.path(project_repo, "results/niche_models")

## memory and processing
#unix::rlimit_as(1e18)
#unix::rlimit_all()
#plan(multisession, workers = detectCores()-5)


# 1. SPECIES DATA AND ENVIRONMENTAL VARIABLES  ####-----------------------------
## 1.1. Species records  ####
## import species records
data <- read.csv(in_records)

## species names with '_'
data$species <- str_replace_all(data$species, " ", "_")

## split records by species in a list
data_list <- split(data, data$species)

## get species names
species <- names(data_list)


## 1.2. Environmental variables  ####
## set scenario names
current <- "1970-2000"

## future projection names
time <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
scenarios <- c(paste0("MPI-ESM1-2-HR_", c("ssp126", "ssp370", "ssp585")))
futures <- apply(expand.grid(time, scenarios), 1, paste, collapse = "/")

## variable names
var_names <- c("pc01", "pc02", "pc03", "pc04", "pc05", "pc06")

## import riverbasin shapefile and set the projection
doce <- vect(riverbasin)
#proj4string(doce) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## import predictor variables
var_pred_all <- list.files(file.path(in_variables, current),
                       pattern = "asc$",
                       full.names = TRUE) %>%
  rast()

## set raster names
names(var_pred_all) <- var_names

## import future projection variables
# paralelize with furrr
var_future <- furrr::future_map(.x = futures,
                                .f = function(x){
                                  # import raster for future variables
                                  var_future <- list.files(file.path(in_variables, x),
                                             pattern = "asc$",
                                             full.names = TRUE) %>%
                                    rast()
                    
                    # set raster names
                    names(var_future) <- var_names
                    
                    # mask and mask future rasters by the riverbasin
                    var_future <- terra::crop(terra::mask(var_future, doce), ext(doce))
                    
                    # assign the 'x' future scenario to an object
                    assign(paste0("var_", x), var_future)
                  })

## set raster names
names(var_future) <- futures


# 2. MODEL CALIBRATION  ####----------------------------------------------------
## create a folder for niche model results if it does not exist
ifelse(dir.exists(out_models), "Results directory already exists!",
       dir.create(out_models, recursive = TRUE))

## for each species
for(sp in species[179:186]){
  
  # set a path to save species' niche models
  model_path <- file.path(out_models, sp)
  
  # create a folder for niche model of each species if it does not exist
  ifelse(dir.exists(model_path), paste0(sp, " directory already exists!"),
         dir.create(model_path))
  
  # presence data
  pres <- data_list[[sp]] %>%
    subset(., .valid == TRUE) # select only valide presences

  # thin species' presences by a minimum distance of 10 km
  thin(loc.data = pres, lat.col = "decimalLatitude", long.col = "decimalLongitude",
       spec.col = "species", thin.par = 10, reps = 50,
       locs.thinned.list.return = FALSE,
       write.files = TRUE, max.files = 1, out.dir = model_path,
       out.base = paste0(sp, "_records"), write.log.file = FALSE, verbose = TRUE)
  
  # import thinned presences
  pres <- read.csv(file.path(model_path, paste0(sp, "_records_thin1.csv")))
  
  # read presence records as a spatial points dataframe
  pres <- cbind(1, pres)
  pres <- pres[,-2]
  colnames(pres) <- c(sp, "lon", "lat")
  
  #coordinates(pres) <- ~lon+lat
  pres <- vect(pres, geom = c("lon", "lat"))
  crs(pres) <- "+proj=longlat +datum=WGS84 +no_defs"

  
  ## 2.1. Create pseudo-absences and background points ####
  # Delimit calibration area as the Minimum Convex Polygon (MCP) with a 300 km buffer
  # Create a MCP including all records
  #train_area <- adehabitatHR::mcp(pres, percent = 100) %>%
  #  spTransform(., CRS("+init=epsg:32724")) %>% # Project MCP shapefile
  #  raster::buffer(., width = 300*1000, dissolve = TRUE) %>% # Create a 200 km buffer
  #  spTransform(., CRS("+proj=longlat +datum=WGS84 +no_defs")) # Back to WGS84 unprojected
  
  train_area <- terra::convHull(pres) %>%
    terra::buffer(., width = 300*1000)
  
  # create a calibration/training area with 500-km buffer from presences
  #train_area <- raster::buffer(pres, width = 500*1000, dissolve = TRUE)
  
  # crop predictor variables by the training area extent
  var_pred <- crop(var_pred_all, ext(train_area))
  
  # pseudo-absences were randomly generated with a 100-km minimum distance of presence records and within the calibration/training area
  # create a mask with a min and max buffer (100 and 500 km) from presences
  msk_pseudo <- terra::buffer(pres, width = 100*1000) %>%
    aggregate(.) %>%
    terra::symdif(., train_area)
#    rgeos::gSymdifference(., train_area) # set the maximum limit
  
  # create pseudo-absences outside the mask
  pseudo <- dismo::randomPoints(raster(terra::mask(var_pred[[1]], msk_pseudo, inverse = FALSE)), n = 1000, ext = extent(raster(var_pred[[1]])), extf = 1)
  colnames(pseudo) <- c("lon", "lat") # Set column names
  pseudo <- as.data.frame((pseudo))

  # Bind presences and pseudo-absences in the same dataframe
  #pres_pseudo <- rbind(records[,c("lon","lat")], pseudo)
  #pres_pseudo$Mincanus <- c(rep(1, nrow(records)), rep(0, nrow(pseudo)))

  # Import presences/pseudo-absences dataframe as spatial points data frame
  pseudo <- cbind(0, pseudo)
  colnames(pseudo) <- c(sp, "lon", "lat")

  # Read absences as a spatial points dataframe
  coordinates(pseudo) <- ~lon+lat
  crs(pseudo) <- "+proj=longlat +datum=WGS84 +no_defs"

  # Create background points within the calibration area
  background <- dismo::randomPoints(raster(terra::mask(var_pred[[1]], train_area)), n = 10000, ext = extent(raster(var_pred[[1]])), extf = 1)



  ## 2.2. Model calibration and projections  ####
  # Presence-only and presence-background methods (i.e. Mahalanobis, Domain, and Maxent) were calibrated in dismo package
  # Presence-absence (i.e. GLM and BRT) were calibrated in biomod2 package
  ## Automate modelling procedures, by calibrating replicates and combining them into an average model for each algorithm through:
  # 2.2.1. Data preparation
  # 2.2.2. Model calibration
  # 2.2.3. Model prediction for current scenario
  # 2.2.4. Model projection for future scenarios
  # 2.2.5. Model evaluation
  
  ## Group each record type by kfold = 10 repetitions
  pres_group <- kfold(pres, 10)
  pseudo_group <- kfold(pseudo, 10)
  back_group <- kfold(background, 10)
  
  rep_eval <- furrr::future_map(.x = 1:10,
                    .f = function(n){
  
                      require(dismo)
    print(paste("Initializing model calibration of", sp, "- partition", n))
  
    ## 2.2.1 Data preparation  ####
    # For presence records
    pres_train <- pres[pres_group != n, ]
    pres_test <- pres[pres_group == n, ]
  
    # For pseudo-absences
    pseudo_train <- pseudo[pseudo_group != n, ]
    pseudo_test <- pseudo[pseudo_group == n, ]
  
    # For background points
    back_train <- background[back_group != n, ]
    back_test <- background[back_group == n, ]
  
    ## Formating data for biomod2
    resp <- rbind(values(pres_train), select(pseudo_train@data, sp))
    xy <- rbind(crds(pres_train), data.frame(x = pseudo_train@coords[,"lon"], y = pseudo_train@coords[,"lat"]))
    
    # Create a directory for biomod outputs
    biomod_path <- file.path(model_path, "biomod", paste0("replicate_", n))
    
    # create a folder for niche model of each species if it does not exist
    ifelse(dir.exists(biomod_path), paste0(biomod_path, " directory already exists!"),
           dir.create(biomod_path, recursive = TRUE))
    
    data <- BIOMOD_FormatingData(resp.var = resp,
                                 expl.var = var_pred,
                                 resp.xy = xy,
                                 resp.name = sp,
                                 dir.name = biomod_path,
                                 na.rm = TRUE)
  
    options <- bm_ModelingOptions(data.type = "binary",
                                  strategy = "default")
    
    ## 2.2.2. Model calibration  ####
    # in dismo package
    bioclim <- bioclim(stack(var_pred), as(pres_train, "Spatial"))
    mahal <- mahal(stack(var_pred), as(pres_train, "Spatial"))
    domain <- domain(stack(var_pred), as(pres_train, "Spatial"))
    maxent <- maxent(stack(var_pred), as(pres_train, "Spatial"), a = back_train,
                     removeDuplicates = TRUE, args = "outputformat=logistic")
    
    
    # in biomod package
    out <- BIOMOD_Modeling(bm.format = data,
                           models = c("GLM","GBM","RF"),
                           bm.options = options,
                           CV.strategy = "random",
                           CV.nb.rep = 1,
                           CV.perc = 1,
                           CV.do.full.models = FALSE,
                           scale.models = TRUE,
                           modeling.id = "biomod")
  
    
    ## 2.2.3. Model prediction for current scenario  ####
    print("Projecting to current scenario")
    
    # in dismo package
    pred_bioclim <- predict(stack(var_pred), bioclim, progress = "")
    pred_mahal <- predict(stack(var_pred), mahal, progress = "")
    pred_domain <- predict(stack(var_pred), domain, progress = "")
    pred_maxent <- predict(stack(var_pred), maxent, progress = "")
    
    # in biomod package
    pred.biomod <- BIOMOD_Projection(bm.mod = out,
                                     new.env = var_pred,
                                     proj.name = "current",
                                     models.chosen = "all",
                                     build.clamping.mask = FALSE,
                                     compress = "xz",
                                     output.format = ".grd",
                                     do.stack = TRUE)
    
    # stack all models calibrated by biomod2
    pred <- stack(paste0(biomod_path, "/", out@sp.name, "/proj_current/proj_current_", out@sp.name, ".grd"))
    
    # stack all raster predictions
    pred <- stack(pred_bioclim, pred_mahal, pred_domain, pred_maxent, pred) %>%
      rast()
    
    # set raster names by model and replicates
    names(pred) <- paste(c("bioclim", "mahal", "domain", "maxent",
                           "GLM", "GBM", "RF"), n, sep = "_")
    
    # rescale and reduce value precision
    # get range values of each model
    range_values <- data.frame(model = names(pred), min = NA, max = NA)
    
    # rescale from 0 to 100 and get min/max for each model
    for(m in names(pred)){
      
      range <- range(values(pred[[m]]), na.rm = TRUE)
      range_values[which(range_values$model == m), c("min", "max")] <- range
      
      values(pred[[m]]) <- scales::rescale(values(pred[[m]]), to = c(0, 100), from = range)
      values(pred[[m]]) <- round(values(pred[[m]]), digits = 2)
      
    }

    # export asc file
    current_path <- file.path(model_path, "current", "replicates")
    
    ifelse(dir.exists(current_path), "Current models directory already exists!",
           dir.create(current_path, recursive = TRUE))
    
    terra::writeRaster(pred, filename = file.path(current_path, paste0(sp, "_current_", names(pred), ".tif")), 
#                format = "ascii", bylayer = TRUE, suffix = names(pred),
overwrite = TRUE)
    
    
    ## 2.2.4. Model projection for future scenarios  ####
    print("Projecting to future scenarios...")
    
    for(x in futures){
      print(paste0("... projecting to ", x))
      
      var_proj <- var_future[[x]]
      
      proj_bioclim <- predict(stack(var_proj), bioclim, progress = "")
      proj_mahal <- predict(stack(var_proj), mahal, progress = "")
      proj_domain <- predict(stack(var_proj), domain, progress = "")
      proj_maxent <- predict(stack(var_proj), maxent, progress = "")
      
      proj.biomod <- BIOMOD_Projection(bm.mod = out,
                                       new.env = var_proj,
                                       proj.name = "future",
                                       models.chosen = "all",
                                       build.clamping.mask = FALSE,
                                       compress = "xz",
                                       output.format = ".grd",
                                       do.stack = TRUE)
      
      # stack all models calibrated by biomod2
      proj <- stack(paste0(biomod_path, "/", out@sp.name, "/proj_future/proj_future_", out@sp.name, ".grd"))
      
      # stack all raster predictions
      proj <- stack(proj_bioclim, proj_mahal, proj_domain, proj_maxent, proj) %>%
        rast()
      
      names(proj) <- paste(c("bioclim", "mahal", "domain", "maxent",
                             "GLM", "GBM", "RF"), n, sep = "_")
      
      # rescale from 0 to 100 and get min/max for each model
      for(m in names(proj)){
        
        range <- range_values[which(range_values$model == m), c("min", "max")]
        
        values(proj[[m]]) <- scales::rescale(values(proj[[m]]), to = c(0, 100), from = c(range$min, range$max))
        values(proj[[m]]) <- round(values(proj[[m]]), digits = 2)
        
      }
      
      # export asc file
      future_path <- file.path(model_path, x, "replicates")
      
      ifelse(dir.exists(future_path), paste0(x, " directory already exists!"),
             dir.create(future_path, recursive = TRUE))
      
      # adjust the projection names
      proj_names <- paste(gsub("/", "_", x, fixed = TRUE), names(proj), sep = "_")
      
      terra::writeRaster(proj, filename = file.path(future_path, paste0(sp, "_", gsub("/", "_", x), "_", names(proj), ".tif")),
#                  format = "ascii", bylayer = TRUE, suffix = names(proj),
overwrite = TRUE)
    }
    
    
    ## 2.3. Model evaluation  ####
    print("Model evaluation")
    
    # remove #s if running loop as for instead of future_map
#    if(n == 1){
      df <- data.frame(model = names(pred), AUC = NA, pAUC = NA,
                       sens = NA, spec = NA, TSS = NA,
                       kappa = NA, spec_sens = NA,
                       no_omission = NA, prevalence = NA,
                       equal_sens_spec = NA, sensitivity = NA)
#      }else{
#        df <- rbind(df, data.frame(model = names(pred), AUC = NA, pAUC = NA,
#                                   sens = NA, spec = NA, TSS = NA,
#                                   kappa = NA, spec_sens = NA,
#                                   no_omission = NA, prevalence = NA,
#                                   equal_sens_spec = NA, sensitivity = NA))
#        }
    
    p <- terra::extract(pred, pres_test)
    a <- terra::extract(pred, vect(pseudo_test))
    p_train <- terra::extract(pred, pres_train)
    a_train <- terra::extract(pred, vect(pseudo_train))
    
    for(i in names(pred)){
      # extract AUC values and threshold values
      eval <- dismo::evaluate(p[,i], a[,i])
      df[which(df$model == i), "AUC"] <- eval@auc
      
      # get the threshold values
      thr <- threshold(eval)
      metrics <- colnames(thr)
      df[which(df$model == i), metrics] <- thr
      
      # Calculate partial AUC
      table <- rbind(cbind(pred = p[,i], obs = 1), cbind(pred = a[,i], obs = 0))
      curve <- roc(as.factor(table[,"obs"]), table[,"pred"])
      pauc <- auc(curve, partial.auc = c(min(curve$sensitivities[curve$sensitivities > 0]), 1), partial.auc.focus = "sens", partial.auc.correct = FALSE)
      df[which(df$model == i), "pAUC"] <- pauc
      
      # Get threshold value that maximizes TSS score
      values <- c(p_train[,i], a_train[,i])
      table <- data.frame(threshold = values, sens = NA, spec = NA, TSS = NA)
      
      for(x in unique(values)){
        binary <- pred[[i]] > x
        
        pres.values <- terra::extract(binary, pres_test)
        pres.values <- na.omit(pres.values)
        true.pres <- length(which(pres.values == 1))
        false.pres <- length(which(pres.values == 0))
        
        abse.values <- terra::extract(binary, vect(pseudo_test))
        abse.values <- na.omit(abse.values)
        true.abse <- length(which(abse.values == 0))
        false.abse <- length(which(abse.values == 1))
        
        sensibility <- true.pres/(true.pres + false.pres)
        specificity <- true.abse/(true.abse + false.abse)
        tss <- (sensibility + specificity) - 1
        
        table[which(table$threshold == x), "sens"] <- sensibility
        table[which(table$threshold == x), "spec"] <- specificity
        table[which(table$threshold == x), "TSS"] <- tss
        }
      
      table <- na.omit(table)
      
      df[which(df$model == i), c("sens", "spec", "TSS")] <- unique(table[which(table$TSS == max(table$TSS)), c("sens", "spec", "TSS")])
      
      
      ## Export binary models
      # select threshold method
      t <- "spec_sens"
      
      # for present
      list <- list.files(file.path(model_path, "current/replicates"),
                         pattern = paste0(i, ".tif$"), full.names = TRUE)

      m <- rast(list)
      
      terra::writeRaster(m > df[which(df$model == i), t],
                  filename = paste0(model_path, "/current/replicates/", sp, "_current_", names(m), "_", t, ".tif"),
                  #format = "ascii",
                  overwrite = TRUE)
      
      # for future
      for(x in futures){
        list <- list.files(file.path(model_path, x, "replicates"),
                           pattern = paste0(i, ".tif$"), full.names = TRUE)
      
        m <- rast(list)
        terra::writeRaster(m > df[which(df$model == i), t],
                    filename = paste0(model_path,"/", x, "/replicates/", sp, "_", gsub("/", "_", x), "_", gsub(".", "-", names(m), fixed = TRUE), "_", t, ".tif"), 
                    #format = "ascii",
                    overwrite = TRUE)
      }
    }
    return(df)
  }
  )
  # bind rows of dataframes in list
  rep_eval <- bind_rows(rep_eval)
  
  # export replicates evaluation results
  write_csv(rep_eval, paste0(model_path, "/", sp, "_replicates_evaluation.csv"))
  
  # remove the directory with biomod2 products
  unlink(file.path(model_path, "biomod"), recursive = TRUE)
}


# 3. FINAL MODELS AND ENSEMBLE MODELS  ####------------------------------------
# generate final models per algorithms and the ensemble model for each species
for(sp in species[75:186]){

  print(paste("Initializing model averaging of", sp))
  
  # set the path to species' niche models
  model_path <- file.path(out_models, sp)
  
  if(length(list.files(model_path)) == 1){
    
    print(paste(sp, "does not have niche models!!"))
    
  }else{
  
  
  # model names
  model_names <- c("bioclim", "mahal", "domain", "maxent", "GLM", "GBM", "RF")
  
  # import evaluation metrics for replicates
  df <- read.csv(file.path(model_path, paste0(sp, "_replicates_evaluation.csv")))
  
  # create a dataframe for final model performance metrics
  df_avg <- data.frame(model = c(model_names, "ensemble"), n = NA,
                       AUC_mean = NA, AUC_sd = NA,
                       pAUC_mean = NA, pAUC_sd = NA,
                       TSS_mean = NA, TSS_sd = NA)
  
  ## 3.1. Final models by averaged replicates  ####
  ## for each scenario...
  for(t in c("current", futures)){
    
    # list rasters for scenario 't'
    list <- list.files(file.path(model_path, t), pattern = ".tif$",
                       recursive = TRUE , full.names = TRUE)
    
    # for each model...
    for(m in model_names){
      
      print(paste("Averaging models of", m, "in", t))
      
      # get good models (TSS >= 0.5)
      good <- subset(df, TSS >= 0.5)
      good <- good[grep(m, good$model),]
      
      # create final model for model 'm' if it has at least one good replicate
      if(nrow(good) == 0){
        
        # give 0 for models without any good replicates
        df_avg[which(df_avg$model == m), "n"] <- 0
        print(paste0(m, " does not have any good model! :("))
        
      } else {
        
        # list of good replicates for model 'm'
        good <- good$model
        
        ### 3.1.1. Average of continuous replicates  ####
        cont <- grep("_spec_sens", list, value = TRUE, invert = TRUE) %>%
          grep(m, ., value = TRUE) %>%
          rast()
        
        final_cont <- mean(cont)
        
        ## 3.1.2. Average of binary replicates  ####
        bin <- grep("_spec_sens", list, value = TRUE, invert = FALSE) %>%
          grep(m, ., value = TRUE) %>%
          rast()
        
        final_bin <- sum(bin) > length(good)/2
        
        ## 3.1.3. Export final models  ####
        # adjust scenario name
        t_name <- gsub("/", "_", t, fixed = TRUE)
        
        # export average continuous model for the algorithm 'm'
        writeRaster(final_cont, paste0(file.path(model_path, t), "/", sp, "_", t_name, "_", m, "_avg.tif"), overwrite = TRUE)
        
        # export average binary model for the algorithm 'm'
        writeRaster(final_bin, paste0(file.path(model_path, t), "/", sp, "_", t_name, "_", m, "_spec_sens.tif"), overwrite = TRUE)
        
      }
      
      ## 3.1.4. Export final model performance  ####
      if(t == "current" & length(good) >= 1){
        
        df_avg[which(df_avg$model == m), "n"] <- length(good)
        df_avg[which(df_avg$model == m), "AUC_mean"] <- mean(df[which(df$model %in% good), "AUC"])
        df_avg[which(df_avg$model == m), "AUC_sd"] <- sd(df[which(df$model %in% good), "AUC"])
        df_avg[which(df_avg$model == m), "pAUC_mean"] <- mean(df[which(df$model %in% good), "pAUC"])
        df_avg[which(df_avg$model == m), "pAUC_sd"] <- sd(df[which(df$model %in% good), "pAUC"])
        df_avg[which(df_avg$model == m), "TSS_mean"] <- mean(df[which(df$model %in% good), "TSS"])
        df_avg[which(df_avg$model == m), "TSS_sd"] <- sd(df[which(df$model %in% good), "TSS"])
        
      }
    }
    
    ## 3.2. Ensemble model by averaged final models  ####
    # list final models for scenario 't'
    list_final <- list.files(file.path(model_path, t), pattern = ".tif$",
                             recursive = FALSE , full.names = TRUE)
    
    print(paste("Ensemble model for", t))
    
    ### 3.2.1. Ensemble model with continuous values  ####
    cont_final <- grep("_spec_sens", list_final, value = TRUE, invert = TRUE) %>%
      rast()
    
    ensemble_cont <- mean(cont_final)
    
    ## 3.1.2. Average of binary replicates  ####
    bin_final <- grep("_spec_sens", list_final, value = TRUE, invert = FALSE) %>%
      rast()
    
    ensemble_bin <- sum(bin_final) > length(names(bin_final))/2
    
    ## 3.1.3. Export ensemble model  ####
    # adjust scenario name
    t_name <- gsub("/", "_", t, fixed = TRUE)
    
    # export average continuous model for the algorithm 'm'
    writeRaster(ensemble_cont, paste0(file.path(model_path, t), "/", sp, "_", t_name, "_ensemble_avg.tif"), overwrite = TRUE)
    
    # export average binary model for the algorithm 'm'
    writeRaster(final_bin, paste0(file.path(model_path, t), "/", sp, "_", t_name, "_ensemble_spec_sens.tif"), overwrite = TRUE)
    
    df_avg[which(df_avg$model == "ensemble"), "n"] <- length(names(bin_final))
    
    df_avg[which(df_avg$model == "ensemble"), "AUC_sd"] <- sd(df_avg$AUC_mean, na.rm = TRUE)
    df_avg[which(df_avg$model == "ensemble"), "AUC_mean"] <- mean(df_avg$AUC_mean, na.rm = TRUE)

    df_avg[which(df_avg$model == "ensemble"), "pAUC_sd"] <- sd(df_avg$pAUC_mean, na.rm = TRUE)
    df_avg[which(df_avg$model == "ensemble"), "pAUC_mean"] <- mean(df_avg$pAUC_mean, na.rm = TRUE)
    
    df_avg[which(df_avg$model == "ensemble"), "TSS_sd"] <- sd(df_avg$TSS_mean, na.rm = TRUE)
    df_avg[which(df_avg$model == "ensemble"), "TSS_mean"] <- mean(df_avg$TSS_mean, na.rm = TRUE)

    
  }

  write_csv(df_avg, paste0(model_path, "/", sp, "_avg_models_evaluation.csv"))
  
  }
}


# 4. AVERAGE MODEL PERFORMANCE  ####--------------------------------------------
names <- list.files(in_models) %>%
  grep(".png|.csv", ., value = TRUE, invert = TRUE)

## import results with model performances
## use only pAUC and TSS mean values among all species
eval <- list.files(in_models, pattern = "avg_models_evaluation.csv", recursive = TRUE, full.names = TRUE) %>%
  grep(paste(names, collapse = "|"), ., value = TRUE) %>% # select only species with migclim results
  read_csv() %>%
  filter(model == "ensemble") %>%
  group_by(model) %>%
  summarise(n_mean = mean(n),
            AUC_avg = mean(AUC_mean), AUC_sd = sd(AUC_mean),
            pAUC_avg = mean(pAUC_mean), pAUC_sd = sd(pAUC_mean),
            TSS_avg = mean(TSS_mean), TSS_sd = sd(TSS_mean))

## get mean value of pAUC and TSS values
write_csv(eval, file.path(out_models, "niche_models_average_performance.csv"))



### 4. FIGURES AND AREA GAIN/LOSS  ####-----------------------------------------
# import riverbasin limits
doce <- vect(riverbasin) %>%
  project("+init=EPSG:4326")

# list enm (only ensemble models) file paths
model_files <- list.files(out_models, pattern = "ensemble_spec_sens.tif", full.names = TRUE, recursive = TRUE)


## 4.1. Area gain/loss in future scenarios  ####
# create a dataframe for the results
area <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(area) <- c("species", "time", "scenario", "gain", "loss", "stable")

# calculate distribution area for each species
for(sp in species[76:186]){
  
  # select only enms of the species
  models <- grep(sp, model_files, value = TRUE)
  
  if(length(models) == 0){
    
    print(paste(sp, "does not have niche models!!"))
    
  }else{
  
  # import initial distribution in current
  init <- grep("current", models, value = TRUE) %>%
    rast() %>%
    crop(., ext(doce)) %>% # crop and mask by the riverbasin limits
    mask(., doce) 
  
  # final distribution in different scenarios 
  for(s in scenarios){
    
    for(t in time){
    
      # import final distribution in 's' scenario and 't' time
    final <- grep(paste(t, s, sep = "_"), models, value = TRUE) %>%
      rast() %>%
      crop(., ext(doce)) %>%
      mask(., doce)
    
    # sum initial distribution (= -1) and final distribution (= +1)
    final[final == 1] <- 2
    
    # for some species, their known geographic range is smaller than the riverbasin
    # in these cases, NA values within the basin is converted to zeros
    if(ext(init) != ext(final)){
      init <- extend(init, doce) %>%
        replace(., is.na(.), 0) %>%
        mask(., final)
    }
    
    sum <- init + final
    
    # calculate gain area (= 2)
    gain <- sum ; gain[gain != 2] <- NA
    cell_size <- terra::cellSize(gain, mask = TRUE)
    gain <- sum(values(cell_size), na.rm = TRUE)
    
    # calculate loss area (= 1)
    loss <- sum ; loss[loss != 1] <- NA
    cell_size <- terra::cellSize(loss, mask = TRUE)
    loss <- sum(values(cell_size), na.rm = TRUE)
    
    # calculate stable area (= 3)
    stable <- sum ; stable[stable != 3] <- NA
    cell_size <- terra::cellSize(stable, mask = TRUE)
    stable <- sum(values(cell_size), na.rm = TRUE)
    
    # resume as a data.frame
    sp_area <- data.frame(species = sp,
                          time = t,
                          scenario = s,
                          gain = gain,
                          loss = loss,
                          stable = stable)
    
    # export to data.frame output
    area <- rbind(area, sp_area)
    
    }
  }
  }
}

# summarize and calculate initial and final distribution areas
data <- area %>% 
  mutate_at(vars(4:6), as.numeric) %>%
  mutate(initial = stable + loss, final = stable + gain) %>%
  mutate(gain_loss = ((final - initial)/initial)*100)

# export
write_csv(data, file.path(out_models, "niche_models_area_results.csv"))

