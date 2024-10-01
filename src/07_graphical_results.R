###
###  GRAPHICAL RESULTS
###
###  DESCRIPTION: 
###
###  Created & Edited by: Alan Braz (@brazagm)
###
###  OBSERVATIONS:
###
###  NEXT TASKS:
###     (1) Revise section 4.
###

# load packages
library(ggpattern)
library(terra)
library(tidyterra)
library(tidyverse)


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories for model results
project_repo <- "/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD"

## inputs
in_models <- file.path(project_repo, "results", "niche_models")
in_migclim <- file.path(project_repo, "results", "migclim")
in_migclim_area <- file.path(in_migclim, "Migclim_area_results.csv")
#in_records <- "./processed_data/02-2_Gbif_records_clean_Revised.csv"
in_riverbasin <- "./data/shape/limite_BHRD.shp"
in_traits <- "./results/Dispersal_dist.csv"
in_munic <- "/home/alan/Documentos/Repositório/GIS/shapes/borders/municipios.shp"


## outputs
out_performance <- "./results/figures/01_enm_performance.png"
out_gainloss <- "./results/figures/02_gain_loss_area_results.png"
out_restoration <- "./results/figures/03_restoration_priority.png"


# 1. DATA PREPARATION  ####-----------------------------------------------------
## 1.1. Import data  ####
## import Doce Riverbasin
doce <- vect(in_riverbasin) %>%
  project(., "+init=EPSG:4326") # reproject from SIRGAS2000 to WGS84

## import dispersal distance values
dist <- read.csv(in_traits) %>%
  mutate(Species = str_replace(Species, " ", "_"))

## 1.2. Species names  ####
## species' name
#species <- list.files(in_models) %>%
#  grep(".csv|.png", ., value = TRUE, invert = TRUE)

## species' name
#species <- list.files(in_models) %>%
#  grep(".csv|.png", ., value = TRUE, invert = TRUE)
omit_spp <- c("Enterolobium_timbouva", "Piptocarpha_angustifolia") # omit because they do not occur within BHRD

spp_names <- read.csv(in_migclim_area) %>%
  filter(!species %in% omit_spp) %>%
  pull(species) %>%
  unique()


## 1.3. Future projection names  ####
time <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
scenarios <- c(paste0("MPI-ESM1-2-HR_", c("ssp126", "ssp370", "ssp585")))
futures <- apply(expand.grid(time, scenarios), 1, paste, collapse = "/")


# 2. MODEL PERFORMANCES  ####---------------------------------------------------
## 2.1. Import results with model performances  ####
## use only pAUC and TSS mean values among all species
eval <- list.files(in_models, pattern = "avg_models_evaluation.csv", recursive = TRUE, full.names = TRUE) %>%
  grep(paste(spp_names, collapse = "|"), ., value = TRUE) %>% # select only species with migclim results
  read_csv() %>%
#  filter(model == "ensemble") %>%
  pivot_longer(c("pAUC_mean", "TSS_mean"), names_to = "metric", values_to = "values") %>%
  select(model, n, metric, values) %>%
  mutate(metric = str_replace(metric, "_mean", ""))

## get the type of algorithm
eval <- eval %>%
  mutate(class = if_else(model %in% c("bioclim", "domain", "mahal"), "presence-only", if_else(model %in% c("GLM", "GBM", "RF"), "presence-absence", if_else(model == "maxent", "presence-background", "ensemble")))) %>%
  mutate(model = str_replace(model, "bioclim", "Bioclim"),
         model = str_replace(model, "domain", "Domain"),
         model = str_replace(model, "mahal", "Mahalanobis"),
         model = str_replace(model, "maxent", "Maxent"),
         model = str_replace(model, "ensemble", "Ensemble"))

## 2.2. Calculate mean value of pAUC and TSS values  ####
## mean, sd and max/min values fo each algorithm
mean <- eval %>%
  group_by(model, metric) %>%
  summarise(n = length(na.omit(values)), mean = mean(na.omit(values)), sd = sd(na.omit(values)), max = max(na.omit(values)), min = min(na.omit(values)))

## 
#summary <- mean %>%
#  group_by(metric, model) %>%
#  summarise(n = length(mean))


## 2.3. Boxplot the performance distribution for each algorithm  ####
## reorder the groups
eval$model <- factor(eval$model , levels = c("Bioclim", "Domain", "Mahalanobis", "GLM", "GBM", "RF", "Maxent", "Ensemble"))
eval$class <- factor(eval$class , levels = c("presence-only", "presence-absence", "presence-background", "ensemble"))

## plot boxplots
ggplot(data = filter(eval, metric == "TSS")) +
  geom_boxplot(aes(x = model, y = values, fill = class), alpha = 1, notch = TRUE) +
  scale_fill_brewer(palette = "Set3", name = "Model method") +
  xlab("Model algorithm") +
  ylab("True Skill Statistic (TSS)") +
  theme_bw()

## export figure
ggsave(out_performance, width = 8, height = 4, dpi = 300, bg = "transparent")

## plot histograms
#ggplot(data = eval, aes(x = values)) +
#  geom_histogram(binwidth = 0.02) +
#  geom_vline(data = mean, aes(xintercept = mean), col = 'red', size = 1, linetype = 2)+
#  xlab("Perfomance metric values") +
#  ylab("Frequency") +
#  theme_bw() +
#  theme(legend.position = "bottom", strip.text = element_text(size = 15), strip.background = element_rect(colour = "black", fill = "white")) +
#  facet_wrap(~metric)

## export figure
#ggsave(paste0("./results/figures/01_enm_performance.png"),
#       width = 8, height = 4, dpi = 300, bg = "transparent")


# 3. AREA GAIN & LOSS FIGURES  ####---------------------------------------------
## 3.1. Variation of species numbers through time periods  ####
# import enm area results
enm <- read.csv(file.path(project_repo, "results", "niche_models", "niche_models_area_results.csv")) %>%
  filter(species %in% spp_names) # select only species with migclim results

# import migclim area results
migclim <- read.csv(file.path(project_repo, "results/migclim/Migclim_area_results.csv")) %>%
  filter(species %in% spp_names)   # select only data from 2081-2100 interval

# join dataframes
data <- left_join(migclim, enm, by = c("species", "scenario", "time"), suffix = c("_migclim", "_enm")) %>%
  select(species, scenario, time, gain_loss_migclim, gain_loss_enm) %>%
  mutate(scenario = str_replace(scenario, "MPI-ESM1-2-HR_ssp", "SSP")) %>%
  rename("migclim" = "gain_loss_migclim", "enm" = "gain_loss_enm")

# prepare dataframe for ggplot by classifying each area result
# and count the number of species for each result
# ... for enm results
dt_enm <- data %>% 
  mutate(estimate = "Potential",
         class = if_else(enm >= 5, "Gain", if_else(enm <= -5, "Loss", "Stable"))) %>%
  group_by(scenario, time, class) %>%
  summarise(estimate = unique(estimate), count = n())

# ... for migclim results
dt_migclim <- data %>% 
  mutate(estimate = "Colonizable",
         class = if_else(migclim >= 5, "Gain", if_else(migclim <= -5, "Loss", "Stable"))) %>%
  group_by(scenario, time, class) %>%
  summarise(estimate = unique(estimate), count = n())

# join both dataframes into once and reorde 'class' and 'estimate'
dt <- bind_rows(dt_enm, dt_migclim) %>%
  mutate(percent = (count/160)*100,
         class = factor(class, levels = c("Gain", "Stable", "Loss")),
         estimate = factor(estimate, levels = c("Potential", "Colonizable")))

# plot the results in facet_grid
ggplot(dt) +
  geom_point(aes(x = time, y = percent, color = class), size = 2) +
  geom_line(aes(x = time, y = percent, group = class, color = class)) +
  scale_color_manual(values = c("#4E84C4", "#FFDB6D", "#D16103"), name = "Distribution area") +
  ylab("Percentage of species (%)") +
  xlab("Time period") +
  theme_bw() +
  facet_grid(estimate ~ scenario)+
  theme(axis.text.x = element_text(angle = -20, vjust = 1, hjust = 0), strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

## export figure
ggsave(out_gainloss, width = 11, height = 5, dpi = 300, bg = "transparent")


## 3.2. Differences between threatened status and endemism  ####
# import iucn status and endemism info
threat <- read.csv("./processed_data/01-1_species_list_revised.csv") %>%
  select(search.str, threat.status, domain) %>%
  setNames(c("species", "status", "domain")) %>%
  mutate(species = str_replace(species, " ", "_"),
         endemism = if_else(domain == "Mata Atlântica", "Endemic", "Not endemic"))

# check the species name that does not match with the species listed in the
# iucn dataframe (based on Flora do Brasil 2020)
diff <- setdiff(unique(data$species), unique(threat$species))

# set the synonyms between gbif and Flora do Brasil names
synonym <- data.frame(gbif = diff, flora = c("Mollinedia_schottiana", "Casearia_commersoniana"))

# replace the flora do brasil by the gbif synonyms
threat$species <- str_replace_all(threat$species, setNames(synonym$gbif, synonym$flora))

# join both dataframes into one
# NA values in threatened status was assumed as Least Concern
threat <- left_join(data, threat, by = "species") %>%
  mutate(status = replace_na(status, "LC")) %>%
  mutate(status = factor(status, levels = c("LC", "NT", "VU", "EN", "CR"))) %>%
  mutate(endemism = factor(endemism, levels = c("Not endemic", "Endemic"))) %>%
  filter(time == "2081-2100") # only 2081-2100 results

# plot boxplots
ggplot(threat) +
  geom_boxplot(aes(x = status, y = migclim, fill = status), alpha = 2, notch = FALSE) +
  scale_fill_grey(start = 1, end = 0) +
#  scale_fill_brewer(palette = "Set3", name = "Conservation status") +
  facet_wrap(~ scenario) +
  xlab("Conservation status") +
  ylab("Area gain or loss (%)") +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

## export figure
ggsave("./results/figures/02_gain_loss_conserv_status.png", width = 8, height = 3, dpi = 300, bg = "transparent")

# plot boxplots
ggplot(threat) +
  geom_boxplot(aes(x = endemism, y = migclim, fill = endemism), alpha = 2, notch = FALSE) +
  scale_fill_grey(start = 1, end = 0.5) +
  #  scale_fill_brewer(palette = "Set3", name = "Conservation status") +
  facet_wrap(~ scenario) +
  xlab("Atlantic Forest") +
  ylab("Area gain or loss (%)") +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

## export figure
ggsave("./results/figures/02_gain_loss_endemism.png", width = 8, height = 3, dpi = 300, bg = "transparent")


# 4. MAPS OF FUTURE SHIFTS ON SPECIES DISTRIBUTION  ####------------------------ 
## 4.1. Colonized/decolonized areas until 2100  ####
## create a folder for niche model results if it does not exist
ifelse(dir.exists("./results/figures/distribution_shifts"), "Results directory already exists!",
       dir.create("./results/figures/distribution_shifts", recursive = TRUE))

map(.x = final_spp,
    .f = function(x){
      
      ### 3.3.1. Import Migclim results  ####
      # import
      list <- list.files(file.path("./results/migclim", x),
                         pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep(paste0(x, "_MPI-ESM1-2-HR_ssp"), ., value = TRUE)
      
      migclim <- rast(list)
      names(migclim) <- paste0(scenarios, "_2081-2100")
      
      
      r <- migclim %>%
        as.data.frame(., xy = TRUE) %>%
        na.omit() %>%
        rename_all(~c("x", "y", "SSP126", "SSP370", "SSP585"))
      
      r_pivot <- r %>%
        pivot_longer(
          c(-x, -y),
          names_to = "scenario",
          values_to = "value"
        )
      
      
      label_colors <- c("red", "light grey", "light green", "dark green", "light grey")
      label_names <- c("decolonized", "NA", "stable", "colonized", "NA")
      
      
      ggplot() +
        geom_raster(data = r_pivot, aes(x = x, y = y, fill = as.character(value))) +
        coord_equal() +
        geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
        facet_wrap(scenario ~ .) +
        scale_fill_manual(name = "", values = label_colors, labels = label_names) +
        theme_void() +
        theme(legend.position = "bottom", strip.text = element_text(size = 15))
      
      ggsave(paste0("./results/figures/distribution_shifts/Migclim_", x, "_distribution_shifts.png"),
             width = 10, height = 4, dpi = 300, bg = "transparent")
      
      return(plot)
      
    }
)


## 4.2. Potential/colonized distribution until 2100  ####
# 
ifelse(dir.exists("./results/figures/potential_range"), "Results directory already exists!",
       dir.create("./results/figures/potential_range", recursive = TRUE))

map(.x = final_spp,
    .f = function(x){
      
      ### 3.3.1. Import Migclim results (= colonized distribution)  ####
      # import
      list <- list.files(file.path("./results/migclim", x),
                         pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      migclim <- rast(list) > 0
      names(migclim) <- apply(expand.grid(scenarios, time), 1, paste, collapse = "/")
      
      migclim <- as.list(migclim)
      names(migclim) <- apply(expand.grid(scenarios, time), 1, paste, collapse = "_")
      
      migclim <- list(SSP126 = rast(list[grep("ssp126", list)]) > 0,
                      SSP370 = rast(list[grep("ssp370", list)]) > 0,
                      SSP585 = rast(list[grep("ssp585", list)]) > 0)
      
      #values(final)[values(final) < 0] <- -1
      #values(final)[values(final) > 1] <- 2
      # %>%
      #  grep("_step_", ., value = TRUE, invert = TRUE) %>%
      #  grep(paste0(x, "_MPI-ESM1-2-HR"), ., value = TRUE)
      
      # migclim <- stack(list) > 0
      
      ### 3.3.2. Import ENM results (= potential distribution)  ####
      list <- list.files(file.path(project_repo, "results/niche_models", x), pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep("ensemble_spec_sens", ., value = TRUE) %>%
        grep("current", ., value = TRUE, invert = TRUE)
      
      enm <- rast(list)
      names(enm) <- apply(expand.grid(scenarios, time), 1, paste, collapse = "/")
      
      
      enm <- list(SSP126 = rast(list[grep("ssp126", list)]),
                  SSP370 = rast(list[grep("ssp370", list)]),
                  SSP585 = rast(list[grep("ssp585", list)]))
      
      ### 3.3.3. Difference between colonized and potential distributions  ####
      r <- map(.x = c("SSP126", "SSP370", "SSP585"),
               .f = function(x){
                 
                 r <- enm[[x]] + crop(migclim[[x]], ext(enm[[x]]))
                 
               }
      )
      
      
      #r <- enm + migclim
      r <- rast(r)
      names <- apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_")
      names <- names[order(names)]
      names(r) <- names
      
      r <- r %>%
        as.data.frame(., xy = TRUE) %>%
        na.omit() %>%
        rename_all(~c("x", "y", names))
      
      r_pivot <- r %>%
        pivot_longer(
          c(-x, -y),
          names_to = "scenario",
          values_to = "value"
        )
      
      if(length(levels(as.factor(r_pivot$value))) > 2){
        label_colors <- c("light grey", "light yellow", "light green")
        label_names <- c("NA", "potential distribution", "colonized distribution")
      }else{
        
        label_colors <- c("light grey", "light green")
        label_names <- c("NA", "colonized distribution")
      }
      ggplot() +
        geom_raster(data = r_pivot, aes(x = x, y = y, fill = as.character(value))) +
        coord_equal() +
        geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
        facet_wrap(scenario ~ .) +
        scale_fill_manual(name = "", values = label_colors, labels = label_names) +
        theme_void() +
        theme(legend.position = "bottom", strip.text = element_text(size = 15))
      
      ggsave(paste0("./results/figures/potential_range/Migclim_", x, "_potential_range.png"),
             width = 10, height = 7, dpi = 300, bg = "white")
      
      return(plot)
      
    }
)


# 5. MAP RESTORATION PRIORITY  ####---------------------------------------------
## 5.1. Map restoration and regeneration priority  ####
map(.x = scenarios,
    .f = function(x){
      
      ### 5.1.1. Map regeneration priority  ####
      ### regeneration priority maps are created based on the sum of the
      ### colonized distribution of all species (i.e., migclim results)
      # import all migclim steps between time periods (steps 120, 220, 320, 420)
      list <- list.files(file.path(in_migclim, spp_names, x), pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      # import steps as binary rasters
#      migclim <- rast(list) > 0
#      names(migclim) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
#      migclim <- as.list(migclim)
#      names(migclim) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "_")
      # import migclim steps for each time period
      migclim <- list(`2021-2040` = rast(list[grep("step_120", list)]),
                      `2041-2060` = rast(list[grep("step_220", list)]),
                      `2061-2080` = rast(list[grep("step_320", list)]),
                      `2081-2100` = rast(list[grep("step_420", list)]))
      
      migclim <- lapply(migclim, function(r){
        return(r > 0 & r < 1000)
      })
      
      sum_migclim <- lapply(migclim, function(i){
        sum_raster <- sum(i)
        return(sum_raster)
      })
      
      # export files
      ifelse(dir.exists(file.path("./results/restoration", x)), "Directory already exists!",
             dir.create(file.path("./results/restoration", x), recursive = TRUE))
      
      writeRaster(sum_migclim$`2021-2040`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2021-2040.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2041-2060`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2041-2060.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2061-2080`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2061-2080.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2081-2100`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2081-2100.tiff"), overwrite = TRUE)
      
      ### 5.1.2. Map restoration priority  ####
      ### restoration priority maps are created based on the difference between
      ### colonized distribution and potential distribution of all species
      ### (i.e., enm - migclim results)
      list <- list.files(file.path(project_repo, "results/niche_models", spp_names), pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep("ensemble_spec_sens", ., value = TRUE) %>%
        grep(x, ., value = TRUE)
      
      enm <- rast(list)
      names(enm) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
      list <- list.files(file.path(in_migclim, spp_names, x),
                         pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      migclim <- rast(list) > 0 & rast(list) < 1000
      names(migclim) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
      #potential <- mask(enm, crop(migclim, extent(enm)), maskvalues = 0, updatevalue = 0, inverse = FALSE)
      potential <- (enm + crop(migclim, ext(enm))) == 1
      names(potential) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
      
      r_list <- list(`2021-2040` = sum(potential[[grep("2021.2040", names(potential))]]),
                     `2041-2060` = sum(potential[[grep("2041.2060", names(potential))]]),
                     `2061-2080` = sum(potential[[grep("2061.2080", names(potential))]]),
                     `2081-2100` = sum(potential[[grep("2081.2100", names(potential))]]))
      
      ### export files
      ifelse(dir.exists(file.path("./results/restoration", x)), "Directory already exists!",
             dir.create(file.path("./results/restoration", x), recursive = TRUE))
      
      writeRaster(r_list$`2021-2040`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2021-2040.tiff"), overwrite = TRUE)
      writeRaster(r_list$`2041-2060`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2041-2060.tiff"), overwrite = TRUE)
      writeRaster(r_list$`2061-2080`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2061-2080.tiff"), overwrite = TRUE)
      writeRaster(r_list$`2081-2100`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2081-2100.tiff"), overwrite = TRUE)
      
    }
)


## 5.2. Mapping continuous indices within BHRD  ####
### 5.2.1. Figure: regeneration priority within BHRD  ####
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("regeneration_priority", ., value = TRUE)

r <- rast(list)
names <- apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_")
names <- names[order(names)]
names(r) <- names

r <- r %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit() %>%
  rename_all(~c("x", "y", names))

r_pivot <- r %>%
  pivot_longer(
    c(-x, -y),
    names_to = "scenario",
    values_to = "value"
  ) %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"))

if(length(levels(as.factor(r_pivot$value))) > 2){
  label_colors <- c("light grey", "light yellow", "light green")
  label_names <- c("NA", "potential distribution", "colonized distribution")
}else{
  
  label_colors <- c("light grey", "light green")
  label_names <- c("NA", "colonized distribution")
}

ggplot() +
  geom_raster(data = r_pivot, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  geom_sf(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_gradientn(colours = c("white", "lightgreen", "darkgreen"), name = "Regeneration priority") +
  #  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

ggsave("./results/figures/03_regeneration_priority.png",
       width = 10, height = 7, dpi = 300, bg = "white")


## 5.2.2. Figure: restoration priority within BHRD  ####
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("restoration_priority", ., value = TRUE)


r <- rast(list)
names <- apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_")
names <- names[order(names)]
names(r) <- names

r <- r %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit() %>%
  rename_all(~c("x", "y", names))

r_pivot <- r %>%
  pivot_longer(
    c(-x, -y),
    names_to = "scenario",
    values_to = "value"
  ) %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"))

if(length(levels(as.factor(r_pivot$value))) > 2){
  label_colors <- c("light grey", "light yellow", "light green")
  label_names <- c("NA", "potential distribution", "colonized distribution")
}else{
  
  label_colors <- c("light grey", "light green")
  label_names <- c("NA", "colonized distribution")
}

ggplot() +
  geom_raster(data = r_pivot, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  geom_sf(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), name = "Restoration priority") +
#  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

ggsave(out_restoration,
       width = 10, height = 7, dpi = 300, bg = "white")


## 5.3. Restoration priority per municipality  ####
### 5.3.1. Priority indices per municipality  ####
# import restoration priority map
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("restoration_priority", ., value = TRUE)
#  grep("2081-2100", ., value = TRUE)

# import as rasters
r_list <- as.list(rast(list))

# assign scenario name to each raster
names(r_list) <- sort(apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_"))

# import shapefile of municipalities within the river basin
munic <- vect(in_munic) %>%
  crop(., ext(r_list[[1]]))

# calculate the priority indices for each municipality
# (1) 'priority_mean' index is calculated by the mean value of all cells within
#    the municipality limits
# (2) 'priority_sum' index is calculated by the sum of all cells within the
#    municipality limits
values <- map(.x = names(r_list),
              .f = function(x){
                
                v <- terra::extract(r_list[[x]], munic) %>%
                  rename(., ID = "ID", value = names(r_list[[x]])) %>%
                  group_by(ID) %>%
                  summarise(priority_mean = mean(value, na.rm = TRUE),
                            priority_sum = sum(value, na.rm = TRUE))
                
                munic$index_sum <- v$priority_sum
                names(munic) <- c(head(names(munic), -1), paste0("Index_", x))
                
                results <- data.frame(scenario = x, municipality = munic$NOMEMUNICP, uf = munic$NOMEUF, v) %>%
                  drop_na() %>%
                  group_by(municipality) %>%
                  summarize(uf = unique(uf),
                            scenario = scenario,
                            priority_mean = mean(priority_mean, na.rm = TRUE),
                            priority_sum = sum(priority_sum, na.rm = TRUE)) %>%
                  arrange(desc(priority_mean))
                
                
                return(results)
                
              })

# join all dataframes and split scenarios and time period
values <- bind_rows(values) %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"))

# export the results
write_csv(values, "./results/restoration/restoration_priority_by_munic.csv")


## 5.3.2. Map restoration priority per municipality  ####
#values <- values %>%
#  pivot_wider(names_from = "scenario", values_from = c("priority_mean", "priority_sum"))

### 5.3.1. Horizontal graphic for municipilaties rank  ####
# import municipalities features from shapefile
#dt <- as.data.frame(munic) %>%
 # rename("municipality" = "NOMEMUNICP") %>%
#  select(municipality)

# join municipalities features and priority indices
#dt <- left_join(values, dt, by = "municipality")

# 
#dt <- dt[which(dt$'priority_mean_SSP585_2081-2100' != 0),] %>%
#  rename(Index = 'priority_mean_SSP585_2081-2100') %>%
#  arrange(desc(Index)) %>%
#  head(., 10)

revals <- values %>%
  arrange(scenario, time, desc(priority_mean)) %>%
  group_by(scenario, time) %>%
  slice(1:10)


#
ggplot(revals, aes(x = priority_mean, y = reorder(municipality, priority_mean))) +
  geom_col() +
  xlab("Restoration priority index") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(size = 15)) +
  facet_grid(scenario ~ time, scales = "free")

ggsave("./results/restoration/Restoration_priority_per_municipality_bars.png",
       width = 10, height = 5, dpi = 300, bg = "transparent")

### 5.3.4. Maps with priority index per municipality  ####
# shapefile with municipilaties within the river basin
munic_bhrd <- crop(munic, doce)
crs(munic_bhrd) <- "+init=EPSG:4326"

# 
vals <- values %>%
  rename("NOMEMUNICP" = "municipality")

data <- list()

for(s in c("SSP126", "SSP370", "SSP585")){
  for(t in time){
    
    v <- filter(vals, scenario == s & time == t)
    
    raster <- merge(munic_bhrd, v, by.x = "NOMEMUNICP") %>%
      rasterize(., rast(., res = 0.005), field = "priority_mean")
    
    data[[paste(s, t, sep = "_")]] <- raster
    
  }
}

data <- lapply(data, function(x){
  df <- as.data.frame(x) %>%
    bind_cols(., crds(x))
  return(df)
  }
  )

data <- data %>%
  bind_rows(., .id = "scenario") %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"), )


ggplot() +
  geom_raster(data = data, aes(x = x, y = y, fill = priority_mean)) +
  coord_equal() +
  geom_spatvector(data = munic_bhrd, color = "white", linewidth = 0.1, fill = NA) +
  geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_gradientn(colours = c("darkblue", "blue", "purple", "yellow", "orange", "red", "darkred"), name = "Regeneration priority mean") +
#  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), name = "Restoration priority") +
#  theme_void() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

ggsave("./results/figures/04_restoration_priority_munic.png",
       width = 10, height = 7, dpi = 300, bg = "white")



