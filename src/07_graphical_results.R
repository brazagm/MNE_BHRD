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
out_gainloss_enm <- "./results/figures/02-1_niche_models_area_results.png"
out_gainloss_migclim <- "./results/figures/02-2_migclim_area_results.png"
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
  mutate(class = if_else(model %in% c("bioclim", "domain", "mahal"), "presence-only", if_else(model %in% c("GLM", "GBM", "RF"), "presence-absence", if_else(model == "maxent", "presence-background", "ensemble"))))

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
eval$model <- factor(eval$model , levels = c("bioclim", "domain", "mahal", "GLM", "GBM", "RF", "maxent", "ensemble"))
eval$class <- factor(eval$class , levels = c("presence-only", "presence-absence", "presence-background", "ensemble"))

## plot boxplots
ggplot(data = filter(eval, metric == "TSS")) +
  geom_boxplot(aes(x = model, y = values, fill = class), alpha = 1, notch = TRUE) +
  scale_fill_brewer(palette = "Set3", name = "Model class") +
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


# 3. AREA GAIN & LOSS BAR PLOTS  ####-------------------------------------------
## 3.1. Gain-loss area graph for niche models (climate-only)  ####
## import area gain and losses results
enm <- read.csv(file.path(project_repo, "results", "niche_models", "niche_models_area_results.csv")) %>%
  filter(time == "2081-2100") %>%  # select only data from 2081-2100 interval
  filter(species %in% spp_names) # select only species with migclim results
#  mutate(species = str_replace(species, "_", " "))

## prepare dataframe for bar plot
## reorder the entries by an unique ID (species name + scenario)
## and classify each entry as area gain or loss
data <- enm %>% 
  arrange(scenario, gain_loss) %>%
  mutate(scenario = str_replace(scenario, "MPI-ESM1-2-HR_ssp", "SSP")) %>%
  mutate(id = paste(species, scenario, sep = "_")) %>%
  mutate(id = factor(id, levels = id),
         class = if_else(gain_loss >= 0, "gain", "loss"))

## plot bars (vertical option)
ggplot(data, aes(x = gain_loss, y = id, fill = class)) +
  geom_col(show.legend = FALSE, width = 1) +
  geom_vline(aes(xintercept = 0), col = 'black', size = 0.2) +
  scale_x_continuous(breaks = c(-100, 0, 100, 250, 500, 1000, 1400)) +
  xlab("Area gain and loss until 2100 (%)") +
  ylab("Species rank") + # do not show Label of Y
  scale_fill_manual(breaks = c("gain", "loss"), values = c("#00798c", "#d1495b")) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 15),  strip.background = element_rect(colour = "black", fill = "white")) +
  facet_wrap(~scenario, scales = "free_y")

## export image
ggsave(out_gainloss_enm, width = 10, height = 4, dpi = 300, bg = "transparent")

## horizontal option
#ggplot(data, aes(x = gain_loss, y = species, fill = class)) +
#  geom_col(show.legend = FALSE, width = 1) +
#  geom_vline(aes(xintercept = 0), col = 'black', size = 0.2) +
#  coord_flip() +
#  scale_x_continuous(breaks = c(-100, 0, 100, 250, 500, 1000, 1400)) +
#  xlab("Area gain and loss until 2100 (%)") +
#  ylab("Species rank") + # do not show Label of Y
#  scale_fill_manual(breaks = c("gain", "loss"), values = c("#00798c", "#d1495b")) +
#  theme_bw() +
#  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
#        panel.grid.major.x = element_blank(), strip.text = element_text(size = 15),  strip.background = element_rect(colour = "black", fill = "white")) +
#  facet_wrap(~scenario, scales = "free_x")

## export image
#ggsave(file.path("results", "figures", paste0("niche_models_area_results.png")), #width = 12, height = 7, dpi = 300, bg = "transparent")


## 3.2. Gain-loss area graph for distribution (with migclim)  ####
## import area gain and losses results from niche models
#enm <- read.csv(file.path(project_repo, "results", "niche_models", "niche_models_area_results.csv")) %>%
#  filter(time == "2081-2100") %>%
#  filter(species %in% spp_names) %>%  # select only data from 2081-2100 interval
#  mutate(species = str_replace(species, "_", " "))

## import area gain and losses results from migclim
migclim <- read.csv(file.path(project_repo, "results/migclim/Migclim_area_results.csv")) %>%
  filter(species %in% spp_names)   # select only data from 2081-2100 interval
#  mutate(species = str_replace(species, "_", " "))

## join dataframes
data <- left_join(migclim, enm, by = c("species", "scenario"), suffix = c("_migclim", "_enm")) %>%
  select(species, scenario, gain_loss_migclim, gain_loss_enm) %>%
#  pivot_longer(c(gain_loss_migclim, gain_loss_enm), names_to = "model", values_to = "values") %>%
  mutate(scenario = str_replace(scenario, "MPI-ESM1-2-HR_ssp", "SSP")) %>%
  rename("migclim" = "gain_loss_migclim", "enm" = "gain_loss_enm")

## prepare dataframe for bar plot
## reorder the entries by an unique ID (species name + scenario)
## and classify each entry as area gain or loss
data <- data %>% 
  arrange(scenario, enm) %>%
  mutate(id = paste(species, scenario)) %>%
  mutate(id = factor(id, levels = id),
         class_migclim = if_else(migclim >= 0, "gain", "loss"),
         class_enm = if_else(enm >= 0, "gain", "loss"))

## vertical option
ggplot(data) +
  geom_col(aes(x = enm, y = id), alpha = 0.7, show.legend = FALSE, width = 1) +
  geom_col(aes(x = migclim, y = id, fill = class_migclim), alpha = 0.5, show.legend = FALSE, width = 1) +
  geom_vline(aes(xintercept = 0), col = 'black', size = 0.2) +
  scale_x_continuous(breaks = c(-100, 0, 100, 250, 500, 1000, 1400)) +
  xlab("Area gain and loss until 2100 (%)") +
  ylab("Species rank") + # do not show Label of Y
  scale_fill_manual(breaks = c("gain", "loss"), values = c("#00798c", "#d1495b")) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        panel.grid.major.y = element_blank(), strip.text = element_text(size = 15),  strip.background = element_rect(colour = "black", fill = "white")) +
  facet_wrap(~scenario, scales = "free_y")

## export figure
ggsave(out_gainloss_migclim, width = 12, height = 5, dpi = 300, bg = "transparent")


# 3. AREA GAIN & LOSS BAR PLOTS  ####-------------------------------------------
migclim



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
      
      migclim <- list(`2021-2040` = rast(list[grep("step_120", list)]) > 0,
                      `2041-2060` = rast(list[grep("step_220", list)]) > 0,
                      `2061-2080` = rast(list[grep("step_320", list)]) > 0,
                      `2081-2100` = rast(list[grep("step_420", list)]) > 0)
      
      sum_migclim <- lapply(migclim, function(i){
        sum_raster <- sum(i)
        return(sum_raster)
      })
      
      ### export files
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
      
      migclim <- rast(list) > 0
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
      
      writeRaster(sum_migclim$`2021-2040`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2021-2040.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2041-2060`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2041-2060.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2061-2080`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2061-2080.tiff"), overwrite = TRUE)
      writeRaster(sum_migclim$`2081-2100`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2081-2100.tiff"), overwrite = TRUE)
      
    }
)


## 5.2.  Figure: restoration priority within BHRD  ####
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
  )

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
  facet_wrap(scenario ~ .) +
  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), name = "Restoration priority") +
  theme_void() +
  theme(legend.position = "bottom", strip.text = element_text(size = 15))

ggsave(out_restoration,
       width = 10, height = 7, dpi = 300, bg = "transparent")


## 5.3. Restoration priority per municipality  ####
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("restoration_priority", ., value = TRUE) %>%
  grep("2081-2100", ., value = TRUE)

r_list <- as.list(rast(list))
names(r_list) <- paste0(c("SSP126", "SSP370", "SSP585"), "_2081-2100")

munic <- vect(in_munic) %>%
  crop(., ext(r_list[[1]]))

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
                  arrange(desc(priority_sum))
                
                
                return(results)
                
              })




values <- bind_rows(values)

write_csv(values, "./results/restoration/restoration_priority_by_munic.csv")


## 5.4. Map restoration priority per municipality  ####
values <- values %>%
  pivot_wider(names_from = "scenario", values_from = c("priority_mean", "priority_sum"))

### 5.4.1. Horizontal graphic for municipilaties rank  ####
dt <- as.data.frame(munic) %>%
  rename("municipality" = "NOMEMUNICP")

dt <- left_join(values, dt, by = "municipality")

dt <- dt[which(dt$'priority_mean_SSP585_2081-2100' != 0),] %>%
  rename(Index = 'priority_mean_SSP585_2081-2100') %>%
  arrange(desc(Index)) %>%
  head(., 10)

ggplot(dt, aes(x = Index, y = reorder(municipality, Index))) +
  geom_col() +
  xlab("Restoration priority index") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(size = 15))

ggsave("./results/restoration/Restoration_priority_per_municipality_bars.png",
       width = 10, height = 5, dpi = 300, bg = "transparent")

### 5.4.2. 
values <- values %>%
  rename("NOMEMUNICP" = "municipality")

merge(munic, values, by.x = "NOMEMUNICP")

dt <- left_join(munic, values, by = "NOMEMUNICP")


doce <- as(doce, "Spatial")
doce <- vect(doce)

dt <- crop(dt, doce)
crs(dt) <- crs(doce)

dt2 <- map(.x = c("priority_mean_SSP126_2081-2100", "priority_mean_SSP370_2081-2100", "priority_mean_SSP585_2081-2100"),
           .f = function(x){
             rasterize(dt, rast(dt, res = 0.005), field = x)
           })

#names(dt2) <- c("SSP126", "SSP370", "SSP585")

xy <- crds(dt2[[1]])

dt2 <- as.data.frame(dt2) %>%
  cbind(xy, .) %>%
  pivot_longer(cols = c("priority_mean_SSP126_2081.2100", "priority_mean_SSP370_2081.2100", "priority_mean_SSP585_2081.2100"), names_to = "variable", values_to = "values")

#munic <- st_set_crs(munic, st_crs(doce))
#m <- sf::st_intersection(munic, doce)

#m <- vect(intersect(munic, doce))

#munic[which(munic$Index_SSP126_2081-2100 > 0),]


ggplot() +
  geom_spatvector(data = dt, aes(fill = `priority_mean_SSP126_2081-2100`), color = "black", linewidth = 0.5) +
  geom_sf(aes(fill = Index_SSP585_2081.2100)) +
  #  coord_equal() +
  #  geom_sf(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  #  facet_wrap(scenario ~ .) +
  scale_fill_gradientn(colours = c("white", "yellow", "red"), name = "Restoration priority") +
  ggtitle("Scenario", subtitle = "SSP585 2081-2100") +
  theme_void() +
  theme(legend.position = "bottom", strip.text = element_text(size = 15))

ggplot() +
  geom_raster(data = dt2, aes(x = x, y = y, fill = values)) +
  coord_equal() +
  geom_spatvector(data = dt, color = "black", linewidth = 0.2, fill = NA) +
  geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_wrap(variable ~ .) +
#  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), name = "Restoration priority") +
  theme_void() +
  theme(legend.position = "bottom", strip.text = element_text(size = 15))

ggsave("./results/restoration/Restoration_priority_munic.png",
       width = 10, height = 7, dpi = 300, bg = "transparent")

