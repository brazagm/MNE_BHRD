###
###  GRAPHICAL RESULTS
###
###  Description: Produce figures, graphics and maps based on the outputs of 
###     other project's scripts.
###
###  Created & Edited by: Alan Braz (@brazagm)
###
###  Observations:
###
###  Next tasks:
###     (1) Revise section 4.
###

# load packages
library(ggpattern)
library(terra)
library(tidyterra)
library(tidyverse)


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories for model results
repo <- "/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD"

## inputs
in_models <- file.path(repo, "results", "niche_models")
in_migclim <- file.path(repo, "results", "migclim")
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
## species to omit because they do not occur within BHRD
omit_spp <- c("Enterolobium_timbouva", "Piptocarpha_angustifolia")

## get species names with migclim results
spp_names <- read.csv(in_migclim_area) %>%
  filter(!species %in% omit_spp) %>%
  pull(species) %>%
  unique()


## 1.3. Future projection names  ####
time <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
scenarios <- c(paste0("MPI-ESM1-2-HR_", c("ssp126", "ssp370", "ssp585")))
futures <- apply(expand.grid(time, scenarios), 1, paste, collapse = "/")


# 2. NICHE MODEL PERFORMANCES  ####---------------------------------------------
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


# 3. AREA GAIN & LOSS RESULTS  ####---------------------------------------------
## 3.1. Organize data  ####
# import enm area results
enm <- read.csv(file.path(repo, "results", "niche_models", "niche_models_area_results.csv")) %>%
  filter(species %in% spp_names)  %>%
  unique() # select only species with migclim results

# import migclim area results
migclim <- read.csv(file.path(repo, "results/migclim/Migclim_area_results.csv")) %>%
  filter(species %in% spp_names) %>% # select only data from 2081-2100 interval
  unique()

# join dataframes
data <- left_join(migclim, enm, by = c("species", "scenario", "time"), suffix = c("_migclim", "_enm")) %>%
  select(species, scenario, time, gain_loss_migclim, gain_loss_enm) %>%
#  mutate(scenario = str_replace(scenario, "MPI-ESM1-2-HR_ssp", "SSP")) %>%
  mutate(scenario = str_replace_all(scenario, pattern = c("MPI-ESM1-2-HR_ssp126" = "Optimistic", "MPI-ESM1-2-HR_ssp370" = "Intermediate", "MPI-ESM1-2-HR_ssp585" = "Pessimistic"))) %>%
  mutate(scenario = factor(scenario, levels = c("Optimistic", "Intermediate", "Pessimistic"))) %>%
  rename("migclim" = "gain_loss_migclim", "enm" = "gain_loss_enm")


## 3.1. Variation of species numbers through time periods  ####
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

# test for differences between groups using Kruskal-Wallis non-parametric test
kruskal_threat <- map(.x = c(levels(threat$scenario)),
                      .f = function(x){
                        sub <- filter(threat, scenario == x)
                        return(kruskal.test(migclim ~ status, data = sub))
                     })

wilcox_threat <- map(.x = c(levels(threat$scenario)),
    .f = function(x){
      
      sub <- filter(threat, scenario == x)
      kruskal.test(migclim ~ status, data = sub)
      
      return(pairwise.wilcox.test(sub$migclim, sub$status, p.adjust.method = "BH"))
      
    })

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

# test for differences between groups using Kruskal-Wallis non-parametric test
kruskal_endemism <- map(.x = c(levels(threat$scenario)),
                      .f = function(x){
                        sub <- filter(threat, scenario == x)
                        return(kruskal.test(migclim ~ endemism, data = sub))
                      })

wilcox_endemism <- map(.x = c(levels(threat$scenario)),
                     .f = function(x){
                       
                       sub <- filter(threat, scenario == x)
                       return(pairwise.wilcox.test(sub$migclim, sub$endemism, p.adjust.method = "BH"))
                       
                     })

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


# 4. SHIFTS ON SPECIES DISTRIBUTION  ####--------------------------------------- 
## 4.1. Colonized/decolonized areas until 2100  ####
## These maps show colonized and decolonized regions from current distribution
## to future distribution in 2100 based on the migclim results for each species.
## Results are shown for each socioeconomic pathway. Note that colonized and
## decolonized areas are defined in differences between current and distribution
## in 2100.
##
## create a folder for niche model results if it does not exist
ifelse(dir.exists("./results/figures/distribution_shifts"), "Results directory already exists!",
       dir.create("./results/figures/distribution_shifts", recursive = TRUE))

## create maps for colonized/decolonized areas for each species
map(.x = spp_names,
    .f = function(x){
      
      ### 4.1.1. Import Migclim results  ####
      # import file names
      list <- list.files(file.path(in_migclim, x),
                         pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep(paste0(x, "_MPI-ESM1-2-HR_ssp"), ., value = TRUE)
      
      # import rasters and give time period name
      migclim <- rast(list)
      names(migclim) <- paste0(scenarios, "_2081-2100")
      
      # raster as dataframe
      r <- migclim %>%
        as.data.frame(., xy = TRUE) %>%
        na.omit() %>%
        rename_all(~c("x", "y", "SSP126", "SSP370", "SSP585"))
      
      # create future distribution for each scenario, i.e., "-1", "0" and
      # "30'000" values are converted in "absences" and "1' and '2' values are
      # converted in "presences"
      future_distr <- r %>%
        pivot_longer(
          c(-x, -y),
          names_to = "scenario",
          values_to = "value"
        ) %>%
        mutate(value = str_replace_all(value, pattern = c("-1" = "0", "30000" = "0", "2" = "1")))
      
      # create a single current distribution shapefile, , i.e., "0" and
      # "30'000" values are converted in "absences" and "1' and '2' values are
      # converted in "presences"
      current_distr <- list.files(file.path(in_models, x, "current"),
                                  pattern = "ensemble_spec_sens.tif", full.names = TRUE) %>%
        rast() %>%
        terra::crop(., ext(doce)) %>%
        terra::mask(., doce) %>%
        subst(., 0, NA) %>%
        terra::as.polygons(aggregate = TRUE, values = TRUE)
      
      # give categorical value for variable 'value' in the polygon vector
      values(current_distr) <- "Presence"

      # set legend colors and names
      label_title <- str_replace(x, "_", " ")
      label_colors <- c("light grey", "#315e26")
      label_names <- c("Absence", "Presence")
      label_pattern_names <- c("Species'
distribution")
      
      # plot distribution shifts for each scenario
      ggplot() +
        geom_raster(data = future_distr, aes(x = x, y = y, fill = as.character(value))) +
        coord_equal() +
        geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
        geom_sf_pattern(data = current_distr, aes(pattern = value), show.legend = TRUE, color = "black", linewidth = 0.2, fill = NA, pattern_size = .05, pattern_density = .12, pattern_spacing = .01, pattern_fill = "black") +
        facet_wrap(scenario ~ .) +
        scale_fill_manual(name = "2081-2100 scenario", values = label_colors, labels = label_names) +
        scale_pattern_manual(name = "Current scenario", values = "stripe", labels = label_pattern_names) +
        ggtitle(label_title) +
        theme_void() +
        theme(plot.title = element_text(face = "bold.italic", size = 20, hjust = 0.5, vjust = 3.5), legend.title = element_text(face = "bold"), legend.position = "bottom", strip.text = element_text(size = 15), legend.spacing.x = unit(2, "cm"), legend.spacing.y = unit(1, "cm")) +
        guides(fill = guide_legend(order = 2, override.aes = list(pattern = "none"), title.position = "top", direction = "vertical", byrow = TRUE), 
               pattern = guide_legend(order = 1, title.position = "top", direction = "vertical"))

      # export png image
      ggsave(paste0("./results/figures/distribution_shifts/Migclim_", x, "_distribution_shifts.png"), width = 11, height = 5, dpi = 300, bg = "white")
      
      # map() return
      return(plot)
      
    }
)


## 4.2. Potential/colonized distribution until 2100  ####
## These maps show potential climatic distribution and colonized distribution
## in 2100 based on, respectively, the enm and migclim results for each species.
## Results are shown for each socioeconomic pathway. Note that potentical and
## colonized areas are defined in differences between enm and migclim results
## in 2081-2100 scenario.
##
## create a folder for niche model results if it does not exist
ifelse(dir.exists("./results/figures/potential_range"), "Results directory already exists!",
       dir.create("./results/figures/potential_range", recursive = TRUE))

## create maps for colonized/decolonized areas for each species
map(.x = spp_names,
    .f = function(x){
      
      ### 4.2.1. Import Migclim results (= colonized distribution)  ####
      # import file paths
      list <- list.files(file.path(in_migclim, x),
                         pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      # import migclim results as raster
      migclim <- list(SSP126 = rast(list[grep("ssp126", list)]) > 0,
                      SSP370 = rast(list[grep("ssp370", list)]) > 0,
                      SSP585 = rast(list[grep("ssp585", list)]) > 0)
      
      
      ### 4.2.2. Import ENM results (= potential distribution)  ####
      # import file paths
      list <- list.files(file.path(repo, "results/niche_models", x), pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep("ensemble_spec_sens", ., value = TRUE) %>%
        grep("current", ., value = TRUE, invert = TRUE)
      
      # import as binary rasters
      enm <- list(SSP126 = rast(list[grep("ssp126", list)]),
                  SSP370 = rast(list[grep("ssp370", list)]),
                  SSP585 = rast(list[grep("ssp585", list)]))
      
      
      ### 4.4.3. Difference between colonized and potential distributions  ####
      # calculate differences between potential and colonized distributions
      # '0' = absence, '1' = potential, '2' = colonized
      r <- map(.x = c("SSP126", "SSP370", "SSP585"),
               .f = function(x){
                 
                 r <- enm[[x]] + crop(migclim[[x]], ext(enm[[x]]))
                 names(r) <- paste(x, c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), sep = "_")
                 return(r)
                 
               }
      ) %>% rast()
      
      # raster as longer dataframe
      r_pivot <- r %>%
        as.data.frame(., xy = TRUE) %>%
        na.omit() %>%
        pivot_longer(c(-x, -y), names_to = "scenario", values_to = "value") %>%
        separate("scenario", into = c("scenario", "time"), sep = "_")
      
      # set legend colors and names
      label_title <- str_replace(x, "_", " ")
      if(length(levels(as.factor(r_pivot$value))) > 2){
        label_colors <- c("light grey", "#FCCF55", "#315e26")
        label_names <- c("", "Potential distribution", "Colonized distribution")
      }else{
        
        label_colors <- c("light grey", "#315e26")
        label_names <- c("", "Colonized distribution")
      }
      
      # plot potential/colonized distribution for each scenario
      ggplot() +
        geom_raster(data = r_pivot, aes(x = x, y = y, fill = as.character(value))) +
        coord_equal() +
        geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
        facet_grid(scenario ~ time, switch = "y") +
        scale_fill_manual(name = "", values = label_colors, labels = label_names) +
        ggtitle(label_title) +
        guides(fill = guide_legend(override.aes = list(fill = c("white", label_colors[-1])))) +
      theme(plot.title = element_text(face = "bold.italic", size = 20, hjust = 0.5, vjust = 3.5), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 15),  strip.background = element_rect(colour = "white", fill = "white"), panel.background = element_blank())
      
      # export png image
      ggsave(paste0("./results/figures/potential_range/Migclim_", x, "_potential_range.png"),
             width = 11, height = 9, dpi = 300, bg = "white")
      
      # map() return
      return(plot)
      
    }
)


# 5. RESTORATION/REGENERATION PRIORITY MAPS  ####---------------------------------------------
## 5.1. Create restoration and regeneration priority rasters  ####
map(.x = scenarios,
    .f = function(x){
      
      ### 5.1.1. Regeneration priority rasters ####
      ### regeneration priority maps are created based on the sum of the
      ### colonized distribution of all species (i.e., migclim results)
      # import all migclim steps between time periods (steps 120, 220, 320, 420)
      list <- list.files(file.path(in_migclim, spp_names, x), pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      # import migclim steps for each time period (note that 'step' names were
      # listed into groups of 'scenario' names)
      migclim <- list(`2021-2040` = rast(list[grep("step_120", list)]),
                      `2041-2060` = rast(list[grep("step_220", list)]),
                      `2061-2080` = rast(list[grep("step_320", list)]),
                      `2081-2100` = rast(list[grep("step_420", list)]))
      
      # binarize rasters (presence/absence in the time scenario)
      # if value is between [0, 30.000] = 1; else value = 0
      # IMPORT: I don't know why max values in raster files reaches only 29'920...
      # so I rounded for 25'000 because of that
      migclim <- lapply(migclim, function(r){
        return(r > 0 & r < 25000)
      })
      
      # sum all present species by cell
      sum_migclim <- lapply(migclim, function(i){
        sum_raster <- sum(i)
        return(sum_raster)
      })
      
      # rescale values (160 is the maximum number of species)
      resc_migclim <- lapply(sum_migclim, function(i){
        values <- scales::rescale(values(i), to = c(0, 1), from = c(0, 160))
        values(i) <- values 
        return(i)
      })
      
      # check if output directory already exists or create it
      ifelse(dir.exists(file.path("./results/restoration", x)), "Directory already exists!",
             dir.create(file.path("./results/restoration", x), recursive = TRUE))
      
      # export files
      writeRaster(resc_migclim$`2021-2040`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2021-2040.tiff"), overwrite = TRUE)
      writeRaster(resc_migclim$`2041-2060`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2041-2060.tiff"), overwrite = TRUE)
      writeRaster(resc_migclim$`2061-2080`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2061-2080.tiff"), overwrite = TRUE)
      writeRaster(resc_migclim$`2081-2100`, paste0("./results/restoration/", x, "/regeneration_priority_", x, "_2081-2100.tiff"), overwrite = TRUE)
      
      
      ### 5.1.2. Restoration priority rasters  ####
      ### restoration priority maps are created based on the difference between
      ### colonized distribution and potential distribution of all species
      ### (i.e., enm - migclim results)
      list <- list.files(file.path(repo, "results/niche_models", spp_names), pattern = ".tif", recursive = TRUE, full.names = TRUE) %>%
        grep("ensemble_spec_sens", ., value = TRUE) %>%
        grep(x, ., value = TRUE)
      
      # import enms and rename them
      enm <- rast(list)
      names(enm) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
      # list migclim results
      list <- list.files(file.path(in_migclim, spp_names, x),
                         pattern = ".asc", recursive = TRUE, full.names = TRUE) %>%
        grep("step_120|step_220|step_320|step_420", ., value = TRUE)
      
      # binarize rasters (presence/absence in the time scenario)
      # if value is between [0, 30.000] = 1; else value = 0
      # IMPORT: I don't know why max values in raster files reaches only 29'920...
      # so I rounded for 25'000 because of that
      migclim <- rast(list) > 0 & rast(list) < 25000
      names(migclim) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/") # rename
      
      # calculate unoccupied potential distribution for each species in each scenario
      potential <- (enm + crop(migclim, ext(enm))) == 1
      names(potential) <- apply(expand.grid(time, spp_names), 1, paste, collapse = "/")
      
      # sum number of species with unoccupied potential distribution for each cell
      r_list <- list(`2021-2040` = sum(potential[[grep("2021.2040", names(potential))]]),
                     `2041-2060` = sum(potential[[grep("2041.2060", names(potential))]]),
                     `2061-2080` = sum(potential[[grep("2061.2080", names(potential))]]),
                     `2081-2100` = sum(potential[[grep("2081.2100", names(potential))]]))
      
      # rescale values (160 is the maximum number of species)
      resc_r <- lapply(r_list, function(i){
        values <- scales::rescale(values(i), to = c(0, 1), from = c(0, 160))
        values(i) <- values 
        return(i)
      })
      
      # check if output directory already exists or create it
      ifelse(dir.exists(file.path("./results/restoration", x)), "Directory already exists!",
             dir.create(file.path("./results/restoration", x), recursive = TRUE))
      
      # export files
      writeRaster(resc_r$`2021-2040`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2021-2040.tiff"), overwrite = TRUE)
      writeRaster(resc_r$`2041-2060`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2041-2060.tiff"), overwrite = TRUE)
      writeRaster(resc_r$`2061-2080`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2061-2080.tiff"), overwrite = TRUE)
      writeRaster(resc_r$`2081-2100`, paste0("./results/restoration/", x, "/restoration_priority_", x, "_2081-2100.tiff"), overwrite = TRUE)
      
    }
)


## 5.2. Continuous restoration/regeneration indices  ####
### 5.2.1. Figure: continuous regeneration priority  ####
# list all regeneration index rasters
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("regeneration_priority", ., value = TRUE)

# import rasters
r <- rast(list)

# rename them
names <- apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_")
names <- names[order(names)]
names(r) <- names

# raster as dataframe
r <- r %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit() %>%
  rename_all(~c("x", "y", names))

# enlongate dataframe
r_pivot <- r %>%
  pivot_longer(
    c(-x, -y),
    names_to = "scenario",
    values_to = "value"
  ) %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"))

# plot
ggplot() +
  geom_raster(data = r_pivot, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  geom_sf(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_stepsn(colours = c("white", "lightgreen", "darkgreen", "#081108"), breaks = c(0, .2, .4, .6, .8, 1), name = "Regeneration priority", limits = c(0,1), labels = c("0", ".20", ".40", ".60", ".80", "1")) +
#  scale_fill_gradientn(colours = c("white", "lightgreen", "darkgreen"), name = "Regeneration priority") +
  #  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

# export results
ggsave("./results/figures/03_regeneration_priority.png",
       width = 10, height = 7, dpi = 300, bg = "white")


## 5.2.2. Figure: continuous restoration priority  ####
# list all restoration index rasters
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("restoration_priority", ., value = TRUE)

# import rasters
r <- rast(list)

# rename them
names <- apply(expand.grid(c("SSP126", "SSP370", "SSP585"), time), 1, paste, collapse = "_")
names <- names[order(names)]
names(r) <- names

# rasters as dataframe
r <- r %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit() %>%
  rename_all(~c("x", "y", names))

# enlongate dataframe
r_pivot <- r %>%
  pivot_longer(
    c(-x, -y),
    names_to = "scenario",
    values_to = "value"
  ) %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"))

# plot
ggplot() +
  geom_raster(data = r_pivot, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  geom_sf(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_stepsn(colours = c("white", "yellow", "orange", "red", "darkred"), breaks = c(0, .05, .1, .15, .2), name = "Restoration priority", limits = c(0,0.2), labels = c("0", ".05", ".10", ".15", ".20")) +
#  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), name = "Restoration priority") +
#  theme_bw() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

# export results
ggsave(out_restoration,
       width = 10, height = 7, dpi = 300, bg = "white")


## 5.3. Restoration priority per municipality  ####
### 5.3.1. Calculate restoration priority indices per municipality  ####
# import restoration priority map
list <- list.files("./results/restoration", pattern = ".tiff", recursive = TRUE, full.names = TRUE) %>%
  grep("restoration_priority", ., value = TRUE)

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


## 5.3.2. Horizontal graphic for municipalities rank   ####
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


### 5.3.3. Figure: maps with priority index per municipality  ####
# shapefile with municipilaties within the river basin
munic_bhrd <- crop(munic, doce)
crs(munic_bhrd) <- "+init=EPSG:4326"

# rename to match column names with 
vals <- values %>%
  rename("NOMEMUNICP" = "municipality")

# create a list for data input
data <- list()

# get index mean for each municipality in shapefile and rasterize it
for(s in c("SSP126", "SSP370", "SSP585")){
  for(t in time){
    
    v <- filter(vals, scenario == s & time == t)
    
    raster <- merge(munic_bhrd, v, by.x = "NOMEMUNICP") %>%
      rasterize(., rast(., res = 0.005), field = "priority_mean")
    
    data[[paste(s, t, sep = "_")]] <- raster
    
  }
}

# raster to a dataframe list
data <- lapply(data, function(x){
  df <- as.data.frame(x) %>%
    bind_cols(., crds(x))
  return(df)
  }
  )

# merge lists into a single dataframe
data <- data %>%
  bind_rows(., .id = "scenario") %>%
  separate_wider_delim(scenario, "_", names = c("scenario", "time"), )

# plot
ggplot() +
  geom_raster(data = data, aes(x = x, y = y, fill = priority_mean)) +
  coord_equal() +
  geom_spatvector(data = munic_bhrd, color = "white", linewidth = 0.1, fill = NA) +
  geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_grid(scenario ~ time, switch = "y") +
  scale_fill_gradientn(colours = c("darkblue", "blue", "purple", "yellow", "orange", "red", "darkred"), name = "Restoration priority mean") +
#  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), name = "Restoration priority") +
#  theme_void() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

# export
ggsave("./results/figures/04_restoration_priority_munic.png",
       width = 10, height = 7, dpi = 300, bg = "white")


### 5.3.4. Figure: overall map with priority index per municipality  ####
# group municipality's index by scenario using the index sum across time
data_group <- data %>%
  group_by(scenario, x, y) %>%
  summarise(priority_sum = mean(priority_mean))

# plot
ggplot() +
  geom_raster(data = data_group, aes(x = x, y = y, fill = priority_sum)) +
  coord_equal() +
  geom_spatvector(data = munic_bhrd, color = "white", linewidth = 0.1, fill = NA) +
  geom_spatvector(data = doce, color = "black", linewidth = 0.5, fill = NA) +
  facet_wrap(scenario ~ .) +
  scale_fill_gradientn(colours = c("darkblue", "blue", "purple", "yellow", "orange", "red", "darkred"), name = "Restoration priority mean") +
  #  scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), name = "Restoration priority") +
  #  theme_void() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", strip.text = element_text(size = 13),  strip.background = element_rect(colour = "black", fill = "white"))

# export
ggsave("./results/figures/04_restoration_priority_munic_summary.png",
       width = 11, height = 5, dpi = 300, bg = "white")


