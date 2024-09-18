###
###  DISPERSAL PARAMETRIZATION
###
###  Description: estimate the dispersal parameters for simulations based on the
###     species' functional traits. 
###
###  Created & Edited by: Alan Braz (@brazagm)
###
###  Observations:
###
###  Next tasks:
###      1. Replace strings with abbreviated species names by their complete names
###

library(AICcmodavg)
library(nlme)
library(tidyverse)


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories for model results
project_repo <- "/home/alan/Documentos/Repositório/Projetos/INMA/MNE_BHRD"

## 1.1. Inputs and outputs  ####
in_spp <- file.path(project_repo, "results/niche_models")
in_list <- "./processed_data/01-1_species_list_revised.csv"
in_traits <- "./data/renova_traits_20240808.csv"


# 1. DISPERSAL PARAMETRIZATION  ####--------------------------------------------
## 1.1. Import and prepare data  ####
## import species names based on GBIF match (names for gathering occurrence records)
spp_names <- list.files(in_spp) %>%
  grep(".png|.csv", ., value = TRUE, invert = TRUE) %>%
  as.data.frame() %>%
  rename("species" = ".") %>%
  mutate(species = str_replace(species, "_", " "))

## import trait information and replace very small seed mass (= 0g) by 0.001g
## species names here area based on Flora do Brasil 2020
data <- read.csv(in_traits) %>%
  mutate(SeedMass_g = replace(SeedMass_g, SeedMass_g == 0.00, 0.001)) %>%
  rename("species" = "search.str") %>%
  mutate(sp_flora = species, .after = "id.flora")

## import and add life form information
form <- read.csv(in_list) %>%
  select(search.str, life.form) %>%
  rename("species" = "search.str", "Life_form" = "life.form") %>%
  unique()

## merge species names and trait measures
traits <- left_join(spp_names, data, by = "species") %>%
  left_join(., form, by = "species") %>%
  unique()

## species names not found in trait dataframe
subset(traits, is.na(sp_flora))

## IMPORTANT NOTE: we identified 6 missing species from original list in GBIF
## occurrence data. These species are correct synonyms and must be included
## manually here.
# species names based on GBIF
gbif_names <- c("Heliocarpus americanus", "Jupunba barnebyana", "Iochroma arborescens", "Piparea dentata", "Mollinedia umbellata", "Phyllostylon brasiliensis")

# species names based on Flora do Brasil 2020
flora_names <- c("Heliocarpus popayanensis", "Abarema barnebyana", "Acnistus arborescens", "Casearia commersoniana", "Mollinedia schottiana", "Phyllostylon brasiliense")

# synonyms dataframe
synonyms <- data.frame(species = gbif_names, sp_flora = flora_names)

# fill empty entries in trait dataframe using synonyms information
for(n in synonyms$species){
  
  name <- synonyms[which(synonyms$species == n), "sp_flora"]
  traits[which(traits$species == n), -c(1, 24)] <- data[which(data$species == name), -3]
  traits[which(traits$species == n), "Life_form"] <- form[which(form$species == name), "Life_form"]
  
  
}



## check again for species with without trait data
subset(traits, is.na(sp_flora))

## how many species has sufficient trait data
with_data <- subset(traits, !is.na(SeedMass_g) & !is.na(Dispersion), species)
nrow(with_data)




#data$SeedMass_g[data$SeedMass_g == 0.00] <- 0.001

# Corrigir o valor NA em Lafoensia glyptocarpa
# Valor aproximado extraído do estudo "CONSERVAÇÃO, GERMINAÇÃO E EFEITOS ALELOPÁTICOS DE Lafoensia Glyptocarpa KOEHNE", tese apresentada por LUIZA PAIVA SILVA DE MORAES a UNIVERSIDADE FEDERAL DE SÃO CARLOS, CENTRO DE CIÊNCIAS BIOLÓGICAS E DA SAÚDE PROGRAMA DE PÓS-GRADUAÇÃO EM ECOLOGIA E RECURSOS NATURAIS
#data[which(data$search.str == "Laf_gly"), "SeedMass_g"] <- 0.02

## 

traits <- traits %>%
  summarise(Species = species,
            GF = Life_form, # categorical growth form
            DS = Dispersion, # categorical dispersal syndrom
            SM = log(SeedMass_g/1000)) %>% # log10 transformed seed mass data (in mg)
  unique()

## see details in the dispeRsal.R instructions for the accepted nomenclature
## for data entry
traits <- traits %>%
  mutate(GF = recode(GF, "Árvore" = "tree", "Arbusto" = "shrub", "Arbusto|Árvore" = "shrub"),
         DS = recode(DS, "zoochoric" = "animal", "Zoochoric" = "animal", "Nonzoochoric" = "wind.none", "anemochoric" = "wind.special", "autochoric" = "ballistic"))

traits$GF <- as.factor(traits$GF)
traits$DS <- as.factor(traits$DS)


## 1.2. Dispersal parameters estimation  ####
load("/home/alan/Dropbox/Alan/INMA/Dados/MNE_BHRD/src/functions/dispeRsal/dispeRsal.rda")
source("./src/functions/dispeRsal/dispeRsal_rev.R")

disper <- dispeRsal2(traits, model = 4, CI = TRUE, random = FALSE, write.result = FALSE)

disper <- left_join(traits, disper$predictions[,-2], by = "Species")

write_csv(disper, "./results/Dispersal_dist.csv")

