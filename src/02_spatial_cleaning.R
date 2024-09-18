###-----------------------------------------------------------------------------
###  Title: SPATIAL CLEANING OF OCCURRENCE RECORDS 
###
###  Description:
###
###  Obs: 
###
###  Created by: Alan Braz (@brazagm) // Jul 2024
###-----------------------------------------------------------------------------

install.packages("CoordinateCleaner")
## load required libraries
library(CoordinateCleaner)
library(rWCVP)
library(terra)
library(tidyverse)

# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## local repositories
repo <- "/home/alan/Documentos/Repositório/GIS/shapes"

## inputs
in_spp <- "./processed_data/01-1_species_list_revised.csv"
in_match <- "./processed_data/01-2_Gbif_taxon_match.csv"
in_records <- "./processed_data/01-3_Gbif_records.csv"
in_world <- file.path(repo, "borders/ne_10m_admin_0_countries.shp")
in_bra <- file.path(repo, "borders/estados.shp")
in_biomes <- file.path(repo, "Biomas_250mil/lm_bioma_250.shp")

## outputs


# 1. IMPORT OCCURRENCE RECORDS  ####--------------------------------------------
## import records dataframe
records <- read.csv(in_records)

## import taxon match
match <- read.csv(in_match)

## import species details
species <- read.csv(in_spp)

#vect(., geom = c("decimalLongitude", "decimalLatitude"))

## get species names from the new and old list
spp_names <- sort(unique(records$species))
spp_list <- sort(unique(species$search.str))

## check the species in list but not in gbif records
miss <- setdiff(spp_list, spp_names)

## get synonyms
synonyms <- filter(match, original_sciname %in% miss) %>%
  select(species, original_sciname) %>%
  unique() %>%
  rename("gbif_species" = "species", "list_species" = "original_sciname")
  
## check species in gbif records but not in the list (or their synonyms)
remove <- setdiff(spp_names, c(spp_list, synonyms$gbif_species))

## remove them from gbif records
records <- filter(records, !species %in% remove)

## get the updated species names
spp_names <- sort(unique(records$species))


# 2. SPATIAL CLEANING OF SPECIES' RECORDS  ####---------------------------------
## set world border limits
wm <- borders("world", colour = "gray50", fill = "white")

## create a check list for the human verification of the automatic cleaning
check <- data.frame(species = spp_names, human_check = NA)


for(sp in spp_names[141:186]){
  
  # get species' records...
  dt <- subset(records, species == sp)
  
  # ... get spatial point records...
  points <- vect(dt, geom = c("decimalLongitude", "decimalLatitude"))
  
  # ... and get the spatial extent of the records
  e <- ext(points)
  
  # is a Brazilian endemic species?
  end <- if_else(subset(species, search.str == sp)$endemism == "Endemica", TRUE, FALSE)
  
  #
  
  native_range <- try({
    native_range <- wcvp_distribution(sp, taxon_rank = "species",
                                      introduced = FALSE, extinct = FALSE, 
                                      location_doubtful = FALSE) %>%
      vect(.)
  }, silent = TRUE)
  
  if(class(native_range) != "try-error"){
  
  native_range <- wcvp_distribution(sp, taxon_rank = "species",
                                    introduced = FALSE, extinct = FALSE, 
                                    location_doubtful = FALSE) %>%
    vect(.)
  
  native_records <- terra::extract(native_range, points) %>%
    select(id.y, LEVEL3_NAM, occurrence_type) %>%
    filter(!is.na(occurrence_type))
  
  native_records <- dt[native_records$id.y, "gbifID"]
  
  dt <- dt %>%
    mutate(.native_range = if_else(gbifID %in% native_records, TRUE, FALSE))
  }else{
    dt <- dt %>%
      mutate(.native_range = TRUE)
  }
  
  flags <- clean_coordinates(x = dt, lon = "decimalLongitude", lat = "decimalLatitude",
                             countries = "countryCode", species = "species",
                             tests = c("capitals", "centroids",
                                       "equal", "institutions", "zeros", "outliers"),
                             outliers_method = "quantile") # most test are on by default

  
  # ensure that records outside Brazil are not included for endemic species
  # according to Flora do Brazil 2020
#  if(end == TRUE){
#    if(length(flags[which(flags$.summary == TRUE & flags$countryCode != "BR"),]$.summary) != 0){
    
#    flags[which(flags$.summary == TRUE & flags$countryCode != "BR"),]$.summary <- FALSE
#    }
#  }
  valid_records <- flags %>%
    mutate(.valid = if_else(.summary == FALSE | .native_range == FALSE, FALSE, TRUE))
  
  map <- ggplot() +
    coord_fixed(xlim = c(e$xmin, e$xmax), ylim = c(e$ymin, e$ymax)) +
    wm +
    geom_point(data = subset(valid_records, .valid == TRUE),
               aes(x = decimalLongitude, y = decimalLatitude),
               colour = "orange",
               size = 1) +
    geom_point(data = subset(valid_records, .valid == FALSE),
               aes(x = decimalLongitude, y = decimalLatitude),
               colour = "darkred",
               size = 1) +
    ggtitle(sp) +
    theme_bw()
  
  print(map)
  
  if(class(native_range) != "try-error"){
    answer <- readline(paste0("Are the automatic spatial cleaning for ", sp," good enough? \n 1: YES \n 0: NO \n"))
  }else{
    answer <- readline(paste0("Are the automatic spatial cleaning for ", sp," good enough? \n IMPORTANT!!!: Data cleaning did not considered the species' geographic range \n 1: YES \n 0: NO \n"))
  }
  
  if(answer == "1"){
    check[which(check$species == sp), "human_check"] <- FALSE
  }else if(answer == "0"){
    check[which(check$species == sp), "human_check"] <- TRUE
    }
  
  if(sp == spp_names[1]){
    clean <- left_join(records, valid_records[, c(1, 51:60)], by = "gbifID")
  }else{
    clean <- rows_patch(clean, valid_records[, c(1, 51:60)], by = "gbifID")
  }
}
write_csv(clean, "./processed_data/02-1_Gbif_records_clean.csv")
write_csv(check, "./processed_data/02-2_check_records.csv")


# 3. REVISE SPECIES WITH DOUBTFUL OCCURRENCE RECORDS  ####----------------------

check_spp <- read.csv("./processed_data/02-2_check_records.csv") %>%
  filter(human_check == TRUE) %>%
  pull(species) %>%
  unique()

clean <- read.csv("./processed_data/02-1_Gbif_records_clean.csv") %>%
  add_column(.obs = NA)

sp <- check_spp[15]

# get species' records...
dt <- subset(clean, species == sp)

# ... get spatial point records...
points <- vect(dt, geom = c("decimalLongitude", "decimalLatitude"))

# ... and get the spatial extent of the records
e <- ext(points)

wm <- borders("world", colour = "gray50", fill = "white")

ggplot() +
  coord_fixed(xlim = c(e$xmin, e$xmax), ylim = c(e$ymin, e$ymax)) +
  wm +
  geom_point(data = subset(dt, .valid == TRUE),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange",
             size = 1) +
  geom_point(data = subset(dt, .valid == FALSE),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 1) +
  ggtitle(sp) +
  theme_bw()



# SPECIES ENDEMIC TO BRAZIL
# endemic to Brazil according to Flora do Brazil
endemic <- c("Albizia polycephala", "Enterolobium glaziovii")
clean <- clean %>%
  mutate(.valid = replace(.valid, species %in% endemic & countryCode != "BR" & .valid == TRUE, FALSE)) %>%
  mutate(.obs = replace(.obs, species %in% endemic, "Human-based decision: all Brazilian records is valid based on endemic status in the Flora do Brazil"))


# IGNORE OUTLIERS
#"Cabralea canjerana"
# "Myrsine coriacea"
omit_otl <- c("Cabralea canjerana", "Myrsine coriacea")
clean <- clean %>%
  mutate(.otl = replace(.otl, species %in% omit_otl, NA)) %>%
  mutate(.summary = replace(.summary, .val == TRUE & .equ == TRUE & .zer == TRUE & .cap == TRUE & .cen == TRUE & .inst == TRUE, TRUE)) %>%
  mutate(.valid = replace(.valid, .native_range == TRUE & .summary == TRUE, TRUE)) %>%
  mutate(.obs = replace(.obs, species %in% omit_otl, "Human-based decision: outlier filter was ignored"))


# IGNORE NATIVE RANGE FROM PLANTS OF THE WORLD

omit_native <- c("Jupunba barnebyana", "Solanum granulosoleprosum", "Stryphnodendron polyphyllum")
clean <- clean %>%
  mutate(.native_range = replace(.native_range, species %in% omit_native, NA)) %>%
#  mutate(.summary = replace(.summary, .val == TRUE & .equ == TRUE & .zer == TRUE & .cap == TRUE & .cen == TRUE & .inst == TRUE, TRUE)) %>%
  mutate(.valid = replace(.valid, species %in% omit_native & .summary == TRUE & .valid == FALSE, TRUE)) %>%
  mutate(.obs = replace(.obs, species %in% omit_native, "Human-based decision: native range was ignored"))


# IS CORRECT... (?)
correct <- c("Centrolobium robustum", "Handroanthus cristatus", "Hortia brasiliana", "Myrciaria glazioviana", "Platypodium elegans", "Schinus terebinthifolia")
clean <- clean %>%
  mutate(.obs = replace(.obs, species %in% correct, "Human-based decision: distribution apparently correct without change"))


#"Phyllostylon brasiliensis" TEM PROBLEMA TAXONOMICO mas dá pra manter
# "Ramisia brasiliensis" IDEM
doubtful <- c("Phyllostylon brasiliensis", "Ramisia brasiliensis")
clean <- clean %>%
  mutate(.obs = replace(.obs, species %in% doubtful, "Human-based decision: taxonomic revision is doubtful"))


# REGISTROS FORA DAS AMÉRICAS
# "Sapindus saponaria"
remove <- -120
clean <- clean %>%
  mutate(.valid = replace(.valid, decimalLongitude < remove & .valid == TRUE, FALSE)) %>%
  mutate(.obs = replace(.obs, decimalLongitude < remove & .valid == TRUE, "Human-based decision: records outside Americas were excluded based on longitude"))

write_csv(clean, "./processed_data/02-3_Gbif_records_clean_human-revised.csv")


# 4. NUMBER OF RECORDS  ####----------------------------------------------------
clean <- read.csv("./processed_data/02-3_Gbif_records_clean_human-revised.csv")

spp_names <- unique(clean$species)

n_records <- clean %>%
  group_by(species) %>%
  summarise(n_total = length(gbifID), n_valid = length(gbifID[.valid == TRUE])) %>%
  arrange(n_valid)

write_csv(n_records, "./results/02-1_number_of_valid_records.csv")















