###-----------------------------------------------------------------------------
###  Title: DATA CLEANING
###
###  Description: 
###    TO CLEAN INCORRECT OCCURRENCE RECORDS FROM A TABLE
###
###  Obs: based on codes from '0_search_flora.R', '1_download_gbif_occurrences.R'
###    and '2_crop_gbif_occurrences.R'
###
###  Created by: Danielle Moreiro (@daniomoreira) // Feb 2021
###  Edited by: Alan Braz (@brazagm) // Set 2021
###-----------------------------------------------------------------------------

install.packages("devtools")
devtools::install_github("liibre/Rocc")

library(dplyr)
library(CoordinateCleaner)
library(flora)
library(magrittr) # for %T>% pipe
library(pacman)
library(parallel)
library(purrr)
library(readr)
library(raster)
library(rgbif) # for occ_download
library(rgdal)
library(Rocc)
library(taxize) # for get_gbifid_
library(tidyverse)


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## inputs
in_gbif_records <- "./processed_data/01-3_Gbif_records.csv"

## outputs
out_gbif_spatial_clean <- "./processed_data/02-1_Gbif_records_clean.csv"


# 1. SPATIAL CLEANING  ####-----------------------------------------------------
gbif <- read.csv(in_gbif_records)

## 1.1. Identify coordinates outside their reported country  ####
# CoordinateCleaner Package - Removes or flags mismatches between geographic coordinates and additional country information.
# usually this information is reliably reported with specimens). Such a mismatch can occur for example, if latitude and longitude are switched.
data_clean <- clean_coordinates(gbif, 
                                lon = "decimalLongitude",
                                lat = "decimalLatitude",)

# check number of mismatches by flag
summary(data_clean)

# check number of records with at least one spatial mismatch
data_clean %>% 
  count(.summary)

# filter only occurrences with no spatial mismatch
gbif_clean <- as.data.frame(data_clean) %>%
  subset(., .summary == "TRUE", select = names(data_clean)[-c(51:60)])

# to write the new table
#write.csv(data_clean3, "./data/registros/spp_Gualaxo/4_gbif_Gualaxo_amesul_clean.csv")


## 1.2. Dataframe cleanning  ####
# leaving only a relevant columns
gbif_clean <- subset(gbif_clean, select = c("gbifID",
                                            "kingdom", "order", "family",
                                            "genus", "species",
                                            "infraspecificEpithet",
                                            "taxonRank", "scientificName",
                                            "countryCode", "locality",
                                            "stateProvince",
                                            "decimalLatitude",
                                            "decimalLongitude",
                                            "elevation",
                                            "day", "month", "year",
                                            "institutionCode", "collectionCode",
                                            "catalogNumber", "recordNumber",
                                            "identifiedBy", "issue")
)

# rename "decimalLongitude" and "decimalLatitude" columns
gbif_clean <- gbif_clean %>%
  rename(lon = decimalLongitude,
         lat = decimalLatitude)


## 1.3. Select only species with more than 30 occurrences points  ####
# checking the frequency of occurrences for each species
n_records <- table(gbif_clean$species) %>%
  sort() %>%
  as.data.frame() %>%
  setNames(c("species", "n")) %>%
  mutate(valid = if_else(n >= 30, TRUE, FALSE))

rm_spp <- n_records[which(n_records$valid == FALSE), "species"] %>%
  as.character()

if(length(rm_spp) > 0){
  gbif_clean <- filter(gbif_clean, !(species %in% rm_spp))
}

# check the final list of species with more than 30 records
spp_list <- unique(gbif_clean$species)

# export the occurrence records data after spatial cleanning
write.csv(gbif_clean, out_gbif_spatial_clean)


# 2. TAXONOMICAL CLEANNING  ####------------------------------------------------
# import spatial cleaned data
data <- read.csv(out_gbif_spatial_clean)

# get species name after spatial cleaning
spp_names <- unique(data$species)

# check species name using 'Rocc' package
check <- check_string(unique(data$species))

# information of species status
table(check$speciesStatus)

# set species status to stay
# see Details in ?check_string for status results
stay <- c("possibly_ok", "variety", "subspecies", "name_w_authors")

# get species names to stay
spp_stay <- unlist(check[which(check$speciesStatus %in% stay), "species"])

# define species status to verify
# see Details in ?check_string for status results
verify <- c("not_Genus_epithet_format")

# get species names to verify
spp_verify <- unlist(check[which(check$speciesStatus %in% verify), "species"])



not_genus <- spp_names %>%
  filter(., check$speciesStatus %in% verify[1])
unique(not_genus$spp) # how many species need to check genus format
length(unique(not_genus$spp))

# To create a new column (check_ok).
check$check_ok <- "out"
head(check)

check$check_ok[check$speciesStatus %in% c(stay)] <- "ok"
check$check_ok[check$speciesStatus %in% "not_Genus_epithet_format"] <- "verify"

table(check$check_ok)

# now merging with object data2
tax_check <- cbind(data2, check[, -1])
head(tax_check)

# exporting data after check
write.csv(tax_check,
          file = "./data/registros/spp_Gualaxo/5_Gualaxo_occ_TaxonCheck.csv",
          row.names = FALSE)

## 2.1. Get synonyms from...  ####
# 1. Flora 2020
# 2. Kew - World Checklist of Vascular Plants
# 3. TNRS - Taxonomic Name Resolution Service

suggest_spp <- suggest_flora(spp_names)

search_spp <- suggest_spp$species %>%
  unique() %>%
  na.omit()

search_df <- data.frame(scientificName_search = search_spp,
                        search_id = 1:length(search_spp))

# 2. checking if the name exists in Brazilian Flora ####
?check_flora

flora_taxa <- list()
for (i in 1:length(search_spp)) {
  message(paste(i, "species"))
  flora_taxa[[i]] <- check_flora(search_spp[i],
                                 get_synonyms = FALSE,
                                 infraspecies = TRUE)
}

length(flora_taxa)

flora_taxa2 <- lapply(flora_taxa, function(x) x[1]$taxon)

flora_df <- bind_rows(flora_taxa2)

head(flora_df)

table(flora_df$taxonomicStatus)


# Changing column name of object search_df: scientificName_search to species
names(search_df)[names(search_df) == "scientificName_search"] <- "species"
head(search_df)

flora_df2 <- left_join(flora_df, search_df, by = "species")

# writing output
write.csv(flora_df2,
          file = "results/04_taxon_data_flora_check.csv",
          row.names = FALSE)

