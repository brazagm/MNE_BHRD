###-----------------------------------------------------------------------------
###  Title: SPECIES LIST, TAXONOMIC CLEANING AND OCCURRENCE RECORDS
###
###  Description: Run a taxonomic cleaning for the species list provided by the
###    Rede Rio Doce (Dec 2023). Only shrub and trees with occurrence records in
###    Espírito Santo and Minas Gerais states were select for analyses. Synonyms
###    was based on the Flora do Brasil 2020. Occurrence records of all species
###    were downloaded from GBIF.org
###
###  Obs: This script was based on codes from '0_search_flora.R',
###     '1_download_gbif_occurrences.R', and '2_crop_gbif_occurrences.R'
###     available at https://github.com/Projeto-BHRD-INMA/MNE_BHRD
###
###  Created by: Alan Braz (@brazagm) // Abr 2024
###-----------------------------------------------------------------------------

install.packages("devtools")
devtools::install_github("liibre/Rocc")

## load required libraries
library(devtools)
library(flora)
library(rgbif)
library(Rocc)
library(taxize)
library(tidyverse)

## install 'Rocc' package from Github repo if necessary
if(!require("Rocc")) devtools::install_github("liibre/Rocc")


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## inputs
in_list <- "./data/lista_espécies_rede_rio_doce.csv"

## outputs
out_list <- "./processed_data/01-1_species_list_revised.csv"
out_gbif_match <- "./processed_data/01-2_Gbif_taxon_match.csv"
out_gbif_data <- "./processed_data/GBIF_data"
out_gbif_records <- "./processed_data/01-3_Gbif_records.csv"


# 1. IMPORT AND EDIT SPECIES LIST  ####-----------------------------------------
## import and select only distinct scientific names
list <- read.csv(in_list) %>%
  select(NOME.CIENTÍFICO) %>%
  distinct()

## check species name using 'Rocc' package
check <- check_string(list$NOME.CIENTÍFICO)

## see how many names for each status
## ?check_string for status description in 'Details'
table(check$speciesStatus)

## get species name for checking
check[which(check$speciesStatus == "conferre"),]

## get species name for indet
check[which(check$speciesStatus == "indet"),]

## set species status to stay...
## see Details in ?check_string for status results
stay <- c("possibly_ok", "variety", "subspecies", "name_w_authors")

## ... and get species names to stay
spp_stay <- unlist(check[which(check$speciesStatus %in% stay), "species"])


# 2. GET SYNONYMS FROM FLORA DO BRASIL 2020  ####-------------------------------
## 2.1. Get synonyms suggestions  ####
## get synonyms of species names using 'Rocc'
suggest_spp <- suggest_flora(spp_stay)

## select only unique names and omit NAs
search_spp <- suggest_spp$species %>%
  unique() %>%
  na.omit()

## create a dataframe for results
search_df <- data.frame(scientificName_search = search_spp,
                        search_id = 1:length(search_spp))


## 2.2. Checking if the name exists in 'Flora do Brasil 2020'  ####
## get species details from 'Flora do Brasil' and select only trees and shrubs
## with terrestrial habitat and species level identification
spp_list <- get.taxa(search_spp, life.form = TRUE, habitat = TRUE, states = TRUE, establishment = TRUE, domain = TRUE, endemism = TRUE) %>%
  subset(., life.form %in% c("Árvore", "Arbusto|Árvore", "Arbusto") &
           taxon.rank == "species" & habitat == "Terrícola")

## select only species with occurrence records in the states of Espírito Santo
## or Minas Gerais AND within the Atlantic Forest or Cerrado
spp_list <- spp_list[grep("BR-ES|BR-MG", spp_list$occurrence),]
spp_list <- spp_list[grep("Cerrado|Mata Atlântica", spp_list$domain),]

## export species list
write_csv(spp_list, out_list)


# 3. GET SPECIES OCCURRENCE RECORDS FROM GBIF  ####-----------------------------
## 3.1. Provide credential for GBIF.org  ####
## the important part here is to use rgbif::occ_download with pred_in and to fill in your gbif credentials.
## fill in your gbif.org credentials. You will need to create an account at gbif if you don't have it.
user <- rstudioapi::askForPassword("Your 'gbif.org' username") # your gbif.org username
pwd <- rstudioapi::askForPassword("Your 'gbif.org' password") # your gbif.org password
email <- rstudioapi::askForPassword("Your 'gbif.org' e-mail") # your email


## 3.2. Get species ID from GBIF.org  ####
## import species list and get only unique names using 'search.str' column
oc <- read.csv(out_list) %>%
  select(search.str) %>%
  unique()

## IMPORTANT NOTE: we identified 9 missing species from original list in GBIF
## occurrence data. These species are correct synonyms and must be included
## manually here.
synonyms <- c("Abarema barnebyana", "Acnistus arborescens", "Casearia commersoniana", "Handroanthus impetiginosus", "Heliocarpus popayanensis", "Mollinedia schottiana", "Phyllostylon brasiliense", "Psidium cattleyanum", "Trema micrantha")


## get taxon keys from GBIF
## gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
gbif_taxon_keys <- oc %>% # for an file with a list of spp names
  pull(search.str) %>% # specify the column from the list
  taxize::get_gbifid_(method = "backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_csv(file = out_gbif_match) %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status != "DOUBTFUL") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys


## 3.3. Get species occurrence records from GBIF.org  ####
## use matched gbif_taxon_keys from above
## !!very important here to use pred_in!!
x <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN')),
  #pred("geometry","POLYGON((-43.86 -17.57, -43.88 -21.49, -39.79 -19.86, -39.46 -17.98, -43.86
  #     -17.57))"),
  #pred("country", "BR"),
  #pred("continent", "South America"),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_gte("year", 1950),
  format = "SIMPLE_CSV",
  user = user, pwd = pwd, email = email
)

## check download processing status
## wait if 'status: running'
occ_download_wait(x[[1]])

## create a file directory for GBIF zipfile
ifelse(dir.exists(out_gbif_data), "GBIF directory already exists!",
       dir.create(out_gbif_data))

## after finish the processing (= 'status: succeeded'),
## download occurrence records from GBIF.org
d <- occ_download_get(x[1], path = out_gbif_data, overwrite = TRUE) %>%
  occ_download_import(path = out_gbif_data)

## export occurrence records
write_csv(d, out_gbif_records)

