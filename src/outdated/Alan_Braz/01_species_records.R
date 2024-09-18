###-----------------------------------------------------------------------------
###  Title: SPECIES RECORDS
###
###  Description: identify and list the species with occurrence in the states
###    of Espírito Santo and Minas Gerais, Brazil, and download occurrence
###    records of species from GBIF database.
###    GBIF tutorial is available at:
###    https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/
###
###  Obs: based on codes from '0_search_flora.R', '1_download_gbif_occurrences.R'
###    and '2_crop_gbif_occurrences.R'
###
###  Created by: Danielle Moreiro (@daniomoreira) // Feb 2021
###  Edited by: Alan Braz (@brazagm) // Set 2021
###-----------------------------------------------------------------------------

install.packages("devtools")
devtools::install_github("liibre/Rocc")

library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(Rocc)
library(taxize) # for get_gbifid_
library(tidyverse)


# 0. INPUTS & OUTPUTS  ####-----------------------------------------------------
## inputs
in_list_gualaxo <- "./data/01_Species_list_Gualaxo.csv"

## outputs
out_list_states <- "./processed_data/01-1_Species_list_states.csv"
out_gbif_match <- "./processed_data/01-2_Gbif_taxon_match.csv"
out_gbif_records <- "./processed_data/01-3_Gbif_records.csv"


# 1. LIST OF SPECIES NAME  ####-------------------------------------------------
## 1.1. Get flora from the states of Espírito Santo and Minas Gerais, Brazil  ####
# search for species name in the List of Species of the Brazilian Flora 2020 database
mg <- search_flora(
  domain = NULL,
  stateProvince = "MG",
  life_form = "Árvore"
)

es <- search_flora(
  domain = NULL,
  stateProvince = "ES",
  life_form = "Árvore"
)

# identify state/province (est)
mg$est <- "MG"
es$est <- "ES"
list <- rbind(mg, es)

# remove duplicated rows
dupl <- duplicated(list[, "id"])
list2 <- list[!dupl, ]

# identify species occurring in both states
id_spp <- list$id[dupl]
list2[which(list2$id %in% id_spp), ]$est <- "MG/ES"

# number of species occurring in both states
sum(dupl)

# export species list
write.csv(list2, out_list_states, row.names = FALSE)

?search_flora

## PS: To get a list of species for the Doce river watersheld (BHRD), we need to (1) use the list "lista2",  (2) get records of the spp from GBIF, (3) then extract them only for the polygon of the BHRD.

# With the new list only for the BHRD, We will repeat the step 2, to get the records of these spp for South America, for niche modeling.


# 2. GET SPECIES OCCURRENCE RECORDS FROM GBIF  ####-----------------------------
## 2.1. Provide credential for GBIF.org  ####
## The important part here is to use rgbif::occ_download with pred_in and to fill in your gbif credentials.
## fill in your gbif.org credentials. You will need to create an account at gbif if you don't have it.
user <- rstudioapi::askForPassword("Your 'gbif.org' username") # your gbif.org username
pwd <- rstudioapi::askForPassword("Your 'gbif.org' password") # your gbif.org password
email <- rstudioapi::askForPassword("Your 'gbif.org' e-mail") # your email


## 2.2. Get species ID from GBIF.org  ####
oc <- read.csv(in_list_gualaxo)
oc$species <- as.factor(oc$species)
names(oc)

gbif_taxon_keys <- oc %>% # for an file with a list of spp names
  pull(species) %>% # specify the column from the list
  taxize::get_gbifid_(method = "backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_csv(path = out_gbif_match) %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys

# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!


## 2.3. Get species occurrence records from GBIF.org  ####
# use matched gbif_taxon_keys from above
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

d <- occ_download_get(x[1]) %>%
  occ_download_import()

write_csv(d, out_gbif_records)


