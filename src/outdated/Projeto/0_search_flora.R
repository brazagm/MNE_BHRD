###-----------------------------------------------------------------------------
###  Title: OCCURRENCE RECORDS
###
###  Description: code for create a list of species with occurrence in the
###    states of Espírito Santo and Minas Gerais, Brazil
###
###  Created by: Danielle Moreiro (@daniomoreira) // 9 Feb 2021
###  Edited by: Alan Braz (@brazagm) // 23 Set 2021
###-----------------------------------------------------------------------------

install.packages("devtools")
devtools::install_github("liibre/Rocc")
library(Rocc)


# 1. LIST OF SPECIES NAME  ####
# Flora from the states of Espírito Santo and Minas Gerais, Brazil
# Search for species name in the List of Species of the Brazilian Flora 2020 database
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

## Identify state/province (est)
mg$est <- "MG"
es$est <- "ES"
list <- rbind(mg, es)

## Remove duplicated rows
dupl <- duplicated(list[, "id"])
list2 <- list[!dupl, ]

## Identify species occurring in both states
id_spp <- list$id[dupl]
list2[which(list2$id %in% id_spp), ]$est <- "MG/ES"

## Number of species occurring in both states
sum(dupl)

## Export species list
write.csv(list2, "./data/registros/lista_mg_es.csv", row.names = FALSE)

?search_flora

## PS: To get a list of species for the Doce river watersheld (BHRD), we need to (1) use the list "lista2",  (2) get records of the spp from GBIF, (3) then extract them only for the polygon of the BHRD.

# With the new list only for the BHRD, We will repeat the step 2, to get the records of these spp for South America, for niche modeling.
