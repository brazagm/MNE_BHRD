library(modleR)
library(dplyr)
library(raster)
library(progress)
library(modleR)
library(foreach)
library(rgeos)
library(doParallel)
library(purrr)
library(furrr)

getwd()

# reading species data, only names
#target_species <- read.csv("./species/lista_species.csv",
#                           stringsAsFactors = FALSE) %>%
#                  pull(sp)

target_species <- c("Acnistus_arborescens",
                    "Piptadenia_gonoacantha",
                    "Lecythis_pisonis"
                    #"Paubrasilia_echinata",
                    #"Albizia_niopoides",
                    #"Handroanthus_heptaphyllus",
                    #"Zeyheria_tuberculosa",
                    #"Eremanthus_erythropappus",
                    #"Balfourodendron_riedelianum",
                    #"Citharexylum_myrianthum",
                    #"Maclura_tinctoriatinctoria",
                    #"Senna_multijugamultijuga",
                    #"Dalbergia_nigra",
                    #"Schizolobium_parahyba",
                    #"Anadenanthera_peregrinafalcata"
)

# reading occurrence table of all species
# in case that your occurrence table have more species than you want to project

occs <- read.csv2("./azure/species/6_Gualaxo_ModleR_out.csv", sep=",",
                  header=TRUE,  dec=".") %>%
  filter(sp %in% target_species)


occs<-occs[,2:4]
data_list <- split(occs, occs$sp)
names(data_list) #check names
species <- names(data_list)
species



wc <- list.files("./azure/raster/pca/pres/",
                 pattern = "tif$",
                 full.names = TRUE) %>%
  stack()


# detectCores() retorna o número de núcleos da sua máquina
# se preferir, pode reduzir em 1 núcleo de uso - 'workers = detectCores() - 1'
plan(multisession, workers = detectCores())

data_list %>%
  furrr::future_map2(.x = .,
              .y = as.list(species),
              ~ setup_sdmdata(species_name = .y,
                              models_dir = "./azure/modler_gualaxo/",
                              occurrences = .x,
                              predictors = wc,
                              seed = 123,
                              partition_type = "crossvalidation",
                              cv_n = 1,
                              cv_partitions = 10,
                              buffer_type = "mean",
                              write_buffer = TRUE,
                              png_sdmdata = TRUE,
                              n_back = 10000,
                              clean_dupl = TRUE,
                              clean_uni = TRUE,
                              clean_nas = TRUE,
                              geo_filt = TRUE,
                              geo_filt_dist = 0.1,
                              select_variables = FALSE
                              )
  )

species %>%
  as.list(.) %>%
  furrr::future_map(~ do_many(species_name = .,
                              predictors = wc,
                              models_dir = "./azure/modler_gualaxo/",
                              project_model = TRUE,
                              proj_data_folder = "./azure/raster/pca/proj/",
                              dismo_threshold = "spec_sens",
                              png_partitions = TRUE,
                              bioclim = TRUE,
                              #brt = TRUE,
                              #glm = TRUE,
                              maxnet = TRUE,
                              svmk = TRUE,
                              #rf = TRUE,
                              #equalize = TRUE,
                              write_bin_cut = TRUE)
  )


for(i in c("present", "mpi_126", "mpi_585")){
species %>%
  as.list(.) %>%
  furrr::future_map(~ final_model(species_name = .,
                                  algorithms = NULL,
                                  consensus_level = 0.5,
                                  models_dir = "./azure/modler_gualaxo/",
                                  mean_th_par = c("spec_sens"),
                                  which_models = c("raw_mean",
                                                   "bin_mean",
                                                   "bin_consensus"),
                                  proj_dir = i,
                                  uncertainty = TRUE,
                                  png_final = TRUE,
                                  overwrite = TRUE)
  )
}

for(i in c("present", "mpi_126", "mpi_585")){
data_list %>%
  furrr::future_map2(.x = .,
                     .y = as.list(species),
                     ~ ensemble_model(species_name = .y,
                                      occurrences = .x,
                                      models_dir = "./azure/modler_gualaxo/",
                                      performance_metric = "TSSmax",
                                      which_ensemble = c("weighted_average", "consensus", "best"),
                                      which_final = c("raw_mean",
                                                      "bin_mean",
                                                      "bin_consensus"),
                                      consensus_level = 0.5,
                                      proj_dir = i,
                                      png_ensemble = TRUE,
                                      uncertainty = TRUE,
                                      overwrite = TRUE)
  )
}
