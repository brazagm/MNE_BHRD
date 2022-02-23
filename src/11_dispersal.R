library(MigClim)
library(remotes)

remotes::install_version("SDMTools", "1.1-221")
install_github("robinengler/MigClim")

doce <- shapefile("/home/alan/Documentos/azure_MNE_BHRD/azure/shapes/bhrd_wgs84_dissolv.shp")
proj4string(doce) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

init <- raster("/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo/Acnistus_arborescens/present/ensemble/Acnistus_arborescens_ensemble_0.5_consensus.tif")

hsuit <- raster("/home/alan/Documentos/azure_MNE_BHRD/azure/modler_gualaxo/Acnistus_arborescens/mpi_585/ensemble/Acnistus_arborescens_TSSmax_ensemble_weighted_average.tif")

init_crop <- crop(mask(init, doce), extent(doce))
hsuit_crop <- crop(mask(hsuit, doce), extent(doce))

MigClim.migrate(iniDist=init_crop, hsMap=hsuit_crop*1000, 
                rcThreshold=0, envChgSteps=1,
                dispSteps=1, dispKernel=c(1.0,1.0), barrier="", barrierType="strong", iniMatAge=1,
                propaguleProd=c(1.0), lddFreq=0.0, lddMinDist=NULL, lddMaxDist=NULL,
                simulName="MigClimTest", replicateNb=1, overWrite=FALSE, testMode=FALSE,
                fullOutput=FALSE, keepTempFiles=FALSE)
