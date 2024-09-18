################################################################################
###                                                                          ###
###             TO CROP SPECIES RECORDS POINTS TO A SHAPEFILE                ###
###                                                                          ###
###                                                                          ###
###                Created by Danielle de O. Moreira                         ###
###                         09 feb 2021                                      ###
###                                                                          ###
################################################################################

library(data.table)
library(raster)
library(rgdal)
library(sp)


# 1. IMPORT OCCURRENCE RECORDS DATA  ####
# read the species records points
# I use the command data.table::fread() because read.delim() or read.csv() are not reading all rows
spp.table <- fread("./data/registros/spp_Gualaxo/2_Gbif_occ_Gualaxo.csv")

# check the table dimensions
dim(spp.table)

# check the first table rows
head(spp.table)

# check the column names
names(spp.table)
View(spp.table)

# check the number of species
unique(spp.table$species)


# 2. OCCURRENCE RECORDS AS SPATIAL POINTS  ####
# Converting data.frame into a SpatialPointsDataFrame for spatial operations
# note that the lon and lat columns are in columns 23 and 22, respectively
# Get long and lat from your data.frame. Make sure that the order is in lon/lat.
xy <- spp.table[, c("decimalLongitude", "decimalLatitude")]

# convert
spp.shp <- SpatialPointsDataFrame(coords = xy, 
                                          data = spp.table,
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))


# 3. GET RECORDS FROM SOUTH AMERICA  ####
# To clip species records to a specific shapefile; in this case, the South America
## Read the shapefile
amesul <- readOGR("./data/shape/amesul.shp") #loading shapefile

# check CRS information 
crs(amesul)
#proj4string(amesul)

## to assign the crs information  
## reproject shp to geographic system WGS 1984
crs(amesul) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#proj4string(amesul) <- CRS("+proj=longlat +datum=WGS84") #to assign a CRS

proj4string(amesul)

## Plot shp
plot(spp.shp, cex = .1)
plot(amesul, add=TRUE)
dev.off()

# Now let's crop the points into the shapefile limits we want to. In this case, South America
spp.shp.ame <- crop(spp.shp, amesul) #cut by file

dim(spp.shp.ame)

##Plot spatial objects
plot(spp.shp.ame, cex = .3)
plot(amesul, add=TRUE)

#Save table of species for South America
write.csv(spp.shp.ame, "./data/registros/spp_Gualaxo/3_gbif_Gualaxo_amesul.csv")


