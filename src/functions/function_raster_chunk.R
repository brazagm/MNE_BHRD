###
###  FUNCTION: RASTER_CHUNK
###
###  Description: create a function to work with large rasters through small
###     chunks. Based on https://strimas.com/post/processing-large-rasters-in-r/?utm_source=pocket_mylist
###
###  Created & Edited by: Alan Braz (@brazagm)
###
###  Observations:
###
###  Next tasks:
###

raster_chunk <- function(list){
  
# First part  ####
#n_values <- ncell(r) * nlayers(r)
# memory in mb
#mem_est <- 8 * n_values / 2^20
#mem_act <- as.integer(object.size(readAll(r))) / 2^20

#canProcessInMemory(r, verbose = TRUE)

# hack raster internal function
#cs_orig <- raster:::.chunk
#cs_hack <- function(x) getOption("rasterChunkSize")
#assignInNamespace(".chunk", cs_hack, ns = "raster")

# use 1 kb chunks
#rasterOptions(chunksize = 1000, todisk = TRUE)
#t_smallchunks <- system.time(calc(r, mean, na.rm = TRUE))

# undo the hack
#assignInNamespace(".chunk", cs_orig, ns = "raster")
#rasterOptions(default = TRUE)


# Second part ####
# file paths
f_in <- list
f_out <- tempfile(fileext = ".tif")

# input and output rasters
r_in <- stack(f_in)
r_out <- raster(r_in)

# blocks
b <- blockSize(r_in)
print(b)

# open files
r_in <- readStart(r_in)
r_out <- writeStart(r_out, filename = f_out)


# loop over blocks
for (i in seq_along(b$row)) {
  # read values for block
  # format is a matrix with rows the cells values and columns the layers
  v <- getValues(r_in, row = b$row[i], nrows = b$nrows[i])
  
  # FUNCTION -------------------------------------------------------------------
#  v[v == 0] <- NA # substitute 0 by NA
#  v[v != value] <- 0
#  v[v == value] <- 1
  
#  f <- (1000*10)/30 # factor multiplier from 30m to 10km
#  v <- raster::aggregate(v, fact = f, fun = mean, na.rm = TRUE)
  # end ------------------------------------------------------------------------
  
  # write to output file
  r_out <- writeValues(r_out, v, b$row[i])
}

# close files
r_out <- writeStop(r_out)
r_in <- readStop(r_in)


return(r_out)
}
