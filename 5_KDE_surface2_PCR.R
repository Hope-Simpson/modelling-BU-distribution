## https://rdrr.io/cran/spatialEco/man/sp.kde.html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~ Generate background data for BU models ~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#  In this script, use the spatialeco package to generate a kernel density surface around occurrence points
#  then extract background points from this density surface
gc()

rm(list=ls())

## Set up the pathways to access to the data 
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_raster <- paste(path_wd, "/Raster/human_predictors/cont", sep = "") 
path_raster_outputs <- paste(path_wd, "/Raster", sep = "") 
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_BG <- paste(path_wd, "/Shapefiles/BG_points", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("sp", "sf", "rgdal", "raster", "fasterize", "st", "tmaptools", "spatialEco", "dplyr")
ipak(packages)

## Import an African map
setwd(path_shapefiles)
Africa_ADM0 <- readOGR(dsn=path_shapefiles, layer="salb")
## 
PCSR <- Africa_ADM0@proj4string

setwd(path_BG)

## read occurrence dataset
AFRO_BU <- file.path(path_shapefiles, "/BU_occurrences/BU_occurrences_PCR.shp")
AFRO_BU <- st_read(AFRO_BU)

names(AFRO_BU)
nrow(AFRO_BU)

# keep only number and  weight
AFRO_BU <- AFRO_BU[,c("number", "weight" )]
AFRO_BU <- st_transform(AFRO_BU, crs = st_crs(Africa_ADM0))
colnames(AFRO_BU) 
AFRO_BU <- AFRO_BU %>% rename(BU = number)
AFRO_BU$BU <- 1
occurrences <- AFRO_BU
nrow(occurrences)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Load any one of the covariates prepared for this modelling exercise
# just to get the correct geographical extent and resolution
setwd(path_raster)
raster.files <- list.files(path_raster, 
                           pattern="*tif$", full.names=TRUE) 

AnPET <- raster((raster.files[1]))
# plot(AnPET)
PCSR <- st_crs(AnPET)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~ Generate background points ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Generate random bakground points with same bias as occurrence data 
# (Phillips et al. 2009 esa; Fitzpatrick et al. 2013, esa; Vollering 2019)

# 3.1. Creating a surface for the study area with a KDE where higher probability values are assigned to areas with more locality points. 

# creat a vector mask of the continent
Mask <- AnPET
values(Mask) <- 1

Mask <- mask(Mask, AnPET)
# plot(Mask)

# convert occurrences to sp 
pres_abs_sp <- as_Spatial(occurrences)
bias <- cellFromXY(Mask, pres_abs_sp)

# . Make sure there are no duplicates in the cell IDs
cells <- unique(sort(bias))

#  Get XY coordinates for the cells, from the mask
kernelXY <- xyFromCell(Mask, cells)

# convert to df
dfxy <- as.data.frame(kernelXY)

# Count number of presences by cell
samps <- as.numeric(table(bias))

# make spatial object
coordinates(dfxy) <- ~x+y
length(dfxy)

# generate KDE raster around points with bandwidth of 150km (arbitrary choice)
pt.kde <- sp.kde(x = dfxy, y = samps, bw = 150000, newdata = Mask, standardize = TRUE, nr = 5000, scale.factor = 100000000)

# resample to get same resolution
KDEr <- resample(pt.kde, Mask)

# save raster
setwd(path_raster_outputs)
writeRaster(KDEr, "KDEraster_PCR.tif", overwrite=TRUE)

# create a buffer 10 km around occurrences (to prevent contamination)
occ_buff <- st_buffer(occurrences, 10000)
buffer10 <- st_geometry(occ_buff)
buffer10 <- st_union(buffer10)
setwd(path_shapefiles)

st_write(buffer10, "buffer10_PCR.shp", delete_layer = T)

# buffer10 <- readOGR(dsn=path_shapefiles, layer="10kbuffer")
# buffer10 <- spTransform(buffer10, PCS)

buffer10 <- st_read("buffer10_PCR.shp")

buffer10 <- st_transform(buffer10, st_crs(AnPET))

## rasterise buffer areas
buffer_raster <- fasterize(buffer10, Mask, field = NULL, fun = "last", background = NA_real_,
                           by = NULL)

# convert nas to 0s and 1s to NA so buffer areas have NA value
 buffer_raster[is.na(buffer_raster)] <- 0
 buffer_raster[raster::values(buffer_raster)==1] <- NA
 buffer_raster[raster::values(buffer_raster)==0] <- 1
 buffer_raster <- resample(buffer_raster, Mask)
# plot(buffer_raster)
setwd(path_raster_outputs)
writeRaster(buffer_raster, "PCR_10kbuffer_raster.tif", format= "GTiff", overwrite=T)

buffer_raster <- raster("PCR_10kbuffer_raster.tif")

# 3.7. Exclude areas outside the study area or within 10km of occurrence points
# plot(st_geometry(occurrences), pch=19, cex=0.5)
# plot(st_geometry(buffer10), add = T)
head(Mask)

## multiply rasters together so sample area is restricted to IN africa, OUT 10k buffer, IN KDE
KDEraster <- KDEr * Mask * buffer_raster

writeRaster(KDEraster, "KDEraster_PCR.tif", format= "GTiff", overwrite=T)
KDEraster <- raster("KDEraster_PCR.tif")
# plot(KDEraster)

# 3.8. Potential background points will be sampled from a list of point data
# extract points from the KDEraster 
KDEpts <- rasterToPoints(KDEraster, spatial = F)
head(KDEpts)

# 3.9.   sample for background data
set.seed(123)
BG <- KDEpts[sample(seq(1:nrow(KDEpts)), size = nrow(occurrences), replace = T,
                    prob = KDEpts[,"KDEraster_PCR"]),1:2]

BG <- as.data.frame(BG) 
# Assign these as negative points
BG$BU <- 0              

BG_sf = st_as_sf(BG, coords = c("x", "y"), crs = st_crs(Africa_ADM0))

plot(KDEraster)
plot(BG_sf, add=T)

setwd(path_BG)
st_write(BG_sf, "KDE_BG_PCR.shp", delete_layer  = TRUE)
nrow(BG_sf)



