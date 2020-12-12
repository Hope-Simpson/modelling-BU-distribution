## https://rdrr.io/cran/spatialEco/man/sp.kde.html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~ Generate background data for BU models ~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#  In this script, use the spatialeco package to generate a kernel density surface around occurrence points
#  then extract background points from this density surface

## Set up the pathways to access to the data 
rm(list = ls())
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_raster <- paste(path_wd, "/Raster/human_predictors/cont", sep = "") 
path_raster_outputs <- paste(path_wd, "/Raster", sep = "") 
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_BG <- paste(path_wd, "/Shapefiles/BG_points", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

setwd(path_wd)
# Load the packages that we need for this project
source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("sp", "sf", "rgdal", "raster", "fasterize", "st", "tmaptools", "spatialEco", "dplyr")
ipak(packages)

## Import an African map
setwd(path_shapefiles)
Africa_ADM0 <- readOGR(dsn=path_shapefiles, layer="salb")

## read occurrence dataset
AFRO_MU <- file.path(path_shapefiles, "/BU_occurrences/env_occurrences.shp")
AFRO_MU <- st_read(AFRO_MU)

names(AFRO_MU)

AFRO_MU <- AFRO_MU[,c("N_postv", "MU_cnfr" )]

colnames(AFRO_MU) 
AFRO_MU <- AFRO_MU %>% rename(MU = N_postv,
                              weight = MU_cnfr)

AFRO_MU$MU <- 1

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

# Generate random bakground points with same bias as occurrence data (Fitzpatrick et al., 2013, Elith 2009)
# Creating a surface for the study area with a KDE where higher probability values are assigned to areas with more occurrence points. 

# create a vector mask of the continent
Mask <- AnPET
values(Mask) <- 1
Mask <- mask(Mask, AnPET)
# plot(Mask)

# convert occurrences to sp 
occurrences <- AFRO_MU
nrow(occurrences)
occurrences <- st_transform(occurrences, crs = PCSR)
pres_abs_sp <- as_Spatial(occurrences)

# Identify which cells have presences
# Make sure there are no duplicates in the cell IDs
bias <- cellFromXY(Mask, pres_abs_sp)
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

# generate KDE raster around points with bandwidth of 150km 
pt.kde <- sp.kde(x = dfxy, y = samps, bw = 150000, newdata = Mask, standardize = TRUE, nr = 5000, 
                 scale.factor = 100000000)

# resample to get same resolution
KDEr <- resample(pt.kde, Mask)

# save raster
setwd(path_raster_outputs)
writeRaster(KDEr, "KDEraster_env.tif", overwrite=TRUE)

# create a buffer 10 km around occurrences (to prevent contamination)
nrow(occurrences)

occ_buff <- st_buffer(occurrences, 10000)
buffer10 <- st_geometry(occ_buff)
buffer10 <- st_union(buffer10)
setwd(path_shapefiles)

st_write(buffer10, "buffer10_env.shp", delete_layer = T)

buffer10 <- st_read("buffer10_env.shp")
buffer10 <- st_transform(buffer10, st_crs(AnPET))

## rasterise buffer areas
buffer_raster <- fasterize(buffer10, Mask, field = NULL, fun = "last", background = NA_real_,
                           by = NULL)

# convert nas to 0s and 1s to NA so buffer areas have NA value
buffer_raster[is.na(buffer_raster)] <- 0
buffer_raster[raster::values(buffer_raster)==1] <- NA
buffer_raster[raster::values(buffer_raster)==0] <- 1
buffer_raster <- resample(buffer_raster, Mask)
plot(buffer_raster)

setwd(path_raster_outputs)
writeRaster(buffer_raster, "env_10kbuffer_raster.tif", format= "GTiff", overwrite=T)

buffer_raster <- raster("env_10kbuffer_raster.tif")

# 3.7. Exclude areas outside the study area or within 10km of occurrence points
# plot(st_geometry(occurrences), pch=19, cex=0.5)
# plot(st_geometry(buffer10), add = T)

## multiply rasters together so sample area is restricted to IN africa, OUT 10k buffer, IN KDE
KDEraster <- KDEr * Mask * buffer_raster

writeRaster(KDEraster, "KDEraster_env.tif", format= "GTiff", overwrite=T)
KDEraster <- raster("KDEraster_env.tif")
plot(KDEraster)

# Potential background points will be sampled from a list of point data
# extract points from the KDEraster 
KDEpts <- rasterToPoints(KDEraster, spatial = F)
head(KDEpts)

# sample for background data
set.seed(123)
BG <- KDEpts[sample(seq(1:nrow(KDEpts)), size = nrow(occurrences), replace = T,
                    prob = KDEpts[,"KDEraster_env"]),1:2]

BG <- as.data.frame(BG) 
# Assign these as negative points
BG$MU <- 0               

nrow(BG)
BG_sf = st_as_sf(BG, coords = c("x", "y"), crs = st_crs(Africa_ADM0))


plot(BG_sf, add=T)

setwd(path_BG)
st_write(BG_sf, "KDE_BG_env.shp", delete_layer  = TRUE)
BG_sf <- st_read("KDE_BG_env.shp")



