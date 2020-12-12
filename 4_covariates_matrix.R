# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~ Extracting selected covariates to matrices ~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# In this script, generate a regular 5x5km matrix containing the values of selected predictors across Africa


rm(list=ls())
gc()

## Set up the pathways to access to the source of data (spatial and others)
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_raster <- paste(path_wd, "/Raster/human_predictors/cont", sep = "") 
path_BG <- paste(path_wd, "/Shapefiles/BG_points", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

# 2. Load the packages that we need for this script
source(paste(path_scripts,"ipak.R", sep = "/"))
packages <- c("tidyverse","reshape2", "sf", "raster","rgdal","proj4")
ipak(packages)

# "stats","readxl","foreign", "dismo", "ggplot2","cowplot", "usdm","kernlab","ks","sm","biomod2","rJava","pROC","pryr", "RStoolbox","vegan","psycho","corrplot", "automap"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~ Data preparation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Load the Africa map with country boundaries
AFRO_map <- file.path(path_shapefiles,"SALB.shp")
AFRO_map <- st_read(AFRO_map)

occurrences_PCR <- readOGR(dsn=paste(path_shapefiles,"/BU_occurrences",sep=""), layer="BU_occurrences_PCR")

# Load the selected predictors (explanatory variables) and extract values for surveyed locations
AFRO_SALB <- as(AFRO_map, 'Spatial') # coerce sf object to spatial object
PCSR <- AFRO_SALB@proj4string

## Rasterize Africa map and then produce a coordinate matrix at 5km resolution.
raster.map.5km <- raster(AFRO_SALB)
res(raster.map.5km) <- 5000 # 5000 m resolution
raster.map.5km <- rasterize(AFRO_SALB, raster.map.5km)
# plot(raster.map.5km)
pred_coords_5km <- rasterToPoints(raster.map.5km)[,c(1, 2)]
pred_coords_5km <- as.data.frame(pred_coords_5km)

setwd(path_raster)

### Get the path to the raster files
raster.files <- list.files(path_raster,  
                           pattern="*tif$", full.names=TRUE) 

stopifnot(length(raster.files)>0)

### Get the raster names
raster.names <- list.files(path_raster, 
                           pattern="*tif$", 
                           full.names=FALSE)

raster.names <- c(unlist(lapply(strsplit(raster.names,"[.]"), FUN=function(x) { x[1] })))   
raster.names

# Generate prediction surface with covariate data at a 5km resolution
## Generate a SpatialPoint object with the coordinates matrix from Africa
pred_5km.sp <- SpatialPoints(pred_coords_5km[,c("x","y")],
                             proj4string = PCSR)

## 6.2. Extract the values from covariates for the Africa matrix at 5km resolution
for(i in 1:length(raster.files)){
  raster <- raster(raster.files[i])
  c <- data.frame(raster::extract(raster, 
                                  pred_5km.sp,
                                  method = "bilinear"))
  colnames(c) <- raster.names[i]
  pred_coords_5km <- cbind(pred_coords_5km,c)
}

names(pred_coords_5km)

# rm(pred_coords_5km)
rm(i, c, raster)

pred_coords_5kmm<- as.matrix(pred_coords_5km)
nrow(pred_coords_5km)

### Save final matrix in working directory
setwd(path_wd)
saveRDS(pred_coords_5km,"covariates_5km.rds")



