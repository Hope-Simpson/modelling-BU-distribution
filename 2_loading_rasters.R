# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~ BU distribution modelling across Africa ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ 2. Preparing covariates matrix              ~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# read and stack potential predictors which have been prepared within an extent which includes most of the data 
# extract the values of predictors to all occurrence locations and save the covariates matrices

.rs.unloadPackage("tidyr")

rm(list=ls())

## 1. Set up the pathways to access to the source of data (spatial and others)
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_raster <- paste(path_wd, "/Raster/clipped", sep = "") 
path_outputs <- paste(path_wd, "/Outputs/covariates", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 
setwd(path_wd)

## 1. Call a function to install and load multiple R packages.
source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("sp", "raster","rgdal","proj4","ggplot2","gridExtra",
              "readxl", "magrittr")
ipak(packages)

#read in datasets
setwd(path_shapefiles)
occurrences_PCR <- readOGR(dsn=paste(path_shapefiles,"/BU_occurrences",sep=""), layer="BU_occurrences_PCR")
env_occurrences <- readOGR(dsn=paste(path_shapefiles,"/BU_occurrences",sep=""), layer="env_occurrences")

# read frame containing most of data 
frame <- readOGR(dsn=paste(path_shapefiles), layer="Frame_WestAfrica")

## 2. Rasterize Framework for extraction of raster datasets
raster.map.frame <- raster(frame)
res(raster.map.frame) <- 5000 # 5000 m resolution
values(raster.map.frame) <- 1 
raster.map.frame <- rasterize(frame, raster.map.frame)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. Working with the covariates that will be used in the analysis

# 3.1 Load the covariates that have been clipped to the West Africa extent
### continuous
setwd(path_raster)
list.files(pattern="*tif$")
raster.files <- list.files(path_raster, 
                           pattern="*tif$", full.names=TRUE) 

r <- raster((raster.files[2]))
PCSR <- r@crs

raster.map.frame.p <- projectRaster(raster.map.frame, crs = PCSR)

rmf <- spTransform(frame, CRS("+proj=longlat +datum=WGS84"))

# 3.2 Iterate through the raster and adjust extension and crop for the framework area, and finally pack them all
#     as a raster Stack object.
 Covariates <- stack()
 for (j in 1:length(raster.files)){
  i <- raster(raster.files[j])
  i[i < 0] <- NA
  i <- extend(i, raster.map.frame.p)
  i <- raster::resample(i, raster.map.frame.p, method = 'bilinear')
  i <- crop(i, raster.map.frame.p)
  Covariates <- stack(Covariates,i)
  }

names(Covariates)

# save raster stack as RDS
setwd(path_raster)
saveRDS(Covariates, file = "stack_cov1.rds")

# read the covariates RDS
Covariates <- readRDS("stack_cov1.rds")
names(Covariates)

######### EXTRACT RASTER VALUES TO POINTS - HUMAN 
BU.cov <-  raster::extract(Covariates, occurrences_PCR)
######### EXTRACT RASTER VALUES TO POINTS - ENVIRONMENTAL
BU.env_Cov <-  raster::extract(Covariates, env_occurrences)

# write covariates matrices to folder
setwd(path_outputs)

write.csv(BU.cov, file = "BU_cov_PCR.csv", row.names = F)
write.csv(BU.env_Cov, file = "MU_cov.csv", row.names = F)

