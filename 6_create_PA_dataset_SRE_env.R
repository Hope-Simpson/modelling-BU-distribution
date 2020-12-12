# https://www.nature.com/articles/nature12060#Sec8

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~ Generate pseudoabsence data for MU models ~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# In this script, select PA points based on evidence consensus score at adm1 and 0 level, 
# oversampling in areas with lower evidence for MU
# outside of area predicted suitable by SRE

rm(list=ls())
gc()

path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_data <- paste(path_wd, "/Tables", sep = "") 
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_BG <- paste(path_wd, "/Shapefiles/BG_points", sep = "") 
path_raster <- paste(path_wd, "/Raster", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("sf", "rgdal", "scales", "tidyverse", "raster", "fasterize","biomod2", "spatialEco" )
ipak(packages)

## Import an  map of Africa
setwd(path_shapefiles)
Africa_ADM0 <- st_read("salb.shp")

## read environmental occurrences dataset
occurrences_env <- st_read("BU_occurrences/env_occurrences.shp")
occurrences_env <- st_transform(occurrences_env, crs = st_crs(Africa_ADM0))

occurrences_env <- occurrences_env[,c("N_postv", "MU_cnfr")]
names(occurrences_env)

occurrences_env <- occurrences_env %>% rename(MU = N_postv,
                              weight = MU_cnfr)

## 
# plot(st_geometry(Africa_ADM0))
# plot(st_geometry(occurrences_env), add = T)

## read adm1 shapefile from GADM
setwd(path_shapefiles)
EC_vector <- st_read("adm1_africa.shp")

names(EC_vector)

# read in results from BU systematic review
setwd(path_data)

adm0evidence <- read.csv("adm0_scores2.csv", stringsAsFactors = F)
adm1evidence <- read.csv("adm1_scores2.csv", stringsAsFactors = F)
adm1_env <- read.csv("adm1_scores_env.csv", stringsAsFactors = F)

names(adm1_env)
adm1_env <- adm1_env[c( "ADM1_CODE", "pct")]
adm1_combined <- left_join(adm1evidence, adm1_env, by = c("gadm1"= "ADM1_CODE" ) )
head(adm1_combined)
adm1_combined[is.na(adm1_combined$pct),]$pct <- 0
# rescale adm1 scores to 0-1
adm1_combined$pct <- scales::rescale(adm1_combined$pct, c(0, 1))

# assign combined score as the max of adm1 BU and MU evidence
adm1_combined$adm1_score_combined <- max(adm1_combined$adm1_score, adm1_combined$pct)
adm1_combined <- adm1_combined[c( "gadm1", "adm1_score_combined")]

EC_vector <- left_join(EC_vector, adm0evidence, by = c("GID_0"= "ISO3" ))
EC_vector <- left_join(EC_vector, adm1evidence, by = c("GID_1"= "gadm1" ))

# rescale adm0 score to 0- 1
EC_vector$adm0_score <- scales::rescale(EC_vector$adm0_score, c(0, 1))

# scale by 0.5
EC_vector$adm0_score <- (EC_vector$adm0_score * 0.5)
nrow(EC_vector[is.na(EC_vector$adm1_score) & !is.na(EC_vector$adm0_score ),])


# If there is  evidence at ADM1 level, assign the ADM1 evidence score
# if no  evidence at adm1 level, combined score is adm0 score scaled by a factor of 0.5
EC_vector[is.na(EC_vector$adm1_score) & !is.na(EC_vector$adm0_score ),]$adm1_score <- EC_vector[is.na(EC_vector$adm1_score) & !is.na(EC_vector$adm0_score ),]$adm0_score

head(EC_vector)

### only keep required fields
EC_vector <- EC_vector[,c(1,2,6)]

names(EC_vector) <- c("GID_0", "GID_1", "cmb_scr", "geometry")

setwd(path_shapefiles)

st_write(EC_vector, "ev_adm1_0_africa_env.shp", delete_layer =  TRUE)
EC_vector <- st_read("ev_adm1_0_africa_env.shp")
# plot(EC_vector)

## rasterise EC 
raster.files <- list.files(paste(path_raster, "/human_predictors/cont", sep = ""), 
                           pattern="*tif$", full.names=TRUE) 

# create a mask
AI <- raster((raster.files[1]))
PCSR <- AI@crs

# creat a vector mask of the continent
Mask <- AI
values(Mask) <- 1
Mask <- mask(Mask, AI)
# Mask <- aggregate(Mask, 5)

# rasterise EC
EC_raster <- fasterize(EC_vector, Mask, field = "cmb_scr", fun = "last", background = NA_real_,
                       by = NULL)

# plot(EC_raster)
setwd(path_raster)

writeRaster(EC_raster, "EC_raster_env.tif", overwrite = T)
EC_raster <- raster("EC_raster_env.tif")

# plot(EC_raster)

dev.off()

# create a buffer 10km around occurrences_env to prevent contamination of PA points with true occurrences 
all_10kbuffer <- st_buffer(occurrences_env, 10000)
all_10kbuffer <- st_geometry(all_10kbuffer)
all_10kbuffer <- st_union(all_10kbuffer)
all_10kbuffer <- st_sf(all_10kbuffer)

setwd(path_shapefiles)
st_write(all_10kbuffer, "buffer10k_env.shp",  delete_layer = TRUE)

all_10kbuffer <- st_read("buffer10k_env.shp")
all_10kbuffer <- st_transform(all_10kbuffer, st_crs(AI))

plot(all_10kbuffer)

## rasterise buffer areas
buffer_raster <- fasterize(all_10kbuffer, Mask, field = NULL, fun = "last", background = NA_real_,
                           by = NULL)

# read in raster that prevents selection close to occurrence locations
setwd(path_raster)
writeRaster(buffer_raster, "10kbuffer_raster_env.tif", overwrite = T)
buffer_raster_10k <- raster("10kbuffer_raster_env.tif")

# invert the raster so the buffers get values of 0 and area outside buffer gets value of 1
buffer_raster_10k_invert <- buffer_raster_10k
buffer_raster_10k_invert[buffer_raster_10k_invert$X10kbuffer_raster_env ==1] <- 0
buffer_raster_10k_invert[is.na(buffer_raster_10k_invert$X10kbuffer_raster_env )] <- 1

# plot(buffer_raster_10k_invert)

setwd(path_wd)
cov.5km <- readRDS("covariates_5km_env.rds")
names(cov.5km)
PCSR <- Mask@crs

##  Generate a raster stack object with the predictors
myExpl <- stack()
for(i in 3:length(cov.5km)){
  raster1 <- rasterFromXYZ(cov.5km[,c(1,2,i)],
                           crs = PCSR)
  myExpl <- stack(myExpl, raster1)
}

names(myExpl)
myRespXY <- as.data.frame(st_coordinates(occurrences_env)) 
myRespXY$occurrence <- 1
myResp <- myRespXY[,3]

# build a raster layer of the response variable (MU)
myResp <- reclassify(subset(myExpl,1,drop=TRUE), c(-Inf,Inf,0))
myResp[cellFromXY(myResp, myRespXY)] <- 1
# plot(myResp)

# Compute SRE 
SRE <- sre(Response = myResp, Explanatory = myExpl, NewData=myExpl, Quant=0.025)

setwd(path_raster)

writeRaster(SRE, "SRE_env.tif", format= "GTiff", overwrite=TRUE)

# visualise results
plot(SRE, main="Perc025")

# adjust SRE raster, the 10k buffer raster and the EC raster to same extent as Mask
SRE_invert <- raster.invert(SRE)
SRE_invert <- extend(SRE_invert, Mask)
SRE_invert <- raster::resample(SRE_invert, Mask, method = 'bilinear')
SRE_invert <-  crop(SRE_invert, Mask)

buffer_raster_10k_invert <- extend(buffer_raster_10k_invert, Mask)
buffer_raster_10k_invert <- raster::resample(buffer_raster_10k_invert, Mask, method = 'bilinear')
buffer_raster_10k_invert <-  crop(buffer_raster_10k_invert, Mask)

EC_raster <- extend(EC_raster, Mask)
EC_raster <- raster::resample(EC_raster, Mask, method = 'bilinear')
EC_raster <-  crop(EC_raster, Mask)

#  Exclude areas outside of Africa (Mask=0), within the SRE envelope (SRE_invert = 0) or within 10km of occurrence points  (buffer_raster_10k_invert = 0)
## PA selection frame is area within Africa, not within 25% CI for BU suitability, not within 10km of a case
# evidence consensus score will be used to determine probability of selection
summary(EC_raster$EC_raster_env)
# EC_raster$EC_raster_env <- round(EC_raster$EC_raster_env, 3)
PA_limit <- Mask * SRE_invert * buffer_raster_10k_invert * EC_raster 

plot(PA_limit)

setwd(path_raster)
writeRaster(PA_limit, "PA_limit_SRE_env.tif", format= "GTiff", overwrite=TRUE)
PA_limit <- raster("PA_limit_SRE_env.tif")

# 3.8. Potential PA points will be sampled from a list of point data
# extract points from the raster 
Pseudo_points <- rasterToPoints(PA_limit, spatial = F)
Pseudo_points <- as.data.frame(Pseudo_points)

# drop points with a PA score of 0 
Pseudo_points <- Pseudo_points[Pseudo_points$PA_limit_SRE_env>0,]
nrow(Pseudo_points)

#generate field of random values from 0-1
set.seed(123)
Pseudo_points$s <- runif(nrow(Pseudo_points))
head(Pseudo_points)

## PA points are those with EC layer score less than random score (Bhatt Nature 2013)
P_absence <- Pseudo_points[which(Pseudo_points$PA_limit_SRE_env<Pseudo_points$s),]
nrow(P_absence)

# 3.9. Now  sample for PA points - equal number to occurrence points
set.seed(123)
pseudoabs <- P_absence[sample(seq(1:nrow(P_absence)), size = nrow(occurrences_env), replace = T),1:4]
names(pseudoabs) <- c("x", "y","weight", "MU") 
pseudoabs <- pseudoabs[,c(1,2,4,3)]
pseudoabs$MU <- 0

pseudoabs$weight <- NA
head(pseudoabs)
PCSR <- st_crs(EC_raster)

# convert to sf
PA_sf_alloc = st_as_sf(pseudoabs, coords = c("x", "y"), crs = st_crs(occurrences_env))
plot(st_geometry(PA_sf_alloc), pch = 20, cex= 0.3, col= "red", add = T)

setwd(path_BG)
st_write(PA_sf_alloc, "pseudoabs_env.shp", delete_layer  = TRUE)
