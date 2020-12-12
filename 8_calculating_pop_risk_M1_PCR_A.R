### Calulcating total population living in areas predicted suitable for BU, MU, and Both, with weighted mean and upper and lower bounds

gc()
rm(list = ls())

.rs.restartR()

# 1. Set up workspace directories to be used during this projects
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_data <- paste(path_wd, "/Tables", sep = "") 
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 
path_raster <- paste(path_wd, "/Raster", sep = "") 
path_model_PCR <- paste(path_wd, "/Outputs/BU_confirmed_cases", sep = "") 
path_model_env <- paste(path_wd, "/Outputs/MU_environmental", sep = "") 

setwd(path_wd)
# 2. Load the packages that we need for this project
source(paste(path_scripts,"ipak.R", sep = "/"))
# List of packages
packages <- c("tidyverse","reshape2","stats","readxl","foreign",
              "sf","dismo","raster","rgdal","proj4","ggplot2",
              "RStoolbox")
ipak(packages)

## 3 Load the African map with country boundaries
AFRO_map <- file.path(path_shapefiles,"salb.shp")
AFRO_map <- st_read(AFRO_map)
st_crs(AFRO_map)

# plot(AFRO_map$geometry)

## 3.1 Load the predicted environmental suitability for BU 
setwd(path_model_PCR)
raster.files <- list.files(path_model_PCR, 
                           pattern="*tif$", full.names=TRUE) 

raster.files
occ_model2_PCR <- raster(("BU_Occurrence.tif"))
LB_model2_PCR <- raster(("BU_Occurrence_LB.tif"))
UB_model2_PCR <- raster(("BU_Occurrence_UB.tif"))

PCS <- crs(occ_model2_PCR)

## 3.3 Load the predicted occurrence for env MU
setwd(path_model_env)
raster.files <- list.files(path_model_env, 
                           pattern="*tif$", full.names=TRUE) 

raster.files
occ_model2_env <- raster(("MU1_Occurrence.tif"))
LB_model2_env <- raster(("MU1_Occurrence_LB.tif"))
UB_model2_env <- raster(("MU1_Occurrence_UB.tif"))

occ_model2_env <- crop(extend(occ_model2_env, occ_model2_PCR),occ_model2_PCR)

occ_model2_env[occ_model2_env < 0] <- 0
# plot(UB_model2_env)

## 3.4 Load the continuous population density for 2020 (WorldPop project)
setwd(path_raster)
raster.files <- list.files(path_raster, 
                           pattern="*tif$", full.names=TRUE) 

PopDens2020 <- raster(("PopDens2020adj.tif"))
PopDens2020 <- projectRaster(PopDens2020, crs = PCS)
PopDens2020[PopDens2020<0] <- 0

#####################
## 4.2 Combine MU occurrence and BU occurrence surface
## Sum up both surface OCC
# first categorise 3 levels
# 2 = BU only, 1= MU only, 3 = both 
occ_model2_PCR[occ_model2_PCR == 1] <- 2
occ_model2_PCR[is.na(occ_model2_PCR)] <- 0
MU.BU.Occ <- occ_model2_PCR + occ_model2_env

setwd(path_model_PCR)
plot(MU.BU.Occ)

writeRaster(MU.BU.Occ, filename = "BU_MU_Occ2.tif", format= "GTiff", overwrite=TRUE)

## Reclassify to (0,1), so that 0. Absence or only 1 present;  2. = BU & MU
MU.BU.Occ[MU.BU.Occ < 3] <- 0
MU.BU.Occ[MU.BU.Occ == 3] <- 1
# plot(MU.BU.Occ)
# plot(AFRO_map$geometry, add = TRUE)

# do the same for LB and UB
MU.BU.LB <- LB_model2_PCR + LB_model2_env
MU.BU.LB[MU.BU.LB < 2] <- 0
MU.BU.LB[MU.BU.LB == 2] <- 1

MU.BU.UB <- UB_model2_PCR + UB_model2_env
MU.BU.UB[MU.BU.UB < 2] <- 0
MU.BU.UB[MU.BU.UB == 2] <- 1

setwd(path_model_PCR)
raster.files <- list.files(path_model_PCR, 
                           pattern="*tif$", full.names=TRUE) 

raster.files
occ_model2_PCR <- raster(("BU_Occurrence.tif"))

#  Extract Population estimates for each category

## Population density for BU
## Align both raster dataset: crop, extend and resample processing

PopDens2020 <- crop(extend(PopDens2020, occ_model2_PCR), occ_model2_PCR)
PopDens2020 <- raster::resample(PopDens2020, occ_model2_PCR, method = "ngb")

# generate raster showing population in each occurence cell
BU.PopDens2020r <- PopDens2020 * occ_model2_PCR
MU.PopDens2020r <- PopDens2020 * occ_model2_env
MU.BU.PopDens2020r <- PopDens2020 * MU.BU.Occ

# 7. Estimate total population by country across Africa
## 7.1 Total Country Population
AFRO_map <- st_transform(AFRO_map, crs = 3857)
AFRO.tb <- AFRO_map %>%
  st_set_geometry(NULL)

Total.Pop2020 <- raster::extract(PopDens2020, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
Total.Pop2020 <- as.data.frame(Total.Pop2020)

names(AFRO.tb)
Total.Pop2020 <- cbind(AFRO.tb[,c(3,7,10)],Total.Pop2020)
colnames(Total.Pop2020)[5] <- c("TotalPop2020") 

## 7.2 Extract population in areas of BU presence at national level
BU.Pop2020e <- raster::extract(BU.PopDens2020r, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.Pop2020 <- as.data.frame(BU.Pop2020e)
BU.Pop2020 <- cbind(AFRO.tb[,c(3,7,10)],BU.Pop2020)
colnames(BU.Pop2020)[5] <- c("BUPop2020") 

## 7.3  Extract population in areas of MU presence 
MU.Pop2020e <- raster::extract(MU.PopDens2020r, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.Pop2020 <- as.data.frame(MU.Pop2020e)
MU.Pop2020 <- cbind(AFRO.tb[,c(3,7,10)],MU.Pop2020)

colnames(MU.Pop2020)[5] <- c("MU.Pop2020") 
head(BU.Pop2020)

## 7.3 Extract population in areas of MU + BU presence
MU.BU.Pop2020e <- raster::extract(MU.BU.PopDens2020r, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.Pop2020 <- as.data.frame(MU.BU.Pop2020e)
MU.BU.Pop2020 <- cbind(AFRO.tb[,c(3,7,10)],MU.BU.Pop2020)

colnames(MU.BU.Pop2020)[5] <- c("MU.BU.Pop2020") 
head(MU.BU.Pop2020)

## 7.5 Assemble datasets and export as a CSV table
Assemble.Pop2020 <- cbind(Total.Pop2020, BU.Pop2020[5], MU.Pop2020[5], MU.BU.Pop2020[5] )

## 7.6 Finally, estimate total population affected  and credible intervals
BU.Pop2020.LBr <- PopDens2020 * LB_model2_PCR
BU.Pop2020.UBr <- PopDens2020 * UB_model2_PCR
BU.Pop2020.WMr <- PopDens2020 * occ_model2_PCR

MU.Pop2020.LBr <- PopDens2020 * LB_model2_env
MU.Pop2020.UBr <- PopDens2020 * UB_model2_env
MU.Pop2020.WMr <- PopDens2020 * occ_model2_env

MU.BU.Pop2020.LBr <- PopDens2020 * MU.BU.LB
MU.BU.Pop2020.UBr <- PopDens2020 * MU.BU.UB
MU.BU.Pop2020.WMr <- PopDens2020 * MU.BU.Occ

## Extract lower bound population living in at-risk areas by country
BU.Pop2020.LBe <- raster::extract(BU.Pop2020.LBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.Pop2020.LB <- as.data.frame(BU.Pop2020.LBe)

# lower bound area
BU.area.LBe <- raster::extract(LB_model2_PCR, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.area.LB <- as.data.frame(BU.area.LBe)

head(BU.Pop2020.LB)
head(BU.Pop2020.LBe)
head(BU.Pop2020.LB)

BU.Pop2020.LB <- cbind(AFRO.tb[,c(3,7,10)], BU.Pop2020.LB, BU.area.LB[2])
colnames(BU.Pop2020.LB)[5] <- c("BUPop2020_LB") 
colnames(BU.Pop2020.LB)[c(6)] <- c("BUarea_LB") 

BU.Pop2020.LB$BUarea_LB <- BU.Pop2020.LB$BUarea_LB * 25

########
# MU lower bound  pop
MU.Pop2020.LBe <- raster::extract(MU.Pop2020.LBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.Pop2020.LB <- as.data.frame(MU.Pop2020.LBe)

# MU lower bound  pop area
MU.area.LBe <- raster::extract(LB_model2_env, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.area.LB <- as.data.frame(MU.area.LBe)

MU.Pop2020.LB <- cbind(AFRO.tb[,c(3,7,10)], MU.Pop2020.LB, MU.area.LB[2])
colnames(MU.Pop2020.LB)[5] <- c("MUPop2020_LB") 
colnames(MU.Pop2020.LB)[6] <- c("MUarea_LB") 

MU.Pop2020.LB$MUarea_LB <- MU.Pop2020.LB$MUarea_LB * 25

#### combined
## Lower Bound population
MU.BU.Pop2020.LBe <- raster::extract(MU.BU.Pop2020.LBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.Pop2020.LB <- as.data.frame(MU.BU.Pop2020.LBe)

# lower bound area
MU.BU.area.LBe <- raster::extract(MU.BU.LB, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.area.LB <- as.data.frame(MU.BU.area.LBe)

MU.BU.Pop2020.LB <- cbind(AFRO.tb[,c(3,7,10)], MU.BU.Pop2020.LB, MU.BU.area.LB[2])
colnames(MU.BU.Pop2020.LB)[5] <- c("MU.BUPop2020_LB") 
colnames(MU.BU.Pop2020.LB)[c(6)] <- c("MU.BUarea_LB") 

MU.BU.Pop2020.LB$MU.BUarea_LB <- MU.BU.Pop2020.LB$MU.BUarea_LB * 25

## Upper Bound
BU.Pop2020.UBe <- raster::extract(BU.Pop2020.UBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.Pop2020.UB <- as.data.frame(BU.Pop2020.UBe)

BU.area.UBe <- raster::extract(UB_model2_PCR, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.area.UB <- as.data.frame(BU.area.UBe)

BU.Pop2020.UB <- cbind(AFRO.tb[,c(3,7,10)],BU.Pop2020.UB, BU.area.UB[2])
colnames(BU.Pop2020.UB)[5] <- c("BUPop2020_UB") 
colnames(BU.Pop2020.UB)[6] <- c("BUarea_UB") 

head(BU.Pop2020)

# multiply each cell by 25 to get area in sqkm
BU.Pop2020.UB$BUarea_UB <- BU.Pop2020.UB$BUarea_UB * 25

###  Upper bound population 
MU.Pop2020.UBe <- raster::extract(MU.Pop2020.UBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.Pop2020.UB <- as.data.frame(MU.Pop2020.UBe)

MU.area.UBe <- raster::extract(UB_model2_env, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.area.UB <- as.data.frame(MU.area.UBe)

MU.Pop2020.UB <- cbind(AFRO.tb[,c(3,7,10)],MU.Pop2020.UB, MU.area.UB[2])
colnames(MU.Pop2020.UB)[5] <- c("MUPop2020_UB") 
colnames(MU.Pop2020.UB)[6] <- c("MUarea_UB") 

MU.Pop2020.UB$MUarea_UB <- MU.Pop2020.UB$MUarea_UB * 25

### mu bu 
MU.BU.Pop2020.UBe <- raster::extract(MU.BU.Pop2020.UBr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.Pop2020.UB <- as.data.frame(MU.BU.Pop2020.UBe)

MU.BU.area.UBe <- raster::extract(MU.BU.UB, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.area.UB <- as.data.frame(MU.BU.area.UBe)

MU.BU.Pop2020.UB <- cbind(AFRO.tb[,c(3,7,10)],MU.BU.Pop2020.UB, MU.BU.area.UB[2])
colnames(MU.BU.Pop2020.UB)[5] <- c("MU.BUPop2020_UB") 
colnames(MU.BU.Pop2020.UB)[6] <- c("MU.BUarea_UB") 

# multiply each cell by 25 to get area in sqkm
MU.BU.Pop2020.UB$MU.BUarea_UB <- MU.BU.Pop2020.UB$MU.BUarea_UB  * 25

## Weighted Mean
BU.Pop2020.WMe <- raster::extract(BU.Pop2020.WMr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.Pop2020.WM <- as.data.frame(BU.Pop2020.WMe)

BU.area.WMe <- raster::extract(occ_model2_PCR, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
BU.area.WM <- as.data.frame(BU.area.WMe)
head(BU.area.WM)
BU.Pop2020.WM <- cbind(AFRO.tb[,c(3,7,10)],BU.Pop2020.WM, BU.area.WM[2])

colnames(BU.Pop2020.WM)[5] <- c("BUPop2020WM") 
colnames(BU.Pop2020.WM)[6] <- c("BUarea_WM") 

BU.Pop2020.WM$BUarea_WM <- BU.Pop2020.WM$BUarea_WM * 25

#### MU 
MU.Pop2020.WMe <- raster::extract(MU.Pop2020.WMr, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.Pop2020.WM <- as.data.frame(MU.Pop2020.WMe)

MU.area.WMe <- raster::extract(occ_model2_env, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.area.WM <- as.data.frame(MU.area.WMe)

MU.Pop2020.WM <- cbind(AFRO.tb[,c(3,7,10)],MU.Pop2020.WM, MU.area.WM[2])
colnames(MU.Pop2020.WM)[5] <- c("MUPop2020") 
colnames(MU.Pop2020.WM)[6] <- c("MUarea_WM") 

MU.Pop2020.WM$MUarea_WM <- MU.Pop2020.WM$MUarea_WM * 25

## combined
MU.BU.Pop2020.WMe <- raster::extract(MU.BU.Pop2020.WMr, AFRO_map, fun=sum,  na.rm=TRUE, df=TRUE)
MU.BU.Pop2020.WM <- as.data.frame(MU.BU.Pop2020.WMe)

MU.BU.area.WMe <- raster::extract(MU.BU.Occ, AFRO_map, fun=sum, na.rm=TRUE, df=TRUE)
MU.BU.area.WM <- as.data.frame(MU.BU.area.WMe)

MU.BU.Pop2020.WM <- cbind(AFRO.tb[,c(3,7,10)],MU.BU.Pop2020.WM, MU.BU.area.WM[2])
colnames(MU.BU.Pop2020.WM)[5] <- c("MU.BU.Pop2020") 
colnames(MU.BU.Pop2020.WM)[6] <- c("MU.BU.area_WM") 

MU.BU.Pop2020.WM$MU.BU.area_WM <- MU.BU.Pop2020.WM$MU.BU.area_WM * 25

head(Assemble.MU.Pop2020)

## Assemble tables
Assemble.BU.Pop2020 <- cbind(Total.Pop2020, BU.Pop2020.WM[c(5,6)], BU.Pop2020.LB[c(5,6)], BU.Pop2020.UB[c(5,6)])
Assemble.MU.Pop2020 <- cbind(Total.Pop2020, MU.Pop2020.WM[c(5,6)], MU.Pop2020.LB[c(5,6)], MU.Pop2020.UB[c(5,6)])
Assemble.MU.BU.Pop2020 <- cbind(Total.Pop2020, MU.BU.Pop2020.WM[c(5,6)], MU.BU.Pop2020.LB[c(5,6)], MU.BU.Pop2020.UB[c(5,6)])
head(Assemble.BU.Pop2020)

# 8. Export CSV tables with final results
setwd(path_data)

write.csv(Assemble.Pop2020, file = "Assemble_Pop2020_PCRA_oct.csv", row.names = FALSE, na="")
write.csv(Assemble.BU.Pop2020, file = "Assemble_BU_Pop2020_PCRA_oct.csv", row.names = FALSE, na="")
write.csv(Assemble.MU.Pop2020, file = "Assemble.MU.Pop2020_PCRA_oct.csv", row.names = FALSE, na="")
write.csv(Assemble.MU.BU.Pop2020, file = "Assemble.MU.BU.Pop2020_PCRA_oct.csv", row.names = FALSE, na="")
getwd()

BU_area <- freq(occ_model2_PCR, useNA='no')
apc <- prod(res(occ_model2_PCR))
BU_area <- cbind(BU_area, area=BU_area[,2] * apc)

BU_areaLB <- freq(LB_model2_PCR, useNA='no')
apc <- prod(res(LB_model2_PCR))
BU_areaLB <- cbind(BU_areaLB, area=BU_areaLB[,2] * apc)

BU_areaUB <- freq(UB_model2_PCR, useNA='no')
apc <- prod(res(UB_model2_PCR))
BU_areaUB <- cbind(BU_areaUB, area=BU_areaUB[,2] * apc)

#####################

MU_area <- freq(occ_model2_env, useNA='no')
apc <- prod(res(occ_model2_env))
MU_area <- cbind(MU_area, area=MU_area[,2] * apc)

BU_areaLB <- freq(LB_model2_PCR, useNA='no')
apc <- prod(res(LB_model2_PCR))
BU_areaLB <- cbind(BU_areaLB, area=BU_areaLB[,2] * apc)

BU_areaUB <- freq(UB_model2_PCR, useNA='no')
apc <- prod(res(UB_model2_PCR))
BU_areaUB <- cbind(BU_areaUB, area=BU_areaUB[,2] * apc)

################################

# Clean up workspace of unnecessary objects and save objects created in the corresponding folder
gc()
save.image(paste0(path_model_PCR, sep = "/", "Population_Estimates.RData"))
