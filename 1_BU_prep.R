# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~ BU distribution modelling across Africa ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ 1. Preparing the data                       ~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

gc()

rm(list=ls())

## 0. Set up the paths to access data 
# set your own file path
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_shapefiles <- paste(path_wd, "/Shapefiles", sep = "") 
path_data <- paste(path_wd, "/Tables", sep = "") 
path_outputs <- paste(path_wd, "/Outputs", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

setwd(path_wd)

## 1. Call a function to install and load multiple R packages.
source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("sp","raster","dismo","maptools","rgdal","proj4","ggplot2","gridExtra",
              "cowplot","biomod2","rJava","pROC","matrixStats","usdm","kernlab",
              "ks","sm","gbm","mgcv","nlme","Metrics","tidyr","readxl","RStoolbox","psych","vegan", "sf")
ipak(packages)

## 2a. Import list of occurrences for BU Human cases
setwd(path_data)

# from literature records (open access)
BU_occurrences <- read.csv("human_occurrences-LIT_2020_23_10.csv", na.strings = "", stringsAsFactors = F)
nrow(BU_occurrences) # 3126

# from national programme records (restricted access)
# BU_occurrences_NPR <- read.csv("human_occurrences-NPR_2020_23_10.csv", na.strings = "", stringsAsFactors = F)
# BU_occurrences_NPR[is.na(BU_occurrences_NPR$max_year),]$max_year <- 2004

# combine into a single df
# BU_occurrences <- rbind(BU_occurrences, BU_occurrences_NPR)

# remove absence records (most are not survey absences, so not reliable)
BU_occurrences <- BU_occurrences[(BU_occurrences$number>0),] #   5490 
nrow(BU_occurrences)

# remove non-georeferenced locations
BU_occurrences <- BU_occurrences[!is.na(BU_occurrences$lat) & !is.na(BU_occurrences$long),] # 3728
nrow(BU_occurrences)

# remove non-reliable georeferenced locations
BU_occurrences <- BU_occurrences[(BU_occurrences$Geo_reliability<5) & !is.na(BU_occurrences$Geo_reliability),] #  3531
nrow(BU_occurrences)

table(BU_occurrences$max_year, useNA = "ifany")

# 3a. Assign contemporariness score
BU_occurrences$year_score <- 1
BU_occurrences[BU_occurrences$max_year<2003,]$year_score <- 0.5
BU_occurrences[BU_occurrences$max_year<1991,]$year_score <- 0.25
# table(BU_occurrences$max_year, BU_occurrences$year_score, useNA = 'ifany')

# 4a. Assign weight to records- mean of year and diagnostic scores
BU_occurrences$weight <- (BU_occurrences$year_score + BU_occurrences$diagnostic_score)/2

# order by weight so that 
# if any points are dropped by deduplicaiton, the point with higher weight will be retained
BU_occurrences <-  BU_occurrences[with(BU_occurrences, order(-weight)), ]

## 5a. Create a spatial layer from the occurrences
BU_occurrences.sp <- SpatialPoints(BU_occurrences[,c("long","lat")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
BU_occurrences.ly <- SpatialPointsDataFrame(BU_occurrences.sp, BU_occurrences)

## 2b. Import list of occurrences for BU environmental detection and animal cases
env_occurrences <- read.csv("environmental_results-LIT_2020_23_09.csv", stringsAsFactors = FALSE) #  406
env_occurrences$lat <- as.numeric(as.character(env_occurrences$lat))
env_occurrences$long <- as.numeric(as.character(env_occurrences$long))
nrow(env_occurrences)

# drop non-confirmed occurrences
env_occurrences <- env_occurrences[(env_occurrences$N_positive>0),] # 326 
nrow(env_occurrences)

env_occurrences <- env_occurrences[(env_occurrences$MU_confirmed == "yes"), ]#  173
nrow(env_occurrences)

# drop non-georeferenced occurrences
env_occurrences <- env_occurrences[!is.na(env_occurrences$lat),] #  151

# drop non-reliably-georeferenced occurrences
env_occurrences <- env_occurrences[(env_occurrences$Geo_reliability<5) & !is.na(env_occurrences$Geo_reliability),] #   133
nrow(env_occurrences)

## 3b. Create a spatial layer from the environmental occurrences
env_occurrences.sp <- SpatialPoints(env_occurrences[,c("long","lat")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
env_occurrences.ly <- SpatialPointsDataFrame(env_occurrences.sp, env_occurrences)

## 6. Import map of African countries
setwd(path_shapefiles)

Africa_ADM0 <- readOGR(dsn=path_shapefiles, layer="salb")
PCS <- Africa_ADM0@proj4string

setwd(path_wd)

## 7. Re-project point data to the projected coordinate system of Africa layer
BU_occurrences.ly <- spTransform(BU_occurrences.ly, PCS)
env_occurrences.ly <- spTransform(env_occurrences.ly, PCS)

# plot(Africa_ADM0)
# plot(BU_occurrences.ly, add = T)

 ## 8. Remove occurrence records that do not fall within mainland
BU_occurrences.ly <- BU_occurrences.ly[Africa_ADM0,] # remove occurrence record  out of Africa mainland
nrow(BU_occurrences.ly)  # 3316

env_occurrences.ly <- env_occurrences.ly[Africa_ADM0,] # remove occurrence record  out of Africa mainland
nrow(env_occurrences.ly)  # 3316

# 9. generate clinical dataset
BU_clinical.ly <- BU_occurrences.ly  
BU_clinical.ly <- BU_clinical.ly[BU_clinical.ly$diagnostic_score>0.25,]
BU_clinical.sp <- as(BU_clinical.ly, "SpatialPoints")
nrow(BU_clinical.ly@data)  # 3316

# 10. generate lab-confirmed dataset
## restrict to PCR and histologically confirmed cases
BU_conf_occurrences.ly <- BU_occurrences.ly[(BU_occurrences.ly$diagnostic_score>0.99),] 
BU_conf_occurrences.sp <- as(BU_conf_occurrences.ly, "SpatialPoints")
nrow(BU_conf_occurrences.ly)  # 1336

## 11. Identify and remove duplicate locations
dups <- duplicated(BU_conf_occurrences.ly@data[,c("long","lat")])
length(which(dups== TRUE)) # 598 duplicated locations (based on coordinates); 742 remain
BU_conf_occurrences.ly <- BU_conf_occurrences.ly[!dups,]
nrow(BU_conf_occurrences.ly)   # 738

rm(dups)

## 11. Identify  and remove duplicate locations from BU_clinical.ly
dups <- duplicated(BU_clinical.ly@data[,c("long","lat")])
length(which(dups== TRUE)) # 1135 duplicated locations (based on coordinates); 2181 remain

BU_clinical.ly <- BU_clinical.ly[!dups,]

nrow(BU_clinical.ly)
table(BU_conf_occurrences.ly$ISO_code3)
table(BU_clinical.ly$ISO_code3)

rm(dups)

env_occurrences.sp <- as(env_occurrences.ly, "SpatialPoints")
nrow(env_occurrences.ly)

## 11. Identify  and remove duplicate locations from env_occurrences.ly
dups <- duplicated(env_occurrences.ly@data[,c("long","lat")])
length(which(dups== TRUE)) # 11 duplicated locations (based on coordinates) 79 remain

env_occurrences.ly <- env_occurrences.ly[!dups,]
dim(env_occurrences.ly) #   unique occurrences

rm(dups)
names(BU_clinical.ly)

table(BU_clinical.ly$max_year)

BU_clinical.ly$max_year <- as.numeric(as.character(BU_clinical.ly$max_year))

BU_conf_occurrences.ly$max_year <- as.numeric(as.character(BU_conf_occurrences.ly$max_year))

## 12. Save individual datasets
writeOGR(BU_clinical.ly, dsn = (paste(path_shapefiles,"/BU_occurrences", sep = "")), layer ="BU_occurrences_all", 
         driver = "ESRI Shapefile", overwrite_layer = T)

nrow(BU_conf_occurrences.ly)
writeOGR(BU_conf_occurrences.ly, dsn = (paste(path_shapefiles,"/BU_occurrences", sep = "")), layer ="BU_occurrences_PCR", 
         driver = "ESRI Shapefile", overwrite_layer = T)

plot(BU_occurrences.ly)

writeOGR(env_occurrences.ly, dsn = (paste(path_shapefiles,"/BU_occurrences", sep = "")), layer ="env_occurrences", 
         driver = "ESRI Shapefile", overwrite_layer = T)


env.table <- env_occurrences.ly@data
confirmed_BU.table <- BU_occurrences.ly@data
all_occ.table <- BU_clinical.ly@data

nrow(BU_occurrences.ly)
nrow(BU_clinical.ly)
nrow(BU_conf_occurrences.ly)

nrow(env_occurrences.ly)

setwd(path_data)

write.csv(all_occ.table, "all_occ.table.csv", row.names = F)

