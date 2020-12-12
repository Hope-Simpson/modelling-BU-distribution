
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~ Identifying correlation between predictors ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# In this script, read in all the candidate predictors at a reduced spatial extent, 
# extract their values at the locations of BU occurrence points, 
# then check for collinearity between the predictors

rm(list=ls())
gc()

## 0. Set up the pathways to access to the source of data (spatial and others)
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_outputs <- paste(path_wd, "/Outputs/covariates", sep = "") 
path_pca <- paste(path_wd, "/PCA_BU", sep = "") 
path_scripts <- paste(path_wd, "/final_scripts", sep = "") 

setwd(path_wd)

## Call a function to install and load multiple R packages.
source(paste(path_scripts,"ipak.R", sep = "/"))
# List of packages
packages <- c("sp","raster","dismo","maptools","rgdal","proj4","ggplot2","gridExtra",
              "cowplot","biomod2","rJava","pROC","matrixStats","usdm","kernlab",
              "ks","sm","gbm","mgcv","nlme","Metrics","tidyr","readxl","RStoolbox","psych","vegan", "sf", "fasterize")
ipak(packages)

## 1. read the covariates matrix produced in script 2
setwd(path_outputs)

cova <- read.csv("BU_cov_PCR.csv")
names(cova)

res <- cor(cova, use = "complete.obs")
res <- round(res, 2)

getwd()
write.csv(res, "BU_correlation_matrix.csv", row.names = T)

# check for any variables >80% correlated with each other
# aridity index correlated with BC 16- retain BC 16 as precipitation has stronger evidence for association with BU
# max T correlated with min T- retain min T as this is expected to be more significant limiting factor

# AnPET_clip	bc14pr_dr_m	bc16pr_wt_q		minT	evi_clip			
# ED_osm_ww1	ED_dam_point1	EucDist_agr_area	EDdefor wetness
#  bc15_pr_ssn  strm1
 

###############################################################
###         NOW MOVE ALL OF THE SELECTED RASTERS            ###
###              AT CONTINENT-SCALE INTO THE                ###
###             human_predictors/cont   folder              ###
###############################################################

##### #########################################################
###            REPEAT FOR ENVIRONMENTAL DATA                ###
###############################################################
setwd(path_outputs)

## covariates extracted for 79 occurrence points
env_cova <- read.csv("MU_cov.csv")

names(env_cova)
res <- cor(env_cova, use = "complete.obs")
res <- round(res, 2)

setwd(path_outputs)
write.csv(res, "MU_correlation_matrix.csv", row.names = T)

# drop aridityind1 bc16pr_wt_q bc15_pr_ssn maxxT  strm1
# retain  bc14pr_dr_m    ED_dam_point1  ED_osm_ww1  EDdefor  EucDist_agr_area  evi_clip  minT

###############################################################
###         NOW MOVE ALL OF THE SELECTED RASTERS            ###
###              AT CONTINENT-SCALE INTO THE                ###
###                   env_cont   folder                     ###
###############################################################
