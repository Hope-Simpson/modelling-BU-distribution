#  https://www.rdocumentation.org/packages/biomod2/versions/3.3-7.1/topics/BIOMOD_Modeling

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~ Modelling Spatial Distribution of Buruli ulcer across Africa ~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## Start Date: 10.09.2019
## End Date: 10.09.201
## Place: London, UK
## Project: Delineating environmental limits for Buruli ulcer across Africa
## Modelling approach: Environmental modelling using occurrence of Buruli ulcer and 
## including pseudoabsence points from countries with evidence consensus for absence
## Background points created randomly around presences to account for potential
## geographical bias


## 0. Set up file paths 
gc()

rm(list=ls())

# 1. Set up workspace directories to be used during this projects
path_shapefiles <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/Shapefiles"
path_data <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/Tables"
path_covars <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_BG <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/Shapefiles/BG_points"
path_maxent <- "C:/Users/Hope/Desktop/max_ent"
path_wd <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
path_func <- "C:/Users/Hope/Dropbox (LSoHaTM)/ESPEN_SCH/R scripts"
path_outputs <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/Outputs/BU_confirmed_cases"
path_scripts <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/final_scripts"

## 1. Call a function to install and load multiple R packages.
source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("tidyverse","reshape2","stats","readxl","foreign",
              "sf","dismo","raster","rgdal","proj4","ggplot2","cowplot",
              "rJava","pROC","pryr","RStoolbox","vegan", "biomod2",
              "psycho","corrplot","spatialEco","maptools","gbm","spatstat", "exactextractr", "sdmvspecies")
ipak(packages)

# 2. Importing data to be used in the models
## 2.1 Load the Africa map with country boundaries
AFRO_map <-  readOGR(dsn=path_shapefiles, layer="salb")

## 2.2 Load the data points: occurrence of BU 
setwd(path_wd)
AFRO_BU <- file.path(path_shapefiles, "/BU_occurrences/BU_occurrences_PCR.shp")
AFRO_BU <- st_read(AFRO_BU)

AFRO_BU <- AFRO_BU[,c("number", "weight")]
AFRO_BU <- AFRO_BU %>% rename(BU = number)
AFRO_BU <- st_transform(AFRO_BU, crs = 3857  )

AFRO_BU$BU <- 1
occurrences <- AFRO_BU
nrow(occurrences)

plot(AFRO_map)
plot(occurrences$geometry, cex = .2, col = "red", add = TRUE)

## 2.3 Load the pseudo-absence points generated for all occurrence locations (6_create_PA_dataset_SRE) 
pseudo_abs <- file.path(path_BG, "pseudoabs_alloc.shp")
pseudo_abs <- st_read(pseudo_abs)
pseudo_abs <- st_transform(pseudo_abs, crs = 3857)
pseudo_abs$type <- "PA"
pseudo_abs$BU <- 0
pseudo_abs <- pseudo_abs[,c(1,2,4,3)]
nrow(pseudo_abs) ## 2181 PA points created (same number as all occurrences)
names(pseudo_abs)

### take random sample of PA points to get number equal to confirmed cases 
set.seed(123)
pseudo_abs <- pseudo_abs[sample(seq(1:nrow(pseudo_abs)), size = nrow(occurrences), replace = F),]
nrow(pseudo_abs) ## 738 PA points created (same number as confirmed occurrences)

plot(pseudo_abs$geometry, cex = .2, col = "blue", add = TRUE)
setwd(path_BG)
st_write(pseudo_abs, "pseudo_abs.SRE.PCR.shp", delete_layer = T)

## 2.4 Load the background points   created in 5_KDE_surface_PCR.R
BG <- file.path(path_BG, "KDE_BG_PCR.shp")    
BG <- st_read(BG)
BG <- st_transform(BG, crs = st_crs(AFRO_map))

nrow(BG) ## 
BG$weight <- 0 # assign weight of 0 to background data
BG$type <- "BG"
plot(BG$geometry, cex = .2, col = "grey", add = TRUE)

# BG <- BG[,c(1,3,4,2)]

names(BG)

# 3. Assemble background and pseudo-absences
PA.BKG.data <- rbind(pseudo_abs, BG)
table(PA.BKG.data$type) # 1380 records
head(PA.BKG.data)

PA.BKG.data <- st_transform(PA.BKG.data, crs = st_crs(AFRO_map))

# 4. Load the covariates to be used in the modelling
# saved as a matrix in 4_covariates matrix
setwd(path_wd)

cov.5km <- readRDS("covariates_5km_alloc.rds")
PCS = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"

# loop through columns in the matrix and generate a raster stack
covariates <- stack()
for(i in 3:length(cov.5km)){
  raster1 <- rasterFromXYZ(cov.5km[,c(1,2,i)],
                           crs = PCS)
  covariates <- stack(covariates, raster1)
}

names(covariates)
# warnings()
rm(i, raster1)

## remove wetness index as it makes little contribution (after running initial modelling step)
covariates <- covariates[[-c(12)]]

## convert distance to dams + agr areas from m to km
covariates[[5]] <- covariates[[5]] /1000
covariates[[7]] <- covariates[[7]] /1000

# convert temp in F to C
# °C = (°F - 32)/1.8
names(covariates)

dev.off()

# 5. Prepare occurrence and background points datasets for being use in the model
## 5.1 Extract predictor values for occurrence locations (check for potential NA)
BU.pred <- raster::extract(covariates, occurrences)
BU.pred <- as.data.frame(BU.pred)

## 5.2 Cross-check for  NA values that can make the model crash.
nl <- which(complete.cases(BU.pred)) # ## 1041 observations complete data
occurrences <- occurrences[nl,] # keep records with complete values
nrow(occurrences)

## 5.3 Extract predictor values for PA/BCK locations (check for potential NA)
PA.BKG.data <- PA.BKG.data[!is.na(PA.BKG.data$BU),]
BU.pred <- raster::extract(covariates, PA.BKG.data)
BU.pred <- as.data.frame(BU.pred)

## 5.4 Cross-check for potential NA values that can make the modelling run crashed out.
nl <- which(complete.cases(BU.pred)) # 2073  PA/BCK complete data
PA.BKG.data <- PA.BKG.data[nl,] # keep records with complete values 
nrow(PA.BKG.data)

# 6 Check that no points are placed in the same pixel.
##  6.1 Identify records which are located a distance lower than 5000m
PA.BKG.data1 <- as(PA.BKG.data, "Spatial")
head(PA.BKG.data1)
i <- zerodist(PA.BKG.data1, zero = 5000)
i <- as.vector(i)
i <- unique(i)
length(i)

## 6.2 Remove those that are closer than 5000m and are background points
j <- which(PA.BKG.data1@data$BU == 0 & PA.BKG.data1@data$type != "BG")
k <- which(i %in% j)
i <- i[k] ## remove records from pseudoabsences located at a distance lower than 5km from either backgrounds or pseudabs

head(PA.BKG.data1)

PA.BKG.data1 <- PA.BKG.data1[-i,] # 2 records with this criteria
nrow(PA.BKG.data1)

# 5. Set up a vector of weights to be used when fitting environmental models
PA.BKG.data$weight <- 0
PA.BKG.data1 <- PA.BKG.data1[,-3]

PA.BKG.data <- st_as_sf(PA.BKG.data1)
AFRO_BU <- rbind(occurrences, PA.BKG.data)

## 5.1 set weight to 1 for occ and 0.5 for PA / BG
AFRO_BU[AFRO_BU$BU == 1,]$weight <- 1
AFRO_BU[AFRO_BU$BU == 0,]$weight <- 0.5

## 5.2 Coerce weight column to vector (to be passed to the Biomod Modelling function)
weights <- as.vector(AFRO_BU$weight)
head(AFRO_BU)
AFRO_BU <- AFRO_BU[1]

# 6. Preparation of final dataset to be used in the modelling
AFRO_BU <- as(AFRO_BU, "Spatial")

plot(covariates[[1]])
plot(AFRO_BU, add = T)


## ~~~~~~~~~~~~~~~~~~~~~~~~~ Run Ecological Modelling using biomod2 package ~~~~~~~~~~~ #
# source(paste(path_scripts,"BIOMOD_FormatingData.R", sep = "/"))

BU.Formatted <- BIOMOD_FormatingData(resp.var = AFRO_BU, 
                                        expl.var = covariates,
                                        resp.name = "BU")

jar <- paste(path_maxent, "maxent.jar", sep='/')

myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = jar))

# 7. Develop a fist model based on the following algorithms: # "GLM","GAM","ANN","MARS", "RF","MAXENT.Phillips"

setwd(path_outputs)
getwd()
set.seed(123)
system.time(BU.model.PCR <- BIOMOD_Modeling(BU.Formatted,
                                                      models = c("GLM","GAM","GBM","ANN","MARS", "RF","MAXENT.Phillips"),
                                                      models.options = myBiomodOption,
                                                      NbRunEval = 20,                  # run at least 20 single models
                                                      DataSplit = 80,                  # splitting training/testing every round
                                                      Prevalence = NULL,  
                                                      Yweights = weights,
                                                      VarImport = 3,
                                                      models.eval.meth = c("TSS","ROC","ACCURACY"),
                                                      SaveObj = TRUE,                  # keep results on hard drive
                                                      rescal.all.models = FALSE,       # 
                                                      do.full.models = FALSE,          # 
                                                      modeling.id = paste("BU","ENM1",sep =".")))

## system elapsed: 4698.65 (78.3 min)

BU.model.PCR # Model failing: none

## 7.1. Extract evaluations from all models run
BU.model.PCR.Eval <- get_evaluations(BU.model.PCR)
dimnames(BU.model.PCR.Eval)

## Visualize performance of the models 
## Call a for modified performance graphs.
source(paste(path_scripts,"models_scores_graph.R", sep = "/"))

setwd(path_outputs)

theme_set(theme_grey())
performance.models1 <- models_scores_graph(BU.model.PCR, 
                                           by = "models", 
                                           metrics = c("ROC","TSS"))

performance.models1 +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size= 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12))

dev.print(file="V1_Performance_models.png", device=png, height=600, width=900)

dev.off()

performance.models2 <- models_scores_graph(BU.model.PCR, 
                                           by = "models", 
                                           metrics = c("ROC","ACCURACY"))
performance.models2 +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size= 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12))

dev.print(file="V2_Performance_models.png", device=png, height=600, width=900)
dev.off()

## 7.2. Extract TSS/ROC/ACCURACY for test data (TSS: True Skill Statistic)

TSS.test.data <- as.data.frame(t(BU.model.PCR.Eval["TSS","Testing.data",,,]))
ROC.test.data <- as.data.frame(t(BU.model.PCR.Eval["ROC","Testing.data",,,]))
PCC.test.data <- as.data.frame(t(BU.model.PCR.Eval["ACCURACY","Testing.data",,,]))

summary(TSS.test.data)
summary(ROC.test.data)
summary(PCC.test.data)

## TSS

TSS.Mean <- apply(TSS.test.data, 2, function (x) mean(x, na.rm = TRUE)) 
TSS.Median <- apply(TSS.test.data, 2, function (x) median(x, na.rm = TRUE))
TSS.1Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.25), type = 7, na.rm = TRUE))
TSS.3Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.75), type = 7, na.rm = TRUE))

TSS <- t(data.frame(TSS.Mean, TSS.Median, TSS.1Q, TSS.3Q))

## ROC

ROC.Mean <- apply(ROC.test.data, 2, function (x) mean(x, na.rm = TRUE)) 
ROC.Median <- apply(ROC.test.data, 2, function (x) median(x, na.rm = TRUE))
ROC.1Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.25), type = 7, na.rm = TRUE))
ROC.3Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.75), type = 7, na.rm = TRUE))

ROC <- t(data.frame(ROC.Mean, ROC.Median, ROC.1Q, ROC.3Q))

## PCC
## Accuracy: PCC (Positive Correcly Classified): % presences which are correctly predicted.

PCC.Mean <- apply(PCC.test.data, 2, function (x) mean(x, na.rm = TRUE)) 
PCC.Median <- apply(PCC.test.data, 2, function (x) median(x, na.rm = TRUE))
PCC.1Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.25), type = 7, na.rm = TRUE))
PCC.3Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.75), type = 7, na.rm = TRUE))

PCC <- t(data.frame(PCC.Mean, PCC.Median, PCC.1Q, PCC.3Q))

## Save the outputs

write.csv(TSS, file = "TSS_BU.csv")
write.csv(ROC, file = "ROC_BU.csv")
write.csv(PCC, file = "PCC_BU.csv")

rm(TSS.Mean, TSS.Median, TSS.1Q, TSS.3Q, ROC.Mean, ROC.Median, ROC.1Q, ROC.3Q,
   PCC.Median, PCC.Mean, PCC.1Q, PCC.3Q)

## 7.3. Get variable importance for the different models

BU.model.PCR.Var <- get_variables_importance(BU.model.PCR)
dimnames(BU.model.PCR.Var)

Variable.Contribution <- t(apply(BU.model.PCR.Var, 2, function (x) rowMeans(x, na.rm = TRUE)))

write.csv(Variable.Contribution, file = "Variable_Contribution_BU_PCR.csv")

## 7.4. Explore variable response curves for RF models
setwd(path_scripts)
## 1. Call a function to plot response curves
source("Response_plot.R")

setwd(path_outputs)

BU.RF <- BIOMOD_LoadModels(BU.model.PCR, models = "RF")

RFsPlot2D <- response.plot2(models = BU.RF, 
                            Data = get_formal_data(BU.model.PCR,"expl.var"), 
                            show.variables = get_formal_data(BU.model.PCR, "expl.var.names"), 
                            do.bivariate = FALSE, 
                            fixed.var.metric = "median", 
                            col = c("blue", "red"),
                            legend = TRUE, 
                            data_species = get_formal_data(BU.model.PCR, "resp.var"))

rm(list = ls(pattern = "BU_AllData")) 


select_models <-  BU.model.PCR@models.computed
# 8. Ensemble Modelling. Come up with a consensus model and ensemble for every type of model
BU.model.PCR.Ensemble <- BIOMOD_EnsembleModeling(modeling.output = BU.model.PCR,
                                                           chosen.models = select_models,
                                                           em.by = "all",
                                                           eval.metric = c("ROC"),
                                                           eval.metric.quality.threshold = c(0.8), 
                                                           prob.mean = T,
                                                           prob.cv = T,
                                                           prob.ci = T,
                                                           prob.ci.alpha = 0.05,
                                                           prob.median = T,
                                                           committee.averaging = T,
                                                           prob.mean.weight = T,
                                                           prob.mean.weight.decay = "proportional")

## 8.1. Extract evaluation of the ensemble model

Ensemble.Algorithms <- c("Mean","Coef.Var","InfCI","SupCI","Median","MCAvg","WeightedMean")

BU.Ensemble.Eval <- get_evaluations(BU.model.PCR.Ensemble)

BU.Ensemble.TSS <- t(sapply(BU.Ensemble.Eval, "[",2,1:4))
rownames(BU.Ensemble.TSS) <- Ensemble.Algorithms

BU.Ensemble.KAPPA <- t(sapply(BU.Ensemble.Eval, "[",1,1:4))
rownames(BU.Ensemble.KAPPA) <- Ensemble.Algorithms

BU.Ensemble.ROC <- t(sapply(BU.Ensemble.Eval, "[",3,1:4))
rownames(BU.Ensemble.ROC) <- Ensemble.Algorithms

cutoff <- (BU.Ensemble.ROC[1,2])/1000

setwd(path_outputs)
write.csv(BU.Ensemble.TSS, file = "BU_Ensemble_TSS.csv")
write.csv(BU.Ensemble.KAPPA, file = "BU_Ensemble_KAPPA.csv")
write.csv(BU.Ensemble.ROC, file = "BU_Ensemble_ROC.csv")

# 9. Final projection (extrapolation) based on the best fitted model. 
## Space projection for the algorithms selected

BU.Proj <- BIOMOD_Projection(modeling.output = BU.model.PCR,
                              new.env = covariates,
                              proj.name = "Final",
                              selected.models = select_models,
                              binary.meth = "ROC",
                              compress = "gzip",
                              build.clamping.mask = FALSE,
                              output.format = ".grd")

BU.projections <- get_predictions(BU.Proj) # Set up a raster stack with the predictions

# 10. Ensemble forescasting combining projections based on models ensemble rules defined at the ensemble modelling step (combine single surfaces based on ensemble parameters)

BU.EF <- BIOMOD_EnsembleForecasting(EM.output = BU.model.PCR.Ensemble,
                                     projection.output = BU.Proj)

BU.Ensemble.final <- get_predictions(BU.EF)

names(BU.Ensemble.final)

# 10. Extract best performing models: mean of probabilities with CI and models committe averaging
# probability raster are provided scaled 0 to 1000, so we need to rescale them 0 to 1
# This is made because float and double raster are heavier (and take up more disk space)

BU.mean <- raster(BU.Ensemble.final, layer=1)/1000 
BU.LB <- raster(BU.Ensemble.final, layer=3)/1000
BU.UB <- raster(BU.Ensemble.final, layer = 4)/1000
BU.median <- raster(BU.Ensemble.final, layer = 5)/1000
BU.MCA <- raster(BU.Ensemble.final, layer = 6)/1000
BU.WM <- raster(BU.Ensemble.final, layer = 7)/1000

getwd()
par(mfrow=c(1,1))

plot(BU.WM)
plot(AFRO_map, add = T)

# 11. Export resulting raster datasets
setwd(path_wd)
setwd(path_outputs)

writeRaster(BU.mean,filename = "BU_Mean.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.LB,filename = "BU_LB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.UB,filename = "BU_UB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.median,filename = "BU_median.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.MCA,filename = "BU_MCA.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.WM,filename = "BU_WM.tif", format= "GTiff", overwrite=TRUE)

# 11. Generate binary surfaces (probable occurrence) from environmental suitability surfaces
## We can generate a binary surfaces using the optimal cutoff estimated from best tradeoff 
## for sensitivity and specificity as defined by the selected indicator of predictive
## performance (Receiver Operating Characteristic, ROC).
getwd()
m <- c(0, cutoff, 0, cutoff, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

setwd(path_outputs)

BU.WM <- raster("BU_WM.tif")
BU.UB <- raster("BU_UB.tif")
BU.LB <- raster("BU_LB.tif")

BU.Occ <- reclassify(BU.WM, rclmat)
BU.Occ <- round(BU.Occ,0)

BU.Occ.LB <- reclassify(BU.LB, rclmat)
BU.Occ.LB <- round(BU.Occ.LB,0)

BU.Occ.UB <- reclassify(BU.UB, rclmat)
BU.Occ.UB <- round(BU.Occ.UB,0)

dev.off()

getwd()
setwd(path_outputs)
## Export resulting raster datasets
writeRaster(BU.Occ,filename = "BU_Occurrence.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.Occ.LB,filename = "BU_Occurrence_LB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(BU.Occ.UB,filename = "BU_Occurrence_UB.tif", format= "GTiff", overwrite=TRUE)

par(mfrow=c(1,3))
plot(BU.Occ.LB)
plot(BU.Occ)
plot(BU.Occ.UB)


### error
BU.CI <- BU.UB - BU.LB 
BU.error <- sdmvspecies::rescale(BU.CI)

par(mfrow=c(1,2))
plot(BU.WM)
plot(BU.error)

setwd(path_outputs)
writeRaster(BU.error, filename = "BU.error_PCR.tif", format= "GTiff", overwrite=TRUE)



# ~~~~~~~~~~~~~ Preparing Marginal Effect Plots for Random Forest based models ~~~~~~~~~~~~~~ ##

## 12. Creating marginal effect plots for covariates modelled using RF

names(RFsPlot2D)
## AnPET
AnPET <- gather(RFsPlot2D$AnPET, 
                Model, 
                Marginal_Effect, 
                BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

AnPET.plot <- ggplot(AnPET, aes(x= AnPET, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(1000,3000)) +
    stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Potential Evapotranspiration (mm/month)") +
  theme_bw() + theme(axis.text=element_text(size=14),
                     axis.title=element_text(size=16))
AnPET.plot


## AnPET
mintemp <- gather(RFsPlot2D$mintemp, 
                Model, 
                Marginal_Effect, 
                BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

mintemp.plot <- ggplot(mintemp, aes(x= mintemp, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
#  scale_x_continuous(limits = c(1000,3000)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Minimum temperature (°C)") +
  theme_bw() + theme(axis.text=element_text(size=14),
                     axis.title=element_text(size=16))
mintemp.plot
names(covariates)

### bioclim_15 = Precipitation Seasonality (Coefficient of Variation)
Prec_ssn <- gather(RFsPlot2D$bc_15_pr_ssn, 
                 Model, 
                 Marginal_Effect, 
                 BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

Prec_ssn.plot <- ggplot(Prec_ssn, aes(x= bc_15_pr_ssn, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0, 100)) +
  
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Precipitation Seasonality (%)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))
Prec_ssn.plot

names(covariates)

### Precipitation of Wettest Quarter
BC_16 <- gather(RFsPlot2D$BC_16_precip_wettst_q, 
                  Model, 
                  Marginal_Effect, 
                  BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

BC_16.plot <- ggplot(BC_16, aes(x= BC_16_precip_wettst_q, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
scale_x_continuous(limits = c(0, 1500)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Precipitation in the Wettest Quarter (mm)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

BC_16.plot

### Precipitation of Driest Quarter
bio17 <- gather(RFsPlot2D$bio17, 
                     Model, 
                     Marginal_Effect, 
                     BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

bio17.plot <- ggplot(bio17, aes(x= bio17, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0, 200)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Precipitation in the Driest Quarter (mm)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

bio17.plot


### Distance to nearest water body (km)
DistWB <- gather(RFsPlot2D$ED_OSM_WBS1, 
                 Model, 
                 Marginal_Effect, 
                 BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

DistWB.plot <- ggplot(DistWB, aes(x= ED_OSM_WBS1, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0, 100)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Distance to nearest waterbody (km)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

DistWB.plot

names(covariates)

### Distance to deforested areas (km)
DistForest <- gather(RFsPlot2D$EucDist_DeforestedAreas_Km, 
                     Model, 
                     Marginal_Effect, 
                     BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

DistForest.plot <- ggplot(DistForest, aes(x= EucDist_DeforestedAreas_Km, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0, 100)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Distance to deforested area (km)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

DistForest.plot

names(covariates)


### Distance to nearest river (km)
Distrivers <- gather(RFsPlot2D$EucDist_km_Rivers1, 
                     Model, 
                     Marginal_Effect, 
                     BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

Distrivers.plot <- ggplot(Distrivers, aes(x= EucDist_km_Rivers1, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,100)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Distance to nearest river or stream (km)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

Distrivers.plot

### Distance to nearest dam (km)
Dist_dams <- gather(RFsPlot2D$dist_dam1, 
                     Model, 
                     Marginal_Effect, 
                     BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

Dist_dams.plot <- ggplot(Dist_dams, aes(x= dist_dam1, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,100)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Distance to nearest dam (km)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

Dist_dams.plot

### Distance to nearest agricultural area (km)
Dist_agr <- gather(RFsPlot2D$eu_dist_agr, 
                     Model, 
                     Marginal_Effect, 
                     BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

Dist_agr.plot <- ggplot(Dist_agr, aes(x= eu_dist_agr, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,100)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Distance to nearest agricultural area (km)") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

Dist_agr.plot

### EVI
EVI <- gather(RFsPlot2D$evi_5k, 
              Model, 
              Marginal_Effect, 
              BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

EVI.plot <- ggplot(EVI, aes(x= evi_5k, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +  
  scale_x_continuous(limits = c(0,0.5)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Enhanced Vegitation Index*") +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))

EVI.plot


names(covariates)

setwd(path_outputs)
## Export a figure with the marginal effect plots for the selected covariates
plot_grid( DistWB.plot,  AnPET.plot, mintemp.plot, Distrivers.plot, BC_16.plot, EVI.plot, Prec_ssn.plot, bio17.plot,   DistForest.plot,    Dist_dams.plot,  Dist_agr.plot, 
           labels = c("A","B","C","D", "E", "F","G","H","I", "J", "K"), ncol = 3, nrow = 4, label_size = 10)

setwd(path_outputs)
dev.print(file="Marginal_Plots_RF.png", device=png, height = 1200, width = 2000)
dev.off()
getwd()


# Clean up workspace of unnecessary objects and save objects created in the corresponding folder
gc()
save.image("BU_ENVM_Modelling.RData")

