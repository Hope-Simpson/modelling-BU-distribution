# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~ Modelling Spatial Distribution of Mycobacterium ulcerans across Africa ~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## Start Date: 10.09.2019
## End Date: 10.11.2020
## Place: London, UK
## Project: Delineating environmental limits for Buruli ulcer across Africa
## Modelling approach: Environmental modelling using occurrence of Buruli ulcer and 
## including pseudoabsence points from countries with evidence consensus for absence
## Background points created randomly around presences to account for potential
## geographical bias

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
path_outputs <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/Outputs/MU_environmental"
path_scripts <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling/final_scripts"

## 1. Call a function to install and load multiple R packages.
setwd(path_wd)
source(paste(path_scripts,"ipak.R", sep = "/"))

# List of packages
packages <- c("tidyverse","reshape2","stats","readxl","foreign",
              "sf","dismo","raster","rgdal","proj4","ggplot2","cowplot",
              "biomod2","rJava","pROC","pryr","RStoolbox","vegan",
              "psycho","corrplot","spatialEco","maptools","gbm","spatstat")
ipak(packages)

# 3. Import data 
## 3.1 Load the Africa map with country boundaries
AFRO_map <-  readOGR(dsn=path_shapefiles, layer="salb")

## 3.2 Load the data points: occurrence of MU 
AFRO_MU <- file.path(path_shapefiles, "/BU_occurrences/env_occurrences.shp")
AFRO_MU <- st_read(AFRO_MU)
AFRO_MU <- st_transform(AFRO_MU, crs = 3857)
names(AFRO_MU)
occurrences <- AFRO_MU[,"geometry"]
head(occurrences) ## 79 occurrences recorded

occurrences$MU <- 1
occurrences$weight <- 1

names(occurrences)
occurrences <- occurrences[,c(2,3,1)]
names(occurrences) <- c("MU", "weight", "geometry")

plot(AFRO_map)
plot(occurrences$geometry, cex = .2, col = "red", add = TRUE)

### Load the pseudo-absences points generated (6_create_PA_dataset_SRE_env script) 
pseudo_abs <- file.path(path_BG, "pseudoabs_env.shp")
pseudo_abs <- st_read(pseudo_abs)
pseudo_abs <- st_transform(pseudo_abs, crs = 3857)
pseudo_abs$type <- "PA"
pseudo_abs$MU <- 0
pseudo_abs <- pseudo_abs[,c(1,2,4,3)]

nrow(pseudo_abs) ## 79 pa points

plot(pseudo_abs$geometry, cex = .2, col = "blue", add = TRUE)

### Load the background points
BG <- file.path(path_BG, "KDE_BG_env.shp")
BG <- st_read(BG)
BG <- st_transform(BG, crs = st_crs(AFRO_map))

nrow(BG) ## 
BG$weight <- 0 # assign weight of 0 to background data
BG$type <- "BG"

BG <- BG[,c(1,3,4,2)]
plot(pseudo_abs$geometry, cex = .2, col = "grey", add = TRUE)

### Assemble background and pseudo-absences
PA.BKG.data <- rbind(pseudo_abs, BG)
names(PA.BKG.data)
nrow(PA.BKG.data) # 158 records
PA.BKG.data <- st_transform(PA.BKG.data, crs = 3857)

## 3.3 Load the covariates to be used in the modelling
path_covars <- "C:/Users/Hope/Dropbox (LSoHaTM)/BU_Modelling"
cov.5km <- file.path(path_covars, "covariates_5km_env.rds")
cov.5km <- readRDS(cov.5km)
PCS = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"

names(cov.5km)

## convert distance to dams + agr areas from m to km
cov.5km$eu_dist_agr <- cov.5km$eu_dist_agr / 1000
cov.5km$dist_dam1 <- cov.5km$dist_dam1 / 1000

covariates <- stack()
for(i in 3:length(cov.5km)){
  raster1 <- rasterFromXYZ(cov.5km[,c(1,2,i)],
                           crs = PCS)
  covariates <- stack(covariates, raster1)
}

names(covariates)

rm(i, raster1)

plot(covariates)

# 5. Prepare occurrence and background points datasets for being use in the model
## 5.1 Extract predictor values for occurrence locations (check for potential NA)
MU.pred <- raster::extract(covariates, occurrences)
MU.pred <- as.data.frame(MU.pred)

## 5.2 Cross-check for  NA values that can make the model crash.
nl <- which(complete.cases(MU.pred)) # ## 687 out of 690 observations completed data
occurrences <- occurrences[nl,] # keep those records with completed values (12 have na)

## 5.3 Extract predictor values for PA/BCK locations (check for potential NA)
MU.pred2 <- raster::extract(covariates, PA.BKG.data)
MU.pred2 <- as.data.frame(MU.pred2)

## 5.4 Cross-check for potential NA values that can make the modelling run crashed out.
nl <- which(complete.cases(MU.pred2)) # 157 PA/BCK completed data
PA.BKG.data <- PA.BKG.data[nl,] # keep those records with completed values 

# 6 Check that no points are placed in the same pixel.
##  6.1 Identify records which are located a distance lower than 5000m
PA.BKG.data1 <- as(PA.BKG.data, "Spatial")
head(PA.BKG.data1)
i <- zerodist(PA.BKG.data1, zero = 5000)
i <- as.vector(i)
i <- unique(i)

## 6.2 Remove those that are closer than 5000m and are background points
j <- which(PA.BKG.data1@data$MU == 0 & PA.BKG.data1@data$type != "BG")
k <- which(i %in% j)
i <- i[k] ## remove records from pseudoabsences located at a distance lower than 5km from either backgrounds or pseudabs
# PA.BKG.data1 <- PA.BKG.data1[-i,] # 

names(PA.BKG.data1)

# 5. Set up a vector of weights to be used when fitting environmental models
PA.BKG.data1$weight <- 0.5
PA.BKG.data1 <- PA.BKG.data1[,-3]

PA.BKG.data <- st_as_sf(PA.BKG.data1)
AFRO_MU <- rbind(occurrences, PA.BKG.data)

## 5.1 Checking that total weight for presences and for PA/BCK sum up the same.
names(AFRO_MU)

sum(AFRO_MU[AFRO_MU$MU == 0,]$weight)
sum(AFRO_MU[AFRO_MU$MU == 1,]$weight)

## 5.2 Coerce weight column to vector (to be passed to the Biomod Modelling function)
# then drop weights from occurrence data
weights <- as.vector(AFRO_MU$weight)
head(AFRO_MU)
AFRO_MU <- AFRO_MU[,-2]

# 6. Preparation of final datasets to be used in the modelling
AFRO_MU1 <- as(AFRO_MU, "Spatial")
nrow(AFRO_MU1)

plot(covariates[[1]])
plot(AFRO_MU1, add = T)

## ~~~~~~~~~~~~~~~~~~~~~~~~~ Run Ecological Modelling using biomod2 package: no BCG points ~~~~~~~~~~~ #
#  Compiling the different elements to be using when computing models
#  https://www.rdocumentation.org/packages/biomod2/versions/3.3-7.1/topics/BIOMOD_Modeling

## define majorlanduses and urbchx as factors
setwd(path_outputs)
MU.Formatted <- BIOMOD_FormatingData(resp.var = AFRO_MU1, 
                                        expl.var = covariates,
                                        resp.name = "BU")

jar <- paste(path_maxent, "maxent.jar", sep='/')

myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = jar))

# 7. Develop a fist model based on the following algorithms: # "GLM","GAM","ANN","SRE","MARS","RF","MAXENT.Phillips"

system.time(MU.model1ENV <- BIOMOD_Modeling(MU.Formatted,
                                            models = c("GLM","GAM","GBM","ANN", "MARS","RF","MAXENT.Phillips"),
                                         models.options = myBiomodOption,
                                         NbRunEval = 20,                  # run at least 20 single models
                                         DataSplit = 80,                  # spliting training/testing every round
                                         Prevalence = NULL,  
                                         Yweights = weights,
                                         VarImport = 3,
                                         models.eval.meth = c("TSS","ROC","ACCURACY"),
                                         SaveObj = TRUE,                  # keep results on hard frive
                                         rescal.all.models = FALSE,       # 
                                         do.full.models = FALSE,          # 
                                         modeling.id = paste("BU","ENM1",sep =".")))

getwd()
## system elapsed: 4698.65 (78.3 min)

MU.model1ENV # Model failing 0

## 7.1. Extract evaluations from all models run
MU.model1ENV.Eval <- get_evaluations(MU.model1ENV)
dimnames(MU.model1ENV.Eval)


## Visualize performance of the models 
## Call a for modified performance graphs.
source(paste(path_scripts,"models_scores_graph.R", sep = "/"))

setwd(path_outputs)
theme_set(theme_grey())
performance.models1 <- models_scores_graph(MU.model1ENV, 
                                           by = "models", 
                                           metrics = c("ROC","TSS"))

performance.models1 +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size= 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12))

dev.print(file="Performance_models.png", device=png,  height=600, width=900)
dev.off()

performance.models2 <- models_scores_graph(MU.model1ENV, 
                                           by = "models", 
                                           metrics = c("ROC","ACCURACY"))
performance.models2 +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size= 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12))

dev.print(file="Performance_models2.png", device=png,  height=600, width=900)
dev.off()

## 7.2. Extract TSS/ROC/ACCURACY for test data (TSS: True Skill Statistic)

TSS.test.data <- as.data.frame(t(MU.model1ENV.Eval["TSS","Testing.data",,,]))
ROC.test.data <- as.data.frame(t(MU.model1ENV.Eval["ROC","Testing.data",,,]))
PCC.test.data <- as.data.frame(t(MU.model1ENV.Eval["ACCURACY","Testing.data",,,]))

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

write.csv(TSS, file = "TSS_MU.csv")
write.csv(ROC, file = "ROC_MU.csv")
write.csv(PCC, file = "PCC_MU.csv")

rm(TSS.Mean, TSS.Median, TSS.1Q, TSS.3Q, ROC.Mean, ROC.Median, ROC.1Q, ROC.3Q,
   PCC.Median, PCC.Mean, PCC.1Q, PCC.3Q)

## 7.3. Get variable importance for the different models

MU.Model1_ENV2.Var <- get_variables_importance(MU.model1ENV)
dimnames(MU.Model1_ENV2.Var)

Variable.Contribution <- t(apply(MU.Model1_ENV2.Var, 2, function (x) rowMeans(x, na.rm = TRUE)))

write.csv(Variable.Contribution, file = "Variable_Contribution_MU_model1.csv")

## 7.4. Explore variable response curves for RF models
# setwd("C:/Users/Hope/Documents")
setwd(path_scripts)
## 1. Call a function to install and load multiple R packages.
source("Response_plot.R")

setwd(path_outputs)

MU1.RF <- BIOMOD_LoadModels(MU.model1ENV, models = "RF")

RFsPlot2D <- response.plot2(models = MU1.RF, 
                            Data = get_formal_data(MU.model1ENV,"expl.var"), 
                            show.variables = get_formal_data(MU.model1ENV, "expl.var.names"), 
                            do.bivariate = FALSE, 
                            fixed.var.metric = "median", 
                            col = c("blue", "red"),
                            legend = TRUE, 
                            data_species = get_formal_data(MU.model1ENV, "resp.var"))



rm(list = ls(pattern = "BU_AllData")) 

# 8. Ensemble Modelling. Come up with a consensus model and ensemble for every type of model
select_models <- MU.model1ENV@models.computed

MU.model1ENV.Ensemble <- BIOMOD_EnsembleModeling(modeling.output = MU.model1ENV,
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

MU1.Ensemble.Eval <- get_evaluations(MU.model1ENV.Ensemble)

MU1.Ensemble.TSS <- t(sapply(MU1.Ensemble.Eval, "[",2,1:4))
rownames(MU1.Ensemble.TSS) <- Ensemble.Algorithms

MU1.Ensemble.KAPPA <- t(sapply(MU1.Ensemble.Eval, "[",1,1:4))
rownames(MU1.Ensemble.KAPPA) <- Ensemble.Algorithms

MU1.Ensemble.ROC <- t(sapply(MU1.Ensemble.Eval, "[",3,1:4))
rownames(MU1.Ensemble.ROC) <- Ensemble.Algorithms

cutoff <- (MU1.Ensemble.ROC[1,2])/1000


write.csv(MU1.Ensemble.TSS, file = "MU1_Ensemble_TSS.csv")
write.csv(MU1.Ensemble.KAPPA, file = "MU1_Ensemble_KAPPA.csv")
write.csv(MU1.Ensemble.ROC, file = "MU1_Ensemble_ROC.csv")

# 9. Final projection (extrapolation) based on the best fitted model. 
## Space projection for the algorithms selected

MU1.Proj <- BIOMOD_Projection(modeling.output = MU.model1ENV,
                              new.env = covariates,
                              proj.name = "Final",
                              selected.models = select_models,
                              binary.meth = "ROC",
                              compress = "gzip",
                              build.clamping.mask = FALSE,
                              output.format = ".grd")

MU1.projections <- get_predictions(MU1.Proj) # Set up a raster stack with the predictions

# 10. Ensemble forescasting combining projections based on models ensemble rules defined at the ensemble modelling step (combine single surfaces based on ensemble parameters)

MU1.EF <- BIOMOD_EnsembleForecasting(EM.output = MU.model1ENV.Ensemble,
                                     projection.output = MU1.Proj)

MU1.Ensemble.final <- get_predictions(MU1.EF)

names(MU1.Ensemble.final)

# 10. Extract best performing models: mean of probabilities with CI and models committe averaging
# probability raster are provided scaled 0 to 1000, so we need to rescale them 0 to 1
# This is made because float and double raster are heavier (and take up more disk space)

MU1.mean <- raster(MU1.Ensemble.final, layer=1)/1000 
MU1.LB <- raster(MU1.Ensemble.final, layer=3)/1000
MU1.UB <- raster(MU1.Ensemble.final, layer = 4)/1000
MU1.median <- raster(MU1.Ensemble.final, layer = 5)/1000
MU1.MCA <- raster(MU1.Ensemble.final, layer = 6)/1000
MU1.WM <- raster(MU1.Ensemble.final, layer = 7)/1000

par(mfrow=c(1,3))

plot(MU1.LB)
plot(MU1.WM)
plot(MU1.UB)

par(mfrow=c(1,1))

plot(MU1.mean)
plot(AFRO_map, add = T)

setwd(path_outputs)
# 11. Export resulting raster datasets
writeRaster(MU1.mean,filename = "MU1_Mean.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.LB,filename = "MU1_LB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.UB,filename = "MU1_UB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.median,filename = "MU1_median.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.MCA,filename = "MU1_MCA.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.WM,filename = "MU1_WM.tif", format= "GTiff", overwrite=TRUE)

# 11. Generate binary surfaces (probable occurrence) from environmental suitability surfaces
## We can generate a binary surfaces using the optimal cutoff estimated from best tradeoff 
## for sensitivity and specificity as defined by the selected indicator of predictive
## performance (Receiver Operating Characteristic, ROC).

m <- c(0, cutoff, 0, cutoff, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

MU1.Occ <- reclassify(MU1.WM, rclmat)
MU1.Occ <- round(MU1.Occ,0)
plot(MU1.Occ)

MU1.Occ.LB <- reclassify(MU1.LB, rclmat)
MU1.Occ.LB <- round(MU1.Occ.LB,0)
plot(MU1.Occ.LB)

MU1.Occ.UB <- reclassify(MU1.UB, rclmat)
MU1.Occ.UB <- round(MU1.Occ.UB,0)
plot(MU1.Occ.UB)

## Export resulting raster datasets
setwd(path_outputs)
writeRaster(MU1.Occ,filename = "MU1_Occurrence.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.Occ.LB,filename = "MU1_Occurrence_LB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(MU1.Occ.UB,filename = "MU1_Occurrence_UB.tif", format= "GTiff", overwrite=TRUE)
getwd()

### error

plot(MU1.UB)

MU.CI <-  MU1.UB - MU1.LB
MU.error <- sdmvspecies::rescale(MU.CI)

par(mfrow=c(1,2))
plot(MU1.WM)
plot(MU.error)

setwd(path_outputs)
writeRaster(MU.error, filename = "MU.error_env_all.tif", format= "GTiff", overwrite=TRUE)
getwd()


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


## min temp
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

# DistWB$ED_OSM_WBS1 <- DistWB$ED_OSM_WBS1 * 1000

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



### Wetness Index1
Wetness_Index1 <- gather(RFsPlot2D$Wetness_Index1, 
                         Model, 
                         Marginal_Effect, 
                         BU_AllData_RUN1_RF:BU_AllData_RUN20_RF)

Wetness_Index1.plot <- ggplot(Wetness_Index1, aes(x= Wetness_Index1, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(9,18.5)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Response Curve") +
  xlab("Wetness Index") +
  theme_bw()+ theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=18))

Wetness_Index1.plot

setwd(path_outputs)
## Export a figure with the marginal effect plots for the selected covariates
plot_grid( DistForest.plot,  AnPET.plot,  Dist_dams.plot, mintemp.plot, bio17.plot,  EVI.plot, Distrivers.plot, DistWB.plot,  Dist_agr.plot,  Wetness_Index1.plot,
           labels = c("A","B","C","D", "E", "F","G","H","I", "J"), ncol = 3, nrow = 4, label_size = 10)

setwd(path_outputs)
dev.print(file="Marginal_Plots_RF.png", device=png, height = 1200, width = 2000)
dev.off()
getwd()







# Clean up workspace of unnecessary objects and save objects created in the corresponding folder
gc()
save.image(paste0(path_outputs, sep = "/", "BU_ENVM_Modelling_v01.RData"))
 
