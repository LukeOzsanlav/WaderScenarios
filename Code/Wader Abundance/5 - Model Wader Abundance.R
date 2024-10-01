##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 02/01/24
## 
## Goal: Model drivers of wader abundance using a Random Forest
##
##------------------------------------------------------##

## Updates:
## - try two step modelling approach, see here: https://link.springer.com/article/10.1007/s00607-023-01224-3
## Alternative Spatial RF package
## https://blasbenito.github.io/spatialRF/
## See this paper for details on spatial RF https://link.springer.com/article/10.1007/s11004-021-09946-w
## - I could round the Lapwing and Redshank estimates to the nearest half pair? This might get rid of the problems of lots of 0.1 pairs everywhere

## Load in required packages
pacman::p_load(here, tidyverse, data.table, randomForest, caret, randomForestSRC, interp, sf, gridExtra)
options(scipen = 100, digits = 4) # set notation
source("Code/Helper functions.R") # get helper functions

## Use `here` package for specifying file paths
here::i_am("Code/Wader Abundance/5 - Model Wader Abundance.R")
library(here)
here()



##----------------------##
#### 0.1 Data read in ####
##----------------------##

## Read in filtered breeding pairs estimates
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")

## Read in data on field characteristics from survey
FieldChar <- read_csv("CleanData/Wader Abundance/1-CleaningScript/Field_Characteristics_Clean.csv")

## read in BWWM field shapefiles
BWWM <- st_read("RawData/BWWM field shapefile/BWWM_fields_27Sep2023.shp") |> select(F_LOC_ID)



##------------------------##
#### 0.2 Join data sets ####
##------------------------##

## Join together the Wader data sets and the field characteristics
WData <- left_join(Waders, FieldChar, by = c("F_LOC_ID", "year"))
summary(WData)



##---------------------------##
#### 0.3 Add extra columns ####
##---------------------------##

## Add additional columns to the data set for modelling
WData <- WData |> 
         mutate(FieldArea = FieldArea/10000, # convert field area to hectares
                CorvDens = MagDens + CrowDens, # sum corvid density
                Reserve = ifelse(RSPB_Reserve == "Y" | LNR == "Y" | NNR == "Y", "Y", "N"), # reserve indicator variable
                RUSH_PERCENT = ifelse(is.na(RUSH_PERCENT)==T, 0, RUSH_PERCENT), # change NA's in rush cover to 0
                Fence_Coverage = ifelse(Fence_Coverage == 0, "N", "Y"), # fence indicator variable
                Lap_Density = est_pairsL/FieldArea, # breeding Lapwing density
                Red_Density = est_pairsR/FieldArea, # breeding Redshank density
                Sni_Density = est_pairsS/FieldArea,)# breeding Snipe density

## Plot histogram of breeding bird densities
hist(WData$Lap_Density, main = "Breeding Lapwing density")
hist(WData$Red_Density, main = "Breeding Redshank density")
hist(WData$Sni_Density, main = "Breeding Snipe density")



##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##



##--------------------------------##
#### 1.0 Lapwing One-step Model ####
##--------------------------------##

##------------------------##
#### Lapwing: Prep Data ####
##------------------------##

## Filter out the Lapwing survey data
Lap <- filter(WData, is.na(est_pairsL)==F)
nrow(Lap); nrow(filter(Lap, est_pairsL > 0)) # total no. of survey fields and total with waders 

## Check if survey columns have NA's
colnames(Lap) # get column names
Lap |> select(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, GROUND_DAMP, STOCK, VEG_HEIGHT, 
              VEG_STRUCTURE, RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE, WaterCoverage, Fence_Coverage) |> 
       is.na() |> summary()

## remove rows with with NAs in these columns
LapF <- Lap |> drop_na(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                       WaterCoverage, GRASSLAND_TYPE, Ave_WiderWater500)

## Select columns wanted for modelling
colnames(LapF)
LapF <- LapF |> 
           select(est_pairsL, S_LOC_ID, Landscape, FieldArea, CorvDens,  Reserve, SSSI, ESS_Wader, CSS_Wader,
                  TALL_BOUNDARY_PERCENT, WaterCoverage, Lap_Density, InterTidal_Distm, GRASSLAND_TYPE, F_LOC_ID, Fence_Coverage,
                  STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                  PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                  PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                  PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                  PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                  PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                  PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000) 

## Set each column as the correct variable type (factor or numeric)
LapF <- LapF |>
            mutate(across(c(est_pairsL, FieldArea, CorvDens, STANDING_WATER_TOTAL_PERCENT,  RUSH_PERCENT, 
                            WaterCoverage, Lap_Density, InterTidal_Distm, TALL_BOUNDARY_PERCENT, 
                            PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                            PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                            PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                            PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                            PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                            PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000), as.numeric),
                   across(c(Reserve, SSSI, ESS_Wader, CSS_Wader, STOCK, VEG_STRUCTURE, GROUND_DAMP, Landscape, GRASSLAND_TYPE, Fence_Coverage), as.factor)) %>% 
           as.data.frame()


## Summary of the number of data points used in the model
nrow(LapF) # total number of feilds in the model
length(unique(LapF$S_LOC_ID))
table(LapF$Landscape) # number of fields per region
LapF |> group_by(Landscape) |> summarise(Sites = length(unique(S_LOC_ID))) # number of sites 




##-------------------------------##
#### Lapwing: Train/Test split ####
##-------------------------------##

## Split the data set into a train-test set
## Choose x percent of different S_LOC_IDs to go into the train set
set.seed(104)
v <- unique(LapF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

## make the split
Lap_train <- LapF |> filter(S_LOC_ID %in% Sampv)
Lap_test <- LapF |> filter(!S_LOC_ID %in% Sampv)


##----------------------##
#### Lapwing: Tune RF ####
##----------------------##

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
# TuneParams <- tune(est_pairsL ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
#                                WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
#                                WaterCoverage*Landscape + 
#                                PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                                PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                                PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                                PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                                PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                                PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                                PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                                PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage, 
#                    data = Lap_train,
#                    mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
#                    nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)
# 
# ## plot the surface from parameter tuning
# TuneParams$optimal # returns optimal parameters from tuning
# plot.tune(TuneParams)



##---------------------##
#### Lapwing: Run RF ####
##---------------------##

## Run random forest regression
set.seed(112)
RF_fit <- rfsrc(est_pairsL ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
                               WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                               WaterCoverage*Landscape + 
                               PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
                               PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
                               PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
                               PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
                               PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
                               PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
                               PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
                               PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
                data = Lap_train, ntree = 2000, mtry = 5, nodesize = 1)
saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF.rds") # Save model for later use if wanted

## run full model with all the data
# set.seed(112)
# RF_fit <- rfsrc(est_pairsL ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
#                                WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
#                                WaterCoverage*Landscape + 
#                                PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                                PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                                PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                                PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                                PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                                PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                                PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                                PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
#                 data = LapF, ntree = 2000, mtry = 5, nodesize = 1)
# saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds") # Save model for later use if wanted

## model summary
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model
RF_fit 
sqrt(mean((Lap_train$est_pairsL - RF_fit$predicted.oob)^2)) # root mean square oob error


## plot variable importance
## VIMP (variable importance) is a technique for estimating the importance of a variable by 
## comparing performance of the estimated model with and without the variable in it.
o <- vimp(RF_fit)
## save plot as png
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_Vimp.png", width = 40, height = 20, units = "cm", res = 1200)
plot(o, ylab = "Vairb")
dev.off()


## Partial Dependence Plot
## These plots display the predicted conditional mean of the outcome as a function of an x variable
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_VarRel1.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("Fence_Coverage", "FieldArea", "Landscape"), partial = TRUE)
dev.off()

png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_VarRel2.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("WaterCoverage", "PropUrban_2000", "PropWetGrass_2000"), partial = TRUE)
dev.off()



# ## find possible interactions, not sure quite how this works
# find.interaction(RF_fit, method = "vimp", nrep = 3)
# 
# ## Create confidence intervals around variable importance (by subsampling) and plot them
# smp_o <- subsample(RF_fit, B=100, subratio= 0.5)
# plot.subsample(smp_o, alpha = 0.05)



##----------------------------------##
#### Lapwing: Predict on test set ####
##----------------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fit, Lap_test)
Lap_test$Predicted <- pred_vals$predicted

## Calculate the root mean square error
sqrt(mean((Lap_test$est_pairsL - Lap_test$Predicted)^2))

## plot the real vs predicted values
ggplot(data = Lap_test, mapping = aes(y= Predicted, x= est_pairsL)) + 
  geom_point(colour = "orange", alpha = 0.3, size = 2) +
  xlab("Real Lapwing Pairs") + ylab("Predicited Lapwing Pairs")+
  geom_smooth(formula = y~x) +
  theme_light() 
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_RealvsPred_TestSet.png", 
       plot = last_plot(), width = 15, height = 15, units = "cm")


## Calculate the predicted number of pairs for each region and all regions combined
## All regions combined first
sum(Lap_test$Predicted)
sum(Lap_test$est_pairsL)

Lap_test2 <- filter(Lap_test, Landscape == "Broads")
sum(Lap_test2$Predicted)
sum(Lap_test2$est_pairsL)
sqrt(mean((Lap_test2$est_pairsL - Lap_test2$Predicted)^2)) # RMSE

Lap_test3 <- filter(Lap_test, Landscape == "Greater Thames")
sum(Lap_test3$Predicted)
sum(Lap_test3$est_pairsL)
sqrt(mean((Lap_test3$est_pairsL - Lap_test3$Predicted)^2)) # RMSE

Lap_test4 <- filter(Lap_test, Landscape == "Somerset Levels and Moors")
sum(Lap_test4$Predicted)
sum(Lap_test4$est_pairsL)
sqrt(mean((Lap_test4$est_pairsL - Lap_test4$Predicted)^2)) # RMSE




## Add on the site total for each site
Lap_test <- Lap_test |>  group_by(S_LOC_ID) |> 
              summarise(SitePopPred = sum(Predicted, na.rm = T),
                        SitePopReal = sum(est_pairsL, na.rm = T)) |> 
              full_join(Lap_test, by = "S_LOC_ID")

## Add the field shapes onto the data set
L_TestFields <- left_join(Lap_test, BWWM, by = "F_LOC_ID") |> st_as_sf()

## filter out the landscape individually
LS_TestFields <- filter(L_TestFields, Landscape == "Somerset Levels and Moors")
LGT_TestFields <- filter(L_TestFields, Landscape == "Greater Thames")
LB_TestFields <- filter(L_TestFields, Landscape == "Broads")

## Calculate the RMSE at the site level
sqrt(mean((L_TestFields$SitePopReal - L_TestFields$SitePopPred)^2)) # RMSE
sqrt(mean((LS_TestFields$SitePopReal - LS_TestFields$SitePopPred)^2)) # RMSE
sqrt(mean((LGT_TestFields$SitePopReal - LGT_TestFields$SitePopPred)^2)) # RMSE
sqrt(mean((LB_TestFields$SitePopReal - LB_TestFields$SitePopPred)^2)) # RMSE

## Plot predicted vs real Lapwing abundance for Somerset
LapSFi1 <- ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Predicted") +
  theme_light()

LapSFi2 <- ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Real") +
  theme_light()

## Arrange plots side-by side and save
LapSFi <- grid.arrange(LapSFi1, LapSFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_Som_FieldComp.png", plot = LapSFi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Lapwing abundance for Broads
LapBFi1 <- ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Predicted") +
  theme_light()

LapBFi2 <-ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Real") +
  theme_light()

## Arrange plots side-by side and save
LapBFi <- grid.arrange(LapBFi1, LapBFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_Broads_FieldComp.png", plot = LapBFi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Lapwing abundance for Greater Thames
## Set the plot extent so that all plots have the same area no matter if they have different landarcels
PlotExt <- coord_sf(xlim = c(st_bbox(LGT_TestFields)[1]+22000, st_bbox(LGT_TestFields)[3]-4000), 
                    ylim = c(st_bbox(LGT_TestFields)[2]-500, st_bbox(LGT_TestFields)[4]-9000), 
                    crs = 27700, expand = FALSE) 
LapGTFi1 <- ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Predicted") +
  PlotExt +
  theme_light()

LapGTFi2 <- ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Real") +
  PlotExt +
  theme_light()

## Arrange plots side-by side and save
LapGTFi <- grid.arrange(LapGTFi1, LapGTFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_GT_FieldComp.png", plot = LapGTFi, width = 40, height = 20, units = "cm")




## Plot predicted vs real Lapwing abundance for Somerset
LapSSi1 <- ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Predicted") +
  theme_light()

LapSSi2 <- ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Real") +
  theme_light()

## Arrange plots side-by side and save
LapSSi <- grid.arrange(LapSSi1, LapSSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_Som_SiteComp.png", plot = LapSSi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Lapwing abundance for Broads
LapBSi1 <- ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Predicted") +
  theme_light()

LapBSi2 <- ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Real") +
  theme_light()

## Arrange plots side-by side and save
LapBSi <- grid.arrange(LapBSi1, LapBSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_Broads_SiteComp.png", plot = LapBSi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Lapwing abundance for Greater Thames
LapGTSi1 <- ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Predicted") +
  PlotExt +
  theme_light()

LapGTSi2 <- ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Real") +
  PlotExt +
  theme_light()

## Arrange plots side-by side and save
LapGTSi <- grid.arrange(LapGTSi1, LapGTSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Lapwing_GT_SiteComp.png", plot = LapGTSi, width = 40, height = 20, units = "cm")




##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##



##---------------------------------##
#### 2.0 Redshank One-step Model ####
##---------------------------------##

##-------------------------##
#### Redshank: Prep Data ####
##-------------------------##

## Filter out the Redshank survey data fro the Thames and Broads
Red <- filter(WData, is.na(est_pairsR)==F) |> filter(Landscape %in% c("Greater Thames", "Broads"))
nrow(Red); nrow(filter(Red, est_pairsR > 0)) # total no. of survey fields and total with waders 

## Check if survey columns have NA's
colnames(Red) # get column names
Red |> select(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, GROUND_DAMP, STOCK, VEG_HEIGHT,
              VEG_STRUCTURE, RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE, WaterCoverage, Fence_Coverage) |> 
       is.na() |> summary()

## remove rows with with NAs in these columns
RedF <- Red |> drop_na(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                       WaterCoverage, GRASSLAND_TYPE, Fence_Coverage, Ave_WiderWater500)


## Select columns wanted for modelling
colnames(RedF)
RedF <- RedF |> 
           select(est_pairsR, S_LOC_ID, Landscape, FieldArea, CorvDens,  Reserve, SSSI, ESS_Wader, CSS_Wader,
                  TALL_BOUNDARY_PERCENT, WaterCoverage, Red_Density, InterTidal_Distm, GRASSLAND_TYPE, F_LOC_ID, Fence_Coverage,
                  STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                  PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                  PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                  PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                  PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                  PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                  PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000) 

## Set each column as the correct variable type (factor or numeric)
RedF <- RedF |>
            mutate(across(c(est_pairsR, FieldArea, CorvDens, STANDING_WATER_TOTAL_PERCENT,  RUSH_PERCENT, 
                            WaterCoverage, Red_Density, InterTidal_Distm, TALL_BOUNDARY_PERCENT, 
                            PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                            PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                            PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                            PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                            PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                            PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000), as.numeric),
                   across(c(Reserve, SSSI, ESS_Wader, CSS_Wader, STOCK, VEG_STRUCTURE, GROUND_DAMP, Landscape, GRASSLAND_TYPE, Fence_Coverage), as.factor)) %>% 
           as.data.frame()


## Summary of the number of data points used in the model
nrow(RedF) # total number of feilds in the model
length(unique(RedF$S_LOC_ID))
table(RedF$Landscape) # number of fields per region
RedF |> group_by(Landscape) |> summarise(Sites = length(unique(S_LOC_ID))) # number of sites 




##--------------------------------##
#### Redshank: Train/Test split ####
##--------------------------------##

## Split the data set into a train-test set
## Choose x percent of different S_LOC_IDs
set.seed(121)
v <- unique(RedF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

Red_train <- RedF |> filter(S_LOC_ID %in% Sampv)
Red_test <- RedF |> filter(!S_LOC_ID %in% Sampv)



##----------------------##
#### Redshank: Tune RF ####
##----------------------##

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
# TuneParams <- tune(est_pairsR ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
#                                WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
#                                WaterCoverage*Landscape + InterTidal_Distm*Landscape +
#                                PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                                PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                                PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                                PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                                PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                                PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                                PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                                PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
#                    data = Red_train,
#                    mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
#                    nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)
# 
# ## plot the surface from parameter tuning
# TuneParams$optimal # returns optimal parameters from tuning
# plot.tune(TuneParams)



##---------------------##
#### Redshank: Run RF ####
##---------------------##

## Run random forest regression
set.seed(1212)
RF_fit <- rfsrc(est_pairsR ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
                               WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                               WaterCoverage*Landscape + InterTidal_Distm*Landscape +
                               PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
                               PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
                               PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
                               PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
                               PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
                               PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
                               PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
                               PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
                data = Red_train, ntree = 2000, mtry = 5, nodesize = 1)
saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF.rds") # Save model for later use if wanted

## run model on full data
# set.seed(1212)
# RF_fit <- rfsrc(est_pairsR ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage +
#                                WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT +
#                                WaterCoverage*Landscape + InterTidal_Distm*Landscape +
#                                PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                                PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                                PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                                PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                                PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                                PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                                PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                                PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
#                 data = RedF, ntree = 2000, mtry = 5, nodesize = 1)
# saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds") # Save model for later use if wanted

## model summary
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model
RF_fit
sqrt(mean((Red_train$est_pairsR - RF_fit$predicted.oob)^2)) # root mean square oob error


## plot variable importance
## VIMP (variable importance) is a technique for estimating the importance of a variable by 
## comparing performance of the estimated model with and without the variable in it.
o <- vimp(RF_fit)
## save plot as png
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_Vimp.png", width = 40, height = 20, units = "cm", res = 1200)
plot(o)
dev.off()


## Partial Dependence Plot
## These plots display the predicted conditional mean of the outcome as a function of an x variable
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_VarRel1.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("FieldArea","Fence_Coverage", "PropUrban_2000"), partial = TRUE)
dev.off()

png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_VarRel2.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("Ave_WiderWater500", "PropWetGrass_1000","PropWetGrass_1500"), partial = TRUE)
dev.off()



# ## find possible interactions, not sure quite how this works
# find.interaction(RF_fit, method = "vimp", nrep = 3)
# 
# ## Create confidence intervals around variable importance (by sub-sampling) and plot them
# smp_o <- subsample(RF_fit, B=100, subratio= 0.5)
# plot.subsample(smp_o, alpha = 0.05)



##----------------------------------##
#### Redshank: Predict on test set ####
##----------------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fit, Red_test)
Red_test$Predicted <- pred_vals$predicted

## Calculate the root mean square error
sqrt(mean((Red_test$est_pairsR - Red_test$Predicted)^2))

## plot the real vs predicted values
ggplot(data = Red_test, mapping = aes(y= Predicted, x= est_pairsR)) + 
  geom_point(colour = "orange", alpha = 0.3, size = 2) +
  xlab("Real Redshank Pairs") + ylab("Predicited Redshank Pairs")+
  geom_smooth(formula = y~x) +
  theme_light() 
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_RealvsPred_TestSet.png", 
       plot = last_plot(), width = 15, height = 15, units = "cm")
sum(Red_test$Predicted) # predicted pairs
sum(Red_test$est_pairsR) # observed pairs



## Calculate the predicted number of pairs for each region and all regions combined
## All regions combined first
sum(Red_test$Predicted)
sum(Red_test$est_pairsR)

Red_test2 <- filter(Red_test, Landscape == "Broads")
sum(Red_test2$Predicted)
sum(Red_test2$est_pairsR)
sqrt(mean((Red_test2$est_pairsR - Red_test2$Predicted)^2)) # RMSE

Red_test3 <- filter(Red_test, Landscape == "Greater Thames")
sum(Red_test3$Predicted)
sum(Red_test3$est_pairsR)
sqrt(mean((Red_test3$est_pairsR - Red_test3$Predicted)^2)) # RMSE



## Add on the site total for each site
Red_test <- Red_test |>  group_by(S_LOC_ID) |> 
              summarise(SitePopPred = sum(Predicted, na.rm = T),
                        SitePopReal = sum(est_pairsR, na.rm = T)) |> 
              full_join(Red_test, by = "S_LOC_ID")

## Add the field shapes onto the data set
R_TestFields <- left_join(Red_test, BWWM, by = "F_LOC_ID") |> st_as_sf()

## filter out the landscape individually
RGT_TestFields <- filter(R_TestFields, Landscape == "Greater Thames")
RB_TestFields <- filter(R_TestFields, Landscape == "Broads")

## Calculate the RMSE at the site level
sqrt(mean((R_TestFields$SitePopReal - R_TestFields$SitePopPred)^2)) # RMSE
sqrt(mean((RGT_TestFields$SitePopReal - RGT_TestFields$SitePopPred)^2)) # RMSE
sqrt(mean((RB_TestFields$SitePopReal - RB_TestFields$SitePopPred)^2)) # RMSE


## Plot predicted vs real Redshank abundance for Broads
RedBFi1 <- ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Predicted") +
  theme_light()

RedBFi2 <-ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsR)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Real") +
  theme_light()

## Arrange plots side-by side and save
RedBFi <- grid.arrange(RedBFi1, RedBFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_Broads_FieldComp.png", plot = RedBFi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Redshank abundance for Greater Thames
## Set the plot extent so that all plots have the same area no matter if they have different landarcels
PlotExt <- coord_sf(xlim = c(st_bbox(RGT_TestFields)[1]+22000, st_bbox(RGT_TestFields)[3]-4000), 
                    ylim = c(st_bbox(RGT_TestFields)[2]-500, st_bbox(RGT_TestFields)[4]-9000), 
                    crs = 27700, expand = FALSE) 
RedGTFi1 <- ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Predicted") +
  PlotExt +
  theme_light()

RedGTFi2 <- ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsR)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Real") +
  PlotExt +
  theme_light()

## Arrange plots side-by side and save
RedGTFi <- grid.arrange(RedGTFi1, RedGTFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_GT_FieldComp.png", plot = RedGTFi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Redshank abundance for Broads
RedBSi1 <- ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Predicted") +
  theme_light()

RedBSi2 <- ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Real") +
  theme_light()

## Arrange plots side-by side and save
RedBSi <- grid.arrange(RedBSi1, RedBSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_Broads_SiteComp.png", plot = RedBSi, width = 40, height = 20, units = "cm")



## Plot predicted vs real Redshank abundance for Greater Thames
RedGTSi1 <- ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Predicted") +
  PlotExt +
  theme_light()

RedGTSi2 <- ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Real") +
  PlotExt +
  theme_light()

## Arrange plots side-by side and save
RedGTSi <- grid.arrange(RedGTSi1, RedGTSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_GT_SiteComp.png", plot = RedGTSi, width = 40, height = 20, units = "cm")





##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##



##------------------------------##
#### 3.0 Snipe One-step Model ####
##------------------------------##

##----------------------##
#### Snipe: Prep Data ####
##----------------------##

## Filter out the Snipe survey data
Snipe <- filter(WData, is.na(est_pairsS)==F)
Snipe <- filter(Snipe, Landscape == "Somerset Levels and Moors")
nrow(Snipe); nrow(filter(Snipe, est_pairsS > 0)) # total no. of survey fields and total with waders 

## Check if survey columns have NA's
colnames(Snipe) # get column names
Snipe |> select(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, GROUND_DAMP, STOCK, VEG_HEIGHT,
              VEG_STRUCTURE, RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE, WaterCoverage, Fence_Coverage) |> 
       is.na() |> summary()

## remove rows with with NAs in these columns
SnipeF <- Snipe |> drop_na(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                       WaterCoverage, GRASSLAND_TYPE, Fence_Coverage, Ave_WiderWater500)

## Select columns wanted for modelling
colnames(SnipeF)
SnipeF <- SnipeF |> 
           select(est_pairsS, S_LOC_ID, Landscape, FieldArea, CorvDens,  Reserve, SSSI, ESS_Wader, CSS_Wader,
                  TALL_BOUNDARY_PERCENT, WaterCoverage, Red_Density, InterTidal_Distm, GRASSLAND_TYPE, F_LOC_ID, Fence_Coverage,
                  STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, Peat_Soil,
                  PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                  PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                  PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                  PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                  PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                  PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000) 

## Set each column as the correct variable type (factor or numeric)
SnipeF <- SnipeF |>
            mutate(across(c(est_pairsS, FieldArea, CorvDens, STANDING_WATER_TOTAL_PERCENT,  RUSH_PERCENT, 
                            WaterCoverage, Red_Density, InterTidal_Distm, TALL_BOUNDARY_PERCENT, 
                            PropWood_500, PropUrban_500, PropWetGrass_500, Ave_WiderWater500,
                            PropWood_750, PropUrban_750, PropWetGrass_750, Ave_WiderWater750,
                            PropWood_1000, PropUrban_1000, PropWetGrass_1000, Ave_WiderWater1000,
                            PropWood_1250, PropUrban_1250, PropWetGrass_1250, Ave_WiderWater1250,
                            PropWood_1500, PropUrban_1500, PropWetGrass_1500, Ave_WiderWater1500,
                            PropWood_2000, PropUrban_2000, PropWetGrass_2000, Ave_WiderWater2000), as.numeric),
                   across(c(Reserve, SSSI, ESS_Wader, CSS_Wader, STOCK, VEG_STRUCTURE, GROUND_DAMP, Landscape, GRASSLAND_TYPE, Fence_Coverage, Peat_Soil), as.factor)) %>% 
           as.data.frame()

SnipeF <- SnipeF |> mutate(est_pairsSW = as.factor(round2(est_pairsS, digits =0)))


## Summary of the number of data points used in the model
nrow(SnipeF) # total number of feilds in the model
length(unique(SnipeF$S_LOC_ID))
table(SnipeF$Landscape) # number of fields per region
SnipeF |> group_by(Landscape) |> summarise(Sites = length(unique(S_LOC_ID))) # number of sites 




##-----------------------------##
#### Snipe: Train/Test split ####
##-----------------------------##

## Split the data set into a train-test set
## Choose x percent of different S_LOC_IDs
set.seed(111)
v <- unique(SnipeF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

Sni_train <- SnipeF |> filter(S_LOC_ID %in% Sampv)
Sni_test <- SnipeF |> filter(!S_LOC_ID %in% Sampv)



##--------------------##
#### Snipe: Tune RF ####
##--------------------##

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
# TuneParams <- tune(est_pairsS  ~ FieldArea + CorvDens + GRASSLAND_TYPE + Fence_Coverage + Peat_Soil +
#                              WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
#                              WaterCoverage*Peat_Soil +
#                              PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                              PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                              PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                              PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                              PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                              PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                              PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                              PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage,
#                    data = Sni_train,
#                    mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
#                    nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)
# 
# ## plot the surface from parameter tuning
# TuneParams$optimal # returns optimal parameters from tuning
# plot.tune(TuneParams)



##-------------------##
#### Snipe: Run RF ####
##-------------------##


## Run random forest regression
set.seed(1212)
RF_fit <- rfsrc(est_pairsS ~ FieldArea + CorvDens + GRASSLAND_TYPE + Fence_Coverage + Peat_Soil +
                             WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                             WaterCoverage*Peat_Soil +
                             PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
                             PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
                             PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
                             PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
                             PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
                             PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
                             PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
                             PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage, 
                data = Sni_train, 
                ntree = 2000, mtry = 4, nodesize = 1)
saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/SnipeRF.rds") # Save model for later use if wanted

## run on full data set
# set.seed(1212)
# RF_fit <- rfsrc(est_pairsS ~ FieldArea + CorvDens + GRASSLAND_TYPE + Fence_Coverage + Peat_Soil +
#                              WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
#                              WaterCoverage*Peat_Soil +
#                              PropWood_500 + PropUrban_500 + PropWetGrass_500 + Ave_WiderWater500 +
#                              PropWood_1000 + PropUrban_1000 + PropWetGrass_1000 + Ave_WiderWater1000 +
#                              PropWood_1500 + PropUrban_1500 + PropWetGrass_1500 + Ave_WiderWater1500 +
#                              PropWood_2000 + PropUrban_2000 + PropWetGrass_2000 + Ave_WiderWater2000 +
#                              PropWood_500*WaterCoverage + PropUrban_500*WaterCoverage + PropWetGrass_500*WaterCoverage + Ave_WiderWater500*WaterCoverage +
#                              PropWood_1000*WaterCoverage + PropUrban_1000*WaterCoverage + PropWetGrass_1000*WaterCoverage + Ave_WiderWater1000*WaterCoverage +
#                              PropWood_1500*WaterCoverage + PropUrban_1500*WaterCoverage + PropWetGrass_1500*WaterCoverage + Ave_WiderWater1500*WaterCoverage +
#                              PropWood_2000*WaterCoverage + PropUrban_2000*WaterCoverage + PropWetGrass_2000*WaterCoverage + Ave_WiderWater2000*WaterCoverage, 
#                 data = SnipeF, 
#                 ntree = 2000, mtry = 4, nodesize = 1)
# saveRDS(RF_fit, "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/SnipeRF_FullData.rds") # Save model for later use if wanted


RF_fit # model summary
sqrt(mean((Sni_train$est_pairsS - RF_fit$predicted.oob)^2)) # root mean square oob error
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model

## plot variable importance
## VIMP (variable importance) is a technique for estimating the importance of a variable by 
## comparing performance of the estimated model with and without the variable in it.
o <- vimp(RF_fit)
## save plot as png
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Snipe_Vimp.png", width = 40, height = 20, units = "cm", res = 1200)
plot(o)
dev.off()



## Partial Dependence Plot
## These plots display the predicted conditional mean of the outcome as a function of an x variable
png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Snipe_VarRel1.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("Ave_WiderWater1000", "Ave_WiderWater1500", "Ave_WiderWater2000"), partial = TRUE)
dev.off()

png(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Snipe_VarRel2.png", width = 20, height = 10, units = "cm", res = 1200)
plot.variable(RF_fit, c("Ave_WiderWater500", "PropWetGrass_500", "PropWood_2000"), partial = TRUE)
dev.off()


# ## Create confidence intervals around variable importance (by subsampling) and plot them
# smp_o <- subsample(RF_fit, B=100, subratio= 0.5)
# plot.subsample(smp_o, alpha = 0.05)
# 
# ## find possible interactions, not sure quite how this works
# find.interaction(RF_fit, method = "vimp", nrep = 3)




##--------------------------------##
#### Snipe: Predict on test set ####
##--------------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fit, Sni_test)
Sni_test$Predicted <- pred_vals$predicted


## Calculate the root mean square error
sqrt(mean((Sni_test$est_pairsS - Sni_test$Predicted)^2))

## plot the real vs predicted values
ggplot(data = Sni_test, mapping = aes(y= Predicted, x= est_pairsS)) + 
  geom_point(colour = "orange", alpha = 0.5, size = 3) +
  xlab("Real Snipe Density") + ylab("Predicited Snipe Density")+
  geom_smooth(formula = y~x) +
  theme_light() 
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Redshank_RealvsPred_TestSet.png", 
       plot = last_plot(), width = 15, height = 15, units = "cm")
sum(as.numeric(as.character(Sni_test$Predicted))) # predicted pairs
sum(Sni_test$est_pairsS) # observed pairs


## Add on the site total for each site
Sni_test <- Sni_test |>  group_by(S_LOC_ID) |> 
              summarise(SitePopPred = sum(Predicted, na.rm = T),
                        SitePopReal = sum(est_pairsS, na.rm = T)) |> 
              full_join(Sni_test, by = "S_LOC_ID")

## Add the field shapes onto the data set
S_TestFields <- left_join(Sni_test, BWWM, by = "F_LOC_ID") |> st_as_sf()

## Calculate the RMSE at the site level
sqrt(mean((S_TestFields$SitePopReal - S_TestFields$SitePopPred)^2)) # RMSE


## Plot predicted vs real Snipe abundance for Somerset
SniSFi1 <- ggplot() + 
  geom_sf(data = S_TestFields, mapping=aes(geometry=geometry, fill = (Predicted)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Snipe Somerset: Predicted") +
  theme_light()

SniSFi2 <- ggplot() + 
  geom_sf(data = S_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsS)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Snipe Somerset: Real") +
  theme_light()

## Arrange plots side-by side and save
SniSFi <- grid.arrange(SniSFi1, SniSFi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Snipe_Som_FieldComp.png", plot = SniSFi, width = 40, height = 20, units = "cm")




## Plot predicted vs real Snipe abundance for Somerset
SniSSi1 <- ggplot() + 
  geom_sf(data = S_TestFields, mapping=aes(geometry=geometry, fill = (SitePopPred)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Snipe Somerset: Predicted") +
  theme_light()

SniSSi2 <- ggplot() + 
  geom_sf(data = S_TestFields, mapping=aes(geometry=geometry, fill = (SitePopReal)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Snipe Somerset: Real") +
  theme_light()

## Arrange plots side-by side and save
SniSSi <- grid.arrange(SniSSi1, SniSSi2, ncol=2)
ggsave(filename = "CleanData/Wader Abundance/5-ModelWaderAbundance/Model Plots/Snipe_Som_SiteComp.png", plot = SniSSi, width = 40, height = 20, units = "cm")




##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##




##--------------------------------##
#### 4.0 Lapwing Two-step Model ####
##--------------------------------##

##--------------------------------##
#### Lapwing2: Train/Test split ####
##--------------------------------##

# ## Split the data set into a train-test set, randomly choose rows
# set.seed(1012)
# inTrain <- createDataPartition(
#   y = LapF$est_pairsL, # the outcome data
#   p = 0.7, # The percentage of data in the training set
#   list = FALSE)
# 
# ## Separate the train and test data set
# LapF <- LapF[ inTrain,]
# Lap_test  <- LapF[-inTrain,]

## Split the data set into a train-test set
## Choose x percent of different S_LOC_IDs
set.seed(3412)
v <- unique(LapF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

## Add column that indicates whether Lapwing where 
LapF <- LapF |> mutate(LapPres = ifelse(est_pairsL > 0, 1, 0),
                       LapPres = as.factor(LapPres))

## split the data into a training/test split
Lap_train <- LapF |> filter(S_LOC_ID %in% Sampv)
Lap_test <- LapF |> filter(!S_LOC_ID %in% Sampv)



##-----------------------##
#### Lapwing2: Tune RF ####
##-----------------------##

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
TuneParams <- tune(Lap_Density ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage +
                             WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + GROUND_DAMP + RUSH_PERCENT + 
                             PropWood_500 + PropWetGrass_500 + PropUrban_500 +
                             WaterCoverage*Landscape +
                             PropWood_500*WaterCoverage + PropWetGrass_500*WaterCoverage + PropUrban_500*WaterCoverage +
                             PropWood_300*WaterCoverage + PropUrban_300*WaterCoverage +  PropWetGrass_300*WaterCoverage +
                             PropWood_900*WaterCoverage + PropUrban_900*WaterCoverage +  PropWetGrass_900*WaterCoverage + 
                             PropWood_1300*WaterCoverage + PropUrban_1300*WaterCoverage +  PropWetGrass_1300*WaterCoverage + 
                             PropWood_1700*WaterCoverage + PropUrban_1700*WaterCoverage +  PropWetGrass_1700*WaterCoverage,
                   data = Lap_train,
                   mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
                   nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)
## returns optimal parameters from tuning
TuneParams$optimal 
## plot the surface from parameter tuning
plot.tune(TuneParams)



##----------------------##
#### Lapwing2: Run RF ####
##----------------------##

## Run random forest regression to predict Lapwing presence/absence
RF_fitBin <- imbalanced(LapPres ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
                               WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                               PropWood_500 + PropWetGrass_500 + PropUrban_500 + 
                               WaterCoverage*Landscape +
                               PropWood_500*WaterCoverage + PropWetGrass_500*WaterCoverage + PropUrban_500*WaterCoverage +
                               PropWood_300*WaterCoverage + PropUrban_300*WaterCoverage +  PropWetGrass_300*WaterCoverage +
                               PropWood_900*WaterCoverage + PropUrban_900*WaterCoverage +  PropWetGrass_900*WaterCoverage + 
                               PropWood_1300*WaterCoverage + PropUrban_1300*WaterCoverage +  PropWetGrass_1300*WaterCoverage + 
                               PropWood_1700*WaterCoverage + PropUrban_1700*WaterCoverage +  PropWetGrass_1700*WaterCoverage,
                data = Lap_train, ntree = 5000, method = "rfq", splitrule = "auc")
get.imbalanced.performance(RF_fitBin)
RF_fitBin

## plot variable importance
o <- vimp(RF_fitBin)
plot(o)

## Run random forest regression to predict Lapwing Density when they are present
## First filter out data from fields with Lapwing
LaptrainPres <- filter(Lap_train, Lap_Density > 0)
RF_fitCont <- rfsrc(Lap_Density ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage + 
                               WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                               PropWood_500 + PropWetGrass_500 + PropUrban_500 + 
                               WaterCoverage*Landscape +
                               PropWood_500*WaterCoverage + PropWetGrass_500*WaterCoverage + PropUrban_500*WaterCoverage +
                               PropWood_300*WaterCoverage + PropUrban_300*WaterCoverage +  PropWetGrass_300*WaterCoverage +
                               PropWood_900*WaterCoverage + PropUrban_900*WaterCoverage +  PropWetGrass_900*WaterCoverage + 
                               PropWood_1300*WaterCoverage + PropUrban_1300*WaterCoverage +  PropWetGrass_1300*WaterCoverage + 
                               PropWood_1700*WaterCoverage + PropUrban_1700*WaterCoverage +  PropWetGrass_1700*WaterCoverage,
                data = LaptrainPres, ntree = 3000, mtry = 4, nodesize = 2)

RF_fitCont # model summary
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model

## plot variable importance
## VIMP (variable importance) is a technique for estimating the importance of a variable by 
## comparing performance of the estimated model with and without the variable in it.
o <- vimp(RF_fitCont)
plot(o)

## plot how the error rate change with the number of trees
plot(RF_fitCont)
find.interaction(RF_fitCont, method = "vimp", nrep = 3)

## Create confidence intervals around variable importance (by subsampling) and plot them
smp_o <- subsample(RF_fitCont, B=100, subratio= 0.5)
plot.subsample(smp_o, alpha = 0.05)


## Partial Dependence Plot
## These plots display the predicted conditional mean of the outcome as a function of an x variable
plot.variable(RF_fitCont, "Fence_Coverage", partial = TRUE)  
plot.variable(RF_fitCont, "WaterCoverage", partial = TRUE)    
plot.variable(RF_fitCont, "CorvDens", partial = TRUE)
plot.variable(RF_fitCont, "Landscape", partial = TRUE)
plot.variable(RF_fitCont, "FieldArea", partial = TRUE)
plot.variable(RF_fitCont, "PropWetGrass_500", partial = TRUE)
plot.variable(RF_fitCont, "PropUrban_500", partial = TRUE)
plot.variable(RF_fitCont, "PropWood_500", partial = TRUE)
plot.variable(RF_fitCont, "TALL_BOUNDARY_PERCENT", partial = TRUE)
plot.variable(RF_fitCont, "VEG_STRUCTURE", partial = TRUE)
plot.variable(RF_fitCont, "RUSH_PERCENT", partial = TRUE)

## find possible interactions, not sure quite how this works
find.interaction(RF_fitCont, method = "vimp", nrep = 3)



##----------------------------------##
#### Lapwing2: Predict on test set ####
##----------------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fitBin, Lap_test)
Lap_test$Predicted <- pred_vals$class
Lap_test$Prob_0 <- pred_vals[["predicted"]][,1]
Lap_test$Prob_1 <- pred_vals[["predicted"]][,2]

## Create confusion matrix, when rounding real and predicted values
example <- confusionMatrix(data= pred_vals$class, 
                           reference = Lap_test$LapPres)
example


## Filter the test set for only the fields that have been classified as 1
Lap_test_Pres <- filter(Lap_test, Predicted == 1)

## Then predict the lapwing density for the fields predicted to have Lapwing
pred_vals2 <- predict(object = RF_fitCont, Lap_test_Pres)

## Assign the predicted abundance to a new column, then multiple the probability of presence and field area
Lap_test_Pres$PredictedAbund <- pred_vals2$predicted
Lap_test_Pres$PredictedAbund <- Lap_test_Pres$PredictedAbund*Lap_test_Pres$FieldArea

## Filter out fields the when surveyed had Lapwing or did not
Absent <- Lap_test_Pres |> filter(LapPres == 0)
Present <- Lap_test_Pres |> filter(LapPres == 1)
## Now calculate the number of pairs observed vs predicted
sum(Absent$PredictedAbund);sum(Present$PredictedAbund)
sum(Lap_test_Pres$PredictedAbund)
sum(Lap_test$est_pairsL)



## Filter the test set for only the fields that have been classified as 1
Lap_test_Abs <- filter(Lap_test, Predicted == 0) |> mutate(PredictedAbund = 0)
Lap_test_All <- rbind(Lap_test_Pres, Lap_test_Abs)

## For each region compare the real number of wader pairs to the predicted number of pairs
Lap_test2 <- filter(Lap_test_All, Landscape == "Broads")
sum(Lap_test2$PredictedAbund)
sum(Lap_test2$est_pairsL)

Lap_test3 <- filter(Lap_test_All, Landscape == "Greater Thames")
sum(Lap_test3$PredictedAbund)
sum(Lap_test3$est_pairsL)

Lap_test4 <- filter(Lap_test_All, Landscape == "Somerset Levels and Moors")
sum(Lap_test4$PredictedAbund)
sum(Lap_test4$est_pairsL)



## Add the field shapes onto the data set
L_TestFields <- left_join(Lap_test_All, BWWM, by = "F_LOC_ID") |> st_as_sf()

## filter out the landscape individually
LS_TestFields <- filter(L_TestFields, Landscape == "Somerset Levels and Moors")
LGT_TestFields <- filter(L_TestFields, Landscape == "Greater Thames")
LB_TestFields <- filter(L_TestFields, Landscape == "Broads")

## Plot predicted vs real Lapwing abundance for Somerset
ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (PredictedAbund)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = LS_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Somerset: Real") +
  theme_light()


## Plot predicted vs real Lapwing abundance for Broads
ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (PredictedAbund)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = LB_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Broads: Real") +
  theme_light()


## Plot predicted vs real Lapwing abundance for Greater Thames
ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (PredictedAbund)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = LGT_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsL)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Lapwing Greater Thames: Real") +
  theme_light()



##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##



##---------------------------------##
#### 5.0 Redshank Two-step Model ####
##---------------------------------##

##-------------------------##
#### Redshank2: Prep Data ####
##-------------------------##

## Filter out the Redshank survey data
Red <- filter(WData, is.na(est_pairsR)==F)
nrow(Red); nrow(filter(Red, est_pairsR > 0)) # total no. of survey fields and total with waders 

## Check if survey columns have NA's
colnames(Red) # get column names
Red |> select(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, GROUND_DAMP, STOCK, VEG_HEIGHT,
              VEG_STRUCTURE, RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE, WaterCoverage, Fence_Coverage) |> 
       is.na() |> summary()

## remove rows with with NAs in these columns
RedF <- Red |> drop_na(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                       WaterCoverage, GRASSLAND_TYPE, Fence_Coverage, WaterEdge)

## Select columns wanted for modelling
colnames(RedF)
RedF <- RedF |> 
           select(est_pairsR, S_LOC_ID, Landscape, FieldArea, CorvDens,  Reserve, SSSI, ESS_Wader, CSS_Wader, 
                  PropWood_500, PropUrban_500, PropIntense_500, PropWetGrass_500, PropSSSI_500, PropReserve_500, 
                  PropStewWa_500, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                  TALL_BOUNDARY_PERCENT, WaterCoverage, Red_Density, InterTidal_Distm, GRASSLAND_TYPE, F_LOC_ID, Fence_Coverage,
                  PropWood_900, PropUrban_900,  PropWetGrass_900, PropWood_1300, PropUrban_1300,  PropWetGrass_1300, 
                  PropWood_1700, PropUrban_1700,  PropWetGrass_1700, PropWood_300, PropUrban_300,  PropWetGrass_300,
                  Ave_WiderWater500, Ave_WiderWater1000, Ave_WiderWater1500, WaterEdge) 

## Set each column as the correct variable type (factor or numeric)
RedF <- RedF |>
            mutate(across(c(est_pairsR, FieldArea, CorvDens, PropWood_500, PropUrban_500, PropIntense_500, PropWetGrass_500, PropSSSI_500, 
                            PropReserve_500, PropStewWa_500, STANDING_WATER_TOTAL_PERCENT,  RUSH_PERCENT, TALL_BOUNDARY_PERCENT, WaterEdge,
                            PropWood_900, PropUrban_900,  PropWetGrass_900, PropWood_1300, PropUrban_1300,  PropWetGrass_1300, 
                            PropWood_1700, PropUrban_1700,  PropWetGrass_1700, PropWood_300, PropUrban_300,  PropWetGrass_300,
                            WaterCoverage, Red_Density, InterTidal_Distm, Ave_WiderWater500, Ave_WiderWater1000, Ave_WiderWater1500), as.numeric),
                   across(c(Reserve, SSSI, ESS_Wader, CSS_Wader, STOCK, VEG_STRUCTURE, GROUND_DAMP, Landscape, GRASSLAND_TYPE, Fence_Coverage), as.factor)) %>% 
           as.data.frame()


##-------------------------------##
#### Redshank2: Train/Test split ####
##-------------------------------##

# ## Split the data set into a train-test set, randomly choose rows
# set.seed(1012)
# inTrain <- createDataPartition(
#   y = RedF$est_pairsL, # the outcome data
#   p = 0.7, # The percentage of data in the training set
#   list = FALSE)
# 
# ## Separate the train and test data set
# RedF <- RedF[ inTrain,]
# Red_test  <- RedF[-inTrain,]

## Split the data set into a train-test set
## Choose x percent of different S_LOC_IDs
set.seed(1012)
v <- unique(RedF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

Red_train <- RedF |> filter(S_LOC_ID %in% Sampv)
Red_test <- RedF |> filter(!S_LOC_ID %in% Sampv)



##----------------------##
#### Redshank2: Tune RF ####
##----------------------##

## Management variables
## FieldArea + Reserve + SSSI + ESS_Wader + CSS_Wader + PropSSSI_500 + PropReserve_500 + PropStewWa_500
## Habitat variables
## FieldArea + CorvDens + STANDING_WATER_TOTAL_PERCENT + STOCK + VEG_STRUCTURE + GROUND_DAMP + RUSH_PERCENT + PropWood_500 + PropWetGrass_500 + PropUrban_500

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
TuneParams <- tune(Red_Density ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE +
                             WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                             PropWood_500 + PropWetGrass_500 + PropUrban_500 + InterTidal_Distm +
                             PropWood_300 + PropUrban_300 +  PropWetGrass_300 +
                             PropWood_900 + PropUrban_900 +  PropWetGrass_900 + 
                             PropWood_1300 + PropUrban_1300 +  PropWetGrass_1300 + 
                             PropWood_1700 + PropUrban_1700 +  PropWetGrass_1700 + 
                             WaterCoverage*Landscape + InterTidal_Distm*Landscape +
                             PropWood_500*WaterCoverage + PropWetGrass_500*WaterCoverage + PropUrban_500*WaterCoverage +
                             PropWood_300**WaterCoverage + PropUrban_300**WaterCoverage +  PropWetGrass_300**WaterCoverage +
                             PropWood_900**WaterCoverage + PropUrban_900**WaterCoverage +  PropWetGrass_900**WaterCoverage + 
                             PropWood_1300**WaterCoverage + PropUrban_1300**WaterCoverage +  PropWetGrass_1300**WaterCoverage + 
                             PropWood_1700**WaterCoverage + PropUrban_1700**WaterCoverage +  PropWetGrass_1700**WaterCoverage +
                             Ave_WiderWater500 + Ave_WiderWater1000 + Ave_WiderWater1500 +
                             Ave_WiderWater500*WaterCoverage + Ave_WiderWater1000*WaterCoverage + Ave_WiderWater1500*WaterCoverage,
                   data = Red_train,
                   mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
                   nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)

TuneParams$optimal # returns optimal parameters from tuning

## plot the surface from parameter tuning
plot.tune(TuneParams)



##---------------------##
#### Redshank2: Run RF ####
##---------------------##

# Red_trainSom <- filter(Red_train, Landscape == "Somerset Levels and Moors")
# Red_trainThames <- filter(Red_train, Landscape == "Greater Thames")
# Red_trainBr <- filter(Red_train, Landscape == "Broads")
# ggplot() + geom_point(data = Red_trainThames, mapping =aes(x =WaterCoverage, y=STANDING_WATER_TOTAL_PERCENT))
# ggplot() + geom_point(data = Red_trainSom, mapping =aes(x =WaterCoverage, y=STANDING_WATER_TOTAL_PERCENT))
# ggplot() + geom_point(data = Red_trainBr, mapping =aes(x =WaterCoverage, y=STANDING_WATER_TOTAL_PERCENT))

#WaterCoverage
#STANDING_WATER_TOTAL_PERCENT
## Run random forest regression
set.seed(1212)
RF_fit <- rfsrc(Red_Density ~ FieldArea + Landscape + CorvDens + GRASSLAND_TYPE + Fence_Coverage +
                             WaterCoverage + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + RUSH_PERCENT + 
                             PropWood_500 + PropWetGrass_500 + PropUrban_500 + InterTidal_Distm +
                             PropWood_300 + PropUrban_300 +  PropWetGrass_300 +
                             PropWood_900 + PropUrban_900 +  PropWetGrass_900 + 
                             PropWood_1300 + PropUrban_1300 +  PropWetGrass_1300 + 
                             PropWood_1700 + PropUrban_1700 +  PropWetGrass_1700 + 
                             WaterCoverage*Landscape + InterTidal_Distm*Landscape + WaterEdge +
                             PropWood_500*WaterCoverage + PropWetGrass_500*WaterCoverage + PropUrban_500*WaterCoverage +
                             PropWood_300*WaterCoverage + PropUrban_300*WaterCoverage +  PropWetGrass_300*WaterCoverage +
                             PropWood_900*WaterCoverage + PropUrban_900*WaterCoverage +  PropWetGrass_900*WaterCoverage + 
                             PropWood_1300*WaterCoverage + PropUrban_1300*WaterCoverage +  PropWetGrass_1300*WaterCoverage + 
                             PropWood_1700*WaterCoverage + PropUrban_1700*WaterCoverage +  PropWetGrass_1700*WaterCoverage +
                             Ave_WiderWater500 + Ave_WiderWater1000 + Ave_WiderWater1500 +
                             Ave_WiderWater500*WaterCoverage + Ave_WiderWater1000*WaterCoverage + Ave_WiderWater1500*WaterCoverage,
                data = Red_train, ntree = 3000, mtry = 3, nodesize = 1)


RF_fit # model summary
sqrt(mean((Red_train$Red_Density - RF_fit$predicted.oob)^2)) # root mean square oob error
o <- vimp(RF_fit)
plot(o)
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model

## plot variable importance
## VIMP (variable importance) is a technique for estimating the importance of a variable by 
## comparing performance of the estimated model with and without the variable in it.
o <- vimp(RF_fit)
plot(o)

## Create confidence intervals around variable importance (by subsampling) and plot them
smp_o <- subsample(RF_fit, B=100, subratio= 0.5)
plot.subsample(smp_o, alpha = 0.05)


## Partial Dependence Plot
## These plots display the predicted conditional mean of the outcome as a function of an x variable
plot.variable(RF_fit, "WaterCoverage", partial = TRUE)    
plot.variable(RF_fit, "CorvDens", partial = TRUE)
plot.variable(RF_fit, "Landscape", partial = TRUE)
plot.variable(RF_fit, "FieldArea", partial = TRUE)
plot.variable(RF_fit, "PropWetGrass_500", partial = TRUE)
plot.variable(RF_fit, "PropUrban_500", partial = TRUE)
plot.variable(RF_fit, "PropWood_500", partial = TRUE)
plot.variable(RF_fit, "TALL_BOUNDARY_PERCENT", partial = TRUE)
plot.variable(RF_fit, "VEG_STRUCTURE", partial = TRUE)
plot.variable(RF_fit, "RUSH_PERCENT", partial = TRUE)


## find possible interactions, not sure quite how this works
find.interaction(RF_fit, method = "vimp", nrep = 3)



##----------------------------------##
#### Redshank2: Predict on test set ####
##----------------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fit, Red_test)
Red_test$Predicted <- pred_vals$predicted

## Calculate the root mean square error
sqrt(mean((Red_test$Red_Density - Red_test$Predicted)^2))

## plot the real vs predicted values
ggplot(data = Red_test, mapping = aes(y= Predicted, x= Red_Density)) + 
  geom_point(colour = "orange", alpha = 0.5, size = 3) +
  xlab("Real Redshank Density") + ylab("Predicited Redshank Density")+
  geom_smooth(formula = y~x) +
  theme_light() 
sum(Red_test$Predicted*Red_test$FieldArea) # predicted pairs
sum(Red_test$est_pairsR) # observed pairs


Red_test2 <- filter(Red_test, Landscape == "Broads")
sum(Red_test2$Predicted*Red_test2$FieldArea)
sum(Red_test2$est_pairsR)

Red_test3 <- filter(Red_test, Landscape == "Greater Thames")
sum(Red_test3$Predicted*Red_test3$FieldArea)
sum(Red_test3$est_pairsR)

Red_test4 <- filter(Red_test, Landscape == "Somerset Levels and Moors")
sum(Red_test4$Predicted*Red_test4$FieldArea)
sum(Red_test4$est_pairsR)

sum(Red_test$Predicted*Red_test$FieldArea)
sum(Red_test$est_pairsR)


## Add the field shapes onto the data set
R_TestFields <- left_join(Red_test, BWWM, by = "F_LOC_ID") |> st_as_sf()

## filter out the landscape individually
RS_TestFields <- filter(R_TestFields, Landscape == "Somerset Levels and Moors")
RGT_TestFields <- filter(R_TestFields, Landscape == "Greater Thames")
RB_TestFields <- filter(R_TestFields, Landscape == "Broads")

## Plot predicted vs real redshank abundance for Somerset
ggplot() + 
  geom_sf(data = RS_TestFields, mapping=aes(geometry=geometry, fill = (Predicted*FieldArea)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Somerset: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = RS_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsR)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Somerset: Real") +
  theme_light()


## Plot predicted vs real redshank abundance for Broads
ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (Predicted*FieldArea)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = RB_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsR)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Broads: Real") +
  theme_light()


## Plot predicted vs real redshank abundance for Greater Thames
ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (Predicted*FieldArea)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Predicted") +
  theme_light()

ggplot() + 
  geom_sf(data = RGT_TestFields, mapping=aes(geometry=geometry, fill = (est_pairsR)), colour = NA) +
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ggtitle("Redshank Greater Thames: Real") +
  theme_light()



##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------------------------------------##

