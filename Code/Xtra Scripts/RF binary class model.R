## Filter out the Lapwing survey data
Lap <- filter(WData, is.na(est_pairsL)==F)
#Lap <- filter(Lap, Landscape == "Broads")
nrow(Lap); nrow(filter(Lap, est_pairsL > 0)) # total no. of survey fields and total with waders 

## Check if survey columns have NA's
colnames(Lap) # get column names
Lap |> select(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, GROUND_DAMP, STOCK, VEG_HEIGHT,
              VEG_STRUCTURE, RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE) |> 
       is.na() |> summary()

## remove rows with with NAs in these columns
LapF <- Lap |> drop_na(TALL_BOUNDARY_PERCENT, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT)

## Select columns wanted for modelling
LapF <- LapF |> 
           select(est_pairsL, Lap_Density, S_LOC_ID, Landscape, FieldArea, CorvDens,  Reserve, SSSI, ESS_Wader, CSS_Wader, 
                  PropWood_500, PropUrban_500, PropIntense_500, PropWetGrass_500, PropSSSI_500, PropReserve_500, 
                  PropStewWa_500, STANDING_WATER_TOTAL_PERCENT, STOCK, VEG_STRUCTURE, GROUND_DAMP, RUSH_PERCENT, 
                  TALL_BOUNDARY_PERCENT) 

## Set each column as the correct variable type (factor or numeric)
LapF <- LapF |>
            mutate(across(c(est_pairsL, Lap_Density, FieldArea, CorvDens, PropWood_500, PropUrban_500, PropIntense_500, PropWetGrass_500, PropSSSI_500, 
                            PropReserve_500, PropStewWa_500, STANDING_WATER_TOTAL_PERCENT,  RUSH_PERCENT, TALL_BOUNDARY_PERCENT), as.numeric),
                   across(c( Reserve, SSSI, ESS_Wader, CSS_Wader, STOCK, VEG_STRUCTURE, GROUND_DAMP, Landscape), as.factor)) %>% 
           as.data.frame()


## Split the data set into a train-test set, randomly choose rows
# set.seed(1012)
# inTrain <- createDataPartition(
#   y = LapF$est_pairsL, # the outcome data
#   p = 0.7, # The percentage of data in the training set
#   list = FALSE)
# 
# LapF <- LapF |> mutate(LapPres = ifelse(est_pairsL > 0, 1, 0),
#                        LapPres = as.factor(LapPres))
# 
# ## Separate the train and test data set
# Laptrain <- LapF[ inTrain,]
# Lap_test  <- LapF[-inTrain,]


set.seed(1012)
v <- unique(LapF$S_LOC_ID)
Sampv <- sample(v, size = round(length(v)*0.7))

LapF <- LapF |> mutate(LapPres = ifelse(est_pairsL > 0, 1, 0),
                       LapPres = as.factor(LapPres))

Laptrain <- LapF |> filter(S_LOC_ID %in% Sampv)
Lap_test <- LapF |> filter(!S_LOC_ID %in% Sampv)



##-------------##
#### Tune RF ####
##-------------##

LaptrainPres <- filter(Laptrain, Lap_Density > 0)

## Hyper parameter tuning ##
## Finds the optimal "mtry" and "nodesize" tuning parameter for a random forest using out-of-sample error.
## increasing "sampsize" and probably "ntreeTry" increase the accuracy of this assessment
## ranks different parameters based on out-of-sample error
TuneParams <- tune(Lap_Density ~ FieldArea + Landscape + CorvDens + 
                                 STANDING_WATER_TOTAL_PERCENT + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + GROUND_DAMP + RUSH_PERCENT + 
                                 Reserve + SSSI + ESS_Wader + CSS_Wader + PropSSSI_500 + PropReserve_500 + PropStewWa_500 +
                                 PropWood_500 + PropWetGrass_500 + PropUrban_500,
                   data = LaptrainPres,
                   mtryStart = 3, ntreeTry = 1000, sampsize = 250, maxIter = 100,
                   nodesizeTry = c(1:14, seq(15, 40, by = 5)), improve = 5e-4)

TuneParams$optimal # returns optimal parameters from tuning

## plot the surface from parameter tuning
plot.tune(TuneParams)



##------------##
#### Run RF ####
##------------##

## Run random forest regression
RF_fit <- imbalanced(LapPres ~ FieldArea + Landscape + CorvDens + 
                                 STANDING_WATER_TOTAL_PERCENT + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + GROUND_DAMP + RUSH_PERCENT + 
                                 Reserve + SSSI + ESS_Wader + CSS_Wader + PropSSSI_500 + PropReserve_500 + PropStewWa_500 +
                                 PropWood_500 + PropWetGrass_500 + PropUrban_500,
                data = Laptrain, ntree = 10000, method = "rfq", splitrule = "auc")
get.imbalanced.performance(RF_fit)
RF_fit

## Run random forest regression
RF_fit2 <- rfsrc(est_pairsL ~ FieldArea + Landscape + CorvDens + 
                                 STANDING_WATER_TOTAL_PERCENT + TALL_BOUNDARY_PERCENT + STOCK + VEG_STRUCTURE + GROUND_DAMP + RUSH_PERCENT + 
                                 Reserve + SSSI + ESS_Wader + CSS_Wader + PropSSSI_500 + PropReserve_500 + PropStewWa_500 +
                                 PropWood_500 + PropWetGrass_500 + PropUrban_500,
                data = LaptrainPres, ntree = 10000, mtry = 10, nodesize = 2)

RF_fit2 # model summary
## Model performance is displayed in terms of out-of-bag (OOB) prediction error. 
## (OOB) Requested performance error = mean square oob error, mean(real-oob predicted)^2
## (OOB) R squared = R squared of RF model
## Since MSE has scale invariance and lacks interpretation, standardized MSE, defined as the MSE divided by the variance of the outcome
## is used and converted to R squared or the percent of variance explained by a random forest model


## plot how the error rate change with the number of trees
plot(RF_fit)

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
plot.variable(RF_fit, "STANDING_WATER_TOTAL_PERCENT", partial = TRUE)    
plot.variable(RF_fit, "CorvDens", partial = TRUE)
plot.variable(RF_fit, "Landscape", partial = TRUE)
plot.variable(RF_fit, "FieldArea", partial = TRUE)
plot.variable(RF_fit, "PropWetGrass_500", partial = TRUE)

## find possible interactions, not sure quite how this works
find.interaction(RF_fit, method = "vimp", nrep = 3)



##-------------------------##
#### Predict on test set ####
##-------------------------##

## Use the RF model to predict the number of pairs for the test data set
pred_vals <- predict(object = RF_fit, Lap_test)
Lap_test$Predicted <- pred_vals$class
Lap_test$Prob_0 <- pred_vals[["predicted"]][,1]
Lap_test$Prob_1 <- pred_vals[["predicted"]][,2]

## Create confusion matrix, when rounding real and predicted values
example <- confusionMatrix(data= pred_vals$class, 
                           reference = Lap_test$LapPres)
example


## Filter the test set for only the fields that have been classified as 1
Lap_test_Pres <- filter(Lap_test, Predicted == 1)
pred_vals2 <- predict(object = RF_fit2, Lap_test_Pres)
Lap_test_Pres$PredictedAbund <- pred_vals2$predicted
Lap_test_Pres$PredictedAbund2 <- Lap_test_Pres$PredictedAbund*Lap_test_Pres$Prob_1
Absent <- Lap_test_Pres |> filter(LapPres == 0)
Present <- Lap_test_Pres |> filter(LapPres == 1)
sum(Absent$PredictedAbund2);sum(Present$PredictedAbund2)
sum(Lap_test_Pres$PredictedAbund2)
sum(Lap_test$est_pairsL)


Lap_test_Pres$PredictedDens <- pred_vals2$predicted
Lap_test_Pres$PredictedAbund <- Lap_test_Pres$PredictedDens*Lap_test_Pres$FieldArea
Lap_test_Pres$PredictedAbund2 <- Lap_test_Pres$PredictedAbund*Lap_test_Pres$Prob_1

Absent <- Lap_test_Pres |> filter(LapPres == 0)
Present <- Lap_test_Pres |> filter(LapPres == 1)
hist(Absent$PredictedAbund2)
hist(Present$PredictedAbund2)
sum(Lap_test_Pres$PredictedAbund2)
sum(Lap_test$est_pairsL)







pred_vals3 <- predict(object = RF_fit2, Lap_test)
Lap_test$PredictedAbund <- pred_vals3$predicted
Absent <- Lap_test |> filter(LapPres == 0)
Present <- Lap_test |> filter(LapPres == 1)
hist(Absent$PredictedAbund)
hist(Present$PredictedAbund)
sum(Lap_test$PredictedAbund)
sum(Lap_test$est_pairsL)





    
