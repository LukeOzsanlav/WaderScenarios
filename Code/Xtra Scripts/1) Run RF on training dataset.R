## Luke Ozsanlav-Harris
## Train a Random Forest alorithm on my full and resmpled Acceleromter training data set
## Check model performance using 10-fold cross validation and vary mtry to find optimum value

## Created: 22/06/2020


## packages required
library(tidyverse)
library(data.table)
library(randomForest)
library(caret)


## 1. Read in data and prepare for RF
## 2. Split dataset into test-training portions
## 3. Tune mtry parameter in Random forest using 10-fold CV
## 4. Test model accuracy on test set




########
## 1. ##
######## Read in data and prepare for RF

## Read in the training data
Full_train <- fread("Training data sets/New/Full_training_set_with_features.csv")
Resamp_train <- fread("Training data sets/New/Resamp_training_set_with_features.csv")

## extract the columns that need to run the random forest
Full <- Full_train[,c(2:50,53)]
Resamp <- Resamp_train[,c(2:50,53)]

## combine the alert and rest behaviours
Full$Behavior <- ifelse(Full$Behavior == "alert" | Full$Behavior == "rest", 
                        "stationary", Full$Behavior)

## set behaviour to a factor
Full$Behavior <- as.factor(Full$Behavior)
Resamp$Behavior <- as.factor(Resamp$Behavior)




########
## 2. ##
######## Split dataset into test-training portions

## Using a CV + hold out appraoch
## Use a train-test split to test the accuracy of the model and 10-fold CV to tune hyperparamters


## First need to create a train-test split
## Create the data partition
set.seed(1012)
inTrain <- createDataPartition(
  y = Full$Behavior, # the outcome data
  p = 0.7, # The percentage of data in thetraining set
  list = FALSE
)

## Seperate the train and test set
training <- Full[ inTrain,]
testing  <- Full[-inTrain,]




########
## 3. ##
######## Tune mtry parameter in Random forest using 10-fold CV

#### using the caret package to tune the RF parameter mtry (only mtry available for tuning) ####

## seeting controls to 3 repeats of 10-fold CV
train_control <- trainControl(method="repeatedcv", number=10, repeats=3, search= "grid")
## try a range of mtry values
tunegrid <- expand.grid(.mtry=c(1:12))
## run random forest
set.seed(1012)
RF_tune <- train(Behavior ~ . , data = training,
            method = "rf", trControl = train_control, tuneGrid=tunegrid, ntree= 5000)
print(RF_tune) # model ouput
plot(RF_tune, main = "Tune mtry parameter with 10 fold CV x3, ntree = 5000, 4 Vars") # plot accuracy vs mtry value

## Output:
## The final value used for the model was mtry = 6
## mtry= 5, accuracy=  0.9882482  , kappa=  0.9827883
## Note little change in accuracy of mtry 3 to 11, only 0.1% diff in accuracy between them



#### using tuning function from RandomForest package to find best mtry value ####

set.seed(1012)
best_mtry <- tuneRF(x = training[,1:49], y =training$Behavior, mtryStart = 2, imporve = 0.00001, stepFactor = 1.5, ntree= 5000)
best_mtry ## get tuning output

## Ouput: 23/06/2020
## Best mtry value was 3 with an OOB error was 1.19%


## Code for defining your own random forest algorithm that can then allow ntree and mtry to both vary
## Found here: https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/


## now run random forest with a specifc mtry value if need be
## if all above approaches are similar then may just go with caret tuning
#set.seed(1012)
#RF <- randomForest(Behavior ~ . , data = training, method = "rf", ntree= 1000, mtry =6, classwt = c(0.3, 0.3, 0.3, 0.1))
#print(RF) # print model summary


## save the trained model for later use
#saveRDS(RF_tune, file = "Random Forest Model/RF_model_full_4behavs_walk.rds") #save as rds
RF_tune <- readRDS("Random Forest Model/RF_model_full4w.rds") # read in the model



########
## 4. ##
######## Test model accuracy on test set

## predict on the testing data set
predClasses <- predict(RF_tune, newdata = testing)

## create confusion matirx to test performance of the model
confusionMatrix(data = predClasses, testing$Behavior)

## Output: 23/06/2020
## Accuracy : 0.9847          
## 95% CI : (0.9767, 0.9905)
## Confusion Matrix:
##              Reference:
## Prediction:   flight graze stationary walk
## flight        520     7          0    0
## graze           5   493          0    3
## stationary      0     0        223    1
## walk            0     0          5  113


## Output: 26/06/2020
## Changed the training data, more stringent on selection of walking behaviour

## Accuracy : 0.992          
## 95% CI : (0.9857, 0.996)
## Confusion Matrix:
##              Reference
## Prediction   flight graze stationary walk
## flight        522     5          0    0
## graze           3   501          0    0
## stationary      0     0        225    1
## walk            0     0          2  111


## Also check which variables were most important for the classification
RFimp <- varImp(RF_tune);RFimp
plot(RFimp)




























########
## x. ##
######## Potentially important inputs for Random forest algorithm ####

##norm.votes	
##If TRUE (default), the final result of votes are expressed as fractions. 
##If FALSE, raw vote counts are returned (useful for combining results from different runs).

##importance	
##Should importance of predictors be assessed?

##strata	
##A (factor) variable that is used for stratified sampling.

##sampsize	
##Size(s) of sample to draw. For classification, if sampsize is a vector of the length the number of strata, 
##then sampling is stratified by strata, and the elements of sampsize indicate the numbers to be drawn from the strata.

##replace	
##Should sampling of cases be done with or without replacement?

##subset	
##an index vector indicating which rows should be used.
##list of rows to be used to train the model, can just be a sibset of the data frame provided
