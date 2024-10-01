## Luke Ozsanlav-Harrs

## Calculate the Acc features for all the non-training acc data
## then use the train Random forest algorithm to classify the remainder of the behaviors

## Created: 23/06/2020

## packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(randomForest)
library(caret)


## 1. Reading in Ornitela Acc data and cleaning
## 2. Preparing data for calculating Acc features
## 3. Calculate Acc features
## 4. Use RF algorithm to predict behaviors




########
## 1. ##
######## Reading in Ornitela Acc data and cleaning

## set wroking directory
#setwd("~/Ornitela data/ALL_ORNI_DATA/Cleaned CSVs")

## read in data
#orn_data <- fread("IRE4_clean.csv")
orn_data <- readRDS(file = "AccData/Ornitela_ACCClean.RDS")
## Try one of the two below ways to read in multi files quicker on the servers
#orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))
#orn_data$V1 <- NULL

## alternate method for reading in csv list
#csv.list <- list.files(pattern = "*.csv")
#orn_data <- setDT(ldply(csv.list, fread, .progress="text"))

## sub-setting out all of the acc data
#orn_data <- filter(orn_data, datatype == "SENSORS")

## needed to prep data timestamp
orn_data$UTC_datetime <- as.POSIXct(as.character(ymd_hms(orn_data$UTC_datetime)), format= "%Y-%m-%d %H:%M:%S", tz= "GMT")
orn_data <- orn_data %>% mutate(month=month(UTC_datetime), year = year(UTC_datetime))




########
## 2. ##
######## Preparing data for calculating Acc features

#### LABELLING ALL THE BURSTS SO THEY HAVE A UNIQUE CODE #### 

## creating a time difference column
orn_data$timedif <- orn_data$UTC_datetime - lag(orn_data$UTC_datetime)

## changing all gaps between bursts to 1 and within burst gaps to 0 
#(problem: sometimes the time diffs go back in time where tags join together)
orn_data$timedif <- ifelse(orn_data$timedif < 0, (orn_data$timedif*-1) , orn_data$timedif) ##solve problem
orn_data$timedif <- ifelse(orn_data$timedif > 1, 1, 0) ##change to all 0 and 1s
orn_data$timedif[1] <- 1 ##first line is an NA so changing to 1 as it starts a burst

## now use cumsum to number the bursts, as there is a 1 at each burst boundary burst will be numbered increasingly in steps of 1
orn_data$burst_no <- cumsum(orn_data$timedif)



#### CHECKING THE ACC DATA IS CLEAN ####

## creating column with burst no and ID
orn_data$burst_no_ID <- paste(orn_data$burst_no, orn_data$device_id, sep = "_")
orn_data$dummy <- 1 ##dummy column

## creating data set with stats on each burst
Acc_data_checks <- orn_data %>% group_by(burst_no_ID) %>% summarise(
  Readings_in_burst = sum(dummy), end_time = max(UTC_datetime), 
  start_time =min(UTC_datetime), deviceID = min(device_id)) 
Acc_data_checks$burst_length_in_secs <- as.numeric(Acc_data_checks$end_time - Acc_data_checks$start_time) ##burst length(time)
Acc_data_checks$burst_freq <- Acc_data_checks$Readings_in_burst/Acc_data_checks$burst_length_in_secs ##burst_freq

## summary of checks
table(Acc_data_checks$Readings_in_burst) ##should all be 30
table(Acc_data_checks$burst_freq) ## should all be 10 and 15
sum(subset(Acc_data_checks, select = c("deviceID", "start_time")) %>% duplicated) ## should be 0
rm(Acc_data_checks) ## remove this from environment




########
## 3. ##
######## Calculate Acc features

## need to convert the acc to g by dividung by 2048 (Only for Orni Data)
orn_data$X <- ((orn_data$acc_x)/2048)
orn_data$Y <- ((orn_data$acc_y)/2048)
orn_data$Z <- ((orn_data$acc_z)/2048)

## streamline the data set so will run slightly quicker
Orn_stream <- subset(orn_data, select = c("burst_no_ID", "X", "Y", "Z"))
rm(orn_data)


#### My own functions to calculate some of the Acc data features ####
## Function 1: norm of a vector
norm_vecf <- function(x){sqrt(sum(x^2))} 

## Function 2 Caulcalate number of line crossings/number of rows
line_crossing_calc <- function(axis1, axis2) {
  
  dummy <- rep(NA, 30)
  dummy <- ifelse((axis1 - axis2) > 0 & lead(axis1 - axis2) < 0, 1, 
                  ifelse((axis1 - axis2) < 0 & lead(axis1 - axis2) > 0, 1, 0))
  return(sum(dummy, na.rm = T)/30)
  
}

## Function 3: Calclate local wave amplitude: diff between min and max every 10 readings
local_wave_amp <- function(axis) {
  
  mean(c(max(axis[1:10])- min(axis[1:10]), 
         max(axis[11:20])- min(axis[11:20]), 
         max(axis[21:30])- min(axis[21:30])))
}



## function to calculate acc features per burst using data.table
Acc_feautures <- function(data) {
  
  # change data to a data.table
  data <- as.data.table(data)
  
  # calculate the dynamic acceleration value for each axis and sum all for ODBA
  output <- data[,.(dba_x = mean(abs(X - mean(X))),
                    dba_y = mean(abs(Y - mean(Y))),
                    dba_z = mean(abs(Z - mean(Z))),
                    odba = (mean(abs(X - mean(X))) + mean(abs(Y - mean(Y))) + mean(abs(Z - mean(Z)))),
                    mean_x = mean(X),
                    mean_y = mean(Y),
                    mean_z = mean(Z),
                    min_x = min(X),
                    min_y = min(Y),
                    min_z = min(Z),
                    max_x = max(X),
                    max_y = max(Y),
                    max_z = max(Z),
                    sd_x = sd(X),
                    sd_y = sd(Y),
                    sd_z = sd(Z),
                    skewness_x = moments::skewness(X),
                    skewness_y = moments::skewness(Y),
                    skewness_z = moments::skewness(Z),
                    kurtosis_x = e1071::kurtosis(X, type = 1),
                    kurtosis_y = e1071::kurtosis(Y, type = 1),
                    kurtosis_z = e1071::kurtosis(Z, type = 1),
                    Q25_x = quantile(X, probs= 0.25),
                    Q25_y = quantile(Y, probs= 0.25),
                    Q25_z = quantile(Z, probs= 0.25),
                    Q50_x = quantile(X, probs= 0.50),
                    Q50_y = quantile(Y, probs= 0.50),
                    Q50_z = quantile(Z, probs= 0.50),
                    Q75_x = quantile(X, probs= 0.75),
                    Q75_y = quantile(Y, probs= 0.75),
                    Q75_z = quantile(Z, probs= 0.75),
                    mean_diffxy = mean(X-Y),
                    mean_diffxz = mean(X-Z),
                    mean_diffyz = mean(Y-Z),
                    sd_diffxy = sd(X-Y),
                    sd_diffxz = sd(X-Z),
                    sd_diffyz = sd(Y-Z),
                    #cor_xy = cor(X, Y, method=c("pearson")),
                    #cor_xz = cor(X, Z, method=c("pearson")),
                    #cor_yz = cor(Y, Z, method=c("pearson")),
                    cov_xy = cov(X, Y, method=c("pearson")),
                    cov_xz = cov(X, Z, method=c("pearson")),
                    cov_yz = cov(Y, Z, method=c("pearson")),
                    norm_vec_x = norm_vecf(X),
                    norm_vec_y = norm_vecf(Y),
                    norm_vec_z = norm_vecf(Z),
                    line_crossingsxy = line_crossing_calc(axis1 = X, axis2 = Y),
                    line_crossingsxz = line_crossing_calc(axis1 = X, axis2 = Z),
                    line_crossingsyz = line_crossing_calc(axis1 = Y, axis2 = Z),
                    wave_amp_x = local_wave_amp(axis = X),
                    wave_amp_y = local_wave_amp(axis = Y),
                    wave_amp_z = local_wave_amp(axis = Z)),
                 by = .(burst_no_ID)]
}


## calculate Acc features for each burst using my function
system.time(output <- Acc_feautures(Orn_stream))
## time for IRE1 file: 905.108 secs

## going to drop rows that have any NA values in
output <- drop_na(output)
# In future may want to not drop rows if it becomes large
# NAs occur due to all values in the z axis in a burst being the same
# Code below sets NAs in Skewness and kurtosis z columns to zero
# Not sure how to deal with NAs in cor_yz and cor_xz columns
summary(output)
output$skewness_z <- ifelse(is.na(output$skewness_z) == TRUE, 0, output$skewness_z)
output$kurtosis_z <- ifelse(is.na(output$kurtosis_z) == TRUE, 0, output$kurtosis_z)




########
## 4. ##
######## Use RF algorithm to predict behaviours

## import that model trainied in another script
#setwd("~/Classification of Acc data/Random Forest Model")
RF <- readRDS("Random Forest Model/RF_model_full_4behavs_walk.rds") # read in the model

## use predict to classify the burst to a behavior
Classified_acc <- output[,1]
system.time(Classified_acc$Behaviour <- predict(RF, newdata = output[, 2:50]))
## time for IRE1 file: 447 secs


## bind back behavior with the rest of the info on that burst, i.e. timestamp
## Extract only one line per burst for both data sets
line_extractor <- function(data){
  
  dups <- subset(data, select = c("burst_no_ID")) %>% duplicated
  sum(dups) ##how many rows are duplicates 
  burst_info <- data %>% filter(!dups) ##now remove these duplicate rows
  ## select the rows I want
  return(subset(burst_info, select = c("device_id", "UTC_datetime", "burst_no_ID")))
  
}

## extract line per burst
burst_info <- line_extractor(orn_data)

## now join on the extra info for each burst
Acc_behaviours <- inner_join(Classified_acc, burst_info, by = "burst_no_ID")
stopifnot(nrow(Acc_behaviours) == nrow(Classified_acc))

## check number of each behaviour
table(Acc_behaviours$Behaviour)

## Ouput: 28/06/2020 for IRE1 file
## Without cor feautre and walks changed    
## flight     graze   stationary     walk 
## 49784     255174     335098      88534 

## With cor and without walk changed    
## flight     graze   stationary     walk 
## 47120     238009     331852     111609



## Read out the data for later use
#setwd("~/Classification of Acc data/Classified datasets")
write.csv(Acc_behaviours, file = "Classified datasets/Classified_ACC.csv", row.names = F)
saveRDS(Acc_behaviours, "Classified datasets/Classified_ACC.RDS")


