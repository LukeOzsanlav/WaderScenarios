##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 12/10/2023
## Goal: Combine and clean BWWM survey data for 2021/22/23 producing,
##       Breeding pair estimate per field and summarize data per field
## 
##------------------------------------------------------##



## Load in required packages
pacman::p_load(here, tidyverse, data.table, sf)
source("Code/Helper functions.R") # get helper functions

## Use `here` package for specifying file paths
here::i_am("Code/Wader Abundance/1 - Cleaning Script.R")
library(here)
here()

## Queries:
## Why some fields only have one entry in the field data but then visit number is 2, 3 or 4??
## Over 1000 fields were only one visit and they were covered and no entry in other habitat column and no corresponding bird data
## WHy are there lots of "B"s in the VEG_HEIGHT column? Should be "S" or "T"



##---------------------##
#### 0. Data Read In ####
##---------------------##

##------------------------------------##
#### 0.1 Read in field habitat data ####
##------------------------------------##

## 2021 field habitat data
Field_21 <- read_csv(here("RawData", "Field data", "2021", "bwwm-visits-field@16-11-2021.csv"))

## 2022 field habitat data
Field_22 <- read_csv(here("RawData", "Field data", "2022", "bwwm-visits-field@13-09-2022.csv"))

## 2023 field habitat data 
## Make start and end times be a character, reading in a time can corrupt data
Field_23_GT <- read_csv(here("RawData", "Field data", "2023", "Greater Thames", "Field habitat data entry Greater Thames 17 08 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Field_23_NB <- read_csv(here("RawData", "Field data", "2023", "Norfolk Broads", "Field habitat data entry Broads 11 09 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Field_23_SL1 <- read_csv(here("RawData", "Field data", "2023", "Somerset Levels", "Field habitat data entry Somerset 17 08 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Field_23_SL2 <- read_csv(here("RawData", "Field data", "2023", "Somerset Levels", "Field habitat data entry WSM and KSM 17 08 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))

## Add origin column to each data set so I can keep track of where they came from
Field_21$Origin <- "Field_21"
Field_22$Origin <- "Field_22"
Field_23_GT$Origin <- "Field_23_GT"
Field_23_NB$Origin <- "Field_23_NB"
Field_23_SL1$Origin <- "Field_23_SL1"
Field_23_SL2$Origin <- "Field_23_SL2"



##----------------------------------##
#### 0.2 Read in bird survey data ####
##----------------------------------##

## 2021 bird survey data
Bird_21 <- read_csv(here("RawData", "Bird data", "2021", "bwwm-obs@16-11-2021.csv"))

## 2022 bird survey data
Bird_22 <- read_csv(here("RawData", "Bird data", "2022", "bwwm-obs@13-09-2022.csv"))

## 2023 bird survey data
## ## Make start and end times be a character, reading in a time can corrupt data
Bird_23_GT <- read_csv(here("RawData", "Bird data", "2023", "Greater Thames", "Bird data entry Greater Thames 17 08 2023 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Bird_23_NB <- read_csv(here("RawData", "Bird data", "2023", "Norfolk Broads", "Bird data entry Broads 11 09 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Bird_23_SL1 <- read_csv(here("RawData", "Bird data", "2023", "Somerset Levels", "Bird data entry Somerset 17 08 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))
Bird_23_SL2 <- read_csv(here("RawData", "Bird data", "2023", "Somerset Levels", "Bird data entry WSM and KSM 17 08 23 FINAL.csv"), col_types = cols(START_TIME = col_character(), END_TIME = col_character()))

## Add origin column to each data set so I can keep track of where they came from
Bird_21$Origin <- "Bird_21"
Bird_22$Origin <- "Bird_22"
Bird_23_GT$Origin <- "Bird_23_GT"
Bird_23_NB$Origin <- "Bird_23_NB"
Bird_23_SL1$Origin <- "Bird_23_SL1"
Bird_23_SL2$Origin <- "Bird_23_SL2"



##---------------------------------##
#### 0.3 Read in field shapefile ####
##---------------------------------##

## Read in the shapefile for all fields
FieldShapes <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp"))

## Check all F_LOC_IDs are unique, should return TRUE
length(unique(FieldShapes$F_LOC_ID)) == nrow(FieldShapes)

## Check no F_LOC_IDs are missing, should all be FALSE
table(is.na(FieldShapes$F_LOC_ID))

## Extract the data to work out F_LOC_ID from S_LOC_ID & FIELD_NUMBER 
FieldConv <- FieldShapes |> 
  select(c(S_LOC_ID, F_LOC_ID, Field_no)) |> 
  st_drop_geometry() |> 
  dplyr::rename("FIELD_NUMBER" = "Field_no") |> 
  mutate(FIELD_NUMBER = as.numeric(FIELD_NUMBER))




##----------------------------##
#### 1. Field Data Cleaning ####
##----------------------------##

##-----------------------------------------##
#### 1.1 Add F_LOC_ID to 2023 field data ####
##-----------------------------------------##  

## Fix column problem in Field_23_SL2 data set
## Many columns that not present in the other 2023 dataset
Field_23_SL2 <- Field_23_SL2 |> 
  select(-c(RUSH_PERCENT...17, SWARD_RANGE, SWARD_HEIGHT_DOMINANCE, TUSSOCK_PERCENT, TUSSOCK_HEIGHT, STANDING_WATER_PERCENT,
            DIS_VEG_HEIGHT, DIS_VEG_PATCH, DIS_VEG_AREA_PERCENT, DIS_VEG_AREA_PERCENT, VEG_HEIGHT)) |> 
  rename(RUSH_PERCENT = RUSH_PERCENT...31)
## check which columns are missing in Field_23_SL2 data set
misscols <- match(colnames(Field_23_SL1), colnames(Field_23_SL2))
colnames(Field_23_SL1)[is.na(misscols) == T]

## Bind together all 2023 data
## NOTE** need to work out why Field_23_SL2 is different and then add it back in
Field_23 <- rbind(Field_23_GT, Field_23_NB, Field_23_SL1) |> plyr::rbind.fill(Field_23_SL2)

## Now add F_LOC_ID to field data, with check for data set length
lpre <- nrow(Field_23)
Field_23 <- left_join(Field_23, FieldConv, by = c("S_LOC_ID", "FIELD_NUMBER"))
stopifnot(nrow(Field_23) == lpre)



##------------------------------##
#### 1.2 Bind all: field data ####
##------------------------------## 

##** Bind 2021 / 2022 data together
## remove this column from the 2021 data to make it match with the 2022 data
Field_21 <- Field_21 |> select(-UniF_Code)
## Bind 2021 and 2022 data sets now
Field21_22 <- rbind(Field_21, Field_22)


##** Bind on 2023
## Work out which colnames to change
test <- match(colnames(Field_23), colnames(Field21_22))
extracols <- Field_23[,(is.na(test)==T)]
colnames(extracols)
glimpse(extracols)

test2 <- match(colnames(Field21_22), colnames(Field_23))
extracols2 <- Field21_22[,(is.na(test2)==T)]
colnames(extracols2)
glimpse(extracols2)

## Change col names in 2021/22 data set to match 2023 dataset
Field21_22 <- Field21_22 |> 
  rename(VISIT_NUMBER =  VISIT,
         STANDING_WATER_TOTAL_PERCENT = STANDING_WATER_PERCENT)
## add 2023 data
Field_All <- plyr::rbind.fill(Field21_22, Field_23)
  


##--------------------------------------##
#### 1.3 General cleaning: field data ####
##--------------------------------------## 

##** Remove any entire rows missing all data
preFilt <- nrow(Field_All) # record rows at start
Field_All <- filter(Field_All, rowSums(is.na(Field_All)) != (ncol(Field_All)-1))
preFilt - nrow(Field_All) # change in row number


##** Remove any rows missing F_LOC_IDs
preFilt <- nrow(Field_All) # record rows at start
table(is.na(Field_All$F_LOC_ID)) # check NAs
## If fine to remove rows missing F_LOC_IDs by filtering them out
Field_All <- filter(Field_All, is.na(F_LOC_ID)==FALSE)
preFilt - nrow(Field_All) # change in row number


##** Remove any row duplicates
## select columns, find duplicates and remove them
dups <- Field_All |> 
  select(F_LOC_ID, VISIT_NUMBER, VISIT_DATE) %>%
  duplicated()
sum(dups) # check number if duplicates
Field_All <- Field_All |> filter(!dups) # remove duplicates


##** Remove rows with key data missing
## fix missing dates in FIELD_SL2 part of data 
Field_All <- Field_All %>% mutate(START_TIME = ifelse(is.na(VISIT_DATE)==T & COVERAGE == "N" & VISIT_NUMBER == 1, 
                                  "06:00", START_TIME),
                                  VISIT_DATE = ifelse(is.na(VISIT_DATE)==T & COVERAGE == "N" & VISIT_NUMBER == 1, 
                                  "20.04.2023", VISIT_DATE))
miss <- complete.cases(Field_All[, c("F_LOC_ID", "VISIT_DATE")])
table(miss) # number rows with missing info
Field_All <- Field_All %>% filter(miss)



##-----------------------------------------------##
#### 1.4 Standardize column syntax: field data ####
##-----------------------------------------------## 

## For GRASSLAND_TYPE is the value is not S or I then change to NA
Field_All <- Field_All |> mutate(GRASSLAND_TYPE = ifelse(GRASSLAND_TYPE == "I" | GRASSLAND_TYPE == "S", GRASSLAND_TYPE, NA))

## ** Fix number of inconsistencies in recording listed below
## Change "P" to "Y in PREDATOR_FENCING 
## Change "S" to "Y" in GROUND_COND_UNKNOWN
## Change "NR" to NAs in VEG_STRUCTURE 
Field_All <- Field_All |> mutate(PREDATOR_FENCING = ifelse(PREDATOR_FENCING == "P", "Y", PREDATOR_FENCING),
                                 GROUND_COND_UNKNOWN = ifelse(GROUND_COND_UNKNOWN == "S", "Y", GROUND_COND_UNKNOWN),
                                 VEG_STRUCTURE = ifelse(VEG_STRUCTURE == "NR", NA, VEG_STRUCTURE))
table(Field_All$PREDATOR_FENCING) # check change
table(Field_All$GROUND_COND_UNKNOWN) # check change
table(Field_All$VEG_STRUCTURE) # check change


##** Error in VISIT_NUMBER, some have wrong number recording
## Will change these with ifelse as there is only three specific mistakes
Field_All <- Field_All %>% 
  mutate(VISIT_NUMBER = ifelse(VISIT_NUMBER == "2 (1)", 2, VISIT_NUMBER),
         VISIT_NUMBER = ifelse(VISIT_NUMBER == "3 (1)", 3, VISIT_NUMBER),
         VISIT_NUMBER = ifelse(VISIT_NUMBER == "s", 2, VISIT_NUMBER))
table(Field_All$VISIT_NUMBER) # check change


##** Some columns have mixed cases so make all upper case
Field_All <- Field_All %>% 
  mutate(across(c(ALL_DITCHES_ASSESSED, VEG_STRUCTURE, COVERAGE), toupper)) |> 
  mutate(ALL_DITCHES_ASSESSED = ifelse(ALL_DITCHES_ASSESSED == "N/A", NA, ALL_DITCHES_ASSESSED)) # NA been input wrongly so change


##** Some column contain a > or < followed by a 1, change all these to 0.5
## Create function that does this change
SymbRem <- function(x){
  Sub <- (substr(x, 0, 1) == "<")
  x[which(Sub == TRUE)] <- 0.5
  Sub2 <- (substr(x, 0, 1) == ">")
  x[which(Sub2 == TRUE)] <- 0.5
  return(x)
}
## apply the function to all the columns that had this error
Field_All <- Field_All %>% 
  mutate(across(c(STANDING_WATER_IN_POOLS_PERCENT, STANDING_WATER_IN_FOOT_DRAINS_PERCENT, 
                  STANDING_WATER_TOTAL_PERCENT, TALL_BOUNDARY_PERCENT, RUSH_PERCENT), SymbRem))
table(Field_All$RUSH_PERCENT) # check change


##** STANDING_WATER_TOTAL_PERCENT & RUSH_PERCENT has "x to y" so change these to actual numbers
Field_All <- Field_All |> 
  mutate(RUSH_PERCENT = case_when(RUSH_PERCENT == "NR" ~ NA,
                                  RUSH_PERCENT == "0 to 5" ~ "2.5",
                                  RUSH_PERCENT == "5 to 10" ~ "7.5",
                                  .default = RUSH_PERCENT),
         RUSH_PERCENT = as.numeric(RUSH_PERCENT),
         STANDING_WATER_TOTAL_PERCENT = case_when(STANDING_WATER_TOTAL_PERCENT == "NR" ~ NA,
                                                  STANDING_WATER_TOTAL_PERCENT == "10 to 20" ~ "15",
                                                  STANDING_WATER_TOTAL_PERCENT == "15 to 20" ~ "17.5",
                                                  STANDING_WATER_TOTAL_PERCENT == "60 to 70" ~ "65",
                                                  STANDING_WATER_TOTAL_PERCENT == "5 to 10" ~ "7.5",
                                                  STANDING_WATER_TOTAL_PERCENT == "70 TO 80" ~ "75",
                                                  STANDING_WATER_TOTAL_PERCENT == "70 to 80" ~ "75",
                                                  STANDING_WATER_TOTAL_PERCENT == "80 to 90" ~ "85",
                                                  .default = STANDING_WATER_TOTAL_PERCENT),
         STANDING_WATER_TOTAL_PERCENT = as.numeric(STANDING_WATER_TOTAL_PERCENT))
table(Field_All$RUSH_PERCENT)# check change
table(Field_All$STANDING_WATER_TOTAL_PERCENT)# check change



##** Some start and end dates have a "<" instead of a ":"
## ## Create function that does this change
SymbRem2 <- function(x){
  Sub <- (substr(x, 3, 3) == "<")
  x[which(Sub == TRUE)] <- paste0(substr(x[which(Sub == TRUE)], 0, 2), ":", substr(x[which(Sub == TRUE)], 4, 5))
  Sub2 <- (substr(x, 3, 3) == ">")
  x[which(Sub2 == TRUE)] <- paste0(substr(x[which(Sub2 == TRUE)], 0, 2), ":", substr(x[which(Sub2 == TRUE)], 4, 5))
  Sub3 <- (substr(x, 2, 2) == ">")
  x[which(Sub3 == TRUE)] <- paste0(substr(x[which(Sub3 == TRUE)], 0, 1), ":", substr(x[which(Sub3 == TRUE)], 3, 4))
  return(x)
}
## apply change to stand and end times
Field_All <- Field_All %>% mutate(across(c(START_TIME, END_TIME), SymbRem2))
table(Field_All$END_TIME); table(Field_All$START_TIME) # check change


##** Visit Dates differ between 2021/22 and 2023 so need to standardize format
## split data so that I can fix the problem
Sub21_22 <- filter(Field_All, Origin %in% c("Field_21", "Field_22"))
Sub_23 <- filter(Field_All, Origin %in% c("Field_23_GT", "Field_23_NB", "Field_23_SL1",  "Field_23_SL2"))

## define list of possible date times in the data set
datetime_formats <- c("dmY HMS", "dmY HM", "Ymd HMS") #specify date & time format 

## format datetime in 2021/22
Sub21_22 <- Sub21_22 |> mutate(VISIT_DATE = parse_date_time(VISIT_DATE, orders=datetime_formats))
table(is.na(Sub21_22$VISIT_DATE)) # check this worked
## Some 2022 data had the wrong data (2021) so fix these
Sub21_22 <- Sub21_22 |> 
  mutate(year = year(VISIT_DATE),
         VISIT_DATE = ifelse(Origin == "Field_22" & year == 2021, paste(VISIT_DATE + years(1)), paste(VISIT_DATE)),
         VISIT_DATE = ymd_hms(VISIT_DATE)) |> 
  select(-year)

## format datetime in 2023
Sub_23 <- Sub_23 |> mutate(VISIT_DATE = parse_date_time(paste(VISIT_DATE, START_TIME), orders=datetime_formats))
table(is.na(Sub_23$VISIT_DATE)) # check this worked
# Check that all years are 2023 and change if not
Sub_23 <- Sub_23 |> 
  mutate(year = year(VISIT_DATE),
         VISIT_DATE = ifelse(!year == 2023, paste0("2023", substr(VISIT_DATE, 5, 19)), paste0(VISIT_DATE)),
         VISIT_DATE = ymd_hms(VISIT_DATE)) |> 
  select(-year)
table(is.na(Sub_23$VISIT_DATE)); table(year(Sub_23$VISIT_DATE)) # check this worked

## bind the two halves of the data back together
Field_All <- rbind(Sub21_22, Sub_23)






##---------------------------##
#### 2. Bird Data Cleaning ####
##---------------------------##

##-----------------------------------------##
#### 2.1 Add F_LOC_ID to 2023 bird data ####
##-----------------------------------------##

## Remove these columns from Somerset data, not needed
Bird_23_SL1 <- Bird_23_SL1 |> select(-`ADD TO WSM DATASET`)
Bird_23_SL2 <- Bird_23_SL2 |> select(-`...18`)

## One of the S_LOC_IDs is wrong so fix this
Bird_23_GT <- Bird_23_GT |> 
  mutate(S_LOC_ID = ifelse(SITE_NAME== "Stoke B" & SITE_CODE == "B_977", "LOC2990268", S_LOC_ID))

## Bind all the 2023 data
Bird_23 <- rbind(Bird_23_GT, Bird_23_NB, Bird_23_SL1, Bird_23_SL2)

## Now add F_LOC_ID to bird survey data, with check for data set length
lpre <- nrow(Bird_23)
Bird_23 <- left_join(Bird_23, FieldConv, by = c("S_LOC_ID", "FIELD_NUMBER"))
stopifnot(nrow(Bird_23) == lpre)



##-----------------------------##
#### 2.2 Bind all: bird data ####
##-----------------------------##

##** Bind 2021 / 2022 bird survey data together
Bird21_22 <- rbind(Bird_21, Bird_22)

##** Bind on 2023
## Work out which colnames to change
test <- match(colnames(Bird_23), colnames(Bird21_22))
extracols <- Bird_23[,(is.na(test)==T)]
colnames(extracols)
glimpse(extracols)

test2 <- match(colnames(Bird21_22), colnames(Bird_23))
extracols2 <- Bird21_22[,(is.na(test2)==T)]
colnames(extracols2)
glimpse(extracols2)

## Change col names in 2021/22 data set to match 2023 dataset
Bird21_22 <- Bird21_22 |> rename(VISIT_NUMBER =  VISIT)
## add 2023 data
Bird_All <- plyr::rbind.fill(Bird21_22, Bird_23)



##-------------------------------------##
#### 2.3 General cleaning: bird data ####
##-------------------------------------##

##** Remove any entire rows missing all data
preFilt <- nrow(Bird_All) # record rows at start
Bird_All <- filter(Bird_All, rowSums(is.na(Bird_All)) != (ncol(Bird_All)-1))
preFilt - nrow(Bird_All) # change in row number


##** Remove any rows missing F_LOC_IDs
preFilt <- nrow(Bird_All) # record rows at start
table(is.na(Bird_All$F_LOC_ID)) # check NAs
## Read out these rows for Rob to check 
Bird_MissFLOC <- filter(Bird_All, is.na(F_LOC_ID)==TRUE)
# write.csv(Bird_MissFLOC, "Bird Data with missing FLOCID.csv")
## If fine to remove rows missing F_LOC_IDs by filtering them out
Bird_All <- filter(Bird_All, is.na(F_LOC_ID)==FALSE)
preFilt - nrow(Bird_All) # change in row number


##** Remove any row duplicates
## select columns, find duplicates and remove them
dups <- Bird_All |> 
  select(F_LOC_ID, VISIT_NUMBER, VISIT_DATE, ENGLISH_NAME) %>%
  duplicated()
sum(dups) # check number if duplicates
Bird_All <- Bird_All |> filter(!dups) # remove duplicates


##** Remove rows with key data missing
miss <- complete.cases(Bird_All[, c("F_LOC_ID", "VISIT_DATE")])
table(miss) # number rows with missing info
Bird_All <- Bird_All %>% filter(miss)



##-----------------------------------------------##
#### 2.4 Standardize column syntax: bird data ####
##-----------------------------------------------## 

##** Change Visit 4B to Visit 4
Bird_All <- Bird_All |> 
  mutate(VISIT_NUMBER = ifelse(VISIT_NUMBER == "4B", "4", VISIT_NUMBER))

##** Remove "2?" from estimates pairs columns
Bird_All <- Bird_All |> 
  mutate(ESTIMATED_PAIRS = ifelse(ESTIMATED_PAIRS == "2?", "2", ESTIMATED_PAIRS))
table(Bird_All$ESTIMATED_PAIRS) # check change


##** Change suspected in YOUNG_PRESENT to Y
Bird_All <- Bird_All |> 
  mutate(YOUNG_PRESENT = ifelse(YOUNG_PRESENT == "suspected", "Y", 
                                ifelse(YOUNG_PRESENT == "(suspected)", "Y", YOUNG_PRESENT)))
table(Bird_All$YOUNG_PRESENT) # check change


##** Some start and end times have a "<" instead of a ":"
## ## Create function that does this change
SymbRem2 <- function(x){
  Sub <- (substr(x, 3, 3) == "<")
  x[which(Sub == TRUE)] <- paste0(substr(x[which(Sub == TRUE)], 0, 2), ":", substr(x[which(Sub == TRUE)], 4, 5))
  Sub2 <- (substr(x, 3, 3) == ">")
  x[which(Sub2 == TRUE)] <- paste0(substr(x[which(Sub2 == TRUE)], 0, 2), ":", substr(x[which(Sub2 == TRUE)], 4, 5))
  Sub3 <- (substr(x, 2, 2) == ">")
  x[which(Sub3 == TRUE)] <- paste0(substr(x[which(Sub3 == TRUE)], 0, 1), ":", substr(x[which(Sub3 == TRUE)], 3, 4))
  return(x)
}
## apply change to stand and end times
Bird_All <- Bird_All %>% mutate(across(c(START_TIME, END_TIME), SymbRem2))
table(Bird_All$END_TIME); table(Bird_All$START_TIME)  # check the fix has worked


##** Visit Dates differ between 2021/22 and 2023 so need to standardize format
##*split data so that I can fix the problem
Sub21_22 <- filter(Bird_All, Origin %in% c("Bird_21", "Bird_22"))
Sub_23 <- filter(Bird_All, Origin %in% c("Bird_23_GT", "Bird_23_NB", "Bird_23_SL1",  "Bird_23_SL2"))

## define list of possible date times in the data set
datetime_formats <- c("dmY HMS", "dmY HM", "Ymd HMS", "dm hM", "Ydm HM", "Ydm hM") #specify date & time format 

## format datetime in 2021/22
Sub21_22 <- Sub21_22 |> 
  mutate(VISIT_DATE = parse_date_time(VISIT_DATE, orders=datetime_formats))
table(is.na(Sub21_22$VISIT_DATE)) # check this worked
## Some 2022 data had the wrong data (2021) so fix these
Sub21_22 <- Sub21_22 |> 
  mutate(year = year(VISIT_DATE),
         VISIT_DATE = ifelse(Origin == "Bird_22" & year == 2021, paste(VISIT_DATE + years(1)), paste(VISIT_DATE)),
         VISIT_DATE = ymd_hms(VISIT_DATE)) |> 
  select(-year)

## format datetime in 2023
##Also correct one date typo
Sub_23 <- Sub_23 |> mutate(VISIT_DATE = ifelse(VISIT_DATE == "15.06", "15.06.2023", VISIT_DATE)) 
Sub_23 <- Sub_23 |> 
  mutate(VISIT_DATE = parse_date_time(paste(VISIT_DATE, START_TIME), orders=datetime_formats))
table(is.na(Sub_23$VISIT_DATE)) # check this worked
## Check that all years are 2023 and change if not
Sub_23 <- Sub_23 |> 
  mutate(year = year(VISIT_DATE),
         VISIT_DATE = ifelse(!year == 2023, paste0("2023", substr(VISIT_DATE, 5, 19)), paste0(VISIT_DATE)),
         VISIT_DATE = ymd_hms(VISIT_DATE)) |> 
  select(-year)
table(is.na(Sub_23$VISIT_DATE)); table(year(Sub_23$VISIT_DATE)) # check this worked

## bind the two halves of the data back together
Bird_All <- rbind(Sub21_22, Sub_23)




##--------------------------------##
#### 3. Combine Bird/Field Data ####
##--------------------------------##

##--------------------##
#### 3.1 Merge Data ####
##--------------------##

## filter just the wader records
Wader_All <- Bird_All |> filter(ENGLISH_NAME %in% c("Lapwing", "Redshank", "Curlew", "Snipe", "Oystercatcher"))


## Add on year column to wader data and fix visit number column
Wader_All <- Wader_All |> 
  mutate(year = year(VISIT_DATE)) 
## Pick columns from the bird survey data that I want to keep for the join
colnames(Wader_All)
Wader_All <- Wader_All |> 
  select(c(F_LOC_ID, S_LOC_ID, VISIT_NUMBER, year, #columns that match field data, used for the join
           ENGLISH_NAME, TOTAL_ADULTS, YOUNG_PRESENT, # columns that don't match and want to keep 
           TOTAL_YOUNG, ESTIMATED_PAIRS, ACTIVITY_CODES, Origin))  # columns that don't match and want to keep 


## Remove a couple of columns form the Field data that are not needed
Field_All2 <- Field_All |> 
    select(-c(SUB_ID, SITE_SUB_ID)) |>  
    mutate(year = year(VISIT_DATE))

## Join bird and field survey data together
Comb_All <- full_join(Wader_All, Field_All2, by = c("F_LOC_ID", "S_LOC_ID", "VISIT_NUMBER", "year"))
nrow(Wader_All)
nrow(Field_All)
nrow(Comb_All)

## Check join has worked properly
## The Field_All data set will be extended by the number of repeats in Wader_All due to multiple wader species being recorded in same field on same visit
DupWad <- Wader_All |> select(F_LOC_ID, S_LOC_ID, VISIT_NUMBER, year) %>% duplicated()
sum(DupWad) + nrow(Field_All) ## should equal length of Comb_All but does not yet
nrow(Comb_All)

## Check which (if any) waders records do not have survey data
set1 <- Comb_All |> filter(is.na(ENGLISH_NAME) == F & is.na(COVERAGE) == T)
nrow(set1) # no records now



##--------------------------------------------------##
#### 3.2 Remove fields surveyed in multiple years ####
##--------------------------------------------------##

## Some fields surveyed in 2021/22 were re-surveyed in 2023
## Remove these earlier surveys and keep the later survey (generally 2023)

## Summarise number of visits per field in each year
Field_Sum <- Comb_All |> 
  mutate(year = year(VISIT_DATE)) |> 
  group_by(F_LOC_ID, year) |> 
  summarise(n = n())

## Work out which fields were surveyed in multiple years (i.e. a repeat survey was carried out)
Ydups <- Field_Sum |> select(F_LOC_ID) %>% duplicated()
sum(Ydups) # check number of duplicates

## filter out all fields that are part of a repeat survey
All_dups <- Field_Sum |> filter(F_LOC_ID %in% unique(Field_Sum$F_LOC_ID[Ydups]))

## For each field with multiple survey years calculate the latest year of survey
Field_max <- Field_Sum |> 
  filter(F_LOC_ID %in% Field_Sum$F_LOC_ID[Ydups]) |> 
  group_by(F_LOC_ID) |> 
  summarise(year = max(year)) |> 
  mutate(keep = "Y") |> 
  full_join(All_dups)  

## Filter out the field years that need to be removed from main data set 
## i.e. those fields with repeat surveys from earlier years 
RepeatField <- Field_max[is.na(Field_max$keep)==T,]

## Remove the repeat field years from the main data set
Comb_All <- Comb_All |> mutate(year = year(VISIT_DATE))
Comb_All2 <- anti_join(Comb_All, RepeatField, by = c("F_LOC_ID", "year"))
nrow(Comb_All) - nrow(Comb_All2) ## how many rows lost; 2446
## Calculate the total number of unique fields were data was removed
losses <- inner_join(Comb_All, RepeatField, by = c("F_LOC_ID", "year"))
length(unique(losses$F_LOC_ID)) # 898 fields



##----------------------------------##
#### 3.3 Extract Somerset Repeats ####
##----------------------------------##

## For Somerset want to keep the earlier survey for fields that were surveyed in multiple years
## Will use these for the Lapwing Pair estimate calculation

## Read in RSPB priority landscapes and extract Somerset landscapes
Lanscap <- st_read(here("RawData", "Priority Landscapes", "EnglandWales_PriorityLandscapes.shp"))
Lanscap <- filter(Lanscap, Att3Value %in% c("Somerset Levels and Moors"))

## Read in the ESAs and extract Somerset ESA
ESAs <- st_read(here("RawData", "ESAs", "Wader_focused_ESAs.shp"))
ESAs <- st_combine(filter(ESAs, NAME == "SOMERSET LEVELS AND MOORS"))

## Combine priority landscapes/ESA and buffer
Somerset_buf <- st_union(Lanscap, ESAs)
Somerset_buf <-st_buffer(Somerset_buf, dist = 5000)
plot(Somerset_buf$geometry)


## For the resurveyed fields add on the sf field shapes so I can work out if they are in Somserset
FieldShapes2 <- FieldShapes |> select(F_LOC_ID)
RepeatField <- left_join(RepeatField, FieldShapes2, by = "F_LOC_ID")


## Extract the resurveyed fields fields that fall within the buffered Somerset landscapes
Somerset_buf <- Somerset_buf |> select(Att3Value)
Overlap <- st_intersection(Somerset_buf, st_as_sf(RepeatField))
rm(FieldShapes2, Somerset_buf, ESAs, Lanscap); gc()

## Get the full field data for each of the Somerset repeats
Overlap <- Overlap |> st_drop_geometry() |> select(F_LOC_ID, year)
Som_Repeats <- inner_join(Comb_All, Overlap, by = c("F_LOC_ID", "year"))




##----------------------------------##
#### 4. Save Cleaned Survey Data ####
##---------------------------------##

## Check whether any bird records are missing a percentage covered
MissCov <- Comb_All2 |> filter(is.na(ENGLISH_NAME)==F & is.na(COVERED_PERCENT)==T)
nrow(MissCov)

## Just before reading out cleaned survey data need to correct the percentage covered column
## If Coverage is 0 and there are birds going to presume this should be 100
Comb_All2 <- Comb_All2 |> 
  mutate(COVERED_PERCENT = ifelse(is.na(ENGLISH_NAME)==F & COVERED_PERCENT == 0, 100, COVERED_PERCENT))

## Now rite out this file
## Later can use this to work out Field level variables for analysis, e.g. hedge coverage, standing water, cattle presence
write_csv(x = Comb_All2, file = here("CleanData", "Wader Abundance", "1-CleaningScript", "Survey_Data_LONG_Clean.csv"))




##----------------------------------##
#### 5. Breeding Pair: Estimation ####
##----------------------------------##

## Add on Somerset Data from repeated years, this will be used to calculate Lapwing abundance only
## First create column that indicated if field is part of repeat survey and whether it is the earlier (1) or later survey (2)
Som_Repeats <- Som_Repeats |> mutate(Somerset = 1)
Comb_All2 <- Comb_All2 |> mutate(Somerset = ifelse(F_LOC_ID %in% unique(Som_Repeats$F_LOC_ID), 2, 0))
Comb_All2 <- rbind(Comb_All2, Som_Repeats)

## Some times have been entered incorrectly, e.g. first visits in April that were in the evening as these should all be morning visits
Comb_All2 <- Comb_All2 |> 
  mutate(VISIT_DATE = fifelse(hour(VISIT_DATE) > 17 & VISIT_NUMBER == 1, (ymd_hms(VISIT_DATE) - hours(12)), VISIT_DATE))


## Summary of how these calculation are currently carried out:
## Step 1. Remove fields with 0% coverage
## Step 2. Calculate corrected number of pairs and adults based on proportion of the field surveyed
## Step 3. Filter data for appropriate species time period
## Step 4. Calculate number of pairs for particular species 
## Step 5. Combine pairs counts with full list of fields so that have all the NAs. 
##    Overall gives the following counts for different situations:
##    - Birds present + survey appropriate period --> Species Calculation
##    - Birds not present + survey appropriate period --> 0
##    - No survey in appropriate period --> NA
##    - Never bird detected + No survey in appropriate period --> NA/0
##    - Single survey + Never bird detected --> NA/0
## Step 6. To pair counts add on Other Habitat, if birds ever observed in field, number of visits with coverage, earliest survey, latest survey 

## Change the name of the data set for this section
Pairs <- Comb_All2 


## Step 1** Remove fields with 0% Coverage
Pairs <- filter(Pairs, !COVERED_PERCENT == 0)
table(Pairs$COVERED_PERCENT) # check the values now


## Step 2** Calculate corrected number of pairs and adults based on proportion of the field surveyed
## NOTE: if less then 50% covered then cap the increase at x2

## Change "NA"s --> 0 in TOTAL_ADULTS/ESTIMATED_PAIRS/TOTAL_YOUNG columns
## These should all be true zeros now columns with 0% coverage are removed
Pairs <- Pairs |> 
  mutate(TOTAL_YOUNG = replace(TOTAL_YOUNG, is.na(TOTAL_YOUNG), 0),
         TOTAL_ADULTS = replace(TOTAL_ADULTS, is.na(TOTAL_ADULTS), 0),
         ESTIMATED_PAIRS = replace(ESTIMATED_PAIRS, is.na(ESTIMATED_PAIRS), 0))

## Calculate the corrected total adults and pairs
Pairs <- Pairs |> 
  mutate(across(c(COVERED_PERCENT, TOTAL_ADULTS, ESTIMATED_PAIRS), as.numeric),
         PropVisCov = COVERED_PERCENT/100,
         PropVisCov = ifelse(PropVisCov < 0.5, 0.5, PropVisCov),
         CORRECT_ADULTS = ifelse(PropVisCov == 0, 0, TOTAL_ADULTS/PropVisCov), 
         CORRECT_PAIRS = ifelse(PropVisCov == 0, 0, ESTIMATED_PAIRS/PropVisCov),
         yday =yday(VISIT_DATE))

## How many fields have less than 100% coverage
Tabs <- table(Pairs$PropVisCov<0.99)
(Tabs[2]/sum(Tabs))*100 # percentage fo fields with less than 100% coverage


## Summary of species-level conversions
# **Lapwing** - find visit at site-level with highest lapwing count between 15 April and 31 May. From this visit, half the number of individuals in each field.
# **Redshank - mean count of individuals in each field between 15 April – 20 May (O'Brien & Smith 1992 Bird Study)
# **Snipe** - EITHER the max number of pairs recorded in ESTIMATED_PAIRS on evening visits 1 May to 24 June
#           - OR if Snipe were seen in May or June, the number of individuals/2 (TOTAL_ADULTS) (whichever was higher)
  # => it wasn't recorded whether adults counted were chipping/drumming or just flushed
  # => SO assumed snipe recorded on evening visits were likely displaying males (classed as later than 18:00 - majority of the 4th visits for snipe) 
  # => those recorded at any other time could have been displaying or flushed
# **Oystercatcher** - half the maximum number of individuals on any visits between 15 April and 31 May
# **Curlew** - mean count of individuals in each field in all visits multiplied by 0.71 and then add 0.1 (optional) ((mean number of individuals x 0.71) + 0.1)

## Create a full list of all fields (with Year) with at least one survey with some level of coverage
## Needed to add fields that will receive NAs for pair estimates as they never had a survey inappropriate period
FieldList <- Pairs |> group_by(F_LOC_ID, year) |> summarise() |> ungroup()

## Before calculations make one data set with earlier repeat surveys and another with later repeat surveys
Pairs_E <- filter(Pairs, Somerset %in% c(0, 1))
Pairs_L <- filter(Pairs, Somerset %in% c(0, 2))


##-----------------------##
#### 5.1 Lapwing Pairs ####
##-----------------------##

## Filter for the correct date range and create lapwing adults column
Pairs_Lap <- Pairs_E |> 
  filter(yday >= 105 & yday <=  151) |> 
  mutate(CORRECT_ADULTS_Lap = ifelse(ENGLISH_NAME == "Lapwing", CORRECT_ADULTS, 0),
         CORRECT_ADULTS_Lap = ifelse(is.na(CORRECT_ADULTS_Lap) == T, 0, CORRECT_ADULTS_Lap))

## Find the site visit with the highest number of lapwings
## For ties choose the site visit where the highest number of fields were covered
SiteLap <- Pairs_Lap |> 
  group_by(S_LOC_ID, year, VISIT_NUMBER) |> 
  summarise(TOT_LAPWING = sum(CORRECT_ADULTS_Lap),
            F_LOCs = n()) |> 
  ungroup() |> group_by(S_LOC_ID) |> 
  filter(TOT_LAPWING == max(TOT_LAPWING)) |> 
  filter(F_LOCs == max(F_LOCs))

## If any sites are tied on total lapwing and number of field visited then remove later visits
SiteLap <- SiteLap[!duplicated(SiteLap$S_LOC_ID), ] |> select(-TOT_LAPWING)

## Use inner join to filter out survey data for the site visits with highest total lapwing number
FiltPairs_Lap <- inner_join(Pairs_Lap, SiteLap, by = c("S_LOC_ID", "year", "VISIT_NUMBER")) 
length(unique(Pairs_Lap$F_LOC_ID)); length(unique(FiltPairs_Lap$F_LOC_ID)) # Check how many fields lost

## Calculate estimated pairs as half number of adults in each field
Lap <- FiltPairs_Lap |> 
  ungroup() |> group_by(F_LOC_ID, year) |> 
  summarise(est_pairsL = max(CORRECT_ADULTS_Lap, na.rm=T)) |> # take max number of adults per field
  mutate(est_pairsL = ifelse(est_pairsL == 0, 0, est_pairsL/2)) # halve none zero adult estimates

## Add on other fields not in lapwing survey period to get fields that should be NAs
nrow(FieldList) - nrow(Lap) # number of fields without survey in appropriate period
Lap <- full_join(FieldList, Lap, by = c("F_LOC_ID", "year"))


##------------------------##
#### 5.2 Redshank Pairs ####
##------------------------##

## Calculate the mean number of Redshank per field across all visits
## Need to include visits to fields where no Redshank where recorded (including those with no birds recorded)
Red <- Pairs_L |> 
  filter(yday >= 105 & yday <=  140 & hour(VISIT_DATE) < 17) |> # select visits in correct time period during the day
  mutate(CORRECT_ADULTS_Red = ifelse(!ENGLISH_NAME == "Redshank" | is.na(ENGLISH_NAME) ==T, 0, CORRECT_ADULTS)) |> # create a column for number of corrected redshank adults
  group_by(F_LOC_ID, year, VISIT_NUMBER) |> 
  summarise(est_pairsR = max(CORRECT_ADULTS_Red, na.rm=T)) |> # get the max Redshank count for each field visit
  ungroup() |> group_by(F_LOC_ID, year) |> 
  summarise(est_pairsR = mean(est_pairsR, na.rm=T)) # get the mean number of Redshank across all field visits

## Add on other Fields not within the redhsnk survey period  as they should be NAs
nrow(FieldList) - nrow(Red) # number of fields without survey in appropriate period
Red <- full_join(FieldList, Red, by = c("F_LOC_ID", "year"))


##---------------------##
#### 5.3 Snipe Pairs ####
##---------------------##

## Filter for the correct date range and create Snipe adults column
## Apply correction to number of pairs here by halving counts during the day
Pairs_Sni <- Pairs_L  |> 
  filter(yday >= 121 & yday <=  175) |> 
  mutate(Eve = ifelse(hour(VISIT_DATE) >= 18, "Evening", "Day"),
         CORRECT_ADULTS_Sni = ifelse(!ENGLISH_NAME == "Snipe" | is.na(ENGLISH_NAME) == T, 0, CORRECT_ADULTS),
         CORRECT_ADULTS_Sni = ifelse(Eve == "Day", CORRECT_ADULTS_Sni/2, CORRECT_ADULTS_Sni)) 

## Find the site visit with the highest number of Snipe (this could be day or night)
## For ties choose the site visit where the highest number of fields were covered
SiteSni <- Pairs_Sni |> 
  group_by(S_LOC_ID, year, VISIT_NUMBER, Eve) |> 
  summarise(TOT_SNIPE = sum(CORRECT_ADULTS_Sni),
            F_LOCs = n()) |> 
  ungroup() |> group_by(S_LOC_ID) |> 
  filter(TOT_SNIPE == max(TOT_SNIPE)) |> 
  filter(F_LOCs == max(F_LOCs))

## If any sites are still tied on total snipe and number of fields visited then remove later visits
SiteSni <- SiteSni[!duplicated(SiteSni$S_LOC_ID), ] |> select(-TOT_SNIPE)

## Use inner join to filter out survey data from site visits with highest total snipe number
FiltPairs_Sni <- inner_join(Pairs_Sni, SiteSni, by = c("S_LOC_ID", "year", "VISIT_NUMBER", "Eve")) 
length(unique(Pairs_Sni$F_LOC_ID)); length(unique(FiltPairs_Sni$F_LOC_ID)) # Check how many fields lost

## Calculate number of corrected snipe per field, halving adults if the survey was from the day
## NOTE: should be 414 repeats due to some repeat fields for Somerset Survey
Sni <- FiltPairs_Sni |> 
  ungroup() |> group_by(F_LOC_ID, year) |> 
  summarise(est_pairsS = max(CORRECT_ADULTS_Sni, na.rm=T)) ## take max number of corrected snip for each field visit

## Add on other Fields not in Snipe survey period to get fields that should be NAs
nrow(FieldList) - nrow(Sni) # number of fields without survey in appropriate period
Sni <- full_join(FieldList, Sni, by = c("F_LOC_ID", "year"))


##-------------------##
#### 5.4 Oyc Pairs ####
##-------------------##

## Calculate the max number of Oystercatchers in each field across all visits 
Oyc <- Pairs_L |> 
  filter(yday >= 105 & yday <=  151) |> 
  mutate(CORRECT_ADULTS_Oyc = ifelse(!ENGLISH_NAME == "Curlew" | is.na(ENGLISH_NAME) ==T, 0, CORRECT_ADULTS),
         CORRECT_ADULTS_Oyc = ifelse(CORRECT_ADULTS_Oyc >= 8, 0, CORRECT_ADULTS_Oyc)) |> # change counts of flocks (deemed to counts greater than 8) to zero
  group_by(F_LOC_ID, year) |> 
  summarise(est_pairsO = max(CORRECT_ADULTS_Oyc, na.rm=T)) |> 
  mutate(est_pairsO = ifelse(est_pairsO == 0, 0, est_pairsO/2))

## Add on other Fields not within period to get fields that should be NAs
nrow(FieldList) - nrow(Oyc) # number of fields without survey in appropriate period
Oyc <- full_join(FieldList, Oyc, by = c("F_LOC_ID", "year"))


##----------------------##
#### 5.5 Curlew Pairs ####
##----------------------##

## Calculate the mean number of Curlew pairs in each field across all visits (then multiply by 0.71)
Cur <- Pairs_L |> 
  filter(hour(VISIT_DATE) < 17) |> # remove night visits
  mutate(CORRECT_ADULTS_Cur = ifelse(!ENGLISH_NAME == "Curlew" | is.na(ENGLISH_NAME) ==T, 0, CORRECT_ADULTS), # create number of curlew column
         CORRECT_ADULTS_Cur = ifelse(CORRECT_ADULTS_Cur > 8, 0, CORRECT_ADULTS_Cur)) |> # change counts of flocks (deemed to be counts greater than 8) to zero
  group_by(S_LOC_ID, F_LOC_ID, year, VISIT_NUMBER) |> 
  summarise(est_pairsC = max(CORRECT_ADULTS_Cur, na.rm=T)) |> # get the max Curlew count for each field visit
  ungroup() |> group_by(S_LOC_ID, F_LOC_ID, year) |> 
  summarise(est_pairsC = mean(est_pairsC, na.rm=T)) |> # calculate mean number of Curlew across all field visits
  mutate(est_pairsC = ifelse(est_pairsC == 0, 0, (est_pairsC*0.71))) # multiple mean count by 0.71

## Need to add on 0.1 to each site
## Do this by dividing 0.1 by the number of occupied fields in a site and adding this amount on to occupied fields
Cur <- Cur |>  
  filter(est_pairsC > 0) |> 
  ungroup() |> group_by(S_LOC_ID) |>
  summarise(Occ_fields = n()) |> 
  full_join(Cur, by = "S_LOC_ID") |> 
  mutate(est_pairsC = fifelse(est_pairsC > 0, (est_pairsC + (0.1/Occ_fields)), est_pairsC)) |> 
  select(-c(Occ_fields, S_LOC_ID))


## Add on other Fields not within period (Should be 0 for Curlews)
nrow(FieldList) - nrow(Cur) # number of fields without survey in appropriate period
Cur <- full_join(FieldList, Cur, by = c("F_LOC_ID", "year"))



##-----------------------##
#### 5.6 Combine Pairs ####
##-----------------------##

## Join together all pair estimates
EstCounts <- full_join(Lap, Oyc) |> full_join(Cur) |> full_join(Red) |> full_join(Sni) 
table(duplicated(EstCounts$F_LOC_ID)) # double check duplicated fields (should be 414 due to Somerset repeats as 2021 fields for Lapwing)


## Make a summary of the field data that I can add onto the breeding pair estimates
## Add on Other Habitat, if birds ever observed in field, any night surveys
FieldSum <- Pairs |> 
            mutate(Unsuit = ifelse(OTHER_HABITAT =="U", 1, 0),
                   Arab = ifelse(OTHER_HABITAT == "A", 1, 0),
                   SuitOther = ifelse(OTHER_HABITAT == "S", 1, 0),
                   Eve = ifelse(hour(VISIT_DATE) >= 18, 1, 0)) |> 
            group_by(F_LOC_ID, year, VISIT_DATE) |> 
            summarise(Unsuit = max(Unsuit),
                      Arab = max(Arab),
                      SuitOther = max(SuitOther),
                      Eve = max(Eve),
                      Somerset = max(Somerset)) |> 
            ungroup() |> group_by(F_LOC_ID, year) |> 
            summarise(N_visits = n(),
                      Prop_unsuit = sum(Unsuit, na.rm = T)/N_visits,
                      Prop_arable = sum(Arab, na.rm = T)/N_visits,
                      Prop_suitother = sum(SuitOther, na.rm = T)/N_visits,
                      N_Eve_visits = sum(Eve),
                      Max_date = max(VISIT_DATE),
                      Min_date = min(VISIT_DATE),
                      Unsuit_Hab = ifelse(Prop_unsuit < 1, 0, 1),
                      Arable_Hab = ifelse(Prop_arable < 1, 0, 1),
                      Other_Hab = ifelse(Prop_suitother < 1, 0, 1),
                      Somerset = max(Somerset)) |> 
           ungroup() |> select(-c(Prop_unsuit, Prop_arable, Prop_suitother)) |> 
           mutate(Somerset = ifelse(Somerset == 0, "All", ifelse(Somerset == 1, "L", "S")))


## Now bind this onto the counts
EstCounts_Final <- full_join(EstCounts, FieldSum, by = c("F_LOC_ID", "year"))
nrow(EstCounts_Final) == nrow(EstCounts) # check join 

## Since for some fields we only examined if there were Lapwing the whole row could be NA
## So remove rows where all wader estimates are NA
EstCounts_Final <- EstCounts_Final |>  
  filter(is.na(est_pairsL) ==F | is.na(est_pairsO) ==F | is.na(est_pairsC) ==F | 
         is.na(est_pairsR) ==F | is.na(est_pairsS) ==F)





##-------------------------------##
#### 6. Breeding Pair: Summary ####
##-------------------------------##

## Summarise the number of individual fields with wader records and the total number of pairs in those fields

## Lapwing summary
LapSum <- EstCounts_Final |> 
  filter(Somerset %in% c("All", "L")) |> 
  filter(est_pairsL > 0) |> 
  group_by(F_LOC_ID, year) |> 
  summarise(est = max(est_pairsL)) |> 
  ungroup() |> group_by(year) |> 
  summarise(fields_occ = n(),
            total_pairs = sum(est))

## Redshank summary
RedSum <- EstCounts_Final |> 
  filter(Somerset %in% c("All", "S")) |> 
  filter(est_pairsR > 0) |> 
  group_by(F_LOC_ID, year) |> 
  summarise(est = max(est_pairsR)) |> 
  ungroup() |> group_by(year) |> 
  summarise(fields_occ = n(),
            total_pairs = sum(est))

## Snipe summary
SniSum <- EstCounts_Final |> 
  filter(Somerset %in% c("All", "S")) |> 
  filter(est_pairsS > 0) |> 
  group_by(F_LOC_ID, year) |> 
  summarise(est = max(est_pairsS)) |> 
  ungroup() |> group_by(year) |> 
  summarise(fields_occ = n(),
            total_pairs = sum(est))




##----------------------------##
#### 7. Breeding Pair: Save ####
##----------------------------##

## Write out this full clean data set
write_csv(x = EstCounts_Final, file = here("CleanData", "Wader Abundance", "1-CleaningScript", "Estimated_Pairs_Clean.csv"))
colnames(EstCounts_Final)



##-------------------------------##
#### 8. Field Data: Amalgamate ####
##-------------------------------##

## Change the name of the data set for this section
Fields <- Comb_All2 

## Remove fields with 0% Coverage
# Fields <- filter(Fields, !COVERED_PERCENT == 0)
# table(Fields$COVERED_PERCENT) # check the values now

## Summaries the non-bird data for each field
## What is the "DITCH_WATER" column

## Summary of each field
colnames(Fields)
FieldSum <- Fields |> 
            mutate(Stock = fifelse((GRZ_CATTLE == "Y" | GRZ_SHEEP == "Y" | GRZ_OTHER == "Y"), 1, 0),
                   Stock = fifelse(is.na(Stock) == T, 0, Stock),
                   OTHER_HABITAT = ifelse(is.na(OTHER_HABITAT) == T & is.na(GRASSLAND_TYPE) == F, "G", OTHER_HABITAT)) |> 
            arrange(VISIT_NUMBER) |> # order by visit number, ensures I can extract standing water from first visit
            group_by(F_LOC_ID, S_LOC_ID, year) |> 
            summarise(GRASSLAND_TYPE = getmode(GRASSLAND_TYPE),
                      OTHER_HABITAT = getmode(OTHER_HABITAT),
                      TALL_BOUNDARY_PERCENT = getmode(TALL_BOUNDARY_PERCENT),
                      STANDING_WATER_TOTAL_PERCENT_vis1 = head(STANDING_WATER_TOTAL_PERCENT, n=1L, na.rm = T), # standing water from first visit
                      STANDING_WATER_TOTAL_PERCENT = mean(STANDING_WATER_TOTAL_PERCENT, na.rm = T),
                      GROUND_DAMP = getmode(GROUND_DAMP),
                      STOCK = max(Stock, na.rm = T), 
                      VEG_HEIGHT = getmode(VEG_HEIGHT), 
                      VEG_STRUCTURE = getmode(VEG_STRUCTURE), 
                      RUSH_PERCENT = mean(RUSH_PERCENT, na.rm = T), 
                      RUSH_DISTRIBUTION = getmode(RUSH_DISTRIBUTION)) |> 
            mutate(RUSH_PERCENT = ifelse(RUSH_PERCENT == "NaN", NA, RUSH_PERCENT),
                   GROUND_DAMP = fifelse(STANDING_WATER_TOTAL_PERCENT > 0 & GROUND_DAMP == "N", "Y", GROUND_DAMP))
            



##-------------------------##
#### 9. Field Data: Save ####
##-------------------------##

## Write out the non-bird data for each field
write_csv(x = FieldSum, file = here("CleanData", "Wader Abundance", "1-CleaningScript", "Field_Characteristics_Clean.csv"))
colnames(FieldSum)












