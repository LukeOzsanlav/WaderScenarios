##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 16/07/2024
## 
## 
## Aim: 
## 1. Calculate average size of an AES agreement and RSPB Reserve
## 2. Calculate the uptake of different wader associated options by region
## 
## 
## 
##------------------------------------------------------##



## Load in packages & helper functions
pacman::p_load(sf, tidyverse) #  terra
options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R")


##------------------------------##
#### 1.0 Average Reserve Size ####
##------------------------------##

## Read in the shapefile of RSPB reserves
RSPB <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")

## Read in priority landscapes
PrLand <- st_read("RawData/Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
PrLand <- filter(PrLand, Att3Value %in% c("Somerset Levels and Moors", "Broads", "Greater Thames"))

## Select just the reserves that fall within the priority landscapes of interest
sub <- st_intersection(RSPB, PrLand)
PrRes <- filter(RSPB, Att3Value %in% sub$Att3Value)


## calculate the area and and standard deviation of reserve sizes
AllSizes <- data.frame(ReserveAve = mean(as.numeric(st_area(PrRes))/10000),
                       ReserveSD = sd(as.numeric(st_area(PrRes))/10000))




##------------------------------------##
#### 2.0 Average AES agreement size ####
##------------------------------------##

## Read in the RPA point parcels for the Somerset regions
PPoints1 <- st_read("RawData/RPA/rpa_parcel_points_Somerset/rpa_parcel_pointsPoint.shp")
PPoints2 <- st_read("RawData/RPA/rpa_parcel_points_East_England/rpa_parcel_pointsPoint.shp")
PPoints <- rbind(PPoints1, PPoints2)
rm(PPoints1, PPoints2); gc()

## Read in the RPA point parcels for the Somerset regions
PointsID <- read_csv("RawData/RPA/RPA_ParcelList_20240424&Suffolk.csv")

## Join together the field parcel points and the customer data from the RPA
Join <- inner_join(PointsID, PPoints, by = join_by(NGC==parcel_ref)) |> st_as_sf() |> select(FARM_ID)
plot(Join$geometry)

## Read in the UKCEH land parcels shapefile
CEHParcels <- st_read("RawData/LandCover/UKCEH_LCvector_2021.shp") |> select(geometry)


# Perform spatial join to find where points fall within polygons
st_crs(CEHParcels)==st_crs(Join) # check crs are the same
join_result <- st_join(CEHParcels, Join, join = st_contains)
join_result <- join_result |> filter(is.na(FARM_ID)==F) # remove fields without a parcel point in them

## remove any duplicate shapes, occurs if there are multiple points in a polygon
dupind <- duplicated(join_result$geometry)
join_result <- join_result[dupind==F, ]

# Group by polygon ID and count the number of points in each
AESSize <- join_result |> 
           group_by(FARM_ID) |> 
           summarise(Size = sum(as.numeric(st_area(geometry)))/10000)

## Now calculate the average and sd of AES agreement size and assign these to the data frame
AllSizes <- AllSizes |>  mutate(AESAve = mean(AESSize$Size),
                                AESsd = sd(AESSize$Size))


## save the final reserve and AES sizes
write_csv(AllSizes, "CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")




##------------------------------##
#### 3.0 AES uptake frequency ####
##------------------------------##

##-- CSS Options --##

## Read in the different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")

## Filter out breeding/wintering wader focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))

## Filter out all the land parcels that have an option for breeding or wintering waders
AllWaderFields <- filter(CSS_Ops, parcref %in% CSS_OpsWa$parcref)


## Read in a shape file of the priority landscapes
Landsc <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |> select(Att3Value)


## Use st_join to extract and label the field points in different priority landscapes
Landsc_Fields <- st_join(AllWaderFields,
                         Landsc,
                         left = F) # not left = inner join

## filter out the points from my landscapes of interest
Landsc_Fields <- filter(Landsc_Fields, Att3Value %in% c("Broads", "Somerset Levels and Moors", "Essex", "North Kent"))
table(Landsc_Fields$Att3Value)

## Create data frame of the different options by region
TableCSS_WaderManage <- as.data.frame((table(Landsc_Fields$opt_desc, Landsc_Fields$Att3Value))) |>  filter(Freq>0)
CSS_Trial <- as.data.frame((table(Landsc_Fields$opt_desc, Landsc_Fields$opt_code))) |>  filter(Freq>0)


AllWader <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp")
## Read in filtered breeding pairs estimates
## And calculate total wader abundance for each field
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
Waders <- Waders |> mutate(Tot_abund = rowSums(dplyr::select(Waders, est_pairsL, est_pairsR, est_pairsS), na.rm = T)) |> 
                    select(F_LOC_ID, Tot_abund) |> group_by(F_LOC_ID) |>
                    summarise(Tot_abund = sum(Tot_abund, na.rm = T))


## Join the Wader Cluster ID code onto the main canvas
WaderAbund <- left_join(AllWader, Waders, by = "F_LOC_ID")

WaderAbund <- filter(WaderAbund, Tot_abund > 0 & is.na(Region)==F) |> select(Tot_abund, Region)

## Use st_join to extract and label the field points in different priority landscapes
Landsc_Fields2 <- st_join(AllWaderFields,
                         WaderAbund,
                         left = F) # not left = inner join

## Create data frame of the different options by region
TableCSS_WaderOcc <- as.data.frame((table(Landsc_Fields2$opt_desc, Landsc_Fields2$Region))) |>  filter(Freq>0)



##-- ESS Options --##

## Read in the ESS options
ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESS <- rbind(ESSWest, ESSEast)
rm(ESSWest, ESSEast)



## Filter out ESS wader focused options
ESS_OpsWa <- filter(ESS, optcode %in% c("HK9", "HK10", "HK11", "HK12", "HK13", "HK14"))
## Filter out all the land parcels that have an option for breeding or wintering waders
ESSWaderFields <- filter(ESS, parcref %in% ESS_OpsWa$parcref)

## Read in a shape file of the priority landscapes
Landsc <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |> select(Att3Value)

## Use st_join to extract and label the field points in different priority landscapes
Landsc_Fields <- st_join(ESSWaderFields,
                         Landsc,
                         left = F) # not left = inner join

## filter out the points from my landscapes of interest
Landsc_Fields <- filter(Landsc_Fields, Att3Value %in% c("Broads", "Somerset Levels and Moors", "Essex", "North Kent"))
table(Landsc_Fields$Att3Value)

## Create data frame of the different options by region
TableESS_WaderManage <- as.data.frame((table(Landsc_Fields$opttitle, Landsc_Fields$optcode, Landsc_Fields$Att3Value))) |>  filter(Freq>0)
Trial2 <- as.data.frame((table(Landsc_Fields$opttitle, Landsc_Fields$optcode))) |>  filter(Freq>0)


## Red in a shapefile of fields that are in my wader clusters
AllWader <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp")
## Read in filtered breeding pairs estimates
## And calculate total wader abundance for each field
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
Waders <- Waders |> mutate(Tot_abund = rowSums(dplyr::select(Waders, est_pairsL, est_pairsR, est_pairsS), na.rm = T)) |> 
                    select(F_LOC_ID, Tot_abund) |> group_by(F_LOC_ID) |>
                    summarise(Tot_abund = sum(Tot_abund, na.rm = T))


## Join the Wader Cluster ID code onto the field polygons
WaderAbund <- left_join(AllWader, Waders, by = "F_LOC_ID")
## Filter out fields that have at least a breeding pair of waders
WaderAbund <- filter(WaderAbund, Tot_abund > 0 & is.na(Region)==F) |> select(Tot_abund, Region)

## Use st_join to extract and label the field points in different priority landscapes
Landsc_Fields2 <- st_join(ESSWaderFields,
                          WaderAbund,
                          left = F) # not left = inner join

## Create data frame of the different options by region
TableESS_WaderOcc <- as.data.frame((table(Landsc_Fields2$opttitle, Landsc_Fields2$Region))) |>  filter(Freq>0)







##----------------------------------##
#### 4.0 Regional Wader AES Spend ####
##----------------------------------##

## create data frame to record costs
AESTable <- data.frame(Landscape = NA,
                       Cost = NA,
                       Area = NA,
                       Fields = NA,
                       Length = 1:4) 

##-------------------------##
## 4.1 North Kent AES Cost ##
##-------------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)
table(duplicated(Canvshp$ParcRef))

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")
table(duplicated(Canv$ParcRef))

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv,
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA,
               LapAbund = NA,
               RedAbund = NA) |>
  rename(FieldArea=ParcArea)


## Add a new column that separates out breeding wader vs wintering wader payment
Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})


##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")

## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
## This includes AES Only fields and reserve fields with AES
WhichAES <- which(is.na(Canv$AES_options) == F)

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){
  
  ## How many years to spread one off costs over
  Yrs <- 15

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
AESTable$Cost[1] <- sum(Canv$AESCostGDP[WhichAES], na.rm=T)
AESTable$Area[1] <- sum(Canv$FieldArea[WhichAES], na.rm=T)
AESTable$Fields[1] <- length(Canv$FieldArea[WhichAES])
AESTable$Landscape[1] <- "North Kent"



##-----------------------##
## 4.2 Somerset AES Cost ##
##-----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)


##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
## This includes AES Only fields and reserve fields with AES
WhichAES <- which(is.na(Canv$AES_options) == F)

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
AESTable$Cost[2] <- as.numeric(sum(Canv$AESCostGDP[WhichAES], na.rm=T))
AESTable$Area[2] <- sum(Canv$FieldArea[WhichAES], na.rm=T)
AESTable$Fields[2] <- length(Canv$FieldArea[WhichAES])
AESTable$Landscape[2] <- "Somerset"



##--------------------------##
## 4.3 Essex Coast AES Cost ##
##--------------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)


##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
## This includes AES Only fields and reserve fields with AES
WhichAES <- which(is.na(Canv$AES_options) == F)

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
AESTable$Cost[3] <- as.numeric(sum(Canv$AESCostGDP[WhichAES], na.rm=T))
AESTable$Area[3] <- sum(Canv$FieldArea[WhichAES], na.rm=T)
AESTable$Fields[3] <- length(Canv$FieldArea[WhichAES])
AESTable$Landscape[3] <- "Essex"




##------------------------------##
## 4.1 Norfolk Broads AES Costs ##
##------------------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)



##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
## This includes AES Only fields and reserve fields with AES
WhichAES <- which(is.na(Canv$AES_options) == F)

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
AESTable$Cost[4] <- as.numeric(sum(Canv$AESCostGDP[WhichAES], na.rm=T))
AESTable$Area[4] <- sum(Canv$FieldArea[WhichAES], na.rm=T)
AESTable$Fields[4] <- length(Canv$FieldArea[WhichAES])
AESTable$Landscape[4] <- "Norfolk"




##-----------------------------##
## 4.5 Summarise the AES costs ##
##-----------------------------##

## see the regional results
print(AESTable)

## Calculate mean value
mean(AESTable$Cost)
mean(AESTable$Area)
mean(AESTable$Fields)

## 10%/20% and 50% expenditure
mean(AESTable$Cost)*0.1
mean(AESTable$Cost)*0.2
mean(AESTable$Cost)*0.5

## Finally write out the table so I can use it in other scripts
write_csv(AESTable, file = "CleanData/Scenarios/4-AnnotateCanvas/Total_AES_Expenditure.csv")





##-------------------------------------##
#### 5.0 Reserve level wader density ####
##-------------------------------------##

##-- North Kent --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)
table(duplicated(Canvshp$ParcRef))

## Read in annotated canvas as csv
CanvNK <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")
table(duplicated(CanvNK$ParcRef))

## Join polygons to annotations
CanvNK <- left_join(CanvNK, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
CanvNK <- mutate(CanvNK, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING"))) |> 
        rename(FieldArea=ParcArea) |> 
        select(ReserveGroup, Category, RSPB, FieldArea, Tot_abund, Landscape)



##-- Somerset Levels --##
 
## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
CanvSom <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.rds")

## Join polygons to annotations
CanvSom <- left_join(CanvSom, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
CanvSom <- mutate(CanvSom, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING"))) |> 
        rename(FieldArea=ParcArea) |> 
        select(ReserveGroup, Category, RSPB, FieldArea, Tot_abund, Landscape)
  
  
  
##-- Essex Coast --##  

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
CanvEs <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.rds")

## Join polygons to annotations
CanvEs <- left_join(CanvEs, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## Read in a data set of Wallasea ISland F_LOC_IDs that are not wet grassland but actually saltmarsh
Wall <- read.csv("RawData/RSPB Reserves/WallaseaIsland_NonWetgrass.csv")

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
CanvEs <- mutate(CanvEs, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               Category = ifelse(F_LOC_ID %in% Wall$F_LOC_ID, "NoOpp", Category),
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               ) |> 
        rename(FieldArea=ParcArea) |> 
        select(ReserveGroup, Category, RSPB, FieldArea, Tot_abund, Landscape)




##-- Norfolk Broads --##


## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
CanvBr <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.rds")

## Join polygons to annotations
CanvBr <- left_join(CanvBr, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
CanvBr <- mutate(CanvBr, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING"))) |> 
        rename(FieldArea=ParcArea)|> 
        select(ReserveGroup, Category, RSPB, FieldArea, Tot_abund, Landscape)


##-- Calculate breeding wader density --##

## Join all the regions together
AllCanv <- rbind(CanvNK, CanvSom, CanvEs, CanvBr)

## Calculate the wader density
WaderDens <- AllCanv |>
             filter(is.na(ReserveGroup)==F) |>
             filter(is.na(Tot_abund) == F) |> 
             filter(Category == "Reserve") |>
             group_by(ReserveGroup, RSPB, Landscape) |>
             summarise(TotArea = sum(FieldArea),
                       TotWader = sum(Tot_abund, na.rm=T),
                       WaderDensity = TotWader/TotArea) |> 
             mutate(Cutoff = ifelse(Landscape == "Somerset Levels and Moors", 0.15, 0.3),
                    ResQual = ifelse(WaderDensity > Cutoff, "HighQ", "LowQ"))


## finally write out the breeding wader densities
write_csv(WaderDens |> st_drop_geometry(), "CleanData/Scenarios/4-AnnotateCanvas/ReserveWaderDensity.csv")



