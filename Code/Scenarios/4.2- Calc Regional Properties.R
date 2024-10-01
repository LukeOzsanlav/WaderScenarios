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









