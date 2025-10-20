##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 21/12/23
## Goal: Add wider landscape attributes to survey fields
## 
## Breakdown:
## 1. Calculate percentage of standing water in semi-natural wet grassland
## 2. Calculate percentage of woodland
## 3. Calculate percentage of urban areas
## 4. Calculate percentage of intensive land use
## 5. Calculate percentage of semi-natural wet grassland
## 6. Calculate percentage SSSI coverage
## 7. Calculate percentage reserve coverage
## 8. Calculate percentage wader AES coverage
## 
##------------------------------------------------------##


## Load in required packages
pacman::p_load(here, tidyverse, data.table, sf, terra)
options(scipen = 100, digits = 4)

## Use `here` package for specifying file paths
here::i_am("Code/Wader Abundance/4 - Add Landscape Attributes.R")
library(here)
here()

## Global Setting: Buffer Distances to trial
BuffDists <- seq(from = 500, to = 2000, by = 500)






##---------------------##
#### 0. Data Read In ####
##---------------------##

## Read in filtered breeding pairs estimates
Waders <- read_csv(here("CleanData", "Wader Abundance", "3-AddFieldAttributes", "Breeding_Pairs_Attribs.csv")) |> select(-geometry)
# Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")

## read in shaepfile for BWWM fields and join to breeding waders pairs
BWWM <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp")) |> select(F_LOC_ID)
Waders <- left_join(Waders, BWWM, by = "F_LOC_ID")
SL <- nrow(Waders) # keep track of data set length throughout

## Read in the UKCEH landcover data
LC <- rast(here("RawData", "LandCover", "gblcm25m2021.tif"))
LC <- LC[[1]]

## Read in RSPB priority landscapes, with greater Thames landscape split in two
SplitLanscape <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |> 
  filter(Att3Value %in% c("Somerset Levels and Moors", "Broads", "Greater Thames", "North Kent", "Essex"))




##------------------------------------------------##
#### 0.2 Function: Landscape values from raster ####
##------------------------------------------------##

## Function that takes the wader survey data (as shape file of the fields) and a raster  
## Then calculates the proportion coverage of the raster in a buffer around the center of each field

## Survey = Survey data as a shape file of each field
## Rast = Raster of variable of interest, should only contain 0 and 1
## Dist = Buffer distances to trial
## Col = Start of the column name to use for the proportion cover column (distance used as suffix)
Prop_Raster_Cov <- function(Survey = Survey, Rast = Rast, Dist = Dist, Col = Col){
  
  ## crop the raster data to the extent of the survey data (helps with speed)
  Rastcrop <- crop(Rast, vect(st_buffer(st_as_sf(Survey), dist = 5000)))
  plot(Rastcrop)
  
  ## Get the centers for each field and then add the buffer
  Centers <- st_centroid(st_as_sf(Survey))
  
  ## loop through each distance
  for(j in 1:length(Dist)){
  
  ## Progress message  
  message("Buffer ", j, " out of ", length(Dist))  
  
  ## Add the buffer on to the field center
  CentBuffs <- st_buffer(Centers, dist = Dist[j])
     
  ## Extract values from raster to fields (weighted average)
  ExtractVals <- extract(Rastcrop, vect(CentBuffs), fun = mean, weights = T, na.rm = T)
  
  ## Add the buffer values to the main wader data set
  Survey <- mutate(Survey, VAL = ExtractVals[,2]) 
  colnames(Survey)[colnames(Survey) == "VAL"] <- paste0(Col, Dist[j])
  
  }
  return(Survey)
} # end of function




##---------------------------------------------------##
#### 0.3 Function: Landscape values from shapefile ####
##---------------------------------------------------##

## Function that takes the wader survey data (as shape file of the fields) and a shape file of a variable of interest
## Then calculates the proportion of each buffered field covered by the shapes of interest

## Survey = Survey data as a shape file of each field
## Shapes = Shape file of variable of interest
## Dist = Buffer distances to trial
## Col = Start of the column name to use for the proportion cover column (distance used as suffix)
Prop_Shapes_Cov <- function(Survey = Survey, Shapes = Shapes, Dist = Dist, Col = Col){
  
  ## crop the shapefile data to the extent of the survey data (helps with speed)
  Shapescrop <- st_crop(Shapes, st_buffer(st_as_sf(Survey), dist = 5000))
  plot(Shapescrop$geometry)
  
  ## Get the centers for each field and then add the buffer
  Centers <- st_centroid(st_as_sf(Survey))
  
  ## loop through each distance
  for(j in 1:length(Dist)){
    
  ## Progress message  
  message("Buffer ", j, " out of ", length(Dist))  
  
  ## Add the buffer on to the field center
  CentBuffs <- st_buffer(Centers, dist = Dist[j])

  ## Work out the area of the shapes that overlaps each buffered survey field
  Overlap_sub <- CentBuffs |> filter(Somerset %in% c("All", "S")) |> 
                 st_intersection(Shapescrop) |> # intersect shapes and survey fields
                 mutate(OverlapArea = st_area(geometry)) |> # calculate area of each intersected polygon
                 select(F_LOC_ID, OverlapArea) |> 
                 st_drop_geometry() |> 
                 group_by(F_LOC_ID) |> # group by field as buffered survey fields may overlap multiple shapes
                 summarise(Tot_Overlap = sum(OverlapArea))

  ## Join the overlap areas for each buffered survey field onto the main fields data
  ## & Work out the proportion of each buffered field overlapped by shapes of interest
  Survey <- left_join(Survey, Overlap_sub, by = "F_LOC_ID") |> 
             mutate(Area = st_area(CentBuffs$geometry),
                    PropOverlap = ifelse(is.na(Tot_Overlap)==T, 0, Tot_Overlap/Area)) |> 
             select(-c(Area, Tot_Overlap)) 
  
  ## change the column Name to one of your choosing
  colnames(Survey)[colnames(Survey)== "PropOverlap"] <- paste0(Col, Dist[j])
  }
  return(Survey)
} # end of function



##-----------------------------------##
#### 1. Add split landscape labels ####
##-----------------------------------## 

## Remove duplicated from the fields data set, this happens as we have different survey fields for Lapwing/Snipe in Somerset
Waderssub <- Waders[duplicated(Waders$F_LOC_ID)==F,]


## Work out if fields fall within any of the priority landscapes
Overlap <- st_intersection(SplitLanscape, st_as_sf(Waderssub))

## Join the priority landscape labels onto that data set with all the fields in
## first streamline the data set from the overlap
Overlap <- Overlap |> 
  select(F_LOC_ID, Att3Value) |> 
  rename(Landscape4 = Att3Value) |> 
  st_drop_geometry()


## Work out if fields fall within any of the buffered priority landscapes
## This is done int two steps here as the boundaries of North kent and Essex overlap when they are buffered
OverlapBuf <- st_intersection(st_buffer(SplitLanscape, dist = 5000), st_as_sf(Waderssub))

## Join the priority landscape labels onto that data set with all the fields in
## first streamline the data set from the overlap
OverlapBuf <- OverlapBuf |> 
  select(F_LOC_ID, Att3Value) |> 
  rename(Landscape4buf = Att3Value) |> 
  st_drop_geometry()


## Join the data sets together, and for row change in join
PreL <- nrow(Waders)
Waders <- left_join(Waders, Overlap, by = "F_LOC_ID")
nrow(Waders) == PreL # check no rows lost


## Join the data sets together, and for row change in join
PreL <- nrow(Waders)
OverlapBuf <- OverlapBuf[duplicated(OverlapBuf$F_LOC_ID)==F,] # remove duplicated from where Essex/North Kent boundaries overlap
Waders <- left_join(Waders, OverlapBuf, by = "F_LOC_ID")
nrow(Waders) == PreL # check no rows lost

## Mutate the 
Waders <- Waders |> 
  mutate(Landscape4 = ifelse(is.na(Landscape4)==T, Landscape4buf, Landscape4))
summary(is.na(Waders$Landscape4)) ## just check that all fields have been given a landscape label







##-----------------------------------------------##
#### 2. Calculate average standing water value ####
##-----------------------------------------------##


## Read in rasters of predicted standing water
SomW <- rast("RawData/PredictedWaterRasters/Som_StandingWaterAll_alt.tif")
NorfolkW <- rast("RawData/PredictedWaterRasters/Broads_StandingWaterAll.tif")
KentW <- rast("RawData/PredictedWaterRasters/NKent_StandingWaterAll.tif")
Essex1W <- rast("RawData/PredictedWaterRasters/Esx1_StandingWaterAll.tif")
Essex2W <- rast("RawData/PredictedWaterRasters/Esx2_StandingWaterAll.tif")

## Mosaic the three Greater Thames rasters
ThamesWat <- mosaic(KentW, Essex1W, Essex2W)
plot(ThamesWat)
rm(KentW, Essex1W, Essex2W); gc()


## Calculate the mean weighted pixel value for each field
## First mask the ditches around each field using field parcels buffered in by 10m
KentCanv <- st_read("CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp") |> st_buffer(dist = -10)
EssexCanv <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp") |> st_buffer(dist = -10)
ThamesCanv <- rbind(KentCanv, EssexCanv)
ThamesWat <- mask(ThamesWat, vect(ThamesCanv)); rm(ThamesCanv, KentCanv, EssexCanv)
WadersThames <- filter(Waders, Landscape == "Greater Thames")
WadersThames <- Prop_Raster_Cov(Survey = WadersThames, Rast = ThamesWat, Dist = BuffDists, Col = "Ave_WiderWater")


## Calculate the mean weighted pixel value for each field
## First mask the ditches around each field using field parcels buffered in by 10m
SomCanv <- st_read("CleanData/Scenarios/1-Starting Canvas/Som_Canvas.shp") |> st_buffer(dist = -10)
SomW <- mask(SomW, vect(SomCanv)); rm(SomCanv)
WadersSom <- filter(Waders, Landscape == "Somerset Levels and Moors")
WadersSom <- Prop_Raster_Cov(Survey = WadersSom, Rast = SomW, Dist = BuffDists, Col = "Ave_WiderWater")


## Calculate the mean weighted pixel value for each field
## First mask the ditches around each field using field parcels buffered in by 10m
BroadCanv <- st_read("CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp") |> st_buffer(dist = -10)
NorfolkW <- mask(NorfolkW, vect(BroadCanv)); rm(BroadCanv)
WadersBroads <- filter(Waders, Landscape == "Broads")
WadersBroads <- Prop_Raster_Cov(Survey = WadersBroads, Rast = NorfolkW, Dist = BuffDists, Col = "Ave_WiderWater")


## Bind all the different regions back together again
Waders <- rbind(WadersThames, WadersSom, WadersBroads)
## Change any NAN's to NA's in the data set
## I visually checked these fields and they all all generally NA's as the field were not lowland wet grassland
## Many were either arable (can see on Google maps) or saltmarsh type habitats
Waders$WaterCoverage <- ifelse(Waders$WaterCoverage == "NaN", NA, Waders$WaterCoverage)
rm(WadersThames, WadersSom, WadersBroads, WaterVals, ThamesWat, NorfolkW, SomW); gc()




##-----------------------------------------##
#### 3. Calculate percentage of woodland ####
##-----------------------------------------##

## Using UKCEH land cover data (Codes 1/2) calculate the percentage woodland around each field

## Reclassify the raster so woodland pixels are 1 and all others habitats are 0
m <- c(0, 2.5, 1,
       2.6, 22, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
WoodLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
plot(WoodLC) # check reclassification has worked


## Now calculate the proportion of woodland cover in each field buffer
Waders <- Prop_Raster_Cov(Survey = Waders, Rast = WoodLC, Dist = BuffDists, Col = "PropWood_")
rm(WoodLC); gc()




##--------------------------------------------##
#### 4. Calculate percentage of urban areas ####
##--------------------------------------------##

## Using UKCEH land cover data (Codes 20/21) calculate the percentage Urban cover around each field

## Reclassify the raster so urban pixels are 1 and all others habitats are 0
m <- c(0, 19.5, 0,
       19.6, 22, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
UrbanLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
plot(UrbanLC) # check reclassification has worked


## Now calculate the proportion of woodland cover in each field buffer
Waders <- Prop_Raster_Cov(Survey = Waders, Rast = UrbanLC, Dist = BuffDists, Col = "PropUrban_")
rm(UrbanLC); gc()




##---------------------------------------------------##
#### 5. Calculate percentage of intensive land use ####
##---------------------------------------------------##

# ## Using UKCEH land cover data (codes 3/4) and NE priority habitats (semi-natural grassland categories)
# ## calculate the percentage intensive land use around each field
# 
# ## Reclassify the raster so intensive land use pixels are 1 and all others habitats are 0
# m <- c(0, 2.5, 0,
#        2.5, 4.5, 1,
#        4.6, 22, 0) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# IntenseLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
# IntenseLC <- crop(IntenseLC, vect(st_buffer(st_as_sf(Waders), dist = 5000))) # crop to extent of survey fields
# plot(IntenseLC) # check reclassification has worked
# 
# 
# 
# ## Read in NE priority habitats
# PriHabs <- st_read(here("RawData", "NEPriorityHabs", "Priority_Habitat_Inventory_England.shp"))
# PriHabs <- st_crop(PriHabs, st_buffer(st_as_sf(Waders), dist = 5000)) # crop to study area bbox
# 
# ## Get tables of the main habitat types present
# table(PriHabs$mainhabs)
# AddHabs <- as.data.frame(table(PriHabs$addhabs))
# 
# ## Filter the priority habitats for semi natural grassland
# ## Create list of main habitats that are semi natural grassland
# MAINHABS <- c("Calaminarian grassland", "Coastal and floodplain grazing marsh", "Coastal and floodplain grazing marsh,Coastal saltmarsh",
#               "Coastal and floodplain grazing marsh,Lowland meadows",  "Good quality semi improved grassland", 
#               "Good quality semi improved grassland,Traditional orchard", "Grass moorland",
#               "Lowland calcareous grassland", "Lowland calcareous grassland,Maritime cliff and slope", 
#               "Lowland dry acid grassland", "Lowland dry acid grassland,Lowland heathland",
#               "Lowland meadows", "Purple moor grass and rush pastures", "Upland calcareous grassland")
# ## Create list of additional habitat codes that are semi natural grassland
# ADDHABS <- c("CALAM", "CFPGM", "GQSIG", "LCGRA", "LDAGR", "LMEAD", "PMGRP")
# 
# ## Create column of just the first 5 characters of the additional habitats column
# PriHabs <- PriHabs |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))
# ## Filter the data set
# SemiNatGrass <- filter(PriHabs, mainhabs %in% MAINHABS |
#                                (mainhabs == "No main habitat but additional habitats present" & addhabs  %in% ADDHABS))
# 
# 
# ## Rasterize the priority habitats layer
# ## If a pixel overlaps with a grassland priority habitat polygon then give the pixel get a value of -1 (no matter how much the pixel and polygon overlap)
# SemiNatGrass$Val <- -1
# SemiNatGrassr <- rasterize(vect(SemiNatGrass), IntenseLC, field = "Val", fun = "min",  background = 0)
# plot(SemiNatGrassr)
# 
# ## Add together the two raster 
# ## If a pixel was arable or improved grassland but also a priority grassland habitat then will get pixel value of 0
# ## If it was only arable or improved grassland then the pixel still has a value of 1
# Intense <- IntenseLC + SemiNatGrassr
# 
# ## recalssify so that all -1 pixels just become 0
# m <- c(-1.5, 0.5, 0,
#        0.6, 1.1, 1) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# Intense <- classify(Intense, rclmat, include.lowest=TRUE) 
# plot(Intense)
# 
# 
# ## Now calculate the proportion of intensive habitat uses in each field buffer
# Waders <- Prop_Raster_Cov(Survey = Waders, Rast = Intense, Dist = BuffDists, Col = "PropIntense_")
# rm(Intense, IntenseLC, SemiNatGrassr, PriHabs); gc()




##-------------------------------------------------------##
#### 6. Calculate percentage of semi-nat wet grassland ####
##-------------------------------------------------------##

# ## Reclassify the UKCEH raster so semi natural grassland pixels pixels are 1 and all others habitats are 0
# m <- c(0, 4.5, 0,
#        4.6, 7.5, 1,
#        7.6, 22, 0) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# SemiNatLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
# SemiNatLC <- crop(SemiNatLC, vect(st_buffer(st_as_sf(Waders), dist = 5000)))
# plot(SemiNatLC) # check reclassification has worked
# 
# 
# ## Rasterize the priority habitats layer
# ## If a pixel overlaps with a grassland priority habitat polygon then give the pixel get a value of 1
# SemiNatGrass$Val <- 1
# SemiNatGrassr <-rasterize(vect(SemiNatGrass), SemiNatLC, field = "Val", fun = "max",  background = 0)
# 
# 
# ## Add together two raster from the UKCEH and NE priority habitats and then reclassify so any value above 1 becomes 1
# SemiNatLC <- SemiNatLC + SemiNatGrassr
# m <- c(-0.5, 0.5, 0,
#        0.6, 2.5, 1) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# SemiNatLC <- classify(SemiNatLC, rclmat, include.lowest=TRUE) # reclassify the raster
# plot(SemiNatLC)
# 
# 
# 
# ## Read in the NATMAP soil vector
# Soil <- st_read(here("RawData", "Soil", "Soilscapes_England_Wales_27700.shp"))
# table(Soil$SOILSCAPE) # get table of soil types
# 
# ## Filter out seasonally wet soil
# ## Create vector of soil types to retains that are wet or peaty
# WetTypes <- c("Blanket bog peat soils", "Fen peat soils", 
#               "Lime-rich loamy and clayey soils with impeded drainage",
#               "Loamy and sandy soils with naturally high groundwater and a peaty surface",
#               "Loamy soils with naturally high groundwater",
#               "Loamy and clayey soils of coastal flats with naturally high groundwater",
#               "Loamy and clayey floodplain soils with naturally high groundwater",
#               "Naturally wet very acid sandy and loamy soils",
#               "Raised bog peat soils", "Very acid loamy upland soils with a wet peaty surface", 
#               "Slightly acid loamy and clayey soils with impeded drainage",
#               "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
#               "Slowly permeable seasonally wet acid loamy and clayey soils", 
#               "Slowly permeable wet very acid upland soils with a peaty surface")
# # filter data set
# WetSoils <- filter(Soil, SOILSCAPE %in% WetTypes)
# # plot(WetSoils$geometry)
# 
# ## Rasterize this wet/peaty soils layer, so if a wet/peaty are overlaps a pixel it gets a value of 1
# WetSoils$Val <- 1
# WetSoilsr <- rasterize(vect(WetSoils), SemiNatLC, field = "Val", fun = "max",  background = 0)
# plot(WetSoilsr)
# 
# 
# ## Add together raster to identify semi natural grassland + wet/peaty soils (pixel values of 2)
# SemiNatLC <- WetSoilsr + SemiNatLC
# ## re calssify raster so all values less then 2 become 0
# m <- c(-0.1, 1.5, 0, 1.6, 2.5, 1) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# SemiNatLC <- classify(SemiNatLC, rclmat, include.lowest=TRUE) # reclassify the raster
# plot(SemiNatLC)

## Read in and mosaic the lowland wet grassland rasters
LowWetGrass <- mosaic(rast("RawData/LowlandWetGrassRasters/Broads_LowWetGrass.tif"),
                      rast("RawData/LowlandWetGrassRasters/Essex_LowWetGrass.tif"),
                      rast("RawData/LowlandWetGrassRasters/NKent_LowWetGrass.tif"),
                      rast("RawData/LowlandWetGrassRasters/Som_LowWetGrass.tif"))

## Reclassify to 1s and 0s
m <- c(-0.1, 0.49, 0, 0.49, 1, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LowWetGrass <- classify(LowWetGrass, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LowWetGrass)

## Now calculate the proportion of semi natural wet grassland in each field buffer
Waders <- Prop_Raster_Cov(Survey = Waders, Rast = LowWetGrass, Dist = BuffDists, Col = "PropWetGrass_")
rm(LowWetGrass); gc()




##-------------------------------------------##
#### 7. Calculate percentage SSSI coverage ####
##-------------------------------------------##

# ## SSSI's 
# SSSI <- st_read(here("RawData", "PAs", "Sites_of_Special_Scientific_Interest_England.shp"))
# 
# ## Use function to calculate the proportion of each buffered field within a SSSI
# Waders <- Prop_Shapes_Cov(Survey = Waders, Shapes = SSSI, Dist = BuffDists, Col = "PropSSSI_")
# rm(SSSI); gc()




##----------------------------------------------##
#### 8. Calculate percentage reserve coverage ####
##----------------------------------------------##

## Read in the reserve outlines
## This is a combined multi-polygon of all RSPB reserves, LNRs and NNRs
Comb <- st_read(here("RawData", "Other Reserves", "All_Reserves_Poly.shp"))

## Use function to calculate the proportion of each buffered field within reserves
Waders <- Prop_Shapes_Cov(Survey = Waders, Shapes = Comb, Dist = BuffDists, Col = "PropReserve_")
rm(RSPB, LNR, NNR, IntComb, Comb); gc()




##------------------------------------------------##
#### 9. Calculate percentage wader AES coverage ####
##------------------------------------------------##

# ## Read in CS Schemes
# #CSS <- st_read(here("RawData", "Stewardship", "Countryside_Stewardship_Scheme_2016_Management_Areas_England.shp"))
# CSS_Ops <- st_read(here("RawData", "Stewardship", "Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp"))
# 
# ## Read in ES Schemes
# ESS <- st_read(here("RawData", "Stewardship", "Environmental_Stewardship_Scheme_Agreements_(England)___Natural_England.shp"))
# ESS <- filter(ESS, appstat == "Live")
# 
# 
# ## ESS COVERAGE ##
# 
# ## Add column to ESS data to indicate if agreements includes options that: 
# ## Directly benefit waders OR general options that could benefit waders
# ## See Table 3 of Rob's BWWM report for what is classified as "wader focused" and "general" options
# ESS <- ESS |> 
#   mutate(Wader_Focus = paste0(ifelse(is.na(str_match(optsonag, "HK9|HK10|HK11|HK12|HK13|HK14|HK19"))==T, "N", "Y")),
#          General_Focus = paste0(ifelse(is.na(str_match(optsonag, "DR|WDC|HK15|HK16|HK17|EK2|OK2|HK2|OHK2|EK3|OK3|HK3|OHK3|EK4|OK4|HK4|OHK4|HK6|HK7|HK8|HL8|HQ3|HQ4|HQ5|HQ6HQ7|HQ8|HR6|EK21|HK21|OHK21"))==T, "N", "Y")))
# 
# ## Calculate proportion of buffered field falling in ESS land With wader focused options
# ESSWa <- filter(ESS, Wader_Focus == "Y")
# Waders <- Prop_Shapes_Cov(Survey = Waders, Shapes = ESSWa, Dist = BuffDists, Col = "PropESSWa_")
# rm(ESS); gc()
# 
# 
# ## CSS COVERAGE ##
# 
# ## Add column to CSS data to indicate if agreements includes options that benefit waders OR general options that could benefit waders
# ## Filter out wader focused options
# CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))
# ## Filter out general options
# #CSS_OpsGe <- filter(CSS_Ops, opt_code %in% c("WN3", "WN4", "GS2", "GS5", "GS6", "GS7", "GS8", "UP2", "WT6", "WT7", "WT8", "WT9", "GS4", "GS13", "GS14", "GS16"))
# rm(CSS_Ops); gc()
# 
# 
# ## Read in the UK CEH land parcels data set
# Parcels <- st_read(here("RawData", "LandCover", "UKCEH_LCvector_2021.shp")) |> select(OBJECTID)
# Parcels <- st_crop(Parcels, st_buffer(st_as_sf(Waders), dist = 5000))
# 
# ## Determine which land parcels that have CSS wader-specific option points within them
# CSS_InterWa <- as.data.frame(st_covers(Parcels, CSS_OpsWa))
# Parcels <- Parcels[CSS_InterWa$row.id,]
# 
# ## Calculate proportion of buffered fields falling in CSS land With wader focused options
# Waders <- Prop_Shapes_Cov(Survey = Waders, Shapes = Parcels, Dist = BuffDists, Col = "PropCSSWa_")
# rm(CSS_InterWa); gc()
# 
# 
# ## ALL STEWARDSHIP COVERAGE ##
# 
# ## COULD DO WITH CHECKING THE JOINING OF POLYGONS HERE
# ## Combine ESS wader focused and CSS wader focused land parcels
# ESSWa <- st_crop(ESSWa, st_buffer(st_as_sf(Waders), dist = 5000)) # crop to speed up processing
# ESSWa <- vect(ESSWa)
# Parcels <- vect(Parcels)
# StewComb <- c(ESSWa, Parcels)
# StewComb <- terra::aggregate(vect(StewComb)) |> st_as_sf()
# 
# ## Calculate proportion of buffered fields falling in CSS land With wader focused options
# Waders <- Prop_Shapes_Cov(Survey = Waders, Shapes = StewComb, Dist = BuffDists, Col = "PropStewWa_")
# rm(StewComb, CSS_OpsWa, ESSWa); gc()




##------------------------------------------------##
#### 10. Save data set with landscape attributes ####
##------------------------------------------------##

## Check i didn't lose any data here
SL == nrow(Waders)

## Save Data
Waders |> st_drop_geometry() |>  select(-geometry) |> 
write_csv(here("CleanData", "Wader Abundance", "4-AddLandscapeAttributes/", "Breeding_Pairs_FullAttrib3.csv"))
write_sf(st_as_sf(Waders), here("CleanData", "Wader Abundance", "4-AddLandscapeAttributes/", "Breeding_Pairs_FullAttrib3.shp"))

