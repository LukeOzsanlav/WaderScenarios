##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 12/12/23
## Goal: Add attributes to survey fields
## 
## Breakdown:
## 1. Calculate Crow/Magpie Density
## 2. Calculate distance to inter tidal habitat
## 3. Merge with predator fence spreadsheet
## 4. Extract water coverage from predicted rasters
## 5. Identify fields in reserves
## 6. Identify fields in a SSSI
## 7. Identify fields in stewardship
##
##
##------------------------------------------------------##

## Queries:
## Are there yearly data sets for AES? Rob's report makes it out that there are?
## The CSS data sets don't quite link up, some options do not overlap land parcels...

## Load in required packages
pacman::p_load(here, tidyverse, data.table, sf, terra, gt)

## Use `here` package for specifying file paths
here::i_am("Code/Wader Abundance/3 - Add Field Attributes.R")
library(here)
here()



##---------------------##
#### 0. Data Read In ####
##---------------------##

## Read in filtered breeding pairs estimates
Waders <- read_csv(here("CleanData", "Wader Abundance", "2-FieldSelectionPlotting", "Breeding_Pairs_Selected.csv"))

## read in shaepfile for BWWM fields and join to breeding waders pairs
BWWM <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp")) |> select(F_LOC_ID)
Waders <- left_join(Waders, BWWM, by = "F_LOC_ID")
SL <- nrow(Waders) # keep track of data set length throughout

## Read in the UKCEH landcover data
LC <- rast(here("RawData", "LandCover", "gblcm25m2021.tif"))
LC <- LC[[1]]

## Read in the UK 1km grid (goes with predator abundance data)
Grid1KM <- st_read(here("RawData", "Gridlines", "OSGB_Grid_1km.shp"))

## Read in the bird densities modeled from the BTO BBS survey
## These are provided as point at the center of each 1km square
BTODens <- read_csv(here("RawData", "BTO Density Data", "model_popn_dens_0709.csv"))

## Read in my spreadsheet of predator fences
Fences <- read_csv(here("RawData", "Predator Fences", "RSPB Predator Fence Locaitons.csv")) |> filter(is.na(F_LOC_ID)==F)


## Read in the various reserve outlines
## RSPB reserves
RSPB <- st_read(here("RawData", "RSPB Reserves", "EnglandWales_RSPBReserves.shp"))

## Local Nature reserves
LNR <- st_read(here("RawData", "Other Reserves", "Local_Nature_Reserves_EnglandPolygon.shp"))

## National Nature Reserves
NNR <- st_read(here("RawData", "Other Reserves", "National_Nature_Reserves_EnglandPolygon.shp"))

## SSSI's 
SSSI <- st_read(here("RawData", "PAs", "Sites_of_Special_Scientific_Interest_England.shp"))

## Read in the Different Stewardship Schemes
CSS <- st_read(here("RawData", "Stewardship", "Countryside_Stewardship_Scheme_2016_Management_Areas_England.shp"))
CSS_Ops <- st_read(here("RawData", "Stewardship", "Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp"))
ESS <- st_read(here("RawData", "Stewardship", "Environmental_Stewardship_Scheme_Agreements_(England)___Natural_England.shp"))



##----------------------------##
#### 1. Crow/Magpie Density ####
##----------------------------##

## streamline the BTO breeding bird density data
BTODens <- BTODens |>  select(GridRef, MGdens, `C.dens`)

## Join the 1km square shapefile with the BTO data
## Each 1km square has a unique code
Grid1KM <- left_join(Grid1KM, BTODens, by = join_by(PLAN_NO == GridRef))

## Convert the vector 1km grid into a raster 
template = rast(vect(Grid1KM), res = 1000) # create temple

## Extract the density of Magpies and Crows for each 1km pixel
Magpie <- rasterize(vect(Grid1KM), template, field = "MGdens") # magpie density
plot(Magpie, col=grDevices::hcl.colors(50, palette = "Sunset"))

Crow <- rasterize(vect(Grid1KM), template, field = "C.dens") # Crow density
plot(Crow, col=grDevices::hcl.colors(50, palette = "Sunset"))
rm(BTODens, Grid1KM, template); gc()



## Calculate the mean weighted Magpie density for each field
## Crop raster to study areas and fill in NA pixels in a buffer
Magpie <- crop(Magpie, vect(st_buffer(st_as_sf(Waders), dist = 1000))) # crop to save space
Magpie <- focal(x = Magpie, w = 3, fun = "mean", na.policy = "only", na.rm = T)

## run Magpie extraction
MagDens <- extract(Magpie, vect(st_as_sf(Waders)), fun = mean, na.rm=TRUE) 
summary(MagDens$focal_mean); hist(MagDens$focal_mean) # checks

## assign density back to main data set
Waders$MagDens <- MagDens$focal_mean
rm(Magpie, MagDens);gc() # clean environment



## Calculate the mean weighted Carrion Crow density for each field
## Crop raster to study areas and fill in NA pixels in a buffer
Crow <- crop(Crow, vect(st_buffer(st_as_sf(Waders), dist = 1000))) # crop to save space
Crow <- focal(x = Crow, w = 3, fun = "mean", na.policy = "only", na.rm = T)

## run Crow extraction
CrowDens <- extract(Crow, vect(st_as_sf(Waders)), fun = mean)
summary(CrowDens$focal_mean); hist(CrowDens$focal_mean) # checks

## assign density back to main data set
Waders$CrowDens <- CrowDens$focal_mean
rm(Crow, CrowDens);gc() # clean environment




##----------------------------------------##
#### 2. Distance to Inter tidal Habitat ####
##----------------------------------------##

## crop the UK CEH data to south of England
LC <- crop(LC, vect(st_buffer(st_as_sf(Waders), dist = 10000)))
plot(LC)

## Reclassify the raster so that inter tidal habitats are 1 and all others habitats are NA
m <- c(0, 12.5, NA,
       12.6, 13.5, 1,
       13.6, 14.5, NA, 
       14.6, 19.5, 1, 
       19.6, 22, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
rc1 <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
plot(rc1) # check reclassification has worked

## Calculate distance from NA pixels (None inter tidal) to none NA pixel (inter tidal habitat)
DistToInter <- distance(rc1)
plot(DistToInter)

## Calculate the mean weighted pixel value for each field
DistsToInter <- extract(DistToInter, vect(st_as_sf(Waders)), fun = mean, weights = T)
hist(DistsToInter$gblcm25m2021_1)

## Assign these weighted means to the main data set
Waders$InterTidal_Distm <- DistsToInter$gblcm25m2021_1

## make a plot to check this worked
WadersBroad <- st_as_sf(Waders) |> filter(Landscape == "Greater Thames")
ggplot() + geom_sf(data = WadersBroad, mapping = aes(geometry = geometry, fill = InterTidal_Distm), colour = NA) + theme_minimal()
rm(DistsToInter, rc1, rclmat, m); gc()




##--------------------------------##
#### 3. Predator Fence Presence ####
##--------------------------------##

## Manipulate predator fence to join with breeding pair estimate
colnames(Fences)
Fences <- Fences |> select(F_LOC_ID, FENCE_COVERAGE, Year) |> rename(Fence_Coverage = FENCE_COVERAGE, year = Year)

## Join to the breeding pair estimates
Waders <- left_join(Waders, Fences, by = c("F_LOC_ID", "year"))
Waders <- Waders |> mutate(Fence_Coverage = ifelse(is.na(Fence_Coverage)==T, 0, Fence_Coverage)) # change the NAs to 0's
gc()



##----------------------------------##
#### 4. Standing water total edge ####
##----------------------------------##

## Read in shape files of water polygons for each region
## This is made by creating a binary water raster and then turning the raster into polygons
NorfolkWS <- st_read("RawData/BinaryWaterRasters/Broads_BinaryWaterAll.shp")
SomWS <- st_read("RawData/BinaryWaterRasters/Som_BinaryWaterAll.shp")
NKentWS <- st_read("RawData/BinaryWaterRasters/NKent_BinaryWaterAll.shp")
EssexWS <- st_read("RawData/BinaryWaterRasters/Essex_BinaryWaterAll.shp")
ThamesWS <- rbind(NKentWS, EssexWS)


## Calculate the perimeter of all water polygons that fall within a field
WadersBroads <- filter(Waders, Landscape == "Broads") |> st_as_sf()
BroadsOver <- st_buffer(WadersBroads, dist = -10) |> 
              st_intersection(NorfolkWS) |> 
              mutate(WaterEdge = perim(vect(geometry))) |> 
              select(F_LOC_ID, WaterEdge) |> 
              st_drop_geometry()

## Summaries the total water edge within each field
BroadsOver <- BroadsOver |> group_by(F_LOC_ID) |> summarise(WaterEdge = sum(WaterEdge))

## Now join total water edge length onto the fields data set, change NAs to 0m edge
WadersBroads <- left_join(WadersBroads, BroadsOver, by = "F_LOC_ID") |> 
                mutate(WaterEdge = ifelse(is.na(WaterEdge)==T, 0, WaterEdge)) 



## Calculate the perimeter of all water polygons that fall within a field
WadersSom <- filter(Waders, Landscape == "Somerset Levels and Moors") |> st_as_sf()
SomOver <- st_buffer(WadersSom, dist = -10) |> 
              st_intersection(SomWS) |> 
              mutate(WaterEdge = perim(vect(geometry))) |> 
              select(F_LOC_ID, WaterEdge) |> 
              st_drop_geometry()

## Summaries the total water edge within each field
SomOver <- SomOver |> group_by(F_LOC_ID) |> summarise(WaterEdge = sum(WaterEdge))

## Now join total water edge length onto the fields data set, change NAs to 0m edge
WadersSom <- left_join(WadersSom, SomOver, by = "F_LOC_ID") |> 
                mutate(WaterEdge = ifelse(is.na(WaterEdge)==T, 0, WaterEdge)) 



## Calculate the perimeter of all water polygons that fall within a field
WadersThames <- filter(Waders, Landscape == "Greater Thames") |> st_as_sf()
ThamesOver <- st_buffer(WadersThames, dist = -10) |> 
              st_intersection(ThamesWS) |> 
              mutate(WaterEdge = perim(vect(geometry))) |> 
              select(F_LOC_ID, WaterEdge) |> 
              st_drop_geometry()

## Summaries the total water edge within each field
ThamesOver <- ThamesOver |> group_by(F_LOC_ID) |> summarise(WaterEdge = sum(WaterEdge))

## Now join total water edge length onto the fields data set, change NAs to 0m edge
WadersThames <- left_join(WadersThames, ThamesOver, by = "F_LOC_ID") |> 
                mutate(WaterEdge = ifelse(is.na(WaterEdge)==T, 0, WaterEdge)) 


## Bind all the different regions back together again
Waders <- rbind(WadersThames, WadersSom, WadersBroads)
rm(NorfolkWS, SomWS, NKentWS, EssexWS, ThamesWS, ThamesOver, SomOver, BroadsOver); gc()




##-------------------------------------##
#### 5. Average standing water value ####
##-------------------------------------##

## Read in rasters of predicted standing water
SomW <- rast("RawData/PredictedWaterRasters/Som_StandingWaterAll_alt.tif")
NorfolkW <- rast("RawData/PredictedWaterRasters/Broads_StandingWaterAll.tif")
KentW <- rast("RawData/PredictedWaterRasters/NKent_StandingWaterAll.tif")
Essex1W <- rast("RawData/PredictedWaterRasters/Esx1_StandingWaterAll.tif")
Essex2W <- rast("RawData/PredictedWaterRasters/Esx2_StandingWaterAll.tif")

## Mosaic the three greater thames rasters
ThamesWat <- mosaic(KentW, Essex1W, Essex2W)
plot(ThamesWat)
rm(KentW, Essex1W, Essex2W); gc()


## Calculate the mean weighted pixel value for each field
WadersThames <- filter(Waders, Landscape == "Greater Thames")
WaterVals <- extract(ThamesWat, vect(st_buffer(st_as_sf(WadersThames), dist = -10)), fun = mean, weights=T, na.rm=TRUE)
hist(WaterVals$mNDWI)

## Assign these weighted means to the main data set
WadersThames$WaterCoverage <- WaterVals$mNDWI


## Calculate the mean weighted pixel value for each field
WadersSom <- filter(Waders, Landscape == "Somerset Levels and Moors")
WaterVals <- extract(SomW, vect(st_buffer(st_as_sf(WadersSom), dist = -10)), fun = mean, weights=T, na.rm=TRUE)
hist(WaterVals$mNDWI)

## Assign these weighted means to the main data set
WadersSom$WaterCoverage <- WaterVals$mNDWI


## Calculate the mean weighted pixel value for each field
WadersBroads <- filter(Waders, Landscape == "Broads")
WaterVals <- extract(NorfolkW, vect(st_buffer(st_as_sf(WadersBroads), dist = -10)), fun = mean, weights=T, na.rm=TRUE)
hist(WaterVals$mNDWI)

## Assign these weighted means to the main data set
WadersBroads$WaterCoverage <- WaterVals$mNDWI


## Bind all the different regions back together again
Waders <- rbind(WadersThames, WadersSom, WadersBroads)
## Change any NAN's to NA's in the data set
## I visually checked these fields and they all all generally NA's as the field were not lowland wet grassland
## Many were either arable (can see on Google maps) or saltmarsh type habitats
Waders$WaterCoverage <- ifelse(Waders$WaterCoverage == "NaN", NA, Waders$WaterCoverage)
rm(WadersThames, WadersSom, WadersBroads, WaterVals, ThamesWat, NorfolkW, SomW); gc()




##------------------------------------------##
#### 6. Identify Status of Survey Fields  ####
##------------------------------------------##

##---------------------------------##
#### 6.1 Function: Overlap ID'er ####
##---------------------------------##

## Fields = polygons of the fields with breeding pair estimates from the three priority Landscapes
## OverlapShapes = polygons of the thing thing that you want to determine if the fields overlap with
## Thresh = the threshold for proportion overlap, above which fields are considered to be in OverlapShapes
## ColName = The name of column used to indicate whether the field is within OverlapShapes
Field_Overlap_IDer <- function(Fields = Fields, OverlapShapes = OverlapShapes, 
                               Thresh = Thresh, ColName = ColName){
  
  ## Crop OverlapShapes to region covered by Fields
  OverlapShapes <- st_crop(OverlapShapes, st_as_sf(Fields))
  plot(OverlapShapes$geometry)
  
  ## Ensure that the fields with breeding pair estimates is an sf object
  Fields <- st_as_sf(Fields)
  class(Fields)
  
  ## Work out which fields are in Reserves and what the size of the overlap is
  Overlap_sub <- Fields |> filter(Somerset %in% c("All", "S")) |> 
                 st_intersection(OverlapShapes) |> 
                 mutate(OverlapArea = st_area(geometry)) |> 
                 select(F_LOC_ID, OverlapArea) |> 
                 st_drop_geometry()
  
  ## Remove any duplicated in this data set, can't work out why duplicates happen, maybe overlapping shapes?
  Dups <- Overlap_sub |> select(F_LOC_ID, OverlapArea) |> duplicated()
  Overlap_sub <- Overlap_sub[Dups==FALSE,]
  Overlap_sub <- Overlap_sub |> group_by(F_LOC_ID) |> summarise(OverlapArea = sum(OverlapArea))

  ## Now join the overlap size onto the main fields data &
  ## Work out the proportion of each field in a reserve
  Fields <- left_join(Fields, Overlap_sub, by = "F_LOC_ID") |> 
             mutate(Area = st_area(geometry),
                    PropOverlap = ifelse(is.na(OverlapArea)==T, 0, OverlapArea/Area), 
                    Inside = ifelse(PropOverlap >= Thresh, "Y", "N")) |> 
             select(-c(Area, OverlapArea, PropOverlap)) 
  
  ## Plot just to check this 
  ggplot() + geom_sf(data = Fields, mapping = aes(geometry = geometry, fill = Inside)) + theme_minimal()
  
  ## change the column Name to one of your choosing
  colnames(Fields)[colnames(Fields)== "Inside"] <- ColName
  
  
  ## return this object
  return(Fields)
  
}



##-------------------------------------------------------##
#### 6.2 Identify Fields that Overlap with Peaty Soils ####
##-------------------------------------------------------##

## Read in Soil data here as it is large and coudl crash R
Soils <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## filter out just the soils that are in some way peaty
Soils <- filter(Soils, SOILSCAPE %in% c("Blanket bog peat soils", "Fen peat soils", "Loamy and sandy soils with naturally high groundwater and a peaty surface",
  "Raised bog peat soils", "Slowly permeable wet very acid upland soils with a peaty surface",
  "Very acid loamy upland soils with a wet peaty surface"))

## Now label fields as peaty of they are more then 50% covered by the peaty soils polygon
Soils <- st_as_sf(terra::aggregate(vect(Soils)))
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = Soils, 
                              Thresh = 0.5, ColName = "Peat_Soil")

## Plot to check that this worked in Somerset
Som <- Waders |> filter(Landscape == "Somerset Levels and Moors")
OverlapSoils <- st_crop(Soils, st_as_sf(Som))
ggplot() + geom_sf(data =Som, aes(fill = Peat_Soil), colour = NA) + geom_sf(data = OverlapSoils, fill = NA, colour = "orange", linewidth = 1) + theme_light()
rm(Soils, OverlapSoils, Som); rm()



##-----------------------------------------------------------##
#### 6.3 Identify Fields that Overlap with Protected Sites ####
##-----------------------------------------------------------##

## RSPB reserve overlap
RSPB <- st_as_sf(terra::aggregate(vect(RSPB)))
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = RSPB, 
                              Thresh = 0.5, ColName = "RSPB_Reserve")
## LNR's
LNR <- st_as_sf(terra::aggregate(vect(LNR)))
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = LNR, 
                              Thresh = 0.5, ColName = "LNR")
## NNR's
NNR <- st_as_sf(terra::aggregate(vect(NNR)))
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = NNR, 
                              Thresh = 0.5, ColName = "NNR")
## SSSI's 
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = SSSI, 
                              Thresh = 0.5, ColName = "SSSI")

## plots to check
DataCheck <- filter(Waders, Landscape == "Greater Thames") # Somerset Levels and Moors # Broads
ggplot() + geom_sf(data = DataCheck, mapping = aes(geometry = geometry, fill = RSPB_Reserve), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck, mapping = aes(geometry = geometry, fill = LNR), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck, mapping = aes(geometry = geometry, fill = NNR), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck, mapping = aes(geometry = geometry, fill = SSSI), colour = NA) + theme_minimal()




##---------------------------------------------------------------##
#### 6.4 Identify Fields that Overlap with Stewardship Schemes ####
##---------------------------------------------------------------##

## Check Types of ESS and CSS
table(CSS$cs_type)
table(ESS$scheme)
table(CSS_Ops$opt_code)


## Identify CSS Higher Tier
CSSH <- filter(CSS, cs_type == "Higher Tier")
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = CSSH, 
                              Thresh = 0.5, ColName = "CSS_Higher")
rm(CSSH);gc()


## Identify CSS Mid Tier
CSSM <- filter(CSS, cs_type == "Mid Tier")
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = CSSM, 
                              Thresh = 0.5, ColName = "CSS_Mid")
rm(CSSM);gc()


## Identify ESS Higher Tier or Entry Level plus Higher Level Stewardship
ESSH <- filter(ESS, scheme %in% c("Entry Level plus Higher Level Stewardship", "Higher Level Stewardship"))
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = ESSH, 
                              Thresh = 0.5, ColName = "ESS_Higher")
rm(ESSH);gc()


## plots to check
DataCheck2 <- filter(Waders, Landscape == "Greater Thames") # Somerset Levels and Moors # Broads
ggplot() + geom_sf(data = DataCheck2, mapping = aes(geometry = geometry, fill = CSS_Higher), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck2, mapping = aes(geometry = geometry, fill = CSS_Mid), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck2, mapping = aes(geometry = geometry, fill = ESS_Higher), colour = NA) + theme_minimal()



## Add column to ESS data to indicate if agreements includes options that benefit waders or general options that could benefit waders
## See Table 3 of Rob's BWWM report for what is classified as "wader focused" and "general" options
ESS <- ESS |> 
  mutate(Wader_Focus = paste0(ifelse(is.na(str_match(optsonag, "HK9|HK10|HK11|HK12|HK13|HK14|HK19"))==T, "N", "Y")),
         General_Focus = paste0(ifelse(is.na(str_match(optsonag, "DR|WDC|HK15|HK16|HK17|EK2|OK2|HK2|OHK2|EK3|OK3|HK3|OHK3|EK4|OK4|HK4|OHK4|HK6|HK7|HK8|HL8|HQ3|HQ4|HQ5|HQ6HQ7|HQ8|HR6|EK21|HK21|OHK21"))==T, "N", "Y")))

## Identify ESS With wader focused options
ESSWa <- filter(ESS, Wader_Focus == "Y")
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = ESSWa, 
                              Thresh = 0.5, ColName = "ESS_Wader")
rm(ESSWa);gc()

## Identify ESS With general options benefiting waders
ESSGen <- filter(ESS, General_Focus == "Y")
Waders <- Field_Overlap_IDer(Fields = Waders, OverlapShapes = ESSGen, 
                              Thresh = 0.5, ColName = "ESS_Gen")
rm(ESSGen);gc()


## plots to check
DataCheck3 <- filter(Waders, Landscape == "Greater Thames") # Somerset Levels and Moors # Broads
ggplot() + geom_sf(data = DataCheck3, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()
ggplot() + geom_sf(data = DataCheck3, mapping = aes(geometry = geometry, fill = ESS_Gen), colour = NA) + theme_minimal()



## Add column to CSS data to indicate if agreements includes options that benefit waders or general options that could benefit waders
## See Table 3 of Rob's BWWM report for what is classified as "wader focused" and "general" options
## Filter out wader focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))
## Filter out general options
CSS_OpsGe <- filter(CSS_Ops, opt_code %in% c("WN3", "WN4", "GS2", "GS5", "GS6", "GS7", "GS8", "UP2", "WT6", "WT7", "WT8", "WT9", "GS4", "GS13", "GS14", "GS16"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
CSS_InterWa <- as.data.frame(st_covers(Waders, CSS_OpsWa))
Waders$CSS_Wader <- "N"
Waders$CSS_Wader[CSS_InterWa$row.id] <- "Y"

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a general CSS point fell within it
CSS_InterGe <- as.data.frame(st_covers(Waders, CSS_OpsGe))
Waders$CSS_Gen <- "N"
Waders$CSS_Gen[CSS_InterGe$row.id] <- "Y"


## plots to check
DataCheck4 <- filter(Waders, Landscape == "Somerset Levels and Moors") # Greater Thames # Broads
CSS_OpsWa2 <- st_crop(CSS_OpsWa, DataCheck4) # crop to size
CSS_OpsGe2 <- st_crop(CSS_OpsGe, DataCheck4) # crop to size
ggplot() + geom_sf(data = DataCheck4, mapping = aes(geometry = geometry, fill = CSS_Wader), colour = NA) + 
           geom_sf(data = CSS_OpsWa2, mapping = aes(geometry = geometry), fill = "purple") + 
           theme_minimal()
ggplot() + geom_sf(data = DataCheck4, mapping = aes(geometry = geometry, fill = CSS_Gen), colour = NA) +
           geom_sf(data = CSS_OpsGe2, mapping = aes(geometry = geometry, fill = "purple")) + 
           theme_minimal()





##-------------------------##
#### 7. Save Survey Data ####
##-------------------------##

## Check that I did not lose any rows in the data set
SL == nrow(Waders)

## Write out the main survey data set
Waders |>  select(-geometry) |> 
write_csv("CleanData/Wader Abundance/AddFieldAttributes-3/Breeding_Pairs_Attribs.csv")
write_sf(st_as_sf(Waders), "CleanData/Wader Abundance/AddFieldAttributes-3/Breeding_Pairs_Attribs.shp")







