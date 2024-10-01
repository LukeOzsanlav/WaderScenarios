##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 22/04/2024
## 
## 
## Aim: Determine the opportunity areas for conservation actions
##      Particularly area for AES/Reserve creation in grassland and arable reversion
## 
##------------------------------------------------------##


## Load in package required
pacman::p_load(sf, terra, tidyverse, exactextractr)
source("Code/Helper functions.R")



##-----------------------------------##
#### 0.1 Read in general data sets ####
##-----------------------------------##

## Read in the priority landscape bounding boxes
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")

## Read in the UKCEH landcover data set @ 25m resolution for 2021/21/23
LC21 <- rast("RawData/LandCover/gblcm25m2021.tif")
LC21 <- LC21[[1]]
LC22 <- rast("RawData/LandCover/gblcm2022_25m.tif")
LC22 <- LC22[[1]]
LC23 <- rast("RawData/LandCover/gblcm2023_25m.tif")
LC23 <- LC23[[1]]

## Read in the Natural England priority habitat inventory
PriHabs <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")

## Read in the NATMAP soil vector
Soil <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## Read in the BWWM field shapefile
BWWM <- st_read("RawData/BWWM field shapefile/BWWM_fields_27Sep2023.shp")



##-----------------------------------------------##
#### 0.2 Function to define conservation areas ####
##-----------------------------------------------##
## This function identifies area of grassland and arable that are potential sites for wader management
## This is based of elevation, soil type and existing habitat type/quality
OppArea <- function(ROI,             # polygons around priority landscape 
                    DEM,             # Digital elevation model for landscape
                    LC21=LC21,       # UK landcover raster for 2021, use same data for each landscape
                    LC22=LC22,       # UK landcover raster for 2022, use same data for each landscape
                    LC23=LC23,       # UK landcover raster for 2023, use same data for each landscape
                    PriHabs=PriHabs, # NE priority habitat shape file for the UK
                    BWWM=BWWM,       # All UK BWWM fields as a shape file
                    Soil=Soil,       # UK NATMAP soil vectors
                    MainHab,         # NE priority habitats categories that could be lowland wet grass
                    AddHabs,         # NE priority habitats categories that could be lowland wet grass (short codes)
                    WetSoils,        # NATMAP soil categories that could be wet soils
                    GrassMask,       # Digitised grassland that needs to be masked from action areas
                    GrassAdd=NULL,   # Digitised extra grassland
                    outpath          # Outpath for arable and grassland opportunity, include landscape name at end of path
                    ){ ## Start of OppArea function
  
  
  ##---------------------##
  #### Landscape Setup ####
  ##---------------------##
  
  ## Get the outline of the priority landscape
  plot(ROI$geometry)
  
  ## Crop the UKCEH land cover map 
  ## This can be used to rasterise polygons 
  LC_ROI <- terra::crop(LC21, vect(st_buffer(ROI$geometry, dist = 1000)))
  LC22ROI <- terra::crop(LC22, vect(st_buffer(ROI$geometry, dist = 1000)))
  LC23ROI <- terra::crop(LC23, vect(st_buffer(ROI$geometry, dist = 1000)))
  plot(LC_ROI, main = "Landcover Map 2021")
  
  ## Rasterize the ROI
  ROIr <- rasterize(vect(st_buffer(ROI$geometry, dist = 50)), LC_ROI, cover = T, background = 0)
  ROIr <- Half_reclass(ROIr)
  plot(ROIr, main = "ROI rasterized")
  
  
  ##---------------------------------##
  #### ID land at suitable heights ####
  ##---------------------------------##
  
  ## Extract all the BWWM fields for that region in the landscape
  BWWMROI <- st_intersection(BWWM, st_buffer(st_as_sf(ROI), dist = 10)) # crop to study area bbox
  plot(BWWMROI$geometry, main = "BWWM inside priority landscape")
  
  ## Mask the elevation raster using the BWWm field shapes
  ElevM <- crop(DEM, vect(BWWMROI$geometry))
  ElevM <- mask(ElevM, vect(BWWMROI$geometry))
  plot(ElevM, main = "BWWM fields elevation")
  gc()
  
  ## Extract the 99.5th quantile for all elevation values in BWWM fields
  quantvals <- quantile(values(ElevM), probs = c(0.995), na.rm = T)
  minheight <- -100
  maxheight <- 2000
  gc()
  
  ## Reclassify the raster so all pixel below the threshold have a value of 1  
  m <- c(minheight-1, as.numeric(quantvals[1]), 1, 
         as.numeric(quantvals[1]), maxheight+1, 0) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  ElevBin <- classify(DEM, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(ElevBin, main = "Suitable Height Land")
  
  rm(DEM, ElevM, BWWMROI, quantvals, minheight, maxheight); gc() # free up space
  
  
  
  ##--------------------------------------##
  #### ID Suitable NE priority habitats ####
  ##--------------------------------------##
  
  ## Crop the NE priority habitats
  PriHabsROI <- st_crop(PriHabs, st_buffer(st_as_sf(ROI), dist = 1000)) # crop to study area bbox
  
  ## Get tables of the main habitat types present
  table(PriHabsROI$mainhabs)
  AddHabs <- as.data.frame(table(PriHabsROI$addhabs))
  
  ## Filter the priority habitats for semi natural grassland
  ## Create list of main habitats that are semi natural grassland
  ## NOTES: Rob included Lowland calcareous grassland and Lowland dry acid grassland
  ##        These are present in the ROI seem far away from the priority landscape and again don't look right on satellite view 
  MAINHABS <- MainHab
  ## Create list of additional habitat codes that are semi natural grassland
  ADDHABS <- AddHabs
  
  ## Create column of just the first 5 characters of the additional habitats column
  PriHabsROI <- PriHabsROI |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))
  ## Filter the data set
  SNPrty <- filter(PriHabsROI, mainhabs %in% MAINHABS |
                                 (mainhabs == "No main habitat but additional habitats present" & addhabs  %in% ADDHABS))
  
  ## Rasterize the priority habitats layer
  ## If a pixel more than 50% overlaps with a grassland priority habitat polygon then give the pixel get a value of 1 
  SNPrtyr <- rasterize(vect(SNPrty), LC_ROI, cover = T, background = 0)
  SNGrass <- Half_reclass(SNPrtyr)
  plot(SNGrass, main = "Suitable NE Priority Habitats")
  
  rm(ADDHABS, MAINHABS, SNPrtyr, SNPrty, AddHabs, PriHabsROI); gc() # free up space
  
  
  
  ##-------------------------------##
  #### ID wet soil using NATMAP  ####
  ##-------------------------------##
  
  ## Crop the NATMAP soil vector
  SoilROI <- st_crop(Soil, st_buffer(st_as_sf(ROI), dist = 1000)) # crop to study area bbox
  table(SoilROI$SOILSCAPE) # get table of soil types
  
  ## Filter out seasonally wet soil
  ## Create vector of soil types to retains that are wet or peaty
  WetTypes <- WetSoils
  
  ## filter data set
  WetSoils <- filter(SoilROI, SOILSCAPE %in% WetTypes)
  ggplot(data =WetSoils) + geom_sf(fill= "blue") + theme_bw() + ggtitle("Wet Soils")
  
  ## Rasterize this wet/peaty soils layer, so if a wet/peaty are overlaps a pixel it gets a value of 1
  WetSoils <- rasterize(vect(WetSoils), LC_ROI, cover = T, background = 0)
  WetSoils <- Half_reclass(WetSoils)
  plot(WetSoils, main = "Wet Soils Raster")
  
  rm(SoilROI, WetTypes); gc() # free up space
  
  
  
  ##-----------------------------------------------##
  #### Combine Semi-nat grass & wet soils layers ####
  ##-----------------------------------------------##
  
  ## Add together raster to identify semi natural grassland + wet/peaty soils (pixel values of 2)
  SNWetGr <- WetSoils + SNGrass
  plot(SNWetGr)
  
  ## re classify raster so all values less then 2 become 0
  m <- c(0, 1.5, 0, 1.5, 2.5, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  SNWetGr <- classify(SNWetGr, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(SNWetGr, main = "Semi nat wet grassland")
  
  
  
  ##--------------------##
  #### ID BWWM fields ####
  ##--------------------##
  
  ## Crop BWWM fields
  BWWMROI <- st_crop(BWWM, st_buffer(st_as_sf(ROI), dist = 1000)) # crop to study area bbox
  
  ## Rasterize the BWWM fields
  ## If a pixel more than 50% overlaps with a BWWM field polygon then give the pixel get a value of 1 
  BWWMROI <- rasterize(vect(BWWMROI), LC_ROI, cover = T, background = 0)
  BWWMROI <- Half_reclass(BWWMROI)
  plot(BWWMROI, main = "BWWM Raster")
  
  
  
  ##----------------------------------------------##
  #### Combine Semi-nat wet grass & BWWM layers ####
  ##----------------------------------------------##
  
  ## Add together the BWWM fields and the Semi-natural wet grassland (ID'd using NE priority habitat, NATMAP and UKCEH)
  SNWetGr <- SNWetGr + BWWMROI
  plot(SNWetGr, main = "NE wet grass plus BWWM fields")
  
  ## Reclassify so that any pixel that was wet grassland from either layer becomes a 1
  m <- c(0, 0.5, 0, 0.5, 2.5, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  SNWetGr <- classify(SNWetGr, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(SNWetGr, main = "All Semi nat wet grassland")
  
  rm(BWWMROI); gc() # free up space
  
  
  
  ##---------------------------------##
  #### Mask using UKCEH Land cover ####
  ##---------------------------------##
  
  ## Reclassify the UKCEH raster so woodland, arable, coastal, water pixels are all given a value of 1
  ## Only Semi-natural grassland pixels have a value of 0
  m <- c(0, 3.5, 1,
         3.5, 8.5, 0,
         8.5, 22, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  
  ## For all three years find none-grassland pixels and then add together all the layers
  OtherLC <- classify(LC_ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(OtherLC) # check reclassification has worked
  
  LC22Other <- classify(LC22ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(LC22Other) # check reclassification has worked
  
  LC23Other <- classify(LC23ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(LC23Other) # check reclassification has worked
  
  ## Add all three years together
  OtherLC <- OtherLC + LC22Other + LC23Other
  plot(OtherLC, main = "All other habitats combined")
  rm(LC22Other, LC23Other); gc()
  
  
  ## Take away the `other habitat` raster from my semi natural wet grassland raster
  SNWetGr2 <- SNWetGr-OtherLC
  plot(SNWetGr2,  main = "Semi nat wet grassland minus none grassland habitats")
  
  ## Reclassify so that only pixels that are definitely semi-nat wet grass remain as a 1
  m <- c(-5, 0.5, 0, 0.5, 1.5, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  SNWetGr2 <- classify(SNWetGr2, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(SNWetGr2,  main = "Semi nat wet grassland minus none grassland habitats cleaned")
  
  ## Reasmple the elevation raster (binary as indicted if right elevation or not), to same resolution as wet grass layer
  ElevBinGr <- terra::resample(ElevBin, SNWetGr2)
  m <- c(0, 0.5, 0, 0.5, 1, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  ElevBinGr <- classify(ElevBinGr, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(ElevBinGr, main = "Suitable Height Binary")
  
  ## Now add together seasonally wet lowland grass and suitable elevation tasters to find suitable habitat
  SNWetGr2 <- ElevBinGr + SNWetGr2
  m <- c(0, 1.5, 0, 1.5, 3, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  SNWetGr2 <- classify(SNWetGr2, rclmat, include.lowest=TRUE) # reclassify the raster
  values(SNWetGr2) <- ifelse(is.na(values(SNWetGr2))==T, 0, values(SNWetGr2))
  plot(SNWetGr2, main = "Semi nat wet grassland at suitable height")
  
  ## Finally mask using a manually digitized layer
  ROI_Mask <- GrassMask
  SNWetGr2 <- mask(SNWetGr2, vect(ROI_Mask), inverse = T, updatevalue = 0)
  plot(SNWetGr2, main = "Semi nat wet grassland at suitable height minus digitized land")
  
  ## Add any bit of grassland that were manually digitised
  if(is.null(GrassAdd)==FALSE){
    
    ## And fill in a missing gap using a manually digitised layer
    GrassAdd2 <- rasterize(vect(GrassAdd), LC_ROI, cover = T, background = 0)
    GrassAdd2 <- Half_reclass(GrassAdd2)
    plot(GrassAdd2)
    SNWetGr2 <- max(SNWetGr2, GrassAdd2, na.rm = T)
    plot(SNWetGr2, main = "Semi nat wet grassland at suitable height add extra land")
    
  }
  
  ## Finally multiply by the ROI raster to remove any fields outside of the priority landscape
  SNWetGr2 <- SNWetGr2*ROIr
  plot(SNWetGr2, main = "Semi nat wet grass, cropped to ROI")
  
  ## Write out the raster to use later
  writeRaster(SNWetGr2, paste0(outpath, "LowWetGrass.tif"), overwrite=TRUE)
  rm(OtherLC, ROI_Mask, ElevBinGr, SNWetGr, GrassMask); gc() # free up space
  
  
  
  ##---------------------------------##
  #### ID All suitable Arable land ####
  ##---------------------------------##
  
  ## Reclassify the UKCEH raster so arable land is all given a value of 1
  ## All other habitats  have a value of 0
  m <- c(0, 2.5, 0,
         2.5, 3.5, 1,
         3.5, 22, 0) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  
  ## For all three years find arable pixels and then add together all the layers
  OtherLC <- classify(LC_ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(OtherLC) # check reclassification has worked
  
  LC22Other <- classify(LC22ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(LC22Other) # check reclassification has worked
  
  LC23Other <- classify(LC23ROI, rclmat, include.lowest=T) # reclassify the raster
  plot(LC23Other) # check reclassification has worked
  
  ## Combine all three years
  OtherLC <- OtherLC + LC22Other + LC23Other
  plot(OtherLC)
  rm(LC22Other, LC23Other); gc()
  
  ## Reclassify so that if a pixel was arable in any year then it becomes a 1
  m <- c(0, 0.5, 0, 0.5, 4, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  ArableLC <- classify(OtherLC, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(ArableLC, main = "Arable fields")
  
  
  ## Now add together seasonally wet soils and Arable land to find arable land on suitable soil types
  ArableSoil <- WetSoils + ArableLC
  m <- c(0, 1.5, 0, 1.5, 3, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  ArableSoil <- classify(ArableSoil, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(ArableSoil, main = "Arable fields on wet soils")
  
  ## Finally re-sample my binary elevation dataset to the same res as the ArableSoil data
  ElevBin2 <- terra::resample(ElevBin, ArableSoil)
  m <- c(0, 0.5, 0, 0.5, 1, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  ElevBin2 <- classify(ElevBin2, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(ElevBin2)
  
  ## Add together Elevation and Arable soil data sets to find Arable land with suitable soil types and at suitable elevations
  SuitArab <- ElevBin2 + ArableSoil
  m <- c(0, 1.5, 0, 1.5, 3, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  SuitArab <- classify(SuitArab, rclmat, include.lowest=TRUE) # reclassify the raster
  values(SuitArab) <- ifelse(is.na(values(SuitArab))==T, 0, values(SuitArab))
  plot(SuitArab, main = "Arable Opportunity Land")
  
  ## Multiply by the ROI raster to remove any fields outside of the priority landscape
  SuitArab <- SuitArab*ROIr
  plot(SuitArab, main = "Arable Opportunity, cropped to ROI")
  
  ## Some bits of grassland where manually added above
  ## Some of these could be arable fields that have been wrongly classifed by UKCEH
  ## SO just apply this extra grassland as a mask here
  if(is.null(GrassAdd)==FALSE){
    
    ## Mask Arable layer using the additonal manually digitized grassland
    SuitArab <- mask(SuitArab, vect(GrassAdd), inverse = T, updatevalue = 0)
    plot(SuitArab, main = "Arable Opportunity, extra grass mask")
    
  }
  
  ## final plot to check
  plot(SNWetGr2+SuitArab, main = "Grass + Arable opporunity (should only be 0's and 1's)")
  
  ## Write out the raster to use later
  writeRaster(SuitArab, paste0(outpath, "ArableSuitable.tif"), overwrite=TRUE)
  rm(ElevBin2, ArableSoil, OtherLC, LC_ROI, LC22ROI, LC23ROI); gc() # free up space
  
  
} ## End of OppArea function





##----------------------------------------##
#### 1.0 Opportunity Areas: Essex Coast ####
##----------------------------------------##

## Get the outline of the priority landscape
Essex <- filter(MyBoxes, Att3Value == "Essex")
plot(Essex$geometry)

## Read in the elevation data for the landscape (plus 1000m buffer)
ThamesDEM <- rast("RawData/LiDAR/Thames_DTM_2m.tif")
EssexDEM <- crop(ThamesDEM, vect(Essex))
rm(ThamesDEM)


## Create list of main habitats that are semi natural grassland
## NOTES: Rob included Lowland calcareous grassland and Lowland dry acid grassland
##        These are present in the Essex seem far away from the priority landscape and again don't look right on satellite view 
EssexMainHab <- c("Coastal and floodplain grazing marsh", "Coastal and floodplain grazing marsh,Lowland meadows",   
                    "Good quality semi improved grassland", 
                    "Lowland meadows",
                    "Purple moor grass and rush pastures")

## Create list of additional habitat codes that are semi natural grassland
EssexAddHabs <- c("CFPGM", "GQSIG", "LMEAD", "PMGRP")

## Filter out seasonally wet soil
## Create vector of soil types to retains that are wet or peaty
EssexWetSoils <- c("Fen peat soils",
                    "Lime-rich loamy and clayey soils with impeded drainage",
                    "Loamy and clayey soils of coastal flats with naturally high groundwater",
                    "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                    "Naturally wet very acid sandy and loamy soils",
                    "Slightly acid loamy and clayey soils with impeded drainage", 
                    "Slowly permeable seasonally wet acid loamy and clayey soils",
                    "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
                    "water")

## Read in digitized layer of areas that are definitely not lowland wet grassland
EssexGrassMask <- st_read("RawData/ManualLandCoverMasks/Essex_NoneGrass_Mask.shp")

## Read in digitized layer of areas that are lowland wet grassland that were missed
EssexGrassAdd <- st_read("RawData/ManualLandCoverMasks/Essex_ExtraGrass.shp")
gc()



# ROI=Essex
# DEM=EssexDEM
# MainHab=EssexMainHab
# AddHabs=EssexAddHabs
# WetSoils=EssexWetSoils
# LC21=LC21
# LC22=LC22
# LC23=LC23
# PriHabs=PriHabs
# BWWM=BWWM
# Soil=Soil
# GrassMask=EssexGrassMask
# GrassAdd=EssexGrassAdd
# outpath="CleanData/Scenarios/3-DefineActionAreas/Essex_"

## Run the function for Essex
OppArea(ROI=Essex, DEM=EssexDEM, 
        MainHab=EssexMainHab, AddHabs=EssexAddHabs, WetSoils=EssexWetSoils, 
        LC21=LC21, LC22=LC22, LC23=LC23, PriHabs=PriHabs, BWWM=BWWM, Soil=Soil,
        GrassMask=EssexGrassMask, GrassAdd=EssexGrassAdd,
        outpath="CleanData/Scenarios/3-DefineActionAreas/Essex_")

## Plot the outputs of the function
EssexGrass <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_LowWetGrass.tif")
plot(EssexGrass, main = "Essex grassland opportunity")
EssexArab <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_ArableSuitable.tif")
plot(EssexArab, main = "Essex arable opportunity")
rm(EssexGrass, EssexArab, EssexDEM, Essex); gc()





##-------------------------------------##
#### 2.0 Opportunity Areas: Somerset ####
##-------------------------------------##

## Get the outline of the priority landscape
Som <- filter(MyBoxes, Att3Value == "Somerset Levels and Moors")
plot(Som$geometry)

## Read in the elevation data for the landscape (plus 1000m buffer)
SomDEM <- rast("RawData/LiDAR/Somerset_DTM_2m.tif")

## Create list of main habitats that are semi natural grassland
## NOTES: Rob included Lowland calcareous grassland and Lowland dry acid grassland
##        These are present in the Som seem far away from the priority landscape and again don't look right on satellite view 
SomMainHab <- c("Coastal and floodplain grazing marsh", "Coastal and floodplain grazing marsh,Lowland meadows",   
              "Good quality semi improved grassland", 
              "Lowland meadows",
              "Purple moor grass and rush pastures")

## Create list of additional habitat codes that are semi natural grassland
SomAddHabs <- c("CFPGM", "GQSIG", "LMEAD", "PMGRP")

## Filter out seasonally wet soil
## Create vector of soil types to retains that are wet or peaty
SomWetSoils <- c("Fen peat soils",
                  "Lime-rich loamy and clayey soils with impeded drainage",
                  "Loamy and clayey floodplain soils with naturally high groundwater",
                  "Loamy and clayey soils of coastal flats with naturally high groundwater ",
                  "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                  "Raised bog peat soils",
                  "Slightly acid loamy and clayey soils with impeded drainage",
                  "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
                  "water")

## Read in digitized layer of areas that are definately not lowland wet grassland
SomGrassMask <- st_read("RawData/ManualLandCoverMasks/Som_NoneGrass_Mask.shp")



## Run the function for Somerset Levels
OppArea(ROI=Som, DEM=SomDEM, 
        LC21=LC21, LC22=LC22, LC23=LC23,
        PriHabs=PriHabs, BWWM=BWWM, Soil=Soil,
        MainHab=SomMainHab, AddHabs=SomAddHabs, WetSoils=SomWetSoils, 
        GrassMask=SomGrassMask, GrassAdd=NULL,
        outpath="CleanData/Scenarios/3-DefineActionAreas/Som_")

## PLot the outputs of the function
SomGrass <- rast("CleanData/Scenarios/3-DefineActionAreas/Som_LowWetGrass.tif")
plot(SomGrass)
SomArab <- rast("CleanData/Scenarios/3-DefineActionAreas/Som_ArableSuitable.tif")
plot(SomArab)
rm(SomGrass, SomArab, SomDEM, Som); gc()




##-----------------------------------##
#### 3.0 Opportunity Areas: Broads ####
##-----------------------------------##

## Get the outline of the priority landscape
Broads <- filter(MyBoxes, Att3Value == "Broads")
plot(Broads$geometry)

## Read in the elevation data for the landscape (plus 1000m buffer)
BroadsDEM <- rast("RawData/LiDAR/Broads_DTM_2m.tif")

## Create list of main habitats that are semi natural grassland
## NOTES: Rob included Lowland calcareous grassland and Lowland dry acid grassland
##        These are present in the Broads seem far away from the priority landscape and again don't look right on satellite view 
BroadsMainHab <- c("Coastal and floodplain grazing marsh", "Coastal and floodplain grazing marsh,Lowland meadows",   
                    "Good quality semi improved grassland", 
                    "Lowland meadows",
                    "Purple moor grass and rush pastures")

## Create list of additional habitat codes that are semi natural grassland
BroadsAddHabs <- c("CFPGM", "GQSIG", "LMEAD", "PMGRP")

## Filter out seasonally wet soil
## Create vector of soil types to retains that are wet or peaty
BroadsWetSoils <- c("Fen peat soils",
                    "Lime-rich loamy and clayey soils with impeded drainage",
                    "Loamy and clayey soils of coastal flats with naturally high groundwater",
                    "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                    "Naturally wet very acid sandy and loamy soils",
                    "Slightly acid loamy and clayey soils with impeded drainage", 
                    "Slowly permeable seasonally wet acid loamy and clayey soils",
                    "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
                    "water")

## Read in digitized layer of areas that are definitely not lowland wet grassland
BroadsGrassMask <- st_read("RawData/ManualLandCoverMasks/Broads_NoneGrass_Mask.shp")

## Read in digitized layer of areas that are lowland wet grassland that were missed
BroadsGrassAdd <- st_read("RawData/ManualLandCoverMasks/Broads_ExtraGrass.shp")



## Run the function for Broads Levels
OppArea(ROI=Broads, DEM=BroadsDEM, 
        LC21=LC21, LC22=LC22, LC23=LC23, PriHabs=PriHabs, BWWM=BWWM, Soil=Soil,
        MainHab=BroadsMainHab, AddHabs=BroadsAddHabs, WetSoils=BroadsWetSoils, 
        GrassMask=BroadsGrassMask, GrassAdd=BroadsGrassAdd,
        outpath="CleanData/Scenarios/3-DefineActionAreas/Broads_")

## Plot the outputs of the function
BroadsGrass <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_LowWetGrass.tif")
plot(BroadsGrass, main = "Broads grassland opportunity")
BroadsArab <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_ArableSuitable.tif")
plot(BroadsArab, main = "Broads arable opportunity")
rm(BroadsGrass, BroadsArab, BroadsDEM, Broads); gc()




##---------------------------------------##
#### 4.0 Opportunity Areas: North Kent ####
##---------------------------------------##

## Get the outline of the priority landscape
NKent <- filter(MyBoxes, Att3Value == "North Kent")
plot(NKent$geometry)

## Read in the elevation data for the landscape (plus 1000m buffer)
ThamesDEM <- rast("RawData/LiDAR/Thames_DTM_2m.tif")
NKentDEM <- crop(ThamesDEM, vect(NKent))
rm(ThamesDEM)


## Create list of main habitats that are semi natural grassland
## NOTES: Rob included Lowland calcareous grassland and Lowland dry acid grassland
##        These are present in the NKent seem far away from the priority landscape and again don't look right on satellite view 
NKentMainHab <- c("Coastal and floodplain grazing marsh",   
                  "Good quality semi improved grassland", 
                  "Lowland meadows")

## Create list of additional habitat codes that are semi natural grassland
NKentAddHabs <- c("CFPGM", "GQSIG", "LMEAD")

## Filter out seasonally wet soil
## Create vector of soil types to retains that are wet or peaty
NKentWetSoils <- c("Loamy and clayey floodplain soils with naturally high groundwater",
                    "Loamy and clayey soils of coastal flats with naturally high groundwater",
                    "Loamy soils with naturally high groundwater",
                    "Slightly acid loamy and clayey soils with impeded drainage",
                    "Slowly permeable seasonally wet acid loamy and clayey soils",
                    "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
                    "water")

## Read in digitized layer of areas that are definitely not lowland wet grassland
NKentGrassMask <- st_read("RawData/ManualLandCoverMasks/NKent_NoneGrass_Mask.shp")

## Read in digitized layer of areas that are lowland wet grassland that were missed
NKentGrassAdd <- st_read("RawData/ManualLandCoverMasks/NKent_ExtraGrass.shp")



## Run the function for North Kent
OppArea(ROI=NKent, DEM=NKentDEM, 
        LC21=LC21, LC22=LC22, LC23=LC23, PriHabs=PriHabs, BWWM=BWWM, Soil=Soil,
        MainHab=NKentMainHab, AddHabs=NKentAddHabs, WetSoils=NKentWetSoils, 
        GrassMask=NKentGrassMask, GrassAdd=NKentGrassAdd,
        outpath="CleanData/Scenarios/3-DefineActionAreas/NKent_")

## Plot the outputs of the function
NKentGrass <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_LowWetGrass.tif")
plot(NKentGrass, main = "North Kent grassland opportunity")
NKentArab <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_ArableSuitable.tif")
plot(NKentArab, main = "North Kent arable opportunity")
rm(NKentGrass, NKentArab, NKentDEM, NKent); gc()



