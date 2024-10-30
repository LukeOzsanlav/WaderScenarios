##----------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 09/05/2024
##
## Aim: Create all guidelines for scenario modelling of the Broads Landscape
## 
##----------------------------------------------------------## 


## Load in packages
pacman::p_load(tidyverse, sf, terra, tidyterra, gstat, stars, automap, spatstat, ggspatial, leastcostpath)
options(scipen=999) # turn off scientific notation

## Load in various helper functions
source("Code/Helper functions.R")




##----------------------------------##
#### 0. Read in general data sets ####
##----------------------------------##

## Read in my priority landscapes boundary and buffer it by 3000m
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
BrOutline <- MyBoxes |> filter(Att3Value == "Broads")
Br <- MyBoxes |> filter(Att3Value == "Broads") |> st_buffer(dist=3000)
rm(MyBoxes) # free up space

## Read in UKCEH landcover and then crop and mask it with the landscape boundary
## This will be used as the canvas to map all other variables across
LCfull <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]]
LC_Br <- crop(LCfull, vect(Br)) |> mask(vect(Br))
rm(LCfull)
plot(LC_Br) # free up space




##-------------------------------##
#### 1. Better only guidelines ####
##-------------------------------##

##*-- Target smallest existing populations --*##

## Read in full BWWM data set that has each wader cluster labelled
Reg_Fields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Broads")

## Read in filtered breeding pairs estimates
## And calculate total wader abundance for each field
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
Waders <- Waders |> mutate(Tot_abund = rowSums(dplyr::select(Waders, est_pairsL, est_pairsR), na.rm = T)) |> 
                    select(F_LOC_ID, Tot_abund) |> group_by(F_LOC_ID) |> 
                    summarise(Tot_abund = sum(Tot_abund, na.rm = T))

## Join pair estimates onto the wader cluster fields
## Summarise the total number of pairs per cluster
PopSize <- left_join(Reg_Fields, Waders, by = "F_LOC_ID") |> st_drop_geometry() |> 
               group_by(ClustGroup) |> 
               summarise(ClustPop = sum(Tot_abund, na.rm = T)) |> 
               left_join(Reg_Fields, ., by = "ClustGroup")
ggplot(data= PopSize) + geom_sf(mapping = aes(geometry=geometry, fill = ClustPop), colour = NA) + theme_light()
rm(Waders, Reg_Fields); gc() # free up space

## Write out this variable for later use
write_csv(PopSize, "CleanData/Guideline Creation/Norfolk/Broads_ClusterPopSize.csv")




##-------------------------------##
#### 2. Bigger only guidelines ####
##-------------------------------##

##*-- Target areas near sites with more breeding waders for expansion --*##

## Loop through each cluster
for(j in 1:length(unique(PopSize$ClustGroup))){

   ## Creat raster of buffered cluster where in pixel inside the buffered cluster recieves the value equal to the number of breeding wader pairs  
   Subbuf <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 2000)), LC_Br, field = "ClustPop", fun = "max",  background = NA)
   
   ## Create the same raster as above but this time do it for the unbuffered cluster (this is so we can calcualte the distance to the cluster)
   Subno <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 1)), LC_Br, field = "ClustPop", fun = "max",  background = NA)
   
   ## Calculate the distance of all pixels to the the un buffered cluster
   ## then mask it so that I only retain pixels that are within 2km of the cluster
   ClusterD <- distance(Subno, unit="m") |> mask(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 2000)))
   ## Now take the inverse scale of the distacne within the 2km buffer, this gives pixels within or right next to the cluster a score of 1 and those 2km away a score of 0
   values(ClusterD) <- Inv_scale_vals(values(ClusterD))

   ## Multiple the raster of the inverse distance to the cluster and the population size of the cluster
   Final <- Subbuf*ClusterD

   ## Now for each loop make sure to take the maximum value for any pixels across all clusters
   if(j==1){AllRast <- Final}else{AllRast <- max(AllRast, Final, na.rm = T)}

}

## With final raster make sure if is scaled and then NA pixels are given a value of 0
values(AllRast) <- scale_vals(values(AllRast))
values(AllRast) <- ifelse(is.na((values(AllRast)))==T, 0, values(AllRast))
plot(AllRast)
ExpandClust <- AllRast

## free up space
rm(AllRast, Final, Subbuf, Subno, ClusterD); gc()




##-----------------------------##
#### 3. More only guidelines ####
##-----------------------------##




##------------------------------------##
#### 4. Bigger/More only guidelines ####
##------------------------------------##

##*-- Target area with silty soils --*##

## Read in the NATMAP soil vector
Soil <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## Create vector of soil types to retains that are silty and wet
## The rule from the stakeholder said silty soils but if the silty soils are freely draining then they will be useless for wader management
table(Soil$SOILSCAPE) # get table of soil types
SiltWetTypes <- c("Lime-rich loamy and clayey soils with impeded drainage",
                  "Loamy soils with naturally high groundwater",
                  "Loamy and clayey soils of coastal flats with naturally high groundwater",
                  "Loamy and clayey floodplain soils with naturally high groundwater",
                  "Naturally wet very acid sandy and loamy soils",
                  "Slightly acid loamy and clayey soils with impeded drainage",
                  "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
                  "Slowly permeable seasonally wet acid loamy and clayey soils")

## filter data set, and assign a Val column for rasterisation
SiltSoils <- filter(Soil, SOILSCAPE %in% SiltWetTypes) #|> mutate(Val = 1)

## Rasterize the silty soils polygons
SiltSoil <- rasterize(vect(SiltSoils), LC_Br, cover = T, background = 0)
## If a pixel is more than 50% covered by silty soils then assign a pixel value of 1, using my function
SiltSoil <- Half_reclass(SiltSoil)
plot(SiltSoil)
rm(Soil, SiltWetTypes, SiltSoils); gc()



##*-- Target areas surrounded by less woodland --*##

## Reclassify land cover raster so that it indicates if pixel is woodland or not
LC_BrWid <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]] |>  crop(vect(Br))
m <- c(0.0, 2.5, 1,
       2.5, 23, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Wood <- classify(LC_BrWid, rclmat, include.lowest=TRUE)
plot(Wood)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
## then calculate inverse scaled values as we want areas with lower tree cover
WoodDens <- focal(Wood, w=41, fun="mean", na.rm = T)
values(WoodDens) <- Inv_scale_vals(values(WoodDens))
plot(WoodDens) # plot smoothed raster
rm(LC_BrWid, Wood); gc() # free up space



##*-- Target areas surrounded by less urban areas --*##

## Reclassify land cover raster so that it indicates if pixel is woodland or not
LC_BrWid <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]] |>  crop(vect(Br))
m <- c(0.0, 19.5, 0,
       19.5, 23, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Urban <- classify(LC_BrWid, rclmat, include.lowest=TRUE)
plot(Urban)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
## then calculate inverse scaled values as we want areas with lower tree cover
UrbanDens <- focal(Urban, w=41, fun="mean", na.rm = T)
values(UrbanDens) <- Inv_scale_vals(values(UrbanDens))
plot(UrbanDens) # plot smoothed raster
rm(LC_BrWid, Urban); gc() # free up space



##*-- Target areas with wintering waterbird CS agreements already --*##

## Read in the Different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")
ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESS <- rbind(ESSWest, ESSEast)
rm(ESSWest, ESSEast)

## Read in the Canvas as this has the field parcel shapes 
Canv <- st_read("CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp") |> select(geometry)

## Filter out winter water birds focused options and breeding waders water options
ESS_OpsWint <- filter(ESS, optcode %in% c("HK10", "HK12", "HK14"))
ESS_OpsBreed <- filter(ESS, optcode %in% c("HK9", "HK11", "HK13"))

## Remove any winter options that also have a 
ESS_OpsWint <- filter(ESS_OpsWint, !parcref %in% ESS_OpsBreed$parcref)

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
ESS_InterWa <- as.data.frame(st_covers(Canv, ESS_OpsBreed))
Canv$ESS_Wader <- "N"
Canv$ESS_Wader[ESS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()


## Filter out winter waterbird focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS10", "GS12"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
CSS_InterWa <- as.data.frame(st_covers(Canv, CSS_OpsWa))
Canv$CSS_Wader <- "N"
Canv$CSS_Wader[CSS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = CSS_Wader), colour = NA) +
           theme_minimal()

## Rasterize this wet/peaty soils layer, so if a wet/peaty are overlaps a pixel it gets a value of 1
Canv <- filter(Canv, CSS_Wader == "Y" | ESS_Wader == "Y")
AESWint <- rasterize(vect(Canv), LC_Br, cover = T, background = 0)
AESWint <- Half_reclass(AESWint)
plot(AESWint, main = "Winter scheme only AES fields")

## free up space
rm(Canv, ESS, ESS_InterWa, CSS_InterWa, CSS_Ops, CSS_OpsWa, ESS_OpsWint, ESS_OpsBreed); gc()



##*-- Target areas that link existing wader sites --*##

## Read in polygons of the wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Broads/Broads_Wader_Sites.shp")
WaderSites <- WaderSites |> st_centroid()
plot(WaderSites$geometry)

## Read in raster of potential lowland wet grassland and suitable arable land
## Then combine so have tiered values
LowGrass <- rast("RawData/LowlandWetGrassRasters/Broads_LowWetGrass.tif")
ArableSuit <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_ArableSuitable.tif") |> resample(LowGrass)
LowGrass <- max(LowGrass*4, ArableSuit*2) + 1
plot(LowGrass)

## Creates a cost surface using the values in the supplied SpatRaster. 
Cond <- create_cs(LowGrass, neighbours = 16)
plot(Cond)

## Calculates Least-cost paths from-everywhere-to-everywhere. This is based on the approach proposed by White and Barber (2012).
LCP <- create_FETE_lcps(x=Cond, locations=WaderSites)
plot(LCP)
LCP$Val <- 1


## No turn the lines that connect polygons into a raster
## Create a base raster with a pixel size of 2km using my base land cover raster
LC_Br2km <- aggregate(LC_Br, fact = 80)

## Rasterize the lines connecting polygons using my 2km raster as a base
Connectr <- rasterize(vect(LCP), LC_Br2km, field = "Val", fun = "sum",  background = 0)
ggplot() + geom_spatraster(data = Connectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)

## Scale the values from the rasterized lines
values(Connectr) <- scale_vals(values(Connectr))
plot(Connectr)

## Finally resample the raster back down to 25m
Connectr <- resample(Connectr, LC_Br)

## Finally turn any NAs back to zeros
values(Connectr) <- ifelse(is.na(values(Connectr))==T, 0, values(Connectr))
plot(Connectr)

## free up space
rm(WaderSites, LowGrass, Cond, LCP, LC_Br2km); gc()






##-------------------------------------##
#### 5. Arable conversion guidelines ####
##-------------------------------------##

##*-- Target arable reversion in isolated patches within grassland --*##

## Read in raster of potential lowland wet grassland 
LowGrass <- rast("RawData/LowlandWetGrassRasters/Broads_LowWetGrass.tif")
plot(LowGrass)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
LowGrassDens <- focal(LowGrass, w=41, fun="mean", na.rm = T)
LowGrassDens <- mask(LowGrassDens, LowGrass) # crop to the extent of the original raster
plot(LowGrassDens) # plot smoothed raster
rm(LowGrass); gc()



##*-- Target arable conversion on areas adjacent to wet grassland --*##

## Rasterize the broads wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Broads/Broads_Wader_Sites.shp")
WaderSites <- rasterize(vect(WaderSites), LC_Br, cover = T, background = 0)
plot(WaderSites)

## Any pixels that overlap with a wader site by more then 50% becomes a 1 and less than 50% a 0
m <- c(0, 0.5, NA,
       0.5, 1, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
WaderSites <- classify(WaderSites, rclmat, include.lowest=TRUE)
plot(WaderSites)

## Calculate distance from NA pixels (other) to none NA pixel (wader sites)
WaderSites <- distance(WaderSites, unit="m")

## Create a raster where only pixels with value of <1000m have a value of 1
m <- c(0, 1000, 1,
       1000, 200000, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
WaderSites2 <- classify(WaderSites, rclmat, include.lowest=TRUE)
plot(WaderSites2)

## Mask pixels only with 1km of wader sites
WaderSitesDist <- terra::mask(WaderSites, WaderSites2) # crop to ROI
values(WaderSitesDist) <- Inv_scale_vals(values(WaderSitesDist)) # inverse scale so low distances have high value
values(WaderSitesDist) <- ifelse(is.na(values(WaderSitesDist))==T, 0, values(WaderSitesDist))
plot(WaderSitesDist)
rm(WaderSites2, WaderSites); gc()



##*-- Target low grade arable land --*##

## Grades of different agricultural land
AgriGrade <- st_read("RawData/AgriGrades/Agricultural_Land_Classification_Provisional_EnglandPolygon.shp")

## Add new column for scoring of grading
AgriGrade <- filter(AgriGrade, !alc_grade == "Urban")
AgriGrade$Grade <- ifelse(AgriGrade$alc_grade %in% c("Grade 1"), 0,
                         ifelse(AgriGrade$alc_grade %in% c("Grade 2"), 1,
                            ifelse(AgriGrade$alc_grade %in% c("Grade 3"), 2,
                                  ifelse(AgriGrade$alc_grade %in% c("Grade 4", "Non Agricultural"), 3, 0))))


## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1
AgriGrades <- rasterize(vect(AgriGrade), LC_Br, field = "Grade", fun = "min", background = 0)
values(AgriGrades) <- scale_vals(values(AgriGrades))
plot(AgriGrades)
rm(AgriGrade); gc()



##*--Target areas that were originally grassland --*##

# ## 2007 data set
# LC2007 <- rast("RawData/LandCover/lcm2007gb25m.tif")[[1]]
# 
# ## crop the UK CEH data to polygon of interest
# LC2007 <- crop(LC2007, vect(st_buffer(st_as_sf(Br), dist = 500)))
# plot(LC2007)


## 2000 data set, need to re download this data set as it will not read in
LC2000 <- rast("RawData/LandCover/LCM2000_GB_25m.tif")[[1]]

## crop the UK CEH data to polygon of interest
LC2000 <- crop(LC2000, vect(st_buffer(st_as_sf(Br), dist = 500)))
plot(LC2000)


## 1990 data set
LC1990<- rast("RawData/LandCover/gb1990lcm25m.tif")[[1]]

## crop the UK CEH data to polygon of interest
LC1990 <- crop(LC1990, vect(st_buffer(st_as_sf(Br), dist = 500)))
plot(LC1990)


# ## 2007 raster, reclassify so grassland is 1 and all other habitats are 0
# m <- c(0, 3.5, 0,
#        3.5, 9.5, 1,
#        9.5, 23, 0) # matrix for re-classification
# rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
# LC2007_Grass <- classify(LC2007, rclmat, include.lowest=TRUE) # reclassify the raster
# plot(LC2007_Grass)


## 2000 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 49, 0,
       49, 92, 1,
       92, 110, 0,
       110, 112, 1, 
       112, 250, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC2000_Grass <- classify(LC2000, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC2000_Grass)


## 1990 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 3.5, 0,
       3.5, 8.5, 1,
       8.5, 23, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC1990_Grass <- classify(LC1990, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC1990_Grass)


## Add together the two grassland layers
All_Grass <- (LC1990_Grass+LC2000_Grass)/2
HistoricGrass <- resample(All_Grass, LC_Br) # resample so that the extent is the same as all my other grading layers
plot(HistoricGrass)
rm(All_Grass, LC1990_Grass, LC2000_Grass, LC1990, LC2000); gc()





##--------------------------------##
#### 6. All strategy guidelines ####
##--------------------------------##

##*-- Target sites in larger areas of continuous wet grassland --*##

## Read in raster of potential lowland wet grassland 
LowGrass <- rast("RawData/LowlandWetGrassRasters/Broads_LowWetGrass.tif")
plot(LowGrass)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
LowGrassDens <- focal(LowGrass, w=41, fun="mean", na.rm = T)
LowGrassDens <- mask(LowGrassDens, LowGrass) # crop to the extent of the original raster
plot(LowGrassDens) # plot smoothed raster
rm(LowGrass); gc()



##*-- Target lowest land within hydro units --*##

## 2m LiDAR raster for the Somerset
Br_Elev <- rast("RawData/LiDAR/Broads_DTM_2m.tif")

## Shapefile of Hydro Units 
Br_Units <- st_read("RawData/HydroUnits/Broads/Broads_HydroUnits.shp")
 
## Add values to Somerset hydro units data set
Br_Units$Region <- "Broads"
Br_Units$id <- 1:nrow(Br_Units)
plot(Br_Units$geometry)

## create an empty list to put all the raster outputs into
List <- vector(mode='list', length= nrow(Br_Units))


## Loop through each hydro unit and extract the LiDAR data and grade it
for(j in 1:nrow(Br_Units)){
  
  message("iteration ", j)
  
  ## filter out the hydro unit outline
  Br_Unit <- vect(Br_Units[j,])
  
  ## Crop the hydro unit from the LiDAR data
  Unit_Lidar <- crop(Br_Elev, Br_Unit, mask = T)
  
  ## Calculate empirical cumulative distribution function
  x <- values(Unit_Lidar)
  y <- ecdf(x)(x)
  
  ## Take inverse of distribution function so lowest land has value of 1 and highest land has value of 1
  values(Unit_Lidar) <- Inv_scale_vals(y)
  plot(Unit_Lidar)
  
  ## add raster to slot in list
  List[[j]] <- Unit_Lidar
  
}


## free up some space first
rm(Br_Elev, Unit_Lidar); gc()

## Combine all the raster into a raster collection
List <- sprc(List)

## now merge that raster collection
Hydro_HeightLow <- terra::merge(List)

## Change all NAs in the raster to zeros
values(Hydro_HeightLow) <- ifelse(is.na(values(Hydro_HeightLow))==T, 0, values(Hydro_HeightLow))
gc()

## re-sample the fine-scale topography raster to my base raster
Hydro_HeightLow <- resample(Hydro_HeightLow, LC_Br)
values(Hydro_HeightLow) <- ifelse(is.na(values(Hydro_HeightLow))==T, 0, values(Hydro_HeightLow))
plot(Hydro_HeightLow, col=grDevices::hcl.colors(50, palette = "Sunset"))
rm(Br_Units, Br_Unit, List); gc() # free up some space



##*-- Target areas within small hydro units --*##

## Shape file of Hydro Units 
## Units with 2 pumps instead of one have Unit == "2Pump"
## I have multiplied the sizes of these units by two to also incorporate the fact
## that another rule proposed favoring units with fewer pumps
Br_Units <- st_read("RawData/HydroUnits/Broads/Broads_HydroUnits.shp") |> 
            mutate(area = as.numeric(st_area(geometry)),
                   Unit = ifelse(is.na(Unit)==T, "1Pump", Unit),
                   area = ifelse(Unit == "2Pump", area*2, area))

## rasterize the unit using the adjusted areas as the values for pixels
UnitSize <- rasterize(vect(Br_Units), LC_Br, field = "area", fun = "max",  background = NA)
plot(UnitSize)

## now inverse scale the raster so the smallest unit has a grading of 1
values(UnitSize) <- Inv_scale_vals(as.numeric(values(UnitSize)))
values(UnitSize) <- ifelse(is.na(values(UnitSize))==T, 0, values(UnitSize))
plot(UnitSize)
rm(Br_Units); gc() # free up space




##*-- Target hydro logical units that have few landowners --*##

## Read in the RPA point parcels for the Somerset regions
PPoints <- st_read("RawData/RPA/rpa_parcel_points_East_England/rpa_parcel_pointsPoint.shp")

## Read in the RPA point parcels for the Somerset regions
PointsID <- fread("RawData/RPA/RPA_ParcelList_20240424&Suffolk.csv")

## Join together the field parcel points and the customer data from the RPA
Join <- inner_join(PointsID, PPoints, by = join_by(NGC==parcel_ref)) |> st_as_sf() |> select(FARM_ID)

## Read in the shapevfile of Hydro Units 
Br_Units <- st_read("RawData/HydroUnits/Broads/Broads_HydroUnits.shp")
Br_Units$Site <- 1:nrow(Br_Units)
Br_Units <- Br_Units |> select(Site)

## Extract to label the parcel points with the hydro unit and then assign these back to the parcel points data set
TestExtr <- extract(vect(Br_Units), vect(Join))
Join$Site <- NA
Join$Site[TestExtr$id.y] <- TestExtr$Site

## Now summarise the number of different farm customers for each of the hydro units
UnitLandowners <- Join |> 
                  filter(is.na(Site) == F) |> 
                  group_by(Site) |> 
                  dplyr::summarise(No_farms = length(unique(FARM_ID))) |> st_drop_geometry()

## Calculate the density of landowners within each hydro unit
Br_Units <- left_join(Br_Units, UnitLandowners, by = "Site") |> 
             mutate(Area = st_area(geometry),
                    LandownerDens = No_farms/as.numeric(Area/1000000)) 

## Rasterize the hydro unit data using the land owner density values
LandownerDens <- rasterize(vect(Br_Units), LC_Br, field = "LandownerDens", fun = "max",  background = NA)

## Inverse scale the values so the units with the lowest landowner denisty have a score closer to 1
values(LandownerDens) <- Inv_scale_vals(values(LandownerDens))
values(LandownerDens) <- ifelse(is.na(values(LandownerDens))==T, 0, values(LandownerDens))
plot(LandownerDens)  
rm(PPoints, PointsID, Join, TestExtr, UnitLandowners, Br_Units); gc()



##*-- Target hydro units with natural variation in topography --*##

## 2m LiDAR raster for the Somerset
Br_Elev <- rast("RawData/LiDAR/Broads_DTM_2m.tif")
plot(Br_Elev)
## Shapefile of Hydro Units 
Br_Units <- st_read("RawData/HydroUnits/Broads/Broads_HydroUnits.shp")
 
## Add values to Somerset hydro units data set
Br_Units$Region <- "Broads"
Br_Units$id <- 1:nrow(Br_Units)
plot(Br_Units$geometry)

## Calculate the standard deviation in elevation for each hydro unit
ElevSD <- terra::extract(x=Br_Elev, y=vect(Br_Units), fun = sd, na.rm= T)

## assign the standard deviation back to the hydro unit polygons
Br_Units$ElevationSD <- ElevSD$file41d44237428b

## Rasterize the hydro unit data using the land owner density values
ElevSD <- rasterize(vect(Br_Units), LC_Br, field = "ElevationSD", fun = "min",  background = NA)

## Inverse scale the values so the units with the lowest landowner denisty have a score closer to 1
values(ElevSD) <- scale_vals(values(ElevSD))
values(ElevSD) <- ifelse(is.na(values(ElevSD))==T, 0, values(ElevSD))
plot(ElevSD)  
rm(Br_Elev, Br_Units); gc()



##*-- Target SSSI areas not currently notified for waders --*##

## Read in file that summarises if SSSI in the Broads are notified specifcally for breeding waders or not
Citation <- read.csv("RawData/Broads/SSSI_Citations.csv")
Citation <- Citation |> filter(Waders == "N")

## Read in a shapefile of all SSSI that I can filter for SSSI without a wader citation in the Broads
SSSI <- st_read("RawData/SSSI/data/Sites_of_Special_Scientific_Interest_England.shp")
SSSI <- SSSI |> filter(sssi_name %in% Citation$Name)
plot(SSSI$geometry)


## Turn the SSSI layer into raster
SSSI_NoNotif <- rasterize(vect(SSSI), LC_Br, cover = T, background = 0)
## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1, using my function
SSSI_NoNotif <- Half_reclass(SSSI_NoNotif)
plot(SSSI_NoNotif)

## free up spcae
rm(SSSI, Citation)



##*-- Target areas where more water will be available in the future --*##

## This uses a map that Andrea Kelly sent me, This map was produced by Figure 8 but I do not have a link to the finished report
## Read in the shapefile that I digitised
WaterAvail <- st_read("RawData/Broads/Surface Availability for Recharge.shp")
plot(WaterAvail)

## Assign values based on the colours of the map
## Green area will get a value of 1 (Water available)
## Yellow areas will be 0.5 (restricted water available)
## Red areas will be 0 (water not available)
WaterAvail <- mutate(WaterAvail, Value = ifelse(Avail=="Red", 0, ifelse(Avail=="Yellow", 0.5, 1)))
## rasterize the shapefile
WaterAvail <- rasterize(vect(WaterAvail), LC_Br, field = "Value", fun = "max",  background = 0)
plot(WaterAvail)







##----------------------##
#### 7. Mask Creation ####
##----------------------##

##*--Avoid priority habitats --*##

## Read in the NE Priority Habitat
NE <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")

## Crop out the Natural England priority habitats for the buffered priority landscape
NE <- st_intersection(NE, Br)

## Get a list of the main habitats and additional habitats
table(NE$mainhabs)
table(NE$addhabs)

## Create vector of habitats that I want to mask out
MAINHABS <- c("Coastal saltmarsh", 
              "Coastal sand dunes", "Coastal sand dunes,Deciduous woodland", "Coastal sand dunes,Reedbeds", 
              "Deciduous woodland", 
              "Lowland dry acid grassland",
              "Lowland fens", "Lowland fens,Reedbeds",
              "Lowland heathland",
              "Mudflats",
              "Reedbeds", "Reedbeds,Coastal saltmarsh",
              "Traditional orchard")

## Create column of just the first 5 characters of the additional habitats column
NE <- NE |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))

## Filter the data set by Main habitat and then if no main habitat retain those that did not have CFPGM in the additional habitat column
NE_Area_Masks <- filter(NE, mainhabs %in% MAINHABS |
                           (mainhabs == "No main habitat but additional habitats present" & !addhabs %in% c("CFPGM", "GQSIG", "LMEAD")))

## plot these masks
ggplot() + geom_sf(data = NE_Area_Masks, mapping = aes(geometry = geometry, fill = mainhabs), colour = NA) + theme_light()

## Rasterize the priotity habitat polygons
NE_Mask <- rasterize(vect(NE_Area_Masks), LC_Br, cover = T, background = 0)
## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1, using my function
NE_Mask <- Half_reclass(NE_Mask)
plot(NE_Mask)
rm(NE, MAINHABS, NE_Area_Masks); gc() # free up space



##*--Avoid scheduled monuments --*##

## Read in polygons of scheduled monuments from National Heritage
Monument <- st_read("RawData/Monuments/National_Heritage_List_for_England/Scheduled_Monuments.shp")

## Crop out the monuments for the buffered priority landscape and buffer then my 20m
Monument <- st_intersection(Monument, Br) |> st_buffer(dist = 20)
plot(Monument$geometry)

## Rasterize the monument polygons
Monument_Mask <- rasterize(vect(Monument), LC_Br, cover = T, background = 0)
## If a pixel is more than 50% covered by a monument then assign a pixel value of 1, using my function
Monument_Mask <- Half_reclass(Monument_Mask)
plot(Monument_Mask)
rm(Monument); gc()



##*--Avoid urban and coastal habitats --*##

## Reclassify raster so that it indicates if pixel is Urban/Sub-urban/Coastal or not
m <- c(0.0, 14.5, 0,
       14.5, 18.5, 1,
       18.5, 19.5, 0, 
       19.5, 22.0, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Urban_mask <- classify(LC_Br, rclmat, include.lowest=TRUE)
plot(Urban_mask)





##-----------------------------------------##
##-----------------------------------------##
#### 8. Define strategy opportunity area ####
##-----------------------------------------##
##-----------------------------------------##



##------------------------------------##
#### 9. Group 1: Combine Guidelines ####
##------------------------------------##

##------------------------##
## Group 1: Combine masks ##
##------------------------##

## Take max value from all mask layers
G1_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
G1_mask <- classify(G1_mask, rclmat, include.lowest=TRUE)
plot(G1_mask)
writeRaster(G1_mask, "CleanData/Guideline Creation/Norfolk/Broads_MasksAll_G1.tif", overwrite=TRUE)

## create a general mask to use in Scenarios later
Gen_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Gen_mask <- classify(Gen_mask, rclmat, include.lowest=TRUE)
plot(Gen_mask)
writeRaster(Gen_mask, "CleanData/Guideline Creation/Norfolk/Broads_MasksAll_Gen.tif", overwrite=TRUE)



##------------------------------------##
## Group 1: Combine better guidelines ##
##------------------------------------##

Bett_G1 <- LowGrassDens + Hydro_HeightLow + UnitSize + LandownerDens + ElevSD + WaterAvail
names(Bett_G1) <- "Better"

ggplot() + geom_spatraster(data=Bett_G1) + labs(fill = "Grading") + scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G1, "CleanData/Guideline Creation/Norfolk/Broads_Better_G1.tif", overwrite=TRUE)


##------------------------------------##
## Group 1: Combine bigger guidelines ##
##------------------------------------##

Big_G1 <- SiltSoil + WoodDens + LowGrassDens + AESWint + Hydro_HeightLow + UnitSize + LandownerDens + ElevSD + WaterAvail + Connectr
names(Big_G1) <- "Bigger"

ggplot() + geom_spatraster(data=Big_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G1, "CleanData/Guideline Creation/Norfolk/Broads_Bigger_G1.tif", overwrite=TRUE) # read out grading



##----------------------------------##
## Group 1: Combine more guidelines ##
##----------------------------------##

More_G1 <- SiltSoil + WoodDens + LowGrassDens + AESWint + Hydro_HeightLow + UnitSize + LandownerDens + ElevSD + WaterAvail + Connectr
names(More_G1) <- "More"

ggplot() + geom_spatraster(data=More_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G1, "CleanData/Guideline Creation/Norfolk/Broads_More_G1.tif", overwrite=TRUE) # read out grading



##---------------------------------------------##
## Group 1: Combine agri-conversion guidelines ##
##---------------------------------------------##

ArableConv_G1 <- Hydro_HeightLow + LowGrassDens + WaderSitesDist + HistoricGrass
names(ArableConv_G1) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableConv_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableConv_G1, "CleanData/Guideline Creation/Norfolk/Broads_ArableRev_G1.tif", overwrite=TRUE)




##-------------------------------------##
#### 10. Group 2: Combine Guidelines ####
##-------------------------------------##

##------------------------##
## Group 2: Combine masks ##
##------------------------##

## Take max value from all mask layers
G2_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
G2_mask <- classify(G2_mask, rclmat, include.lowest=TRUE)
plot(G2_mask)
writeRaster(G2_mask, "CleanData/Guideline Creation/Norfolk/Broads_MasksAll_G2.tif", overwrite=TRUE)


##------------------------------------##
## Group 2: Combine better guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Bett_G2 <- Hydro_HeightLow + SSSI_NoNotif + WaterAvail
names(Bett_G2) <- "Better"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Bett_G2) + labs(fill = "Grading") + scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G2, "CleanData/Guideline Creation/Norfolk/Broads_Better_G2.tif", overwrite=TRUE)



##------------------------------------##
## Group 2: Combine bigger guidelines ##
##------------------------------------##

Big_G2 <- SiltSoil + WoodDens + UrbanDens + Hydro_HeightLow + LandownerDens + ExpandClust + SSSI_NoNotif + WaterAvail
names(Big_G2) <- "Bigger"

ggplot() + geom_spatraster(data=Big_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G2, "CleanData/Guideline Creation/Norfolk/Broads_Bigger_G2.tif", overwrite=TRUE) # read out grading



##----------------------------------##
## Group 2: Combine more guidelines ##
##----------------------------------##

More_G2 <- SiltSoil + WoodDens + UrbanDens + Hydro_HeightLow + LandownerDens + SSSI_NoNotif + WaterAvail
names(More_G2) <- "More"

ggplot() + geom_spatraster(data= More_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G2, "CleanData/Guideline Creation/Norfolk/Broads_More_G2.tif", overwrite=TRUE) # read out grading



##----------------------------------------------##
## Group 2: Combine agri-conversion guidelines ##
##---------------------------------------------##

## Add together all rules and assign name to layer
ArableConv_G2 <- Hydro_HeightLow + AgriGrades + ExpandClust
names(ArableConv_G2) <- "ArableConv"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=ArableConv_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableConv_G2, "CleanData/Guideline Creation/Norfolk/Broads_ArableRev_G2.tif", overwrite=TRUE)





##-------------------------------------------------##
#### 11.0 Plot graded maps with opportunity areas ####
##-------------------------------------------------##

##---------------------------------##
## 11.1 Prepare Landscape Canvas ####
##---------------------------------##


##-- Label lowland wet grassland fields --##

## Read in the Broads canvas with all land parcel polygons
CanvasGr <- st_read("CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp")
plot(CanvasGr$geometry)

## Read in a raster icreaetd of lowland wet grassland extent in the Broads
LowGrass <- rast("RawData/LowlandWetGrassRasters/Broads_LowWetGrass.tif")
plot(LowGrass)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(LowGrass, vect(CanvasGr), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 25% covered by lowland wet grassland pixels
CanvasGr <- CanvasGr |>
          mutate(WetGrassProp = Scores$layer) |>
          filter(WetGrassProp > 0.50) |> select(-WetGrassProp)
plot(CanvasGr$geometry)
rm(Scores); gc()

## Crop the canvas to just the piority landscape
CanvasGr <- CanvasGr |> st_intersection(BrOutline)
plot(CanvasGr$geometry)



##-- Label arable opportunity fields --##

## Read in the Broads canvas with all land parcel polygons
CanvasAr <- st_read("CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp")
plot(CanvasAr$geometry)

## Read in a raster I created of lowland wet grassland extent in the Broads
SuitArable <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_ArableSuitable.tif")
plot(SuitArable)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(SuitArable, vect(CanvasAr), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 50% covered by lowland wet grassland pixels
CanvasAr <- CanvasAr |>
          mutate(SuitArabProp = Scores[,2]) |>
          filter(SuitArabProp > 0.50) |> select(-SuitArabProp)
plot(CanvasAr$geometry)
rm(Scores); gc()

## Crop the canvas to just the priority landscape
CanvasAr <- CanvasAr |> st_intersection(BrOutline)
plot(CanvasAr$geometry)



##-- Extract the grades/masks for each of the Lawton strategies --##

## Code to read in the rules if needed
G1_mask <- rast("CleanData/Guideline Creation/Norfolk/Broads_MasksAll_G1.tif")
Bett_G1 <- rast("CleanData/Guideline Creation/Norfolk/Broads_Better_G1.tif")
Big_G1 <- rast("CleanData/Guideline Creation/Norfolk/Broads_Bigger_G1.tif")
More_G1 <- rast("CleanData/Guideline Creation/Norfolk/Broads_More_G1.tif")
ArableConv_G1 <- rast("CleanData/Guideline Creation/Norfolk/Broads_ArableRev_G1.tif")

G2_mask <- rast("CleanData/Guideline Creation/Norfolk/Broads_MasksAll_G2.tif")
Bett_G2 <- rast("CleanData/Guideline Creation/Norfolk/Broads_Better_G2.tif")
Big_G2 <- rast("CleanData/Guideline Creation/Norfolk/Broads_Bigger_G2.tif")
More_G2 <- rast("CleanData/Guideline Creation/Norfolk/Broads_More_G2.tif")
ArableConv_G2 <- rast("CleanData/Guideline Creation/Norfolk/Broads_ArableRev_G2.tif")


## More strategy
MoreGrades <- extract(More_G1, CanvasGr, fun = mean, na.rm = T)
MoreGrades2 <- extract(More_G2, CanvasGr, fun = mean, na.rm = T)
CanvasGr$MoreGrade_G1 <- MoreGrades$More
CanvasGr$MoreGrade_G2 <- MoreGrades2$More

## Bigger strategy
BiggerGrades <- extract(Big_G1, CanvasGr, fun = mean, na.rm = T)
BiggerGrades2 <- extract(Big_G2, CanvasGr, fun = mean, na.rm = T)
CanvasGr$BigGrade_G1 <- BiggerGrades$Bigger
CanvasGr$BigGrade_G2 <- BiggerGrades2$Bigger

## Better Strategy
BetterGrades <- extract(Bett_G1, CanvasGr, fun = mean, na.rm = T)
BetterGrades2 <- extract(Bett_G2, CanvasGr, fun = mean, na.rm = T)
CanvasGr$BetterGrade_G1 <- BetterGrades$Better
CanvasGr$BetterGrade_G2 <- BetterGrades2$Better

## Arable reversion Strategy
ArableGrades <- extract(ArableConv_G1, CanvasAr, fun = mean, na.rm = T)
ArableGrades2 <- extract(ArableConv_G2, CanvasAr, fun = mean, na.rm = T)
CanvasAr$ArableGrade_G1 <- ArableGrades$ArableConv
CanvasAr$ArableGrade_G2 <- ArableGrades2$ArableConv


## Add on cluster size to the Better grading
PopSize <- PopSize |> select(ClustPop, F_LOC_ID) |> mutate(ClustPop = Inv_scale_vals(ClustPop)) # extract population cluster sizes
CanvasGr <- left_join(CanvasGr, PopSize, by = "F_LOC_ID") # add these pop sizes to the fields (if they are in a cluster)
CanvasAr <- left_join(CanvasAr, PopSize, by = "F_LOC_ID") # add these pop sizes to the fields (if they are in a cluster)
## Do the final addition to work out grade
CanvasGr <- CanvasGr |> mutate(ClustPop2 = ifelse(is.na(ClustPop)==T, 0, ClustPop), 
                           BetterGrade_G1 = BetterGrade_G1+ClustPop2,
                           BetterGrade_G2 = BetterGrade_G2+ClustPop2) |> 
                    select(-ClustPop2) # Add inverse scaled clust pop size to the better grading 

                    
## Which fields need to be masked
## First re-label raster to pixels to be masked have value of 0
Masks <- G1_mask
values(Masks) <- ifelse(is.na(values(Masks))==T, 0, values(Masks))
Masks2 <- G2_mask
values(Masks2) <- ifelse(is.na(values(Masks2))==T, 0, values(Masks2))

## Now extract the average value for each land parcel, fully masked fields will have a score of 1
OppScore <- extract(Masks, CanvasGr, fun = mean, na.rm = T)
OppScore2 <- extract(Masks2, CanvasGr, fun = mean, na.rm = T)
CanvasGr$Mask_G1 <- OppScore$layer
CanvasGr$Mask_G2 <- OppScore2$layer

## re do this for the arable canvas
OppScore <- extract(Masks, CanvasAr, fun = mean, na.rm = T)
OppScore2 <- extract(Masks2, CanvasAr, fun = mean, na.rm = T)
CanvasAr$Mask_G1 <- OppScore$layer
CanvasAr$Mask_G2 <- OppScore2$layer



##-----------------------------##
## 11.2 General plot styling ####
##-----------------------------##

## Set a general theme so that all pots look similar
GeneralThemeing <- theme(
          axis.title.x = element_text(size = 17),
          #axis.text.x = element_text(hjust=0.7, angle = 45),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 17),
          axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
          text = element_text(color = "#2D2D2E"), 
          panel.grid = element_line(color = "#ebebe5", linewidth = 0.2),
          panel.background = element_rect(fill = "#f5f5f2", color = NA))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(CanvasGr)[1]-20, ext(CanvasGr)[2]+20), ylim = c(ext(CanvasGr)[3]-20, ext(CanvasGr)[4]+20), 
                    crs = 27700, expand = FALSE) 



##-------------------------------------##
## 11.3 Group 1: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasGrBetter <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Broads: Better Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBetter)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G1_Better.png", width = 20, height = 20, units = "cm")



##-------------------------------------##
## 11.4 Group 1: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Broads: Bigger Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBig)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G1_Bigger.png", width = 20, height = 20, units = "cm")



##-----------------------------------##
## 11.5 Group 1: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Broads: More Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrMore)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G1_More.png", width = 20, height = 20, units = "cm")




##-----------------------------------------------##
## 11.6 Group 1: Plot arable reversion grading ####
##-----------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG1 <- filter(CanvasAr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG1, mapping=aes(geometry=geometry, fill = ArableGrade_G1), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Broads: Arable Reversion G1") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG1)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G1_ArableRev.png", width = 20, height = 20, units = "cm")




##---------------------------------------##
## 11.7 Group 1: Plot Map of Landscape ####
##---------------------------------------##

## Read in priority landscape boundary
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
Broads <- MyBoxes |> filter(Att3Value == "Broads") 

## Read in RSPB reserve outlines
Reserves <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")
ReservesBr <- st_crop(Reserves, (Br |> st_buffer(dist=500)))

## make plot
ggplot() + 
  ## add polygons
  geom_sf(data=Broads, mapping=aes(geometry=geometry), fill = "lightblue", colour = NA) + 
  geom_sf(data=ReservesBr, mapping=aes(geometry=geometry), fill = "red", alpha = 0.5, colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## Add North arrow and scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  ## set labels
  ggtitle("Norfolk Broads") +
  ## set them
  theme_light() + 
  GeneralThemeing +
  theme(legend.position = "none")

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_LandscapeMap.png", width = 20, height = 20, units = "cm")




##-------------------------------------##
## 11.8 Group 2: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasGrBetter <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Broads: Better Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBetter)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G2_Better.png", width = 20, height = 20, units = "cm")
      


##-------------------------------------##
## 11.9 Group 2: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Broads: Bigger Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBig)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G2_Bigger.png", width = 20, height = 20, units = "cm")



##-----------------------------------##
## 11.10 Group 2: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Broads: More Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrMore)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G2_More.png", width = 20, height = 20, units = "cm")



##------------------------------------------------##
## 11.11 Group 2: Plot arable reversion grading ####
##------------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG2 <- filter(CanvasAr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=BrOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG2, mapping=aes(geometry=geometry, fill = ArableGrade_G2), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Norfolk Broads: Arable Reversion G2") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG2)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Broads_G2_ArableRev.png", width = 20, height = 20, units = "cm")



