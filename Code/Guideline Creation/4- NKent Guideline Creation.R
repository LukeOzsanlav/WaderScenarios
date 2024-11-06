##----------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 09/05/2024
##
## Aim: Create all guidelines for scenario modelling of the North Kent Landscape
## 
##----------------------------------------------------------## 


## Load in packages
pacman::p_load(tidyverse, data.table, sf, terra, tidyterra, leastcostpath, ggspatial)
options(scipen=999) # turn off scientific notation

## Load in various helper functions
source("Code/Helper functions.R")



##----------------------------------##
#### 0. Read in general data sets ####
##----------------------------------##

## Read in my priority landscapes boundary and buffer it by 3000m
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
NKOutline <- MyBoxes |> filter(Att3Value == "North Kent")
NK <- MyBoxes |> filter(Att3Value == "North Kent") |> st_buffer(dist=3000)
rm(MyBoxes) # free up space

## Read in UKCEH landcover and then crop and mask it with the landscape boundary
## This will be used as the canvas to map all other variables across
LCfull <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]]
LC_NK <- crop(LCfull, vect(NK)) |> mask(vect(NK))
rm(LCfull)
plot(LC_NK) # free up space




##-------------------------------##
#### 1. Better only guidelines ####
##-------------------------------##

##*-- Target smallest existing populations --*##

## Read in full BWWM data set that has each wader cluster labelled
Reg_Fields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Kent") 
InterFields <- Reg_Fields |> st_intersection(NKOutline)
Reg_Fields <- Reg_Fields |> filter(F_LOC_ID %in% InterFields$F_LOC_ID)

## Read in filtered breeding pairs estimates
## And calculate total wader abundance for each field
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
Waders <- Waders |> mutate(Tot_abund = rowSums(dplyr::select(Waders, est_pairsL, est_pairsR, est_pairsS), na.rm = T)) |> 
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
write_csv(PopSize, "CleanData/Guideline Creation/Kent/Kent_ClusterPopSize.csv")



##*-- Avoid RSPB reserves and Elmley --*##

## filter out the wader sites with more than 50 pairs
LargeSites <- PopSize |> filter(ClustPop > 50) |> st_as_sf()

## RSPB reserves, crop to large wader sites to see which RSPB reserves I want
RSPB <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp") 
Inter <- RSPB |> st_intersection(LargeSites)
#plot(Inter$geometry); table(Inter$Att3Value)
RSPB <- filter(RSPB, Att3Value %in% c("CLIFFE POOLS", "HIGHAM MARSHES", "NORTHWARD HILL", "SHORNE MARSHES", "GREAT BELLS FARM")) |> select(geometry)

## National Nature Reserves, crop to large wader sites
NNR <- st_read("RawData/Other Reserves/National_Nature_Reserves_EnglandPolygon.shp")
Inter <- NNR |> st_intersection(LargeSites)
#plot(Inter$geometry); table(Inter$nnr_name)
NNR <- filter(NNR, nnr_name == "Elmley") |> select(geometry)

## Rasterize the selected RSPB reserves and NNR's
Shapes <- rbind(RSPB, NNR) |> st_buffer(dist =10)
AvoidRSPB <- rasterize(vect(Shapes), LC_NK, cover = T, background = 0)
values(AvoidRSPB) <- ifelse(values(AvoidRSPB)==0, 1, 0)
plot(AvoidRSPB)





##-------------------------------##
#### 2. Bigger only guidelines ####
##-------------------------------##







##-----------------------------##
#### 3. More only guidelines ####
##-----------------------------##

#*-- Target areas that link most isolated sites to other sites --*##

## Read in polygons of the wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Kent/Kent_Wader_Sites.shp")
WaderSites <- WaderSites |> st_crop(NK) |> st_centroid()
plot(WaderSites$geometry)

## Read in raster of potential lowland wet grassland and suitable arable land
## Then combine so have tiered values
LowGrass <- rast("RawData/LowlandWetGrassRasters/NKent_LowWetGrass.tif")
ArableSuit <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_ArableSuitable.tif") |> resample(LowGrass)
LowGrass <- max(LowGrass*4, ArableSuit*2) + 1
plot(LowGrass)

## Creates a cost surface using the values in the supplied SpatRaster.
Cond <- create_cs(LowGrass, neighbours = 16)
plot(Cond)

## Calculates Least-cost paths from-everywhere-to-everywhere. This is based on the approach proposed by White and Barber (2012).
LCP <- create_FETE_lcps(x=Cond, locations=WaderSites)
plot(LCP)
LCP$Val <- 1

## For each wader cluster node
## Calcualte the total path distacen of all connection and then join these distacnes
## back onto all of the least cost paths
LCP_sum <- LCP |> mutate(Path_length = st_length(geometry)) |> 
                  group_by(origin_ID) |> 
                  summarise(TotLength= sum(as.numeric(Path_length))) |> 
                  st_drop_geometry() |> 
                  left_join(LCP, by = "origin_ID") |> 
                  st_as_sf()

## No turn the lines that connect polygons into a raster
## Create a base raster with a pixel size of 2km using my base land cover raster
LC_NK2km <- aggregate(LC_NK, fact = 80)

## Rasterize the lines connecting polygons using my 2km raster as a base
IsoConnectr <- rasterize(vect(LCP_sum), LC_NK2km, field = "TotLength", fun = "sum",  background = 0)
ggplot() + geom_spatraster(data = IsoConnectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)

## Scale the values from the rasterized lines
values(IsoConnectr) <- scale_vals(values(IsoConnectr))
plot(IsoConnectr)

## Finally resample the raster back down to 25m
IsoConnectr <- resample(IsoConnectr, LC_NK)

## Finally turn any NAs back to zeros
values(IsoConnectr) <- ifelse(is.na(values(IsoConnectr))==T, 0, values(IsoConnectr))
plot(IsoConnectr)

## free up space
rm(WaderSites, LowGrass, Cond, LCP, LC_NK2km); gc()




##-------------------------------##
#### 4. Bigger/More guidelines ####
##-------------------------------##

##*-- Target areas near/far away sites with more breeding waders --*##

## The sentiment in the wokrshop was to expand the larger sites
## I am going to use a cut off of 50 pairs to classify a site as having more then 50 pairs
LargeSites <- PopSize |> filter(ClustPop > 50)
ggplot(data= LargeSites) + geom_sf(mapping = aes(geometry=geometry, fill = ClustPop), colour = NA) + theme_light()

## Rasterize the large sites shape file
LargeSites <- rasterize(vect(st_as_sf(LargeSites)), LC_NK, cover = T, background = NA)
plot(LargeSites)

## Calculate distances to these large sites and then mask to the priority landscape
SiteDist <- distance(LargeSites) |> mask(NKOutline)
plot(SiteDist)

## Finally take the inverse scale of the values and turn all NA's into 0's
InvSiteDist <- SiteDist
values(InvSiteDist) <- Inv_scale_vals(values(InvSiteDist))
values(InvSiteDist) <- ifelse(is.na(values(InvSiteDist))==T, 0, values(InvSiteDist))
plot(InvSiteDist)

## Finally take the inverse scale of the values and turn all NA's into 0's
values(SiteDist) <- scale_vals(values(SiteDist))
values(SiteDist) <- ifelse(is.na(values(SiteDist))==T, 0, values(SiteDist))
plot(SiteDist)


## free up space
rm(LargeSites); gc()



##*-- Target areas surrounded by less urban areas --*##

## Reclassify land cover raster so that it indicates if pixel is woodland or not
LC_NKWid <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]] |>  crop(vect(NK))
m <- c(0.0, 19.5, 0,
       19.5, 23, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Urban <- classify(LC_NKWid, rclmat, include.lowest=TRUE)
plot(Urban)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
## then calculate inverse scaled values as we want areas with lower tree cover
UrbanDens <- focal(Urban, w=41, fun="mean", na.rm = T)
values(UrbanDens) <- Inv_scale_vals(values(UrbanDens))
values(UrbanDens) <- ifelse(is.na(values(UrbanDens))==T, 0, values(UrbanDens))
plot(UrbanDens) # plot smoothed raster
rm(LC_NKWid, Urban); gc() # free up space



##*-- Target areas with wintering waterbird CS agreements already --*##

## Read in the Different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")
ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESS <- rbind(ESSWest, ESSEast)
rm(ESSWest, ESSEast)

## Read in the Canvas as this has the field parcel shapes 
Canv <- st_read("CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp") |> select(geometry)

## Filter out winter water birds focused options and breeding waders water options
ESS_OpsWa <- filter(ESS, optcode %in% c("HK9", "HK11", "HK13", "HK10", "HK12", "HK14"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
ESS_InterWa <- as.data.frame(st_covers(Canv, ESS_OpsWa))
Canv$ESS_Wader <- "N"
Canv$ESS_Wader[ESS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()


## Filter out winter waterbird or breeding focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
CSS_InterWa <- as.data.frame(st_covers(Canv, CSS_OpsWa))
Canv$CSS_Wader <- "N"
Canv$CSS_Wader[CSS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = CSS_Wader), colour = NA) +
           theme_minimal()

## Rasterize this AES field polygons
Canv <- filter(Canv, CSS_Wader == "Y" | ESS_Wader == "Y")
AESWader <- rasterize(vect(Canv), LC_NK, cover = T, background = 0)
AESWader <- Half_reclass(AESWader)
plot(AESWader, main = "Wader AES fields")

## free up space
rm(Canv, ESS, ESS_InterWa, CSS_InterWa, CSS_Ops, CSS_OpsWa, ESS_OpsWa); gc()





##-------------------------------------##
#### 5. Arable conversion guidelines ####
##-------------------------------------##

##*-- None --*##






##--------------------------------##
#### 6. All strategy guidelines ####
##--------------------------------##

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
AgriGrades <- rasterize(vect(AgriGrade), LC_NK, field = "Grade", fun = "min", background = 0)
values(AgriGrades) <- scale_vals(values(AgriGrades))
plot(AgriGrades)
rm(AgriGrade); gc()



##*-- Target areas less likely to be targeted for salt marsh creation --*##

## Read in the two shape file layer, one for opportunities and one for 
SaltOpp <- st_read("RawData/Saltmarsh Creation Sites/SustainableShores_all_potential_sites.shp") |> st_transform(crs = st_crs(LC_NK))
SaltPr <- st_read("RawData/Saltmarsh Creation Sites/SustainableShores_priority_sites.shp") |> st_transform(crs = st_crs(LC_NK))

## rasterize the opportunity areas with a value of 0.5
SaltOpp$Val <- 0.5
SaltOpp <- rasterize(vect(SaltOpp), LC_NK, field = "Val", fun = "max",  background = 1)
plot(SaltOpp)

## rasterize the opportunity areas with a value of 1
SaltPr$Val <- 0
SaltPr <- rasterize(vect(SaltPr), LC_NK, field = "Val", fun = "max",  background = 1)
plot(SaltPr)

## Take the max value of the two rasters
SaltmarshOpp <- min(SaltOpp, SaltPr)
plot(SaltmarshOpp)

## free up space
rm(SaltOpp, SaltPr); gc()



##*-- Target areas that link existing wader sites --*##

## Read in polygons of the wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Kent/Kent_Wader_Sites.shp")
WaderSites <- WaderSites |> st_crop(NK) |> st_centroid()
plot(WaderSites$geometry)

## Read in raster of potential lowland wet grassland and suitable arable land
## Then combine so have tiered values
LowGrass <- rast("RawData/LowlandWetGrassRasters/NKent_LowWetGrass.tif")
ArableSuit <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_ArableSuitable.tif") |> resample(LowGrass)
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
LC_NK2km <- aggregate(LC_NK, fact = 80)

## Rasterize the lines connecting polygons using my 2km raster as a base
Connectr <- rasterize(vect(LCP), LC_NK2km, field = "Val", fun = "sum",  background = 0)
ggplot() + geom_spatraster(data = Connectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)

## Scale the values from the rasterized lines
values(Connectr) <- scale_vals(values(Connectr))
plot(Connectr)

## Finally resample the raster back down to 25m
Connectr <- resample(Connectr, LC_NK)

## Finally turn any NAs back to zeros
values(Connectr) <- ifelse(is.na(values(Connectr))==T, 0, values(Connectr))
plot(Connectr)

## free up space
rm(WaderSites, LowGrass, Cond, LCP, LC_NK2km); gc()



##*--Target areas that were originally grassland --*##

## 2007 data set
LC2007 <- rast("RawData/LandCover/lcm2007gb25m.tif")[[1]]

## crop the UK CEH data to polygon of interest
LC2007 <- crop(LC2007, vect(st_buffer(st_as_sf(NK), dist = 500)))
plot(LC2007)


## 2000 data set, need to re download this data set as it will not read in
LC2000 <- rast("RawData/LandCover/LCM2000_GB_25m.tif")[[1]]

## crop the UK CEH data to polygon of interest
LC2000 <- crop(LC2000, vect(st_buffer(st_as_sf(NK), dist = 500)))
plot(LC2000)


## 1990 data set
LC1990<- rast("RawData/LandCover/gb1990lcm25m.tif")[[1]]

## crop the UK CEH data to polygon of interest
LC1990 <- crop(LC1990, vect(st_buffer(st_as_sf(NK), dist = 500)))
plot(LC1990)


## 2007 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 3.5, 0,
       3.5, 9.5, 1,
       9.5, 23, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC2007_Grass <- classify(LC2007, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC2007_Grass)


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
All_Grass <- (LC1990_Grass+LC2000_Grass+LC2007_Grass)/3
HistoricGrass <- resample(All_Grass, LC_NK) # resample so that the extent is the same as all my other grading layers
plot(HistoricGrass)
rm(All_Grass, LC1990_Grass, LC2000_Grass, LC1990, LC2000); gc()



##*-- Avoid areas near large heronries --*##

## Read in the shapefile of Heronries
Heronries <- st_read("RawData/Thames/North Kent Heronries/Kent Heronries.shp") |> st_centroid()

## Buffer by 2000
## Could not find a strong reference for this 
HeronriesBuff <- Heronries |> st_buffer(dist = 2000)

## Rasterize the colonies
Heronries <- rasterize(vect(Heronries), LC_NK, cover = T, background = NA)

## calculate the distance of all pixles from the colony centre
## Then mask to just pixls within 2km of colony centre
HeronriesDist <- distance(Heronries) |> mask(HeronriesBuff)
plot(HeronriesDist)

## scale the values and give all NAs a value of 1
values(HeronriesDist) <- scale_vals(values(HeronriesDist))
values(HeronriesDist) <- ifelse(is.na(values(HeronriesDist))==T, 1, values(HeronriesDist))
plot(HeronriesDist)



##*-- Target areas with more water in future --*##

## Read in data set on water abstraction availability
WaterAbst <- st_read("RawData/Thames/Water Availability/Resource_Reliability_pct_of_timePolygon.shp")
WaterAbst$ID <- 1:nrow(WaterAbst) # giev unique row ID

## Filter out the polygons that at least partly fall within the priority landscape
Sub <- st_intersection(WaterAbst, NK)
WaterAbst <- WaterAbst |> filter(ID %in% Sub$ID)
ggplot() + geom_sf(data = WaterAbst, mapping=aes(geometry=geometry, fill = resavail)) # quick plot

## Change the character % time water is available to a number, this matches the conversion done for Somerset
WaterAbst <- WaterAbst |> 
  mutate(PctVal = case_when(resavail == "less than 30%" ~ 15,
                            resavail == "at least 30%" ~ 30,
                            resavail == "at least 50%" ~ 50,
                            resavail == "at least 70%" ~ 70,
                            resavail == "at least 95%" ~ 95,
                            .default = 15))
ggplot() + geom_sf(data = WaterAbst, mapping=aes(geometry=geometry, fill = PctVal)) # quick plot


## Rasterize the polygons data set and give the pixels the values of the % of time water is available for abstraction
WaterAbstr <- rasterize(vect(WaterAbst), LC_NK, field = "PctVal", fun = "min",  background = 15)
## Scale values and change any NAs to zeros
values(WaterAbstr) <- scale_vals(values(WaterAbstr))
values(WaterAbstr) <- ifelse(values(WaterAbstr) == "NaN", 0, values(WaterAbstr))
plot(WaterAbstr)

## free up space
rm(WaterAbst); gc()



##*-- Target low height land within hydro units --*##

## Read in Thames Lidar then mask to the Kent region
Kent_Elev <- rast("RawData/LiDAR/Thames_DTM_2m.tif") |> crop(vect(NKOutline)) |> mask(vect(NKOutline))
## aggregate the Lidar up to 4mx4m resolution, this speeds up processing time
Kent_Elev <- aggregate(Kent_Elev, fact = 2)
plot(Kent_Elev)

## Read in a shape file of Hydro Units then extract the ones at least partly in the priority landscape
Kent_Units <- st_read("RawData/HydroUnits/Thames/Thames_Hydrog_Units.shp")
Kent_Units$id <- 1:nrow(Kent_Units)
Kent_Inter <- Kent_Units |> st_intersection(NKOutline)
Kent_Units <- Kent_Units |> filter(id %in% Kent_Inter$id)


## Loop through each hydro unit and extract the LiDAR data and grade it using a empirical cumulative distribution function
## create an empty list to put all the raster outputs into
List <- vector(mode='list', length= nrow(Kent_Units))

## Start the loop
for(j in 1:nrow(Kent_Units)){
  
  ## message to console
  message("iteration ", j)
  
  ## filter out the hydro unit outline
  Kent_Unit <- vect(Kent_Units[j,])
  
  ## Crop the hydro unit from the LiDAR data
  Unit_Lidar <- crop(Kent_Elev, Kent_Unit, mask = T)
  plot(Unit_Lidar)
  
  ## Calculate empirical cumulative distribution function
  x <- values(Unit_Lidar)
  y <- ecdf(x)(x)
  
  ## Take inverse of distribution function so lowest land has value of 1 and highest land has value of 1
  values(Unit_Lidar) <- Inv_scale_vals(y)
  plot(Unit_Lidar)
  
  ## add raster to slot in list
  List[[j]] <- Unit_Lidar
  
}

## Combine all the raster into a raster collection
List <- sprc(List)

## now merge that raster collection
Hydro_Height <- terra::merge(List)
plot(Hydro_Height, col=grDevices::hcl.colors(50, palette = "Sunset"))



## Now for all pixels that are in the priority landscape but not within a hydro unit
## We need to assign these values for height grading

## Firstly set all the pixels outside the priority landscape to 0 (from NA)
## Rasterize the outline of the priority landscape with a bit of a buffer around it
NKOutline2 <- NKOutline |> st_buffer(dist = 25) |> mutate(ID=1)
Landscape <- rasterize(vect(NKOutline2), Hydro_Height, field = "ID", fun = "max",  background = NA)
plot(Landscape)
## Change all NAs in the raster to zeros that are not within the priority landscape
values(Hydro_Height) <- ifelse(is.na(values(Hydro_Height))==T & is.na(values(Landscape))==T, 0, values(Hydro_Height))
plot(Hydro_Height, col=grDevices::hcl.colors(50, palette = "Sunset"))
gc()


## Now need to create an ecdf for the pixels that were in a hyrdo unit
## Create raster of all the hydro units
UnitRast <- rasterize(vect(Kent_Units), Hydro_Height, field = "id", fun = "max",  background = NA)
ElevVect <- ifelse(is.na(values(UnitRast))==F, values(Kent_Elev), NA)
ElevVect <- ElevVect[is.na(ElevVect)==FALSE]
func <- ecdf(ElevVect)


## Now fill in the remaining NAs with the graded values
## Do this by calcualting the inverse of there ecdf values
Hydro_Height <- resample(Hydro_Height, Kent_Elev)
values(Hydro_Height) <- ifelse(is.na(values(Hydro_Height))==T, Inv_scale_vals(func(values(Kent_Elev))), values(Hydro_Height))
plot(Hydro_Height, col=grDevices::hcl.colors(50, palette = "Sunset"))


## resample the fine-scale topography raster to my base raster so I can add it to all the other graded rasters
Hydro_Height <- resample(Hydro_Height, LC_NK)
values(Hydro_Height) <- ifelse(is.na(values(Hydro_Height))==T, 0, values(Hydro_Height))
plot(Hydro_Height, col=grDevices::hcl.colors(50, palette = "Sunset"))
rm(Kent_Elev, Kent_Units, Kent_Unit, List, Unit_Lidar, NKOutline2, Landscape, UnitRast, ElevVect, func); gc() # free up some space






##----------------------##
#### 7. Mask Creation ####
##----------------------##


##*--Avoid priority habitats --*##

## Read in the NE Priority Habitat
NE <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")

## Crop out the Natural England priority habitats for the buffered priority landscape
NE <- st_intersection(NE, NK)

## Get a list of the main habitats and additional habitats
table(NE$mainhabs)
table(NE$addhabs)

## Create vector of habitats that I want to mask out
MAINHABS <- c("Coastal saltmarsh",
              "Coastal vegetated shingle",
              "Deciduous woodland", "Deciduous woodland,Maritime cliff and slope",
              "Lowland calcareous grassland",
              "Maritime cliff and slope", "Maritime cliff and slope,Coastal saltmarsh",
              "Mudflats",
              "Reedbeds", "Reedbeds,Coastal saltmarsh",
              "Saline lagoons",
              "Traditional orchard")

## Create column of just the first 5 characters of the additional habitats column
NE <- NE |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))

## Filter the data set by Main habitat and then if no main habitat retain those that did not have CFPGM in the additional habitat column
NE_Area_Masks <- filter(NE, mainhabs %in% MAINHABS |
                           (mainhabs == "No main habitat but additional habitats present" & !addhabs %in% c("CFPGM", "GQSIG", "LMEAD")))

## plot these masks
ggplot() + geom_sf(data = NE_Area_Masks, mapping = aes(geometry = geometry, fill = mainhabs), colour = NA) + theme_light()

## Rasterize the priotity habitat polygons
NE_Mask <- rasterize(vect(NE_Area_Masks), LC_NK, cover = T, background = 0)
## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1, using my function
NE_Mask <- Half_reclass(NE_Mask)
plot(NE_Mask)
rm(NE, MAINHABS, NE_Area_Masks); gc() # free up space



##*--Avoid scheduled monuments --*##

## Read in polygons of scheduled monuments from National Heritage
Monument <- st_read("RawData/Monuments/National_Heritage_List_for_England/Scheduled_Monuments.shp")

## Crop out the monuments for the buffered priority landscape
Monument <- st_intersection(Monument, NK) |> st_buffer(dist = 20)
plot(Monument$geometry)

## Rasterize the monument polygons
Monument_Mask <- rasterize(vect(Monument), LC_NK, cover = T, background = 0)
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
Urban_mask <- classify(LC_NK, rclmat, include.lowest=TRUE)
plot(Urban_mask)




##-----------------------------------------##
#### 8. Define strategy opportunity area ####
##-----------------------------------------##

##*-- Better Strategy --*##

## Read in full BWWM data set that has each wader cluster labelled
ClustFields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Essex") |> 
              st_buffer(dist=25) |> 
              st_crop(NK)


## Rasterize the wader clusters
Better_Opp <- rasterize(vect(ClustFields), LC_NK, cover = T, background = 0)
## If a pixel is more than 50% covered by a monument then assign a pixel value of 1, using my function
Better_Opp <- Half_reclass(Better_Opp)
plot(Better_Opp)
m <- c(0.0, 0.5, NA,
       0.5, 1.5, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Better_Oppr <- classify(Better_Opp, rclmat, include.lowest=TRUE)
plot(Better_Oppr)



##*-- Bigger/More Strategy --*##

## Just take the inverse of the Better strategy clusters to get the opportunity area for bigger/more
BigMore_Opp <- Better_Opp
values(BigMore_Opp) <- Inv_scale_vals(values(BigMore_Opp))
m <- c(0.0, 0.5, NA,
       0.5, 1.5, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
BigMore_Opp <- classify(BigMore_Opp, rclmat, include.lowest=TRUE)
BigMore_Opp <- mask(BigMore_Opp , vect(NK))
plot(BigMore_Opp)



##*-- Agricultural Reversion --*##

## Reclassify raster so that it indicates if pixel is Urban/Sub-urban/Coastal or not
m <- c(0.0, 2.5, NA,
       2.5, 3.5, 1,
       3.5, 22, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Arable_Opp <- classify(LC_NK, rclmat, include.lowest=TRUE)
Arable_Opp <- mask(Arable_Opp , vect(NK))
plot(Arable_Opp)










##------------------------------------##
#### 9. Group 1: Combine Guidelines ####
##------------------------------------##

##------------------------##
## Group 1: Combine masks ##
##------------------------##

G1_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
G1_mask <- classify(G1_mask, rclmat, include.lowest=TRUE)
plot(G1_mask)
writeRaster(G1_mask, "CleanData/Guideline Creation/Kent/Kent_MasksAll_G1.tif", overwrite=TRUE)

## create a general mask to use in Scenarios later
Gen_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Gen_mask <- classify(Gen_mask, rclmat, include.lowest=TRUE)
plot(Gen_mask)
writeRaster(Gen_mask, "CleanData/Guideline Creation/Kent/Kent_MasksAll_Gen.tif", overwrite=TRUE)



##------------------------------------##
## Group 1: Combine better guidelines ##
##------------------------------------##

Bett_G1 <- InvSiteDist + AvoidRSPB + HeronriesDist + WaterAbstr
names(Bett_G1) <- "Better"

ggplot() + geom_spatraster(data=Bett_G1 |> mask(NKOutline)) + labs(fill = "Grading") + scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G1, "CleanData/Guideline Creation/Kent/Kent_Better_G1.tif", overwrite=TRUE)



##------------------------------------##
## Group 1: Combine bigger guidelines ##
##------------------------------------##

Big_G1 <- InvSiteDist + UrbanDens + AgriGrades + SaltmarshOpp + Connectr + HistoricGrass + HeronriesDist + WaterAbstr
names(Big_G1) <- "Bigger"

ggplot() + geom_spatraster(data=Big_G1 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G1, "CleanData/Guideline Creation/Kent/Kent_Bigger_G1.tif", overwrite=TRUE)



##----------------------------------##
## Group 1: Combine more guidelines ##
##----------------------------------##

More_G1 <- UrbanDens + AgriGrades + SaltmarshOpp + IsoConnectr + HistoricGrass + HeronriesDist + WaterAbstr
names(More_G1) <- "More"

ggplot() + geom_spatraster(data=More_G1 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G1, "CleanData/Guideline Creation/Kent/Kent_More_G1.tif", overwrite=TRUE)



##----------------------------------------------------##
## Group 1: Combine agri-conversion Bigger guidelines ##
##----------------------------------------------------##

ArableBig_G1 <- Big_G1
names(ArableBig_G1) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableBig_G1 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableBig_G1, "CleanData/Guideline Creation/Kent/Kent_ArableBig_G1.tif", overwrite=TRUE)


##--------------------------------------------------##
## Group 1: Combine agri-conversion More guidelines ##
##--------------------------------------------------##

ArableMore_G1 <- More_G1
names(ArableMore_G1) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableMore_G1 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableMore_G1, "CleanData/Guideline Creation/Kent/Kent_ArableMore_G1.tif", overwrite=TRUE)





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
writeRaster(G2_mask, "CleanData/Guideline Creation/Kent/Kent_MasksAll_G2.tif", overwrite=TRUE)



##------------------------------------##
## Group 2: Combine better guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Bett_G2 <- Hydro_Height
names(Bett_G2) <- "Better"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Bett_G2 |> mask(NKOutline)) + labs(fill = "Grading") + scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G2, "CleanData/Guideline Creation/Kent/Kent_Better_G2.tif", overwrite=TRUE)



##------------------------------------##
## Group 2: Combine bigger guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Big_G2 <- InvSiteDist + UrbanDens + AESWader + Hydro_Height
names(Big_G2) <- "Bigger"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Big_G2 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G2, "CleanData/Guideline Creation/Kent/Kent_Bigger_G2.tif", overwrite=TRUE)



##----------------------------------##
## Group 2: Combine more guidelines ##
##----------------------------------##

## Add together all rules and assign name to layer
More_G2 <- SiteDist + UrbanDens + AESWader + Hydro_Height
names(More_G2) <- "More"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=More_G2 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G2, "CleanData/Guideline Creation/Kent/Kent_More_G2.tif", overwrite=TRUE)


##----------------------------------------------------##
## Group 2: Combine agri-conversion Bigger guidelines ##
##----------------------------------------------------##

ArableBig_G2 <- Big_G2
names(ArableBig_G2) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableBig_G2 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableBig_G2, "CleanData/Guideline Creation/Kent/Kent_ArableBig_G2.tif", overwrite=TRUE)


##--------------------------------------------------##
## Group 2: Combine agri-conversion More guidelines ##
##--------------------------------------------------##

ArableMore_G2 <- More_G2
names(ArableMore_G2) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableMore_G2 |> mask(NKOutline)) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableMore_G2, "CleanData/Guideline Creation/Kent/Kent_ArableMore_G2.tif", overwrite=TRUE)





##-------------------------------------------------##
#### 11.0 Plot graded maps with opportunity areas ####
##-------------------------------------------------##

##---------------------------------##
## 11.1 Prepare Landscape Canvas ####
##---------------------------------##


##-- Label lowland wet grassland fields --##

## Read in the somerset canvas with all land parcel polygons
CanvasGr <- st_read("CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp")
plot(CanvasGr$geometry)

## Read in a raster icreaetd of lowland wet grassland extent in the somerset levels
LowGrass <- rast("RawData/LowlandWetGrassRasters/NKent_LowWetGrass.tif")
plot(LowGrass)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(LowGrass, vect(CanvasGr), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 25% covered by lowland wet grassland pixels
CanvasGr <- CanvasGr |>
          mutate(WetGrassProp = Scores$layer) |>
          filter(WetGrassProp > 0.50) |> select(-WetGrassProp)
plot(CanvasGr$geometry)


##-- Label arable opportunity fields --##

## Read in the Broads canvas with all land parcel polygons
CanvasAr <- st_read("CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp")
plot(CanvasAr$geometry)

## Read in a raster I created of lowland wet grassland extent in the Broads
SuitArable <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_ArableSuitable.tif")
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
CanvasAr <- CanvasAr |> st_intersection(NKOutline)
plot(CanvasAr$geometry)


##-- Extract the grades/masks for each of the Lawton strategies --##

## Code to read in the rules if needed
# G1_mask <- rast("CleanData/Guideline Creation/Kent/Kent_MasksAll_G1.tif")
# Bett_G1 <- rast("CleanData/Guideline Creation/Kent/Kent_Better_G1.tif")
# Big_G1 <- rast("CleanData/Guideline Creation/Kent/Kent_Bigger_G1.tif")
# More_G1 <- rast("CleanData/Guideline Creation/Kent/Kent_More_G1.tif")
# ArableConv_G1 <- rast("CleanData/Guideline Creation/Kent/Kent_ArableRev_G1.tif")
# 
# G2_mask <- rast("CleanData/Guideline Creation/Kent/Kent_MasksAll_G2.tif")
# Bett_G2 <- rast("CleanData/Guideline Creation/Kent/Kent_Better_G2.tif")
# Big_G2 <- rast("CleanData/Guideline Creation/Kent/Kent_Bigger_G2.tif")
# More_G2 <- rast("CleanData/Guideline Creation/Kent/Kent_More_G2.tif")
# ArableConv_G2 <- rast("CleanData/Guideline Creation/Kent/Kent_ArableRev_G2.tif")


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
ArableBig1 <- extract(ArableBig_G1, CanvasAr, fun = mean, na.rm = T)
ArableMore1 <- extract(ArableMore_G1, CanvasAr, fun = mean, na.rm = T)
ArableBig2 <- extract(ArableBig_G2, CanvasAr, fun = mean, na.rm = T)
ArableMore2 <- extract(ArableMore_G2, CanvasAr, fun = mean, na.rm = T)
CanvasAr$ArableBig_G1 <- ArableBig1$ArableConv
CanvasAr$ArableBig_G2 <- ArableBig2$ArableConv
CanvasAr$ArableMore_G1 <- ArableMore1$ArableConv
CanvasAr$ArableMore_G2 <- ArableMore2$ArableConv


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
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("North Kent: Better (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBetter)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G1_Better.png", width = 20, height = 12, units = "cm")



##-------------------------------------##
## 11.4 Group 1: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("North Kent: Bigger (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBig)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G1_Bigger.png", width = 20, height = 12, units = "cm")



##-----------------------------------##
## 11.5 Group 1: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("North Kent: More (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrMore)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G1_More.png", width = 20, height = 12, units = "cm")



##-----------------------------------------##
## 11.6 Group 1: Plot arable big grading ####
##-----------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG1 <- filter(CanvasAr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG1, mapping=aes(geometry=geometry, fill = ArableBig_G1), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("North Kent: Arable Reversion for Bigger (Conservationists)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG1)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G1_ArableBig.png", width = 20, height = 12, units = "cm")



##------------------------------------------##
## 11.7 Group 1: Plot arable more grading ####
##------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG1 <- filter(CanvasAr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG1, mapping=aes(geometry=geometry, fill = ArableMore_G1), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("North Kent: Arable Reversion for More (Conservationists)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG1)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G1_ArableMore.png", width = 20, height = 12, units = "cm")






##---------------------------------------##
## 11.8 Group 1: Plot Map of Landscape ####
##---------------------------------------##

## Read in priority landscape boundary
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
Kent <- MyBoxes |> filter(Att3Value == "North Kent") 

## Read in RSPB reserve outlines
Reserves <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")
ReservesEss <- st_crop(Reserves, (NK |> st_buffer(dist=500)))

## make plot
ggplot() + 
  ## add polygons
  geom_sf(data=Kent, mapping=aes(geometry=geometry), fill = "lightblue", colour = NA) + 
  geom_sf(data=ReservesEss, mapping=aes(geometry=geometry), fill = "red", alpha = 0.5, colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## Add North arrow and scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  ## set labels
  ggtitle("North Kent") +
  ## set them
  theme_light() + 
  GeneralThemeing +
  theme(legend.position = "none")

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_LandscapeMap.png", width = 20, height = 12, units = "cm")




##-------------------------------------##
## 11.9 Group 2: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasGrBetter <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("North Kent: Better (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBetter)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G2_Better.png", width = 20, height = 12, units = "cm")
      


##-------------------------------------##
## 11.10 Group 2: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("North Kent: Bigger (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrBig)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G2_Bigger.png", width = 20, height = 12, units = "cm")




##------------------------------------##
## 11.11 Group 2: Plot more grading ####
##------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("North Kent: More (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing
rm(CanvasGrMore)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G2_More.png", width = 20, height = 12, units = "cm")




##------------------------------------------##
## 11.12 Group 2: Plot arable big grading ####
##------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG2 <- filter(CanvasAr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG2, mapping=aes(geometry=geometry, fill = ArableBig_G2), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("North Kent: Arable Reversion for Bigger (Land Managers)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG2)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G2_ArableBig.png", width = 20, height = 12, units = "cm")




##-------------------------------------------##
## 11.13 Group 2: Plot arable more grading ####
##-------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG2 <- filter(CanvasAr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## add landscape outline
  geom_sf(data=NKOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG2, mapping=aes(geometry=geometry, fill = ArableMore_G2), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("North Kent: Arable Reversion for More (Land Managers)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG2)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/NorthKent_G2_ArableMore.png", width = 20, height = 12, units = "cm")

