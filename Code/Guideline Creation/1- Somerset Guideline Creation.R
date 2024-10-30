##----------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 09/05/2024
##
## Aim: Create all guidelines for scenario modelling of the Somerset Landscape
## 
##----------------------------------------------------------## 


## Load in packages
pacman::p_load(tidyverse, data.table, sf, terra, tidyterra, gstat, stars, automap, spatstat, ggspatial, leastcostpath)
options(scipen=999) # turn off scientific notation

## Load in various helper functions
source("Code/Helper functions.R")



##----------------------------------##
#### 0. Read in general data sets ####
##----------------------------------##

## Read in my priority landscapes boundary and buffer it by 3000m
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
SomOutline <- MyBoxes |> filter(Att3Value == "Somerset Levels and Moors")
Som <- MyBoxes |> filter(Att3Value == "Somerset Levels and Moors") |> st_buffer(dist=3000)
rm(MyBoxes) # free up space

## Read in UKCEH landcover and then crop and mask it with the landscape boundary
## This will be used as the canvas to map all other variables across
LCfull <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]]
LC_Som <- crop(LCfull, vect(Som)) |> mask(vect(Som))
rm(LCfull)
plot(LC_Som) # free up space




##-------------------------------##
#### 1. Better only guidelines ####
##-------------------------------##

##*-- Target smallest existing populations --*##

## Read in full BWWM data set that has each wader cluster labelled
Reg_Fields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Somerset")

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
write_csv(PopSize, "CleanData/Guideline Creation/Somerset/Som_ClusterPopSize.csv")





##-------------------------------##
#### 2. Bigger only guidelines ####
##-------------------------------##

##*-- Target areas within landscape recovery block --*##

LR <- st_read("RawData/Somerset Levels/LR Boundary/LR boundary zone for consultation.shp") |> 
      st_transform(crs = st_crs(Som)) |> 
      mutate(Val = 1)

## Rasterize the landscape recovery polygon
LRr <- rasterize(vect(LR), LC_Som, field = "Val", fun = "max",  background = 0)
ggplot() + geom_spatraster(data=LRr) + scale_fill_viridis_c() + theme_light()



##*-- Target areas near sites with more breeding waders for expansion --*##

## Loop through each cluster
for(j in 1:length(unique(PopSize$ClustGroup))){

   ## Creat raster of buffered cluster where in pixel inside the buffered cluster recieves the value equal to the number of breeding wader pairs  
   Subbuf <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 2000)), LC_Som, field = "ClustPop", fun = "max",  background = NA)
   
   ## Create the same raster as above but this time do it for the unbuffered cluster (this is so we can calcualte the distance to the cluster)
   Subno <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 1)), LC_Som, field = "ClustPop", fun = "max",  background = NA)
   
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

##*-- Target sites that historically had breeding waders --*##

## Read in data set of digistised wader records from 1976
## These have been digisised so that a single point contains the total number of Redshank, Lapwing and Snipe on a moor by moor basis
## I also random added on absences to the data set in the hashed area on the abundance maps provided in the report
## Report: An ornithological survey of the Somerset Levels 1976-66, P. Round. 1978, Repory by Wessex Water Authority and RSPB. 
Historic <- st_read("RawData/Somerset Levels/Historic Wader Abund/Combined_Wader_Abundance 1977.shp") |> 
            st_transform(crs = st_crs(Som))  
Historic$PointNum <- 1:length(Historic)

## Read in the polygons of the wader sites from the cluster analysis
Reg_Fields <- st_read("CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Sites.shp")

## Extract just the historic counts that fall outside of existing wader sites
Inter <- st_intersection(Historic, Reg_Fields)
Historic <- filter(Historic, !PointNum %in% c(Inter$PointNum))
plot(Historic$geometry)


## Carry out ordinary kriging to interpolate the points into a continuous raster
## Guidance on this webpage: https://geobgu.xyz/r/spatial-interpolation-of-point-data.html#ordinary-kriging

## Fit a variogram model to an empirical variogram using function from automap
## The function chooses the best fitting type of model, and also fine tunes its parameters automatically
v_mod <- autofitVariogram(Abundance ~ 1, Historic)
plot(v_mod)

## The variogram model can then be passed to the gstat function, and we can carry on with the Ordinary Kriging interpolation
gsta <- gstat(formula = Abundance ~ 1, model = v_mod$var_model, data = Historic)
Canv <- LC_Som
values(Canv) <- 1
zPred <- predict(gsta, st_as_stars(Canv))

## Subset the predicted values attribute and rename it
HistAbund <- zPred['var1.pred',,]
names(HistAbund) = 'abund'

## Convert the stars output to a terra object raster and scale the values
HistAbund <- as(HistAbund, "SpatRaster")
values(HistAbund) <- scale_vals(values(HistAbund))

## plot the scales final raster along with the points used to do the krigin
ggplot() + geom_spatraster(data=HistAbund) + scale_fill_viridis_c() + theme_light()
rm(Historic, v_mod, gsta, zPred, Canv); gc()





##-------------------------------##
#### 4. Bigger/More guidelines ####
##-------------------------------##


##*-- Target fields with lower public footfall --*##

## Read in the shapefiles of the paths from the OrVal data set
## PID links that path lines to the attribute data
Paths <- st_read("RawData/ORVal/Paths/paths_england.shp")

## Read in the attribute data for the access points and streamline to retain just visitor numbers
Simple <- read_csv("RawData/ORVal/Somerset Simple Data/ORVal Somerset Site Info Paths1.csv")
Simple <- Simple |> select(pid, vis)

## Crop the paths to just the region of interest
Paths <- st_intersection(Paths, Som) |> select(PID)
plot(Paths$geometry)

## Join paths with attributes data, a few paths are lost as they did not have attribute data
PathsAtt <- inner_join(Paths, Simple, by = c("PID"= "pid"))
ggplot() + geom_sf(data = Som) +  geom_sf(data=st_buffer(PathsAtt, dist = 500), mapping=aes(colour=vis)) + theme_light() # 


## Calculate the distance between all pixels and the nearest path
PathDist <- distance(LC_Som, vect(PathsAtt), rasterize=T)
PathDist <- mask(PathDist, vect(Som))
plot(PathDist)

## Only retain path distances within 500m of a path
## This was the distance that the below paper found that human presence could impact breeding waders
## https://onlinelibrary.wiley.com/doi/10.1111/j.1474-919X.2008.00889.x
PathsAttBuf <- st_buffer(PathsAtt, dist = 500)
PathsAttBuf <- st_combine(PathsAttBuf)
PathDist <- mask(PathDist, vect(PathsAttBuf))
plot(PathDist)

## For each path buffered by 500m calcualte the max visitor value for each pixel
PathsAttBuf <- st_buffer(PathsAtt, dist = 500)
for(j in 1:nrow(PathsAttBuf)){
  message(j)
  PathsAttBuf1 <- PathsAttBuf[j,]
  PathValsr <- rasterize(vect(PathsAttBuf1), LC_Som, field = "vis", fun = "max",  background = 0)
  if(j==1){AllPathValsr <- PathValsr}else{AllPathValsr <- max(AllPathValsr, PathValsr)}
}
plot(AllPathValsr)

## Take inverse of distance to paths so pixels closer to path have higher values
values(PathDist) <- Inv_scale_vals(values(PathDist))
plot(PathDist)

## Now multiple the inverse distance to path by the max visitor value to create an index of disturbance
DisurbanceIndex <- PathDist*AllPathValsr
plot(DisurbanceIndex)

## finally take inverse and change all NAs to 1
## Inverse insures that all low value disturbed pixels have a higher grade
## All NA pixels need to be 1 as these are not in proximity to a path
values(DisurbanceIndex) <- Inv_scale_vals(values(DisurbanceIndex))
values(DisurbanceIndex) <- ifelse(is.na(values(DisurbanceIndex))==T, 1, values(DisurbanceIndex))
plot(DisurbanceIndex)
rm(PathDist, PathsAttBuf, PathsAttBuf1, AllPathValsr, PathValsr, PathsAtt, Paths, Simple); gc()




##*-- Target areas that link existing wader sites --*##

## Read in polygons of the wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Sites.shp")
WaderSites <- WaderSites |> st_centroid()
plot(WaderSites$geometry)

## Read in raster of potential lowland wet grassland and suitable arable land
## Then combine so have tiered values
LowGrass <- rast("RawData/LowlandWetGrassRasters/Som_LowWetGrass.tif")
ArableSuit <- rast("CleanData/Scenarios/3-DefineActionAreas/Som_ArableSuitable.tif") |> resample(LowGrass)
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
LC_Som2km <- aggregate(LC_Som, fact = 80)

## Rasterize the lines connecting polygons using my 2km raster as a base
Connectr <- rasterize(vect(LCP), LC_Som2km, field = "Val", fun = "sum",  background = 0)
ggplot() + geom_spatraster(data = Connectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)

## Scale the values from the rasterized lines
values(Connectr) <- scale_vals(values(Connectr))
plot(Connectr)

## Finally resample the raster back down to 25m
Connectr <- resample(Connectr, LC_Som)

## Finally turn any NAs back to zeros
values(Connectr) <- ifelse(is.na(values(Connectr))==T, 0, values(Connectr))
plot(Connectr)

## free up space
rm(WaderSites, LowGrass, Cond, LCP, LC_Som2km); gc()

## OLD VERSION OF SITE LINKAGES BELOW
# ## Read in polygons of the wader sites
# WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Sites.shp")
# plot(WaderSites$geometry)
# 
# ## From a section above get the population size for each of the wader clusters
# PopSize2 <- PopSize |> 
#   group_by(ClustGroup) |> 
#   summarise(ClustPop = max(ClustPop, nam.rm = T)) |> 
#   as.data.frame()
# 
# ## Now add the cluster populaiton sizes onto the polygons
# WaderSites <- left_join(WaderSites, PopSize2, by = "ClustGroup")
# 
# 
# ## Calculate the shortest distance between each of the wader site polygons
# # Initialize an empty list to store the results
# results <- list()
# 
# # Loop through all pairs of polygons
# for (i in 1:(nrow(WaderSites) - 1)) {
#   for (j in (i+1):nrow(WaderSites)) {
#     # Check if the pair (i, j) has already been processed
#     pair <- paste(sort(c(WaderSites$ClustGroup[i], WaderSites$ClustGroup[j])), collapse = "_")
#     if (!(pair %in% names(results))) {
#       # Multiply the pair of values and store the result
#       result <- st_nearest_points(WaderSites[i,], WaderSites[j,])
#       Pop_size <- WaderSites$ClustPop[i] + WaderSites$ClustPop[j]
#       
#       results[[pair]] <- list(Line = result, Pop = Pop_size)
# 
#     }
#   }
# }
# 
# ## Bind the results
# library(plyr)
# results2 <- ldply (results, data.frame)
# results2 <-  st_as_sf(results2)
# ggplot(data = results2) + geom_sf()
# 
# ## No turn the lines that connect polygons into a raster, I am going to weight the lines by the size of the two populations they connect
# ## Create a base raster with a pixel size of 1km using my base land cover raster
# LC_Som500 <- aggregate(LC_Som, fact = 40)
# 
# ## Rasterize the lines connecting polygons using my 1km raster as a base
# Connectr <- rasterize(vect(results2), LC_Som500, field = "Pop", fun = "sum",  background = 0)
# ggplot() + geom_spatraster(data = Connectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)
# 
# ## Scale the values from the rasterized lines
# values(Connectr) <- scale_vals(values(Connectr))
# plot(Connectr)
# 
# ## Finally resample the raster back down to 25m and mask any pixels that are part of wader site alread
# Connectr <- resample(Connectr, LC_Som)
# Connectr <- mask(Connectr, vect(WaderSites), inverse=T)
# ## Finally turn any NAs back to zeros
# values(Connectr) <- ifelse(is.na(values(Connectr))==T, 0, values(Connectr))
# plot(Connectr)
# rm(WaderSites, PopSize2, pair, result, results, results2, LC_Som500);gc()




##*-- Target hydrological units with more gradually sloping boundaries --*##

# ## 2m LiDAR raster for the Somerset
# Som_Elev <- rast("RawData/LiDAR/Somerset_DTM_2m.tif")
# 
# ## Shapefile of Hydro Units 
# Som_Units <- st_read("RawData/HydroUnits/Somerset/Main_Management_Areas.shp")
# ## Add row values to Somerset hydro units data set
# Som_Units$id <- 1:nrow(Som_Units)
# plot(Som_Units$geometry)
# ## Add columns to add different surface area value for the unit
# Som_Units$Area2D <- NA
# Som_Units$Area3D <- NA
# 
# 
# ## Loop through each hydro unit and extract the 2D and 3D surface areas
# for(j in 1:nrow(Som_Units)){
#   
#   message("iteration ", j)
#   
#   ## filter out the hydro unit outline
#   Som_Unit <- vect(Som_Units[j,])
#   
#   ## Crop the hydro unit from the LiDAR data
#   Unit_Lidar <- crop(Som_Elev, Som_Unit, mask = T)
#   plot(Unit_Lidar)
#   
#   ## Calculate the cell sizes of all pixels in the hydro unit
#   area <- cellSize(Unit_Lidar) 
#   area <- mask(area, Unit_Lidar) # mask to the area of the hydrounit
#   plot(area)
#   
#   ## Calcualte the slop of eaxh pixel in the hydro unit
#   slope <- terrain(Unit_Lidar, v="slope", neighbors=8, unit = "radians")
#   plot(slope)
#   
#   ## Simple calculation to work out the 3D area of the pixel
#   adjarea <- area / cos(slope)
#   area <- mask(area, adjarea) # crop 2D pixel area to that of 3D pixel area so calculating area over same number of pixels
#   plot(adjarea)
#   
#   ## Assign 2D and 3D surface areas to each hydro unit
#   Som_Units$Area2D[j] <- sum(values(area), na.rm = T)
#   Som_Units$Area3D[j] <- sum(values(adjarea), na.rm = T)
#   
# }
# 
# ## plot the data set and save it under a new object
# ggplot() + geom_sf(data = Som_Units, mapping = aes(geometry=geometry, fill = (Area3D/Area2D)))
# Som_UnitsBowl <- Som_Units
# rm(Som_Elev, Som_Units, Unit_Lidar, area, slope, adjarea); gc() # free up space



##*-- Target areas within existing Moorland Associations --*##

## Read in the two shapefiles of the farming clusters
FarmClust1 <- st_read("RawData/Somerset Levels/Farming Clusters/Current MA Bounds/All.shp") |> select(geometry)
FarmClust2 <- st_read("RawData/Somerset Levels/Farming Clusters/LAPWDP MA/LAPWDP_All.shp") |> select(geometry)
FarmClust1 <- st_transform(FarmClust1, crs = st_crs(FarmClust2))
FarmClust1 <- rbind(FarmClust1, FarmClust2); rm(FarmClust2)

## Rasterize the landscape recovery polygon
FarmClusts <- rasterize(vect(FarmClust1), LC_Som, cover = T, background = 0)
FarmClusts <- Half_reclass(FarmClusts)
ggplot() + geom_spatraster(data=FarmClusts) + scale_fill_viridis_c() + theme_light()
rm(FarmClust1); gc()



##*-- Target areas further away from urban areas --*##

## Reclassify raster so that it indicates if pixel is Urban/Sub-urban or not
m <- c(0.0, 20.5, NA,
       20.5, 22.0, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Urban <- classify(LC_Som, rclmat, include.lowest=TRUE)
plot(Urban)

## Calculate distance from NA pixels (other habitats) to none NA pixel (Urban areas)
DistToUrban <- distance(Urban)
DistToUrban <- terra::mask(DistToUrban, LC_Som) # crop to ROI

## Re scale values so largest distance to urban = 1
values(DistToUrban) <- scale_vals(values(DistToUrban))
plot(DistToUrban)
rm(Urban); gc()



##*-- Target SSSI area --*##

## Read in shapefile of SSSI areas
SSSI <- st_read("RawData/SSSI/data/Sites_of_Special_Scientific_Interest_England.shp")
SSSI <- st_crop(SSSI, ext(LC_Som)) # crop to save space

## Rasterize the landscape recovery polygon
SSSI <- rasterize(vect(SSSI), LC_Som, cover = T, background = 0)
SSSI <- Half_reclass(SSSI)
ggplot() + geom_spatraster(data=SSSI) + scale_fill_viridis_c() + theme_light()



##*-- Target hydro logical units that have few landowners --*##

## Read in the RPA point parcels for the Somerset regions
PPoints <- st_read("RawData/RPA/rpa_parcel_points_Somerset/rpa_parcel_pointsPoint.shp")

## Read in the RPA point parcels for the Somerset regions
PointsID <- fread("RawData/RPA/RPA_ParcelList_20240424&Suffolk.csv")

## Join together the field parcel points and the cutomer data from the RPA
Join <- inner_join(PointsID, PPoints, by = join_by(NGC==parcel_ref)) |> st_as_sf() |> select(FARM_ID)
plot(Join$geometry)

## Read in the shapevfile of Hydro Units 
Som_Units <- st_read("RawData/HydroUnits/Somerset/Main_Management_Areas.shp")
Som_Units$id <- NULL

## Extract to label the parcel points with he hydro unit and then assign these back to the parcel points data set
TestExtr <- extract(vect(Som_Units), vect(Join))
Join$Site <- TestExtr$Site

## Now summarise the number of different farm customers for each of the hydro units
UnitLandowners <- Join |> 
                  filter(is.na(Site) == F) |> 
                  group_by(Site) |> 
                  dplyr::summarise(No_farms = length(unique(FARM_ID))) |> st_drop_geometry()

## Calculate the density of landowners within each hydro unit
Som_Units <- left_join(Som_Units, UnitLandowners, by = "Site") |> 
             mutate(Area = st_area(geometry),
                    LandownerDens = No_farms/as.numeric(Area/1000000)) 

## Rasterize the hydro unit data using the land owner density values
LandownerDens <- rasterize(vect(Som_Units), LC_Som, field = "LandownerDens", fun = "max",  background = NA)

## Inverse scale the values so the units with the lowest landowner denisty have a score closer to 1
values(LandownerDens) <- Inv_scale_vals(values(LandownerDens))
values(LandownerDens) <- ifelse(is.na(values(LandownerDens))==T, 0, values(LandownerDens))
plot(LandownerDens)  
rm(PPoints, PointsID, Join, TestExtr, UnitLandowners, Som_Units); gc()




##-------------------------------------##
#### 5. Arable conversion guidelines ####
##-------------------------------------##

##*-- Target arable conversion on land where Maize is grown --*##

Crops <- st_read("RawData/LandCover/Land Cover plus Crops 2021/LC_Crops21.shp")
Crops <- filter(Crops, crop_code == "ma")

## Crop out the crop polygons for the buffered priority landscape
Crops <- st_intersection(Crops, Som)

## Rasterize the maize polygons
CropsMaize <- rasterize(vect(Crops), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by silty soils then assign a pixel value of 1, using my function
CropsMaize <- Half_reclass(CropsMaize)
plot(CropsMaize)
rm(Crops); gc()



##*-- Target arable reversion near to raised water level areas --*##

## Read in the shapefile of the rasied water levels
RWLA <- st_read("RawData/Somerset Levels/RWLAshapefiles/RWLAS.shp")
RWLA <- st_transform(RWLA, crs = st_crs(LC_Som))
plot(RWLA$geometry)

## Rasterize the rasied water levels polygons
RWLA <- rasterize(vect(RWLA), LC_Som, cover = T, background = 0)
RWLA <- Half_reclass(RWLA)
plot(RWLA)

## Reclassify raised water level raster so 0 pixels become NA
m <- c(0, 0.5, NA,
       0.5, 1, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
RWLA2 <- classify(RWLA, rclmat, include.lowest=TRUE)
plot(RWLA2)

## Calculate distance from NA pixels (other) to none NA pixel (raised water levels)
DistToRWLA <- distance(RWLA2, unit="m")
m <- c(0, 1000, 1,
       1000, 200000, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
DistToRWLA2 <- classify(DistToRWLA, rclmat, include.lowest=TRUE)
DistToRWLA <- terra::mask(DistToRWLA, DistToRWLA2) # crop to ROI
values(DistToRWLA) <- Inv_scale_vals(values(DistToRWLA))
values(DistToRWLA) <- ifelse(is.na(values(DistToRWLA))==T, 0, values(DistToRWLA))
plot(DistToRWLA)
rm(RWLA, RWLA2, DistToRWLA2); gc()



##*-- Target arable reversion in isolated patches within grassland --*##

## Read in raster of potential lowland wet grassland 
LowGrass <- rast("RawData/LowlandWetGrassRasters/Som_LowWetGrass.tif")
plot(LowGrass)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
LowGrassDens <- focal(LowGrass, w=41, fun="mean", na.rm = T)
LowGrassDens <- mask(LowGrassDens, LowGrass) # crop to the extent of the original raster
plot(LowGrassDens) # plot smoothed raster
#rm(LowGrass); gc()



##*-- Target arable conversion on areas adjacent to wet grassland --*##

## Reclassify raised water level raster so 0 pixels become NA
m <- c(0, 0.5, NA,
       0.5, 1, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LowGrass <- classify(LowGrass, rclmat, include.lowest=TRUE)
plot(LowGrass)

## Calculate distance from NA pixels (other) to none NA pixel (lowland wet grassland)
LowGrass <- distance(LowGrass, unit="m")
## Create a raster where only pixels with value of <1000m have a value of 1
m <- c(0, 1000, 1,
       1000, 200000, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LowGrass2 <- classify(LowGrass, rclmat, include.lowest=TRUE)
plot(LowGrass2)
## Mask pixles only with 1km of lowland wet grassland
LowGrassDist <- terra::mask(LowGrass, LowGrass2) # crop to ROI
values(LowGrassDist) <- Inv_scale_vals(values(LowGrassDist)) # inverse scale so low distances have high value
values(LowGrassDist) <- ifelse(is.na(values(LowGrassDist))==T, 0, values(LowGrassDist))
plot(LowGrassDist)
rm(LowGrass2, LowGrass); gc()



##*-- Target arable conversion on peat soils --*##

## Read in NATMAP vector data
Soil <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## Crop the NATMAP soil vector
Soil <- st_crop(Soil, Som) # crop to study area bbox
table(Soil$SOILSCAPE) # get table of soil types

## filter out just the soils that are in some way peaty
Soil <- Soil |> 
  filter(SOILSCAPE %in% c("Blanket bog peat soils", "Fen peat soils", "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                          "Raised bog peat soils", "Slowly permeable wet very acid upland soils with a peaty surface",
                          "Very acid loamy upland soils with a wet peaty surface"))
ggplot(data =Soil) + geom_sf(fill= "blue") + theme_bw() + ggtitle("Peaty Soils")

## Rasterize this wet/peaty soils layer, so if a wet/peaty are overlaps a pixel it gets a value of 1
PeatSoils <- rasterize(vect(Soil), LC_Som, cover = T, background = 0)
PeatSoils <- Half_reclass(PeatSoils)
plot(PeatSoils, main = "Wet Soils Raster")
rm(Soil); gc() # free up space








##--------------------------------##
#### 6. All strategy guidelines ####
##--------------------------------##


##*-- Target areas with more water in future --*##

## Read in data set on water abstraction availability
WaterAbst <- st_read("RawData/Somerset Levels/Water Abstraction/All_Somerset_Catchments.shp")
ggplot() + geom_sf(data = WaterAbst, mapping=aes(geometry=geometry, fill = RsrcAvl)) # quick plot

## Rasterize the polygons data set and give the pixels the values of the % of time water is available for abstraction
WaterAbstr <- rasterize(vect(WaterAbst), LC_Som, field = "RsrcAvl", fun = "mean",  background = 0)
## Scale values and change any NAs to zeros
values(WaterAbstr) <- scale_vals(values(WaterAbstr))
values(WaterAbstr) <- ifelse(values(WaterAbstr) == "NaN", 0, values(WaterAbstr))
plot(WaterAbstr)
rm(WaterAbst); gc()




##*-- Target silty areas for Lapwing management --*##

## Read in the NATMAP soil vector
Soil <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## Create vector of soil types to retains that are silty/wet
table(Soil$SOILSCAPE) # get table of soil types
WetTypes <- c("Lime-rich loamy and clayey soils with impeded drainage",
              "Loamy soils with naturally high groundwater",
              "Loamy and clayey soils of coastal flats with naturally high groundwater",
              "Loamy and clayey floodplain soils with naturally high groundwater",
              "Naturally wet very acid sandy and loamy soils",
              "Slightly acid loamy and clayey soils with impeded drainage",
              "Slowly permeable seasonally wet slightly acid but base-rich loamy and clayey soils",
              "Slowly permeable seasonally wet acid loamy and clayey soils")

## filter data set, and assign a Val column for rasterisation
SiltSoils <- filter(Soil, SOILSCAPE %in% WetTypes) #|> mutate(Val = 1)

## Rasterize the silty soils polygons
SiltSoil <- rasterize(vect(SiltSoils), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by silty soils then assign a pixel value of 1, using my function
SiltSoil <- Half_reclass(SiltSoil)
plot(SiltSoil)
rm(Soil, WetTypes); gc()



##*-- Target areas not part of floodplain/flood storage area on Curry moor --*##

## Read in digitized hydro logical units polygons
HydroUnits <- st_read("RawData/HydroUnits/Somerset/Main_Management_Areas.shp")
## Extract the polygon for Curry Moor
FloodCurry <- filter(HydroUnits, Site == "Currymoor") 

## Rasterize the Curry moor polygons
FloodCurry <- rasterize(vect(FloodCurry), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by Curry moor then assign a pixel value of 1, using my function
FloodCurry <- Half_reclass(FloodCurry)
values(FloodCurry) <- Inv_scale_vals(values(FloodCurry))
plot(FloodCurry)
rm(HydroUnits); gc()



##*-- Target areas away from other priority habitat expansion areas --*##

## Read in Network maps for three key priority habitat inthe Somerset Levels
## Reedbed, Lowland raised bog and lowland fen
Reed <- st_read("RawData/Somerset Levels/Habitat_Networks_Individual_England_Reedbeds/Habitat_Networks_Individual_England_ReedbedsPolygon.shp") |> 
        filter(class %in% c("Fragmentation Action Zone", "Habitat Restoration-Creation", "Network Enhancement Zone 1"))
Fen <- st_read("RawData/Somerset Levels/Habitat_Networks_Individual_England_Lowland_Fen/Habitat_Networks_Individual_England_Lowland_FenPolygon.shp") |> 
        filter(class %in% c("Fragmentation Action Zone", "Habitat Restoration-Creation", "Network Enhancement Zone 1"))
Bog <- st_read("RawData/Somerset Levels/Habitat_Networks_Individual_England_Lowland_Raised_Bog/Habitat_Networks_Individual_England_Lowland_Raised_BogPolygon.shp") |> 
       filter(class %in% c("Fragmentation Action Zone", "Habitat Restoration-Creation", "Network Enhancement Zone 1"))

## Filter out these classes as the top priority for restoration
Hab1 <- bind_rows(Reed, Fen, Bog) |> filter(class %in% c("Fragmentation Action Zone", "Habitat Restoration-Creation"))
## This class is a secondary priority so will be used but given a half score
Hab0_5 <- bind_rows(Reed, Fen, Bog) |> filter(class == "Network Enhancement Zone 1")

## Rasterize the Hab1 polygons
Hab1 <- rasterize(vect(Hab1), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by Hab1 then assign a pixel value of 1, using my function
Hab1 <- Half_reclass(Hab1)
plot(Hab1)

## Rasterize the Hab0_5 polygons
Hab0_5 <- rasterize(vect(Hab0_5), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by Hab0_5 then assign a pixel value of 1, using my function
Hab0_5 <- Half_reclass(Hab0_5)
plot(Hab0_5)

## Add together the two polygons, make sure the max value is one and then inverse scale the variable
## This will ensure all areas where priority habitat expansion is not highlighted will have a value of 1
Priority_HabExp <- Hab1 + Hab0_5/2
values(Priority_HabExp) <- ifelse(values(Priority_HabExp) > 1, 1, values(Priority_HabExp))
values(Priority_HabExp) <- Inv_scale_vals(values(Priority_HabExp))
plot(Priority_HabExp)
rm(Reed, Fen, Bog, Hab1, Hab0_5); gc()



##*-- Target mid height land within hydro units --*##

## 2m LiDAR raster for the Somerset
Som_Elev <- rast("RawData/LiDAR/Somerset_DTM_2m.tif")

## Shapefile of Hydro Units 
Som_Units <- st_read("RawData/HydroUnits/Somerset/Main_Management_Areas.shp")
 
## Add values to Somerset hydro units data set
Som_Units$Region <- "Somerset"
Som_Units$id <- 1:nrow(Som_Units)
plot(Som_Units$geometry)

## create an empty list to put all the raster outputs into
List <- vector(mode='list', length= nrow(Som_Units))


## Loop through each hydro unit and extract the LiDAR data and grade it
for(j in 1:nrow(Som_Units)){
  
  message("iteration ", j)
  
  ## filter out the hydro unit outline
  Som_Unit <- vect(Som_Units[j,])
  
  ## Crop the hydro unit from the LiDAR data
  Unit_Lidar <- crop(Som_Elev, Som_Unit, mask = T)
  
  ## Grade each pixel in the hydrological unit, first extract pixels
  x <- values(Unit_Lidar)
  ## empirical cumulative distribution function
  y <- ecdf(x)(x)
  ## Now scale the elevation values so that values nearer the 50th quantile are given a value of 1 and the highest and lowest areas a value of 1
  values(Unit_Lidar) <- Inv_scale_vals(abs(y-0.5))
  
  ## add raster to slot in list
  List[[j]] <- Unit_Lidar
  
}


## free up some space first
rm(Som_Elev, Unit_Lidar); gc()

## Combine all the raster into a raster collection
List <- sprc(List)

## now merge that raster collection
Hydro_Height <- terra::merge(List)

## Change all NAs in the raster to zeros
values(Hydro_Height) <- ifelse(is.na(values(Hydro_Height))==T, 0, values(Hydro_Height))
gc()

## resample the fine-scale topography raster to my base raster
Hydro_Height <- resample(Hydro_Height, LC_Som)
values(Hydro_Height) <- ifelse(is.na(values(Hydro_Height))==T, 0, values(Hydro_Height))
plot(Hydro_Height, col=grDevices::hcl.colors(50, palette = "Sunset"))
rm(Som_Units, Som_Unit, List); gc() # free up some space




##*-- Target lowest land within hydro units --*##

## 2m LiDAR raster for the Somerset
Som_Elev <- rast("RawData/LiDAR/Somerset_DTM_2m.tif")

## Shapefile of Hydro Units 
Som_Units <- st_read("RawData/HydroUnits/Somerset/Main_Management_Areas.shp")
 
## Add values to Somerset hydro units data set
Som_Units$Region <- "Somerset"
Som_Units$id <- 1:nrow(Som_Units)
plot(Som_Units$geometry)

## create an empty list to put all the raster outputs into
List <- vector(mode='list', length= nrow(Som_Units))


## Loop through each hydro unit and extract the LiDAR data and grade it
for(j in 1:nrow(Som_Units)){
  
  message("iteration ", j)
  
  ## filter out the hydro unit outline
  Som_Unit <- vect(Som_Units[j,])
  
  ## Crop the hydro unit from the LiDAR data
  Unit_Lidar <- crop(Som_Elev, Som_Unit, mask = T)
  
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
rm(Som_Elev, Unit_Lidar); gc()

## Combine all the raster into a raster collection
List <- sprc(List)

## now merge that raster collection
Hydro_HeightLow <- terra::merge(List)

## Change all NAs in the raster to zeros
values(Hydro_HeightLow) <- ifelse(is.na(values(Hydro_HeightLow))==T, 0, values(Hydro_HeightLow))
gc()

## re-sample the fine-scale topography raster to my base raster
Hydro_HeightLow <- resample(Hydro_HeightLow, LC_Som)
values(Hydro_HeightLow) <- ifelse(is.na(values(Hydro_HeightLow))==T, 0, values(Hydro_HeightLow))
plot(Hydro_HeightLow, col=grDevices::hcl.colors(50, palette = "Sunset"))
rm(Som_Units, Som_Unit, List); gc() # free up some space




##*-- Target raised water level areas --*##

## Read in the shapefile of farming clusters
RWLA <- st_read("RawData/Somerset Levels/RWLAshapefiles/RWLAS.shp")
RWLA <- st_transform(RWLA, crs = st_crs(LC_Som))
plot(RWLA$geometry)

## Rasterize the rasied water levels polygons
RWLA <- rasterize(vect(RWLA), LC_Som, cover = T, background = 0)
RWLA <- Half_reclass(RWLA)
ggplot() + geom_spatraster(data=RWLA) + scale_fill_viridis_c() + theme_light()



##*-- Target areas with CS agreements already --*##

## Read in the Different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")
ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESS <- rbind(ESSWest, ESSEast)
rm(ESSWest, ESSEast)

## Read in the Canvas as this has the field parcel shapes 
Canv <- st_read("CleanData/Scenarios/1-Starting Canvas/Som_Canvas.shp") |> select(geometry)

## Filter out wader focused options
ESS_OpsWa <- filter(ESS, optcode %in% c("HK9", "HK10", "HK11", "HK12", "HK13", "HK14", "HK19"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
ESS_InterWa <- as.data.frame(st_covers(Canv, ESS_OpsWa))
Canv$ESS_Wader <- "N"
Canv$ESS_Wader[ESS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()

## Add column to CSS data to indicate if agreements includes options that benefit waders or general options that could benefit waders
## See Appenidx 1 of Rob's BWWM report for what is classified as "wader focused" and "general" options
## Filter out wader focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))

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
AESFields <- rasterize(vect(Canv), LC_Som, cover = T, background = 0)
AESFields <- Half_reclass(AESFields)
plot(AESFields, main = "AES fields")

## free up space
rm(Canv, ESS, ESS_InterWa, ESS_OpsWa, CSS_InterWa, CSS_Ops, CSS_OpsWa); gc()




##*-- Target periphery of peat soils edges --*##

## Read in NATMAP vector data
Soil <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")

## Crop the NATMAP soil vector
Soil <- st_crop(Soil, Som) # crop to study area bbox
table(Soil$SOILSCAPE) # get table of soil types

## filter out just the soils that are in some way peaty
SoilNP <- Soil |> 
  filter(!SOILSCAPE %in% c("Blanket bog peat soils", "Fen peat soils", "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                          "Raised bog peat soils", "Slowly permeable wet very acid upland soils with a peaty surface",
                          "Very acid loamy upland soils with a wet peaty surface", "water"))
ggplot(data =SoilNP) + geom_sf(fill= "blue") + theme_bw() + ggtitle("None-Peaty Soils")

SoilP <- Soil |> 
  filter(SOILSCAPE %in% c("Blanket bog peat soils", "Fen peat soils", "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                          "Raised bog peat soils", "Slowly permeable wet very acid upland soils with a peaty surface",
                          "Very acid loamy upland soils with a wet peaty surface"))
ggplot(data =SoilP) + geom_sf(fill= "blue") + theme_bw() + ggtitle("Peaty Soils")
rm(Soil); gc()

## Rasterize none-peaty soils
NonePeatSoils <- rasterize(vect(st_buffer(SoilNP, dist = 10)), LC_Som, cover = T, background = 0)
NonePeatSoils <- Half_reclass(NonePeatSoils)
plot(NonePeatSoils, main = "None Peaty Soils")
rm(SoilNP); gc()

## Reclassify the raster so that none peaty soils are 1 and peat soils are NA
m <- c(0, 0.5, NA,
       0.5, 1, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
NonePeatSoils <- classify(NonePeatSoils, rclmat, include.lowest=TRUE) # reclassify the raster
plot(NonePeatSoils, main = "None peaty soils") # check reclassification has worked

## Calculate distance from NA pixels (peat soil) to none NA pixel (none peat soils)
## then crop to just priority landscape
PeatEdge <- distance(NonePeatSoils)
PeatEdge <- mask(PeatEdge, vect(SoilP)) |> mask(SomOutline)
plot(PeatEdge, main = "Distance to peat soil edge")
rm(NonePeatSoils, SoilP); gc()

## Now any pixels that are peaty soils take the inverse scaled distance to none peaty soils
values(PeatEdge) <- ifelse(values(PeatEdge)==0, NA, values(PeatEdge))
values(PeatEdge) <- Inv_scale_vals(values(PeatEdge))
values(PeatEdge) <- ifelse(is.na(values(PeatEdge))==T, 0, values(PeatEdge))
plot(PeatEdge, main = "Inverse Scaled Distance to peat soil edge")









##----------------------##
#### 7. Mask Creation ####
##----------------------##


##*--Avoid priority habitats --*##

## Read in the NE Priority Habitat
NE <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")

## Crop out the Natural England priority habitats for the buffered priority landscape
NE <- st_intersection(NE, Som)

## Get a list of the main habitats and additional habitats
table(NE$mainhabs)
table(NE$addhabs)

## Create vector of habitats that I want to mask out
MAINHABS <- c("Deciduous woodland", "Deciduous woodland,Lowland raised bog",
              "Lowland calcareous grassland", 
              "Lowland dry acid grassland", "Lowland dry acid grassland,Lowland heathland",
              #"Lowland fens", "Lowland fens,Reedbeds",
              "Lowland heathland",
              "Lowland raised bog",
              "Reedbeds",
              "Traditional orchard")

## Create column of just the first 5 characters of the additional habitats column
NE <- NE |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))

## Filter the data set by Main habitat and then if no main habitat retain those that did not have CFPGM in the additional habitat column
NE_Area_Masks <- filter(NE, mainhabs %in% MAINHABS |
                           (mainhabs == "No main habitat but additional habitats present" & !addhabs == "CFPGM" ))

## plot these masks
ggplot() + geom_sf(data = NE_Area_Masks, mapping = aes(geometry = geometry, fill = mainhabs), colour = NA) + theme_light()

# ## Put a 50m buffer around all these priority habitat
# NE_Area_Masks <- st_buffer(NE_Area_Masks, dist = 50)

## Rasterize the priotity habitat polygons
NE_Mask <- rasterize(vect(NE_Area_Masks), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1, using my function
NE_Mask <- Half_reclass(NE_Mask)
plot(NE_Mask)
rm(NE, MAINHABS, NE_Area_Masks); gc() # free up space



##*--Avoid peat extraction sites --*##

## Read in peat extraction sites that I digitized from aerial images
PeatSites <- st_read("RawData/Somerset Levels/Peat Extraction Sites/Somerset Peat Extraction Sites.shp")

## Rasterize the peat extraction sites polygons
PeatExtr_Mask <- rasterize(vect(PeatSites), LC_Som, cover = T, background = 0)
## If a pixel is more than 50% covered by a peat extraction site then assign a pixel value of 1, using my function
PeatExtr_Mask <- Half_reclass(PeatExtr_Mask)
plot(PeatExtr_Mask)
rm(PeatSites); gc() # free up space



##*--Avoid scheduled monuments --*##

## Read in polygons of scheduled monuments from National Heritage
Monument <- st_read("RawData/Monuments/National_Heritage_List_for_England/Scheduled_Monuments.shp")

## Crop out the monuments for the buffered priority landscape
Monument <- st_intersection(Monument, Som) |> st_buffer(dist = 20)
plot(Monument$geometry)

## Rasterize the monument polygons
Monument_Mask <- rasterize(vect(Monument), LC_Som, cover = T, background = 0)
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
Urban_mask <- classify(LC_Som, rclmat, include.lowest=TRUE)
plot(Urban_mask)




##-----------------------------------------##
#### 8. Define strategy opportunity area ####
##-----------------------------------------##

##*-- Better Strategy --*##

## Read in full BWWM data set that has each wader cluster labelled
ClustFields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Somerset") |> 
              st_buffer(dist=25)

## Rasterize the wader clusters
Better_Opp <- rasterize(vect(ClustFields), LC_Som, cover = T, background = 0)
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
BigMore_Opp <- mask(BigMore_Opp , vect(Som))
plot(BigMore_Opp)



##*-- Agricultural Reversion --*##

## Reclassify raster so that it indicates if pixel is Urban/Sub-urban/Coastal or not
m <- c(0.0, 2.5, NA,
       2.5, 3.5, 1,
       3.5, 22, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Arable_Opp <- classify(LC_Som, rclmat, include.lowest=TRUE)
Arable_Opp <- mask(Arable_Opp , vect(Som))
plot(Arable_Opp)







##------------------------------------##
#### 9. Group 1: Combine Guidelines ####
##------------------------------------##

##------------------------##
## Group 1: Combine masks ##
##------------------------##

G1_mask <- max(NE_Mask, PeatExtr_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
G1_mask <- classify(G1_mask, rclmat, include.lowest=TRUE)
plot(G1_mask)
writeRaster(G1_mask, "CleanData/Guideline Creation/Somerset/Som_MasksAll_G1.tif", overwrite=TRUE)

## create a general mask to use in Scenarios later
Gen_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Gen_mask <- classify(Gen_mask, rclmat, include.lowest=TRUE)
plot(Gen_mask)
writeRaster(Gen_mask, "CleanData/Guideline Creation/Somerset/Som_MasksAll_Gen.tif", overwrite=TRUE)


##------------------------------------##
## Group 1: Combine better guidelines ##
##------------------------------------##

Bett_G1 <- Hydro_Height + WaterAbstr + Priority_HabExp
names(Bett_G1) <- "Better"
#Bett_G1 <- Bett_G1*G1_mask*Better_Oppr
ggplot() + geom_spatraster(data=Bett_G1) + labs(fill = "Grading") +
  scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G1, "CleanData/Guideline Creation/Somerset/Som_Better_G1.tif", overwrite=TRUE)


##------------------------------------##
## Group 1: Combine bigger guidelines ##
##------------------------------------##

Big_G1 <- Hydro_Height + WaterAbstr + LRr + (SiltSoil/2) + FloodCurry + Priority_HabExp + DisurbanceIndex + Connectr
names(Big_G1) <- "Bigger"
#Big_G1 <- Big_G1*G1_mask*BigMore_Opp
ggplot() + geom_spatraster(data=Big_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G1, "CleanData/Guideline Creation/Somerset/Som_Bigger_G1.tif", overwrite=TRUE)


##----------------------------------##
## Group 1: Combine more guidelines ##
##----------------------------------##

More_G1 <- Hydro_Height + WaterAbstr + HistAbund + (SiltSoil/2) + FloodCurry + Priority_HabExp + DisurbanceIndex + Connectr
names(More_G1) <- "More"
#More_G1 <- More_G1*G1_mask*BigMore_Opp
ggplot() + geom_spatraster(data=More_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G1, "CleanData/Guideline Creation/Somerset/Som_More_G1.tif", overwrite=TRUE)

## Option to carry out sclaing at this stage of the analysis
# More_G1 <- c(Hydro_Height, HistAbund ,SiltSoil ,FloodCurry ,Priority_HabExp)
# More_G1 <- More_G1*G1_mask*BigMore_Opp
# plot(More_G1)
# 
# ## Now scale all of the raster in the stack
# nx <- minmax(More_G1)    
# rn <- (More_G1 - nx[1,]) / (nx[2,] - nx[1,])
# plot(rn)
# ggplot() + geom_spatraster(data=More_G1[[3]]) + scale_fill_viridis_c() + theme_light()
# rn2 <- sum(rn)
# ggplot() + geom_spatraster(data=rn2) + scale_fill_viridis_c() + theme_light()


##---------------------------------------------##
## Group 1: Combine agri-conversion guidelines ##
##---------------------------------------------##

ArableConv_G1 <- CropsMaize 
names(ArableConv_G1) <- "ArableConv"
#ArableConv_G1 <- ArableConv_G1*G1_mask*Arable_Opp
ggplot() + geom_spatraster(data=ArableConv_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableConv_G1, "CleanData/Guideline Creation/Somerset/Som_ArableRev_G1.tif", overwrite=TRUE)





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
writeRaster(G2_mask, "CleanData/Guideline Creation/Somerset/Som_MasksAll_G2.tif", overwrite=TRUE)


##------------------------------------##
## Group 2: Combine better guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Bett_G2 <- Hydro_Height + WaterAbstr + RWLA
names(Bett_G2) <- "Better"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Bett_G2) + labs(fill = "Grading") +
  scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G2, "CleanData/Guideline Creation/Somerset/Som_Better_G2.tif", overwrite=TRUE)


##------------------------------------##
## Group 2: Combine bigger guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Big_G2 <- Hydro_Height + FloodCurry + HistAbund + WaterAbstr + FarmClusts + DistToUrban + RWLA + SSSI + LandownerDens + ExpandClust
names(Big_G2) <- "Bigger"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Big_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G2, "CleanData/Guideline Creation/Somerset/Som_Bigger_G2.tif", overwrite=TRUE)


##----------------------------------##
## Group 2: Combine more guidelines ##
##----------------------------------##

## Add together all rules and assign name to layer
More_G2 <- Hydro_Height + Connectr + FloodCurry + HistAbund + WaterAbstr + FarmClusts + DistToUrban + RWLA + SSSI + LandownerDens
names(More_G2) <- "More"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=More_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G2, "CleanData/Guideline Creation/Somerset/Som_More_G2.tif", overwrite=TRUE)


##----------------------------------------------##
## Group 2: Combine agri-conversion guidelines ##
##---------------------------------------------##

## Add together all rules and assign name to layer
ArableConv_G2 <- CropsMaize + (LowGrassDens/2) + (LowGrassDist/2) + DistToRWLA 
names(ArableConv_G2) <- "ArableConv"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=ArableConv_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableConv_G2, "CleanData/Guideline Creation/Somerset/Som_ArableRev_G2.tif", overwrite=TRUE)





##-------------------------------------##
#### 11. Group 3: Combine Guidelines ####
##-------------------------------------##

##------------------------##
## Group 3: Combine masks ##
##------------------------##

G3_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
G3_mask <- classify(G3_mask, rclmat, include.lowest=TRUE)
plot(G3_mask)
writeRaster(G3_mask, "CleanData/Guideline Creation/Somerset/Som_MasksAll_G3.tif", overwrite=TRUE)


##------------------------------------##
## Group 3: Combine better guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Bett_G3 <- WaterAbstr + DisurbanceIndex + RWLA + AESFields + (SiltSoil/2) + Hydro_HeightLow + PeatEdge
names(Bett_G3) <- "Better"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Bett_G3) + labs(fill = "Grading") +
  scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G3, "CleanData/Guideline Creation/Somerset/Som_Better_G3.tif", overwrite=TRUE)


##------------------------------------##
## Group 3: Combine bigger guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Big_G3 <- HistAbund + WaterAbstr + DisurbanceIndex + RWLA + AESFields + (SiltSoil/2) + Hydro_HeightLow + PeatEdge
names(Big_G3) <- "Bigger"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Big_G3) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G3, "CleanData/Guideline Creation/Somerset/Som_Bigger_G3.tif", overwrite=TRUE)


##----------------------------------##
## Group 3: Combine more guidelines ##
##----------------------------------##

## Add together all rules and assign name to layer

More_G3 <- HistAbund  + WaterAbstr + DisurbanceIndex + RWLA + AESFields + (SiltSoil/2) + Hydro_HeightLow + PeatEdge
names(More_G3) <- "More"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=More_G3) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G3, "CleanData/Guideline Creation/Somerset/Som_More_G3.tif", overwrite=TRUE)


##----------------------------------------------##
## Group 3: Combine agri-conversion guidelines ##
##---------------------------------------------##

## Add together all rules and assign name to layer
ArableConv_G3 <- CropsMaize + PeatSoils
names(ArableConv_G3) <- "ArableConv"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=ArableConv_G3) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableConv_G3, "CleanData/Guideline Creation/Somerset/Som_ArableRev_G3.tif", overwrite=TRUE)







##-------------------------------------------------##
#### 12.0 Plot graded maps with opportunity areas ####
##-------------------------------------------------##

##---------------------------------##
## 12.1 Prepare Landscape Canvas ####
##---------------------------------##

## Read in the somerset canvas with all land parcel polygons
Canvas <- st_read("CleanData/Scenarios/1-Starting Canvas/Som_Canvas.shp")
plot(Canvas$geometry)


##-- Label lowland wet grassland fields --##

## Read in a raster icreaetd of lowland wet grassland extent in the somerset levels
LowGrass <- rast("RawData/LowlandWetGrassRasters/Som_LowWetGrass.tif")
plot(LowGrass)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(LowGrass, vect(Canvas), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 25% covered by lowland wet grassland pixels
Canvas <- Canvas |>
          mutate(WetGrassProp = Scores$layer) |>
          filter(WetGrassProp > 0.50) |> select(-WetGrassProp)
plot(Canvas$geometry)

## Crop the canvas to just the piority landscape
Canvas <- Canvas |> st_intersection(SomOutline)
plot(Canvas$geometry)


##-- Extract the grades/masks for each of the Lawton strategies --##

## More strategy
MoreGrades <- extract(More_G1, Canvas, fun = mean, na.rm = T)
MoreGrades2 <- extract(More_G2, Canvas, fun = mean, na.rm = T)
MoreGrades3 <- extract(More_G3, Canvas, fun = mean, na.rm = T)
Canvas$MoreGrade_G1 <- MoreGrades$More
Canvas$MoreGrade_G2 <- MoreGrades2$More
Canvas$MoreGrade_G3 <- MoreGrades3$More

## Bigger strategy
BiggerGrades <- extract(Big_G1, Canvas, fun = mean, na.rm = T)
BiggerGrades2 <- extract(Big_G2, Canvas, fun = mean, na.rm = T)
BiggerGrades3 <- extract(Big_G3, Canvas, fun = mean, na.rm = T)
Canvas$BigGrade_G1 <- BiggerGrades$Bigger
Canvas$BigGrade_G2 <- BiggerGrades2$Bigger
Canvas$BigGrade_G3 <- BiggerGrades3$Bigger

## Better Strategy
BetterGrades <- extract(Bett_G1, Canvas, fun = mean, na.rm = T)
BetterGrades2 <- extract(Bett_G2, Canvas, fun = mean, na.rm = T)
BetterGrades3 <- extract(Bett_G3, Canvas, fun = mean, na.rm = T)
Canvas$BetterGrade_G1 <- BetterGrades$Better
Canvas$BetterGrade_G2 <- BetterGrades2$Better
Canvas$BetterGrade_G3 <- BetterGrades3$Better


## Add on cluster size to the Better grading
PopSize <- PopSize |> select(ClustPop, F_LOC_ID) |> mutate(ClustPop = Inv_scale_vals(ClustPop)) # extract population cluster sizes
Canvas <- left_join(Canvas, PopSize, by = "F_LOC_ID") # add these pop sizes to the fields (if they are in a cluster)
## Do the final addition to work out grade
Canvas <- Canvas |> mutate(ClustPop2 = ifelse(is.na(ClustPop)==T, 0, ClustPop), 
                           BetterGrade_G1 = BetterGrade_G1+ClustPop2,
                           BetterGrade_G2 = BetterGrade_G2+ClustPop2,
                           BetterGrade_G3 = BetterGrade_G3+ClustPop2) |> 
                    select(-ClustPop2) # Add inverse scaled clust pop size to the better grading 

                    
## Which fields need to be masked
## First re-label raster to pixels to be masked have value of 0
SomMasks <- G1_mask
values(SomMasks) <- ifelse(is.na(values(SomMasks))==T, 0, values(SomMasks))
SomMasks2 <- G2_mask
values(SomMasks2) <- ifelse(is.na(values(SomMasks2))==T, 0, values(SomMasks2))
SomMasks3 <- G3_mask
values(SomMasks3) <- ifelse(is.na(values(SomMasks3))==T, 0, values(SomMasks3))

## Now extract the average value for each land parcel, fully masked fields will have a score of 1
OppScore <- extract(SomMasks, Canvas, fun = mean, na.rm = T)
OppScore2 <- extract(SomMasks2, Canvas, fun = mean, na.rm = T)
OppScore3 <- extract(SomMasks3, Canvas, fun = mean, na.rm = T)
Canvas$Mask_G1 <- OppScore$layer
Canvas$Mask_G2 <- OppScore2$layer
Canvas$Mask_G3 <- OppScore3$layer



##-----------------------------##
## 12.2 General plot styling ####
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

## Set the plot extent so that all plots have the same area no matter if they have different landarcels
PlotExt <- coord_sf(xlim = c(ext(Canvas)[1]-20, ext(Canvas)[2]+20), ylim = c(ext(Canvas)[3]-20, ext(Canvas)[4]+20), 
                    crs = 27700, expand = FALSE) 



##-------------------------------------##
## 12.3 Group 1: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasBetter <- filter(Canvas, (Mask_G1 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Somerset Levels: Better Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G1_Better.png", width = 20, height = 20, units = "cm")
      


##-------------------------------------##
## 12.4 Group 1: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasBig <- filter(Canvas, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBig, mapping=aes(geometry=geometry, fill = BigGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Somerset Levels: Bigger Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G1_Bigger.png", width = 20, height = 20, units = "cm")



##-----------------------------------##
## 12.5 Group 1: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasMore <- filter(Canvas, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasMore, mapping=aes(geometry=geometry, fill = MoreGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Somerset Levels: More Grading G1") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G1_More.png", width = 20, height = 20, units = "cm")



##---------------------------------------##
## 12.6 Group 1: Plot Map of Landscape ####
##---------------------------------------##

## Read in priority landscape boundary
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
Som <- MyBoxes |> filter(Att3Value == "Somerset Levels and Moors") 

## Read in RSPB reserve outlines
Reserves <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")
ReservesSom <- st_crop(Reserves, (Som |> st_buffer(dist=500)))

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=Som, mapping=aes(geometry=geometry), fill = "lightblue", colour = NA) + 
  geom_sf(data=ReservesSom, mapping=aes(geometry=geometry), fill = "red", alpha = 0.5, colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## Add North arrow and scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  # annotation_north_arrow(location = "br", which_north = "true",
  #                        pad_x = unit(0.13, "in"), pad_y = unit(0.25, "in"),
  #                        style = north_arrow_orienteering,
  #                        height = unit(0.7, "cm"), width = unit(0.5, "cm"),) +
  annotate("text", x = 334849, y = 122571, label = "West Sedgemoor") +
  annotate("text", x = 338849, y = 137571, label = "Greylake") +
  ## set labels
  ggtitle("Somerset Levels") +
  ## set them
  theme_light() + 
  GeneralThemeing +
  theme(legend.position = "none")

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_LandscapeMap.png", width = 20, height = 20, units = "cm")




##-------------------------------------##
## 12.7 Group 2: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasBetter <- filter(Canvas, (Mask_G2 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Somerset Levels: Better Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G2_Better.png", width = 20, height = 20, units = "cm")
      


##-------------------------------------##
## 12.8 Group 2: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasBig <- filter(Canvas, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBig, mapping=aes(geometry=geometry, fill = BigGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Somerset Levels: Bigger Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G2_Bigger.png", width = 20, height = 20, units = "cm")



##-----------------------------------##
## 12.9 Group 2: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasMore <- filter(Canvas, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasMore, mapping=aes(geometry=geometry, fill = MoreGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Somerset Levels: More Grading G2") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G2_More.png", width = 20, height = 20, units = "cm")



##--------------------------------------##
## 12.10 Group 3: Plot better grading ####
##--------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasBetter <- filter(Canvas, (Mask_G3 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G3), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Somerset Levels: Better Grading G3") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G3_Better.png", width = 20, height = 20, units = "cm")
      


##--------------------------------------##
## 12.11 Group 3: Plot bigger grading ####
##--------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasBig <- filter(Canvas, (Mask_G3 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasBig, mapping=aes(geometry=geometry, fill = BigGrade_G3), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Somerset Levels: Bigger Grading G3") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G3_Bigger.png", width = 20, height = 20, units = "cm")



##------------------------------------##
## 12.12 Group 3: Plot more grading ####
##------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasMore <- filter(Canvas, (Mask_G3 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() + 
  ## add landscape outline
  geom_sf(data=SomOutline, mapping=aes(geometry=geometry), colour = "grey", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasMore, mapping=aes(geometry=geometry, fill = MoreGrade_G3), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Somerset Levels: More Grading G3") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Somerset_G3_More.png", width = 20, height = 20, units = "cm")







