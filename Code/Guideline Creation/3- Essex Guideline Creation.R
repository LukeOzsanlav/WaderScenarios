##----------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 09/05/2024
##
## Aim: Create all guidelines for scenario modelling of the Essex Landscape
## 
##----------------------------------------------------------## 


## Load in packages
pacman::p_load(tidyverse, sf, terra, tidyterra, gstat, stars, ggnewscale, basemaps,
               automap, spatstat, ggspatial, leastcostpath, data.table)
options(scipen=999) # turn off scientific notation

## Load in various helper functions
source("Code/Helper functions.R")



##----------------------------------##
#### 0. Read in general data sets ####
##----------------------------------##

## Read in my priority landscapes boundary and buffer it by 3000m
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
EssOutline <- MyBoxes |> filter(Att3Value == "Essex")
Ess <- MyBoxes |> filter(Att3Value == "Essex") |> st_buffer(dist=3000)
rm(MyBoxes) # free up space

## Read in UKCEH landcover and then crop and mask it with the landscape boundary
## This will be used as the canvas to map all other variables across
LCfull <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]]
LC_Ess <- crop(LCfull, vect(Ess)) |> mask(vect(Ess))
rm(LCfull)
plot(LC_Ess) # free up space




##-------------------------------##
#### 1. Better only guidelines ####
##-------------------------------##

##*-- Target smallest existing populations --*##

## Read in full BWWM data set that has each wader cluster labelled
Reg_Fields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Essex")

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

## crop to just cluster in the Essex priority landscapes
Popinter <- st_intersection(st_as_sf(PopSize), EssOutline)
PopSize <- filter(PopSize, F_LOC_ID %in% Popinter$F_LOC_ID)
ggplot(data= PopSize) + geom_sf(mapping = aes(geometry=geometry, fill = ClustPop), colour = NA) + theme_light()

## Write out this variable for later use
write_csv(PopSize, "CleanData/Guideline Creation/Essex/Ess_ClusterPopSize.csv")





##-------------------------------##
#### 2. Bigger only guidelines ####
##-------------------------------##

##*-- Target areas near sites with more breeding waders for expansion --*##

## Loop through each cluster
for(j in 1:length(unique(PopSize$ClustGroup))){

   ## Creat raster of buffered cluster where in pixel inside the buffered cluster recieves the value equal to the number of breeding wader pairs  
   Subbuf <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 2000)), LC_Ess, field = "ClustPop", fun = "max",  background = NA)
   
   ## Create the same raster as above but this time do it for the unbuffered cluster (this is so we can calcualte the distance to the cluster)
   Subno <- rasterize(vect(st_buffer(st_as_sf(filter(PopSize, ClustGroup==j)), dist= 1)), LC_Ess, field = "ClustPop", fun = "max",  background = NA)
   
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

##*-- Target large site creation on the Dengie peninsula --*##

## Read in shapefile with Dengie polygon in
ThamesBlocks <- st_read("RawData/Thames/Thames Regional Polygons.shp")
Dengie <- filter(ThamesBlocks, Area == "Dengie")

## Turn the Dengie block into raster
Dengie <- rasterize(vect(Dengie), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by the polygon then assign a pixel value of 1, using my function
Dengie <- Half_reclass(Dengie)
plot(Dengie)

## free up space
rm(ThamesBlocks); gc()



##*-- Target larger site creation further from existing populations --*##

## Create the same raster from the polygons of wader clusters (this is so we can calculate the distance to the cluster)
Subno <- rasterize(vect(st_as_sf(PopSize)), LC_Ess, cover = T, background = NA)
plot(Subno)

## Calculate the distance of all pixels to the nearest cluster
## then mask it so that I only retain pixels that are within the priority landscape
ClusterD <- distance(Subno, unit="m") |> mask(vect(Ess))
values(ClusterD) <- scale_vals(values(ClusterD))
values(ClusterD) <- ifelse(is.na(values(ClusterD))==T, 0, values(ClusterD))
plot(ClusterD)

## free up space 
rm(Subno); gc()




##*-- Target areas that link existing wader sites --*##

## Read in polygons of the wader sites
WaderSites <- st_read("CleanData/Scenarios/2-DefineWaderSites/Essex/Essex_Wader_Sites.shp")
WaderSites <- WaderSites |> st_centroid()
plot(WaderSites$geometry)

## Read in raster of potential lowland wet grassland and suitable arable land
## Then combine so have tiered values
LowGrass <- rast("RawData/LowlandWetGrassRasters/Essex_LowWetGrass.tif")
ArableSuit <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_ArableSuitable.tif") |> resample(LowGrass)
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
LC_Ess2km <- aggregate(LC_Ess, fact = 80)

## Rasterize the lines connecting polygons using my 2km raster as a base
Connectr <- rasterize(vect(LCP), LC_Ess2km, field = "Val", fun = "sum",  background = 0)
ggplot() + geom_spatraster(data = Connectr) + geom_sf(data = WaderSites, colour = "red", fill = NA)

## Scale the values from the rasterized lines
values(Connectr) <- scale_vals(values(Connectr))
plot(Connectr)

## Finally resample the raster back down to 25m
Connectr <- resample(Connectr, LC_Ess)

## Finally turn any NAs back to zeros
values(Connectr) <- ifelse(is.na(values(Connectr))==T, 0, values(Connectr))
plot(Connectr)

## free up space
rm(WaderSites, LowGrass, Cond, LCP, LC_Ess2km); gc()



##*-- Target area D --*##

## Read in the shapefile of the digitised Tjames areas, These were drwan by wrokshop participants
Areas <- st_read("RawData/Thames/Thames Regional Polygons.shp")
plot(Areas$geometry)

## filter out the region of interest
AreaD <- filter(Areas, Area == "Region D")

## Turn AreaD block into raster
AreaD <- rasterize(vect(AreaD), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by the polygon then assign a pixel value of 1, using my function
AreaD <- Half_reclass(AreaD)
plot(AreaD)

## free up space
rm(Areas); gc()







##-------------------------------##
#### 4. Bigger/More guidelines ####
##-------------------------------##

##*-- Target areas surrounded by less urban areas --*##

## Reclassify land cover raster so that it indicates if pixel is woodland or not
LC_EssWid <- rast("RawData/LandCover/gblcm25m2021.tif")[[1]] |>  crop(vect(Ess))
m <- c(0.0, 19.5, 0,
       19.5, 23, 1) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Urban <- classify(LC_EssWid, rclmat, include.lowest=TRUE)
plot(Urban)

## smooth the raster with focal window of 41 (9*25m = 1.025km)
## then calculate inverse scaled values as we want areas with lower tree cover
UrbanDens <- focal(Urban, w=41, fun="mean", na.rm = T)
values(UrbanDens) <- Inv_scale_vals(values(UrbanDens))
values(UrbanDens) <- ifelse(is.na(values(UrbanDens))==T, 0, values(UrbanDens))
plot(UrbanDens) # plot smoothed raster
rm(LC_EssWid, Urban); gc() # free up space



##*-- Target SSSI area --*##

## Read in shapefile of SSSI areas
SSSI <- st_read("RawData/SSSI/data/Sites_of_Special_Scientific_Interest_England.shp")
SSSI <- st_crop(SSSI, ext(LC_Ess)) # crop to save space

## Rasterize the SSSI polygons
SSSI <- rasterize(vect(SSSI), LC_Ess, cover = T, background = 0)
SSSI <- Half_reclass(SSSI)
values(SSSI) <- Inv_scale_vals(values(SSSI))
InvSSSI <- SSSI
ggplot() + geom_spatraster(data=SSSI) + scale_fill_viridis_c() + theme_light()

## free up space
rm(SSSI); gc()



##*-- Target areas with more water in future --*##

## Read in data set on water abstraction availability
WaterAbst <- st_read("RawData/Thames/Water Availability/Resource_Reliability_pct_of_timePolygon.shp")
WaterAbst$ID <- 1:nrow(WaterAbst) # giev unique row ID

## Filter out the polygons that at least partly fall within the priority landscape
Sub <- st_intersection(WaterAbst, Ess)
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
WaterAbstr <- rasterize(vect(WaterAbst), LC_Ess, field = "PctVal", fun = "min",  background = 0)
## Scale values and change any NAs to zeros
values(WaterAbstr) <- scale_vals(values(WaterAbstr))
values(WaterAbstr) <- ifelse(values(WaterAbstr) == "NaN", 0, values(WaterAbstr))
plot(WaterAbstr)

## free up space
rm(WaterAbst); gc()



##*-- Target areas less likely to be targeted for salt marsh creation --*##

## Read in the two shape file layer, one for opportunities and one for 
SaltOpp <- st_read("RawData/Saltmarsh Creation Sites/SustainableShores_all_potential_sites.shp") |> st_transform(crs = st_crs(LC_Ess))
SaltPr <- st_read("RawData/Saltmarsh Creation Sites/SustainableShores_priority_sites.shp") |> st_transform(crs = st_crs(LC_Ess))

## rasterize the opportunity areas with a value of 0.5
SaltOpp$Val <- 0.5
SaltOpp <- rasterize(vect(SaltOpp), LC_Ess, field = "Val", fun = "max",  background = 1)
plot(SaltOpp)

## rasterize the opportunity areas with a value of 1
SaltPr$Val <- 0
SaltPr <- rasterize(vect(SaltPr), LC_Ess, field = "Val", fun = "max",  background = 1)
plot(SaltPr)

## Take the max value of the two rasters
SaltmarshOpp <- min(SaltOpp, SaltPr)
plot(SaltmarshOpp)

## free up space
rm(SaltOpp, SaltPr); gc()



##*-- Target area A and B --*##

## Read in the shapefile of the digitised Tjames areas, These were drwan by wrokshop participants
Areas <- st_read("RawData/Thames/Thames Regional Polygons.shp")
plot(Areas$geometry)

## filter out the region of interest
AreaAB <- filter(Areas, Area %in% c("Region A", "Region B1", "Region B2"))

## Turn AreaAB block into raster
AreaAB <- rasterize(vect(AreaAB), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by the polygon then assign a pixel value of 1, using my function
AreaAB <- Half_reclass(AreaAB)
plot(AreaAB)

## free up space
rm(Areas); gc()



##*-- Target areas with wintering waterbird CS agreements already --*##

## Read in the Different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")
ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
ESS <- rbind(ESSWest, ESSEast)
rm(ESSWest, ESSEast)

## Read in the Canvas as this has the field parcel shapes 
Canv <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp") |> select(geometry)

## Filter out winter water birds focused options and breeding waders water options
ESS_OpsWa <- filter(ESS, optcode %in% c("HK9", "HK11", "HK13", "HK10", "HK12", "HK14", "DR","WDC","HK15","HK16","HK17","EK2","OK2","HK2","OHK2","EK3","OK3",
                                        "HK3","OHK3","EK4","OK4","HK4","OHK4","HK6","HK7","HK8","HL8","HQ3","HQ4","HQ5","HQ6HQ7","HQ8","HR6","EK21","HK21","OHK21"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
ESS_InterWa <- as.data.frame(st_covers(Canv, ESS_OpsWa))
Canv$ESS_Wader <- "N"
Canv$ESS_Wader[ESS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()


## Filter out winter waterbird or breeding focused options
CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12", "WN3", "WN4", "GS2", "GS5", "GS6", "GS7", 
                                             "GS8", "UP2", "WT6", "WT7", "WT8", "WT9", "GS4", "GS13", "GS14", "GS16"))

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
AESWader <- rasterize(vect(Canv), LC_Ess, cover = T, background = 0)
AESWader <- Half_reclass(AESWader)
plot(AESWader, main = "Wader AES fields")

## free up space
rm(Canv, ESS, ESS_InterWa, CSS_InterWa, CSS_Ops, CSS_OpsWa, ESS_OpsWa); gc()



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
SiltSoil <- rasterize(vect(SiltSoils), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by silty soils then assign a pixel value of 1, using my function
SiltSoil <- Half_reclass(SiltSoil)
plot(SiltSoil)
rm(Soil, SiltWetTypes, SiltSoils); gc()



##*-- Target areas away from the coastal footpath --*##

## Read in the coastal path shapefile
CoastPath <- st_read("RawData/England_Coast_Path_Route/England_Coast_Path_RouteLine.shp")

## Extract the section in the prioirty landscape
EssexPath <- st_crop(CoastPath, Ess)
plot(EssexPath$geometry)

## Buffer by 500m becasue...
## Black-tailed Godwit paper found effect of disturbance up to 500m away
## https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1474-919X.2008.00889.x
EssexPath <- st_buffer(EssexPath, dist = 500)
plot(EssexPath$geometry)


## Rasterize the buffered coastal path
## Need to take inverse so that areas within buffer have value of 0 and all other areas have avlue of 1
EssexPath <- rasterize(vect(EssexPath), LC_Ess, cover = T, background = 0)
EssexPath <- Half_reclass(EssexPath)
values(EssexPath) <- Inv_scale_vals(values(EssexPath))
plot(EssexPath, main = "Coatal path buffered")



##*-- Target areas further away from open tips/landfills --*##

## Read in the shapefile that I edited, this has landfill sites that:
## had the following type `Waste Landfilling; >10 T/D With Capacity >25,000T Excluding Inert Waste -  5.2 A(1) a)` 
## and where the status was `effective` and they looked active on satellite view
Landfill <- st_read("RawData/Landfill Active/Edited Active None Inert Landfill Sites.shp") |> 
             st_crop(st_buffer(EssOutline, dist = 20000))
plot(Landfill$geometry)

## extend the landcover raster so that we can rasterize landfill sites that are outside the priority landscape
Ext_LC <-  extend(LC_Ess, y = 1000)
plot(Ext_LC)

## now rasterize the landfill sites
Landfillr <- rasterize(vect(Landfill), Ext_LC, cover = T, background = NA)
plot(Landfillr)

## calcualte the distance to the neareast landfill sites
LandfillDist <- distance(Landfillr)
plot(LandfillDist)

## Now reample this larger raster back to to the extent of the the Essex land cover raster
LandfillDist <- resample(LandfillDist, LC_Ess) 
LandfillDist <- mask(LandfillDist, EssOutline) # mask to retain just priotiy landscape
plot(LandfillDist)

## Scale values and turn NAs to zero
values(LandfillDist) <- scale_vals(values(LandfillDist))
values(LandfillDist) <- ifelse(is.na(values(LandfillDist))==T, 0, values(LandfillDist))
plot(LandfillDist)

## free up space
rm(Landfill, Ext_LC); gc()



##*-- Target fields with lower public footfall --*##

# ## Read in the shapefiles of the paths from the OrVal data set
# ## PID links that path lines to the attribute data
# Paths <- st_read("RawData/ORVal/Paths/paths_england.shp")
# 
# ## Read in the attribute data for the access points and streamline to retain just visitor numbers
# Simple <- read_csv("RawData/ORVal/Thames Simple Data/ORVal Fri, 09 Feb 2024 13_56_25 GMT.csv")
# Simple <- Simple |> select(pid, vis)
# 
# ## Crop the paths to just the region of interest
# Paths <- st_intersection(Paths, Ess) |> select(PID)
# plot(Paths$geometry)
# 
# ## Join paths with attributes data, a few paths are lost as they did not have attribute data
# PathsAtt <- inner_join(Paths, Simple, by = c("PID"= "pid"))
# ggplot() + geom_sf(data = Ess) +  geom_sf(data=st_buffer(PathsAtt, dist = 500), mapping=aes(colour=vis)) + theme_light() # 
# 
# 
# ## Calculate the distance between all pixels and the nearest path
# PathDist <- distance(LC_Ess, vect(PathsAtt), rasterize=T)
# PathDist <- mask(PathDist, vect(Ess))
# plot(PathDist)
# 
# ## Only retain path distances within 500m of a path
# ## This was the distance that the below paper found that human presence could impact breeding waders
# ## https://onlinelibrary.wiley.com/doi/10.1111/j.1474-919X.2008.00889.x
# PathsAttBuf <- st_buffer(PathsAtt, dist = 500)
# PathsAttBuf <- st_combine(PathsAttBuf)
# PathDist <- mask(PathDist, vect(PathsAttBuf))
# plot(PathDist)
# 
# ## For each path buffered by 500m calcualte the max visitor value for each pixel
# PathsAttBuf <- st_buffer(PathsAtt, dist = 500)
# for(j in 1:nrow(PathsAttBuf)){
#   message(j)
#   PathsAttBuf1 <- PathsAttBuf[j,]
#   PathValsr <- rasterize(vect(PathsAttBuf1), LC_Ess, field = "vis", fun = "max",  background = 0)
#   if(j==1){AllPathValsr <- PathValsr}else{AllPathValsr <- max(AllPathValsr, PathValsr)}
# }
# plot(AllPathValsr)
# 
# ## Take inverse of distance to paths so pixels closer to path have higher values
# values(PathDist) <- Inv_scale_vals(values(PathDist))
# plot(PathDist)
# 
# ## Now multiple the inverse distance to path by the max visitor value to create an index of disturbance
# DisurbanceIndex <- PathDist*AllPathValsr
# plot(DisurbanceIndex)
# 
# ## finally take inverse and change all NAs to 1
# ## Inverse insures that all low value disturbed pixels have a higher grade
# ## All NA pixels need to be 1 as these are not in proximity to a path
# values(DisurbanceIndex) <- Inv_scale_vals(values(DisurbanceIndex))
# values(DisurbanceIndex) <- ifelse(is.na(values(DisurbanceIndex))==T, 1, values(DisurbanceIndex))
# plot(DisurbanceIndex)
# rm(PathDist, PathsAttBuf, PathsAttBuf1, AllPathValsr, PathValsr, PathsAtt, Paths, Simple); gc()





##-------------------------------------##
#### 5. Arable conversion guidelines ####
##-------------------------------------##

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
AgriGrades <- rasterize(vect(AgriGrade), LC_Ess, field = "Grade", fun = "min", background = 0)
values(AgriGrades) <- scale_vals(values(AgriGrades))
plot(AgriGrades)
rm(AgriGrade); gc()


##*-- Avoid Area C --*##

## Read in the shapefile of the digitised Tjames areas, These were drwan by wrokshop participants
Areas <- st_read("RawData/Thames/Thames Regional Polygons.shp")
plot(Areas$geometry)

## filter out the region of interest
AreaC <- filter(Areas, Area %in% c("Region C"))

## Turn AreaC block into raster
AreaC <- rasterize(vect(AreaC), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by the polygon then assign a pixel value of -1, using my function
AreaC <- Half_reclass(AreaC)
values(AreaC) <- Inv_scale_vals(values(AreaC)) 
plot(AreaC)

## free up space
rm(Areas); gc()



##*-- Target areas closer to the sea wall --*##

## Read in a shapefile of the coastline
## I have added a section of additional coastline that wasn't included in South Essex along th Thames
Shoreline <- st_read("RawData/Shoreline Management/Shoreline_Management_Plan_MappingLine.shp")
Shoreline <- st_crop(Shoreline, st_buffer(Ess, dist = 5000))
plot(Shoreline$geometry)

## Calculate distance of all pixels from the coastline and then mask by the priority landscape boundary
distance_raster <- distance(LC_Ess, vect(Shoreline))
distance_raster <- mask(distance_raster, vect(EssOutline))
plot(distance_raster)

## Now scale the distance to coast values
## For any pixels over 5km from the coast just reduce these down to 5000
## Areas on the seawall should have a value of 1
CoastDist <- distance_raster
values(CoastDist) <- ifelse(values(CoastDist)> 5000, 5000, values(CoastDist))
values(CoastDist) <- Inv_scale_vals(values(CoastDist))
values(CoastDist) <- ifelse(is.na(values(CoastDist))==T, 0, values(CoastDist))
plot(CoastDist)

## free up space
rm(Shoreline, distance_raster); gc()



##*-- Target areas near fields with arable Lapwing plots --*##

## Read in the Different Stewardship Schemes
CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")

## Read in the Canvas as this has the field parcel shapes 
Canv <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp") |> select(geometry)

## Filter out winter waterbird or breeding focused options
CSS_LapPlot <- filter(CSS_Ops, opt_code %in% c("AB5"))

## Since these CSS options are points, determine if a points falls within one of the survey fields
## And label each field whether a wader specific CSS point fell within it
CSS_InterWa <- as.data.frame(st_covers(Canv, CSS_LapPlot))
Canv$CSS_Wader <- "N"
Canv$CSS_Wader[CSS_InterWa$row.id] <- "Y"

## plots to check
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = CSS_Wader), colour = NA) +
           theme_minimal()

## Rasterize this AES field polygons
Canv <- filter(Canv, CSS_Wader == "Y")
PlotFields <- rasterize(vect(Canv), LC_Ess, cover = T, background = NA)
plot(PlotFields, main = "Lapwing Plots")

## calcualte distance to Lapwing plot within the priority landscapes
## Change distacnes over 4km back to 4km
## This is based off this paper where most local breeding dispersals are within 4km https://ornisfennica.journal.fi/article/view/133716/82266
DistLapPlot <- distance(PlotFields) |> mask(EssOutline)
values(DistLapPlot) <- ifelse(values(DistLapPlot) > 4000, 4000, values(DistLapPlot))

## Take inverse scaled values of raster 
values(DistLapPlot) <- Inv_scale_vals(values(DistLapPlot))
values(DistLapPlot) <- ifelse(is.na(values(DistLapPlot))==T, 0, values(DistLapPlot))
plot(DistLapPlot)


## free up space
rm(Canv, CSS_InterWa, CSS_Ops, CSS_LapPlot); gc()



##--* Target arable land with poor soil health --*##

## Going to use susceptibility to erosion as a proxy for soil health

## Read in water and wind erosion data set
WindEr <- rast("RawData/Soil Erosion/ILSWE_UK_wind.tif") |> project(y=vect(EssOutline)) |> crop(vect(EssOutline)) |> mask(vect(EssOutline))
plot(WindEr)
WaterEr <- rast("RawData/Soil Erosion/RUSLE_2016_UK_water.tif") |> project(y=vect(EssOutline)) |> crop(vect(EssOutline)) |> mask(vect(EssOutline))
plot(WaterEr)


## Turn the wind erosion raster from a categorical raster to a numeric raster
WindEr <- as.numeric(WindEr)
values(WindEr) <- ifelse(values(WindEr)==max(values(WindEr), na.rm=T), NA, values(WindEr)) ## some NAs take on a large value so remove these
values(WindEr) <- Inv_scale_vals(values(WindEr)) ## take inverse, so areas with high erosion have highest score
plot(WindEr)

## scale the values in the water erosion raster then reamsple so I can add it to the wind erosion raster
values(WaterEr) <- scale_vals(values(WaterEr))
WaterEr <- resample(WaterEr, WindEr)
plot(WaterEr)

## Add together the two raster layers, the higher values are the most at risk from erosion
ErosionIndex <- WindEr+WaterEr
plot(ErosionIndex)


## Created smoothed rasters to fill in some of the NAs in the raster
## do this will a 3x3 square and 5x5 square
ErosionIndexSmooth3 <- focal(x=ErosionIndex, w=3, fun="mean", na.rm=T)
plot(ErosionIndexSmooth3)

ErosionIndexSmooth5 <- focal(x=ErosionIndex, w=5, fun="mean", na.rm=T)
plot(ErosionIndexSmooth5)

## Replace any NAs in the original raster with the 3x3 smoothed raster first 
## The replace any remaining NAs using the 5x5 raster
values(ErosionIndex) <- ifelse(is.na(values(ErosionIndex))==T, values(ErosionIndexSmooth3), values(ErosionIndex))
plot(ErosionIndex)
values(ErosionIndex) <- ifelse(is.na(values(ErosionIndex))==T, values(ErosionIndexSmooth5), values(ErosionIndex))
plot(ErosionIndex)

## finally simple scale the values
values(ErosionIndex) <- scale_vals(values(ErosionIndex)) 

## Finally re-sample this coarse raster to the same resolution as my base raster
## and change NAs to 0's
ErosionIndex <- resample(ErosionIndex, LC_Ess)
values(ErosionIndex) <- ifelse(is.na(values(ErosionIndex))==T, 0, values(ErosionIndex))
plot(ErosionIndex)

## free up space
rm(WindEr, WaterEr, ErosionIndexSmooth3, ErosionIndexSmooth5); gc()




##--------------------------------##
#### 6. All strategy guidelines ####
##--------------------------------##

##*-- Target lowest land within hydro units --*##

## 2m LiDAR raster for the Essex
## and crop then mask to the priority landscape
Ess_Elev <- rast("RawData/LiDAR/Thames_DTM_2m.tif") |> crop(vect(Ess)) |> mask(vect(Ess))
Ess_Elev <- aggregate(Ess_Elev, fact = 4)
plot(Ess_Elev)

## If any values are over 30m in elevation then set these to 30m 
## This is done as any values of 30m are not suitable and should all ultimately have a value of 0
gc()
values(Ess_Elev) <- ifelse(values(Ess_Elev)>30, 30, values(Ess_Elev))
plot(Ess_Elev)
gc()

## Take the inverse scale of the elevation, so the lowest areas are given a value of 1 and anywhere 30m or over 0
values(Ess_Elev) <- Inv_scale_vals(values(Ess_Elev))
plot(Ess_Elev)

## Re-sample the raster so it has the same resolution as the landcover raster
Ess_Elev <- resample(Ess_Elev, LC_Ess)
values(Ess_Elev) <- ifelse(is.na(values(Ess_Elev))==T, 0, values(Ess_Elev)) # change NAs to 0
plot(Ess_Elev)




##*-- Target sites with access to streams or boreholes --*##

## Only able to find data on water courses so going to use that
## Currently just have the OS watercourse data set but applied for the higher res UKCEH data
Stream <- st_read("RawData/Open Rivers/oprvrs_essh_gb/data/WatercourseLink.shp") |> 
             st_crop(st_buffer(EssOutline, dist = 20000)) 
plot(Stream$geometry)

## Read in the UK coastline shapefile
## Then crop the streams data set so that only streams actually on land can be used
UK <- st_read("RawData/UK_Coastline/Coastline Higher res/infuse_gb_2011_clipped.shp") |>  
  st_crop(st_buffer(EssOutline, dist = 20000)) |> 
  st_buffer(dist=-25)
Stream <- st_intersection(Stream, UK)
plot(Stream$geometry)


## extend the landcover raster so that we can rasterize landfill sites that are outside the priority landscape
Ext_LC <-  extend(LC_Ess, y = 500)
plot(Ext_LC)

## now rasterize the landfill sites
Stream <- st_buffer(Stream, dist=10) 
Streamr <- rasterize(vect(Stream), Ext_LC, cover = T, background = NA)
plot(Streamr)

## calcualte the distance to the neareast stream
StreamDist <- distance(Streamr)
plot(StreamDist)

## Now reample this larger raster back to to the extent of the the Essex land cover raster
StreamDist <- resample(StreamDist, LC_Ess) 
StreamDist <- mask(StreamDist, EssOutline) # mask to retain just priotiy landscape
plot(StreamDist)

## Scale values and turn NAs to zero
values(StreamDist) <- Inv_scale_vals(values(StreamDist))
values(StreamDist) <- ifelse(is.na(values(StreamDist))==T, 0, values(StreamDist))
plot(StreamDist)

## free up space
rm(Stream, Ext_LC); gc()




##*-- Target areas at lower risk from tidal inundation --*##

## Read in a shape file of the coastline
## I have added a section of additional coastline that wasn't included in South Essex along the Thames
Shoreline <- st_read("RawData/Shoreline Management/Shoreline_Management_Plan_MappingLine.shp")
Shoreline <- st_crop(Shoreline, st_buffer(Ess, dist = 5000))
plot(Shoreline$geometry)

## Now add a new column to the shoreline management plan
## If there is no action (NAI) or management realignment (MR) then these areas are at risk
## If the the polyci is hold the line (HTL) then this is a good and gets the higher score
Shoreline <- Shoreline |> 
             mutate(Policy50 = ifelse(policy_50 == "MR" | policy_50 == "NAI" , 2, 1))


## For each of the the two options classified above calcualte the distacen from each pixel to the nearest bit of coastline
Dist_raster1 <- distance(LC_Ess, vect(Shoreline |> filter(Policy50==1)))
plot(Dist_raster1)
Dist_raster2 <- distance(LC_Ess, vect(Shoreline |> filter(Policy50==2)))
plot(Dist_raster2)

## now created a raster were the pixel value corresponds to the closet policy option in distance
MinDistr <- Dist_raster1
values(MinDistr) <- ifelse(values(Dist_raster1) < values(Dist_raster2), 2, 1)
plot(MinDistr)
rm(Dist_raster1, Dist_raster2);gc()


## Use LiDAR data to extract pixel within 3m of sea level as these may be at risk from inundation
## 2m LiDAR raster for the Essex
## and crop then mask to the priority landscape
Ess_El <- rast("RawData/LiDAR/Thames_DTM_2m.tif") |> crop(vect(Ess)) |> mask(vect(Ess))
Ess_El <- aggregate(Ess_El, fact = 4) ## aggregate to save space
plot(Ess_El)

## If any values are over 3m in elevation then set these to NA 
values(Ess_El) <- ifelse(values(Ess_El)>3, NA, 1)
Ess_El <- resample(Ess_El, MinDistr) ## resample so it can be used as a mask
plot(Ess_El)
gc()


## Now mask my neareast policy raster using the elevation 3m raster
IndRisk <- crop(MinDistr, Ess_El) 
IndRisk <- mask(IndRisk, Ess_El)
values(IndRisk) <- values(IndRisk) -1
values(IndRisk) <- ifelse(is.na(values(IndRisk)) ==T, 1, values(IndRisk))
plot(IndRisk)


## free up space
rm(MinDistr, Ess_El); gc()






##----------------------##
#### 7. Mask Creation ####
##----------------------##


##*--Avoid priority habitats --*##

## Read in the NE Priority Habitat
NE <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")

## Crop out the Natural England priority habitats for the buffered priority landscape
NE <- st_intersection(NE, Ess)

## Get a list of the main habitats and additional habitats
table(NE$mainhabs)
table(NE$addhabs)

## Create vector of habitats that I want to mask out
MAINHABS <- c("Coastal saltmarsh", "Coastal saltmarsh,Saline lagoons",
              "Coastal vegetated shingle",
              "Deciduous woodland", 
              "Lowland calcareous grassland",
              "Lowland dry acid grassland", 
              #"Lowland fens", "Lowland fens,Reedbeds",
              "Lowland heathland",
              "Lowland raised bog",
              "Maritime cliff and slope", "Maritime cliff and slope,Coastal saltmarsh",
              "Mudflats",
              "Reedbeds", "Reedbeds,Coastal saltmarsh",
              "Saline lagoons",
              "Traditional orchard")

## Create column of just the first 5 characters of the additional habitats column
NE <- NE |>  mutate(addhabs = substr(addhabs, start = 1, stop = 5))

## Filter the data set by Main habitat and then if no main habitat retain those that did not have CFPGM in the additional habitat column
NE_Area_Masks <- filter(NE, mainhabs %in% MAINHABS |
                           (mainhabs == "No main habitat but additional habitats present" & !addhabs %in% c("CFPGM", "GQSIG", "LMEAD") ))

## plot these masks
ggplot() + geom_sf(data = NE_Area_Masks, mapping = aes(geometry = geometry, fill = mainhabs), colour = NA) + theme_light()

## Rasterize the priotity habitat polygons
NE_Mask <- rasterize(vect(NE_Area_Masks), LC_Ess, cover = T, background = 0)
## If a pixel is more than 50% covered by the priority habitat then assign a pixel value of 1, using my function
NE_Mask <- Half_reclass(NE_Mask)
plot(NE_Mask)
rm(NE, MAINHABS, NE_Area_Masks); gc() # free up space



##*--Avoid scheduled monuments --*##

## Read in polygons of scheduled monuments from National Heritage
Monument <- st_read("RawData/Monuments/National_Heritage_List_for_England/Scheduled_Monuments.shp")

## Crop out the monuments for the buffered priority landscape
Monument <- st_intersection(Monument, Ess) |> st_buffer(dist = 20)
plot(Monument$geometry)

## Rasterize the monument polygons
Monument_Mask <- rasterize(vect(Monument), LC_Ess, cover = T, background = 0)
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
Urban_mask <- classify(LC_Ess, rclmat, include.lowest=TRUE)
plot(Urban_mask)




##-----------------------------------------##
#### 8. Define strategy opportunity area ####
##-----------------------------------------##

##*-- Better Strategy --*##

## Read in full BWWM data set that has each wader cluster labelled
ClustFields <- st_read("CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp") |> 
              filter(Region == "Essex") |> 
              st_buffer(dist=25)

## Rasterize the wader clusters
Better_Opp <- rasterize(vect(ClustFields), LC_Ess, cover = T, background = 0)
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
BigMore_Opp <- mask(BigMore_Opp , vect(Ess))
plot(BigMore_Opp)



##*-- Agricultural Reversion --*##

## Reclassify raster so that it indicates if pixel is Urban/Sub-urban/Coastal or not
m <- c(0.0, 2.5, NA,
       2.5, 3.5, 1,
       3.5, 22, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Arable_Opp <- classify(LC_Ess, rclmat, include.lowest=TRUE)
Arable_Opp <- mask(Arable_Opp , vect(Ess))
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
writeRaster(G1_mask, "CleanData/Guideline Creation/Essex/Ess_MasksAll_G1.tif", overwrite=TRUE)

## create a general mask to use in Scenarios later
Gen_mask <- max(NE_Mask, Monument_Mask, Urban_mask, na.rm=T)
m <- c(0, 0.5, 1, 0.5, 1.5, NA) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
Gen_mask <- classify(Gen_mask, rclmat, include.lowest=TRUE)
plot(Gen_mask)
writeRaster(Gen_mask, "CleanData/Guideline Creation/Essex/Ess_MasksAll_Gen.tif", overwrite=TRUE)



##------------------------------------##
## Group 1: Combine better guidelines ##
##------------------------------------##

Bett_G1 <- WaterAbstr + IndRisk
names(Bett_G1) <- "Better"

ggplot() + geom_spatraster(data=Bett_G1) + labs(fill = "Grading") +
  scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G1, "CleanData/Guideline Creation/Essex/Ess_Better_G1.tif", overwrite=TRUE)



##------------------------------------##
## Group 1: Combine bigger guidelines ##
##------------------------------------##

Big_G1 <- ExpandClust + UrbanDens + InvSSSI + WaterAbstr + SaltmarshOpp + EssexPath + LandfillDist + IndRisk
names(Big_G1) <- "Bigger"

ggplot() + geom_spatraster(data=Big_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G1, "CleanData/Guideline Creation/Essex/Ess_Bigger_G1.tif", overwrite=TRUE)



##----------------------------------##
## Group 1: Combine more guidelines ##
##----------------------------------##

More_G1 <- WaterAbstr + IndRisk + UrbanDens + InvSSSI + Dengie + ClusterD + SaltmarshOpp + EssexPath + LandfillDist 
names(More_G1) <- "More"

ggplot() + geom_spatraster(data=More_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G1, "CleanData/Guideline Creation/Essex/Ess_More_G1.tif", overwrite=TRUE)


##----------------------------------------------------##
## Group 1: Combine agri-conversion Bigger guidelines ##
##----------------------------------------------------##

ArableBig_G1 <- AgriGrades + DistLapPlot + ErosionIndex + Big_G1
names(ArableBig_G1) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableBig_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableBig_G1, "CleanData/Guideline Creation/Essex/Ess_ArableBig_G1.tif", overwrite=TRUE)


##--------------------------------------------------##
## Group 1: Combine agri-conversion More guidelines ##
##--------------------------------------------------##

ArableMore_G1 <- AgriGrades + DistLapPlot + ErosionIndex + More_G1
names(ArableMore_G1) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableMore_G1) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableMore_G1, "CleanData/Guideline Creation/Essex/Ess_ArableMore_G1.tif", overwrite=TRUE)





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
writeRaster(G2_mask, "CleanData/Guideline Creation/Essex/Ess_MasksAll_G2.tif", overwrite=TRUE)



##------------------------------------##
## Group 2: Combine better guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Bett_G2 <- Ess_Elev + StreamDist + IndRisk
names(Bett_G2) <- "Better"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Bett_G2) + labs(fill = "Grading") +
  scale_fill_viridis_c(na.value = "lightgrey") + theme_light() 
writeRaster(Bett_G2, "CleanData/Guideline Creation/Essex/Ess_Better_G2.tif", overwrite=TRUE)



##------------------------------------##
## Group 2: Combine bigger guidelines ##
##------------------------------------##

## Add together all rules and assign name to layer
Big_G2 <- UrbanDens + AreaAB + Ess_Elev + AESWader + SiltSoil + StreamDist + IndRisk
names(Big_G2) <- "Bigger"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=Big_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(Big_G2, "CleanData/Guideline Creation/Essex/Ess_Bigger_G2.tif", overwrite=TRUE)



##----------------------------------##
## Group 2: Combine more guidelines ##
##----------------------------------##

## Add together all rules and assign name to layer
More_G2 <- UrbanDens + Connectr + AreaD + AreaAB + Ess_Elev + AESWader + SiltSoil + StreamDist + IndRisk
names(More_G2) <- "More"

## Plot combined grading layer and save as a tif file
ggplot() + geom_spatraster(data=More_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(More_G2, "CleanData/Guideline Creation/Essex/Ess_More_G2.tif", overwrite=TRUE)


##----------------------------------------------------##
## Group 2: Combine agri-conversion Bigger guidelines ##
##----------------------------------------------------##

ArableBig_G2 <- AgriGrades + AreaC + CoastDist + Big_G2
names(ArableBig_G2) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableBig_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableBig_G2, "CleanData/Guideline Creation/Essex/Ess_ArableBig_G2.tif", overwrite=TRUE)


##--------------------------------------------------##
## Group 2: Combine agri-conversion More guidelines ##
##--------------------------------------------------##

ArableMore_G2 <- AgriGrades + AreaC + CoastDist + More_G2
names(ArableMore_G2) <- "ArableConv"

ggplot() + geom_spatraster(data=ArableMore_G2) + scale_fill_viridis_c(na.value = "lightgrey") + theme_light()
writeRaster(ArableMore_G2, "CleanData/Guideline Creation/Essex/Ess_ArableMore_G2.tif", overwrite=TRUE)





##-------------------------------------------------##
#### 11.0 Plot graded maps with opportunity areas ####
##-------------------------------------------------##

##---------------------------------##
## 11.1 Prepare Landscape Canvas ####
##---------------------------------##


##-- Label lowland wet grassland fields --##

## Read in the essex canvas with all land parcel polygons
CanvasGr <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp")
plot(CanvasGr$geometry)

## Read in a raster I created of lowland wet grassland extent in the Essex
LowGrass <- rast("RawData/LowlandWetGrassRasters/Essex_LowWetGrass.tif")
plot(LowGrass)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(LowGrass, vect(CanvasGr), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 50% covered by lowland wet grassland pixels
CanvasGr <- CanvasGr |>
          mutate(WetGrassProp = Scores$layer) |>
          filter(WetGrassProp > 0.50) |> select(-WetGrassProp)
plot(CanvasGr$geometry)

## Crop the canvas to just the piority landscape
CanvasGr <- CanvasGr |> st_intersection(EssOutline)
plot(CanvasGr$geometry)



##-- Label arable opportunity fields --##

## Read in the essex canvas with all land parcel polygons
CanvasAr <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp")
plot(CanvasAr$geometry)

## Read in a raster I created of lowland wet grassland extent in the Essex
SuitArable <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_ArableSuitable.tif")
plot(SuitArable)

## Extract the average score from the lowland wet grassland raster (1=wet grass. 0=other habitats)
Scores <- terra::extract(SuitArable, vect(CanvasAr), fun = mean, na.rm = T)

## Classify land parcels as wet grass land if it was more than 50% covered by lowland wet grassland pixels
CanvasAr <- CanvasAr |>
          mutate(SuitArabProp = Scores[,2]) |>
          filter(SuitArabProp > 0.50) |> select(-SuitArabProp)
plot(CanvasAr$geometry)

## Crop the canvas to just the piority landscape
CanvasAr <- CanvasAr |> st_intersection(EssOutline)
plot(CanvasAr$geometry)




##-- Extract the grades/masks for each of the Lawton strategies --##

## Code to read in the rules if needed

# G1_mask <- rast("CleanData/Guideline Creation/Essex/Ess_MasksAll_G1.tif")
# Bett_G1 <- rast("CleanData/Guideline Creation/Essex/Ess_Better_G1.tif")
# Big_G1 <- rast("CleanData/Guideline Creation/Essex/Ess_Bigger_G1.tif")
# More_G1 <- rast("CleanData/Guideline Creation/Essex/Ess_More_G1.tif")
# ArableConv_G1 <- rast("CleanData/Guideline Creation/Essex/Ess_ArableRev_G1.tif")
# 
# G2_mask <- rast("CleanData/Guideline Creation/Essex/Ess_MasksAll_G2.tif")
# Bett_G2 <- rast("CleanData/Guideline Creation/Essex/Ess_Better_G2.tif")
# Big_G2 <- rast("CleanData/Guideline Creation/Essex/Ess_Bigger_G2.tif")
# More_G2 <- rast("CleanData/Guideline Creation/Essex/Ess_More_G2.tif")
# ArableConv_G2 <- rast("CleanData/Guideline Creation/Essex/Ess_ArableRev_G2.tif")


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

## Set the plot extent so that all plots have the same area no matter if they have different landarcels
PlotExt <- coord_sf(xlim = c(ext(CanvasAr)[1]-20, ext(CanvasAr)[2]+20), ylim = c(ext(CanvasAr)[3]-20, ext(CanvasAr)[4]+20), 
                    crs = 27700, expand = FALSE) 

## Read in an outline of the UK to put in the background of the maps
Coast <- st_read("RawData/UK_Coastline/UK_Coatline.shp")
## crop to area just around the Broads priority landscape
Coast <- st_transform(Coast, crs = st_crs(EssOutline)) |> st_crop(EssOutline |> st_buffer(dist = 2000))


##-------------------------------------##
## 11.3 Group 1: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasGrBetter <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Essex Coast: Better (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G1_Better.png", width = 20, height = 13, units = "cm")
      


##-------------------------------------##
## 11.4 Group 1: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Essex Coast: Bigger (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G1_Bigger.png", width = 20, height = 13, units = "cm")



##-----------------------------------##
## 11.5 Group 1: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G1), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Essex Coast: More (Conservationists)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G1_More.png", width = 20, height = 13, units = "cm")



##-----------------------------------------##
## 11.6 Group 1: Plot arable big grading ####
##-----------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG1 <- filter(CanvasAr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG1, mapping=aes(geometry=geometry, fill = ArableBig_G1), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Essex Coast: Arable Reversion for Bigger (Conservationists)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG1)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G1_ArableBig.png", width = 20, height = 13, units = "cm")



##------------------------------------------##
## 11.7 Group 1: Plot arable more grading ####
##------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG1 <- filter(CanvasAr, (Mask_G1 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG1, mapping=aes(geometry=geometry, fill = ArableMore_G1), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Essex Coast: Arable Reversion for More (Conservationists)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG1)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G1_ArableMore.png", width = 20, height = 13, units = "cm")



##---------------------------------------##
## 11.7 Group 1: Plot Map of Landscape ####
##---------------------------------------##

## Read in priority landscape boundary
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
Essex <- MyBoxes |> filter(Att3Value == "Essex") 

## Read in RSPB reserve outlines
Reserves <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")
ReservesEss <- st_crop(Reserves, (Ess |> st_buffer(dist=500)))


## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArr <- filter(CanvasAr, (Mask_G1 > 0.5)) |> mutate(Habb = "#ffd166") |> select(Habb)
CanvasGrr <- filter(CanvasGr, (Mask_G1 > 0.5)) |> mutate(Habb = "#06d6a0") |> select(Habb)
Canv <- rbind(CanvasArr, CanvasGrr)

## make plot
ggplot() +
  ## Add in the OS basemap
  basemap_gglayer(Essex |> st_buffer(dist = 4000), map_service = "osm", map_type = "streets") +
  scale_fill_identity() +
  coord_sf(expand = FALSE) +
  
  ## add landscape outline
  geom_sf(data=Essex |> st_transform(crs=3857), mapping=aes(geometry=geometry), colour = "#273746", fill = NA, linewidth = 0.5) +
  coord_sf(expand = FALSE) + 
  
  ## add polygons
  new_scale_fill() +
  geom_sf(data=Canv |> st_transform(crs=3857), mapping=aes(geometry=geometry, fill = Habb), colour = NA, alpha =0.8) +
  coord_sf(expand = FALSE) + 
  scale_fill_manual(name = "Opportunity\nHabitat",   # Change legend title
                    values = c("#06d6a0", "#ffd166"),
                    labels = c("Grassland", "Arable")) +

  ## Add North arrow and scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  ## set labels
  ggtitle("Essex Coast: Opportunity Habitat") +
  ## set them
  theme(legend.position = "right", 
        axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
        text = element_text(color = "#2D2D2E"), 
        panel.grid = element_line(color = "#ebebe5", linewidth = 0.2),
        panel.background = element_rect(fill = "#f5f5f2", color = NA),
        axis.title = element_blank()) 

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_OpportunityHabitatMap.png", width = 20, height = 20, units = "cm")



## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArr <- filter(CanvasAr, (Mask_G1 > 0.5)) |> mutate(Habb = "Arable") |> select(Habb, ClustPop)
CanvasGrr <- filter(CanvasGr, (Mask_G1 > 0.5)) |> mutate(Habb = "Grassland") |> select(Habb, ClustPop)
Canv <- rbind(CanvasGrr, CanvasArr) |> mutate(ClustPop= ifelse(is.na(ClustPop)==T, "#ef476f", "#118ab2"))

## make plot
ggplot() +
  ## Add in the OS basemap
  basemap_gglayer(Essex |> st_buffer(dist = 4000), map_service = "osm", map_type = "streets") +
  scale_fill_identity() +
  coord_sf(expand = FALSE) +

  ## add landscape outline
  geom_sf(data=Essex |> st_transform(crs=3857), mapping=aes(geometry=geometry), colour = "#273746", fill = NA, linewidth = 0.5) +
  coord_sf(expand = FALSE) + 
  
  ## add field polygons
  new_scale_fill() +
  geom_sf(data=Canv |> st_transform(crs=3857), mapping=aes(geometry=geometry, fill = ClustPop), colour = NA, alpha =0.75) +
  coord_sf(expand = FALSE) + 
  scale_fill_manual(name = "Strategy",   # Change legend title
                    values = c("#ef476f", "#118ab2"),
                    labels = c("Better", "Bigger/More")) +
  
  ## Add North arrow and scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.05, "in")) +
  ## set labels
  ggtitle("Essex Coast: Targetting Strategy") +
  ## set them
  theme(legend.position = "right", 
        axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
        text = element_text(color = "#2D2D2E"), 
        panel.grid = element_line(color = "#ebebe5", linewidth = 0.2),
        panel.background = element_rect(fill = "#f5f5f2", color = NA),
        axis.title = element_blank()) 

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_LawtonPrincipleMap.png", width = 20, height = 20, units = "cm")






##-------------------------------------##
## 11.8 Group 2: Plot better grading ####
##-------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasGrBetter <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==F)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrBetter, mapping=aes(geometry=geometry, fill = BetterGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Better Grading") +
  ggtitle("Essex Coast: Better (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G2_Better.png", width = 20, height = 13, units = "cm")
      


##-------------------------------------##
## 11.9 Group 2: Plot bigger grading ####
##-------------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrBig <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrBig, mapping=aes(geometry=geometry, fill = BigGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Bigger Grading") +
  ggtitle("Essex Coast: Bigger (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G2_Bigger.png", width = 20, height = 13, units = "cm")



##-----------------------------------##
## 11.10 Group 2: Plot more grading ####
##-----------------------------------##

## filter out data that is not 50% covered by mask and is not in a cluster
CanvasGrMore <- filter(CanvasGr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) + 
  ## add polygons
  geom_sf(data=CanvasGrMore, mapping=aes(geometry=geometry, fill = MoreGrade_G2), colour = NA) + 
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") + 
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "More Grading") +
  ggtitle("Essex Coast: More (Land Managers)") +
  ## set them
  theme_light() + 
  GeneralThemeing

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G2_More.png", width = 20, height = 13, units = "cm")




##------------------------------------------##
## 11.11 Group 2: Plot arable big grading ####
##------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG2 <- filter(CanvasAr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG2, mapping=aes(geometry=geometry, fill = ArableBig_G2), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Essex Coast: Arable Reversion for Bigger (Land Managers)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG2)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G2_ArableBig.png", width = 20, height = 13, units = "cm")



##-------------------------------------------##
## 11.12 Group 2: Plot arable more grading ####
##-------------------------------------------##

## Plot the better grading
## Filter out fields withing population clusters and not masked
CanvasArG2 <- filter(CanvasAr, (Mask_G2 > 0.5) & is.na(ClustPop)==T)

## make plot
ggplot() +
  ## Add coastline
  geom_sf(data=Coast, mapping = aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  ## add landscape outline
  geom_sf(data=EssOutline, mapping=aes(geometry=geometry), colour = "#273746", fill = NA) +
  ## add polygons
  geom_sf(data=CanvasArG2, mapping=aes(geometry=geometry, fill = ArableMore_G2), colour = NA) +
  ## give a viridis fill to shapes
  scale_fill_viridis_c(na.value = "lightgrey") +
  ## Set plot extent so all plots have the same extent
  PlotExt +
  ## set labels
  labs(fill = "Arable Grading") +
  ggtitle("Essex Coast: Arable Reversion for More (Land Managers)") +
  ## set them
  theme_light() +
  GeneralThemeing
rm(CanvasArG2)

## save plot as png
ggsave(filename = "CleanData/Guideline Creation/Plots/Essex_G2_ArableMore.png", width = 20, height = 13, units = "cm")
