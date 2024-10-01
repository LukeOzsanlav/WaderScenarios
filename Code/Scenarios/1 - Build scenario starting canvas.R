##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 22/04/2024
## 
## 
## Aim: Creating starting point for scenario modelling by
##      combining BWWM survey fields with UKCEH land parcels
##      into a single canvas of land parcels
## 
##------------------------------------------------------##

## Load in packages
pacman::p_load(tidyverse, data.table, sf, terra)
options(scipen=999) # turn off scientific notation



##------------------------------##
## F. Function to create canvas ##
##------------------------------##

## START OF FUNCTION
Canvas_Maker <- function(OverlapBWWM=OverlapBWWM, OverlapCEH=OverlapCEH, 
                         Min_Parcel=Min_Parcel, RatioCutoff=0.1, Outpath=Outpath){
  
  
  ##-----------------------------------------##
  #### F.1 Intersect land parcel data sets ####
  ##-----------------------------------------##
  
  ## First remove any fragement parcels
  OverlapBWWM <- OverlapBWWM |> 
                 mutate(Area = st_area(geometry), # calculate area
                        Perim = st_length(st_cast(geometry,"MULTILINESTRING")), # calculate perimeter
                        RatioArPr = Perim/Area) |>  # ratio or perimeter vs area
                 filter(as.numeric(Area) >= Min_Parcel) |> 
                 filter(as.numeric(RatioArPr) < RatioCutoff) |> 
                 select(-c(Area, Perim, RatioArPr))
  OverlapCEH <- OverlapCEH |> 
                 mutate(Area = st_area(geometry), # calculate area
                        Perim = st_length(st_cast(geometry,"MULTILINESTRING")), # calculate perimeter
                        RatioArPr = Perim/Area) |>  # ratio or perimeter vs area
                 filter(as.numeric(Area) >= Min_Parcel) |> 
                 filter(as.numeric(RatioArPr) < RatioCutoff) |> 
                 select(-c(Area, Perim, RatioArPr))
  
  ## Overlap the UKCEH land parcels with the BWWM field polygons
  ## To find and return areas where polygons from the two spatial data sets overlap
  Inter <- st_intersection(OverlapBWWM, OverlapCEH) # order does not matter as it computes shared geometries
  Inter <- st_collection_extract(Inter, "POLYGON") # remove any points/lines where shapes just touch/cross
  glimpse(Inter)
  # plot(Inter$geometry) (for debugging)
  
  
  ## Calculate the size of overlapping area for each UKCEH land parcel
  Inter2 <- Inter
  Inter2$Area <- st_area(Inter2$geometry) # calculate the area of each intersection
  Inter2 <- Inter2 |>  select(OBJECTID, Area) |> st_drop_geometry() |> rename(Area_Inter = Area) # streamline data set
  ## Group by UKCEH land parcel ID as one parcel can overlap with multiple BWWM fields, making distinct polygons
  Inter2 <- Inter2 |> group_by(OBJECTID) |> summarise(Area_Inter = sum(as.numeric(Area_Inter)))
  
  
  ## Join overlapping area onto UKCEH parcels to calculate the proportion of parcel overlapping a BWWM field
  OverlapCEH$Area <- st_area(OverlapCEH$geometry) # calculate area of each UKCEH parcel first
  OverlapCEH <- left_join(OverlapCEH, Inter2, by = "OBJECTID") # join UKCEH parcels and overlapping area
  OverlapCEH <- OverlapCEH |>  
    mutate(Area_Inter = ifelse(is.na(Area_Inter)==T, 0, Area_Inter),
           Prop_Overlap = as.numeric(Area_Inter)/as.numeric(Area), 
           Area_NoneOverlap = as.numeric(Area)-as.numeric(Area_Inter))
  hist(OverlapCEH$Prop_Overlap) 
  ## Plot to make sure proportion overlap looks sensible (for debugging)
  # ggplot() + geom_sf(data=OverlapCEH, mapping = aes(geometry=geometry, fill=Prop_Overlap), colour=NA) + theme_light()
  
  
  
  ##----------------------------------------##
  #### F.2 Filter/crop UKCEH land parcels ####
  ##----------------------------------------##
  
  ## Filter the UKCEH parcels depending on how much the parcels overlap with a BWWM field
  No_overlap <- filter(OverlapCEH, Prop_Overlap ==0)
  # plot(No_overlap$geometry)
  ## Filter out parcels that are partly overlapped by a BWWM field so they can be trimmed
  Part_overlap <- filter(OverlapCEH, Prop_Overlap > 0 & Area_NoneOverlap > Min_Parcel)
  # plot(Part_overlap$geometry)
  ## Mainly overlapped fields are removed from the UKCEH data set as there will be a BWWM field that covers pretty much the same area
  Full_overlap <- filter(OverlapCEH, Area_NoneOverlap < Min_Parcel) # 4525
  # plot(Full_overlap$geometry)
  
  
  ## Loop through each of the UKCEH parcels that overlapped a BWWM field by less then 90%
  ## For each UKCEH overlapping field extract the BWWM field it intersected with
  ## Then can extract the part of the UKCEH parcel that does not overlap with any BWWM field
  ## Streamline some data sets to speed up the loop
  OverlapBWWM <- OverlapBWWM |> select(F_LOC_ID)
  OverlapCEH <- OverlapCEH |> select(OBJECTID, X_mode)
  
  ## Start loop
  for(j in 1:nrow(Part_overlap)){
    
    ## Message of loop progress
    message(j,  " out of ", nrow(Part_overlap))
    skip_to_next <- FALSE # to skip loop if intersection fails
    
    ## Extract UKCEH land parcel using ID code
    IDj <- Part_overlap$OBJECTID[j]
  
    ## Filter our the UKCEH land parcel
    UKField <- filter(OverlapCEH, OBJECTID ==IDj)
    ## Filter out the BWWM fields that overlap with the chosen UKCEH land parcel
    BWWMField <- filter(OverlapBWWM, F_LOC_ID %in% filter(Inter, OBJECTID ==IDj)$F_LOC_ID)
  
    ## Intersect two sets of fields, if this fails then it just skips to the next iteration of the loop
    tryCatch({InterLoop <- st_difference(UKField, st_union(BWWMField))}, error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next } 
  
    ## Extract just the polygons, if needed
    if(!st_geometry_type(InterLoop)=="POLYGON"){InterLoop <- st_collection_extract(InterLoop, "POLYGON")}
  
    ## Split up all mulit-polygons into individual polygons
    InterLoop <- st_cast(InterLoop, "POLYGON")
    #ggplot() + geom_sf(data = InterLoop, fill = "blue") + geom_sf(data = BWWMField) 
    
    ## Remove small polygons
    InterLoop <- InterLoop[as.numeric(st_area(InterLoop))>Min_Parcel,]
    
    ## Join together shapes sequentially
    if(j == 1){AllNoneOver <- InterLoop}else{AllNoneOver <- bind_rows(AllNoneOver, InterLoop)}
  }
  
  # ggplot() +
  #   geom_sf(data =InterLoop, mapping = aes(geometry=geometry, fill = 1:nrow(InterLoop)), alpha = 0.6) +
  #   geom_sf(data =NoneOverlap, mapping = aes(geometry=geometry), colour = "red", alpha = 0.6) +
  #   geom_sf(data =BWWMField, mapping = aes(geometry=geometry), fill = "green", alpha = 0.6) 
  
  ## plot the output (for deugging)
  # plot(AllNoneOver$geometry)
  
  ## Calculate the area of each portion of a UKCEH land parcel that was extracted
  AllNoneOver <- AllNoneOver |> 
    mutate(Cutarea = st_area(geometry), # calculate area
           Cutperim = st_length(st_cast(geometry,"MULTILINESTRING")), # calculate perimeter
           RatioArPr = Cutperim/Cutarea, # ratio or perimeter vs area
           UKCEH_crop = 1) |> 
    filter(as.numeric(RatioArPr) < RatioCutoff) # remove object with high ration, generally fragments of field edges
  paste0("Min parcel size UKCEH ", min(AllNoneOver$Cutarea))
  
  
  
  ##-------------------------------------------##
  #### F.3 Create combined landscape parcels ####
  ##-------------------------------------------##
  
  ## Now to join all the three data sets to make a full canvas
  No_overlap <- No_overlap |> select(OBJECTID, X_mode)
  OverlapBWWM <- OverlapBWWM |> select(F_LOC_ID)
  AllNoneOver <- AllNoneOver |> select(OBJECTID, X_mode, UKCEH_crop) ##### **CHECK HERE** #### may end up with some OBJECTIDs being repeated, could fix this with a group_by and st_uniotn after the loop
  AllParcels <- bind_rows(No_overlap, OverlapBWWM) |> bind_rows(AllNoneOver)
  
  ## Add indicator column to show whether polygon is an original BWWM field or not
  AllParcels <- AllParcels |> mutate(BWWM = ifelse(is.na(F_LOC_ID)==T, 0, 1))
  write_sf(AllParcels, Outpath)
  
} # END OF FUNCTION




##-----------------------------------##
#### 1.0 Read in spatial data sets ####
##-----------------------------------##

## Read in priority landscape polygons
Pr <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |>  
  select(Att3Value) |>  rename(Landscape = Att3Value)
Pr <- filter(Pr, Landscape %in% c("Broads", "Somerset Levels and Moors", "North Kent", "Essex")) 
plot(Pr$geometry)
## Buffer priority landscapes
Pr <- st_buffer(Pr, dist = 3000)
plot(Pr$geometry)

## Read in BWWM field polygons
BWWM <- st_read("RawData/BWWM field shapefile/BWWM_fields_27Sep2023.shp")
BWWM <- BWWM |> select(F_LOC_ID)

## Read in the UKCEH land parcels data set from 2021
UKC <- st_read("RawData/LandCover/UKCEH_LCvector_2021.shp") |> select(OBJECTID, X_mode)
glimpse(UKC)




##--------------------------------------------##
#### 2.0 Crop out fields for each landscape ####
##--------------------------------------------##

## Work out which fields fall within any of the buffered priority landscapes
OverlapBWWM <- st_join(BWWM, Pr)

## Work out which fields fall within any of the buffered priority landscapes
OverlapCEH <- st_join(UKC, Pr)
rm(UKC, BWWM); gc()


## Filter out Broads fields for BWWM and UKCEH
BroadsBWWM <- filter(OverlapBWWM , Landscape == "Broads")
plot(BroadsBWWM$geometry)
BroadsCEH <- filter(OverlapCEH , Landscape == "Broads")
plot(BroadsCEH$geometry)

## Filter out Broads fields for BWWM and UKCEH
NKentBWWM <- filter(OverlapBWWM , Landscape == "North Kent")
plot(NKentBWWM$geometry)
NKentCEH <- filter(OverlapCEH , Landscape == "North Kent")
plot(NKentCEH$geometry)

## Filter out Broads fields for BWWM and UKCEH
EssexBWWM <- filter(OverlapBWWM , Landscape == "Essex")
plot(EssexBWWM$geometry)
EssexCEH <- filter(OverlapCEH , Landscape == "Essex")
plot(EssexCEH$geometry)

## Filter out Broads fields for BWWM and UKCEH
SomBWWM <- filter(OverlapBWWM , Landscape == "Somerset Levels and Moors")
plot(SomBWWM$geometry)
SomCEH <- filter(OverlapCEH , Landscape == "Somerset Levels and Moors")
plot(SomCEH$geometry)
rm(OverlapBWWM, OverlapCEH); gc() # free up space




##------------------------------------------##
#### 3.0 Create canvas for each landscape ####
##------------------------------------------##

## Run function for each landscape
Canvas_Maker(OverlapBWWM=BroadsBWWM, OverlapCEH=BroadsCEH, RatioCutoff = 0.06,
             Min_Parcel=2500, Outpath="CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp")

Canvas_Maker(OverlapBWWM=NKentBWWM, OverlapCEH=NKentCEH, RatioCutoff = 0.06,
             Min_Parcel=2500, Outpath="CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp")

Canvas_Maker(OverlapBWWM=SomBWWM, OverlapCEH=SomCEH, RatioCutoff = 0.06,
             Min_Parcel=2500, Outpath="CleanData/Scenarios/1-Starting Canvas/Som_Canvas.shp")

Canvas_Maker(OverlapBWWM=EssexBWWM, OverlapCEH=EssexCEH, RatioCutoff = 0.06,
             Min_Parcel=5000, Outpath="CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp")


