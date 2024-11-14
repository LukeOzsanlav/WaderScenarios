##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 23/04/2024
## 
## Aim: Annotate scenario canvas with important spatial data sets
## 
##------------------------------------------------------##

## Need to go pack and sort out a few columns
## Mainly their is some mismatch between Reserve/AES classification and arable/grass opportunity column

## Load in packages & helper functions
pacman::p_load(tidyverse, data.table, sf, terra, exactextractr, ggpattern)
options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R")




##---------------------------------##
#### 0.1 Define Worker Functions ####
##---------------------------------##

## Determine if small polygons fall within a larger polygon 
## Can manipulate the percentage of small polygon has to be covered by the larger polygon to be considered "overlapping"

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
  Overlap_sub <- Fields |> 
                 st_intersection(OverlapShapes) |> 
                 mutate(OverlapArea = st_area(geometry)) |> 
                 select(ParcRef, OverlapArea) |> 
                 st_drop_geometry()
  
  ## Remove any duplicated in this data set, can't work out why duplicates happen, maybe overlapping shapes?
  Dups <- Overlap_sub |> select(ParcRef, OverlapArea) |> duplicated()
  Overlap_sub <- Overlap_sub[Dups==FALSE,]
  Overlap_sub <- Overlap_sub |> group_by(ParcRef) |> summarise(OverlapArea = sum(OverlapArea))

  ## Now join the overlap size onto the main fields data &
  ## Work out the proportion of each field in a reserve
  Fields <- left_join(Fields, Overlap_sub, by = "ParcRef") |> 
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


## Function that takes land parcel polygons and a raster  
## Then calculates the average value of the raster within in a buffer around the center of each field

## Survey = land parcels data as a shape file of each field
## Rast = Raster of variable of interest, should only contain 0 and 1
## Dist = Buffer distances to trial
## Col = Start of the column name to use for the proportion cover column (distance used as suffix)
Buff_RastVal <- function(Parcels = Parcels, Rast = Rast, Dist = Dist, Col = Col){
  
  ## crop the raster data to the extent of the survey data (helps with speed)
  Rastcrop <- crop(Rast, vect(st_buffer(st_as_sf(Parcels), dist = 2050)))
  plot(Rastcrop)
  
  ## Get the centers for each field and then add the buffer
  CentBuffs <- st_centroid(st_as_sf(Parcels)) |> st_buffer(dist = Dist)
     
  ## Extract values from raster to fields (weighted average)
  
  ExtractVals <- exact_extract(x = Rastcrop,
                               y = st_as_sf(Parcels) |> st_centroid() |> st_buffer(dist = Dist),
                               fun = "mean")
  extrac
  ## Add the buffer values to the main wader data set
  Parcels <- mutate(Parcels, VAL = ExtractVals[,2]) 
  colnames(Parcels)[colnames(Parcels) == "VAL"] <- paste0(Col, Dist)
  
  return(Parcels)
} # end of function


## **Code for checking function, can delete when all finalsies
# Canv = EssexCanvas
# AllClust = EsWaderClust
# GrassOpp = EssexGrassOpp
# ArableOpp = EssexArableOpp
#  PeatSoilsClass = NULL
# WatRast = EsW
# IntertidalDist=TRUE
# Label = EsLabel
# guidepath = Espaths
# outpath = EsOut
# RunGrades=FALSE

# Canv = KentCanvas
# AllClust = NKWaderClust
# GrassOpp = NKentGrassOpp
# ArableOpp = NKentArableOpp
# PeatSoilsClass = NULL
# WatRast = NKW
# IntertidalDist=TRUE
# Label = NKLabel
# guidepath = NKpaths
# outpath = NKOut
# RunGrades=TRUE

# Canv = Canvas
# AllClust = WaderClust
# GrassOpp = SomGrassOpp
# ArableOpp = SomArableOpp
# PeatSoilsClass = SomPeatClass
# WatRast = SomW
# Label = SomLabel
# guidepath = Sompaths
# outpath = SomOut
# RunGrades=TRUE


## START of AnnotateCanv ##
## This is the function to annotate the canvas for any given priority landscape
## These canvas form the starting point of the scenario modelling
AnnotateCanv <- function(Canv, ## Blank un-labelled canvas for each priority landscape created in `1 - Build scenario starting canvas.R`
                         AllClust, ## All wader clusters for the priority landscapes, created in `2 - Define wader sites.R`
                         GrassOpp, ## lowland grassland opportunity area for a priority landscape, created in `3 - Define conservation action areas.R`
                         ArableOpp, ## arable opportunity area for a priority landscape, created in `3 - Define conservation action areas.R`
                         PeatSoilsClass=NULL, ## If the landscape needs to be labelled with peat soil presence then provide a vector of peaty soil types in the NATMAP vector
                         IntertidalDist=FALSE, ## If TRUE then calculate the distance to the nearest patch of inter tidal habitat for each field
                         WatRast, ## Raster of standing water coverage, this is created in a different project- see `Standing Water Classification`
                         Label, ## Label for the priority landscape
                         guidepath, ## Part of file path to read in the stakeholder guidelines e.g. "Somerset/Som"
                         outpath, # This is the file path for the final .csv and .shp to be read out, start of filepath is "CleanData/Scenarios/4-AnnotateCanvas/"
                         RunGrades=FALSE # TRUE/FALSE if want to calculate stakeholder gradings for priority landscape or not
                         ){ 
  
  ##-----------------------------##
  #### F1.0 Add Parcel Details ####
  ##-----------------------------##
  
  ## assign an area and a unique reference ID to each parcel
  Canv <- Canv |> 
    mutate(ParcArea = st_area(geometry),
           ParcRef = 1:nrow(Canv))
  
  
  
  
  ##----------------------------##
  #### F1.1 Label Wader Sites ####
  ##----------------------------##
  
  ## shape file of wader sites for region of interest
  plot(AllClust$geometry)
  
  ## Intersect canvas and wader site polygons and work out proportion of land parcel in the cluster
  Inter <- st_intersection(Canv, AllClust) |> 
           mutate(InterArea = st_area(geometry),
                  OverlapProp = as.numeric(InterArea)/as.numeric(ParcArea)) 
  ## pot the intersection
  # ggplot(data = Inter) + geom_sf(aes(fill = as.numeric(OverlapProp)))
  
  ## If a field is only 50% covered by a site them remove it from the site
  Inter <- filter(Inter, OverlapProp >= 0.5) |> 
           select(ParcRef, ClustGroup) |> 
           st_drop_geometry()
  
  ## Join the Wader Cluster ID code onto the main canvas
  Canv2 <- left_join(Canv, Inter, by = "ParcRef")
  stopifnot(nrow(Canv2)==nrow(Canv))
  Canv <- Canv2; rm(Canv2)
  
  ## visualize the cluster ID numbers
  # ggplot() + 
  #    geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = as.character(ClustGroup))) +
  #    theme_light()
  
  ## free up space
  rm(Inter); gc()
  
  
  
  
  ##--------------------------------##
  #### F1.2 Label Wader Abundance ####
  ##--------------------------------##
  
  ## Read in filtered breeding pairs estimates
  ## And calculate total wader abundance for each field
  Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
  Waders <- Waders |> mutate(Tot_abund = rowSums(dplyr::select(Waders, est_pairsL, est_pairsR, est_pairsS), na.rm = T)) |> 
                      select(F_LOC_ID, Tot_abund) |> group_by(F_LOC_ID) |>
                      summarise(Tot_abund = sum(Tot_abund, na.rm = T))
  
  
  ## Join the Wader Cluster ID code onto the main canvas
  Canv2 <- left_join(Canv, Waders, by = "F_LOC_ID")
  stopifnot(nrow(Canv2)==nrow(Canv))
  Canv <- Canv2; rm(Canv2)
  
  ## free up space
  rm(Waders); gc()
  
  
  
  
  ##--------------------------------------------##
  #### F1.3 Calculate Distance to Wader Sites ####
  ##--------------------------------------------##
  
  ## Now add on the distance from a field's centroid to the nearest wader sites
  ## This will be 0m for fields in wader sites
  ## First create index of the nearest wader site for each field centroid
  NFeat <- st_nearest_feature(st_centroid(Canv), AllClust)
  
  ## Now use this index to calculate the distance between the field centorid and its corresponding nearest wader site
  DistMatrix <- st_distance(st_centroid(Canv), AllClust[c(NFeat),], by_element = T)
  
  ## Assign the distances to a column
  Canv$ClustDist <- as.numeric(DistMatrix)
  
  ## plot to visualise the distance
  # ggplot() + 
  #    geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = (ClustDist)), colour = NA) +
  #    theme_light()
  
  ## free up space
  rm(NFeat, DistMatrix); gc()
  
  
  
  
  ##----------------------------------##
  #### F1.4 Label Opportunity Areas ####
  ##----------------------------------##
  
  ## Convert the water raster into a binary raster
  WatRast2 <- WatRast
  values(WatRast2) <- ifelse(values(WatRast2) =="NaN", 0, 1)

  ## Extract values from opportunity raster to field parcel canvas
  GrVals <- exact_extract(x = GrassOpp, y = Canv, fun = "mean")
  ArVals <- exact_extract(x = ArableOpp, y = Canv, fun = "mean")
  WetGrassVals <- exact_extract(x = WatRast2, y = Canv, fun = "mean")

  ## Label whether fields are part of grassland or arable opportunity areas or not
  Canv <- Canv |> mutate(GrassOpp = GrVals,
                         ArableOpp = ArVals,
                         WiderWetGrass = WetGrassVals,
                         GrassOpp = ifelse(GrassOpp == "NaN", 0, GrassOpp), 
                         ArableOpp = ifelse(ArableOpp == "NaN", 0, ArableOpp),
                         WiderWetGrass = ifelse(WiderWetGrass == "NaN", 0, WiderWetGrass),
                         GrassOpp = ifelse(GrassOpp >=0.5, 1, 0), 
                         ArableOpp = ifelse(ArableOpp >=0.5, 1, 0), 
                         WiderWetGrass = ifelse(WiderWetGrass >=0.5, 1, 0),
                         Opp = ifelse(GrassOpp ==1, "Grass",
                                      ifelse(ArableOpp ==1 , "Arable", "None")))
  
  # ## visualize the extractions
  # ggplot() + 
  #    geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = (GrassOpp))) +
  #    theme_light()
  # 
  # ## visualize the extractions
  # ggplot() + 
  #    geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = (ArableOpp))) +
  #    theme_light()
  # 
  # ## visualize the extractions
  # ggplot() + 
  #    geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = (Opp))) +
  #    theme_light()
  
  ## free up space
  rm(GrassOpp, ArableOpp, GrVals, ArVals, WatRast2); gc()
  
  
  
  
  ##------------------------------##
  #### F1.5 Label Reserve Areas ####
  ##------------------------------##
  
  ## Read in the various reserve outlines
  ## RSPB reserves
  RSPB <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp")
  ## rename reserve name column and remove none-wetgrassland reserves, done visually via QGIS
  RSPB <- RSPB |> rename(ReserveName=Att4Value) |> filter(!ReserveName == "BLEAN WOODS")
  RSPB <- RSPB |> group_by(ReserveName) |>  summarize(geometry = st_union(geometry)) |> mutate(RSPBArea = as.numeric(st_area(geometry)))

  ## National Nature Reserves
  NNR <- st_read("RawData/Other Reserves/National_Nature_Reserves_EnglandPolygon.shp")
  ## rename reserve name column and remove none-wet grassland reserves, done visually via QGIS
  NNR <- NNR |> filter(!nnr_name %in% c("Leigh", "Dengie", "Colne Estuary", "Blackwater Estuary", "Ant Broads and Marshes", "Winterton Dunes", "Rodney Stoke", "Ebbor Gorge"))
  NNR <- NNR |> group_by(nnr_name) |>  summarize(geometry = st_union(geometry)) |> mutate(NNRArea = as.numeric(st_area(geometry)))
  
  ## Overlap the RSPB and NNR polygons
  Inter <- st_intersection(RSPB, NNR) |> mutate(OverlapArea = as.numeric(st_area(geometry)),
                                                OverlapProp = OverlapArea/NNRArea)
  
  ## Remove the NNR polygons if they are overlapped by a RSPB reserve
  InterMost <- filter(Inter, OverlapProp >= 0.10)
  NNR <- NNR |> filter(!nnr_name %in% InterMost$nnr_name)
  
  
  
  ## RSPB reserve overlap
  RSPBagg <- st_as_sf(terra::aggregate(vect(RSPB)))
  Canv <- Field_Overlap_IDer(Fields = Canv, OverlapShapes = RSPBagg, 
                             Thresh = 0.5, ColName = "RSPB")
  ## NNR's
  NNRagg <- st_as_sf(terra::aggregate(vect(NNR)))
  Canv <- Field_Overlap_IDer(Fields = Canv, OverlapShapes = NNRagg, 
                             Thresh = 0.5, ColName = "NNR")
  
  ## Add on the reserve names
  Canv <- st_join(Canv, RSPB["ReserveName"], join = st_intersects)
  Canv <- st_join(Canv, NNR["nnr_name"], join = st_intersects)
  
  ## remove any duplicates, these come about if a field is overlapped by multiple reserves
  Canv <- Canv[duplicated(Canv$ParcRef)==F,]

  
  
  ## Add extra column to indicate if a field is of a reserve of any type
  Canv <- Canv |> mutate(Reserve = ifelse(RSPB == "Y" | NNR == "Y", "Y", "N"),
                         ReserveGroup = ifelse(NNR == "Y", nnr_name, ReserveName),
                         ReserveGroup = ifelse(Reserve == "N", NA, ReserveGroup)) |> 
                  select(-c(nnr_name, ReserveName))
  
  ## Add on the size of the reserve group
  RSPB <- RSPB |> rename(ReserveGroup=ReserveName, GroupArea=RSPBArea) |> st_drop_geometry()
  NNR <- NNR |> rename(ReserveGroup=nnr_name, GroupArea=NNRArea) |> st_drop_geometry()
  ReserveSizes <- rbind(RSPB, NNR)
  Canv <- left_join(Canv, ReserveSizes, by = "ReserveGroup")
  
  ## plots to check
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = RSPB), colour = NA) + theme_minimal()
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = NNR), colour = NA) + theme_minimal()
  
  ## free up space
  rm(RSPBagg, NNRagg, RSPB, NNR); gc()
  
  
  
  
  ##---------------------------------##
  #### F1.6 Label Wader AES fields ####
  ##---------------------------------##
  
  ## Read in the Different Stewardship Schemes
  CSS_Ops <- st_read("RawData/Stewardship/Countryside_Stewardship_Scheme_Management_Options_(England)___Natural_England.shp")
  ESSWest <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_WestEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
  ESSEast <- st_read("RawData/Stewardship/Environmental_Stewardship_Scheme_Options_EastEngland/Environmental_Stewardship_Scheme_Options_EnglandPoint.shp")
  ESS <- rbind(ESSWest, ESSEast)
  rm(ESSWest, ESSEast)

  ## Crop the extents of the stewardship options to speed up processing speed
  CSS_Ops <- st_crop(CSS_Ops, ext(Canv))
  ESS <- st_crop(ESS, ext(Canv))

  
  ## Filter out wader focused options
  ESS_OpsWa <- filter(ESS, optcode %in% c("HK9", "HK10", "HK11", "HK12", "HK13", "HK14", "HK19"))

  ## Since these CSS options are points, determine if a points falls within one of the survey fields
  ## And label each field whether a wader specific CSS point fell within it
  ESS_InterWa <- as.data.frame(st_covers(Canv, ESS_OpsWa))
  Canv$ESS_Wader <- "N"
  Canv$ESS_Wader[ESS_InterWa$row.id] <- "Y"
  
  
  ## Filter out all the land parcels that have an option for breeding or wintering waders
  ESSWaderFields <- filter(ESS, parcref %in% ESS_OpsWa$parcref)
  
  ## Read in the costings spreadsheet, that gives the options codes that we are going to work out the payment
  ESS_Costs <- read.csv("RawData/AES Costings/ESS_Cost_Sheet.csv") |> 
               filter(!CSS_Code == "SP1") # remove this one as is only for small fields
  
  ## Join the AES costs spreadsheet onto the AES option points
  ## This means that each ESS option has a equivalent CSS option
  ESSWaderFields <- left_join(ESSWaderFields, ESS_Costs, by = c("optcode" = "ESS_Code"))

  ## Create a list column of all the ESS options of interest for each parcel
  Field_ESSOpts <- ESSWaderFields |> 
                   filter(optcode %in% ESS_Costs$ESS_Code) |> 
                   group_by(parcref) |> 
                   summarise(AES_Options = list(unique(CSS_Code)))
  
  ## No assign these CSS options to the correct field in my canvas
  ESS_InterOpts <- as.data.frame(st_covers(Canv, Field_ESSOpts))
  Canv$AES_options <- list(NA)
  Canv$AES_options[ESS_InterOpts$row.id] <- Field_ESSOpts$AES_Options[ESS_InterOpts$col.id]


  ## plots to check
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Wader), colour = NA) + theme_minimal()
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = ESS_Gen), colour = NA) + theme_minimal()
  
  
  
  
  ## Add column to CSS data to indicate if agreements includes options that benefit waders or general options that could benefit waders
  ## See Appenidx 1 of Rob's BWWM report for what is classified as "wader focused" and "general" options
  ## Filter out wader focused options
  CSS_OpsWa <- filter(CSS_Ops, opt_code %in% c("GS9", "GS10", "GS11", "GS12"))

  ## Since these CSS options are points, determine if a points falls within one of the survey fields
  ## And label each field whether a wader specific CSS point fell within it
  CSS_InterWa <- as.data.frame(st_covers(Canv, CSS_OpsWa))
  Canv$CSS_Wader <- "N"
  Canv$CSS_Wader[CSS_InterWa$row.id] <- "Y"
  
  
  ## Filter out all the land parcels that have an option for breeding or wintering waders
  CSSWaderFields <- filter(CSS_Ops, parcref %in% CSS_OpsWa$parcref)
  
  ## Read in the costings spreadsheet, that gives the options codes that we are going to work out the payment
  CSS_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv") |> 
               filter(!CSS_Code == "SP1") # remove this one as is only for small fields
  
  ## Create a list column of all the CSS options of interest for each parcel
  Field_CSSOpts <- CSSWaderFields |> 
                   filter(opt_code %in% CSS_Costs$CSS_Code) |> 
                   group_by(parcref) |> 
                   summarise(AES_Options = list(opt_code))
  
  ## No assign these CSS options to the correct field in my canvas
  ## Save the options in the same column as the ESS options
  ## I have already converted the ESS option codes to CSS codes
  CSS_InterOpts <- as.data.frame(st_covers(Canv, Field_CSSOpts))
  Canv$AES_options[CSS_InterOpts$row.id] <- Field_CSSOpts$AES_Options[CSS_InterOpts$col.id]
    
  
  ## plots to check
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = CSS_Wader), colour = NA) + 
  #            theme_minimal()
  # ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = CSS_Gen), colour = NA) +
  #            theme_minimal()
  
  ## free up space
  rm(ESS, CSS_InterWa, CSS_Ops, CSS_OpsWa, ESS_OpsWa, ESS_InterOpts, CSS_InterOpts, CSSWaderFields,
     ESS_InterWa, ESSWaderFields, CSS_Costs, ESS_Costs, Field_ESSOpts, Field_CSSOpts); gc()
  
  
  
  
  ##-------------------------------------------------##
  #### F1.7 Categorize Land Management/Opportunity ####
  ##-------------------------------------------------##
  
  ## Firstly, I need to correct some Land management or opportunity classification of some fields
  ## This happens perhaps when fields get wrongly classified as arable from UKCEH
  ## I am going to presume that any reserve only parcels that are arable opportunity should actually be grass opportunity
  ## But based off looking at maps I am going to say that any fields that are CSS or ESS and arable opportunity are not in fact part of the AES agreement
  ## Looking at the Broads this seems to be because large parts of the ESS outline are actually arable and not lowland grassland
  if(!Label == "Essex"){
    
      Canv <- Canv |> mutate(GrassOpp = ifelse(ArableOpp == 1 & Reserve == "Y", 1, GrassOpp),
                         ArableOpp = ifelse(ArableOpp == 1 & Reserve == "Y", 0, ArableOpp),
                         CSS_Wader = ifelse(ArableOpp == 1 & CSS_Wader == "Y", "N", CSS_Wader),
                         ESS_Wader = ifelse(ArableOpp == 1 & ESS_Wader == "Y", "N", ESS_Wader)) 
  }
  
  ## Different function for Essex
  ## I visually looked and if a reserve is classified as arable then it is normally a mistake in Kent and Broads
  ## But in Essex most of the reserve fields that are Arable are actually arable or none-wet grassland so turn these from reserve
  ## Back to just arable opportunity and Reserve as "N"
  if(Label == "Essex"){
    
      Canv <- Canv |> mutate(Reserve = ifelse(ArableOpp == 1 & Reserve == "Y", "N", Reserve),
                             CSS_Wader = ifelse(ArableOpp == 1 & CSS_Wader == "Y", "N", CSS_Wader),
                             ESS_Wader = ifelse(ArableOpp == 1 & ESS_Wader == "Y", "N", ESS_Wader)) 
  }
  
  ## Ensure the Opp column is correct, i.e. if it is wader management then it 
  Canv <- Canv |> mutate(Opp = ifelse(Reserve == "Y" | ESS_Wader == "Y" | CSS_Wader == "Y", "Grass", Opp))

  
  ## Add a category variable to look at where different land categories lie within the landscape
  Canv <- Canv |> mutate(Category= case_when(
                                Reserve == "Y" & GrassOpp == 1 ~ "Reserve",
                                (CSS_Wader == "Y" | ESS_Wader  == "Y") & Reserve == "N" & GrassOpp == 1  ~ "AES Only",
                                GrassOpp == 1 & CSS_Wader == "N" & ESS_Wader  == "N" & Reserve == "N" ~ "Grass Opp",
                                ArableOpp == 1 & CSS_Wader == "N" & ESS_Wader  == "N" & Reserve == "N" ~ "Arable Opp",
                                GrassOpp == 0 | ArableOpp == 0 ~ "NoOpp",
                                .default = "NoOpp"))
  
  ## Plot the land categorization classes
  ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
    scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "lightgrey", "#4C77DD")) +
             theme_light()
  
  ## save the plot
  ggsave(paste0("CleanData/Scenarios/4-AnnotateCanvas/Plots/", Label, "_Categories.png"), width = 20, height = 20, units = "cm")
  
  
  ## Plot the land categorization classes
  ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
             scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "lightgrey", "#4C77DD")) +
             geom_sf(data = AllClust, mapping = aes(geometry = geometry, colour = ""), fill = NA) +
             scale_colour_manual(values = c("black")) +
             labs(colour = "Wader Sites") +
             theme_light()
  
  ## save the plot
  ggsave(paste0("CleanData/Scenarios/4-AnnotateCanvas/Plots/", Label, "_Categories&WaderSites.png"), width = 20, height = 20, units = "cm")
  
  
  
  ##------------------------------------------##
  #### F1.8 Label with Stakeholder Gradings ####
  ##------------------------------------------##
  
  if(RunGrades==TRUE){
    
    ## Read in all the gradings and change the layer names to something meaningful
    Arable_G1 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_ArableBig_G1.tif")) 
    names(Arable_G1) <- "Arable_G1"
    Arable_G2 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_ArableBig_G2.tif"))
    names(Arable_G2) <- "Arable_G2"
    Better_G1 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Better_G1.tif"))
    names(Better_G1) <- "Better_G1"
    Better_G2 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Better_G2.tif"))
    names(Better_G2) <- "Better_G2"
    Bigger_G1 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Bigger_G1.tif"))
    names(Bigger_G1) <- "Bigger_G1"
    Bigger_G2 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Bigger_G2.tif"))
    names(Bigger_G2) <- "Bigger_G2"
    More_G1 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_More_G1.tif"))
    names(More_G1) <- "More_G1"
    More_G2 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_More_G2.tif"))
    names(More_G2) <- "More_G2"
    Mask_G1 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_MasksAll_G1.tif"))
    names(Mask_G1) <- "Mask_G1"
    Mask_G2 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_MasksAll_G2.tif"))
    names(Mask_G2) <- "Mask_G2"
    
    ## Alter mask layers to that mask ==1 and none maksk = 0
    values(Mask_G1) <- ifelse(values(Mask_G1)== "NaN", 1, 0)
    values(Mask_G2) <- ifelse(values(Mask_G2)== "NaN", 1, 0)
    
    ## For Somerset need to add on a third group that was not done for the other three regions
    if(Label == "Somerset Levels and Moors"){
      
    Arable_G3 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_ArableBig_G3.tif")) 
    names(Arable_G3) <- "Arable_G3"
    Better_G3 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Better_G3.tif"))
    names(Better_G3) <- "Better_G3"
    Bigger_G3 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_Bigger_G3.tif"))
    names(Bigger_G3) <- "Bigger_G3"
    More_G3 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_More_G3.tif"))
    names(More_G3) <- "More_G3"
    Mask_G3 <- rast(paste0("CleanData/Guideline Creation/", guidepath,"_MasksAll_G3.tif"))
    names(Mask_G3) <- "Mask_G3"
    
    ## Alter mask layers to that mask ==1 and none maksk = 0
    values(Mask_G3) <- ifelse(values(Mask_G3)== "NaN", 1, 0)
      
    }
    

    ## Stack all the raster layers so that we can to the extraction in one go
    if(Label == "Somerset Levels and Moors"){
      Grades <- c(Arable_G1, Arable_G2, Better_G1, Better_G2, Bigger_G1, Bigger_G2, More_G1, More_G2, Mask_G1, Mask_G2,
                  Arable_G3, Better_G3, Bigger_G3, More_G3, Mask_G3)
      rm(Arable_G1, Arable_G2, Better_G1, Bigger_G1, Bigger_G2, More_G1, More_G2, Better_G2,
         Mask_G1, Mask_G2, Arable_G3, Better_G3, Bigger_G3, More_G3, Mask_G3);gc()
    }
    if(!Label == "Somerset Levels and Moors"){
      Grades <- c(Arable_G1, Arable_G2, Better_G1, Better_G2, Bigger_G1, Bigger_G2, More_G1, More_G2, Mask_G1, Mask_G2)
      rm(Arable_G1, Arable_G2, Better_G1, Bigger_G1, Bigger_G2, More_G1, More_G2, Better_G2, Mask_G1, Mask_G2);gc()
    }
    
    ## Extract all the gradings to all the land parcels with the canvas
    Vals <- exact_extract(x = Grades, y = Canv, fun = "mean")
    names(Vals) <- substring(names(Vals), 6) # remove mean from the start of all column names
    rm(Grades); gc()
    
    ## Now add on the extracted values back to the canvas
    Canv <- cbind(Canv, Vals)
    
    
    ## One extra grading needs to be added on to the better rule (cluster size)
    ## Read in data set that gives population size for each cluster
    PopSize <- read_csv(paste0("CleanData/Guideline Creation/", guidepath, "_ClusterPopSize.csv"))
    PopSize <- PopSize |> select(ClustPop, F_LOC_ID) |> mutate(ClustPop = Inv_scale_vals(ClustPop)) # extract population cluster sizes
    Canv <- left_join(Canv, PopSize, by = "F_LOC_ID") # add these pop sizes to the fields (if they are in a cluster)
    
    ## Do the final addition to work out grade
    # Canv$ClustPop <- ifelse(is.na(Canv$ClustPop)==T, 0, Canv$ClustPop)
    # Canv$Better_G1 <- Canv$Better_G1 + Canv$ClustPop
    # Canv$Better_G2 <- Canv$Better_G2 + Canv$ClustPop
    Canv <- Canv |> dplyr::mutate(ClustPop = ifelse(is.na(ClustPop)==T, 0, ClustPop),
                                  Better_G1 = Better_G1 + ClustPop,
                                  Better_G2 = Better_G2 + ClustPop)

    
    if(Label == "Somerset Levels and Moors"){Canv <- Canv |>  dplyr::mutate(Better_G3 = Better_G3 + ClustPop) |> select(-ClustPop)}
    if(!Label == "Somerset Levels and Moors"){Canv <- Canv |> select(-ClustPop)}
    
    ## Plot the land categorization classes
    ArableOpp <- filter(Canv, ArableOpp == 1)
    ggplot() + geom_sf(data = ArableOpp, mapping = aes(geometry = geometry, fill = Arable_G2), colour = NA) +
               theme_minimal()
    
    GrOpp <- filter(Canv, GrassOpp == 1)
    ggplot() + geom_sf(data = GrOpp, mapping = aes(geometry = geometry, fill = Bigger_G2), colour = NA) +
               theme_minimal()
    
    ## free up space
    rm(Vals, PopSize, ArableOpp, GrOpp); gc()
    
  }
  
  
  
  
  ##------------------------------------------##
  #### F1.9 Annotate BWWM Survey Variables  ####
  ##------------------------------------------##
  
  ## Firstly going to split up the data into Grassland opportunity and other areas
  ## I will only label all of the grassland opportunity fields with information
  CanvGr <- filter(Canv, Category %in% c("Grass Opp", "Reserve", "AES Only"))
  CanvOther <- filter(Canv, !Category %in% c("Grass Opp", "Reserve", "AES Only"))
  ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
             theme_minimal()
  
  ## Read in the field survey data and select the variables that were included in the model
  FieldChar <- read_csv("CleanData/Wader Abundance/1-CleaningScript/Field_Characteristics_Clean.csv")
  
  ## Remove duplicated in field characteristics (this occurred due to repeat surveys in Somerset)
  ## Also just extract columns needed for modelling
  Dups <- duplicated(FieldChar$F_LOC_ID)
  FieldChar <- FieldChar[Dups==F,] |> select(F_LOC_ID, GRASSLAND_TYPE, TALL_BOUNDARY_PERCENT, STOCK, VEG_STRUCTURE, RUSH_PERCENT)
  
  ## Join the survey data onto the Canvas and then do a visualization
  CanvGr2 <- left_join(CanvGr, FieldChar, by = "F_LOC_ID")
  stopifnot(nrow(CanvGr2)==nrow(CanvGr)) # check join
  CanvGr <- CanvGr2; rm(CanvGr2)
  
  ## Visualize survey data
  #CanvGrTest <- filter(CanvGr, is.na(RUSH_PERCENT)==T)
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = RUSH_PERCENT), colour = NA) +
  #            theme_minimal()
  rm(FieldChar)
  
  
  ## This function takes a polygon data set of the landscape canvas
  ## It also takes a column name, in this case one from the BWWM survey
  ## For Any rows with missing habitat survey data it fills these in with data from other fields
  ## If there is a field within 1000m with habitat values then it uses the values from the closest field
  ## If no fields with habitat values are nearby then sample from all other fields that fall in the same category (reserve/AES only/other)
  FillBlankHab <- function(Data=Data, ColName = ColName){
    
    
    ## 
    Can2 <- Data 
    
    ## Seperate the reserve parcels without habitat values to the ones with habitat values
    ## Use Reserve== "Y" top capture any reserve parcels that might be just outside the priority landscape
    noHab <- which(is.na(Can2[,colnames(Can2)==ColName])==T & Can2$Category == "Reserve")
    Hab <- which(is.na(Can2[,colnames(Can2)==ColName])==F & Can2$Reserve== "Y")
    
    # Find the nearest polygon in Can2[Hab,] for each polygon in Can2[noHab,]
    nearest_indices <- st_nearest_feature(Can2[noHab,], Can2[Hab,])
    
    # Calculate the distances to the nearest polygons
    distances <- st_distance(Can2[noHab,], Can2[Hab,][nearest_indices, ], by_element = TRUE)
    
    ## If the nearest parcel is more than 1000m away then make it an NA
    distacneedit <- ifelse(as.numeric(distances)>1000, NA, nearest_indices)
    
    ## Now assign the new habitat values to the correct rows that are not too far from their nearest neighbor
    Data[noHab, colnames(Data)%in% ColName] <- Can2[Hab,][distacneedit, colnames(Can2)%in% ColName] |> st_drop_geometry()
    
    ## Remove geometry as can mess up sub-setting by columns
    Can2 <- Data |> st_drop_geometry()
    
    ## For Reserve parcels missing a habitat value sample a habitat value from the same category
    InterS <- sample(x= na.omit((Can2[Can2$Reserve== "Y",])[, colnames(Can2)==ColName]), 
                     size = nrow(Can2[is.na(Can2[,colnames(Can2)==ColName])==T & Can2$Reserve== "Y",]), replace = T)
    ## Assign the sampled value to the rows with the missing habitat values
    Can2[is.na(Can2[,colnames(Can2)==ColName])==T & Can2$Reserve== "Y" , colnames(Can2)==ColName] <- InterS
    
    ## Assign habitat values to missing rows
    Data[colnames(Data)==ColName] <- Can2[colnames(Can2)==ColName]
    
    
    ##
    
        
    ## 
    Can2 <- Data 
    
    ## Seperate the AES only parcels without habitat values to the ones with habitat values
    ## Capture any AES only parcels that might be just outside the priority landscape for Hab set
    noHab <- which(is.na(Can2[,colnames(Can2)==ColName])==T & Can2$Category == "AES Only")
    Hab <- which(is.na(CanvGr[,colnames(CanvGr)==ColName])==F & (CanvGr$CSS_Wader == "Y" | CanvGr$ESS_Wader  == "Y") & CanvGr$Reserve== "N")
    
    # Find the nearest polygon in Can2[Hab,] for each polygon in Can2[noHab,]
    nearest_indices <- st_nearest_feature(Can2[noHab,], Can2[Hab,])
    
    # Calculate the distances to the nearest polygons
    distances <- st_distance(Can2[noHab,], Can2[Hab,][nearest_indices, ], by_element = TRUE)
    
    ## If the nearest parcel is more than 1000m away then make it an NA
    distacneedit <- ifelse(as.numeric(distances)>1000, NA, nearest_indices)
    
    ## Now assign the new habitat values to the correct rows that are not too far from their nearest neighbour
    Data[noHab, colnames(Data)%in% ColName] <- Can2[Hab,][distacneedit, colnames(Can2)%in% ColName] |> st_drop_geometry()

    ## Remove geometry as can mess up sub-setting by columns
    Can2 <- Data |> st_drop_geometry()
    
    ## For AES only parcels missing a habitat value sample a habitat value from the same category
    InterS <- sample(x= na.omit((Can2[(Can2$CSS_Wader == "Y" | Can2$ESS_Wader  == "Y") & Can2$Reserve== "N" ,])[, colnames(Can2)==ColName]), 
                     size = nrow(Can2[is.na(Can2[,colnames(Can2)==ColName])==T & (Can2$CSS_Wader == "Y" | Can2$ESS_Wader  == "Y") & Can2$Reserve== "N",]), replace = T)
    ## Assign the sampled value to the rows with the missing habitat values
    Can2[is.na(Can2[,colnames(Can2)==ColName])==T & (Can2$CSS_Wader == "Y" | Can2$ESS_Wader  == "Y") & Can2$Reserve== "N" , colnames(Can2)==ColName] <- InterS
    
    ## Assign habitat values to missing rows
    Data[colnames(Data)==ColName] <- Can2[colnames(Can2)==ColName]
    
    
    ##
    
    
    ## 
    Can2 <- Data 
    
    ## Seperate the grassland opportunity parcels without habitat values to the ones with habitat values
    ## Capture any grassland opportunity parcels that might be just outside the priority landscape for Hab set
    noHab <- which(is.na(CanvGr[,colnames(CanvGr)==ColName])==T & CanvGr$Category == "Grass Opp")
    Hab <- which(is.na(CanvGr[,colnames(CanvGr)==ColName])==F & CanvGr$CSS_Wader == "N" & CanvGr$ESS_Wader  == "N" & CanvGr$Reserve == "N")
    
    # Find the nearest polygon in Can2[Hab,] for each polygon in Can2[noHab,]
    nearest_indices <- st_nearest_feature(Can2[noHab,], Can2[Hab,])
    
    # Calculate the distances to the nearest polygons
    distances <- st_distance(Can2[noHab,], Can2[Hab,][nearest_indices, ], by_element = TRUE)
    
    ## If the nearest parcel is more than 1000m away then make it an NA
    distacneedit <- ifelse(as.numeric(distances)>1000, NA, nearest_indices)
    
    ## Now assign the new habitat values to the correct rows that are not too far from their nearest neighbour
    Can2[noHab, colnames(Can2)%in% ColName] <- Can2[Hab,][distacneedit, colnames(Can2)%in% ColName] |> st_drop_geometry()
    
    ## Remove geometry as can mess up sub-setting by columns
    Can2 <- Data |> st_drop_geometry()
    
    ## For parcels that do not have AES or reserve and have missing a habitat value sample a habitat value from the same category
    InterS <- sample(x= na.omit((Can2[(Can2$CSS_Wader == "N" & Can2$ESS_Wader  == "N" & Can2$Reserve == "N"),])[, colnames(Can2)==ColName]), 
                     size = nrow(Can2[is.na(Can2[,colnames(Can2)==ColName])==T & Can2$CSS_Wader == "N" & Can2$ESS_Wader  == "N" & Can2$Reserve == "N",]), replace = T)
    ## Assign the sampled value to the rows with the missing habitat values
    Can2[is.na(Can2[,colnames(Can2)==ColName])==T & (Can2$CSS_Wader == "N" & Can2$ESS_Wader  == "N" & Can2$Reserve == "N") , colnames(Can2)==ColName] <- InterS
    
    ## Assign habitat values to missing rows
    Data[colnames(Data)==ColName] <- Can2[colnames(Can2)==ColName]
  
    
    
    ## Return the data set
    return(Data)} # end of function
  
  
  ## Plot all of the survey variables across grassland opportunity areas
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = TALL_BOUNDARY_PERCENT), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = GRASSLAND_TYPE), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = VEG_STRUCTURE), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = RUSH_PERCENT), colour = NA) +
  #            theme_minimal()
  
  
  ## Run the function for all survey variables
  CanvGr <- FillBlankHab(Data=CanvGr, ColName = "GRASSLAND_TYPE")
  CanvGr <- FillBlankHab(Data=CanvGr, ColName = "TALL_BOUNDARY_PERCENT")
  CanvGr <- FillBlankHab(Data=CanvGr, ColName = "STOCK")
  CanvGr <- FillBlankHab(Data=CanvGr, ColName = "VEG_STRUCTURE")
  CanvGr <- FillBlankHab(Data=CanvGr, ColName = "RUSH_PERCENT")
  
  
  ## Plot all of the survey variables across grassland opportunity areas
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = TALL_BOUNDARY_PERCENT), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = GRASSLAND_TYPE), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = VEG_STRUCTURE), colour = NA) +
  #            theme_minimal()
  # ggplot() + geom_sf(data = CanvGr, mapping = aes(geometry = geometry, fill = RUSH_PERCENT), colour = NA) +
  #            theme_minimal()
  
  ## Re-join the two halves of the data set back together
  Canv <- plyr::rbind.fill(CanvGr, CanvOther) 
  
  ## free up space
  rm(CanvGr, CanvOther); gc()
  
  
  
  
  ##---------------------------------------------------##
  #### F1.10 Annotate canvas with Crow/Magpie Density ####
  ##---------------------------------------------------##
  
  ## Split up opportunity areas (arable and grassland) from none-opportunity parcels
  CanvOpp <- filter(Canv, Category %in% c("Grass Opp", "Reserve", "AES Only", "Arable Opp"))
  CanvOther <- filter(Canv, !Category %in% c("Grass Opp", "Reserve", "AES Only", "Arable Opp"))
  
  
  ## Read in the UK 1km grid (goes with predator abundance data)
  Grid1KM <- st_read("RawData/Gridlines/OSGB_Grid_1km.shp")
  
  ## Read in the bird densities modeled from the BTO BBS survey
  ## These are provided as point at the center of each 1km square
  BTODens <- read_csv("RawData/BTO Density Data/model_popn_dens_0709.csv")
  
  ## streamline the BTO breeding bird density data
  BTODens <- BTODens |>  select(GridRef, MGdens, `C.dens`)
  
  ## Join the 1km square shapefile with the BTO data
  ## Each 1km square has a unique code
  Grid1KM <- left_join(Grid1KM, BTODens, by = join_by(PLAN_NO == GridRef))
  
  ## Convert the vector 1km grid into a raster 
  template = rast(vect(Grid1KM), res = 1000) # create temple
  
  
  ## Extract the density of Magpies and Crows for each 1km pixel
  Magpie <- rasterize(vect(Grid1KM), template, field = "MGdens") # magpie density
  # plot(Magpie, col=grDevices::hcl.colors(50, palette = "Sunset"))
  
  Crow <- rasterize(vect(Grid1KM), template, field = "C.dens") # Crow density
  # plot(Crow, col=grDevices::hcl.colors(50, palette = "Sunset"))
  rm(BTODens, Grid1KM, template); gc()
  
  
  ## Calculate the mean weighted Magpie density for each field
  ## Crop raster to study areas and fill in NA pixels in a buffer
  Magpie <- crop(Magpie, vect(st_buffer(st_as_sf(CanvOpp), dist = 1000))) # crop to save space
  Magpie <- focal(x = Magpie, w = 3, fun = "mean", na.policy = "only", na.rm = T)
  
  ## run Magpie extraction
  MagDens <- exact_extract(x = Magpie, y = st_as_sf(CanvOpp), fun = "mean")
  summary(MagDens); hist(MagDens) # checks
  
  ## assign density back to main data set
  CanvOpp$MagDens <- MagDens
  rm(Magpie, MagDens);gc() # clean environment
  
  
  ## Calculate the mean weighted Carrion Crow density for each field
  ## Crop raster to study areas and fill in NA pixels in a buffer
  Crow <- crop(Crow, vect(st_buffer(st_as_sf(CanvOpp), dist = 1000))) # crop to save space
  Crow <- focal(x = Crow, w = 3, fun = "mean", na.policy = "only", na.rm = T)
  
  ## run Crow extraction
  CrowDens <- exact_extract(x = Crow, y = st_as_sf(CanvOpp), fun = "mean")
  summary(CrowDens); hist(CrowDens) # checks
  
  ## assign density back to main data set
  CanvOpp$CrowDens <- CrowDens
  rm(Crow, CrowDens);gc() # clean environment
  
  ## Calculate the overall Corvid density
  CanvOpp <- CanvOpp |> 
            mutate(CorvDens = CrowDens + MagDens) |> 
            select(-c(CrowDens, MagDens))
  
  
  
  
  ##-------------------------------##
  #### F1.11 Intertidal Distance ####
  ##-------------------------------##
  
  ## Only run if user has set IntertidalDist to TRUE
  if(IntertidalDist==TRUE){
    
    ## Read in the UKCEH landcover data
    LC <- rast("RawData/LandCover/gblcm25m2021.tif")
    LC <- LC[[1]]
    
    ## crop the UK CEH data to south of England
    LC <- crop(LC, vect(st_as_sf(CanvOpp)))
    
    ## Reclassify the raster so that inter tidal habitats are 1 and all others habitats are NA
    m <- c(0, 12.5, NA,
           12.6, 13.5, 1,
           13.6, 14.5, NA, 
           14.6, 19.5, 1, 
           19.6, 22, NA) # matrix for re-classification
    rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
    rc1 <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
    # plot(rc1) # check reclassification has worked
    
    ## Calculate distance from NA pixels (None inter tidal) to none NA pixel (inter tidal habitat)
    DistToInter <- distance(rc1)
    # plot(DistToInter, main = "Distance to Intertidal (m)")
    
    ## run Magpie extraction
    InterDist <- exact_extract(x = DistToInter, y = st_as_sf(CanvOpp), fun = "mean")
    summary(InterDist); hist(InterDist) # checks
    
    ## assign density back to main data set
    CanvOpp$InterTidal_Distm <- InterDist
    rm(rc1, LC, InterDist, DistToInter);gc() # clean environment
    
  }


  
  
  ##-----------------------------------##
  #### F1.12 Predator Fence Presence ####
  ##-----------------------------------##
  
  ## Read in my spreadsheet of predator fences
  Fences <- read_csv("RawData/Predator Fences/RSPB Predator Fence Locaitons.csv") |> filter(is.na(F_LOC_ID)==F)
  
  ## Remove any duplicates
  IDX <- duplicated(Fences$F_LOC_ID)
  Fences <- Fences[IDX==F,]
  
  ## Manipulate predator fence to join with breeding pair estimate
  Fences <- Fences |> group_by(F_LOC_ID) |> summarise(Fence_Coverage = max(FENCE_COVERAGE))
  
  ## Join to the breeding pair estimates
  CanvOpp2 <- left_join(CanvOpp, Fences, by = c("F_LOC_ID"))
  stopifnot(nrow(CanvOpp2)==nrow(CanvOpp))
  CanvOpp <- CanvOpp2; rm(CanvOpp2)
  
  ## Change column to binary 1/0 column (1=fence)
  CanvOpp <- CanvOpp |> mutate(Fence_Coverage = ifelse(is.na(Fence_Coverage)==T, 0, Fence_Coverage)) # change the NAs to 0's
  
  ## free up space
  rm(Fences); gc()
  
  
  
  
  ##--------------------------------------------##
  #### F1.13 Identify Fields with Peaty Soils ####
  ##--------------------------------------------##
  
  if(is.null(PeatSoilsClass)==F){
    
    ## Read in Soil data here as it is large and coudl crash R
    Soils <- st_read("RawData/Soil/Soilscapes_England_Wales_27700.shp")
    
    ## filter out just the soils that are in some way peaty
    Soils <- filter(Soils, SOILSCAPE %in% PeatSoilsClass)
    
    ## Now label fields as peaty of they are more then 50% covered by the peaty soils polygon
    Soils <- st_as_sf(terra::aggregate(vect(Soils)))
    CanvOpp <- Field_Overlap_IDer(Fields = CanvOpp, OverlapShapes = Soils, 
                                  Thresh = 0.5, ColName = "Peat_Soil")
    
    ## Plot to check that this worked in Somerset
    # ggplot() + geom_sf(data = CanvOpp, aes(fill = Peat_Soil), colour = NA) + theme_light()
    rm(Soils); rm()
    
  }
  
  
  
  
  ##------------------------------------------------##
  #### F1.14 Calculate Landscape-scale Wood/Urban ####
  ##------------------------------------------------##
  
  ## Read in the UKCEH landcover data
  LC <- rast("RawData/LandCover/gblcm25m2021.tif")
  LC <- LC[[1]]
  
  
  
  ##-- Calculate percentage of woodland --##
  
  ## Using UKCEH land cover data (Codes 1/2) calculate the percentage woodland around each field
  ## Reclassify the raster so woodland pixels are 1 and all others habitats are 0
  m <- c(0, 2.5, 1,
         2.6, 22, 0) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  WoodLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
  plot(WoodLC) # check reclassification has worked
  
  ## Crop the binary raster to save on processing time
  WoodLC <- crop(WoodLC, vect(st_buffer(st_as_sf(CanvOpp), dist = 2050)))
  plot(WoodLC)
  
  ## Now calculate the proportion of woodland cover in each field buffer
  CanvOpp$PropWood_500 <- exact_extract(x = WoodLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 500), fun = "mean")
  CanvOpp$PropWood_1000 <- exact_extract(x = WoodLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 1000), fun = "mean")
  CanvOpp$PropWood_1500 <- exact_extract(x = WoodLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 1500), fun = "mean")
  CanvOpp$PropWood_2000 <- exact_extract(x = WoodLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 2000), fun = "mean")
  rm(WoodLC); gc()
  
  
  
  ##--  Calculate percentage of urban areas --##
  
  ## Using UKCEH land cover data (Codes 20/21) calculate the percentage Urban cover around each field
  ## Reclassify the raster so urban pixels are 1 and all others habitats are 0
  m <- c(0, 19.5, 0,
         19.6, 22, 1) # matrix for re-classification
  rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
  UrbanLC <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
  # plot(UrbanLC) # check reclassification has worked
  
  ## Crop the binary raster to save on processing time
  UrbanLC <- crop(UrbanLC, vect(st_buffer(st_as_sf(CanvOpp), dist = 2050)))
  # plot(UrbanLC)
  
  ## Now calculate the proportion of woodland cover in each field buffer
  CanvOpp$PropUrban_500 <- exact_extract(x = UrbanLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 500), fun = "mean")
  CanvOpp$PropUrban_1000 <- exact_extract(x = UrbanLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 1000), fun = "mean")
  CanvOpp$PropUrban_1500 <- exact_extract(x = UrbanLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 1500), fun = "mean")
  CanvOpp$PropUrban_2000 <- exact_extract(x = UrbanLC, y = st_as_sf(CanvOpp) |> st_centroid() |> st_buffer(dist = 2000), fun = "mean")
  rm(UrbanLC); gc()
  
  ## Re-join the two halves of the data set back together
  Canv <- plyr::rbind.fill(CanvOpp, CanvOther) 
  
  
  
  
  ##----------------------------------------##
  #### F1.15 Average standing water value ####
  ##----------------------------------------##
  
  ## I want to label standing water in all the fields that are within my opportunity area and are
  ## Grassland opp, AES or reserve. I also want to label those extra grassland fields that are just outside my priority landscape
  ## but are still potentially lowland wet grassland
  PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp") | Canv$WiderWetGrass == 1)

  ## Calculate the mean weighted pixel value for each field
  # WaterVals <- exact_extract(x = WatRast, y = st_as_sf(st_buffer(st_as_sf(CanvOpp[PrIndex, ]), dist = -10)), fun = "mean") this faster package doesn't work here
  # y = st_as_sf(st_buffer(st_as_sf(CanvOpp), dist = -10))
  # table(is.na(sf::st_dimension(y) == 2))
  WaterVals <- terra::extract(WatRast, vect(st_buffer(st_as_sf(Canv[PrIndex, ]), dist = -10)), fun = mean, na.rm=TRUE)
  hist(WaterVals$mNDWI)
  
  ## Assign these weighted means to the main data set
  Canv$WaterCoverage <- NA
  Canv$WaterCoverage[PrIndex] <- WaterVals$mNDWI
  
  ## Change NaNs to NA
  Canv$WaterCoverage <- ifelse(Canv$WaterCoverage == "NaN", NA, Canv$WaterCoverage)
  rm(WaterVals, WatRast); gc()
  
  ## plot the extracted water values
  # ggplot() + geom_sf(data = st_as_sf(Canv), aes(fill = WaterCoverage), colour = NA) + theme_light()
  
  
  
  
  ##-------------------------------------------------##
  #### F1.16 Calculate Landscape-scale Water/Grass ####
  ##-------------------------------------------------##
  
  ##-- Calculate average standing water value across landscape --##
  
  ## Crop the UKCEH landcover to Somerset to use as a base raster
  LC_Som <- crop(LC, vect(st_buffer(st_as_sf(Canv), dist = 2050)))
  plot(LC_Som)
  
  ## Makse sure the canvas is a shape file
  Canv <- Canv |> st_as_sf()
  
  ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
  ## The base raster is just the UKCEH land cover data set at 25m res
  ## First create index of rows for fields that are suitable
  Suit <- which(Canv$Opp == "Grass" | Canv$WiderWetGrass == 1) # use this as some areas masked out are still functional for waders just not possible for change
  WaterRast <- rasterize_polygons(Canv[Suit, ], LC_Som, min_coverage = 0.25)
  
  ### Now assign the pixels a water coverage value from the correct field
  values(WaterRast) <- Canv[Suit,]$WaterCoverage[values(WaterRast)]
  plot(WaterRast)
  
  ## Get the index of the rows that I want to calculate landscape variables for
  ## These will be the ones I want to predict wader into
  PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
  
  ## Update each of the landscape scale water coverage variables
  ## The fun = "Mean does the following:  the mean cell value, weighted by the fraction of each cell that is covered by the polygon
  ## See: https://isciences.gitlab.io/exactextractr/reference/exact_extract.html
  Canv$Ave_WiderWater500 <- NA
  Canv$Ave_WiderWater500[PrIndex] <- exact_extract(x = WaterRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 500)), 
                                                     fun = "mean")
  
  Canv$Ave_WiderWater1000 <- NA
  Canv$Ave_WiderWater1000[PrIndex] <- exact_extract(x = WaterRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 1000)), 
                                                     fun = "mean")
  
  Canv$Ave_WiderWater1500 <- NA
  Canv$Ave_WiderWater1500[PrIndex] <- exact_extract(x = WaterRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 1500)), 
                                                     fun = "mean")
  
  Canv$Ave_WiderWater2000 <- NA
  Canv$Ave_WiderWater2000[PrIndex] <- exact_extract(x = WaterRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 2000)), 
                                                     fun = "mean")
  
  
  
  ##-- Calculate percentage of semi-nat wet grassland --##
  
  ## Do the same thing for lowland wet grassland
  Suit <- which(Canv$Opp == "Grass" | Canv$WiderWetGrass == 1) # use this as some areas masked out are still functional for waders just not possible for change
  GrassRast <- rasterize_polygons(Canv[Suit, ], LC_Som, min_coverage = 0.5)
  
  ### Now assign the pixels a water coverage value from the correct field
  values(GrassRast) <- ifelse(is.na(values(GrassRast)==T), 0, 1)
  # plot(GrassRast)
  
  ## Get the index of the rows that I want to calculate landscape variables for
  ## These will be the ones I want to predict wader into
  PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
  
  ## Update each of the landscape scale water coverage variables
  ## The fun = "Mean does the following:  the mean cell value, weighted by the fraction of each cell that is covered by the polygon
  ## See: https://isciences.gitlab.io/exactextractr/reference/exact_extract.html
  Canv$PropWetGrass_500 <- NA
  Canv$PropWetGrass_500[PrIndex] <- exact_extract(x = GrassRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 500)), 
                                                     fun = "mean")
  
  Canv$PropWetGrass_1000 <- NA
  Canv$PropWetGrass_1000[PrIndex] <- exact_extract(x = GrassRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 1000)), 
                                                     fun = "mean")
  
  Canv$PropWetGrass_1500 <- NA
  Canv$PropWetGrass_1500[PrIndex] <- exact_extract(x = GrassRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 1500)), 
                                                     fun = "mean")
  
  Canv$PropWetGrass_2000 <- NA
  Canv$PropWetGrass_2000[PrIndex] <- exact_extract(x = GrassRast, 
                                                     y = (Canv[PrIndex, ] |> st_centroid() |> st_buffer(dist = 2000)), 
                                                     fun = "mean")
  
  
  
  ##--------------------------------------------##
  #### F1.17 Calculate Purchase/Forgone Costs ####
  ##--------------------------------------------##

  ## First need to give each parcel an agricultural grading, cine the costs could vary based on the grading
  ## Grades of different agricultural land
  AgriGrade <- st_read("RawData/AgriGrades/Agricultural_Land_Classification_Provisional_EnglandPolygon.shp")
  
  # ## Add new column for scoring of grading, this allows me to rasterize the grading polygons
  AgriGrade <- AgriGrade |> 
          mutate(Grade = case_when(alc_grade == "Grade 1" ~ 1,
                                    alc_grade == "Grade 2" ~ 2,
                                    alc_grade == "Grade 3" ~ 3,
                                    alc_grade == "Grade 4" ~ 4,
                                    alc_grade == "Grade 5" ~ 5,
                                    .default = NA))
  
  ## rasterize the agri grading polygons
  AgriGrades <- rasterize(vect(AgriGrade), LC_Som, field = "Grade", fun = "min", background = NA)
  plot(AgriGrades)
  
  ## If a field is covered by the raster then take the minimum gradings, only do this for field that could beocme reserve in a scenario
  AllIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp", "Arable Opp"))
  Canv$AgriGrade <- NA
  Canv$AgriGrade[AllIndex] <- exact_extract(x = AgriGrades, 
                                           y = Canv[AllIndex, ], 
                                           fun = "min")
  
  ## If a parcel did not receive a grading and it is an opportunity field then just give it  grade of 3 (the middle)
  Canv$AgriGrade[AllIndex] <- ifelse(is.na(Canv$AgriGrade[AllIndex])==T, 3, Canv$AgriGrade[AllIndex])
  
  ## plot the output for now to check
  ggplot(data=Canv[AllIndex,]) + geom_sf(mapping=aes(geometry=geometry, fill=AgriGrade))
  
  ## Calculate the rough purchase cost of the land parcel
  ## this is based on the region, habitat and grading
  ## Source here: https://rural.struttandparker.com/article/english-estates-farmland-market-review-winter-2023-24/
  Purchase <- read.csv("RawData/AES Costings/Land purchase costs.csv")
  
  ## select out the correct rows for the region of chouce
  if(Label=="Somerset Levels and Moors"){P1 <- Purchase |> filter(Region == "South West")}
  if(Label=="North Kent"){P1 <- Purchase |> filter(Region == "South East")}
  if(Label=="Essex"){P1 <- Purchase |> filter(Region == "South East")}
  if(Label=="Broads"){P1 <- Purchase |> filter(Region == "East England")}
  
  ## Now for different land categories and grading assign a land cost
  Canv <- Canv |> 
          mutate(Hectares = as.numeric(m2_to_ha(ParcArea)),
                 PurchaseGDP = case_when(Category == "AES Only" & AgriGrade <=2 ~ P1$Q75[P1$Land_Type=="Pasture"],
                                         Category == "AES Only" & AgriGrade ==3 ~ P1$Q50[P1$Land_Type=="Pasture"],
                                         Category == "AES Only" & AgriGrade >=4 ~ P1$Q25[P1$Land_Type=="Pasture"],
                                         Category == "Grass Opp" & AgriGrade <=2 ~ P1$Q75[P1$Land_Type=="Pasture"],
                                         Category == "Grass Opp" & AgriGrade ==3 ~ P1$Q50[P1$Land_Type=="Pasture"],
                                         Category == "Grass Opp" & AgriGrade >=4 ~ P1$Q25[P1$Land_Type=="Pasture"],
                                         Category == "Arable Opp" & AgriGrade <=2 ~ P1$Q75[P1$Land_Type=="Arable"],
                                         Category == "Arable Opp" & AgriGrade ==3 ~ P1$Q50[P1$Land_Type=="Arable"],
                                         Category == "Arable Opp" & AgriGrade >=4 ~ P1$Q25[P1$Land_Type=="Arable"],
                                         Category == "Reserve" ~ 0,
                                         .default = NA),
                 PurchaseGDP = PurchaseGDP*Hectares) |> select(-Hectares)

  
  
  ## Calculate the income foregone, this is based on the land grading, current habitat and AES agreements
  ## Read in the data sheet with the income foregone for various scenarios
  FG <- read.csv("RawData/AES Costings/Income foregone summary.csv")
  
  ## First need to label which columns are already in breeding or wintering wader schemes
  ## This in schemes already will have lower rates of income foregone
  WhichAES <- which(Canv$Category == "AES Only")
  Canv$BreedScheme <- NA; Canv$WintScheme <- NA
  Canv$BreedScheme[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% c("GS9", "GS11"))})
  Canv$WintScheme[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% c("GS10", "GS12"))})

  ## Now calcualte the income foregone for different AES schemes
  ## For agricualtural land of varying qualities
  ## Might want to also vary the income forgone based on gradings for Category == "Grass Opp"
  Canv <- Canv |> 
        mutate(Hectares = as.numeric(m2_to_ha(ParcArea)),
               ForegoneGDP = case_when(Category == "AES Only" & AgriGrade <=2 ~ FG$Foregone[FG$Land_Type=="Lowland Ley + Permanent Pasture + Dairy"],
                                       Category == "AES Only" & AgriGrade ==3 ~ FG$Foregone[FG$Land_Type=="Lowland Permanent Pasture + Dairy"],
                                       Category == "AES Only" & AgriGrade >=4 ~ FG$Foregone[FG$Land_Type=="Lowland Permanent Pasture"],
                                       Category == "Grass Opp" & AgriGrade <=2 ~ FG$Foregone[FG$Land_Type=="Lowland Ley + Permanent Pasture + Dairy"],
                                       Category == "Grass Opp" & AgriGrade ==3 ~ FG$Foregone[FG$Land_Type=="Lowland Permanent Pasture + Dairy"],
                                       Category == "Grass Opp" & AgriGrade >=4 ~ FG$Foregone[FG$Land_Type=="Lowland Permanent Pasture"],
                                       Category == "Arable Opp" & AgriGrade <=2 ~ FG$Foregone[FG$Land_Type=="High yield arable"],
                                       Category == "Arable Opp" & AgriGrade ==3 ~ FG$Foregone[FG$Land_Type=="Medium yield arable"],
                                       Category == "Arable Opp" & AgriGrade >=4 ~ FG$Foregone[FG$Land_Type=="Low yield arable"],
                                       Category == "Reserve" ~ 0,
                                       .default = NA),
               ForegoneGDP = ForegoneGDP*Hectares) |> select(-c(Hectares, BreedScheme, WintScheme))
  
  
  rm(AgriGrade, AllIndex, FG, WhichAES, Purchase, P1); gc()
    
  
  
  ##---------------------------##
  #### F1.18 Read out Canvas ####
  ##---------------------------##
  
  ## Assign the landscape variable to the data set
  Canv$Landscape <- Label
  
  ## Read out canvas as both a shape file and a rds file to preserve the list column of AES options
  Canv |> st_drop_geometry() |> 
          write_rds(paste0("CleanData/Scenarios/4-AnnotateCanvas/", outpath, ".rds"))
          #write_csv(paste0("CleanData/Scenarios/4-AnnotateCanvas/", outpath, ".csv"))
  
  Canv |> select(-AES_options) |> st_as_sf() |> 
  write_sf(paste0("CleanData/Scenarios/4-AnnotateCanvas/", outpath, ".shp"), overwrite = T, append=FALSE)
  
  
  
} #### END of AnnotateCanv ####





##-------------------------------##
#### Annotate: Somerset Levels ####
##-------------------------------##

## Read in canvas from script 1
Canvas <- st_read("CleanData/Scenarios/1-Starting Canvas/Som_Canvas.shp")

## Read in a shape file of wader sites for region of interest
WaderClust <- st_read("CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Sites.shp")
plot(WaderClust$geometry)

## Read in raster data set that define opportunity area for grassland and arable for ROI
SomGrassOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Som_LowWetGrass.tif")
SomArableOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Som_ArableSuitable.tif")

## Read in rasters of predicted standing water for ROI
SomW <- rast("RawData/PredictedWaterRasters/Som_StandingWaterAll_alt.tif")

## Define part of the path to read in the stakeholder gradings
Sompaths = "Somerset/Som"

## Read in the NATMAP soil vector categories that are peaty
SomPeatClass <- c("Blanket bog peat soils", "Fen peat soils", "Loamy and sandy soils with naturally high groundwater and a peaty surface",
                    "Raised bog peat soils", "Slowly permeable wet very acid upland soils with a peaty surface",
                    "Very acid loamy upland soils with a wet peaty surface")

## Assign the landscape name so the canvas can be labelled
SomLabel <- "Somerset Levels and Moors"

## Read out canvas as both a shape file and csv file
SomOut <- "Som_AnnotatedCanv"


## Run function for Somerset
AnnotateCanv(Canv = Canvas, AllClust = WaderClust,
             GrassOpp = SomGrassOpp, ArableOpp = SomArableOpp,
             PeatSoilsClass = SomPeatClass, WatRast = SomW, IntertidalDist=FALSE,
             Label = SomLabel, guidepath = Sompaths, outpath = SomOut, RunGrades=TRUE)

## free up space
rm(Canvas, WaderClust, SomGrassOpp, SomArableOpp, SomW, SomPeatClass, Sompaths, SomLabel, SomOut); gc()




##--------------------------##
#### Annotate: North Kent ####
##--------------------------##

## Read in canvas from script 1
KentCanvas <- st_read("CleanData/Scenarios/1-Starting Canvas/NKent_Canvas.shp")

## Read in a shape file of wader sites for region of interest
NKWaderClust <- st_read("CleanData/Scenarios/2-DefineWaderSites/Kent/Kent_Wader_Sites.shp")
plot(NKWaderClust$geometry)

## Read in raster data set that define opportunity area for grassland and arable for ROI
NKentGrassOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_LowWetGrass.tif")
NKentArableOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/NKent_ArableSuitable.tif")

## Read in rasters of predicted standing water for ROI
NKW <- rast("RawData/PredictedWaterRasters/NKent_StandingWaterAll.tif")

## Define part of the path to read in the stakeholder gradings
## start of file path is "CleanData/Guideline Creation/"
NKpaths = "Kent/Kent"

## Assign the landscape name so the canvas can be labelled
NKLabel <- "North Kent"

## Read out canvas as both a shape file and csv file
NKOut <- "NKent_AnnotatedCanv"

## Run function for North Kent
AnnotateCanv(Canv = KentCanvas, AllClust = NKWaderClust,
             GrassOpp = NKentGrassOpp, ArableOpp = NKentArableOpp,
             PeatSoilsClass = NULL, WatRast = NKW, IntertidalDist=TRUE,
             Label = NKLabel, guidepath = NKpaths, outpath = NKOut, RunGrades=TRUE)

ALLNK <- NKentGrassOpp + NKentArableOpp 
plot(ALLNK)

## free up space
rm(KentCanvas, NKWaderClust, NKentGrassOpp, NKentArableOpp, NKW, NKpaths, NKLabel, NKOut); gc()






##---------------------------##
#### Annotate: Essex Coast ####
##---------------------------##

## Read in canvas from script 1
EssexCanvas <- st_read("CleanData/Scenarios/1-Starting Canvas/Essex_Canvas.shp")

## Read in a shape file of wader sites for region of interest
EsWaderClust <- st_read("CleanData/Scenarios/2-DefineWaderSites/Essex/Essex_Wader_Sites.shp")
plot(EsWaderClust$geometry)

## Read in raster data set that define opportunity area for grassland and arable for ROI
EssexGrassOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_LowWetGrass.tif")
EssexArableOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Essex_ArableSuitable.tif")

## Read in rasters of predicted standing water for ROI
EsW1 <- rast("RawData/PredictedWaterRasters/Esx1_StandingWaterAll.tif")
EsW2 <- rast("RawData/PredictedWaterRasters/Esx2_StandingWaterAll.tif")
EsW <- mosaic(EsW1, EsW2)
plot(EsW)
rm(EsW1, EsW2);gc()
  
## Define part of the path to read in the stakeholder gradings
## start of file path is "CleanData/Guideline Creation/"
Espaths = "Essex/Ess"

## Assign the landscape name so the canvas can be labelled
EsLabel <- "Essex"

## Read out canvas as both a shape file and csv file
EsOut <- "Essex_AnnotatedCanv"

## Run function for Essex
AnnotateCanv(Canv = EssexCanvas, AllClust = EsWaderClust,
             GrassOpp = EssexGrassOpp, ArableOpp = EssexArableOpp,
             PeatSoilsClass = NULL, WatRast = EsW, IntertidalDist=TRUE,
             Label = EsLabel, guidepath = Espaths, outpath = EsOut, RunGrades=TRUE)

AllEssex <- EssexGrassOpp+EssexArableOpp 
plot(AllEssex)

## free up space
rm(EssexCanvas, EsWaderClust, EssexGrassOpp, EssexArableOpp, EsW, Espaths, EsLabel, EsOut); gc()




##------------------------------##
#### Annotate: Norfolk Broads ####
##------------------------------##

## Read in canvas from script 1
BroadsCanvas <- st_read("CleanData/Scenarios/1-Starting Canvas/Broads_Canvas.shp")

## Read in a shape file of wader sites for region of interest
BrWaderClust <- st_read("CleanData/Scenarios/2-DefineWaderSites/Broads/Broads_Wader_Sites.shp")
plot(BrWaderClust$geometry)

## Read in raster data set that define opportunity area for grassland and arable for ROI
BroadsGrassOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_LowWetGrass.tif")
BroadsArableOpp <- rast("CleanData/Scenarios/3-DefineActionAreas/Broads_ArableSuitable.tif")

## Read in rasters of predicted standing water for ROI
BrW <- rast("RawData/PredictedWaterRasters/Broads_StandingWaterAll.tif")
  
## Define part of the path to read in the stakeholder gradings
## start of file path is "CleanData/Guideline Creation/"
Brpaths = "Norfolk/Broads"

## Assign the landscape name so the canvas can be labelled
BrLabel <- "Broads"

## Read out canvas as both a shape file and csv file
BrOut <- "Broads_AnnotatedCanv"

## Run function for Broads
AnnotateCanv(Canv = BroadsCanvas, AllClust = BrWaderClust,
             GrassOpp = BroadsGrassOpp, ArableOpp = BroadsArableOpp,
             PeatSoilsClass = NULL, WatRast = BrW, IntertidalDist=TRUE,
             Label = BrLabel, guidepath = Brpaths, outpath = BrOut, RunGrades=TRUE)

## Read in raster data set that define opportunity area for grassland and arable for ROI
ALLB <- BroadsGrassOpp + BroadsArableOpp 
plot(ALLB)


