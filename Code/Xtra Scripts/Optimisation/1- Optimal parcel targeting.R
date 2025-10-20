##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 22/04/2024
## 
## 
## Aim: Determine which land parcels should be improved 
##      in order to cause the largest increases in waders
## 
##------------------------------------------------------##


## Load in packages & helper functions
pacman::p_load(tidyverse, # used
               sf, # used
               terra, # used
               randomForestSRC, # used
               exactextractr, # used
               ggpubr, # used
               truncnorm,
               doFuture) # used

options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R")



##---------------------------------##
## F 1 START OptimiseTarg Function ##
##---------------------------------##

# Code for testing function
# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# NewCat= "Reserve"
# Plus = F
# Outpath = "CleanData/Scenarios/6-OpimalTargeting/Kent_OptimalTargeting"
# FenceParcels = FALSE
# LoopStart = 1
# LoopStop = 100


## Defining function arguments
# Canvas= landscape canvas as a shape file variables need to match those in the RF model
# SniModel= RF model for predicting Snipe abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LapModel= RF model for predicting Lapwing abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# RedModel= RF model for predicting Redshank abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LandCov= Landcover data raster covering same areas as Canvas. This is used to create as raster out of Canvas to calculate landscape variables
# FenceParcels = Should the selected parcel have a fenc round it (TRUE) or be unfences (FALSE)
# NewCat= Field parcel categories that are used to `improve` the chosen opportunity fields to (can use c("AES Only"), or c("Reserve"))
# Plus= (TRUE/FALSE) should habitat plus improvement be undertaken? Where opportunity fields are improved to the standard of fields with waders for the chosen management (NewCat)
# Outpath = Outpath for the log file
OptimiseTarg <- function(Canvas, 
                          SniModel, 
                          LapModel,
                          RedModel, 
                          LandCov,
                          FenceParcels,
                          NewCat, 
                          Plus, 
                          LoopStart = NA,
                          LoopStop = NA,
                          Outpath){

  

  ##---------------------------##
  #### F 1.1 Scenario set up ####
  ##---------------------------##
  
  ## Remove all the fields that are not needed during the scenario modelling
  ## This is to try and increase the speed
  ## Any fields with no opportunity or that were masked are removed unless they are wider wet grassland adjacent to the priority landscape
  Canvas <- filter(Canvas, Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve") | WiderWetGrass == 1)
  ggplot(data = Canvas) + geom_sf(mapping = aes(fill=Category)) + theme_light()
  
  ## Create data frame to record wader abundance at the end of each segment
  ## Also record the starting abundance in a seperate columns
  if(LoopStart == 1 | is.na(LoopStart) ==T){  
    ScTracker <- data.frame(ParcRef = Canvas$ParcRef,
                            SnipeNew = NA,
                            #Snipe_unfenced = NA,
                            BaseSnipe = ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T)),
                            LapwingNew = NA,
                            #Lapwing_unfenced = NA,
                            BaseLapwing = ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T)),
                            RedshankNew = NA,
                            #Redshank_unfenced = NA,
                            BaseRedshank = ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T)))
    }else{ScTracker <- read_csv(paste0(Outpath, ".csv"))}

  
  

  ## Create second species abundance columns
  Canvas$SnipeAbundAlt <- Canvas$SnipeAbund
  Canvas$LapAbundAlt <- Canvas$LapAbund
  Canvas$RedAbundAlt <- Canvas$RedAbund

  ## Pre-compute buffer for all fields so I just have to do it once
  CanvCants <- Canvas |> st_centroid() |> select(geometry)
  Parc500 <- CanvCants |> st_buffer(dist = 500)
  Parc1000 <- CanvCants |> st_buffer(dist = 1000)
  Parc1500 <- CanvCants |> st_buffer(dist = 1500)
  Parc2000 <- CanvCants |> st_buffer(dist = 2000)
  
  ## Store a copy of the unaltered canvas here
  ## Can then call on this to reset the canvas at the start of every loop
  MASTERCanv <- Canvas
  
  
  
  
  ##------------------------------------##
  #### F 1.2 Summarise target habitat ####
  ##------------------------------------##
  
  ## Get the row numbers of fields with the management option that new opportunity fields will become
  ## for plus sampling remove fields that do not have waders in
  if(Plus==F){IndNewHabs <- which(Canvas$Category %in% NewCat)}
  if(Plus==T){IndNewHabs <- which(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)}
  
  ## define which within habitat categories I am going to update
  ## also update the categopry for that row
  ## For Arable also update the TALL_BOUNDARY_PERCENT, normally for grassland fields I leave this the same but in arable it needs to be updates
  HabCats <- c("Category", "GrassOpp", "GRASSLAND_TYPE", "WaterCoverage", "STOCK", "VEG_STRUCTURE", "RUSH_PERCENT", "TALL_BOUNDARY_PERCENT")
  # if(any(OppCat %in% "Arable Opp")){HabCats <- c(HabCats, "TALL_BOUNDARY_PERCENT")}
  
  ## Extract the field for the category that I want to sample from
  ## And just the colums that I am going to update in the scenario
  ImproveHabs <- Canvas[IndNewHabs, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
  glimpse(ImproveHabs)
  
  ## Summarise the habitat values for the category
  ## These will be used to update the habitat for all fields in tern
  SumHabs <- ImproveHabs |> 
             summarise(Category = getmode(Category),
                       GrassOpp = getmode(GrassOpp),
                       GRASSLAND_TYPE = getmode(GRASSLAND_TYPE),
                       WaterCoverage = mean(WaterCoverage),
                       STOCK = getmode(STOCK),
                       VEG_STRUCTURE = getmode(VEG_STRUCTURE),
                       RUSH_PERCENT = mean(RUSH_PERCENT),
                       Fence_Coverage = getmode(VEG_STRUCTURE),
                       TALL_BOUNDARY_PERCENT = mean(TALL_BOUNDARY_PERCENT))
  
  ## Work out what column order my summary data frame needs to match that of the Canvas
  ReIndex <- match(colnames(Canvas[1, colnames(Canvas) %in% HabCats] |> st_drop_geometry()), colnames(SumHabs))
  
  ## Re-order, the summarised habitats so that the column order matches that of the main canvas
  SumHabs <- SumHabs[, ReIndex]
  
  
  
  
  ##------------------------------------------------##
  #### F 1.3 Loop through each opportunity field ####
  ##------------------------------------------------##
  
  ## Code to allow the loop to be started and stooped in specific loops
  if(is.na(LoopStop)==TRUE){LoopMax <- nrow(Canvas)}
  if(is.na(LoopStop)==FALSE){LoopMax <- LoopStop}
  if(is.na(LoopStart)==TRUE){LoopMin <- 1}
  if(is.na(LoopStart)==FALSE){LoopMin <- LoopStart}
  
  
  
  ##-- LOOP STARTS HERE --##
  
  ## For each of my random sample fields update the habitat within the fields and wider landscape water
  ## Then re-calculate the wader abundance variable
  for(i in LoopMin:LoopMax){ # 
    
    ## send message to console
    message("Running for field ", i, " out of ", nrow(Canvas))
    
    ## If the category is "NoOpp" then skip this field as we do not want to calculate the abundance
    if(Canvas$Category[i] == "NoOpp"){ next }
    
    ## save the original category of the field as this will get updated
    ## This will be used to update the landscape scale wet grassland value
    OrigCat <- Canvas$Category[i]
        

    
    ##----------------------------------------##        
    #### F 1.4 Update within field habitats ####
    ##----------------------------------------##
    
    ## send message to console
    message("Updating field habitat variables...")

    ## Assign new habitat values to the opportunity fields
    Canvas[i, colnames(Canvas) %in% HabCats] <- SumHabs
    
    ## Decide whether the selected parcel should be fenced or not
    ## this is based on the value of a argument in the function
    if(FenceParcels == TRUE){ Canvas$Fence_Coverage[i] <- "Y" }
    if(FenceParcels == FALSE){ Canvas$Fence_Coverage[i] <- "N" }
    

    
    ##------------------------------------##    
    #### F 1.5 Update landscape habitat ####
    ##------------------------------------##
    
    ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
    ## The base raster is just the UKCEH land cover data set at 25m res
    ## First create index of rows for fields that are suitable
    Suit <- which(Canvas$GrassOpp == 1 | Canvas$WiderWetGrass == 1)
    WaterRast <- rasterize_polygons(Canvas[Suit, ], LandCov, min_coverage = 0.25)
    
    ### Now assign the pixels a water coverage value from the correct field
    values(WaterRast) <- Canvas[Suit,]$WaterCoverage[values(WaterRast)]
    # plot(WaterRast)
    
    ## Work out which fields need to be updated, only a portion of fields will need to be updated as not all fields
    ## will be within 2000m of a field that has had it's water value changed
    UpdateParcels <- st_intersection(Canvas, (Canvas[i,] |> st_buffer(dist = 2050) |> st_union()))
    
    ## label those fields that need to be updated
    ## These are fields within 2000m of any field that have changed category
    ## and those fields that are reserve, AES only or grass
    ## NOTE might need to filter UpdateParcels
    Canvas$UpdateScape <- NA
    Canvas$UpdateScape <- ifelse(Canvas$ParcRef %in% c(UpdateParcels$ParcRef) & Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp"), 1, 0)
    
    ## get row number of fields that need to be updated
    IndUpdate <- which(Canvas$UpdateScape ==1)
    
    ## Update each of the landscape scale water coverage variables
    ## Assign updates values only to the rows identified above
    ## Calculate Centroid and then buffer any fields where landscape scale water needs updating
    ## Then extract the average water value from each buffered field and assign it back to the main data set
    ## The fun = "Mean does the following:  the mean cell value, weighted by the fraction of each cell that is covered by the polygon
    ## See: https://isciences.gitlab.io/exactextractr/reference/exact_extract.html
    Canvas$Ave_WiderWater500[IndUpdate] <- exact_extract(x = WaterRast, y = Parc500[IndUpdate, ], fun = "mean")
    
    Canvas$Ave_WiderWater1000[IndUpdate] <- exact_extract(x = WaterRast, y = Parc1000[IndUpdate, ], fun = "mean")
    
    Canvas$Ave_WiderWater1500[IndUpdate] <- exact_extract(x = WaterRast, y = Parc1500[IndUpdate, ], fun = "mean")
    
    Canvas$Ave_WiderWater2000[IndUpdate] <- exact_extract(x = WaterRast, y = Parc2000[IndUpdate, ], fun = "mean")
    
    
    ## Option to update extent of lowland wet grassland
    if(OrigCat == "Arable Opp"){
      
      ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
      ## The base raster is just the UKCEH land cover data set at 25m res
      ## First create index of rows for fields that are suitable
      Suit <- which(Canvas$GrassOpp == 1 | Canvas$WiderWetGrass == 1)
      GrassRast <- rasterize_polygons(Canvas[Suit, ], LandCov, min_coverage = 0.5)
      
      ### Now assign the pixels as either low wet grass (1) or not (0)
      values(GrassRast) <- ifelse(is.na(values(GrassRast)==T), 0, 1)
      # plot(GrassRast)
      
      ## Update each of the landscape scale water coverage variables
      ## Assign updates values only to the rows identified above
      ## Calculate Centroid and then buffer any fields where landscape scale water needs updating
      ## Then extract the average water value from each buffered field and assign it back to the main data set
      ## The fun = "Mean does the following:  the mean cell value, weighted by the fraction of each cell that is covered by the polygon
      ## See: https://isciences.gitlab.io/exactextractr/reference/exact_extract.html
      
      Canvas$PropWetGrass_500[IndUpdate] <- exact_extract(x = GrassRast, y = Parc500[IndUpdate, ], fun = "mean")
      
      Canvas$PropWetGrass_1000[IndUpdate] <- exact_extract(x = GrassRast, y = Parc1000[IndUpdate, ], fun = "mean")
      
      Canvas$PropWetGrass_1500[IndUpdate] <- exact_extract(x = GrassRast, y = Parc1500[IndUpdate, ], fun = "mean")
      
      Canvas$PropWetGrass_2000[IndUpdate] <- exact_extract(x = GrassRast, y = Parc2000[IndUpdate, ], fun = "mean")
    }
    

    ##-----------------------------------##
    #### F 1.6 Predict Wader abundance ####
    ##-----------------------------------##
    
    ## send message to console
    message("Updating Wader Abundance")
    
    ## Get the index of the rows that I want to predict wader abundance into
    Index <- which(Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp") & Canvas$UpdateScape ==1)
    
    ## Predict wader abundance and then assign it to the correct abundance column
    ## If there is no model object provided then do not predict for that species
    ## Predict for Snipe
    if(is.null(SniModel)==F){
        set.seed(1212)
        Pred <- (predict.rfsrc(object = SniModel, 
                              newdata = (Canvas[Index, ] |> select(SniModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                              jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$SnipeAbund[Index]))stop("Prediction error: incorrect number of Snipe fields predicited")
        Canvas$SnipeAbund[Index] <- Pred
        ScTracker$SnipeNew[i] <- sum(Canvas$SnipeAbund, na.rm = T)
        } # number of Snipe in the landscape
      
    
    ## Predict for Lapwing
    if(is.null(LapModel)==F){
        set.seed(1212)
        Pred <- (predict.rfsrc(object = LapModel, 
                                newdata = (Canvas[Index, ] |> select(LapModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$LapAbund[Index]))stop("Prediction error: incorrect number of Lapwing fields predicited")
        Canvas$LapAbund[Index] <- Pred
        ScTracker$LapwingNew[i] <- sum(Canvas$LapAbund, na.rm = T)
        } # number of Lapwing in the landscape 
      
    
    ## Predict for Redshank
    if(is.null(RedModel)==F){
        set.seed(1212)
        Pred <- (predict.rfsrc(object = RedModel, 
                                newdata = (Canvas[Index, ] |> select(RedModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$RedAbund[Index]))stop("Prediction error: incorrect number of Redshank fields predicited")
        Canvas$RedAbund[Index] <- Pred
        ScTracker$RedshankNew[i] <- sum(Canvas$RedAbund, na.rm = T)
        } # number of Lapwing in the landscape 
  
    
    
    
    ##------------------------------------------##
    #### F 1.7 Save current iteration outputs ####
    ##------------------------------------------##
    
    ## Read out shape file of current scenario state
    ## Might want to only save a subset of columns that might be important for plotting
    ## Alternatively could just save the final output as a shapefile
    write_csv(ScTracker, paste0(Outpath, ".csv")) # save the tracker log
    
    ## return the tracker data set
    glimpse(ScTracker[i,])
    
    ## Then reset the canvas at this point, this will mean that the in each round of the i loop we go back to the starting canvas
    ## and do not start from where we left off from on the last loop
    Canvas <- MASTERCanv
    
    
  } # End of segment loop (i)
  
} 






##-------------------------------##
##-------------------------------##
#### 1. *North Kent Scenarios* ####
##-------------------------------##
##-------------------------------##

##----------------------##
#### 1.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)
table(duplicated(Canvshp$ParcRef))

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")
table(duplicated(Canv$ParcRef))

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)

## Visualize canvas categories
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
           theme_minimal()


## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds")




##---------------------------------##
#### 1.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------------##
#### 1.3 Run optimization function ####
##-----------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# NewCat= "Reserve"
# Plus = F
# Outpath = "CleanData/Scenarios/6-OpimalTargeting/Kent_OptimalTargeting"
# FenceParcels = FALSE

## Run function for North Kent
OptimiseTarg(Canvas=Canv,
             SniModel=NULL, LapModel=LapMod, RedModel=RedMod,
             LandCov=LC_crop, NewCat= "Reserve", Plus = T, FenceParcels = FALSE,
             LoopStart = 1,
             Outpath = "CleanData/Optimisation/1-OpimalTargeting/Kent_OptimalTargeting")




##---------------------------------------------##
#### 1.4 Visualise outputs from optimization ####
##---------------------------------------------##

## read in the outputs csv file from the optimisaiton function
OptOut <- read_csv("CleanData/Optimisation/1-OpimalTargeting/Kent_OptimalTargeting.csv")

## calculate the change in wader abundance for each field
OptOut <- OptOut |> 
          mutate(SnipeChange = SnipeNew-BaseSnipe,
                 LapwingChange = LapwingNew-BaseLapwing,
                 RedshankChange = RedshankNew-BaseRedshank) |> 
          select(ParcRef, SnipeChange, LapwingChange, RedshankChange )


## Join the output file to the main canvas
CanvUp <- left_join(Canv, OptOut, by = "ParcRef")


## Plot how wader abundance would change overall given improvements at the field level
NK1 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() + 
    geom_sf(mapping = aes(geometry = geometry, fill = (LapwingChange/FieldArea)), colour = NA) +
    scale_fill_viridis_c(name = "Change Lapwing/ha", option = "turbo") +
    ggtitle("North Kent: Change in Lapwing Abundance") +
    theme_light()

ggsave(plot=NK1, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/North Kent Change in Lapwing Abundance.png", units = "in", height = 10, width = 13)
    
    
NK2 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() +  
    geom_sf(mapping = aes(geometry = geometry, fill = (RedshankChange/FieldArea)), colour = NA) +
    scale_fill_viridis_c(name = "Change Redshank/ha", option = "turbo")
    ggtitle("North Kent: Change in Redshank Abundance") +
    theme_light()

ggsave(plot=NK2, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/North Kent Change in Redshank Abundance.png", units = "in", height = 10, width = 13)






##----------------------------##
##----------------------------##
#### 2 *Somerset Scenarios* ####
##----------------------------##
##----------------------------##

##----------------------##
#### 2.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.csv")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = ParcArea/10000, # match RF variable format
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)

## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#            theme_minimal()



## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_Som <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_Som) # free up space

## Read in wader models
SniMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/SnipeRF.rds")
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF.rds")




##---------------------------------##
#### 2.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
# TEST <- Canv[PrIndex, ] |> select(SniMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()

## Predict abundance for Snipe
set.seed(1212)
Canv$SnipeAbund[PrIndex] <- (predict.rfsrc(object = SniMod, 
                             newdata = (Canv[PrIndex, ] |> select(SniMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$SnipeAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------------##
#### 2.3 Run optimization function ####
##-----------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# NewCat= "Reserve"
# Plus = F
# Outpath = "CleanData/Optimisation/1-OpimalTargeting/Kent_OptimalTargeting"
# FenceParcels = FALSE

## Run function for Somerset
OptimiseTarg(Canvas=Canv,
             SniModel=SniMod, LapModel=LapMod, RedModel=NULL,
             LandCov=LC_Som, NewCat= "Reserve", Plus = F, FenceParcels = FALSE, 
             LoopStop = 2000, LoopStart = 1001,
             Outpath = "CleanData/Optimisation/1-OpimalTargeting/Som_OptimalTargeting")




##---------------------------------------------##
#### 2.4 Visualise outputs from optimization ####
##---------------------------------------------##

## read in the outputs csv file from the optimisaiton function
OptOut <- read_csv("CleanData/Optimisation/1-OpimalTargeting/Som_OptimalTargeting.csv")

## calculate the change in wader abundance for each field
OptOut <- OptOut |> 
          mutate(SnipeChange = SnipeNew-BaseSnipe,
                 LapwingChange = LapwingNew-BaseLapwing,
                 RedshankChange = RedshankNew-BaseRedshank) |> 
          select(ParcRef, SnipeChange, LapwingChange, RedshankChange )


## Join the output file to the main canvas
CanvUp <- left_join(Canv, OptOut, by = "ParcRef")


## Plot how wader abundance would change overall given improvements at the field level
Som1 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() + 
    geom_sf(mapping = aes(geometry = geometry, fill = LapwingChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Somerset: Change in Lapwing Abundance") +
    theme_light()

ggsave(plot=Som1, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Somerset Change in Lapwing Abundance.png", units = "in", height = 10, width = 13)
    
    
Som2 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() +  
    geom_sf(mapping = aes(geometry = geometry, fill = SnipeChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Somerset: Change in Snipe Abundance") +
    theme_light()

ggsave(plot=Som2, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Somerset Change in Snipe Abundance.png", units = "in", height = 10, width = 13)







##-----------------------------------##
##-----------------------------------##
#### 3. *Norfolk Broads Scenarios* ####
##-----------------------------------##
##-----------------------------------##

##----------------------##
#### 3.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.csv")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = ParcArea/10000, # match RF variable format
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)

## Visualize canvas categories
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
           theme_minimal()



## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF.rds")




##---------------------------------##
#### 3.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  
## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------------##
#### 3.3 Run optimization function ####
##-----------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# NewCat= "Reserve"
# Plus = F
# Outpath = "CleanData/Optimisation/1-OpimalTargeting/Kent_OptimalTargeting"
# FenceParcels = FALSE

## Run function for North Kent
OptimiseTarg(Canvas=Canv,
             SniModel=NULL, LapModel=LapMod, RedModel=RedMod,
             LandCov=LC_crop, NewCat= "Reserve", Plus = F, FenceParcels = FALSE,
             Outpath = "CleanData/Optimisation/1-OpimalTargeting/Broads_OptimalTargeting")




##---------------------------------------------##
#### 3.4 Visualise outputs from optimization ####
##---------------------------------------------##

## read in the outputs csv file from the optimisaiton function
OptOut <- read_csv("CleanData/Optimisation/1-OpimalTargeting/Broads_OptimalTargeting.csv")

## calculate the change in wader abundance for each field
OptOut <- OptOut |> 
          mutate(SnipeChange = SnipeNew-BaseSnipe,
                 LapwingChange = LapwingNew-BaseLapwing,
                 RedshankChange = RedshankNew-BaseRedshank) |> 
          select(ParcRef, SnipeChange, LapwingChange, RedshankChange )


## Join the output file to the main canvas
CanvUp <- left_join(Canv, OptOut, by = "ParcRef")


## Plot how wader abundance would change overall given improvements at the field level
Br1 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() + 
    geom_sf(mapping = aes(geometry = geometry, fill = LapwingChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Broads: Change in Lapwing Abundance") +
    theme_light()

ggsave(plot=Br1, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Broads Change in Lapwing Abundance.png", units = "in", height = 10, width = 13)
    
    
Br2 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() +  
    geom_sf(mapping = aes(geometry = geometry, fill = RedshankChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Broads: Change in Redshank Abundance") +
    theme_light()

ggsave(plot=Br2, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Broads Change in Redshank Abundance.png", units = "in", height = 10, width = 13)






##--------------------------------##
##--------------------------------##
#### 4. *Essex Coast Scenarios* ####
##--------------------------------##
##--------------------------------##

##----------------------##
#### 4.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.csv")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = ParcArea/10000, # match RF variable format
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)

## Visualize canvas categories
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
           theme_minimal()



## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF.rds")




##---------------------------------##
#### 4.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------------##
#### 4.3 Run optimization function ####
##-----------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# NewCat= "Reserve"
# Plus = F
# Outpath = "CleanData/Optimisation/1-OpimalTargeting/Kent_OptimalTargeting"
# FenceParcels = FALSE

## Run function for North Kent
OptimiseTarg(Canvas=Canv,
             SniModel=NULL, LapModel=LapMod, RedModel=RedMod,
             LandCov=LC_crop, NewCat= "Reserve", Plus = F, FenceParcels = FALSE,
             Outpath = "CleanData/Optimisation/1-OpimalTargeting/Essex_OptimalTargeting")




##---------------------------------------------##
#### 4.4 Visualise outputs from optimization ####
##---------------------------------------------##

## read in the outputs csv file from the optimisaiton function
OptOut <- read_csv("CleanData/Optimisation/1-OpimalTargeting/Essex_OptimalTargeting.csv")

## calculate the change in wader abundance for each field
OptOut <- OptOut |> 
          mutate(SnipeChange = SnipeNew-BaseSnipe,
                 LapwingChange = LapwingNew-BaseLapwing,
                 RedshankChange = RedshankNew-BaseRedshank) |> 
          select(ParcRef, SnipeChange, LapwingChange, RedshankChange )


## Join the output file to the main canvas
CanvUp <- left_join(Canv, OptOut, by = "ParcRef")


## Plot how wader abundance would change overall given improvements at the field level
Ess1 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() + 
    geom_sf(mapping = aes(geometry = geometry, fill = LapwingChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Essex: Change in Lapwing Abundance") +
    theme_light()

ggsave(plot=Ess1, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Essex Change in Lapwing Abundance.png", units = "in", height = 10, width = 13)
    
    
Ess2 <- CanvUp |> filter(Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve")) |> 
    ggplot() +  
    geom_sf(mapping = aes(geometry = geometry, fill = RedshankChange), colour = NA) +
    scale_fill_viridis_c()
    ggtitle("Essex: Change in Redshank Abundance") +
    theme_light()

ggsave(plot=Ess2, filename= "CleanData/Optimisation/1-OpimalTargeting/Plots/Essex Change in Redshank Abundance.png", units = "in", height = 10, width = 13)


