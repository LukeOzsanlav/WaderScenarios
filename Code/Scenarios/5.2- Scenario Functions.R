## Load in packages & helper functions
pacman::p_load(tidyverse, # used
               data.table, #??
               sf, # used
               terra, # used
               randomForestSRC, # used
               exactextractr, # used
               tictoc, # used for timing scenario creation function
               ggpubr, # used
               truncnorm) # used

options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R") # additional simple functions for use in the CreateScenario function



## TO DO:
## - Need to work out correction for reserve --> reserve+ creation (current rough calculation)
  




##---------------------------------##
## Start of Small Helper Functions ##
##---------------------------------##

##-----------------------------------------------##

# Define function used for finding neighbours around existing cluster
# Function to get neighbors within 100m buffer of current cluster polygons, it excludes those already in the cluster
# cluster =  the field polygons present in the cluster up until this point in the cluster expansion
# all_polygons = the polygons considered for cluster expansion (small subset of the entire canvas to speed up processing)
# clustbuff = the distance in meters around a cluster to consider for expansion, if no polygons within this distance from `cluster` then no neighbors found
get_neighbors <- function(cluster, all_polygons, clustbuff) {

   ## Create buffer around cluster
   cluster_buffer <- st_buffer(st_union(cluster), dist = clustbuff)  # 100 meters buffer
   ## Extract all fields within this buffer
   intersects <- st_intersects(all_polygons, cluster_buffer, sparse = FALSE)
   potential_neighbors <- all_polygons[apply(intersects, 1, any), ]
  
   # Exclude polygons already in the cluster
   cluster_ids <- cluster$ParcRef  # Assuming each polygon has a unique 'ID' column
   neighbors <- potential_neighbors %>% filter(!ParcRef %in% cluster_ids)
   # plot(neighbors["FieldArea"]) # check
  
  return(neighbors) # return the potential neighbors

}



##-----------------------------------------------##


## Calculate the costs reserve management for a land parcel based off Malcolm's economies of scale project
## This takes into account the size of the land parcel and the total size of the reserve to which it belongs
## TotalArea = total area of the site in ha
## ParcelArea = size of the specific land parcel in ha
## InflationAdjust = Adjust the final cost for the field by inflation from Malcolm's 2011/12 estimate

Manage_Costs <- function(TotalArea, ParcelArea, InflationAdjust=TRUE){
  
  ## ha to m2 conversion factor
  ha_to_m = 10000
  
  ## change the total area so they do not go beyond these values,
  ## Malcolm's regression does not go beyond these values
  TotalArea <- ifelse(TotalArea < 50, 50, 
                      ifelse(TotalArea > 1500, 1500, TotalArea))
  
  ## calculate the log of the cost per m2 for the area provided
  LogCost = 4.04 + (0.18 + 0.34)*(log(TotalArea*ha_to_m))
  
  ## Calculate the cost per ha
  OutCost = (exp(LogCost)/(TotalArea*10000))*10000
  
  ## Calculate the total cost for the land parcel
  OutCost = OutCost*ParcelArea
  
  ## Adjust this costs for inflation
  ## This converts from 2011/12 costs to 2023/24 costs 
  if(InflationAdjust==TRUE){OutCost = OutCost*(108.249/82.346)}

  return(OutCost)
  
}



##-----------------------------------------------##



## Calculate the costs per ha of reserve management based off Malcolm's economies of scale project
## ParcelArea = size of the specific land parcel in ha
## Cost_ha = cost of reserve creation per ha
## Reversion_ha = 
## StartCategory =
## StartHab =

CrCost <- read_csv("RawData/AES Costings/Reserve_creation_costs.csv")

Create_Costs <- function(ParcelArea, StartCategory, StartHab = "grass",
                         CreateCost = CrCost$Cost[CrCost$Element =="Creation costs"],
                         ReversionCost = CrCost$Cost[CrCost$Element =="Arable reversion"],
                         Correction = CrCost$Cost[CrCost$Element =="Reserve correction"]){
  
  ## ha to m2 conversion factor
  ha_to_m = 10000
  
  ## calculate the log of the total cost for the area provided
  TotCost = ParcelArea*CreateCost
  
  ## Add the values in for this at a later date
  TotCost = ifelse(StartHab == "arable",  TotCost+(ParcelArea*ReversionCost), TotCost)
  
  ## If the field was already a reserve then cange the creation costs to 0
  TotCost = ifelse(StartCategory == "Reserve", TotCost*Correction, TotCost)
  
  
  return(TotCost)
  
}



##-----------------------------------------------##



## Calculate the costs per ha of reserve management based off Malcolm's economies of scale project
## ParcelArea = size of the specific land parcel in ha
## Cost_ha = cost of reserve creation per ha

CrCost <- read_csv("RawData/AES Costings/Reserve_creation_costs.csv")

SimpFenceCost <- function(ParcelArea, Cost_m=CrCost$Cost[CrCost$Element =="Fencing"]){

  ## ha to m2 conversion factor
  ha_to_m = 10000
  
  ## calculate the perimeter of the area if it were a a rectangle with a ratio of length to width of 1.5
  w = sqrt((ParcelArea*ha_to_m/1.5))
  TotPerim = 5*w
  
  
  ## Multiply the fence perimeter by the cost per m
  TotCost = TotPerim*Cost_m
  
  return(TotCost)
  
}



##-----------------------------------------------##



## Given an annotate canvas for a priority landscape work out the total opportunity area for all scenarios
## START of AreaForScenarios function
AreaForScenarios <- function(Canvas, ## the canvas for the priority landscape, created in `4 - Annotate scenario canvas.R`
                             ScenList ## list of different scenario types that are going to by run
                             ){

  ## Create a data frame to put all different opportunity area sizes into
  ScenAreas <- data.frame(Landscape = Canvas$Landscape[1:nrow(ScenList)],
                          Strategy = NA,
                          OppCat= NA,
                          NewCat=  NA,
                          Plus = NA,
                          OppArea = NA)

  ## run through each row in the data frame of scenarios
  for(j in 1:nrow(ScenList)){

    ##Join together opportunity and replacement categories for CreateScenario() function
    OppCatChoice = c(paste0(if(is.na(ScenList$OppCat1[j])==F){paste0(ScenList$OppCat1[j])}),
                     paste0(if(is.na(ScenList$OppCat2[j])==F){paste0(ScenList$OppCat2[j])}),
                     paste0(if(is.na(ScenList$OppCat3[j])==F){paste0(ScenList$OppCat3[j])}),
                     paste0(if(is.na(ScenList$OppCat4[j])==F){paste0(ScenList$OppCat4[j])}))

    NewCatChoice = c(paste0(if(is.na(ScenList$NewCat1[j])==F){paste0(ScenList$NewCat1[j])}),
                     paste0(if(is.na(ScenList$NewCat2[j])==F){paste0(ScenList$NewCat2[j])}))

    Strategy = ScenList$Strategy[j]
    Plus = ScenList$Plus[j]

    ## Get the index of rows representing the opportunity area for this scenario
    ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
    if(Strategy %in% c("Big", "More") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice) }
    if(Strategy %in% c("Better") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice) }

    ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders
    ## (i.e. AES fields with waders will not be targeted and those without will be)
    if(Strategy %in% c("Big", "More") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice & !(Canvas$Category %in% NewCatChoice & Canvas$Tot_abund > 0)) }
    if(Strategy %in% c("Better") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice & !(Canvas$Category %in% NewCatChoice & Canvas$Tot_abund > 0)) }

    ## Calculate the total field Area for my opportunity Area
    ScenAreas$OppArea[j] <- sum(Canvas[IndOpp,]$FieldArea)*ScenList$PropOppArea[j]

    ScenAreas$Strategy[j] <- Strategy
    ScenAreas$OppCat[j] <- paste0(paste0(if(is.na(ScenList$OppCat1[j])==F){paste0(ScenList$OppCat1[j], " ")}),
                                  paste0(if(is.na(ScenList$OppCat2[j])==F){paste0(ScenList$OppCat2[j], " ")}),
                                  paste0(if(is.na(ScenList$OppCat3[j])==F){paste0(ScenList$OppCat3[j], " ")}),
                                  paste0(if(is.na(ScenList$OppCat4[j])==F){paste0(ScenList$OppCat4[j])}))
    ScenAreas$NewCat[j] <- paste0(paste0(if(is.na(ScenList$NewCat1[j])==F){paste0(ScenList$NewCat1[j], " ")}),
                                  paste0(if(is.na(ScenList$NewCat2[j])==F){paste0(ScenList$NewCat2[j], " ")}))
    ScenAreas$Plus[j] <- Plus

  } ## End of j loop

  return(ScenAreas)

}




##-----------------------------------------------##





##-------------------------------------##
## FUNC START: CreateScenario Function ##
##-------------------------------------##

## Defining function arguments
# Canvas= landscape canvas as a shape file variables need to match those in the RF model
# ScenType = What type of scenario should by run, can only pick one at a time. options are "random", "cluster", "stakeholder", "stakeholder"
# SniModel= RF model for predicting Snipe abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LapModel= RF model for predicting Lapwing abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# RedModel= RF model for predicting Redshank abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LandCov= Landcover data raster covering same areas as Canvas. This is used to create as raster out of Canvas to calculate landscape variables
# N_sets= number of sets to split up the total opportunity area into, each set is run and then wader abundance calculated
# Additive = should additive sampling be carried out (TRUE/FALSE). Additive sampling essentially resets between rounds of N_sets and if Additive=FALSE management is sequentially added to the canvas
# Unfenced = should all parcel created be un-fenced (FALSE- the default). Or should they be fenced at roughly the same proportion as the category being drawn from in the landscape (TRUE)
# OppCat= Field parcel categories that are defined as the opportunity area (e.g. c("Grass Opp"), or c("Arable Opp"))
# NewCat= Field parcel categories that are used to `improve` the chosen opportunity fields to (can use c("AES Only"), or c("Reserve") or both together)
# Plus= (TRUE/FALSE) should habitat plus improvement be undertaken? Where opportunity fields are improved to the standard of fields with waders for the chosen management (NewCat)
# Strategy= Which strategy do you want to use to guide random sampling, options are "Better", "Big", "More"
# PropOppArea= This variable gives the proportion of the opportunity area that is used for a given strategy. Since Bigger and More use the same opportunity area then it has to be divided up. 
# OppAreaUsed = This is the total area in hectares of the opportunity area that can be altered in any scenario. For multiplicative samples this is the total area when adding up all loops of N_sets. For additive sampling provide a vector of length N_sets. Each element in the vector is an area in hectares that should be altered in each round of sampling (allows the Area of each sample to vary) 
# ClustMean= if ScenType=="cluster" then this determines the average size of clusters to create (if available land allows)
# ClustSD= if ScenType=="cluster" then this determines the standard deviation of ClustMean as cluster sizes are drawn from a truncated normal distribution (lower bound =0)
# ClustBuf= if ScenType=="cluster" then this determines the distance around a cluster to look for additional fields to add to the cluster up until cluster size is reached
# StakeGroup= if ScenType == "stakeholder" then write down which group's grading to use (either "G1", "G2" or "G3")
# Outpath = Outpath for the log file and canvas
# SaveCanv = (TRUE/FALSE) should the canvas be saved in between each round of sampling (number of rounds = N_sets)

CreateScenario <- function(Canvas, 
                           ScenType, 
                           SniModel, 
                           LapModel,
                           RedModel, 
                           LandCov,
                           N_sets, 
                           Additive,
                           Unfenced=FALSE,
                           OppCat, 
                           NewCat, 
                           Plus, 
                           Strategy, 
                           PropOppArea=0.5, # currently defunct
                           OppAreaUsed, 
                           ClustBuf= NULL,
                           ClustMean= NULL, 
                           ClustSD= NULL, 
                           StakeGroup = NULL,
                           Outpath, 
                           SaveCanv){

  
  ##---------------------------##  
  ##---------------------------##
  #### F 1.1 Scenario set up ####
  ##---------------------------##
  ##---------------------------##
  
  ##-- SET UP CANVAS --##
  
  ## Remove all the fields that are not needed during the scenario modelling
  ## This is to try and increase the speed
  ## Any fields with no opportunity or that were masked are removed unless they are wider wet grassland adjacent to the priority landscape
  Canvas <- filter(Canvas, Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve") | WiderWetGrass == 1)
  # plot(Canvas$geometry)
  
  ## The number segments to split scenario modeling up into
  ## For each segment habitat variables are updated and wader abundance is calculated
  N_sets <- N_sets
  
  
  ##-- CREATE TRACKER DOCUMENT --##
  
  ## Create data frame to record wader abundance at the end of each segment 
  ## Also record the starting abundance in a separate columns
  ScTracker <- data.frame(Landscape = Canvas$Landscape[1],
                          ScenType =  ScenType,
                          Strategy = Strategy,
                          OppCat= paste(OppCat,collapse="&"),
                          NewCat=  paste(NewCat,collapse="&"), 
                          Plus = Plus, 
                          Additive = Additive,
                          Unfenced = Unfenced,
                          ClustMean = ClustMean,
                          ClustActual = NA,
                          Segment = 1:(N_sets+1), 
                          SegmentArea = NA,
                          BaseAESCost = sum(Canvas$AESCostGDP[Canvas$Category=="AES Only"], na.rm=T),
                          NewAESCost = NA,
                          BaseResMaintCost = sum(Canvas$ResMaintGDP[Canvas$Category=="Reserve"], na.rm=T),
                          NewResMaintCost = NA,
                          NewResCreatCost = NA,
                          FencingCost = NA,
                          PurchaseCost = NA,
                          ForegoneCost = NA,
                          Snipe = NA,
                          Snipe_unfenced = NA, 
                          BaseSnipe = ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T)),
                          Lapwing = NA, 
                          Lapwing_unfenced = NA, 
                          BaseLapwing = ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T)),
                          Redshank = NA,
                          Redshank_unfenced = NA, 
                          BaseRedshank = ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T)),
                          SegTime = NA)
  
  ## Set the first row as zero area and baseline abundance
  ScTracker$Snipe[1] <- ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T))
  ScTracker$Lapwing[1] <- ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T))
  ScTracker$Redshank[1] <- ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T))
  ScTracker$SegmentArea[1] <- 0
  
  
  
  ##-- DEFINE OPPORTUNITY AREA --##
  
  ## Get the index of rows representing the opportunity area for this scenario
  ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
  if(Strategy %in% c("Big", "More") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat) }
  if(Strategy %in% c("Better") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat) }
  
  ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders 
  ## (i.e. AES fields with waders will not be targeted and those without will be)
  if(Strategy %in% c("Big", "More") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat & !(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)) }
  if(Strategy %in% c("Better") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat & !(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)) }
  
  ## Calculate the total field Area for my opportunity Area
  # TotArea <- sum(Canvas[IndOpp,]$FieldArea)*PropOppArea
  ## Edited here: now the user defines an area in hectares to alter rather than a percentage of the total area
  TotArea <- sum(OppAreaUsed)
  
  ## Split up that area into segments
  if(Additive==TRUE){AreaSeg <- OppAreaUsed}else{AreaSeg <- rep(TotArea/N_sets, times = N_sets)}


  
  
  ##---------------------------------------##  
  ##---------------------------------------##
  ## Label fields in each scenario segment ##
  ##---------------------------------------##
  ##---------------------------------------##

  ##---------------------------------------------##  
  #### F 1.2 Label fields for cluster strategy ####
  ##---------------------------------------------##  
  
  
  ##-- CLUSTERED FIELD LABELLING --##
  
  if(ScenType =="cluster" | ScenType == "stakeholder2"){
    
    ## Message to update console
    message("Labelling segments with ", ScenType, "...")
    
    ## If using ScenType "stakeholder2" then will need to retain some extra columns along the way, if using cluster just set this to Null
    if(ScenType == "stakeholder2"){StakeCols <- c(paste0("Bigger_", StakeGroup), paste0("Better_", StakeGroup), 
                                                  paste0("More_", StakeGroup), paste0("Arable_", StakeGroup))}else{StakeCols <- NULL}
    
    
    
    ##-- PREPARE SUB-CANVAS  --##
    
    ## filter out the fields that are within the opportunity area for this scenario
    OppFields <- Canvas[IndOpp,] |> select(all_of(c("FieldArea", "ParcRef", "ClustDist", "ReserveGroup", StakeCols))) # could filter the column of OppFields here as later on just need ParcRef & RandSamp
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    
    ## Start a counter that will label each cluster sequentially
    ## And create columns to put counter in and a column that gives the distance to the center of the cluster
    ClusterCount <- 1
    OppFields <- OppFields |> mutate(ClustID=NA,
                                     CentDist=NA)


    
    ##-- CREATE POLYGONS OF ALL RESERVES (ADDITVE=F ONLY) --##
    
    ## Create a combined set of outlines for all RSPB reserves, need this if creating reserve
    ## If Additive = TRUE then create this inside the loop so that it can be reset after every loop
    if(any(NewCat %in% "Reserve") & Additive==FALSE){
      
      ReserveOutlines <- Canvas |> filter(Reserve=="Y") |> group_by(ReserveGroup) |>  
                          summarize(geometry = st_union(geometry)) }
    
    
    
    ##-- CREAST LIST FOR LOOP RECORDING --##
    
    ## Create an empty list where I can put the parcel references that belong to each round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    listClustID <- vector(mode = "list", length = N_sets)
    listReserveGroup <- vector(mode = "list", length = N_sets)
    listPosit <- vector(mode = "list", length = N_sets)
    
    
    
    ##-- START LOOP --##
    
    ## Label opportunity fields for each round of sampling
    for(k in 1:N_sets){ 
      
      ## Message for progress of cluster sampling
      message("Creating set of clusters ", k , " out of ", N_sets)
      

      
      ##-- SELECT SUB-SET OF COLUMNS --##
      
      ## Filter out the fields that are still in the opportunity area
      ## For each iteration of the loop this will gradually decrease the size of the data set I am working with as more fields become part of a cluster
      # OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols)))
      if(Additive==FALSE){OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(all_of(c("FieldArea", "ParcRef", "ReserveGroup", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols)))}
      if(Additive==TRUE){OppFieldsLoop <- OppFields |> select(all_of(c("FieldArea", "ParcRef", "ReserveGroup", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols))) |> 
                                                       mutate(ClustID=NA, CentDist=NA, RandSamp=NA)}
      
      
      
      ##-- CREATE RESERVE POLYGONS --##
      
      ## Create a combined set of outlines for all RSPB reserves, need this if creating reserve
      ## Create inside the loop for Additive = TRUE as we want to reset the canvas at the start of very k for loop
      if(any(NewCat %in% "Reserve") & Additive==TRUE){
      
        ReserveOutlines <- Canvas |> filter(Reserve=="Y") |> group_by(ReserveGroup) |>  
                            summarize(geometry = st_union(geometry)) }
      
      
      
      ##-- SET LOOP CONTROLS AND OUTPUTS --##
      
      ## set the total area cluster for this round of the loop to zero
      ## This only resets outside of all the while loops
      TotClustArea <- 0
      
      ## On first loop set the shortfall to 0 or always return to zero if Additve==T
      if(k==1 | Additive==TRUE){NextGapExtra <- 0}
      
      ## Counter for while loop
      ## And a list to store the cluster in if the loop has to run multiple times
      WhileCount = 1
      ClustStore <- vector(mode = "list")
      
    
      
      ##-- 1st WHILE LOOP START --##
      
      ## While total area of all clusters created in this iteration of the while function is < AreaSeg, keep on creating more clusters
      ## Need to add on NextGapExtra to make sure sampling isn't too much or too little (depending on if previous loops over/undersampled)
      while(TotClustArea < (AreaSeg[k]+NextGapExtra)) {
        
      
      ##-- SET THRESHOLD AREA --##  
    
      # Generate a random max threshold in hectares for current cluster
      threshold_area <- rtruncnorm(n=1, a=0, mean = ClustMean, sd = ClustSD)  # Example: between 10 and 50 hectares
      
      
      
      ##-- CONTROLS TO MODERATE THRESOLD AREA --##

      ## If doing additive sample where the canvas resets after every full sampling round
      ## And creating large cluster (hence > 200) Then just set the target to the value in Area seg as will just create one cluster per round of sampling
      if(Additive==TRUE & ClustMean > 200){threshold_area <- AreaSeg[k]}
      
      #  1. Stops too much land being changed at final step of the loop by lowering threshold_area to how much area is still required to convert
      if(k==N_sets & (TotClustArea+threshold_area) > (AreaSeg[k]+NextGapExtra) & Additive==FALSE){threshold_area <- (AreaSeg[k]+NextGapExtra)-TotClustArea} 
      #  2. If on the first iteration of the while loop the threshold is exceeded with a single cluster then lower the cluster size to the threshold
      if(TotClustArea==0 & threshold_area > (AreaSeg[k]+NextGapExtra)){
          message("Cluster size larger than enitre segment, lower cluster size") # warning to lower cluster size
          threshold_area <- AreaSeg[k]+NextGapExtra
        } 
      #  3. If the cluster will exceeds the sub threshold before the final loop by more then half the mean cluster size then break the loop and start next sampling round
      #     If it does not exceed it by more than half then create the cluster. (Should try make samples more even in size)
      if((TotClustArea+threshold_area) > (AreaSeg[k]+NextGapExtra+(ClustMean*0.49)) & Additive==FALSE & !k==N_sets) break
      #  4. For any iteration in the Adiitive = TRUE loops we do not want the TotClustArea to exceed the AreaSeg
      #     SO if this happens then shift it back down
      if((TotClustArea+threshold_area) > (AreaSeg[k]+NextGapExtra) & Additive==TRUE){threshold_area <- (AreaSeg[k]+NextGapExtra)-TotClustArea} 
    
      
      
      ##-- SAMPLE CLUSTER NUCLEUS AND INTIATE CLUSTER --##
      
      # Select an initial polygon as the initiation point of the cluster ( this will depend upon the strategy chosen)
      # If ScenType =="cluster" weight the random draw by the distance to cluster (more) or inverse distance to cluster (Bigger)
      OppFieldsLoop2 <- OppFieldsLoop |> filter(is.na(RandSamp)) # remove fields already in cluster from earlier iterations of this while loop
      
      if(nrow(OppFieldsLoop2) ==0) break # break out of while loop if there is no more rows left (this applies if PropOppArea=1)
      
      if(Strategy == "Big" & ScenType =="cluster"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$InvClustDist) 
      }
      if(Strategy == "More" & ScenType =="cluster"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$PClustDist) 
      }
      if(Strategy == "Better" & ScenType =="cluster"){ initial_polygon <- sample_n(OppFieldsLoop2, 1) 
      }
      
      # If ScenType == "stakeholder2"  weight the random draw by the stakeholder gradings
      if(Strategy == "Big" & ScenType == "stakeholder2" & !any(OppCat %in% "Arable Opp")){  
        ColIn <- colnames(OppFieldsLoop2) == paste0("Bigger_", StakeGroup)
        initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = (OppFieldsLoop2[,ColIn] |> st_drop_geometry() |> unlist()))
      }
      if(Strategy == "More" & ScenType == "stakeholder2" & !any(OppCat %in% "Arable Opp")){  
        ColIn <- colnames(OppFieldsLoop2) == paste0("More_", StakeGroup)
        initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = (OppFieldsLoop2[,ColIn] |> st_drop_geometry() |> unlist()))
      }
      if(Strategy == "Better" & ScenType == "stakeholder2" & !any(OppCat %in% "Arable Opp")){ 
        ColIn <- colnames(OppFieldsLoop2) == paste0("Better_", StakeGroup)
        initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = (OppFieldsLoop2[,ColIn] |> st_drop_geometry() |> unlist())) 
      }
      if(ScenType == "stakeholder2" & any(OppCat %in% "Arable Opp")){ 
        ColIn <- colnames(OppFieldsLoop2) == paste0("Arable_", StakeGroup)
        initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = (OppFieldsLoop2[,ColIn] |> st_drop_geometry() |> unlist())+0.000000001) # add small amount in case there are some 0 gradings
        }
      
      # Create the initial cluster
      cluster <- initial_polygon
      cluster_area <- sum(cluster$FieldArea)  # Area in hectares
      
      # Create data set of potential fields for this cluster, should speed up sampling
      PotentialRefs <- (OppFieldsLoop |> filter(is.na(RandSamp)) |> st_crop(st_buffer(cluster, dist = 7500)))$ParcRef 
      ClustOppFields <- filter(OppFieldsLoop, ParcRef %in% c(PotentialRefs))
      # plot(ClustOppFields$geometry)
      
      
      
      ##-- 2nd WHILE LOOP START --##

      # Expand the cluster while the cluster size is below the max threshold (threshold_area)
      while(cluster_area < threshold_area) {
        
        
        ##-- FIND NEIGHBOURS --##
        
        ## Retrieve potential neighbors
        neighbors <- get_neighbors(cluster, ClustOppFields, clustbuff = ClustBuf)
        # plot(neighbors["FieldArea"]) # check
    
        # Break if no more neighbors are found (Could add option to try a larger buffer area, would need to add argument to get_neighbors())
        if(nrow(neighbors) == 0) break
        
        
        ##-- SELECT NEIGHBOURS --##
        
        # If all the chosen neighbors don't take the cluster over the threshold then just add them all
        # If threshold exceeded by taking all field then take closest ones sequentially until threshold is exceed
        NeighbourArea <- sum(neighbors$FieldArea) # total area of selected fields
        if(NeighbourArea < (threshold_area-cluster_area)){selected_neighbor <- neighbors
        }else{
          ## Calculate the distance between potential neighbors and cluster center, then order with the closest fields first
          neighbors <- neighbors |> mutate(CentDist = as.numeric(st_distance(initial_polygon, neighbors))) |> 
                                    arrange(CentDist)
          ## Select the closest neighbors sequentially until the total area exceeds the threshold_area for the entire cluster
          indT <- which(cumsum((neighbors$FieldArea)) < (threshold_area-cluster_area)) # rows of fields that go just up to threshold
          if(length(indT)>0){selected_neighbor <- neighbors[c(indT, (max(indT)+1)),]} # select fields chosen above and add one extra field to just take over threshold
          if(length(indT) == 0){selected_neighbor <- neighbors[1,]} # just use first fields if none selected to take it just over threshold
          
        }
        
        
        ##-- UPDATE CLUSTER --##
        
        # Add the selected neighbour(s) to the cluster
        cluster <- rbind(cluster, selected_neighbor)
        # plot(cluster["FieldArea"]) # check
        
        # Update the cluster area
        cluster_area <- sum(cluster$FieldArea)  # Area in hectares
    
        
      } ## 2nd WHILE LOOP END -- stop expanding single cluster

      # Plot the final cluster
      # plot(st_geometry(cluster), border = 'blue', col = alpha('blue', 0.5))

      
      
      ##-- CONTROLS TO ENSURE ONLY COMPLETED CLUSTERS USED --##
      
      ## if the full cluster was not made then at first do not use it
      ## Store cluster that are too small in a list, hoping that a cluster of the right size can be made
      ## If after 10 number of runs there have been no cluster large enough created then just use the largest one created so far
      if(cluster_area < threshold_area & WhileCount <10){
        ## Store cluster if it is too small, add 1 to loop counter and go on to next iteration of while loop
        ClustStore[[WhileCount]] <- cluster |> st_drop_geometry()
        WhileCount = WhileCount+1 
        next 
        } 
      if(cluster_area < threshold_area & WhileCount >=10){
        ## If loop run lots of times then calculate the cluster with the largest area and choose that as the cluster to go forward with
        MAXClust <- which.max(do.call(rbind, lapply(ClustStore, function(i)colSums(i['FieldArea']))))
        cluster <-  ClustStore[[MAXClust]]
      }
      
      
      
      ##-- RECORD CLUSTER PROPERTIES IN `OppFieldsLoop` --##
      
      ## Label the fields in my opportunity canvas that are part of the cluster in this round
      ## Can use the loop counter to label the random sample number
      ## Then also add on a label for each cluster and the distance of each field to the center of the cluster
      cluster <- cluster |>  arrange(ParcRef) # arrange by ParcRef as this is how main Canvas is ordered
      clusterInd <- which(OppFieldsLoop$ParcRef %in% cluster$ParcRef)
      OppFieldsLoop$RandSamp[clusterInd] <- k
      OppFieldsLoop$ClustID[clusterInd] <- ClusterCount; ClusterCount <- ClusterCount+1

      
      
      ##-- ASSIGN RESERVE CLUSTERS TO EXSITING RESERVE OR NEW RESERVE --##
      
      ## Work out if this cluster will expand an existing reserve or not
      ## If it does then give it the same Reserve Group ID as the reserve it is near
      ## If not then just give it a unique code
      if(any(NewCat %in% "Reserve")){
        
        ## Calculate if any existing reserves are near to my newly created reserve
        BuffClust <- OppFieldsLoop[clusterInd,] |> st_buffer(dist=250) |> summarize(geometry = st_union(geometry))
        # plot(BuffClust$geometry)
        
        ## Test for overlap new reserve and existing reserves
        Overlap <- st_overlaps(BuffClust, ReserveOutlines)
        
        ## if no overlap then just give the cluster a unique grouping code that can not be confused with others
        if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) == 0L){
          OppFieldsLoop$ReserveGroup[clusterInd] <- paste("Group", ClusterCount-1, "_", k) }
        
        ## If there is overlap then give the new reserve cluster the same ReserveGroup name as the reserve that it overlaps with
        if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) > 0L){
          
          if(length(Overlap[[1]])==1){ResNumber <- Overlap[[1]]}
          if(length(Overlap[[1]])>1){ResNumber <- Overlap[[1]][1]}

          
          OppFieldsLoop$ReserveGroup[clusterInd] <- paste0(ReserveOutlines$ReserveGroup[ResNumber])

          ## Update the reserves geometry so that the new cluster become part of the reserve it overlaps with
          ReserveOutlines <- rbind(ReserveOutlines,
                                   data.frame(ReserveGroup = paste0(ReserveOutlines$ReserveGroup[ResNumber]),
                                              geometry= st_union(OppFieldsLoop$geometry[clusterInd])) |> st_as_sf())
        }
        
      }
      
      
      
      ##-- 1st WHILE LOOP CONTROLS --##
      
      ## Calculate the total area of all clusters in this round of for loop
      TotClustArea <- sum((OppFieldsLoop |> filter(RandSamp==k))$FieldArea)
      
      ## Now that a cluster has been chosen I can reset the while loop counter for the creation of the next cluster
      WhileCount = 1
      ClustStore <- vector(mode = "list")
      
        
    } ## 1st WHILE LOOP END -- stop creating clusters
      
      
      
    ##-- RECORD CLUSTER PROPOERTIES IN LISTS --##  
      
    ## Now I have created a series of clusters of the suitable total size
    ## I need to now label these using the list i created earlier
    ## Create a RandSamp column to find all the fields in this round of cluster sampling
    clusterAll <- which(OppFields$ParcRef %in% (OppFieldsLoop |> filter(RandSamp==k))$ParcRef)
    OppFields$RandSamp[clusterAll] <- k
  
    ## Assign the parcel refs chosen in this round of sampling to the list
    ## Also assign the ClustID to a list so each ParcRef can be matched to a ClustID number
    listSamples[[k]] <- (OppFieldsLoop |> filter(RandSamp==k))$ParcRef
    listClustID[[k]] <- (OppFieldsLoop |> filter(RandSamp==k))$ClustID
    
    ## When reserve is create I also want to record what reserve the cluster is part of and how far away each field is from the centee of the cluster
    if(any(NewCat %in% "Reserve")){
        
      ## If creating reserve then also assign what ReserveGroup each clster belongs to  
      listReserveGroup[[k]] <- (OppFieldsLoop |> filter(RandSamp==k))$ReserveGroup
      
      ## Create cluster ID column so that i can group the data set
      OppFields$ClustID[clusterAll] <- (OppFieldsLoop |> filter(RandSamp==k))$ClustID
      
      ## If want to create reserve and AES also calculate whether fields are in the center or periphary of clusters and assign these to a list
      FieldPost <- OppFieldsLoop |> 
                    filter(RandSamp==k) |> 
                    group_by(ClustID) |> 
                    mutate(CentDist = as.numeric(st_distance(st_centroid(st_union(geometry)), geometry))) |> 
                    ungroup()
  
      ## Assign distances to a list
      listPosit[[k]] <- FieldPost$CentDist
    
    }

    
    
    ##-- K LOOP CONTROLS --##
    
    ## Calculate the shortfall in FieldArea for this round of for loop, can add this on next time (stops last sample being very large/and carryover shortfall to next iteration of loop)
    if(Additive==FALSE){NextGapExtra <- (AreaSeg[k]*k) - sum((OppFields |> filter(is.na(RandSamp)==F))$FieldArea)}
      
  }# end of loop k 
    
    ## see the total area of each set of fields
    # OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))
    # ggplot(data = OppFields) + geom_sf(mapping = aes(geometry = geometry, fill = CentDist), colour = NA)
    
  ## Select just the unique field ID column and random set number for a join
  if(Additive==FALSE){OppFields <- OppFields |> select(ParcRef, RandSamp, ClustID, CentDist) |> st_drop_geometry()}
    
        
  } # end of cluster strategy labeling

  
  
  
  ##--------------------------------------------##    
  #### F 1.3 Label fields for random strategy ####
  ##--------------------------------------------##    

  ##-- RANDOM FIELD LABELLING --##
  
  if(ScenType == "random"){
    
    ## Message to update console
    message("Labelling segments randomly...")
    
    
    
    ##-- PREPARE TRIMMED CANVAS --##
    
    ## filter out the fields that are within the opportunity area for this scenario
    OppFields <- Canvas[IndOpp,] |> select(c("FieldArea", "ParcRef", "ClustDist")) |> select(all_of(c("FieldArea", "ParcRef", "ClustDist")))
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    
    
    
    ##-- ORDER OPPORTUNITY FIELDS --##
    
    # sample the rows of OppFields, this mixes up the order and allows for random sampling
    # for bigger weight sampling by inverse distance to wader population
    # for more weight sampling by distance to wader population
    set.seed(1212)
    if(Strategy == "Big"){  my_samp <- sample(1:nrow(OppFields), replace = FALSE, 
                                                     prob = OppFields$InvClustDist) }
    if(Strategy == "More"){  my_samp <- sample(1:nrow(OppFields),  replace = FALSE, 
                                                     prob = OppFields$PClustDist) }
    if(Strategy == "Better"){ my_samp <- sample(1:nrow(OppFields), replace = FALSE) }
    
    
    
    ##-- CREATE LIST/RESERVE OUTLINES FOR LOOP --##
    
    ## Create an empty list where I can put the parcel references that belong to eacg round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    listClustID <- vector(mode = "list", length = N_sets)
    listReserveGroup <- vector(mode = "list", length = N_sets)
    listPosit <- vector(mode = "list", length = N_sets)
    
    ## Create a combined set of outlines for all RSPB reserves, need this if creating reserve
    ## If Additive = TRUE then create this inside the loop so that it can be reset after every loop
    if(any(NewCat %in% "Reserve") & Additive==FALSE){
      
      ReserveOutlines <- Canvas |> filter(Reserve=="Y") |> group_by(ReserveGroup) |>  
                          summarize(geometry = st_union(geometry)) }
    
    
    
    ##-- START LABELLING LOOP --##
    
    ## Label opportunity fields for each round of sampling
    for(j in 1:N_sets){
      
      ##-- CREATE RESERVE OUTLINES --##
      
      ## Create a combined set of outlines for all RSPB reserves, need this if creating reserve
      ## If Additive = TRUE then create this inside the loop so that it can be reset after every loop
      ## If Additive = FALSE want to to sequentially add to these reservegroup outlines
      if(any(NewCat %in% "Reserve") & Additive==TRUE){
      
      ReserveOutlines <- Canvas |> filter(Reserve=="Y") |> group_by(ReserveGroup) |>  
                          summarize(geometry = st_union(geometry)) }
      
      ## On first loop set the shortfall to 0 and also reset if for each round of Addtive sampling
      if(j==1 | Additive==TRUE){NextGapExtra <- 0}
      
      
      
      ##-- SELECT FIELDS FOR RANDOM SAMPLE --##
      
      ## return the sampled rows with cumulative sum < condition
      sampRows <- my_samp[which(cumsum(OppFields$FieldArea[my_samp]) < (AreaSeg[j]+NextGapExtra))] 
      
      ## assign the parcel refs in this sample to the list
      listSamples[[j]] <- OppFields$ParcRef[sampRows]
      OppFields$RandSamp[sampRows] <- j # label the sample with loop counter number, used to keep track of shortfall

      
      
      ##-- ASSIGN RANDOM SAMPLES TO CLUSTERS - AES ONLY --##
      
      ## Here create clusters from the randomly sampled fields
      ## This is only to allow me to calcualte an average cluster size across all loop runs. 
      if(any(NewCat %in% "AES Only")){
        
        ## Buffer and then join the fields in this sample
        BuffClust <- OppFields[sampRows,] |> st_buffer(dist=250) |> summarize(geometry = st_union(geometry))
        
        ## split the multi-polygons created above into polygons and add a cluster ID column
        BuffClust_Sep <- st_cast(BuffClust, "POLYGON") |> mutate(ClustID = 1:length(geometry))
        
        ## Join the cluster ID onto the fields
        OppFieldsJoin <- st_join(OppFields[sampRows,], BuffClust_Sep)
        
        ## Assign the ClustID to a list so each parc ref can be matched to a cluster ID  
        ## NOTE: these cluster IDs are not unique between loops  
        listClustID[[j]] <- OppFieldsJoin$ClustID[match(OppFieldsJoin$ParcRef, OppFields$ParcRef[sampRows])] 
        
      }
      
      
      ##-- ASSIGN RANDOM SAMPLES TO CLUSTERS - RESERVES --##
      
      ## Work out if the drawn field form cluster or expand existing reserves
      ## If it does then give it the same Reserve Group ID as the reserve it is near
      ## If not then just give it a unique code of the small cluster if forms a part of (if any)
      if(any(NewCat %in% "Reserve")){
        
        ## Buffer and then join the fields in this sample
        BuffClust <- OppFields[sampRows,] |> st_buffer(dist=250) |> summarize(geometry = st_union(geometry))
        
        ## split the multi-polygons created above into polygons and add a cluster ID column
        BuffClust_Sep <- st_cast(BuffClust, "POLYGON") |> mutate(ClustID = 1:length(geometry))
        
        ## Join the cluster ID onto the fields
        OppFieldsJoin <- st_join(OppFields[sampRows,], BuffClust_Sep)
        
        ## Now test for overlap between clusters created above and existing reserves
        Overlap <- st_overlaps(BuffClust_Sep, ReserveOutlines)
        
        ## start a cluster counter
        ClusterCount <- 1
        
        ## Now loop through each of the cluster that I created above (these are from fields that were randomly sampled)
        ## Basically checking if the cluster overlap with an existing reserve
        ## For additive=FALSE it also tests for overlap with cluster created in previous iterations of the loop
        for(z in 1:length(Overlap)){
          
          ## if no overlap then just give the cluster a unique grouping code that can not be confused with others
          if(is.integer(Overlap[[z]]) && length(Overlap[[z]]) == 0L){
            
            OppFieldsJoin$ReserveGroup[OppFieldsJoin$ClustID==z] <- paste("Group", ClusterCount-1, "_", j)
            ClusterCount <- ClusterCount+1 # add one to the cluster counter
          }
          
          ## If there is overlap then give the new reserve cluster the same ReserveGroup name as the reserve that it overlaps with
          if(is.integer(Overlap[[z]]) && length(Overlap[[z]]) > 0L){
            
            if(length(Overlap[[z]])==1){ResNumber <- Overlap[[z]]}
            if(length(Overlap[[z]])>1){ResNumber <- Overlap[[z]][1]}
  
            OppFieldsJoin$ReserveGroup[OppFieldsJoin$ClustID==z] <- paste0(ReserveOutlines$ReserveGroup[ResNumber])
          }
          
        } # end of z loop

        
        
      ##-- UPDATE RESERVE OUTLINES --##
        
      ## Update the reserves geometry so that the new cluster become part of the reserve it overlaps with
      ## Also add on the new cluster that I created in this iteraiton of the loop
      ## Only do this for Additve=FALSE
      if(Additive==FALSE){
      ReserveOutlines <- OppFieldsJoin |> select(ReserveGroup) |> 
                     rbind(ReserveOutlines) |> 
                     group_by(ReserveGroup) |>  
                     summarize(geometry = st_union(geometry))}
        
      
        
      ##-- ASSIGN CLUSTER DETAILS TO LIST --##
          
      ## Assign the ClustID to a list so each parc ref can be matched to a cluster ID  
      ## NOTE: these cluster IDs are not unique between loops      
      listClustID[[j]] <- OppFieldsJoin$ClustID[match(OppFieldsJoin$ParcRef, OppFields$ParcRef[sampRows])] 
      
      ## Assign what ReserveGroup each cluster belongs to a list
      listReserveGroup[[j]] <- OppFieldsJoin$ReserveGroup[match(OppFieldsJoin$ParcRef, OppFields$ParcRef[sampRows])]
      
      ## Calculate the distance of each field from the centre of its cluster
      FieldPost <- OppFieldsJoin |>  
                      group_by(ClustID) |> 
                      mutate(CentDist = as.numeric(st_distance(st_centroid(st_union(geometry)), geometry))) |> 
                      ungroup() 
      
      ## Assign distances to a list
      listPosit[[j]] <- FieldPost$CentDist[match(FieldPost$ParcRef, OppFields$ParcRef[sampRows])]
        
      }
      
      
      
      ##-- UPDATE MY_SAMP FOR ADDITVE TRUE/FALSE --##
      
      ## Calculate the shortfall in FieldArea, can add this on next time (stops last sample being very large/and carryover shortfall to next iteration of loop)
      if(Additive==FALSE){NextGapExtra <- ((AreaSeg[j]*j) - sum(OppFields$FieldArea[is.na(OppFields$RandSamp)==F]))}
      
      ## update the sample of rows to remove the rows number just sampled, this create none-overlapping sequential samples
      if(Additive==FALSE){my_samp <- my_samp[(length(sampRows)+1):length(my_samp)]}
      
      ## recreate my samp with a different order for additive sample so samples could overlap
      if(Additive==TRUE){
          if(Strategy == "Big"){  my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = OppFields$InvClustDist) }
          if(Strategy == "More"){  my_samp <- sample(1:nrow(OppFields),  replace = FALSE, prob = OppFields$PClustDist) }
          if(Strategy == "Better"){ my_samp <- sample(1:nrow(OppFields), replace = FALSE) }
       }
      
    } # end of sampling loop (j)
    
    
    ## see the total area of each set of fields
    if(Additive==FALSE){OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))} 
    
  } # end of "random" strategy labeling
  
  
  
  
  ##--------------------------------------------------##    
  #### F 1.4 Label fields for stakeholder1 strategy ####
  ##--------------------------------------------------##    
  
  
  ##-- RANDOM STAKEHOLDER FIELD LABELLING --##
  ## **NOT CURRENTLY IN USE**
  
  if(ScenType == "stakeholder1"){
    
    ## Message to update console
    message("Labelling segments with stakeholder gradings")
    
    ## filter out the fields that are within the opportunity area for this scenario
    ## and select the columns that are needed for the sampling
    OppFields <- Canvas[IndOpp,] |> select(c("FieldArea", "ParcRef", "ClustDist", paste0("Bigger_", StakeGroup), paste0("Better_", StakeGroup), 
                                             paste0("More_", StakeGroup), paste0("Arable_", StakeGroup))) 
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                   InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                   PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    
    
    # sample the rows of OppFields, this mixes up the order and allows for random sampling weighted by stakeholder gradings
    if(Strategy == "Big" & !any(OppCat %in% "Arable Opp")){ 
      Gr <- OppFields[,colnames(OppFields)==paste0("Bigger_", StakeGroup)] |> st_drop_geometry()
      my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1]))
      }
    if(Strategy == "More" & !any(OppCat %in% "Arable Opp")){  
      Gr <- OppFields[,colnames(OppFields)==paste0("More_", StakeGroup)] |> st_drop_geometry()
      my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])) 
      }
    if(Strategy == "Better" & !any(OppCat %in% "Arable Opp")){ 
      Gr <- OppFields[,colnames(OppFields)==paste0("Better_", StakeGroup)] |> st_drop_geometry()
      my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])) 
      }
    if(any(OppCat %in% "Arable Opp")){ 
      Gr <- OppFields[,colnames(OppFields)==paste0("Arable_", StakeGroup)] |> st_drop_geometry()
      my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])+0.000000001) 
      }
    
    ## Create an empty list where I can put the parcel references that belong to each round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    
    ## Label opportunity fields for each round of sampling
    for(j in 1:N_sets){
      
      ## On first loop set the shortfall to 0 and also reset if for each round of Addtive sampling
      if(j==1 | Additive==TRUE){NextGapExtra <- 0}
      
      ## return the sampled rows with cumulative sum < condition
      sampRows <- my_samp[which(cumsum(OppFields$FieldArea[my_samp]) < (AreaSeg[j]+NextGapExtra))] 
      
      ## assign the parcel refs in this sample to the list
      listSamples[[j]] <- OppFields$ParcRef[sampRows]
      OppFields$RandSamp[sampRows] <- j # label the sample with loop counter number, used to keep track of shortfall
      
      ## Calculate the shortfall in FieldArea, can add this on next time (stops last sample being very large/and carryover shortfall to next iteration of loop)
      if(Additive==FALSE){NextGapExtra <- ((AreaSeg[j]*j) - sum(OppFields$FieldArea[is.na(OppFields$RandSamp)==F]))}
      
      ## update the sample of rows to remove the rows number just sampled, this create none-overlapping sequential samples
      if(Additive==FALSE){my_samp <- my_samp[(length(sampRows)+1):length(my_samp)]}
      
      ## recreate my samp with a different order for additive sample so samples could overlap, need to do this using gradings again
      if(Additive==TRUE){
          if(Strategy == "Big" & !any(OppCat %in% "Arable Opp")){ 
            Gr <- OppFields[,colnames(OppFields)==paste0("Bigger_", StakeGroup)] |> st_drop_geometry()
            my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1]))
            }
          if(Strategy == "More" & !any(OppCat %in% "Arable Opp")){  
            Gr <- OppFields[,colnames(OppFields)==paste0("More_", StakeGroup)] |> st_drop_geometry()
            my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])) 
            }
          if(Strategy == "Better" & !any(OppCat %in% "Arable Opp")){ 
            Gr <- OppFields[,colnames(OppFields)==paste0("Better_", StakeGroup)] |> st_drop_geometry()
            my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])) 
            }
          if(any(OppCat %in% "Arable Opp")){ 
            Gr <- OppFields[,colnames(OppFields)==paste0("Arable_", StakeGroup)] |> st_drop_geometry()
            my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = unlist(Gr[,1])) 
            }
         }
      } # end of sampling loop (j)
    
    ## see the total area of each set of fields
    if(Additive==FALSE){OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))} 
    
  } # end of "stakeholder" strategy labeling
  
  
  
  
  ##-------------------------------------##
  ##-------------------------------------##
  ## Loop each segment and update canvas ##
  ##-------------------------------------##
  ##-------------------------------------##
  
  
  ##-- CREATE BUFFERED CANVAS --##
  
  ## Pre-compute buffer for all fields so I just have to do it once
  CanvCants <- Canvas |> st_centroid() |> select(geometry)
  Parc500 <- CanvCants |> st_buffer(dist = 500)
  Parc1000 <- CanvCants |> st_buffer(dist = 1000)
  Parc1500 <- CanvCants |> st_buffer(dist = 1500)
  Parc2000 <- CanvCants |> st_buffer(dist = 2000)
  
  
  
  ##-- ADD COLUMNS TO CANVAS --##
  
  ## Add a column to the data set that keeps track of which fields have been altered during the scenario modelling procedure
  Canvas$Upgraded <- NA
  Canvas$CategoryOrig <- Canvas$Category
  
  ## If you want to also calculate the abundance if no areas are fenced then create columns
  ## for the un-fenced abundance
  if(Unfenced==TRUE){  Canvas <- Canvas |> mutate(Fence_CoverageOrig = Fence_Coverage,
                                                  SnipeAbundUn = SnipeAbund,
                                                  RedAbundUn = RedAbund,
                                                  LapAbundUn = LapAbund) }

  
  ##-- CREATE MASTER COPY --##
  
  ## Store a copy of the unaltered canvas here
  ## Can then call on this to reset the canvas at the start of every loop
  MASTERCanv <- Canvas

  
  
  ##-- START MAIN LOOP THROUGH SAMPLES --##
  
  ## For each of my random sample fields update the habitat within the fields and wider landscape water
  ## Then re-calculate the wader abundance variable
  for(i in 1:N_sets){ 
        
    ## send message to console
    message("Running scenario segment part ", i)
    
    ## start timer for segment
    tic()
    
    
    
    
    ##----------------------------------------##        
    #### F 1.5 Update within field habitats ####
    ##----------------------------------------##
    
    ## send message to console
    message("Updating habitat variables")
    
    
    
    ##-- IDENTIFY FIELDS IN THE CURRENT SET --##
            
    ## Get the row numbers of the opportunity fields in the current random sampling group
    ## Update here so that sampled fields are drawn from a list of parcel refs and not from a column
    # SetInd <- which(Canvas$RandSamp==i)
    SetInd <- which(Canvas$ParcRef %in% listSamples[[i]])
    
    ## Now label the fields in this loop number with the loop count to show that they have been updated from their original state
    Canvas$Upgraded[SetInd] <- i
    
    
    
    ##-- CHOOSE HABITAT COLUMNS TO UPDATE --##
    
    ## define which within habitat categories I am going to update
    ## also update the categopry for that row
    ## For Arable also update the TALL_BOUNDARY_PERCENT, normally for grassland fields I leave this the same but in arable it needs to be updates
    HabCats <- c("Category", "GrassOpp", "GRASSLAND_TYPE", "WaterCoverage", "STOCK", "VEG_STRUCTURE", "RUSH_PERCENT", "Fence_Coverage", "AES_options")
    if(any(OppCat %in% "Arable Opp")){HabCats <- c(HabCats, "TALL_BOUNDARY_PERCENT")}

    ## Assign the total area of this sampling set to the data set that tracks progress
    ScTracker$SegmentArea[i+1] <- sum(Canvas$FieldArea[SetInd])
  
    
    
    ##-- CHOOSE FIELDS FOR NEW HABTIAT SAMPLES --##
    
    ## Get the row numbers of fields with the management option that new opportunity fields will become
    ## for plus sampling remove fields that do not have waders in
    if(Plus==F){IndNewHabs <- which(Canvas$Category %in% NewCat)}
    if(Plus==T){IndNewHabs <- which(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)}
    
    
    
    ##-- UPDATE WITHIN FIELD HABITAT --##
    
    ## Sample a set of fields from the chosen management option equal to the number of fields in the random opportunty sample
    ## If the set if chosen management fields is smaller then the opportunity fields that are going to be changed 
    ## then use replace = TRUE or there will not be enough options in the random sample
    if(length(NewCat)==1){
        
        # does sampling need to be done with replacement of not (depends on required length of sample)
        if(length(IndNewHabs)>length(SetInd)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} 
        NewHabs <- sample(x = IndNewHabs, size = length(SetInd), replace = ReplaceSet)
        
        ## Now assign new habitat values to the opportunity fields
        Canvas[SetInd, colnames(Canvas) %in% HabCats] <- Canvas[NewHabs, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
        print(table(Canvas$Category))
        
        ## If creating reserve also alter the reserve column
        if(any(NewCat %in% "Reserve")){ Canvas$Reserve[SetInd] <- "Y" }
        
    }
    
    ## **NOT CURRENTLY IN USE**
    if(length(NewCat)==2){
      
      if(ScenType == "random" | ScenType == "stakeholder1"){
        
      ## How many of each management option can we sample from
      INH1 <- IndNewHabs[Canvas$Category[IndNewHabs]=="Reserve"]
      INH2 <- IndNewHabs[Canvas$Category[IndNewHabs]=="AES Only"]
      
      ## take half the sample required from management strategy 1
      ## take floor to deal with odd numbers
      if(length(INH1) > floor(length(SetInd)/2)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} # does sampling need to be done with replacement?
      NewHabs1 <- sample(x = INH1, size = floor(length(SetInd)/2), replace = ReplaceSet)
      
      ## take half the sample required from management strategy 2
      ## take ceiling to deal with odd numbers (not floor above)
      if(length(INH2) > ceiling(length(SetInd)/2)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} # does sampling need to be done with replacement?
      NewHabs2 <- sample(x = INH2, size = ceiling(length(SetInd)/2), replace = ReplaceSet)
      
      ## Join together two numbers from two strategies
      NewHabs <- c(NewHabs1, NewHabs2)
      
      ## Now assign new habitat values to the opportunity fields
      Canvas[SetInd, colnames(Canvas) %in% HabCats] <- Canvas[NewHabs, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
      ResInd <- which(Canvas$ParcRef %in% listSamples[[i]] & Category == "Reserve")
      Canvas$Reserve[ResInd,] <- "Y" # also update reserve column
      print(table(Canvas$Category))
        
      }
      if(ScenType =="cluster" | ScenType == "stakeholder2"){
        
              
      ## How many of each management option can we sample from
      INHR <- IndNewHabs[Canvas$Category[IndNewHabs]=="Reserve"]
      INHA <- IndNewHabs[Canvas$Category[IndNewHabs]=="AES Only"]
      
      ## Return which rows are central and which one peripheral for this sampling set
      SetCent <- which(Canvas$ParcRef %in% listSamples[[i]][listPosit[[i]]== "Cent"])
      SetPeriph <- which(Canvas$ParcRef %in% listSamples[[i]][listPosit[[i]]== "Periph"])

      ## take half the sample required from management strategy 1
      ## take floor to deal with odd numbers
      if(length(INHR) > length(SetCent)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} # does sampling need to be done with replacement?
      NewHabsC <- sample(x = INHR, size = length(SetCent), replace = ReplaceSet)
      
      ## take half the sample required from management strategy 2
      ## take ceiling to deal with odd numbers (not floor above)
      if(length(INHA) > length(SetPeriph)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} # does sampling need to be done with replacement?
      NewHabsP <- sample(x = INHA, size = length(SetPeriph), replace = ReplaceSet)
      
      ## Now assign new habitat values to the opportunity fields, need to be done differently for cluster option with reserve and AES being created
      Canvas[SetCent, colnames(Canvas) %in% HabCats] <- Canvas[NewHabsC, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
      Canvas$Reserve[SetCent] <- "Y" # also update reserve column
      Canvas[SetPeriph, colnames(Canvas) %in% HabCats] <- Canvas[NewHabsP, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
      print(table(Canvas$Category))
        
      }
    }


    
    ##------------------------------------##    
    #### F 1.6 Update landscape habitat ####
    ##------------------------------------##
  
    ##-- RASTERIZE CANVAS STANDING WATER --##
    
    ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
    ## The base raster is just the UKCEH land cover data set at 25m res
    ## First create index of rows for fields that are suitable
    Suit <- which(Canvas$GrassOpp == 1 | Canvas$WiderWetGrass == 1)
    # Suit <- which(Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp"))
    WaterRast <- rasterize_polygons(Canvas[Suit, ], LandCov, min_coverage = 0.25)
    
    ### Now assign the pixels a water coverage value from the correct field
    values(WaterRast) <- Canvas[Suit,]$WaterCoverage[values(WaterRast)]
    # plot(WaterRast)
    
    
    
    ##-- WHICH FIELDS NEED UPDATEDING --##
    
    ## Work out which fields need to be updated, only a portion of fields will need to be updated as not all fields
    ## will be within 2000m of a field that has had it's water value changed
    UpdateParcels <- st_intersection(Canvas, (Canvas[SetInd,] |> st_buffer(dist = 2050) |> st_union()))
    
    ## label those fields that need to be updated
    ## These are fields within 2000m of any field that have changed category
    ## and those fields that are reserve, AES only or grass
    ## NOTE might need to filter UpdateParcels
    Canvas$UpdateScape <- NA
    Canvas$UpdateScape <- ifelse(Canvas$ParcRef %in% c(UpdateParcels$ParcRef) & Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp"), 1, 0)
    
    ## get row number of fields that need to be updated
    IndUpdate <- which(Canvas$UpdateScape ==1)
    
    
    
    ##-- UPDATE LANDSCAPE STANDING WATER --##
    
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
    
    
    
    ##-- UPDATE LANDSCAPE WET GRASSLAND --##
    
    ## Option to update extent of lowland wet grassland
    if(any(OppCat %in% "Arable Opp")){
      
      ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
      ## The base raster is just the UKCEH land cover data set at 25m res
      ## First create index of rows for fields that are suitable
      Suit <- which(Canvas$GrassOpp == 1 | Canvas$WiderWetGrass == 1)
      # Suit <- which(Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp"))
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
    #### F 1.7 Predict Wader abundance ####
    ##-----------------------------------##
    
    ## send message to console
    message("Updating Wader Abundance")
    
    ## Get the index of the rows that I want to predict wader abundance into
    Index <- which(Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp") & Canvas$UpdateScape ==1)
    
    
    
    ##-- UPDATE ABUNDANCE (PROPORTIONAL FENCED NEW PARCELS) --##
    
    ## Predict wader abundance and then assign it to the correct abundance column
    ## If there is no model object provided then do not predict for that species
    ## Predict for Snipe
    if(is.null(SniModel)==F){
        Pred <- (predict.rfsrc(object = SniModel, 
                              newdata = (Canvas[Index, ] |> select(SniModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                              jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$SnipeAbund[Index]))stop("Prediction error: incorrect number of Snipe fields predicited")
        Canvas$SnipeAbund[Index] <- Pred
        ScTracker$Snipe[i+1] <- sum(Canvas$SnipeAbund, na.rm = T)
        } # number of Snipe in the landscape
      
    
    ## Predict for Lapwing
    if(is.null(LapModel)==F){
        Pred <- (predict.rfsrc(object = LapModel, 
                                newdata = (Canvas[Index, ] |> select(LapModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$LapAbund[Index]))stop("Prediction error: incorrect number of Lapwing fields predicited")
        Canvas$LapAbund[Index] <- Pred
        ScTracker$Lapwing[i+1] <- sum(Canvas$LapAbund, na.rm = T)
        } # number of Lapwing in the landscape 
      
    
    ## Predict for Redshank
    if(is.null(RedModel)==F){
        Pred <- (predict.rfsrc(object = RedModel, 
                                newdata = (Canvas[Index, ] |> select(RedModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$RedAbund[Index]))stop("Prediction error: incorrect number of Redshank fields predicited")
        Canvas$RedAbund[Index] <- Pred
        ScTracker$Redshank[i+1] <- sum(Canvas$RedAbund, na.rm = T)
        } # number of Lapwing in the landscape 
    
    
    
    
    ##-- UPDATE ABUNDANCE (UNFENCED NEW PARCELS) --##
    
    ## If Unfenced==TRUE then also calculate the wader abundance if no fences were used in the landscape
    if(Unfenced==TRUE){
      
      if(is.null(SniModel)==F & any(NewCat %in% "Reserve")){
        Pred <- (predict.rfsrc(object = SniModel, 
                               newdata = (Canvas[Index, ] |> mutate(Fence_Coverage = Fence_CoverageOrig) |>  
                                            select(SniModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame() ), 
                               jitt=FALSE))$predicted
        if(!length(Pred)==length(Canvas$SnipeAbundUn[Index]))stop("Prediction error: incorrect number of Snipe fields predicited")
        Canvas$SnipeAbundUn[Index] <- Pred
        ScTracker$Snipe_unfenced[i+1] <- sum(Canvas$SnipeAbundUn, na.rm = T)
        } # number of Snipe in the landscape
      
    
      ## Predict for Lapwing
      if(is.null(LapModel)==F & any(NewCat %in% "Reserve")){
          Pred <- (predict.rfsrc(object = LapModel, 
                                  newdata = (Canvas[Index, ] |> mutate(Fence_Coverage = Fence_CoverageOrig) |> 
                                               select(LapModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame() ), 
                                  jitt=FALSE))$predicted
          if(!length(Pred)==length(Canvas$LapAbundUn[Index]))stop("Prediction error: incorrect number of Lapwing fields predicited")
          Canvas$LapAbundUn[Index] <- Pred
          ScTracker$Lapwing_unfenced[i+1] <- sum(Canvas$LapAbundUn, na.rm = T)
          } # number of Lapwing in the landscape 
        
      
      ## Predict for Redshank
      if(is.null(RedModel)==F & any(NewCat %in% "Reserve")){
          Pred <- (predict.rfsrc(object = RedModel, 
                                  newdata = (Canvas[Index, ] |> mutate(Fence_Coverage = Fence_CoverageOrig) |> 
                                               select(RedModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame() ), 
                                  jitt=FALSE))$predicted
          if(!length(Pred)==length(Canvas$RedAbundUn[Index]))stop("Prediction error: incorrect number of Redshank fields predicited")
          Canvas$RedAbundUn[Index] <- Pred
          ScTracker$Redshank_unfenced[i+1] <- sum(Canvas$RedAbundUn, na.rm = T)
          } # number of Lapwing in the landscape 

    }

    
    
    ##---------------------------##
    #### F 1.8 Update costings ####
    ##---------------------------##

    ##-- AES ONLY COSTING --##  
     
    ## If creating AES only land parcels then create the costings
    if(any(NewCat %in% "AES Only")){
      
      ## Read in a table with the costing for different AES options
      AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
      SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
      AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")
      
      ## Work out which land parcels should be updated, that are in this sample and have been converted to AES
      ## Need this step in case creating AES only land and reserve
      AESUpd <- which(Canvas$ParcRef %in% listSamples[[i]])
      AES_Check <- which(Canvas$Category[AESUpd] == "AES Only")
      AESUpd <- AESUpd[AES_Check]


      ## Now using the list column of the different AES schemes this function looks through all the AES codes
      ## Then it works out the total payment for those fields based on the AES payment rates and the area of the field
      ## This is only done for fields that are being updates to AES in this round of sampling
      for(j in 1:nrow(AES_Costs)){
        
        ## create columns used to track AES payment 
        Canvas$Add <- NA
        Canvas$AddAmount <- 0
        
        ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
        Canvas$Add[AESUpd] <- lapply(Canvas$AES_options[AESUpd], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})
        
        ## Calculate the cost for that field of the current AES scheme
        ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
        Canvas[AESUpd,] <- Canvas[AESUpd,] %>%
                   mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                .default = 0),
                          AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                          AESCostGDP = AESCostGDP + AddAmount)
        
        ## remove the column I created 
        if(j == nrow(AES_Costs)){Canvas <- Canvas |> select(-c(Add, AddAmount))}
      }
      
      ## Calculate the total cost of all the AES payments
      ScTracker$NewAESCost[i+1] <- sum(Canvas$AESCostGDP[Canvas$Category=="AES Only"], na.rm=T)
      
    }
    
    
    
    ##-- RESERVE COSTING --## 
    
    ## If creating Reserve land parcels then calculate the costings
    if(any(NewCat %in% "Reserve")){
      

      
      ##-- UPDATE RESERVE AREAS --##
      
      ## Calculate the size of all existing reserves (needs to be done as some existing reserves larger then the shapes in the canvas)
      ## Can then add on size of newly created reserves
      OrigAreas <- Canvas |> 
                   filter(is.na(ReserveGroup)==F) |> st_drop_geometry() |> 
                   group_by(ReserveGroup) |> 
                   summarize(OrigArea = max(GroupArea, na.rm=T)) 

      ## For any parcels that has been turned into a reserve in this loop
      ## I want to assign the ReserveGroup ID, this could be a new code OR a reserve name if an existing reserve was expanded
      Canvas$ReserveGroup[SetInd] <- listReserveGroup[[i]]
      ## For new reserve that was deemed to be an extension of an old reserve
      ## Add on the size of the original reserve in the `GroupArea` column
      Canvas$GroupArea[SetInd] <- OrigAreas$OrigArea[match(Canvas$ReserveGroup[SetInd], OrigAreas$ReserveGroup)]
      
      
      ## For any newly created reserves, calculate the size of all the land parcels
      RGarea <- Canvas[SetInd,] |> 
                group_by(ReserveGroup) |> 
                summarize(geometry = st_union(geometry)) |> 
                mutate(AreaAdd = as.numeric(m2_to_ha(st_area(geometry)))) |> 
                st_drop_geometry()
      
      ## Now join the additional areas onto the main canvas
      ## Add the new areas onto the existing reserves size
      ## For completely new resereves only the size of the newly created area is used
      Canvas <- left_join(Canvas, RGarea, by = "ReserveGroup") |> 
                rowwise() %>% 
                mutate(GroupArea = sum(GroupArea, AreaAdd, na.rm=T)) |> 
                select(-AreaAdd)
      
      
      
      ##-- UPDATE MANGEMENT COSTS --##
      
      ## Work out which rows are reserves 
      WhichReserve <- which(Canvas$Category == "Reserve")
      
      ## Work out the costs of ongoing management of the field per hectare
      Canvas$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canvas$GroupArea[WhichReserve], 
                                                       ParcelArea = Canvas$FieldArea[WhichReserve],
                                                       InflationAdjust=T)
      
      ## Calculate the total cost of managing all reserve fields
      ## This includes all the fields that are reserve at this point in the scenario function
      ScTracker$NewResMaintCost[i+1] <- sum(Canvas$ResMaintGDP[Canvas$Category=="Reserve"], na.rm=T)
      
      
      
      ##-- CALCULATE CREATION COSTS --##          
      
      ## Work out the costs of site creation for each parcel that has been created
      ## if the parcels changed to reserve were originally arable then add on an extra creation cost
      Canvas$ResCreatGDP[SetInd] <- Create_Costs(ParcelArea= Canvas$FieldArea[SetInd], 
                                                 StartHab = ifelse(Canvas$ArableOpp[SetInd] == 1, "arable","grass"), 
                                                 StartCategory = Canvas$CategoryOrig[SetInd])
      
      ## This calculate the total cost of doing creation up to this point in the scenario creation
      ## For the additive = FALSE scenarios this will sum up all the creation costs across all loops
      ScTracker$NewResCreatCost[i+1] <- sum(Canvas$ResCreatGDP[Canvas$Category=="Reserve"], na.rm=T)

      
      
      ##-- CALCULATE FENCING COSTS --## 
      
      ## First just select the parcel that have been altered in this round of the i loop
      ## This reduces the size of the data set to speed up processing a bit as well
      Set_i <- Canvas[SetInd,] |> select(c(ParcRef, Fence_Coverage, FieldArea))
      
      ## for the cluster add on their cluster ID and the distance of the each field to the geometric center of the cluster
      Set_i$Cluster_ID <- listClustID[[i]]
      Set_i$ClusterCentDist <-listPosit[[i]]
      
      ## plot to check this worked
      # ggplot(data=Set_i)+geom_sf(mapping = aes(geometry=geometry, fill=Fence_Coverage))
      
      ## For each cluster of fields calculate the total hecterage that the fence needs to cover
      Fence_per_Clust <- Set_i |> 
                         group_by(Cluster_ID) |> st_drop_geometry() |> 
                         arrange(ClusterCentDist) |> 
                         summarise(TotalFence = sum(ifelse(Fence_Coverage=="Y", FieldArea, 0))) 
      
      ## Finally calculate the cost to fence this area if it were a rectange with length = 1.5*width
      FenceCost <- Fence_per_Clust |> 
                    mutate(CostGDP = SimpFenceCost(ParcelArea=TotalFence))
      
      ## For additive=FALSE want to include the costs of all the fencing in previous loops
      if(Additive==T){FinalFenceCost <- sum(FenceCost$CostGDP, na.rm=T)}
      if(Additive==F){
        if(i==1){FinalFenceCost <- sum(FenceCost$CostGDP, na.rm=T)}else{
          FinalFenceCost <- FinalFenceCost + sum(FenceCost$CostGDP, na.rm=T) }
        }
      
      
      ## Now assign the total fencing costs to the tracker document
      ScTracker$FencingCost[i+1] <- FinalFenceCost
      
      
      
      ##-- CALCULATE PURCHASE COSTS --##
      
      ## Calculate the cost of purchasing all of the land
      ## For additive=FALSE this includes the costs of all the land in previous loops
      if(Additive==T){PurchaseCost <- sum(Canvas$PurchaseGDP[SetInd], na.rm = T)}
      if(Additive==F){
        
        ## extracting purchase costs like this means that all pirchase costs from previous iteration of the 
        ## i loop are also included
        AllUpgrade <- which(is.na(Canvas$Upgraded)==F)
        PurchaseCost <- sum(Canvas$PurchaseGDP[AllUpgrade], na.rm = T)
        
        }
      
      ## Assign the purchase cost to the tracker
      ScTracker$PurchaseCost[i+1] <- PurchaseCost

                          
      
      ##-- CALCULATE INCOME FOREGONE COSTS --##
      
      ## Calculate the income foregone for turning all of this land into effectively reserve
      ## For now this presumes for now that the the entire gross margin (minus savings for arable reversion)
      ## is the income foregone, i.e. there is no gross margin for the reserve. 
      if(Additive==T){ ForegoneCost <- sum(Canvas$ForegoneGDP[SetInd], na.rm = T) }
      if(Additive==F){
        
        ## extracting purchase costs like this means that all pirchase costs from previous iteration of the 
        ## i loop are also included
        AllUpgrade <- which(is.na(Canvas$Upgraded)==F)
        ForegoneCost <- sum(Canvas$ForegoneGDP[AllUpgrade], na.rm = T)
        
        }

      ## Assign the costs of income foregone to the tracker
      ScTracker$ForegoneCost[i+1] <- ForegoneCost

    }
    
    
    
    
    ##------------------------------------------##
    #### F 1.9 Save current iteration outputs ####
    ##------------------------------------------##
    
    ##-- SAVE CANVAS AFTER EVERY LOOP --##
    
    ## First calculate the average cluster size for this round of sampling
    # Create a data frame to combine cluster ID and field parcel areas
    DaTa <- data.frame(Area = Canvas$FieldArea[SetInd], Groups = listClustID[[i]])

    # Sum the values by cluster ID and then take the mean total size across groups
    ScTracker$ClustActual[i+1] <- mean(aggregate(Area ~ Groups, DaTa, sum)$Area)
    
    ## Read out shape file of current scenario state
    ## Might want to only save a subset of columns that might be important for plotting
    ## Alternatively could just save the final output as a shapefile
    if(SaveCanv==T){write_sf(Canvas, paste0(Outpath, "Canvas_", Strategy, "_", paste(NewCat,collapse="&"), if(Plus==T){paste0("_plus")}, if(any(OppCat %in% "Arable Opp")){"_Arable"}, i, ".shp"), append=FALSE)}

    ## Update how long it took to run this segment
    toc(log = TRUE, quiet = TRUE)
    ScTracker$SegTime[i+1] <- paste0(tic.log(format = T))
    tic.clearlog()
    write_csv(ScTracker, paste0(Outpath, 
                              "Track_", 
                              if(Additive==TRUE){paste0("Add_")},
                              if(Additive==FALSE){paste0("Mult_")},
                              Strategy, "_", 
                              paste(NewCat,collapse="&"), 
                              if(Plus==T){paste0("_plus")}, 
                              if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
    
    ## If doing the additive testing, i.e. only can look at additive effects
    ## Then reset the canvas at this point, this will mean that the in each round of the i loop we go back to the starting canvas
    ## and do not start from where we left off from on the last loop
    if(Additive==TRUE){Canvas <- MASTERCanv}
    if(Additive==FALSE){Canvas <- Canvas}
    
    
  } # End of segment loop (i)
  
  
  ##-- FINAL SAVE CANVAS AFTER ALL LOOPS --##
    
  ## Add on column for cumulative sum of area
  if(Additive==FALSE){ ScTracker <- mutate(ScTracker, CumArea = cumsum(SegmentArea)) }
  
  ## return the tracker data set
  print(ScTracker)
    
  ## finally write out the file once again
  write_csv(ScTracker, paste0(Outpath, 
                              "Track_", 
                              if(Additive==TRUE){paste0("Add_")},
                              if(Additive==FALSE){paste0("Mult_")},
                              Strategy, "_", 
                              paste(NewCat,collapse="&"), 
                              if(Plus==T){paste0("_plus")}, 
                              if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
  
} ## FUNCTION `CreateScenario` END




##-----------------------------------------------##




inpath="CleanData/Scenarios/5-ScenarioCreation/Kent/"
outpath="CleanData/Scenarios/5-ScenarioCreation/Kent/Plots/"
species="LapRed"

## This function takes the scenario modelling output files from the CreateScenario function and 
## plots the outputs to plot a whole series of plots for a write up
PlotScenario <- function(inpath,
                         outpath,
                         species){
  
  ## Here set the number of years to spread one off costs over
  ## This includes Fencing costs and Reserve creation costs
  Yr <- 15
  
  
  ##--------------------------------##
  ##--------------------------------##
  ##  Plot Multiplicative Scenarios ##
  ##--------------------------------##
  ##--------------------------------##
  
  # Identify all scenario tracker files that were saved
  files <- c(dir(paste0(inpath, "rand/"), pattern  = "Track_Mult", full.names  = T),
             dir(paste0(inpath, "clustlarge/"), pattern  = "Track_Mult", full.names  = T),
             dir(paste0(inpath, "clustsmall/"), pattern  = "Track_Mult", full.names  = T))
  
  ## Read in all these files and bind them together
  AllScn <- files %>%
              map(fread) %>% # read in all files into a list
              bind_rows() %>% # bind all the rows
              mutate(ScenType = ifelse(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T), "clusterlarge", ScenType),
                     ScenType = ifelse(ScenType == "cluster" & ClustMean == min(ClustMean, na.rm = T), "clustersmall", ScenType),
                     ChangeSnipe = Snipe-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeCosts = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                     PlotCat = paste0(NewCat, ifelse(Plus==T, "+", "-"), " ", ifelse(OppCat== "Arable Opp", "Arable", ""))) 
  
  
  ## Calculate the cost per pairs depending on whether it is Snipe or Lapwing/Redshank
  if(species == "LapRed"){ AllScn <- AllScn |> mutate(PairCost = ChangeWaders/ChangeCosts) }
  if(species == "Snipe"){ AllScn <- AllScn |> mutate(ChangeWaders = ChangeSnipe,
                                                     PairCost = ChangeSnipe/ChangeCosts) }
  
  
  
  ##------------------------------------------------------##
  ##  Multiplicative: Abundance Change vs Cumulative Area ##
  ##------------------------------------------------------##
  
  ## Plot the scenario outputs, facet by strategy and plot category
  Mult1 <- ggplot(data = AllScn, mapping = aes(x=CumArea, y=ChangeWaders, colour = ScenType, group = ScenType)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(name = "Scenario Type",   # Change legend title
                        labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small",
                                   "random" ="Random"),  # Change legend labels
                        values = c("clusterlarge"="#540d6e", "clustersmall"="#ee4266",
                                   "random"="#ffd23f")) + # Change legend colors
  xlab("Cumulative Area Converted (ha)") + ylab("Change in Wader Abundance (pairs)") +
  labs(colour = "Scenario Type") +
  facet_grid(PlotCat~Strategy, scales = "free_x")+
  theme_light() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "top")
  
  ## save the plot
  ggsave(plot=Mult1, filename= paste0(outpath, "Mult_Pairs_vs_CumArea.png"), units = "in", height = 10, width = 13)
  
  
  
  
  ##----------------------------------------------------------##
  ##  Multiplicative: Cost per Pair Change vs Cumulative Area ##
  ##----------------------------------------------------------##
  
  ## Plot the scenario outputs, facet by strategy and plot category
  Mult2 <- ggplot(data = AllScnSub, mapping = aes(x=CumArea, y=PairCost*100000, colour = ScenType, group = ScenType)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(name = "Scenario Type",   # Change legend title
                        labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small",
                                   "random" ="Random"),  # Change legend labels
                        values = c("clusterlarge"="#540d6e", "clustersmall"="#ee4266",
                                   "random"="#ffd23f")) + # Change legend colors
  xlab("Cumulative Area Converted (ha)") + ylab("Breeding Wader Pairs/ £100,000") +
  labs(colour = "Scenario Type") +
  facet_grid(PlotCat~Strategy, scales = "free_x") +
  theme_light() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "top")
  
  ## save the plot
  ggsave(plot=Mult2, filename= paste0(outpath, "Mult_PairCost_vs_CumArea.png"), units = "in", height = 10, width = 13)
  
  
  
  
  ##--------------------------##
  ##--------------------------##
  ##  Plot Additive Scenarios ##
  ##--------------------------##
  ##--------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  files <- c(dir(paste0(inpath, "rand/"), pattern  = "Track_Add", full.names  = T),
             dir(paste0(inpath, "clustlarge/"), pattern  = "Track_Add", full.names  = T),
             dir(paste0(inpath, "clustsmall/"), pattern  = "Track_Add", full.names  = T))
  
  
  ## Read in all these files and bind them together
  AllScn2 <- files %>%
              map(fread) %>% # read in all files into a list
              bind_rows() %>% # bind all the rows
              mutate(ScenType = ifelse(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T), "clusterlarge", ScenType),
                     ScenType = ifelse(ScenType == "cluster" & ClustMean == min(ClustMean, na.rm = T), "clustersmall", ScenType),
                     ChangeSnipe = Snipe-BaseSnipe,
                     ChangeSnipenoF = Snipe_unfenced-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeWadersnoF = (Lapwing_unfenced+Redshank_unfenced) - (BaseLapwing+BaseRedshank),
                     ChangeCostsFG = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                     ChangeCostsPU = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+(PurchaseCost/Yr))-BaseResMaintCost, 0)),
                     ChangeCostsNoF = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                     PlotCat = paste0(NewCat, ifelse(Plus==T, "+", "-"), " ", ifelse(OppCat== "Arable Opp", "Arable", ""))) |>
              filter(!SegmentArea == 0) 
  
  
  ## Calculate the cost per pairs depending on whether it is Snipe or Lapwing/Redshank
  if(species == "LapRed"){ AllScn2 <- AllScn2 |> mutate(PairCost = ChangeWaders/ChangeCostsFG,
                                                       PairCostPU = ChangeWaders/ChangeCostsPU,
                                                       PairCostNoF = ChangeWadersnoF/ChangeCostsNoF,
                                                       Waders_100Ha = (ChangeWaders/SegmentArea)*100) }
  
  if(species == "Snipe"){ AllScn2 <- AllScn2 |> mutate(PairCost = ChangeSnipe/ChangeCostsFG,
                                                      PairCostPU = ChangeSnipe/ChangeCostsPU,
                                                      PairCostNoF = ChangeSnipenoF/ChangeCostsNoF,
                                                      Waders_100Ha = (ChangeSnipe/SegmentArea)*100) }
  
  
  ## Set a general theme for the bar plots
  BarPlotTheme <- theme_light() +
                  theme(panel.grid.minor = element_blank(),
                        strip.text = element_text(size = 13, face = "bold"),
                        axis.title = element_text(size = 15),
                        axis.text.y = element_text(size = 13),
                        axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5),
                        legend.title = element_text(size = 13, face = "bold"),
                        legend.text = element_text(size = 13),
                        legend.position = "top")
  
  
  
  ##-----------------------------------------##
  ##  Additive: Pairs per 100 ha vs Category ##
  ##-----------------------------------------##
  
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSum <- AllScn2 |>
    mutate(PlotCatFull = paste0(Strategy, " & ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Waders_100Ha = mean(Waders_100Ha))
  
  
  ## Create Bar plot
  ggplot(ScenSum, aes(x= PlotCatFull, y= Waders_100Ha)) +
    geom_col(aes(fill = ScenType), width=0.55, position=position_dodge(0.55)) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Scenario Type",   # Change legend title
                        labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                        values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
    ylab("Change in Breeding Wader Pairs/ 100 ha") +
    xlab("Scenario Category") +
    labs(fill = "Targeting Strategy") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pairha_vs_Category.png"), units = "in", height = 9, width = 11)
  
  
  
  
  ##-------------------------------------------##
  ##  Additive: Pairs per £100,000 vs Category ##
  ##-------------------------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSumCost <- AllScn2 |>
    mutate(PlotCatFull = paste0(Strategy, " & ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost))
  
  
  ## Create Bar plot
  ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost*100000)) +
    geom_col(aes(fill = ScenType), width=0.55, position=position_dodge(0.55)) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Scenario Type",   # Change legend title
                        labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                        values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
    ylab("Breeding Wader Pairs/ £100,000") +
    xlab("Scenario Category") +
    labs(fill = "Targeting Strategy") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_Category.png"), units = "in", height = 9, width = 12)
  
  
  
  
  ##------------------------------------------------------##
  ##  Additive: Pairs per 100 ha vs Summarized Categories ##
  ##------------------------------------------------------##
  
  ## Summarise the data so there are less rows
  ## Essentially the data is summarized to by taking the average across the different targeting strategies
  ScenSum2 <- ScenSum |>
            filter(!PlotCatFull %in% c("Big & Reserve from Arable", "More & Reserve from Arable")) |> 
            group_by(Strategy, Plus, NewCat) |>
            summarise(Waders_100Ha = mean(Waders_100Ha)) |>
            mutate(PlotCatOther = paste0(NewCat, " ", ifelse(Plus==T, "+", "-")),
                   PlotCatFull= paste0(Strategy, " & ", PlotCatOther),
                   PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling"))
  
  ## Create bar plot
  ## Much simpler than above and mainly focusses on showing the differences between AES and reserve for the different Lawton Principles
  ggplot(ScenSum2, aes(x= Strategy, y= Waders_100Ha)) +
    geom_col(aes(fill = NewCat), width=0.55, position=position_dodge(0.55)) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Scenario Type",   # Change legend title
                      labels = c("AES Only"= "AES Only", "Reserve"= "Nature Reserve"),  # Change legend labels
                      values = c("AES Only"="#F076A5", "Reserve"= "#76A5F0")) +
    ylab("Change in Breeding Wader Pairs/ 100 ha") +
    xlab("Lawton Principle") +
    labs(fill = "Management Used") +
    BarPlotTheme
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pairha_vs_SumCategory.png"), units = "in", height = 9, width = 11)
  
  
  
  
  ##--------------------------------------------------------##
  ##  Additive: Pairs per £100,000 vs Summarized Categories ##
  ##--------------------------------------------------------##
  
  ## Summarise the data so there are less rows
  ## Essentially the data is summarized to by taking the average across the different targeting strategies
  ScenSumCost2 <- ScenSumCost |>
              filter(!PlotCatFull %in% c("Big & Reserve from Arable", "More & Reserve from Arable")) |> 
            group_by(Strategy, Plus, NewCat) |>
            summarise(Ave_PairCost = mean(Ave_PairCost)) |>
            mutate(PlotCatOther = paste0(NewCat, " ", ifelse(Plus==T, "+", "-")),
                   PlotCatFull= paste0(Strategy, " & ", PlotCatOther),
                   PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling"))
  
  ## Create bar plot
  ## Much simpler than above and mainly focuses on showing the differences between AES and reserve for the different Lawton Principles
  ggplot(ScenSumCost2, aes(x= Strategy, y= Ave_PairCost*100000)) +
    geom_col(aes(fill = NewCat), width=0.55, position=position_dodge(0.55)) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Scenario Type",   # Change legend title
                      labels = c("AES Only"= "AES Only", "Reserve"= "Nature Reserve"),  # Change legend labels
                      values = c("AES Only"="#F076A5", "Reserve"= "#76A5F0")) +
    ylab("Breeding Wader Pairs/ £100,000") +
    xlab("Lawton Principle") +
    labs(fill = "Management Used") +
    BarPlotTheme
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_SumCategory.png"), units = "in", height = 9, width = 11)
  
  
  
  
  ##--------------------------------------##
  ##--------------------------------------##
  ##  Plot Scenarios with/without Fencing ##
  ##--------------------------------------##
  ##--------------------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ## Calculate the change in pairs per unit cost for un-fenced reserve
  ScenSumCost_NoFence <- AllScn2 |>
  filter(!OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCostNoF)) |> 
    ungroup() |> 
    mutate(Fencing = "Reserve Unfenced\n(grassland)") |> 
    filter(NewCat == "Reserve")
  
  ## Calculate the change in pairs per unit cost for fenced reserve
  ScenSumCost_Fence <- AllScn2 |>
    filter(!OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Fencing = "Reserve Fenced\n(grassland)") |> 
    filter(NewCat == "Reserve")
  
  ## Calculate the change in pairs per unit cost for AES only land
  ScenSumCost_AES <- AllScn2 |>
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Fencing = "AES Only") |> 
    filter(NewCat == "AES Only")
  
  
  
  ## Calculate the change in pairs per unit cost for un-fenced reserve from arable land
  ScenSumCost_NoFArable <- AllScn2 |>
    filter(OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCostNoF)) |> 
    ungroup() |> 
    mutate(Fencing = "Reserve Unfenced\n(arable)") |> 
    filter(NewCat == "Reserve")
  
  ## Calculate the change in pairs per unit cost for fenced reserve from arable land
  ScenSumCost_FenceArable <- AllScn2 |>
    filter(OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Fencing = "Reserve Fenced\n(arable)") |> 
    filter(NewCat == "Reserve")
  
  
  ## Combine all data sets for plotting
  FenceComp <- rbind(ScenSumCost_NoFence, ScenSumCost_Fence, ScenSumCost_AES, ScenSumCost_NoFArable, ScenSumCost_FenceArable)
  
  
  
  ## Create Bar plot
  ggplot(FenceComp, aes(x= Strategy, y= Ave_PairCost*100000)) +
    geom_col(aes(fill = Fencing), width=0.55, position=position_dodge(0.6, preserve = "single")) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Created Habitat\n(orginal habitat)",   # Change legend title
                      values = c("#073b4c", "#ffd166", "#06d6a0",  "#ef476f", "#118ab2")) +
    ylab("Breeding Wader Pairs/ £100,000") +
    xlab("Lawton Principle") +
    labs(fill = "Targeting Strategy") +
    guides(fill=guide_legend(nrow=2, byrow=F)) +
    BarPlotTheme
  
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_FencingEffect.png"), units = "in", height = 9, width = 11)
  
  
  
  
  
  ##------------------------------------------------------##
  ##------------------------------------------------------##
  ##  Plot Scenarios for Land Purchase vs Income Foregone ##
  ##------------------------------------------------------##
  ##------------------------------------------------------##
  
  ## Calculate change in breeding pairs for reserve creation on grassland and pay income foregone
  ScenSumCost_FG <- AllScn2 |>
    filter(!OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Created = "Reserve- Grassland Income Foregone") |> 
    filter(NewCat == "Reserve")
  
  ## Calculate change in breeding pairs for reserve creation on grassland and purchase land
  ScenSumCost_PU <- AllScn2 |>
    filter(!OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCostPU)) |> 
    ungroup() |> 
    mutate(Created = "Reserve- Purchase Grassland") |> 
    filter(NewCat == "Reserve")
  
  
  ## Calculate change in breeding pairs for reserve creation on arable land and pay income foregone
  ScenSumCost_FGArable <- AllScn2 |>
    filter(OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Created = "Reserve- Arable Income Foregone") |> 
    filter(NewCat == "Reserve")
  
  ## Calculate change in breeding pairs for reserve creation on arable land and purchase land
  ScenSumCost_PUArable <- AllScn2 |>
    filter(OppCat == "Arable Opp") |> 
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCostPU)) |> 
    ungroup() |> 
    mutate(Created = "Reserve- Purchase Arable") |> 
    filter(NewCat == "Reserve")
  
  
  ## Calculate change in breeding pairs for AES creation
  ScenSumCost_AES <- AllScn2 |>
    mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
           PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
    summarise(Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(Created = "AES Only") |> 
    filter(NewCat == "AES Only")
  
  
  ## Combine all data sets for plotting
  CostComp <- rbind(ScenSumCost_FG, ScenSumCost_PU, ScenSumCost_AES, ScenSumCost_FGArable, ScenSumCost_PUArable)
  
  
  ## Create Bar plot
  ggplot(CostComp, aes(x= Strategy, y= Ave_PairCost*100000)) +
    geom_col(aes(fill = Created), width=0.55, position=position_dodge(0.6, preserve = "single")) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Created Habitat-\n(land costs)",   # Change legend title
                      values = c("#073b4c", "#ffd166", "#06d6a0",  "#ef476f", "#118ab2")) +
    ylab("Breeding Wader Pairs/ £100,000") +
    xlab("Lawton Principle") +
    labs(fill = "Targeting Strategy") +
    guides(fill=guide_legend(nrow=2, byrow=F)) +
    BarPlotTheme
  
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_ForegonePurchase.png"), units = "in", height = 9, width = 12)
  
  
  
  
  ##-------------------------------------##
  ##-------------------------------------##
  ##  Plot Scenario Average Cluster Size ##
  ##-------------------------------------##
  ##-------------------------------------##  
  
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSum <- AllScn2 |>
    mutate(PlotCatFull = paste0(Strategy, " & ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
    group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
    summarise(ClustActual= mean(ClustActual))
  
  ## Read in the 
  Sizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")

  ## Create Bar plot
  ggplot(ScenSum, aes(x= PlotCatFull, y= ClustActual)) +
    geom_col(aes(fill = ScenType), width=0.55, position=position_dodge(0.55)) +
    facet_wrap(~PlusFull) +
    scale_fill_manual(name = "Scenario Type",   # Change legend title
                        labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                        values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
    geom_hline(yintercept=Sizes$ReserveAve, linetype="dashed", color = "grey", linewidth = 1.2) +
    geom_hline(yintercept=Sizes$AESAve, linetype="dashed", color = "grey", linewidth = 1.2) +
    ylab("Average cluster size/ ha") +
    xlab("Scenario Category") +
    labs(fill = "Targeting Strategy") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_AveClusterSizeha_vs_Category.png"), units = "in", height = 9, width = 11)
  
  
  
  
  
  
  ##------------------------------------##
  ##------------------------------------##
  ##  Plot Scenarios Lapwing & Redshank ##
  ##------------------------------------##
  ##------------------------------------##
  
  if(species == "LapRed"){
    
    ## Read in all these files and bind them together
    ## Here just calculating the change in breeding pairs per unit cost for Lapwing/Redshank separately 
    AllScn2 <- files %>%
                map(fread) %>% # read in all files into a list
                bind_rows() %>% # bind all the rows
                mutate(ScenType = ifelse(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T), "clusterlarge", ScenType),
                       ScenType = ifelse(ScenType == "cluster" & ClustMean == min(ClustMean, na.rm = T), "clustersmall", ScenType),
                       ChangeSnipe = Snipe-BaseSnipe,
                       ChangeLapwing = Lapwing-BaseLapwing,
                       ChangeRedshank = Redshank-BaseRedshank,
                       ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                       ChangeWadersnoF = (Lapwing_unfenced+Redshank_unfenced) - (BaseLapwing+BaseRedshank),
                       ChangeCostsFG = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                            ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                       PairCostLap = ChangeLapwing/ChangeCostsFG,
                       PairCostRed = ChangeRedshank/ChangeCostsFG,
                       Plus2 = ifelse(Plus==T, "+", "-"),
                       Arable2 = ifelse(OppCat== "Arable Opp", "Arable", ""),
                       PlotCat = paste0(NewCat, Plus2, " ", Arable2)) |>
                       select(-c(Plus2, Arable2)) |>
                filter(!SegmentArea == 0) |> filter(!NewCat == "Reserve&AES Only")
    
    
    
    ## Calculate the change in Lapwing pairs per unit cost
    ScenSumLap <- AllScn2 |>
      mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
             PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
             PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
      group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
      summarise(Ave_PairCost = mean(PairCostLap)) |> 
      ungroup() |> 
      mutate(Species = "Lapwing")
    
    ## Calculate the change in Redshank pairs per unit cost
    ScenSumRed <- AllScn2 |>
      mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
             PlotCatFull = paste0(Strategy, " & ", NewCat, Arable2),
             PlusFull = ifelse(Plus==FALSE, "Normal Sampling", "Plus Sampling")) |>
      group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
      summarise(Ave_PairCost = mean(PairCostRed)) |> 
      ungroup() |> 
      mutate(Species = "Redshank")
      
    
    ## Combine the data sets for Lapwing and Redshank
    SpeciesComb <- rbind(ScenSumLap, ScenSumRed)
    
    
    ## Create Bar plot
    ggplot(SpeciesComb, aes(x= PlotCatFull, y= Ave_PairCost*100000)) +
      geom_col(aes(fill = Species), width=0.55, position=position_dodge(0.6)) +
      facet_wrap(~PlusFull) +
      scale_fill_manual(name = "Species",   # Change legend title
                        values = c("#6464f4", "#f46464")) +
      ylab("Breeding Wader Pairs/ £100,000") +
      xlab("Scenario Category") +
      BarPlotTheme + 
      theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
    
    
    ## save the plot
    ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_LapwingRedshankComparison.png"), units = "in", height = 9, width = 11)
  
  }

}













# ## CODE TO TEST SCENARIO FUNCTION FOR NORTH KENT
# 
# ##---------------------------##
# ##---------------------------##
# ## 2. *North Kent Scenarios* ##
# ##---------------------------##
# ##---------------------------##
# 
# ##------------------##
# ## 2.1 Read in data ##
# ##------------------##
# 
# ## Read in canvas polygons
# Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)
# table(duplicated(Canvshp$ParcRef))
# 
# ## Read in annotated canvas as csv
# Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")
# table(duplicated(Canv$ParcRef))
# 
# 
# 
# ## Join polygons to annotations
# Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
# rm(Canvshp)
# 
# ## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
# ## This will remove any fields that are not suitable for wader management
# Canv <- mutate(Canv,
#                Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
#                WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
#                Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
#                ParcArea = m2_to_ha(ParcArea), # match RF variable format
#                GroupArea = m2_to_ha(GroupArea), # convert to hectares
#                Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
#                # initiate columns for wader abundance
#                SnipeAbund = NA,
#                LapAbund = NA,
#                RedAbund = NA) |>
#         rename(FieldArea=ParcArea)
# 
# ## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#            theme_minimal()
# 
# 
# 
# ## Read in the UKCEH landcover data as a base map to rasterize my fields with
# LC <- rast("RawData/LandCover/gblcm25m2021.tif")
# LC <- LC[[1]]
# LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
# plot(LC_crop) # free up space
# rm(LC); gc()
# 
# ## Read in wader models
# LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")
# RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds")
# 
# 
# 
# ### NOTE: I could put the next two sections just inside the function itself,would streamline the whole thing a bit ###
# 
# ##-----------------------------##
# ## 2.2 Calc starting abundance ##
# ##-----------------------------##
# 
# ## Get the index of the rows that I want to predict wader abundance into
# PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
# 
# ## Predict abundance for Snipe
# set.seed(1212)
# Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod,
#                              newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()),
#                              jitt=FALSE))$predicted
# sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
# 
# 
# ## Predict abundance for Lapwing
# set.seed(1212)
# Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod,
#                           newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()),
#                           jitt=FALSE))$predicted
# sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape
# 
# 
# 
# ##-------------------------##
# ## 2.3 Calc starting costs ##
# ##-------------------------##
# 
# ## Calculate the costs of the AES agreements in the canvas
# 
# ## Read in the AES scheme
# ## Remove the small field supplement as not going to use this option for now
# AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
# SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
# AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")
# 
# 
# ## Initiate a column for creating the costings
# Canv$AESCostGDP <- 0
# Canv$ResMaintGDP <- 0
# Canv$ResCreatGDP <- 0
# 
# ## Work out which rows in the data set are AES fields
# WhichAES <- which(Canv$Category == "AES Only")
# 
# ## Now using the list column of the different AES schemes this function looks through all the AES codes
# ## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
# for(j in 1:nrow(AES_Costs)){
# 
#   ## create columns used to track AES payment
#   Canv$Add <- NA
#   Canv$AddAmount <- 0
# 
#   ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
#   Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})
# 
#   ## Calculate the cost for that field of the current AES scheme
#   ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
#   ## Also if the field is less then 1ha add on a small field supplement
#   Canv[WhichAES,] <- Canv[WhichAES,] %>%
#                      mutate(AddAmount = ifelse(Add==T & AES_Costs$Unit[j]=="ha", FieldArea*AES_Costs$Cost[j],
#                                               ifelse(Add==T & AES_Costs$Unit[j]=="m", Perim*AES_Costs$Cost[j], 0)),
#                             AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
#                             AESCostGDP = AESCostGDP + AddAmount)
# 
#   if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
# }
# 
# ## Calculate the costs in £1000's of pounds
# sum(Canv$AESCostGDP[WhichAES], na.rm=T)/1000
# 
# 
# 
# ## Calculate the costs of the reserve management in the canvas
# 
# ## Work out which rows are reserves
# WhichReserve <- which(Canv$Category == "Reserve")
# 
# ## Work out the costs of ongoing management of the field per hectare
# Canv$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canv$GroupArea[WhichReserve],
#                                                ParcelArea = Canv$FieldArea[WhichReserve],
#                                                InflationAdjust=T)
# 
# ## Calculate the costs of reserve management in £1000's of pounds
# sum(Canv$ResMaintGDP[WhichReserve], na.rm=T)/1000
# 
# 
# 
# 
# ##-------------------------##
# ## 2.4 Batch Run Scenarios ##
# ##-------------------------##
# 
# 
# ## Read in the AES scheme
# ## Remove the small field supplement as not going to use this option for now
# # CostSheet <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv") |>
# #              filter(!CSS_Code == "SP1")
# #
# # AES_Costs
# 
# 
# ## Read in spreadsheet with setting for all scenarios I want to run
# Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")
# 
# ## Read in spreadsheet with sizes for clusters
# clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
# Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == min(Scen_settings$ClustMean, na.rm=T), clustsizes$AESAve, Scen_settings$ClustMean)
# Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == max(Scen_settings$ClustMean, na.rm=T), clustsizes$ReserveAve, Scen_settings$ClustMean)
# 
# ## For Broads I want to remove large cluster from Better scenarios
# ## For now also do not try and run any stakeholder scenarios as gradings not yet calculated
# Scen_settingsNK<- Scen_settings #|> filter(!(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T) & Strategy == "Better"))
# Scen_settingsNK$OppCat4 <- NA
# 
# ## Run the batch function for additive
# # BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
# #           Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/", SaveCanv = FALSE,
# #           N_sets = 5, runsetting = Scen_settingsNK, zstart = 1, Additive = TRUE)
# 
# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/"
# SaveCanv = FALSE
# N_sets = 5
# runsetting = Scen_settingsNK
# zstart = 1
# Additive = TRUE
# 
# 
# 
# ##-----------------------------##
# ## F 2.1 Calc Opp Area Targets ##
# ##-----------------------------##
# 
# ## run `AreaForScenarios` function to calculate the size of the opportunity areas for each scenario
# ## only need to do this for random as it is the same for clustered appraoches
# AreasList <- AreaForScenarios(Canvas=Canvas, ScenList=(runsetting |> filter(ScenType == "random")))
# 
# ## Calculate the minimum opportunity area for each strategy/lawton principle
# AreasSum <- AreasList |>
#              group_by(Strategy) |>
#              summarise(MinOppArea = min(OppArea, na.rm = T))
# 
# ## Now join the MinOppArea onto the full list of scenarios
# runsetting <- left_join(runsetting, AreasSum, by = "Strategy")
# 
# ## When additive = TRUE I need a vector of areas to test that are the same length as N_sets
# ## I could just repeat the same area but I though it would be good to introduce a bit of variance
# ## For now I will use the cluster sizes from the large clusters
# if(Additive==TRUE){
#     clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
#     threshold_areaList <- rtruncnorm(n=N_sets, a=0, mean = clustsizes$ReserveAve, sd = max(runsetting$ClustSD, na.rm=T))
#   }
# 
# 
# 
# ##---------------------------##
# ## F 2.2 Loop over Scenarios ##
# ##---------------------------##
# 
# z=1
# 
# ##Join together opportunity and replacement categories for CreateScenario() function
# OppCatChoice = c(paste0(if(is.na(runsetting$OppCat1[z])==F){paste0(runsetting$OppCat1[z])}),
#                  paste0(if(is.na(runsetting$OppCat2[z])==F){paste0(runsetting$OppCat2[z])}),
#                  paste0(if(is.na(runsetting$OppCat3[z])==F){paste0(runsetting$OppCat3[z])}),
#                  paste0(if(is.na(runsetting$OppCat4[z])==F){paste0(runsetting$OppCat4[z])}))
# 
# NewCatChoice = c(paste0(if(is.na(runsetting$NewCat1[z])==F){paste0(runsetting$NewCat1[z])}),
#                  paste0(if(is.na(runsetting$NewCat2[z])==F){paste0(runsetting$NewCat2[z])}))
# 
# ## create the file path for output
# ## Create a column earlier on to differentiate small vs large cluster
# Typ <- runsetting$ScenType[z]
# FinalOutPath = paste0(Outpath, if(Typ == "random"){paste0("rand", "/")},
#                                if(Typ == "cluster" & runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustlarge/")},
#                                if(Typ == "cluster" & !runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustsmall/")},
#                                if(Typ == "stakeholder2" | Typ == "stakeholder1"){paste0("stake", runsetting$StakeGroup[z], "/")})
# 
# ## Define what will be the OppAreaUsed, this is either a single number of the total area across all N_sets
# ## Or a vecotr the same length of N_sets that dictate how much management is deployed on each sampling loop
# if(Additive==TRUE){OppAreaChoice <- threshold_areaList}
# if(Additive==FALSE){OppAreaChoice <- runsetting$MinOppArea[z]}
# 
# # Create directory if needs be
# dir.create(FinalOutPath, showWarnings = F)
# 
# 
# # ## Run all of the scenarios
# # CreateScenario(Canvas=Canvas, SniModel=SniModel, LapModel=LapModel, RedModel=RedModel, LandCov=LandCov, N_sets=N_sets, Additive = Additive,
# #                ScenType = runsetting$ScenType[z], Strategy = runsetting$Strategy[z], StakeGroup = runsetting$StakeGroup[z],
# #                OppCat= OppCatChoice, NewCat= NewCatChoice, Plus = runsetting$Plus[z],  OppAreaUsed = OppAreaChoice,
# #                ClustBuf= runsetting$ClustBuf[z], ClustMean= runsetting$ClustMean[z], ClustSD= runsetting$ClustSD[z],
# #                Outpath = FinalOutPath, SaveCanv = SaveCanv)
# Canvas=Canvas
# SniModel=SniModel
# LapModel=LapModel
# RedModel=RedModel
# LandCov=LandCov
# N_sets=N_sets
# Additive = Additive
# ScenType = runsetting$ScenType[z]
# Strategy = runsetting$Strategy[z]
# StakeGroup = runsetting$StakeGroup[z]
# OppCat= OppCatChoice
# NewCat= NewCatChoice
# Plus = runsetting$Plus[z]
# OppAreaUsed = OppAreaChoice
# ClustBuf= runsetting$ClustBuf[z]
# ClustMean= runsetting$ClustMean[z]
# ClustSD= runsetting$ClustSD[z]
# Outpath = FinalOutPath
# SaveCanv = SaveCanv
# Unfenced=FALSE