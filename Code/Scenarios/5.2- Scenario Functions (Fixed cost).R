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
   cluster_buffer <- st_buffer(st_union(cluster), dist = clustbuff) 
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

Create_Costs <- function(ParcelArea, StartCategory="grass", StartHab = "grass",
                         CreateCost = CrCost$Cost[CrCost$Element =="Creation costs"],
                         ReversionCost = CrCost$Cost[CrCost$Element =="Arable reversion"],
                         Correction = CrCost$Cost[CrCost$Element =="Reserve correction"]){
  
  ## ha to m2 conversion factor
  ha_to_m = 10000
  
  ## calculate the log of the total cost for the area provided
  TotCost = ParcelArea*CreateCost
  
  ## Add the values in for this at a later date
  TotCost = ifelse(StartHab == "arable",  TotCost+(ParcelArea*ReversionCost), TotCost)
  
  ## If the field was already a reserve then change the creation costs
  if(StartCategory == "Reserve"){
    TotCost <- TotCost*(Correction/100)
  }
  
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
                          OppCat1 = NA,
                          OppCat2 = NA,
                          NewCat=  NA, 
                          OppArea = NA)
  
  ## run through each row in the data frame of scenarios
  for(j in 1:nrow(ScenList)){
    
    ##Join together opportunity and replacement categories for CreateScenario() function
    OppCatChoice = c(paste0(if(is.na(ScenList$OppCat1[j])==F){paste0(ScenList$OppCat1[j])}),
                     paste0(if(is.na(ScenList$OppCat2[j])==F){paste0(ScenList$OppCat2[j])}))
    NewCatChoice = paste0(ScenList$NewCat1[j])
  
    Strategy = ScenList$Strategy[j]
    
    ## Get the index of rows representing the opportunity area for this scenario
    ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
    if(Strategy %in% c("Big", "More") & NewCatChoice == "AES Only"){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice) | 
                                                                             (is.na(Canvas$ClustGroup)==T & Canvas$Category %in% NewCatChoice & Canvas$BrAES == F)) }
    if(Strategy %in% c("Better") & NewCatChoice == "AES Only"){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice) |
                                                                          (is.na(Canvas$ClustGroup)==F & Canvas$Category %in% NewCatChoice & Canvas$BrAES == F)) }
    
    
    ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders 
    ## (i.e. AES fields with waders will not be targeted and those without will be)
    if(Strategy %in% c("Big", "More") & NewCatChoice == "Reserve"){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice) | 
                                                                              (is.na(Canvas$ClustGroup)==T & Canvas$Category %in% NewCatChoice & Canvas$RSPB == "N")) }
    if(Strategy %in% c("Better") & NewCatChoice == "Reserve"){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice) | 
                                                                         (is.na(Canvas$ClustGroup)==F & Canvas$Category %in% NewCatChoice & Canvas$RSPB == "N")) }
    
    ## Only sample the appropriate Arable fields if doing arable reversion
    if(Strategy %in% c("Big", "More") & NewCatChoice == "Reserve" & 
       any(OppCatChoice %in% "Arable Opp")){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice)) }
    if(Strategy %in% c("Better") & NewCatChoice == "Reserve" & 
       any(OppCatChoice %in% "Arable Opp")){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice)) }
    
    
    ## Calculate the total field Area for my opportunity Area
    ScenAreas$OppArea[j] <- sum(Canvas[IndOpp,]$FieldArea)
    
    ScenAreas$Strategy[j] <- Strategy
    ScenAreas$OppCat1[j] <- ScenList$OppCat1[j]
    ScenAreas$OppCat2[j] <- ScenList$OppCat2[j]
    ScenAreas$NewCat[j] <- paste0(NewCatChoice)
    
  } ## End of j loop
  
  return(ScenAreas)
  
} 





##-----------------------------------------------##

# Canvas=Canvas
# SniModel=SniModel
# LapModel=LapModel
# RedModel=RedModel
# LandCov=LandCov
# N_sets=N_sets
# Budget = runsetting$Budget[z]
# ClustBuf= ClustBuf
# OppCat= OppCatChoice
# NewCat= NewCatChoice
# Strategy = runsetting$Strategy[z]
# CostSpread = CostSpread
# SaveCanv = F
# Outpath = Outpath



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
                           ScenType, # currently defunct 
                           SniModel, 
                           LapModel,
                           RedModel, 
                           LandCov,
                           N_sets, 
                           Budget,
                           ClustBuf= 200, 
                           OppCat, 
                           NewCat, 
                           Strategy, 
                           CostSpread = 20,
                           MinClust,
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
  message("Going to run ", N_sets, " repeats")
  
  
  
  ##-- CREATE TRACKER DOCUMENT --##
  
  ## Create data frame to record wader abundance at the end of each segment 
  ## Also record the starting abundance in a separate columns
  ScTracker <- data.frame(Landscape = Canvas$Landscape[1],
                          #ScenType =  ScenType,
                          Strategy = Strategy,
                          OppCat= paste(OppCat,collapse="&"),
                          NewCat=  paste(NewCat,collapse="&"), 
                          Budget = Budget,
                          BudgetReached = NA,
                          CostSpread = paste0(CostSpread, "yrs"),

                          AveClustArea = NA,
                          AveClustWider = NA, # for reserve creation this is the adjusted cluster area, when adding area of newly created reserve and connected existing reserve
                          ClustNum = NA,
                          Segment = 1:(N_sets+1), 
                          SegmentArea = NA,
                          
                          BaseAESCost = sum(Canvas$AESCostGDP[Canvas$Category=="AES Only"], na.rm=T),
                          NewAESCost = NA,
                          BaseResMaintCost = sum(Canvas$ResMaintGDP[Canvas$Category=="Reserve"], na.rm=T),
                          NewResOverallCost = NA,
                          NewResMaintCost = NA,
                          NewResCreatCost = NA,
                          NewResFencingCost = NA,
                          NewResPurchCost = NA,
                          
                          Snipe = NA,
                          BaseSnipe = ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T)),
                          Lapwing = NA, 
                          BaseLapwing = ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T)),
                          Redshank = NA,
                          BaseRedshank = ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T)),
                          SegTime = NA)
  
  ## Set the first row as zero area and baseline abundance
  ScTracker$Snipe[1] <- ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T))
  ScTracker$Lapwing[1] <- ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T))
  ScTracker$Redshank[1] <- ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T))
  ScTracker$SegmentArea[1] <- 0
  
  
  
  ##-- DEFINE OPPORTUNITY AREA --##
  
  ## Get the index of rows representing the opportunity area for scenarios involving AES creation
  ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
  ## Grassland opportunity fields are selected and AES fields that do not have a breeding wader option
  if(Strategy %in% c("Big", "More") & NewCat == "AES Only"){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat) | 
                                                                             (is.na(Canvas$ClustGroup)==T & Canvas$Category %in% NewCat & Canvas$BrAES == F)) }
  if(Strategy %in% c("Better") & NewCat == "AES Only"){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat) |
                                                                        (is.na(Canvas$ClustGroup)==F & Canvas$Category %in% NewCat & Canvas$BrAES == F)) }
  
  
  ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders 
  ## (i.e. AES fields with waders will not be targeted and those without will be)
  if(Strategy %in% c("Big", "More") & NewCat == "Reserve"){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat) | 
                                                                            (is.na(Canvas$ClustGroup)==T & Canvas$Category %in% NewCat & Canvas$ResQual == "LowQ")) }
  if(Strategy %in% c("Better") & NewCat == "Reserve"){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat) | 
                                                                       (is.na(Canvas$ClustGroup)==F & Canvas$Category %in% NewCat & Canvas$ResQual == "LowQ ")) }
  
  ## Only sample the appropriate Arable fields if doing arable reversion
  if(Strategy %in% c("Big", "More") & NewCat == "Reserve" & 
     any(OppCat %in% "Arable Opp")){ IndOpp <- which((is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat)) }
  if(Strategy %in% c("Better") & NewCat == "Reserve" & 
     any(OppCat %in% "Arable Opp")){ IndOpp <- which((is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat)) }
  

  ##-- DEFINE FIELDS TO USE AS SAMPLES OF NEW HABITAT --##
  
  ## These are the fields from which I will sample costs and habtiat values
  ## Pre-compute the AES rate per ha and trim down to two columns to speed up the loop
  AESBr_fields <- filter(Canvas, Category == "AES Only" & BrAES == T) |> 
                  mutate(AES_Rate = AESCostGDP/FieldArea) |> 
                  select(ParcRef, AES_Rate)
  
  QReserve_fields <- filter(Canvas, Category == "Reserve" & ResQual == "HighQ") |> 
                  mutate(AES_Rate = AESCostGDP/FieldArea) |> 
                  select(ParcRef, AES_Rate, Fence_Coverage)
  
  

  
  ##---------------------------------------##  
  ##---------------------------------------##
  ## Label fields in each scenario segment ##
  ##---------------------------------------##
  ##---------------------------------------##

  ## For each Opportunity field randomly choose a suitable breeding wader AES field to pair it with
  ## This provides the costs for the field but also the habitat values that the field will receive later on 
  

  ##-------------------------------------------##  
  #### F 1.2 Label fields for AES deployment ####
  ##-------------------------------------------##  
  
  
  ##-- FIELD LABELLING --##
  
  if(NewCat == "AES Only"){
    
    ## Message to update console
    message("Building AES creation scenarios...")
    
    
    
    ##-- PREPARE SUB-CANVAS  --##
    
    ## filter out the fields that are within the opportunity area for this scenario
    MoneyCols <- c("AESCostGDP") # cols with costing
    OppFields <- Canvas[IndOpp,] |> select(all_of(c("FieldArea", "ParcRef", "ClustDist", MoneyCols))) # could filter the column of OppFields here as later on just need ParcRef & RandSamp
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     ClustID=NA,
                                     CentDist=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)),
                                     CurrentAESRate = AESCostGDP/FieldArea)


    
    ##-- CREAST LIST FOR LOOP RECORDING --##
    
    ## Create an empty list where I can put the parcel references that belong to each round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    listClustID <- vector(mode = "list", length = N_sets)
    listPairParcref <- vector(mode = "list", length = N_sets)
    listReserveGroup <- vector(mode = "list", length = N_sets)
    listPosit <- vector(mode = "list", length = N_sets)
    
    
    
    ##-- START LOOP --##
    
    ## Label opportunity fields for each round of sampling
    for(k in 1:N_sets){ 
      
      ## Message for progress of cluster sampling
      message("Creating scenario ", k , " out of ", N_sets)
      

      ##-- SELECT SUB-SET OF COLUMNS --##
      
      ## Filter out the fields that are still in the opportunity area
      ## For each iteration of the loop this will gradually decrease the size of the data set I am working with as more fields become part of a cluster
      OppFieldsLoop <- OppFields |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", "CurrentAESRate", MoneyCols))) |> 
                                                       mutate(ClustID=NA, CentDist=NA, RandSamp=NA, PairRefs = NA)
      
      
      ##-- ASSIGN PAIR FIELD AND CALCULATE COSTS --##
      
      ## Sample fields to the same length as the opportunity fields
      PairRefs = sample_n(AESBr_fields, nrow(OppFieldsLoop), replace = TRUE)
      
      ## Now calculate the change in costs for a parcel as it may already have a wintering wader AES scheme
      OppFieldsLoop$PairRefs <- PairRefs$ParcRef
      OppFieldsLoop$NewAESRate <- PairRefs$AES_Rate
      OppFieldsLoop <- OppFieldsLoop |> 
                       mutate(ChangeRate = NewAESRate-CurrentAESRate,
                              ChangeRate = ifelse(ChangeRate < 129, 129, ChangeRate),
                              AESCostGDP = ChangeRate*FieldArea) |> 
                       select(-c(CurrentAESRate, ChangeRate, NewAESRate))

      
      
      
      ##-- SET LOOP CONTROLS AND OUTPUTS --##
      
      ## set the total area cluster for this round of the loop to zero
      ## This only resets outside of all the while loops
      TotalCost <- 0
      
      ## On first loop set the shortfall to 0 or always return to zero if Additive==T
      Spent <- 0
      
      ## Counter for while loop
      ## And a list to store the cluster in if the loop has to run multiple times
      WhileCount = 1
      SmallClustCount = 1
      ClustStore <- vector(mode = "list")
      
    
      
      ##-- 1st WHILE LOOP START --##
      
      ## While total area of all clusters created in this iteration of the while function is < AreaSeg, keep on creating more clusters
      ## Need to add on NextGapExtra to make sure sampling isn't too much or too little (depending on if previous loops over/undersampled)
      while(TotalCost < Budget) {
      
      
      ##-- SAMPLE CLUSTER NUCLEUS AND INTIATE CLUSTER --##
      
      # Select an initial polygon as the initiation point of the cluster ( this will depend upon the strategy chosen)
      # If ScenType =="cluster" weight the random draw by the distance to cluster (more) or inverse distance to cluster (Bigger)
      if(WhileCount==1){OppFieldsLoop2 <- OppFieldsLoop}
      if(WhileCount>1){OppFieldsLoop2 <- OppFieldsLoop |> filter(!ParcRef %in% All_Clusters$ParcRef)} # remove fields already in cluster from earlier iterations of this while loop
      
      # break out of while loop if there is no more rows left and record this in tracker document
      if(nrow(OppFieldsLoop2) ==0){
        message("All opportunity fields used")
        ScTracker$BudgetReached[k+1] <- "Y"
        break}  
      
      #message("Sampling starting point")
      if(Strategy == "Big"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$InvClustDist) 
      }
      if(Strategy == "More"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$PClustDist) 
      }
      if(Strategy == "Better"){ initial_polygon <- sample_n(OppFieldsLoop2, 1) 
      }
        
        
        
      ##-- START SECOND WHILE LOOP --##
      ## Keep on expanding the cluster until I run out of neighbours or the budget is exceeded
      ## If budget is exceeded then break the loop from inside
      neighborsL <- 1
      cluster <- initial_polygon
      Fullcluster_cost <- sum(cluster$AESCostGDP)
      
      while(neighborsL > 0){

      ##-- FIND NEIGHBOURS --##
      #message("finding neighbours")
      # Create data set of potential fields for this cluster, should speed up sampling
      PotentialRefs <- (OppFieldsLoop2 |> st_crop(st_buffer(cluster, dist = 10000)))$ParcRef 
      PotentialRefs <- PotentialRefs[PotentialRefs %in% cluster$ParcRef == FALSE] # remove fields already in the expanding cluster
      ClustOppFields <- filter(OppFieldsLoop2, ParcRef %in% c(PotentialRefs))
      if(nrow(ClustOppFields)==0){
        message("Can't find any neighbours")
        if(nrow(cluster)==1){ CheckCost <- Fullcluster_cost + Spent }
        break
      }
      # plot(ClustOppFields$geometry)
      
      
      ## Retrieve potential neighbors
      neighbors <- get_neighbors(cluster, ClustOppFields, clustbuff = ClustBuf)
      neighborsL <- nrow(neighbors)
      # plot(neighbors["FieldArea"]) # check
      
      ## Calculate the area of the all the neighbors and initial field
      CluserSize <- sum(neighbors$FieldArea) + sum(cluster$FieldArea)
      
      
      
      ##-- BUDGET CHECK --##
      
      Fullcluster <- rbind(neighbors, cluster)
      
      ## This is the average annual cost of the AES
      Fullcluster_cost <- sum(Fullcluster$AESCostGDP)
      
      ## Calculate the total expenditure, including the cost of previous reserve creation
      CheckCost <- Fullcluster_cost + Spent
      
      
      
      ##-- IF ENTIRE NEW CLUSTER IS BELOW BUDGET --##
      
      if(CheckCost < Budget){
        
        ## if under budget then just 
        cluster <- Fullcluster
        
        ## Assign the final cluster cost as well
        Finalcluster_cost <- Fullcluster_cost
        
      }
        
        
      
      ##-- IF ENTIRE NEW CLUSTER IS ABOVE BUDGET --##
      
      ## If I have gone over budget then just add fields sequentially until
      if(CheckCost >= Budget){
        
      message("Budget exceeded, reducing size of cluster") 
        
      ## Calculate the distance between potential neighbors and cluster center, then order with the closest fields first
      Fullcluster <- Fullcluster |> mutate(CentDist = as.numeric(st_distance(initial_polygon, Fullcluster))) |> 
                                    arrange(CentDist)
        
      ## calculate the budget remaining
      RemainingCost <- Budget - Spent
        
      ## Select the closest neighbors sequentially until the total area exceeds the threshold_area for the entire cluster
      indT <- which(cumsum((Fullcluster$AESCostGDP)) < (RemainingCost)) # rows of fields that go just up to threshold but not over threshold
      if(length(indT)>0){cluster <- Fullcluster[indT,]} # select fields chosen above and add one extra field to just take over threshold
      if(length(indT) == 0){} # just use first fields if none selected to take it just over threshold
        
        }
        

      ## If the budget has been reached then break out of the while loop
      ## The break above only breaks out of the j for loop
      if(CheckCost>=Budget){break}
        
      } # END OF 2nd WHILE LOOP
      
      
      
      ##-- CHECK CLUSTER IS LARGE ENOUGH --##
      
      ## If the cluster is below a threshold then do not use it
      ## Only do this if the budget is not reached as to reach the budget might require a small cluster
      ## calculate the size of the cluster
      ClustArea <- sum(cluster$FieldArea)
      
      ## if the cluster does not exceed our threshold and not tried 5 times
      if(ClustArea < MinClust & CheckCost < Budget & SmallClustCount < 5){
        
        ## Store cluster if it is too small, add 1 to loop counter and go on to next iteration of while loop
        ClustStore[[SmallClustCount]] <- cluster 
        SmallClustCount = SmallClustCount+1 
        next } 
      ## if the cluster is not large enough but I have tried 5 times
      if(ClustArea < MinClust & CheckCost < Budget & SmallClustCount >=5){
        ## If loop run lots of times then calculate the cluster with the largest area and choose that as the cluster to go forward with

        ClustStore[[SmallClustCount]] <- cluster 
        MAXClust <- which.max(do.call(rbind, lapply(ClustStore, function(i)colSums((i |> st_drop_geometry())['FieldArea']))))
        cluster <-  ClustStore[[MAXClust]]
        ClustStore <- vector(mode = "list")
        SmallClustCount = 1
      }
      

      ##-- SORT OUT CHOSEN FINAL CLUSTER --##
      
      message("Storing a cluster")
      
      ## Now I have my final cluster, either the full cluster or the one trimmed down to size
      ## First give it a cluster ID, this is just a loop counter for the while loop
      cluster$ClustID <- WhileCount
      
      ## Now assign which fields were chosen and but in the OppFieldsLoop data set
      cluster <- cluster |>  arrange(ParcRef) # arrange by ParcRef as this is how main Canvas is ordered
      clusterInd <- which(OppFieldsLoop$ParcRef %in% cluster$ParcRef)
      OppFieldsLoop$RandSamp[clusterInd] <- k
        
      
      ## Assign the final cluster to a larger data set of all clusters
      if(WhileCount==1){All_Clusters <- cluster}
      if(WhileCount>1){All_Clusters <- rbind(All_Clusters, cluster)}
      
      ## Add one on to the loop counter as this records how many cluster I have made
      WhileCount <- WhileCount+1
      
      ## Reset cluster store and loop counter for small cluster
      ClustStore <- vector(mode = "list")
      SmallClustCount = 1
        
      ## update how much has been spent
      Spent <- Finalcluster_cost + Spent
      
      
      ## assign the check cost to the total cost
      ## This will break loop if I went over the limit
      ## If i did not go over the limit then it will go around to the next iteration and assign the total spent to TotalCost
      message("Check cost total: £", round(CheckCost/1000), "k")
      TotalCost <- CheckCost


    } ## 1st WHILE LOOP END -- stop creating clusters
      
      
      
    ##-- RECORD CLUSTER PROPOERTIES IN LISTS --##  

  
    ## Assign the parcel refs chosen in this round of sampling to the list
    ## Also assign the ClustID to a list so each ParcRef can be matched to a ClustID number
    All_Clusters <- All_Clusters |>  arrange(ParcRef)
    listSamples[[k]] <- All_Clusters$ParcRef
    listClustID[[k]] <- All_Clusters$ClustID
    listPairParcref[[k]] <- All_Clusters$PairRefs
    
    ## Assign stats about the random sample I have drawn to the document
    ScTracker$NewAESCost[k+1] <- sum(All_Clusters$AESCostGDP)
    ScTracker$SegmentArea[k+1] <- sum(All_Clusters$FieldArea)
    ScTracker$AveClustArea[k+1] <- mean(aggregate(FieldArea ~ ClustID, All_Clusters, sum)$FieldArea) ## This might need improving as multiple cluster might actually be the same
    ScTracker$ClustNum[k+1] <- WhileCount-1
    
      
  }# end of loop k 
    
        
  } # end of AES scenario sampling

  
  
  
  
  ##-----------------------------------------------##  
  #### F 1.3 Label fields for Reserve deployment ####
  ##-----------------------------------------------##   
  
  
  ##-- FIELD LABELLING --##
  
  if(NewCat == "Reserve"){
    
    ## Message to update console
    message("Building Reserve creation scenarios...")

    
    ##-- PREPARE SUB-CANVAS  --##

    ## filter out the fields that are within the opportunity area for this scenario
    MoneyCols <- c("ForegoneGDP", "PurchaseGDP", "ResMaintGDP", "ResCreatGDP")
    OppFields <- Canvas[IndOpp,] |> select(all_of(c("FieldArea", "ParcRef", "ClustDist", "ReserveGroup", "GroupArea", "ArableOpp", MoneyCols))) # could filter the column of OppFields here as later on just need ParcRef & RandSamp
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     ClustID=NA,
                                     CentDist=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    # plot(OppFields["ReserveGroup"])

    
    
    ##-- CREAST LIST FOR LOOP RECORDING --##
    
    ## Create an empty list where I can put the parcel references that belong to each round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    listClustID <- vector(mode = "list", length = N_sets)
    listPairParcref <- vector(mode = "list", length = N_sets)
    listReserveGroup <- vector(mode = "list", length = N_sets)
    listPosit <- vector(mode = "list", length = N_sets)
    
    
    
    ##-- START LOOP --##
    
    ## Label opportunity fields for each round of sampling
    for(k in 1:N_sets){ 
      
      ## Message for progress of cluster sampling
      message("Creating scenario ", k , " out of ", N_sets)
      

      
      ##-- SELECT SUB-SET OF COLUMNS --##
      
      ## Filter out the fields that are still in the opportunity area
      ## For each iteration of the loop this will gradually decrease the size of the data set I am working with as more fields become part of a cluster
      # OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols)))
      OppFieldsLoop <- OppFields |> select(all_of(c("FieldArea", "ParcRef", "ReserveGroup", "GroupArea", "RandSamp", "InvClustDist", 
                                                    "PClustDist", "ClustID", "CentDist", "ArableOpp", MoneyCols))) |> 
                                    mutate(ClustID=NA, CentDist=NA, RandSamp=NA, PairRefs = NA)
      
            
            
      ##-- ASSIGN PAIR FIELD FOR EACH OPPORTUNITY FIELD --##
      
      ## sample parcel ID for pair fields
      PairRefs = sample_n(QReserve_fields, nrow(OppFieldsLoop), replace = TRUE)

      ## Assign the parcel ID and the fence indicator variable to my main data set
      OppFieldsLoop$FenceCov <- PairRefs$Fence_Coverage
      OppFieldsLoop$PairRefs <- PairRefs$ParcRef

      
      
  
      ##-- SET WHILE LOOP CONTROLS AND OUTPUTS --##
      
      ## set the total area cluster for this round of the loop to zero
      ## This only resets outside of all the while loops
      TotalCost = 0
      
      ## On first loop set the shortfall to 0 or always return to zero if Additve==T
      Spent = 0
      FenceSpent = 0
      
      ## Counter for while loop
      ## And a list to store the cluster in if the loop has to run multiple times
      WhileCount = 1
      SmallClustCount = 1
      ClusterCount = 1
      ClustStore = vector(mode = "list")
      CostStore = vector(mode = "list")
      FenceStore = vector(mode = "list")

      
    
      
      ##-- 1st WHILE LOOP START --##
      
      ## UPDATE
      while(TotalCost < Budget) {
        
        
        
      ##-- CREATE POLYGONS OF ALL RESERVES --##
    
      ## Create a combined set of outlines for all RSPB reserves, need this if creating reserve
      ReserveFields <- Canvas |> filter(Reserve=="Y") |> select(ParcRef, FieldArea) 
      ReserveOutlines <- Canvas |> filter(Reserve=="Y") |> group_by(ReserveGroup) |>  
                                   summarize(GroupArea= max(GroupArea),
                                             geometry = st_union(geometry)) 
      # plot(ReserveOutlines)
      
      
      
      ##-- SAMPLE CLUSTER NUCLEUS AND INTITIATE CLUSTER --##
      
      # Select an initial polygon as the initiation point of the cluster ( this will depend upon the strategy chosen)
      # If ScenType =="cluster" weight the random draw by the distance to cluster (more) or inverse distance to cluster (Bigger)
      if(WhileCount==1){OppFieldsLoop2 <- OppFieldsLoop}
      if(WhileCount>1){OppFieldsLoop2 <- OppFieldsLoop |> filter(!ParcRef %in% All_Clusters$ParcRef)} # remove fields already in cluster from earlier iterations of this while loop
      
      # break out of while loop if there is no more rows left and record this in tracker document
      if(nrow(OppFieldsLoop2) ==0){
        message("All opportunity fields used")
        ScTracker$BudgetReached[k+1] <- "Y"
        break}  
      
      if(Strategy == "Big"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$InvClustDist) 
      }
      if(Strategy == "More"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$PClustDist) 
      }
      if(Strategy == "Better"){ initial_polygon <- sample_n(OppFieldsLoop2, 1) 
      }
  
      
      ##-- START SECOND WHILE LOOP --##
      ## Keep on expanding the cluster until I run out of neighbours or the bidget is exceeded
      ## If budget is exceeded then break the loop from inside
      neighborsL <- 1
      cluster <- initial_polygon
      
      while(neighborsL > 0){
        
      ##-- FIND NEIGHBOURS --##
      
      # Create data set of potential fields for this cluster, should speed up sampling
      # OppFieldsLoop2 has already had field selected from previously created cluster removed
      PotentialRefs <- (OppFieldsLoop2 |> st_crop(st_buffer(cluster, dist = 10000)))$ParcRef 
      PotentialRefs <- PotentialRefs[PotentialRefs %in% cluster$ParcRef == FALSE] # remove fields already in the expanding cluster
      ClustOppFields <- filter(OppFieldsLoop2, ParcRef %in% c(PotentialRefs))
      # plot(ClustOppFields$geometry)
      
      ## Retrieve potential neighbors
      neighbors <- get_neighbors(cluster, ClustOppFields, clustbuff = ClustBuf)
      neighborsL <- nrow(neighbors)
      # plot(neighbors["FieldArea"]) # check
      
      ## Calculate the area of the all the neighbors and initial field
      CluserSize <- sum(neighbors$FieldArea) + sum(cluster$FieldArea)
      # cluster <- rbind(initial_polygon, neighbors)   
      
      
      ##-- BUDGET CHECK --##
      
      Fullcluster <- rbind(neighbors, cluster)
      
      ## Work out if this cluster will expand an existing reserve or not
      ## If it does then give it the same Reserve Group ID as the reserve it is near
      ## If not then just give it a unique code
        
      ## Calculate if any existing reserves are near to my newly created reserve
      BuffClust <- Fullcluster |> st_buffer(dist=ClustBuf) |> summarize(geometry = st_union(geometry))
      # plot(BuffClust$geometry)
      
      ## Test for overlap new reserve and existing reserves
      Overlap <- st_overlaps(BuffClust, ReserveOutlines)
      
      
      ## if no overlap then just give the cluster a unique grouping code that can not be confused with others
      if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) == 0L){
        Fullcluster$ReserveGroup <- paste("Group", WhileCount, "_", k)
        Fullcluster$GroupArea <- sum(Fullcluster$FieldArea)
        
        }
        
      ## If there is overlap then give the new reserve cluster the same ReserveGroup name as the reserve that it overlaps with
      if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) > 0L){
        
        if(length(Overlap[[1]])==1){ResNumber <- Overlap[[1]]}
        if(length(Overlap[[1]])>1){ResNumber <- Overlap[[1]][1]}
        
        ## Do they have any fields in common, if they do then need to subtract this from total cluster area
        OverlapFields <- filter(ReserveFields, ParcRef %in% Fullcluster$ParcRef)

        ## Assign a reserve group code and calcualte the total area of the cluster
        Fullcluster$ReserveGroup <- paste0(ReserveOutlines$ReserveGroup[ResNumber])
        Fullcluster$GroupArea <- ReserveOutlines$GroupArea[ResNumber] + sum(Fullcluster$FieldArea) - sum(OverlapFields$FieldArea)

      }

      
      ## Work out the costs of ongoing management of the field per hectare
      Fullcluster$ResMaintGDP <- Manage_Costs(TotalArea = Fullcluster$GroupArea, 
                                                  ParcelArea = Fullcluster$FieldArea,
                                                  InflationAdjust=T)
     
      ## Work out the costs of site creation for each parcel that has been created
      ## if the parcels changed to reserve were originally arable then add on an extra creation cost
      ## StartCategory argument not used here
      Fullcluster$ResCreatGDP <- (Create_Costs(ParcelArea= Fullcluster$FieldArea, 
                                               StartHab = ifelse(Fullcluster$ArableOpp == 1, "arable","grass")))/CostSpread
      
      ## First just select the parcel that have been altered in this round of the i loop
      ## This reduces the size of the data set to speed up processing a bit as well
      FenceIndex <- which(Fullcluster$FenceCov == "Y")
      
      ## Finally calculate the cost to fence this area if it were a rectangle with length = 1.5*width
      ## Divide by 15 years here as we are estimating that the full fence will need replacing every 15 years
      FenceCost <- SimpFenceCost(ParcelArea= sum(Fullcluster$FieldArea[FenceIndex]))/15
      
      ## This is the average annual cost over the next 20 years
      Fullcluster_cost <- sum(Fullcluster$ResMaintGDP) + sum(Fullcluster$ResCreatGDP) + FenceCost + sum(Fullcluster$PurchaseGDP)/CostSpread
      
      
      ## Calculate the total expenditure, including the cost of previous reserve creation
      CheckCost <- Fullcluster_cost + Spent
      
      
      ## I have created and costed a cluster
      ## If it reaches the budget then need to trim
      ## If it does not reach budget then go to next round of while loop to add to cluster
      
      
      ##-- IF ENTIRE NEW CLUSTER IS BELOW BUDGET --##

      if(CheckCost < Budget){

        ## if under budget then just 
        cluster <- Fullcluster
        
        ## Assign the final cluster cost as well
        Finalcluster_cost <- Fullcluster_cost
        Finalfencecost <- FenceCost
        
      }
        
        
      
      ##-- IF ENTIRE NEW CLUSTER IS ABOVE BUDGET --##
      
      ## If I have gone over budget then just add fields sequentially until
      if(CheckCost >= Budget){
        
      print("Budget exceeded, reducing size of cluster") 
        
      ## Calculate the distance between potential neighbors and cluster center, then order with the closest fields first
      neighbors <- neighbors |> mutate(CentDist = as.numeric(st_distance(initial_polygon, neighbors))) |> 
                                arrange(CentDist)
      
      ## Now loop through all the neighbors and sequentially add them too to the cluster
      ## Create a data set to record all the costs
      Notebook <- data.frame(Run = 0:nrow(neighbors), Cost = NA, FenceCost = NA)
      TrimStore <- vector(mode = "list")

      
      for(j in 0:nrow(neighbors)){
        
        ## create the cluster for this round of the loop
        if(j==0){Testcluster <- cluster}
        if(j>0){Testcluster <- rbind(cluster, neighbors[1:j,])}
        
        
        ## Now calculate the cost of the cluster
        ## Work out if this cluster will expand an existing reserve or not
        ## If it does then give it the same Reserve Group ID as the reserve it is near
        ## If not then just give it a unique code
          
        ## Calculate if any existing reserves are near to my newly created reserve
        TestBuffClust <- Testcluster |> st_buffer(dist=ClustBuf) |> summarize(geometry = st_union(geometry))
        # plot(BuffClust$geometry)
        
        ## Test for overlap new reserve and existing reserves
        Overlap <- st_overlaps(TestBuffClust, ReserveOutlines)
        
        
        ## if no overlap then just give the cluster a unique grouping code that can not be confused with others
        if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) == 0L){
          Testcluster$ReserveGroup <- paste("Group", WhileCount, "_", k)
          Testcluster$GroupArea <- sum(Testcluster$FieldArea)
          
          }
        
        ## If there is overlap then give the new reserve cluster the same ReserveGroup name as the reserve that it overlaps with
        if(is.integer(Overlap[[1]]) && length(Overlap[[1]]) > 0L){
          
          if(length(Overlap[[1]])==1){ResNumber <- Overlap[[1]]}
          if(length(Overlap[[1]])>1){ResNumber <- Overlap[[1]][1]}
          
          ## Do they have any fields in common, if they do then need to subtract this from total cluster area
          OverlapFields <- filter(ReserveFields, ParcRef %in% Testcluster$ParcRef)

          ## Assign a reserve group code and calcualte the total area of the cluster
          Testcluster$ReserveGroup <- paste0(ReserveOutlines$ReserveGroup[ResNumber])
          Testcluster$GroupArea <- ReserveOutlines$GroupArea[ResNumber] + sum(Testcluster$FieldArea) - sum(OverlapFields$FieldArea)

        }
  
        
        ## Work out the costs of ongoing management of the field per hectare
        Testcluster$ResMaintGDP <- Manage_Costs(TotalArea = Testcluster$GroupArea, 
                                                    ParcelArea = Testcluster$FieldArea,
                                                    InflationAdjust=T)
       
        ## Work out the costs of site creation for each parcel that has been created
        ## if the parcels changed to reserve were originally arable then add on an extra creation cost
        ## StartCategory argument not used here
        Testcluster$ResCreatGDP <- (Create_Costs(ParcelArea= Testcluster$FieldArea, 
                                                     StartHab = ifelse(Testcluster$ArableOpp == 1, "arable","grass")))/CostSpread
        
        ## First just select the parcel that have been altered in this round of the i loop
        ## This reduces the size of the data set to speed up processing a bit as well
        FenceIndex <- which(Testcluster$FenceCov == "Y")
        
        ## Finally calculate the cost to fence this area if it were a rectangle with length = 1.5*width
        ## Divide by 15 years here as we are estimating that the full fence will need replacing every 15 years
        FenceCost <- SimpFenceCost(ParcelArea= sum(Testcluster$FieldArea[FenceIndex]))/15
        
        ## This is the average annual cost over the next 20 years
        Testcluster_cost <- sum(Testcluster$ResMaintGDP) + sum(Testcluster$ResCreatGDP) + FenceCost + sum(Testcluster$PurchaseGDP)/CostSpread
        Notebook$Cost[j+1] <- Testcluster_cost
        Notebook$FenceCost[j+1] <- FenceCost
        TrimStore[[j+1]] <- Testcluster # store in j+1 as can not store when j=0
        
        ## Calculate the total expenditure, including the cost of previous reserve creation
        TestCost <- Testcluster_cost + Spent
        
        
        ## If the budget is not exceeded then try adding on another field, at some point this loop will break the budget
        if(TestCost<Budget){next}
        
        ## the budget is exceeded then I not the tipping point
        ## Can now pick the appropriate sized reserve and assign it to Finalcluster
        if(TestCost>=Budget){
          
          ## Return the last cluster trialed, NOTE cluster are stored in slot j+1
          if(j==0){cluster <- cluster}
          if(j>0){cluster <- TrimStore[[j]]}
          
          ## Assign the final cluster cost as well
          ## use j-1 here as this will make sure I am just under the budget, the AES strategies always come in just under
          if(j==0){
           Finalcluster_cost <- Notebook$Cost[1]
           Finalfencecost <- Notebook$FenceCost[1]
          }
          if(j>0){
           Finalcluster_cost <- Notebook$Cost[j]
           Finalfencecost <- Notebook$FenceCost[j]}

          
          ## Now I have reached the budget break out of the loop
          break
          
          }
        
        } # END OF j LOOP: trimming clusters over budget
        
        }
      
      ## If the budget has been reached then break out of the while loop
      ## The break above only breaks out of the j for loop
      if(CheckCost>=Budget){break}
        
      } # END OF 2nd WHILE LOOP
      
      
      ##-- CHECK CLUSTER IS LARGE ENOUGH --##
      
      ## If the cluster is below a threshold then do not use it
      ## Only do this if the budget is not reached as to reach the budget might require a small cluster
      ## calculate the size of the cluster
      ClustArea <- sum(cluster$FieldArea)
      
      ## if the cluster does not exceed our threshold and not tried 5 times
      if(ClustArea < MinClust & CheckCost < Budget & SmallClustCount < 5){
        message("Cluster too small storing it")
        ## Store cluster if it is too small, add 1 to loop counter and go on to next iteration of while loop
        ClustStore[[SmallClustCount]] <- cluster 
        CostStore[[SmallClustCount]] <- Finalcluster_cost
        FenceStore[[SmallClustCount]] <- Finalfencecost
        SmallClustCount = SmallClustCount+1 
        next } 
      
      ## if the cluster is not large enough but I have tried 5 times
      if(ClustArea < MinClust & CheckCost < Budget & SmallClustCount >=5){
        ## If loop run lots of times then calculate the cluster with the largest area and choose that as the cluster to go forward with
        ## Store the values for the current iteraiton
        ClustStore[[SmallClustCount]] <- cluster 
        CostStore[[SmallClustCount]] <- Finalcluster_cost
        FenceStore[[SmallClustCount]] <- Finalfencecost
        
        ## select the largest cluster in terms of area
        message("Maxed out small cluster, picking biggest")
        MAXClust <- which.max(do.call(rbind, lapply(ClustStore, function(i)colSums((i |> st_drop_geometry())['FieldArea']))))
        
        ## Retrieve the largest cluster and its total cost and fence cost
        cluster <-  ClustStore[[MAXClust]]
        Finalcluster_cost <- CostStore[[MAXClust]]
        Finalfencecost <- FenceStore[[MAXClust]]
        CostStore = vector(mode = "list")
        FenceStore = vector(mode = "list")
        ClustStore = vector(mode = "list")
        SmallClustCount = 1
      }
        
  
      ##-- SORT OUT CHOSEN FINAL CLUSTER --##
      
      message("Storing a cluster")
      
      ## Now I have my final cluster, either the full cluster or the one trimmed down to size
      ## First give it a cluster ID, this is just a loop counter for the while loop
      cluster$ClustID <- WhileCount
      
      ## Now assign which fields were chosen and but in the OppFieldsLoop data set
      cluster <- cluster |>  arrange(ParcRef) # arrange by ParcRef as this is how main Canvas is ordered
      clusterInd <- which(OppFieldsLoop$ParcRef %in% cluster$ParcRef)
      OppFieldsLoop$RandSamp[clusterInd] <- k
        
      
      ## Assign the final cluster to a larger data set of all clusters
      if(WhileCount==1){All_Clusters <- cluster}
      if(WhileCount>1){All_Clusters <- rbind(All_Clusters, cluster)}
      
      ## Add one on to the loop counter as this records how many cluster I have made
      WhileCount <- WhileCount+1
        
      ## update how much has been spent
      stopifnot(is.na(Finalcluster_cost)==F)
      message("Updating spent... ", Finalcluster_cost, " + ", Spent)
      Spent <- Finalcluster_cost + Spent
      FenceSpent <- Finalfencecost + FenceSpent
      
      ## Reset cluster store and loop counter for small cluster
      ClustStore <- vector(mode = "list")
      CostStore = vector(mode = "list")
      FenceStore = vector(mode = "list")
      SmallClustCount = 1

      
      ## assign the check cost to the total cost
      ## This will break loop if I went over the limit
      ## If i did not go over the limit then it will go around to the next iteration and assign the total spent to TotalCost
      message("Check cost total: £", round(CheckCost/1000), "k")
      TotalCost <- CheckCost
      

    } ## 1st WHILE LOOP END -- stop creating clusters
      
      
      
      
    ##-- RECORD CLUSTER PROPERTIES IN LISTS --##  

    ## Assign the parcel refs chosen in this round of sampling to the list
    ## Also assign the ClustID to a list so each ParcRef can be matched to a ClustID number

    # plot(All_Clusters["ClustID"])
    All_Clusters <- All_Clusters |>  arrange(ParcRef)
    listSamples[[k]] <- All_Clusters$ParcRef
    listClustID[[k]] <- All_Clusters$ClustID
    listPairParcref[[k]] <- All_Clusters$PairRefs
    
    ## Assign stats about the random sample I have drawn to the document
    ScTracker$NewResOverallCost[k+1] <- Spent
    ScTracker$NewResFencingCost[k+1] <- FenceSpent
    ScTracker$NewResMaintCost[k+1] <- sum(All_Clusters$ResMaintGDP)
    ScTracker$NewResCreatCost[k+1] <- sum(All_Clusters$ResCreatGDP)
    ScTracker$NewResPurchCost[k+1] <- sum(All_Clusters$PurchaseGDP)/CostSpread
    ScTracker$SegmentArea[k+1] <- sum(All_Clusters$FieldArea)
    ScTracker$AveClustArea[k+1] <- mean(aggregate(FieldArea ~ ClustID, All_Clusters, sum)$FieldArea) ## This might need improving as multiple cluster might actually be the same
    ScTracker$AveClustWider[k+1] <- mean(aggregate(GroupArea ~ ClustID, All_Clusters, max)$GroupArea)
    ScTracker$ClustNum[k+1] <- WhileCount-1
    
      
  }# end of loop k 
    
        
  } # end of Reserve scenario sampling
  
  
  
  

  
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
    #### F 1.4 Update within field habitats ####
    ##----------------------------------------##
    
    ## send message to console
    message("Updating habitat variables")
    
    
    ##-- IDENTIFY FIELDS IN THE CURRENT SET --##
            
    ## Get the row numbers of the opportunity fields in the current random sampling group
    ## Update here so that sampled fields are drawn from a list of parcel refs and not from a column
    SetInd <- which(Canvas$ParcRef %in% listSamples[[i]])
    
    ## Now label the fields in this loop number with the loop count to show that they have been updated from their original state
    Canvas$Upgraded[SetInd] <- i
    
    
    
    ##-- CHOOSE HABITAT COLUMNS TO UPDATE --##
    
    ## define which within habitat categories I am going to update, also update the category for that row
    ## For Arable also update the TALL_BOUNDARY_PERCENT, normally for grassland fields I leave this the same but in arable it needs to be updates
    HabCats <- c("Category", "GrassOpp", "GRASSLAND_TYPE", "WaterCoverage", "STOCK", "VEG_STRUCTURE", "RUSH_PERCENT", "Fence_Coverage", "AES_options")
    if(any(OppCat %in% "Arable Opp")){HabCats <- c(HabCats, "TALL_BOUNDARY_PERCENT")}


    
    ##-- UPDATE WITHIN FIELD HABITAT --##
      
    ## Get the row numbers of the fields from which I am going to sample habitat from
    ## These were created earlier as each field that is being converted needed to have a cost
    NewHabs <- match(listPairParcref[[i]], Canvas$ParcRef)
   
    ## Now update the within field habitat columns
    Canvas[SetInd, colnames(Canvas) %in% HabCats] <- Canvas[NewHabs, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
    print(table(Canvas$Category)) 
    
    ## If creating reserve also alter the reserve column
    if(any(NewCat %in% "Reserve")){ Canvas$Reserve[SetInd] <- "Y" }

    
    
    
    ##------------------------------------##    
    #### F 1.5 Update landscape habitat ####
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
    #### F 1.6 Predict Wader abundance ####
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
    
    
    
    ##------------------------------------------##
    #### F 1.7 Save current iteration outputs ####
    ##------------------------------------------##
    
    ##-- SAVE TRACKER DOC AFTER EVERY LOOP --##
    
    ## First calculate the average cluster size for this round of sampling
    # Create a data frame to combine cluster ID and field parcel areas
    DaTa <- data.frame(Area = Canvas$FieldArea[SetInd], Groups = listClustID[[i]])
    
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
                              Strategy, "_",
                              paste0(round(Budget/1000),"k"), 
                              paste(NewCat,collapse="&"), 
                              if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
    
    ## If doing the additive testing, i.e. only can look at additive effects
    ## Then reset the canvas at this point, this will mean that the in each round of the i loop we go back to the starting canvas
    Canvas <- MASTERCanv

    
    
  } # End of segment loop (i)
  
  
  
  ##-- FINAL SAVE TRACKER DOC AFTER ALL LOOPS --##
  
  ## return the tracker data set
  print(ScTracker)
    
  ## finally write out the file once again
  write_csv(ScTracker, paste0(Outpath, 
                              "Track_", 
                              Strategy, "_", 
                              paste0(round(Budget/1000),"k"), 
                              paste(NewCat,collapse="&"), 
                              if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
  
  ##-- SAVE LISTS OF SAMPLES --##
  
  write_rds(listSamples, paste0(Outpath, "Samples/",
                                "ParcRefs_", 
                                Strategy, "_", 
                                paste0(round(Budget/1000),"k"), 
                                paste(NewCat,collapse="&"), 
                                if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".rds"))
  
  write_rds(listClustID, paste0(Outpath, "Samples/",
                                "ClustID_", 
                                Strategy, "_", 
                                paste0(round(Budget/1000),"k"), 
                                paste(NewCat,collapse="&"), 
                                if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".rds"))
    
  
  
} ## FUNCTION `CreateScenario` END




##-----------------------------------------------##




# inpath="CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/"
# outpath="CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/Plots/"
# species="LapRed"

## This function takes the scenario modelling output files from the CreateScenario function and 
## plots the outputs to plot a whole series of plots for a write up
PlotScenario <- function(inpath,
                         outpath,
                         species){
  

  ##-----------------##
  ##-----------------##
  ##  Plot Scenarios ##
  ##-----------------##
  ##-----------------##
  
  dir.create(outpath, showWarnings = F)
  
  # Identify all scenario tracker files that were saved
  files <- dir(paste0(inpath), pattern  = "Track_", full.names  = T)

  ## Read in all these files and bind them together
  AllScn <- files %>%
              map(fread) %>% # read in all files into a list
              bind_rows() %>% # bind all the rows
              mutate(ChangeSnipe = Snipe-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeCosts = ifelse(NewCat=="AES Only", NewAESCost,
                                          ifelse(NewCat=="Reserve", NewResOverallCost, 0)),
                     PlotCat = paste0(NewCat, " ", ifelse(OppCat== "Arable Opp", "Arable", ""))) |>
              filter(!SegmentArea == 0) 
  
  ## Calculate the cost per pairs depending on whether it is Snipe or Lapwing/Redshank
  if(species == "LapRed"){ AllScn <- AllScn |> mutate(PairCost = ChangeWaders/(ChangeCosts/Budget),
                                                      PairCost100k = ChangeWaders/(ChangeCosts/100000),
                                                      Waders_100Ha = (ChangeWaders/SegmentArea)*100) }
  
  if(species == "Snipe"){ AllScn <- AllScn |> mutate(ChangeWaders = ChangeSnipe,
                                                     PairCost = ChangeSnipe/(ChangeCosts/Budget),
                                                     PairCost100k = ChangeWaders/(ChangeCosts/100000),
                                                     Waders_100Ha = (ChangeWaders/SegmentArea)*100) }
  
  
  ## Set a general theme for the bar plots
  BarPlotTheme <- theme_light() +
                theme(panel.grid.minor = element_blank(),
                      strip.text = element_text(size = 17, face = "bold"),
                      axis.title = element_text(size = 18),
                      axis.text.y = element_text(size = 16),
                      axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5),
                      legend.title = element_text(size = 17, face = "bold"),
                      legend.text = element_text(size = 16),
                      legend.position = "top")
  
  
  
  ##-------------------------------##
  ##  Pairs per 100 ha vs Strategy ##
  ##-------------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSum <- AllScn |>
    mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""))) |>
    group_by(PlotCatFull, NewCat, Strategy) |>
    summarise(StEr = sd(Waders_100Ha)/sqrt(n()),
              t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
              MeanWaders_100Ha = mean(Waders_100Ha)) |> 
    ungroup() |> 
    mutate(LowerCI = MeanWaders_100Ha - (t_score * StEr),
           UpperCI = MeanWaders_100Ha + (t_score * StEr))
 

  ## Create Bar plot
  ggplot(ScenSum, aes(x= PlotCatFull, y= MeanWaders_100Ha)) +
    geom_col(aes(fill = Strategy), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(.9), colour = "#424949") +
    scale_fill_manual(name = "Strategy",   # Change legend title
                        labels = c("Better"= "Better", "More"= "More", "Big" ="Bigger"),  # Change legend labels
                        values = c("Better"="#A5F076", "More"="#76A5F0", "Big"="#F076A5")) +
    ylab("Change in Breeding Wader Pairs/ 100 ha") +
    xlab("Scenario Category") +
    labs(fill = "Targeting Strategy") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Pairha_vs_AmalgBudget.png"), units = "in", height = 9, width = 11)
  
  
  
  
  ##------------------------------------##
  ##  Pairs per budget cost vs Strategy ##
  ##------------------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSumCost <- AllScn |>
    mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""))) |>
    group_by(PlotCatFull, NewCat, Strategy) |>
    summarise(StEr = sd(PairCost)/sqrt(n()),
              t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
              Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(LowerCI = Ave_PairCost - (t_score * StEr),
           UpperCI = Ave_PairCost + (t_score * StEr))
  
  
  ## Create Bar plot
  ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost)) +
    geom_col(aes(fill = Strategy), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(.9), colour = "#424949") +
    scale_fill_manual(name = "Strategy",   # Change legend title
                        labels = c("Better"= "Better", "More"= "More", "Big" ="Bigger"),  # Change legend labels
                        values = c("Better"="#A5F076", "More"="#76A5F0", "Big"="#F076A5")) +
    ylab("Breeding Wader Pairs/Budget Cost") +
    xlab("Scenario Category") +
    labs(fill = "Targeting Strategy") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "PairBudget_vs_AmalgBudget.png"), units = "in", height = 9, width = 12)
  
  
  
   
  ##-------------------------------------##
  ## Pairs per 100 ha vs Strategy/Budget ##
  ##-------------------------------------##
  
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSum <- AllScn |>
    mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                           ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                  paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
    group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
    summarise(StEr = sd(Waders_100Ha)/sqrt(n()),
              t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
              MeanWaders_100Ha = mean(Waders_100Ha)) |> 
    ungroup() |> 
    mutate(LowerCI = MeanWaders_100Ha - (t_score * StEr),
           UpperCI = MeanWaders_100Ha + (t_score * StEr),
           BudgetW = fct_reorder(BudgetW, Budget))
 
  ## Create Bar plot
  ggplot(ScenSum, aes(x= PlotCatFull, y= MeanWaders_100Ha, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Change in Breeding Wader Pairs/ 100 ha") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Pairha_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  
  
  
  
  ##-------------------------------------------##
  ##  Pairs per budget cost vs Strategy/Budget ##
  ##-------------------------------------------##
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSumCost <- AllScn |>
    mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                           ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                  paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
    group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
    summarise(StEr = sd(PairCost)/sqrt(n()),
              t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
              Ave_PairCost = mean(PairCost)) |> 
    ungroup() |> 
    mutate(LowerCI = Ave_PairCost - (t_score * StEr),
           UpperCI = Ave_PairCost + (t_score * StEr),
           BudgetW = fct_reorder(BudgetW, Budget))
 
  ## Create Bar plot
  ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Change in Breeding Wader Pairs/Budget cost") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "PairBudget_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  

  
  ##-------------------------------------##
  ##  Pairs per £100k vs Strategy/Budget ##
  ##-------------------------------------##

  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSumCost <- AllScn |>
    mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
           BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                           ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                  paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
    group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
    summarise(StEr = sd(PairCost100k)/sqrt(n()),
              t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
              Ave_PairCostK = mean(PairCost100k)) |> 
    ungroup() |> 
    mutate(LowerCI = Ave_PairCostK - (t_score * StEr),
           UpperCI = Ave_PairCostK + (t_score * StEr),
           BudgetW = fct_reorder(BudgetW, Budget))
 
  ## Create Bar plot
  ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCostK, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Change in Breeding Wader Pairs/£100,000") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "Pair£100k_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  
  
  

  ##----------------------------------------##
  ## Total Scenario Area vs Strategy/Budget ##
  ##----------------------------------------##

  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ScenSize <- AllScn |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(SegmentArea)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            AveSegmentArea = mean(SegmentArea)) |> 
  ungroup() |> 
  mutate(LowerCI = AveSegmentArea - (t_score * StEr),
         UpperCI = AveSegmentArea + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

    
  ## Create Bar plot
  ggplot(ScenSize, aes(x= PlotCatFull, y= AveSegmentArea, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Area altered in scenario/ha") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "TotalArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  
  
  
  ##----------------------------------------##
  ## Total Cluster Area vs Strategy/Budget ##
  ##----------------------------------------##

  ## First create a data set that can be used to create a bar plot that covers all scenarios
  ClustSize <- AllScn |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(AveClustArea)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            AveAveClustArea = mean(AveClustArea)) |> 
  ungroup() |> 
  mutate(LowerCI = AveAveClustArea - (t_score * StEr),
         UpperCI = AveAveClustArea + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

    
  ## Create Bar plot
  ggplot(ClustSize, aes(x= PlotCatFull, y= AveAveClustArea, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Cluster area/ha") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "ClusterArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  
  
  
  ##----------------------------------------##
  ## Wider Cluster Area vs Strategy/Budget ##
  ##----------------------------------------##

  ## First create a data set that can be used to create a bar plot that covers all scenarios
  WClustSize <- AllScn |>
  filter(NewCat == "Reserve") |> 
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(AveClustWider)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            AveAveClustWider = mean(AveClustWider)) |> 
  ungroup() |> 
  mutate(LowerCI = AveAveClustWider - (t_score * StEr),
         UpperCI = AveAveClustWider + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

    
  ## Create Bar plot
  ggplot(WClustSize, aes(x= PlotCatFull, y= AveAveClustWider, group = BudgetW)) +
    geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
    geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
    scale_fill_manual(name = "Budget",
                        values = c("#A5F076", "#76A5F0", "#F076A5")) +
    ylab("Wider Cluster area/ha") +
    xlab("Scenario Category") +
    labs(fill = "Budget") +
    BarPlotTheme + 
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))
  
  ## save the plot
  ggsave(plot=last_plot(), filename= paste0(outpath, "ClusterArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
  

}













# ## CODE TO TEST SCENARIO FUNCTION FOR NORTH KENT ####
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
#   rename(FieldArea=ParcArea)
# 
# ## Just do this here so that i do not upgrade Elmley
# Canv$RSPB <- ifelse(Canv$ReserveGroup == "Elmley", "Y", Canv$RSPB)
# 
# 
# ## Add a new column that separates out breeding wader vs wintering wader payment
# Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})
# 
# 
# ## Need to add column hat differentiates high quality vs lower quality reserves
# 
# 
# ## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#   theme_minimal()
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
#                                          newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()),
#                                          jitt=FALSE))$predicted
# sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
# 
# 
# ## Predict abundance for Lapwing
# set.seed(1212)
# Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod,
#                                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()),
#                                          jitt=FALSE))$predicted
# sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape
# 
# 
# 
# ##-----------------------------##
# ## 2.3 Calc starting costs ##
# ##-----------------------------##
# 
# ##-- Calculate the costs of the AES agreements in the canvas --##
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
# 
# ## Work out which rows in the data set are AES fields
# WhichAES <- which(Canv$Category == "AES Only")
# 
# ## Now using the list column of the different AES schemes this function looks through all the AES codes
# ## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
# for(j in 1:nrow(AES_Costs)){
#   
#   ## How many years to spread one off costs over
#   Yrs <- 15
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
#                      mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
#                                                   Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
#                                                   Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
#                                                   Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
#                                                   .default = 0),
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
# ##-- Calculate the costs of the reserve management in the canvas --##
# 
# ## Initiate a column for creating the costings
# Canv$ResMaintGDP <- 0
# Canv$ResCreatGDP <- 0
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
#   group_by(Strategy) |>
#   summarise(MinOppArea = min(OppArea, na.rm = T))
# 
# ## Now join the MinOppArea onto the full list of scenarios
# runsetting <- left_join(runsetting, AreasSum, by = "Strategy")
# 
# ## When additive = TRUE I need a vector of areas to test that are the same length as N_sets
# ## I could just repeat the same area but I though it would be good to introduce a bit of variance
# ## For now I will use the cluster sizes from the large clusters
# if(Additive==TRUE){
#   clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
#   threshold_areaList <- rtruncnorm(n=N_sets, a=0, mean = clustsizes$ReserveAve, sd = max(runsetting$ClustSD, na.rm=T))
# }
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
#                       if(Typ == "cluster" & runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustlarge/")},
#                       if(Typ == "cluster" & !runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustsmall/")},
#                       if(Typ == "stakeholder2" | Typ == "stakeholder1"){paste0("stake", runsetting$StakeGroup[z], "/")})
# 
# ## Define what will be the OppAreaUsed, this is either a single number of the total area across all N_sets
# ## Or a vecotr the same length of N_sets that dictate how much management is deployed on each sampling loop
# if(Additive==TRUE){OppAreaChoice <- threshold_areaList}
# if(Additive==FALSE){OppAreaChoice <- runsetting$MinOppArea[z]}
# 
# # Create directory if needs be
# dir.create(FinalOutPath, showWarnings = F)


# CreateScenario(Canvas=Canvas, SniModel=SniModel, LapModel=LapModel, RedModel=RedModel, LandCov=LandCov,
#                N_sets=N_sets, Budget = 300000, ClustBuf= 200, OppCat= c("Arable Opp"),
#                NewCat= "Reserve", Strategy = runsetting$Strategy[z],
#                CostSpread = 20, SaveCanv = F, Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/")

# Canvas=Canvas
# SniModel=SniModel
# LapModel=LapModel
# RedModel=RedModel
# LandCov=LandCov
# N_sets=N_sets
# Budget = 300000
# ClustBuf= 200
# OppCat= c("AES Only", "Grass Opp")
# NewCat= "Reserve"
# Strategy = runsetting$Strategy[z]
# CostSpread = 20
# SaveCanv = F
# Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/"
# 
# 
# "Arable Opp"

