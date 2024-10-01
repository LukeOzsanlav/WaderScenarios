
##-----------------------------------##
## F 1 START CreateScenario Function ##
##-----------------------------------##

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

## Defining function arguments
# Canvas= landscape canvas as a shape file variables need to match those in the RF model
# ScenType = What type of scenario should by run, can only pick one at a time. options are "random", "cluster", "stakeholder", "stakeholder"
# SniModel= RF model for predicting Snipe abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LapModel= RF model for predicting Lapwing abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# RedModel= RF model for predicting Redshank abundance (if don't predict want to predict this species use =NULL), make sure canvas has same columns as xvars in model
# LandCov= Landcpover data raster covering same areas as Canvas. This is used to create as raster out of Canvas to calculate landscape variables
# N_sets= number of sets to split up the total opportunity area into, each set is run and then wader abundance calculated
# Additive = should additive sampling be carried out (TRUE/FALSE). Additive sampling essentially resets between rounds of N_sets and if Additive=FALSE managemnt is sequentially added to the canvas
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
                           Additive=FALSE,
                           OppCat, 
                           NewCat, 
                           Plus, 
                           Strategy, 
                           PropOppArea=0.5, # currently defunct
                           OppAreaUsed, # new addition
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
  
  ## Remove all the fields that are not needed during the scenario modelling
  ## This is to try and increase the speed
  ## Any fields with no opportunity or that were masked are removed unless they are wider wet grassland adjacent to the priority landscape
  Canvas <- filter(Canvas, Category %in% c("AES Only", "Arable Opp",  "Grass Opp", "Reserve") | WiderWetGrass == 1)
  plot(Canvas$geometry)
  
  ## The number segments to split scenario modeling up into
  ## For each segment habitat variables are updated and wader abundance is calculated
  N_sets <- N_sets
  
  ## Create data frame to record wader abundance at the end of each segment 
  ## Also record the starting abundance in a seperate columns
  ScTracker <- data.frame(Landscape = Canvas$Landscape[1],
                          ScenType =  ScenType,
                          Strategy = Strategy,
                          OppCat= paste(OppCat,collapse="&"),
                          NewCat=  paste(NewCat,collapse="&"), 
                          Plus = Plus, 
                          # PropOppArea = PropOppArea, 
                          ClustMean = ClustMean,
                          Segment = 1:(N_sets+1), 
                          SegmentArea = NA,
                          SegTime = NA,
                          Snipe = NA, 
                          BaseSnipe = ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T)),
                          Lapwing = NA, 
                          BaseLapwing = ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T)),
                          Redshank = NA,
                          BaseRedshank = ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T)))
  
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
  
  if(ScenType =="cluster" | ScenType == "stakeholder2"){
    
    ## Message to update console
    message("Labelling segments with ", ScenType, "...")
    
    ## If using ScenType "stakeholder2" then will need to retain some extra columns along the way, if using cluster just set this to Null
    if(ScenType == "stakeholder2"){StakeCols <- c(paste0("Bigger_", StakeGroup), paste0("Better_", StakeGroup), 
                                                  paste0("More_", StakeGroup), paste0("Arable_", StakeGroup))}else{StakeCols <- NULL}
    
    
    ## filter out the fields that are within the opportunity area for this scenario
    OppFields <- Canvas[IndOpp,] |> select(all_of(c("FieldArea", "ParcRef", "ClustDist", StakeCols))) # could filter the column of OppFields here as later on just need ParcRef & RandSamp
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    
    ## Start a counter that will label each cluster sequentially
    ## And create columns to put counter in and a column that gives the distance to the center of the cluster
    ClusterCount <- 1
    OppFields <- OppFields |> mutate(ClustID=NA,
                                     CentDist=NA)

    # Define function used for finding neighbours around existing cluster
    # Function to get neighbors within 100m buffer of current cluster polygons, it excludes those already in the cluster
    get_neighbors <- function(cluster, all_polygons, clustbuff) {
  
       ## 50m around buffer cluster
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
    
    ## Create an empty list where I can put the parcel references that belong to each round of sampling
    listSamples <- vector(mode = "list", length = N_sets)
    listClustID <- vector(mode = "list", length = N_sets)
    listPosit <- vector(mode = "list", length = N_sets)
    
    ## Label opportunity fields for each round of sampling
    for(k in 1:N_sets){ 
      
      ## Message for progress of cluster sampling
      message("Creating set of clusters ", k , " out of ", N_sets)
      
      ## set the total area cluster for this round of the loop to zero
      ## This only resets outside of all the while loops
      TotClustArea <- 0
      
      ## On first loop set the shortfall to 0
      if(k==1 | Additive==TRUE){NextGapExtra <- 0}
      ## On the final loop set the shortfall as a large number so all fields are sampled
      ## REMOVED as currently no scenarios sample all the opportunity area
      # if(k==N_sets & PropOppArea == 1){NextGapExtra <- sum(OppFields$FieldArea)+1}
      
      ## Filter out the fields that are still in the opportunity area
      ## For each iteration of the loop this will gradually decrease the size of the data set I am working with as more fields become part of a cluster
      # OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols)))
      if(Additive==FALSE){OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols)))}
      if(Additive==TRUE){OppFieldsLoop <- OppFields |> select(all_of(c("FieldArea", "ParcRef", "RandSamp", "InvClustDist", "PClustDist", "ClustID", "CentDist", StakeCols))) |> 
                                                       mutate(ClustID=NA, CentDist=NA, RandSamp=NA)}
      
      ## Counter for while loop
      ## And a list to store the cluster in if the loop has to run multiple times
      WhileCount = 1
      ClustStore <- vector(mode = "list")
      
    
      ##-- 1st WHILE LOOP START --##
      ## While total area of all clusters created in this iteration of the while function is < AreaSeg, keep on creating more clusters
      ## Need to add on NextGapExtra to make sure sampling isn't too much or too little (depending on if previous loops over/undersampled)
      while(TotClustArea < (AreaSeg[k]+NextGapExtra)) {
        
    
      # Generate a random max threshold in hectares for current cluster
      threshold_area <- rtruncnorm(n=1, a=0, mean = ClustMean, sd = ClustSD)  # Example: between 10 and 50 hectares
      
      ## **INSERT COMMENTS....
      if(Additive==TRUE & ClustMean > 200){threshold_area <- AreaSeg[k]}

      
      ## **CONTROLS TO MODERATE THRESOLD AREA
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
        initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = (OppFieldsLoop2[,ColIn] |> st_drop_geometry() |> unlist())+0.000000001) # add small amount incase there are some 0 gradings
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
        
        
        ## Retrieve potential neighbors
        neighbors <- get_neighbors(cluster, ClustOppFields, clustbuff = ClustBuf)
        # plot(neighbors["FieldArea"]) # check
    
        # Break if no more neighbors are found (Could add option to try a larger buffer area, would need to add argument to get_neighbors())
        if(nrow(neighbors) == 0) break
        
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
        # plot(selected_neighbor["FieldArea"]) # check
        
        # Add the selected neighbour(s) to the cluster
        cluster <- rbind(cluster, selected_neighbor)
        # plot(cluster["FieldArea"]) # check
        
        # Update the cluster area
        cluster_area <- sum(cluster$FieldArea)  # Area in hectares
    
        
      } ## 2nd WHILE LOOP END -- stop expanding single cluster

      # Plot the final cluster
      # plot(st_geometry(cluster), border = 'blue', col = alpha('blue', 0.5))

      
      ## **CONTROLS TO ENSURE ONLY COMPLETED CLUSTERS USED
      ## if the full cluster was not made then at first do not use it
      ## Store cluster that are too small in a list, hoping that a cluster of the right size can be made
      ## If after 10 number of runs there have been no cluster large enough created then just use the largest one created so far
      if(cluster_area < threshold_area & WhileCount <10){
        ## Sotre cluster if it is too small, add 1 to loop counter and go on to next iteraiton of while loop
        ClustStore[[WhileCount]] <- cluster |> st_drop_geometry()
        WhileCount = WhileCount+1 
        next 
        } 
      if(cluster_area < threshold_area & WhileCount >=10){
        ## If loop run lots of times then calculate the cluster with the largest area and choose that as the cluster to go forward with
        MAXClust <- which.max(do.call(rbind, lapply(ClustStore, function(i)colSums(i['FieldArea']))))
        cluster <-  ClustStore[[MAXClust]]
      }
      
      ## Label the fields in my opportunity canvas that are part of the cluster in this round
      ## Can use the loop counter to label the random sample number
      ## Then also add on a label for each cluster and the distance of each field to the centre of the cluster
      cluster <- cluster |>  arrange(ParcRef) # add distance to center of cluster # |> mutate(CentDist = as.numeric(st_distance(st_centroid(st_union(geometry)), geometry))) 
      clusterInd <- which(OppFieldsLoop$ParcRef %in% cluster$ParcRef)
      OppFieldsLoop$RandSamp[clusterInd] <- k
      OppFieldsLoop$ClustID[clusterInd] <- ClusterCount; ClusterCount <- ClusterCount+1
      # OppFieldsLoop$CentDist[clusterInd] <- cluster$CentDist
      
      ## Calculate the total area of all clusters in this round of for loop
      TotClustArea <- sum((OppFieldsLoop |> filter(RandSamp==k))$FieldArea)
      
      ## Now that a cluster has been chosen I can reset the while loop counter for the creation of the next cluster
      WhileCount = 1
      ClustStore <- vector(mode = "list")
      
        
      } ## 1st WHILE LOOP END -- stop creating clusters
      
      
      ## Now I have created a series of clusters of the suitable total size
      ## I need to now label these in the full data set and not just the smaller data set used in the loop
      clusterAll <- which(OppFields$ParcRef %in% (OppFieldsLoop |> filter(RandSamp==k))$ParcRef)
      OppFields$RandSamp[clusterAll] <- k
      OppFields$ClustID[clusterAll] <- (OppFieldsLoop |> filter(RandSamp==k))$ClustID
      # OppFields$CentDist[clusterAll] <- (OppFieldsLoop |> filter(RandSamp==k))$CentDist
    
      ## Calculate the shortfall in FieldArea for this round of for loop, can add this on next time (stops last sample being very large/and carryover shortfall to next iteration of loop)
      if(Additive==FALSE){NextGapExtra <- (AreaSeg[k]*k) - sum((OppFields |> filter(is.na(RandSamp)==F))$FieldArea)}
      
      ## Assign the parcel refs chosen in this round of sampling to the list
      ## Also assign the ClustID to a list so each parc ref can be matched to a cluster ID
      listSamples[[k]] <- (OppFieldsLoop |> filter(RandSamp==k))$ParcRef
      listClustID[[k]] <- (OppFieldsLoop |> filter(RandSamp==k))$ClustID
      ## If want to create reserve and AES also calculate whether fields are in the center or periphary of clusters and assign these to a list
      if(length(NewCat)==2){
        
          FieldPost <- OppFieldsLoop |> 
                        filter(RandSamp==k) |> 
                        group_by(ClustID) |> 
                        mutate(CentDist = as.numeric(st_distance(st_centroid(st_union(geometry)), geometry)),
                               ClustPosit = ifelse(CentDist >= median(CentDist), "Periph", "Cent")) |> 
                        ungroup()
      
          listPosit[[k]] <- FieldPost$ClustPosit
      }
      
    }# end of loop k 
    
    ## see the total area of each set of fields
    ## I can hash this out if I want to save on speed for a final run 
    # OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))
    # ggplot(data = OppFields) + geom_sf(mapping = aes(geometry = geometry, fill = CentDist), colour = NA)
    
    ## Select just the unique field ID column and random set number for a join
    if(Additive==FALSE){OppFields <- OppFields |> select(ParcRef, RandSamp, ClustID, CentDist) |> st_drop_geometry()}
    
    # ## Join set numbers onto my main data set
    # Canvas2 <- left_join(as.data.frame(Canvas), OppFields, by = "ParcRef")
    # stopifnot(nrow(Canvas) == nrow(Canvas2)) # check nothing went wrong during join
    # Canvas <- Canvas2 |> st_as_sf(); rm(Canvas2) # revert back to original name
    
    ## Here based off the distance to the center of the cluster work out which fields are in the center and which on the periphery
    ## Only need to do this if the scenario involves creating reserves and AES only sites
    # if(length(NewCat)==2){
    #   Canvas <- Canvas |> 
    #             group_by(ClustID) |> 
    #             mutate(ClustPosit = ifelse(CentDist > median(CentDist), "Periph", "Cent")) |> 
    #             ungroup()
    #   }
        
  } # end of cluster strategy labeling
  
  
  
  ##--------------------------------------------##    
  #### F 1.3 Label fields for random strategy ####
  ##--------------------------------------------##    
  
  if(ScenType == "random"){
    
    ## Message to update console
    message("Labelling segments randomly...")
    
    ## filter out the fields that are within the opportunity area for this scenario
    OppFields <- Canvas[IndOpp,] |> select(c("FieldArea", "ParcRef", "ClustDist")) |> select(all_of(c("FieldArea", "ParcRef", "ClustDist")))
  
    ## Add on columns for sample numbering and columns for scaled distance to nearest cluster and the inverse
    OppFields <- OppFields |> mutate(RandSamp=NA,
                                     InvClustDist = (1/(ClustDist+10))/mean(1/(ClustDist+10)),
                                     PClustDist= ((ClustDist+10))/mean((ClustDist+10)))
    
    # sample the rows of OppFields, this mixes up the order and allows for random sampling
    # for bigger weight sampling by inverse distance to wader population
    # for more weight sampling by distance to wader population
    set.seed(1212)
    if(Strategy == "Big"){  my_samp <- sample(1:nrow(OppFields), replace = FALSE, 
                                                     prob = OppFields$InvClustDist) }
    if(Strategy == "More"){  my_samp <- sample(1:nrow(OppFields),  replace = FALSE, 
                                                     prob = OppFields$PClustDist) }
    if(Strategy == "Better"){ my_samp <- sample(1:nrow(OppFields), replace = FALSE) }
    
    ## Create an empty list where I can put the parcel references that belong to eacg round of sampling
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
      
      ## recreate my samp with a different order for additive sample so samples could overlap
      if(Additive==TRUE){
          if(Strategy == "Big"){  my_samp <- sample(1:nrow(OppFields), replace = FALSE, prob = OppFields$InvClustDist) }
          if(Strategy == "More"){  my_samp <- sample(1:nrow(OppFields),  replace = FALSE, prob = OppFields$PClustDist) }
          if(Strategy == "Better"){ my_samp <- sample(1:nrow(OppFields), replace = FALSE) }
          }
      } # end of sampling loop (j)
    
    
    ## see the total area of each set of fields
    if(Additive==FALSE){OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))} 
    
    ## Don't need this anymore as use `listSamples` to extract samples and not the RandSamp column
    # ## Select just the unique field ID column and random set number for a join
    # OppFields <- OppFields |> select(ParcRef, RandSamp) |> st_drop_geometry()
    # 
    # ## Join set numbers onto my main data set
    # Canvas2 <- left_join(as.data.frame(Canvas), OppFields, by = "ParcRef")
    # stopifnot(nrow(Canvas) == nrow(Canvas2)) # check nothing went wrong during join
    # Canvas <- Canvas2 |> st_as_sf(); rm(Canvas2) # revert back to original name
    # table(Canvas$RandSamp)
    
  } # end of "random" strategy labeling
  
  
  
  ##--------------------------------------------------##    
  #### F 1.4 Label fields for stakeholder1 strategy ####
  ##--------------------------------------------------##    
  
  if(ScenType == "stakeholder1"){
    
    ##** If I just want to do the stakeholder scenarios at the field level then I can mainly reuse this code
    ##** I would order OppFields by the grading (highest first) and then run the which(cumsum()) section on the ordered data set, taking the highest grades first
    ##** I could also maybe add in another weighting for distance from cluster.e.g. more grading x  PClustDist
    ##** Alternatively I could copy the approach above and provide a smoothed raster of grading (smoothed by cluster size)
    ##** Then in this line of code       if(Strategy == "Big"){  initial_polygon <- sample_n(OppFieldsLoop2, 1, weight = OppFieldsLoop2$InvClustDist) }
    ##** I would just change OppFieldsLoop2$InvClustDist for my bigger grade or multiple the grade by InvClustDist
    
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

    ## Old version which essentially just samples in order, I could actually use this for the mutliplicative sampling
    # # depending on the strategy being used order OppFields by a different grading
    # set.seed(1212)
    # if(Strategy == "Big" & !any(OppCat %in% "Arable Opp")){  
    #   OppFieldsOrd <- OppFields |> dplyr::arrange(across(paste0("Bigger_", StakeGroup), desc))
    # }
    # if(Strategy == "More" & !any(OppCat %in% "Arable Opp")){  
    #   OppFieldsOrd <- OppFields |> dplyr::arrange(across(paste0("More_", StakeGroup), desc))
    # }
    # if(Strategy == "Better" & !any(OppCat %in% "Arable Opp")){ 
    #   OppFieldsOrd <- OppFields |> dplyr::arrange(across(paste0("Better_", StakeGroup), desc))
    # }
    # if(any(OppCat %in% "Arable Opp")){ 
    #   OppFieldsOrd <- OppFields |> dplyr::arrange(across(paste0("Arable_", StakeGroup), desc))
    # }
    # 
    # 
    # ## Create an empty list where I can put the parcel references that belong to eacg round of sampling
    # listSamples <- vector(mode = "list", length = N_sets)
    # 
    # ## Starting point, it is a counter that keep up to date with what fields have already been labelled in a cluster
    # Start_point <- 1
    # 
    # ## Label opportunity fields for each round of sampling
    # for(j in 1:N_sets){
    # 
    #   ## On first loop set the shortfall to 0 and also reset if for each round of Addtive sampling
    #   if(j==1 | Additive==TRUE){NextGapExtra <- 0}
    # 
    #   ## return the number of sampled rows with cumulative sum < condition
    #   sampRows <- length(which(cumsum(OppFieldsOrd$FieldArea[Start_point:nrow(OppFieldsOrd)]) < (AreaSeg+NextGapExtra)))
    #   OppFieldsOrd$RandSamp[Start_point:(Start_point+sampRows-1)] <- j # label the sample with loop counter number
    # 
    #   ## Calculate the shortfall in FieldArea, can add this on next time (stops last sample being very large/and carryover shorfall to next iteraiton of loop)
    #   NextGapExtra <- ((AreaSeg*j) - sum(OppFieldsOrd$FieldArea[is.na(OppFieldsOrd$RandSamp)==F]))
    #   
    #   ## update the next start point, this will start at the first NA field
    #   Start_point <- Start_point+sampRows
    #   
    #   } # end of sampling loop (j)
    
    # ## see the total area of each set of fields
    # OppFieldsOrd |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))
    # 
    # ## Select just the unique field ID column and random set number for a join
    # OppFieldsOrd <- OppFieldsOrd |> select(ParcRef, RandSamp) |> st_drop_geometry()
    # 
    # ## Join set numbers onto my main data set
    # Canvas2 <- left_join(as.data.frame(Canvas), OppFieldsOrd, by = "ParcRef")
    # stopifnot(nrow(Canvas) == nrow(Canvas2)) # check nothing went wrong during join
    # Canvas <- Canvas2 |> st_as_sf(); rm(Canvas2) # revert back to original name
    # # ggplot(data = Canvas) + geom_sf(mapping = aes(geometry = geometry, fill = RandSamp), colour = NA) + scale_fill_viridis_c()
    
    
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
  
  ## Pre-compute buffer for all fields so I just have to do it once
  CanvCants <- Canvas |> st_centroid() |> select(geometry)
  Parc500 <- CanvCants |> st_buffer(dist = 500)
  Parc1000 <- CanvCants |> st_buffer(dist = 1000)
  Parc1500 <- CanvCants |> st_buffer(dist = 1500)
  Parc2000 <- CanvCants |> st_buffer(dist = 2000)
  
  ## Store a copy of the unaltered canvas here
  ## Can then call on this to reset the canvas at the start of everyloop
  MASTERCanv <- Canvas

  
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
    
            
    ## Get the row numbers of the opportunity fields in the current random sampling group
    ## Update here so that sampled fields are drawn from a list of parcel refs and not from a column
    # SetInd <- which(Canvas$RandSamp==i)
    SetInd <- which(Canvas$ParcRef %in% listSamples[[i]])
    
    ## define which within habitat categories I am going to update
    ## also update the categopry for that row
    HabCats <- c("Category", "GrassOpp", "Opp", "GRASSLAND_TYPE", "WaterCoverage", "STOCK", "VEG_STRUCTURE", "RUSH_PERCENT", "Fence_Coverage")
    
    ## Assign the total area of this sampling set to the data set that tracks progress
    ScTracker$SegmentArea[i+1] <- sum(Canvas$FieldArea[SetInd])
  
    ## Get the row numbers of fields with the management option that new opportunity fields will become
    ## for plus sampling remove fields that do not have waders in
    if(Plus==F){IndNewHabs <- which(Canvas$Category %in% NewCat)}
    if(Plus==T){IndNewHabs <- which(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)}
    
    ## Sample a set of fields from the chosen management option equal to the number of fields in the random opportunty sample
    ## If the set if chosen management fields is smaller then the opportunity fields that are going to be changed 
    ## then use replace = TRUE or there will not be enough options in the random sample
    set.seed(1212)
    if(length(NewCat)==1){
        
        # does sampling need to be done with replacement of not (depends on required length of sample)
        if(length(IndNewHabs)>length(SetInd)){ReplaceSet = FALSE}else{ReplaceSet = TRUE} 
        NewHabs <- sample(x = IndNewHabs, size = length(SetInd), replace = ReplaceSet)
        
        ## Now assign new habitat values to the opportunity fields
        Canvas[SetInd, colnames(Canvas) %in% HabCats] <- Canvas[NewHabs, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
        print(table(Canvas$Category))
        
      }
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
      Canvas[SetPeriph, colnames(Canvas) %in% HabCats] <- Canvas[NewHabsP, colnames(Canvas) %in% HabCats] |> st_drop_geometry()
      print(table(Canvas$Category))
        
      }
    }


    
    ##------------------------------------##    
    #### F 1.6 Update landscape habitat ####
    ##------------------------------------##
        
    ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
    ## The base raster is just the UKCEH land cover data set at 25m res
    ## First create index of rows for fields that are suitable
    Suit <- which(Canvas$Opp == "Grass" | Canvas$WiderWetGrass == 1)
    # Suit <- which(Canvas$Category %in% c("Reserve", "AES Only", "Grass Opp"))
    WaterRast <- rasterize_polygons(Canvas[Suit, ], LandCov, min_coverage = 0.25)
    
    ### Now assign the pixels a water coverage value from the correct field
    values(WaterRast) <- Canvas[Suit,]$WaterCoverage[values(WaterRast)]
    # plot(WaterRast)
    
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
    if(any(OppCat %in% "Arable Opp")){
      
      ## Rasterize my field polygons, returns the index of the polygons that overlaps most of the pixel
      ## The base raster is just the UKCEH land cover data set at 25m res
      ## First create index of rows for fields that are suitable
      Suit <- which(Canvas$Opp == "Grass" | Canvas$WiderWetGrass == 1)
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
    
    ## Predict wader abundance and then assign it to the correct abundance column
    ## If there is no model object provided then do not predict for that species
    ## Predict for Snipe
    if(is.null(SniModel)==F){
        Canvas$SnipeAbund[Index] <- (predict.rfsrc(object = SniModel, 
                              newdata = (Canvas[Index, ] |> select(SniModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                              jitt=FALSE))$predicted
        ScTracker$Snipe[i+1] <- sum(Canvas$SnipeAbund, na.rm = T)
        } # number of Snipe in the landscape
      
    
    ## Predict for Lapwing
    if(is.null(LapModel)==F){
        Canvas$LapAbund[Index] <- (predict.rfsrc(object = LapModel, 
                                  newdata = (Canvas[Index, ] |> select(LapModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                  jitt=FALSE))$predicted
        ScTracker$Lapwing[i+1] <- sum(Canvas$LapAbund, na.rm = T)
        } # number of Lapwing in the landscape 
      
    
    ## Predict for Redshank
    if(is.null(RedModel)==F){
        Canvas$RedAbund[Index] <- (predict.rfsrc(object = RedModel, 
                                  newdata = (Canvas[Index, ] |> select(RedModel[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                                  jitt=FALSE))$predicted
        ScTracker$Redshank[i+1] <- sum(Canvas$RedAbund, na.rm = T)
        } # number of Lapwing in the landscape 
  
    
    
    ##------------------------------------------##
    ##------------------------------------------##
    #### F 1.8 Save current iteration outputs ####
    ##------------------------------------------##
    ##------------------------------------------##
    
    ## Read out shape file of current scenario state
    ## Might want to only save a subset of columns that might be important for plotting
    ## Alternatively could just save the final output as a shapefile
    if(SaveCanv==T){write_sf(Canvas, paste0(Outpath, "Canvas_", Strategy, "_", paste(NewCat,collapse="&"), if(Plus==T){paste0("_plus")}, if(any(OppCat %in% "Arable Opp")){"_Arable"}, i, ".shp"), append=FALSE)}

    ## Update how long it took to run this segment
    toc(log = TRUE, quiet = TRUE)
    ScTracker$SegTime[i+1] <- paste0(tic.log(format = T))
    tic.clearlog()
    write_csv(ScTracker, paste0(Outpath, "Tracker_", Strategy, "_", paste(NewCat,collapse="&"), if(Plus==T){paste0("_plus")}, if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
    
    ## return the tracker data set
    print(ScTracker)
    
    ## If doing the additive testing, i.e. only can look at additive effects
    ## Then reset the canvas at this point, this will mean that the in each round of the i loop we go back to the starting canvas
    ## and do not start from where we left off from on the last loop
    if(Additive==TRUE){Canvas <- MASTERCanv}
    if(Additive==FALSE){Canvas <- Canvas}
    
    
  } # End of segment loop (i)
    
  ## Add on column for cumulative sum of area
  ScTracker <- mutate(ScTracker, CumArea = cumsum(SegmentArea))
  write_csv(ScTracker, paste0(Outpath, 
                              "Tracker_", 
                              if(Additive==TRUE){paste0("Additive")},
                              if(Additive==FALSE){paste0("Multiplicative")},
                              Strategy, "_", 
                              paste(NewCat,collapse="&"), 
                              if(Plus==T){paste0("_plus")}, 
                              if(any(OppCat %in% "Arable Opp")){"_Arable"}, ".csv")) # save the tracker log
  
} 
#### F 1 END ####



