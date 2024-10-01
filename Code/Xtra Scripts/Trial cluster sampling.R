# ## Read in the various reserve outlines
# ## RSPB reserves
# RSPB <- st_read("RawData/RSPB Reserves/EnglandWales_RSPBReserves.shp") |> st_crop(Canv)
# plot(RSPB$geometry, border = 'blue', col = alpha('blue', 0.5))
# mean(st_area(RSPB$geometry)/10000)
# ## Local Nature reserves
# LNR <- st_read("RawData/Other Reserves/Local_Nature_Reserves_EnglandPolygon.shp") |> st_crop(Canv)
# ## National Nature Reserves
# NNR <- st_read("RawData/Other Reserves/National_Nature_Reserves_EnglandPolygon.shp") |> st_crop(Canv)
# plot(NNR$geometry,  col = NNR$nnr_name)
# ggplot() + geom_sf(data=NNR, mapping = aes(geometry=geometry, fill = nnr_name))



# Load the necessary libraries
library(sf)
library(terra)
library(dplyr)

## set variables for following code
Canvas=Canv
SniModel=SniMod
LapModel=LapMod
RedModel=NULL
N_sets=5
OppCat= c("Grass Opp","AES Only")
NewCat= c("Reserve")
Plus = F
Strategy = "More"
PropOppArea = 1/3
Outpath = "CleanData/Scenarios/5-SomScenario/Rand/"
SaveCanv = F



####-- Scenario set up --####

## The number segments to split scenario modeling up into
## For each segment habitat variables are updated and wader abundance is calculated
N_sets <- N_sets

## Create data frame to record wader abundance at the end of each segment 
## Also record the starting abundance in a seperate columns
ScTracker <- data.frame(Landscape = Canvas$Landscape[1],
                        Segment = 1:(N_sets+1), 
                        SegmentArea = NA,
                        SegTime = NA,
                        Snipe = NA, 
                        BaseSnipe = ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T)),
                        Lapwing = NA, 
                        BaseLapwing = ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T)),
                        Redshank = NA,
                        BaseRedshank = ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T)),
                        OppCat= paste(OppCat,collapse="&"),
                        NewCat=  paste(NewCat,collapse="&"), 
                        Plus = Plus, 
                        Strategy = Strategy, 
                        PropOppArea = PropOppArea)

## Set the first row as zero area and baseline abundance
ScTracker$Snipe[1] <- ifelse(sum(Canvas$SnipeAbund, na.rm = T)==0, NA, sum(Canvas$SnipeAbund, na.rm = T))
ScTracker$Lapwing[1] <- ifelse(sum(Canvas$LapAbund, na.rm = T)==0, NA, sum(Canvas$LapAbund, na.rm = T))
ScTracker$Redshank[1] <- ifelse(sum(Canvas$RedAbund, na.rm = T)==0, NA, sum(Canvas$RedAbund, na.rm = T))
ScTracker$SegmentArea[1] <- 0

  ## Get the index of rows representing the opportunity area for this scenario
  ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
  if(Strategy %in% c("Big", "More") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat) }
  if(Strategy %in% c("Better") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat) }
  
  ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders 
  ## (i.e. AES fields with waders will not be targeted and those without will be)
  if(Strategy %in% c("Big", "More") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCat & !(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)) }
  if(Strategy %in% c("Better") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCat & !(Canvas$Category %in% NewCat & Canvas$Tot_abund > 0)) }
  
  ## Calculate the total field Area for my opportunity Area
  TotArea <- sum(Canvas[IndOpp,]$FieldArea)*PropOppArea
  
  ## Split up that area into segments
  AreaSeg <- TotArea/N_sets
  
  ## filter out the fields that are within the opportunity area for this scenario
  OppFields <- Canvas[IndOpp,]
  OppFields$RandSamp <- NA # column to ID fields in different segments
  
  OppFields$InvClustDist <- (1/(OppFields$ClustDist+10))/mean(1/(OppFields$ClustDist+10)) ## **NEW code could update original function
  OppFields$PClustDist <-((OppFields$ClustDist+10))/mean((OppFields$ClustDist+10))  ## **NEW code
  plot(OppFields$geometry)
  
  # sample the rows of OppFields, this mixes up the order and allows for random sampling
  # for bigger weight sampling by inverse distance to wader population
  # for more weight sampling by distance to wader population
  # set.seed(1212)
  # if(Strategy == "Big"){  my_samp <- sample(1:nrow(OppFields), replace = FALSE, 
  #                                                  prob = (1/(OppFields$ClustDist+10))/mean(1/(OppFields$ClustDist+10)))   }
  # if(Strategy == "More"){  my_samp <- sample(1:nrow(OppFields),  replace = FALSE, 
  #                                                  prob = ((OppFields$ClustDist+10))/mean((OppFields$ClustDist+10)))   }
  # if(Strategy == "Better"){ my_samp <- sample(1:nrow(OppFields), replace = FALSE) }

## CODE IS SAME AS BEFORE UP TO HERE ##


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

  
## Label opportunity fields for each round of sampling
for(j in 1:N_sets){ 
  
  
  ## Message for progress of cluster sampling
  message("Createing set of clusters ", j , " out of ", N_sets)
  
  ## set the total area cluster for this round of the loop to zero
  ## This only resets outside of all the while loops
  TotClustArea <- 0
  
  ## On first loop set the shortfall to 0
  if(j==1){NextGapExtra <- 0}
  ## On the final loop set the shortfall as a large number so all fields are sampled
  if(j==N_sets & PropOppArea == 1){NextGapExtra <- sum(OppFields$FieldArea)+1}
  
  ## Filter out the fields that are still in the opportunity area
  ## For each iteration of the loop this will gradually decrease the size of the data set I am working with
  OppFieldsLoop <- OppFields |> filter(is.na(RandSamp)) |> select(c("RandSamp", "InvClustDist", "FieldArea", "ParcRef"))
  
  
  ##-- 1st WHILE LOOP START --##
  ## While total area of all clusters created in this iteration of the while function is < AreaSeg, keep on creating more clusters
  ## Need to add on NextGapExtra to make sure sampling isn't too much or too little (depending on if previous loops over/undersampled)
  while(TotClustArea < (AreaSeg+NextGapExtra)) {
    

  # Generate a random max threshold in hectares for current cluster
  # Currently 345 is average RSPB reserve in Somerset
  # 69 is average farm size in SW England: https://assets.publishing.service.gov.uk/media/6627dac8838212a903a7e6d9/regional-profiles-stats-region-south_west-25apr24.pdf
  threshold_area <- rtruncnorm(n=1, a=0, mean = 345, sd = 50)  # Example: between 10 and 50 hectares
  
  ## If this cluster exceeds the threshold then break the while loop so that it is not created
  ## For the last loop allow the last cluster to go over the threshold
  if((TotClustArea+threshold_area) > (AreaSeg+NextGapExtra) & !j==N_sets) break

  # Select an initial polygon as the initiation point of the cluster
  # Weight the random draw by the inverse distance to cluster i.e. largest values closest to cluster
  # remove fields already in cluster from earlier iterations of this while loop
  OppFields2 <- OppFieldsLoop |> filter(is.na(RandSamp)) # remove any fields already part of a cluster
  initial_polygon <- sample_n(OppFields2, 1, weight = OppFields2$InvClustDist) # or PClustDist
  
  # Create the initial cluster
  cluster <- initial_polygon
  cluster_area <- sum(cluster$FieldArea)  # Area in hectares
  
  # Create data set of potential fields for this cluster
  PotentialRefs <- (OppFieldsLoop |> filter(is.na(RandSamp)) |> st_crop(st_buffer(cluster, dist = 7500)))$ParcRef 
  ClustOppFields <- filter(OppFieldsLoop, ParcRef %in% c(PotentialRefs))
  # plot(ClustOppFields$geometry)
  
  
  
  ##-- 2nd WHILE LOOP START --##
  message("making a cluster")
  # Expand the cluster while the cluster size is below the max threshold (threshold_area)
  while(cluster_area < threshold_area) {
    
    
    ## Retrieve potential neighbors
    neighbors <- get_neighbors(cluster, ClustOppFields, clustbuff = 100)
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
      selected_neighbor <- neighbors[which(cumsum((neighbors$FieldArea)) < (threshold_area-cluster_area)),] |> select(-CentDist)
      if(nrow(selected_neighbor) == 0) break # break if no selected_neighbor, means a single fields takes you over the threshold 
      # if(nrow(selected_neighbor) == 0){selected_neighbor <- neighbors[1,] |> select(-CentDist)} # just use first fields if none selected
      
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
  
  ## Label the fields in my opportunity canvas that are part of the cluster in this round
  ## Can use the loop counter to do the labeling
  ## Go back to full data set here to do the labelling
  clusterInd <- which(OppFieldsLoop$ParcRef %in% cluster$ParcRef)
  OppFieldsLoop$RandSamp[clusterInd] <- j
  
  ## Calculate the total area of all clusters in this round of for loop
  TotClustArea <- sum((OppFieldsLoop |> filter(RandSamp==j))$FieldArea)
  
    
  } ## 1st WHILE LOOP END -- stop expanding single cluster
  
  
  ## Now I have created a series of clusters of the suitable total size
  ## I need to now label these in the full data set and not just the smaller data set used in the loop
  clusterAll <- which(OppFields$ParcRef %in% (OppFieldsLoop |> filter(RandSamp==j))$ParcRef)
  OppFields$RandSamp[clusterAll] <- j

  ## Calculate the shortfall in FieldArea for this round of for loop, can add this on next time (stops last sample being very large/and carryover shortfall to next iteration of loop)
  NextGapExtra <- (AreaSeg*j) - sum((OppFields |> filter(is.na(RandSamp)==F))$FieldArea)
  
  
} # end of sampling loop (j)


## Plot the final product with clusters colored according the the round of the loop they were created
ggplot(data = OppFields) + geom_sf(mapping = aes(geometry = geometry, fill = RandSamp), colour = NA)


## see the total area of each set of fields
OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea)) 


"1644.542  1243.405  1271.711  1143.550  1456.729 13033.899"

  

   
  #   
  # 
  # ####-- Label fields in each segment of scenario --####
  # 
  # ## Message to update console
  # message("Labelling segments....")
  # 
  # 
  # ## Label opportunity fields for each round of sampling
  # for(j in 1:N_sets){
  #   
  #   ## On first loop set the shortfall to 0
  #   if(j==1){NextGapExtra <- 0}
  #   ## On the final loop set the shortfall as a large number so all fields are sampled
  #   if(j==N_sets & Strategy == "Better"){NextGapExtra <- sum(OppFields$FieldArea)+1}
  #   
  #   ## return the sampled rows with cumulative sum < condition
  #   sampRows <- my_samp[which(cumsum(OppFields$FieldArea[my_samp]) < (AreaSeg+NextGapExtra))]  
  #   OppFields$RandSamp[sampRows] <- j # lable the sample with loop counter number
  #   
  #   ## Calculate the shortfall in FieldArea, can add this on next time (stops last sample being very large/and carryover shorfall to next iteraiton of loop)
  #   NextGapExtra <- ((AreaSeg*j) - sum(OppFields$FieldArea[is.na(OppFields$RandSamp)==F]))
  #   
  #   ## update the sample of rows to remove the rows number just sampled
  #   my_samp <- my_samp[(length(sampRows)+1):length(my_samp)]} # end of sampling loop (j)
  # 
  # 
  # ## see the total area of each set of fields
  # OppFields |> group_by(RandSamp) |> summarise(Total = sum(FieldArea))
  # 
  # ## Select just the unique field ID column and random set number for a join
  # OppFields <- OppFields |> select(ParcRef, RandSamp) |> st_drop_geometry()
  # 
  # ## Join set numbers onto my main data set
  # Canvas2 <- left_join(as.data.frame(Canvas), OppFields, by = "ParcRef")
  # stopifnot(nrow(Canvas) == nrow(Canvas2)) # check nothing went wrong during join
  # Canvas <- Canvas2 |> st_as_sf(); rm(Canvas2) # revert back to original name
  