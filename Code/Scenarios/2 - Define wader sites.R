##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 15/03/24
## Goal: Define breeding wader sites in each region
##
##
##------------------------------------------------------##

## Overview of Workflow
## 1. K-means clustering within landscapes
## 2. Kernel density estimates within clusters
## 3. Remove small cluster (based on both # occupied fields | # breeding pairs)
## 4. Join close adjacent clusters


## Load in required packages
pacman::p_load(here, sf, tidyverse, data.table, factoextra, amt)
options(scipen = 100, digits = 4) # set notation
source("Code/Helper functions.R") # get helper functions

## Use `here` package for specifying file paths
here::i_am("Code/Scenarios/2 - Define wader sites.R")
library(here)
here()




##----------------------##
#### 0.1 Data read in ####
##----------------------##

## Read in filtered breeding pairs estimates
Waders <- read_csv("CleanData/Wader Abundance/4-AddLandscapeAttributes/Breeding_Pairs_FullAttrib3.csv")
colnames(Waders)

## Read in shapefile of all BWWM fields and then rename to make a data set to label fields in different clusters
Fields <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp")) |> select(F_LOC_ID)
FieldsClusts <- Fields |> mutate(OrigArea = st_area(geometry))


##--------------------------##
#### 0.2 Define Functions ####
##--------------------------##

## k-means clustering loop function
## function that loops through different values of k and calculates within group sum of squares
wssplot <- function(data, vline=7, nc=25, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares",
       cex.lab=1.5, cex.axis=1.25)
  abline(v = vline, col="darkgrey", lwd=3, lty=2)
  wss
}




##-------------------------------##
#### 1.0 Define Somerset Sites ####
##-------------------------------##

##------------------------------------##
#### 1.1 Extract fields with waders ####
##------------------------------------##

## Filter out the field data for just Somerset
Somerset <- filter(Waders, Landscape == "Somerset Levels and Moors")
FieldLength <- sqrt(median(Somerset$FieldArea))

## Work out which fields have at least one pair of the target species (Snipe, Redshank, Lapwing)
Somerset <- Somerset %>%  
  mutate(Tot_abund = rowSums(dplyr::select(., est_pairsL, est_pairsR, est_pairsS), na.rm = T),
         Wader_Pres = ifelse(Tot_abund>0, "Y", ifelse(is.na(Tot_abund)== T, "N", "N"))) |> 
  select(-Tot_abund)

## Filter out the fields with at least one pair of target wader species present
Somerset_Waders <- filter(Somerset, Wader_Pres == "Y")

## plot the fields that have at least one pair of waders
ggplot() +
  geom_point(data = Somerset_Waders, aes(x = FieldX, y= FieldY)) +
  theme_light()

## Save spatial data set to indicate which fields have waders
SomersetPres <- Somerset |> select(F_LOC_ID, year, est_pairsL, est_pairsR, est_pairsS, N_visits, Wader_Pres)
SomersetPres <- inner_join(SomersetPres, Fields, by = "F_LOC_ID")
write_sf(SomersetPres, "CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Present.shp")



##----------------------------##
#### 1.2 K-means clustering ####
##----------------------------##

## Get the X and Y coordinates for the fields with at least one pair of waders
Somerset_WadersFields <- Somerset_Waders |> select(FieldX, FieldY)

## Apply function to Somerset field centers
png(filename= "CleanData/Scenarios/2-DefineWaderSites/Example Pics/Som_K_vs_SumofSquares.png", units = "px", height = 475, width = 600)
wssplot(data = Somerset_WadersFields, vline=7)
dev.off()

## Run k-means clustering with chosen values for k
kmeans_out <- kmeans(x = Somerset_WadersFields, centers = 7, iter.max = 25, nstart = 25)
Somerset_Waders$cluster <- kmeans_out$cluster

## Plot the clusters using ggplot and the factoextra function fviz_cluster()
ggplot() +
  geom_point(data = Somerset_Waders, aes(x = FieldX, y= FieldY, colour = as.factor(cluster))) +
  theme_light()

fviz_cluster(kmeans_out, data = Somerset_WadersFields,
             xlab = "X coordinate",
             ylab = "Y coordinate",
             legend.title = "Cluster ID",
             main = " ",
             labelsize =0,
             ggtheme = theme_light(),
             font.x = 18, font.y = 18,
             font.legend = 15, font.tickslab = 15)
# save the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/2-DefineWaderSites/Example Pics/Som_KmeansClusts.png", units = "in", height = 8, width = 11)
  

## select the columns from the main data set that I want to keep
SomersetCLust <- Somerset_Waders |> select(F_LOC_ID, est_pairsL, est_pairsR, est_pairsS, cluster, FieldX, FieldY)




##------------------------------------##
#### 1.3 Calculate KDE's by cluster ####
##------------------------------------##

## Use `amt` functions to calculate KDE for each cluster
SomersetCLust <- as.data.frame(SomersetCLust)
Somtrk <- make_track(tbl = SomersetCLust, .x = FieldX, .y = FieldY,
                     ID = cluster, crs = 27700)

## grouping factor to calculate KDE
groupfac <- "ID"

## make template raster for KDE's
trast1 <- make_trast(Somtrk, res = 25) # res for rastet that kde's calculated on

## make KDE per grouping factor
## use average field with as bandwidth for kernel density estimate
hr_Som <- Somtrk %>%  
  nest(data = -groupfac) %>%  
  mutate(kde = map(data, hr_kde, trast = trast1, levels = 0.95, h = c(FieldLength, FieldLength)))



##---------------------##
#### 1.4. Plot KDE's ####
##---------------------##

## Extract the polygon for the 95% KDE's from each cluster
for(j in 1:nrow(hr_Som)){
  
  if(j == 1){Areas <- hr_isopleths(hr_Som$kde[[j]])}else{Areas <- rbind(Areas,hr_isopleths(hr_Som$kde[[j]]))}
  
}

## make plot of all 95% KDE's
ggplot() + 
    
    ## Add home ranges
    geom_sf(data = Areas, aes(fill = "blue"), alpha = 0.5) +
    scale_fill_manual(name = expression(""), values =c("blue", "blue"), labels = c("95% kde")) +
    
    ## change labels
    labs(x = "Longitude", y = "Latitude") +
    
    # set theme
    theme_light() +
    theme(legend.box.background=element_rect(colour = "#BDC3C7"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text=element_text(colour="black"),
          axis.title.x = element_text(size = 17),
          axis.text.x = element_text(hjust=0.7, size = 13,),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 17),
          axis.text.y = element_text(size = 13, hjust=0.7, angle=90, vjust=0.3))

# save the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/2-DefineWaderSites/Example Pics/Som_KDE_clusters.png", units = "in", height = 9, width = 10)
  
  

##-------------------------------##
#### 1.5 Remove small clusters ####
##-------------------------------##

## Separate mutliple polygons into distinct polygons (some 95% KDE's were not single polygons)
for(p in 1:nrow(Areas)){
  
  Inter <- as.data.frame(st_cast(st_sfc(Areas$geometry[p]), "POLYGON"))
  Inter$ClustID <- hr_Som$ID[p]
  
  if(p == 1){Areas_poly <- Inter}else{Areas_poly <- rbind(Areas_poly, Inter)}

}

## Assign each polygon a distinct number
Areas_poly <- st_as_sf(Areas_poly)
Areas_poly$ShapeID <- 1:nrow(Areas_poly)

## Calculate the centroid for each field with waders
SomersetPres_Cords <- filter(SomersetPres, Wader_Pres == "Y")
Somerset_Cent <- st_centroid(st_as_sf(SomersetPres_Cords))

## find points within polygons
Field_in_Clust <- st_join(Somerset_Cent, st_as_sf(Areas_poly), join = st_within)

## Calculate the number of fields with waders and wader pairs within each shape
## Remove shapes that do not have at least three occupied fields or three pairs of any species
Areas_poly <- Field_in_Clust |> 
  group_by(ShapeID) |> 
  summarise(n_fields = n(),
            n_lap = sum(est_pairsL, na.rm = T),
            n_red = sum(est_pairsR, na.rm = T),
            n_snipe = sum(est_pairsS, na.rm = T),
            n_wad = sum(n_lap, n_red, n_snipe)) |> 
  st_drop_geometry() |> 
  filter(n_fields >= 3 | n_wad >= 3) |> 
  select(ShapeID, n_wad) |> 
  inner_join(Areas_poly, by = "ShapeID") |> 
  st_as_sf()

## plot the different groupings
ggplot(data = Areas_poly) + geom_sf(aes(fill = n_wad)) + theme_light()



##-----------------------------------------##
#### 1.6 Join close/overlapping clusters ####
##-----------------------------------------##

## Loop through each cluster/shape and group those that are closer than Field Length value
counter <- 1 # counter to keep track of the most recent cluster group ID number
Areas_poly$ClustGroup <- NA # initiate column where cluster group ID's go

for(j in 1:nrow(Areas_poly)){
  
  message("iteration ", j)
  
  ## Calculate distances to all other features
  index <- st_distance(Areas_poly[j,], Areas_poly[-j,])
  
  mindist <- as.numeric(min(index))
  
  ## if no cluster within the threshold distance then assign current cluster next cluster number
  if(mindist > 2*FieldLength){Areas_poly$ClustGroup[j] <- counter; counter <- counter + 1}
  
  ## if there are cluster within the threshold distance then they must be numbered together
  if(mindist < 2*FieldLength){
    
    ## Check if any cluster identified in this loop already have a cluster ID
    PrevClust <- max(Areas_poly$ClustGroup[which(as.numeric(index) < 2*FieldLength)], na.rm = T)
    
    ## If there was already a ClustGroup ID for any of the clusters below the threshold distance then use that number to label the current clusters
    ## If there was no ClustGroup ID use the next cluster number to label this group of clusters
    if(PrevClust == -Inf){ClustNum <- counter; counter <- counter + 1}else(ClustNum <-PrevClust)
    
    ## Label all shapes in the current cluster with the appropriate cluster number
    Areas_poly$ClustGroup[which(as.numeric(index) < 2*FieldLength)] <- ClustNum
    Areas_poly$ClustGroup[j] <- ClustNum
    
  }
  
}

## plot the different groupings
ggplot(data = Areas_poly) + 
  geom_sf(aes(fill = as.factor(ClustGroup))) +




## Now join the Cluster Groups together into single polygons
Areas_polys <- Areas_poly |> group_by(ClustGroup) |> summarise()
ggplot(data = Areas_polys) + geom_sf(aes(fill = as.factor(ClustGroup))) +
    labs(x = "Longitude", y = "Latitude", fill = "Cluster ID") +
    theme_light() +
    theme(legend.box.background=element_rect(colour = "#BDC3C7"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text=element_text(colour="black"),
          axis.title.x = element_text(size = 17),
          axis.text.x = element_text(hjust=0.7, size = 13,),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 17),
          axis.text.y = element_text(size = 13, hjust=0.7, angle=90, vjust=0.3))

# save the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/2-DefineWaderSites/Example Pics/Som_Final_clusters.png", units = "in", height = 9, width = 10)
  

## Save this shape file of the wader sites
write_sf(Areas_polys, "CleanData/Scenarios/2-DefineWaderSites/Somerset/Somerset_Wader_Sites.shp")



##----------------------------------------##
#### 1.7 Extract fields in each cluster ####
##----------------------------------------##

## Intersect the wader site boundaries and the BWWM fields shapefile
## Calculate the proportion of each field in the wader site
SomInter <- st_intersection(Areas_polys, FieldsClusts) |> 
            mutate(InterArea = st_area(geometry),
                   OverlapProp = as.numeric(InterArea)/as.numeric(OrigArea), 
                   Region = "Somerset") 

ggplot(data = SomInter) + geom_sf(aes(fill = as.factor(ClustGroup)))

## If a field is only 50% covered by a site them remove it from the site
SomInter <- filter(SomInter, OverlapProp >= 0.5) |> 
            select(F_LOC_ID, ClustGroup, Region) |> 
            st_drop_geometry()




##-----------------------------##
#### 2.0 Define Broads Sites ####
##-----------------------------##


##------------------------------------##
#### 2.1 Extract fields with waders ####
##------------------------------------##

## Filter out the field data for just Broads
## FieldX filter removes fields added in buffered priority landscape
Broads <- filter(Waders, Landscape == "Broads" & FieldX > 630000)
FieldLengthB <- sqrt(median(Broads$FieldArea))

## Work out which fields have at least one pair of the target species (Snipe, Redshank, Lapwing)
Broads <- Broads %>%  
  mutate(Tot_abund = rowSums(dplyr::select(., est_pairsL, est_pairsR), na.rm = T),
         Wader_Pres = ifelse(Tot_abund>0, "Y", ifelse(is.na(Tot_abund)== T, "N", "N"))) |> 
  select(-Tot_abund)

## Filter out the fields with at least one pair of target wader species present
Broads_Waders <- filter(Broads, Wader_Pres == "Y")

## plot the fields that have at least one pair of waders
ggplot() +
  geom_point(data = Broads_Waders, aes(x = FieldX, y= FieldY)) +
  theme_light()

## Save spatial data set to indicate which fields have waders
BroadsPres <- Broads |> select(F_LOC_ID, year, est_pairsL, est_pairsR, N_visits, Wader_Pres)
BroadsPres <- inner_join(BroadsPres, Fields, by = "F_LOC_ID")
write_sf(BroadsPres, "CleanData/Scenarios/2-DefineWaderSites/Broads/Broads_Wader_Present.shp")



##----------------------------##
#### 2.2 K-means clustering ####
##----------------------------##

## Get the X and Y coordinates for the fields with at least one pair of waders
Broads_WadersFields <- Broads_Waders |> select(FieldX, FieldY)

## Apply function to Broads field centers
wssplot(data = Broads_WadersFields)

## Run k-means clustering with chosen values for k
kmeans_outB <- kmeans(x = Broads_WadersFields, centers = 7, iter.max = 25, nstart = 25)
Broads_Waders$cluster <- kmeans_outB$cluster

## Plot the clusters using ggplot and the factoextra function fviz_cluster()
ggplot() +
  geom_point(data = Broads_Waders, aes(x = FieldX, y= FieldY, colour = as.factor(cluster))) +
  theme_light()

fviz_cluster(kmeans_outB, data = Broads_WadersFields)

## select the columns from the main data set that I want to keep
BroadsCLust <- Broads_Waders |> select(F_LOC_ID, est_pairsL, est_pairsR, cluster, FieldX, FieldY)




##------------------------------------##
#### 2.3 Calculate KDE's by cluster ####
##------------------------------------##

## Use `amt` functions to calculate KDE for each cluster
BroadsCLust <- as.data.frame(BroadsCLust)
Broadtrk <- make_track(tbl = BroadsCLust, .x = FieldX, .y = FieldY,
                     ID = cluster, crs = 27700)

## grouping factor to calculate KDE
groupfac <- "ID"

## make template raster for KDE's
trastB <- make_trast(Broadtrk, res = 25) # res for rastet that kde's calculated on

## make KDE per grouping factor
## use average field with as bandwidth for kernel density estimate
hr_Broad <- Broadtrk %>%  
  nest(data = -groupfac) %>%  
  mutate(kde = map(data, hr_kde, trast = trastB, levels = 0.95, h = c(FieldLengthB, FieldLengthB)))



##--------------------##
#### 2.4 Plot KDE's ####
##--------------------##

## Extract the polygon for the 95% KDE's from each cluster
for(j in 1:nrow(hr_Broad)){
  
  if(j == 1){AreasBr <- hr_isopleths(hr_Broad$kde[[j]])}else{AreasBr <- rbind(AreasBr,hr_isopleths(hr_Broad$kde[[j]]))}
  
}

## make plot of all 95% KDE's
ggplot() + 
    
    ## Add home ranges
    geom_sf(data = AreasBr, aes(fill = "blue"), alpha = 0.5) +
    scale_fill_manual(name = expression(""), values =c("blue", "blue"), labels = c("95% kde")) +
    
    ## change labels
    labs(x = "Longitude", y = "Latitude") +
    
    # set theme
    theme_light() +
    theme(legend.box.background=element_rect(colour = "#BDC3C7"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text=element_text(colour="black"),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(hjust=0.7),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
          axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

  

##-------------------------------##
#### 2.5 Remove small clusters ####
##-------------------------------##

## Separate mutliple polygons into distinct polygons (some 95% KDE's were not single polygons)
for(p in 1:nrow(AreasBr)){
  
  Inter <- as.data.frame(st_cast(st_sfc(AreasBr$geometry[p]), "POLYGON"))
  Inter$ClustID <- hr_Som$ID[p]
  
  if(p == 1){AreasBr_poly <- Inter}else{AreasBr_poly <- rbind(AreasBr_poly, Inter)}

}


## Assign each polygon a distinct number
AreasBr_poly <- st_as_sf(AreasBr_poly)
AreasBr_poly$ShapeID <- 1:nrow(AreasBr_poly)

## Calculate the centroid for each field with waders
BroadsPres_Cords <- filter(BroadsPres, Wader_Pres == "Y")
Broads_Cent <- st_centroid(st_as_sf(BroadsPres_Cords))

## find points within polygons
Field_in_Clust <- st_join(Broads_Cent, st_as_sf(AreasBr_poly), join = st_within)

## Calculate the number of fields with waders and wader pairs within each shape
## Remove shapes that do not have at least three occupied fields or three pairs of any species
AreasBr_poly <- Field_in_Clust |> 
  group_by(ShapeID) |> 
  summarise(n_fields = n(),
            n_lap = sum(est_pairsL, na.rm = T),
            n_red = sum(est_pairsR, na.rm = T),
            n_wad = sum(n_lap, n_red)) |> 
  st_drop_geometry() |> 
  filter(n_fields >= 3 | n_wad >= 3) |> 
  select(ShapeID, n_wad) |> 
  inner_join(AreasBr_poly, by = "ShapeID") |> 
  st_as_sf()

## plot the different groupings
ggplot(data = AreasBr_poly) + geom_sf(aes(fill = n_wad)) + theme_light()




##-----------------------------------------##
#### 2.6 Join close/overlapping clusters ####
##-----------------------------------------##

## Loop through each cluster/shape and group those that are closer than Field Length value
counter <- 1 # counter to keep track of the most recent cluster group ID number
AreasBr_poly$ClustGroup <- NA # initiate column where cluster group ID's go

for(j in 1:nrow(AreasBr_poly)){
  
  message("iteration ", j)
  
  ## Calculate distances to all other features
  index <- st_distance(AreasBr_poly[j,], AreasBr_poly[-j,])
  
  mindist <- as.numeric(min(index))
  
  ## if no cluster within the threshold distance then assign current cluster next cluster number
  if(mindist > 2*FieldLengthB){AreasBr_poly$ClustGroup[j] <- counter; counter <- counter + 1}
  
  ## if there are cluster within the threshold distance then they must be numbered together
  if(mindist < 2*FieldLengthB){
    
    ## Check if any cluster identified in this loop already have a cluster ID
    PrevClust <- max(AreasBr_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthB)], na.rm = T)
    
    ## If there was already a ClustGroup ID for any of the clusters below the threshold distance then use that number to label the current clusters
    ## If there was no ClustGroup ID use the next cluster number to label this group of clusters
    if(PrevClust == -Inf){ClustNum <- counter; counter <- counter + 1}else(ClustNum <-PrevClust)
    
    ## Label all shapes in the current cluster with the appropriate cluster number
    AreasBr_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthB)] <- ClustNum
    AreasBr_poly$ClustGroup[j] <- ClustNum
    
  }
  
}

## plot the different groupings
ggplot(data = AreasBr_poly) + geom_sf(aes(fill = as.factor(ClustGroup)))

## Now join the Cluster Groups together into single polygons
AreasBr_polys <- AreasBr_poly |> group_by(ClustGroup) |> summarise()
ggplot(data = AreasBr_polys) + geom_sf(aes(fill = as.factor(ClustGroup))) + theme_light()

## Save this shape file of the wader sites
write_sf(AreasBr_polys, "CleanData/Scenarios/2-DefineWaderSites/Broads/Broads_Wader_Sites.shp")



##----------------------------------------##
#### 2.7 Extract fields in each cluster ####
##----------------------------------------##

## Intersect the wader site boundaries and the BWWM fields shapefile
## Calculate the proportion of each field in the wader site
BroadInter <- st_intersection(AreasBr_polys, FieldsClusts) |> 
            mutate(InterArea = st_area(geometry),
                   OverlapProp = as.numeric(InterArea)/as.numeric(OrigArea), 
                   Region = "Broads") 

ggplot(data = BroadInter) + geom_sf(aes(fill = as.factor(ClustGroup)))

## If a field is only 50% covered by a site them remove it from the site
BroadInter <- filter(BroadInter, OverlapProp >= 0.5) |> 
            select(F_LOC_ID, ClustGroup, Region) |> 
            st_drop_geometry()




##----------------------------##
#### 3.0 Define Essex Sites ####
##----------------------------##

##------------------------------------##
#### 3.1 Extract fields with waders ####
##------------------------------------##

## Filter out the field data for just Broads
## FieldX filter removes fields added in buffered priority landscape
Thames <- filter(Waders, Landscape == "Greater Thames")
Essex <- filter(Thames, FieldY > 183000)
FieldLengthE <-sqrt(median(Essex$FieldArea))


## Work out which fields have at least one pair of the target species (Snipe, Redshank, Lapwing)
Essex <- Essex %>%  
  mutate(Tot_abund = rowSums(dplyr::select(., est_pairsL, est_pairsR), na.rm = T),
         Wader_Pres = ifelse(Tot_abund>0, "Y", ifelse(is.na(Tot_abund)== T, "N", "N"))) |> 
  select(-Tot_abund)

## Filter out the fields with at least one pair of target wader species present
Essex_Waders <- filter(Essex, Wader_Pres == "Y")

## plot the fields that have at least one pair of waders
ggplot() +
  geom_point(data = Essex_Waders, aes(x = FieldX, y= FieldY)) +
  theme_light()

## Save spatial data set to indicate which fields have waders
EssexPres <- Essex |> select(F_LOC_ID, year, est_pairsL, est_pairsR, N_visits, Wader_Pres)
EssexPres <- inner_join(EssexPres, Fields, by = "F_LOC_ID")
write_sf(EssexPres, "CleanData/Scenarios/2-DefineWaderSites/Essex/Essex_Wader_Present.shp")



##----------------------------##
#### 3.2 K-means clustering ####
##----------------------------##

## Get the X and Y coordinates for the fields with at least one pair of waders
Essex_WadersFields <- Essex_Waders |> select(FieldX, FieldY)

## Apply function to Essex field centers
wssplot(data = Essex_WadersFields)

## Run k-means clustering with chosen values for k
kmeans_outE <- kmeans(x = Essex_WadersFields, centers = 7, iter.max = 25, nstart = 25)
Essex_Waders$cluster <- kmeans_outE$cluster

## Plot the clusters using ggplot and the factoextra function fviz_cluster()
ggplot() +
  geom_point(data = Essex_Waders, aes(x = FieldX, y= FieldY, colour = as.factor(cluster))) +
  theme_light()

fviz_cluster(kmeans_outE, data = Essex_WadersFields)

## select the columns from the main data set that I want to keep
EssexCLust <- Essex_Waders |> select(F_LOC_ID, est_pairsL, est_pairsR, cluster, FieldX, FieldY)




##------------------------------------##
#### 3.3 Calculate KDE's by cluster ####
##------------------------------------##

## Use `amt` functions to calculate KDE for each cluster
EssexCLust <- as.data.frame(EssexCLust)
Essextrk <- make_track(tbl = EssexCLust, .x = FieldX, .y = FieldY,
                     ID = cluster, crs = 27700)

## grouping factor to calculate KDE
groupfac <- "ID"

## make template raster for KDE's
trastE <- make_trast(Essextrk, res = 25) # res for rastet that kde's calculated on

## make KDE per grouping factor
## use average field with as bandwidth for kernel density estimate
hr_Essex <- Essextrk %>%  
  nest(data = -groupfac) %>%  
  mutate(kde = map(data, hr_kde, trast = trastE, levels = 0.95, h = c(FieldLengthE, FieldLengthE)))



##--------------------##
#### 3.4 Plot KDE's ####
##--------------------##

## Extract the polygon for the 95% KDE's from each cluster
for(j in 1:nrow(hr_Essex)){
  
  if(j == 1){AreasEs <- hr_isopleths(hr_Essex$kde[[j]])}else{AreasEs <- rbind(AreasEs,hr_isopleths(hr_Essex$kde[[j]]))}
  
}

## make plot of all 95% KDE's
ggplot() + 
    
    ## Add home ranges
    geom_sf(data = AreasEs, aes(fill = "blue"), alpha = 0.5) +
    scale_fill_manual(name = expression(""), values =c("blue", "blue"), labels = c("95% kde")) +
    
    ## change labels
    labs(x = "Longitude", y = "Latitude") +
    
    # set theme
    theme_light() +
    theme(legend.box.background=element_rect(colour = "#BDC3C7"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text=element_text(colour="black"),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(hjust=0.7),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
          axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

  

##-------------------------------##
#### 3.5 Remove small clusters ####
##-------------------------------##

## Separate mutliple polygons into distinct polygons (some 95% KDE's were not single polygons)
for(p in 1:nrow(AreasEs)){
  
  Inter <- as.data.frame(st_cast(st_sfc(AreasEs$geometry[p]), "POLYGON"))
  Inter$ClustID <- hr_Som$ID[p]
  
  if(p == 1){AreasEs_poly <- Inter}else{AreasEs_poly <- rbind(AreasEs_poly, Inter)}

}


## Assign each polygon a distinct number
AreasEs_poly <- st_as_sf(AreasEs_poly)
AreasEs_poly$ShapeID <- 1:nrow(AreasEs_poly)

## Calculate the centroid for each field with waders
EssexPres_Cords <- filter(EssexPres, Wader_Pres == "Y")
Essex_Cent <- st_centroid(st_as_sf(EssexPres_Cords))

## find points within polygons
Field_in_Clust <- st_join(Essex_Cent, st_as_sf(AreasEs_poly), join = st_within)

## Calculate the number of fields with waders and wader pairs within each shape
## Remove shapes that do not have at least three occupied fields or three pairs of any species
AreasEs_poly <- Field_in_Clust |> 
  group_by(ShapeID) |> 
  summarise(n_fields = n(),
            n_lap = sum(est_pairsL, na.rm = T),
            n_red = sum(est_pairsR, na.rm = T),
            n_wad = sum(n_lap, n_red)) |> 
  st_drop_geometry() |> 
  filter(n_fields >= 3 | n_wad >= 3) |> 
  select(ShapeID, n_wad) |> 
  inner_join(AreasEs_poly, by = "ShapeID") |> 
  st_as_sf()

## plot the different groupings
ggplot(data = AreasEs_poly) + geom_sf(aes(fill = n_wad)) + theme_light()




##-----------------------------------------##
#### 3.6 Join close/overlapping clusters ####
##-----------------------------------------##

## Loop through each cluster/shape and group those that are closer than Field Length value
counter <- 1 # counter to keep track of the most recent cluster group ID number
AreasEs_poly$ClustGroup <- NA # initiate column where cluster group ID's go

for(j in 1:nrow(AreasEs_poly)){
  
  message("iteration ", j)
  
  ## Calculate distances to all other features
  index <- st_distance(AreasEs_poly[j,], AreasEs_poly[-j,])
  
  mindist <- as.numeric(min(index))
  
  ## if no cluster within the threshold distance then assign current cluster next cluster number
  if(mindist > 2*FieldLengthE){AreasEs_poly$ClustGroup[j] <- counter; counter <- counter + 1}
  
  ## if there are cluster within the threshold distance then they must be numbered together
  if(mindist < 2*FieldLengthE){
    
    ## Check if any cluster identified in this loop already have a cluster ID
    PrevClust <- max(AreasEs_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthE)], na.rm = T)
    
    ## If there was already a ClustGroup ID for any of the clusters below the threshold distance then use that number to label the current clusters
    ## If there was no ClustGroup ID use the next cluster number to label this group of clusters
    if(PrevClust == -Inf){ClustNum <- counter; counter <- counter + 1}else(ClustNum <-PrevClust)
    
    ## Label all shapes in the current cluster with the appropriate cluster number
    AreasEs_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthE)] <- ClustNum
    AreasEs_poly$ClustGroup[j] <- ClustNum
    
  }
  
}

## plot the different groupings
ggplot(data = AreasEs_poly) + geom_sf(aes(fill = as.factor(ClustGroup)))

## Now join the Cluster Groups together into single polygons
AreasEs_polys <- AreasEs_poly |> group_by(ClustGroup) |> summarise()
ggplot(data = AreasEs_polys) + geom_sf(aes(fill = as.factor(ClustGroup))) + theme_light()

## Save this shape file of the wader sites
write_sf(AreasEs_polys, "CleanData/Scenarios/2-DefineWaderSites/Essex/Essex_Wader_Sites.shp")



##----------------------------------------##
#### 3.7 Extract fields in each cluster ####
##----------------------------------------##

## Intersect the wader site boundaries and the BWWM fields shapefile
## Calculate the proportion of each field in the wader site
EssexInter <- st_intersection(AreasEs_polys, FieldsClusts) |> 
            mutate(InterArea = st_area(geometry),
                   OverlapProp = as.numeric(InterArea)/as.numeric(OrigArea), 
                   Region = "Essex") 

ggplot(data = EssexInter) + geom_sf(aes(fill = as.factor(ClustGroup)))

## If a field is only 50% covered by a site them remove it from the site
EssexInter <- filter(EssexInter, OverlapProp >= 0.5) |> 
            select(F_LOC_ID, ClustGroup, Region) |> 
            st_drop_geometry()




##---------------------------##
#### 4.0 Define Kent Sites ####
##---------------------------##

##------------------------------------##
#### 4.1 Extract fields with waders ####
##------------------------------------##

## Filter out the field data for just Broads
## FieldX filter removes fields added in buffered priority landscape
Thames <- filter(Waders, Landscape == "Greater Thames")
Kent <- filter(Thames, FieldY < 183000)
FieldLengthK <-sqrt(median(Kent$FieldArea))


## Work out which fields have at least one pair of the target species (Snipe, Redshank, Lapwing)
Kent <- Kent %>%  
  mutate(Tot_abund = rowSums(dplyr::select(., est_pairsL, est_pairsR), na.rm = T),
         Wader_Pres = ifelse(Tot_abund>0, "Y", ifelse(is.na(Tot_abund)== T, "N", "N"))) |> 
  select(-Tot_abund)

## Filter out the fields with at least one pair of target wader species present
Kent_Waders <- filter(Kent, Wader_Pres == "Y")

## plot the fields that have at least one pair of waders
ggplot() +
  geom_point(data = Kent_Waders, aes(x = FieldX, y= FieldY)) +
  theme_light()

## Save spatial data set to indicate which fields have waders
KentPres <- Kent |> select(F_LOC_ID, year, est_pairsL, est_pairsR, N_visits, Wader_Pres)
KentPres <- inner_join(KentPres, Fields, by = "F_LOC_ID")
write_sf(KentPres, "CleanData/Scenarios/2-DefineWaderSites/Kent/Kent_Wader_Present.shp")



##----------------------------##
#### 4.2 K-means clustering ####
##----------------------------##

## Get the X and Y coordinates for the fields with at least one pair of waders
Kent_WadersFields <- Kent_Waders |> select(FieldX, FieldY)

## Apply function to Kent field centers
wssplot(data = Kent_WadersFields)

## Run k-means clustering with chosen values for k
kmeans_outK <- kmeans(x = Kent_WadersFields, centers = 5, iter.max = 25, nstart = 25)
Kent_Waders$cluster <- kmeans_outK$cluster

## Plot the clusters using ggplot and the factoextra function fviz_cluster()
ggplot() +
  geom_point(data = Kent_Waders, aes(x = FieldX, y= FieldY, colour = as.factor(cluster))) +
  theme_light()

fviz_cluster(kmeans_outK, data = Kent_WadersFields)

## select the columns from the main data set that I want to keep
KentCLust <- Kent_Waders |> select(F_LOC_ID, est_pairsL, est_pairsR, cluster, FieldX, FieldY)




##------------------------------------##
#### 4.3 Calculate KDE's by cluster ####
##------------------------------------##

## Use `amt` functions to calculate KDE for each cluster
KentCLust <- as.data.frame(KentCLust)
Kenttrk <- make_track(tbl = KentCLust, .x = FieldX, .y = FieldY,
                     ID = cluster, crs = 27700)

## grouping factor to calculate KDE
groupfac <- "ID"

## make template raster for KDE's
trastK <- make_trast(Kenttrk, res = 25) # res for rastet that kde's calculated on

## make KDE per grouping factor
## use average field with as bandwidth for kernel density estimate
hr_Kent <- Kenttrk %>%  
  nest(data = -groupfac) %>%  
  mutate(kde = map(data, hr_kde, trast = trastK, levels = 0.95, h = c(FieldLengthK, FieldLengthK)))



##--------------------##
#### 4.4 Plot KDE's ####
##--------------------##

## Extract the polygon for the 95% KDE's from each cluster
for(j in 1:nrow(hr_Kent)){
  
  if(j == 1){AreasKe <- hr_isopleths(hr_Kent$kde[[j]])}else{AreasKe <- rbind(AreasKe,hr_isopleths(hr_Kent$kde[[j]]))}
  
}

## make plot of all 95% KDE's
ggplot() + 
    
    ## Add home ranges
    geom_sf(data = AreasKe, aes(fill = "blue"), alpha = 0.5) +
    scale_fill_manual(name = expression(""), values =c("blue", "blue"), labels = c("95% kde")) +
    
    ## change labels
    labs(x = "Longitude", y = "Latitude") +
    
    # set theme
    theme_light() +
    theme(legend.box.background=element_rect(colour = "#BDC3C7"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text=element_text(colour="black"),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(hjust=0.7),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
          axis.text.y = element_text(hjust=0.7, angle=90, vjust=0.3))

  

##-------------------------------##
#### 4.5 Remove small clusters ####
##-------------------------------##

## Separate mutliple polygons into distinct polygons (some 95% KDE's were not single polygons)
for(p in 1:nrow(AreasKe)){
  
  Inter <- as.data.frame(st_cast(st_sfc(AreasKe$geometry[p]), "POLYGON"))
  Inter$ClustID <- hr_Som$ID[p]
  
  if(p == 1){AreasKe_poly <- Inter}else{AreasKe_poly <- rbind(AreasKe_poly, Inter)}

}


## Assign each polygon a distinct number
AreasKe_poly <- st_as_sf(AreasKe_poly)
AreasKe_poly$ShapeID <- 1:nrow(AreasKe_poly)

## Calculate the centroid for each field with waders
KentPres_Cords <- filter(KentPres, Wader_Pres == "Y")
Kent_Cent <- st_centroid(st_as_sf(KentPres_Cords))

## find points within polygons
Field_in_Clust <- st_join(Kent_Cent, st_as_sf(AreasKe_poly), join = st_within)

## Calculate the number of fields with waders and wader pairs within each shape
## Remove shapes that do not have at least three occupied fields or three pairs of any species
AreasKe_poly <- Field_in_Clust |> 
  group_by(ShapeID) |> 
  summarise(n_fields = n(),
            n_lap = sum(est_pairsL, na.rm = T),
            n_red = sum(est_pairsR, na.rm = T),
            n_wad = sum(n_lap, n_red)) |> 
  st_drop_geometry() |> 
  filter(n_fields >= 3 | n_wad >= 3) |> 
  select(ShapeID, n_wad) |> 
  inner_join(AreasKe_poly, by = "ShapeID") |> 
  st_as_sf()

## plot the different groupings
ggplot(data = AreasKe_poly) + geom_sf(aes(fill = n_wad)) + theme_light()



##-----------------------------------------##
#### 4.6 Join close/overlapping clusters ####
##-----------------------------------------##

## Loop through each cluster/shape and group those that are closer than Field Length value
counter <- 1 # counter to keep track of the most recent cluster group ID number
AreasKe_poly$ClustGroup <- NA # initiate column where cluster group ID's go

for(j in 1:nrow(AreasKe_poly)){
  
  message("iteration ", j)
  
  ## Calculate distances to all other features
  index <- st_distance(AreasKe_poly[j,], AreasKe_poly[-j,])
  
  mindist <- as.numeric(min(index))
  
  ## if no cluster within the threshold distance then assign current cluster next cluster number
  if(mindist > 2*FieldLengthK){AreasKe_poly$ClustGroup[j] <- counter; counter <- counter + 1}
  
  ## if there are cluster within the threshold distance then they must be numbered together
  if(mindist < 2*FieldLengthK){
    
    ## Check if any cluster identified in this loop already have a cluster ID
    PrevClust <- max(AreasKe_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthK)], na.rm = T)
    
    ## If there was already a ClustGroup ID for any of the clusters below the threshold distance then use that number to label the current clusters
    ## If there was no ClustGroup ID use the next cluster number to label this group of clusters
    if(PrevClust == -Inf){ClustNum <- counter; counter <- counter + 1}else(ClustNum <-PrevClust)
    
    ## Label all shapes in the current cluster with the appropriate cluster number
    AreasKe_poly$ClustGroup[which(as.numeric(index) < 2*FieldLengthK)] <- ClustNum
    AreasKe_poly$ClustGroup[j] <- ClustNum
    
  }
  
}

## plot the different groupings
ggplot(data = AreasKe_poly) + geom_sf(aes(fill = as.factor(ClustGroup)))

## Now join the Cluster Groups together into single polygons
AreasKe_polys <- AreasKe_poly |> group_by(ClustGroup) |> summarise()
ggplot(data = AreasKe_polys) + geom_sf(aes(fill = as.factor(ClustGroup))) + theme_light()

## Save this shape file of the wader sites
write_sf(AreasKe_polys, "CleanData/Scenarios/2-DefineWaderSites/Kent/Kent_Wader_Sites.shp")



##----------------------------------------##
#### 4.7 Extract fields in each cluster ####
##----------------------------------------##

## Intersect the wader site boundaries and the BWWM fields shapefile
## Calculate the proportion of each field in the wader site
KentInter <- st_intersection(AreasKe_polys, FieldsClusts) |> 
            mutate(InterArea = st_area(geometry),
                   OverlapProp = as.numeric(InterArea)/as.numeric(OrigArea), 
                   Region = "Kent") 

ggplot(data = KentInter) + geom_sf(aes(fill = as.factor(ClustGroup)))

## If a field is only 50% covered by a site them remove it from the site
KentInter <- filter(KentInter, OverlapProp >= 0.5) |> 
            select(F_LOC_ID, ClustGroup, Region) |> 
            st_drop_geometry()






##-----------------------------------##
#### 5.0 Join all regional cluster ####
##-----------------------------------##

## Bind all the regional data sets of fields within a wader site
All_Inter <- rbind(KentInter, EssexInter, SomInter, BroadInter)

## Join the cluster IDs to the shapefile of all BWWM fields
FieldsClustsAll <- left_join(FieldsClusts, All_Inter, by = "F_LOC_ID")

## save the data set as a shapefile
write_sf(FieldsClustsAll, "CleanData/Scenarios/2-DefineWaderSites/All Regions/All_BWWM_Fields_Clusters.shp")












