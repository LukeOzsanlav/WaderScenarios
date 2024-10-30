##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 17/10/2023
## Goal: Refine data set for statistical analysis
## 
## Breakdown:
##  1. Select priority landscape fields 
##  2. Add spatial properties of each field 
##  3. Remove non-grassland fields
##  4. Plots of breeding pairs
##  
##------------------------------------------------------##


## Load in required packages
pacman::p_load(here, tidyverse, data.table, sf, terra, gt)
source("Code/Helper functions.R") # get helper functions

## Use `here` package for specifying file paths
here::i_am("Code/Wader Abundance/2 - Field Selection + Plotting.R")
library(here)
here()

## Queries: 
## When removing non-grassland fields there are some fields that are removed as they are salt marsh but they are RSPB reserves (e.g. Wallasea Island)
## Some fields that are unsuitable are actually grassland, not obvious what to do with these




##---------------------##
#### 0. Data Read In ####
##---------------------##

##-----------------------------##
#### 0.1 Read in survey data ####
##-----------------------------##

## read in cleaned survey data output in "Cleaning Script.R"
Wader <- read_csv(here("CleanData", "Wader Abundance", "1-CleaningScript", "Estimated_Pairs_clean.csv"))




####----------------------------##
#### 0.2 Read in spatial data ####
##------------------------------##

## Read in high res world UK coastline (for plotting to give context)
UK <- st_read(here("RawData", "UK_Coastline", "UK_Coatline.shp"))

## Read in the field shapes GIS layer
FieldShapes <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp"))

## Read in RSPB priority landscapes
Lanscap <- st_read(here("RawData", "Priority Landscapes", "EnglandWales_PriorityLandscapes.shp"))
Lanscap <- filter(Lanscap, Att3Value %in% c("Somerset Levels and Moors", "Greater Thames", 
                                            "Broads", "Suffolk Coast")) # Filter out the landscapes I am focusing on

## Read in RSPB priority landscapes, with greater Thames landscape split in two
SplitLanscape <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |> 
                        filter(Att3Value %in% c("Somerset Levels and Moors", "Broads", "Greater Thames", "North Kent", "Essex"))

## Read in the UKCEH landcover data
LC <- rast(here("RawData", "LandCover", "gblcm25m2021.tif"))
LC <- LC[[1]]




##-----------------------------------##
#### 1. Spatial Data Manipulations ####
##-----------------------------------##

##---------------------------##
#### 1.1 Plot spatial data ####
##---------------------------##

## plot field shapes
p1 <- ggplot() + geom_sf(data = FieldShapes, mapping = aes(geometry = geometry), fill = "orange") +
         geom_sf(data = UK, mapping = aes(geometry = geometry), fill = NA) +
         theme_light()

## plot priority landscapes
p2 <- ggplot() + geom_sf(data = Lanscap, mapping = aes(geometry = geometry, fill = Att3Value), colour = NA) +
         geom_sf(data = UK, mapping = aes(geometry = geometry), fill = NA) +
         coord_sf(xlim = c(327378.1-10000, 653935.4+10000), ylim = c(120152.9-10000, 330111.0+10000), crs = 27700,
                  expand = FALSE) +
         theme_light()




##----------------------------------##
#### 1.2 Create Landscape Buffers ####
##----------------------------------##

##** Combine and buffer Priority Landscapes and ESAs
## Create bounds for Somerset area
Somerset_Prior <- filter(Lanscap, Att3Value == "Somerset Levels and Moors")
Somerset_buf <-st_buffer(Somerset_Prior, dist = 5000)
plot(Somerset_buf$geometry)

## Create bounds for Norfolk Broads
Broads_Prior <- filter(Lanscap, Att3Value == "Broads")
Broads_buf <-st_buffer(Broads_Prior, dist = 5000)
plot(Broads_buf$geometry)

## Create bounds for Greater Thames
Thames_Prior <- filter(Lanscap, Att3Value == "Greater Thames")
Thames_buf <-st_buffer(Thames_Prior, dist = 5000)
plot(Thames_buf$geometry)

## Now bind all of these shapes together in one data set
Lanscap_buf <- rbind(Thames_buf, Broads_buf, Somerset_buf)
gc()




##------------------------------------##
#### 1.3 Calculate field dimensions ####
##------------------------------------##

## Subset field shape file to a few important columns
Fields <- FieldShapes |> select(Field_no, AES, F_LOC_ID, S_LOC_ID, REG, X2023_field)

## Calculate the field area and the field centroids
Fields <- Fields |> 
  mutate(FieldArea = as.numeric(st_area(geometry)),
         FieldX = st_coordinates(st_centroid(geometry))[,1],
         FieldY = st_coordinates(st_centroid(geometry))[,2])




##----------------------------------------------------##
#### 1.4 Identify fields within landscape boundaries ####
##----------------------------------------------------##

## Work out if fields fall within any of the buffered priority landscapes
Overlap <- st_intersection(Lanscap_buf, Fields)

## Join the priority landscape labels onto that data set with all the fields in
## first streamline the data set from the overlap
Overlap <- Overlap |> 
  select(F_LOC_ID, Att3Value) |> 
  rename(Landscape = Att3Value) |> 
  st_drop_geometry()

## Join the data sets together, and for row change in join
PreL <- nrow(Fields)
Fields <- full_join(Fields, Overlap, by = "F_LOC_ID")
nrow(Fields) == PreL # check no rows lost


## VISUALISE THE INTERSECTION ##

## Filter the fields identified in the Somerset region and the Somerset landscape outline
SomersetFields <- filter(Fields, Landscape == "Somerset Levels and Moors")

## plot the priority region and fields inside in orange and fields not inside in grey
Int1 <- ggplot() + 
           geom_sf(data = Somerset_buf, mapping = aes(geometry = geometry), fill = "lightgrey", colour = NA) +
           geom_sf(data = FieldShapes, mapping = aes(geometry = geometry), fill = "darkgrey") +
           geom_sf(data = SomersetFields, mapping = aes(geometry = geometry), fill = "orange", colour = NA) +
           geom_sf(data = Somerset_Prior, mapping = aes(geometry = geometry), fill = NA, colour = "blue") +
           coord_sf(xlim = c(st_bbox(Somerset_buf)[1]-1000, st_bbox(Somerset_buf)[3]+1000), 
                    ylim = c(st_bbox(Somerset_buf)[2]-1000, st_bbox(Somerset_buf)[4]+1000), crs = 27700,
                    expand = TRUE) +
           theme_light()


## Filter the fields identified in the Broads region and the Broads landscape outline
BroadsFields <- filter(Fields, Landscape == "Broads")

## plot the priority region and fields inside in orange and fields not inside in grey
Int2 <- ggplot() + 
           geom_sf(data = Broads_buf, mapping = aes(geometry = geometry), fill = "lightgrey", colour = NA) +
           geom_sf(data = FieldShapes, mapping = aes(geometry = geometry), fill = "darkgrey") +
           geom_sf(data = BroadsFields, mapping = aes(geometry = geometry), fill = "orange", colour = NA) +
           geom_sf(data = Broads_Prior, mapping = aes(geometry = geometry), fill = NA, colour = "blue") +
           coord_sf(xlim = c(st_bbox(Broads_buf)[1]-1000, st_bbox(Broads_buf)[3]+1000), 
                    ylim = c(st_bbox(Broads_buf)[2]-1000, st_bbox(Broads_buf)[4]+1000), crs = 27700,
                    expand = TRUE) +
           theme_light()


## Filter the fields identified in the Thames region and the Thames landscape outline
ThamesFields <- filter(Fields, Landscape == "Greater Thames")

## plot the priority region and fields inside in orange and fields not inside in grey
Int3 <- ggplot() +
           geom_sf(data = Thames_buf, mapping = aes(geometry = geometry), fill = "lightgrey", colour = NA) +
           geom_sf(data = FieldShapes, mapping = aes(geometry = geometry), fill = "darkgrey") +
           geom_sf(data = ThamesFields, mapping = aes(geometry = geometry), fill = "orange", colour = NA) +
           geom_sf(data = Thames_Prior, mapping = aes(geometry = geometry), fill = NA, colour = "blue") +
           coord_sf(xlim = c(st_bbox(Thames_buf)[1]-1000, st_bbox(Thames_buf)[3]+1000), 
                    ylim = c(st_bbox(Thames_buf)[2]-1000, st_bbox(Thames_buf)[4]+1000), crs = 27700,
                    expand = TRUE) +
           theme_light()




##-----------------------------------------------------##
#### 2.1 Join Pair Estimates and Spatial Calculations ####
##-----------------------------------------------------##

## Trim Fields data set for join
FieldSub <- Fields |> select(F_LOC_ID, FieldArea, FieldX, FieldY, Landscape)

## Join Spatial calculation and wader survey data
PreL <- nrow(Wader)
Wader <- left_join(Wader, FieldSub, by = "F_LOC_ID")
nrow(Wader) == PreL



##-------------------------##
#### 2.2 Summary of data ####
##-------------------------##

## filter out the fields in priority landscapes
PL_data <- filter(Wader, is.na(Landscape)==F)
nrow(PL_data) # the number of fields across all landscapes
sum(PL_data$FieldArea)/10000 # total survey area in hectares

## summarise the data by landscape
PL_data |> 
  group_by(Landscape) |> 
  summarise(TotArea = sum(FieldArea)/10000,
            TotFields = n())

## calculate the number of fields with only one visit
Visits <- table(PL_data$N_visits)
Visits[1]/sum(Visits)




##------------------------------##
#### 3. Field Selection Step  ####
##------------------------------##

##** CRITERIA 1: Priority Landscape filter ##
## Remove any fields that are not in the buffered priority landscape
PFields <- filter(Wader, is.na(Landscape)==F)
StartLength <- nrow(PFields) # record how many fields before filtering


## Work out dominant habitat in each field and the proportion of each field that is grassland
## Firstly, crop the data set to the south of England
LC <- crop(LC, vect(Lanscap_buf))
plot(LC)

## Function to calculate the mode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]} 

## extract the dominant habitat from the UKCEH raster (using get mode function)
MainHab <- extract(LC, vect(PFields$geometry), fun = getmode)

## Add these dominant habitat categories to the main data set
PFields$DomHab <- MainHab$gblcm25m2021_1
table(PFields$DomHab)

## Calculate proportion of each that is grassland
## Firstly, reclassify raster so grassland is 1 and non-grassland is 0
## NOTE: Fen/Swamp on satellite looks like wet grassland
m <- c(0, 3.5, 0,
       3.6, 8.5, 1,
       8.6, 21.5, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
rc1 <- classify(LC, rclmat, include.lowest=TRUE) # reclassify the raster
plot(rc1)

## Calculate the mean weighted pixel value for each field
PropGrass <- extract(rc1, vect(PFields$geometry), fun = mean, weights = T)
hist(PropGrass$gblcm25m2021_1)

## Assign these weighted means to the main data set
PFields$PropGrass <- PropGrass$gblcm25m2021_1



##** CRITERIA 2: Remove Arable ##
## remove any fields that were assigned as Arable during every survey visit
PFields <- filter(PFields, !Arable_Hab == 1)

##** CRITERIA 3: Remove unsuitable non-grassland ##
## remove any fields that were assigned as none-grassland from the UKCEH land cover maps
## AND that were identified as unsuitable or suitable other in every survey visit
SelectFields <- filter(PFields, !(DomHab %in% c(1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 20, 21) & (Unsuit_Hab == 1 | Other_Hab == 1)))


## Check size of current data set
nrow(PFields) # number of fields left
StartLength-nrow(PFields) # number of fields removed
((StartLength-nrow(PFields))/StartLength)*100 # % of fields removed
PFields |> 
  group_by(Landscape) |> 
  summarise(TotArea = sum(FieldArea)/10000,
            TotFields = n())



##--------------------------------------------------------##
#### 4. Save Breeding Pair data for Priority Landscapes ####
##--------------------------------------------------------##

## First, Add column whether waders are present
SurveySum <- SelectFields %>%  
  mutate(WaderTot = rowSums(dplyr::select(., est_pairsL, est_pairsO, est_pairsC, est_pairsR, est_pairsS), na.rm = T),
         wader_occ = ifelse(WaderTot>0, "Y", ifelse(is.na(WaderTot)== T, "N", "N"))) |> 
  select(-WaderTot)

## Write out this file here for plotting in QGIS
SurveySum |> select(-geometry) |> 
write_csv(here("CleanData", "Wader Abundance", "2-FieldSelectionPlotting", "Breeding_Pairs_Selected.csv"))
write_sf(st_as_sf(SurveySum), here("CleanData", "Wader Abundance", "2-FieldSelectionPlotting", "Breeding_Pairs_Selected.shp"))




##---------------------------------------------------##
#### 5. Plots of Wader Abundance and Survey Effort ####
##---------------------------------------------------##

##-------------------------------##
#### 5.1 Prepare data for plot ####
##-------------------------------##

##** Make heat maps of:
## year, whether surveyed, number of visits, wader distribution
## Need to do a bit of data prep first

## Join to field shapes to work out fields not surveyed
SurveyCov <- SelectFields |> 
  filter(Somerset %in% c("All", "S")) |> 
  select(F_LOC_ID) |>  
  mutate(Covered = "Y") |> 
  full_join(Fields, by = "F_LOC_ID") |> 
  mutate(Covered = ifelse(is.na(Covered)==T, "N", Covered))


            
##------------------------##
#### 5.2 Somerset plots ####
##------------------------##

## Filter out the data for Somerset
## Remember there will be different datasets for snipe and Lapwing
SomersetSurv_Lap <- filter(SurveySum, Landscape == "Somerset Levels and Moors" & Somerset %in% c("All", "L"))
SomersetSurv_Sni <- filter(SurveySum, Landscape == "Somerset Levels and Moors" & Somerset %in% c("All", "S"))
SomersetCov <- filter(SurveyCov, Landscape == "Somerset Levels and Moors")
Somerset <- filter(Lanscap, Att3Value == "Somerset Levels and Moors")

## Create a base plot that I can then add field data on top of
SL_Base <- ggplot() + geom_sf(data = Somerset, mapping = aes(geometry = geometry), fill = "#f5f5f2", colour = "NA") +
           geom_sf(data = Somerset, mapping = aes(geometry = geometry), fill = NA, colour = "orange") +
           coord_sf(xlim = c(st_bbox(SomersetSurv_Sni$geometry)[1]-1000, st_bbox(SomersetSurv_Sni$geometry)[3]+1000), 
                    ylim = c(st_bbox(SomersetSurv_Sni$geometry)[2]-1000, st_bbox(SomersetSurv_Sni$geometry)[4]+1000), crs = 27700,
                    expand = F) +
           theme_minimal() +
           theme(text = element_text(family = "Karla", color = "#2D2D2E"), # change the font
            legend.position = c(0.9,0.1), # Move the legend to the bottom left corner
            legend.margin = margin(t = 0, unit = "cm"), # Set top of legend margin to zero
            legend.title = element_text(size = 13), # Adjust legend text size
            legend.text = element_text(size = 10, color = "#878585"), # Adjust legend text size
            legend.key.width = unit(1, "cm"), # Change legend colorbar width
            legend.key.height = unit(0.5, "cm"), # Change legend colorbar height
            legend.background = element_blank(), # Remove legend background
            axis.title = element_blank(), # Remove axis titles
            panel.grid = element_line(color = "#ebebe5", linewidth = 0.2), # Change lat/long graticule line colors
            panel.background = element_rect(fill = "white", color = NA)) +  # Change plot panel background color
            annotate(geom = "text", x = st_bbox(SomersetSurv_Sni$geometry)[1]-1000, y = st_bbox(SomersetSurv_Sni$geometry)[4],
                     label = "Somerset", family = "Karla", color = "#2D2D2E",
                     size = 9, hjust = 0)


## Plot the number of visits to each field
## create colour pallet so values are easy to tell apart 
my_colors <- RColorBrewer::brewer.pal(6, "Blues")[c(2,3,5,6)]

## create the plot for number of visits per field
SL_Base + geom_sf(data = SomersetSurv_Sni, mapping = aes(geometry = geometry, fill = as.character(N_visits)), colour = NA) +
          scale_fill_manual(values = my_colors, labels = c("1", "2","3", "4"), name = "# Visits") +
          annotate(geom = "text", x = st_bbox(SomersetSurv_Sni$geometry)[1]-1000, y = st_bbox(SomersetSurv_Sni$geometry)[4]-1400,
                   label = "Visit number", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Somerset/Somerset_Visits.png", width = 18, height = 20, units = "cm")



## Plot estimated number of wader pairs
## Lapwing breeding pairs estimates
SomersetSurv_Lap <- SomersetSurv_Lap |> 
  mutate(est_pairsL = ifelse(est_pairsL == 0, NA, est_pairsL))
SL_Base + geom_sf(data = SomersetSurv_Lap, mapping = aes(geometry = geometry, fill = est_pairsL), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Lapwing Pairs")  +
          annotate(geom = "text", x = st_bbox(SomersetSurv_Lap$geometry)[1]-1000, y = st_bbox(SomersetSurv_Lap$geometry)[4]-1400,
                   label = "Breeding Lapwing", family = "Karla", color = "#566573", size = 5, hjust = 0)

ggsave("Plots/Somerset/Somerset_Lapwing.png", width = 18, height = 20, units = "cm")

## Snipe breeding pairs estimates
SomersetSurv_Sni <- SomersetSurv_Sni |> 
  mutate(est_pairsS = ifelse(est_pairsS == 0, NA, est_pairsS))
SL_Base + geom_sf(data = SomersetSurv_Sni, mapping = aes(geometry = geometry, fill = est_pairsS), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Snipe Pairs") +
          annotate(geom = "text", x = st_bbox(SomersetSurv_Sni$geometry)[1]-1000, y = st_bbox(SomersetSurv_Sni$geometry)[4]-1400,
                   label = "Breeding Snipe", family = "Karla", color = "#566573", size = 5, hjust = 0)

ggsave("Plots/Somerset/Somerset_Snipe.png", width = 18, height = 20, units = "cm")


## Plot Whether field were surveyed or not
SL_Base + geom_sf(data = SomersetCov, mapping = aes(geometry = geometry, fill = as.character(Covered)), colour = NA) +
          scale_fill_manual(values = c("#9947F7", "#47F799"), labels = c("N", "Y"), name = "Surveyed?") +
          annotate(geom = "text", x = st_bbox(SomersetSurv_Sni$geometry)[1]-1000, y = st_bbox(SomersetSurv_Sni$geometry)[4]-1400,
                   label = "Coverage", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Somerset/Somerset_Coverage.png", width = 18, height = 20, units = "cm")




##------------------------------##
#### 5.3 Norfolk Broads plots ####
##------------------------------##

## Filter out the data for the Broads
BroadsSurv <- filter(SurveySum, Landscape == "Broads")
BroadsCov <- filter(SurveyCov, Landscape == "Broads")
Broads <- filter(Lanscap, Att3Value == "Broads")

## Create a base plot that I can then add field data on top of
NB_Base <- ggplot() + geom_sf(data = Broads, mapping = aes(geometry = geometry), fill = "#f5f5f2", colour = "NA") +
           geom_sf(data = Broads, mapping = aes(geometry = geometry), fill = NA, colour = "orange") +
           coord_sf(xlim = c(st_bbox(BroadsSurv$geometry)[1]-2000, st_bbox(BroadsSurv$geometry)[3]+2000), 
                    ylim = c(st_bbox(BroadsSurv$geometry)[2]-2000, st_bbox(BroadsSurv$geometry)[4]+2000), crs = 27700,
                    expand = F) +
           theme_minimal() +
           theme(text = element_text(family = "Karla", color = "#2D2D2E"), # change the font
            legend.position = c(0.14,0.15), # Move the legend to the bottom left corner
            legend.margin = margin(t = 0, unit = "cm"), # Set top of legend margin to zero
            legend.title = element_text(size = 13), # Adjust legend text size
            legend.text = element_text(size = 10, color = "#878585"), # Adjust legend text size
            legend.key.width = unit(1, "cm"), # Change legend colorbar width
            legend.key.height = unit(0.5, "cm"), # Change legend colorbar height
            legend.background = element_blank(), # Remove legend background
            axis.title = element_blank(), # Remove axis titles
            panel.grid = element_line(color = "#ebebe5", linewidth = 0.2), # Change lat/long graticule line colors
            panel.background = element_rect(fill = "white", color = NA)) +  # Change plot panel background color
            annotate(geom = "text", x = st_bbox(BroadsSurv$geometry)[1], y = st_bbox(BroadsSurv$geometry)[4]-500,
                     label = "Broads", family = "Karla", color = "#2D2D2E",
                     size = 9, hjust = 0)



## Plot the number of visits to each field
## create colour pallet so values are easy to tell apart 
my_colors <- RColorBrewer::brewer.pal(6, "Blues")[c(2,3,5,6)]

## create the plot for number of visits per field
NB_Base + geom_sf(data = BroadsSurv, mapping = aes(geometry = geometry, fill = as.character(N_visits)), colour = NA) +
          scale_fill_manual(values = my_colors, labels = c("1", "2","3", "4"), name = "# Visits") +
          annotate(geom = "text", x = st_bbox(BroadsSurv$geometry)[1], y = st_bbox(BroadsSurv$geometry)[4]-3500,
                   label = "Visit number", family = "Karla", color = "#566573", size = 6, hjust = 0) +
          theme(legend.position = c(0.1,0.1))

ggsave("Plots/Broads/Broads_Visits.png", width = 18, height = 20, units = "cm")


## Plot estimated number of wader pairs
## Lapwing breeding pairs estimates
BroadsSurv <- BroadsSurv |> 
  mutate(est_pairsL = ifelse(est_pairsL == 0, NA, est_pairsL))
NB_Base + geom_sf(data = BroadsSurv, mapping = aes(geometry = geometry, fill = est_pairsL), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Lapwing Pairs")  +
          annotate(geom = "text", x = st_bbox(BroadsSurv$geometry)[1], y = st_bbox(BroadsSurv$geometry)[4]-3500,
                   label = "Breeding Lapwing", family = "Karla", color = "#566573", size = 4.5, hjust = 0)

ggsave("Plots/Broads/Broads_Lapwing.png", width = 18, height = 20, units = "cm")

## Redshank breeding pairs estimates
BroadsSurv <- BroadsSurv |> 
  mutate(est_pairsR = ifelse(est_pairsR == 0, NA, est_pairsR))
NB_Base + geom_sf(data = BroadsSurv, mapping = aes(geometry = geometry, fill = est_pairsR), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Redshank Pairs") +
          annotate(geom = "text", x = st_bbox(BroadsSurv$geometry)[1], y = st_bbox(BroadsSurv$geometry)[4]-3500,
                   label = "Breeding Redshank", family = "Karla", color = "#566573", size = 4, hjust = 0)

ggsave("Plots/Broads/Broads_Redshank.png", width = 18, height = 20, units = "cm")


## Plot Whether field were surveyed or not
NB_Base + geom_sf(data = BroadsCov, mapping = aes(geometry = geometry, fill = as.character(Covered)), colour = NA) +
          scale_fill_manual(values = c("#9947F7", "#47F799"), labels = c("N", "Y"), name = "Surveyed?") +
          annotate(geom = "text", x = st_bbox(BroadsSurv$geometry)[1], y = st_bbox(BroadsSurv$geometry)[4]-3500,
                   label = "Coverage", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Broads/Broads_Coverage.png", width = 18, height = 20, units = "cm")




##------------------------------##
#### 5.4 Greater Thames plots ####
##------------------------------##

## Filter out the data for the Broads
ThamesSurv <- filter(SurveySum, Landscape == "Greater Thames")
ThamesCov <- filter(SurveyCov, Landscape == "Greater Thames")
Thames <- filter(Lanscap, Att3Value == "Greater Thames")

## Create a base plot that I can then add field data on top of
GT_Base <- ggplot() + geom_sf(data = Thames, mapping = aes(geometry = geometry), fill = "#f5f5f2", colour = "NA") +
           geom_sf(data = Thames, mapping = aes(geometry = geometry), fill = NA, colour = "orange") +
           coord_sf(xlim = c(st_bbox(ThamesSurv$geometry)[1]-250, st_bbox(ThamesSurv$geometry)[3]+250), 
                    ylim = c(st_bbox(ThamesSurv$geometry)[2]-250, st_bbox(ThamesSurv$geometry)[4]+250), crs = 27700,
                    expand = F) +
           theme_minimal() +
           theme(text = element_text(family = "Karla", color = "#2D2D2E"), # change the font
            legend.position = c(0.9,0.15), # Move the legend to the bottom left corner
            legend.margin = margin(t = 0, unit = "cm"), # Set top of legend margin to zero
            legend.title = element_text(size = 13), # Adjust legend text size
            legend.text = element_text(size = 10, color = "#878585"), # Adjust legend text size
            legend.key.width = unit(1, "cm"), # Change legend colorbar width
            legend.key.height = unit(0.5, "cm"), # Change legend colorbar height
            legend.background = element_blank(), # Remove legend background
            axis.title = element_blank(), # Remove axis titles
            panel.grid = element_line(color = "#ebebe5", size = 0.2), # Change lat/long graticule line colors
            panel.background = element_rect(fill = "white", color = NA)) +  # Change plot panel background color
            annotate(geom = "text", x = st_bbox(ThamesSurv$geometry)[1], y = st_bbox(ThamesSurv$geometry)[4]-1000,
                     label = "Greater Thames", family = "Karla", color = "#2D2D2E",
                     size = 9, hjust = 0)



## Plot the number of visits to each field
## create colour pallet so values are easy to tell apart 
my_colors <- RColorBrewer::brewer.pal(6, "Blues")[c(2,3,5,6)]

## create the plot for number of visits per field
GT_Base + geom_sf(data = ThamesSurv, mapping = aes(geometry = geometry, fill = as.character(N_visits)), colour = NA) +
          scale_fill_manual(values = my_colors, labels = c("1", "2","3", "4"), name = "# Visits") +
          annotate(geom = "text", x = st_bbox(ThamesSurv$geometry)[1]+200, y = st_bbox(ThamesSurv$geometry)[4]-5000,
                   label = "Visit number", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Thames/Thames_Visits.png", width = 18, height = 20, units = "cm")


## Plot estimated number of wader pairs
## Lapwing breeding pairs estimates
ThamesSurv <- ThamesSurv |> 
  mutate(est_pairsL = ifelse(est_pairsL == 0, NA, est_pairsL))
GT_Base + geom_sf(data = ThamesSurv, mapping = aes(geometry = geometry, fill = est_pairsL), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Lapwing Pairs")  +
          annotate(geom = "text", x = st_bbox(ThamesSurv$geometry)[1]+200, y = st_bbox(ThamesSurv$geometry)[4]-5000,
                   label = "Breeding Lapwing", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Thames/Thames_Lapwing.png", width = 18, height = 20, units = "cm")

## Redshank breeding pairs estimates
ThamesSurv <- ThamesSurv |> 
  mutate(est_pairsR = ifelse(est_pairsR == 0, NA, est_pairsR))
GT_Base + geom_sf(data = ThamesSurv, mapping = aes(geometry = geometry, fill = est_pairsR), colour = NA) +
          scale_fill_gradient(low = "#A6D1FB", high = "darkblue", na.value = "#E5E7E9") +
          labs(fill = "Redshank Pairs") +
          annotate(geom = "text", x = st_bbox(ThamesSurv$geometry)[1]+200, y = st_bbox(ThamesSurv$geometry)[4]-5000,
                   label = "Breeding Redshank", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Thames/Thames_Redshank.png", width = 18, height = 20, units = "cm")


## Plot Whether field were surveyed or not
GT_Base + geom_sf(data = ThamesCov, mapping = aes(geometry = geometry, fill = as.character(Covered)), colour = NA) +
          scale_fill_manual(values = c("#9947F7", "#47F799"), labels = c("N", "Y"), name = "Surveyed?") +
          annotate(geom = "text", x = st_bbox(ThamesSurv$geometry)[1], y = st_bbox(ThamesSurv$geometry)[4]-5000,
                   label = "Coverage", family = "Karla", color = "#566573", size = 6, hjust = 0)

ggsave("Plots/Thames/Thames_Coverage.png", width = 18, height = 20, units = "cm")





##------------------------------------------##
#### 6. Summaries Regional Survey Outputs ####
##------------------------------------------##

##** Summary of the number of visits per surveyed field in each landscape
VisitSum <- SurveySum |> 
  filter(is.na(Landscape)==F & Somerset %in% c("All", "S")) |> 
  mutate(Visit1 = ifelse(N_visits == 1, 1, 0),
         Visit2 = ifelse(N_visits == 2, 1, 0),
         Visit3 = ifelse(N_visits == 3, 1, 0),
         Visit4 = ifelse(N_visits == 4, 1, 0)) |> 
  group_by(Landscape) |> 
  summarise(`One Visit` = sum(Visit1),
            `Two Visits` = sum(Visit2),
            `Three Visits` = sum(Visit3),
            `Four Visits` = sum(Visit4))

## Make into a table and save
VisitSum |> 
  gt() |> 
  tab_header(title = md("**Visits per Field**")) |> 
  opt_stylize(style = 3, color = "blue") |> 
  gtsave(filename = "Plots/Summary_Tables/Number_of_Visits.png")


##** Survey of the number of fields surveyed and those not
CovSum <- SurveyCov |> 
  filter(is.na(Landscape)==F) |> 
  mutate(Yes = ifelse(Covered == "Y", 1, 0),
         No = ifelse(Covered == "N", 1, 0)) |> 
  group_by(Landscape) |>   
  summarise(Covered = sum(Yes),
            `Not Covered` = sum(No))

## Make into a table and save
CovSum |> 
  gt() |> 
  tab_header(title = md("**Survey Coverage**")) |> 
  opt_stylize(style = 3, color = "green") |> 
  gtsave(filename = "Plots/Summary_Tables/Survey_coverage.png")
  

##** The total number of wader pairs in each landscape
WaderSum <- SurveySum |> 
  filter(is.na(Landscape)==F) |> 
  group_by(Landscape) |> 
  summarise(Lapwing = sum(est_pairsL, na.rm = T), 
            Redshank = sum(est_pairsR, na.rm = T),
            Snipe = sum(est_pairsS, na.rm = T),
            Oyc = sum(est_pairsO, na.rm = T),
            Curlew = sum(est_pairsC, na.rm = T))
  
## Make into a table and save
WaderSum |> 
  select(-Oyc) |> 
  mutate(across(c(Lapwing, Redshank, Snipe, Curlew),\(x) round2(x, digits =0))) |> 
  gt() |> 
  tab_header(title = md("**Breeding Pair Estimates**")) |> 
  opt_stylize(style = 3, color = "blue") |> 
  gtsave(filename = "Plots/Summary_Tables/Breeding_Pairs.png")




