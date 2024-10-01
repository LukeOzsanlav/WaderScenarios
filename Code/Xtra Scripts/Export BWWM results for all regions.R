##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
##  
##------------------------------------------------------##


## Load in required packages
pacman::p_load(here, tidyverse, data.table, sf, terra, gt)
source("Code/Helper functions.R") # get helper functions

## Use `here` package for specifying file paths
here::i_am("Code/2 - Field Selection + Plotting.R")
library(here)
here()




##---------------------##
#### 0. Data Read In ####
##---------------------##

##-----------------------------##
#### 0.1 Read in survey data ####
##-----------------------------##

## read in cleaned survey data output in "Cleaning Script.R"
Wader <- read_csv(here("CleanData", "Script 1", "Estimated_Pairs_clean.csv"))

## Optional Step to filter all date + the additional lapwing fields for Somerset 
#Wader <- filter(Wader, Somerset %in% c("All", "L"))



####----------------------------##
#### 0.2 Read in spatial data ####
##------------------------------##

## Read in high res world UK coastline (for plotting to give context)
UK <- st_read(here("RawData", "UK_Coastline", "UK_Coatline.shp"))

## Read in the field shapes GIS layer
FieldShapes <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp"))

## Read in RSPB priority landscapes
Lanscap <- st_read(here("RawData", "Priority Landscapes", "EnglandWales_PriorityLandscapes.shp"))
## Filter out the landscapes I am focusing on
Lanscap <- filter(Lanscap, Att3Value %in% c("Somerset Levels and Moors", "Greater Thames", "Broads", "Suffolk Coast"))

## Read in the ESAs
ESAs <- st_read(here("RawData", "ESAs", "Wader_focused_ESAs.shp"))

## Read in the UKCEH landcover data
LC <- rast(here("RawData", "LandCover", "gblcm25m2021.tif"))
LC <- LC[[1]]



##-----------------------------------##
#### 1. Spatial Data Manipulations ####
##-----------------------------------##

##---------------------------##
#### 1.1 Plot spatial data ####
##---------------------------##

# ## plot field shapes
# p1 <- ggplot() + geom_sf(data = FieldShapes, mapping = aes(geometry = geometry), fill = "orange") +
#          geom_sf(data = UK, mapping = aes(geometry = geometry), fill = NA) +
#          theme_light()
# 
# ## plot priority landscapes
# p2 <- ggplot() + geom_sf(data = Lanscap, mapping = aes(geometry = geometry, fill = Att3Value), colour = NA) +
#          geom_sf(data = UK, mapping = aes(geometry = geometry), fill = NA) +
#          coord_sf(xlim = c(327378.1-10000, 653935.4+10000), ylim = c(120152.9-10000, 330111.0+10000), crs = 27700,
#                   expand = FALSE) +
#          theme_light()
# 
# ## plot ESAs
# p3 <- ggplot() + geom_sf(data = ESAs, mapping = aes(geometry = geometry, fill = NAME), colour = NA) +
#          geom_sf(data = UK, mapping = aes(geometry = geometry), fill = NA) +
#          coord_sf(xlim = c(327378.1-10000, 653935.4+10000), ylim = c(120152.9-10000, 330111.0+10000), crs = 27700,
#                   expand = FALSE) +
#          theme_light()



##----------------------------------##
#### 1.2 Create Landscape Buffers ####
##----------------------------------##

##** Combine and buffer Priority Landscapes and ESAs
## Create bounds for Somerset area
Somerset_Prior <- filter(Lanscap, Att3Value == "Somerset Levels and Moors")
#Somerset_ESA <- st_combine(filter(ESAs, NAME == "SOMERSET LEVELS AND MOORS"))
Somerset_buf <- Somerset_Prior
Somerset_buf <-st_buffer(Somerset_buf, dist = 500)
plot(Somerset_buf$geometry)

## Create bounds for Norfolk Broads
Broads_Prior <- filter(Lanscap, Att3Value == "Broads")
#Broads_ESA <- st_combine(filter(ESAs, NAME == "BROADS"))
Broads_buf <- Broads_Prior
Broads_buf <-st_buffer(Broads_buf, dist = 500)
plot(Broads_buf$geometry)

## Create bounds for Greater Thames
Thames_Prior <- filter(Lanscap, Att3Value == "Greater Thames")
#Thames_ESA <- st_combine(filter(ESAs, NAME %in% c("ESSEX COAST", "NORTH KENT MARSHES")))
Thames_buf <- Thames_Prior
plot(Thames_buf$geometry)
Thames_buf <-st_buffer(Thames_buf, dist = 500)
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



##----------------------------------------##
#### 1.4 Extract fields from Landscapes ####
##----------------------------------------##

## Work out if fields fall within any of the priority landscapes
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





##-----------------------------------------------------##
#### 2. Join Pair Estimates and Spatial Calculations ####
##-----------------------------------------------------##

## Trim Fields data set for join
FieldSub <- Fields |> select(F_LOC_ID, FieldArea, FieldX, FieldY, Landscape)

## Join Spatial calculation and wader survey data
PreL <- nrow(Wader)
Wader <- left_join(Wader, FieldSub, by = "F_LOC_ID")
nrow(Wader) == PreL




##-------------------------------------##
#### 3. Remove Non-grassland Fields  ####
##-------------------------------------##

## Remove any fields that are not in a priority landscape buffer
PFields <- filter(Wader, is.na(Landscape)==F)

## Remove fields that are the not grassland based off the UK CEH data set
## Firstly, crop the data set to the south of England
LC <- crop(LC, vect(Lanscap_buf))
plot(LC)

## Calculate the dominant habitat class for each field
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]} # function to calculate the mode

## extract the dominant habitat from the UKCEH raster
MainHab <- extract(LC, vect(PFields$geometry), fun = getmode)

## Add these habitat to main data set
PFields$DomHab <- MainHab$gblcm25m2021_1
table(PFields$DomHab)


## Calculate proportion of each that is grassland
## Reclassify raster so grassland is 1 and non-grassland is 0
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

## select columns to send to Damon and give DomHab habitat column meaningful names
PFields <- PFields |>  select(c(F_LOC_ID, year, Landscape, est_pairsL, est_pairsR, est_pairsS, est_pairsC, FieldArea, DomHab, PropGrass, Unsuit_Hab, geometry)) |> 
  mutate(DomHab = ifelse(DomHab < 3, "Woodland",
                          ifelse(DomHab == 3, "Arable",
                                  ifelse(DomHab == 4, "Improved Grass",
                                          ifelse(DomHab == 5, "Neutral Grass", 
                                                  ifelse(DomHab == 8, "Fen",
                                                         ifelse(DomHab == 13, "Water",
                                                                ifelse(DomHab == 14, "Water",
                                                                       ifelse(DomHab == 19, "Saltmarsh",
                                                                              ifelse(DomHab >= 20, "Urban", NA))))))))))



##---------------------------------##
#### 4. Add on extra field data  ####
##---------------------------------##

## First, Add column whether waders are present
PFields <- PFields %>%  
  mutate(WaderTot = rowSums(dplyr::select(., est_pairsL, est_pairsC, est_pairsR, est_pairsS), na.rm = T),
         Wader_occ = ifelse(WaderTot>0, "Y", ifelse(is.na(WaderTot)== T, "N", "N"))) |> 
  select(-WaderTot)

## Read in the habitat data from the field surveys
FieldHab <- read_csv(here("CleanData", "Script 1", "Field_Characteristics_Clean.csv"))

## select the columns needed from the field survey data
FieldHab <- FieldHab |> select(c(F_LOC_ID, year, STANDING_WATER_TOTAL_PERCENT,
                                 TALL_BOUNDARY_PERCENT, VEG_HEIGHT, VEG_STRUCTURE,
                                 RUSH_PERCENT, RUSH_DISTRIBUTION, GRASSLAND_TYPE, OTHER_HABITAT))

## Add this onto the bird surveys data
PFields <- left_join(PFields, FieldHab, by = c("F_LOC_ID", "year"))
PFields <- PFields |> mutate(OTHER_HABITAT = ifelse(OTHER_HABITAT == "G", NA, OTHER_HABITAT))

## write out this file
st_write(PFields, here("CleanData", "All regions", "All_regions_BWWM_results.shp"))
#st_write(PFields, here("CleanData", "All regions", "All_regions_BWWM_results_Lapwing.shp"))



##------------------------------------------##
#### 6. Summaries Regional Survey Outputs ####
##------------------------------------------##

##** The total number of wader pairs in each landscape
WaderSum <- PFields |> 
  group_by(Landscape) |> 
  summarise(Lapwing = sum(est_pairsL, na.rm = T), 
            Redshank = sum(est_pairsR, na.rm = T),
            Snipe = sum(est_pairsS, na.rm = T),
            Curlew = sum(est_pairsC, na.rm = T))
  
## Make into a table and save
WaderSum |> 
  mutate(across(c(Lapwing, Redshank, Snipe, Curlew),\(x) round2(x, digits =0))) |> 
  gt() |> 
  tab_header(title = md("**Breeding Pair Estimates**")) |> 
  opt_stylize(style = 3, color = "blue") |> 
  gtsave(filename = "CleanData/All regions/Breeding_Pairs.png")




