## Load in packages
pacman::p_load(tidyverse, data.table, sf, terra)
options(scipen=999) # turn off scientific notation


## Read in priority landscape polygons
Pr <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp") |>  
  select(Att3Value) |>  rename(Landscape = Att3Value)
Pr <- filter(Pr, Landscape %in% c("Broads", "Somerset Levels and Moors", "North Kent", "Essex")) 
plot(Pr$geometry)
## Buffer priority landscapes
Pr <- st_buffer(Pr, dist = 10000)
plot(Pr$geometry)



## Read in the UKCEH land parcels data set from 2021
UKC <- st_read("RawData/LandCover/UKCEH_LCvector_2021.shp") |> select(OBJECTID, X_mode)
glimpse(UKC)

## Work out which fields fall within any of the buffered priority landscapes
OverlapCEH <- st_join(UKC, Pr)
glimpse(OverlapCEH)
OverlapCEH <- OverlapCEH |> filter(is.na(Landscape) ==F)
write_sf(OverlapCEH, "RawData/LandCover/Crop_UKCEK_LCvector/UKCEH_LCvector_2021.shp")
rm(UKC, BWWM); gc()




## Read in the UKCEH land parcels data set from 2021
NEHab <- st_read("RawData/NEPriorityHabs/Priority_Habitat_Inventory_England.shp")
glimpse(NEHab)

## Work out which fields fall within any of the buffered priority landscapes
OverlapNEHab <- st_join(NEHab, Pr)
#glimpse(OverlapNEHab)
OverlapNEHab <- OverlapNEHab |> filter(is.na(Landscape) ==F) |> select(-Landscape)
write_sf(OverlapNEHab, "RawData/NEPriorityHabs/Crop/Priority_Habitat_Inventory_England.shp")





Rast <- rast("RawData/LandCover/gblcm2023_25m.tif")

## crop the raster data to the extent of the survey data (helps with speed)
Rastcrop <- crop(Rast, vect(Pr))
#plot(Rastcrop)

terra::writeRaster(Rastcrop, "RawData/LandCover/Crops/gblcm2023_25m.tif")


