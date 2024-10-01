## Created 27/03/2024
## Aim: Create derived data set needed for various guidelines
## Currently this script takes multiple years of UKCEH data and works out field which current arable fields were historically grassland 


pacman::p_load(tidyverse, sf, terra)
here::i_am("Code/7- Spatially Map Guidelines.R")
library(here)
here()




##
#### Read in Datasets ####
##

## Load in priority landscapes
PL <- st_read("RawData/Priority Landscapes/EnglandWales_PriorityLandscapes.shp")
plot(PL)
Broads <- filter(PL, Att3Value == "Broads")


## Load in Landcover data sets

## 2021 data set
LC2021 <- rast(here("RawData", "LandCover", "gblcm25m2021.tif"))
LC2021 <- LC2021[[1]]
## crop the UK CEH data to polgon of interesr
LC2021 <- crop(LC2021, vect(st_buffer(st_as_sf(Broads), dist = 500)))
plot(LC2021)
crs(LC2021)


## 2007 data set
LC2007 <- rast(here("RawData", "LandCover", "lcm2007gb25m.tif"))
LC2007 <- LC2007[[1]]
## crop the UK CEH data to polgon of interesr
LC2007 <- crop(LC2007, vect(st_buffer(st_as_sf(Broads), dist = 500)))
plot(LC2007)
crs(LC2007)

## 2000 data set, need to re download this data set as it will not read in
LC2000 <- rast(here("RawData", "LandCover", "LCM2000_GB_25m.tif"))
LC2000 <- LC2000[[1]]
## crop the UK CEH data to polgon of interesr
LC2000 <- crop(LC2000, vect(st_buffer(st_as_sf(Broads), dist = 500)))
plot(LC2000)
unique(values(LC2000))


## 1990 data set
LC1990<- rast(here("RawData", "LandCover", "gb1990lcm25m.tif"))
LC1990 <- LC1990[[1]]
## crop the UK CEH data to polgon of interesr
LC1990 <- crop(LC1990, vect(st_buffer(st_as_sf(Broads), dist = 500)))
plot(LC1990)
crs(LC1990)



##
#### reclassify Rasters
##


## 2021 raster, reclassify so Arable is 1 and all other habitats are 0
m <- c(0, 2.5, 0,
       2.5, 3.5, 1,
       3.5, 22, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC2021_Arable <- classify(LC2021, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC2021_Arable)


## 2007 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 3.5, 0,
       3.5, 9.5, 1,
       9.5, 23, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC2007_Grass <- classify(LC2007, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC2007_Grass)


## 2000 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 49, 0,
       49, 92, 1,
       92, 110, 0,
       110, 112, 1, 
       112, 250, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC2000_Grass <- classify(LC2000, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC2000_Grass)


## 1990 raster, reclassify so grassland is 1 and all other habitats are 0
m <- c(0, 3.5, 0,
       3.5, 8.5, 1,
       8.5, 23, 0) # matrix for re-classification
rclmat <- matrix(m, ncol=3, byrow=TRUE) # make matrix
LC1990_Grass <- classify(LC1990, rclmat, include.lowest=TRUE) # reclassify the raster
plot(LC1990_Grass)



##
#### Combine Rasters ####
##

## Add together the two grassland layers
All_Grass <- LC2007_Grass+LC1990_Grass+LC2000_Grass
plot(All_Grass)

## Multiple together the arable and grassland layers
## ANy pixel with a value of 1/2 was Arable in 2021 but grassland in 1 or 2 years previously
ArableHist <- LC2021_Arable*All_Grass
plot(ArableHist)




