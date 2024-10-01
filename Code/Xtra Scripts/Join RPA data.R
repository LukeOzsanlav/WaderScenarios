pacman::p_load(tidyverse, data.table, sf, terra)


## Read in the csv from the RPA that contains FarmIDs for each parcel ID
RPA_Farms <- fread("RawData/RPA/RPA_ParcelList_20240424&Suffolk.csv")

## Read in the parcel IDs from from area of interest
## https://www.data.gov.uk/dataset/ed6f8d0e-df4f-4027-b4e3-ee43d7cebf35/rpa-parcel-points-england#:~:text=Parcel%20Points%20is%20a%20simple,or%20attach%20their%20own%20attributes.
RPA_parcels <- st_read("RawData/RPA/rpa_parcel_points_East_England/rpa_parcel_pointsPoint.shp")
RPA_parcels2 <- st_read("RawData/RPA/rpa_parcel_points_Somerset/rpa_parcel_pointsPoint.shp")
RPA_parcels <- bind_rows(RPA_parcels2, RPA_parcels) # bind together East England and SOmerset

## These are the columns that match up bteween the data sets
length(unique(RPA_Farms$NGC))
length(unique(RPA_parcels$parcel_ref))

## Can join the two data sest using either the NGC or ID columns
## joining using NGC and parcel_ref returns more points so perhaps use this for now
Join <- inner_join(RPA_Farms, RPA_parcels, by = join_by(NGC==parcel_ref))
plot(Join$geometry)
write_sf(Join, "CleanData/Inter_farms.shp")

