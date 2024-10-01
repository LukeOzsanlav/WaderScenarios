pacman::p_load(mapdeck, sf, tidyverse, here, geojsonsf)


## set mapbeck token
key <- 'pk.eyJ1Ijoiam96c2FsMSIsImEiOiJjbHdnazA1ZTMwNGwxMnFybXZoOG8wd3ZsIn0.dfHN-YRHHczxvDafxD14vQ'    ## put your own token here
mapdeck(token = key)
set_token(key)


## Read in filtered breeding pairs estimates
Waders <- read_csv(here("CleanData", "Script 3", "Breeding_Pairs_Attribs.csv")) |> select(-geometry)

## read in shaepfile for BWWM fields and join to breeding waders pairs
BWWM <- st_read(here("RawData", "BWWM field shapefile", "BWWM_fields_27Sep2023.shp")) |> select(F_LOC_ID)
Waders <- left_join(Waders, BWWM, by = "F_LOC_ID")
Waders <- st_as_sf(Waders)


som <- filter(Waders, Landscape == "Somerset Levels and Moors") |> st_transform(crs = 4326)
som <- mutate(som, est_pairsR10 = est_pairsR*100)



mapdeck(token = key, style = mapdeck_style("dark")) %>%
  add_polygon(
    data = som
    , layer = "polygon_layer"
    , fill_colour = "est_pairsR"
    , elevation = "est_pairsR10"
  )


table(Waders$Landscape)
Br <- filter(Waders, Landscape == "Broads") |> st_transform(crs = 4326)
Br <- mutate(Br, est_pairsR10 = est_pairsR*100)



mapdeck(token = key, style = mapdeck_style("dark")) %>%
  add_polygon(
    data = Br
    , layer = "polygon_layer"
    , fill_colour = "est_pairsR"
    , elevation = "est_pairsR10"
  )
