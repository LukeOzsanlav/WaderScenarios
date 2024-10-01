## Luke Ozsanlav-Harris

## Combine all classified acc datasets and create an Rdata file
## Do some check along the way too

library(tidyverse)
library(data.table)

## read in all the classified dats sets
Data <-list.files(path = "Classified datasets") %>% map_df(~fread(paste0("Classified datasets/", .)))

