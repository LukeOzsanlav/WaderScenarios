##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 28/08/2024
## 
## 
## Aim: Plot packages for manuscript on scenario modelling
## 
## 
##------------------------------------------------------##

## Read in packages and custom function
pacman::p_load(tidyverse, sf, terra, tidyterra, ggspatial, ggpubr, tidyterra)
options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R")
source("Code/Scenarios/5.2- Scenario Functions.R")


##---------------------------------------##
#### 0. Read in data sets for plotting ####
##---------------------------------------##

## Read in priority landscape boundary
MyBoxes <- st_read("RawData/Priority Landscapes/Split Priority Landscapes/EnglandWales_PriorityLandscapes.shp")

## Read in outline of UK
UKCoast <- st_read("RawData/UK_Coastline/UK_Coatline.shp")
UKCoastHighRes <- st_read("RawData/UK_Coastline/Countries Boundaries/UK_Country_Boundaries.shp")



##-- READ IN & PREPARE CANVAS --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
        rename(FieldArea=ParcArea)






##-----------------------------##
## 0.2 General plot styling ####
##-----------------------------##

## Set a general theme so that all pots look similar
GeneralThemeing <- theme(
          axis.title.x = element_text(size = 14),
          #axis.text.x = element_text(hjust=0.7, angle = 45),
          axis.title.y = element_text(angle=90, vjust = 0.4, size = 14),
          axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
          text = element_text(color = "#2D2D2E"), 
          panel.grid = element_line(color = "#ebebe5", linewidth = 0.2),
          panel.background = element_rect(fill = "#f5f5f2", color = NA),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

## Set a general theme so that all pots look similar
GeneralThemeing2 <- theme(
  axis.title.x = element_text(size = 11),
  #axis.text.x = element_text(hjust=0.7, angle = 45),
  axis.title.y = element_text(angle=90, vjust = 0.4, size = 11),
  axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
  text = element_text(color = "#2D2D2E"), 
  panel.grid = element_line(color = "#ebebe5", linewidth = 0.2),
  panel.background = element_rect(fill = "#f5f5f2", color = NA),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14))





##-----------------------------------##
#### 1.0 Plot categories in canvas ####
##-----------------------------------##

##-- READ IN & PREPARE BROADS CANVAS --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
  rename(FieldArea=ParcArea)

## Filter out the Broads Landscape boundary
Broads <- MyBoxes |> filter(Att3Value == "Broads") 
## crop to area just around the Broads priority landscape
BroadsCoast <- st_transform(UKCoast, crs = st_crs(Broads)) |> st_crop(Broads |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Broads)[1]-2000, ext(Broads)[2]+2000), ylim = c(ext(Broads)[3]-2000, ext(Broads)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## filter the canvas so that I just have land that is in scope
CanvCat <- filter(Canv, !Category %in% c("NoOpp", "Masked Hab")) |> 
           mutate(Category = ifelse(Category== "AES Only", "Wader-AES", Category))

## make plot of the different categories in the landscape
Categories <- ggplot() +

  ## Add the coastline
  geom_sf(data=BroadsCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Broads, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
  scale_fill_manual(values = c("#DDB24C", "#4CDD6A", "#4C77DD", "#DD4CC0"), name="Opportunity\nLand Type") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  theme(legend.position = "bottom") +
  GeneralThemeing2


## Add a column to indicate the opportunity areas for different Lawton principles
CanvCat$WaderSite <- ifelse(is.na(CanvCat$ClustGroup)==F, "Better", "Bigger/More")

## Plot the opportunity areas for different Lawton principles
Lawton <- ggplot() +

  ## Add the coastline
  geom_sf(data=BroadsCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Broads, mapping=aes(geometry=geometry),  colour = "#273746", fill = NA) + 
  #scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = WaderSite), colour = NA) +
  scale_fill_manual(values = c("#e9cd50", "#506ce9"), name="Location\nStrategy") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  #annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  theme(legend.position = "bottom") +
  GeneralThemeing2


## Arrange the above two plot together
P1 <- ggarrange(Categories, Lawton, ncol = 2, labels="auto")

## save the arranged plot
ggsave(plot = P1, filename = "CleanData/Paper Plots/Canvas Examps/Canvas_Combo_Broads.png", width = 35, height = 20, units = "cm")



##-- READ IN & PREPARE ESEEX CANVAS --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
  rename(FieldArea=ParcArea)

## Filter out the Essex Landscape boundary
Essex <- MyBoxes |> filter(Att3Value == "Essex") 

## crop to area just around the Essex priority landscape
EssexCoast <- st_transform(UKCoast, crs = st_crs(Essex)) |> st_crop(Essex |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Essex)[1]-2000, ext(Essex)[2]+2000), ylim = c(ext(Essex)[3]-2000, ext(Essex)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## filter the canvas so that I just have land that is in scope
CanvCat <- filter(Canv, !Category %in% c("NoOpp", "Masked Hab")) |> 
  mutate(Category = ifelse(Category== "AES Only", "Wader-AES", Category))

## make plot of the different categories in the landscape
Categories <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=EssexCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Essex, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
  scale_fill_manual(values = c("#DDB24C", "#4CDD6A", "#4C77DD", "#DD4CC0"), name="Opportunity\nLand Type") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Add a column to indicate the opportunity areas for different Lawton principles
CanvCat$WaderSite <- ifelse(is.na(CanvCat$ClustGroup)==F, "Better", "Bigger/More")

## Plot the opportunity areas for different Lawton principles
Lawton <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=EssexCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Essex, mapping=aes(geometry=geometry),  colour = "#273746", fill = NA) + 
  #scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = WaderSite), colour = NA) +
  scale_fill_manual(values = c("#e9cd50", "#506ce9"), name="Location\nStrategy") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  #annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Arrange the above two plot together
P2 <- ggarrange(Categories, Lawton, ncol = 1, labels="auto")

## save the arranged plot
ggsave(plot = P2, filename = "CleanData/Paper Plots/Canvas Examps/Canvas_Combo_Essex.png", width = 20, height = 23, units = "cm")




##-- READ IN & PREPARE NORTH KENT CANVAS --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
  rename(FieldArea=ParcArea)

## Filter out the NorthKent Landscape boundary
NorthKent <- MyBoxes |> filter(Att3Value == "North Kent") 

## crop to area just around the NorthKent priority landscape
NorthKentCoast <- st_transform(UKCoast, crs = st_crs(NorthKent)) |> st_crop(NorthKent |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(NorthKent)[1]-2000, ext(NorthKent)[2]+2000), ylim = c(ext(NorthKent)[3]-2000, ext(NorthKent)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## filter the canvas so that I just have land that is in scope
CanvCat <- filter(Canv, !Category %in% c("NoOpp", "Masked Hab")) |> 
  mutate(Category = ifelse(Category== "AES Only", "Wader-AES", Category))

## make plot of the different categories in the landscape
Categories <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=NorthKentCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=NorthKent, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
  scale_fill_manual(values = c("#DDB24C", "#4CDD6A", "#4C77DD", "#DD4CC0"), name="Opportunity\nLand Type") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Add a column to indicate the opportunity areas for different Lawton principles
CanvCat$WaderSite <- ifelse(is.na(CanvCat$ClustGroup)==F, "Better", "Bigger/More")

## Plot the opportunity areas for different Lawton principles
Lawton <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=NorthKentCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=NorthKent, mapping=aes(geometry=geometry),  colour = "#273746", fill = NA) + 
  #scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = WaderSite), colour = NA) +
  scale_fill_manual(values = c("#e9cd50", "#506ce9"), name="Location\nStrategy") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  #annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Arrange the above two plot together
P3 <- ggarrange(Categories, Lawton, ncol = 1, labels="auto")

## save the arranged plot
ggsave(plot = P3, filename = "CleanData/Paper Plots/Canvas Examps/Canvas_Combo_NorthKent.png", width = 22, height = 18, units = "cm")



##-- READ IN & PREPARE SOMERSET CANVAS --##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.rds")

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               WaterCoverage = ifelse(is.na(WaterCoverage)==T, 0, WaterCoverage), # Any fields that did not get a water coverage value just get a zero
               Fence_Coverage = ifelse(Fence_Coverage == 0 | is.na(Fence_Coverage) == T, "N", "Y"), # match RF variable format
               ParcArea = m2_to_ha(ParcArea), # match RF variable format
               GroupArea = m2_to_ha(GroupArea), # convert to hectares
               Perim = st_length(st_cast(geometry,"MULTILINESTRING")),
               # initiate columns for wader abundance
               SnipeAbund = NA, 
               LapAbund = NA,
               RedAbund = NA) |> 
  rename(FieldArea=ParcArea)

## Filter out the Somerset Landscape boundary
Somerset <- MyBoxes |> filter(Att3Value == "Somerset Levels and Moors") 

## crop to area just around the Somerset priority landscape
SomersetCoast <- st_transform(UKCoast, crs = st_crs(Somerset)) |> st_crop(Somerset |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Somerset)[1]-2000, ext(Somerset)[2]+2000), ylim = c(ext(Somerset)[3]-2000, ext(Somerset)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## filter the canvas so that I just have land that is in scope
CanvCat <- filter(Canv, !Category %in% c("NoOpp", "Masked Hab")) |> 
  mutate(Category = ifelse(Category== "AES Only", "Wader-AES", Category))

## make plot of the different categories in the landscape
Categories <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=SomersetCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Somerset, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
  scale_fill_manual(values = c("#DDB24C", "#4CDD6A", "#4C77DD", "#DD4CC0"), name="Opportunity\nLand Type") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Add a column to indicate the opportunity areas for different Lawton principles
CanvCat$WaderSite <- ifelse(is.na(CanvCat$ClustGroup)==F, "Better", "Bigger/More")

## Plot the opportunity areas for different Lawton principles
Lawton <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=SomersetCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Somerset, mapping=aes(geometry=geometry),  colour = "#273746", fill = NA) + 
  #scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = WaderSite), colour = NA) +
  scale_fill_manual(values = c("#e9cd50", "#506ce9"), name="Location\nStrategy") +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  #annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2


## Arrange the above two plot together
P4 <- ggarrange(Categories, Lawton, ncol = 1, labels="auto")

## save the arranged plot
ggsave(plot = P4, filename = "CleanData/Paper Plots/Canvas Examps/Canvas_Combo_Somerset.png", width = 22, height = 26, units = "cm")










##-----------------------------------------##
#### 2.0 Plot scenario opportunity areas ####
##-----------------------------------------##

## Create a new column that indicated the opportunity areas for a scenario creating AES on unmanned grassland
CanvCat <- CanvCat |> 
           mutate(AES_Site = case_when(
                                Category == "Reserve"  ~ "No Opportunity",
                                Category == "AES Only"  ~ "No Opportunity",
                                Category == "Grass Opp" & is.na(CanvCat$ClustGroup)==T ~ "AES creation outside wader sites",
                                Category == "Grass Opp" & is.na(CanvCat$ClustGroup)==F~ "AES creation within wader sites",
                                .default = "No Opportunity"))


## make plot showing the opportunity areas for a scenario creating AES on unmanned grassland for all Lawton principles
AESOpportunity <- ggplot() +

  ## Add the coastline
  geom_sf(data=BroadsCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Broads, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name= "Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = AES_Site), colour = NA) +
  scale_fill_manual(values = c("#64ed90", "#ed64c1", "darkgrey"), name= "AES Opportunity Areas", 
                    labels=c("AES creation\noutside wader sites",
                             "AES creation\nwithin wader sites",
                             "No Opportunity")) +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing +
  guides(fill = guide_legend(order = 2), col = guide_legend(order = 1))




## Create a new column that indicated the opportunity areas for a scenario creating reserve on unmanned grassland or AES only land
CanvCat <- CanvCat |> 
           mutate(Res_Site = case_when(
                                Category == "Reserve"  ~ "No Opportunity",
                                Category %in% c("Grass Opp", "AES Only") & is.na(CanvCat$ClustGroup)==T ~ "Reserve creation outside wader sites",
                                Category %in% c("Grass Opp", "AES Only") & is.na(CanvCat$ClustGroup)==F~ "Reserve creation within wader sites",
                                .default = "No Opportunity"))

## make plot showing the opportunity areas for a scenario creating reserve on unmanned grassland or AES only land
ResOpportunity <- ggplot() +

  ## Add the coastline
  geom_sf(data=BroadsCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) + 
  
  ## Add the landscape boundary
  geom_sf(data=Broads, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name= "Landscape\nBoundary") +
  
  ## Add the canvas categories
  geom_sf(data = CanvCat, mapping = aes(geometry = geometry, fill = Res_Site), colour = NA) +
  scale_fill_manual(values = c("darkgrey", "#64ed90",  "#ed64c1"), name= "Reserve Opportunity\nAreas", 
                    labels=c("No Opportunity",
                             "Reserve creation\noutside wader sites",
                             "Reserve creation\nwithin wader sites")) +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing 





## Arrange the above two plot together
P2 <- ggarrange(AESOpportunity, ResOpportunity, ncol = 2, labels="auto", widths = c(1, 0.972))

## save the arranged plot
ggsave(plot = P2, filename = "CleanData/Paper Plots/Canvas Examps/Opportunity_Combo_Broads.png", width = 40, height = 20, units = "cm")








##---------------------------------------##
#### 3.0 Plot Reserve Management Curve ####
##---------------------------------------##


## Using my manage costs calculator function work out the costs of managing reserves of different sizes
## For original price and then adjusted for inflation
MaC1 <- data.frame(TotalHa = seq(from=50, to=1500, by = 10),
                   CostHa = Manage_Costs(TotalArea = seq(from=50, to=1500, by = 10), ParcelArea = 1, InflationAdjust=FALSE),
                   Year = "2011/12")
MaC2 <- data.frame(TotalHa = seq(from=50, to=1500, by = 10),
                   CostHa = Manage_Costs(TotalArea = seq(from=50, to=1500, by = 10), ParcelArea = 1, InflationAdjust=TRUE),
                   Year = "2023/24")

## combine data sets
MaC <- rbind(MaC1, MaC2)

## Create the plot
P3 = ggplot() +
  geom_line(data=MaC, mapping = aes(x= TotalHa, y= CostHa, group=Year, colour = Year), linewidth = 1.5) +
  scale_colour_manual(values = c("#4a88ea", "#eaac4a")) +
  xlab("Total area of land managed/ ha") +
  ylab("Cost (£-GDP) of management/ ha") +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top")

## save the arranged plot
ggsave(plot = P3, filename = "CleanData/Paper Plots/ManageCosts/Costs_of_management_ha.png", width = 25, height = 20, units = "cm")






##-----------------------------------##
#### 4.0 Plot Landscape Boundaries ####
##-----------------------------------##

## Filter out the landscapes of interest
Scapes <- filter(MyBoxes, Att3Value %in% c("Broads", "Essex", "North Kent", "Somerset Levels and Moors"))
Scapes$Att3Value <- ifelse(Scapes$Att3Value =="Somerset Levels and Moors", "Somerset Levels", Scapes$Att3Value)
Scapes$Att3Value <- ifelse(Scapes$Att3Value =="Broads", "Norfolk Broads", Scapes$Att3Value)
ScapeExt <- ext(Scapes)

## create the plot
LP1 <- ggplot() + 
  ## add coastline
 geom_sf(data = UKCoastHighRes, aes(geometry = geometry), color = "#ffffff", fill = "#CBCACA") +
 # add map of landscapes
 geom_sf(data = Scapes, aes(geometry = geometry, fill = Att3Value), color = NA) + 
 scale_fill_manual(name = "Landscape:",
                   values = c("#06d6a0", "#ffd166", "#26547c", "#ef476f")) +
 # set plot limits
 coord_sf(xlim = c(ScapeExt[1]-198000, ScapeExt[2]+15000), ylim = c(ScapeExt[3]-115000, ScapeExt[4]+80000), crs = 27700,
          expand = FALSE) +
 # add labels
 labs(x = "Longitude", y = "Latitude", fill = "Landscape:") +
 # add styling
 theme_light() +
 theme(legend.position = "bottom",
       legend.title = element_text(size = 13),
       legend.text = element_text(size = 12),
       axis.title.x = element_text(size = 15),
       axis.text.x = element_text(hjust=0.7, angle = 45),
       axis.title.y = element_text(angle=90, vjust = 0.4, size = 15),
       axis.text.y = element_text(hjust=0.7,angle=45,vjust=0.3),
       text = element_text(family = "Karla", color = "#2D2D2E"), 
       panel.grid = element_line(color = "#d9d9d4", size = 0.3),
       panel.background = element_rect(fill = "#f5f5f2", color = NA)) 


## Save the plot
ggsave("CleanData/Paper Plots/AllLandscapePlots.png", 
       plot = LP1,  width = 20, height = 20, units = "cm")


## create the plot
LP1 <- ggplot() + 
  ## add coastline
 geom_sf(data = UKCoastHighRes, aes(geometry = geometry), color = "#ffffff", fill = "#CBCACA") +
 # add map of landscapes
 geom_sf(data = Scapes, aes(geometry = geometry, fill = Att3Value), color = NA) + 
 scale_fill_manual(name = "Landscape:",
                   values = c("#06d6a0", "#ffd166", "#26547c", "#ef476f")) +
guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
 # set plot limits
 coord_sf(xlim = c(ScapeExt[1]-198000, ScapeExt[2]+10000), ylim = c(ScapeExt[3]-115000, ScapeExt[4]+80000), crs = 27700,
          expand = FALSE) +
 # add labels
 labs(x = "Longitude", y = "Latitude", fill = "Landscape:") +
 # add styling
 theme_light() +
 theme(legend.position = "bottom",
       legend.title = element_text(size = 28, face="bold"),
       legend.text = element_text(size = 24),
       axis.title = element_blank(),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       text = element_text(family = "Karla", color = "#2D2D2E"), 
       panel.grid = element_line(color = "#d9d9d4", size = 0.3),
       panel.background = element_rect(fill = "#f5f5f2", color = NA)) 


## Save the plot
ggsave("CleanData/Paper Plots/AllLandscapePlotsLARGETEXT.png", 
       plot = LP1,  width = 20, height = 20, units = "cm")





##-----------------------------------------##
#### 5.0 Plot standing water predictions ####
##-----------------------------------------##

## |- North kent ----

## Read in the raster of North Kent
NKWater <- rast("RawData/BinaryWaterRasters/NKent_BinaryWaterAll.tif")

## Filter out the landscapes of interest
NorthKent <- MyBoxes |> filter(Att3Value == "North Kent") 

## crop to area just around the NorthKent priority landscape
NorthKentCoast <- st_transform(UKCoast, crs = st_crs(NorthKent)) |> st_crop(NorthKent |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(NorthKent)[1]-2000, ext(NorthKent)[2]+2000), ylim = c(ext(NorthKent)[3]-2000, ext(NorthKent)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## Create plot of standing water cover
NKWa <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=NorthKentCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  
  ## Add the landscape boundary
  geom_sf(data=NorthKent, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the water raster
  geom_spatraster(data = as.factor(NKWater)) +
  scale_fill_viridis_d(name = "Pixel Classification:", option = "D", na.value = NA, 
                       labels = c("Dry grassland", "Standing water", "")) +
  guides(fill = guide_legend(na.translate = FALSE)) +

  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2 +
  theme(legend.position = "top")


## Save the plot
ggsave("CleanData/Paper Plots/NorthKent_StandingWater_Predicted.png", 
       plot = NKWa,  width = 22, height = 12, units = "cm")


## |- Somerset ----

## Read in the raster of North Kent
SomWater <- rast("RawData/BinaryWaterRasters/Som_BinaryWaterAll.tif")

## Filter out the landscapes of interest
Somer <- MyBoxes |> filter(Att3Value == "Somerset Levels and Moors") 

## crop to area just around the Somer priority landscape
SomerCoast <- st_transform(UKCoast, crs = st_crs(Somer)) |> st_crop(Somer |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Somer)[1]-2000, ext(Somer)[2]+2000), ylim = c(ext(Somer)[3]-2000, ext(Somer)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## Create plot of standing water cover
SomWa <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=SomerCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  
  ## Add the landscape boundary
  geom_sf(data=Somer, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the water raster
  geom_spatraster(data = as.factor(SomWater)) +
  scale_fill_viridis_d(name = "Pixel Classification:", option = "D", na.value = NA, 
                       labels = c("Dry grassland", "Standing water", "")) +
  guides(fill = guide_legend(na.translate = FALSE)) +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2 +
  theme(legend.position = "top")


## Save the plot
ggsave("CleanData/Paper Plots/Somerset_StandingWater_Predicted.png", 
       plot = SomWa,  width = 20, height = 20, units = "cm")



## |- Broads ----

## Read in the raster of North Kent
BroWater <- rast("RawData/BinaryWaterRasters/Broads_BinaryWaterAll.tif")

## Filter out the landscapes of interest
Broads <- MyBoxes |> filter(Att3Value == "Broads") 

## crop to area just around the Broads priority landscape
BroadsCoast <- st_transform(UKCoast, crs = st_crs(Broads)) |> st_crop(Broads |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Broads)[1]-2000, ext(Broads)[2]+2000), ylim = c(ext(Broads)[3]-2000, ext(Broads)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## Create plot of standing water cover
SomWa <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=BroadsCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  
  ## Add the landscape boundary
  geom_sf(data=Broads, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the water raster
  geom_spatraster(data = as.factor(BroWater)) +
  scale_fill_viridis_d(name = "Pixel Classification:", option = "D", na.value = NA, 
                       labels = c("Dry grassland", "Standing water", "")) +
  guides(fill = guide_legend(na.translate = FALSE)) +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2 +
  theme(legend.position = "top")


## Save the plot
ggsave("CleanData/Paper Plots/Broads_StandingWater_Predicted.png", 
       plot = SomWa,  width = 20, height = 20, units = "cm")



## |- Essex ----

## Read in the raster of North Kent
EssWater <- rast("RawData/BinaryWaterRasters/Essex_BinaryWaterAll.tif")

## Filter out the landscapes of interest
Essex <- MyBoxes |> filter(Att3Value == "Essex") 

## crop to area just around the Essex priority landscape
EssexCoast <- st_transform(UKCoast, crs = st_crs(Essex)) |> st_crop(Essex |> st_buffer(dist = 2000))

## Set the plot extent so that all plots have the same area no matter if they have different land parcels
PlotExt <- coord_sf(xlim = c(ext(Essex)[1]-2000, ext(Essex)[2]+2000), ylim = c(ext(Essex)[3]-2000, ext(Essex)[4]+2000), 
                    crs = 27700, expand = FALSE) 

## Create plot of standing water cover
EssWa <- ggplot() +
  
  ## Add the coastline
  geom_sf(data=EssexCoast, mapping=aes(geometry=geometry), fill = "lightgrey", colour = NA) +
  
  ## Add the landscape boundary
  geom_sf(data=Essex, mapping=aes(geometry=geometry, colour = "#273746"), fill = "NA", ) + 
  scale_color_manual(values = c("#273746"), labels = c(""), name="Landscape\nBoundary") +
  
  ## Add the water raster
  geom_spatraster(data = as.factor(EssWater)) +
  scale_fill_viridis_d(name = "Pixel Classification:", option = "D", na.value = NA, 
                       labels = c("Dry grassland", "Standing water", "")) +
  guides(fill = guide_legend(na.translate = FALSE)) +
  
  ## Set plot extent so all plots have the same extent
  PlotExt +
  
  ## Add x and y axis labels
  ylab("Latitude") +
  xlab("Longitude") +
  
  ## Add scale bar
  annotation_scale(location = "br", line_width = unit(0.25, "cm"), height = unit(0.1, "cm"), pad_y = unit(0.1, "in")) +
  
  ## set theme
  theme_light() + 
  GeneralThemeing2 +
  theme(legend.position = "top")


## Save the plot
ggsave("CleanData/Paper Plots/Essex_StandingWater_Predicted.png", 
       plot = EssWa,  width = 21, height = 18, units = "cm")

