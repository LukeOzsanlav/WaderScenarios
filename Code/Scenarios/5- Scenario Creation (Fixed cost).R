##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## 23/04/2024
## 
## 
## Aim: Calculate base-level values for canvas i.e. Wader abundance and costs
## 
## 
##------------------------------------------------------##

## Load in packages & helper functions
pacman::p_load(tidyverse, # used
               data.table, #??
               sf, # used
               terra, # used
               randomForestSRC, # used
               exactextractr, # used
               tictoc, # used for timing scenario creation function
               ggpubr, # used
               truncnorm) # used

options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R")
source("Code/Scenarios/5.2- Scenario Functions (Fixed cost).R")


## Overview of Random Scenario modelling
## - **Read in annotated canvas
## - **Apply the final masks to each field
## - **Predict wader abundance across all lowland grass fields
## - **Calculate the total area for each strategy (AES only on grass, Reserve only on grass, Reserve only on arable)
## - **Divide this area by 3 for bigger/more strategies
## - **Calculate 20/40/60/80/100% of area for each strategy to reach
## - **Randomly draw fields from opportunity till percentage of area surpassed (some kind of sample while function?)
## - **Assign these fields habitat values from corresponding strategy (E.g. from AES only or reserve fields)
## - **Calculate all fields that fall within 2km of fields that have changed
## - **Recalculate variable landscape variables for each field in buffer
## - **Predict the number of waders in the whole canvas and save this
## - **Save the important variables of the canvas so that I could plot changes across 20/40/60/80/100%



##-------------------------------##
## F 1 START BatchBake Function  ##
##-------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/"
# SaveCanv = FALSE
# N_sets = 2
# runsetting = Scen_settingsNK
# zstart = 1
# CostSpread = 20
# ClustBuf=200


## This function take the CreatScenario function and batch runs them over a spreadsheet with scenario settings
BatchBake <- function(Canvas, 
                      SniModel, 
                      LapModel,  
                      RedModel,
                      LandCov, 
                      Outpath, 
                      N_sets, 
                      SaveCanv, 
                      runsetting, 
                      CostSpread,
                      ClustBuf,
                      zstart){
  
  ## read in average sizes for AES and resevres
  clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
  
  ## Loop through each row in the setting doc as these will all be separate scenarios
  for(z in zstart:nrow(runsetting)){
    
    ## Update to console on progress
    message("Running scenario ", z, " out of ", nrow(runsetting))
    

    ##Join together opportunity categories for CreateScenario() function
    OppCatChoice = c(paste0(if(is.na(runsetting$OppCat1[z])==F){paste0(runsetting$OppCat1[z])}), 
                     paste0(if(is.na(runsetting$OppCat2[z])==F){paste0(runsetting$OppCat2[z])}))
    
    ## set the new category that select parcels become                 
    NewCatChoice = paste0(runsetting$NewCat1[z])
    
    ## Set what the minimum cluster size to aim for should be
    ClustMin <- ifelse(NewCatChoice =="Reserve", clustsizes$ReserveAve, clustsizes$AESAve)
    
    # Create directory if needs be
    dir.create(Outpath, showWarnings = F)
    
    ## Run all of the scenarios
    CreateScenario(Canvas=Canvas, SniModel=SniModel, LapModel=LapModel, RedModel=RedModel, LandCov=LandCov,
                   N_sets=N_sets, Budget = runsetting$Budget[z], ClustBuf= ClustBuf, 
                   OppCat= OppCatChoice, NewCat= NewCatChoice, Strategy = runsetting$Strategy[z],
                   CostSpread = CostSpread, SaveCanv = F, Outpath = Outpath, MinClust = ClustMin)
  
  }
} 
## F 1 END ##






##----------------------------##
##----------------------------##
#### 1 *Somerset Scenarios* ####
##----------------------------##
##----------------------------##

##----------------------##
#### 1.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Som_AnnotatedCanv.rds")
table(duplicated(Canv$ParcRef))

## Read in csv of the wader densities
WaderDens <- read.csv("CleanData/Scenarios/4-AnnotateCanvas/ReserveWaderDensity.csv") |> select(ReserveGroup, ResQual)


## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## Join the reserve quality
Canv <- left_join(Canv, WaderDens, by = "ReserveGroup") |> st_as_sf()


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

## Add a new column that separates out breeding wader vs wintering wader payment
Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})

## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#            theme_minimal()
            
## Get the areas of different categories
Canv |> st_drop_geometry() |> group_by(Category) |> summarise(TotArea = sum(FieldArea))


## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_Som <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_Som) # free up space

## Read in wader models
SniMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/SnipeRF_FullData.rds")




##---------------------------------##
#### 1.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
# TEST <- Canv[PrIndex, ] |> select(SniMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()

## Predict abundance for Snipe
set.seed(1212)
Canv$SnipeAbund[PrIndex] <- (predict.rfsrc(object = SniMod, 
                             newdata = (Canv[PrIndex, ] |> select(SniMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$SnipeAbund, na.rm = T) # number of Snipe in the landscape




##-----------------------------##
#### 1.3 Calc starting costs ####
##-----------------------------##

##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
WhichAES <- which(Canv$Category == "AES Only")

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){
  
  ## How many years to spread one off costs over
  Yrs <- 20

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
sum(Canv$AESCostGDP[WhichAES], na.rm=T)/1000



##-- Calculate the costs of the reserve management in the canvas --##

## Initiate a column for creating the costings
Canv$ResMaintGDP <- 0
Canv$ResCreatGDP <- 0

## Work out which rows are reserves
WhichReserve <- which(Canv$Category == "Reserve")

## Work out the costs of ongoing management of the field per hectare
Canv$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canv$GroupArea[WhichReserve],
                                               ParcelArea = Canv$FieldArea[WhichReserve],
                                               InflationAdjust=T)

## Calculate the costs of reserve management in £1000's of pounds
sum(Canv$ResMaintGDP[WhichReserve], na.rm=T)/1000




##------------------------------------##
#### 1.4 Calc Opportunity Area Size ####
##------------------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")
lowset <- Scen_settings |> filter(Budget == "low")

## run the function to calculate the size of the opportunity areas for each scenario
SomAreas <- AreaForScenarios(Canvas=Canv, ScenList=lowset)




##-----------------------------##
#### 1.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")

## Read in spreadsheet with AES costs for each region
GovSpend <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Total_AES_Expenditure.csv")
AveSpend <- mean(GovSpend$Cost)
Scen_settingsSom <- Scen_settings |> mutate(Budget= case_when(Budget=="low" ~ AveSpend*0.1,
                                          Budget=="medium" ~ AveSpend*0.2,
                                          Budget=="high" ~ AveSpend*0.3,
                                          .default = 0))
## Run the batch function 
BatchBake(Canvas=Canv, SniModel=SniMod, LapModel=NULL,  RedModel=NULL, LandCov=LC_Som,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Som/SetCost/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsSom, zstart = 1, CostSpread = 20, ClustBuf=200)




##---------------------------------##
#### 1.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Som/SetCost/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Som/SetCost/Plots/",
             species="Snipe")





# ##------------------------------##
# #### 1.7 Plot Scenario Canvas ####
# ##------------------------------##
# 
# ## Create a function that will read in the canvas's for a scenario and plot them
# PlotProg <- function(infolder, outfolder, scenario){
#   
#   ## list the different scenarios run in folder
#   ## just extract the first canvas as I will also need to read in the third and fifth canvas
#   paths = dir(infolder, pattern     = "1.shp$", full.names  = T)
#   
#   # Remove the end of the file paths so I can use it do read in the first third and fifth canvas
#   updated_filepaths = lapply(paths, function(filepath){substr(filepath, 1, nchar(filepath) - 5)})
# 
#   ## Get the file names instead of the full file paths this time
#   filenames = dir(infolder, pattern     = "1.shp$", full.names  = F)
#   ## remove the fist part of the file path
#   updated_filenames = lapply(filenames, function(filepath){substr(filepath, 8, nchar(filepath) - 5)})
#   
#   ## Loop through each canvas type
#   for(t in 1:length(updated_filepaths)){
#     
#   ## send message to console
#   message("Canvas ", t, " out of ", length(updated_filepaths))  
#     
#   ## Create the start of a title for each plot
#   title = gsub("_", " ", updated_filenames[t])
#   
#   
#   ## Plot canvas after first round
#   Start <- Canv|> ggplot() + 
#     geom_sf(mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#     ggtitle(paste0(title, ": Start")) +
#     scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "darkgrey", "lightgrey", "#4C77DD")) +
#     theme_light()
#   
#   ## Plot canvas after first round
#   P1 <- st_read(paste0(updated_filepaths[t], "1.shp")) |> ggplot() + 
#     geom_sf(mapping = aes(geometry = geometry, fill = Categry), colour = NA) +
#     ggtitle(paste0(title, ": Round 1")) +
#     scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "darkgrey", "lightgrey", "#4C77DD")) +
#     theme_light()
#   
#   ## Plot canvas after third round
#   P3 <- st_read(paste0(updated_filepaths[t], "3.shp")) |> ggplot() + 
#     geom_sf(mapping = aes(geometry = geometry, fill = Categry), colour = NA) +
#     ggtitle(paste0(title, ": Round 3")) +
#     scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "darkgrey", "lightgrey", "#4C77DD")) +
#     theme_light()
#   
#   ## Plot canvas after final round
#   P5 <- st_read(paste0(updated_filepaths[t], "5.shp")) |> ggplot() + 
#     geom_sf(mapping = aes(geometry = geometry, fill = Categry), colour = NA) +
#     ggtitle(paste0(title, ": Round 5")) +
#     scale_fill_manual(values = c("#DD4CC0", "#DDB24C", "#4CDD6A", "darkgrey", "lightgrey", "#4C77DD")) +
#     theme_light()
#   
#   
#   ## arrange the three plots in a row
#   FullPlot <- ggarrange(Start, P1, P3, P5, nrow = 2, ncol =2, common.legend = TRUE, legend = "bottom")
#   
#   ## save the plot as a png
#   plotout = paste0(outfolder, scenario, "_", updated_filenames[t], ".png")
#   ggsave(plot = FullPlot, filename = plotout, width = 26, height = 26, units = "cm")
#     
#     
#   }
#   
# } 
# ## END ##
# 
# 
# ## Run the function for all the canvas boards that I created
# PlotProg(infolder = "CleanData/Scenarios/5-ScenarioCreation/Som/random/",
#          outfolder = "CleanData/Scenarios/5-ScenarioCreation/Som/Plots/",
#          scenario = "random")
# 
# PlotProg(infolder = "CleanData/Scenarios/5-ScenarioCreation/Som/clusterlarge/",
#          outfolder = "CleanData/Scenarios/5-ScenarioCreation/Som/Plots/",
#          scenario = "cluster-large")
# 
# PlotProg(infolder = "CleanData/Scenarios/5-ScenarioCreation/Som/clustersmall/",
#          outfolder = "CleanData/Scenarios/5-ScenarioCreation/Som/Plots/",
#          scenario = "cluster-small")
# 
# PlotProg(infolder = "CleanData/Scenarios/5-ScenarioCreation/Som/stakeholder2G1/",
#          outfolder = "CleanData/Scenarios/5-ScenarioCreation/Som/Plots/",
#          scenario = "stakeholder-G1")







##-------------------------------##
##-------------------------------##
#### 2. *North Kent Scenarios* ####
##-------------------------------##
##-------------------------------##

##----------------------##
#### 2.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.shp") |> select(ParcRef)
table(duplicated(Canvshp$ParcRef))

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/NKent_AnnotatedCanv.rds")
table(duplicated(Canv$ParcRef))

## Read in csv of the wader densities
WaderDens <- read.csv("CleanData/Scenarios/4-AnnotateCanvas/ReserveWaderDensity.csv") |> select(ReserveGroup, ResQual)

## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## Join the the reserve quality
Canv <- left_join(Canv, WaderDens, by = "ReserveGroup") |> st_as_sf()

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

## Add a new column that separates out breeding wader vs wintering wader payment
Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})

## Visualize canvas categories
ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
           theme_minimal()

## Get the areas of different categories
Canv |> st_drop_geometry() |> group_by(Category) |> summarise(TotArea = sum(FieldArea))



## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds")





##---------------------------------##
#### 2.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------##
#### 2.3 Calc starting costs ####
##-----------------------------##

##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
WhichAES <- which(Canv$Category == "AES Only")

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){
  
  ## How many years to spread one off costs over
  Yrs <- 20

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
sum(Canv$AESCostGDP[WhichAES], na.rm=T)/1000



##-- Calculate the costs of the reserve management in the canvas --##

## Initiate a column for creating the costings
Canv$ResMaintGDP <- 0
Canv$ResCreatGDP <- 0

## Work out which rows are reserves
WhichReserve <- which(Canv$Category == "Reserve")

## Work out the costs of ongoing management of the field per hectare
Canv$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canv$GroupArea[WhichReserve],
                                               ParcelArea = Canv$FieldArea[WhichReserve],
                                               InflationAdjust=T)

## Calculate the costs of reserve management in £1000's of pounds
sum(Canv$ResMaintGDP[WhichReserve], na.rm=T)/1000




##------------------------------------##
#### 2.4 Calc Opportunity Area Size ####
##------------------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")
lowset <- Scen_settings |> filter(Budget == "low")

## run the function to calculate the size of the opportunity areas for each scenario
NKAreas <- AreaForScenarios(Canvas=Canv, ScenList=lowset)



##-----------------------------##
#### 2.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")

## Read in spreadsheet with AES costs for each region
GovSpend <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Total_AES_Expenditure.csv")
AveSpend <- mean(GovSpend$Cost)
Scen_settingsNK <- Scen_settings |> mutate(Budget= case_when(Budget=="low" ~ AveSpend*0.1,
                                          Budget=="medium" ~ AveSpend*0.2,
                                          Budget=="high" ~ AveSpend*0.3,
                                          .default = 0))
## Run the batch function 
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsNK, zstart = 1, CostSpread = 20, ClustBuf=200)





##---------------------------------##
#### 2.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/Plots/",
             species="LapRed")







##--------------------------------##
##--------------------------------##
#### 3. *Essex Coast Scenarios* ####
##--------------------------------##
##--------------------------------##

##----------------------##
#### 3.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Essex_AnnotatedCanv.rds")
table(duplicated(Canvshp$ParcRef))

## Read in csv of the wader densities
WaderDens <- read.csv("CleanData/Scenarios/4-AnnotateCanvas/ReserveWaderDensity.csv") |> select(ReserveGroup, ResQual)

## Read in a data set of Wallasea ISland F_LOC_IDs that are not wet grassland but actually saltmarsh
Wall <- read.csv("RawData/RSPB Reserves/WallaseaIsland_NonWetgrass.csv")


## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## Join the the reserve quality
Canv <- left_join(Canv, WaderDens, by = "ReserveGroup") |> st_as_sf()


## If any parcels are current reserve, AES-only or opportunity then apply some simple additional masks
## This will remove any fields that are not suitable for wader management
Canv <- mutate(Canv, 
               Category = ifelse(Mask_G1 > 0.5 & !Category == "NoOpp", "Masked Hab", Category), # apply a mask to canvas (priority habs, monuments...etc)
               Category = ifelse(F_LOC_ID %in% Wall$F_LOC_ID, "NoOpp", Category),
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

## Add a new column that separates out breeding wader vs wintering wader payment
Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})

## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#            theme_minimal()
            
## Get the areas of different categories
Canv |> st_drop_geometry() |> group_by(Category) |> summarise(TotArea = sum(FieldArea))


## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds")




##---------------------------------##
#### 3.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------##
#### 3.3 Calc starting costs ####
##-----------------------------##

##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
WhichAES <- which(Canv$Category == "AES Only")

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){
  
  ## How many years to spread one off costs over
  Yrs <- 20

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
sum(Canv$AESCostGDP[WhichAES], na.rm=T)/1000



##-- Calculate the costs of the reserve management in the canvas --##

## Initiate a column for creating the costings
Canv$ResMaintGDP <- 0
Canv$ResCreatGDP <- 0

## Work out which rows are reserves
WhichReserve <- which(Canv$Category == "Reserve")

## Work out the costs of ongoing management of the field per hectare
Canv$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canv$GroupArea[WhichReserve],
                                               ParcelArea = Canv$FieldArea[WhichReserve],
                                               InflationAdjust=T)

## Calculate the costs of reserve management in £1000's of pounds
sum(Canv$ResMaintGDP[WhichReserve], na.rm=T)/1000




##------------------------------------##
#### 3.4 Calc Opportunity Area Size ####
##------------------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")
lowset <- Scen_settings |> filter(Budget == "low")

## run the function to calculate the size of the opportunity areas for each scenario
EsAreas <- AreaForScenarios(Canvas=Canv, ScenList=lowset)




##-----------------------------##
#### 3.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")

## Read in spreadsheet with AES costs for each region
GovSpend <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Total_AES_Expenditure.csv")
AveSpend <- mean(GovSpend$Cost)
Scen_settingsEs <- Scen_settings |> mutate(Budget= case_when(Budget=="low" ~ AveSpend*0.1,
                                          Budget=="medium" ~ AveSpend*0.2,
                                          Budget=="high" ~ AveSpend*0.3,
                                          .default = 0))
## Run the batch function 
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Essex/SetCost/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsEs, zstart = 1, CostSpread = 20, ClustBuf=200)



##---------------------------------##
#### 3.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Essex/SetCost/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Essex/SetCost/Plots/",
             species="LapRed")





##-----------------------------------##
##-----------------------------------##
#### 4. *Norfolk Broads Scenarios* ####
##-----------------------------------##
##-----------------------------------##

##----------------------##
#### 4.1 Read in data ####
##----------------------##

## Read in canvas polygons
Canvshp <- st_read("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.shp") |> select(ParcRef)

## Read in annotated canvas as csv
Canv <- readRDS("CleanData/Scenarios/4-AnnotateCanvas/Broads_AnnotatedCanv.rds")

## Read in csv of the wader densities
WaderDens <- read.csv("CleanData/Scenarios/4-AnnotateCanvas/ReserveWaderDensity.csv") |> select(ReserveGroup, ResQual)


## Join polygons to annotations
Canv <- left_join(Canv, Canvshp, by = "ParcRef") |> st_as_sf()
rm(Canvshp)

## Join the the reserve quality
Canv <- left_join(Canv, WaderDens, by = "ReserveGroup") |> st_as_sf()


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

## Add a new column that separates out breeding wader vs wintering wader payment
Canv$BrAES <- lapply(Canv$AES_options, FUN=function(x){any(x %in% c("GS9", "GS11"))})

## Visualize canvas categories
# ggplot() + geom_sf(data = Canv, mapping = aes(geometry = geometry, fill = Category), colour = NA) +
#            theme_minimal()
            
## Get the areas of different categories
Canv |> st_drop_geometry() |> group_by(Category) |> summarise(TotArea = sum(FieldArea))

## Read in the UKCEH landcover data as a base map to rasterize my fields with
LC <- rast("RawData/LandCover/gblcm25m2021.tif")
LC <- LC[[1]]
LC_crop <- crop(LC, vect(Canv)) #|> mask(vect(Som))
plot(LC_crop) # free up space
rm(LC); gc()

## Read in wader models
LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")
RedMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/RedRF_FullData.rds")




##---------------------------------##
#### 4.2 Calc starting abundance ####
##---------------------------------##

## Get the index of the rows that I want to predict wader abundance into
PrIndex <- which(Canv$Category %in% c("Reserve", "AES Only", "Grass Opp"))
# TEST <- Canv[PrIndex, ] |> select(SniMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()

## Predict abundance for Snipe
set.seed(1212)
Canv$RedAbund[PrIndex] <- (predict.rfsrc(object = RedMod, 
                             newdata = (Canv[PrIndex, ] |> select(RedMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                             jitt=FALSE))$predicted
sum(Canv$RedAbund, na.rm = T) # number of Snipe in the landscape
  

## Predict abundance for Lapwing
set.seed(1212)
Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
                          newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
                          jitt=FALSE))$predicted
sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




##-----------------------------##
#### 3.3 Calc starting costs ####
##-----------------------------##

##-- Calculate the costs of the AES agreements in the canvas --##

## Calculate the costs of the AES agreements in the canvas

## Read in the AES scheme
## Remove the small field supplement as not going to use this option for now
AES_Costs <- read.csv("RawData/AES Costings/CSS_Cost_Sheet.csv")
SmallSup <-  AES_Costs |> filter(CSS_Code == "SP1")
AES_Costs <-  AES_Costs |> filter(!CSS_Code == "SP1")


## Initiate a column for creating the costings
Canv$AESCostGDP <- 0

## Work out which rows in the data set are AES fields
WhichAES <- which(Canv$Category == "AES Only")

## Now using the list column of the different AES schemes this function looks through all the AES codes
## Then ti works out the total payment for those fields based on the AES payment rates and the area of the field
for(j in 1:nrow(AES_Costs)){
  
  ## How many years to spread one off costs over
  Yrs <- 20

  ## create columns used to track AES payment
  Canv$Add <- NA
  Canv$AddAmount <- 0

  ## Work out which fields have the the current AES scheme or not (TRUE/FALSE)
  Canv$Add[WhichAES] <- lapply(Canv$AES_options[WhichAES], FUN=function(x){any(x %in% AES_Costs$CSS_Code[j])})

  ## Calculate the cost for that field of the current AES scheme
  ## Add this costs onto any existing costs from other schemes calculated earlier in the loop
  ## Also if the field is less then 1ha add on a small field supplement
  Canv[WhichAES,] <- Canv[WhichAES,] %>%
                     mutate(AddAmount = case_when(Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="Y" ~ FieldArea*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="m" & AES_Costs$Annual[j]=="Y" ~ as.numeric(Perim)*AES_Costs$Cost[j],
                                                  Add==T & AES_Costs$Unit[j]=="ha" & AES_Costs$Annual[j]=="N" ~ FieldArea*(AES_Costs$Cost[j]/Yrs),
                                                  Add==T & AES_Costs$Unit[j]=="item" & AES_Costs$Annual[j]=="N" ~ ((2*AES_Costs$Cost[j])/Yrs),
                                                  .default = 0),
                            AddAmount = ifelse(FieldArea<1, AddAmount+(FieldArea*SmallSup$Cost), AddAmount),
                            AESCostGDP = AESCostGDP + AddAmount)

  if(j == nrow(AES_Costs)){Canv <- Canv |> select(-c(Add, AddAmount))}
}

## Calculate the costs in £1000's of pounds
sum(Canv$AESCostGDP[WhichAES], na.rm=T)/1000



##-- Calculate the costs of the reserve management in the canvas --##

## Initiate a column for creating the costings
Canv$ResMaintGDP <- 0
Canv$ResCreatGDP <- 0

## Work out which rows are reserves
WhichReserve <- which(Canv$Category == "Reserve")

## Work out the costs of ongoing management of the field per hectare
Canv$ResMaintGDP[WhichReserve] <- Manage_Costs(TotalArea = Canv$GroupArea[WhichReserve],
                                               ParcelArea = Canv$FieldArea[WhichReserve],
                                               InflationAdjust=T)

## Calculate the costs of reserve management in £1000's of pounds
sum(Canv$ResMaintGDP[WhichReserve], na.rm=T)/1000




##------------------------------------##
#### 4.4 Calc Opportunity Area Size ####
##------------------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")
lowset <- Scen_settings |> filter(Budget == "low")

## run the function to calculate the size of the opportunity areas for each scenario
NBAreas <- AreaForScenarios(Canvas=Canv, ScenList=lowset)




##-----------------------------##
#### 4.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_fixedcost.csv")

## Read in spreadsheet with AES costs for each region
GovSpend <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Total_AES_Expenditure.csv")
AveSpend <- mean(GovSpend$Cost)
Scen_settingsBr <- Scen_settings |> mutate(Budget= case_when(Budget=="low" ~ AveSpend*0.1,
                                          Budget=="medium" ~ AveSpend*0.2,
                                          Budget=="high" ~ AveSpend*0.3,
                                          .default = 0))
## Run the batch function 
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Broads/SetCost/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsBr, zstart = 1, CostSpread = 20, ClustBuf=200)




##---------------------------------##
#### 4.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Broads/SetCost/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Broads/SetCost/Plots",
             species="LapRed")






##------------------------------##
##------------------------------##
#### 5. *Combine All Regions* ####
##------------------------------##
##------------------------------##


## Create empty list to put the output data for each region into
ListSet <- list()

## Create a vector of input paths for each region
Paths <- c("CleanData/Scenarios/5-ScenarioCreation/Broads/SetCost/",
           "CleanData/Scenarios/5-ScenarioCreation/Essex/SetCost/",
           "CleanData/Scenarios/5-ScenarioCreation/Kent/SetCost/",
           "CleanData/Scenarios/5-ScenarioCreation/Som/SetCost/")

## Create vector of the species that are used for each region
Species <- c("LapRed", "LapRed", "LapRed", "Snipe")


## Loop through each region and calculate a series of columns needed for plotting
for(j in 1:length(Paths)){ 
  
  ## number of years to spread one off costs over
  # Yr=15
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  files <- dir(paste0(Paths[j]), pattern  = "Track_", full.names  = T)
  
  
  ## Read in all these files and bind them together
  AllScn <- files %>%
              map(fread) %>% # read in all files into a list
              bind_rows() %>% # bind all the rows
              mutate(BaseWaderTot = rowSums(dplyr::select(., BaseSnipe, BaseLapwing, BaseRedshank), na.rm = T),
                     ChangeSnipe = Snipe-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeCosts = ifelse(NewCat=="AES Only", NewAESCost,
                                          ifelse(NewCat=="Reserve", NewResOverallCost, 0)),
                     PlotCat = paste0(NewCat, " ", ifelse(OppCat== "Arable Opp", "Arable", ""))) |>
              filter(!SegmentArea == 0) 

  
  ## Calculate the cost per pairs depending on whether it is Snipe or Lapwing/Redshank
  if(Species[j] == "LapRed"){ ListSet[[j]] <- AllScn |> mutate(PairCost = ChangeWaders/(ChangeCosts/Budget),
                                                                PairCost100k = ChangeWaders/(ChangeCosts/100000),
                                                                Waders_100Ha = (ChangeWaders/SegmentArea)*100,
                                                                Perc_Waders_100Ha = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(SegmentArea))*100,
                                                                Perc_PairCost = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(ChangeCosts)),
                                                                Perc_PairCostSt = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(ChangeCosts/Budget))) }
  
  if(Species[j] == "Snipe"){ ListSet[[j]] <- AllScn |> mutate(ChangeWaders = ChangeSnipe,
                                                               PairCost = ChangeSnipe/(ChangeCosts/Budget),
                                                               PairCost100k = ChangeWaders/(ChangeCosts/100000),
                                                               Waders_100Ha = (ChangeWaders/SegmentArea)*100,
                                                               Perc_Waders_100Ha = (((ChangeSnipe/(BaseSnipe))*100)/(SegmentArea))*100,
                                                               Perc_PairCost = (((ChangeSnipe/(BaseSnipe))*100)/(ChangeCosts)),
                                                               Perc_PairCostSt = (((ChangeSnipe/(BaseSnipe))*100)/(ChangeCosts/Budget))) }

}



## Now finally combine all the scenario outputs across regions
Set <- do.call("rbind", ListSet)

## set the start of the outpath for all of my plots
outpath <- "CleanData/Scenarios/5-ScenarioCreation/Combined/SetCost/"


## Set a general theme for the bar plots
## Set a general theme for the bar plots
BarPlotTheme <- theme_light() +
                theme(panel.grid.minor = element_blank(),
                      strip.text = element_text(size = 17, face = "bold"),
                      axis.title = element_text(size = 18),
                      axis.text.y = element_text(size = 16),
                      axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5),
                      legend.title = element_text(size = 17, face = "bold"),
                      legend.text = element_text(size = 16),
                      legend.position = "top")



##-------------------------------##
##  Pairs per 100 ha vs Strategy ##
##-------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSum <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""))) |>
  group_by(PlotCatFull, NewCat, Strategy) |>
  summarise(StEr = sd(Waders_100Ha)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            MeanWaders_100Ha = mean(Waders_100Ha)) |> 
  ungroup() |> 
  mutate(LowerCI = MeanWaders_100Ha - (t_score * StEr),
         UpperCI = MeanWaders_100Ha + (t_score * StEr))


## Create Bar plot
ggplot(ScenSum, aes(x= PlotCatFull, y= MeanWaders_100Ha)) +
  geom_col(aes(fill = Strategy), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(.9), colour = "#424949") +
  scale_fill_manual(name = "Strategy",   # Change legend title
                      labels = c("Better"= "Better", "More"= "More", "Big" ="Bigger"),  # Change legend labels
                      values = c("Better"="#A5F076", "More"="#76A5F0", "Big"="#F076A5")) +
  ylab("Change in Breeding Wader Pairs/ 100 ha") +
  xlab("Scenario Category") +
  labs(fill = "Targeting Strategy") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Pairha_vs_AmalgBudget.png"), units = "in", height = 9, width = 11)




##------------------------------------##
##  Pairs per budget cost vs Strategy ##
##------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""))) |>
  group_by(PlotCatFull, NewCat, Strategy) |>
  summarise(StEr = sd(PairCost)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(PairCost)) |> 
  ungroup() |> 
  mutate(LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr))


## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost)) +
  geom_col(aes(fill = Strategy), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(.9), colour = "#424949") +
  scale_fill_manual(name = "Strategy",   # Change legend title
                      labels = c("Better"= "Better", "More"= "More", "Big" ="Bigger"),  # Change legend labels
                      values = c("Better"="#A5F076", "More"="#76A5F0", "Big"="#F076A5")) +
  ylab("Breeding Wader Pairs/Budget Cost") +
  xlab("Scenario Category") +
  labs(fill = "Targeting Strategy") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "PairBudget_vs_AmalgBudget.png"), units = "in", height = 9, width = 12)



 
##-------------------------------------##
## Pairs per 100 ha vs Strategy/Budget ##
##-------------------------------------##


## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSum <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(Waders_100Ha)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            MeanWaders_100Ha = mean(Waders_100Ha)) |> 
  ungroup() |> 
  mutate(LowerCI = MeanWaders_100Ha - (t_score * StEr),
         UpperCI = MeanWaders_100Ha + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

## Create Bar plot
ggplot(ScenSum, aes(x= PlotCatFull, y= MeanWaders_100Ha, group = BudgetW)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Change in Breeding Wader Pairs/ 100 ha") +
  xlab("Scenario Category") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Pairha_vs_SplitBudget.png"), units = "in", height = 9, width = 11)




##-------------------------------------------##
##  Pairs per budget cost vs Strategy/Budget ##
##-------------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(PairCost)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(PairCost)) |> 
  ungroup() |> 
  mutate(LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 2.5, 5.5),  # approximate starting positions
  xmax = c(2.5, 5.5, 8.5)   # approximate ending positions
)

## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost, group = BudgetW)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.2, 0.5)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Change in Breeding Wader Pairs / Budget used") +
  xlab("Mechanism") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)) +
  # Add labels
  annotate("text", x = 1.5, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 15,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 4, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 15,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 7, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 15,
           label = "More", fontface = "bold", size = 5) +
  scale_x_discrete(labels = c(
    "Better for AES Only"  = "AES only",
    "Better for Reserve"  = "Reserve\n(grassland conversion)",
    "Big for AES Only" = "AES only",
    "Big for Reserve" = "Reserve\n(grassland conversion)",
    "Big for Reserve\nfrom Arable" = "Reserve\n(arable reversion)",
    "More for AES Only"  = "AES only",
    "More for Reserve"  = "Reserve\n(grassland conversion)",
    "More for Reserve\nfrom Arable" = "Reserve\n(arable reversion)"))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "PairBudget_vs_SplitBudget.png"), units = "in", height = 9, width = 11)





##-------------------------------------##
##  Pairs per £100k vs Strategy/Budget ##
##-------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                         ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(PairCost100k)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCostK = mean(PairCost100k)) |> 
  ungroup() |> 
  mutate(LowerCI = Ave_PairCostK - (t_score * StEr),
         UpperCI = Ave_PairCostK + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCostK, group = BudgetW)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Change in Breeding Wader Pairs/£100,000") +
  xlab("Scenario Category") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 0.5))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Pair£100k_vs_SplitBudget.png"), units = "in", height = 9, width = 11)




##----------------------------------------##
## Total Scenario Area vs Strategy/Budget ##
##----------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSize <- Set |>
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                       ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                              paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
summarise(StEr = sd(SegmentArea)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          AveSegmentArea = mean(SegmentArea)) |> 
ungroup() |> 
mutate(LowerCI = AveSegmentArea - (t_score * StEr),
       UpperCI = AveSegmentArea + (t_score * StEr),
       BudgetW = fct_reorder(BudgetW, Budget))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 2.5, 5.5),  # approximate starting positions
  xmax = c(2.5, 5.5, 8.5)   # approximate ending positions
)


## Create Bar plot
ggplot(ScenSize, aes(x= PlotCatFull, y= AveSegmentArea, group = BudgetW)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.2, 0.5)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Area altered in scenario / ha") +
  xlab("Mechanism") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)) +
  # Add labels
  annotate("text", x = 1.5, y = max(ScenSize$UpperCI, na.rm = TRUE) + 40,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 4, y = max(ScenSize$UpperCI, na.rm = TRUE) + 40,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 7, y = max(ScenSize$UpperCI, na.rm = TRUE) + 40,
           label = "More", fontface = "bold", size = 5) +
  scale_x_discrete(labels = c(
    "Better for AES Only"  = "AES only",
    "Better for Reserve"  = "Reserve\n(grassland conversion)",
    "Big for AES Only" = "AES only",
    "Big for Reserve" = "Reserve\n(grassland conversion)",
    "Big for Reserve\nfrom Arable" = "Reserve\n(arable reversion)",
    "More for AES Only"  = "AES only",
    "More for Reserve"  = "Reserve\n(grassland conversion)",
    "More for Reserve\nfrom Arable" = "Reserve\n(arable reversion)"))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "TotalArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)



##----------------------------------------##
## Total Cluster Area vs Strategy/Budget ##
##----------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ClustSize <- Set |>
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                       ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                              paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
summarise(StEr = sd(AveClustArea)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          AveAveClustArea = mean(AveClustArea)) |> 
ungroup() |> 
mutate(LowerCI = AveAveClustArea - (t_score * StEr),
       UpperCI = AveAveClustArea + (t_score * StEr),
       BudgetW = fct_reorder(BudgetW, Budget))

# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 2.5, 5.5),  # approximate starting positions
  xmax = c(2.5, 5.5, 8.5)   # approximate ending positions
)

  
## Create Bar plot
ggplot(ClustSize, aes(x= PlotCatFull, y= AveAveClustArea, group = BudgetW)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.2, 0.5)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Cluster area / ha") +
  xlab("Mechanism") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)) +
  # Add labels
  annotate("text", x = 1.5, y = max(ClustSize$UpperCI, na.rm = TRUE) + 10,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 4, y = max(ClustSize$UpperCI, na.rm = TRUE) + 10,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 7, y = max(ClustSize$UpperCI, na.rm = TRUE) + 10,
           label = "More", fontface = "bold", size = 5) +
  scale_x_discrete(labels = c(
    "Better for AES Only"  = "AES only",
    "Better for Reserve"  = "Reserve\n(grassland conversion)",
    "Big for AES Only" = "AES only",
    "Big for Reserve" = "Reserve\n(grassland conversion)",
    "Big for Reserve\nfrom Arable" = "Reserve\n(arable reversion)",
    "More for AES Only"  = "AES only",
    "More for Reserve"  = "Reserve\n(grassland conversion)",
    "More for Reserve\nfrom Arable" = "Reserve\n(arable reversion)"))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "ClusterArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)



##----------------------------------------##
## Wider Cluster Area vs Strategy/Budget ##
##----------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
WClustSize <- Set |>
filter(NewCat == "Reserve") |> 
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                       ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                              paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
summarise(StEr = sd(AveClustWider)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          AveAveClustWider = mean(AveClustWider)) |> 
ungroup() |> 
mutate(LowerCI = AveAveClustWider - (t_score * StEr),
       UpperCI = AveAveClustWider + (t_score * StEr),
       BudgetW = fct_reorder(BudgetW, Budget))

  
## Create Bar plot
ggplot(WClustSize, aes(x= PlotCatFull, y= AveAveClustWider, group = BudgetW)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                      values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("Wider Cluster area/ha") +
  xlab("Scenario Category") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "ClusterArea_vs_SplitBudget.png"), units = "in", height = 9, width = 11)




##---------------------------------------------------##
##  Percentage Pairs per £100,000 vs Strategy/Budget ##
##---------------------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                          ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                 paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(Perc_PairCost)) |> 
  ungroup() |> 
  mutate(LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 2.5, 5.5),  # approximate starting positions
  xmax = c(2.5, 5.5, 8.5)   # approximate ending positions
)

## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost*100000, group = BudgetW)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.2, 0.5)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                    values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("% Change in Breeding Wader Pairs per £100,000") +
  xlab("Mechanism") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)) +
  # Add labels
  annotate("text", x = 1.5, y = max(ScenSumCost$UpperCI*100000, na.rm = TRUE) + 0.75,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 4, y = max(ScenSumCost$UpperCI*100000, na.rm = TRUE) + 0.75,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 7, y = max(ScenSumCost$UpperCI*100000, na.rm = TRUE) + 0.75,
           label = "More", fontface = "bold", size = 5) +
  scale_x_discrete(labels = c(
    "Better for AES Only"  = "AES only",
    "Better for Reserve"  = "Reserve\n(grassland conversion)",
    "Big for AES Only" = "AES only",
    "Big for Reserve" = "Reserve\n(grassland conversion)",
    "Big for Reserve\nfrom Arable" = "Reserve\n(arable reversion)",
    "More for AES Only"  = "AES only",
    "More for Reserve"  = "Reserve\n(grassland conversion)",
    "More for Reserve\nfrom Arable" = "Reserve\n(arable reversion)"))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "PercPair£100000_vs_SplitBudget.png"), units = "in", height = 9, width = 11)




##------------------------------------------------------##
##  Percentage Pairs per Budget Used vs Strategy/Budget ##
##------------------------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
  mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
         BudgetW = ifelse(Budget==max(Budget), paste0("High: £", round(max(Budget)/1000), "k"),
                          ifelse(Budget==min(Budget), paste0("Low: £", round(min(Budget)/1000), "k"),
                                 paste0("Medium: £", round(median(Budget)/1000), "k")))) |>
  group_by(PlotCatFull, NewCat, Strategy, Budget, BudgetW) |>
  summarise(StEr = sd(Perc_PairCostSt)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(Perc_PairCostSt)) |> 
  ungroup() |> 
  mutate(LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr),
         BudgetW = fct_reorder(BudgetW, Budget))

# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 2.5, 5.5),  # approximate starting positions
  xmax = c(2.5, 5.5, 8.5)   # approximate ending positions
)

## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Ave_PairCost, group = BudgetW)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.2, 0.5)) +
  geom_col(aes(fill = BudgetW), width=0.55, position=position_dodge(0.55)) +
  geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=.2, position=position_dodge(0.55), colour = "#424949") +
  scale_fill_manual(name = "Budget",
                    values = c("#A5F076", "#76A5F0", "#F076A5")) +
  ylab("% Change in Breeding Wader Pairs / Budget used") +
  xlab("Mechanism") +
  labs(fill = "Budget") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)) +
  # Add labels
  annotate("text", x = 1.5, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 5,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 4, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 5,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 7, y = max(ScenSumCost$UpperCI, na.rm = TRUE) + 5,
           label = "More", fontface = "bold", size = 5) +
  scale_x_discrete(labels = c(
    "Better for AES Only"  = "AES only",
    "Better for Reserve"  = "Reserve\n(grassland conversion)",
    "Big for AES Only" = "AES only",
    "Big for Reserve" = "Reserve\n(grassland conversion)",
    "Big for Reserve\nfrom Arable" = "Reserve\n(arable reversion)",
    "More for AES Only"  = "AES only",
    "More for Reserve"  = "Reserve\n(grassland conversion)",
    "More for Reserve\nfrom Arable" = "Reserve\n(arable reversion)"))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "PercPairBudget_vs_SplitBudget.png"), units = "in", height = 9, width = 11)
