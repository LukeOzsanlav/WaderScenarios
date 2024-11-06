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
source("Code/Scenarios/5.2- Scenario Functions.R")

## Set yrs variables
## This controls how many years fixed costs are spread out over
Yrs=20


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


## THINGS TO DO:
## - Maybe change data in the RF model so landscape standing water and low wet grass is calculated in new way. 






##-------------------------------------##
## F 1 START AreaForScenarios Function ##
##-------------------------------------##

## Given an annotate canvas for a priority landscape work out the total opportunity area for all scenarios
## START of AreaForScenarios function
AreaForScenarios <- function(Canvas, ## the canvas for the priority landscape, created in `4 - Annotate scenario canvas.R`
                             ScenList ## list of different scenario types that are going to by run
                             ){
  
  ## Create a data frame to put all different opportunity area sizes into
  ScenAreas <- data.frame(Landscape = Canvas$Landscape[1:nrow(ScenList)],
                          Strategy = NA,
                          OppCat= NA,
                          NewCat=  NA, 
                          Plus = NA, 
                          OppArea = NA)
  
  ## run through each row in the data frame of scenarios
  for(j in 1:nrow(ScenList)){
    
    ##Join together opportunity and replacement categories for CreateScenario() function
    OppCatChoice = c(paste0(if(is.na(ScenList$OppCat1[j])==F){paste0(ScenList$OppCat1[j])}),
                     paste0(if(is.na(ScenList$OppCat2[j])==F){paste0(ScenList$OppCat2[j])}),
                     paste0(if(is.na(ScenList$OppCat3[j])==F){paste0(ScenList$OppCat3[j])}),
                     paste0(if(is.na(ScenList$OppCat4[j])==F){paste0(ScenList$OppCat4[j])}))
                     
    NewCatChoice = c(paste0(if(is.na(ScenList$NewCat1[j])==F){paste0(ScenList$NewCat1[j])}), 
                     paste0(if(is.na(ScenList$NewCat2[j])==F){paste0(ScenList$NewCat2[j])}))
  
    Strategy = ScenList$Strategy[j]
    Plus = ScenList$Plus[j]
    
    ## Get the index of rows representing the opportunity area for this scenario
    ## This is decided based upon whether fields are inside or outside of wader clusters and the category of the field
    if(Strategy %in% c("Big", "More") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice) }
    if(Strategy %in% c("Better") & Plus ==F){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice) }
    
    ## For the plus strategies we also remove any fields that are the same category as the new management type and already have waders 
    ## (i.e. AES fields with waders will not be targeted and those without will be)
    if(Strategy %in% c("Big", "More") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==T & Canvas$Category %in% OppCatChoice & !(Canvas$Category %in% NewCatChoice & Canvas$Tot_abund > 0)) }
    if(Strategy %in% c("Better") & Plus ==T){ IndOpp <- which(is.na(Canvas$ClustGroup)==F & Canvas$Category %in% OppCatChoice & !(Canvas$Category %in% NewCatChoice & Canvas$Tot_abund > 0)) }
    
    ## Calculate the total field Area for my opportunity Area
    ScenAreas$OppArea[j] <- sum(Canvas[IndOpp,]$FieldArea)*ScenList$PropOppArea[j]
    
    ScenAreas$Strategy[j] <- Strategy
    ScenAreas$OppCat[j] <- paste0(paste0(if(is.na(ScenList$OppCat1[j])==F){paste0(ScenList$OppCat1[j], " ")}), 
                                  paste0(if(is.na(ScenList$OppCat2[j])==F){paste0(ScenList$OppCat2[j], " ")}), 
                                  paste0(if(is.na(ScenList$OppCat3[j])==F){paste0(ScenList$OppCat3[j], " ")}), 
                                  paste0(if(is.na(ScenList$OppCat4[j])==F){paste0(ScenList$OppCat4[j])}))
    ScenAreas$NewCat[j] <- paste0(paste0(if(is.na(ScenList$NewCat1[j])==F){paste0(ScenList$NewCat1[j], " ")}), 
                                  paste0(if(is.na(ScenList$NewCat2[j])==F){paste0(ScenList$NewCat2[j], " ")}))
    ScenAreas$Plus[j] <- Plus
    
  } ## End of j loop
  
  return(ScenAreas)
  
} 
## F 1 END ##









##-------------------------------##
## F 2 START BatchBake Function  ##
##-------------------------------##

# Canvas=Canv
# SniModel=NULL
# LapModel=LapMod
# RedModel=RedMod
# LandCov=LC_crop
# Outpath = "CleanData/Scenarios/5-ScenarioCreation/Broads/"
# SaveCanv = FALSE
# N_sets = 5
# runsetting = Scen_settingsBr
# zstart = 1
# Additive = TRUE


## This function take the CreatScenario function and batch runs them over a spreadsheet with scenario settings
BatchBake <- function(Canvas, SniModel, LapModel,  RedModel, LandCov, Outpath, 
                      N_sets, Additive, SaveCanv, runsetting, Unfenced=FALSE, zstart=1, zend=nrow(runsetting)){
  
  ##-----------------------------##        
  ## F 2.1 Calc Opp Area Targets ##
  ##-----------------------------##
  
  message("Preparing to run scenarios...")
  
  ## run `AreaForScenarios` function to calculate the size of the opportunity areas for each scenario
  ## only need to do this for random as it is the same for clustered appraoches
  AreasList <- AreaForScenarios(Canvas=Canvas, ScenList=(runsetting |> filter(ScenType == "random")))
  
  ## Calculate the minimum opportunity area for each strategy/lawton principle
  AreasSum <- AreasList |> 
               group_by(Strategy) |> 
               summarise(MinOppArea = min(OppArea, na.rm = T))
  
  ## Now join the MinOppArea onto the full list of scenarios
  runsetting <- left_join(runsetting, AreasSum, by = "Strategy") 
  
  ## When additive = TRUE I need a vector of areas to test that are the same length as N_sets
  ## I could just repeat the same area but I though it would be good to introduce a bit of variance
  ## For now I will use the cluster sizes from the large clusters
  if(Additive==TRUE){
      clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
      threshold_areaList <- rtruncnorm(n=N_sets, a=0, mean = clustsizes$ReserveAve, sd = max(runsetting$ClustSD, na.rm=T))
    }

  
  
  ##---------------------------##        
  ## F 2.2 Loop over Scenarios ##
  ##---------------------------##

  ## Loop through each row in the setting doc as these will all be separate scenarios
  for(z in zstart:zend){
    
    ## Update to console on progress
    message("Running scenario ", z, " out of ", zend)
    

    ##Join together opportunity and replacement categories for CreateScenario() function
    OppCatChoice = c(paste0(if(is.na(runsetting$OppCat1[z])==F){paste0(runsetting$OppCat1[z])}), 
                     paste0(if(is.na(runsetting$OppCat2[z])==F){paste0(runsetting$OppCat2[z])}),
                     paste0(if(is.na(runsetting$OppCat3[z])==F){paste0(runsetting$OppCat3[z])}),
                     paste0(if(is.na(runsetting$OppCat4[z])==F){paste0(runsetting$OppCat4[z])}))
                     
    NewCatChoice = c(paste0(if(is.na(runsetting$NewCat1[z])==F){paste0(runsetting$NewCat1[z])}), 
                     paste0(if(is.na(runsetting$NewCat2[z])==F){paste0(runsetting$NewCat2[z])}))
    
    ## create the file path for output
    ## Create a column earlier on to differentiate small vs large cluster
    Typ <- runsetting$ScenType[z]
    FinalOutPath = paste0(Outpath, if(Typ == "random"){paste0("rand", "/")}, 
                                   if(Typ == "cluster" & runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustlarge/")},
                                   if(Typ == "cluster" & !runsetting$ClustMean[z] == max(runsetting$ClustMean, na.rm =T)){paste0("clustsmall/")},
                                   if(Typ == "stakeholder2" | Typ == "stakeholder1"){paste0("stake", runsetting$StakeGroup[z], "/")})
    
    ## Define what will be the OppAreaUsed, this is either a single number of the total area across all N_sets
    ## Or a vecotr the same length of N_sets that dictate how much management is deployed on each sampling loop
    if(Additive==TRUE){OppAreaChoice <- threshold_areaList}
    if(Additive==FALSE){OppAreaChoice <- runsetting$MinOppArea[z]}
    
    # Create directory if needs be
    dir.create(FinalOutPath, showWarnings = F)
    
    ## Run all of the scenarios
    CreateScenario(Canvas=Canvas, SniModel=SniModel, LapModel=LapModel, RedModel=RedModel, LandCov=LandCov, N_sets=N_sets, Additive = Additive,
                   ScenType = runsetting$ScenType[z], Strategy = runsetting$Strategy[z], StakeGroup = runsetting$StakeGroup[z],
                   OppCat= OppCatChoice, NewCat= NewCatChoice, Plus = runsetting$Plus[z],  OppAreaUsed = OppAreaChoice,
                   ClustBuf= runsetting$ClustBuf[z], ClustMean= runsetting$ClustMean[z], ClustSD= runsetting$ClustSD[z],
                   Outpath = FinalOutPath, SaveCanv = SaveCanv, Unfenced=Unfenced)
  
  }
} 
## F 2 END ##






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
# LapMod <- readRDS("CleanData/Wader Abundance/5-ModelWaderAbundance/Model Objects/LapRF_FullData.rds")





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
  

## Predict abundance for Lapwing
# set.seed(1212)
# Canv$LapAbund[PrIndex] <- (predict.rfsrc(object = LapMod, 
#                           newdata = (Canv[PrIndex, ] |> select(LapMod[["xvar.names"]]) |> st_drop_geometry() |> as.data.frame()), 
#                           jitt=FALSE))$predicted
# sum(Canv$LapAbund, na.rm = T) # number of Lapwing in the landscape 




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
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")
Ran_scen <- Scen_settings |> filter(ScenType == "random")
Ran_scen$OppCat4 <- NA

Scen_settings2 <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_masterv2.csv")
Ran_scen2 <- Scen_settings2 |> filter(ScenType == "random")

## run the function to calculate the size of the opportunity areas for each scenario
SomAreas <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen)
SomAreas2 <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen2)




##-----------------------------##
#### 1.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")

## Read in spreadsheet with sizes for clusters
clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == min(Scen_settings$ClustMean, na.rm=T), clustsizes$AESAve, Scen_settings$ClustMean)
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == max(Scen_settings$ClustMean, na.rm=T), clustsizes$ReserveAve, Scen_settings$ClustMean)


## For Somerset I may want to remove the arable reversion scenarios
Scen_settingsSom <- filter(Scen_settings, !OppCat1 == "Arable Opp") #|> filter(!(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T) & Strategy == "Better"))
Scen_settingsSom$OppCat4 <- NA


## Run the batch function for additive
BatchBake(Canvas=Canv, SniModel=SniMod, LapModel=NULL,  RedModel=NULL, LandCov=LC_Som,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Som/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsSom, zstart = 1, Additive = TRUE, Unfenced=TRUE)

## Run the batch function for multiplicative
# BatchBake(Canvas=Canv, SniModel=SniMod, LapModel=NULL,  RedModel=NULL, LandCov=LC_Som,
#           Outpath = "CleanData/Scenarios/5-ScenarioCreation/Som/", SaveCanv = FALSE,
#           N_sets = 5, runsetting = Scen_settingsSom, zstart = 1, Additive = FALSE, Unfenced=TRUE)





##---------------------------------##
#### 1.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Som/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Som/Plots/",
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
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")
Ran_scen <- Scen_settings |> filter(ScenType == "random")
Ran_scen$OppCat4 <- NA

Scen_settings2 <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_masterv2.csv")
Ran_scen2 <- Scen_settings2 |> filter(ScenType == "random")

## run the function to calculate the size of the opportunity areas for each scenario
NKAreas <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen)
NKAreas2 <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen2)




##-----------------------------##
#### 2.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")

## Read in spreadsheet with sizes for clusters
clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == min(Scen_settings$ClustMean, na.rm=T), clustsizes$AESAve, Scen_settings$ClustMean)
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == max(Scen_settings$ClustMean, na.rm=T), clustsizes$ReserveAve, Scen_settings$ClustMean)

## For Broads I want to remove large cluster from Better scenarios
## For now also do not try and run any stakeholder scenarios as gradings not yet calculated
Scen_settingsNK<- Scen_settings #|> filter(!(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T) & Strategy == "Better"))
Scen_settingsNK$OppCat4 <- NA

## Run the batch function for additive
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsNK, zstart = 1, Additive = TRUE, Unfenced=TRUE)

## Run the batch function for multiplicative
# BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
#           Outpath = "CleanData/Scenarios/5-ScenarioCreation/Kent/", SaveCanv = FALSE,
#           N_sets = 5, runsetting = Scen_settingsNK, zstart = 1, Additive = FALSE, Unfenced=TRUE)




##---------------------------------##
#### 2.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Kent/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Kent/Plots/",
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
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")
Ran_scen <- Scen_settings |> filter(ScenType == "random")
Ran_scen$OppCat4 <- NA

Scen_settings2 <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_masterv2.csv")
Ran_scen2 <- Scen_settings2 |> filter(ScenType == "random")

## run the function to calculate the size of the opportunity areas for each scenario
EsAreas <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen)
EsAreas2 <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen2)




##-----------------------------##
#### 3.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")

## Read in spreadsheet with sizes for clusters
clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == min(Scen_settings$ClustMean, na.rm=T), clustsizes$AESAve, Scen_settings$ClustMean)
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == max(Scen_settings$ClustMean, na.rm=T), clustsizes$ReserveAve, Scen_settings$ClustMean)

## For Broads I want to remove large cluster from Better scenarios
## For now also do not try and run any stakeholder scenarios as gradings not yet calculated
Scen_settingsEs<- Scen_settings #|> filter(!(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T) & Strategy == "Better"))
Scen_settingsEs$OppCat4 <- NA


## Run the batch function for additive
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Essex/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsEs, zstart = 1, Additive = TRUE, Unfenced=TRUE)

## Run the batch function for multiplicative
# BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
#           Outpath = "CleanData/Scenarios/5-ScenarioCreation/Essex/", SaveCanv = FALSE,
#           N_sets = 5, runsetting = Scen_settingsEs, zstart = 1, Additive = FALSE, Unfenced=TRUE)





##---------------------------------##
#### 3.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Essex/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Essex/Plots/",
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
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")
Ran_scen <- Scen_settings |> filter(ScenType == "random")
Ran_scen$OppCat4 <- NA

Scen_settings2 <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_masterv2.csv")
Ran_scen2 <- Scen_settings2 |> filter(ScenType == "random")

## run the function to calculate the size of the opportunity areas for each scenario
NBAreas <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen)
NBAreas2 <- AreaForScenarios(Canvas=Canv, ScenList=Ran_scen2)




##-----------------------------##
#### 4.5 Batch Run Scenarios ####
##-----------------------------##

## Read in spreadsheet with setting for all scenarios I want to run
Scen_settings <- read_csv("CleanData/Scenarios/5-ScenarioCreation/All_scenarios_settings_master.csv")

## Read in spreadsheet with sizes for clusters
clustsizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == min(Scen_settings$ClustMean, na.rm=T), clustsizes$AESAve, Scen_settings$ClustMean)
Scen_settings$ClustMean <- ifelse(Scen_settings$ClustMean == max(Scen_settings$ClustMean, na.rm=T), clustsizes$ReserveAve, Scen_settings$ClustMean)

## For Broads I want to remove large cluster from Better scenarios
## For now also do not try and run any stakeholder scenarios as gradings not yet calculated
Scen_settingsBr <- Scen_settings #|> filter(!(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T) & Strategy == "Better"))
Scen_settingsBr$OppCat4 <- NA

## Run the batch function for additive
BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
          Outpath = "CleanData/Scenarios/5-ScenarioCreation/Broads/", SaveCanv = FALSE,
          N_sets = 20, runsetting = Scen_settingsBr, zstart = 1, Additive = TRUE, Unfenced=TRUE)

## Run the batch function for multiplicative
# BatchBake(Canvas=Canv, SniModel=NULL, LapModel=LapMod,  RedModel=RedMod, LandCov=LC_crop,
#           Outpath = "CleanData/Scenarios/5-ScenarioCreation/Broads/", SaveCanv = FALSE,
#           N_sets = 5, runsetting = Scen_settingsBr, zstart = 1, Additive = FALSE, Unfenced=TRUE)





##---------------------------------##
#### 4.6 Plot Scenario Abundance ####
##---------------------------------##

## Plot the outputs from scenario modelling
PlotScenario(inpath="CleanData/Scenarios/5-ScenarioCreation/Broads/",
             outpath="CleanData/Scenarios/5-ScenarioCreation/Broads/Plots/",
             species="LapRed")






##------------------------------##
##------------------------------##
#### 5. *Combine All Regions* ####
##------------------------------##
##------------------------------##

## Create empty list to put the output data for each region into
ListSet <- list()

## Create a vector of input paths for each region
Paths <- c("CleanData/Scenarios/5-ScenarioCreation/Broads/",
           "CleanData/Scenarios/5-ScenarioCreation/Essex/",
           "CleanData/Scenarios/5-ScenarioCreation/Kent/",
           "CleanData/Scenarios/5-ScenarioCreation/Som/")

## Create vector of the species that are used for each region
Species <- c("LapRed", "LapRed", "LapRed", "Snipe")

## file path for all of the plots made here
outpath <- "CleanData/Scenarios/5-ScenarioCreation/Combined/"


## Loop through each region and calculate a series of columns needed for plotting
for(j in 1:length(Paths)){ 
  
  ## First create a data set that can be used to create a bar plot that covers all scenarios
  files <- c(dir(paste0(Paths[j], "rand/"), pattern  = "Track_Add", full.names  = T),
             dir(paste0(Paths[j], "clustlarge/"), pattern  = "Track_Add", full.names  = T),
             dir(paste0(Paths[j], "clustsmall/"), pattern  = "Track_Add", full.names  = T))
  
  
  ## Read in all these files and bind them together
  AllScn2 <- files %>%
              map(fread) %>% # read in all files into a list
              bind_rows() %>% # bind all the rows
              mutate(ScenType = ifelse(ScenType == "cluster" & ClustMean == max(ClustMean, na.rm = T), "clusterlarge", ScenType),
                     ScenType = ifelse(ScenType == "cluster" & ClustMean == min(ClustMean, na.rm = T), "clustersmall", ScenType),
                     ChangeSnipe = Snipe-BaseSnipe,
                     ChangeSnipenoF = Snipe_unfenced-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeWadersnoF = (Lapwing_unfenced+Redshank_unfenced) - (BaseLapwing+BaseRedshank),
                     ChangeCostsFG = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yrs)+(FencingCost/15)+ForegoneCost)-BaseResMaintCost, 0)),
                     ChangeCostsPU = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yrs)+(FencingCost/15)+(PurchaseCost/Yrs))-BaseResMaintCost, 0)),
                     ChangeCostsNoF = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yrs)+ForegoneCost)-BaseResMaintCost, 0)),
                     PlotCat = paste0(NewCat, ifelse(Plus==T, "+", "-"), " ", ifelse(OppCat== "Arable Opp", "Arable", ""))) |>
              filter(!SegmentArea == 0) 
              
  ## ratio of income forgone to
  mean(AllScn2$PurchaseCost/AllScn2$ForegoneCost, na.rm=T)
  
  ## Calculate the cost per pairs depending on whether it is Snipe or Lapwing/Redshank
  ## Here I calculate the percentage change in breeding pairs per 100ha
  ## Then I calculate the percentage breeding pairs per pound (£), for different costings
  if(Species[j] == "LapRed"){ ListSet[[j]] <- AllScn2 |> mutate(Perc_Waders_100Ha = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(SegmentArea))*100,
                                                                 Perc_PairCost = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(ChangeCostsFG)),
                                                                 Perc_PairCostLap = (((ChangeLapwing/(BaseLapwing))*100)/(ChangeCostsFG)),
                                                                 Perc_PairCostRed = (((ChangeRedshank/(BaseRedshank))*100)/(ChangeCostsFG)),
                                                                 Perc_PairCostPU = (((ChangeWaders/(BaseLapwing+BaseRedshank))*100)/(ChangeCostsPU)),
                                                                 Perc_PairCostNoF = (((ChangeWadersnoF/(BaseLapwing+BaseRedshank))*100)/(ChangeCostsNoF))) }
  
            
  if(Species[j] == "Snipe"){ ListSet[[j]] <- AllScn2 |> mutate(Perc_Waders_100Ha = (((ChangeSnipe/(BaseSnipe))*100)/(SegmentArea))*100,
                                                                Perc_PairCost = (((ChangeSnipe/(BaseSnipe))*100)/(ChangeCostsFG)),
                                                                Perc_PairCostLap = NA,
                                                                Perc_PairCostRed = NA,
                                                                Perc_PairCostPU = (((ChangeSnipe/(BaseSnipe))*100)/(ChangeCostsPU)),
                                                                Perc_PairCostNoF = (((ChangeSnipenoF/(BaseSnipe))*100)/(ChangeCostsNoF))) }

}


## Now finally combine all the scenario outputs across regions
#ListSet[[1]] <- ListSet[[1]] |> mutate(ClustActual=NA) ## remove after I have rerun the scenarios
Set <- do.call("rbind", ListSet)


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


##-----------------------------------------##
##  Additive: Pairs per 100 ha vs Category ##
##-----------------------------------------##


## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSum <- Set |>
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_Waders_100Ha)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Perc_Waders_100Ha = mean(Perc_Waders_100Ha)) |>  
ungroup() |> 
mutate(LowerCI = Perc_Waders_100Ha - (t_score * StEr),
       UpperCI = Perc_Waders_100Ha + (t_score * StEr))


## Create Bar plot
ggplot(ScenSum, aes(x= PlotCatFull, y= Perc_Waders_100Ha, fill = ScenType)) +
geom_col(width=0.55, position=position_dodge(width=0.55)) +
geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=0.2, position=position_dodge(width=0.55), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Scenario Type",   # Change legend title
                    labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                    values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
ylab("Change in Breeding Wader Pairs/ 100 ha") +
xlab("Scenario Category") +
labs(fill = "Targeting Strategy") +
BarPlotTheme + 
theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pairha_vs_Category.png"), units = "in", height = 9, width = 11)




##-------------------------------------------##
##  Additive: Pairs per £100,000 vs Category ##
##-------------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSumCost <- Set |>
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Perc_PairCost = mean(Perc_PairCost)) |>     
ungroup() |> 
mutate(LowerCI = Perc_PairCost - (t_score * StEr),
       UpperCI = Perc_PairCost + (t_score * StEr))


## Create Bar plot
ggplot(ScenSumCost, aes(x= PlotCatFull, y= Perc_PairCost*100000, fill = ScenType)) +
geom_col(width=0.55, position=position_dodge(0.55)) +
geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=0.2, position=position_dodge(width=0.55), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Scenario Type",   # Change legend title
                    labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                    values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
ylab("Breeding Wader Pairs/ £100,000") +
xlab("Scenario Category") +
labs(fill = "Targeting Strategy") +
BarPlotTheme + 
theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))


## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_Category.png"), units = "in", height = 9, width = 12)




##------------------------------------------------------##
##  Additive: Pairs per 100 ha vs Summarized Categories ##
##------------------------------------------------------##

## Summarise the data so there are less rows
## Essentially the data is summarized to by taking the average across the different targeting strategies
ScenSum2 <- Set |>
        filter(!OppCat %in% c("Arable Opp")) |> 
        group_by(Strategy, Plus, NewCat) |>
        summarise(StEr = sd(Perc_Waders_100Ha)/sqrt(n()),
                  t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
                  Perc_Waders_100Ha = mean(Perc_Waders_100Ha)) |>
        ungroup() |>
        mutate(PlotCatOther = paste0(NewCat, " ", ifelse(Plus==T, "+", "-")),
               PlotCatFull= paste0(Strategy, " for ", PlotCatOther),
               PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality"),
               LowerCI = Perc_Waders_100Ha - (t_score * StEr),
               UpperCI = Perc_Waders_100Ha + (t_score * StEr)) 

## Create bar plot
## Much simpler than above and mainly focusses on showing the differences between AES and reserve for the different Lawton Principles
ggplot(ScenSum2, aes(x= Strategy, y= Perc_Waders_100Ha, fill = NewCat)) +
geom_col(width=0.55, position=position_dodge(0.55)) +
geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=0.2, position=position_dodge(width=0.55), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Scenario Type",   # Change legend title
                  labels = c("AES Only"= "AES Only", "Reserve"= "Nature Reserve"),  # Change legend labels
                  values = c("AES Only"="#F076A5", "Reserve"= "#76A5F0")) +
ylab("Change in Breeding Wader Pairs/ 100 ha") +
xlab("Lawton Principle") +
labs(fill = "Management Used") +
BarPlotTheme

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pairha_vs_SumCategory.png"), units = "in", height = 9, width = 11)




##--------------------------------------------------------##
##  Additive: Pairs per £100,000 vs Summarized Categories ##
##--------------------------------------------------------##

## Summarise the data so there are less rows
## Essentially the data is summarized to by taking the average across the different targeting strategies
ScenSumCost2 <- Set |>
          filter(!OppCat %in% c("Arable Opp")) |> 
        group_by(Strategy, Plus, NewCat) |>
        summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
                  t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
                  Perc_PairCost = mean(Perc_PairCost)) |>
        mutate(PlotCatOther = paste0(NewCat, " ", ifelse(Plus==T, "+", "-")),
               PlotCatFull= paste0(Strategy, " for ", PlotCatOther),
               PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality"),
               LowerCI = Perc_PairCost - (t_score * StEr),
               UpperCI = Perc_PairCost + (t_score * StEr))

## Create bar plot
## Much simpler than above and mainly focuses on showing the differences between AES and reserve for the different Lawton Principles
ggplot(ScenSumCost2, aes(x= Strategy, y= Perc_PairCost*100000, fill = NewCat)) +
geom_col(width=0.55, position=position_dodge(0.55)) +
geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=0.2, position=position_dodge(width=0.55), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Scenario Type",   # Change legend title
                  labels = c("AES Only"= "AES Only", "Reserve"= "Nature Reserve"),  # Change legend labels
                  values = c("AES Only"="#F076A5", "Reserve"= "#76A5F0")) +
ylab("Breeding Wader Pairs/ £100,000") +
xlab("Lawton Principle") +
labs(fill = "Management Used") +
BarPlotTheme

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_SumCategory.png"), units = "in", height = 9, width = 11)




##--------------------------------------##
##--------------------------------------##
##  Plot Scenarios with/without Fencing ##
##--------------------------------------##
##--------------------------------------##

## First create a data set that can be used to create a bar plot that covers all scenarios
## Calculate the change in pairs per unit cost for un-fenced reserve
ScenSumCost_NoFence <- Set |>
filter(!OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCostNoF)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCostNoF)) |> 
ungroup() |> 
mutate(Fencing = "Reserve Unfenced\n(grassland)",
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")

## Calculate the change in pairs per unit cost for fenced reserve
ScenSumCost_Fence <- Set |>
filter(!OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Fencing = "Reserve Fenced\n(grassland)",
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")

## Calculate the change in pairs per unit cost for AES only land
ScenSumCost_AES <- Set |>
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Fencing = "AES Only",
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "AES Only")



## Calculate the change in pairs per unit cost for un-fenced reserve from arable land
ScenSumCost_NoFArable <- Set |>
filter(OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCostNoF)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCostNoF)) |> 
ungroup() |> 
mutate(Fencing = "Reserve Unfenced\n(arable)",
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")

## Calculate the change in pairs per unit cost for fenced reserve from arable land
ScenSumCost_FenceArable <- Set |>
filter(OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Fencing = "Reserve Fenced\n(arable)",
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")


## Combine all data sets for plotting
FenceComp <- rbind(ScenSumCost_NoFence, ScenSumCost_Fence, ScenSumCost_AES, ScenSumCost_NoFArable, ScenSumCost_FenceArable)


## Create Bar plot
ggplot(FenceComp, aes(x= Strategy, y= Ave_PairCost*100000, fill = Fencing)) +
geom_col(width=0.55, position=position_dodge(0.6, preserve = "single")) +
geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=0.2, position=position_dodge(0.6, preserve = "single"), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Created Habitat\n(orginal habitat)",   # Change legend title
                  values = c("#073b4c", "#ffd166", "#06d6a0",  "#ef476f", "#118ab2")) +
ylab("Breeding Wader Pairs/ £100,000") +
xlab("Lawton Principle") +
labs(fill = "Targeting Strategy") +
guides(fill=guide_legend(nrow=2, byrow=F)) +
BarPlotTheme


## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_FencingEffect.png"), units = "in", height = 9, width = 11)





##------------------------------------------------------##
##------------------------------------------------------##
##  Plot Scenarios for Land Purchase vs Income Foregone ##
##------------------------------------------------------##
##------------------------------------------------------##

## Calculate change in breeding pairs for reserve creation on grassland and pay income foregone
ScenSumCost_FG <- Set |>
filter(!OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Created = "Reserve- Grassland Income Foregone",           
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")

## Calculate change in breeding pairs for reserve creation on grassland and purchase land
ScenSumCost_PU <- Set |>
filter(!OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCostPU)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCostPU)) |> 
ungroup() |> 
mutate(Created = "Reserve- Purchase Grassland",           
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")


## Calculate change in breeding pairs for reserve creation on arable land and pay income foregone
ScenSumCost_FGArable <- Set |>
filter(OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Created = "Reserve- Arable Income Foregone",           
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")

## Calculate change in breeding pairs for reserve creation on arable land and purchase land
ScenSumCost_PUArable <- Set |>
filter(OppCat == "Arable Opp") |> 
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCostPU)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCostPU)) |> 
ungroup() |> 
mutate(Created = "Reserve- Purchase Arable",           
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "Reserve")


## Calculate change in breeding pairs for AES creation
ScenSumCost_AES <- Set |>
mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
       PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(Perc_PairCost)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          Ave_PairCost = mean(Perc_PairCost)) |> 
ungroup() |> 
mutate(Created = "AES Only",           
       LowerCI = Ave_PairCost - (t_score * StEr),
       UpperCI = Ave_PairCost + (t_score * StEr)) |> 
filter(NewCat == "AES Only")


## Combine all data sets for plotting
CostComp <- rbind(ScenSumCost_FG, ScenSumCost_PU, ScenSumCost_AES, ScenSumCost_FGArable, ScenSumCost_PUArable)


## Create Bar plot
ggplot(CostComp, aes(x= Strategy, y= Ave_PairCost*100000, fill = Created)) +
geom_col(width=0.55, position=position_dodge(0.6, preserve = "single")) +
geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=0.2, position=position_dodge(0.6, preserve = "single"), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Created Habitat-\n(land costs)",   # Change legend title
                  values = c("#073b4c", "#ffd166", "#06d6a0",  "#ef476f", "#118ab2")) +
ylab("Breeding Wader Pairs/ £100,000") +
xlab("Lawton Principle") +
labs(fill = "Targeting Strategy") +
guides(fill=guide_legend(nrow=2, byrow=F)) +
BarPlotTheme


## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_ForegonePurchase.png"), units = "in", height = 9, width = 12)




##-------------------------------------##
##-------------------------------------##
##  Plot Scenario Average Cluster Size ##
##-------------------------------------##
##-------------------------------------##  


## First create a data set that can be used to create a bar plot that covers all scenarios
ScenSum <- Set |>
mutate(PlotCatFull = paste0(Strategy, " for ", NewCat, ifelse(OppCat== "Arable Opp", "\nfrom Arable", "")),
       PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
group_by(PlotCatFull, ScenType, PlusFull, Plus, NewCat, Strategy) |>
summarise(StEr = sd(ClustActual)/sqrt(n()),
          t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
          ClustActual= mean(ClustActual)) |> 
ungroup() |> 
mutate(LowerCI = ClustActual - (t_score * StEr),
       UpperCI = ClustActual + (t_score * StEr))

## Read in the 
Sizes <- read_csv("CleanData/Scenarios/4-AnnotateCanvas/Average_reserve_AES_size.csv")

## Create Bar plot
ggplot(ScenSum, aes(x= PlotCatFull, y= ClustActual, fill = ScenType)) +
geom_col(width=0.55, position=position_dodge(0.55)) +
geom_errorbar(aes(ymin= LowerCI, ymax= UpperCI), width=0.2, position=position_dodge(0.55), colour = "#6a6b6b") +
facet_wrap(~PlusFull) +
scale_fill_manual(name = "Scenario Type",   # Change legend title
                    labels = c("clusterlarge"= "Cluster-large", "clustersmall"= "Cluster-small", "random" ="Random"),  # Change legend labels
                    values = c("clusterlarge"="#A5F076", "clustersmall"="#76A5F0", "random"="#F076A5")) +
geom_hline(yintercept=Sizes$ReserveAve, linetype="dashed", color = "grey", linewidth = 1.2) +
geom_hline(yintercept=Sizes$AESAve, linetype="dashed", color = "grey", linewidth = 1.2) +
ylab("Average cluster size/ ha") +
xlab("Scenario Category") +
labs(fill = "Targeting Strategy") +
BarPlotTheme + 
theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))

## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_AveClusterSizeha_vs_Category.png"), units = "in", height = 9, width = 11)



##------------------------------------##
##------------------------------------##
##  Plot Scenarios Lapwing & Redshank ##
##------------------------------------##
##------------------------------------##

## Calculate the change in Lapwing pairs per unit cost
ScenSumLap <- Set |>
  filter(!Landscape == "Somerset Levels and Moors") |> 
  mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
         PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
         PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
  group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
  summarise(StEr = sd(Perc_PairCostLap)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(Perc_PairCostLap)) |> 
  ungroup() |> 
  mutate(Species = "Lapwing",
         LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr))

## Calculate the change in Redshank pairs per unit cost
ScenSumRed <- Set |>
  filter(!Landscape == "Somerset Levels and Moors") |> 
  mutate(Arable2 = ifelse(OppCat== "Arable Opp", "\nfrom Arable", ""),
         PlotCatFull = paste0(Strategy, " for ", NewCat, Arable2),
         PlusFull = ifelse(Plus==FALSE, "Average-Quality", "High-Quality")) |>
  group_by(PlotCatFull, PlusFull, Plus, NewCat, Strategy) |>
  summarise(StEr = sd(Perc_PairCostRed)/sqrt(n()),
            t_score = qt(p=0.05/2, df=(n()-1), lower.tail=F),
            Ave_PairCost = mean(Perc_PairCostRed)) |> 
  ungroup() |> 
  mutate(Species = "Redshank",
         LowerCI = Ave_PairCost - (t_score * StEr),
         UpperCI = Ave_PairCost + (t_score * StEr))
  

## Combine the data sets for Lapwing and Redshank
SpeciesComb <- rbind(ScenSumLap, ScenSumRed)


## Create Bar plot
ggplot(SpeciesComb, aes(x= PlotCatFull, y= Ave_PairCost*100000, fill = Species)) +
  geom_col(width=0.55, position=position_dodge(0.6)) +
  geom_errorbar(aes(ymin= LowerCI*100000, ymax= UpperCI*100000), width=0.2, position=position_dodge(0.6), colour = "#6a6b6b") +
  facet_wrap(~PlusFull) +
  scale_fill_manual(name = "Species",   # Change legend title
                    values = c("#6464f4", "#f46464")) +
  ylab("Breeding Wader Pairs/ £100,000") +
  xlab("Scenario Category") +
  BarPlotTheme + 
  theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1))


## save the plot
ggsave(plot=last_plot(), filename= paste0(outpath, "Additive_Pair£_vs_LapwingRedshankComparison.png"), units = "in", height = 9, width = 11)
  


