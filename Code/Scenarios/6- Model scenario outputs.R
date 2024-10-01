##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 11/09/2024
## 
## Aim: Model abundance differences between scenarios
## 
##------------------------------------------------------##

## Load in packages & helper functions
pacman::p_load(tidyverse, # used
               data.table,
               sf) # used

options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R") # additional simple functions 






##------------------------------##
#### 1. Read in Scenario data ####
##------------------------------##

## Create empty list to put the output data for each region into
ListSet <- list()

## Create a vector of input paths for each region
Paths <- c("CleanData/Scenarios/5-ScenarioCreation/Broads/",
           "CleanData/Scenarios/5-ScenarioCreation/Essex/",
           #"CleanData/Scenarios/5-ScenarioCreation/Kent/",
           "CleanData/Scenarios/5-ScenarioCreation/Som/")

## Create vector of the species that are used for each region
Species <- c("LapRed", "LapRed", "LapRed", "Snipe")


## Loop through each region and calculate a series of columns needed for plotting
for(j in 1:length(Paths)){ 
  
  ## number of years to spread one off costs over
  Yr=15
  
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
                     WaderTot = rowSums(dplyr::select(., Snipe, Lapwing, Redshank), na.rm = T),
                     BaseWaderTot = rowSums(dplyr::select(., BaseSnipe, BaseLapwing, BaseRedshank), na.rm = T),
                     ChangeSnipe = Snipe-BaseSnipe,
                     ChangeSnipenoF = Snipe_unfenced-BaseSnipe,
                     ChangeLapwing = Lapwing-BaseLapwing,
                     ChangeRedshank = Redshank-BaseRedshank,
                     ChangeWaders = (Lapwing+Redshank) - (BaseLapwing+BaseRedshank),
                     ChangeWadersnoF = (Lapwing_unfenced+Redshank_unfenced) - (BaseLapwing+BaseRedshank),
                     ChangeCostsFG = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                     ChangeCostsPU = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/Yr)+(PurchaseCost/Yr))-BaseResMaintCost, 0)),
                     ChangeCostsNoF = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+ForegoneCost)-BaseResMaintCost, 0)),
                     PlotCat = paste0(NewCat, ifelse(Plus==T, "+", "-"), " ", ifelse(OppCat== "Arable Opp", "Arable", "")),
                     Plus = ifelse(Plus==T, "Plus", "NonePlus")) |>
              filter(!SegmentArea == 0) 

  
  
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
Set <- do.call("rbind", ListSet)






##----------------------------------##
#### 2. Model Change in Abundance ####
##----------------------------------##

## Load in a packages used for modelling and model checks
library(glmmTMB)
library(DHARMa)
library(lme4)
library(performance)
library(emmeans)
library(effects)

## Look at the structure of the data set
glimpse(Set)



##-- Organize data for models --##

## Select columns for model and calculate change in abundance
SetD <- Set |> select(WaderTot, BaseWaderTot, ChangeCostsFG, ScenType, Strategy, NewCat, Plus, Landscape, OppCat) |> 
  mutate(AbChange = WaderTot-BaseWaderTot) |> 
  filter(!OppCat == "Arable Opp")


## Set the right format for each column
SetD <- SetD |>
        mutate(across(c(WaderTot, ChangeCostsFG, AbChange), as.numeric),
                   across(c(ScenType, BaseWaderTot, Strategy, NewCat, Plus, Landscape), as.factor)) %>% 
        as.data.frame()
summary(is.na(SetD)) # check for NAs
summary(SetD) # summary of data

## Check for data point across different factor levels
with(SetD, 
  as.data.frame(table(ScenType, Strategy, NewCat, Plus, Landscape)))
scale(SetD$WaderTot)




##-- Running Model --##

## baseline: random cluster, better
Mod1 <- glmmTMB(AbChange ~ BaseWaderTot*ChangeCostsFG + ScenType + Strategy + NewCat + Plus,
                data = SetD,
                family = t_family(link = "identity")) # Chose t-family as data has long tail

summary(Mod1)
confint(Mod1)

glimpse(SetD)
## baseline: random cluster, better, non-plus
Mod2 <- glmmTMB(AbChange ~ BaseWaderTot + ScenType + Strategy + NewCat + Plus + 
                           ScenType*Strategy*Plus*NewCat,
                data = SetD,
                #offset = scale(ChangeCostsFG),
                family = t_family(link = "identity"))

summary(Mod2)
confint(Mod2)


noise.emm <- emmeans(Mod2, ~ ScenType*Strategy*Plus*NewCat)

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")

joint_tests(Mod2)



#, lines=list(multiline=TRUE)


##-- Model Diagnostics --##

## Check model using DHARMa
simulationOutput <- simulateResiduals(fittedModel = Mod1, plot = F)
plot(simulationOutput)

## plot the residuals against other predictors
par(mfrow=c(3,2))
plotResiduals(simulationOutput, form = Set$BaseWaderTot) # not violated
plotResiduals(simulationOutput, form = Set$ChangeCostsFG) # not violated
plotResiduals(simulationOutput, form = Set$ScenType) # not violated
plotResiduals(simulationOutput, form = Set$Strategy) # not violated
plotResiduals(simulationOutput, form = Set$NewCat) # not violated
plotResiduals(simulationOutput, form = Set$Plus) # not violated

## check model using performance package
par(mfrow=c(1,1))
check_model(Mod1)










