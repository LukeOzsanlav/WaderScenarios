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
                                                      Waders_100Ha = (ChangeWaders/SegmentArea)*100) }
  
  if(Species[j] == "Snipe"){ ListSet[[j]] <- AllScn |> mutate(ChangeWaders = ChangeSnipe,
                                                     PairCost = ChangeSnipe/(ChangeCosts/Budget),
                                                     PairCost100k = ChangeWaders/(ChangeCosts/100000),
                                                     Waders_100Ha = (ChangeWaders/SegmentArea)*100) }

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
SetD <- Set |> select(PairCost, Landscape, BaseWaderTot, Budget, Strategy, NewCat, OppCat) |> 
  mutate(Mechanism = ifelse(OppCat == "Arable Opp", paste0(NewCat, "_", OppCat), paste0(NewCat)))


## Set the right format for each column
SetD <- SetD |>
        mutate(across(c(PairCost), as.numeric),
               across(c(Landscape, BaseWaderTot, Budget, Strategy, NewCat, OppCat, Mechanism), as.factor)) %>% 
        as.data.frame()
summary(is.na(SetD)) # check for NAs
summary(SetD) # summary of data

## Check for data point across different factor levels
# with(SetD, 
#   as.data.frame(table(ScenType, Strategy, NewCat, Plus, Landscape)))
# scale(SetD$WaderTot)




##-- Running Model --##

## baseline: random cluster, better
Mod1 <- glmmTMB(PairCost ~ Landscape + Budget + Strategy + Mechanism + 
                  Landscape*Strategy + Budget*Strategy + Landscape*Mechanism + Budget*Mechanism,
                data = SetD,
                family = t_family(link = "identity")) # Chose t-family as data has long tail

summary(Mod1)
confint(Mod1)

glimpse(SetD)


(Emm_ScenType <- emmeans(Mod1, "Landscape"))
emmeans(Mod1, ~ Landscape) |> plot()

(Emm_Strategy <- emmeans(Mod1, "Budget"))
emmeans(Mod1, ~ Budget) |> plot()
 
(Emm_NewCat <- emmeans(Mod1, "Strategy"))
emmeans(Mod1, ~ Strategy) |> plot()

(Emm_Plus <- emmeans(Mod1, "Mechanism"))
emmeans(Mod1, ~ Mechanism) |> plot()



noise.emm <- emmeans(Mod1, ~ Landscape*Strategy, type = "response")

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")

noise.emm <- emmeans(Mod1, ~ Budget*Strategy)

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")

noise.emm <- emmeans(Mod1, ~ Budget*Mechanism)

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")

noise.emm <- emmeans(Mod1, ~ Landscape*Mechanism)

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")

noise.emm <- emmeans(Mod1, ~ Budget + Mechanism + Landscape + Strategy)

contrast(noise.emm, "pairwise", simple = "each", combine = T, adjust = "mvt")




emm_s.t <- emmeans(Mod1, pairwise ~  Strategy | Landscape)
emm_s.t

emm_s.t <- emmeans(Mod1, pairwise ~  Strategy | Budget)
emm_s.t 

emm_t.t <- emmeans(Mod1, pairwise ~ Mechanism | Landscape)
emm_t.t 

emm_t.t <- emmeans(Mod1, pairwise ~ Mechanism | Budget)
emm_t.t 



## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(Mod1))
rownames(Ests) <- c("Intercept", "Region: Broads", "Region: Thames", "Targetting: Small sites", "Targetting: Random ", "Location: Bigger", 
                    "Location: More", "Mechanism: Reserve", "Mechanism: Higher quality")

## make a forest plot fo the model estimates
Ests <- Ests %>%
        as.data.frame() |> 
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))


## Wexford forest plot
ggplot(Ests) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y=Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y=Param, x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#B2BABB", "#E74C3C")) +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 





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










