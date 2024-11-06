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
           "CleanData/Scenarios/5-ScenarioCreation/Kent/",
           "CleanData/Scenarios/5-ScenarioCreation/Som/")

## Create vector of the species that are used for each region
Species <- c("LapRed", "LapRed", "LapRed", "Snipe")


## Loop through each region and calculate a series of columns needed for plotting
for(j in 1:length(Paths)){ 
  
  ## number of years to spread one off costs over
  Yr=20
  
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
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/15)+ForegoneCost)-BaseResMaintCost, 0)),
                     ChangeCostsPU = ifelse(NewCat=="AES Only", NewAESCost - BaseAESCost,
                                          ifelse(NewCat=="Reserve", (NewResMaintCost+(NewResCreatCost/Yr)+(FencingCost/15)+(PurchaseCost/Yr))-BaseResMaintCost, 0)),
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




##-- Organize data for models --##

## Look at the structure of the data set
glimpse(Set)

## Select columns for model and calculate change in abundance
SetD <- Set |> select(WaderTot, BaseWaderTot, ChangeCostsFG, ScenType, Strategy, NewCat, Plus, Landscape, OppCat) |> 
  mutate(AbChange = WaderTot-BaseWaderTot,
         NewCat = ifelse(OppCat == "Arable Opp", paste0(NewCat, "_Arable"), NewCat))


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

## Re-level certain factors for better plotting
SetD <- SetD |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random"))





##-- Running Model --##

## baseline: random cluster, better
Mod1 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ Landscape + ScenType + Strategy + NewCat + Plus,
                data = SetD,
                family = t_family(link = "identity")) # Chose t-family as data has long tail

summary(Mod1)
confint(Mod1)


##-- Plot model outputs --##


## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(Mod1))
rownames(Ests) <- c("Intercept", "Region: Broads", "Region: Essex", "Region: North Kent", "Targetting: Large Sites", "Targetting: Small Sites ", "Location: Bigger", 
                    "Location: More", "Mechanism: Reserve", "Mechanism: Reserve from Arable", "Quality: High-quality")
Ests <- Ests[!rownames(Ests)=="Intercept",] # remove the intercept
Ests <- Ests |> as.data.frame()

Ests$Pos <- c(8, 9, 10, 3, 4, 1, 2, 5, 6, 7)
Ests <- Ests[order(Ests$Pos), ]


## make a forest plot fo the model estimates
Ests <- Ests %>%
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))


## this vector might be useful for other plots/analyses
level_order <- c(Ests$Param)


## Create forest plot
ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#E74C3C")) +
        theme(panel.grid.minor.y = element_blank(),
              #panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 

## save the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_ForestPlot.png", units = "in", height = 8, width = 11)
  



##-- Comparisons --##

glimpse(SetD)
emmeans(Mod1, specs = pairwise ~ ScenType + Strategy)
emm <- emmeans(Mod1, ~ ScenType + Strategy)
simp <- pairs(emm, simple = "each")
pairs

(Emm_ScenType <- emmeans(Mod1, "ScenType"))
confint(pairs(Emm_ScenType))
emmeans(Mod1, ~ ScenType) |> plot()

(Emm_Strategy <- emmeans(Mod1, "Strategy"))
confint(pairs(Emm_Strategy))
emmeans(Mod1, ~ Strategy) |> plot()
 
(Emm_NewCat <- emmeans(Mod1, "NewCat"))
pairs(Emm_NewCat)
emmeans(Mod1, ~ NewCat) |> plot()

(Emm_Plus <- emmeans(Mod1, "Plus"))
pairs(Emm_Plus)
emmeans(Mod1, ~ Plus) |> plot()




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






##---------------------------------##
#### 3. Model Regional Abundance ####
##---------------------------------##


##-- Model for Somerset --##

## baseline: random cluster, better
Som <- SetD |> filter(Landscape == "Somerset Levels and Moors")
ModS <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ ScenType + Strategy + NewCat + Plus,
                data = Som,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModS)

(Emm_ScenType <- emmeans(ModS, "ScenType"))
confint(pairs(Emm_ScenType))
emmeans(ModS, ~ ScenType) |> plot()

(Emm_Strategy <- emmeans(ModS, "Strategy"))
confint(pairs(Emm_Strategy))
emmeans(ModS, ~ Strategy) |> plot()
 
(Emm_NewCat <- emmeans(ModS, "NewCat"))
confint(pairs(Emm_NewCat))
emmeans(ModS, ~ NewCat) |> plot()

(Emm_Plus <- emmeans(ModS, "Plus"))
confint(pairs(Emm_Plus))
emmeans(ModS, ~ Plus) |> plot()



## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(ModS))
## make a forest plot fo the model estimates
Ests <- Ests %>%
        as.data.frame() |> 
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))
## Create forest plot
ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= Param, x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#E74C3C")) +
        theme(panel.grid.minor.y = element_blank(),
              #panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 




##-- Model for Broads --##

## baseline: random cluster, better
Broads <- SetD |> filter(Landscape == "Broads")
ModB <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ ScenType + Strategy + NewCat + Plus,
                data = Broads,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModB)

(Emm_ScenType <- emmeans(ModB, "ScenType"))
confint(pairs(Emm_ScenType))
emmeans(ModB, ~ ScenType) |> plot()

(Emm_Strategy <- emmeans(ModB, "Strategy"))
confint(pairs(Emm_Strategy))
emmeans(ModB, ~ Strategy) |> plot()
 
(Emm_NewCat <- emmeans(ModB, "NewCat"))
confint(pairs(Emm_NewCat))
emmeans(ModB, ~ NewCat) |> plot()

(Emm_Plus <- emmeans(ModB, "Plus"))
confint(pairs(Emm_Plus))
emmeans(ModB, ~ Plus) |> plot()


## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(ModB))
## make a forest plot fo the model estimates
Ests <- Ests %>%
        as.data.frame() |> 
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))
## Create forest plot
ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= Param, x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#E74C3C")) +
        theme(panel.grid.minor.y = element_blank(),
              #panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 




##-- Model for North Kent --##

## baseline: random cluster, better
NKent <- SetD |> filter(Landscape == "North Kent")
ModNK <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ ScenType + Strategy + NewCat + Plus,
                data = NKent,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModNK)

(Emm_ScenType <- emmeans(ModNK, "ScenType"))
confint(pairs(Emm_ScenType))
emmeans(ModNK, ~ ScenType) |> plot()

(Emm_Strategy <- emmeans(ModNK, "Strategy"))
confint(pairs(Emm_Strategy))
emmeans(ModNK, ~ Strategy) |> plot()
 
(Emm_NewCat <- emmeans(ModNK, "NewCat"))
confint(pairs(Emm_NewCat))
emmeans(ModNK, ~ NewCat) |> plot()

(Emm_Plus <- emmeans(ModNK, "Plus"))
pairs(Emm_Plus)
emmeans(ModNK, ~ Plus) |> plot()


## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(ModNK))
## make a forest plot fo the model estimates
Ests <- Ests %>%
        as.data.frame() |> 
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))
## Create forest plot
ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= Param, x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#E74C3C", "grey")) +
        theme(panel.grid.minor.y = element_blank(),
              #panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 




##-- Model for Essex --##

## baseline: random cluster, better
Essex <- SetD |> filter(Landscape == "Essex")
ModEs <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ ScenType + Strategy + NewCat + Plus,
                data = Essex,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModEs)

(Emm_ScenType <- emmeans(ModEs, "ScenType"))
confint(pairs(Emm_ScenType))
emmeans(ModEs, ~ ScenType) |> plot()

(Emm_Strategy <- emmeans(ModEs, "Strategy"))
confint(pairs(Emm_Strategy))
emmeans(ModEs, ~ Strategy) |> plot()
 
(Emm_NewCat <- emmeans(ModEs, "NewCat"))
pairs(Emm_NewCat)
emmeans(ModEs, ~ NewCat) |> plot()

(Emm_Plus <- emmeans(ModEs, "Plus"))
pairs(Emm_Plus)
emmeans(ModEs, ~ Plus) |> plot()


## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(ModEs))
## make a forest plot fo the model estimates
Ests <- Ests %>%
        as.data.frame() |> 
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))
## Create forest plot
ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= Param, x= Estimate, colour = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        scale_color_manual(values=c("#E74C3C", "grey")) +
        theme(panel.grid.minor.y = element_blank(),
              #panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 
