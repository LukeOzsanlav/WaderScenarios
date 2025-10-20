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
               sf, 
               ggpubr) # used

options(scipen=999) # turn off scientific notation
source("Code/Helper functions.R") # additional simple functions 

## Emmeans vignette:
##- https://rvlenth.github.io/emmeans/articles/AQuickStart.html



## Set a general theme for the bar plots
BarPlotTheme <- theme_light() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        plot.title=element_text(size=18, face = "bold", hjust = 0.5))



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






##-----------------------------##
#### 2. Model Abundance-Cost ####
##-----------------------------##

## Load in a packages used for modelling and model checks
library(glmmTMB)
library(DHARMa)
library(lme4)
library(performance)
library(emmeans)
library(effects)




## |- Organize data for models ----

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
         ScenType= fct_relevel(ScenType, "random"),
         Strategy= fct_relevel(Strategy, c("More", "Better", "Big")))




## |- Running Model ----

# ## baseline: random cluster, better
# Mod1 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ Landscape + ScenType + Strategy + NewCat + Plus,
#                 data = SetD,
#                 family = t_family(link = "identity")) # Chose t-family as data has long tail
# 
# summary(Mod1)
# confint(Mod1)


## Run model
Mod1 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  Landscape + ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = SetD,
                family = t_family(link = "identity")) # Chose t-family as data has long tail

summary(Mod1)
confint(Mod1)
# Strategy*NewCat




## |-  Plot model outputs ----


## Forest plot from model
## get the estimates and confidence intervals from the models
Ests <- cbind(confint(Mod1))
rownames(Ests) <- c("Intercept", "Region (Broads)", "Region (Essex)", "Region (North Kent)", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                    "Location (Bigger)", "Mechanism (Reserve)", "Mechanism (Reserve - arable)", "Quality (High)",
                    "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                    "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                    "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                    "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                    "Scale (Large site) * Mechanism (Reserve - arable)", "Scale (Small sites) * Mechanism (Reserve - arable)",
                    "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                    "Mechanism (Reserve) * Quality (High)", "Mechanism (Reserve - arable) * Quality (High)")
Ests <- Ests |> as.data.frame()
Ests$Pos <- c(1, 9, 10, 11, 4, 5, 2, 3, 6, 7, 8, 12:25) # define an order for the variables
Ests <- Ests[order(Ests$Pos), ]


## make a forest plot fo the model estimates
Ests <- Ests %>%
        mutate(Param = rownames(Ests),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#FF2A00", "#FFC3B8"))


## this vector might be useful for other plots/analyses
level_order <- c(Ests$Param)


## Create forest plot
FP <- ggplot(Ests) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
        theme_bw() +
        ylab("") +
        xlab("Coefficient Estimate") +
        ggtitle("Cost") +
        scale_color_manual(values=c("#FFC3B8", "#FF2A00")) +
        theme(panel.grid.minor.y = element_blank(),
              plot.title=element_text(size=20, face = "bold", hjust = 0.5),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=18), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=12), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FP, filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_ForestPlot_Cost.png", units = "in", height = 11, width = 11)
  


## |-  Plot EMMs ----

## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(Mod1, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(Mod1, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(Mod1, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(Mod1, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="Reserve"]/Emm3$emmean[Emm3$Var=="AES Only"] ," times more effective in terms of cost than wader AES ")
message("high-quality management producing, on average, ", Emm4$emmean[2]/Emm4$emmean[1]," times as many waders than standard-quality management for the same investment (cost and land)")


## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Mechanism: Reserve\nfrom Arable", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 9, 7, 8, 11, 10)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
        mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
EP <- ggplot(Emms) +
        #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
        geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
        scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 17, 10, 10))) +
        theme_bw() +
        ylab("") +
        xlab("Estimated Marginal Mean") +
        ggtitle("All Regions: Cost") +
        xlim(0, 26) +
        theme(panel.grid.minor.y = element_blank(),
              plot.title=element_text(size=20, face = "bold", hjust = 0.5),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=18),
              legend.title=element_text(size=14),
              axis.text=element_text(size=16),
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none")

## save the plot
ggsave(plot=EP, filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_EstMargMeans_Cost.png", units = "in", height = 8, width = 11)




## |- Plot Interactions ----

## re set the factors for the model
SetDNew <- SetD |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
Mod1 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  Landscape + ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = SetDNew,
                family = t_family(link = "identity"))

## create interaction style plot for estimated marginal means
testp <- emmip(Mod1, ScenType ~ Strategy | NewCat | Plus, 
              CIs = T, tlab = "Scale", ylab= "Estimated marginal mean\n(Breeding pairs per £100,000)",
              xlab = "Location",
              CIarg = c(linewidth = 0.6))
testp + 
  theme_bw() +
  ggtitle("Cost") +
  scale_color_manual(values=c("#FFD8B8", "#FF9640", "#A34900"),
                     labels = c("Field-scale", "Several small sites", "Single large site")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid.minor.x = element_blank()) 


## read out the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/6-ModelOutputs/Cost_Interaction_Vis.png", units = "in", height = 8, width = 11)




## Get the emmeans across the four response variables
Costemmean <- emmeans(Mod1, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenCost <- Costemmean |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)

# Create the base plot
FullCost <- ggplot(ScenCost, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.55, position = position_dodge(0.55, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.55, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = max(ScenCost$asymp.UCL, na.rm = TRUE) + 2,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 5, y = max(ScenCost$asymp.UCL, na.rm = TRUE) + 2,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 8, y = max(ScenCost$asymp.UCL, na.rm = TRUE) + 2,
           label = "More", fontface = "bold", size = 5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    # labels = c("AES Only" = "AES only",
                    #            "Reserve" = "Reserve creation (grassland conversion)",
                    #            "Reserve from Arable" = "Reserve creation (arable reversion)"),
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0",
                               "Reserve creation\n(arable reversion)" = "#DDB24C")) +
  ylab("Change in Breeding Pairs per £100,000") +
  xlab("Scale") +
  #ylim(-1, 10) +
  labs(fill = "Targeting Strategy") +
  BarPlotTheme +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1)) +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field-scale",
    "Better & Several small sites" = "Several small sites",
    "Better & Single large site" = "Single large site",
    "Bigger & Field-scale" = "Field-scale",
    "Bigger & Several small sites" = "Several small sites",
    "Bigger & Single large site" = "Single large site",
    "More & Field-scale" = "Field-scale",
    "More & Several small sites" = "Several small sites",
    "More & Single large site" = "Single large site"))

## read out the plot
ggsave(plot=FullCost, filename= "CleanData/Scenarios/6-ModelOutputs/Cost_Interaction_VisFULL.png", units = "in", height = 8, width = 11)





## |-  Comparisons ----

glimpse(SetD)
emmeans(Mod1, specs = pairwise ~ ScenType + Strategy)
emm <- emmeans(Mod1, ~ ScenType + Strategy)
simp <- pairs(emm, simple = "each")

Emm_ScenType <-emmeans(Mod1, "ScenType")
confint(pairs(Emm_ScenType, adjust = "mvt"))
emmeans(Mod1, ~ ScenType) |> plot()

Emm_Strategy <- emmeans(Mod1, "Strategy")
confint(pairs(Emm_Strategy, adjust = "mvt"))
emmeans(Mod1, ~ Strategy) |> plot()
 
Emm_NewCat <- emmeans(Mod1, "NewCat")
confint(pairs(Emm_NewCat, adjust = "mvt"))
emmeans(Mod1, ~ NewCat) |> plot()

Emm_Plus <- emmeans(Mod1, "Plus")
confint(pairs(Emm_Plus))
emmeans(Mod1, ~ Plus) |> plot()




## |-  Model Diagnostics ----

# ## Check model using DHARMa
# simulationOutput <- simulateResiduals(fittedModel = Mod1, plot = F)
# plot(simulationOutput)
# 
# ## plot the residuals against other predictors
# par(mfrow=c(3,2))
# plotResiduals(simulationOutput, form = Set$BaseWaderTot) # not violated
# plotResiduals(simulationOutput, form = Set$ChangeCostsFG) # not violated
# plotResiduals(simulationOutput, form = Set$ScenType) # not violated
# plotResiduals(simulationOutput, form = Set$Strategy) # not violated
# plotResiduals(simulationOutput, form = Set$NewCat) # not violated
# plotResiduals(simulationOutput, form = Set$Plus) # not violated
# 
# ## check model using performance package
# par(mfrow=c(1,1))
# check_model(Mod1)





##-----------------------------##
#### 3. Model Abundance-Area ####
##-----------------------------##

## |- Organize data for models ----

## Look at the structure of the data set
glimpse(Set)

## Select columns for model and calculate change in abundance
SetA <- Set |> select(WaderTot, BaseWaderTot, ChangeCostsFG, ScenType, Strategy, NewCat, Plus, Landscape, OppCat, SegmentArea) |> 
  mutate(AbChange = WaderTot-BaseWaderTot,
         NewCat = ifelse(OppCat == "Arable Opp", paste0(NewCat, "_Arable"), NewCat))


## Set the right format for each column
SetA <- SetA |>
  mutate(across(c(WaderTot, ChangeCostsFG, AbChange, SegmentArea), as.numeric),
         across(c(ScenType, BaseWaderTot, Strategy, NewCat, Plus, Landscape), as.factor)) %>% 
  as.data.frame()
summary(is.na(SetA)) # check for NAs
summary(SetA) # summary of data

## Check for data point across different factor levels
with(SetA, 
     as.data.frame(table(ScenType, Strategy, NewCat, Plus, Landscape)))
# scale(SetA$WaderTot)

## Re-level certain factors for better plotting
SetA <- SetA |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random"),
         Strategy= fct_relevel(Strategy, c("More", "Better", "Big")))




## |- Running model ----

## Run model
ModA <- glmmTMB(((AbChange/SegmentArea)*100) ~ 
                  Landscape + ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = SetA,
                family = t_family(link = "identity")) # Chose t-family as data has long tail

summary(ModA)
confint(ModA)


## |-  Plot model outputs ----

## Forest plot from model
## get the estimates and confidence intervals from the models
EstsA <- cbind(confint(ModA))
rownames(EstsA) <- c("Intercept", "Region (Broads)", "Region (Essex)", "Region (North Kent)", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                    "Location (Bigger)", "Mechanism (Reserve)", "Mechanism (Reserve - arable)", "Quality (High)",
                    "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                    "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                    "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                    "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                    "Scale (Large site) * Mechanism (Reserve - arable)", "Scale (Small sites) * Mechanism (Reserve - arable)",
                    "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                    "Mechanism (Reserve) * Quality (High)", "Mechanism (Reserve - arable) * Quality (High)")
EstsA <- EstsA |> as.data.frame()
EstsA$Pos <- c(1, 9, 10, 11, 4, 5, 2, 3, 6, 7, 8, 12:25) # define an order for the variables
EstsA <- EstsA[order(EstsA$Pos), ]


## make a forest plot fo the model estimates
EstsA <- EstsA %>%
  mutate(Param = rownames(EstsA),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#0040FF", "#B8C9FF"))


## this vector might be useful for other plots/analyses
level_order <- c(EstsA$Param)


## Create forest plot
FPA <- ggplot(EstsA) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
  theme_bw() +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Area") +
  scale_color_manual(values=c("#0040FF", "#B8C9FF")) +
  theme(panel.grid.minor.y = element_blank(),
        #panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FPA, filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_ForestPlot_Area.png", units = "in", height = 11, width = 11)



## Combine the two plots above
EstCombo <- ggarrange(FP, FPA + theme(axis.text.y=element_blank()),
                       ncol=2, common.legend = F, widths = c(2, 1.25))

## save the plot
ggsave(plot=EstCombo, filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_ForestPlot_Combo.png", units = "in", height = 11, width = 16)




## |-  Plot EMMs ----

## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(ModA, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(ModA, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(ModA, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(ModA, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="Reserve"]/Emm3$emmean[Emm3$Var=="AES Only"] ," times more effective in terms of land area than wader AES ")
message("high-quality management producing, on average, ", Emm4$emmean[2]/Emm4$emmean[1]," times as many waders than standard-quality management for the same investment (cost and land)")


## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Mechanism: Reserve\nfrom Arable", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 9, 7, 8, 11, 10)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
  mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
EP <- ggplot(Emms) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
  scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 17, 10, 10))) +
  theme_bw() +
  ylab("") +
  xlab("Estimated Marginal Mean") +
  ggtitle("All Regions: Area") +
  xlim(0, 40) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        axis.text=element_text(size=16),
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

## save the plot
ggsave(plot=EP, filename= "CleanData/Scenarios/6-ModelOutputs/Parameter_EstMargMeans_Area.png", units = "in", height = 8, width = 11)




## |- Plot Interactions ----


## re set the factors for the model
SetANew <- SetA |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
ModA2 <- glmmTMB(((AbChange/SegmentArea)*100) ~ 
                  Landscape + ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = SetANew,
                family = t_family(link = "identity"))

## create interaction style plot for estimated marginal means
testA <- emmip(ModA2, ScenType ~ Strategy | NewCat | Plus, 
               CIs = T, tlab = "Scale", ylab= "Estimated marginal mean\n(Breeding pairs per 100ha)",
               xlab = "Location",
               CIarg = c(linewidth = 0.6))
testA + 
  theme_bw() +
  ggtitle("Area") +
  scale_color_manual(values=c("#FFD8B8", "#FF9640", "#A34900"),
                     labels = c("Field-scale", "Several small sites", "Single large site")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid.minor.x = element_blank()) 


## read out the plot
ggsave(plot=last_plot(), filename= "CleanData/Scenarios/6-ModelOutputs/Area_Interaction_Vis.png", units = "in", height = 8, width = 11)




## Get the emmeans across the four response variables
Areaemmean <- emmeans(ModA2, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenArea <- Areaemmean |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)

# Create the base plot
FullArea <- ggplot(ScenArea, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.55, position = position_dodge(0.55, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.55, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = max(ScenArea$asymp.UCL, na.rm = TRUE) + 2,
           label = "Better", fontface = "bold", size = 5) +
  annotate("text", x = 5, y = max(ScenArea$asymp.UCL, na.rm = TRUE) + 2,
           label = "Bigger", fontface = "bold", size = 5) +
  annotate("text", x = 8, y = max(ScenArea$asymp.UCL, na.rm = TRUE) + 2,
           label = "More", fontface = "bold", size = 5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    # labels = c("AES Only" = "AES only",
                    #            "Reserve" = "Reserve creation (grassland conversion)",
                    #            "Reserve from Arable" = "Reserve creation (arable reversion)"),
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0",
                               "Reserve creation\n(arable reversion)" = "#DDB24C")) +
  ylab("Change in Breeding Pairs per 100ha") +
  xlab("Scale") +
  #ylim(-1, 10) +
  labs(fill = "Targeting Strategy") +
  BarPlotTheme +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1)) +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field-scale",
    "Better & Several small sites" = "Several small sites",
    "Better & Single large site" = "Single large site",
    "Bigger & Field-scale" = "Field-scale",
    "Bigger & Several small sites" = "Several small sites",
    "Bigger & Single large site" = "Single large site",
    "More & Field-scale" = "Field-scale",
    "More & Several small sites" = "Several small sites",
    "More & Single large site" = "Single large site"))

## read out the plot
ggsave(plot=FullArea, filename= "CleanData/Scenarios/6-ModelOutputs/Area_Interaction_VisFULL.png", units = "in", height = 8, width = 11)


## Combine the two plots above
FullCombo <- ggarrange(FullCost, FullArea,
                       nrow=2, common.legend = T)

## save the plot
ggsave(plot=FullCombo, filename= "CleanData/Scenarios/6-ModelOutputs/Combo_Interaction_VisFULL.png", units = "in", height = 18, width = 14)





## |-  Comparisons ----

glimpse(SetA)
emmeans(ModA, specs = pairwise ~ ScenType + Strategy)
emm <- emmeans(ModA, ~ ScenType + Strategy)
simp <- pairs(emm, simple = "each")

Emm_ScenType <-emmeans(ModA, "ScenType")
confint(pairs(Emm_ScenType, adjust = "mvt"))
emmeans(ModA, ~ ScenType) |> plot()

Emm_Strategy <- emmeans(ModA, "Strategy")
confint(pairs(Emm_Strategy, adjust = "mvt"))
emmeans(ModA, ~ Strategy) |> plot()

Emm_NewCat <- emmeans(ModA, "NewCat")
confint(pairs(Emm_NewCat, adjust = "mvt"))
emmeans(ModA, ~ NewCat) |> plot()

Emm_Plus <- emmeans(ModA, "Plus")
confint(pairs(Emm_Plus))
emmeans(ModA, ~ Plus) |> plot()




## |-  Model Diagnostics ----

# ## Check model using DHARMa
# simulationOutput <- simulateResiduals(fittedModel = ModA, plot = F)
# plot(simulationOutput)
# 
# ## plot the residuals against other predictors
# par(mfrow=c(3,2))
# plotResiduals(simulationOutput, form = Set$BaseWaderTot) # not violated
# plotResiduals(simulationOutput, form = Set$ChangeCostsFG) # not violated
# plotResiduals(simulationOutput, form = Set$ScenType) # not violated
# plotResiduals(simulationOutput, form = Set$Strategy) # not violated
# plotResiduals(simulationOutput, form = Set$NewCat) # not violated
# plotResiduals(simulationOutput, form = Set$Plus) # not violated
# 
# ## check model using performance package
# par(mfrow=c(1,1))
# check_model(ModA)







##---------------------------------##
#### 3. Model Regional Abundance ####
##---------------------------------##

## Set a general theme for the bar plots
BarPlotThemeR <- theme_light() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top",
        plot.title=element_text(size=14, face = "bold", hjust = 0),
        plot.margin = margin(t = 0, b = 0))





## |-  Model for Essex ----

## baseline: random cluster, better
Essex <- SetD |> filter(Landscape == "Essex")
ModEs <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                   ScenType + Strategy + NewCat + Plus +
                   ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                   Strategy*Plus  + Plus*NewCat,
                 data = Essex,
                 family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModEs)


##  Plot model outputs 

## Forest plot from model
## get the estimates and confidence intervals from the models
EstEss <- cbind(confint(ModEs))
rownames(EstEss) <- c("Intercept", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                      "Location (Bigger)", "Mechanism (Reserve)", "Mechanism (Reserve - arable)", "Quality (High)",
                      "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                      "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                      "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                      "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                      "Scale (Large site) * Mechanism (Reserve - arable)", "Scale (Small sites) * Mechanism (Reserve - arable)",
                      "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                      "Mechanism (Reserve) * Quality (High)", "Mechanism (Reserve - arable) * Quality (High)")
EstEss <- EstEss |> as.data.frame()
EstEss$Pos <- c(1, 4, 5, 2, 3, 6, 7, 8, 9:22) # define an order for the variables
EstEss <- EstEss[order(EstEss$Pos), ]


## make a forest plot fo the model estimates
EstEss <- EstEss %>%
  mutate(Param = rownames(EstEss),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#FF2A00", "#FFC3B8"))


## this vector might be useful for other plots/analyses
level_order <- c(EstEss$Param)


## Create forest plot
FEss <- ggplot(EstEss) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
  theme_bw() +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Essex: Cost") +
  scale_color_manual(values=c("#FF2A00", "#FFC3B8")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FEss, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Essex_Parameter_ForestPlot_Cost.png", units = "in", height = 11, width = 11)



## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(ModEs, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(ModEs, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(ModEs, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(ModEs, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="Reserve"]/Emm3$emmean[Emm3$Var=="AES Only"]," times more effective in terms of land area than wader AES ")
message("high-quality management producing, on average, ", Emm4$emmean[2]/Emm4$emmean[1]," times as many waders than standard-quality management for the same investment (cost and land)")



## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Mechanism: Reserve\nfrom Arable", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 9, 7, 8, 11, 10)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
  mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
EsEP <- ggplot(Emms) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
  scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 17, 10, 10))) +
  theme_bw() +
  ylab("") +
  xlab("Estimated Marginal Mean") +
  ggtitle("Essex: Cost") +
  xlim(0, 48) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        axis.text=element_text(size=16),
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

## save the plot
ggsave(plot=EsEP, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Essex_Parameter_EstMargMeans_Cost.png", units = "in", height = 8, width = 11)





## re set the factors for the model
SetEsNew <- Essex |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
ModEs2 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                    ScenType + Strategy + NewCat + Plus +
                    ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                    Strategy*Plus  + Plus*NewCat,
                  data = SetEsNew,
                  family = t_family(link = "identity"))


## Get the emmeans across the four response variables
CostemmEs <- emmeans(ModEs2, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenCostEss <- CostemmEs |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)


# Create the base plot
FullCostEss <- ggplot(ScenCostEss, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.75, position = position_dodge(0.75, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.75, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Better", fontface = "bold", size = 3.5) +
  annotate("text", x = 5, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Bigger", fontface = "bold", size = 3.5) +
  annotate("text", x = 8, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "More", fontface = "bold", size = 3.5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0",
                               "Reserve creation\n(arable reversion)" = "#DDB24C")) +
  ylab("Change in Breeding Pairs per £100,000") +
  xlab("Scale") +
  ggtitle("Essex") +
  labs(fill = "Targeting Strategy") +
  ylim(-17, ( max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)) +
  BarPlotThemeR +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field",
    "Better & Several small sites" = "Small",
    "Better & Single large site" = "Large",
    "Bigger & Field-scale" = "Field",
    "Bigger & Several small sites" = "Small",
    "Bigger & Single large site" = "Large",
    "More & Field-scale" = "Field",
    "More & Several small sites" = "Small",
    "More & Single large site" = "Large"))


## read out the plot
ggsave(plot=FullCostEss, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Essex_Cost_Interaction_Bar.png", units = "in", height = 8, width = 11)



## |-  Model for Somerset ----

## baseline: random cluster, better
Som <- SetD |> filter(Landscape == "Somerset Levels and Moors")
ModS <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = Som,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModS)


##  Plot model outputs 

## Forest plot from model
## get the estimates and confidence intervals from the models
EstSom <- cbind(confint(ModS))
rownames(EstSom) <- c("Intercept", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                    "Location (Bigger)", "Mechanism (Reserve)", "Quality (High)",
                    "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                    "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                    "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                    "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                    "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                    "Mechanism (Reserve) * Quality (High)")
EstSom <- EstSom |> as.data.frame()
EstSom$Pos <- c(1, 4, 5, 2, 3, 6, 7, 8:18) # define an order for the variables
EstSom <- EstSom[order(EstSom$Pos), ]


## make a forest plot fo the model estimates
EstSom <- EstSom %>%
  mutate(Param = rownames(EstSom),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#FF2A00", "#FFC3B8"))


## this vector might be useful for other plots/analyses
level_order <- c(EstSom$Param)


## Create forest plot
FSom <- ggplot(EstSom) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
  theme_bw() +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Somerset Levels: Cost") +
  scale_color_manual(values=c("#FF2A00", "#FFC3B8")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FSom, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Som_Parameter_ForestPlot_Cost.png", units = "in", height = 11, width = 11)



## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(ModS, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(ModS, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(ModS, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(ModS, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="AES Only"]/Emm3$emmean[Emm3$Var=="Reserve"] ," times more effective in terms of land area than wader AES ")

## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 8, 7, 10, 9)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
  mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
SomEP <- ggplot(Emms) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
  scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 10, 10))) +
  theme_bw() +
  ylab("") +
  xlab("Estimated Marginal Mean") +
  ggtitle("Somerset Levels: Cost") +
  xlim(0, 10) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        axis.text=element_text(size=16),
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

## save the plot
ggsave(plot=SomEP, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Som_Parameter_EstMargMeans_Cost.png", units = "in", height = 8, width = 11)



## re set the factors for the model
SetSomNew <- Som |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
ModS2 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = SetSomNew,
                family = t_family(link = "identity"))


## Get the emmeans across the four response variables
CostemmSom <- emmeans(ModS2, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenCostSom <- CostemmSom |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)


# Create the base plot
FullCostSom <- ggplot(ScenCostSom, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.75, position = position_dodge(0.75, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.75, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = (max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)/2,
           label = "Better", fontface = "bold", size = 3.5) +
  annotate("text", x = 5, y = (max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)/2,
           label = "Bigger", fontface = "bold", size = 3.5) +
  annotate("text", x = 8, y = (max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)/2,
           label = "More", fontface = "bold", size = 3.5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    # labels = c("AES Only" = "AES only",
                    #            "Reserve" = "Reserve creation (grassland conversion)",
                    #            "Reserve from Arable" = "Reserve creation (arable reversion)"),
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0")) +
  ylab("Change in Breeding Pairs per £100,000") +
  xlab("Scale") +
  ggtitle("Somerset Levels") +
  labs(fill = "Targeting Strategy") +
  ylim(-17/2, ( max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)/2) +
  BarPlotThemeR +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field",
    "Better & Several small sites" = "Small",
    "Better & Single large site" = "Large",
    "Bigger & Field-scale" = "Field",
    "Bigger & Several small sites" = "Small",
    "Bigger & Single large site" = "Large",
    "More & Field-scale" = "Field",
    "More & Several small sites" = "Small",
    "More & Single large site" = "Large"))

## read out the plot
ggsave(plot=FullCostSom, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Som_Cost_Interaction_Bar.png", units = "in", height = 8, width = 11)






## |-  Model for Broads ----

## baseline: random cluster, better
Broads <- SetD |> filter(Landscape == "Broads")
ModB <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = Broads,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModB)


##  Plot model outputs 

## Forest plot from model
## get the estimates and confidence intervals from the models
EstBro <- cbind(confint(ModB))
rownames(EstBro) <- c("Intercept", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                      "Location (Bigger)", "Mechanism (Reserve)", "Mechanism (Reserve - arable)", "Quality (High)",
                      "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                      "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                      "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                      "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                      "Scale (Large site) * Mechanism (Reserve - arable)", "Scale (Small sites) * Mechanism (Reserve - arable)",
                      "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                      "Mechanism (Reserve) * Quality (High)", "Mechanism (Reserve - arable) * Quality (High)")
EstBro <- EstBro |> as.data.frame()
EstBro$Pos <- c(1, 4, 5, 2, 3, 6, 7, 8, 9:22) # define an order for the variables
EstBro <- EstBro[order(EstBro$Pos), ]


## make a forest plot fo the model estimates
EstBro <- EstBro %>%
  mutate(Param = rownames(EstBro),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#FF2A00", "#FFC3B8"))


## this vector might be useful for other plots/analyses
level_order <- c(EstBro$Param)


## Create forest plot
FBro <- ggplot(EstBro) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
  theme_bw() +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Broads: Cost") +
  scale_color_manual(values=c("#FF2A00", "#FFC3B8")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FBro, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Broads_Parameter_ForestPlot_Cost.png", units = "in", height = 11, width = 11)



## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(ModB, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(ModB, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(ModB, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(ModB, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="Reserve"]/Emm3$emmean[Emm3$Var=="AES Only"]," times more effective in terms of land area than wader AES ")
message("high-quality management producing, on average, ", Emm4$emmean[2]/Emm4$emmean[1]," times as many waders than standard-quality management for the same investment (cost and land)")


## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Mechanism: Reserve\nfrom Arable", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 9, 7, 8, 11, 10)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
  mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
BroEP <- ggplot(Emms) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
  scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 17, 10, 10))) +
  theme_bw() +
  ylab("") +
  xlab("Estimated Marginal Mean") +
  ggtitle("Broads: Cost") +
  xlim(0, 48) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        axis.text=element_text(size=16),
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

## save the plot
ggsave(plot=BroEP, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Broads_Parameter_EstMargMeans_Cost.png", units = "in", height = 8, width = 11)


## re set the factors for the model
SetBrNew <- Broads |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
ModBr2 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                   ScenType + Strategy + NewCat + Plus +
                   ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                   Strategy*Plus  + Plus*NewCat,
                 data = SetBrNew,
                 family = t_family(link = "identity"))


## Get the emmeans across the four response variables
CostemmBro <- emmeans(ModBr2, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenCostBro <- CostemmBro |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)


# Create the base plot
FullCostBro <- ggplot(ScenCostBro, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.75, position = position_dodge(0.75, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.75, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Better", fontface = "bold", size = 3.5) +
  annotate("text", x = 5, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Bigger", fontface = "bold", size = 3.5) +
  annotate("text", x = 8, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "More", fontface = "bold", size = 3.5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0",
                               "Reserve creation\n(arable reversion)" = "#DDB24C")) +
  ylab("Change in Breeding Pairs per £100,000") +
  xlab("Scale") +
  ggtitle("Broads") +
  labs(fill = "Targeting Strategy") +
  ylim(-17, ( max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)) +
  BarPlotThemeR +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field",
    "Better & Several small sites" = "Small",
    "Better & Single large site" = "Large",
    "Bigger & Field-scale" = "Field",
    "Bigger & Several small sites" = "Small",
    "Bigger & Single large site" = "Large",
    "More & Field-scale" = "Field",
    "More & Several small sites" = "Small",
    "More & Single large site" = "Large"))

## read out the plot
ggsave(plot=FullCostBro, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/Broads_Cost_Interaction_Bar.png", units = "in", height = 8, width = 11)







## |-  Model for North Kent ----

## baseline: random cluster, better
NKent <- SetD |> filter(Landscape == "North Kent")
ModNK <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                  ScenType + Strategy + NewCat + Plus +
                  ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                  Strategy*Plus  + Plus*NewCat,
                data = NKent,
                family = t_family(link = "identity")) # Chose t-family as data has long tail
summary(ModNK)


##  Plot model outputs 

## Forest plot from model
## get the estimates and confidence intervals from the models
EstNK <- cbind(confint(ModNK))
rownames(EstNK) <- c("Intercept", "Scale (Large site)", "Scale (Small sites)", "Location (Better)", 
                      "Location (Bigger)", "Mechanism (Reserve)", "Mechanism (Reserve - arable)", "Quality (High)",
                      "Scale (Large site) * Location (Better)", "Scale (Small sites) * Location (Better)",
                      "Scale (Large site) * Location (Bigger)", "Scale (Small sites) * Location (Bigger)",
                      "Scale (Large site) * Quality (High)", "Scale (Small sites) * Quality (High)",
                      "Scale (Large site) * Mechanism (Reserve)", "Scale (Small sites) * Mechanism (Reserve)",
                      "Scale (Large site) * Mechanism (Reserve - arable)", "Scale (Small sites) * Mechanism (Reserve - arable)",
                      "Location (Better) * Quality (High)", "Location (Bigger) * Quality (High)",
                      "Mechanism (Reserve) * Quality (High)", "Mechanism (Reserve - arable) * Quality (High)")
EstNK <- EstNK |> as.data.frame()
EstNK$Pos <- c(1, 4, 5, 2, 3, 6, 7, 8, 9:22) # define an order for the variables
EstNK <- EstNK[order(EstNK$Pos), ]


## make a forest plot fo the model estimates
EstNK <- EstNK %>%
  mutate(Param = rownames(EstNK),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "#FF2A00", "#FFC3B8"))


## this vector might be useful for other plots/analyses
level_order <- c(EstNK$Param)


## Create forest plot
FNK <- ggplot(EstNK) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param, level = rev(level_order)), xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param, level = rev(level_order)), x= Estimate, colour = Sig), size = 3) +
  theme_bw() +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("North Kent: Cost") +
  scale_color_manual(values=c("#FF2A00", "#FFC3B8")) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

## save the plot
## Current reference level is for
## Somerset levels, Location - More, Quality - Average, Scale - random field, Mechanism - AES
ggsave(plot=FNK, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/NKent_Parameter_ForestPlot_Cost.png", units = "in", height = 11, width = 11)


## Retrieve all the Emmeans
Emm1 <- as.data.frame(emmeans(ModNK, "ScenType"))
colnames(Emm1)[1] <- "Var"
Emm2 <- as.data.frame(emmeans(ModNK, "Strategy"))
colnames(Emm2)[1] <- "Var"
Emm3 <- as.data.frame(emmeans(ModNK, "NewCat"))
colnames(Emm3)[1] <- "Var"
Emm4 <- as.data.frame(emmeans(ModNK, "Plus"))
colnames(Emm4)[1] <- "Var"

## comparison of reserve from grassland and AES
message("reserve creation through grassland conversion was ", Emm3$emmean[Emm3$Var=="Reserve"]/Emm3$emmean[Emm3$Var=="AES Only"] ," times more effective in terms of land area than wader AES ")
message("high-quality management producing, on average, ", Emm4$emmean[2]/Emm4$emmean[1]," times as many waders than standard-quality management for the same investment (cost and land)")



## bind them all together
Emms <- rbind(Emm1, Emm2, Emm3, Emm4)
Emms$Var <- c("Scale: Random", "Scale: Large Site", "Scale: Small Sites" , "Location: More", "Location: Better",
              "Location: Bigger", "Mechanism: AES",  "Mechanism: Reserve", "Mechanism: Reserve\nfrom Arable", "Quality: Average-quality", "Quality: High-quality")

## Define the order I want
Emms$Pos <- c(6, 4, 5, 3, 1, 2, 9, 7, 8, 11, 10)
Emms <- Emms[order(Emms$Pos), ]

## make a forest plot fo the model estimates
Emms <- Emms %>%
  mutate(Param2 = Var)


## this vector might be useful for other plots/analyses
level_order2 <- c(Emms$Param)


## Create plot
NKEP <- ggplot(Emms) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y= factor(Param2, level = rev(level_order2)), xmin= asymp.LCL, xmax=asymp.UCL), height = 0.2, linewidth =0.5) +
  geom_point(aes(y= factor(Param2, level = rev(level_order2)), x= emmean, shape = factor(Param2, level = rev(level_order2))), colour = "#E74C3C", size = 3.5) +
  scale_shape_manual(values = rev(c(15, 15, 15, 8, 8, 8, 17, 17, 17, 10, 10))) +
  theme_bw() +
  ylab("") +
  xlab("Estimated Marginal Mean") +
  ggtitle("North Kent: Cost") +
  xlim(0, 48) +
  theme(panel.grid.minor.y = element_blank(),
        plot.title=element_text(size=20, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        axis.text=element_text(size=16),
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

## save the plot
ggsave(plot=NKEP, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/NKent_Parameter_EstMargMeans_Cost.png", units = "in", height = 8, width = 11)



## re set the factors for the model
SetNKNew <- NKent |> 
  mutate(Landscape= fct_relevel(Landscape, "Somerset Levels and Moors"), 
         ScenType= fct_relevel(ScenType, "random",  "clustersmall", "clusterlarge"),
         Strategy=as.character(Strategy),
         Strategy = ifelse(Strategy=="Big", "Bigger", Strategy),
         Plus = ifelse(Plus== "NonePlus", "Quality: Average", "Quality: High"),
         Strategy= fct_relevel(Strategy, c("Better", "More", "Bigger")),
         NewCat = case_when(NewCat == "AES Only" ~ "AES only",
                            NewCat  == "Reserve" ~ "Reserve creation\n(grassland conversion)",  
                            NewCat  == "Reserve_Arable" ~ "Reserve creation\n(arable reversion)"))

## Run model
ModNK2 <- glmmTMB((AbChange/(ChangeCostsFG/100000)) ~ 
                    ScenType + Strategy + NewCat + Plus +
                    ScenType*Strategy + ScenType*Plus + ScenType*NewCat +
                    Strategy*Plus  + Plus*NewCat,
                  data = SetNKNew,
                  family = t_family(link = "identity"))


## Get the emmeans across the four response variables
CostemmNK <- emmeans(ModNK2, ~ ScenType * Strategy * NewCat * Plus) |> as.data.frame()

## alter the data set
ScenCostNK <- CostemmNK |>
  mutate(ScenType = case_when(ScenType == "random" ~ "Field-scale",
                              ScenType  == "clustersmall" ~ "Several small sites",  
                              ScenType  == "clusterlarge" ~ "Single large site"),
         PlotXCat = paste0(Strategy, " & ", ScenType)) |> 
  filter(!(NewCat == "Reserve creation\n(arable reversion)" & Strategy == "Better"))


# Define x-axis positions for groupings
group_ranges <- data.frame(
  category = c("Better", "Bigger", "More"),
  xmin = c(0.5, 3.5, 6.5),  # approximate starting positions
  xmax = c(3.5, 6.5, 9.5)   # approximate ending positions
)


# Create the base plot
FullCostNK <- ggplot(ScenCostNK, aes(x = PlotXCat, y = emmean, fill = NewCat)) +
  # Add grey rectangles
  geom_rect(data = group_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, 
            fill = "grey", 
            alpha = c(0, 0.15, 0.4, 0, 0.15, 0.4)) +  # Increasing darkness
  ## add bars
  geom_col(width = 0.75, position = position_dodge(0.75, preserve = "single")) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                position = position_dodge(width = 0.75, preserve = "single"), 
                colour = "#6a6b6b") +
  
  # Add labels
  annotate("text", x = 2, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Better", fontface = "bold", size = 3.5) +
  annotate("text", x = 5, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "Bigger", fontface = "bold", size = 3.5) +
  annotate("text", x = 8, y = max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2,
           label = "More", fontface = "bold", size = 3.5) +
  
  # Continue with your existing layers
  facet_wrap(~Plus) +
  scale_fill_manual(name = "Mechanism",
                    values = c("AES only" = "#F076A5",
                               "Reserve creation\n(grassland conversion)" = "#76A5F0",
                               "Reserve creation\n(arable reversion)" = "#DDB24C")) +
  ylab("Change in Breeding Pairs per £100,000") +
  xlab("Scale") +
  ggtitle("North Kent") +
  labs(fill = "Targeting Strategy") +
  ylim(-17, ( max(ScenCostEss$asymp.UCL, na.rm = TRUE) + 2)) +
  BarPlotThemeR +
  scale_x_discrete(labels = c(
    "Better & Field-scale" = "Field",
    "Better & Several small sites" = "Small",
    "Better & Single large site" = "Large",
    "Bigger & Field-scale" = "Field",
    "Bigger & Several small sites" = "Small",
    "Bigger & Single large site" = "Large",
    "More & Field-scale" = "Field",
    "More & Several small sites" = "Small",
    "More & Single large site" = "Large"))

## read out the plot
ggsave(plot=FullCostNK, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/NKent_Cost_Interaction_Bar.png", units = "in", height = 8, width = 11)




## arrange all the plots into one
RegCat <- ggarrange(FullCostBro, FullCostSom + ylab("") + xlab(""), FullCostNK + ylab("") + xlab(""), FullCostEss + ylab("") + xlab(""), common.legend = T, 
                    nrow= 2, ncol = 2, legend = "top",
                    align='hv',
                    hjust=-4.5, vjust=3)

## save the plot
ggsave(plot = RegCat, filename= "CleanData/Scenarios/6-ModelOutputs/Regional/AllRegions_Cost_Interaction_Bar.png", units = "cm", height = 20.3, width = 21.6)

