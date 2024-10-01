

## Load in required packages
pacman::p_load(here, tidyverse, data.table, coda, stars, spAbundance, bayesplot, MCMCvis, caret)
options(scipen = 100, digits = 4) # set notation
source("Code/Helper functions.R") # get helper functions
set.seed(1234) # makes results reproducble



##----------------------##
#### 0.1 Data read in ####
##----------------------##

## Read in filtered breeding pairs estimates
Waders <- read_csv("CleanData/Script 4/Breeding_Pairs_FullAttrib.csv")
 
## Read in data on field characteristics from survey
FieldChar <- read_csv("CleanData/Script 1/Field_Characteristics_Clean.csv")



##------------------------##
#### 0.2 Join data sets ####
##------------------------##

## Join together the Wader data sets and the field characteristics
WData <- left_join(Waders, FieldChar, by = c("F_LOC_ID", "year"))
summary(WData)



##---------------------------##
#### 0.3 Add extra columns ####
##---------------------------##

## Add additional columns to the data set for modelling
WData <- WData |> 
         mutate(FieldArea = FieldArea/10000, # convert field area to hectares
                CorvDens = MagDens + CrowDens, # sum corvid denisty
                Reserve = ifelse(RSPB_Reserve == "Y" | LNR == "Y" | NNR == "Y", "Y", "N"), # reserve indicator variable
                RUSH_PERCENT = ifelse(is.na(RUSH_PERCENT)==T, 0, RUSH_PERCENT), # change NAs in rush cover to 0
                Fence_Coverage = ifelse(Fence_Coverage == 0, "N", "Y"), # fence indicator variable
                Lap_Density = est_pairsL/FieldArea, # breeding Lapwing density
                Red_Density = est_pairsR/FieldArea, # breeding Redshank density
                Sni_Density = est_pairsS/FieldArea,)# breeding Snipe density

## Plot histogram of breeding bird densities
hist(WData$est_pairsL, main = "Breeding Lapwing density")
hist(WData$est_pairsR, main = "Breeding Redshank density")
hist(WData$Sni_Density, main = "Breeding Snipe density")


##
#### Sort out list structure for spAbundacne models ####
##

## filter out just the Lapwing data for the Broads
LapBroad <- filter(WData, is.na(est_pairsL)==F)
colnames(LapBroad)
LapBroad <- mutate(LapBroad, est_pairsL = round2(est_pairsL, digits = 0),
                   S_LOC_ID = as.numeric(as.factor(LapBroad$S_LOC_ID)), 
                   Landscape = as.numeric(as.factor(LapBroad$Landscape)))

## Split the data set into a train-test set, randomly choose rows
# set.seed(1012)
# inTrain <- createDataPartition(
#   y = LapBroad$est_pairsL, # the outcome data
#   p = 0.15, # The percentage of data in the training set
#   list = FALSE)
# 
# ## Separate out the subset of data
# LapBroad <- LapBroad[ inTrain,]


## Know make my data into a list for modelling
LapBroad <- LapBroad |> rename(X = FieldX, Y = FieldY)
hist(LapBroad$est_pairsL)
Lap_list <- list(y = LapBroad$est_pairsL,
                 covs = LapBroad |> select(PropWood_500, PropIntense_500, PropWetGrass_500,
                                           ESS_Wader, CSS_Wader, Reserve, FieldArea, Landscape) |> as.data.frame(), 
                 coords = cbind(X=LapBroad$X, Y=LapBroad$Y))
str(Lap_list)





## This is the formal for the model
## Scaling the continuous variables will help with model convergence
Lap_formula <- ~ scale(PropWood_500) + scale(PropIntense_500) + scale(PropWetGrass_500) +
                 ESS_Wader + CSS_Wader + Reserve + 
                 scale(PropWood_500)*CSS_Wader + scale(PropIntense_500)*CSS_Wader + scale(PropWetGrass_500)*CSS_Wader +
                 scale(PropWood_500)*Reserve + scale(PropIntense_500)*Reserve + scale(PropWetGrass_500)*Reserve +
                 (1|Landscape) + offset(FieldArea)

## Specify the specific family we will use to model the data
Lap_family <- 'NB'

# Pair-wise distances between all sites
dist.mat <- dist(Lap_list$coords)
# Exponential covariance model
cov.model <- 'exponential'
n.neighbors <- 15
search.type <- 'cb'

min.dist <- min(dist.mat)
max.dist <- max(dist.mat)
priors <- list(beta.normal = list(mean = 25, var = 75), # tightened this prior based off first run
               kappa.unif = c(0, 100),
               sigma.sq.ig = c(2, 1),
               phi.unif = c(0.0001, 0.05), # tightened this prior based off first run
               sigma.sq.mu.ig = list(0.1, 0.1))


## In this modelling approach, we break up the total number of MCMC samples into a set of “batches”, 
## where each batch has a specific number of MCMC samples. 
## Thus, we must specify the total number of batches (n.batch) 
## as well as the number of MCMC samples each batch contains (batch.length) when specifying the function arguments. 
batch.length <- 25
n.batch <- 1000
# Total number of MCMC samples per chain
batch.length * n.batch

## Parameters for the MCMC 
tuning <- list(beta = 0.5, kappa = 0.5, beta.star = 0.5, w = 0.5, phi = 0.5)
# accept.rate = 0.43 by default, so we do not specify it.
n.burn <- 15000
n.thin <- 20
n.chains <- 3

## Run spatial abundance model
out.sp <- spAbund(formula = Lap_formula,
                  data = Lap_list,
                  #inits = inits,
                  priors = priors,
                  n.batch = n.batch,
                  batch.length = batch.length,
                  tuning = tuning,
                  cov.model = cov.model,
                  NNGP = TRUE,
                  n.neighbors = n.neighbors,
                  search.type = search.type,
                  n.omp.threads = 20,
                  n.report = 200,
                  family = 'NB',
                  verbose = TRUE,
                  n.burn = n.burn,
                  n.thin = n.thin,
                  n.chains = n.chains)
waicAbund(out.sp)


# ## Parameters for the MCMC 
# priors2 <- list(beta.normal = list(mean = 0, var = 100),
#                kappa.unif = c(0, 100),
#                sigma.sq.mu.ig = list(0.1, 0.1))
# tuning2 <- list(beta = 0.5, kappa = 0.5, beta.star = 0.5)
# out <- abund(formula = Lap_formula,
#              data = Lap_list,
#              #inits = inits,
#              priors = priors2,
#              n.batch = n.batch,
#              batch.length = batch.length,
#              tuning = tuning2,
#              n.omp.threads = 1,
#              n.report = 200,
#              family = 'NB',
#              verbose = TRUE,
#              n.burn = n.burn,
#              n.thin = n.thin,
#              n.chains = n.chains)
# 
# waicAbund(out)



summary(out.sp) # Regression coefficients
plot(out.sp, param = 'beta', density = FALSE) # MCMC chain traces
MCMCplot(out.sp$beta.samples, ref_ovl = T, ci = c(50, 92)) # plot covariate model estimates

## covariate model estimates using the bayesplot package
posterior <- as.array(out.sp)
dimnames(posterior)
color_scheme_set("red")
mcmc_intervals(posterior[["beta.samples"]]) # , pars = c("I", "drat", "am", "sigma")
mcmc_areas(
  posterior[["beta.samples"]], 
  prob = 0.92, # 92% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

## A posterior predictive check using a Freeman-Tukey test statistic, and summarize it with a Bayesian p-value
ppc.out.sp <- ppcAbund(out.sp, fit.stat = 'freeman-tukey', group = 0)
summary(ppc.out.sp)

## plot of the fitted values versus the trues to get a better sense of how our model performed
# Extract fitted values
y.rep.samples <- fitted(out.sp)
# Get means of fitted values
y.rep.means <- apply(y.rep.samples, 2, mean)
# Simple plot of True vs. fitted values
plot(Lap_list$y, y.rep.means, pch = 19, xlab = 'True', ylab = 'Fitted')
abline(0, 1)

color_scheme_set("brightblue")
ppc_dens_overlay(Lap_list$y, as.matrix(t(y.rep.means))) 







