#=========================================================================
# TREESCAPES: Participatory scenario planning: future land use scenario
#    
# AREA: NORTH PENNINES & DALES PRIORITY LANDSCAPE
#
# GROUP: BLUE
#=========================================================================

#============================================
# 1. set working directory & create functions
#============================================
setwd("D:/GIS")
source("setup.R")

# Function for calculating the distance of non-NA pixels to focal land cover, e.g. for calculating the distance from opportunity area (non-NA) to nearest woodland
dist_to_focal <- function(r, mask = NULL){
  # Target cells (e.g. woodland)
  target_cells <- which(r[] != 0)
  target_xy <- cbind(xyFromCell(r, target_cells))
  target_tree <- createTree(coordinates(target_xy))
  # Focal cells - non-NA 
  if(is.null(mask)){
    focal_cells <- which(!is.na(r[]))
  } else{
    focal_cells <- which(!is.na(mask[]))
  }
  focal_xy <- cbind(xyFromCell(r, focal_cells))
  ids <- knnLookup(target_tree, newdat= coordinates(focal_xy), k = 1)
  r[focal_cells] <- sqrt((focal_xy[, 1]-target_xy[ids, 1])^2 + (focal_xy[,2]-target_xy[ids, 2])^2)
  return(r)
}

# Function for getting object from string
get2 <- function(x){eval(parse(text = x))}

#======================================
# 2. CREATE OPPORTUNITY AREAS FOR EACH INTERVENTION
# 2a Woodland (riparian and upland gills)
# 2b Semi-natural grassland expansion
# 2c Edge broadleaves (standards on boundaries)
# 2d Remove conifer plantations (on peat soils - to bog)
# 2e Remove conifer plantations (on non peat soils - to broadleaved)
# 2f Scrub (riparian and upland gills)
# 2g Orchards
# 2h Peatland restoration
#======================================

## Read Landscape-specific raster stacks --------------
rast_NPD <- rast("rast/rast2_NPD.tif")
# Read LCM classes
lcm_classes <- read_csv("lcm_classes_edit.csv")
#rast_NPD$lcm <- rast_NPD$lcm+1 #check whether still need to do this
#levels(rast_NPD$lcm_edit) <- lcm_classes
levels(rast_NPD$landscape) <- "NPD"


#==========================================================
# 2a Mixed Woodland (riparian and upland gills - lower elevation)
#==========================================================
wood <- rast_NPD$landscape == "NPD" #actual landscape not buffer
wood[is.na(rast_NPD$rivers_buf)] <- NA
wood[!rast_NPD$lcm_edit %in% c("Improved grassland", "Acid grassland", "Calcareous grassland", "Neutral grassland", "Heather", "Heather grassland")] <- NA  #NA for not these land cover types
wood[!is.na(rast_NPD$peat)] <- NA # not in areas of peat
wood[rast_NPD$p_hab %in% c("Calaminarian grassland", "Limestone pavement", "Purple moor grass and rush pastures", "Upland calcareous grassland", "Upland flushes, fens and swamps", "Upland hay meadow")] <- NA    
wood[rast_NPD$waders > 10] <- NA
wood[!is.na(rast_NPD$bgrouse)] <- NA
wood[rast_NPD$treecover > 80] <- NA
wood[!is.na(rast_NPD$moorline)] <- NA

names(wood) <- "op"
plot(wood, col = "black")

#====================================
# 2b Semi natural grassland expansion
#====================================
grass <- rast_NPD$landscape == "NPD" #actual landscape not buffer
grass[!rast_NPD$lcm_edit %in% c("Improved grassland")] <- NA  #NA for not these land cover types
buff <- rast_NPD$landscape == "NPD" 
buff[!rast_NPD$p_hab %in% c("Upland hay meadow")] <- NA  #NA for not these land cover types
buff <- buff %>% 
  as.polygons()
buff <- buffer(buff, 1000) 
buff_rast <- rasterize(buff, rast_NPD)
grass[is.na(buff_rast)] <-NA # areas not in hay meadow buffer

names(grass) <- "op"

plot(grass, col = "black")

#===============================================
# 2c Edge broadleaves (standards on boundaries)
#===============================================
tree <- rast_NPD$landscape == "NPD" #actual landscape not buffer
tree[rast_NPD$boundary == 0] <- NA
tree[!rast_NPD$lcm_edit %in% c("Acid grassland", "Calcareous grassland", "Neutral grassland")] <- NA  #NA for not these land cover types
tree[rast_NPD$treecover > 10] <- NA

names(tree) <- "op"
plot(tree, col = "black")

#============================================
# 2d Remove conifer plantations (not peat)
#============================================

conifer <- rast_NPD$landscape == "NPD" #actual landscape not buffer
conifer[!rast_NPD$lcm_edit %in% c("Coniferous woodland")] <- NA
conifer[!is.na(rast_NPD$peat)] <- NA

names(conifer) <- "op"
plot(conifer, col = "black")

#============================================
# 2e Remove conifer plantations (peat)
#============================================

conifer_peat <- rast_NPD$landscape == "NPD" #actual landscape not buffer
conifer_peat[!rast_NPD$lcm_edit %in% c("Coniferous woodland")] <- NA
conifer_peat[is.na(rast_NPD$peat)] <- NA

names(conifer_peat) <- "op"
plot(conifer_peat, col = "black")

#================
# 2f Scrub
#================

uplandgill <- rast_NPD$landscape == "NPD"
uplandgill[is.na(rast_NPD$moorline)] <-NA # within the moorline
uplandgill[is.na(rast_NPD$rivers_buf)] <-NA # within river 50m buff
uplandgill[rast_NPD$slope < 5] <- NA

brackenslope <- rast_NPD$landscape == "NPD"
brackenslope[is.na(rast_NPD$bracken)] <- NA

brackengill <- merge(brackenslope, uplandgill)

scrub <- rast_NPD$landscape == "NPD" #actual landscape not buffer
scrub[is.na(brackengill)] <- NA
scrub[!rast_NPD$lcm_edit %in% c("Improved grassland", "Acid grassland", "Calcareous grassland", "Neutral grassland", "Heather", "Heather grassland")] <- NA  #NA for not these land cover types
scrub[!is.na(rast_NPD$peat)] <- NA # not in areas of peat
scrub[rast_NPD$p_hab %in% c("Calaminarian grassland", "Limestone pavement", "Purple moor grass and rush pastures", "Upland calcareous grassland", "Upland flushes, fens and swamps", "Upland hay meadow")] <- NA    
scrub[rast_NPD$waders > 10] <- NA
scrub[!is.na(rast_NPD$bgrouse)] <- NA
scrub[rast_NPD$treecover > 30] <- NA

names(scrub) <- "op"
plot(scrub, col = "black")

#================================
# 2g Silvopasture (fruit trees)
#================================
orchard <- rast_NPD$landscape == "NPD" #actual landscape not buffer
orchard[rast_NPD$p_hab %in% c("Calaminarian grassland", "Limestone pavement", "Purple moor grass and rush pastures", "Upland calcareous grassland", "Upland flushes, fens and swamps", "Upland hay meadow", "Traditional orchard")] <- NA    
orchard[!rast_NPD$lcm_edit %in% c("Improved grassland")] <- NA  #NA for not these land cover types
orchard[rast_NPD$waders > 10] <- NA
orchard[!is.na(rast_NPD$bgrouse)] <- NA
orchard[!is.na(rast_NPD$moorline)] <- NA
orchard[rast_NPD$treecover > 10] <- NA
orchard[!is.na(rast_NPD$peat)] <- NA # not in areas of peat

names(orchard) <- "op"
plot(orchard, col = "black")

#==============================
# 2h Peatland restoration
#==============================

peat <- rast_NPD$landscape == "NPD" #actual landscape not buffer
peat[is.na(rast_NPD$peat)] <- NA 
peat[!rast_NPD$lcm_edit %in% c("Improved grassland", "Acid grassland", "Calcareous grassland", "Neutral grassland", "Heather", "Heather grassland", "Degraded bog")] <- NA  #NA for not these land cover types
names(peat) <- "op"
plot(peat, col = "black")

#===========================================
# 3. Ranking methods of intervention creation
#===========================================

#including near to existing hay meadow
## Rank order
# Calculate distance to hay meadow
grass$hay_meadow_dist <- dist_to_focal(r = rast_NPD$p_hab == "Upland hay meadow", mask = grass$op)
grass$hay_meadow_dist[is.na(grass$op)] <- NA 
grass$tree <- rast_NPD$treecover
grass$tree[is.na(grass$op)] <- NA
# Calculate rank order 
grass$order <- (tibble(dist = grass$hay_meadow_dist[], tree = grass$tree[]) %>% mutate(rownumber = 1:nrow(.)) %>% arrange(dist, tree) %>% mutate(ordernew = 1:nrow(.)) %>% arrange(rownumber))$ordernew 
grass$order[is.na(grass$op)] <- NA


#calculate scrub/woodland highest ALC first
scrub$alc <- rast_NPD$alc
scrub$alc[is.na(scrub$op)] <- NA
scrub$tree <- rast_NPD$treecover
scrub$tree[is.na(scrub$op)] <- NA
scrub$order <- (tibble(alc = scrub$alc[], tree = scrub$tree[]) %>% mutate(rownumber = 1:nrow(.)) %>% arrange(-alc, tree) %>% mutate(ordernew = 1:nrow(.)) %>% arrange(rownumber))$ordernew 
scrub$order[is.na(scrub$op)] <- NA

wood$alc <- rast_NPD$alc
wood$alc[is.na(wood$op)] <- NA
wood$tree <- rast_NPD$treecover
wood$tree[is.na(wood$op)] <- NA
wood$order <- (tibble(alc = wood$alc[], tree = wood$tree[]) %>% mutate(rownumber = 1:nrow(.)) %>% arrange(-alc, tree) %>% mutate(ordernew = 1:nrow(.)) %>% arrange(rownumber))$ordernew 
wood$order[is.na(wood$op)] <- NA

# rank order for everything else near low tree cover

orchard$tree <- rast_NPD$treecover
orchard$tree[is.na(orchard$op)] <- NA
orchard$order <- rank(orchard$tree[] * orchard$op[], na.last = "keep", ties.method = "random") # lowest tree cover ranking

tree$tree <- rast_NPD$treecover
tree$tree[is.na(tree$op)] <- NA
tree$order <- rank(tree$tree[] * tree$op[], na.last = "keep", ties.method = "random") # lowest tree cover ranking

# tree$order_rand <- rank(tree$op[], na.last = "keep", ties.method = "random") #random ranking
conifer$order <- rank(conifer$op[], na.last = "keep", ties.method = "random") #random ranking
peat$order <- rank(peat$op[], na.last = "keep", ties.method = "random") #random ranking
conifer_peat$order <- rank(conifer_peat$op[], na.last = "keep", ties.method = "random") #random ranking

#==========================
# 4. Land change scenarios 
#==========================
# write opportunity rasters, clear files and load back in to save memory
writeRaster(wood, file = "rast/data/NPD_blue/wood_opp.tif", overwrite = TRUE)
writeRaster(scrub, file = "rast/data/NPD_blue/scrub_opp.tif", overwrite=TRUE)
writeRaster(tree, file = "rast/data/NPD_blue/tree_opp.tif", overwrite=TRUE)
writeRaster(conifer, file = "rast/data/NPD_blue/conifer_opp.tif", overwrite=TRUE)
writeRaster(peat, file = "rast/data/NPD_blue/peat_opp.tif", overwrite=TRUE)
writeRaster(conifer_peat, file = "rast/data/NPD_blue/conifer_peat_opp.tif", overwrite=TRUE)
writeRaster(orchard, file = "rast/data/NPD_blue/orchard_opp.tif", overwrite=TRUE)
writeRaster(grass, file = "rast/data/NPD_blue/grass_opp.tif", overwrite=TRUE)

#read back in opportunity areas
wood <- rast("rast/data/NPD_blue/wood_opp.tif")
scrub <- rast("rast/data/NPD_blue/scrub_opp.tif")
tree <- rast("rast/data/NPD_blue/tree_opp.tif")
conifer <- rast("rast/data/NPD_blue/conifer_opp.tif")
peat <- rast("rast/data/NPD_blue/peat_opp.tif")
conifer_peat <- rast("rast/data/NPD_blue/conifer_peat_opp.tif")
orchard <- rast("rast/data/NPD_blue/orchard_opp.tif")
grass <- rast("rast/data/NPD_blue/grass_opp.tif")


n_woodop <- sum(!is.na(wood$op[]))
n_scrubop <- sum(!is.na(scrub$op[]))
n_treeop <- sum(!is.na(tree$op[]))
n_coniferop <- sum(!is.na(conifer$op[]))
n_conifer_peatop <- sum(!is.na(conifer_peat$op[]))
n_peatop <- sum(!is.na(peat$op[]))
n_orchardop <- sum(!is.na(orchard$op[]))
n_grassop <- sum(!is.na(grass$op[]))


#load in original raster again 
rast_NPD <- rast("rast/rast2_NPD.tif")
#lcm_classes <- read_csv("lcm_classes_edit.csv")
#levels(rast_NPD$lcm_edit) <- lcm_classes

levels(rast_NPD$landscape) <- "NPD"

## Update LCM
# Create function
update_lcm_fun <- function(lcm_in, order, target, lcm_to){
  # Initialise new LCM
  lcm_new <- lcm_in     
  lcm_stay <- rast_NPD$lcm_edit
  lcm_stay[is.na(order)] <- NA
  
  if(!is.na(lcm_to)) lcm_new[order <= target] <- lcm_to
  #tcd_new[order <= target] <- tcd_to  #if order <= target is TRUE then tcd will change
  
  if(is.na(lcm_to)) lcm_new[order <= target] <- lcm_stay
  
  return(lcm_new)
}


# Create list of scenario parameters
params <- tibble(landscape = "NPD",
                 order_name_wood = c("alc"),
                 lcm_in = "rast_NPD$lcm_edit",
                 order_wood = c("wood$order"),
                 lcm_to_wood = 26) %>% 
  crossing(target_prop_wood = c(0.25, 0.5, 0.75, 1)) %>%
  mutate(target_prop_grass = target_prop_wood,       
         order_name_grass = "near_hay",
         order_grass = "grass$order",
         lcm_to_grass = 28) %>% 
  mutate(target_prop_conifer = target_prop_wood,       
         order_name_conifer = "rand",
         order_conifer = "conifer$order",
         lcm_to_conifer = 1) %>% 
  mutate(target_prop_tree = target_prop_wood,       
         order_name_tree = "rand",
         order_tree = "tree$order",
         lcm_to_tree = NA) %>% 
  mutate(target_prop_conifer_peat = target_prop_wood,       
         order_name_conifer_peat = "rand",
         order_conifer_peat = "conifer_peat$order",
         lcm_to_conifer_peat = 11) %>%
  mutate(target_prop_scrub = target_prop_wood,       
         order_name_scrub = "alc",
         order_scrub = "scrub$order",
         lcm_to_scrub = 25) %>%
  mutate(target_prop_peat = target_prop_wood,       
         order_name_peat = "rand",
         order_peat = "peat$order",
         lcm_to_peat = 11) %>%
  mutate(target_prop_orchard = target_prop_wood,       
         order_name_orchard = "rand",
         order_orchard = "orchard$order",
         lcm_to_orchard = 27) %>%
  mutate(target_wood = round(n_woodop * target_prop_wood),
         target_grass = round(n_grassop * target_prop_grass),
         target_conifer = round(n_coniferop * target_prop_conifer),
         target_conifer_peat = round(n_conifer_peatop * target_prop_conifer_peat),
         target_peat = round(n_peatop * target_prop_peat),
         target_tree = round(n_treeop * target_prop_tree),
         target_orchard = round(n_orchardop * target_prop_orchard),
         target_scrub = round(n_scrubop * target_prop_scrub)) %>% 
  distinct()


# Create a list of new landcover rasters    
lcm_new <- as.list(rep(NA, nrow(params)))

# Loop through params, update LCM, updates lcm with woodland first, then edits that output with other interventions
# if other mitigation measures are added then need to consider order (e.g. woodland last after grassland)
for(i in 1:length(lcm_new)){
  lcm_new[[i]] <- update_lcm_fun(get2(params[i, ]$lcm_in), get2(params[i, ]$order_peat), params[i, ]$target_peat, params[i, ]$lcm_to_peat)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_orchard), params[i, ]$target_orchard, params[i, ]$lcm_to_orchard)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_scrub), params[i, ]$target_scrub, params[i, ]$lcm_to_scrub)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_conifer_peat), params[i, ]$target_conifer_peat, params[i, ]$lcm_to_conifer_peat)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_conifer), params[i, ]$target_conifer, params[i, ]$lcm_to_conifer)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_tree), params[i, ]$target_tree, params[i, ]$lcm_to_tree)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_grass), params[i, ]$target_grass, params[i, ]$lcm_to_grass)
  lcm_new[[i]] <- update_lcm_fun(lcm_new[[i]],             get2(params[i, ]$order_wood), params[i, ]$target_wood, params[i, ]$lcm_to_wood)
}

# Stack

lcm_new <- rast(lcm_new)
lcm_new[is.na(rast_NPD$landscape)] <- NA
writeRaster(lcm_new, file = "rast/data/NPD_blue/lcm_new.tif", overwrite = T)
#lcm_new <- rast("rast/data/NPD_blue/lcm_new.tif")

#==========================
# 5. Tree density change 
#==========================


#update tree canopy density function
update_tcd_fun <- function(tcd_in, order, target, tcd_to){
  # Initialise new TCD
  tcd_new <- tcd_in     
  tcd_stay <- rast_NPD$treecover
  tcd_stay[is.na(order)] <- NA
  
  if(!is.na(tcd_to)) tcd_new[order <= target] <- tcd_to
  #tcd_new[order <= target] <- tcd_to  #if order <= target is TRUE then tcd will change
  
  if(is.na(tcd_to)) tcd_new[order <= target] <- tcd_stay
  
  return(tcd_new)
}


# Create list of scenario parameters

params_trees <- tibble(landscape = "NPD",
                 order_name_wood = c("alc"),
                 tcd_in = "rast_NPD$treecover",
                 order_wood = c("wood$order"),
                 tcd_to_wood = 80) %>% 
  crossing(target_prop_wood = c(0.25, 0.5, 0.75, 1)) %>%
  mutate(target_prop_grass = target_prop_wood,       
         order_name_grass = "near_hay",
         order_grass = "grass$order",
         tcd_to_grass = NA) %>% 
  mutate(target_prop_conifer = target_prop_wood,       
         order_name_conifer = "rand",
         order_conifer = "conifer$order",
         tcd_to_conifer = 80) %>% 
  mutate(target_prop_conifer_peat = target_prop_wood,       
         order_name_conifer_peat = "rand",
         order_conifer_peat = "conifer_peat$order",
         tcd_to_conifer_peat = 0) %>%
  mutate(target_prop_scrub = target_prop_wood,       
         order_name_scrub = "alc",
         order_scrub = "scrub$order",
         tcd_to_scrub = 30) %>%
  mutate(target_prop_peat = target_prop_wood,       
         order_name_peat = "rand",
         order_peat = "peat$order",
         tcd_to_peat = NA) %>%
  mutate(target_prop_tree = target_prop_wood,       
         order_name_tree = "rand",
         order_tree = "tree$order",
         tcd_to_tree = 10) %>%
  mutate(target_prop_orchard = target_prop_wood,       
         order_name_orchard = "rand",
         order_orchard = "orchard$order",
         tcd_to_orchard = 10) %>%
  mutate(target_wood = round(n_woodop * target_prop_wood),
         target_grass = round(n_grassop * target_prop_grass),
         target_conifer = round(n_coniferop * target_prop_conifer),
         target_conifer_peat = round(n_conifer_peatop * target_prop_conifer_peat),
         target_peat = round(n_peatop * target_prop_peat),
         target_orchard = round(n_orchardop * target_prop_orchard),
         target_tree = round(n_treeop * target_prop_tree),
         target_scrub = round(n_scrubop * target_prop_scrub)) %>% 
  distinct()



# Create a list of new tree density rasters    
tcd_new <- as.list(rep(NA, nrow(params_trees)))

# Loop through params, update TCD with standards first, then the other interventions in the same order as the LCM change
for(i in 1:length(tcd_new)){
 tcd_new[[i]] <- update_tcd_fun(get2(params_trees[i, ]$tcd_in), get2(params_trees[i, ]$order_peat), params_trees[i, ]$target_peat, params_trees[i, ]$tcd_to_peat)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_orchard), params_trees[i, ]$target_orchard, params_trees[i, ]$tcd_to_orchard)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_scrub), params_trees[i, ]$target_scrub, params_trees[i, ]$tcd_to_scrub)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_conifer_peat), params_trees[i, ]$target_conifer_peat, params_trees[i, ]$tcd_to_conifer_peat)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_conifer), params_trees[i, ]$target_conifer, params_trees[i, ]$tcd_to_conifer)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_tree), params_trees[i, ]$target_tree, params_trees[i, ]$tcd_to_tree)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_grass), params_trees[i, ]$target_grass, params_trees[i, ]$tcd_to_grass)
 tcd_new[[i]] <- update_tcd_fun(tcd_new[[i]],             get2(params_trees[i, ]$order_wood), params_trees[i, ]$target_wood, params_trees[i, ]$tcd_to_wood)
}


# Stack
tcd_new <- rast(tcd_new)
tcd_new[is.na(rast_NPD$landscape)] <- NA
rast_NPD$treecover[is.na(rast_NPD$landscape)] <- NA
writeRaster(tcd_new, file = "rast/data/NPD_blue/tcd_new.tif", overwrite = T)

#calculate total canopy cover in new and present tcd maps

global(tcd_new, sum, na.rm=T)
global(rast_NPD$treecover, sum, na.rm=T) #divide by 10000 to get hectares

#tcd_new<- rast("rast/data/NPD_blue/tcd_new.tif")
names(tcd_new[[1]]) <- "25%"
names(tcd_new[[2]]) <- "50%"
names(tcd_new[[3]]) <- "75%"
names(tcd_new[[4]]) <- "100%"
png("rast/data/NPD_blue/tcd_new.png", width = 1200, height = 1000, pointsize = 20)
plot(tcd_new)
dev.off()

#===============
# 6. Plots
#===============

#plot opportunity areas
plot(peat$op, col = "brown")
plot(orchard$op, col = "lightblue", add =T)
plot(scrub$op, col = "green", add =T)
plot(conifer_peat$op, col = "orange", add =T)
plot(conifer$op, col = "darkgreen", add = T)
plot(tree$op, col = "red", add = T)
plot(grass$op, col = "blue", add = T)
plot(wood$op, col = "purple", add = T)

lcm_new[rast_NPD$lcm_edit == lcm_new] <- NA
names(lcm_new[[1]]) <- "25%"
names(lcm_new[[2]]) <- "50%"
names(lcm_new[[3]]) <- "75%"
names(lcm_new[[4]]) <- "100%"
png("rast/data/NPD_blue/lcm_new.png", width = 1400, height = 1000, pointsize = 20)
plot(lcm_new, col = c("darkgreen", "purple", "yellow", "lightgreen", "red", "blue"))
dev.off()


# Tidy
lcm_new[is.na(rast_NPD$landscape)] <- NA
rast_NPD$lcm_edit[is.na(rast_NPD$landscape)] <- NA



## Crosstabulate areas  #how many pixels under each land cover
# Initialise
areas <- as.list(rep(NA, nrow(params)))
# Loop through each scenario
for(i in 1:length(areas)){
  
  areas[[i]] <- crosstab(c(rast_NPD$lcm_edit, lcm_new[[i]]), long = TRUE, useNA = TRUE) %>% 
    as_tibble() %>% 
    set_names(c("lcm_from", "lcm_to", "freq")) %>% 
    mutate(area = freq * 0.01, # 10m pixel = 0.01 ha
           lcm_to = lcm_to,
           lcm_from = lcm_from) %>% # Fix factor levels 
    select(-freq) %>% 
    arrange(-area) %>%
    filter(!is.na(lcm_from)) %>% 
    crossing(params %>% 
               slice(i) %>% 
               select(landscape, target_prop_scrub
               ))
}

# Bind
areas <- bind_rows(areas)

current_land <- areas %>% 
  rename("Ambition" = target_prop_scrub) %>% 
  group_by(lcm_from, Ambition) %>% 
  summarise(total_area = sum(area)) 

write.csv(current_land, "rast/data/NPD_blue/current_land_cover_NPD.csv")

future_land <- areas %>% 
  rename("Ambition" = target_prop_scrub) %>% 
  group_by(lcm_to, Ambition) %>% 
  summarise(total_area = sum(area)) %>% 
  mutate(Group = "Blue")

write.csv(future_land, "rast/data/NPD_blue/future_land_cover_NPD.csv")

