pacman::p_load(terra, tidyterra, tidyverse)


## Read in Somerset opportunity map
R <- rast("Outputs/Somerset_Opportunity_Map.tif")
# plot masked and smoother raster
ggplot() +
  geom_spatraster(data = R) + 
  scale_fill_viridis_c(na.value = "transparent") + theme_light()

## smooth the raster with focal windoe of 9 (9*200m = 1.8km)
Rsmooth <- focal(R, w=9, fun="mean", na.rm = T)
plot(Rsmooth) # plot smoothed raster
Rsmooth <- mask(Rsmooth, R) # crop to the extent of the orginal raster

# plot masked and smoother raster
ggplot() +
  geom_spatraster(data = Rsmooth) + 
  scale_fill_viridis_c(na.value = "transparent") + theme_light()



## Find the maximum value within the 13-cell neighborhood of each cell
f <- function(X) max(X, na.rm=TRUE)
ww <- matrix(1, nrow=13, ncol=13) ## Weight matrix for cells in moving window
localmax <- focal(R, fun=f, w=ww, expand=TRUE, fillvalue=NA)

## Does each cell have the maximum value in its neighborhood?
r2 <- R==localmax
plot(r2)

## Get x-y coordinates of those cells that are local maxima
maxXY <- xyFromCell(r2, which(values(r2)==1))
head(maxXY)

## create mask for opportunity map to see what values these peaks had in the opportunity maps
r3 <- r2
values(r3) <- ifelse(values(r3)==FALSE, NA, values(r3))
Rspeak <- mask(R, r3)
ggplot() +
  geom_spatraster(data = Rspeak) + 
  scale_fill_viridis_c(na.value = "transparent") + theme_light()

## mask any peaks that have a value less than 10
values(Rspeak) <- ifelse(values(Rspeak) < 10, NA, values(Rspeak))
plot(Rspeak)

## order the none-NA pixels, so give peaks an order
Rspeak$order <- (tibble(layer = Rspeak$layer[]) %>% mutate(rownumber = 1:nrow(.)) %>% arrange(-layer) %>% mutate(ordernew = 1:nrow(.)) %>% arrange(rownumber))$ordernew
Rspeak <- mask(Rspeak, Rspeak$layer)
plot(Rspeak$order)

# plot the peaks using the tidyterra package
ggplot() +
  geom_spatraster(data = Rspeak, aes(fill = order)) + 
  scale_fill_viridis_c(na.value = "transparent") + theme_light()
