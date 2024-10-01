##------------------------------------------------------##
##
## Luke Ozsanlav-Harris
## Created: 04/01/2023
## Helper Function for this Project
## 
##------------------------------------------------------##


## Function to get mode of a vector (v), excluding NAs
getmode <- function(v) {
   uniqv <- unique(v)[ !is.na(unique(v)) ]
   uniqv[which.max(tabulate(match(v, uniqv)))]} # function to calculate the mode



##-----------------------------------------------##



## nice little wrapper for plotting results of parameter tuning in the randomForestSRC package
plot.tune <- function(o, linear = TRUE) {
    x <- o$results[,1]
    y <- o$results[,2]
    z <- o$results[,3]
    so <- interp(x=x, y=y, z=z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]
    y0 <- y[idx]
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette =
                     colorRampPalette(c("yellow", "red")),
                   xlab = "nodesize",
                   ylab = "mtry",
                   main = "error rate for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                                points(x,y,pch=16,cex=.25)})
}



##-----------------------------------------------##



## rounding function that rounds up .5 values
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}



##-----------------------------------------------##



## Function to scale a variable within the bounds of 0 to 1
## The largest value in the vector will become 1 and the smallest 0
scale_vals <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}




## Function to inverse scale a variable within the bounds of 0 to 1
## The largest value in the vector will become 0 and the smallest 1
Inv_scale_vals <- function(x){ abs((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))-1) }



##-----------------------------------------------##



## Take a SpatRast and reclassify it so that any values below 0.5 are assigned to 0 
## and any values above 0.5 are assigned a value of 1
Half_reclass <- function(x){
  
 m = c(0, 0.5, 0, 0.5, 1, 1) # matrix for re-classification
 rclmat = matrix(m, ncol=3, byrow=TRUE) # make matrix
 x = classify(x, rclmat, include.lowest=TRUE) 
 return(x) 
  
}



##-----------------------------------------------##


## Function: Overlap ID'er ##

## Fields = polygons of the fields with breeding pair estimates from the three priority Landscapes
## OverlapShapes = polygons of the thing thing that you want to determine if the fields overlap with
## Thresh = the threshold for proportion overlap, above which fields are considered to be in OverlapShapes
## ColName = The name of column used to indicate whether the field is within OverlapShapes
Field_Overlap_IDer <- function(Fields = Fields, OverlapShapes = OverlapShapes, 
                               Thresh = Thresh, ColName = ColName){
  
  ## Crop OverlapShapes to region covered by Fields
  OverlapShapes <- st_crop(OverlapShapes, st_as_sf(Fields))
  plot(OverlapShapes$geometry)
  
  ## Ensure that the fields with breeding pair estimates is an sf object
  Fields <- st_as_sf(Fields)
  class(Fields)
  
  ## Work out which fields are in Reserves and what the size of the overlap is
  Overlap_sub <- Fields |> filter(Somerset %in% c("All", "S")) |> 
                 st_intersection(OverlapShapes) |> 
                 mutate(OverlapArea = st_area(geometry)) |> 
                 select(F_LOC_ID, OverlapArea) |> 
                 st_drop_geometry()
  
  ## Remove any duplicated in this data set, can't work out why duplicates happen, maybe overlapping shapes?
  Dups <- Overlap_sub |> select(F_LOC_ID, OverlapArea) |> duplicated()
  Overlap_sub <- Overlap_sub[Dups==FALSE,]
  Overlap_sub <- Overlap_sub |> group_by(F_LOC_ID) |> summarise(OverlapArea = sum(OverlapArea))

  ## Now join the overlap size onto the main fields data &
  ## Work out the proportion of each field in a reserve
  Fields <- left_join(Fields, Overlap_sub, by = "F_LOC_ID") |> 
             mutate(Area = st_area(geometry),
                    PropOverlap = ifelse(is.na(OverlapArea)==T, 0, OverlapArea/Area), 
                    Inside = ifelse(PropOverlap >= Thresh, "Y", "N")) |> 
             select(-c(Area, OverlapArea, PropOverlap)) 
  
  ## Plot just to check this 
  ggplot() + geom_sf(data = Fields, mapping = aes(geometry = geometry, fill = Inside)) + theme_minimal()
  
  ## change the column Name to one of your choosing
  colnames(Fields)[colnames(Fields)== "Inside"] <- ColName
  
  
  ## return this object
  return(Fields)
  
}




##-----------------------------------------------##

## Convert form m2 to hectares

m2_to_ha <- function(x){y <- as.numeric(x/10000)
                        return(y)}















