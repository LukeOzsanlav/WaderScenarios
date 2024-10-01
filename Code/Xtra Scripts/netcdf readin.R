## Luke Ozsanlav-Harris
## Read in the Future flows climate data set
## Create indices of the time series that show regional water availability

pacman::p_load(tidyverse, sf, terra, ncdf4)

## For 20 years of data for the whole of the UK for precip and evapotranpiration
## There are nearly 13 billion values and then there are 11 different simulations..... (143 billion)
## But for each pixel I want the monthly net water gain so this will decrease the final data size by a factor of 24

PrecipFile <- "RawData/FutureFlows_ClimateData/FF-HadRM3-afgcx-APr-2010-2039.nc"



our_nc_data <- nc_open(PrecipFile)
print(our_nc_data)


attributes(our_nc_data$dim)
# Get latitude and longitude with the ncvar_get function and store each into their own object:
Northing <- ncvar_get(our_nc_data, "Northing")
nNorthing <- dim(Northing) #to check it matches the metadata: 23
Easting <- ncvar_get(our_nc_data, "Easting")
nEasting <- dim(Easting) #to check, should be 24
# Check your lat lon dimensions match the information in the metadata we explored before:
print(c(nEasting, nNorthing))
# Get the time variable. Remember: our metadata said our time units are in seconds since 1981-01-01 00:00:00, so you will not see a recognizable date and time format, but a big number like "457185600". We will take care of this later
time <- ncvar_get(our_nc_data, "Time")
head(time) # just to have a look at the numbers
tunits <- ncatt_get(our_nc_data, "Time", "units") #check units
nt <- dim(time) #should be 2622



attributes(our_nc_data$var)
#get the variable in "matrix slices"
AP_array <- ncvar_get(our_nc_data, "avlble_pr") 

fillvalue <- ncatt_get(our_nc_data, "avlble_pr", "_FillValue")
dim(AP_array) #to check; this should give you 24 23 2622
#right away let's replace the nc FillValues with NAs
AP_array[AP_array==fillvalue$value] <- NA
AP_array



v3      <- our_nc_data$var[[1]]
varsize <- v3$varsize
ndims   <- v3$ndims
nt      <- varsize[ndims]
i=1000
#for( i in 1:nt ) {
	# Initialize start and count to read one timestep of the variable.
	start <- rep(1,ndims)	# begin with start=(1,1,1,...,1)
	start[ndims] <- i	# change to start=(1,1,1,...,i) to read timestep i
	count <- varsize	# begin w/count=(nx,ny,nz,...,nt), reads entire var
	count[ndims] <- 20	# change to count=(nx,ny,nz,...,1) to read 1 tstep
	data3 <- ncvar_get( our_nc_data, v3, start=start, count=count, raw_datavals =T)
  print(data3)
  data3[250,601,10]
	# Now read in the value of the timelike dimension
	timeval <- ncvar_get( our_nc_data, v3$dim[[ndims]]$name, start=i, count=1 )

	print(paste("Data for",v3$name,"at timestep",i,
		"(time value =",timeval,v3$dim[[ndims]]$units,")"))
	
#	}

test <-	terra::rast(PrecipFile)
test <- aggregate(test[[1:30]], by = c(rep(1, 30)), fun = sum)
plot(test)





araw <- tidync(PrecipFile) %>%
         hyper_filter(Time = Time > 32870) %>%
  hyper_array()

araw[["avlble_pr"]]




