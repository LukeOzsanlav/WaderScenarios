
# Note for Rob 18,05,22 --------------------------------------------------
# This is code I used to format and manipulate the original survey data files provided by BTO. 
# There were separate files for bird obs, habitat data, site level data etc which meant these needed to be combined. 
# Also only positive counts were included in most cases (i.e. no absence values) and multiple species other than the 5 core waders were included so needed subsetting out. 
# Estimated numbers of pairs also needed to be calculated following the standard OB methods.
# Also some code in here used to help with covariate creation (TWI, Landscape wader abundance)
# This code is mostly for info only as I've already run it for the 2021 data: analysis datasets are already sorted and saved in the data folder. 
# You shouldn't need to do this again for the 2021 data hopefully!
# I'm including it in case you end up needing to follow the same process for the 2022 data once it's in.


#### Packages work in R version 3.5.2 

#install.packages("glmmADMB",repos="http://glmmadmb.r-forge.r-project.org/repos")
library("plyr")
library("reshape")
library("reshape2")
library("lme4")
library("lmerTest")
library("R2admb")
library("glmmADMB")
library("bbmle")    # for comparing AICs
library("plotrix")
library("car")
library("rnrfa")
library("sp")
library("raster")
library("rgeos")
library("ggplot2")
library("glmmTMB") ## Note - this package MASKS glmmADMB objects - do not load if I want to run a glmmADMB model
library("boot")
library("emmeans")
library("multcomp")
library("DescTools")
library("rgdal")
library("snow")
library("FSA")
library("sf")
library("viridis")
library("magick")
#library("exactextractr")
#library("Rcpp")
#library("terra")
#install.packages("velox")
#library(car)
library("multcompView")
library("MuMIn")




###### FIELD-VISIT-LEVEL DATA 

## THIS IS ALL DONE AND UP TO DATE - DATA EXPORTED AS FIELD_VISIT_WADER_DAT_MASTER



##### STEP 1: add F_LOC_IDs and S_LOC_IDs to GIS data ####
## so it can be matched with the wader observations and habitat data. 
# SKIP THIS STEP ONCE IT HAS BEEN RUN ONCE AND THE NEW FIELD LAYER GIS CREATED


### Luke - during this first stage you will resolve any fields with missing uniq_codes. If none are missing then you can skip this stage

#bwwm_fieldsF4GIS <- read.csv("bwwm-fieldsF4GIS.csv") # field metadata

temp.meta.data <- read.csv("../../Data/BWWM 2022/Field code cross refs 2022/Fields_2022_GIS_tab.csv") # Read in a CSV of you BWWM field GIS

# remove leading zeros

temp.meta.data$Uniq_code_zeros_ex <- sub("^0+", "", temp.meta.data$Uniq_code)

length(unique(temp.meta.data$Uniq_code))
length(unique(temp.meta.data$Uniq_code_zeros_ex))
length(unique(temp.meta.data$FID))

## 17 fields are missing a Uni_code 

###### For the 17 fields that dont have a Uniq_code, these can be created from the site code and field number

which(temp.meta.data$Uniq_code == " ") ## NOTE 1 - 11 of the 17 fields have no infomation in the uni_code column

n_occur <- data.frame(table(temp.meta.data$Uniq_code)) ## NOTE 2 - whilst the remaining 6 are duplicated across multiple fields
n_occur[n_occur$Freq > 1,]

##### For the three sites where Greg C altered the field boundaries in 2022, calculate a new Uni_code using the site and field code

temp.meta.data[temp.meta.data$Code_2020 %in% c("7540", "B_2870", "10160"),]$Uniq_code_zeros_ex <- paste(temp.meta.data[temp.meta.data$Code_2020 %in% c("7540", "B_2870", "10160"),]$Code_2020, temp.meta.data[temp.meta.data$Code_2020 %in% c("7540", "B_2870", "10160"),]$Field_no, sep = "_")

length(unique(temp.meta.data$Uniq_code_zeros_ex)) # perfect, 3840 uniq_codes


###### Export data and reimport

#write.csv(temp.meta.data , "../../Data/BWWM 2022/Field code cross refs 2022/Fields_2022_GIS_tabV2.csv")

temp.meta.data <- read.csv("../../Data/BWWM 2022/Field code cross refs 2022/Fields_2022_GIS_tabV2.csv")



###### Use Lucy's 2021 cross referncing data to get F_LOC_IDs and S_LOC_IDs for most of the fields. Note, for this matching process to work, I need to use the 'Uniq_codes' with the leading zeros removed


##### Luke - you can skip this step, every field in your GIS should have a unique F_LOC_ID (but worth checking)


field.cross.refs <- read.csv("../../Data/BWWM 2021/Field code cross-referencing/FIELD ID CROSS_REF.csv")

head(field.cross.refs)

temp.meta.data$F_LOC_ID <- field.cross.refs$X[match(temp.meta.data$Uniq_code_zeros_ex, field.cross.refs$X.2)] 
temp.meta.data$S_LOC_ID <- field.cross.refs$X.3[match(temp.meta.data$Uniq_code_zeros_ex, field.cross.refs$X.2)] 

nrow(temp.meta.data) # 38430
length(unique(temp.meta.data$F_LOC_ID)) # 38412 - 18 fields dont have an F LOC ID

temp.meta.data[which(is.na(temp.meta.data$F_LOC_ID)),]$Uniq_code_zeros_ex ## all the NAs are from the the three sites that Greg adjusted in 2022

missing.F.and.S.LOC.IDs <- temp.meta.data[which(is.na(temp.meta.data$F_LOC_ID)),]$Uniq_code_zeros_ex


###### Mannually add the F_LOC_IDs to these fields using the submitted habitat data (see: C:\Data\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\RSPB\BWWM\BWWM project handover\Data\BWWM 2022\Original data exports from BTO 2022\bwwm-visits-field@13-09-2022)


bwwm_fields <- read.csv("../../Data/BWWM 2022/Original data exports from BTO 2022/bwwm-fields@13-09-2022.csv") #  field list (fields visted in 2022)
  
bwwm_fields$Unif_code <- paste(bwwm_fields$SITE_CODE, bwwm_fields$FIELD_NUMBER, sep = "_")
  
missing.F.and.S.LOC.IDs <- data.frame(UniF_code = missing.F.and.S.LOC.IDs, F_LOC_ID = NA, S_LOC_ID = NA)

missing.F.and.S.LOC.IDs$F_LOC_ID <- bwwm_fields$F_LOC_ID[match(missing.F.and.S.LOC.IDs$UniF_code, bwwm_fields$Unif_code)]
missing.F.and.S.LOC.IDs$S_LOC_ID <- bwwm_fields$S_LOC_ID[match(missing.F.and.S.LOC.IDs$UniF_code, bwwm_fields$Unif_code)]

##### Add the missing F_LOC_IDs and S_LOC_IDs back into the main data frame

temp.meta.data[which(is.na(temp.meta.data$F_LOC_ID)),]$F_LOC_ID <-  missing.F.and.S.LOC.IDs$F_LOC_ID
temp.meta.data[which(is.na(temp.meta.data$S_LOC_ID)),]$S_LOC_ID <-  missing.F.and.S.LOC.IDs$S_LOC_ID


#### Check that every field has a unique F_LOC_ID

nrow(temp.meta.data) # 38430
length(unique(temp.meta.data$F_LOC_ID)) # Good! every field has a unique F_LOC_ID
length(unique(temp.meta.data$S_LOC_ID))


#### Next, export this table to ARC GIS so I can add the correct F_LOC_IDs and S_LOC_IDs to every field. When in ARC GIS, link the two tables by their FID (but check this has worked!!)
names(temp.meta.data)
head(temp.meta.data)
#write.csv(temp.meta.data[c(2,3,4,157,158,159,160)] , "../../Data/BWWM 2022/Field code cross refs 2022/Fields_2022_GIS_tabV3.csv")



## Checks - do the UniF_codes from this years field Excel web export (downloaded by Greg) match up with the UniF_codes from my new 2022 field GIS table


table(bwwm_fields$F_LOC_ID %in% temp.meta.data$F_LOC_ID) ### all 14112 fields that were included in Gregs web export for this year are in the 2022 GIS

bwwm_fields$Unif_code_zeros_ex <- sub("^0+", "", bwwm_fields$Unif_code)

identical(bwwm_fields$Unif_code_zeros_ex, temp.meta.data$Uniq_code_zeros_ex[match(bwwm_fields$F_LOC_ID, temp.meta.data$F_LOC_ID)]) ## perfect - this shows that my GIS data F_LOC_IDs (which were created using a mix of last years (2021) cross ref data and some of this years data from the three sites that had chnaged fields), can be used to look up Unif_codes in an idependant dataset (the field data export the Greg downloaded) via F_LOC_IDs. In summary, this shows that my updated 2022 field-level GIS file should now contain the correct F_LOC_IDs and S_LOC_IDs 



## Finish - new GIS saved here: C:\Data\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\RSPB\BWWM\BWWM project handover\Maps\Shapefiles\2022 BWWM\BWWM_Fields_19Aug2022_with_FLOCIDs





#### STEP 2: Read in and process 2022 avian data ####

################# Luke, you will probably need to repeat this step for EACH survey year, 2021, 2022 and 2023


bwwm_fields_visit <- read.csv("../../Data/BWWM 2022/Original data exports from BTO 2022/bwwm-visits-field@13-09-2022.csv") # visit level field habitat data 
bwwm_obs <- read.csv("../../Data/BWWM 2022/Original data exports from BTO 2022/bwwm-obs@13-09-2022.csv") # wader and other bird spp observations

#bwwm_AES <- read.csv("BWWM2021_AES_Options_per_field_#3.csv") # IGNORE AES DATA, which I will get seperately from Gavin
bwwm_field_meta <- shapefile("../../Maps/Shapefiles/2022 BWWM/BWWM_Fields_19Aug2022_with_FLOCIDs.shp") ## Shapefile - the BWWM field shapefile, with F_LOC_IDs and S_LOC_IDs added 

### Get field area from GIS

for (i in 1:nrow(bwwm_field_meta)){
  bwwm_field_meta@data$Area_RH[i] <- gArea(bwwm_field_meta[i,]) / 10000}
bwwm_field_meta@data$Area_RH <- round(as.numeric(bwwm_field_meta@data$Area_RH), 5)

bwwm_field_meta <- bwwm_field_meta@data # convert shapefile to a dataframe (to save disk space)


#foreign::read.dbf("C:\\Users\\Lucymason\\Documents\\LRM\\BWWM survey 2020\\Maps\\BWWM_Fields_13Oct2021_with_FFLOCID.dbf") # data from GIS table ## Lucy's old code for reading in the 2021 shapefiles



# investigate species records 
# with(bwwm_obs, aggregate(TOTAL_ADULTS, list(ENGLISH_NAME), sum))

# subset all wader species records

ALLwaders <- list("Avocet", "Bar-tailed Godwit", "Black-tailed Godwit", "Black-winged Stilt", "Common Sandpiper", "Curlew", "Dunlin", "Golden Plover", "Great Knot", "Great Snipe", "Green Sandpiper", "Greenshank", "Grey Phalarope", "Grey Plover", "Jack Snipe", "Knot", "Lapwing", "Little Ringed Plover", "Oystercatcher", "Pacific Golden Plover", "Pectoral Sandpiper", "Redshank", "Ruff", "Snipe", "Spotted Redshank", "Stone-curlew", "Temminck's Stint", "Whimbrel", "Wood Sandpiper", "Woodcock")

bwwm_obs_ALLWaders <- bwwm_obs[!is.na(match(bwwm_obs$ENGLISH_NAME, ALLwaders)),] 

# subset all 4 KEY WADER records
KEYwaders <- list("Lapwing", "Redshank", "Curlew", "Snipe")

bwwm_obs_KEYWaders <- bwwm_obs[!is.na(match(bwwm_obs$ENGLISH_NAME, KEYwaders)),] 

# There is enough data to look at Oystercatcher as well if we wanted to
KEYwadersOYC <-list("Lapwing", "Redshank", "Curlew", "Snipe", "Oystercatcher")

# THERE IS A ROW PER SPECIES PER FIELD-VISIT

bwwm_obs_KEYWadersOYC <- bwwm_obs[!is.na(match(bwwm_obs$ENGLISH_NAME, KEYwadersOYC)),] 


# merge field level metadata, habitat and wader observation data into one visit-level dataset
# retain OYC for now

head(bwwm_fields_visit)

field_wader_dat <- merge(bwwm_obs_KEYWadersOYC, bwwm_fields_visit[,-c(1,2,4,5,6,9)], by = c("F_LOC_ID", "VISIT"))  # merges without duplicate columns

write.csv(field_wader_dat, "field_wader_data.csv")


table(field_wader_dat$CBC_CODE)

# CU  L.  OC  RK  SN 
# 47 270 143  85  57

# There are 270 field X Visits with lapwing records => can we assume that where there are lapwing but no other species recorded, the count for the other species = 0? YES I think we need to assume this => so can add in zero counts for the other species for these fields


############################## How many fields held each of the five species


length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "L.",]$F_LOC_ID)) # 167
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "RK",]$F_LOC_ID)) # 58
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "OC",]$F_LOC_ID)) # 93
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "CU",]$F_LOC_ID)) # 32
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "SN",]$F_LOC_ID)) # 49





## BREAK

#field_wader_dat <-  read.csv("field_wader_data.csv")


# merge in data from GIS dbf file: field area, previous survey coverage, reserve info and some old AES info

# match("Area_ha",names(bwwm_field_meta))
# match("Cov_2002",names(bwwm_field_meta))
# match("Cov_1982",names(bwwm_field_meta))
# match("Cov_200910",names(bwwm_field_meta))
# match("RSPB_reser",names(bwwm_field_meta))
# match("NNR",names(bwwm_field_meta))
# match("LNR",names(bwwm_field_meta))
# match("oldCSS_sta",names(bwwm_field_meta))
# match("oldESA_sta",names(bwwm_field_meta))
# match("oldESA",names(bwwm_field_meta))

# match(c("Area_ha","Cov_2002","Cov_1982","Cov_200910","RSPB_reser","NNR","LNR","oldCSS_sta","oldESA_sta","oldESA"),names(bwwm_field_meta))

# field_wader_dat <- merge(field_wader_dat,
# 	bwwm_field_meta[, match(c("F_LOC_ID","Area_ha","Cov_2002","Cov_1982","Cov_200910","RSPB_reser","NNR","LNR","oldCSS_sta","oldESA_sta","oldESA"), names(bwwm_field_meta))], 
# 		by = c("F_LOC_ID"))  # merges with only necessary columns

# 654 rows are lost or not copied across from the field metadata here: which ones are missing? This occurs whether F_LOC_ID or UniF_Code are used as the merge column

# field_wader_dat <- merge(field_wader_dat,
# bwwm_field_meta[, match(c("UniF_Code","Area_ha","Cov_2002","Cov_1982","Cov_200910","RSPB_reser","NNR","LNR","oldCSS_sta","oldESA_sta","oldESA"), names(bwwm_field_meta))],
# 	by = c("UniF_Code"))  # merges with only necessary columns

# bwwm_field_meta$matchvar <- bwwm_field_meta$F_LOC_ID %in% field_wader_dat$F_LOC_ID
# length(bwwm_field_meta[bwwm_field_meta$matchvar == "TRUE", "F_LOC_ID"]) # match in 3419 cases (fields)
# 
# length(unique(field_wader_dat$F_LOC_ID)) # 3608 fields needed - where are the missing fields??

# closer inspection shows that In the field data (exports of wader data and habitat data) some of the field identifiers in UniF_Code have a zero in front of them, e.g. 0870_C_3. In the GIS layer BWWM_Fields_13Oct2021_with_SFLOCID the UniF_Code has NO ZERO in front of it i.e. 870_C_3. This has caused some issues in matching up the field level data with the metadata from the GIS dataset

# remove the unecessary zeros from the field codes in the field_wader_dat dataset

#field_wader_dat$UniF_Code = stringr::str_remove(field_wader_dat$UniF_Code, "^0")

# merge datasets:

head(bwwm_field_meta)

field_wader_dat <- merge(field_wader_dat,
                         bwwm_field_meta[, match(c("F_LOC_ID", "Uniq_code_","Area_RH","Cov_2002","Cov_1982","Cov_200910","RSPB_reser","NNR","LNR","oldCSS_sta","oldESA_sta","oldESA", "Field_id"), names(bwwm_field_meta))],
                         by = c("F_LOC_ID"))  # merges with only necessary columns




# Field_id is the field code from 2009/10 survey - though NOTE, some of the blank cells probably SHOULD have a field_ID. Check this later in the script

# No need to merge in AES data. I will do this later with the bespoke export that Gavin has sent me.

MASTER_field_wader_dat <- field_wader_dat 

write.csv(MASTER_field_wader_dat, "field_visit_wader_dat_MASTER.csv") 

# bwwm_AES$matchvar <- bwwm_AES$UniF_Code %in% field_wader_dat$UniF_Code
# length(bwwm_AES[bwwm_AES$matchvar == "TRUE", "UniF_Code"]) # match in 1451 cases (fields)
# length(bwwm_AES[bwwm_AES$matchvar == "FALSE", "UniF_Code"])
# 
# length(unique(field_wader_dat$UniF_Code)) # 3608 fields needed
# 
# 
# field_wader_dat$matchvar <- field_wader_dat$UniF_Code %in% bwwm_AES$UniF_Code
# length(field_wader_dat[field_wader_dat$matchvar == "TRUE", "UniF_Code"]) # match in 3907 cases (fields)
# length(field_wader_dat[field_wader_dat$matchvar == "FALSE", "UniF_Code"]) # no match in 6176 => these add to 10083 (ok)
# 
# bwwm_AES$duplicated <- duplicated(bwwm_AES$UniF_Code)
# # moved to excel => several duplicated UniF_Codes with different AES option data associated with them
# emailed Greg 2nd November


##### AES history info and options present etc all done in Excel, see BWWM MASTER FIELD AES METADATA 




###### FIELD-LEVEL WADER PAIR ESTIMATES #######


#  NEED TO LOAD IN THIS DATA:

field_wader_dat <-  read.csv("field_visit_wader_dat_MASTER.csv")  # NOTE this is visit level data
head(field_wader_dat)

# histograms of counts: any non-breeding flocks

hist(field_wader_dat[field_wader_dat$ENGLISH_NAME == "Lapwing", "TOTAL_ADULTS"])

library(ggplot2)

ggplot(field_wader_dat, aes(x = TOTAL_ADULTS)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(ENGLISH_NAME ~ VISIT)

with(field_wader_dat, aggregate(TOTAL_ADULTS, list(visit = VISIT, species = ENGLISH_NAME), max))

# visit       species  x
# 1      1        Curlew  8
# 2      2        Curlew  3
# 3      3        Curlew  4
# 4      1       Lapwing 55
# 5      2       Lapwing 43
# 6      3       Lapwing 55
# 7      4       Lapwing  8
# 8      1 Oystercatcher  6
# 9      2 Oystercatcher  8
# 10     3 Oystercatcher  7
# 11     4 Oystercatcher  3
# 12     1      Redshank 12
# 13     2      Redshank  5
# 14     3      Redshank 14
# 15     4      Redshank  2
# 16     1         Snipe  8
# 17     2         Snipe  6
# 18     3         Snipe 11
# 19     4         Snipe  4

lattice::xyplot(TOTAL_ADULTS ~ ESTIMATED_PAIRS | VISIT * ENGLISH_NAME,
                data = field_wader_dat[field_wader_dat$ENGLISH_NAME == "Lapwing",],
                layout = c(2,2))


lattice::xyplot(TOTAL_ADULTS ~ ESTIMATED_PAIRS | VISIT * ENGLISH_NAME,
                data = field_wader_dat[field_wader_dat$ENGLISH_NAME == "Redshank",],
                layout = c(2,2))


lattice::xyplot(TOTAL_ADULTS ~ ESTIMATED_PAIRS | VISIT * ENGLISH_NAME,
                data = field_wader_dat[field_wader_dat$ENGLISH_NAME == "Curlew",],
                layout = c(2,2))


lattice::xyplot(TOTAL_ADULTS ~ ESTIMATED_PAIRS | VISIT * ENGLISH_NAME,
                data = field_wader_dat[field_wader_dat$ENGLISH_NAME == "Snipe",],
                layout = c(2,2))


lattice::xyplot(TOTAL_ADULTS ~ ESTIMATED_PAIRS | VISIT * ENGLISH_NAME,
                data = field_wader_dat[field_wader_dat$ENGLISH_NAME == "Oystercatcher",],
                layout = c(2,2))



# for [species] records; for i in [field ID vector]; calculated number of pairs across visits meeting these criteria:

# Lapwing - half the maximum number of individuals counted on any visit between 15 April and 31 May
# Curlew - the number of pairs derived from the following formula ((maximum number of individuals x 0.71) + 0.1)
# Redshank - mean count of individuals between 24 April and 31 May.
# Oystercatcher - half the maximum number of individuals on any visits between 15 April and 31 May
# # Snipe - it wasn't recorded in the data whether adult counts were chipping/drumming or just flushed
# => SO assumed snipe recorded on evening visits (classed as later than 18:00 - majority of the 4th extra visits for snipe fitted within this period) were likely displaying males, while those recorded at any other time could have been displaying or flushed
# => for pair estimate took 
# EITHER the max number of pairs recorded in PAIR COUNT 	(ESTIMATED_PAIRS) on evening visits 1May to 24June
# OR if Snipe were seen in May or June, the number of individuals/2 (FROM TOTAL_ADULTS)
# (whichever was higher)


#  "CU" "L." "OC" "RK" "SN"

# convert visit dates to R readable format
field_wader_dat$rDATE <- as.Date(field_wader_dat$VISIT_DATE, "%d/%m/%Y")
# [originally in the data download some of the dates were recorded as 2020 or 1921 => these corrected to 2021 in the raw data files prior to R import]

################ There area some 2021 records in the dataset, investigate this

field_wader_dat$YEAR <- substr(field_wader_dat$rDATE, 1,4) 
table(field_wader_dat$YEAR) ## 48 2021 records.... (vs 554 2022 records)

unique(field_wader_dat[field_wader_dat$YEAR == "2021",]$S_LOC_ID)

######## On closer inspection it seems that the 2021 data is a typo. This is because the TWO sites (and their fields) with 2021 data also have 2022 data.  

######## Convert the 2021 data to 2022

field_wader_dat[field_wader_dat$YEAR == "2021",]$rDATE <- paste("2022",substr(field_wader_dat[field_wader_dat$YEAR == "2021",]$rDATE, 5,10),sep="")


###### Correct pair counts for each visit based on proportion of the field surveyed
# Correction of pair calculation for fields NOT 100% covered
# - for each visit need to correct the # individuals based on % of field covered, so if 2 lapwing recorded for a visit where 50% of field covered, # lapwing to use for that visit will be DOUBLE
# also need to do the same for estimated pairs (used for snipe calculations)
# 	- calculation would be x / proportion surveyed
# 	- then use THESE CORRECTED VALUES 
# 	- [then the density calculation would be based on the full field area, i.e. no correction needed at that stage]
# will need to note that this makes an assumption that the # individuals recorded (and the habitat surveyed) was representative of the whole field

field_wader_dat$PropVisCov <- field_wader_dat$COVERED_PERCENT / 100

######## One field claims to have 0% coverage, but on the other two visits it had 100% coverage. Change this probable typo 0 to 100% 

field_wader_dat$PropVisCov[441] <- 1


field_wader_dat$CORRECTED_ADULTS <- with(field_wader_dat, TOTAL_ADULTS / PropVisCov)

field_wader_dat$CORRECTED_PAIRS <- with(field_wader_dat, ESTIMATED_PAIRS / PropVisCov)

#field_wader_dat[field_wader_dat$F_LOC_ID == "LOC3379636",]


# use the CORRECTED ADULT counts to calculate estimated pairs below:

# test <- field_wader_dat[1:20,]
# halfmax <- function(x) max(x)/2
# est_test <- aggregate(test$TOTAL_ADULTS, list(UniF_Code = test$UniF_Code), halfmax)
# names(est_test)[names(est_test) == 'x'] <- "pairs"

### LAPWING ####
halfmax <- function(x) max(x)/2
est_pairsL <- with(field_wader_dat[field_wader_dat$CBC_CODE == "L." & field_wader_dat$rDATE >= "2022-04-15" & field_wader_dat$rDATE <= "2022-05-31",], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))
names(est_pairsL)[names(est_pairsL) == 'x'] <- "Lpairs"

### OYSTERCATCHER ####
halfmax <- function(x) max(x)/2
est_pairsOC <- with(field_wader_dat[field_wader_dat$CBC_CODE == "OC" & field_wader_dat$rDATE >= "2022-04-15" & field_wader_dat$rDATE <= "2022-05-31",], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))
names(est_pairsOC)[names(est_pairsOC) == 'x'] <- "OCpairs"

### CURLEW (on any visit, assumed) ####
cufun <- function(x) (max(x) * 0.71) + 0.1
est_pairsCU <- with(field_wader_dat[field_wader_dat$CBC_CODE == "CU",], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), cufun))
names(est_pairsCU)[names(est_pairsCU) == 'x'] <- "CUpairs"

### REDSHANK ####
# mean count of individuals between 24 April and 31 May
est_pairsRK <- with(field_wader_dat[field_wader_dat$CBC_CODE == "RK" & field_wader_dat$rDATE >= "2022-04-24" & field_wader_dat$rDATE <= "2022-05-31",], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), mean))
names(est_pairsRK)[names(est_pairsRK) == 'x'] <- "RKpairs"


### Snipe calculations ###

# Rob Hawkes, note to self. In late 2022 I adjusted the Snipe pair estimation code due to some errors in the original. The below code for estimating Snipe pairs was extracted from my 'Final pre analysis manips' R script. i.e. this is the code I ended up using in my finla BWWM analyses.

# BEFORE DO THIS NEED TO HAVE:
# - split the visit_date column into date and time separately
# - create a column called "Evening" in the field_wader_dat spreadsheet by hand
# - Evening = Yes for any field visit conducted at or after 18:00 (irrespective of visit number)
# have re-uploaded field_wader_dat from file

######### Calculate which observation were evening observations in R (much easier than Excel)

field_wader_dat$Vis.Time <- substr(field_wader_dat$VISIT_DATE, 12,19) 
field_wader_dat$Vis.Hour <- as.numeric(substr(field_wader_dat$Vis.Time, 1,2)) 
table(field_wader_dat$Vis.Hour)
field_wader_dat$TimeOfDay <- ifelse(field_wader_dat$Vis.Hour >= 18, "Evening", "Day")
evefields <- unique(field_wader_dat[field_wader_dat$TimeOfDay == "Evening", "F_LOC_ID"])  # Just 9 fields with evening visits in 2022




#### proceed with my revised calculations for 2022 ####
## (i.e. counting day counts where evening counts = zero, and counting max adults in the evening instead of the estimated number of pairs)

field_wader_dat$Postive_SN_record <- NA ### flag up data rows that show a positive Snipe record
for(i in 1:nrow(field_wader_dat)){
  field_wader_dat$Postive_SN_record[i] <- ifelse(field_wader_dat$CBC_CODE[i] == "SN" & field_wader_dat$CORRECTED_ADULTS[i] > 0, "Yes","No")}

evefields.with.SN <- unique(field_wader_dat[field_wader_dat$TimeOfDay == "Evening" & field_wader_dat$Postive_SN_record == "Yes", "F_LOC_ID"])  # 5 fields with evening visits also recorded Snipe

est_pairsSNeve <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & field_wader_dat$TimeOfDay == "Evening" & !is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), max))

est_pairsSNother <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))

est_pairsSN_alt <- rbind(est_pairsSNeve, est_pairsSNother) ## NOTE - the zeros in this data set are cases where the observer noted zero for a Snipe evening count (and thus isnt picked up in the 'evefields.with.SN' vector), and no Snipe were seen on a subsequent visit. So, after a lot of manual checking ,I can confirm that these zeros are fine to keep in.

names(est_pairsSN_alt)[names(est_pairsSN_alt) == 'x'] <- "SNpairs_alt"
sum(est_pairsSN_alt$SNpairs_alt) ## 24.30556 in 2022 (for 2021, this should add up to 456.79 Snipe pairs across the 17,053 fields that Lucy Mason originally selected for inclusion in the analysis, but the figure is 530 something for ALL the 2021 fields where SN were observed)

length(est_pairsSN_alt$F_LOC_ID)
length(unique(est_pairsSN_alt$F_LOC_ID))



### Old (incorrect Snipe code, you can ignore this Luke)

# 
# ### SNIPE ####
# 
# # BEFORE DO THIS NEED TO HAVE:
# # - split the visit_date column into date and time separately
# # - create a column called "Evening" in the field_wader_dat spreadsheet by hand
# # - Evening = Yes for any field visit conducted at or after 18:00 (irrespective of visit number)
# # have re-uploaded field_wader_dat from file
# 
# ######### Calculate which observation were evening observations in R (much easier than Excel)
# 
# field_wader_dat$Vis.Time <- substr(field_wader_dat$VISIT_DATE, 12,19) 
# field_wader_dat$Vis.Hour <- as.numeric(substr(field_wader_dat$Vis.Time, 1,2)) 
# table(field_wader_dat$Vis.Hour)
# 
# field_wader_dat$TimeOfDay <- ifelse(field_wader_dat$Vis.Hour >= 18, "Evening", "Day")
# 
# evefields <- unique(field_wader_dat[field_wader_dat$TimeOfDay == "Evening", "F_LOC_ID"])  # Just 9 fields with evening visits in 2022
# 
# 
# #  where evefields match UniF_Code, aggregate by taking the max of the evening survey estimated pairs (where Evening == Yes)
# #  where evefields do not match UniF_Code, aggregate by taking the max individuals / 2
# 
# est_pairsSNeve <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & field_wader_dat$TimeOfDay == "Evening" & !is.na(match(field_wader_dat$F_LOC_ID, evefields)),], aggregate(CORRECTED_PAIRS, list(F_LOC_ID = F_LOC_ID), max))
# 
# est_pairsSNother <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & is.na(match(field_wader_dat$F_LOC_ID, evefields)),], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))
# 
# est_pairsSN <- rbind(est_pairsSNeve, est_pairsSNother)
# names(est_pairsSN)[names(est_pairsSN) == 'x'] <- "SNpairs"
# 
# 
# 
# ##### There is an oddiety with the Snipe calculation. If an evening visit was taken but no Snipe were recorded, then any May/June day visits are dicarded. This means that one field with Snipe during a day visit (but not during the evening visit) has been removed. 
# 
# 
# ###### Alternative Snipe calculations - only considering evening visits when there IS Snipe
# 
# 
# evefields.with.SN <- unique(field_wader_dat[field_wader_dat$TimeOfDay == "Evening" & field_wader_dat$CBC_CODE== "SN", "F_LOC_ID"])  # Just 5 fields with evening visits also recorded Snipe
# 
# 
# est_pairsSNeve <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & field_wader_dat$TimeOfDay == "Evening" & !is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_PAIRS, list(F_LOC_ID = F_LOC_ID), max))
# 
# est_pairsSNother <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-05-01" & field_wader_dat$rDATE <= "2022-06-24" & is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))
# 
# 
# est_pairsSN_alt <- rbind(est_pairsSNeve, est_pairsSNother)
# names(est_pairsSN_alt)[names(est_pairsSN_alt) == 'x'] <- "SNpairs_alt"
# 
# 
# 
# ###### Second alternative Snipe calculation - same as above (only consider evening visits when there IS Snipe) AND extend the lower Snipe observation cut-off to mid April (April 15, as per Lapwing and Oystercatcher). This will be known as the 'lenient' Snipe metric. 
# 
# est_pairsSNeve <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-04-15" & field_wader_dat$rDATE <= "2022-06-24" & field_wader_dat$TimeOfDay == "Evening" & !is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_PAIRS, list(F_LOC_ID = F_LOC_ID), max))
# 
# est_pairsSNother <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2022-04-15" & field_wader_dat$rDATE <= "2022-06-24" & is.na(match(field_wader_dat$F_LOC_ID, evefields.with.SN)),], aggregate(CORRECTED_ADULTS, list(F_LOC_ID = F_LOC_ID), halfmax))
# 
# 
# est_pairsSN_lenient <- rbind(est_pairsSNeve, est_pairsSNother)
# names(est_pairsSN_lenient)[names(est_pairsSN_lenient) == 'x'] <- "SNpairs_lenient"
# 
# 
# 
# 
# 
# 
# 
# # # testing:
# # test <- with(field_wader_dat[field_wader_dat$CBC_CODE == "SN" & field_wader_dat$rDATE >= "2021-05-01" & field_wader_dat$rDATE <= "2021-06-24" ,], aggregate(TOTAL_ADULTS, list(UniF_Code = UniF_Code), halfmax))
# # test[!is.na(match(test$UniF_Code, evefields)),]
# # merge(est_pairsSNeve, test[!is.na(match(test$UniF_Code, evefields)),], by = c("UniF_Code"), all.x = TRUE) 
# 


# MERGE these wader estimates together into one dataset

MergeWaderPAIRS <- function(x, y){
  df <- merge(x, y, by= "F_LOC_ID", all.x= TRUE, all.y= TRUE)
  return(df) }

waderpairs2022 <- Reduce(MergeWaderPAIRS, list(est_pairsL, est_pairsRK, est_pairsCU, est_pairsOC, est_pairsSN_alt)) 


############### Check data
'%!in%' <- function(x,y)!('%in%'(x,y))


### Lapwing

nrow(field_wader_dat[field_wader_dat$CBC_CODE == "L.",]) # 270 field x visits have Lapwing in 2022
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "L.",]$F_LOC_ID)) # 167 fields have Lapwing in 2022
nrow(waderpairs2022[!is.na(waderpairs2022$Lpairs),]) # 150 fields have Lapwing WHEN all the qualifying factors are taken into account (e.g. visit date range)

unique(field_wader_dat[field_wader_dat$CBC_CODE == "L." & field_wader_dat$F_LOC_ID %!in%  waderpairs2022[!is.na(waderpairs2022$Lpairs),]$F_LOC_ID,]$F_LOC_ID) ## this confirms that 17 'missing' fields were not included in the final Lapwing totals because they failed to qualify for the Lapwing inclusion cretria

### Redshank

nrow(field_wader_dat[field_wader_dat$CBC_CODE == "RK",]) # 85 field x visits have Redshank in 2022
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "RK",]$F_LOC_ID)) # 58 fields have Redshank in 2022
nrow(waderpairs2022[!is.na(waderpairs2022$RKpairs),]) # 45 fields have Redshank WHEN all the qualifying factors are taken into account (e.g. visit date range)

unique(field_wader_dat[field_wader_dat$CBC_CODE == "RK" & field_wader_dat$F_LOC_ID %!in%  waderpairs2022[!is.na(waderpairs2022$RKpairs),]$F_LOC_ID,]$F_LOC_ID) ## this confirms that 13 'missing' fields were not included in the final Redshank totals because they failed to qualify for the Lapwing inclusion cretria

### Curlew

nrow(field_wader_dat[field_wader_dat$CBC_CODE == "CU",]) # 47 field x visits have Curlew in 2022
length(unique(field_wader_dat[field_wader_dat$CBC_CODE == "CU",]$F_LOC_ID)) # 32 fields have Curlew in 2022
nrow(waderpairs2022[!is.na(waderpairs2022$CUpairs),]) # 32 fields have Curlew WHEN all the qualifying factors are taken into account (e.g. visit date range)

unique(field_wader_dat[field_wader_dat$CBC_CODE == "CU" & field_wader_dat$F_LOC_ID %!in%  waderpairs2022[!is.na(waderpairs2022$CUpairs),]$F_LOC_ID,]$F_LOC_ID) ## this confirms that 0 'missing' fields were inlcuded in the final Curlew totals

### Snipe


## Add your own checks for Snipe

# There are lots of fields with just one species recorded => can we assume that where there are lapwing but no other species recorded, the count for the other species = 0? YES I think we need to assume this => so can add in zero counts for the other species for these fields
# i.e. where data - NA here this should be 0


waderpairs2022$Lpairs2022 <- ifelse(is.na(waderpairs2022$Lpairs), 0, waderpairs2022$Lpairs)
waderpairs2022$RKpairs2022 <- ifelse(is.na(waderpairs2022$RKpairs), 0, waderpairs2022$RKpairs)
waderpairs2022$OCpairs2022 <- ifelse(is.na(waderpairs2022$OCpairs), 0, waderpairs2022$OCpairs)
waderpairs2022$CUpairs2022 <- ifelse(is.na(waderpairs2022$CUpairs), 0, waderpairs2022$CUpairs)
waderpairs2022$SNpairs2022 <- ifelse(is.na(waderpairs2022$SNpairs), 0, waderpairs2022$SNpairs)
waderpairs2022$SNpairs2022_alt <- ifelse(is.na(waderpairs2022$SNpairs_alt), 0, waderpairs2022$SNpairs_alt)
waderpairs2022$SNpairs2022_lenient <- ifelse(is.na(waderpairs2022$SNpairs_lenient), 0, waderpairs2022$SNpairs_lenient)


write.csv(waderpairs2022, "wader pairs 2022.csv")



#### ADDING ZERO COUNTS TO FIELDS WITH NO WADERS ####
## The above code focuses on those fields that were monitored that held waders. 
## I now need to add in the zeros (waders absent), and these come from two sources. 
## 1) fields that were monitored but held no waders, and 
## 2) fields that were not monitored (the site OR field was deemed unsuitable) and assumed to hold no waders



#### START WITH SOME CHECKS ####

#################### First, estbalish which fields had a 2022 habitat survey, irrespective of wether the field was surveyed or not (this show which sites WHERE visited (i.e. deemed suitable))


bwwm_field_meta$Hab_survey_2022 <- bwwm_field_meta$F_LOC_ID %in% bwwm_fields_visit$F_LOC_ID

### how many sites have at least field with 2022 habitat data?

length(unique(bwwm_field_meta[bwwm_field_meta$Hab_survey_2022 == TRUE,]$S_LOC_ID)) # 277


#################### Next, establish which fields were deemed unsuitable in 2022 (and therefore the field-level wader counts are assumed to be zero)

### Read in Greg's summary of the site level data, LUKE - for you, this file is saved in 'Field data/2022' (for 2022) AND 'Field data/2021' (for 2021). There is no site summary for the 2023 data because every site we took on WAS covered

site.summary.22 <- read.csv("../../Data/BWWM 2022/Original data exports from BTO 2022/BWWM Allocations & Data Entry@13-09-2022.csv")


site.summary.22$covered.22 <- site.summary.22$N_VISIT_DATA > 0 # this sites were 'covered' during 2022


## Check, how many of the covered but unsuitable sites have habitat data?

site.summary.22$unsuit.site <- site.summary.22$UNSUITABLE %in% c("Ac","AD","Bf","Gl","Hd","HX","Ld","LX","NS","Pw","Rd","Sg","SI","Sm","ZZ")

table(site.summary.22[site.summary.22$unsuit.site == TRUE,]$S_LOC_ID %in% bwwm_field_meta[bwwm_field_meta$Hab_survey_2022 == TRUE,]$S_LOC_ID) ## of 240 unsuitable sites, only 8 have have habitat data


## Check, how many of the covered but suitable sites have 2022 habitat data?

site.summary.22$visted.site <- site.summary.22$N_VISIT_DATA > 0 & site.summary.22$UNSUITABLE == ""

site.summary.22$visted.site <- ifelse(is.na(site.summary.22$visted.site) == TRUE, FALSE, site.summary.22$visted.site)

table(site.summary.22[site.summary.22$visted.site == TRUE,]$S_LOC_ID %in% bwwm_field_meta[bwwm_field_meta$Hab_survey_2022 == TRUE,]$S_LOC_ID) ## 269/273 of 'suitable sites' have habitat data


## Explore the 4 'suitable' sites that DONT have habitat data

site.summary.22[site.summary.22$visted.site == TRUE,]$S_LOC_ID[which(site.summary.22[site.summary.22$visted.site == TRUE,]$S_LOC_ID %in% bwwm_field_meta[bwwm_field_meta$Hab_survey_2022 == TRUE,]$S_LOC_ID == FALSE)]

# "LOC2990482" "LOC2991245" "LOC2988645" "LOC3309447"



## IN SUMMARRY - of the the 513 sites that are were 'covered' in 2022, 273 were visted and deemed suitable and 240 were not visited and deemed unsuitable. In total, 277 of these sites have habitat data (269 from sites that were visited and deemed suitable, and 8 from thats that werent visited and deemed 'unsuitable')

######## Final check, are any of the covered sites included in 2021 data (i.e. after I have applied all the field exclusion steps)??


#### Luke, for you, I would skip this step for now and check later one which fields have been monitored in multiple years. There should be minimal duplication between 2021 and 2022, but a fair few 2021/22 fields were resurveyed in 2023. As previously discussed, where this has happened, take the 2023. The only time we might want to make an exception to this rule is for West Sedgemoor - Here we DEFINATELY want to take the 2023 surveys for Snipe, but we may want to take the 2021 surveys for Lapwing. This is because the 2021 West Sedgemoor visit timings in 2021 favour Lapwing, whilst the visit timings in 2023 favour Snipe.

processed.2021.data <- read.csv("../../Data/Covariates/Covariate_R_project/Outputs/Covariate_df/2021_fields_master_df.csv")


nrow(processed.2021.data[processed.2021.data$exclude.field == "No",]) # 13750 fields included in the 2021 data

table(unique(processed.2021.data[processed.2021.data$exclude.field == "No",]$Site_id) %in% na.omit(site.summary.22[site.summary.22$covered.22 == TRUE,]$S_LOC_ID)) ## one site is present in both 2021 and 2022 - investigate.....

# "LOC3309447" is the site in question


site.summary.22[site.summary.22$S_LOC_ID == "LOC3309447",] ## I have no idea why this site was visted in 2021 and 2022


bwwm_fields_visit[bwwm_fields_visit$S_LOC_ID == "LOC3309447",] ## NOTE - this is one of the four sites that were apparently visted in 2022 by contained no habitat (visit) data. Ok to proceed 



#### Add 2022 field coverage infomation ####
## (and site suitability infomation) to the 38,000 + master field dataframe

bwwm_field_meta$Covered_22 <-bwwm_field_meta$S_LOC_ID %in% na.omit(site.summary.22[site.summary.22$covered.22 == TRUE,]$S_LOC_ID)
bwwm_field_meta$Suitable_22 <-bwwm_field_meta$S_LOC_ID %in% na.omit(site.summary.22[site.summary.22$visted.site == TRUE,]$S_LOC_ID)
bwwm_field_meta$Unsuit_22 <-bwwm_field_meta$S_LOC_ID %in% na.omit(site.summary.22[site.summary.22$unsuit.site == TRUE,]$S_LOC_ID)

###### Check that all fields with 2022 wader data are registered as covered (and ideally suitable)  for 2022

table(waderpairs2022$F_LOC_ID %in% bwwm_field_meta[bwwm_field_meta$Suitable_22 == TRUE,]$F_LOC_ID) # yes


length(unique(bwwm_field_meta[bwwm_field_meta$Covered_22 == TRUE,]$S_LOC_ID)) ## 513 sites, which is correct


### Luke, ressume for here



#### Create a new master 2022 dataframe for the occupancy and density analyses ####
## that inlcudes the relevant meta data for all fields with waders and all fields without waders 
## (with both suitable and unsuitable sites inlcuded). 
## Note, this dataset will be subject to further field exclusions in the Covariate project


MASTER_2022_analysis_dataset <- data.frame(F_LOC_ID = bwwm_field_meta[bwwm_field_meta$Covered_22 == TRUE,]$F_LOC_ID, S_LOC_ID = NA, Area_ha = NA, Lpairs2022 = NA, RKpairs2022 = NA, OCpairs2022 = NA, CUpairs2022 = NA, SNpairs2022 = NA, SNpairs2022_alt = NA, SNpairs2022_lenient = NA, Lpa2022 = NA, RKpa2022 = NA, OCpa2022 = NA, CUpa2022 = NA, SNpa2022 = NA, SNpa2022_alt = NA, SNpa2022_lenient = NA, AnySppPA = NA, AnySppPA_SNlenient = NA, RESERVE = "TBC", SSSI = "TBC", GO_REGION = NA, GO_REGION_CODE = NA, REGION_COMBINED = NA, Field.cov.22 = NA, Field.suitable.22 = NA, Field.unsuitable.22 = NA)

length(unique(MASTER_2022_analysis_dataset$F_LOC_ID)) # 5251 fields are currently in the 2022 dataset




#### Add data to the 2022 master dataframe for the occupancy and density analyses ####


MASTER_2022_analysis_dataset$S_LOC_ID <- bwwm_field_meta$S_LOC_ID[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta$F_LOC_ID)] 

length(unique(MASTER_2022_analysis_dataset$S_LOC_ID)) # 513 sites are in the df for 2022

MASTER_2022_analysis_dataset$Area_ha <- bwwm_field_meta$Area_RH[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta$F_LOC_ID)] 

MASTER_2022_analysis_dataset$Lpairs2022<- waderpairs2022$Lpairs2022[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$RKpairs2022<- waderpairs2022$RKpairs2022[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$OCpairs2022<- waderpairs2022$OCpairs2022[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$CUpairs2022<- waderpairs2022$CUpairs2022[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$SNpairs2022<- waderpairs2022$SNpairs2022[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$SNpairs2022_alt<- waderpairs2022$SNpairs2022_alt[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 

MASTER_2022_analysis_dataset$SNpairs2022_lenient<- waderpairs2022$SNpairs2022_lenient[match(MASTER_2022_analysis_dataset$F_LOC_ID, waderpairs2022$F_LOC_ID)] 


############## Convert NA fields to zeros

MASTER_2022_analysis_dataset$Lpairs2022<-ifelse(is.na(MASTER_2022_analysis_dataset$Lpairs2022),0, MASTER_2022_analysis_dataset$Lpairs2022)
MASTER_2022_analysis_dataset$RKpairs2022<-ifelse(is.na(MASTER_2022_analysis_dataset$RKpairs2022),0, MASTER_2022_analysis_dataset$RKpairs2022)
MASTER_2022_analysis_dataset$OCpairs2022<-ifelse(is.na(MASTER_2022_analysis_dataset$OCpairs2022),0, MASTER_2022_analysis_dataset$OCpairs2022)
MASTER_2022_analysis_dataset$CUpairs2022<-ifelse(is.na(MASTER_2022_analysis_dataset$CUpairs2022),0, MASTER_2022_analysis_dataset$CUpairs2022)
MASTER_2022_analysis_dataset$SNpairs2022<-ifelse(is.na(MASTER_2022_analysis_dataset$SNpairs2022),0, MASTER_2022_analysis_dataset$SNpairs2022)
MASTER_2022_analysis_dataset$SNpairs2022_alt<-ifelse(is.na(MASTER_2022_analysis_dataset$SNpairs2022_alt),0, MASTER_2022_analysis_dataset$SNpairs2022_alt)
MASTER_2022_analysis_dataset$SNpairs2022_lenient<-ifelse(is.na(MASTER_2022_analysis_dataset$SNpairs2022_lenient),0, MASTER_2022_analysis_dataset$SNpairs2022_lenient)

############## Species-specific - presence and absence 

MASTER_2022_analysis_dataset$Lpa2022<-ifelse(MASTER_2022_analysis_dataset$Lpairs2022 >0, 1,0)
MASTER_2022_analysis_dataset$RKpa2022<-ifelse(MASTER_2022_analysis_dataset$RKpairs2022 >0, 1,0)
MASTER_2022_analysis_dataset$OCpa2022<-ifelse(MASTER_2022_analysis_dataset$OCpairs2022 >0, 1,0)
MASTER_2022_analysis_dataset$CUpa2022<-ifelse(MASTER_2022_analysis_dataset$CUpairs2022 >0, 1,0)
MASTER_2022_analysis_dataset$SNpa2022<-ifelse(MASTER_2022_analysis_dataset$SNpairs2022 >0, 1,0)
MASTER_2022_analysis_dataset$SNpa2022_alt<-ifelse(MASTER_2022_analysis_dataset$SNpairs2022_alt >0, 1,0)
MASTER_2022_analysis_dataset$SNpa2022_lenient<-ifelse(MASTER_2022_analysis_dataset$SNpairs2022_lenien >0, 1,0)

############# Any wader - presence and absence

MASTER_2022_analysis_dataset$AnySppPA <- MASTER_2022_analysis_dataset$Lpa2022 + MASTER_2022_analysis_dataset$RKpa2022 + MASTER_2022_analysis_dataset$OCpa2022 + MASTER_2022_analysis_dataset$CUpa2022 + MASTER_2022_analysis_dataset$SNpa2022

MASTER_2022_analysis_dataset$AnySppPA_SNlenient <- MASTER_2022_analysis_dataset$Lpa2022 + MASTER_2022_analysis_dataset$RKpa2022 + MASTER_2022_analysis_dataset$OCpa2022 + MASTER_2022_analysis_dataset$CUpa2022 + MASTER_2022_analysis_dataset$SNpa2022_lenient

MASTER_2022_analysis_dataset$AnySppPA <- ifelse(MASTER_2022_analysis_dataset$AnySppPA >0,1,0)
MASTER_2022_analysis_dataset$AnySppPA_SNlenient  <- ifelse(MASTER_2022_analysis_dataset$AnySppPA_SNlenient  >0,1,0)



########### region information (less important for you Luke, as you will be flagging and extracting Somerset, Greater Thames, Suffolk Coast, and Broads data in a later step)

### region infomation is not availible on the 2022 'bwwm_field_meta' df (which is based on the 2022 field-level shapefile). Use the 2021 field-level shapefile data to get region codes for the 2022 data

#bwwm_field_meta_2021 <- shapefile("../../Maps/Shapefiles/2021 BWWM/BWWM_Fields_16Nov2021_with_FLOCID_GOR_COV2021_FINAL.shp") ## Shapefile - the 2021 field-level shapefile

#bwwm_field_meta_2021 <- bwwm_field_meta_2021@data # convert shapefile to dataframe to save space

MASTER_2022_analysis_dataset$GO_REGION <- bwwm_field_meta_2021$REGION_2[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta_2021$F_LOC_ID_2)] 

MASTER_2022_analysis_dataset$GO_REGION_CODE <- bwwm_field_meta_2021$REGCODE[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta_2021$F_LOC_ID_2)] 

### add region infomation for 2022 fields, from three sites, that are missing from the 2021 data set

MASTER_2022_analysis_dataset$GO_REGION <- ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2988391", "East of England", ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2991836", "East of England", ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2989266", "Yorkshire and The Humber", MASTER_2022_analysis_dataset$GO_REGION )))
  
MASTER_2022_analysis_dataset$GO_REGION_CODE <- ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2988391", "EA", ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2991836", "EA", ifelse(MASTER_2022_analysis_dataset$S_LOC_ID == "LOC2989266", "YK", MASTER_2022_analysis_dataset$GO_REGION_CODE)))
  

MASTER_2022_analysis_dataset[which(is.na(MASTER_2022_analysis_dataset$GO_REGION)),]$GO_REGION # GOOD, no NAs
MASTER_2022_analysis_dataset[which(is.na(MASTER_2022_analysis_dataset$GO_REGION_CODE)),]$GO_REGION_CODE # GOOD, no NAs



### Add site suitability infomation (unsuitable and suitable sites) - unsuitable sites were not visited, suitable sites were visted

MASTER_2022_analysis_dataset$Field.cov.22 <- bwwm_field_meta$Covered_22[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta$F_LOC_ID)] 
MASTER_2022_analysis_dataset$Field.suitable.22 <- bwwm_field_meta$Suitable_22[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta$F_LOC_ID)] 
MASTER_2022_analysis_dataset$Field.unsuitable.22<- bwwm_field_meta$Unsuit_22[match(MASTER_2022_analysis_dataset$F_LOC_ID, bwwm_field_meta$F_LOC_ID)]   

table(MASTER_2022_analysis_dataset$Field.cov.22)
table(MASTER_2022_analysis_dataset$Field.suitable.22)
table(MASTER_2022_analysis_dataset$Field.unsuitable.22)

nrow(MASTER_2022_analysis_dataset) ## there are currently 5251 fields for 2022 (though only 2927 were actually visted for breeding waders)


################### EXPORT MASTER DATAFRAME FOR THE OCCPANCY AND ABUNDANCE ANALYSIS 

write.csv(MASTER_2022_analysis_dataset, "../../Data/Analysis datasets USE THESE/2022/MASTER_2022_analysis_dataset.csv")



## NOTE - at this point in Lucy's 2021 code she excluded some inapprioiate fields from the dataframe (based on habitat, e.g. arable fields). However, now that I have a more rigerious exclusion pipeline in  '2022 covariate project', there is no need to filter the data at this stage (i.e the xxx fields will be filtered down in the covariate project)






#### ADD 1) 2022 WADER DATA AND 2) REGIONAL DATA ATTRIBUTES ONTO THE YOUR GIS SHAPEFILE ####

#### Import the 2022 field-level GIS afresh (Luke, for you, thsi will be your BWWM shapefile)

bwwm_field_2022_GIS <- shapefile("../../Maps/Shapefiles/2022 BWWM/BWWM_Fields_19Aug2022_with_FLOCIDs.shp")

#### Add abundance data

bwwm_field_2022_GIS@data$Lprs22 <- MASTER_2022_analysis_dataset$Lpairs2022[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]
bwwm_field_2022_GIS@data$RKprs22 <- MASTER_2022_analysis_dataset$RKpairs2022[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]
bwwm_field_2022_GIS@data$CUprs22 <- MASTER_2022_analysis_dataset$CUpairs2022[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]
bwwm_field_2022_GIS@data$OCprs22 <- MASTER_2022_analysis_dataset$OCpairs2022[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]
bwwm_field_2022_GIS@data$SNprs22 <- MASTER_2022_analysis_dataset$SNpairs2022[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]
bwwm_field_2022_GIS@data$SNprs22Len <- MASTER_2022_analysis_dataset$SNpairs2022_lenient[match(bwwm_field_2022_GIS@data$F_LOC_ID, MASTER_2022_analysis_dataset$F_LOC_ID)]

bwwm_field_2022_GIS@data$Lprs22 <- ifelse(is.na(bwwm_field_2022_GIS@data$Lprs22),0,bwwm_field_2022_GIS@data$Lprs22)
bwwm_field_2022_GIS@data$RKprs22 <- ifelse(is.na(bwwm_field_2022_GIS@data$RKprs22),0,bwwm_field_2022_GIS@data$RKprs22)
bwwm_field_2022_GIS@data$CUprs22 <- ifelse(is.na(bwwm_field_2022_GIS@data$CUprs22),0,bwwm_field_2022_GIS@data$CUprs22)
bwwm_field_2022_GIS@data$OCprs22 <- ifelse(is.na(bwwm_field_2022_GIS@data$OCprs22),0,bwwm_field_2022_GIS@data$OCprs22)
bwwm_field_2022_GIS@data$SNprs22 <- ifelse(is.na(bwwm_field_2022_GIS@data$SNprs22),0,bwwm_field_2022_GIS@data$SNprs22)
bwwm_field_2022_GIS@data$SNprs22Len <- ifelse(is.na(bwwm_field_2022_GIS@data$SNprs22Len),0,bwwm_field_2022_GIS@data$SNprs22Len)


#### Add region data

bwwm_field_2022_GIS@data$REG <- bwwm_field_meta_2021$REGION_2[match(bwwm_field_2022_GIS@data$F_LOC_ID, bwwm_field_meta_2021$F_LOC_ID_2)] 

bwwm_field_2022_GIS@data$REGCODE <- bwwm_field_meta_2021$REGCODE[match(bwwm_field_2022_GIS@data$F_LOC_ID, bwwm_field_meta_2021$F_LOC_ID_2)] 

bwwm_field_2022_GIS@data$REG <- ifelse(bwwm_field_2022_GIS@data$S_LOC_ID == "LOC2988391", "East of England", ifelse(bwwm_field_2022_GIS@data$S_LOC_ID  == "LOC2991836", "East of England", ifelse(bwwm_field_2022_GIS@data$S_LOC_ID  == "LOC2989266", "Yorkshire and The Humber", bwwm_field_2022_GIS@data$REG)))

bwwm_field_2022_GIS@data$REGCODE <- ifelse(bwwm_field_2022_GIS@data$S_LOC_ID == "LOC2988391", "EA", ifelse(bwwm_field_2022_GIS@data$S_LOC_ID == "LOC2991836", "EA", ifelse(bwwm_field_2022_GIS@data$S_LOC_ID == "LOC2989266", "YK", bwwm_field_2022_GIS@data$REGCODE)))


bwwm_field_2022_GIS@data[which(is.na(bwwm_field_2022_GIS@data$REG)),]$REG # GOOD, no NAs
bwwm_field_2022_GIS@data[which(is.na(bwwm_field_2022_GIS@data$REGCODE)),]$REGCODE # GOOD, no NAs



#### Export the new shapefile

writeOGR(bwwm_field_2022_GIS, dsn = "../../Maps/Shapefiles/2022 BWWM", layer = 'BWWM_Fields_19Aug2022_with_FLOCIDs_22DAT_REG', driver = "ESRI Shapefile")


bwwm_field_2022_GIS <- NULL


 



#### LUKE - you can ignore all the below code, this is only needed if you ever end up bolting on historical BWWM data ####

