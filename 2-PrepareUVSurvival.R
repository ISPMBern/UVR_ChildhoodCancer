
#==========================================================================================#
#
# Author: Christian Kreis & Antonella Mazzei
#
# Created: 30 November 2021
#
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#
# Prepare data for survival analysis of UV exposure and childhood hematological cancer risk
#
#==========================================================================================#

rm(list=ls(all=TRUE))

library(dplyr)
library(Hmisc)
library(parallel)
library(raster)
library(survival)

Sys.setlocale("LC_TIME", "english")

options(scipen=10000)
options(stringsAsFactors=FALSE)

# working directory
.wd <- "Y:/ENVEPI/RESTRICTED/temp/p_UV"

# file paths
.dat <- "R/data"

##########################
# Set file paths directory
#.dat  <- "R/data"
#.out  <- "R/graphres"
#.outt <- "R/textres"
#.date <- Sys.Date() 

#==========================================================================================#
# Set global parameters and functions
#==========================================================================================#

.t0 <- Sys.time()

#PARAMETERS TO CHANGE

time.exposure <- "cumulative"
#time.exposure <- "atentry"

lag <- 0

#==========================================================================================#
# Import
#==========================================================================================#

# SNC Childhood Risk Set #
#------------------------#

# import childhood SNC time to event data set
snc <- readRDS(file.path(.wd,.dat,paste0("SNCriskset_cum_exposure_",lag,"y_lag_2021-06-25.rds")))

# UV Climatological Data #
#------------------------#

# import monthly 11-15h UVI data
uvr <- readRDS(file.path(.wd,.dat,"UVR_climatological_data_monthly_11_15.rds"))

#==========================================================================================#
# Preparation
#==========================================================================================#

# create flag for first row of each individual to make calculations easier
snc$first <- as.numeric(NA)
snc$first[snc$age.t1.update==0] <- 1

table(snc$first) #316,837 individuals with lag 0
table(snc$Case[snc$first==1]) #3137 cases with lag 0

#==========================================================================================#
# Create new variables
#==========================================================================================#

## Diagnosis groups ##
 
snc$DxMainGroup <- snc$DxICCC3MainGroup_FINAL
snc$Dx_groups[snc$Case==1] <- snc$DxMainGroup[snc$Case==1]
snc$Dx_groups[snc$Dx_groups >= 4] <- 4
snc$Dx_groups <- factor(snc$Dx_groups, 
                        levels = c(1,2,3,4), 
                        labels = c("Leukemia","Lymphoma","CNS tumors","Other malignant tumor"))

## Categorize hematological malignancies ##

# Main groups
snc$HM       <- as.numeric(snc$Case & snc$DxMainGroup %in% c(1,2))
snc$Leukemia <- as.numeric(snc$Case & snc$DxMainGroup %in% 1)
snc$Lymphoma <- as.numeric(snc$Case & snc$DxMainGroup %in% 2)

# Sub-groups
snc$ALL      <- as.numeric(snc$Case & snc$DxICCC3E %in% 11)
snc$AML      <- as.numeric(snc$Case & snc$DxICCC3E %in% 12)
snc$HL       <- as.numeric(snc$Case & snc$DxICCC3E %in% 21)
snc$NHL      <- as.numeric(snc$Case & snc$DxICCC3E %in% 22)

# categorize year of birth
snc$yob_c <- cut(as.numeric(snc$yob), breaks=seq(1975,2015,5), right=FALSE, include.lowest=TRUE)

## Year of entry at cohort ##
snc$yearEntry <- as.numeric(format(snc$first.gc, '%Y'))
snc$yearEntry.cat <- snc$yearEntry
snc$yearEntry.cat[snc$yearEntry > 2010] <- "2011-2014"
snc$yearEntry.cat <- factor(snc$yearEntry.cat)

# Add BFS major areas (Grossregionen)

snc$region <- as.character(NA)

snc$region[snc$cantonDx %in% c("GE","VD","VS")] <- "1"
snc$region[snc$cantonDx %in% c("BE","FR","JU","NE","SO")] <- "2"
snc$region[snc$cantonDx %in% c("AG","BL","BS")] <- "3"
snc$region[snc$cantonDx %in% c("ZH")] <- "4"
snc$region[snc$cantonDx %in% c("AR","AI","GL","GR","SG","SH","TG")] <- "5"
snc$region[snc$cantonDx %in% c("LU","NW","OW","SZ","UR","ZG")] <- "6"
snc$region[snc$cantonDx %in% c("TI")] <- "7"

snc$region <- factor(snc$region, labels=c("R?gion l?manique","Espace Mittelland",
                     "Northeastern Switzerland","Zurich","Eastern Switzerland",
                     "Central Switzerland","Ticino"))

#==========================================================================================#
# Calculate UV radiation exposure variables
#==========================================================================================#

## Extract UV climatological data for SNC time-to-event data address points ##
#----------------------------------------------------------------------------#

# extract coordinates of SNC time-to-event data set
XY <- snc[,c("geox","geoy")]

# parcelling unit to run data extraction using parallel computing
.n <- 10000

# split XY-coordinates into subsets for parallel computing
r <- nrow(XY)
if (r>.n) XY <- split(XY,rep(1:ceiling(r/.n),each=.n,length.out=r)) else XY <- list(XY)

# create stacked raster file of monthly climatological UVI data
UVR <- brick(lapply(month.name, function(x) rasterFromXYZ(uvr[uvr$month==x,c("Longitude","Latitude","Mean_UV")],crs=CRS("+init=epsg:4150"),digits=3)))
warnings()

# convert LV03 XY-coordinates to WGS84 longitude/latitude
source(file.path(.wd,"R/scripts/functions/convertLV03toWGS84.R"))

cl <- makeCluster(detectCores(), methods=FALSE)
clusterExport(cl, c("UVR","LonLat"))
   UVI <- parLapply(cl, LonLat, function(x) raster::extract(UVR, x))
stopCluster(cl)

# Unlist the processed data batches
UVI <- do.call(rbind,UVI)

# Rename variables
dimnames(UVI)[[2]] <- paste("meanUVI",month.abb,sep=".")

if(nrow(UVI)!=r) warning("Mismatch between UVI extractions and XY data points")

## Prepare UV climatological data for analysis ##
#-----------------------------------------------#

snc <- cbind(snc,UVI)

# calculate annual mean UV Index value
snc$meanUVI <- rowMeans(snc[,names(snc) %in% paste("meanUVI",month.abb,sep=".")], na.rm=TRUE)

# categorize UVI data into quintiles
snc$meanUVI_q <- cut(snc$meanUVI, breaks=quantile(snc$meanUVI[snc$HM %in% 0],
                     probs=seq(0,1,0.2), na.rm=TRUE), include.lowest=TRUE)
snc$meanUVI.Jul_q <- cut(snc$meanUVI.Jul, breaks=quantile(snc$meanUVI.Jul[snc$HM %in% 0],
                     probs=seq(0,1,0.2), na.rm=TRUE), include.lowest=TRUE)

# create variable of UVI quintile-category mean
snc$meanUVI_cat_mean <- as.numeric(by(snc$meanUVI,snc$meanUVI_q,mean))[snc$meanUVI_q]

snc$meanUVI.Jul_cat_mean <- as.numeric(by(snc$meanUVI.Jul,snc$meanUVI.Jul_q,mean))[snc$meanUVI.Jul_q]

#==========================================================================================#
# Calculate ionizing radiation exposure variables
#==========================================================================================#

snc <- snc[order(snc$Group,snc$Case,snc$sncidNUM,snc$age.t1),]
snc$ID.group.record <- factor(paste0(snc$Group,snc$sncidNUM))

# caluclate Cs at entry
snc$Csyear <- snc$yearEntry-1989
snc$CsEntry <- ifelse(snc$Csyear>=-3, snc$Cs137*exp(-0.05*snc$Csyear), snc$Cs137*exp(-0.05*100000))

# calculate total IR at entry
snc$totalIRentry <- rowSums(snc[,c("terrestrial","cosmic","CsEntry")])

# calculate total cum terr, cosmic and Cs separately per individual
snc$totalTerr <- ave(snc$terrestrial.cum, by=snc$ID.group.record, FUN=sum)
snc$totalCosm <- ave(snc$cosmic.cum, by=snc$ID.group.record, FUN=sum)
snc$totalCs   <- ave(snc$Cs.cum, by=snc$ID.group.record, FUN=sum)

#==========================================================================================#
# Format factor variables for tables
#==========================================================================================#

## Sex
snc$sex <- factor(snc$sex, levels = c(0,1), labels = c("Male","Female"))

## Year of birth
snc$yob <- as.numeric(snc$yob)
snc$yob.cat <- cut(snc$yob, breaks = c(1975,1979,1989,1999,2009,2014),
                   labels = c("1975-1979","1980-1989","1990-1999","2000-2009","2010-2014"),include.lowest = T)

## Degree of urbanization
snc$urban.cat <- factor(snc$urban, levels = c(1,2,3), labels = c("Urban","Peri-urban","Rural"))

## Swiss-SEP
snc$ssep_q <- factor(snc$ssep_q, levels = c(1,2,3,4,5), labels = c("1st quintile (low SEP)", "2nd quintile", "3rd quintile", "4th quintile", "5th quintile (high SEP)"))

## NO2
snc$NO2Dx[snc$NO2Dx<0] <- 0

#==========================================================================================#
# Export
#==========================================================================================#

saveRDS(snc, file.path(.wd,.dat,"UVR_CHM_SNCriskset.rds"))

#------------------------------------------------------------------------------------------#

Sys.time() - .t0

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#

