
#==========================================================================================#
#
# Author: Christian Kreis
#
# Created: 2 November 2021
#
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#
# Extract predictors and prepare data set for analysis
#
#==========================================================================================#

rm(list=ls(all=TRUE))

library(parallel)
library(raster)

Sys.setlocale("LC_TIME", "C")

options(scipen=10000)
options(stringsAsFactors=FALSE)

# Working directory
.wd <- "/Volumes/FS/_ISPM/ENVEPI/RESTRICTED"

# File paths
.dat <- "temp/test_snc_prep/snc/R/data/SNC.TimetoEvent.snc3_90_00_16_new_full_20200824"
.proj <- "temp/p_UV"

#==========================================================================================#
# Set global parameters and functions
#==========================================================================================#

.t0 <- Sys.time()

# data parcelling unit to run data extraction using parallel computing
.n <- 10000

#==========================================================================================#
# Import
#==========================================================================================#

# Time to event data #
#--------------------#

# import CH1903/LV03 coordinates of childhood SNC time to event data set
XY <- read.table(file=file.path(.wd,.dat,"SNC.TimetoEvent.xy.txt"), header=TRUE)

# UV Climatological Data #
#------------------------#

# import monthly 11-15h UVI data
uvr <- readRDS(file.path(.wd,.proj,"R/data/UVR_climatological_data_monthly_11_15.rds"))

#==========================================================================================#
# Preparation
#==========================================================================================#

# Time to event data #
#--------------------#

# split XY-coordinates into subsets for parallel computing
r <- nrow(XY)
if (r>.n) XY <- split(XY,rep(1:ceiling(r/.n),each=.n,length.out=r)) else XY <- list(XY)

#------------------------------------------------------------------------------------------#

# UV Climatological Data #
#------------------------#

# create stacked raster file of monthly climatological UVI data
UVR <- brick(lapply(month.name, function(x) rasterFromXYZ(uvr[uvr$month==x,c("Longitude","Latitude","Mean_UV")],crs=CRS("+init=epsg:4150"),digits=3)))
warnings()

#==========================================================================================#
# Extraction
#==========================================================================================#

# UV Climatological Data #
#------------------------#

# convert LV03 XY-coordinates to WGS84 longitude/latitude
source(file.path(.wd,.proj,"R/scripts/functions/convertLV03toWGS84.R"))

cl <- makeCluster(detectCores(), methods=FALSE)
clusterExport(cl, c("UVR","LonLat"))
   UVI <- parLapply(cl, LonLat, function(x) raster::extract(UVR, x))
stopCluster(cl)

# Unlist the processed data batches
UVI <- do.call(rbind,UVI)

if(nrow(UVI)!=r) warning("Mismatch between UVI extractions and XY data points")

rm(list=ls()[!ls() %in% c("XY","UVI")])

#==========================================================================================#
# Export
#==========================================================================================#

# Combine XY-coordinates with exposure data #
#-------------------------------------------#

# Unlist the XY data batches
XY <- do.call(rbind,XY)

# Merge XY-coordinates withe exposure data
XYUVI <- cbind(XY,UVI)

saveRDS(XYUVI, file.path(.wd,.proj,"R/data/XYUVI_11-15.rds"))

#------------------------------------------------------------------------------------------#

Sys.time() - .t0

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#




