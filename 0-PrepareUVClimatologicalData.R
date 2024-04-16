
#==========================================================================================#
#
# Author: Christian Kreis
#
# Created: 19 October 2021
#
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#
# Prepare UV climatological data for analysis
#
#==========================================================================================#

rm(list=ls(all=TRUE))

Sys.setlocale("LC_TIME", "C")

options(scipen=10000)
options(stringsAsFactors=FALSE)

# Working directory
.wd <- "Y:/ENVEPI/RESTRICTED"

# File paths
.orig <- "origdata/climatological_data"
.proj <- "temp/p_UV"

#==========================================================================================#
# Set global parameters and functions
#==========================================================================================#

.t0 <- Sys.time()

# Global Parameters
#------------------

# Functions
#----------

# Convert WGS84 coordinates to LV95/CH1903+ coordinates
source(file.path(.wd,.proj,"R/scripts/functions/convertwgs84tolv.R"))

#==========================================================================================#
# Data Import
#==========================================================================================#

# UV Climatological Data #
#------------------------#

uvr <- read.csv(file.path(.wd,.proj,.orig,"Bern_mean_month.csv"))
uvr_11_15 <- read.csv(file.path(.wd,.proj,.orig,"Bern_mean_month_11-15.csv"))

#==========================================================================================#
# Data Preparation
#==========================================================================================#

# Convert WGS84 coordinates to LV95/CH1903 coordinates
# X=East-West; Y=North-South; for compatibility with BFS shape files

uvr$X <- WGS.to.CH.y(uvr$Latitude,uvr$Longitude)
uvr$Y <- WGS.to.CH.x(uvr$Latitude,uvr$Longitude)

uvr$X <- uvr$X + 2000000
uvr$Y <- uvr$Y + 1000000

#------------------------------------------------------------------------------------------#

uvr_11_15$X <- WGS.to.CH.y(uvr_11_15$Latitude,uvr_11_15$Longitude)
uvr_11_15$Y <- WGS.to.CH.x(uvr_11_15$Latitude,uvr_11_15$Longitude)

uvr_11_15$X <- uvr_11_15$X + 2000000
uvr_11_15$Y <- uvr_11_15$Y + 1000000

#------------------------------------------------------------------------------------------#

# Add month names

uvr$month <- month.name[uvr$Month]
uvr_11_15$month <- month.name[uvr_11_15$Month]

#==========================================================================================#
# Data Export
#==========================================================================================#

saveRDS(uvr, file=file.path(.wd,.proj,"R/data/UVR_climatological_data_monthly.rds"))
saveRDS(uvr_11_15, file=file.path(.wd,.proj,"R/data/UVR_climatological_data_monthly_11_15.rds"))

#------------------------------------------------------------------------------------------#

Sys.time() - .t0

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#

