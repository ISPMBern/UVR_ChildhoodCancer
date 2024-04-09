
#==========================================================================================#
#
# Author: Christian Kreis
#
# Created: 19 October 2021
#
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#
# Import and visualize UV climatological data for exploratory analysis
#
#==========================================================================================#

rm(list=ls(all=TRUE))

library(ggplot2)
library(ggspatial)
library(maptools)
library(rgdal)
library(raster)

Sys.setlocale("LC_TIME", "C")

sessionInfo()

options(scipen=10000)
options(stringsAsFactors=FALSE)

# Working directory
.wd <- "/Volumes/FS/_ISPM/ENVEPI/RESTRICTED"

# File paths
.proj <- "temp/p_UV"

#==========================================================================================#
# Set global parameters and functions
#==========================================================================================#

.t0 <- Sys.time()

.plotting <- "screen"
.plotting <- "PDF"

#==========================================================================================#
# Data Import
#==========================================================================================#

# Shapefiles
#-----------

# Cantons
.cantons <- readOGR(file.path(.wd,"GIS/origdata/PoliticalBoundaries/gg2021/ggg_2021-LV95/shp"), layer="g1k21")

# Communes
#.communes <-  readOGR(file.path(.wd,"GIS/origdata/PoliticalBoundaries/gg2021/ggg_2021-LV95/shp"), layer="g1g21_01012021")

# UV Climatological Data #
#------------------------#

uvr <- readRDS(file.path(.wd,.proj,"R/data/UVR_climatological_data_monthly_11_15.rds"))

#==========================================================================================#
# Data Preparation
#==========================================================================================#

# Shapefiles
#-----------

# Cantons WGS 84
.cantonswgs84 <- spTransform(.cantons,CRS("+init=epsg:4150"))

# Communes WGS 84
#.communeswgs84 <- spTransform(.communes,CRS("+init=epsg:4150"))

# UV Climatological Data #
#------------------------#

# Create stacked raster file of monthly UV climatological data

UVR <- brick(lapply(month.name, function(x) rasterFromXYZ(uvr[uvr$month==x,c("Longitude","Latitude","Mean_UV")],crs=CRS("+init=epsg:4150"),digits=3)))

#==========================================================================================#
# Data Visualization
#==========================================================================================#

if (.plotting=="PDF") {
  pdf(file=file.path(.wd,.proj,paste0("R/graphres/DataExploration/Monthly_UVIndex_11_15.pdf")),
      width=35/2.54, height=16/2.54)
} else {
  quartz(width=35/2.54, height=16/2.54)
}

layout(matrix(1:12,nrow=3,byrow=TRUE))
par(mai=c(0.3,0.4,0.2,0.5), col.axis="grey70")

lapply(1:12, function(i) {
   plot(UVR[[i]],main=month.name[i],las=1)
   plot(.cantonswgs84,add=TRUE)
   }
)

if (.plotting=="PDF") dev.off()

#------------------------------------------------------------------------------------------#

# function to plot UV Index maps

plot.uvi <- function(data) {

ggplot(uvi) + 
   geom_raster(aes(x, y, fill = uvi)) +
   coord_equal() +
# various color schemes
#   scale_fill_gradientn(colours=terrain.colors(100)) +
#   scale_fill_gradientn(colours=plasma(100)) +
  scale_fill_gradientn(colours=rev(rainbow(10,start=0,end=2/3)),name = "UV Index") +
  layer_spatial(.cantonswgs84, fill='transparent',col='grey40') +

# grey scales
   # scale_fill_gradientn(colours=gray.colors(20, start=0, end=1), name = "UV Index") +
   # layer_spatial(.cantonswgs84, fill='transparent',col='black') +

   theme(
     plot.margin=margin(t=-40,r=0,b=-40,l=0),
     axis.title.x=element_blank(),
     axis.title.y=element_blank()
   ) +
   annotation_scale(
      width_hint = 0.13,
      style      = 'bar',
      location   = 'bl',
      pad_x      = unit(0.7, "cm"),
      pad_y      = unit(0.5, "cm"),
      height     = unit(0.15, "cm")
      ) +
   annotation_north_arrow(
      location = 'tr',
      height   = unit(1.0, "cm"),
      width    = unit(1.1, "cm"),
      pad_x    = unit(0.5, "cm"),
      pad_y    = unit(0.7, "cm"),
      style = north_arrow_fancy_orienteering
   )
}

# UV Map July
#------------

uvi <- UVR[[7]]

coords <- xyFromCell(uvi, seq_len(ncell(uvi)))
uvi <- getValues(uvi)
uvi <- as.data.frame(cbind(coords, uvi))

if (.plotting=="PDF") {
  pdf(file=file.path(.wd,.proj,paste0("R/graphres/DisplayElements/Map_UVIndex_11_15_July.pdf")),
      width=16/2.54, height=9.5/2.54)
} else {
  quartz(width=16/2.54, height=9.5/2.54)
}

plot.uvi(uvi)

if (.plotting=="PDF") dev.off()

# UV Map Annual Mean
#-------------------

# Calculate annual mean UVI
uvi <- calc(UVR,fun=mean)

coords <- xyFromCell(uvi, seq_len(ncell(uvi)))
uvi <- getValues(uvi)
uvi <- as.data.frame(cbind(coords, uvi))

if (.plotting=="PDF") {
  pdf(file=file.path(.wd,.proj,paste0("R/graphres/DisplayElements/Map_UVIndex_11_15_Year.pdf")),
      width=16/2.54, height=9.5/2.54)
} else {
  quartz(width=16/2.54, height=9.5/2.54)
}

plot.uvi(uvi)

if (.plotting=="PDF") dev.off()

Sys.time() - .t0

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#




