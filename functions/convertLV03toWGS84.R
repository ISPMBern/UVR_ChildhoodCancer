
#==========================================================================================#
#
# Created: 22 March 2018
#
# Convert LV03 XY-coordinates to WGS84 Longitude and Latitude
#
#==========================================================================================#

# import Swisstopo functions to convert LV03 x,y-coordindates to WGS84 Longitude-Latitude
source(file.path(.wd,.proj,"R/scripts/functions/convertwgs84tolv.R"))

# convert LV03 x,y-coordindates of SNC data points to WGS84 latlon
cl <- makeCluster(detectCores(), methods=FALSE)
clusterExport(cl, c("XY","CH.to.WGS.lng","CH.to.WGS.lat"))
   LonLat <- parLapply(cl, XY, function(i) {temp <- t(apply(i, 1, 
      function(x) c(CH.to.WGS.lng(x[1],x[2]),CH.to.WGS.lat(x[1],x[2]))))
         dimnames(temp)[[2]] <- c("lon","lat")
         return(temp)
      }
   )
stopCluster(cl)

#==========================================================================================#

