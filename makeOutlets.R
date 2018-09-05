###############################################################################
## Script to find outlets and modify borders
###############################################################################
###############################################################################
## Libraries
library(rgdal)
library(raster)
source("helpers.R")

## Files
dem.r = raster("./Data/af30c.nc")
## Calculate slope (prior to buffer)
dem.r = extend(dem.r, c(1,1), value = 1e6)
mask.r = dem.r == 1e6

###############################################################################
## Estimate cell areas
area.r = area(dem.r)
## Calculate slope
slope.r = terrain(dem.r, "slope", unit="degrees")
## Find outlet
outlet.r = focal(dem.r, w=matrix(1,3,3), binOutlet, pad=TRUE, padValue=1e6)
mb <- modBorder(outlet.r, dem.r, mask.r)
dem2.r <- mb$dem
mask2.r <- mb$mask
###############################################################################

###############################################################################
## Function to find outlet cells
## Defined as cells next to the border mask, with no lower elevation neighbors
binOutlet <- function (x) {
  outlet = FALSE
  if (max(x) == 1e6) { ## Are we next to an edge?
    if (x[5] != 1e6) { ## Is the cell an edge cell?
      if (which.min(x) == 5) {
        outlet = TRUE
      }
    } 
  } 
}

###############################################################################
## Function to convert radians to degrees
rad2deg <- function(x) {
  x*180/pi
}

###############################################################################
## Function to convert radians to degrees
deg2rad <- function(x) {
  x/180*pi
}

###############################################################################
## Function to calculate great circle distances
gcDist <- function(lon1, lat1, lon2, lat2, r=6378) {
  ## Degree conversion
  lon1 <- lon1 * pi/180
  lon2 <- lon2 * pi/180
  lat1 <- lat1 * pi/180
  lat2 <- lat2 * pi/180
  ## Central angle
  ca <- acos((sin(lat1)*sin(lat2)) + 
      (cos(lat1)*cos(lat2) * cos(abs(lon1-lon2))))
  d <- r * ca
  return(d)
}

