###############################################################################
## Version of the SWAM model using hydroCA to move water around
##
## Ref: Guidolin et al. (2016). A weighted cellular automata 2D inundation model 
## for rapid flood analysis Env. Mod. Soft., 84, 378-394
##
###############################################################################
###############################################################################
## Test with hydroshed parameters
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
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]
delt = 60*60 ## Daily integration
cellem = 1000 ## Approximately 1000m cell centers

###############################################################################
## Make up precip grid (mm/day)
pptconst = 20
# ppt.r = setValues(dem.r, 0)
# ppt.r[1,1] <- pptconst
ppt.r = setValues(dem.r, pptconst)
evap.r = setValues(dem.r, 0)
runoff.r = setValues(dem.r, pptconst)
baseflow.r = setValues(dem.r, 0)

###############################################################################
## Grid to record total inflow and outflow from each timestep
itot.r = setValues(dem.r, 0)
otot.r = setValues(dem.r, 0)

###############################################################################
## Estimate cell areas
area.r = area(dem.r)
## Calculate slope
slope.r = terrain(dem.r, "slope", unit="degrees")
## Find outlet
outlet.r = focal(dem.r, w=matrix(1,3,3), binOutlet, pad=TRUE, padValue=1e6)
###############################################################################

###############################################################################
## Convert to matrices
dem = as.matrix(dem.r)
cella = as.matrix(area.r)
mask = as.matrix(mask.r)
ppt = as.matrix(ppt.r)
evap = as.matrix(evap.r)
runoff = as.matrix(runoff.r)
baseflow = as.matrix(baseflow.r)
itot = as.matrix(itot.r)
otot = as.matrix(otot.r)
wse = as.matrix(dem.r)
###############################################################################

for (i in 1:240) {
  print(i)
  sim.out = swamCA_1t(gridx, gridy, dem, mask, cella, 
                      ppt, evap, runoff, baseflow,
                      wse, otot, itot, delt,
                      mannN=0.05, cellem=cellem,
                      tolwd=0.0001, tolslope=0.001)
  wse = sim.out$wse
  print(sum(sim.out$wse-sim.out$dem))
  ## Convert wse to raster for plotting
  wse.r = setValues(dem.r, matrix(sim.out$wse - sim.out$dem, 
                                  nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
  plot(wse.r)
  
}

bclake = readOGR("./Data/929m lake level.kml")
plot(bclake, add=TRUE)
