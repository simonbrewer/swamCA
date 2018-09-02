###############################################################################
## Version of the SWAM model using hydroCA to move water around
##
## Ref: Guidolin et al. (2016). A weighted cellular automata 2D inundation model 
## for rapid flood analysis Env. Mod. Soft., 84, 378-394
##
###############################################################################

###############################################################################
## Libraries
library(raster)
dyn.load("./src/swamCA.so")

swamCA_1t <- function(gridx, gridy, dem, mask,
                      ppt, evap, runoff, baseflow,
                      wse, otot, itot, delt,  
                      mannN=0.05, cellem=50, 
                      tolwd=0.0001, tolslope=0.001) {
  
  simcf = .Fortran("swamca_1t",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   mask = as.integer(mask), 
                   ppt = as.double(ppt),
                   evap = as.double(evap),
                   runoff = as.double(runoff),
                   baseflow = as.double(baseflow),
                   wse = as.double(wse),
                   otot = as.double(dem), itot = as.double(dem),
                   dt = as.double(delt), 
                   mannn = as.double(mannN), cellx = as.double(cellem),
                   cellem = as.double(cellem), cella = as.double(cellem*cellem),
                   tolwd = as.double(tolwd), tolslope = as.double(tolslope))
  return(simcf)
  
}

## Files
dem.r = raster("dem1.nc")
dem.r = extend(dem.r, c(1,1), value = 1e6)
mask.r = dem.r == 1e6
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]
delt = 60

###############################################################################
## Make up precip grid (mm/hr)
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
## Convert to matrices
dem = as.matrix(dem.r)
mask = as.matrix(mask.r)
ppt = as.matrix(ppt.r)
evap = as.matrix(evap.r)
runoff = as.matrix(runoff.r)
baseflow = as.matrix(baseflow.r)
itot = as.matrix(itot.r)
otot = as.matrix(otot.r)
wse = as.matrix(dem.r)

for (i in 1:2) {
  print(i)
  sim.out = swamCA_1t(gridx, gridy, dem, mask,
                      ppt, evap, runoff, baseflow,
                      wse, otot, itot, delt,
                      mannN=0.05, cellem=50,
                      tolwd=0.0001, tolslope=0.001)
  wse = sim.out$wse
  print(sum(sim.out$wse-sim.out$dem))
  ## Convert wse to raster for plotting
  wse.r = setValues(dem.r, matrix(sim.out$wse - sim.out$dem, nrow=22, ncol=32))
  plot(wse.r)
  
}

