###############################################################################
## Helper functions
dyn.load("./src/swamCA.so")

swamCA_1t <- function(gridx, gridy, dem, mask, cella,
                      ppt, evap, runoff, baseflow,
                      wse, otot, itot, delt,  
                      mannN=0.05, cellem=50, 
                      tolwd=0.0001, tolslope=0.001) {
  
  simcf = .Fortran("swamca_1t",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   mask = as.integer(mask), 
                   cella = as.double(cella), 
                   ppt = as.double(ppt),
                   evap = as.double(evap),
                   runoff = as.double(runoff),
                   baseflow = as.double(baseflow),
                   wse = as.double(wse),
                   otot = as.double(dem), itot = as.double(dem),
                   dt = as.double(delt), 
                   mannn = as.double(mannN), cellx = as.double(cellem),
                   cellem = as.double(cellem), 
                   tolwd = as.double(tolwd), tolslope = as.double(tolslope))
  return(simcf)
  
}

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

###############################################################################
## Function to lower border next to outlets
modBorder <- function(outlet, dem, mask) {
  oID <- Which(outlet==1, cells=TRUE)
  if (length(oID) > 0) {
    for (i in 1:length(oID)) {
      cen.crds = xyFromCell(mask, x[i])
      ngb.rc = adjacent(mask, x[i], directions = 8)
      ngb.crds = xyFromCell(mask, ngb.rc[,2])
      ngb.vals = extract(mask, ngb.rc[,2])
      ngb.crds <- ngb.crds[which(ngb.vals==1),]
      ngb.rc <- ngb.rc[which(ngb.vals==1),]
      ngb.dist = distCosine(cen.crds, ngb.crds)
      
      bID <- which.min(ngb.dist)
      dem[ngb.rc[bID,2]] <- dem[ngb.rc[bID,2]]*-1
      mask[ngb.rc[bID,2]] <- mask[ngb.rc[bID,2]]*-1
      
    }
    
  }
  return(list(dem=dem,mask=mask))
}