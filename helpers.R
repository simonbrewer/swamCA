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

