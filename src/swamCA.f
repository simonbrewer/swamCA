!-------------------------------------------------------------------------------
! SWAMCA: Version of the surface water area model that uses a weighted 2D 
!         cellular automata lake water distribution
!
! Includes code and equations from:
!
! Coe (1997). Simulating continental surface waters: an application to 
! Holocene Northern Africa. J. Clim. 10, 1680-1689
!
! Coe (1998). A linked global model of hydrological processes:
! simulation of modern rivers, lakes and wetlands. JGR, 103, 8885-8899
!
! Guidolin et al. (2016). A weighted cellular automata 2D inundation model
! for rapid flood analysis Env. Mod. Soft., 84, 378-394
!
! Ver. 0.1 Takes codes from wca2d.f routine (basic flood model) and
! updates time and linear reservoir code
!-------------------------------------------------------------------------------

      subroutine swamca_1t( m, n, dem,  
     >                      ppt, evap, runoff, baseflow,
     >                      wse, otot, itot, dt, nullcell, 
     >                      mannn, cellx, cellem, cella,
     >                      tolwd, tolslope)
      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      double precision basin( m, n ) ! Elevation values (m)
      double precision ppt( m, n ) ! PPT grid values (mm)
      double precision evap( m, n ) ! PPT grid values (mm)
      double precision runoff( m, n ) ! PPT grid values (mm)
      double precision baseflow( m, n ) ! PPT grid values (mm)
      double precision dt ! Current time step (s)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)
      double precision tolwd ! Minimum water depth for flow to occur
      double precision tolslope ! Minimum slope for flow to occur

      ! Output variables
      double precision wse( m, n ) ! water surface elevation (m)
      double precision fluxes( m, n, 8 ) ! Fluxes to neighboring cells

      ! Internal variables
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)

      !-------------------------------------------------------------------------
      ! Main loop through grid cells

      call calc_flows( m, n, dem, wse, otot, itot,
     >                 fluxes, dt, nullcell, mannn,
     >                 cellx, cellem, cella)

      ! call update_depth( m, n, ppt, wse, itot, otot, dt, cella )

      end

      !-------------------------------------------------------------------------
      ! Subroutine to calculate flows between all cells
      subroutine calc_flows( m, n, dem, wse, otot, itot,
     >                       fluxes, dt, nullcell, mannn,
     >                       cellx, cellem, cella)

      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      double precision dt ! Current time step (s)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)

      ! Output variables
      double precision wse( m, n ) ! PPT grid values (mm)
      double precision fluxes( m, n, 8 ) ! Flux of water between cells

      ! Internal variables
      integer i,j,k
      double precision taut ! Threshold for flow to occur
      integer offx(8),offy(8)
      double precision dl0i(8),dl0i1(8)
      double precision dV0i(8)
      double precision dVmin,dVmax,dVtot
      double precision wi(8) ! Weights
      double precision wi0 ! Weights
      integer maxid(1) ! Weights
      double precision d0,sd0g,nd0g ! Water depth (m)
      double precision vmax ! Maximum velocity (m/s)
      double precision im ! Max. volume to neighbor (m3)
      double precision g ! Gravitational acceleration (m/s-2)
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)
      double precision itotdt

      ! Set up parameters
      parameter(taut = 1e-16)
      parameter(g = 9.81) ! (m/s-2)
      ! offx = (/-1,+1,0,0/) ! Offsets for von Neumann neighborhood
      ! offy = (/0,0,-1,+1/)
      offx = (/-1,0,+1,-1,+1,-1,0,+1/) ! Offsets for Moore neighborhood
      offy = (/-1,-1,-1,0,0,+1,+1,+1/)

      !write(*,*) offx

      !-------------------------------------------------------------------------
      ! Main loop through grid cells

      ! Reset all inflows to zero
      itot(:,:) = 0
      ! Reset all outflows to zero
      otot(:,:) = 0

      do 10 i=2,(m-1)

      do 20 j=2,(n-1)

      if (dem(i,j).ne.nullcell) then ! Check for barrier cells
              ! write(*,*) i,j,dem(i,j),ppt(i,j),itot(i,j),wse(i,j)
      if (wse(i,j).gt.dem(i,j)) then ! Check for standing water in cell
      ! Get neighbors
      do 30 k=1,8
      ! Estimate difference in level
      dl0i(k) = ( wse(i,j) + tau ) - wse((i+offx(k)),(j+offy(k)))
      !write(*,*) 1,k,dl0i(k)

      ! Set negative level differences to 0
      !dl0i(k) = dmax1( myz, dl0i1(k) )
      if (dl0i(k).le.0.0) then
        dl0i(k)=0.0
      endif
      !write(*,*) 2,k,dl0i(k)
      dV0i(k) = max( 0.0, dl0i(k)) * cella ! Equation 2
      !write(*,*) 3,k,dl0i(k),dV0i(k)

30    continue

      ! Volume changes
      dVmax = maxval(dV0i)
      dVmin = minval(dV0i, mask=dV0i.gt.0.0)
      dVtot = sum(dV0i)

      endif ! Standing water check

      endif ! Barrier cell check

20    continue
10    continue

      end


