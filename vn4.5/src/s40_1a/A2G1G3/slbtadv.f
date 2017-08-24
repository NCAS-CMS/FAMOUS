C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C Use, duplication or disclosure of this code is subject to the
C restrictions as set forth in the contract.
C
C                Meteorological Office
C                London Road
C                BRACKNELL
C                Berkshire UK
C                RG12 2SZ
C 
C If no contract has been raised with this copy of the code, the use,
C duplication or disclosure of it is strictly prohibited.  Permission
C to do so must first be obtained in writing from the Head of Numerical
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C
C*LL
CLL   SUBROUTINE SLAB_TEMP_ADVECT
CLL   ---------------------------
CLL
CLL   SUBROUTINE TO ADVECT SLAB OCEAN TEMPERATURE
CLL   TEMPERATURE ON ARAKAWA B GRID USING UPSTREAM
CLL   DIFFERENCING FOR SLAB OCEAN MODEL.
CLL
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   WRITTEN BY R.E.CARNELL (05/07/94)
CLL
CLL  MODEL            MODIFICATION HISTORY FROM INSERTION INTO UM 3.4:
CLL VERSION  DATE
CLL   4.0           Added vertical SST advection. R Carnell(J Crossley)
CLL
CLL   ADHERES TO THE STANDARDS OF DOCUMENTATION PAPER 4, VERSION 6.
CLLEND---------------------------------------------------------------
C*L
      subroutine slab_temp_advect(
     & L1,u_field,landmask
     &,icols,jrows,jrowsm1,lglobal
     &,delta_lat,delta_long,timestep,a,dz1
     &,usea,vsea,tmask,opensea
     &,cos_p_latitude,sec_p_latitude
     &,cos_u_latitude,sin_u_latitude
     &,slabtemp,wtsfc,wtbase
     & )
C
      implicit none
C
      integer
     & L1                            ! in points in p field
     &,u_field                       ! in points in u field
     &,icols                         ! in number of columns E-W
     &,jrows                         ! in number of rows N-S
     &,jrowsm1                       ! in number of rows N-S - 1
      logical
     & lglobal                       ! in true for a global model
     &,landmask(icols,jrows)         ! in land-sea mask (p-grid)
     &,opensea(icols,jrows)          ! in true if box is open sea
      real
     & delta_lat                     ! in EW grid spacing in degrees
     &,delta_long                    ! in NS grid spacing in degrees
     &,timestep                      ! in timestep in seconds
     &,a                             ! in radius of earth in metres
     &,dz1                           ! in slab ocean thickness (m)
     &,usea(icols,jrowsm1)           ! in zonal surface current (M/S)
     &,vsea(icols,jrowsm1)           ! in meri surface current (M/S)
     &,tmask(icols,jrows)            !in mask 1.0 for opensea 0.0 la/ice
      real
     & cos_p_latitude(icols,jrows)   ! in cos(latitude) on p grid
     &,cos_u_latitude(icols,jrowsm1) ! in cos(latitude) on uv grid
     &,sin_u_latitude(icols,jrowsm1) ! in sin(latitude) on uv grid
     &,sec_p_latitude(icols,jrowsm1) ! in sec(latitude) on p grid
     &,slabtemp(icols,jrows)         ! inout temperatue of slab ocean C
     &,wtsfc(icols,jrows)            ! inout w * surface slab temp
     &,wtbase(icols,jrows)           ! inout w * base slab temp
C
C Variables local to this subroutine are now defined
C
      EXTERNAL ZONM
C
      real
     & slabt_init(icols,jrows)       ! initial slab temperature
     &,wbase(icols,jrows)            ! vertical velocity at base
     &,wsfc                          ! vertical velocity at surface =0
     &,divuv(icols,jrows)            ! divergence of u and v
     &,um                            ! usea at i,j
     &,ump                           ! usea at i+1,j
     &,vm                            ! vsea at i,j
     &,vmp                           ! vsea at i,j+i
     &,wm                            ! wsea at surface
     &,wmp                           ! wsea at base
     &,tinx                          ! advection across facex
     &,tinxp                         ! advection across facexp
     &,tiny                          ! advection across facey
     &,tinyp                         ! advection across faceyp
     &,tinz                          ! advection across face top
     &,tinzp                         ! advection across face bottom
     &,tinzsum                       ! sum of vertical advection comps
     &,tinzwt                        ! total weight of upwelling points
     &,tiny_pole                     ! advection across facey at pole
     &,tinyp_pole                    ! advection across faceyp at pole
     &,area                          ! grid box area
     &,fractarea                     ! area*0.125
     &,facex                         ! length of west side of box
     &,facexp                        ! length of east side of box
     &,facey                         ! length of north side of box
     &,faceyp                        ! length of south side of box
     &,wtbpwts                       ! wtbase + wtsfc for conservation
     &,p_levels                      ! no of p levels = 1
     &,latitude_step_inverse         ! 1 / latitude increment
     &,longitude_step_inverse        ! 1 / longitude increment
     &,smask(icols,jrows)            ! Sea mask (p-grid) for zonal mean
     &,s_pmass(icols,jrows)          ! dummy wgt for surf(p_grid)
     &,z_slabtemp(jrows)             ! zonal mean slab temp
     &,dtarea                        ! timestep/area
     &,dtdz                          ! timestep/slab depth
     &,cosplat                       ! cos_p_latitude
     &,coslimit                      ! cos(latitude limit for vertadv)
      parameter ( coslimit = 0.642 ) ! limit is 50 degrees
C
      integer
     & icolsm1                       ! icols - 1
     &,icolsm2                       ! icols - 2
     &,jrowsm2                       ! jrows - 2
     &,i,j,i1,i2,ii                  ! loop counts
     &,spts(jrows)                   ! No of sea points/row
C
      LOGICAL
     & lspts(jrows)                  ! Marks rows with sea pts
C*
C start executable code
C
      icolsm1    = icols-1
      icolsm2    = icols-2
      jrowsm2    = jrows-2
      tinyp_pole = 0.0
      tiny_pole  = 0.0
      p_levels   = 1
      tinzsum    = 0.0
      tinzwt     = 0.0
      dtdz       = timestep / dz1
      latitude_step_inverse  = 1. / delta_lat
      longitude_step_inverse = 1. / delta_long
C
C 1. Calculate zonal mean temperature
C
C 1.1 Set up masks for weighted sums
C
      DO j=1,jrows
        DO i=1,icols
          IF (.not. LANDMASK(I,j)) THEN
            SMASK(I,j) = 1.0
           else
            SMASK(I,j) = 0.0
          ENDIF
        END DO
      END DO
C
C 1.2 Calculate no of sea points on row-by-row basis
C     and set logical array to denote active sea rows
C
        DO j=1,jROWS
          SPTS(j) = 0
          DO i=1,icols
            SPTS(j) = SPTS(j) + SMASK(I,j)
          end do
          if (spts(j) .gt.0) then
            LSPTS(j) = .true.
           else
            LSPTS(j) = .false.
          endif
        end do
C
C 1.3 Set dummy weighting for surface variables to one
C
      DO j=1,jrows
        DO i=1,icols
          S_PMASS(i,j) = 1.0
        end do
      end do
C
C 1.4 Mass weighted zonal mean slabtemp on p grid
C
      call zonm(slabtemp,z_slabtemp
     &,smask,s_pmass,lspts,icols,jrows)
C
C 2. Calculate vertical velocity
C
C     with boundary condition wmp=0 at surface
      wsfc = 0.0
C
C 2.1 Divergence of horizontal velocity
C
      call div_calc(usea,vsea,u_field,L1,p_levels,
     &  icols,sec_p_latitude,cos_u_latitude,
     &  latitude_step_inverse,longitude_step_inverse,1,divuv)
C
C 2.2 Vertical velocity = -dz1 * div (u)
C         wbase positive downwards
C
      do j=1,jrowsm2
        do i=1,icols
          wbase(i,j+1) = -1. * dz1 * divuv(i,j+1)
        end do
      end do
C
C 3. Copy initial slab temperature to workspace
C
      do j=1,jrows
        do i=1,icols
          slabt_init(i,j) = slabtemp(i,j)
        end do
      end do
C
C 4. Conservation of vertical advection
C
C  Make vertical advection of slab temperature conservative.
C  As wtbase and wtsfc are calculated add to get total sum,
C  to conserve this sum needs to be zero.
C  Total made zero by adding difference onto upwelling points.
C  This is done using 5 iterations.
C  Vertical advection only applied between latitude limit (50 N/S)
C
C 4.1 Calculate wtbase and wtsfc and total vertical advection
C
      do j=1,jrowsm2
        cosplat = cos_p_latitude(1,j+1)
C
        do i=1,icols
          i1 = i+1
          i2 = i+2
          if (i1.gt.icols) i1 = i1-icols
          if (i2.gt.icols) i2 = i2-icols
C
          if ( tmask(i1,j+1) .gt. 0.1 ) then
C
            wm  = wsfc
            wmp = wbase(i1,j+1)
C
            if ( cosplat .ge. coslimit ) then
C
C Downwelling advection terms
              if (wmp .gt. 0.) then
                wtsfc(i1,j+1)   = -wmp*slabt_init(i1,j+1)
                wtbase(i1,j+1)  = 0.0
                tinzsum         = tinzsum +wtsfc(i1,j+1)*cosplat
              endif
C
C No vertical advection at these points
              if (wmp .eq. 0.) then
                wtbase(i1,j+1)  = 0.0
                wtsfc(i1,j+1)   = 0.0
              endif
C
C Upwelling advection terms
              if (wmp.lt.0.) then
                wtbase(i1,j+1)  = -wmp*z_slabtemp(j+1)
                wtsfc(i1,j+1)   = 0.0
                tinzwt          = tinzwt  + 1.0*cosplat
                tinzsum         = tinzsum  + wtbase(i1,j+1)*cosplat
              endif
C
             else
C These points are polewards of area of vertical advection.
              wtsfc(i1,j+1)     = 0.0
              wtbase(i1,j+1)    = 0.0
            endif
          endif
        end do
      end do
C
C 4.2 tinsum = weighted sum of vertical advection terms
C      want to make this zero
C      so add it to each upwelling point
C      and divide by weights of upwelling points
C
      wtbpwts =  tinzsum / tinzwt
C
C 4.3 Make tinzsum zero using 5 iterations
C
      do ii=1,5
      tinzwt  = 0.0
      tinzsum = 0.0
      cosplat = 0.0
C
      do j=1,jrowsm2
        cosplat = cos_p_latitude(1,j+1)
C
        do i=1,icols
          i1 = i+1
          i2 = i+2
          if (i1.gt.icols) i1 = i1-icols
          if (i2.gt.icols) i2 = i2-icols
C
          if ( tmask(i1,j+1) .gt. 0.1 ) then
            wmp = wbase(i1,j+1)
C
            if ( cosplat .ge. coslimit ) then
C
C Change upwelling points by wtbpwts
C Add upwelling to total and calculate total weight
              if (wmp.lt.0.) then
                wtbase(i1,j+1)  = wtbase(i1,j+1) - wtbpwts
                tinzwt          = tinzwt  + 1.0 * cosplat
                tinzsum         = tinzsum + wtbase(i1,j+1)*cosplat
              endif
C
C Add downwelling to total
              if (wmp.gt.0.) then
                tinzsum         = tinzsum  + wtsfc(i1,j+1)*cosplat
              endif
C
            endif
          endif
        enddo
      enddo
C
C Recalculate difference to add to upweeling points
      wtbpwts = tinzsum  / tinzwt
      enddo
C
C 5. Loop over grid advecting slab temperature
C
      cosplat =0.0
      facex   = ABS(a * delta_lat)
      facexp  = facex
C
C 5.1 All but polar rows
C
      do j=1,jrowsm2
        area      = ABS(a**2 * delta_long * (sin_u_latitude(1,j+1)
     &               -sin_u_latitude(1,j)))
        fractarea = area * 0.125
        dtarea    = timestep / area
        facey     = ABS(a*cos_u_latitude(1,j)*delta_long)
        faceyp    = ABS(a*cos_u_latitude(1,j+1)*delta_long)
        cosplat   = cos_p_latitude(1,j+1)
C
        do i=1,icols
          i1 = i+1
          i2 = i+2
          if (i1.gt.icols) i1 = i1-icols
          if (i2.gt.icols) i2 = i2-icols
C
          if ( tmask(i1,j+1) .gt. 0.1 ) then
C
            um  = usea(i,j)
            vm  = vsea(i,j)
            ump = usea(i1,j)
            vmp = vsea(i,j+1)
            wm  = wsfc
            wmp = wbase(i1,j+1)
C
C Vertical components
            if (wmp .gt. 0.) then
              tinzp = wtsfc(i1,j+1)
            endif
            if (wmp .eq. 0.) then
              tinzp = 0.0
            endif
            if (wmp .lt. 0.) then
              tinzp = wtbase(i1,j+1)
            endif
            tinz  = wm
C
C Horizontal components
            if (um.ge.0.)  tinx  = um*slabt_init(i,j+1)*tmask(i,j+1)
            if (um.lt.0.)  tinx  = um*slabt_init(i1,j+1)*tmask(i1,j+1)
            if (vm.ge.0.)  tiny  =-vm*slabt_init(i1,j+1)*tmask(i1,j+1)
            if (vm.lt.0.)  tiny  =-vm*slabt_init(i1,j)*tmask(i1,j)
            if (ump.ge.0.) tinxp =-ump*slabt_init(i1,j+1)*tmask(i1,j+1)
            if (ump.lt.0.) tinxp =-ump*slabt_init(i2,j+1)*tmask(i2,j+1)
            if (vmp.ge.0.) tinyp = vmp*slabt_init(i1,j+2)*tmask(i1,j+2)
            if (vmp.lt.0.) tinyp = vmp*slabt_init(i1,j+1)*tmask(i1,j+1)
C
C Polar rows
            if (j.eq.1)       tiny_pole = tiny_pole
     &                            -tiny*facey*timestep/fractarea
            if (j.eq.jrowsm2) tinyp_pole = tinyp_pole
     &                            -tinyp*faceyp*timestep/fractarea
C
C Advection equation
            slabtemp(i1,j+1) =  slabtemp(i1,j+1)
     &      + ( tinx*facex + tinxp*facexp
     &      +   tiny*facey + tinyp*faceyp ) * dtarea
     &      + ( tinz + tinzp ) * dtdz
C
          endif
        end do
      end do
C
C 5.2 Special code for polar rows.
C
      tiny_pole  = tiny_pole/icols
      tinyp_pole = tinyp_pole/icols
C
C First row
      do i=1,icols
        i1 = i+1
        i2 = i+2
        if (i1.gt.icols) i1 = i1-icols
        if (i2.gt.icols) i2 = i2-icols
        if (tmask(i1,1).gt.0.1) then
          slabtemp(i1,1) = slabtemp(i1,1) + tiny_pole
        endif
      end do
C Last row
      do i=1,icols
        i1 = i+1
        i2 = i+2
        if (i1.gt.icols) i1 = i1-icols
        if (i2.gt.icols) i2 = i2-icols
        if (tmask(i1,jrows).gt.0.1) then
          slabtemp(i1,jrows) = slabtemp(i1,jrows) + tinyp_pole
        endif
      end do
C
      return
      end
