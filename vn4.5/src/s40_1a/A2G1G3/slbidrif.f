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
!+ Slab routine controlling simple ice advection scheme
!
! Subroutine Interface:
      SUBROUTINE slab_icedrift(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     & l1,l2,icols,jrows,jrowsm1,landmask
     &,lglobal,aicemin,amxnorth,amxsouth,ah
     &,delta_lat,delta_long,base_lat,timestep
     &,cos_u_latitude,cos_p_latitude
     &,sec_p_latitude
     &,sin_u_latitude
     &,wsx,wsy
     &,ucurrent,vcurrent
     &,aice,hice,hsnow
     &,icy,newice
     &,hinc_diff,hinc_adv,hsinc_adv,areas
     & )

!
! Description:
!   Sea ice dynamics subroutine for slab model, controls
!   simple advection scheme, calling slab_ice_advect to
!   advect the sea ice using surface currents, and slab_diff
!   to diffuse the ice depth using del2 diffusion.
!
! Method:
! This routine first sets up landsea masks and work arrays
! for ice depth, ice fraction and snow depth with zeroes at
! land points for use by the advection routine. U and v current
! components are then interpolated from Arakawa B grid to a C
! grid and zeroed where they would advect ice into a square with
! depth >= 4 metres. This parameterises ice rheology crudely but
! effectively. SLAB_ICE_ADVECT is called to perform upwind advection
! of ice depth, fraction and snow depth on the C grid. The mask is
! then extended so that although ice could advect beyond the previous
! ice edge, it is not allowed to diffuse out from the ice edge.
! Del squared diffusion (based on the horizontal diffusion in ocean
! deck TRACER) is then applied to ice depth. Logical NEWICE is altered
! to account for the extension of ice area by advection.
!
! Current Code Owner: J.F.Thomson
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.4      24/3/94   Original code. J Thomson
! 4.0                Add extra diagnostics and alter criteria for
!                    advection cut-off so it DOES check downstream
!                    ice thickness. J.F.Crossley
!LL  4.4   04/08/97  Add missing ARGOINDX to various argument lists.
!LL                  D. Robinson.
!LL  4.5   03/09/98  Corrected a bug. ucurrent and vcurrent were
!LL                  originally being updated by this subroutine. Now
!LL                  use variables local to this subroutine.
!LL                  C. D. Hewitt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:              P40
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

C*L---------------COMDECK C_SLAB----------------------------------------
C PARAMETERS REQUIRED BY SLAB OCEAN MODEL AND NOT DEFINED ELSEWHERE
C
C CONRATIO IS THE RATIO OF THERMAL CONDUCTIVITIES OF ICE AND SNOW
C (DIMENSIONLESS)
C RHOCP IS THE VOLUMETRIC HEAT CAPACITY OF SEA WATER (J K-1 M-3)
C RHOICE IS THE DENSITY OF ICE (KG M-3)
C RHOSNOW IS THE DENSITY OF SNOW (KG M-3)
C NB ** RHOSNOW is also defined in the common deck C_SOILH, which
C cannot be used in the slab routines as it contains a duplicate
C definition of RHO_WATER, which is also found in C_DENSTY **
C ** It should be noted that the value of RHOSNOW defined here matches
C    the value defined in C_SOIL_H, but differs from that currently
C    used in the ocean GCM (300 Kg m-3)
C
       REAL CONRATIO,RHOCP,RHOICE,RHOSNOW
C
       PARAMETER(CONRATIO=6.5656)
       PARAMETER(RHOCP=4.04E6)
       PARAMETER(RHOICE=900.0)
       PARAMETER(RHOSNOW=250.0)
C
C*----------------------------------------------------------------------

! Subroutine arguments
!========================== COMDECK TYPOINDX ===========================
!
! Description:
!
!       This comdeck contains all the indices and row-wise loop
!       control variables required by the ocean MPP code.
!
! History:
!
!=======================================================================

      INTEGER J_1     ! Local value of loop control for J = 1, n
     &,       J_2     !   "     "    "   "     "     "  J = 2, n
     &,       J_3     !   "     "    "   "     "     "  J = 3, n
     &,       J_JMT   !   "     "    "   "     "     "  J = n, JMT
     &,       J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
     &,       J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
     &,       J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
     &,       JST     ! First row this process considers (no halo)
     &,       JFIN    ! Last   "    "     "     "        "
     &,       J_FROM_LOC
     &,       J_TO_LOC
     &,       JMT_GLOBAL       ! Global value of JMT
     &,       JMTM1_GLOBAL     ! Global value of JMT - 1
     &,       JMTM2_GLOBAL     ! Global value of JMT - 2
     &,       JMTP1_GLOBAL     ! Global value of JMT + 1
     &,       J_OFFSET         ! Global value of JST - 1
     &,       O_MYPE           ! MYPE for ocean arg lists
     &,       O_EW_HALO        ! EW HALO for ocean arg lists
     &,       O_NS_HALO        ! NS HALO for ocean arg lists
     &,       J_PE_JSTM1
     &,       J_PE_JSTM2
     &,       J_PE_JFINP1
     &,       J_PE_JFINP2
     &,       O_NPROC
     &,       imout(4),jmout(4)! i,j indices for pts in Med outflow
     &,       J_PE_IND_MED(4)  ! no for each PE in Med outflow
     &,       NMEDLEV          ! no of levels for Med outflow
     &,      lev_med   ! level at which deep advective Med outflow
c                        exits the Mediterranean
     &,      lev_hud   ! level at which deep advective flow
c                        enters the Hudson Bay
     &,      imout_hud(4)  ! zonal index for Hudson Bay outflow
     &,      jmout_hud(4)  ! merid index for Hudson Bay outflow
     &,      J_PE_IND_HUD(4)  ! PE's involved in HB outflow
     &,      med_topflow   ! last level for which there is inflow to 
C                          ! Mediterranean




!   Scalar arguments with intent(in):
      INTEGER l1             ! length of data
      INTEGER l2             ! length of data to be updated.
      INTEGER icols          ! number of columns.
      INTEGER jrows          ! number of rows.
      INTEGER jrowsm1        ! number of rows minus 1.
      LOGICAL lglobal        ! true if model is global.
      REAL aicemin           ! minimum ice fraction.
      REAL amxnorth          ! maximum ice fraction (arctic).
      REAL amxsouth          ! maximum ice fraction (antarctic).
      REAL delta_lat         ! meridional grid spacing in degrees.
      REAL delta_long        ! zonal grid spacing in degrees.
      REAL base_lat          ! base latitude in degrees.
      REAL timestep          ! slab timestep in seconds.
      REAL ah                ! diffusion coeff for sea ice.

!   Array  arguments with intent(in):
      LOGICAL landmask(icols,jrows)     ! mask true at land points.
      REAL cos_p_latitude(icols,jrows)  ! cos of latitude on p grid.
      REAL cos_u_latitude(icols,jrowsm1)! cos of latitude on u grid.
      REAL sec_p_latitude(icols,jrows)  ! sec of latitude on p grid.
      REAL sin_u_latitude(icols,jrowsm1)! sin of latitude on u grid.
      REAL wsx(icols,jrowsm1)           ! zonal wind stress.
      REAL wsy(icols,jrowsm1)           ! meridional wind stress.
      REAL ucurrent(icols,jrowsm1)      ! zonal surface current.
      REAL vcurrent(icols,jrowsm1)      ! meridional surface current

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      LOGICAL icy(icols,jrows)   ! true ocean points, aice>.001
      LOGICAL newice(icols,jrows)! true points where ice forming.
      REAL aice(icols,jrows)     ! fractional ice concentration.
      REAL hice(icols,jrows)     ! ice depth avg over grid square (m)
      REAL hsnow(icols,jrows)    ! snow depth over ice fract (m)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL uice(icols,jrowsm1)   ! zonal current for advection.
      REAL vice(icols,jrowsm1)   ! meridional current for advectn.
      REAL hinc_diff(icols,jrows)! ice depth inc due to diffusion.
      REAL hinc_adv(icols,jrows) ! ice depth inc due to advection.
      REAL hsinc_adv(icols,jrows)! snow depth inc. * aice(advection).
      REAL areas(icols,jrows)    ! area of grid squares.

! Local parameters:

! Local scalars:
      integer
     & i,j              ! loop counters
     &,icolsm1          ! number of tracer columns minus 1
     &,jrowsby2,jby2p1
      real
     &        zero      ! 0.0
     &,       ahdt      ! ah*timestep
     &,       uv        ! workspace scalar.
     &,       cu        ! workspace scalar.
     &,       cv        ! workspace scalar.
     &,       cumask    ! workspace scalar.
     &,       cvmask    ! workspace scalar.
     &,       dlat_rad  ! grid spacing in radians.
     &,       dlon_rad  ! grid spacing in radians.

! Local dynamic arrays:
      logical
     & ocean(icols,jrows) ! true for non-land points on p grid
     &,dmask(icols,jrows) ! mask for diffusion.
      REAL    amx(jrows)  ! max ice fraction as function of rows.
      real
     & hmask(icols,jrows) ! 1.0 for land 0.0 for sea points.
     &,umask(icols,jrowsm1) ! 1.0 for uv land 0.0 for sea.
      real
     & ucurrent_c(icols,jrowsm1)! u current on C grid h pts
     &,vcurrent_c(icols,jrowsm1)! v current on C grid h pts
     &,ucurrent_l(icols,jrowsm1) ! local copy of u current on P grid
     &,vcurrent_l(icols,jrowsm1) ! local copy of v current on P grid
     &,aice_old(icols,jrows)    ! initial ice fraction
     &,aice_work(icols,jrows)   ! ice fraction (no mdi)
     &,hice_old(icols,jrows)    ! initial ice depth
     &,hice_work(icols,jrows)   ! ice depth (no mdi)
     &,hice_cu(icols,jrowsm1)   ! ice depth on c grid u points
     &,hice_cv(icols,jrowsm1)   ! ice depth on c grid v points
     &,hsnow_old(icols,jrows)   ! initial snow depth
     &,hsnow_work(icols,jrows)  ! snow depth (no mdi)
     &,diffus(icols,jrows)      ! ice depth increments due to diffusion

! Function & Subroutine calls:
      External uv_to_cu,uv_to_cv,slab_ice_advect,slabdiff

!- End of header

! initialise various constants.
      icolsm1 = icols-1
      zero    = 0.000E+00
      ahdt    = ah * timestep
      dlat_rad= delta_lat * pi_over_180
      dlon_rad= delta_long * pi_over_180

! First set up land sea and ice-free sea masks
      do j = 1,jrows
        do i = 1,icols
          ocean(i,j)    =  .not.landmask(i,j)
          dmask(i,j)    =  ocean(i,j)
          hmask(i,j)    =  0.0
          if (ocean(i,j)) hmask(i,j) = 1.0
          aice_work(i,j) = aice(i,j)
          if (aice_work(i,j).eq.rmdi) aice_work(i,j)=zero
          hice_work(i,j) = hice(i,j)
          if (hice_work(i,j).eq.rmdi) hice_work(i,j)=zero
          hsnow_work(i,j) = hsnow(i,j)
          if (hsnow_work(i,j).eq.rmdi) hsnow_work(i,j)=zero
        end do
      end do
      jrowsby2=jrows/2
      jby2p1=jrowsby2+1
      do j=1,jrowsby2
        amx(j)=amxnorth
      end do
      do j=jby2p1,jrows
        amx(j)=amxsouth
      end do

! Calculate Arakawa B grid ice velocity mask.
      do j = 1,jrowsm1
        do i = 1,icolsm1
          umask(i,j)    =  1.0
          uv = hmask(i,j)+hmask(i+1,j)+hmask(i,j+1)+hmask(i+1,j+1)
          if (uv.lt.3.5)  umask(i,j)    =  0.0
        end do
        if (lglobal) then
          umask(icols,j)  =  1.0
          uv = hmask(icols,j)+hmask(icols,j+1)+hmask(1,j)+hmask(1,j+1)
          if (uv.lt.3.5)    umask(icols,j) = 0.0
        else
          umask(icols,j) = 0.0     ! what should i do here ?
        endif
      end do
      do j=1,jrowsm1
        do i=1,icols
          ucurrent_l(i,j) = ucurrent(i,j)*umask(i,j)
          vcurrent_l(i,j) = vcurrent(i,j)*umask(i,j)
        end do
      end do

! Interpolate currents to C grid.
      call uv_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     ucurrent_l,ucurrent_c,jrowsm1,icols)

      call uv_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vcurrent_l,vcurrent_c,jrowsm1,icols)

! Copy initial ice depths and snow depths to workspace
      do j=1,jrows
        do i=1,icols
          aice_old(i,j) = aice(i,j)
          hice_old(i,j) = hice(i,j)
          hsnow_old(i,j) = hsnow(i,j)
          diffus(i,j) = 0.0
        end do
      end do

! Zero currents if ice depth >= 4 metres and ice flowing to thicker area
      do j=1,jrowsm1
        do i=1,icolsm1
          if (ucurrent_c(i,j).gt.0.0
     &    .and.hice(i+1,j+1).gt.4.0)
     &    ucurrent_c(i,j)=0.0
          if (ucurrent_c(i,j).le.0.0
     &    .and.hice(i,j+1).gt.4.0)
     &    ucurrent_c(i,j)=0.0
          if (vcurrent_c(i,j).gt.0.0
     &    .and.hice(i+1,j).gt.4.0)
     &    vcurrent_c(i,j)=0.0
          if (vcurrent_c(i,j).le.0.0
     &    .and.hice(i+1,j+1).gt.4.0)
     &    vcurrent_c(i,j)=0.0
        end do
        if (lglobal) then
          if (ucurrent_c(icols,j) .gt. 0.0 .and.
     &    hice(1,j+1).gt.4.0)
     &    ucurrent_c(icols,j)=0.0
          if (ucurrent_c(icols,j) .lt. 0.0 .and.
     &    hice(icols,j+1).gt.4.0)
     &    ucurrent_c(icols,j)=0.0
          if (vcurrent_c(icols,j) .gt. 0.0 .and.
     &    hice(1,j).gt.4.0)
     &    vcurrent_c(icols,j)=0.0
          if (vcurrent_c(icols,j) .le. 0.0 .and.
     &    hice(1,j+1).gt.4.0)
     &    vcurrent_c(icols,j)=0.0
        else
          ucurrent_c(icols,j)=ucurrent_c(icolsm1,j)
          vcurrent_c(icols,j)=vcurrent_c(icolsm1,j)
        end if
      end do
      do i=1,icols
        vcurrent_c(i,1) = 0.0
      end do

! Call ice_advect to advect ice thickness, compactness and snow depth.
      call slab_ice_advect(
     & icols,jrows,jrowsm1,lglobal,dlat_rad,dlon_rad,timestep,a
     &,ucurrent_c,vcurrent_c,hmask,umask,icy,cos_p_latitude
     &,cos_u_latitude,sin_u_latitude,aice_work,hice_work,hsnow_work
     &,areas
     &          )

! Calculate diffusion increments using slabdiff (as used for slabtemp)
! Extend dmask to be zero over all ice-free areas to prevent diffusion
! out from ice edge and initialise increment array.
      do j=1,jrows
        do i=1,icols
         if (ocean(i,j)) then
           if (aice_work(i,j).eq.zero) dmask(i,j)=.false.
           hinc_diff(i,j)=hice_work(i,j)
           hinc_adv(i,j) =hice_work(i,j)-hice_old(i,j)
           hsinc_adv(i,j)=hsnow_work(i,j)*aice_work(i,j)-hsnow_old(i,j)
     &                    *aice_old(i,j)
          endif
        end do
      end do

! call diffusion subroutine
      call slabdiff ( hice_work
     &,dmask
     &,l1,l2
     &,jrows
     &,icols
     &,ahdt
     &,delta_long,delta_lat,base_lat
     &,cos_p_latitude,cos_u_latitude,sec_p_latitude
     &                )

! Calculate increment in ice depth due to diffusion.
      do j=1,jrows
        do i=1,icols
         if (ocean(i,j)) hinc_diff(i,j)=hice_work(i,j)-hinc_diff(i,j)
        end do
      end do

! Adjust ice fractions greater than the max, or less than the min.
! Also adjust snow depth accordingly and reset icy and newice.
      do j=1,jrows
        do i=1,icols
          if (aice_work(i,j).gt.amx(j)) then
            hsnow_work(i,j) = hsnow_work(i,j)*aice_work(i,j)/amx(j)
            aice_work(i,j)  = amx(j)
          elseif (aice_work(i,j).gt.zero.and.aice_work(i,j).lt.aicemin
     &            .and. ocean(i,j) ) then
            hsnow_work(i,j) = hsnow_work(i,j)*aice_work(i,j)/aicemin
            aice_work(i,j)  = aicemin
          endif
          icy(i,j) = (aice_work(i,j).gt.zero)
        end do
      end do

! Deal with boxes where ice has advected over open ocean.
      do j=1,jrows
        do i=1,icols
          if (icy(i,j)) newice(i,j)=.false.
        end do
      end do

! copy work variables to primary space
      do j=1,jrows
        do i=1,icols
          if (aice_work(i,j).gt.zero) aice(i,j)=aice_work(i,j)
          if (hice_work(i,j).gt.zero) hice(i,j)=hice_work(i,j)
          if (hsnow_work(i,j).gt.zero) hsnow(i,j)=hsnow_work(i,j)
        end do
      end do

      return
      end
