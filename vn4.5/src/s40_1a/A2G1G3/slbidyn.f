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
!+ Main cavitating fluid sea ice dynamic subroutine for slab ocean.
!
! Subroutine Interface:
      SUBROUTINE slab_icedyn(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
! Input arguments
     & icols,jrows,jrowsm1,lglobal
     &,delta_lat,delta_long,timestep
     &,amxsouth,amxnorth,aicemin
     &,Pstar_ice_strength,kappa_ice_strength,cdw
     &,tol_ifree,nmax_ifree,weight_ifree,tol_icav,nmax_icav
! Inout arguments
     &,landmask
     &,cos_u_latitude,cos_p_latitude
     &,sec_p_latitude
     &,sin_u_latitude
     &,coriolis
     &,wsx,wsy,ucurrent,vcurrent
     &,aice,hice,hsnow
     &,icy,newice
     &,opensea
     &,uice,vice
! Output arguments
     &,pmax,pressure
     & )

      IMPLICIT NONE
!
! Description:
!   This routine calculates u and v components of ice velocity
! on Arakawa C grid velocity points and advects ice thickness,
! ice fraction, and snow depth using these velocities and upstream
! advection.
!
! Method:
! Various constants are initialised, and real land-sea masks are
! calculated for the P grid and the C grid velocity points. Primary
! variables ( ice fraction, thickness and snow depth) are copied into
! workspace arrays with zero for land points. The ice strength field
! is calculated (this is a function of ice depth and fraction) and the
! internal ice pressure and various work arrays are zeroed.
!   Currents and velocities are interpolated to the p grid using
! service routines and drag coefficients and the coriolis parameter
! are calculated on the p grid.These are then interpolated to the C
! grid where stresses and dynamics coefficients are calculated.
! SLAB_ICEFREEC is called to calculate free drift velocities given a
! zero pressure field and SLAB_ICECAVRX is called to correct these
! velocities for the effects of internal ice stress. Drag coeff.s and
! dynamics coeff.s are then recalculated and SLAB_ICEFREEC is called
! again, this time with the pressure field produced by SLAB_ICECAVRX,
! and SLAB_ICECAVRX is called to produce final velocities. SLAB_ICE_
! _ADVECT is used to advect the primary variables using an upwind
! scheme. Ice fractions (and corresponding snow depths) outside the
! specified max/min are adjusted, and the logicals ICY and NEWICE
! are reset to reflect the altered ice edge.
!
! Current Code Owner: J.F.Thomson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   3.4     6/94     Original code. J.F.Thomson
!   4.0              Correct external statement. J.F.Crossley
!LL  4.4   04/08/97  Add missing ARGOINDX to various argument lists.
!LL                  D. Robinson.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:              P40
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
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
      INTEGER icols                      ! number of columns
      INTEGER jrows                      ! number of rows
      INTEGER jrowsm1                    ! number of rows minus 1
      LOGICAL lglobal                    ! true for global model
      REAL    delta_lat                  ! meridional grid spacing(deg)
      REAL    delta_long                 ! zonal grid spacing (deg)
      REAL    timestep                   ! slab timestep (seconds)
      REAL    amxsouth                   ! min leads (southern hemi.)
      REAL    amxnorth                   ! min leads (northern hemi.)
      REAL    aicemin                    ! min ice fraction
      REAL    Pstar_ice_strength         ! Parameter in ice strength
      REAL    kappa_ice_strength         ! Parameter in ice strength
      REAL    cdw                        ! Quadratic water drag coeff
      REAL    tol_ifree                  ! tolerance in free drift calc
      INTEGER nmax_ifree                 ! max iterations in free drift
      REAL    weight_ifree               ! underrelaxation weight
      REAL    tol_icav                   ! tolerance in cav fluid calc
      INTEGER nmax_icav                  ! max iterations in cav fluid

!   Array  arguments with intent(in):
      LOGICAL landmask(icols,jrows)        ! true=land points
      REAL    cos_u_latitude(icols,jrowsm1)! cos(lat) u grid
      REAL    cos_p_latitude(icols,jrows)  ! cos(lat) p grid
      REAL    sec_p_latitude(icols,jrows)  ! sec(lat) p grid
      REAL    sin_u_latitude(icols,jrowsm1)! sin(lat) u grid
      REAL    coriolis(icols,jrows)        ! 2 omega sin(lat) p grid
      REAL    wsx(icols,jrowsm1)           ! x wind stress (N/m2)
      REAL    wsy(icols,jrowsm1)           ! y wind stress (N/m2)
      REAL    ucurrent(icols,jrowsm1)      ! x current (m/s)
      REAL    vcurrent(icols,jrowsm1)      ! y current (m/s)

!   Array  arguments with intent(InOut):
      REAL    aice(icols,jrows)     ! fractional ice concentration.
      REAL    hice(icols,jrows)     ! depth avged over grid square (m)
      REAL    hsnow(icols,jrows)    ! snow depth over ice fract only (m)
      REAL    uice(icols,jrowsm1)   ! zonal sea ice velocity.
      REAL    vice(icols,jrowsm1)   ! meridional sea ice velocity.
      LOGICAL icy(icols,jrows)      ! true for ocean pts with aice>.001
      LOGICAL newice(icols,jrows)   ! true for pts where ice is forming
      LOGICAL opensea(icols,jrows)  ! true for ocean pts with no ice

!   Array  arguments with intent(out):
      REAL    pmax(icols,jrows)     ! ice strength (max pressure)
      REAL    pressure(icols,jrows) ! internal ice pressure

! Local parameters:
      REAL    zero                  ! 0.0
      PARAMETER ( zero=0.0 )

! Local scalars:
      INTEGER i         ! loop counter
      INTEGER j         ! loop counter
      INTEGER i1        ! loop counter
      INTEGER j1        ! loop counter
      INTEGER im1       ! loop counter
      INTEGER ip1       ! loop counter
      INTEGER jm1       ! loop counter
      INTEGER jp1       ! loop counter
      INTEGER jrowsby2  ! jrows/2
      INTEGER icolsm1   ! number of columns minus 1
      REAL    uv                    ! workspace scalar.
      REAL    cu                    ! workspace scalar.
      REAL    cv                    ! workspace scalar.
      REAL    dlat_rad              ! meridional grid spacing (radians)
      REAL    dlon_rad              ! zonal grid spacing (radians)

! Local dynamic arrays:
      LOGICAL ccalc(icols,jrows)    ! true= C grid calcs needed
      LOGICAL cavrow(jrows)         ! true= C grid calcs needed on row
      LOGICAL icyrow(jrows)         ! true= ice in this row
      LOGICAL ocrow(jrows)          ! true= ocean points on this row
      REAL psi(jrows)               ! constant water turning angle.
      REAL amx(jrows)               ! maximum ice fraction
      REAL hmask(icols,jrows)       ! 1.0 for land 0.0 for sea points.
      REAL umask(icols,jrowsm1)     ! 1.0 for uv land 0.0 for sea.
      REAL umc(icols,jrowsm1)       ! 1.0 for cu land 0.0 for sea.
      REAL vmc(icols,jrowsm1)       ! 1.0 for cv land 0.0 for sea.
      REAL aice_old(icols,jrows)    ! initial ice fraction field.
      REAL hice_old(icols,jrows)    ! initial ice thickness field.
      REAL hsno_old(icols,jrows)    ! initial snow depth field.
      REAL aice_work(icols,jrows)   ! ice fraction (zero over land)
      REAL hice_work(icols,jrows)   ! ice thickness (zero over land)
      REAL hsnow_work(icols,jrows)  ! snow depth (zero over land)
      REAL uiceh(icols,jrows)       ! x vel comp. on P pts
      REAL viceh(icols,jrows)       ! y vel comp. on P pts
      REAL aice_uv(icols,jrowsm1)   ! ice fraction field on B uv pts
      REAL aice_cu(icols,jrowsm1)   ! ice fraction field on C u pts
      REAL aice_cv(icols,jrowsm1)   ! ice fraction field on C v pts
      REAL wsx_cu(icols,jrowsm1)    ! zonal wind stress on C u pts
      REAL wsy_cv(icols,jrowsm1)    ! merid. wind stress on C v pts
      REAL ucurrent_h(icols,jrows)  ! u current on P pts
      REAL vcurrent_h(icols,jrows)  ! v current on P pts
      REAL uw_cu(icols,jrowsm1)     ! u current on C u pts
      REAL uw_cv(icols,jrowsm1)     ! u current on C v pts
      REAL vw_cu(icols,jrowsm1)     ! v current on C u pts
      REAL vw_cv(icols,jrowsm1)     ! v current on C v pts
      REAL cwstar_h(icols,jrows)    ! drag coeff. on P points
!                                     cdw * rho_water * mod(Uw-U)
      REAL cwstar_cu(icols,jrowsm1) ! drag coeff. on C u pts
      REAL cwstar_cv(icols,jrowsm1) ! drag coeff. on C v pts
      REAL x_stress(icols,jrowsm1)  ! zonal stress on ice for C grid
      REAL y_stress(icols,jrowsm1)  ! merid. stress on ice for C grid
      REAL bh(icols,jrows)          ! mf + cwstar * sin(psi) on P pts
     &                              ! mf is mean ice mass * coriolis

      REAL bx(icols,jrowsm1)        ! bh for C u pts
      REAL by(icols,jrowsm1)        ! bh for C v pts
      REAL ax(icols,jrowsm1)        ! cwstar * sin(psi) for C u pts
      REAL ay(icols,jrowsm1)        ! cwstar * sin(psi) for C v pts

! Function & Subroutine calls:
      External TIMER,SLAB_ICEFREEC,SLAB_ICECAVRX,SLAB_ICE_ADVECT
     &,H_TO_CU,H_TO_CV,UV_TO_H,UV_TO_CU,UV_TO_CV
     &,CU_TO_H,CV_TO_H

!- End of header
! start executable code
! initialise various constants.
      icolsm1  = icols-1
      dlat_rad = delta_lat * pi_over_180
      dlon_rad = delta_long * pi_over_180
! psi is the angle between ice-ocean stress and velocity of ice
! relative to the ocean, assumed constant but different in sign
! between hemispheres.
! psi initialised as 0.4363 i.e. 25 degrees expressed in radians.
      do j=1,jrows
        do i=1,icols
          psi(j) = 0.4363
          if (sin_u_latitude(1,j).lt.0.0) psi(j) = -0.4363
        end do
      end do

! Set up land sea and ice-free sea masks
! Set up amx (max ice fraction as function of latitude)
! Copy ice fraction, thickness and snow depth to workspace
! setting land values to zero.
      do j = 1,jrows
        do i = 1,icols
          hmask(i,j)      =  0.0
          if (.not.landmask(i,j))     hmask(i,j)       = 1.0
          aice_work(i,j)  = aice(i,j)
          if (aice_work(i,j).eq.rmdi) aice_work(i,j)   = zero
          hice_work(i,j)  = hice(i,j)
          if (hice_work(i,j).eq.rmdi) hice_work(i,j)   = zero
          hsnow_work(i,j) = hsnow(i,j)
          if (hsnow_work(i,j).eq.rmdi) hsnow_work(i,j) = zero
        end do
      end do
      jrowsby2 = jrows/2
      do j=1,jrowsby2
        amx(j) = amxnorth
      end do
      do j=jrowsby2+1,jrows
        amx(j) = amxsouth
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

! Ccalc controls where C grid velocity calculations occur.
      do i = 1,icols
        ip1=i+1
        if (ip1.gt.icols) ip1=ip1-icols
        im1=i-1
        if (im1.eq.0) im1=icols
        ccalc(i,1)    =  ( ( icy(im1,1) .or. icy(i,1)
     &                   .or. icy(ip1,1) )
     &                   .and.  ( hmask(i,1) .gt. 0.5 ) )
        ccalc(i,jrows)=  ( ( icy(i,jrowsm1) .or. icy(im1,jrows)
     &                   .or. icy(i,jrows) .or. icy(ip1,jrows) )
     &                   .and.  ( hmask(i,jrows) .gt. 0.5 ) )
      end do
      do j = 2,jrowsm1
        do i = 1,icols
          ip1=i+1
          if (ip1.gt.icols) ip1=ip1-icols
          im1=i-1
          if (im1.eq.0) im1=icols
          ccalc(i,j)    =  ( ( icy(i,j-1) .or. icy(im1,j) .or. icy(i,j)
     &                     .or. icy(ip1,j) .or. icy(ip1,j+1) )
     &                     .and.  ( hmask(i,j) .gt. 0.5 ) )
        end do
        if (lglobal) then
          ccalc(1,j)    = ( ( icy(1,j-1) .or. icy(icols,j) .or. icy(1,j)
     &                     .or. icy(2,j) .or. icy(2,j+1) )
     &                     .and.  ( hmask(1,j) .gt. 0.5 ) )
          ccalc(icols,j)= ((icy(icols,j-1).or.icy(icolsm1,j)
     &                     .or.icy(icols,j).or.icy(1,j).or.icy(1,j+1))
     &                     .and.  ( hmask(icols,j) .gt. 0.5 ) )
        else
          ccalc(1,j)    = ( ( icy(1,j-1) .or. icy(1,j)
     &                     .or. icy(2,j) .or. icy(2,j+1) )
     &                     .and.  ( hmask(1,j) .gt. 0.5 ) )
          ccalc(icols,j)= ( ( icy(icols,j-1).or.icy(icolsm1,j)
     &                     .or.icy(icols,j) )
     &                     .and.  ( hmask(icols,j) .gt. 0.5 ) )
        endif

      end do

! Calculate C grid masks and row masks for icecavrx.
      do j=1,jrowsm1
        do i=1,icolsm1
          umc(i,j)= 1.0
          cu      = hmask(i,j+1) + hmask(i+1,j+1)
          if ( cu .lt. 1.5 ) umc(i,j) = 0.0
          vmc(i,j)= 1.0
          cv      =  hmask(i+1,j) + hmask(i+1,j+1)
          if ( cv .lt. 1.5 ) vmc(i,j) = 0.0
        end do
        if (lglobal) then
          umc(icols,j)= 1.0
          cu      = hmask(icols,j+1) + hmask(1,j+1)
          if ( cu .lt. 1.5 ) umc(icols,j) = 0.0
          vmc(icols,j)= 1.0
          cv      =  hmask(1,j) + hmask(1,j+1)
          if ( cv .lt. 1.5 ) vmc(icols,j) = 0.0
        else
          umc(icols,j) = zero
          vmc(icols,j) = zero
        endif
      end do
      do j = 1,jrows
        ocrow(j)    =  .false.
        icyrow(j)   =  .false.
        cavrow(j)   =  .false.
        do i=1,icols
          if (icy(i,j))   icyrow(j) = .true.
          if (.not.landmask(i,j)) ocrow(j)  = .true.
        end do
      end do
      do j = 1,jrows
        jp1= j+1
        jm1= j-1
        if (j.eq.jrows) jp1=j
        if (j.eq.1)     jm1=j
        cavrow(j)   =  ocrow(j) .and.
     &  (icyrow(jm1).or.icyrow(j).or.icyrow(jp1))
      end do

! Calculate ice strength, pmax and zero pressure field.
      do j=1,jrows
        do i=1,icols
          pmax(i,j) = pstar_ice_strength * hice_work(i,j)
     &                * exp(-kappa_ice_strength*(1-aice_work(i,j)))
          pressure(i,j) = zero
          aice_old(i,j) = aice_work(i,j)
          hice_old(i,j) = hice_work(i,j)
          hsno_old(i,j) = hsnow_work(i,j)
        end do
      end do
      do j=1,jrowsm1
        do i=1,icols
          ucurrent_h(i,j)= zero
          vcurrent_h(i,j)= zero
          uiceh(i,j)     = zero
          viceh(i,j)     = zero
          uw_cu(i,j)     = zero
          uw_cv(i,j)     = zero
          vw_cu(i,j)     = zero
          vw_cv(i,j)     = zero
          wsx_cu(i,j)    = zero
          wsy_cv(i,j)    = zero
        end do
      end do
C
C Interpolate currents and velocities to mass grid.
C
      call uv_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     uice,uiceh,icols,jrows,jrowsm1)
      call uv_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vice,viceh,icols,jrows,jrowsm1)
      call uv_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     ucurrent,ucurrent_h,icols,jrows,jrowsm1)
      call uv_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vcurrent,vcurrent_h,icols,jrows,jrowsm1)
C
C Calculate drag coefficients and coriolis parameter on mass grid.
C
      do j=1,jrows
        do i=1,icols
          cwstar_h(i,j) = rho_water * cdw * sqrt ( (ucurrent_h(i,j)
     &    -uiceh(i,j))**2 + (vcurrent_h(i,j)-viceh(i,j))**2 )
     &    * aice_work(i,j)
          bh(i,j) = ( hice_work(i,j)*rhoice + hsnow_work(i,j)
     &    *aice_work(i,j)*rhosnow ) * coriolis(i,j)
          if ( cwstar_h(i,j) .lt. 0.25 ) cwstar_h(i,j) = 0.25
        end do
      end do
C
C Interpolate drag coeffs and coriolis param. to C grid u and v points
C
      call h_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     cwstar_h,cwstar_cu,jrows,jrowsm1,icols)
      call h_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     cwstar_h,cwstar_cv,jrows,jrowsm1,icols)
      call h_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     bh,bx,jrows,jrowsm1,icols)
      call h_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     bh,by,jrows,jrowsm1,icols)
C
C Interpolate currents from B grid uv points to C grid u and v points.
C ( Also wind stress components. )
C
      call uv_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     ucurrent,uw_cu,jrowsm1,icols)
      call uv_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vcurrent,vw_cu,jrowsm1,icols)
      call uv_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     ucurrent,uw_cv,jrowsm1,icols)
      call uv_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vcurrent,vw_cv,jrowsm1,icols)
      call uv_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     wsx,wsx_cu,jrowsm1,icols)
      call uv_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     wsy,wsy_cv,jrowsm1,icols)
C
C Calculate coeffs ax ay bx by and forcing x_stress y_stress on C grid.
C
      do j=1,jrowsm1
        do i=1,icols
          ax(i,j) = cwstar_cu(i,j)*cos(psi(j))*umc(i,j)
          ay(i,j) = cwstar_cv(i,j)*cos(psi(j))*vmc(i,j)
          bx(i,j) = bx(i,j) + cwstar_cu(i,j)*sin(psi(j))*umc(i,j)
          by(i,j) = by(i,j) + cwstar_cv(i,j)*sin(psi(j))*vmc(i,j)
          x_stress(i,j) = wsx_cu(i,j)
     &              + cwstar_cu(i,j)*(-sin(psi(j))*vw_cu(i,j)
     &              + cos(psi(j))*uw_cu(i,j) )*umc(i,j)
          y_stress(i,j) = wsy_cv(i,j)
     &              + cwstar_cv(i,j)*(sin(psi(j))*uw_cv(i,j)
     &              + cos(psi(j))*vw_cv(i,j) )*vmc(i,j)
        end do
      end do
! Iterative C grid 'free drift' velocity calculations
      call slab_icefreec(
     &              icols,jrows,jrowsm1,lglobal
     &             ,tol_ifree,nmax_ifree,weight_ifree
     &             ,dlat_rad,dlon_rad,cos_p_latitude,umask,ccalc
     &             ,umc,vmc,pressure
     &             ,ax,ay,bx,by,x_stress,y_stress
     &             ,uice,vice
     &             )

! Cavitating Fluid Correction Scheme
      call slab_icecavrx(
     &              icols,jrows,jrowsm1,dlat_rad,dlon_rad
     &             ,tol_icav,nmax_icav
     &             ,cos_p_latitude,cos_u_latitude,sec_p_latitude
     &             ,hmask,umask,umc,vmc,ccalc,cavrow
     &             ,ax,ay,bx,by,pmax
     &             ,uice,vice,pressure
     &             )

! Recalculate drag coefficients, forcing, A and B coefficients using
! mean of initial velocities and corrected velocities.
!
! First interpolate ice velocity to C grid mass points
      call cu_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     uice,uiceh,icols,jrows,jrowsm1)
      call cv_to_h(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vice,viceh,icols,jrows,jrowsm1)

! Calculate drag coefficients and coriolis parameter on mass grid.
      do j=1,jrows
        do i=1,icols
          cwstar_h(i,j) = cwstar_h(i,j)*0.5 +
     &             0.5 * rho_water * cdw * sqrt ( (ucurrent_h(i,j)
     &    -uiceh(i,j))**2 + (vcurrent_h(i,j)-viceh(i,j))**2 )
     &    * aice(i,j)
          if ( cwstar_h(i,j) .lt. 0.25 ) cwstar_h(i,j) = 0.25
        end do
      end do

! Interpolate drag coeffs and coriolis param. to C grid u and v points
      call h_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     cwstar_h,cwstar_cu,jrows,jrowsm1,icols)
      call h_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     cwstar_h,cwstar_cv,jrows,jrowsm1,icols)
      call h_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     bh,bx,jrows,jrowsm1,icols)
      call h_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     bh,by,jrows,jrowsm1,icols)

! Calculate coeffs ax ay bx by and forcing x_stress y_stress on C grid.
      do j=1,jrowsm1
        do i=1,icols
          ax(i,j) = cwstar_cu(i,j)*cos(psi(j))*umc(i,j)
          ay(i,j) = cwstar_cv(i,j)*cos(psi(j))*vmc(i,j)
          bx(i,j) = bx(i,j) + cwstar_cu(i,j)*sin(psi(j))*umc(i,j)
          by(i,j) = by(i,j) + cwstar_cv(i,j)*sin(psi(j))*vmc(i,j)
          x_stress(i,j) = wsx_cu(i,j)
     &              + cwstar_cu(i,j)*(-sin(psi(j))*vw_cu(i,j)
     &              + cos(psi(j))*uw_cu(i,j) )*umc(i,j)
          y_stress(i,j) = wsy_cv(i,j)
     &              + cwstar_cv(i,j)*(sin(psi(j))*uw_cv(i,j)
     &              + cos(psi(j))*vw_cv(i,j) )*vmc(i,j)
        end do
      end do

! Iterative C grid 'free drift' velocity calculations
! Using pressure field from first half timestep.
      call slab_icefreec(
     &              icols,jrows,jrowsm1,lglobal
     &             ,tol_ifree,nmax_ifree,weight_ifree
     &             ,dlat_rad,dlon_rad,cos_p_latitude,umask,ccalc
     &             ,umc,vmc,pressure
     &             ,ax,ay,bx,by,x_stress,y_stress
     &             ,uice,vice
     &             )

! Cavitating Fluid Correction Scheme
      call slab_icecavrx(
     &              icols,jrows,jrowsm1,dlat_rad,dlon_rad
     &             ,tol_icav,nmax_icav
     &             ,cos_p_latitude,cos_u_latitude,sec_p_latitude
     &             ,hmask,umask,umc,vmc,ccalc,cavrow
     &             ,ax,ay,bx,by,pmax
     &             ,uice,vice,pressure
     &             )

! Call ice_advect to advect ice thickness, compactness and snow depth.
      call slab_ice_advect(
     &                icols,jrows,jrowsm1,lglobal
     &               ,dlat_rad,dlon_rad,timestep,a
     &               ,uice,vice,hmask,umask,icy
     &               ,cos_p_latitude,cos_u_latitude,sin_u_latitude
     &               ,aice_work,hice_work,hsnow_work
     &               )

! Adjust ice fractions greater than the max or less than the min.
! Also adjust snow depth accordingly and reset icy and newice.
! And copy work arrays into main arrays.
      do j=1,jrows
        do i=1,icols
          if (aice_work(i,j).gt.amx(j)) then
            hsnow(i,j) = hsnow_work(i,j)*aice_work(i,j)/amx(j)
            aice(i,j)  = amx(j)
          elseif ( (aice_work(i,j).gt.zero)
     &            .and. (.not.landmask(i,j))) then
            if ( aice_work(i,j).lt.aicemin ) then
              hsnow(i,j) = hsnow_work(i,j)*aice_work(i,j)/aicemin
              aice(i,j)  = aicemin
            else
              hsnow(i,j) = hsnow_work(i,j)
              aice(i,j)  = aice_work(i,j)
            endif
          endif
          icy(i,j) = (aice(i,j).gt.zero)
          if (icy(i,j)) then
            newice(i,j)=.false.
            hice(i,j)  = hice_work(i,j)
          endif
        enddo
      enddo

      return
      end
