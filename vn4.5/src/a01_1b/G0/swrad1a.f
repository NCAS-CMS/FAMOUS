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
CLL  Subroutine SWRAD --------------------------------------------------
CLL
CLL  Its main function is to gather the input data for daylit points and
CLL  pass them to SWMAST, the top-level routine for P234, the
CLL  plug-compatible interaction of shortwave radiation with the
CLL  atmosphere, and to scatter the output back.  It may return fluxes
CLL  at all layer boundaries, or heating rates produced by differencing
CLL  the fluxes (plus the surface flux); it can also deal with
CLL  shortwave diagnostics.
CLL  Before SWRAD is called, SWLKIN (in deck SWTRAN) must be CALLed to
CLL                                     initialize LUT
CLL  If *DEF IBM is set, the code is standard FORTRAN 77 except for
CLL  having ! comments (it then sets the "vector length" to be 1) but
CLL  otherwise it includes CRAY automatic arrays also.
CLL
CLL   Option CCLD3, set if A01_2A is chosen, combines the
CLL   layer clouds so that at each point the plug-compatible SW only
CLL   sees one layer cloud of each "type" ("high", "medium" and "low")
CLL   - as well as the convective tower, of course.  The boundaries
CLL   between these types are defined in terms of eta, and the model eta
CLL   values passed in are used to convert these to layer numbers.
CLL   These layer cloud amounts reduced to 3 layers are also made
CLL   available as diagnostics.                 William Ingram  8/10/92
CLL  Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.4   31/8/94     nupdate *IFs replaced by FORTRAN IFs  (W Ingram)
CLL   3.4   31/8/94       Compiling system directives added  (W Ingram)
CLL   4.0  1/2/95    SWOUT zeroed at ALL night points for safety.  WJI
!LL   4.0   28/9/95  FOCWWIL COMDECK now subroutine CALL. (A Bushell)
CLL   4.1   17.1.95  NSSSB1 is renamed NSSB1 and is now set over all
CLL                  points instead of over sea surface only; needed
CLL                  over land to derive downward SW in band 1 (TDSSB1)
CLL                  which provides photosynthetically active radiation
CLL                  for vegetation model in Section 3.  This is added
CLL                  to the SWOUT array as an 'extra level', without
CLL                  Zenith Angle adjustment, to enable use in all
CLL                  physics timesteps.             (R A Betts)
CLL                                                       
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled.
CLL                       S.J.Swarbrick
!     4.4   01/7/97  Allow replacement of FOCWWIL parametrization by
!                    direct ratio of prognostic cloud ice to liquid
!                    in _LAYER_ cloud calculations. Note that FOCWWIL is
!                    also used to partition CONVECTIVE cloud and thus is
!                    retained.  (A C Bushell)
CLL
CLL  It technically conforms with standard A of UMDP 4 (version 3,
CLL  07/9/90), but makes many assumptions about STASH structure, and is
CLL  not plug-compatible.
CLL
CLL  It is part of component P233 (ancillary calculations for the
CLL  shortwave radiation), which is in task P23 (radiation).  It also
CLL  performs some of the functions of D23 (radiation diagnostics).
CLL
CLL  Offline documentation (where appropriate) is in UMDP 23.
CLLEND --------------------------------------------------------------
C*L
      SUBROUTINE SWRAD (H2OIN, CO2, O3IN, PSTIN, ABIN, BBIN, LCAIN,
     &     LCW1IN, LCW2IN, CCAIN, CCWPIN, CCBIN, CCTIN, SALIIN, SAOSIN,
     &     AICE, COSZIN, LIT, LAND, LIST, TAC, SCS, LUT, PTS,
     &     OSDIA, OSON, CSOSDI, CSOSON, NSSB1, NSS1ON, TDSS, TDSSON,
     &     CSSSD, CSSSDO, CSSSU,  CSSSUO, LCASW, LCASWO, CCASW,  CCASWO,
     &     LCAAR, LCAARO, LCAARL, LCAARB, LCAAF, LCAAFO, LCAAFL, LCAAFB,
     &     CCAAR, CCAARO, CCAARB, CCAAF, CCAAFO, CCAAFB, TCASW,  TCASWO,
     &     CREFF,CREFFO,LREFF,LREFFO,CVAMT,CVAMTO,LRAMT,LRAMTO,
     &     CWPAJ,CWPAJON,MICRO,
     &     LCA3L, LCA3ON, LCLD3,
     &     L_CLOUD_WATER_PARTITION,
     &     NLIT, NDO, NLEVS, NCLDS, NWET, NOZONE, L1,                   
     &                                           NTSWIN,  SWSEA,  SWOUT)
      EXTERNAL SWMAST, SWDTCA, LSP_FOCWWIL
C*
!     EITHER
!       Use temperature dependent focwwil for convection but calculate
!       ratio in layer cloud from prognostic cloud ice produce as part
!       of large-scale precipitation scheme 3A, OR
!       Use the subroutine LSP_FOCWWIL (from Section 4) consistently to
C     ! derive cloud radiative properties and precipitation amount,
C     ! taking into account that cloud does not freeze as soon as it is
C     ! cooled below the freezing point of bulk water.  The release of
C     ! latent heat of fusion (not a major term) is done differently in
C     ! order to allow energy conservation (UMDP 29).  This is the
C     ! reason for two layer cloud water contents being passed in and
C     ! then combined and differently split.
C
C     !   Dimensions:
      INTEGER NBANDS         ! Number of spectral bands in the shortwave
      PARAMETER (NBANDS=4)   ! This run uses the standard set of
C shortwave bands as described by Slingo (May 1989, J.Atmos.Sci., 46,
C 10, 1419-1427) or UMDP 23.
      INTEGER NGASES
C Number of absorbing gases treated in the shortwave
      PARAMETER (NGASES=3)    !  Standard set is water vapour, ozone
C                             !  and carbon dioxide.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER NLKUPS          ! DIMENSION OF LOOK-UP TABLES
      PARAMETER (NLKUPS=50)
      INTEGER!, INTENT(IN) ::
     &     L1,                       ! Number of points in input arrays
     &     NDO,                      ! Number of points to be treated
     &     NLIT,                     ! Number of them to be sunlit
     &     NLEVS,                    ! Number of levels
     &     NCLDS,                    ! Number of possibly cloudy levels
     &     NWET,                     ! Number of levels with water vapor
     &     NOZONE                    ! Number of levels with ozone
C     !  Physical inputs:
      REAL!, INTENT(IN) ::
     &     H2OIN(L1,NWET), CO2,          ! Mass mixing ratios of
     &     O3IN(L1,NOZONE),              !         absorbing gases
     &     PSTIN(L1),                    ! Surface pressure
     &     ABIN(NLEVS+1), BBIN(NLEVS+1), ! As and Bs at layer boundaries
     &     LCAIN(L1,1/(NCLDS+1)+NCLDS),  ! Layer cloud fractional cover
     &     LCW1IN(L1,1/(NCLDS+1)+NCLDS), ! Layer cloud frozen and liquid
     &     LCW2IN(L1,1/(NCLDS+1)+NCLDS), !               water contents
C     !   These are specific cloud water contents, mass per unit mass,
C     !               and, as explained above, only their sum is used.
     &     CCAIN(L1),                    ! Convective Cloud Amount
     &     CCWPIN(L1),                   !      and condensed water path
     &     SALIIN(L1),                   ! (True) Surface Albedo for
     &     SAOSIN(L1,2),                 !  land & ice, for & open sea
     &     COSZIN(L1),                   ! Mean (cos solar zenith angle)
C                                        !     while the point is sunlit
     &     LIT(L1),                      ! Fraction the point is sunlit
     &     TAC(L1,NLEVS),                ! Atmospheric temperatures
     &     AICE(L1),                     ! Sea-ice fraction
     &     SCS,                          ! Solar Constant Scaling factor
C     ! - inverse-square factor which multiplies the solar constant to
C     ! get the normal solar irradiance at this day's earth-sun distance
     &     LUT(NLKUPS,NTRANS,NGASES,2),
C     ! Look-up tables of transmissivities for each gas and of
C     ! differences of their successive elements.
     &     PTS                           ! Time interval at which the
C     ! increments to be returned are to be added in ("physics
C     ! timestep").  The time interval over which they are valid
C     ! ("shortwave timestep") is not used directly here, but as an
C     ! input to the astronomy code it affects COSZIN, LIT and LIST.
      INTEGER!, INTENT(IN) ::
     &     LIST(NLIT),               ! List of the NLIT sunlit points
     &     CCBIN(L1),                ! Convective cloud base & top,
     &     CCTIN(L1)                 ! layer boundaries counting up from
C                                    !                 the surface as 1
C     !  Control quantities:
      LOGICAL!, INTENT(IN) ::
     &     LAND(L1)                  ! Land/sea mask (.TRUE. for land)
     &     , OSON, CSOSON            ! Are OSDIA & CSOSDI wanted ?
     &     , NSS1ON, TDSSON          ! And are NSSSB1 and TDSS ?
     &     , CREFFO, LREFFO          ! And are CREFF and LREFF...
     &     , CVAMTO, LRAMTO          ! ... and CVAMT and LRAMT ?
     &     , CWPAJON                 ! Is CWP O/P wanted?
     &     , MICRO                   ! Is microphysics code activated?
     &     , LCA3ON                  !  And is LCA3L ?
C     ! Note that if LCLD3, LCA3L is needed to calculate TCASW & so
C     !  will be calculated whenever TCASWO or LCA3ON - so space must
C     !  then be available (via "implied diagnostics" in the std UM).
     &     , CSSSDO, CSSSUO          !       & are CSSSD & CSSSU,
     &     , LCASWO, CCASWO          !             LCASW & CCASW,
     &     , LCAARO, LCAAFO          !             LCAAR & LCAAF,
     &     , CCAARO, CCAAFO          !             CCAAR & CCAAF,
     &     , TCASWO                  !                 & TCASW ?
     &     , LCAARL(NCLDS),  LCAARB(NBANDS), LCAAFL(NCLDS)
     &     , LCAAFB(NBANDS), CCAARB(NBANDS), CCAAFB(NBANDS)
C     !  If L/C CAA R/F are wanted, for which (levels and) bands ?
     &    , LCLD3
     &    , L_CLOUD_WATER_PARTITION
C     !  And outputs:
      REAL!, INTENT(OUT) ::
     &     SWOUT(L1,NLEVS+2),        ! This is filled by SWMAST with the
C     !  normalized net downward shortwave flux at all layer boundaries.
C     !  SWRAD multiplies them by the normal incoming insolation to give
C     !  dimensioned fluxes (still not the actual fluxes as the cosz
C     !  term is not put in here) and differences them in the vertical
C     !  to give SW heating rates (except for the cosz) in each
C     !  atmospheric layer, leaving a surface net downward SW flux in
C     !  the first level for use in the surface scheme.  It also
C     !  modifies the latter so that it refers to land-and-ice only (the
C     !  surfaces dealt with in the atmospheric model), being the value
C     !  over that surface (except the cosz) times the fraction of the
C     !  grid-box covered by land or sea-ice.
C     !    The 'level' NLEV+2 holds NSSB1 without Zenith Angle
C     !  adjustment,for use in physics timesteps in RAD_CTL and CLD_CTL
     &     SWSEA(L1)                 ! The net downward SW flux over
C     !  open sea.  SWMAST returns this normalized and SWRAD converts
C     !  it into an actual flux with weighting by the open sea fraction
C     !  (so that it can be added to the corresponding land-and-ice
C     !  term to give the overall net downward SW flux.)
     &     , NTSWIN(L1)            !  Net SW absorption by the planet
     &     , OSDIA(L1)                ! Diagnosed actual and clear-sky
     &     , CSOSDI(L1)               !            outgoing solar at toa
     &     , CSSSD(L1)                ! Clear-sky total downward &
     &     , CSSSU(L1)                !   upward SW flux at the surface
     &     , LCASW(L1,NCLDS)          ! Layer/Convective Cloud Amount
     &     , CCASW(L1)                !    in SW (zero at night points)
     &     , LCAAR(L1,*)              ! Layer/Convective Cloud Amount *
     &     , LCAAF(L1,*)              !    Albedo to diRect and diFfuse
     &     , CCAAR(L1,*)              !    light (set to zero at night
     &     , CCAAF(L1,*)              !    points)
     &     , TCASW(L1)                !   Total cloud amount in SW
C     ! (i.e. fraction of the grid-box with cloud at some level)
     &     , NSSB1(L1)
C     !   Net downward SW flux at the surface in band 1
     &     , TDSS(L1)
C     !   Total downward SW flux at the surface (multiply-reflected
C     !   light being multiply counted).
     &     , TDSSB1(L1)
C     !   Total downward SW flux at surface in band 1
     &     , CREFF(L1)                ! Convective cloud rE * cld amount
     &     , LREFF(L1,NCLDS)          ! Layer cloud rE * cld amount
     &     , CVAMT(L1)                ! Convective cloud amount in SWRAD
     &     , LRAMT(L1,NCLDS)          ! Layer cloud amount in SWRAD
     &     , CWPAJ(L1,NCLDS)          ! Lyr cld CWP for 3-cld scheme
     &     , LCA3L(L1,NCLDS)          ! Diagnostic of layer cloud amount
C     ! restricted to 3 layers, calculated at all points on SW timesteps
C*
C     !  Constants:
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
!
! Description:
!
!  Contains various parameters for the SW effective radius
!  parametrisation, defined for land and sea areas.
!
!  NTOT is the number concentration (m-3) of activated CCN;
!  KPARAM is the ratio of volume mean radius to effective radius;
!  DCONRE is the effective radius (m) for deep convective clouds.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!
      REAL NTOT_LAND,
     &     NTOT_SEA,
     &     KPARAM_LAND,
     &     KPARAM_SEA,
     &     DCONRE_LAND,
     &     DCONRE_SEA

      PARAMETER ( NTOT_LAND = 6.0E08,
     &            NTOT_SEA = 1.5E08,
     &            KPARAM_LAND = 0.67,
     &            KPARAM_SEA = 0.80,
     &            DCONRE_LAND = 9.5E-06,
     &            DCONRE_SEA = 13.5E-06 )
      REAL CPBYG                        ! Helps convert fluxes to
      PARAMETER ( CPBYG = CP / G )      !                 heating rates
      REAL SC                           !  Solar constant
      PARAMETER ( SC = 1365. )
      REAL REICE,                       !  rE for ice cloud
     &     RELIQ,                       !  rE for liquid cloud
     &     DRE                          !  RELIQ-REICE
      PARAMETER ( REICE = 30.E-6, RELIQ = 7.E-6, DRE = RELIQ - REICE )
      REAL COSMIN                       ! Minimum value for COSZ, to
      PARAMETER ( COSMIN = 1.E-4 )      !     avoid underflow in SWCLOP
C     !  Local variables:
      REAL NSI,                         ! Normal Solar Irradiance
     &     TEMPOR,                      ! Temporary store
     &     DACON1, DBCON1,              ! Conversion factors for turning
C     ! fluxes into increments - difference of As and Bs across the
C     ! current layer, times CPBYG and divided by the timestep.
     &     DACON2, DBCON2               ! Conversion factors for turning
C     ! mixing ratio into pathlength - difference of As and Bs across
C     ! the current layer, divided by g.
      REAL DCONRE,   ! Cloud droplet rE for deep convective clouds.
     &     SCONRE,   !   "      "    "   " shallow   "         "  .
     &     NTOT,     ! Total CCN concentration (/m**3).
     &     KPARAM,   ! k parameter (=rV/rE).
     &     PCCTOP,   ! Convective cloud top pressure.
     &     PCCBOT,   !      "       "   base    "   .
     &     LCMMR,    ! Layer cloud mass mixing ratio (kg/kg).
     &     LWC,      ! Cloud liquid water content (kg/m**3).
     &     RHOAIR,   ! Local density of air (kg/m**3).
     &     DELTAZ,   ! Thickness of convective cloud (m).
     &     PRESS1,   ! Pressure at bottom...
     &     PRESS2,   !        ...and top of layer boundaries.
     &     TAU,      ! Area-averaged optical depth.
     &     L1AJ      ! Cloud amount dummy-variable.
C
C*L
CL    !  Dynamically allocated workspace:
C     !  3*NDO+ NLIT*(3*NCLDS+NWET+NOZONE+4*NBANDS+8) +2*(NLEVS+1)
      REAL H2OGI(NLIT,NWET),             ! Gathered and inverted inputs:
     &     O3GI(NLIT,NOZONE),            ! just as the corresponding
     &     PSTGI(NLIT),                  ! ...IN arrays, except that the
     &     ABGI(NLEVS+1), BBGI(NLEVS+1), ! two LCW arrays are combined,
     &     LCAGI(NLIT,NCLDS),            ! since the ice/liquid split is
     &     LCWPGI(NLIT,NCLDS),           ! done differently for
     &     CCAGI(NLIT),                  ! radiation and precipitation
     &     CCWPGI(NLIT),                 ! than for latent heat release,
     &     COSZGI(NLIT),                 ! and also converted from cloud
     &                                   ! water content to path.
     &     SAGI(NLIT,NBANDS,2),          ! Gathered surface albedos for
     &     SAOSGI(NLIT,NBANDS,2)         ! each band, for the whole
C     ! grid-box and open sea only (for SWMAST to calculate SWSEA with)
C
      INTEGER CCBGI(NLIT),               ! Convective cloud base & top,
     &     CCTGI(NLIT)                   ! layers counting down from the
C                                        !                top layer as 1
     &     , INDEX(NDO)
C     !  Index for maximum(input)/only(used) cloud cover for a "type"
C     !  (This, and MAXCLD below, are dimensioned NDO rather than NLIT
C     !          because full field size is used if LCA3L is wanted.)
      REAL CRE(NLIT),                    ! Equivalent radii calculated
     &     LRE(NLIT,NCLDS),              ! as functions of temperature.
     &     LAYERE(NLIT,NCLDS),           ! Liquid-only rE
     &     CWPAJGI(NLIT,NCLDS),          ! CWP gathered & inverted
     &     MAXCLD(NDO),                  !  Maximum cloud cover & total
     &     TOTCWC(NLIT),                 !    water content for a "type"
     &     IITOA(NDO)                    ! Incoming Insolation at the
C                                        !         Top Of the Atmosphere
C*
      INTEGER LEVEL, J,                  ! Loopers over level and point
     &     BAND,                         !                    and band.
     &     OFFSET,                       ! Index for diagnostics SWRAD
C     ! returns (potentially) compressed, allowing just the bands or
C     ! level-and-band combinations wanted to be allocated and set.
     &     DIRDIF,                       !    and direct/diffuse albedos
     &     TYPE,                         !       & cloud "type" (H/M/L)
     &     RANGE(3,2),                   ! The range of level numbers
C     !  (counting down from the highest potentially cloudy level) for
C     !  the 3 cloud "types" - i.e. the RANGE(n,1)th to RANGE(n,2)th
C     !  potentially cloudy levels are assigned to the nth cloud type.
C     !  The values are set by comparing model eta values with BOUNDS.
     &     FSTLEV,                       ! The equivalent of RANGE for
     &     LSTLEV,                       !  a particular cloud type, but
C                                        !  counting up from the surface
     &     NCLEAR,                       ! NLEVS-NCLDS
     &     NNIGHT,                       ! NDO-NLIT
     &     NLP1B2                        ! (NLEVS+1)/2
      REAL BOUNDS(2),                    ! Eta values that define where
C     ! cloud changes from "high" to "medium", & from "medium" to "low"
     &     ETA,                          ! Eta at the layer boundary
C     !                                  !    currently being checked
     &     ETALST                        !       & the previous one
     &     , FOCWWIL
!       Local value of Fraction Of Cloud Water Which Is Liquid
     &     , TFOC
!       and the cloud temperature used to calculate it.
      LOGICAL SET                        ! Has RANGE been set yet ?
      DATA BOUNDS / .37, .79 /
      DATA SET / .FALSE. /
      SAVE RANGE, SET                    ! SET must be specified too as
C     !   FORTRAN requires a variable initialized by a DATA statement to
C     !   have the SAVE attribute only if its value has not changed.
      IF (MICRO) THEN

C   Zero effective radius arrays if diagnostics requested:
        IF (CREFFO) THEN
          DO II=1, NDO
            CREFF(II) = 0.0
          END DO
        END IF
        IF (LREFFO) THEN
          DO JJ=1, NCLDS
            DO II=1, NDO
              LREFF(II,JJ) = 0.0
            END DO
          END DO
        END IF
C   Zero Cloud-Amount-In-SWRAD arrays if diagnostics requested:
        IF (CVAMTO) THEN
          DO II=1, NDO
            CVAMT(II) = 0.0
          END DO
        END IF
        IF (LRAMTO) THEN
          DO JJ=1, NCLDS
            DO II=1, NDO
              LRAMT(II,JJ) = 0.0
            END DO
          END DO
        END IF
C   Zero Layer-Cloud-CWP-In-SWRAD arrays if diagnostics requested:
        IF (CWPAJON) THEN
          DO JJ=1, NCLDS
            DO II=1, NDO
              CWPAJ(II,JJ)=0.0
            END DO
          END DO
        END IF

      END IF

CL
CL    !  Section 1 - invert and gather input data for SWMAST
CL       ~~~~~~~~~
CL    !  As & Bs of course only need inverting:
Cfpp$ NoConcur L
      DO 11 LEVEL=1, NLEVS+1
       ABGI(LEVEL) = ABIN(NLEVS+2-LEVEL)
       BBGI(LEVEL) = BBIN(NLEVS+2-LEVEL)
   11 CONTINUE
      NCLEAR = NLEVS - NCLDS
C
CL    ! &, if LCLD3 is on, the first time into the routine, find where
CL    ! cloud type boundaries will lie in terms of the numbering of this
CL    !  run's eta levels:
C
      IF ( LCLD3 .AND. .NOT. SET ) THEN
        RANGE(1,1) = 1
        LEVEL = NCLEAR + 1
        DO J=1, 2
  101     ETA = BBGI(LEVEL) + ABGI(LEVEL) / PREF
          IF ( ETA .LT. BOUNDS(J) ) THEN
             LEVEL  = LEVEL + 1
             ETALST = ETA
C            ! This assumes the vertical resolution is not too crude in
C            !    the troposphere - but it would have to be rather worse
C            !    even than the old 11-layer Cyber climate model.
             GO TO 101
           ELSE
C            ! This has found the first layer boundary below BOUNDS -
C            !   is this or the previous one closer ?
             IF ( BOUNDS(J)-ETALST .LT. ETA-BOUNDS(J) ) LEVEL = LEVEL-1
             RANGE(J+1,1) = LEVEL - NCLEAR
             RANGE(J,2)   = RANGE(J+1,1) - 1
          ENDIF
        ENDDO
        RANGE(3,2) = NCLDS
        SET = .TRUE.
      ENDIF
C
C
CL    !  while single-level or no-level data would just need gathering
C     !  - except that convective cloud rE must be calculated from the
C     !  temperature of the highest layer the cloud extends into, and
C     !  convective cloud base and top must be altered to count from the
C     !  top down and to refer to layer centres rather than layer
C     !  boundaries, and constrained to have a valid value (where CCA=0,
!        P27 does not set CCB or CCT.)  MAXCLD is used as temporary
!        storage for the gathered temperature input to ROCWWIP (also
!        used later by the microphsyics option), and CRE for the output.
      DO J=1, NLIT
       PSTGI(J) = PSTIN(LIST(J))
       CCAGI(J) = CCAIN(LIST(J))
       CCWPGI(J)= CCWPIN(LIST(J))
C      !  Conversion of CCWP here omitted for the time being.
       COSZGI(J)= COSZIN(LIST(J))
       IF ( COSZGI(J) .LT. COSMIN )  COSZGI(J) = COSMIN
       CCTGI(J) = NLEVS+2 - CCTIN(LIST(J))
       IF (  CCTGI(J) .GT. NLEVS  .OR.  CCTGI(J) .LE. NCLEAR  )
     &     CCTGI(J) = NCLEAR + 1
       CCBGI(J) = NLEVS+1 - CCBIN(LIST(J))
       IF (  CCBGI(J) .GT. NLEVS  .OR.  CCBGI(J) .LE. NCLEAR  )
     &     CCBGI(J) = NLEVS
C      ! CCTGI (where it was defined) was indexed similarly to TAC, but
C      !  we would have to subtract 1 to get the temperature at the
C      !  layer centre BELOW the layer boundary indicated by CCT.  To
C      !  be sure we do not access outside the valid range, we must
C      !  actually use CCTGI, which makes it a little less clear.
       MAXCLD(J) = TAC(LIST(J),NLEVS+1-CCTGI(J))
      END DO
      CALL LSP_FOCWWIL (MAXCLD, NLIT, CRE)
      DO J=1, NLIT
       TFOC = MAXCLD(J)
       FOCWWIL = CRE(J)
      IF (MICRO) THEN

       IF (LAND(LIST(J))) THEN
         DCONRE = DCONRE_LAND       ! Continental clouds.
         KPARAM = KPARAM_LAND
         NTOT = NTOT_LAND
       ELSE
         DCONRE = DCONRE_SEA        ! Maritime clouds.
         KPARAM = KPARAM_SEA
         NTOT = NTOT_SEA
       END IF
       IF (CCAGI(J).LE.0.0) THEN
         CRE(J)=0.0                    ! Set rE to zero for no cloud.
       ELSE
         PCCTOP=ABIN(CCTIN(LIST(J)))+BBIN(CCTIN(LIST(J)))*PSTGI(J)
         PCCBOT=ABIN(CCBIN(LIST(J)))+BBIN(CCBIN(LIST(J)))*PSTGI(J)
         DELTAZ=(R*TFOC/G)*ALOG(PCCBOT/PCCTOP)
         IF (DELTAZ .LT. 500.0) THEN             ! Shallow convection.
           LWC=(CCWPGI(J)/DELTAZ)
           SCONRE=(3.0*LWC/(4.0*PI*RHO_WATER*KPARAM*NTOT))**(1.0/3.0)
           CRE(J)=REICE+(SCONRE-REICE)*FOCWWIL
C         Set safe rE limits (for SWCLOP):
           IF (CRE(J).LT.0.35E-06) CRE(J)=0.35E-06
           IF (CRE(J).GT.37.0E-06) CRE(J)=37.0E-06
         ELSE
           CRE(J)=REICE+(DCONRE-REICE)*FOCWWIL   ! Deep convection.
         END IF
       END IF
       IF (CREFFO) CREFF(LIST(J))=CRE(J) * CCAGI(J) * 1000000.0
       IF (CVAMTO) CVAMT(LIST(J))=CCAGI(J) * 1000000.0

      ELSE

       CRE(J) = REICE + DRE * FOCWWIL

      END IF

      ENDDO
C
CL    !  Water is gathered and inverted at NWET levels:
      DO 14 LEVEL=1, NWET
Cfpp$  Select(CONCUR)
       DO 14 J=1, NLIT
        H2OGI(J,LEVEL) = H2OIN(LIST(J),NWET+1-LEVEL)
   14 CONTINUE
C
CL    !  and ozone at NOZONE...
      DO 15 LEVEL=1, NOZONE
Cfpp$  Select(CONCUR)
       DO 15 J=1, NLIT
        O3GI(J,LEVEL) = O3IN(LIST(J),NOZONE+1-LEVEL)
   15 CONTINUE
C
CL    !  Layer cloud data are gathered and inverted at NCLDS levels.
C     !  rE is calculated as for convective cloud,
C     !  and also QL & QF are added together.
      DO 16 LEVEL=1, NCLDS
       DACON2 = ( ABIN(NCLDS+1-LEVEL) - ABIN(NCLDS+2-LEVEL) ) / G
       DBCON2 = ( BBIN(NCLDS+1-LEVEL) - BBIN(NCLDS+2-LEVEL) ) / G
Cfpp$  Select(CONCUR)
       DO J=1, NLIT
        LCAGI(J,LEVEL) = LCAIN(LIST(J),NCLDS+1-LEVEL)
        MAXCLD(J) = TAC(LIST(J),NCLDS+1-LEVEL)
       END DO
       IF (L_CLOUD_WATER_PARTITION)  THEN
!   calculate proportion of liquid water focwwil as ratio qcl/(qcl+qcf)
         DO J=1, NLIT
           IF (LCAGI(J,LEVEL) .GT. 0.) THEN
             LRE(J,LEVEL) = LCW1IN(LIST(J),NCLDS+1-LEVEL) /
     &     (LCW1IN(LIST(J),NCLDS+1-LEVEL)+LCW2IN(LIST(J),NCLDS+1-LEVEL))
           ELSE
!          Arbitrary number: makes it safe & vectorizable
             LRE(J,LEVEL) = 0.0
           ENDIF
         END DO
       ELSE
!   set proportion of liquid water focwwil from parametrized function
         CALL LSP_FOCWWIL (MAXCLD, NLIT, LRE(1,LEVEL))
       ENDIF
!
       DO J=1, NLIT
        TFOC = MAXCLD(J)
        FOCWWIL = LRE(J,LEVEL)
      IF (MICRO) THEN

        IF (LAND(LIST(J))) THEN
          KPARAM = KPARAM_LAND
          NTOT = NTOT_LAND
        ELSE
          KPARAM = KPARAM_SEA
          NTOT = NTOT_SEA
        END IF
        LCMMR = ( LCW1IN(LIST(J), NCLDS+1-LEVEL)
     &          + LCW2IN(LIST(J), NCLDS+1-LEVEL) )
        IF (LCAGI(J,LEVEL) .GT. 0.0) THEN
          LCMMR = LCMMR / LCAGI(J,LEVEL)
          PRESS1=ABIN(NCLDS+1-LEVEL)+BBIN(NCLDS+1-LEVEL)*PSTGI(J)
          PRESS2=ABIN(NCLDS+2-LEVEL)+BBIN(NCLDS+2-LEVEL)*PSTGI(J)
          RHOAIR=(EXP((ALOG(PRESS1)+ALOG(PRESS2))/2.0)) / (R*TFOC)
          LWC=LCMMR * RHOAIR
          IF (LEVEL .GE. RANGE(3,1)) THEN                   ! Low cloud
            LAYERE(J,LEVEL)=(6.0*LWC/(4.0*PI*RHO_WATER*KPARAM*NTOT))
     &                                                      **(1.0/3.0)
          ELSE
            LAYERE(J,LEVEL)=(3.0*LWC/(4.0*PI*RHO_WATER*KPARAM*NTOT))
     &                                                      **(1.0/3.0)
          END IF
          LRE(J,LEVEL)=REICE+(LAYERE(J,LEVEL)-REICE)*FOCWWIL
C               Set safe rE limits (for SWCLOP):
          IF (LRE(J,LEVEL).LT.0.35E-06) LRE(J,LEVEL)=0.35E-06
          IF (LRE(J,LEVEL).GT.37.0E-06) LRE(J,LEVEL)=37.0E-06
        ELSE
          LRE(J,LEVEL)=0.0
          LAYERE(J,LEVEL)=0.0
        END IF

      ELSE

        LRE(J,LEVEL) = REICE + DRE * FOCWWIL

      END IF

        LCWPGI(J,LEVEL) = ( DACON2 + DBCON2 * PSTGI(J) ) *
     & ( LCW1IN(LIST(J),NCLDS+1-LEVEL) + LCW2IN(LIST(J),NCLDS+1-LEVEL) )
        IF ( ( .NOT. LCLD3 ) .AND. LCAGI(J,LEVEL) .GT. 0. )
     &         LCWPGI(J,LEVEL)= LCWPGI(J,LEVEL) / LCAGI(J,LEVEL)
       END DO
   16 CONTINUE
CL    ! If the option to combine layer clouds into 3 layers is on, do so
      IF ( LCLD3 ) THEN
C
CL    ! Now, find which layer holds most cloud of each "type":
C     !  (The loops over TYPE, and over LEVEL inside them, are from the
C     !  top down, as usual for loops involving TYPE or ..GI arrays.)
C
      DO TYPE=1, 3
Cfpp$   Select(CONCUR)
        DO J=1, NLIT
          TOTCWC(J) = LCWPGI(J,RANGE(TYPE,1))
          MAXCLD(J) = LCAGI(J,RANGE(TYPE,1))
          INDEX(J)  = RANGE(TYPE,1)
        ENDDO
        DO LEVEL=RANGE(TYPE,1)+1, RANGE(TYPE,2)
Cfpp$     Select(CONCUR)
          DO 161 J=1, NLIT
            TOTCWC(J) = TOTCWC(J) + LCWPGI(J,LEVEL)
            IF ( MAXCLD(J) .LT. LCAGI(J,LEVEL) ) THEN
              MAXCLD(J) = LCAGI(J,LEVEL)
              INDEX(J)  = LEVEL
            ENDIF
  161     CONTINUE ! Next J
        ENDDO                                  ! Next LEVEL
C
CL      !  and use it to set the values in the array passed to SWMAST:
C
C       !  We have the level of maximum cover for each type in the input
C       !   data, which will be the only one left non-zero.  Its CWC is
C       !   set to the sum of the CWC in all the levels of that "type"
C       !   (this sum being done on the grid-box means, which will then
C       !   be converted to an in-cloud value using the selected
C       !   (maximum) cloud amount).  The other levels' CWC and the rE
C       !   are not altered.
        DO LEVEL=RANGE(TYPE,1), RANGE(TYPE,2)
Cfpp$     Select(CONCUR)
          DO 162 J=1, NLIT
            IF ( LEVEL .EQ. INDEX(J) ) THEN
                IF ( LCAGI(J,LEVEL) .GT. 0. )
     &               TOTCWC(J) = TOTCWC(J) / LCAGI(J,LEVEL)
               LCWPGI(J,LEVEL) = TOTCWC(J)
               IF (MICRO) CWPAJGI(J,LEVEL) = LCWPGI(J,LEVEL)
             ELSE
               LCAGI(J,LEVEL)  = 0.
               IF (MICRO) CWPAJGI(J,LEVEL) = 0.0
            ENDIF
  162     CONTINUE                            ! Next J
        ENDDO                                 ! Next LEVEL
      ENDDO                                   ! Next TYPE
C
C     !  If wanted repeat the reduction-to-three-cloud-layers, but now
C     !  for all points & the other way up, for the diagnostic LCA3L.
C     ! This must be done if this diagnostic is wanted in its own right
C     !   or if TCASW is, as the latter is calculated from it.
C     !  (The loop over TYPE is still from the top down, but the loops
C     !  over LEVEL are now from the bottom up, to match how the clouds
C     !  are input and the output has to be output.)
C
      IF ( LCA3ON .OR. TCASWO ) THEN
        DO TYPE=1, 3
          FSTLEV = NCLDS + 1 - RANGE(TYPE,2)
          LSTLEV = NCLDS + 1 - RANGE(TYPE,1)
Cfpp$     Select(CONCUR)
          DO J=1, NDO
            MAXCLD(J) = LCAIN(J,FSTLEV)
            INDEX(J)  = FSTLEV
          ENDDO
          DO LEVEL=FSTLEV+1, LSTLEV
Cfpp$       Select(CONCUR)
            DO 163 J=1, NDO
              IF ( MAXCLD(J) .LT. LCAIN(J,LEVEL) ) THEN
                MAXCLD(J) = LCAIN(J,LEVEL)
                INDEX(J) = LEVEL
              ENDIF
  163       CONTINUE                             ! Next J
          ENDDO                                  ! Next LEVEL
          DO LEVEL=FSTLEV, LSTLEV
Cfpp$       Select(CONCUR)
            DO 164 J=1, NDO
              IF ( LEVEL .EQ. INDEX(J) ) THEN
                 LCA3L(J,LEVEL) = MAXCLD(J)
               ELSE
                 LCA3L(J,LEVEL) = 0.
              ENDIF
  164       CONTINUE                            ! Next J
          ENDDO                                 ! Next LEVEL
        ENDDO                                   ! Next TYPE
      END IF
      ENDIF   !  LCLD3
C
      IF (MICRO) THEN

        DO II=1,NCLDS
          DO JJ=1,NLIT
          L1AJ=LCAGI(JJ,II)
          IF (L1AJ .GT. 0.0) THEN
            TAU=(1.5*CWPAJGI(JJ,II)/(1000.0*LRE(JJ,II)))*L1AJ
          ELSE
            TAU=0.0
          END IF
          IF (TAU .LT. 5.0) L1AJ = 0.0
            IF (LREFFO) THEN
              LREFF(LIST(JJ),NCLDS+1-II) = LAYERE(JJ,II)*L1AJ*1.0E06
            END IF
            IF (LRAMTO) THEN
              LRAMT(LIST(JJ),NCLDS+1-II) = L1AJ * 1.0E06
            END IF
            IF (CWPAJON) THEN
              CWPAJ(LIST(JJ),NCLDS+1-II) = CWPAJGI(JJ,II) * L1AJ
            END IF
          ENDDO
        ENDDO

      END IF
C
CL    !  Gathering the clear-sky surface albedos, multiple copies are
CL    !  needed as P234 code expects band-dependent ones, which P233
CL    !                                          does not yet produce.
      DO 171 DIRDIF=1, 2
Cfpp$  Select(CONCUR)
       DO 17 J=1, NLIT
        SAOSGI(J,1,DIRDIF) = SAOSIN(LIST(J),DIRDIF)
        SAGI(J,1,DIRDIF)   = SALIIN(LIST(J))
        IF ( .NOT. LAND(LIST(J)) )
     &   SAGI(J,1,DIRDIF) = SAGI(J,1,DIRDIF) * AICE(LIST(J)) +
     &                    SAOSGI(J,1,DIRDIF) * ( 1.-AICE(LIST(J)) )
   17  CONTINUE
       DO 171 BAND=2, NBANDS
Cfpp$   Select(CONCUR)
        DO 171 J=1, NLIT
         SAGI(J,BAND,DIRDIF)   = SAGI(J,1,DIRDIF)
         SAOSGI(J,BAND,DIRDIF) = SAOSGI(J,1,DIRDIF)
  171 CONTINUE
C
C     !  Diagnose cloud-if-sunlit if wanted:
C
      IF ( CCASWO ) THEN
        DO J=1, NDO
          CCASW(J) = 0.0
        END DO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          CCASW(LIST(J)) = CCAGI(J)
        END DO
      END IF
      IF ( LCASWO ) THEN
        DO LEVEL=1, NCLDS
Cfpp$     Select(Concur)
          DO J=1, NDO
            LCASW(J,LEVEL) = 0.0
          END DO
CDir$     IVDep
Cfpp$     NoConcur L
          DO J=1, NLIT
            LCASW(LIST(J),LEVEL) = LCAGI(J,NCLDS+1-LEVEL)
          END DO
        END DO
      END IF
C
CL    !  Set NNIGHT, the number of night points to be treated by this
CL    !                                                  CALL to SWRAD
      NNIGHT=NDO-NLIT
C
CL
CL    ! Section 2 - CALL SWMAST
CL      ~~~~~~~~~
      CALL SWMAST (H2OGI, CO2, O3GI, PSTGI, ABGI, BBGI, LCAGI, LCWPGI,
     &     LRE, CCAGI, CCWPGI, CRE, CCBGI, CCTGI, COSZGI,
     &     SAGI, SAOSGI, LUT,
     &     CSOSDI(1+NNIGHT), CSOSON, NSSB1(1+NNIGHT), NSS1ON,
     &     TDSS(1+NNIGHT), TDSSON,
     &     CSSSD(1+NNIGHT),   CSSSDO, CSSSU(1+NNIGHT), CSSSUO,
     &     LCAAR(1+NNIGHT,1), LCAARO, LCAARL, LCAARB,
     &     LCAAF(1+NNIGHT,1), LCAAFO, LCAAFL, LCAAFB,
     &     CCAAR(1+NNIGHT,1), CCAARO, CCAARB,
     &     CCAAF(1+NNIGHT,1), CCAAFO, CCAAFB,
     &     NLIT, NLEVS, NCLDS,                                          
     &     NWET, NOZONE, NLIT, L1, SWSEA(1+NNIGHT), SWOUT(1+NNIGHT,1) )
C
C
CL    ! Also, zero areas of SWOUT & SWSEA that will not be set by SWMAST
C
C     !        (They are multiplied, here or in the control routines,
C     !   by the mean cosz for each physics timestep, i.e. zero at night
C     !   points, but this would fail if a word were not a valid real.)
C     !
      IF ( NDO.GT.NLIT ) THEN
        DO 20 LEVEL=1, NLEVS+2
Cfpp$    Select(CONCUR)
         DO 20 J=1, NNIGHT
          SWOUT(J,LEVEL) = 0.
   20   CONTINUE
        DO J=1, NNIGHT
          SWSEA(J) = 0.
        ENDDO
      ENDIF
C
C
CL    !  Section 3 - convert normalized net downward flux to atmospheric
CL    !  ~~~~~~~~~   heating rates and surface actual net downward flux
C
CL    !  Set up normalized-to-actual flux conversion factors:
CL    !  the incoming insolation at the top of the atmosphere
C
      NSI = SC * SCS
      DO 31 J=1, NDO
       IITOA(J) = NSI * COSZIN(J) * LIT(J)
   31 CONTINUE
C
CL    ! and set COSZGI to the same for daylit points
C
      DO 32 J=1, NLIT
        COSZGI(J) = IITOA(LIST(J))
   32 CONTINUE
C
CL    !  Fill NTSWIN:
C
      DO J=1, NDO
        NTSWIN(J) = 0.
      ENDDO
C
CDir$   IVDep
Cfpp$   NoConcur L
        DO 323 J=1, NLIT
        NTSWIN(LIST(J)) = COSZGI(J) * SWOUT(J+NNIGHT,1)
  323   CONTINUE
C
C     !  Before flux-differencing, diagnose outgoing solar if wanted :
C
      IF ( OSON ) THEN
        DO J=1, NDO
          OSDIA(J) = 0.
        ENDDO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          OSDIA(LIST(J)) = COSZGI(J) * ( 1. - SWOUT(J+NNIGHT,1) )
        ENDDO
      ENDIF
CL
CL    !  and if CSOSDI is wanted, scatter it back and convert it from
CL    !  normalized to actual flux:
CL
      IF ( CSOSON ) THEN
        DO J=1, NNIGHT
          CSOSDI(J) = 0.
        ENDDO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          CSOSDI(LIST(J)) = CSOSDI(J+NNIGHT)
        ENDDO
        DO J=1, NDO
          CSOSDI(J) = IITOA(J) * CSOSDI(J)
        ENDDO
      ENDIF
C
CL    !  Scatter NSSB1 back and convert from normalized to actual flux
C     !     (including multiplication by open-sea fraction), and set to
C     !     zero over land:
C
      IF( NSS1ON) THEN
        DO J=1, NNIGHT
          NSSB1(J) = 0.
        ENDDO
CDir$ IVDep
Cfpp$ NoConcur L
        DO J=1, NLIT
          NSSB1(LIST(J)) = NSSB1(J+NNIGHT)
        ENDDO
C Set NSSB1 over both land and sea surface
        DO J=1, NDO
          IF ( LAND(J) ) THEN
            NSSB1(J) = IITOA(J) * NSSB1(J)
          ELSE
            NSSB1(J) = IITOA(J) * ( 1. - AICE(J) ) * NSSB1(J)
          ENDIF

C Find total downward SW flux in band 1
          TDSSB1(J) = NSSB1(J) / (1.0 - SALIIN(J))
C (Albedo should never equal 1.0)
C Store TDSSB1 without zenith angle adjustment in SWOUT
          IF(IITOA(J).NE.0.0) THEN
            SWOUT(J,NLEVS+2) = TDSSB1(J) / (COSZIN(J) * LIT(J))
          ENDIF
        ENDDO                   ! NDO

      ELSE                      ! NSS1ON is false
C Photosynthetically active radiation not required, but initialise to
C  zero to avoid possible problems accessing uninitialised data later.
        DO J=1,NDO
           SWOUT(J,NLEVS+2) = 0.0
        ENDDO                   ! NDO 

      ENDIF                     ! NSS1ON
C
CL    !  Scatter TDSS back and convert from normalized to actual flux:
C
      IF ( TDSSON ) THEN
        DO J=1, NNIGHT
          TDSS(J) = 0.
        ENDDO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          TDSS(LIST(J)) = TDSS(J+NNIGHT)
        ENDDO
        DO J=1, NDO
          TDSS(J) = IITOA(J) * TDSS(J)
        ENDDO
      ENDIF
C
CL    !  And the same for CSSSD and CSSSU:
C
      IF ( CSSSDO ) THEN
        DO J=1, NNIGHT
          CSSSD(J) = 0.
        ENDDO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          CSSSD(LIST(J)) = CSSSD(J+NNIGHT)
        ENDDO
        DO J=1, NDO
          CSSSD(J) = IITOA(J) * CSSSD(J)
        ENDDO
      ENDIF
      IF ( CSSSUO ) THEN
        DO J=1, NNIGHT
          CSSSU(J) = 0.
        ENDDO
CDir$   IVDep
Cfpp$   NoConcur L
        DO J=1, NLIT
          CSSSU(LIST(J)) = CSSSU(J+NNIGHT)
        ENDDO
        DO J=1, NDO
          CSSSU(J) = IITOA(J) * CSSSU(J)
        ENDDO
      ENDIF
C
CL    !  and cloud albedo diagnostics:
C
      IF ( LCAARO ) THEN
        OFFSET = 1
        DO 338 BAND=1, NBANDS
          DO 338 LEVEL=1, NCLDS
            IF ( LCAARL(LEVEL) .AND. LCAARB(BAND) ) THEN
CDir$         IVDep
Cfpp$         NoConcur L
              DO J=1, NLIT
                LCAAR(LIST(J),OFFSET) = LCAAR(J+NNIGHT,OFFSET)
              ENDDO
CDir$         IVDep
              DO J=1, NDO
                IF ( LIT(J) .EQ. 0. ) LCAAR(J,OFFSET) = 0.
              ENDDO
              OFFSET = OFFSET + 1
            ENDIF
  338   CONTINUE
      ENDIF
      IF ( LCAAFO ) THEN
        OFFSET = 1
        DO 337 BAND=1, NBANDS
          DO 337 LEVEL=1, NCLDS
            IF ( LCAAFL(LEVEL) .AND. LCAAFB(BAND) ) THEN
CDir$         IVDep
Cfpp$         NoConcur L
              DO J=1, NLIT
                LCAAF(LIST(J),OFFSET) = LCAAF(J+NNIGHT,OFFSET)
              ENDDO
CDir$         IVDep
              DO J=1, NDO
                IF ( LIT(J) .EQ. 0. ) LCAAF(J,OFFSET) = 0.
              ENDDO
              OFFSET = OFFSET + 1
            ENDIF
  337   CONTINUE
      ENDIF
      IF ( CCAARO ) THEN
        OFFSET = 1
        DO 336 BAND=1, NBANDS
          IF ( CCAARB(BAND) ) THEN
CDir$       IVDep
Cfpp$       NoConcur L
            DO J=1, NLIT
              CCAAR(LIST(J),OFFSET) = CCAAR(J+NNIGHT,OFFSET)
            ENDDO
CDir$       IVDep
            DO J=1, NDO
              IF ( LIT(J) .EQ. 0. ) CCAAR(J,OFFSET) = 0.
            ENDDO
            OFFSET = OFFSET + 1
          ENDIF
  336   CONTINUE
      ENDIF
      IF ( CCAAFO ) THEN
        OFFSET = 1
        DO 335 BAND=1, NBANDS
          IF ( CCAAFB(BAND) ) THEN
CDir$       IVDep
Cfpp$       NoConcur L
            DO J=1, NLIT
              CCAAF(LIST(J),OFFSET) = CCAAF(J+NNIGHT,OFFSET)
            ENDDO
CDir$       IVDep
            DO J=1, NDO
              IF ( LIT(J) .EQ. 0. ) CCAAF(J,OFFSET) = 0.
            ENDDO
            OFFSET = OFFSET + 1
          ENDIF
  335   CONTINUE
      ENDIF
C
CL    !  Invert SWOUT and scatter it and SWSEA back
C
CDir$ IVDep
Cfpp$ NoConcur L
      DO 33 J=1, NLIT
        SWSEA(LIST(J)) = SWSEA(J+NNIGHT)
   33 CONTINUE
      NLP1B2=(NLEVS+1)/2
CIf this were NLEVS/2+1, could omit special case (do (twice) as general)
      DO 34 LEVEL=1, NLP1B2
CDir$  IVDep
Cfpp$  NoConcur L
       DO 34 J=1, NLIT
        TEMPOR = SWOUT(J+NNIGHT,LEVEL)
        SWOUT(LIST(J),LEVEL) = SWOUT(J+NNIGHT,NLEVS+2-LEVEL)
        SWOUT(LIST(J),NLEVS+2-LEVEL) = TEMPOR
   34 CONTINUE
      IF ( NLEVS/2*2 .EQ. NLEVS ) THEN      ! Middle level: scatter only
CDir$   IVDep
Cfpp$   NoConcur L
        DO 35 J=1, NLIT
         SWOUT(LIST(J),LEVEL) = SWOUT(J+NNIGHT,LEVEL)
   35   CONTINUE
      ENDIF
C
CL    !  If wanted, diagnose total cloud amount as seen by the SW:
C
      IF ( TCASWO ) THEN
        IF ( LCLD3 ) THEN
           CALL SWDTCA (LCA3L, CCAIN, NCLDS, L1, NDO, TCASW)
         ELSE
           CALL SWDTCA (LCAIN, CCAIN, NCLDS, L1, NDO, TCASW)
        ENDIF
      ENDIF
C
CL    !  Convert fluxes to increments (Eq 1.1), and also put NSI in
C     !   - but omit cosz term (we could multiply by IITOA to get values
C     !   averaged over the whole SW timestep, but this is omitted so
C     !   that the control code can multiply by the correct mean cosz
C     !   for each physics timestep).  Also zero the heating rates for
C     !   night points in the later part of the scattered-back vector
C     !   - these should be multiplied by cosz=0 before being added in,
C     !   but there is the possibility of rounding-error-sized cosz
C     !   (from when the sun sets just as the timestep starts, or rises
C     !   just as it finishes) not being calculated consistently on some
C     !   machines, so it is safest to zero them in case, rather than
C     !   leave in the values for some day point which would then be
C     !   added in multiplied by a (very small) cosz to give (very
C     !   small) spurious and batching-dependent heating.
C
      DO 37 LEVEL=NLEVS, 1, -1
       DACON1 = ( ABIN(LEVEL) - ABIN(LEVEL+1) ) * CPBYG / ( PTS * NSI )
       DBCON1 = ( BBIN(LEVEL) - BBIN(LEVEL+1) ) * CPBYG / ( PTS * NSI )
       DO 38 J=1, NDO
        SWOUT(J,LEVEL+1) = ( SWOUT(J,LEVEL+1) - SWOUT(J,LEVEL) )
     &                                 / ( DACON1 + PSTIN(J) * DBCON1 )
   38  CONTINUE
       DO J=NNIGHT+1, NDO
        IF ( IITOA(J) .EQ. 0. ) SWOUT(J,LEVEL+1) = 0.
       ENDDO
   37 CONTINUE
C
CL    ! Finally, subtract the open-sea contribution from the total
CL    !  net downward surface flux to leave the land-and-sea-ice
CL    !  contribution, and convert both from normalized fluxes to
CL    !  dimensioned ones - they did not get multiplied by NSI as the
CL    !  atmospheric heating rates have just been.  The term to be used
CL    !  over land or sea-ice is not multiplied by the cos(solar zenith
CL    !  angle) term because this will be done for each physics
CL    !  timestep in the control routines (though again it is set to
CL    !  zero at night points), but SWSEA and NSSB1 are.
C
      DO 39 J=1, NDO
       IF ( LAND(J) ) THEN
          SWSEA(J) = 0.
        ELSE
          SWSEA(J)   = SWSEA(J) * ( 1.-AICE(J) )
          SWOUT(J,1) = SWOUT(J,1) - SWSEA(J)
          SWSEA(J)   = IITOA(J) * SWSEA(J)
       ENDIF
       SWOUT(J,1) = SWOUT(J,1) * NSI
   39 CONTINUE
      DO J=NNIGHT+1, NDO
        IF ( IITOA(J) .EQ. 0. ) SWOUT(J,1) = 0.
      ENDDO
C
      RETURN
      END
