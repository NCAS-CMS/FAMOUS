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
CLL Subroutine LWRAD    ----------------------------------------------
CLL
CLL Purpose :
CLL  It is the top-level routine in component P232. It CALLs LWMAST to
CLL  produce longwave fluxes (after setting convective cloud base and
CLL  top to safe values) and then differences these fluxes and returns
CLL  timestep-by timestep increments.  It will diagnose Outgoing
CLL  Longwave Radiation (OLR) if requested.
CLL  Before LWRAD is called, LWLKIN (in deck LWTRAN) must be CALLed to
CLL  initialize LUT.
CLL
CLL        Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.4    31/10/94  Stephanie Woodward
CLL Added 4 arguments (mmr's of minor gases) for compatibility with
CLL new LWRAD1C. They are declared but never used.
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions replaced by T3E lib functions.
CLL                       S.J.Swarbrick
!     4.4  12/03/97  Argument L_CLOUD_WATER_PARTITION passed into
!                    cloud calculation code. A Bushell
CLL
CLL Programming standard :
CLL  It technically conforms with standard A of UMDP 4 (version 2,
CLL  18/1/90)
CLL  If UPDATE *DEF CRAY is off, the code is standard FORTRAN 77
CLL  except for the addition of ! comments (and then sets the "vector
CLL  length" to 1) but if not it includes CRAY automatic arrays also.
CLL
CLL Logical components covered : P232, D23
CLL  It is the top-level routine in component P232.
CLL  P232 (longwave radiation), D23 (radiation diagnostics).
CLL
CLL Project task : P23
CLL
CLL External documentation:  (where appropriate) is in UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE LWRAD(H2O, CO2, O3, N2OMMR, CH4MMR, C11MMR, C12MMR,
     &      TAC, PEXNER, TSTAR, PSTAR, AB,
     &     BB, AC, BC, AICE, LCA, LCCWC1, LCCWC2, CCA, CCCWP, CCB, CCT,
     &     LAND, PTS, LUT,
     &     TCADIA, TCAON, CSOLRD, CSOLON, SFDN, SFDNON, CSSFDN, CSSDON,
     &     L_CLOUD_WATER_PARTITION,
     &     L2, NLEVS, NCLDS,                                            
     &     NWET, NOZONE, L1,                     OLR,  LWSEA,  LWOUT)
C
CL   External Routines called
      EXTERNAL LWMAST                ! Top-level of the LW physics
     &     , LWDCSF                  ! Diagnoses clear-sky fraction
C     !  Dimensions:
      INTEGER IT              ! Dimension of look-up tables
!  Dimension of Lookup Tables
      PARAMETER (IT=25)
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
! longwave bands as described by Slingo and Wilderspin
! (April 1986, Quart.J.R.Met.Soc., 112, 472, 371-386) or UMDP 23.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=4)    !  Standard set is water vapour line and
C continuum (which have to be treated separately because the pathlength
C scaling is different), ozone and carbon dioxide.
C     !  Array dimensions must be constants in FORTRAN:
      INTEGER!, INTENT(IN) ::
     &     L2,                       ! Number of points to be treated
     &     NLEVS,                    ! Number of levels
     &     NCLDS,                    ! Number of possibly cloudy levels
     &     L1,                       ! Full field dimension
     &     NWET,                     ! Number of levels with moisture
     &     NOZONE                    ! Number of levels with ozone
      REAL!, INTENT(IN) ::
     &     PSTAR(L1),                ! Surface pressure
     &     AB(NLEVS+1), BB(NLEVS+1), ! As and Bs at layer boundaries
     &     AC(NLEVS), BC(NLEVS),     ! As and Bs at layer centres
     &     H2O(L1,NWET), CO2,        ! Mixing ratios of the three
     &     O3(L1,NOZONE),            !               absorbing gases
     &     N2OMMR,                   ! mmrs for minor gases
     &     CH4MMR,                   ! not used in this version
     &     C11MMR,                   ! but included for compatibility
     &     C12MMR,                   ! with 1c
     &     TAC(L1,NLEVS),            ! Temperature at layer centres
     &     PEXNER(L1,NLEVS+1),       ! Exner function @ layer boundaries
     &     TSTAR(L1),                ! Surface temperature
     &     LUT(IT,NTRANS,NGASES,2),  ! Look-up tables for LWTRAN
     &     AICE(L1),                 ! Sea-ice fraction
     &     LCCWC1(L1,1/(NCLDS+1)+NCLDS), LCCWC2(L1,1/(NCLDS+1)+NCLDS),
C     ! Layer cloud condensed water contents (specific contents, mass
C     ! per unit mass).  Only the sum of these two fields is used.
     &     LCA(L1,1/(NCLDS+1)+NCLDS),! Layer cloud fractional cover
     &     CCCWP(L1),                ! Convective cloud fractional cover
     &     CCA(L1),                  !          and condensed water path
     &     PTS                       ! Time interval that increments to
C     ! be returned are to be added in at ("physics timestep").  The
C     ! interval over which they are used ("longwave timestep") has no
C     ! effect on the longwave code: it just sets how often it is CALLed
      INTEGER!, INTENT(IN) ::
     &     CCB(L1), CCT(L1)          ! Convective cloud base & top
      LOGICAL!, INTENT(IN) ::
     &     LAND(L1)                  ! Land/sea mask (.TRUE. for land)
     &     , CSOLON                  !  Is CSOLRD wanted ?
     &     , TCAON                   !                   & is TCADIA ?
     &     , SFDNON                  !                     & is SFDN ?
     &     , CSSDON                  !                    & is CSSFDN ?
     &     , L_CLOUD_WATER_PARTITION
      REAL!, INTENT(OUT) ::
     &     LWOUT(L1,NLEVS+1),        ! This is filled by LWMAST with the
C     !  net downward longwave flux at all layer boundaries.  LWRAD
C     !  converts these to atmospheric heating rates, leaving only the
C     !  surface fluxes in the first level.
     &     LWSEA(L1)                 ! Then it uses numbers LWPLAN has
C     !  put into LWSEA to separate out the total surface flux over the
C     !  grid-box into the open-ocean and solid-surface (sea-ice or
C     !  land) contributions and returns these in LWSEA and the first
C     !  level of LWOUT respectively.
     &     , OLR(L1)                 !  Outgoing Longwave Radiation
     &     , CSOLRD(L1)              ! and its clear-sky equivalent
     &     , TCADIA(L2)              ! Total Cloud Amount diagnostic
     &     , SFDN(L2)                ! Surface flux down diagnostic
     &     , CSSFDN(L1)              ! and its clear-sky equivalent
C*
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

      REAL CPBYG
      PARAMETER ( CPBYG = CP / G )
      REAL DACON, DBCON              ! Conversion factors for turning
C     ! fluxes into increments - difference of As and Bs across the
C     ! current layer, times CPBYG and divided by the timestep.
      INTEGER LEVEL, J               ! Loopers over level and point
      LOGICAL SFDNCA                 !  Is SFDN to be calculated by
C     !  LWMAST ? - set if either SFDNON or CSSDON is, because SFDN is
C     !  needed to find CSSSDN even if not wanted for its own sake.
C     !  Space for SFDN is assigned by the "implied diagnostics"
C     !  arrangements in that case, but SFDNON is only set if it is
C     !  wanted for its own sake.
C
CL    Section 1 - correct input data where necessary
CL    ---------
C
CL    ! Restrict convective cloud base and top to their physical range.
C
      DO 10 J=1, L2
       IF ( CCB(J) .GT. NCLDS .OR. CCB(J) .LT. 1 ) CCB(J) = 1
       IF ( CCT(J) .GT. (NCLDS+1) .OR. CCT(J) .LT. 2 ) CCT(J) = NCLDS+1
   10 CONTINUE
C
      SFDNCA = SFDNON .OR. CSSDON
C
CL    Section 2 - CALL LWMAST
CL    ---------
C
      CALL LWMAST (H2O, CO2, O3, TAC, PEXNER, TSTAR, PSTAR, AB, BB,
     &     AC, BC, AICE, LCA, LCCWC1, LCCWC2, CCA, CCCWP, CCB, CCT, LUT,
     &     CSOLRD, CSOLON, SFDN, SFDNCA, CSSFDN, CSSDON,
     &     L_CLOUD_WATER_PARTITION,
     &     L2, NLEVS, NCLDS,                                            
     &     NWET, NOZONE, L1,                           LWSEA,   LWOUT)
C
CL    Section 3 - convert fluxes to increments
CL    ---------
C
CL    !  but first copy the top layer into OLR:
C
      DO J=1, L2
        OLR(J) = - LWOUT(J,NLEVS+1)
      ENDDO
C
CL    !  Convert fluxes to increments (Eq 1.1) within the atmosphere,
CL    !  leaving the surface net downward flux at the beginning of LWOUT
C
      DO 30 LEVEL=NLEVS, 1, -1
       DACON = ( AB(LEVEL) - AB(LEVEL+1) ) * CPBYG / PTS
       DBCON = ( BB(LEVEL) - BB(LEVEL+1) ) * CPBYG / PTS
       DO 33 J=1, L2
        LWOUT(J,LEVEL+1) = ( LWOUT(J,LEVEL+1) - LWOUT(J,LEVEL) )
     &                             / ( DACON + PSTAR(J) * DBCON )
   33  CONTINUE
   30 CONTINUE
C
CL    ! Separate the contributions for a solid surface (land or sea-ice)
CL    ! to be added in by the model's surface scheme, from those for
CL    ! open sea, to be used in the ocean model.  Initially LWOUT
CL    ! has the actual box-mean flux and LWSEA has the difference
CL    ! between upward surface fluxes for open-sea and sea-ice.
CL    ! The values that will never be used (open-sea value at land
CL    ! points and solid-surface values at ice-free sea points) are
CL    ! zeroed so that the 2 fields will sum to the actual box-mean flux
      DO 35 J=1, L2
       IF (LAND(J)) THEN
          LWSEA(J) = 0.
        ELSE IF ( AICE(J) .EQ. 0. ) THEN
          LWSEA(J) = LWOUT(J,1)
          LWOUT(J,1) = 0.
        ELSE
C         ! Overall, LWOUT(,1) = AICE * ( LWOUT(,1) + (1.-AICE)*LWSEA )
          LWSEA(J) = (1.-AICE(J)) * ( LWOUT(J,1) - AICE(J)*LWSEA(J) )
          LWOUT(J,1) = LWOUT(J,1) - LWSEA(J)
       ENDIF
   35 CONTINUE
C
CL    Section 5 - CALL LWDCSF to calculate clear-sky fraction
CL    ---------             and put away total cloud amount if requested
C
      IF ( TCAON ) THEN
        IF ( NCLDS .GT. 0 ) THEN
           CALL LWDCSF (LCA, CCA, CCB, CCT, NCLDS, L1, L2, TCADIA)
           DO J=1, L2
             TCADIA(J) = 1. - TCADIA(J)
           ENDDO
         ELSE
           DO J=1, L2
             TCADIA(J) = 0.
           ENDDO
        ENDIF
      ENDIF
C
      RETURN
      END
