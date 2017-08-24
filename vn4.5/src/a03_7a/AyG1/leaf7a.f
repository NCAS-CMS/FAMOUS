C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
! Written by Peter Cox (February 1996)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
!**********************************************************************
      SUBROUTINE LEAF_C3 (LAND_FIELD,VEG_PTS,VEG_INDEX,FT
     &,                   DQ,APAR,TL,CA,OA,PSTAR,FSMC
     &,                   GL,AL,CI,RD)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points
!                                 !    on the land grid.
     &,FT                         ! IN Plant functional type.

      REAL
     & DQ(LAND_FIELD)             ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)           ! IN Absorbed PAR (W/m2)
     &,TL(LAND_FIELD)             ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)             ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)             ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_FIELD)          ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_FIELD)           ! IN Soil water factor.
     &,GL(LAND_FIELD)             ! OUT Leaf conductnace for H2O (m/s).
     &,AL(LAND_FIELD)             ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,RD(LAND_FIELD)             ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_FIELD)             ! OUT Internal CO2 pressure (Pa).
     &,ACR(LAND_FIELD)            ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,B1(LAND_FIELD)             !
     &,B2(LAND_FIELD)             !
     &,B3(LAND_FIELD)             ! WORK Coefficients of the quadratic.
     &,CCP(LAND_FIELD)            ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_FIELD)           ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_FIELD)          ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)          ! WORK Leaf conductnace for CO2 (m/s).
     &,KC(LAND_FIELD)             ! WORK Michaelis constant for CO2 (Pa)
     &,KO(LAND_FIELD)             ! WORK Michaelis constant for O2 (Pa).
     &,QTENF(LAND_FIELD)          ! WORK Q10 function.
     &,TAU(LAND_FIELD)            ! WORK CO2/O2 specificity ratio.
     &,TDEGC(LAND_FIELD)          ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)            ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)          ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).
     &,WL(LAND_FIELD)             ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)          ! WORK Carboxylation,
     &,WLITE(LAND_FIELD)          !      Light, and
     &,WEXPT(LAND_FIELD)          !      export limited gross
!                                 !      photosynthetic rates
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_FIELD)             ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).

      INTEGER
     & CLOS_INDEX(LAND_FIELD)     ! WORK Index of land points
!                                 !      with closed stomata.
     &,CLOS_PTS                   ! WORK Number of land points
!                                 !      with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)     ! WORK Index of land points
!                                 !      with open stomata.
     &,OPEN_PTS                   ! WORK Number of land points
!                                 !      with open stomata.

      INTEGER
     & J,L                        ! WORK Loop counters.

      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice
C-----------------------------------------------------------------------
C Functional Type dependent parameters
C-----------------------------------------------------------------------
      INTEGER
     + C3(NPFT)                   ! 1 for C3 Plants, 0 for C4 Plants.

      REAL
     + ALPHA(NPFT)                ! Quantum efficiency
C                                 ! (mol CO2/mol PAR photons).
     +,A_WL(NPFT)                 ! Allometric coefficient relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index (kg C/m2).
     +,A_WS(NPFT)                 ! Woody biomass as a multiple of
C                                 ! live stem biomass.
     +,B_WL(NPFT)                 ! Allometric exponent relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index.
     +,DGL_DM(NPFT)               ! Rate of change of leaf turnover
C                                 ! rate with moisture availability.
     +,DGL_DT(NPFT)               ! Rate of change of leaf turnover
C                                 ! rate with temperature (/K)
     +,DQCRIT(NPFT)               ! Critical humidity deficit
C                                 ! (kg H2O/kg air).
     +,ETA_SL(NPFT)               ! Live stemwood coefficient
C                                 ! (kg C/m/LAI).
     +,FSMC_OF(NPFT)              ! Moisture availability below
C                                 ! which leaves are dropped.
     +,F0(NPFT)                   ! CI/CA for DQ = 0.
     +,GLMIN(NPFT)                ! Minimum leaf conductance for H2O
     +,G_AREA(NPFT)               ! Disturbance rate (/360days).
     +,G_GROW(NPFT)               ! Rate of leaf growth (/360days).
     +,G_LEAF_0(NPFT)             ! Minimum turnover rate for leaves
!                                 ! (/360days).
     +,G_ROOT(NPFT)               ! Turnover rate for root biomass
!                                 ! (/360days).
     +,G_WOOD(NPFT)               ! Turnover rate for woody biomass
!                                 ! (/360days).
     +,KPAR(NPFT)                 ! PAR Extinction coefficient
C                                 ! (m2 leaf/m2 ground).
     +,LAI_MAX(NPFT)              ! Maximum projected LAI.
     +,LAI_MIN(NPFT)              ! Minimum projected LAI.
     +,NL0(NPFT)                  ! Top leaf nitrogen concentration
C                                 ! (kg N/kg C).
     +,NR_NL(NPFT)                ! Ratio of root nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,NS_NL(NPFT)                ! Ratio of stem nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,OMEGA(NPFT)                ! Leaf scattering coefficient
C                                 ! for PAR.
     +,R_GROW(NPFT)               ! Growth respiration fraction.
     +,SIGL(NPFT)                 ! Specific density of leaf carbon
C                                 ! (kg C/m2 leaf).
     +,TLEAF_OF(NPFT)             ! Temperature below which leaves are
C                                 ! dropped.
     +,TLOW(NPFT)                 ! Lower temperature for
C                                 ! photosynthesis (deg C)
     +,TUPP(NPFT)                 ! Upper temperature for
C                                 ! photosynthesis (deg C)

C----------------------------------------------------------------------
C                        BT     NT    C3G    C4G     S
C----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0,     1 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040,  0.08 /
      DATA A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10 /               
      DATA A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 /
      DATA B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667 /
      DATA DGL_DM  /   10.0,  10.0,   0.0,   0.0,  10.0 /               
      DATA DGL_DT  /    9.0,   0.0,   0.0,   0.0,   0.0 /
      DATA DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100 /
      DATA ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01 /
      DATA F0      /  0.875, 0.875, 0.900, 0.800, 0.900 /
      DATA FSMC_OF /   0.20,  0.20,  0.20,  0.20,  0.20 /               
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6 /
      DATA G_AREA  /  0.004, 0.004,  0.10,  0.10,  0.05 /
      DATA G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00 /
      DATA G_LEAF_0/   0.20,  0.20,  0.20,  0.20,  0.20 /               
      DATA G_ROOT  /   0.20,  0.20,  0.20,  0.20,  0.20 /
      DATA G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05 /               
      DATA KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50 /
      DATA LAI_MAX /  10.00, 10.00,  4.00,  4.00,  6.00 /
      DATA LAI_MIN /   4.00,  4.00,  1.00,  1.00,  1.00 /
      DATA NL0     /  0.050, 0.030, 0.060, 0.030, 0.030 /
      DATA NR_NL   /   1.00,  1.00,  1.00,  1.00,  1.00 /               
      DATA NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10 /
      DATA OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 /
      DATA R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA SIGL    / 0.0300,0.1000,0.0250,0.0500,0.0500 /
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,223.15 /               
      DATA TLOW    /    0.0,  -5.0,   0.0,  13.0,   0.0 /               
      DATA TUPP    /   36.0,  31.0,  36.0,  45.0,  36.0 /               

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL
     & BETA1, BETA2   ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4     ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4 ! Constant relating VCMAX and leaf N
!                     ! from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                     ! - assuming dry matter is 40% carbon by mass)
!                     ! and Jacobs 1994:
!                     ! C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                     ! (mol/m2/s)
     &,R              ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
!                     ! resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93
     &,          FDC3 = 0.015,  FDC4 = 0.025
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3
     &,          R = 8.3144  , RATIO = 1.6
     &,          ZERODEGC = 273.15)

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
        IF (FSMC(L).EQ.0.0 .OR. DQ(L).GE.DQCRIT(FT)
     &                     .OR. APAR(L).EQ.0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFFC3 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
        CCP(L) = 0.5 * OA(L) / TAU(L)


        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * QTENF(L)

      ENDDO

CDIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFFC3 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * QTENF(L)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
        ACR(L) = APAR(L) / 2.19E5

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
        KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))

        WCARB(L) = VCM(L) * (CI(L) - CCP(L))
     &           / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))

        WLITE(L) = ALPHA(FT) * ACR(L) * (CI(L) - CCP(L))
     &           / (CI(L) + 2 * CCP(L))

        WEXPT(L) = 0.5 * VCM(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------

CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
        AL(L) = (WL(L) - RD(L)) * FSMC(L)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = RATIO * GLCO2(L)

      ENDDO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L).LE.GLMIN(FT) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT)
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      RETURN
      END

!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1991) C4 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
! Written by Peter Cox (February 1996)
!**********************************************************************
      SUBROUTINE LEAF_C4 (LAND_FIELD,VEG_PTS,VEG_INDEX,FT
     &,                   DQ,APAR,TL,CA,OA,PSTAR,FSMC
     &,                   GL,AL,CI,RD)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points
!                                 !    on the land grid.
     &,FT                         ! IN Plant functional type.

      REAL
     & DQ(LAND_FIELD)             ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)           ! IN Absorbed PAR (W/m2)
     &,TL(LAND_FIELD)             ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)             ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)             ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_FIELD)          ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_FIELD)           ! IN Soil water factor.
     &,GL(LAND_FIELD)             ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_FIELD)             ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,RD(LAND_FIELD)             ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_FIELD)             ! OUT Internal CO2 pressure (Pa).
     &,ACR(LAND_FIELD)            ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,B1(LAND_FIELD)             !
     &,B2(LAND_FIELD)             !
     &,B3(LAND_FIELD)             ! WORK Coefficients of the quadratic.
     &,CCP(LAND_FIELD)            ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_FIELD)           ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_FIELD)          ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)          ! WORK Leaf conductance for CO2 (m/s).
     &,QTENF(LAND_FIELD)          ! WORK Q10 function.
     &,TDEGC(LAND_FIELD)          ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)            ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)          ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).
     &,WL(LAND_FIELD)             ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)          ! WORK Carboxylation,
     &,WLITE(LAND_FIELD)          !      Light, and
     &,WEXPT(LAND_FIELD)          !      export limited gross
!                                 !      photosynthetic rates
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_FIELD)             ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).

      INTEGER
     & CLOS_INDEX(LAND_FIELD)     ! WORK Index of land points
!                                 !      with closed stomata.
     &,CLOS_PTS                   ! WORK Number of land points
!                                 !      with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)     ! WORK Index of land points
!                                 !      with open stomata.
     &,OPEN_PTS                   ! WORK Number of land points
!                                 !      with open stomata.
      INTEGER
     & J,L                        ! WORK Loop counters.

      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice
C-----------------------------------------------------------------------
C Functional Type dependent parameters
C-----------------------------------------------------------------------
      INTEGER
     + C3(NPFT)                   ! 1 for C3 Plants, 0 for C4 Plants.

      REAL
     + ALPHA(NPFT)                ! Quantum efficiency
C                                 ! (mol CO2/mol PAR photons).
     +,A_WL(NPFT)                 ! Allometric coefficient relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index (kg C/m2).
     +,A_WS(NPFT)                 ! Woody biomass as a multiple of
C                                 ! live stem biomass.
     +,B_WL(NPFT)                 ! Allometric exponent relating
C                                 ! the target woody biomass to
C                                 ! the leaf area index.
     +,DGL_DM(NPFT)               ! Rate of change of leaf turnover
C                                 ! rate with moisture availability.
     +,DGL_DT(NPFT)               ! Rate of change of leaf turnover
C                                 ! rate with temperature (/K)
     +,DQCRIT(NPFT)               ! Critical humidity deficit
C                                 ! (kg H2O/kg air).
     +,ETA_SL(NPFT)               ! Live stemwood coefficient
C                                 ! (kg C/m/LAI).
     +,FSMC_OF(NPFT)              ! Moisture availability below
C                                 ! which leaves are dropped.
     +,F0(NPFT)                   ! CI/CA for DQ = 0.
     +,GLMIN(NPFT)                ! Minimum leaf conductance for H2O
     +,G_AREA(NPFT)               ! Disturbance rate (/360days).
     +,G_GROW(NPFT)               ! Rate of leaf growth (/360days).
     +,G_LEAF_0(NPFT)             ! Minimum turnover rate for leaves
!                                 ! (/360days).
     +,G_ROOT(NPFT)               ! Turnover rate for root biomass
!                                 ! (/360days).
     +,G_WOOD(NPFT)               ! Turnover rate for woody biomass
!                                 ! (/360days).
     +,KPAR(NPFT)                 ! PAR Extinction coefficient
C                                 ! (m2 leaf/m2 ground).
     +,LAI_MAX(NPFT)              ! Maximum projected LAI.
     +,LAI_MIN(NPFT)              ! Minimum projected LAI.
     +,NL0(NPFT)                  ! Top leaf nitrogen concentration
C                                 ! (kg N/kg C).
     +,NR_NL(NPFT)                ! Ratio of root nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,NS_NL(NPFT)                ! Ratio of stem nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,OMEGA(NPFT)                ! Leaf scattering coefficient
C                                 ! for PAR.
     +,R_GROW(NPFT)               ! Growth respiration fraction.
     +,SIGL(NPFT)                 ! Specific density of leaf carbon
C                                 ! (kg C/m2 leaf).
     +,TLEAF_OF(NPFT)             ! Temperature below which leaves are
C                                 ! dropped.
     +,TLOW(NPFT)                 ! Lower temperature for
C                                 ! photosynthesis (deg C)
     +,TUPP(NPFT)                 ! Upper temperature for
C                                 ! photosynthesis (deg C)

C----------------------------------------------------------------------
C                        BT     NT    C3G    C4G     S
C----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0,     1 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040,  0.08 /
      DATA A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10 /               
      DATA A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 /
      DATA B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667 /
      DATA DGL_DM  /   10.0,  10.0,   0.0,   0.0,  10.0 /               
      DATA DGL_DT  /    9.0,   0.0,   0.0,   0.0,   0.0 /
      DATA DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100 /
      DATA ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01 /
      DATA F0      /  0.875, 0.875, 0.900, 0.800, 0.900 /
      DATA FSMC_OF /   0.20,  0.20,  0.20,  0.20,  0.20 /               
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6 /
      DATA G_AREA  /  0.004, 0.004,  0.10,  0.10,  0.05 /
      DATA G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00 /
      DATA G_LEAF_0/   0.20,  0.20,  0.20,  0.20,  0.20 /               
      DATA G_ROOT  /   0.20,  0.20,  0.20,  0.20,  0.20 /
      DATA G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05 /               
      DATA KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50 /
      DATA LAI_MAX /  10.00, 10.00,  4.00,  4.00,  6.00 /
      DATA LAI_MIN /   4.00,  4.00,  1.00,  1.00,  1.00 /
      DATA NL0     /  0.050, 0.030, 0.060, 0.030, 0.030 /
      DATA NR_NL   /   1.00,  1.00,  1.00,  1.00,  1.00 /               
      DATA NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10 /
      DATA OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 /
      DATA R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA SIGL    / 0.0300,0.1000,0.0250,0.0500,0.0500 /
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,223.15 /               
      DATA TLOW    /    0.0,  -5.0,   0.0,  13.0,   0.0 /               
      DATA TUPP    /   36.0,  31.0,  36.0,  45.0,  36.0 /               

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL
     & BETA1, BETA2   ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4     ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4 ! Constant relating VCMAX and leaf N
!                     ! from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                     ! - assuming dry matter is 40% carbon by mass)
!                     ! and Jacobs 1994:
!                     ! C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                     ! (mol/m2/s)
     &,R              ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
!                     ! resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93
     &,          FDC3 = 0.015,  FDC4 = 0.025
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3
     &,          R = 8.3144  , RATIO = 1.6
     &,          ZERODEGC = 273.15)

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
        IF (FSMC(L).EQ.0.0 .OR. DQ(L).GE.DQCRIT(FT)
     &                     .OR. APAR(L).EQ.0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFFC4 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        CCP(L) = 0.0

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * QTENF(L)

      ENDDO

CDIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFFC4 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * QTENF(L)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
        ACR(L) = APAR(L) / 2.19E5

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        WCARB(L) = VCM(L)

        WLITE(L) = ALPHA(FT) * ACR(L)

        WEXPT(L) = 20000.0 * VCM(L) * CI(L) / PSTAR(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------

CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
        AL(L) = (WL(L) - RD(L)) * FSMC(L)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = GLCO2(L) * RATIO

      ENDDO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L).LE.GLMIN(FT) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT)
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      RETURN
      END
