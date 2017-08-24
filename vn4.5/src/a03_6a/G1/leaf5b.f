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

!!!********************************************************************
!!! Calculates the leaf resistance and net photosynthesis using:
!!!  (i) Collatz et al. (1992) C3 photosynthesis model
!!! (ii) Jacobs (1994) CI/CA closure.
!!!
!!! Written by Peter Cox (February 1996)
! 4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!********************************************************************
      SUBROUTINE LEAF_C3 (LAND_PTS
     &,                   LAND_INDEX,P1
     &,                   LAND_FIELD,P_FIELD
     &,                   VEG_PTS,VEG_INDEX
     &,                   FT,DQ,APAR,TL,CA,OA,PSTAR
     &,                   NL0,FSMC
     &,                   GL,AL,CI,RD,LTIMER)

      IMPLICIT NONE

      INTEGER
     & LAND_PTS               ! IN Number of land points to be
!                                  processed.
     &,LAND_FIELD             ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)   ! IN Index of land points.
     &,P1                     ! IN First P point to be processed.
     &,P_FIELD                ! IN Total number of P points.
     &,VEG_PTS                ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)  ! IN Index of vegetated points
!                                  on the land grid.

      INTEGER
     & FT(LAND_FIELD)         ! IN Plant functional type.

      REAL
     & DQ(LAND_FIELD)         ! IN Canopy level specific humidity
!                                  deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)       ! IN Absorbed PAR (W/m2)
     &,TL(P_FIELD)            ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)         ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)         ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(P_FIELD)         ! IN Atmospheric pressure (Pa).
     &,NL0(LAND_FIELD)        ! IN Leaf nitrogen conncentration
!                                 (kg N/kg C).
     &,FSMC(LAND_FIELD)       ! IN Soil water factor.

! OUTPUT
      REAL
     & GL(LAND_FIELD)         ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_FIELD)         ! OUT Net Leaf photosynthesis
!                                   (mol CO2/m2/s).
     &,RD(LAND_FIELD)         ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_FIELD)         ! OUT Internal CO2 pressure (Pa).


! WORKspace
      REAL
     & ACR(LAND_FIELD)        ! WORK Absorbed PAR (mol photons/m2/s).
     &,B1(LAND_FIELD)         ! WORK Coefficients of the quadratic.
     &,B2(LAND_FIELD)         !                    "
     &,B3(LAND_FIELD)         !                    "
     &,CCP(LAND_FIELD)        ! WORK Photorespiratory compensatory
!                                  point (mol/m3).
     &,CONV(LAND_FIELD)       ! WORK Factor for converting mol/m3
!                                  into Pa (J/mol).
     &,DENOM(LAND_FIELD)      ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)      ! WORK Leaf conductance for CO2 (m/s).
     &,KC(LAND_FIELD)         ! WORK Michaelis constant for CO2 (Pa)
     &,KO(LAND_FIELD)         ! WORK Michaelis constant for O2 (Pa).
     &,QTENF(LAND_FIELD)      ! WORK Q10 function.
     &,TAU(LAND_FIELD)        ! WORK CO2/O2 specificity ratio.
     &,TDEGC(LAND_FIELD)      ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)        ! WORK Maximum rate of carboxylation
!                                    of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)      ! WORK Maximum rate of carboxylation
!                                    of Rubisco - without the
!                                    temperature factor
!                                    (mol CO2/m2/s).
     &,WL(LAND_FIELD)         ! WORK Gross leaf phtosynthesis
!                                    (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)      ! WORK Carboxylation,
     &,WLITE(LAND_FIELD)      !      Light, and
     &,WEXPT(LAND_FIELD)      !      export limited gross
!                                    photosynthetic rates
!                                    (mol CO2/m2/s).
     &,WP(LAND_FIELD)         ! WORK Smoothed minimum of
!                                    Carboxylation and Light
!                                    limited gross photosynthesis
!                                    (mol CO2/m2/s).

      INTEGER
     & CLOS_INDEX(LAND_FIELD) ! WORK Index of land points
!                                    with closed stomata.
     &,CLOS_PTS               ! WORK Number of land points
!                                    with closed stomata.
     &,OPEN_INDEX(LAND_FIELD) ! WORK Index of land points
!                                    with open stomata.
     &,OPEN_PTS               ! WORK Number of land points
!                                    with open stomata.

      LOGICAL
     & LTIMER

      INTEGER
     & I,J,L                  ! WORK Loop counters.

!-----------------------------------------------------------------------
! Functional Type dependent parameters
!-----------------------------------------------------------------------
      INTEGER
     & C3(4)                  ! 1 for C3 Plants, 0 for C4 Plants.

      REAL
     & ALPHA(4)               ! Quantum efficiency
!                               (mol CO2/mol PAR photons).
     &,DQCRIT(4)              ! Critical humidity deficit
!                               (kg H2O/kg air).
     &,F0(4)                  ! CI/CA for DQ = 0.
     &,GLMIN(4)               ! Minimum leaf conductance for H2O
!                               (m/s).
!----------------------------------------------------------------------
!                        BT     NT    C3G    C4G
!----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040 /
      DATA DQCRIT  /  0.090, 0.060, 0.150, 0.075 /
      DATA F0      /  0.875, 0.875, 0.925, 0.800 /
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6 /

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL
     & BETA1, BETA2   ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4     ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4 ! Constant relating VCMAX and leaf N
!                       from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                       - assuming dry matter is 40% carbon by mass)
!                       and Jacobs 1994:
!                       C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                       (mol/m2/s)
     &,R              ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
!                       resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93
     &,          FDC3 = 0.015,  FDC4 = 0.025
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3
     &,          R = 8.3144  , RATIO = 1.6
     &,          ZERODEGC = 273.15)

      IF (LTIMER) THEN
        CALL TIMER('LEAF_C3  ',103)
      ENDIF

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
!----------------------------------------------------------------------
! Only carry out calculations for points with C3 plants
!----------------------------------------------------------------------
        IF (C3(FT(L)).EQ.1) THEN

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
          IF (FSMC(L).EQ.0.0 .OR. DQ(L).GE.DQCRIT(FT(L))
     &                       .OR. APAR(L).LE.0.0) THEN
            CLOS_PTS = CLOS_PTS + 1
            CLOS_INDEX(CLOS_PTS) = J
          ELSE
            OPEN_PTS = OPEN_PTS + 1
            OPEN_INDEX(OPEN_PTS) = J
          ENDIF

        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))
        I = LAND_INDEX(L)

        VCMAX(L) = NEFFC3 * NL0(L)
        TDEGC(L) = TL(I) - ZERODEGC

        TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
        CCP(L) = 0.5 * OA(L) / TAU(L)


        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - 36.0)))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * VCM(L)

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(I)

      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))
        I = LAND_INDEX(L)

        VCMAX(L) = NEFFC3 * NL0(L)
        TDEGC(L) = TL(I) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - 36.0)))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * VCM(L)

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(I)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT(L))
     &        * (1 - DQ(L) / DQCRIT(FT(L))) + CCP(L)

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
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
        KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))

        WCARB(L) = VCM(L) * (CI(L) - CCP(L))
     &           / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))

        WLITE(L) = ALPHA(FT(L)) * ACR(L) * (CI(L) - CCP(L))
     &           / (CI(L) + 2 * CCP(L))

        WEXPT(L) = 0.5 * VCM(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))
     .         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))
     .         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
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
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L).LE.GLMIN(FT(L)) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT(L))
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT(L))
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('LEAF_C3  ',104)
      ENDIF

      RETURN
      END

!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1991) C4 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
! Written by Peter Cox (February 1996)
!**********************************************************************
      SUBROUTINE LEAF_C4 (LAND_PTS
     &,                   LAND_INDEX,P1
     &,                   LAND_FIELD,P_FIELD
     &,                   VEG_PTS,VEG_INDEX
     &,                   FT,DQ,APAR,TL,CA,OA,PSTAR
     &,                   NL0,FSMC
     &,                   GL,AL,CI,RD,LTIMER)

      IMPLICIT NONE

      INTEGER
     & LAND_PTS               ! IN Number of land points to be
!                                  processed.
     &,LAND_FIELD             ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)   ! IN Index of land points.
     &,P1                     ! IN First P point to be processed.
     &,P_FIELD                ! IN Total number of P points.
     &,VEG_PTS                ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)  ! IN Index of vegetated points
!                                  on the land grid.

      INTEGER
     & FT(LAND_FIELD)         ! IN Plant functional type.

      REAL
     & DQ(LAND_FIELD)         ! IN Canopy level specific humidity
!                                  deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)       ! IN Absorbed PAR (W/m2)
     &,TL(P_FIELD)            ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)         ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)         ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(P_FIELD)         ! IN Atmospheric pressure (Pa).
     &,NL0(LAND_FIELD)        ! IN Leaf nitrogen conncentration
!                                  (kg N/kg C).
     &,FSMC(LAND_FIELD)       ! IN Soil water factor.

! OUTPUT
      REAL
     & GL(LAND_FIELD)         ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_FIELD)         ! OUT Net Leaf photosynthesis
!                                   (mol CO2/m2/s).
     &,RD(LAND_FIELD)         ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_FIELD)         ! OUT Internal CO2 pressure (Pa).

!WORKspace
      REAL
     & ACR(LAND_FIELD)        ! WORK Absorbed PAR (mol photons/m2/s).
     &,B1(LAND_FIELD)         ! WORK Coefficients of the quadratic.
     &,B2(LAND_FIELD)         !                   "
     &,B3(LAND_FIELD)         !                   "
     &,CCP(LAND_FIELD)        ! WORK Photorespiratory compensatory
!                                    point (mol/m3).
     &,CONV(LAND_FIELD)       ! WORK Factor for converting mol/m3
!                                    into Pa (J/mol).
     &,DENOM(LAND_FIELD)      ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)      ! WORK Leaf conductance for CO2 (m/s).
     &,QTENF(LAND_FIELD)      ! WORK Q10 function.
     &,TDEGC(LAND_FIELD)      ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)        ! WORK Maximum rate of carboxylation
!                                    of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)      ! WORK Maximum rate of carboxylation
!                                    of Rubisco - without the
!                                    temperature factor
!                                    (mol CO2/m2/s).
     &,WL(LAND_FIELD)         ! WORK Gross leaf phtosynthesis
!                                    (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)      ! WORK Carboxylation,
     &,WLITE(LAND_FIELD)      !      Light, and
     &,WEXPT(LAND_FIELD)      !      export limited gross
!                                    photosynthetic rates
!                                    (mol CO2/m2/s).
     &,WP(LAND_FIELD)         ! WORK Smoothed minimum of Carboxylation
!                                    and Light limited gross
!                                    photosynthesis (mol CO2/m2/s).

      INTEGER
     & CLOS_INDEX(LAND_FIELD) ! WORK Index of land points
!                                    with closed stomata.
     &,CLOS_PTS               ! WORK Number of land points
!                                    with closed stomata.
     &,OPEN_INDEX(LAND_FIELD) ! WORK Index of land points
!                                    with open stomata.
     &,OPEN_PTS               ! WORK Number of land points
!                                    with open stomata.

      LOGICAL
     & LTIMER

      INTEGER
     & I,J,L                  ! WORK Loop counters.

!-----------------------------------------------------------------------
! Functional Type dependent parameters
!-----------------------------------------------------------------------
      INTEGER
     & C3(4)                  ! 1 for C3 Plants, 0 for C4 Plants.

      REAL
     & ALPHA(4)               ! Quantum efficiency
!                               (mol CO2/mol PAR photons).
     &,DQCRIT(4)              ! Critical humidity deficit
!                               (kg H2O/kg air).
     &,F0(4)                  ! CI/CA for DQ = 0.
     &,GLMIN(4)               ! Minimum leaf conductance for H2O
!                               (m/s).
!----------------------------------------------------------------------
!                        BT     NT    C3G    C4G
!----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040 /
      DATA DQCRIT  /  0.090, 0.060, 0.150, 0.075 /
      DATA F0      /  0.875, 0.875, 0.925, 0.800 /
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6 /

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL
     & BETA1, BETA2   ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4     ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4 ! Constant relating VCMAX and leaf N
!                       from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                       - assuming dry matter is 40% carbon by mass)
!                       and Jacobs 1994:
!                       C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                       (mol/m2/s)
     &,R              ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
!                     ! resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93
     &,          FDC3 = 0.015,  FDC4 = 0.025
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3
     &,          R = 8.3144  , RATIO = 1.6
     &,          ZERODEGC = 273.15)

      IF (LTIMER) THEN
        CALL TIMER('LEAF_C4  ',103)
      ENDIF

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
!----------------------------------------------------------------------
! Only carry out calculations for points with C4 plants
!----------------------------------------------------------------------
        IF (C3(FT(L)).EQ.0) THEN

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
          IF (FSMC(L).EQ.0.0 .OR. DQ(L).GE.DQCRIT(FT(L))
     &                       .OR. APAR(L).LE.0.0) THEN
            CLOS_PTS = CLOS_PTS + 1
            CLOS_INDEX(CLOS_PTS) = J
          ELSE
            OPEN_PTS = OPEN_PTS + 1
            OPEN_INDEX(OPEN_PTS) = J
          ENDIF

        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))
        I = LAND_INDEX(L)

        VCMAX(L) = NEFFC4 * NL0(L)
        TDEGC(L) = TL(I) - ZERODEGC

        CCP(L) = 0.0

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - 36.0)))
     &           * (1 + EXP (0.3 * (13.0 - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * VCM(L)

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(I)

      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))
        I = LAND_INDEX(L)

        VCMAX(L) = NEFFC4 * NL0(L)
        TDEGC(L) = TL(I) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - 36.0)))
     &           * (1 + EXP (0.3 * (13.0 - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * VCM(L)

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(I)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT(L))
     &        * (1 - DQ(L) / DQCRIT(FT(L))) + CCP(L)

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
        I = LAND_INDEX(L)

        WCARB(L) = VCM(L)

        WLITE(L) = ALPHA(FT(L)) * ACR(L)

        WEXPT(L) = 20000.0 * VCM(L) * CI(L) / PSTAR(I)

      ENDDO

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))
     .         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))
     .         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
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
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L).LE.GLMIN(FT(L)) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT(L))
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT(L))
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('LEAF_C4  ',104)
      ENDIF

      RETURN
      END
