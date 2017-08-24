C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C Routine to calculate the bulk stomatal resistance and the canopy
C CO2 fluxes
C
C Written by Peter Cox (Nov 1995)
C**********************************************************************
      SUBROUTINE SF_STOM  (LAND_PTS,LAND_INDEX,P1,P_POINTS
     +,                    FT,CO2,HT,IPAR,LAI,NL0,PSTAR
     +,                    Q1,RA,ROOT,TSTAR,V_CRIT,V_ROOT,V_WILT
     +,                    VEGF,GPP,NPP,RESP_P,GC,LTIMER,FSMC)


      IMPLICIT NONE

      INTEGER
     + LAND_PTS                   ! IN Number of land points to be
C                                 !    processed.
     +,LAND_INDEX(LAND_PTS)       ! IN Index of land points.
     +,P1                         ! IN First P point to be processed.
     +,P_POINTS                   ! IN Number of P points to be
C                                 !    processed.

      INTEGER
     + FT(LAND_PTS)               ! IN Plant functional type.

      REAL
     + CO2                        ! IN Atmospheric CO2 concentration
C                                 !    (kg CO2/kg air).
     +,HT(LAND_PTS)               ! IN Canopy height (m).
     +,IPAR(P_POINTS)             ! IN Incident PAR (W/m2).
     +,LAI(LAND_PTS)              ! IN Leaf area index.
     +,NL0(LAND_PTS)              ! IN Nitrogen concentration
C                                 !    of top leaf (kg N/kg C).
     +,PSTAR(P_POINTS)            ! IN Surface pressure (Pa).
     +,Q1(P_POINTS)               ! IN Specific humidity of level 1
C                                 !    (kg H2O/kg air).
     +,RA(P_POINTS)               ! IN Aerodynamic resistance (s/m).
     +,ROOT(LAND_PTS)             ! IN Root biomass (kg C/m2).
     +,TSTAR(P_POINTS)            ! IN Surface temperature (K).
     +,V_CRIT(LAND_PTS)           ! IN Volumetric soil moisture
C                                 !    concentration above which
C                                 !    stomata are not sensitive
C                                 !    to soil water (m3 H2O/m3 soil).
     +,V_ROOT(LAND_PTS)           ! IN Volumetric soil moisture
C                                 !    concentration in the rootzone
C                                 !    (m3 H2O/m3 soil).
     +,V_WILT(LAND_PTS)           ! IN Volumetric soil moisture
C                                 !    concentration below which
C                                 !    stomata close (m3 H2O/m3 soil).
     +,VEGF(LAND_PTS)             ! IN Vegetated fraction.
     +,GPP(LAND_PTS)              ! OUT Gross Primary Productivity
C                                 !     (kg C/m2/s).
     +,NPP(LAND_PTS)              ! OUT Net Primary Productivity
C                                 !     (kg C/m2/s).
     +,RESP_P(LAND_PTS)           ! OUT Plant respiration rate
C                                 !     (kg C/m2/sec).
     +,GC(LAND_PTS)               ! INOUT Canopy resistance to H2O
C                                 !       (m/s).

      REAL
     + ANETC(LAND_PTS)            ! WORK Net canopy photosynthesis
C                                 !     (mol CO2/m2/s).
     +,CO2C(LAND_PTS)             ! WORK Canopy level CO2 concentration
C                                 !      (kg CO2/kg air).
     +,CI(LAND_PTS)               ! WORK Internal CO2 pressure (Pa).
     +,DQ(P_POINTS)               ! WORK Specific humidity deficit
C                                 !      (kg H2O/kg air).
     +,DQC(LAND_PTS)              ! WORK Canopy level specific humidity
C                                 !      deficit (kg H2O/kg air).
     +,FPAR(LAND_PTS)             ! WORK PAR absorption factor.
     +,FSMC(LAND_PTS)             ! OUT Soil water factor.
     +,NL(LAND_PTS)               ! WORK Mean leaf nitrogen
C                                 !      concentration (kg N/kg C).
     +,N_LEAF(LAND_PTS)           ! WORK Nitrogen contents of the leaf,
     +,N_ROOT(LAND_PTS)           !      root,
     +,N_STEM(LAND_PTS)           !      and stem (kg N/m2).
     +,QS(P_POINTS)               ! WORK Saturated specific humidity
C                                 !      (kg H2O/kg air).
     +,RA_RC(LAND_PTS)            ! WORK Ratio of aerodynamic resistance
C                                 !      to canopy resistance.
     +,RDC(LAND_PTS)              ! WORK Canopy dark respiration,
C                                 !      without soil water dependence
C                                 !      (mol CO2/m2/s).
     +,RESP_P_G(LAND_PTS)         ! WORK Plant growth respiration rate
C                                 !      (kg C/m2/sec).
     +,RESP_P_M(LAND_PTS)         ! WORK Plant maintenance respiration
C                                 !      rate (kg C/m2/sec).
     +,RHOSTAR(P_POINTS)          ! WORK Surface air density (kg/m3).

      LOGICAL
     + LTIMER

      INTEGER
     + I,J,K,L                    ! WORK Loop counters.
     +,VEG_PTS                    ! WORK Number of vegetated points.
     +,VEG_INDEX(LAND_PTS)        ! WORK Index of vegetated points
C                                 !      on the land grid.

C-----------------------------------------------------------------------
C Parameters
C-----------------------------------------------------------------------
      REAL
     + RAIR                       ! Gas constant for dry air (J/kg/K).
      PARAMETER (RAIR = 287.05)

      REAL
     + O2                         ! Atmospheric concentration of
C                                 ! oxygen (kg O2/kg air).
      PARAMETER (O2 = 0.23)

      INTEGER
     + ITER                       ! Number of iterations to
C                                 ! determine the canopy climate.
      PARAMETER (ITER = 2)

C-----------------------------------------------------------------------
C Functional Type dependent parameters
C-----------------------------------------------------------------------
      REAL
     + ETA_SL(4)                  ! Live stemwood coefficient
C                                 ! (kg C/m/LAI).
     +,KPAR(4)                    ! PAR Extinction coefficient.
     +,NR_NL(4)                   ! Ratio of root nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,NS_NL(4)                   ! Ratio of stem nitrogen
C                                 ! concentration to leaf
C                                 ! nitrogen concentration.
     +,R_GROW(4)                  ! Growth respiration fraction.
     +,SIGL(4)                    ! Specific leaf density
C                                 ! (kg C/projected LAI).
C----------------------------------------------------------------------
C                       BT    NT   C3G   C4G
C----------------------------------------------------------------------
      DATA ETA_SL  /  0.01, 0.01, 0.01, 0.01 /  ! Friend et al. (1995)
      DATA KPAR    /  0.50, 0.50, 0.50, 0.50 /  ! Friend et al. (1995)
      DATA NR_NL   /  1.00, 1.00, 1.00, 1.00 /  !
      DATA NS_NL   /  0.04, 0.10, 1.00, 1.00 /  !
      DATA R_GROW  /  0.25, 0.25, 0.25, 0.25 /  ! Bonan (1995)
      DATA SIGL    /  0.04, 0.10, 0.04, 0.04 /  ! Schulze et al. (1994)


      IF (LTIMER) THEN
        CALL TIMER('SFSTOM  ',103)
      ENDIF

C-----------------------------------------------------------------------
C Create index of vegetated points on the land grid
C-----------------------------------------------------------------------
      VEG_PTS=0
      DO L=1,LAND_PTS
        I = LAND_INDEX(L)-(P1-1)
        IF (VEGF(L).GT.0.0 .AND. LAI(L).GT.0.0) THEN
          VEG_PTS = VEG_PTS + 1
          VEG_INDEX(VEG_PTS) = L
        ELSE
          GC(L) = 0.0
          NPP(L) = 0.0
          GPP(L) = 0.0
          RESP_P(L) = 0.0
        ENDIF
      ENDDO

C-----------------------------------------------------------------------
C Calculate the surface to level 1 humidity deficit and the surface
C density of the air
C-----------------------------------------------------------------------
      CALL QSAT(QS,TSTAR,PSTAR,P_POINTS)
      DO I=1,P_POINTS
        DQ(I) = MAX(0.0,(QS(I) - Q1(I)))
        RHOSTAR(I) = PSTAR(I) / (RAIR * TSTAR(I))
      ENDDO

C-----------------------------------------------------------------------
C Calculate the soil water factor (Cox and Huntingford, 1995)
C-----------------------------------------------------------------------
      DO L=1,LAND_PTS

        IF (V_ROOT(L) .GT. V_CRIT(L)) THEN
          FSMC(L) = 1.0
        ELSEIF (V_ROOT(L) .LE. V_WILT(L)) THEN
          FSMC(L) = 0.0
        ELSE
          FSMC(L) = (V_ROOT(L) - V_WILT(L))
     &            / (V_CRIT(L) - V_WILT(L))
        ENDIF

      ENDDO

C-----------------------------------------------------------------------
C Calculate the PAR absorption factor
C-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        FPAR(L) = (1 - EXP(-KPAR(FT(L))*LAI(L))) / KPAR(FT(L))

      ENDDO


C-----------------------------------------------------------------------
C Iterate to ensure that the canopy humidity deficit is consistent with
C the H2O flux. Ignore the (small) difference between the canopy and
C reference level CO2 concentration. Intially set the canopy humidity
C deficit using the previous value of GC.
C-----------------------------------------------------------------------
      DO K=1,ITER

C-----------------------------------------------------------------------
C Diagnose the canopy level humidity deficit and CO2 concentration
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          I = LAND_INDEX(L) - (P1-1)

          RA_RC(L) = RA(I) * GC(L)
          DQC(L) = DQ(I) / (1 + RA_RC(L))
          CO2C(L) = CO2

        ENDDO

C-----------------------------------------------------------------------
C Call CANOPY to calculate the canopy resistance and photosynthesis
C-----------------------------------------------------------------------
        CALL CANOPY (LAND_PTS,LAND_INDEX,P1,P_POINTS
     +,              VEG_PTS,VEG_INDEX
     +,              FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR
     +,              NL0,FPAR,FSMC,LAI
     +,              GC,ANETC,CI,RDC,LTIMER)

      ENDDO ! LOOP OVER ITER

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

C-----------------------------------------------------------------------
C Calculate the mean leaf nitrogen concentration assuming perfect
C light acclimation
C-----------------------------------------------------------------------
        NL(L) = (FPAR(L) / LAI(L)) * NL0(L)

C-----------------------------------------------------------------------
C Calculate the total nitrogen content of the leaf, root and stem
C-----------------------------------------------------------------------
        N_LEAF(L) = NL(L) * SIGL(FT(L)) * LAI(L)
        N_ROOT(L) = NR_NL(FT(L)) * NL(L) * ROOT(L)
        N_STEM(L) = NS_NL(FT(L)) * NL(L)
     &            * ETA_SL(FT(L)) * HT(L) * LAI(L)

C-----------------------------------------------------------------------
C Calculate the Gross Primary Productivity and the plant maintenance
C respiration rate in kg C/m2/sec
C-----------------------------------------------------------------------
        GPP(L) = 12.0E-3 * (ANETC(L) + RDC(L)*FSMC(L))
        RESP_P_M(L) = 12.0E-3 * RDC(L)
     &     * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)

C-----------------------------------------------------------------------
C Calculate the total plant respiration and the Net Primary Productivity
C-----------------------------------------------------------------------
        RESP_P_G(L) = R_GROW(FT(L)) * (GPP(L) - RESP_P_M(L))
        RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
        NPP(L) = GPP(L) - RESP_P(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SFSTOM  ',104)
      ENDIF

      RETURN
      END

C***********************************************************************
C Calculates the canopy resistance, net photosynthesis and transpiration
C by scaling-up the leaf level response using the "Big-Leaf" approach
C of Sellers et al. (1994)
C
C Written by Peter Cox (May 1995)
C***********************************************************************
      SUBROUTINE CANOPY (LAND_PTS,LAND_INDEX,P1,P_POINTS
     +,                  VEG_PTS,VEG_INDEX
     +,                  FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR,NL0
     +,                  FPAR,FSMC,LAI
     +,                  GC,ANETC,CI,RDC,LTIMER)

      IMPLICIT NONE

      INTEGER
     + LAND_PTS                   ! IN Number of land points to be
C                                 !    processed.
     +,LAND_INDEX(LAND_PTS)       ! IN Index of land points.
     +,P1                         ! IN First P point to be processed.
     +,P_POINTS                   ! IN Number of P points to be
C                                 !    processed.
     +,VEG_PTS                    ! IN Number of vegetated points.
     +,VEG_INDEX(LAND_PTS)        ! IN Index of vegetated points
C                                 !    on the land grid.

      INTEGER
     + FT(LAND_PTS)               ! IN Plant functional type.

      REAL
     + CO2C(LAND_PTS)             ! IN Canopy level CO2 concentration
C                                 !    (kg CO2/kg air).
     +,DQC(LAND_PTS)              ! IN Canopy level specific humidity
C                                 !    deficit (kg H2O/kg air).
     +,O2                         ! IN Atmospheric O2 concentration
C                                 !    (kg O2/kg air).
     +,PSTAR(P_POINTS)            ! IN Surface pressure (Pa).
     +,IPAR(P_POINTS)             ! IN Incident PAR (W/m2).
     +,TSTAR(P_POINTS)            ! IN Surface temperature (K).
     +,NL0(LAND_PTS)              ! IN Nitrogen concentration of
C                                 !    top leaf (kg N/kg C).
     +,FPAR(LAND_PTS)             ! IN PAR absorption factor.
     +,FSMC(LAND_PTS)             ! IN Soil water factor.
     +,LAI(LAND_PTS)              ! IN Leaf area index
C                                 !    (m2 leaf/m2 ground).


      REAL
     + ANETC(LAND_PTS)            ! OUT Net canopy photosynthesis
C                                 !     (mol CO2/m2/s).
     +,CI(LAND_PTS)               ! OUT Internal CO2 concentration
C                                 !     (mol CO2/m3).
     +,GC(LAND_PTS)               ! OUT Canopy conductance for H2O
C                                 !     (m/s).
     +,RDC(LAND_PTS)              ! OUT Canopy dark respiration
C                                 !     (mol CO2/m2/s).
     +,ANETL(LAND_PTS)            ! WORK Net leaf photosynthesis
C                                 !      (mol CO2/m2/s/LAI).
     +,APAR(LAND_PTS)             ! WORK PAR absorbed by the top leaf
C                                 !      (W/m2).
     +,CA(LAND_PTS)               ! WORK Canopy level CO2 pressure
C                                 !      (Pa).
     +,DQM(LAND_PTS)              ! WORK Canopy level humidity
C                                 !      deficit (mol H2O/m3).
     +,GL(LAND_PTS)               ! WORK Leaf conductance for H2O
C                                 !      (m/s).
     +,OA(LAND_PTS)               ! WORK Atmospheric O2 pressure
C                                 !      (Pa).
     +,RD(LAND_PTS)               ! WORK Dark respiration of top leaf
C                                 !      (mol CO2/m2/s).

      LOGICAL
     + LTIMER

      INTEGER
     + I,J,L                      ! WORK Loop counters.

C-----------------------------------------------------------------------
C Functional Type dependent parameters
C-----------------------------------------------------------------------
      REAL
     + OMEGA(4)                   ! Leaf scattering coefficient for PAR.
C-----------------------------------------------------------------------
C                       BT    NT   C3G   C4G
C-----------------------------------------------------------------------
      DATA OMEGA   /  0.15, 0.15, 0.15, 0.17 /  ! Sellers et al. (1994)

C-----------------------------------------------------------------------
C Parameters
C-----------------------------------------------------------------------
      REAL
     + R                          ! Gas constant (J/K/mol)
      PARAMETER (R = 8.3144)

      REAL
     + EPSILON                    ! Ratio of molecular weights of water
C                                 ! and dry air.
     +,EPCO2                      ! Ratio of molecular weights of CO2
C                                 ! and dry air.
     +,EPO2                       ! Ratio of molecular weights of O2
C                                 ! and dry air.
      PARAMETER (EPSILON = 0.62198, EPCO2 = 1.5194, EPO2 = 1.106)

      IF (LTIMER) THEN
        CALL TIMER('CANOPY  ',103)
      ENDIF

C-----------------------------------------------------------------------
C Calculate the atmospheric pressures of CO2 and O2
C-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L) - (P1-1)

        CA(L) = CO2C(L) / EPCO2 * PSTAR(I)
        OA(L) = O2 / EPO2 * PSTAR(I)
        DQM(L) = DQC(L) / EPSILON * PSTAR(I) / (R * TSTAR(I))

C-----------------------------------------------------------------------
C Calculate the PAR absorbed by the top leaf
C-----------------------------------------------------------------------
        APAR(L) = (1 - OMEGA(FT(L))) * IPAR(I)

      ENDDO

C-----------------------------------------------------------------------
C Call the leaf level model for the top leaf of the C3 and C4 plants
C-----------------------------------------------------------------------
      CALL LEAF_C3 (LAND_PTS,LAND_INDEX,P1,P_POINTS
     +,             VEG_PTS,VEG_INDEX
     +,             FT,DQC,APAR,TSTAR,CA,OA,PSTAR
     +,             NL0,FSMC
     +,             GL,ANETL,CI,RD,LTIMER)

      CALL LEAF_C4 (LAND_PTS,LAND_INDEX,P1,P_POINTS
     +,             VEG_PTS,VEG_INDEX
     +,             FT,DQC,APAR,TSTAR,CA,OA,PSTAR
     +,             NL0,FSMC
     +,             GL,ANETL,CI,RD,LTIMER)

C-----------------------------------------------------------------------
C Scale-up to the canopy level
C-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        ANETC(L) = ANETL(L) * FPAR(L)
        GC(L) = FPAR(L) * GL(L)
        RDC(L) = RD(L) * FPAR(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('CANOPY  ',104)
      ENDIF

      RETURN

      END
