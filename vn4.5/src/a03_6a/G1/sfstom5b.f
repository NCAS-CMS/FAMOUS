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
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes
!
! Written by Peter Cox (Nov 1995)
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!**********************************************************************
      SUBROUTINE SF_STOM  (LAND_PTS,LAND_FIELD,LAND_MASK,P1,LAND1
     &,                   LAND_INDEX
     &,                   P_POINTS,P_FIELD
     &,                   FT,CO2,HT,IPAR,LAI,NL0,PSTAR
     &,                   Q1,RA,ROOT,TSTAR,V_CRIT,V_ROOT,V_WILT
     &,                   VEGF,GPP,NPP,RESP_P,GC,LTIMER,FSMC)


      IMPLICIT NONE

      INTEGER
     & LAND_PTS             ! IN Number of land points to be
!                                processed.
     &,LAND_FIELD           ! IN Total number of land points
     &,P_FIELD              ! IN Total number of p points
     &,LAND_INDEX(LAND_FIELD)
!                           ! IN Index of land points.
     &,P1                   ! IN First P point to be processed.
     &,LAND1                ! IN First P point to be processed.
     &,P_POINTS             ! IN Number of P points to be processed.

      LOGICAL
     & LAND_MASK(P_POINTS)  ! IN .TRUE. for land points

      INTEGER
     & FT(LAND_FIELD)       ! IN Plant functional type.

      REAL
     & CO2                  ! IN Atmospheric CO2 concentration
!                                (kg CO2/kg air).
     &,HT(LAND_FIELD)       ! IN Canopy height (m).
     &,IPAR(P_FIELD)        ! IN Incident PAR (W/m2).
     &,LAI(LAND_FIELD)      ! IN Leaf area index.
     &,NL0(LAND_FIELD)      ! IN Nitrogen concentration of top leaf
!                                (kg N/kg C).
     &,PSTAR(P_FIELD)       ! IN Surface pressure (Pa).
     &,Q1(P_FIELD)          ! IN Specific humidity of level 1
!                                (kg H2O/kg air).
     &,RA(P_FIELD)          ! IN Aerodynamic resistance (s/m).
     &,ROOT(LAND_FIELD)     ! IN Root biomass (kg C/m2).
     &,TSTAR(P_FIELD)       ! IN Surface temperature (K).
     &,V_CRIT(LAND_FIELD)   ! IN Volumetric soil moisture concentration
!                                above which stomata are not sensitive
!                                to soil water (m3 H2O/m3 soil).
     &,V_ROOT(LAND_FIELD)   ! IN Volumetric soil moisture concentration
!                                in the rootzone (m3 H2O/m3 soil).
     &,V_WILT(LAND_FIELD)   ! IN Volumetric soil moisture concentration
!                                below which stomata close
!                                (m3 H2O/m3 soil).
     &,VEGF(LAND_FIELD)     ! IN Vegetated fraction.


! OUTPUT
      REAL
     & GPP(LAND_FIELD)      ! OUT Gross Primary Productivity
!                                 (kg C/m2/s).
     &,NPP(LAND_FIELD)      ! OUT Net Primary Productivity (kg C/m2/s).
     &,RESP_P(LAND_FIELD)   ! OUT Plant respiration rate (kg C/m2/sec).
     &,GC(LAND_FIELD)       ! INOUT Canopy resistance to H2O (m/s).


! WORK
      REAL
     & ANETC(LAND_FIELD)    ! WORK Net canopy photosynthesis
!                                  (mol CO2/m2/s).
     &,CO2C(LAND_FIELD)     ! WORK Canopy level CO2 concentration
!                                  (kg CO2/kg air).
     &,CI(LAND_FIELD)       ! WORK Internal CO2 pressure (Pa).
     &,DQ(P_FIELD)          ! WORK Specific humidity deficit
!                                  (kg H2O/kg air).
     &,DQC(LAND_FIELD)      ! WORK Canopy level specific humidity
!                                  deficit (kg H2O/kg air).
     &,FPAR(LAND_FIELD)     ! WORK PAR absorption factor.
     &,FSMC(LAND_FIELD)     ! OUT Soil water factor.
     &,NL(LAND_FIELD)       ! WORK Mean leaf nitrogen
!                                  concentration (kg N/kg C).
     &,N_LEAF(LAND_FIELD)   ! WORK Nitrogen contents of the leaf,
     &,N_ROOT(LAND_FIELD)   !      root,
     &,N_STEM(LAND_FIELD)   !      and stem (kg N/m2).
     &,QS(P_FIELD)          ! WORK Saturated specific humidity
!                                  (kg H2O/kg air).
     &,RA_RC(LAND_FIELD)    ! WORK Ratio of aerodynamic resistance
!                                  to canopy resistance.
     &,RDC(LAND_FIELD)      ! WORK Canopy dark respiration, without
!                                  soil water dependence (mol CO2/m2/s).
     &,RESP_P_G(LAND_FIELD) ! WORK Plant growth respiration rate
!                                  (kg C/m2/sec).
     &,RESP_P_M(LAND_FIELD) ! WORK Plant maintenance respiration
!                                  rate (kg C/m2/sec).
     &,RHOSTAR(P_FIELD)     ! WORK Surface air density (kg/m3).

      LOGICAL
     & LTIMER

      INTEGER
     & I,J,K,L              ! WORK Loop counters.
     &,VEG_PTS              ! WORK Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)! WORK Index of vegetated points
!                                  on the land grid.

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL
     & RAIR                 ! Gas constant for dry air (J/kg/K).
      PARAMETER (RAIR = 287.05)

      REAL
     & O2                   ! Atmospheric concentration of
!                             oxygen (kg O2/kg air).
      PARAMETER (O2 = 0.23)

      INTEGER
     & ITER                 ! Number of iterations to
!                             determine the canopy climate.
      PARAMETER (ITER = 2)

!-----------------------------------------------------------------------
! Functional Type dependent parameters
!-----------------------------------------------------------------------
      REAL
     & ETA_SL(4)            ! Live stemwood coefficient
!                             (kg C/m/LAI).
     &,KPAR(4)              ! PAR Extinction coefficient.
     &,NR_NL(4)             ! Ratio of root nitrogen
!                             concentration to leaf
!                             nitrogen concentration.
     &,NS_NL(4)             ! Ratio of stem nitrogen
!                             concentration to leaf
!                             nitrogen concentration.
     &,R_GROW(4)            ! Growth respiration fraction.
     &,SIGL(4)              ! Specific leaf density
!                             (kg C/projected LAI).

!----------------------------------------------------------------------
!                       BT    NT   C3G   C4G
!----------------------------------------------------------------------
      DATA ETA_SL  /  0.01, 0.01, 0.01, 0.01 /  ! Friend et al. (1995)
      DATA KPAR    /  0.50, 0.50, 0.50, 0.50 /  ! Friend et al. (1995)
      DATA NR_NL   /  1.00, 1.00, 1.00, 1.00 /  !
      DATA NS_NL   /  0.04, 0.10, 1.00, 1.00 /  !
      DATA R_GROW  /  0.25, 0.25, 0.25, 0.25 /  ! Bonan (1995)
      DATA SIGL    /  0.04, 0.10, 0.04, 0.04 /  ! Schulze et al. (1994)


      IF (LTIMER) THEN
        CALL TIMER('SFSTOM  ',103)
      ENDIF

!-----------------------------------------------------------------------
! Create index of vegetated points on the land grid
!-----------------------------------------------------------------------

      VEG_PTS=0
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        IF (VEGF(L).GT.0.0 .AND. LAI(L).GT.0.0) THEN
          VEG_PTS = VEG_PTS + 1
          VEG_INDEX(VEG_PTS) = L
        ELSE
          GC(L) = 0.0
          NPP(L) = 0.0
          GPP(L) = 0.0
          RESP_P(L) = 0.0
        ENDIF
        ENDDO ! Loop over land-points

!-----------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit and the surface
! density of the air
!-----------------------------------------------------------------------

      CALL QSAT(QS(P1),TSTAR(P1),PSTAR(P1),P_POINTS)

      DO I=P1,P1+P_POINTS-1
        DQ(I) = MAX(0.0,(QS(I) - Q1(I)))
        RHOSTAR(I) = PSTAR(I) / (RAIR * TSTAR(I))
      ENDDO

!-----------------------------------------------------------------------
! Calculate the soil water factor (Cox and Huntingford, 1995)
!-----------------------------------------------------------------------

      DO L=LAND1,LAND1+LAND_PTS-1

        IF (V_ROOT(L) .GT. V_CRIT(L)) THEN
          FSMC(L) = 1.0
        ELSEIF (V_ROOT(L) .LE. V_WILT(L)) THEN
          FSMC(L) = 0.0
        ELSE
          FSMC(L) = (V_ROOT(L) - V_WILT(L))
     &            / (V_CRIT(L) - V_WILT(L))
        ENDIF

      ENDDO

!-----------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        FPAR(L) = (1 - EXP(-KPAR(FT(L))*LAI(L))) / KPAR(FT(L))

      ENDDO


!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Intially set the canopy humidity
! deficit using the previous value of GC.
!-----------------------------------------------------------------------
      DO K=1,ITER

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          I = LAND_INDEX(L)

          RA_RC(L) = RA(I) * GC(L)
          DQC(L) = DQ(I) / (1 + RA_RC(L))
          CO2C(L) = CO2

        ENDDO

!-----------------------------------------------------------------------
! Call CANOPY to calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------

        CALL CANOPY (LAND_PTS,LAND_FIELD,P1
     &,              LAND_INDEX
     &,              P_POINTS,P_FIELD
     &,              VEG_PTS,VEG_INDEX
     &,              FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR
     &,              NL0,FPAR,FSMC,LAI
     &,              GC,ANETC,CI,RDC,LTIMER)

      ENDDO ! LOOP OVER ITER

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!-----------------------------------------------------------------------
! Calculate the mean leaf nitrogen concentration assuming perfect
! light acclimation
!-----------------------------------------------------------------------
        NL(L) = (FPAR(L) / LAI(L)) * NL0(L)

!-----------------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------

        N_LEAF(L) = NL(L) * SIGL(FT(L)) * LAI(L)
        N_ROOT(L) = NR_NL(FT(L)) * NL(L) * ROOT(L)
        N_STEM(L) = NS_NL(FT(L)) * NL(L)
     &            * ETA_SL(FT(L)) * HT(L) * LAI(L)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity and the plant maintenance
! respiration rate in kg C/m2/sec
!-----------------------------------------------------------------------

        GPP(L) = 12.0E-3 * (ANETC(L) + RDC(L)*FSMC(L))
        RESP_P_M(L) = 12.0E-3 * RDC(L)
     &     * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------

        RESP_P_G(L) = R_GROW(FT(L)) * (GPP(L) - RESP_P_M(L))
        RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
        NPP(L) = GPP(L) - RESP_P(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SFSTOM  ',104)
      ENDIF

      RETURN
      END

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)
!
! Written by Peter Cox (May 1995)
!***********************************************************************

      SUBROUTINE CANOPY (LAND_PTS,LAND_FIELD,P1
     &,                  LAND_INDEX
     &,                  P_POINTS
     &,                  P_FIELD,VEG_PTS,VEG_INDEX
     &,                  FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR,NL0
     &,                  FPAR,FSMC,LAI
     &,                  GC,ANETC,CI,RDC,LTIMER)

      IMPLICIT NONE

      INTEGER
     & LAND_PTS               ! IN Number of land points to be processed
     &,LAND_FIELD             ! IN Total number of land points
     &,P_FIELD                ! IN Total number of p points
     &,LAND_INDEX(LAND_FIELD) ! IN Index of land points.
     &,P1                     ! IN First P point to be processed.
     &,P_POINTS               ! IN Number of P points to be processed.
     &,VEG_PTS                ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)  ! IN Index of vegetated points
!                                  on the land grid.

      INTEGER
     & FT(LAND_FIELD)         ! IN Plant functional type.

      REAL
     & CO2C(LAND_FIELD)       ! IN Canopy level CO2 concentration
!                                  (kg CO2/kg air).
     &,DQC(LAND_FIELD)        ! IN Canopy level specific humidity
!                                  deficit (kg H2O/kg air).
     &,O2                     ! IN Atmospheric O2 concentration
!                                  (kg O2/kg air).
     &,PSTAR(P_FIELD)         ! IN Surface pressure (Pa).
     &,IPAR(P_FIELD)          ! IN Incident PAR (W/m2).
     &,TSTAR(P_FIELD)         ! IN Surface temperature (K).
     &,NL0(LAND_FIELD)        ! IN Nitrogen concentration of
!                                  top leaf (kg N/kg C).
     &,FPAR(LAND_FIELD)       ! IN PAR absorption factor.
     &,FSMC(LAND_FIELD)       ! IN Soil water factor.
     &,LAI(LAND_FIELD)        ! IN Leaf area index (m2 leaf/m2 ground).


      REAL
     & ANETC(LAND_FIELD)      ! OUT Net canopy photosynthesis
!                                  (mol CO2/m2/s).
     &,CI(LAND_FIELD)         ! OUT Internal CO2 concentration
!                                   (mol CO2/m3).
     &,GC(LAND_FIELD)         ! OUT Canopy conductance for H2O (m/s).
     &,RDC(LAND_FIELD)        ! OUT Canopy dark respiration
!                                   (mol CO2/m2/s).

! WORK
      REAL
     & ANETL(LAND_FIELD)      ! WORK Net leaf photosynthesis
!                                    (mol CO2/m2/s/LAI).
     &,APAR(LAND_FIELD)       ! WORK PAR absorbed by the top leaf (W/m2)
     &,CA(LAND_FIELD)         ! WORK Canopy level CO2 pressure (Pa).
     &,DQM(LAND_FIELD)        ! WORK Canopy level humidity
!                                    deficit (mol H2O/m3).
     &,GL(LAND_FIELD)         ! WORK Leaf conductance for H2O (m/s).
     &,OA(LAND_FIELD)         ! WORK Atmospheric O2 pressure (Pa).
     &,RD(LAND_FIELD)         ! WORK Dark respiration of top leaf
!                                    (mol CO2/m2/s).

      LOGICAL
     & LTIMER

      INTEGER
     & I,J,L                  ! WORK Loop counters.

!-----------------------------------------------------------------------
! Functional Type dependent parameters
!-----------------------------------------------------------------------
      REAL
     & OMEGA(4)               ! Leaf scattering coefficient for PAR.
!-----------------------------------------------------------------------
!                       BT    NT   C3G   C4G
!-----------------------------------------------------------------------
      DATA OMEGA   /  0.15, 0.15, 0.15, 0.17 /  ! Sellers et al. (1994)

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL
     & R                      ! Gas constant (J/K/mol)
      PARAMETER (R = 8.3144)

      REAL
     & EPSILON                ! Ratio of molecular weights of water
!                               and dry air.
     &,EPCO2                  ! Ratio of molecular weights of CO2
!                               and dry air.
     &,EPO2                   ! Ratio of molecular weights of O2
!                               and dry air.
      PARAMETER (EPSILON = 0.62198, EPCO2 = 1.5194, EPO2 = 1.106)

      IF (LTIMER) THEN
        CALL TIMER('CANOPY  ',103)
      ENDIF

!-----------------------------------------------------------------------
! Calculate the atmospheric pressures of CO2 and O2
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)

        CA(L) = CO2C(L) / EPCO2 * PSTAR(I)
        OA(L) = O2 / EPO2 * PSTAR(I)
        DQM(L) = DQC(L) / EPSILON * PSTAR(I) / (R * TSTAR(I))

!-----------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf
!-----------------------------------------------------------------------
        APAR(L) = (1 - OMEGA(FT(L))) * IPAR(I)

      ENDDO

!-----------------------------------------------------------------------
! Call the leaf level model for the top leaf of the C3 and C4 plants
!-----------------------------------------------------------------------

      CALL LEAF_C3 (LAND_PTS
     &,             LAND_INDEX,P1
     &,             LAND_FIELD,P_FIELD
     &,             VEG_PTS,VEG_INDEX
     &,             FT,DQC,APAR,TSTAR,CA,OA,PSTAR
     &,             NL0,FSMC
     &,             GL,ANETL,CI,RD,LTIMER)


      CALL LEAF_C4 (LAND_PTS
     &,             LAND_INDEX,P1
     &,             LAND_FIELD,P_FIELD
     &,             VEG_PTS,VEG_INDEX
     &,             FT,DQC,APAR,TSTAR,CA,OA,PSTAR
     &,             NL0,FSMC
     &,             GL,ANETL,CI,RD,LTIMER)

!-----------------------------------------------------------------------
! Scale-up to the canopy level
!-----------------------------------------------------------------------

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
