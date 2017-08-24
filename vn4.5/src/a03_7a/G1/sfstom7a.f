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
! Written by Peter Cox (June 1997)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
!**********************************************************************
      SUBROUTINE SF_STOM  (LAND_FIELD,LAND_INDEX,P_FIELD
     &,                    VEG_PTS,VEG_INDEX
     &,                    FT,CO2,CO2_3D,CO2_DIM,L_CO2_INTERACTIVE
     &,                    FSMC,HT,IPAR,LAI,PSTAR
     &,                    Q1,RA,TSTAR
     &,                    GPP,NPP,RESP_P,RESP_W,GC)


      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,LAND_INDEX(LAND_FIELD)     ! IN Index of land points on the
!                                 !    P-grid.
     &,P_FIELD                    ! IN Total number of P-gridpoints.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points
!                                 !    on the land grid.
     &,CO2_DIM           ! dimension of CO2 field

      INTEGER
     & FT                         ! IN Plant functional type.
      LOGICAL L_CO2_INTERACTIVE   ! switch for 3D CO2 field

      REAL
     & CO2                        ! IN Atmospheric CO2 concentration
     &,CO2_3D(CO2_DIM)            ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
     &,FSMC(LAND_FIELD)           ! IN Soil water factor.
     &,HT(LAND_FIELD)             ! IN Canopy height (m).
     &,IPAR(P_FIELD)              ! IN Incident PAR (W/m2).
     &,LAI(LAND_FIELD)            ! IN Leaf area index.
     &,PSTAR(LAND_FIELD)          ! IN Surface pressure (Pa).
     &,Q1(P_FIELD)                ! IN Specific humidity of level 1
!                                 !    (kg H2O/kg air).
     &,RA(LAND_FIELD)             ! IN Aerodynamic resistance (s/m).
     &,TSTAR(LAND_FIELD)          ! IN Surface temperature (K).
     &,GPP(LAND_FIELD)            ! OUT Gross Primary Productivity
!                                 !     (kg C/m2/s).
     &,NPP(LAND_FIELD)            ! OUT Net Primary Productivity
!                                 !     (kg C/m2/s).
     &,RESP_P(LAND_FIELD)         ! OUT Plant respiration rate
!                                 !     (kg C/m2/sec).
     &,RESP_W(LAND_FIELD)         ! OUT Wood respiration rate
!                                 !     (kg C/m2/sec).
     &,GC(LAND_FIELD)             ! INOUT Canopy resistance to H2O
!                                 !       (m/s).

      REAL
     & ANETC(LAND_FIELD)          ! WORK Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CO2C(LAND_FIELD)           ! WORK Canopy level CO2 concentration
!                                 !      (kg CO2/kg air).
     &,CI(LAND_FIELD)             ! WORK Internal CO2 pressure (Pa).
     &,DQ(LAND_FIELD)             ! WORK Specific humidity deficit
!                                 !      (kg H2O/kg air).
     &,DQC(LAND_FIELD)            ! WORK Canopy level specific humidity
!                                 !      deficit (kg H2O/kg air).
     &,FPAR(LAND_FIELD)           ! WORK PAR absorption factor.
     &,LAI_BAL(LAND_FIELD)        ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,NL(LAND_FIELD)             ! WORK Mean leaf nitrogen
!                                 !      concentration (kg N/kg C).
     &,NL_BAL(LAND_FIELD)         ! WORK Mean leaf nitrogen
!                                 !      concentration in balanced
!                                 !      growth state (kg N/kg C).
     &,N_LEAF(LAND_FIELD)         ! WORK Nitrogen contents of the leaf,
     &,N_ROOT(LAND_FIELD)         !      root,
     &,N_STEM(LAND_FIELD)         !      and stem (kg N/m2).
     &,QS(LAND_FIELD)             ! WORK Saturated specific humidity
!                                 !      (kg H2O/kg air).
     &,RA_RC(LAND_FIELD)          ! WORK Ratio of aerodynamic resistance
!                                 !      to canopy resistance.
     &,RDC(LAND_FIELD)            ! WORK Canopy dark respiration,
!                                 !      without soil water dependence
!                                 !      (mol CO2/m2/s).
     &,RESP_P_G(LAND_FIELD)       ! WORK Plant growth respiration rate
!                                 !      (kg C/m2/sec).
     &,RESP_P_M(LAND_FIELD)       ! WORK Plant maintenance respiration
!                                 !      rate (kg C/m2/sec).
     &,ROOT(LAND_FIELD)           ! WORK Root carbon (kg C/m2).

      INTEGER
     & I,J,K,L                    ! WORK Loop counters.

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL
     & RAIR                       ! Gas constant for dry air (J/kg/K).
      PARAMETER (RAIR = 287.05)

      REAL
     & O2                         ! Atmospheric concentration of
!                                 ! oxygen (kg O2/kg air).
      PARAMETER (O2 = 0.23)

      INTEGER
     & ITER                       ! Number of iterations to
!                                 ! determine the canopy climate.
      PARAMETER (ITER = 1)

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

!-----------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit and the surface
! density of the air
!-----------------------------------------------------------------------
      CALL QSAT(QS,TSTAR,PSTAR,LAND_FIELD)
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)
        DQ(L) = MAX(0.0,(QS(L) - Q1(I)))
      ENDDO

!-----------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        FPAR(L) = (1 - EXP(-KPAR(FT)*LAI(L))) / KPAR(FT)

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
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          RA_RC(L) = RA(L) * GC(L)
          DQC(L) = DQ(L) / (1 + RA_RC(L))
        ENDDO
      IF (L_CO2_INTERACTIVE) THEN
!  use full 3D CO2 field
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          I = LAND_INDEX(L)
          CO2C(L) = CO2_3D(I)
        ENDDO
      ELSE
!  just use single CO2_MMR value
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          CO2C(L) = CO2
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Call CANOPY to calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
        CALL CANOPY (LAND_FIELD,LAND_INDEX,P_FIELD
     &,              VEG_PTS,VEG_INDEX
     &,              FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR
     &,              FPAR,FSMC,LAI
     &,              GC,ANETC,CI,RDC)

      ENDDO

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!-----------------------------------------------------------------------
! Assume that root biomass is equal to balanced growth leaf biomass
!-----------------------------------------------------------------------
        LAI_BAL(L) = (A_WS(FT)*ETA_SL(FT)*HT(L)/A_WL(FT))
     &             **(1.0/(B_WL(FT)-1))
        ROOT(L) = SIGL(FT) * LAI_BAL(L)

!-----------------------------------------------------------------------
! Calculate the actual and balanced mean leaf nitrogen concentration
! assuming perfect light acclimation
!-----------------------------------------------------------------------
        NL(L) = (FPAR(L) / LAI(L)) * NL0(FT)
        NL_BAL(L) = (1 - EXP(-KPAR(FT)*LAI_BAL(L)))
     &            / (KPAR(FT)*LAI_BAL(L)) * NL0(FT)

!-----------------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
        N_LEAF(L) = NL(L) * SIGL(FT) * LAI(L)
        N_ROOT(L) = NR_NL(FT) * NL_BAL(L) * ROOT(L)
        N_STEM(L) = NS_NL(FT) * NL_BAL(L) * ETA_SL(FT) * HT(L) * LAI(L)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity, the plant maintenance
! respiration rate, and the wood maintenance respiration rate
! in kg C/m2/sec
!-----------------------------------------------------------------------
        GPP(L) = 12.0E-3 * (ANETC(L) + RDC(L)*FSMC(L))
        RESP_P_M(L) = 12.0E-3 * RDC(L)
     &     * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)
        RESP_W(L) = 12.0E-3 * RDC(L) * N_STEM(L) / N_LEAF(L)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
        RESP_P_G(L) = R_GROW(FT) * (GPP(L) - RESP_P_M(L))
        RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
        NPP(L) = GPP(L) - RESP_P(L)

      ENDDO

      RETURN
      END

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)
!
! Written by Peter Cox (May 1995)
!***********************************************************************
      SUBROUTINE CANOPY (LAND_FIELD,LAND_INDEX,P_FIELD
     &,                  VEG_PTS,VEG_INDEX
     &,                  FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR
     &,                  FPAR,FSMC,LAI
     &,                  GC,ANETC,CI,RDC)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,LAND_INDEX(LAND_FIELD)     ! IN Index of land points on the
!                                 !    P-grid.
     &,P_FIELD                    ! IN Total number of P-gridpoints.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points
!                                 !    on the land grid.

      INTEGER
     & FT                         ! IN Plant functional type.

      REAL
     & CO2C(LAND_FIELD)           ! IN Canopy level CO2 concentration
!                                 !    (kg CO2/kg air).
     &,DQC(LAND_FIELD)            ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,O2                         ! IN Atmospheric O2 concentration
!                                 !    (kg O2/kg air).
     &,PSTAR(LAND_FIELD)          ! IN Surface pressure (Pa).
     &,IPAR(P_FIELD)              ! IN Incident PAR (W/m2).
     &,TSTAR(LAND_FIELD)          ! IN Surface temperature (K).
     &,FPAR(LAND_FIELD)           ! IN PAR absorption factor.
     &,FSMC(LAND_FIELD)           ! IN Soil water factor.
     &,LAI(LAND_FIELD)            ! IN Leaf area index
!                                 !    (m2 leaf/m2 ground).


      REAL
     & ANETC(LAND_FIELD)          ! OUT Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CI(LAND_FIELD)             ! OUT Internal CO2 concentration
!                                 !     (mol CO2/m3).
     &,GC(LAND_FIELD)             ! OUT Canopy conductance for H2O
!                                 !     (m/s).
     &,RDC(LAND_FIELD)            ! OUT Canopy dark respiration
!                                 !     (mol CO2/m2/s).
     &,ANETL(LAND_FIELD)          ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
     &,APAR(LAND_FIELD)           ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
     &,CA(LAND_FIELD)             ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
     &,DQM(LAND_FIELD)            ! WORK Canopy level humidity
!                                 !      deficit (mol H2O/m3).
     &,GL(LAND_FIELD)             ! WORK Leaf conductance for H2O
!                                 !      (m/s).
     &,OA(LAND_FIELD)             ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
     &,RD(LAND_FIELD)             ! WORK Dark respiration of top leaf
!                                 !      (mol CO2/m2/s).

      INTEGER
     & I,J,L                      ! WORK Loop counters.

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

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL
     & R                          ! Gas constant (J/K/mol)
      PARAMETER (R = 8.3144)
C*L------------------COMDECK CCARBON------------------------------------
C carbon cycle and vegetation parameters
      REAL
     & M_CO2                      ! molecular weight of CO2
     &,M_AIR                      ! molecular weight of dry air
     &,EPSILON                    ! Ratio of molecular weights of water
!                                 !  and dry air.
     &,EPCO2                      ! Ratio of molecular weights of CO2
!                                 !  and dry air.
     &,EPO2                       ! Ratio of molecular weights of O2
!                                 !  and dry air.
     &,CO2CONV_A2O                ! conversion factor for atmos to
!                                 !  ocean passing of CO2 (mmr to ppmv)
     &,CO2CONV_O2A                ! conversion factor for ocean to
!                                 !  atmos passing of CO2 flux
!                                 !  (mol C/m2/yr to Kg CO2/m2/s)

      PARAMETER (M_AIR=28.966, EPCO2=1.5194, M_CO2=M_AIR*EPCO2,
     &           EPSILON = 0.62198, EPO2 = 1.106)

      PARAMETER (CO2CONV_A2O = M_AIR * 1E6 / M_CO2,
     &           CO2CONV_O2A = M_CO2 * 1e-3 / (360.0 * 24.0 * 3600.0))
C*----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate the atmospheric pressures of CO2 and O2
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)

        CA(L) = CO2C(L) / EPCO2 * PSTAR(L)
        OA(L) = O2 / EPO2 * PSTAR(L)
        DQM(L) = DQC(L) / EPSILON * PSTAR(L) / (R * TSTAR(L))

!-----------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf
!-----------------------------------------------------------------------
        APAR(L) = (1 - OMEGA(FT)) * IPAR(I)

      ENDDO

!-----------------------------------------------------------------------
! Call the leaf level model for the top leaf of the C3 and C4 plants
!-----------------------------------------------------------------------

      IF ( C3(FT) .EQ. 1 ) THEN

        CALL LEAF_C3 (LAND_FIELD,VEG_PTS,VEG_INDEX,FT
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC
     &,               GL,ANETL,CI,RD)

      ELSE

        CALL LEAF_C4 (LAND_FIELD,VEG_PTS,VEG_INDEX,FT
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC
     &,               GL,ANETL,CI,RD)

      ENDIF

!-----------------------------------------------------------------------
! Scale-up to the canopy level
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        ANETC(L) = ANETL(L) * FPAR(L)
        GC(L) = FPAR(L) * GL(L)
        RDC(L) = RD(L) * FPAR(L)

      ENDDO

      RETURN

      END
