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
!!! Subroutine TRIFFID ------------------------------------------------
!!!
!!!                     Top-down
!!!                     Representation of
!!!                     Interactive
!!!                     Foliage and
!!!                     Flora
!!!                     Including
!!!                     Dynamics
!!!
!!! Purpose : Simulates changes in vegetation structure, areal
!!!           coverage and the carbon contents of vegetation and soil.
!!!           can be used to advance these variables dynamically
!!!           (GAMMA=1/TIMESTEP) or to iterate towards  equilibrium
!!!           (GAMMA --> 0.0, FORW=1.0).
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE TRIFFID (LAND_FIELD,TRIF_PTS,TRIF_INDEX,FORW,GAMMA    
     &,                   FRAC_VS,G_ANTH,G_LEAF,NPP,RESP_S,RESP_W
     &,                   CS,FRAC,HT,LAI,C_VEG,CV,LIT_C,LIT_C_T)


      IMPLICIT NONE

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

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,TRIF_PTS                   ! IN Number of points on which 
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_FIELD)     ! IN Indices of land points on 
!                                 !    which TRIFFID may operate
     &,L,N,T                      ! WORK Loop counters

      REAL
     & FORW                       ! IN Forward timestep weighting.
     &,FRAC_VS(LAND_FIELD)        ! IN Total fraction of gridbox 
!                                 !    covered by veg or soil.
     &,GAMMA                      ! IN Inverse timestep (/360days).
     &,G_ANTH(LAND_FIELD)         ! IN Anthropogenic disturbance rate
C                                 !    (/360days).
     &,G_LEAF(LAND_FIELD,NPFT)    ! IN Turnover rate for leaf and
C                                 !    fine root biomass (/360days).
     &,NPP(LAND_FIELD,NPFT)       ! INOUT Net primary productivity
C                                 !       (kg C/m2/360days).
     &,RESP_S(LAND_FIELD)         ! INOUT Soil respiration 
C                                 !       (kg C/m2/360days).
     &,RESP_W(LAND_FIELD,NPFT)    ! INOUT Wood maintenance respiration
C                                 !       (kg C/m2/360days).
     &,CS(LAND_FIELD)             ! INOUT Soil carbon (kg C/m2).
     &,FRAC(LAND_FIELD,NTYPE)     ! INOUT Fractional cover of each
C                                 !       Functional Type.
     &,HT(LAND_FIELD,NPFT)        ! INOUT Vegetation height (m).
     &,LAI(LAND_FIELD,NPFT)       ! INOUT Leaf area index.
     &,C_VEG(LAND_FIELD,NPFT)     ! OUT Total carbon content of
C                                 !     the vegetation (kg C/m2).
     &,CV(LAND_FIELD)             ! OUT Gridbox mean vegetation
C                                 !     carbon (kg C/m2).
     &,LIT_C(LAND_FIELD,NPFT)     ! OUT Carbon Litter (kg C/m2/360days).
     &,LIT_C_T(LAND_FIELD)        ! OUT Gridbox mean carbon litter
C                                 !     (kg C/m2/360days).

      REAL
     & DCVEG(LAND_FIELD,NPFT)     ! WORK Change in vegetation carbon
C                                 !      during the timestep
C                                 !      (kg C/m2/timestep).
     &,DFRAC(LAND_FIELD,NPFT)     ! WORK Change in areal fraction
C                                 !      during the timestep
C                                 !      (/timestep).
     &,LAI_BAL(LAND_FIELD,NPFT)   ! WORK Leaf area index in balanced
C                                 !      growth state.
     &,LEAF(LAND_FIELD,NPFT)      ! WORK Leaf biomass (kg C/m2).
     &,PC_S(LAND_FIELD,NPFT)      ! WORK Net carbon flux available
C                                 !      for spreading 
C                                 !      (kg C/m2/yr).
     &,PHEN(LAND_FIELD,NPFT)      ! WORK Phenological state.
     &,ROOT(LAND_FIELD,NPFT)      ! WORK Root biomass (kg C/m2).
     &,WOOD(LAND_FIELD,NPFT)      ! WORK Woody biomass (kg C/m2).

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

C----------------------------------------------------------------------
C Loop through Functional Types
C----------------------------------------------------------------------
      DO N=1,NPFT

C----------------------------------------------------------------------
C Loop through TRIFFID points   
C----------------------------------------------------------------------
        DO T=1,TRIF_PTS
          L=TRIF_INDEX(T) 

C----------------------------------------------------------------------
C Diagnose the balanced-growth leaf area index and the associated leaf,
C wood, root and total vegetation carbon
C----------------------------------------------------------------------
          LAI_BAL(L,N) = (A_WS(N)*ETA_SL(N)*HT(L,N)
     &              /A_WL(N))**(1.0/(B_WL(N)-1))
          LEAF(L,N) = SIGL(N)*LAI_BAL(L,N)
          ROOT(L,N) = LEAF(L,N)
          WOOD(L,N) = A_WL(N)*(LAI_BAL(L,N)**B_WL(N))
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)
C----------------------------------------------------------------------
C Diagnose the phenological state
C----------------------------------------------------------------------
          PHEN(L,N) = LAI(L,N)/LAI_BAL(L,N)

        ENDDO

C----------------------------------------------------------------------
C Update vegetation carbon contents
C----------------------------------------------------------------------
        CALL VEGCARB (LAND_FIELD,TRIF_PTS,TRIF_INDEX,N,FORW   
     &,               GAMMA,G_LEAF(1,N),NPP(1,N),RESP_W(1,N)
     &,               LEAF(1,N),ROOT(1,N),WOOD(1,N)
     &,               DCVEG(1,N),PC_S(1,N))

      ENDDO

C-----------------------------------------------------------------------
C Diagnose the new value of Canopy Height, Leaf Area Index and Total
C Vegetation Carbon
C-----------------------------------------------------------------------
      DO N=1,NPFT

        DO T=1,TRIF_PTS
          L=TRIF_INDEX(T) 

          HT(L,N) = WOOD(L,N) / (A_WS(N) * ETA_SL(N))
     &            * (A_WL(N)/WOOD(L,N))**(1.0/B_WL(N))
          LAI_BAL(L,N) = LEAF(L,N) / SIGL(N)
          LAI(L,N) = PHEN(L,N) * LAI_BAL(L,N)
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)

        ENDDO

      ENDDO

C----------------------------------------------------------------------
C Update the areal coverage of each functional type
C----------------------------------------------------------------------
      CALL LOTKA (LAND_FIELD,TRIF_PTS,TRIF_INDEX
     &,           C_VEG,FORW,FRAC_VS,GAMMA,G_ANTH,LAI_BAL,PC_S
     &,           FRAC,DFRAC)

C----------------------------------------------------------------------
C Diagnose the litterfall from the carbon balance of each vegetation
C type (assumes explicit update).
C----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 

        LIT_C_T(L) = 0.0

        DO N=1,NPFT
          LIT_C(L,N) = NPP(L,N)-GAMMA*(C_VEG(L,N)*FRAC(L,N)
     &     -(C_VEG(L,N)-DCVEG(L,N))*(FRAC(L,N)-DFRAC(L,N)))/FRAC(L,N)
          LIT_C_T(L) = LIT_C_T(L)+FRAC(L,N)*LIT_C(L,N)
        ENDDO
      ENDDO

C----------------------------------------------------------------------
C Call SOIL_C to update the soil carbon content
C----------------------------------------------------------------------
      CALL SOILCARB (LAND_FIELD,TRIF_PTS,TRIF_INDEX   
     &,              FORW,GAMMA,LIT_C_T,RESP_S,CS)

C----------------------------------------------------------------------
C Diagnose the gridbox mean vegetation carbon
C----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 
        CV(L) = 0.0
        DO N=1,NPFT
          CV(L) = CV(L)+FRAC(L,N)*C_VEG(L,N)
        ENDDO
      ENDDO

      RETURN
      END
