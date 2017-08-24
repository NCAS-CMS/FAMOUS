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
!!! Subroutine VEGCARB ------------------------------------------------
!!!
!!! Purpose : Updates carbon contents of the vegetation.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!
!!!END ----------------------------------------------------------------
       SUBROUTINE VEGCARB (LAND_FIELD,TRIF_PTS,TRIF_INDEX   
     &,                    N,FORW,GAMMA,G_LEAF,NPP,RESP_W
     &,                    LEAF,ROOT,WOOD,DCVEG,PC_S)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,TRIF_PTS                   ! IN Number of points on which 
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_FIELD)     ! IN Indices of land points on 
!                                 !    which TRIFFID may operate
     &,N                          ! IN Plant functional type.
     &,L,T                        ! WORK Loop counters

      REAL
     & FORW                       ! IN Forward timestep weighting.
     &,GAMMA                      ! IN Inverse timestep (/360days).
     &,G_LEAF(LAND_FIELD)         ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
     &,NPP(LAND_FIELD)            ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
     &,RESP_W(LAND_FIELD)         ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
     &,LEAF(LAND_FIELD)           ! INOUT Leaf biomass (kg C/m2).
     &,ROOT(LAND_FIELD)           ! INOUT Root biomass (kg C/m2).
     &,WOOD(LAND_FIELD)           ! INOUT Woody biomass (kg C/m2).
     &,DCVEG(LAND_FIELD)          ! OUT Change in vegetation carbon
C                                 !     during the timestep
C                                 !     (kg C/m2/timestep).
     &,PC_S(LAND_FIELD)           ! OUT Net carbon flux available
!                                 !     for spreading (kg C/m2/360days).
     &,DFPAR_DLAI                 ! WORK Rate of change of FPAR
C                                 !      with leaf area index.
     &,DLAI                       ! WORK Increment to the leaf area
C                                 !      index.
     &,DLAMG_DLAI,DLIT_DLAI       ! WORK Required for calculation
C                                 !      of the equilibrium increments.
     &,DNPP_DLAI(LAND_FIELD)      ! WORK Rate of change of NPP
C                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
     &,DPC_DLAI(LAND_FIELD)       ! WORK Rate of change of PC
C                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
     &,DPCG_DLAI(LAND_FIELD)      ! WORK Rate of change of PC_G
C                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
     &,DRESPW_DLAI                ! WORK Rate of change of RESP_W
C                                 !      with leaf area index
     &,FPAR                       ! WORK PAR interception factor.
     &,LAI(LAND_FIELD)            ! WORK Leaf area index.
     &,LAMBDA_G                   ! WORK Fraction of NPP available
C                                 !      for spreading.
     &,LIT_C_L(LAND_FIELD)        ! WORK Local rate of Carbon Litter
!                                 !      production (kg C/m2/360days).
     &,PC(LAND_FIELD)             ! WORK Net carbon flux available
!                                 !      to vegetation (kg C/m2/360days)
     &,PC_G(LAND_FIELD)           ! WORK Net carbon flux available
!                                 !      for growth (kg C/m2/360days).

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


      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 

        LAI(L) = LEAF(L)/SIGL(N)
C----------------------------------------------------------------------
C Calculate the local production rate for carbon litter
C----------------------------------------------------------------------
        LIT_C_L(L) = G_LEAF(L)*LEAF(L)+G_ROOT(N)*ROOT(L)
     &               + G_WOOD(N)*WOOD(L)

C----------------------------------------------------------------------
C Diagnose the net local carbon flux into the vegetation
C----------------------------------------------------------------------
        PC(L) = NPP(L) - LIT_C_L(L)

C----------------------------------------------------------------------
C Variables required for the implicit and equilibrium calculations
C----------------------------------------------------------------------
        DLIT_DLAI = (G_LEAF(L)*LEAF(L)+G_ROOT(N)*ROOT(L))/LAI(L)
     &            + B_WL(N)*G_WOOD(N)*WOOD(L)/LAI(L)

        FPAR = (1 - EXP(-KPAR(N)*LAI(L)))/KPAR(N)
        DFPAR_DLAI = EXP(-KPAR(N)*LAI(L))

        DNPP_DLAI(L) = NPP(L)*DFPAR_DLAI/FPAR
     &               + (1-R_GROW(N))*RESP_W(L)
     &               *(DFPAR_DLAI/FPAR-B_WL(N)/LAI(L))

        LAMBDA_G = 1 - (LAI(L) - LAI_MIN(N))
     &                /(LAI_MAX(N) - LAI_MIN(N))
        DLAMG_DLAI = -1.0/(LAI_MAX(N) - LAI_MIN(N))

        PC_G(L) = LAMBDA_G * NPP(L) - LIT_C_L(L)
        DPCG_DLAI(L) = LAMBDA_G*DNPP_DLAI(L)
     &               + DLAMG_DLAI*NPP(L)
     &               - DLIT_DLAI
        DPC_DLAI(L) = DNPP_DLAI(L) - DLIT_DLAI

      ENDDO

C----------------------------------------------------------------------
C Update vegetation carbon contents
C----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 
        DCVEG(L) = LEAF(L)+ROOT(L)+WOOD(L)
      ENDDO

      CALL GROWTH (LAND_FIELD,TRIF_PTS,TRIF_INDEX   
     &,            N,DPCG_DLAI,FORW,GAMMA,PC_G,LEAF,ROOT,WOOD) 

      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 
        DCVEG(L) = LEAF(L)+ROOT(L)+WOOD(L)-DCVEG(L)
      ENDDO

C----------------------------------------------------------------------
C Diagnose the carbon available for spreading and apply implicit
C corrections to the driving fluxes.
C----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T) 
        DLAI = LEAF(L)/SIGL(N) - LAI(L)
        PC_S(L) = PC(L) + FORW*DPC_DLAI(L)*DLAI - DCVEG(L)*GAMMA

        FPAR = (1 - EXP(-KPAR(N)*LAI(L)))/KPAR(N)
        DFPAR_DLAI = EXP(-KPAR(N)*LAI(L))
        DRESPW_DLAI = RESP_W(L)*B_WL(N)/LAI(L)

        NPP(L) = NPP(L) + FORW*DNPP_DLAI(L)*DLAI
        RESP_W(L) = RESP_W(L) + FORW*DRESPW_DLAI*DLAI
      ENDDO

      RETURN
      END
