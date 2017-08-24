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
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.
!
! Written by Peter Cox (June 1997)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
! 4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!**********************************************************************
      SUBROUTINE PHYSIOL (LAND_FIELD,LAND_PTS,LAND1
     &,                   LAND_INDEX
     &,                   P_FIELD,NSHYD,TILE_PTS,TILE_INDEX
     &,                   CO2,CO2_3D,CO2_DIM,L_CO2_INTERACTIVE
     &,                   CS,FRAC,HT,IPAR,LAI,PSTAR,Q1
     &,                   STHU,TIMESTEP,TSOIL,TSTAR_TILE
     &,                   V_CRIT,V_SAT,V_WILT,WIND,Z0V,Z1
     &,                   G_LEAF,GS,GS_TILE,GPP,GPP_FT,NPP,NPP_FT
     &,                   RESP_P,RESP_P_FT,RESP_S,RESP_W_FT,SMCT,WT_EXT)

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
     &,LAND_PTS                   ! IN Number of land points to be
!                                 !    processed.
     &,LAND1                      ! IN First land point to be
!                                 !    processed.
     &,LAND_INDEX(LAND_FIELD)     ! IN Index of land points on the
!                                 !    P-grid.
     &,P_FIELD                    ! IN Total number of P-gridpoints.
     &,CO2_DIM           ! dimension of CO2 field
     &,NSHYD                      ! IN Number of soil moisture
!                                 !    levels.
     &,TILE_PTS(NTYPE)            ! IN Number of land points which
!                                 !    include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE)
!                                 ! IN Indices of land points which
!                                 !    include the nth surface type.
      LOGICAL L_CO2_INTERACTIVE   ! switch for 3D CO2 field

      REAL
     & CO2                        ! IN Atmospheric CO2 concentration
     &,CO2_3D(CO2_DIM)            ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
     &,CS(LAND_FIELD)             ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_FIELD,NTYPE)     ! IN Tile fractions.
     &,HT(LAND_FIELD,NPFT)        ! IN Canopy height (m).
     &,IPAR(P_FIELD)              ! IN Incident PAR (W/m2).
     &,LAI(LAND_FIELD,NPFT)       ! IN Leaf area index.
     &,PSTAR(P_FIELD)             ! IN Surface pressure (Pa).
     &,Q1(P_FIELD)                ! IN Specific humidity at level 1
!                                 !    (kg H2O/kg air).
     &,STHU(LAND_FIELD,NSHYD)     ! IN Soil moisture content in each
!                                 !    layer as a fraction of saturation
     &,TIMESTEP                   ! IN Timestep (s).
     &,TSOIL(LAND_FIELD)          ! IN Soil temperature (K).
     &,TSTAR_TILE(LAND_FIELD,NTYPE)
!                                 ! IN Tile surface temperatures (K).
     &,V_CRIT(LAND_FIELD)         ! IN Volumetric soil moisture
!                                 !    concentration above which
!                                 !    stomata are not sensitive
!                                 !    to soil water (m3 H2O/m3 soil).
     &,V_SAT(LAND_FIELD)          ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
     &,V_WILT(LAND_FIELD)         ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).
     &,WIND(P_FIELD)              ! IN Windspeed (m/s).
     &,Z0V(LAND_FIELD,NTYPE)      ! IN Tile roughness lengths (m).
     &,Z1(P_FIELD)                ! IN Windspeed reference height(m).
     &,GS(LAND_FIELD)             ! INOUT Gridbox mean surface
!                                 !       conductance (m/s).

      REAL
     & G_LEAF(LAND_FIELD,NPFT)    ! OUT Leaf turnover rate (/360days).
     &,GS_TILE(LAND_FIELD,NTYPE)  ! OUT Surface conductance for
!                                 !     land tiles (m/s).
     &,GPP(LAND_FIELD)            ! OUT Gridbox mean gross primary
!                                 !     productivity (kg C/m2/s).
     &,GPP_FT(LAND_FIELD,NPFT)    ! OUT Gross primary productivity
!                                 !     (kg C/m2/s).
     &,NPP(LAND_FIELD)            ! OUT Gridbox mean net primary
!                                 !     productivity (kg C/m2/s).
     &,NPP_FT(LAND_FIELD,NPFT)    ! OUT Net primary productivity
!                                 !     (kg C/m2/s).
     &,RESP_P(LAND_FIELD)         ! OUT Gridbox mean plant respiration
!                                 !     (kg C/m2/s).
     &,RESP_P_FT(LAND_FIELD,NPFT) ! OUT Plant respiration (kg C/m2/s).
     &,RESP_S(LAND_FIELD)         ! OUT Soil respiration (kg C/m2/s).
     &,RESP_W_FT(LAND_FIELD,NPFT) ! OUT Wood maintenance respiration
!                                 !     (kg C/m2/s).
     &,SMCT(LAND_FIELD)           ! OUT Available moisture in the
!                                 !     soil profile (mm).
     &,WT_EXT(LAND_FIELD,NSHYD)   ! OUT Fraction of evapotranspiration
!                                 !     which is extracted from each
!                                 !     soil layer.

      REAL
     & F_ROOT(NSHYD)              ! WORK Fraction of roots in each soil
!                                 !      layer.
     &,FSMC(LAND_FIELD)           ! WORK Moisture availability factor.
     &,PSTAR_LAND(LAND_FIELD)     ! WORK Surface pressure (Pa).
     &,RA(LAND_FIELD)             ! WORK Aerodynamic resistance (s/m).
     &,RIB(P_FIELD)               ! WORK Bulk Richardson Number.

      INTEGER
     & I,J,K,L,N                  ! Loop indices

C-----------------------------------------------------------------------
C Surface parameters for each Plant Functional Type
C-----------------------------------------------------------------------
      REAL
     + ALBSNC_MAX(NPFT)           ! Snow-covered albedo for large LAI.
     +,ALBSNC_MIN(NPFT)           ! Snow-covered albedo for zero LAI.   
     +,ALBSNF_MAX(NPFT)           ! Snow-free albedo for large LAI.
     +,DZ0V_DH(NPFT)              ! Rate of change of vegetation
C                                 ! roughness length with height.
     +,CATCH0(NPFT)               ! Minimum canopy capacity (kg/m2).
     +,DCATCH_DLAI(NPFT)          ! Rate of change of canopy capacity
C                                 ! with LAI.
     +,INFIL_F(NPFT)              ! Infiltration enhancement factor.
     +,KEXT(NPFT)                 ! Light extinction coefficient.
     +,ROOTD_FT(NPFT)             ! Rootdepth (m).
C----------------------------------------------------------------------
C                           BT    NT   C3G   C4G    S
C----------------------------------------------------------------------
      DATA ALBSNC_MAX  /  0.15, 0.15, 0.60, 0.60, 0.40 /                
      DATA ALBSNC_MIN  /  0.30, 0.30, 0.80, 0.80, 0.80 /                
      DATA ALBSNF_MAX  /  0.10, 0.10, 0.20, 0.20, 0.20 /
      DATA DZ0V_DH     /  0.05, 0.05, 0.10, 0.10, 0.10 /
      DATA CATCH0      /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA DCATCH_DLAI /  0.05, 0.05, 0.05, 0.05, 0.05 /
      DATA INFIL_F     /  4.00, 4.00, 2.00, 2.00, 2.00 /
      DATA KEXT        /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA ROOTD_FT    /  3.00, 1.00, 0.50, 0.50, 0.50 /                
      REAL
     + ALBSNC_NVG(NNVG)           ! Snow-covered albedo.
     +,ALBSNF_NVG(NNVG)           ! Snow-free albedo.
     +,CATCH_NVG(NNVG)            ! Canopy capacity (kg/m2).
     +,GS_NVG(NNVG)               ! Surface conductance (m/s).
     +,INFIL_NVG(NNVG)            ! Infiltration enhancement factor.
     +,ROOTD_NVG(NNVG)            ! Rootdepth (m).
     +,Z0_NVG(NNVG)               ! Roughness length (m).
C----------------------------------------------------------------------
C                         Urban  Water  Soil   Ice
C----------------------------------------------------------------------
      DATA ALBSNC_NVG  /  0.40,  0.80,  0.80,  0.80 /
      DATA ALBSNF_NVG  /  0.18,  0.06, -1.00,  0.75 /
      DATA CATCH_NVG   /  0.50,  0.00,  0.50,  0.00 /
      DATA GS_NVG      /  5E-3,   1E2,  1E-2,  0.00 /
      DATA INFIL_NVG   /  0.10,  1.00,  0.50,  0.00 /
      DATA ROOTD_NVG   /  0.50,  1.00,  0.10,  0.00 /
      DATA Z0_NVG      /  1.50,  3E-4,  3E-4,  1E-4 /
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
      REAL
     + DZSOIL(4)               ! Soil layer thicknesses (m).
      DATA DZSOIL /0.100, 0.250, 0.650, 2.000 /
C-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      DO K=1,NSHYD
        F_ROOT(K)=0.0
        DO L=1,LAND_FIELD
          WT_EXT(L,K)=0.0
        ENDDO
      ENDDO

      DO L=1,LAND_FIELD
        I = LAND_INDEX(L)
        PSTAR_LAND(L) = PSTAR(I)
      ENDDO

      DO I=1,P_FIELD
        RIB(I)=0.0
      ENDDO

      DO N=1,NPFT
        DO L=1,LAND_FIELD
          G_LEAF(L,N)=0.0
          GPP_FT(L,N)=0.0
          NPP_FT(L,N)=0.0
          RESP_P_FT(L,N)=0.0
          RESP_W_FT(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NTYPE
        DO L=1,LAND_FIELD
          GS_TILE(L,N)=GS(L)
        ENDDO
      ENDDO

      DO L=1,LAND_FIELD
        GPP(L)=0.0
        NPP(L)=0.0
        RESP_P(L)=0.0
        RESP_S(L)=0.0
        SMCT(L)=0.0
        GS(L)=0.0
        FSMC(L)=0.0
        RA(L)=0.0
      ENDDO

!-----------------------------------------------------------------------
! Loop over Plant Functional Types to calculate the available moisture
! and the values of canopy conductance, the carbon fluxes and the leaf
! turnover rate
!-----------------------------------------------------------------------
      DO N=1,NPFT

        CALL ROOT_FRAC(NSHYD,DZSOIL,ROOTD_FT(N),F_ROOT)

        CALL SMC_EXT (LAND_FIELD,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)
     &,               F_ROOT,FRAC(1,N),STHU,V_CRIT,V_SAT,V_WILT
     &,               WT_EXT,FSMC)

        CALL RAERO (LAND_FIELD,LAND_INDEX,P_FIELD
     &,             TILE_PTS(N),TILE_INDEX(1,N)
     &,             RIB,WIND,Z0V(1,N),Z0V(1,N),Z1,RA)

        CALL SF_STOM (LAND_FIELD,LAND_INDEX,P_FIELD
     &,               TILE_PTS(N),TILE_INDEX(1,N),N
     &,               CO2,CO2_3D,CO2_DIM,L_CO2_INTERACTIVE
     &,               FSMC,HT(1,N),IPAR,LAI(1,N),PSTAR_LAND
     &,               Q1,RA,TSTAR_TILE(1,N)
     &,               GPP_FT(1,N),NPP_FT(1,N),RESP_P_FT(1,N)
     &,               RESP_W_FT(1,N),GS_TILE(1,N))

        CALL LEAF_LIT (LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N)
     &,                N,FSMC,TSTAR_TILE(1,N),G_LEAF(1,N))

      ENDDO

!----------------------------------------------------------------------
! Loop over non-vegetated surface types to calculate the available
! moisture and the surface conductance. Land-ice (tile NTYPE) excluded.
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE-1

        CALL ROOT_FRAC(NSHYD,DZSOIL,ROOTD_NVG(N-NPFT),F_ROOT)

        CALL SMC_EXT (LAND_FIELD,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)
     &,               F_ROOT,FRAC(1,N),STHU,V_CRIT,V_SAT,V_WILT
     &,               WT_EXT,FSMC)

        DO J=1,TILE_PTS(N)
          L=TILE_INDEX(J,N)
          GS_TILE(L,N)=FSMC(L)*GS_NVG(N-NPFT)
        ENDDO

      ENDDO

!----------------------------------------------------------------------
! Calculate the rate of soil respiration
!----------------------------------------------------------------------
      CALL MICROBE (LAND_FIELD,LAND_PTS,LAND1
     &,             CS,STHU,V_SAT,V_WILT,TSOIL,RESP_S)

!----------------------------------------------------------------------
! Form gridbox mean values
!----------------------------------------------------------------------
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L=TILE_INDEX(J,N)
          GS(L)=GS(L)+FRAC(L,N)*GS_TILE(L,N)
        ENDDO
      ENDDO

      DO N=1,NPFT
        DO J=1,TILE_PTS(N)
          L=TILE_INDEX(J,N)

          GPP(L)=GPP(L)+FRAC(L,N)*GPP_FT(L,N)
          NPP(L)=NPP(L)+FRAC(L,N)*NPP_FT(L,N)
          RESP_P(L)=RESP_P(L)+FRAC(L,N)*RESP_P_FT(L,N)

        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Diagnose the available moisture in the soil profile
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO L=LAND1,LAND1+LAND_PTS-1
          SMCT(L) = SMCT(L) + MAX( 0. ,
     &             RHO_WATER*DZSOIL(N)*(STHU(L,N)*V_SAT(L)-V_WILT(L)))
        ENDDO
      ENDDO

      RETURN
      END
