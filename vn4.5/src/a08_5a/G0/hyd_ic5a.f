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
!   Subroutine HYD_INTCTL------------------------------------------
!
!   Level 3 control routine
!
!   Purpose: Calls HYDROL to calculate and add hydrology increments.
!            Avoids the need for *IF DEF around calls to different
!            versions of HYDROL.
!            Multilayer Version.
!
!   Written for the CRAY YMP
!
!      Author C.Bunton      Reviewer J.Lean 4/10/94
!
!   Modification history:
!
!   version  Date
!   4.1     5/6/96  New deck.
!   4.4   29/10/97  Modified for prognostic snow albedo scheme
!   
!
!   Programming standard : unified model documentation paper No 3
!
!   System components covered : P25
!
!   System task : P0
!
!   Documentation: Unified Model documentation paper P0
!                  version No 11 dated (26/11/90)
!  END -----------------------------------------------------------------
!   Arguments

      SUBROUTINE HYD_INTCTL(
     &     LAND,LICE_PTS,LICE_INDEX,ST_LEVELS,SM_LEVELS,
     &     SOIL_PTS,SOIL_INDEX,
     &     B_EXP,CAN_CPY,CANOPY_EVAPORATION,
     &     CONV_RAIN,CONV_SNOW,EXT,
     &     HCAP,HCON,INFIL_FAC,
     &     LAYER_DEPTH,LS_RAIN,LS_SNOW,
     &     ROOTD,SATCON,SATHH,SNOW_SUBLIMATION,
     &     SOILB,SOIL_EVAPORATION,SURF_HT_FLUX,
     &     VFRAC,VSAT,VWILT,TIMESTEP,
     &     CAN_WCNT,RGRAIN,L_SNOW_ALBEDO,SNODEP,STHF,STHU,
     &     TSTAR,T_DEEP_SOIL,
     &     INFIL,STF_HF_SNOW_MELT,
     &     HF_SNOW_MELT,SMC,SMCL,
     &     SNOW_MELT,STF_HF_SNOMLT_SUB,
     &     SNOMLT_SUB_HTF,
     &     STF_SUB_SURF_ROFF,SUB_SURF_ROFF,SURF_ROFF,
     &     TOT_TFALL
! Additional arguments for 7A hydrology (MOSES II) 
     &    ,TILE_PTS,TILE_INDEX
     &    ,CAN_CPY_TILE,CAN_EVAP_TILE
     &    ,FRAC,SNOW_FRAC,SOIL_SURF_HTF,SNOW_SURF_HTF,TSTAR_SNOW
     &    ,CAN_WCNT_TILE,TSNOW
C LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE

!----------------------------------------------------------------------
! Some of the following variables are 'dummy' as they are used by
! other HYDROL versions
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! IN Variables
!----------------------------------------------------------------------
      INTEGER
     &   LAND,                    ! IN No. land only points
     &   LICE_PTS,                ! IN No. of land-ice points
     &   LICE_INDEX(LAND),        ! IN Index of land-ice pts.
!                                 !    on land grid
     &   SM_LEVELS,               ! IN No. of soil moisture levels
     &   ST_LEVELS,               ! IN No. deep soil temp. levels
     &   SOIL_PTS,                ! IN No. of soil points
!                                 !    excludes land-ice points
     &   SOIL_INDEX(LAND)         ! IN Index of soil points on land
!                                 !    grid excludes land-ice pts.
      REAL
     &  B_EXP(LAND),              ! IN Exponent used in calculation
!                                 !    of soil water suction and
!                                 !    hydraulic conductivity
     &  CAN_CPY(LAND),            ! IN Canopy Capacity (Kg/m2)
     &  CANOPY_EVAPORATION(LAND), ! IN Canopy evaporation (Kg/m2/s)
     &  CONV_RAIN(LAND),          ! IN Convective rain (Kg/m2/s)
     &  CONV_SNOW(LAND),          ! IN Convective snow(Kg/m2/s)
     &  EXT(LAND,SM_LEVELS),      ! IN Extraction of water from each
!                                 !    soil layer (kg/m2/s)
     &  HCAP(LAND),               ! IN Soil heat capacity (J/K/m3)
     &  HCON(LAND),               ! IN Soil thermal conductivity (W/M/K)
     &  INFIL_FAC(LAND),          ! IN Infiltration enhancement factor
     &  LAYER_DEPTH(SM_LEVELS),   ! IN Layer depths
     &  LS_RAIN(LAND),            ! IN large scale rain (Kg/m2/s)
     &  LS_SNOW(LAND),            ! IN large scale snow (Kg/m2/s)
     &  ROOTD(LAND),              ! IN Root depth (m)
     &  SATCON(LAND),             ! IN Saturated hydraulic conductivity
!                                 !    (Kg/m2/s)
     &  SATHH(LAND),              ! IN Soil water suction (m)
     &  SNOW_SUBLIMATION(LAND),   ! IN Sublimation of snow (Kg/m2/s)
     &  SOILB,                    ! IN Dummy
     &  SOIL_EVAPORATION(LAND),   ! IN Surface evaporation (Kg/m2/s)
     &  SURF_HT_FLUX(LAND),       ! IN Net downward surface heat flux
!                                 !    (W/m2)
     &  TSTAR(LAND),              ! IN Surface temperature (K) 
     &  VFRAC(LAND),              ! IN Vegetated fraction
     &  VSAT(LAND),               ! IN Volumetric soil moisture content
!                                 !    at saturation (m3/m3 soil)
     &  VWILT(LAND),              ! IN Volumetric soil moisture at
!                                 !     wilting point (m3/m3 soil)
     &  TIMESTEP                  ! IN seconds between
!                                 !    physics routines call

      LOGICAL
     &   L_SNOW_ALBEDO,           ! IN Flag for prognostic snow albedo 
     &   STF_HF_SNOW_MELT,        ! IN Stash flag for snow melt heat
!                                 !    flux
     &   STF_HF_SNOMLT_SUB,       ! IN Stash flag for sub-surface
!                                 !    snow melt heat flux
     &   STF_SUB_SURF_ROFF        ! IN Stash flag for sub-surface runoff
!---------------------------------------------------------------------
! INOUT Variables
!----------------------------------------------------------------------
      REAL
     &  CAN_WCNT(LAND),           ! INOUT Canopy water content (Kg/m2)
     &  RGRAIN(LAND),             ! INOUT Snow grain size (microns)
     &  SNODEP(LAND),             ! INOUT Snow depth (Kg of water)
     &  STHF(LAND,SM_LEVELS),     ! INOUT Frozen soil moisture content
!                                 !       of each layer as a fraction
!                                 !       of saturation
     &  STHU(LAND,SM_LEVELS),     ! INOUT UNfrozen soil moisture cont
!                                 !       of each layer as a fraction
!                                 !       of saturation
     &  T_DEEP_SOIL(LAND,ST_LEVELS) ! INOUT Deep soil temps. (K)
!---------------------------------------------------------------------
! OUT Variables
!----------------------------------------------------------------------
      REAL
     &  INFIL(LAND),              ! OUT Maximum surface infiltration
!                                 !     rate (kg/m2/s)
     &  HF_SNOW_MELT(LAND),       ! OUT Total snow melt heat flux (W/m2)
     &  SMC(LAND),                ! OUT Available soil moisture
!                                       (Kg/m2/s)
     &  SMCL(LAND,SM_LEVELS),     ! OUT Soil moisture content of each
!                                 !     layer (Kg/Kg)

     &  SNOW_MELT(LAND),          ! OUT Snow melt (Kg/m2/s)
     &  SNOMLT_SUB_HTF(LAND),     ! OUT Sub-surface snow melt heat flux
!                                 !    (W/m2)
     &  SUB_SURF_ROFF(LAND),      ! OUT Subsurface runoff (Kg/m2/s)
     &  SURF_ROFF(LAND),          ! OUT Surface runoff (Kg/m2/s)
     &  TOT_TFALL(LAND)           ! OUT Total throughfall (Kg/m2/s)
!
! Additional arguments for 7A hydrology (MOSES II)
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
     &   TILE_PTS(NTYPE),
     &   TILE_INDEX(NTYPE)
      REAL
     &   CAN_CPY_TILE(NTYPE-1),
     &   CAN_EVAP_TILE(NTYPE-1),
     &   CAN_WCNT_TILE(NTYPE-1),
     &   FRAC(NTYPE),
     &   SNOW_FRAC,
     &   SOIL_SURF_HTF,
     &   SNOW_SURF_HTF,
     &   TSTAR_SNOW,
     &   TSNOW

!
      LOGICAL LTIMER             ! Logical switch for TIMER diags
!
!      CHARACTER*80
!     &       CMESSAGE     ! Error message if return code >0

!  External subroutines called

      EXTERNAL
     &      HYDROL


!--------------- SECTION 8 HYDROLOGY -----------------------------------

!  8.2 Call HYDROL to calculate and add hydrology increments

! As ST_LEVELS = SM_LEVELS in the MOSES scheme (HYDROL5A) only
! SM_LEVELS is input to HYDROL5A to avoid changing the variable NSHYD
! in HYDROL and below.
!
      CALL HYDROL(
     &     LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,
     &     LAND,SM_LEVELS,
     &     B_EXP,CAN_CPY,
     &     CONV_RAIN,CONV_SNOW,CANOPY_EVAPORATION,
     &     EXT,HCAP,HCON,INFIL_FAC,LS_RAIN,LS_SNOW,
     &     ROOTD,SATCON,SATHH,SNOW_SUBLIMATION,
     &     SURF_HT_FLUX,TIMESTEP,
     &     VFRAC,VSAT,VWILT,
     &     TSTAR,RGRAIN,L_SNOW_ALBEDO,
     &     CAN_WCNT,
     &     HF_SNOW_MELT,
     &     STF_HF_SNOW_MELT,SMCL,
     &     STHF,STHU,SNODEP,T_DEEP_SOIL,
     &     INFIL,SMC,SNOW_MELT,
     &     SNOMLT_SUB_HTF,
     &     STF_SUB_SURF_ROFF,
     &     SUB_SURF_ROFF,SURF_ROFF,
     &     TOT_TFALL
C LOGICAL LTIMER
     +,LTIMER
     +)


      RETURN
      END
