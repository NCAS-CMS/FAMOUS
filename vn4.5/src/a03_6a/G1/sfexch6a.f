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
!!!   SUBROUTINE SF_EXCH------------------------------------------------
!!!
!!!  Purpose: Calculate coefficients of turbulent exchange between
!!!           the surface and the lowest atmospheric layer, and
!!!           "explicit" fluxes between the surface and this layer.
!!!
!!!  Suitable for Single Column use.
!!!
!!!          Canopy evaporation made implicit
!!!     with respect to canopy water content (requiring TIMESTEP to be
!!!     passed in).
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4   10/09/95   New deck    R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE SF_EXCH (
     & P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,N_TYPES
     &,P1,LAND1
     &,LAND_INDEX,GATHER
     &,P_1,TILE_FRAC
     &,CANOPY,CATCH,CO2
     &,SM_LEVELS,DZSOIL,HCONS,F_TYPE
     &,HT,LAI,PAR,GPP,NPP,RESP_P
     &,ICE_FRACT,LAND_MASK,LYING_SNOW,PSTAR,Q_1
     &,QCF_1,QCL_1,RADNET_C,GC,RESIST
     &,ROOTD,SMC,SMVCCL,SMVCWT
     &,T_1,TIMESTEP,TI,TS1,TSTAR_GB
     &,TSTAR_TILE,U_1,V_1,U_0,V_0
     &,V_ROOT,V_SOIL,VFRAC
     &,Z0V_GB,Z0V,SIL_OROG,HO2R2_OROG,ZH
     &,Z1_UV,Z1_TQ,CANCAP,Z0MSEA,ALPHA1_GB,ALPHA1,ASHTF          
     &,BQ1_GB,BT1_GB,CD,CH
     &,FQW_1,FQW1_GB,FTL_1,FTL1_GB
     &,EPOT,EPOT_GB,FSMC,FSMC_GB
     &,E_SEA,H_SEA,FRACA,RESFS,F_SE
     &,RESFT,RESFT_GB,RHOKE,RHOKH_1,RHOKH_1_GB
     &,RHOKM_1_GB,RHOKPM,RHOKPM_GB,RHOKPM_POT,RHOKPM_POT_GB
     &,RIB_GB,RIB,TL_1,VSHR,Z0H_T,Z0M_T,Z0M_EFF_T,Z0M_EFF
     &,H_BLEND_OROG,H_BLEND,T1_SD,Q1_SD,TV1_SD,U_S,FB_SURF
     &,RHO_CD_MODV1,WIND_BLEND_FACTOR,HEAT_BLEND_FACTOR
     &,CDR10M,CHR1P5M,CER1P5M,FME
     &,SU10,SV10,SQ1P5,ST1P5,SFME
     &,RHO_ARESIST,ARESIST,RESIST_B,NRML
     &,L_Z0_OROG,L_RMBL,LTIMER
     &)

      IMPLICIT NONE

!  Input variables.  All fields are on P grid except where noted.
!  Fxxx in a comment indicates the file from which the data are taken.


!       GENERAL NOTES ABOUT GRID-DEFINITION INPUT VARIABLES.
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  For global data :-

!  An Arakawa B-grid is assumed in which each pole is represented by a
!  row of P-grid points. Entire fields of P-grid values are taken as
!  input, but  the two polemost rows are (a) not updated, in the case
!  of INOUT fields, or (b) set to zero, in the case of OUT fields.

!  If defined variable IBM is selected then land point calculations are
!  performed using the array LAND_INDEX to select land points. But note
!  that elements of LAND_INDEX define land points on the full field
!  (ie including polar rows).

!  Entire fields of UV-grid values are taken as input, but the two
!  polemost rows are (a) not updated, in the case of INOUT fields, or
!  (b) set to zero, in the case of OUT fields.

!  For limited-area data :-

!  The above applies, but for "polar rows", etc., read "rows at the
!  north and south boundaries of the area", etc.  E.g. if you want to
!  do calculations in UV-rows n to m inclusive, the input data will be
!  on P-rows n-1 to m+1, and UV-rows n-1 to m+1.  P-rows n to m will
!  then be updated.  Land specific variables are processed as for global
!  data.

!  For both cases, the following equalities apply amongst the input
!  grid-definition variables :-

!            P_POINTS = P_ROWS * ROW_LENGTH
!            U_POINTS = U_ROWS * ROW_LENGTH
!              U_ROWS = P_ROWS + 1
!            LAND_PTS <= P_POINTS

!  NB: All this has severe implications for batching/macrotasking;
!      effectively it can't be done on a shared-memory machine without
!      either rewriting this routine or using expensive synchronizations
!      (or other messy and/or undesirable subterfuges).


      LOGICAL LTIMER

      INTEGER                !    Variables defining grid.
     & P_POINTS              ! IN Number of P-grid points to be
!                               processed.
     &,P_FIELD               ! IN Total number of P-grid points.
     &,P1                    ! IN First P-point to be processed.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.
     &,LAND_FIELD            ! IN Total number of land points.
     &,N_TYPES               ! IN Number of tiles per land point.
     &,LAND_INDEX(LAND_FIELD)! IN Index for compressed land point array;
!                               ith element holds position in the FULL
!                               field of the ith land pt to be processed

      LOGICAL
     & GATHER                ! IN If true then leads variables are comp-
!                               ressed for sea-ice calculations. This
!                               saves duplicating calculations if there
!                               are a relatively few of sea-ice points.
!                               Set to false for a limited area run
!                               with a high proportion of sea-ice.

! Extra variables for the interactive stomatal resistance model

      INTEGER
     & SM_LEVELS             ! IN Number of soil moisture levels.
     &,F_TYPE(LAND_FIELD,N_TYPES)
!                              IN Plant functional type:
!                                1 - Broadleaf Tree
!                                2 - Needleleaf Tree
!                                3 - C3 Grass
!                                4 - C4 Grass
      REAL
     & CANOPY(LAND_FIELD)    ! IN Surface water (kg per sq metre). F642.
     &,CATCH(LAND_FIELD,N_TYPES)
!                              IN Surface capacity (max. surface water)
!                               (kg per sq metre).  F6416.
     &,CO2                   ! IN CO2 Mass Mixing Ratio
     &,DZSOIL(SM_LEVELS)     ! IN Soil layer thicknesses (m)
     &,HCONS(LAND_FIELD)     ! IN Soil thermal conductivity (W/m/K).
     &,HO2R2_OROG(LAND_FIELD)! IN Peak to trough height of unresolved
!                               orography devided by 2SQRT(2) (m).
     &,HT(LAND_FIELD,N_TYPES)! IN Canopy height (m).
     &,ICE_FRACT(P_FIELD)    ! IN Fraction of gridbox which is sea-ice.
     &,LAI(LAND_FIELD,N_TYPES)!IN Leaf area index.
     &,LYING_SNOW(P_FIELD)   ! IN Lying snow amount (kg per sq metre).
     &,PAR(P_FIELD)          ! IN Photosynthetically active radiation
!                               (W/m2).
     &,PSTAR(P_FIELD)        ! IN Surface pressure (Pascals).
     &,P_1(P_FIELD)          ! IN pressure lowest atmospheric
     &,Q_1(P_FIELD)          ! IN Specific humidity for lowest
!                               atmospheric layer (kg water per kg air).
     &,QCF_1(P_FIELD)        ! IN Cloud ice for lowest atmospheric layer
!                               (kg water per kg air).
     &,QCL_1(P_FIELD)        ! IN Cloud liquid water for lowest atm
!                               layer (kg water per kg air).
     &,RESIST(LAND_FIELD,N_TYPES)
!                            ! IN "Stomatal" resistance to evaporation
!                               (seconds per metre).  F6415.
     &,ROOTD(LAND_FIELD,N_TYPES)
!                              IN "Root depth" (metres).  F6412.
     &,SIL_OROG(LAND_FIELD)  ! IN Silhouette area of unresolved
!                               orography per unit horizontal area
     &,SMC(LAND_FIELD,N_TYPES)!IN Soil moisture content (kg/m2). F621.
     &,SMVCCL(LAND_FIELD)    ! IN Critical volumetric SMC (cubic metres
!                               per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_FIELD)    ! IN Volumetric wilting point (cubic m of
!                               water per cubic m of soil).  F6231.

!    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
!          moisture content.                            ~~~~~~~~~

      REAL                   !    (Split to avoid > 19 continuations.)
     & T_1(P_FIELD)          ! IN Temperature for lowest atmospheric
!                               layer (Kelvin).
     &,TILE_FRAC(P_FIELD,N_TYPES)
!                              IN Fractional coverage for each tile
     &,TIMESTEP              ! IN Timestep in seconds for EPDT calc.
     &,TI(P_FIELD)           ! IN Temperature of sea-ice surface layer
!                               (Kelvin)
     &,TL_1(P_FIELD)         ! IN Liquid/frozen water temperature for
!                               lowest atmospheric layer (K).
     &,TS1(LAND_FIELD)       ! IN Temperature of top soil layer (K)
     &,TSTAR_TILE(P_FIELD,N_TYPES)
!                              IN Tile surface temperature (K).
     &,TSTAR_GB(P_FIELD)     ! IN Mean gridbox surface temperature (K).
     &,U_0(P_FIELD)          ! IN West-to-east component of ocean
!                               surface current (m/s; ASSUMED zero over
!                               land). UV grid.  F615.
     &,U_1(P_FIELD)          ! IN West-to-east wind component for lowest
!                               atmospheric layer (m/s).  On UV grid.
     &,V_0(P_FIELD)          ! IN South-to-north component of ocean
!                               surface current (m/s; ASSUMED zero over
!                               land). UV grid.  F616.
     &,V_1(P_FIELD)          ! IN South-to-north wind component for
!                               lowest atm. layer (m/s).  On UV grid.
     &,V_ROOT(LAND_FIELD,N_TYPES)
!                            ! IN Volumetric soil moisture concentration
!                               in the rootzone (m3 H2O/m3 soil).
     &,V_SOIL(LAND_FIELD)    ! IN Volumetric soil moisture concentration
!                               in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_FIELD,N_TYPES)
!                            ! IN Vegetation fraction.
     &,Z0V(P_FIELD,N_TYPES)  ! IN Tile vegetative roughness length (m).
     &,Z0V_GB(P_FIELD)       ! IN Gridbox veg. roughness length (m).
     &,Z1_UV(P_FIELD)        ! IN Height of lowest uv level (m).
     &,Z1_TQ(P_FIELD)        ! IN Height of lowest tq level (m).
!                               Note, if the grid used is staggered in
!                               the vertical, Z1_UV and Z1_TQ can be
!                               different.
     &,ZH(P_FIELD)           ! IN Height of top of boundary layer (m).

      LOGICAL
     & LAND_MASK(P_FIELD)    ! IN .TRUE. for land; .FALSE. elsewhere.
!                               F60.
     &,SU10                  ! IN STASH flag for 10-metre W wind.
     &,SV10                  ! IN STASH flag for 10-metre S wind.
     &,SQ1P5                 ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                 ! IN STASH flag for 1.5-metre temperature.
     &,SFME                  ! IN STASH flag for wind mixing energy flux
     &,L_RMBL                ! IN T to use rapidly mixing boundary
!                               scheme in IMPL_CAL
     &,L_Z0_OROG             ! IN .TRUE. to use orographic roughness.

!  Modified (INOUT) variables.

      REAL
     & CANCAP(P_FIELD,N_TYPES)! INOUT Volumetric heat capacity of
C                            !       vegetation canopy (J/Kg/m3).
     &,RADNET_C(P_FIELD,N_TYPES) ! INOUT Adjusted net radiation for 
C                            !          vegetation over land (W/m2).
     &,Z0MSEA(P_FIELD)       ! INOUT Sea-surface roughness length for   
!                                momentum (m).  F617.
     &,GC(LAND_FIELD,N_TYPES)! INOUT "Stomatal" conductance to
!                                evaporation (m/s).

!  Output variables.
!
      REAL
     & ALPHA1(P_FIELD,N_TYPES)!OUT Gradients of saturated specific
!                                humidity with respect to temperature
!                                between the bottom model layer and tile
!                                surface
     &,ALPHA1_GB(P_FIELD)    ! OUT Gradient of saturated specific
!                                humidity with respect to temperature
!                                between the bottom model layer and the
!                                mean surface
     &,ASHTF(P_FIELD)        ! OUT Coefficient to calculate surface
!                                heat flux into soil or sea-ice (W/m2/K)
     &,BQ1_GB(P_FIELD)       ! OUT A buoyancy parameter for lowest atm
!                                level ("beta-q twiddle").
     &,BT1_GB(P_FIELD)       ! OUT A buoyancy parameter for lowest atm
!                                level ("beta-T twiddle").
     &,CD(P_FIELD)           ! OUT Bulk transfer coefficient for
!                                momentum.
     &,CH(P_FIELD)           ! OUT Bulk transfer coefficient for heat
!                                and/or moisture.
     &,CDR10M(P_FIELD)       ! OUT Reqd for calculation of 10m wind
!                                (u & v).
!                                NBB: This is output on the UV-grid, but
!                                with the first and last rows set to a
!                                "missing data indicator".
!                                Sea-ice leads ignored. See 3.D.7 below.
     &,CHR1P5M(P_FIELD)      ! OUT Reqd for calculation of 1.5m temp.
!                                Sea-ice leads ignored. See 3.D.7 below.
     &,CER1P5M(P_FIELD)      ! OUT Reqd for calculation of 1.5m sp
!                                humidity. Sea-ice leads ignored.
!                                See 3.D.7 below.
     &,E_SEA(P_FIELD)        ! OUT Evaporation from sea times leads
!                                fraction (kg/m2/s). Zero over land.
     &,FME(P_FIELD)          ! OUT Wind mixing energy flux (Watts/sq m).
     &,F_SE(P_FIELD,N_TYPES) ! OUT Fraction of the evapotranspiration
!                                which is bare soil evaporation.
     &,EPOT(P_FIELD,N_TYPES) ! OUT "Explicit" potential evaporation
!                                on P-grid (kg/m2/s).
     &,EPOT_GB(P_FIELD)      ! OUT "Explicit" potential evaporation
!                                on P-grid (kg/m2/s)
!                                for whole grid box.
     &,FSMC(LAND_FIELD,N_TYPES)
!                              OUT soil moisture availability.
     &,FSMC_GB(LAND_FIELD)     
!                              OUT soil moisture availability
!                                  for whole grid box.
     &,FQW_1(P_FIELD,N_TYPES)! OUT "Explicit" surface flux of QW (i.e.
!                                 evaporation), on P-grid (kg/m2/s).
     &,FQW1_GB(P_FIELD)      ! OUT "Explicit" surface flux of QW (i.e.
!                                evaporation), on P-grid (kg/m2/s). for
!                                whole grid-box
     &,FTL_1(P_FIELD,N_TYPES)! OUT "Explicit" surface flux of TL = H/CP.
!                                (sensible heat / CP).
     &,FTL1_GB(P_FIELD)      ! OUT "Explicit" surface flux of TL = H/CP.
!                                (sensible heat / CP). grid-box mean
     &,FRACA(P_FIELD,N_TYPES)! OUT Fraction of surface moisture flux
!                                with only aerodynamic resistance.
     &,GPP(LAND_FIELD,N_TYPES)!OUT Gross Primary Productivity
!                               (kg C/m2/s).
     &,H_BLEND(P_FIELD)      ! OUT Blending height for tiles
     &,H_BLEND_OROG(P_FIELD) ! OUT Blending height for orographic
!                                roughness
     &,H_SEA(P_FIELD)        ! OUT Surface sensible heat flux over sea
!                                times leads fraction (W/m2).
!                                Zero over land.
     &,NPP(LAND_FIELD,N_TYPES)!OUT Net Primary Productivity (kg C/m2/s).
     &,Q1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                                fluctuations of surface layer
!                                specific humidity (kg/kg).
     &,RESFS_GB(P_FIELD)     ! OUT Combined soil, stomatal and
!                                aerodynamic resistance factor =
!                                PSIS/(1+RS/RA) for fraction (1-FRACA)
     &,RESFT_GB(P_FIELD)     ! OUT Total resistance factor
!                                FRACA+(1-FRACA)*RESFS.
     &,RESP_P(LAND_FIELD,N_TYPES)
!                            ! OUT Plant respiration rate (kg C/m2/s).
     &,RIB_GB(P_FIELD)       ! OUT Mean bulk Richardson number for
!                                lowest layer
     &,T1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                                fluctuations of surface layer
!                                temperature (K).
     &,TV1_SD(P_FIELD)       ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     virtual temperature (K).
     &,U_S(P_FIELD)          ! OUT Surface friction velocity (m/s)
     &,FB_SURF(P_FIELD)      ! OUT Surface flux buoyancy over density
!                            !     (m^2/s^3)
!
     &,VSHR(P_FIELD)         ! OUT Magnitude of surface-to-lowest-level
!                                wind
     &,Z0H(P_FIELD)          ! OUT Roughness length for heat & moisture
     &,Z0M(P_FIELD)          ! OUT Roughness length for momentum (m).
     &,Z0M_EFF(P_FIELD)      ! OUT Effective roughness length for
!                                momentum
     &,RHO_ARESIST(P_FIELD)  ! OUT, RHOSTAR*CD_STD*VSHR  for SCYCLE
     &,ARESIST(P_FIELD)      ! OUT, 1/(CD_STD*VSHR)      for SCYCLE
     &,RESIST_B(P_FIELD)     ! OUT, (1/CH-1/CD_STD)/VSHR for SCYCLE


! Surface exchange coefficients;passed to subroutine IMPL_CAL
      REAL
     & RHO_CD_MODV1(P_FIELD) ! OUT rhostar*cD*vshr before horizontal
!                                interpolation output as a diagnostic.
     &,RHOKE_GB(P_FIELD)     ! OUT For FQW
     &,RHOKH_1(P_FIELD,N_TYPES)
!                            ! OUT For FTL
     &,RHOKH_1_GB(P_FIELD)   ! OUT For FTL
     &,RHOKM_1_GB(P_FIELD)   ! OUT For momentum. NB: This is output on
!                                UV-grid, but with the first and last
!                                rows set to a "missing data indicator".
     &,RHOKPM_GB(P_FIELD)    ! OUT Mixing coefficient for Penman-
!                                Monteith scheme
     &,RHOKPM_POT(P_FIELD,N_TYPES)
!                              OUT Surface exchange coeff. for
!                                potential evaporation.
     &,RHOKPM_POT_GB(P_FIELD)! OUT Surface exchange coeff. for
!                                potential evaporation
!                                for whole grid box.

      INTEGER
     & NRML(P_FIELD)         ! OUT 1 if surface layer unstable, else 0.

!  Symbolic constants ------------------------------------------------

!   (1) UM-wide common parameters.

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------


!   (2) Boundary Layer local parameters.

!!!-----------COMDECK C_CHARNK FOR SUBROUTINE SF_EXCH----------
! CHARNOCK is a constant in the Charnock formula for sea-surface
!          roughness length for momentum (Z0MSEA).
      REAL CHARNOCK

      PARAMETER(CHARNOCK = 0.011)
!!----------------------------------------------------------------------
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C*L-----------COMDECK C_HT_M FOR SUBROUTINE SF_EXCH----------
C Z10M  = height of 10m level for diagnostic calculations (m).
C Z1P5M = height of 1.5m level for diagnostic calculations (m).
      REAL Z10M,Z1P5M

      PARAMETER(Z10M  = 10.0,
     &          Z1P5M = 1.5)
C*----------------------------------------------------------------------
!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0HSEA = 4.0E-5,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
!!----------------------------------------------------------------------
      REAL    RI_CRIT   ! Critical Richardson number, where Z0M_EFF=Z0M.
!                       ! Linear interpolation between RIB=0 and RI_CRIT
                                                                       
      REAL    OROG_DRAG_PARAM    ! Tunable parameter in calculation of
!                                ! Effective roughness length for 
!                                ! momentum                     
      PARAMETER(                                                
     & RI_CRIT=0.5,                                              
     & OROG_DRAG_PARAM=0.3)                                         
!*----------------------------------------------------------------------
C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
      REAL KAPPAI
      PARAMETER (
     + KAPPAI=2.09          ! Thermal conductivity of sea-ice (W per
C                           ! m per K).
     +)
      REAL DE
      PARAMETER (
     + DE = 0.1             ! Effective thickness of sea-ice surface
C                           ! layer (m).
     +)
C-----------------------------------------------------------------------
C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------


!   (3) Derived local parameters.

      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP,H_BLEND_MIN,H_BLEND_MAX

      PARAMETER (
     & ETAR=1./(1.-EPSILON)  ! Used in calc of buoyancy parameter BETAC.
     &,GRCP=G/CP             ! Used in calc of dT across surface layer.
     &,LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LFRCP=LF/CP           ! Freezing-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     &,H_BLEND_MIN=0.0       ! Minimum blending height.
     &,H_BLEND_MAX=1000.0    ! Maximum blending height (m).
     &)


!   External subprograms called.

      EXTERNAL SF_ROUGH,SF_LBEST,SF_RIB,FCDCH,QSAT,SFL_INT,SF_RESIST,
     &         SF_FLUX,SF_STOM,TIMER


!   Define local storage.

!   (a) Workspace.

!  Workspace ---------------------------------------------------------
      REAL
     & BQ_1(P_FIELD,N_TYPES)!A buoyancy parameter for lowest atm level
!                                ("beta-q twiddle").
     &,BT_1(P_FIELD,N_TYPES)!A buoyancy parameter for lowest atm level.
!                                ("beta-T twiddle").
     &,CD_LEAD(P_FIELD)     ! Bulk transfer coefficient for momentum
!                              over sea-ice leads.Missing data over non
!                              sea-ice points.(Temporary store for
!                              Z0MIZ)
     &,CD_MIZ(P_FIELD)      ! Bulk transfer coefficient for momentum
!                              over the sea-ice Marginal Ice Zone.
!                              Missing data indicator over non sea-ice.
     &,CD_STD_T(P_FIELD,N_TYPES)
!                             Local drag coefficient for
!                              calculation of interpolation coefficients
     &,CD_STD(P_FIELD)    ! Local drag coefficient for
!                         !  calculation of interpolation coefficients
     &,CD_T(P_FIELD,N_TYPES)! Drag coefficient on tile
     &,CH_LEAD(P_FIELD)     ! Bulk transfer coefficient for heat and
!                              or moisture over sea ice leads.
!                              Missing data indicator over non sea-ice.
     &,CH_MIZ(P_FIELD)      ! Bulk transfer coefficient for heat and
!                              or moisture over the Marginal Ice Zone.
!                              Missing data indicator over non sea-ice.
     &,CH_T(P_FIELD,N_TYPES)! Transfer coefficient for heat and
!                              moisture on tile
     &,DQ(P_FIELD,N_TYPES)  ! Sp humidity difference between surface
!                              and lowest atmospheric level (Q1 - Q*).
!                              Holds value over sea-ice where ICE_FRACT
!                              >0 i.e. Leads contribution not included.
     &,DQ_LEAD(P_FIELD)     ! DQ for leads fraction of gridsquare.
!                              Missing data indicator over non sea-ice.
     &,DTEMP(P_FIELD,N_TYPES)!Liquid/ice static energy difference
!                              between surface and lowest atmospheric
!                              level, divided by CP (a modified
!                              temperature difference).
!                              Holds value over sea-ice where ICE_FRACT
!                              >0 i.e. Leads contribution not included.
     &,DTEMP_LEAD(P_FIELD)  ! DTEMP for leads fraction of gridsquare.
!                              Missing data indicator over non sea-ice.
     &,EPDT(P_FIELD)        ! "Potential" Evaporation * Timestep
     &,HEAT_BLEND_FACTOR(P_FIELD)
!                             used in estimation of heat and
!                              moisture at blending height
     &,NL0(LAND_FIELD)      ! Nitrogen concentration of the top leaf
!                            (kg N/kg C).
     &,PSIS(P_FIELD,N_TYPES)! Soil moisture availability factor.
     &,PSTAR_ICE(P_FIELD)   ! Surface pressure over sea ice (Pa).
     &,Q_BLEND(P_FIELD)     ! Estimate of blending height Q
     &,QS_BLEND(P_FIELD)    ! Sat. specific humidity
!                              qsat(TL_BLEND,PSTAR)
     &,QW_BLEND(P_FIELD)    ! Estimate of blending height Q
     &,QS1(P_FIELD)         ! Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSL(P_FIELD)         ! Saturated sp humidity at liquid/ice
!                              temperature and pressure of lowest
!                              atmospheric level.
     &,QSTAR_GB(P_FIELD)    ! Gridbox mean QSTAR
     &,QSTAR(P_FIELD)       ! Surface saturated sp humidity. Holds
!                              value over sea-ice where ICE_FRACT > 0.
!                              i.e. Leads contribution not included.
     &,QSTAR_LEAD(P_FIELD)  ! QSTAR for sea-ice leads.
!                              Missing data indicator over non sea-ice.
     &,RA(P_FIELD)          ! Aerodynamic resistance.
     &,RESFS(P_FIELD,N_TYPES)!Combined soil, stomatal and aerodynamic
!                              resistance factor = PSIS/(1+RS/RA) for
!                              fraction (1-FRACA)
     &,RESFT(P_FIELD,N_TYPES)!Total resistance factor
!                              FRACA+(1-FRACA)*RESFS.
     &,RHOKE(P_FIELD,N_TYPES)!For FQW
     &,RHOKM_1(P_FIELD,N_TYPES)
!                             RHOKM for tile
     &,RHOKPM(P_FIELD,N_TYPES)
!                             Mixing coefficient
     &,RHOSTAR(P_FIELD,N_TYPES)
!                             Surface air density in kg per cubic metre.
     &,RHOSTAR_GB(P_FIELD)  ! Surface air density in kg per cubic metre.
     &,DB_GB(P_FIELD)       ! Gridbox mean buoyancy difference.
     &,DB_LEAD(P_FIELD)     ! Buoyancy difference for lead part of grdbx
     &,DB(P_FIELD,N_TYPES)  ! Buoyancy difference for surface tile
     &,RIB_LEAD(P_FIELD)    ! Bulk Richardson no. for sea-ice leads at
!                             lowest layer. At non sea-ice points holds
!                             RIB for FCDCH calculation, then set to
!                             to missing data indicator.
     &,RIB(P_FIELD,N_TYPES) ! Bulk Richardson no. for surface tile
     &,ROOT(LAND_FIELD)     ! Root biomass (kg C/m2).
     &,T_BLEND(P_FIELD)     ! Estimate of blending height T
     &,TL_BLEND(P_FIELD)    ! Estimate of blending height TL
     &,TSTAR_NL(P_FIELD)    ! TSTAR No Leads: surface temperature
!                              over sea-ice fraction of gridsquare.
!                              =TSTAR over non sea-ice points.
     &,U_BLEND(P_FIELD)     ! Estimate of blending height U
     &,V_BLEND(P_FIELD)     ! Estimate of blending height V
     &,WIND_BLEND_FACTOR(P_FIELD)
!                              used in estimation of winds at
!                              blending height
     &,WIND_PROFILE_FACTOR(P_FIELD,N_TYPES)
!                              For transforming effective surface
!                              transfer coefficients to those excluding
!                              form drag.
     &,RECIP_L_MO(P_FIELD,N_TYPES)
!                           ! Reciprocal of the Monin-Obukhov length.
     &,V_S(P_FIELD,N_TYPES) ! Surface scaling velocity (friction velocit
!                           ! modified with convective turbulence
!                           ! velocity) including orographic form drag
!                           ! effects.
     &,V_S_STD(P_FIELD,N_TYPES)
!                           ! Surface scaling velocity (friction velocit
!                           ! modified with convective turbulence
!                           ! velocity) excluding orographic form drag
!                           ! effects.
     &,V_S_LEAD(P_FIELD)    ! Surface scaling velocity (friction velocit
!                           ! modified with convective turbulence veloci
!                           ! for leads part of sea gridbox.
     &,Z0HS(P_FIELD)        ! Roughness length for heat and moisture
!                              transport over sea.
     &,Z0M_EFF_T(P_FIELD,N_TYPES)
!                             Effective roughness length for momentum
     &,Z0H_T(P_FIELD,N_TYPES)!Tile roughness length for heat and
!                              moisture
     &,Z0M_T(P_FIELD,N_TYPES)!Local tileroughness length for momentum

!  Workspace (reqd for compression).
      INTEGER
     & SICE_INDEX(P_FIELD)   ! Index vector for gather to sea-ice points

      LOGICAL ITEST(P_FIELD) !Used as 'logical' for compression.


!   (b) Scalars.

      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,ITILE       ! Loop counter (tile index).
     &,J,K         ! Offset counter within I-loop.
     &,L,N         ! Loop counter (land point field index).
     &,NSICE       ! Number of sea-ice points.
     &,SI          ! Loop counter (sea-ice field index).
      REAL
     & TAU         ! Magnitude of surface wind stress over sea.
     &,W_S_CUBED   ! Cube of surface layer free convective scaling
!                  ! velocity
     &,W_M         ! Turbulent velocity scale for surface layer

      LOGICAL
     & L_LAND      ! a logical

! Extra work variables for the canopy (stomatal) conductance model.
      LOGICAL
     & INT_STOM              ! T for interactive stomatal resistance.
      PARAMETER (INT_STOM=.TRUE.)


!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',3)
      ENDIF


!-----------------------------------------------------------------------
!!  1.  Construct SICE_INDEX for compression onto sea points in
!!      sea-ice leads calculations.
!-----------------------------------------------------------------------

        DO I=P1,P1+P_POINTS-1
          ITEST(I) = .FALSE.
          IF (ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I))
     &      ITEST(I) = .TRUE.
        ENDDO

        NSICE = 0
        DO I=P1,P1+P_POINTS-1
          IF(ITEST(I))THEN
            NSICE = NSICE + 1
            SICE_INDEX(NSICE) = I
          END IF
        ENDDO

!-----------------------------------------------------------------------
!!  2.  Calculate QSAT values required later and components of ocean
!!      current.
!!       Done here to avoid loop splitting.
!!       QSTAR 'borrowed' to store P at level 1 (just this once).
!!       PSIS 'borrowed' to store leads and non sea-ice surface temp.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!!  2.1 IF (GATHER) THEN
!!       Calculate temperatures and pressures for QSAT calculations.
!!       Calculate QSAT values. For sea-ice points, separate values
!!       are required for the leads (QSTAR_LEAD) and sea-ice (QSTAR)
!!       fractions respectively. QSTAR_LEAD = missing data, elsewhere.
!!       Use RS to store compressed PSTAR for this section only.
!!       NB Unlike QSTAR, TSTAR values at sea-ice points are gridsq.
!!       means and so include the leads contribution.
!!      ELSE
!!       As above with QSTAR_LEAD done on full field.
!!      ENDIF
!-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=P1,P1+P_POINTS-1
          TSTAR_NL(I) = TSTAR_GB(I)
          QSTAR_LEAD(I) = 1.0E30                ! Missing data indicator
        ENDDO
        IF (NSICE.GT.0) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
          DO SI = 1,NSICE
            I = SICE_INDEX(SI)

            TSTAR_NL(I) = (TSTAR_GB(I)-(1.0-ICE_FRACT(I)) *TFS)
     &                    / ICE_FRACT(I)                     ! P2430.1
            PSIS(SI,1) = TFS
            PSTAR_ICE(SI) = PSTAR(I)
          ENDDO
        ENDIF

        CALL QSAT(QSL(P1),TL_1(P1),P_1(P1),P_POINTS)
        CALL QSAT(QS1(P1),TL_1(P1),PSTAR(P1),P_POINTS)

        CALL QSAT(QSTAR(P1),TSTAR_NL(P1),PSTAR(P1),P_POINTS)
!            ...values at sea-ice points contain ice contribution only
        IF (NSICE.GT.0) CALL QSAT(QSTAR_LEAD,PSIS,PSTAR_ICE,NSICE)
!            ...values at sea-ice points only

        CALL QSAT(QSTAR_GB(P1),TSTAR_GB(P1),PSTAR(P1),P_POINTS)
!            ...values at sea-ice points gb-average

      ELSE

        DO I=P1,P1+P_POINTS-1
          TSTAR_NL(I) = TSTAR_GB(I)
! Set to missing data at non sea-ice points after QSAT.
          PSIS(I,1) = TSTAR_GB(I)
          IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
            TSTAR_NL(I) = (TSTAR_GB(I)-(1.0-ICE_FRACT(I)) *TFS)
     &                / ICE_FRACT(I)                          ! P2430.1
            PSIS(I,1) = TFS
          ENDIF
        ENDDO
        CALL QSAT(QSL(P1),TL_1(P1),P_1(P1),P_POINTS)
        CALL QSAT(QS1(P1),TL_1(P1),PSTAR(P1),P_POINTS)

        CALL QSAT(QSTAR(P1),TSTAR_NL(P1),PSTAR(P1),P_POINTS)
!          ...values at sea-ice points contain ice contribution only

        IF (NSICE.GT.0)
     &       CALL QSAT(QSTAR_LEAD(P1),PSIS(P1,1),PSTAR(P1),P_POINTS)
!          ...values at sea-ice points contain leads contribution only

        CALL QSAT(QSTAR_GB(P1),TSTAR_GB(P1),PSTAR(P1),P_POINTS)
!            ...values at sea-ice points gb-average

        DO I=P1,P1+P_POINTS-1
          IF ( .NOT.(ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) )
     &      QSTAR_LEAD(I) = 1.0E30
        ENDDO
      ENDIF                ! End of IF (GATHER) THEN... ELSE.

!-----------------------------------------------------------------------
!!  2.2  Reset aggregated quantities
!-----------------------------------------------------------------------

      DO I=1,P_FIELD
        RHO_ARESIST(I) = 0.0
        ARESIST(I) = 0.0
        RESIST_B(I) = 0.0
        EPOT_GB(I) = 0.0
        FQW1_GB(I)=0.0
        FTL1_GB(I)=0.0
        RIB_GB(I)=0.0
        DB_GB(I)=0.0
        BT1_GB(I)=0.0
        BQ1_GB(I)=0.0
        RESFS_GB(I)=0.0
        RESFT_GB(I)=0.0
        ALPHA1_GB(I) = 0.0
        RHOKE_GB(I) = 0.0
        CD(I)=0.0
        CD_STD(I)=0.0
        CH(I)=0.0
        RHOKH_1_GB(I) = 0.0
        RHOKM_1_GB(I) = 0.0
        RHOKPM_GB(I) = 0.0
        RHOKPM_POT_GB(I) = 0.0
        T1_SD(I)=0.0
        Q1_SD(I)=0.0
        TV1_SD(I)=0.0
        RHOSTAR_GB(I)=0.0
        NRML(I) = 0

        DO ITILE=1,N_TYPES
           DB(I,ITILE)=0.0
           RIB(I,ITILE)=0.0
        ENDDO
      ENDDO

      DO L=1,LAND_FIELD
        FSMC_GB(L) = 0.0
      ENDDO

!-----------------------------------------------------------------------
!!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!  3.1 Calculate neutral roughness lengths and blending height for
!!      surface
!-----------------------------------------------------------------------

! Grid box mean value for estimating model values at bending height

      L_LAND=.FALSE.  ! Calc over all points)

      CALL SF_ROUGH (
     & P_FIELD,P_POINTS,LAND_FIELD,LAND_PTS,LAND_MASK,L_LAND,P1,LAND1,
     & LAND_INDEX,
     & L_Z0_OROG,Z1_UV,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V_GB,SIL_OROG,HO2R2_OROG,RIB_GB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR(1,1),H_BLEND_OROG,CD_LEAD,Z0HS,
     & LTIMER)


! Estimate model values at blending height from neutral profile

      CALL SF_LBEST (
     & P_POINTS,P_FIELD,P1,H_BLEND_OROG,
     & QCL_1,QCF_1,QSTAR_GB,Q_1,TSTAR_GB,T_1,U_1,V_1,
     & Z0M_EFF,Z0H,Z0M,Z1_UV,Z1_TQ,H_BLEND,HEAT_BLEND_FACTOR,
     & Q_BLEND,QW_BLEND,T_BLEND,TL_BLEND,U_BLEND,V_BLEND,
     & WIND_BLEND_FACTOR,LTIMER
     & )

! Calc. QSAT at blending height
      CALL QSAT(QS_BLEND(P1),TL_BLEND(P1),PSTAR(P1),P_POINTS)


! Calc QSTAR_no_leads and store in QSTAR_GB
      CALL QSAT(QSTAR_GB(P1),TSTAR_NL(P1),PSTAR(P1),P_POINTS)


! Start of loop over tiles
      DO ITILE=1,N_TYPES

!-----------------------------------------------------------------------
!  3.1.1 Tile roughnesses
!-----------------------------------------------------------------------

! Only calculate roughnesses for sea points once

        L_LAND=.FALSE.

        CALL SF_ROUGH (
     &   P_FIELD,P_POINTS,LAND_FIELD,LAND_PTS,LAND_MASK,L_LAND,P1,LAND1,
     &   LAND_INDEX,
     &   L_Z0_OROG,Z1_UV,Z0MSEA,ICE_FRACT,
     &   LYING_SNOW,Z0V(1,ITILE),SIL_OROG,HO2R2_OROG,RIB(1,ITILE),
     &   Z0M_EFF_T(1,ITILE),Z0M_T(1,ITILE),Z0H_T(1,ITILE),
     &   WIND_PROFILE_FACTOR(1,ITILE),H_BLEND_OROG,CD_LEAD,Z0HS,
     &   LTIMER
     &   )

!-----------------------------------------------------------------------
!!  3.2 Calculate buoyancy parameters and bulk Richardson number for
!!      the lowest model level.
!-----------------------------------------------------------------------


! Tile temperature passed to sf_rib through tstar_nl
        DO I=P1,P1+P_POINTS-1
          IF ( LAND_MASK(I) ) TSTAR_NL(I)=TSTAR_TILE(I,ITILE)
        ENDDO


        CALL QSAT(QSTAR(P1),TSTAR_NL(P1),PSTAR(P1),P_POINTS)


! qstar over sea-ice doesn not include leads

        DO I=P1,P1+P_POINTS-1
          IF(.NOT.LAND_MASK(I)) QSTAR(I)=QSTAR_GB(I)
        ENDDO

        CALL SF_RIB (
     &   P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,LAND_MASK,L_LAND,INT_STOM,
     &   P1,LAND1,
     &   GATHER,LAND_INDEX,
     &   NSICE,SICE_INDEX,ICE_FRACT,Q_BLEND,QW_BLEND,QCL_1,QCF_1,
     &   T_BLEND,TL_BLEND,QSL,QSTAR,QSTAR_LEAD,
     &   QS_BLEND,TSTAR_NL,Z1_TQ,Z1_UV,Z0M_EFF_T(1,ITILE),
     &   Z0M_T(1,ITILE),Z0H_T(1,ITILE),Z0HS,Z0MSEA,
     &   WIND_PROFILE_FACTOR(1,ITILE),U_BLEND,U_0,V_BLEND,V_0,
     &   ROOTD(1,ITILE),SMVCCL,SMVCWT,SMC(1,ITILE),VFRAC(1,ITILE),
     &   V_SOIL,CANOPY,CATCH(1,ITILE),
     &   LYING_SNOW,GC(1,ITILE),RESIST(1,ITILE),
     &   DB(1,ITILE),DB_LEAD,RIB(1,ITILE),RIB_LEAD,PSIS(1,ITILE),VSHR,
     &   ALPHA1(1,ITILE),BT_1(1,ITILE),BQ_1(1,ITILE),
     &   FRACA(1,ITILE),RESFS(1,ITILE),
     &   DQ(1,ITILE),DQ_LEAD,DTEMP(1,ITILE),DTEMP_LEAD,LTIMER
     &   )

!-----------------------------------------------------------------------
!!  3.3 Calculate stability corrected effective roughness length.
!!  Simple linear interpolation when RIB between 0 and RIB_CRIT (>0) for
!!  form drag term.
!-----------------------------------------------------------------------


! Stability correction only applies to land points
        L_LAND = .TRUE.

        CALL SF_ROUGH (
     &   P_FIELD,P_POINTS,LAND_FIELD,LAND_PTS,LAND_MASK,L_LAND,P1,LAND1,
     &   LAND_INDEX,
     &   L_Z0_OROG,Z1_UV,Z0MSEA,ICE_FRACT,
     &   LYING_SNOW,Z0V(1,ITILE),SIL_OROG,HO2R2_OROG,RIB(1,ITILE),
     &   Z0M_EFF_T(1,ITILE),Z0M_T(1,ITILE),Z0H_T(1,ITILE),
     &   WIND_PROFILE_FACTOR(1,ITILE),H_BLEND_OROG,CD_LEAD,Z0HS,
     &   LTIMER
     &   )

      ENDDO ! n_types


      DO ITILE=1,N_TYPES

! Calculate 'mean' richardson number for mean roughness lengths
        DO I=P1,P1+P_POINTS-1

          IF (.NOT.LAND_MASK(I)) THEN
            RIB(I,ITILE) = RIB(I,1)
            DB(I,ITILE) = DB(I,1)
          ENDIF

          RIB_GB(I) = RIB_GB(I) + RIB(I,ITILE) * TILE_FRAC(I,ITILE)
          DB_GB(I) = DB_GB(I) + DB(I,ITILE) * TILE_FRAC(I,ITILE)

        ENDDO !End of p_point loop


        DO I=P1,P1+P_POINTS-1
          IF (.NOT. LAND_MASK(I).AND.ITILE.GT.1) THEN
             PSIS(I,ITILE)=PSIS(I,1)
             DQ(I,ITILE)=DQ(I,1)
             DTEMP(I,ITILE)=DTEMP(I,1)
             FRACA(I,ITILE)=FRACA(I,1)
             RESFS(I,ITILE)=RESFS(I,1)
             BT_1(I,ITILE)=BT_1(I,1)
             BQ_1(I,ITILE)=BQ_1(I,1)
             ALPHA1(I,ITILE)=ALPHA1(I,1)
             Z0M_EFF_T(I,ITILE) = Z0M_EFF_T(I,1)
             Z0M_T(I,ITILE) = Z0M_T(I,1)
             Z0H_T(I,ITILE) = Z0H_T(I,1)
             WIND_PROFILE_FACTOR(I,ITILE) = WIND_PROFILE_FACTOR(I,1)
          ENDIF
        ENDDO !P_POINTS
      ENDDO !loop over tiles

! stability correction for grid box roughness lengths

      CALL SF_ROUGH (
     & P_FIELD,P_POINTS,LAND_FIELD,LAND_PTS,LAND_MASK,L_LAND,P1,LAND1,
     & LAND_INDEX,
     & L_Z0_OROG,Z1_UV,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V_GB,SIL_OROG,HO2R2_OROG,RIB_GB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR(1,1),H_BLEND_OROG,CD_LEAD,Z0HS,
     & LTIMER)


!-----------------------------------------------------------------------
!!  3.4 Calculate CD, CH via routine FCDCH.
!!  Calculate CD_MIZ,CH_MIZ,CD_LEAD,CH_LEAD on full field then set
!!  non sea-ice points to missing data (contain nonsense after FCDCH)
!   Unlike the QSAT calculations above, arrays are not compressed to
!   sea-ice points for FCDCH. This is because it would require extra
!   work space and initial tests showed that with with the extra
!   compression calculations required no time was saved.
!   NB CD_LEAD stores Z0MIZ for calculation of CD_MIZ,CH_MIZ.
!-----------------------------------------------------------------------

      L_LAND=.FALSE.

      CALL FCDCH(P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,DB_GB,VSHR,
     &           CD_LEAD,CD_LEAD,ZH,Z1_UV,Z1_TQ,
     &           WIND_PROFILE_FACTOR(1,1),
     &           CD_MIZ,CH_MIZ,CD_STD_T(1,1),V_S(1,1),V_S_STD(1,1),
     &           RECIP_L_MO(1,1),LTIMER)
!                                           ! Marginal Ice Zone.P2430.9
!
      CALL FCDCH(P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,DB_LEAD,VSHR,
     &           Z0MSEA,Z0HS,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR(1,1),
     &           CD_LEAD,CH_LEAD,CD_STD_T(1,1),V_S_LEAD,V_S_STD(1,1),
     &           RECIP_L_MO(1,1),LTIMER)
!                                           ! Sea-ice leads.P2430.8

      DO ITILE=1,N_TYPES

        IF (ITILE.EQ.1) THEN
          L_LAND=.FALSE.
        ELSE
          L_LAND=.TRUE.
        ENDIF

        CALL FCDCH(P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,DB(1,ITILE),
     &             VSHR,Z0M_EFF_T(1,ITILE),Z0H_T(1,ITILE),ZH,
     &             Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR(1,ITILE),
     &             CD_T(1,ITILE),CH_T(1,ITILE),CD_STD_T(1,ITILE),
     &             V_S(1,ITILE),V_S_STD(1,ITILE),RECIP_L_MO(1,ITILE),
     &             LTIMER)


        DO I=P1,P1+P_POINTS-1

          IF (.NOT.LAND_MASK(I).AND.ITILE.GT.1) THEN
             CD_T(I,ITILE)=CD_T(I,1)
             CH_T(I,ITILE)=CH_T(I,1)
             CD_STD_T(I,ITILE)=CD_STD_T(I,1)
          ENDIF
        ENDDO  ! loop over P-points

      ENDDO ! loop over tiles


      DO I=P1,P1+P_POINTS-1
!       IF ( an ordinary sea points (no sea-ice) or a land point)
        IF (.NOT.(ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) ) THEN
          CD_MIZ(I) = 1.E30
          CH_MIZ(I) = 1.E30
          CD_LEAD(I) = 1.E30
          CH_LEAD(I) = 1.E30
          RIB_LEAD(I) = 1.E30
        ENDIF
      ENDDO


!-----------------------------------------------------------------------
!!  4.  Loop round gridpoints to be processed, performing calculations
!!      AFTER call to FCDCH which necessitates splitting of loop.
!-----------------------------------------------------------------------

      DO ITILE=1,N_TYPES

!-----------------------------------------------------------------------
! 4.1 If the interactive surface resistance is requested call SF_STOM
!-----------------------------------------------------------------------

      IF (INT_STOM) THEN

!-----------------------------------------------------------------------
! Calculate the aerodynamic resistance
!-----------------------------------------------------------------------
        DO I=P1,P1+P_POINTS-1
          RA(I) = 1.0 / CH_T(I,ITILE)
        ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
!-----------------------------------------------------------------------
! For mesoscale model release assume uniform functional types and top
! leaf nitrogen concentrations. Assume that (fine) root biomass is
! equal to leaf biomass.
!-----------------------------------------------------------------------
          NL0(L) = 50.0E-3
          ROOT(L) = 0.04 * LAI(L,ITILE)

        ENDDO ! Loop over land-points


        IF(LAND_PTS.GT.0) THEN    ! Omit if no land points
          CALL SF_STOM  (
     &     LAND_PTS,LAND_FIELD,LAND_MASK,P1,LAND1,
     &     LAND_INDEX,
     &     P_POINTS,P_FIELD,
     &     F_TYPE(1,ITILE),CO2,HT(1,ITILE),PAR,LAI(1,ITILE),
     &     NL0,PSTAR,Q_1,RA,ROOT,TSTAR_TILE(1,ITILE),SMVCCL,
     &     V_ROOT(1,ITILE),SMVCWT,VFRAC(1,ITILE),GPP(1,ITILE),
     &     NPP(1,ITILE),RESP_P(1,ITILE),
     &     GC(1,ITILE),LTIMER,FSMC(1,ITILE))
        ENDIF                     ! End test on land points


!-----------------------------------------------------------------------
! Initialise gridbox mean carbon fluxes on uncalculated points
!-----------------------------------------------------------------------
      IF(LAND_FIELD.GT.0) THEN
        DO L=1,LAND1-1
          GPP(L,ITILE)=0.
          NPP(L,ITILE)=0.
          RESP_P(L,ITILE)=0.
        ENDDO
        DO L=LAND_PTS+LAND1,LAND_FIELD
          GPP(L,ITILE)=0.
          NPP(L,ITILE)=0.
          RESP_P(L,ITILE)=0.
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Convert carbon fluxes to gridbox mean values
!-----------------------------------------------------------------------


        DO L = LAND1,LAND1+LAND_PTS-1

            GPP(L,ITILE) = VFRAC(L,ITILE) * GPP(L,ITILE)
            NPP(L,ITILE) = VFRAC(L,ITILE) * NPP(L,ITILE)
            RESP_P(L,ITILE) = VFRAC(L,ITILE) * RESP_P(L,ITILE)

        ENDDO ! Loop over land-points

      ENDIF  ! INT_STOM


!-----------------------------------------------------------------------
!!  4.2 Recalculate RESFS using "true" CH and EPDT

!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
          EPDT(I) = -PSTAR(I)/(R*TSTAR_TILE(I,ITILE))*CH_T(I,ITILE)*
     &                         DQ(I,ITILE)*TIMESTEP

      ENDDO ! Loop over land-points


        CALL SF_RESIST (
     &   P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,LAND_MASK,INT_STOM,
     &   P1,LAND1,
     &   LAND_INDEX,
     &   ROOTD(1,ITILE),SMVCCL,SMVCWT,SMC(1,ITILE),V_SOIL,
     &   VFRAC(1,ITILE),CANOPY,CATCH(1,ITILE),DQ(1,ITILE),EPDT,
     &   LYING_SNOW,GC(1,ITILE),RESIST(1,ITILE),CH_T(1,ITILE),
     &   PSIS(1,ITILE),FRACA(1,ITILE),RESFS(1,ITILE),F_SE(1,ITILE),
     &   RESFT(1,ITILE),LTIMER
     &   )

      ENDDO ! loop over tiles


!-----------------------------------------------------------------------
!!  4.D Call SFL_INT to calculate CDR10M, CHR1P5M and CER1P5M -
!!      interpolation coefficients used in SF_EVAP and IMPL_CAL to
!!      calculate screen temperature, specific humidity and 10m winds.
!-----------------------------------------------------------------------

      IF (SU10 .OR. SV10 .OR. SQ1P5 .OR. ST1P5) THEN

!sjtemp        ITILE=3 ! short grass tile

         ITILE=1  ! single tile mode only

        CALL SFL_INT (
     &  P_POINTS,P_FIELD,P1,
     &  Z0M_EFF_T(1,ITILE),Z0H_T(1,ITILE),CD_T(1,ITILE),CH_T(1,ITILE),
     &  Z0M_T(1,ITILE),CD_STD_T(1,ITILE),
     &  RESFT(1,ITILE),RECIP_L_MO(1,ITILE),
     &  V_S(1,ITILE),V_S_STD(1,ITILE),
     &  CDR10M,CHR1P5M,CER1P5M,
     &  SU10,SV10,ST1P5,SQ1P5,
     &  LTIMER
     & )

      ENDIF

!-----------------------------------------------------------------------
!!  4.2 Now that diagnostic calculations are over, update sea ice CD
!!      and CH to their correct values (i.e. gridsquare means).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I).LT. 0.7 ) THEN
            CD_T(I,1) = ( ICE_FRACT(I)*CD_MIZ(I) +
     &                (0.7-ICE_FRACT(I))*CD_LEAD(I) ) / 0.7  ! P2430.5
            CH_T(I,1) = ( ICE_FRACT(I)*CH_MIZ(I) +
     &                (0.7-ICE_FRACT(I))*CH_LEAD(I) ) / 0.7  ! P2430.4
            CD_STD_T(I,1)=CD_T(I,1)  ! for SCYCLE: no orog. over sea+ice
          ELSE
            CD_T(I,1) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +
     &                (ICE_FRACT(I)-0.7)*CD_T(I,1) ) / 0.3     ! P2430.7
            CH_T(I,1) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +
     &              (ICE_FRACT(I)-0.7)*CH_T(I,1) ) / 0.3       ! P2430.7
            CD_STD_T(I,1)=CD_T(I,1)  ! for SCYCLE: no orog. over sea+ice
          ENDIF
        ENDIF

      ENDDO !loop over points for sea ice


      DO ITILE=1,N_TYPES
        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!!  4.3 Calculate the surface exchange coefficients RHOK(*).
!-----------------------------------------------------------------------

          RHOSTAR(I,ITILE) = PSTAR(I) / ( R*TSTAR_TILE(I,ITILE) )
!                        ... surface air density from ideal gas equation

          RHOKM_1(I,ITILE) = RHOSTAR(I,ITILE) * CD_T(I,ITILE)
                                                            ! P243.124
          RHOKH_1(I,ITILE) = RHOSTAR(I,ITILE) * CH_T(I,ITILE)
                                                            ! P243.125
          RHOKE(I,ITILE) = RESFT(I,ITILE) * RHOKH_1(I,ITILE)

!  Calculate resistances for use in Sulphur Cycle
!  (Note that CD_STD, CH and VSHR should never = 0)
          RHO_ARESIST(I) = RHO_ARESIST(I) + TILE_FRAC(I,ITILE) *
     &                 (RHOSTAR(I,ITILE) * CD_STD_T(I,ITILE))

          ARESIST(I) = ARESIST(I) + TILE_FRAC(I,ITILE) /
     &                              CD_STD_T(I,ITILE)

          RESIST_B(I)= RESIST_B(I) + TILE_FRAC(I,ITILE)*
     &                (CD_STD_T(I,ITILE)/CH_T(I,ITILE) - 1.0) /
     &                 CD_STD_T(I,ITILE)

!     RHOSTAR * CD * VSHR stored for diagnostic output before
!     horizontal interpolation.

        ENDDO ! loop over p-points
      ENDDO ! n_types


      DO ITILE=1,N_TYPES
        IF(ITILE.EQ.1) THEN
          L_LAND=.FALSE.
        ELSE
          L_LAND=.TRUE.
        ENDIF

        CALL SF_FLUX (
     &   P_POINTS,P_FIELD,LAND_PTS,LAND_FIELD,LAND_MASK,L_LAND,P1,LAND1,
     &   LAND_INDEX,
     &   ALPHA1(1,ITILE),DQ(1,ITILE),DQ_LEAD,DTEMP(1,ITILE),DTEMP_LEAD,
     &   DZSOIL,HCONS,ICE_FRACT,
     &   LYING_SNOW,QS_BLEND,QW_BLEND,RADNET_C(1,ITILE),RESFT(1,ITILE),
     &   RHOKE(1,ITILE),RHOKH_1(1,ITILE),TI,TL_BLEND,TS1,
     &   Z0H_T(1,ITILE),Z0M_EFF_T(1,ITILE),Z1_TQ,Z1_UV,
     &   ASHTF,E_SEA,EPOT(1,ITILE),FQW_1(1,ITILE),FTL_1(1,ITILE),H_SEA,
     &   RHOKPM(1,ITILE),RHOKPM_POT(1,ITILE),LTIMER 
     &,  TSTAR_TILE(1,ITILE),VFRAC(1,ITILE),TIMESTEP,CANCAP(1,ITILE)
     &   )

      ENDDO ! n_types

      DO ITILE=1,N_TYPES
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
! average fluxes, resistances and other things

            FTL1_GB(I)=FTL1_GB(I)+FTL_1(I,ITILE)*TILE_FRAC(I,ITILE)
            FQW1_GB(I)=FQW1_GB(I)+FQW_1(I,ITILE)*TILE_FRAC(I,ITILE)
            EPOT_GB(I)=EPOT_GB(I)+EPOT(I,ITILE)*TILE_FRAC(I,ITILE)

            RESFS_GB(I) = RESFS_GB(I) +
     &                    TILE_FRAC(I,ITILE) * RESFS(I,ITILE)
            RESFT_GB(I) = RESFT_GB(I) +
     &                    TILE_FRAC(I,ITILE) * RESFT(I,ITILE)

            RHOKH_1_GB(I) = RHOKH_1_GB(I) +
     &                      RHOKH_1(I,ITILE) * TILE_FRAC(I,ITILE)
            RHOKM_1_GB(I) = RHOKM_1_GB(I) +
     &                      RHOKM_1(I,ITILE) * TILE_FRAC(I,ITILE)
            RHOKE_GB(I) = RHOKE_GB(I) +
     &                    RHOKE(I,ITILE) * TILE_FRAC(I,ITILE)
            RHOKPM_GB(I) = RHOKPM_GB(I) +
     &                     RHOKPM(I,ITILE) * TILE_FRAC(I,ITILE)
            RHOKPM_POT_GB(I) = RHOKPM_POT_GB(I) +
     &                     RHOKPM_POT(I,ITILE) * TILE_FRAC(I,ITILE)
            FSMC_GB(L) = FSMC_GB(L) +
     &                     FSMC(L,ITILE) * TILE_FRAC(I,ITILE)

            ALPHA1_GB(I) = ALPHA1_GB(I) +
     &                     ALPHA1(I,ITILE) * TILE_FRAC(I,ITILE)

            CD(I) = CD(I) + CD_T(I,ITILE) * TILE_FRAC(I,ITILE)

            CD_STD(I) = CD_STD(I) +
     &                   CD_STD_T(I,ITILE) * TILE_FRAC(I,ITILE)

            CH(I) = CH(I) + CH_T(I,ITILE) * TILE_FRAC(I,ITILE)

            BT1_GB(I) = BT1_GB(I) + BT_1(I,ITILE) * TILE_FRAC(I,ITILE)
            BQ1_GB(I) = BQ1_GB(I) + BQ_1(I,ITILE) * TILE_FRAC(I,ITILE)
            RHOSTAR_GB(I) = PSTAR(I) / ( R*TSTAR_GB(I) )
!                        ... surface air density from ideal gas equation


        ENDDO ! Loop over land-points

      ENDDO ! loop over tiles


      DO I=P1,P1+P_POINTS-1
        IF(.NOT.LAND_MASK(I)) THEN
          FTL1_GB(I) = FTL_1(I,1)
          FQW1_GB(I) = FQW_1(I,1)
          EPOT_GB(I) = EPOT(I,1)

          RESFS_GB(I) = RESFS(I,1)
          RESFT_GB(I) = RESFT(I,1)

          RHOKH_1_GB(I) = RHOKH_1(I,1)
          RHOKM_1_GB(I) = RHOKM_1(I,1)
          RHOKE_GB(I) = RHOKE(I,1)
          RHOKPM_GB(I) = RHOKPM(I,1)
          RHOKPM_POT_GB(I) = RHOKPM_POT(I,1)

          ALPHA1_GB(I) = ALPHA1(I,1)

          CD(I) = CD_T(I,1)
          CD_STD(I) = CD_STD_T(I,1)
          CH(I) = CH_T(I,1)

          BT1_GB(I) = BT_1(I,1)
          BQ1_GB(I) = BQ_1(I,1)
          RHOSTAR_GB(I) = RHOSTAR(I,1)

        ENDIF

        RHO_CD_MODV1(I) = RHOKM_1_GB(I) ! diagnostic required for VAR

      ENDDO

!-----------------------------------------------------------------------
!!  4.4   Calculate the standard deviations of layer 1 turbulent
!!        fluctuations of temperature and humidity using approximate
!!        formulae from first order closure.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1

        U_S(I) = SQRT(CD(I) * VSHR(I))
        FB_SURF(I) = G * ( BT1_GB(I)*FTL1_GB(I) +
     &                     BQ1_GB(I)*FQW1_GB(I) ) / RHOSTAR_GB(I)

        W_S_CUBED = 75.0 * FB_SURF(I)
C       ! 75.0 = 2.5 * height above the surface of 30 m
C       !---------------------------------------------------------------
C       ! Only calculate standard deviations for unstable surface layers
C       !---------------------------------------------------------------
        IF (W_S_CUBED .GT. 0.0) THEN
          W_M  = ( W_S_CUBED + U_S(I) * U_S(I) * U_S(I) ) ** (1.0/3.0)
          T1_SD(I) = 1.93 * FTL1_GB(I) / (RHOSTAR_GB(I) * W_M)
          Q1_SD(I) = 1.93 * FQW1_GB(I) / (RHOSTAR_GB(I) * W_M)
          TV1_SD(I) = T_1(I) *
     &                ( 1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I) ) *
     &                ( BT1_GB(I)*T1_SD(I) + BQ1_GB(I)*Q1_SD(I) )
          T1_SD(I) = MAX ( 0.0 , T1_SD(I) )
          Q1_SD(I) = MAX ( 0.0 , Q1_SD(I) )
          IF (TV1_SD(I) .LE. 0.0) THEN
            TV1_SD(I) = 0.0
            T1_SD(I) = 0.0
            Q1_SD(I) = 0.0
          ENDIF
        ELSE
          T1_SD(I) = 0.0
          Q1_SD(I) = 0.0
          TV1_SD(I) = 0.0
        ENDIF
!-----------------------------------------------------------------------
!!  4.5 For diagnostic output calculate the dimensionless surface
!!      transfer coefficients.
!----------------------------------------------------------------------
        CD(I) = CD(I) / VSHR(I)
        CH(I) = CH(I) / VSHR(I)
!
      ENDDO

!-----------------------------------------------------------------------
!!  4.6 For sea points, calculate the wind mixing energy flux and the
!!      sea-surface roughness length on the P-grid, using time-level n
!!      quantities.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1

        IF (.NOT.LAND_MASK(I)) THEN
          TAU = RHOSTAR_GB(I) * V_S(I,1) * V_S(I,1)
          IF (ICE_FRACT(I) .GT. 0.0)
     &      TAU = RHOSTAR_GB(I) * V_S_LEAD(I) * V_S_LEAD(I)
          IF (SFME) FME(I) = (1.0-ICE_FRACT(I)) * TAU * SQRT(TAU/RHOSEA)
!                                                             ! P243.96
          Z0MSEA(I) = 1.54E-6 / SQRT(TAU / RHOSTAR_GB(I)) +
     &                (CHARNOCK/G) * (TAU / RHOSTAR_GB(I))
!                                                  ... (S.Smith formula)
        ENDIF ! of IF (.NOT. LAND_MASK), land-points done in next loop.
      ENDDO ! Loop over points for sections 4.2 - 4.6
      DO L=LAND1,LAND1+LAND_PTS-1
      I = LAND_INDEX(L)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  4.7 Set Z0MSEA to Z0V, FME to zero for land points.
!   (Former because UM uses same storage for Z0V
!   and Z0MSEA.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Z0MSEA(I) = Z0V_GB(I)

      IF (SFME) FME(I) = 0.0

      ENDDO ! Loop over points for section 4.7

      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',4)
      ENDIF

      RETURN
      END
