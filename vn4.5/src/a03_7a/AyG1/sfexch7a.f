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
!!!  Model            Modification history:
!!! version  Date
!!!  4.3   17/11/95   New deck      Simon Jackson
!!!  4.4    16/7/97   Version for MOSES II tile model.  Richard Essery
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!  4.5    17/11/98  Introduce Z0H_Z0M and initialise FTL_TILE and
!!!                   RIB_TILE on all tiles at all points. Richard Betts
!!!
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version 2, dated 18/1/90.
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!                 See especially sub-section (ix).
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE SF_EXCH (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,LAND_INDEX,
     & TILE_INDEX,TILE_PTS,
     & BQ_1,BT_1,CANOPY,CATCH,DZSOIL,GC,HCONS,HO2R2_OROG,
     & ICE_FRACT,LYING_SNOW,PSTAR,P_1,QW_1,RADNET,RADNET_SNOW,SIL_OROG,
     & SMVCST,TILE_FRAC,TIMESTEP,TL_1,TI,TS1,TSNOW,TSTAR_TILE,TSTAR,
     & VSHR,Z0_TILE,Z0_SF_GB,Z1_UV,Z1_TQ,LAND_MASK,
     & SU10,SV10,SQ1P5,ST1P5,SFME,LTIMER,L_Z0_OROG,Z0MSEA,
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_SNOW,CD,CH,CDR10M,
     & CHR1P5M,CHR1P5M_SICE,E_SEA,FME,FQW_1,FQW_TILE,FQW_ICE,
     & FTL_1,FTL_TILE,FTL_ICE,FRACA,H_BLEND_OROG,H_SEA,
     & Q1_SD,RESFS,RESFT,RIB,RIB_TILE,T1_SD,Z0M_EFF,
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,RHO_ARESIST,ARESIST,RESIST_B,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & RHO_CD_MODV1,RHOKH_1,RHOKH_1_SICE,RHOKM_1,RHOKPM,RHOKPM_SICE,
     & NRML
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS              ! IN Number of P-grid points to be
!                            !    processed.
     &,P_FIELD               ! IN Total number of P-grid points.
     &,P1                    ! IN First P-point to be processed.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.
     &,LAND_FIELD            ! IN Total number of land points.
     &,NTYPE                 ! IN Number of tiles per land point.
     &,LAND_INDEX(P_FIELD)   ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD,NTYPE)
!                            ! IN Index of tile points.
     &,TILE_PTS(NTYPE)       ! IN Number of tile points.

      REAL
     & BQ_1(P_FIELD)         ! IN A buoyancy parameter for lowest atm
!                            !    level ("beta-q twiddle").
     &,BT_1(P_FIELD)         ! IN A buoyancy parameter for lowest atm
!                            !    level ("beta-T twiddle").
     &,CANOPY(LAND_FIELD,NTYPE-1)
!                            ! IN Surface water for land tiles
!                            !    (kg/m2).
     &,CATCH(LAND_FIELD,NTYPE-1)
!                            ! IN Surface capacity (max. surface water)
!                            !    of snow-free land tiles (kg/m2).
     &,DZSOIL                ! IN Soil or land-ice surface layer
!                            !    thickness (m).
     &,GC(LAND_FIELD,NTYPE)  ! IN "Stomatal" conductance to evaporation
!                            !    for land tiles (m/s).
     &,HCONS(LAND_FIELD)     ! IN Soil thermal conductivity including
!                            !    effects of water and ice (W/m/K).
     &,HO2R2_OROG(LAND_FIELD)! IN Peak to trough height of unresolved
!                            !    orography divided by 2SQRT(2) (m).
     &,ICE_FRACT(P_FIELD)    ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_FIELD)   ! IN Lying snow amount (kg per sq metre).
     &,PSTAR(P_FIELD)        ! IN Surface pressure (Pascals).
     &,P_1(P_FIELD)          ! IN Level 1 atmospheric pressure.
     &,QW_1(P_FIELD)         ! IN Total water content of lowest
!                            !    atmospheric layer (kg per kg air).
     &,RADNET(P_FIELD)       ! IN Net surface radiation over snow-free
!                            !    land or sea-ice (W/m2)
     &,RADNET_SNOW(P_FIELD)  ! IN Net surface radiation over snow or
!                            !    land-ice (W/m2)
     &,SIL_OROG(LAND_FIELD)  ! IN Silhouette area of unresolved
!                            !    orography per unit horizontal area
     &,SMVCST(LAND_FIELD)    ! IN Volumetric saturation point
!                            !    - zero at land-ice points.
     &,TILE_FRAC(LAND_FIELD,NTYPE)
!                            ! IN Tile fractions.
     &,TIMESTEP              ! IN Timestep in seconds for EPDT calc.
     &,TL_1(P_FIELD)         ! IN Liquid/frozen water temperature for
!                            !    lowest atmospheric layer (K).
     &,TI(P_FIELD)           ! IN Temperature of sea-ice surface layer
!                            !    (K)
     &,TS1(LAND_FIELD)       ! IN Temperature of top soil or land-ice
!                            !    layer (K)
     &,TSNOW(LAND_FIELD)     ! IN Temperature of surface snow layer (K)
!                            !    = TS1 at land-ice points.
     &,TSTAR_TILE(LAND_FIELD,NTYPE)
!                            ! IN Tile surface temperatures (K).
     &,TSTAR(P_FIELD)        ! IN Gridbox mean surface temperature (K).
     &,VSHR(P_FIELD)         ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear
     &,Z0_TILE(LAND_FIELD,NTYPE)
!                            ! IN Tile roughness lengths (m).
     &,Z0_SF_GB(P_FIELD)     ! IN Snow-free GBM roughness length (m).
     &,Z1_UV(P_FIELD)        ! IN Height of lowest uv level (m).
     &,Z1_TQ(P_FIELD)        ! IN Height of lowest tq level (m).
!                            !    Note, if the grid used is staggered in
!                            !    the vertical, Z1_UV and Z1_TQ can be
!                            !    different.

      LOGICAL
     & LAND_MASK(P_FIELD)    ! IN .TRUE. for land; .FALSE. elsewhere.
     &,SU10                  ! IN STASH flag for 10-metre W wind.
     &,SV10                  ! IN STASH flag for 10-metre S wind.
     &,SQ1P5                 ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                 ! IN STASH flag for 1.5-metre temperature.
     &,SFME                  ! IN STASH flag for wind mixing energy flux
     &,LTIMER                ! IN Logical for TIMER.
     &,L_Z0_OROG             ! IN .TRUE. to use orographic roughness.

!  Modified (INOUT) variables.

      REAL
     & Z0MSEA(P_FIELD)       ! INOUT Sea-surface roughness length for
!                            !       momentum (m).  F617.

!  Output variables.
!
      REAL
     & ALPHA1(LAND_FIELD,NTYPE)
!                            ! OUT Gradients of saturated specific
!                            !     humidity with respect to temperature
!                            !     between the bottom model layer and
!                            !     tile surface
     &,ALPHA1_SICE(P_FIELD)  ! OUT ALPHA1 for sea-ice.
     &,ASHTF(P_FIELD)        ! OUT Coefficient to calculate surface heat
!                            !     flux into soil or sea-ice (W/m2/K)
     &,ASHTF_SNOW(P_FIELD)   ! OUT Coefficient to calculate surface heat
!                            !     flux into snow (W/m2/K)
     &,CD(P_FIELD)           ! OUT Bulk transfer coefficient for
!                            !      momentum.
     &,CH(P_FIELD)           ! OUT Bulk transfer coefficient for heat
!                            !     and/or moisture.
     &,CDR10M(P_FIELD)       ! OUT Reqd for calculation of 10m wind
!                            !     (u & v).
!                            !     NBB: This is output on the UV-grid,
!                            !     but with the first and last rows set
!                            !     to a "missing data indicator".
!                            !     Sea-ice leads ignored.
     &,CHR1P5M(LAND_FIELD,NTYPE)
!                            ! OUT Reqd for calculation of 1.5m temp for
!                            !     land tiles.
     &,CHR1P5M_SICE(P_FIELD) ! OUT CHR1P5M for sea and sea-ice
!                            !     (leads ignored).
     &,E_SEA(P_FIELD)        ! OUT Evaporation from sea times leads
!                            !     fraction (kg/m2/s). Zero over land.
     &,FME(P_FIELD)          ! OUT Wind mixing energy flux (Watts/sq m).
     &,FQW_1(P_FIELD)        ! OUT "Explicit" surface flux of QW (i.e.
!                            !     evaporation), on P-grid (kg/m2/s).
!                            !     for whole grid-box
     &,FQW_TILE(LAND_FIELD,NTYPE)
!                            ! OUT Local FQW_1 for land tiles.
     &,FQW_ICE(P_FIELD)      ! OUT GBM FQW_1 for sea-ice.
     &,FTL_1(P_FIELD)        ! OUT "Explicit" surface flux of TL = H/CP.
!                            !     (sensible heat / CP). grid-box mean
     &,FTL_TILE(LAND_FIELD,NTYPE)
!                            ! OUT Local FTL_1 for land tiles.
     &,FTL_ICE(P_FIELD)      ! OUT GBM FTL_1 for sea-ice.
     &,FRACA(LAND_FIELD,NTYPE-1)
!                            ! OUT Fraction of surface moisture flux
!                            !     with only aerodynamic resistance
!                            !     for snow-free land tiles.
     &,H_BLEND_OROG(P_FIELD) ! OUT Blending height for orographic
!                            !     roughness
     &,H_SEA(P_FIELD)        ! OUT Surface sensible heat flux over sea
!                            !     times leads fraction (W/m2).
!                            !     Zero over land.
     &,Q1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     specific humidity (kg/kg).
     &,RESFS(LAND_FIELD,NTYPE-1)
!                            ! OUT Combined soil, stomatal and
!                            !     aerodynamic resistance factor for
!                            !     fraction 1-FRACA of snow-free tiles
     &,RESFT(LAND_FIELD,NTYPE)
!                            ! OUT Total resistance factor
!                            !     FRACA+(1-FRACA)*RESFS for snow-free
!                            !     tiles, 1 for snow.
     &,RIB(P_FIELD)          ! OUT Mean bulk Richardson number for
!                            !     lowest layer
     &,RIB_TILE(LAND_FIELD,NTYPE)
!                            ! OUT RIB for land tiles.
     &,T1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     temperature (K).
     &,Z0M_EFF(P_FIELD)      ! OUT Effective roughness length for
!                            !     momentum
     &,Z0H(P_FIELD)          ! OUT Roughness length for heat
!                            !     and moisture
     &,Z0H_TILE(LAND_FIELD,NTYPE)
!                            ! OUT Tile roughness lengths for heat
!                            !     and moisture
     &,Z0M(P_FIELD)          ! OUT Roughness length for momentum
     &,Z0M_TILE(LAND_FIELD,NTYPE)
!                            ! OUT Tile roughness lengths for momentum
     &,RHO_ARESIST(P_FIELD)  ! OUT RHOSTAR*CD_STD*VSHR  for SCYCLE
     &,ARESIST(P_FIELD)      ! OUT 1/(CD_STD*VSHR)      for SCYCLE
     &,RESIST_B(P_FIELD)     ! OUT (1/CH-1/CD_STD)/VSHR for SCYCLE
     &,RHO_ARESIST_TILE(LAND_FIELD,NTYPE)
!                            ! OUT RHOSTAR*CD_STD*VSHR on land tiles
     &,ARESIST_TILE(LAND_FIELD,NTYPE)
!                            ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_FIELD,NTYPE)
!                            ! OUT (1/CH-1/CD_STD)/VSHR on land tiles

! Surface exchange coefficients;passed to subroutine IMPL_CAL
      REAL
     & RHO_CD_MODV1(P_FIELD) ! OUT rhostar*cD*vshr before horizontal
!                            !     interpolation output as a diagnostic.
     &,RHOKH_1(LAND_FIELD,NTYPE)
!                            ! OUT Surface exchange coefficient for land
!                            !     tiles.
     &,RHOKH_1_SICE(P_FIELD) ! OUT Surface exchange coefficient for sea
!                            !     or sea-ice.
     &,RHOKM_1(P_FIELD)      ! OUT For momentum. NB: This is output on
!                            !     UV-grid, but with the first and last
!                            !     rows set to "missing data indicator".
     &,RHOKPM(LAND_FIELD,NTYPE)
!                            ! OUT Mixing coefficient for land tiles.
     &,RHOKPM_SICE(P_FIELD)  ! OUT Mixing coefficient for sea-ice.

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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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


! Derived local parameters.
      REAL LS
      PARAMETER (
     & LS=LF+LC            ! Latent heat of sublimation.
     & )

!   (2) Boundary Layer local parameters.

      REAL
     + LB                         ! Blending height (m).
      PARAMETER (LB = 550.0)
C*L-----------COMDECK C_CHARNK FOR SUBROUTINE SF_EXCH----------
C CHARNOCK is a constant in the Charnock formula for sea-surface
C          roughness length for momentum (Z0MSEA).
      REAL CHARNOCK

      PARAMETER(CHARNOCK = 0.012)
C*----------------------------------------------------------------------
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
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
C*L-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
C Z0FSEA = roughness length for free convective heat and moisture
C          transport over the sea (m).
C Z0HSEA = roughness length for free heat and moisture transport
C          over the sea (m).
C Z0MIZ  = roughness length for heat, moisture and momentum over
C          the Marginal Ice Zone (m).
C Z0SICE = roughness length for heat, moisture and momentum over
C          sea-ice (m).
      REAL Z0FSEA,Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0FSEA = 1.3E-3,
     &          Z0HSEA = 1.0E-4,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
C*----------------------------------------------------------------------
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

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
      REAL
     + Z0H_Z0M(9)                 ! Ratio of roughness length for heat
C                                 ! to roughness length for momentum.
C----------------------------------------------------------------------
C                         BT   NT   C3G  C4G  Shr  Urb  Wat  Soil Ice
C----------------------------------------------------------------------
      DATA Z0H_Z0M     /  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 /

      REAL H_BLEND_MIN
      PARAMETER (
     & H_BLEND_MIN=0.0       ! Minimum blending height.
     &)

!   External subprograms called.

      EXTERNAL SF_OROG,SF_OROG_GB,QSAT,SFL_INT,SF_RESIST,TIMER,
     &         STDEV1_SEA,STDEV1_LAND,SF_RIB_SEA,SF_RIB_LAND,
     &         FCDCH_SEA,FCDCH_LAND,SF_FLUX_SEA,SF_FLUX_LAND

!   Define local storage.

!   (a) Workspace.

      REAL
     & QS1(P_FIELD)                ! Sat. specific humidity
!                                  ! qsat(TL_1,PSTAR)
     &,RHOSTAR(P_FIELD)            ! Surface air density

!  Workspace for sea and sea-ice leads
      REAL
     & CD_SEA(P_FIELD)             ! Drag coefficient
     &,CH_SEA(P_FIELD)             ! Transfer coefficient for heat and
!                                  ! moisture
     &,QSTAR_SEA(P_FIELD)          ! Surface saturated sp humidity
     &,RIB_SEA(P_FIELD)            ! Bulk Richardson number
     &,TSTAR_SEA(P_FIELD)          ! Surface temperature
     &,Z0F_SEA(P_FIELD)            ! Roughness length for free-convec.
!                                  ! heat and moisture transport
     &,Z0H_SEA(P_FIELD)            ! Roughness length for heat and
!                                  ! moisture transport

!  Workspace for sea-ice and marginal ice zone
      REAL
     & CD_ICE(P_FIELD)             ! Drag coefficient
     &,CD_MIZ(P_FIELD)             ! Drag coefficient
     &,CH_ICE(P_FIELD)             ! Transfer coefficient for heat and
!                                  ! moisture
     &,CH_MIZ(P_FIELD)             ! Transfer coefficient for heat and
!                                  ! moisture
     &,QSTAR_ICE(P_FIELD)          ! Surface saturated sp humidity
     &,RIB_ICE(P_FIELD)            ! Bulk Richardson number
     &,RIB_MIZ(P_FIELD)            ! Bulk Richardson number
     &,TSTAR_ICE(P_FIELD)          ! Surface temperature
     &,Z0_ICE(P_FIELD)             ! Roughness length.
     &,Z0_MIZ(P_FIELD)             ! Roughness length.
      INTEGER
     & SICE_INDEX(P_FIELD)         ! Index of sea-ice points
     &,NSICE                       ! Number of sea-ice points.

!  Workspace for land tiles
      REAL
     & CD_STD(LAND_FIELD,NTYPE)    ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
     &,CD_TILE(LAND_FIELD,NTYPE)   ! Drag coefficient
     &,CH_TILE(LAND_FIELD,NTYPE)   ! Transfer coefficient for heat and
!                                  ! moisture
     &,CHN(LAND_FIELD)             ! Neutral value of CH.
     &,DQ(LAND_FIELD)              ! Sp humidity difference between
!                                  ! surface and lowest atmospheric lev
     &,EPDT(LAND_FIELD)            ! "Potential" Evaporation * Timestep
     &,PSTAR_LAND(LAND_FIELD)      ! Surface pressure for land points.
     &,QSTAR_TILE(LAND_FIELD,NTYPE)! Surface saturated sp humidity.
     &,RHOKM_1_TILE(LAND_FIELD,NTYPE)
!                                  ! Momentum exchange coefficient.
     &,WIND_PROFILE_FACTOR(LAND_FIELD,NTYPE)
!                                  ! For transforming effective surface
!                                  ! transfer coefficients to those
!                                  ! excluding form drag.
     &,Z0_GB(LAND_FIELD)           ! GBM roughness length including snow
     &,Z0M_EFF_TILE(LAND_FIELD,NTYPE)
!                                  ! Effective momentum roughness length
     &,Z0F_TILE(LAND_FIELD,NTYPE)  !Roughness length for free convective
!                                  ! heat and moisture transport

!   (b) Scalars.

      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Loop counter (tile field index).
     &,L           ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).
      REAL
     & TAU         ! Magnitude of surface wind stress over sea.
     &,ZETAM       ! Temporary in calculation of CHN.
     &,ZETAH       ! Temporary in calculation of CHN.
     &,ZETA1       ! Work space
     &,Z0          ! yet more workspace

      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  0. Initialise FTL_TILE and RIB_TILE on all tiles at all points,
!!     to allow STASH to process these as diagnostics.
!-----------------------------------------------------------------------
      DO N=1,NTYPE
        DO L=1,LAND_FIELD
          FTL_TILE(L,N) = 0.0
          RIB_TILE(L,N) = 0.0
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!  1. Index array for sea-ice
!-----------------------------------------------------------------------

      NSICE = 0
      DO I=P1,P1+P_POINTS-1
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          NSICE = NSICE + 1
          SICE_INDEX(NSICE) = I
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
!!  2.  Calculate QSAT values required later.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        TSTAR_SEA(I) = TSTAR(I)
        TSTAR_ICE(I) = TSTAR(I)
        RHOSTAR(I) = PSTAR(I) / ( R*TSTAR(I) )
!                        ... surface air density from ideal gas equation
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          TSTAR_ICE(I) = ( TSTAR(I) - (1.0-ICE_FRACT(I))*TFS )
     &                    / ICE_FRACT(I)                       ! P2430.1
          TSTAR_SEA(I) = TFS
        ENDIF
      ENDDO
      CALL QSAT(QS1(P1),TL_1(P1),PSTAR(P1),P_POINTS)
      CALL QSAT(QSTAR_SEA(P1),TSTAR_SEA(P1),PSTAR(P1),P_POINTS)
      CALL QSAT(QSTAR_ICE(P1),TSTAR_ICE(P1),PSTAR(P1),P_POINTS)
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        PSTAR_LAND(L) = PSTAR(I)
      ENDDO
      DO N=1,NTYPE
        CALL QSAT(QSTAR_TILE(LAND1,N),TSTAR_TILE(LAND1,N),
     &            PSTAR_LAND(LAND1),LAND_PTS)
      ENDDO

!-----------------------------------------------------------------------
!!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!  3.1 Calculate neutral roughness lengths
!-----------------------------------------------------------------------

! Sea, sea-ice leads, sea-ice and marginal ice zone
      DO I=P1,P1+P_POINTS-1
        Z0H_SEA(I) = Z0HSEA
        Z0F_SEA(I) = Z0FSEA
        Z0_MIZ(I) = Z0MIZ
        Z0_ICE(I) = Z0SICE
        RIB_SEA(I) = 0.
        RIB_ICE(I) = 0.
      ENDDO

! Land tiles
! Z0_TILE contains the appropriate value for land-ice points, but has to
! be modified for snow-cover on non-land-ice points
      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          Z0M_TILE(L,N) = Z0_TILE(L,N)
          IF ( N.EQ.NTYPE .AND. SMVCST(L).NE.0. ) THEN
            I = LAND_INDEX(L)
            Z0 = Z0_SF_GB(I) - 4.0E-4*LYING_SNOW(I)/TILE_FRAC(L,N)
            ZETA1 = MIN( 5.0E-4 , Z0_SF_GB(I) )
            Z0M_TILE(L,N) = MAX( ZETA1 , Z0 )
          ENDIF
          Z0H_TILE(L,N) = Z0H_Z0M(N)*Z0M_TILE(L,N)
          Z0F_TILE(L,N) = Z0H_Z0M(N)*Z0M_TILE(L,N)
          RIB_TILE(L,N) = 0.
        ENDDO
      ENDDO

      DO N=1,NTYPE
        CALL SF_OROG (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   L_Z0_OROG,LTIMER,
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG,Z0M_TILE(1,N),Z1_UV,
     &   WIND_PROFILE_FACTOR(1,N),Z0M_EFF_TILE(1,N)
     &   )
      ENDDO

!-----------------------------------------------------------------------
! Calculate RESFT with neutral CH and EPDT=0 for use in calculation
! of Richardson number. RESFT=1 for snow.
!-----------------------------------------------------------------------

! Snow-free land tiles
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          ZETAM = LOG ( (Z1_UV(I) + Z0M_TILE(L,N))/Z0M_TILE(L,N) )
          ZETAH = LOG ( (Z1_TQ(I) + Z0M_TILE(L,N))/Z0H_TILE(L,N) )
          CHN(L) = (VKMAN/ZETAH)*(VKMAN/ZETAM)*WIND_PROFILE_FACTOR(L,N)
          DQ(L) = QW_1(I) - QSTAR_TILE(L,N)
          EPDT(L) = 0.0
        ENDDO
        CALL SF_RESIST (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   CANOPY(1,N),CATCH(1,N),CHN,DQ,EPDT,GC(1,N),VSHR,
     &   FRACA(1,N),RESFS(1,N),RESFT(1,N),LTIMER
     &   )
      ENDDO

! Snow and land-ice tile
      DO J=1,TILE_PTS(NTYPE)
        L = TILE_INDEX(J,NTYPE)
        RESFT(L,NTYPE) = 1.
      ENDDO

!-----------------------------------------------------------------------
!!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------

! Sea, sea-ice and sea-ice leads
      CALL SF_RIB_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK,NSICE,SICE_INDEX,
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_ICE,
     & TSTAR_SEA,VSHR,Z0_ICE,Z0H_SEA,Z0_ICE,Z0MSEA,Z1_TQ,Z1_UV,
     & RIB_SEA,RIB_ICE,LTIMER
     & )

! Land tiles
      DO N=1,NTYPE
        CALL SF_RIB_LAND (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   BQ_1,BT_1,QSTAR_TILE(1,N),QW_1,RESFT(1,N),TL_1,
     &   TSTAR_TILE(1,N),VSHR,Z0H_TILE(1,N),Z0M_TILE(1,N),Z1_TQ,Z1_UV,
     &   RIB_TILE(1,N),LTIMER
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  3.3 Calculate stability corrected effective roughness length.
!!  Stability correction only applies to land points.
!-----------------------------------------------------------------------

      DO N=1,NTYPE
        CALL SF_OROG (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   L_Z0_OROG,LTIMER,
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG,Z0M_TILE(1,N),Z1_UV,
     &   WIND_PROFILE_FACTOR(1,N),Z0M_EFF_TILE(1,N)
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  3.4 Calculate CD, CH via routine FCDCH.
!-----------------------------------------------------------------------

! Sea-ice
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,LAND_MASK,
     &               RIB_ICE,Z0_ICE,Z0_ICE,Z0_ICE,Z1_UV,Z1_TQ,
     &               CD_ICE,CH_ICE,LTIMER)

! Marginal Ice Zone
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,LAND_MASK,
     &               RIB_ICE,Z0_MIZ,Z0_MIZ,Z0_MIZ,Z1_UV,Z1_TQ,
     &               CD_MIZ,CH_MIZ,LTIMER)

! Sea and sea-ice leads
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,LAND_MASK,
     &               RIB_SEA,Z0MSEA,Z0H_SEA,Z0F_SEA,Z1_UV,Z1_TQ,
     &               CD_SEA,CH_SEA,LTIMER)

! Land tiles
      DO N=1,NTYPE
        CALL FCDCH_LAND (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N),LAND_INDEX,
     &   RIB_TILE(1,N),WIND_PROFILE_FACTOR(1,N),
     &   Z0M_EFF_TILE(1,N),Z0H_TILE(1,N),Z0F_TILE(1,N),Z1_UV,Z1_TQ,
     &   CD_TILE(1,N),CH_TILE(1,N),CD_STD(1,N),LTIMER
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  4.1 Recalculate RESFT using "true" CH and EPDT for snow-free land
!!      tiles
!-----------------------------------------------------------------------

      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          DQ(L) = QW_1(I) - QSTAR_TILE(L,N)
          EPDT(L) = - RHOSTAR(I)*CH_TILE(L,N)*VSHR(I)*DQ(L)*TIMESTEP
        ENDDO
        CALL SF_RESIST (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   CANOPY(1,N),CATCH(1,N),CH_TILE(1,N),DQ,EPDT,GC(1,N),VSHR,
     &   FRACA(1,N),RESFS(1,N),RESFT(1,N),LTIMER
     &   )
      ENDDO

!-----------------------------------------------------------------------
! Calculate gridbox-means of transfer coefficients.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        CD(I) = 0.
        CH(I) = 0.

! Sea and sea-ice
        IF ( .NOT.LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I) .LT. 0.7 ) THEN
            CD(I) = ( ICE_FRACT(I)*CD_MIZ(I) +
     &                  (0.7-ICE_FRACT(I))*CD_SEA(I) ) / 0.7   ! P2430.5
            CH(I) = ( ICE_FRACT(I)*CH_MIZ(I) +
     &                (0.7-ICE_FRACT(I))*CH_SEA(I) ) / 0.7     ! P2430.4
          ELSE
            CD(I) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +
     &                  (ICE_FRACT(I)-0.7)*CD_ICE(I) ) / 0.3   ! P2430.7
            CH(I) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +
     &                  (ICE_FRACT(I)-0.7)*CH_ICE(I) ) / 0.3   ! P2430.7
          ENDIF
        ENDIF

      ENDDO

! Land tiles
      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          CD(I) = CD(I) + TILE_FRAC(L,N)*CD_TILE(L,N)
          CH(I) = CH(I) + TILE_FRAC(L,N)*CH_TILE(L,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!  4.3 Calculate the surface exchange coefficients RHOK(*) and
!       resistances for use in Sulphur Cycle
!       (Note that CD_STD, CH and VSHR should never = 0)
!     RHOSTAR * CD * VSHR stored for diagnostic output before
!     horizontal interpolation.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        RHO_ARESIST(I) = 0.
        ARESIST(I) = 0.
        RESIST_B(I) = 0.
        RHOKM_1(I) = 0.

! Sea and sea-ice
        IF ( .NOT.LAND_MASK(I) ) THEN
          RHOKM_1(I) = RHOSTAR(I)*CD(I)*VSHR(I)               ! P243.124
          RHOKH_1_SICE(I) = RHOSTAR(I) * CH(I) * VSHR(I)      ! P243.125
          RHO_ARESIST(I) = RHOSTAR(I) * CD(I) * VSHR(I)
          ARESIST(I) =  1. / (CD(I) * VSHR(I))
          RESIST_B(I)= (CD(I)/CH(I) - 1.0) * ARESIST(I)
        ENDIF

      ENDDO

! Land tiles
      DO N=1,NTYPE
        DO L=LAND1,LAND1+LAND_PTS-1
          RHO_ARESIST_TILE(L,N) = 0.
          ARESIST_TILE(L,N) = 0.
          RESIST_B_TILE(L,N) = 0.
        ENDDO
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          RHOKM_1_TILE(L,N) = RHOSTAR(I)*CD_TILE(L,N)*VSHR(I) ! P243.124
          RHOKM_1(I) = RHOKM_1(I) + TILE_FRAC(L,N)*RHOKM_1_TILE(L,N)
          RHOKH_1(L,N) = RHOSTAR(I)*CH_TILE(L,N)*VSHR(I)      ! P243.125
          RHO_ARESIST_TILE(L,N) = RHOSTAR(I) * CD_STD(L,N) * VSHR(I)
          ARESIST_TILE(L,N) = 1. / ( CD_STD(L,N) * VSHR(I) )
          RESIST_B_TILE(L,N) = ( CD_STD(L,N)/CH_TILE(L,N) - 1.0 ) *
     &                                                 ARESIST_TILE(L,N)
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
        RHO_CD_MODV1(I) = RHOKM_1(I)      ! diagnostic required for VAR
      ENDDO

!-----------------------------------------------------------------------
!!  Calculate local and gridbox-average surface fluxes of heat and
!!  moisture. Parameters for snow tile depend on whether or not a land
!!  point has permanent ice cover.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        FTL_1(I) = 0.
        FQW_1(I) = 0.
        ASHTF(I) = 2 * KAPPAI / DE
      ENDDO

      DO N=1,NTYPE
        DO L = LAND1,LAND1+LAND_PTS-1
          FTL_TILE(L,N) = 0.
          FQW_TILE(L,N) = 0.
        ENDDO
      ENDDO

      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        ASHTF(I) = 2.0 * HCONS(L) / DZSOIL
        ASHTF_SNOW(I) = ASHTF(I)
        IF ( SMVCST(L).NE.0. ) THEN
          ASHTF_SNOW(I) = 2.0 * SNOW_HCON / DEFF_SNOW
        ENDIF
      ENDDO

! Sea and sea-ice
      CALL SF_FLUX_SEA (
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,LAND_MASK,
     & ASHTF,ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1_SICE,
     & TI,TL_1,TSTAR_ICE,TSTAR_SEA,Z0_ICE,Z0_ICE,Z0H_SEA,Z0MSEA,Z1_TQ,
     & ALPHA1_SICE,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,RHOKPM_SICE,
     & LTIMER
     & )

! Snow-free land tiles
      DO N=1,NTYPE-1
        CALL SF_FLUX_LAND (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   ASHTF,LC,QS1,QSTAR_TILE(1,N),QW_1,RADNET,RESFT(1,N),
     &   RHOKH_1(1,N),TILE_FRAC(1,N),TL_1,TS1,TSTAR_TILE(1,N),
     &   Z0H_TILE(1,N),Z0M_EFF_TILE(1,N),Z1_TQ,
     &   FQW_1,FTL_1,
     &   ALPHA1(1,N),FQW_TILE(1,N),FTL_TILE(1,N),RHOKPM(1,N),LTIMER
     &   )
      ENDDO

! Snow and land-ice tile
      N=NTYPE
      CALL SF_FLUX_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     & ASHTF_SNOW,LS,QS1,QSTAR_TILE(1,N),QW_1,RADNET_SNOW,RESFT(1,N),
     & RHOKH_1(1,N),TILE_FRAC(1,N),TL_1,TSNOW,TSTAR_TILE(1,N),
     & Z0H_TILE(1,N),Z0M_EFF_TILE(1,N),Z1_TQ,
     & FQW_1,FTL_1,
     & ALPHA1(1,N),FQW_TILE(1,N),FTL_TILE(1,N),RHOKPM(1,N),LTIMER
     & )

!-----------------------------------------------------------------------
!!  4.4   Calculate the standard deviations of layer 1 turbulent
!!        fluctuations of temperature and humidity using approximate
!!        formulae from first order closure.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        Q1_SD(I) = 0.
        T1_SD(I) = 0.
      ENDDO

! Sea and sea-ice
      CALL STDEV1_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK,
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_1,RHOSTAR,VSHR,
     & Z0MSEA,Z0_ICE,Z1_TQ,
     & Q1_SD,T1_SD,LTIMER
     & )

! Land tiles
      DO N=1,NTYPE
        CALL STDEV1_LAND (
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),LAND_INDEX,TILE_INDEX(1,N),
     &   BQ_1,BT_1,FQW_TILE(1,N),FTL_TILE(1,N),RHOKM_1_TILE(1,N),
     &   RHOSTAR,VSHR,Z0M_TILE(1,N),Z1_TQ,
     &   Q1_SD,T1_SD,LTIMER
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  4.5 Set indicator for unstable suface layer (buoyancy flux +ve.).
!-----------------------------------------------------------------------
! Set to 0 - rapidly mixing boundary layer not available with MOSES II

      DO I=P1,P1+P_POINTS-1
        NRML(I) = 0
      ENDDO

!-----------------------------------------------------------------------
!!  4.6 For sea points, calculate the wind mixing energy flux and the
!!      sea-surface roughness length on the P-grid, using time-level n
!!      quantities.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1

        IF (SFME) FME(I) = 0.0
        IF (.NOT.LAND_MASK(I)) THEN
          TAU = RHOKM_1(I) * VSHR(I)                         ! P243.130
          IF (ICE_FRACT(I) .GT. 0.0)
     &      TAU = RHOSTAR(I) * CD_SEA(I) * VSHR(I) * VSHR(I)

          IF (SFME) FME(I) = (1.0-ICE_FRACT(I)) * TAU * SQRT(TAU/RHOSEA)
!                                                             ! P243.96
          Z0MSEA(I) = MAX ( Z0HSEA ,
     &                      (CHARNOCK/G) * (TAU / RHOSTAR(I)) )
!                                         ... P243.B6 (Charnock formula)
!                      TAU/RHOSTAR is "mod VS squared", see eqn P243.131
        ENDIF

      ENDDO

!-----------------------------------------------------------------------
! Calculate effective roughness lengths, orographic blending heights
! and gridbox-average Richardson numbers.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        RIB(I) = 0.
        Z0M_EFF(I) = 1.

! Sea and sea-ice (leads ignored)
        IF ( .NOT.LAND_MASK(I) ) THEN
          H_BLEND_OROG(I) = H_BLEND_MIN
          RIB(I) = RIB_SEA(I)
          Z0M_EFF(I) = Z0MSEA(I)
          Z0M(I) = Z0MSEA(I)
          Z0H(I) = Z0HSEA
          IF ( ICE_FRACT(I) .GT. 0. ) THEN
            RIB(I) = RIB_ICE(I)
            Z0M_EFF(I) = Z0_ICE(I)
            Z0M(I) = Z0_ICE(I)
            Z0H(I) = Z0_ICE(I)
          ENDIF
        ENDIF

      ENDDO

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          RIB(I) = RIB(I) + TILE_FRAC(L,N)*RIB_TILE(L,N)
        ENDDO
      ENDDO


      DO L = LAND1,LAND1+LAND_PTS-1
        Z0_GB(L) = Z0_SF_GB(LAND_INDEX(L))
      ENDDO
      DO J=1,TILE_PTS(NTYPE)
        L = TILE_INDEX(J,NTYPE)
        Z0 = TILE_FRAC(L,NTYPE) / ( LOG(LB/Z0M_TILE(L,NTYPE))**2 ) +
     &               (1. - TILE_FRAC(L,NTYPE)) / ( LOG(LB/Z0_GB(L))**2 )
        Z0_GB(L) = LB * EXP( - SQRT(1./Z0) )
      ENDDO

      CALL SF_OROG_GB(
     & P_FIELD,P1,P_POINTS,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,
     & LAND_MASK,L_Z0_OROG,HO2R2_OROG,RIB,SIL_OROG,Z0_GB,Z1_UV,
     & H_BLEND_OROG,Z0M_EFF,LTIMER
     & )

!-----------------------------------------------------------------------
!! Call SFL_INT to calculate CDR10M and CHR1P5M - interpolation coeffs
!! used to calculate screen temperature, humidity and 10m winds.
!-----------------------------------------------------------------------

      IF (SU10 .OR. SV10 .OR. SQ1P5 .OR. ST1P5) THEN

! Sea and sea-ice (leads ignored)
        DO I=P1,P1+P_POINTS-1
          CDR10M(I) =0.
          IF ( .NOT.LAND_MASK(I) .AND. ICE_FRACT(I).GT.0. ) THEN
            CD_SEA(I) = CD_ICE(I)
            CH_SEA(I) = CH_ICE(I)
            Z0H_SEA(I) = Z0_ICE(I)
            Z0F_SEA(I) = Z0_ICE(I)
          ENDIF
        ENDDO

        CALL SFL_INT_SEA (
     &   P_POINTS,P_FIELD,P1,
     &   CD_SEA,CH_SEA,RIB,Z0M_EFF,Z0H_SEA,Z0F_SEA,Z1_UV,
     &   LAND_MASK,SU10,SV10,ST1P5,SQ1P5,LTIMER,
     &   CDR10M,CHR1P5M_SICE
     &   )

! Land tiles
        DO N=1,NTYPE
          CALL SFL_INT_LAND (
     &     P_FIELD,LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N),LAND_INDEX,
     &     CD_STD(1,N),CD_TILE(1,N),CH_TILE(1,N),RIB_TILE(1,N),
     &     TILE_FRAC(1,N),WIND_PROFILE_FACTOR(1,N),Z0M_TILE(1,N),
     &     Z0M_EFF_TILE(1,N),Z0H_TILE(1,N),Z0F_TILE(1,N),Z1_UV,
     &     SU10,SV10,ST1P5,SQ1P5,LTIMER,
     &     CDR10M,CHR1P5M(1,N)
     &     )
        ENDDO

      ENDIF

      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',4)
      ENDIF

      RETURN
      END
