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
!!!  SUBROUTINE BDY_LAYR-----------------------------------------------
!!!
!!!  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
!!!           between (a) surface and atmosphere, (b) atmospheric levels
!!!           within the boundary layer, and/or the effects of these
!!!           fluxes on the primary model variables.  The flux of heat
!!!           into and through the soil is also modelled.  Numerous
!!!           related diagnostics are also calculated.
!!!           F E Hewer, July 1990: removed call to LS_CLD.
!!!    This version passes out liquid/frozen water temperature in
!!!    array "T" (TL), and total water content in array "Q" (QW).
!!!    These may be converted to T and Q respectively by calling
!!!    the large scale cloud routine, LS_CLD.
!!!            F E Hewer, August 1990: land point data stored
!!!    on land points only.
!!!    Arrays whose elements may contain values over both sea and land
!!!    points are compressed onto land points for land calculations if
!!!    defined variable IBM is NOT selected. RHOKM,RHOKH redefined as
!!!    workspace.
!!!
!!!
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.3  7/2/97     New deck. S Jackson
!!!   4.4 25/6/97     Modified for MOSES II tile model. R Essery
!!!   4.4 25/6/97     Move grid definitions up to BL_INTCT.  R.A.Betts
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!   4.5  7/5/98     Set TSTAR, SNOW_SURF_HTF and SOIL_SURF_HTF to 0
!!!                   at all land points, to avoid problems of
!!!                   non-initialised data.  R.A.Betts
!!!   4.5 21/5/98     Add optional error check for negative surface
!!!                   temperature.  R.A.Betts
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_LAYR (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,LAND_FIELD,
     & P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH,
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,LAND1,LAND_PTS,U_POINTS,U1,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN soil/vegetation/land surface data :
     & LAND_INDEX,
     & LAND_MASK,L_Z0_OROG,
     & NTYPE,TILE_INDEX,TILE_PTS,SM_LEVELS,
     & CANOPY,CATCH,GC,HCON,HO2R2_OROG,LYING_SNOW,
     & SIL_OROG_LAND,SMC,SMVCST,STHF,STHU,
     & TILE_FRAC,WT_EXT,Z0_SF_GB,Z0_TILE,

! IN sea/sea-ice data :
     & DI,ICE_FRACT,U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & PSTAR,RADNET,RADNET_SNOW,TIMESTEP,VSHR,
     & L_RMBL,L_BL_LSPICE,L_MOM,L_NEG_TSTAR,

! IN STASH flags :-
     & SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,

! INOUT data :
     & Q,T,T_SOIL,TSNOW,TI,TSTAR,TSTAR_TILE,U,V,Z0MSEA,

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,ECAN,E_SEA,ESOIL_TILE,FQW,
     & FTL,FTL_TILE,H_SEA,RHOKH,RHOKM_UV,
     & RIB,RIB_TILE,SEA_ICE_HTF,SURF_HT_FLUX,TAUX,TAUY,

! OUT diagnostic requiring STASH flags :
     & FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,
     & Q1P5M,T1P5M,U10M,V10M,

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & NRML,

! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,RHO_KM,

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ESOIL,EXT,SNOWMELT,ZH,
     & SOIL_SURF_HTF,SNOW_SURF_HTF,
     & T1_SD,Q1_SD,ERROR,

! LOGICAL LTIMER
     & LTIMER
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid
!                                  !    (for dimensioning only).
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
     &,LAND_FIELD                  ! IN No.of land points in whole grid.
     &,P_ROWS                      ! IN No. of P-rows in whole grid
!                                  !    (for dimensioning only).
     &,FIRST_ROW                   ! IN First row of data to be treated,
!                                  !    referred to P-grid.
     &,N_ROWS                      ! IN No. of rows of data to be
!                                  !    treated, referred to P-grid.
     &,ROW_LENGTH                  ! IN No. of points in one row.
     &,N_P_ROWS   ! IN No of P-rows being processed.
     &,N_U_ROWS   ! IN No of UV-rows being processed.
     &,P_POINTS   ! IN No of P-points being processed.
     &,P1         ! IN First P-point to be processed.
     &,LAND1      ! IN First land-point to be processed.
!                 !       1 <= LAND1 <= LAND_FIELD
     &,LAND_PTS   ! IN No of land points being processed.
     &,U_POINTS   ! IN No of UV-points being processed.
     &,U1         ! IN First UV-point to be processed.

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                  !    allowed. Assumed <= 30 for dim-
!                                  !    ensioning GAMMA in common deck
!                                  !    C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                    ! IN Total no. of vertical levels in
!                                  !    the model atmosphere.
      REAL
     & AK(P_LEVELS)                ! IN Hybrid 'A' for all levels.
     &,BK(P_LEVELS)                ! IN Hybrid 'B' for all levels.
     &,AKH(P_LEVELS+1)             ! IN Hybrid 'A' for layer interfaces.
     &,BKH(P_LEVELS+1)             ! IN Hybrid 'B' for layer interfaces.
     &,DELTA_AK(P_LEVELS)          ! IN Difference of hybrid 'A' across
!                                  !    layers (K-1/2 to K+1/2).
!                                  !    NB: Upper minus lower.
     &,DELTA_BK(P_LEVELS)          ! IN Difference of hybrid 'B' across
!                                  !     layers (K-1/2 to K+1/2).
!                                  !     NB: Upper minus lower.
     &,EXNER(P_FIELD,BL_LEVELS+1)  ! IN Exner function.  EXNER(,K) is
!                                  !    value for LOWER BOUNDARY of
!                                  !    level K.

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL
     & LAND_MASK(P_FIELD)          ! IN T if land, F elsewhere.
     &,L_Z0_OROG                   ! IN T to use orog.roughness
!                                  !    treatment in SF_EXCH

      INTEGER
     & LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in P_FIELD is the Ith
!                                  !    land point.

      INTEGER
     & SM_LEVELS                   ! IN No. of soil moisture levels
     &,NTYPE                       ! IN No. of land tiles
     &,TILE_INDEX(LAND_FIELD,NTYPE)! IN Index of tile points
     &,TILE_PTS(NTYPE)             ! IN Number of tile points

      REAL
     & CANOPY(LAND_FIELD,NTYPE-1)  ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_FIELD,NTYPE-1)   ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,GC(LAND_FIELD,NTYPE)        ! IN "Stomatal" conductance to
!                                  !     evaporation for land tiles
!                                  !     (m/s).
     &,HCON(LAND_FIELD)            ! IN Soil thermal conductivity
!                                  !    (W/m/K).
     &,LYING_SNOW(P_FIELD)         ! IN Lying snow (kg/sq m).
!                                     Must be global for coupled model,
!                                     ie dimension P_FIELD not
!                                     LAND_FIELD
     &,SMC(LAND_FIELD)             ! IN Available soil moisture (kg/m2).
     &,SMVCST(LAND_FIELD)          ! IN Volumetric saturation point
!                                  !    (m3/m3 of soil).
     &,STHF(LAND_FIELD,SM_LEVELS)  ! IN Frozen soil moisture content of
!                                  !    each layer as a fraction of
!                                  !    saturation.
     &,STHU(LAND_FIELD,SM_LEVELS)  ! IN Unfrozen soil moisture content
!                                  !    of each layer as a fraction of
!                                  !    saturation.
     &,TILE_FRAC(LAND_FIELD,NTYPE) ! IN Tile fractions including
!                                  ! snow cover in the ice tile.
     &,WT_EXT(LAND_FIELD,SM_LEVELS)! IN Fraction of evapotranspiration
!                                  !    extracted from each soil layer.
     &,Z0_TILE(LAND_FIELD,NTYPE)   ! IN Tile roughness lengths (m).
     &,Z0_SF_GB(P_FIELD)           ! IN GBM roughness length for
!                                  !    snow-free land (m).
     &,SIL_OROG_LAND(LAND_FIELD)   ! IN Silhouette area of unresolved
!                                  !    orography per unit horizontal
!                                  !    area on land points only.
     &,HO2R2_OROG(LAND_FIELD)      ! IN Standard Deviation of orography.
!                                  !    equivilent to peak to trough
!                                  !    height of unresolved orography
!                                  !    divided by 2SQRT(2) on land
!                                  !    points only (m)

! (d) Sea/sea-ice data.

      REAL
     & DI(P_FIELD)                 ! IN "Equivalent thickness" of
!                                  !     sea-ice(m).
     &,ICE_FRACT(P_FIELD)          ! IN Fraction of gridbox covered by
!                                  !     sea-ice (decimal fraction).
     &,U_0(U_FIELD)                ! IN W'ly component of surface
!                                  !    current (m/s).
     &,V_0(U_FIELD)                ! IN S'ly component of surface
!                                  !    current (m/s).

! (e) Cloud data.

      REAL
     & CF(P_FIELD,BL_LEVELS)       ! IN Cloud fraction (decimal).
     &,QCF(P_FIELD,BL_LEVELS)      ! IN Cloud ice (kg per kg air)
     &,QCL(P_FIELD,BL_LEVELS)      ! IN Cloud liquid water (kg
!                                  !    per kg air).
     &,CCA(P_FIELD)                ! IN Convective Cloud Amount
!                                  !    (decimal)

      INTEGER
     & CCB(P_FIELD)                ! IN Convective Cloud Base
     &,CCT(P_FIELD)                ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL
     & PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).
     &,RADNET(P_FIELD)             ! IN Surface net radiation for sea-
!                                  !    ice or snow-free land (W/sq m).
     &,RADNET_SNOW(P_FIELD)        ! IN Snow surface net radiation.
     &,TIMESTEP                    ! IN Timestep (seconds).
     &,VSHR(P_FIELD)               ! IN Magnitude of surface-to-lowest
!                                  !    atm level wind shear (m per s).

      LOGICAL
     & LTIMER                      ! IN Logical switch for TIMER diags
     &,L_RMBL                      ! IN T to use rapidly mixing boundary
!                                  !    scheme
!                                  !    - not available in MOSES II
     &,L_BL_LSPICE                 ! IN Use if 3A large scale precip
     &,L_MOM                       ! IN Switch for convective momentum
!                                  !    transport.
     &,L_NEG_TSTAR                ! IN Switch for -ve TSTAR error check

!  STASH flags :-

      LOGICAL
     & SFME    ! IN Flag for FME (q.v.).
     &,SIMLT   ! IN Flag for SICE_MLT_HTF (q.v.)
     &,SMLT    ! IN Flag for SNOMLT_SURF_HTF (q.v.)
     &,SLH     ! IN Flag for LATENT_HEAT (q.v.)
     &,SQ1P5   ! IN Flag for Q1P5M (q.v.)
     &,ST1P5   ! IN Flag for T1P5M (q.v.)
     &,SU10    ! IN Flag for U10M (q.v.)
     &,SV10    ! IN Flag for V10M (q.v.)

!  In/outs :-

      REAL
     & Q(P_FIELD,BL_LEVELS)        ! IN  Specific humidity ( kg/kg air).
!                                  ! OUT Total water content (QW)
!                                  !     (kg/kg air).
     &,T(P_FIELD,BL_LEVELS)        ! IN  Atmospheric temperature (K).
!                                  ! OUT Liquid/frozen water
!                                  !     temperature (TL) (K).
     &,T_SOIL(LAND_FIELD,SM_LEVELS)! INOUT Soil temperatures (K).
     &,TI(P_FIELD)                 ! INOUT Sea-ice surface layer
!                                  !       temperature (K).
     &,TSNOW(LAND_FIELD)           ! INOUT Snow surface layer
!                                  !       temperature (K).
!                                  !       =T_SOIL(*,1) for land-ice
     &,TSTAR(P_FIELD)              ! INOUT GBM surface temperature (K).
     &,TSTAR_TILE(LAND_FIELD,NTYPE)! INOUT Surface tile temperatures
     &,U(U_FIELD,BL_LEVELS)        ! INOUT W'ly wind component (m/s)
     &,V(U_FIELD,BL_LEVELS)        ! INOUT S'ly wind component (m/s)
     &,Z0MSEA(P_FIELD)             ! INOUT Sea-surface roughness
!                                  !       length for momentum (m).

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL
     & CD(P_FIELD)                 ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     momentum.
     &,CH(P_FIELD)                 ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     heat and/or moisture.
     &,ECAN(P_FIELD)               ! OUT Gridbox mean evaporation from
!                                  !     canopy/surface store (kg/m2/s).
!                                  !     Zero over sea.
     &,E_SEA(P_FIELD)              ! OUT Evaporation from sea times
!                                  !     leads fraction. Zero over land.
!                                  !     (kg per square metre per sec).
     &,ESOIL_TILE(LAND_FIELD,NTYPE-1)
!                                  ! OUT ESOIL for snow-free land tiles
     &,FQW(P_FIELD,BL_LEVELS)      ! OUT Moisture flux between layers
!                                  !     (kg per square metre per sec).
!                                  !     FQW(,1) is total water flux
!                                  !     from surface, 'E'.
     &,FTL(P_FIELD,BL_LEVELS)      ! OUT FTL(,K) contains net turbulent
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the
!                                  !     surface sensible heat, H.(W/m2)
     &,FTL_TILE(LAND_FIELD,NTYPE)  ! OUT Surface FTL for land tiles
     &,H_SEA(P_FIELD)              ! OUT Surface sensible heat flux over
!                                  !     sea times leads fraction (W/m2)
     &,RHOKH(P_FIELD,BL_LEVELS)    ! OUT Exchange coeffs for moisture.
     &,RHOKM_UV(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                                  !     momentum (on UV-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RIB(P_FIELD)                ! OUT Mean bulk Richardson number for
!                                  !     lowest layer.
     &,RIB_TILE(LAND_FIELD,NTYPE)  ! OUT RIB for land tiles.
     &,SEA_ICE_HTF(P_FIELD)        ! OUT Heat flux through sea-ice
!                                  !     (W/m2, positive downwards).
     &,SURF_HT_FLUX(P_FIELD)       ! OUT Net downward heat flux at
!                                  !     surface over land or sea-ice
!                                  !     fraction of gridbox (W/m2).
     &,TAUX(U_FIELD,BL_LEVELS)     ! OUT W'ly component of surface wind
!                                  !     stress (N/sq m). (On UV-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUY(U_FIELD,BL_LEVELS)     ! OUT S'ly component of surface wind
!                                  !     stress (N/sq m).  On UV-grid;
!                                  !     comments as per TAUX.
     &,RHO_CD_MODV1(P_FIELD)       ! OUT Surface air density * drag coef
!                                  !     *mod(v1 - v0) before interp
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT Air density * turbulent mixing
!                                  !     coefficient for momentum before
!                                  !     interpolation.
     &,RHO_ARESIST(P_FIELD)        ! OUT RHOSTAR*CD_STD*VSHR for Sulphur
!                                  !     cycle
     &,ARESIST(P_FIELD)            ! OUT 1/(CD_STD*VSHR) for Sulphur
!                                  !     cycle
     &,RESIST_B(P_FIELD)           ! OUT (1/CH-1/(CD_STD)/VSHR for
!                                  !     Sulphur cycle
     &,RHO_ARESIST_TILE(LAND_FIELD,NTYPE)
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
     &,ARESIST_TILE(LAND_FIELD,NTYPE)
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_FIELD,NTYPE)
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles

      INTEGER
     & NRML(P_FIELD)               ! OUT Number of model layers in the
!                                  !     Rapidly Mixing Layer; set to
!                                  !     zero in SF_EXCH for MOSES II.

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL
     & FME(P_FIELD)                ! OUT Wind mixing "power" (W/m2).
     &,SICE_MLT_HTF(P_FIELD)       ! OUT Heat flux due to melting of
!                                  !     sea-ice (Watts per sq metre).
     &,SNOMLT_SURF_HTF(P_FIELD)    ! OUT Heat flux required for surface
!                                  !     melting of snow (W/m2).
     &,LATENT_HEAT(P_FIELD)        ! OUT Surface latent heat flux, +ve
!                                  !     upwards (Watts per sq m).
     &,Q1P5M(P_FIELD)              ! OUT Q at 1.5 m (kg water / kg air).
     &,T1P5M(P_FIELD)              ! OUT T at 1.5 m (K).
     &,U10M(U_FIELD)               ! OUT U at 10 m (m per s).
     &,V10M(U_FIELD)               ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL
     & EI(P_FIELD)                 ! OUT Sublimation from lying snow or
!                                  !     sea-ice (kg/m2/s).
     &,ECAN_TILE(LAND_FIELD,NTYPE-1)! OUT ECAN for snow-free land tiles
     &,ESOIL(P_FIELD)              ! OUT Surface evapotranspiration
!                                  !     from soil moisture store
!                                  !     (kg/m2/s).
     &,EXT(LAND_FIELD,SM_LEVELS)   ! OUT Extraction of water from each
!                                  !     soil layer (kg/m2/s).
     &,SOIL_SURF_HTF(LAND_FIELD)   ! OUT Net downward heat flux at
!                                  !     snow-free land surface (W/m2).
     &,SNOW_SURF_HTF(LAND_FIELD)   ! OUT Net downward heat flux at
!                                  !     snow surface (W/m2).
     &,SNOWMELT(P_FIELD)           ! OUT Snowmelt (kg/m2/s).
     &,ZH(P_FIELD)                 ! OUT Height above surface of top of
!                                  !     boundary layer (metres).
     &,T1_SD(P_FIELD)              ! OUT Standard deviation of turbulent
!                                  !     fluctuations of layer 1 temp;
!                                  !     used in initiating convection.
     &,Q1_SD(P_FIELD)              ! OUT Standard deviation of turbulent
!                                  !     flucs of layer 1 humidity;
!                                  !     used in initiating convection.
      INTEGER
     & ERROR          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL Z,HEAT_CON,SF_EXCH,BOUY_TQ,BTQ_INT,
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,SF_MELT,
     & IM_CAL_UV,SCREEN_TQ
      EXTERNAL TIMER
      EXTERNAL UV_TO_P,P_TO_UV

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

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
C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
      REAL
     + DZSOIL(4)               ! Soil layer thicknesses (m).
      DATA DZSOIL /0.100, 0.250, 0.650, 2.000 /
C-----------------------------------------------------------------------
! MPP Common block
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the MPP-UM
!
!   Two sets of parameters are set up -
!     i)  for the MPP-UM itself.
!     ii) for the interface to the Message Passing Software.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.1      27/1/96  New comdeck based on first section of
!                     old PARVARS.   P.Burton
!   4.2      21/11/96 Add new field type parameter and
!                     magic number used in addressing to indicate
!                     if a calculation is for local data, or data
!                     on the dump on disk (ie. global data)  P.Burton
!   4.2      18/11/96 Moved MaxFieldSize to comdeck AMAXSIZE and
!                     removed Maxbuf.  P.Burton
!   4.2      18/7/96  Removed some unused variables      P.Burton
!   4.4      11/07/97 Reduced MAXPROC to 256 to save memory  P.Burton
!
! ---------------------- PARAMETERS ---------------------
!
! =======================================================
! Parameters needed for the MPP-UM
! =======================================================

      INTEGER   Ndim_max        ! maximum number of spatial dimensions
      PARAMETER (Ndim_max = 3 ) ! 3d data


      INTEGER
     &   fld_type_p           ! indicates a grid on P points
     &,  fld_type_u           ! indicates a grid on U points
     &,  fld_type_unknown     ! indicates a non-standard grid.
      PARAMETER (
     &   fld_type_p=1
     &,  fld_type_u=2
     &,  fld_type_unknown=-1)

      INTEGER
     &   local_data
     &,  global_dump_data
      PARAMETER (
     &   local_data=1        ! Used in addressing to indicate if
     &,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================


      INTEGER
     &   Maxproc              ! Max number of processors
      PARAMETER (
     &   MAXPROC = 256)

      INTEGER
     &   PNorth       ! North processor address in the neighbour array
     &,  PEast        ! East  processor address in the neighbour array
     &,  PSouth       ! South processor address in the neighbour array
     &,  PWest        ! West  processor address in the neighbour array
     &,  NoDomain     ! Value in neighbour array if the domain has
     &                !  no neighbor in this direction. Otherwise
     &                !  the value will be the tid of the neighbor
      PARAMETER (
     &   PNorth   = 1
     &,  PEast    = 2
     &,  PSouth   = 3
     &,  PWest    = 4
     &,  NoDomain = -1)

      INTEGER
     &   BC_STATIC            ! Static boundary conditions
     &,  BC_CYCLIC            ! Cyclic boundary conditions
      PARAMETER (
     &   BC_STATIC = 1
     &,  BC_CYCLIC = 2)

! ---------------------- End of comdeck PARPARM ---------------------
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the MPP-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for MPP-UM
!                                                      P.Burton
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER
     &   first_comp_pe       ! top left pe in LPG
     &,  last_comp_pe        ! bottom right pe in LPG
     &,  current_decomp_type ! current decomposition type
     &,  Offx                ! halo size in EW direction
     &,  Offy                ! halo size in NS direction
     &,  glsize(Ndim_max)    ! global data size
     &,  lasize(Ndim_max)    ! local data size
     &,  blsizep(Ndim_max)   ! personal p data area
     &,  blsizeu(Ndim_max)   ! personal u data area
     &,  datastart(Ndim_max) ! position of personal data in global data
     &                       !   (in terms of standard Fortran array
     &                       !    notation)
     &,  gridsize(Ndim_max)  ! size of the LPG in each dimension
     &,  gridpos(Ndim_max)   ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

      LOGICAL
     &    atbase             ! process at the bottom of the LPG
     &,   attop              ! process at the top of the LPG
     &,   atleft             ! process at the left of the LPG
     &,   atright            ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

      COMMON /UM_PARVAR/
     &                  first_comp_pe,last_comp_pe
     &,                 current_decomp_type,Offx, Offy
     &,                 glsize,lasize,blsizep,blsizeu
     &,                 datastart,gridsize,gridpos
     &,                 atbase,attop,atleft,atright

! =======================================================
! Common block for the Message Passing Software
! =======================================================

      INTEGER
     &  bound(Ndim_max)           ! type of boundary (cyclic or static)
     &                            !  in each direction
     &, g_lasize(Ndim_max,0:maxproc)
!                                 ! global copy of local data size
     &, g_blsizep(Ndim_max,0:maxproc)
!                                 ! global copy of personal p data area
     &, g_blsizeu(Ndim_max,0:maxproc)
!                                 ! global copy of personal u data area
     &, g_datastart(Ndim_max,0:maxproc)
!                                 ! global copy of datastart
     &, g_gridpos(Ndim_max,0:maxproc)
!                                 ! global copy of gridpos
     &, nproc                     ! number of processors in current
!                                 ! decomposition
     &, nproc_max                 ! maximum number of processors
     &, nproc_x                   ! number of processors in x-direction
     &, nproc_y                   ! number of processors in y-direction
     &, mype                      ! number of this processor
     &                            !  (starting from 0)
     &, neighbour(4)              ! array with the tids of the four
     &                            ! neighbours in the horizontal plane
     &, gc_proc_row_group         ! GID for procs along a proc row
     &, gc_proc_col_group         ! GID for procs along a proc col
     &, gc_all_proc_group         ! GID for all procs

      COMMON /MP_PARVAR/
     &                  bound
     &,                 g_lasize,g_blsizep,g_blsizeu
     &,                 g_datastart,g_gridpos
     &,                 nproc,nproc_max,nproc_x,nproc_y,mype
     &,                 neighbour,gc_proc_row_group
     &,                 gc_proc_col_group, gc_all_proc_group



! ---------------------- End of comdeck PARCOMM -----------------------
! --------------------- End of comdeck PARVARS ---------------------

! Derived local parameters.

      REAL LCRCP,LS,LSRCP

      PARAMETER (
     & LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     &  )

!-----------------------------------------------------------------------

!  Workspace :-

      REAL
     & ALPHA1(LAND_FIELD,NTYPE) ! Mean gradient of saturated
!                               ! specific humidity with respect to
!                               ! temperature between the bottom model
!                               ! layer and tile surfaces
     &,ALPHA1_SICE(P_FIELD)     ! ALPHA1 for sea-ice.
     &,ASHTF(P_FIELD)           ! Coefficient to calculate surface
!                               ! heat flux into soil or sea-ice.
     &,ASHTF_SNOW(P_FIELD)      ! ASHTF for snow or land-ice.
     &,ASURF(P_FIELD)           ! Reciprocal areal heat capacity
!                               ! of sea-ice surface layer (K m**2 / J).
     &,BF(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta F tilde)
     &,BQ(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta q tilde).
     &,BT(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta T tilde).
     &,DELTAP(P_FIELD,BL_LEVELS)! Difference in pressure between levels
     &,DELTAP_UV(P_FIELD,BL_LEVELS)
!                               ! Difference in pressure between levels
!                               ! on UV points
     &,DTRDZ(P_FIELD,BL_LEVELS) ! -g.dt/dp for model layers.
     &,DTRDZ_UV(U_FIELD,BL_LEVELS)
!                               ! -g.dt/dp for model wind layers.
     &,DTRDZ_RML(P_FIELD)       ! -g.dt/dp for the rapidly
!                               !  mixing layer.
     &,DZL(P_FIELD,BL_LEVELS)   ! DZL(,K) is depth in m of layer
!                               ! K, i.e. distance from boundary
!                               ! K-1/2 to boundary K+1/2.
     &,DU(U_FIELD,BL_LEVELS)    ! BL increment to u wind foeld
     &,DV(U_FIELD,BL_LEVELS)    ! BL increment to v wind foeld
     &,DU_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to u wind field
     &,DV_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to v wind field
     &,DTL_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to TL field
     &,DQW_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to QW field
     &,FQW_TILE(LAND_FIELD,NTYPE)! Surface FQW for land tiles
     &,FQW_ICE(P_FIELD)         ! Surface FQW for sea-ice
     &,FTL_ICE(P_FIELD)         ! Surface FTL for sea-ice
     &,FRACA(LAND_FIELD,NTYPE-1)! Fraction of surface moisture flux
!                               ! with only aerodynamic resistance
!                               ! for snow-free land tiles.
     &,HCONS(LAND_FIELD)        ! Soil thermal conductivity including
!                               ! the effects of water and ice (W/m2)
     &,QW(P_FIELD,BL_LEVELS)    ! Total water content, but
!                               ! replaced by specific humidity
!                               ! in LS_CLD.
     &,P(P_FIELD,BL_LEVELS)     ! Pressure at model levels
     &,RDZ(P_FIELD,BL_LEVELS)   ! RDZ(,1) is the reciprocal of the
!                               ! height of level 1, i.e. of the
!                               ! middle of layer 1.  For K > 1,
!                               ! RDZ(,K) is the reciprocal
!                               ! of the vertical distance
!                               ! from level K-1 to level K.
     &,RDZUV(U_FIELD,BL_LEVELS) ! RDZ (K > 1) on UV-grid.
!                               ! Comments as per RHOKM (RDZUV).
     &,RESFS(LAND_FIELD,NTYPE-1)! Combined soil, stomatal
!                               ! and aerodynamic resistance
!                               ! factor for fraction (1-FRACA) of
!                               ! snow-free land tiles.
     &,RESFT(LAND_FIELD,NTYPE)  ! Total resistance factor.
!                               ! FRACA+(1-FRACA)*RESFS for snow-free
!                               ! land, 1 for snow.
     &,RHO(P_FIELD,BL_LEVELS)   ! Density of model layer
     &,RHOKH_TILE(LAND_FIELD,NTYPE)
!                               ! Surface exchange coefficients
!                               ! for land tiles
     &,RHOKH_SICE(P_FIELD)      ! Surface exchange coefficients
!                               ! for sea and sea-ice
     &,RHOKM(P_FIELD,BL_LEVELS) ! Exchange coefficients for
!                               ! momentum on P-grid
     &,RHOKPM(LAND_FIELD,NTYPE) ! Land surface exchange coeff.
     &,RHOKPM_SICE(P_FIELD)     ! Sea-ice surface exchange coeff.
     &,TL(P_FIELD,BL_LEVELS)    ! Ice/liquid water temperature,
!                               ! but replaced by T in LS_CLD.
     &,TV(P_FIELD,BL_LEVELS)    ! Virtual temp
     &,U_P(P_FIELD,BL_LEVELS)   ! U on P-grid.
     &,V_P(P_FIELD,BL_LEVELS)   ! V on P-grid.
     &,ZLB(P_FIELD,0:BL_LEVELS) ! ZLB(,K) is the height of the
!                               ! upper boundary of layer K
!                               ! ( = 0.0 for "K=0").
       REAL
     & Z1(P_FIELD)              ! Height of lowest level (i.e.
!                               ! height of middle of lowest
!                               ! layer).
     &,H_BLEND_OROG(P_FIELD)    ! Blending height used as part of
!                               ! effective roughness scheme
     &,Z0H(P_FIELD)             ! Roughness length for heat and
!                               ! moisture (m).
     &,Z0H_TILE(LAND_FIELD,NTYPE)
!                               ! Tile roughness lengths for heat and
!                               ! moisture (m).
     &,Z0M(P_FIELD)             ! Roughness length for momentum (m).
     &,Z0M_TILE(LAND_FIELD,NTYPE)
!                               ! Tile roughness lengths for momentum.
     &,Z0M_EFF(P_FIELD)         ! Effective grid-box roughness
!                               ! length for momentum
     &,CDR10M(P_FIELD)          ! Ratio of CD's reqd for calculation
!                               ! of 10 m wind. On P-grid
     &,CDR10M_UV(U_FIELD)       ! Ratio of CD's reqd for calculation
!                               ! of 10 m wind. On UV-grid; comments as
!                               ! per RHOKM.
     &,CHR1P5M(LAND_FIELD,NTYPE)! Ratio of coefffs for calculation of
!                               ! 1.5m temp for land tiles.
     &,CHR1P5M_SICE(P_FIELD)    ! CHR1P5M for sea and sea-ice
!                               ! (leads ignored).


!  Local scalars :-

      REAL
     & WK         ! LOCAL 0.5 * DZL(I,K) * RDZ(I,K)
     &,WKM1       ! LOCAL 0.5 * DZL(I,K-1) * RDZ(I,K)

      INTEGER
     & I,J,L      ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).
     &,N          ! LOCAL Loop counter (tile index).

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF
      ERROR = 0

!-----------------------------------------------------------------------
!! 1.  Perform calculations in what the documentation describes as
!!     subroutine Z_DZ.  In fact, a separate subroutine isn't used.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 1.1 Initialise ZLB(,0) (to zero, of course, this being the height
!!     of the surface above the surface).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
      ENDDO

!-----------------------------------------------------------------------
!! 1.2 Calculate layer depths and heights, and construct wind fields on
!!     P-grid.  This involves calling subroutines Z and UV_TO_P.
!!     Virtual temperature is also calculated, as a by-product.
!-----------------------------------------------------------------------
!  NB RDZ  TEMPORARILY used to return DELTA_Z_LOWER, the lower half
!     layer thickness

      DO K=1,BL_LEVELS
        CALL Z(P_POINTS,EXNER(P1,K),EXNER(P1,K+1),PSTAR(P1),
     &    AKH(K),BKH(K),Q(P1,K),QCF(P1,K),
     &    QCL(P1,K),T(P1,K),ZLB(P1,K-1),TV(P1,K),
     &    ZLB(P1,K),DZL(P1,K),RDZ(P1,K),LTIMER)

        CALL UV_TO_P(U(U1,K),U_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V(U1,K),V_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)

! DZL can contain incorrect data in halos, so call SWAPBOUNDS.
      CALL SWAPBOUNDS(DZL(P1,1),ROW_LENGTH,N_U_ROWS,1,0,BL_LEVELS)

! du_nt 'borrowed to store dzl on uv grid
        CALL P_TO_UV (DZL(P1,K),DU_NT(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

      ENDDO

! set pressure array.
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          P(I,K) = AK(K) + BK(K)*PSTAR(I)

! These will be used in new dynamics scheme - currently unused
          DTL_NT(I,K)=0.0
          DQW_NT(I,K)=0.0

        ENDDO

      ENDDO  ! end of loop over bl_levels

      DO K=BL_LEVELS,2,-1

        DO I=P1,P1+P_POINTS-1
          RDZ(I,K)=1.0/(RDZ(I,K)+(DZL(I,K-1)-RDZ(I,K-1)))
          DELTAP(I,K)=DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)

          DTRDZ(I,K) = -G * TIMESTEP/ DELTAP(I,K)
!     &                  (DELTA_AK(K) + PSTAR(I)*DELTA_BK(K))
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
        Z1(I)=RDZ(I,1)
        RDZ(I,1)=1.0/RDZ(I,1)
        DELTAP(I,1)=DELTA_AK(1) + PSTAR(I)*DELTA_BK(1)
        DTRDZ(I,1) = -G * TIMESTEP/DELTAP(I,1)
!     &                  (DELTA_AK(1) + PSTAR(I)*DELTA_BK(1))
      ENDDO

      DO K=1,BL_LEVELS

! Calculate RDZUV here

        IF(K.GE.2)THEN

          DO I=U1+ROW_LENGTH,U1-ROW_LENGTH+U_POINTS-1
            RDZUV(I,K) = 2.0 / ( DU_NT(I,K) + DU_NT(I,K-1) )
          ENDDO

!-----------------------------------------------------------------------
! 1.3 Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------

      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
      ENDIF

        ENDIF   ! K .ge. 2

! Calculate DTRDZ_UV here.

!        CALL P_TO_UV (DTRDZ(P1,K),DTRDZ_UV(U1+ROW_LENGTH,K),
!     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        CALL P_TO_UV (DELTAP(P1,K),DELTAP_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DTRDZ_UV(I,K) = -G * TIMESTEP / DELTAP_UV(I,K)
        ENDDO


      ENDDO ! loop over bl_levels

! "borrowed" du_nt reset to zero
! Non turbulent increments for new dynamics scheme (currently not used)
        DO K=1,BL_LEVELS
          DO I=1,U_FIELD
            DU_NT(I,K) =0.0
            DV_NT(I,K) =0.0
          ENDDO
        ENDDO

      IF (LAND_FIELD.GT.0) THEN    ! Omit if no land points

!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.
!-----------------------------------------------------------------------
        CALL HEAT_CON (LAND_FIELD,HCON,STHU,STHF,SMVCST,HCONS,LTIMER)

      ENDIF                     ! End test on land points

!-----------------------------------------------------------------------
!! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          QW(I,K) = Q(I,K) + QCL(I,K) + QCF(I,K)              ! P243.10
          TL(I,K) = T(I,K) - LCRCP*QCL(I,K) - LSRCP*QCF(I,K)  ! P243.9
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! Calculate buoyancy parameters BT and BQ.
!-----------------------------------------------------------------------
      CALL BOUY_TQ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,CF,Q,QCF,QCL,T,TL
     &,BT,BQ,BF,L_BL_LSPICE,LTIMER
     & )

!-----------------------------------------------------------------------
!! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!!     (P243a, routine SF_EXCH).
!!     Wind mixing "power" and some values required for other, later,
!!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------

      CALL SF_EXCH (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,LAND_INDEX,
     & TILE_INDEX,TILE_PTS,
     & BQ(1,1),BT(1,1),CANOPY,CATCH,DZSOIL(1),GC,HCONS,HO2R2_OROG,
     & ICE_FRACT,LYING_SNOW,PSTAR,P(1,1),QW(1,1),RADNET,RADNET_SNOW,
     & SIL_OROG_LAND,SMVCST,TILE_FRAC,TIMESTEP,
     & TL(1,1),TI,T_SOIL(1,1),TSNOW,TSTAR_TILE,TSTAR,
     & VSHR,Z0_TILE,Z0_SF_GB,Z1,Z1,
     & LAND_MASK,SU10,SV10,SQ1P5,ST1P5,SFME,LTIMER,L_Z0_OROG,Z0MSEA,
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_SNOW,CD,CH,CDR10M,CHR1P5M,
     & CHR1P5M_SICE,E_SEA,FME,FQW(1,1),FQW_TILE,FQW_ICE,
     & FTL(1,1),FTL_TILE,FTL_ICE,FRACA,H_BLEND_OROG,H_SEA,
     & Q1_SD,RESFS,RESFT,RIB,RIB_TILE,T1_SD,Z0M_EFF,
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,RHO_ARESIST,ARESIST,RESIST_B,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & RHO_CD_MODV1,RHOKH_TILE,RHOKH_SICE,RHOKM(1,1),RHOKPM,RHOKPM_SICE,
     & NRML
     & )

!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!      Interpolate BT and BQ to interface between layers.
!-----------------------------------------------------------------------

      CALL BTQ_INT (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,BQ,BT,BF,DZL,RDZ,QW,QCF,TL
     &,L_BL_LSPICE,LTIMER
     &  )

!-----------------------------------------------------------------------
!! 5.3  Calculate the diffusion coefficients Km and Kh.
!-----------------------------------------------------------------------

! Repeat of KMKH calculation, could be passed in from KMKH.

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)

! Calculate rho at K-1/2, from P243.111 :-
          RHO(I,K) =
     &     ( AKH(K) + BKH(K)*PSTAR(I) )    ! Pressure at K-1/2, P243.112
     &     /                               ! divided by ...
     &     ( R *                           ! R times ...
     &     ( TV(I,K-1)*WK + TV(I,K)*WKM1 ) ! TV at K-1/2, from P243.113
     &     )
        ENDDO
      ENDDO

      CALL KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & TIMESTEP,P,CCA,BT,BQ,BF,CF,DZL,DTRDZ,
     & RDZ,U_P,V_P,FTL,FQW,
     & RHO,Z0M_EFF,ZLB(1,0),H_BLEND_OROG,
     & QW,QCF,RHOKM,RHO_KM(1,2),RHOKH,TL,ZH,
     & CCB,CCT,L_MOM,
     & NRML,L_BL_LSPICE,LTIMER
     & )

!-----------------------------------------------------------------------
!! 5.4 Interpolate RHOKM's and CDR10M to uv points ready for the
!!     calculation of the explcit fluxes TAU_X and TAU_Y at levels
!!     above the surface.
!-----------------------------------------------------------------------

! RHOKM(*,1) contains duff data in halos. The P_TO_UV can interpolate
! this into the real data, so first we must update east/west halos

      CALL SWAPBOUNDS(RHOKM(P1,1),ROW_LENGTH,N_U_ROWS,1,0,1)
      CALL SWAPBOUNDS(RHOKM(1,2),ROW_LENGTH,
     &                U_FIELD/ROW_LENGTH,1,1,BL_LEVELS-1)

      DO K=1,BL_LEVELS

        CALL P_TO_UV (RHOKM(P1,K),RHOKM_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
      ENDIF

      ENDDO ! loop over bl_levels

        IF (SU10. OR. SV10)THEN

        CALL P_TO_UV (CDR10M(P1),CDR10M_UV(U1+ROW_LENGTH),P_POINTS,
     &     U_POINTS,ROW_LENGTH,N_P_ROWS)
!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------
        IF (attop) THEN
          DO I=U1,U1+ROW_LENGTH-1
            CDR10M_UV(I) = 1.0E30
          ENDDO
        ENDIF

        IF (atbase) THEN
          DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1
            CDR10M_UV(I) = 1.0E30
          ENDDO
        ENDIF

        ENDIF

!-----------------------------------------------------------------------
!! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------

      CALL EX_FLUX_TQ (
     &  P_POINTS,P_FIELD,P1,BL_LEVELS
     &, TL,QW,RDZ,FTL,FQW,RHOKH
     &, LTIMER
     &  )

!-----------------------------------------------------------------------
!! 5.6 Calculation of explicit fluxes of U and V.
!-----------------------------------------------------------------------

      CALL EX_FLUX_UV ( ! For U
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, U,U_0,RDZUV(1,2),RHOKM_UV,TAUX
     &, LTIMER
     &  )

      CALL EX_FLUX_UV ( ! For V
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, V,V_0,RDZUV(1,2),RHOKM_UV,TAUY
     &, LTIMER
     &  )

!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1 + (N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH -1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
      ENDIF
      ENDDO

!-----------------------------------------------------------------------
!! 6.  "Implicit" calculation of increments for TL and QW
!-----------------------------------------------------------------------

      CALL IM_CAL_TQ (
     &  P_FIELD,P1,P_POINTS,BL_LEVELS,LAND_FIELD,LAND_INDEX,NTYPE,
     &  TILE_INDEX,TILE_PTS,LAND_MASK,LTIMER,
     &  ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_SNOW,DTL_NT,DQW_NT,DTRDZ,
     &  ICE_FRACT,RDZ,RESFT,RHOKH(1,2),
     &  RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,TILE_FRAC,
     &  FQW,FQW_ICE,FQW_TILE,E_SEA,
     &  FTL,FTL_ICE,FTL_TILE,H_SEA,QW,TL
     &  )

!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
Cfpp$ Select(CONCUR)
        DO  I=P1,P1+P_POINTS-1
          FTL(I,K) = FTL(I,K)*CP
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
        FTL_ICE(I) = CP*FTL_ICE(I)
      ENDDO

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          FTL_TILE(L,N) = CP*FTL_TILE(L,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!   Sea-ice (P241, routine SICE_HTF).
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        IF ( .NOT.LAND_MASK(I) )
     &    SURF_HT_FLUX(I) = RADNET(I) - LS*FQW_ICE(I) - FTL_ICE(I)
      ENDDO

      CALL SICE_HTF(
     & ASHTF,DI,ICE_FRACT,SURF_HT_FLUX,TIMESTEP,
     & LAND_MASK,P_FIELD,P_POINTS,P1,TI,TSTAR,ASURF,
     & SEA_ICE_HTF,LTIMER
     &)

!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN
        DO L=LAND1,LAND1+LAND_PTS-1
          IF (T_SOIL(L,1).LT.0) THEN
            ERROR = 1
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'
            WRITE(6,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
            WRITE(6,*) 'LAND POINT ',L
          ENDIF
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!!   Diagnose the land surface temperature (previously in SOIL_HTF)
!-----------------------------------------------------------------------

      DO N=1,NTYPE
        DO L=LAND1,LAND1+LAND_PTS-1
          TSTAR_TILE(L,N) = T_SOIL(L,1)
        ENDDO
      ENDDO

      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          TSTAR_TILE(L,N) = TSTAR_TILE(L,N) +
     &    ( RADNET(I) - LC*FQW_TILE(L,N) - FTL_TILE(L,N) ) / ASHTF(I)
        ENDDO
      ENDDO

      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        TSTAR_TILE(L,N) = TSNOW(L) + ( RADNET_SNOW(I) - LS*FQW_TILE(L,N)
     &                                 - FTL_TILE(L,N) ) / ASHTF_SNOW(I)
      ENDDO

!-----------------------------------------------------------------------
!! 7.  Surface evaporation components and updating of surface
!!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
      CALL SF_EVAP (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,
     & LAND_INDEX,TILE_INDEX,TILE_PTS,SM_LEVELS,LTIMER,
     & ASHTF,ASHTF_SNOW,CANOPY,DTRDZ(1,1),FRACA,LYING_SNOW,RESFS,
     & RESFT,RHOKH_TILE,TILE_FRAC,SMC,WT_EXT,TIMESTEP,
     & FQW(1,1),FQW_TILE,FTL(1,1),FTL_TILE,QW(1,1),TL(1,1),TSTAR_TILE,
     & ECAN,ECAN_TILE,ESOIL,ESOIL_TILE,EXT
     & )

!-----------------------------------------------------------------------
!!     Surface melting of snow.
!!     Melting of sea-ice.
!-----------------------------------------------------------------------
      N = NTYPE
      CALL SF_MELT (
     & P_POINTS,P_FIELD,P1,LAND_FIELD,LAND_INDEX,
     & TILE_INDEX(1,N),TILE_PTS(N),LAND_MASK,LTIMER,SIMLT,SMLT,
     & ALPHA1(1,N),ALPHA1_SICE,ASHTF,ASHTF_SNOW,DTRDZ(1,1),ICE_FRACT,
     & LYING_SNOW,RHOKH_TILE(1,N),RHOKH_SICE,TILE_FRAC(1,N),TIMESTEP,
     & FQW(1,1),FQW_ICE,FQW_TILE(1,N),FTL(1,1),FTL_TILE(1,N),
     & QW(1,1),TL(1,1),TSTAR,TSTAR_TILE(1,N),TI,
     & EI,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT
     & )

!-----------------------------------------------------------------------
!!     Specific humidity and temperature at 1.5 metres.
!-----------------------------------------------------------------------
      CALL SCREEN_TQ (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,
     & LAND_INDEX,TILE_INDEX,TILE_PTS,LAND_MASK,
     & SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QW(1,1),RESFT,
     & TILE_FRAC,TL(1,1),TSTAR,TSTAR_TILE,
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z1,
     & Q1P5M,T1P5M
     & )

!7.1 Copy T and Q from workspace to INOUT space.

      DO K=1,BL_LEVELS
Cfpp$  Select(CONCUR)
        DO I=P1,P1+P_POINTS-1
          T(I,K)=TL(I,K)
          Q(I,K)=QW(I,K)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!     Gridbox-mean surface temperature and net surface heat fluxes
!-----------------------------------------------------------------------
      DO L=1,LAND_FIELD
        I = LAND_INDEX(L)
          TSTAR(I) = 0.
          SNOW_SURF_HTF(L) = 0.
          SOIL_SURF_HTF(L) = 0.
      ENDDO

      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          SOIL_SURF_HTF(L) = SOIL_SURF_HTF(L) + TILE_FRAC(L,N) *
     &                    (RADNET(I) - LC*FQW_TILE(L,N) - FTL_TILE(L,N))
          TSTAR(I) = TSTAR(I) + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
        ENDDO
      ENDDO

      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        SNOW_SURF_HTF(L) = TILE_FRAC(L,N) *
     &               (RADNET_SNOW(I) - LS*FQW_TILE(L,N) - FTL_TILE(L,N))
     &                     - LF*SNOWMELT(I)
        TSTAR(I) = TSTAR(I) + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
      ENDDO


      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        SURF_HT_FLUX(I) = SOIL_SURF_HTF(L) + SNOW_SURF_HTF(L) 
      ENDDO

      DO I=P1,P1+P_POINTS-1
        IF ( .NOT.LAND_MASK(I) )
     &    SURF_HT_FLUX(I) = RADNET(I) - LS*FQW_ICE(I) - FTL_ICE(I)
      ENDDO

!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          IF (TSTAR(I).LT.0) THEN
            ERROR = 1
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'
            WRITE(6,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',L
          ENDIF
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!! 8 "Implicit" calculation of increments for U and V.
!-----------------------------------------------------------------------

      CALL IM_CAL_UV (  ! For U
     & U_FIELD,U1
     &,U_POINTS,BL_LEVELS,ROW_LENGTH
     &,GAMMA
     &,RHOKM_UV(1,2)
     &,U,U_0,TIMESTEP
     &,RHOKM_UV(1,1),DU_NT,DU
     &,DTRDZ_UV,RDZUV(1,2),TAUX
     &,LTIMER
     &)

      CALL IM_CAL_UV (  ! For V
     & U_FIELD,U1
     &,U_POINTS,BL_LEVELS,ROW_LENGTH
     &,GAMMA
     &,RHOKM_UV(1,2)
     &,V,V_0,TIMESTEP
     &,RHOKM_UV(1,1),DV_NT,DV
     &,DTRDZ_UV,RDZUV(1,2),TAUY
     &,LTIMER
     & )

!----------------------------------------------------------------------
!! 8.1 Update U_V.
!----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          U(I,K) = U(I,K) + DU(I,K)
          V(I,K) = V(I,K) + DV(I,K)
        ENDDO
      ENDDO

! U component of 10m wind
      IF (SU10)THEN
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          U10M(I) = (U(I,1) -U_0(I))*CDR10M_UV(I) + U_0(I)
        ENDDO
      ENDIF

! V component of 10m wind
      IF (SV10)THEN
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          V10M(I) = (V(I,1) -V_0(I))*CDR10M_UV(I) + V_0(I)
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!! 9.  Calculate surface latent heat flux.
!-----------------------------------------------------------------------

      IF (SLH) THEN
        DO I=P1,P1+P_POINTS-1
          LATENT_HEAT(I) = LC*FQW(I,1) + LF*EI(I)
        ENDDO
      ENDIF

  999  CONTINUE  ! Branch for error exit.

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
