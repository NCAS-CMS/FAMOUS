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
!!!
!!!  Suitable for single column use - activate *IF definition IBM.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.4   10/09/97  New deck.   R.N.B.Smith
!!!   4.5   Jul. 98   Kill the IBM specific lines. (JCThil)
!!!
!!! Programming standard : unified model documentation paper No 3
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
     & P_FIELD,U_FIELD,N_TYPES,LAND_FIELD,
     & P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN soil/vegetation/land surface data :
     & LAND_MASK,GATHER,LAND_INDEX,
     & ST_LEVELS,SM_LEVELS,TILE_FRAC,HT_TILE,CANOPY,
     & CATCH,CATCH_TILE,HCON,
     & LYING_SNOW,RESIST,RESIST_TILE,ROOTD,ROOTD_TILE,
     & SMVCCL,SMVCST,SMVCWT,STHF,STHU,
     & VFRAC_TILE,Z0V,Z0V_TILE,SIL_OROG_LAND,L_Z0_OROG,HO2R2_OROG,
     & LAI_TILE,

! IN sea/sea-ice data :
     & DI,ICE_FRACT,U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,RAD_HR,RADHR_DIM1,
     & TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,

! INOUT data :
     & Q,GC,T,T_SOIL,TI,TSTAR,TSTAR_TILE,U,V,Z0MSEA,

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,E_SEA,EPOT,ETRAN,FQW,FQW_TILE,FSMC,FTL,FTL_TILE,
     & H_SEA,RHOKH,RHOKM_UV,
     & RIB_GB,RIB,SEA_ICE_HTF,SURF_HT_FLUX_GB,TAUX,TAUY,VSHR,ZHT,
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,    

! OUT diagnostic requiring STASH flags :
     & FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,
     & Q1P5M,T1P5M,U10M,V10M,

! (IN) STASH flags :-
     & SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,
     & NRML,

! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,RHO_KM,

! OUT data required elsewhere in UM system :
     & ECAN,EI,ES_GB,EXT,SNOWMELT,ZH,
     & GPP,NPP,RESP_P,
     & T1_SD,Q1_SD,ERROR,

! LOGICAL LTIMER
     & LTIMER
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid
!                                     (for dimensioning only).
     &,RADHR_DIM1                  ! IN Dimension of Radiative heating
!                                  !    rate (P_FIELD but used for
!                                  !    dynamic allocation)
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
!                                     (Checked for consistency with
!                                     P_FIELD and P_ROWS; there must
!                                     be 1 less UV than P row.)
     &,N_TYPES                     ! IN number of land tiles
     &,LAND_FIELD                  ! IN No.of land points in whole grid.
!                                     (Checked for consistency with
!                                     P_FIELD )
     &,P_ROWS                      ! IN No. of P-rows in whole grid
!                                     (for dimensioning only).
     &,FIRST_ROW                   ! IN First row of data to be treated,
!                                     referred to P-grid (must be > 1
!                                     since "polar" rows are never
!                                     treated).
     &,N_ROWS                      ! IN No. of rows of data to be
!                                     treated, referred to P-grid.
!                                     FIRST_ROW+N_ROWS-1 must be less
!                                     than P_ROWS, since "polar" rows
!                                     are never treated.
     &,ROW_LENGTH                  ! IN No. of points in one row.
!                                     (Checked for consistency with
!                                     P_FIELD and N_ROWS.)

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                    ! IN Total no. of vertical levels in
!                                       the model atmosphere.
      REAL
     & AK(P_LEVELS)                ! IN Hybrid 'A' for all levels.
     &,BK(P_LEVELS)                ! IN Hybrid 'B' for all levels.
     &,AKH(P_LEVELS+1)             ! IN Hybrid 'A' for layer interfaces.
     &,BKH(P_LEVELS+1)             ! IN Hybrid 'B' for layer interfaces.
     &,DELTA_AK(P_LEVELS)          ! IN Difference of hybrid 'A' across
!                                     layers (K-1/2 to K+1/2).
!                                     NB: Upper minus lower.
     &,DELTA_BK(P_LEVELS)          ! IN Difference of hybrid 'B' across
!                                     layers (K-1/2 to K+1/2).
!                                     NB: Upper minus lower.
     &,EXNER(P_FIELD,BL_LEVELS+1)  ! IN Exner function.  EXNER(,K) is
!                                     value for LOWER BOUNDARY of
!                                     level K.

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL
     & LAND_MASK(P_FIELD)          ! IN T if land, F elsewhere.
     &,L_Z0_OROG                   ! IN T to use orog.roughness
!                                     treatment in SF_EXCH
     &,L_RMBL                      ! IN T to use rapidly mixing boundary
!                                     scheme in IMPL_CAL
     &,L_BL_LSPICE                 ! IN True if 3A large-scale ppn
!                                       scheme is used.
     &,L_MOM                       ! IN Switch for convective momentum
!                                  !    transport.
     &,GATHER                      ! IN T if gather to sea-ice points
!                                     in SF_EXCH. Saves a lot of un-
!                                     necessary calculations if there
!                                     are relatively few sea-ice points

      INTEGER
     & LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

      INTEGER
     & ST_LEVELS                   ! IN No. of deep soil temp. levels
     &,SM_LEVELS                   ! IN No. of soil moisture levels

      REAL
     & CANOPY(LAND_FIELD)          ! IN Surface/canopy water (kg/m2)
     &,CATCH(LAND_FIELD)           ! IN Surface/canopy water capacity
!                                     (kg/m2).
     &,CATCH_TILE(LAND_FIELD,N_TYPES)
!                                    IN Surface/canopy water capacity
!                                     (kg per sq m).
     &,HCON(LAND_FIELD)            ! IN Soil thermal conductivity
!                                     (W/m/K).
     &,HT_TILE(LAND_FIELD,N_TYPES) ! IN Canopy height (m)
     &,LAI_TILE(LAND_FIELD,N_TYPES)! IN Leaf area index.
     &,LYING_SNOW(P_FIELD)         ! IN Lying snow (kg/sq m).
!                                     Must be global for coupled model,
!                                     ie dimension P_FIELD not
!                                     LAND_FIELD
     &,RESIST(LAND_FIELD)          ! IN "Stomatal" resistance to
!                                     evaporation (seconds per metre).
     &,RESIST_TILE(LAND_FIELD,N_TYPES)
!                                    IN "Stomatal" resistance to
!                                     evaporation (seconds per metre).
     &,ROOTD(LAND_FIELD)           ! IN Depth of active soil layer
!                                     ("root depth") (metres).
     &,ROOTD_TILE(LAND_FIELD,N_TYPES)
!                                    IN Depth of active soil layer
!                                     ("root depth") (metres).
     &,SMVCCL(LAND_FIELD)          ! IN Critical volumetric SMC (m3/m3
!                                     of soil).
     &,SMVCST(LAND_FIELD)          ! IN Volumetric saturation point
!                                     (m3/m3 of soil).
     &,SMVCWT(LAND_FIELD)          ! IN Volumetric wilting point (m3/m3
!                                     of soil).
     &,STHF(LAND_FIELD,SM_LEVELS)  ! IN Frozen soil moisture content of
!                                     each layer as a fraction of
!                                     saturation.
     &,STHU(LAND_FIELD,SM_LEVELS)  ! IN Unfrozen soil moisture content
!                                     of each layer as a fraction of
!                                     saturation.
     &,TILE_FRAC(P_FIELD,N_TYPES)  ! IN fractional coverage for each
!                                     surface tile
     &,VFRAC_TILE(LAND_FIELD,N_TYPES)
!                                  ! IN Vegetation fraction.
     &,Z0V(P_FIELD)                ! IN Vegetative roughness length (m).
!                                     NB:UM uses same storage for Z0MSEA
!                                     so for sea points this is INOUT.
     &,Z0V_TILE(P_FIELD,N_TYPES)   ! IN Vegetative roughness length (m)
!                                     for surface tile
     &,SIL_OROG_LAND(LAND_FIELD)   ! IN Silhouette area of unresolved
!                                     orography per unit horizontal area
!                                     on land points only.
     &,HO2R2_OROG(LAND_FIELD)      ! IN Standard Deviation of orography.
!                                     equivilent to peak to trough
!                                     height of unresolved orography
!                                     devided by 2SQRT(2) on land
!                                     points only (m)

! (d) Sea/sea-ice data.

      REAL
     & DI(P_FIELD)                 ! IN "Equivalent thickness" of
!                                     sea-ice(m).
     &,ICE_FRACT(P_FIELD)          ! IN Fraction of gridbox covered by
!                                     sea-ice (decimal fraction).
     &,U_0(U_FIELD)                ! IN W'ly component of surface
!                                     current (m/s).
     &,V_0(U_FIELD)                ! IN S'ly component of surface
!                                     current (m/s).

! (e) Cloud data.

      REAL
     & CF(P_FIELD,BL_LEVELS)       ! IN Cloud fraction (decimal).
     &,QCF(P_FIELD,BL_LEVELS)      ! IN Cloud ice (kg per kg air)
     &,QCL(P_FIELD,BL_LEVELS)      ! IN Cloud liquid water (kg
!                                     per kg air).
     &,CCA(P_FIELD)                ! IN Convective Cloud Amount
!                                     (decimal)

      INTEGER
     & CCB(P_FIELD)                ! IN Convective Cloud Base
     &,CCT(P_FIELD)                ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL
     & CO2_MMR                     ! IN CO2 Mass Mixing Ratio
     &,PHOTOSYNTH_ACT_RAD(P_FIELD) ! IN Net downward shortwave radiation
!                                     in band 1 (w/m2).
     &,PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).
     &,RAD_HR(RADHR_DIM1,BL_LEVELS)! IN Radiative heating rate (K/s).
     &,RADNET(P_FIELD)             ! IN Surface net radiation (W/sq m,
!                                     positive downwards).
     &,TIMESTEP                    ! IN Timestep (seconds).

      LOGICAL LTIMER               ! Logical switch for TIMER diags

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
     & GC(LAND_FIELD,N_TYPES)      ! INOUT "Stomatal" conductance to
!                                      evaporation (m/s).
     &,Q(P_FIELD,BL_LEVELS)        ! INOUT Input:specific humidity
!                                      ( kg/kg air).
!                                      Output:total water content
!                                      (Q)(kg/Kg air).
     &,T(P_FIELD,BL_LEVELS)        ! INOUT Input:atmospheric temp(K)
!                                      Output:liquid/frozen water
!                                      temperature (TL) (K)
     &,T_SOIL(LAND_FIELD,SM_LEVELS)! INOUT Soil temperatures (K).
     &,TI(P_FIELD)                 ! INOUT Sea-ice surface layer
!                                      temperature (K).
     &,TSTAR(P_FIELD)              ! INOUT Surface temperature (K).
     &,TSTAR_TILE(P_FIELD,N_TYPES) ! INOUT Surface tile temperature
     &,U(U_FIELD,BL_LEVELS)        ! INOUT W'ly wind component (m/s)
     &,V(U_FIELD,BL_LEVELS)        ! INOUT S'ly wind component (m/s)
     &,Z0MSEA(P_FIELD)             ! INOUT Sea-surface roughness
!                                      length for momentum (m).
!                                      NB: same storage is used
!                                      for Z0V, so the intent is
!                                      IN for land points.

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL
     & CD(P_FIELD)                 ! OUT Turbulent surface exchange
!                                     (bulk transfer) coefficient for
!                                     momentum.
     &,CH(P_FIELD)                 ! OUT Turbulent surface exchange
!                                     (bulk transfer) coefficient for
!                                     heat and/or moisture.
     &,E_SEA(P_FIELD)              ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
     &,EPOT_TILE(P_FIELD,N_TYPES)  ! WORK potential evaporation
!                                     over tile (kg/m2/s).
     &,EPOT(P_FIELD)               ! OUT potential evaporation (kg/m2/s)
     &,FQW(P_FIELD,BL_LEVELS)      ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FQW_TILE(P_FIELD,N_TYPES)   ! OUT surface tile moisture flux
     &,FSMC(LAND_FIELD)            ! OUT soil moisture availability.
     &,FSMC_TILE(LAND_FIELD,N_TYPES)
                                   ! WORK soil moisture availability
!                                     over tile.
     &,FTL(P_FIELD,BL_LEVELS)      ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,FTL_TILE(P_FIELD,N_TYPES)   ! OUT surface tile heat flux
     &,H_SEA(P_FIELD)              ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
     &,RHOKH(P_FIELD,BL_LEVELS)    ! OUT Exchange coeffs for moisture.
     &,RHOKM_UV(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                                     momentum (on UV-grid, with 1st
!                                     and last rows undefined (or, at
!                                     present, set to "missing data"))
     &,RIB(P_FIELD,N_TYPES)        ! OUT Tile bulk Richardson number for
!                                     lowest layer.
     &,RIB_GB(P_FIELD)             ! OUT Mean bulk Richardson number for
!                                     lowest layer.
     &,SEA_ICE_HTF(P_FIELD)        ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).
     &,SURF_HT_FLUX_GB(P_FIELD)    ! OUT Net downward heat flux at
!                                     surface over land or sea-ice
!                                     fraction of gridbox (W/m2).
     &,TAUX(U_FIELD,BL_LEVELS)     ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(U_FIELD,BL_LEVELS)     ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
     &,VSHR(P_FIELD)               ! OUT Magnitude of surface-to-lowest
!                                     atm level wind shear (m per s).
     &,ZHT(P_FIELD)                ! OUT Height below which there may be
!                                  !     turbulent mixing (m).
     &,BL_TYPE_1(P_FIELD)          ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_2(P_FIELD)          ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_3(P_FIELD)          ! OUT Indicator set to 1.0 if well 
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_4(P_FIELD)          ! OUT Indicator set to 1.0 if 
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_5(P_FIELD)          ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
     &,BL_TYPE_6(P_FIELD)          ! OUT Indicator set to 1.0 if a 
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,RHO_CD_MODV1(P_FIELD)       ! OUT Surface air density * drag coef
!                                     *mod(v1 - v0) before interpolation
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT Air density * turbulent mixing
!                                     coefficient for momentum before
!                                     interpolation.
     &,RHO_ARESIST(P_FIELD)        ! OUT RHOSTAR*CD_STD*VSHR for SULPHUR
!                                     cycle
     &,ARESIST(P_FIELD)            ! OUT 1/(CD_STD*VSHR) for Sulphur
!                                     cycle
     &,RESIST_B(P_FIELD)           ! OUT (1/CH-1/(CD_STD)/VSHR for
!                                     Sulphur cycle

      INTEGER
     & NRML(P_FIELD)               ! OUT Number of model layers in the
!                                     Rapidly Mixing Layer; diagnosed
!                                     in SF_EXCH and KMKH and used in
!                                     IMPL_CAL, SF_EVAP and TR_MIX.

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL
     & FME(P_FIELD)              ! OUT Wind mixing "power" (W per sq m).
     &,SICE_MLT_HTF(P_FIELD)     ! OUT Heat flux due to melting of sea-
!                                   ice (Watts per sq metre).
     &,SNOMLT_SURF_HTF(P_FIELD)  ! OUT Heat flux required for surface
!                                   melting of snow (W/m2).
     &,LATENT_HEAT(P_FIELD)      ! OUT Surface latent heat flux, +ve
!                                   upwards (Watts per sq m).
     &,Q1P5M(P_FIELD)            ! OUT Q at 1.5 m (kg water per kg air).
     &,T1P5M(P_FIELD)            ! OUT T at 1.5 m (K).
     &,U10M(U_FIELD)             ! OUT U at 10 m (m per s).
     &,V10M(U_FIELD)             ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL
     & EI(P_FIELD)               ! OUT Sublimation from lying snow or
!                                   sea-ice (kg/m2/s).
     &,ECAN(P_FIELD)             ! OUT Gridbox mean evaporation from
!                                   canopy/surface store (kg/m2/s).
!                                   Zero over sea.
     &,ES_GB(P_FIELD)            ! OUT Surface evapotranspiration
!                                   through a resistance which is not
!                                   entirely aerodynamic i.e. "soil
!                                   evaporation".  Always non-negative.
!                                   (kg/m2/s).
     &,ETRAN(P_FIELD,N_TYPES)    ! OUT Transpiration (kg/m2/s).
     &,EXT(LAND_FIELD,SM_LEVELS) ! OUT Extraction of water from each
!                                   soil layer (kg/m2/s).
     &,GPP(LAND_FIELD,N_TYPES)   ! OUT Gross primary productivity
!                                   (kg C/m2/s).
     &,NPP(LAND_FIELD,N_TYPES)   ! OUT Net primary productivity
!                                   (kg C/m2/s).
     &,RESP_P(LAND_FIELD,N_TYPES)! OUT Plant respiration (kg C/m2/s).
     &,SNOWMELT(P_FIELD)         ! OUT Snowmelt (kg/m2/s).
     &,ZH(P_FIELD)               ! INOUT Height above surface of top of
!                                   boundary layer (metres).
     &,T1_SD(P_FIELD)            ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 temperature;
!                                   for use in initiating convection.
     &,Q1_SD(P_FIELD)            ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 humidity;
!                                   for use in initiating convection.
      INTEGER
     & ERROR          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL Z,HEAT_CON,SMC_ROOT,SF_EXCH,BOUY_TQ,BTQ_INT,
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,
     & IM_CAL_UV
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
     & A_DQSDT(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_DQSDTM(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,ALPHA1(P_FIELD,N_TYPES)  ! Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
     &,ALPHA1_GB(P_FIELD)       ! Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces
     &,A_QS(P_FIELD,BL_LEVELS)  ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_QSM(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,ASHTF(P_FIELD)           ! Coefficient to calculate surface
!                                 heat flux into soil or sea-ice.
     &,ASURF(P_FIELD)           ! Reciprocal areal heat capacity
!                                 of soil layer or sea-ice
!                                 surface layer (K m**2 / J).
     &,BQ(P_FIELD,BL_LEVELS)    ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BQ_CLD(P_FIELD,BL_LEVELS)! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BQM(P_FIELD,BL_LEVELS)   ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BQM_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,BT(P_FIELD,BL_LEVELS)    ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BT_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BTM(P_FIELD,BL_LEVELS)   ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BTM_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,DB(P_FIELD,2:BL_LEVELS)
!                               ! Buoyancy jump across layer interface.
     &,DELTAP(P_FIELD,BL_LEVELS)! Difference in pressure between levels
     &,DELTAP_UV(P_FIELD,BL_LEVELS)
!                                 Difference in pressure between levels
!                                 on UV points
     &,DQSDT(P_FIELD,BL_LEVELS) ! Derivative of q_SAT w.r.t. T
     &,DQW_1(P_FIELD)           ! Increment for QW(,1).
     &,DTRDZ(P_FIELD,BL_LEVELS) ! -g.dt/dp for model layers.
     &,DTRDZ_UV(U_FIELD,BL_LEVELS)
!                                 -g.dt/dp for model wind layers.
     &,DTRDZ_RML(P_FIELD)       ! -g.dt/dp for the rapidly
!                                 mixing layer.
     &,DZL(P_FIELD,BL_LEVELS)   ! DZL(,K) is depth in m of layer
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,DU(U_FIELD,BL_LEVELS)    ! BL increment to u wind foeld
     &,DV(U_FIELD,BL_LEVELS)    ! BL increment to v wind foeld
     &,DU_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to u wind field
     &,DV_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to v wind field
     &,DTL_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to TL field
     &,DQW_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to QW field
     &,ES(P_FIELD,N_TYPES)      ! Surface evapotranspiration
!                                 through a resistance which is not
!                                 entirely aerodynamic i.e. "soil
!                                 evaporation".  Always non-negative.
!                                 (kg/m2/s).
     &,ESOIL(P_FIELD,N_TYPES)   ! Evaporation from bare soil (kg/m2
     &,FB_SURF(P_FIELD)         ! Surface flux buoyancy over density
!                               ! (m^2/s^3)
!
     &,FRACA(P_FIELD,N_TYPES)   ! Fraction of surface moisture flux
!                                 with only aerodynamic resistance.
     &,F_SE(P_FIELD,N_TYPES)    ! Fraction of the evapotranspiration
!                                 which is bare soil evaporation.
     &,GRAD_Q_ADJ(P_FIELD)      ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,GRAD_T_ADJ(P_FIELD)      ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,HEAT_BLEND_FACTOR(P_FIELD,N_TYPES)
!                                 Blending factor used as part of
!                                 tile scheme
     &,HCONS(LAND_FIELD)        ! Soil thermal conductivity includi
!                                 the effects of water and ice (W/m
     &,QW(P_FIELD,BL_LEVELS)    ! Total water content, but
!                                 replaced by specific humidity
!                                 in LS_CLD.
     &,P(P_FIELD,BL_LEVELS)     ! P(*,K) is pressure at full level k.
     &,P_HALF(P_FIELD,BL_LEVELS)! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,Z_FULL(P_FIELD,BL_LEVELS)! Z_FULL(*,K) is height of full level k.
     &,Z_HALF(P_FIELD,BL_LEVELS)! Z_HALF(*,K) is height of half level
!                               ! k-1/2.
     &,Z_UV(P_FIELD,BL_LEVELS)  ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(P_FIELD,BL_LEVELS)  ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,RDZ(P_FIELD,BL_LEVELS)   ! RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,RDZUV(U_FIELD,BL_LEVELS) !  RDZ (K > 1) on UV-grid.
!                                  Comments as per RHOKM (RDZUV).
     &,RESFS(P_FIELD,N_TYPES)   ! Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
     &,RESFT_TILE(P_FIELD,N_TYPES)
!                                 Total resistance factor for tile
!                                 FRACA+(1-FRACA)*RESFS.
     &,RESFT(P_FIELD)           ! Mean total resistance factor
!                                 FRACA+(1-FRACA)*RESFS.
     &,RHO_FULL(P_FIELD,BL_LEVELS)
!                               ! RHO_FULL(*,K) is the density at full
!                               ! model level k.
     &,RHO_HALF(P_FIELD,BL_LEVELS)
!                               ! RHO_HALF(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_UV(P_FIELD,BL_LEVELS)
!                               ! RHO_UV(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_TQ(P_FIELD,BL_LEVELS)
!                               ! RHO_TQ(*,K) is the density at half
!                               ! level k+1/2.
     &,RHOKE(P_FIELD,N_TYPES)   ! Surface exchange coefficient for FQW
     &,RHOKH_TILE(P_FIELD,N_TYPES)
!                                 Tile surface exchange coefficients
!                                 for heat
     &,RHOKHZ(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
     &,RHOKH_TOP(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of heat and
!                               ! moisture.
     &,RHOKM(P_FIELD,BL_LEVELS) ! Turbulent mixing coefficient for
!                                 momentum on P-grid.
     &,RHOKMZ(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing
!                                 coefficient for momentum.
     &,RHOKM_TOP(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of momentum.
     &,RHOKPM(P_FIELD)          ! Surface exchange coefficient.
     &,RHOKPM_POT(P_FIELD)      ! WORK Surface exchange coeff. for
!                                 potential evaporation.
     &,RHOKPM_POT_TILE(P_FIELD,N_TYPES)
                                ! WORK Tile surface exchange coeff.
!                                 for potential evaporaiotn.
     &,RHOKPM_TILE(P_FIELD,N_TYPES)
!                                 Surface exchange coefficient.
     &,SMC(LAND_FIELD,N_TYPES)  ! Soil moisture content in root depth
!                                  (kg/m2).
     &,SURF_HT_FLUX(P_FIELD,N_TYPES)
!                                 Net downward heat flux at surface
!                                 over land or sea-ice fraction of
!                                 gridbox (W/m2).
     &,TL(P_FIELD,BL_LEVELS)    ! Ice/liquid water temperature,
!                                 but replaced by T in LS_CLD.
     &,TV(P_FIELD,BL_LEVELS)    ! Virtual temp
     &,TV1_SD(P_FIELD)          ! Standard deviation of turbulent
!                               ! fluctuations of surface layer
!                               ! virtual temperature (K).
     &,U_P(P_FIELD,BL_LEVELS)   ! U on P-grid.
     &,U_0_P(P_FIELD)           ! U_0 on P-grid.
     &,U_S(P_FIELD)             ! Surface friction velocity (m/s)
     &,V_P(P_FIELD,BL_LEVELS)   ! V on P-grid.
     &,V_0_P(P_FIELD)           ! V_0 on P-grid.
     &,V_ROOT(LAND_FIELD,N_TYPES)! Volumetric soil moisture
!                                  concentration in the rootzone
!                                  (m3 H2O/m3 soil).
     &,V_SOIL(LAND_FIELD)       ! Volumetric soil moisture
!                                 concentration in the top
!                                 soil layer (m3 H2O/m3 soil).
     &,WIND_BLEND_FACTOR(P_FIELD,N_TYPES)
!                                 Blending factor used as part of
!                                 tile scheme
     &,WT_EXT(LAND_FIELD,SM_LEVELS)
!                                 Fraction of transpiration which is
!                                 extracted from each soil layer.
     &,ZLB(P_FIELD,0:BL_LEVELS) ! ZLB(,K) is the height of the
!                                 upper boundary of layer K
!                                 ( = 0.0 for "K=0").
       REAL
     & Z0H(P_FIELD,N_TYPES)     ! Roughness length for heat and
!                                 moisture.
     &,Z0M(P_FIELD,N_TYPES)     ! Roughness length for momentum.
     &,Z1(P_FIELD)              ! Height of lowest level (i.e.
!                                 height of middle of lowest
!                                 layer).
     &,H_BLEND_OROG(P_FIELD)    ! Blending height used as part of
!                                 effective roughness scheme
     &,H_BLEND(P_FIELD)         ! Blending height for tiles
     &,Z0M_EFF_GB(P_FIELD)      ! Effective grid-box roughness
!                                 length for momentum
     &,Z0M_EFF(P_FIELD,N_TYPES) ! Effective tile roughness length
!                                 for momentum
     &,Z_LCL(P_FIELD)           ! Height of lifting condensation level.
     &,CDR10M(P_FIELD)          ! Ratio of CD's reqd for calculation
!                                 of 10 m wind. On P-grid
     &,CDR10M_UV(U_FIELD)       ! Ratio of CD's reqd for calculation
!                                 of 10 m wind. On UV-grid; comments as
!                                 per RHOKM.
     &,CER1P5M(P_FIELD)         ! Ratio of coefficients reqd for
!                                 calculation of 1.5 m Q.
     &,CHR1P5M(P_FIELD)         ! Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
!
!   Variables for Vegetation Thermal Canopy
!
      REAL
     + CANCAP(P_FIELD,N_TYPES)    ! WORK Volumetric heat capacity of
!                                 !      vegetation canopy (J/Kg/m3).
     +,RADNET_C(P_FIELD,N_TYPES)  ! WORK Adjusted net radiation for
!                                 !      vegetation canopy over land
!                                 !      (W/m2).

      INTEGER
     & F_TYPE(LAND_FIELD,N_TYPES)! Plant functional type:
                                 !       1 - Broadleaf Tree
                                 !       2 - Needleleaf Tree
                                 !       3 - C3 Grass
                                 !       4 - C4 Grass
      INTEGER
     & NTML(P_FIELD)            ! Number of model levels in the
!                                 turbulently mixed layer.
     &,NTDSC(P_FIELD)           ! Top level for turbulent mixing in
!                               ! cloud layer.
      LOGICAL
     & CUMULUS(P_FIELD)         ! Logical switch for cumulus in the b.l.
     &,UNSTABLE(P_FIELD)        ! Logical switch for unstable 
!                                 surface layer.
     &,DSC(P_FIELD)             ! Flag set if decoupled stratocumulus
!                               ! layer found.

!  Local scalars :-

      REAL
     & WK         ! LOCAL 0.5 * DZL(I,K) * RDZ(I,K)
     &,WKM1       ! LOCAL 0.5 * DZL(I,K-1) * RDZ(I,K)

      INTEGER
     & I,J,L      ! LOCAL Loop counter (horizontal field index).
     &,ITILE      ! LOCAL Loopy counter (tile index).
     &,N          ! LOCAL Loop counter (soil levels)
     &,K          ! LOCAL Loop counter (vertical level index).
     &,N_P_ROWS   ! LOCAL No of P-rows being processed.
     &,N_U_ROWS   ! LOCAL No of UV-rows being processed.
     &,P_POINTS   ! LOCAL No of P-points being processed.
     &,P1         ! LOCAL First P-point to be processed.
     &,LAND1      ! LOCAL First land-point to be processed.
!                         1 <= LAND1 <= LAND_FIELD
     &,LAND_PTS   ! LOCAL No of land points being processed.
     &,U_POINTS   ! LOCAL No of UV-points being processed.
     &,U1         ! LOCAL First UV-point to be processed.

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF
      ERROR = 0
C-----------------------------------------------------------------------
C Initialise RADNET_C to be the same as RADNET over all points
C-----------------------------------------------------------------------
      DO ITILE=1,N_TYPES
        DO I=1,P_FIELD
          RADNET_C(I,ITILE) = RADNET(I)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
!!    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
!!    most UV-rows are used only for interpolation and are not updated.
!-----------------------------------------------------------------------

      IF ( BL_LEVELS.LT.1 .OR. ST_LEVELS.LT.1 .OR. SM_LEVELS.LT.1
     & .OR. P_ROWS.LT.3 ) THEN
        ERROR = 1
        GOTO999
      ELSEIF ( U_FIELD .NE. P_ROWS*ROW_LENGTH ) THEN
        ERROR = 2
        GOTO999
      ELSEIF ( P_FIELD .NE. P_ROWS*ROW_LENGTH ) THEN
        ERROR = 3
        GOTO999
      ELSEIF ( FIRST_ROW.LE.1 .OR. FIRST_ROW.GE.P_ROWS ) THEN
        ERROR = 4
        GOTO999
      ELSEIF ( N_ROWS.LE.0 ) THEN
        ERROR = 5
        GOTO999
      ELSEIF ( (FIRST_ROW+N_ROWS-1) .GT. P_ROWS ) THEN
        ERROR = 6
        GOTO999
      ELSEIF ( LAND_FIELD.GT.P_FIELD ) THEN
        ERROR = 7
        GOTO999
      ENDIF

!-----------------------------------------------------------------------
!!    Set pointers, etc.
!-----------------------------------------------------------------------

      N_P_ROWS = N_ROWS
      N_U_ROWS = N_ROWS + 1

      P_POINTS = N_P_ROWS * ROW_LENGTH
      U_POINTS = N_U_ROWS * ROW_LENGTH

      P1 = 1 + (FIRST_ROW-1)*ROW_LENGTH
      U1 = 1 + (FIRST_ROW-2)*ROW_LENGTH

!-----------------------------------------------------------------------
!!    Set compressed land point pointers.
!-----------------------------------------------------------------------

      LAND1=0
      DO I=1,P1+P_POINTS-1
        IF (LAND_INDEX(I).GE.P1) THEN
          LAND1 = I
          GOTO2
        ENDIF
      ENDDO
   2  CONTINUE

      LAND_PTS=0
      DO I=P1,P1+P_POINTS-1
        IF (LAND_MASK(I)) LAND_PTS = LAND_PTS + 1
      ENDDO


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
      ENDDO
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          Z_FULL(I,K) = ZLB(I,K) - 0.5 * DZL(I,K)
          Z_HALF(I,K) = ZLB(I,K-1)
          Z_UV(I,K) = ZLB(I,K-1)
          Z_TQ(I,K) = ZLB(I,K)
       ENDDO
      ENDDO
      DO K=1,BL_LEVELS

        CALL UV_TO_P(U(U1,K),U_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V(U1,K),V_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)


! du_nt 'borrowed to store dzl on uv grid
        CALL P_TO_UV (DZL(P1,K),DU_NT(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

      ENDDO

        CALL UV_TO_P(U_0(U1),U_0_P(P1),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V_0(U1),V_0_P(P1),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)


! set pressure array.
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          P(I,K) = AK(K) + BK(K)*PSTAR(I)
          P_HALF(I,K) = AKH(K) + BKH(K)*PSTAR(I)

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


!!----------------------------------------------------------------------
!! 2. Diagnose the plant functional types at each location.
!! Assume : Broadleaf Trees if rootdepth > 0.8m
!          C3 Grass        if rootdepth < 0.8m
!-----------------------------------------------------------------------
      DO ITILE=1,N_TYPES
        DO L=1,LAND_FIELD
          IF (ROOTD_TILE(L,ITILE).GT.0.8) THEN
            F_TYPE(L,ITILE)=1
          ELSE
            F_TYPE(L,ITILE)=3
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.
!-----------------------------------------------------------------------

      IF(LAND_FIELD.GT.0) THEN    ! Omit if no land points
        CALL HEAT_CON (LAND_FIELD,HCON,STHU,STHF,SMVCST,HCONS,LTIMER)

!-----------------------------------------------------------------------
! Calculate the soil moisture in the root zone.
!-----------------------------------------------------------------------

        DO ITILE=1,N_TYPES
          CALL SMC_ROOT (LAND_FIELD,SM_LEVELS,F_TYPE(1,ITILE),DZSOIL,
     &                   ROOTD_TILE(1,ITILE),
     &                   STHU,VFRAC_TILE(1,ITILE),SMVCST,SMVCWT,
     &                   SMC(1,ITILE),V_ROOT(1,ITILE),V_SOIL,WT_EXT,
     &                   LTIMER)
        ENDDO

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
!! 3.  Calls to SICE_HTF and SOIL_HTF now after IMPL_CAL
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!!     (P243a, routine SF_EXCH).
!!     Wind mixing "power" and some values required for other, later,
!!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------


! Set lots of things to zero

      DO ITILE=1,N_TYPES
        DO I=1,P_FIELD
          ETRAN(I,ITILE)=0.0
          ALPHA1(I,ITILE)=0.0
          FQW_TILE(I,ITILE)=0.0
          FTL_TILE(I,ITILE)=0.0
          FRACA(I,ITILE)=0.0
          RESFS(I,ITILE)=0.0
          RESFT_TILE(I,ITILE)=0.0
          RHOKH_TILE(I,ITILE)=0.0
          RHOKPM_TILE(I,ITILE)=0.0
          Z0H(I,ITILE)=0.0
          Z0M_EFF(I,ITILE)=0.0
          WIND_BLEND_FACTOR(I,ITILE)=0.0
          HEAT_BLEND_FACTOR(I,ITILE)=0.0

          IF(.NOT. LAND_MASK(I)) TILE_FRAC(I,ITILE)=0.0

          tstar_tile(i,itile)=tstar(i)  ! temporary for single tile only

        ENDDO
      ENDDO

      DO N=1,SM_LEVELS
        DO I=LAND1,LAND1+LAND_PTS-1
          EXT(I,N)=0.0
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
!         IF(.NOT. LAND_MASK(I)) TILE_FRAC(I,1)=1.0
         TILE_FRAC(I,1)=1.0  ! hard wired for single tile only

         SURF_HT_FLUX_GB(I)=0.0
         ES_GB(I)=0.0

      ENDDO



      CALL SF_EXCH (
     & P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,N_TYPES,
     & P1,LAND1,
     & LAND_INDEX,GATHER,
     & P(1,1),TILE_FRAC,
     & CANOPY,CATCH_TILE,CO2_MMR,
     & SM_LEVELS,DZSOIL,HCONS,F_TYPE,
     & HT_TILE,LAI_TILE,PHOTOSYNTH_ACT_RAD,GPP,NPP,RESP_P,
     & ICE_FRACT,LAND_MASK,LYING_SNOW,PSTAR,Q(1,1),
     & QCF(1,1),QCL(1,1),RADNET_C,GC,RESIST_TILE,
     & ROOTD_TILE,SMC,SMVCCL,SMVCWT,
     & T(1,1),TIMESTEP,TI,T_SOIL(1,1),TSTAR,
     & TSTAR_TILE,U_P(1,1),V_P(1,1),U_0_P,V_0_P,
     & V_ROOT,V_SOIL,VFRAC_TILE,
     & Z0V,Z0V_TILE,SIL_OROG_LAND,HO2R2_OROG,ZH,
     & Z1,Z1,CANCAP,Z0MSEA,ALPHA1_GB,ALPHA1,ASHTF,              
     & BQ(1,1),BT(1,1),CD,CH,
     & FQW_TILE,FQW(1,1),FTL_TILE,FTL(1,1),
     & EPOT_TILE,EPOT,FSMC_TILE,FSMC,
     & E_SEA,H_SEA,FRACA,RESFS,F_SE,
     & RESFT_TILE,RESFT,RHOKE,RHOKH_TILE,
     & RHOKH,RHOKM,RHOKPM_TILE,RHOKPM,RHOKPM_POT_TILE,RHOKPM_POT,
     & RIB_GB,RIB,TL(1,1),VSHR,Z0H,Z0M,Z0M_EFF,Z0M_EFF_GB,
     & H_BLEND_OROG,H_BLEND,T1_SD,Q1_SD,TV1_SD,U_S,FB_SURF,
     & RHO_CD_MODV1,WIND_BLEND_FACTOR,HEAT_BLEND_FACTOR,
     & CDR10M,CHR1P5M,CER1P5M,FME,
     & SU10,SV10,SQ1P5,ST1P5,SFME,
     & RHO_ARESIST,ARESIST,RESIST_B,NRML,
     & L_Z0_OROG,L_RMBL,LTIMER
     &)


!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 5.1  Calculate bouyancy parameters BT and BQ.
!-----------------------------------------------------------------------

      CALL BOUY_TQ (
     & P_FIELD,P1
     &,P_POINTS,BL_LEVELS
     &,P,T,Q,QCF,QCL
     &,BT,BQ,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,LTIMER
     &  )


!-----------------------------------------------------------------------
!! 5.2  Interpolate BT and BQ to half levels.
!-----------------------------------------------------------------------

      CALL BTQ_INT (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,DZL,RDZ,BQ,BT,BQ_CLD,BT_CLD,A_QS,A_DQSDT
     &,BQM,BTM,BQM_CLD,BTM_CLD,A_QSM,A_DQSDTM
     &,LTIMER
     &  )


!-----------------------------------------------------------------------
!! 5.3  Calculate the diffusion coefficients Km and Kh.
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHO_FULL(I,K) =
     &     ( AK(K) + BK(K)*PSTAR(I) )      ! Pressure at K
     &     /                               ! divided by ...
     &     ( R * TV(I,K) )                 ! R times TV at K
        ENDDO
      ENDDO
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)
          RHO_HALF(I,K) = WK*RHO_FULL(I,K-1) + WKM1*RHO_FULL(I,K)
        ENDDO
      ENDDO
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHO_UV(I,K) = RHO_HALF(I,K)
        ENDDO
      ENDDO
      DO K=1,BL_LEVELS-1
        DO I=P1,P1+P_POINTS-1
          RHO_TQ(I,K) = RHO_HALF(I,K+1)
        ENDDO
      ENDDO
      DO I=P1,P1+P_POINTS-1
        RHO_HALF(I,1) = RHO_FULL(I,1)
        RHO_UV(I,1) = RHO_FULL(I,1)
        RHO_TQ(I,BL_LEVELS) = RHO_FULL(I,BL_LEVELS)
      ENDDO

      CALL KMKHZ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & P,P_HALF,T,Q,QCL,QCF,BT,BQ,CF,DZL,
     & RDZ,DELTAP,FTL,FQW,
     & Z0M_EFF_GB,Z_FULL,Z_HALF,Z_UV,Z_TQ,U_S,FB_SURF,
     & QW,RHOKMZ(1,2),DB(1,2),RHOKHZ(1,2),TL,ZH,TV1_SD,T1_SD,Q1_SD,
     & NTML,GRAD_T_ADJ,GRAD_Q_ADJ,
     & BTM,BQM,DQSDT,BTM_CLD,BQM_CLD,A_QSM,A_DQSDTM,RHO_TQ,RHO_UV,
     & RAD_HR,RADHR_DIM1,CUMULUS,Z_LCL,RHOKM_TOP(1,2),RHOKH_TOP(1,2),
     & ZHT,BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,
     & UNSTABLE,NTDSC,DSC,
     & LTIMER
     & )

      CALL EX_COEF (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,CCB,CCT,NTML,L_MOM
     &,CCA,DZL,RDZ,DB(1,2),U_P,V_P
     &,RHO_HALF,ZH,Z_HALF,Z0M,H_BLEND_OROG
     &,CUMULUS,Z_LCL
     &,RHOKM,RHOKH
     &,LTIMER
     & )

      CALL KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,RHOKM,RHO_KM(1,2),RHOKH
     &,RHOKMZ(1,2),RHOKHZ(1,2)
     &,NTML,CUMULUS,RHOKM_TOP(1,2),RHOKH_TOP(1,2)
     &,UNSTABLE,NTDSC,DSC
     &,LTIMER
     & )

!
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

! CDR10M contains incorrect data in halos. The P_TO_UV can interpolate
! this into the real data, so first we must update east/west halos.
      CALL SWAPBOUNDS(CDR10M(P1),ROW_LENGTH,N_U_ROWS,1,0,1)


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

      IF (L_BL_LSPICE) THEN

        DO K = 1,BL_LEVELS
          DO I = P1,P1+P_POINTS-1
            QW(I,K) = Q(I,K) + QCL(I,K)
            TL(I,K) = T(I,K) - LCRCP * QCL(I,K)
          ENDDO
        ENDDO

      ENDIF

!-----------------------------------------------------------------------
!! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------


      CALL EX_FLUX_TQ (
     &  P_POINTS,P_FIELD,P1,BL_LEVELS
     &, TL,QW,RDZ,FTL,FQW,RHOKH
     &, RHOKHZ(1,2)
     &, GRAD_T_ADJ,GRAD_Q_ADJ
     &, NTML
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
     & P_FIELD,P1
     &,LAND_INDEX
     &,LAND_PTS,LAND1
     &,P_POINTS,BL_LEVELS,N_TYPES,TILE_FRAC
     &,ALPHA1_GB,ALPHA1,ASHTF
     &,DTRDZ,DTRDZ_RML,RHOKH(1,2),RDZ
     &,ICE_FRACT,LYING_SNOW,RADNET_C,RESFT_TILE,RHOKPM_TILE          
     &,RHOKPM_POT_TILE
     &,TIMESTEP,LAND_MASK
     &,EPOT,EPOT_TILE
     &,FQW,FQW_TILE,FTL,FTL_TILE,E_SEA,H_SEA,DQW_NT,QW        
     &,GAMMA,RHOKE,RHOKH(1,1),DTL_NT,TL
     &,SURF_HT_FLUX,NRML
     &,LTIMER
     &)


!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
!      over sea points, so copy this to remaining tiles
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
Cfpp$ Select(CONCUR)
        DO  I=P1,P1+P_POINTS-1
          FTL(I,K) = FTL(I,K)*CP
        ENDDO
      ENDDO

      DO ITILE=1,N_TYPES
        DO I=P1,P1+P_POINTS-1
          IF(LAND_MASK(I)) THEN
            FTL_TILE(I,ITILE) = FTL_TILE(I,ITILE)*CP
          ELSE
            FTL_TILE(I,ITILE) = FTL(I,1)
            FQW_TILE(I,ITILE) = FQW_TILE(I,1)
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!  Diagnose surface temperature and increment sub-surface temperatures
!!  for land and sea-ice.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!   Sea-ice (P241, routine SICE_HTF).
!-----------------------------------------------------------------------

      CALL SICE_HTF(
     & ASHTF,DI,ICE_FRACT,SURF_HT_FLUX(1,1),TIMESTEP,
     & LAND_MASK,P_FIELD,P_POINTS,P1,TI,TSTAR,ASURF,
     & SEA_ICE_HTF,LTIMER
     &)

!-----------------------------------------------------------------------
!!   Diagnose the land surface temperature (previously in SOIL_HTF)
!-----------------------------------------------------------------------


      DO I=LAND1,LAND1+LAND_PTS-1
        J = LAND_INDEX(I)
        TSTAR(J)=0.0
      ENDDO


      DO ITILE=1,N_TYPES
        DO J=P1,P1+P_POINTS-1
          IF (.NOT. LAND_MASK(J)) TSTAR_TILE(J,ITILE)=TSTAR(J)
        ENDDO

        DO I=LAND1,LAND1+LAND_PTS-1
          J = LAND_INDEX(I)
          TSTAR_TILE(J,ITILE) = T_SOIL(I,1) + SURF_HT_FLUX(J,ITILE)
     &                                       / ASHTF(J)
          TSTAR(J)=TSTAR(J)+TSTAR_TILE(J,ITILE)*TILE_FRAC(J,ITILE)
        ENDDO
      ENDDO ! tile loop

!-----------------------------------------------------------------------
!! 7.  Surface evaporation components and updating of surface
!!     temperature (P245, routine SF_EVAP).
!!     The following diagnostics are also calculated, as requested :-
!!     Heat flux due to melting of sea-ice; specific humidity at 1.5
!!     metres; temperature at 1.5 metres.
!-----------------------------------------------------------------------

      CALL SF_EVAP (
     & P_FIELD,P1,N_TYPES,LAND_FIELD,LAND1,GAMMA,
     & P_POINTS,BL_LEVELS,LAND_MASK,LAND_PTS,LAND_INDEX,
     & TILE_FRAC,ALPHA1,ASURF,ASHTF,CANOPY,CATCH_TILE,
     & DTRDZ,DTRDZ_RML,E_SEA,FRACA,
     & ICE_FRACT,NRML,RHOKH_TILE,SMC,TIMESTEP,CER1P5M,CHR1P5M,
     & PSTAR,RESFS,RESFT_TILE,Z0M,Z0H,SQ1P5,ST1P5,SIMLT,SMLT,
     & FTL,FTL_TILE,FQW,FQW_TILE,LYING_SNOW,QW,SURF_HT_FLUX,
     & TL,TSTAR_TILE,TSTAR,TI,ECAN,ES,EI,
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,
     & H_BLEND,HEAT_BLEND_FACTOR,QCL(1,1),QCF(1,1),Z1,
     & Q1P5M,T1P5M,LTIMER
     & )


!7.1 Copy T and Q from workspace to INOUT space.

      DO K=1,BL_LEVELS
Cfpp$  Select(CONCUR)
        DO I=P1,P1+P_POINTS-1
          T(I,K)=TL(I,K)
          Q(I,K)=QW(I,K)
        ENDDO
      ENDDO

C-----------------------------------------------------------------------
C Diagnose the true value of the surface soil heat flux over land points
C-----------------------------------------------------------------------
      DO ITILE=1,N_TYPES
        DO I=LAND1,LAND1+LAND_PTS-1
          J = LAND_INDEX(I)
          SURF_HT_FLUX(J,ITILE) = SURF_HT_FLUX(J,ITILE)
     +                          - CANCAP(J,ITILE) *
     +            (TSTAR_TILE(J,ITILE) - T_SOIL(I,1)) / TIMESTEP
        ENDDO
      ENDDO
      DO ITILE=1,N_TYPES
        DO I=P1,P1+P_POINTS-1
          SURF_HT_FLUX_GB(I) = SURF_HT_FLUX_GB(I) + TILE_FRAC(I,ITILE)*
     &                                           SURF_HT_FLUX(I,ITILE)
        ENDDO
      ENDDO

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
!! 9.  Calculate diagnostics
!  9.1 Surface latent heat flux.
!-----------------------------------------------------------------------

      IF (SLH) THEN
        DO I=P1,P1+P_POINTS-1
          LATENT_HEAT(I) = LC*FQW(I,1) + LF*EI(I)
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
! 9.2 Diagnose the soil evaporation, the transpiration and the water
!     extracted from each soil layer
!-----------------------------------------------------------------------
      DO ITILE=1,N_TYPES
        DO N=1,SM_LEVELS
          DO I=LAND1,LAND1+LAND_PTS-1
            J = LAND_INDEX(I)
            EXT(I,N)=EXT(I,N) + WT_EXT(I,N) * (1-F_SE(J,ITILE))*
     &                          ES(J,ITILE) * TILE_FRAC(J,ITILE)

          ENDDO ! land_points
        ENDDO ! sm_levels

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=LAND1,LAND1+LAND_PTS-1
          J = LAND_INDEX(I)
          ESOIL(J,ITILE)=F_SE(J,ITILE)*ES(J,ITILE)
          ETRAN(J,ITILE)=(1-F_SE(J,ITILE))*ES(J,ITILE)
          EXT(I,1)=EXT(I,1)+ESOIL(J,ITILE)*TILE_FRAC(J,ITILE)
          ES_GB(J)=ES_GB(J)+ES(J,ITILE)*TILE_FRAC(J,ITILE)
        ENDDO
      ENDDO ! Tile loop

!-----------------------------------------------------------------------
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        RHOKH(I,1) = GAMMA(1)*RHOKH(I,1)
      ENDDO
      DO K = 2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHOKH(I,K) = GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
        ENDDO
      ENDDO

  999  CONTINUE  ! Branch for error exit.

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
