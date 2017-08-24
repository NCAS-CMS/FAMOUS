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
!!! Subroutine BL_INTCT -------------------------------------------
!!!
!!! Purpose : Intermediate control level to call requested version of
!!!           BDY_LAYR with the appropriate arguments.
!!!
!!! Level 3 control routine
!!! version for CRAY YMP
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.3   2/2/97  new deck. S Jackson
!!!  4.4  25/6/97  Modified for MOSES II tile model. R Essery
!!!  4.4  18/09/97 Extra arguments RAD_HR and RADHR_DIM1 for BDYLYR6A
!!!                                                     Cyndy Bunton
!!!  4.4  24/11/97 Move grid definitions up from BDY_LAYR.  R.A.Betts
!!!  4.5   24/4/98 New diagnostics ZHT and BL_TYPE_1 to _6 for 
!!!                BDYLYR6A                             R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!  4.5  24/06/98 Output TILE_FRAC as diagnostic.  R.A.Betts
!!!  4.5  07/09/98 Output GPP_FT and RESP_P_FT as diagnostics.
!!!                                                     Richard Betts
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!! System components covered : P24
!!!
!!! System task : P0
!!!
!!!END -----------------------------------------------------------------
!    Arguments :-
      SUBROUTINE BL_INTCT (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,LAND_FIELD,LAND_FIELD_TRIF,NPFT_TRIF,
     & P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN soil/vegetation/land surface data :
     & LAND_MASK,GATHER,LAND_INDEX,
     & ST_LEVELS,SM_LEVELS,CANHT,CANOPY,CATCH,HCAP,
     & HCON,LAI,LAYER_DEPTH,
     & LYING_SNOW,RESIST,ROOTD,SMC,SMVCCL,SMVCST,SMVCWT,
     & VFRAC,Z0V,SIL_OROG_LAND,L_Z0_OROG,
     & HO2R2_OROG,

! IN sea/sea-ice data :
     & DI,ICE_FRACT,U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,
     & CCA,CCB,CCT,

! IN everything not covered so far :
     & RAD_HR,RADHR_DIM1,
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,
     & TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN,

! INOUT data :
     & GS,Q,STHF,STHU,T,T_DEEP_SOIL,TI,TSTAR,U,V,Z0MSEA,

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,E_SEA,ETRAN,FQW,FTL,GPP,H_SEA,
     & NPP,RESP_P,RHOKH,RHOKM,RIB,SEA_ICE_HTF,
     & TAUX,TAUY,VSHR,ZHT,
     & EPOT,FSMC,

! OUT diagnostic requiring STASH flags :
     & FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,
     & Q1P5M,T1P5M,U10M,V10M,
! (IN) STASH flags :-
     & SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,
     & NRML,

! OUT data required for 4D_VAR :
     & RHO_CD_MODV1,RHO_KM,

! OUT data required elsewhere in UM system :
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,
     & ECAN,EI,ES,EXT,SNOWMELT,
     & SURF_HT_FLUX,ZH,T1_SD,Q1_SD,ERROR,

! Additional arguments for 7A boundary layer (MOSES II)
! IN
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,
     & CANHT_FT,CANOPY_TILE,CATCH_TILE,CS,LAI_FT,
     & FRAC,SNOW_FRAC,RAD_NO_SNOW,RAD_SNOW,TSNOW,Z0V_TILE,
     & CO2_3D,CO2_DIM,L_CO2_INTERACTIVE,
! INOUT
     & TSTAR_TILE,
     & G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,
! OUT
     & ECAN_TILE,ESOIL_TILE,FTL_TILE,
     & G_LEAF,GPP_FT,NPP_FT,RESP_P_FT,RESP_S,RESP_W_FT,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & RIB_TILE,SNOW_SURF_HTF,SOIL_SURF_HTF,
     & TILE_INDEX,TILE_PTS,TILE_FRAC,

! LOGICAL LTIMER
     & LTIMER
     &)
      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency in BDY_LAYR.
!    All dimensions set to 1 for single column model.

      INTEGER
     & P_FIELD                   ! IN No. of P-points in whole grid
!                                !    (for dimensioning only).
     &,U_FIELD                   ! IN No. of UV-points in whole grid.
     &,LAND_FIELD                ! IN No. of land points in whole grid.
     &,LAND_FIELD_TRIF           ! IN For dimensioning land fields
!                                !    available only with TRIFFID
!                                !    Set to LAND_FIELD when TRIFFID on,
!                                !    set to 1 when TRIFFID off.
     &,NPFT_TRIF                 ! IN For dimensioning PFT fields
!                                !    available only with TRIFFID
!                                !    Set to NPFT when TRIFFID on,
!                                !    set to 1 when TRIFFID off.
     &,P_ROWS                    ! IN No. of P-rows in whole grid
!                                !    (for dimensioning only).
     &,FIRST_ROW                 ! IN First row of data to be treated,
!                                !    referred to P-grid.
     &,N_ROWS                    ! IN No. of rows of data to be
!                                !    treated, referred to P-grid.
     &,ROW_LENGTH                ! IN No. of points in one row.

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                 ! IN Max. no. of "boundary" levels
!                                !    allowed.Assumed <= 30 for dim-
!                                !    sioning of GAMMA in common deck
!                                !    C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                  ! IN Total no. of vertical levels in
!                                !    the model atmosphere.
     &,RADHR_DIM1                ! IN Dimension for RAD_HR
      REAL
     & AK(P_LEVELS)              ! IN Hybrid 'A' for all levels.
     &,BK(P_LEVELS)              ! IN Hybrid 'B' for all levels.
     &,AKH(P_LEVELS+1)           ! IN Hybrid 'A' for layer interfaces.
     &,BKH(P_LEVELS+1)           ! IN Hybrid 'B' for layer interfaces.
     &,DELTA_AK(P_LEVELS)        ! IN Difference of hybrid 'A' across
!                                !    layers (K-1/2 to K+1/2).
!                                !    NB: Upper minus lower.
     &,DELTA_BK(P_LEVELS)        ! IN Difference of hybrid 'B' across
!                                !    layers (K-1/2 to K+1/2).
!                                !    NB: Upper minus lower.
     &,EXNER(P_FIELD,BL_LEVELS+1)! IN Exner function.  EXNER(,K) is
!                                !    value for LOWER BOUNDARY of
!                                !    level K.

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL
     & LAND_MASK(P_FIELD)        ! IN T if land, F elsewhere.
     &,L_CO2_INTERACTIVE
     &,L_Z0_OROG                 ! IN T to use simple orog.roughness
!                                !    treatment in SF_EXCH
     &,GATHER                    ! IN T if gather to sea-ice points
!                                !    in SF_EXCH. Saves a lot of un-
!                                !    necessary calculations if there
!                                !    are relatively few sea-ice points

      INTEGER
     & LAND_INDEX(P_FIELD)       ! IN LAND_INDEX(I)=J => the Jth
!                                !    point in P_FIELD is the Ith
!                                !    land point.

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
     & ST_LEVELS                 ! IN No. of deep soil temp. levels
     &,SM_LEVELS                 ! IN No. of soil moisture levels
     &,CO2_DIM          ! number of points in CO2 field.

      REAL
     & CANHT_FT(LAND_FIELD,NPFT) ! IN Canopy height (m)
     &,CANOPY_TILE(LAND_FIELD,NTYPE-1)
!                                ! IN Surface/canopy water for snow-free
!                                !    land tiles (kg per sq m)
     &,CATCH_TILE(LAND_FIELD,NTYPE-1)
!                                ! IN Surface/canopy water capacity of
!                                !    snow-free land tiles (kg per sq m)
     &,CS(LAND_FIELD)            ! IN Soil carbon (kg C/m2).
     &,HCON(LAND_FIELD)          ! IN Soil thermal conductivity (W/m/K).
     &,HO2R2_OROG(LAND_FIELD)    ! IN Dummy used only in version 3A.
     &,LAI_FT(LAND_FIELD,NPFT)   ! IN Leaf area index
     &,LYING_SNOW(P_FIELD)       ! IN Lying snow (kg per sq m).
!                                !    Must be global for coupled model,
!                                !    ie dimension P_FIELD not LAND_FIEL
     &,SIL_OROG_LAND(LAND_FIELD) ! IN Silhouette area of unresolved
!                                !    orography per unit horizontal area
!                                !    on land points only.
     &,SMVCCL(LAND_FIELD)        ! IN Critical volumetric SMC (cubic m
!                                !    per cubic m of soil).
     &,SMVCST(LAND_FIELD)        ! IN Volumetric saturation point (cubic
!                                !    per cubic m of soil).
     &,SMVCWT(LAND_FIELD)        ! IN Volumetric wilting point (cubic m
!                                !    per cubic m of soil).
     &,STHF(LAND_FIELD,SM_LEVELS)! IN Frozen soil moisture content of
!                                !    each layer as a fraction of
!                                !    saturation.
     &,STHU(LAND_FIELD,SM_LEVELS)! IN Unfrozen soil moisture content of
!                                !    each layer as a fraction of
!                                !    saturation.
     &,FRAC(LAND_FIELD,NTYPE)    ! IN Tile fracs excluding snow cover
     &,SNOW_FRAC(LAND_FIELD)     ! IN Snow fraction.
     &,TSNOW(LAND_FIELD)         ! IN Snow surface layer temp. (K).
     &,Z0V(P_FIELD)              ! IN GBM snow-free roughness length (m)
!                                !    NB:UM uses same storage for Z0MSEA
!                                !    so for sea points this is INOUT.
     &,Z0V_TILE(LAND_FIELD,NTYPE)! IN Tile roughness lengths (m).

! (d) Sea/sea-ice data.

      REAL
     & DI(P_FIELD)               ! IN "Equivalent thickness" of sea-ice
!                                !    (m).
     &,ICE_FRACT(P_FIELD)        ! IN Fraction of gridbox covered by
!                                !    sea-ice (decimal fraction).
     &,U_0(U_FIELD)              ! IN W'ly component of surface current
!                                !    (metres per second).
     &,V_0(U_FIELD)              ! IN S'ly component of surface current
!                                !    (metres per second).

! (e) Cloud data.

      REAL
     & CF(P_FIELD,BL_LEVELS)     ! IN Cloud fraction (decimal).
     &,QCF(P_FIELD,BL_LEVELS)    ! IN Cloud ice (kg per kg air)
     &,QCL(P_FIELD,BL_LEVELS)    ! IN Cloud liquid water (kg/kg air).
     &,CCA(P_FIELD)              ! IN Convective Cloud Amount (decimal).

      INTEGER
     & CCB(P_FIELD)              ! IN Convective Cloud Base
     &,CCT(P_FIELD)              ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL
     & CO2_MMR                   ! IN CO2 Mass Mixing Ratio
     &,CO2_3D(CO2_DIM)  ! 3D CO2 field if required.
     &,PHOTOSYNTH_ACT_RAD(P_FIELD)! IN Net downward shortwave radiation
!                                !     in band 1 (w/m2).
     &,PSTAR(P_FIELD)            ! IN Surface pressure (Pascals).
     +,RADNET(P_FIELD)           ! IN Surface net radiation (W/sq m,
C                                !    positive downwards).
     &,RAD_NO_SNOW(P_FIELD)      ! IN Surface net radiation, snow-free
!                                !    fraction of gridbox.
     &,RAD_SNOW(P_FIELD)         ! IN Surface net radiation, snow-
!                                !    covered fraction of gridbox.
     +,RAD_HR(RADHR_DIM1,BL_LEVELS)
!                                ! IN Radiative heating rates
!                                !    - not used in A03_7A.
     &,TIMESTEP                  ! IN Timestep (seconds).

      LOGICAL
     & LTIMER                    ! IN Logical switch for TIMER diags
     &,L_RMBL                    ! IN T to use rapidly mixing
!                                !    boundary scheme in IMPL_CAL
     &,L_BL_LSPICE               ! IN Use if 3A large scale precip
     &,L_MOM                     ! IN Switch for convective momentum
!                                !    transport.
     &,L_PHENOL                  ! IN Indicates whether phenology in use
     &,L_TRIFFID                 ! IN Indicates whether TRIFFID in use.
     &,L_NEG_TSTAR              ! IN Switch for -ve TSTAR error check

!  STASH flags :-

      LOGICAL
     & SFME    ! IN Flag for FME (q.v.).
     &,SMLT    ! IN Flag for SICE_MLT_HTF (q.v.)
     &,SIMLT   ! IN Flag
     &,SLH     ! IN Flag for LATENT_HEAT (q.v.)
     &,SQ1P5   ! IN Flag for Q1P5M (q.v.)
     &,ST1P5   ! IN Flag for T1P5M (q.v.)
     &,SU10    ! IN Flag for U10M (q.v.)
     &,SV10    ! IN Flag for V10M (q.v.)

!  In/outs :-

      REAL
     & GS(LAND_FIELD)            ! INOUT "Stomatal" conductance to
!                                !       evaporation (m/s).
     &,Q(P_FIELD,BL_LEVELS)      ! INOUT Input:specific humidity
!                                !       ( kg water per kg air).
!                                !      Output:total water content
!                                !      (Q)(kg water per kg air).
     &,T(P_FIELD,BL_LEVELS)      ! INOUT Input:atmospheric temp(K)
!                                !      Output:liquid/frozen water
!                                !       temperature (TL) (K)
     &,T_DEEP_SOIL(LAND_FIELD,ST_LEVELS)
!                                ! INOUT Deep soil temperatures (K).
     &,TI(P_FIELD)               ! INOUT Sea-ice surface layer
!                                !       temperature (K)
     &,TSTAR(P_FIELD)            ! INOUT Surface temperature (K).
     &,TSTAR_TILE(LAND_FIELD,NTYPE)
!                                ! INOUT Surface tile temperature
     &,U(U_FIELD,BL_LEVELS)      ! INOUT W'ly wind component (m/s).
     &,V(U_FIELD,BL_LEVELS)      ! INOUT S'ly wind component (m/s).
     &,Z0MSEA(P_FIELD)           ! INOUT Sea-surface roughness
!                                !       length for momentum (m).
!                                !       NB: same storage is used
!                                !       for Z0V, so the intent is
!                                !       IN for land points.

!  Accumulation prognostics for PHENOLOGY and TRIFFID.
!  NPP_FT_ACC, RESP_W_FT_ACC and RESP_S_ACC are only allocated D1 space
!  when TRIFFID is in use, so their dimensions here are set accordingly.

      REAL
     & G_LEAF_ACC(LAND_FIELD,NPFT)          ! INOUT Accumulated G_LEAF
     &,NPP_FT_ACC(LAND_FIELD_TRIF,NPFT_TRIF)! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_FIELD_TRIF,NPFT_TRIF) ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_FIELD_TRIF)          ! INOUT Accumulated RESP_S


!  Outputs :-

!-1 Diagnostic (or effectively so - includes coupled model requisites):-

      INTEGER
     & TILE_INDEX(LAND_FIELD,NTYPE)
!                               ! OUT Index of tile points.
     &,TILE_PTS(NTYPE)          ! OUT Number of tile points.

!  (a) Calculated anyway (use STASH space from higher level) :-

      REAL
     & CD(P_FIELD)              ! OUT Turbulent surface exchange (bulk
!                               !     transfer) coefficient for
!                               !     momentum.
     &,CH(P_FIELD)              ! OUT Turbulent surface exchange (bulk
!                               !     transfer) coefficient for heat
!                               !     and/or moisture.
     &,ECAN(P_FIELD)            ! OUT Gridbox mean evaporation from
!                               !     canopy / surface store (kg/m2/s).
!                               !     Zero over sea.
     &,E_SEA(P_FIELD)           ! OUT Evaporation from sea times leads
!                               !     fraction. Zero over land.
!                               !     (kg per square metre per sec).
     &,EPOT(P_FIELD)            ! Dummy.
     &,ESOIL_TILE(LAND_FIELD,NTYPE-1)
                                ! OUT ES for snow-free land tiles
     &,FQW(P_FIELD,BL_LEVELS)   ! OUT Moisture flux between layers
!                               !     (kg per square metre per sec).
!                               !     FQW(,1) is total water flux
!                               !     from surface, 'E'.
     &,FSMC(LAND_FIELD)         ! Dummy.
     &,FTL(P_FIELD,BL_LEVELS)   ! OUT FTL(,K) contains net turbulent
!                               !     sensible heat flux into layer K
!                               !     from below; so FTL(,1) is the
!                               !     surface sensible heat, H.  (W/m2)
     &,FTL_TILE(LAND_FIELD,NTYPE)
!                               ! OUT Surface FTL for land tiles
     &,G_LEAF(LAND_FIELD,NPFT)  ! OUT Leaf turnover rate (/360days).
     &,GPP(LAND_FIELD)          ! OUT Gross primary productivity
!                               !     (kg C/m2/s).
     &,GPP_FT(LAND_FIELD,NPFT)  ! OUT Gross primary productivity
!                               !     on PFTs (kg C/m2/s).
     &,H_SEA(P_FIELD)           ! OUT Surface sensible heat flux over
!                               !     sea times leads fraction. (W/m2)
     &,NPP(LAND_FIELD)          ! OUT Net primary productivity
!                               !      (kg C/m2/s).
     &,NPP_FT(LAND_FIELD,NPFT)  ! OUT Net primary productivity
!                               !     (kg C/m2/s).
     &,RESP_P(LAND_FIELD)       ! OUT Plant respiration (kg C/m2/s).
     &,RESP_P_FT(LAND_FIELD,NPFT) ! OUT Plant respiration on PFTs
!                                 !     (kg C/m2/s).

     &,RESP_S(LAND_FIELD)       ! OUT Soil respiration (kg C/m2/s).
     &,RESP_W_FT(LAND_FIELD,NPFT)! OUT Wood maintenance respiration
!                                !     (kg C/m2/s).
     &,RHOKH(P_FIELD,BL_LEVELS) ! OUT Exchange coeffs for moisture.
     &,RHOKM(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                               !     momentum (on UV-grid, with 1st
!                               !     and last rows undefined (or, at
!                               !     present, set to "missing data")).
     &,RIB(P_FIELD)             ! OUT Bulk Richardson number for lowest
!                               !     layer.
     &,RIB_TILE(LAND_FIELD,NTYPE)! OUT RIB for land tiles.
     &,SEA_ICE_HTF(P_FIELD)     ! OUT Heat flux through sea-ice (W per
!                               !     sq m, positive downwards).
     &,SMC(LAND_FIELD)          ! OUT Available moisture in the
!                               !     soil profile (mm).
     &,SURF_HT_FLUX(P_FIELD)    ! OUT Net downward heat flux at surface
!                               !     over land or sea-ice fraction of
!                               !     gridbox (W/m2)
     &,TAUX(U_FIELD,BL_LEVELS)  ! OUT W'ly component of surface wind
!                               !     stress (N/sq m).(On UV-grid with
!                               !     first and last rows undefined or
!                               !     at present, set to 'missing data'
     &,TAUY(U_FIELD,BL_LEVELS)  ! OUT S'ly component of surface wind
!                               !     stress (N/sq m).  On UV-grid;
!                               !     comments as per TAUX.
     &,TILE_FRAC(LAND_FIELD,NTYPE)
!                               ! OUT Tile fractions adjusted for snow.
!                               !     1 to NTYPE-1: snow-free fraction.
!                               !     NTYPE:land-ice plus snow fraction.
     &,VSHR(P_FIELD)            ! OUT Magnitude of surface-to-lowest
!                               !     atm level wind shear (m per s).
     &,RHO_CD_MODV1(P_FIELD)    ! OUT Surface air density * drag coef.
!                               ! mod(v1 - v0) before interpolation.
     &,RHO_KM(P_FIELD,2:BL_LEVELS)! OUT Air density * turbulent mixing
!                               ! coef. for momentum before
     &,RHO_ARESIST(P_FIELD)     ! OUT, RHOSTAR*CD_STD*VSHR for SCYCLE
     &,ARESIST(P_FIELD)         ! OUT, 1/(CD_STD*VSHR)    for SCYCLE
     &,RESIST_B(P_FIELD)        ! OUT,(1/CH-1/CD_STD)/VSHR for SCYCLE
     &,RHO_ARESIST_TILE(LAND_FIELD,NTYPE)
!                               ! OUT RHOSTAR*CD_STD*VSHR on land tiles
     &,ARESIST_TILE(LAND_FIELD,NTYPE)
!                               ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_FIELD,NTYPE)
!                               ! OUT (1/CH-1/CD_STD)/VSHR on land tiles
!
      INTEGER
     & NRML(P_FIELD)            ! OUT Number of model layers in the
!                               !     Rapidly Mixing Layer; diagnosed
!                               !     in SF_EXCH and KMKH and used in
!                               !     IMPL_CAL, SF_EVAP and TR_MIX.

! (b) Not passed between lower-level routines (not in workspace at this
!     level) :-

      REAL
     & FME(P_FIELD)             ! OUT Wind mixing "power" (W per sq m).
     &,SICE_MLT_HTF(P_FIELD)    ! OUT Heat flux due to melting of sea-
!                               !     ice (Watts per sq metre).
     &,SNOMLT_SURF_HTF(P_FIELD)
     &,LATENT_HEAT(P_FIELD)     ! OUT Surface latent heat flux, +ve
!                               !     upwards (Watts per sq m).
     &,Q1P5M(P_FIELD)           ! OUT Q at 1.5 m (kg water per kg air).
     &,T1P5M(P_FIELD)           ! OUT T at 1.5 m (K).
     &,U10M(U_FIELD)            ! OUT U at 10 m (m per s).
     &,V10M(U_FIELD)            ! OUT V at 10 m (m per s).
     &,ZHT(P_FIELD)              ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_1(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_2(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_3(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_4(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_5(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_6(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)


!-2 Genuinely output, needed by other atmospheric routines :-

      REAL
     & ECAN_TILE(LAND_FIELD,NTYPE-1)
                      ! OUT ECAN for snow-free land tiles
     &,EI(P_FIELD)    ! OUT Sublimation from lying snow or sea-ice
!                     !     (kg/m2/s).
     &,ES(P_FIELD)    ! OUT Surface evapotranspiration from soil
!                     !     moisture store (kg/m2/s).
     &,EXT(LAND_FIELD,SM_LEVELS)
!                     ! OUT Extraction of water from each soil layer
!                     !     (kg/m2/s).
     &,SNOWMELT(P_FIELD)
!                     ! OUT Snowmelt (kg/m/s).
     &,SNOW_SURF_HTF(LAND_FIELD)
!                     ! OUT Net downward heat flux at
!                     !     snow surface (W/m2).
     &,SOIL_SURF_HTF(LAND_FIELD)
!                     ! OUT Net downward heat flux at
!                     !     snow-free land surface (W/m2).
     &,ZH(P_FIELD)    ! OUT Height above surface of top of boundary
!                     !     layer (metres).
     &,T1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
!                     !     of layer 1 temperature; for use in
!                     !     initiating convection.
     &,Q1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
!                     !     of layer 1 humidity; for use in initiating
!                     !     convection.

      INTEGER
     & ERROR          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;

! Local variables
      REAL
     & GS_TILE(LAND_FIELD,NTYPE)! LOCAL Surface conductance for
!                               !       land tiles
     &,WT_EXT(LAND_FIELD,SM_LEVELS)
!                               ! LOCAL Fraction of evapotranspiration
!                               !       which is extracted from each
!                               !       soil layer.
     &,Z1(P_FIELD)              ! LOCAL Height of lowest level (m).
      INTEGER
     & TILE_INDEX_S(LAND_FIELD,NTYPE)
!                               ! LOCAL Index for TILE_FRAC.
     &,TILE_PTS_S(NTYPE)        ! LOCAL Number of points for TILE_FRAC.

      INTEGER
     & I                        ! LOCAL P-point index
     &,L                        ! LOCAL Land point index
     &,N                        ! LOCAL Tile index
     &,N_P_ROWS                 ! LOCAL No of P-rows being processed.
     &,N_U_ROWS                 ! LOCAL No of UV-rows being processed.
     &,P_POINTS                 ! LOCAL No of P-points being processed.
     &,P1                       ! LOCAL First P-point to be processed.
     &,LAST_POINT               ! LOCAL Last P-point to be processed.
     &,LAND1                    ! LOCAL First land-point to be processed
!                               !       1 <= LAND1 <= LAND_FIELD
     &,LAND_PTS                 ! LOCAL No of land points processed.
     &,U_POINTS                 ! LOCAL No of UV-points being processed.
     &,U1                       ! LOCAL First UV-point to be processed.

      REAL
     & SECS_PER_360DAYS         ! LOCAL Number of seconds in 360 days

      PARAMETER(SECS_PER_360DAYS=31104000.0)



! Dummy variables not used by MOSES II
      REAL
     & CANHT(LAND_FIELD)
     &,CANOPY(LAND_FIELD)
     &,CATCH(LAND_FIELD)
     &,ETRAN(P_FIELD)
     &,HCAP(LAND_FIELD)
     &,LAI(LAND_FIELD)
     &,LAYER_DEPTH(SM_LEVELS)
     &,RESIST(LAND_FIELD)
     &,ROOTD(LAND_FIELD)
     &,VFRAC(LAND_FIELD)
      LOGICAL
     & L_MIXLEN

! External subroutines called
      EXTERNAL
     & TILEPTS     ! Calculates number of points occupied by each
!                  ! tile and their indices on the land field
     &,VSHR_Z1     ! Calculates level 1 windspeed and height.
     &,PHYSIOL     ! Models plant physiology
     &,BDY_LAYR    ! Models surface fluxes and boundary layer processes

!-----------------------------------------------------------------------
!! 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
!!    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
!!    most UV-rows are used only for interpolation and are not updated.
!-----------------------------------------------------------------------

      IF ( BL_LEVELS.LT.1 .OR. SM_LEVELS.LT.1 .OR. P_ROWS.LT.3 ) THEN
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

      LAST_POINT = P1 + P_POINTS - 1

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
! Call TILEPTS to calculate TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------
      CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &             FRAC,TILE_PTS,TILE_INDEX)



!-----------------------------------------------------------------------
! Call MOSES II physiology routine to calculate surface conductances
! and carbon fluxes.
! VSHR_Z1 provides level 1 windspeed and height.
!-----------------------------------------------------------------------

      CALL VSHR_Z1 (
     & P_FIELD,U_FIELD,LTIMER,
     & N_ROWS,FIRST_ROW,ROW_LENGTH,
     & AKH,BKH,EXNER,PSTAR,Q,QCF,QCL,T,U,V,U_0,V_0,
     & VSHR,Z1
     & )


      CALL PHYSIOL (
     & LAND_FIELD,LAND_PTS,LAND1,
     & LAND_INDEX,
     & P_FIELD,SM_LEVELS,TILE_PTS,TILE_INDEX,
     & CO2_MMR,CO2_3D,CO2_DIM,L_CO2_INTERACTIVE,
     & CS,FRAC,CANHT_FT,PHOTOSYNTH_ACT_RAD,
     & LAI_FT,PSTAR,Q,STHU,TIMESTEP,T_DEEP_SOIL,TSTAR_TILE,
     & SMVCCL,SMVCST,SMVCWT,VSHR,Z0V_TILE,Z1,
     & G_LEAF,GS,GS_TILE,GPP,GPP_FT,NPP,NPP_FT,
     & RESP_P,RESP_P_FT,RESP_S,RESP_W_FT,SMC,WT_EXT
     & )

!----------------------------------------------------------------------
! Increment accumulation of leaf turnover rate.
! This is required for leaf phenology and/or TRIFFID, either of
! which can be enabled independently of the other.
!----------------------------------------------------------------------
      IF (L_PHENOL.OR.L_TRIFFID) THEN
        DO N=1,NPFT
          DO L=LAND1,LAND1+LAND_PTS-1
            G_LEAF_ACC(L,N) = G_LEAF_ACC(L,N) +
     &      G_LEAF(L,N)*(TIMESTEP/SECS_PER_360DAYS)
          ENDDO
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Increment accumulation prognostics for TRIFFID
!----------------------------------------------------------------------
      IF (L_TRIFFID) THEN
        DO N=1,NPFT
          DO L=LAND1,LAND1+LAND_PTS-1
            NPP_FT_ACC(L,N) = NPP_FT_ACC(L,N) + NPP_FT(L,N)*TIMESTEP
            RESP_W_FT_ACC(L,N) = RESP_W_FT_ACC(L,N)
     &                                      + RESP_W_FT(L,N)*TIMESTEP
          ENDDO
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          RESP_S_ACC(L) = RESP_S_ACC(L) + RESP_S(L)*TIMESTEP
        ENDDO
      ENDIF


!-----------------------------------------------------------------------
! Calculate modified snow-free tile fractions for all but the ice tile
!-----------------------------------------------------------------------
      DO N=1,NTYPE-1
        DO L=1,LAND_FIELD
          TILE_FRAC(L,N) = (1. - SNOW_FRAC(L))*FRAC(L,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate the areal fraction of an "ice plus snow" tile by adding the
! snow-covered fractions of all other tiles onto the areal fraction of
! the land-ice tile
!-----------------------------------------------------------------------
      N = NTYPE
      DO L=1,LAND_FIELD
        TILE_FRAC(L,N) = FRAC(L,N) + SNOW_FRAC(L)*(1-FRAC(L,N))
      ENDDO

!-----------------------------------------------------------------------
! Call TILEPTS to calculate TILE_PTS_S and TILE_INDEX_S
!-----------------------------------------------------------------------
      CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &             TILE_FRAC,TILE_PTS_S,TILE_INDEX_S)

!-----------------------------------------------------------------------
! Call boundary layer routine carrying-out tile calculations on
! snow-modified tiles
!-----------------------------------------------------------------------

      CALL BDY_LAYR (

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
     & NTYPE,TILE_INDEX_S,TILE_PTS_S,SM_LEVELS,
     & CANOPY_TILE,CATCH_TILE,GS_TILE,HCON,HO2R2_OROG,LYING_SNOW,
     & SIL_OROG_LAND,SMC,SMVCST,STHF,STHU,
     & TILE_FRAC,WT_EXT,Z0V,Z0V_TILE,

! IN sea/sea-ice data :
     & DI,ICE_FRACT,U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & PSTAR,RAD_NO_SNOW,RAD_SNOW,TIMESTEP,VSHR,
     & L_RMBL,L_BL_LSPICE,L_MOM,L_NEG_TSTAR,

! IN STASH flags :-
     & SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,

! INOUT data :
     & Q,T,T_DEEP_SOIL,TSNOW,TI,TSTAR,TSTAR_TILE,
     & U,V,Z0MSEA,

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,ECAN,E_SEA,ESOIL_TILE,FQW,
     & FTL,FTL_TILE,H_SEA,RHOKH,RHOKM,
     & RIB,RIB_TILE,SEA_ICE_HTF,SURF_HT_FLUX,TAUX,TAUY,

! OUT diagnostic requiring STASH flags :
     & FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,
     & Q1P5M,T1P5M,U10M,V10M,

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & NRML,

! OUT data required for 4D_VAR :
     & RHO_CD_MODV1,RHO_KM,

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,ZH,
     & SOIL_SURF_HTF,SNOW_SURF_HTF,
     & T1_SD,Q1_SD,ERROR,

! LOGICAL LTIMER
     & LTIMER
     & )

  999  CONTINUE  ! Branch for error exit.

      RETURN
      END

