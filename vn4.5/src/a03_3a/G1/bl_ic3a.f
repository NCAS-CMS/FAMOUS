CLL Subroutine BL_INTCT -------------------------------------------
CLL
CLL Purpose : Intermediate control level to call requested version of
CLL           BDY_LAYR with the appropriate arguments.
CLL
CLL Level 3 control routine
CLL version for CRAY YMP
CLL
CLL C.Bunton    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  4.0  22/05/95  Original code to introduce intermediate control
CLL                 level-all arguments for every version
CLL                 are passed from BL_CTL and then this routine
CLL                 calls appropriate version of BDY_LAYR with a
CLL                 subset of these
CLL                 (Outlined in proposal by S. Foreman 22/08/94)
CLL                 C.Bunton
CLL  4.0  24/4/95   Diagnostics for 4D-Var project added
CLL                           Simon Jackson
CLL  4.1  07/05/96  3 surface resistance coefficients to argument list
CLL                 for use in Sulphur Cycle      M. Woodage.
CLL  4.1  20/5/96   Correct stash flags for sea-ice melt and snow melt
CLL                 heat fluxes      Simon Jackson
CLL  4.3  04/02/97  Logical switches L_MOM and L_MIXLEN passed down
CLL                 to BDY_LYR
CLL                                                     R.N.B.Smith
CLL  4.4  08/09/97 L_BL_LSPICE specifies mixed phase precipitation
CLL                scheme.                        D.Wilson
CLL
CLL
!!!  4.4  18/09/97 Extra arguments RAD_HR and RADHR_DIM1 for BDYLYR6A
!!!                                                     Cyndy Bunton
!!!  4.5   24/4/98 New diagnostics ZHT and BL_TYPE_1 to _6 for 
!!!                BDYLYR6A                             R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!  4.5  27/04/98 Add dummy arguments LAND_FIELD_TRIF, NPFT_TRIF,
!!!                GPP_FT, RESP_P_FT, TILE_FRAC, L_PHENOL and L_TRIFFID
!!!                and change TILE_INDEX and TILE_PTS from IN to OUT,
!!!                consistent with changes to BL_CTL required for
!!!                TRIFFID.
!!!                                                     Richard Betts
CLL Programming standard : unified model documentation paper No 3
CLL
CLL System components covered : P24
CLL
CLL System task : P0
CLL
CLLEND -----------------------------------------------------------------
C    Arguments :-
      SUBROUTINE BL_INTCT(

C IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,LAND_FIELD,LAND_FIELD_TRIF,NPFT_TRIF,
     + P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH,

C IN values defining vertical grid of model atmosphere :
     + BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     + EXNER,

C IN soil/vegetation/land surface data :
     + LAND_MASK,GATHER,LAND_INDEX,
     + ST_LEVELS,SM_LEVELS,CANHT,CANOPY,CATCH,HCAP,
     + HCON,LAI,LAYER_DEPTH,
     + LYING_SNOW,RESIST,ROOTD,SMC,SMVCCL,SMVCST,SMVCWT,
     + VFRAC,Z0V,SIL_OROG_LAND,L_Z0_OROG,
     + HO2R2_OROG,

C IN sea/sea-ice data :
     + DI,ICE_FRACT,U_0,V_0,

C IN cloud data :
     + CF,QCF,QCL,
     + CCA,CCB,CCT,

C IN everything not covered so far :
     & RAD_HR,RADHR_DIM1,
     + CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,TIMESTEP,
     + L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN,

C INOUT data :
     + GS,Q,STHF,STHU,T,T_DEEP_SOIL,TI,TSTAR,U,V,Z0MSEA,

C OUT Diagnostic not requiring STASH flags :
     + CD,CH,E_SEA,ETRAN,FQW,FTL,GPP,H_SEA,
     + NPP,RESP_P,RHOKH,RHOKM,RIB,SEA_ICE_HTF,
     + TAUX,TAUY,VSHR,ZHT,
     + EPOT,FSMC,

C OUT diagnostic requiring STASH flags :
     + FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,
     + Q1P5M,T1P5M,U10M,V10M,
C (IN) STASH flags :-
     + SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,

C OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,
     & NRML,

C OUT data required for 4D_VAR :
     & RHO_CD_MODV1,RHO_KM,

C OUT data required elsewhere in UM system :
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,
     + ECAN,EI,ES,EXT,SNOWMELT,
     + SURF_HT_FLUX,ZH,T1_SD,Q1_SD,ERROR,

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

C LOGICAL LTIMER
     + LTIMER
     +)
      IMPLICIT NONE
C
C  Inputs :-
C
C (a) Defining horizontal grid and subset thereof to be processed.
C
      INTEGER
     + P_FIELD                    ! IN No. of P-points in whole grid
C                                 !    (for dimensioning only).
     +,RADHR_DIM1                 ! IN Dimension for RAD_HR 
     +,U_FIELD                    ! IN No. of UV-points in whole grid.
C                                 !    (Checked for consistency with
C                                 !    P_FIELD and P_ROWS; there must
C                                 !    be 1 less UV than P row.)
     +,LAND_FIELD                 ! IN No.of land points in whole grid.
C                                 !    (Checked for consistency with
C                                 !    P_FIELD )
     +,P_ROWS                     ! IN No. of P-rows in whole grid
C                                 !    (for dimensioning only).
     +,FIRST_ROW                  ! IN First row of data to be treated,
C                                 !    referred to P-grid (must be > 1
C                                 !    since "polar" rows are never
C                                 !    treated).
     +,N_ROWS                     ! IN No. of rows of data to be
C                                 !    treated, referred to P-grid.
C                                 !    FIRST_ROW+N_ROWS-1 must be less
C                                 !    than P_ROWS, since "polar" rows
C                                 !    are never treated.
     +,ROW_LENGTH                 ! IN No. of points in one row.
C                                 !    (Checked for consistency with
C                                 !    P_FIELD and N_ROWS.)
C
C (b) Defining vertical grid of model atmosphere.
C
      INTEGER
     + BL_LEVELS                  ! IN Max. no. of "boundary" levels
C                                 !    allowed.Assumed <= 30 for dim-
C                                 !    sioning of GAMMA in common deck
C                                 !    C_GAMMA used in SF_EXCH and KMKH
     +,P_LEVELS                   ! IN Total no. of vertical levels in
C                                 !    the model atmosphere.
      REAL
     + AK(P_LEVELS)                ! IN Hybrid 'A' for all levels.
     +,BK(P_LEVELS)                ! IN Hybrid 'B' for all levels.
     +,AKH(P_LEVELS+1)             ! IN Hybrid 'A' for layer interfaces.
     +,BKH(P_LEVELS+1)             ! IN Hybrid 'B' for layer interfaces.
     +,DELTA_AK(P_LEVELS)          ! IN Difference of hybrid 'A' across
C                                  !    layers (K-1/2 to K+1/2).
C                                  !    NB: Upper minus lower.
     +,DELTA_BK(P_LEVELS)          ! IN Difference of hybrid 'B' across
C                                  !    layers (K-1/2 to K+1/2).
C                                  !    NB: Upper minus lower.
     +,EXNER(P_FIELD,BL_LEVELS+1)  ! IN Exner function.  EXNER(,K) is
C                                  !    value for LOWER BOUNDARY of
C                                  !    level K.
C
C (c) Soil/vegetation/land surface parameters (mostly constant).
C
      LOGICAL
     + LAND_MASK(P_FIELD)        ! IN T if land, F elsewhere.
     &,L_CO2_INTERACTIVE   ! dummy CO2 variable - not used with A03_3A.
     +,GATHER                    ! IN T if gather to sea-ice points
C                                !    in SF_EXCH. Saves a lot of un-
C                                !    necessary calculations if there
C                                !    are relatively few sea-ice points
     +,L_RMBL                    ! IN T to use rapidly mixing
C     !    boundary scheme in IMPL_CAL          AJS   
     +,L_BL_LSPICE           ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
     +,L_Z0_OROG                 ! IN T to use simple orog.roughness
     &,L_MOM                     ! IN Switch for convective momentum
C                                !    transport.
     &,L_MIXLEN                  ! IN Switch for reducing the turbulent
C                                !    mixing length above the top of the
C                                !    boundary layer.
C                                !    treatment in SF_EXCH
      INTEGER
     + LAND_INDEX(P_FIELD)       ! IN LAND_INDEX(I)=J => the Jth
C                                !    point in P_FIELD is the Ith
C                                !    land point.
      INTEGER
     + ST_LEVELS                 ! IN No. of deep soil temp levls
     +,SM_LEVELS                 ! IN No. of soil moisture levls
     &,CO2_DIM       ! dummy CO2 variable - field not used with A03_3A.
C                                !
      REAL
     + CANHT(LAND_FIELD)         ! IN Canopy height (m)
     +,CANOPY(LAND_FIELD)        ! IN Surface/canopy water (kg per sq m)
     +,CATCH(LAND_FIELD)         ! IN Surface/canopy water capacity
C                                !    (kg per sq m).
     +,HCAP(LAND_FIELD)          ! IN Soil heat capacity (J/K/m**3)
     +,HCON(LAND_FIELD)          ! IN Soil thermal conductivity (W/m/K).
     +,HO2R2_OROG(LAND_FIELD)    ! IN Dummy used only in version 3A.

     +,LAI(LAND_FIELD)           ! IN Leaf area index
     +,LAYER_DEPTH(SM_LEVELS)    !    Dummy Variable not used
C                                !    in this version
     +,LYING_SNOW(P_FIELD)       ! IN Lying snow (kg per sq m).
C                                !    Must be global for coupled model,
C                                !    ie dimension P_FIELD not LAND_FIEL
     +,RESIST(LAND_FIELD)        ! IN "Stomatal" resistance to
C                                !    evaporation (seconds per metre).
     +,ROOTD(LAND_FIELD)         ! IN Depth of active soil layer ("root
C                                !    depth") (metres).
     +,SMC(LAND_FIELD)           ! IN Soil moisture content (kg / sq m).
     +,SMVCCL(LAND_FIELD)        ! IN Critical volumetric SMC (cubic m
C                                !    per cubic m of soil).
     +,SMVCST(LAND_FIELD)        ! IN Volumetric saturation point (cubic
C                                !    per cubic m of soil).
     +,SMVCWT(LAND_FIELD)        ! IN Volumetric wilting point (cubic m
C                                !    per cubic m of soil).
     +,VFRAC(LAND_FIELD)         ! IN Vegetation fraction.
     +,Z0V(P_FIELD)              ! IN Vegetative roughness length (m).
C                                !    NB:UM uses same storage for Z0MSEA
C                                !    so for sea points this is INOUT.
     +,SIL_OROG_LAND(LAND_FIELD) ! IN Dummy used only in version 3A.
C
C (d) Sea/sea-ice data.
C
      REAL
     + DI(P_FIELD)               ! IN "Equivalent thickness" of sea-ice
C                                !    (m).
     +,ICE_FRACT(P_FIELD)        ! IN Fraction of gridbox covered by
C                                !    sea-ice (decimal fraction).
     +,U_0(U_FIELD)              ! IN W'ly component of surface current
C                                !    (metres per second).
     +,V_0(U_FIELD)              ! IN S'ly component of surface current
C                                !    (metres per second).
C
C (e) Cloud data.
C
      REAL
     + CF(P_FIELD,BL_LEVELS)     ! IN Cloud fraction (decimal).
     +,QCF(P_FIELD,BL_LEVELS)    ! IN Cloud ice (kg per kg air)
     +,QCL(P_FIELD,BL_LEVELS)    ! IN Cloud liquid water (kg
C                                !       per kg air).
     +,CCA(P_FIELD)              ! IN Convective Cloud Amount (decimal).
      INTEGER
     + CCB(P_FIELD)              ! IN Convective Cloud Base
     +,CCT(P_FIELD)              ! IN Convective Cloud Top
C
C (f) Atmospheric + any other data not covered so far, incl control.
C
      REAL
     + CO2_MMR                   ! IN CO2 Mass Mixing Ratio
     &,CO2_3D        ! dummy CO2 variable - field not used with A03_3A.
     +,PHOTOSYNTH_ACT_RAD(P_FIELD) ! IN Net downward shortwave radiation
C                                !    in band 1 (w/m2).
     +,PSTAR(P_FIELD)            ! IN Surface pressure (Pascals).
     +,RADNET(P_FIELD)           ! IN Surface net radiation (W/sq m,
C                                !    positive downwards).
     &,RAD_HR(RADHR_DIM1,BL_LEVELS)
!                                 ! IN Radiative heating rates 
!                                 !    - not used in A03_3A.
     +,TIMESTEP                  ! IN Timestep (seconds).
C
      LOGICAL LTIMER             ! Logical switch for TIMER diags
C
C  STASH flags :-
C
      LOGICAL
     + SFME    ! IN Flag for FME (q.v.).
     +,SIMLT   ! IN Flag for SICE_MLT_HTF (q.v.)
     +,SMLT    ! IN Flag for SNOWMLT_SURF_HTF (q.v.)
     +,SLH     ! IN Flag for LATENT_HEAT (q.v.)
     +,SQ1P5   ! IN Flag for Q1P5M (q.v.)
     +,ST1P5   ! IN Flag for T1P5M (q.v.)
     +,SU10    ! IN Flag for U10M (q.v.)
     +,SV10    ! IN Flag for V10M (q.v.)
C
C  In/outs :-
C
      REAL
     + GS(LAND_FIELD)                  ! INOUT "Stomatal" conductance to
C                                      !       evaporation (m/s).
     +,Q(P_FIELD,BL_LEVELS)            ! INOUT Input:specific humidity
C                                      !       ( kg water per kg air).
C                                      !      Output:total water content
C                                      !      (Q)(kg water per kg air).
     +,STHF(LAND_FIELD,SM_LEVELS)      ! INOUT Frozen soil moisture
C                                      !       content of each layer
C                                      !       as a fraction of
C                                      !       saturation.
     +,STHU(LAND_FIELD,SM_LEVELS)      ! INOUT UNfrozen soil moisture
C                                      !       content of each layer
C                                      !       as a fraction of
C                                      !       saturation.
     +,T(P_FIELD,BL_LEVELS)            ! INOUT Input:atmospheric temp(K)
C                                      !      Output:liquid/frozen water
C                                      !       temperature (TL) (K)
     +,T_DEEP_SOIL(LAND_FIELD,ST_LEVELS) ! INOUT Deep soil temperatures
C                                        !       (K).
     +,TI(P_FIELD)                     ! INOUT Sea-ice surface layer
C                                      ! temperature (K)
     +,TSTAR(P_FIELD)                  ! INOUT Surface temperature
C                                      !       (= top soil layer
C                                      !       temperature) (K).
     +,U(U_FIELD,BL_LEVELS)            ! INOUT W'ly wind component
C                                      !       (metres per second).
     +,V(U_FIELD,BL_LEVELS)            ! INOUT S'ly wind component
C                                      !       (metres per second).
     +,Z0MSEA(P_FIELD)                 ! INOUT Sea-surface roughness
C                                      !       length for momentum (m).
C                                      !       NB: same storage is used
C                                      !       for Z0V, so the intent is
C                                      !       IN for land points.
C
C  Outputs :-
C
C-1 Diagnostic (or effectively so - includes coupled model requisites):-
C
C  (a) Calculated anyway (use STASH space from higher level) :-
C
      REAL
     & BL_TYPE_1(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_2(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_3(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_4(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_5(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     &,BL_TYPE_6(P_FIELD)        ! OUT Dummy (diagnostics for BDYLYR6A)
     +,CD(P_FIELD)               ! OUT Turbulent surface exchange (bulk 
C                                !     transfer) coefficient for
C                                !     momentum.
     +,CH(P_FIELD)               ! OUT Turbulent surface exchange (bulk
C                                !     transfer) coefficient for heat
C                                !     and/or moisture.
     +,E_SEA(P_FIELD)            ! OUT Evaporation from sea times leads
C                                !     fraction. Zero over land.
C                                !     (kg per square metre per sec).
     +,EPOT(P_FIELD)             ! Dummy.
     +,ETRAN(P_FIELD)            ! OUT Transpiration (kg/m2/s).
     +,EXT(LAND_FIELD,SM_LEVELS)
C                                ! OUT Extraction of water from
C                                !     each soil layer (kg/m2/s)
     +,FQW(P_FIELD,BL_LEVELS)    ! OUT Moisture flux between layers
C                                !     (kg per square metre per sec).
C                                !     FQW(,1) is total water flux
C                                !     from surface, 'E'.
     +,FSMC(LAND_FIELD)          ! Dummy.
     +,FTL(P_FIELD,BL_LEVELS)    ! OUT FTL(,K) contains net turbulent
C                                !     sensible heat flux into layer K
C                                !     from below; so FTL(,1) is the
C                                !     surface sensible heat, H.  (W/m2)
     +,GPP(LAND_FIELD)           ! OUT Gross primary productivity
C                                !     (kg C/m2/s).
     +,H_SEA(P_FIELD)            ! OUT Surface sensible heat flux over
C                                !     sea times leads fraction. (W/m2)
     +,NPP(LAND_FIELD)           ! OUT Net primary productivity
C                               !      (kg C/m2/s).
     +,RESP_P(LAND_FIELD)        ! OUT Plant respiration (kg C/m2/s).

     +,RHOKH(P_FIELD,BL_LEVELS)  ! OUT Exchange coeffs for moisture.
C                                !     Surface:out of SF_EXCH containing
C                                !     GAMMA(1)*RHOKH,after IMPL_CAL
C                                !     contains only RHOKH.
C                                !     Above surface:out of KMKH cont-
C                                !     aining GAMMA(1)*RHOKH(,1)*RDZ(,1)
     +,RHOKM(U_FIELD,BL_LEVELS)  ! OUT Exchange coefficients for
C                                !     momentum (on UV-grid, with 1st
C                                !     and last rows undefined (or, at
C                                !     present, set to "missing data")).
C                                !     Surface:out of SF_EXCH containing
C                                !     GAMMA(1)*RHOKH,after IMPL_CAL
C                                !     contains only RHOKH.
C                                !     Above surface:out of KMKH cont-
C                                !     aining GAMMA(1)*RHOKH(,1)*RDZ(,1)
     +,RIB(P_FIELD)              ! OUT Bulk Richardson number for lowest
C                                !     layer.
     +,SEA_ICE_HTF(P_FIELD)      ! OUT Heat flux through sea-ice (W per
C                                !     sq m, positive downwards).
     +,SOIL_HT_FLUX(P_FIELD)     ! OUT Dummy
     +,SURF_HT_FLUX(P_FIELD)     ! OUT Net downward heat flux at surface
!                                !     over land or sea-ice fraction of
!                                !     gridbox (W/m2)
     +,TAUX(U_FIELD,BL_LEVELS)   ! OUT W'ly component of surface wind
C                                !     stress (N/sq m).(On UV-grid with
C                                !     first and last rows undefined or
C                                !     at present, set to 'missing data'
     +,TAUY(U_FIELD,BL_LEVELS)   ! OUT S'ly component of surface wind
C                                !     stress (N/sq m).  On UV-grid;
C                                !     comments as per TAUX.
     +,VSHR(P_FIELD)             ! OUT Magnitude of surface-to-lowest
C                                !     atm level wind shear (m per s).
C
     &,ZHT(P_FIELD)              ! OUT Dummy (diagnostics for BDYLYR6A)
     +,RHO_CD_MODV1(P_FIELD)     ! OUT Surface air density * drag coef.
C                                ! mod(v1 - v0) before interpolation.

     +,RHO_KM(P_FIELD,2:BL_LEVELS)! OUT Air density * turbulent mixing
C                                ! coef. for momentum before
     +,RHO_ARESIST(P_FIELD)       ! OUT, RHOSTAR*CD_STD*VSHR for SCYCLE
     +,ARESIST(P_FIELD)           ! OUT, 1/(CD_STD*VSHR)    for SCYCLE
     +,RESIST_B(P_FIELD)          ! OUT,(1/CH-1/CD_STD)/VSHR for SCYCLE
C
      INTEGER
     + NRML(P_FIELD)             ! OUT Number of model layers in the
C                                !     Rapidly Mixing Layer; diagnosed
C                                !     in SF_EXCH and KMKH and used in
C                                !     IMPL_CAL, SF_EVAP and TR_MIX.
C
C  (b) Not passed between lower-level routines (not in workspace at this
C      level) :-
C
      REAL
     + FME(P_FIELD)             ! OUT Wind mixing "power" (W per sq m).
     +,SICE_MLT_HTF(P_FIELD)    ! OUT Heat flux due to melting of sea-
C                               !     ice (Watts per sq metre).
     +,SNOMLT_SURF_HTF(P_FIELD)
     +,LATENT_HEAT(P_FIELD)     ! OUT Surface latent heat flux, +ve
C                               !     upwards (Watts per sq m).
     +,Q1P5M(P_FIELD)           ! OUT Q at 1.5 m (kg water per kg air).
     +,T1P5M(P_FIELD)           ! OUT T at 1.5 m (K).
     +,U10M(U_FIELD)            ! OUT U at 10 m (m per s).
     +,V10M(U_FIELD)            ! OUT V at 10 m (m per s).

 
C
C-2 Genuinely output, needed by other atmospheric routines :-
C
      REAL
     + ECAN(P_FIELD)  ! OUT Gridbox mean evaporation from canopy/surface
C                     !     store (kg per sq m per s).  Zero over sea.
     +,EI(P_FIELD)    ! OUT Sublimation from lying snow or sea-ice (kg
C                     !     per sq m per sec).
     +,ES(P_FIELD)    ! OUT Surface evapotranspiration through a
C                     !     resistance which is not entirely aerodynamic
C                     !     i.e. "soil evaporation".  Always non-
C                     !     negative.  Kg per sq m per sec.
     +,SNOWMELT(P_FIELD) ! OUT Snowmelt (kg/m/s)
     +,ZH(P_FIELD)    ! OUT Height above surface of top of boundary
C                     !     layer (metres).
     &,T1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
C                     !     of layer 1 temperature; for use in
C                     !     initiating convection.
     &,Q1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
C                     !     of layer 1 humidity; for use in initiating
C                     !     convection.
! Additional arguments for 7A boundary layer (MOSES II)
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
     & TILE_INDEX(NTYPE)
     &,TILE_PTS(NTYPE)
     &,LAND_FIELD_TRIF
     &,NPFT_TRIF
      REAL
     & ARESIST_TILE(NTYPE)
     &,CANHT_FT(NPFT)
     &,CANOPY_TILE(NTYPE-1)
     &,CATCH_TILE(NTYPE-1)
     &,CS
     &,ECAN_TILE(NTYPE-1)
     &,ESOIL_TILE(NTYPE-1)
     &,FRAC(NTYPE)
     &,FTL_TILE(NTYPE)
     &,G_LEAF(NPFT)
     &,G_LEAF_ACC(NPFT)
     &,GPP_FT(NPFT)
     &,LAI_FT(NPFT)
     &,NPP_FT(NPFT)
     &,NPP_FT_ACC(NPFT)
     &,RAD_NO_SNOW(P_FIELD)
     &,RAD_SNOW(P_FIELD)
     &,RESIST_B_TILE(NTYPE)
     &,RESP_P_FT(NPFT)
     &,RESP_S
     &,RESP_S_ACC
     &,RESP_W_FT(NPFT) 
     &,RESP_W_FT_ACC(NPFT)
     &,RHO_ARESIST_TILE(NTYPE)
     &,RIB_TILE(NTYPE)
     &,SNOW_FRAC
     &,SNOW_SURF_HTF
     &,SOIL_SURF_HTF
     &,TILE_FRAC(NPFT)
     &,TSNOW
     &,TSTAR_TILE(NTYPE)
     &,Z0V_TILE(NTYPE)
      LOGICAL
     & L_PHENOL
     &,L_TRIFFID
     &,L_NEG_TSTAR

      INTEGER
     + ERROR          ! OUT 0 - AOK;
     +,I
C                     !     1 to 7  - bad grid definition detected;
C                     !     11 - error in SF_EXCH;
C                     !     21 - error in KMKH;
C                     !     31 - error in IMPL_CAL;
C                     !     41 - error in SOIL_HTF.
C*----------------------------------------------------------------------

      CALL BDY_LAYR (

C IN values defining field dimensions and subset to be processed :
     + P_FIELD,U_FIELD,LAND_FIELD,P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH

C IN values defining vertical grid of model atmosphere :
     +,BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK
     +,EXNER

C IN soil/vegetation/land surface data :
     +,LAND_MASK,GATHER,LAND_INDEX,ST_LEVELS
     +,SM_LEVELS,CANOPY,CATCH,HCAP,HCON,LAYER_DEPTH
     +,LYING_SNOW,RESIST,ROOTD,SMC,SMVCCL,SMVCWT,Z0V,SIL_OROG_LAND
     +,L_Z0_OROG,HO2R2_OROG

C IN sea/sea-ice data :
     +,DI,ICE_FRACT,U_0,V_0

C IN cloud data :
     +,CF,QCF,QCL
     +,CCA,CCB,CCT

C IN everything not covered so far :
     +,PSTAR,RADNET,TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN

C INOUT data :
     +,Q,T,T_DEEP_SOIL,TSTAR,U,V,Z0MSEA

C OUT Diagnostic not requiring STASH flags :
     +,CD,CH,E_SEA,FQW,FTL,H_SEA,RHOKH,RHOKM,RIB
     +,SEA_ICE_HTF,SOIL_HT_FLUX,TAUX,TAUY,VSHR

C OUT diagnostic requiring STASH flags :
     +,FME,SICE_MLT_HTF,LATENT_HEAT,Q1P5M,T1P5M,U10M,V10M
C (IN) STASH flags :-
     +,SFME,SIMLT,SLH,SQ1P5,ST1P5,SU10,SV10

C OUT data required for tracer mixing :
     &,NRML

C OUT data required for 4D-VAR :
     &,RHO_CD_MODV1,RHO_KM

C OUT data required elsewhere in UM system :
     +,ECAN,EI,ES,ZH,T1_SD,Q1_SD,ERROR

C LOGICAL LTIMER
     +,LTIMER
     +)

      RETURN
      END

