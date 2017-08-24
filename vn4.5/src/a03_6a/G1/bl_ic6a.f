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
!!!  4.4    10/09/97  New deck.  R.N.B.Smith
!!!  4.5   24/04/98   New diagnostics ZHT and BL_TYPE_1 to _6.
!!!                   R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!  4.5  27/04/98 Add dummy arguments LAND_FIELD_TRIF, NPFT_TRIF,
!!!                GPP_FT, RESP_P_FT, TILE_FRAC_OUT, L_PHENOL and
!!!                L_TRIFFID and change TILE_INDEX and TILE_PTS from
!!!                IN to OUT, consistent with changes to BL_CTL
!!!                required for TRIFFID.   Richard Betts                
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!! System components covered : P24
!!!
!!! System task : P0
!!!
!!!END -----------------------------------------------------------------
!    Arguments :-
      SUBROUTINE BL_INTCT(

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
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,TIMESTEP,
     & L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN,

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

! Additional arguments for A03_7A (MOSES-II) - not used in A03_6A.
! IN
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,
     & CANHT_FT,CANOPY_TILE,CATCH_TILE,CS,LAI_FT,
     & FRAC,SNOW_FRAC,RAD_NO_SNOW,RAD_SNOW,TSNOW,Z0V_TILE,
     & CO2_3D,CO2_DIM,L_CO2_INTERACTIVE,
! INOUT
     & TSTAR_TILE_DUMMY,
     & G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,
! OUT
     & ECAN_TILE,ESOIL_TILE,FTL_TILE_DUMMY,
     & G_LEAF,GPP_FT,NPP_FT,RESP_P_FT,RESP_S,RESP_W_FT,
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,
     & RIB_TILE,SNOW_SURF_HTF,SOIL_SURF_HTF,
     & TILE_INDEX,TILE_PTS,TILE_FRAC_OUT,

! LOGICAL LTIMER
     & LTIMER
     &)
      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER
     & P_FIELD                   ! IN No. of P-points in whole grid
!                                !    (for dimensioning only).
     &,RADHR_DIM1                ! IN Dimension of Radiative heating
!                                !    rate (P_FIELD but used for
!                                !    dynamic allocation)
     &,U_FIELD                   ! IN No. of UV-points in whole grid.
!                                !    (Checked for consistency with
!                                !    P_FIELD and P_ROWS; there must
!                                !    be 1 less UV than P row.)
     &,LAND_FIELD                ! IN No.of land points in whole grid.
!                                !    (Checked for consistency with
!                                !    P_FIELD )
     &,P_ROWS                    ! IN No. of P-rows in whole grid
!                                !    (for dimensioning only).
     &,FIRST_ROW                 ! IN First row of data to be treated,
!                                !    referred to P-grid (must be > 1
!                                !    since "polar" rows are never
!                                !    treated).
     &,N_ROWS                    ! IN No. of rows of data to be
!                                !    treated, referred to P-grid.
!                                !    FIRST_ROW+N_ROWS-1 must be less
!                                !    than P_ROWS, since "polar" rows
!                                !    are never treated.
     &,ROW_LENGTH                ! IN No. of points in one row.
!                                !    (Checked for consistency with
!                                !    P_FIELD and N_ROWS.)

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                 ! IN Max. no. of "boundary" levels
!                                !    allowed.Assumed <= 30 for dim-
!                                !    sioning of GAMMA in common deck
!                                !    C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                  ! IN Total no. of vertical levels in
!                                !    the model atmosphere.
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
     &,L_CO2_INTERACTIVE   ! dummy CO2 variable - not used with A03_6A.
     &,GATHER                    ! IN T if gather to sea-ice points
!                                !    in SF_EXCH. Saves a lot of un-
!                                !    necessary calculations if there
!                                !    are relatively few sea-ice points
     &,L_RMBL                    ! IN T to use rapidly mixing
!                                !    boundary scheme in IMPL_CAL
     &,L_BL_LSPICE               ! IN True if 3A large-scale ppn
!                                !    scheme is used.
     &,L_MOM                     ! IN Switch for convective momentum
!                                !    transport.
     &,L_MIXLEN                  ! IN Switch for reducing the turbulent
!                                !    mixing length above the top of the
!                                !    boundary layer - dummy for A03_5B
     &,L_Z0_OROG                 ! IN T to use simple orog.roughness
!                                !    treatment in SF_EXCH
      INTEGER
     & LAND_INDEX(P_FIELD)       ! IN LAND_INDEX(I)=J => the Jth
!                                !    point in P_FIELD is the Ith
!                                !    land point.
      INTEGER
     & ST_LEVELS                 ! IN No. of deep soil temp. levels
     &,SM_LEVELS                 ! IN No. of soil moisture levels
     &,CO2_DIM       ! dummy CO2 variable - field not used with A03_6A.

      REAL
     & CANHT(LAND_FIELD)         ! IN Canopy height (m)
     &,CANOPY(LAND_FIELD)        ! IN Surface/canopy water (kg per sq m)
     &,CATCH(LAND_FIELD)         ! IN Surface/canopy water capacity
!                                !    (kg per sq m).
     &,HCAP(LAND_FIELD)          ! IN Soil heat capacity (J/K/m**3)
     &,HCON(LAND_FIELD)          ! IN Soil thermal conductivity (W/m/K).
     &,HO2R2_OROG(LAND_FIELD)    ! IN Dummy used only in version 3A.

     &,LAI(LAND_FIELD)           ! IN Leaf area index
     &,LAYER_DEPTH(SM_LEVELS)    !    Dummy Variable not used
!                                !    in this version
     &,LYING_SNOW(P_FIELD)       ! IN Lying snow (kg per sq m).
!                                !    Must be global for coupled model,
!                                !    ie dimension P_FIELD not LAND_FIEL
     &,RESIST(LAND_FIELD)        ! IN "Stomatal" resistance to
!                                !    evaporation (seconds per metre).
     &,ROOTD(LAND_FIELD)         ! IN Depth of active soil layer ("root
!                                !    depth") (metres).
     &,SIL_OROG_LAND(LAND_FIELD) ! IN Silhouette area of unresolved
!                                     orography per unit horizontal area
!                                     on land points only.
     &,SMC(LAND_FIELD)           ! IN Soil moisture content (kg / sq m).
     &,SMVCCL(LAND_FIELD)        ! IN Critical volumetric SMC (cubic m
!                                !    per cubic m of soil).
     &,SMVCST(LAND_FIELD)        ! IN Volumetric saturation point (cubic
!                                !    per cubic m of soil).
     &,SMVCWT(LAND_FIELD)        ! IN Volumetric wilting point (cubic m
!                                !    per cubic m of soil).
     &,TILE_FRAC(P_FIELD)        ! IN Dummy
     &,VFRAC(LAND_FIELD)         ! IN Vegetation fraction.
     &,Z0V(P_FIELD)              ! IN Vegetative roughness length (m).
!                                !    NB:UM uses same storage for Z0MSEA
!                                !    so for sea points this is INOUT.
!
! Additional arrays for A03_7A (MOSES-II) - not used in A03_6A.
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
     &,FTL_TILE_DUMMY(NTYPE)
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
     &,TILE_FRAC_OUT(NPFT)
     &,TSNOW
     &,TSTAR_TILE_DUMMY(NTYPE)
     &,Z0V_TILE(NTYPE)
      LOGICAL
     & L_PHENOL
     &,L_TRIFFID
     &,L_NEG_TSTAR

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
     &,CO2_3D        ! dummy CO2 variable - field not used with A03_6A.
     &,PHOTOSYNTH_ACT_RAD(P_FIELD) ! IN Net downward shortwave radiation
!                                !    in band 1 (w/m2).
     &,PSTAR(P_FIELD)            ! IN Surface pressure (Pascals).
     &,RADNET(P_FIELD)           ! IN Surface net radiation (W/sq m,
!                                !    positive downwards).
     &,RAD_HR(RADHR_DIM1,BL_LEVELS)
!                                ! IN Radiative heating rate (K/s).
     &,TIMESTEP                  ! IN Timestep (seconds).

      LOGICAL LTIMER             ! Logical switch for TIMER diags

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
     &,STHF(LAND_FIELD,SM_LEVELS)! INOUT Frozen soil moisture
!                                !       content of each layer
!                                !       as a fraction of saturation.
     &,STHU(LAND_FIELD,SM_LEVELS)! INOUT UNfrozen soil moisture
!                                !       content of each layer
!                                !       as a fraction of saturation.
     &,T(P_FIELD,BL_LEVELS)      ! INOUT Input:atmospheric temp(K)
!                                !      Output:liquid/frozen water
!                                !       temperature (TL) (K)
     &,T_DEEP_SOIL(LAND_FIELD,ST_LEVELS)
!                                ! INOUT Deep soil temperatures (K).
     &,TI(P_FIELD)               ! INOUT Sea-ice surface layer
!                                ! temperature (K)
     &,TSTAR(P_FIELD)            ! INOUT Surface temperature
!                                !       (= top soil layer
!                                !       temperature) (K).
     &,TSTAR_TILE(P_FIELD)       ! INOUT Dummy
     &,U(U_FIELD,BL_LEVELS)      ! INOUT W'ly wind component (m/s).
     &,V(U_FIELD,BL_LEVELS)      ! INOUT S'ly wind component (m/s).
     &,Z0MSEA(P_FIELD)           ! INOUT Sea-surface roughness
!                                !       length for momentum (m).
!                                !       NB: same storage is used
!                                !       for Z0V, so the intent is
!                                !       IN for land points.

!  Outputs :-

!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-

      REAL
     & BL_TYPE_1(P_FIELD)        ! OUT Indicator set to 1.0 if stable
!                                !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_2(P_FIELD)        ! OUT Indicator set to 1.0 if Sc over
!                                !     stable surface layer diagnosed,
!                                !     0.0 otherwise.
     &,BL_TYPE_3(P_FIELD)        ! OUT Indicator set to 1.0 if well 
!                                !     mixed b.l. diagnosed,
!                                !     0.0 otherwise.
     &,BL_TYPE_4(P_FIELD)        ! OUT Indicator set to 1.0 if 
!                                !     decoupled Sc layer (not over
!                                !     cumulus) diagnosed,
!                                !     0.0 otherwise.
     &,BL_TYPE_5(P_FIELD)        ! OUT Indicator set to 1.0 if
!                                !     decoupled Sc layer over cumulus
!                                      diagnosed, 0.0 otherwise.
     &,BL_TYPE_6(P_FIELD)        ! OUT Indicator set to 1.0 if a 
!                                !     cumulus capped b.l. diagnosed,
!                                !     0.0 otherwise.                   
     &,CD(P_FIELD)              ! OUT Turbulent surface exchange (bulk  
!                               !     transfer) coefficient for
!                               !     momentum.
     &,CH(P_FIELD)              ! OUT Turbulent surface exchange (bulk
!                               !     transfer) coefficient for heat
!                               !     and/or moisture.
     &,E_SEA(P_FIELD)           ! OUT Evaporation from sea times leads
!                               !     fraction. Zero over land.
!                               !     (kg per square metre per sec).
     &,EPOT(P_FIELD)            ! OUT potential evaporation-rate
!                               !     (kg/m2/s).
     &,ETRAN(P_FIELD)           ! OUT Transpiration (kg/m2/s).
     &,EXT(LAND_FIELD,SM_LEVELS)! OUT Extraction of water from
!                               !     each soil layer (kg/m2/s)
     &,FQW(P_FIELD,BL_LEVELS)   ! OUT Moisture flux between layers
!                               !     (kg per square metre per sec).
!                               !     FQW(,1) is total water flux
!                               !     from surface, 'E'.
     &,FSMC(LAND_FIELD)         ! OUT soil moisture availability.
     &,FTL_TILE(P_FIELD)        ! OUT Dummy
     &,FTL(P_FIELD,BL_LEVELS)   ! OUT FTL(,K) contains net turbulent
!                               !     sensible heat flux into layer K
!                               !     from below; so FTL(,1) is the
!                               !     surface sensible heat, H.  (W/m2)
     &,FQW_TILE(P_FIELD)        ! OUT Dummy
     &,GPP(LAND_FIELD)          ! OUT Gross primary productivity
!                               !     (kg C/m2/s).
     &,H_SEA(P_FIELD)           ! OUT Surface sensible heat flux over
!                               !     sea times leads fraction. (W/m2)
     &,NPP(LAND_FIELD)          ! OUT Net primary productivity
!                               !      (kg C/m2/s).
     &,RESP_P(LAND_FIELD)       ! OUT Plant respiration (kg C/m2/s).
     &,RHOKH(P_FIELD,BL_LEVELS) ! OUT Exchange coeffs for moisture.
     &,RHOKM(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                               !     momentum (on UV-grid, with 1st
!                               !     and last rows undefined (or, at
!                               !     present, set to "missing data")).
     &,RIB(P_FIELD)             ! OUT Bulk Richardson number for lowest
!                               !     layer.
     &,RIB_GB(P_FIELD)          ! OUT Dummy
     &,SEA_ICE_HTF(P_FIELD)     ! OUT Heat flux through sea-ice (W per
!                               !     sq m, positive downwards).
     &,SOIL_HT_FLUX(P_FIELD)    ! OUT Dummy
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
     &,VSHR(P_FIELD)            ! OUT Magnitude of surface-to-lowest
!                               !     atm level wind shear (m per s).
     &,ZHT(P_FIELD)              ! OUT Height below which there may be
!                                !     turbulent mixing (m).
     &,RHO_CD_MODV1(P_FIELD)    ! OUT Surface air density * drag coef.
!                               ! mod(v1 - v0) before interpolation.
     &,RHO_KM(P_FIELD,2:BL_LEVELS)! OUT Air density * turbulent mixing
!                               ! coef. for momentum before
     &,RHO_ARESIST(P_FIELD)     ! OUT, RHOSTAR*CD_STD*VSHR for SCYCLE
     &,ARESIST(P_FIELD)         ! OUT, 1/(CD_STD*VSHR)    for SCYCLE
     &,RESIST_B(P_FIELD)        ! OUT,(1/CH-1/CD_STD)/VSHR for SCYCLE
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

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL
     & ECAN(P_FIELD)  ! OUT Gridbox mean evaporation from canopy/surface
!                     !     store (kg per sq m per s).  Zero over sea.
     &,EI(P_FIELD)    ! OUT Sublimation from lying snow or sea-ice (kg
!                     !     per sq m per sec).
     &,ES(P_FIELD)    ! OUT Surface evapotranspiration through a
!                     !     resistance which is not entirely aerodynamic
!                     !     i.e. "soil evaporation".  Always non-
!                     !     negative.  Kg per sq m per sec.
     &,SNOWMELT(P_FIELD) ! OUT Snowmelt (kg/m/s)
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
     &,I
!                     !     1 to 7  - bad grid definition detected;


C-----------------------------------------------------------------------

      INTEGER N_TYPES
      PARAMETER (N_TYPES=1)
C  Number tiles per land gridbox
C-----------------------------------------------------------------------

!*----------------------------------------------------------------------


      CALL BDY_LAYR (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,N_TYPES,LAND_FIELD,
     & P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN soil/vegetation/land surface data :
     & LAND_MASK,GATHER,LAND_INDEX,
     & ST_LEVELS,SM_LEVELS,TILE_FRAC,CANHT,CANOPY,
     & CATCH,CATCH,HCON,
     & LYING_SNOW,RESIST,RESIST,ROOTD,ROOTD,
     & SMVCCL,SMVCST,SMVCWT,STHF,STHU,
     & VFRAC,Z0V,Z0V,SIL_OROG_LAND,L_Z0_OROG,HO2R2_OROG,
     & LAI,

! IN sea/sea-ice data :
     & DI,ICE_FRACT,U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,RAD_HR,RADHR_DIM1,
     & TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,

! INOUT data :
     & Q,GS,T,T_DEEP_SOIL,TI,TSTAR,TSTAR_TILE,U,V,Z0MSEA,

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,E_SEA,EPOT,ETRAN,FQW,FQW_TILE,FSMC,FTL,FTL_TILE,
     & H_SEA,RHOKH,RHOKM,
     & RIB_GB,RIB,SEA_ICE_HTF,SURF_HT_FLUX,TAUX,TAUY,VSHR,ZHT,
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,
 

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
     & ECAN,EI,ES,EXT,SNOWMELT,ZH,
     & GPP,NPP,RESP_P,
     & T1_SD,Q1_SD,ERROR,

! LOGICAL LTIMER
     & LTIMER
     &)
      RETURN
      END

