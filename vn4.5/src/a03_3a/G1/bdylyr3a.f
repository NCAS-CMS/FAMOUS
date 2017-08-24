C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C
C*LL  SUBROUTINE BDY_LAYR-----------------------------------------------
CLL
CLL  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
CLL           between (a) surface and atmosphere, (b) atmospheric levels
CLL           within the boundary layer, and/or the effects of these
CLL           fluxes on the primary model variables.  The flux of heat
CLL           into and through the soil is also modelled.  Numerous
CLL           related diagnostics are also calculated.
CLL
CLL  If the point for a single column run is a sea point then the
CLL    land mask ensures that correct calculations are done; LAND_FIELD
CLL    dimensions land point arrays and so its value is irrelevent in
CLL    this case. An error is generated if LAND_FIELD is not one in the
CLL    land point case.
CLL           F E Hewer, July 1990: removed call to LS_CLD.
CLL    This version passes out liquid/frozen water temperature in
CLL    array "T" (TL), and total water content in array "Q" (QW).
CLL    These may be converted to T and Q respectively by calling
CLL    the large scale cloud routine, LS_CLD.
CLL            F E Hewer, August 1990: land point data stored
CLL    on land points only (dimension: LAND_FIELD, arrays:CANOPY, CATCH
CLL    HCAP, HCON, RESIST, ROOTDEP, SMC, SMVCCL, SMVCWT, T_DEEP_SOIL)
CLL    Arrays whose elements may contain values over both sea and land
CLL    points are compressed onto land points for land calculations if
CLL    defined variable IBM is NOT selected. RHOKM,RHOKH redefined as
CLL    workspace.
CLL
CLL  Suitable for single column use
CLL
CLL F.Hewer     <- programmer of some or all of previous code or changes
CLL C.Wilson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  05/01/95   rhostar*cD*modv1 and rho*Km output via argument
CLL                   list for STASHing.                R.N.B.Smith
CLL
CLL   4.0  30/12/94   References to z0h_eff removed for reformulated
CLL                   surface transfer coefficients for scalars.
CLL                                                   R.N.B.Smith
!     4.1   23/05/96   MPP code: Change updateable area  P.Burton
CLL   4.1   5/6/96    DS_LEVELS chnaged to ST_LEVELS and SM_LEVELS
CLL                   C.Bunton
CLL
CLL   4.1  08/05/96   Logical L_RMBL included to select rapidly mixing
CLL                   boundary layer scheme        Simon Jackson
CLL
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL
CLL   4.3  04/02/97   Logical switches L_MOM and L_MIXLEN passed down
CLL                   to KHKH and thence EXCOEF.
CLL                                                     R.N.B.Smith
CLL   4.4  08/09/97   L_BL_LSPICE specifies mixed phase precipitation
CLL                   scheme                     D.Wilson 
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version ?, dated ?.
CLL
CLL  System component covered: P24.
CLL
CLL  Project task:
CLL
CLL  Documentation: UMDP 24.
CLL
CLL---------------------------------------------------------------------
C*
C*L---------------------------------------------------------------------
C    Arguments :-
      SUBROUTINE BDY_LAYR (

C IN values defining field dimensions and subset to be processed :
     + P_FIELD,U_FIELD,LAND_FIELD,P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH

C IN values defining vertical grid of model atmosphere :
     +,BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK
     +,EXNER

C IN soil/vegetation/land surface data :
     +,LAND_MASK,GATHER,LAND_INDEX
     +,ST_LEVELS,SM_LEVELS,CANOPY,CATCH,HCAP,HCON,LAYER_DEPTH
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
     +,SFME,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10

C OUT data required for tracer mixing :
     &,NRML

C OUT data required for 4D-VAR :
     &,RHO_CD_MODV1,RHO_KM

C OUT data required elsewhere in UM system :
     +,ECAN,EI,ES,ZH,T1_SD,Q1_SD,ERROR

C LOGICAL LTIMER
     +,LTIMER
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
     +,GATHER                    ! IN T if gather to sea-ice points
C                                !    in SF_EXCH. Saves a lot of un-
C                                !    necessary calculations if there
C                                !    are relatively few sea-ice points
     +,L_Z0_OROG                 ! IN T to use orog.roughness
C                                !    treatment in SF_EXCH
     +,L_RMBL                    ! IN T to use rapidly mixing boundary
C                                !    scheme in IMPL_CAL
     &,L_BL_LSPICE           ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
     &,L_MOM                     ! IN Switch for convective momentum
C                                !    transport.
     &,L_MIXLEN                  ! IN Switch for reducing the turbulent
C                                !    mixing length above the top of the
C                                !    boundary layer.
C
      INTEGER
     + LAND_INDEX(P_FIELD)       ! IN LAND_INDEX(I)=J => the Jth
C                                !    point in P_FIELD is the Ith
C                                !    land point.
      INTEGER
     + ST_LEVELS                 ! IN No. of deep soil levels
     +,SM_LEVELS                 ! IN No. of soil moisture levels
C                                !
      REAL
     + CANOPY(LAND_FIELD)        ! IN Surface/canopy water (kg per sq m)
     +,CATCH(LAND_FIELD)         ! IN Surface/canopy water capacity
C                                !    (kg per sq m).
     +,HCAP(LAND_FIELD)          ! IN Soil heat capacity (J / K / m**3).
     +,HCON(LAND_FIELD)          ! IN Soil thermal conductivity (W/m/K).
     +,LAYER_DEPTH(ST_LEVELS+1) ! IN Depths of soil layers, as multiple
C                                !    of depth of top layer, counting
C                                !    from top down and including first
C                                !    layer (i.e. 1st value must be 1).
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
     +,SMVCWT(LAND_FIELD)        ! IN Volumetric wilting point (cubic m
C                                !    per cubic m of soil).
     +,Z0V(P_FIELD)              ! IN Vegetative roughness length (m).
C                                !    NB:UM uses same storage for Z0MSEA
C                                !    so for sea points this is INOUT.
     +,SIL_OROG_LAND(LAND_FIELD) ! IN Silhouette area of unresolved
C                                !    orography per unit horizontal area
C                                !    on land points only.
     +,HO2R2_OROG(LAND_FIELD)    ! IN Standard Deviation of orography.
C                                !    equivilent to peak to trough
C                                !    height of unresolved orography
C                                !    devided by 2SQRT(2) on land
C                                !    points only (m)
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
     + PSTAR(P_FIELD)            ! IN Surface pressure (Pascals).
     +,RADNET(P_FIELD)           ! IN Surface net radiation (W/sq m,
C                                !    positive downwards).
     +,TIMESTEP                  ! IN Timestep (seconds).
C
      LOGICAL LTIMER             ! Logical switch for TIMER diags
C
C  STASH flags :-
C
      LOGICAL
     + SFME    ! IN Flag for FME (q.v.).
     +,SMLT    ! IN Flag for SICE_MLT_HTF (q.v.)
     +,SLH     ! IN Flag for LATENT_HEAT (q.v.)
     +,SQ1P5   ! IN Flag for Q1P5M (q.v.)
     +,ST1P5   ! IN Flag for T1P5M (q.v.)
     +,SU10    ! IN Flag for U10M (q.v.)
     +,SV10    ! IN Flag for V10M (q.v.)
C
C  In/outs :-
C
      REAL
     + Q(P_FIELD,BL_LEVELS)            ! INOUT Input:specific humidity
C                                      !       ( kg water per kg air).
C                                      !      Output:total water content
C                                      !      (Q)(kg water per kg air).
     +,T(P_FIELD,BL_LEVELS)            ! INOUT Input:atmospheric temp(K)
C                                      !      Output:liquid/frozen water
C                                      !       temperature (TL) (K)
     +,T_DEEP_SOIL(LAND_FIELD,ST_LEVELS) ! INOUT Deep soil temperatures
C                                        !       (K).
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
     + CD(P_FIELD)               ! OUT Turbulent surface exchange (bulk
C                                !     transfer) coefficient for
C                                !     momentum.
     +,CH(P_FIELD)               ! OUT Turbulent surface exchange (bulk
C                                !     transfer) coefficient for heat
C                                !     and/or moisture.
     +,E_SEA(P_FIELD)            ! OUT Evaporation from sea times leads
C                                !     fraction. Zero over land.
C                                !     (kg per square metre per sec).
     +,FQW(P_FIELD,BL_LEVELS)    ! OUT Moisture flux between layers
C                                !     (kg per square metre per sec).
C                                !     FQW(,1) is total water flux
C                                !     from surface, 'E'.
     +,FTL(P_FIELD,BL_LEVELS)    ! OUT FTL(,K) contains net turbulent
C                                !     sensible heat flux into layer K
C                                !     from below; so FTL(,1) is the
C                                !     surface sensible heat, H.  (W/m2)
     +,H_SEA(P_FIELD)            ! OUT Surface sensible heat flux over
C                                !     sea times leads fraction. (W/m2)
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
     +,SOIL_HT_FLUX(P_FIELD)     ! OUT Heat flux from soil layer 1 to
C                                !     soil layer 2, i.e. from surface
C                                !     to deep soil layer 1, i.e. +ve
C                                !     downwards (Watts per sq m).
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
     &,RHO_CD_MODV1(P_FIELD)     ! OUT Surface air density * drag coef.*
C                                !     mod(v1 - v0) before interpolation
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT Air density * turbulent mixing
C                                !     coefficient for momentum before
C                                      interpolation.
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
     +,ZH(P_FIELD)    ! OUT Height above surface of top of boundary
C                     !     layer (metres).
     &,T1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
C                     !     of layer 1 temperature; for use in
C                     !     initiating convection.
     &,Q1_SD(P_FIELD) ! OUT Standard deviation of turbulent fluctuations
C                     !     of layer 1 humidity; for use in initiating
C                     !     convection.
      INTEGER
     + ERROR          ! OUT 0 - AOK;
C                     !     1 to 7  - bad grid definition detected;
C                     !     11 - error in SF_EXCH;
C                     !     21 - error in KMKH;
C                     !     31 - error in IMPL_CAL;
C                     !     41 - error in SOIL_HTF.
C*----------------------------------------------------------------------
C*L---------------------------------------------------------------------
C  External routines called :-
C
      EXTERNAL Z,SICE_HTF,SOIL_HTF,SF_EXCH,KMKH,IMPL_CAL,SF_EVAP
      EXTERNAL TIMER
      EXTERNAL UV_TO_P
C*----------------------------------------------------------------------
C*L---------------------------------------------------------------------
C   Symbolic constants (parameters) reqd in top-level routine :-
C
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*----------------------------------------------------------------------
C
C  Workspace :-
C
      REAL
     + DZL(P_FIELD,BL_LEVELS)     ! WORK DZL(,K) is depth in m of layer
C                                 !      K, i.e. distance from boundary
C                                 !      K-1/2 to boundary K+1/2.
     +,FQW_LEAD(P_FIELD)          ! WORK FQW for leads fractn of grid sq
     +,FTL_LEAD(P_FIELD)          ! WORK FTL for leads fractn of grid sq
     +,QW(P_FIELD,BL_LEVELS)      ! WORK Total water content, but
C                                 !      replaced by specific humidity
C                                 !      in LS_CLD.
     +,RDZ(P_FIELD,BL_LEVELS)     ! WORK RDZ(,1) is the reciprocal
C                                 !      of the height of level 1, i.e.
C                                 !      of the middle of layer 1.  For
C                                 !      K > 1, RDZ(,K) is the
C                                 !      reciprocal of the vertical
C                                 !      distance from level K-1 to
C                                 !      level K.
     +,TL(P_FIELD,BL_LEVELS)      ! WORK Ice/liquid water temperature,
C                                 !      but replaced by T in LS_CLD.
     +,TSTAR_NL(P_FIELD)          ! WORK TSTAR No Leads (K).
     +,TV_RDZUV(P_FIELD,BL_LEVELS) ! WORK Virtual temp at start (TV).
C                                  !      RDZ (K > 1) on UV-grid.
C                                  !      Comments as per RHOKM (RDZUV).
     +,U_P(P_FIELD,BL_LEVELS)     ! WORK U on P-grid.
     +,V_P(P_FIELD,BL_LEVELS)     ! WORK V on P-grid.
     +,ZLB(P_FIELD,0:BL_LEVELS)   ! WORK ZLB(,K) is the height of the
C                                 !      upper boundary of layer K
C                                 !      ( = 0.0 for "K=0").
      REAL
     + ASOIL_1(P_FIELD)           ! WORK Soil thermodynamic coefficient.
     +,BQ_1(P_FIELD)              ! WORK A buoyancy parameter for the
C                                 !      lowest atmospheric level.
     +,BT_1(P_FIELD)              ! WORK A buoyancy parameter for the
C                                 !      lowest atmospheric level.
     &,BF_1(P_FIELD)
!        WORK   A bouyancy parameter for the lowest atmospheric level
     +,DQW_1(P_FIELD)             ! WORK Increment for QW(,1).
     +,DQW_RML(P_FIELD)           ! WORK Increment for QW in the
C                                 !      Rapidly Mixing Layer.
     +,DTRDZ(P_FIELD,BL_LEVELS)   ! WORK -g.dt/dp for model layers.
     +,DTRDZ_RML(P_FIELD)         ! WORK -g.dt/dp for the rapidly
C                                 !      mixing layer.
     +,DTSTAR(P_FIELD)            ! WORK Increment for TSTAR.
     +,EA(P_FIELD)                ! WORK Surface evaporation with only
C                                 !      aerodynamic resistance (+ve),
C                                 !      or condensation (-ve),
C                                 !      averaged over gridbox.
     +,ESL(P_FIELD)               ! WORK ES (q.v.) without fractional
C                                 !      weighting factor FRACS.  "L"
C                                 !      is for "local".
     +,RHOKE(P_FIELD)             ! WORK Surface exchange coefficient
C                                 !      for FQW.
     +,RHOKEA(P_FIELD)            ! WORK Surface exchange coefficient
C                                 !      for EA.
     +,RHOKES(P_FIELD)            ! WORK Surface exchange coefficient
C                                 !      for ES.
     +,RHOKESL(P_FIELD)           ! WORK Surface exchange coefficient
C                                 !      for ESL.
     +,Z0H(P_FIELD)               ! WORK Roughness length for heat and
C                                 !      moisture.
     +,Z0M(P_FIELD)               ! WORK Roughness length for momentum.
     +,Z1(P_FIELD)                ! WORK Height of lowest level (i.e.
C                                 !      height of middle of lowest
C                                 !      layer).
     +,H_BLEND(P_FIELD)           ! WORK Blending height used as part of
C                                 !      effective roughness scheme
     +,Z0M_EFF(P_FIELD)           ! WORK Effective roughness length for
C                                 !      momentum
      REAL
     + CDR10M(U_FIELD)            ! WORK Ratio of CD's reqd for
C                                 !      calculation of 10 m wind.
C                                 !      On UV-grid; comments as per
C                                 !      RHOKM.
     +,CER1P5M(P_FIELD)           ! WORK Ratio of coefficients reqd for
C                                 !      calculation of 1.5 m Q.
     +,CHR1P5M(P_FIELD)           ! WORK Ratio of coefficients reqd for
C                                 !      calculation of 1.5 m T.
C
C  Local scalars :-
C
      INTEGER
     + ERR        ! LOCAL Return codes from lower-level routines.
     +,I          ! LOCAL Loop counter (horizontal field index).
     +,K          ! LOCAL Loop counter (vertical level index).
     +,N_P_ROWS   ! LOCAL No of P-rows being processed.
     +,N_U_ROWS   ! LOCAL No of UV-rows being processed.
     +,P_POINTS   ! LOCAL No of P-points being processed.
     +,P1         ! LOCAL First P-point to be processed.
     +,LAND1      ! LOCAL First land-point to be processed.
C                 !           1 <= LAND1 <= LAND_FIELD
     +,LAND_PTS   ! LOCAL No of land points being processed.
     +,SEA_FIELD  ! LOCAL No of sea points in the full field.
     +,SEA_PTS    ! LOCAL No of sea points being processed.
     +,U_POINTS   ! LOCAL No of UV-points being processed.
     +,U1         ! LOCAL First UV-point to be processed.
      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF
      ERROR = 0
C
C-----------------------------------------------------------------------
CL 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
CL    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
CL    most UV-rows are used only for interpolation and are not updated.
C-----------------------------------------------------------------------
C
      IF ( BL_LEVELS.LT.1 .OR. ST_LEVELS.LT.1 
     &.OR. P_ROWS.LT.3 ) THEN
        ERROR = 1
        GOTO999
      ELSEIF ( U_FIELD .NE. (P_ROWS)*ROW_LENGTH ) THEN
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
C
C-----------------------------------------------------------------------
CL    Set pointers, etc.
C-----------------------------------------------------------------------
C
      N_P_ROWS=N_ROWS
      N_U_ROWS=N_ROWS+1

      P_POINTS=N_P_ROWS*ROW_LENGTH
      U_POINTS=N_U_ROWS*ROW_LENGTH

      P1=1+(FIRST_ROW-1)*ROW_LENGTH
      U1=1+(FIRST_ROW-2)*ROW_LENGTH
C
C-----------------------------------------------------------------------
CL    Set compressed land point pointers.
C-----------------------------------------------------------------------
C
      LAND1=0
      DO 1 I=1,P1+P_POINTS-1
        IF (LAND_INDEX(I).GE.P1) THEN
          LAND1 = I
          GOTO2
        ENDIF
   1  CONTINUE
   2  CONTINUE
      LAND_PTS=0
      DO 3 I=P1,P1+P_POINTS-1
        IF (LAND_MASK(I)) LAND_PTS = LAND_PTS + 1
   3  CONTINUE
C
C-----------------------------------------------------------------------
CL 1.  Perform calculations in what the documentation describes as
CL     subroutine Z_DZ.  In fact, a separate subroutine isn't used.
C-----------------------------------------------------------------------
C
CL 1.1 Initialise ZLB(,0) (to zero, of course, this being the height
CL     of the surface above the surface).
C
      DO 4 I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
    4 CONTINUE
C
CL 1.2 Calculate layer depths and heights, and construct wind fields on
CL     P-grid.  This involves calling subroutines Z and UV_TO_P.
CL     Virtual temperature is also calculated, as a by-product.
C
C NB RDZ  TEMPORARILY used to return DELTA_Z_LOWER, the lower half layer
C    thickness
C
      DO 5 K=1,BL_LEVELS
        CALL Z(P_POINTS,EXNER(P1,K),EXNER(P1,K+1),PSTAR(P1),
     +    AKH(K),BKH(K),Q(P1,K),QCF(P1,K),
     +    QCL(P1,K),T(P1,K),ZLB(P1,K-1),TV_RDZUV(P1,K),
     +    ZLB(P1,K),DZL(P1,K),RDZ(P1,K),LTIMER)
        CALL UV_TO_P(U(U1,K),U_P(P1,K),
     +               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V(U1,K),V_P(P1,K),
     +               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
    5 CONTINUE
      DO 61 K=BL_LEVELS,2,-1
        DO 62 I=P1,P1+P_POINTS-1
          RDZ(I,K)=1.0/(RDZ(I,K)+(DZL(I,K-1)-RDZ(I,K-1)))
   62   CONTINUE
   61 CONTINUE
      DO 6 I=P1,P1+P_POINTS-1
        Z1(I)=RDZ(I,1)
        RDZ(I,1)=1.0/RDZ(I,1)
    6 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.  Heat flux through sea-ice (P241, routine SICE_HTF).
C-----------------------------------------------------------------------
C
      CALL SICE_HTF(
     + DI(P1),ICE_FRACT(P1),LAND_MASK(P1),TSTAR(P1),P_POINTS,
     + SEA_ICE_HTF(P1),LTIMER
     +)
C
C-----------------------------------------------------------------------
CL 3.  Heat fluxes through soil (P242, routine SOIL_HTF).
C-----------------------------------------------------------------------
C
      CALL SOIL_HTF(
     + HCAP,HCON,LAYER_DEPTH,LYING_SNOW,TSTAR,LAND_MASK,TIMESTEP,
     + P_FIELD,LAND_FIELD,P_POINTS,P1,
     + LAND_PTS,LAND1,LAND_INDEX,
     + ST_LEVELS+1,T_DEEP_SOIL,ASOIL_1,SOIL_HT_FLUX
     +,LTIMER)
C
C-----------------------------------------------------------------------
CL 4.  Surface turbulent exchange coefficients and "explicit" fluxes
CL     (P243a, routine SF_EXCH).
CL     Wind mixing "power" and some values required for other, later,
CL     diagnostic calculations, are also evaluated if requested.
C-----------------------------------------------------------------------
C
      CALL SF_EXCH (
     + P_POINTS,LAND_PTS,U_POINTS,ROW_LENGTH,N_P_ROWS,N_U_ROWS,
     + LAND_INDEX(LAND1),P1,GATHER,
     + AK(1),BK(1),TIMESTEP,CANOPY(LAND1),CATCH(LAND1),CF(P1,1),
     + ICE_FRACT(P1),LAND_MASK(P1),PSTAR(P1),Q(P1,1),QCF(P1,1),
     + QCL(P1,1),RESIST(LAND1),ROOTD(LAND1),SMC(LAND1),
     + SMVCCL(LAND1),SMVCWT(LAND1),LYING_SNOW(P1),T(P1,1),TSTAR(P1),
     + U(U1,1),V(U1,1),U_P(P1,1),V_P(P1,1),U_0(U1),V_0(U1),
     + Z0V(P1),SIL_OROG_LAND(LAND1),Z1(P1),Z0MSEA(P1),HO2R2_OROG(LAND1),
     & BQ_1(P1),BT_1(P1),BF_1(P1),CD(P1),CH(P1),
     + EA(P1),ES(P1),ESL(P1),FQW(P1,1),FQW_LEAD(P1),FTL(P1,1),
     + FTL_LEAD(P1),TAUX(U1,1),TAUY(U1,1),QW(P1,1),RHOKE(P1),
     + RHOKEA(P1),RHOKES(P1),RHOKESL(P1),RHOKH(P1,1),RHOKM(U1,1),
     + RIB(P1),TL(P1,1),TSTAR_NL(P1),VSHR(P1),Z0H(P1),Z0M(P1),
     + Z0M_EFF(P1),H_BLEND(P1),T1_SD(P1),Q1_SD(P1),
     & RHO_CD_MODV1(P1),CDR10M(U1),CHR1P5M(P1),CER1P5M(P1),FME(P1),
     + SU10,SV10,SQ1P5,ST1P5,SFME,
     + NRML(P1),L_Z0_OROG,L_RMBL,L_BL_LSPICE,ERR,LTIMER
     +)
      IF (ERR.GT.0) THEN
        ERROR = ERR + 10
        GOTO999
      ENDIF
C
C-----------------------------------------------------------------------
CL 5.  Turbulent exchange coefficients and "explicit" fluxes between
CL     model layers in the boundary layer (P243b, routine KMKH).
C-----------------------------------------------------------------------
C
      CALL KMKH (
     + P_FIELD,U_FIELD,P1,U1,
     + P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS,N_U_ROWS,BL_LEVELS,
     + TIMESTEP,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,CCA,BQ_1,BT_1,BF_1,
     & CF,DZL,
     + PSTAR,Q,QCF,QCL,RDZ,T,TV_RDZUV,
     + U,U_P,V,V_P,Z0M_EFF,ZLB(1,0),H_BLEND,
     + FQW,FTL,TAUX,TAUY,QW,
     + RHOKM,RHOKH,TL,ZH,TV_RDZUV(1,2),RHO_KM(1,2),
     + CCB,CCT,L_MOM,L_MIXLEN,
     & L_BL_LSPICE,
     + NRML,ERR,LTIMER
     +)
      IF (ERR.GT.0) THEN
        ERROR = ERR + 20
        GOTO999
      ENDIF
C
C-----------------------------------------------------------------------
CL 6.  "Implicit" calculation of increments for TSTAR and atmospheric
CL     boundary layer variables (P244, routine IMPL_CAL).
CL     10-metre wind components are also diagnosed if requested.
C-----------------------------------------------------------------------
C
      CALL IMPL_CAL(
     + P_FIELD,U_FIELD,P1,U1,
     + P_POINTS,U_POINTS,BL_LEVELS,ROW_LENGTH,N_P_ROWS,N_U_ROWS,
     + CDR10M,U_0,V_0,SU10,SV10,
     + ASOIL_1,DELTA_AK,DELTA_BK,FQW_LEAD,FTL_LEAD,
     + RHOKE,RHOKEA,RHOKES,RHOKESL,
     + RHOKH(1,2),RHOKM(1,2),
     + ICE_FRACT,LAND_MASK,LYING_SNOW,PSTAR,RADNET,
     + SEA_ICE_HTF,SOIL_HT_FLUX,TIMESTEP,TSTAR_NL,
     + EA,ES,ESL,FQW,FTL,QW,RHOKH(1,1),RHOKM(1,1),TSTAR,TL,
     + U,V,DQW_1,DQW_RML,DTRDZ,DTRDZ_RML,DTSTAR,
     + E_SEA,H_SEA,TAUX,TAUY,
     + U10M,V10M,
     + NRML,ERR,LTIMER
     +)
      IF (ERR.GT.0) THEN
        ERROR = ERR + 30
        GOTO999
      ENDIF
C
CL 6.1 Convert FTL to sensible heat flux in Watts per square metre.
C
      DO 7 K=1,BL_LEVELS
Cfpp$ Select(CONCUR)
        DO 71 I=P1,P1+P_POINTS-1
          FTL(I,K) = FTL(I,K)*CP
   71   CONTINUE
    7 CONTINUE
C
C-----------------------------------------------------------------------
CL 7.  Surface evaporation components and updating of surface
CL     temperature (P245, routine SF_EVAP).
CL     The following diagnostics are also calculated, as requested :-
CL     Heat flux due to melting of sea-ice; specific humidity at 1.5
CL     metres; temperature at 1.5 metres.
C-----------------------------------------------------------------------
C
      CALL SF_EVAP (
     + P_FIELD,P1,LAND_FIELD,LAND1,
     + P_POINTS,BL_LEVELS,LAND_MASK,LAND_PTS,LAND_INDEX,
     + ASOIL_1,CANOPY,CATCH,DTRDZ,DTRDZ_RML,
     + EA,ESL,SOIL_HT_FLUX,E_SEA,
     + ICE_FRACT,LYING_SNOW,RADNET,SMC,TIMESTEP,
     + CER1P5M,CHR1P5M,PSTAR,Z1,Z0M_EFF,Z0H,
     + SQ1P5,ST1P5,SMLT,
     + NRML,
     + RHOKH,DQW_1,DQW_RML,DTSTAR,FTL,FQW,ES,QW,TL,TSTAR,
     + ECAN,EI
     +,SICE_MLT_HTF,Q1P5M,T1P5M,LTIMER
     +)
C
CL 7.1 Copy T and Q from workspace to INOUT space.
C
      DO K=1,BL_LEVELS
Cfpp$  Select(CONCUR)
        DO I=P1,P1+P_POINTS-1
          T(I,K)=TL(I,K)
          Q(I,K)=QW(I,K)
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
CL 8.  Calculate surface latent heat flux diagnostic.
C-----------------------------------------------------------------------
C
      IF (SLH) THEN
        DO 9 I=P1,P1+P_POINTS-1
          LATENT_HEAT(I) = LC*FQW(I,1) + LF*EI(I)
    9   CONTINUE
      ENDIF
  999  CONTINUE  ! Branch for error exit.

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
