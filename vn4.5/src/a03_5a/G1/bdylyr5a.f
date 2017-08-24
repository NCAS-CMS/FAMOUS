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
C*LL  SUBROUTINE BDY_LAYR-----------------------------------------------
CLL
CLL  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
CLL           between (a) surface and atmosphere, (b) atmospheric levels
CLL           within the boundary layer, and/or the effects of these
CLL           fluxes on the primary model variables.  The flux of heat
CLL           into and through the soil is also modelled.  Numerous
CLL           related diagnostics are also calculated.
CLL           F E Hewer, July 1990: removed call to LS_CLD.
CLL    This version passes out liquid/frozen water temperature in
CLL    array "T" (TL), and total water content in array "Q" (QW).
CLL    These may be converted to T and Q respectively by calling
CLL    the large scale cloud routine, LS_CLD.
CLL            F E Hewer, August 1990: land point data stored
CLL    on land points only (dimension: LAND_FIELD, arrays:CANOPY, CATCH
CLL    HCAP, HCON, RESIST, ROOTDEP, SMC, SMVCCL, SMVCWT, T_SOIL)
CLL    Arrays whose elements may contain values over both sea and land
CLL    points are compressed onto land points for land calculations if
CLL    defined variable IBM is NOT selected. RHOKM,RHOKH redefined as
CLL    workspace.
CLL
CLL  Suitable for single column use
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   4.1  5/6/96     New deck. C.Bunton
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL   4.3  04/02/97   Logical switches L_MOM and L_MIXLEN passed down
CLL                   to KHKH and thence EXCOEF.
CLL                                                     R.N.B.Smith
CLL   4.3  15/05/97   By-pass calls to HEAT_CON and SMC_ROOT when land
CLL                   points=0 to prevent occasional failures with
CLL                   MPP. R.Rawlins.
CLL
CLL
CLL   4.3   28/04/97  Some fields not fully initialised.
CLL   4.4   16/10/97  Minor initialisation bug. S.D.Mullerworth
CLL                   SD Mullerworth
CLL
CLL   4.4   08/09/97  L_BL_LSPICE specifies mixed phase precipitation
CLL                   scheme.                  D.Wilson 
CLL   4.5   Jul. 98  Kill the IBM specific lines. (JCThil)
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
     + P_FIELD,U_FIELD,LAND_FIELD
     +,P_ROWS,FIRST_ROW,N_ROWS,ROW_LENGTH

C IN values defining vertical grid of model atmosphere :
     +,BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK
     +,EXNER

C IN soil/vegetation/land surface data :
     +,LAND_MASK,GATHER,LAND_INDEX
     +,ST_LEVELS,SM_LEVELS,CANOPY,CATCH,HCON
     +,LYING_SNOW,RESIST,ROOTD,SMVCCL,SMVCST,SMVCWT
     +,STHF,STHU,VFRAC,Z0V,SIL_OROG_LAND,L_Z0_OROG,HO2R2_OROG
     +,HT,LAI

C IN sea/sea-ice data :
     +,DI,ICE_FRACT,U_0,V_0

C IN cloud data :
     +,CF,QCF,QCL
     +,CCA,CCB,CCT

C IN everything not covered so far :
     +,CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RADNET,TIMESTEP
     +,L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN

C INOUT data :
     +,Q,GC,T,T_SOIL,TI,TSTAR,U,V,Z0MSEA

C OUT Diagnostic not requiring STASH flags :
     &,CD,CH,E_SEA,EPOT,ETRAN,FQW,FSMC,FTL,H_SEA,RHOKH,RHOKM,RIB
     +,SEA_ICE_HTF,SURF_HT_FLUX
     +,TAUX,TAUY,VSHR

C OUT diagnostic requiring STASH flags :
     +,FME,SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT
     +,Q1P5M,T1P5M,U10M,V10M
C (IN) STASH flags :-
     +,SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10

C OUT data required for tracer mixing :
     &,RHO_ARESIST,ARESIST,RESIST_B
     &,NRML

C OUT data required for 4D-VAR :
     &,RHO_CD_MODV1,RHO_KM

C OUT data required elsewhere in UM system :
     +,ECAN,EI,ES,EXT,SNOWMELT,ZH
     +,GPP,NPP,RESP_P
     +,T1_SD,Q1_SD,ERROR
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
     +,L_RMBL                    ! IN T to use rapidly mixing boundary
C                                !    scheme in IMPL_CAL
     &,L_BL_LSPICE           ! IN 
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
     +,L_Z0_OROG                 ! IN T to use orog.roughness
C                                !    treatment in SF_EXCH
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
     + ST_LEVELS                 ! IN No. of deep soil temp. levels
     +,SM_LEVELS                 ! IN No. of soil moisture levels
      REAL
     + CANOPY(LAND_FIELD)        ! IN Surface/canopy water (kg per sq m)
     +,CATCH(LAND_FIELD)         ! IN Surface/canopy water capacity
C                                !    (kg per sq m).
C                                !    Must be global for coupled model,
C                                !    ie dimension P_FIELD not LAND_FIEL
     +,HCON(LAND_FIELD)          ! IN Soil thermal conductivity excludin
C                                !    the effects of water and ice (W/m/
     +,HT(LAND_FIELD)            ! IN Canopy height (m)
     +,LAI(LAND_FIELD)           ! IN Leaf area index.
     +,LYING_SNOW(P_FIELD)       ! IN Lying snow (kg/sq m).
     +,RESIST(LAND_FIELD)        ! IN Fixed surface resistance to
C                                !    evaporation (s/m).
     +,ROOTD(LAND_FIELD)         ! IN Depth of active soil layer ("root
C                                !    depth") (metres).
     +,SMVCCL(LAND_FIELD)        ! IN Critical volumetric SMC (cubic m
C                                !    per cubic m of soil).
     +,SMVCST(LAND_FIELD)        ! IN Volumetric saturation point (cubic
C                                !    per cubic m of soil).
     +,SMVCWT(LAND_FIELD)        ! IN Volumetric wilting point (cubic m
C                                !    per cubic m of soil).
     +,STHF(LAND_FIELD,SM_LEVELS)! IN Frozen soil moisture content of
C                                !    each layer as a fraction of
C                                !    saturation.
     +,STHU(LAND_FIELD,SM_LEVELS)! IN Unfrozen soil moisture content of
C                                !    each layer as a fraction of
C                                !    saturation.
     +,VFRAC(LAND_FIELD)         ! IN Vegetation fraction.
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
     + CO2_MMR                   ! IN CO2 Mass Mixing Ratio
     +,PHOTOSYNTH_ACT_RAD(P_FIELD) ! IN Net downward shortwave radiation
C                                !    in band 1 (w/m2).
     +,PSTAR(P_FIELD)            ! IN Surface pressure (Pascals).
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
     +,SIMLT   ! IN Flag for SICE_MLT_HTF (q.v.)
     +,SMLT    ! IN Flag for SNOMLT_SURF_HTF (q.v.)
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
     +,GC(LAND_FIELD)                  ! INOUT "Stomatal" conductance to
C                                      !       evaporation (m/s).
     +,T(P_FIELD,BL_LEVELS)            ! INOUT Input:atmospheric temp(K)
C                                      !      Output:liquid/frozen water
C                                      !       temperature (TL) (K)
     +,T_SOIL(LAND_FIELD,ST_LEVELS)    ! INOUT Soil temperatures (K).
     +,TI(P_FIELD)                     ! INOUT Sea-ice surface layer
C                                      !       temperature (K).
     +,TSTAR(P_FIELD)                  ! INOUT Surface temperature (K).
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
     &,EPOT(P_FIELD)             ! OUT potential evaporation (kg/m2/s).
     +,FQW(P_FIELD,BL_LEVELS)    ! OUT Moisture flux between layers
C                                !     (kg per square metre per sec).
C                                !     FQW(,1) is total water flux
C                                !     from surface, 'E'.
     &,FSMC(LAND_FIELD)          ! OUT soil moisture availability.
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
     +,SURF_HT_FLUX(P_FIELD)     ! OUT Net downward heat flux at surface
C                                !     over land or sea-ice fraction of
C                                !     gridbox (W/m2).
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
C                                  !     coefficient for momentum before
C                                  !     interpolation.
     &,RHO_ARESIST(P_FIELD)      ! OUT RHOSTAR*CD_STD*VSHR for SULPHUR c
     &,ARESIST(P_FIELD)          ! OUT 1/(CD_STD*VSHR) for Sulphur cycle
     &,RESIST_B(P_FIELD)         ! OUT (1/CH-1/(CD_STD)/VSHR for Sulpur 

      INTEGER
     & NRML(P_FIELD)             ! OUT Number of model layers in the
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
     +,SNOMLT_SURF_HTF(P_FIELD) ! OUT Heat flux required for surface
C                               !     melting of snow (W/m2).
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
     + EI(P_FIELD)    ! OUT Sublimation from lying snow or sea-ice
C                     !     (kg per sq m per sec).
     +,ECAN(P_FIELD)  ! OUT Gridbox mean evaporation from canopy/surface
C                     !     store (kg per sq m per s).  Zero over sea.
     +,ES(P_FIELD)    ! OUT Surface evapotranspiration through a
C                     !     resistance which is not entirely aerodynamic
C                     !     i.e. "soil evaporation".  Always non-
C                     !     negative.  Kg per sq m per sec.
     +,ETRAN(P_FIELD) ! OUT Transpiration (kg/m2/s).
     +,EXT(LAND_FIELD,SM_LEVELS)
C                     ! OUT Extraction of water from each soil layer
C                     !     (kg/m2/s).
     +,GPP(LAND_FIELD)! OUT Gross primary productivity (kg C/m2/s).
     +,NPP(LAND_FIELD)! OUT Net primary productivity (kg C/m2/s).
     +,RESP_P(LAND_FIELD)
C                     ! OUT Plant respiration (kg C/m2/s).
     +,SNOWMELT(P_FIELD) ! OUT Snowmelt (kg/m2/s).
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
C*----------------------------------------------------------------------
C*L---------------------------------------------------------------------
C  External routines called :-
C
      EXTERNAL Z,SICE_HTF,SF_EXCH,KMKH,IMPL_CAL,SF_EVAP
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
      REAL
     + DZSOIL(4)               ! Soil layer thicknesses (m).
      DATA DZSOIL /0.100, 0.250, 0.650, 2.000 /
C-----------------------------------------------------------------------
C*----------------------------------------------------------------------
C
C  Workspace :-
C
      REAL
     + ALPHA1(P_FIELD)          ! WORK Gradient of saturated
C                               !      specific humidity with
C                               !      respect to temperature between
C                               !      the bottom model layer and the
C                               !      surface
     +,ASHTF(P_FIELD)           ! WORK Coefficient to calculate surface
C                               !      heat flux into soil or sea-ice.
     +,ASURF(P_FIELD)           ! WORK Reciprocal areal heat capacity
C                               !      of soil layer or sea-ice
C                               !      surface layer (K m**2 / J).
     +,BQ_1(P_FIELD)            ! WORK A buoyancy parameter for the
C                               !      lowest atmospheric level.
     +,BT_1(P_FIELD)            ! WORK A buoyancy parameter for the
C                               !      lowest atmospheric level.
     &,BF_1(P_FIELD)            
!        WORK   A bouyancy parameter for the lowest atmospheric level
     +,DQW_1(P_FIELD)           ! WORK Increment for QW(,1).
     +,DTRDZ(P_FIELD,BL_LEVELS) ! WORK -g.dt/dp for model layers.
     +,DTRDZ_RML(P_FIELD)       ! WORK -g.dt/dp for the rapidly
C                               !      mixing layer.
     +,DZL(P_FIELD,BL_LEVELS)   ! WORK DZL(,K) is depth in m of layer
C                               !      K, i.e. distance from boundary
C                               !      K-1/2 to boundary K+1/2.
     +,ESOIL(P_FIELD)           ! WORK Evaporation from bare soil (kg/m2
     +,FRACA(P_FIELD)           ! WORK Fraction of surface moisture flux
C                               !      with only aerodynamic resistance.
     +,F_SE(P_FIELD)            ! WORK Fraction of the evapotranspiratio
C                               !      which is bare soil evaporation.
     +,HCONS(LAND_FIELD)        ! WORK Soil thermal conductivity includi
C                               !      the effects of water and ice (W/m
     +,QW(P_FIELD,BL_LEVELS)    ! WORK Total water content, but
C                               !      replaced by specific humidity
C                               !      in LS_CLD.
     +,RESFS(P_FIELD)           ! WORK Combined soil, stomatal
C                               !      and aerodynamic resistance
C                               !      factor = PSIS/(1+RS/RA) for
C                               !      fraction (1-FRACA)
     +,RESFT(P_FIELD)           ! WORK Total resistance factor
C                               !      FRACA+(1-FRACA)*RESFS.
     +,RHOKE(P_FIELD)           ! WORK Surface exchange coefficient for
C                               !      FQW.
     +,RHOKPM(P_FIELD)          ! WORK Surface exchange coeff.
     &,RHOKPM_POT(P_FIELD)      ! WORK Surface exchange coeff. for 
!                                      potential evaporation
     +,RDZ(P_FIELD,BL_LEVELS)   ! WORK RDZ(,1) is the reciprocal
C                               !      of the height of level 1, i.e.
C                               !      of the middle of layer 1.  For
C                               !      K > 1, RDZ(,K) is the
C                               !      reciprocal of the vertical
C                               !      distance from level K-1 to
C                               !      level K.
     +,SMC(LAND_FIELD)          ! WORK Soil moisture content (kg/m2).
     +,TL(P_FIELD,BL_LEVELS)    ! WORK Ice/liquid water temperature,
C                               !      but replaced by T in LS_CLD.
     +,TV_RDZUV(P_FIELD,BL_LEVELS) ! WORK Virtual temp at start (TV).
C                                  !      RDZ (K > 1) on UV-grid.
C                                  !      Comments as per RHOKM (RDZUV).
     +,U_P(P_FIELD,BL_LEVELS)   ! WORK U on P-grid.
     +,V_P(P_FIELD,BL_LEVELS)   ! WORK V on P-grid.
     +,V_ROOT(LAND_FIELD)       ! WORK Volumetric soil moisture
C                               !      concentration in the rootzone
C                               !      (m3 H2O/m3 soil).
     +,V_SOIL(LAND_FIELD)       ! WORK Volumetric soil moisture
C                               !      concentration in the top
C                               !      soil layer (m3 H2O/m3 soil).
     +,WT_EXT(LAND_FIELD,SM_LEVELS)! WORK Fraction of transpiration whic
C                               !      extracted from each soil layer.
     +,ZLB(P_FIELD,0:BL_LEVELS) ! WORK ZLB(,K) is the height of the
C                               !      upper boundary of layer K
C                               !      ( = 0.0 for "K=0").
       REAL
     + Z0H(P_FIELD)               ! WORK Roughness length for heat and
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
C   Variables for Vegetation Thermal Canopy
C
      REAL
     + CANCAP(P_FIELD)            ! WORK Volumetric heat capacity of
C                                 !      vegetation canopy (J/Kg/m3).
     +,RADNET_C(P_FIELD)          ! WORK Adjusted net radiation for 
C                                 !      vegetation canopy over land
C                                 !      (W/m2).

      INTEGER
     + F_TYPE(LAND_FIELD)         ! WORK Plant functional type:
C                                 !       1 - Broadleaf Tree
C                                 !       2 - Needleleaf Tree
C                                 !       3 - C3 Grass
C                                 !       4 - C4 Grass
C
C  Local scalars :-
C
      INTEGER
     + ERR        ! LOCAL Return codes from lower-level routines.
     +,I,J,L      ! LOCAL Loop counter (horizontal field index).
     +,K,N        ! LOCAL Loop counter (vertical level index).
     +,N_P_ROWS   ! LOCAL No of P-rows being processed.
     +,N_U_ROWS   ! LOCAL No of UV-rows being processed.
     +,P_POINTS   ! LOCAL No of P-points being processed.
     +,P1         ! LOCAL First P-point to be processed.
     +,LAND1      ! LOCAL First land-point to be processed.
C                 !           1 <= LAND1 <= LAND_FIELD
     +,LAND_PTS   ! LOCAL No of land points being processed.
     +,U_POINTS   ! LOCAL No of UV-points being processed.
     +,U1         ! LOCAL First UV-point to be processed.

C-----------------------------------------------------------------------
C Functional Type dependent parameters
C-----------------------------------------------------------------------
      INTEGER
     + R_LAYERS(4)   ! Number of soil layers from which
                     ! water can be extracted
C-----------------------------------------------------------------------
C                       BT    NT   C3G   C4G
C-----------------------------------------------------------------------
      DATA R_LAYERS/     4,    4,    3,    3 /

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF
      ERROR = 0
C-----------------------------------------------------------------------
C Initialise RADNET_C to be the same as RADNET over all points
C-----------------------------------------------------------------------
      DO I=1,P_FIELD
        RADNET_C(I) = RADNET(I)
      ENDDO

C
C-----------------------------------------------------------------------
CL 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
CL    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
CL    most UV-rows are used only for interpolation and are not updated.
C-----------------------------------------------------------------------
C
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

!-----------------------------------------------------------------------
! Diagnose the plant functional types at each location.
! Assume : Broadleaf Trees if rootdepth > 0.8m
!          C3 Grass        if rootdepth < 0.8m
!-----------------------------------------------------------------------
      DO L=1,LAND_FIELD
        IF (ROOTD(L).GT.0.8) THEN
          F_TYPE(L)=1
        ELSE
          F_TYPE(L)=3
        ENDIF
      ENDDO

      IF(LAND_FIELD.GT.0) THEN    ! Omit if no land points
!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.
!-----------------------------------------------------------------------
      CALL HEAT_CON (LAND_FIELD,HCON,STHU,STHF,SMVCST,HCONS,LTIMER)

!-----------------------------------------------------------------------
! Calculate the soil moisture in the root zone.
!-----------------------------------------------------------------------
      CALL SMC_ROOT (LAND_FIELD,SM_LEVELS,F_TYPE,DZSOIL,ROOTD,STHU,   
     &         VFRAC,SMVCST,SMVCWT,SMC,V_ROOT,V_SOIL,WT_EXT,LTIMER)     
      ENDIF                     ! End test on land points

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
C
C-----------------------------------------------------------------------
CL 3.  Calls to SICE_HTF now after IMPL_CAL
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
CL 4.  Surface turbulent exchange coefficients and "explicit" fluxes
CL     (P243a, routine SF_EXCH).
CL     Wind mixing "power" and some values required for other, later,
CL     diagnostic calculations, are also evaluated if requested.
C-----------------------------------------------------------------------
C

      IF(LAND_FIELD.GT.0) THEN
C Initialise any uncalculated points
      DO I=1,LAND1
        GPP(I)=0.
        NPP(I)=0.
        RESP_P(I)=0.
        FSMC(I)=0.
      ENDDO
      DO I=LAND_PTS+LAND1-1,LAND_FIELD
        GPP(I)=0.
        NPP(I)=0.
        RESP_P(I)=0.
        FSMC(I)=0.
      ENDDO
      ENDIF                        ! if land points exist

      CALL SF_EXCH (
     + P_POINTS,LAND_PTS,U_POINTS,ROW_LENGTH,N_P_ROWS,N_U_ROWS,
     + LAND_INDEX(LAND1),P1,GATHER,
     + AK(1),BK(1),
     + CANOPY(LAND1),CATCH(LAND1),CO2_MMR,CF(P1,1),
     + SM_LEVELS,DZSOIL(1),HCONS(LAND1),F_TYPE(LAND1),
     + HT(LAND1),LAI(LAND1),PHOTOSYNTH_ACT_RAD(P1),GPP(LAND1),
     + NPP(LAND1),RESP_P(LAND1),
     + ICE_FRACT(P1),LAND_MASK(P1),LYING_SNOW(P1),PSTAR(P1),Q(P1,1),
     + QCF(P1,1),QCL(P1,1),RADNET_C(P1),GC(LAND1),RESIST(LAND1),
     + ROOTD(LAND1),SMC(LAND1),SMVCCL(LAND1),SMVCWT(LAND1),
     + T(P1,1),TIMESTEP,TI(P1),T_SOIL(LAND1,1),TSTAR(P1),
     + U(U1,1),V(U1,1),U_P(P1,1),V_P(P1,1),U_0(U1),V_0(U1),
     + V_ROOT(LAND1),V_SOIL(LAND1),VFRAC(LAND1),
     + Z0V(P1),SIL_OROG_LAND(LAND1),Z1(P1),
     + CANCAP(P1),Z0MSEA(P1),HO2R2_OROG(LAND1), 
     & ALPHA1(P1),ASHTF(P1),BQ_1(P1),BT_1(P1),BF_1(P1),CD(P1),CH(P1),
     & EPOT(P1),FQW(P1,1),FSMC(LAND1),FTL(P1,1),E_SEA(P1),H_SEA(P1),
     + TAUX(U1,1),TAUY(U1,1),QW(P1,1),FRACA(P1),RESFS(P1),F_SE(P1),
     & RESFT(P1),RHOKE(P1),RHOKH(P1,1),RHOKM(U1,1),
     & RHOKPM(P1),RHOKPM_POT(P1),
     + RIB(P1),TL(P1,1),VSHR(P1),Z0H(P1),Z0M(P1),
     + Z0M_EFF(P1),H_BLEND(P1),T1_SD(P1),Q1_SD(P1),
     + RHO_CD_MODV1(P1),CDR10M(U1),CHR1P5M(P1),CER1P5M(P1),FME(P1),
     + SU10,SV10,SQ1P5,ST1P5,SFME,
     + RHO_ARESIST(P1),ARESIST(P1),RESIST_B(P1),
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
     + ALPHA1,ASHTF,CDR10M,DELTA_AK,DELTA_BK,
     + RHOKH(1,2),RHOKM(1,2),
     + ICE_FRACT,LYING_SNOW,PSTAR,RADNET_C,RESFT,RHOKPM,           
     & RHOKPM_POT,
     + U_0,V_0,TIMESTEP,LAND_MASK,SU10,SV10,
     & EPOT,FQW,FTL,E_SEA,H_SEA,QW,
     + RHOKE,RHOKH(1,1),RHOKM(1,1),TL,U,V,
     + DTRDZ,DTRDZ_RML,TAUX,TAUY,SURF_HT_FLUX,U10M,V10M,NRML,
     + ERR,LTIMER
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
C    Diagnose surface temperature and increment sub-surface temperatures
C    for land and sea-ice.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL   Sea-ice (P241, routine SICE_HTF).
C-----------------------------------------------------------------------
C
      CALL SICE_HTF(
     + ASHTF(P1),DI(P1),ICE_FRACT(P1),SURF_HT_FLUX(P1),TIMESTEP,
     + LAND_MASK(P1),P_POINTS,TI(P1),TSTAR(P1),ASURF(P1),
     + SEA_ICE_HTF(P1),LTIMER
     +)
C
C-----------------------------------------------------------------------
CL   Diagnose the land surface temperature (previously in SOIL_HTF)
C-----------------------------------------------------------------------
C
      DO I=LAND1,LAND1+LAND_PTS-1
        J = LAND_INDEX(I)
        TSTAR(J) = T_SOIL(I,1) + SURF_HT_FLUX(J) / ASHTF(J)
      ENDDO
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
     + ALPHA1,ASURF,ASHTF,CANOPY,CATCH,
     + DTRDZ,DTRDZ_RML,E_SEA,FRACA,
     + ICE_FRACT,NRML,RHOKH,SMC,TIMESTEP,CER1P5M,CHR1P5M,
     + PSTAR,RESFS,RESFT,Z1,Z0M,Z0H,SQ1P5,ST1P5,SIMLT,SMLT,
     + FTL,FQW,LYING_SNOW,QW,SURF_HT_FLUX,TL,TSTAR,TI,
     + ECAN,ES,EI,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,
     + Q1P5M,T1P5M,LTIMER
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

!-----------------------------------------------------------------------
! Diagnose the soil evaporation, the transpiration and the water
! extracted from each soil layer
!-----------------------------------------------------------------------
      DO N=1,SM_LEVELS
        DO I=LAND1,LAND1+LAND_PTS-1
          J = LAND_INDEX(I)
          EXT(I,N)=WT_EXT(I,N)*(1-F_SE(J))*ES(J)
        ENDDO
      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
C Initialise ETRAN otherwise sea points remain uninitialised
      DO I=1,P_FIELD
        ETRAN(I)=0.
      ENDDO
      DO I=LAND1,LAND1+LAND_PTS-1
        J = LAND_INDEX(I)
          ESOIL(J)=F_SE(J)*ES(J)
          ETRAN(J)=(1-F_SE(J))*ES(J)
          EXT(I,1)=EXT(I,1)+ESOIL(J)
      ENDDO

C-----------------------------------------------------------------------
C Diagnose the true value of the surface soil heat flux over land points
C-----------------------------------------------------------------------
      DO I=LAND1,LAND1+LAND_PTS-1
        J = LAND_INDEX(I)
        SURF_HT_FLUX(J) = SURF_HT_FLUX(J) - CANCAP(J) *
     +          ( TSTAR(J) - T_SOIL(I,1) ) / TIMESTEP
      ENDDO


      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
