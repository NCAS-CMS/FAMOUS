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
C*LL  SUBROUTINE IMPL_CAL ----------------------------------------------
CLL
CLL  Purpose: Calculate increments for surface temperature and
CLL           atmospheric variables in the boundary layer, using an
CLL           implicit numerical scheme.  The tridiagonal matrices are
CLL           inverted using simple Gaussian elimination.  See written
CLL           documentation for a very detailed explanation.
CLL
CLL  Suitable for single column use; activate *IF definition IBM.
CLL
CLL  Model           Modification history from model version 3.0
CLL version  Date
CLL  3.1    12/01/93 Alternative, more complete implicit numerical
CLL                  scheme made available under *IF DEF,A03_2C.
CLL  3.3  07/07/93  Only versions 1B & 2C carried forward to UM3.3
CLL  3.4  06/06/94  DEF TIMER replaced by LOGICAL LTIMER
CLL                  Argument LTIMER added
CLL                                                 S.J.Swarbrick
!     3.5    9/5/95   MPP code: Change updateable area  P.Burton
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
!LL   4.3  14/01/97   MPP code : Corrected setting of polar rows
!LL                                                     P.Burton
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL
CLL  FH, RNBS  <- Programmers of some or all of previous code or changes
CLL
CLL  Programming standard: UM Documentation Paper No 4, Version 2,
CLL                        dated 18/1/90
CLL
CLL  System component covered: P244
CLL
CLL  Project task: P24
CLL
CLL  Documentation: UM Documentation Paper No 24.
CLL
C*----------------------------------------------------------------------
C*L  Arguments :-
      SUBROUTINE IMPL_CAL (
     + P_FIELD,U_FIELD,P1,U1,
     + P_POINTS,U_POINTS,BL_LEVELS,ROW_LENGTH,P_ROWS,U_ROWS
     +,CDR10M,U0,V0,SU10,SV10
     +,ASOIL_1,DELTA_AK,DELTA_BK,FQW_LEAD,FTL_LEAD
     +,GAMMA_RHOKE,GAMMA_RHOKEA,GAMMA_RHOKES,GAMMA_RHOKESL
     +,GAMMA_RHOKH_RDZ,GAMMA_RHOKM_RDZUV
     +,ICE_FRACT,LAND_MASK,LYING_SNOW,PSTAR,RADNET
     +,SEA_ICE_HTF,SOIL_HT_FLUX,TIMESTEP,TSTAR_NL
     +,EA,ES,ESL,FQW,FTL,QW,GAMMA_RHOKH_1,GAMMA_RHOKM_1,TSTAR,TL,U,V
     +,DQW_1,DQW_RML,DTRDZ,DTRDZ_RML,DTSTAR,E_SEA,H_SEA,TAUX,TAUY
     +,U10M,V10M
     +,NRML
     +,ERROR,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + P_FIELD                     ! IN No. of points in P-grid.
     +,U_FIELD                     ! IN No. of points in UV-grid.
     +,P1                          ! IN First point to be processed in
C                                  !    P-grid.
     +,U1                          ! IN First point to be processed in
C                                  !    UV-grid.
     +,P_POINTS                    ! IN Number of P-grid points to be
C                                  !    processed.
     +,U_POINTS                    ! IN Number of UV-grid points.
     +,BL_LEVELS                   ! IN No. of atmospheric levels for
C                                  !    which boundary layer fluxes are
C                                  !    calculated.
     +,ROW_LENGTH                  ! IN No. of points in latitude row.
     +,P_ROWS                      ! IN No. of P-rows of data to be
C                                  !    processed.
     +,U_ROWS                      ! IN No. of UV-rows of data to be
C                                  !    processed.
      REAL
     + CDR10M(U_FIELD)             ! IN Used in calc of 10m wind - from
C                                  !    P243 (routine SF_EXCH).  First
C                                  !    and last rows are "missing data"
C                                  !    and not used.  UV-grid.
     +,U0(U_FIELD)                 ! IN Westerly component of surface
C                                  !    current (m/s; 0 over land). UVG.
     +,V0(U_FIELD)                 ! IN Southerly component of surface
C                                  !    current (m/s; 0 over land). UVG.
      REAL
     + ASOIL_1(P_FIELD)            ! IN Soil thermodynamic coefficient
C                                  !    from P242 (SOIL_HTF).  Sq m K
C                                  !    per Joule * timestep.
     +,DELTA_AK(BL_LEVELS)         ! IN Difference of hybrid 'A' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     +,DELTA_BK(BL_LEVELS)         ! IN Difference of hybrid 'B' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     +,FQW_LEAD(P_FIELD)           ! IN "Explicit" surface flux of
C                                  !     QW (i.e. evaporation), on
C                                  !     P-grid for leads fraction
C                                  !     of gridsquare (kg/m2/s).
C                                  !     Missing data at non-sea-ice pts
     +,FTL_LEAD(P_FIELD)           ! IN "Explicit surface flux of TL
C                                  !     = H/CP (sensible heat/Cp) for
C                                  !     leads fraction of gridsquare.
C                                  !     Missing data at non-sea-ice pts
     +,GAMMA_RHOKE(P_FIELD)        ! IN Surface exchange coeff. for FQW,
C                                  !    =GAMMA(1)*RHOKE from SF_EXCH.
     +,GAMMA_RHOKEA(P_FIELD)       ! IN Surface exchange coeff. for EA,
C                                  !    =GAMMA(1)*RHOKEA from SF_EXCH.
     +,GAMMA_RHOKES(P_FIELD)       ! IN Surface exchange coeff. for ES,
C                                  !    =GAMMA(1)*RHOKES from SF_EXCH.
     +,GAMMA_RHOKESL(P_FIELD)      ! IN Surface exchange coeff. for ESL,
C                                  !    =GAMMA(1)*RHOKESL from SF_EXCH.
     +,GAMMA_RHOKH_RDZ(P_FIELD,2:BL_LEVELS)
C                                  ! IN Exchange coeff for FTL above
C                                  !    surface, =GAMMA(K)*RHOKH(,K)
C                                  !    *RDZ(K) for K>=2 (from KMKH).
     +,GAMMA_RHOKM_RDZUV(U_FIELD,2:BL_LEVELS)
C                                  ! IN Exchange coefficients for
C                                  !    momentum, on UV-grid with
C                                  !    first and last rows ignored.
C                                  !     =GAMMA(K)*RHOKM(,K)*RDZUV(,K)
C                                  !    for K>=2 (from KMKH).
      REAL                         ! Split to avoid > 19 continuations.
     + ICE_FRACT(P_FIELD)          ! IN Fraction of grid-box which is
C                                  !    sea-ice (decimal fraction).
     +,LYING_SNOW(P_FIELD)         ! IN Lying snow (kg per sq m ie "mm")
     +,PSTAR(P_FIELD)              ! IN Surface pressure (Pa).
     +,RADNET(P_FIELD)             ! IN Area weighted ice component
C                                  !    of surface net radiation.
C                                  !    (+ve downwards, W per sq m)
     +,SEA_ICE_HTF(P_FIELD)        ! IN Heat flux through sea-ice (W per
C                                  !    sq m, +ve downwards).  From P241
C                                  !    (routine SICE_HTF).
     +,SOIL_HT_FLUX(P_FIELD)       ! IN Heat flux from soil layer 1 to
C                                  !    soil layer 2 (deep soil layer 1)
C                                  !    i.e. +ve is downwards (W per sq
C                                  !    m). From P242 (SOIL_HTF).
     +,TIMESTEP                    ! IN Timestep in seconds.
     +,TSTAR_NL(P_FIELD)           ! IN TSTAR No Leads: surface temp-
C                                  !    erature of the ice at sea-ice
C                                  !    points and =TSTAR elsewhere.(K)
      LOGICAL
     + LAND_MASK(P_FIELD)          ! IN T for land, F elsewhere.
     +,SU10                        ! IN STASH flag for 10-metre W wind.
     +,SV10                        ! IN STASH flag for 10-metre S wind.
C
C  Next 5 arrays are all IN as "explicit" fluxes from P243 (SF_EXCH and
C  possibly KMKH), and OUT as "implicit" fluxes.
C
      REAL
     + EA(P_FIELD)                 ! INOUT Surface evaporation with only
C                                  !       aerodynamic resistance (+ve),
C                                  !       or condensation (-ve),
C                                  !       averaged over gridbox.
C                                  !       Kg per sq m per sec.
     +,ES(P_FIELD)                 ! INOUT Surface evapotranspiration
C                                  !       (through a not-entirely
C                                  !       aerodynamic resistance).
C                                  !       Always non-negative.  Kg per
C                                  !       sq m per sec.
     +,ESL(P_FIELD)                ! INOUT ES without fractional
C                                  !       weighting factor FRACS.  "L"
C                                  !       is for "local".  Kg/sq m/sec.
     +,FQW(P_FIELD,BL_LEVELS)      ! INOUT Flux of QW (ie., for surface,
C                                  !       total evaporation). Kg/sq m/s
     +,FTL(P_FIELD,BL_LEVELS)      ! INOUT Flux of TL (ie., for surface,
C                                  !       H/Cp where H is sensible heat
C                                  !       in W per sq m).
     +,QW(P_FIELD,BL_LEVELS)       ! INOUT Total water content (kg per
C                                  !       kg air).  From P243.
     +,GAMMA_RHOKH_1(P_FIELD)      ! IN    Surface exchange coeffs for
C                                  !       FTL, =GAMMA(1)*RHOKH(,1)
C                                  !       from P243 (SF_EXCH).
C                                  !   OUT =RHOKH_1 to satisfy STASH.
     +,GAMMA_RHOKM_1(U_FIELD)      ! IN    Surface exchange coeffs for
C                                  !       momentum, on UV-grid
C                                  !       with first and last rows
C                                  !       ignored. =GAMMA(1)*RHOKM(,1)
C                                  !       from P243 (SF_EXCH).
C                                  !   OUT =RHOKM_1 to satisfy STASH.
     +,TSTAR(P_FIELD)              ! INOUT Gridbox mean surface
C                                  !       temperature.
     +,TL(P_FIELD,BL_LEVELS)       ! INOUT Liquid/frozen water
C                                  !       temperature (K).  From P243.
     +,U(U_FIELD,BL_LEVELS)        ! INOUT Westerly wind (m/s).  On UV-
C                                  !       grid, 1st & last rows unused.
     +,V(U_FIELD,BL_LEVELS)        ! INOUT Southerly wind (m/s).  On UV-
C                                  !       grid, 1st & last rows unused.
      REAL
     + DQW_1(P_FIELD)              ! OUT Increment for QW(,1) (needed
C                                  !     again in P245).
     +,DQW_RML(P_FIELD)            ! OUT Rapidly mixing layer increment
C                                  !     to QW (needed again in P245).
     +,DTRDZ(P_FIELD,BL_LEVELS)    ! OUT -g.dt/dp for bottom BL_LEVELS
C                                  !     model layers (needed in P245).
     +,DTRDZ_RML(P_FIELD)          ! OUT -g.dt/dp for the rapidly
C                                  !     mixing layer (needed in P245).
     +,DTSTAR(P_FIELD)             ! OUT Surface temperature increment
C                                  !     (needed again in P245).
     +,E_SEA(P_FIELD)              ! OUT Evaporation from sea times
C                                  !     leads fraction. Zero over land
C                                  !     (kg per square metre per sec).
     +,H_SEA(P_FIELD)              ! OUT Surface sensible heat flux
C                                  !     over sea times leads fraction.
C                                  !     Zero over land. ( W/m2)
     +,TAUX(U_FIELD,BL_LEVELS)     ! INOUT x-component of turbulent
C                                  !       stress at levels k-1/2;
C                                  !     eg. TAUX(,1) is surface stress.
C                                  !     UV-grid, 1st and last rows set
C                                  !     to "missing data". (N/sq m)
     +,TAUY(U_FIELD,BL_LEVELS)     ! INOUT y-component of turbulent
C                                  !       stress at levels k-1/2;
C                                  !     eg. TAUY(,1) is surface stress.
C                                  !     UV-grid, 1st and last rows set
C                                  !     to "missing data". (N/sq m)
     +,U10M(U_FIELD)               ! OUT Westerly wind at 10m (m/s).
C                                  !     1st & last rows "missing data".
     +,V10M(U_FIELD)               ! OUT Southerly wind at 10m (m/s).
C                                  !     1st & last rows "missing data".
      INTEGER
     + NRML(P_FIELD)               ! IN The number of model layers
C                                  !    in the unstable rapidly mixing
C                                  !    layer. Zero if surface layer
C                                  !    is stable.
     +,ERROR                       ! OUT 1 if bad arguments, else 0.
C*
C*L  External references :-
      EXTERNAL QSAT
      EXTERNAL P_TO_UV
      EXTERNAL TIMER
C*
C*L  Local and other symbolic constants :-
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------
      REAL LS
      PARAMETER (
     + LS=LC+LF     ! Latent heat of sublimation (J per kg).
     +)
C*
C*L Workspace :-
C   6*BL_LEVELS + 2 blocks of real workspace are required.
      REAL
     + ALPHAS(P_FIELD)             ! Partial dQsat/dT* from Clausius-
C                                  ! Clapeyron relation.
     +,AQ_AM(U_FIELD,BL_LEVELS)    ! As AQ: "Q" elements of matrix in
C                                  ! eqn P244.79 (modified during
C                                  ! Gaussian elimination process).
C                                  ! As AM: elements of matrix in eqn
C                                  ! P244.80 (also get modified).
     +,AQ_RML(U_FIELD)             ! Matrix element for humidity in
C                                  ! rapidly mixing layer. Then briefly
C                                  ! used for DELTAP on the UV grid.
     +,AT_ATQ(P_FIELD,BL_LEVELS)   ! Elements in atmospheric T rows of
C                                  ! matrix in eqn P244.79 (modified
C                                  ! during Gaussian elimination).
     +,AT_RML(P_FIELD)             ! Matrix element for temperature in
C                                  ! rapidly mixing layer.
     +,ATSTAR(P_FIELD)             ! Element in surface temperature row
C                                  ! of matrix in eqn P244.79 (modified
C                                  ! during Gaussian elimination).
C                                  ! Also used for saturated sp hum, and
C                                  ! in interpolations to UV-grid.
     +,DELTAP(U_FIELD,BL_LEVELS)   ! Vertical pressure difference across
C                                  ! hybrid layers (upper minus lower)
C                                  ! (Pa).
     +,DELTAP_RML(U_FIELD)         ! Vertical pressure difference across
C                                  ! the rapidly mixing layer (Pa).
     +,DQW_DU(U_FIELD,BL_LEVELS)   ! As DQW: delta QW elements of vector
C                                  ! on RHS, then LHS, of eqn P244.79.
C                                  ! As DU: delta U elements of vector
C                                  ! on RHS, then LHS, of eqn P244.80.
     +,DTL_DV(U_FIELD,BL_LEVELS)   ! As DTL: delta TL (for atmosphere)
C                                  ! elements of vector on RHS, then
C                                  ! LHS, of eqn P244.79.
C                                  ! As DV: delta V elements of vector
C                                  ! on RHS, then LHS, of eqn P244.80.
     +,DTL_RML(U_FIELD)            ! Delta TL for rapidly mixing layer.
     +,DTRDZ_UV(U_FIELD,BL_LEVELS) ! -g.dt/dp for model layers
C                                  ! interpolated to the UV-grid.
     +,FQW_ICE(P_FIELD)            ! "Explicit" surface flux of QW for
C                                  ! sea-ice fraction of gridsquare.
     +,FTL_ICE(P_FIELD)            ! "Explicit" surface flux of TL for
C                                  ! sea-ice fraction of gridsquare.
C*
C  Local scalars :-
      REAL
     + CTQ      ! Matrix element in P244.??, for local increments to rml
     +,CM       ! Matrix element in eqn P244.80.
     +,CQ       ! Matrix element in "Q" row in eqn P244.79.
     +,CQ_RML   ! As above but for rapidly mixing layer increment.
     +,CT       ! Matrix element in "T" row in eqn P244.79.
     +,CT_RML   ! As above but for rapidly mixing layer increment.
     +,RBTQ     ! Reciprocal of B P244.??, for local increments to rml
     +,RBM      ! Reciprocal of BM(') (eqns P244.81, 85, 89).
     +,RBQ      ! Reciprocal of BQ(') (eqns P244.98, 101, 104).
     +,RBQ_RML  ! As above but for the rapidly mixing layer increment.
     +,RBT      ! Reciprocal of BT' (eqns P244.107, 110, 113).
     +,RBT_RML  ! As above but for the rapidly mixing layer increment.
     +,DELTDQ   ! P244.126 times GAMMA(1) (factor required in eqns
C               ! P244.122-4).
     +,DTRICE   ! = DTSTAR/ICE_FRACT at sea-ice points.
     +,DTIMEG   ! TIMESTEP * G (used in several places).
     +,LAT_HEAT ! Latent heat (either LC or LS) in section 2.2.
      INTEGER
     + BLM1     ! BL_LEVELS minus 1.
     +,NRMLP1   ! NRML plus 1.
     +,I        ! Loop counter (horizontal field index).
     +,J        ! Offset version of I.
     +,K        ! Loop counter (vertical index).
     +,KM1      ! K minus 1.
     +,KP1      ! K plus 1.
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
C
C-----------------------------------------------------------------------
CL  0.  Check that the scalars input to define the grid are consistent.
C       See comments to routine SF_EXCH for details.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
      CALL TIMER('IMPLCAL ',3)
      ENDIF
      ERROR=0
      IF (                                                          
     +     U_POINTS .NE. (U_ROWS*ROW_LENGTH) .OR.
     +     P_POINTS .NE. (P_ROWS*ROW_LENGTH) )  THEN
        ERROR=1
        GOTO999
      ENDIF
      DTIMEG = TIMESTEP * G
      BLM1 = BL_LEVELS-1
C
CL----------------------------------------------------------------------
CL (A) Calculations on P-grid.
CL----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL 1.  Calculate pressure across layer (in hybrid coordinates), DELTAP,
CL     and then -gdt/dP = dt/rho*dz for use throughout section (A)
C-----------------------------------------------------------------------
C
      DO 1 K=1,BL_LEVELS
        DO 11 I=P1,P1+P_POINTS-1
          DELTAP(I,K) = DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)
          DTRDZ(I,K) = -DTIMEG / DELTAP(I,K)
   11   CONTINUE
   1  CONTINUE
C
C-----------------------------------------------------------------------
CL 2.  Calculate implicit T and Q increments due to local mixing within
CL     the rapidly mixing layer (where it exists).
CL     The surface fluxes FTL(I,1), FQW(I,1) are used for calculating
CL     the rapidly mixing layer (rml) increments but not here.
CL     Therefore the matrix equation we must solve to find the implicit
CL     T and Q increments due to local mixing within the rml does not
CL     have a "surface" row and we can solve for the T and Q increments
CL     for K = 1 to NRML simultaneously.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL 2.1 Start 'upward sweep' with lowest model layer, which will be the
CL     bottom of the rapidly mixing layer (rml) if it exists.
C-----------------------------------------------------------------------
C
      DO 21 I=P1,P1+P_POINTS-1
        IF (NRML(I) .GE. 2) THEN
C
C  "Explicit" increments due to local mixing within the rml.
C  P244.49/31 but surface flux used in rml increment calculations.
C
          DQW_DU(I,1) = -DTRDZ(I,1) * FQW(I,2)
          DTL_DV(I,1) = -DTRDZ(I,1) * FTL(I,2)
C
C  Define matrix elements (CTQ always zero for this case).
C
          AT_ATQ(I,1) = -DTRDZ(I,1) * GAMMA_RHOKH_RDZ(I,2)    ! P244.28
          RBTQ = 1.0 / ( 1.0 - AT_ATQ(I,1) )    ! Reciprocal of P244.110
C
C  Now start Gaussian elimination
C
          DQW_DU(I,1) = RBTQ * DQW_DU(I,1)                    ! P244.102
          DTL_DV(I,1) = RBTQ * DTL_DV(I,1)                    ! P244.111
          AT_ATQ(I,1) = RBTQ * AT_ATQ(I,1)                    ! P244.112
C
C  Start calculating DELTAP_RML. Mid-level depths added in 2.2 below.
C
          DELTAP_RML(I) = DELTAP(I,1)
        ELSE ! No rapidly mixing layer calculations
          DTRDZ_RML(I) = 1.E30
          DQW_RML(I) = 1.E30
          DTL_RML(I) = 1.E30
          AQ_RML(I) = 1.E30
          AT_RML(I) = 1.E30
          DELTAP_RML(I) = 1.E30
        ENDIF
   21 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.2 Continue upward sweep through middle of the rapidly mixing layer
CL     (if it exists) and to its top. NB NRML is always <= BLM1.
C-----------------------------------------------------------------------
C
      DO 22 K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO 221 I=P1,P1+P_POINTS-1
C
C   If in the top rapidly mixing layer then do not include flux at its
C   top in the calculation ie FQW(I,NRML+1) and FTL(I,NRML+1) are not
C   included here; they are taken account of in the non-local mixing
C   through the "rapidly mixing layer".
C
          IF ( K .EQ. NRML(I) ) THEN
C
C   Add final DELTAP contribution to DELTAP_RML and then calculate
C   DTRDZ_RML.  Lower level contributions added in 2.1 and below.
C
            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
            DTRDZ_RML(I) = -DTIMEG / DELTAP_RML(I)
C
C  "Explicit" flux divergence across layer giving explicit
C  increment due to the local mixing at the top of rml.
C
            DQW_DU(I,K) = DTRDZ(I,K) * FQW(I,K)
            DTL_DV(I,K) = DTRDZ(I,K) * FTL(I,K)
C
C  Define matrix elements (A always zero for this case).
C
            CTQ = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)           ! P244.36
            RBTQ = 1.0 / ( 1.0  - CTQ*( 1.0 + AT_ATQ(I,KM1) ) )
C                                             ... Reciprocal of P244.113
C  Now start Gaussian elimination
C
            DQW_DU(I,K) = RBTQ * ( DQW_DU(I,K) - CTQ*DQW_DU(I,KM1) )
            DTL_DV(I,K) = RBTQ * ( DTL_DV(I,K) - CTQ*DTL_DV(I,KM1) )
          ELSEIF (K .LT. NRML(I)) THEN
C
C  Add layer depths to form total rml depth.
C
            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
C
C  "Explicit" flux divergence across layer giving explicit
C  increment due to the local mixing.P244.54/38
C
            DQW_DU(I,K) = -DTRDZ(I,K) * ( FQW(I,KP1) - FQW(I,K) )
            DTL_DV(I,K) = -DTRDZ(I,K) * ( FTL(I,KP1) - FTL(I,K) )
C
C  Define matrix elements.
C
            AT_ATQ(I,K) = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,KP1) ! P244.35
            CTQ = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)           ! P244.36
            RBTQ = 1.0 / ( 1.0  - AT_ATQ(I,K)
     &                          - CTQ * ( 1.0 + AT_ATQ(I,KM1) ) )
C                                             ... Reciprocal of P244.113
C  Now start Gaussian elimination
C
            DQW_DU(I,K) = RBTQ * ( DQW_DU(I,K) - CTQ*DQW_DU(I,KM1) )
            DTL_DV(I,K) = RBTQ * ( DTL_DV(I,K) - CTQ*DTL_DV(I,KM1) )
            AT_ATQ(I,K) = RBTQ * AT_ATQ(I,K)                 ! P244.115
          ENDIF
  221   CONTINUE
   22 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.3 Downward sweep through matrix. Add implicit increments due to
CL     local mixing within the rapidly mixing layer.  Update fluxes of
CL     heat and moisture at the surface and the top-of-rml using
CL     local mixing increments for layers 1 and NRML respectively.
C-----------------------------------------------------------------------
C
      DO 23 K=BLM1,1,-1
        KP1 = K + 1
        DO 231 I=P1,P1+P_POINTS-1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRML(I))) THEN
            QW(I,K) = QW(I,K) + DQW_DU(I,K)                   ! P244.128
            TL(I,K) = TL(I,K) + DTL_DV(I,K)                   ! P244.127
            FQW(I,KP1) = FQW(I,KP1)
     &                    + GAMMA_RHOKH_RDZ(I,KP1)*DQW_DU(I,K)
            FTL(I,KP1) = FTL(I,KP1)
     &                    + GAMMA_RHOKH_RDZ(I,KP1)*DTL_DV(I,K)
          ELSEIF ((NRML(I) .GE. 2) .AND. (K .LT. NRML(I))) THEN
            DQW_DU(I,K) = DQW_DU(I,K)
     &                           - AT_ATQ(I,K)*DQW_DU(I,KP1)  ! P244.???
            DTL_DV(I,K) = DTL_DV(I,K)
     &                           - AT_ATQ(I,K)*DTL_DV(I,KP1)  ! P244.???
            QW(I,K) = QW(I,K) + DQW_DU(I,K)                   ! P244.128
            TL(I,K) = TL(I,K) + DTL_DV(I,K)                   ! P244.127
          ENDIF
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. 1)) THEN
            IF (LAND_MASK(I)) THEN
              FTL(I,1) = FTL(I,1) - GAMMA_RHOKH_1(I) * DTL_DV(I,1)
              EA(I) = EA(I) - GAMMA_RHOKEA(I)*DQW_DU(I,1)     ! P244.122
              ESL(I) = ESL(I) - GAMMA_RHOKESL(I)*DQW_DU(I,1)  ! P244.124
              ES(I) = ES(I) - GAMMA_RHOKES(I)*DQW_DU(I,1)     ! P244.123
              FQW(I,1) = EA(I) + ES(I)                        ! P244.125
            ELSEIF (ICE_FRACT(I) .GT. 0.0) THEN
              FQW_ICE(I) = FQW(I,1) - FQW_LEAD(I)
     &                      - GAMMA_RHOKE(I) * ICE_FRACT(I)*DQW_DU(I,1)
              FTL_ICE(I) = FTL(I,1) - FTL_LEAD(I)
     &                      - GAMMA_RHOKH_1(I)*ICE_FRACT(I)*DTL_DV(I,1)
              FQW_LEAD(I) = FQW_LEAD(I) - GAMMA_RHOKE(I)
     &                        * ( 1.0 - ICE_FRACT(I) ) * DQW_DU(I,1)
              FTL_LEAD(I) = FTL_LEAD(I) - GAMMA_RHOKH_1(I)
     &                        * ( 1.0 - ICE_FRACT(I) ) * DTL_DV(I,1)
              FTL(I,1) = FTL_ICE(I) + FTL_LEAD(I)
              FQW(I,1) = FQW_ICE(I) + FQW_LEAD(I)
            ELSE ! ordinary sea point
              FTL(I,1) = FTL(I,1) - GAMMA_RHOKH_1(I)*DTL_DV(I,1)
              FQW(I,1) = FQW(I,1) - GAMMA_RHOKE(I)*DQW_DU(I,1)
            ENDIF
          ENDIF
  231   CONTINUE
   23 CONTINUE
C
C-----------------------------------------------------------------------
CL 3.  Calculate those matrix and vector elements on the LHS of eqn
CL     P244.79 which are to do with implicit solution of the moisture
CL     transport problem at the surface, above the rml (if it exists)
CL     and between all levels if it does not.
CL     Begin with "upward sweep" through lower half of matrix).
C-----------------------------------------------------------------------
CL 3.1 Row of matrix applying to QW transport into top "boundary"
CL     layer of model atmosphere.
C-----------------------------------------------------------------------
C
      DO 31 I=P1,P1+P_POINTS-1
        DQW_DU(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FQW(I,BL_LEVELS)
C                                                            ... P244.58
C
        AQ_AM(I,BL_LEVELS) = -DTRDZ(I,BL_LEVELS)               ! P244.56
     +                        * GAMMA_RHOKH_RDZ(I,BL_LEVELS)
C
        RBQ = 1.0 / ( 1.0 - AQ_AM(I,BL_LEVELS) )    ! Reciprocal P244.98
        DQW_DU(I,BL_LEVELS) = RBQ * DQW_DU(I,BL_LEVELS)        ! P244.99
        AQ_AM(I,BL_LEVELS) = RBQ * AQ_AM(I,BL_LEVELS)         ! P244.100
   31 CONTINUE
C
C-----------------------------------------------------------------------
CL 3.2 Rows of matrix applying to "middle of boundary layer" model
CL     layers, i.e. all but the topmost and bottom layers.
C-----------------------------------------------------------------------
C
      DO 32 K=BLM1,2,-1
        KP1=K+1
        DO 321 I=P1,P1+P_POINTS-1
C
C  "Explicit" flux divergence across layer giving explicit QW increment.
C
          IF ( K .GT. NRML(I) ) THEN
            DQW_DU(I,K) = -DTRDZ(I,K) * ( FQW(I,KP1) - FQW(I,K) )
C                                                            ... P244.54
C
            CQ = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,KP1)          ! P244.52
            AQ_AM(I,K) = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)    ! P244.51
            RBQ = 1.0 / ( 1.0 - AQ_AM(I,K) - CQ*( 1.0 + AQ_AM(I,KP1) ) )
C                       1                       2                    2 1
C                                             ... reciprocal of P244.101
C
            DQW_DU(I,K) = RBQ * ( DQW_DU(I,K) - CQ*DQW_DU(I,KP1) )
C                                                           ... P244.102
C
            AQ_AM(I,K) = RBQ * AQ_AM(I,K)                     ! P244.103
          ENDIF
  321   CONTINUE
   32 CONTINUE
C
C-----------------------------------------------------------------------
CL 3.3 Bottom model layer QW row of matrix equation.
C-----------------------------------------------------------------------
C
      DO 33 I=P1,P1+P_POINTS-1
        IF ( NRML(I) .GE. 2 ) THEN
C
C-----------------------------------------------------------------------
CL 3.3.1 Start calculating rapidly mixing layer increments.
C-----------------------------------------------------------------------
C
          NRMLP1 = NRML(I) + 1
C
C  "Explicit" QW increment for the rapidly mixing layer.
C
          DQW_RML(I) = -DTRDZ_RML(I) * ( FQW(I,NRMLP1) - FQW(I,1) )
C
C  Define coefficients A,B,C, for implicit calculations.
C
          AQ_RML(I) = -DTRDZ_RML(I) * GAMMA_RHOKE(I)
          CQ_RML = -DTRDZ_RML(I) * GAMMA_RHOKH_RDZ(I,NRMLP1)
          RBQ_RML = 1.0 / ( 1.0 - AQ_RML(I)
     &                      - CQ_RML * ( 1.0 + AQ_AM(I,NRMLP1) ) )
          DQW_RML(I) = RBQ_RML * ( DQW_RML(I)
     &                                 - CQ_RML * DQW_DU(I,NRMLP1) )
          AQ_RML(I) = RBQ_RML * AQ_RML(I)
        ELSE
C
C  "Explicit" increment for QW(1) when there is no rapidly mixing
C  layer or when it is only one model layer in depth.
C
          DQW_DU(I,1) = -DTRDZ(I,1) * ( FQW(I,2) - FQW(I,1) )  ! P244.49
          AQ_AM(I,1) = -DTRDZ(I,1) * GAMMA_RHOKE(I)            ! P244.46
          CQ = -DTRDZ(I,1) * GAMMA_RHOKH_RDZ(I,2)              ! P244.47
          RBQ = 1.0 / ( 1.0 - AQ_AM(I,1) - CQ*( 1.0 + AQ_AM(I,2) ) )
C                     1                       2                  2 1
C                                             ... reciprocal of P244.104
C
          DQW_DU(I,1) = RBQ * ( DQW_DU(I,1) - CQ*DQW_DU(I,2) )! P244.105
          AQ_AM(I,1) = RBQ * AQ_AM(I,1)                       ! P244.106
        ENDIF
C
C-----------------------------------------------------------------------
CL 4.  Calculate those matrix and vector elements on the LHS of eqn
CL     P244.79 which are to do with implicit solution of the heat
CL     transport problem (i.e. for surface temperature and ice/liquid
CL     water temperatures in atmospheric boundary layer), and begin the
CL     solution algorithm (perform "upward sweep" through upper half of
CL     matrix).
C-----------------------------------------------------------------------
CL 4.1 Surface temperature row of matrix (different treatments for
CL     different surface types).
C      ALPHAS used as a temporary store for estimate of TSTAR.
C-----------------------------------------------------------------------
C
        IF (LAND_MASK(I)) THEN
          IF (LYING_SNOW(I).GT.0.0) THEN
            LAT_HEAT = LS
          ELSE
            LAT_HEAT = LC
          ENDIF
          DTSTAR(I) = ASOIL_1(I) * ( RADNET(I)
     &                               - CP*FTL(I,1) - LAT_HEAT*FQW(I,1)
     &                               - SOIL_HT_FLUX(I) )
C           ... P244.15 (Explicit surface temperature increment;
C                        note that ASOIL_1 contains the TIMESTEP factor)
C
          ALPHAS(I) = TSTAR_NL(I) + DTSTAR(I)
C           ...(Estimate of TSTAR at timelevel n+1 based on the
C               explicit surface temperature increment)
C
        ELSEIF (ICE_FRACT(I) .GT. 0.0) THEN
          FTL_ICE(I) = FTL(I,1) - FTL_LEAD(I)
          FQW_ICE(I) = FQW(I,1) - FQW_LEAD(I)
          DTSTAR(I) = TIMESTEP * AI * ( RADNET(I)
     &                              - CP*FTL_ICE(I) - LS*FQW_ICE(I)
     &                              - SEA_ICE_HTF(I) )
C         ... P244.19 (Gridbox mean explicit surface temp. increment)
C
          ALPHAS(I) = TSTAR_NL(I) + DTSTAR(I)/ICE_FRACT(I)
C          ...(Estimate of T*(I), surface temperature of sea-ice, at
C              timelevel n+1 based on the explicit surface
C              temperature increment.)
C
        ELSE   ! Sea without sea-ice.
          DTSTAR(I) = 0.0   ! Surface temperature not updated.
          ALPHAS(I) = TSTAR_NL(I) + DTSTAR(I)
C          ...(Estimate of TSTAR at timelevel n+1 based on the
C              explicit surface temperature increment)
C
        ENDIF
   33 CONTINUE
      CALL QSAT (ATSTAR(P1),TSTAR_NL(P1),PSTAR(P1),P_POINTS)
C       ...Qsat for T*n (T*n(I) at sea-ice points)  temporarily
C          stored in ATSTAR.
C
      CALL QSAT (AT_ATQ(P1,1),ALPHAS(P1),PSTAR(P1),P_POINTS)
C       ...Qsat for (T*n + explicit deltaT*), or (T*n(I) +
C          explicit deltaT*(I) at sea-ice points), temporarily
C          storedin AT_ATQ(*,1).
C
      DO 41 I = P1,P1+P_POINTS-1
        IF (LAND_MASK(I)) THEN
          IF (LYING_SNOW(I).GT.0.0) THEN
            LAT_HEAT = LS
          ELSE
            LAT_HEAT = LC
          ENDIF
          IF (DTSTAR(I) .GE. -2.0 .AND. DTSTAR(I) .LE. 2.0) THEN
            ALPHAS(I)=(EPSILON*LAT_HEAT*ATSTAR(I))
     +                  / (R*TSTAR_NL(I)*TSTAR_NL(I))
C          ...P244.13a(Clausius-Clapeyron formula for partial dQsat/dT*
C                      used for small explicit T* increments.)
C
          ELSE
            ALPHAS(I) = ( AT_ATQ(I,1) - ATSTAR(I) ) / DTSTAR(I)
C           ...P244.13b(partial dQsat/dT* calculated as a finite
C               difference for large explicit T* increments.)
C
          ENDIF
          ATSTAR(I) = -ASOIL_1(I) * CP * GAMMA_RHOKH_1(I)     ! P244.16
          CT = -GAMMA_RHOKE(I) * ASOIL_1(I) * LAT_HEAT
C                   ... P244.17 (Note: ASOIL_1 includes factor TIMESTEP)
C
          IF ( NRML(I) .GE. 2 ) THEN
            RBT = 1.0 / ( 1.0 - ATSTAR(I) - CT * ALPHAS(I)
     &                                      * ( 1.0 + AQ_RML(I) ) )
            DTSTAR(I) = RBT * ( DTSTAR(I) - CT * DQW_RML(I) )
          ELSE
            RBT = 1.0 / ( 1.0 - ATSTAR(I) - CT * ALPHAS(I)
     &                                         * ( 1.0 + AQ_AM(I,1) ) )
C                       1                        2                  2 1
C                                             ... Reciprocal of P244.107
C
            DTSTAR(I) = RBT * ( DTSTAR(I) - CT * DQW_DU(I,1) )
C                                          ... P244.108 (giving DTSTAR')
C
          ENDIF
          ATSTAR(I) = RBT * ATSTAR(I)                         ! P244.109
        ELSEIF (ICE_FRACT(I).GT.0.0) THEN
          DTRICE = DTSTAR(I)/ICE_FRACT(I)
          IF (DTRICE .GE. -2.0 .AND. DTRICE .LE. 2.0) THEN
            ALPHAS(I) = (EPSILON*LS*ATSTAR(I))
     +                  / (R*TSTAR_NL(I)*TSTAR_NL(I))
C          ...P244.13a (Clausius-Clapeyron formula for partial
C             dQsat/dT*(I) used for small explicit T*(I) increments.)
C
          ELSE
            ALPHAS(I) = ( AT_ATQ(I,1) - ATSTAR(I) ) / DTRICE
C           ...P244.13b (partial dQsat/dT*(I) calculated as a finite
C               difference used for large explicit T*(I) increments.)
C
          ENDIF
          ATSTAR(I) = -GAMMA_RHOKH_1(I) * TIMESTEP * AI * CP   ! P244.20
          CT = -GAMMA_RHOKE(I) * TIMESTEP * AI * LS            ! P244.21
          IF ( NRML(I) .GE. 2 ) THEN
            RBT = 1.0 / ( 1.0 -ATSTAR(I) - CT * ALPHAS(I) *
     &                          ( 1.0 + ICE_FRACT(I) * AQ_RML(I) ) )
            DTSTAR(I) = RBT * ( DTSTAR(I) - ICE_FRACT(I)
     &                                       * CT * DQW_RML(I) )
          ELSE
            RBT = 1.0 / ( 1.0 - ATSTAR(I) - CT * ALPHAS(I)
     &                            * ( 1.0 + ICE_FRACT(I)*AQ_AM(I,1) ) )
C                       1           2                               2 1
C                                             ... Reciprocal of P2440.14
C
            DTSTAR(I) = RBT * ( DTSTAR(I) - ICE_FRACT(I)*CT*DQW_DU(I,1))
C                                           ... P2440.15, giving DTSTAR'
C
          ENDIF
          ATSTAR(I) = RBT * ATSTAR(I) * ICE_FRACT(I)          ! P2440.16
        ELSE               ! Sea with no sea-ice.
          ALPHAS(I) = (EPSILON*LC*ATSTAR(I))
     +                  / (R*TSTAR_NL(I)*TSTAR_NL(I))
C          ...P244.13a(Clausius-Clapeyron formula for partial dQsat/dT*)
          ATSTAR(I) = 0.0
        ENDIF
C-----------------------------------------------------------------------
CL 4.2 Lowest atmospheric layer TL row of matrix.
C-----------------------------------------------------------------------
        IF (NRML(I) .GE. 2) THEN
          NRMLP1 = NRML(I) + 1
C
C  "Explicit" rapidly mixing layer increment for TL.
C
          DTL_RML(I) = -DTRDZ_RML(I) * ( FTL(I,NRMLP1) - FTL(I,1) )
          AT_RML(I) = -DTRDZ_RML(I) * GAMMA_RHOKH_RDZ(I,NRMLP1)
          CT_RML = -DTRDZ_RML(I) * GAMMA_RHOKH_1(I)
          RBT_RML = 1.0 / ( 1.0 - AT_RML(I) - CT_RML
     &                                         * ( 1.0 + ATSTAR(I) ) )
          DTL_RML(I) = RBT_RML * ( DTL_RML(I)
     &                                 - CT_RML * DTSTAR(I) )
          AT_RML(I) = RBT_RML * AT_RML(I)
        ELSE
C
C  "Explicit" increment to TL(1) when there is no rapidly mixing layer
C  or it does not extend beyond the bottom model layer.
C
          DTL_DV(I,1) = -DTRDZ(I,1) * ( FTL(I,2) - FTL(I,1) )  ! P244.31
          CT = -DTRDZ(I,1) * GAMMA_RHOKH_1(I)                  ! P244.29
          AT_ATQ(I,1) = -DTRDZ(I,1) * GAMMA_RHOKH_RDZ(I,2)     ! P244.28
          RBT = 1.0 / ( 1.0 - AT_ATQ(I,1) - CT*( 1.0 + ATSTAR(I) ) )
C                     1                        2                 2 1
C                                             ... Reciprocal of P244.110
C
          DTL_DV(I,1) = RBT * ( DTL_DV(I,1) - CT*DTSTAR(I) )  ! P244.111
          AT_ATQ(I,1) = RBT * AT_ATQ(I,1)                     ! P244.112
        ENDIF
   41 CONTINUE
C-----------------------------------------------------------------------
CL 4.3 Rows of matrix applying to TL transport into model layers in the
CL     "middle" of the "boundary" layer, i.e. all but the bottom layer
CL     and the top "boundary" layer.
C-----------------------------------------------------------------------
      DO 43 K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO 431 I=P1,P1+P_POINTS-1
C
C   "Explicit" flux divergence across layer giving explicit TL increment
C   due to mixing above rml if it exists or for all levels if it does
C   not.
C
          NRMLP1 = NRML(I) + 1
          IF (K .GT. NRML(I)) THEN
            DTL_DV(I,K) = -DTRDZ(I,K) * ( FTL(I,KP1) - FTL(I,K) )
C                                                            ... P244.38
C
            AT_ATQ(I,K) = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,KP1) ! P244.35
            CT = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)            ! P244.36
            IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
              RBT = 1.0 / ( 1.0 - AT_ATQ(I,K) - CT*( 1.0 + AT_RML(I) ) )
              DTL_DV(I,K) = RBT * ( DTL_DV(I,K) - CT*DTL_RML(I) )
            ELSE
              RBT = 1.0 / ( 1.0 - AT_ATQ(I,K)
     &                          - CT*( 1.0 + AT_ATQ(I,KM1) ) )
C                         1          2                     2 1
C                                             ... Reciprocal of P244.113
C
              DTL_DV(I,K) = RBT * ( DTL_DV(I,K) - CT*DTL_DV(I,KM1) )
C                                                           ... P244.114
C
            ENDIF
            AT_ATQ(I,K) = RBT * AT_ATQ(I,K)                   ! P244.115
          ENDIF
  431   CONTINUE
   43 CONTINUE
C-----------------------------------------------------------------------
CL 4.4 Top "boundary" layer TL row of matrix.  TL for this layer can
CL     then be, and is, updated.
C-----------------------------------------------------------------------
      DO 44 I=P1,P1+P_POINTS-1
        DTL_DV(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FTL(I,BL_LEVELS)
C                                                            ... P244.42
C
        CT = -DTRDZ(I,BL_LEVELS) * GAMMA_RHOKH_RDZ(I,BL_LEVELS)! P244.40
        IF (NRML(I) .EQ. BLM1) THEN
          RBT = 1.0 / ( 1.0 - CT*( 1.0 + AT_RML(I) ) )
          DTL_DV(I,BL_LEVELS) = RBT * ( DTL_DV(I,BL_LEVELS)
     +                                  - CT*DTL_RML(I) )
        ELSE
         RBT = 1.0 / ( 1.0 - CT*( 1.0 + AT_ATQ(I,BLM1) ) )
C                    1          2                      2 1
C                                             ... Reciprocal of P244.116
C
          DTL_DV(I,BL_LEVELS) = RBT * ( DTL_DV(I,BL_LEVELS)   ! P244.117
     +                                  - CT*DTL_DV(I,BLM1) )
        ENDIF
        TL(I,BL_LEVELS) = TL(I,BL_LEVELS) + DTL_DV(I,BL_LEVELS)
C                                                           ... P244.127
C
   44 CONTINUE
C
C-----------------------------------------------------------------------
CL 5.  "Downward sweep" through whole matrix.  TL, QW and TSTAR are
CL     updated when the final implicit increments have been calculated.
C-----------------------------------------------------------------------
CL 5.1 Remaining TL rows of matrix and add implicit increments above
CL     the rml or at all levels if it is less than two layers thick.
C-----------------------------------------------------------------------
C
      DO 51 K=BLM1,1,-1
        DO 511 I=P1,P1+P_POINTS-1
          IF ( (K .GT. NRML(I)) .OR. (NRML(I) .LT. 2) ) THEN
            DTL_DV(I,K) = DTL_DV(I,K) - AT_ATQ(I,K)*DTL_DV(I,K+1)
                                                              ! P244.118
            TL(I,K) = TL(I,K) + DTL_DV(I,K)                   ! P244.127
          ENDIF
  511   CONTINUE
   51 CONTINUE
      DO 52 I=P1,P1+P_POINTS-1
C
C-----------------------------------------------------------------------
CL 5.2 Surface temperature (TSTAR) row of matrix and rapidly mixing
CL     layer increments.
C-----------------------------------------------------------------------
        IF ( NRML(I) .GE. 2 ) THEN
          NRMLP1 = NRML(I) + 1
          DTL_RML(I) = DTL_RML(I) - AT_RML(I) * DTL_DV(I,NRMLP1)
          DTSTAR(I) = DTSTAR(I) - ATSTAR(I) * DTL_RML(I)
          DQW_RML(I) = DQW_RML(I)
     &                      - AQ_RML(I) * ALPHAS(I) * DTSTAR(I)
          TL(I,1) = TL(I,1) + DTL_RML(I)
          QW(I,1) = QW(I,1) + DQW_RML(I)
        ELSE
          DTSTAR(I) = DTSTAR(I) - ATSTAR(I)*DTL_DV(I,1)       ! P244.119
C
C-----------------------------------------------------------------------
CL 5.3 Lowest-level QW row of matrix; local mixing where there is no rml
C-----------------------------------------------------------------------

          DQW_DU(I,1) = DQW_DU(I,1) - AQ_AM(I,1)*ALPHAS(I)*DTSTAR(I)
C                                                           ... P244.120
          QW(I,1) = QW(I,1) + DQW_DU(I,1)     ! P244.128
          DQW_1(I) = DQW_DU(I,1)
        ENDIF
        TSTAR(I) = TSTAR(I) + DTSTAR(I)
   52 CONTINUE
C
C-----------------------------------------------------------------------
CL 5.4 Remaining QW rows of matrix + updating of QW's.
CL     Add implicit increments due to mixing above rml or at all levels
CL     if there it does not exist.
C-----------------------------------------------------------------------
      DO 54 K=2,BL_LEVELS
        DO 541 I=P1,P1+P_POINTS-1
          IF (K .GT. NRML(I)) THEN
            IF ((NRML(I) .GE. 2) .AND. (K-1 .EQ. NRML(I))) THEN
              DQW_DU(I,K) = DQW_DU(I,K) - AQ_AM(I,K)*DQW_RML(I)
            ELSE
              DQW_DU(I,K) = DQW_DU(I,K) - AQ_AM(I,K)*DQW_DU(I,K-1)
C                                                           ... P244.121
            ENDIF
            QW(I,K) = QW(I,K) + DQW_DU(I,K)                  ! P244.128
          ELSE
C
C  Add the increments due to rapid mixing if in the rapidly mixing layer
C
            TL(I,K) = TL(I,K) + DTL_RML(I)
            QW(I,K) = QW(I,K) + DQW_RML(I)
          ENDIF
  541   CONTINUE
   54 CONTINUE
C
C-----------------------------------------------------------------------
CL 6.  Calculate final implicit fluxes of heat and moisture.
C-----------------------------------------------------------------------
CL 6.1 Surface fluxes for the 3 surface types: land, sea-ice, ordinary
CL     sea. Pass out the value of RHOKH(,1) in GAMMA_RHOKH_1 to satisfy
CL     STASH GAMMA_RHOKH_RDZ will contain precisely that on output.
C-----------------------------------------------------------------------
C
      DO 61 I=P1,P1+P_POINTS-1
        IF (LAND_MASK(I)) THEN
          FTL_ICE(I) = 0.0
          H_SEA(I) = 0.0
          FQW_ICE(I) = 0.0
          E_SEA(I) = 0.0
          IF ( NRML(I) .GE. 2 ) THEN
            FTL(I,1) = FTL(I,1) - GAMMA_RHOKH_1(I) *
     &                             ( DTL_RML(I) - DTSTAR(I) )
            DELTDQ = DQW_RML(I) - ALPHAS(I) * DTSTAR(I)
          ELSE
            FTL(I,1) = FTL(I,1) - GAMMA_RHOKH_1(I) *
     &                             ( DTL_DV(I,1) - DTSTAR(I) )
C                 ... P244.11 (for FTL = H/Cp, rather than for H itself)
C
            DELTDQ = DQW_DU(I,1) - ALPHAS(I)*DTSTAR(I)        ! P244.126
          ENDIF
          EA(I) = EA(I) - GAMMA_RHOKEA(I)*DELTDQ              ! P244.122
          ES(I) = ES(I) - GAMMA_RHOKES(I)*DELTDQ              ! P244.123
          ESL(I) = ESL(I) - GAMMA_RHOKESL(I)*DELTDQ           ! P244.124
          FQW(I,1) = EA(I) + ES(I)                            ! P244.125
        ELSEIF (ICE_FRACT(I).GT.0.0) THEN
          IF ( NRML(I) .GE. 2 ) THEN
            H_SEA(I) = CP * ( FTL_LEAD(I) - GAMMA_RHOKH_1(I)
     &                      * ( 1.0 - ICE_FRACT(I) ) * DTL_RML(I) )
            FTL_ICE(I) = FTL_ICE(I) - GAMMA_RHOKH_1(I)
     &                    * ( ICE_FRACT(I)*DTL_RML(I) - DTSTAR(I) )
            E_SEA(I) = FQW_LEAD(I) - GAMMA_RHOKE(I)
     &               * (1.0 - ICE_FRACT(I)) * DQW_RML(I)
            DELTDQ = ICE_FRACT(I)*DQW_RML(I) - ALPHAS(I)*DTSTAR(I)
          ELSE
            H_SEA(I) = CP * ( FTL_LEAD(I) - GAMMA_RHOKH_1(I)
     &                     * ( 1.0 - ICE_FRACT(I) ) * DTL_DV(I,1) )
            FTL_ICE(I) = FTL_ICE(I) - GAMMA_RHOKH_1(I)
     &             * ( ICE_FRACT(I)*DTL_DV(I,1) - DTSTAR(I) ) ! P2440.1
            E_SEA(I) = FQW_LEAD(I) - GAMMA_RHOKE(I)
     &               * (1.0 - ICE_FRACT(I)) * DQW_DU(I,1)     ! P2440.2
            DELTDQ = ICE_FRACT(I)*DQW_DU(I,1) - ALPHAS(I)*DTSTAR(I)
          ENDIF
          FQW_ICE(I) = FQW_ICE(I) - GAMMA_RHOKE(I) * DELTDQ
          FTL(I,1) = FTL_ICE(I) + H_SEA(I)/CP
          FQW(I,1) = FQW_ICE(I) + E_SEA(I)
          EA(I) = FQW(I,1)
          ES(I) = 0.0
          ESL(I) = 0.0
        ELSE  ! ordinary sea point
          IF ( NRML(I) .GE. 2) THEN
            H_SEA(I) = CP*( FTL(I,1) - GAMMA_RHOKH_1(I)*DTL_RML(I) )
            E_SEA(I) = FQW(I,1) - GAMMA_RHOKE(I)*DQW_RML(I)
          ELSE
            H_SEA(I) = CP * ( FTL(I,1) - GAMMA_RHOKH_1(I)*DTL_DV(I,1) )
            E_SEA(I) = FQW(I,1) - GAMMA_RHOKE(I)*DQW_DU(I,1)
          ENDIF
          FTL_ICE(I) = 0.0
          FQW_ICE(I) = 0.0
          FTL(I,1) = H_SEA(I) / CP
          FQW(I,1) = E_SEA(I)
          EA(I) = FQW(I,1)
          ES(I) = 0.0
          ESL(I) = 0.0
        ENDIF
        GAMMA_RHOKH_1(I) = GAMMA_RHOKH_1(I) / GAMMA(1)
   61 CONTINUE
C
C-----------------------------------------------------------------------
CL 6.2 Fluxes at layer interfaces above the surface.
C-----------------------------------------------------------------------
      DO 62 K=2,BL_LEVELS
        KM1 = K-1
        DO 621 I=P1,P1+P_POINTS-1
C
C  Calculate and store fluxes due to local mixing.
C  FTL(local mixing) stored in array AT,
C  FQW(local mixing) stored in array AQ_AM.
C
          NRMLP1 = NRML(I) + 1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
            AT_ATQ(I,K) = FTL(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( DTL_DV(I,K) - DTL_RML(I) )
            AQ_AM(I,K) = FQW(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( DQW_DU(I,K) - DQW_RML(I) )
          ELSE
            AT_ATQ(I,K) = FTL(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( DTL_DV(I,K) - DTL_DV(I,KM1) )
            AQ_AM(I,K) = FQW(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( DQW_DU(I,K) - DQW_DU(I,KM1) )
          ENDIF
C
C  Now calculate the implicit fluxes including both local mixing and
C  if appropriate also the fluxes due to rapid mixing through layers.
C
          IF ( K .EQ. 2 ) THEN
            IF ( NRML(I) .GE. 2 ) THEN
              FTL(I,K) = AT_ATQ(I,K)
     &                    + FTL(I,KM1) - DTL_RML(I) / DTRDZ(I,KM1)
              FQW(I,K) = AQ_AM(I,K)
     &                    + FQW(I,KM1) - DQW_RML(I) / DTRDZ(I,KM1)
            ELSE
              FTL(I,K) = AT_ATQ(I,K)
              FQW(I,K) = AQ_AM(I,K)
            ENDIF
          ELSEIF ( K .LE. NRML(I) ) THEN
            FTL(I,K) = AT_ATQ(I,K) - AT_ATQ(I,KM1)
     &                    + FTL(I,KM1) - DTL_RML(I) / DTRDZ(I,KM1)
            FQW(I,K) = AQ_AM(I,K) - AQ_AM(I,KM1)
     &                    + FQW(I,KM1) - DQW_RML(I) / DTRDZ(I,KM1)
          ELSE
            FTL(I,K) = AT_ATQ(I,K)
            FQW(I,K) = AQ_AM(I,K)
          ENDIF
  621   CONTINUE
   62 CONTINUE
C
CL----------------------------------------------------------------------
CL (B) Calculations on UV-grid.
CL----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL 7.  Solve matrix equation P244.80 for implicit increments to U and V.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  7.0 Interpolate to UV-grid the pressure changes across the model
C      layers
C-----------------------------------------------------------------------
      DO 70 K=1,BL_LEVELS
        CALL P_TO_UV(DELTAP(P1,K),AQ_RML(U1+ROW_LENGTH),
     +     P_POINTS,P_POINTS-ROW_LENGTH,ROW_LENGTH,P_ROWS)
        DO 701 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DELTAP(I,K) = AQ_RML(I)
  701   CONTINUE
   70 CONTINUE
C
C-----------------------------------------------------------------------
CL 7.1 Initial calculations and "upward sweep".
CL (a) "Surface" fluxes.
C-----------------------------------------------------------------------
C
      DO 71 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
        DTRDZ_UV(I,1) = -DTIMEG / DELTAP(I,1)
        DELTAP_RML(I) = 0.0
C
C  "Explicit" increments to U(1) and V(1) when there is no rapidly
C  mixing layer or it does not extend beyond the bottom model layer.
C
        DQW_DU(I,1) = DTRDZ_UV(I,1) * ( TAUX(I,2) - TAUX(I,1) )
C                                                            ... P244.67
C
        DTL_DV(I,1) = DTRDZ_UV(I,1) * ( TAUY(I,2) - TAUY(I,1) )
C                                                            ... P244.67
C
        CM = -DTRDZ_UV(I,1) * GAMMA_RHOKM_1(I)               ! P244.66
        AQ_AM(I,1) = -DTRDZ_UV(I,1) * GAMMA_RHOKM_RDZUV(I,2)   ! P244.64
        RBM = 1.0 / ( 1.0 - AQ_AM(I,1) - CM )    ! Reciprocal of P244.81
        DQW_DU(I,1) = RBM * DQW_DU(I,1)                        ! P244.82
        DTL_DV(I,1) = RBM * DTL_DV(I,1)                        ! P244.83
        AQ_AM(I,1) = RBM * AQ_AM(I,1)                          ! P244.84
   71 CONTINUE
C
C
C-----------------------------------------------------------------------
CL (b) Fluxes at (or rows representing) layer interfaces above the
CL     surface but below the top of the boundary layer.
C-----------------------------------------------------------------------
      DO 72 K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO 721 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DTRDZ_UV(I,K) = -DTIMEG / DELTAP(I,K)
            DQW_DU(I,K) = DTRDZ_UV(I,K) * ( TAUX(I,KP1) - TAUX(I,K) )
C                                                            ... P244.74
            DTL_DV(I,K) = DTRDZ_UV(I,K) * ( TAUY(I,KP1) - TAUY(I,K) )
C                                                            ... P244.74
            AQ_AM(I,K) = -DTRDZ_UV(I,K) * GAMMA_RHOKM_RDZUV(I,KP1)
C                                                            ... P244.71
          CM = -DTRDZ_UV(I,K) * GAMMA_RHOKM_RDZUV(I,K)         ! P244.72
            RBM = 1.0 / ( 1.0 - AQ_AM(I,K) -CM*( 1.0 + AQ_AM(I,KM1) ) )
C                     1                      2                    2 1
C                                              ... Reciprocal of P244.85
C
            DQW_DU(I,K) = RBM * ( DQW_DU(I,K) - CM*DQW_DU(I,KM1) )
C                                                            ... P244.86
C
            DTL_DV(I,K) = RBM * ( DTL_DV(I,K) - CM*DTL_DV(I,KM1) )
C                                                            ... P244.87
C
          AQ_AM(I,K) = RBM * AQ_AM(I,K)                        ! P244.88
  721   CONTINUE
   72 CONTINUE
C-----------------------------------------------------------------------
CL (c) Top "boundary" layer; also increment U and V here, as implicit
CL     flux for this layer is got from "upward sweep" alone.
C-----------------------------------------------------------------------
      DO 73 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
        DTRDZ_UV(I,BL_LEVELS) = -DTIMEG / DELTAP(I,BL_LEVELS)
        DQW_DU(I,BL_LEVELS) = -DTRDZ_UV(I,BL_LEVELS) * TAUX(I,BL_LEVELS)
C                                                            ... P244.78
C
        DTL_DV(I,BL_LEVELS) = -DTRDZ_UV(I,BL_LEVELS) * TAUY(I,BL_LEVELS)
C                                                            ... P244.78
C
        CM = -DTRDZ_UV(I,BL_LEVELS) * GAMMA_RHOKM_RDZUV(I,BL_LEVELS)
C                                                            ... P244.76
        RBM = 1.0 / ( 1.0 - CM*( 1.0 + AQ_AM(I,BLM1) ) )
C                     1          2                     2 1
C                                              ... Reciprocal of P244.89
C
        DQW_DU(I,BL_LEVELS) = RBM * ( DQW_DU(I,BL_LEVELS)    ! P244.90
     +                                    - CM*DQW_DU(I,BLM1) )
        DTL_DV(I,BL_LEVELS) = RBM * ( DTL_DV(I,BL_LEVELS)    ! P244.91
     +                                    - CM*DTL_DV(I,BLM1) )
        U(I,BL_LEVELS) = U(I,BL_LEVELS) + DQW_DU(I,BL_LEVELS)  ! P244.94
        V(I,BL_LEVELS) = V(I,BL_LEVELS) + DTL_DV(I,BL_LEVELS)  ! P244.95
   73 CONTINUE

C-----------------------------------------------------------------------
CL 7.4 Complete solution of matrix equation by performing "downward
CL     sweep", then update U and V.
C-----------------------------------------------------------------------
      DO 74 K=BLM1,1,-1
        KP1 = K+1
        DO 741 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DQW_DU(I,K) = DQW_DU(I,K) - AQ_AM(I,K)*DQW_DU(I,KP1) ! P244.92
          DTL_DV(I,K) = DTL_DV(I,K) - AQ_AM(I,K)*DTL_DV(I,KP1) ! P244.93
          U(I,K) = U(I,K) + DQW_DU(I,K)                        ! P244.94
          V(I,K) = V(I,K) + DTL_DV(I,K)                        ! P244.95
  741   CONTINUE
   74 CONTINUE
C
C-----------------------------------------------------------------------
CL 8.  Essentially diagnostic calculations, though some values (i.e. the
CL     surface wind stresses) are required by the coupled version of the
CL     model.
C-----------------------------------------------------------------------
CL 8.1 Surface wind stress components.
CL     Pass out the value of RHOKM(,1) in GAMMA_RHOKM_1 to satisfy STASH
CL     GAMMA_RHOKM_RDZ will contain precisely that on output.
C-----------------------------------------------------------------------
C
      DO 81 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
        TAUX(I,1) = TAUX(I,1) + GAMMA_RHOKM_1(I)*DQW_DU(I,1)
C                                                            ... P244.61
C
        TAUY(I,1) = TAUY(I,1) + GAMMA_RHOKM_1(I)*DTL_DV(I,1)
C                                                            ... P244.61
C
        GAMMA_RHOKM_1(I) = GAMMA_RHOKM_1(I) / GAMMA(1)
   81 CONTINUE
C
C  Set extreme rows to "missing data indicator".
C
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          TAUX(I,1)=1.E30
          TAUY(I,1)=1.E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1 + (U_ROWS-1)*ROW_LENGTH , U1 + U_ROWS*ROW_LENGTH - 1
          TAUX(I,1)=1.E30
          TAUY(I,1)=1.E30
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
CL 8.2 Wind stress components at layer interfaces above the surface.
C-----------------------------------------------------------------------
      DO 82 K=2,BL_LEVELS
        KM1 = K-1
        DO 821 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          AQ_AM(I,K) = TAUX(I,K) + GAMMA_RHOKM_RDZUV(I,K)    ! P244.61
     +                        *( DQW_DU(I,K) - DQW_DU(I,KM1) )
          AT_ATQ(I,K) = TAUY(I,K) + GAMMA_RHOKM_RDZUV(I,K)   ! P244.61
     +                        *( DTL_DV(I,K) - DTL_DV(I,KM1) )
          TAUX(I,K) = AQ_AM(I,K)
          TAUY(I,K) = AT_ATQ(I,K)
  821   CONTINUE
C
C  Set extreme rows to "missing data indicator".
C
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1 + (U_ROWS-1)*ROW_LENGTH , U1 + U_ROWS*ROW_LENGTH - 1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
      ENDIF
   82 CONTINUE
C
C-----------------------------------------------------------------------
CL 8.3 Wind components at 10 metres above the surface.
C-----------------------------------------------------------------------
      IF (SU10 .OR. SV10) THEN
CMIC$ DO ALL VECTOR SHARED(U_POINTS, ROW_LENGTH, U1, SU10, SV10, U_FIELD
CMIC$1   , U0, CDR10M, BL_LEVELS, U, U10M, V0, V, V10M) PRIVATE(I)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 83 I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          IF (SU10)
     +      U10M(I) = U0(I) + CDR10M(I)*( U(I,1) - U0(I) )     ! P244.96
          IF (SV10)
     +      V10M(I) = V0(I) + CDR10M(I)*( V(I,1) - V0(I) )     ! P244.97
   83   CONTINUE
C
C  Set extreme rows to "missing data indicator".
C
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          IF (SU10) THEN
            U10M(I)=1.E30
          ENDIF
          IF (SV10) THEN
            V10M(I)=1.E30
          ENDIF
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1 + (U_ROWS-1)*ROW_LENGTH , U1 + U_ROWS*ROW_LENGTH - 1
          IF (SU10) THEN
            U10M(I)=1.E30
          ENDIF
          IF (SV10) THEN
            V10M(I)=1.E30
          ENDIF
        ENDDO
      ENDIF
      ENDIF
  999 CONTINUE   ! Branch for error exit.
      IF (LTIMER) THEN
      CALL TIMER('IMPLCAL ',4)
      ENDIF
      RETURN
      END
