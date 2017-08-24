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
C*LL  SUBROUTINE KMKH---------------------------------------------------
CLL
CLL  Purpose: To calculate the turbulent mixing coefficients KM and KH
CLL           and then the "explicit" turbulent fluxes for layer
CLL           interfaces 1+1/2 to BL_LEVELS-1/2.
CLL
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  05/01/95   rho*Km before interpolation output via argument
CLL                   list.                            R.N.B.Smith
!LL   4.1  29/05/96   Remove U_ROWS check and add swapbounds for MPP
!LL                   code.  P.Burton
CLL
CLL   4.1  27/03/96   Correct declaration of ZLB.  S Jackson
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
!LL  4.2  08/01/97  Add SWAPBOUNDS for RHOKM.  RTHBarnes.
!LL   4.3  14/01/97   MPP code : Corrected setting of polar rows
!LL                                                     P.Burton
!LL  4.3  18/03/97  Protect swapbounds by *IF DEF,MPP.  RTHBarnes.
CLL   4.3  04/02/97   Number of model layers in turbulently mixed layer
CLL                   diagnosed and stored in NTML(*).
CLL                   NTML also passed to EX_COEF.
CLL                   Logical switches passed to EX_COEF.
CLL                                                        R.N.B.Smith
!LL  4.4  08/09/97  L_BL_LSPICE specifies mixed phase precipitation
!LL                 scheme                         D.Wilson
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  System component covered: Part of P243.
CLL
CLL  Documentation: UMDP No.24
CLL
CLL---------------------------------------------------------------------
C*
C*L  Arguments :-
      SUBROUTINE KMKH (
     + P_FIELD,U_FIELD,P1,U1,
     + P_POINTS,U_POINTS,ROW_LENGTH,P_ROWS,U_ROWS,BL_LEVELS
     +,TIMESTEP,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,CCA,BQ_1,BT_1,BF_1
     &,CF,DZL
     +,PSTAR,Q,QCF,QCL,RDZ,T,TV,U,U_P,V,V_P,Z0M,ZLB,H_BLEND
     +,FQW,FTL,TAUX,TAUY,QW,RHOKM,RHOKH,TL,ZH
     +,RDZUV,RHO_KM
     +,CCB,CCT,L_MOM,L_MIXLEN
     &,L_BL_LSPICE
     +,NRML
     +,ERROR,LTIMER
     +)

      IMPLICIT NONE
      LOGICAL LTIMER  ! IN Flag for TIMER diagnostics
     &,L_BL_LSPICE         ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
C
      LOGICAL
     & L_MOM       ! IN Switch for convective momentum transport.
     &,L_MIXLEN    ! IN Switch for reducing the turbulent mixing
C                  !    length above the top of the boundary layer.
C
      INTEGER
     + P_FIELD     ! IN No. of P-grid points in whole field.
     +,U_FIELD     ! IN No. of UV-grid points in whole field.
     +,P1          ! IN First P-grid point to be processed.
     +,U1          ! IN First UV-grid point to be processed.
     +,P_POINTS    ! IN No. of P-grid points to be processed.
     +,U_POINTS    ! IN No. of UV-grid points to be processed.
     +,ROW_LENGTH  ! IN No. of points in latitude row (inclusive of
C                  !    endpoints for limited area model).
     +,P_ROWS      ! IN No. of P-rows of data to be processed.
     +,U_ROWS      ! IN No. of UV-rows of data to be processed.
     +,BL_LEVELS   ! IN No. of atmospheric levels for which boundary
C                  !    layer fluxes are calculated.  Assumed <=30
C                  !    for dimensioning GAMMA() in common deck C_GAMMA
      REAL
     + TIMESTEP                  ! IN Timestep in seconds.
     +,AK(BL_LEVELS)             ! IN Hybrid "A" for boundary levels.
     +,BK(BL_LEVELS)             ! IN Hybrid "B" for boundary levels.
     +,AKH(BL_LEVELS)            ! IN Hybrid "A" for layer interfaces.
C                                !    AKH(K) is value for lower boundary
C                                !    of layer K.
     +,BKH(BL_LEVELS)            ! IN Hybrid "B" for layer interfaces.
C                                !    BKH(K) is value for lower boundary
C                                !    of layer K.
     +,DELTA_AK(BL_LEVELS)       ! IN Hybrid A-difference across layer.
     +,DELTA_BK(BL_LEVELS)       ! IN Hybrid B-difference across layer.
     +,CCA(P_FIELD)              ! IN Convective Cloud Amount.
     +,BQ_1(P_FIELD)             ! IN A buoyancy parameter for lowest
C                                !    atmospheric level (from routine
C                                !    SF_EXCH).
     +,BT_1(P_FIELD)             ! IN A buoyancy parameter for lowest
C                                !    atmospheric level (from routine
C                                !    SF_EXCH).
     &,BF_1(P_FIELD)             
!        IN A buoyancy parameter for lowest
!           atmospheric level (from routine SF_EXCH).
     +,CF(P_FIELD,BL_LEVELS)     ! IN Cloud fractions for boundary levs.
     +,DZL(P_FIELD,BL_LEVELS)    ! IN Layer depths (m).  DZL(,K) is the
C                                !    distance from layer boundary K-1/2
C                                !    to layer boundary K+1/2.  For K=1
C                                !    the lower boundary is the surface.
     +,PSTAR(P_FIELD)            ! IN Surface pressure (Pa).
     +,Q(P_FIELD,BL_LEVELS)      ! IN Sp humidity (kg water per kg air).
     +,QCF(P_FIELD,BL_LEVELS)    ! IN Cloud ice (kg per kg air).
     +,QCL(P_FIELD,BL_LEVELS)    ! IN Cloud liq water (kg per kg air).
      REAL                       ! Split to keep continuation cards < 19
     + RDZ(P_FIELD,BL_LEVELS)    ! IN Reciprocal of distance between
C                                !    hybrid levels (m-1).  1/RDZ(,K) is
C                                !    the vertical distance from level
C                                !    K-1 to level K, except that for
C                                !    K=1 it is just the height of the
C                                !    lowest atmospheric level.
     +,T(P_FIELD,BL_LEVELS)      ! IN Temperature (K).
     +,TV(P_FIELD,BL_LEVELS)     ! IN Virtual temperature (K) - from
C                                !    SUBROUTINE Z.
     +,U(U_FIELD,BL_LEVELS)      ! IN Westerly wind component on UV-grid
C                                !    (m per s).
     +,U_P(P_FIELD,BL_LEVELS)    ! IN Westerly wind component on P-grid
C                                !    (m per s).
     +,V(U_FIELD,BL_LEVELS)      ! IN Southerly wind component on
C                                !    UV-grid (m per s).
     +,V_P(P_FIELD,BL_LEVELS)    ! IN Southerly wind component on P-grid
C                                !    (m per s).
     +,Z0M(P_FIELD)              ! IN Roughness length for momentum (m).
     +,ZLB(P_FIELD,BL_LEVELS)  ! IN ZLB(,K) is height above surface of
C                                !    lower boundary of layer K (m).
     +,H_BLEND(P_FIELD)          ! IN Blending height for effective
C                                !    roughness scheme passed through
C                                !    EX_COEF
      REAL                       ! Note: for all these INOUT arrays,
C                                !       apart from ZH, level 1 is IN
C                                !       (though not always used), and
C                                !       the other levels are all OUT.
     + FQW(P_FIELD,BL_LEVELS)    ! INOUT "Explicit" flux of QW (i.e.
C                                !       evaporation) from layer below,
C                                !       on P-grid (kg per sq m per s).
     +,FTL(P_FIELD,BL_LEVELS)    ! INOUT "Explicit" flux of TL = H/CP
C                                !       (sensible heat/CP) from layer
C                                !       below, on P-grid.
     +,TAUX(U_FIELD,BL_LEVELS)   ! INOUT "Explicit" x-component of
C                                !       turbulent stress at level K-1/2
C                                !       On UV-grid with first and last
C                                !       rows set to "missing data".
     +,TAUY(U_FIELD,BL_LEVELS)   ! INOUT "Explicit" y-component of
C                                !       turbulent stress at level K-1/2
C                                !       On UV-grid with first and last
C                                !       rows set to "missing data".
     +,QW(P_FIELD,BL_LEVELS)     ! INOUT Total water content (kg per kg
C                                !       air).
     +,RHOKM(U_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
C                                !       exchange coefficient for
C                                !       momentum, on UV-grid with first
C                                !       and last rows set to "missing
C                                !       data".Output as GAMMA*RHOKM*
C                                !       RDZUV for P244 (IMPL_CAL).
     +,RHOKH(P_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
C                                !       exchange coefficient for FTL,
C                                !       on P-grid.Output as GAMMA*
C                                !       *RHOKH*RDZ for P244(IMPL_CAL)
     +,TL(P_FIELD,BL_LEVELS)     ! INOUT Liquid/frozen water temperature
C                                !       (K).
     +,ZH(P_FIELD)               ! INOUT Boundary layer height (m).
      INTEGER
     + CCB(P_FIELD)              ! IN  Convective Cloud Base.
     +,CCT(P_FIELD)              ! IN  Convective Cloud Top.
     +,NRML(P_FIELD)             ! INOUT Number of model layers in the
C                                !       unstable Rapidly Mixing Layer.
     +,ERROR                     ! OUT 1 if bad arguments; 0 otherwise.
      REAL
     + RDZUV(U_FIELD,2:BL_LEVELS)  ! OUT RDZ on UV-grid, with first and
C                                  !     last rows set to "missing
C                                  !     data".
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT RHO * KM before interpolation
C                                  !     to UV-grid.
C*
C*L---------------------------------------------------------------------
C    External references :-
      EXTERNAL P_TO_UV
      EXTERNAL TIMER
      EXTERNAL QSAT,QSAT_WAT,EX_COEF
C*
C*L---------------------------------------------------------------------
C    Local and other symbolic constants :-
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     + ETAR=1.0/(1.0-EPSILON)   ! Used in buoyancy parameter BETAC.
     +,GRCP=G/CP                ! Used in DZTL, FTL calculations.
     +,LCRCP=LC/CP              ! Latent heat of condensation / CP.
     +,LFRCP=LF/CP              ! Latent heat of fusion / CP.
     +,LS=LC+LF                 ! Latent heat of sublimation.
     +,LSRCP=LS/CP              ! Latent heat of sublimation / CP.
     +)
C*
C
C  Define local storage.
C
C  (a) Workspace. 6*BL_LEVELS-1 blocks of real workspace are required
C      plus 1 block of logical.
C
C
      LOGICAL
     + TOPBL(P_FIELD)            ! Flag set when top of boundary layer
C                                ! is reached.
      REAL
     + BQ(P_FIELD,2:BL_LEVELS)   ! A buoyancy parameter (beta q tilde).
     +,BT(P_FIELD,2:BL_LEVELS)   ! A buoyancy parameter (beta T tilde).
C                                ! Re-used for DZLUV(,2:BL_LEVELS).
     &,BF(P_FIELD,2:BL_LEVELS)   
!        A buoyancy parameter (beta F tilde).
     +,RI(P_FIELD,2:BL_LEVELS)   ! Richardson number for lower interface
C                                ! of layer.
     +,RIM(P_FIELD,2:BL_LEVELS)  ! Modified Richardson number for lower
C                                ! interface of layer.
     +,DELTAP(P_FIELD,BL_LEVELS) ! Pressure difference across layer.
     +,DTRDZ(P_FIELD,BL_LEVELS)  ! -g.dt/dp for model layers.
     +,DELTAP_RML(P_FIELD)       ! Pressure difference across the
C                                ! rapidly mixing layer (if it exists).
     +,QSL(P_FIELD)              ! Saturated sp humidity at pressure and
C                                ! liquid/frozen water temperature of
C                                ! successive levels.
C                                ! Re-used for RHOKM.
      INTEGER
     & NTML(P_FIELD)             ! No. of turbulently mixed model levels
C
C
C  (b) Scalars.
C
      REAL
     + AL      ! Temporary in calculation of buoyancy parameters.
     +,ALPHAL  ! Temporary in calculation of buoyancy parameters.
     +,BETAC   ! Temporary in calculation of buoyancy parameters.
     +,BETAQ   ! Temporary in calculation of buoyancy parameters.
     +,BETAT   ! Temporary in calculation of buoyancy parameters.
     +,BQM     ! Temporary in calculation of Richardson number RI(,K)
     +,BTM     ! Temporary in calculation of Richardson number RI(,K)
     &,BFM     
!        Temporary in calculation of Richardson number RI(,K)
     +,DZB     ! Temporary in calculation of Richardson number RI(,K).
     +,DZQW    ! Difference of QW between "current" and "lower" levels.
     &,DZQI    
!        Difference of QI between "current" and "lower" levels.
     +,DZTL    ! Liquid/ice static energy difference between adjacent
C              ! levels.
     +,DZU     ! Westerly wind shear between adjacent levels.
     +,DZV     ! Southerly wind shear between adjacent levels.
     +,DVMOD2  ! Square of modulus of wind shear between adjacent levels
     +,WK      ! Temporary in calculation of RHO.
     +,WKM1    ! Temporary in calculation of RHO.
     +,DTRDZ_RML   ! -g.dt/dp for the rapidly mixing layer.
     +,DTL_RML     ! Explicit TL increment for the rapidly mixing layer.
     +,DQW_RML     ! Explicit QW increment for the rapidly mixing layer.
     +,DTL_RMLP1   ! Explicit TL increment for the model layer above
C                  ! the rapidly mixing layer.
     +,DQW_RMLP1   ! Explicit QW increment for the model layer above
C                  ! the rapidly mixing layer.
     +,RIT         ! Temporary calculation of the modified Richardson
C                  ! number at the interface between the
C                  ! rapidly mixing layer and the model layer above it.
      INTEGER
     + I       ! Loop counter (horizontal field index).
     +,J       ! Offset counter in certain I loops.
     +,K       ! Loop counter (vertical level index).
     +,KM1     ! K-1.
     +,KP1     ! K-1
     +,MBL     ! Maximum number of model layers allowed in the rapidly
C              ! mixing layer; set to BL_LEVELS-1.
     +,NRMLP1  ! NRML + 1
     +,NRMLP2  ! NRML + 2
     +,IT_COUNTER  ! Iteration loop counter.
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
        CALL TIMER('KMKH    ',3)
      ENDIF
      ERROR=0
      IF (                                                          
     +     U_POINTS .NE. (U_ROWS*ROW_LENGTH) .OR.
     +     P_POINTS .NE. (P_ROWS*ROW_LENGTH) )  THEN
        ERROR=1
        GOTO9
      ENDIF
C
C  Set MBL, "maximum number of boundary levels" for the purposes of
C  boundary layer height calculation.
C
      MBL = BL_LEVELS - 1
C
C-----------------------------------------------------------------------
CL 1.  First loop round "boundary" levels.
C-----------------------------------------------------------------------
C
      DO 1 K=2,BL_LEVELS
C
C-----------------------------------------------------------------------
CL 1.1 Calculate saturated specific humidity at pressure and ice/liquid
CL     temperature of current level (TL, which is calculated here).
C      Store pressure temporarily in BQ(*,K).
C-----------------------------------------------------------------------
C
        DO 11 I=P1,P1+P_POINTS-1
          BQ(I,K) = AK(K) + BK(K)*PSTAR(I)
          IF (L_BL_LSPICE) THEN
            TL(I,K) = T(I,K) - LCRCP*QCL(I,K)                   ! P243.9
          ELSE
            TL(I,K) = T(I,K) - LCRCP*QCL(I,K) - LSRCP*QCF(I,K)  ! P243.9
          ENDIF
   11   CONTINUE
        IF (L_BL_LSPICE) THEN
          CALL QSAT_WAT(QSL(P1),TL(P1,K),BQ(P1,K),P_POINTS)
        ELSE
          CALL QSAT(QSL(P1),TL(P1,K),BQ(P1,K),P_POINTS)
        ENDIF
C
C-----------------------------------------------------------------------
CL 1.2 Calculate buoyancy parameters BT and BQ, required for the
CL     calculation of Richardson numbers above the surface.
C-----------------------------------------------------------------------
C
        DO 12 I=P1,P1+P_POINTS-1
C
CL 1.2.1 Calculate total water content (QW) at current level.
C
          IF (L_BL_LSPICE) THEN
            QW(I,K) = Q(I,K) + QCL(I,K)                        ! P243.10
          ELSE
            QW(I,K) = Q(I,K) + QCL(I,K) +QCF(I,K)              ! P243.10
          ENDIF
          BETAT = 1.0/T(I,K)                         ! P243.19 (1st eqn)
          BETAQ = C_VIRTUAL/(1.0 +C_VIRTUAL*Q(I,K) -QCL(I,K) -QCF(I,K))
C                                                  ... P243.19 (2nd eqn)
C
          IF (TL(I,K).GT.TM.OR.L_BL_LSPICE) THEN
            ALPHAL = (EPSILON * LC * QSL(I)) / ( R * TL(I,K) * TL(I,K) )
C                 ... P243.20 (Clausius-Clapeyron) for TL above freezing
C
            AL = 1.0 / (1.0 + LCRCP*ALPHAL)
C                                      ... P243.21 for TL above freezing
C
            BETAC = CF(I,K) * AL * ( LCRCP*BETAT - ETAR*BETAQ )
C                            ... P243.19 (3rd eqn) for TL above freezing
C
          ELSE
            ALPHAL = (EPSILON * LS * QSL(I)) / ( R * TL(I,K) * TL(I,K) )
C                 ... P243.20 (Clausius-Clapeyron) for TL below freezing
C
            AL = 1.0 / (1.0 + LSRCP*ALPHAL)
C                                      ... P243.21 for TL below freezing
C
            BETAC = CF(I,K) * AL * ( LSRCP*BETAT - ETAR*BETAQ )
C                            ... P243.19 (3rd eqn) for TL below freezing
C
          ENDIF
          BT(I,K) = BETAT - ALPHAL*BETAC             ! P243.18 (1st eqn)
          BQ(I,K) = BETAQ + BETAC                    ! P243.18 (2nd eqn)
          BF(I,K) = BETAQ*EPSILON*ETAR               ! P243.18 (2nd eqn)
!
          IF (L_BL_LSPICE) THEN
            TL(I,K) = T(I,K) - LCRCP*QCL(I,K)                   ! P243.9
          ENDIF
!
   12   CONTINUE
    1 CONTINUE
C
C-----------------------------------------------------------------------
CL 1.3 Initialise flag for having reached top of boundary layer
CL     and also the number of turbulently mixed layers.
C-----------------------------------------------------------------------
C
      DO 13 I=P1,P1+P_POINTS-1
        TOPBL(I) = .FALSE.
        NTML(I) = 1
   13 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.  Second loop round "boundary" levels.
C-----------------------------------------------------------------------
C
      DO 2 K=2,BL_LEVELS
        KM1 = K-1
        DO 21 I=P1,P1+P_POINTS-1
C
C-----------------------------------------------------------------------
CL 2.1 Calculate wind shear and other "jumps" between "current" and
CL     previous (lower) level.
C-----------------------------------------------------------------------
C
          DZU = U_P(I,K) - U_P(I,KM1)
          DZV = V_P(I,K) - V_P(I,KM1)
          DVMOD2 = MAX ( 1.0E-12 , DZU*DZU + DZV*DZV ) ! Used in P243.C1
          DZTL = TL(I,K) - TL(I,KM1) + GRCP/RDZ(I,K)   ! Used in P243.C2
          DZQW = QW(I,K) - QW(I,KM1)                   ! Used in P243.C2
          DZQI = QCF(I,K) - QCF(I,KM1)                 ! Used in P243.C2
C
C-----------------------------------------------------------------------
CL 2.2 Calculate buoyancy parameters BT, BQ, DZB at interface between
CL     "current" and previous levels (i.e. at level K-1/2, if current
CL     level is level K).
C-----------------------------------------------------------------------
C
          WKM1 = 0.5 * DZL(I,KM1) * RDZ(I,K)         ! P243.C5 (2nd eqn)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)             ! P243.C5 (1st eqn)
          IF (K.EQ.2) THEN
            BTM = WK*BT(I,K) + WKM1*BT_1(I)          ! P243.C3 (1st eqn)
            BQM = WK*BQ(I,K) + WKM1*BQ_1(I)          ! P243.C3 (2nd eqn)
            BFM = WK*BF(I,K) + WKM1*BF_1(I)          ! P243.C3 (2nd eqn)
          ELSE
            BTM = WK*BT(I,K) + WKM1*BT(I,KM1)        ! P243.C3 (1st eqn)
            BQM = WK*BQ(I,K) + WKM1*BQ(I,KM1)        ! P243.C3 (2nd eqn)
            BFM = WK*BF(I,K) + WKM1*BF(I,KM1)        ! P243.C3 (2nd eqn)
          ENDIF
          IF (L_BL_LSPICE) THEN
            DZB = BTM*DZTL + BQM*DZQW - BFM*DZQI       ! P243.C2 (no g)
          ELSE
            DZB = BTM*DZTL + BQM*DZQW                  ! P243.C2 (no g)
          ENDIF
C
C  For rationale of next IF block, see the discussion in the last
C  paragraph of Appendix C of the P243 documentation.
C
          IF (DZB.GT.0.0) THEN
            IF (K.EQ.2) THEN
              BTM = 0.5 * ( BT(I,K) + BT_1(I) )
              BQM = 0.5 * ( BQ(I,K) + BQ_1(I) )
              BFM = 0.5 * ( BF(I,K) + BF_1(I) )
            ELSE
              BTM = 0.5 * ( BT(I,K) + BT(I,KM1) )
              BQM = 0.5 * ( BQ(I,K) + BQ(I,KM1) )
              BFM = 0.5 * ( BF(I,K) + BF(I,KM1) )
            ENDIF
          IF (L_BL_LSPICE) THEN
            DZB = BTM*DZTL + BQM*DZQW - BFM*DZQI
          ELSE
            DZB = BTM*DZTL + BQM*DZQW
          ENDIF
          ENDIF
          IF (K.EQ.2) THEN
            BT_1(I) = BTM
            BQ_1(I) = BQM
            BF_1(I) = BFM
          ELSE
            BT(I,KM1) = BTM
            BQ(I,KM1) = BQM
            BF(I,KM1) = BFM
          ENDIF
C
C-----------------------------------------------------------------------
CL 2.3 Calculate Richardson number Ri at the same interface.
C-----------------------------------------------------------------------
C
          RI(I,K) = G * DZB / ( RDZ(I,K) * DVMOD2 )
C
C-----------------------------------------------------------------------
CL 2.4 If either a stable layer (Ri > 1) or the maximum boundary layer
CL     height has been reached, set boundary layer height (ZH) to
CL     the height of the lower boundary of the current layer and set
CL     the number of rapidly mixing layers if the surface layer is
CL     unstable (as determined in subroutine SF_EXCH).
C-----------------------------------------------------------------------
C
          IF ( .NOT.TOPBL(I) .AND. (RI(I,K).GT.1.0 .OR. K.GT.MBL) ) THEN
            TOPBL(I) = .TRUE.
            ZH(I) = ZLB(I,K)
            NTML(I) = K-1
            IF ( NRML(I) .GT. 0 ) NRML(I) = K-1
          ENDIF
   21   CONTINUE
    2 CONTINUE
C
C-----------------------------------------------------------------------
CL 3.  Subroutine EX_COEF.
C-----------------------------------------------------------------------
C
      CALL EX_COEF (
     + P_FIELD,U_FIELD,P1,P_POINTS,BL_LEVELS,
     + CCB,CCT,NTML,L_MOM,L_MIXLEN,
     + AKH,BKH,CCA,DZL,PSTAR,RDZ,RI,TV,U_P,V_P,ZH,ZLB,Z0M,H_BLEND,
     + RHOKM,RHOKH,LTIMER
     +)
C
C-----------------------------------------------------------------------
CL 4.  Calculate "explicit" fluxes of TL and QW.
C-----------------------------------------------------------------------
C
      DO 4 K=2,BL_LEVELS
        KM1 = K-1
C-----------------------------------------------------------------------
CL 4.1 "Explicit" fluxes of TL and QW, on P-grid.
C-----------------------------------------------------------------------
        DO 41 I=P1,P1+P_POINTS-1
          FTL(I,K) = -RHOKH(I,K) *
     +      ( ( ( TL(I,K) - TL(I,KM1) ) * RDZ(I,K) ) + GRCP )
C           1 2 3                     3            2        1
          FQW(I,K) = -RHOKH(I,K) * ( QW(I,K) - QW(I,KM1) ) * RDZ(I,K)
   41   CONTINUE
    4 CONTINUE
C
C-----------------------------------------------------------------------
CL 4.2 Use explicit fluxes to calculate a modified Richardson number
CL     at the top of the rapidly mixing layer (if it exists); if this
CL     indicates that the r.m.l. can deepen due to heat and/or moisture
CL     input from the surface then increase NRML(I).
C-----------------------------------------------------------------------
C
C Initialise pressure difference, dp, for the rapidly mixing layer.
C
      DO I = P1,P1+P_POINTS-1
        DELTAP_RML(I) = 0.0
      ENDDO ! Loop over p-points
C
C Calculate pressure differences, dp, and -g.dt/dp for model layers
C and dp for the rapidly mixing layer.
C
      DO 420 K = 1,BL_LEVELS
        DO I = P1,P1+P_POINTS-1
          DELTAP(I,K) = DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)
          DTRDZ(I,K) = -G * TIMESTEP / DELTAP(I,K)
          IF (K .LE. NRML(I)) DELTAP_RML(I)
     &                           = DELTAP_RML(I) + DELTAP(I,K)
          IF (K .GE. 2) RIM(I,K) = RI(I,K)
        ENDDO ! Loop over p-points
 420  CONTINUE ! Loop over bl-levels
C
C-----------------------------------------------------------------------
CL 4.2.1 Iterate BL_LEVELS-2 times; this allows an initial r.m.l. of
CL       1 model layer to deepen to BL_LEVELS-1 model layers a layer at
CL       a time in each iteration.  Some of these iterations may be
CL       redundant if in fact the r.m.l. cannot deepen.
C-----------------------------------------------------------------------
      DO 421 IT_COUNTER = 1,BL_LEVELS-2
C
        DO I = P1,P1+P_POINTS-1
C         !-------------------------------------------------------------
C         ! Only check to see if the r.m.l. can deepen if it exists
C         ! in the first place and has less than its maximum depth.
C         !-------------------------------------------------------------
          IF ((NRML(I) .GE. 1) .AND. (NRML(I) .LE. BL_LEVELS-2)) THEN
            NRMLP1 = NRML(I) + 1
            NRMLP2 = NRML(I) + 2
            DTRDZ_RML = -G * TIMESTEP / DELTAP_RML(I)
C           !-----------------------------------------------------------
C           ! "Explicit" rapidly mixing layer increments to TL and QW.
C           !-----------------------------------------------------------
            DTL_RML = -DTRDZ_RML * ( FTL(I,NRMLP1) - FTL(I,1) )
            DQW_RML = -DTRDZ_RML * ( FQW(I,NRMLP1) - FQW(I,1) )
C           !-----------------------------------------------------------
C           ! "Explicit" increments to TL and QW in the model layer
C           ! above the rapidly mixing layer.
C           !-----------------------------------------------------------
            DTL_RMLP1 = -DTRDZ(I,NRMLP1) *
     &                           ( FTL(I,NRMLP2) - FTL(I,NRMLP1) )
            DQW_RMLP1 = -DTRDZ(I,NRMLP1) *
     &                           ( FQW(I,NRMLP2) - FQW(I,NRMLP1) )
            DZU = U_P(I,NRMLP1) - U_P(I,NRML(I))
            DZV = V_P(I,NRMLP1) - V_P(I,NRML(I))
            DVMOD2 = MAX (1.0E-12 , DZU*DZU + DZV*DZV)
C           !-----------------------------------------------------------
C           ! Calculate a modified Richardson number for the interface
C           ! between the rapidly mixing layer and the model layer above
C           !-----------------------------------------------------------
            IF (NRML(I) .EQ. 1) THEN
              RIT = RI(I,NRMLP1) + ( G/(RDZ(I,NRMLP1) * DVMOD2) ) *
     &                     ( BT_1(I) * ( DTL_RMLP1 - DTL_RML )
     &                     + BQ_1(I) * ( DQW_RMLP1 - DQW_RML ) )
            ELSE ! NRML(I) .GT. 1
              RIT = RI(I,NRMLP1) + ( G/(RDZ(I,NRMLP1) * DVMOD2) ) *
     &                ( BT(I,NRML(I)) * ( DTL_RMLP1 - DTL_RML )
     &                + BQ(I,NRML(I)) * ( DQW_RMLP1 - DQW_RML ) )
            ENDIF
            IF ( RIT .LE. 1.0 ) THEN
C             !---------------------------------------------------------
C             ! Deepen the rapidly mixing layer by one model layer
C             !---------------------------------------------------------
              NRML(I) = NRMLP1
              DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,NRMLP1)
              ZH(I) = ZLB(I,NRMLP2)
              IF (RIT .LT. RI(I,NRMLP1)) RIM(I,NRMLP1) = RIT
            ENDIF
          ENDIF ! NRML(I) .GE. 1 and .LE. BL_LEVELS-2
        ENDDO ! Loop over p-points
 421  CONTINUE! Loop over iterations
C
C-----------------------------------------------------------------------
CL 4.2.2 Adjust the Richardson number at the top of the rapidly mixing
CL       layer- the 'DZ' in the expression (P243.C1) is set to 100.0 m
CL       rather than the model level separation (if the latter is
CL       greater than 100 m).  This is a simple way of adjusting for
CL       inaccuracies in calculating gradients at inversions when the
CL       model's vertical resolution is coarse.
C-----------------------------------------------------------------------
C
      DO 422 K = 2,BL_LEVELS
        KM1 = K-1
        DO I = P1,P1+P_POINTS-1
          IF ( (KM1 .EQ. NRML(I)) .AND. (100.0*RDZ(I,K) .LT. 1.0) )
     &       RIM(I,K) = 100.0*RDZ(I,K)*RIM(I,K)
        ENDDO ! Loop over p-points
 422  CONTINUE! Loop over bl-levels
C
C-----------------------------------------------------------------------
CL 5.  Subroutine EX_COEF using adjusted Richardson Number, RIM
C-----------------------------------------------------------------------
C
      CALL EX_COEF (
     + P_FIELD,U_FIELD,P1,P_POINTS,BL_LEVELS,
     + CCB,CCT,NTML,L_MOM,L_MIXLEN,
     + AKH,BKH,CCA,DZL,PSTAR,RDZ,RIM,TV,U_P,V_P,ZH,ZLB,Z0M,H_BLEND,
     + RHOKM,RHOKH,LTIMER
     +)
C
C-----------------------------------------------------------------------
CL 6.  Recalculate "explicit" fluxes of TL and QW, see section 4. above.
C-----------------------------------------------------------------------
C
      DO 6 K=2,BL_LEVELS
        KM1 = K-1
C-----------------------------------------------------------------------
CL 6.1 "Explicit" fluxes of TL and QW, on P-grid.
CL       Overwrite exchange coefficient with GAMMA*RHOKH*RDZ to avoid
CL       unecessary multiplication in P244 (IMPL_CAL). RHOKH only used
CL       in this combination hearafter.
CL      RHOKM stored here for diagnostic output (N.B. before interpoln.)
C-----------------------------------------------------------------------
        DO 61 I=P1,P1+P_POINTS-1
          FTL(I,K) = -RHOKH(I,K) *
     +      ( ( ( TL(I,K) - TL(I,KM1) ) * RDZ(I,K) ) + GRCP )
C           1 2 3                     3            2        1
          FQW(I,K) = -RHOKH(I,K) * ( QW(I,K) - QW(I,KM1) ) * RDZ(I,K)
          RHOKH(I,K) = GAMMA(K) * RHOKH(I,K) * RDZ(I,K)
          RHO_KM(I,K) = RHOKM(I,K)
   61   CONTINUE
    6 CONTINUE
C
C
C-----------------------------------------------------------------------
CL 7. Calculate "explicit" turbulent stress components.
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
CL 7.1 "Explicit" turbulent stress components, on UV-grid.
CL       Overwrite exchange coefficient with GAMMA*RHOKM*RDZUV to avoid
CL       unecessary multiplication in P244 (IMPL_CAL). RHOKM only used
CL       in this combination hearafter.
C-----------------------------------------------------------------------
C
C  Store DZL(,1) on UV-grid, in BT_1.
C
      CALL P_TO_UV(DZL(P1,1),BT_1(U1+ROW_LENGTH),P_POINTS,
     +             P_POINTS,ROW_LENGTH,P_ROWS)
C
C  Need to update haloes of RHOKM before calling P_TO_UV
      CALL SWAPBOUNDS(RHOKM(1,2),ROW_LENGTH,
     & U_FIELD/ROW_LENGTH,1,1,BL_LEVELS-1)

      DO 7 K = 2,BL_LEVELS
        KM1 = K-1
C
C  Store DZL(,K) on UV-grid in BT(,K), interpolate RHOKM to UV-grid via
C  QSL.  Note addressing of array elements; the first and last rows
C  of genuinely UV-grid arrays are never accessed here.  They get set
C  to "missing data" (1.0E30) before being passed out of the routine.
C
        CALL P_TO_UV (DZL(P1,K),BT(U1+ROW_LENGTH,K),
     +     P_POINTS,P_POINTS,ROW_LENGTH,P_ROWS)
C
        CALL P_TO_UV (RHOKM(P1,K),QSL(U1+ROW_LENGTH),
     +     P_POINTS,P_POINTS,ROW_LENGTH,P_ROWS)
C
        DO 71 I=U1+ROW_LENGTH,U1-ROW_LENGTH+U_POINTS-1
          RHOKM(I,K) = QSL(I)
          IF (K.EQ.2) THEN
            RDZUV(I,K) = 2.0 / ( BT(I,K) + BT_1(I) )
          ELSE
            RDZUV(I,K) = 2.0 / ( BT(I,K) + BT(I,KM1) )
          ENDIF
          TAUX(I,K) = RHOKM(I,K) * ( U(I,K) - U(I,KM1) ) * RDZUV(I,K)
          TAUY(I,K) = RHOKM(I,K) * ( V(I,K) - V(I,KM1) ) * RDZUV(I,K)
          RHOKM(I,K) = GAMMA(K) * RHOKM(I,K) * RDZUV(I,K)
   71   CONTINUE
C
C-----------------------------------------------------------------------
CL 8.  Set first and last rows of UV-grid variables to "missing data".
C-----------------------------------------------------------------------
C
      IF (attop) THEN
        DO I=U1,U1+ROW_LENGTH-1
          RHOKM(I,K) = 1.0E30
          TAUX(I,K)  = 1.0E30
          TAUY(I,K)  = 1.0E30
          RDZUV(I,K) = 1.0E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= U1 + (U_ROWS-1)*ROW_LENGTH , U1 + U_ROWS*ROW_LENGTH - 1
          RHOKM(I,K) = 1.0E30
          TAUX(I,K)  = 1.0E30
          TAUY(I,K)  = 1.0E30
          RDZUV(I,K) = 1.0E30
        ENDDO
      ENDIF
    7 CONTINUE
C
C
    9 CONTINUE  ! Branch for error exit.
C
      IF (LTIMER) THEN
        CALL TIMER('KMKH    ',4)
      ENDIF
C
      RETURN
      END
