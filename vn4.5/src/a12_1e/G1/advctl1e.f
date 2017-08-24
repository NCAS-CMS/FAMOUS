C ******************************COPYRIGHT******************************
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
C
CLL   SUBROUTINE ADV_CTL --------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES THE RIGHT-HAND SIDES OF EQUATIONS (40) TO
CLL              (42) REPRESENTING THE MASS WEIGHTED FIELDS AFTER
CLL              ADVECTION AND THE ADDITION OF THE CORIOLIS TERM DUE
CLL              TO VERTICAL MOTION. THE SPATIAL DIFFERENCING SCHEME
CLL              (35) TO (38) IS USED. ONE MORE PRESSURE ROW THAN
CLL              VELOCITY ROW IS UPDATED. DIVERGENCE DAMPS VELOCITY
CLL              FIELDS AS DESCRIBED IN SECTION 3.4 OF DOCUMENTATION
CLL              PAPER NO. 10
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   WAS VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Required for new interface to THQADV
!LL                   Version 1E not bit reprod with 1C
CLL   4.4    04/08/97  Optimisation for T3E  D.Salmond
CLL
CLL
CLL
CLL   PROGRAMMING STANDARD:
CLL
CLL   LOGIACL COMPONENTS COVERED: P12
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION:        THE EQUATIONS USED ARE (35) TO (46) AND
CLL                         SECTION 3.4 IN UNIFIED MODEL DOCUMENTATION
CLL                         NO. 10  M.J.P. CULLEN, T.DAVIES AND
CLL                         M.H. MAWSON VERSION 17, DATED 11/02/91.
CLLEND-------------------------------------------------------------
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE ADV_CTL
     1              (THETAL,QT,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,U,V,
     &              COS_U_LATITUDE,COS_P_LATITUDE,
     2              SEC_P_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,
     3              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     4              LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5              SOUTHERN_FILTERED_P_ROW,Q_LEVELS,
     6              U_FIELD,P_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     7              P_LEVELS,SEC_U_LATITUDE,F1,F2,AK,BK,KD,AKH,BKH,
     8              COS_U_LONGITUDE,SIN_U_LONGITUDE,TRIGS,IFAX,
     9              FILTER_WAVE_NUMBER_P_ROWS,OMEGA,QCL,QCF,P_EXNER,
     &              LLINTS,LWHITBROM,
     &              L_TRACER_THETAL_QT,NSWEEP,L_SUPERBEE)

      IMPLICIT NONE

! All TYPFLDPT arguments are intent IN
! Comdeck TYPFLDPT
! Variables which point to useful positions in a horizontal field

      INTEGER
     &  FIRST_ROW        ! First updatable row on field
     &, TOP_ROW_START    ! First point of north-pole (global) or
!                        ! Northern (LAM) row
!                        ! for processors not at top of LPG, this
!                        ! is the first point of valid data
!                        ! (ie. Northern halo).
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
!                        ! for processors not at base of LPG, this
!                        ! is the start of the last row of valid data
!                        ! (ie. Southern halo).
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
! For the last three variables, these indexes are the start points
! and end points of "local" data - ie. missing the top and bottom
! halo regions.
     &, FIRST_VALID_PT   ! first valid point of data on field
     &, LAST_P_VALID_PT  ! last valid point of data on field
     &, LAST_U_VALID_PT  ! last valid point of data on field
     &, VALID_P_ROWS     ! number of valid rows of P data
     &, VALID_U_ROWS     ! number of valid rows of U data
     &, START_POINT_NO_HALO
!                        ! first non-polar point of field (misses
!                        ! halo for MPP code)
     &, START_POINT_INC_HALO
!                        ! first non-polar point of field (includes
!                        ! halo for MPP code)
     &, END_P_POINT_NO_HALO
!                        ! last non-polar point of P field (misses
!                        ! halo for MPP code)
     &, END_P_POINT_INC_HALO
!                        ! last non-polar point of P field (includes
!                        ! halo for MPP code)
     &, END_U_POINT_NO_HALO
!                        ! last non-polar point of U field (misses
!                        ! halo for MPP code)
     &, END_U_POINT_INC_HALO
!                        ! last non-polar point of U field (includes
!                        ! halo for MPP code)
     &, FIRST_ROW_PT     ! first data point along a row
     &, LAST_ROW_PT      ! last data point along a row
! For the last two variables, these indexes are the start and
! end points along a row of the "local" data - ie. missing out
! the east and west halos
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!

     &, MY_PROC_ID         ! my processor id
     &, NP_PROC_ID         ! processor number of North Pole Processor
     &, SP_PROC_ID         ! processor number of South Pole Processor
     &, GC_ALL_GROUP       ! group id of group of all processors
     &, GC_ROW_GROUP       ! group id of group of all processors on this
!                          ! processor row
     &, GC_COL_GROUP       ! group id of group of all processors on this
!                          ! processor column
     &, N_PROCS            ! total number of processors

     &, EW_Halo            ! Halo size in the EW direction
     &, NS_Halo            ! Halo size in the NS direction

     &, halo_4th           ! halo size for 4th order calculations
     &, extra_EW_Halo      ! extra halo size required for 4th order
     &, extra_NS_Halo      ! extra halo size required for 4th order
     &, LOCAL_ROW_LENGTH   ! size of local row
     &, FIRST_GLOBAL_ROW_NUMBER
!                          ! First row number on Global Grid    

! Variables which indicate if special operations are required at the
! edges.
      LOGICAL
     &  at_top_of_LPG    ! Logical variables indicating if this
     &, at_right_of_LPG  ! processor is at the edge of the Logical
     &, at_base_of_LPG   ! Processor Grid and should process its edge
     &, at_left_of_LPG   ! data differently.

! End of comdeck TYPFLDPT
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

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN    NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *                     ! MOVING TOWARDS THE EQUATOR.
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN.
     *                     ! MOVING TOWARDS SOUTHPOLE.
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.
     *, NSWEEP(glsize(2),P_LEVELS) !IN No.of EW sweeps for all rows.
     *                  ! NUMBER OF EAST_WEST TIMESTEPS NEEDED FOR
     *                  ! EACH LATITUDE WHEN USING TRACER ADVECTION.
     *, FIRST_POINT     !
     *, POINTS          !

      INTEGER
     &  FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
!       LAST WAVE NUMBER NOT TO BE CHOPPED ON A P ROW
      LOGICAL
     &  L_SUPERBEE             ! FORM OF LIMITER USED IN TRACER
     &                         ! ADVECTION
     & ,L_TRACER_THETAL_QT     ! LOGICAL TRUE IF USING TRACER
     &                         ! ADVECTION FOR THETAL & QT
       INTEGER
     &  P_POINTS_UPDATE
     & ,START_U_REQUIRED
     & ,P_POINTS_REQUIRED
     & ,U_POINTS_REQUIRED

     & ,LLINTS              !Logical switch for linear TS
     & ,LWHITBROM           !Log swch for White & Bromley terms

      REAL
     * U_MEAN(U_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                           !   FROM ADJUSTMENT STEP
     *,V_MEAN(U_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                           !   * COS(LAT) FROM ADJUSTMENT STEP
     *,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *,PSTAR(P_FIELD)            !IN PSTAR FIELD AT NEW TIME-LEVEL
     *,PSTAR_OLD(P_FIELD)        !IN PSTAR AT PREVIOUS TIME-LEVEL
     *,RS(P_FIELD,P_LEVELS)      !IN RS FIELD
     *,TRIGS(ROW_LENGTH)        !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                          ! IN FILTERING.
     *,QCL(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCL.
     *,QCF(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCF.
     *,P_EXNER(P_FIELD,P_LEVELS+1) !IN. PRIMARY ARRAY FOR P_EXNER.

      REAL
     * U(U_FIELD,P_LEVELS)       !INOUT U FIELD, MASS-WEIGHTED ON OUT.
     *,V(U_FIELD,P_LEVELS)       !INOUT V FIELD, MASS-WEIGHTED ON OUT.
     *,THETAL(P_FIELD,P_LEVELS)  !INOUT THETAL FIELD
     *,QT(P_FIELD,Q_LEVELS)      !INOUT QT FIELD.

      REAL
     * DELTA_AK(P_LEVELS)     !IN    LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)     !IN    LAYER THICKNESS
     *,AK(P_LEVELS)           !IN    FIRST TERM IN HYBRID CO-ORDS.
     *,BK(P_LEVELS)           !IN    SECOND TERM IN HYBRID CO-ORDS.
     *,AKH(P_LEVELS+1)        !IN    AK AT HALF LEVELS
     *,BKH(P_LEVELS+1)        !IN    BK AT HALF LEVELS
     &,COS_P_LATITUDE(P_FIELD) !IN  COS_LAT AT P_POINTS (2D ARRAY)
     *,SEC_P_LATITUDE(P_FIELD) !IN  1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,COS_U_LATITUDE(U_FIELD) !IN  COS(LAT) AT U POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN  1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,COS_U_LONGITUDE(ROW_LENGTH) !IN COS(LONGITUDE) AT U-POINTS.
     *,SIN_U_LONGITUDE(ROW_LENGTH) !IN SIN(LONGITUDE) AT U-POINTS.
     *,LONGITUDE_STEP_INVERSE !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE  !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP     !IN
     *,NU_BASIC               !IN STANDARD NU TERM FOR MODEL RUN.
     *,F1(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F2(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,KD(P_LEVELS)           !IN DIVERGENCE DAMPING COEFFICIENTS

      REAL
     * OMEGA(U_FIELD,P_LEVELS) !OUT TRUE VERTICAL VELOCITY
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C  DEFINE LOCAL ARRAYS: 3 ARE REQUIRED
      REAL
     * WORK1(U_FIELD)         ! GENERAL WORKSPACE
     *,WORK2(P_FIELD)         ! GENERAL WORKSPACE
     &,   OMEGA_P(P_FIELD)    ! HOLDS OMEGA AT P POINTS.

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,K,K1

      INTEGER X_FIELD  ! 1 IF 2ND ORDER ELSE U_FIELD IF 4TH ORDER

C   REAL SCALARS
      REAL
     &  CONST1,LC_LF,TIMESTEP
     &  ,PK, PK1       ! Pressure at half levels
     &  ,P_EXNER_FULL  ! Exner pressure at full model level

C LOGICAL VARIABLE
      LOGICAL
     * L_SECOND                ! SET TO TRUE IF NU_BASIC EQUAL TO ZERO
C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL TH_Q_ADV,UV_ADV,P_TO_UV,DIV_DAMP
     &         ,TRAC_ADV,TRAC_VERT_ADV,UV_TO_P,POLAR
C*---------------------------------------------------------------------

      INTEGER extended_address(P_FIELD)

CL    COMDECK C_THADV HOLDS PHYSICAL CONSTANTS REQUIRED BY ROUTINE
CL    TH_ADV.
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
CL    END OF COMDECK C_THADV.
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------



      IF (NU_BASIC .NE. 0.0) THEN
! Calculate the mapping between points on the normal horizontal
! field, and points in the extended field (with double halos for
! the fourth order code)
! Logic: extended_address=old_address
!         + ROW_LENGTH*extra_NS_Halo
!           -> extra halo row at top of field
!         + (row_number+1)*2*extra_EW_Halo
!           -> two extra halo points for each preceeding row
!         + extra_EW_Halo   -> extra halo point at start of this row
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          extended_address(I)=I + ROW_LENGTH*extra_NS_Halo +
     &      (((I-1)/ROW_LENGTH)+extra_NS_Halo)*2*extra_EW_Halo +
     &      extra_EW_Halo
        ENDDO
      ENDIF
CL    MAXIMUM VECTOR LENGTH ASSUMED IS U_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
C****************************************************************
C         INTEGERS AND VARIABLES NEEDED WHEN USING
C         TRACER ADVECTION OF THETAL & QT
C*****************************************************************
      IF(L_TRACER_THETAL_QT)THEN
       LC_LF=LC + LF
       P_POINTS_UPDATE=upd_P_ROWS*ROW_LENGTH
       START_U_REQUIRED  = START_POINT_NO_HALO-ROW_LENGTH
       P_POINTS_REQUIRED = (upd_P_ROWS+2)*ROW_LENGTH
       U_POINTS_REQUIRED = (upd_U_ROWS+2)*ROW_LENGTH
      ENDIF
CL---------------------------------------------------------------------
CL    SECTION 1.     INTERPOLATE FIELDS ONTO U GRID.
CL---------------------------------------------------------------------

      IF(NU_BASIC.EQ.0.) THEN
        L_SECOND=.TRUE.
      X_FIELD=1
      ELSE
        L_SECOND=.FALSE.
      X_FIELD=U_FIELD
      END IF
! Initialise arrays WORK1 & WORK2
      DO I = 1,P_FIELD
        WORK1(I) = 1.0
        WORK2(I) = 1.0
      END DO

C----------------------------------------------------------------------
CL    SECTION 1.1    INTERPOLATE PSTAR ONTO U GRID.
C----------------------------------------------------------------------

      CALL P_TO_UV(PSTAR,WORK1,P_FIELD,U_FIELD,ROW_LENGTH,tot_P_ROWS)

C----------------------------------------------------------------------
CL    SECTION 1.2    INTERPOLATE PSTAR_OLD ONTO U GRID.
C----------------------------------------------------------------------

      CALL P_TO_UV(PSTAR_OLD,WORK2,P_FIELD,U_FIELD,ROW_LENGTH,
     &             tot_P_ROWS)

! Update the halos of WORK1 and WORK2
      CALL SWAPBOUNDS(WORK1,ROW_LENGTH,tot_P_ROWS,EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(WORK2,ROW_LENGTH,tot_P_ROWS,EW_Halo,NS_Halo,1)

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALL DIV_DAMP TO PERFORM DIVERGENCE DAMPING.
CL---------------------------------------------------------------------

C PSTAR_OLD ON U GRID IS HELD IN WORK2.

      CALL DIV_DAMP(U,V,RS,SEC_U_LATITUDE,WORK2,COS_U_LATITUDE,KD,
     *              LONGITUDE_STEP_INVERSE,LATITUDE_STEP_INVERSE,
     *              P_FIELD,U_FIELD,ROW_LENGTH,P_LEVELS,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *              BKH,ADVECTION_TIMESTEP,DELTA_AK,DELTA_BK,
     *              COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE)

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALL UV_ADV TO ADVECT U AND V.
CL---------------------------------------------------------------------

C PSTAR ON U GRID IS HELD IN WORK1.
C PSTAR_OLD ON U GRID IS HELD IN WORK2.

      CALL UV_ADV (U,V,WORK2,WORK1,U_MEAN,V_MEAN,
     *             SEC_U_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,AK,
     *             BK,F1,F2,LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,
     *             NU_BASIC,LONGITUDE_STEP_INVERSE,U_FIELD,P_FIELD,
     *             ROW_LENGTH,P_LEVELS,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *             COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE,
     &             AKH,BKH,OMEGA,L_SECOND,LLINTS,
     &            extended_address,
     &             LWHITBROM,X_FIELD)
! Update the halos for the OMEGA array
      CALL SWAPBOUNDS(OMEGA,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)

! U and V are not swapped here, but in ATM_DYN, after the call to
! MASS_UWT which spoils the halo.

CL
CL---------------------------------------------------------------------
CL    SECTION 4.     CALL TH_Q_ADV TO ADVECT THETAL AND
CL                   QT USING STANDARD HEUN ADVECTION.
CL                   IF USING TRACER ADVECTION FOR THETAL & QT
CL                   THEN CALL APPROPRIATE TRACER ROUTINES.
CL---------------------------------------------------------------------
      IF(.NOT.L_TRACER_THETAL_QT)THEN
CL---------------------------------------------------------------
C    SECTION 4.1  HEUN ADVVECTION SCHEME
C
CL----------------------------------------------------------------

      CALL TH_Q_ADV (THETAL,QT,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,
     *            SEC_P_LATITUDE,
     *            ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,LATITUDE_STEP_INVERSE
     *            ,ADVECTION_TIMESTEP,NU_BASIC,LONGITUDE_STEP_INVERSE,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW,
     *            P_LEVELS,U_FIELD,P_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *            TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,SEC_U_LATITUDE,
     *            AKH,BKH,QCL,QCF,P_EXNER,OMEGA,
     &            Q_LEVELS,AK,BK,L_SECOND,
     &            extended_address,
     &            LWHITBROM)

! Update the halos for the THETAL array
      CALL SWAPBOUNDS(THETAL,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)

! Update the halos for the QT array
      CALL SWAPBOUNDS(QT,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,Q_LEVELS)
      ELSE
CL---------------------------------------------------------------
C    SECTION 4.2  TRACER ADVECTION OF THETAL AND QT
C
CL----------------------------------------------------------------
      DO K=1,P_LEVELS
        CALL TRAC_ADV(THETAL(1,K),NSWEEP(1,K),U_MEAN(1,K),V_MEAN(1,K),
     &                U_FIELD,P_FIELD,ADVECTION_TIMESTEP,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     &                SEC_P_LATITUDE,COS_P_LATITUDE,RS(1,K),
     &                PSTAR_OLD,DELTA_AK(K),DELTA_BK(K),
     &                LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                L_SUPERBEE)
      END DO

C  Set temperature flux through lower boundary to zero
      DO I=1,P_FIELD
        WORK2(I)=0.
      END DO

      FIRST_POINT = START_POINT_NO_HALO
      POINTS = upd_P_ROWS * ROW_LENGTH

      TIMESTEP=ADVECTION_TIMESTEP
      CONST1=R/(CP*CP)*TIMESTEP
      CALL TRAC_VERT_ADV(THETAL,ETADOT_MEAN,PSTAR,P_FIELD,
     &                   TIMESTEP,1,P_LEVELS,FIRST_POINT,
     &                   POINTS,P_LEVELS,1,P_LEVELS,RS,AK,BK,DELTA_AK,
     &                   DELTA_BK,WORK2,L_TRACER_THETAL_QT,L_SUPERBEE)
C ---------------------------------------------------------------------
CL                   INTERPOLATE OMEGA TO P GRID AND CALCULATE
CL                   REMAINING TERM IN ADVECTION EQUATION.
CL                   CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD.
C ---------------------------------------------------------------------

          DO 110 K=1,P_LEVELS

          CALL UV_TO_P(OMEGA(START_U_REQUIRED,K),
     &                 OMEGA_P(START_POINT_NO_HALO),U_POINTS_REQUIRED,
     &                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

C TOTAL MASS-WEIGHTED HORIZONTAL AND VERTICAL INCREMENTS ARE CALCULATED
C SEPARATELY.

          IF(K.LT.Q_LEVELS+1) THEN
            DO  I = FIRST_POINT,FIRST_POINT+POINTS-1

              PK  = AKH(K+1)+ BKH(K+1)*PSTAR(I)
              PK1 = AKH(K) + BKH(K)*PSTAR(I)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I,K+1),P_EXNER(I,K),PK,PK1,KAPPA)

              WORK2(I) =
     &                   -(LC*QCL(I,K)+LC_LF*QCF(I,K))*CONST1*
     &                    OMEGA_P(I)/((AK(K)+BK(K)*PSTAR(I))
     &                    *(P_EXNER_FULL)*
     &        RS(I,K)*RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)))
              THETAL(I,K) =THETAL(I,K)+WORK2(I)
            END DO
          END IF

CL END LOOP OVER P_LEVELS+1
 110  CONTINUE

! Update the halos for the THETAL array
      CALL SWAPBOUNDS(THETAL,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)


      DO K=1,Q_LEVELS
        CALL TRAC_ADV(QT(1,K),NSWEEP(1,K),U_MEAN(1,K),V_MEAN(1,K),
     &                U_FIELD,P_FIELD,ADVECTION_TIMESTEP,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     &                SEC_P_LATITUDE,COS_P_LATITUDE,RS(1,K),
     &                PSTAR_OLD,DELTA_AK(K),DELTA_BK(K),
     &                LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                L_SUPERBEE)
      END DO

C  Set moisture flux through lower boundary to zero
      DO I=1,P_FIELD
        WORK2(I)=0.
      END DO

! Values of FIRST_POINT and POINTS
! should be unaltered from those set for Thetal

      CALL TRAC_VERT_ADV(QT,ETADOT_MEAN,PSTAR,P_FIELD,
     &                   TIMESTEP,1,Q_LEVELS,FIRST_POINT,
     &                   POINTS,P_LEVELS,1,Q_LEVELS,RS,AK,BK,DELTA_AK,
     &                   DELTA_BK,WORK2,L_TRACER_THETAL_QT,L_SUPERBEE)

C     END DO

! Update the halos for the QT array
      CALL SWAPBOUNDS(QT,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,Q_LEVELS)
      ENDIF ! L_TRACER_THETAL_QT

CL END OF ROUTINE ADV_CTL

      RETURN
      END
