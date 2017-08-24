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
CLL  SUBROUTINE DYN_DIAG---------------------------------------
CLL
CLL PURPOSE: Calculate various diagnostics required for operational
CLL          and climate oputput. U and V compnts of wind on P levels
CLL          , clear air turbulence,
CLL          and MAXIMUM WIND Compnts and LEVEL.
CLL          and Potential Vorticty on isentropic surfaces
CLL
CLL JH, DR, RS  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.0   30/12/92  Two references to the former deckname DYNDIAG1
CLL                in CMESSAGE changed to subroutine name DYN_DIAG. MJH
CLL   3.1   25/01/93  Include test diagnostic (a simple analytic
CLL                function), items 231,232,233,234. R. Rawlins
CLL   3.1   14/01/93  Add routines to calculate potential vorticity on
CLL                   a pressure surface and theta on a pv surface.
CLL   3.4   26/05/94  Argument LLINTS added and passed to CALC_PV,
CLL                            CALC_PV_P, THETA_PV    S.J.Swarbrick
CLL  4.1   31/05/96     The number of v points to be processed on a
CLL                     C grid differs from u by row_length. u,v
CLL                     dimensioned separately in call to WLLTOEQ.
CLL                     Requirement for VAR.
CLL                     Author I.Edmond       Reviewer D. Goddard
!LL   4.2   08/01/97  Initialise PUV array to remove any NaNs in
!LL                   halo regions.                     P.Burton
!LL  4.4  09/04/97 : Add new diagnositics 235 qw, 236 heavyside 
!LL                  function (1 if pressure surface above land
!LL                  zero if below)        
!LL                  and 237 total kinetic energy per unit area. 
!LL       30/07/97 : Also geopotential height on u grid, 238 Z, 
!LL                  239 uZ, 240  VZ.
!LL       19/08/97 : 241 mountain torque per unit area.
!LL                  R. A. Stratton.                                   
!LL  4.5  15/04/98 Added start-end arguments to V_INT, V_INT_T and
!LL                V_INT_Z routines, and also to a lot of loops
!LL                over fields. Consequently, NS halos of diagnostics
!LL                are not set in this routine - instead STASHWORK is
!LL                initialised in ST_DIAG1.  S.D.Mullerworth
!LL       23/09/98   Allow 50m winds to be above second model level
!LL                                                        P.Burton
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  System components covered : D16D
CLL
CLL  Project task:
CLL
CLL  External documentation:
!LL            Unified Model Documentation Paper no D4                  
!LL            describes the diagnostics u*v, v*T etc.                  
CLL
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE DYN_DIAG(
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
C   primary data in
     &     PSTAR,U,V,Q,
     &     THETA,OROG,P_EXNER_HALF,PSTAR_OLD,
C   primary data constants
     &     U_ROWS,P_ROWS,ROW_LENGTH,P_LEVELS,Q_LEVELS,P_FIELD,
     &     U_FIELD,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     &     NMOST_LAT,WMOST_LONG,NS_SPACE,EW_SPACE,PHI_POLE,
     &     LAMBDA_POLE,SEC_U_LATITUDE,ROTATE_UV,ROTATE_MAX_UV,
     &     ELF,ETA_MATRIX_INV,MATRIX_P_O,LATITUDE_STEP_INVERSE,
     &     LONGITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,SEC_P_LATITUDE,
     &     COS_U_LATITUDE,F3,FORECAST_HRS,
C  Required Theta values
     &     DESIRED_THETA,PV_PRESS,DESIRED_PV,REQ_THETA_PV_LEVS,
     & n_levels,
C   Required pressures
     &     UCOMP_PRESS,VCOMP_PRESS,CAT_PROB_PRESS,T_PRESS,W_PRESS,
     &     Q_PRESS,TESTD_PRESS,HEAVY_PRESS,Z_PRESS,
C   Required Model levels
     &     TESTD_MODEL,
C   Indices for product fields
     &     UV_IND,UT_IND,VT_IND,T2_IND,U2_IND,V2_IND,WT_IND,WU_IND,
     &     WV_IND,QU_IND,QV_IND,QW_IND,UZ_IND,VZ_IND,
C   DIAGNOSTICS OUT
     &     UCOMP_P,VCOMP_P,MAX_CAT_PROB,MAX_CAT_LEVEL,
     &     CAT_PROB_SINGLE,MAX_WIND_HEIGHT,
     &     MAX_WIND_ICAO_HEIGHT,MAX_WIND_PRESSURE ,UCOMP_MAX_WIND,
     &     VCOMP_MAX_WIND,CAT_PROB_MEAN,UCOMP50_WIND,VCOMP50_WIND,
     &     POTN_VORT_THETA,
     &     UV_P,T_P,UT_P,VT_P,T2_P,U2_P,V2_P,W_P,WT_P,WU_P,WV_P,Q_P,
     &     UQ_P,VQ_P,
     &     POTN_VORT_ON_P,THETA_ON_PV,
     &     TESTDIAG1,TESTDIAG2,TESTDIAG3,TESTDIAG4,
     &     WQ_P,HEAVYSIDE_P,TOTAL_KE,                                   
     &     Z_P,UZ_P,VZ_P,M_TORQUE, 
C   diagnostic lengths
     &     UCOMP_P_LEVS,VCOMP_P_LEVS,CAT_PROB_LEVS,POTN_VORT_THETA_LEVS,
     &     POTN_VORT_P_LEVS,THETA_PV_LEVS,THETA_PV_P_LEVS,
     &     UV_P_LEVS,T_P_LEVS,
     &     UT_P_LEVS,VT_P_LEVS,T2_P_LEVS,U2_P_LEVS,V2_P_LEVS,W_P_LEVS,
     &     WT_P_LEVS,WU_P_LEVS,WV_P_LEVS,Q_P_LEVS,QU_P_LEVS,QV_P_LEVS,
     &     TESTD_P_LEVS,TESTD_M_LEVS, QW_P_LEVS,HEAVY_P_LEVS,           
     &     Z_P_LEVS,UZ_P_LEVS,VZ_P_LEVS,
C   diagnostic logical indicators
     &     QUCOMP_P,              QVCOMP_P,          QMAX_CAT_PROB,
     &     QMAX_CAT_LEVEL,        QCAT_PROB_SINGLE,  QMAX_WIND_HEIGHT,
     &     QMAX_WIND_ICAO_HEIGHT, QMAX_WIND_PRESSURE,QUCOMP_MAX_WIND,
     &     QVCOMP_MAX_WIND,       QCAT_PROB_MEAN,    QUCOMP50_WIND,
     &     QVCOMP50_WIND, QPOTN_VORT_THETA,
     &     QUV_P, QT_P, QUT_P, QVT_P, QT2_P,
     &     QU2_P, QV2_P, QW_P, QWT_P, QWU_P, QWV_P, QQ_P, QUQ_P, QVQ_P,
     &     QPOTN_VORT_PRESS, QTHETA_ON_PV,
     &     QDIA1,QDIA2,QDIA3,QDIA4,
     &     QWQ_P,QHEAVY_P,QTOTAL_KE,QZ_P,QUZ_P,QVZ_P,Q_MT,Z_REF,
C   diagnostic rerun code and message
     &     ICODE,CMESSAGE,
C   Logical switch LLINTS - passed to other routines
     &     LLINTS)
C
      IMPLICIT NONE
      LOGICAL  LLINTS
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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

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
      INTEGER
     *  P_FIELD            !IN    1ST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD            !IN    1ST DIMENSION OF FIELD OF U,V
     *, U_ROWS             !IN    NUMBER OF ROWS FOR U,V FIELD
     *, P_ROWS             !IN    NUMBER OF ROWS FOR P,T FIELD
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, LEVELS             !IN    NUMBER OF MODEL LEVELS
     *, P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS
     *, Q_LEVELS           !IN    NUMBER OF WET LEVELS
     *, MATRIX_P_O         !IN Order of polynomial used in calculation
     *                     !   of ETA_HALF inverse matrix
     *, FORECAST_HRS       !IN    FORECAST HOURS AFTER ANALYSIS T+nn
     *, ICODE              ! RETURN CODE      :    IRET=0   NORMAL EXIT
      INTEGER
     *  UCOMP_P_LEVS       !IN    NO OF LEVS ON WHICH TO INTERP U_P
     *, VCOMP_P_LEVS       !IN    NO OF LEVS ON WHICH TO INTERP V_P
     *, CAT_PROB_LEVS      !IN    NO OF LEVS ON WHICH TO CALC/INTERP CAT
     & ,n_levels      ! number of levels for dtheta/dp spline
     &, POTN_VORT_THETA_LEVS !IN  no of theta levs on which to calc pv
     &, POTN_VORT_P_LEVS   !IN    no of p levs on which to calc pv
     &, THETA_PV_LEVS      !IN    no of pv levs on which to calc theta
     &, THETA_PV_P_LEVS    !IN    no of p levs to calculate pv at,which
     &                     !     are then used to calculate theta on p.
     *, UV_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP UV_P
     *, T_P_LEVS           !IN    NO OF LEVS ON WHICH TO INTERP T_P
     *, UT_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP UT_P
     *, VT_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP VT_P
     *, T2_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP T2_P
     *, U2_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP U2_P
     *, V2_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP V2_P
     *, W_P_LEVS           !IN    NO OF LEVS ON WHICH TO INTERP W_P
     *, WT_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP WT_P
     *, WU_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP WU_P
     *, WV_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP WV_P
     *, Q_P_LEVS           !IN    NO OF LEVS ON WHICH TO INTERP Q_P
     *, QU_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP UQ_P
     *, QV_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP VQ_P
     *, QW_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP WQ_P    
     *, TESTD_P_LEVS       !IN    NO OF PRESS LEVELS TO CALC TESTDIAG
     *, TESTD_M_LEVS       !IN    NO OF MODEL LEVELS TO CALC TESTDIAG
     &, HEAVY_P_LEVS       !IN    NO OF PRESS levels to calc HEAVYSIDE  
     &, Z_P_LEVS           !IN    No of press levels for geopotential
     &, UZ_P_LEVS          !IN    No of press levels for UZ
     &, VZ_P_LEVS          !IN    No of press levels for VZ
     &, Z_REF              !IN   Level of model used to calculate PMSL
  
      INTEGER
     *  UV_IND(UV_P_LEVS,2) !IN  index for pressure levels in u and v
     *, UT_IND(UT_P_LEVS,2) !IN  index for pressure levels in U and T
     *, VT_IND(VT_P_LEVS,2) !IN  index for pressure levels in V and T
     *, T2_IND(T2_P_LEVS)   !IN  index for pressure levels in T
     *, U2_IND(U2_P_LEVS)   !IN  index for pressure levels in U
     *, V2_IND(V2_P_LEVS)   !IN  index for pressure levels in V
     *, WT_IND(WT_P_LEVS,2) !IN  index for pressure levels in W and T
     *, WU_IND(WU_P_LEVS,2) !IN  index for pressure levels in W and U
     *, WV_IND(WV_P_LEVS,2) !IN  index for pressure levels in W and V
     *, QU_IND(QU_P_LEVS,2) !IN  index for pressure levels in q and U
     *, QV_IND(QV_P_LEVS,2) !IN  index for pressure levels in Q and V
     *, QW_IND(QW_P_LEVS,2) !IN  index for pressure levels in Q and W   
     *, UZ_IND(UZ_P_LEVS,2) !IN  index for pressure levels in u and z  
     *, VZ_IND(VZ_P_LEVS,2) !IN  index for pressure levels in v and z  
C
      LOGICAL
     * QUCOMP_P      !IN  LOGICAL FLAG FOR PRESS INTER U COMPONENTS
     *,QVCOMP_P      !IN     "     "    "    "     "   V     "
     *,QMAX_CAT_PROB     !IN  "      "     "     MAXIMUM CAT PROBABILITY
     *,QMAX_CAT_LEVEL    !IN  "      "     "     MAX CAT PROB LEVEL
     *,QCAT_PROB_SINGLE  !IN  "      "     "     CAT PROB ON PRESS SFCE
     &,QMAX_WIND_HEIGHT !IN  "     "    "    MAX_WIND HEIGHT
     &,QMAX_WIND_ICAO_HEIGHT !IN   "    "    MAX_WIND  ICAO HEIGHT
     &,QMAX_WIND_PRESSURE !IN      "    "    MAX_WIND  PRESSURE
     &,QUCOMP_MAX_WIND  ! IN "     "    "    UCOMP_MAX _WIND
     &,QVCOMP_MAX_WIND  ! IN "     "    "    VCOMP_MAX _WIND
     *,QCAT_PROB_MEAN   ! IN  "      "     "  CAT PROB MEAN 300/250/200m
     &,ELF              ! IN  True if ELF i.e rotated LAM grid
     &,ROTATE_UV        ! IN  True if winds to be rotated
     &,ROTATE_MAX_UV    ! IN  True if winds to be rotated
     &,QUCOMP50_WIND    ! IN  Logical flag for QUCOMP50
     &,QVCOMP50_WIND    ! IN     "      "   "     "
     &,QPOTN_VORT_THETA ! IN     "      "   "     computing pv on theta
     &,QPOTN_VORT_PRESS ! IN     "      "   "     computing pv on
     &                  !                         pressure.
     &,QTHETA_ON_PV     ! IN     "      "   "     computing theta on PV
     &,QUV_P            ! IN     "      "   "  UV on pressure levels
     &,QT_P             ! IN     "      "   "  T  on pressure levels
     &,QUT_P            ! IN     "      "   "  UT on pressure levels
     &,QVT_P            ! IN     "      "   "  VT on pressure levels
     &,QT2_P            ! IN     "      "   "  T2 on pressure levels
     &,QU2_P            ! IN     "      "   "  U2 on pressure levels
     &,QV2_P            ! IN     "      "   "  V2 on pressure levels
     &,Qw_P             ! IN     "      "   "  w  on pressure levels
     &,QwT_P            ! IN     "      "   "  wT on pressure levels
     &,QwU_P            ! IN     "      "   "  wU on pressure levels
     &,QwV_P            ! IN     "      "   "  wV on pressure levels
     &,QQ_P             ! IN     "      "   "  q  on pressure levels
     &,Quq_P            ! IN     "      "   "  uq on pressure levels
     &,Qvq_P            ! IN     "      "   "  vq on pressure levels
     &,Qwq_P            ! IN     "      "   "  wq on pressure levels    
     &,QDIA1            ! IN     "      "   "  test diagnostic 1
     &,QDIA2            ! IN     "      "   "  test diagnostic 2
     &,QDIA3            ! IN     "      "   "  test diagnostic 3
     &,QDIA4            ! IN     "      "   "  test diagnostic 4
     &,QHEAVY_P         ! IN     "      "   "  heavyside function p lev 
     &,QTOTAL_KE        ! IN     "      "   "  Total KE                 
     &,QZ_P             ! IN     "      "   "  Z on pressure levels
     &,QUZ_P            ! IN     "      "   "  UZ on pressure levels
     &,QVZ_P            ! IN     "      "   "  VZ on pressure levels
     &,Q_MT             ! IN      mountain torque per unit area
C
      CHARACTER CMESSAGE*(*)

      REAL
     * PSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR PSTAR FIELD
     *,PSTAR_OLD(P_FIELD)     !IN    Pstar before dynamics.
     *,OROG(P_FIELD)          !IN    PRIMARY MODEL OROGRAPHY
     *,P_EXNER_HALF(P_FIELD,P_LEVELS+1) !IN  EXNER PRESS ON 1/2 LVLS
     *,THETA(P_FIELD,P_LEVELS)!IN PRIMARY MODEL ARRAY FOR THETA FIELD
     *,U(U_FIELD,P_LEVELS)    !INT PRIMARY MODEL ARRAY FOR U FIELD
     *,V(U_FIELD,P_LEVELS)    !IN PRIMARY MODEL ARRAY FOR V FIELD
     *,Q(P_FIELD,Q_LEVELS)    !IN PRIMARY MODEL ARRAY FOR HUMIDITY
C   DIAGNOSTICS OUT
     *,UCOMP_P(U_FIELD,UCOMP_P_LEVS)  !OUT  UCOMP ON ANY PRESSURE SFCE
     *,VCOMP_P(U_FIELD,VCOMP_P_LEVS)  !OUT  VCOMP ON ANY PRESSURE SFCE
     *,UCOMP50_WIND(U_FIELD)       !OUT 50 M wind zonal cmpnt.
     *,VCOMP50_WIND(U_FIELD)       !OUT 50 M wind zonal cmpnt.
     *,MAX_CAT_PROB(U_FIELD) !OUT MAX CAT PROB FROM LEVELS 300/250/200mb
     *,MAX_CAT_LEVEL(U_FIELD)!OUT LEVEL OF MAX CAT PROB
     *,CAT_PROB_SINGLE(U_FIELD,CAT_PROB_LEVS)!OUT CAT PROB ON PRESS SFCE
     *,MAX_WIND_HEIGHT(U_FIELD) !OUT HEIGHT LEVEL OF MAX WIND
     *,MAX_WIND_ICAO_HEIGHT(U_FIELD) !OUT ICAO HEIGHT LEVEL OF MAX WIND
     *,MAX_WIND_PRESSURE(U_FIELD) !OUT PRESSURE LEVEL OF MAX WIND
     *,UCOMP_MAX_WIND(U_FIELD) !OUT U COMPONENT OF  MAX WIND
     *,VCOMP_MAX_WIND(U_FIELD) !OUT V COMPONENT OF  MAX WIND
     *,CAT_PROB_MEAN(U_FIELD)!OUT CAT PROB MEAN OVER LEVELS 300/250/200m
     &,POTN_VORT_THETA(P_FIELD,POTN_VORT_THETA_LEVS) !OUT pv on theta
     &,POTN_VORT_ON_P(P_FIELD,POTN_VORT_P_LEVS) !OUT pv on pressure
     &,THETA_ON_PV(P_FIELD,THETA_PV_LEVS) !OUT Pot. temp on a pv surface
     &,e_levels(n_levels)      ! Model half-levels for dtheta/dp.
     &,dthe_dph(p_field,n_levels)   ! dtheta/dp for potential vorticity
     &                             !    on half-levels.
     *,UV_P(U_FIELD,UV_P_LEVS)     ! UV on pressure levels, wind grid
     *,T_P(U_FIELD,T_P_LEVS)       ! T  on pressure levels, wind grid
     *,UT_P(U_FIELD,UT_P_LEVS)     ! UT on pressure levels, wind grid
     *,VT_P(U_FIELD,VT_P_LEVS)     ! VT on pressure levels, wind grid
     *,T2_P(U_FIELD,T2_P_LEVS)     ! T2 on pressure levels, wind grid
     *,U2_P(U_FIELD,U2_P_LEVS)     ! U2 on pressure levels, wind grid
     *,V2_P(U_FIELD,V2_P_LEVS)     ! V2 on pressure levels, wind grid
     *,W_P(U_FIELD,W_P_LEVS)       ! w  on pressure levels, wind grid
     *,WT_P(U_FIELD,WT_P_LEVS)     ! wT on pressure levels, wind grid
     *,WU_P(U_FIELD,WU_P_LEVS)     ! wU on pressure levels, wind grid
     *,WV_P(U_FIELD,WV_P_LEVS)     ! wV on pressure levels, wind grid
     *,Q_P(U_FIELD,Q_P_LEVS)       ! q  on pressure levels, wind grid
     *,UQ_P(U_FIELD,QU_P_LEVS)     ! qu on pressure levels, wind grid
     *,VQ_P(U_FIELD,QV_P_LEVS)     ! qv on pressure levels, wind grid
     *,WQ_P(U_FIELD,QW_P_LEVS)     ! qw on pressure levels, wind grid   
     *,TESTDIAG1(U_FIELD)             ! OUT Diag 1 single lev, u grid
     *,TESTDIAG2(P_FIELD)             ! OUT Diag 2 single lev, p grid
     *,TESTDIAG3(P_FIELD,TESTD_P_LEVS)! OUT Diag 3 press levs, p grid
     *,TESTDIAG4(P_FIELD,TESTD_M_LEVS)! OUT Diag 4 model levs, p grid
     &,HEAVYSIDE_P(U_FIELD,HEAVY_P_LEVS) ! OUT heavyside on p levs      
     &,TOTAL_KE(U_FIELD)           ! total KE per unit area, u grid     
     &,Z_P(U_FIELD,Z_P_LEVS)       ! z on pressure levels, u grid
     &,UZ_P(U_FIELD,UZ_P_LEVS)     ! Uz on pressure levels, u grid
     &,VZ_P(U_FIELD,VZ_P_LEVS)     ! Vz on pressure levels, u grid
     &,M_TORQUE(U_FIELD)           ! mountain torque per unit area, u
C
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,
      REAL
     * AKH(P_LEVELS+1)       !IN    LAYER THICKNESS
     *,BKH(P_LEVELS+1)       !IN    LAYER THICKNESS
     *,AK (P_LEVELS)         !IN    VALUE AT LAYER CENTRE
     *,BK (P_LEVELS)         !IN    VALUE AT LAYER CENTRE
     *,DELTA_AK (P_LEVELS)   !IN
     *,DELTA_BK (P_LEVELS)   !IN
     *,NMOST_LAT             !Northern most latitude of grid
     *,WMOST_LONG            !Western most longitude
     *,EW_SPACE              !Delta longitude
     *,NS_SPACE              !Delta latitude
     *,PHI_POLE              !Latitude of the pseudo pole
     *,LAMBDA_POLE           !Longitude of the pseudo pole
     *,SEC_U_LATITUDE(U_FIELD)!IN 1/COS(LAT) AT U POINTS
     *,ETA_MATRIX_INV(MATRIX_P_O,MATRIX_P_O,P_LEVELS)!IN Inverse matrix
     *                                               !   of ETA_HALF
     *,UCOMP_PRESS(UCOMP_P_LEVS)    !IN Required pressure surface
     *,VCOMP_PRESS(VCOMP_P_LEVS)    !IN Required pressure surface
     *,CAT_PROB_PRESS(CAT_PROB_LEVS)!IN     "       "       "
     &,DESIRED_THETA(POTN_VORT_THETA_LEVS) !IN required theta surfaces
     &,PV_PRESS(POTN_VORT_P_LEVS)   !IN required pressure surfaces
     &,DESIRED_PV(THETA_PV_LEVS)    !IN required pv surfaces
     &,REQ_THETA_PV_LEVS(THETA_PV_P_LEVS) !IN required p surfaces.
     *,T_PRESS(T_P_LEVS)            !IN     "       "       "
     *,W_PRESS(W_P_LEVS)            !IN     "       "       "
     *,Q_PRESS(Q_P_LEVS)            !IN     "       "       "
     *,TESTD_PRESS(TESTD_P_LEVS)    !IN Required pressures for test diag
     &,HEAVY_PRESS(HEAVY_P_LEVS)    !IN Required pressures for heavyside
     &,Z_PRESS(Z_P_LEVS)            !IN Required pressures for geopot
     *,TESTD_MODEL(TESTD_M_LEVS)    !IN Required model lvs for test diag
     *,LATITUDE_STEP_INVERSE        !IN 1/latitude increment
     *,LONGITUDE_STEP_INVERSE       !IN 1/longitude increment
     *,ADVECTION_TIMESTEP           !IN advection timestep
     *,SEC_P_LATITUDE(P_FIELD)      !IN 1/cos(lat) p points
     *,COS_U_LATITUDE(U_FIELD)      !IN cos(lat) u points
     &,F3(U_FIELD)                  !IN    Coriolis term.

C Local variables

      LOGICAL
     &  found_levels  ! TRUE if level search is successful

      INTEGER
     &  level1,level2  ! Model levels either side of 50m

C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS:
C   REAL ARRAYS REQUIRED AT FULL FIELD LENGTH
C   1 INTEGER INDEX ARRAY
C
C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL TROP,V_INT,WINDMAX,ICAO_HT,V_INT_T,OMEGA_DIAG,
     &         P_TO_UV,CAT,CALC_PV,CALC_PV_P,THETA_PV
     &       ,V_INT_ZH,V_INT_Z
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS+1) * ROWLENGTH
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      LOGICAL
     *  TEST          !     NUMBER OF P POINTS NEEDED
      REAL
     *  PHI_STAR(P_FIELD)
     *, P(P_FIELD,P_LEVELS)   ! PRESSURE ARRAY
     *, PUV(U_FIELD,P_LEVELS) ! PRESSURE ARRAY ON U,V POINTS
     *, PZ(P_FIELD)           ! PRESSURE SURFACE ON WHICH RESULTS REQD
     *, WORK1(U_FIELD)        ! Work array
     *, WORK5(P_FIELD)        ! Work array
     *, ETA1,ETA2,ETA50,C1,C2 ! Used in the calculation of 50 M winds
     *, OMEGA(P_FIELD,P_LEVELS)   ! Omega array
     &, THETA_ON_P(P_FIELD,POTN_VORT_P_LEVS) !holds Pot. temperature on
     &                                       !on a pressure surface.
     & ,model_half_height(p_field,p_levels+1) !heights on model half lev
     & ,PSTAR_UV(U_FIELD)     ! pstar on uv grid.                       
     & ,FACTOR                ! factor for KE calculation               
     & ,PLEV                  ! pressure level for Heavyside calculation
C
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      INTEGER    K,I,II,IK,J,LEVEL! LOOP COUNTERS IN ROUTINE
     *, T_REF    ! reference level for below surface T extrapolation.
     &, U_FLD_VALID ! Set to No of points in U-field excluding
                    ! unused rows and halos
C
           ICODE=0
           U_FLD_VALID=LAST_U_FLD_PT-FIRST_FLD_PT+1

C-----------------------------------------------------------------------
CL    Calculate variables required by various subroutines
C-----------------------------------------------------------------------
      IF(QUCOMP_P.OR.QVCOMP_P.OR.QCAT_PROB_SINGLE.OR.QCAT_PROB_MEAN.OR.
     &   QMAX_CAT_PROB.OR.QMAX_CAT_LEVEL.OR.QT_P.OR.QW_P)THEN
        DO K=1,P_LEVELS
          DO I=1,P_FIELD
            P(I,K)=AK(K)+BK(K)*PSTAR(I)
          ENDDO
        ENDDO
! QAN fix : Initialise unused rows
        IF (at_base_of_LPG) THEN
          DO K=1,P_LEVELS
            DO I=LAST_U_VALID_PT,U_FIELD
              PUV(I,K)=0.0
            ENDDO
          ENDDO
        ENDIF
        DO K=1,P_LEVELS
          CALL P_TO_UV(P(1,K),PUV(1,K),P_FIELD,U_FIELD,ROW_LENGTH,
     &       P_ROWS)
        ENDDO
        CALL SWAPBOUNDS(PUV,ROW_LENGTH,U_ROWS,
     &    EW_Halo,NS_Halo,P_LEVELS)
      ENDIF
C-----------------------------------------------------------------------
      IF(QUCOMP_P.OR.QVCOMP_P.OR.QZ_P) THEN 
        DO I=1,P_FIELD
          PHI_STAR(I)=OROG(I)*G
        ENDDO
C
CL------------------Interpolate U cmpnt of wind onto Pressure ---------

          IF(QUCOMP_P) THEN
            DO  K=1,UCOMP_P_LEVS
              DO  I=FIRST_FLD_PT,LAST_U_FLD_PT
                PZ(I)=UCOMP_PRESS(K)*100.0   ! convert to Pascals
              ENDDO
              CALL V_INT(PUV,PZ,U,UCOMP_P(1,K),U_FIELD,P_LEVELS,
     &        UCOMP_MAX_WIND,MAX_WIND_PRESSURE,.FALSE.
     &        ,FIRST_FLD_PT,LAST_U_FLD_PT)
            ENDDO  ! Levels loop
          ENDIF
CL------------------Interpolate V cmpnt of wind onto Pressure ---------

          IF(QVCOMP_P) THEN
            DO  K=1,VCOMP_P_LEVS
              DO  I=FIRST_FLD_PT,LAST_U_FLD_PT
                PZ(I)=VCOMP_PRESS(K)*100.0   ! convert to Pascals
              ENDDO
              CALL V_INT(PUV,PZ,V,VCOMP_P(1,K),U_FIELD,P_LEVELS,
     &        VCOMP_MAX_WIND,MAX_WIND_PRESSURE,.FALSE.
     &        ,FIRST_FLD_PT,LAST_U_FLD_PT)
            ENDDO  ! Levels loop
          ENDIF
      ENDIF     ! End of QUCOMP or QVCOMP
CL------------------Calculate the maximum wind-------------------------

      IF (QUCOMP_MAX_WIND.AND.
     *    QVCOMP_MAX_WIND.AND.
     *    QMAX_WIND_PRESSURE) THEN
          CALL WINDMAX(
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
     &     PSTAR,U,V,U_ROWS,P_ROWS,ROW_LENGTH,P_LEVELS,
     *     P_FIELD,U_FIELD,AK,BK,AKH,BKH,ETA_MATRIX_INV,MATRIX_P_O,
     *     UCOMP_MAX_WIND,VCOMP_MAX_WIND,MAX_WIND_PRESSURE)
CL---   ICAO HT of the MAX WIND pressure  ?
        IF(QMAX_WIND_ICAO_HEIGHT) THEN
          CALL ICAO_HT(MAX_WIND_PRESSURE(FIRST_FLD_PT),U_FLD_VALID
     &      ,MAX_WIND_ICAO_HEIGHT(FIRST_FLD_PT))
        ENDIF
      ELSEIF(QUCOMP_MAX_WIND.NEQV.QVCOMP_MAX_WIND.OR.QUCOMP_MAX_WIND.
     *   NEQV.QMAX_WIND_PRESSURE.OR.QVCOMP_MAX_WIND.NEQV.
     *   QMAX_WIND_PRESSURE)THEN
        WRITE(6,*)' Subroutine WINDMAX not called - U & VCOMP_MAX_WIND'
        WRITE(6,*)'   and MAX_WIND_PRESSURE all must be selected'
      ENDIF  ! Top IF block for maxwind


C-----------------------------------------------------------------------
CL    Section 15 Item 205 CAT PROBABILITY at pressure levels
C-----------------------------------------------------------------------
      IF(QCAT_PROB_SINGLE.AND.QMAX_WIND_PRESSURE.AND.QUCOMP_MAX_WIND.
     *   AND.QVCOMP_MAX_WIND)THEN
        DO K=1,CAT_PROB_LEVS
          DO I=FIRST_VALID_PT,LAST_U_VALID_PT
            PZ(I)=CAT_PROB_PRESS(K)*100.0      ! Convert to pascals
          ENDDO
          CALL CAT(
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
     &      glsize(2),
     &      U,V,PUV,PSTAR,PZ,MAX_WIND_PRESSURE,
     &      CAT_PROB_SINGLE(1,K),P_FIELD,U_FIELD,P_LEVELS,ROW_LENGTH,
     &      P_ROWS,SEC_U_LATITUDE,AK,BK,EW_SPACE,NS_SPACE)
        ENDDO
      ELSEIF(QCAT_PROB_SINGLE.AND.(.NOT.(QMAX_WIND_PRESSURE.AND.
     *  QUCOMP_MAX_WIND.AND.QVCOMP_MAX_WIND)))THEN
        WRITE(6,*)' Subroutine CAT not called - PRESSURE, U & VCOMP of'
      WRITE(6,*)' MAX WIND must be selected as well as CAT_PROB_SINGLE'
      ENDIF
C-----------------------------------------------------------------------
CL    Section 15 Item 211 MEAN CAT PROBABILITY over levels 300,250,200mb
C-----------------------------------------------------------------------
      IF(QCAT_PROB_MEAN.AND.QMAX_WIND_PRESSURE.AND.QUCOMP_MAX_WIND.
     *   AND.QVCOMP_MAX_WIND)THEN
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          CAT_PROB_MEAN(I)=0.0
        ENDDO
C-----------------------------------------------------------------------
CL    Call CAT for three levels - CAT probability output to WORK1
C-----------------------------------------------------------------------
        DO K=1,3
          IF(K.EQ.1)PZ(FIRST_VALID_PT)=30000.0
          IF(K.EQ.2)PZ(FIRST_VALID_PT)=25000.0
          IF(K.EQ.3)PZ(FIRST_VALID_PT)=20000.0
          DO I=FIRST_VALID_PT+1,LAST_U_VALID_PT
            PZ(I)=PZ(FIRST_VALID_PT)
          ENDDO
          CALL CAT(
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
     &      glsize(2),
     &      U,V,PUV,PSTAR,PZ,MAX_WIND_PRESSURE,WORK1,P_FIELD,
     &      U_FIELD,P_LEVELS,ROW_LENGTH,P_ROWS,SEC_U_LATITUDE,AK,BK,
     &      EW_SPACE,NS_SPACE)
          DO I=FIRST_FLD_PT,LAST_U_FLD_PT
            CAT_PROB_MEAN(I)=CAT_PROB_MEAN(I)+WORK1(I)
          ENDDO
        ENDDO
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          CAT_PROB_MEAN(I)=CAT_PROB_MEAN(I)/3.0
        ENDDO
      ELSEIF(QCAT_PROB_MEAN.AND.(.NOT.(QMAX_WIND_PRESSURE.AND.
     *  QUCOMP_MAX_WIND.AND.QVCOMP_MAX_WIND)))THEN
        WRITE(6,*)' Subroutine CAT not called - PRESSURE, U & VCOMP of'
        WRITE(6,*)' MAX WIND must be selected as well as CAT_PROB_MEAN'
      ENDIF
C-----------------------------------------------------------------------
CL    Section 15 Items 203/204 MAXIMUM CAT PROBABILITY AND LEVEL
C-----------------------------------------------------------------------
      IF(QMAX_CAT_PROB.AND.QMAX_CAT_LEVEL.AND.QMAX_WIND_PRESSURE.AND.
     *   QUCOMP_MAX_WIND.AND.QVCOMP_MAX_WIND)THEN
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          MAX_CAT_PROB(I)=0.0
          MAX_CAT_LEVEL(I)=30000.0
        ENDDO
C-----------------------------------------------------------------------
CL    Call CAT for three levels - CAT probability output to WORK1
C-----------------------------------------------------------------------
        DO K=1,3
          IF(K.EQ.1)PZ(FIRST_VALID_PT)=30000.0
          IF(K.EQ.2)PZ(FIRST_VALID_PT)=25000.0
          IF(K.EQ.3)PZ(FIRST_VALID_PT)=20000.0
          DO I=FIRST_VALID_PT+1,LAST_U_VALID_PT
            PZ(I)=PZ(FIRST_VALID_PT)
          ENDDO
          CALL CAT(
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
     &      glsize(2),
     &      U,V,PUV,PSTAR,PZ,MAX_WIND_PRESSURE,WORK1,P_FIELD,
     &      U_FIELD,P_LEVELS,ROW_LENGTH,P_ROWS,SEC_U_LATITUDE,AK,BK,
     &      EW_SPACE,NS_SPACE)
          IF(K.EQ.1)PZ(1)=30000.0
          IF(K.EQ.2)PZ(1)=25000.0
          IF(K.EQ.3)PZ(1)=20000.0
          DO I=FIRST_FLD_PT,LAST_U_FLD_PT
            IF(MAX_CAT_PROB(I).LT.WORK1(I))THEN
              MAX_CAT_PROB(I)=WORK1(I)
              MAX_CAT_LEVEL(I)=PZ(1)
            ENDIF
          ENDDO
        ENDDO
      ELSEIF(QMAX_CAT_PROB.NEQV.QMAX_CAT_LEVEL)THEN
        WRITE(6,*)' Subroutine CAT not called - both MAX_CAT_PROB and'
        WRITE(6,*)' MAX_CAT_LEVEL must be selected'
      ELSEIF(QMAX_CAT_PROB.NEQV.QMAX_CAT_LEVEL)THEN
        WRITE(6,*)' Subroutine CAT not called - MAX_CAT_PROB and'
        WRITE(6,*)'       MAX_CAT_LEVEL must both be selected'
      ELSEIF(QMAX_CAT_PROB.AND.QMAX_CAT_LEVEL.AND.
     * (.NOT.(QMAX_WIND_PRESSURE.AND.QUCOMP_MAX_WIND.AND.
     * QVCOMP_MAX_WIND)))THEN
        WRITE(6,*)' Subroutine CAT not called - PRESSURE, U & VCOMP of'
        WRITE(6,*)' MAX WIND must be selected as well as MAX_CAT_PROB '
        WRITE(6,*)'                                    & MAX_CAT_LEVEL'
      ENDIF
C-----------------------------------------------------------------------
CL    Section 15 Items 212/213 50 M U and V components.
C-----------------------------------------------------------------------
C=====================================================================C
C     50 METRE WINDS                                                  C
C     USE U(50)=C1*U(ETA2)+C2(ETA1)                                   C
C     C1=LOG(ETA50/ETA1)/LOG(ETA2/ETA1) = 0.135                       C
C     C2=LOG(ETA2/ETA50)/LOG(ETA2/ETA1) = 0.865                       C
C     ETA50=0.994 IE CORRESPONDS TO Z=50M AND TBAR=283K and assumes   C
C     the first 5 levels are sigma levels.                            C
C=====================================================================C
C     First check between which levels the 50 M level lies (ETA=0.994)
      IF(QUCOMP50_WIND.OR.QVCOMP50_WIND) THEN
        ETA50=0.994
      found_levels=.FALSE.
      K=1
      DO WHILE ((.NOT. found_levels) .AND. (K .LT. P_LEVELS))
        level1=K
        level2=K+1
        ETA1=AK(level1)/PREF+BK(level1)
        ETA2=AK(level2)/PREF+BK(level2)

        IF ((ETA1 .GE. ETA50) .AND. (ETA2 .LE. ETA50))
     &    found_levels=.TRUE.

        K=K+1
      ENDDO

      IF (.NOT. found_levels) THEN
        ICODE=1
        CMESSAGE='DYN_DIAG: Error in calculating 50 M winds'
        RETURN
      ENDIF
        IF(ETA1.LT.ETA50) THEN
          ICODE=1
          CMESSAGE='DYN_DIAG: Error in calculating 50 M winds'
          RETURN
        ENDIF
        IF(ETA2.GT.ETA50) THEN
          ICODE=1
          CMESSAGE='DYN_DIAG: Error in calculating 50 M winds'
          RETURN
        ENDIF
        C1=ALOG(ETA50/ETA1)/ALOG(ETA2/ETA1)
        C2=ALOG(ETA2/ETA50)/ALOG(ETA2/ETA1)
      ENDIF
      IF(QUCOMP50_WIND) THEN
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          UCOMP50_WIND(I)=C1*U(I,level2)+C2*U(I,level1)
        ENDDO
      ENDIF
      IF(QVCOMP50_WIND) THEN
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          VCOMP50_WIND(I)=C1*V(I,level2)+C2*V(I,level1)
        ENDDO
      ENDIF
C want to compute pv for some theta level
      if(qpotn_vort_theta.or.qpotn_vort_press.or.qtheta_on_pv)then
      n_levels=p_levels-1
      call dthe_dp(pstar,theta,p_field,p_levels
     2          ,ak,bk,akh,bkh,n_levels
     3          ,e_levels,dthe_dph)
      endif
      IF (QPOTN_VORT_THETA) THEN
        DO I = 1,POTN_VORT_THETA_LEVS
           CALL CALC_PV
     1                 (PSTAR,THETA,U,V,P_FIELD,U_FIELD,P_LEVELS,
     2                  ROW_LENGTH,
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
     &                  RMDI,AK,BK,DESIRED_THETA(I),F3,
     &                  e_levels,n_levels,dthe_dph,
     3                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     4                  COS_U_LATITUDE,SEC_P_LATITUDE,
     5                  POTN_VORT_THETA(1,I),LLINTS)
        ENDDO
      ENDIF

C want to compute pv for some pressure level
      IF (QPOTN_VORT_PRESS) THEN
        DO I = 1,POTN_VORT_P_LEVS
           CALL CALC_PV_P
     1                 (PSTAR,THETA,U,V,P_FIELD,U_FIELD,P_LEVELS,
     2                  ROW_LENGTH,
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
     &                  RMDI,AK,BK,PV_PRESS(I),F3,
     & e_levels,n_levels,dthe_dph,
     3                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     4                  COS_U_LATITUDE,SEC_P_LATITUDE,
     5                  POTN_VORT_ON_P(1,I),THETA_ON_P(1,I),LLINTS)
        ENDDO
      END IF
C want to compute theta on some pv surfaces.
C the loop over the surfaces is contained inside the subroutine.
      IF (QTHETA_ON_PV) THEN
          CALL THETA_PV(
     1                  PSTAR,THETA,U,V,P_FIELD,U_FIELD,P_LEVELS,
     2                  ROW_LENGTH,
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
     & RMDI,AK,BK,F3,
     & e_levels,n_levels,dthe_dph,
     3                  THETA_PV_LEVS,DESIRED_PV,
     4                  THETA_PV_P_LEVS,REQ_THETA_PV_LEVS,
     5                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                  COS_U_LATITUDE,SEC_P_LATITUDE,
     7                  THETA_ON_PV,LLINTS)
      ENDIF

C-----------------------------------------------------------------------
CL  Section 15 item 215,  UV on pressure levels
CL  Only possible if U and V have been requested on pressure levels
CL  required for UV
CL
      IF (QUV_P) THEN
           do K=1,UV_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               UV_P(I,K) = UCOMP_P(I,UV_IND(K,1))*VCOMP_P(I,UV_IND(K,2))
             ENDDO
           ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 216,  T on pressure levels on the wind grid
CL   T is first calculated on p-grid for the pressure level and
CL   then interpolated to the u-grid.
CL
      IF (QT_P) THEN
        T_REF=2       !used in vertical interpolation
        DO K=1,T_P_LEVS
          DO I=FIRST_FLD_PT,P_FIELD
            PZ(I)=T_PRESS(K)*100.0    ! convert to pascals
          ENDDO
          CALL V_INT_T(WORK5,PZ,P(1,T_REF),PSTAR,P_EXNER_HALF
     &      ,THETA,P_FIELD,P_LEVELS,T_REF,AKH,BKH
     &      ,FIRST_FLD_PT,P_FIELD)
          CALL P_TO_UV(WORK5(FIRST_FLD_PT),T_P(FIRST_FLD_PT,K)
     &      ,P_FIELD-FIRST_FLD_PT+1,U_FIELD-FIRST_FLD_PT+1
     &      ,ROW_LENGTH,P_LAST_ROW)
        ENDDO
C-----------------------------------------------------------------------
CL  Section 15 item 217,  UT on pressure levels
CL
        IF (QUT_P) THEN
           do K=1,UT_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               UT_P(I,K) = UCOMP_P(I,UT_IND(K,1))*T_P(I,UT_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 218,  VT on pressure levels
CL
        IF (QVT_P) THEN
           do K=1,VT_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               VT_P(I,K) = VCOMP_P(I,VT_IND(K,1))*T_P(I,VT_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 219,  T**2 on pressure levels
CL
        IF (QT2_P) THEN
           DO K=1,T2_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               T2_P(I,K) = T_P(I,T2_IND(K))*T_P(I,T2_IND(K))
             ENDDO
           ENDDO
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 220,  U2 on pressure levels
CL  Only possible if U has been requested on the same pressure levels
CL
      IF (QU2_P) THEN
           DO K=1,U2_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               U2_P(I,K) = UCOMP_P(I,U2_IND(K))*UCOMP_P(I,U2_IND(K))
             ENDDO
           ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 221,  v2 on pressure levels
CL  Only possible if v has been requested on the same pressure levels
CL
      IF (QV2_P) THEN
           DO K=1,V2_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               V2_P(I,K) = VCOMP_P(I,V2_IND(K))*VCOMP_P(I,V2_IND(K))
             ENDDO
           ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 222,  w on pressure levels and wind grid
CL
CL
      IF (QW_P) THEN
        CALL OMEGA_DIAG(
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
     &                  U,V,OMEGA,SEC_P_LATITUDE,COS_U_LATITUDE,
     1                  PSTAR,PSTAR_OLD,DELTA_AK,DELTA_BK,
     2                  AK,BK,AKH,BKH,U_FIELD,P_FIELD,P_LEVELS,
     3                  ROW_LENGTH,LATITUDE_STEP_INVERSE,
     4                  LONGITUDE_STEP_INVERSE,ADVECTION_TIMESTEP)

C  omega returned at u points on all p_levels
C   Interpolation of omega to pressure levels
       DO K=1,W_P_LEVS
         DO I=FIRST_FLD_PT,LAST_U_FLD_PT
           PZ(I) = W_PRESS(K)*100.0   ! convert to Pascals
         ENDDO
         CALL V_INT(PUV,PZ,OMEGA,W_P(1,K),U_FIELD,P_LEVELS,WORK1,WORK5,
     &                       .FALSE.,FIRST_FLD_PT,LAST_U_FLD_PT)
       ENDDO

C-----------------------------------------------------------------------
CL  Section 15 item 223-225,  wT, wu, wv on pressure levels
CL  Only possible if w and T have been requested on the same pressure
CL levels
        IF (QWT_P) THEN
           do K=1,WT_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               WT_P(I,K) = W_P(I,WT_IND(K,1))*T_P(I,WT_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
        IF (QWU_P) THEN
           do K=1,WU_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               WU_P(I,K) = W_P(I,WU_IND(K,1))*UCOMP_P(I,WU_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
        IF (QWV_P) THEN
           do K=1,WV_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               WV_P(I,K) = W_P(I,WV_IND(K,1))*VCOMP_P(I,WV_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
CL  Section 15 item 226-228,  q, qu, qv on pressure levels
CL             item 235       qw on pressure levels                     
CL
      IF (QQ_P) THEN
        DO K=1,Q_P_LEVS
          DO I=FIRST_FLD_PT,P_FIELD
           PZ(I) = Q_PRESS(K)*100.0   ! convert to Pascals
          ENDDO
          CALL V_INT(P,PZ,Q,WORK5,P_FIELD,Q_LEVELS,WORK1,WORK1,.FALSE.
     &      ,FIRST_FLD_PT,P_FIELD)
          CALL P_TO_UV(WORK5(FIRST_FLD_PT),Q_P(FIRST_FLD_PT,K)
     &      ,P_FIELD-FIRST_FLD_PT+1,U_FIELD-FIRST_FLD_PT+1,ROW_LENGTH
     &      ,P_LAST_ROW)
        ENDDO

        IF (QuQ_P) THEN
           do K=1,QU_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               UQ_P(I,K) = Q_P(I,QU_IND(K,1))*UCOMP_P(I,QU_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
        IF (QVQ_P) THEN
           do K=1,QV_P_LEVS
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               VQ_P(I,K) = Q_P(I,QV_IND(K,1))*VCOMP_P(I,QV_IND(K,2))
             ENDDO
           ENDDO
        ENDIF
! added later therefore jump in stashcode                               
        IF (QWQ_P) THEN  
           do K=1,QW_P_LEVS  
             DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               WQ_P(I,K) = Q_P(I,QW_IND(K,1))*W_P(I,QW_IND(K,2))  
             ENDDO 
           ENDDO  
        ENDIF   
      ENDIF
C ---------------------------------------------------------------------
CL  Section 15 items 231,232,233,234 test diagnostics
CL   231 single level on u grid, 232 single level on p grid
CL   233 press levels on p grid, 234 model levels on p grid
CL
      IF (QDIA1.OR.QDIA2.OR.QDIA3.OR.QDIA4) THEN
        CALL TESTDIAG(
     1  P_FIELD,U_FIELD,P_ROWS,U_ROWS,ROW_LENGTH,EW_SPACE,NS_SPACE,
     2  NMOST_LAT,WMOST_LONG,ELF,PHI_POLE,LAMBDA_POLE,
     3  TESTD_PRESS,TESTD_P_LEVS,
     4  TESTD_MODEL,TESTD_M_LEVS,FORECAST_HRS,
     5  TESTDIAG1,TESTDIAG2,TESTDIAG3,TESTDIAG4,
     6  QDIA1,QDIA2,QDIA3,QDIA4)
C
      ENDIF                                                             
! --------------------------------------------------------------------- 
!L Items 236 and 237 both require pstar on the uv grid                  
!L                                                                      
      IF (QHEAVY_P.or.QTOTAL_KE) THEN                                   
        CALL P_TO_UV(PSTAR,PSTAR_UV,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)  
                                                                        
!L Section 15 item 236 Heavyside function on pressure levels for        
!L                     u-grid.                                          
!L The Heavyside function is defined as 1.0 if the pressure level       
!L  is above the surface (i.e. pstar) and 0.0 if below. A time mean of  
!L  this will give information on the fraction of time a pressure       
!L  level is above the land or sea surface.                             
                                                                        
        IF (QHEAVY_P) THEN                                              
          DO K=1,HEAVY_P_LEVS                                           
            PLEV=HEAVY_PRESS(K)*100.   ! pressure in Pascals            
            DO I=FIRST_FLD_PT,LAST_U_FLD_PT
              IF (PSTAR_UV(I).LT.PLEV) THEN                             
                 HEAVYSIDE_P(I,K)=0.0                                   
              ELSE                                                      
                 HEAVYSIDE_P(I,K)=1.0                                   
              ENDIF                                                     
            ENDDO                                                       
          ENDDO                                                         
        ENDIF                                                           
! --------------------------------------------------------------------- 
!L Section 15 item 237 Total kinetic energy in a column u-grid          
!L                                                                      
!L  KE = SUM [0.5/g (u*u + v*v)] dp over model levels.                  
!L                                                                      
!L  Output scaled by 1.0e-6 to prevent accuracy problems                
                                                                        
        IF (QTOTAL_KE) THEN                                             
          FACTOR=0.5*1.0e-6/g                                           
          DO I=FIRST_FLD_PT,LAST_U_FLD_PT
            TOTAL_KE(I)=0.0                                             
          ENDDO                                                         
          DO K=1,P_LEVELS                                               
            DO I=FIRST_FLD_PT,LAST_U_FLD_PT
              TOTAL_KE(I)=TOTAL_KE(I) -                                 
     &                         FACTOR*(U(I,K)*U(I,K)+V(I,K)*V(I,K))     
     &                        *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I))    
            ENDDO                                                       
          ENDDO                                                         
        ENDIF                                                           
      ENDIF                                                             
! --------------------------------------------------------------------- 
!L Section 15 item 238  Geopotential on pressure levels
      IF (QZ_P) THEN
        CALL V_INT_ZH(P_EXNER_HALF,THETA,Q,PHI_STAR,MODEL_HALF_HEIGHT,
     &                P_FIELD,P_LEVELS,Q_LEVELS)
        DO K=1,Z_P_LEVS
          DO I=1,P_FIELD
            PZ(I)=Z_PRESS(k)*100.0   ! convert to pascals
          ENDDO
          CALL V_INT_Z(PZ,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,
     &    MODEL_HALF_HEIGHT,WORK5,P_FIELD,p_LEVELS,Q_LEVELS,
     &    Z_REF,AKH,BKH,FIRST_FLD_PT,P_FIELD)

! put on u grid
          CALL P_TO_UV(WORK5(FIRST_FLD_PT),Z_P(FIRST_FLD_PT,k)
     &      ,P_FIELD-FIRST_FLD_PT+1,U_FIELD-FIRST_FLD_PT+1
     &      ,ROW_LENGTH,P_LAST_ROW)
        ENDDO

!L Section 15 item 239 U*Z

        IF (QUZ_P) THEN
          DO K=1,UZ_P_LEVS
            DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               UZ_P(I,K) = UCOMP_P(I,UZ_IND(k,1))*Z_P(I,UZ_IND(K,2))
            ENDDO
          ENDDO
        ENDIF 


!L Section 15 item 240 V*Z

        IF (QVZ_P) THEN
          DO K=1,VZ_P_LEVS
            DO I=FIRST_FLD_PT,LAST_U_FLD_PT
               VZ_P(I,K) = VCOMP_P(I,VZ_IND(k,1))*Z_P(I,VZ_IND(K,2))
            ENDDO
          ENDDO
        ENDIF 
      ENDIF

!L Section 15 item 241 mountain torque per unit area
!L
!L  a*orog* E-W pressure gradient =  a* orography * dp /(a dlong)

      IF (Q_MT) THEN
       FACTOR=LONGITUDE_STEP_INVERSE
! wrong at row ends if non MPP code
        DO i=FIRST_FLD_PT,LAST_U_FLD_PT
         M_TORQUE(I)=0.25*factor*( (orog(i)+orog(i+1))*
     &              (pstar(i+1)-pstar(i))
     &           + (orog(i+row_length)+orog(i+row_length+1))*
     &             (pstar(i+row_length+1)-pstar(i+row_length)))
        ENDDO 
      ENDIF
C ---------------------------------------------------------------------

      RETURN
      END
