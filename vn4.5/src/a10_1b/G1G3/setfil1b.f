C ******************************COPYRIGHT******************************
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
C
CLL  SUBROUTINE SET_FIL     -----------------------------------------
CLL
CLL  PURPOSE: CALCULATES THE WAVE-NUMBER ON EACH ROW AT WHICH FOURIER
CLL           FILTERING SHOULD BEGIN. AS FOURIER CHOPPING IS BEING USED
CLL           THIS WAVE-NUMBER ONLY AND NO WEIGHTS ARE RETURNED. THE
CLL           ROW IN EACH HEMISPHERE ON WHICH FILTERING ENDS IS ALSO
CLL           RETURNED.
CLL  N.B.     CODE WILL ALWAYS FILTER 1 ROW IN NORTHERN HEMISPHERE.
CLL           THIS DOES NOT IMPLY THAT SOME WAVE NUMBERS ON THAT ROW
CLL           ARE CHOPPED, IT IS MERELY A WAY OF AVOIDING ALLOCATING A
CLL           ZERO LENGTH ARRAY IN FILTER WHICH CAUSES A FAILURE.
CLL  NOT SUITABLE FOR I.B.M USE.
CLL
CLL  WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.2   25/10/96  New deck for HADCM2-specific section A10_1B,
!LL                   as SETFIL1A but with reintroduced error in
!LL                   computation of SH filter row.       T.Johns
!LL   4.3   14/04/97  Print warning if filtering entire area.  T Johns
!LL                   Correct HADCM2-specific code to use global row.
!LL   4.4   08/08/97  Remove FFT filter data common block,
!LL                   Add code for chunking optimisation
!LL                                                         P.Burton
!LL   4.5   24/02/98  Catch error when whole model is filtered.
!LL                                                 Paul Burton
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL  SYSTEM COMPONENTS COVERED: P141
CLL
CLL  SYSTEM TASK: P1
CLL
CLL  DOCUMENTATION:        EQUATIONS (50) TO (52) IN SECTION 3.5
CLL                        OF UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE SET_FIL
     1                  (U,C,ADJUSTMENT_TIMESTEP,ADVECTION_TIMESTEP,
     2                   SEC_P_LATITUDE,FILTER_WAVE_NUMBER_P_ROWS,
     3                   FILTER_WAVE_NUMBER_U_ROWS,
     4                   LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5                   SOUTHERN_FILTERED_P_ROW,P_FIELD,U_FIELD,
     6                   P_LEVELS,ROW_LENGTH,
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
     7                   FILTERING_SAFETY_FACTOR,TWO_D_GRID_CORRECTION)

      IMPLICIT NONE

! All FLDPTR arguments are intent IN
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
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF MODEL LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS ON A ROW.
     *, NORTHERN_FILTERED_P_ROW !IN LAST ROW, MOVING EQUATORWARDS,
     *                          ! IN NORTHERN HEMISPHERE
     *                          ! ON WHICH FILTERING IS PERFORMED.
     *, SOUTHERN_FILTERED_P_ROW !IN LAST ROW, MOVING EQUATORWARDS,
     *                          ! IN SOUTHERN HEMISPHERE
     *                          ! ON WHICH FILTERING IS PERFORMED.

      INTEGER
     &  FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &      ! LARGEST WAVE NUMBER ON EACH P ROW WHICH IS NOT FILTERED
     &, FILTER_WAVE_NUMBER_U_ROWS(GLOBAL_U_FIELD/GLOBAL_ROW_LENGTH)
     &      ! LARGEST WAVE NUMBER ON EACH U ROW WHICH IS NOT FILTERED


      REAL
     * U(U_FIELD,P_LEVELS)     !IN ZONAL WIND COMPONENT.
     *,SEC_P_LATITUDE(P_FIELD) !IN 1/(COS(LAT)) AT P POINTS.
     *,ADJUSTMENT_TIMESTEP     !IN
     *,ADVECTION_TIMESTEP      !IN
     *,C                       !IN EXTENAL GRAVITY WAVE SPEED.
     *,LONGITUDE_STEP_INVERSE  !IN 1/(DELTA LAMDA)
     *,TWO_D_GRID_CORRECTION(P_FIELD/ROW_LENGTH) !IN TERM REPRESENTING
     *                         ! 2-D EFFECT OF WAVES COMPARED TO 1D
     *                         ! REPRESENTATION USED HERE.
     *,FILTERING_SAFETY_FACTOR !IN TERM REPRESENTING SAFETY MARGIN TO
     *                         ! ADD TO GRID_CORRECTION
C*---------------------------------------------------------------------
! Common blocks and parameters for MPP code
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
!LL 4.2      16/08/96 Removed filt_wave_no_common variable.
!LL                                                        P.Burton
!LL 4.4      08/08/97 Removed sarr_len and rarr_len arrays
!LL                   Added filt_level variable
!LL                   Increased parameters to maximum likely values
!LL                                               P.Burton

! Called by SET_FIL and FILTER - used to communicate the decomposition
! of data in the ffts
! NB : Comdeck PARVARS must be *CALLed before this comdeck.

      INTEGER MAX_ROW_LEN,MAX_ROWS,MAX_LEVELS,MAX_ROWS_TO_FILTER
      PARAMETER(MAX_ROW_LEN=500,MAX_ROWS=500,MAX_LEVELS=60)
      PARAMETER(MAX_ROWS_TO_FILTER=0.5*MAX_ROWS*MAX_LEVELS)

! Common block for communication between SETFILT and FILTER
! We set up these arrays:
! filt_send_map(7,n_items_to_send,fld_type) - contains information about
! the rows of data that this processor has to send off to be filtered
!
! filt_recv_map(7,n_items_to_send,fld_type) - contains information about
! the rows of data that this processor receives to be filtered
!
! filt_info(row_number,fld_type) - contains information about the
! rows of data that this processor will be filtering
!

! Two sets of everything - one for P_FIELDs and one for U_FIELDs.

      INTEGER
     & south_filt_p_row  ! southern filtered p row


      REAL
     & global_trigs(MAX_ROW_LEN) ! global version of TRIGS array

      INTEGER
     & fft_rows(2)   ! total number of rows I will fft

      LOGICAL
     &  filter_off ! set to true if no filtering to be done (usually
!                    indicates an error has occurred

      INTEGER filt_smap_len, filt_rmap_len
      PARAMETER (filt_smap_len = MAX_ROWS_TO_FILTER,
     &           filt_rmap_len = MAX_ROWS_TO_FILTER)
      INTEGER filt_send_map(7,filt_smap_len,2),
     &        filt_recv_map(7,filt_rmap_len,2),
     &        n_items_to_send(2), n_items_to_recv(2),
     &        filt_info(MAX_ROWS_TO_FILTER,2),
     &        filt_level(MAX_ROWS_TO_FILTER,2),
     &        filt_send_start(filt_smap_len,2),
     &        filt_recv_start(filt_smap_len,2),
     &        filt_send_max(filt_smap_len,2),
     &        filt_recv_max(filt_smap_len,2)
      COMMON /PAR_FFT/ south_filt_p_row,
     &                 global_trigs,
     &                 filt_send_map, filt_recv_map,
     &                 n_items_to_send, n_items_to_recv, filt_info,
     &                 filt_level,
     &                 fft_rows, filter_off,
     &                 filt_send_start, filt_recv_start,
     &                 filt_send_max, filt_recv_max

! End COMDECK PARFFTS

CDIR$ FIXED
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GC - General Communication primitives package. For use on
C multiprocessor shared memory and message passing systems.
C
C
C LICENSING TERMS
C
C  GC is provided free of charge. Unless otherwise agreed with SINTEF,
C  use and redistribution in source and binary forms are permitted
C  provided that
C
C      (1) source distributions retain all comments appearing within
C          this file header, and
C
C      (2) distributions including binaries display the following
C          acknowledgement:
C
C              "This product includes software developed by SINTEF.",
C
C          in the documentation or other materials provided with the
C          distribution and in all advertising materials mentioning
C          features or use of this software.
C
C  The name of SINTEF may not be used to endorse or promote products
C  derived from this software without specific prior written
C  permission.  SINTEF disclaims any warranty that this software will
C  be fit for any specific purposes. In no event shall SINTEF be liable
C  for any loss of performance or for indirect or consequential damage
C  or direct or indirect injury of any kind. In no case shall SINTEF
C  be liable for any representation or warranty make to any third party
C  by the users of this software.
C
C
C Fortran header file. PLEASE use the parameter variables in user
C routines calling GC and NOT the numeric values. The latter are
C subject to change without further notice.
C
C---------------------------------------------- ------------------------
C $Id: gpb2f402,v 1.10 1996/11/28 20:36:24 t11pb Exp $
C (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

C    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
C                    value to be set from a shell variable.
C                      Author: Bob Carruthers  Cray Research.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     GC general options
      INTEGER GC_OK, GC_FAIL, GC_NONE, GC_ANY, GC_DONTCARE,
     $     GC_SHM_DIR, GC_SHM_GET, GC_SHM_PUT, GC_USE_GET, GC_USE_PUT
     &   , GC_NAM_TIMEOUT, GC_SHM_SAFE
      PARAMETER (GC_OK         =     0)
      PARAMETER (GC_FAIL       =    -1)
      PARAMETER (GC_NONE       =     0)
      PARAMETER (GC_ANY        =    -1)
      PARAMETER (GC_DONTCARE   =    -1)
      PARAMETER (GC_SHM_DIR    =     1)
      PARAMETER (GC_SHM_SAFE   =     2)
      PARAMETER (GC_NAM_TIMEOUT=     4)
      PARAMETER (GC_SHM_GET    = -9999)
      PARAMETER (GC_SHM_PUT    = -9998)
      PARAMETER (GC_USE_GET    = -9999)
      PARAMETER (GC_USE_PUT    = -9998)

C     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

C     GC groups (GCG) support
      INTEGER GC_ALLGROUP, GCG_ALL
      PARAMETER (GC_ALLGROUP = 0)
      PARAMETER (GCG_ALL = GC_ALLGROUP)

C     GC groups (GCG) functions
      INTEGER GCG_ME

C     GC reserved message tags
      INTEGER GC_MTAG_LOW, GC_MTAG_HIGH
      PARAMETER (GC_MTAG_LOW   = 999999901)
      PARAMETER (GC_MTAG_HIGH  = 999999999)

C     GCG_RALLETOALLE index parameters
      INTEGER S_DESTINATION_PE, S_BASE_ADDRESS_IN_SEND_ARRAY,
     $     S_NUMBER_OF_ELEMENTS_IN_ITEM, S_STRIDE_IN_SEND_ARRAY,
     $     S_ELEMENT_LENGTH, S_BASE_ADDRESS_IN_RECV_ARRAY,
     $     S_STRIDE_IN_RECV_ARRAY
      PARAMETER (S_DESTINATION_PE = 1)
      PARAMETER (S_BASE_ADDRESS_IN_SEND_ARRAY = 2)
      PARAMETER (S_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (S_STRIDE_IN_SEND_ARRAY = 4)
      PARAMETER (S_ELEMENT_LENGTH = 5)
      PARAMETER (S_BASE_ADDRESS_IN_RECV_ARRAY = 6)
      PARAMETER (S_STRIDE_IN_RECV_ARRAY = 7)

      INTEGER R_SOURCE_PE, R_BASE_ADDRESS_IN_RECV_ARRAY,
     $     R_NUMBER_OF_ELEMENTS_IN_ITEM, R_STRIDE_IN_RECV_ARRAY,
     $     R_ELEMENT_LENGTH, R_BASE_ADDRESS_IN_SEND_ARRAY,
     $     R_STRIDE_IN_SEND_ARRAY
      PARAMETER (R_SOURCE_PE = 1)
      PARAMETER (R_BASE_ADDRESS_IN_RECV_ARRAY = 2)
      PARAMETER (R_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (R_STRIDE_IN_RECV_ARRAY = 4)
      PARAMETER (R_ELEMENT_LENGTH = 5)
      PARAMETER (R_BASE_ADDRESS_IN_SEND_ARRAY = 6)
      PARAMETER (R_STRIDE_IN_SEND_ARRAY = 7)

C*L  DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C   DEFINE LOCAL ARRAYS: NONE ARE REQUIRED

C*---------------------------------------------------------------------
C   DEFINE LOCAL VARIABLES

      REAL
     *  COURANT         !HOLDS MAXIMUM COURANT NUMBER.
     *, ADVECTION_COURANT ! HOLDS ADVECTION COURANT NUMBER.
     *, ADJUSTMENT_COURANT ! HOLDS ADJUSTMENT COURANT NUMBER.
     *, C_TERM          ! HOLDS TERM WHICH MULTIPLIES EPSI_J IN
     *                  ! EQUATION 50.
     *, C_CONST         ! HOLDS CONSTANT PART OF C_TERM

      REAL
     *  SCALAR1
     *, SCALAR2
     *, ALPHA           ! HOLDS DELTA LAMDA / COS(LAT)
      REAL
     &  U_MAX_INITIAL(P_FIELD/ROW_LENGTH)
! Holds maximum windspeed along each U row
     &,    U_MAX(P_FIELD/ROW_LENGTH)  ! Holds maximum windspeed
                                      ! along each P row
C TWO DIFFERENT NAMES NEEDED WHERE ONE WOULD SUFFICE TO STRICTLY CONFORM
C WITH CRAY'S CASE CONSTRUCT. SEE ALSO CUT_OFF_N AND CUT_OFF_S.

C   COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,K
     &, CUT_OFF         ! HOLDS SMALLEST WAVE NUMBER FILTERED.
     *, P_ROWS          ! HOLDS NUMBER OF ROWS OF DATA ON P_GRID
     *, HALF_P_ROWS     ! HOLDS HALF NUMBER OF ROWS OF DATA ON P_GRID
     *, ROW             ! LOOP COUNTER WHICH IS THE ROW ON WHICH
     *                  ! CALCULATION IS TAKING PLACE.
      INTEGER
     &  ROW_START,ROW_END     ! Start and end points for loop over rows
     &, global_row            ! global row number of local row
     &, local_row             ! local row number of global row
     &, glob_filterable_rows  ! Global number of rows for filtering
     &, approx_ffts_per_proc  ! Rough estimate of number of ffts/proc
     &, remainder_ffts        ! "Spare" ffts to be added on
     &, filt_row_number       ! filterable row number
     &, filt_row_index        ! index for a (row,level) to be fft'd
     &, change_pt             ! used in calculation of decomposition
     &, proc_row              ! Current processor row
     &, proc_dest             ! Processor to send a row to
     &, n_send                ! No. of fft data blocks I will send
     &, n_recv                ! No. of fft data blocks I will receive
     &, fld_type              ! indicator of field type (P or U)
     &, north_filt_row        ! local vn for fld_type
     &, south_filt_row        ! local vn for fld_type
     &, n_glob_rows           ! number of global rows for fld_type
     &, oproc_dest            ! destination of last previous row
     &, sourcepe              ! processor this row comes from
     &, row_number            ! row number of this row
     &, fft_row_len           ! row length of fft data (glsize(1)+2)
     &, J                     ! loop counter
     &, info                  ! return code for comms


      LOGICAL
     &  sending               ! Am I sending this row?

      LOGICAL called_before
      DATA called_before/.FALSE./
      INTEGER old_northern_filt,old_southern_filt
      SAVE called_before,old_northern_filt,old_southern_filt



C*L  EXTERNAL SUBROUTINE CALLS:------------------------------------
C NO EXTERNAL SUBROUTINE CALLS
C*---------------------------------------------------------------------
CL CALL COMDECK TO GET CONSTANTS USED.

CLL COMDECK C_SETFIL HOLDS CONSTANTS FOR ROUTINE SET_FIL
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C END OF COMDECK C_SETFIL

CL  MAXIMUM VECTOR LENGTH ASSUMED IS ROW_LENGTH
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL


      IF (.NOT. called_before) THEN
        called_before=.TRUE.
        old_northern_filt=-1
        old_southern_filt=-1
      ENDIF


      C_CONST = C* ADJUSTMENT_TIMESTEP * LONGITUDE_STEP_INVERSE / A
      C_CONST = C_CONST * C_CONST
      ALPHA=1.0/LONGITUDE_STEP_INVERSE
      P_ROWS=P_FIELD/ROW_LENGTH
      HALF_P_ROWS=glsize(2)/2    ! half number of global rows

! Caculate the local maximum U along each row

      DO ROW=1,Offy
        U_MAX_INITIAL(ROW)=0.0
        U_MAX_INITIAL(P_ROWS+1-ROW)=0.0
      ENDDO
      DO ROW=1+Offy,P_ROWS-Offy  ! loop over local rows, missing halos
        SCALAR1=0.0
        DO K=1,P_LEVELS
          DO I=(ROW-1)*ROW_LENGTH+1+Offx,ROW*ROW_LENGTH-Offx
            SCALAR1=MAX(ABS(U(I,K)),SCALAR1)
          ENDDO
        ENDDO  ! loop over levels
        U_MAX_INITIAL(ROW)=SCALAR1
      ENDDO  ! loop over rows

! We only have the local maximum so far. We must now do a maximum
! along each global row.

       CALL GCG_RMAX(P_ROWS-2*Offy,gc_proc_row_group,info,
     &               U_MAX_INITIAL(1+Offy))

! We need to know the values of the rows above and below us, so
! fill in the halos with this information
      CALL SWAPBOUNDS(U_MAX_INITIAL,1,P_ROWS,0,Offy,1)

! Now calculate the maximum "wind" along each P_ROW, by taking the
! maximum wind of the U_ROW either side of it:
! MAX_WIND(P_ROW)=MAX(U_MAX(P_ROW),U_MAX(P_ROW-1))

      IF (attop) THEN
        ROW_START=Offy+2
      ELSE
        ROW_START=Offy+1
      ENDIF
      IF (atbase) THEN
        ROW_END= P_ROWS-Offy-1
      ELSE
        ROW_END=P_ROWS-Offy
      ENDIF

      DO ROW=ROW_START,ROW_END  ! loop over local rows
! Deliberately reintroduce error in computation of SH filter row
! Must use the global rather than local row number to decide if MPP 
        IF (ROW+datastart(2)-Offy-1 .GE. HALF_P_ROWS+2) THEN
! This is the error; 0.0 implies only one row is considered in SH
          U_MAX(ROW)=MAX(U_MAX_INITIAL(ROW),0.0)
        ELSE
! This is the correct code; two adjacent rows are considered in NH
          U_MAX(ROW)=
     &      MAX(U_MAX_INITIAL(ROW),U_MAX_INITIAL(ROW-1))
        ENDIF
! and scale U_MAX
        U_MAX(ROW)=U_MAX(ROW)*ADVECTION_TIMESTEP*
     &    LONGITUDE_STEP_INVERSE*SEC_P_LATITUDE(ROW*ROW_LENGTH)/A
     &    *(TWO_D_GRID_CORRECTION(ROW)+FILTERING_SAFETY_FACTOR)
! and square it
        U_MAX(ROW)=U_MAX(ROW)*U_MAX(ROW)
      ENDDO

! Now loop round each row, setting the filter_wave_number, and
! calculating NORTHERN/SOUTHERN_FILTERED_P_ROW.
! These values are an initial guess, and will have to be MAX/MINed in
! the North/South direction in order to get the correct global values.

      NORTHERN_FILTERED_P_ROW=2
      SOUTHERN_FILTERED_P_ROW=glsize(2)
      CUT_OFF=glsize(1)/2

      DO ROW=ROW_START,ROW_END  ! loop over local rows
        global_row=ROW+datastart(2)-Offy-1  ! global row number
        C_TERM=SEC_P_LATITUDE(ROW*ROW_LENGTH)*
     &    SEC_P_LATITUDE(ROW*ROW_LENGTH)*C_CONST*
     &    TWO_D_GRID_CORRECTION(ROW)*
     &    TWO_D_GRID_CORRECTION(ROW)

        DO I=2,CUT_OFF
          SCALAR1= SIN(ALPHA*(I+1))  ! EPSI_V
          SCALAR2 = SIN(I*ALPHA*.5)  ! EPSI_J
! First term in eqn 50 - Advection Courant Number
          ADVECTION_COURANT = U_MAX(ROW)*SCALAR1*SCALAR1
! Second term in eqn 50 - Adjustment Courant Number
          ADJUSTMENT_COURANT = C_TERM*SCALAR2*SCALAR2
          COURANT = MAX(ADJUSTMENT_COURANT,ADVECTION_COURANT)
          IF (COURANT .GT. 1.0) THEN  ! violated stability criteria

! If Northern hemisphere:
            IF (global_row .LE. HALF_P_ROWS+1)
     &        NORTHERN_FILTERED_P_ROW=global_row
! If Southern hemisphere:
            IF ((global_row .GE. HALF_P_ROWS+2) .AND.
     &          (SOUTHERN_FILTERED_P_ROW .EQ. glsize(2)))
     &        SOUTHERN_FILTERED_P_ROW=global_row

            FILTER_WAVE_NUMBER_P_ROWS(global_row) = I-1

            GOTO 100
          ENDIF ! if stability criteria violated
        ENDDO ! loop over wave-numbers
! If stability criteria not violated, set FILTER_WAVE_NUMBER
! to CUT_OFF
        FILTER_WAVE_NUMBER_P_ROWS(global_row) = CUT_OFF

 100    CONTINUE
      ENDDO ! loop over local rows

! Work out the correct global values for FILTERED_P_ROWs

      CALL GC_IMAX(1,nproc,info,NORTHERN_FILTERED_P_ROW)
      CALL GC_IMIN(1,nproc,info,SOUTHERN_FILTERED_P_ROW)

! Keep a copy of SOUTHERN_FILTERED_P_ROW - we'll use this to decide
! if the field we're filtering is a U_FIELD or P_FIELD
      south_filt_p_row=SOUTHERN_FILTERED_P_ROW

! We have to enforce the condition that the FILTER_WAVE_NUMBER_P_ROWS
! cannot increase as we move towards the pole.
! This is easiest if all the data is available on each processor

! First, broadcast my local part of the filt_wave_number array

      DO I=0,nproc-1,nproc_x
        CALL GC_IBCAST(I,g_blsizep(2,I),I,nproc,info,
     &                 FILTER_WAVE_NUMBER_P_ROWS(g_datastart(2,I)))
      ENDDO

! Northern hemisphere first:
      DO ROW=HALF_P_ROWS,2,-1
        FILTER_WAVE_NUMBER_P_ROWS(ROW) =
     &    MIN(FILTER_WAVE_NUMBER_P_ROWS(ROW),
     &        FILTER_WAVE_NUMBER_P_ROWS(ROW+1))
      ENDDO

! And similarly for the Southern hemisphere
      DO ROW=HALF_P_ROWS+2,glsize(2)-1
        FILTER_WAVE_NUMBER_P_ROWS(ROW) =
     &    MIN(FILTER_WAVE_NUMBER_P_ROWS(ROW),
     &        FILTER_WAVE_NUMBER_P_ROWS(ROW-1))
      ENDDO

! And set the values at the poles:
      FILTER_WAVE_NUMBER_P_ROWS(1) = 0
      FILTER_WAVE_NUMBER_P_ROWS(glsize(2)) = 0

! Make the U_FIELD version
      DO I=1,HALF_P_ROWS
        FILTER_WAVE_NUMBER_U_ROWS(I) = FILTER_WAVE_NUMBER_P_ROWS(I)
        FILTER_WAVE_NUMBER_U_ROWS(glsize(2)-I) =
     &    FILTER_WAVE_NUMBER_P_ROWS(glsize(2)+1-I)
      ENDDO

      IF ((NORTHERN_FILTERED_P_ROW .EQ. HALF_P_ROWS+1) .AND.
     &    (SOUTHERN_FILTERED_P_ROW .EQ. HALF_P_ROWS+2)) THEN
        WRITE(6,*) '*** SETFIL : WARNING *** ',
     &             'Filtering entire model domain.'
        WRITE(6,*) 'Timestep should be reduced.'
        WRITE(6,*) 'Filtering has been switched off.'
        filter_off=.TRUE.
        GOTO 9999
      ENDIF

! --------------------------------------------------------------------

      IF ((old_northern_filt .NE. NORTHERN_FILTERED_P_ROW) .OR.
     &    (old_southern_filt .NE. SOUTHERN_FILTERED_P_ROW)) THEN
        old_northern_filt = NORTHERN_FILTERED_P_ROW
        old_southern_filt = SOUTHERN_FILTERED_P_ROW


!     WRITE(6,*) 'SETFIL : Updating filtering area to ',
!    &           '2 - > ',NORTHERN_FILTERED_P_ROW,' and ',
!    &           SOUTHERN_FILTERED_P_ROW,' -> ',glsize(2)-1
!     WRITE(6,*) 'SETFIL: Filtering ',
!    &  REAL((glsize(2)-2)-
!    &  (SOUTHERN_FILTERED_P_ROW-NORTHERN_FILTERED_P_ROW-1))/
!    &  REAL(glsize(2)-2),' % of grid (not counting polar rows)'
      fft_row_len=glsize(1)+2
      DO fld_type=1,2  ! 1=P_FIELD and 2=U_FIELD

        north_filt_row=NORTHERN_FILTERED_P_ROW
        IF (fld_type .EQ. 1) THEN
          south_filt_row=SOUTHERN_FILTERED_P_ROW
          n_glob_rows=glsize(2)
        ELSE
          south_filt_row=SOUTHERN_FILTERED_P_ROW-1
          n_glob_rows=glsize(2)-1
        ENDIF

! Set up all the information required for decomposing the ffts
! The information is passed to the filter routine via the
! PAR_FFT common block

        glob_filterable_rows=n_glob_rows-(south_filt_row-
     &                       north_filt_row-1)-2
! The "-2" arises because we don't want to filter the polar
! rows.

        approx_ffts_per_proc=(glob_filterable_rows*P_LEVELS)/nproc
        remainder_ffts=(glob_filterable_rows*P_LEVELS)-
     &                  approx_ffts_per_proc*nproc

        filter_off=.FALSE.
        IF (REAL(glob_filterable_rows*P_LEVELS)/REAL(nproc) .GT.
     &    MAX_ROWS_TO_FILTER) THEN
          WRITE(6,*) 'Error : Too many rows to filter'
          WRITE(6,*) 'Increase MAX_ROWS_TO_FILTER in comdeck ',
     &                       'PARFFTS.'
          filter_off=.TRUE.
          GOTO 9999
        ENDIF
! approx_ffts_per_proc is the number of rows each processor would have
! to fft, if the total number of fft rows divided exactly between
! the total number of processors.
! remainder_ffts is the number of fft rows left unassigned, if each
! processor takes on approx_ffts_per_proc rows to do.

! The rows are distributed like this:
! The first "remainder_ffts" processors are given
! approx_ffts_per_proc+1 rows each to fft
! All the other processors are given approx_ffts_per_proc rows each.
! The indexing of the filterable rows is like this:
! index  global_row_number  level
! 1      2                  1
! 2      2                  2
! :      2                  :
! :      2                  maxlevel
! :      3                  1
! etc.

        n_send=0      ! number of rows I will send
        n_recv=0   ! number of processor rows I will receive from
        fft_rows(fld_type)=0    ! total number of rows I will fft
        proc_dest = -1
        change_pt=(approx_ffts_per_proc+1)*remainder_ffts


        DO ROW=2,n_glob_rows-1  ! loop over global rows, missing poles

          IF ((ROW .LE. north_filt_row) .OR.
     &        (ROW .GE. south_filt_row)) THEN
! This row needs to be filtered

            IF (ROW .LE. north_filt_row) THEN
            ! Northern hemisphere
              filt_row_number=ROW-1
            ELSE
            ! Southern hemisphere
              filt_row_number=ROW-
     &                       (south_filt_row-
     &                        north_filt_row-1)-1
            ENDIF

! Is this row of data on my processor?
            sending=.FALSE.
            IF ((ROW .GE. datastart(2)) .AND.
     &          (ROW .LE. datastart(2)+blsizep(2)-1)) sending=.TRUE.

! Which processor row does it live on?
            DO I=0,nproc_y-1  ! loop over processor row
              local_row=ROW-g_datastart(2,I*nproc_x)+1
              IF ((local_row .LE. g_blsizep(2,I*nproc_x)) .AND.
     &            (local_row .GT. 0)) proc_row=I
            ENDDO

            DO K=1,P_LEVELS
              filt_row_index=(K+(filt_row_number-1)*P_LEVELS)

! Which processor is it going to?
              oproc_dest = proc_dest
              IF (filt_row_index .LE. change_pt) THEN
! It's a processor with approx_ffts_per_proc+1 rows to do
                proc_dest=(filt_row_index-1)/(approx_ffts_per_proc+1)
              ELSE
! It's a processor with approx_ffts_per_proc rows to do
                proc_dest=remainder_ffts+
     &            (filt_row_index-1-change_pt)/approx_ffts_per_proc
              ENDIF

! Am I sending it?
              IF (sending) THEN
                 if (proc_dest .eq. oproc_dest .and. k .ge. 2) then
                    filt_send_max(n_send,fld_type) =
     &                   filt_send_max(n_send,fld_type) + 1
                 else
                    n_send = n_send + 1
                    filt_send_start(n_send,fld_type) = K
                    filt_send_max(n_send,fld_type) = 1
                    filt_send_map(S_DESTINATION_PE,n_send,fld_type) =
     &                proc_dest

                    filt_send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,
     &                            n_send,fld_type) =
     &                (K-1)*p_field +
     &                (ROW-datastart(2)+offy)*lasize(1) + offx + 1

                    filt_send_map(S_STRIDE_IN_SEND_ARRAY,
     &                            n_send,fld_type) =
     &                p_field

                    filt_send_map(S_ELEMENT_LENGTH,n_send,fld_type) =
     &                blsizep(1)
! Calculate the number of this fft within the receiving processor
! We do this by calculating the first fft the receiving processor does
! and taking this number away from the filt_row_index
                    IF (proc_dest .LE. remainder_ffts-1) THEN

                       filt_send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,
     &                               n_send,fld_type)=
     &                      (filt_row_index
     &                      +1-(proc_dest*(approx_ffts_per_proc+1)+1)-1)
     &                      *fft_row_len + datastart(1)
                    ELSE

                       filt_send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,
     &                               n_send,fld_type)=
     &                      (filt_row_index+
     &                      1-((remainder_ffts*(approx_ffts_per_proc+1))
     &                      +(proc_dest-remainder_ffts)
     &                      *approx_ffts_per_proc+1)
     &                      -1)*fft_row_len + datastart(1)
                    ENDIF

                    filt_send_map(S_STRIDE_IN_RECV_ARRAY,
     &                            n_send,fld_type) =
     &                fft_row_len

                 endif
              ENDIF ! IF sending

! Am I receiving it?
              IF (proc_dest .EQ. mype) THEN  ! Is this my row to fft?
                fft_rows(fld_type)=fft_rows(fld_type)+1
                filt_info(fft_rows(fld_type),fld_type) = row

                if (proc_dest .eq. oproc_dest .and. k .ge. 2) then
                   do j = 0,nproc_x-1
                      filt_recv_max(n_recv-j,fld_type) =
     &                     filt_recv_max(n_recv-j,fld_type) + 1
                   enddo
                else
                   DO J=0,nproc_x-1 ! loop over processor row

                      sourcepe=proc_row*nproc_x+J
                      row_number=
     &                     ROW-g_datastart(2,proc_row*nproc_x)+Offy+1
! The following loop could be replaced with a strided receive

                      n_recv = n_recv + 1
                      filt_recv_start(n_recv,fld_type) = K
                      filt_recv_max(n_recv,fld_type) = 1
                      filt_recv_map(R_SOURCE_PE,n_recv,fld_type) =
     &                  sourcepe

                      filt_recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,
     &                              n_recv,fld_type) =
     &                  (fft_rows(fld_type)-1)*fft_row_len +
     &                  g_datastart(1,sourcepe)

                      filt_recv_map(R_STRIDE_IN_RECV_ARRAY,
     &                              n_recv,fld_type) =
     &                  fft_row_len

                      filt_recv_map(R_ELEMENT_LENGTH,
     &                              n_recv,fld_type) =
     &                  g_blsizep(1,sourcepe)

                      filt_recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,
     &                              n_recv,fld_type) =
     &                  (K-1)*g_lasize(1,sourcepe)*g_lasize(2,sourcepe)
     &                  + 1 + offx + (row_number-1)*
     &                  g_lasize(1,sourcepe)

                      filt_recv_map(R_STRIDE_IN_SEND_ARRAY,
     &                              n_recv,fld_type) =
     &                  g_lasize(1,sourcepe)*g_lasize(2,sourcepe)

                   ENDDO        ! J
                endif

              ENDIF  ! IF my row to fft

            ENDDO  ! loop over levels (K)
          ENDIF  ! IF this row needs to be filtered
        ENDDO  ! loop over rows (ROW)

      n_items_to_send(fld_type) = n_send                               
      n_items_to_recv(fld_type) = n_recv                               

      ENDDO  ! loop over field types (fld_type)
      ENDIF  ! IF we need to update fft decomposition

!
!      WRITE(6,*) 'Processor ',mype
!      WRITE(6,*) 'Sending ',n_items_to_send,' chunks.'
!      WRITE(6,*) 'Receiving ',n_items_to_recv,' chunks.'
!       WRITE(6,*) 'Doing ',fft_rows,' ffts.'
!
!       write(6,*) 'filter wave nos:'
!       do i=1,glsize(2)
!          write(6,*) i,filt_wave_number(i)
!       enddo

 9999   CONTINUE


! END OF ROUTINE SET_FIL

      RETURN
      END
