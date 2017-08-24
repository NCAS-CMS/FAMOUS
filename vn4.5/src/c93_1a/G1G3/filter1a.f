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
CLL  SUBROUTINE FILTER     -----------------------------------------
CLL
CLL  PURPOSE: FOURIER DAMPS ALL WAVES AFTER THE WAVE-NUMBER HELD IN
CLL           FILTER_WAVE_NUMBER FOR ALL ROWS WHERE FILTERING IS
CLL           NECESSARY. USES THE ECMWF FFT ROUTINES TO CALCULATE THE
CLL           FOURIER PARTS OF THE CODE.
CLL  NOT SUITABLE FOR I.B.M USE.
CLL
CLL  WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
!   4.1      09/05/96 Added MPP code.   P.Burton
!    4.2    18/10/96  New name for group of processors in all_to_all
!                     P.Burton
!LL 4.2      16/08/96 Added TYPLDPT variables.
!LL                   Made FILTER_WAVE_NUMBER globally sized.
!LL                   Removal of identical deck FILTER2A, requires
!LL                   *DEFs to be changed.                P.Burton
!LL 4.2      17/10/96 Changed arguments to ALLTOALL for GCC v1.1
!LL                   P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
!LL 4.4      08/08/97 Remove filter data common block, adding new
!LL                   subroutine MPP_FILTER.
!LL                   Add chunking optimisation.
!LL                   Add T3E optimised communications.
!LL                                                       P.Burton
!LL 4.5      22/06/98 Remove redundant divides.
!LL                     Author: Bob Carruthers, Cray Research
CLL
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                        STANDARD B. VERSION 2, DATED 18/01/90
CLL  SYSTEM COMPONENTS COVERED: 142
CLL  SYSTEM TASK: P1
CLL  DOCUMENTATION:        EQUATIONS (53) AND (54) IN SECTION 3.5
CLL                        OF UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 11, DATED 10/10/90.
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE FILTER
     1                 (FIELD,FIELD_LENGTH,LEVELS,FILTER_SPACE,
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
     &                  FILTER_WAVE_NUMBER,TRIGS,IFAX,
     3                  NORTHERN_FILTERED_ROW,
     4                  SOUTHERN_FILTERED_ROW)

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

      INTEGER
     *  FIELD_LENGTH       !IN HORIZONTAL DIMENSION OF FIELD TO BE
     *                     !   FILTERED.
     *, LEVELS             !IN NUMBER OF MODEL LEVELS IN FIELD.
     *, ROW_LENGTH         !IN NUMBER OF POINTS ON A ROW.
     *, NORTHERN_FILTERED_ROW !IN LAST ROW, MOVING EQUATORWARDS,
     *                          ! IN NORTHERN HEMISPHERE
     *                          ! ON WHICH FILTERING IS PERFORMED.
     *, SOUTHERN_FILTERED_ROW !IN LAST ROW, MOVING EQUATORWARDS,
     *                          ! IN SOUTHERN HEMISPHERE
     *                          ! ON WHICH FILTERING IS PERFORMED.
     *, FILTER_SPACE       !IN HORIZONTAL DIMENSION OF ARRAY NEEDED TO
     *                     ! HOLD DATA TO BE FILTERED TO PASS TO FFT'S.
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED IN FFT'S
     *, FILTER_WAVE_NUMBER(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &     ! LAST WAVE NUMBER ON EACH ROW WHICH IS NOT FILTERED

      REAL
     * FIELD(FIELD_LENGTH,LEVELS) !INOUT HOLDS FIELD TO BE FILTERED.

      REAL
     * TRIGS(ROW_LENGTH)    !IN HOLDS TRIGONOMETRIC TERMS USED IN FFT'S
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
C   DEFINE LOCAL ARRAYS: 1 IS REQUIRED

C*---------------------------------------------------------------------
C   DEFINE LOCAL VARIABLES

C   COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J ,IJ
     &, LOT                 ! NUMBER OF DATA VECTORS PASSED TO FFT'S.
     &, JUMP                ! NUMBER OF STORAGE LOCATIONS BETWEEN THE
     &                      ! START OF CONSECUTIVE DATA VECTOR.
     &, INCREMENT           ! NUMBER OF STORAGE LOCATIONS BETWEEN EACH
     &                      ! ELEMENT OF THE SAME DATA VECTOR, 1, IF
     &                      ! CONSECUTIVE.
     &, FFTSIGN             ! PARAMETER DETERMINING WHETHER SPECTRAL
     &                      ! TO GRID-POINT (FFTSIGN = 1) OR GRID-POINT
     &                      ! TO SPECTRAL (FFTSIGN = -1) FFT'S ARE
     &                      ! REQUIRED.
     &, fld_type,fft_row_len
     &, info,max_field_length

C*L  EXTERNAL SUBROUTINE CALLS:------------------------------------
      EXTERNAL FOURIER
C*---------------------------------------------------------------------


      IF (filter_off) THEN  ! if filter has been switched off for
!                             some reason
        WRITE(6,*) 'FILTERING SWITCHED OFF.'
        WRITE(6,*) 'See earlier in output for error message'
        GOTO 9999
      ENDIF


! Determine what field type this is, and set variables accordingly
      fld_type=south_filt_p_row-SOUTHERN_FILTERED_ROW+1
!             =1 for p_field and =2 for u_field

! Find the number of levels involved in the data transpose

      DO I = 1, n_items_to_send(fld_type)
        filt_send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,i,fld_type) =
     &    MAX( MIN(filt_send_max(i,fld_type),
     &             LEVELS - filt_send_start(i,fld_type) + 1),
     &         0)
      ENDDO
      DO I = 1, n_items_to_recv(fld_type)
        filt_recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,i,fld_type) =
     &    MAX( MIN(filt_recv_max(i,fld_type),
     &             LEVELS - filt_recv_start(i,fld_type) + 1),
     &         0)
      ENDDO


      fft_row_len=glsize(1)+2
! Call routine to perform redistribution of data and filtering

! Find the maximum field size across processors

      max_field_length=FIELD_LENGTH*LEVELS
      CALL GC_IMAX(1,nproc,info,max_field_length)
      
      CALL MPP_FILTER(FIELD,FIELD_LENGTH,LEVELS,max_field_length,
     &                fft_rows(fld_type),fft_row_len,fld_type,
     &                GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH,
     &                FILTER_WAVE_NUMBER,IFAX)
 9999 CONTINUE

CL END OF ROUTINE FILTER

      RETURN
      END


      SUBROUTINE MPP_FILTER
     &                     (FIELD,FIELD_LENGTH,LEVELS,dummy_len,
     &                      rows_to_fft,fft_row_len,fld_type,
     &                      global_rows,
     &                      filter_wave_number,ifax)

      IMPLICIT NONE

      INTEGER
     &  FIELD_LENGTH  ! IN : horizontal size of FIELD array
     &, LEVELS        ! IN : number of levels in FIELD
     &, dummy_len     ! IN : size for dummy_FIELD
     &, rows_to_fft   ! IN : number of rows I will fft
     &, fft_row_len   ! IN : size of rows that I will fft
     &, fld_type      ! IN : field type (P or U) of FIELD
     &, global_rows   ! IN : number of rows in global field
     &, filter_wave_number(global_rows)
!                       IN : last wave number on each global row
!                            which is not filtered
     &, ifax(10)      ! IN : factors of row length used in ffts

      REAL
     &  FIELD(FIELD_LENGTH,LEVELS)  ! IN/OUT : data to be filtered


! Dynamic array for putting data to be filtered

      REAL
     &  FILTER_DATA(fft_row_len*rows_to_fft)
cdir$ cache_align filter_data

! Comdecks/ commonblocks/parameters
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


! The following variables are used to make the MPP filter more
! efficient. The FILTER_DATA array may contain rows of data for
! levels LEVELS+1 -> P_LEVELS containing null data, but which
! was filtered. Now, the chunk_start(i) array points to chunks of
! rows (being chunk_row(i) rows in length) which require filtering
! n_chunks says how many chunks there are for filtering

      INTEGER
     &  chunk_start(rows_to_fft)
     &, chunk_rows(rows_to_fft)
     &, n_chunks,chunk

! Local scalars

      INTEGER
     &  increment,jump,fftsign,lot,row_number,filt_wave_no
     &, i,j

      INTEGER
     &  info,flag

! 1.0 Set up chunk arrays describing where the filterable data
!     in filter_data will be.

      n_chunks=1
      chunk_start(1)=-1
      chunk_rows(1)=0

      DO i=1,rows_to_fft

        IF ((filt_level(i,fld_type) .GT. LEVELS) .AND.
     &      (chunk_rows(n_chunks) .NE. 0)) THEN
          n_chunks=n_chunks+1
          chunk_start(n_chunks)=-1
          chunk_rows(n_chunks)=0
        ENDIF

        IF (filt_level(i,fld_type) .LE. LEVELS) THEN
          IF (chunk_start(n_chunks) .EQ. -1)
     &      chunk_start(n_chunks)=i
          chunk_rows(n_chunks)=chunk_rows(n_chunks)+1
        ENDIF

      ENDDO

      IF (chunk_start(n_chunks) .EQ. -1) THEN
        n_chunks=n_chunks-1
      ENDIF

! 2.0 Move the data to the filter_data array


      flag=GC_NONE

      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)
      info=GC_NONE

      CALL gcg_ralltoalle(field, filt_send_map(1,1,fld_type),
     &                    n_items_to_send(fld_type),
     &                    FIELD_LENGTH*LEVELS,
     &                    filter_data, filt_recv_map(1,1,fld_type),
     &                    n_items_to_recv(fld_type),
     &                    fft_row_len*rows_to_fft,
     &                    gc_all_proc_group, flag, info)


! 3.0 Move to wave space

      INCREMENT=1
      JUMP=fft_row_len
      FFTSIGN=-1

      DO chunk=1,n_chunks

        LOT=chunk_rows(chunk)

        CALL FOURIER(FILTER_DATA(1+(chunk_start(chunk)-1)*fft_row_len),
     &               global_trigs,
     &               IFAX,INCREMENT,JUMP,glsize(1),LOT,FFTSIGN)
      ENDDO

! 4.0 Perform the truncation

      DO chunk=1,n_chunks
        DO I=chunk_start(chunk),chunk_start(chunk)+chunk_rows(chunk)-1

          row_number=filt_info(I,fld_type)
          filt_wave_no=FILTER_WAVE_NUMBER(row_number)

          DO J=3+(filt_wave_no)*2,fft_row_len

            FILTER_DATA((I-1)*fft_row_len+J)=
     &                       FILTER_DATA((I-1)*fft_row_len+J)*
     &                       (filt_wave_no/REAL(((J-1)/2)*2))**2

          ENDDO ! J : loop along fft row
        ENDDO ! I : loop over fft rows in chunk
      ENDDO ! chunk : loop over chunks of rows

! 5.0 Move back to grid-point space

      FFTSIGN=1

      DO chunk=1,n_chunks

        LOT=chunk_rows(chunk)

        CALL FOURIER(FILTER_DATA(1+(chunk_start(chunk)-1)*fft_row_len),
     &               global_trigs,
     &               IFAX,INCREMENT,JUMP,glsize(1),LOT,FFTSIGN)
      ENDDO

! 6.0 And finally, move the data back to the FIELD array


      flag=GC_NONE

      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_GET,info)
      info=GC_NONE
      CALL gcg_ralltoalle(filter_data, filt_recv_map(1,1,fld_type),
     &                   n_items_to_recv(fld_type),                    
     &                    fft_row_len*rows_to_fft,
     &                    field, filt_send_map(1,1,fld_type),
     &                    n_items_to_send(fld_type),
     &                    FIELD_LENGTH*LEVELS,
     &                    gc_all_proc_group, flag, info)


      RETURN
      END
