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
CLL   SUBROUTINE SET_TRAC ---------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES NUMBER OF EAST-WEST SWEEPS OF HORIZONTAL
CLL              ADVECTION REQUIRED ON EACH ROW TO MAINTAIN A CFL
CLL              NUMBER LESS THAN 0.5.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP, CRAY T3E, and Workstations.     
CLL
CLL   WRITTEN BY M.H. MAWSON
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed (was used to
CLL                   switch on ISAMAX,ISAMIN - for PVP systems only).
CLL                                   S.J.Swarbrick
CLL  4.2  15/08/96  MPP code added. Loop structure modified. RTHBarnes.
!LL  4.3  17/03/97  Remove print statement.  RTHBarnes.
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P123
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION: U.M. Doc. Paper 11, by M.J.P. Cullen
CLL
CLLEND-----------------------------------------------------------------

C
C*L   ARGUMENTS:-------------------------------------------------------

      SUBROUTINE SET_TRAC
     &                   (TRACER_EW_SWEEPS,U,P_FIELD,U_FIELD,
     &                    P_LEVELS,ROW_LENGTH,
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
     &                    LONGITUDE_STEP_INVERSE,  
     &                    SEC_P_LATITUDE,ADVECTION_TIMESTEP,
     &                    PSTAR,DELTA_AK,DELTA_BK,RS)

      IMPLICIT NONE

      INTEGER
     & P_FIELD             !IN DIMENSION OF FIELDS ON PRESSURE GRID.
     &,U_FIELD             !IN DIMENSION OF FIELDS ON VELOCITY GRID.
     &,ROW_LENGTH          !IN NUMBER OF POINTS PER ROW.
     &,P_LEVELS            !IN NUMBER OF PRESSURE LEVELS.
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

      INTEGER
     & TRACER_EW_SWEEPS(glsize(2),P_LEVELS) ! OUT. 
!  Number of East-West sweeps required for each row (of global field)

      REAL
     & U(U_FIELD,P_LEVELS)        !IN ADVECTING U FIELD, MASS-WEIGHTED.
     &,SEC_P_LATITUDE(P_FIELD)    !IN 1/COS(LAT) AT P POINTS
     &,RS(P_FIELD,P_LEVELS)       !IN EFFECTIVE RADIUS OF EARTH

      REAL
     & LONGITUDE_STEP_INVERSE     !IN 1/(DELTA LAMDA)
     &,ADVECTION_TIMESTEP         !IN
     &,PSTAR(P_FIELD)             !IN
     &,DELTA_AK(P_LEVELS)         !IN
     &,DELTA_BK(P_LEVELS)         !IN

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS:  4 ARE REQUIRED

      REAL
     &  U_WORK(U_FIELD)
     &, COURANT_NUMBER(P_FIELD)
     &, PSTAR_UV(U_FIELD)
     &, RS_UV(U_FIELD)
      INTEGER
     & info   ! Return code from GCom routines.
     &,LOCAL_EW_SWEEPS(P_LEVELS,P_FIELD/ROW_LENGTH) ! for this PE
     &,ALL_EW_SWEEPS(P_LEVELS,glsize(2)) ! to hold values from all PEs
! N.B. so that data are contiguous for inter-PE message passing routines
! GCG_IMAX & GC_IBCAST, these arrays are declared (levels,rows)
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     &  I,J,K           ! Do loop counters.
     &, I_START
     &, I_MAX
     &, P_ROWS
     &, global_row   ! row number in global array ALL_EW_SWEEPS
     &, HALF_P_ROWS  ! half of total no. of rows in global array

      REAL
     &  MAX_COURANT

C*L   EXTERNAL SUBROUTINE CALLS:------------------------------------

      EXTERNAL
     & P_TO_UV

C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS U_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 0.     INITIALISATION
CL---------------------------------------------------------------------

      P_ROWS = P_FIELD/ROW_LENGTH
      HALF_P_ROWS=glsize(2)/2    ! half total number of rows 
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE NUMBER OF SWEEPS REQUIRED ON EACH ROW.
CL---------------------------------------------------------------------

CL    Interpolate pressure field to velocity points.

      CALL P_TO_UV(PSTAR,PSTAR_UV,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)

CL    Loop over all levels.

      DO K=1,P_LEVELS
! Initialise polar/North- & South-most rows.
        ALL_EW_SWEEPS(K,1) = 1
        ALL_EW_SWEEPS(K,glsize(2)) = 1
        LOCAL_EW_SWEEPS(K,1) = 1
        LOCAL_EW_SWEEPS(K,P_ROWS) = 1
        IF (at_top_of_LPG) LOCAL_EW_SWEEPS(K,2) = 1
        IF (at_base_of_LPG) LOCAL_EW_SWEEPS(K,P_ROWS-1) = 1

CL    Interpolate RS field to velocity points.

      CALL P_TO_UV(RS(1,K),RS_UV,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)

CL    Remove mass-weight from U.

!!!        DO I=1,U_FIELD                                               
        DO  I = FIRST_VALID_PT,LAST_U_FLD_PT
          U_WORK(I) = U(I,K)/
     &                (RS_UV(I)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I)))
        END DO

CL    Calculate Courant number on each interior P_ROW.

!!!        DO I=ROW_LENGTH+1,P_FIELD-ROW_LENGTH                         
        DO  I = START_POINT_NO_HALO,END_P_POINT_NO_HALO
          COURANT_NUMBER(I) = .5*(U_WORK(I) + U_WORK(I-ROW_LENGTH))
     &                        *ADVECTION_TIMESTEP*SEC_P_LATITUDE(I)
     &                        *LONGITUDE_STEP_INVERSE/RS(I,K)
        END DO

CL    Loop over all rows.

!!!        DO J=2,P_ROWS-1                                              
        DO  J = FIRST_ROW,P_LAST_ROW
CL    Calculate maximum absolute courant number.

          MAX_COURANT = 0.0
          DO  I = (J-1)*ROW_LENGTH+Offx+1,J*ROW_LENGTH-Offx
            MAX_COURANT = MAX(ABS(COURANT_NUMBER(I)),MAX_COURANT)
          END DO
CL    Set number of sweeps so that maximum courant number on each row
Cl    is less than 0.25

          LOCAL_EW_SWEEPS(K,J) = 1 + 4*MAX_COURANT

CL    End loop over rows.
        END DO
CL    End loop over levels.                                             
        END DO                                                          

! Find max value of local_ew_sweeps for each row along all processors
!  in group gc_proc_row_group
      CALL GCG_IMAX((P_ROWS-2*Offy)*P_LEVELS,gc_proc_row_group,info,
     &               LOCAL_EW_SWEEPS(1,1+Offy))
! Copy to correct place in global array all_ew_sweeps
      DO  J = FIRST_ROW,P_LAST_ROW
        global_row = J+datastart(2)-Offy-1
        DO  K = 1,P_LEVELS
          ALL_EW_SWEEPS(K,global_row) = LOCAL_EW_SWEEPS(K,J)
        END DO
      END DO
! Broadcast section of global array from this processor to all others
! - only needs to be done by one processor per row.
      DO  I = 0,nproc-1,nproc_x
        CALL GC_IBCAST(I,P_LEVELS*g_blsizep(2,I),I,nproc,info,
     &                 ALL_EW_SWEEPS(1,g_datastart(2,I)))
      END DO

CL    Loop over all levels.   
      DO K=1,P_LEVELS                                                   
                                                                        
CL    Make number of sweeps in each hemisphere monotonic increasing
Cl    as you go towards the pole.
        DO  J = HALF_P_ROWS,2,-1
          IF (ALL_EW_SWEEPS(K,J) .lt. ALL_EW_SWEEPS(K,J+1)) THEN
            ALL_EW_SWEEPS(K,J) = ALL_EW_SWEEPS(K,J+1)
          END IF
        END DO
        DO  J = HALF_P_ROWS+1,glsize(2)-1
          IF (ALL_EW_SWEEPS(K,J) .lt. ALL_EW_SWEEPS(K,J-1)) THEN
            ALL_EW_SWEEPS(K,J) = ALL_EW_SWEEPS(K,J-1)
          END IF
        END DO
! Initialise North- & South-most values
        ALL_EW_SWEEPS(K,1) = 1
        ALL_EW_SWEEPS(K,glsize(2)) = 1
! Copy from local global array all_ew_sweeps to array tracer_ew_sweeps
! that will be passed to TRAC_ADV. 
      DO  J = 1,glsize(2)
        TRACER_EW_SWEEPS(J,K) = ALL_EW_SWEEPS(K,J)
      END DO

CL    End loop over levels.
      END DO

CL    END OF ROUTINE SET_TRAC

      RETURN
      END
