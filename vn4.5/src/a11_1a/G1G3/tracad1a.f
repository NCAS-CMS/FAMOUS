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
CLL   SUBROUTINE TRAC_ADV ---------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES ADVECTION INCREMENTS TO A FIELD AT A
CLL              SINGLE MODEL LEVEL USING A POSITIVE DEFINITE SCHEME.
CLL              IN CALCULATING THE INCREMENTS THE TEST FOR THE
CLL              DIRECTION OF THE WIND HAS BEEN REVERSED TO TAKE INTO
CLL              ACCOUNT THE CHANGE IN SIGN INTRODUCED BY MASS
CLL              WEIGHTING.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL
CLL M.MAWSON    <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.2  20/08/96  MPP code added.  RTHBarnes.
CLL  4.3  17/03/97  Corrections to MPP code.  RTHBarnes.
CLL       WARNING   Owing to compiler optimisation differences,
CLL       non-MPP & MPP runs will only bit compare if this deck
CLL       is compiled with minimum optimisation (-Oscalar0).
!LL  4.4  05/09/97  Ensure halos OK both at start and end of routine.
!LL                 S.D.Mullerworth
!LL  4.5  23/06/98  Single PE optimisations
!LL                 D.Salmond, B.Carruthers and P.Burton
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
      SUBROUTINE TRAC_ADV
     &                   (FIELD,N_SWEEP,U_MEAN,V_MEAN,U_FIELD,P_FIELD,
     &                    ADVECTION_TIMESTEP,ROW_LENGTH,
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
     &                    SEC_P_LATITUDE,COS_P_LATITUDE,RS,PSTAR,
     &                    DELTA_AK,DELTA_BK,LATITUDE_STEP_INVERSE,   
     &                    LONGITUDE_STEP_INVERSE,L_SUPERBEE)

      IMPLICIT NONE

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
      LOGICAL
     & L_SUPERBEE          !IN True then use SUPERBEE limiter,
     &                     !   False then use VAN LEER limiter.

      INTEGER
     & P_FIELD             !IN DIMENSION OF FIELDS ON PRESSURE GRID.
     &,U_FIELD             !IN DIMENSION OF FIELDS ON VELOCITY GRID.
     &,ROW_LENGTH          !IN NUMBER OF POINTS PER ROW.
     &,N_SWEEP(glsize(2))  ! Number of sweeps to be done East-West
!             ! for each row in full domain (needed for MAX_SWEEPS)

      REAL
     & U_MEAN(U_FIELD)        !IN ADVECTING U FIELD, MASS-WEIGHTED.
     &,V_MEAN(U_FIELD)        !IN ADVECTING V FIELD, MASS-WEIGHTED.
     &,FIELD(P_FIELD)         !IN FIELD TO BE ADVECTED.
     &,ADVECTION_TIMESTEP     !IN

      REAL
     & LONGITUDE_STEP_INVERSE     !IN 1/(DELTA LAMDA)
     &,LATITUDE_STEP_INVERSE      !IN 1/(DELTA PHI)
     &,SEC_P_LATITUDE(P_FIELD)    !IN 1/COS(LAT) AT P POINTS
     &,COS_P_LATITUDE(P_FIELD)    !IN COS(LAT) AT P POINTS
     &,RS(P_FIELD)                !IN RS_FIELD
     &,PSTAR(P_FIELD)             !IN
     &,DELTA_AK                   !IN
     &,DELTA_BK                   !IN

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 15 ARE REQUIRED

      REAL
     & FLUX_DELTA_T(P_FIELD) ! FLUX * ADVECTION TIMESTEP
     &,B1(P_FIELD)           ! ARGUMENT OF B_TERM
     &,B2(P_FIELD)           ! ARGUMENT OF B_TERM
     &,B_TERM(P_FIELD)       !
     &,COURANT(P_FIELD)      ! COURANT NUMBER
     &,ABS_COURANT(P_FIELD)  ! ABSOLUTE VALUE OF COURANT NUMBER
     &,COURANT_MW(P_FIELD)   ! MASS WEIGHTED COURANT NUMBER
     &,RS_SQUARED_DELTAP(P_FIELD) ! MASS * RADIUS OF EARTH
     &,MW(P_FIELD)           ! MASS WEIGHTING ASSOCIATED WITH
     &                       ! GRID BOX BOUNDARY FLUXES
     &,MW_RECIP(P_FIELD)     ! 1./MW
     &,RS_SQUARED_DELTAP_RECIP(P_FIELD) ! HOLDS 1./RS_SQUARED_DELTAP
     &,FIELD_INC(P_FIELD)    ! HOLDS INCREMENT TO FIELD.
     &,B_SWITCH(P_FIELD)     ! Entropy condition switch.
     &,SHIFT_N(ROW_LENGTH,2)   ! Local copy of polar rows for 180 deg
     &,B2_SHIFT_N(ROW_LENGTH)  ! rotational shift by GCG_RVECSHIFT
     &,SHIFT_S(ROW_LENGTH,2)   ! Local copy of polar rows for 180 deg
     &,B2_SHIFT_S(ROW_LENGTH)  ! rotational shift by GCG_RVECSHIFT

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     &  P_POINTS_UPDATE ! NUMBER OF P POINTS TO BE UPDATED
     & ,START_P_UPDATE  ! FIRST P POINT TO BE UPDATED
     & ,END_P_UPDATE    ! LAST P POINT TO BE UPDATED
     & ,ROWS            ! NUMBER OF ROWS TO BE UPDATED
     & ,I,J,L           ! Do loop counters.
     & ,I_SWEEP         ! holds east-west sweep number
     & ,MAX_SWEEPS      ! holds maximum number of east-west sweeps
     & ,I_CNTL          ! Control variable for do loop inside do while.
     & ,START_P(2)      ! Start point for loop with I_CNTL =1 or 2
     & ,END_P(2)        ! End point for loop with I_CNTL=1 or 2
     & ,U_ROWS          ! Number of u rows
     & ,N_HEMI          ! Number of hemispheres east-west advection
     &                  ! being performed in.
     & ,J1,J2,J3        ! for deciding which rows still need E-W sweeps
     & ,info            ! return code for GCom routines
     & ,HALF_RL         ! ROW_LENGTH/2 for GCG_RVECSHIFT function.

      REAL
     & NORTH_POLE_INC
     &,SOUTH_POLE_INC
     &,ROW_LENGTH_RECIP     ! 1/ROW_LENGTH
     &,R_SWEEP              ! 1/N_SWEEP(J-1)

      INTEGER II
      REAL real_ROW_LENGTH
C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 0.     INITIALISATION
CL---------------------------------------------------------------------

      ROWS             = P_FIELD/ROW_LENGTH
      P_POINTS_UPDATE  = (ROWS-2) * ROW_LENGTH
      IF (at_top_of_LPG) P_POINTS_UPDATE = P_POINTS_UPDATE-ROW_LENGTH
      IF (at_base_of_LPG) P_POINTS_UPDATE = P_POINTS_UPDATE-ROW_LENGTH
      START_P_UPDATE   = (FIRST_ROW-1) * ROW_LENGTH + 1
      END_P_UPDATE     = START_P_UPDATE + P_POINTS_UPDATE - 1
      I_SWEEP          = 1
      MAX_SWEEPS       = MAX(N_SWEEP(2),N_SWEEP(glsize(2)-1))
      ROW_LENGTH_RECIP = 1./GLOBAL_ROW_LENGTH 
! Ensure halos are OK before starting
      real_ROW_LENGTH=ROW_LENGTH

CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE FIELD INCREMENTS FOR U ADVECTION
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE COURANT NUMBER
C----------------------------------------------------------------------

C DIVIDE BY N_SWEEP TO ENSURE COURANT NUMBER LESS THAN 1.

      DO  J = datastart(2)+FIRST_ROW-1,datastart(2)+P_LAST_ROW-1
        J1 = (J-datastart(2))*ROW_LENGTH
! Because of unused 1st row of global fields.
          R_SWEEP = 1.0/N_SWEEP(J-1)
        DO I = 1,ROW_LENGTH   
          COURANT_MW(J1+I) = 0.5*(U_MEAN(J1+I-ROW_LENGTH)+U_MEAN(J1+I))
     &                * ADVECTION_TIMESTEP * LONGITUDE_STEP_INVERSE 
     &                * SEC_P_LATITUDE(J1+I) * R_SWEEP 
        END DO
      END DO                                                            

      DO I=START_P_UPDATE-ROW_LENGTH,END_P_UPDATE+ROW_LENGTH            
!      DO I=1,P_FIELD                                                   
        RS_SQUARED_DELTAP(I) = RS(I)*RS(I)*(DELTA_AK+DELTA_BK*PSTAR(I))
        RS_SQUARED_DELTAP_RECIP(I) = 1./RS_SQUARED_DELTAP(I)
      END DO

      DO I = START_P_UPDATE,END_P_UPDATE
        MW(I)=0.5*(RS_SQUARED_DELTAP(I)+RS_SQUARED_DELTAP(I+1))
        II=MIN(0.0,SIGN(1.0,COURANT_MW(I)))
        COURANT(I) = COURANT_MW(I)*RS_SQUARED_DELTAP_RECIP(I-II)
      END DO

! Do EW Swapbounds to fill in right halo points.
      CALL SWAPBOUNDS(MW,ROW_LENGTH,ROWS,EW_Halo,0,1) !single level
      CALL SWAPBOUNDS(COURANT,ROW_LENGTH,ROWS,EW_Halo,0,1) 
      MW(START_P_UPDATE-1)=1.
      MW(END_P_UPDATE+1)=1.
      MW_RECIP(START_P_UPDATE-1) = 1.
      MW_RECIP(END_P_UPDATE+1) = 1.
      COURANT(START_P_UPDATE-1)=0.
      COURANT(END_P_UPDATE+1)=0.

C SET ABSOLUTE VALUE OF COURANT NUMBER. THIS LOOP IS SEPARATE   
C TO LOOPS CALCULATING COURANT NUMBER SINCE INCLUDING IT THERE
C PREVENTS TOTAL OPTIMISATION.

      DO I= START_P_UPDATE,END_P_UPDATE
        ABS_COURANT(I) = ABS(COURANT(I))
        MW_RECIP(I) = 1./MW(I)
      END DO

      ABS_COURANT(START_P_UPDATE-1) = COURANT(START_P_UPDATE-1)
      ABS_COURANT(END_P_UPDATE+1) = COURANT(END_P_UPDATE+1)

CL    PERFORM N_SWEEPS OF ADVECTION ON EACH ROW.
CL    LOOP OVER NUMBER OF SWEEPS REQUIRED.

      DO I_SWEEP = 1,MAX_SWEEPS
        N_HEMI = 1

! Loop over non-halo rows to decide which points need computed this 
!  sweep. 
C IF FIRST SWEEP THEN PERFORM ADVECTION OVER ALL POINTS.
        IF(I_SWEEP.EQ.1) THEN
          START_P(1) = START_P_UPDATE
          END_P(1) = END_P_UPDATE
        ELSE
! Initialise
          J1 = 0
          J2 = P_LAST_ROW
          J3 = 0
          START_P(1) = 0

          DO  J = FIRST_ROW,P_LAST_ROW
            I = J+datastart(2)-Offy-1
!  above line because N_SWEEP is global array of values for all rows
            IF (J1 .eq. 0 .and. N_SWEEP(I) .ge. I_SWEEP)  J1 = J
            IF (J1 .ne. 0 .and. N_SWEEP(I) .lt. I_SWEEP .and.
     &          J2 .eq. P_LAST_ROW)  J2 = J-1
            IF (J3 .eq. 0 .and. N_SWEEP(I) .ge. I_SWEEP .and.
     &          J2 .ne. P_LAST_ROW)  J3 = J
          END DO
          IF (J1 .gt. 0) THEN
            START_P(1) = START_POINT_NO_HALO+(J1-FIRST_ROW)*ROW_LENGTH
            END_P(1) = END_P_POINT_NO_HALO-(P_LAST_ROW-J2)*ROW_LENGTH
          END IF
          IF (J3 .gt. 0) THEN
            START_P(2) = START_POINT_NO_HALO+(J3-FIRST_ROW)*ROW_LENGTH
            END_P(2) = END_P_POINT_NO_HALO
            N_HEMI = 2
          END IF
        END IF ! I_SWEEP

CL    LOOP OVER NUMBER OF HEMISPHERES.

! There may be no work for some processors to do,
        IF (START_P(1) .gt. 0) THEN
        DO I_CNTL = 1,N_HEMI
          START_P_UPDATE = START_P(I_CNTL)
          END_P_UPDATE = END_P(I_CNTL)

C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE FLUX_DELTA_T AND B1
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS

          DO I=START_P_UPDATE,END_P_UPDATE
            FLUX_DELTA_T(I) = (FIELD(I+1) - FIELD(I)) *
     &                         COURANT_MW(I)
            FIELD_INC(I) = 0.0
            B1(I)=FLUX_DELTA_T(I)*0.5*(1.0-ABS_COURANT(I))
          END DO

! but all processors must call SWAPBOUNDS together,
! so also end possible loop over hemispheres.
        END DO 
        END IF
! Do EW Swapbounds to fill in right halo points.
      CALL SWAPBOUNDS(FLUX_DELTA_T(1),ROW_LENGTH,ROWS,EW_Halo,0,1)
        IF (START_P(1) .gt. 0) THEN
        DO I_CNTL = 1,N_HEMI   
          START_P_UPDATE = START_P(I_CNTL) 
          END_P_UPDATE = END_P(I_CNTL)        

C----------------------------------------------------------------------
CL    SECTION 1.3    CALCULATE B1 AND B2
C----------------------------------------------------------------------


C GLOBAL MODEL.
C LOOP OVER ALL POINTS

          FLUX_DELTA_T(START_P_UPDATE-1)=0.
          FLUX_DELTA_T(END_P_UPDATE+1)=0.
          DO I=START_P_UPDATE,END_P_UPDATE
            II=SIGN(1.0,COURANT_MW(I))
            B2(I) = FLUX_DELTA_T(I+II)*0.5*(MW(I)*MW_RECIP(I+II)-
     &              ABS_COURANT(I+II))
            B_SWITCH(I) = SIGN(1.,COURANT(I)*COURANT(I+II))
          END DO


C----------------------------------------------------------------------
CL    SECTION 1.4    CALCULATE B_TERM
C----------------------------------------------------------------------

          IF (L_SUPERBEE) THEN

CL    SUPERBEE LIMITER.

            DO I=START_P_UPDATE,END_P_UPDATE
              IF(ABS(B2(I)).GT.1.0E-8) THEN
                B_SWITCH(I) = B_SWITCH(I)*B1(I)/B2(I)
              IF (B_SWITCH(I).GT.0.5.AND.
     &                 B_SWITCH(I).LT.2.0) THEN
                B_TERM(I) = B2(I) * MAX(B_SWITCH(I),1.0)
              ELSE IF (B_SWITCH(I).LE.0.0) THEN
                B_TERM(I) = 0.0
              ELSE
                B_TERM(I) = 2.0 * B2(I) * MIN(B_SWITCH(I),1.0)
              END IF
              ELSE
                B_SWITCH(I) = 0.
                  B_TERM(I) = 0.0
              END IF
            END DO


          ELSE

CL    VAN LEER LIMITER.

C LOOP OVER ALL POINTS
            DO I=START_P_UPDATE,END_P_UPDATE
              B_TERM(I) = 0.0
              IF (B1(I)*B2(I)*B_SWITCH(I).GT.0.0)
     &           B_TERM(I) = 2.0*B1(I)*B2(I)*B_SWITCH(I)/
     &                       (B1(I)+B2(I)*B_SWITCH(I))
            END DO

          END IF

! All processors must call SWAPBOUNDS together,
! so also end possible loop over hemispheres.
        END DO 
        END IF
! Do EW Swapbounds to fill in left halo points.
      CALL SWAPBOUNDS(B_TERM(1),ROW_LENGTH,ROWS,EW_Halo,0,1)
        IF (START_P(1) .gt. 0) THEN
        DO I_CNTL = 1,N_HEMI   
          START_P_UPDATE = START_P(I_CNTL) 
          END_P_UPDATE = END_P(I_CNTL)        
                                                                        
C----------------------------------------------------------------------
CL    SECTION 1.5    CALCULATE INCREMENTS TO FIELD
C----------------------------------------------------------------------

          DO I=START_P_UPDATE,END_P_UPDATE
            FIELD_INC(I) = FIELD_INC(I) - B_TERM(I)
            IF (COURANT_MW(I).GE.0.0)
     &        FIELD_INC(I)= FIELD_INC(I)+ 2.*B_TERM(I)- FLUX_DELTA_T(I)
          END DO

          DO I=START_P_UPDATE,END_P_UPDATE-1
            FIELD_INC(I+1) = FIELD_INC(I+1)-B_TERM(I)
            IF (COURANT_MW(I).LT.0.0)
     &         FIELD_INC(I+1) = FIELD_INC(I+1) + 2.*B_TERM(I)-
     &                                          FLUX_DELTA_T(I)
          END DO


C----------------------------------------------------------------------
CL    SECTION 1.6    UPDATE FIELD
C----------------------------------------------------------------------

C UPDATE MASS WEIGHTING FIELDS
          COURANT_MW(START_P_UPDATE-1) = 0.
          DO I=START_P_UPDATE,END_P_UPDATE
            RS_SQUARED_DELTAP(I) = RS_SQUARED_DELTAP(I) +
     &        (COURANT_MW(I-1)  - COURANT_MW(I))
            RS_SQUARED_DELTAP_RECIP(I) = 1./RS_SQUARED_DELTAP(I)
          END DO

C ADD INCREMENTS TO FIELD
          DO I=START_P_UPDATE,END_P_UPDATE
            FIELD(I) = FIELD(I)+FIELD_INC(I)*RS_SQUARED_DELTAP_RECIP(I)
          END DO

CL    END INNER LOOP OVER NUMBER OF HEMISPHERES.
        END DO
        END IF ! START_P(1) > 0
! Do EW Swapbounds to fill in left & right halo points.
      IF (I_SWEEP .lt. MAX_SWEEPS) THEN
        CALL SWAPBOUNDS(FIELD,ROW_LENGTH,ROWS,EW_Halo,0,1)
      END IF

CL    END LOOP OVER NUMBER OF SWEEPS
      END DO

! Swap all halo points to ensure updated arrays are fully up-to-date
!  for North-South sweep.
      CALL SWAPBOUNDS(FIELD,ROW_LENGTH,ROWS,EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(RS_SQUARED_DELTAP,ROW_LENGTH,ROWS,
     &                EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(RS_SQUARED_DELTAP_RECIP,ROW_LENGTH,ROWS,
     &                EW_Halo,NS_Halo,1)
!!! Might be better to recompute RS_SQUARED_DELTAP_RECIP haloes.
CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE FIELD INCREMENTS FOR V ADVECTION
CL---------------------------------------------------------------------

      DO I=1,P_FIELD
        FIELD_INC(I)=0.0
      END DO

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE COURANT NUMBER
C----------------------------------------------------------------------

      DO  I = FIRST_VALID_PT+1,LAST_U_VALID_PT  
!      DO I=2,U_FIELD 
        COURANT_MW(I) = 0.5*(V_MEAN(I)+V_MEAN(I-1))*ADVECTION_TIMESTEP*
     &                    LATITUDE_STEP_INVERSE
      END DO

! Do EW Swapbounds to fill in west halo points.
      CALL SWAPBOUNDS(COURANT_MW,ROW_LENGTH,ROWS,EW_Halo,0,1)

      DO  I = FIRST_VALID_PT,LAST_U_FLD_PT 
!      DO I=1,U_FIELD 
        MW(I)=0.5*(RS_SQUARED_DELTAP(I)*COS_P_LATITUDE(I)+              
     &     RS_SQUARED_DELTAP(I+ROW_LENGTH)*COS_P_LATITUDE(I+ROW_LENGTH))
! Split this loop to try and retain MPP & nonMPP bit comparison
        COURANT(I) = COURANT_MW(I)*RS_SQUARED_DELTAP_RECIP(I)
        COURANT(I) = COURANT(I)*SEC_P_LATITUDE(I)
        IF (COURANT_MW(I).GT.0.) THEN
          COURANT(I) = COURANT_MW(I)*SEC_P_LATITUDE(I+ROW_LENGTH)
     &                 *RS_SQUARED_DELTAP_RECIP(I+ROW_LENGTH)
        ENDIF
      END DO

C ABSOLUTE VALUE OF COURANT NUMBER CALCULATED IN THIS LOOP AS PUTTING
C IT IN PREVIOUS LOOP PREVENTS FULL OPTIMISATION.

! Do NS Swapbounds to fill in south halo points.
      CALL SWAPBOUNDS(MW,ROW_LENGTH,ROWS,0,NS_Halo,1) !single level
      CALL SWAPBOUNDS(COURANT,ROW_LENGTH,ROWS,0,NS_Halo,1) 

      DO  I = FIRST_VALID_PT,LAST_U_VALID_PT 
        ABS_COURANT(I) = ABS(COURANT(I))
        MW_RECIP(I) = 1./MW(I)
      END DO

C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE FLUX_DELTA_T AND B1
C----------------------------------------------------------------------

      DO  I = FIRST_VALID_PT,LAST_U_FLD_PT 
!      DO I=1,U_FIELD   
        FLUX_DELTA_T(I) = COURANT_MW(I) * (FIELD(I)-FIELD(I+ROW_LENGTH))
      END DO

! Do NS Swapbounds to fill in south halo points.
      CALL SWAPBOUNDS(FLUX_DELTA_T(1),ROW_LENGTH,ROWS,0,NS_Halo,1)
      DO  I = FIRST_VALID_PT,LAST_U_FLD_PT 
        B1(I) = FLUX_DELTA_T(I)*0.5*(1.0-ABS_COURANT(I))                
      END DO                                                            
                                                                        
C----------------------------------------------------------------------
CL    SECTION 2.3    CALCULATE B1 AND B2
C----------------------------------------------------------------------

C CALCULATE FLUXES AT VELOCITY POINTS.
C FIRST LOOP OVER ALL POINTS NOT AT NORTHERN OR SOUTHERN BOUNDARY.

      real_row_length=row_length
cdir$ unroll
      DO  I = START_POINT_NO_HALO,END_U_POINT_NO_HALO
        ii=nint(sign(real_row_length, courant_mw(i)))
        B2(I) = FLUX_DELTA_T(I-II)*0.5*(MW(I)
     &           *MW_RECIP(I-II) - ABS_COURANT(I-II))
      END DO
c
      DO  I = START_POINT_NO_HALO,END_U_POINT_NO_HALO
        ii=nint(sign(real_row_length, courant_mw(i)))
        B_SWITCH(I) = SIGN(1.,COURANT(I)*COURANT(I-II))
      END DO

      IF (at_top_of_LPG) THEN
! Needs values over top of pole on another processor
      HALF_RL = GLOBAL_ROW_LENGTH/2
! Copy North polar row into copy arrays
      DO  I = 1,ROW_LENGTH
        SHIFT_N(I,1) = MW(TOP_ROW_START+I-1)
        SHIFT_N(I,2) = COURANT(TOP_ROW_START+I-1)
      END DO
! Rotate these arrays by half the global row length to get values on
!  opposite side of the pole.
      CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*Offx,1+Offx,2,
     &                   HALF_RL,.TRUE.,SHIFT_N,GC_ROW_GROUP,info)  
      DO  I = 1,ROW_LENGTH
        B2_SHIFT_N(I) = 
     &    FLUX_DELTA_T(TOP_ROW_START+I-1)*0.5*(SHIFT_N(I,1)
     &   *MW_RECIP(TOP_ROW_START+I-1)-ABS_COURANT(TOP_ROW_START+I-1))
        B_SWITCH(TOP_ROW_START+I-1) =
     & SIGN(1.0,COURANT(TOP_ROW_START+I-1)*SHIFT_N(I,2))
      END DO
      CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*Offx,1+Offx,1,
     &                   HALF_RL,.TRUE.,B2_SHIFT_N,GC_ROW_GROUP,info)  
      DO  I = 1,ROW_LENGTH
        B2(TOP_ROW_START+I-1) = B2_SHIFT_N(I)
      END DO
      DO  I = TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
        IF (COURANT_MW(I) .lt. 0.0) THEN
          B2(I) = FLUX_DELTA_T(I+ROW_LENGTH)*0.5*(MW(I) 
     &          *MW_RECIP(I+ROW_LENGTH) - ABS_COURANT(I+ROW_LENGTH))  
          B_SWITCH(I) = SIGN(1.,COURANT(I)*COURANT(I+ROW_LENGTH)) 
        END IF 
      END DO
      END IF ! at_top_of_LPG

      IF (at_base_of_LPG) THEN
      DO  I = U_BOT_ROW_START,LAST_U_VALID_PT
        B2(I) = FLUX_DELTA_T(I-ROW_LENGTH)*0.5*(MW(I)
     &          *MW_RECIP(I-ROW_LENGTH)-ABS_COURANT(I-ROW_LENGTH))
        B_SWITCH(I) = SIGN(1.0,COURANT(I)*COURANT(I-ROW_LENGTH))
      END DO
      HALF_RL = GLOBAL_ROW_LENGTH/2
      DO  I = 1,ROW_LENGTH
        SHIFT_S(I,1) = MW(U_BOT_ROW_START+I-1)
        SHIFT_S(I,2) = COURANT(U_BOT_ROW_START+I-1)
      END DO
      CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*Offx,1+Offx,2,
     &                   HALF_RL,.TRUE.,SHIFT_S,GC_ROW_GROUP,info)  
      DO  I = 1,ROW_LENGTH
        B2_SHIFT_S(I) = 
     &    FLUX_DELTA_T(U_BOT_ROW_START+I-1)*0.5*(SHIFT_S(I,1)
     & *MW_RECIP(U_BOT_ROW_START+I-1)-ABS_COURANT(U_BOT_ROW_START+I-1))
        IF (COURANT_MW(U_BOT_ROW_START+I-1) .lt. 0.0) THEN
        B_SWITCH(U_BOT_ROW_START+I-1) = SIGN(1.0,
     &   COURANT(U_BOT_ROW_START+I-1)*SHIFT_S(I,2))
        END IF
      END DO
      CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*Offx,1+Offx,1,
     &                   HALF_RL,.TRUE.,B2_SHIFT_S,GC_ROW_GROUP,info)  
      DO  I = 1,ROW_LENGTH
        IF (COURANT_MW(U_BOT_ROW_START+I-1) .lt. 0.0) THEN
        B2(U_BOT_ROW_START+I-1) = B2_SHIFT_S(I)
        END IF
      END DO
      END IF ! at_base_of LPG


C----------------------------------------------------------------------
CL    SECTION 2.4    CALCULATE B_TERM
C----------------------------------------------------------------------

      IF (L_SUPERBEE) THEN

CL    SUPERBEE LIMITER.

        DO  I = FIRST_FLD_PT,LAST_U_FLD_PT
!        DO I=1,U_FIELD    
          IF(ABS(B2(I)).GT.1.0E-8) THEN
            B_SWITCH(I) = B_SWITCH(I)*B1(I)/B2(I)
            IF (B_SWITCH(I).GT.0.5.AND.B_SWITCH(I).LT.2.0) THEN
              B_TERM(I) = B2(I) * MAX(B_SWITCH(I),1.0)
            ELSE IF (B_SWITCH(I).LE.0.0) THEN
              B_TERM(I) = 0.0
            ELSE
              B_TERM(I) = 2.0 * B2(I) * MIN(B_SWITCH(I),1.0)
            END IF
          ELSE
            B_SWITCH(I) = 0.
            B_TERM(I) = 0.0
          END IF
        END DO

      ELSE

CL    VAN LEER LIMITER.

C LOOP OVER ALL POINTS
        DO  I = FIRST_FLD_PT,LAST_U_FLD_PT
!        DO I=1,U_FIELD  
          B_TERM(I) = 0.0
          IF (B1(I)*B2(I)*B_SWITCH(I).GT.0.0)
     &           B_TERM(I) = 2.0*B1(I)*B2(I)*B_SWITCH(I)/
     &                       (B1(I)+B2(I)*B_SWITCH(I))
        END DO

      END IF

! Do NS Swapbounds to fill in halo points.
      CALL SWAPBOUNDS(B_TERM(1),ROW_LENGTH,ROWS,0,NS_Halo,1)
                                                                        
C----------------------------------------------------------------------
CL    SECTION 2.5    CALCULATE INCREMENTS TO FIELD
CL---------------------------------------------------------------------

CDIR$ IVDEP
      DO  I = FIRST_FLD_PT,LAST_U_FLD_PT
!      DO I=1,U_FIELD 
        FIELD_INC(I) = FIELD_INC(I) - B_TERM(I)
        IF (COURANT_MW(I).LT.0.0)
     &    FIELD_INC(I) = FIELD_INC(I) - FLUX_DELTA_T(I) +2.*B_TERM(I)
      END DO
CDIR$ IVDEP
      DO  I = FIRST_VALID_PT,LAST_U_FLD_PT
!      DO I=1,U_FIELD 
        FIELD_INC(I+ROW_LENGTH) = FIELD_INC(I+ROW_LENGTH) -
     &                              B_TERM(I)
        IF (COURANT_MW(I).GE.0.0)
     &    FIELD_INC(I+ROW_LENGTH) = FIELD_INC(I+ROW_LENGTH) -
     &                              FLUX_DELTA_T(I) + 2.*B_TERM(I)
      END DO


C----------------------------------------------------------------------
CL    SECTION 2.6    CALCULATE POLAR INCREMENTS
CL---------------------------------------------------------------------

C CALCULATE AVERAGE POLAR INCREMENT
      NORTH_POLE_INC = 0.0
      SOUTH_POLE_INC = 0.0

! Use reproducible vector sum of points on polar rows.
      IF (at_top_of_LPG) THEN
        CALL GCG_RVECSUMR(P_FIELD,ROW_LENGTH-2*EW_Halo,
     &                    TOP_ROW_START+EW_Halo,1,
     &                    FIELD_INC,GC_ROW_GROUP,info,NORTH_POLE_INC)
      END IF
! Because values are summed in reverse order in non-MPP loop
! SOUTH_POLE_INC will probably not bit-compare
      IF (at_base_of_LPG) THEN
        CALL GCG_RVECSUMR(P_FIELD,ROW_LENGTH-2*EW_Halo,
     &                    P_BOT_ROW_START+EW_Halo,1,
     &                    FIELD_INC,GC_ROW_GROUP,info,SOUTH_POLE_INC)
      END IF

      NORTH_POLE_INC = NORTH_POLE_INC * ROW_LENGTH_RECIP * 2.0
      SOUTH_POLE_INC = SOUTH_POLE_INC * ROW_LENGTH_RECIP * 2.0

      IF (at_top_of_LPG) THEN
        DO  I = TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
          FIELD_INC(I) =  NORTH_POLE_INC 
        END DO
      END IF
      IF (at_base_of_LPG) THEN
        DO  J = P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
          FIELD_INC(J) = SOUTH_POLE_INC   
        END DO
      END IF

C----------------------------------------------------------------------
CL    SECTION 2.7    UPDATE FIELD
CL---------------------------------------------------------------------

C UPDATE MASS WEIGHTING
      DO  I = START_POINT_NO_HALO,END_P_POINT_NO_HALO
!      DO I=ROW_LENGTH+1,P_FIELD-ROW_LENGTH  
        RS_SQUARED_DELTAP(I) = RS_SQUARED_DELTAP(I) +
     &        (COURANT_MW(I)-COURANT_MW(I-ROW_LENGTH))*SEC_P_LATITUDE(I)
      END DO

C     POLAR VALUES
      NORTH_POLE_INC = 0.
      SOUTH_POLE_INC = 0.
! Use reproducible vector sum of points on polar rows.
      IF (at_top_of_LPG) THEN
        CALL GCG_RVECSUMR(P_FIELD,ROW_LENGTH-2*EW_Halo,
     &                    TOP_ROW_START+EW_Halo,1,
     &                    COURANT_MW,GC_ROW_GROUP,info,NORTH_POLE_INC)
      END IF
      IF (at_base_of_LPG) THEN
        CALL GCG_RVECSUMR(P_FIELD,ROW_LENGTH-2*EW_Halo,
     &                    U_BOT_ROW_START+EW_Halo,1,
     &                    COURANT_MW,GC_ROW_GROUP,info,SOUTH_POLE_INC)
      SOUTH_POLE_INC = -SOUTH_POLE_INC
      END IF

      NORTH_POLE_INC = NORTH_POLE_INC*ROW_LENGTH_RECIP*
     &                   SEC_P_LATITUDE(TOP_ROW_START)*2
      SOUTH_POLE_INC = SOUTH_POLE_INC*ROW_LENGTH_RECIP*
     &                   SEC_P_LATITUDE(P_BOT_ROW_START)*2
      IF (at_top_of_LPG) THEN
        DO  I = TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
          RS_SQUARED_DELTAP(I) = RS_SQUARED_DELTAP(I) + NORTH_POLE_INC  
        END DO
      END IF
      IF (at_base_of_LPG) THEN
        DO  J = P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
          RS_SQUARED_DELTAP(J) = RS_SQUARED_DELTAP(J) + SOUTH_POLE_INC 
        END DO
      END IF

C ADD INCREMENTS TO FIELD
      DO  I = FIRST_FLD_PT,LAST_P_FLD_PT
!      DO I=1,P_FIELD   
        FIELD(I) = FIELD(I)+
     &     FIELD_INC(I)*SEC_P_LATITUDE(I) / RS_SQUARED_DELTAP(I)
      END DO

CL    END OF ROUTINE TRAC_ADV

      RETURN
      END
