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
CLL   SUBROUTINE SETFILT ------------------------------------------
CLL
CLL   PURPOSE:
CLL   SUBROUTINE 'SETFILT' - COMPUTES FACTORS OF N & TRIGONOMETRIC
CLL   FUNCTIONS REQUIRED BY FOURIER, FFT99 & FFT991
CLL   UNIFIED MODEL RE-WRITE OF ECMWF ROUTINE SET99
CLL   ALSO CALCULATES TWO_D_GRID_CORRECTION FACTORS USED IN SET_FIL
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   REWRITTEN TO UNIFIED MODEL PROGRAMMING STANDARDS FROM ECMWF
CLL   CODE BY M.H.MAWSON; ORIGINAL CODE WRITER C. TEMPERTON
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
!     4.1     19/06/95  Added MPP code and argument ROWS  P.Burton
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:   UNIFIED MODEL DOCUMENTATION PAPER NUMBER 10
CLL                    M.J.P. CULLEN, T. DAVIES AND M.H. MAWSON
CLL                    VERSION 8 DATED 1/05/90.
CLLEND-------------------------------------------------------------

C*L   ARGUMENT LIST

      SUBROUTINE SETFILT(TRIGS,IFAX,N,ROWS,TWO_D_GRID_CORRECTION,
     &                   P_FIELD, COS_P_LATITUDE)

      IMPLICIT NONE

      INTEGER
     *        IFAX(10) !OUT HOLDS FACTORS OF N

      INTEGER
     *        N       !IN NUMBER OF POINTS OF DATA ON A SLICE
     *,       ROWS   !IN NUMBER OF ROWS IN P_FIELD
     *,       P_FIELD !IN HOLDS NUMBER OF POINTS IN THE FIELD.

      REAL
     *     COS_P_LATITUDE(P_FIELD) !IN COSINE OF LATITUDE AT P POINTS.

      REAL
     *     TRIGS(N)   !OUT HOLDS TRIGONOMETRIC FUNCTIONS NEEDED BY
     *                !    FOURIER, FFT99 AND FFT991
     *,    TWO_D_GRID_CORRECTION(ROWS) !OUT HOLDS FACTOR ON
     *                ! EACH ROW WHICH MODIFIES THE STABILITY
     *                ! RELATION EQUATION (50) SO THAT IT TAKES INTO
     *                ! ACCOUNT THE 2-D NATURE OF THE GRID.
C*   --------------------------------------------------------------
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


C*L   WORKSPACE. 2 ARRAYS ARE REQUIRED ----------------------------
      INTEGER
     *        JFAX(10) ! HOLDS FACTORS OF N BEFORE THEY ARE ORDERED
     *                 ! AND STORED IN IFAX
     *        ,LFAX(8) ! HOLDS LIST OF VALID FACTORS. LAST VALUE IS
     *                 ! MINUS PENULTIMATE ONE. THIS IS USED TO
     *                 ! SEND CODE TO THE ERROR MESSAGE AS N CONTAINS
     *                 ! AN INVALID FACTOR.
C*   --------------------------------------------------------------

C*L   NO EXTERNAL ROUTINES  ---------------------------------------
C*  ---------------------------------------------------------------

C LOCAL VARIABLES
      REAL
     *  DEL    ! HOLDS ANGLE IN RADIANS BETWEEN POINTS ON LATITUDE
     *         ! CIRCLE
     * ,ANGLE  ! HOLDS ANGLE IN RADIANS AT A POINT ON THE LATITUDE
     *         ! CIRCLE

      INTEGER
     *  K      ! LOOP COUNTER
     * ,I      ! LOOP COUNTER
     * ,NU     ! HOLDS CURRENT FACTORISED VALUE OF N
     * ,IFAC   ! HOLDS CURRENT FACTOR BEING TESTED FOR
     * ,L      ! HOLDS CURRENT POSITION IN LFAX ARRAY
     * ,NFAX   ! HOLDS TOTAL NUMBER OF FACTORS OF N

CL
CL ----------------------------------------------------------------
CL    SECTION 1.  INITIALISATION.
CL ----------------------------------------------------------------

      LFAX(1)=6
      LFAX(2)=8
      LFAX(3)=5
      LFAX(4)=4
      LFAX(5)=3
      LFAX(6)=2
      LFAX(7)=1
      LFAX(8)=-1

CL
CL ----------------------------------------------------------------
CL    SECTION 2.  SET TRIGONOMETRIC FUNCTIONS
CL ----------------------------------------------------------------

      DEL=4.0*ASIN(1.0)/N
      DO 200 K=0,N/2-1
        ANGLE=K*DEL
! Set global_trigs array - carried by COMMON info FILTER
        global_trigs(2*K+1)=COS(ANGLE)
        global_trigs(2*K+2)=SIN(ANGLE)
 200  CONTINUE

CL
CL ----------------------------------------------------------------
CL    SECTION 3. FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
CL               STORE FACTORS IN DESCENDING ORDER.
CL ----------------------------------------------------------------

C LOOK FOR SIXES FIRST AND STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
C K HOLDS NUMBER OF FACTORS FOUND. L IS USED TO MOVE THROUGH LIST OF
C VALID FACTORS.
      K=0
      L=1
  300 CONTINUE
      IF (MOD(NU,IFAC).EQ.0) THEN
C IF IFAC IS A FACTOR OF N.
        K=K+1
        JFAX(K)=IFAC
        IF (IFAC.EQ.8.AND.K.NE.1) THEN
          JFAX(1)=8
          JFAX(K)=6
        ENDIF
        NU=NU/IFAC
C IF N FACTORISES COMPLETELY JUMP PAST ERROR MESSAGE
        IF(NU.EQ.1) GO TO 310
C IF FACTOR IS NOT AN 8 SEE IF IT APPEARS MORE THAN ONCE.
        IF(IFAC.NE.8) GO TO 300
      ENDIF
      L=L+1
      IFAC=LFAX(L)
C IF NOT COMPLETELY FACTORISED GO BACK TO TOP OF PROCEDURE
      IF(IFAC.GT.1) GO TO 300

C ILLEGAL FACTOR IN N. PRINT ERROR MESSAGE THEN JUMP TO END OF
C ROUTINE.

      WRITE(6,1300) N
 1300 FORMAT(' ERROR IN SETFILT. N = ',I4,' CONTAINS AN ILLEGAL FACTOR')
      GO TO 330

C NOW REVERSE ORDER OF FACTORS SO THAT THEY ARE IN ASCENDING ORDER

  310 CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 320 I=1,NFAX
        IFAX(NFAX+2-I)=JFAX(I)
 320  CONTINUE
      IFAX(10)=N
 330  CONTINUE

CL
CL ----------------------------------------------------------------
CL    SECTION 4. CALCULATE 2-D GRID CORRECTION TERM.
CL               SEE DOC. PAPER 10 SECTION 3.5.
CL ----------------------------------------------------------------

CL CALCULATION NOT PERFORMED ON POLAR ROWS.
      DO I=1+Offy,ROWS-Offy
        TWO_D_GRID_CORRECTION(I) =
     &   (1.+4.*(glsize(2)-1.)*(glsize(2)-1.)/(glsize(1)*glsize(1))*
     &   COS_P_LATITUDE(I*lasize(1))*COS_P_LATITUDE(I*lasize(1)))**0.5
      ENDDO

      IF (attop) TWO_D_GRID_CORRECTION(1+Offy)=0.0
      IF (atbase) TWO_D_GRID_CORRECTION(ROWS-Offy)=0.0

CL    END OF ROUTINE SETFILT
      RETURN
      END
