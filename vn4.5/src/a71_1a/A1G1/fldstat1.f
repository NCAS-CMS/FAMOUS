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
CLL  SUBROUTINES FLDSTAT AND FLDDIAG ---------------------------------
CLL
CLL  PURPOSE:
CLL   CALCULATE VALUES OF increments of T,RH,U,V between timesteps
CLL                FLDDIAG:
CLL   PRINT VALUES OF max,min increments of T,RH,U,V between timesteps
CLL
CLL  MODIFIED VERSION OF FLDDIAG FOR CRAY Y-MP BASED ON
CLL  EARLIER ROUTINE BY S.BELL WRITTEN BY F. RAWLINS
CLL
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL RR / DR     <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL 4.2  8 Jan 97  Changes for MPP. Gather domains from each pe to
CLL                provide full global fields on pe0 to be written
CLL                and read from disk files (1 per pe). R.Rawlins
CLL 4.3 15 May 97  Correction to 4.2 change: V increments against V
CLL                instead of U. Correct RH label. R.Rawlins
CLL 4.4 28 Aug 97  Change method of I/O from Fortran unformatted to
Cll                C buffer streams with portable I/O, thus freeing
!LL 4.5 13/01/98   Replace reference to IOVARS comdeck to ATM_LSM
!LL                                                      P.Burton
CLL 4.5 25 Mar 98  Change formatting of printed diagnostics to cater
CLL                for 10**7 points in horizontal field (from 10**5):
CLL                needed for new op. resolution. Rick Rawlins
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 5,
CLL  VERSION 4, DATED 31/05/90
CLL
CLL  SYSTEM TASK: increment diagnostics  D67
CLL
CLL  DOCUMENTATION:        None
CLL
CLLEND-------------------------------------------------------------

      SUBROUTINE FLDSTAT (NUP,NRP,ILENP,JLENP,LENTHP,LENUVP,KSTEP,
     +                    AK,BK,AKH,BKH,P_EXNER,
     +                    PSTAR,TH,Q,U,V,
     &                    LTHETA,PRFLD_STEP,PRFLD_FIRST,PRFLD_LAST,
     &                    NDEV_FLD,LEN_FLD_FILENAME,FLD_FILENAME)
C
C FLDSTAT   GETS STATS FOR MEANS/MAX/MIN OF PROGNOSTIC VARIABLES
C           PLUS STATS ON CHANGE SINCE LAST TIMESTEP
C           IT CALLS FLDDIAG AND ALSO DOES I/O TO UNIT NDEV_FLD 
C*
      IMPLICIT NONE

      EXTERNAL FLDDIAG,QSAT
C
C*L  ARGUMENTS:---------------------------------------------------

      INTEGER
     +      NUP,                  ! (IN) TOTAL NUMBER OF LEVELS
     +      NRP,                  ! (IN) NUMBER OF WET LEVELS
     +      ILENP,                ! (IN) NUMBER OF POINTS ON ROW
     +      JLENP,                ! (IN) NUMBER OF ROWS
     +      LENTHP,               ! (IN) NUMBER OF POINTS IN MASS FIELD
     +      LENUVP,               ! (IN) NUMBER OF POINTS IN WIND FIELD
     +      KSTEP,                ! (IN) CURRENT MODEL TIMESTEP
     +      PRFLD_STEP,           ! (IN) STEP INTERVAL FOR PRINTING
     +      PRFLD_FIRST,          ! (IN) FIRST STEP    FOR PRINTING
     &      PRFLD_LAST,           ! (IN) LAST STEP     FOR PRINTING
     &      NDEV_FLD,             ! (IN) OUTPUT DEVICE NUMBER
     &      LEN_FLD_FILENAME      ! (IN) Filename length of NDEV_FLD

      CHARACTER*80 FLD_FILENAME   ! (IN) Filename of NDEV_FLD file 
 
      LOGICAL
     +     LTHETA                      ! (IN) THETA OR TEMPERATURE
      REAL
     +     AK(NUP),BK(NUP),            ! (IN) HYBRID CO-ORDS - full levs
     +     AKH(NUP+1),BKH(NUP+1),      ! (IN) HYBRID CO-ORDS - 1/2 levs
     +     P_EXNER(LENTHP,NUP+1),      ! (IN) EXNER PRESSURE
     +     PSTAR(LENTHP),              ! (IN) PROG VARIABLE PSTAR
     +     TH   (LENTHP,NUP),          ! (IN) THETA  (LTHETA=.T. OR
C                                        TEMPERATURE (LTHETA=.F.)
     +     Q    (LENTHP,NRP),          ! (IN) PROG VARIABLE Q
     +     U    (LENUVP,NUP),          ! (IN) PROG VARIABLE U
     +     V    (LENUVP,NUP)           ! (IN) PROG VARIABLE V

! DECOMPTP comdeck
!
! Description
!
! Magic numbers indicating decomposition types.
! These numbers are used to index the arrays defined in the
! DECOMPDB comdeck, and are required as an argument to
! the CHANGE_DECOMPOSITION subroutine.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 4.3       17/02/97  Added new ocean decomposition decomp_nowrap_ocean
!                     which does not contain extra wrap points at
!                     start and end of row.                  P.Burton

! Magic Numbers indicating decomposition types

      INTEGER
     &  max_decomps            ! maximum number of decompositions
     &, decomp_unset           ! no decomposition selected
     &, decomp_standard_atmos  ! standard 2D atmosphere
!                              ! decomposition
     &, decomp_standard_ocean  ! standard 1D ocean decomposition
     &, decomp_nowrap_ocean    ! 1D ocean without extra wrap-around
!                              ! points at ends of each row

      PARAMETER (
     &  max_decomps=3
     &, decomp_unset=-1
     &, decomp_standard_atmos=1
     &, decomp_standard_ocean=2
     &, decomp_nowrap_ocean=3)

! End of DECOMPTP comdeck
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
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton

      INTEGER

     &  ROW_LENGTH_MAX  ! Maximum row length
     &, P_ROWS_MAX      ! Maximum number of rows
     &, HORIZ_DIM_MAX   ! MAX(ROW_LENGTH_MAX,P_ROWS_MAX)
     &, HALO_MAX        ! Maximum MPP halo width
     &, P_LEVELS_MAX    ! Maximum number of total levels
     &, Q_LEVELS_MAX    ! Maximum number of wet levels

      PARAMETER ( ROW_LENGTH_MAX = 432
     &,           P_ROWS_MAX = 325
     &,           HORIZ_DIM_MAX = 432
     &,           HALO_MAX = 2  ! fourth order double width halo
     &,           P_LEVELS_MAX = 42
     &,           Q_LEVELS_MAX = 42)

! Derived sizes

      INTEGER
     &  Max2DFieldSize
     &, Max3DFieldSize
     &, MaxHaloSize

      PARAMETER (
     &  Max2DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX
     &, Max3DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX*P_LEVELS_MAX
     &, MaxHaloSize = HORIZ_DIM_MAX*HALO_MAX
     & )

      INTEGER
     &  MAX_SHMEM_COMMON_SIZE,
     &  MaxFieldSize,
     &  MaxFieldSizeMes                                                 
      PARAMETER ( MAX_SHMEM_COMMON_SIZE = 3000000 ,
     &            MaxFieldSize   = MAX_SHMEM_COMMON_SIZE/4 ,
     &            MaxFieldSizeMes= MAX_SHMEM_COMMON_SIZE/6 )
!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!
! History:
!   Model    Date     Modification history
!   version
!   4.5      12/01/98 New comdeck created.                P.Burton
!

      LOGICAL
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

      COMMON /Atmos_LSM_Common/
     &  atmos_landmask
     &, atmos_landmask_local
     &, atmos_number_of_landpts

CDIR$ CACHE_ALIGN /Atmos_LSM_Common/

! End of comdeck ATM_LSM
C*L --------------------- Comdeck: CENVIR   ----------------------------
C
C    Purpose: COMDECK defining Character enviroment variables used
C             by portable IO to open and close files
C
C    Author : R A Stratton      Date : 22/10/92
C
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
CLL
CLL 3.1     15/01/93  Increase no. of unit nos. from 1-99  to 1-199
CLL                   Dummy names have been set up temporarily for
CLL                   files 104-119. R.Rawlins
CLL
CLL 3.3     09/03/94  Separate data statements into COMDECK
CLL                   CENVIRDT. Also includes mods originally
CLL                   in RB221193 : Add source terms at unit no.110
CLL                   P.Burton and R.T.H Barnes
CLL

C    Vn3.0  12/02/93 - Environment variables PERTURB and TRANSP put in
C                      positions 37 and 97 respectively in character
C                      array FT_ENVIRON, and the appropriate character
C                      lengths put in LEN_FT_ENVIR. C. S. Douglas
C
C  Type declarations
C
      CHARACTER*8 FT_ENVIRON(199)  ! Array holding enviroment variables
C                                   for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
C


C
C Common Blocks for character and integer arrays
C
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
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

 
C*L  WORKSPACE USAGE:-------------------------------------------------
C DYNAMIC SPACE FOR LAST TIMESTEP PROGNOSTIC VARIABLES
      REAL WORKPTR(LENTHP),WORKUV(LENUVP),
     +     P    (LENTHP),             ! WORKSPACE FOR PRESSURE
     +     T    (LENTHP,NUP),         ! WORKSPACE FOR TEMPERATURE
     +     RH   (LENTHP,NRP)          ! WORKSPACE FOR RELATIVE HUMIDITY
      REAL WORK_FULL(glsize(1)*glsize(2))
     &    ,WORKPTR_FULL(glsize(1)*glsize(2))
      INTEGER
     &     gather_pe
     &    ,info       ! return code for MPP gather


      INTEGER LEV,                    ! LEVEL COUNTER
     &        I,                      ! POINT COUNTER  
     &        ICODE,                  ! ERROR RETURN CODE FROM I/O
     &        LEN_IO,                 ! I/O LENGTH RETURNED FROM I/O
     &        IPOS                    ! I/O POINTER  

      REAL PLEV,PLEVP1                ! Pressures at half levels
                                      ! LEV and LEV+1
      REAL P_EXNER_FULL               ! Exner pressure at full model
                                      ! levels.
      REAL A_IO                       ! Error return from buffer i/o 

      LOGICAL FIRST                   ! FIRST TIME THROUGH CODE?
      DATA FIRST /.TRUE./
      SAVE FIRST 

      gather_pe=0     ! only PE 0 for MPP gathering

C  (  THE LOGICAL DEVICE IS OPENED AT TOP LEVEL IN ROUTINE INITIAL)
 
CL
CL GET TEMPERATURE OF EACH LEVEL
CL
      DO LEV=1,NUP
      IF(LTHETA) THEN
         DO I=1,LENTHP
         PLEVP1 = AKH(LEV+1) + BKH(LEV+1)*PSTAR(I)
         PLEV   = AKH(LEV)   + BKH(LEV)  *PSTAR(I)
         P_EXNER_FULL = P_EXNER_C
     +   (P_EXNER(I,LEV+1),P_EXNER(I,LEV),PLEVP1,PLEV,KAPPA)
         T(I,LEV) = TH(I,LEV) * P_EXNER_FULL
         ENDDO
      ELSE
         DO I=1,LENTHP
         T(I,LEV)=TH(I,LEV)
         ENDDO
      ENDIF
      ENDDO
CL
CL GET RELATIVE HUMIDITY FOR EACH LEVEL
CL
      DO LEV=1,NRP
        DO I=1,LENTHP
        P(I)=AK(LEV) + BK(LEV)*PSTAR(I)
        ENDDO
        CALL QSAT(WORKPTR,T(1,LEV),P,LENTHP)
        DO I=1,LENTHP
        RH(I,LEV)=Q(I,LEV)/WORKPTR(I)*100.0
        ENDDO
      ENDDO

      IF(FIRST)THEN
         FIRST=.FALSE.
         WRITE(6,*) ' First call to FLDSTAT at step ',KSTEP 
      ELSE
CL
CL READ PREVIOUS TIMESTEP AND CALL FLDDIAG TO GET STATS
CL
       IF(MOD(KSTEP-PRFLD_FIRST,PRFLD_STEP).EQ.0) THEN
          IPOS=0                ! Point to start of file
          CALL SETPOS_SINGLE(NDEV_FLD,IPOS,ICODE)
          CALL BUFFIN_SINGLE(NDEV_FLD,WORKPTR,LENTHP,LEN_IO,A_IO)
 
           CALL GATHER_FIELD(PSTAR,WORK_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of PSTAR'
           ENDIF
           CALL GATHER_FIELD(WORKPTR,WORKPTR_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of PSTAR work'
           ENDIF

           IF(mype.eq.gather_pe) THEN
              CALL FLDDIAG(WORK_FULL,WORKPTR_FULL,KSTEP,
     &            glsize(1)*glsize(2),     1,' PSTAR ')
           ENDIF    ! test on gather PE
          DO 10 LEV=1,NUP
           CALL BUFFIN_SINGLE(NDEV_FLD,WORKPTR,LENTHP,LEN_IO,A_IO)
           CALL GATHER_FIELD(T(1,LEV),WORK_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of T'
           ENDIF
           CALL GATHER_FIELD(WORKPTR,WORKPTR_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of T work'
           ENDIF

           IF(mype.eq.gather_pe) THEN
              CALL FLDDIAG(WORK_FULL,WORKPTR_FULL,KSTEP,
     &            glsize(1)*glsize(2),   LEV,' T     ')
           ENDIF    ! test on gather PE
10        CONTINUE
          DO 11 LEV=1,NRP
           CALL BUFFIN_SINGLE(NDEV_FLD,WORKPTR,LENTHP,LEN_IO,A_IO) 
           CALL GATHER_FIELD(RH(1,LEV),WORK_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of RH'
           ENDIF
           CALL GATHER_FIELD(WORKPTR,WORKPTR_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2),
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of RH work'
           ENDIF

           IF(mype.eq.gather_pe) THEN
              CALL FLDDIAG(WORK_FULL,WORKPTR_FULL,KSTEP,
     &            glsize(1)*glsize(2),   LEV,' RH    ')
           ENDIF    ! test on gather PE
11        CONTINUE
          DO 12 LEV=1,NUP
           CALL BUFFIN_SINGLE(NDEV_FLD,WORKUV,LENUVP,LEN_IO,A_IO)
           CALL GATHER_FIELD(U(1,LEV),WORK_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2)-1,
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of U'
           ENDIF
           CALL GATHER_FIELD(WORKUV ,WORKPTR_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2)-1,
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of U work'
           ENDIF

           IF(mype.eq.gather_pe) THEN
              CALL FLDDIAG(WORK_FULL,WORKPTR_FULL,KSTEP,
     &            glsize(1)*(glsize(2)-1),   LEV,' U     ')
           ENDIF    ! test on gather PE
12        CONTINUE
          DO 13 LEV=1,NUP
           CALL BUFFIN_SINGLE(NDEV_FLD,WORKUV,LENUVP,LEN_IO,A_IO)
           CALL GATHER_FIELD(V(1,LEV),WORK_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2)-1,
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of V'
           ENDIF
           CALL GATHER_FIELD(WORKUV ,WORKPTR_FULL,
     &          lasize(1),lasize(2),glsize(1),glsize(2)-1,
     &          gather_pe,GC_ALL_PROC_GROUP,info)
           IF(info.NE.0) THEN      ! Check return code
              write(6,*) 'FLDSTAT1: Error in GATHER_FIELD of V work'
           ENDIF

           IF(mype.eq.gather_pe) THEN
              CALL FLDDIAG(WORK_FULL,WORKPTR_FULL,KSTEP,
     &            glsize(1)*(glsize(2)-1),   LEV,' V     ')
           ENDIF    ! test on gather PE
13        CONTINUE
       ENDIF
      ENDIF
CL
CL CLOSE DEVICE IF LAST TIMESTEP FOR DIAGNOSTIC (AND DELETE) 
CL
      IF(KSTEP.EQ.PRFLD_LAST) THEN
C   Close and delete explicit file name
         CALL CLOSE_SINGLE(NDEV_FLD,FLD_FILENAME,
     &                    LEN_FLD_FILENAME,1,1,ICODE)
      ELSE
CL
CL SAVE THIS TIMESTEP TO TMP DISK FILE USING (UNIT NDEV_FLD)
CL
         IPOS=0                ! Point to start of file
         CALL SETPOS_SINGLE(NDEV_FLD,IPOS,ICODE)
         CALL BUFFOUT_SINGLE(NDEV_FLD,PSTAR,LENTHP,LEN_IO,A_IO)
 
         DO LEV=1,NUP
          CALL BUFFOUT_SINGLE(NDEV_FLD,T (1,LEV),LENTHP,LEN_IO,A_IO)
         ENDDO ! LEV
         DO LEV=1,NRP 
          CALL BUFFOUT_SINGLE(NDEV_FLD,RH (1,LEV),LENTHP,LEN_IO,A_IO)
         ENDDO ! LEV  
         DO LEV=1,NUP
          CALL BUFFOUT_SINGLE(NDEV_FLD,U (1,LEV),LENUVP,LEN_IO,A_IO)
         ENDDO ! LEV  
         DO LEV=1,NUP  
          CALL BUFFOUT_SINGLE(NDEV_FLD,V (1,LEV),LENUVP,LEN_IO,A_IO)
         ENDDO ! LEV
      ENDIF

      RETURN
      END
      SUBROUTINE FLDDIAG(THIS,LAST,KSTEP,LENP,LEV,TITLE)
C
C CALC MAX MIN MEAN OF FIELD 'THIS'
C &    MAX MIN MEAN AND RMS OF FIELD 'THIS' MINUS 'LAST'
C NO AREA WEIGHTING IS APPLIED
C LOCATION OF MAX/MIN IS ALSO PRINTED

CLL Modification
CLL vn3.3  22/11/93 : Arrays THIS and LAST were declared before LENP(N.F
      IMPLICIT NONE

C
C*L  ARGUMENTS:---------------------------------------------------

      INTEGER
     +        KSTEP,            ! (IN) CURRENT TIMESTEP NO.
     +        LENP,             ! (IN) FIELD LENGTH
     +        LEV               ! (IN) MODEL LEVEL

      REAL
     +     THIS(LENP),          ! (IN) CURRENT FIELD
     +     LAST(LENP)           ! (IN) PREVIOUS FIELD
      CHARACTER*6 TITLE         ! (IN) FIELD TITLE

C
C DYNAMIC SPACE
C
      REAL
     +     DIFF(LENP),
     +     AMAX,AMIN,DMAX,DMIN,AMEAN,DMEAN,DRMS

      INTEGER
     +        IPT,                     ! POINT COUNTER
     +        IAMAX,IAMIN,IDMAX,IDMIN  ! FIELD MAX, MIN NO.

      AMAX=THIS(1)
      IAMAX=0
      DO 10 IPT=2,LENP
      IF(THIS(IPT).GT.AMAX)THEN
        AMAX=THIS(IPT)
        IAMAX=IPT
      ENDIF
10    CONTINUE

      AMIN=THIS(1)
      IAMIN=0
      DO 11 IPT=2,LENP
      IF(THIS(IPT).LT.AMIN)THEN
        AMIN=THIS(IPT)
        IAMIN=IPT
      ENDIF
11    CONTINUE

      DO 12 IPT=1,LENP
      DIFF(IPT)=THIS(IPT)-LAST(IPT)
12    CONTINUE

      DMAX=DIFF(1)
      IDMAX=0
      DO 13 IPT=2,LENP
      IF(DIFF(IPT).GT.DMAX)THEN
        DMAX=DIFF(IPT)
        IDMAX=IPT
      ENDIF
13    CONTINUE

      DMIN=DIFF(1)
      IDMIN=0
      DO 14 IPT=2,LENP
      IF(DIFF(IPT).LT.DMIN)THEN
        DMIN=DIFF(IPT)
        IDMIN=IPT
      ENDIF
14    CONTINUE

      AMEAN=THIS(1)
      DO 15 IPT=2,LENP
      AMEAN=THIS(IPT)+AMEAN
15    CONTINUE
      AMEAN=AMEAN/LENP

      DMEAN=DIFF(1)
      DO 16 IPT=2,LENP
      DMEAN=DIFF(IPT)+DMEAN
16    CONTINUE
      DMEAN=DMEAN/LENP

      DO 17 IPT=1,LENP
      DIFF(IPT)=DIFF(IPT)*DIFF(IPT)
17    CONTINUE

      DRMS=DIFF(1)
      DO 18 IPT=2,LENP
      DRMS=DIFF(IPT)+DRMS
18    CONTINUE
      DRMS=DRMS/LENP
      DRMS=SQRT(DRMS)
      IF(TITLE.EQ.' PSTAR ')THEN
C CONVERT TO MB
      AMAX=AMAX*.01
      AMIN=AMIN*.01
      AMEAN=AMEAN*.01
      DMAX=DMAX*.01
      DMIN=DMIN*.01
      DMEAN=DMEAN*.01
      DRMS=DRMS*.01
      WRITE(6,*)' STEP TITLE   LEV ',
     *          ' AMAX   IAMAX  AMIN    IAMIN  AMEAN  ', 
     *          ' DMAX   IDMAX    DMIN   IDMIN   DMEAN  DRMS '
      ENDIF
      WRITE(6,60)KSTEP,TITLE,LEV,
     *           AMAX,IAMAX,AMIN,IAMIN,AMEAN,
     *           DMAX,IDMAX,DMIN,IDMIN,DMEAN,DRMS
60    FORMAT(1X,I4,1X,A6,1X,I4,1X,
     *       F6.1,1X,I7,1X,F6.1,1X,I7,1X,F6.1,1X,
     *       F6.2,1X,I7,1X,F7.2,1X,I7,1X,F6.2,1X,F6.2)
      RETURN
      END
