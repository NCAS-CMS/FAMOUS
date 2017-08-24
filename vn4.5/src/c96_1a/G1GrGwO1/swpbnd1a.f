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
C
!+ Parallel UM: Updates halo areas
!
! Subroutine interface:
      SUBROUTINE SWAPBOUNDS(FIELD,X_SIZE,Y_SIZE,X_OFF,Y_OFF,N_LEVELS)

      IMPLICIT NONE
!
! Description:
! This routine fills the halo areas (of size X_OFF in the x dimension
! and Y_OFF in the y dimension) of the first N_LEVELS of the array
! FIELD with the appropriate data from adjacent processors.
! If *DEF,GLOBAL is set, a east-west wrap around of data will
! occur.
!
! Method:
! Data to be sent to adjacent processors is packed into a sending
! array (which is on a COMMON block for T3D shmem useage), and
! sent to the relevant processors using GC_RSEND. The data is
! received by GC_RRECV, and then unpacked into the relevant halo
! area.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    9/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + R.Skaalin
!    4.1    18/3/96  Changed message ids to avoid clashes.
!                    Buffer logic revised. Removed timer   P.Burton
!    4.2    28/11/96 Added comdeck AMAXSIZE - reqd for BUFFERS
!    4.2    17/10/96  Use SETOPT to force shmem comms type.  P.Burton
!    4.3    24/02/97  Replace EW comms with single send/receive pairs
!                     if less than 3 PEs in the EW direction.
!                                                            P.Burton
!    4.3    18/03/97  Update N & S superpolar rows with polar values,
!                     to stop drift to large -ve temps in these rows.
!                                                       R.T.H.Barnes.
!
! Subroutine Arguments:

      INTEGER
     &   X_SIZE       ! IN  : X dimension of field (inc. halos)
     &,  Y_SIZE       ! IN  : Y dimension of field (inc. halos)
     &,  X_OFF        ! IN  : X halo size
     &,  Y_OFF        ! IN  : Y halo size
     & , N_LEVELS     ! IN  : Number of levels to be swapped

      REAL FIELD(X_SIZE*Y_SIZE,N_LEVELS)
!                     ! IN/OUT : Field to take place in
!                     !          boundary data exchange.

! Parameters and Common blocks

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
!====================== COMDECK BUFFERS ========================
!    Description:
!       This COMDECK defines buffer space used in various
!    Communication routines. Buffer space must be held in
!    COMMON for use with CRAY SHMEM.
!    NB: Requires PARVARS to be *CALLed first for MAXBUF

!    NB: Requires AMAXSIZE to be *CALLed first for maxFieldSize
!

      REAL
     &  buf1(MaxHaloSize*P_LEVELS_MAX)
     &, buf2(MaxHaloSize*P_LEVELS_MAX)
     &, buf3(MaxHaloSize*P_LEVELS_MAX)
     &, buf4(MaxHaloSize*P_LEVELS_MAX)

      COMMON /Halo_Buffer_Common/
     &  buf1,buf2,buf3,buf4

! Local variables

      INTEGER i,j,k, ioff, isize, jsize, info

! ------------------------------------------------------------------
      CALL GC_SSYNC(nproc,info)

      IF (Y_OFF .GT. 0) THEN
!       Do North/South communication

!       Send to Northen neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PNorth) .NE. NoDomain) THEN
          ioff = X_SIZE*Y_OFF+X_OFF
          isize = (X_SIZE-(2*X_OFF))*N_LEVELS
          DO j = 1, Y_OFF
            DO k = 1, N_LEVELS
              DO i = 1, (X_SIZE-(2*X_OFF))
                buf1((j-1)*isize+(k-1)*(X_SIZE-(2*X_OFF))+i)
     &             = FIELD(ioff+i,k)
              ENDDO
            ENDDO
            ioff = ioff + X_SIZE
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(1001, isize*Y_OFF, neighbour(PNorth),
     &                  info, buf3, buf1)
        ELSE
! Copy values to Northmost row (outside domain)
          DO k = 1,N_LEVELS
            DO j = 1,Y_OFF
              DO i = 1,X_SIZE
                FIELD((j-1)*x_size+i,k) = FIELD(y_off*x_size+i,k)
              END DO
            END DO
          END DO
        ENDIF

!       Send to Southern neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PSouth) .NE. NoDomain) THEN
          ioff = X_SIZE*Y_SIZE-2*Y_OFF*X_SIZE+X_OFF
          isize = (X_SIZE-(2*X_OFF))*N_LEVELS
          DO j = 1, Y_OFF
            DO k = 1, N_LEVELS
              DO i = 1, (X_SIZE-(2*X_OFF))
                buf2((j-1)*isize+(k-1)*(X_SIZE-(2*X_OFF))+i) =
     &            FIELD(ioff+i,k)
              ENDDO
            ENDDO
            ioff = ioff + X_SIZE
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(2002, isize*Y_OFF, neighbour(PSouth),
     &                  info, buf4, buf2)
        ELSE
! Copy values to Southmost row (outside domain)
          DO k = 1,N_LEVELS
            DO j = 1,Y_OFF
              DO i = 1,X_SIZE
                FIELD((y_size-j)*x_size+i,k) =
     &            FIELD((y_size-y_off-1)*x_size+i,k)
              END DO
            END DO
          END DO
        ENDIF

!       Synchronize before receiving

        CALL GC_SSYNC(nproc,info)

!       Receive from Southern neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PSouth) .NE. NoDomain) THEN
          ioff = X_SIZE*(Y_SIZE-Y_OFF)+X_OFF
          isize = (X_SIZE-(2*X_OFF))*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(1001, isize*Y_OFF, neighbour(PSouth),
     &                  info, buf3, buf1)
          DO j = 1, Y_OFF
            DO k = 1, N_LEVELS
              DO i = 1, (X_SIZE-(2*X_OFF))
                FIELD(ioff+i,k) =
     &            buf3((j-1)*isize+(k-1)*(X_SIZE-(2*X_OFF))+i)
              ENDDO
            ENDDO
            ioff = ioff + X_SIZE
          ENDDO
        ENDIF

!       Receive from Northen neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PNorth) .NE. NoDomain) THEN
          ioff = X_OFF
          isize = (X_SIZE-(2*X_OFF))*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(2002, isize*Y_OFF, neighbour(PNorth),
     &                  info, buf4, buf2)
          DO j = 1, Y_OFF
            DO k = 1, N_LEVELS
              DO i = 1, (X_SIZE-(2*X_OFF))
                FIELD(ioff+i,k) =
     &            buf4((j-1)*isize+(k-1)*(X_SIZE-(2*X_OFF))+i)
              ENDDO
            ENDDO
            ioff = ioff + X_SIZE
          ENDDO
        ENDIF

      ENDIF  ! should we do north/south communication?

      IF (X_OFF .GT. 0) THEN

!       Do East/West communication

!       Send to Western neighbour
        IF (gridsize(1) .GT. 2) THEN
!         If there are more than two processors in the EW direction
!         it is safe to do two sends (E+W) before we receive our
!         halos - as the two processors we are sending to must
!         be different. (SHMEM_NAM only allows one outstanding
!         send to each processor).

!       A full column (Y_SIZE) is sent, so corner points are included


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PWest) .NE. NoDomain) THEN
          ioff = X_OFF
          jsize = Y_SIZE*N_LEVELS
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                buf3((i-1)*jsize+(k-1)*Y_SIZE+j) =
     &            FIELD((j-1)*X_SIZE+ioff+1,k)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(3003, jsize*X_OFF, neighbour(PWest),
     &                  info, buf1, buf3)
        ENDIF

!       Send to Eastern neighbour.
!       A full column (Y_SIZE) is sent, so corner points are included


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PEast) .NE. NoDomain) THEN
          ioff = X_SIZE - 2*X_OFF
          jsize = Y_SIZE*N_LEVELS
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                buf4((i-1)*jsize+(k-1)*Y_SIZE+j) =
     &            FIELD((j-1)*X_SIZE+ioff+1,k)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(4004, jsize*X_OFF, neighbour(PEast),
     &                  info, buf2, buf4)
        ENDIF

!       Synchronize before receiving

        CALL GC_SSYNC(NPROC,INFO)

!       Receive from Eastern neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PEast) .NE. NoDomain) THEN
          ioff = X_SIZE - X_OFF
          jsize = Y_SIZE*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(3003, jsize*X_OFF, neighbour(PEast),
     &                  info, buf1, buf3)
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                FIELD((j-1)*X_SIZE+ioff+1,k) =
     &            buf1((i-1)*jsize+(k-1)*Y_SIZE+j)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
        ENDIF

!       Receive from Western neighbour


          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PWest) .NE. NoDomain) THEN
          ioff = 0
          jsize = Y_SIZE*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(4004, jsize*X_OFF, neighbour(PWest),
     &    info, buf2, buf4)
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                FIELD((j-1)*X_SIZE+ioff+1,k) =
     &            buf2((i-1)*jsize+(k-1)*Y_SIZE+j)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
        ENDIF
        ELSE  ! One or two processors in EW direction
!         Do the same communication as before - but as two
!         seperate send/receive phases

!       Send to Western neighbour
!       A full column (Y_SIZE) is sent, so corner points are included


        CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PWest) .NE. NoDomain) THEN
          ioff = X_OFF
          jsize = Y_SIZE*N_LEVELS
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                buf3((i-1)*jsize+(k-1)*Y_SIZE+j) =
     &            FIELD((j-1)*X_SIZE+ioff+1,k)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(3003, jsize*X_OFF, neighbour(PWest),
     &                  info, buf1, buf3)
        ENDIF

!       Synchronize before receiving

        CALL GC_SSYNC(NPROC,INFO)

!       Receive from Eastern neighbour


        CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PEast) .NE. NoDomain) THEN
          ioff = X_SIZE - X_OFF
          jsize = Y_SIZE*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(3003, jsize*X_OFF, neighbour(PEast),
     &                  info, buf1, buf3)
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                FIELD((j-1)*X_SIZE+ioff+1,k) =
     &            buf1((i-1)*jsize+(k-1)*Y_SIZE+j)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
        ENDIF

!       Send to Eastern neighbour.
!       A full column (Y_SIZE) is sent, so corner points are included


        CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PEast) .NE. NoDomain) THEN
          ioff = X_SIZE - 2*X_OFF
          jsize = Y_SIZE*N_LEVELS
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                buf4((i-1)*jsize+(k-1)*Y_SIZE+j) =
     &            FIELD((j-1)*X_SIZE+ioff+1,k)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
          info=GC_NONE
          CALL GC_RSEND(4004, jsize*X_OFF, neighbour(PEast),
     &                  info, buf2, buf4)
        ENDIF

!       Synchronize before receiving

        CALL GC_SSYNC(NPROC,INFO)

!       Receive from Western neighbour

        CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! use shmem_put

        IF (neighbour(PWest) .NE. NoDomain) THEN
          ioff = 0
          jsize = Y_SIZE*N_LEVELS
          info=GC_NONE
          CALL GC_RRECV(4004, jsize*X_OFF, neighbour(PWest),
     &    info, buf2, buf4)
          DO i = 1, X_OFF
            DO k = 1, N_LEVELS
              DO j = 1, Y_SIZE
                FIELD((j-1)*X_SIZE+ioff+1,k) =
     &            buf2((i-1)*jsize+(k-1)*Y_SIZE+j)
              ENDDO
            ENDDO
            ioff = ioff + 1
          ENDDO
        ENDIF

        ENDIF ! check for number of PEs in EW direction

      ENDIF  ! should we be doing east/west communications ?

!     CALL GC_SSYNC(nproc,info)
      RETURN
      END
