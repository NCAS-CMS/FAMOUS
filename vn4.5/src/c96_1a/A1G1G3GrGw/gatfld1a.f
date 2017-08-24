C *****************************COPYRIGHT******************************
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
!+ Gathers a field from many processors to one processor
!
! Subroutine Interface:
      SUBROUTINE GATHER_FIELD(LOCAL_FIELD,GLOBAL_FIELD,
     &                        LOCAL_ROW_LEN,LOCAL_ROWS,
     &                        GLOBAL_ROW_LEN,GLOBAL_ROWS,
     &                        GATHER_PE,PROC_GROUP,
     &                        INFO)

      IMPLICIT NONE

!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    22/1/96   New DECK created for the Parallel Unified
!                     Model. P.Burton
!    4.2    17/10/96  Modify send/receive maps and change args to
!                     alltoall for GCOM/GCG v1.1     P.Burton
!    4.4    06/08/97  Recalculate maps if decomposition has changed
!                                                       Paul Burton
!
! Subroutine Arguments:

      INTEGER
     &  LOCAL_ROW_LEN    ! IN length of rows in local part of field
     &, LOCAL_ROWS       ! IN number of rows in local part of field
     &, GLOBAL_ROW_LEN   ! IN length of rows in global field
     &, GLOBAL_ROWS      ! IN number of rows in global field
     &, GATHER_PE        ! IN processor to gather global field to
     &, PROC_GROUP       ! IN group ID of processors involved here
     &, INFO             ! OUT return code from comms

      REAL
     &  LOCAL_FIELD(LOCAL_ROW_LEN*LOCAL_ROWS)
!                        ! IN local part of field
     &, GLOBAL_FIELD(GLOBAL_ROW_LEN*GLOBAL_ROWS)
!                        ! OUT (on PE GATHER_PE) global field

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

! Local variables

      INTEGER
     &   send_map(7,1)
     &,  receive_map(7,MAXPROC)
     &,  n_mess_to_rec

      INTEGER
     &  old_GLOBAL_ROW_LEN    ! value on last call
     &, old_GLOBAL_ROWS       ! value on last call
     &, old_PROC_GROUP        ! value on last call
     &, old_GATHER_PE         ! value on last call
     &, old_DECOMP            ! value on last call

      SAVE send_map,receive_map,n_mess_to_rec,
     &     old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,
     &     old_GATHER_PE,old_DECOMP
      DATA old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,
     &     old_GATHER_PE,old_DECOMP
     &   / -1234, -1234, -1234, -1234, -1234/

      INTEGER
     &  fld_type
     &, iproc
     &, flag

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

      IF ((GLOBAL_ROW_LEN .NE. old_GLOBAL_ROW_LEN) .OR.
     &    (GLOBAL_ROWS    .NE. old_GLOBAL_ROWS   ) .OR.
     &    (PROC_GROUP     .NE. old_PROC_GROUP    ) .OR.
     &    (GATHER_PE     .NE. old_GATHER_PE    ) .OR.
     &    (current_decomp_type .NE. old_DECOMP  )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

! 1.0 Find the type of field (P or U) being done

        IF (GLOBAL_ROWS .EQ. glsize(2)) THEN
          fld_type=fld_type_p
        ELSEIF (GLOBAL_ROWS .EQ. glsize(2)-1) THEN
          fld_type=fld_type_u
        ELSE
          WRITE(6,*) 'Bad field type in GATHER_FIELD'
          info=-1
          GOTO 9999
        ENDIF


! 2.0 Set up send map

        send_map(S_DESTINATION_PE,1) = GATHER_PE
!       processor to send to

        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =
     &    Offy*LOCAL_ROW_LEN+1+Offx
!       first data to send

        IF (atbase) THEN
          IF (fld_type .EQ. fld_type_p) THEN
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=LOCAL_ROWS-2*Offy
!           number of rows

          ELSE
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=LOCAL_ROWS-2*Offy-1

!           One less row at the bottom of a U field
          ENDIF
        ELSE
          send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1) = LOCAL_ROWS-2*Offy
!         number of rows

        ENDIF
        send_map(S_STRIDE_IN_SEND_ARRAY,1) = LOCAL_ROW_LEN
!       stride between row starts

        send_map(S_ELEMENT_LENGTH,1) = LOCAL_ROW_LEN-2*Offx
!       length of local row minus halos

        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =
     &    datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
!       start position in global data of this local data

        send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN
!       stride between rows in global data


! 3.0 Set up the receive map (for PE GATHER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

        n_mess_to_rec=0

        IF (mype .EQ. GATHER_PE) THEN
          DO iproc=0,nproc-1
            receive_map(R_SOURCE_PE,iproc+1) = iproc

            receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =
     &        g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsize(1)

            IF (fld_type .EQ. fld_type_p) THEN
              receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =
     &          g_blsizep(2,iproc)

            ELSE
              receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =
     &          g_blsizeu(2,iproc)
            ENDIF
            receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) =
     &        GLOBAL_ROW_LEN

            receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizep(1,iproc)

            receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =
     &        Offy*g_lasize(1,iproc)+Offx+1

            receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) =
     &        g_lasize(1,iproc)

          ENDDO
          n_mess_to_rec=nproc
        ENDIF

        old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
        old_GLOBAL_ROWS=GLOBAL_ROWS
        old_PROC_GROUP=PROC_GROUP
        old_GATHER_PE=GATHER_PE
        old_DECOMP=current_decomp_type

      ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

      flag=0  ! This is currently ignored at GCG v1.1

      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! gather operation
      info=GC_NONE

      CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,1,
     &                    LOCAL_ROW_LEN*LOCAL_ROWS,
     &                    GLOBAL_FIELD,receive_map,n_mess_to_rec,
     &                    GLOBAL_ROW_LEN*GLOBAL_ROWS,
     &                    PROC_GROUP,flag,info)

 9999 CONTINUE

      RETURN
      END
