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
!+ Parallel UM: Perform data decomposition for ocean model
!
! Subroutine Interface:
      SUBROUTINE DECOMPOSE_OCEAN(global_row_len,global_n_rows,
     &                           tot_levels,
     &                           nproc_EW, nproc_NS,
     &                           local_row_len,local_n_rows,
     &                           l_ocyclic)
      IMPLICIT NONE
!
! Desciption:
! This routine currently performs a 1D North-South decomposition on the
! ocean model. nproc_EW is currently ignored.
! The decomposition strategy is much the same as the atmosphere's -
! First try and divide the rows equally between processors, and then
! distribute any left over rows to the processors, starting from the
! top.
!
! Method:
! The local data sizes are calculated and sotred in the COMMON block
! DECOMPDB. The boundary conditions are set (cyclic in East/West
! direction if *DEF,GLOBAL
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/8/96  New deck created for MPP ocean model.  P.Burton
!  4.3      17/02/97 Added initialisation of new ocean decompositon
!                    decomp_nowrap_ocean - which does not include
!                    the wrap around points on the ends of rows.
!                    This requires passing in the l_ocyclic variable
!                    to indicate if these points are present.
!                                                         P.Burton
!
! Subroutine arguments:

      INTEGER

     &  global_row_len  ! IN :  number of E-W points of entire model
     &, global_n_rows   ! IN :  number of N-S points of entire model
     &, tot_levels      ! IN :  total number of levels
     &, nproc_EW        ! IN :  number of processors to decompose E-W
     &, nproc_NS        ! IN :  number of processors to decompose N-S
     &, local_row_len   ! OUT : local number of E-W points
     &, local_n_rows    ! OUT : local number of N-S points
!                       ! local_row_len and local_n_rows include
!                       ! any halos

      LOGICAL

     &  l_ocyclic       ! IN : true if extra wrap points are present
!                       !      at the start/ends of rows

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
! DECOMPDB comdeck
!
! Description:
!
! DECOMPDB comdeck (Decomposition Database) contains information
! describing the various decompositions used by the MPP-UM
! The CHANGE_DECOMPOSITION subroutine can be used to select
! a particular decomposition (which copies the appropriate
! decomposition information into the PARVARS common block).
!
! Requires comdeck PARVARS to be *CALLed before it.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton

! Common blocks containing information about each decomposition
! (For description of variables see the PARVARS comdeck)

      INTEGER
     &  decomp_db_bound(Ndim_max,max_decomps)
     &, decomp_db_glsize(Ndim_max,max_decomps)
     &, decomp_db_gridsize(Ndim_max,max_decomps)
     &, decomp_db_g_lasize(Ndim_max,0:maxproc,max_decomps)
     &, decomp_db_g_blsizep(Ndim_max,0:maxproc,max_decomps)
     &, decomp_db_g_blsizeu(Ndim_max,0:maxproc,max_decomps)
     &, decomp_db_g_datastart(Ndim_max,0:maxproc,max_decomps)
     &, decomp_db_g_gridpos(Ndim_max,0:maxproc,max_decomps)
     &, decomp_db_halosize(Ndim_max,max_decomps)
     &, decomp_db_neighbour(4,max_decomps)
     &, decomp_db_first_comp_pe(max_decomps)
     &, decomp_db_last_comp_pe(max_decomps)
     &, decomp_db_nproc(max_decomps)
     &, decomp_db_gc_proc_row_group(max_decomps)
     &, decomp_db_gc_proc_col_group(max_decomps)
     &, decomp_db_gc_all_proc_group(max_decomps)

      LOGICAL
     &  decomp_db_set(max_decomps)  ! indicates if a decomposition
!                                   ! has been initialised

      COMMON /DECOMP_DATABASE/
     &  decomp_db_bound , decomp_db_glsize
     &, decomp_db_g_lasize , decomp_db_gridsize
     &, decomp_db_g_blsizep , decomp_db_g_blsizeu
     &, decomp_db_g_datastart , decomp_db_g_gridpos
     &, decomp_db_halosize , decomp_db_neighbour
     &, decomp_db_first_comp_pe , decomp_db_last_comp_pe
     &, decomp_db_nproc
     &, decomp_db_gc_proc_row_group , decomp_db_gc_proc_col_group
     &, decomp_db_gc_all_proc_group
     &, decomp_db_set

! End of DECOMPDB comdeck

! Local variables
      INTEGER iproc,irest,jrest,info,in_ocn_decomp,idim

! ------------------------------------------------------------------

! Halo Sizes

      decomp_db_halosize(1,decomp_standard_ocean) = 0
      decomp_db_halosize(2,decomp_standard_ocean) = 1
      decomp_db_halosize(3,decomp_standard_ocean) = 0

! Size of global data

      decomp_db_glsize(1,decomp_standard_ocean) = global_row_len
      decomp_db_glsize(2,decomp_standard_ocean) = global_n_rows
      decomp_db_glsize(3,decomp_standard_ocean) = tot_levels

! Make sure there's actually enough work for all the processors to do

      IF (nproc_NS .GT. global_n_rows) THEN
        IF (mype .EQ. 0) THEN
          WRITE(6,*) 'Warning : Ocean model has more processors than ',
     &               'rows. Reducing nproc_y to ',global_n_rows
        ENDIF
        nproc_NS=global_n_rows
      ENDIF

      decomp_db_nproc(decomp_standard_ocean)=nproc_NS
      decomp_db_first_comp_pe(decomp_standard_ocean) = 0
      decomp_db_last_comp_pe(decomp_standard_ocean) =
     &  decomp_db_nproc(decomp_standard_ocean)-1

! Set the size of the Logical Processor Grid (LPG)

      decomp_db_gridsize(1,decomp_standard_ocean) = nproc_EW  ! =1
      decomp_db_gridsize(2,decomp_standard_ocean) = nproc_NS
      decomp_db_gridsize(3,decomp_standard_ocean) = 1

! Calculate processor specific information.

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
!       ! Loop over all processors in this decomposition

! NB : Although the decomposition is currently only N-S, all
! the code is included to allow an E-W decomposition too.
! All that is required is to supply nproc_NS > 1.

! Calculate the position in the LPG:
        decomp_db_g_gridpos(3,iproc,decomp_standard_ocean) = 0
        decomp_db_g_gridpos(2,iproc,decomp_standard_ocean) =
     &    iproc / decomp_db_gridsize(1,decomp_standard_ocean)
        decomp_db_g_gridpos(1,iproc,decomp_standard_ocean) =
     &    iproc - decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*
     &            decomp_db_gridsize(1,decomp_standard_ocean)

! Calculate the local data sizes for processor iproc

! East-West decomposition

        decomp_db_g_blsizep(1,iproc,decomp_standard_ocean) =
     &    decomp_db_glsize(1,decomp_standard_ocean) /
     &    decomp_db_gridsize(1,decomp_standard_ocean)
        irest = decomp_db_glsize(1,decomp_standard_ocean)-
     &          decomp_db_g_blsizep(1,iproc,decomp_standard_ocean)*
     &          decomp_db_gridsize(1,decomp_standard_ocean)
        decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =
     &    decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)*
     &    decomp_db_g_blsizep(1,iproc,decomp_standard_ocean) + 1

        IF (decomp_db_g_gridpos(1,iproc,decomp_standard_ocean) .LT.
     &      irest) THEN
          decomp_db_g_blsizep(1,iproc,decomp_standard_ocean) =
     &      decomp_db_g_blsizep(1,iproc,decomp_standard_ocean)+1
          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +
     &      decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)
        ELSE
          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +
     &      irest
        ENDIF

        decomp_db_g_lasize(1,iproc,decomp_standard_ocean)=
     &    decomp_db_g_blsizep(1,iproc,decomp_standard_ocean) +
     &    2*decomp_db_halosize(1,decomp_standard_ocean)

! North-South decomposition

        decomp_db_g_blsizep(2,iproc,decomp_standard_ocean) =
     &    decomp_db_glsize(2,decomp_standard_ocean) /
     &    decomp_db_gridsize(2,decomp_standard_ocean)
        jrest = decomp_db_glsize(2,decomp_standard_ocean)-
     &          decomp_db_g_blsizep(2,iproc,decomp_standard_ocean)*
     &          decomp_db_gridsize(2,decomp_standard_ocean)
        decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =
     &    decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*
     &    decomp_db_g_blsizep(2,iproc,decomp_standard_ocean) + 1

        IF (decomp_db_g_gridpos(2,iproc,decomp_standard_ocean) .LT.
     &      jrest) THEN
          decomp_db_g_blsizep(2,iproc,decomp_standard_ocean) =
     &      decomp_db_g_blsizep(2,iproc,decomp_standard_ocean)+1
          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +
     &      decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)
        ELSE
          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +
     &      jrest
        ENDIF

        decomp_db_g_lasize(2,iproc,decomp_standard_ocean)=
     &    decomp_db_g_blsizep(2,iproc,decomp_standard_ocean) +
     &    2*decomp_db_halosize(2,decomp_standard_ocean)

! No decomposition in the vertical

        decomp_db_g_datastart(3,iproc,decomp_standard_ocean) = 1
        decomp_db_g_blsizep(3,iproc,decomp_standard_ocean) =
     &    tot_levels
        decomp_db_g_lasize(3,iproc,decomp_standard_ocean) =
     &    tot_levels

! One less U row at bottom

        decomp_db_g_blsizeu(1,iproc,decomp_standard_ocean) =
     &    decomp_db_g_blsizep(1,iproc,decomp_standard_ocean)
        IF (  decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)
     &  .EQ. (decomp_db_gridsize(2,decomp_standard_ocean)-1))
     &  THEN
          decomp_db_g_blsizeu(2,iproc,decomp_standard_ocean) =
     &    decomp_db_g_blsizep(2,iproc,decomp_standard_ocean) - 1
        ELSE
          decomp_db_g_blsizeu(2,iproc,decomp_standard_ocean) =
     &    decomp_db_g_blsizep(2,iproc,decomp_standard_ocean)
        ENDIF
        decomp_db_g_blsizeu(3,iproc,decomp_standard_ocean) =
     &    decomp_db_g_blsizep(3,iproc,decomp_standard_ocean)

      ENDDO  ! loop over processors

! Set up the boundary types

      decomp_db_bound(1,decomp_standard_ocean) = BC_STATIC
!       ! No East-West wrap around
      decomp_db_bound(2,decomp_standard_ocean) = BC_STATIC
!       ! No North-South wrap around
      decomp_db_bound(3,decomp_standard_ocean) = BC_STATIC
!       ! No vertical wrap around

! And set up the neighbour array

      CALL SET_NEIGHBOUR(
     &  decomp_standard_ocean)

! Set up the GCOM groups

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_standard_ocean) .EQ. 1)
     & THEN
       decomp_db_gc_proc_row_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,
     &    decomp_db_g_gridpos(2,mype,decomp_standard_ocean),
     &    info,
     &    decomp_db_gc_proc_row_group(decomp_standard_ocean))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_standard_ocean) .EQ. 1)
     & THEN
        decomp_db_gc_proc_col_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,
     &    decomp_db_g_gridpos(1,mype,decomp_standard_ocean),
     &    info,
     &    decomp_db_gc_proc_col_group(decomp_standard_ocean))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_standard_ocean) .EQ. nproc_max)
     & THEN
        decomp_db_gc_all_proc_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        IF ((mype .GE. decomp_db_first_comp_pe(decomp_standard_ocean))
     &    .AND.
     &     (mype .LE. decomp_db_last_comp_pe(decomp_standard_ocean) ))
     &  THEN
          in_ocn_decomp=1
        ELSE
          in_ocn_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_ocn_decomp,info,
     &    decomp_db_gc_all_proc_group(decomp_standard_ocean))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_standard_ocean)=.TRUE.

! Initialise decomp_nowrap_ocean which doesn't contain extra wrap
! points at start and end of each row
! Mostly it is a straight copy of the original ocean decomposition

      DO idim=1,Ndim_max
        decomp_db_bound(idim,decomp_nowrap_ocean)=
     &    decomp_db_bound(idim,decomp_standard_ocean)
        decomp_db_glsize(idim,decomp_nowrap_ocean)=
     &    decomp_db_glsize(idim,decomp_standard_ocean)
        decomp_db_gridsize(idim,decomp_nowrap_ocean)=
     &    decomp_db_gridsize(idim,decomp_standard_ocean)
        decomp_db_halosize(idim,decomp_nowrap_ocean)=
     &    decomp_db_halosize(idim,decomp_standard_ocean)
      ENDDO

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
        DO idim=1,Ndim_max
          decomp_db_g_lasize(idim,iproc,decomp_nowrap_ocean)=
     &      decomp_db_g_lasize(idim,iproc,decomp_standard_ocean)
          decomp_db_g_blsizep(idim,iproc,decomp_nowrap_ocean)=
     &      decomp_db_g_blsizep(idim,iproc,decomp_standard_ocean)
          decomp_db_g_blsizeu(idim,iproc,decomp_nowrap_ocean)=
     &      decomp_db_g_blsizeu(idim,iproc,decomp_standard_ocean)
          decomp_db_g_datastart(idim,iproc,decomp_nowrap_ocean)=
     &      decomp_db_g_datastart(idim,iproc,decomp_standard_ocean)
          decomp_db_g_gridpos(idim,iproc,decomp_nowrap_ocean)=
     &      decomp_db_g_gridpos(idim,iproc,decomp_standard_ocean)
        ENDDO
      ENDDO

      DO idim=1,4
        decomp_db_neighbour(idim,decomp_nowrap_ocean)=
     &    decomp_db_neighbour(idim,decomp_standard_ocean)
      ENDDO

      decomp_db_first_comp_pe(decomp_nowrap_ocean)=
     &  decomp_db_first_comp_pe(decomp_standard_ocean)
      decomp_db_last_comp_pe(decomp_nowrap_ocean)=
     &  decomp_db_last_comp_pe(decomp_standard_ocean)
      decomp_db_nproc(decomp_nowrap_ocean)=
     &  decomp_db_nproc(decomp_standard_ocean)
      decomp_db_gc_proc_row_group(decomp_nowrap_ocean)=
     &  decomp_db_gc_proc_row_group(decomp_standard_ocean)
      decomp_db_gc_proc_col_group(decomp_nowrap_ocean)=
     &  decomp_db_gc_proc_col_group(decomp_standard_ocean)
      decomp_db_gc_all_proc_group(decomp_nowrap_ocean)=
     &  decomp_db_gc_all_proc_group(decomp_standard_ocean)

      IF (l_ocyclic) THEN
! Make modifications to the decompositions to remove the point at
! the beginning and end of each row
        decomp_db_glsize(1,decomp_nowrap_ocean)=
     &    decomp_db_glsize(1,decomp_nowrap_ocean)-2

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),
     &         decomp_db_last_comp_pe(decomp_standard_ocean)

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)
     &        .EQ. 0) THEN  ! this processor at left of LPG

            decomp_db_g_lasize(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_lasize(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_blsizep(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_blsizep(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_blsizeu(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_blsizeu(1,iproc,decomp_nowrap_ocean)-1

          ELSE  ! processor not at left of LPG

            decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)-1

          ENDIF

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)
     &       .EQ. (decomp_db_gridsize(1,decomp_nowrap_ocean)-1) )
     &    THEN  ! this processor at right of LPG

            decomp_db_g_lasize(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_lasize(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_blsizep(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_blsizep(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_blsizeu(1,iproc,decomp_nowrap_ocean)=
     &        decomp_db_g_blsizeu(1,iproc,decomp_nowrap_ocean)-1

          ENDIF

        ENDDO

      ENDIF

! Finally, indicate this decomposition has been initialised

      decomp_db_set(decomp_nowrap_ocean)=.TRUE.
! And return the new horizontal dimensions

      local_row_len=decomp_db_g_lasize(1,mype,decomp_standard_ocean)
      local_n_rows=decomp_db_g_lasize(2,mype,decomp_standard_ocean)

      RETURN
      END

