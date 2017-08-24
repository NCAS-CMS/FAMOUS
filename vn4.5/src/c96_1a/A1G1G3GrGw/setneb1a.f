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
!+ Parallel UM: Sets up array of tids of adjacent processes.
!
! Subroutine interface:
      SUBROUTINE SET_NEIGHBOUR(decomp_type)

      IMPLICIT NONE

!
! Description:
! This routine finds the tids of the North, South, East and West
! neighbouring processes. It takes account of the boundary
! condition in each dimension (X and Y) which can be either:
! cyclic : wrap around
! static : no wrap around
!
! Method:
! The tid of each neighbouring process is calculated (taking into
! account the relevant boundary conditions) and placed in the
! neighbour array.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    4/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + R.Skaalin
!    4.2    21/08/96 Updated to use DECOMPDB variables rather than
!                    PARVARS. Only argument required is the
!                    decomposition type.              P.Burton

!
! Subroutine Arguments:

      INTEGER decomp_type  ! Decomposition type to update neighbours
!                          ! for

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

! ------------------------------------------------------------------


C Set Northen neighbour

      IF ( decomp_db_g_gridpos(2,mype,decomp_type) .GT. 0) THEN
!       This processor is not at the top of the LPG
        decomp_db_neighbour(PNorth,decomp_type) =
     &    mype - decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type) .EQ. BC_CYCLIC) THEN
!       This processor at the top of LPG, and has cyclic BCs.
        decomp_db_neighbour(PNorth,decomp_type) =
     &    mype + decomp_db_nproc(decomp_type) -
     &    decomp_db_gridsize(1,decomp_type)
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PNorth,decomp_type) =NoDomain
      ENDIF

C Set Southern neighbour

      IF ( decomp_db_g_gridpos(2,mype,decomp_type) .LT.
     &    (decomp_db_gridsize(2,decomp_type)-1) ) THEN
!       This processor is not at the bottom of the LPG
        decomp_db_neighbour(PSouth,decomp_type) =
     &    mype + decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type) .EQ. BC_CYCLIC) THEN
!       This processor at the bottom of LPG, and has cyclic BCs.
        decomp_db_neighbour(PSouth,decomp_type) =
     &    mype - decomp_db_nproc(decomp_type) +
     &    decomp_db_gridsize(1,decomp_type)
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PSouth,decomp_type) =NoDomain
      ENDIF

C Set Western neighbour

      IF ( decomp_db_g_gridpos(1,mype,decomp_type) .GT. 0) THEN
!       This processor is not at the left of the LPG
        decomp_db_neighbour(PWest,decomp_type) =
     &    mype - 1
      ELSEIF (decomp_db_bound(1,decomp_type) .EQ. BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
        decomp_db_neighbour(PWest,decomp_type) =
     &    mype + decomp_db_gridsize(1,decomp_type) - 1
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PWest,decomp_type) =NoDomain
      ENDIF

C Set Eastern neighbour
      IF ( decomp_db_g_gridpos(1,mype,decomp_type) .LT.
     &    (decomp_db_gridsize(1,decomp_type)-1) ) THEN
!       This processor is not at the right of the LPG
        decomp_db_neighbour(PEast,decomp_type) =
     &    mype + 1
      ELSEIF (decomp_db_bound(1,decomp_type) .EQ. BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
        decomp_db_neighbour(PEast,decomp_type) =
     &    mype - decomp_db_gridsize(1,decomp_type) + 1
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PEast,decomp_type) =NoDomain
      ENDIF

      RETURN
      END
