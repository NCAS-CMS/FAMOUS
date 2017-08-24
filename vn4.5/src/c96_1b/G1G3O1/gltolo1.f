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
!+ Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_SUBDOMAIN: converts global subdomain boundaries
!                            to local subdomain boundaries
! GLOBAL_TO_LOCAL_RC: converts global row,column co-ordinates to
!                     processor co-ordinates plus local
!                     co-ordinates within the processor.
!
! Subroutine Interface:
      SUBROUTINE GLOBAL_TO_LOCAL_SUBDOMAIN(
     &                              L_include_halosEW,
     &                              L_include_halosNS,
     &                              grid_code,procid,
     &                              global_north_in,global_east_in,
     &                              global_south_in,global_west_in,
     &                              local_north,local_east,
     &                              local_south,local_west)
      IMPLICIT NONE
!
! Description:
! Takes a global definition of a subdomain region (in terms of
! model gridpoints) and translates it into local numbers.
! This effectively means local co-ordinates of the region of the
! subdomain which intersects with this processor's area.
!
! Method:
! Use the datastart variable in PARVARS to see if the requested
! subdomain intersects with this processor's area, if it does
! then use datastart to convert to local co-ordinate and do a bit
! of logic using MAX and MIN to ensure the local co-ordinates
! actually lie within the local area  Then make any corrections
! necessary to account for a subdomain which crosses over the
! 0 longitude line. Finally, if L_include_halos is set to
! .TRUE. - include any relevant halo regions.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      03/09/96 New deck created for MPP code.  P.Burton
!  4.3      13/03/97 Various bug fixes               P.Burton
!  4.4      12/06/97 Another bug fix                 P.Burton
!
! Subroutine arguments:

      LOGICAL
     &  L_include_halosEW  ! IN : include East-West halos in local
!                          !      region if set to .TRUE.
     &, L_include_halosNS  ! IN : include North-South halos in local
!                          !      region if set to .TRUE.
      INTEGER

     &  grid_code        ! IN : STASH grid type of field
     &, procid           ! IN : processor to produce result for
     &, global_north_in  ! IN : northern boundary of global subdomain
     &, global_east_in   ! IN : eastern boundary of global subdomain
     &, global_south_in  ! IN : southern boundary of global subdomain
     &, global_west_in   ! IN : western boundary of global subdomain

     &, local_north   ! OUT : northern boundary of local subdomain
     &, local_east    ! OUT : eastern boundary of local subdomain
     &, local_south   ! OUT : southern boundary of local subdomain
     &, local_west    ! OUT : western boundary of local subdomain

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
CLL  Comdeck: STERR ----------------------------------------------------
CLL
CLL  Purpose: PARAMETER names for STASH processing error codes;
CLL           fatal errors have positive codes, warnings negative.
CLL
CLL  Author:   S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for MPP code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D70
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
C
C Warning codes
C
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
C
C Error codes
C
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)


! Local variables
      INTEGER
! Copies of the input arguments, that can be modified for
! wrap-around calculations
     &  global_north,global_east,global_south,global_west
     &, fld_type  ! is field on P or U grid?
     &, row_len_nh    ! row length when halos are removed
     &, nrows_nh      ! number of rows when halos are removed
     &, first_global_pt_EW ! global point number of first and last
     &, last_global_pt_EW  ! local points in local area
     &, first_global_pt_NS ! in the East-West and
     &, last_global_pt_NS  ! North-South directions

      LOGICAL
! Logicals indicating if this processor contains part of a
! subdomain
     &  NS_intersect,EW_intersect
     &, wrap ! set to .TRUE. if the subdomain passes over the
!            ! the 0 degree longitude line
     &, fullfield ! if the field is NOT a subdomain

      INTEGER GET_FLD_TYPE  ! function

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

      global_north=global_north_in
      global_east=global_east_in
      global_south=global_south_in
      global_west=global_west_in

! Find out if the data is on a mass or velocity grid

      fld_type=GET_FLD_TYPE(grid_code)

      IF (fld_type .EQ. fld_type_unknown) THEN
        WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',
     &    'field with gridtype code ',grid_code
        WRITE(6,*) 'Unable to process this field.'
        local_north=st_no_data
        local_south=st_no_data
        local_east=st_no_data
        local_west=st_no_data
        GOTO 9999
      ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

      fullfield= ((( global_west .EQ. 1 ) .AND.
     &             ( global_east .EQ. glsize(1)) .AND.
     &             ( global_north .EQ. 1 )) .AND.
     &            (((fld_type .EQ. fld_type_p) .AND.
     &              (global_south .EQ. glsize(2))) .OR.
     &             ((fld_type .EQ. fld_type_u) .AND.
     &              (global_south .EQ. glsize(2)-1))))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

      IF (fullfield) THEN

        IF (L_include_halosNS) THEN
          local_north=1
          local_south=g_lasize(2,procid)
        ELSE
          local_north=1+Offy
          local_south=g_lasize(2,procid)-Offy
        ENDIF
        IF (L_include_halosEW) THEN
          local_west=1
          local_east=g_lasize(1,procid)
        ELSE
          local_west=1+Offx
          local_east=g_lasize(1,procid)-Offx
        ENDIF

      ELSE ! a subdomain requires some careful analysis:

        row_len_nh=g_blsizep(1,procid)
        IF (fld_type .EQ. fld_type_p) THEN
          nrows_nh=g_blsizep(2,procid)
        ELSE
          nrows_nh=g_blsizeu(2,procid)
        ENDIF

! Set up variables giving the global point numbers of the
! start and end of this processor's subdomain

        first_global_pt_EW=g_datastart(1,procid)
        last_global_pt_EW=first_global_pt_EW+row_len_nh-1

        first_global_pt_NS=g_datastart(2,procid)
        last_global_pt_NS=first_global_pt_NS+nrows_nh-1

! If global_east is greater than the global row length, this
! indicates a wrap around - but not in the format this code
! wants - where it expects a wrap around to be indicated by
! the east column being less than the west column.

        IF (global_east .LT. global_west) THEN
          wrap=.TRUE.
        ELSEIF (global_east .GT. glsize(1)) THEN
          wrap=.TRUE.
          global_east=global_east-glsize(1)
        ELSE
          wrap=.FALSE.
        ENDIF

        EW_intersect =
     &    (( .NOT. wrap) .AND.
     &     ((global_east .GE. first_global_pt_EW) .AND.
     &      (global_west .LE. last_global_pt_EW)))
     &    .OR.
     &    ((wrap) .AND.
     &     ((global_west .LE. last_global_pt_EW) .OR.
     &      (global_east .GE. first_global_pt_EW)))

        NS_intersect =
     &    ((global_south .GE. first_global_pt_NS) .AND.
     &     (global_north .LE. last_global_pt_NS))

        IF (NS_intersect) THEN

          IF ((global_north .GE. first_global_pt_NS) .AND.
     &        (global_north .LE. last_global_pt_NS)) THEN
! This processor contains the NS start of the subarea
            local_north=global_north-first_global_pt_NS+Offy+1
          ELSE
! This processor is below the start of the subarea
            local_north=1+Offy
          ENDIF

          IF ((global_south .GE. first_global_pt_NS) .AND.
     &        (global_south .LE. last_global_pt_NS)) THEN
! This processor contains the NS end of the subarea
            local_south=global_south-first_global_pt_NS+Offy+1
          ELSE
! This processor is above the end of the subarea
            local_south=Offy+nrows_nh
          ENDIF

        ELSE

          local_north=st_no_data
          local_south=st_no_data

        ENDIF

        IF (EW_intersect) THEN

          IF ((global_west .GE. first_global_pt_EW) .AND.
     &        (global_west .LE. last_global_pt_EW)) THEN
! This processor contains the EW start of the subarea
            local_west=global_west-first_global_pt_EW+Offx+1
          ELSE
! This processor is to the right of the start of the subarea
            local_west=1+Offx
          ENDIF

          IF ((global_east .GE. first_global_pt_EW) .AND.
     &        (global_east .LE. last_global_pt_EW)) THEN
! This processor contains the EW end of the subarea
            local_east=global_east-first_global_pt_EW+Offx+1
          ELSE
! This processor is to the left of the end of the subarea
            local_east=Offx+row_len_nh
          ENDIF

        ELSE

          local_east=st_no_data
          local_west=st_no_data

        ENDIF

      ENDIF ! is this a fullfield?

 9999 CONTINUE

      RETURN
      END

! Subroutine Interface:
      SUBROUTINE GLOBAL_TO_LOCAL_RC(grid_code,
     &                              global_column_in , global_row,
     &                              processor_x , processor_y,
     &                              local_column, local_row)

      IMPLICIT NONE
!
! Description:
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.
!
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      17 /09/96 New deck created for MPP code.  P.Burton
!  4.3      13/03/97  Various bug fixes               P.Burton
!  4.4      18/06/97  Check that row number is valid      P.Burton
!           06/10/97  Set correct row length and n_rows
!                     in dowhile loop.                    P.Burton
!
! Subroutine arguments:

      INTEGER
     &  grid_code          ! IN : STASH grid type code
     &, global_column_in   ! IN : global column number
     &, global_row         ! IN : global row number
     &, processor_x        ! OUT : processor X (EW) co-ordinate
!                          !       (0->nproc_x)
     &, processor_y        ! OUT : processor Y (NS) co-ordinate
!                               (0->nproc_y)
     &, local_column       ! OUT : local column number on processor
     &, local_row          ! OUT : local row number on processor

! Parameters and COMMON blocks
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
CLL  Comdeck: STERR ----------------------------------------------------
CLL
CLL  Purpose: PARAMETER names for STASH processing error codes;
CLL           fatal errors have positive codes, warnings negative.
CLL
CLL  Author:   S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for MPP code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D70
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
C
C Warning codes
C
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
C
C Error codes
C
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)


! Local variables

      INTEGER
     &  global_column ! modified version of global_column_in which
!                     ! takes account of domains wrapping over
!                     ! 0 degree longitude
     &, fld_type      ! field stored on P grid or U grid?
     &, row_len_nh,nrows_nh  ! row_len and n_rows when halos removed
     &, proc  ! loop counter for loop over processors
! global column and row numbers delimiting a processors area
     &, start_col,end_col,start_row,end_row

      INTEGER GET_FLD_TYPE  ! function

! ------------------------------------------------------------------

! Find out if the data is on a mass or velocity grid

      fld_type=GET_FLD_TYPE(grid_code)

      IF (fld_type .EQ. fld_type_unknown) THEN
        WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',
     &    'field with gridtype code ',grid_code
        WRITE(6,*) 'Unable to process this field.'
        processor_x=st_no_data
        processor_y=st_no_data
        local_column=st_no_data
        local_row=st_no_data
        GOTO 9999
      ENDIF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

      IF (global_column_in .GT. glsize(1)) THEN
        global_column=MOD(global_column_in+1,glsize(1))-1
      ELSE
        global_column=global_column_in
      ENDIF

      IF ((global_column .LT. 1) .OR.
     &    (global_row .LT. 1) .OR.
     &    (global_row .GT. glsize(2))) THEN

        WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',
     &  'impossible global row/column co-ordinates ',
     &  'row: ',global_row,' column: ',global_column

        processor_x=st_no_data
        processor_y=st_no_data
        local_column=st_no_data
        local_row=st_no_data

      ENDIF

! Make a first guess at the processor co-ordinates

      processor_x=MIN(global_column/(glsize(1)/gridsize(1)),
     &                nproc_x-1)
      processor_y=MIN(global_row/(glsize(2)/gridsize(2)),
     &                nproc_y-1)

      proc=processor_x+processor_y*gridsize(1)

      row_len_nh=g_blsizep(1,proc)
      IF (fld_type .EQ. fld_type_p) THEN
        nrows_nh=g_blsizep(2,proc)
      ELSE
        nrows_nh=g_blsizeu(2,proc)
      ENDIF

      start_col=g_datastart(1,proc)
      end_col=start_col+row_len_nh-1
      start_row=g_datastart(2,proc)
      end_row=start_row+nrows_nh-1

! Now iterate around these processors until we hit the right one

      DO WHILE
     &    (((global_column .LT. start_col) .OR.
     &      (global_column .GT. end_col  ))
     &   .OR.
     &     ((global_row .LT. start_row) .OR.
     &      (global_row .GT. end_row)))


        IF (global_column .LT. start_col) THEN
          processor_x=processor_x-1
        ELSEIF (global_column .GT. end_col) THEN
          processor_x=processor_x+1
        ENDIF

        IF (global_row .LT. start_row) THEN
          processor_y=processor_y-1
        ELSEIF (global_row .GT. end_row) THEN
          processor_y=processor_y+1
        ENDIF

        proc=processor_x+processor_y*gridsize(1)

      row_len_nh=g_blsizep(1,proc)
      IF (fld_type .EQ. fld_type_p) THEN
        nrows_nh=g_blsizep(2,proc)
      ELSE
        nrows_nh=g_blsizeu(2,proc)
      ENDIF
        start_col=g_datastart(1,proc)
        end_col=start_col+row_len_nh-1
        start_row=g_datastart(2,proc)
        end_row=start_row+nrows_nh-1

      ENDDO

! Now we have the processor co-ordinates, we can calculate the
! local co-ordinates.

      local_column=Offx+global_column-start_col+1
      local_row=Offy+global_row-start_row+1

 9999 CONTINUE

      RETURN
      END

! Function Interface
      INTEGER FUNCTION GET_FLD_TYPE (grid_type_code)

      IMPLICIT NONE

!
! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/11/96 New deck created for MPP code.  P.Burton
!
! Subroutine arguments:

      INTEGER
     &  grid_type_code     ! IN : STASH grid type code

! Parameters
CLL  Comdeck: CPPXREF --------------------------------------------------
CLL
CLL  Purpose: Holds PARAMETERs describing structure of PP_XREF file,
CLL           and some values for valid entries.
CLL
CLL  Author    Dr T Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
CLL                  1.Removes the limit on primary STASH item numbers.
CLL                  2.Removes the assumption that (section,item)
CLL                    defines the sub-model.
CLL                  3.Thus allows for user-prognostics.
CLL                  Add a PPXREF record for model number.
CLL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
CLL  3.5   10/3/94   Sub-Models project:
CLL                 List of PPXREF addressing codes augmented, in order
CLL                 to include all of the pre_STASH master information
CLL                 in the new PPXREF file.
CLL                 PPXREF_CODELEN increased to 38.
CLL                 PPXREF_IDLEN deleted - no longer relevant.
CLL                   S.J.Swarbrick
CLL  4.1   June 96  Wave model parameters included.
CLL                 ppx_ address parameters adjusted to allow for 
CLL                  reading option code as 4x5 digit groups.
CLL                   S.J.Swarbrick  
CLL
CLL  Logical components covered: C40
CLL
C-----------------------------------------------------------------------
C Primary file record definition
      INTEGER
     *       PPXREF_IDLEN,PPXREF_CHARLEN,PPXREF_CODELEN
     *      ,PPXREF_PACK_PROFS
      PARAMETER(
     *       PPXREF_IDLEN=2,               ! length of id in a record
! WARNING: PPXREF_CHARLEN must be an exact multiple of 4
!                         to avoid overwriting
     *       PPXREF_CHARLEN=36,            ! total length of characters
     *       PPXREF_PACK_PROFS=10,         ! number of packing profiles
     *       PPXREF_CODELEN=40)            ! total length of codes      
C Derived file record sizes
      INTEGER
     *       PPX_CHARWORD,PPX_RECORDLEN
      PARAMETER(
C            Assume that an integer is at least 4 bytes long.
C            This wastes some space and memory on 8 byte machines.
     *       PPX_CHARWORD=((PPXREF_CHARLEN+3)/4), ! i.e., ppx_charword=9
     *       PPX_RECORDLEN=
     *           PPX_CHARWORD+PPXREF_CODELEN)  ! read buffer record len
C
C-----------------------------------------------------------------------
C Addressing codes within PPXREF
      INTEGER
     &       ppx_model_number   ,
     &       ppx_section_number ,ppx_item_number    ,
     &       ppx_version_mask   ,ppx_space_code     ,
     &       ppx_timavail_code  ,ppx_grid_type      ,
     &       ppx_lv_code        ,ppx_lb_code        ,
     &       ppx_lt_code        ,ppx_lev_flag       ,
     &       ppx_opt_code       ,ppx_pt_code        ,
     &       ppx_pf_code        ,ppx_pl_code        ,
     &       ppx_ptr_code       ,ppx_lbvc_code      ,
     &       ppx_dump_packing   ,ppx_rotate_code    ,
     &       ppx_field_code     ,ppx_user_code      ,
     &       ppx_meto8_levelcode,ppx_meto8_fieldcode,
     &       ppx_cf_levelcode   ,ppx_cf_fieldcode   ,
     &       ppx_base_level     ,ppx_top_level      ,
     &       ppx_ref_LBVC_code  ,ppx_data_type      ,
     &       ppx_packing_acc    ,ppx_pack_acc
      PARAMETER(
     &       ppx_model_number   = 1,  ! Model number address
     &       ppx_section_number = 2,  ! Section number address
     &       ppx_item_number    = 3,  ! Item number address
     &       ppx_version_mask   = 4,  ! Version mask address
     &       ppx_space_code     = 5,  ! Space code address
     &       ppx_timavail_code  = 6,  ! Time availability code address
     &       ppx_grid_type      = 7,  ! Grid type code address
     &       ppx_lv_code        = 8,  ! Level type code address
     &       ppx_lb_code        = 9,  ! First level code address
     &       ppx_lt_code        =10,  ! Last level code address
     &       ppx_lev_flag       =11,  ! Level compression flag address
     &       ppx_opt_code       =12,  ! Sectional option code address
     &       ppx_pt_code        =16,  ! Pseudo dimension type address   
     &       ppx_pf_code        =17,  ! First pseudo dim code address   
     &       ppx_pl_code        =18,  ! Last pseudo dim code address    
     &       ppx_ptr_code       =19,  ! Section 0 point-back code addres
     &       ppx_dump_packing   =20,  ! Dump packing code address       
     &       ppx_lbvc_code      =21,  ! PP LBVC code address            
     &       ppx_rotate_code    =22,  ! Rotation code address           
     &       ppx_field_code     =23,  ! PP field code address           
     &       ppx_user_code      =24,  ! User code address               
     &       ppx_meto8_levelcode=25,  ! CF level code address           
     &       ppx_meto8_fieldcode=26,  ! CF field code address           
     &       ppx_cf_levelcode   =25,                                    
     &       ppx_cf_fieldcode   =26,                                    
     &       ppx_base_level     =27,  ! Base level code address         
     &       ppx_top_level      =28,  ! Top level code address          
     &       ppx_ref_lbvc_code  =29,  ! Ref level LBVC code address     
     &       ppx_data_type      =30,  ! Data type code address          
     &       ppx_packing_acc    =31,  ! Packing accuracy code address (1
     &       ppx_pack_acc       =31)                                    
C
C Valid grid type codes
      INTEGER
     &       ppx_atm_nonstd,ppx_atm_tall,ppx_atm_tland,ppx_atm_tsea,
     &       ppx_atm_uall,ppx_atm_uland,ppx_atm_usea,ppx_atm_compressed,
     &       ppx_atm_ozone,ppx_atm_tzonal,ppx_atm_uzonal,ppx_atm_rim,
     &       ppx_atm_tmerid,ppx_atm_umerid,ppx_atm_scalar,
     &       ppx_atm_cuall,ppx_atm_cvall,
     &       ppx_ocn_nonstd,ppx_ocn_tall,ppx_ocn_tcomp,ppx_ocn_tfield,
     &       ppx_ocn_uall,ppx_ocn_ucomp,ppx_ocn_ufield,
     &       ppx_ocn_tzonal,ppx_ocn_uzonal,ppx_ocn_tmerid,
     &       ppx_ocn_umerid,ppx_ocn_scalar,ppx_ocn_rim,
     &       ppx_ocn_cuall,ppx_ocn_cvall,    
     &       ppx_wam_all,ppx_wam_sea,ppx_wam_rim
C Valid rotation type codes
      INTEGER
     &       ppx_unrotated,ppx_elf_rotated
C Valid level type codes
      INTEGER
     &       ppx_full_level,ppx_half_level
C Valid data type codes
      INTEGER
     &       ppx_type_real,ppx_type_int,ppx_type_log
C Valid meto8 level type codes
      INTEGER
     &       ppx_meto8_surf
C Valid dump packing codes
      INTEGER
     &       ppx_pack_off,ppx_pack_32,ppx_pack_wgdos,ppx_pack_cfi1
C
C
C
C
      PARAMETER(
     &       ppx_atm_nonstd=0,      ! Non-standard atmos grid
     &       ppx_atm_tall=1,        ! All T points (atmos)
     &       ppx_atm_tland=2,       ! Land-only T points (atmos)
     &       ppx_atm_tsea=3,        ! Sea-only T points (atmos)
     &       ppx_atm_tzonal=4,      ! Zonal field at T points (atmos)
     &       ppx_atm_tmerid=5,      ! Merid field at T points (atmos)
     &       ppx_atm_uall=11,       ! All u points (atmos)
     &       ppx_atm_uland=12,      ! Land-only u points (atmos)
     &       ppx_atm_usea=13,       ! Sea-only u points (atmos)
     &       ppx_atm_uzonal=14,     ! Zonal field at u points (atmos)
     &       ppx_atm_umerid=15,     ! Merid field at u points (atmos)
     &       ppx_atm_scalar=17,     ! Scalar (atmos)
     &       ppx_atm_cuall=18,      ! All C-grid (u) points (atmos)
     &       ppx_atm_cvall=19,      ! All C-grid (v) points (atmos)
     &       ppx_atm_compressed=21, ! Compressed land points (atmos)
     &       ppx_atm_ozone=22,      ! Field on ozone grid (atmos)
     &       ppx_atm_rim=25,        ! Rim type field (LAM BCs atmos)
     &       ppx_ocn_nonstd=30,     ! Non-standard ocean grid
     &       ppx_ocn_tcomp=31,      ! Compressed T points (ocean)
     &       ppx_ocn_ucomp=32,      ! Compressed u points (ocean)
     &       ppx_ocn_tall=36,       ! All T points incl. cyclic (ocean)
     &       ppx_ocn_uall=37,       ! All u points incl. cyclic (ocean)
     &       ppx_ocn_cuall=38,      ! All C-grid (u) points (ocean)
     &       ppx_ocn_cvall=39,      ! All C-grid (v) points (ocean)
     &       ppx_ocn_tfield=41,     ! All non-cyclic T points (ocean)
     &       ppx_ocn_ufield=42,     ! All non-cyclic u points (ocean)
     &       ppx_ocn_tzonal=43,     ! Zonal n-c field at T points(ocean)
     &       ppx_ocn_uzonal=44,     ! Zonal n-c field at u points(ocean)
     &       ppx_ocn_tmerid=45,     ! Merid n-c field at T points(ocean)
     &       ppx_ocn_umerid=46,     ! Merid n-c field at u points(ocean)
     &       ppx_ocn_scalar=47,     ! Scalar (ocean)
     &       ppx_ocn_rim=51,        ! Rim type field (LAM BCs ocean)    
     &       ppx_wam_all=60,        ! All points (wave model)
     &       ppx_wam_sea=62,        ! Sea points only (wave model)
     &       ppx_wam_rim=65)        ! Rim type field (LAM BCs wave)
C
      PARAMETER(
     &       ppx_unrotated=0,       ! Unrotated output field
     &       ppx_elf_rotated=1)     ! Rotated ELF field
C
      PARAMETER(
     &       ppx_full_level=1,      ! Model full level
     &       ppx_half_level=2)      ! Model half level
C
      PARAMETER(
     &       ppx_type_real=1,       ! Real data type
     &       ppx_type_int=2,        ! Integer data type
     &       ppx_type_log=3)        ! Logical data type
C
      PARAMETER(
     &       ppx_meto8_surf=9999)   ! MetO8 surface type code
C
      PARAMETER(
     &       ppx_pack_off=0,        ! Field not packed (ie. 64 bit)
     &       ppx_pack_32=-1,        ! Field packed to 32 bit in dump
     &       ppx_pack_wgdos=1,      ! Field packed by WGDOS method
     &       ppx_pack_cfi1=11)      ! Field packed using CFI1 (ocean)
C
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

      IF (((grid_type_code .GE. ppx_atm_tall)    .AND.
     &     (grid_type_code .LE. ppx_atm_tmerid)) .OR.
     &     (grid_type_code .EQ. ppx_atm_cuall)   .OR.
     &     (grid_type_code .EQ. ppx_atm_ozone)   .OR.
     &     (grid_type_code .EQ. ppx_atm_compressed) .OR.
     &     (grid_type_code .EQ. ppx_ocn_tall)    .OR.
     &     (grid_type_code .EQ. ppx_ocn_cuall)  .OR.
     &     (grid_type_code .EQ. ppx_ocn_tfield)  .OR.
     &     (grid_type_code .EQ. ppx_ocn_tzonal)  .OR.
     &     (grid_type_code .EQ. ppx_ocn_tmerid)) THEN
        GET_FLD_TYPE=fld_type_p
      ELSEIF
     &   (((grid_type_code .GE. ppx_atm_uall)    .AND.
     &     (grid_type_code .LE. ppx_atm_umerid)) .OR.
     &     (grid_type_code .EQ. ppx_atm_cvall)   .OR.
     &     (grid_type_code .EQ. ppx_ocn_uall)    .OR.
     &     (grid_type_code .EQ. ppx_ocn_cvall)   .OR.
     &     (grid_type_code .EQ. ppx_ocn_ufield)  .OR.
     &     (grid_type_code .EQ. ppx_ocn_uzonal)  .OR.
     &     (grid_type_code .EQ. ppx_ocn_umerid)) THEN
        GET_FLD_TYPE=fld_type_u
      ELSE
        GET_FLD_TYPE=fld_type_unknown
      ENDIF

      RETURN

      END

