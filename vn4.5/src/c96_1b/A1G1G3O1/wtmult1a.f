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
!+ Parallel UM interface to BUFFOUT
!
! Subroutine Interface:
      SUBROUTINE WRITE_MULTI(NFT,D1,ISIZE,LEN_IO,LOCAL_LEN,
     &                       IOSTAT,LOOKUP,FIXHD12,COMPBUF,
     &                       ADDR_INFO,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Unified Model. It is used where each process must write out a
!  local section of a global field.
!
! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + D.Salmond
!    4.1    18/3/96   Simplified communications    P.Burton
!    4.2    18/11/96  Added *CALL AMAXSIZE for IOVARS comdeck
!                     Added atmos_ prefix to landmask fields P.Burton
!    4.2    18/9/96   Modify send/receive maps and change args to
!                     alltoall for GCOM/GCG v1.1     P.Burton
!    4.2    18/10/96  New name for group of processors in scatter_field
!                     P.Burton
!    4.3    18/03/97  Added D1_ADDR argument, and rewrote field
!                     recognition tests. Now handles diagnostic
!                     fields too.                         P.Burton
!    4.4    13/06/97  Use GENERAL_GATHER_FIELD for collecting data
!                     and change decomposition depending on model
!                     type of field being read in            P.Burton
!    4.5    13/01/98  Replace SHMEM COMMON block with dynamic array
!                                                          P.Burton
!
! Subroutine Arguments:

      INTEGER
     &  NFT          !   IN : FORTRAN unit number
     & ,ISIZE        !   IN : no. of words to write out
     & ,LEN_IO       !  OUT : no. of words written out
     & ,LOCAL_LEN    !  OUT : size of the local field written out
     & ,LOOKUP(64)   !   IN : LOOKUP header from dump
     & ,FIXHD12      !   IN : 12th. element of fixed length header
     &               !        required for packing fields

! Required for dimensioning ADDR_INFO
C     COMDECK D1_ADDR
C     Information for accessing D1 addressing array
C     Number of items of info needed for each object and maximum
C     number of objects in D1 -
      INTEGER
     &  D1_LIST_LEN

C Number of items of information in D1 addressing array
      PARAMETER(
     &  D1_LIST_LEN=16
     &  )
C     Names of items in D1 addressing array
      INTEGER
     &  d1_object_type   ! Prognostic, Diagnostic, Secondary or other
     &  ,d1_imodl        ! Internal model id
     &  ,d1_section      ! Section
     &  ,d1_item         ! Item
     &  ,d1_address      ! Address in D1
     &  ,d1_length       ! Record length
     &  ,d1_grid_type    ! Grid type
     &  ,d1_no_levels    ! Number of levels
     &  ,d1_stlist_no    ! Stash list number for diags. -1 for progs
     &  ,d1_lookup_ptr   ! Pointer to dump header lookup table
     &  ,d1_north_code   ! Northern row address
     &  ,d1_south_code   ! Southern row address
     &  ,d1_east_code    ! Eastern row address
     &  ,d1_west_code    ! Western row address
     &  ,d1_gridpoint_code ! gridpoint info address
     &  ,d1_proc_no_code ! Processing Code address

C Codes for items in D1 array. Update D1_LIST_LEN above if items added
      PARAMETER(
     &  d1_object_type=1
     &  ,d1_imodl=2
     &  ,d1_section=3
     &  ,d1_item=4
     &  ,d1_address=5
     &  ,d1_length=6
     &  ,d1_grid_type=7
     &  ,d1_no_levels=8
     &  ,d1_stlist_no=9
     &  ,d1_lookup_ptr=10
     &  ,d1_north_code=11
     &  ,d1_south_code=12
     &  ,d1_east_code=13
     &  ,d1_west_code=14
     &  ,d1_gridpoint_code=15
     &  ,d1_proc_no_code=16
     &  )

C     Types of items for d1_type
      INTEGER
     &  prognostic
     &  ,diagnostic
     &  ,secondary
     &  ,other

      PARAMETER(
     &  prognostic=0
     &  ,diagnostic=1
     &  ,secondary=2
     &  ,other=3
     &  )


      INTEGER
     &  ADDR_INFO(D1_LIST_LEN)   ! IN addressing info about field
      REAL
     &  IOSTAT       !  OUT : Return code
     & ,D1(*)        !   IN : Array to write out
     & ,COMPBUF(*)   !   IN : Workspace for compressing field

      CHARACTER*(80)
     &  CMESSAGE     !  OUT : Error message

! Parameters and Common blocks

C*L------------------ COMDECK LOOKADD ----------------------------------
CLL
CLL Purpose : Contains information about the format
CLL           of the PP header
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
CLL                   BRSVD2 to BHULEV and definitions for BRLEV and
CLL                   BHRLEV. Corresponding changes made to STWORK1A
CLL                   and PPHEAD1A. (Andrew Brady)   
CLL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
CLL
CLL Programming standard :
CLL
CLL Logical components covered : F092
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
C
      INTEGER
C Validity time
     &       LBYR,       ! Year
     &       LBMON,      ! Month
     &       LBDAT,      ! Day of month
     &       LBHR,       ! Hour
     &       LBMIN,      ! Minute
     &       LBDAY       ! Day number

C Data time

      INTEGER
     &       LBYRD,      ! Year
     &       LBMOND,     ! Month
     &       LBDATD,     ! Day of month
     &       LBHRD,      ! Hour
     &       LBMIND,     ! Minute
     &       LBDAYD      ! Day number

      INTEGER
     &       LBTIM,      ! Time indicator
     &       LBFT,       ! Forcast period (hours)
     &       LBLREC,     ! Length of data record
     &       LBCODE,     ! Grid type code
     &       LBHEM,      ! Hemisphere indicator
     &       LBROW,      ! Number of rows in grid
     &       LBNPT,      ! Number of points per row
     &       LBEXT,      ! Length of extra data
     &       LBPACK,     ! Packing method indicator
     &       LBREL       ! Header release number

      INTEGER
     &       LBFC,       ! Field code
     &       LBCFC,      ! Second field code
     &       LBPROC,     ! Processing code
     &       LBVC,       ! Vertical coordinate type
     &       LBRVC,      ! Coordinate type for reference level
     &       LBEXP,      ! Experiment number
     &       LBEGIN,     ! Start record
     &       LBNREC,     ! No of records-Direct access only
     &       LBPROJ,     ! Met-O-8 projection number
     &       LBTYP,      ! Met-O-8 field type
     &       LBLEV,      ! Met-O-8 level code
     &       LBRSVD1,    ! Reserved for future PP-package use
     &       LBRSVD2,    ! Reserved for future PP-package use
     &       LBRSVD3,    ! Reserved for future PP-package use
     &       LBRSVD4,    ! Reserved for future PP-package use
     &       LBSRCE      ! =1111 to indicate following apply to UM
      INTEGER
     &       DATA_TYPE,  ! Indicator for real/int or timeseries
     &       NADDR,      ! Start address in DATA_REAL or DATA_INT
     &       LBUSER3,    ! Free for user-defined function   
     &       ITEM_CODE,  ! Stash code
     &       LBPLEV,     ! Pseudo-level indicator (if defined)
     &       LBUSER6,    ! Free for user-defined function
     &       MODEL_CODE ! internal model identifier
      INTEGER
     &       BULEV,      ! Upper level boundary (Bk for ATMOS)
     &       BHULEV,     ! Upper level boundary (Ak for ATMOS)   
     &       BRSVD3,     ! Reserved for future PP-package use
     &       BRSVD4,     ! Reserved for future PP-package use
     &       BDATUM,     ! Datum value
     &       BACC,       ! (Packed fields) Packing accuracy
     &       BLEV,       ! Level
     &       BRLEV,      ! Lower level boundary (Bk for ATMOS)   
     &       BHLEV,      ! (Hybrid levels) A-level of value
     &       BHRLEV,     ! Lower level boundary (Ak for ATMOS)   
     &       BPLAT,      ! Real latitude of 'pseudo' N Pole
     &       BPLON,      ! Real longitude of 'pseudo' N Pole
     &       BGOR,       ! Grid orientation
     &       BZY,        ! Zeroth latitude
     &       BDY,        ! Latitude interval
     &       BZX,        ! Zeroth longitude
     &       BDX,        ! Longitude interval
     &       BMDI,       ! Missing data indicator
     &       BMKS        ! M,K,S scaling factor

C Mapping of MPP_LOOKUP; analogous to mapping in PP header

      INTEGER
     &       P_NADDR,    ! Address on local PE
     &       P_LBLREC    ! Local length of record

      PARAMETER (
     &       P_NADDR=1,
     &       P_LBLREC=2)
C*----------------------------------------------------------------------
C NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
C ITEM_CODE is the location in PP header for a code defined as
C           (section number)*1000+item number
C DATA_TYPE is the location in the PP header defining data as REAL or
C           INTEGER.
C LBNPT is the location defining the number of points per row
C
      PARAMETER(
C Validity time
     &       LBYR=1,
     &       LBMON=2,
     &       LBDAT=3,
     &       LBHR=4,
     &       LBMIN=5,
     &       LBDAY=6,
C Data time
     &       LBYRD=7,
     &       LBMOND=8,
     &       LBDATD=9,
     &       LBHRD=10,
     &       LBMIND=11,
     &       LBDAYD=12)

      PARAMETER (
     &       LBTIM=13,
     &       LBFT=14,
     &       LBLREC=15,
     &       LBCODE=16,
     &       LBHEM=17,
     &       LBROW=18,
     &       LBNPT=19,
     &       LBEXT=20,
     &       LBPACK=21,
     &       LBREL=22,
     &       LBFC=23,
     &       LBCFC=24,
     &       LBPROC=25,
     &       LBVC=26,
     &       LBRVC=27)

      PARAMETER (
     &       LBEXP=28,
     &       LBEGIN=29,
     &       LBNREC=30,
     &       LBPROJ=31,
     &       LBTYP=32,
     &       LBLEV=33,
     &       LBRSVD1=34,
     &       LBRSVD2=35,
     &       LBRSVD3=36,
     &       LBRSVD4=37,
     &       LBSRCE=38,
     &       DATA_TYPE=39,
     &       NADDR=40,
     &       LBUSER3=41,    
     &       ITEM_CODE=42,
     &       LBPLEV=43,
     &       LBUSER6=44,
     &       MODEL_CODE=45)

      PARAMETER (
     &       BULEV=46,
     &       BHULEV=47, 
     &       BRSVD3=48,
     &       BRSVD4=49,
     &       BDATUM=50,
     &       BACC=51,
     &       BLEV=52,
     &       BRLEV=53,
     &       BHLEV=54,
     &       BHRLEV=55,
     &       BPLAT=56,
     &       BPLON=57,
     &       BGOR=58,
     &       BZY=59,
     &       BDY=60,
     &       BZX=61,
     &       BDX=62,
     &       BMDI=63,
     &       BMKS=64)

C
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
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      INTEGER
     *   A_IM,ATMOS_IM        ! Atmosphere internal model
     *  ,O_IM,OCEAN_IM        ! Ocean      internal model
     *  ,S_IM, SLAB_IM        ! Slab       internal model
     *  ,W_IM, WAVE_IM        ! Wave       internal model
     *  ,I_IM,SEAICE_IM       ! Sea-ice    internal model
     *  ,N_IM,NATMOS_IM       ! New dynamics (Charney-Phillips grid)
!                               atmosphere internal model
!
      PARAMETER(
     *   A_IM=1,ATMOS_IM=1       ! Atmosphere internal model
     *  ,O_IM=2,OCEAN_IM=2       ! Ocean      internal model
     *  ,S_IM=3, SLAB_IM=3       ! Slab       internal model
     *  ,W_IM=4, WAVE_IM=4       ! Wave       internal model
     *  ,I_IM=5,SEAICE_IM=5      ! Sea-ice    internal model
     *  ,N_IM=6,NATMOS_IM=6      ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
      INTEGER
     *   A_SM,ATMOS_SM        ! Atmosphere submodel partition
     *  ,O_SM,OCEAN_SM        ! Ocean      submodel partition
     *  ,W_SM, WAVE_SM        ! Wave       submodel partition
     *  ,N_SM,NATMOS_SM       ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
!
      PARAMETER(
     *   A_SM=1,ATMOS_SM=1    ! Atmosphere submodel partition
     *  ,O_SM=2,OCEAN_SM=2    ! Ocean      submodel partition
     *  ,W_SM=4, WAVE_SM=4    ! Wave       submodel partition
     *  ,N_SM=6,NATMOS_SM=6   ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
C
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

      REAL buf(ISIZE*2)  ! Buffer for holding data to be written.
!                        ! Factor of two is incase the data is
!                        ! packed (ISIZE is the size on disk)
CDIR$ CACHE_ALIGN buf
      INTEGER
     &  ICODE ! return code from GENERAL_GATHER_FIELD
     &, I      ! loop counter

      INTEGER
     &  orig_decomp  ! decomposition on entry
     &, new_decomp   ! decomposition to change to

      EXTERNAL CHANGE_DECOMPOSITION,GENERAL_GATHER_FIELD,
     &         BUFFOUT_SINGLE,
     &         PACK21,EXPAND21,P21BITS
      INTEGER P21BITS

! ------------------------------------------------------------------
      IOSTAT=-1.0
      ICODE=-1
      LEN_IO=ISIZE
      LOCAL_LEN=0

      orig_decomp=current_decomp_type
      new_decomp=orig_decomp

      IF ((ADDR_INFO(d1_imodl) .EQ. ATMOS_IM) .AND.
     &    (orig_decomp .NE. decomp_standard_atmos)) THEN

        new_decomp=decomp_standard_atmos

      ELSEIF ((ADDR_INFO(d1_imodl) .EQ. OCEAN_IM) .AND.
     &        (ADDR_INFO(d1_object_type) .EQ. prognostic) .AND.
     &        (orig_decomp .NE. decomp_standard_ocean)) THEN

        new_decomp=decomp_standard_ocean

      ELSEIF ((ADDR_INFO(d1_imodl) .EQ. OCEAN_IM) .AND.
     &        (ADDR_INFO(d1_object_type) .NE. prognostic) .AND.
     &        (orig_decomp .NE. decomp_nowrap_ocean)) THEN

        new_decomp=decomp_nowrap_ocean

      ENDIF

      IF (new_decomp .NE. orig_decomp) THEN

        icode=0
        CALL CHANGE_DECOMPOSITION(new_decomp,icode)

        IF (icode .NE. 0) THEN
          IF (mype .EQ. 0) THEN
            WRITE(6,*) 'ERROR : WRITE_MULTI'
            WRITE(6,*) 'Failed to change decomposition to ',new_decomp
            WRITE(6,*) 'Field M,S,I ',
     &                  ADDR_INFO(d1_imodl),
     &                  ADDR_INFO(d1_section),
     &                  ADDR_INFO(d1_item)
          ENDIF
          IOSTAT=-100
          CMESSAGE='WRITE_MULTI : Failed to change decomposition'
          GOTO 9999
        ENDIF

      ENDIF

! Gather the field from the local D1 array to buf

      ICODE=0

      CALL GENERAL_GATHER_FIELD(
     &  D1,buf,LOCAL_LEN,LOOKUP(LBLREC),
     &  ADDR_INFO,0,
     &  ICODE,CMESSAGE)

      IF (ICODE .EQ. 1) THEN
        WRITE(6,*) 'WRITE_MULTI: Field number ',LOOKUP(ITEM_CODE),
     &             'with dimensions ', LOOKUP(LBNPT),' x ',
     &             LOOKUP(LBROW),' and gridtype ',
     &             ADDR_INFO(d1_grid_type),
     &             'was unrecognized and not written out.'
        IOSTAT=1.0
        CMESSAGE='Unrecognized field on write'
        GOTO 9999
      ELSEIF (ICODE .NE. 0) THEN
        IOSTAT=2.0
        GOTO 9999
      ENDIF

! ------------------------------------------------------------------
! ------------------------------------------------------------------
! And finally the code to write the global field in array buf
! out to disk.

      IF (mype .EQ. 0) THEN      
!       Does this field need to be compressed?
        IF(MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
          IF(LOOKUP(DATA_TYPE) .EQ. 1) THEN
            CALL PACK21(LOOKUP(LBLREC),buf,
     &                  COMPBUF,P21BITS(FIXHD12))
          ENDIF
        ELSE ! no compression required - just do a copy
          DO i=1,LOOKUP(LBLREC) 
            COMPBUF(i)=buf(i)
          ENDDO
        ENDIF
        
! Now write out the global field
      
        CALL buffout_single(NFT,COMPBUF,ISIZE,LEN_IO,IOSTAT)
      
      ENDIF ! am I PE 0 ?

!     CALL GC_GSYNC(nproc,info)
      
! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

      IF(MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
        IF(LOOKUP(DATA_TYPE) .EQ. 1) THEN
          CALL PACK21(LOCAL_LEN,D1,COMPBUF,P21BITS(FIXHD12))
          CALL EXPAND21(LOCAL_LEN,COMPBUF,D1,P21BITS(FIXHD12))
        ENDIF
      ENDIF

      IF (new_decomp .NE. orig_decomp) THEN  ! change back

        icode=0
        CALL CHANGE_DECOMPOSITION(orig_decomp,icode)

      ENDIF
 9999 CONTINUE  ! point to jump to if there is a failure
     
      RETURN
      END
