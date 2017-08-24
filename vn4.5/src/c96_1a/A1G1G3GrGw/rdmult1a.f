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
!+ Parallel UM interface to BUFFIN
!
! Subroutine Interface:
      SUBROUTINE READ_MULTI(NFT,D1,ISIZE,LEN_IO,LOCAL_LEN,
     &                      IOSTAT,LOOKUP,FIXHD12,
     &                      ADDR_INFO,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  Unified Model. It is used where each process must read in a
!  local section of a global field.
!
! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.
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
!    4.2    17/9/96   Modify send/receive maps and change args to
!                     alltoall for GCOM/GCG v1.1     P.Burton
!    4.2    12/11/96  Detects non-constant PSTAR on pole rows
!                     P.Burton
!    4.2    18/10/96  New name for group of processors in scatter_field
!                     Allowed reading of LBC files
!                     P.Burton
!    4.3    18/03/97  Added D1_ADDR argument, and rewrote field
!                     recognition tests. Now handles diagnostic
!                     fields too.                      P.Burton
!                     Fix for reading in packed LBCS   R.Barnes
!    4.4    13/06/97  Use GENERAL_SCATTER_FIELD for distributing data
!                     and change decomposition depending on model
!                     type of field being read in            P.Burton
!    4.4    25/04/97  Changes to read well-formed records if the
!                     input dumpfile is in that format (almost PP file
!                     format)
!                       Author: Bob Carruthers, Cray Research
!    4.5    13/01/98  Replace SHMEM COMMON block with dynamic array
!                                                          P.Burton
!    4.5    14/05/98  Corrected error messages to say READ_MULTI
!                                                       P.Burton
!
! Subroutine Arguments:

      INTEGER
     &  NFT          !   IN : FORTRAN unit number
     & ,ISIZE        !   IN : no. of words to be read in (global field)
     & ,LEN_IO       !  OUT : no. of words read in (global field)
     & ,LOCAL_LEN    !  OUT : no. of words in local field
     & ,LOOKUP(64)   !   IN : LOOKUP header from dump
     & ,FIXHD12      !   IN : 12th element of fixed length header
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

      CHARACTER*(80)
     &  CMESSAGE     !  OUT : Error message

      REAL
     &  IOSTAT       !  OUT : Return code
     & ,D1(*)        !  OUT : Array to read data in to

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


! Local variables

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
      INTEGER
     &  ICODE  ! return code from GENERAL_SCATTER_FIELD
     &, info   ! GCOM return code
     &, I      ! loop counter

      INTEGER
     &  orig_decomp  ! decomposition on entry
     &, new_decomp   ! decomposition to change to
     &, io_ret_codes(2) ! return code from I/O
      REAL pstar_pole_val

      INTEGER pstar_const
      COMMON /RDMULT_SHMEM_COM/ pstar_const

      REAL buf(ISIZE*2)  ! buffer for reading the field into
CDIR$ CACHE_ALIGN buf

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
            WRITE(6,*) 'ERROR : READ_MULTI'
            WRITE(6,*) 'Failed to change decomposition to ',new_decomp
            WRITE(6,*) 'Field M,S,I ',
     &                  ADDR_INFO(d1_imodl),
     &                  ADDR_INFO(d1_section),
     &                  ADDR_INFO(d1_item)
          ENDIF
          IOSTAT=-100
          CMESSAGE='READ_MULTI : Failed to change decomposition'
          GOTO 9999
        ENDIF

      ENDIF
      
! First thing to do is to read the field in to PE 0

      IF (mype .EQ. 0) THEN
      
! If lateral boundary data are 32-bit packed,
!  need to halve buffin length,
        IF (LOOKUP(LBHEM) .EQ. 99 .AND.
     &      MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
          len_io=(isize+1)/2
          CALL BUFFIN_SINGLE(NFT,buf,(ISIZE+1)/2,LEN_IO,IOSTAT)
!   and double LEN_IO to satisfy tests against ISIZE.
          LEN_IO = 2*LEN_IO
        ELSE
          CALL BUFFIN_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
        END IF
        
!       Has the data been read in OK?
        IF (.NOT.((IOSTAT .NE. -1.0) .OR. (LEN_IO .NE. ISIZE))) THEN
        
! We must check to see if it is a 32 bit field on disk, and if
! so, expand it before doing anything with it.
          IF (MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
            IF (LOOKUP(DATA_TYPE) .EQ. 1) THEN
! For special case of lateral boundary data, the length
! is given by ISIZE.
              IF (LOOKUP(LBHEM) .EQ. 99) THEN
                CALL EXPAND32B( ISIZE , buf, FIXHD12 )
              ELSE
              CALL EXPAND32B( LOOKUP(LBLREC) , buf, FIXHD12 )
              ENDIF
            ELSE
              IOSTAT=100
            ENDIF
          ENDIF
        ENDIF
        
        io_ret_codes(1)=LEN_IO
        io_ret_codes(2)=IOSTAT
      ENDIF  ! IF (mype .EQ. 0)

! Broadcast the error codes, so we can check if anything's wrong      
      CALL GC_IBCAST(333,2,0,nproc,info,io_ret_codes)
      LEN_IO=io_ret_codes(1)
      IOSTAT=io_ret_codes(2)
      
      IF ((IOSTAT .NE. -1.0) .OR. (LEN_IO .NE. ISIZE)) GOTO 9999
              
! Check the field for non-constant PSTAR at poles
      IF ((ADDR_INFO(d1_imodl) .EQ. ATMOS_IM) .AND.
     &    (ADDR_INFO(d1_object_type) .EQ. prognostic) .AND.
     &    (ADDR_INFO(d1_section) .EQ. 0) .AND.
     &    (ADDR_INFO(d1_item) .EQ. 1)) THEN  ! this is PSTAR

        pstar_const=0
        IF (mype .EQ. 0) THEN

          pstar_pole_val=buf(1)
          DO I=2,glsize(1)  ! NP row
            IF ( buf(I) .NE. pstar_pole_val ) pstar_const=1
          ENDDO

          pstar_pole_val=buf((glsize(2)-1)*glsize(1) + 1)
          DO I=2,glsize(1)  ! SP row
            IF ( buf(I+ (glsize(2)-1)*glsize(1)) .NE. pstar_pole_val )
     &        pstar_const=1
          ENDDO

        ENDIF ! if mype .eq. 0

        CALL GC_IMAX(1,nproc,info,pstar_const)

        IF (pstar_const .EQ. 1) THEN
          WRITE(6,*) 'Non constant polar PSTAR found in dump'
          WRITE(6,*) 'Model run aborted'
          IOSTAT=1.0
          CMESSAGE='Non constant polar PSTAR found in dump'
          GOTO 9999
        ENDIF

      ENDIF  ! is this a pstar field

! Now we can distribute it out to the other processes

! For atmosphere zonal ozone fields - set to zonal grid type

      IF (( ADDR_INFO(d1_grid_type) .EQ. ppx_atm_ozone) .AND.
     &    ( LOOKUP(LBNPT) .EQ. 1)) THEN

        ADDR_INFO(d1_grid_type)=ppx_atm_tzonal

      ENDIF

! Now decompose the field in buf to the local D1 arrays

      ICODE=0

      CALL GENERAL_SCATTER_FIELD(
     &  D1,buf,LOCAL_LEN,LOOKUP(LBLREC),
     &  ADDR_INFO,0,
     &  ICODE,CMESSAGE)

      IF (ICODE .EQ. 1) THEN
        WRITE(6,*) 'READ_MULTI: Field number ',LOOKUP(ITEM_CODE),
     &             'with dimensions ', LOOKUP(LBNPT),' x ',
     &             LOOKUP(LBROW),' and gridtype ',
     &             ADDR_INFO(d1_grid_type),
     &             'was unrecognized and not read in.'
        IOSTAT=1.0
        CMESSAGE='Unrecognized field on read'
      ELSEIF (ICODE .NE. 0) THEN
        IOSTAT=2.0
      ENDIF

      IF (new_decomp .NE. orig_decomp) THEN  ! change back

        icode=0
        CALL CHANGE_DECOMPOSITION(orig_decomp,icode)

      ENDIF
      
9999  CONTINUE
      
      RETURN
      END
