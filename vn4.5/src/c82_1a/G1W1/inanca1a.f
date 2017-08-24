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
CLL Subroutine INANCILA
CLL
CLL Purpose : Takes as input,the code defining the frequency of update
CLL           of ancillary fields as set by the user interface.
CLL           Converts them into a list of numbers of timesteps after
CLL           which each field must be updated, and calculates the
CLL           frequency with which this list must be interrogated.
CLL           Where the update interval is in months or years,
CLL           the check will be carried out each day. The physical
CLL           files required are also determined by input code,
CLL           and the headers and lookup tables are read into
CLL           the arguments FIXHD,INTHD,LOOKUP which are in
CLL           COMMON/ANCILHDA/ of calling routine INANCCTL.
CLL           Indexes for each possible ancillary field are set up in
CLL           COMMON/IXANCILA/
CLL
CLL Level 2 Control routine for CRAY YMP
CLL
CLL CW, DR      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  22/02/93  Changes to add 2 SLAB fields (STASH items 178,179)
CLL                  - to be updated from existing atmosphere files.
CLL   3.3  22/11/93  Add aerosol ancillary fields.  R T H Barnes.
CLL   3.3  21/12/93  Fix put in to prevent array 'out of bounds'
CLL                  problem in section 1.6. Problem to be investigated
CLL                  for 3.4 D. Robinson.
CLL   3.4  16/06/94  DEF CAL360 replaced by LOGICAL LCAL360
CLL                                                   S.J.Swarbrick
CLL  3.4  05/09/94  Add murk and user ancillary fields.  RTHBarnes.
CLL   3.4  22/06/94  Array 'out of bounds' problem solved. D. Robinson
CLL  3.4   11/10/94   Part of modset which sorts out some handling
CLL                   of unset data by recon_dump.
CLL                   Necessary to port model to a T3D.
CLL                   Author D.M. Goddard
CLL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
CLL  3.5  24/07/95  Check fields for updating have valid address. RTHB
!    4.0  01/09/95  Add diagnostic information to output about
!                   ozone ancillary fields and test correct ozone
!                   data provided.  D. Goddard & D. Robinson
CLL  4.0  10/10/95  Set LOOKUP(45) in ancillary files. D. Robinson.
CLL
CLL  4.0  29/09/95  Need extra rewind of namelist file. RTHBarnes.
CLL  4.0  05/08/95  Temporary solution to get round problem of
CLL                 no. of soil moisture levels being hard-wired
CLL                 to no. of deep soil temperature levels
CLL                 This causes a problem with introduction of
CLL                 Penman-Monteith BL code at 4.0 - use if test
CLL                 on number of deep soil temperature
CLL                 levels which is set to 4 for Penman-Monteith code
CLL                 (set to 3 for all other BL options)
CLL                 Permanent solution suggested for 4.1
CLL                 search on C**** for comments
CLL                 J Smith
CLL  4.0  06/01/96  SI array received for two internal models (atmos
CLL                 and slab) in argument list. Hardwire processing of
CLL                 slab ancillary field (item code 177) to use
CLL                 SI_SLAB. D. Robinson
!    4.1  03/05/96  Use READHEAD to read in ancillary file headers.
!                   D. Robinson
!    4.1  18/06/96  Changes to cope with changes in STASH addressing
!                   Author D.M. Goddard.
CLL  4.1  22/05/96  Call new CANC* comdecks. Use new arrays in
CLL                 CANCFLDA. Cater for new sulphur ancillary files.
CLL                 Remove hardwired fix for slab ancillary fields
CLL                 introduced at 4.0 D. Robinson.
!LL  4.4  28/07/97  Add LAMIPII to namelist for special updating of 
!LL                 ice in AMIP II runs. R A Stratton
CLL  4.4  16/09/97  Set number of headers for multi-pseudo-level 
CLL                 ancillary fields for surface and vegetation types.
CLL                                              Richard Betts
!LL  4.4  09/09/97  New namelist UPANCA for updating information.
!LL                 D. Robinson.
!LL  4.4  10/09/97  Check calendar indicator in Anc File. D Robinson.
!    4.5  22/10/98  Set LEVELS array for new user hulti-layer 
!                   ancillary fields
!                   Author D.M Goddard
!LL  4.5  19/01/98  Remove SOIL_VARS and VEG_VARS. D. Robinson.
!LL  4.5  05/05/98  Improve error message for missing files. R. Rawlins
CLL
CLL System components covered : C710
CLL
CLL System task : C7
CLL
CLL Documentation : Unified Model Documentation Paper No C7
CLL                 Version No 4  dated 15/06/90
CLLEND
      SUBROUTINE INANCILA(LEN_FIXHD,LEN_INTHD,LEN_REALHD,  !Intent (In)
     &                    LEN1_LEVDEPC,LEN2_LEVDEPC,
     &                    FIXHD,INTHD,REALHD,LOOKUP,
     &                    A_FIXHD,A_REALHD,A_LEVDEPC,
     &                    NDATASETS,NLOOKUPS,FTNANCIL,
     &                    LOOKUP_START,LEN1_LOOKUP,ROW_LENGTH,
     &                    P_ROWS,U_ROWS,P_LEVELS,
     &                    TR_LEVELS,ST_LEVELS,SM_LEVELS,
     &                    OZONE_LEVELS,TITLE,
     &                    SI_ATMOS,SI_SLAB,SILEN,
     &                    ANCILLARY_STEPS,STEPS_PER_HR,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &                    ICODE,CMESSAGE,LCAL360)         ! Intent (Out)


      IMPLICIT NONE

      LOGICAL LCAL360  ! Logical switch for 360-day calendar

      INTEGER
     &        LEN_FIXHD,       ! Length of header blocks in ancillary
C                              ! data sets
     &        LEN_INTHD,       !
     &        LEN_REALHD,      !
     &        LEN1_LEVDEPC,    ! Dimension of LEVDEPC in model
     &      LEN2_LEVDEPC
     &     ,ANCILLARY_STEPS,
     &        STEPS_PER_HR


      INTEGER
     &        NDATASETS,       ! No of physical files
     &        NLOOKUPS,        ! No of lookups required(set by User I.)
     &                IOUNIT,
     &        FTNANCIL(NDATASETS), ! Fortran nos of physical files
     &        LOOKUP_START(NDATASETS),!start of each individual lookup
C                                     !in overall LOOKUP array
     &        LEN1_LOOKUP,     ! Length of PP header
     &        ROW_LENGTH,      ! Atmosphere model dimensions
     &        P_ROWS,          ! No. of rows for pressure-type variables
     &        U_ROWS,          ! No. of rows for wind-type variables
     &        P_LEVELS,        ! No. of pressure levels
     &        TR_LEVELS,       ! No. of tracer levels
     &        FILE_LEVELS,     ! Number of levels of data in files
C                              ! contining multi-level data.
     &        ST_LEVELS,       ! No. of soil temperature levels 
     &        SM_LEVELS,       ! No. of soil moisture levels
     &      OZONE_LEVELS

!      For atmos only runs SI_SLAB is a copy of SI_ATMOS
!      SI_SLAB is only used in SLAB runs.

     &       ,SILEN             ! Length for SI_ATMOS/SLAB arrays
     &       ,SI_ATMOS(SILEN)   ! ) STASHin addresses of atmos and
     &       ,SI_SLAB(SILEN)    ! ) slab ancillary fields.
      CHARACTER*80 TITLE(NDATASETS) ! Titles of each dataset

      INTEGER
     &        FIXHD(LEN_FIXHD,NDATASETS),! Overall Fixed header array
     &        A_FIXHD(LEN_FIXHD), ! Fixed header for Dump
     &        INTHD(LEN_INTHD,NDATASETS),! Overall Integer header array
     &        LOOKUP(LEN1_LOOKUP,NLOOKUPS),!Overall Lookup array
     &        ICODE            ! Return code =0 Normal Exit
C                              !             >0 Error

      REAL
     &      REALHD(LEN_REALHD,NDATASETS),!
     &      A_REALHD(LEN_REALHD),!
     &      A_LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC),
     &      LEVDEPC(P_LEVELS*4)! Space to hold level dependent constants
C                              ! from data set

      CHARACTER*100
     &        CMESSAGE         ! Out error message if I>0

! Comdecks:----------------------------------------------------------
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
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

!
!  2. Maximum internal model/submodel array sizes for this version.
!
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of 
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.         
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      INTEGER
     * N_INTERNAL_MODEL_MAX      ! Max no. of internal models
     *,N_SUBMODEL_PARTITION_MAX  ! Max no. of submodel dump partitions
     *,INTERNAL_ID_MAX           ! Max value of internal model id
     *,SUBMODEL_ID_MAX           ! Max value of submodel dump id

      PARAMETER(
     * N_INTERNAL_MODEL_MAX=4,                                          
     * N_SUBMODEL_PARTITION_MAX=4,
     * INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX,
     * SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX)
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER
     * N_INTERNAL_MODEL          ! No. of internal models
     *,N_SUBMODEL_PARTITION      ! No. of submodel partitions
     *,INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX) ! Internal models
     *,SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX) ! Submodel identifier
     *                           ! for each internal model in list
     &,SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX) ! Submodel number for
!                                  each submodel id
!
! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION
     *,INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM
!
!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
      INTEGER
     * N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)  ! No of internal models in
!              each submodel partition indexed by sm identifier
     *,SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)    ! List of
!              submodel partition identifiers
     *,SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)  ! Submodel partition
!              identifier indexed by internal model identifier
     *,INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)      ! Sequence number of
!              internal model indexed by internal model identifier:
!              required to map from id to STASH internal model sequence
      LOGICAL
     * LAST_IM_IN_SM(INTERNAL_ID_MAX)      ! Last internal model within
!                                a submodel partition if .TRUE.,
!                                indexed by internal model id.
! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,
     *     INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,
     *     N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,
     *     SUBMODEL_PARTITION_INDEX,
     *     INTERNAL_MODEL_INDEX,
     *     LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      LOGICAL
     * new_im   ! new internal model next group of timesteps if .true.
     *,new_sm   ! new submodel dump  next group of timesteps if .true.

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT                                    
C
C
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
! COMDECK PPXLOOK
! Description:
!
!   Declares ppxref look-up arrays used by the UM and associated
!    arrays and parameters.
!   Comdecks CSUBMODL,CPPXREF must be *CALLed before this         
!    comdeck
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       May. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.0       Dec. 95   Replace dynamic dim ppxRecs with 
!                     NUM_DIAG_MAX in PPXC   N. Farnon
! 4.1       July 96   *CALL VERSION introduced - NUM_DIAG_MAX made
!                      equal to NDIAGP.
!                     NUM_USR_DIAG_MAX increased from 200 to 300
!                      (just in case).    
! 4.4       03/11/97  Removed MKPPXRF *DEF references. K Rogers
! 4.4       04/11/97  Changed -RECON def line to allow for other small
!                     execs which had used the RECON def. K Rogers
!
! Declarations:

! Global parameters:
! COMDECK VERSION
! Description:                                                          
!   STASH parameter definitions
!                                                                       
! Current code owner: S.J.Swarbrick                                     
!                                                                       
! History:                                                              
! Version   Date      Comment                                           
! -------   ----      -------                                           
! 3.5       Mar. 95   Original code.  S.J.Swarbrick                     
! 4.0                                 S.J.Swarbrick
! 4.1       Apr. 96   Rationalise MDI  S.J.Swarbrick
!  4.1  29/02/96  Increase OUTFILE_E.  RTHBarnes.
!  4.2  27/11/96  MPP code : Increase NELEMP   P.Burton
!  4.3  04/06/97  Increase NELEMP for D1 addressing S.D.Mullerworth
!  4.4  04/06/97  Increase NELEMP for sampling offset. S.D.Mullerworth 
!  4.5  28/01/98  Increade NELEMP for MPP code.   P.Burton
!  4.5  18/09/98  Modify name of VERSION common block to stop potential
!                 clashes with Fortran variable names          P.Burton
!  4.5  30/09/98  Increase NRECDP from 600 to 800. D. Robinson.
!                                                                       
! Declarations:                                                         
      INTEGER   NSECTP           !Max. no. of STASH sections            
      PARAMETER(NSECTP=99)       ! per internal model (44 in practise)  
      INTEGER   NITEMP           !Max. no. of STASH items per           
      PARAMETER(NITEMP=512)      ! section                              
      INTEGER   NRECDP           !Max. no. of STASH list records        
      PARAMETER(NRECDP=800)      ! (prognostic + diagnostic)
      INTEGER   NTIMEP           !Max. no. of output times tables       
      PARAMETER(NTIMEP=100)      ! in STASHC                            
      INTEGER   NPROFTP          !Max. no. of time profiles 
      PARAMETER(NPROFTP=100)     ! in STASHC                            
      INTEGER   NPROFDP          !Max. no. of domain profiles/levels    
      PARAMETER(NPROFDP=100)     ! lists in STASHC (used for both)      
      INTEGER   NTimSerP         !Max. total no. of time series 
      PARAMETER(NTimSerP=1500)    ! in STASHC
      INTEGER   tsdp             !Max. no. time series per
      PARAMETER(tsdp=250)        ! domain profile 
      INTEGER   NPROFUP          !Max. no. of useage profiles           
      PARAMETER(NPROFUP=40)      ! in STASHC                            
      INTEGER   NLEVP            !Max. no. of levels in a               
      PARAMETER(NLEVP=50)        ! levels list                          
      INTEGER   NPSLEVP          !Max. no. of pseudo levels in a        
      PARAMETER(NPSLEVP=40)      ! pseudo levels list                   
      INTEGER   NPSLISTP         !Max. no. of pseudo levels lists       
      PARAMETER(NPSLISTP=40)     ! in STASHC                            
      INTEGER   NDIAGP           !Max. no. non-blank records in         
      PARAMETER(NDIAGP=1800)     ! PPXREF file  
      INTEGER   NDIAGPM          !Same as NRECDP                        
      PARAMETER(NDIAGPM=NRECDP)  ! (will be tidied)                     
      INTEGER   NELEMP           !No. of elements in a ppxref record    
      PARAMETER(NELEMP=33)
      INTEGER   NLEVP_S
      PARAMETER(NLEVP_S=NLEVP*6+1)
      INTEGER   NLEVLSTSP
      PARAMETER(NLEVLSTSP=NPROFDP)                                      
      INTEGER   NMEANP           !No. of meaning periods 
      PARAMETER(NMEANP=4)
! OUTFILE_S, OUTFILE_L and OUTFILE_E must be consistent with
! NUNITS and NUNITS_LEN in comdeck CHSUNITS.
      INTEGER   OUTFILE_S        !Range of                              
      PARAMETER(OUTFILE_S=20)    ! output file                          
      INTEGER   OUTFILE_E        ! numbers                              
      PARAMETER(OUTFILE_E=149)   ! 
      INTEGER   OUTFILE_L
      PARAMETER(OUTFILE_L=OUTFILE_E-OUTFILE_S+1)
!Global scalar:    
      CHARACTER*55 STASH_SET     !Names of stasets files                
!Common block:                                                          
      COMMON/common_VERSION/ STASH_SET
C-----------------------------------------------------------------

! No. of STASH items per section
      INTEGER      PPXREF_ITEMS
        PARAMETER (PPXREF_ITEMS    =NITEMP) 
! No. of STASH sections per internal model
      INTEGER      PPXREF_SECTIONS
        PARAMETER (PPXREF_SECTIONS =NSECTP-55)    
! Max. number of non-null records in ppxref file (>1200) 
      INTEGER      NUM_DIAG_MAX
        PARAMETER (NUM_DIAG_MAX    =NDIAGP)
! Max. number of user-defined ppxref records allowed
      INTEGER      NUM_USR_DIAG_MAX
        PARAMETER (NUM_USR_DIAG_MAX=300)

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
      INTEGER      ppxRecs

! Global arrays:
! ppxref look-up arrays
      INTEGER   PPXI(ppxRecs,PPXREF_CODELEN)
      CHARACTER PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN)
! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
      INTEGER   PPXI_U(NUM_USR_DIAG_MAX,PPXREF_CODELEN)
      CHARACTER PPXC_U(NUM_USR_DIAG_MAX,PPXREF_CHARLEN)
! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
      CHARACTER OriginFlag(NUM_DIAG_MAX)
! Array of indices to identify which ppxref record corresponds to
!   any given row of PPXI, PPXC
      INTEGER   RowIndex(NUM_DIAG_MAX)
! Pointer array for PPXI, PPXC arrays  
      INTEGER PPXPTR
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS)

! Common block:
      COMMON/PPX_INT/ RowIndex,PPXI_U
      COMMON/PPX_CHA/ OriginFlag,PPXC_U
! - End --------------------------------------------------------------
! COMDECK MODEL     
! Description:
!   Defines model-dependent quantities used by data addressing
!     and STASH.
!
!   Comdecks CSUBMODL and VERSION must be *CALLed before this comdeck
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.1       Apr. 96   Generalisation and incorporation of
!                      wave model     S.J.Swarbrick
! 4.2       28/11/96  MPP code : Added variables
!                                global_LPRIM and global_LDUMP to
!                                provide relevant information for
!                                 global (dump)data.  P.Burton
! 4.5       29/07/98  Remove redundant code. Processing for Boundary
!                     files moved to INTF_CTL. D. Robinson.
!
! Declarations:
! Imported global variables:
!

! Global parameters:

      INTEGER      AASSETS
      PARAMETER   (AASSETS=9)
      INTEGER      MEAD_TYPES
      PARAMETER   (MEAD_TYPES=4)

      REAL         H_A_EWSPACE ,H_A_NSSPACE      
      REAL         H_A_FIRSTLAT,H_A_FIRSTLONG 
      REAL         H_A_POLELAT ,H_A_POLELONG     

      INTEGER      H_A_GROUP
      INTEGER      H_OROG_ROUGH
      INTEGER      A_ASSMGRPS
      INTEGER      NUM_PVPR

      LOGICAL      A_RECON
      LOGICAL      H_OROG_GRAD
      LOGICAL      ATMODS
      LOGICAL      CMODS
      LOGICAL      LMESO

      LOGICAL      TRACER_A    (0:29)
      LOGICAL      AASSET   (AASSETS)
      INTEGER      MAX_AOBS
      PARAMETER   (MAX_AOBS=100)
      INTEGER      AOBINC   (MAX_AOBS)
      INTEGER      AOBGRP   (MAX_AOBS)
      INTEGER      AASPF    (AASSETS)
      INTEGER      AASPL    (AASSETS)
      INTEGER      RUN_TARGET_END( 6)

      COMMON/MODELA/
     & H_A_EWSPACE    ,H_A_NSSPACE       ,H_A_FIRSTLAT
     &,H_A_FIRSTLONG  ,H_A_POLELAT       ,H_A_POLELONG
     &,A_ASSMGRPS     ,NUM_PVPR
     &,A_RECON        ,H_A_GROUP     
     &,H_OROG_GRAD    ,ATMODS            ,CMODS
     &,LMESO
     &,TRACER_A       ,AASSET            ,AASPF
     &,AASPL             

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)     
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)     
!Total data length for diagnostic flds for each submodel data partition 
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX) 
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model          
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX) 
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX) 
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)    
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX) 
!Total length of extra space for each submod. data part. 
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)     
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/
     & LPRIM       ,LDUMP    ,LSECD  ,LWORK    ,NHEAD   ,      
     & LEXTRA      ,LPRIM_O2 ,              
     & LPrimIM     ,LDumpIM  ,LSecdIM,NHeadSub ,
     & ITEM_MAX_REQ,
     & ITEM_MAX_ALL,
     & NSERBLK_S   ,NSERREC_S,NLEVL_S,NMAXLEV_S,NPSLISTS_S,    
     & NMAXPSL_S   ,LSTUSER  ,NRECS_S,NTIMES_S ,NHEAD_FILE       
     &, global_LPRIM , global_LPrimIM
     &, global_LDUMP , global_LDumpIM
      CHARACTER*1  H_ATMOS
      CHARACTER*1  H_OCEAN
      CHARACTER*1  H_SLAB
      CHARACTER*1  H_WAVE 
      CHARACTER*1  H_FLOOR
      CHARACTER*1  H_STRAT
      CHARACTER*1  H_SLAB_CAL
      CHARACTER*1  H_TOTEM
      CHARACTER*1  H_GLOBAL(N_INTERNAL_MODEL_MAX         )   
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/
     & H_ATMOS ,H_OCEAN   ,H_SLAB ,H_WAVE ,                  
     & H_GLOBAL,H_SLAB_CAL,H_TOTEM,H_FLOOR,H_STRAT   

      COMMON/HVERS/ H_VERS

      REAL H_O_EWSPACE ,H_O_NSSPACE
      REAL H_O_FIRSTLAT,H_O_FIRSTLONG
      REAL H_O_POLELAT ,H_O_POLELONG

      INTEGER H_O_PTSPROW
      INTEGER N_COMP_O
      INTEGER H_NSIDEIMTO       ,H_NSIDEJMTO
      INTEGER SEAICE_TYPE       ,OCEAN_BASINS

      LOGICAL COX_Z,COX_Y,COX_P,COX_L,COX_PMSL
      LOGICAL COX_O,COX_X
      LOGICAL COX_1234
      LOGICAL COX_LCASE_I
      LOGICAL COX_LCASE_C,COX_OCARB
      LOGICAL TRACER_O(0:18)

      CHARACTER*1 O_ASSM_FIELDS(6)

      COMMON/MODELO/
     & H_O_EWSPACE       ,H_O_NSSPACE      ,H_O_FIRSTLAT,
     & H_O_FIRSTLONG     ,H_O_POLELAT      ,H_O_POLELONG,
     & H_O_PTSPROW       ,N_COMP_O         ,  
     & H_NSIDEIMTO       ,H_NSIDEJMTO      ,
     & SEAICE_TYPE       ,OCEAN_BASINS     ,
     & COX_Z,COX_Y,COX_P ,COX_L,COX_1234   ,COX_PMSL    ,    
     & COX_O,COX_X,
     & COX_LCASE_I,
     & COX_LCASE_C,       COX_OCARB        ,
     & TRACER_O   ,       O_ASSM_FIELDS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER

      REAL    H_W_EWSPACE ,H_W_NSSPACE              
      REAL    H_W_FIRSTLAT,H_W_FIRSTLONG

      COMMON/MODELW/ 
     &              H_W_EWSPACE ,H_W_NSSPACE,              
     &              H_W_FIRSTLAT,H_W_FIRSTLONG

! Variables read in by namelist and used in SETMODL 
      INTEGER      OCAAA   ,OCAAO   ,OCAAW  
      INTEGER      NROWSO  ,NCOLSO  ,NLEVSO             
      INTEGER      NROWSW  ,NCOLSW
      INTEGER      NWTRAIN
      REAL         EWSPACEA,NSSPACEA    
      REAL         EWSPACEO,NSSPACEO 
      REAL         EWSPACEW,NSSPACEW
      REAL         FRSTLATA,FRSTLONA     
      REAL         FRSTLATO,FRSTLONO        
      REAL         FRSTLATW,FRSTLONW
  
      LOGICAL      ZonAvOzone
      INTEGER      IVDF
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      LWINC
      INTEGER      NECF(50)                                             
      INTEGER      OASFLDID(4) 
      INTEGER      OASLEV(6) ! dimensioned by max no of O-Assm groups
      INTEGER      OBAS
      INTEGER      OBS
      INTEGER      OCALB
      INTEGER      OCBOHaney
      INTEGER      OICE
      INTEGER      OIDYN
      INTEGER      OMP(4)
      REAL         POLELATA
      REAL         POLELONA
      REAL         POLELATO
      REAL         POLELONO
      INTEGER      PSA
      INTEGER      StLevGWdrag
      INTEGER      SWBND
      INTEGER      SWINC
      INTEGER      TCA(29)
      INTEGER      TCO(29)
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/   
     &RUN_TARGET_END,
     &OCAAA      ,EWSPACEA ,POLELATA ,FRSTLATA,LATS    ,           
     &            NSSPACEA ,POLELONA ,FRSTLONA,LONS    ,

     &OCAAO      ,EWSPACEO ,POLELATO ,FRSTLATO,NCOLSO  ,NLEVSO ,      
     &            NSSPACEO ,POLELONO ,FRSTLONO,NROWSO  ,

     &OCAAW      ,EWSPACEW ,          FRSTLATW,NCOLSW  ,
     &            NSSPACEW ,          FRSTLONW,NROWSW  ,NWTRAIN,

     &SWBND      ,LWBND    ,SWINC    ,LWINC   ,
     &ZonAvOzone ,AOBINC   ,
     &StLevGWdrag,AOBGRP   ,
     &BotVDiffLev, TopVDiffLev,
     &OCALB      ,TCA      ,
     &OIDYN      ,OBAS     ,OCBOHaney,OBS     ,OICE    ,IVDF ,
     &PSA        ,NECF     ,                                            
     &OASLEV     ,TCO      ,                                    
     &OMP        ,OASFLDID    

      CHARACTER*2  ATMOS_SR(0:NSECTP)
      CHARACTER*2  OCEAN_SR(0:NSECTP)
      CHARACTER*2  SLAB_SR (0:NSECTP)
      CHARACTER*2  WAVE_SR (0:NSECTP)
      CHARACTER*2  INDEP_SR(0:NSECTP)

      CHARACTER*1  BSPMSL
      CHARACTER*1  CCEW
      CHARACTER*1  FLOOR
      CHARACTER*1  IDO
      CHARACTER*1  LOSSM
      CHARACTER*1  MLMO
      CHARACTER*1  OAFLD(4)
      CHARACTER*1  OCARB
      CHARACTER*1  OROGR
      CHARACTER*1  OSFC
      CHARACTER*1  SCAL
      CHARACTER*1  SSTAnom
      CHARACTER*1  SWMCR
      CHARACTER*1  TOTAE
      CHARACTER*1  TOTEM
      CHARACTER*1  UPD175
      CHARACTER*1  MESO

      COMMON/STSHCHAR/   
     &BSPMSL, CCEW, INDEP_SR,
     &                   FLOOR  ,IDO   ,LOSSM   ,ATMOS_SR,      
     &MLMO    ,OAFLD    ,OCARB  ,OROGR ,OSFC    ,OCEAN_SR,     
     &SCAL    ,SSTAnom  ,SWMCR  ,TOTAE ,TOTEM   ,SLAB_SR ,   
     &UPD175  ,MESO     ,                        WAVE_SR           



      NAMELIST/STSHCOMP/
     &RUN_TARGET_END,
     &INDEP_SR    ,ATMOS_SR    ,OCEAN_SR ,SLAB_SR ,WAVE_SR,   
     &OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA,LATS   , 
     &             NSSPACEA    ,POLELONA ,FRSTLONA,LONS   ,
     &OCAAO       ,EWSPACEO    ,POLELATO ,FRSTLATO,NCOLSO ,NLEVSO   ,   
     &             NSSPACEO    ,POLELONO ,FRSTLONO,NROWSO ,
     &OCAAW       ,EWSPACEW    ,          FRSTLATW,NCOLSW ,
     &             NSSPACEW    ,          FRSTLONW,NROWSW ,
     &SWBND       ,LWBND       ,SWINC    ,LWINC   ,OROGR  ,            
     &ZonAvOzone  ,SWMCR       ,MESO     ,        
     &StLevGWdrag ,BotVDiffLev, TopVDiffLev,
     &OCALB       ,FLOOR       ,AOBINC   ,TOTAE   ,TOTEM  ,TCA  ,      
     &SSTAnom     ,SCAL        ,AOBGRP   ,        
     &NECF        ,BSPMSL      ,CCEW              ,UPD175 ,         
     &OIDYN       ,OBAS        ,OCBOHaney,OBS     ,OICE   ,IVDF ,IDO,   
     &OCARB       ,MLMO        ,PSA      ,OSFC    ,                     
     &LOSSM       ,OASLEV      ,OAFLD    ,TCO     ,                     
     &OMP         ,OASFLDID    
                           







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
C*L--------------------COMDECK  CANCILA ---------------------------
!
! Purpose : Contains index blocks for control of update of
!           ancillary fields.
!
! System component F0171
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  3.4   23/06/94  Update comments. D. Robinson
!  3.4   04/10/94  Increase NANCIL_FIELDS from 43 to 71. RTHBarnes
!  4.1   22/05/96  Move NANCIL_FIELDS to comdeck CANCMAXA. D. Robinson.
!  4.4   28/07/97  Add LAMIPII to common block. R A Stratton
!
! -------------------------------------------------------------------
C
C*L--------------------COMDECK  CANCMAXA  ---------------------------
!
! Purpose : Store maximum total no of atmosphere/slab ancillary fields.
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.1   20/05/96  New comdeck. Increase from 71 to 81. D. Robinson.
!  4.3   18/3/97  Increase from 81 to 82.  William Ingram.
!  4.4   16/9/97  Increase from 82 to 89.  Richard Betts.
!  4.5   05/05/98 Increase from 89 to 109. D.M. Goddard
!
! -------------------------------------------------------------------
! Type Declarations

      INTEGER NANCIL_FIELDS  ! No of Atmosphere & Slab Ancillary fields
      PARAMETER (NANCIL_FIELDS = 109) 

! Type Declarations

      INTEGER
     &  FILEANCIL,       ! File number associated with ancillary fields
     &  NLOOKUP,         ! Position of ancillary field in lookup tables.
     &  LOOKUP_STEP,     ! Interval between PP Headers refering to
C                        ! to the same ancillary fields at diferent time
     &  LEVELS,          ! Number of levels of data in each ancillary
C                        ! field.
     &  STASHANCIL,      ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block


      COMMON/IXANCILA/ FILEANCIL(NANCIL_FIELDS),
     &           NLOOKUP(NANCIL_FIELDS),
     &           LOOKUP_STEP(NANCIL_FIELDS),
     &           LEVELS(NANCIL_FIELDS),
     &           STASHANCIL(NANCIL_FIELDS),
     &           D1_ANCILADD(NANCIL_FIELDS)

C*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL
     &         UPDATE
     &,      L_SSTANOM          ! Indicator if SST anom to be formed
     &                          ! (RECON=T) or used (-DEF,RECON)
     & ,     LAMIPII            ! True if special AMIP II updating 

      INTEGER  FIELDCODE,
     &         STEPS
C*----------------------------------------------------------------------
      COMMON/CTANCILA/
     &                L_SSTANOM,LAMIPII,
     &         FIELDCODE(2,NANCIL_FIELDS),
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS)

!========================== COMDECK CSENARIO ===========================
!
!  Comdeck & COMMON for histories/scenarios of climate change forcings.
!
! Release   Changes:
!
!   4.3    New, in AWI1F403, author William Ingram, reviewer Cath Senior
!   4.5    Code for extra (H)(C)FCs added. J. M. Edwards 18/05/98
!
!=======================================================================

      INTEGER NCLMFCGS,          ! Number of such scenarios, made up of:
     &        NWMGHG,            ! Number of well-mixed greenhouse gases
     &        NSULPAT,           ! Number of sulphate loading patterns
     &        LENSCEN            ! Maximum length of scenarios
      PARAMETER ( NWMGHG = 9, NSULPAT = 2, LENSCEN = 50,
     &            NCLMFCGS = NWMGHG + NSULPAT )
      INTEGER S_CO2, S_CH4, S_N2O, S_CFC11, S_CFC12, S_SO4
     &        , S_CFC113, S_HCFC22, S_HFC125, S_HFC134A
!     !  Indices indicating which scenario corresponds to which forcing:
      PARAMETER ( S_CO2 = 1, S_CH4 = 2, S_N2O = 3, S_CFC11 = 4,
     &          S_CFC12 = 5,  S_SO4 = 6,     S_CFC113 = 8,
     &          S_HCFC22 = 9, S_HFC125 = 10, S_HFC134A = 11 )
!     !  Carbon dioxide (CO2), methane (CH4), nitrous oxide (N2O),
!     !  trichlorofluoromethane (CCl3F, "CFC-11"),
!     !  dichlorodifluoromethane (CCl2F2, "CFC-12"), and then the first
!     !  HadCM2-style anthropogenic sulphate loading pattern - these
!     !  come at the end as their number in principle may vary.
      INTEGER CLIM_FCG_YEARS(LENSCEN,NCLMFCGS),
!     ! Years at which a rate or level is specified
     &        CLIM_FCG_NYEARS(NCLMFCGS)
!     ! Number of such years, for each forcing
      REAL    CLIM_FCG_LEVLS(LENSCEN,NCLMFCGS),
     &        CLIM_FCG_RATES(LENSCEN,NCLMFCGS)
!     ! Values, or rates of increase, for the designated years.
!     !  See GAS_CALC (in Section 70) or the umui panels for details.
      COMMON / CSENARIO /  CLIM_FCG_NYEARS, CLIM_FCG_YEARS,
     &                     CLIM_FCG_LEVLS,  CLIM_FCG_RATES
      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
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


! Comdecks for ancillary files/fields.
C*L--------------------COMDECK  CANCFLDA ---------------------------
!
! Purpose : List of Ancillary Fields - Atmosphere and Slab
!           Stash Codes, Model Codes and Logical file Numbers
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.1   20/05/96  New comdeck. D. Robinson.
!  4.3   18/03/97  Add sulphate aerosol forcing.  William Ingram.
!  4.4   16/09/97  Add surface type initial fractions, initial leaf
!                  area index, initial canopy height, initial canopy
!                  conductance, vegetation disturbed fraction, soil 
!                  albedo and soil carbon content.  Richard Betts.
!  4.5   22/04/98  Add ammonia emission, surface soot emission and 
!                  elevated soot emission (68-70). R.Rawlins
!  4.5   05/05/98  Removal of Van-Genuchten 'B' parameter and
!                  increase number of user prognostics from 4 to 20
!                  Author D.M. Goddard 
!
! -------------------------------------------------------------------
! Type Declarations

      INTEGER
     &  ITEM_CODES_ANCIL(NANCIL_FIELDS)   ! Stash Codes
     &, MODEL_CODES_ANCIL(NANCIL_FIELDS)  ! Model Codes
     &, ANCIL_FILE_NO(NANCIL_FIELDS)      ! Logical file numbers

      DATA ITEM_CODES_ANCIL /
     &  30,  33,  34,  35,  36,  37,  60,  21,  23,  20,    !  1-10
     &  40,  41,  42,  43,  44,  45,  46,  47,  50,  51,    ! 11-20
     &  52,  53,  54,  55,  56,  26,  31,  24,  32,  28,    ! 21-30
     &  29,  93, 177,  19,  48,   9, 178, 179,  58,  59,    ! 31-40
     &  88,  87,  85,  57,  90,  17,  18, 301, 302, 303,    ! 41-50
     & 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,    ! 51-60
     & 314, 315, 316, 317, 318, 319, 320, 127, 128, 129,    ! 61-70
     &   0, 121, 122, 123, 124, 125, 126, 251, 207, 208,    ! 71-80
     & 209, 160, 216, 217, 218, 213, 219, 220, 223, 321,    ! 81-90
     & 322, 323, 324, 325, 326, 327, 328, 329, 330, 331,    ! 91-100
     & 332, 333, 334, 335, 336, 337, 338, 339, 340          ! 101-110
     & /

      DATA MODEL_CODES_ANCIL /
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    !  1-10
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 11-20
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 21-30
     &   1,   1,   3,   1,   1,   1,   3,   3,   1,   1,    ! 31-40
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 41-50
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 51-60
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 61-70
     &   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 71-80
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 81-90
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    ! 91-100
     &   1,   1,   1,   1,   1,   1,   1,   1,   1          ! 101-110
     & /

      DATA ANCIL_FILE_NO /
     &   9,  10,  10,  10,  10,  10,   1,   2,   2,   3,    !  1-10
     &   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,    ! 11-20
     &   5,   5,   5,   5,   5,   5,   7,   6,   7,   8,    ! 21-30
     &   8,   9,  11,  10,   4,   2,   6,   7,  12,  12,    ! 31-40
     &  13,  13,  13,  14,  14,  10,  10,  15,  15,  15,    ! 41-50
     &  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,    ! 51-60
     &  15,  15,  15,  15,  15,  15,  15,  12,  23,  23,    ! 61-70
     &   0,  17,  18,  18,  18,  18,  12,  24,   4,   5,    ! 71-80
     &   5,  19,  20,  21,  21,  21,  22,   4,   4,  16,    ! 81-90
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,    ! 91-100
     &  16,  16,  16,  16,  16,  16,  16,  16,  16          ! 101-110
     & /
 
CL External subroutines called:

      EXTERNAL
     &        FILE_OPEN,
     &        READ_FLH, READHEAD, SETPOS

CL Namelist input

      NAMELIST/ANCILCTA/L_SSTANOM,LAMIPII

!     UPANCA Namelist
      INTEGER
     &   ANC_REF_NO       ! Ancil Ref. No : See comdeck CANCFLDA
     &  ,PERIOD           ! Period of Updating Interval (Y/M/D/H)
     &  ,INTERVAL         ! Updating Interval

      NAMELIST /UPANCA/ ANC_REF_NO,PERIOD,INTERVAL

C Local Variables

      INTEGER
     &        I,               !
     &        ITEM,            !
     &        J,               !
     &        J1,              !
     &        K,               !
     &        LEN_IO,          !
     &        LOOKUPS,         !
     &        NFTIN,           ! Current FTN number for ancillary data
     &        START_BLOCK,     !
     &        JSOIL_DEPTHS     ! Pointer for soil thicknesses in level
C                              ! dep.consts.
     &       ,STASH_CODE       ! Stash item code
     &       ,NREC_A,NREC_S    ! No of atmos & slab records
     &       ,STASH_ADDR       ! Stash address
     &       ,DUMMY            !                                       
     &       ,N_ANC_UPD        ! No of ancillaries to be updated
      DATA DUMMY /1/                                                    

      CHARACTER*8 CPERIOD      ! PERIOD in characters.
      LOGICAL
     &        LFILE            !

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2)) .GT. (1.E-6*ABS(P1+P2)))

CL Internal Structure

      ICODE=0
      CMESSAGE=' '
      IOUNIT=0

C
CL  1.  Initialisation for atmosphere model

      DO I=1,NANCIL_FIELDS
        FILEANCIL(I) =ANCIL_FILE_NO(I)
        STASHANCIL(I)=ITEM_CODES_ANCIL(I)
      ENDDO


! Set default values

      L_SSTANOM=.FALSE.
      LAMIPII=.FALSE.

CL  Read in control information from namelist

        REWIND 5
      READ(5,ANCILCTA)

!     Initialise FIELDCODE from Namelist UPANCA
      N_ANC_UPD = 0
      DO I=1,NANCIL_FIELDS
        READ (5,UPANCA,ERR=101,END=101)
        FIELDCODE(1,ANC_REF_NO) = PERIOD
        FIELDCODE(2,ANC_REF_NO) = INTERVAL
        N_ANC_UPD = N_ANC_UPD+1
      ENDDO

 101  CONTINUE
      WRITE (6,*) ' '
      WRITE (6,*) N_ANC_UPD,' Atmos & Slab Ancillaries to be updated.'
      DO I=1,NANCIL_FIELDS
        IF (FIELDCODE(1,I).GT.0) THEN
        IF (FIELDCODE(1,I).EQ.1) CPERIOD=' Years'
        IF (FIELDCODE(1,I).EQ.2) CPERIOD=' Months'
        IF (FIELDCODE(1,I).EQ.3) CPERIOD=' Days'
        IF (FIELDCODE(1,I).EQ.4) CPERIOD=' Hours'
        WRITE (6,*) 'Anc Ref No ',I,' Stash code ',ITEM_CODES_ANCIL(I),
     &  ' Interval ',FIELDCODE(2,I),CPERIOD
        ENDIF
      ENDDO
      WRITE (6,*) ' '

! Check that ancillary field has valid address (>1) before proceding
!  to try and update it.  If not, switch off updating via FIELDCODE.
      DO I=1,NANCIL_FIELDS
        if (model_codes_ancil(i).eq.slab_im) then
          stash_addr = si_slab(stashancil(i))
         else
          stash_addr = si_atmos(stashancil(i))
        endif
        IF (stash_addr .le. 1) THEN
          IF (FIELDCODE(1,I).gt.0) THEN
           WRITE(6,*)' INANCILA: update requested for item ',i,
     &     ' STASHcode ',stashancil(i),' but prognostic address not set'
            WRITE(6,*)' FIELDCODE values reset to zeroes'
            FIELDCODE(1,I) = 0
            FIELDCODE(2,I) = 0
          END IF
        END IF
      END DO

CL  1.1 Set number of steps after which each ancillary field is updated
C       Zero is used for fields not to be updated

      DO I=1,NANCIL_FIELDS
        STEPS(I)=0
        IF (FIELDCODE(1,I).EQ.4)THEN
          STEPS(I)=FIELDCODE(2,I)*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I).EQ.3) THEN
          STEPS(I)=FIELDCODE(2,I)*24*STEPS_PER_HR
        END IF

      IF (LCAL360) THEN
        IF (FIELDCODE(1,I).EQ.2) THEN
          STEPS(I)=FIELDCODE(2,I)*30*24*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I).EQ.1) THEN
          STEPS(I)=FIELDCODE(2,I)*360*24*STEPS_PER_HR
        END IF
      ELSE
C Gregorian calender:
C If update interval is months or years, test each day. Further testing
C done in REPLANCA.

        IF (FIELDCODE(1,I).EQ.1.OR.FIELDCODE(1,I).EQ.2)THEN
         STEPS(I)=24*STEPS_PER_HR
        END IF
      END IF

      END DO

CL  1.2 Set master number of steps ANCILLARY_STEPS at which
CL      individual switches are tested.

C   Find first active field

      DO I=1,NANCIL_FIELDS
        IF (STEPS(I).GT.0) THEN
          ANCILLARY_STEPS=STEPS(I)
          GOTO 121
        END IF
      END DO

C No above fields found

      ANCILLARY_STEPS=0

      GOTO 900
121   ITEM=I

CL      Set ANCILLARY_STEPS to lowest common denominater of
CL      frequencies for active fields

      DO I=ITEM+1,NANCIL_FIELDS
        IF (STEPS(I).LT.ANCILLARY_STEPS
     *      .AND. STEPS(I).GT.0) THEN
          IF (MOD(ANCILLARY_STEPS,STEPS(I)).EQ.0) THEN
            ANCILLARY_STEPS=STEPS(I)
          ELSE
            J1=STEPS(I)-1
            DO J=J1,1,-1
              IF ((MOD(ANCILLARY_STEPS,J).EQ.0).AND.
     &           (MOD(STEPS(I),J).EQ.0)) THEN
                 GOTO 124
              ENDIF
            END DO
124         ANCILLARY_STEPS = J
          END IF
        END IF
      END DO

CL 1.2.4 Sea surface temperature must be updated when sea ice is update

      IF (STEPS(27).GT.0.AND.STEPS(28).LE.0) THEN
         STEPS(28)=1
      END IF


CL 1.3 Set number of headers for each ancillary field

      DO I=1,NANCIL_FIELDS
        LEVELS(I)=1
C   Multilayer hydrology
        IF(I.EQ.36)LEVELS(I)=SM_LEVELS
C   Multilayer aerosols
        IF(I.GE.41.AND.I.LE.43) LEVELS(I)=TR_LEVELS
C   Multilayer murk concentration and source
        IF(I.GE.44.AND.I.LE.45) LEVELS(I)=P_LEVELS
C   Multilayer user ancillaries
        IF(I.GE.90.AND.I.LE.109) LEVELS(I)=P_LEVELS
!   Multi-level ancillaries for sulphur cycle
        IF (I.EQ.72) LEVELS(I) = P_LEVELS
        IF (I.EQ.73) LEVELS(I) = P_LEVELS
        IF (I.EQ.74) LEVELS(I) = P_LEVELS
        IF (I.EQ.75) LEVELS(I) = P_LEVELS
        IF (I.EQ.76) LEVELS(I) = P_LEVELS
        IF (I.EQ.82) LEVELS(I) = NSULPAT
        IF (I.EQ.83) LEVELS(I) = NTYPE
        IF (I.EQ.84) LEVELS(I) = NPFT
        IF (I.EQ.85) LEVELS(I) = NPFT
      END DO

      LEVELS(7)=OZONE_LEVELS
      LEVELS(10)=ST_LEVELS


CL 1.4 Read headers

      LOOKUPS=0

      DO I=1,NDATASETS

C  Initialise LOOKUP_START (=0 implies file I not required)
        LOOKUP_START(I)=0

CL Check whether each physical file is needed

        LFILE=.FALSE.
        DO 141 J=1,NANCIL_FIELDS


          IF (FILEANCIL(J).EQ.I.AND.STEPS(J).GT.0) THEN


            LFILE=.TRUE.
          END IF
141     CONTINUE

        IF(LFILE) THEN

      WRITE(6,*) ' '
      WRITE(6,*) ' Ancillary data file ',I,', unit no ',FTNANCIL(I),
     &           ', ',TITLE(I)

CL Read headers for physical files required

          NFTIN=FTNANCIL(I)

CL 1.4.1 Buffer in fixed length header record


        CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),
     &                 LEN_FT_ENVIR(NFTIN),0,0,ICODE)
        IF(ICODE.NE.0)THEN
          CMESSAGE='INANCLA: Error opening file'
          write(6,*) 'INANCILA: Error opening file on unit ',NFTIN,
     &               ' accessed from env.var.: ',FT_ENVIRON(NFTIN)
          RETURN
        ENDIF
        CALL SETPOS(NFTIN,0,ICODE)

C       Read in fixed header to get array dimensions
        CALL READ_FLH(NFTIN,FIXHD(1,I),LEN_FIXHD,ICODE,CMESSAGE)
        IF (ICODE.GT.0) THEN
          WRITE (6,*) ' Error in reading fixed header for file ',I
          GO TO 9999   !  Return
        ENDIF

C       Check for negative dimensions
        IF (FIXHD(101,I).LE.0) FIXHD(101,I)=1
        IF (FIXHD(106,I).LE.0) FIXHD(106,I)=1
        IF (FIXHD(111,I).LE.0) FIXHD(111,I)=1
        IF (FIXHD(112,I).LE.0) FIXHD(112,I)=1
        IF (FIXHD(151,I).LE.0) FIXHD(151,I)=1
        IF (FIXHD(152,I).LE.0) FIXHD(152,I)=1
        IF (FIXHD(161,I).LE.0) FIXHD(161,I)=1

C Set start position of boundary fields for file                       
        LOOKUP_START(I)=LOOKUPS+1 
                                     
        IF (LOOKUPS+FIXHD(152,I).GT.NLOOKUPS) THEN              
          WRITE (6,*) 'No room in LOOKUP table for Ancillary File ',I   
          CMESSAGE='INANCILA: Insufficient space for LOOKUP headers'   
          ICODE=14                                              
          GO TO 9999   !  Return
        END IF                                                      

        CALL SETPOS(NFTIN,0,ICODE)
        IF (ICODE.GT.0) THEN
          WRITE (6,*) ' ERROR in SETPOS called from INANCA1A'
          WRITE (6,*) ' SETPOS attempted with Unit No ',NFTIN
          CMESSAGE = 'INANCA1A : ERROR in SETPOS'
          GO TO 9999    !   Return
        ENDIF
       
        CALL READHEAD(NFTIN,
     &                FIXHD(1,I),LEN_FIXHD,
     &                INTHD(1,I),FIXHD(101,I),
     &                REALHD(1,I),FIXHD(106,I),
     &                LEVDEPC,FIXHD(111,I),FIXHD(112,I),
     &                DUMMY,DUMMY,DUMMY,
     &                DUMMY,DUMMY,DUMMY,
     &                DUMMY,DUMMY,DUMMY,
     &                DUMMY,DUMMY,
     &                DUMMY,DUMMY,
     &                DUMMY,DUMMY,
     &                DUMMY,DUMMY,
     &                DUMMY,DUMMY,
     &                LOOKUP(1,LOOKUPS+1),FIXHD(151,I),FIXHD(152,I),
     &                FIXHD(161,I),
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &                START_BLOCK,ICODE,CMESSAGE)

        IF (ICODE.GT.0) THEN
           WRITE(6,*) 'ERROR in READHEAD for Ancillary File ',I
           WRITE(6,*) 'Unit Number ',NFTIN
           GO TO 9999   !   Return
        ENDIF

!     Check calendar indicator
        IF ((     LCAL360 .and. FIXHD(8,I).NE.2) .or.
     &      (.not.LCAL360 .and. FIXHD(8,I).NE.1) ) THEN   
          ICODE=100+I
          CMESSAGE='INANCILA : Wrong calendar set in Ancillary File'
          WRITE (6,*) ' ******** Error in INANCILA ********'
          WRITE (6,*) ' Wrong calendar setting in Ancillary File ',I
          IF (LCAL360) THEN
            WRITE (6,*) ' Model run is set up for 360 day calendar.'
            WRITE (6,*) ' Ancillary File is for 365 day calendar.'
          ELSE
            WRITE (6,*) ' Model run is set up for 365 day calendar.'
            WRITE (6,*) ' Ancillary File is for 360 day calendar.'
          ENDIF
          WRITE (6,*) ' Rerun with correct ancillary file.'
          GO TO 9999   !  Return
        ENDIF

        FILE_LEVELS=1

        IF(I.EQ.1) THEN
          FILE_LEVELS=OZONE_LEVELS
        ELSE IF(I.EQ.2) THEN
          FILE_LEVELS=SM_LEVELS
C This is the maximum value that might be present on the ancillary
C file if it includes soil moisture in layers; otherwise only single
C level data is present and PR_FIXHD will not check value since
C FIXHD(110) will be zero
        ELSE IF(I.EQ.3) THEN
            FILE_LEVELS=ST_LEVELS
        ELSE IF(I.EQ.13) THEN   ! for multilevel aerosols
            FILE_LEVELS=TR_LEVELS
        ELSE IF(I.EQ.14.or.I.EQ.16) THEN   ! for murk and user ancil.
            FILE_LEVELS=P_LEVELS
        ELSE IF(I.EQ.17.or.I.EQ.18) THEN
!           multi-level sulphur cycle ancillary files.
            FILE_LEVELS=P_LEVELS
        END IF


CL 1.4.2 Buffer in integer constants

           IF(FIXHD(100,I).GT.0) THEN

C Check for error in file pointers

C Check validity of integer data and print out information
C All files except ozone should contain full fields

            IF(INTHD(6,I).NE.ROW_LENGTH) THEN
C Ozone may contain zonal mean data
              IF(I.NE.1.OR.INTHD(6,I).NE.1) THEN
                ICODE=4
                CMESSAGE='INANCILA:integer header error'
                WRITE(6,*) ' INTHD(6) : ',INTHD(6,I),' ??'
                RETURN
              END IF
            END IF

            IF(INTHD(7,I).NE.P_ROWS.AND.(I.EQ.9.AND.INTHD
     &        (7,I).NE.U_ROWS)) THEN
              ICODE=5
              CMESSAGE='INANCILA:integer header error'
              WRITE(6,*) ' INTHD(7) : ',INTHD(7,I),' ??'
              RETURN
            END IF

            IF (I.EQ.1) THEN   !  Ozone file
              WRITE (6,*) ' '
              IF (INTHD(6,I).EQ.1)THEN
                WRITE (6,*) ' OZONE file contains zonal mean data for ',
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ELSEIF (INTHD(6,I).EQ.ROW_LENGTH)THEN
                WRITE (6,*) ' OZONE file contains full fields for ',
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ENDIF
! Check that correct ozone file has been provided.
              IF (ZonAvOzone) THEN
                IF (INTHD(6,I).NE.1) THEN
                  WRITE (6,*) ' Zonal Ozone Data is expected',
     &            ' for 1 point x ',P_ROWS,' rows'
                  ICODE = 51
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ELSE
                IF (INTHD(6,I).NE.ROW_LENGTH) THEN
                  WRITE (6,*) ' Ozone Data is expected for ',
     &            ROW_LENGTH,' points x ',P_ROWS,' rows.'
                  ICODE = 52
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ENDIF
            ENDIF

          END IF

CL 1.4.3 Buffer in real constants

          IF(FIXHD(105,I).GT.0) THEN

C Check validity of real header and print out information

           DO J=1,6
             IF(REALHD(J,I).GT.(A_REALHD(J)+0.1).OR.
     &         REALHD(J,I).LT.(A_REALHD(J)-0.1))THEN
             IF(I.NE.1.OR.(J.NE.1.AND.J.NE.4))THEN
               WRITE(6,*)(REALHD(K,I),K=1,6),(A_REALHD(K),K=1,6)
               ICODE=8
               CMESSAGE='INANCILA: REAL header Error.'
               RETURN
             END IF
             END IF
           END DO

         END IF

CL 1.4.4 Buffer in level dependent constants if required
C        Not retained in model after initial check

         IF(FIXHD(110,I).GT.0) THEN

CL Only files 1 (Ozone), and 3 (Soil temperature)should contain multi
CL level data. File 2 (Soil moisture,snow depth,fractional snow time
CL and soil moisture in layers) may possibly also have multi level data.
CL FILES 13,14,16 (aerosols, murkiness, user ancil.) may also have
CL  multi level data.

CL If ozone file, check against model levels

           IF(I.EQ.1) THEN
             DO J=1,OZONE_LEVELS
               DO J1=1,4
                 IF(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC
     &                   (J+P_LEVELS-OZONE_LEVELS,J1))) THEN
      WRITE(6,*)'Error in level dependent constants:Level=',J
                   WRITE(6,*)'Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC
     &                   (J+P_LEVELS-OZONE_LEVELS,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=11
               CMESSAGE='INANCILA: error in LEVDEPC.'
                   RETURN
                 END IF
               END DO
             END DO

           ELSE IF (I.EQ.2.OR.I.EQ.3) THEN

             IF (A_FIXHD(12).LE.0) THEN
               ICODE = 121
               CMESSAGE = 'INANCA1A : FIXHD(12) not set in A_FIXHD'
               WRITE (6,*) ' FIXHD(12) not set in A_FIXHD.'
               WRITE (6,*) ' Run reconfiguration program to set.'
               GO TO 9999    !   Return
             ELSEIF (A_FIXHD(12).LT.305) THEN
               JSOIL_DEPTHS = 13
             ELSE
               JSOIL_DEPTHS = 6
             ENDIF

C soil moisture levels
C
C  If deep soil temperatures or multilayer soil moistures, check
C  against model soil level/layer depths.
C
             IF (I.EQ.2) THEN
               DO J=1,SM_LEVELS
                 IF (LNER(LEVDEPC(J),A_LEVDEPC(J,JSOIL_DEPTHS))) THEN 
                   ICODE=12
                   CMESSAGE='INANCILA: error in LEVDEPC.'
                  RETURN
                 END IF
               END DO
             END IF
             IF (I.EQ.3) THEN
               DO J=1,ST_LEVELS
C Penman-Monteith BL version
                 IF(ST_LEVELS.EQ.SM_LEVELS)THEN
                   IF (LNER(LEVDEPC(J),A_LEVDEPC(J,JSOIL_DEPTHS)))THEN
                     ICODE=12
                     CMESSAGE='INANCILA: error in LEVDEPC.'
                     RETURN
                   ENDIF
C All other BL versions
                 ELSE
                   IF (LNER(LEVDEPC(J),A_LEVDEPC(J+1,JSOIL_DEPTHS)))THEN
                     ICODE=12
                     CMESSAGE='INANCILA: error in LEVDEPC.'
                     RETURN
                   END IF
                 END IF
               END DO
             END IF

CL If aerosol file, check against model levels

           ELSE IF (I.EQ.13) THEN
             DO J=1,TR_LEVELS
               DO J1=1,4
                 IF(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC
     &                   (J,J1))) THEN
      WRITE(6,*)'Error in level dependent constants:Level=',J
                   WRITE(6,*)'Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC(J,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=16
               CMESSAGE='INANCILA: error in LEVDEPC.'
                   RETURN
                 END IF
               END DO
             END DO

CL If murk or user ancillary file, check against model levels

           ELSE IF (I.EQ.14.or.I.EQ.16) THEN
             DO J=1,P_LEVELS
               DO J1=1,4
                 IF(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC
     &                   (J,J1))) THEN
      WRITE(6,*)'Error in level dependent constants:Level=',J
                   WRITE(6,*)'Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC(J,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=16
               CMESSAGE='INANCILA: error in LEVDEPC.'
                   RETURN
                 END IF
               END DO
             END DO

           END IF

         END IF

CL 1.4.5 Buffer in lookup table
C Set start position of boundary fields for file

         IF(FIXHD(150,I).GT.0) THEN


           NREC_A = 0
           NREC_S = 0
           DO J = 1,FIXHD(152,I)
             IF (LOOKUP(MODEL_CODE,LOOKUPS+J) .eq. 0 .or.
     &           LOOKUP(MODEL_CODE,LOOKUPS+J) .eq. imdi) THEN
               STASH_CODE = LOOKUP(ITEM_CODE,LOOKUPS+J)
               IF ((STASH_CODE.GE.177 .and. STASH_CODE.LE.179) .or.
     &             (STASH_CODE.GE.210 .and. STASH_CODE.LE.212)) THEN
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = slab_im
                 NREC_S = NREC_S+1
               ELSE
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = atmos_im
                 NREC_A = NREC_A+1
               END IF
             END IF
           END DO
           IF (NREC_A.GT.0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_A,
     &       ' records set to atmos_im in ancillary file ',I
           ENDIF
           IF (NREC_S.GT.0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_S,
     &       ' records set to slab_im in ancillary file ',I
           ENDIF

         END IF

         LOOKUPS=LOOKUPS+FIXHD(152,I)

       ELSE

CL  If file not required, zero fixed length header
         DO J=1,LEN_FIXHD
      FIXHD(J,I)=0
         END DO

         LOOKUP_START(I)=LOOKUPS+1
       END IF

      END DO

CL 1.5 Set positions in main data blocks


      DO I=1,NANCIL_FIELDS
        IF (MODEL_CODES_ANCIL(I).EQ.SLAB_IM) THEN
          D1_ANCILADD(I)=SI_SLAB(STASHANCIL(I))
        ELSE
          D1_ANCILADD(I)=SI_ATMOS(STASHANCIL(I))
        ENDIF
      ENDDO

CL 1.51 If a request is made to update a field, ensure that space for
CL     that field has been allocted in D1.

      DO I=1,NANCIL_FIELDS
        IF((FIELDCODE(1,I).GT.0).AND.(D1_ANCILADD(I).LE.1)) THEN
          WRITE(6,*)' An address in D1 has not been set for ancillary
     & field number ',I
          ICODE=30
          CMESSAGE='INANCILA: updating for ancillary field is requested
     & but no space has been allocated in D1'
          RETURN
        ENDIF
      END DO


CL 1.52 Reset target STASH codes for fields mapped onto other fields
CL      (eg. SLAB ref SST/ice depth are obtained from normal SST and
CL           ice depth ancillary files as used in atmos-only mode).

      STASHANCIL(37)=STASHANCIL(28)   ! SLAB ref SST    == Normal SST
      STASHANCIL(38)=STASHANCIL(29)   ! SLAB ref icedep == Normal icedep

CL 1.6 Set positions of data

      DO I=1,NANCIL_FIELDS
      NLOOKUP(I) =0
      LOOKUP_STEP(I)=0

C If LOOKUP_START=0 for file FILEANCIL(I), no fields required.
        IF (LOOKUP_START(FILEANCIL(I)).GT.0) THEN

        DO J=LOOKUP_START(FILEANCIL(I)),LOOKUPS

          IF (LOOKUP(ITEM_CODE,J).EQ.STASHANCIL(I)) THEN
            NLOOKUP(I)=J-LOOKUP_START(FILEANCIL(I))+1
            GOTO 161
          END IF

        END DO

C Find second occurence of data to set LOOKUP_STEP

161     LOOKUP_STEP(I)=0


        IF(J.LT.LOOKUPS) THEN

          DO J1=J+LEVELS(I),LOOKUPS
            IF (LOOKUP(ITEM_CODE,J1).EQ.STASHANCIL(I)) THEN
              LOOKUP_STEP(I)=J1-NLOOKUP(I)-LOOKUP_START(FILEANCIL(I))+1
              GOTO 164
            END IF
          END DO
164      CONTINUE
        END IF

        END IF

      END DO

CL SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
CL fractional time fields

      LEVELS(27)=2
      LEVELS(9)=2

 900  CONTINUE
 9999 CONTINUE
      RETURN
      END
