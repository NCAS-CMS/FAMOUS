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
CLL  SUBROUTINE READHEAD---------------------------------------
CLL
CLL AD, DR      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.2  06/05/93    Skip call to CHKLOOK if PP type file
CLL                    Author: A. Dickinson    Reviewer: D. Richardson
CLL
CLL  3.1   22/12/92     Allow use by ancillary field headers
CLL                     Author A. Dickinson    Reviewer C. Wilson
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL  3.2   12/05/93     Adapt to read prognostic fields only.
CLL                     Author D. Robinson     Reviewer A. Dickinson
CLL   3.5    28/03/95 MPP code: New code for parallel I/O
CLL                                              P.Burton
!     3.5  21/06/95  Set lookup(45) if initial dump pre version 3.5
!                    Author D.M.Goddard    Reviewer S Swarbrick
!     4.0  06/10/95  Set variable MODEL for all diagnostics in dump
!                    Author D.M. Goddard
!     4.1  23/05/96  Removed resetting of FIXHD(161) for MPP code
!                    P.Burton
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.

CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  Logical component: R30
CLL
CLL  System task: F3
CLL
CLL  Purpose: Reads in model dump header records on unit NFTIN and
CLL           checks model and dump dimensions for consistency.
CLL
CLL  Documentation: Unified Model Documentation Paper No F3
CLL                 Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE READHEAD(NFTIN,FIXHD,LEN_FIXHD,       ! Intent (In)
     &                    INTHD,LEN_INTHD,
     &                    REALHD,LEN_REALHD,
     &                    LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,
     &                    ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,
     &                    COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,
     &                    FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,
     &                    EXTCNST,LEN_EXTCNST,
     &                    DUMPHIST,LEN_DUMPHIST,
     &                    CFI1,LEN_CFI1,
     &                    CFI2,LEN_CFI2,
     &                    CFI3,LEN_CFI3,
     &                    LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,
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
     &                    START_BLOCK,ICODE,CMESSAGE)   ! Intent (Out)

      IMPLICIT NONE

      INTEGER
     * NFTIN         !IN Unit no of dump
     *,LEN_FIXHD     !IN Length of fixed length header
     *,LEN_INTHD     !IN Length of integer header
     *,LEN_REALHD    !IN Length of real header
     *,LEN1_LEVDEPC  !IN 1st dim of level dep consts
     *,LEN2_LEVDEPC  !IN 2ndt dim of level dep consts
     *,LEN1_ROWDEPC  !IN 1st dim of row dep consts
     *,LEN2_ROWDEPC  !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC  !IN 1st dim of column dep consts
     &,LEN2_COLDEPC  !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC  !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC  !IN 2nd dim of field dep consts
     &,LEN_EXTCNST   !IN Length of extra constants
     &,LEN_DUMPHIST  !IN Length of history block
     &,LEN_CFI1      !IN Length of comp field index 1
     &,LEN_CFI2      !IN Length of comp field index 2
     &,LEN_CFI3      !IN Length of comp field index 3
     &,LEN1_LOOKUP   !IN 1st dim of lookup
     &,LEN2_LOOKUP   !IN 2nd dim of lookup

      INTEGER
     * LEN_DATA       !IN Length of model data
     *,START_BLOCK    !OUT Pointer to position of each block.
     *                !Should point to start of model data block on exit
     *,ICODE          !OUT Return code; successful=0
     *                !                 error > 0

      CHARACTER*(80)
     * CMESSAGE       !OUT Error message if ICODE > 0

      INTEGER
     * FIXHD(LEN_FIXHD) !IN Fixed length header
     *,INTHD(LEN_INTHD) !IN Integer header
     *,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) !IN PP lookup tables
     *,CFI1(LEN_CFI1+1) !IN Compressed field index no 1
     *,CFI2(LEN_CFI2+1) !IN Compressed field index no 2
     *,CFI3(LEN_CFI3+1) !IN Compressed field index no 3

      REAL
     & REALHD(LEN_REALHD) !IN Real header
     &,LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC) !IN Lev dep consts
     &,ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC) !IN Row dep consts
     &,COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC) !IN Col dep consts
     &,FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC) !IN Field dep consts
     &,EXTCNST(LEN_EXTCNST+1)   !IN Extra constants
     &,DUMPHIST(LEN_DUMPHIST+1) !IN History block

C Local arrays:------------------------------------------------
C None
C -------------------------------------------------------------
C External subroutines called:---------------------------------
      EXTERNAL IOERROR,POSERROR,PR_FIXHD,CHK_LOOK,BUFFIN
C*-------------------------------------------------------------
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
C Local variables:---------------------------------------------
      INTEGER K
      INTEGER LEN_IO
      INTEGER FIXHD_152    !  Original value of FIXHD(152)
      LOGICAL L_A_DUMP
      LOGICAL L_O_DUMP
      REAL A
C -------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

CL 1. Buffer in fixed length header record

      CALL BUFFIN(NFTIN,FIXHD(1),LEN_FIXHD,LEN_IO,A)


C Check for I/O errors
      IF(A.NE.-1.0.OR.LEN_IO.NE.LEN_FIXHD)THEN
        CALL IOERROR('buffer in of fixed length header',A,LEN_IO
     *               ,LEN_FIXHD)
        CMESSAGE='READHEAD: I/O error'
        ICODE=1
        RETURN
      ENDIF

      START_BLOCK=LEN_FIXHD+1

      FIXHD_152 = FIXHD(152)    !  Store original value

C     Test if atmos dump read in
      L_A_DUMP = FIXHD(5).EQ.1 .AND. FIXHD(2).EQ.1
     *           . AND . LEN_DATA.NE.IMDI

C     Test if ocean dump read in
      L_O_DUMP = FIXHD(5).EQ.1 .AND. FIXHD(2).EQ.2
     *           . AND . LEN_DATA.NE.IMDI

      IF (L_A_DUMP .OR. L_O_DUMP) THEN
        IF (FIXHD(152).NE.LEN2_LOOKUP) THEN
CXX       WRITE (6,*) 'FIXHD(152) being reset from ',FIXHD(152),' to ',
CXX  *    LEN2_LOOKUP
          FIXHD(152) = LEN2_LOOKUP
        ENDIF
      ENDIF

C Check validity of data and print out fixed header information

      IF (mype .EQ. 0) THEN
      CALL PR_FIXHD(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC
     *,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC
     *,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1
     *,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA
     *,ICODE,CMESSAGE)

      IF(ICODE.GT.0)RETURN

      ENDIF
CL 2. Buffer in integer constants

      IF(FIXHD(100).GT.0)THEN

C Check for error in file pointers
       IF(FIXHD(100).NE.START_BLOCK)THEN
        CALL POSERROR('integer constants',START_BLOCK,100,FIXHD(100))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=2
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,INTHD(1),FIXHD(101),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(101))THEN
        CALL IOERROR('buffer in of integer constants',A,LEN_IO,
     *               FIXHD(101))
        CMESSAGE='READHEAD: I/O error'
        ICODE=3
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(101)

      ENDIF

CL 3. Buffer in real constants

      IF(FIXHD(105).GT.0)THEN

C Check for error in file pointers
       IF(FIXHD(105).NE.START_BLOCK)THEN
        CALL POSERROR('real constants',START_BLOCK,105,FIXHD(105))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=4
        RETURN
       ENDIF

C Check for I/O errors
      CALL BUFFIN(NFTIN,REALHD(1),FIXHD(106),LEN_IO,A)

       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(106))THEN
        CALL IOERROR('buffer in of real constants',A,LEN_IO,
     *                FIXHD(106))
        CMESSAGE='READHEAD: I/O error'
        ICODE=5
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(106)


      ENDIF

CL 4. Buffer in level dependent constants

      IF(FIXHD(110).GT.0.AND.LEN1_LEVDEPC.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(110).NE.START_BLOCK)THEN
        CALL POSERROR('level dependent constants',
     *  START_BLOCK,110,FIXHD(110))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=6
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,LEVDEPC(1),FIXHD(111)*FIXHD(112),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(111)*FIXHD(112))THEN
        CALL IOERROR('buffer in of level dependent constants',A,LEN_IO,
     *               FIXHD(111)*FIXHD(112))
        CMESSAGE='READHEAD: I/O error'
        ICODE=7
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(111)*FIXHD(112)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(111)*FIXHD(112)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 5. Buffer in row dependent constants

      IF(FIXHD(115).GT.0.AND.LEN1_ROWDEPC.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(115).NE.START_BLOCK)THEN
        CALL POSERROR('row dependent constants',
     *  START_BLOCK,115,FIXHD(115))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=8
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,ROWDEPC(1),FIXHD(116)*FIXHD(117),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(116)*FIXHD(117))THEN
        CALL IOERROR('buffer in of row dependent constants',A,LEN_IO,
     *                FIXHD(116)*FIXHD(117))
        CMESSAGE='READHEAD: I/O error'
        ICODE=9
        RETURN
      ENDIF


       START_BLOCK=START_BLOCK+FIXHD(116)*FIXHD(117)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(116)*FIXHD(117)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 6. Buffer in column dependent constants

      IF(FIXHD(120).GT.0.AND.LEN1_COLDEPC.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(120).NE.START_BLOCK)THEN
        CALL POSERROR('column dependent constants',
     *  START_BLOCK,120,FIXHD(120))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=10
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,COLDEPC(1),FIXHD(121)*FIXHD(122),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(121)*FIXHD(122))THEN
        CALL IOERROR('buffer in of column dependent constants',A,LEN_IO,
     *               FIXHD(121)*FIXHD(122))
        CMESSAGE='READHEAD: I/O error'
        ICODE=11
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(121)*FIXHD(122)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(121)*FIXHD(122)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 7. Buffer in constants stored as fields

      IF(FIXHD(125).GT.0.AND.LEN1_FLDDEPC.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(125).NE.START_BLOCK)THEN
        CALL POSERROR('fields of constants',
     *  START_BLOCK,125,FIXHD(125))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=12
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,FLDDEPC(1),FIXHD(126)*FIXHD(127),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(126)*FIXHD(127))THEN
        CALL IOERROR('buffer in of field dependent constants',A,LEN_IO,
     *               FIXHD(126)*FIXHD(127))
        CMESSAGE='READHEAD: I/O error'
        ICODE=13
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(126)*FIXHD(127)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(126)*FIXHD(127)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 8. Buffer in extra constants

      IF(FIXHD(130).GT.0.AND.LEN_EXTCNST.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(130).NE.START_BLOCK)THEN
        CALL POSERROR('extra constants',
     *  START_BLOCK,130,FIXHD(130))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=14
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,EXTCNST(1),FIXHD(131),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(131))THEN
        CALL IOERROR('buffer in extra constants',A,LEN_IO,
     *               FIXHD(131))
        CMESSAGE='READHEAD: I/O error'
        ICODE=15
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(131)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' EXTRA CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(131)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 9. Buffer in temporary history block

      IF(FIXHD(135).GT.0.AND.LEN_DUMPHIST.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(135).NE.START_BLOCK)THEN
        CALL POSERROR('history',
     *  START_BLOCK,136,FIXHD(136))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=16
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,DUMPHIST(1),FIXHD(136),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(136))THEN
        CALL IOERROR('buffer in of history file',A,LEN_IO,
     *               FIXHD(136))
        CMESSAGE='READHEAD: I/O error'
        ICODE=17
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(136)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(136)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 10. Buffer in compressed field index1

      IF(FIXHD(140).GT.0.AND.LEN_CFI2.NE.0)THEN

C Check for error in file pointers

       IF(FIXHD(140).NE.START_BLOCK)THEN
        CALL POSERROR('compressed field index1',
     *  START_BLOCK,140,FIXHD(140))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=18
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,CFI1(1),FIXHD(141),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(141))THEN
        CALL IOERROR('buffer in of compressed index1',A,LEN_IO,
     *               FIXHD(141))
        CMESSAGE='READHEAD: I/O error'
        ICODE=19
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(141)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(141)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 11. Buffer in compressed field index2

      IF(FIXHD(142).GT.0.AND.LEN_CFI2.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(142).NE.START_BLOCK)THEN
        CALL POSERROR('compressed field index2',
     *  START_BLOCK,142,FIXHD(142))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=20
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,CFI2(1),FIXHD(143),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(143))THEN
       CALL IOERROR('buffer in of compressed index2',A,LEN_IO,
     *               FIXHD(143))
        CMESSAGE='READHEAD: I/O error'
        ICODE=21
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(143)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(143)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 12. Buffer in compressed field index3

      IF(FIXHD(144).GT.0.AND.LEN_CFI3.NE.0)THEN

C Check for error in file pointers
       IF(FIXHD(144).NE.START_BLOCK)THEN
        CALL POSERROR('compressed field index3',
     *  START_BLOCK,144,FIXHD(144))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=22
        RETURN
       ENDIF

      CALL BUFFIN(NFTIN,CFI3(1),FIXHD(145),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(145))THEN
        CALL IOERROR('buffer in of compressed index3',A,LEN_IO,
     *               FIXHD(145))
        CMESSAGE='READHEAD: I/O error'
        ICODE=23
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(145)

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(145)

      ENDIF ! mype .EQ. 0
      ENDIF

CL 13. Buffer in lookup table

      IF(FIXHD(150).GT.0)THEN

C Supress checking if not full dump
      IF(LEN_DUMPHIST.NE.0)THEN
C Check for error in file pointers
       IF(FIXHD(150).NE.START_BLOCK)THEN
        CALL POSERROR('lookup table',
     *  START_BLOCK,150,FIXHD(150))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=24
        RETURN
       ENDIF
      ENDIF

C Move to start of Look Up Table
      CALL SETPOS(NFTIN,FIXHD(150)-1,ICODE)

C Read in fields from LOOKUP table
      CALL BUFFIN(NFTIN,LOOKUP(1,1),FIXHD(151)*FIXHD(152),LEN_IO,A)

C Check for I/O errors
       IF(A.NE.-1.0.OR.LEN_IO.NE.FIXHD(151)*FIXHD(152))THEN
        CALL IOERROR('buffer in of lookup table',A,LEN_IO,
     *               FIXHD(151)*FIXHD(152))
        CMESSAGE='READHEAD: I/O error'
        ICODE=25
        RETURN
       ENDIF

C Point to start of data section ( Use original FIXHD(152) )
       START_BLOCK=START_BLOCK+FIXHD(151)*FIXHD_152

       IF (mype .EQ. 0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LOOKUP TABLE'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(151)*FIXHD(152)

       IF (FIXHD(152).LT.FIXHD_152) THEN
         WRITE(6,'('' '')')
         WRITE(6,'('' '',I6,'' Entries in Look Up Table.'')') FIXHD_152
         WRITE(6,'('' '',I6,'' Entries read in.'')') FIXHD(152)
       ENDIF

      ENDIF ! mype .EQ. 0
!---------------------------------------------------------------
! Reset LOOKUP(45) if not set
!---------------------------------------------------------------

        DO K=1,LEN2_LOOKUP
          IF(LOOKUP(45,K).EQ.0.OR.LOOKUP(45,K).EQ.IMDI)THEN 

!Section 0: Prognostic fields.
          IF(LOOKUP(42,K).LE.100.OR.
     &      (LOOKUP(42,K).GE.200.AND.LOOKUP(42,K).LE.205))THEN
            LOOKUP(45,K)=1

          ELSE IF((LOOKUP(42,K).GT.100.AND.LOOKUP(42,K).LE.176).OR.
     &            (LOOKUP(42,K).GE.180.AND.LOOKUP(42,K).LT.200))THEN
            LOOKUP(45,K)=2

          ELSE IF((LOOKUP(42,K).GE.177.AND.LOOKUP(42,K).LE.179).OR.
     &            (LOOKUP(42,K).GE.210.AND.LOOKUP(42,K).LE.212))THEN
            LOOKUP(45,K)=3
          
! Sections 1 - 99: Diagnostic fields
          ELSE IF(LOOKUP(42,K).GE.1000.AND.LOOKUP(42,K).LE.29999)THEN
            IF((LOOKUP(42,K).GE.21177.AND.LOOKUP(42,K).LE.21179).OR.
     &         (LOOKUP(42,K).GE.21225.AND.LOOKUP(42,K).LE.21227).OR.
     &         (LOOKUP(42,K).GE.22177.AND.LOOKUP(42,K).LE.22179).OR.
     &         (LOOKUP(42,K).GE.22225.AND.LOOKUP(42,K).LE.22227).OR.
     &         (LOOKUP(42,K).GE.23177.AND.LOOKUP(42,K).LE.23179).OR.
     &         (LOOKUP(42,K).GE.23225.AND.LOOKUP(42,K).LE.23227).OR.
     &         (LOOKUP(42,K).GE.24177.AND.LOOKUP(42,K).LE.24179).OR.
     &         (LOOKUP(42,K).GE.24225.AND.LOOKUP(42,K).LE.24227))THEN
              LOOKUP(45,K)=3        !Slab diagnostic

            ELSE
              LOOKUP(45,K)=1        !Atmosphere diagnostic

            END IF

          ELSE IF(LOOKUP(42,K).GE.30000.AND.LOOKUP(42,K).LE.99999)THEN
            IF(LOOKUP(42,K).GE.40000.AND.LOOKUP(42,K).LE.40999)THEN
              LOOKUP(45,K)=3        !Slab diagnostic

            ELSE
              LOOKUP(45,K)=2        !Ocean diagnostic

            END IF

          ELSE
            WRITE(6,*) 'WARNING: User defined field found - ',
     &                 'STASH code : ', LOOKUP(42,K)
            WRITE(6,*) ' Internal model number can not be defined.'
            WRITE(6,*) ' Setting internal model number to atmosphere.'
            LOOKUP(45,K)=1

          ENDIF  

        ENDIF 

      ENDDO 
C---------------------------------------------------------------
C  Reset LOOKUP headers if dump created earlier than vn2.8
C---------------------------------------------------------------

      IF(FIXHD(12).LT.208)THEN
        CALL NEWPACK(LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP)
      ENDIF

C Check LOOKUP for consistency with PARAMETER statements
      IF(LOOKUP(LBNREC,1).EQ.0 . OR.
C        Prog lookups in dump before vn3.2:
     *  (LOOKUP(LBNREC,1).EQ.IMDI. AND. FIXHD(12).LE.301)) THEN
        IF(LEN_DATA.NE.IMDI)THEN
      CALL CHK_LOOK(FIXHD,LOOKUP,LEN1_LOOKUP,LEN_DATA,
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
     &              ICODE,CMESSAGE) 
        ENDIF
      ENDIF

      ENDIF

      RETURN
      END
