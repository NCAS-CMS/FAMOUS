C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL Subroutine INTF_OUT -----------------------------------------------
C
CLL  Purpose: To open boundary files
CLL           Ocean has 4 on Fortran unit numbers 100-103
CLL           Atmos has 8 on Fortran unit numbers 140-147
CLL
CLL  Model            Modification history from model version 4.5
CLL version  Date
CLL  4.5   3/09/98    New deck added M.J.Bell
CLL
CLLEND ---------------------------------------------------------------
       subroutine intf_out (
C*L----------------- COMDECK ADUMLEN --------------------------------
C       Array dimensions
     +  LEN_FIXHD, LEN_INTHD, LEN_REALHD,
     +  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,
     +  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC, LEN2_FLDDEPC,
     +  LEN_EXTCNST,  LEN_DUMPHIST,
     +  LEN_CFI1, LEN_CFI2, LEN_CFI3,
     +  LEN1_LOOKUP, LEN2_LOOKUP,
     &  MPP_LEN1_LOOKUP,
     &  len_ixsts, len_spsts,
C*----------------------------------------------------------------------
C*L----------------- COMDECK AINFLEN --------------------------------
C       Array dimensions
     # N_INTF, INTF_LOOKUPS, PP_LEN_REALHD, PP_LEN_INTHD,
     # MAX_INTF_P_LEVELS, INTF_LEN2_LEVDEPC,
     # TOT_LEN_INTF_P, TOT_LEN_INTF_U, NPTS_U_FIELD,
C*----------------------------------------------------------------------
C===========================COMDECK ARGDUM==========================
     &FIXHD, INTHD, CFI1, CFI2, CFI3, REALHD, LEVDEPC,
     &ROWDEPC, COLDEPC, FLDDEPC, EXTCNST, DUMPHIST,
! PP lookup headers and Ocean stash array + index with lengths
     &LOOKUP,
     &MPP_LOOKUP,
     &ixsts, spsts,
C========================END OF COMDECK ARGDUM======================
CLL----------Headers for  interface data sets ----------------
     & FIXHD_INTF, INTHD_INTF, LOOKUP_INTF,
     & REALHD_INTF,LEVDEPC_INTF,
C*L   Interpolation constants for atmosphere interface data sets.
     &  P_INDEX_B_L, P_INDEX_B_R,  U_INDEX_B_L,  U_INDEX_B_R,
     & P_WEIGHT_T_R, P_WEIGHT_B_L, P_WEIGHT_B_R, P_WEIGHT_T_L,
     & U_WEIGHT_T_R, U_WEIGHT_B_L, U_WEIGHT_B_R, U_WEIGHT_T_L,

C*L   Rotation coefficients for atmosphere interface data sets
     &  COEFF1, COEFF2, COEFF3, COEFF4, COEFF5, COEFF6,

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
     & NFTOUT, JINTF, im, mype,
     & INTF_PACK, INTFWIDTH, LEN_INTF_P, LEN_INTF_U,
     & len_intf_data, item_intf, len_bdy_flds,
     & dump_lookup_intf, intf_data, icode, cmessage )
!---------------------------------------

      implicit none

C*L----------------- COMDECK CDUMLEN --------------------------------
C       Array dimensions
       INTEGER
     +  LEN_FIXHD, LEN_INTHD, LEN_REALHD,
     +  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,
     +  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC, LEN2_FLDDEPC,
     +  LEN_EXTCNST,  LEN_DUMPHIST,
     +  LEN_CFI1, LEN_CFI2, LEN_CFI3,
     +  LEN1_LOOKUP, LEN2_LOOKUP,
     &  MPP_LEN1_LOOKUP,
     &  len_ixsts, len_spsts
C*----------------------------------------------------------------------
C*L----------------- COMDECK CINFLEN --------------------------------
C       Array dimensions
      INTEGER
     # N_INTF, INTF_LOOKUPS, PP_LEN_REALHD, PP_LEN_INTHD,
     # MAX_INTF_P_LEVELS, INTF_LEN2_LEVDEPC,
     # TOT_LEN_INTF_P, TOT_LEN_INTF_U, NPTS_U_FIELD
C*----------------------------------------------------------------------
CL =================== COMDECK TYPDUM ================================
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL                           to be called in the same module.
CL --------------- Dump headers  -----------------
      INTEGER
     &FIXHD(LEN_FIXHD),                       ! fixed length header
     &INTHD(LEN_INTHD),                       ! integer header
     &CFI1(LEN_CFI1+1),                       ! compress field index
     &CFI2(LEN_CFI2+1),                       ! compress field index
     &CFI3(LEN_CFI3+1)                        ! compress field index

      REAL
     &REALHD(LEN_REALHD),                     ! real header
     &LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC+1), ! level  dep const
     &ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC+1), ! row    dep const
     &COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC+1), ! column dep const
     &FLDDEPC(LEN1_FLDDEPC*LEN2_FLDDEPC+1), ! field  dep const
     &EXTCNST(LEN_EXTCNST+1),                 ! extra constants
     &DUMPHIST(LEN_DUMPHIST+1)                  ! temporary hist file

CL --------------- PP headers ---------------------------
      INTEGER
     &LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)         ! IN/OUT: lookup heads
     &, MPP_LOOKUP(MPP_LEN1_LOOKUP,LEN2_LOOKUP)
     &, ixsts(len_ixsts)                      ! stash index array

      INTEGER
     &  spsts(len_spsts)                      ! stash array

CL
CL
CL
C*L---Headers for interface data sets
       INTEGER
     &  FIXHD_INTF(LEN_FIXHD,N_INTF)        ! Fixed header
     & ,INTHD_INTF(PP_LEN_INTHD,N_INTF)     ! Integer header
     & ,LOOKUP_INTF(LEN1_LOOKUP,INTF_LOOKUPS,N_INTF)  ! Lookups

       REAL
     &  REALHD_INTF(PP_LEN_REALHD,N_INTF)   ! Real header

C !!!! check next line !!!

     & ,LEVDEPC_INTF(MAX_INTF_P_LEVELS,INTF_LEN2_LEVDEPC,N_INTF)

C*L---Interpolation constants for interface data sets.
       INTEGER
C                                    Index of corner in source grid box:
     &  P_INDEX_B_L(TOT_LEN_INTF_P),    ! Bottom left  ( p grid)
     &  P_INDEX_B_R(TOT_LEN_INTF_P),    ! Bottom right ( p grid)
     &  U_INDEX_B_L(TOT_LEN_INTF_U),    ! Bottom left  ( u grid)
     &  U_INDEX_B_R(TOT_LEN_INTF_U)     ! Bottom right ( u grid)
       REAL
C                                    Weight applied to value at:
     &  P_WEIGHT_T_R(TOT_LEN_INTF_P),   ! Top    right (p grid)
     &  P_WEIGHT_B_L(TOT_LEN_INTF_P),   ! Bottom left  (p grid)
     &  P_WEIGHT_B_R(TOT_LEN_INTF_P),   ! Bottom right (p grid)
     &  P_WEIGHT_T_L(TOT_LEN_INTF_P),   ! Top    left  (p grid)
     &  U_WEIGHT_T_R(TOT_LEN_INTF_U),   ! Top    right (u grid)
     &  U_WEIGHT_B_L(TOT_LEN_INTF_U),   ! Bottom left  (u grid)
     &  U_WEIGHT_B_R(TOT_LEN_INTF_U),   ! Bottom right (u grid)
     &  U_WEIGHT_T_L(TOT_LEN_INTF_U)    ! Top    left  (u grid)

C*L---Rotation coefficients for atmosphere interface data sets
       REAL
     &  COEFF1(TOT_LEN_INTF_U),
     &  COEFF2(TOT_LEN_INTF_U),
     &  COEFF3(NPTS_U_FIELD),
     &  COEFF4(NPTS_U_FIELD),
     &  COEFF5(TOT_LEN_INTF_U),
     &  COEFF6(TOT_LEN_INTF_U)


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

      integer
     &    NFTOUT,         ! unit to write to
     &    JINTF,          ! number of this boundary file
     &    im,             ! internal model identifier
     &    mype,           ! number of "my" processor
     &    INTF_PACK(N_INTF),       ! Packing Indicator for data
     *    INTFWIDTH(N_INTF),       ! Width of interface zone
     &    LEN_INTF_P(N_INTF),      ! Length of interface p field
     &    LEN_INTF_U(N_INTF),      ! Length of interface u field
     &    len_intf_data,   ! length of field of data to output
     &    item_intf(INTF_LOOKUPS),    ! stash item numbers of fields
     &    len_bdy_flds(INTF_LOOKUPS),  ! length of interface fields
     &    dump_lookup_intf(INTF_LOOKUPS) ! numbers of corresponding
                                          ! dump lookup tables

      real intf_data( len_intf_data )  ! boundary data for output
      integer icode
      character*256 cmessage

!-----------------------------------------------------------
C*L================ COMDECK CMAXSIZE ==========================
C   Description:
C     This COMDECK contains maximum sizes for dimensioning arrays
C   of model constants whose sizes are configuration dependent. This
C   allows constants to be read in from a NAMELIST file and maintain
C   the flexibility of dynamic allocation for primary variables. The
C   maximum sizes should agree with the maximum sizes implicit in the
C   front-end User Interface.
C
CLL
CLL  Model            Modification history:
CLL version  Date
CLL 3.2   26/03/93  New COMDECK. Author R.Rawlins
CLL  3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
CLL                 in macro-tasked calls to SWRAD, LWRAD & CONVECT.
CLL                 Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
CLL  3.5  22/05/95  Add MAX_N_INTF. D. Robinson
CLL  4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.

CLL
C
C
C     MAX_N_INTF/MAX_N_INTF_A to be sorted out in next version
      INTEGER  MAX_N_INTF     ! Max no. of interface areas
        PARAMETER (MAX_N_INTF =  8  )


C*L --------------------- Comdeck: CHSUNITS -------------------------
CLL
CLL Purpose: COMDECK defining the number of i/o units
CLL
CLL  Author : R A Stratton
CLL
CLL  Model            Modification history:
CLL version  date
CLL   3.1  03/02/93   Introduced at version 3.1
CLL   4.1  21/02/96   Increase no.of i/o units to accommodate wave 
CLL                   sub-model.  RTHBarnes.
CLL
CLL Project task:
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLL
CLLEND---------------------------------------------------------------
C
C*L Type declarations
C
      INTEGER NUNITS          ! No. of I/O units
      INTEGER NUNITS_LEN      ! length of most unit no arrays
!
!     These values must be consistent with OUTFILE_S, OUTFILE_L
!     and OUTFILE_E in comdeck VERSION.
      PARAMETER(NUNITS=149)
      PARAMETER(NUNITS_LEN=NUNITS-19)
C
C   The above parameter statements must not be altered without 
C   considering the effect on the following HISTORY COMDECKs 
C    CHISTO, CLFHIST and IHISTO.    
C  This comdeck must always preceed the above history file comdecks.
C   New file environment variable names may need to be added to     
C    CLFHIST and/or CENVIRDT (usually both) depending on manner of I/O.
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     & MODEL_DATA_TIME(6),     ! Array containing model data time
!                              ! (Same as MODEL_BASIS_TIME/MODEL
!                              !  ANALYSIS_HRS depending whether
!                              !  before/after assimilation)
     & RUN_MEANCTL_RESTART,    ! Indicator for next mean period
!                              ! to be processed
     & RUN_INDIC_OP            ! Indicator of operational run type
C
      INTEGER
     & RUN_RESUBMIT_TARGET(6), ! Final target date for the run
C
     & FT_LASTFIELD(20:NUNITS) ! Last field written/read per FT unit
C
C
C History Common Block for overall model integers variables.
C
      COMMON /IHISTO/
     & MODEL_DATA_TIME,
     & RUN_MEANCTL_RESTART, RUN_INDIC_OP,
     & RUN_RESUBMIT_TARGET, FT_LASTFIELD
C
      NAMELIST /NLIHISTO/
     & MODEL_DATA_TIME,
     & RUN_MEANCTL_RESTART, RUN_INDIC_OP,
     & RUN_RESUBMIT_TARGET, FT_LASTFIELD
C
C*L --------------------- Comdeck: CLFHIST  ----------------------------
CLL
CLL  Purpose: COMDECK defining unit numbers relevant to history file
CLL           and variables used to hold the logical to physical
CLL           file associations made within the model
CLL
CLL  Author : A. Sangster
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLL                  Version 5  18/6/90
CLL
CLL  Model             Modification history from model version 3.0
CLL version  Date
CLL
CLL  3.4  30/09/94  Add files MURKFILE,OUSRANCL,OUSRMULT at 109,113,114
CLL  3.4  05/09/94  Add files USRANCIL,USRMULTI at unit nos. 111,112.
CLL
CLL  3.3  22/11/93  Add file SOURCES at unit number 110. R.T.H.Barnes.
CLL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
CLL  Vn3.0  12/02/93 - Variables PERTURB and TRANSP equivalenced to unit
CLL                    numbers 37, and 97 respectively. C.S. Douglas
CLL  3.4  1/8/94     Revised Obs file specification: Stuart Bell
CLL  3.5  01/05/95  Sub-models stage 1: History/control files. RTHBarnes
!    4.0  22/09/95  Added units for Spectral data for Radiation scheme.
!                                        (J. M. Edwards)
CLL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1  26/02/96  Associate new env. variables SO2NATEM and CHEMOXID
!                   with unit nos. 115 & 116. Rename SOURCES to 
!                   SULPEMIS. D. Robinson.
!  4.3   18/3/97  Add aerosol forcings of climate change.  Will Ingram
!  4.4   4/7/97   Add ANLINCR  Chris Jones/Stuart Bell
CLL  4.4   12/9/97  Associate ancillary file EVs for initial surface 
CLL                 type fracs, initial vegetation state and vegetation 
CLL                 disturbance with unit no.s 135-137 R. Betts
CLL  4.4  17/10/97  Associate env var. CACHED with Unit 138. D Robinson 
CLL  4.5  22/04/98  Add new ancillary file for soot emissions:  
CLL                 SOOTEMIS - in I/O unit 139. R.Rawlins
CLL  4.5  29/07/98  Add new variables ALABCOU5/6/7/8. D. Robinson.
CLL  4.5  17/08/98  Add new variables OLABCOU1/2/3/4. Remove
CLL                 OLABCOUT. D. Robinson.
CLL
CLL  Type declarations
CLL
CLL
CLL  Logical Filenames used in the model
CLL
      CHARACTER*80 HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,
     *             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,
     *             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,
     2             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,
     *             AOPSUM4,AOMEAN,SSU,
     3             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,
     *             SICEIN,PERTURB,MASK,
     4             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3,
     *             AOPSTMP4,
     5             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,
     &             SWSPECTD,
     6             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,
     &             OBS01,OBS02,OBS03,OBS04,OBS05,
     &             OBS06,OBS07,OBS08,OBS09,OBS10,
     8             LWSPECTD,WAVEOUT,SURGEOUT,MESOUT,STRATOUT,WFOUT,     
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT,
     *             CURNTOUT,ALABCIN,OROG,OLABCIN,OCNDEPTH,CURNTIN,
     *             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND
     &             ,TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,
     *             OUSRANCL,OUSRMULT,MURKFILE,
     *             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4
     &            ,ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8
     &            ,OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4
     &            ,ANLINCR
     &            ,WINITIAL,WSTART,WRESTART,WAVANL,WAVANCIN  
     &            ,SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB
     &            ,CACHED,SOOTEMIS
     &            ,CO2EMITS
C
      CHARACTER*80 MODEL_FT_UNIT ! Array holding FORTRAN unit file
C                                ! associations details for each unit
C
      INTEGER
     *        MCTL_UNIT,         ! Master control namelist file unit
     *        ICTL_UNIT,         ! Interim control namelist file unit
     *        PHIST_UNIT,        ! Permanent history file unit
     *        IHIST_UNIT,        ! Interim history file unit
     *        THIST_UNIT,        ! Temporary history file unit
     *        FTXX_UNIT,         ! Logical/physical file associations
     *        HKFILE_UNIT        ! Operational houskeeping file unit
C*
C  Parameters specifying unit numbers relevant to control/history tasks
C
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(PHIST_UNIT =10)
      PARAMETER(IHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)
      PARAMETER(FTXX_UNIT  =13)
!
! Namelist of all permissible logical files.
!
      NAMELIST / NLCFILES /
     &             HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,
     &             AOPSUM4,AOMEAN,SSU,
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,
     &             SICEIN,PERTURB,MASK,
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3,
     &             AOPSTMP4,
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,
     &             SWSPECTD,
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,
     &             OBS01,OBS02,OBS03,OBS04,OBS05,
     &             OBS06,OBS07,OBS08,OBS09,OBS10,
     &             LWSPECTD,WAVEOUT,SURGEOUT,MESOUT,STRATOUT,WFOUT,     
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT,
     &             CURNTOUT,ALABCIN,OROG,OLABCIN,OCNDEPTH,CURNTIN,
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND
     &             ,TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,
     &             OUSRANCL,OUSRMULT,MURKFILE,
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4
     &            ,ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8
     &            ,OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4
     &            ,ANLINCR
     &            ,WINITIAL,WSTART,WRESTART,WAVANL,WAVANCIN 
     &            ,SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB
     &            ,CACHED,SOOTEMIS
     &            ,CO2EMITS
C
C Common block definition
C
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
C
C  Equivalence logical filenames within array MODEL_FT_UNIT
C
      EQUIVALENCE
     *(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ),
     *(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ),
     *(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ),
     *(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ),
     *(ICTL      ,MODEL_FT_UNIT(9)  ),(PHIST      ,MODEL_FT_UNIT(10) ),
     *(IHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ),
     *(FTXX      ,MODEL_FT_UNIT(13) ),
     *(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ),
     *(AOTRANS   ,MODEL_FT_UNIT(17) ),(ASWAP      ,MODEL_FT_UNIT(18) ),
     *(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ),
     *(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ),
     *(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ),
     *(AOPSUM3   ,MODEL_FT_UNIT(25) )
C
      EQUIVALENCE
     *(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ),
     *(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ),
     *(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ),
     *(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ),
     *(VEGTYPE   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ),
     *(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ),
     *(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ),
     *(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ),
     *(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ),
     *(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ),
     *(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ),
     *(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ),
     *(WFIN      ,MODEL_FT_UNIT(50) )
C
      EQUIVALENCE
     *(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ),
     *(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ),
     *                                (FLUXCORR   ,MODEL_FT_UNIT(56) ),
     *(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), 
     *(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ),
     *(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ),
     *(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ),
     *(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ),
     *(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ),
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ),
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ),
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) )
C
      EQUIVALENCE
     &(OBS06     ,MODEL_FT_UNIT(75) ),(OBS07      ,MODEL_FT_UNIT(76) ),
     &(OBS08     ,MODEL_FT_UNIT(77) ),(OBS09      ,MODEL_FT_UNIT(78) ),
     &(OBS10     ,MODEL_FT_UNIT(79) ),(LWSPECTD   ,MODEL_FT_UNIT(80) ), 
     *(WAVEOUT   ,MODEL_FT_UNIT(81) ),(SURGEOUT   ,MODEL_FT_UNIT(82) ),
     *(MESOUT    ,MODEL_FT_UNIT(83) ),(STRATOUT   ,MODEL_FT_UNIT(84) ),
     *(WFOUT     ,MODEL_FT_UNIT(85) ),(HFLUXOUT   ,MODEL_FT_UNIT(86) ),
     *(PMEOUT    ,MODEL_FT_UNIT(87) ),(ICEFOUT    ,MODEL_FT_UNIT(88) ),
     &(MOSOUT    ,MODEL_FT_UNIT(89) ),
     *(SSTOUT    ,MODEL_FT_UNIT(91) ),(SICEOUT    ,MODEL_FT_UNIT(92) ),
     *(CURNTOUT  ,MODEL_FT_UNIT(93) ),(FLXCROUT   ,MODEL_FT_UNIT(94) ),
     *(ALABCIN   ,MODEL_FT_UNIT(95) ),(OROG       ,MODEL_FT_UNIT(96) ),
     *(TRANSP    ,MODEL_FT_UNIT(97) ),(OLABCIN    ,MODEL_FT_UNIT(98) ),
     *(OCNDEPTH  ,MODEL_FT_UNIT(99) ),
     &(OLABCOU1  ,MODEL_FT_UNIT(100)),(OLABCOU2   ,MODEL_FT_UNIT(101)),
     &(OLABCOU3  ,MODEL_FT_UNIT(102)),(OLABCOU4   ,MODEL_FT_UNIT(103)),
     &(ANLINCR   ,MODEL_FT_UNIT(108)),(MURKFILE   ,MODEL_FT_UNIT(109)),
     &(SULPEMIS  ,MODEL_FT_UNIT(110)),(USRANCIL   ,MODEL_FT_UNIT(111)),
     *(USRMULTI  ,MODEL_FT_UNIT(112)),(OUSRANCL   ,MODEL_FT_UNIT(113)),
     *(OUSRMULT  ,MODEL_FT_UNIT(114)),(SO2NATEM   ,MODEL_FT_UNIT(115)),
     &(CHEMOXID  ,MODEL_FT_UNIT(116)),(AEROFCG    ,MODEL_FT_UNIT(117)),
     *(CO2EMITS  ,MODEL_FT_UNIT(118)),
     *(WINITIAL  ,MODEL_FT_UNIT(130)),(WSTART     ,MODEL_FT_UNIT(131)), 
     *(WRESTART  ,MODEL_FT_UNIT(132)),(WAVANL     ,MODEL_FT_UNIT(133)), 
     *(WAVANCIN  ,MODEL_FT_UNIT(134)),(FRACINIT   ,MODEL_FT_UNIT(135)),
     *(VEGINIT   ,MODEL_FT_UNIT(136)),(DISTURB    ,MODEL_FT_UNIT(137)),
     &(CACHED    ,MODEL_FT_UNIT(138)),(SOOTEMIS   ,MODEL_FT_UNIT(139)),
     &(ALABCOU1  ,MODEL_FT_UNIT(140)),(ALABCOU2   ,MODEL_FT_UNIT(141)),
     &(ALABCOU3  ,MODEL_FT_UNIT(142)),(ALABCOU4   ,MODEL_FT_UNIT(143)),
     &(ALABCOU5  ,MODEL_FT_UNIT(144)),(ALABCOU6   ,MODEL_FT_UNIT(145)),
     &(ALABCOU7  ,MODEL_FT_UNIT(146)),(ALABCOU8   ,MODEL_FT_UNIT(147)) 
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
! ----------------------- Comdeck: CNTLALL  ----------------------------
! Description: COMDECK defining Control variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  25/10/95  Add user switch CONTROL_RESUBMIT. RTHBarnes
!  4.4  28/07/97  Add user switch LCLIMREALYR. M Gallani
!  4.4  11/10/97  Add logical switch L_AO_D1_MEMORY. D. Robinson. 
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     &        MODEL_BASIS_TIME(6),! Array holding original data time
!                                 ! (prior to assimilation)
     &        MODEL_ANALYSIS_HRS, ! Model analysis time in hours since
!                                 ! Basis Time
     &        MODEL_HRS_PER_GROUP,! No. of hours in coupling period
     &        NCPU,               ! No of CPUs assigned to the program
     &        ANCIL_REFTIME(6),   ! Ref. time for updating ancillaries
     &        FT_PLOTSEL(60:69),  ! interval for plotting pp file
     &        RUN_TARGET_END(6),   ! Target end time for this run
     &        RUN_RESUBMIT_INC(6), ! Increment to be added on each
!                                  ! resubmission of the job.
     &   PP_LEN2_LOOK(20:NUNITS),  ! Number of field headers reserved
!                                  !  for non-mean PPfiles on each unit
     &   PP_PACK_CODE(20:NUNITS),  ! Internally defined PP packing code
     &   FT_STEPS(20:NUNITS),   ! Frequency of initialisation of FTunit
     &   FT_FIRSTSTEP(20:NUNITS)   ! ... starting at step number ..
!
      LOGICAL
     &       LATMOSNEXT,LOCEANNEXT,  ! Flags to select atmosphere/ocean
     &       LPP,                    ! Activate PPCTL
     &       LPP_SELECT(20:NUNITS),  ! Activate PP init'sation on unit
     &       LDUMP,                  ! Activate DUMPCTL
     &       LMEAN,                  ! Activate MEANCTL
     &       LHISTORY,               ! Update TEMP history file
     &       LPRINT,                 ! Activate PRINTCTL
     &       LINTERFACE,             ! Activate GEN_INTF
     &       LEXIT,                  ! Activate EXITCHEK
     &       LJOBRELEASE,            ! Activate JOBCTL
     &       LMEANPR(4),             ! Select printed diags from means
     &       LANCILLARY,             ! Activate UP_ANCIL
     &       LBOUNDARY,              ! Activate UP_BOUND
     &       LASSIMILATION,          ! Activate assimilation
     &       LCAL360,                ! 360-day calendar
     &       LTIMER                  ! Activate detailed TIMER routine
     &      ,L_AO_D1_MEMORY  ! T : D1 copied to memory for AO coupling
     &      ,LCLIMREALYR             ! Real-period climate means

      CHARACTER*4  EXPT_ID          ! Unique alphanumeric serial number
!                                   ! associated with model
!                                   ! (Non-Operational expts)
!                                   !
!                                   ! Operational run name
!                                   ! (Operational expts)
      CHARACTER*8  EXPT_ALIAS       ! Non unique user defined expt name
      CHARACTER*1  JOB_ID           ! Unique alphanumeric job identifier
!                                   ! used for networking
      CHARACTER*4  EXPT_ID_IN       ! Experiment ID of driving model if
!                                   ! limited-area run
      CHARACTER*4  JOB_ID_IN        ! Job ID of driving model if
!                                   ! limited-area run
      CHARACTER*14 MODEL_STATUS     ! Operational or NonOperational
      CHARACTER*14 MODEL_ASSIM_MODE ! Atmosphere,Ocean,Coupled or None
      CHARACTER*17 TIME_CONVENTION  ! Relative, Timestep, Absolute_long,
!                                    Absolute_standard or Absolute_short
      CHARACTER*1  FT_WSSEND(60:69) ! "Y" if file to be sent to HP
!
      CHARACTER*1 TYPE_LETTER_1(20:NUNITS) ! File type letter #1
      CHARACTER*1 TYPE_LETTER_2(20:NUNITS) ! File type letter #2
      CHARACTER*1 TYPE_LETTER_3(20:NUNITS) ! File type letter #3
!
      CHARACTER*1  FT_INPUT (20:NUNITS) ! "Y" if input file on unit
      CHARACTER*1  FT_OUTPUT(20:NUNITS) ! "Y" if output file on unit
      CHARACTER*1  FT_SELECT(20:NUNITS) ! "Y" if file selected for post
!                                          processing request.
      CHARACTER*1  FT_ARCHSEL(20:NUNITS) ! "Y" if file to be archived.
!
      CHARACTER*10 RUN_ASSIM_MODE      ! cf MODEL_ASSIM_MODE (Oper use)
      CHARACTER*1  CONTROL_RESUBMIT    ! User flag for auto resubmit

      NAMELIST / NLSTCALL /
     & MODEL_BASIS_TIME, MODEL_ANALYSIS_HRS,
     & MODEL_HRS_PER_GROUP,
     & NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,
     & RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,
     & FT_STEPS, FT_FIRSTSTEP,
     & LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,
     & LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,
     & LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,
     & LCAL360, LTIMER, L_AO_D1_MEMORY,
     & LCLIMREALYR,
     & EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,
     & EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,
     & TIME_CONVENTION, FT_WSSEND,
     & TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,
     & FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,
     & RUN_ASSIM_MODE, CONTROL_RESUBMIT

      COMMON / CNTLCALL /
     & MODEL_BASIS_TIME, MODEL_ANALYSIS_HRS,
     & MODEL_HRS_PER_GROUP,
     & NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,
     & RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,
     & FT_STEPS, FT_FIRSTSTEP,
     & LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,
     & LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,
     & LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,
     & LCAL360, LTIMER, L_AO_D1_MEMORY,
     & LCLIMREALYR,
     & EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,
     & EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,
     & TIME_CONVENTION, FT_WSSEND,
     & TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,
     & FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,
     & RUN_ASSIM_MODE, CONTROL_RESUBMIT
! ----------------------- Comdeck: CNTLGEN  ----------------------------
! Description: COMDECK defining Control variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0   3/11/95  Move character array MEANSim to the end of the 
!                 common block to ensure that it starts correctly on a
!                 word boundary. [No problem is apparent on the Cray  
!                 if N_INTERNAL_MODEL_MAX is an even no.]
!                 Rick Rawlins                                    
!  4.1  03/04/96  Add new array DUMP_PACKim. D. Robinson
!  4.5  10/11/98  Increase number of dumps allowed at irregular 
!                 timesteps from 10 to 40: Move lengths into
!                 CNTLGEN. R Rawlins
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     & DUMPTIMES_LEN1   ! Max no. of irregular times for dumps
     &,PRINTFREQ_LEN1   ! No. of areas of zonal mean prints
     &,MEANFREQ_LEN1    ! No. of time intervals for climate meaning
     &,JOBREL_LEN1      ! Max no. of irregular times for job release

      PARAMETER(
     & DUMPTIMES_LEN1 = 40 
     &,PRINTFREQ_LEN1 = 5
     &,MEANFREQ_LEN1  = 4
     &,JOBREL_LEN1    = 10
     &) 
      INTEGER
     & STEPS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
     &,SECS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
     &,EXITFREQim(N_INTERNAL_MODEL_MAX)  ! Number of advection
!                               timesteps between checks for model exit
     &,DUMPFREQim(N_INTERNAL_MODEL_MAX)  ! Number of steps between
!                                              atmosphere restart dumps
     &,ARCHDUMP_FREQim(N_INTERNAL_MODEL_MAX)  ! Archiving frequency
!                                                   for atmos dumps
     &,DUMPTIMESim(DUMPTIMES_LEN1,N_INTERNAL_MODEL_MAX) ! Timesteps 
!            (from start of run) at which restart dumps are written
     &,MEANFREQim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Indicators 
!             for mean dump frequency
     &,MEANARCHim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Switches 
!             for mean dump arch.
     &,PPSELECTim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !PP field 
!             selectors 
     &,ARCHPPSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)!Switches 
!             for pp field archive
     &,PLOTSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)  !Switches 
!             for chart plotting
     &,PP_LEN2_MEANim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Number of
!             field headers to reserve for internal model mean PPfiles
     &,MEAN_REFTIMEim(6,N_INTERNAL_MODEL_MAX)  ! Reference time for
!                                                production of means
     &,PRINTFREQim(PRINTFREQ_LEN1,N_INTERNAL_MODEL_MAX) ! Indicators   
!             of zonal mean print frequency
     &,JOBREL_STEPim(JOBREL_LEN1,N_INTERNAL_MODEL_MAX)  ! Step numbers 
!             at which to release user-specified scripts 
     &,ARCHDUMP_OFFSETim(N_INTERNAL_MODEL_MAX)!Offset for dump archiving
     &,FT_MEANim(N_INTERNAL_MODEL_MAX)     ! Unit reserved for mean PPs
     &,DUMP_PACKim(N_INTERNAL_MODEL_MAX)  ! Packing indicator for dumps
      CHARACTER*1  MEANWSim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) 
!                                    "Y" if mean file to be sent to HP
      LOGICAL    LLBOUTim(N_INTERNAL_MODEL_MAX)  ! Lateral b.c.'s
     &          ,LANCILim(N_INTERNAL_MODEL_MAX)  ! Ancillary files
C
      NAMELIST / NLSTCGEN /
     & STEPS_PER_PERIODim, SECS_PER_PERIODim,
     & EXITFREQim, DUMPFREQim,
     & ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,
     & ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,
     & PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim,
     & FT_MEANim,
     & DUMP_PACKim,
     & MEANWSim, LLBOUTim, LANCILim
C
      COMMON / CNTLCGEN /
     & STEPS_PER_PERIODim, SECS_PER_PERIODim,
     & EXITFREQim, DUMPFREQim,
     & ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,
     & ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,
     & PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim,
     & FT_MEANim,
     & DUMP_PACKim,
     &  LLBOUTim, LANCILim,
     &  MEANWSim


CLL  Comdeck: CTIME ----------------------------------------------------
CLL
CLL  Purpose: Derived model time/step information including start/end
CLL           step numbers and frequencies (in steps) of interface field
CLL           generation, boundary field updating, ancillary field
CLL           updating; and assimilation start/end times.
CLL           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
CLL           Also contains current time/date information, current
CLL           step number (echoed in history file) and steps-per-group.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.1   13/02/93  Dimension arrays A_INTERFACE_STEPS/FSTEP/LSTEP
CLL                   D. Robinson
CLL   3.3  01/02/94  Add BASIS_TIME_DAYS to BASIS_TIME_SECS for revised
CLL                  (32-bit portable) model clock calculations. TCJ
CLL  3.4  13/12/94  Change COMMOM name from CTIME to CTIMED to satisfy
CLL                 DEC alpha compiler for portability.  N.Farnon.
CLL  3.5  12/04/95  Stage 1 submodel changes: move to dimensioning
CLL                 arrays by internal model. R.Rawlins
CLL  4.4  06/10/97  Data time of IAU dump added. Adam Clayton.
CLL  4.5  21/08/98  Remove redundant code. D. Robinson.
CLL
CLL Programming standard :
CLL
CLL  Logical components covered: C0
CLL
CLL Project task :
CLL
CLL External documentation: Unified Model documentation paper No:
CLL                         Version:
CLL
CLLEND -----------------------------------------------------------------
C
      INTEGER
     1     I_YEAR,                 ! Current model time (years)
     2     I_MONTH,                ! Current model time (months)
     3     I_DAY,                  ! Current model time (days)
     4     I_HOUR,                 ! Current model time (hours)
     5     I_MINUTE,               ! Current model time (minutes)
     6     I_SECOND,               ! Current model time (seconds)
     7     I_DAY_NUMBER,           ! Current model time (day no)
     8     PREVIOUS_TIME(7),       ! Model time at previous step
     9     DATA_MINUS_BASIS_HRS,   ! Data time - basis time (hours)
     A     IAU_DATA_TIME(6)        ! Data time of IAU dump.
      INTEGER
     &       BASIS_TIME_DAYS,     ! Integral no of days to basis time
     3       BASIS_TIME_SECS,     ! No of seconds-in-day at basis time
     4       FORECAST_HRS         ! Hours since Data Time (ie T+nn)
      INTEGER
     H       O_CLM_FIRSTSTEP,     ! First } step for ocean climate
     I       O_CLM_LASTSTEP       ! Last  } increments
C
      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,
     1               I_DAY_NUMBER,PREVIOUS_TIME,
     &               BASIS_TIME_DAYS,BASIS_TIME_SECS,
     &               FORECAST_HRS,DATA_MINUS_BASIS_HRS,
     &               IAU_DATA_TIME,
     C               O_CLM_FIRSTSTEP,   O_CLM_LASTSTEP

      INTEGER
     * STEPim(INTERNAL_ID_MAX)            ! Step no since basis time
     *,GROUPim(INTERNAL_ID_MAX)           ! Number of steps per group
     *,TARGET_END_STEPim(INTERNAL_ID_MAX) ! Finish step number this run

      REAL
     & SECS_PER_STEPim(INTERNAL_ID_MAX)   ! Timestep length in secs

      INTEGER
     * INTERFACE_STEPSim(MAX_N_INTF,INTERNAL_ID_MAX)     ! Frequency of
!                              ! interface field generation in steps
     *,INTERFACE_FSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)     ! Start steps
!                              ! for interface field generation
     *,INTERFACE_LSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)     ! End   steps
!                              ! for interface field generation
     *,BOUNDARY_STEPSim(INTERNAL_ID_MAX)                 ! Frequency of
!                              ! updating boundary fields in steps
     *,BNDARY_OFFSETim(INTERNAL_ID_MAX)!  No of steps from boundary data
!                              ! prior to basis time to model basis time
     *,ANCILLARY_STEPSim(INTERNAL_ID_MAX) ! Lowest frequency for
!                              ! updating of ancillary fields in steps
     *,ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX) ! Start steps for assimilation
     *,ASSIM_STEPSim(INTERNAL_ID_MAX)     ! Number of assimilation
!                              ! steps to analysis
     *,ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)! Number of assimilation
!                              ! steps after analysis
      COMMON/CTIMEE/
     & STEPim,GROUPim,TARGET_END_STEPim
     &,INTERFACE_STEPSim
     &,INTERFACE_FSTEPim
     &,INTERFACE_LSTEPim
     &,BOUNDARY_STEPSim
     &,BNDARY_OFFSETim
     &,ANCILLARY_STEPSim
     &,ASSIM_FIRSTSTEPim
     &,ASSIM_STEPSim
     &,ASSIM_EXTRASTEPSim
     &,SECS_PER_STEPim
!
!====================== COMDECK CNTL_IO ========================
! Description:
!
!     Defines the sector size for well-formed transfers on Cray
!     Research systems.  Disk addresses must start on a sector
!     boundary, and transfers must be a number of sectors.  Disk
!     word addresses start at 0.
!
!     On the T3E, well-formed transfers must also start on a
!     cache-line boundary in memory.
!
!   4.3    30/04/97  New deck       B. Carruthers, Cray Research
!   4.4    27/10/97  Remove DATA statement. C.P. Jones
!
C
      INTEGER UM_SECTOR_SIZE    ! Sector size on disk for I/O
C
      COMMON / CNTL_IO / UM_SECTOR_SIZE

CL Local variables

      integer
     &  LEN_PPNAME,  ! length of pp file name
     &  NTIME,       ! number of field in output file
     &  lookup_start,! start location to write lookup table
     &  LEN_IO,
     &  disk_address,
     *  start_addr,
     &  len_data,       !
     &  var,            ! loop index for variable
     &  i,              ! loop index
     &  N1,             ! local packing index
     &  disk_length,
     &  iaddr,
     &  j,              ! loop index
     &  data_start,     ! start location for writing interface data

     &  EXPPXI,     ! function
     &  P21BITS     ! function
       EXTERNAL EXPPXI, P21BITS

      integer SEC,YY,MM,DD,HR,MN,SS,DAY_NO

      real a_io

      LOGICAL LPACK_32B,    ! pack as 32 bit numbers
     &        LPACK_PPXREF  !

      CHARACTER*80 STRING         ! work array
      CHARACTER*14 PPNAME         ! boundary output filename
!-------------------------------------------------------------

CL 0. Miscellaneous Preliminaries

      LPACK_32B = INTF_PACK(JINTF).EQ.1
      LPACK_PPXREF = INTF_PACK(JINTF).EQ.2

CL 1.0   Open file; determine where to write new data

CL     Open boundary output file if reinitialised during run

      IF (FT_STEPS(NFTOUT).GT.0) THEN
        STRING = MODEL_FT_UNIT(NFTOUT)
        PPNAME = STRING(18:31)
        LEN_PPNAME = LEN(PPNAME)
        CALL FILE_OPEN(NFTOUT,PPNAME,LEN_PPNAME,1,1,ICODE)
        IF (ICODE.NE.0) THEN
          CMESSAGE="INTF_OUT: Error opening preassigned boundary file"
          GO TO 999   !  Return
        ENDIF
      ENDIF

C      Determine position where to Buffer out data to
      NTIME=FT_LASTFIELD(NFTOUT)+1

CL 2.  Set up headers

CL 2.1 Fixed length header
      FIXHD_INTF(152,JINTF) = INTF_LOOKUPS*NTIME
      FIXHD_INTF(161,JINTF) = LEN_INTF_DATA*NTIME

CL 2.2 Integer Constants
      INTHD_INTF(3,JINTF) = NTIME

CL 2.3 LOOKUP Table

C  2.3.1   Determine position in LOOKUP table
      LOOKUP_START=FIXHD_INTF(150,JINTF) +
     &             FIXHD_INTF(151,JINTF)*INTF_LOOKUPS*(NTIME-1) - 1

C 2.3.2  For well-formed I/O re-read the last lookup
C       table on disk to find disk_address
C       also set initial start address

      if(ntime.ne.1) then
        call setpos(nftout, lookup_start-len1_lookup, icode)
        call buffin(nftout, lookup_intf(1, 1, jintf), len1_lookup,
     &   len_io, a_io)

c--check for errors
        if(a_io.ne.-1.0 .or. len_io.ne.len1_lookup) then
          call ioerror('intf_out: Buffer in of Last Lookup Header',
     &     a_io, len_io, len1_lookup)
          cmessage=' intf_out: I/O Error on reading last lookup'
          icode=5
          goto 999
        endif

c--compute the new disk address from the last address and length
        disk_address=lookup_intf(lbegin, 1, jintf)+
     &               lookup_intf(lbnrec, 1, jintf)

      else     ! ntime

        disk_address=fixhd_intf(160, jintf)-1
      endif     ! ntime

c--round this disk address to ensure we start on a sector boundary
      disk_address=((disk_address+um_sector_size-1)/
     & um_sector_size)*um_sector_size

C - start address (not used by well formed I/O ?)
      START_ADDR = FIXHD_INTF(161,JINTF)-LEN_INTF_DATA+1

C 2.3.3  Check that there is enough space for this entry in LOOKUP table

      IF (FIXHD_INTF(150,JINTF)+
     &    FIXHD_INTF(151,JINTF)*FIXHD_INTF(152,JINTF).GT.
     &   FIXHD_INTF(160,JINTF)) THEN
        CMESSAGE=' INTF_OUT: Insufficient space for headers in boundary
     &                       dataset.'
        ICODE=1
        GO TO 999   !  Return
      ENDIF

C 2.3.5 Set validity times

      SEC = STEPim(im) * SECS_PER_PERIODim(im) /
     &      STEPS_PER_PERIODim(im)

      CALL SEC2TIME(0,SEC,BASIS_TIME_DAYS,BASIS_TIME_SECS,
     &                  YY,MM,DD,HR,MN,SS,DAY_NO,LCAL360)

      DO VAR = 1,  INTF_LOOKUPS

C 2.3.6 Initialise lookup tables (with values from dump lookup tables)

        DO I=1,LEN1_LOOKUP
          LOOKUP_INTF(I,VAR,JINTF)=LOOKUP(I,dump_lookup_intf(var))
        ENDDO

C 2.3.7 Set times in lookup tables

        LOOKUP_INTF(LBYR ,VAR,JINTF) = YY
        LOOKUP_INTF(LBMON,VAR,JINTF) = MM
        LOOKUP_INTF(LBDAT,VAR,JINTF) = DD
        LOOKUP_INTF(LBHR ,VAR,JINTF) = HR
        LOOKUP_INTF(LBMIN,VAR,JINTF) = MN
        LOOKUP_INTF(LBDAY,VAR,JINTF) = DAY_NO

        LOOKUP_INTF(LBYRD ,VAR,JINTF) = FIXHD(21)
        LOOKUP_INTF(LBMOND,VAR,JINTF) = FIXHD(22)
        LOOKUP_INTF(LBDATD,VAR,JINTF) = FIXHD(23)
        LOOKUP_INTF(LBHRD ,VAR,JINTF) = FIXHD(24)
        LOOKUP_INTF(LBMIND,VAR,JINTF) = FIXHD(25)
        LOOKUP_INTF(LBDAYD,VAR,JINTF) = FIXHD(27)

C  2.3.8 Set the length of the field in LOOKUP table
C (simpler than in original atmosphere code)  !! CHECK THIS !!

        LOOKUP_INTF(LBLREC,VAR,JINTF) = len_bdy_flds(var)

C 2.3.9 Set packing info
        N1 = 0   !  Data not packed
        IF (LPACK_32B) N1 = 2  ! Data packed as 32 bits
        IF (LPACK_PPXREF) THEN
          N1 = EXPPXI(im,0,item_intf,ppx_dump_packing,
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
     &                 icode,cmessage)
          if (icode .gt. 0) then
             write(6,*) 'exppxi failed in intf_out'
             go to 999
          end if
        END IF
        LOOKUP_INTF(LBPACK,VAR,JINTF)= N1

C 2.3.10 Store the disk address; and calculate for next field
        lookup_intf(lbegin, var, jintf)=disk_address

c--fetch the data field length, allowing for packing
        if(mod(lookup_intf(lbpack, var, jintf), 10).eq.2) then
          disk_length=(lookup_intf(lblrec, var, jintf)+1)/2
        else
          disk_length=lookup_intf(lblrec, var, jintf)
        endif

c--store the rounded-up length
C NB !! This length is not checked to fit sectors !!
        lookup_intf(lbnrec, var, jintf)=disk_length

c--update the disk address
        disk_address=disk_address+disk_length

C 2.3.11 Set other elements in the lookup table

C grid code ; should be 1 for all variables (correction at 4.4)
C !!! ideally should be 101 for rotated grid
        LOOKUP_INTF(LBCODE,VAR,JINTF)=1

        LOOKUP_INTF(LBHEM,VAR,JINTF)=99
        LOOKUP_INTF(LBROW,VAR,JINTF)=INTFWIDTH(JINTF)

C numbers of rows & columns;  var=2 or 3 is not suitable for ocean
        LOOKUP_INTF(LBNPT,VAR,JINTF) =
     &           LEN_INTF_P(JINTF)/INTFWIDTH(JINTF)
        IF (  ( IM .EQ. 1. .AND. (VAR.EQ.2.OR.VAR.EQ.3) )
     &   .OR. ( IM .EQ. 2. .AND. (VAR.EQ.3.OR.VAR.EQ.4) ) ) THEN
          LOOKUP_INTF(LBNPT,VAR,JINTF) =
     &           LEN_INTF_U(JINTF)/INTFWIDTH(JINTF)
        END IF

        LOOKUP_INTF(LBLEV,VAR,JINTF)=-1
        LOOKUP_INTF(NADDR,VAR,JINTF) = START_ADDR
        START_ADDR = START_ADDR + LOOKUP_INTF(LBLREC,VAR,JINTF)

      END DO  ! VAR

CL 3. Pack data as required

        IADDR = 1
        LEN_DATA = 0

        DO VAR = 1,INTF_LOOKUPS
          IF (MOD(LOOKUP_INTF(LBPACK,VAR,JINTF),10).EQ.2) THEN
CL 3.1 Pack this data field

            IF (mype .EQ. 0) THEN
            CALL PACK21(LOOKUP_INTF(LBLREC,VAR,JINTF),
     &                  INTF_DATA(IADDR),INTF_DATA(LEN_DATA+1),
     &                  P21BITS(FIXHD_INTF(12,JINTF)))
            ENDIF

c--the (+1) in the expression below is unnecessary, since
c  LBC data is composed of two rows NS and two rows EW, and
c  thus always has an even number of data points.  If this
c  is not true, then READFLDS will either get the data one
c  out downwards if the (+1) is omitted, or one word upwards
c  if the (+1) is added.  In other words, the packing will
c  cause either one word to be omitted or one word added in
c  the data after the read.  This is because READFLDS reads
c  and converts the whole LBC record at one go, rather than
c  as a series of separate records.
            LEN_DATA = LEN_DATA+(LOOKUP_INTF(LBLREC,VAR,JINTF)+1)/2

c--check that we are not packing an odd nuber of words
            if((lookup_intf(lblrec,var,jintf)/2)*2 .ne.
     &       lookup_intf(lblrec,var,jintf)) then
              write(6,7734) lookup_intf(lblrec,var,jintf)
7734          format(/'LBC Data contains ',i10,' Words, which is',
     &         ' an Odd Number which is not allowed for 32-bit',
     &         ' Packing')
            endif

          ELSE        !    LOOKUP_INTF(LBPACK..

CL 3.2 Copy unpacked data to new location if necessary
            IF (LEN_DATA+1.LT.IADDR) THEN
            IF (mype .EQ. 0) THEN
              DO J = 1,LOOKUP_INTF(LBLREC,VAR,JINTF)
                INTF_DATA(LEN_DATA+J) = INTF_DATA(IADDR+J-1)
              ENDDO
            ENDIF
            ENDIF
            LEN_DATA = LEN_DATA+LOOKUP_INTF(LBLREC,VAR,JINTF)

          ENDIF

          IADDR = IADDR+LOOKUP_INTF(LBLREC,VAR,JINTF)
        ENDDO          ! VAR

CL 4.0 Write out headers/data

CL 4.1 Fixed length header

        IADDR = 0
        CALL SETPOS (NFTOUT,IADDR,ICODE)
        CALL BUFFOUT(NFTOUT,FIXHD_INTF(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)

C Check for I/O Errors

        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.LEN_FIXHD) THEN
          CALL IOERROR('buffer out of fixed length header',A_IO,LEN_IO,
     &                  LEN_FIXHD)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=2
          GO TO 999   !  Return
        END IF

CL 4.2 Integer constants

        CALL BUFFOUT (NFTOUT,INTHD_INTF(1,JINTF),
     &                PP_LEN_INTHD,LEN_IO,A_IO)

C Check for I/O Errors

        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.PP_LEN_INTHD) THEN
          CALL IOERROR('buffer out of integer header',A_IO,LEN_IO,
     &                  PP_LEN_INTHD)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=3
          GO TO 999   !  Return
        END IF

CL 4.3 PP headers in LOOKUP table
        CALL SETPOS(NFTOUT,LOOKUP_START,ICODE)
        CALL BUFFOUT(NFTOUT,LOOKUP_INTF(1,1,JINTF),
     &               LEN1_LOOKUP*INTF_LOOKUPS,LEN_IO,A_IO)

C Check for I/O Errors

        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.LEN1_LOOKUP*INTF_LOOKUPS) THEN
          CALL IOERROR('buffer out of PP header',A_IO,LEN_IO,
     &                  LEN1_LOOKUP*INTF_LOOKUPS)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=4
          GO TO 999   !  Return
        END IF

CL 4.4 Interface data
C       Determine position in data section

        DATA_START =
     &   lookup_intf(lbegin, 1, jintf)
c--round this disk length to a multiple of the sector size
        len_data=((len_data+um_sector_size-1)/
     &    um_sector_size)*um_sector_size
        CALL SETPOS(NFTOUT,DATA_START,ICODE)
        CALL BUFFOUT(NFTOUT,INTF_DATA(1),LEN_DATA,LEN_IO,A_IO)

C Check for I/O Errors

        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.LEN_DATA) THEN
          CALL IOERROR('buffer out of boundary data',A_IO,LEN_IO,
     &                  LEN_DATA)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=51
          GO TO 999   !  Return
        END IF


CL 5.    Close boundary output file if reinitialised during run
      IF (FT_STEPS(NFTOUT).GT.0) THEN
        LEN_PPNAME=LEN(PPNAME)
        CALL FILE_CLOSE(NFTOUT,PPNAME,LEN_PPNAME,1,0,ICODE)
      END IF

CL 6.  Update FT_LASTFIELD
      FT_LASTFIELD(NFTOUT) = FT_LASTFIELD(NFTOUT) + 1

 999  RETURN
      END
