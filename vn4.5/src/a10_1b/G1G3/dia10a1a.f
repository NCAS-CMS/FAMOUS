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
CLL   SUBROUTINE DIAG10_A --------------------------------------------
CLL
CLL  PURPOSE: Calculate diagnostics from section 10 before call to
CLL           THETL_QT.
CLL
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL 3.2 26/07/93 CHANGE DIMENSION OF SF TO INCLUDE (0:NITEMS,  R.RAWLINS
CLL 3.4 29/04/94 : Correct calculations of temperature from theta
CLL                (was T=theta/p_exner now T=theta*p_exner) R Stratton
CLL     11/10/94 : Correct calls to COPYDIAG_3D. R A Stratton
!   4.2 25/04/95 : Scale many of the output fields by 1.0e-6 to avoid
!                  problems with partial sums of the field reaching
!                  numbers too big to be packed to 32 bits (ie > 1.e9)
!LL 4.3     11/02/97 Added ARGFLDPT and ARGPPX arguments   P.Burton
!LL 4.4  10/09/97 : Correct error introduced by MPP
!LL 4.5  28/10/98   Introduce Single Column Model. J-C Thil.
CLL
CLL   Programming standard: U M DOC  Paper NO. 4,
CLL
CLL   Logical components covered : D3111
CLL
CLL   Project task: P1
CLL
CLL   External documentation: U.M. Doc. Paper 10. Appendix 3.
CLL
CLLEND---------------------------------------------------------------

C*L  ARGUMENTS:------------------------------------------------------

      SUBROUTINE DIAG10_A(
     &                    PSTAR,PSTAR_OLD,U_ADJ,V_ADJ,Q,ETADOT,
     &                    THETA,P_EXNER,RS,SEC_U_LATITUDE,
     &                    ROW_LENGTH,P_LEVELS,Q_LEVELS,P_FIELD,
     &                    U_FIELD,AK,BK,AKH,BKH,ADVECTION_TIMESTEP,
     &                    FIRST_POINT,LAST_POINT,
     &                    NSECTS,NITEMS,TOTITEMS,NUM_STASH_LEVELS,
     &                    NUM_LEVEL_LISTS,LEN_STLIST,STASHLEN,SF,
     &                    STINDEX,STLIST,SI,STASH_LEVELS,STASHWORK,
     &                    FIELD,WORK_LENGTH,
     &                    im_ident,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
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
     &                    ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     &  P_FIELD            !IN  1ST DIMENSION OF FIELD OF PSTAR
     &, U_FIELD            !IN  1ST DIMENSION OF FIELD OF U,V
     &, ROW_LENGTH         !IN  NUMBER OF POINTS PER ROW
     &, P_LEVELS           !IN  NUMBER OF PRESSURE LEVELS
     &, Q_LEVELS           !IN  NUMBER OF WET LEVELS
     &, FIRST_POINT        !IN  FIRST POINT OUTPUT REQUIRED FOR.
     &, LAST_POINT         !IN  LAST POINT OUTPUT REQUIRED FOR.
     &, WORK_LENGTH        !IN  SIZE OF DYNAMICALLY ALLOCATED WORKSPACE

      INTEGER
     &  im_ident           !IN : Internal model indent

! Comdeck TYPFLDPT
! Variables which point to useful positions in a horizontal field

      INTEGER
     &  FIRST_ROW        ! First updatable row on field
     &, TOP_ROW_START    ! First point of north-pole (global) or
!                        ! Northern (LAM) row
!                        ! for processors not at top of LPG, this
!                        ! is the first point of valid data
!                        ! (ie. Northern halo).
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
!                        ! for processors not at base of LPG, this
!                        ! is the start of the last row of valid data
!                        ! (ie. Southern halo).
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
! For the last three variables, these indexes are the start points
! and end points of "local" data - ie. missing the top and bottom
! halo regions.
     &, FIRST_VALID_PT   ! first valid point of data on field
     &, LAST_P_VALID_PT  ! last valid point of data on field
     &, LAST_U_VALID_PT  ! last valid point of data on field
     &, VALID_P_ROWS     ! number of valid rows of P data
     &, VALID_U_ROWS     ! number of valid rows of U data
     &, START_POINT_NO_HALO
!                        ! first non-polar point of field (misses
!                        ! halo for MPP code)
     &, START_POINT_INC_HALO
!                        ! first non-polar point of field (includes
!                        ! halo for MPP code)
     &, END_P_POINT_NO_HALO
!                        ! last non-polar point of P field (misses
!                        ! halo for MPP code)
     &, END_P_POINT_INC_HALO
!                        ! last non-polar point of P field (includes
!                        ! halo for MPP code)
     &, END_U_POINT_NO_HALO
!                        ! last non-polar point of U field (misses
!                        ! halo for MPP code)
     &, END_U_POINT_INC_HALO
!                        ! last non-polar point of U field (includes
!                        ! halo for MPP code)
     &, FIRST_ROW_PT     ! first data point along a row
     &, LAST_ROW_PT      ! last data point along a row
! For the last two variables, these indexes are the start and
! end points along a row of the "local" data - ie. missing out
! the east and west halos
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!

     &, MY_PROC_ID         ! my processor id
     &, NP_PROC_ID         ! processor number of North Pole Processor
     &, SP_PROC_ID         ! processor number of South Pole Processor
     &, GC_ALL_GROUP       ! group id of group of all processors
     &, GC_ROW_GROUP       ! group id of group of all processors on this
!                          ! processor row
     &, GC_COL_GROUP       ! group id of group of all processors on this
!                          ! processor column
     &, N_PROCS            ! total number of processors

     &, EW_Halo            ! Halo size in the EW direction
     &, NS_Halo            ! Halo size in the NS direction

     &, halo_4th           ! halo size for 4th order calculations
     &, extra_EW_Halo      ! extra halo size required for 4th order
     &, extra_NS_Halo      ! extra halo size required for 4th order
     &, LOCAL_ROW_LENGTH   ! size of local row
     &, FIRST_GLOBAL_ROW_NUMBER
!                          ! First row number on Global Grid    

! Variables which indicate if special operations are required at the
! edges.
      LOGICAL
     &  at_top_of_LPG    ! Logical variables indicating if this
     &, at_right_of_LPG  ! processor is at the edge of the Logical
     &, at_base_of_LPG   ! Processor Grid and should process its edge
     &, at_left_of_LPG   ! data differently.

! End of comdeck TYPFLDPT
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
      INTEGER
     &  ICODE              !OUT RETURN CODE. NON-ZERO IF ERROR-DETECTED

      CHARACTER
     &  CMESSAGE*(*)       !OUT ERROR MESSAGE

C INPUT DATA

      REAL
     &  PSTAR(P_FIELD)          !IN PRIMARY MODEL ARRAY FOR PSTAR FIELD
     &, PSTAR_OLD(P_FIELD)      !IN PSTAR FIELD AT PREVIOUS TIMESTEP.
     &, P_EXNER(P_FIELD,P_LEVELS+1) !IN  EXNER PRESS ON 1/2 LVLS
     &, THETA(P_FIELD,P_LEVELS) !IN PRIMARY MODEL ARRAY FOR THETA FIELD
     &, U_ADJ(U_FIELD,P_LEVELS) !IN MEAN U OVER ADJUSTMENT STEPS
     &, V_ADJ(U_FIELD,P_LEVELS) !IN MEAN V OVER ADJUSTMENT STEPS
     &, Q(P_FIELD,Q_LEVELS)     !IN PRIMARY MODEL ARRAY FOR HUMIDITY
     &, RS(P_FIELD,P_LEVELS)    !IN EFFECTIVE RADIUS OF EARTH.
     &, ETADOT(P_FIELD,P_LEVELS)!IN VERTICAL VELOCITY.
     &, SEC_U_LATITUDE(U_FIELD) !IN 1./(COS(LAT)) AT U POINTS.

      REAL
     &  AKH(P_LEVELS+1)         !IN  LAYER THICKNESS
     &, BKH(P_LEVELS+1)         !IN  LAYER THICKNESS
     &, AK (P_LEVELS)           !IN  VALUE AT LAYER CENTRE
     &, BK (P_LEVELS)           !IN  VALUE AT LAYER CENTRE
     &, ADVECTION_TIMESTEP      !IN  ADVECTION TIMESTEP.

      REAL
     &  FIELD(P_FIELD*P_LEVELS) ! WORK-SPACE FOR OUTPUT FIELD

C STASH REQUIREMENTS.

      INTEGER
     &  NSECTS             !IN NO OF PROCESSING SECTIONS (MASTER PCRS)
     &, NITEMS             !IN MAX NO OF STASH ITEMS IN A SECTION
     &, TOTITEMS           !IN MAX NO OF TOTAL STASH ITEMS
     &, NUM_STASH_LEVELS   !IN MAX NUMBER OF LEVELS IN A LEVELS LIST
     &, NUM_LEVEL_LISTS    !IN MAX NUMBER OF LEVELS LIST
     &, LEN_STLIST         !IN LENGTH OF LIST OF ITEMS FROM STASH
     &, STASHLEN           !IN SIZE OF STASHWORK

      INTEGER
     &  STINDEX(2,NITEMS,0:NSECTS)    !IN
     &, STLIST(LEN_STLIST,TOTITEMS)   !IN
     &, SI(NITEMS,0:NSECTS)           !IN
     &, STASH_LEVELS(NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS) !IN

      LOGICAL
     &  SF(0:NITEMS,0:NSECTS)        !IN

      REAL
     &  STASHWORK(STASHLEN) !INOUT. WORK SPACE HOLDING STASH OUTPUT.

C*--------------------------------------------------------------------

C*L   DEFINE LOCAL ARRAYS AND VARIABLES USED IN THIS ROUTINE----------
C DEFINE LOCAL ARRAYS: 2 ARE REQUIRED.
      REAL
     &  VELOCITY(WORK_LENGTH)      ! WORK-SPACE FOR INTERPOLATED
     &                             ! WIND FIELD.
      INTEGER
     &  FIRST_U,FIRST_P       ! first point for COPYDIAG for U & P grids
     & ,LAST_U,LAST_P         ! last point for COPYDIAG for U & P grids

      LOGICAL
     &  LIST(P_LEVELS)

C*--------------------------------------------------------------------

C DEFINE LOCAL VARIABLES
      REAL
     &  RECIP_TIMESTEP
     &, EARTH_RADIUS_INVERSE
     &, SCALAR
     &, SCALAR_A
     &, SCALAR_B
     &, PKP1,PK              !  Pressures at half levels k+1 and k
     &, P_EXNER_FULL         !  Exner Pressure at full model level
     &, TEMP_I,TEMP_IP1      !  Temperatures at points/rows i and i+1
     & ,FACTOR                 ! scaling factor

      INTEGER
     &  I,K,K1,LEVEL

CLL   COMDECK C_DG10_1

CL    Define local variables for STASH processing in routines
CL    DIAG10_A and DIAG10_B. Values are assigned in comdeck C_DG10_2

      LOGICAL
     &  L_UADJ_DP        !TRUE IF DIAGNOSTIC REQUIRED.
     &, L_VADJ_DP        !           "
     &, L_EFF_RADIUS     !           "
     &, L_ETADOT         !           "
     &, L_PRESS_TEND     !           "
     &, L_GEOPOTENTIAL   !           "
     &, L_UADJ_T_DP      !           "
     &, L_VADJ_T_DP      !           "
     &, L_UADJ_Q_DP      !           "
     &, L_VADJ_Q_DP      !           "

      LOGICAL
     &  L_UADJ_TL_DP      !TRUE IF DIAGNOSTIC REQUIRED.
     &, L_VADJ_TL_DP      !           "
     &, L_UADJ_QT_DP      !           "
     &, L_VADJ_QT_DP      !           "
     &, L_UADJ_U_DP       !           "
     &, L_VADJ_U_DP       !           "
     &, L_UADJ_V_DP       !           "
     &, L_VADJ_V_DP       !           "
     &, L_UADJ_GEOPOT_DP  !           "
     &, L_VADJ_GEOPOT_DP  !           "
     &, L_UADJ_ENERGY_DP  !           "
     &, L_VADJ_ENERGY_DP  !           "

      INTEGER
     &  INDEX_UADJ_DP     ! LEVELS LIST INDEX ADDRESS
     &, INDEX_VADJ_DP     !           "
     &, INDEX_EFF_RADIUS  !           "
     &, INDEX_ETADOT      !           "
     &, INDEX_UADJ_T_DP   !           "
     &, INDEX_VADJ_T_DP   !           "
     &, INDEX_UADJ_Q_DP   !           "
     &, INDEX_VADJ_Q_DP   !           "

      INTEGER
     &  INDEX_UADJ_TL_DP                  ! LEVELS LIST INDEX ADDRESS
     &, INDEX_VADJ_TL_DP                  !           "
     &, INDEX_UADJ_QT_DP                  !           "
     &, INDEX_VADJ_QT_DP                  !           "
     &, INDEX_UADJ_U_DP                   !           "
     &, INDEX_VADJ_U_DP                   !           "
     &, INDEX_UADJ_V_DP                   !           "
     &, INDEX_VADJ_V_DP                   !           "
     &, INDEX_UADJ_GEOPOT_DP              !           "
     &, INDEX_VADJ_GEOPOT_DP              !           "
     &, INDEX_UADJ_ENERGY_DP              !           "
     &, INDEX_VADJ_ENERGY_DP              !           "

      INTEGER
     &  LOC_UADJ_DP        !LOCATION IN STASHWORK OF OUTPUT.
     &, LOC_VADJ_DP        !           "
     &, LOC_EFF_RADIUS     !           "
     &, LOC_ETADOT         !           "
     &, LOC_PRESS_TEND     !           "
     &, LOC_UADJ_T_DP      !           "
     &, LOC_VADJ_T_DP      !           "
     &, LOC_UADJ_Q_DP      !           "
     &, LOC_VADJ_Q_DP      !           "

      INTEGER
     &  LOC_UADJ_TL_DP         !LOCATION IN STASHWORK OF OUTPUT.
     &, LOC_VADJ_TL_DP         !           "
     &, LOC_UADJ_QT_DP         !           "
     &, LOC_VADJ_QT_DP         !           "
     &, LOC_UADJ_U_DP          !           "
     &, LOC_VADJ_U_DP          !           "
     &, LOC_UADJ_V_DP          !           "
     &, LOC_VADJ_V_DP          !           "
     &, LOC_UADJ_GEOPOT_DP     !           "
     &, LOC_VADJ_GEOPOT_DP     !           "
     &, LOC_UADJ_ENERGY_DP     !           "
     &, LOC_VADJ_ENERGY_DP     !           "
     &, LOC_GEOPOTENTIAL       !           "

CLL   End of comdeck C_DG10_1

C     Get UM constants
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------

C*L   EXTERNAL SUBROUTINES CALLED ------------------------------------
      EXTERNAL COPYDIAG_3D,SET_LEVELS_LIST,COPYDIAG
C*--------------------------------------------------------------------

C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------

C     Comdeck C_DG10_2 initialises local variables defined in C_DG10_1
CLL   COMDECK C_DG10_2

CL    Set up pointers and switches for STASH processing in routines
CL    DIAG10_A and DIAG10_B. Variables defined in comdeck C_DG10_1

      L_UADJ_DP         = SF(201,10)
      L_VADJ_DP         = SF(202,10)
      L_EFF_RADIUS      = SF(203,10)
      L_ETADOT          = SF(204,10)
      L_PRESS_TEND      = SF(205,10)
      L_GEOPOTENTIAL    = SF(206,10)
      L_UADJ_T_DP       = SF(207,10)
      L_VADJ_T_DP       = SF(208,10)
      L_UADJ_Q_DP       = SF(209,10)
      L_VADJ_Q_DP       = SF(210,10)
      L_UADJ_TL_DP      = SF(211,10)
      L_VADJ_TL_DP      = SF(212,10)
      L_UADJ_QT_DP      = SF(213,10)
      L_VADJ_QT_DP      = SF(214,10)
      L_UADJ_U_DP       = SF(215,10)
      L_VADJ_U_DP       = SF(216,10)
      L_UADJ_V_DP       = SF(217,10)
      L_VADJ_V_DP       = SF(218,10)
      L_UADJ_GEOPOT_DP  = SF(219,10)
      L_VADJ_GEOPOT_DP  = SF(220,10)
      L_UADJ_ENERGY_DP  = SF(221,10)
      L_VADJ_ENERGY_DP  = SF(222,10)

      INDEX_UADJ_DP             = STINDEX(1,201,10)
      INDEX_VADJ_DP             = STINDEX(1,202,10)
      INDEX_EFF_RADIUS          = STINDEX(1,203,10)
      INDEX_ETADOT              = STINDEX(1,204,10)
      INDEX_UADJ_T_DP           = STINDEX(1,207,10)
      INDEX_VADJ_T_DP           = STINDEX(1,208,10)
      INDEX_UADJ_Q_DP           = STINDEX(1,209,10)
      INDEX_VADJ_Q_DP           = STINDEX(1,210,10)
      INDEX_UADJ_TL_DP          = STINDEX(1,211,10)
      INDEX_VADJ_TL_DP          = STINDEX(1,212,10)
      INDEX_UADJ_QT_DP          = STINDEX(1,213,10)
      INDEX_VADJ_QT_DP          = STINDEX(1,214,10)
      INDEX_UADJ_U_DP           = STINDEX(1,215,10)
      INDEX_VADJ_U_DP           = STINDEX(1,216,10)
      INDEX_UADJ_V_DP           = STINDEX(1,217,10)
      INDEX_VADJ_V_DP           = STINDEX(1,218,10)
      INDEX_UADJ_GEOPOT_DP      = STINDEX(1,219,10)
      INDEX_VADJ_GEOPOT_DP      = STINDEX(1,220,10)
      INDEX_UADJ_ENERGY_DP      = STINDEX(1,221,10)
      INDEX_VADJ_ENERGY_DP      = STINDEX(1,222,10)

      LOC_UADJ_DP        = SI(201,10)
      LOC_VADJ_DP        = SI(202,10)
      LOC_EFF_RADIUS     = SI(203,10)
      LOC_ETADOT         = SI(204,10)
      LOC_PRESS_TEND     = SI(205,10)
      LOC_GEOPOTENTIAL   = SI(206,10)
      LOC_UADJ_T_DP      = SI(207,10)
      LOC_VADJ_T_DP      = SI(208,10)
      LOC_UADJ_Q_DP      = SI(209,10)
      LOC_VADJ_Q_DP      = SI(210,10)
      LOC_UADJ_TL_DP     = SI(211,10)
      LOC_VADJ_TL_DP     = SI(212,10)
      LOC_UADJ_QT_DP     = SI(213,10)
      LOC_VADJ_QT_DP     = SI(214,10)
      LOC_UADJ_U_DP      = SI(215,10)
      LOC_VADJ_U_DP      = SI(216,10)
      LOC_UADJ_V_DP      = SI(217,10)
      LOC_VADJ_V_DP      = SI(218,10)
      LOC_UADJ_GEOPOT_DP = SI(219,10)
      LOC_VADJ_GEOPOT_DP = SI(220,10)
      LOC_UADJ_ENERGY_DP = SI(221,10)
      LOC_VADJ_ENERGY_DP = SI(222,10)

CLL   End of comdeck C_DG10_2

CL--------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL--------------------------------------------------------------------
      FIRST_U = FIRST_FLD_PT
      FIRST_P = FIRST_FLD_PT
      LAST_U  = LAST_U_FLD_PT
      LAST_P  = LAST_P_FLD_PT

CL -------------------------------------------------------------------
CL SECTION 1.  DIAGNOSTICS INVOLVING MEAN U OVER ADJUSTMENT STEP.
CL -------------------------------------------------------------------

      EARTH_RADIUS_INVERSE = 1./A
      FACTOR=1.0e-6

C --------------------------------------------------------------------
CL SECTION 1.1 MEAN PRESSURE WEIGHTED U OVER ADJUSTMENT STEPS.
C --------------------------------------------------------------------

      IF (L_UADJ_DP) THEN
CL   REMOVE RADIUS OF EARTH FROM U FIELD.
C MINUS SIGN SETS DELTA P TO POSITIVE VALUE.
       DO 110 K=1,P_LEVELS
          K1 = (K-1)*U_FIELD
          DO I=FIRST_U,LAST_U
            FIELD(K1+I) = -U_ADJ(I,K)*EARTH_RADIUS_INVERSE
          END DO
 110   CONTINUE

        CALL COPYDIAG_3D (STASHWORK(LOC_UADJ_DP),FIELD,FIRST_U,
     &                    LAST_U,U_FIELD,ROW_LENGTH,P_LEVELS,
     &                    STLIST(1,INDEX_UADJ_DP),LEN_STLIST,
     &                    STASH_LEVELS,NUM_STASH_LEVELS+1,
     &                    im_ident,10,201,
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
     &                    ICODE,CMESSAGE)
        IF(ICODE.GT.0) THEN
          RETURN
        END IF
      END IF

CL CHECK TO SEE IF ANY U DIAGNOSTICS REQUESTED WHICH NEED U_ADJ TO
CL BE INTERPOLATED.

      IF(L_UADJ_T_DP.OR.L_UADJ_Q_DP) THEN

C --------------------------------------------------------------------
CL SECTION 1.2 INTERPOLATE U TO C-GRID U POINTS.
C --------------------------------------------------------------------

C MINUS SIGN SETS DELTA P TO POSITIVE VALUE.
        DO 120 K=1,P_LEVELS
          K1 = (K-1)*P_FIELD
          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            VELOCITY(K1+I) = -.5*(U_ADJ(I,K) + U_ADJ(I-ROW_LENGTH,K))
     &                                      *EARTH_RADIUS_INVERSE
          END DO
C SET POLAR VALUES EQUAL TO VALUE ON ADJACENT ROW.

          IF (at_top_of_LPG) THEN
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              VELOCITY(K1+I) = -U_ADJ(I,K)*EARTH_RADIUS_INVERSE
            ENDDO
          ENDIF

          IF (at_base_of_LPG) THEN
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              VELOCITY(K1+I) = -U_ADJ(I-ROW_LENGTH,K)*
     &                          EARTH_RADIUS_INVERSE
            ENDDO
          ENDIF
 120    CONTINUE

C --------------------------------------------------------------------
CL SECTION 1.3 MEAN PRESSURE WEIGHTED U OVER ADJUSTMENT STEPS
CL                      * TEMPERATURE.
C --------------------------------------------------------------------

        IF (L_UADJ_T_DP) THEN
          DO 130 K=1,P_LEVELS
            K1 = (K-1)*P_FIELD
            DO I=FIRST_P,LAST_P-1

              PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I)
              PK   = AKH(K)   + BKH(K)  *PSTAR(I)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I,K+1),P_EXNER(I,K),PKP1,PK,KAPPA)
              TEMP_I = THETA(I,K) * P_EXNER_FULL

              PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I+1)
              PK   = AKH(K)   + BKH(K)  *PSTAR(I+1)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I+1,K+1),P_EXNER(I+1,K),PKP1,PK,KAPPA)
              TEMP_IP1 = THETA(I+1,K) * P_EXNER_FULL

              FIELD(K1+I) = VELOCITY(K1+I) * 0.5 * (TEMP_I + TEMP_IP1)
     &                       *factor

            END DO
! Set last point of field (halo) to a valid number
            FIELD(K1+LAST_P)=FIELD(K1+LAST_P-1)
 130    CONTINUE

          CALL COPYDIAG_3D (STASHWORK(LOC_UADJ_T_DP),FIELD,FIRST_P,
     &                      LAST_P,P_FIELD,ROW_LENGTH,P_LEVELS,
     &                      STLIST(1,INDEX_UADJ_T_DP),LEN_STLIST,
     &                      STASH_LEVELS,
     &                      NUM_STASH_LEVELS+1,
     &                      im_ident,10,207,
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
     &                      ICODE,CMESSAGE)
          IF(ICODE.GT.0) THEN
            RETURN
          END IF
        END IF

C --------------------------------------------------------------------
CL SECTION 1.4 MEAN PRESSURE WEIGHTED U OVER ADJUSTMENT STEPS
CL                      * HUMIDITY.
C --------------------------------------------------------------------

        IF (L_UADJ_Q_DP) THEN
          DO 140 K=1,Q_LEVELS
            K1 = (K-1)*P_FIELD
            DO I=FIRST_P,LAST_P-1
              FIELD(K1+I) = VELOCITY(K1+I)*.5* (Q(I,K)+Q(I+1,K))
            END DO
! Set last point of field (halo) to a valid number
            FIELD(K1+LAST_P)=FIELD(K1+LAST_P-1)
 140      CONTINUE

          CALL COPYDIAG_3D (STASHWORK(LOC_UADJ_Q_DP),FIELD,FIRST_P,
     &                      LAST_P,P_FIELD,ROW_LENGTH,Q_LEVELS,
     &                      STLIST(1,INDEX_UADJ_Q_DP),LEN_STLIST,
     &                      STASH_LEVELS,
     &                      NUM_STASH_LEVELS+1,
     &                      im_ident,10,209,
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
     &                      ICODE,CMESSAGE)
          IF(ICODE.GT.0) THEN
            RETURN
          END IF
        END IF

C END IF FOR U DIAGNOSTICS
      END IF

CL -------------------------------------------------------------------
CL SECTION 2.  DIAGNOSTICS INVOLVING MEAN V OVER ADJUSTMENT STEP.
CL -------------------------------------------------------------------

C --------------------------------------------------------------------
CL SECTION 2.1 MEAN PRESSURE WEIGHTED V OVER ADJUSTMENT STEPS.
C --------------------------------------------------------------------

      IF (L_VADJ_DP) THEN
CL   REMOVE RADIUS OF EARTH * COSINE OF LATITUDE FROM V FIELD.
C MINUS SIGN SETS DELTA P TO POSITIVE VALUE.

        DO 210 K=1,P_LEVELS
          K1 = (K-1)*U_FIELD
          DO I=FIRST_U,LAST_U
            FIELD(K1+I) = -V_ADJ(I,K)*EARTH_RADIUS_INVERSE
     &                                       *SEC_U_LATITUDE(I)
          END DO
 210    CONTINUE

        CALL COPYDIAG_3D(STASHWORK(LOC_VADJ_DP),FIELD,FIRST_U,
     &                    LAST_U,U_FIELD,ROW_LENGTH,P_LEVELS,
     &                    STLIST(1,INDEX_VADJ_DP),LEN_STLIST,
     &                    STASH_LEVELS,
     &                    NUM_STASH_LEVELS+1,
     &                    im_ident,10,202,
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
     &                    ICODE,CMESSAGE)
        IF(ICODE.GT.0) THEN
          RETURN
        END IF
      END IF

CL CHECK TO SEE IF ANY V DIAGNOSTICS REQUESTED WHICH NEED V_ADJ TO
CL BE INTERPOLATED.

      IF(L_VADJ_T_DP.OR.L_VADJ_Q_DP) THEN

C --------------------------------------------------------------------
CL SECTION 2.2 INTERPOLATE V TO C-GRID V POINTS.
C --------------------------------------------------------------------

C MINUS SIGN SETS DELTA P TO POSITIVE VALUE.
        DO 220 K=1,P_LEVELS
          K1 = (K-1)*U_FIELD
          DO I=FIRST_U+1,LAST_U
            VELOCITY(K1+I)= -.5*(V_ADJ(I,K)*SEC_U_LATITUDE(I)
     &                                      +V_ADJ(I-1,K)
     &                                      *SEC_U_LATITUDE(I-1))
     &                                      *EARTH_RADIUS_INVERSE
          END DO
! Set first point of field (halo) to a valid number
          VELOCITY(K1+FIRST_U)= VELOCITY(K1+FIRST_U+1)
 220    CONTINUE

C --------------------------------------------------------------------
CL SECTION 2.3 MEAN PRESSURE WEIGHTED V OVER ADJUSTMENT STEPS
CL                      * TEMPERATURE.
C --------------------------------------------------------------------

        IF (L_VADJ_T_DP) THEN
          DO 230 K=1,P_LEVELS
            K1 = (K-1)*U_FIELD
            DO I=FIRST_U,LAST_U

              PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I)
              PK   = AKH(K)   + BKH(K)  *PSTAR(I)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I,K+1),P_EXNER(I,K),PKP1,PK,KAPPA)
              TEMP_I = THETA(I,K) * P_EXNER_FULL

              PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I+ROW_LENGTH)
              PK   = AKH(K)   + BKH(K)  *PSTAR(I+ROW_LENGTH)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I+ROW_LENGTH,K+1),P_EXNER(I+ROW_LENGTH,K),
     &         PKP1,PK,KAPPA)
              TEMP_IP1 = THETA(I+ROW_LENGTH,K) * P_EXNER_FULL

              FIELD(K1+I) = VELOCITY(K1+I) * 0.5 * (TEMP_I + TEMP_IP1)
     &                      * FACTOR

            END DO
 230      CONTINUE

          CALL COPYDIAG_3D (STASHWORK(LOC_VADJ_T_DP),FIELD,FIRST_U,
     &                      LAST_U,U_FIELD,ROW_LENGTH,P_LEVELS,
     &                      STLIST(1,INDEX_VADJ_T_DP),LEN_STLIST,
     &                      STASH_LEVELS,
     &                      NUM_STASH_LEVELS+1,
     &                      im_ident,10,208,
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
     &                      ICODE,CMESSAGE)
          IF(ICODE.GT.0) THEN
            RETURN
          END IF
        END IF

C --------------------------------------------------------------------
CL SECTION 2.4 MEAN PRESSURE WEIGHTED V OVER ADJUSTMENT STEPS
CL                      * HUMIDITY.
C --------------------------------------------------------------------

        IF (L_VADJ_Q_DP) THEN
          DO 240 K=1,Q_LEVELS
            K1 = (K-1)*U_FIELD
            DO I=FIRST_U,LAST_U
              FIELD(K1+I) = VELOCITY(K1+I)*.5*
     &                                 (Q(I,K)+Q(I+ROW_LENGTH,K))
            END DO
 240      CONTINUE

          CALL COPYDIAG_3D (STASHWORK(LOC_VADJ_Q_DP),FIELD,FIRST_U,
     &                      LAST_U,U_FIELD,ROW_LENGTH,Q_LEVELS,
     &                      STLIST(1,INDEX_VADJ_Q_DP),LEN_STLIST,
     &                      STASH_LEVELS,
     &                      NUM_STASH_LEVELS+1,
     &                      im_ident,10,210,
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
     &                      ICODE,CMESSAGE)
          IF(ICODE.GT.0) THEN
            RETURN
          END IF
        END IF

C END IF FOR V DIAGNOSTICS
      END IF

CL -------------------------------------------------------------------
CL SECTION 3.  DIAGNOSTICS NOT INVOLVING MEAN HORIZONTAL VELOCITIES.
CL -------------------------------------------------------------------

C --------------------------------------------------------------------
CL SECTION 3.1 EFFECTIVE EARTH RADIUS AT MODEL LEVELS.
C --------------------------------------------------------------------

      IF (L_EFF_RADIUS) THEN
        CALL COPYDIAG_3D (STASHWORK(LOC_EFF_RADIUS),RS,FIRST_POINT,
     &                    LAST_POINT,P_FIELD,ROW_LENGTH,P_LEVELS,
     &                    STLIST(1,INDEX_EFF_RADIUS),LEN_STLIST,
     &                    STASH_LEVELS,
     &                    NUM_STASH_LEVELS+1,
     &                    im_ident,10,203,
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
     &                    ICODE,CMESSAGE)
        IF(ICODE.GT.0) THEN
          RETURN
        END IF
      END IF

C --------------------------------------------------------------------
CL SECTION 3.2 MEAN ETADOT IN ADJUSTMENT STEPS.
C --------------------------------------------------------------------

      IF (L_ETADOT) THEN
        CALL COPYDIAG_3D(STASHWORK(LOC_ETADOT),ETADOT,FIRST_POINT,
     &                   LAST_POINT,P_FIELD,ROW_LENGTH,P_LEVELS,
     &                   STLIST(1,INDEX_ETADOT),LEN_STLIST,
     &                   STASH_LEVELS,
     &                   NUM_STASH_LEVELS+1,
     &                   im_ident,10,204,
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
     &                   ICODE,CMESSAGE)
        IF(ICODE.GT.0) THEN
          RETURN
        END IF

CL CALL SET_LEVELS_LIST TO DETERMINE WHICH LEVELS OUTPUT ARRAY WAS
CL REQUESTED ON.

        CALL SET_LEVELS_LIST(P_LEVELS,LEN_STLIST,
     &                       STLIST(1,INDEX_ETADOT),
     &                       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,
     &                       CMESSAGE)
        K=0
CL CHECK TO SEE IF LEVEL 1 WAS REQUESTED AS THIS NEEDS SPECIAL TREATMENT
        IF(LIST(1)) THEN
          K=K+1
CL IF NOT STRATOSPHERIC MODEL THEN SET OUTPUT ETADOT FIELD TO ZERO.
          DO I=0,P_FIELD-1
            STASHWORK(LOC_ETADOT+I) = 0.
          END DO
        END IF
CL NOW REMOVE MASS-WEIGHT FROM ALL OTHER REQUESTED LEVELS.
        DO LEVEL=2,P_LEVELS
          IF(LIST(LEVEL)) THEN
C SCALAR HOLDS DELTA ETA / RADIUS OF EARTH SQUARED.
            SCALAR=((AK(LEVEL)-AK(LEVEL-1))/PREF
     &             +(BK(LEVEL)-BK(LEVEL-1)))
     &             /(A*A)
C SCALAR_A HOLDS DIFFERENCE IN AK PART OF DP
            SCALAR_A= AK(LEVEL)-AK(LEVEL-1)
C SCALAR_B HOLDS DIFFERENCE IN BK PART OF DP
            SCALAR_B= BK(LEVEL)-BK(LEVEL-1)
C REMOVE A*A*DP/DETA FROM ETADOT FIELD.
            DO I=FIRST_P-1,LAST_P-1
              STASHWORK(LOC_ETADOT+K*P_FIELD+I) =
     &                               STASHWORK(LOC_ETADOT+K*P_FIELD+I)
     &                               *SCALAR
     &                               /(SCALAR_A+SCALAR_B*PSTAR(I+1))
            END DO
            K=K+1
          END IF
        END DO
      END IF

C --------------------------------------------------------------------
CL SECTION 3.3 SURFACE PRESSURE TENDENCY.
C --------------------------------------------------------------------

      IF (L_PRESS_TEND) THEN
        RECIP_TIMESTEP = 1./ADVECTION_TIMESTEP
        DO I=FIRST_P,LAST_P
          FIELD(I) = (PSTAR(I) - PSTAR_OLD(I))*RECIP_TIMESTEP
        END DO
        CALL COPYDIAG(STASHWORK(LOC_PRESS_TEND),FIELD,FIRST_POINT,
     &   LAST_POINT,P_FIELD,ROW_LENGTH,
     &   im_ident,10,205,
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
     &   ICODE,CMESSAGE)

        IF (ICODE .GT. 0) RETURN
      ENDIF

C --------------------------------------------------------------------
CL SECTION 3.4 GEOPOTENTIAL.
CL             THIS DIAGNOSTIC ACCUMULATED IN ADJ_CTL OVER ALL P_LEVELS
CL             USED IN DIAG10_B TO CALCULATE ENERGY.
CL             ALREADY HELD IN STASHWORK.
C --------------------------------------------------------------------

CL    END OF ROUTINE DIAG10_A

      RETURN
      END
