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
CLL Subroutine REPLANCA ---------------------------------------------
CLL
CLL Purpose:  Updates ancillary fields as requested in FIELDCODE array.
CLL   Tests whether update is required for each field, allowing for
CLL   dependencies between fields. Uses LOOKUP array to find data for
CLL   appropriate time, reads a record and checks for current data
CLL   type. Reads second record if time interpolation required. Updates
CLL   the field. Under DEF RECON, the interface to the routine is
CLL   modified for use in the reconfiguration rather than the model.
CLL   Under DEF CAL360 the 360 day rather than the Gregorian calender
CLL   is used.
CLL
CLL Level 2 control routine for CRAY YMP
CLL
CLL C.Wilson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  22/02/93  Changes to allow updating of SLAB ref SST and ice
CLL                 ancillary fields (items 178,179) from SST/ice files.
CLL                 Correct bug if SST updated but not ice fraction.
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.2  15/04/93  Remove misleading warning messages if no time
CLL                  interpolation of SST
CLL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
CLL                  elapsed times in days & secs, for portability. TCJ
CLL  3.3  22/11/93  Source term and aerosol ancillary fields added.
CLL   3.3  17/11/93 Initializing Integer UPDATE_MONTHS (N.Farnon)
CLL   3.3  08/12/93  Extra argument for READFLDS. D. Robinson
CLL   3.4  17/06/94  DEF CAL360 replaced by LOGICAL LCAL360
CLL                  Argument LCAL360 passed to SEC2TIME, TIME2SEC
CLL                                                   S.J.Swarbrick
CLL  3.4  20/07/94  Improve time interpolation by using reference time
CLL                  for ancillary updating.   R.T.H.Barnes.
CLL  3.4  05/09/94  Add murk and user ancillary fields.  R.T.H.Barnes.
CLL   3.4  18/05/94  Allow prognostic slabtemp under sea ice.J Thomson
CLL   3.4  31/8/94   More error trapping.        (William Ingram)
CLL  4.0  06/09/95  Only print time interpolation diagnostics when it
CLL                 is really done.  RTHBarnes.
CLL   4.0  08/09/95   Allow time interpolation of ozone fields and
CLL                   cater for zonal/full fields. D. Robinson
CLL   4.0  29/11/95   Set land points to zero for surface currents
CLL                   and Heat convergence fields. D. Robinson
!     4.1  18/06/96   Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
CLL   4.1  22/05/96   Replace list of ancillary fields with call to
CLL                   comdeck CANCLSTA. D. Robinson.
CLL   4.2  08/11/96   Initialise fields to ensure haloes contain data
CLL                   for time interpolation for MPP runs. D. Robinson
!     4.1  16/12/96   Check ancillary files for non-constant polar rows,
!                     in reconfiguration only. Correct if LPOLARCHK=T
!                     Author D.M. Goddard
!LL   4.4  07/07/97   Alter SST and ice updating for AMIPII runs
!LL                   R A Stratton
!     4.4  13/11/97   Ancilary fields 72 - 89 no longer used for
!                     for user defined ancillaries. Code altered to
!                     ensure correct treatment of these fields.
!                     Author D.M. Goddard
CLL   4.4  25/07/97   (Reconfiguration only). Prevent failure when 
CLL                   non-constant polar values for ancillary files are 
CLL                   corrected.                         R. Rawlins
CLL   4.5  22/04/98   Add control of new NH3, soot aerosol emission
CLL                   ancillary fields. Plus minor message changes.
CLL                   R.Rawlins
!     4.5  22/10/98   Increase number of user ancillary fields by
!                     deleting existing four fields 68 - 72 and
!                     adding twenty to end 90 - 109
!                     Author D.M. Goddard
!LL   4.5  22/01/98   Correct level of second field read in for time
!LL                   interpolation.  D. Robinson
CLL
CLL Programing standard : UMDP no 3, version no 2, dated 07/09/90
CLL
CLL Logical component covered : C71
CLL
CLL System task : C7
CLL
CLL   External Documentation: UMDP no C7
CLL
CLLEND-------------------------------------------------------------

       SUBROUTINE REPLANCA(I_YEAR,I_MONTH,I_DAY,I_HOUR,
     &                     I_MINUTE,I_SECOND,I_DAY_NUMBER,
     &                     ANCIL_REFTIME,OFFSET_STEPS,
     &                     P_FIELD,P_ROWS,U_FIELD,D1,LAND,
     &                     A_STEP,LAND_FIELD,STEPS_PER_HR,
     &                     ICE_FRACTION,TSTAR,TSTAR_ANOM,
     &                     NS_SPACE,FIRST_LAT,
     &                     LEN1_LOOKUP,LEN_FIXHD,LEN_INTHD,
     &                     LEN_REALHD,LEN_D1,FIXHD,INTHD,REALHD,
     &                     LOOKUP,RLOOKUP,FTNANCIL,LOOKUP_START,
     &                     NDATASETS,NLOOKUPS,
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
     &                     ICODE,CMESSAGE,LCAL360)        ! Intent Out

      IMPLICIT NONE

      LOGICAL LCAL360

      INTEGER
     &       I_YEAR,            ! Curent Model Time
     &       I_MONTH,           !   "      "     "
     &       I_DAY,             !   "      "     "
     &       I_HOUR,            !   "      "     "
     &       I_MINUTE,          !   "      "     "
     &       I_SECOND,          !   "      "     "
     &       I_DAY_NUMBER,
     &       ANCIL_REFTIME(6),  ! Reference time for ancillary updating
     &       OFFSET_STEPS,      ! Offset in timesteps of ref. from basis


     &       A_STEP,LAND_FIELD,STEPS_PER_HR,


     &       P_FIELD,           ! Size of horizontal fields
     &       P_ROWS,            !
     &       U_FIELD,           !   "  "      "         "
     &       NDATASETS,         ! Number of ancillary datasets
     &       NLOOKUPS,          ! Number of lookup tables
     &       LEN_D1             ! Size of primary data array

      INTEGER
     &       LEN1_LOOKUP,       ! First dimension of lookup table
     &       LEN_FIXHD,         ! Length of headers in data sets
     &       LEN_INTHD,         !
     &       LEN_REALHD, !
     &       FIXHD(LEN_FIXHD,NDATASETS),  ! Data set headers
     &       INTHD(LEN_INTHD,NDATASETS),  !
     &       LOOKUP(LEN1_LOOKUP,NLOOKUPS),! Data set lookup tables
     &       FTNANCIL(NDATASETS),         ! FTN numbers of data sets
     &       LOOKUP_START(NDATASETS)      ! Start of lookup tables
C                                         ! referring to each data set.

      REAL
     &       D1(LEN_D1),        !INOUT  Primary data array used to hold
C                               !       all fields except TSTAR and
C                               !       ICE_FRACTION
     &       ICE_FRACTION(P_FIELD), !INOUT  Ice fraction, updated if
C                                   !       requested
     &       TSTAR(P_FIELD),    !INOUT  TSTAR, updated if requested
     &       TSTAR_ANOM(P_FIELD),!INOUT  SST anomaly,formed in recon;
     &                           !       added if requested in model run
     &       REALHD(LEN_REALHD,NDATASETS),
     &       RLOOKUP(LEN1_LOOKUP,NLOOKUPS)
     &       ,NS_SPACE         ! NS latitude spacing
     &       ,FIRST_LAT        ! latitude of first gridpoint

      LOGICAL
     &       LAND(P_FIELD)      ! Land sea mask

      INTEGER
     &       ICODE     ! Return code
     &      ,IOUNIT       !OUT I/O unit passed out in RECON mode

      CHARACTER*(80)
     &       CMESSAGE  ! Error message
C*
! Comdecks:------------------------------------------------------------
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
C*L -----------------COMDECK PHYSCONS----------------------------------
C
C   Purpose : contains physical constants required by the whole of the
C             model. It is made up of individual COMDECKS for sets of
C             of related constants, each routine can access one or
C             several of these COMDECKS seperately
C   System Component : F07
C   System task : Z
C  END
C*----------------------------------------------------------------------
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
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

C*L------------------COMDECK C_KT_FT-----------------------------------
      REAL KT2MS,    ! Knots to m/s conversion
     &     FT2M      ! Feet to meters conversion

      PARAMETER(
     & KT2MS=1852.0/3600.0,
     & FT2M =0.3048)
C*----------------------------------------------------------------------
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
C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------

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

C*L   Subroutines called;

      EXTERNAL
     &       TIME2SEC,
     &       READFLDS,
     &       T_INT,

     &       TO_LAND_POINTS,
     &       SEC2TIME,TIME_DF,

     &       T_INT_C

C*
C*L   Local real arrays

      REAL
     &       ANCIL1(P_FIELD),   ! Buffers to hold values of ancillary
C                               ! data for time interpolation.
     &       ANCIL2(P_FIELD),   !
     &       ANCIL_DATA(P_FIELD),! Field of ancillary data held prior
C                               ! to selective updating.
     &       SNOW_CHANGE(P_FIELD),! Fractional time of change of
C                               ! snow cover
     &       ICE_EXTENT(P_FIELD,2),! Fractional time of change
C                               ! of ice cover
     &       PRES_VALUE(P_FIELD) ! Prescribed value of data when
C                               ! controlling field is zero.
     &,      NO_ICE_EXTENT(P_FIELD)! Indicator for no sea ice
C                               ! =0 if ice cover
C*
C     Local variables

      INTEGER
     &       I,                 !
     &       I1,                !
     &       I2,                !
     &       I3,
     &       ID,                !
     &       IM,                !
     &       IY,                !
     &       K,
     &       FIELD,             ! Current field number.
     &       FILE               !

      INTEGER
     &       INTERVAL,          ! Interval between data times
     &       STEP,              ! Number of data times skipped.
     &       MONTHS,            ! Used in calculation of position
C                               ! of data required.
     &       HOURS,             !
     &       PERIOD,            ! Period of periodic data
     &       START_MONTH,       !
     &       LEVEL,             !
     &       NFTIN,             ! Current FTN number for ancillary field
     &       ANCIL_REF_DAYS,    ! Ancil.reference time in whole days
     &       ANCIL_REF_SECS,    ! Ancil.reference time in extra seconds
     &       DAY,SEC,           ! Times relative to reference time
     &       DAY1,SEC1,         ! Times relative to reference time
     &       INCR_SEC,          ! Increment in sec
     &       LEN
     &      ,IEND             
     &      ,II,ROW_LENGTH,J
      INTEGER
     &       I_YEAR1,            ! Copy of Curent Model Time year
     &       I_MONTH1,           !   "      "     "          month
     &       I_DAY1,             !   "      "     "          day
     &       I_HOUR1             !   "      "     "          hour

      INTEGER
     &       UPDATE_MONTHS      ! update frequency (months) if Gregorian
      LOGICAL
     &       LGREG_MONTHLY      ! True for Gregorian monthly updating
C
C *IF -DEF,CAL360
C
      INTEGER
     &       I_YEAR_BASIS,            ! Basis Model Time
     &       I_MONTH_BASIS,           !   "     "     "
     &       I_DAY_BASIS,             !   "     "     "
     &       I_HOUR_BASIS,            !   "     "     "
     &       I_MINUTE_BASIS,          !   "     "     "
     &       I_SECOND_BASIS,          !   "     "     "
     &       I_DAY_NUMBER_BASIS
C
C *ENDIF
C

      INTEGER
     &       I_YEAR_REF,              ! Reference Time
     &       I_MONTH_REF,             !    "       "
     &       I_DAY_REF,               !    "       "
     &       I_HOUR_REF,              !    "       "
     &       I_MINUTE_REF,            !    "       "
     &       I_SECOND_REF             !    "       "


      LOGICAL
     &       LINTERPOLATE,      ! Indicates whether time
C                               ! interpolation needed.
     &       LT_INT_C,          ! Indicates use of controlled time
C                               ! interpolation
     &       LMISMATCH,         ! Used in header checks
     &       LICE_FRACTION,     !
     &       LSNOW_DEPTH,       !
     &       SINGLE_TIME,       ! Indicates that only one time is
C                               ! available in data set
     &       PERIODIC,          ! Data set is periodic
     &       REGULAR            ! Interval between data times in
C                               ! dataset is regular in model timesteps.
     &       ,LICE_DEPTH

      REAL
     &       ZERO,            !
     &       TIME1,           ! Times if data used in time interpolation
     &       TIME2,           !
     &       TIME       !Target time for time interpolation
     &       ,LAT_P          ! latitude of point

CL    Internal structure

CL  List of Atmosphere & Slab Ancillary fields.
C*L--------------------COMDECK  CANCLSTA ---------------------------
!
! Purpose : Cross-Reference List of Ancillary Fields
!           Atmosphere and Slab. Ancillary Reference numbers,
!           Stash Codes, Model Codes and Logical file Numbers
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.1   20/05/96  New comdeck. Add new fields 72-81. D. Robinson.
!  4.3   18/3/97  Add new field 82.  William Ingram.
!  4.4   16/9/97  Add new fields 83-89.  Richard Betts.
!  4.5  22/04/98  Add new fields 68-70. R.Rawlins
!  4.5  23/03/98  Remove fields 68-71 and 78, and add fields 90-109.
!                 D.M. Goddard
!
! -------------------------------------------------------------------
!   Column A  : Ancillary Reference Number
!          B  : Internal Model Number
!          C  : Stash Item code
!          D  : Logical file number
!          E  : Field Description
!
!   A  B    C   D  E
!   ----------------------------------------------------------------
!   1  1   30   9  Land Sea Mask
!   2  1   33  10  Orography
!   3  1   34  10  Orographic Variance
!   4  1   35  10  Orographic gradient XX
!   5  1   36  10  Orographic gradient XY
!   6  1   37  10  Orographic gradient YY
!   7  1   60   1  Ozone
!   8  1   21   2  Soil Moisture
!   9  1   23   2  Snow Depth
!  10  1   20   3  Deep Soil Temperature
!  11  1   40   4  Vol SMC at Wilting
!  12  1   41   4  Vol SMC at Critical Point
!  13  1   42   4  Vol SMC at Field Capacity
!  14  1   43   4  Vol SMC at Saturation
!  15  1   44   4  Saturated Soil Conductivity
!  16  1   45   4  Eagleson's Exponent
!  17  1   46   4  Thermal Capacity
!  18  1   47   4  Thermal Conductivity
!  19  1   50   5  Vegetation Fraction
!  20  1   51   5  Root Depth
!  21  1   52   5  Snow Free Surface Albedo
!  22  1   53   5  Deep Snow Surface Albedo
!  23  1   54   5  Surface Resistance to Evaporation
!  24  1   55   5  Surface Capacity
!  25  1   56   5  Infiltration Factor
!  26  1   26   5  Surface Roughness (vegetation)
!  27  1   31   7  Sea-ice Fraction
!  28  1   24   6  Sea-surface Temperature
!  29  1   32   7  Sea-ice Thickness
!  30  1   28   8  Surface Currents : u-component
!  31  1   29   8  Surface Currents : v-component
!  32  1   93   9  Runoff coastal outflow point
!  33  3  177  11  Heat convergence (slab model)
!  34  1   19  10  Orographic roughness
!  35  1   48   4  Saturated soil water suction
!  36  1    9   2  Soil moisture in layers
!  37  3  178   6  SLAB reference SST climatology
!  38  3  179   7  SLAB reference GBM ice depth climatology
!  39  1   58  12  Sulphur dioxide emission
!  40  1   59  12  Dimethyl sulphur emission
!  41  1   88  13  Sulphate aerosol mass mixing ratio
!  42  1   87  13  Sulphuric acid aerosol mass mixing ratio
!  43  1   85  13  Soot aerosol mass mixing ratio
!  44  1   57  14  Multi-level murk source term emission
!  45  1   90  14  Multi-level murk concentration
!  46  1   17  10  Silhouette area of orography (orog. roughness)
!  47  1   18  10  Peak to trough height (for orog. roughness scheme)
!  48  1  301  15  User ancillary field 1
!  49  1  302  15  User ancillary field 2
!  50  1  303  15  User ancillary field 3
!  51  1  304  15  User ancillary field 4
!  52  1  305  15  User ancillary field 5
!  53  1  306  15  User ancillary field 6
!  54  1  307  15  User ancillary field 7
!  55  1  308  15  User ancillary field 8
!  56  1  309  15  User ancillary field 9
!  57  1  310  15  User ancillary field 10
!  58  1  311  15  User ancillary field 11
!  59  1  312  15  User ancillary field 12
!  60  1  313  15  User ancillary field 13
!  61  1  314  15  User ancillary field 14
!  62  1  315  15  User ancillary field 15
!  63  1  316  15  User ancillary field 16
!  64  1  317  15  User ancillary field 17
!  65  1  318  15  User ancillary field 18
!  66  1  319  15  User ancillary field 19
!  67  1  320  15  User ancillary field 20
!  68  1  127  12  NH3 (ammonia) aerosol emission
!  69  1  128  23  Surface fresh soot aerosol emission
!  70  1  129  23  Elevated fresh soot aerosol emission
!  71              Not used
!  72  1  121  17  Natural Sulphur dioxide emissions
!  73  1  122  18  OH concentrations
!  74  1  123  18  HO2 concentrations
!  75  1  124  18  H2O2 concentrations
!  76  1  125  18  Ozone (CHEM) concentrations
!  77  1  126  12  Sulphur dioxide high level emission
!  78  1  251  24  Surface CO2 emissions
!  79  1  207   4  Clapp-Hornberger parameter
!  80  1  208   5  Leaf area index of vegetated fraction
!  81  1  209   5  Canopy height of vegetated fraction
!  82  1  160  19  Aerosol data for radiative forcing of climate change
!  83  1  216  20  Initial fractions of surface types
!  84  1  217  21  Initial leaf area index of plant functional types
!  85  1  218  21  Initial canopy height of plant functional types
!  86  1  213  21  Initial gridbox mean canopy conductance
!  87  1  219  22  Fraction of vegetation subject to disturbance
!  88  1  220   4  Snow free albedo of bare soil
!  89  1  223   4  Soil carbon content
!  90  1  321  16  User ancillary multi 1
!  91  1  322  16  User ancillary multi 2
!  92  1  323  16  User ancillary multi 3
!  93  1  324  16  User ancillary multi 4
!  94  1  325  16  User ancillary multi 5
!  95  1  326  16  User ancillary multi 6
!  96  1  327  16  User ancillary multi 7
!  97  1  328  16  User ancillary multi 8
!  98  1  329  16  User ancillary multi 9
!  99  1  330  16  User ancillary multi 10
! 100  1  331  16  User ancillary multi 11
! 101  1  332  16  User ancillary multi 12
! 102  1  333  16  User ancillary multi 13
! 103  1  334  16  User ancillary multi 14
! 104  1  335  16  User ancillary multi 15
! 105  1  336  16  User ancillary multi 16
! 106  1  337  16  User ancillary multi 17
! 107  1  338  16  User ancillary multi 18
! 108  1  339  16  User ancillary multi 19
! 109  1  340  16  User ancillary multi 20
!  ------------------------------------------------------------------

CL  1.  Initialisation for atmosphere


      ICODE=0 
      IOUNIT=0
      UPDATE_MONTHS=0
      INCR_SEC = 0

!     Initialise ANCIL1/2. Includes Halos for MPP runs.
      DO I=1,P_FIELD
        ANCIL1(I)=0.0
        ANCIL2(I)=0.0
      ENDDO
CL  1.1 Set logical UPDATE for each ancillary field independently

      DO FIELD=1,NANCIL_FIELDS


        UPDATE(FIELD)=.FALSE.
        IF(STEPS(FIELD).NE.0) THEN
C         UPDATE(FIELD)=MOD(A_STEP,STEPS(FIELD)).EQ.0
          UPDATE(FIELD)=(MOD(A_STEP+OFFSET_STEPS,STEPS(FIELD)).EQ.0
     &                   .OR.A_STEP.EQ.0)
     &                    .AND.FIELDCODE(1,FIELD).GT.0
     &                     .AND.D1_ANCILADD(FIELD).GT.1
        END IF

CL  1.05 Copy ancillary updating reference time to local variables
      I_YEAR_REF   = ANCIL_REFTIME(1)
      I_MONTH_REF  = ANCIL_REFTIME(2)
      I_DAY_REF    = ANCIL_REFTIME(3)
      I_HOUR_REF   = ANCIL_REFTIME(4)
      I_MINUTE_REF = ANCIL_REFTIME(5)
      I_SECOND_REF = ANCIL_REFTIME(6)
CL       and convert to reference days & secs
            CALL TIME2SEC(I_YEAR_REF,I_MONTH_REF,I_DAY_REF,
     &                    I_HOUR_REF,I_MINUTE_REF,I_SECOND_REF,
     &                    0,0,ANCIL_REF_DAYS,ANCIL_REF_SECS,LCAL360)

C
      IF (.NOT. LCAL360) THEN

CL  1.11 Set logical UPDATE for Gregorian calender updates at monthly
CL       or yearly intervals. NB STEPS value set to 1 day in INANCILA
        IF(FIELDCODE(1,FIELD).EQ.1.OR.FIELDCODE(1,FIELD).EQ.2) THEN
          MONTHS=I_MONTH+I_YEAR*12-(I_MONTH_REF+I_YEAR_REF*12)
          UPDATE_MONTHS= FIELDCODE(2,FIELD)*
     &     ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          UPDATE(FIELD)=MOD(MONTHS,UPDATE_MONTHS).EQ.0.AND.I_DAY.EQ.1
        END IF
      END IF !  (.NOT.LCAL360)
C


      END DO

CL 1.2 Allow for dependencies between fields
C Sea surface temperature must be updated when sea ice is updated

      UPDATE(28)=UPDATE(27).OR.UPDATE(28)

C Both surface current components must be updated together

      UPDATE(30)=UPDATE(30).OR.UPDATE(31)
      UPDATE(31)=UPDATE(30)

CL Select method of time interpolation for SST. The interpolation
CL allows for sea ice if ice data is available at the same times
CL as the temperature data. Otherwise linear interpolation is used.

      LT_INT_C=.TRUE.

      IF(UPDATE(28)) THEN
      IF(FIXHD(10,FILEANCIL(27)).EQ.0) LT_INT_C=.FALSE.
        IF(LT_INT_C) THEN
        DO I=21,41
          IF(FIXHD(I,FILEANCIL(27)).NE.FIXHD(I,
     &      FILEANCIL(28))) THEN
            LT_INT_C=.FALSE.
            WRITE(6,*)' WARNING:controlled time interpolation for SST',
     &      ' not available: Mismatch in SST and SEA-ICE ancillary data'
     &     ,' times in FIXED HEADER'
            WRITE(6,*)' position=',I,' SEA-ICE=',FIXHD(I,FILEANCIL(27))
            WRITE(6,*)' position=',I,' SST=',FIXHD(I,FILEANCIL(28))
          END IF
        END DO
        ENDIF
      END IF



CL Loop over ancillary fields(atmosphere)

      DO FIELD=1,NANCIL_FIELDS

        LICE_DEPTH=field.eq.29  ! required for LAMIPII

      IF (UPDATE(FIELD)) THEN  ! (1st level IF)
        FILE=FILEANCIL(FIELD)
        NFTIN=FTNANCIL(FILE)

       IF(LICE_DEPTH.AND.LAMIPII) THEN

! Uses ice fraction set earlier in field loop.
! WARNING this will fail if the order of ancillary fields is ever 
! changed so that ice-depth preceeds ice fraction
! Note : For complete sea ice cover
!        Arctic ice depth    = 2m
!        Antarctic ice depth = 1m
! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
! This results in similar values to those from runs using ancillary 
! files containing ice depths set to 1 or 2m.  

          ROW_LENGTH=P_FIELD/P_ROWS
          DO I=1,P_ROWS
! work out latitude in radians
            LAT_P=FIRST_LAT-NS_SPACE*(I+datastart(2)-Offy-1)
            DO J=1,ROW_LENGTH
              II=J+(I-1)*ROW_LENGTH
              ANCIL_DATA(II)=0.0    
              IF (ICE_FRACTION(II).gt.0.0) THEN
                IF (LAT_P.GT.0.0) THEN   ! Arctic ice depth 
                  ANCIL_DATA(II)=2.*ICE_FRACTION(II)
                ELSE                     ! Antarctic ice depth
                  ANCIL_DATA(II)=1.*ICE_FRACTION(II)
                ENDIF  
              ENDIF
            ENDDO  
          ENDDO 
!L     Sea ice thickness
!L       Update over all sea points (all sea ice points are the only    
!L       ones strictly required, but this cannot be determined easily)
                                                                        
          DO I=1,P_FIELD                                                
            IF(.NOT.LAND(I)) THEN                                       
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)                  
            END IF                                                      
          END DO                                                        
       ELSE  
C     Update required for field

        WRITE(6,*)'REPLANCA: UPDATE REQUIRED FOR FIELD',FIELD

          IF ( FIXHD(10,FILE) .LT. 0 .OR. FIXHD(10,FILE) .GT. 2 ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in fixed header(10) of ancillary
     & file                           '
            RETURN
          ENDIF

CL    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE).EQ.0

CL    Set default values for time interpolation

        LINTERPOLATE=.TRUE.
        IF(SINGLE_TIME) THEN
          LINTERPOLATE=.FALSE.
        END IF

        IF (FIELD.GT.9 .AND. FIELD.LT.27) THEN
          LINTERPOLATE=.FALSE.
        END IF

CL 2.1 Find position of input record

CL    Default settings of search parameters if only one time present

        IF(SINGLE_TIME) THEN
          STEP=0
        ELSE


          LGREG_MONTHLY=.FALSE.
C
      IF (.NOT. LCAL360) THEN
          IF(FIELDCODE(1,FIELD).EQ.1.OR.FIELDCODE(1,FIELD).EQ.2) THEN
            LGREG_MONTHLY=.TRUE.
            UPDATE_MONTHS= FIELDCODE(2,FIELD)*
     &      ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          END IF
      END IF
C


          PERIODIC=FIXHD(10,FILE).EQ.2
          REGULAR=.TRUE.

C
      IF (.NOT. LCAL360) THEN
          REGULAR=FIXHD(35,FILE).EQ.0.AND.FIXHD(36,FILE).
     &    EQ.0
C i.e. data at intervals of days/hours & non-periodic
          IF(PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE).EQ.0
C i.e. data at intervals of hours & periodic
      END IF
C

C         Error checking on time information.

          IF ( FIXHD(35,FILE) .LT. 0 .OR.
     &         FIXHD(36,FILE) .LT. 0 .OR. FIXHD(36,FILE) .GT. 12 .OR.
     & REGULAR .AND. ( FIXHD(37,FILE) .LT. 0 .OR. FIXHD(37,FILE) .GT. 31
     &  .OR. FIXHD(38,FILE) .LT. 0 .OR. FIXHD(38,FILE) .GT. 24 ) ) THEN
C           FIXHD(39-40) are not used by REPLANCA.
C           FIXHD(35-37) have already been used if not CAL360.
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in validity time interval given
     &in ancillary file (FIXHD(35-38))'
            RETURN
          ENDIF

          IF ( FIXHD(21,FILE) .LT. 0 .AND. .NOT. PERIODIC
     &  .OR. .NOT. ( REGULAR .AND. PERIODIC ) .AND.
C    !  If it is REGULAR & PERIODIC more detailed check is applied below
     &     ( FIXHD(22,FILE) .LT. 0 .OR. FIXHD(22,FILE) .GT. 12 .OR.
     &       FIXHD(23,FILE) .LT. 0 .OR. FIXHD(23,FILE) .GT. 31 .OR.
     &       FIXHD(24,FILE) .LT. 0 .OR. FIXHD(24,FILE) .GT. 24 .OR.
     &       FIXHD(25,FILE) .LT. 0 .OR. FIXHD(25,FILE) .GT. 60 .OR.
     &       FIXHD(26,FILE) .LT. 0 .OR. FIXHD(26,FILE) .GT. 60 ) ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in first validity time given in
     & ancillary file (FIXHD(21-26))  '
            RETURN
          ENDIF

          IF(.NOT.PERIODIC) THEN

CL            If data taken from full time series of input data.

            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR
     &                    ,I_MINUTE,I_SECOND
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC
     &                    ,LCAL360)


CL Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (A_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR
     &                    ,I_MINUTE,I_SECOND
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2).EQ.0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (A_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


            IF(REGULAR) THEN
CL 2.1.1  Standard cases:360 day calender;
CL 2.1.1  or Gregorian calendar with
CL        interval between data times in days or hours
CL        updating interval may be regular in model timesteps,
CL        or (LGREG_MONTHLY=T) irregular in model timesteps,

              HOURS=SEC/3600+DAY*24
CL FInd time(in hours) of first ancillary data on file
              CALL TIME2SEC(FIXHD(21,FILE),FIXHD(22,FILE),
     &                   FIXHD(23,FILE),FIXHD(24,FILE),
     &                   FIXHD(25,FILE),FIXHD(26,FILE),
     &                   ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &                   LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

              IF(HOURS.LT.0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF

CL FInd interval(in hours) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

C Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL).EQ.0) THEN
                LINTERPOLATE=.FALSE.
              END IF

              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE

CL 2.1.2 Gregorian calender;ancillary data interval is in months or
CL       years,which is irregular in model timesteps.
!L original code is inaccurate for this section - corrected code under
!L LAMIPII makes use of dates in lookup headers 
!L For a real calendar year the mid-point of each month is different
!L in terms of its hour and day. The old inaccurate method assumes
!L the hour and day are taken from the fixhd values. These are only 
!L usually correct for the first month on the ancillary file. 


CL Adjust YMD time to middle of updating interval

              I_YEAR1=I_YEAR
              I_MONTH1=I_MONTH
              I_DAY1=I_DAY
              I_HOUR1=I_HOUR
              CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,
     &                     I_YEAR,I_MONTH,I_DAY,
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,
     &                     LCAL360)


CL FInd interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
C  Check for time within month
           IF (LAMIPII) THEN   ! corrected code uses pp header
              STEP=MONTHS/INTERVAL                              
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
! Check against day and hour of actual lookup header not first field 
              IF((I_DAY*24+I_HOUR).LT.
     &           (LOOKUP(3,I1)*24+LOOKUP(4,I1))) THEN   
                MONTHS=MONTHS-1                                         
              END IF                                                    
           ELSE              ! old less accurate code uses FIXHD
              IF((I_DAY*24+I_HOUR).LT.                                  
     *           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
           ENDIF ! LAMIPII  

              IF(MONTHS.LT.0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF


CL Adjust YMD time back to start of updating interval

              I_YEAR=I_YEAR1
              I_MONTH=I_MONTH1
              I_DAY=I_DAY1
              I_HOUR=I_HOUR1



              STEP=MONTHS/INTERVAL

           IF (LAMIPII) THEN       ! corrected code
              TIME=REAL(SEC)/3600+REAL(DAY*24)   
! correct calculation of dates uses lookup table dates not fixhd date
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
              I_YEAR1=lookup(1,i1)
              I_MONTH1=lookup(2,i1)     
              I_DAY1=lookup(3,i1)
              I_HOUR1=lookup(4,i1)
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,
     &              FIXHD(25,FILE),FIXHD(26,FILE),                 
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,         
     &              LCAL360)                                       
              TIME1=REAL(SEC)/3600+REAL(DAY*24) 
! I1+1 correct pointer to next field as only one field in ancil file    
              I_YEAR1=lookup(1,i1+1)
              I_MONTH1=lookup(2,i1+1)     
              I_DAY1=lookup(3,i1+1)
              I_HOUR1=lookup(4,i1+1)     
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,
     &              FIXHD(25,FILE),FIXHD(26,FILE),                 
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,         
     &              LCAL360)                                       
              TIME2=REAL(SEC)/3600+REAL(DAY*24)                         

           ELSE   ! LAMIPII test - old inaccurate code using FIXHD
C NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
C Calculate data times for time interpolation
              TIME=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+FIXHD(22,FILE)-1)/12
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),
     &              FIXHD(25,FILE),FIXHD(26,FILE),
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+INTERVAL+FIXHD(22,FILE)-1)/12
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),
     &              FIXHD(25,FILE),FIXHD(26,FILE),
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF     ! end LAMIPII test

C Do not interpolate in time if data time exactly matches model time

              IF(TIME.EQ.TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF ! End of REGULAR/not REGULAR

          ELSE  ! PERIODIC data

CL 2.2   If data is taken from ancillary periodic data.

            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR,
     &                     I_MINUTE,I_SECOND,
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &                     LCAL360)


CL Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (A_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR
     &                    ,I_MINUTE,I_SECOND
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2).EQ.0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (A_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


CL Adjust YMD time to middle of updating interval

            I_YEAR1=I_YEAR
            I_MONTH1=I_MONTH
            I_DAY1=I_DAY
            I_HOUR1=I_HOUR
            CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,
     &                     I_YEAR,I_MONTH,I_DAY,
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,
     &                     LCAL360)



            IF (REGULAR) THEN
CL 2.2.1 Standard cases:1) 360 day calender, with allowed periods of
CL       1 day, 1 month or 1 year;
CL
CL       2) Gregorian calender with update in hours,and period of
CL       data 1 day.
CL
CL       For both updating interval and number of
CL       data times to be skipped in data set calculated in hours.

              HOURS=SEC/3600+DAY*24
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

              PERIOD=INTHD(3,FILE)*INTERVAL

CL   Do not allow non-standard periods
      IF (LCAL360) THEN
              IF(PERIOD.NE.8640.AND.PERIOD.NE.720.AND.PERIOD.NE.24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ELSE
              IF(PERIOD.NE.24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ENDIF
              IF(PERIOD.EQ.24)THEN
C Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                IF(I_HOUR.LT.FIXHD(24,FILE)) HOURS=HOURS+24

              ELSE IF(PERIOD.EQ.720)THEN
C Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                IF((I_DAY*24+I_HOUR).LT.
     &             (FIXHD(23,FILE)*24+FIXHD(24,FILE)))
     &           HOURS=HOURS+720

              ELSE IF(PERIOD.EQ.8640)THEN
C Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                IF((I_MONTH*720+I_DAY*24+I_HOUR).LT.
     &          (FIXHD(22,FILE)*720+FIXHD(23,FILE)*24+FIXHD(24,FILE)))
     &           HOURS=HOURS+8640

              END IF

              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),
     &                     FIXHD(25,FILE),FIXHD(26,FILE),
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &                     LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

C Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL).EQ.0) THEN
                LINTERPOLATE=.FALSE.
              END IF
              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE  ! non regular case

CL 2.2.2 Gregorian calender,and data interval is in months,
CL       period is 1 year
CL       Updating interval and number of data times to be skipped
CL       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
              PERIOD=INTHD(3,FILE)*INTERVAL
              IF(PERIOD.NE.12)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
!  Difference between date now (month) & first date ancil file (month)
              MONTHS=I_MONTH-FIXHD(22,FILE)


           IF (LAMIPII) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which 
! contains values for first field on ancillary file only.
             step=months/INTERVAL    
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*step
             I1=I2+LOOKUP_START(FILE)-1
!  Check for time within month - using ppheader information             
          IF((I_DAY*24+I_HOUR).LT.(lookup(3,i1)*24+lookup(4,i1))) THEN 
                MONTHS=MONTHS-1    
          END IF                   
             IF(MONTHS.LT.0) THEN 
                MONTHS=MONTHS+12   
             END IF               
! recalculate STEP
              STEP=MONTHS/INTERVAL                         
! NB INTERVAL may be > 1 month                             
              MONTHS=STEP*INTERVAL                         
              IY=I_YEAR                                                 
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1               
              IF(IM.GT.I_MONTH) IY=IY-1                          
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
              CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),
     &              FIXHD(25,FILE),FIXHD(26,FILE),                 
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)      
              TIME1=REAL(SEC)/3600+REAL(DAY*24)                         
!  Calculate  TIME2 for second ancillary data time          
!  set IY correctly for time interpolation calculations     
              IY=I_YEAR                                     
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1 
              IF(IM.LT.I_MONTH) IY=IY+1                     
              I1=(IM-1)/INTERVAL                         
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*I1
              I1=I2+LOOKUP_START(FILE)-1
              CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),
     &              FIXHD(25,FILE),FIXHD(26,FILE),                 
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)      
              TIME2=REAL(SEC)/3600+REAL(DAY*24)                         

           ELSE   ! original code inaccurate use of FIXHD dates
C  Check for time within month
              IF((I_DAY*24+I_HOUR).LT.
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
              IF(MONTHS.LT.0) THEN
                MONTHS=MONTHS+12
              END IF

              STEP=MONTHS/INTERVAL
C NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
C  Calculate TIME1 for first ancillary data time
C  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IF(IM.GT.I_MONTH) IY=IY-1
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),
     &              FIXHD(25,FILE),FIXHD(26,FILE),
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
C  Calculate  TIME2 for second ancillary data time
C  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IF(IM.LT.I_MONTH) IY=IY+1
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),
     &              FIXHD(25,FILE),FIXHD(26,FILE),
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF  ! end LAMIPII test

C Do not interpolate in time if data time exactly matches model time

              IF(TIME.EQ.TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF  ! regular/non-regular


CL Adjust YMD time back to start of updating interval

            I_YEAR=I_YEAR1
            I_MONTH=I_MONTH1
            I_DAY=I_DAY1
            I_HOUR=I_HOUR1


          ENDIF  ! non-periodic/periodic

        IF (LINTERPOLATE) THEN
        WRITE(6,*)' REPLANCA - time interpolation for field ',field
        WRITE(6,*)' time,time1,time2 ',time,time1,time2
        WRITE(6,*)' hours,int,period ',hours,interval,period
        END IF

        END IF ! singletime/non-singletime

CL 2.3   Check STASH Code

        I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP

        I1=LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)-1)

        LMISMATCH=.FALSE.
        WRITE(6,*)' Information used in checking ancillary data set:',
     *  ' position of lookup table in dataset:',I2
        WRITE(6,*)' Position of first lookup table referring to ',
     *  'data type ',NLOOKUP(FIELD)
        WRITE(6,*)' Interval between lookup tables referring to data ',
     *  'type ', LOOKUP_STEP(FIELD),' Number of steps', STEP
        WRITE(6,*)' STASH code in dataset ',I1,
     *  '  STASH code requested ',STASHANCIL(FIELD)
        WRITE(6,*)'''Start'' position of lookup tables for dataset ',
     *  'in overall lookup array ' ,LOOKUP_START(FILE)

        IF(I1.NE.STASHANCIL(FIELD)) THEN
        WRITE(6,*)I1,STASHANCIL(FIELD),FIELD
          LMISMATCH=.TRUE.
        END IF

CL Error exit if checks fail

        IF(LMISMATCH) THEN
          ICODE=200+FIELD
         CMESSAGE='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
         RETURN
        END IF

        IF(LINTERPOLATE.AND..NOT.SINGLE_TIME) THEN
CL Check time interpolation factors
          IF(TIME.LT.TIME1.OR.TIME.GT.TIME2) THEN
           WRITE(6,*)' Information used in interpolation/replacement:'
           WRITE(6,*)' Time of first data=', TIME1
           WRITE(6,*)' Validity Time for update=', TIME
           WRITE(6,*)' Time of second data=', TIME2

           ICODE=500+FIELD
           CMESSAGE='REPLANCA: TIME INTERPOLATION ERROR'
           RETURN
          END IF
        END IF

CL 3   Loop over levels of ancillary data for field I
CL Reset pointer for dataset


CL Includes loop over X and Y components of surface currents

         LICE_FRACTION=FIELD.EQ.27
         LSNOW_DEPTH=FIELD.EQ.9
         LICE_DEPTH=FIELD.EQ.29

        DO 30 LEVEL=1,LEVELS(FIELD)

CL Do not go through loop for ice edge or snow edge

        IF((LICE_FRACTION.OR.LSNOW_DEPTH).AND.LEVEL.EQ.2) THEN
          GOTO 30
        END IF

CL 3.1 Read data for single level of ancillary field.

        IF(.NOT.LICE_FRACTION) THEN
! AMIPII case ice depth field not read from ancillary file
         IF(.NOT.(LICE_DEPTH.and.LAMIPII)) THEN
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),
     &                  LEN1_LOOKUP,ANCIL1,P_FIELD,FIXHD(1,FILE),
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
     &                  ICODE,CMESSAGE)

         ENDIF
          IF(ICODE.NE.0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

        ELSE

CL If ice-fraction,read fractional time field as well
CL       UNLESS IT IS A SINGLE TIME FIELD
CL If snow-depth,read fractional time field as well only if time
CL interpolation required.

      IF(.NOT.SINGLE_TIME.and..NOT.LAMIPII) THEN                        
         IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)).EQ.38) THEN
          CALL READFLDS(NFTIN,2,I2,LOOKUP(1,LOOKUP_START(FILE)),
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),
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
     &                  ICODE,CMESSAGE)
          IF(ICODE.NE.0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE CHANGE DATA MISSING'
            RETURN
         END IF
        ELSE    ! single time or LAMIPII - ie no time change field
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),
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
     &                  ICODE,CMESSAGE)
          IF(ICODE.NE.0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA: I/O ERROR'
            RETURN
          ENDIF
        END IF
      ENDIF

        IF(LSNOW_DEPTH.AND.LINTERPOLATE) THEN
      IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)).EQ.27) THEN

           CALL READFLDS(NFTIN,1,I2+1,LOOKUP(1,LOOKUP_START(FILE)),
     &                   LEN1_LOOKUP,SNOW_CHANGE,P_FIELD,FIXHD(1,FILE),
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
          IF(ICODE.NE.0)THEN
             ICODE=FIELD+100
             IOUNIT=NFTIN
            CMESSAGE='REPLANCA :I/O ERROR '
             RETURN
           END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
           CMESSAGE='REPLANCA :SNOW CHANGE DATA MISSING'
           RETURN
         END IF
        END IF

CL If sea surface temperature or other ice fields, read ice fraction
CL and fractional time field if not already pressent and if required
CL by time interpolation.  Similar if SLAB ref SST or ice depth needed.

        IF(FIELD.EQ.29.OR.(FIELD.EQ.28.AND.LT_INT_C).OR.
     &     FIELD.EQ.38.OR.(FIELD.EQ.37.AND.LT_INT_C))
     &    THEN

         IF(.NOT.UPDATE(27)) THEN
          I3 = NLOOKUP(27) + LOOKUP_STEP(27)*STEP + LOOKUP_START(
     &       FILEANCIL(27))
          IF ( LOOKUP(ITEM_CODE,I3) .EQ. 38 ) THEN

            CALL READFLDS(FTNANCIL(FILEANCIL(27)),2,
     &                    NLOOKUP(27)+LOOKUP_STEP(27)*STEP,
     &                    LOOKUP(1,LOOKUP_START(FILEANCIL(27))),
     &                    LEN1_LOOKUP,ICE_EXTENT,
     &                    P_FIELD,FIXHD(1,FILEANCIL(27)),
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
          IF(ICODE.NE.0)THEN
              ICODE=FIELD+100
              IOUNIT=NFTIN
             CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF
          IF ( RLOOKUP(BMDI,I3-1) .NE. RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti
     &mes of sea-ice chge not standard'
            RETURN
          ENDIF


          ELSE
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE FIELD DATA MISSING'
            RETURN
          END IF
         END IF
        END IF

CL 3.3 If time interpolation required, read second record

        IF(LINTERPOLATE) THEN

          I1=I2+ LOOKUP_STEP(FIELD)
          IF(I1.LE.FIXHD(152,FILE)) THEN

! AMIP II and ice depth don't read in ice depth field
          IF (.NOT.(LAMIPII.and.LICE_DEPTH)) THEN

            CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),
     &                    LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),
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
          ENDIF
          IF(ICODE.NE.0)THEN
              ICODE=FIELD+300
              IOUNIT=NFTIN
              CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF

          ELSE !end of data on file

CL  If end of data has been reached go back to the start.If data is
CL  periodic.
CL  Otherwise cancel time interpolation

            IF(PERIODIC) THEN

              I1 = NLOOKUP(FIELD) + LEVEL - 1

              CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),
     &                      LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),
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
          IF(ICODE.NE.0)THEN
                ICODE=FIELD+300
                IOUNIT=NFTIN
               CMESSAGE='REPLANCA :I/O ERROR '
                RETURN
              END IF
            ELSE
              LINTERPOLATE=.FALSE.
            END IF
          END IF! End of position on file test

          ICODE=0
        END IF ! End LINTERPOLATE

CL 3.4 Perform time interpolation

        IF(LINTERPOLATE) THEN

          ZERO=0.0

CL Select appropriate time interpolation for each field
C  Snowdepth: set equal to zero if no snow cover

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              PRES_VALUE(I)=ZERO
            END DO

C For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
C  which was read in from position I2+1.
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I2) .NE. RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti
     &mes of snow change non-standard '
            RETURN
          ENDIF

            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &           TIME,P_FIELD,SNOW_CHANGE,ANCIL1,PRES_VALUE)

C Ice fraction: ice depth set equal to zero if no ice

          ELSE IF(FIELD.EQ.27.OR.FIELD.EQ.29.OR.FIELD.EQ.38) THEN
            IF(FIELD.EQ.27) THEN
C For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
C  which was read in from position I1+1
          IF(.NOT.LAMIPII) THEN
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1) .NE. RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti
     &mes of sea-ice chge non-standard'
            RETURN
          ENDIF
          ENDIF    

             IF (LAMIPII) THEN
! linear uncontrolled time interpolation
              CALL T_INT (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA, 
     &             TIME,P_FIELD)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary 
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

              DO I=1,P_FIELD
                IF (ANCIL_DATA(I).lt.0.3) ANCIL_DATA(I)=0.0
                IF (ANCIL_DATA(I).gt.1.0) ANCIL_DATA(I)=1.0
              ENDDO

             ELSE       ! non AMIPII option
              DO I=1,P_FIELD                                            
                PRES_VALUE(I)=0
              END DO
               
              CALL T_INT_C (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

             ENDIF     ! end AMIPII test

            ELSE IF (FIELD.EQ.29.OR.FIELD.EQ.38) THEN

              DO I=1,P_FIELD
                PRES_VALUE(I)=0
              END DO

              CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)


            END IF


C Sea surface temperature, set equal to TFS if ice present

          ELSE IF ((FIELD.EQ.28.OR.FIELD.EQ.37).AND.LT_INT_C) THEN
           IF (LAMIPII) THEN

            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,          
     &              TIME,P_FIELD)
! remove any T below TFS
            DO I=1,P_FIELD
              IF (ANCIL_DATA(i).LT.TFS)  ANCIL_DATA(I)=TFS              
            ENDDO                                                

           ELSE     ! non AMIPII option
                                                                        
            DO I=1,P_FIELD                                              
                PRES_VALUE(I)=TFS

C Set no_ice_extent indicator for controlled SST interpolation
                IF(ICE_EXTENT(I,1).EQ.0) THEN
                  NO_ICE_EXTENT(I)=1.0
                ELSE
                  NO_ICE_EXTENT(I)=0.0
                ENDIF
            END DO

            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &           TIME,P_FIELD,ICE_EXTENT(1,2),NO_ICE_EXTENT,PRES_VALUE)

           ENDIF   ! end AMIPII test
C Otherwise linear interpolation in time, unless missing data indicator
C present at either time.

          ELSE

C Time interpolation checks the data against the standard missing data
C   indicator - check that the field is labelled as using the same one.
C  (It is to have the right I1 here that I3 is used above.)
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1) .NE. RMDI .OR.
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1) .NE. RMDI ) THEN
            WRITE (6, *) 'LOOKUPS:',
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1),
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Missing data indicator in lookup of an
     &cillary file is non-standard    '
            RETURN
          ENDIF

          LEN=P_FIELD
CL  Ozone, test for zonal mean or full field
          IF(FIELD.EQ.7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1).EQ.1) THEN
              LEN=P_ROWS
            END IF
          END IF

            CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &                 TIME,LEN)

          END IF ! End Lsnow_depth

C If no interpolation, copy data into final array

        ELSE ! no interpolation
         IF(LICE_FRACTION) THEN
          IF (LAMIPII) THEN
          DO I=1,P_FIELD
             
          ANCIL_DATA(I)=ICE_EXTENT(I,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary 
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

             IF (ANCIL_DATA(I).lt.0.3) ANCIL_DATA(I)=0.0
             IF (ANCIL_DATA(I).gt.1.0) ANCIL_DATA(I)=1.0
      
          ENDDO
          ELSE           ! non AMIP II option
            DO I=1,P_FIELD   
             ANCIL_DATA(I)=ICE_EXTENT(I,1)
            ENDDO
          ENDIF           ! end of AMIPII test 
         ELSE IF (LAMIPII.AND.FIELD.EQ.28) THEN
          DO I=1,P_FIELD      
            ANCIL_DATA(I)=ANCIL1(I)    
            IF (ANCIL_DATA(I).lt.TFS) ANCIL_DATA(I)=TFS 
          ENDDO 
         ELSE
          DO I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
          END DO
         ENDIF
        END IF !End interpolate/no interpolate

CL 3.5 Updating action for each field at each level
CL     Fields replaced except that Sea Surface Temperature may be
CL     incremented. Take apropriate action for each field.

        IF(FIELD.LE.2.OR.FIELD.EQ.7.OR.FIELD.EQ.39.OR.FIELD.EQ.40.
     &  OR.FIELD.EQ.41.OR.FIELD.EQ.42.OR.FIELD.EQ.43.
     &  OR.FIELD.EQ.44.OR.FIELD.EQ.45.   ! multi-level murk
     &  OR.(FIELD.GE.68.AND.FIELD.LE.70). !NH3,soot aerosol emissions
     &  OR.(FIELD.GE.72.AND.FIELD.LE.77). !Sulphur cycle
     &  OR.FIELD.EQ.78.                   !CO2 EMISSIONS
     &  OR.FIELD.EQ.82.                   !HADCM2 sulphate aerosol
     &  OR.(FIELD.GE.90.AND.FIELD.LE.109) !multi-level user ancillaries
     &  )THEN

CL 3.5.0 Updates at all points

          LEN=P_FIELD
CL  Ozone, test for zonal mean or full field
          IF(FIELD.EQ.7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1).EQ.1) THEN
              LEN=P_ROWS
            END IF
          END IF

          DO I=1,LEN
            D1(D1_ANCILADD(FIELD)+I-1+(LEVEL-1)*LEN)=ANCIL_DATA(I)
          END DO

CL 3.5.1 Updates over all land points

        ELSEIF((FIELD.GT.2.AND.FIELD.LT.7).
     &   OR.(FIELD.GT.7.AND.FIELD.LT.27).
     &   OR.(FIELD.EQ.32).OR.(FIELD.GE.34.AND.FIELD.LE.36).OR.
     &  (FIELD.GE.48.AND.FIELD.LE.67). ! single level user ancillaries
     &  OR.(FIELD.GE.46.AND.FIELD.LE.47).      !Orographic roughness
     &  OR.(FIELD.GE.79.AND.FIELD.LE.81).      !MOSES-I
     &  OR.(FIELD.GE.83.AND.FIELD.LE.89)) THEN !MOSES-II


CL If not reconfiguration, set snowdepth values at all land points
CL Reset TSTAR to TM if snow cover present

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              IF(LAND(I)) THEN
                D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
                IF(TSTAR(I).GT.TM.AND.ANCIL_DATA(I).GT.0.0) THEN
                  TSTAR(I)=TM
                END IF
              END IF
            END DO

CL Set all other fields , which are stored at land points only

          ELSE
            CALL TO_LAND_POINTS(ANCIL_DATA,D1(D1_ANCILADD(FIELD)+
     &        (LEVEL-1)*LAND_FIELD),LAND,P_FIELD,I)
          END IF



CL 3.5.2 Ice fraction

        ELSE IF(FIELD.EQ.27) THEN
          DO I=1,P_FIELD
            ICE_FRACTION(I)=0.
            IF (.NOT.LAND(I)) THEN
              ICE_FRACTION(I)=ANCIL_DATA(I)
            END IF
          END DO

CL Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

          DO I=1,P_FIELD
            IF(ICE_FRACTION(I).GT.0.0) THEN
              TSTAR(I)=AMIN1(TSTAR(I),TFS)
            ENDIF
          END DO

CL 3.5.3 Sea surface temperatures for atmosphere, allow fields to be
CL       incremented rather than replaced

        ELSE IF (FIELD.EQ.28) THEN


          DO I=1,P_FIELD
            IF (.NOT.LAND(I).AND.ICE_FRACTION(I).EQ.0.0) THEN
              IF(L_SSTANOM) THEN
                TSTAR(I)=ANCIL_DATA(I)+TSTAR_ANOM(I)
              ELSE
                TSTAR(I)=ANCIL_DATA(I)
              END IF
            END IF
          END DO

CL 3.5.3.1 Reference SSTs for SLAB model

        ELSE IF (FIELD.EQ.37) THEN

          DO I=1,P_FIELD
            IF (.NOT.LAND(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1)=TSTAR(I)
            ENDIF
          END DO

CL 3.5.4 Sea ice thickness/Reference seaice thickness for SLAB
CL       Update over all sea points (all sea ice points are the only
CL       ones strictly required, but this cannot be determined easily)

        ELSE IF (FIELD.EQ.29.OR.FIELD.EQ.38) THEN

          DO I=1,P_FIELD
            IF(.NOT.LAND(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            END IF
          END DO

CL 3.5.5 Surface currents

        ELSE IF (FIELD.EQ.30.OR.FIELD.EQ.31) THEN
          DO I=1,U_FIELD
            IF(.NOT.LAND(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1) = 0.0
            END IF
          END DO

CL 3.5.6 Heat convergence (slab model)
CL       Update over all non-land points

        ELSE IF (FIELD.EQ.33) THEN

          DO I=1,P_FIELD
            IF(.NOT.LAND(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1) = 0.0
            END IF
          END DO

        ELSE

        WRITE(6,*)' REPLANCA: ERROR - FIELD ',FIELD,
     &  ' omitted from update block'

        END IF !End tests on FIELD numbers

CL End loop over levels

      I2=I2+1

 30   CONTINUE

CL End loop over ancillary fields (atmosphere)
       ENDIF ! LAMIPII and ice depth test

      END IF    ! End UPDATE(field) test     level 1 IF


      END DO

900   RETURN
      END

