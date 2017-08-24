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
CLL SUBROUTINE REPLANCW
CLL
CLL Purpose:  Updates ancillary fields as requested in FIELDCODE array.
CLL   Tests whether update is required for each field, allowing for
CLL   dependencies between fields. Uses LOOKUP array to find data for
CLL   appropriate time, reads a record and checks for current data
CLL   type. Reads second record if time interpolation required. Updates
CLL   the field. Under DEF RECON, the interface to the routine is
CLL   modified for use in the reconfiguration rather than the model.
CLL   At present there is no reconfiguration for the wave sub-model.
CLL   Under DEF CAL360 the 360 day rather than the Gregorian calender
CLL   is used.
CLL
CLL  Model            Modification history
CLL version  Date
CLL  4.1  08/05/96  New routine for wave sub-model.  RTHBarnes.
CLL
CLL PROGRAMMING STANDARD: UMDP NO 3 VERSION NO 2, DATED 07/09/90
CLL
CLL SYSTEM COMPONENTS COVERED: C71
CLL
CLL SYSTEM TASK C7
CLL
CLL  External documentation: UMDP C7
CLL
CLLEND

      SUBROUTINE REPLANCW
     &(I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,
     & ANCIL_REFTIME,OFFSET_STEPS,IMT,JMT,D1,
CCC *IF -DEF,RECON
     & W_STEP,W_STEPS_P_P,W_SECS_P_P,
CCC *ENDIF
     & LEN1_LOOKUP,
     & LEN_FIXHD,
     & LEN_INTHD,
     & LEN_REALHD,
     & LEN_D1,
     & FIXHD,
     & INTHD,
     & REALHD,
     & LOOKUP,
     & RLOOKUP,
     & FTNANCIL,
     & LOOKUP_START,
     & NDATASETS,
     & NLOOKUPS,
CCC *IF DEF,RECON
CCC      & IOUNIT,
CCC *ENDIF
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
     & ICODE,CMESSAGE,LCAL360)


      IMPLICIT NONE

      LOGICAL LCAL360
C Include COMDECKS
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
C*L--------------------COMDECK  CANCILO ---------------------------
C
C Contains PARAMETERS, headers, and index blocks for control of update
C of ancillary fields.
C System component F0171 Parameters determined by model version and
C size parameters.
! History:
! Version  Date     Comment
! -------  ----     -------
!  4.1  08/05/96  New deck for wave sub-model.  RTHBarnes
C
      INTEGER
     &  NANCIL_FIELDS,   ! Maximum total number of ancillary fields
     &  FILEANCIL,       ! File number associated with ancillary fields
     &  NLOOKUP,         ! Position of given ancillary field in lookup
C                        ! tables
     &  LOOKUP_STEP,     ! Interval between PP Headers refering to
C                        ! to the same ancillary fields at diferent time
     &  LEVELS,          ! Number of levels of data in each ancillary
C                        ! field (Set by INANCILW )
     &  STASHANCIL,      ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block

C PARAMETER statement, fixed by model version

      PARAMETER
     &          (NANCIL_FIELDS=2) ! Hard-wired

C*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL
     &         UPDATE
      INTEGER  FIELDCODE,
     &         STEPS
C*----------------------------------------------------------------------
      COMMON/CTANCILW/
     &         FIELDCODE(2,NANCIL_FIELDS),
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS)
      COMMON/IXANCILW/ FILEANCIL(NANCIL_FIELDS),
     &           NLOOKUP(NANCIL_FIELDS),
     &           LOOKUP_STEP(NANCIL_FIELDS),
     &           LEVELS(NANCIL_FIELDS),
     &           STASHANCIL(NANCIL_FIELDS),
     &           D1_ANCILADD(NANCIL_FIELDS)

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


      INTEGER
     &       I_YEAR,            ! Curent Model Time
     &       I_MONTH,           !   "      "     "
     &       I_DAY,             !   "      "     "
     &       I_HOUR,            !   "      "     "
     &       I_MINUTE,          !   "      "     "
     &       I_SECOND,          !   "      "     "
     &       I_DAY_NUMBER,      !
     &       ANCIL_REFTIME(6),  ! Reference time for ancillary updating
     &       OFFSET_STEPS,      ! Offset in timesteps of ref. from basis
CCC *IF -DEF,RECON
     &       W_STEP,            !
     &       W_STEPS_P_P,       ! steps per period
     &       W_SECS_P_P,        ! seconds per period
CCC *ENDIF
CCC  &       BASIS_TIME_DAYS,   ! Model basis time in whole days
CCC  &       BASIS_TIME_SECS,   ! Model basis time in extra seconds
     &       NDATASETS,         ! Number of ancillary datasets
     &       NLOOKUPS,          ! Number of lookup tables
     &       LEN_D1             ! Size of primary data array
     &      ,IMT                ! Zonal dimension of arrays
     &      ,JMT                ! Meridional dimension of arrays

      INTEGER
     &       LEN1_LOOKUP,       ! First dimension of lookup table
     &       LEN_FIXHD,         ! Length of headers in data sets
     &       LEN_INTHD,
     &       LEN_REALHD,
     &       FIXHD(LEN_FIXHD,NDATASETS),  ! Data set headers
     &       INTHD(LEN_INTHD,NDATASETS),  !
     &       FTNANCIL(NDATASETS),         ! FTN numbers of data sets
     &       LOOKUP_START(NDATASETS),     ! Start of lookup tables
C                                         ! referring to data set
     &       LOOKUP(LEN1_LOOKUP,NLOOKUPS) ! Data set lookup tables

      REAL
     &       D1(LEN_D1),        ! Primary data array
     &       REALHD(LEN_REALHD,NDATASETS),
     &       RLOOKUP(LEN1_LOOKUP,NLOOKUPS)

      INTEGER
     &       I_AO,     ! Sub-model indicator = 1 Atmosphere
C                                = 2 Ocean   = 4 Wave
     &       ICODE,    ! Return code
     &       IOUNIT    ! OUT  I/O unit passed out in recon mode

      CHARACTER*(80)
     &       CMESSAGE  ! Error message


C*L   Subroutines called;
      EXTERNAL
     &       TIME2SEC,
     &       READFLDS,
CCC *IF -DEF,RECON
     &       SEC2TIME,TIME_DF,
CCC *ENDIF
     &       T_INT

C*L   Local integer arrays
      REAL
     &       ANCIL1(IMT*JMT),     ! Buffers to hold values of ancillary
C                                 ! data for time interpolation.
     &       ANCIL2(IMT*JMT),     !
     &       ANCIL_DATA(IMT*JMT)  ! Field of ancillary data held prior
C                                 ! to selective updating.

C     Local variables
      INTEGER
     &       I,                 !
     &       JADDR,
     &       J,                 !
     &       I1,       ! used for checking stash item code in sec 2.3
     &       I2,       ! ptr to 1st field; calculated in sec 2.3
     &       I1LEV,    ! ptr to 2nd field for this level; used in sec 3
     &       I2LEV,    ! ptr to 1st field for this level; used in sec 3
     &       ID,                !
     &       IM,                !
     &       IY,                !
     &       FIELD,             ! Current field number.
     &       FILE               !
      INTEGER
     &       INTERVAL,          ! Interval between data times.
     &       STEP,              ! Number of data times skipped
     &       MONTHS,            ! Used in calculations of position of
C                               ! data required
     &       HOURS,
     &       DAYS,SECONDS,      ! Times used in intermediate calculation
     &       PERIOD,
     &       START_MONTH,       !
     &       NFTIN,             ! Current FTN number for ancillary field
     &       ANCIL_REF_DAYS,    ! Ancil.reference time in whole days
     &       ANCIL_REF_SECS,    ! Ancil.reference time in extra seconds
     &       DAY,SEC,           ! Times relative to basis time
     &       DAY1,SEC1,         ! Times relative to basis time
     &       INCR_SEC,          ! Increment in sec
     &       LEN_IO
     &      ,LEVEL              ! loop index for level number
     &      ,LEN_FLD            ! length of (single level) field
     &      ,LEN_FLD_ACC        ! accumulated length of fields on
C                                 previous levels
     &      ,POS_STRT           ! start position in D1 array

CCC *IF -DEF,RECON
      INTEGER
     &       I_YEAR1,            ! Copy of Curent Model Time year
     &       I_MONTH1,           !   "      "     "          month
     &       I_DAY1,             !   "      "     "          day
     &       I_HOUR1,            !   "      "     "          hour
     &       I_MINUTE1,          !   "      "     "          minute
     &       I_SECOND1           !   "      "     "          second

      INTEGER
     &       UPDATE_MONTHS      ! update frequency (months) if Gregorian
      LOGICAL
     &       LGREG_MONTHLY      ! True for Gregorian monthly updating
      INTEGER
     &       I_YEAR_BASIS,            ! Basis Model Time
     &       I_MONTH_BASIS,           !   "     "     "
     &       I_DAY_BASIS,             !   "     "     "
     &       I_HOUR_BASIS,            !   "     "     "
     &       I_MINUTE_BASIS,          !   "     "     "
     &       I_SECOND_BASIS,          !   "     "     "
     &       I_DAY_NUMBER_BASIS
      INTEGER
     &       I_YEAR_REF,              ! Reference Time
     &       I_MONTH_REF,             !    "       "
     &       I_DAY_REF,               !    "       "
     &       I_HOUR_REF,              !    "       "
     &       I_MINUTE_REF,            !    "       "
     &       I_SECOND_REF             !    "       "
CCC *ENDIF


      LOGICAL
     &       LINTERPOLATE,      ! Indicates whether time
C                               ! interpolation needed.
     &       LMISMATCH,         ! Used in header chacks
     &       LICE_DEPTH,        ! Number of data times skipped
     &       SINGLE_TIME,       ! Indicates that only one time is
C                               ! available in data set
     &       PERIODIC,          ! Data set is periodic
     &       REGULAR            ! Interval between data times in
C                               ! dataset is regular in model timesteps.

      REAL
     &       A_IO               ! Used in check for I/O errors
     &,      TIME               ! Target time for time interpolation
     &,      TIME1              ! Times if data used in time interpn.
     &,      TIME2              !   "


CL  1.  Initialisation for ocean

C     Ocean fields
C     FIELD=  1 u-component wind
C             2 v-component wind

      IOUNIT=0
      INCR_SEC = 0
CL  1.1 Set logical switches for each ancillary field independently

      DO FIELD=1,NANCIL_FIELDS

CCC *IF -DEF,RECON

        UPDATE(FIELD)=.FALSE.
        IF(STEPS(FIELD).NE.0) THEN
C         UPDATE(FIELD)=MOD(W_STEP,STEPS(FIELD)).EQ.0
          UPDATE(FIELD)=(MOD(W_STEP+OFFSET_STEPS,STEPS(FIELD)).EQ.0
     &                   .OR.W_STEP.EQ.0)
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

      IF (.NOT. LCAL360) THEN

CL  1.11 Set logical UPDATE for Gregorian calender updates at monthly
CL       or yearly intervals. NB STEPS value set to 1 day in INANCILO
        IF(FIELDCODE(1,FIELD).EQ.1.OR.FIELDCODE(1,FIELD).EQ.2) THEN
          MONTHS=I_MONTH+I_YEAR*12-(I_MONTH_REF+I_YEAR_REF*12)
          UPDATE_MONTHS= FIELDCODE(2,FIELD)*
     &     ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          UPDATE(FIELD)=MOD(MONTHS,UPDATE_MONTHS).EQ.0.AND.I_DAY.EQ.1
        END IF

      END IF
CCC *ELSE
CCC         UPDATE(FIELD)=FIELDCODE(FIELD).GT.0
CCC *ENDIF

      END DO


CL Loop over ancillary fields (wave)

      DO FIELD=1,NANCIL_FIELDS


      IF (UPDATE(FIELD)) THEN  ! (1st level IF block)
        FILE=FILEANCIL(FIELD)
        NFTIN=FTNANCIL(FILE)

C     Update required for field

        WRITE(6,*)'REPLANCW: UPDATE REQUIRED FOR FIELD',FIELD

CL    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE).EQ.0

CL    Set default values for time interpolation

        LINTERPOLATE=.TRUE.
        IF(SINGLE_TIME) THEN
          LINTERPOLATE=.FALSE.
        END IF


CL 2.1 Find position of input record

CL    Default settings of search parameters if only one time present

        IF(SINGLE_TIME) THEN
          STEP=0
        ELSE

CCC *IF -DEF,RECON

          UPDATE_MONTHS=0
          LGREG_MONTHLY=.FALSE.
      IF (.NOT. LCAL360) THEN

          IF(FIELDCODE(1,FIELD).EQ.1.OR.FIELDCODE(1,FIELD).EQ.2) THEN
            LGREG_MONTHLY=.TRUE.
            UPDATE_MONTHS= FIELDCODE(2,FIELD)*
     &      ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          END IF

      END IF
CCC *ENDIF

          PERIODIC=FIXHD(10,FILE).EQ.2
          REGULAR=.TRUE.
      IF (.NOT. LCAL360) THEN

          REGULAR=FIXHD(35,FILE).EQ.0.AND.FIXHD(36,FILE).
     &    EQ.0
C i.e. data at intervals of days/hours & non-periodic
          IF(PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE).EQ.0
C i.e. data at intervals of hours & periodic

      END IF
          IF(.NOT.PERIODIC) THEN

CL            If data taken from full time series of input data.

            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR
     &                    ,I_MINUTE,I_SECOND
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC
     &                    ,LCAL360)

CCC *IF -DEF,RECON

CL Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*W_SECS_P_P/(W_STEPS_P_P*2)

C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (W_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
      INCR_SEC = -W_SECS_P_P*MOD(OFFSET_STEPS,STEPS(FIELD))/W_STEPS_P_P
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
              IF (W_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
      INCR_SEC = -W_SECS_P_P*MOD(OFFSET_STEPS,STEPS(FIELD))/W_STEPS_P_P
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF

CCC *ENDIF

            IF(REGULAR) THEN
CL 2.1.1  Standard cases:360 day calender;
CL 2.1.1  or Gregorian calendar with
CL        interval between data times in days or hours
CL        updating interval may be regular in model timesteps,
CL        or (LGREG_MONTHLY=T) irregular in model timesteps,

              DAYS   =DAY
              SECONDS=SEC
CL FInd time(in seconds) of first ancillary data on file
              CALL TIME2SEC(FIXHD(21,FILE),FIXHD(22,FILE),
     &                   FIXHD(23,FILE),FIXHD(24,FILE),
     &                   FIXHD(25,FILE),FIXHD(26,FILE),
     &                   ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &                   LCAL360)
              DAYS   =DAYS   -DAY
              SECONDS=SECONDS-SEC
              SECONDS=SECONDS+86400*DAYS

              IF(SECONDS.LT.0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCW: Current time precedes start time of data'
                RETURN
              END IF

CL FInd interval(in seconds) between ancillary data on file
              INTERVAL=(FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE))*3600+
     &                FIXHD(39,FILE)*60+FIXHD(40,FILE)

C Do not interpolate in time if data time exactly matches model time

              IF(MOD(SECONDS,INTERVAL).EQ.0) THEN
                LINTERPOLATE=.FALSE.
              END IF

              STEP=SECONDS/INTERVAL
              TIME=REAL(SECONDS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE

CL 2.1.2 Gregorian calender;ancillary data interval is in months or
CL       years,which is irregular in model timesteps.

CCC *IF -DEF,RECON

CL Adjust YMD time to middle of updating interval

              I_YEAR1=I_YEAR
              I_MONTH1=I_MONTH
              I_DAY1=I_DAY
              I_HOUR1=I_HOUR
              CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,
     &                     I_YEAR,I_MONTH,I_DAY,
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,
     &                     LCAL360)

CCC *ENDIF

CL Find interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
C  Check for time within month
              IF((I_DAY*24+I_HOUR).LT.
     *           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF

              IF(MONTHS.LT.0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCW: Current time precedes start time of data'
                RETURN
              END IF

CCC *IF -DEF,RECON

CL Adjust YMD time back to start of updating interval

              I_YEAR=I_YEAR1
              I_MONTH=I_MONTH1
              I_DAY=I_DAY1
              I_HOUR=I_HOUR1

CCC *ENDIF


              STEP=MONTHS/INTERVAL
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
              SEC=SEC+STEPS(FIELD)*W_SECS_P_P/(W_STEPS_P_P*2)

C  If start-up, adjust for offset of reference time from initial time,
C  & update with values for half a period before first standard update.
              IF (W_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
      INCR_SEC = -W_SECS_P_P*MOD(OFFSET_STEPS,STEPS(FIELD))/W_STEPS_P_P
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
              IF (W_STEP.EQ.0) THEN
                DAY1 = DAY
                SEC1 = SEC
      INCR_SEC = -W_SECS_P_P*MOD(OFFSET_STEPS,STEPS(FIELD))/W_STEPS_P_P
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


CL Adjust YMD time to middle of updating interval

            I_YEAR1=I_YEAR
            I_MONTH1=I_MONTH
            I_DAY1=I_DAY
            I_HOUR1=I_HOUR
            I_MINUTE1=I_MINUTE
            I_SECOND1=I_SECOND
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
CL       data times to be skipped in data set calculated in seconds.

              DAYS   =DAY
              SECONDS=SEC
              INTERVAL=(FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE))*3600+
     &                FIXHD(39,FILE)*60+FIXHD(40,FILE)

C             PERIOD=INTHD(3,FILE)*INTERVAL
              PERIOD=(FIXHD(28,FILE)-FIXHD(21,FILE))*8640*3600
     &             + (FIXHD(29,FILE)-FIXHD(22,FILE))*720*3600
     &             + (FIXHD(30,FILE)-FIXHD(23,FILE))*24*3600
     &             + (FIXHD(31,FILE)-FIXHD(24,FILE))*3600
     &             + (FIXHD(32,FILE)-FIXHD(25,FILE))*60
     &             + (FIXHD(33,FILE)-FIXHD(26,FILE))+INTERVAL
              PERIOD=PERIOD/3600

CL   Do not allow non-standard periods

      IF (LCAL360) THEN

              IF(PERIOD.NE.8640.AND.PERIOD.NE.720.AND.PERIOD.NE.24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCW: Non-standard period for periodic data'
                RETURN
              ENDIF
      ELSE
              IF(PERIOD.NE.24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCW: Non-standard period for periodic data'
                RETURN
              ENDIF
      END IF
              IF(PERIOD.EQ.24)THEN
C Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                IF((I_HOUR*3600+I_MINUTE*60+I_SECOND).LT.(FIXHD(24,FILE)
     &             *3600+FIXHD(25,FILE)*60+FIXHD(26,FILE)))
     &           DAYS=DAYS+1

              ELSE IF(PERIOD.EQ.720)THEN
C Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                IF((I_DAY*24*3600+I_HOUR*3600+I_MINUTE*60+I_SECOND).LT.
     &             (FIXHD(23,FILE)*24*3600+FIXHD(24,FILE)*3600+
     &              FIXHD(25,FILE)*60+FIXHD(26,FILE)))
     &           DAYS=DAYS+30

              ELSE IF(PERIOD.EQ.8640)THEN
C Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                IF((I_MONTH*720*3600+I_DAY*24*3600+I_HOUR*3600+
     &             I_MINUTE*60+I_SECOND).LT.(FIXHD(22,FILE)*720*3600+
     &             FIXHD(23,FILE)*24*3600+FIXHD(24,FILE)*3600+
     &              FIXHD(25,FILE)*60+FIXHD(26,FILE)))
     &           DAYS=DAYS+360

              END IF

              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),
     &                     FIXHD(25,FILE),FIXHD(26,FILE),
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,
     &                     LCAL360)
              DAYS   =DAYS   -DAY
              SECONDS=SECONDS-SEC
              SECONDS=SECONDS+86400*DAYS

C Do not interpolate in time if data time exactly matches model time

              IF(MOD(SECONDS,INTERVAL).EQ.0) THEN
                LINTERPOLATE=.FALSE.
              END IF
              STEP=SECONDS/INTERVAL
              TIME=REAL(SECONDS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE  ! non regular case

CL 2.2.2 Gregorian calender,and data interval is in months,
CL       period is 1 year
CL       Updating interval and number of data times to be skipped
CL       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
C             PERIOD=INTHD(3,FILE)*INTERVAL
              PERIOD=(FIXHD(28,FILE)-FIXHD(21,FILE))*12
     &             + (FIXHD(29,FILE)-FIXHD(22,FILE))+INTERVAL
              IF(PERIOD.NE.12)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCW: Non-standard period for periodic data'
                RETURN
              ENDIF

              MONTHS=I_MONTH-FIXHD(22,FILE)
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
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
C  Calculate  TIME2 for second ancillary data time
C  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IF(IM.LT.I_MONTH) IY=IY+1
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),
     &              FIXHD(25,FILE),FIXHD(26,FILE),
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)

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
            I_MINUTE=I_MINUTE1
            I_SECOND=I_SECOND1


          ENDIF  ! non-periodic/periodic

CCC *IF -DEF,RECON
        IF (LINTERPOLATE) THEN
        WRITE(6,*)' REPLANCW - time interpolation for field ',field
        WRITE(6,*)' time,time1,time2 ',time,time1,time2
        WRITE(6,*)' seconds,int,period ',seconds,interval,period
        END IF
CCC *ENDIF

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
         CMESSAGE='REPLANCO: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
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
           CMESSAGE='REPLANCW: TIME INTERPOLATION ERROR'
           RETURN
          END IF
        END IF

CL  3.  Extract ancillary data for field I  and transfer to D1 array

C set accumulated length of fields input to zero
        LEN_FLD_ACC = 0

C start loop over levels (ends near end of routine)
        DO LEVEL=1,LEVELS(FIELD)

         I2LEV = I2 + LEVEL - 1

CL 3.0 Determine number of values in each field to be input

C  LEN_FLD is length of first field to be read
        LEN_FLD = LOOKUP(LBLREC, LOOKUP_START(FILE) - 1 + I2LEV )

        IF ( LEN_FLD .GT. IMT*JMT ) THEN
          ICODE = 1
          WRITE(6,*) ' length of ancillary field longer than allowed'
          WRITE(6,*) 'LEN_FLD, IMT*JMT = ', LEN_FLD, IMT*JMT
          CMESSAGE='REPLANCW : field length error'
          GO TO 900
        END IF


CL 3.1  Read data for single level of ancillary field.

        CALL READFLDS(NFTIN,1,I2LEV,LOOKUP(1,LOOKUP_START(FILE)),
     &       LEN1_LOOKUP,ANCIL1,LEN_FLD,FIXHD(1,FILE),
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
     &       ICODE,CMESSAGE)
        IF(ICODE.GT.0) THEN
          ICODE=FIELD+100
          CMESSAGE='REPLANCW :I/O error'
          RETURN
        END IF

CLL 3.2 If time interpolation required, read second record
        IF(LINTERPOLATE) THEN

          I1LEV=I2LEV+LOOKUP_STEP(FIELD)

          IF (I1LEV.LE.FIXHD(152,FILE)) THEN
            CALL READFLDS(NFTIN,1,I1LEV,LOOKUP(1,LOOKUP_START(FILE)),
     &        LEN1_LOOKUP,ANCIL2,LEN_FLD,FIXHD(1,FILE),
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
     &        ICODE,CMESSAGE)
              IF(ICODE.GT.0) THEN
                ICODE=FIELD+100
                IOUNIT=NFTIN
                CMESSAGE='REPLANCW :I/O error'
                RETURN
              END IF

          ELSE  ! end of data on file

CL  If end of data has been reached go back to the start if periodic
CL     otherwise cancel time interpolation

            IF(PERIODIC) THEN

C find number of field to access
              I1LEV=NLOOKUP(FIELD) + LEVEL - 1

              CALL READFLDS(NFTIN,1,I1LEV,LOOKUP(1,LOOKUP_START(FILE)),
     &        LEN1_LOOKUP,ANCIL2,LEN_FLD,FIXHD(1,FILE),
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
     &        ICODE,CMESSAGE)
              IF(ICODE.GT.0) THEN
                ICODE=FIELD+300
                IOUNIT=NFTIN
                CMESSAGE='REPLANCW :I/O error'
                RETURN
              END IF
            ELSE
              LINTERPOLATE=.FALSE.
            END IF
          END IF  ! end of position of file test
          ICODE=0
        END IF  ! end of LINTERPOLATE

CL 3.3 Set number of rows of data  (no longer required)

CL 3.4 Perform time interpolation

        IF(LINTERPOLATE) THEN

C Linear interpolation in time, unless missing data indicator
C present at either time.

          CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,
     &             TIME,LEN_FLD)

C If no interpolation, copy data into final array

        ELSE
          DO I=1,LEN_FLD
            ANCIL_DATA(I)=ANCIL1(I)
          END DO
        END IF

CL 3.5 Updating action for each field at each level
CL     Fields replaced.

        POS_STRT = D1_ANCILADD(FIELD) + LEN_FLD_ACC
        DO I=1,LEN_FLD
           D1(POS_STRT+I-1)=ANCIL_DATA(I)
        END DO

CL End loop over levels

        LEN_FLD_ACC = LEN_FLD_ACC + LEN_FLD

       END DO   ! LEVEL  loop ends

CL End loop over ancillary fields (ocean)

      END IF  ! End UPDATE(FIELD) test : 1st level IF block

      END DO

 900  RETURN
      END
