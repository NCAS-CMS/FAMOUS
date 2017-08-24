C******************************COPYRIGHT******************************
C(c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C     Use, duplication or disclosure of this code is subject to the
C     restrictions as set forth in the contract.
C
C     Meteorological Office
C     London Road
C     BRACKNELL
C     Berkshire UK
C     RG12 2SZ
C
CIf no contract has been raised with this copy of the code, the use,
Cduplication or disclosure of it is strictly prohibited.  Permission
Cto do so must first be obtained in writing from the Head of Numerical
CModelling at the above address.
C******************************COPYRIGHT******************************
C
CLL   Routine: INIT_Z_PTR -------------------------------------------
CLL
CLL   Purpose: Initialises address pointers needed by OASIS_STEP when
CLL   coupling the UM with an external model connected by OASIS.
CLL   The fields accessed by the coupler are requiered to be have a
CLL   STASH code and be stored into the D1 array.
CLL
CLL   Tested under compiler:   cft77
CLL   Tested under OS version: UNICOS 9.0.4 (C90)
CLL
CLL  Author:   JC Thil.
CLL
CLL  Code version no: 1.0         Date: 10 Oct 1996
CLL
CLL  Model            Modification history :
CLL  version  date
!LL  4.5     13/01/98 Removed unused AMAXSIZE and IOVARS   P.Burton
CLL
CLL
CLL
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered:
CLL
CLL  Project task:
CLL
CLL  External documentation:
CLL
CLL
CLL  ----------------------------------------------------------------
C*L  Interface and arguments: ---------------------------------------
C
      subroutine ini_z_ptr (
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
! History:                                              
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for MPP.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins                               
     &        D1_ADDR,D1,LD1,ID1, ! IN/OUT:Addressing of D1 & D1 array

C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


C Dump headers
     &A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,
     &A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,
! PP lookup headers and Atmos stash array + index with lengths
     &A_LOOKUP,
     &A_MPP_LOOKUP,
     &a_ixsts, a_spsts,


! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add murk and user ancillary pointers. RTHBarnes
!  4.1  04/12/95  Add pointers JSTHU and JSTHF. J.Smith
!  4.1  26/04/96  Add pointers for Sulphur Cycle variables (12)  MJW
!  4.3   18/3/97  And for HadCM2 sulphate loading patterns.  Will Ingram
!  4.4   05/8/97  And for Conv. cloud amt on model levs. Julie Gregory
!  4.5  04/03/98   Remove pointer SO2_HILEM (add to CARGPT_ATMOS)
!                  Add 1 pointers for NH3 in S Cycle
!                  Add 3 pointers for Soot              M. Woodage
!  4.5  08/05/98   Add 16 new pointers for User Anc.    D. Goddard
!  4.5  13/05/98   Add pointer for RHcrit variable.     S. Cusack
!  4.5  15/07/98   Add pointers for new 3D CO2 array.   C.D.Jones
!  4.5  17/08/98   Remove JSOIL_FLDS and JVEG_FLDS      D. Robinson
C Argument list.
C Pointers for ATMOSPHERE model variables. Configuration dependent.
C
C Addresses in D1 array of primary variables and 'extra' space
C  variable (Exner pressures)
C  Array  variables (dimensions are resolution dependent)
     &       JU, JV, JTHETA, JQ, JQCL, JQCF, J_DEEP_SOIL_TEMP,  JSMCL,
     &       JOZONE, JTRACER, JP_EXNER,
     &       JSO4, JH2SO4, JSOOT, JMURK, JMURK_SOURCE,
     &       JUSER_MULT1, JUSER_MULT2, JUSER_MULT3, JUSER_MULT4,
     &       JUSER_MULT5, JUSER_MULT6, JUSER_MULT7, JUSER_MULT8,
     &       JUSER_MULT9, JUSER_MULT10, JUSER_MULT11, JUSER_MULT12,
     &       JUSER_MULT13, JUSER_MULT14, JUSER_MULT15, JUSER_MULT16,
     &       JUSER_MULT17, JUSER_MULT18, JUSER_MULT19, JUSER_MULT20,
     &       JSTHU, JSTHF,
     &       JSO2,JDMS,JSO4_AITKEN,JSO4_ACCU,JSO4_DISS,JH2O2,
     &       JSO2_NATEM,JOH,JHO2,JH2O2_LIMIT,JO3_CHEM,
     &       JHadCM2_SO4,JCCA,JRHC,JNH3,
     &       JSOOT_NEW,JSOOT_AGD,JSOOT_CLD,JCO2,


     &  internal_model,
     &                      ICODE,CMESSAGE )
C
      IMPLICIT NONE
C
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

C Define Parameters:
      INTEGER  MAX_P_LEVELS     ! Maximum no. of p levels
        PARAMETER (MAX_P_LEVELS = 99  )
      INTEGER  MAX_REQ_THPV_LEVS  ! Max no. of levels for pvort output
        PARAMETER (MAX_REQ_THPV_LEVS = MAX_P_LEVELS )
      INTEGER  MAX_ADJ_TSL      ! Max A_ADJSTEPS
        PARAMETER (MAX_ADJ_TSL  = 10  )
      INTEGER  MAX_N_INTF_A     ! Max no. of atmos interface areas
        PARAMETER (MAX_N_INTF_A =  8  )
      INTEGER  MAX_INTF_LEVELS  ! Max no. of atmos interface levels
        PARAMETER (MAX_INTF_LEVELS = MAX_P_LEVELS )
      INTEGER  MAX_NO_OF_SEGS   ! Maximum number of physics segments
        PARAMETER (MAX_NO_OF_SEGS = 200  )
C     MAX_N_INTF/MAX_N_INTF_A to be sorted out in next version
      INTEGER  MAX_N_INTF     ! Max no. of interface areas
        PARAMETER (MAX_N_INTF =  8  )


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
C*L================ COMDECK TYPSIZE ===========================
C   Description:
C     This COMDECK contains sizes needed for dynamic allocation of
C   main data arrays within the model. Sizes read in from the user
C   interface via NAMELISTs are passed by /COMMON/. Other control
C   sizes that are fundamental in the definition of data structures
C   are assigned by PARAMETER statements.
C
CLL
CLL  Model            Modification history
CLL version  Date
CLL 3.2   30/03/93  New COMDECK created to expedite dynamic allocation
CLL                 of memory. R.Rawlins
CLL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
CLL                  1.Removes the limit on primary STASH item numbers.
CLL                  2.Removes the assumption that (section,item)
CLL                    defines the sub-model.
CLL                  3.Thus allows for user-prognostics.
CLL  3.4    4/07/94  Reduce LEN_DUMPHIST from 2048->0 so that the
CLL                  temporary history file is effectively removed
CLL                  from the dump. R. Rawlins.
!    3.5    MAR. 95  Sub-Models project                                
!                    LEN_STLIST increased to 28 - to allow for         
!                    "internal model identifier".                       
!    4.0    AUG. 95  Introduced S_LEN2_LOOKUP, S_LEN_DATA, 
!                               S_PROG_LOOKUP, S_PROG_LEN
!                       S.J.Swarbrick                                   
!  3.5  17/07/95  Remove ADJ_TIME_SMOOTHING_LENGTH. RTHBarnes.
CLL  4.1  12/03/96  Introduce Wave sub-model.  RTHBarnes.
!  4.1  12/01/96  Add global versions of LBC lengths (the standard
!                 lengths are just local) and U point version
!                 of LENRIMA               MPP code       P.Burton
!  4.2  28/11/96  Increase LEN_STLIST for extra MPP variables
!                 Add global_A_LEN_DATA variable in COMMON
!  4.2  11/10/96  Enable atmos-ocean coupling for MPP.
!                 (1): Coupled fields. Add 'global' sizes for dynamic
!                 allocation of interpolation arrays. R.Rawlins
!  4.2  11/10/96  Enable atmos-ocean coupling for MPP.
!                 (2): Swap D1 memory. Add variables _LEN_D1 to pass 
!                 into coupling routines.           R.Rawlins
!                 Introduce O_LEN_DUALDATA.         S.Ineson
!  4.4  23/09/97  Increment LEN_STLIST due to addition of st_offset_code
!                 S.D. Mullerworth
!  4.4  14/07/97  Add global versions of ocean LBC lengths P.Burton/SI
!  4.4  29/09/97  Add number of levels convective cloud is stored on
!                 N_CCA_LEV                         J.Gregory
!  4.4  29/09/97  New common block to store lengths of Stash Auxillary
!                 arrays and associated index arrays. D. Robinson.
!   4.5  12/09/97 Added LENRIMO_U for ocean boundary velocity 
!                 fields.    C.G. Jones
!  4.5  23/01/98  Increase LEN_STLIST to 33
!  4.5  04/08/98  Add U_FIELD_INTFA. D. Robinson. 
!  4.5  15/04/98  Add new common block MPP_LANDPTS. D. Robinson.
!  4.5  19/01/98  Remove SOIL_VARS and VEG_VARS. D. Robinson.
C All sizes
C Not dependent on sub-model
C     DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
C ATMOS START
C Main sizes of fields for each submodel
C Grid-related sizes for ATMOSPHERE submodel.
      INTEGER
     &       ROW_LENGTH,          ! IN: No of points per row
     &       P_ROWS,              ! IN: No of p-rows
     &       P_LEVELS,            ! IN: No of p-levels
     &       LAND_FIELD           ! IN: No of land points in field
C Physics-related sizes for ATMOSPHERE submodel
      INTEGER
     &       Q_LEVELS,            ! IN: No of moist-levels
     &       CLOUD_LEVELS,        ! IN: No of cloud-levels
     &       ST_LEVELS,           ! IN: No of soil temperature levels
     &       SM_LEVELS,           ! IN: No of soil moisture levels
     &       BL_LEVELS,           ! IN: No of boundary-layer-levels
     &       OZONE_LEVELS         ! IN: No of ozone-levels
     &      ,P_FIELD_CONV         ! IN: field size for conv.incr.copies
C                                 ! 1 if A_CONV_STEP=1, else P_FIELD
C Dynamics-related sizes for ATMOSPHERE submodel
      INTEGER
     &       TR_LEVELS,           ! IN: No of tracer-levels
     &       TR_VARS              ! IN: No of passive tracers
C Dynamics output diagnostic-related sizes for ATMOSPHERE submodel
      INTEGER
     &       THETA_PV_P_LEVS      ! IN: No of levels requested for pvort
C Assimilation-related sizes for ATMOSPHERE submodel  
      INTEGER N_AOBS              ! IN: No. of atmos observation types
C Grid related sizes for data structure
C Data structure sizes for ATMOSPHERE submodel
      INTEGER
     &       A_PROG_LOOKUP,       ! IN: No of prognostic fields
     &       A_PROG_LEN,          ! IN: Total length of prog fields
     &       A_LEN_INTHD,         ! IN: Length of INTEGER header
     &       A_LEN_REALHD,        ! IN: Length of REAL header
     &       A_LEN2_LEVDEPC,      ! IN: No of LEVEL-dependent arrays
     &       A_LEN2_ROWDEPC,      ! IN: No of ROW-dependent arrays
     &       A_LEN2_COLDEPC,      ! IN: No of COLUMN-dependent arrays
     &       A_LEN2_FLDDEPC,      ! IN: No of FIELD arrays
     &       A_LEN_EXTCNST,       ! IN: No of EXTRA scalar constants
     &       A_LEN_CFI1,          ! IN: Length of compressed fld index 1
     &       A_LEN_CFI2,          ! IN: Length of compressed fld index 2
     &       A_LEN_CFI3           ! IN: Length of compressed fld index 3
C ATMOS END
C Data structure sizes for SLAB submodel                                
      INTEGER
     &       S_PROG_LOOKUP,       !IN: No of prognostic fields, SLAB   
     &       S_PROG_LEN           !IN: Tot len of prog fields, SLAB
C SLAB END
C
C OCEAN START
C This *CALL contains TYPOCBAS and COMOCPAR:
C     COMDECK TYPOCPAR
C     ----------------
C   History:         
C   Version   Date     Comment   
C   -------   ----     -------     
C     4.4   15.06.97   Add free surface scalar R.Lenton 
C     COMDECK TYPOCBAS
C     ----------------
C Physics-related sizes for OCEAN submodel
      INTEGER
     &       NT                  ! IN: No of ocean tracers (inc T,S)
C Grid related sizes for OCEAN model
      INTEGER
     &       IMT,                ! IN: No of points per row (incl wrap)
     &       JMT,                ! IN: No of tracer rows
     &       KM                  ! IN: No of tracer levels
C
C      Copies of basic dimensioning variables for OCEAN submodel
C
      INTEGER
     + NT_UI     ! Copy of NT
     +,IMT_UI    ! Copy of IMT
     +,JMT_UI    ! Copy of JMT
     +,KM_UI     ! Copy of KM
CL* COMDECK TYPOASZ; sizes for dynamic allocation of ocean assim.
      INTEGER JO_MAX_OBS_VAL !max number of values in OBS array
     *,JO_LEN_COV            !length of climate/covariances array
     *,JO_MAX_COLS_C     !max number of columns in climate grid
     *,JO_MAX_ROWS_C     !max number of rows    in climate grid
     *,JO_MAX_LEVS_C     !max number of levels  in climate grid
C
      PARAMETER (
     * JO_MAX_OBS_VAL = 1
     *,JO_LEN_COV = 1
     *,JO_MAX_COLS_C = 1
     *,JO_MAX_ROWS_C = 1
     *,JO_MAX_LEVS_C = 1
     *                    )
C
C
C Grid related sizes for OCEAN model
      INTEGER
     &       LSEG,               ! IN: Max no of sets of start/end
C                                      indices for vorticity
     &       NISLE,              ! IN: No of islands
     &       ISEGM,              ! IN: Max no of island segments per box
     &       O_LEN_COMPRESSED,   ! IN: No of ocean points in 3D field
     &       LSEGC               ! IN: No of island basins for mead calc
     &      ,LSEGFS              ! IN: No of start/end indicies for
C                                !     the free surface solution
C Fourier filtering for OCEAN submodel
      INTEGER
     &       LSEGF,    ! IN: max. no of sets of indices for filtering
     &       JFRST,    ! IN: first J row of T to be filtered
     &       JFT0,     ! IN: filtering is done on T with a low
C pass cut off to make the zonal dimension of the box filtered
C effectively the same as that of the boxes on row JFT0
     &       JFT1,     ! IN: last J row of T in SH to be filtered
     &       JFT2,     ! IN: first J row of T in NH to be filtered
     &       JFU0,     ! IN: same function as JFT0 but for U,V
     &       JFU1,     ! IN: last J row of U,V in SH to be filtered
     &       JFU2      ! IN: first J row of U,V in NH to be filtered
C Variables derived from those above
      INTEGER
     &       IMU,      ! IN: total number of U,V grid boxes zonally
     &       IMTP1,    ! IN: IMT+1
     &       IMTM1,    ! IN: IMT-1
     &       IMTM2,    ! IN: IMT-2
     &       IMUM1,    ! IN: IMU-1
     &       IMUM2,    ! IN: IMU-2
     &       JMTP1,    ! IN: JMT+1
     &       JMTM1,    ! IN: JMT-1
     &       JMTM2,    ! IN: JMT-2
     &       JSCAN,    ! IN: JMTM2+1
     &       KMP1,     ! IN: KM+1
     &       KMP2,     ! IN: KM+2
     &       KMM1,     ! IN: KM-1
     &       NSLAB,    ! IN: no of words in one slab
     &       JSKPT,    ! IN: no of rows of T and U,V not filtered in
     &       JSKPU,    ! IN: low and mid latitudes + 1
     &       NJTBFT,   ! IN: no of J rows to be filtered on T
     &       NJTBFU,   ! IN: no of J rows to be filtered on U,V
     &       IMTKM,    ! IN: IMT*KM
     &       NTMIN2    ! IN: maximum of NT or 2
      INTEGER
     &       IMTD2,    ! IN: IMT/2
     &       LQMSUM,   ! IN: IMTD2*(IMT-IMTD2)
     &       LHSUM,    ! IN: IMT*IMTP1/2
     &       IMTX8,    ! IN: IMT*8
     &       IMTIMT    ! IN: IMT*IMT  
      INTEGER
     &       IMROT,    ! X dimension for Coriolis array
     &       JMROT,    ! Y dimension for Coriolis array
     &       IMBC,     ! No of columns in boundary field array
     &       JMBC,     ! No of rows in boundary field array
     &       KMBC,     ! No of levels in boundary field array
     &       NTBC,     ! No of tracers in boundary field array
     &       JMMD,     ! No of rows for mead diagnostic basin indices
     &       LDIV      ! No of divisions mead basin indices
C Grid-related switches for OCEAN submodel
      LOGICAL
     &       CYCLIC_OCEAN,        ! IN: TRUE if CYCLIC E-W boundary
     &       GLOBAL_OCEAN,        ! IN: TRUE if global domain
     &       INVERT_OCEAN         ! IN: TRUE if ocean grid
C                                 !          NS-inverted cf atmos
      PARAMETER
     &      (INVERT_OCEAN=.TRUE.)
C User interface limit for tracers
      INTEGER
     &       O_MAX_TRACERS        ! IN: Max no. tracers in STASHMASTER
      PARAMETER
     &      (O_MAX_TRACERS=20)
C============================= COMDECK COMOCPAR ======================
C
      COMMON /COMOCPAR/ GLOBAL_OCEAN, CYCLIC_OCEAN
     * ,LSEG,NISLE,ISEGM,O_LEN_COMPRESSED,LSEGC,LSEGFS,LSEGF,JFRST,JFT0
     * ,JFT1,JFT2,JFU0,JFU1,JFU2,IMU,IMTP1,IMTM1,IMTM2  
     * ,IMUM1,IMUM2,JMTP1,JMTM1,JMTM2,JSCAN,KMP1,KMP2,KMM1,NSLAB
     * ,JSKPT,JSKPU,NJTBFT,NJTBFU,IMTKM,NTMIN2      
     * ,IMTD2,LQMSUM,LHSUM,IMTX8,IMTIMT                 
     * ,IMROT,JMROT,IMBC,JMBC,KMBC,NTBC,JMMD,LDIV
C
C=====================================================================
C
      INTEGER
     & O_PROG_LOOKUP,O_PROG_LEN,
     & O_LEN_CFI1,O_LEN_CFI2,O_LEN_CFI3,
     & O_LEN_INTHD,O_LEN_REALHD,O_LEN2_LEVDEPC,
     & O_LEN2_ROWDEPC,O_LEN2_COLDEPC,O_LEN2_FLDDEPC,
     & O_LEN_EXTCNST
C OCEAN END

C WAVE SUB-MODEL START                                                  
      INTEGER 
     & NANG,NFRE, ! no.of directions and frequencies of energy spectrum
     & NGX,NGY,   ! length of rows and no.of rows in wave grid
     & NBLO,NIBLO,NOVER,      ! no. & size of blocks, and overlap
     & W_SEA_POINTS,          ! no.of sea points
     & NBLC,NIBLC,NBLD,NIBLD, ! try to remove these later
     & W_PROG_LOOKUP,W_PROG_LEN,                                        
     & W_LEN_CFI1,W_LEN_CFI2,W_LEN_CFI3,                                
     & W_LEN_INTHD,W_LEN_REALHD,W_LEN2_LEVDEPC,                         
     & W_LEN2_ROWDEPC,W_LEN2_COLDEPC,W_LEN2_FLDDEPC,                    
     & W_LEN_EXTCNST 
      LOGICAL   GLOBAL_WAVE   ! true if wave model global  
C WAVE END                                                              

C ATMOS START
C Data structure sizes for ATMOSPHERE ANCILLARY file control routines
      INTEGER
     &       NANCIL_LOOKUPSA      ! IN: Max no of fields to be read
C ATMOS END
C OCEAN START
C Data structure sizes for OCEAN ANCILLARY file control routines
      INTEGER
     &       NANCIL_LOOKUPSO      ! IN: Max no of fields to be read
C OCEAN END

C WAVE START                                                            
C Data structure sizes for WAVE ANCILLARY file control routines         
      INTEGER                                                           
     &       NANCIL_LOOKUPSW      ! IN: Max no of fields to be read     
C WAVE END                                                              
                                                                        
C ATMOS START
C Data structure sizes for ATMOSPHERE INTERFACE file control routines
      INTEGER
     &  N_INTF_A,          ! No of atmosphere interface areas
     &  MAX_INTF_P_LEVELS, ! Max no of model levels in all areas
     &  TOT_LEN_INTFA_P,   ! Total length of interface p grids.
     &  TOT_LEN_INTFA_U    ! Total length of interface u grids.
     & ,U_FIELD_INTFA      ! Length of Model U field (= U_FIELD)
C ATMOS END

                                                                        
      INTEGER
     &  N_INTF_O            ! No of ocean interface areas 

C WAVE START                                                            
C Data structure sizes for WAVE INTERFACE file control routines         
      INTEGER                                                           
     &  N_INTF_W           ! No of atmosphere interface areas           
!     &  ,MAX_INTF_P_LEVELS, ! Max no of model levels in all areas
!     &  TOT_LEN_INTFA_P,   ! Total length of interface p grids.   
!     &  TOT_LEN_INTFA_U    ! Total length of interface u grids.  
C WAVE END                                                              
                                                                        
C ATMOS START
C Data structure sizes for ATMOSPHERE BOUNDARY file control routines
      INTEGER
     &       RIMWIDTHA,           ! IN: No of points width in rim fields
     &       NRIM_TIMESA,         ! IN: Max no of timelevels in rim flds
     &       NFLOOR_TIMESA        ! IN: Max no of t-levs in lwr bdy flds
C ATMOS END
C OCEAN START
C Data structure sizes for OCEAN BOUNDARY file control routines
      INTEGER
     &       RIMWIDTHO,           ! IN: No of points width in rim fields
     &       NRIM_TIMESO          ! IN: Max no of timelevels in rim flds
C OCEAN END
C WAVE START  
C Data structure sizes for WAVE BOUNDARY file control routines   
      INTEGER   
     &       RIMWIDTHW,           ! IN: No of points width in rim fields
     &       NRIM_TIMESW          ! IN: Max no of timelevels in rim flds
C WAVE END 
C Data structure sizes for ATMOS & OCEAN BOUNDARY file control routines
      INTEGER
     &       FLOORFLDSA       ! IN: Total no of lower bndry fields (A)
C
C Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      INTEGER
     &       PP_LEN_INTHD,        ! IN: Length of PP file integer header
     &       PP_LEN_REALHD        ! IN: Length of PP file real    header
C
C
C
C OCEAN sizes common blocks
C============================= COMDECK COMOCBAS ======================
C
      COMMON /COMOCBAS/ NT_UI, IMT_UI, JMT_UI, KM_UI
C
C=====================================================================
C Other sizes passed from namelist into common blocks
      COMMON/NLSIZES/
     & ROW_LENGTH,P_ROWS,LAND_FIELD,P_LEVELS,Q_LEVELS,
     & CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,
     & OZONE_LEVELS,TR_VARS,

     & P_FIELD_CONV,

     & THETA_PV_P_LEVS, N_AOBS, 

     & A_PROG_LOOKUP,A_PROG_LEN,
     & A_LEN_INTHD,A_LEN_REALHD,
     & A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,
     & A_LEN2_FLDDEPC,A_LEN_EXTCNST,
     & A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,

     & S_PROG_LOOKUP, S_PROG_LEN,
   
     & O_PROG_LOOKUP,O_PROG_LEN,
     & O_LEN_CFI1,O_LEN_CFI2,O_LEN_CFI3,
     & O_LEN_INTHD,O_LEN_REALHD,O_LEN2_LEVDEPC,
     & O_LEN2_ROWDEPC,O_LEN2_COLDEPC,O_LEN2_FLDDEPC,
     & O_LEN_EXTCNST,

     & NANG,NFRE,NGX,NGY,NBLO,NIBLO,NOVER,W_SEA_POINTS,
     & NBLC,NIBLC,NBLD,NIBLD, ! try to remove these later
     & W_PROG_LOOKUP,W_PROG_LEN,                                        
     & W_LEN_CFI1,W_LEN_CFI2,W_LEN_CFI3,                                
     & W_LEN_INTHD,W_LEN_REALHD,W_LEN2_LEVDEPC,                         
     & W_LEN2_ROWDEPC,W_LEN2_COLDEPC,W_LEN2_FLDDEPC,                    
     & W_LEN_EXTCNST,                                                   

     & NANCIL_LOOKUPSA,NANCIL_LOOKUPSO,NANCIL_LOOKUPSW, 
                                                                        
     & N_INTF_A,MAX_INTF_P_LEVELS, TOT_LEN_INTFA_P,
     & TOT_LEN_INTFA_U, U_FIELD_INTFA,

     &  N_INTF_O,

     & N_INTF_W,

     & RIMWIDTHA, NRIM_TIMESA,
     & FLOORFLDSA,NFLOOR_TIMESA,
     & RIMWIDTHO, NRIM_TIMESO,         
     & RIMWIDTHW, NRIM_TIMESW,  

     & PP_LEN_INTHD,PP_LEN_REALHD,
     & GLOBAL_WAVE 

C----------------------------------------------------------------------
C     DATA IN STASHC#x MEMBER OF THE JOB LIBRARY

C ATMOS START
C Data structure sizes for ATMOSPHERE submodel (configuration dependent)
      INTEGER
     &       A_LEN2_LOOKUP,       ! IN: Total no of fields (incl diags)
     &       A_LEN_DATA,          ! IN: Total no of words of data
     &       A_LEN_D1             ! IN: Total no of words in atmos D1   
C ATMOS END
C SLAB START                                                            
C Data structure sizes for SLAB       submodel (config dependent)
      INTEGER                                                           
     &       S_LEN2_LOOKUP,       !IN: Tot no of fields (incl diags) 
     &       S_LEN_DATA           !IN: Tot no of words of data       
C SLAB END                                                              
C OCEAN START
C Data structure sizes for OCEAN      submodel (configuration dependent)
      INTEGER
     &       O_LEN2_LOOKUP,       ! IN: Total no of fields (incl diags)
     &       O_LEN_DATA,          ! IN: Total no of words of data
     &       O_LEN_DUALDATA,      ! IN: Words of data at 2 time levels
     &       O_LEN_D1             ! IN: Total no of words in ocean D1
C OCEAN END
C WAVE START                                                            
C Data structure sizes for WAVE       submodel (configuration dependent)
      INTEGER                                                           
     &       W_LEN2_LOOKUP,       ! IN: Total no of fields (incl diags) 
     &       W_LEN_DATA,          ! IN: Total no of words of data
     &       W_LEN_D1             ! IN: Total no of words in atmos D1
C WAVE END                                                              
C Size of main data array for this configuration
      INTEGER
     &       LEN_TOT,             ! IN: Length of D1 array
     &       N_OBJ_D1_MAX         ! IN: No of objects in D1 array
      INTEGER
     &       NSECTS,              ! IN: Max no of diagnostic sections
     &       N_REQ_ITEMS,         ! IN: Max item number in any section
     &       NITEMS,              ! IN: No of distinct items requested
     &       N_PPXRECS,           ! IN: No of PP_XREF records this run
     &       TOTITEMS,            ! IN: Total no of processing requests
     &       NSTTIMS,             ! IN: Max no of STASHtimes in a table
     &       NSTTABL,             ! IN: No of STASHtimes tables
     &       NUM_STASH_LEVELS,    ! IN: Max no of levels in a levelslist
     &       NUM_LEVEL_LISTS,     ! IN: No of levels lists
     &       NUM_STASH_PSEUDO,    ! IN: Max no of pseudo-levs in a list
     &       NUM_PSEUDO_LISTS,    ! IN: No of pseudo-level lists
     &       NSTASH_SERIES_BLOCK, ! IN: No of blocks of timeseries recds
     &       NSTASH_SERIES_RECORDS! IN: Total no of timeseries records

      COMMON/STSIZES/
     &        S_LEN2_LOOKUP,S_LEN_DATA,
     &        A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                        
     &        O_LEN2_LOOKUP,O_LEN_DATA,O_LEN_DUALDATA,O_LEN_D1,         
     &        W_LEN2_LOOKUP,W_LEN_DATA,W_LEN_D1,
     &        LEN_TOT,N_OBJ_D1_MAX,
     &        NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,
     &        NSTTABL,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,
     &        NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS,
     &        NSTTIMS,NSTASH_SERIES_BLOCK,
     &        NSTASH_SERIES_RECORDS

      INTEGER
     &  global_A_LEN_DATA ! global (ie. dump version) of A_LEN_DATA
     &, global_O_LEN_DATA ! global (ie. dump version) of O_LEN_DATA

      COMMON /MPP_STSIZES_extra/
     &       global_A_LEN_DATA,global_O_LEN_DATA


! Sizes of Stash Auxillary Arrays and associated index arrays
!     Initialised in UMINDEX and UMINDEX_A/O/W
      INTEGER LEN_A_IXSTS, LEN_A_SPSTS
      INTEGER LEN_O_IXSTS, LEN_O_SPSTS
      INTEGER LEN_W_IXSTS, LEN_W_SPSTS

      COMMON /DSIZE_STS/ 
     &  LEN_A_IXSTS, LEN_A_SPSTS
     &, LEN_O_IXSTS, LEN_O_SPSTS
     &, LEN_W_IXSTS, LEN_W_SPSTS


!     From 4.5, the number of land points is computed for each
!     PE before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      integer global_land_field    !  Global no of land points
      integer local_land_field     !  Local no of land points
      common /mpp_landpts/ global_land_field,local_land_field

C----------------------------------------------------------------------
C     EXTRA VARIABLES NOT PASSED THROUGH USER INTERFACE
C
C     : FUNDAMENTAL DATA SIZES :
CL   Fundamental parameter  sizes of data structure
C Sizes applicable to all configurations (HISTORY FILE)
      INTEGER
     &       LEN_DUMPHIST         ! IN: Length of history file in dump
      PARAMETER(
     &       LEN_DUMPHIST =    0)
C Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      INTEGER
     &       LEN_FIXHD,           ! IN: Length of dump fixed header
     &       MPP_LEN1_LOOKUP,
     &       LEN1_LOOKUP          ! IN: Size of a single LOOKUP header
      PARAMETER(
     &       LEN_FIXHD    = 256,
     &       MPP_LEN1_LOOKUP= 2,
     &       LEN1_LOOKUP  = 64 )
C Sizes applicable to all configurations (STASH)
      INTEGER
     &       LEN_STLIST,          ! IN: No of items per STASHlist record
     &       TIME_SERIES_REC_LEN  ! IN: No of items per timeseries recd
      PARAMETER(
     &       LEN_STLIST   = 33,
     &       TIME_SERIES_REC_LEN = 9)
      INTEGER
     &        INTF_LEN2_LEVDEPC  ! 1st dim of interface out lev dep cons
      COMMON/DSIZE/
     &        INTF_LEN2_LEVDEPC
C     : SUB-MODEL SIZES        :
C      OUTSIDE *DEF BECAUSE OF THE FUNDAMENTAL ASSUMPTION THAT THE FIRST
C      SECTION.
      INTEGER
     &       MOS_MASK_LEN         ! IN: Size of bit mask for MOS
      COMMON/DSIZE_AO/
     &        MOS_MASK_LEN
C     : SUB-MODEL ATMOSPHERE   :
      INTEGER
C Data structure sizes derived from grid size
     &       A_LEN1_LEVDEPC,      ! IN: 1st dim of level  dep const
     &       A_LEN1_ROWDEPC,      ! IN: 1st dim of row    dep const
     &       A_LEN1_COLDEPC,      ! IN: 1st dim of column dep const
     &       A_LEN1_FLDDEPC       ! IN: 1st dim of field  dep const
C Data structure sizes for ATMOSPHERE INTERFACE file control routines
      INTEGER
     &       INTF_LOOKUPSA        ! No of interface lookups.
      COMMON/DSIZE_A/
     &        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,
     &        A_LEN1_COLDEPC,
     &        INTF_LOOKUPSA
C     : SUB-MODEL ATMOSPHERE   : DERIVED SIZES
C Parameters derived from model grid/levels. Arakawa B-grid
      INTEGER
     &       P_FIELD,             ! IN: No of p-points in field
     &       U_ROWS,              ! IN: No of uv-rows
     &       U_FIELD              ! IN: No of uv-points in field
     &,      N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/
     &        P_FIELD,U_FIELD,U_ROWS,N_CCA_LEV

C     : BOUNDARY UPDATING      : DERIVED VALUES
      INTEGER
     & LENRIMA,LENRIMO,  ! No.of pts.in horiz.strip round bdy. (A&O)
     & LENRIMO_U,        ! No. of points in for velocity fields  
     & LENRIMA_U,        ! Similarly for atmosphere U points
     & RIMFLDSA,RIMFLDSO,! Total no.of fields in lateral bdy.d/s. (A&O)
     & BOUNDFLDS,        ! Total no.of indep.updated groups of bdy.flds.
     & RIM_LOOKUPSA,     ! Total no.of PP headers describing bdy.data(A)
     & RIM_LOOKUPSO,     ! Total no.of PP headers describing bdy.data(O)
     & BOUND_LOOKUPSA,   ! Total no.of PP headers describing fields (A)
     & BOUND_LOOKUPSO,   ! Total no.of PP headers describing fields (O)
     & BOUND_LOOKUPSW,   ! Total no.of PP headers describing fields (W) 
     & LENRIMDATA_A,     ! Length of lat.bdy.data for a single time (A)
     & LENRIMDATA_O      ! Length of lat.bdy.data for a single time (O)
      COMMON/DRSIZ_BO/
     & LENRIMA,LENRIMO,LENRIMA_U,LENRIMO_U,RIMFLDSA,RIMFLDSO,BOUNDFLDS,
     & RIM_LOOKUPSA,RIM_LOOKUPSO,BOUND_LOOKUPSA,BOUND_LOOKUPSO,
     & BOUND_LOOKUPSW,LENRIMDATA_A,LENRIMDATA_O 
! The above variables all refer to local data sizes. For the MPP code
! we also require a few global lengths to dimension arrays with
! before the data is distributed to local processors
      INTEGER
     &  global_LENRIMA      ! global version of LENRIMA
     &, global_LENRIMDATA_A ! global version of LENRIMDATA_A
     &, global_LENRIMO ! global version of LENRIMO
     &, global_LENRIMDATA_O ! global version of LENRIMDATA_O
     
      COMMON /MPP_global_DRSIZE_BO/
     & global_LENRIMA,global_LENRIMDATA_A
     &, global_LENRIMO,global_LENRIMDATA_O

! History:                                              
! Version  Date    Comment
!  4.2  11/10/96   Enable atmos-ocean coupling for MPP.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins 
!  4.5  18/09/98   Modified name of COMMON block to stop clash with
!                  identical Fortran vairable name         P.Burton
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
C     Common block containing the ALT_N_SUBMODEL_PARTITION variables
C     COMDECK CALTSUBM:
C     COMDECK TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX
C     in CSUBMODL. However, they are not always called in the same
C     decks and in the right order. Therefore, copy the values to
C     another comdeck and *CALL it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION
      INTEGER ALT_N_SUBMODEL_PARTITION_MAX

      PARAMETER(ALT_N_SUBMODEL_PARTITION_MAX=4)

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
CL                           to be called in the same module.
      REAL     D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

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


C     D1 addressing array and number of objects in each submodel
      INTEGER
     &  D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,ALT_N_SUBMODEL_PARTITION)

      INTEGER
     &  NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
CL
CL COMDECKS TYPSIZE and CSUBMODL must be *CALLed before this comdeck
CL
CL
C Applicable to all configurations (except MOS variables)
C STASH related variables for describing output requests and space
C management.
CLL
CLL   AUTHOR            Rick Rawlins
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.2             Code creation for Dynamic allocation
CLL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
CLL                  1.Removes the limit on primary STASH item numbers.
CLL                  2.Removes the assumption that (section,item)
CLL                    defines the sub-model.
CLL                  3.Thus allows for user-prognostics.
CLL   3.5  Apr. 95   Sub-Models project.
CLL                  Dimensioning of various STASH arrays altered in
CLL                  accordance with internal model separation scheme.
CLL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
CLL                  longer required.
CLL                  S.J.Swarbrick
CLL
C
CC This *CALL is needed to get ppxref_codelen to dimension PP_XREF
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
C
C Scalars defining sizes in STASH used for defining local array
C dimensions at a lower level.
      INTEGER
     &       MAX_STASH_LEVS,  ! Maximum no of output levels for any diag
     &       PP_LEN2_LOOKUP,  ! Maximum no of LOOKUPs needed in STWORK
     &       MOS_OUTPUT_LENGTH
      COMMON/CARGST/
     &       MAX_STASH_LEVS,  ! Maximum no of output levels for any diag
     &       PP_LEN2_LOOKUP,  ! Maximum no of LOOKUPs needed in STWORK
     &       MOS_OUTPUT_LENGTH
C
      LOGICAL
     &       SF(0:NITEMS,0:NSECTS) ! STASHflag (.TRUE. for processing
C                                  ! this timestep). SF(0,IS) .FALSE.
C                                  ! if no flags on for section IS.
      INTEGER
! STASH list index
     &       STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL),

! List of STASH output requests
     &       STLIST (LEN_STLIST,TOTITEMS),

! Address of item from generating plug compatible routine (often workspa
     &       SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL),

! STASH times tables
     &       STTABL (NSTTIMS,NSTTABL),

! Length of STASH workspace required in each section
     &       STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          ),
     &       PPINDEX            (  NITEMS,N_INTERNAL_MODEL          ),
     &       STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS ),
     &       STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS),
     &       STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS),
     &       STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK),
     &       MOS_MASK(MOS_MASK_LEN)
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL                           to be called in the same module.
!LL  Model            Modification history
!LL version  Date
!LL   4.1    21/03/96 Added arrays to hold local lengths and addresses
!LL                   for MPP code
!LL
CL --------------- Dump headers (atmosphere)-------------
      INTEGER
C                                                 ! IN/OUT:
     &A_FIXHD(LEN_FIXHD),                         ! fixed length header
     &A_INTHD(A_LEN_INTHD),                       ! integer header
     &A_CFI1(A_LEN_CFI1+1),                       ! compress field index
     &A_CFI2(A_LEN_CFI2+1),                       ! compress field index
     &A_CFI3(A_LEN_CFI3+1)                        ! compress field index

      REAL
C                                                 ! IN/OUT:
     &A_REALHD(A_LEN_REALHD),                     ! real header
     &A_LEVDEPC(A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1), ! level  dep const
     &A_ROWDEPC(A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1), ! row    dep const
     &A_COLDEPC(A_LEN1_COLDEPC*A_LEN2_COLDEPC+1), ! column dep const
     &A_FLDDEPC(A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1), ! field  dep const
     &A_EXTCNST(A_LEN_EXTCNST+1),                 ! extra constants
     &A_DUMPHIST(LEN_DUMPHIST+1)                  ! temporary hist file

CL --------------- PP headers ---------------------------
      INTEGER
     &A_LOOKUP(LEN1_LOOKUP,A_LEN2_LOOKUP)         ! IN/OUT: lookup heads
     &, A_MPP_LOOKUP(MPP_LEN1_LOOKUP,A_LEN2_LOOKUP)
     &, a_ixsts(len_a_ixsts)                      ! stash index array

      REAL
     &  a_spsts(len_a_spsts)                      ! atmos stash array

CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL                           to be called in the same module.
!LL  Model            Modification history
!LL version  Date
!LL   4.1    21/03/96 Added arrays to hold local lengths and addresses
!LL                   for MPP code
!LL
CL --------------- Dump headers (ocean) -----------------
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL                           to be called in the same module.
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add pointers to new slab prognostics. J F Thomson.
!  3.5   19/05/95  Remove pointers JK1,JK2,JEXPK1,JEXPK2,JKDA,JKDF
!                  and JRHCRIT. D. Robinson
!  4.1   13/10/95  Add pointers for new soil moisture fraction
!                  prognostics,canopy conductance and
!                  vegetation J.Smith
!  4.1   26/04/96  Add pointers for Sulphur Cycle (12)   MJWoodage
!  4.3    18/3/97  Add pointers for HadCM2 sulphate loading patterns
!                                                   William Ingram
!  4.4   05/09/97  Add pointer for net energy flux prognostic 
!                  S.D.Mullerworth
!  4.4   05/08/97  Add pointer for conv. cld amt on model levs. JMG   
!  4.4  10/09/97   Added pointers for snow grain size and snow soot
!                  content used in prognostic snow albedo scheme
!                                                        R. Essery
!  4.4    16/9/97  Add pointers for new vegetation and land surface
!                  prognostics.                       Richard Betts
!  4.5    1/07/98  Add pointer for ocean CO2 flux. C.D.Jones
!  4.5    19/01/98 Replace JVEG_FLDS and JSOIL_FLDS with
!                  individual pointers. D. Robinson
!  4.5    04/03/98 Add 2 pointers for NH3 in S Cycle     M. Woodage
!                  Add 5 pointers for soot               M. Woodage
!                  Add pointer SO2_HILEM to CARGPT_ATMOS M. Woodage
!  4.5    08/05/98 Add 16 pointers for User Anc.         D. Goddard
!  4.5    13/05/98 Add pointer for RHcrit variable.      S. Cusack
!  4.5    15/07/98 Add pointers for new 3D CO2 array.    C.D.Jones
C Type definition.
C Pointers for ATMOSPHERE model variables. Configuration dependent.
C Addresses in D1 array of primary variables and 'extra' space
C  variable (Exner pressures)
C
C        Array  variables (dimensions are resolution dependent)
      INTEGER
     &       JU(P_LEVELS),
     &       JV(P_LEVELS),
     &       JTHETA(P_LEVELS),
     &       JQ(Q_LEVELS),
     &       JQCL(Q_LEVELS),
     &       JQCF(Q_LEVELS),
     &       JCCA(N_CCA_LEV), ! conv cld amt on model levs.
     &       JRHC(Q_LEVELS),

     &       J_DEEP_SOIL_TEMP(ST_LEVELS),
     &       JSMCL(SM_LEVELS),       !soil moisture content in layers
     &       JSTHU(SM_LEVELS),       ! unfrozen soil moisture fraction
     &       JSTHF(SM_LEVELS),       ! frozen soil moisture fraction

     &       JOZONE(OZONE_LEVELS),
     &       JTRACER(TR_LEVELS,TR_VARS+1),

     &       JMURK_SOURCE(P_LEVELS),  ! multi-level murk source
     &       JMURK(P_LEVELS),     ! multi-level murk concentration
     &       JSO4(TR_LEVELS),     ! (ammonium) sulphate aerosol
     &       JH2SO4(TR_LEVELS),   ! sulphuric acid aerosol
     &       JSOOT(TR_LEVELS),    ! soot aerosol
     &       JSO2(P_LEVELS),         ! sulphur dioxide gas
     &       JDMS(P_LEVELS),         ! dimethyl sulphide gas
     &       JSO4_AITKEN(P_LEVELS),  ! Aitken mode sulphate aerosol
     &       JSO4_ACCU(P_LEVELS),    ! accumulation mode sulphate aer
     &       JSO4_DISS(P_LEVELS),    ! dissloved  sulphate aerosol
     &       JH2O2(P_LEVELS),        ! hydrogen peroxide mmr
     &       JNH3(P_LEVELS),         ! ammonia gas mmr
     &       JSOOT_NEW(P_LEVELS),    ! fresh soot mmr
     &       JSOOT_AGD(P_LEVELS),    ! aged soot mmr
     &       JSOOT_CLD(P_LEVELS),    ! soot in cloud mmr
     &       JSO2_NATEM(P_LEVELS),   ! natural SO2 emissions
     &       JOH(P_LEVELS),          ! hydroxyl radical ancillary
     &       JHO2(P_LEVELS),         ! hydrogen dioxide ancillary
     &       JH2O2_LIMIT(P_LEVELS),  ! limiting H2O2 ancillary
     &       JO3_CHEM(P_LEVELS),     ! ozone for chemistry ancillary
     &       JHadCM2_SO4(2),         ! HadCM2 sulphate loading patterns
     &       JCO2(P_LEVELS),         ! 3D CO2 FIELD
     &       JUSER_MULT1(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT2(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT3(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT4(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT5(P_LEVELS),  ! multi-level user ancillary
     &       JUSER_MULT6(P_LEVELS),  ! multi-level user ancillary
     &       JUSER_MULT7(P_LEVELS),  ! multi-level user ancillary
     &       JUSER_MULT8(P_LEVELS),  ! multi-level user ancillary
     &       JUSER_MULT9(P_LEVELS),  ! multi-level user ancillary
     &       JUSER_MULT10(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT11(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT12(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT13(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT14(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT15(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT16(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT17(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT18(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT19(P_LEVELS), ! multi-level user ancillary
     &       JUSER_MULT20(P_LEVELS), ! multi-level user ancillary
     &       JP_EXNER(P_LEVELS+1)      ! Exner pressure

C        Scalar variables
      INTEGER
     &       JPSTAR,                 ! surface pressure
     &       JSMC,                   ! soil moisture content
     &       JCANOPY_WATER,
     &       JSNSOOT,                ! snow soot content
     &       JRGRAIN,                ! snow grain size
     &       JSNODEP,                ! snow depth
     &       JTSTAR,                 ! surface temperature
     &       JTI,                    ! Sea-ice temperature (RE 27/7/94)
     &       JTSTAR_ANOM,
     &       JZH,                    ! boundary layer depth
     &       JZ0,                    ! roughness length
     &       JLAND,                  ! land sea mask
     &       JICE_FRACTION,
     &       JICE_THICKNESS,
     &       JTCLIM,
     &       JHCLIM,
     &       JICE_EDGE,
     &       JSAT_SOILW_SUCTION,     ! saturated soil water suction
     &       JLAI,                   ! Gridbox mean leaf area index
     &       JCANHT,                 ! Gridbox mean canopy height
     &       JFRAC_TYP,              ! Fractions of surface types
     &       JLAI_PFT,               ! LAI of plant functional types
     &       JCANHT_PFT,             ! Canopy hght of plant func types
     &       JGS,                    ! Gridbox mean canopy conductance
     &       JDISTURB,               ! Disturbed fraction of vegetation
     &       JVOL_SMC_WILT,          ! vol smc at wilting
     &       JVOL_SMC_CRIT,          ! vol smc at critical point
     &       JVOL_SMC_FCAP,          ! vol smc at field capacity
     &       JVOL_SMC_SAT,           ! vol smc at saturation
     &       JSAT_SOIL_COND,         ! saturated soil conductivity
     &       JEAGLE_EXP,             ! eagle's exponent
     &       JTHERM_CAP,             ! thermal capacity
     &       JTHERM_COND,            ! thermal conductivity
     &       JCLAPP_HORN,            ! clapp hornberger B coeff
     &       JVEG_FRAC,              ! vegetation fraction
     &       JROOT_DEPTH,            ! root depth
     &       JSFA,                   ! snow free albedo
     &       JMDSA,                  ! cold deep snow albedo
     &       JSURF_RESIST,           ! surafce resistance
     &       JSURF_CAP,              ! surface capacity
     &       JINFILT,                ! infiltration factor
     &       JSOIL_ALB,              ! Snow-free albedo of bare soil
     &       JSOIL_CARB,             ! Soil carbon content
     &       JNPP_PFT_ACC,           ! Accumulated NPP on PFTs
     &       JG_LF_PFT_ACC,          ! Accum. leaf turnover rate PFTs
     &       JG_PHLF_PFT_ACC,        ! Accumulated phenological leaf
C                                    ! turnover rate on PFTs
     &       JRSP_W_PFT_ACC,         ! Accum. wood respiration on PFTs
     &       JRSP_S_ACC,             ! Accumulated soil respiration 
     &       JTSNOW,                 ! Snow surface layer temperature
     &       JCAN_WATER_NIT,         ! Canopy water content on non-ice 
C                                    ! tiles
     &       JCATCH_NIT,             ! Canopy capacity on non-ice tiles
     &       JTSTAR_TYP,             ! Surface temperature on tiles
     &       JZ0_TYP,                ! Surface roughness on tiles
     &       JOROG,          ! orographic height
     &       JOROG_SD,       ! standard deviation of orography
     &       JOROG_Z0,       ! orographic roughness length (old version)
     &       JOROG_SIL,      ! silhouette area of orography
     &       JOROG_HO2,      ! peak to trough height/(2*sqrt2)

     &       JU_SEA,          ! Surface current (u component)
     &       JV_SEA,          ! Surface current (v component)

     &       JTSLAB,          ! Temperature of slab ocean.
     &       JUICE,           ! X component of ice velocity.
     &       JVICE,           ! Y component of ice velocity.

     &       JCCB,            ! convective cloud base
     &       JCCT,            ! convective cloud top
     &       JCCLWP,          ! convective cloud liquid water path
     &       JSO2_EM,         ! sulphur dioxide emission
     &       JDMS_EM,         ! dimethyl sulphur emission
     &       JSO2_HILEM,      ! high level SO2 emissions
     &       JNH3_EM,         ! ammonia gas surface emiss 
     &       JSOOT_EM,        ! fresh soot surface emissions
     &       JSOOT_HILEM,     ! fresh soot high lev emissions
     &       JNET_FLUX        ! Net energy flux

      INTEGER
     &       JUSER_ANC1,      ! user ancillary field 1
     &       JUSER_ANC2,      ! user ancillary field 2
     &       JUSER_ANC3,      ! user ancillary field 3
     &       JUSER_ANC4,      ! user ancillary field 4
     &       JUSER_ANC5,      ! user ancillary field 5
     &       JUSER_ANC6,      ! user ancillary field 6
     &       JUSER_ANC7,      ! user ancillary field 7
     &       JUSER_ANC8,      ! user ancillary field 8
     &       JUSER_ANC9,      ! user ancillary field 9
     &       JUSER_ANC10,     ! user ancillary field 10
     &       JUSER_ANC11,     ! user ancillary field 11
     &       JUSER_ANC12,     ! user ancillary field 12
     &       JUSER_ANC13,     ! user ancillary field 13
     &       JUSER_ANC14,     ! user ancillary field 14
     &       JUSER_ANC15,     ! user ancillary field 15
     &       JUSER_ANC16,     ! user ancillary field 16
     &       JUSER_ANC17,     ! user ancillary field 17
     &       JUSER_ANC18,     ! user ancillary field 18
     &       JUSER_ANC19,     ! user ancillary field 19
     &       JUSER_ANC20,     ! user ancillary field 20

     &       JRIM,              ! Lateral boundary update fields
     &       JRIM_TENDENCY,     ! Lateral boundary tendencies

     &       JOROG_TENDENCY,    ! Orographic tendencies
     &       JOROG_SD_TENDENCY, ! Orographic variable tendency
     &       JOROG_GRAD_XX,
     &       JOROG_GRAD_XY,
     &       JOROG_GRAD_YY
     &,    J_CO2FLUX     ! Ocean CO2 flux (Kg CO2/m2/s1)
     &,    J_CO2_EMITS      ! Surface CO2 emissions (Kg CO2/m2/s1)

C Addresses in D1 array of primary variables: scalars
      COMMON/CARGPT_ATMOS/
     &  JPSTAR, JSMC,  JCANOPY_WATER, JSNODEP, JTSTAR, JTI,
     &  JSNSOOT, JRGRAIN,
     &  JTSTAR_ANOM, JZH, JZ0, JLAND, JICE_FRACTION,
     &  JGS, JCANHT, JLAI,
     &  JICE_THICKNESS, JTCLIM, JHCLIM, JICE_EDGE, JSAT_SOILW_SUCTION,
     &  JVOL_SMC_WILT, JVOL_SMC_CRIT, JVOL_SMC_FCAP, JVOL_SMC_SAT,
     &  JSAT_SOIL_COND, JEAGLE_EXP, JTHERM_CAP, JTHERM_COND, 
     &  JVEG_FRAC, JROOT_DEPTH, JSFA, JMDSA, JSURF_RESIST,
     &  JSURF_CAP, JINFILT, JCLAPP_HORN,
     &  JOROG, JOROG_SD, JOROG_Z0, JU_SEA, JV_SEA,
     &  JTSLAB,JUICE,JVICE,
     &  JCCB, JCCT, JCCLWP, JSO2_EM, JDMS_EM,
     &  JSO2_HILEM, JNH3_EM, JSOOT_EM, JSOOT_HILEM,
     &  JOROG_SIL, JOROG_HO2, JNET_FLUX,
     &  JUSER_ANC1, JUSER_ANC2, JUSER_ANC3, JUSER_ANC4, JUSER_ANC5,
     &  JUSER_ANC6, JUSER_ANC7, JUSER_ANC8, JUSER_ANC9, JUSER_ANC10,
     &  JUSER_ANC11, JUSER_ANC12, JUSER_ANC13, JUSER_ANC14, JUSER_ANC15,
     &  JUSER_ANC16, JUSER_ANC17, JUSER_ANC18, JUSER_ANC19, JUSER_ANC20,
     &  JRIM, JRIM_TENDENCY, JOROG_TENDENCY, JOROG_SD_TENDENCY,
     &  JOROG_GRAD_XX, JOROG_GRAD_XY, JOROG_GRAD_YY, JFRAC_TYP,
     &  JLAI_PFT, JCANHT_PFT, JDISTURB, JSOIL_ALB, JSOIL_CARB, 
     &  JNPP_PFT_ACC, JG_LF_PFT_ACC, JG_PHLF_PFT_ACC, JRSP_W_PFT_ACC, 
     &  JRSP_S_ACC, JTSNOW, JCAN_WATER_NIT, JCATCH_NIT, JTSTAR_TYP, 
     &  JZ0_TYP
     &,    J_CO2FLUX
     &,    J_CO2_EMITS

C Pointers for ATMOSPHERE model constants. Scalars only.
C Addresses in level dependent constants array.
      INTEGER
     &       JAK,             ! Mid layer values defining
     &       JBK,             ! Hybrid coordinates
     &       JDELTA_AK,       ! Layer
     &       JDELTA_BK,       ! thickness
     &       JTHETA_REF,      ! Reference temperature profile for
C                             ! split-explicit time integrations
     &       JSOIL_THICKNESS, ! Thickness of deep soil layers
     &       JFILTER_WAVE_NUMBER_P_ROWS, ! holds last wave number ^ on p
     &       JFILTER_WAVE_NUMBER_U_ROWS, ! not to be chopped      ^ on u
     &       JNSWEEP          ! No. of E-W sweeps/row (tracer advection)
C Pointers for ATMOSPHERE model constants. Scalars only.
      COMMON/CARGPT_ATMOS/
C Addresses in level dependent constants array.
     &       JAK, JBK, JDELTA_AK, JDELTA_BK, JTHETA_REF,
     &       JSOIL_THICKNESS,
C Addresses in row   dependent constants array.
     &       JFILTER_WAVE_NUMBER_P_ROWS, JFILTER_WAVE_NUMBER_U_ROWS,
     &       JNSWEEP
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL                           to be called in the same module.
CLL History:
CLL Version  Date  Comment
CLL  3.4   18/5/94 Remove sea ice flux correction. J F Thomson.
CLL  4.4   15/06/97 introduce pointers for the free surface solution.
CLL                                                         R.Lenton
CLL  4.5  04/08/97 Add pointers for ocean boundary data in D1 array
CLL                                                C.G. Jones
CLL  4.5    1/07/98 Add new pointer to atmospheric CO2 field. C.D.Jones
C
      integer
     &  internal_model

      INTEGER ICODE             ! OUT - Error return code
      CHARACTER*(*) CMESSAGE    ! OUT - Error return message


C
C -------------------------------------------------------------------
C

C
C*-------------------------------------------------------------------

C
C     Common blocks
C
CLL  Comdeck: CAOPTR -------------------------------------------------
CLL
CLL  Purpose: Holds address pointers for atmosphere-to-ocean coupling
CLL           fields required by SWAP_A2O and SWAP_O2A.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   T.C.Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.5    1/07/98  Add pointers for interactive CO2 coupling
CLL                       fields.                  C.D.Jones
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: C100
CLL
CLL  Project task: C0
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C0 - The top-level control system
CLL
CLL  -------------------------------------------------------------------
C
CL 1. Pointers needed by SWAP_A2O and SWAP_O2A -------------------------

      INTEGER
     *    JO_SST,            ! Sea-surface temperature on ocean grid
     *    JO_UCURR           ! Surface zonal current on ocean grid
C
      COMMON /AOPTR/
     &    JO_SST,JO_UCURR

CL 2. Pointers needed by SWAP_A2O -------------------------------------

      INTEGER
     &    JA_TAUX,           ! Surface x-windstress on atmos grid
     &    JO_TAUX,           ! Surface x-windstress on ocean grid
     &    JA_TAUY,           ! Surface y-windstress on atmos grid
     &    JO_TAUY,           ! Surface y-windstress on ocean grid
     &    JA_WINDMIX,        ! Windmixing power on atmos grid
     &    JO_WINDMIX,        ! Windmixing power on ocean grid
     &    JA_SOLAR,          ! Net downward SW at surf on atmos grid
     &    JA_BLUE,           ! Net blueband SW at surf on atmos grid
     &    JO_BLUE,           ! Net blueband SW at surf on ocean grid
     &    JA_EVAP,           ! Net evaporation over sea on atmos grid
     &    JA_LONGWAVE,       ! Net downward LW at surf on atmos grid
     &    JA_SENSIBLE,       ! Sensible heat flux over sea on atmos grid
     &    JO_HEATFLUX,       ! Non penetrative heatflux on ocean grid
     &    JA_LSSNOW,         ! Large-scale snowfall rate on atmos grid
     &    JA_CVSNOW,         ! Convective snowfall rate on atmos grid
     &    JA_LSRAIN,         ! Large-scale rainfall rate on atmos grid
     &    JA_CVRAIN,         ! Convective rainfall rate on atmos grid
     &    JO_PMINUSE,        ! Precipitation-evaporation on ocean grid
     &    JA_SLOWRUNOFF,     ! Slow (sub-surface) runoff on atmos grid
     &    JA_FASTRUNOFF,     ! Fast (surface) runoff on atmos grid
     &    JA_OCENTPTS,       ! Ocean entry point index to atmos landpts
     &    JO_RIVEROUT        ! Total river outflow on ocean grid
     &,   JA_co2        ! atmos level 1 co2 conc.
     &,   JO_co2
     &,   JA_co2flux    ! ocean co2 flux.
     &,   JO_co2flux
C
      COMMON /A2OPTR/
     &    JA_TAUX,JO_TAUX,JA_TAUY,JO_TAUY,JA_WINDMIX,JO_WINDMIX,
     &    JA_SOLAR,JA_BLUE,JO_BLUE,JA_EVAP,JA_LONGWAVE,JA_SENSIBLE,
     &    JO_HEATFLUX,JA_LSSNOW,JA_CVSNOW,JA_LSRAIN,JA_CVRAIN,JO_PMINUSE
     &,   JA_SLOWRUNOFF,JA_FASTRUNOFF,JA_OCENTPTS,JO_RIVEROUT
     &,   JA_co2, JO_co2, JA_co2flux, JO_co2flux

CL 3. Pointers needed by SWAP_O2A -------------------------------------

      INTEGER
     &    JO_TSTAR,          ! Surface temperature on ocean grid
     &    JA_TSTAR,          ! Surface temperature on atmos grid
     &    JA_UCURR,          ! Surface zonal current on atmos grid
     *    JO_VCURR,          ! Surface merid current on ocean grid
     &    JA_VCURR           ! Surface merid current on atmos grid
C
      COMMON /O2APTR/
     &    JO_TSTAR,JA_TSTAR,JA_UCURR,JO_VCURR,JA_VCURR
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
CLL  Comdeck: STPARAM --------------------------------------------------
CLL
CLL  Purpose: Meaningful PARAMETER names for STASH processing routines.
CLL           Both a long name and short name have been declared, to
CLL           reduce the existence of "magic" numbers in STASH.
CLL           Format is that first the address of the item is declare in
CLL           both long and short form. example is;
CLL             integer st_item_code,s_item  !Item number (declaration)
CLL             parameter(st_item_code=3,s_item=3)
CLL
CLL  Author:   S.Tett             Date:           22 January 1991
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.5    Mar. 95  Sub-models project.
CLL                   st_model_code=28 added to STLIST addresses
CLL                                   S.J.Swarbrick
!LL   4.2    27/11/96 MPP code: Added new stlist "magic numbers" :
!LL                   st_dump_output_length, st_dump_output_addr
!LL                                                       P.Burton
!LL   4.4    23/09/97 Add st_offset_code to the STASH list
!LL                   S.D. Mullerworth
!    4.4  02/12/96 Time mean timeseries added R A Stratton.             
!    4.5  23/01/98 Added new stlist magic number
!                  st_dump_level_output_length
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
CLLEND--------------------------------------------------------------
C
         integer st_model_code,s_modl ! Internal model number address
         parameter(st_model_code=28,s_modl=28)

         integer st_sect_no_code,s_sect ! Section Number address
     &          ,st_sect_code
         parameter(st_sect_no_code=2,s_sect=2,st_sect_code=2)

         integer st_item_code,s_item  !Item number address
         parameter(st_item_code=1,s_item=1)

         integer st_proc_no_code,s_proc ! Processing Code address
         parameter(st_proc_no_code=3,s_proc=3)

CL subsidiary codes for st_proc_no_code now

         integer st_replace_code
         parameter(st_replace_code=1)

         integer st_accum_code
         parameter(st_accum_code=2)

         integer st_time_mean_code
         parameter(st_time_mean_code=3)

         integer st_time_series_code
         parameter(st_time_series_code=4)

         integer st_max_code
         parameter(st_max_code=5)

         integer st_min_code
         parameter(st_min_code=6)

         integer st_append_traj_code
         parameter(st_append_traj_code=7)

         integer st_time_series_mean                                    
         parameter(st_time_series_mean=8)                               
                                                                        
         integer st_variance_code
         parameter(st_variance_code=9)                                  

         integer st_freq_code,s_freq ! Frequency (Input & output) addres
         parameter(st_freq_code=4,s_freq=4)

         integer st_offset_code,s_offs ! Offset for sampling
         parameter(st_offset_code=30,s_offs=30)

         integer st_start_time_code,s_times ! start timestep address
         parameter(st_start_time_code=5,s_times=5)

         integer st_end_time_code,s_timee ! end timestep address
         parameter(st_end_time_code=6,s_timee=6)

         integer st_period_code,s_period ! period in timesteps address
         parameter(st_period_code=7,s_period=7)

         integer st_infinite_time        ! infinite end/period value
         parameter(st_infinite_time=-1)

         integer st_end_of_list          ! end-of-list marker in times
         parameter(st_end_of_list=-1)

C ---------------------------- grid point stuff
         integer st_gridpoint_code,s_grid ! gridpoint info address
         parameter(st_gridpoint_code=8,s_grid=8)
CL now subsid grid point stuff
         integer stash_null_mask_code,s_nomask ! no masking done
         parameter(stash_null_mask_code=1,s_nomask=1)

         integer stash_land_mask_code,s_lndms ! land mask conds
         parameter(stash_land_mask_code=2,s_lndms=2)

         integer stash_sea_mask_code,s_seams  ! sea mask code
         parameter(stash_sea_mask_code=3,s_seams =3)

CL processing options

         integer block_size ! size of block for gridpoint code
         parameter(block_size=10)

         integer extract_top ! max code for vertical mean subroutine
         integer extract_base ! base codes for vertical mean subroutine
         parameter(extract_base=block_size*0)
         parameter(extract_top=block_size*1)

         integer vert_mean_top ! max code for vertical mean subroutine
         integer vert_mean_base ! base codes for vertical mean subroutin
         parameter(vert_mean_base=block_size*1)
         parameter(vert_mean_top=block_size*2)

         integer zonal_mean_top ! max code for zonal mean subroutine
         integer zonal_mean_base ! base codes for zonal mean subroutine
         parameter(zonal_mean_base=block_size*2)
         parameter(zonal_mean_top=block_size*3)

         integer merid_mean_top ! max code for meridional mean subroutin
         integer merid_mean_base ! base codes for meridional mean subrou
         parameter(merid_mean_base=block_size*3)
         parameter(merid_mean_top=block_size*4)

         integer field_mean_top ! max code for field mean subroutine
         integer field_mean_base ! base codes for field mean subroutine
         parameter(field_mean_base=block_size*4)
         parameter(field_mean_top=block_size*5)

         integer global_mean_top ! max code for global mean subroutine
         integer global_mean_base ! base codes for global mean subroutin
         parameter(global_mean_base=block_size*5)
         parameter(global_mean_top=block_size*6)

CL Weighting

         integer st_weight_code,s_weight ! weighting info address
         parameter(st_weight_code=9,s_weight=9)

         integer stash_weight_null_code,s_noweight ! value of null weigh
         parameter(stash_weight_null_code=0,s_noweight=0)

         integer stash_weight_area_code,s_areaweight ! value of area wei
         parameter(stash_weight_area_code=1,s_areaweight=1)

         integer stash_weight_volume_code,s_volweight
         parameter(stash_weight_volume_code=2,s_volweight=2)

         integer stash_weight_mass_code,s_massweight ! value of mass wei
         parameter(stash_weight_mass_code=3,s_massweight=3)

CL Domain definition

         integer st_north_code,s_north ! northern row address
         parameter(st_north_code=12,s_north=12)

         integer st_south_code,s_south ! southern row address
         parameter(st_south_code=13,s_south =13)

         integer st_west_code,s_west ! western column address
         parameter(st_west_code=14,s_west=14)

         integer st_east_code,s_east ! eastern row address
         parameter(st_east_code=15,s_east =15)

CL Levels

         integer st_input_bottom,s_bottom ! input bottom level address
         parameter(st_input_bottom=10,s_bottom =10)

         integer  st_special_code,s_special ! special code
         parameter(st_special_code=100,s_special=100)

         integer st_input_top,s_top          ! input top level address
         parameter(st_input_top=11,s_top=11)

         integer st_output_bottom,s_outbot   ! output bottom level addre
         parameter(st_output_bottom=21,s_outbot=21)

         integer st_output_top,s_outtop      ! output top level address
         parameter(st_output_top=22,s_outtop=22)

         integer st_model_level_code,s_model
         parameter(st_model_level_code=1,s_model=1)

         integer st_pressure_level_code,s_press ! code for pressure leve
         parameter( st_pressure_level_code=2,s_press=2)

         integer st_height_level_code,s_height ! code for height levels
         parameter(st_height_level_code=3,s_height=3)

         integer st_input_code,s_input               ! input code addres
         parameter(st_input_code=16,s_input=16)

         integer st_input_length,s_length ! input length of diagnostic
         parameter(st_input_length=17,s_length=17)             ! address

         integer st_output_code,s_output ! output code address
         parameter(st_output_code=18,s_output=18)

C Pointer to D1 addressing information
         integer st_position_in_d1,st_d1pos ! Pos of item in D1 for 
         parameter(st_position_in_d1=29,st_d1pos=29) ! relevant submodel

C Output destination options

         integer st_dump,st_secondary
         parameter(st_dump=1,st_secondary=2)

         integer st_output_length,s_outlen ! output length of diagnostic
         parameter(st_output_length=19,s_outlen=19)           ! address
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

         integer st_lookup_ptr       ! ptr to dump lookup header address
         parameter(st_lookup_ptr=23)

         integer st_series_ptr ! ptr into stash_series where control dat
         parameter(st_series_ptr=24)                            ! addres

CL subsid stuff for time series
         integer series_grid_type
         parameter(series_grid_type=1)

         integer series_grid_code
         parameter(series_grid_code=0)

         integer series_long_code
         parameter(series_long_code=1)

         integer series_size
         parameter(series_size=2)

         integer series_proc_code
         parameter(series_proc_code=3)

         integer series_north
         parameter(series_north=4)

         integer series_south
         parameter(series_south=5)

         integer series_west
         parameter(series_west=6)

         integer series_east
         parameter(series_east=7)

         integer series_list_start
         parameter(series_list_start=8)

         integer series_list_end
         parameter(series_list_end=9)

         integer record_size
         parameter(record_size=9)

C Miscellaneous parameters

         integer st_macrotag   ! system/user tag field in stlist address
         parameter(st_macrotag=25)

C Pseudo-level list pointers

         integer st_pseudo_in        ! pseudo-levels input list address
         parameter(st_pseudo_in=26)

         integer st_pseudo_out       ! pseudo-levels output list address
         parameter(st_pseudo_out=27)

C Internal horizontal gridtype codes common to all diagnostics

         integer st_tp_grid,st_uv_grid, ! T-p grid, u-v grid
     &           st_cu_grid,st_cv_grid, ! C-grid (u point, v point)
     &           st_zt_grid,st_zu_grid, ! Zonal T-grid, u-grid
     &           st_mt_grid,st_mu_grid, ! Meridional T-grid, u-grid
     &           st_scalar              ! Scalar (ie. single value)
         parameter(st_tp_grid=1,
     &             st_uv_grid=2,
     &             st_cu_grid=3,
     &             st_cv_grid=4,
     &             st_zt_grid=5,
     &             st_zu_grid=6,
     &             st_mt_grid=7,
     &             st_mu_grid=8,
     &             st_scalar=9)
CLL   Comdeck COASIS : -----------------------------------------------
CLL
CLL   Common declarations and variables of the routines of the oasis
CLL   module.
CLL
CLL   Tested under compiler:   cft77
CLL   Tested under OS version: UNICOS 9.0.4 (C90)
CLL
CLL   Author:   JC Thil.
CLL
CLL   Code version no: 1.0         Date: 18 Nov 1996
CLL
CLL   Model            Modification history from model version 4.1:
CLL   version  date
CLL
CLL
CLL
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered:
CLL
CLL  Project task:
CLL
CLL  External documentation:
CLL

C*********************************************************************
C declatarions for the atmosphere model :
C*********************************************************************
C     This into a common deck to be included in any OASIS routine
C     which requires it.

      integer
     &  G_ROW_LENGTH           ! Global (atmos) row length
     &  ,G_P_ROWS               ! Global (atmos) p  rows
     &  ,G_U_ROWS               ! Global (atmos) uv rows
     &  ,gather_pe              ! Processor for gathering

      integer nfield

      integer nulgr,nulma,nulsu ! unit no of the grid, mask, surf
                                ! files
      character*255 cficgr, cficma, cficsu ! filenames of the grid,
                                ! mask, surf files.
      character*255 coasis_in   ! oasis input filename.

      integer                   ! items of the field locator array.
     &  istash, lon, lat, msk, srf, grd, direction, exc_frequency
     &  , exc_basis
      parameter(
     &  istash = 1, lon = 2, lat = 3, msk = 4, srf = 5
     &  ,grd = 6, direction = 7, exc_frequency = 8
     &  ,exc_basis = 9 )

! Ditto as above, but for the Zinput array :
      integer                   ! items of the field locator array.
     &  Zistash, Zgrd, Zdirection, Zexc_frequency
     &  ,Zexc_basis
      parameter(
     &  Zistash = 1
     &  ,Zgrd = 2, Zdirection = 3, Zexc_frequency = 4
     &  ,Zexc_basis = 5 )

      integer
     &  MaxCouplingField        ! max number of coupling fields.
     &  , NoCouplingField       ! number of coupling fields.
     &  , NbItem                ! Nb of items of the field locator
                                ! array.
     &  , ZNbItem               ! Nb of items of the input file.
      parameter(  MaxCouplingField = 100 ) ! dimension of
                                           ! FieldLocator.
      parameter(  NbItem = 9 )
      parameter( ZNbItem = 5 )

      character*8
     &  ZInput(ZNbItem, MaxCouplingField)
     &  ,FieldLocator(NbItem,MaxCouplingField)

      integer D1_Zptr(MaxCouplingField)
      integer FieldSize(MaxCouplingField) ! size of each coupling
                                          ! field
! Declaration of the pointers on the atmosphere D1.
      integer   ptr_solar, ptr_blue, ptr_longwave, ptr_sensible,
     &  ptr_evap, ptr_snowls, ptr_snowconv, ptr_rainls, ptr_rainconv,
     &  ptr_ice, ptr_pminus, ptr_heat_flux,
     &  ptr_snowfall, ptr_sublimation_accumul, ptr_sublimation_inst,
     &  ptr_slowrunoff, ptr_fastrunoff, ptr_ocentpts,
     &  ptr_runoff


      integer nulou             ! unit for verbose file.

C     cdfile : alias filename for pipe (char string)
C     cdpipe  : symbolic pipe name (char string)
      character*8
     &  cdfile(MaxCouplingField),
     &  cdpipe(MaxCouplingField)
      common / COM_OASIS /
     &  gather_pe, g_row_length, g_p_rows, g_u_rows,
     &  ptr_solar, ptr_blue, ptr_longwave, ptr_sensible,
     &  ptr_evap, ptr_snowls, ptr_snowconv, ptr_rainls, ptr_rainconv,
     &  ptr_ice, ptr_pminus, ptr_heat_flux,
     &  ptr_snowfall, ptr_sublimation_accumul, ptr_sublimation_inst,
     &  ptr_slowrunoff, ptr_fastrunoff,
     &  ptr_runoff, ptr_ocentpts,
     &  NoCouplingField,
     &  FieldLocator,D1_Zptr, FieldSize,
     &  nulou,
     &  nfield,
     &  cdfile, cdpipe

C     end of common deck
C********************************************************************
C     ! The list below describes the status of the atmosphere model
C     ! Exporting and importing fields to the ocean. A symetrical list
C     ! will have to be drawn up for the Ocean part of the UM when
C     ! coupling both UM atmos and Ocean.


      data  cficgr   / "grids" /
      data  cficma   / "masks" /
      data  cficsu   / "areas" /
      data  coasis_in / "namoasis_atm" /

C     end of namelist
C********************************************************************

      integer
     &  ii,j,i                  ! working indexes


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
C
C  Subroutines called
C
      EXTERNAL FINDPTR
!    &  , FINDLOOKPTR
C
C     Local variables
C

      integer
     &  process_code,           ! processing code
     &  freq_code,              ! frequency code
     &  start,end,period,       ! start, end and period step
     &  gridpt_code,weight_code,! gridpt and weighting codes
     &  bottom_level,top_level, ! bottom and top input level
     &  grid_n,grid_s,grid_w,grid_e, ! grid corner definitions
     &  stashmacro_tag          ! stashmacro tag number

      integer
     &  StashCode,              ! integer describing the stash code of
     &  item, section           ! the item code and section code of
                                ! the field.

      integer
     &  im_ident                ! Internal Model Identifier
     &  ,im_index               ! Internal Model Index in Stash arrays

c*-------------------------------------------------------------------
cl-------------------------------------------------------------------
cl 0.  set grid definition information (undefined as search is on
cl     stashmacro tag number)
cl
      process_code=imdi
      freq_code=imdi
      start=imdi
      end=imdi
      period=imdi
      gridpt_code=imdi
      weight_code=imdi
      bottom_level=imdi
      top_level=imdi
      grid_n=imdi
      grid_s=imdi
      grid_e=imdi
      grid_w=imdi

c  set up internal model identifier and stash index
      im_index = internal_model_index(internal_model)

      if (internal_model .eq. atmos_im) then
cl--------------------------------------------------------------------
cl      atmosphere -> ocean (tag=10)
        stashmacro_tag = 10

cl--------------------------------------------------------------------
cl 1.  get address for each field from its stash section/item code
cl     and stashmacro tag if a diagnostic, or from its primary pointer
cl     if prognostic or ancillary field

        im_ident  = internal_model
        im_index  = internal_model_index(im_ident)

        do i = 1, NoCouplingField
C         read the stash code from the  user namelist
          read(FieldLocator(istash, i),'(i8)') StashCode
          section = StashCode / 1000
          item = StashCode - section * 1000

C         Fields in primary space : extracted from pointers set
C         elsewhere into the set_atm_ptr.
          if (StashCode .eq. 00024) then
C         SST:
            D1_Zptr(i) = jtstar
          elseif (StashCode .eq. 00028) then
C           U surface current :
            D1_Zptr(i) = ju_sea
          elseif (StashCode .eq. 00029) then
C           V surface current :
            D1_Zptr(i) = jv_sea
          elseif (StashCode .eq. 00032) then
C           Ice depth :
            D1_Zptr(i) = jice_thickness
          elseif (StashCode .eq. 00023) then
C           Snow depth :
            D1_Zptr(i) = jsnodep
          elseif (StashCode .eq. 00031) then
C           Seaice fraction :
!*IF DEF,SEAICE
            D1_Zptr(i) = JICE_FRACTION
            ptr_ice    = JICE_FRACTION
!*ENDIF
          else
C           Fields in secondary space have their pointers
C           extracted with findptr :
            call findptr(internal_model, section, item,
     &        process_code,freq_code,start,end,period,
     &        gridpt_code,weight_code,
     &        bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &        stashmacro_tag,imdi,D1_Zptr(i),
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &        icode,cmessage )
            if (icode.lt.0) goto 999
          endif
        enddo

!  Setup of the dedicated coupling fields to couple the atmosphere
!  with the current UM-ocean.
C  Net integrated downward solar on atmos grid
        call findptr(internal_model, 1, 203,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_solar,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Net downward blueband solar on atmos grid
        call findptr(internal_model, 1, 204,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_blue,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Net downward longwave on atmos grid
        call findptr(internal_model, 2, 203,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_longwave,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Sensible heat on atmos grid, area mean over open sea
        call findptr(internal_model, 3, 228,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_sensible,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Surface evaporation over sea weighted by fractional leads
        call findptr(internal_model, 3, 232,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_evap,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Large-scale snowfall rate on atmos grid
C       Convective snowfall rate on atmos grid
        call findptr(internal_model, 4, 204,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_snowls,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Convective snowfall rate on atmos grid
        call findptr(internal_model, 5, 206,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_snowconv,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Large-scale rainfall rate on atmos grid
        call findptr(internal_model, 4, 203,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_rainls,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Convective rainfall rate on atmos grid
        call findptr(internal_model, 5, 205,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_rainconv,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       SLOW runoff on atmos grid
        call findptr(internal_model, 8, 205,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_slowrunoff,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       FAST runoff on atmos grid
        call findptr(internal_model, 8, 204,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_fastrunoff,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999
C       Sublimation in atm D1 :
        call findptr(internal_model, 3, 231,
     &    process_code,freq_code,start,end,period,
     &    gridpt_code,weight_code,
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,
     &    stashmacro_tag,imdi,ptr_sublimation_accumul,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
C Applicable to all configurations
C STASH related variables for describing output requests and space
C management.
! vn3.5 (Apr. 95)  Sub-models project   S.J.Swarbrick
!                  PPXREF, INDEX_PPXREF removed

     &       SF, STINDEX, STLIST, SI, STTABL, STASH_MAXLEN,
     &       PPINDEX,  STASH_LEVELS, STASH_PSEUDO_LEVELS,
     &       STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,


     &    icode,cmessage )
        if (icode.lt.0) goto 999

C       Runoff coastal outflow point (pointer):
        ptr_ocentpts = si(93,0,im_index)

      elseif (internal_model .eq. ocean_im) then
      endif

 999  continue
      if(icode.ne.0) then
        write(nulou,*) cmessage,icode
      endif
      return
      end
