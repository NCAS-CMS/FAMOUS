C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! Calls routines to initialize veg parameters and accumulated C fluxes
!
! Subroutine Interface:
      SUBROUTINE INIT_VEG(A_STEP,
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
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
C Constants for dynamics routines
     &    AK_TO_THE_KAPPA, BK_TO_THE_KAPPA, AKH, BKH, AKH_TO_THE_KAPPA,
     &    BKH_TO_THE_KAPPA, COS_U_LATITUDE, COS_P_LATITUDE,
     &    SEC_U_LATITUDE, SEC_P_LATITUDE, TAN_U_LATITUDE, SIN_LONGITUDE,
     &    COS_LONGITUDE, TRUE_LONGITUDE, F1, F2, F3, F3_P,
     &    TRIGS, TWO_D_GRID_CORRECTION, ETA_MATRIX_INV,
C Constants for physics routines
     &    SOILB, LAND_LIST,
     &    SIN_U_LATITUDE,
     &              ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Initializes vegetation parameters from fractions of surface types
!   and initializes accumulated carbon fluxes to zero if a new TRIFFID
!   calling period is starting.
!
! Method:
!   Calls routine SPARM to initialize vegetation parameters.
!   Calls routine INIT_ACC to initialize accumulated carbon fluxes.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4   10/10/97   Original code.  Richard Betts
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

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
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
CL This COMDECK needs COMDECK TYPSIZE *CALLed first
CL This COMDECK needs COMDECK CMAXSIZE *CALLed first
CL                           to be called in the same module.
CL CMAXSIZE should be called first.
CL
C Constants for ATMOSPHERE model.
C Constants for routines independent of resolution.
C*L================ COMDECK CCONSTS ===========================
C   Description:
C     This COMDECK contains declarations for derived constants within
C   the atmospheric model. Where necessary PARAMETERS are defined to
C   dimension these constants. All constants are placed in the common
C   block CDERIVED.
C   COMDECK CMAXSIZE must be called first
C
C   The derived constants are calculated in the routine SETCONA1.
C
CLL PA, WI      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.2   26/03/93  Remove resolution dependent variables for dynamic
CLL                 allocation. R.Rawlins
!   4.0   20/07/95  Sizes of tables expanded for version 3A
!                   of the radiation in sections 1 and 2. J.M. Edwards
C Logical component: F011
C

C Define Parameters:
      INTEGER  NUM_CLOUD_TYPES    ! No of cloud types ie low/med/high
        PARAMETER (NUM_CLOUD_TYPES = 3)
      INTEGER  LEN_SW_TABLES       ! Size of look-up tables for SW
        PARAMETER (LEN_SW_TABLES = 1)                           
!          This declaration is needed when the comdeck is
!          used outside section 1.
      INTEGER  LEN_LW_TABLES       ! Size of look up tables for LW
C                                  ! radiation.
C
      INTEGER  MATRIX_POLY_ORDER      ! Order of polynomial used in
      PARAMETER (MATRIX_POLY_ORDER=5) ! calculation
C
        PARAMETER (LEN_LW_TABLES = 1 )
C     !  If there are no LW look-up tables - as in the ECMWF-like
C     !    treatment of gaseous LW effects, or if the LW is being
C     !    completely omitted - the definitions are left in the code
C     !    though the space is not needed, so the size is best set to 1.

C Declare derived constants:
      INTEGER  LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER  LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER  MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER  MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER  HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER  HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      REAL     SW_TABLES          ! Tables for SW radiation
      REAL     LW_TABLES          ! Tables for LW radiation
      REAL     ETA_SPLIT          ! Eta values of layer bndry cld types
C                                 ! (high/med/low

      LOGICAL  ELF                ! T if atmosphere model on LAM grid
      INTEGER  IFAX               ! factor for FFT (10)
C Constants for dynamics routines independent of resolution but
C  dependent on choice of dynamics parameter.
      INTEGER ADJ_TIME_SMOOTHING_WEIGHT
      REAL    ADJ_TIME_SMOOTHING_COEFF
C Constants for dynamics output independent of resolution but
C  dependent on choice of levels for output.
      REAL    REQ_THETA_PV_LEVS
C*----------------------------------------------------------------------
      COMMON /CDERIVED/
     &  SW_TABLES(LEN_SW_TABLES),      LW_TABLES(LEN_LW_TABLES),
     &  ETA_SPLIT(NUM_CLOUD_TYPES+1),  LOW_BOT_LEVEL,  LOW_TOP_LEVEL,
     &  MED_BOT_LEVEL, MED_TOP_LEVEL,  HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,
     &  ADJ_TIME_SMOOTHING_WEIGHT(MAX_ADJ_TSL),ADJ_TIME_SMOOTHING_COEFF,
     &  ELF,IFAX(10),
     &  REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)
C
C Derived for dynamics routines. Array sizes are resolution dependent.
C                                     ! of ETA_HALF inverse matrix.
      REAL AK_TO_THE_KAPPA(P_LEVELS)  ! (a/pref)**(r/cp) at full levels
      REAL BK_TO_THE_KAPPA(P_LEVELS)  ! (b/pref)**(r/cp) at full levels
      REAL AKH(P_LEVELS+1)            !  value of a at half levels
      REAL BKH(P_LEVELS+1)            !  value of b at half levels
      REAL AKH_TO_THE_KAPPA(P_LEVELS+1) ! (a/pref)**(r/cp) at hf levels
      REAL BKH_TO_THE_KAPPA(P_LEVELS+1) ! (b/pref)**(r/cp) at hf levels
      REAL COS_U_LATITUDE(U_FIELD)    ! cos(lat) at u points (2-d array)
      REAL COS_P_LATITUDE(P_FIELD)    ! cos(lat) at p points (2-d array)
      REAL SEC_U_LATITUDE(U_FIELD)    ! 1/cos(lat) at u points (2-d )
      REAL SEC_P_LATITUDE(P_FIELD)    ! 1/cos(lat) at p points (2-d )
      REAL SIN_U_LATITUDE(U_FIELD)    ! sin(lat) at u points (2-d array)
      REAL TAN_U_LATITUDE(U_FIELD)    ! tan(lat) at u points (2-d array)
      REAL SIN_LONGITUDE(ROW_LENGTH)  ! sin(long)
      REAL COS_LONGITUDE(ROW_LENGTH)  ! cos(long)
      REAL TRUE_LONGITUDE(P_FIELD)    ! true longitude of point
      REAL F1(U_FIELD),F2(U_FIELD),F3(U_FIELD)   ! Coriolis terms on
C                          ! on u points. These are defined in UMDP 10.
      REAL F3_P(P_FIELD)              ! Coriolis term F3 on p points
      REAL TRIGS(ROW_LENGTH)          ! sin and cos for FFT routines

      REAL TWO_D_GRID_CORRECTION(P_FIELD/ROW_LENGTH) ! factor used in
C                                     ! filtering see UMDP 10.
      REAL ETA_MATRIX_INV (MATRIX_POLY_ORDER,MATRIX_POLY_ORDER,P_LEVELS)
C Constants for physics  routines. Array sizes are resolution dependent.
      REAL    SOILB(LAND_FIELD)       ! SOILB parameter
      INTEGER LAND_LIST(LAND_FIELD)   ! List of land points
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
CLL  Comdeck: CCONTROL -------------------------------------------------
CLL
CLL  Purpose: COMMON block for top level switches and 2nd level switches
CLL           needed by the top level (C0) and 2nd level routines, but
CLL           not held in the history COMMON block.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1    8/02/93 : Changed 99 to NUNITS for i/o. Note this comdeck
CLL                    must always be called after CHSUNITS so that
CLL                    NUNITS is defined.
CLL  3.1  15/02/93  Add L_Z0_OROG orographic roughness switch. R.Barnes.
CLL  3.3  09/07/93  Add L_CONVECT =F to add saved convection increments,
CLL                               =T to call conv.scheme.  R.T.H.Barnes.
CLL   3.3    13/12/93   Insert switches for half timestep
CLL                     dynamics. A.S.Lawless
CLL   3.4 23/08/94  Add switch for local -ve q correction R.A.Stratton.
CLL   3.4    16/06/94   COMMON block DEFLOGIC inserted - declares
CLL                     logical switches for control and other
CLL                     purposes - most of these have replaced *DEFs
CLL                                                     S.J.Swarbrick
CLL   3.4    1/8/94  Add control for assimilation mode S Bell
CLL  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
CLL                 contents. Control expanded and subdivided by
CLL                 overall, generic and specific categories. RTHBarnes
CCL  4.1  23/02/96  Extend for new wave sub-model.  RTHBarnes.
CLL
CLL Logical components covered :
CLL
CLL External documentation: Unified Model documentation paper No
CLL                         Version
CLL
CLLEND ---------------------------------------------------------------

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


! ----------------------- Comdeck: CNTLATM  ----------------------------
! Description: COMDECK defining Control variables for the Atmosphere
!              internal model, and its runtime constants.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  17/08/95  New variables for long physics timestep,
!                 and tracer advection.  RTHBarnes.
!  4.0  7/11/95  Logical switches for convective momentum transports
!                and CAPE closure added to namelists. Pete Inness.
!  4.1  8/5/96   Logical switch for rapidly mixing boundary layer
!  4.1 28/05/96  New control switches added for Sulphur Chemistry
!                Cycle and Hydrology Schemes. D Robinson & D. Goddard.
!  4.3  18/3/97  Flag to indicate if the HadCM2 approximate treatment
!                of sulphate aerosol is being used.       William Ingram
!  4.3 14/04/97  New control switch L_OLD_PWTS for old polar geometric
!                weights, needed for HADCM2.    T Johns.
!  4.3 03/02/97  Logical L_MIXLEN for mixing in the boundary layer
!                  S Jackson
!  4.3 03/02/97  Logical switches L_XSCOMP and L_SDXS for convection
!                scheme  S Jackson
!  4.4 4/7/97    Add control for IAU  Chris Jones/Stuart Bell 
!  4.4 05/09/97  Logical LFLUX_RESET to indicate when net flux field
!                needs initialising to 0. S.D. Mullerworth
!  4.4 17/09/97  Logical switches L_CCW and L_3D_CCA added to enable
!                use of anvil package/3D conv. cloud amt. J.Gregory
!  4.4 08/09/97 Logical switches L_LSPICE, L_BL_LSPICE and L_LSPICE_BDY
!               for mixed phase precipitation.
!                                                       D.Wilson
!  4.4 10/10/97  Logical switch L_SNOW_ALBEDO.  Richard Essery   
!  4.4 10/10/97  Logical switches L_VEG_FRACS, L_TRIFFID, L_PHENOL, 
!                L_NRUN_MID_TRIF and L_TRIF_EQ for veg.  Richard Betts
!  4.5   1/07/98  Add logicals to control interactive CO2 use. C.D.Jones
!   4.5  28/04/98  Add logicals for NH3 and SOOT variables and emiss
!                                                           M Woodage
!  4.5 21/08/98  Logical switch l_ssice_albedo.  Jonathan Gregory
!  4.5 20/05/98  Logical switch L_NEG_TSTAR for negative surface 
!                temperature error check.  Richard Betts  
!  4.5 19/11/98  Add PHENOL_PERIOD and TRIFFID_PERIOD, moved from
!                NLSTCATM.  Richard Betts  
!  4.5 19/05/98  Logical switch L_PHASE_LIM for HADAM3 physics in
!                optimised convection scheme.       Julie Gregory
!  4.5 26/06/98  Logical switches L_RHCPT, L_CLD_AREA for new
!                Section 9 parametrizations.             S. Cusack
!  4.5 22/10/98  Remove redundant switch LMULTIL_HYDROL
!                Author D.M. Goddard
!  4.5 05/06/98  Add Logical switch L_VINT_TP.  D Robinson.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Parameter declarations
      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
!
!   Type declarations
!
      INTEGER H_SWBANDS,      ! Number of shortwave radiation bands
     &        H_LWBANDS,      ! Number of longwave radiation bands
     &        A_ADJSTEPS,         ! No. of adjustment timesteps per
!                                 ! advection step
     &        A_SW_RADSTEP,       ! Number of advection steps per
!                                 ! shortwave radiation step
     &        A_LW_RADSTEP,       ! Number of advection steps per
!                                 ! longwave radiation step
     &        A_SW_SEGMENTS,      ! No of batches used in shortwave code
     &        A_LW_SEGMENTS,      ! No of batches used in longwave code
     &        A_CONV_STEP,        ! No of advection timesteps between
!                                 ! calls to convection scheme
     &        A_CONVECT_SEGMENTS, ! No of batches in convection code
     &        A_NSET_FILTER,      ! No of advection steps after which
!                                 ! filtering wavenumber checked
     &        A_ENERGYSTEPS,      ! Number of advection steps after
!                                 ! which energy adjustment performed
     &        A_ASSIM_START_HR,   ! Time at which data assimilation
!                                 ! starts (Hours after Basis Time)
     &        A_ASSIM_END_HR      ! Time at which data assimilation
!                                 ! ends (Hours after Basis Time)
     &       ,T_IAU_START, T_IAU_END  ! IAU before and after
!
     &       ,A_SWEEPS_DYN ! No.of dynamics sweeps per physics timestep
     &       ,CALL_CHEM_FREQ ! Frequency of calls to CHEM_CTL
     &       ,PHENOL_PERIOD ! Update frequency for leaf phenology (days)
     &       ,TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      LOGICAL
     1       L_SW_RADIATE,           ! Activate SW radiation
     2       L_LW_RADIATE,           ! Activate LW radiation
     3       LADD_RADINCS,           ! Both SW and LW radiation active
     &       L_H2_SULPH,             ! HadCM2 approximate sulphate on
     4       L_SET_FILTER,           ! Recalculate filtering in dynamics
     5       LDAY,                   ! End-of-day
     6       LEXPAND_OZONE,          ! Convert zonal mean ozone to field
     6       L_COMPRESS_LAND,        ! Compress land points in physics
     *       L_NEG_THETA,            ! Test for -ve theta in dynamics
     *       L_NEG_PSTAR,            ! Test for -ve P* in dynamics
     *       L_NEG_QT,               ! Test for -ve QT over layer
     &       L_NEG_TSTAR,           ! Test for -ve surface temperature.
     *       L_FIELD_FLT,            ! Apply field filtering (atmos)
     *       L_SUPERBEE,             ! Superbee(T) or Van Leer(F)
     *                               ! limiter for tracer advection
     7       LENERGY,                ! Recalculate energy drift
     &       LFLUX_RESET,            ! Reset net flux field
     *       L_Z0_OROG,              ! T to use orog.roughness code
     *       L_RMBL,                 ! T to use rapid mixing BL code
     &       L_MIXLEN,               ! T to make mixing length above BL
     &                               ! top independent of BL depth
     *       L_CONVECT               ! T call conv.scheme, F add incrs.
     *      ,L_HALF_TIMESTEP_DYN     ! T if wind threshold exceeded
     *      ,L_HALF_TIMESTEP_DIV     ! T if diverg. threshold exceeded
     &      ,L_QT_POS_LOCAL          ! Apply -ve q correction locally.
     &      ,L_TRACER_THETAL_QT      ! T if using tracer advection
!                                    !  for thetal & qt
     &      ,L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR !Switches for assm mode
     &      ,L_OLD_PWTS              ! T if using old polar weights
     &      ,L_VINT_TP               ! T: Use V_INT_TP to output Temp
                                     ! on model levels.
C
      LOGICAL                       ! Logical switches for:
     &   LFROUDE        ,           !  Limit max grav wave amp
     &   LGWLINP        ,           !  Linear grav wave stress prof
     &   LLINTS         ,           !  Linear TS approx
     &   LWHITBROM      ,           !  White & Bromley terms
     &   LEMCORR        ,           !  Energy & mass correction
     &   LMICROPHY      ,           !  Microphysics in sw rad scheme
     &   L_MURK         , !           :Total aerosol field
     &   L_MURK_ADVECT  , !           :Aerosol advection
     &   L_MURK_SOURCE  , !Bndry      :Aerosol source & sink terms
     &   L_MURK_BDRY    , !Layer      :UK Mes bndry model
     &   L_BL_TRACER_MIX, !model      :Bndry layer tracer mixing
     &   L_MOM,                     !  convective momentum mixing
     &   L_CAPE,                    !  CAPE closure for convection      
     &   L_SDXS,                    ! Convective excess from turbulent
!                                   !  fluctuations
     &   L_XSCOMP                   ! Environmental compensation for 
!                                   !  parcel excess
     &  ,L_3D_CCA                   ! Use 3D conv cloud amount
     &  ,L_CCW                      ! Rain not inc. in conv water path
     &  ,L_PHASE_LIM                ! Select 3B physics for A05_3C
     &  ,L_CLOUD_DEEP               ! Depth criterion applied for anvils
     &  ,LSINGLE_HYDROL             ! Single level hydrology
     &  ,LMOSES                     ! MOSES hydrology only
     &  ,L_SNOW_ALBEDO              ! Prognostic snow albedo
     &  ,l_ssice_albedo             ! Sea-ice albedo affected by snow
     &  ,L_SULPC_SO2                ! Sulphur Cycle : SO2 MMR included 
     &  ,L_SULPC_DMS                ! Sulphur Cycle : DMS MMR included
     &  ,L_SULPC_OZONE              ! Sulphur Cycle : Ozone included
     &  ,L_SO2_SURFEM               ! SO2 Surface Emissions
     &  ,L_SO2_HILEM                ! SO2 High Level Emissions
     &  ,L_SO2_NATEM                ! SO2 Natural Emissions
     &  ,L_DMS_EM                   ! DMS Emissions
     &  ,L_SULPC_NH3           ! S Cycle : NH3 included
     &  ,L_NH3_EM              ! S Cycle : NH3 emiss included
     &  ,L_SOOT                ! Soot included  
     &  ,L_SOOT_SUREM          ! surface Soot emiss included
     &  ,L_SOOT_HILEM          ! elevated Soot emiss included
     &  ,L_USE_SOOT_DIRECT     ! direct radiative effects of soot
     &  ,L_USE_SULPC_DIRECT   !\Use SO4 aerosol from sulphur cycle for
     &  ,L_USE_SULPC_INDIRECT_SW !direct/indirect effect in radiation,
     &  ,L_USE_SULPC_INDIRECT_LW !the latter for both SW and LW.
     &  ,L_CLIMAT_AEROSOL           ! Switch for climatological
!                                   ! aerosols in the radiation.        
     &  ,L_RHCPT                     ! controls the use of new RHcrit
                                     ! parametrization, vn 2B of Sec 9
     &  ,L_CLD_AREA                  ! controls cloud area parametriz.
     &  ,L_IAU_DIAG                 ! controls IAU diagnostics
     &  ,L_IAU                      ! controls IAU calls
     &  ,L_IAU_RAMP                 ! controls IAU weights
     &  ,L_VEG_FRACS                ! Switch for vegetation fractions
     &  ,L_TRIFFID                  ! Switch for interactive veg model
     &  ,L_PHENOL                   ! Switch for leaf phenology
     &  ,L_NRUN_MID_TRIF            ! Switch for starting NRUN mid-way 
C                                   ! through a TRIFFID calling period
     &  ,L_TRIF_EQ                  ! Switch for running TRIFFID in
C                                   ! equilibrium mode
     &      ,L_CO2_INTERACTIVE      ! interactive 3D CO2 field for
                                    !  use with carbon cycle model
     &      ,L_CO2_EMITS            ! include surface emissions
      LOGICAL L_LSPICE              ! New cloud/precip microphysics
     &,       L_BL_LSPICE           ! Full boundary layer treatment
!                                     with new cloud/precip scheme
     &,       L_LSPICE_BDY          ! QCF present in lateral boundaries
!
      CHARACTER*5 A_ASSIM_MODE     ! Switch for BG/AC/3DVAR/4DVAR assm
!
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section versions
!
      NAMELIST / NLSTCATM /
     & L_VEG_FRACS, L_TRIFFID, L_PHENOL, L_NRUN_MID_TRIF, L_TRIF_EQ,
     & PHENOL_PERIOD, TRIFFID_PERIOD,
     & H_SWBANDS, H_LWBANDS, A_ADJSTEPS,
     & A_SW_RADSTEP, A_LW_RADSTEP, A_SW_SEGMENTS, A_LW_SEGMENTS,
     & A_CONV_STEP, A_CONVECT_SEGMENTS, A_NSET_FILTER, A_ENERGYSTEPS,
     & A_ASSIM_START_HR, A_ASSIM_END_HR, A_SWEEPS_DYN, L_H2_SULPH,
     & L_SW_RADIATE, L_LW_RADIATE, LADD_RADINCS, L_SET_FILTER,
     & LDAY, LEXPAND_OZONE, L_COMPRESS_LAND,
     & L_NEG_THETA, L_NEG_PSTAR, L_NEG_QT, L_NEG_TSTAR, L_FIELD_FLT,
     & L_SUPERBEE, LENERGY, L_Z0_OROG, L_RMBL, L_MIXLEN, L_CONVECT,
     & L_HALF_TIMESTEP_DYN, L_HALF_TIMESTEP_DIV, L_QT_POS_LOCAL,
     & L_TRACER_THETAL_QT,LFLUX_RESET,
!    & L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR, in COMMON but not in NAMELIST
     & L_OLD_PWTS, L_VINT_TP,
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR,
     & LMICROPHY, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE,
     & L_MURK_BDRY, L_BL_TRACER_MIX, L_MOM, L_CAPE, L_SDXS, L_XSCOMP,
     & L_3D_CCA, L_CCW, L_PHASE_LIM, L_CLOUD_DEEP,
     & LSINGLE_HYDROL, LMOSES,
     & L_CO2_INTERACTIVE, L_CO2_EMITS,
     & L_SNOW_ALBEDO,
     & l_ssice_albedo,
     &  L_CLIMAT_AEROSOL,
     &  L_RHCPT,L_CLD_AREA,
     & L_LSPICE,L_BL_LSPICE,L_LSPICE_BDY,
     & L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE,
     & L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM,
     & L_SULPC_NH3,L_NH3_EM,L_SOOT,L_SOOT_SUREM,L_SOOT_HILEM,
     & L_USE_SOOT_DIRECT,
     & L_IAU,L_IAU_RAMP,L_IAU_DIAG,
     & T_IAU_START, T_IAU_END,
     & L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT_SW,
     & L_USE_SULPC_INDIRECT_LW, CALL_CHEM_FREQ,
     & A_ASSIM_MODE, H_SECT

      COMMON / CNTLCATM /
     & L_VEG_FRACS, L_TRIFFID, L_PHENOL, L_NRUN_MID_TRIF, L_TRIF_EQ,
     & PHENOL_PERIOD, TRIFFID_PERIOD,
     & H_SWBANDS, H_LWBANDS, A_ADJSTEPS, L_H2_SULPH,
     & A_SW_RADSTEP, A_LW_RADSTEP, A_SW_SEGMENTS, A_LW_SEGMENTS,
     & A_CONV_STEP, A_CONVECT_SEGMENTS, A_NSET_FILTER, A_ENERGYSTEPS,
     & A_ASSIM_START_HR, A_ASSIM_END_HR, A_SWEEPS_DYN,
     & L_SW_RADIATE, L_LW_RADIATE, LADD_RADINCS, L_SET_FILTER,
     & LDAY, LEXPAND_OZONE, L_COMPRESS_LAND,
     & L_NEG_THETA, L_NEG_PSTAR, L_NEG_QT, L_NEG_TSTAR, L_FIELD_FLT,
     & L_SUPERBEE, LENERGY, L_Z0_OROG, L_RMBL, L_MIXLEN, L_CONVECT,
     & L_HALF_TIMESTEP_DYN, L_HALF_TIMESTEP_DIV, L_QT_POS_LOCAL,
     & L_TRACER_THETAL_QT,LFLUX_RESET,
     & L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR,
     & L_OLD_PWTS, L_VINT_TP,
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR,
     & LMICROPHY, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE,
     & L_MURK_BDRY, L_BL_TRACER_MIX, L_MOM, L_CAPE, L_SDXS, L_XSCOMP,
     & L_3D_CCA, L_CCW, L_PHASE_LIM, L_CLOUD_DEEP,
     & LSINGLE_HYDROL, LMOSES,
     & L_CO2_INTERACTIVE, L_CO2_EMITS,
     & L_SNOW_ALBEDO,
     & l_ssice_albedo,
     &  L_CLIMAT_AEROSOL,
     &  L_RHCPT,L_CLD_AREA,
     & L_LSPICE,L_BL_LSPICE,L_LSPICE_BDY,
     & L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE,
     & L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM,
     & L_SULPC_NH3,L_NH3_EM,L_SOOT,L_SOOT_SUREM,L_SOOT_HILEM,
     & L_USE_SOOT_DIRECT,

     & L_IAU,L_IAU_RAMP,L_IAU_DIAG,
     & T_IAU_START, T_IAU_END,
     & L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT_SW,
     & L_USE_SULPC_INDIRECT_LW, CALL_CHEM_FREQ,
     & A_ASSIM_MODE, H_SECT
! ----------------------- Comdeck: CNTLOCN  ----------------------------
! Description: COMDECK defining Control variables for the Ocean
!              internal model.
!   This comdeck contains logical variables which are used on the
!   control of certain sections of Ocean model code
!   They replace the previous method of controlling code using *IF DEFs.
!
! Author : R.T.H.Barnes & R.Hill
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  29/05/96  include L_OZVRT     M. J. Bell 
!  4.3  8.11.96   include L_SLOPEMAX and L_COXCNVC  JMG
!  4.4  11.08.97   Remove L_OCHEQUB.    R. Hill 
!    4.4  10/09/97  Remove all references to SKIPLAND code. R.Hill
!  4.4  8.07.97   include L_FLUXD      R.Lenton
!  4.5  3.11.98    include L_OBIMOM       M. Roberts
!  4.5  3.11.98   include L_OMEDADV and L_OHUDOUT  (new outflow param)
!                 M. Roberts
!  4.5  10.11.98  New logicals: L_OISOMOM, L_OISOGM L_OISOGMSKEW
!                 L_OBIHARMGM and L_OVISHADCM4
!  4.5  3.9.98    Changes for HADCM4 sea-ice. Cresswell and Gregory
!  4.5   1/07/98  Add logical to control interactive CO2 use. C.D.Jones
CLL   4.5 G.J.Rickard include L_OFULARGE (full Large scheme),
CLL                   L_OPANDP (choice of vertical mixing),
CLL                   L_OSTATEC (density calculation choice),
CLL                   L_OUSTARWME (ustar calculation).
CLL
!  4.5  7.8.97    Removed old ocean boundary logicals L_OBGILLS,
!                 L_OBGILLN, L_OSTEVNS, L_OSTEVS and L_BOUNDSO.  
!                 Added in new logicals L_OBDY_NORTH to L_OBDY_STREAM.
!                 M.J. Bell
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     &        O_CLM_START_HR,     ! Time ocean climate increments start
     &        O_CLM_END_HR,       ! Time ocean climate increments end
     &        O_INT_CLM_INC,      ! # ocean steps  } climate incs.
     &        O_INT_ANA_STP,      ! # between      } analysis steps
!
     &        O_INT_EVO_BTS,      ! # ocean steps between fwd evolution
!                                     of bathys and tesacs
     &        O_INT_VRY_BTS,      ! # ocean steps between re-calculation
!                         of future bathys and tesacs valid at this hour
     &        O_INT_WTS_ACC,      ! # ocean steps betwn accumulating wts
!
     &        O_INT_OBS_FRSH,     ! # ocean  } reading new OBS files
     &        O_INT_OBS_OUT,      ! # steps  } outputting new OBS files
     &        O_INT_OBS_STR,      ! # between} caching OBS array
     &        O_INT_FLD_STR,      ! #        } caching model fields
     &        O_ASSIM_START_HR,   ! Time at which data assimilation
!                                 ! starts (Hours after Basis Time)
     &        O_ASSIM_END_HR      ! Time at which data assimilation
!                                 ! ends (Hours after Basis Time)
!
      LOGICAL L_FLUXCORR   ! Heat & water flux correction
     &,       L_OGLOBAL    ! Global ocean
     &,       L_ICEFREEDR  ! Free Drift Sea Ice model
     &,       L_ICESIMPLE  ! Simple Advection Sea Ice model
     &,       L_HADCM4O2I  ! HADCM4 version of ocean-to-ice heat flux
     &,       L_IHANEY     ! Haney Forcing Ice
     &,       L_OADGHR2    ! Ocean assimilation diagnostics
     &,       L_OBDY_NORTH   ! Update northern lateral boundary   
     &,       L_OBDY_SOUTH   ! Update southern lateral boundary
     &,       L_OBDY_EAST    ! Update eastern lateral boundary
     &,       L_OBDY_WEST    ! Update western lateral boundary
     &,       L_OGILL_LBCS   ! Use the Gill boundary scheme
     &,       L_OFRS_LBCS    ! Use the FRS boundary scheme
     &,       L_OSTVNS_LBCS  ! Use the Stevens boundary scheme
     &,       L_OBDY_TRACER  ! Update the tracers
     &,       L_OBDY_UV      ! Update the velocities
     &,       L_OBDY_STREAM  ! Update the stream functions
     &,       L_OBDY_ICE     ! Update ice fields (snow, aice, hice)
     &,       L_OBIOLOGY   ! Effect of phytoplankton on carbon cycle
     &,       L_OCARB14    ! Calculate atmospheric C12/C14 ratio
     &,       L_OCARBON    ! Carbon cycle model
     &,       L_OCNASSM    ! Activate ocean assimilation
     &,       L_OCYCLIC    ! Cyclic boundary conditions
     &,       L_OFILTER    ! Fourier filtering for high latitudes
     &,       L_OFREESFC   ! Use free surface conditions
     &,       L_FLUXD
     &,       L_OHANEY     ! Haney Forcing heat/fresh water fluxes
     &,       L_OHMEAD     ! Mead tracer transport diagnostics
     &,       L_OICECOUP   ! Coupled model with Sea Ice
     &,       L_OIMPADDF   ! Crank-Nicholson vert. advn-diffn scheme
     &,       L_OIMPDIF    ! CN vertical diffusion scheme
     &,       L_OISLANDS   ! Include Island Routines
     &,       L_OISOPYC    ! Isopycnal diffusion scheme
     &,       L_OLATVISC   ! Latitude dependent viscosity
     &,       L_OLISS      ! Liss & Merlivat wind mixing of tracers
     &,       L_OMIXLAY    ! Wind mixing of tracers-mixed layer scheme
     &,       L_ONOCLIN    ! Barotropic solution
     &,       L_ONOPOLO    ! No sea ice at North Pole
     &,       L_OPENBC     ! Read in lateral boundary fields
     &,       L_OPSEUDIC   ! Pseudo-ice routine
     &,       L_ORICHARD   ! Evaluate & use Richardson No.
     &,       L_OROTATE    ! Coriolis force calculation
     &,       L_OSOLAR     ! Calc solar penetration for given water type
     &,       L_OSOLARAL   ! Calc sol. pen. - simplified layer structure
     &,       L_OSYMM      ! Symmetric boundary conditions
     &,       L_OVARYT     ! Varying time step with depth
     &,       L_RIVERS     ! River run-off routines
     &,       L_SEAICE     ! Include Sea Ice model
     &,       L_TRANGRID   ! Spatial interp. in coupled model
     &,       L_OCONJ     ! Whether to use conjugate gradient solver
     &,       L_UPWIND     ! Upwind differencing for tracer advection
     &,       L_OPRINT     ! Whether to print incidental ocean info
     &,       L_OSTVEW     !\
     &,       L_OPMSL      ! \
     &,       L_OTIDAL     !  \___ All for use with O. Alves free
     &,       L_OFOURW     !  /    surface modifications at V4.0
     &,       L_ODELPLUS   ! /
     &,       L_OTROPIC    !/
     &,       L_OISOMOM
     &,       L_OISOGMSKEW
     &,       L_OISOGM
     &,       L_OBIHARMGM
     &,       L_OVISHADCM4
     &,       L_OMEDOUT    ! Mediterranean outflow - 288*144 and 96*73
                           !   grids only - uses hardwired gridpoint nos
     &,       L_OCONVROUS  ! Roussenov convective adjustment
     &,       L_OEXTRAP    ! Extrapolation of vertical density gradients
     &,       L_OISOPYCGM  ! Gent and McWilliams eddy parametrisation.
     &,       L_OISOTAPER  ! Tapering of isopycnal diffusion
     &,       L_OVISBECK    ! Visbeck scheme
     &,       L_OQLARGE     ! Quadratic Large scheme
     &,       L_OFULARGE   ! FULL LARGE SCHEME
     &,       L_OPANDP     ! RI-DEPENDENT VERT MIX SCHEMES
     &,       L_OSTATEC    ! DENSITY CHOICE FOR RI-CALC
     &,       L_OUSTARWME  ! WME OR WSTRESS TO FIND USTAR

     &,       L_OZVRT      ! barotropic vorticity diagnostic switch
                           ! set by OCN_FOR_STEP (not in namelist) 
     &,       L_SLOPEMAX   ! Selects SLOPE_MAX isopycnal diffusion
     &,       L_COXCNVC    ! Selects original Cox convection scheme
     &,       L_OCOMP    ! Land pnts compressed from dump (3d fields)
     &,       L_OMEDADV
     &,       L_OHUDOUT
     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA
     &      ,L_CO2O_INTERACTIVE     ! interactive 3D CO2 field for
                                    !  use with carbon cycle model
     &,       L_OBIMOM  ! biharmonic momentum diffusion
! *IF DEF,OCNASSM
!    Additions to CCONTROL for ocean assimilation
!
      LOGICAL
     &       LAS_CLM_INC,    ! make increments to relax to climate
     &       LAS_ADD_INC,    ! add or subtract analysis increments
     &       LAS_ANA_STP,    ! calculate analysis increments
     &       LAS_EVO_BTS,    ! evolve bathy and tesac obs 1 step
     &       LAS_VRY_BTS,    ! estimate bathys and tesacs at this hour
     &       LAS_WTS_ACC,    ! evolve accumulated weights
     &       LAS_OBS_FRSH,   ! to refresh main OBS data set
     &       LAS_OBS_OUT,    ! output ACOBS file for incremented obs
     &       LAS_FLD_STR,    ! output model fields to cache store
     &       LAS_OBS_STR     ! output obs to cache store
! *ENDIF OCNASSM


      NAMELIST / NLSTCOCN /
     & O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,
     & O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,
     & O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,
     & O_ASSIM_START_HR, O_ASSIM_END_HR, L_FLUXCORR, L_OGLOBAL,
     & L_ICEFREEDR, L_ICESIMPLE, L_IHANEY, L_HADCM4O2I, L_OADGHR2,
     & L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,
     & L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,
     & L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,
     & L_OBIOLOGY, L_OCARB14, L_OCARBON, L_OCNASSM,
     & L_OCYCLIC, L_OFILTER, L_OFREESFC, L_FLUXD,
     & L_OHANEY, L_OHMEAD, L_OICECOUP,
     & L_OIMPADDF, L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC,
     & L_OLISS, L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC, L_OPSEUDIC,
     & L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,        
     & L_OSYMM, L_OVARYT, L_RIVERS, L_SEAICE, L_OCONJ,
     & L_TRANGRID, L_UPWIND, L_OSTVEW, L_OPMSL, L_OTIDAL, L_OPRINT,  
     & L_OFOURW, L_ODELPLUS, L_OTROPIC
     &, L_OISOMOM,L_OISOGMSKEW,L_OISOGM,L_OBIHARMGM,L_OVISHADCM4
     &, L_OMEDOUT
     &, L_OCONVROUS
     &, L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER
     &, L_OVISBECK
     &,L_OBIMOM
     &, L_OQLARGE
     &,L_OCOMP
     &, L_OMEDADV,L_OHUDOUT
     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA 
     &, L_CO2O_INTERACTIVE
     &, L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME
! *IF DEF,OCNASSM
     &,L_SLOPEMAX,L_COXCNVC
!        additions for control of ocean assimilation
     &              ,LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP
     &              ,LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC
     &              ,LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR
! *ENDIF OCNASSM

      COMMON / CNTLCOCN /

     & O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,
     & O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,
     & O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,
     & O_ASSIM_START_HR, O_ASSIM_END_HR, L_FLUXCORR, L_OGLOBAL,
     & L_ICEFREEDR, L_ICESIMPLE, L_IHANEY, L_HADCM4O2I, L_OADGHR2,
     & L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,
     & L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,
     & L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,
     & L_OBIOLOGY, L_OCARB14, L_OCARBON, L_OCNASSM, 
     & L_OCYCLIC, L_OFILTER, L_OFREESFC, L_FLUXD,
     & L_OHANEY, L_OHMEAD, L_OICECOUP,
     & L_OIMPADDF, L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC,
     & L_OLISS, L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC, L_OPSEUDIC,
     & L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,               
     & L_OSYMM, L_OVARYT, L_RIVERS, L_SEAICE, L_OCONJ, 
     & L_TRANGRID, L_UPWIND, L_OSTVEW, L_OPMSL, L_OTIDAL, L_OPRINT,
     & L_OFOURW, L_ODELPLUS, L_OTROPIC, L_OZVRT 
     &, L_OMEDOUT,L_OISOMOM,L_OISOGMSKEW,L_OISOGM
     &, L_OBIHARMGM,L_OVISHADCM4
     &, L_OCONVROUS
     &, L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER
     &, L_OVISBECK
     &,L_OBIMOM

     &, L_OQLARGE
     &,L_OCOMP
     &, L_OMEDADV,L_OHUDOUT

     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA 
     &, L_CO2O_INTERACTIVE
     &, L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME
! *IF DEF,OCNASSM
     &,L_SLOPEMAX,L_COXCNVC
!        additions for control of ocean assimilation
     &              ,LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP
     &              ,LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC
     &              ,LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR
! *ENDIF OCNASSM

! ----------------------- Comdeck: CNTLSLB  ----------------------------
! Description: COMDECK defining Control variables for the Slab
!              internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      LOGICAL                       ! Logical switches for:
     &   L_THERM        ,  !         :Coupled model ice thermodynamics
     &   L_IDYN         ,  !Slab     :Cavitating fluid ice dynamics
     &   L_IDRIF        ,  !model    :Simple ice depth advection
     &   L_SLBADV          !

      NAMELIST / NLSTCSLB /
     & L_THERM, L_IDYN, L_IDRIF, L_SLBADV

      COMMON / CNTLCSLB /
     & L_THERM, L_IDYN, L_IDRIF, L_SLBADV
! ----------------------- Comdeck: CNTLWAV  ----------------------------
! Description: COMDECK defining Control variables for the Wave
!              internal model, and its runtime constants (if any).
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  4.1  23/02/96  New comdeck for wave sub-model.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Parameter declarations
!??   INTEGER MAXSECTS            ! Max. no. of code sections
!??   PARAMETER (MAXSECTS=99)
!
!   Type declarations
!
      INTEGER
     & W_N_SRCE,    ! no.of source timesteps per propagation timestep
     & W_ISHALLOW,  ! 1 for deep, otherwise shallow water
     & W_IREFRACT,  ! refraction options, 0 = none,
!                   ! 1 = depth, 2 = depth & current
     & W_ICASE,     ! 1 for spherical propagation, otherwise Cartesian
     & W_IPER       ! 1 for , otherwise

      LOGICAL
     & L_WAVASSM    ! True if assimilation requested
!
      CHARACTER*5 W_ASSIM_MODE     ! Switch for BG/AC/3DVAR/4DVAR assm
!
!     CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section versions
!
      NAMELIST / NLSTCWAV /
     & W_N_SRCE, W_ISHALLOW, W_IREFRACT, W_ICASE, W_IPER, L_WAVASSM,
     & W_ASSIM_MODE
      COMMON / CNTLCWAV /
     & W_N_SRCE, W_ISHALLOW, W_IREFRACT, W_ICASE, W_IPER, L_WAVASSM,
     & W_ASSIM_MODE

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
!
! Description: This comdeck defines the Run Time Constansts
!              used in the internal models.
!              Currently set up for Atmosphere model. Other
!              internal models to be catered for in future releases.
!
! Current Code Owner: D. Robinson
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.5   28/03/95  Original code. D. Robinson
!   4.0   17/08/95  Add pressure_altitude.  RTHBarnes.
!   4.1   12/03/96  Add SO2_HIGH_LEVEL.  D. Robinson.
!   4.3    18/3/97  Replace three CO2 variables with one - if it is not 
!                    constant CSENARIO will hold the info.  Will Ingram
!   4.3     07/03/97 Add KAY_LEE_GWAVE    S.Webster
!   4.4   26/09/97  Add MPARWTR ,ANVIL_FACTOR & TOWER_FACTOR parameters
!                   for convection.                         J.M.Gregory
!   4.4   15/10/97  Add PHENOL_PERIOD and TRIFFID_PERIOD.   R.A.Betts
!   4.5   07/07/98  Add SOOT_HIGH_LEVEL             M Woodage
!   4.5   19/11/98  Move PHENOL_PERIOD and TRIFFID_PERIOD to NLSTCATM.  
!                   R.A.Betts

! Type Declarations :

C CMAXSIZE must be called before *CALLing this comdeck.

      INTEGER
     &  START_LEVEL_GWDRAG       ! Starting level for GW drag test
     & ,BOTTOM_VDIF_LEVEL        ! Bottom level for vertical diffusion
     & ,TOP_VDIF_LEVEL           ! Top    level for vertical diffusion
     & ,TRAC_ADV_LEVELS          ! No of levels which tracers advected
     & ,DIFF_EXP(MAX_P_LEVELS)   ! Exponent for theta,u,v diffusion
     & ,DIFF_EXP_Q(MAX_P_LEVELS) ! Exponent for q diffusion
     & ,SO2_HIGH_LEVEL ! Model level for SO2 (High level) emissions
     & ,SOOT_HIGH_LEVEL  ! Model level for high level soot emiss

      REAL
     &  LATITUDE_BAND        ! Latitude at which vert diffusion starts
     & ,VERTICAL_DIFFUSION   ! Vertical diffusion coefficient
     & ,CO2_MMR              ! CO2 concentration (if constant)
     & ,GRAV_WAVE_SPEED      ! Gravity wave speed
     & ,FILTERING_SAFETY_FACTOR ! Safety factor used for fourier
!                               ! filtering, see UMDP 10
     & ,WIND_LIMIT           ! Wind speed threshold for top level
     & ,DIV_LIMIT            ! Divergence threshold for top level
     & ,KAY_GWAVE            ! Surface stress constant for GW drag
     & ,KAY_LEE_GWAVE        ! Surface stress constant for GW 
     &                       ! trapped lee wave drag   
     & ,NU_BASIC             ! Advection term
     & ,CW_LAND              ! ) Threshold cloud liquid water content
     & ,CW_SEA               ! ) over land/sea for conversion to ppn
!                            !   (kg water per m**3)
     & ,PRESSURE_ALTITUDE    ! for steep slope test in diffusion code
     & ,DIFF_COEFF(MAX_P_LEVELS)   ! Diffusion coeffs for Theta, u & v
     & ,DIFF_COEFF_Q(MAX_P_LEVELS) ! Diffusion coeffs for q.
     & ,DIV_DAMP_COEFF_ASSM(MAX_P_LEVELS)
!                            ! Divergence damping coeffs-assimilation
     & ,DIV_DAMP_COEFF_FC(MAX_P_LEVELS)
!                            ! Divergence damping coeffs-forecast
     & ,RHCRIT(MAX_P_LEVELS) ! Critical RH for layer cloud formation
     & ,MPARWTR              ! Reservoir of conv. cloud water left in
!                            ! layer whilst rest is rained out (kg/kg)
     & ,ANVIL_FACTOR         ! Determines vertical distribution of
!                            ! convective cloud amount
     & ,TOWER_FACTOR         ! Determines vertical distribution of
!                            ! convective cloud amount
     & ,UD_FACTOR            ! Updraught factor used in calculation
!                            ! of conv. water path if L_CCW = true.

! COMMON block for Atmosphere Run Time Constants :

      COMMON /CRUNTIMA/
     &  LATITUDE_BAND,  VERTICAL_DIFFUSION,
     &  CO2_MMR,  GRAV_WAVE_SPEED,  KAY_GWAVE, KAY_LEE_GWAVE,
     &  NU_BASIC,  FILTERING_SAFETY_FACTOR,  START_LEVEL_GWDRAG,
     &  BOTTOM_VDIF_LEVEL,  TOP_VDIF_LEVEL, PRESSURE_ALTITUDE,
     &  DIFF_COEFF,  DIFF_COEFF_Q,  DIFF_EXP,  DIFF_EXP_Q,
     &  DIV_DAMP_COEFF_ASSM,  DIV_DAMP_COEFF_FC,  RHCRIT,
     &  TRAC_ADV_LEVELS,  CW_SEA,  CW_LAND,  WIND_LIMIT,  DIV_LIMIT
     & ,SOOT_HIGH_LEVEL
     & ,SO2_HIGH_LEVEL
     & ,MPARWTR, ANVIL_FACTOR, TOWER_FACTOR, UD_FACTOR                  

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

      INTEGER
     & A_STEP             ! IN Current timestep in atmosphere model

      INTEGER
     & FIRST_POINT         ! LOCAL First P-point to be processed.
     &,LAST_POINT          ! LOCAL Last P-point to be processed.
     &,LAND1               ! LOCAL First land point to be processed.
     &,LAND_PTS            ! LOCAL Number of land point to be processed.
     &,TILE_PTS(NTYPE)              ! LOCAL Number of land points which
C                                   !       include the nth surface type
     &,TILE_INDEX(LAND_FIELD,NTYPE) ! LOCAL Indices of land points which
C                                   !       include the nth surface type
     &,NSTEP_TRIF                   ! LOCAL Number of atmospheric
C                                   !       timesteps between calls to
C                                   !       TRIFFID.
     &,I                            ! LOCAL Loop counter for all points
     &,L                            ! LOCAL Loop counter for land points
     &,N          ! ** TEMPORARY ** loop counter for types

      REAL
     & Z0_LAND(LAND_FIELD)          ! LOCAL Z0 on land points

      INTEGER ICODE                 ! LOCAL Internal return code
      CHARACTER*80 CMESSAGE         ! LOCAL Internal error message

      EXTERNAL
     & INIT_MIN
     &,TILEPTS
     &,SPARM
     &,INIT_ACC

!-----------------------------------------------------------------------
! Set the TYPFLDPT variables, grid pointers and indices.
!-----------------------------------------------------------------------
! Comdeck SETFLDPT
! Set up all the TYPFLDPT variables

! Pointers to useful positions in the horizontal field

! FIRST_ROW is the first updatable row of data - missing out the
!   North Pole row (global) or Northern-most row (LAM)
! TOP_ROW_START is the first point containing valid data. For a
!   processor at the top of the LPG, this means the North Pole row
!   (global) or Northern-most row (LAM), for all other processors
!   it means the first point of the Northern halo
! START_POINT_NO/INC_HALO is the first non-polar point which does not
! or does lie on a halo.

      IF (attop) THEN
        FIRST_ROW=Offy+2
        TOP_ROW_START=(FIRST_ROW-2)*ROW_LENGTH + 1
        START_POINT_NO_HALO=(FIRST_ROW-1)*ROW_LENGTH+1
        START_POINT_INC_HALO=START_POINT_NO_HALO
      ELSE
        FIRST_ROW=Offy+1
        TOP_ROW_START=1
        START_POINT_NO_HALO=Offy*ROW_LENGTH+1
        START_POINT_INC_HALO=1
      ENDIF

! P/U_LAST_ROW is the last updatable row of data - missing out the
!   South Pole row (global)  or Southern-most row (LAM)
! P/U_BOT_ROW_START is the first point of the last row containing
!   valid data. For a processor at the base of the LPG, this means
!   the South Pole row (global) or Southern-most row (LAM), for all
!   other processors it means the first point of the Southern halo.
! LAST_P/U_FLD_PT is the last point in the field containing local
!   rows - ie. it does not include the Southern halo.
! END_P/U_POINT_NO/INC_HALO is the last point of data before the
! polar row which either does not or does include the halo

      IF (atbase) THEN
        P_LAST_ROW=P_ROWS-Offy-1
        U_LAST_ROW=P_ROWS-Offy-2
        P_BOT_ROW_START=P_LAST_ROW*ROW_LENGTH + 1
        U_BOT_ROW_START=U_LAST_ROW*ROW_LENGTH + 1
        LAST_P_FLD_PT=(P_ROWS-Offy)*ROW_LENGTH
        LAST_U_FLD_PT=(P_ROWS-Offy-1)*ROW_LENGTH
        END_P_POINT_NO_HALO=P_LAST_ROW*ROW_LENGTH
        END_U_POINT_NO_HALO=U_LAST_ROW*ROW_LENGTH
        END_P_POINT_INC_HALO=END_P_POINT_NO_HALO
        END_U_POINT_INC_HALO=END_U_POINT_NO_HALO
      ELSE
        P_LAST_ROW=P_ROWS-Offy
        U_LAST_ROW=P_LAST_ROW
        P_BOT_ROW_START=(P_ROWS-1)*ROW_LENGTH + 1
        U_BOT_ROW_START=P_BOT_ROW_START
        LAST_P_FLD_PT=(P_ROWS-Offy)*ROW_LENGTH
        LAST_U_FLD_PT=LAST_P_FLD_PT
        END_P_POINT_NO_HALO=P_LAST_ROW*ROW_LENGTH
        END_U_POINT_NO_HALO=U_LAST_ROW*ROW_LENGTH
        END_P_POINT_INC_HALO=P_FIELD
        END_U_POINT_INC_HALO=U_FIELD
      ENDIF

! FIRST_FLD_PT is the first point in the field containing local
!   rows - ie. it does not include the Northern halo.
      FIRST_FLD_PT=Offy*ROW_LENGTH+1

      FIRST_ROW_PT=1+Offx  ! first local data point along a row
!                          ! (ie. first non-halo point)
      LAST_ROW_PT=ROW_LENGTH-Offx  ! last data point on row
!                                  ! (ie. last non-halo point)

! and the size of a global row and global fields

      GLOBAL_ROW_LENGTH=glsize(1)
      GLOBAL_P_FIELD=glsize(1)*glsize(2)
      GLOBAL_U_FIELD=glsize(1)*(glsize(2)-1)

      upd_P_ROWS=P_LAST_ROW-FIRST_ROW+1
      upd_U_ROWS=U_LAST_ROW-FIRST_ROW+1

      FIRST_VALID_PT=TOP_ROW_START ! first point on field that must
!                                  ! contain valid (sensible) data
      LAST_P_VALID_PT=P_BOT_ROW_START-1+ROW_LENGTH
      LAST_U_VALID_PT=U_BOT_ROW_START-1+ROW_LENGTH
!                                  ! last point on field that must
!                                  ! contain valid (sensible) data

      tot_P_ROWS=P_FIELD/ROW_LENGTH
      tot_U_ROWS=U_FIELD/ROW_LENGTH

      VALID_P_ROWS=(LAST_P_VALID_PT-FIRST_VALID_PT+1)/ROW_LENGTH
      VALID_U_ROWS=(LAST_U_VALID_PT-FIRST_VALID_PT+1)/ROW_LENGTH

      MY_PROC_ID=mype
      NP_PROC_ID=first_comp_pe
      SP_PROC_ID=last_comp_pe

      GC_ALL_GROUP=0  ! hardwired number, copied from GCOM include file
      GC_ROW_GROUP=gc_proc_row_group
      GC_COL_GROUP=gc_proc_col_group

      N_PROCS=nproc

      EW_Halo = Offx
      NS_Halo = Offy

      halo_4th=2
      extra_EW_Halo = halo_4th-EW_Halo
      extra_NS_Halo = halo_4th-NS_Halo

      LOCAL_ROW_LENGTH=ROW_LENGTH
      FIRST_GLOBAL_ROW_NUMBER = datastart(2)

      at_top_of_LPG=attop
      at_right_of_LPG=atright
      at_base_of_LPG=atbase
      at_left_of_LPG=atleft
! End of comdeck SETFLDPT
      FIRST_POINT=START_POINT_NO_HALO
      LAST_POINT=END_P_POINT_INC_HALO
      LAND1 = 1
      LAND_PTS = 0
      DO L=1,LAND_FIELD
        IF ( LAND_LIST(L) .LT. FIRST_POINT ) THEN
          LAND1 = LAND1 + 1
        ELSEIF ( LAND_LIST(L) .LE. LAST_POINT ) THEN
          LAND_PTS = LAND_PTS + 1
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! If TRIFFID on, call INIT_MIN to ensure PFT fractions are GE minimum
! fraction except where vegetation excluded by ice, water or urban
!-----------------------------------------------------------------------
      IF (L_TRIFFID) THEN
        CALL INIT_MIN(LAND_FIELD,LAND1,LAND_PTS,D1(JFRAC_TYP),
     &                D1(JSOIL_CARB))
      ENDIF

C-----------------------------------------------------------------------
C Call TILEPTS to initialise TILE_PTS and TILE_INDEX
C-----------------------------------------------------------------------
      CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &             D1(JFRAC_TYP),TILE_PTS,TILE_INDEX)

C-----------------------------------------------------------------------
C Initialise tiled and gridbox mean vegetation parameters
C-----------------------------------------------------------------------
      WRITE(6,*) 'INITVEG: CALLING SPARM'
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,TILE_PTS,TILE_INDEX,
     &            D1(JSOIL_ALB),D1(JFRAC_TYP),D1(JCANHT_PFT),
     &            D1(JLAI_PFT),D1(JMDSA),D1(JSFA),D1(JCATCH_NIT),
     &            Z0_LAND,D1(JZ0_TYP))

C-----------------------------------------------------------------------
C Copy Z0 from land field to full field
C-----------------------------------------------------------------------
      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_LIST(L)
        D1(JZ0+I-1)=Z0_LAND(L)
      ENDDO

      IF (L_TRIFFID) THEN
C-----------------------------------------------------------------------
C If this is an NRUN and re-start from mid-way through a TRIFFID calling
C period has not been requested: (i) initialise accumulation prognostics
C to zero, (ii) set TRIFFID_PERIOD in integer header, and
C (iii) initialise ASTEPS_SINCE_TRIFFID integer header to zero.
C If mid-period restart is requested then leave the accumulated fields
C unchanged, and if a new calling period is specified then reset
C calling period header the new value provided that the number of
C atmosphere timesteps since the last call to TRIFFID does not exceed
C the new calling period .
C A_INTHD(22) holds TRIFFID_PERIOD in days.
C A_INTHD(23) holds the number of atmosphere timesteps since the last
C call to TRIFFID.
C-----------------------------------------------------------------------
        IF (A_STEP.EQ.0) THEN
          IF (L_NRUN_MID_TRIF) THEN

            IF (TRIFFID_PERIOD.NE.A_INTHD(22)) THEN
              NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/
     &        SECS_PER_STEPim(atmos_im))

              IF (A_INTHD(23).GT.NSTEP_TRIF) THEN
                WRITE(6,*) '**ERROR IN TRIFFID** YOU HAVE SELECTED TO'
                WRITE(6,*) 'START MID-WAY THROUGH A TRIFFID CALLING'
                WRITE(6,*) 'PERIOD BUT YOUR INITIAL DUMP CONTAINS'
                WRITE(6,*) 'PROGNOSTICS ACCUMULATED OVER A PERIOD'
                WRITE(6,*) 'LONGER THAN THE NEW CALLING PERIOD'

              ELSE

                A_INTHD(22)=TRIFFID_PERIOD

              ENDIF
            ENDIF

          ELSE

            CALL INIT_ACC(LAND_FIELD,
     &                    D1(JNPP_PFT_ACC),
     &                    D1(JG_PHLF_PFT_ACC),D1(JRSP_W_PFT_ACC),
     &                    D1(JRSP_S_ACC),ICODE,CMESSAGE)

            A_INTHD(22)=TRIFFID_PERIOD
            A_INTHD(23)=0

          ENDIF
        ENDIF

      ELSE IF (L_PHENOL) THEN

C-----------------------------------------------------------------------
C Initialise accumulated leaf turnover rate to zero
C-----------------------------------------------------------------------
      DO L = 1,LAND_FIELD
        D1(JG_LF_PFT_ACC+L-1) = 0.0
      ENDDO


      ENDIF

      RETURN
      END
