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
CLL  Routine: EXITPROC -------------------------------------------------
CLL
CLL  Purpose: Tidies up at the end of the run, and takes certain
CLL           actions in the case of model failure (as indicated by
CLL           ICODE input)
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   T.C.Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: C0
CLL
CLL  Project task: C0
CLL
CLL  External documentation: On-line UM document C0 - The top-level
CLL                          control system
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE EXITPROC(ICODE,CMESSAGE)
      IMPLICIT NONE
      INTEGER ICODE            ! INOUT - Error code from model
      CHARACTER*(80) CMESSAGE   ! IN    - Error message from model
C
C*----------------------------------------------------------------------
C  Common blocks
C
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
C*L --------------------- Comdeck: CHISTORY ----------------------------
CLL
CLL  Purpose: COMMON block for history data needed by top level (C0)
CLL           routines, and passed from run to run.  Mostly set by
CLL           the User Interface.
CLL
CLL           Note that CHISTORY *CALLs ALL individual history comdecks
CLL
CLL  Author : A. Sangster
CLL
CLL  Model            Modification history
CLL version  Date
CLL  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
CLL                 contents.  RTHBarnes.
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLLEND----------------------------------------------------------------
C*
CCC   *CALL CHSUNITS
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
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  18/04/96  Add RUN_IN for qxhistreport.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*10 RUN_HIST_TYPE       ! Type of history file
      CHARACTER*8  RUN_TYPE            ! Type of run
      CHARACTER*14 RUN_COMPCODE        ! Run completion code
      CHARACTER*14 RUN_LAST_MEAN       ! Last mean dump created by run
C APPREARS UNUSED                      ! for pp fields
      CHARACTER*1  RUN_MEANS_TO_DO     ! Flag indicating the run stopped
C                                      ! before creating next mean dump
      CHARACTER*1  RUN_OCEAN_FIRST     ! Flag set to true if ocean to be
C                                      ! run first
      CHARACTER*8  RUN_JOB_NAME        ! Jobname this run
      CHARACTER*5  RUN_ID              ! Expt./Job id for this run
      CHARACTER*1  RUN_RESUBMIT        ! Flag controlling auto resubmit
      CHARACTER*12 RUN_RESUBMIT_Q      ! Job queue to which resubmit run
      CHARACTER*20 RUN_RESUBMIT_TIME   ! Time at which run resubmits
      CHARACTER*6  RUN_RESUBMIT_CPU    ! Time limit for resubmitted job
      CHARACTER*6  RUN_RESUBMIT_MEMORY ! Resubmitted job's memory limit
      CHARACTER*2  RUN_RESUBMIT_PRTY   ! Resubmitted job intra q prty
      CHARACTER*8  RUN_RESUBMIT_JOBNAME! Resubmitted jobname
      CHARACTER*1  FT_ACTIVE(20:NUNITS) ! "Y" if file partly written
C
C
C History Common Block for overall model character variables.
C
      COMMON /CHISTO/
     & RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,
     & RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID, 
     & RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,
     & RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE
C
      NAMELIST /NLCHISTO/
     & RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,
     & RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID, 
     & RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,
     & RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
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
     & LENGTH(N_INTERNAL_MODEL_MAX) ! No. of tsteps completed this run
     &,ACTUAL_ENDT(6,N_INTERNAL_MODEL_MAX) ! Model end time this run
!     These 2 appears to be purely diagnostic, and not really used.

     &,H_STEPim(N_INTERNAL_MODEL_MAX) ! History block copy of A/O_STEP
!                                   ! held in COMDECK CTIME
     &,H_GROUPim(N_INTERNAL_MODEL_MAX) ! No of steps in coupling period
     &,MEAN_OFFSETim(N_INTERNAL_MODEL_MAX) ! No of means activated
     &,OFFSET_DUMPSim(N_INTERNAL_MODEL_MAX) ! Offset between
!                 MEAN_REFTIME and model basis time  (in model dumps)
     &,MEAN_NUMBERim(N_INTERNAL_MODEL_MAX)  ! No of mean periods chosen
     &,RUN_MEANCTL_INDICim(4,N_INTERNAL_MODEL_MAX) ! Indicators used to
!    correct logical units are used for atmos/ocean partial sum dump I/O
!
C
C History Common Block for generic model integer variables.
C
      COMMON /IHISTG/
     & H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim
C
      NAMELIST /NLIHISTG/
     & H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim
C
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.4  30/05/97  Added vars LASTATMim, CURRATMim, LASTDMPim.  K Rogers
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*14 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
      CHARACTER*80 RESTARTim(N_INTERNAL_MODEL_MAX) !current restart dump
      CHARACTER*14 SAFEDMPim(N_INTERNAL_MODEL_MAX) 
! Name of old safe restart dump
      CHARACTER*14 NEWSAFEim(N_INTERNAL_MODEL_MAX)  
! Name of new safe restart dump
      CHARACTER*14 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos restart dump
!                                                  ! until ocean dump
      CHARACTER*14 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of 
!                                                  ! current atmos 
!                                                  ! restart dump
      CHARACTER*14 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos/ocean dumps
!                                                  ! until meaning done

!
!
! History Common Block for generic model characters variables.
!
      COMMON /CHISTG/
     & END_DUMPim, RESTARTim
     & , SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim
C
      NAMELIST /NLCHISTG/
     & END_DUMPim, RESTARTim
     & , SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim
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
C
C
C  Subroutines called
C
CL----------------------------------------------------------------------
CL 1. If fatal error occurred in main body of model, suppress resubmit
CL    switch to prevent model resubmit (if activated)
CL
      IF (ICODE.GT.0) THEN
        RUN_RESUBMIT="N"
      ENDIF
CL
CL 1.1  Reset error code
CL
      ICODE=0
 999  CONTINUE
CL----------------------------------------------------------------------
CL 2. Close named pipe unit used for communication with server
CL
      CLOSE(8)
C
      RETURN
CL----------------------------------------------------------------------
      END
