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
!+ Reads History & control files; also interim for CRUN and HK for Op.
!
! Subroutine Interface:
      SUBROUTINE UM_SETUP(ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Reads History and control namelist files for all runs.
!   For CRUNs also reads interim control file of values changed by user.
!   For operational runs reads housekeeping file.
!
! Method: as above
!
! Current Code Owner: RTHBarnes.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.5  04/05/95  Sub-models stage 1: New routine.  RTHBarnes.
!  4.0  06/09/95  Add Timer stop section.  RTHBarnes.
!  4.0  06/12/95  Check env.var. TYPE for CRUN as well as HSTEPim. RTHB
!  4.1  14/05/96  Change GETENV to more portable FORT_GET_ENV  P.Burton
!  4.2  11/09/96  Remove EXTERNAL GET_ENV statement   P.Burton
!  4.5  23/10/98  Introduce Single Column Model. J-C Thil.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: C0
! System Task:              C0
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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


! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
      Integer ICODE ! Error code
      Character*80 CMESSAGE ! Error message
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
      INTEGER I              ! Loop variable
      INTEGER UNITICTL       ! Unit no. for Interim Control file
      LOGICAL CRUN           ! T if continuation run
      INTEGER J              ! dummy for use with FORT_GET_ENV
      CHARACTER*4  TYPE      ! Type of run - NRUN or CRUN

! Local dynamic arrays:

! Function & Subroutine calls:
      External TIMER, READHIST, READCNTL, READHK, EREPORT, ABORT
     & ,FORT_GET_ENV

!- End of header
!
! ----------------------------------------------------------------------
!  0. Start Timer running
!
      CALL TIMER('UM_SETUP',3)
      ICODE=0
! ----------------------------------------------------------------------
!  1. Read History file.
!
      CALL READHIST ( IHIST_UNIT,ICODE,CMESSAGE )
      IF(ICODE .GT. 0) GOTO 999

! ----------------------------------------------------------------------
!  2. Read Control file on standard input.
!
      CALL READCNTL ( 5,ICODE,CMESSAGE )
      IF(ICODE .GT. 0) GOTO 999

! ----------------------------------------------------------------------
!  3. Read Interim Control file for CRUNs only.
!
      CALL FORT_GET_ENV('TYPE',4,TYPE,4,J)
      CRUN = .false.
      DO  I = 1,N_INTERNAL_MODEL_MAX
        IF (H_STEPim(I) .gt. 0) THEN
          CRUN = .true.
          GO TO 300
        END IF
      END DO
  300 CONTINUE
      IF (CRUN .or. TYPE.eq.'CRUN') THEN
        WRITE(6,*)' UMSETUP; CRUN, read CONTCNTL from unit 9'
!  Open and read Interim Control file
        UNITICTL = 9
        CALL READCNTL ( UNITICTL,ICODE,CMESSAGE )
        IF(ICODE .GT. 0) GOTO 999
        CLOSE (UNITICTL)
      END IF

! ----------------------------------------------------------------------
!  4. Read Housekeeping file for operational runs only.
!
      IF(MODEL_STATUS .EQ. 'Operational') THEN
        CALL READHK(HKFILE_UNIT,ICODE,CMESSAGE)
        IF(ICODE .GT. 0) GOTO 999
      END IF

! ----------------------------------------------------------------------
!  5. Output any error/warning message and stop if error.
!
  999 CONTINUE
      IF(ICODE .NE. 0) CALL EREPORT(ICODE,CMESSAGE)
      IF(ICODE .GT. 0) CALL ABORT

! ----------------------------------------------------------------------
!  6. Stop Timer for this routine.
!
      CALL TIMER('UM_SETUP',4)

      RETURN
      END
