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
CLL  Subroutine STWORK -------------------------------------------------
CLL
CLL  Purpose: Processes all the STASHlist entries for a particular
CLL item and section number after a timestep.  Raw input diagnostics to
CLL STWORK come either from D1, the main data array, or STASH_WORK, a
CLL work array dimensioned by the control level and passed in.  Each
CLL field is spatially processed, temporally processed and then output.
CLL The output destination can either be to an address in D1 or to a PP
CLL fieldsfile on a given unit.  In either case a PP-type LOOKUP header
CLL is generated to describe the contents of the output field.  Now
CLL handles atmosphere or ocean diagnostics according to arguments
CLL passed by calling routine.
CLL
CLL  Author:   T.Johns/S.Tett
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.0   31/12/92  Reference to deck name removed from message
CLL                   when deckname STWORK2  changed to STWORK1A  MJH
CLL   3.1   19/01/93  Enforce processing thro' SPATIAL loop if the input
CLL                   and output pseudo-levels lists differ (bugfix).
CLL   3.1   24/02/93  Correct outstanding problems with timeseries (ST).
CLL         16/03/93  Correct coding errors for ppx_ocn GR types (SI).
CLL   3.1   15/01/93  Add NUNITS to argument list to cope with increases
CLL                   in i/o units for 'C' portable code. R.Rawlins
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.2   06/04/93  Correct base_level in no-SPATIAL-processing case
CLL                   when OCEAN diagnostic on input levels list.  Trap
CLL                   error return codes from STOCGT. Bugfix for packing
CLL                   of timeseries (extraw_hdr not extraw). Initialise
CLL                   extraw_hdr (MOS case). Remove rogue comma.(TCJ/ST)
CLL   3.2   02/07/93  Add arguments needed in STOCGT owing to dynamic
CLL                   allocation control level changes. R.T.H.Barnes.
CLL   3.3   30/03/94  Pass lenout dimension to MULTI_SPAIAL and remove
CLL                   redundant addr_out addressing on ppfield; correct
CLL                   bug in SPATIAL processing of ocean fields on
CLL                   input levels lists.  Author Tim Johns.
CLL   3.3   17/11/93 Declaring MODEL_FT_UNIT(NUNITS) after NUNITS
CLL                  has been defined. (N.Farnon)
CLL   3.3   17/09/93  Set LOGICAL lmasswt to denote that level-by-level
CLL                   mass-weights exist and pass to (MULTI)SPATIAL.
CLL                   Also correct error in passing addressed field to
CLL                   SPATIAL in vertical/global mean case (TCJ).
CLL  3.4   04/08/94  No packing indicator change from -26 to -99  PJS
CLL   3.4   07/10/94 : Move call to grib_file so that pp headers are
CLL                   set correctly R A Stratton.
CLL   3.4   15/06/94  Force soft failure if number of fields output to
CLL                   a PP file exceeds the number of pre-allocated
CLL                   headers. R.Rawlins
CLL   3.4   22/07/94  Corrects uninitialised argument to PP_HEAD when
CLL                   stashing data to D1
CLL   3.4   13/01/94  Correct bug in spatial processing loop - cater for
CLL                   C-U and C-V grid types.           T.Johns
CLL   4.0   26/07/95  Allow integer or logical fields to be passed thro'
CLL                   STASH provided no spatial processing done. T.Johns
CLL   3.5   07/04/95  Changes for stage 1 of submodels project. K Rogers
CLL                   In STOCGT remove ARGSIZE and pass in:
CLL                   JMT=T_ROWS,IMT=ROW_LENGTH,KM=P_LEVELS or
CLL                   JMT=T_ROWS,IMTM2=ROW_LENGTH,KM=P_LEVELS
CLL   3.5  24/03/95    Changed OPEN to FILE_OPEN and
CLL                    CLOSE to FILE_CLOSE    P.Burton
CLL   4.0   12/09/95  Remove CMESSAGE='' after call to PP_HEAD.
CLL                   Call PP_HEAD with AK, BK, AHH, BKH, T_levels and
CLL                   JJ (level list value) and without RUN_INDIC_OP
CLL                   (Andrew Brady)  
!LL   4.0  06/09/95  Pass in superarray from higher levels. K Rogers
!     4.0   10/03/95  Allow alternative packing method for grib output.
!                     R A Stratton
CLL   4.0   15/06/95  Reduce i/o required by STASH. Original code read
CLL                   the entire set of pp lookup tables (at least
CLL                   PP_LEN2_LOOK, usually 4096) records, then re-wrote
CLL                   the same set for each call to STASH. This is
CLL                   replaced here by only reading the last pp header
CLL                   written to and subsequently only writing out the
CLL                   next set of pp headers serviced by the STASH call.
CLL                   Rick Rawlins.
!     4.1    Apr. 96  Rationalise *CALLs     S.J.Swarbrick
CLL   4.1   31/05/96  Add calls to STWVGT for wave model. K Rogers
!     4.1   03/04/96  New argument DUMP_PACK : Use to control
!                     packing in dumps. D. Robinson
!     4.2   10/09/96  MPP modifications to STASH           P.Burton
!     4.3   13/3/97   Further MPP STASH modifications      P.Burton
!LL  4.3   30/04/97  Added code to use UM_SECTOR_SIZE to make transfers
!LL                  well-formed.
!LL                  B. Carruthers  Cray Research.
!     4.3   07/05/97   Correct error in STASH-time processed
!                      diagnostics  affecting Ocean and Wave
!                      models only                   M.Carter
!     4.3   17/4/97   Parallelisation of COEX              D.Salmond
!     4.4   25/11/96  New option - mean timeseries. R A Stratton.       
!     4.4   25/06/97   Set PPHORIZ_OUT and NROWS/COLS_OUT correctly
!                      when outputing timeseries.          P.Burton
!LL   4.4  28/07/97  Changes connected with allowing reinitialisation of
!LL                  PP files on Gregorian month boundaries. M.Gallani
!LL   4.4   16/06/97  Add processing after the write, so
!LL                   that all the processors know the answer
!LL                     Author: Bob Carruthers, Cray Rsearch.
!     4.4   12/06/97  MPP improvements                 P.Burton
!     4.4    7/10/97  Force explicit fail of model when attempting to 
!                     write to a file on a re-initialised stream before
!                     it has been opened. R.Rawlins
!LL   4.5   13/01/98  Added global_LENOUT argument, and replace
!LL                   SHMEM COMMON blocks with dynamic arrays
!LL                   for buf and buf3 arrays.         P.Burton
!LL   4.5    18/09/98  Corrected non-standard FORMAT statments
!LL                                                  P.Burton
!LL  4.5   28/05/98  Code for parallel processing in COEX Packing
!LL                    Author: Paul Burton & Bob Carruthers
CLL  4.5  28/02/96  Flush buffer for non-reinitializable files to try
CLL                 to avoid problems with continuation runs following
CLL                 hard failures.  RTHBarnes.
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered : C3, C4, C8
CLL
CLL  Project task: C4
CLL
CLL  External documentation : UMDP no C4
CLL
C*L  Interface and arguments: ------------------------------------------
      SUBROUTINE STWORK (
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
     &  D1,LEN_TOT,STASH_WORK,STASH_WORK_LEN,LENOUT,
     &  global_LENOUT,
     &  len_ocwork,
     &  IS,IM,ILSTART,ILEND,STEP,steps_per_period,secs_per_period,
     &  previous_time,
     &  STLIST,LEN_STLIST,TOTITEMS,SI,NSECTS,NITEMS,
     &  STASH_LEVELS,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,
     &  STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS,
     &  MAX_STASH_LEVS,STTABL,NSTTIMS,
     &  NSTTABL,STASH_SERIES,STASH_SERIES_LEN,
     &  stash_series_rec_len,stash_series_index,stash_ser_index_size,
     &  MOS_MASK,MOS_MASK_LEN,MOS_OUTPUT_LENGTH,
     &  PP_PACK_CODE,MODEL_FT_UNIT,FT_STEPS,FT_FIRSTSTEP,
     &  FIXHD,INTHD,REALHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,
     &  LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,
     &  LOOKUP,RLOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,
     &  PP_LEN2_LOOKUP,NUNITS,PP_LEN2_LOOK,
     &  LCYCLIC,
     &  T_rows,U_ROWS,ROW_LENGTH,T_field,U_FIELD,T_levels,
     &  FCST_PRD,RUN_INDIC_OP,ELF,FT_LASTFIELD,
     &  sm_ident,im_ident,dump_pack,
     &  stsuparrlen, stsuparr, istsuparr, sa_idx, sa_idxlen,
     &  ICODE,CMESSAGE)
C
      IMPLICIT NONE
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
      INTEGER
     *  TOTITEMS           !IN    MAX NO OF ITEMS IN STASHLIST
     *, NSECTS             !IN    MAX NO OF SECTIONS
     *, NITEMS             !IN    MAX NO OF ITEMS IN A SECTION
     *, LEN_TOT            !IN    LENGTH OF REAL DATA ARRAY D1
     *, len_ocwork         !IN    Length of ocean work array OCWORK
     *, LEN_FIXHD          !IN    LENGTH OF FIXED CONSTANTS
     *, LEN_INTHD          !IN    LENGTH OF INTEGER CONSTANTS
     *, LEN_REALHD         !IN   LENGTH OF REAL CONSTANTS
     *, LEN1_LEVDEPC       !IN First dimension of LEVDEPC
     *, LEN2_LEVDEPC       !IN Second dimension of LEVDEPC
     *, NUM_STASH_LEVELS   !IN   DIMENSION OF STASH_LEVELS
     *, NUM_LEVEL_LISTS    !IN   DIMENSION OF STASH_LEVELS
     *, NUM_STASH_PSEUDO   !IN   Maximum no of pseudo-levels in a list
     *, NUM_PSEUDO_LISTS   !IN   Number of pseudo-level lists
     *, MAX_STASH_LEVS     !IN   Max no of output levels for any diag
     *, LEN1_LOOKUP        !IN   First dimension of LOOKUP/IPPLOOK
     *, LEN2_LOOKUP        !IN   Second dimension of LOOKUP
     *, PP_LEN2_LOOKUP     !IN   Largest poss. value in PP_LEN2_LOOK
     *, NUNITS             !IN   Max i/o FT unit no
     *, PP_LEN2_LOOK(20:NUNITS)!IN   Individual PP_LEN2_LOOKs per unit
     *, PP_PACK_CODE(20:NUNITS)!IN   Packing code per unit
     *, FT_LASTFIELD(20:NUNITS)!IN   Current write posn in each PP file
     *, FT_STEPS(20:NUNITS)    !IN   File reinitialisation freq per unit
     *, FT_FIRSTSTEP(20:NUNITS)!IN   First step file initialised
     *, NSTTIMS            !IN   Number of times against to test
     *, NSTTABL            !IN   Number of STASH timetables
     *, NUM_WORDS          !IN    Number of 64 Bit words to hold DATA
     *, sm_ident           !IN    Submodel identifier
     *, im_ident           !IN    Internal model identifier
     &, dump_pack          !IN    Packing Indicator for Dump
     &, sa_idxlen          !IN    Superarray index length
     &, sa_idx(sa_idxlen)  !IN    Superarray index
     &, stsuparrlen        !IN    Superarray index length
     &, istsuparr(stsuparrlen)!IN Integer superarray    
      CHARACTER*80
     *  MODEL_FT_UNIT(NUNITS)  !IN   Current table of file associations

      INTEGER
     *  FIXHD(LEN_FIXHD)   !IN    ARRAY OF FIXED CONSTANTS
     *, PP_FIXHD(LEN_FIXHD)!IN    Array of fixed constants from PP file
     *, INTHD(LEN_INTHD)   !IN    ARRAY OF integer CONSTANTS
     *, ILSTART            !IN    START OF LOOP OVER ENTRIES
     *, ILEND              !IN    END OF LOOP OVER ENTRIES
     *, IS                 !IN    SECTION NUMBERS
     *, IM                 !IN    ITEM NUMBER
     *, STEP               !IN    MODEL STEP NUMBER
     *, steps_per_period   !IN    No of steps in defining period
     *, secs_per_period    !IN    No of secs in period (define timestep)
     *, previous_time(7)   !IN    Time at start of current step
     *, LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Integer LOOKUP headers
     *, RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Real version of LOOKUP
     *, ICODE              !OUT   RETURN CODE FROM ROUTINE
C
      INTEGER
     *  LENOUT              !IN     Length of largest workfield needed
     &, global_LENOUT       !IN     Output length of largest field
     *, T_field             !IN     NO OF TEMP/PRESS POINTS
     *, U_FIELD             !IN     NO OF U,V POINTS
     *, ROW_LENGTH          !IN     No of points per row
     *, U_ROWS              !IN     No of U,V, rows
     *, T_rows              !IN     No of PRESS/TEMP rows
     *, T_levels            !IN     No of model Press/Temp levels
     *, STASH_WORK_LEN      !IN     LENGTH of STASH_WORK
     &, MOS_MASK_LEN        !IN     Size of MOS_MASK array
     *, MOS_OUTPUT_LENGTH   !IN     No of MOS data pts extracted
     *, LEN_STLIST          !IN     No of entries in STASHlist
     *, STLIST(LEN_STLIST,TOTITEMS) !IN STASHLIST
     *, SI(NITEMS,0:NSECTS,N_INTERNAL_MODEL) !IN     STASH IN ADDRESS
     *, STTABL(NSTTIMS,NSTTABL)!IN  STASH TIME TABLES
     *, STASH_LEVELS(NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS)
     *, STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
     &, MOS_MASK(MOS_MASK_LEN) ! IN mask used to output data on MOS grid
     *, FCST_PRD            !IN     Forecast period
     *, RUN_INDIC_OP        !IN     Operational Run indicator (ITAB)
C STASH timeseries information
     *, STASH_SERIES_LEN     ! IN no of STASH SERIES records
     *, stash_series_rec_len ! IN length of each record
     *, STASH_SERIES(stash_series_rec_len,STASH_SERIES_LEN)
C                            ! IN individual sample records
     *, stash_ser_index_size ! IN no. of index blocks
     *, stash_series_index(2,stash_ser_index_size)
C                            ! IN index block (1=start, 2=no of records)
     *, EXPPXI               ! Function to extract ppxref info
     *, im_index             ! Internal model index number
     &, N1                   ! Packing Indicator for Lookup(21)
C
      CHARACTER*(80) CMESSAGE ! OUT MESSAGE FROM ROUTINE
C
C
      LOGICAL
     *  LCYCLIC       !IN TRUE if cyclic EW BCs
     *, LAND(T_field) !IN land sea mask
     *, ELF           !IN True if the input grid is rotated Equatorial
C

      REAL
     *  D1(LEN_TOT)                 !IN  REAL DATA ARRAY
     *, REALHD(LEN_REALHD)          !IN  ARRAY OF REAL CONSTANTS
     *, LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC+1)    !IN level dep constants
     *, STASH_WORK(STASH_WORK_LEN)  !IN    INPUT work array to STASH
     &, stsuparr(stsuparrlen)       !IN  Real superarray    

C
C*----------------------------------------------------------------------
C
C Common blocks and PARAMETERs
C
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
CLL  Comdeck: STERR ----------------------------------------------------
CLL
CLL  Purpose: PARAMETER names for STASH processing error codes;
CLL           fatal errors have positive codes, warnings negative.
CLL
CLL  Author:   S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for MPP code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
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
C
C Warning codes
C
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
C
C Error codes
C
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)

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
C
C Subroutines called
C
      EXTERNAL STOCGT,SPATIAL,MULTI_SPATIAL,STLEVELS,TEMPORAL,
     &         PP_FILE,PP_HEAD,BUFFIN,BUFFOUT,EXPPXI,STWVGT
     & ,GRIB_FILE,IOERROR,SETPOS,FILE_OPEN,FILE_CLOSE,FLUSH_BUFFER
C
C Local variables
C
      REAL
     * PPFIELD(LENOUT)     ! Main internal work array
     &, OCWORK(len_ocwork) ! Extra work array needed by ocean
     &                     ! Also used by wave model in STWVGT.
C                          ! - reduced to minimum length in atmos case
     *, AK_LEV(max_stash_levs) ! The values of A at output levels
     *, BK_LEV(max_stash_levs) ! The values of B at output levels
     *, LEVEL(max_stash_levs)  ! The levels of the data as REAL nos
     *, sample_prd         ! Sampling period in hours for means, etc
     *, A_IO               ! The output code from the unit command.
C

! I/O buffer workspace arrays - used for holding the data to be
! written to disk.

      REAL
     &  buf(global_LENOUT)
     &, buf3(global_LENOUT)

CDIR$ CACHE_ALIGN buf,buf3

      INTEGER
     *  IPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP) ! INTEGER LOOKUP TABLE
     *, N_ROWS_OUT         ! No of rows used for a diagnostic
     *, N_COLS_OUT         ! No of cols used PPHORIZ=N_ROWS*N_COLS_
     *, SROW_IN,SROW_OUT   ! North, South, East and West
     *, NROW_IN,NROW_OUT   ! subdomain limits in the horizontal sense
     *, WCOL_IN,WCOL_OUT   ! corresponding to the subarea before being
     *, ECOL_IN,ECOL_OUT   ! processed (IN) and after processing (OUT)
     *, LEV_IN_ADDR        ! The num of pts skipped in the Input
     *, GR                 ! Grid point code
     *, LBPROC_COMP(14)    ! Array of 1/0 denoting LBPROC components
     *, UNITPP             ! PPinit number (also used in PP_FILE)
     *, LENBUF             ! PPHORIZ_OUT rnd to 512 words (used PP)
     *, COMP_ACCRCY        ! PACKING ACCURACY IN POWER OF 2
     *, PPHORIZ_OUT        ! No of points in the output field
     *, PPHORIZ_IN         ! No of points in the input field
     *, IWA                ! Record no used inSETPOS
     *, LEN_BUF_WORDS      ! Number of 64 Bit words (rounded to 512)
     *, NUM_LEVS_OUT       ! Number of output levels
     *, NUM_LEVS_IN        ! Number of input levels
     *, NI                 ! Number of the INPUT STASH_LIST
     *, NO                 ! Number of the OUTPUT STASH_LIST
     *, INDX1
     *, INDEX_LEV(MAX_STASH_LEVS)    ! Index used to relate input and
C                                    ! output levels
     *, level_list(MAX_STASH_LEVS)   ! model level for each output level
     *, pseudo_level(MAX_STASH_LEVS) ! pseudo-level at each output level
     *, lv                           ! LV code for field from PP_XREF
     *, samples                      ! no of samples (timeseries/trajec)
     *, icurrll_dump_ptr             ! pointer to mother record LOOKUP
     *, start_time(7)                ! start time for PPheader in
C                                    ! timeseries/time mean/accumulation
     *, no_records           ! no of records processed by multi_spatial
     *, record_start         ! the start record for multi_spatial
     &, PEXNER               ! Exner Pressure
     &, PSTAR                ! Surface pressure
C
cdir$ cache_align ipplook
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

      INTEGER
     &        JL,II,IL,JJ,IT,
     &        ntab,        ! Number of the STASH TIMES table
     &        IKOUNT,      ! Local  counter
     &        POINTS,      ! No of points in a field
     &        POSN,        ! Local  temp variable
     &        KL,          ! Local  counter
     &        ML,          ! Local  counter
     *        LEN_IO,      ! The length of data transferred.
     *        ILPREV,      !The counter of the first of a pair of STLIST
     *        ILCURR,      ! The current value of IL
     *        IWL,         ! The word address of the LOOKUP table
     *        LBVCL,       ! Vertical coordinate code in LOOKUP table
     *        ICURRLL,     ! Current position in the PP Lookup table
     *        INT_LEVEL,   ! Integer copy of level
     *        DUMP_PACKING,  ! Copy of packing indicator in dump LOOKUP
     *        DUMP_DATA_TYPE,! Copy of data type in dump LOOKUP
     *        I              ! loop variable
     *       ,J            ! Level indicator used in call to GRIB_FILE
     *       ,PACKING_TYPE ! 0 No packing, 1 for WGDOS, 3 for GRIB
     &       ,GRIB_PACKING ! ppxref profile number used to determine
                           ! grib packing method
      INTEGER VX,VY,VZ     ! SIZES OF ARRAYS.
     &,       st_grid      ! Horizontal grid type, (eg. T-p or u-v)
     *,LEN_PPNAME
      INTEGER INPUT_CODE   ! VALUE OF INPUT_CODE IN STASHLIST
      INTEGER ADDR         ! ADDRESS OF STASH VARIABLE IN EITHER DUMP OR
      INTEGER ADDR_OUT     ! ADDRESS OF SPATIALLY PROCESSED FIELD
      INTEGER ELAP_TIME    ! NO OF TIMESTEPS ELAPSED IN PERIOD.
      INTEGER SERIES_PTR   ! THE ADDRESS IN STASH_SERIES WHERE DOMIN INF
      INTEGER INDEX_SIZE   ! THE NUMBER OF LEVELS IN THE INDEX.
      INTEGER BASE_LEVEL   ! Base model level needed for mass weighting
      INTEGER top_level    ! Top model level for 3D ocean decompress
      INTEGER base_level0,top_level0 ! Ref base/top levels in levs loop
      INTEGER what_proc    ! What kind of processing will be done
      INTEGER what_mean    ! What kind of meaning will be done
      INTEGER output_code  ! Output destination code from STLIST
      INTEGER expected_len         ! expected length of output field
      INTEGER ocnlev_bottom  ! first ocean level diagnosed
      LOGICAL
     &        S_F            ! TRUE for items requiring processing
     &,       OCEAN          ! TRUE if processing an ocean diagnostic
     &,       LLPROC         ! TRUE if spatial processing required
     &,       lnullproc      ! TRUE if null processing indicated
     &,       lfullfield     ! TRUE if output field on full horiz domain
     &,       lmasswt        ! TRUE if level-by-level mass-weights exist
     &,       start_step     ! TRUE at start of time period (TEMPORAL)
     &,       end_step       ! TRUE at end of time period (TEMPORAL)
     &,       MOS            ! TRUE if MOS output is required
     &,       PACKING        ! TRUE if packing required
     &,       GRIB_OUT     ! TRUE if output to be in GRIB code.
     &,       ROTATE         ! TRUE if input data to be rotated
      CHARACTER*14
     &        PPNAME           ! PPfile name decoded from MODEL_FT_UNIT
      CHARACTER*80
     &        STRING         ! PPfile name decoded from MODEL_FT_UNIT

      integer expected_extra ! expected length of extra data
      INTEGER extraw        ! number of extra words this timestep
      INTEGER extraw_hdr    ! number of extra words for the header
      INTEGER data_type_code ! ppx_data_type code copied from PPX file
      INTEGER rotatecode   ! code for rotated grid
      INTEGER NT_DIM         ! Number of tracers
      INTEGER pointer_dummy  ! dummy pointer variable for ocean
       REAL RPPLOOK(64)

      INTEGER
! local versions of the subdomain limits
     &  local_NROW_OUT,local_SROW_OUT,local_ECOL_OUT,local_WCOL_OUT
     &, local_NROW_IN,local_SROW_IN,local_ECOL_IN,local_WCOL_IN
! global versions of the total size of output
     &, global_N_ROWS_OUT,global_N_COLS_OUT, global_PPHORIZ_OUT
! MOS variables
     &, global_NROWS
     &, info ! return variable from GCOM


      INTEGER
     & grid_type_code  ! grid type of field being processed

CL----------------------------------------------------------------------
CL 0. Initialise: set constants relating to input grid type and size
CL
CL 0.1  Set up internal grid type st_grid and input field size
CL      according to the master GR code for the diagnostic
CL

! Get internal model index
      im_index = internal_model_index(im_ident)
        NT_DIM = (sa_idx(2) - sa_idx(1))/2

      if (im_ident .eq. ocean_im) then
        OCEAN = .true.
      else
        OCEAN = .false.
      endif

      if ((im_ident .eq. atmos_im) .or. (im_ident .eq. slab_im)) then
        pexner = stsuparr(sa_idx(7))
        pstar  = stsuparr(sa_idx(8))
      else
        pointer_dummy = 1
        pexner = pointer_dummy
        pstar  = pointer_dummy
      endif


! Get PP_XREF gridtype code
      GR = EXPPXI( im_ident, IS, IM, ppx_grid_type,
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
     &             icode, cmessage)
      grid_type_code=GR

      IF (GR.EQ.ppx_atm_tall.OR.GR.EQ.ppx_atm_tland.OR.
     &    GR.EQ.ppx_atm_tsea) THEN
C Atmosphere data on T-grid
        st_grid=st_tp_grid
        PPHORIZ_IN=T_field
      ELSEIF (GR.EQ.ppx_atm_uall.OR.GR.EQ.ppx_atm_uland.OR.
     &        GR.EQ.ppx_atm_usea) THEN
C Atmosphere data on U-grid
        st_grid=st_uv_grid
        PPHORIZ_IN=U_FIELD
      ELSEIF(GR.EQ.ppx_atm_compressed) THEN
C Atmosphere data on T-grid (compressed)
        st_grid=st_tp_grid
        PPHORIZ_IN=T_field
      ELSEIF(GR.EQ.ppx_atm_cuall) THEN
C Atmosphere data on C-grid (u-points)
        st_grid=st_cu_grid
        PPHORIZ_IN=T_field
      ELSEIF(GR.EQ.ppx_atm_cvall) THEN
C Atmosphere data on C-grid (v-points)
        st_grid=st_cv_grid
        PPHORIZ_IN=U_field
      ELSEIF(GR.EQ.ppx_atm_tzonal) THEN
C Atmosphere zonal data on T-grid
        st_grid=st_zt_grid
        PPHORIZ_IN=T_rows
      ELSEIF(GR.EQ.ppx_atm_uzonal) THEN
C Atmosphere zonal data on u-grid
        st_grid=st_zu_grid
        PPHORIZ_IN=u_rows
      ELSEIF(GR.EQ.ppx_atm_tmerid) THEN
C Atmosphere meridional data on T-grid
        st_grid=st_mt_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR.EQ.ppx_atm_umerid) THEN
C Atmosphere meridional data on u-grid
        st_grid=st_mu_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR.EQ.ppx_atm_scalar) THEN
C Atmosphere scalar
        st_grid=st_scalar
        PPHORIZ_IN=1
      ELSEIF (GR.EQ.ppx_ocn_tcomp.OR.GR.EQ.ppx_ocn_tall.OR.
     &        GR.EQ.ppx_ocn_tfield) THEN
C Ocean data on T-grid (compressed/uncompressed)
        st_grid=st_tp_grid
        PPHORIZ_IN=T_field
      ELSEIF (GR.EQ.ppx_ocn_ucomp.OR.GR.EQ.ppx_ocn_uall.OR.
     &        GR.EQ.ppx_ocn_ufield) THEN
C Ocean data on U-grid (compressed/uncompressed)
        st_grid=st_uv_grid
        PPHORIZ_IN=U_FIELD
      ELSEIF(GR.EQ.ppx_ocn_cuall) THEN
C Ocean data on C-grid (u-points)
        st_grid=st_cu_grid
        PPHORIZ_IN=T_field
      ELSEIF(GR.EQ.ppx_ocn_cvall) THEN
C Ocean data on C-grid (v-points)
        st_grid=st_cv_grid
        PPHORIZ_IN=U_field
      ELSEIF(GR.EQ.ppx_ocn_tzonal) THEN
C Ocean zonal data on T-grid
        st_grid=st_zt_grid
        PPHORIZ_IN=T_rows
      ELSEIF(GR.EQ.ppx_ocn_uzonal) THEN
C Ocean zonal data on u-grid
        st_grid=st_zu_grid
        PPHORIZ_IN=u_rows
      ELSEIF(GR.EQ.ppx_ocn_tmerid) THEN
C Ocean meridional data on T-grid
        st_grid=st_mt_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR.EQ.ppx_ocn_umerid) THEN
C Ocean meridional data on u-grid
        st_grid=st_mu_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR.EQ.ppx_ocn_scalar) THEN
C Ocean scalar
        st_grid=st_scalar
        PPHORIZ_IN=1
      ELSE
C Unknown grid type
        ICODE=1
        CMESSAGE='STWORK   : Unknown grid type found in PP_XREF'
        GOTO 999
      ENDIF
CL
CL 0.2 Set up ROTATE to flag fields which are rotated (eg. ELF winds)
CL     (this is used to set alternative fieldcodes in PPHEAD)
CL

      rotatecode = EXPPXI( im_ident, IS, IM, ppx_rotate_code,
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
     &             icode, cmessage)

      IF (rotatecode .EQ. ppx_elf_rotated .AND. ELF)
     &  THEN
        ROTATE=.TRUE.
      ELSE
        ROTATE=.FALSE.
      ENDIF
CL----------------------------------------------------------------------
CL 1. Loop over entries with this section/item number
CL
      DO 200 IL=ILSTART,ILEND  !loop over num entries for each item/sec
      extraw=0 ! no extra data by default
  104 FORMAT('  MAIN LOOP IN STWORK IL/IS/IM ARE',I6,2X,I4,2X,I4)
C Set MOS flag if output PP unit indicates output to MOS file
      MOS=(STLIST(st_output_code,IL).EQ.-89)
CL
CL 1.1 Set up S_F which has to be set for each STASHLIST entry. The
CL     STASHFLAG is be set for a particular ITEM/SECTION the S_F for
CL     the STASHLIST entry.
CL
      S_F=.FALSE.
      IF (STLIST(st_freq_code,IL).EQ.1.AND.
     &    STEP.GE.STLIST(st_start_time_code,IL).AND.
     &   (STEP.LE.STLIST(st_end_time_code,IL).OR.
     &    STLIST(st_end_time_code,IL).EQ.st_infinite_time)) THEN
C       ... if required every step between start and end
        S_F=.TRUE.
      ELSEIF(STLIST(st_freq_code,IL).LT.0) THEN
C       ... if required at specified times and this is one of them
        NTAB=-STLIST(st_freq_code,IL)
        DO 220 IT=1,NSTTIMS
          IF (STTABL(IT,NTAB).EQ.st_end_of_list) GOTO 230
          IF (STEP.EQ.STTABL(IT,NTAB)) S_F=.TRUE.
  220   CONTINUE
      ELSEIF (STLIST(st_freq_code,IL).gt.0) THEN
        IF   (MOD((STEP-STLIST(st_start_time_code,IL)),
     &             STLIST(st_freq_code,IL)).EQ.0.AND.
     &        STEP.GE.STLIST(st_start_time_code,IL).AND.
     &       (STEP.LE.STLIST(st_end_time_code,IL).OR.
     &        STLIST(st_end_time_code,IL).EQ.st_infinite_time))
C       ... if required every N timesteps and this is one of them
     &  S_F=.TRUE.
      ENDIF
  230 CONTINUE
C
C  S_F now set - Start of IF (S_F) block .......
C
      IF(S_F) THEN
CL
CL 1.2 Find number of input and output levels and relative positions
CL     and set up levels and pseudo-levels arrays for PPheaders.
CL     Set indicator lmasswt if level-by-level mass weighting possible
CL     - only currently available with atmosphere model full levels.
CL
! special case of mean timeseries leave ilcurr pointing to il           
                                                                        
        ilcurr=il   ! The current STASHlist entry IL
        IF (STLIST(st_input_code,IL).LT.0.and.                          
     &      STLIST(st_proc_no_code,IL).ne.st_time_series_mean) THEN 
          ilcurr=-STLIST(st_input_code,il) ! points to prev entry
        ENDIF
C

C Get PP_XREF lbvc code
      lbvcl = EXPPXI( im_ident, IS, IM, ppx_lbvc_code,
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
     &             icode, cmessage)

C Get PP_XREF lv code
      lv = EXPPXI( im_ident, IS, IM, ppx_lv_code,
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
     &             icode, cmessage)

C
        IF (lv.EQ.ppx_full_level.AND.im_ident .NE. ocean_im) THEN
          lmasswt=.TRUE.
        ELSE
          lmasswt=.FALSE.
        ENDIF
C
        IF (lv.EQ.ppx_half_level) THEN
          CALL STLEVELS(stlist(1,ilcurr),len_stlist,
     +      stash_levels,num_stash_levels,num_level_lists,
     +      stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,
     +      max_stash_levs,num_levs_in,num_levs_out,index_size,
     +      index_lev,level_list,
     &      lbvcl,stsuparr(sa_idx(3)),stsuparr(sa_idx(4)), ! akh,bkh
     &      level,pseudo_level,ak_lev,bk_lev,
     +      icode,cmessage)
        ELSE
          CALL STLEVELS(stlist(1,ilcurr),len_stlist,
     +      stash_levels,num_stash_levels,num_level_lists,
     +      stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,
     +      max_stash_levs,num_levs_in,num_levs_out,index_size,
     +      index_lev,level_list,
     &      lbvcl,stsuparr(sa_idx(3)),stsuparr(sa_idx(4)), ! akh,bkh
     &      level,pseudo_level,ak_lev,bk_lev,
     +      icode,cmessage)
        ENDIF
        IF (icode.gt.0) goto 999
        vz=num_levs_in
CL
CL 1.3 Find the horizontal dimensions for the output grid
CL     and the input field subdomain limits in which processing happens.
CL
        WCOL_IN= STLIST(st_west_code,ILCURR)    ! Input subdomain limits
        ECOL_IN= STLIST(st_east_code,ILCURR)
        NROW_IN= STLIST(st_north_code,ILCURR)
        SROW_IN= STLIST(st_south_code,ILCURR)
C
        IF(MOS) THEN
CL
CL 1.3.1 MOS output uses hard-wired MOS_OUTPUT_LENGTH, samples,extraw=0
CL
          PPHORIZ_OUT=MOS_OUTPUT_LENGTH
          NROW_OUT=1   ! Used in PP_HEAD , so set to preserve orig grid
          WCOL_OUT=1   !  "    "    "       "  "   "    "      "     "
          N_ROWS_OUT=1
          N_COLS_OUT=MOS_OUTPUT_LENGTH
          samples=0
          extraw_hdr=0
        ELSE
CL
CL 1.3.2 Other output types need to calculate lengths in detail
CL       (ie. number of rows, columns, and horizontal field size)
CL       according to processing options
CL

! Calculate local versions of the subdomain limits and area

          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(
     &      .TRUE.,.TRUE.,
     &      grid_type_code,mype,
     &      NROW_IN,ECOL_IN,SROW_IN,WCOL_IN,
     &      local_NROW_IN,local_ECOL_IN,
     &      local_SROW_IN,local_WCOL_IN)

          what_proc=stlist(st_proc_no_code,ilcurr)
          what_mean=(stlist(st_gridpoint_code,ilcurr)/block_size)*
     +       block_size
          samples=0         ! Initialise value for non-timeseries
CL
CL 1.3.2.1 Time series or trajectory processing
CL
          IF (what_proc.eq.st_time_series_code.or.
     &        what_proc.eq.st_append_traj_code.or.        
     &        what_proc.eq.st_time_series_mean) THEN                    
CL
CL 1.3.2.2 Compute number of samples in period for timeseries for
CL         input to PP_HEAD.
CL         No of output rows and cols are set to the no of points
CL         in each time sample and number of time samples in the
CL         period spanned by the output field, respectively.
CL
            samples=stlist(st_period_code,ilcurr)/
     &             stlist(st_freq_code,ilcurr)
            wcol_out=1
            ecol_out=samples
            nrow_out=1
            srow_out=stlist(st_output_length,ilcurr)/samples
            local_WCOL_OUT=1
            local_ECOL_OUT=samples
            local_NROW_OUT=1
            local_SROW_OUT=stlist(st_output_length,ilcurr)/samples
CL 1.3.2.3 Multi spatial processing of some other type (not supported)
          ELSEIF (stlist(st_series_ptr,ilcurr).ne.0) THEN
            ICODE=1323
            CMESSAGE='STWORK   : Illegal timeseries processing selected'
            GOTO 999
CL 1.3.2.4 Primary record requesting an extract
          ELSEIF (what_mean.eq.extract_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN
            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN
CL 1.3.2.5 Primary record requesting a vertical mean
          ELSEIF (what_mean.eq.vert_mean_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN
            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN
CL 1.3.2.6 Primary record requesting a zonal mean
          ELSEIF (what_mean.eq.zonal_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN
            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN
CL 1.3.2.7 Primary record requesting a meridional mean
          ELSEIF (what_mean.eq.merid_mean_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= 1
            SROW_OUT= 1
            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= 1
            local_SROW_OUT= 1
CL 1.3.2.8 Primary record requesting a global mean
          ELSEIF (what_mean.eq.global_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= 1
            SROW_OUT= 1
            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= 1
            local_SROW_OUT= 1
CL 1.3.2.9 Primary record requesting a field mean
          ELSEIF (what_mean.eq.field_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= 1
            SROW_OUT= 1
            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= 1
            local_SROW_OUT= 1
CL 1.3.2.10 Error trap for unknown request
          ELSE         ! Invalid option
            icode=st_unknown
            write(cmessage,87)'unknown option in setup',what_mean
            goto 999 ! jump to error return
          ENDIF
CL
CL 1.3.3 Compute expected length. This differs from total output length
CL       when data is appended from multiple timesteps into the same
CL       field, being output_length/number_of_appends in this case.
CL
          IF (stlist(st_output_code,il).ge.0.and.
     &        (what_proc.eq.st_time_series_code.or.                     
     &        what_proc.eq.st_append_traj_code.or.                      
     &        what_proc.eq.st_time_series_mean)) THEN   
          series_ptr=stlist(st_series_ptr,il) !set up ptr to stashseries
            expected_extra=(stash_series_index(2,series_ptr)+1)*6
            extraw_hdr=expected_extra
            expected_len=((stlist(st_output_length,ilcurr)
     &        -expected_extra)*
     &      stlist(st_freq_code,ilcurr))/stlist(st_period_code,ilcurr)
          ELSE
            expected_len=stlist(st_output_length,ilcurr)
            expected_extra=0 ! no extra data for non timeseries stuff
            extraw_hdr=0
          ENDIF
CL
CL 1.3.6 Compute number of rows and columns and field size for output
CL       - first adjust easternmost column if field wraps EW
CL
          IF (WCOL_IN .GT.ECOL_IN .AND.LCYCLIC)
     &      ECOL_IN =ECOL_IN + glsize(1)
          IF (WCOL_OUT.GT.ECOL_OUT.AND.LCYCLIC)
     &      ECOL_OUT=ECOL_OUT + glsize(1)
C

          IF (local_WCOL_OUT .GT. local_ECOL_OUT)
     &      local_ECOL_OUT=local_ECOL_OUT+ROW_LENGTH-2*Offx

          N_ROWS_OUT = local_SROW_OUT - local_NROW_OUT + 1
          N_COLS_OUT = local_ECOL_OUT - local_WCOL_OUT + 1
          global_N_ROWS_OUT = SROW_OUT - NROW_OUT + 1
          global_N_COLS_OUT = ECOL_OUT - WCOL_OUT + 1

          PPHORIZ_OUT= N_ROWS_OUT*N_COLS_OUT
          global_PPHORIZ_OUT=global_N_ROWS_OUT*global_N_COLS_OUT

        ENDIF     !  End of MOS IF block
CL
CL 1.4 Check to see if any processing is required.
CL     Set flag LLPROC if some SPATIAL processing indicated.
CL     NB: If input and output bottom levels differ (or the input and
CL         output pseudo-levels lists differ), level-by-level
CL         processing in the SPATIAL loop IS required.
CL         MULTI-SPATIAL processing is always required for timeseries.
CL
        lfullfield=((st_grid.EQ.st_tp_grid .OR. st_grid.EQ.st_cu_grid)
     &       .AND.  stlist(st_west_code,il).eq.1.and.
     &              stlist(st_east_code,il).eq.glsize(1).and.
     &              stlist(st_north_code,il).eq.1.and.
     &              stlist(st_south_code,il).eq.glsize(2)) .OR.
     &             ((st_grid.EQ.st_uv_grid .OR. st_grid.EQ.st_cv_grid)
     &       .AND.  stlist(st_west_code,il).eq.1.and.
     &              stlist(st_east_code,il).eq.glsize(1).and.
     &              stlist(st_north_code,il).eq.1.and.
     &              stlist(st_south_code,il).eq.glsize(2)-1) .OR.
     &             ((st_grid.EQ.st_zt_grid)
     &       .AND.  stlist(st_north_code,il).eq.1.and.
     &              stlist(st_south_code,il).eq.glsize(2)) .OR.
     &             ((st_grid.EQ.st_zu_grid)
     &       .AND.  stlist(st_north_code,il).eq.1.and.
     &              stlist(st_south_code,il).eq.glsize(2)-1) .OR.
     &             ((st_grid.EQ.st_mt_grid .OR. st_grid.EQ.st_mu_grid)
     &       .AND.  stlist(st_west_code,il).eq.1.and.
     &              stlist(st_east_code,il).eq.glsize(1)) .OR.
     &             (st_grid.EQ.st_scalar)
        lnullproc= lfullfield .AND.
     &             (stlist(st_input_bottom,il).eq.
     &                stlist(st_output_bottom,il)) .and.
     &             (stlist(st_pseudo_in,il).eq.
     &                stlist(st_pseudo_out,il)) .and.
     &             (stlist(st_gridpoint_code,il).eq.
     &               (extract_base+stash_null_mask_code) .and.
     &              stlist(st_weight_code,IL).eq.
     &                stash_weight_null_code )
        IF (STLIST(st_series_ptr,IL).GT.0) THEN
          lnullproc=.FALSE.     ! Timeseries always requires processing
        ENDIF
C  LLPROC must be false for MOS output, output from a prev STLIST
C  or simple extraction of full field with no weighting
        IF (MOS .or. (STLIST(st_input_code,IL).LT.0.and.                
     &      STLIST(st_proc_no_code,IL).ne.st_time_series_mean).or.
     &      lnullproc) THEN          
          LLPROC=.FALSE.
        ELSE
          LLPROC=.TRUE.
        ENDIF
CL
CL 1.5 Check that no spatial processing is requested if the input field
CL     is of integer or logical type -- these types of fields can be
CL     passed directly through STASH, for example for coupling purposes,
CL     but no arithmetic is allowed at present.
CL
        data_type_code=EXPPXI(im_ident,IS,IM,ppx_data_type,
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
     &                        icode,cmessage)
C
        IF ((       data_type_code .EQ.ppx_type_int    .OR.
     &              data_type_code .EQ.ppx_type_log)   .AND.
     &       .NOT. lnullproc)  THEN
          ICODE=st_not_supported
          CMESSAGE='STWORK  : Spatial processing of INT/LOGICAL illegal'
          GOTO 999
        ENDIF
CL----------------------------------------------------------------------
CL 2. Perform spatial processing (loop over output levels)
CL
        IF (LLPROC) THEN   ! Processing is required
          input_code=stlist(st_input_code,il)
CL make sure no volume processing asked for as not supported
CL this will need adding at some point
          IF (stlist(st_weight_code,il).eq.stash_weight_volume_code)
     +      THEN
            icode=st_not_supported
            cmessage='STWORK  : volume processing not supported'
            goto 999
          ENDIF
C Work out vx,vy (depends on kind of grid)
          IF (st_grid.EQ.st_tp_grid.OR.st_grid.EQ.st_cu_grid) THEN
            vx=row_length            ! T_rows x row_length
            vy=T_rows                ! worth of data in input field
          ELSEIF (st_grid.EQ.st_uv_grid.OR.st_grid.EQ.st_cv_grid) THEN
            vx=row_length            ! u_rows x row_length
            vy=u_rows                ! worth of data in input field
          ELSEIF (st_grid.EQ.st_zt_grid) THEN
            vx=1                     ! T_rows x 1
            vy=T_rows                ! worth of data in input field
          ELSEIF (st_grid.EQ.st_zu_grid) THEN
            vx=1                     ! u_rows x 1
            vy=u_rows                ! worth of data in input field
          ELSEIF (st_grid.EQ.st_mt_grid .OR. st_grid.EQ.st_mu_grid) THEN
            vx=row_length            ! 1 x row_length
            vy=1                     ! worth of data in input field
          ELSEIF (st_grid.EQ.st_scalar) THEN
            vx=1                     ! 1 x 1
            vy=1                     ! worth of data in input field
          ENDIF
CL Work out if this is the first timestep in a timeseries.
CL This is required so that the extra data can be generated
C
          series_ptr=stlist(st_series_ptr,il)
          IF (series_ptr.gt.0) THEN ! multi spatial processing reqd.
CL recompute expected sizes
            elap_time=step-stlist(st_start_time_code,il)
            elap_time=mod(elap_time,stlist(st_period_code,il))
            start_step=(elap_time.eq.0)
            IF (start_step) THEN
              expected_len=stlist(st_output_length,ilcurr)
            ELSE
              expected_extra=0  ! reset to zero as no extra data
              expected_len=((stlist(st_output_length,ilcurr)
     &          -((stash_series_index(2,series_ptr)+1)*6))*
     &        stlist(st_freq_code,ilcurr))/stlist(st_period_code,ilcurr)
            ENDIF
CL
CL 2.1 Timeseries extraction section follows
CL
            no_records=stash_series_index(2,series_ptr)
            record_start=stash_series_index(1,series_ptr)
CL
CL 2.1.0 Strip/decompress ocean fields using STOCGT if not already
CL       processed - output result to ocwork
CL       processed - output result to ocwork
CL       Similarly decompress wave fields using STWVGT.
CL
            IF (im_ident .eq. ocean_im) THEN
              base_level=stlist(st_input_bottom,il)
              top_level=stlist(st_input_top,il)
              ocnlev_bottom=base_level
              IF (base_level.LT.0.OR.base_level.EQ.st_special_code)THEN
                base_level=1
                top_level =1
                ocnlev_bottom=1
              ENDIF
              IF (input_code.eq.0) THEN
                CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            d1,im,is,base_level,top_level,                        
     &            ocnlev_bottom,rmdi,ocwork,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                IF (icode.GT.0) GOTO 999
              ELSEIF (input_code.eq.1) THEN
                CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            stash_work,im,is,base_level,top_level,
     &            ocnlev_bottom,rmdi,ocwork,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                IF (icode.GT.0) GOTO 999
              ENDIF

            ELSEIF (im_ident .eq. wave_im) THEN

              IF (input_code.eq.0) THEN
                CALL STWVGT(
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
     &            t_rows,row_length,im_ident,sm_ident,
     &            d1,im,is,rmdi,ocwork,vx,vy,si,
     &            istsuparr(sa_idx(1)),               ! land-sea mask
     &            ICODE,CMESSAGE)
                IF (icode.GT.0) GOTO 999

              ELSEIF (input_code.eq.1) THEN
                CALL STWVGT(
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
     &            t_rows,row_length,im_ident,sm_ident,
     &            stash_work,im,is,rmdi,ocwork,vx,vy,si,
     &            istsuparr(sa_idx(1)),               ! land-sea mask
     &            ICODE,CMESSAGE)
                IF (icode.GT.0) GOTO 999
              ENDIF
            ENDIF
CL
CL 2.1.1 Process a primary field from D1 (timeseries)
CL
            IF (input_code.eq.0) THEN
              IF (im_ident .eq.ocean_im) THEN
                CALL MULTI_SPATIAL(ocwork,
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt,
     +            pphoriz_out,num_levs_out,
     &            d1(pexner), d1(pstar),
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude
     &            stsuparr(sa_idx(11)),       ! land
     +            row_length,T_rows,u_rows,T_levels,
     +            ppfield,lenout,
     +            rmdi,stlist(1,il),len_stlist,
     +            stash_series(1,record_start),
     +            stash_series_rec_len,no_records,
     +            index_size,index_lev,level_list,
     +            start_step,extraw,n_rows_out,n_cols_out,
     +            realhd,len_realhd,inthd,len_inthd,ocean,
     +            icode,cmessage)
              ELSE
                 addr=si(im,is,im_index)
                CALL MULTI_SPATIAL(d1(addr),
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt,
     +            pphoriz_out,num_levs_out,
     &            d1(pexner), d1(pstar),
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude
     &            stsuparr(sa_idx(11)),       ! land
     +            row_length,T_rows,u_rows,T_levels,
     +            ppfield,lenout,
     +            rmdi,stlist(1,il),len_stlist,
     +            stash_series(1,record_start),
     +            stash_series_rec_len,no_records,
     +            index_size,index_lev,level_list,
     +            start_step,extraw,n_rows_out,n_cols_out,
     +            realhd,len_realhd,inthd,len_inthd,ocean,
     +            icode,cmessage)
              ENDIF
CL
CL 2.1.2 Process a field from STASHWORK (timeseries)
CL
            ELSEIF (input_code.eq.1) THEN
              IF (im_ident .eq.ocean_im) THEN
                CALL MULTI_SPATIAL(ocwork,
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt,
     +            pphoriz_out,num_levs_out,
     &            d1(pexner), d1(pstar),
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude
     &            stsuparr(sa_idx(11)),       ! land
     +            row_length,T_rows,u_rows,T_levels,
     +            ppfield,lenout,
     +            rmdi,stlist(1,il),len_stlist,
     +            stash_series(1,record_start),
     +            stash_series_rec_len,no_records,
     +            index_size,index_lev,level_list,
     +            start_step,extraw,n_rows_out,n_cols_out,
     +            realhd,len_realhd,inthd,len_inthd,ocean,
     +            icode,cmessage)
              ELSE
                 addr=si(im,is,im_index)
                CALL MULTI_SPATIAL(stash_work(addr),
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt,
     +            pphoriz_out,num_levs_out,
     &            d1(pexner), d1(pstar),
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude
     &            stsuparr(sa_idx(11)),       ! land
     +            row_length,T_rows,u_rows,T_levels,
     +            ppfield,lenout,
     +            rmdi,stlist(1,il),len_stlist,
     +            stash_series(1,record_start),
     +            stash_series_rec_len,no_records,
     +            index_size,index_lev,level_list,
     +            start_step,extraw,n_rows_out,n_cols_out,
     +            realhd,len_realhd,inthd,len_inthd,ocean,
     +            icode,cmessage)
              ENDIF
            ELSEIF (input_code.lt.0) then
CL
CL 2.1.3 Process a field from previously STASHed position in D1
CL     (currently unsupported since diagnostic-of-diagnostic)
CL
             IF (what_proc.eq.st_time_series_mean) THEN  
! special case of mean timeseries                                       
              IF (im_ident .eq.ocean_im) THEN                           
! Warning never tested for ocean case                                   
                CALL MULTI_SPATIAL(ocwork,                              
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt, 
     +            pphoriz_out,num_levs_out,                             
     &            d1(pexner), d1(pstar),                                
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)    
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)    
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude          
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude          
     &            stsuparr(sa_idx(11)),       ! land                    
     +            row_length,T_rows,u_rows,T_levels,                    
     +            ppfield,lenout,                                       
     +            rmdi,stlist(1,il),len_stlist,                         
     +            stash_series(1,record_start),                         
     +            stash_series_rec_len,no_records,                      
     +            index_size,index_lev,level_list,                      
     +            start_step,extraw,n_rows_out,n_cols_out,              
     +            realhd,len_realhd,inthd,len_inthd,ocean,              
     +            icode,cmessage)                                       
             ELSE                                                       
!   Mother record                                                       
                 ILPREV=-stlist(st_input_code,IL)                       
! address of mother record in D1                                        
                 addr=stlist(20,ILPREV)                                 
                CALL MULTI_SPATIAL(D1(addr),                            
     +            vx,vy,vz,grid_type_code,st_grid,lcyclic,lmasswt,      
     +            pphoriz_out,num_levs_out,                             
     &            d1(pexner), d1(pstar),                                
     &            stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)    
     &            stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)    
     &            stsuparr(sa_idx(9)),        ! cos_p_latitude          
     &            stsuparr(sa_idx(10)),       ! cos_u_latitude          
     &            stsuparr(sa_idx(11)),       ! land                    
     +            row_length,T_rows,u_rows,T_levels,                    
     +            ppfield,lenout,                                       
     +            rmdi,stlist(1,il),len_stlist,                         
     +            stash_series(1,record_start),                         
     +            stash_series_rec_len,no_records,                      
     +            index_size,index_lev,level_list,                      
     +            start_step,extraw,n_rows_out,n_cols_out,              
     +            realhd,len_realhd,inthd,len_inthd,ocean,              
     +            icode,cmessage)                                       
              ENDIF                                                     
             ELSE                                                       
              icode=st_not_supported
              cmessage='STWORK1  : diag-of-diagnostic unsupported'
              goto 999 ! jump to error return
             ENDIF                                                      
            ELSE
              icode=st_unknown
              write(cmessage,87)'unknown input option',input_code
              goto 999
            ENDIF
            if (icode.ne.0) goto 999 ! error exit
            IF (start_step.AND.(extraw.NE.expected_extra)) THEN
              icode=st_bad_array_param
              write(cmessage,89) extraw,expected_extra
 89           FORMAT('STWORK : Inconsistent length for extra data ',
     &          i8,1x,i8)
              goto 999
            endif
C
CL         pphoriz_out has been computed by multi_spatial
CL         it is the size of the output field
          ELSE  ! do "normal" spatial processing
CL
CL 2.2 Standard spatial processing section follows
CL
CL In ocean case, a possibly 3D decompression is required (to top_level)
CL
CL If multi-level processing (ie. vertical, global mean) is performed by
CL SPATIAL, the input field is passed in with the original start address
CL but if single-level processing is done by SPATIAL, the field is
CL passed in with an address pointing to the single level required.
CL
            base_level0=stlist(st_input_bottom,il)
            what_proc=stlist(st_gridpoint_code,il)
            IF (im_ident .eq.ocean_im) THEN
              IF ((what_proc.LT.vert_mean_top .AND.
     &             what_proc.GT.vert_mean_base) .OR.
     &            (what_proc.LT.global_mean_top .AND.
     &             what_proc.GT.global_mean_base)) THEN
                top_level0=stlist(st_input_top,il)
              ELSE
                top_level0=base_level0
              ENDIF
            ENDIF
C
            addr_out=1                   ! Initialise output address
C
            DO kl=1,num_levs_out         ! --- Start of levels loop ---
C Work out model level if model level range, otherwise set to 1
              IF (base_level0.LT.0.OR.base_level0.EQ.st_special_code)
     &        THEN
                base_level=1
                IF (im_ident .eq.ocean_im) THEN
                  top_level =1
                  ocnlev_bottom=1
                ENDIF
              ELSE
                base_level=base_level0+index_lev(kl)-1
                IF (im_ident .eq.ocean_im) THEN
                  top_level=top_level0+index_lev(kl)-1
                  ocnlev_bottom=stlist(st_input_bottom,il)
                ENDIF
              ENDIF
CL
CL 2.2.0 Strip/decompress ocean fields using STOCGT if not already
CL       processed - output result to ocwork
CL       NB: If 3D global mean is required, perform 3D decompress
CL           between base_level and top_level.
CL
              IF (im_ident .eq. ocean_im) THEN
                IF (input_code.eq.0) THEN
                  CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            d1,im,is,base_level,top_level,                        
     &            ocnlev_bottom,rmdi,ocwork,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                  IF (icode.GT.0) GOTO 999
                ELSEIF (input_code.eq.1) THEN
                  CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            stash_work,im,is,base_level,top_level,
     &            ocnlev_bottom,rmdi,ocwork,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                  IF (icode.GT.0) GOTO 999
                ENDIF

              ELSEIF (im_ident .eq. wave_im) THEN

                IF (input_code.eq.0) THEN
                  CALL STWVGT(
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
     &              t_rows,row_length,im_ident,sm_ident,
     &              d1,im,is,rmdi,ocwork,vx,vy,si,
     &              istsuparr(sa_idx(1)),             ! land-sea mask
     &              ICODE,CMESSAGE)
                  IF (icode.GT.0) GOTO 999

                ELSEIF (input_code.eq.1) THEN
                  CALL STWVGT(
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
     &              t_rows,row_length,im_ident,sm_ident,
     &              stash_work,im,is,rmdi,ocwork,vx,vy,si,
     &              istsuparr(sa_idx(1)),             ! land-sea mask
     &              ICODE,CMESSAGE)
                  IF (icode.GT.0) GOTO 999

                ENDIF
              ENDIF
CL
CL 2.2.1 Process a primary field from D1 (or ocwork if OCEAN field)
CL
              IF (input_code.eq.0) THEN
                if (im_ident .eq. ocean_im) then
                  CALL SPATIAL(ocwork,vx,vy,vz,
     +                grid_type_code,st_grid,lcyclic,lmasswt,
     +                n_cols_out,n_rows_out,base_level,
     +                level_list,index_lev,index_size,
     &                d1(pexner), d1(pstar),
     &                stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &                stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &                stsuparr(sa_idx(9)),        ! cos_p_latitude
     &                stsuparr(sa_idx(10)),       ! cos_u_latitude
     &                stsuparr(sa_idx(11)),       ! land
     +                row_length,T_rows,u_rows,T_levels,
     +                ppfield(addr_out),pphoriz_out,
     +                stlist(1,il),len_stlist,rmdi,
     +                icode,cmessage)
                ELSE
                  IF ((what_proc.LT.vert_mean_top .AND.
     &                 what_proc.GT.vert_mean_base) .OR.
     &                (what_proc.LT.global_mean_top .AND.
     &                 what_proc.GT.global_mean_base)) THEN
                    addr=si(im,is,im_index)
                  ELSE
                    addr=si(im,is,im_index)+
     &                   (index_lev(kl)-1)*pphoriz_in
                  ENDIF
                  IF (addr.lt.1.or.addr.gt.len_tot) THEN
                    icode=st_bad_address
                    cmessage='STWORK  : D1 address out of bounds'
                    goto 999
                  ENDIF
                  CALL SPATIAL(d1(addr),vx,vy,vz,
     +                grid_type_code,st_grid,lcyclic,lmasswt,
     +                n_cols_out,n_rows_out,base_level,
     +                level_list,index_lev,index_size,
     &                d1(pexner), d1(pstar),
     &                stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &                stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &                stsuparr(sa_idx(9)),        ! cos_p_latitude
     &                stsuparr(sa_idx(10)),       ! cos_u_latitude
     &                stsuparr(sa_idx(11)),       ! land
     +                row_length,T_rows,u_rows,T_levels,
     +                ppfield(addr_out),pphoriz_out,
     +                stlist(1,il),len_stlist,rmdi,
     +                icode,cmessage)
                ENDIF
CL
CL 2.2.2 Process a field from STASHWORK (or ocwork if OCEAN)
CL
              ELSEIF (input_code.eq.1) THEN
                IF (im_ident .eq. ocean_im) THEN
                  CALL SPATIAL(ocwork,vx,vy,vz,
     +                grid_type_code,st_grid,lcyclic,lmasswt,
     +                n_cols_out,n_rows_out,base_level,
     +                level_list,index_lev,index_size,
     &                d1(pexner), d1(pstar),
     &                stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &                stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &                stsuparr(sa_idx(9)),        ! cos_p_latitude
     &                stsuparr(sa_idx(10)),       ! cos_u_latitude
     &                stsuparr(sa_idx(11)),       ! land
     +                row_length,T_rows,u_rows,T_levels,
     +                ppfield(addr_out),pphoriz_out,
     +                stlist(1,il),len_stlist,rmdi,
     +                icode,cmessage)
                ELSE
                  IF ((what_proc.LT.vert_mean_top .AND.
     &                 what_proc.GT.vert_mean_base) .OR.
     &                (what_proc.LT.global_mean_top .AND.
     &                 what_proc.GT.global_mean_base)) THEN
                    addr=si(im,is,im_index)
                  ELSE
                    addr=si(im,is,im_index)+
     &                   (index_lev(kl)-1)*pphoriz_in
                  ENDIF
                  IF (addr.lt.1.or.addr.gt.stash_work_len) THEN
                    icode=st_bad_address
                    cmessage='STWORK  : STASHWORK addr out of bounds'
                    goto 999
                  ENDIF
                  CALL SPATIAL(stash_work(addr),vx,vy,vz,
     +              grid_type_code,st_grid,lcyclic,lmasswt,
     +              n_cols_out,n_rows_out,base_level,
     +              level_list,index_lev,index_size,
     &              d1(pexner), d1(pstar),
     &              stsuparr(sa_idx(5)),        ! a_levdepc(jdelta_ak)
     &              stsuparr(sa_idx(6)),        ! a_levdepc(jdelta_bk)
     &              stsuparr(sa_idx(9)),        ! cos_p_latitude
     &              stsuparr(sa_idx(10)),       ! cos_u_latitude
     &              stsuparr(sa_idx(11)),       ! land
     +              row_length,T_rows,u_rows,T_levels,
     +              ppfield(addr_out),pphoriz_out,
     +              stlist(1,il),len_stlist,rmdi,
     +              icode,cmessage)
                ENDIF
              ELSEIF (input_code.lt.0) THEN
CL
CL 2.2.3 Process a field from previously STASHed position in D1
CL     (currently unsupported since diagnostic-of-diagnostic)
CL
                icode=st_not_supported
                cmessage='STWORK1  : diag-of-diagnostic unsupported'
                goto 999
              ELSE
                icode=st_unknown
                write(cmessage,87)'unknown input option',input_code
87              format('STWORK1 : >>FATAL ERROR <<',a,1x,i5)
                goto 999
              ENDIF
C
              IF (icode.gt.0) goto 999 ! Trap error
C
CL compute pphoriz_out
CL pphoriz_out is the size of the output vector
CL we should not be doing timeseries processing here.
CL
CL NOTE: n_cols_out and n_rows_out should agree with values calculated
CL       before, but are not checked for consistency.
C
              pphoriz_out=n_cols_out*n_rows_out
              addr_out=addr_out+pphoriz_out ! increment output address
            ENDDO                      ! --- End of levels loop ---
C
          ENDIF         ! End of multi-spatial/spatial IF block
C
          IF (icode.gt.0) goto 999     ! Trap processing error
CL
CL 2.3 Set length of output field and check against expected length
CL
CL check that extrawords are the same as the expected number of extrawor
! Calculate size of global pphoriz_out - the size on disk

          IF (what_proc .eq. st_time_series_code .OR.
     &        what_proc .eq. st_time_series_mean) THEN      
            global_pphoriz_out=pphoriz_out
            global_n_rows_out=n_rows_out
            global_n_cols_out=n_cols_out
          ELSE
            CALL STASH_GET_GLOBAL_SIZE(
     &       stlist(st_north_code,il) , stlist(st_east_code,il),
     &       stlist(st_south_code,il) , stlist(st_west_code,il),
     &       1,
     &       STLIST(st_gridpoint_code,il) , STLIST(st_proc_no_code,il),
     &       global_pphoriz_out,
     &       ICODE, CMESSAGE)

            IF (icode .ne. 0) goto 999

          ENDIF

          IF (pphoriz_out*num_levs_out.ne.expected_len) THEN
            icode=st_bad_array_param
            write(cmessage,88) pphoriz_out*num_levs_out,expected_len
 88         FORMAT('STWORK   : Inconsistent length for output field ',
     &             i8,1x,i8)
            goto 999
          ENDIF
        ELSE
CL----------------------------------------------------------------------
CL 3. No SPATIAL processing - extract output field by direct copy
CL
          IF(MOS) THEN ! Mos OUTPUT
CL
CL 3.1 MOS output.
CL
C Establish which grid the data is on ie wind/pts temp/pts ocean etc
            IF(st_grid.EQ.st_tp_grid.OR.st_grid.EQ.st_cu_grid) THEN
              POINTS=T_field
              global_NROWS=glsize(2)
            ELSE
              POINTS=U_FIELD
              global_NROWS=glsize(2)-1
            ENDIF
            IF(PPHORIZ_IN.NE.POINTS) THEN
              ICODE=1
              CMESSAGE='STWORK  : MOS input must be a complete field'
              GOTO 999
            ENDIF
            IKOUNT=0

! Send the required field to PE 0 (buf), and pack into PPFIELD
            DO KL=1,NUM_LEVS_OUT

! Gather the field into buf1 array on PE 0

              IF(STLIST(st_input_code,IL).LT.0) THEN ! previous entry

                II=-STLIST(st_input_code,IL)
                LEV_IN_ADDR=(KL-1)*PPHORIZ_IN

                CALL GATHER_FIELD(
     &            D1(STLIST(st_output_addr,II)+LEV_IN_ADDR),buf,
     &            ROW_LENGTH,T_ROWS,
     &            glsize(1),global_NROWS,
     &            0,gc_all_proc_group,
     &            info)

              ELSE ! Current STASH list
                LEV_IN_ADDR=(INDEX_LEV(KL)-1)*PPHORIZ_IN

                IF (STLIST(st_input_code,IL).EQ.0) THEN ! DATA in D1

                  CALL GATHER_FIELD(
     &            D1(SI(IM,IS,im_index)+LEV_IN_ADDR),buf,
     &            ROW_LENGTH,T_ROWS,
     &            glsize(1),global_NROWS,
     &            0,gc_all_proc_group,
     &            info)

                ELSE ! DATA in STWORK

                  CALL GATHER_FIELD(
     &            STASH_WORK(SI(IM,IS,im_index)+LEV_IN_ADDR),buf,
     &            ROW_LENGTH,T_ROWS,
     &            glsize(1),global_NROWS,
     &            0,gc_all_proc_group,
     &            info)

                ENDIF
              ENDIF

! PE 0 must now pack the field using the mask in MOS_MASK
! from the buf array into PPFIELD

              DO JL=1,glsize(1)*global_NROWS
                IF (MOS_MASK(JL).EQ.1) THEN
                  IKOUNT=IKOUNT+1

                  IF (mype .EQ. 0)
     &              PPFIELD(IKOUNT)=buf(JL)

                ENDIF
              ENDDO

            ENDDO  ! KL : loop over output levels

            global_PPHORIZ_OUT=PPHORIZ_OUT
            global_N_ROWS_OUT=N_ROWS_OUT
            global_N_COLS_OUT=N_COLS_OUT


              IF(IKOUNT/NUM_LEVS_OUT.NE.PPHORIZ_OUT)THEN
                WRITE(6,*)'MOS_OUTPUT_LENGTH  ',MOS_OUTPUT_LENGTH
                WRITE(6,*)'IKOUNT             ',IKOUNT
                ICODE=1
                CMESSAGE='STWORK  MOS_OUTPUT_LENGTH not = to MOS_MASK'
                GOTO 999
              ENDIF

          ELSE
CL
CL 3.2 Other output - determine input source
CL
            input_code=STLIST(st_input_code,IL)
CL
CL 3.2.1 Ocean fields need to be stripped/decompressed using STOCGT
CL       if not already processed, and output to ppfield array
CL       except when input is already STASHed in D1 (input_code.lt.0)
CL
            IF ( (im_ident .eq. ocean_im).AND.(input_code.eq.0) ) THEN
              base_level=stlist(st_input_bottom,il)
C Set base_level to 1 for special level fields
                IF (base_level.eq.st_special_code) base_level=1
                  CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            d1,im,is,base_level,base_level+num_levs_out-1,
     &            base_level,rmdi,ppfield,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                IF (icode.GT.0) GOTO 999
            ELSEIF ((im_ident.eq.ocean_im).AND.(input_code.eq.1)) THEN
                base_level=stlist(st_input_bottom,il)
C Set base_level to 1 for special level fields or input levels lists
                IF (base_level.LT.0.OR.base_level.EQ.st_special_code)
     &            base_level=1
                CALL STOCGT (
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
     &            t_rows,row_length,t_levels,
     &            im_ident,sm_ident,
     &            stash_work,im,is,base_level,
     &            base_level+num_levs_out-1,
     &            base_level,rmdi,ppfield,vx,vy,vz,
     &            nt_dim,si,istsuparr(sa_idx(1)),     ! joc_tracer 
     &            istsuparr(sa_idx(2)),               ! joc_u
     &            istsuparr(sa_idx(3)),               ! joc_v
     &            istsuparr(sa_idx(6)),               ! o_cfi1,
     &            istsuparr(sa_idx(7)),               ! o_cfi2
     &            istsuparr(sa_idx(8)),               ! o_cfi3
     &            istsuparr(sa_idx(4)),               ! joc_no_seapts
     &            istsuparr(sa_idx(5)),               ! joc_no_segs
     &            icode,cmessage)
                IF (icode.GT.0) GOTO 999
CL
CL 3.2.1.1 Wave fields need to be stripped/decompressed using STWVGT 
CL       if not already processed, and output to ppfield array
CL       except when input is already STASHed in D1 (input_code.lt.0)
CL
            ELSEIF ((im_ident.eq.wave_im).AND.(input_code.eq.0)) THEN
                CALL STWVGT(
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
     &            t_rows,row_length,im_ident,sm_ident,
     &            d1,im,is,rmdi,ppfield,vx,vy,si,
     &            istsuparr(sa_idx(1)),               ! land-sea mask
     &            ICODE,CMESSAGE)
                IF (icode.GT.0) GOTO 999

            ELSEIF ((im_ident.eq.wave_im).AND.(input_code.eq.1)) THEN
                CALL STWVGT(
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
     &            t_rows,row_length,im_ident,sm_ident,
     &            stash_work,im,is,rmdi,ppfield,vx,vy,si,
     &            istsuparr(sa_idx(1)),               ! land-sea mask
     &            ICODE,CMESSAGE)
                IF (icode.GT.0) GOTO 999

CL
CL 3.2.2 Other fields are simply copied
CL
            ELSEIF (input_code.eq.0) THEN
C Simple extraction with no weighting from primary field in D1
C except for those needing special extraction on funny grids
C (Ocean, Wave)
              addr=si(im,is,im_index)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=D1(addr+JL-1)
              ENDDO
            ELSEIF (input_code.eq.1) THEN
C Simple extraction with no weighting from STASH_WORK
C except for those needing special extraction on funny grids
C (Ocean, Wave)
              addr=si(im,is,im_index)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=STASH_WORK(addr+JL-1)
              ENDDO
            ELSEIF (input_code.lt.0) THEN
C Previously STASHed entry in D1
C for all sub-models as diagnostic D1 is always on a proper grid.

              addr=STLIST(st_output_addr,-input_code)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=D1(addr+JL-1)
              ENDDO
            ELSE
C Illegal input code
              ICODE=st_unknown
              CMESSAGE='STWORK   : Unknown input code encountered'
              GOTO 999
            ENDIF
          ENDIF    ! End of IF for test if MOS output
        ENDIF      ! End of LLPROC IF BLOCK    ************************

CL-----------------------------------------------------------------
CL 4. OUTPUT section.
CL
CL    The data is in PPFIELD with a length LENOUT.
CL    The horizontal field size PPHORIZ_OUT and number of output levels
CL    NUM_LEVS_OUT were calculated in section 1.
CL    Output option depends on the STLIST code.
CL
CL 4.0 Find mother STASHlist record if necessary.
CL
C   Packing_type not set from PP_FILE in the ELSE IF part of block.
            PACKING_TYPE = 0  ! Default is unpacked.
        IF(STLIST(st_input_code,IL).LT.0) THEN ! Second of two STLIST
          ILPREV=-STLIST(st_input_code,IL)
        ELSE
          ilprev=IL ! no daughter record
        ENDIF
CL
CL 4.0.1 Set up LBPROC sub-components based on STASH processing info.
CL
        DO JJ=1,14
          LBPROC_COMP(JJ)=0
        ENDDO
C
        IF(STLIST(st_gridpoint_code,ilprev).GE.zonal_mean_base)
     &    LBPROC_COMP(7)=1
C
        IF((STLIST(st_gridpoint_code,ilprev).GE.vert_mean_base .AND.
     &      STLIST(st_gridpoint_code,ilprev).LT.vert_mean_top) .OR.
     &    (STLIST(st_gridpoint_code,ilprev).GE.global_mean_base .AND.
     &     STLIST(st_gridpoint_code,ilprev).LT.global_mean_top))
     &    LBPROC_COMP(12)=1
C
        IF((STLIST(st_proc_no_code,ilprev).EQ.st_accum_code) .OR.
     &     (STLIST(st_proc_no_code,ilprev).EQ.st_time_mean_code).OR. 
     &     (STLIST(st_proc_no_code,ilprev).EQ.st_time_series_mean)) 
     &    LBPROC_COMP(8)=1
C
        IF(STLIST(st_proc_no_code,ilprev).EQ.st_min_code)
     &    LBPROC_COMP(13)=1
C
        IF(STLIST(st_proc_no_code,ilprev).EQ.st_max_code)
     &    LBPROC_COMP(14)=1
C
        output_code=STLIST(st_output_code,IL)
CL
CL 4.1 OUTPUT to PPfile
CL
        IF (output_code.LT.0) THEN                  ! PP Output
C Find appropriate dump header if a daughter record
          IF (il.ne.ilcurr) then
            icurrll_dump_ptr=stlist(st_lookup_ptr,ilcurr)
          ENDIF
CL
CL 4.1.0 Determine output PP unit and associated filename; OPEN file
CL
C If preattached files are used the file is left open by PPCTL following
C the initial OPEN; if reinitialised files are used the unit must be
C OPENed and CLOSEd explicitly every time it is used.
C
          UNITPP=-output_code
          IF (FT_STEPS(UNITPP).NE.0) THEN ! Filename generated by model
! Check if re-initialised file stream has been opened yet 
            IF(STEP.LT.FT_FIRSTSTEP(UNITPP)) THEN ! File stream opened?
               ICODE=1
               CMESSAGE='STWORK  : Re-initialised file not yet created'
               write(6,*)
     &          'STWORK  : FATAL ERROR. Attempt to write to ',
     &          're-initialised file stream before file first opened:'
               write(6,*)
     &          '        : Check that output on unit ',unitpp,' is not',
     &          ' requested before first initialisation of output file:'
               write(6,*) 
     &          '        :  See UMUI window (Initialisation of PP file',
     &          's) accessed from (Post Processing) from (Submodel ',
     &          'independent).' 
               GO TO 999
             ENDIF                                ! File stream opened? 
            STRING=MODEL_FT_UNIT(UNITPP)
            DO JJ=80,1,-1
              IF (STRING(JJ:JJ).EQ.'/') GOTO 411
            ENDDO
            ICODE=1
            CMESSAGE='STWORK  : Illegal output PPfile name'
            GOTO 999
 411        CONTINUE
            IF (JJ.GT.66) THEN
              PPNAME=STRING(JJ+1:80)
            ELSE
              PPNAME=STRING(JJ+1:JJ+14)
            ENDIF
            LEN_PPNAME=LEN(PPNAME)
             CALL FILE_OPEN(UNITPP,PPNAME,LEN_PPNAME,1,1,ICODE)
            IF(ICODE.NE.0)GOTO990
          ENDIF
CL
CL 4.1.1 Read in the pp fixed-length header
CL
          IWA=0
          CALL SETPOS(UNITPP,IWA,ICODE)
          CALL BUFFIN(UNITPP,PP_FIXHD(1),LEN_FIXHD,LEN_IO,A_IO)
          IF(A_IO.NE.-1.0.OR.LEN_IO.NE.LEN_FIXHD) THEN
            CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,
     &                LEN_FIXHD)
            CMESSAGE='STWORK  : I/O error - PP fixed length header'
            ICODE=1
            RETURN
          ENDIF
CL
CL 4.1.2 Find the first available pp lookup record.
CL
          ICURRLL=FT_LASTFIELD(UNITPP) ! Position of the last field
          ICURRLL=ICURRLL+1            ! Position of the next field
CL
CL 4.1.3 Find the first available position for the next data record(s)
CL       by reading last pp lookup record.
CL
          IWL= PP_FIXHD(150)-1  ! NOTE for BUFFIN I/O the start address
C                               ! is zero for word 1. This is pointer
C                               ! to start of lookups.
          IF(ICURRLL.EQ.1) THEN      ! First record
            IWA=PP_FIXHD(160)-1   ! Pointer to start of data

          ELSE

C  Point to start of last pp lookup record and read
            CALL SETPOS(UNITPP,IWL+(ICURRLL-2)*LEN1_LOOKUP,ICODE)
            CALL BUFFIN (UNITPP,IPPLOOK(1,ICURRLL-1),
     &                   LEN1_LOOKUP,LEN_IO,A_IO)
C
            IF(A_IO.NE.-1.0.OR.LEN_IO.NE.PP_FIXHD(151))
     &        THEN
              CALL IOERROR('Buffer in LOOKUP table       ',A_IO,LEN_IO,
     &                    PP_FIXHD(151))
               CMESSAGE='STWORK  : I/O error - PP LOOKUP table       '
               ICODE=2
               RETURN
            ENDIF
C  Pointer to next available data location in output file
            IWA=IPPLOOK(LBEGIN,ICURRLL-1)+IPPLOOK(LBNREC,ICURRLL-1)

          ENDIF                     ! Test on first record
CL
CL 4.1.4 If a daughter record is being processed then recover
CL         size information from dump LOOKUP header referenced by
CL         mother record (unless MOS output)
CL
          IF (il.NE.ilcurr .AND. .NOT.MOS) THEN
            extraw_hdr=lookup(lbext,icurrll_dump_ptr)
            global_pphoriz_out=lookup(lblrec,icurrll_dump_ptr)
            global_n_rows_out=lookup(lbrow,icurrll_dump_ptr)
            global_n_cols_out=lookup(lbnpt,icurrll_dump_ptr)
            IF (what_proc.eq.st_time_series_mean) then
! As work is done on only PE 0 and copy to buf 3 uses pphoriz_out 
! this must be reset.
              pphoriz_out=global_pphoriz_out
              n_rows_out=global_n_rows_out
              n_cols_out=global_n_cols_out
            endif
            IF (what_proc .EQ. st_time_series_code) THEN
              pphoriz_out=global_pphoriz_out
              n_rows_out=global_n_rows_out
              n_cols_out=global_n_cols_out
            ENDIF
          ENDIF
CL
CL 4.1.5 Check PP_PACK_CODE for GRIB output. Set GRIB flag and reset
CL       PP_PACK_CODE to give packing profile.
       IF(PP_PACK_CODE(UNITPP).GE.100)THEN
         GRIB_OUT=.TRUE.
         PP_PACK_CODE(UNITPP)=PP_PACK_CODE(UNITPP)-100
         GRIB_PACKING=PP_PACK_CODE(UNITPP)
       ELSE
         GRIB_OUT=.FALSE.
       ENDIF
CL
CL 4.1.6 Set packing accuracy for output data field and buffer length
CL       Multiple packing profiles are held in PP_XREF and chosen on a
CL       per-unit basis through PP_PACK_CODE.  Profile 0 means unpacked.
CL       If the field has any extra data switch off packing.
CL
          IF (PP_PACK_CODE(UNITPP).EQ.0.OR.extraw_hdr.NE.0) THEN
            PACKING=.FALSE.
            COMP_ACCRCY=-99
          ELSE
            PACKING=.TRUE.
      comp_accrcy= EXPPXI( im_ident, IS, IM, 
     &                     ppx_packing_acc+PP_PACK_CODE(UNITPP)-1,
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
     &             icode, cmessage)
          ENDIF
          IF(GRIB_OUT) THEN  ! reset packing code
            PP_PACK_CODE(UNITPP)=PP_PACK_CODE(UNITPP)+100
          ENDIF
C
          LENBUF=((global_PPHORIZ_OUT+um_sector_size-1)/um_sector_size)*
     2     um_sector_size
!        ! Output length before pack
CL
CL 4.2 Select routine to output data using logical GRIB.
CL     If data to be output in grib code then call GRIB_FILE
CL     If data to be output in PP   code then call PP_FILE
CL
        DO  II=1,NUM_LEVS_OUT           ! --- Start of levels loop ---

! Gather together distributed field to PE 0
! Distributed data is in PP_FIELD, gathered data will be
! in the buf array

          IF ( MOS .OR. (what_proc.eq.st_time_series_code) .OR.
     &    (what_proc.eq.st_time_series_mean)  ) THEN   
! if it's either MOS or timeseries output - just copy on PE 0

            IF (mype .EQ. 0) THEN
              DO I=1,PPHORIZ_OUT
                buf(I)=PPFIELD(I+(II-1)*PPHORIZ_OUT)
              ENDDO
            ENDIF

          ELSE ! not MOS or timeseries output

            CALL STASH_GATHER_FIELD (
     &        PPFIELD(1+(II-1)*PPHORIZ_OUT) , buf,
     &        PPHORIZ_OUT , global_PPHORIZ_OUT , 1,
     &        stlist(st_north_code,il) ,
     &        stlist(st_west_code,il)+global_N_COLS_OUT-1,
     &        stlist(st_north_code,il)+global_N_ROWS_OUT-1,
     &        stlist(st_west_code,il),
     &        grid_type_code,0,.TRUE.,
     &        ICODE,CMESSAGE)

            IF (ICODE .GT. 0) THEN
              WRITE(6,*) 'Error occured in STASH while gathering ',
     &                   'data for output.'
              GOTO 999
            ENDIF

          ENDIF ! IF (MOS)

!  Reset index for PP_HEAD if there is a levels list of hybrid levels
          IF (stlist(st_output_bottom,il).lt.0.and.lbvcl.eq.9) THEN
            JJ = level_list(II)
!  or a range of model levels, as they may not be consecutive,
          ELSE IF (stlist(st_output_bottom,il).ge.1.and.
     &             stlist(st_output_top,il).le.T_levels) THEN
            JJ = level_list(II)
!  otherwise use of level index is alright
          ELSE
            JJ = II
          END IF
C
C    Check that PP output file has sufficient headers pre-allocated
       IF(ICURRLL.GT.PP_FIXHD(152)) THEN
          ICODE=4
          WRITE(6,*) 'ERROR detected in routine STWORK: stop model'
          WRITE(6,*) ': No. of output fields (=',ICURRLL,')',
     &          ' exceeds no. of reserved PP headers for unit ',UNITPP
          CMESSAGE='STWORK   : NO. OF FIELDS EXCEEDS RESERVED HEADERS'
          GOTO 999
       ENDIF     ! end  no. of pp fields check
C
       IF(GRIB_OUT) THEN
!
! NOTE cannot pack data into grib before pphead correctly setup
!
         NUM_WORDS = -99   ! ie unset before call to pp_head
         PACKING_TYPE=3
       ELSE
CL Pack data into PP code.
          IF (mype .EQ. 0) THEN

          CALL PP_FILE(buf,
     1      LENBUF,NUM_WORDS,RMDI,COMP_ACCRCY,
     2      global_PPHORIZ_OUT,UNITPP,IWA,
     &      global_N_COLS_OUT,global_N_ROWS_OUT,
     3     PACKING,PACKING_TYPE,ICODE,CMESSAGE)
          ENDIF ! (IF mype.eq.0)
! Make sure all processors get the return code

          CALL GC_IBCAST(101,1,0,nproc,info,ICODE)

C Num_words is the no of 64 bit words required
          IF(ICODE.gt.0) THEN
            CMESSAGE='STWORK  : Error in PP_FILE'
            GOTO 999
          ENDIF
       ENDIF
          LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*
     2     um_sector_size ! No of words output
CL
CL 4.2.1 Set STASH processing codes and sampling period for PPheader
CL
          GR=STLIST(st_gridpoint_code,ILPREV)! Grid point code
C Any time-processed field has a (non-zero) sample_prd set -
C this will be translated by PP_HEAD into an LBTIM subcode IB of 2
          sample_prd=0.0
          IF(STLIST(st_proc_no_code,ILPREV).GT.st_replace_code) THEN
            sample_prd=REAL(STLIST(st_freq_code,ILPREV)*secs_per_period)
     &                /REAL(steps_per_period*3600)
          ENDIF
CL
CL 4.2.2 Verification time comes from fixhd(28), current time fixhd(21)
CL       2 cases that require consideration here:
CL
CL      (1) this record is not a daughter record.
CL          in which case, set start_step=.true., verif time from fixhd
CL          present time also from fixhd
CL
CL      (2) this record IS a daughter record.
CL          in which case, will need to retreive info on start_time
CL          from dump
CL
          start_step=.true.
          IF (il.eq.ilcurr) THEN ! not daughter record
            CALL PP_HEAD(
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
     *        im_ident,FIXHD,INTHD,REALHD,
     1        LEN_FIXHD,LEN_INTHD,LEN_REALHD,
     2        IM,IS,GR,lfullfield,
     3        level(II),pseudo_level(II),
     4        samples,start_step,fixhd(28),fixhd(21),LEN1_LOOKUP,
     5        extraw_hdr,IPPLOOK(1,ICURRLL),IPPLOOK(1,ICURRLL),
     6        global_N_COLS_OUT,NUM_WORDS,LEN_BUF_WORDS,
     7        global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,
     7        lbproc_comp,sample_prd,
     8        FCST_PRD,COMP_ACCRCY,PACKING_TYPE,
     &        st_grid,IWA,stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),
     &        T_levels,JJ,ROTATE,ELF,
     A        OCEAN,LEVDEPC,LEN1_LEVDEPC,
     B        ICODE,CMESSAGE)
          ELSE ! daughter record so start time is in dump
C set up start_time from data in LOOKUP(lbyr,icurrll_dump_ptr)
            start_time(1)=LOOKUP(lbyr,icurrll_dump_ptr)
            start_time(2)=LOOKUP(lbmon,icurrll_dump_ptr)
            start_time(3)=LOOKUP(lbdat,icurrll_dump_ptr)
            start_time(4)=LOOKUP(lbhr,icurrll_dump_ptr)
            start_time(5)=LOOKUP(lbmin,icurrll_dump_ptr)
            start_time(7)=LOOKUP(lbday,icurrll_dump_ptr)
            CALL PP_HEAD(
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
     *        im_ident,FIXHD,INTHD,REALHD,
     1        LEN_FIXHD,LEN_INTHD,LEN_REALHD,
     2        IM,IS,GR,lfullfield,
     3        level(II),pseudo_level(II),
     4        samples,start_step,start_time,
     4        fixhd(28),LEN1_LOOKUP,
     5        extraw_hdr,IPPLOOK(1,ICURRLL),IPPLOOK(1,ICURRLL),
     6        global_N_COLS_OUT,NUM_WORDS,LEN_BUF_WORDS,
     7        global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,
     7        lbproc_comp,sample_prd,
     8        FCST_PRD,COMP_ACCRCY,PACKING_TYPE,
     &        st_grid,IWA,stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),T_levels,JJ,
     A        ROTATE,ELF,OCEAN,LEVDEPC,LEN1_LEVDEPC,
     B        ICODE,CMESSAGE)
           ENDIF ! end if block over daughter/mother record

           IF(ICODE.gt.0) GOTO 999 ! An error has occured

           IF(GRIB_OUT) THEN
!  Now safe to call grib coder as pphead correctly set apart from
!  length of data
CL Pack data into grib code
             IF (mype .EQ. 0) THEN
               CALL GRIB_FILE(LEN1_LOOKUP,PP_LEN2_LOOKUP,
     &                  IPPLOOK,IPPLOOK,ICURRLL,
     &                  buf,global_PPHORIZ_OUT,
     &                  LENBUF,NUM_WORDS,UNITPP,IWA,GRIB_PACKING,
     &                  ICODE,CMESSAGE)
             ENDIF  ! (IF mype.eq.0)

! Make sure all processors get the return code

          CALL GC_IBCAST(101,1,0,nproc,info,ICODE)
             IF(ICODE.gt.0)THEN
               CMESSAGE='STWORK  : Error in GRIB_FILE'
               GOTO 990
             ENDIF

           ENDIF       ! end of grib_out

           IWA=IPPLOOK(LBEGIN,ICURRLL)+IPPLOOK(LBNREC,ICURRLL)
           ICURRLL=ICURRLL+1  ! Update the counter
           icurrll_dump_ptr=icurrll_dump_ptr+1
C strictly only needs doing if a daughter record
C
         ENDDO                           ! --- End of levels loop ---
         FT_LASTFIELD(UNITPP)=ICURRLL-1  ! Position of the last field
CL
CL 4.3 Write out the pp lookup table IPPLOOK (RPPLOOK) and close file
CL     if reinitialisation possible.
CL     Otherwise force last buffer to be written to file to avoid
CL     problems with continuation runs following hard failures.
CL
CL      Write out only pp lookups for records written this call.
CL      Point to the start of the current lookup record.
          IWA = IWL + (ICURRLL-NUM_LEVS_OUT-1)*LEN1_LOOKUP
          CALL SETPOS(UNITPP,IWA,ICODE)
          CALL BUFFOUT (UNITPP,IPPLOOK(1,ICURRLL-NUM_LEVS_OUT),
     *             LEN1_LOOKUP*NUM_LEVS_OUT,LEN_IO,A_IO)
C
          IF (A_IO.NE.-1.0.OR.LEN_IO.NE.PP_FIXHD(151)*NUM_LEVS_OUT)
     &      THEN
            CALL IOERROR('Buffer out LOOKUP table      ',A_IO,LEN_IO,
     &                    PP_FIXHD(151))
            CMESSAGE='STWORK  : I/O error - PP LOOKUP table       '
            ICODE=3
            RETURN
          ENDIF
            IF (FT_STEPS(UNITPP).NE.0) THEN
              LEN_PPNAME=LEN(PPNAME)
              CALL FILE_CLOSE(UNITPP,PPNAME,LEN_PPNAME,1,0,ICODE)
            ELSE
              CALL FLUSH_BUFFER(UNITPP,ICODE)
              IF(ICODE.NE.0) THEN
                CMESSAGE='STWORK: Problem flushing buffer'
                ICODE=10
                RETURN
              ENDIF

            ENDIF
C
        ELSEIF (output_code.EQ.st_dump.OR.output_code.EQ.st_secondary)
     *  THEN
CL
CL 4.4 OUTPUT to dump or secondary D1 space - this implies some
CL     TEMPORAL processing possibly.  If destination is secondary D1
CL     space, there will be no associated LOOKUP header to update.
CL
C Length is calculated from STASHlist
C NB: Full field length must be coded here, even for partial timeseries
C
          NUM_WORDS=STLIST(st_dump_output_length,IL)/NUM_LEVS_OUT
          ICURRLL=STLIST(st_lookup_ptr,IL) ! Location of DUMP header
C
          DO  II=1,NUM_LEVS_OUT          ! --- Start of levels loop ---
!  Reset index for PP_HEAD if there is a levels list of hybrid levels
            IF (stlist(st_output_bottom,il).lt.0.and.lbvcl.eq.9) THEN
              JJ = level_list(II)
!  or a range of model levels, as they may not be consecutive,
            ELSE IF (stlist(st_output_bottom,il).ge.1.and.
     &               stlist(st_output_top,il).le.T_levels) THEN
              JJ = level_list(II)
!  otherwise use of level index is alright
            ELSE
              JJ = II
            END IF
            addr=stlist(st_output_addr,il) ! start address
            IF (what_proc.eq.st_time_series_code.or.      
     &          what_proc.eq.st_time_series_mean) THEN  
                                                                        
CL
CL 4.4.1 Timeseries addresses are incremented according to timestep
CL
              IF (stlist(st_freq_code,il).lt.1) THEN
                icode=st_not_supported
                cmessage=
     +              'STWORK  : STASHtime for timeseries not supported'
                goto 999 ! got an error so jump to return
              ENDIF
              elap_time=step-stlist(st_start_time_code,il)
              elap_time=(mod(elap_time,stlist(st_period_code,il)))/
     +          stlist(st_freq_code,il)
              addr=addr+(elap_time*pphoriz_out)
CL on the first time step of a timeseries processing
CL pphoriz_out is the length of the entire output vector
CL including extra data -- on other timesteps it is
CL the length of a single record (data for just one timestep)
            ENDIF
CL
CL 4.4.2 TEMPORAL processing from ppfield array to D1
CL
            CALL TEMPORAL(ppfield(1+(ii-1)*pphoriz_out),
     +        d1(addr+(ii-1)*pphoriz_out),pphoriz_out,extraw,
     +        stlist(1,il),len_stlist,OCEAN,step,
     +        icode,cmessage,start_step,rmdi)

            IF (icode.gt.0) goto 999
CL
CL 4.4.3 Set up LOOKUP header if destination is main part of D1
CL
            IF (output_code.EQ.st_dump) THEN
CL
CL 4.4.3 Set other information for input to PPHEAD
CL
              GR=STLIST(st_gridpoint_code,ILPREV)! Grid point code
C Any time-processed field has a (non-zero) sample_prd set -
C this will be translated by PP_HEAD into an LBTIM subcode IB of 2
              sample_prd=0.0
              IF (STLIST(st_proc_no_code,ILPREV).GT.st_replace_code)
     *        THEN
                sample_prd=REAL(STLIST(st_freq_code,ILPREV)*
     &                     secs_per_period)/REAL(steps_per_period*3600)
              ENDIF
C Address of whole field is calculated from STASHlist
              IWA=STLIST(st_dump_output_addr,IL)+(II-1)*NUM_WORDS
CL
CL 4.4.4 Call PPHEAD to set LOOKUP header for field STASHed to D1.
CL       Here pass previous_time as well as start_step from temporal
CL       if start_step is true then start time will be updated.
CL       Value of end time is unimportant as that is handled properly
CL       when data is written out to pp file.
CL       Note that LBNREC is hardwired to 0 and so too is BACC.
CL
            CALL PP_HEAD(
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
     *        im_ident,FIXHD,INTHD,REALHD,
     1          LEN_FIXHD,LEN_INTHD,LEN_REALHD,
     2          IM,IS,GR,lfullfield,
     3        level(II),pseudo_level(II),
     4          samples,start_step,previous_time,fixhd(28),LEN1_LOOKUP,
     5          extraw_hdr,LOOKUP(1,ICURRLL),RLOOKUP(1,ICURRLL),
     6          global_N_COLS_OUT,NUM_WORDS,0,
     7          global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,
     7          lbproc_comp,sample_prd,
     8          FCST_PRD,0,PACKING_TYPE,
     &        st_grid,IWA,stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),T_levels,JJ,
     A          ROTATE,ELF,OCEAN,LEVDEPC,LEN1_LEVDEPC,
     B          ICODE,CMESSAGE)
C
              IF(ICODE.gt.0) goto 999 ! An error has occured
 
C Only (optionally) pack fields if no extra words of data
              IF (extraw_hdr .eq. 0) THEN
                LOOKUP(LBPACK,ICURRLL) =
     &            EXPPXI( im_ident, IS, IM, ppx_dump_packing,
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
     &             icode, cmessage)
                IF (DUMP_PACK.eq.3 ) THEN
!                 Override packing indicator from PPXREF
                  N1 = 0   !   No packing
                  LOOKUP(LBPACK,ICURRLL) =
     &            (LOOKUP(LBPACK,ICURRLL)/10)*10 + N1
                ENDIF
              ELSE
                LOOKUP(LBPACK,ICURRLL)=0
              ENDIF
C Set data type (REAL/INTEGER) from PP_XREF (-ve for timeseries)
              IF (STLIST(st_series_ptr,ilprev).GT.0.OR.
     &   STLIST(st_proc_no_code,ilprev).eq.st_time_series_mean) THEN
                LOOKUP(DATA_TYPE,ICURRLL) =
     &           -EXPPXI( im_ident, IS, IM, ppx_data_type,
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
     &             icode, cmessage)
              ELSE
                LOOKUP(DATA_TYPE,ICURRLL) =
     &            EXPPXI( im_ident, IS, IM, ppx_data_type,
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
     &             icode, cmessage)
              ENDIF
              ICURRLL=ICURRLL+1 ! Update the counter for the next field
C
            ENDIF               ! End of IF output_code=dump
C
          ENDDO                           ! --- End of levels loop ---
C
        ELSE
          ICODE=9
          CMESSAGE='STWORK  : Illegal output destination in STLIST'
          GOTO 999
        ENDIF    ! End of STLIST output destination IF block
C
      ENDIF      ! END OF S_F IF Block ---------------------------------
CL
CL 5. End of loop over STASHlist entries - RETURN to calling routine
CL
  200 CONTINUE
C
  999 CONTINUE
      RETURN
CL
CL 9. IO error exits
CL
  990 WRITE(CMESSAGE,'("STWORK  : Error opening output PP file on unit "
     &                ,I2)') UNITPP
      RETURN
      END
