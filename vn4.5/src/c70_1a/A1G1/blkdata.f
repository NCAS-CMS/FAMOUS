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
CLL Block Data Subprogram : BLKDATA
CLL
CLL Purpose : Holds DATA for any variables that are in common blocks,
CLL           so that they are initialised in only one place.
CLL
CLL Written by P.Burton
CLL
CLL Model vn.  Date    Modification history from vn3.3
CLL  3.4       1/8/94   Add C_VARCTL,C_VARCDT     Stuart Bell
CLL  4.4      24/10/97  Add CNTL_IO and CNTLIODT C.P. Jones
CLL  4.5      14/07/98  Add FLDOP to *IF DEF at top (A Van der Wal)
CLL

      BLOCK DATA BLKDATA

C*L --------------------- Comdeck: CENVIR   ----------------------------
C
C    Purpose: COMDECK defining Character enviroment variables used
C             by portable IO to open and close files
C
C    Author : R A Stratton      Date : 22/10/92
C
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
CLL
CLL 3.1     15/01/93  Increase no. of unit nos. from 1-99  to 1-199
CLL                   Dummy names have been set up temporarily for
CLL                   files 104-119. R.Rawlins
CLL
CLL 3.3     09/03/94  Separate data statements into COMDECK
CLL                   CENVIRDT. Also includes mods originally
CLL                   in RB221193 : Add source terms at unit no.110
CLL                   P.Burton and R.T.H Barnes
CLL

C    Vn3.0  12/02/93 - Environment variables PERTURB and TRANSP put in
C                      positions 37 and 97 respectively in character
C                      array FT_ENVIRON, and the appropriate character
C                      lengths put in LEN_FT_ENVIR. C. S. Douglas
C
C  Type declarations
C
      CHARACTER*8 FT_ENVIRON(199)  ! Array holding enviroment variables
C                                   for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
C


C
C Common Blocks for character and integer arrays
C
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
C
!
! Description:
!   Control for VAR UM Processing
!
! Current Code Owner: Stuart Bell
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.4   1/8/94  Original code. Stuart Bell
!   4.0   23/5/95 Add ObTypes.   Stuart Bell
!   4.0    6/6/95 Add namelist control for interpolation method. S. Bell
!   4.1   4/01/96 Increase NumModelVars and add ItemLevels. S. Bell
!   4.1  11/12/97 Increase size of ObTypes. S. Bell
!   4.5  24/2/98  Changes to Allow OBS_FORMAT=3 (directories). S. Bell
! Declarations:

! Global parameters:
      INTEGER ObsUnitNum              !Unit Num of first Obs file
      PARAMETER ( ObsUnitNum = 70 )
      INTEGER CxUnitNum               !Unit Num of first Cx file
      PARAMETER ( CxUnitNum = 120 )
      INTEGER  MaxObTypes             !Max Num. VAR Obs Types
      PARAMETER ( MaxObTypes = 10 )
      INTEGER  MaxVarFiles            !Max Num. VAR Obs Files
      PARAMETER ( MaxVarFiles = 10 )
      INTEGER NumModelVars    ! Num. of extra sect 0 items used by var
      PARAMETER ( NumModelVars = 12 )
      CHARACTER*256 NameUsedFile(MaxObTypes*MaxVarFiles)!VarObs path
      INTEGER  LenUsedFile(MaxObTypes*MaxVarFiles)!VarObs pathlength

! Global scalars:
        LOGICAL Cx              !switch: generate Cx files
        LOGICAL DiagVar         !switch: calc/print VAR diagNumstics
        INTEGER NumVarFiles     !Num of VAROBS files
        INTEGER ModeHorizInterp !Mode for horizontal Interpolation
        INTEGER OBS_FORMAT     ! Format of OBS Input 2=files,3=directory
        INTEGER NumUsedFiles

! Global dynamic arrays:
      CHARACTER*16 ObTypes(MaxObTypes) !ObType names
      INTEGER SectionIn(NumModelVars) !Input Section Number
      INTEGER ItemIn(NumModelVars)    !Input Item Number
      INTEGER ItemOut(NumModelVars)   !Output Section Number
      INTEGER ItemLevels(NumModelVars) !Number of levels of data to copy
      REAL    XOffset(NumModelVars)   !X offset from A grid
      REAL    YOffset(NumModelVars)   !Y offset from A grid

! COMMON blocks:
      COMMON /C_VARCTL/ Cx, DiagVar, NumVarFiles,
     &          SectionIn, ItemIn, ItemOut, ItemLevels, 
     &          XOffset, YOffset, ModeHorizInterp,
     &          ObTypes
     &          ,OBS_FORMAT,NumUsedFiles,NameUsedFile,LenUsedFile

!- End of COMDECK C_VARCTL
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

C-----------------------------------------------------------------------
CLCOMDECK ACPARM
CL--------------
CL -  PARAMETERS USED FOR DIMENSIONING SMALL PERMANENT ARRAYS
CL -  VALUES ARE MAXIMUM LIKELY TO AVOID RECOMPILE ON RESOLUTION CHANGE
CL -  THE ACTUAL DIMENSIONS ARE PASSED AS ARGUMENTS
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton

      INTEGER

     &  ROW_LENGTH_MAX  ! Maximum row length
     &, P_ROWS_MAX      ! Maximum number of rows
     &, HORIZ_DIM_MAX   ! MAX(ROW_LENGTH_MAX,P_ROWS_MAX)
     &, HALO_MAX        ! Maximum MPP halo width
     &, P_LEVELS_MAX    ! Maximum number of total levels
     &, Q_LEVELS_MAX    ! Maximum number of wet levels

      PARAMETER ( ROW_LENGTH_MAX = 432
     &,           P_ROWS_MAX = 325
     &,           HORIZ_DIM_MAX = 432
     &,           HALO_MAX = 2  ! fourth order double width halo
     &,           P_LEVELS_MAX = 42
     &,           Q_LEVELS_MAX = 42)

! Derived sizes

      INTEGER
     &  Max2DFieldSize
     &, Max3DFieldSize
     &, MaxHaloSize

      PARAMETER (
     &  Max2DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX
     &, Max3DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX*P_LEVELS_MAX
     &, MaxHaloSize = HORIZ_DIM_MAX*HALO_MAX
     & )

      INTEGER
     &  MAX_SHMEM_COMMON_SIZE,
     &  MaxFieldSize,
     &  MaxFieldSizeMes                                                 
      PARAMETER ( MAX_SHMEM_COMMON_SIZE = 3000000 ,
     &            MaxFieldSize   = MAX_SHMEM_COMMON_SIZE/4 ,
     &            MaxFieldSizeMes= MAX_SHMEM_COMMON_SIZE/6 )
      INTEGER NOBTYPMX
      PARAMETER ( NOBTYPMX       =  61 )
C-----------------------------------------------------------------------
CLCOMDECK COMOBS
! Description                                                  
!   This comdeck provides parameters relating to observations in 
!   atmospheric AC assimilation (see DOCOBS for detail)
!                                                                   
!   History:                                                         
!   Model    Date     Modification history                       
!  version                                                              
!   4.4      3/7/97   Add MAX_NUM_ACOB_FILES,PER_FILE_TNDVMAX S. Bell
!   4.5      5/8/98   Increase USED_FILES size S. Bell
CL-------------------------------------------------------------------   
      INTEGER NDATAVMX, NOBLEVMX
      INTEGER MAX_NUM_ACOB_FILES
      INTEGER NUM_OB_FILE_TYPES
      PARAMETER (NDATAVMX = 6+3*P_LEVELS_MAX)
      PARAMETER (NOBLEVMX = P_LEVELS_MAX+1)
      PARAMETER (MAX_NUM_ACOB_FILES=100)
      PARAMETER (NUM_OB_FILE_TYPES = 10)
      INTEGER NOBTYP,NDVHDR,MAXNLEV1,
     + OBSTYP(NOBTYPMX),NOBLEV(NOBTYPMX),NDATAV(NOBTYPMX),
     + NERLEV1(NOBTYPMX),NOBS(NOBTYPMX),OBLEVTYP(NOBTYPMX),
     + MDISPOBT(NOBTYPMX),OBS_NO_ST(NOBTYPMX),
     + OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN

      REAL           MISSD,
     +               OBLEVELS(NOBLEVMX,NOBTYPMX),
     +               OBLAYERB(NOBLEVMX+1,NOBTYPMX),
     +               TIMEINT,TIMENEXT,
     +               OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W

      INTEGER PER_FILE_TNDVMAX(MAX_NUM_ACOB_FILES)

      CHARACTER*30 OB_FILE_TYPE(NUM_OB_FILE_TYPES)

      CHARACTER*256 USED_FILES(MAX_NUM_ACOB_FILES)
      INTEGER FILENAME_LEN(MAX_NUM_ACOB_FILES)
      INTEGER NUM_USED_FILES
                                                                        
      COMMON /COMOBS/ NOBTYP,NDVHDR,MAXNLEV1,
     + OBSTYP,NOBLEV,NDATAV,NERLEV1,NOBS,OBLEVTYP,
     + MDISPOBT,OBS_NO_ST,MISSD,OBLEVELS,OBLAYERB,
     + OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN,
     + TIMEINT, TIMENEXT,
     + OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W, PER_FILE_TNDVMAX,
     + OB_FILE_TYPE,USED_FILES,FILENAME_LEN,NUM_USED_FILES
C-----------------------------------------------------------------------

! History:
! Version  Date  Comment
!  3.4   18/5/94 Correct misspelling of CURNTIN and change length
!                to match. J F Thomson
C*L --------------------- Comdeck: CENVIRDT ---------------------------
C
C    Purpose: Data statements for character enviroment variables used
C             by portable IO to open/close files (links with CENVIR)
C
C    Author : R A Stratton      Date : 22/10/92
C
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  3.4   1/8/94     Revised Obs file spec + CX files: Stuart Bell
!    4.0  22/9/95     Units for spectral data added. J. M. Edwards
CLL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1 26/02/96     New env. variables for sulphur ancillary files.
!                     Rename SOURCES to SULPEMIS. D. Robinson.
!      4.3  18/3/97   Add aerosol fcgs of climate change. William Ingram
!    4.4  4/7/97    Add ANLINCR at 108. Chris Jones/Stuart Bell
CLL  4.4 12/9/97      New ancillary file environment variables for 
CLL                   initial surface type fracs, initial vegetation 
CLL                   state and vegetation disturbance.
CLL                                                  R.A.Betts
!    4.4 28/08/97     Move CACHED from logical unit no 3 to 138 in
!                     order to release a Fortran unit [for OASIS].
!                     R.Rawlins
!    4.5 22/04/98     Add new ancillary file for CO2 emissions:
!                     CO2EMITS - in I/O unit 118. Chris Jones.
!    4.5 22/04/98     Add new ancillary file for soot emissions:
!                     SOOTEMIS - in I/O unit 139. R.Rawlins             
!    4.5 29/07/98     Move ALABCOU1/2/3/4 from 101-103 to 140-143.
!                     Add ALABCOU5/6/7/8 to 144-147. D. Robinson. 
!    4.5 17/08/98     Remove OLABCOUT from Unit 90. Add OLABCOU1/2/3/4
!                     to 101-103. D. Robinson.
CLL
CLL DATA statements for COMDECK CENVIR

      DATA FT_ENVIRON/
     &  'PPXREF  ','PPXREFU ','        ','STASHCTL','        ', !  1- 5
     &  '        ','OUTPUT2 ','        ','        ','XHIST   ',
     &  'IHIST   ','THIST   ','        ','ERRFLAG ','CACHE1  ',
     &  'CACHE2  ','AOTRANS ','ASWAP   ','OSWAP   ','AINITIAL',
     &  'ASTART  ','        ','APSUM1  ','APSTMP1 ','        ',
     &  '        ','AOMEAN  ','ATMANL  ','        ','OZONE   ',
     &  'SMCSNOWD','DSOILTMP','SOILTYPE','VEGTYPE ','SSTIN   ',
     &  'SICEIN  ','PERTURB ','CURNTIN ','        ','OINITIAL',
     &  'OSTART  ','        ','OPSUM1  ','OPSTMP1 ','        ',
     &  '        ','OCNANL  ','ATRACER ','OTRACER ','WFIN    ',
     &  'HFLUXIN ','PMEIN   ','ICEFIN  ','AIRTMP  ','SALINITY',
     &  'FLUXCORR','SWSPECTD','BAS_IND ','SLABHCON','PP0     ',         
     &  'PP1     ','PP2     ','PP3     ','PP4     ','PP5     ',
     &  'PP6     ','PP7     ','PP8     ','PP9     ','OBS01   ',
     &  'OBS02   ','OBS03   ','OBS04   ','OBS05   ','OBS06   ',
     &  'OBS07   ','OBS08   ','OBS09   ','OBS10   ','LWSPECTD',         
     &  'WAVEOUT ','SURGEOUT','MESOUT  ','STRATOUT','WFOUT   ',
     &  'HFLUXOUT','PMEOUT  ','ICFOUT  ','MOSOUT  ','FILE90  ',
     &  'SSTOUT  ','SICEOUT ','CURNOUT ','        ','ALABCIN ',
     &  'OROG    ','TRANSP  ','OLABCIN ','OCNDEPTH',
     &  'OLABCOU1','OLABCOU2','OLABCOU3','OLABCOU4','FILE104 ',
     &  'FILE105 ','FILE106 ','FILE107 ','ANLINCR ','MURKFILE',
     &  'SULPEMIS','USRANCIL','USRMULTI','OUSRANCL','OUSRMULT',
     &  'SO2NATEM','CHEMOXID','AEROFCG ','CO2EMITS','FILE119 ',
     &  'CX01    ','CX02    ','CX03    ','CX04    ','CX05    ',
     &  'CX06    ','CX07    ','CX08    ','CX09    ','CX10    ',
     &  'WINITIAL','WSTART  ','WRESTART','WAVANL  ','WAVANCIN', 
     &  'FRACINIT','VEGINIT ','DISTURB ','CACHED  ','SOOTEMIS',!135-139
     &  'ALABCOU1','ALABCOU2','ALABCOU3','ALABCOU4','ALABCOU5',!140-144
     &  'ALABCOU6','ALABCOU7','ALABCOU8','        ','        ',!145-150
     &          50*'        '
     & /
C
      DATA LEN_FT_ENVIR/6,7,0,8,0,0,7,0,0,5,5,5,0,7,6,6,7,5,5,8,!  1- 20

     &                  6,0,6,7,0,0,6,6,0,5,8,8,8,7,5,6,7,7,0,8,
     &                  6,0,6,7,0,0,6,7,7,4,7,5,6,6,8,8,8,7,8,3,        
     &                  3,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,8,
     &                  7,8,6,8,5, 8,6,6,6,6,       ! 81-90
     &                  6,7,7,0,7, 4,6,7,8,         ! 91-99
     &                  8,8,8,8,7, 7,7,7,7,8,       ! 100-109
     &                  8,8,8,8,8, 8,8,7,8,7,       ! 110-119
     &                  4,4,4,4,4, 4,4,4,4,4,       ! 120-129
     &                  8,6,8,6,8, 8,7,7,6,8,       ! 130-139
     &                  8,8,8,8,8, 8,8,8,0,0,       ! 140-149
     &                  50*0/                       ! 150-199

C End of COMDECK CENVIRDT

C
!
! Description:
!   Data Statements for C_VARCTL
!
! Current Code Owner: Stuart Bell
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.4   1/8/94    New    Stuart Bell
!   4.0   23/5/95   Add ObTypes.  Stuart Bell
!   4.0    9/6/95   correct ItemOut.  Stuart Bell
!   4.1   4/01/96   add extra BL coeffs from section 3.  Stuart Bell
!   4.4   11/12/97  Revise ObTypes names.  Stuart Bell
!   4.5   24/02/98  Revise ObTypes Array again.  Stuart Bell

      DATA      SectionIn  / 6*0          ,  3,  3,  3,  3,  3,  3 /
      DATA      ItemIn     / 1,2,3,4,10,33,225,226,236,245,256,257 /
      DATA      ItemOut    / 1,2,3,4,10,33,281,282,283,284,285,286 /
      DATA      ItemLevels / 6*-1         ,5*1                ,  4 /
      DATA      ObTypes   / 'Surface         ','Scatwind        ',
     &                      'Satwind         ','Aircraft        ',
     &                      'Sonde           ','Sat120          ',
     &                      'GLOSS           ','Test            ',
     &                      'MOPS            ','None            '/      
!- End of COMDECK C_VARCDT
!====================== COMDECK CNTLIODT ========================
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
!   4.4    24/10/97  New deck       C.P. Jones
!   4.5    02/10/98  Increase from 512 to 2048. D. Robinson.
!
C
      DATA UM_SECTOR_SIZE/2048/
!========================== COMDECK PARVARDT ====================
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     PARCOMM.
!
!   History:
!
!   4.4    20/11/98  New comdeck.  A. Van der Wal
!
C
      DATA current_decomp_type/-1/  ! set the initial decomposition
C                                   ! to an "unset" value
C
! ---------------------- End of comdeck PARVARDT ----------------------
!========================== COMDECK DECOMPDT ====================
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     DECOMPDB.
!
!   History:
!
!   4.4    20/11/98  New comdeck.  A. Van der Wal
!
C
      DATA decomp_db_set / max_decomps * .FALSE. /
C
! ---------------------- End of comdeck DECOMPDT ----------------------
!========================== COMDECK COMOBSDT ====================
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     COMOBS.
!
!   History:
!
!   4.6    31/08/99  New comdeck.  A. Van der Wal
!
C
      DATA OB_FILE_TYPE/
     + 'Surface                       ',
     + 'Sonde                         ',
     + 'Aircraft                      ',
     + 'Sat120                        ',
     + 'Sat500                        ',
     + 'GLOSS                         ',
     + 'Satwind                       ',
     + 'Scatwind                      ',
     + 'MOPS                          ',
     + 'Test                          '/
! ---------------------- End of comdeck COMOBSDT ----------------------

      END
