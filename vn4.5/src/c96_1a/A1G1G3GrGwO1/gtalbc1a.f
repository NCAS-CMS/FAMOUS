C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Gathers in atmosphere partial LBCs from boundaries to create
!  global LBC
!
! Subroutine Interface
      SUBROUTINE GATHER_ATMOS_LBCS(FULL_LBC,FULL_LENRIMDATA,
     &                             PART_LBC,PART_LENRIMDATA,
     &                             GATHER_PE,
     &                             ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Description:
! Gathers atmosphere LBCs from relevant processors and assembles
! a global LBC on processor GATHER_PE
!
! Method:
! The code loops around over all the processors.  For each
! processor there is then a loop over the four boundaries
! (North, East, South and West) and if the processor is on
! that particular boundary region, the appropriate data
! is sent from that processor's PART_LBC array to
! GATHER_PE's FULL_LBC array.
! The data is transferred via two arrays, buf (the data to be
! sent) and buf_expand (the data being received) which are both
! on COMMON, to make the comms work with CRAY shmem correctly.
!
! The data structure of the PART_LBC array is a little different
! from the "standard" LBC data structure as described in
! documentation paper C7. The East and West boundarys are both
! dimensioned to have P_ROWS of data on all grids. Not all
! of the rows are used (rows are used starting from the top, so
! some of the lower rows may not contain meaningful data on
! some processors)
!
! Current code owner : Paul Burton
!
! History
!  Model    Date       Modification history from model version 4.1
!  version
!    4.1    3/1/96     New Deck for MPP code   P.Burton
!    4.2    18/11/96   Tidy up use of COMMON blocks.  P.Burton
!    4.2    23/10/96   Use GC_SETOPT to use PUTs under GCOM_shmem
!                      P.Burton
!    4.4    08/12/97   Extend to cope with QCF prognostic if mixed
!                      phase precip scheme in use.  R.T.H.Barnes.
!    4.5    13/01/98   Removed SHMEM COMMON block and replaced by
!                      dynamic arrays.                   P.Burton
!    4.5    26/05/98   Corrected routine name in error message
!                                                    P. Burton
!    4.5    27/08/98   Corrected indexing for tracer variable
!                                                    P.Burton
!    4.5    15/09/98   Replaced L_LSPICE by L_LSPICE_BDY to correct
!                      sizes for non-mixed phase boundary condition
!                      fields in mixed phase dump.   R.Rawlins
!
! Subroutine Arguments:

      INTEGER
     &  FULL_LENRIMDATA   ! IN size of the FULL_LBC array
     &, PART_LENRIMDATA   ! IN size of the PART_LBC array
     &, GATHER_PE         ! IN processor to scatter data from
     &, ICODE             ! OUT error code

      CHARACTER*(80)
     &  CMESSAGE          ! OUT error message

      REAL
     &  FULL_LBC(FULL_LENRIMDATA)  ! IN full LBC (only on GATHER-_PE)
     &, PART_LBC(PART_LENRIMDATA)  ! OUT local part of LBC

! Parameters and COMMON

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
CDIR$ FIXED
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GC - General Communication primitives package. For use on
C multiprocessor shared memory and message passing systems.
C
C
C LICENSING TERMS
C
C  GC is provided free of charge. Unless otherwise agreed with SINTEF,
C  use and redistribution in source and binary forms are permitted
C  provided that
C
C      (1) source distributions retain all comments appearing within
C          this file header, and
C
C      (2) distributions including binaries display the following
C          acknowledgement:
C
C              "This product includes software developed by SINTEF.",
C
C          in the documentation or other materials provided with the
C          distribution and in all advertising materials mentioning
C          features or use of this software.
C
C  The name of SINTEF may not be used to endorse or promote products
C  derived from this software without specific prior written
C  permission.  SINTEF disclaims any warranty that this software will
C  be fit for any specific purposes. In no event shall SINTEF be liable
C  for any loss of performance or for indirect or consequential damage
C  or direct or indirect injury of any kind. In no case shall SINTEF
C  be liable for any representation or warranty make to any third party
C  by the users of this software.
C
C
C Fortran header file. PLEASE use the parameter variables in user
C routines calling GC and NOT the numeric values. The latter are
C subject to change without further notice.
C
C---------------------------------------------- ------------------------
C $Id: gpb2f402,v 1.10 1996/11/28 20:36:24 t11pb Exp $
C (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

C    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
C                    value to be set from a shell variable.
C                      Author: Bob Carruthers  Cray Research.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     GC general options
      INTEGER GC_OK, GC_FAIL, GC_NONE, GC_ANY, GC_DONTCARE,
     $     GC_SHM_DIR, GC_SHM_GET, GC_SHM_PUT, GC_USE_GET, GC_USE_PUT
     &   , GC_NAM_TIMEOUT, GC_SHM_SAFE
      PARAMETER (GC_OK         =     0)
      PARAMETER (GC_FAIL       =    -1)
      PARAMETER (GC_NONE       =     0)
      PARAMETER (GC_ANY        =    -1)
      PARAMETER (GC_DONTCARE   =    -1)
      PARAMETER (GC_SHM_DIR    =     1)
      PARAMETER (GC_SHM_SAFE   =     2)
      PARAMETER (GC_NAM_TIMEOUT=     4)
      PARAMETER (GC_SHM_GET    = -9999)
      PARAMETER (GC_SHM_PUT    = -9998)
      PARAMETER (GC_USE_GET    = -9999)
      PARAMETER (GC_USE_PUT    = -9998)

C     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

C     GC groups (GCG) support
      INTEGER GC_ALLGROUP, GCG_ALL
      PARAMETER (GC_ALLGROUP = 0)
      PARAMETER (GCG_ALL = GC_ALLGROUP)

C     GC groups (GCG) functions
      INTEGER GCG_ME

C     GC reserved message tags
      INTEGER GC_MTAG_LOW, GC_MTAG_HIGH
      PARAMETER (GC_MTAG_LOW   = 999999901)
      PARAMETER (GC_MTAG_HIGH  = 999999999)

C     GCG_RALLETOALLE index parameters
      INTEGER S_DESTINATION_PE, S_BASE_ADDRESS_IN_SEND_ARRAY,
     $     S_NUMBER_OF_ELEMENTS_IN_ITEM, S_STRIDE_IN_SEND_ARRAY,
     $     S_ELEMENT_LENGTH, S_BASE_ADDRESS_IN_RECV_ARRAY,
     $     S_STRIDE_IN_RECV_ARRAY
      PARAMETER (S_DESTINATION_PE = 1)
      PARAMETER (S_BASE_ADDRESS_IN_SEND_ARRAY = 2)
      PARAMETER (S_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (S_STRIDE_IN_SEND_ARRAY = 4)
      PARAMETER (S_ELEMENT_LENGTH = 5)
      PARAMETER (S_BASE_ADDRESS_IN_RECV_ARRAY = 6)
      PARAMETER (S_STRIDE_IN_RECV_ARRAY = 7)

      INTEGER R_SOURCE_PE, R_BASE_ADDRESS_IN_RECV_ARRAY,
     $     R_NUMBER_OF_ELEMENTS_IN_ITEM, R_STRIDE_IN_RECV_ARRAY,
     $     R_ELEMENT_LENGTH, R_BASE_ADDRESS_IN_SEND_ARRAY,
     $     R_STRIDE_IN_SEND_ARRAY
      PARAMETER (R_SOURCE_PE = 1)
      PARAMETER (R_BASE_ADDRESS_IN_RECV_ARRAY = 2)
      PARAMETER (R_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (R_STRIDE_IN_RECV_ARRAY = 4)
      PARAMETER (R_ELEMENT_LENGTH = 5)
      PARAMETER (R_BASE_ADDRESS_IN_SEND_ARRAY = 6)
      PARAMETER (R_STRIDE_IN_SEND_ARRAY = 7)

      REAL
     &  local_buf(FULL_LENRIMDATA)
     &, buf_expand(FULL_LENRIMDATA)

! TYPESIZE contains information about the tracers, and also
! P_LEVELS, Q_LEVELS and RIMWIDTHA which we require
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
C Grid related sizes for COUPLING between ATMOS and OCEAN submodels
C [For MPP, sizes are 'global' values over all PEs.] 
      INTEGER
     &  AOCPL_IMT                ! Ocean rowlength
     & ,AOCPL_JMT                ! Ocean no. of rows
     & ,AOCPL_ROW_LENGTH         ! Atmos rowlength
     & ,AOCPL_P_ROWS             ! Atmos no. of p rows
      COMMON/SIZE_AOCPL/AOCPL_IMT,AOCPL_JMT,
     &                  AOCPL_ROW_LENGTH,AOCPL_P_ROWS

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
C     : SUB-MODEL OCEAN        :
C Data structure sizes derived from grid size
      INTEGER
     &       O_LEN1_LEVDEPC,      ! IN: 1st dim of level  dep const
     &       O_LEN1_ROWDEPC,      ! IN: 1st dim of row    dep const
     &       O_LEN1_COLDEPC,      ! IN: 1st dim of column dep const
     &       O_LEN1_FLDDEPC       ! IN: 1st dim of field  dep const
C Data structure sizes for OCEAN INTERFACE file control routines        
      INTEGER
     &       INTF_LEN2_LEVDEPC_O, ! 2nd dimension of level dep consts
     &       INTF_LOOKUPSO,       ! No of interface lookups (ocean)
     &       MAX_INTF_P_LEVELS_O, ! Max no of model levels in all areas
     &       TOT_LEN_INTFO_P,     ! { Total length of single level  
     &       TOT_LEN_INTFO_U,     ! { interface fields; p & u grids
     &       NPTS_U_FIELD_O       ! No of points in ocean vely field
      COMMON/DSIZE_O/
     &        O_LEN1_LEVDEPC,O_LEN1_FLDDEPC,O_LEN1_ROWDEPC,
     &        O_LEN1_COLDEPC,
     &       INTF_LEN2_LEVDEPC_O, INTF_LOOKUPSO, MAX_INTF_P_LEVELS_O, 
     &       TOT_LEN_INTFO_P, TOT_LEN_INTFO_U, NPTS_U_FIELD_O 
C     : SUB MODEL OCEAN
C  are held in TYPOCPAR

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

! CNTLATM contains L_LSPICE to say whether LBC file needs QCF
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

! Local variables

      INTEGER
     &  global_ROW_LENGTH ! length of row on global P grid
     &, global_P_ROWS     ! number of rows on global P grid
     &, local_ROW_LENGTH  ! length of row on local P grid (no halos)
     &, local_U_ROW_LENGTH ! length of row on local U grid (no halos)
     &, local_P_ROWS      ! number of rows on local P grid (no halos)
     &, global_P_LENRIM   ! size of boundary data for single
!                         !   level on P grid for global data
     &, global_U_LENRIM   ! size of boundary data for single
!                         !   level on U grid for global data
     &, local_P_LENRIM    ! size of boundary data for single
!                         !    level on P grid for local data
     &, local_U_LENRIM    ! size of boundary data for single
!                         !    level on U grid for local data

     &, global_FIRST_SIDE_ROW  ! first row for side LBC section
     &, global_LAST_P_SIDE_ROW ! last row for side LBC section
     &, global_LAST_U_SIDE_ROW ! last row for side LBC section

     &, global_START_ROW       ! first row in global data in boundary
     &, local_START_ROW        ! first row in local data in boundary
     &, global_START_POINT     ! first point in row in global boundary
     &, local_START_POINT      ! first point in row in local boundary
     &, N_P_ROWS               ! number of rows on P grid in boundary
     &, N_P_POINTS             ! number of points in row on P grid
     &, N_U_ROWS               ! number of rows on U grid in boundary
     &, N_U_POINTS             ! number of points in row on U grid
     &, global_P_ROW_LEN       ! length of row on P grid in boundary
     &, global_U_ROW_LEN       ! length of row on U grid in boundary
     &, local_P_ROW_LEN        ! length of row on P grid in boundary
     &, local_U_ROW_LEN        ! length of row on U grid in boundary
     &, global_P_bound_off     ! offset of a P grid boundary
     &, local_P_bound_off      ! offset of a P grid boundary
     &, global_U_bound_off     ! offset of a U grid boundary
     &, local_U_bound_off      ! offset of a U grid boundary


     &, global_PSTAR_START      ! Start addresses of variables
     &, global_U_START          ! in the global LBC data
     &, global_V_START          !
     &, global_THETA_START      !
     &, global_Q_START          !
     &, global_TR_START         ! (Tracer variables)
     &, global_QCF_START        ! cloud ice if mixed phase scheme
     &, global_P_EAST_DATA_off  ! offset of East P data in global LBC
     &, global_U_EAST_DATA_off  ! offset of East U data in global LBC
     &, global_P_SOUTH_DATA_off ! offset of South P data in global LBC
     &, global_U_SOUTH_DATA_off ! offset of South U data in global LBC
     &, global_P_WEST_DATA_off  ! offset of West P data in global LBC
     &, global_U_WEST_DATA_off  ! offset of West U data in global LBC

     &, local_PSTAR_START       ! Start addresses of variables
     &, local_U_START           ! in the local LBC data
     &, local_V_START           !
     &, local_THETA_START       !
     &, local_Q_START           !
     &, local_TR_START          ! (Tracer variables)
     &, local_QCF_START         ! cloud ice if mixed phase scheme
     &, local_P_EAST_DATA_off   ! offset of East P data in local LBC
     &, local_U_EAST_DATA_off   ! offset of East U data in local LBC
     &, local_P_SOUTH_DATA_off  ! offset of South P data in local LBC
     &, local_U_SOUTH_DATA_off  ! offset of South U data in local LBC
     &, local_P_WEST_DATA_off   ! offset of West P data in local LBC
     &, local_U_WEST_DATA_off   ! offset of West U data in local LBC

      INTEGER
     &  iproc  ! processor index
     &, bound_type ! loop index for loop over boundary types
     &, data_size  ! size of message to be sent
     &, buf_pt     ! pointer to position in buffer array
     &, ROW        ! loop index indicating row number
     &, POINT      ! loop index indicating point along row
     &, LEVEL      ! loop index indicating level
     &, TVAR       ! loop index indicating tracer variable
     &, info       ! return code from communications routines

      INTEGER  ! magic numbers for boundary types
     &  bt_top   ! top boundary
     &, bt_right ! right boundary
     &, bt_base  ! bottom boundary
     &, bt_left  ! left boundary
      PARAMETER(bt_top=1,bt_right=2,bt_base=3,bt_left=4)

      LOGICAL
     &  iproc_at_left  ! processor at left of LPG
     &, iproc_at_right ! processor at right of LPG
     &, iproc_at_top   ! processor at top of LPG
     &, iproc_at_base  ! processor at base of LPG

!--------------------------------------------------------------------
!--------------------------------------------------------------------

! 1.0 Set up sizes and addresses for global data

! Set up some sizes for the global data

      global_ROW_LENGTH=glsize(1)
      global_P_ROWS=glsize(2)
      global_P_LENRIM=(global_ROW_LENGTH+global_P_ROWS-2*RIMWIDTHA)*
     &              2*RIMWIDTHA
      global_U_LENRIM=global_P_LENRIM-4*RIMWIDTHA

! Set up some addresses

      global_PSTAR_START=1
      global_U_START=global_P_LENRIM+1
      global_V_START=global_U_START+global_U_LENRIM*P_LEVELS
      global_THETA_START=global_V_START+global_U_LENRIM*P_LEVELS
      global_Q_START=global_THETA_START+global_P_LENRIM*P_LEVELS
      global_TR_START=global_Q_START+global_P_LENRIM*Q_LEVELS
      global_QCF_START=
     & global_TR_START+global_P_LENRIM*TR_LEVELS*TR_VARS

      global_P_EAST_DATA_off=global_ROW_LENGTH*RIMWIDTHA
      global_U_EAST_DATA_off=(global_ROW_LENGTH-1)*RIMWIDTHA
      global_P_SOUTH_DATA_off=global_P_EAST_DATA_off+
     &                        (global_P_ROWS-2*RIMWIDTHA)*RIMWIDTHA
      global_U_SOUTH_DATA_off=global_U_EAST_DATA_off+
     &                        (global_P_ROWS-2*RIMWIDTHA-1)*RIMWIDTHA
      global_P_WEST_DATA_off=global_P_SOUTH_DATA_off+
     &                       global_ROW_LENGTH*RIMWIDTHA
      global_U_WEST_DATA_off=global_U_SOUTH_DATA_off+
     &                       (global_ROW_LENGTH-1)*RIMWIDTHA

      global_FIRST_SIDE_ROW=RIMWIDTHA+1
      global_LAST_P_SIDE_ROW=global_P_ROWS-RIMWIDTHA
      global_LAST_U_SIDE_ROW=global_LAST_P_SIDE_ROW-1

!--------------------------------------------------------------------

! 2.0 Loop over processors

      DO iproc=first_comp_pe,last_comp_pe  ! loop over all processors


!   2.1 Set up logicals, sizes and addresses for local data

! set up logicals indicating position on LPG
        iproc_at_left=.FALSE.
        iproc_at_right=.FALSE.
        iproc_at_top=.FALSE.
        iproc_at_base=.FALSE.
        IF (g_gridpos(1,iproc) .EQ. 0) iproc_at_left=.TRUE.
        IF (g_gridpos(1,iproc) .EQ. nproc_x-1) iproc_at_right=.TRUE.
        IF (g_gridpos(2,iproc) .EQ. 0) iproc_at_top=.TRUE.
        IF (g_gridpos(2,iproc) .EQ. nproc_y-1) iproc_at_base=.TRUE.

! Set up the local data for this processor
        local_ROW_LENGTH=g_blsizep(1,iproc)  ! no halos
        IF (iproc_at_right) THEN
!        ! This processor at right of LPG so one less point on
!        ! U grid
          local_U_ROW_LENGTH=local_ROW_LENGTH-1
        ELSE
          local_U_ROW_LENGTH=local_ROW_LENGTH
        ENDIF

        local_P_ROWS=g_blsizep(2,iproc)      ! again, no halos
        local_P_LENRIM=(local_ROW_LENGTH+local_P_ROWS)*
     &                  2*RIMWIDTHA
        local_U_LENRIM=(local_U_ROW_LENGTH+local_P_ROWS)*
     &                  2*RIMWIDTHA

! Set up some addresses
        local_PSTAR_START=1
        local_U_START=local_P_LENRIM+1
        local_V_START=local_U_START+local_U_LENRIM*P_LEVELS
        local_THETA_START=local_V_START+local_U_LENRIM*P_LEVELS
        local_Q_START=local_THETA_START+local_P_LENRIM*P_LEVELS
        local_TR_START=local_Q_START+local_P_LENRIM*Q_LEVELS
        local_QCF_START=
     &   local_TR_START+local_P_LENRIM*TR_LEVELS*TR_VARS

        local_P_EAST_DATA_off=local_ROW_LENGTH*RIMWIDTHA
        local_U_EAST_DATA_off=local_U_ROW_LENGTH*RIMWIDTHA
        local_P_SOUTH_DATA_off=local_P_EAST_DATA_off+
     &                         local_P_ROWS*RIMWIDTHA
        local_U_SOUTH_DATA_off=local_U_EAST_DATA_off+
     &                         local_P_ROWS*RIMWIDTHA
        local_P_WEST_DATA_off=local_P_SOUTH_DATA_off+
     &                        local_ROW_LENGTH*RIMWIDTHA
        local_U_WEST_DATA_off=local_U_SOUTH_DATA_off+
     &                        local_U_ROW_LENGTH*RIMWIDTHA


!--------------------------------------------------------------------

!   2.2 Loop over boundaries: North, East, South and West

        DO bound_type=bt_top,bt_left  ! loop over all boundaries

          IF (((bound_type .EQ. bt_top) .AND. (iproc_at_top)) .OR.
     &       ((bound_type .EQ. bt_right) .AND. (iproc_at_right)) .OR.
     &        ((bound_type .EQ. bt_base) .AND. (iproc_at_base)) .OR.
     &        ((bound_type .EQ. bt_left) .AND. (iproc_at_left))) THEN
!             Processor iproc has a boundary of type bound_type


!      2.2.1 Set up data pointers and sizes for this boundary

            IF ((bound_type .EQ. bt_top) .OR.   ! What type of
     &          (bound_type .EQ. bt_base)) THEN ! boundary is it?

!             Northern or Southern boundary

              global_START_ROW=1
              local_START_ROW=1

              global_START_POINT=g_datastart(1,iproc)
              local_START_POINT=1

              N_P_ROWS=RIMWIDTHA
              N_P_POINTS=local_ROW_LENGTH
              N_U_ROWS=RIMWIDTHA
              N_U_POINTS=local_U_ROW_LENGTH

              global_P_ROW_LEN=global_ROW_LENGTH
              global_U_ROW_LEN=global_ROW_LENGTH-1
              local_P_ROW_LEN=local_ROW_LENGTH
              local_U_ROW_LEN=local_U_ROW_LENGTH

              IF (bound_type .EQ. bt_top) THEN  ! Northern boundary
                global_P_bound_off=0
                local_P_bound_off=0
                global_U_bound_off=0
                local_U_bound_off=0
              ELSE  ! Southern boundary
                global_P_bound_off=global_P_SOUTH_DATA_off
                local_P_bound_off=local_P_SOUTH_DATA_off
                global_U_bound_off=global_U_SOUTH_DATA_off
                local_U_bound_off=local_U_SOUTH_DATA_off
              ENDIF

            ELSE  ! Eastern or Western boundary

              global_START_ROW=MAX(g_datastart(2,iproc),
     &                             global_FIRST_SIDE_ROW)-
     &                         RIMWIDTHA
              local_START_ROW=1

              global_START_POINT=1
              local_START_POINT=1

              N_P_ROWS=MIN(g_datastart(2,iproc)+local_P_ROWS-1,
     &                     global_LAST_P_SIDE_ROW) -
     &                 (global_START_ROW+RIMWIDTHA) + 1
              N_P_POINTS=RIMWIDTHA
              N_U_ROWS=MIN(g_datastart(2,iproc)+local_P_ROWS-1,
     &                     global_LAST_U_SIDE_ROW) -
     &                 (global_START_ROW+RIMWIDTHA) + 1
              N_U_POINTS=RIMWIDTHA

              global_P_ROW_LEN=RIMWIDTHA
              global_U_ROW_LEN=RIMWIDTHA
              local_P_ROW_LEN=RIMWIDTHA
              local_U_ROW_LEN=RIMWIDTHA

              IF (bound_type .EQ. bt_right) THEN ! Eastern boundary
                global_P_bound_off=global_P_EAST_DATA_off
                local_P_bound_off=local_P_EAST_DATA_off
                global_U_bound_off=global_U_EAST_DATA_off
                local_U_bound_off=local_U_EAST_DATA_off
              ELSE ! Western boundary
                global_P_bound_off=global_P_WEST_DATA_off
                local_P_bound_off=local_P_WEST_DATA_off
                global_U_bound_off=global_U_WEST_DATA_off
                local_U_bound_off=local_U_WEST_DATA_off
              ENDIF

            ENDIF ! What type of boundary is it?

            data_size=(1+P_LEVELS+Q_LEVELS+(TR_VARS*TR_LEVELS))*
     &                   N_P_ROWS*N_P_POINTS +
     &                 2*P_LEVELS*N_U_ROWS*N_U_POINTS
            if (L_LSPICE_BDY) then ! Mixed phase boundary conds. in
              data_size=data_size+Q_LEVELS*N_P_ROWS*N_P_POINTS
            end if
!           the size of the data to be sent from processor iproc

!           Check the buffer is big enough
            IF (FULL_LENRIMDATA .LT.
     &          data_size) THEN
              WRITE(6,*) 'ERROR Buffer not big enough in GATHER_LBCS'
              WRITE(6,*) 'Buffer size is ',FULL_LENRIMDATA
              WRITE(6,*) 'Required size is ',data_size
              ICODE=1
              CMESSAGE='GATHER_LBCS BUFFER TOO SMALL'
              GOTO 9999
            ENDIF


!      2.2.2 Pack all the data for this boundary into the buf array

              CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! gather
              info=GC_NONE

            IF (mype .EQ. iproc) THEN  ! I'm processor iproc

              buf_pt=1

!             --- PSTAR ---

              DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                DO POINT=local_START_POINT,
     &                   local_START_POINT+N_P_POINTS-1
                  local_buf(buf_pt)=PART_LBC(local_PSTAR_START-1+
     &                                  local_P_bound_off+
     &                                  POINT+
     &                                  (ROW-1)*local_P_ROW_LEN)
                  buf_pt=buf_pt+1
                ENDDO
              ENDDO


!             --- U and V winds ---

              DO LEVEL=1,P_LEVELS
                DO ROW=local_START_ROW,local_START_ROW+N_U_ROWS-1
                  DO POINT=local_START_POINT,
     &                     local_START_POINT+N_U_POINTS-1
                    local_buf(buf_pt)=PART_LBC(local_U_START-1+
     &                                   local_U_bound_off+
     &                                   POINT+
     &                                   (ROW-1)*local_U_ROW_LEN+
     &                                   (LEVEL-1)*local_U_LENRIM)
                    buf_pt=buf_pt+1
                    local_buf(buf_pt)=PART_LBC(local_V_START-1+
     &                                   local_U_bound_off+
     &                                   POINT+
     &                                   (ROW-1)*local_U_ROW_LEN+
     &                                   (LEVEL-1)*local_U_LENRIM)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO


!             --- Theta ---

              DO LEVEL=1,P_LEVELS
                DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                  DO POINT=local_START_POINT,
     &                     local_START_POINT+N_P_POINTS-1
                    local_buf(buf_pt)=PART_LBC(local_THETA_START-1+
     &                                   local_P_bound_off+
     &                                   POINT+
     &                                   (ROW-1)*local_P_ROW_LEN+
     &                                   (LEVEL-1)*local_P_LENRIM)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO


!             --- Q ---

              DO LEVEL=1,Q_LEVELS
                DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                  DO POINT=local_START_POINT,
     &                     local_START_POINT+N_P_POINTS-1
                    local_buf(buf_pt)=PART_LBC(local_Q_START-1+
     &                                   local_P_bound_off+
     &                                   POINT+
     &                                   (ROW-1)*local_P_ROW_LEN+
     &                                   (LEVEL-1)*local_P_LENRIM)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO

!             --- Tracer Variables ---

              DO TVAR=1,TR_VARS
                DO LEVEL=1,TR_LEVELS
                  DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                    DO POINT=local_START_POINT,
     &                       local_START_POINT+N_P_POINTS-1
                      local_buf(buf_pt)=PART_LBC(local_TR_START-1+
     &                              local_P_bound_off+
     &                              POINT+
     &                              (ROW-1)*local_P_ROW_LEN+
     &                              (LEVEL-1)*local_P_LENRIM+
     &                              (TVAR-1)*TR_LEVELS*local_P_LENRIM)
                      buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO

            IF (L_LSPICE_BDY) THEN ! Mixed phase boundary conds. in
!             --- QCF ---

              DO LEVEL=1,Q_LEVELS
                DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                  DO POINT=local_START_POINT,
     &                     local_START_POINT+N_P_POINTS-1
                    local_buf(buf_pt)=PART_LBC(local_QCF_START-1+
     &                                   local_P_bound_off+
     &                                   POINT+
     &                                   (ROW-1)*local_P_ROW_LEN+
     &                                   (LEVEL-1)*local_P_LENRIM)
                   buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO
            END IF


!      2.2.3 Send local_buf array to processor GATHER_PE
!            into array buf_expand

              CALL GC_RSEND(iproc+1000*bound_type,data_size,
     &                      GATHER_PE,info,buf_expand,local_buf)

            ENDIF  ! If I'm processor iproc

            CALL GC_SSYNC(nproc,info)



              CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! gather
              info=GC_NONE
            IF (mype .EQ. GATHER_PE) THEN
!             I'm the processor who will gather together all the local
!             LBCs and assemble them into a global version

              CALL GC_RRECV(iproc+1000*bound_type,data_size,
     &                      iproc,info,buf_expand,local_buf)

!      2.2.5 Unpack buf_expand into the FULL_LBC array

              buf_pt=1

!             --- PSTAR ---

              DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                DO POINT=global_START_POINT,
     &                   global_START_POINT+N_P_POINTS-1
                  FULL_LBC(global_PSTAR_START-1+
     &                     global_P_bound_off+
     &                     POINT+
     &                     (ROW-1)*global_P_ROW_LEN)=
     &              buf_expand(buf_pt)
                  buf_pt=buf_pt+1
                ENDDO
              ENDDO


!             --- U and V winds ---

              DO LEVEL=1,P_LEVELS
                DO ROW=global_START_ROW,global_START_ROW+N_U_ROWS-1
                  DO POINT=global_START_POINT,
     &                     global_START_POINT+N_U_POINTS-1
                    FULL_LBC(global_U_START-1+
     &                       global_U_bound_off+
     &                       POINT+
     &                       (ROW-1)*global_U_ROW_LEN+
     &                       (LEVEL-1)*global_U_LENRIM)=
     &                buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                    FULL_LBC(global_V_START-1+
     &                       global_U_bound_off+
     &                       POINT+
     &                       (ROW-1)*global_U_ROW_LEN+
     &                       (LEVEL-1)*global_U_LENRIM)=
     &                buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO


!             -- Theta ---

              DO LEVEL=1,P_LEVELS
                DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                  DO POINT=global_START_POINT,
     &                     global_START_POINT+N_P_POINTS-1
                    FULL_LBC(global_THETA_START-1+
     &                       global_P_bound_off+
     &                       POINT+
     &                       (ROW-1)*global_P_ROW_LEN+
     &                       (LEVEL-1)*global_P_LENRIM)=
     &                buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO


!           --- Q ---

              DO LEVEL=1,Q_LEVELS
                DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                  DO POINT=global_START_POINT,
     &                     global_START_POINT+N_P_POINTS-1
                    FULL_LBC(global_Q_START-1+
     &                       global_P_bound_off+
     &                       POINT+
     &                       (ROW-1)*global_P_ROW_LEN+
     &                       (LEVEL-1)*global_P_LENRIM)=
     &                buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO


!             --- Tracer Variables ---

              DO TVAR=1,TR_VARS
                DO LEVEL=1,TR_LEVELS
                  DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                    DO POINT=global_START_POINT,
     &                       global_START_POINT+N_P_POINTS-1
                      FULL_LBC(global_TR_START-1+
     &                         global_P_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_P_ROW_LEN+
     &                         (LEVEL-1)*global_P_LENRIM+
     &                         (TVAR-1)*TR_LEVELS*global_P_LENRIM)=
     &                   buf_expand(buf_pt)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO

            IF (L_LSPICE_BDY) THEN ! Mixed phase boundary conds. in
!           --- QCF ---

              DO LEVEL=1,Q_LEVELS
                DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                  DO POINT=global_START_POINT,
     &                     global_START_POINT+N_P_POINTS-1
                    FULL_LBC(global_QCF_START-1+
     &                       global_P_bound_off+
     &                       POINT+
     &                       (ROW-1)*global_P_ROW_LEN+
     &                       (LEVEL-1)*global_P_LENRIM)=
     &                buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
              ENDDO
            END IF

            ENDIF  ! If I'm processor GATHER_PE

            CALL GC_SSYNC(nproc,info)

          ENDIF ! If this processor has boundary type bound_type

        ENDDO ! bound_type: loop over boundary types

      ENDDO ! iproc : loop over processors


 9999 CONTINUE  ! point to jump to if failure

      RETURN
      END
