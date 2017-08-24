C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Gathers in ocean partial LBCs from boundaries to create
!  global LBC
!
! Subroutine Interface
      SUBROUTINE GATHER_OCEAN_LBCS(FULL_LBC,FULL_LENRIMDATA,
     &                             PART_LBC,PART_LENRIMDATA,
     &                             GATHER_PE,
     &                             ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Description:
! Gathers ocean LBCs from relevant processors and assembles
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
! The East and West boundaries of PART_LBC are dimensioned to
! have the same number of rows of data as the prognostic fields to
! which they relate (with halo rows excluded).
!
! History
!  Model    Date       Modification history from model version 4.5
!  version
!   4.5   17/06/98    New Deck for MPP code. S.Ineson,M.Bell,P.Burton
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


! Local variables

      INTEGER
     &  global_ROW_LENGTH ! length of row on global P grid
     &, global_P_ROWS     ! number of rows on global P grid
     &, local_ROW_LENGTH  ! length of row on local P grid (no halos)
     &, local_U_ROW_LENGTH ! length of row on local U grid (no halos)
     &, local_P_ROWS      ! number of rows on local P grid (no halos)
     &, local_U_ROWS      ! number of rows on local U grid (no halos)
     &, global_P_LENRIM   ! size of boundary data for single
!                         !   level on P grid for global data
     &, global_U_LENRIM   ! size of boundary data for single
!                         !   level on U grid for global data
     &, local_P_LENRIM    ! size of boundary data for single
!                         !    level on P grid for local data
     &, local_U_LENRIM    ! size of boundary data for single
!                         !    level on U grid for local data

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


     &, global_TRACER_START     ! Start addresses of variables
     &, global_U_START          ! in the global LBC data
     &, global_V_START          !
     &, global_STREAM_START     !
     &, global_ICE_START
     &, global_P_EAST_DATA_off  ! offset of East P data in global LBC
     &, global_U_EAST_DATA_off  ! offset of East U data in global LBC
     &, global_P_SOUTH_DATA_off ! offset of South P data in global LBC
     &, global_U_SOUTH_DATA_off ! offset of South U data in global LBC
     &, global_P_WEST_DATA_off  ! offset of West P data in global LBC
     &, global_U_WEST_DATA_off  ! offset of West U data in global LBC

     &, local_TRACER_START      ! Start addresses of variables
     &, local_U_START           ! in the local LBC data
     &, local_V_START           !
     &, local_STREAM_START      !
     &, local_ICE_START
     &, local_P_EAST_DATA_off   ! offset of East P data in local LBC
     &, local_U_EAST_DATA_off   ! offset of East U data in local LBC
     &, local_P_SOUTH_DATA_off  ! offset of South P data in local LBC
     &, local_U_SOUTH_DATA_off  ! offset of South U data in local LBC
     &, local_P_WEST_DATA_off   ! offset of West P data in local LBC
     &, local_U_WEST_DATA_off   ! offset of West U data in local LBC

     &, offset,offset_p,offset_u   ! used in calculation of offsets
     &, numside_colso           !no of active boundary columns
     &, numside_rowso           !no of active boundary rows

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

! 1.0 Set up sizes and addresses for global data

! Set up some sizes

      global_ROW_LENGTH=glsize(1)
      global_P_ROWS=glsize(2)

! Find no of active rows, columns and set addresses for grid location

      numside_rowso = 0
      numside_colso = 0

      IF ( L_OBDY_NORTH) THEN
        offset_p = global_ROW_LENGTH * RIMWIDTHO
        offset_u = (global_ROW_LENGTH-1) * RIMWIDTHO
        numside_rowso = numside_rowso + rimwidtho
      ELSE
        offset_p = 0
        offset_u = 0
      ENDIF

      IF ( L_OBDY_EAST) THEN
        global_P_EAST_DATA_off = offset_p
        global_U_EAST_DATA_off = offset_u
        offset_p = offset_p + global_P_ROWS * RIMWIDTHO
        offset_u = offset_u + (global_P_ROWS-1) * RIMWIDTHO
        numside_colso = numside_colso + rimwidtho
      ELSE
        global_P_EAST_DATA_off = -1
        global_U_EAST_DATA_off = -1
      ENDIF

      IF ( L_OBDY_SOUTH) THEN
        global_P_SOUTH_DATA_off = offset_p
        global_U_SOUTH_DATA_off = offset_u
        offset_p = offset_p + global_ROW_LENGTH * RIMWIDTHO
        offset_u = offset_u + (global_ROW_LENGTH-1) * RIMWIDTHO
        numside_rowso = numside_rowso + rimwidtho
      ELSE
        global_P_SOUTH_DATA_off = -1
        global_U_SOUTH_DATA_off = -1
      ENDIF

      IF ( L_OBDY_WEST) THEN
        global_P_WEST_DATA_off = offset_p
        global_U_WEST_DATA_off = offset_u
        numside_colso = numside_colso + rimwidtho
      ELSE
        global_P_WEST_DATA_off = -1
        global_U_WEST_DATA_off = -1
      ENDIF

! Set up LENRIM sizes

      global_P_LENRIM = global_ROW_LENGTH*numside_rowso
     &                     + global_P_ROWS* numside_colso

      global_U_LENRIM = (global_ROW_LENGTH-1) *numside_rowso
     &                     + (global_P_ROWS-1) * numside_colso

! Set up addresses for variable types

      IF (L_OBDY_TRACER) THEN
        global_TRACER_START = 1
        offset = global_P_LENRIM * KM_UI * NT_UI
      ELSE
        global_TRACER_START = -1
        offset = 0
      ENDIF

      IF (L_OBDY_UV) THEN
        global_U_START = offset + 1
        global_V_START = offset + global_U_LENRIM * KM_UI + 1
        offset = offset + 2 * global_U_LENRIM * KM_UI
      ELSE
        global_U_START = -1
        global_V_START = -1
      ENDIF

      IF (L_OBDY_STREAM) THEN
        global_STREAM_START = offset + 1
        offset = offset + 2 * global_P_LENRIM
      ELSE
        global_STREAM_START = -1
      ENDIF

      IF (L_OBDY_ICE) THEN
        global_ICE_START = offset + 1
      ELSE
        global_ICE_START = -1
      ENDIF

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
          local_U_ROW_LENGTH=local_ROW_LENGTH-1
        ELSE
          local_U_ROW_LENGTH=local_ROW_LENGTH
        ENDIF

        local_P_ROWS=g_blsizep(2,iproc)      ! again, no halos
        IF (iproc_at_base) THEN
          local_U_ROWS=local_P_ROWS-1
        ELSE
          local_U_ROWS=local_P_ROWS
        ENDIF

        local_P_LENRIM=local_ROW_LENGTH * numside_rowso
     &                     +local_P_ROWS * numside_colso
        local_U_LENRIM=local_U_ROW_LENGTH * numside_rowso
     &                     +local_U_ROWS * numside_colso

! Set up some addresses for variable types

        IF (L_OBDY_TRACER) THEN
          local_TRACER_START = 1
          offset = local_P_LENRIM * KM_UI
        ELSE
          local_TRACER_START = -1
          offset = 0
        ENDIF

        IF (L_OBDY_UV) THEN
          local_U_START = offset + 1
          local_V_START = offset + local_U_LENRIM * KM_UI + 1
          offset = offset + 2 * local_U_LENRIM * KM_UI
        ELSE
          local_U_START = -1
          local_V_START = -1
        ENDIF

        IF (L_OBDY_STREAM) THEN
          local_STREAM_START = offset + 1
          offset = offset + 2 * local_P_LENRIM
        ELSE
          local_STREAM_START = -1
        ENDIF

        IF (L_OBDY_ICE) THEN
          local_ICE_START = offset + 1
          offset = offset + 3 * local_P_LENRIM
        ELSE
          local_ICE_START = -1
        ENDIF

! Set up some addresses for grid location

        IF ( L_OBDY_NORTH) THEN
          offset_p = local_ROW_LENGTH * RIMWIDTHO
          offset_u = local_U_ROW_LENGTH * RIMWIDTHO
        ELSE
          offset_p = 0
          offset_u = 0
        ENDIF

        IF ( L_OBDY_EAST) THEN
          local_P_EAST_DATA_off = offset_p
          local_U_EAST_DATA_off = offset_u
          offset_p = offset_p + local_P_ROWS * RIMWIDTHO
          offset_u = offset_u + local_U_ROWS * RIMWIDTHO
        ELSE
          local_P_EAST_DATA_off = -1
          local_U_EAST_DATA_off = -1
        ENDIF

        IF ( L_OBDY_SOUTH) THEN
          local_P_SOUTH_DATA_off = offset_p
          local_U_SOUTH_DATA_off = offset_u
          offset_p = local_ROW_LENGTH * RIMWIDTHO
          offset_u = local_U_ROW_LENGTH * RIMWIDTHO
        ELSE
          local_P_SOUTH_DATA_off = -1
          local_U_SOUTH_DATA_off = -1
        ENDIF

        IF ( L_OBDY_WEST) THEN
          local_P_WEST_DATA_off = offset_p
          local_U_WEST_DATA_off = offset_u
        ELSE
          local_P_WEST_DATA_off = -1
          local_U_WEST_DATA_off = -1
        ENDIF

!-------------------------------------------------------------------

!   2.2 Loop over boundaries: North, East, South and West

!   For the ocean model the order of the boundary data order is N,E,S,W.
!   Howvever, row 1 of the ocean model is the the southern most row.
!   Hence, the northern bdy data, bt_top ==> processors, iproc_at_base
!          the southern bdy data,bt_base ==> processors, iproc_at_top

        DO bound_type=bt_top,bt_left  ! loop over all boundaries

          IF (((bound_type .EQ. bt_top) .AND. (iproc_at_base)
     &              .AND. (L_OBDY_NORTH))    .OR.
     &       ((bound_type .EQ. bt_right) .AND. (iproc_at_right)
     &              .AND. (L_OBDY_EAST))      .OR.
     &       ((bound_type .EQ. bt_base) .AND. (iproc_at_top)
     &              .AND. (L_OBDY_SOUTH))      .OR.
     &       ((bound_type .EQ. bt_left) .AND. (iproc_at_left)
     &              .AND. (L_OBDY_WEST)))             THEN
!             Processor iproc has a boundary of type bound_type


!      2.2.1 Set up data pointers and sizes for this boundary

            IF ((bound_type .EQ. bt_top) .OR.   ! What type of
     &          (bound_type .EQ. bt_base)) THEN ! boundary is it?

!             Northern or Southern boundary

              global_START_ROW=1
              local_START_ROW=1

              global_START_POINT=g_datastart(1,iproc)
              local_START_POINT=1

              N_P_ROWS=RIMWIDTHO
              N_P_POINTS=local_ROW_LENGTH
              N_U_ROWS=RIMWIDTHO
              N_U_POINTS=local_U_ROW_LENGTH

              global_P_ROW_LEN=global_ROW_LENGTH
              local_P_ROW_LEN=local_ROW_LENGTH
              global_U_ROW_LEN=global_ROW_LENGTH-1
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

              global_START_ROW=g_datastart(2,iproc)
              local_START_ROW=1

              global_START_POINT=1
              local_START_POINT=1

              N_P_ROWS=local_P_ROWS
              N_P_POINTS=RIMWIDTHO
              N_U_ROWS=local_U_ROWS
              N_U_POINTS=RIMWIDTHO

              global_P_ROW_LEN=RIMWIDTHO
              local_P_ROW_LEN=RIMWIDTHO
              global_U_ROW_LEN=RIMWIDTHO
              local_U_ROW_LEN=RIMWIDTHO

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

!  Find size of the data to be sent to processor iproc

            data_size=0
            IF (L_OBDY_TRACER) data_size =
     &           data_size + KM_UI * NT_UI * N_P_ROWS * N_P_POINTS
            IF (L_OBDY_UV) data_size =
     &           data_size + 2 * KM_UI * N_U_ROWS * N_U_POINTS
            IF (L_OBDY_STREAM) data_size =
     &           data_size + 2 * N_P_ROWS * N_P_POINTS
            IF (L_OBDY_ICE) data_size =
     &           data_size + 3 * N_P_ROWS * N_P_POINTS

!           Check the buffer is big enough
            IF (FULL_LENRIMDATA .LT.
     &          data_size) THEN
              WRITE(6,*) 'ERROR Buffer not big enough in GATHER_LBCS'
              WRITE(6,*) 'Buffer size is ',FULL_LENRIMDATA
              WRITE(6,*) 'Required size is ',data_size
              ICODE=1
              CMESSAGE='SCATTER_LBCS BUFFER TOO SMALL'
              GOTO 9999
            ENDIF


!      2.2.2 Pack all the data for this boundary into the buf array

              CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)  ! gather
              info=GC_NONE

            IF (mype .EQ. iproc) THEN  ! I'm processor iproc

              buf_pt=1

!             --- TRACERS ---

              IF (L_OBDY_TRACER) THEN

              DO TVAR=1,NT_UI
                DO LEVEL=1,KM_UI
                  DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                    DO POINT=local_START_POINT,
     &                       local_START_POINT+N_P_POINTS-1
                      local_buf(buf_pt)=PART_LBC(local_TRACER_START-1+
     &                              local_P_bound_off+
     &                              POINT+
     &                              (ROW-1)*local_P_ROW_LEN+
     &                              (LEVEL-1)*local_P_LENRIM+
     &                              (TVAR-1)*KM_UI*local_P_LENRIM)
                      buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO

               ENDIF      !  L_OBDY_TRACER

!             --- U,V CURRENTS ---

               IF (L_OBDY_UV) THEN

                 DO LEVEL=1,KM_UI
                   DO ROW=local_START_ROW,local_START_ROW+N_U_ROWS-1
                     DO POINT=local_START_POINT,
     &                       local_START_POINT+N_U_POINTS-1
                       local_buf(buf_pt)=PART_LBC(local_U_START-1+
     &                              local_U_bound_off+
     &                              POINT+
     &                              (ROW-1)*local_U_ROW_LEN+
     &                               (LEVEL-1)*local_U_LENRIM)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO

                 DO LEVEL=1,KM_UI
                   DO ROW=local_START_ROW,local_START_ROW+N_U_ROWS-1
                     DO POINT=local_START_POINT,
     &                       local_START_POINT+N_U_POINTS-1
                       local_buf(buf_pt)=PART_LBC(local_V_START-1+
     &                              local_U_bound_off+
     &                              POINT+
     &                              (ROW-1)*local_U_ROW_LEN+
     &                              (LEVEL-1)*local_U_LENRIM)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO

               ENDIF      !  L_OBDY_UV

!             --- STREAMFUNCTION ---

               IF (L_OBDY_STREAM) THEN

                DO TVAR=1,2
                 DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                   DO POINT=local_START_POINT,
     &                       local_START_POINT+N_P_POINTS-1
                     local_buf(buf_pt)=PART_LBC(local_STREAM_START-1+
     &                              local_P_bound_off+
     &                              POINT+
     &                              (ROW-1)*local_P_ROW_LEN+
     &                              (TVAR-1)*local_P_LENRIM)
                     buf_pt=buf_pt+1
                   ENDDO
                 ENDDO
                ENDDO

               ENDIF      !  L_OBDY_STREAM

!             --- SEA ICE ---

              IF (L_OBDY_ICE) THEN

               DO TVAR=1,3
                DO ROW=local_START_ROW,local_START_ROW+N_P_ROWS-1
                  DO POINT=local_START_POINT,
     &                 local_START_POINT+N_P_POINTS-1
                    local_buf(buf_pt)= PART_LBC(local_ICE_START-1+
     &                             local_P_bound_off+
     &                             POINT+
     &                            (ROW-1)*local_P_ROW_LEN+
     &                            (TVAR-1)*local_P_LENRIM)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
               ENDDO

              ENDIF     !  L_OBDY_ICE

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

!             --- TRACERS ---

              IF (L_OBDY_TRACER) THEN

              DO TVAR=1,NT_UI
                DO LEVEL=1,KM_UI
                  DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                    DO POINT=global_START_POINT,
     &                       global_START_POINT+N_P_POINTS-1
                      FULL_LBC(global_TRACER_START-1+
     &                         global_P_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_P_ROW_LEN+
     &                         (LEVEL-1)*global_P_LENRIM+
     &                         (TVAR-1)*KM_UI*global_P_LENRIM)=
     &                   buf_expand(buf_pt)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO

               ENDIF    ! L_OBDY_TRACER

!             --- U,V CURRENTS ---

              IF (L_OBDY_UV) THEN

                DO LEVEL=1,KM_UI
                  DO ROW=global_START_ROW,global_START_ROW+N_U_ROWS-1
                    DO POINT=global_START_POINT,
     &                       global_START_POINT+N_U_POINTS-1
                      FULL_LBC(global_U_START-1+
     &                         global_U_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_U_ROW_LEN+
     &                         (LEVEL-1)*global_U_LENRIM)=
     &                   buf_expand(buf_pt)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO

                DO LEVEL=1,KM_UI
                  DO ROW=global_START_ROW,global_START_ROW+N_U_ROWS-1
                    DO POINT=global_START_POINT,
     &                       global_START_POINT+N_U_POINTS-1
                      FULL_LBC(global_V_START-1+
     &                         global_U_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_U_ROW_LEN+
     &                         (LEVEL-1)*global_U_LENRIM)=
     &                   buf_expand(buf_pt)
                       buf_pt=buf_pt+1
                     ENDDO
                   ENDDO
                 ENDDO

               ENDIF    ! L_OBDY_UV

!             --- STREAMFUNCTION ---

              IF (L_OBDY_STREAM) THEN

               DO TVAR = 1,2
                DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                  DO POINT=global_START_POINT,
     &                       global_START_POINT+N_P_POINTS-1
                    FULL_LBC(global_STREAM_START-1+
     &                         global_P_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_P_ROW_LEN+
     &                         (TVAR-1)*global_P_LENRIM)=
     &                   buf_expand(buf_pt)
                     buf_pt=buf_pt+1
                   ENDDO
                 ENDDO
               ENDDO

               ENDIF    ! L_OBDY_STREAM

!             --- SEA ICE ---

              IF (L_OBDY_ICE) THEN

               DO TVAR=1,3
                DO ROW=global_START_ROW,global_START_ROW+N_P_ROWS-1
                  DO POINT=global_START_POINT,
     &                 global_START_POINT+N_P_POINTS-1
                    FULL_LBC(global_ICE_START-1+
     &                         global_P_bound_off+
     &                         POINT+
     &                         (ROW-1)*global_P_ROW_LEN+
     &                         (TVAR-1)*global_P_LENRIM)=
     &                   buf_expand(buf_pt)
                    buf_pt=buf_pt+1
                  ENDDO
                ENDDO
               ENDDO

              ENDIF      !  L_OBDY_ICE

            ENDIF  ! If I'm processor GATHER_PE

            CALL GC_SSYNC(nproc,info)

          ENDIF ! If this processor has boundary type bound_type

        ENDDO ! bound_type: loop over boundary types

      ENDDO ! iproc : loop over processors


 9999 CONTINUE  ! point to jump to if failure

      RETURN
      END
