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
! Subroutine ATMOS_ANG_MOM----------------------------------------
! Description:
!  Routine to calculate the 3 components of atmospheric angular
!  momentum for the wind term and the mass term separately.
!
! Method: The 3 components of angular momentum are defined below
!         as integrals over pressure of
!
! W1=  {u*sin(lat)cos(lon)-v*sin(lon)} r**3 cos(lat) dp/g
! W2=  {u*sin(lat)sin(lon)+v*cos(lon)} r**3 cos(lat) dp/g
! W3= -u*cos(lat) r**3 cos(lat) dp/g
!
! M1=  {r*omega*cos(lat)sin(lat)cos(lon)} r**3 cos(lat) dp/g
! M2=  {r*omega*cos(lat)sin(lat)sin(lon)} r**3 cos(lat) dp/g
! M3= -{r*omega*cos(lat)}cos(lat) r**3 cos(lat) dp/g
!
!
! Current Code Owner: R A Stratton
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.0     11/11/94  Original code. R A Stratton
!  4.4     02/07/97  Added ARGFLDPT args and MPP code    P.Burton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: ?
! System Task:              ?
! --------------------------------------------------------------------
      SUBROUTINE ATMOS_ANG_MOM(
     &         P_FIELD,U_FIELD,P_ROWS,ROW_LENGTH,P_LEVS,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     &         EW_SPACE,NS_SPACE,FIRST_LAT,FIRST_LONG,
     &         PSTAR,U,V,RS,COS_U_LATITUDE,DELTA_AK,DELTA_BK,
     &         L_AMM1,L_AMM2,L_AMM3,L_AMW1,L_AMW2,L_AMW3,
     &         AMM1,AMM2,AMM3,AMW1,AMW2,AMW3)
      IMPLICIT NONE
! Declarations:
!
! Global variables

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
C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
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


! Subroutine arguments
      INTEGER
     &   P_FIELD      ! IN : length of p grid
     &  ,U_FIELD      ! IN : length of u grid
     &  ,P_ROWS       ! IN : number of rows p grid
     &  ,ROW_LENGTH   ! IN : length of row
     &  ,P_LEVS       ! IN : number of model levels

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
      REAL
     &   EW_SPACE     ! IN : East west grid spacing in degrees
     &  ,NS_SPACE     ! IN : North South grid spacing in degrees
     &  ,FIRST_LAT    ! IN : first latitude
     &  ,FIRST_LONG   ! IN : first longitude
      REAL
     &  PSTAR(P_FIELD)     ! IN : pstar
     & ,U(U_field,P_LEVS)  ! IN : u
     & ,V(U_field,P_LEVS)  ! IN : V
     & ,RS(P_field,P_LEVS) ! IN : effect radius of atmosphere
     & ,COS_U_LATITUDE(U_FIELD) ! IN : cos (lat) u-grid
     & ,DELTA_AK(P_LEVS)   ! IN : layer akh(k+1)-akh(k)
     & ,DELTA_BK(P_LEVS)   ! IN : layer bkh(k+1)-bkh(k)

      LOGICAL
     &  L_AMM1        ! IN : true if field required
     & ,L_AMM2        ! IN : true if field required
     & ,L_AMM3        ! IN : true if field required
     & ,L_AMW1        ! IN : true if field required
     & ,L_AMW2        ! IN : true if field required
     & ,L_AMW3        ! IN : true if field required

      REAL
     &  AMM1(U_FIELD)      ! OUT: 1st com of angular momemtum mass term
     & ,AMM2(U_FIELD)      ! OUT: 2nd com of angular momemtum mass term
     & ,AMM3(U_FIELD)      ! OUT: 3rd com of angular momemtum mass term
     & ,AMW1(U_FIELD)      ! OUT: 1st com of angular momemtum wind term
     & ,AMW2(U_FIELD)      ! OUT: 2nd com of angular momemtum wind term
     & ,AMW3(U_FIELD)      ! OUT: 3rd com of angular momemtum wind term

! -------------------------------------------------------------------
! Local variables:

      REAL
     &   R3DP                  ! r**3 dp/g
     &  ,DP                    !  dp
     &  ,ROCOS                 ! r omega cos(lat)
     &  ,FACTOR                ! scaling factor /g
     &  ,COS_LONG
     &  ,SIN_LONG
     &  ,SIN_LAT
      REAL
     &   RS_U(U_FIELD)         ! effective radius on u grid
     &  ,PSTAR_U(U_FIELD)      ! pstar on u grid
     &  ,COSSQ(U_FIELD)        ! cos**2
     &  ,LONGITUDE(U_FIELD)    ! longitude
     &  ,LATITUDE(U_FIELD)     ! latitude
     &  ,SLCP(U_FIELD)         ! sin(lon)cos(lat)
     &  ,CLCP(U_FIELD)         ! cos(lon)cos(lat)
     &  ,SPSLCP(U_FIELD)       ! sin(lat)sin(lon)cos(lat)
     &  ,SPCLCP(U_FIELD)       ! sin(lat)cos(lon)cos(lat)

      INTEGER
     &   I,J,K,II         ! loop counters

! Function & Subroutine calls:
      External    p_to_uv

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
! ------------------------------------------------------------------
! constants
      FACTOR=1.e-24/g

! Calculate longitude  & latitude
      DO I=1,P_ROWS-1
        DO J=1,ROW_LENGTH
          II=J+(I-1)*ROW_LENGTH
          LONGITUDE(II)=(FIRST_LONG+EW_SPACE*
     &                  ((J+datastart(1)-Offx-1)-0.5))*PI_OVER_180
          LATITUDE(II)=(FIRST_LAT-NS_SPACE*
     &                  ((I+datastart(2)-Offy-1)-0.5))*PI_OVER_180
        ENDDO
      ENDDO

! calculate pstar on u grid

      CALL P_TO_UV(PSTAR,PSTAR_U,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)

      DO I=1,U_FIELD
!   Intialise output arrays
        IF (L_AMW1)  AMW1(I)=0.0
        IF (L_AMW2)  AMW2(I)=0.0
        IF (L_AMW3)  AMW3(I)=0.0
        IF (L_AMM1)  AMM1(I)=0.0
        IF (L_AMM2)  AMM2(I)=0.0
        IF (L_AMM3)  AMM3(I)=0.0
      ENDDO
      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        COSSQ(I)=COS_U_LATITUDE(I)*COS_U_LATITUDE(I)
      ENDDO
! calculate cos , sin etc
      IF (L_AMM1.OR.L_AMM2.OR.L_AMW1.OR.L_AMW2) THEN
       DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        COS_LONG=COS(LONGITUDE(I))
        SIN_LONG=SIN(LONGITUDE(I))
        SIN_LAT=SIN(LATITUDE(I))
        spclcp(i)=SIN_LAT*COS_LONG*COS_U_LATITUDE(I)
        spslcp(i)=SIN_LAT*SIN_LONG*COS_U_LATITUDE(I)
        clcp(i)=COS_LONG*COS_U_LATITUDE(I)
        slcp(i)=SIN_LONG*COS_U_LATITUDE(I)
       ENDDO
      ENDIF

! integrate momemtum over p

      DO K=1,P_LEVS     ! loop over model levels
        CALL P_TO_UV(RS(1,K),RS_U,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)

        DO I=FIRST_FLD_PT,LAST_U_FLD_PT

          DP=DELTA_AK(K) + DELTA_BK(K)*PSTAR_U(I)
          R3DP=(RS_U(I)**3)*DP*FACTOR
          ROCOS=OMEGA*RS_U(I)*COS_U_LATITUDE(I)

          IF (L_AMW1) AMW1(I)=AMW1(I)
     &       +(U(I,K)*SPCLCP(I)-V(I,K)*SLCP(I))*R3DP

          IF (L_AMW2) AMW2(I)=AMW2(I)
     &       +(U(I,K)*SPSLCP(I)+V(I,K)*CLCP(I))*R3DP

          IF (L_AMW3) AMW3(I)=AMW3(I) - U(I,K)*COSSQ(I)*R3DP

          IF (L_AMM1) AMM1(I)=AMM1(I) + ROCOS*SPCLCP(I)*R3DP
          IF (L_AMM2) AMM2(I)=AMM2(I) + ROCOS*SPSLCP(I)*R3DP
          IF (L_AMM3) AMM3(I)=AMM3(I) - ROCOS*COSSQ(I)*R3DP
        ENDDO           ! end loop over gridpoints
      ENDDO             ! end loop over model levels

      RETURN
      END
