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
CLL  SUBROUTINES DUCT, STRETCH, SMOOTH, LUBATT--------------------------
CLL
CLL                  PURPOSES:
CLL    DUCT: TO CALCULTE RADIO DUCT INTENSITY AND HEIGHT
CLL           SUITABLE FOR ROTATED GRIDS
CLL
CLL  STRETCH: CALCULATES SEA TEMPERATURE OVER THE OPEN OCEAN WITH
CLL           COASTAL LAND POINTS, GIVEN THE SEA TEMPERATURE FROM
CLL           ADJACENT SEA POINTS
CLL
CLL  SMOOTH:  CALCULATES WEIGHTED SMOOTHED VALUES FOR ALL INTERIOR
CLL           GRID POINTS
CLL
CLL  STRETCH, SMOOTH AND LUBATT  SUITABLE FOR ROTATED GRIDS
CLL
CLL  LUBATT:  CALCULATES EVAPORATION DUCT HEIGHT AND INTENSITY.
CLL
CLL
CLL                 DUCT:
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL                 STRETCH:
CLL P.Smith     <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 2, DATED 18/01/90
CLL
CLL Logical components covered : D484
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE DUCT
     1 (PSTAR,TSTAR,THETA,Q,U,V,AK,BK,AKH,BKH,LAND,
     3  DUCT_HEIGHT,MAX_WAVELENGTH,P_EXNER_HALF,
     4  P_LEVELS,Q_LEVELS,ROW_LENGTH,P_ROWS,U_ROWS,P_FIELD,U_FIELD)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
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

C This gives R,CP and KAPPA=R/CP
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C This gives EPSILON
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C This gives ZERODEGC
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
C This gives RMDI and IMDI
C-----------------------------------------------------------------------
      INTEGER
     *  P_LEVELS       ! IN  No of pressure levels
     *, Q_LEVELS       ! IN  No of humidity levels
     *, ROW_LENGTH     ! IN  No of points in row
     *, P_ROWS         ! IN  No of rows in pressure grid
     *, U_ROWS         ! IN  No of rows in (u,v) grid
     *, P_FIELD        ! IN  No of points in pressure grid
     *, U_FIELD        ! IN  No of points in (u,v) grid
C-----------------------------------------------------------------------
      REAL
     *  PSTAR(P_FIELD) ! IN  Pressure at the surface of the earth
     *, TSTAR(P_FIELD) ! IN  Temperature at the surface of the earth
     *, THETA(P_FIELD,P_LEVELS)! IN  Potential temperature
     *, P_EXNER_HALF(P_FIELD,P_LEVELS+1)! IN Exner pressure at half levs
     *, U(U_FIELD,P_LEVELS)    ! IN  Easterly component of wind
     *, V(U_FIELD,P_LEVELS)    ! IN  Northerly component of wind
     *, Q(P_FIELD,Q_LEVELS)    ! IN  Specific humidity
     *, DUCT_HEIGHT(P_FIELD)   ! OUT evaporation duct height in metres
     *, MAX_WAVELENGTH(P_FIELD)! OUT Maximum wavelength propagated by
     *                         ! duct (duct intensity) in metres
      REAL
     *  AK(P_LEVELS)     !IN  Hybrid Coords. A and B values for
     *, BK(P_LEVELS)     !IN  model full levels
     *, AKH(P_LEVELS+1)  !IN  Hybrid Coords. A and B values for
     *, BKH(P_LEVELS+1)  !IN  model half levels.
C-----------------------------------------------------------------------
      LOGICAL
     *  LAND(P_FIELD)  ! IN  TRUE if land; FALSE if sea
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL STRETCH,SMOOTH,LUBATT
C*----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C   LOCAL CONSTANTS
C-----------------------------------------------------------------------
      REAL
     *  RT         ! Used in the calculation of dew point
     *, R_OVER_M   ! gas constant/molecular weight of water
     *, RLO        ! latent heat of evaporation at 0 deg C
     *, RLE        ! rate of change of latent heat with temp at 0 deg C
      PARAMETER(RT=3.66E-3,
     &          R_OVER_M=461.5,
     &          RLO=2.5E6,
     &          RLE=-2.73E3)
C-----------------------------------------------------------------------
C   LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  I          ! Loop counter
     *, I1,J       ! Row number/position in row for WRITE statement
     *, ROW_NO,COL_NO  ! temporary
      REAL
     *  P_ON_LVL_1       ! Pressure on level 1
     *, P1,P2            ! Pressure on half levels 1 and 2
     *, P_EXNER_FULL     ! Exner pressure on full model level
     *, T_ON_LVL_1       ! Temperature on level 1
     *, HUMIDITY         ! Specific humidity on level 1
     *, VAPOUR_PRESSURE
     *, R1 ! Latent heat of evaporation
     *, R2
     *, T_DEW(P_FIELD)   ! dew point temperature
     *, T_DRY(P_FIELD)   ! dry bulb temperature
     *, WSPEED(P_FIELD)  ! windspeed of points on theta grid
     *, T_SEA(P_FIELD)   ! sea surface temperature
     *, SMOOTHED_T_SEA(P_FIELD) ! smoothed sea surface temperature
     *, SMOOTHED_T_DEW(P_FIELD) ! smoothed dewpoint temperature
     *, SMOOTHED_T_DRY(P_FIELD) ! smoothed dry bulb temperature
     *, SMOOTHED_WSPEED(P_FIELD)! smoothed windspeed
     *, TEMP(P_FIELD)
C
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
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------


C---    Put temperature to degree C-------------------------------------
      DO I=1,P_FIELD
        TEMP(I)=TSTAR(I)-ZERODEGC
      ENDDO
      DO 444 I=1,P_FIELD
C-----------------------------------------------------------------------
CL  2. Calculate dew point temperature
C-----------------------------------------------------------------------
        P_ON_LVL_1=AK(1)+BK(1)*PSTAR(I)
        HUMIDITY=Q(I,1)    ! Units KG/KG also press in Pascals
        IF(HUMIDITY.LE.1.0E-3) HUMIDITY=1.0E-3
        VAPOUR_PRESSURE=(HUMIDITY*PSTAR(I)/100.)/(EPSILON+HUMIDITY)
        VAPOUR_PRESSURE=VAPOUR_PRESSURE/6.11
        P1 = AKH(1) + BKH(1)*PSTAR(I)
        P2 = AKH(2) + BKH(2)*PSTAR(I)
        P_EXNER_FULL = P_EXNER_C
     +  (P_EXNER_HALF(I,2),P_EXNER_HALF(I,1),P2,P1,KAPPA)
        T_ON_LVL_1 = THETA(I,1) * P_EXNER_FULL
        R1=RLO+RLE*(T_ON_LVL_1-ZERODEGC)  ! Latent heat for new temp
        R2=RT-LOG(VAPOUR_PRESSURE)*R_OVER_M/R1
        T_DEW(I)=1.0/R2
C-----------------------------------------------------------------------
CL  3. Store  dry bulb temperature
C-----------------------------------------------------------------------
        T_DRY(I)=T_ON_LVL_1
 444  CONTINUE
C-----------------------------------------------------------------------
CL  4. Calculate wind speed
C-----------------------------------------------------------------------
        DO I=1,U_FIELD
          TEMP(I)=SQRT(U(I,1)*U(I,1)+V(I,1)*V(I,1))
        ENDDO
        CALL UV_TO_P_FULL(TEMP,WSPEED,U_FIELD,P_FIELD,ROW_LENGTH,P_ROWS)
C-----------------------------------------------------------------------
CL  5. Calculate sea surface temperature
C-----------------------------------------------------------------------
      CALL STRETCH(TSTAR,T_SEA,LAND,ROW_LENGTH,P_ROWS,P_FIELD)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CL  6. Smooth sea surface temp, dry bulb temp, dew point temp and
CL     wind speed fields
C-----------------------------------------------------------------------
      CALL SWAPBOUNDS(T_SEA,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(T_DRY,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(T_DEW,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(WSPEED,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      WRITE(6,1011)
 1011 FORMAT(' Successfully reached the start of subroutine SMOOTH')
      CALL SMOOTH(T_SEA,SMOOTHED_T_SEA,ROW_LENGTH,P_ROWS,P_FIELD)
      CALL SMOOTH(T_DRY,SMOOTHED_T_DRY,ROW_LENGTH,P_ROWS,P_FIELD)
      CALL SMOOTH(T_DEW,SMOOTHED_T_DEW,ROW_LENGTH,P_ROWS,P_FIELD)
      CALL SMOOTH(WSPEED,SMOOTHED_WSPEED,ROW_LENGTH,P_ROWS,P_FIELD)
! Ensure all smoothed variables have halos filled.

      CALL SWAPBOUNDS(SMOOTHED_T_SEA,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(SMOOTHED_T_DRY,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(SMOOTHED_T_DEW,ROW_LENGTH,P_ROWS,Offx,Offy,1)
      CALL SWAPBOUNDS(SMOOTHED_WSPEED,ROW_LENGTH,P_ROWS,Offx,Offy,1)

C-----------------------------------------------------------------------
CL  7. Calculate DUCT height using smoothed parameters.
CL     Set DUCT height missing over land.
C-----------------------------------------------------------------------
      DO 460 I=1,P_FIELD
        IF (LAND(I)) THEN
          DUCT_HEIGHT(I)=RMDI
C         DUCT_HEIGHT(I)=0.
          MAX_WAVELENGTH(I)=RMDI
C         MAX_WAVELENGTH(I)=0.
        ELSE
          CALL LUBATT(SMOOTHED_WSPEED(I),SMOOTHED_T_DRY(I)-ZERODEGC,
     &    SMOOTHED_T_SEA(I)-ZERODEGC,SMOOTHED_T_DEW(I)-ZERODEGC,
     &    DUCT_HEIGHT(I),MAX_WAVELENGTH(I))
          J=MOD(I,ROW_LENGTH)
          I1=(I-J)/ROW_LENGTH+1
        ENDIF
 460  CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE STRETCH
     & (TSTAR,T_SEA,LAND,ROW_LENGTH,P_ROWS,P_FIELD)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER
     *  ROW_LENGTH        ! IN  No of columns of longitude
     *, P_ROWS            ! IN  No of rows of latitude
     *, P_FIELD           ! IN  No of points in pressure grid
C-----------------------------------------------------------------------
      REAL
     *  TSTAR(P_FIELD) ! IN  Temperature at the surface of the earth
     *, T_SEA(P_FIELD) ! OUT Sea temperature over open ocean
C-----------------------------------------------------------------------
      LOGICAL
     *  LAND(P_FIELD)     ! IN  TRUE if land; FALSE if sea
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
! Copies of input arrays but with extra large EW halos to allow the
! fourth order "stretching" operation

      REAL TSTAR_COPY((ROW_LENGTH+2)*P_ROWS)
      LOGICAL LAND_COPY((ROW_LENGTH+2)*P_ROWS)

C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
C     NONE
C*----------------------------------------------------------------------
C
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
C-----------------------------------------------------------------------
C   LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  I,J            ! Loop counters
     &, point_src,point_set   ! points of source data and to be set
     &, offset ! distance of POINT_SET from POINT_SRC
     &, max_off ! maximum distance allowed from POINT_SRC

      LOGICAL
     &  end

      INTEGER
     &  extended_ROW_LENGTH
     &, extended_P_FIELD
     &, extended_point_src
     &, extended_J


      max_off=2  ! Go up to two points to look for a value to use

! Copy input arrays into _COPY arrays which contain larger EW halos

      extended_ROW_LENGTH=ROW_LENGTH+2
      extended_P_FIELD=extended_ROW_LENGTH*P_ROWS

      CALL COPY_FIELD(TSTAR,TSTAR_COPY,
     &                P_FIELD,extended_P_FIELD,
     &                ROW_LENGTH,P_ROWS,1,
     &                Offx,Offy,
     &                2,Offy,.TRUE.)

      CALL COPY_FIELD(LAND,LAND_COPY,
     &                P_FIELD,extended_P_FIELD,
     &                ROW_LENGTH,P_ROWS,1,
     &                Offx,Offy,
     &                2,Offy,.TRUE.)


      DO I=1+Offy,P_ROWS-Offy ! miss halos

        DO J=1+Offx,ROW_LENGTH-Offx  ! miss halos

          extended_J=J+1

! Transfer sea temperatures 1 or 2 grid lengths westward onto
! land points

          offset=0
          end=.FALSE.

          point_set=(I-1)*ROW_LENGTH+J

          DO WHILE ((offset .LE. max_off) .AND. (.NOT. end))

            point_src=(I-1)*extended_ROW_LENGTH+
     &                MOD(extended_J+offset-1,extended_ROW_LENGTH)+1

            IF (.NOT. LAND_COPY(point_src)) THEN
              T_SEA(point_set)=TSTAR_COPY(point_src)
              end=.TRUE.
            ENDIF

            offset=offset+1

          ENDDO

          IF (.NOT. end) THEN
            T_SEA(point_set)=TSTAR(point_set)
          ENDIF

! Transfer sea temperatures 1 or 2 grid lengths eastwars onto
! land points

          IF (LAND(point_set)) THEN

            offset=0
            end=.FALSE.

            DO WHILE ((offset .LE. max_off) .AND. (.NOT. end))


              point_src=(I-1)*extended_ROW_LENGTH+
     &                  MOD(extended_J-offset+extended_ROW_LENGTH-1,
     &                      extended_ROW_LENGTH)+1

              IF (.NOT. (LAND_COPY(point_src))) THEN
                T_SEA(point_set)=TSTAR_COPY(point_src)
                end=.TRUE.
              ENDIF

              offset=offset+1

            ENDDO

          ENDIF

        ENDDO ! J : loop along row

      ENDDO ! I : loop over rows
       RETURN
       END
C-----------------------------------------------------------------------
C
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE SMOOTH
     & (A,ABAR,ROW_LENGTH,P_ROWS,P_FIELD)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER
     *  ROW_LENGTH    ! IN  No of points in a row
     *, P_ROWS        ! IN  No of rows of theta grid
     *, P_FIELD       ! IN  No of points in pressure field
C-----------------------------------------------------------------------
      REAL
     *  A(P_FIELD)    ! IN  Array to be smoothed
     *, ABAR(P_FIELD) ! OUT Smoothed array
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
C     NONE
C*----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C   LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  IW1,IW2,IW3,IW4,IW5,IW6  ! Weight constants
     *, I,J                      ! Loop counters
     *, POINT                    ! Point being smoothed
      INTEGER I_start,I_end

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
C-----------------------------------------------------------------------
CL  1. Set weight constants
C-----------------------------------------------------------------------
      IW1=1
      IW2=2
      IW3=8
      IW4=IW3+4*IW2+4*IW1
C-----------------------------------------------------------------------
CL  2. Smooth interior points of grid
C-----------------------------------------------------------------------
      IF (attop) THEN
        I_start=Offy+2
      ELSE
        I_start=Offy+1
      ENDIF
      IF (atbase) THEN
        I_end=P_ROWS-Offy-1
      ELSE
        I_end=P_ROWS-Offy
      ENDIF

      DO 300 I=I_start,I_end
        DO 200 J=2,ROW_LENGTH-1
          POINT=(I-1)*ROW_LENGTH+J
          ABAR(POINT)=(IW1*A(POINT-ROW_LENGTH-1)+IW2*A(POINT-ROW_LENGTH)
     1     +IW1*A(POINT-ROW_LENGTH+1)+IW2*A(POINT-1)+IW3*A(POINT)
     2     +IW2*A(POINT+1)+IW1*A(POINT+ROW_LENGTH-1)
     3     +IW2*A(POINT+ROW_LENGTH)+IW1*A(POINT+ROW_LENGTH+1))/IW4
 200    CONTINUE
 300  CONTINUE
C-----------------------------------------------------------------------
CL  3. Smooth top row (excluding the two end points)
C-----------------------------------------------------------------------
      IW6=IW3+3*IW2+2*IW1
      IF (attop) THEN
      I=Offy+1
      DO 400 J=2,ROW_LENGTH-1
        POINT=(I-1)*ROW_LENGTH+J
        ABAR(POINT)=(IW2*A(POINT-1)+IW3*A(POINT)+IW2*A(POINT+1)
     1   +IW1*A(POINT+ROW_LENGTH-1)+IW2*A(POINT+ROW_LENGTH)
     2   +IW1*A(POINT+ROW_LENGTH+1))/IW6
  400 CONTINUE
      ENDIF
C-----------------------------------------------------------------------
CL  4. Smooth bottom row (excluding the two end points)
C-----------------------------------------------------------------------
      IF (atbase) THEN
      I=P_ROWS-Offy
      DO 500 J=2,ROW_LENGTH-1
        POINT=(I-1)*ROW_LENGTH+J
        ABAR(POINT)=(IW1*A(POINT-ROW_LENGTH-1)+IW2*A(POINT-ROW_LENGTH)
     1   +IW1*A(POINT-ROW_LENGTH+1)+IW2*A(POINT-1)+IW3*A(POINT)
     2   +IW2*A(POINT+1))/IW6
  500 CONTINUE
      ENDIF
C-----------------------------------------------------------------------
CL  5. Smooth four corner grid points
C-----------------------------------------------------------------------
      IW5=IW3+2*IW2+IW1
      IF (attop .AND. atleft) THEN
      POINT=Offy*ROW_LENGTH+Offx+1
      ABAR(POINT)=(IW3*A(POINT)+IW2*A(POINT+1)+IW2*A(POINT+ROW_LENGTH)
     1 +IW1*A(POINT+ROW_LENGTH+1))/IW5
      ENDIF

      IF (attop .AND. atright) THEN
      POINT=(Offy+1)*ROW_LENGTH-Offx
      ABAR(POINT)=(IW2*A(POINT-1)+IW3*A(POINT)
     1 +IW1*A(POINT+ROW_LENGTH-1)+IW2*A(POINT+ROW_LENGTH))/IW5
      ENDIF

      IF (atbase .AND. atleft) THEN
      POINT=(P_ROWS-Offy-1)*ROW_LENGTH+Offx+1
      ABAR(POINT)=(IW2*A(POINT-ROW_LENGTH)+IW1*A(POINT-ROW_LENGTH+1)
     1 +IW3*A(POINT)+IW2*A(POINT+1))/IW5
      ENDIF

      IF (atbase .AND. atright) THEN
      POINT=(P_ROWS-Offy)*ROW_LENGTH-Offx
      ABAR(POINT)=(IW1*A(POINT-ROW_LENGTH-1)+IW2*A(POINT-ROW_LENGTH)
     1 +IW2*A(POINT-1)+IW3*A(POINT))/IW5
      ENDIF
C-----------------------------------------------------------------------
CL  6. Smooth left hand column
C-----------------------------------------------------------------------
      IF (atleft) THEN

      J=Offx+1
      IF (attop) THEN
        I_start=Offy+2
      ELSE
        I_start=Offy+1
      ENDIF

      IF (atbase) THEN
        I_end=P_ROWS-Offy-1
      ELSE
        I_end=P_ROWS-Offy
      ENDIF

      DO 600 I=I_start,I_end
        POINT=(I-1)*ROW_LENGTH+J
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
  600 CONTINUE
      ENDIF
C-----------------------------------------------------------------------
CL  7. Smooth right hand column
C-----------------------------------------------------------------------
      IF (atright) THEN

      J=ROW_LENGTH-Offx
      IF (attop) THEN
        I_start=Offy+2
      ELSE
        I_start=Offy+1
      ENDIF

      IF (atbase) THEN
        I_end=P_ROWS-Offy-1
      ELSE
        I_end=P_ROWS-Offy
      ENDIF

      DO 700 I=I_start,I_end
        POINT=(I-1)*ROW_LENGTH+J
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
  700 CONTINUE
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE LUBATT
     1 (WINDSPEED,T_DRY,T_SEA,T_DEW,DUCT_HEIGHT,
     2  MAX_WAVELENGTH)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      REAL
     *  WINDSPEED       ! IN  wind speed in m/s
     *, T_DRY           ! IN  dry bulb temperature in deg C
     *, T_SEA           ! IN  sea surface temperature in deg C
     *, T_DEW           ! IN  dew point temperature in deg C
     *, DUCT_HEIGHT     ! OUT evaporation duct height in metres
     *, MAX_WAVELENGTH  ! OUT Maximum wavelength propagated by duct
     *                  !     (duct intensity) in metres
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
      INTEGER
     *  VSEA(6)  ! values of sea temperature in table
     *, VWIND(6) ! values of windspeed in table
     *, VDPD(6)  ! values of sea temp depression in table
     *, VSTD(6)  ! values of dewpoint depression in table
     *, IH(6,6,6,6)
     *, IW(6,6,6,6)
     *, H(16)    ! values of IH for upper and lower points of four
     *           !    dimensional space formed by intervals
     *, W(16)    ! As H but for IW
C-----------------------------------------------------------------------
      REAL
     *  H2(8)    ! Heights found by linear interpolation for change in
     *           !                                  sea temperature
     *, H3(4)    ! "       "    "    "    "  "  " " wind speed
     *, H4(2)    ! "       "    "    "    "  "  " " dewpoint depression
     *, W2(8)    ! Wavelength found by linear interpolation for change
     *           !                               in sea temperature
     *, W3(4)    ! "       "    "    "    "  "  " " wind speed
     *, W4(2)    ! "       "    "    "    "  "  " " dewpoint depression
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
C     NONE
C*----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C   LOCAL CONSTANTS
C-----------------------------------------------------------------------
      INTEGER
     *  NSEA  ! No of values of sea temperature in table
     *, NWIN  ! No of values of windspeed in table
     *, NDPD  ! No of values of dewpoint depression in table
     *, NSTD  ! No of values of sea temp depression in table
      PARAMETER(NSEA=6,NWIN=6,NDPD=6,NSTD=6)
C-----------------------------------------------------------------------
C   LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  NS       ! lower boundary of sea temperature interval
     *, NW       ! lower boundary of wind interval
     *, NT       ! lower boundary of sea temp depression interval
     *, ND       ! lower boundary of dewpoint depression interval
     *, I,J,K    ! Loop counters for DATA statements
      REAL
     *  STD      ! sea temp depression  = T_DRY - T_SEA
     *, DPD      ! dewpoint depression  = T_DRY - T_DEW
     *, RS       ! =(T_SEA-VSEA(NS))/(VSEA(NS+1)-VSEA(NS))
     *, RW       ! =(WINDSPEED-VWIND(NW))/(VWIND(NW+1)-VWIND(NW))
     *, RT       ! =(STD-VSTD(NT))/(VSTD(NT+1)-VSTD(NT))
     *, RD       ! =(DPD-VDPD(ND))/(VDPD(ND+1)-VDPD(ND))
C-----------------------------------------------------------------------
CL  1. DATA statements for six values of sea temp, wind speed,
CL     dewpoint depression and sea temp depression
C-----------------------------------------------------------------------
      DATA VSEA/-2,8,16,22,26,28/
      DATA VWIND/1,3,5,8,14,20/
      DATA VDPD/0,2,4,6,8,10/
      DATA VSTD/-12,-8,-4,0,2,4/
C-----------------------------------------------------------------------
CL  2. DATA initialisation for IH(6,6,6,6) array
C-----------------------------------------------------------------------
C     first parameter 6 sea temperature depression values.
C     second parameter 6 dewpoint depression values.
C     third parameter 6 wind speed values.
C     fourth parameter 6 sea temperature values.
C
C T_SEA -2
C WINDSPEED 1,3,5
C
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,1),I=1,6),J=1,6),K=1,3)/
     *         1,1,1,1,0,0,
     *         1,1,1,7,30,30,
     *         1,1,1,10,30,30,
     *         1,1,2,4,30,30,
     *         1,1,2,5,30,30,
     *         1,1,2,5,30,30,
C
     *         1,1,1,1,0,0,
     *         1,1,2,2,3,1,
     *         1,2,2,3,30,30,
     *         1,2,2,5,30,30,
     *         1,2,3,6,30,30,
     *         1,1,3,6,23,30,
C
     *         1,1,1,1,0,0,
     *         1,2,2,2,1,1,
     *         2,2,2,3,5,30,
     *         1,2,3,5,30,30,
     *         1,3,3,6,30,30,
     *         1,2,4,7,30,30/
C-----------------------------------------------------------------------
C T_SEA -2
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,1),I=1,6),J=1,6),K=4,6)/
     *         1,1,1,1,0,0,
     *         2,2,2,2,1,1,
     *         2,3,3,3,3,3,
     *         2,3,4,5,7,12,
     *         1,4,4,6,12,30,
     *         1,3,5,7,20,30,
C
     *         1,1,1,1,0,0,
     *         2,3,3,2,1,1,
     *         3,3,4,4,3,3,
     *         3,4,5,5,6,6,
     *         3,5,6,7,8,9,
     *         2,5,6,8,10,13,
C
     *         1,1,1,1,0,0,
     *         2,3,3,2,1,1,
     *         3,4,4,4,4,3,
     *         4,5,5,5,6,5,
     *         4,5,6,7,7,8,
     *         3,6,7,8,9,10/
C-----------------------------------------------------------------------
C T_SEA 8
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,2),I=1,6),J=1,6),K=1,3)/
     *         1,1,1,1,0,0,
     *         2,2,2,3,30,0,
     *         2,2,2,5,30,30,
     *         2,2,3,5,30,30,
     *         2,2,3,5,30,30,
     *         2,3,3,5,30,30,
C
     *         2,2,2,1,0,0,
     *         3,3,3,3,3,0,
     *         3,3,3,6,30,30,
     *         4,4,4,7,30,30,
     *         4,4,4,7,30,30,
     *         4,4,5,8,30,30,
C
     *         2,3,2,1,0,0,
     *         4,4,4,3,1,0,
     *         4,4,4,6,30,30,
     *         5,5,5,7,30,30,
     *         5,5,6,8,30,30,
     *         5,6,6,10,30,30/
C-----------------------------------------------------------------------
C T_SEA 8
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,2),I=1,6),J=1,6),K=4,6)/
     *         3,4,3,1,0,0,
     *         5,5,4,3,1,0,
     *         6,6,6,6,6,3,
     *         6,6,7,8,18,30,
     *         6,7,7,10,30,30,
     *         7,7,8,11,30,30,
C
     *         4,5,4,1,0,0,
     *         7,7,6,3,1,0,
     *         8,8,7,6,5,3,
     *         9,9,9,9,9,8,
     *         9,10,10,11,14,17,
     *         10,10,11,13,18,30,
C
     *         4,5,4,1,0,0,
     *         9,8,7,4,1,0,
     *         10,9,9,7,5,3,
     *         10,11,11,9,9,7,
     *         11,12,12,12,13,12,
     *         12,13,13,14,16,17/
C-----------------------------------------------------------------------
C T_SEA 16
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,3),I=1,6),J=1,6),K=1,3)/
     *         2,2,2,1,0,0,
     *         3,3,3,4,30,0,
     *         3,3,3,4,30,30,
     *         3,3,3,5,30,30,
     *         3,3,4,5,30,30,
     *         3,4,4,6,30,30,
C
     *         4,4,3,1,0,0,
     *         5,4,4,4,2,0,
     *         5,5,5,6,30,30,
     *         5,5,5,8,30,30,
     *         6,6,6,8,30,30,
     *         6,6,6,9,30,30,
C
     *         5,5,4,1,0,0,
     *         6,6,5,5,1,0,
     *         7,6,6,7,30,30,
     *         7,7,7,9,30,30,
     *         7,8,8,11,30,30,
     *         8,8,9,12,30,30/
C-----------------------------------------------------------------------
C T_SEA 16
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,3),I=1,6),J=1,6),K=4,6)/
     *         5,6,5,1,0,0,
     *         8,8,7,5,1,0,
     *         9,8,8,8,11,3,
     *         9,9,9,11,30,30,
     *         10,10,10,13,30,30,
     *         10,11,11,15,30,30,
C
     *         8,9,6,1,0,0,
     *         11,10,9,5,1,0,
     *         12,12,11,9,7,3,
     *         13,13,13,13,15,13,
     *         14,14,15,16,23,30,
     *         15,15,16,18,30,30,
C
     *         9,10,7,1,0,0,
     *         14,13,10,5,1,0,
     *         15,15,13,10,7,3,
     *         17,16,16,14,13,10,
     *         18,18,18,17,19,19,
     *         19,19,19,20,25,29/
C-----------------------------------------------------------------------
C T_SEA 22
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,4),I=1,6),J=1,6),K=1,3)/
     *         3,3,2,1,0,0,
     *         3,3,3,4,30,0,
     *         4,4,4,5,30,30,
     *         4,4,4,5,30,30,
     *         4,4,5,6,30,30,
     *         4,4,5,6,30,30,
C
     *         5,5,4,1,0,0,
     *         6,5,5,5,3,0,
     *         6,6,6,7,30,30,
     *         7,7,7,8,30,30,
     *         7,7,7,9,30,30,
     *         8,8,8,10,30,30,
C
     *         7,6,5,1,0,0,
     *         8,7,6,6,1,0,
     *         8,8,8,9,30,30,
     *         9,9,9,11,30,30,
     *         9,9,10,12,30,30,
     *         10,10,10,13,30,30/
C-----------------------------------------------------------------------
C T_SEA 22
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,4),I=1,6),J=1,6),K=4,6)/
     *         9,8,6,1,0,0,
     *         10,10,8,6,1,0,
     *         11,11,10,10,19,4,
     *         12,12,12,13,30,30,
     *         13,13,13,15,30,30,
     *         13,13,14,17,30,30,
C
     *         11,11,8,1,0,0,
     *         15,14,11,6,1,0,
     *         16,15,14,12,10,3,
     *         17,17,16,17,20,18,
     *         18,18,18,20,30,30,
     *         19,19,20,22,30,30,
C
     *         14,14,9,1,0,0,
     *         19,17,14,7,2,0,
     *         20,19,17,13,9,3,
     *         22,21,20,18,18,13,
     *         23,23,23,22,26,27,
     *         24,24,25,26,30,30/
C-----------------------------------------------------------------------
C T_SEA 26
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,5),I=1,6),J=1,6),K=1,3)/
     *         4,3,3,1,0,0,
     *         4,4,4,4,30,0,
     *         4,4,4,5,30,30,
     *         4,4,5,6,30,30,
     *         5,5,5,6,30,30,
     *         5,5,6,7,15,30,
C
     *         6,5,4,1,0,0,
     *         7,6,5,6,3,0,
     *         7,7,7,8,30,30,
     *         8,8,7,9,30,30,
     *         8,8,8,10,30,30,
     *         9,9,9,11,22,30,
C
     *         8,7,5,1,0,0,
     *         9,8,7,6,1,0,
     *         10,9,9,10,30,30,
     *         10,10,10,12,30,30,
     *         11,11,11,13,30,30,
     *         11,11,12,14,25,30/
C-----------------------------------------------------------------------
C T_SEA 26
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,5),I=1,6),J=1,6),K=4,6)/
     *         11,9,7,1,0,0,
     *         12,11,10,7,1,0,
     *         13,12,12,12,27,4,
     *         14,14,13,15,30,30,
     *         15,15,15,17,30,30,
     *         15,15,16,19,25,30,
C
     *         14,13,9,1,0,0,
     *         17,16,13,8,1,0,
     *         19,18,16,14,12,4,
     *         20,20,19,19,25,25,
     *         21,21,21,22,30,30,
     *         22,22,23,25,27,30,
C
     *         18,17,11,1,0,0,
     *         22,20,16,8,2,0,
     *         24,23,20,15,11,4,
     *         26,25,23,21,21,16,
     *         27,27,26,26,30,30,
     *         28,28,29,30,28,30/
C-----------------------------------------------------------------------
C T_SEA 28
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,6),I=1,6),J=1,6),K=1,3)/
     *         4,3,3,1,0,0,
     *         4,4,4,4,30,0,
     *         4,4,4,5,30,30,
     *         5,5,5,6,30,30,
     *         5,5,5,6,13,30,
     *         5,5,6,7,16,30,
C
     *         7,6,4,1,0,0,
     *         7,7,6,6,3,0,
     *         8,8,7,8,30,30,
     *         9,8,8,9,30,30,
     *         9,9,9,10,20,30,
     *         9,9,9,11,24,30,
C
     *         9,8,6,1,0,0,
     *         10,9,8,6,1,0,
     *         10,10,9,10,30,30,
     *         11,11,11,12,30,30,
     *         12,12,12,13,22,30,
     *         12,12,12,15,27,30/
C-----------------------------------------------------------------------
C T_SEA 28
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IH(I,J,K,6),I=1,6),J=1,6),K=4,6)/
     *         12,10,7,1,0,0,
     *         13,12,10,7,1,0,
     *         14,13,12,12,30,4,
     *         15,14,14,15,30,30,
     *         16,16,16,18,30,30,
     *         16,16,17,20,27,30,
C
     *         17,14,10,1,0,0,
     *         19,17,14,8,1,0,
     *         20,19,18,15,13,4,
     *         22,21,20,20,27,30,
     *         23,22,22,24,30,30,
     *         24,24,24,27,29,30,
C
     *         20,18,12,1,0,0,
     *         24,21,17,9,2,0,
     *         26,24,22,16,12,4,
     *         28,27,25,23,23,18,
     *         29,29,28,28,30,30,
     *         30,30,30,30,30,30/
C-----------------------------------------------------------------------
CL  3. DATA initialisation for IW(6,6,6,6) array
C-----------------------------------------------------------------------
C T_SEA -2
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,1),I=1,6),J=1,6),K=1,3)/
     *         0,0,0,0,0,0,
     *         0,1,1,2,41,47,
     *         1,1,1,4,51,73,
     *         1,1,2,3,51,77,
     *         1,1,2,4,51,77,
     *         1,1,2,5,50,76,
C
     *         0,0,0,0,0,0,
     *         1,1,1,1,1,0,
     *         1,1,1,2,21,26,
     *         1,1,2,4,28,37,
     *         1,2,3,5,32,41,
     *         1,1,3,5,25,41,
C
     *         0,0,0,0,0,0,
     *         1,1,1,1,0,0,
     *         1,1,1,2,2,13,
     *         1,1,2,4,17,22,
     *         1,2,2,5,22,27,
     *         1,2,4,6,25,31/
C-----------------------------------------------------------------------
C T_SEA -2
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,1),I=1,6),J=1,6),K=4,6)/
     *         0,0,0,0,0,0,
     *         1,1,1,1,0,0,
     *         1,2,2,2,2,1,
     *         1,2,3,4,4,6,
     *         1,3,3,5,9,21,
     *         1,2,4,6,16,25,
C
     *         0,0,0,0,0,0,
     *         1,2,2,1,0,0,
     *         2,2,3,2,2,1,
     *         2,3,4,4,4,4,
     *         2,4,5,6,6,7,
     *         1,4,5,7,8,11,
C
     *         0,0,0,0,0,0,
     *         1,2,2,1,0,0,
     *         2,2,2,2,2,1,
     *         2,3,4,3,4,3,
     *         3,4,5,6,5,6,
     *         2,5,6,7,8,8/
C-----------------------------------------------------------------------
C T_SEA 8
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,2),I=1,6),J=1,6),K=1,3)/
     *         1,1,1,0,0,0,
     *         2,2,2,5,41,0,
     *         2,2,2,4,62,73,
     *         2,2,3,5,61,83,
     *         2,2,4,6,60,82,
     *         2,4,4,6,60,81,
C
     *         2,2,1,0,0,0,
     *         3,3,3,2,1,0,
     *         3,3,3,5,27,27,
     *         5,5,5,7,35,42,
     *         5,5,5,8,41,52,
     *         5,5,6,10,45,54,
C
     *         2,3,1,0,0,0,
     *         4,4,3,2,0,0,
     *         4,4,4,5,17,13,
     *         6,6,6,7,25,26,
     *         6,6,7,9,31,33,
     *         6,8,8,12,35,38/
C-----------------------------------------------------------------------
C T_SEA 8
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,2),I=1,6),J=1,6),K=4,6)/
     *         3,3,2,0,0,0,
     *         5,5,3,2,0,0,
     *         7,6,6,5,4,1,
     *         7,7,8,8,15,21,
     *         7,9,8,11,30,28,
     *         9,9,10,13,34,33,
C
     *         4,4,2,0,0,0,
     *         7,7,5,2,0,0,
     *         9,8,7,5,3,1,
     *         10,10,10,9,8,6,
     *         11,12,12,12,14,16,
     *         13,13,14,15,20,32,
C
     *         4,4,2,0,0,0,
     *         9,7,6,2,0,0,
     *         11,9,9,6,3,2,
     *         11,12,12,9,8,5,
     *         13,14,14,13,13,11,
     *         15,16,16,17,18,18/
C-----------------------------------------------------------------------
C T_SEA 16
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,3),I=1,6),J=1,6),K=1,3)/
     *         3,2,2,0,0,0,
     *         4,4,3,3,41,0,
     *         4,4,4,4,73,75,
     *         5,4,4,6,72,90,
     *         5,5,6,7,71,89,
     *         5,7,7,9,70,89,
C
     *         5,5,3,0,0,0,
     *         7,5,5,3,1,0,
     *         8,7,7,6,32,28,
     *         8,8,7,10,41,48,
     *         10,10,9,11,47,60,
     *         11,10,10,14,51,68,
C
     *         7,6,4,0,0,0,
     *         9,8,6,4,0,0,
     *         11,9,8,7,22,15,
     *         11,11,10,11,32,31,
     *         12,13,12,15,38,39,
     *         14,14,15,18,43,45/
C-----------------------------------------------------------------------
C T_SEA 16
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,3),I=1,6),J=1,6),K=4,6)/
     *         7,7,4,0,0,0,
     *         12,10,8,4,0,0,
     *         14,11,10,8,8,1,
     *         14,14,13,13,31,26,
     *         17,16,15,18,38,34,
     *         17,19,18,23,43,40,
C
     *         11,10,5,0,0,0,
     *         16,13,10,4,0,0,
     *         18,17,14,9,6,2,
     *         21,20,18,16,16,11,
     *         23,23,23,22,29,34,
     *         26,25,26,27,42,40,
C
     *         12,11,6,0,0,0,
     *         20,17,11,4,0,0,
     *         23,21,16,10,6,2,
     *         27,34,22,17,14,9,
     *         30,29,27,23,24,21,
     *         33,32,31,30,35,38/
C-----------------------------------------------------------------------
C T_SEA 22
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,4),I=1,6),J=1,6),K=1,3)/
     *         5,4,2,0,0,0,
     *         5,5,4,3,42,0,
     *         7,7,6,6,84,78,
     *         7,7,7,7,83,99,
     *         8,7,9,10,82,97,
     *         8,8,10,11,81,96,
C
     *         8,7,4,0,0,0,
     *         11,8,7,4,1,0,
     *         11,10,9,8,36,29,
     *         14,13,12,12,46,53,
     *         14,14,13,15,51,66,
     *         17,16,16,18,53,75,
C
     *         12,8,5,0,0,0,
     *         14,11,8,5,0,0,
     *         15,14,12,11,27,16,
     *         18,17,15,16,37,35,
     *         18,18,18,19,45,44,
     *         21,20,19,23,50,50/
C-----------------------------------------------------------------------
C T_SEA 22
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,4),I=1,6),J=1,6),K=4,6)/
     *         15,11,6,0,0,0,
     *         18,16,10,5,0,0,
     *         21,19,15,12,17,2,
     *         24,22,20,19,37,30,
     *         26,25,24,24,44,40,
     *         27,26,27,30,50,47,
C
     *         18,15,8,0,0,0,
     *         26,22,14,5,0,0,
     *         30,26,21,14,9,2,
     *         33,31,27,24,24,18,
     *         36,35,32,32,44,39,
     *         40,38,38,38,50,46,
C
     *         23,19,9,0,0,0,
     *         33,27,18,6,1,0,
     *         37,32,25,15,8,2,
     *         43,38,33,25,22,13,
     *         46,44,41,35,38,35,
     *         50,48,48,45,50,46/
C-----------------------------------------------------------------------
C T_SEA 26
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,5),I=1,6),J=1,6),K=1,3)/
     *         7,5,4,0,0,0,
     *         7,7,6,4,43,0,
     *         8,7,7,7,92,80,
     *         8,8,9,10,91,105,
     *         11,10,10,11,89,104,
     *         11,11,13,14,27,102,
C
     *         11,8,5,0,0,0,
     *         14,11,7,6,1,0,
     *         15,14,12,11,38,30,
     *         18,17,13,14,48,57,
     *         18,18,16,18,52,70,
     *         21,21,19,22,40,79,
C
     *         15,11,6,0,0,0,
     *         18,14,10,6,0,0,
     *         21,17,15,13,30,17,
     *         22,21,19,19,41,37,
     *         25,24,22,23,49,47,
     *         26,25,26,27,46,54/
C-----------------------------------------------------------------------
C T_SEA 26
C WINDSPEED 8,14,20
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,5),I=1,6),J=1,6),K=4,6)/
     *         20,14,8,0,0,0,
     *         24,20,15,7,0,0,
     *         27,23,20,15,26,2,
     *         31,29,24,24,41,33,
     *         34,33,30,30,49,44,
     *         35,34,34,37,45,51,
C
     *         26,20,10,0,0,0,
     *         34,28,19,7,0,0,
     *         40,35,27,18,12,2,
     *         44,41,35,30,34,27,
     *         48,45,42,39,49,43,
     *         52,50,49,49,49,51,
C
     *         33,26,12,0,0,0,
     *         43,35,23,7,1,0,
     *         50,44,33,19,11,2,
     *         57,51,43,32,29,17,
     *         61,58,52,46,49,43,
     *         66,63,62,59,51,50/
C-----------------------------------------------------------------------
C T_SEA 28
C WINDSPEED 1,3,5
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,6),I=1,6),J=1,6),K=1,3)/
     *         7,5,4,0,0,0,
     *         8,7,6,4,43,0,
     *         8,8,7,7,96,81,
     *         11,10,10,10,95,109,
     *         11,11,10,11,22,107,
     *         11,11,13,15,31,106,
C
     *         14,10,5,0,0,0,
     *         15,13,9,6,1,0,
     *         18,16,13,11,40,31,
     *         21,18,16,15,49,59,
     *         22,21,19,19,34,72,
     *         23,22,20,23,46,80,
C
     *         18,13,7,0,0,0,
     *         21,17,12,6,0,0,
     *         22,21,16,14,31,18,
     *         26,24,22,20,43,39,
     *         29,28,26,24,37,49,
     *         30,29,27,31,52,55/
C-----------------------------------------------------------------------
C T_SEA 28
C WINDSPEED 8,14,20,
C-----------------------------------------------------------------------
      DATA (((IW(I,J,K,6),I=1,6),J=1,6),K=4,6)/
     *         23,17,9,0,0,0,
     *         27,22,15,7,0,0,
     *         31,27,21,16,31,2,
     *         35,31,27,25,43,34,
     *         39,37,34,34,51,46,
     *         40,38,38,41,52,54,
C
     *         33,23,12,0,0,0,
     *         40,32,21,8,0,0,
     *         44,39,32,20,14,2,
     *         51,46,39,32,38,34,
     *         55,50,47,45,51,45,
     *         60,57,54,55,55,53,
C
     *         39,30,14,0,0,0,
     *         50,39,26,9,1,0,
     *         58,49,39,22,13,2,
     *         65,58,49,37,33,20,
     *         70,66,59,51,51,45,
     *         75,72,67,61,57,53/
C-----------------------------------------------------------------------
CL  4. Calculate sea temp depression and dewpoint depression
C-----------------------------------------------------------------------
      STD = T_DRY - T_SEA
      DPD = T_DRY - T_DEW
C-----------------------------------------------------------------------
CL  5. Check that sea temperature lies in acceptable range
CL     Find lower boundary value of sea temperature interval
CL     for interpolation
C-----------------------------------------------------------------------
      IF(T_SEA.LT.VSEA(1)) T_SEA = VSEA(1)
      IF(T_SEA.GT.VSEA(NSEA)) T_SEA = VSEA(NSEA)
      NS = 0
  100 NS = NS + 1
      IF(T_SEA.GT.VSEA(NS+1)) GO TO 100
C-----------------------------------------------------------------------
CL  6. Check that wind speed lies in acceptable range
CL     Find lower boundary value of wind speed interval
CL     for interpolation
C-----------------------------------------------------------------------
      IF(WINDSPEED.LT.VWIND(1)) WINDSPEED = VWIND(1)
      IF(WINDSPEED.GT.VWIND(NWIN)) WINDSPEED = VWIND(NWIN)
      NW = 0
  120 NW = NW + 1
      IF(WINDSPEED.GT.VWIND(NW+1)) GO TO 120
C-----------------------------------------------------------------------
CL  7. Check that sea temperature depression lies in acceptable range
CL     Find lower boundary value of sea temperature depression
CL     for interpolation
C-----------------------------------------------------------------------
      IF(STD.LT.VSTD(1)) STD = VSTD(1)
      IF(STD.GT.VSTD(NSTD)) THEN
        STD=VSTD(NSTD)
      ENDIF
      NT = 0
  140 NT = NT + 1
      IF(STD.GT.VSTD(NT+1)) GO TO 140
C-----------------------------------------------------------------------
CL  8. Check that dew point depression lies in acceptable range
CL     Find lower boundary value of dew point depression
CL     for interpolation
C-----------------------------------------------------------------------
      IF(DPD.LT.VDPD(1)) DPD = VDPD(1)
      IF(DPD.GT.VDPD(NDPD)) DPD = VDPD(NDPD)
      ND = 0
  160 ND = ND + 1
      IF(DPD.GT.VDPD(ND+1)) GO TO 160
C-----------------------------------------------------------------------
CL  9. Four dimensional interpolation of duct height
C-----------------------------------------------------------------------
      H(1) = IH(NT,ND,NW,NS)
      H(2) = IH(NT,ND,NW,NS+1)
C
      H(3) = IH(NT,ND,NW+1,NS)
      H(4) = IH(NT,ND,NW+1,NS+1)
C
      H(5) = IH(NT,ND+1,NW,NS)
      H(6) = IH(NT,ND+1,NW,NS+1)
      H(7) = IH(NT,ND+1,NW+1,NS)
      H(8) = IH(NT,ND+1,NW+1,NS+1)
C
      H(9) = IH(NT+1,ND,NW,NS)
      H(10) =IH(NT+1,ND,NW,NS+1)
      H(11) =IH(NT+1,ND,NW+1,NS)
      H(12) =IH(NT+1,ND,NW+1,NS+1)
      H(13) =IH(NT+1,ND+1,NW,NS)
      H(14) =IH(NT+1,ND+1,NW,NS+1)
      H(15) =IH(NT+1,ND+1,NW+1,NS)
      H(16) =IH(NT+1,ND+1,NW+1,NS+1)
C-----------------------------------------------------------------------
C INTERPOLATE HEIGHT FOR CHANGE IN SEA TEMPERATURE
C-----------------------------------------------------------------------
      RS = (T_SEA - VSEA(NS))/(VSEA(NS+1) - VSEA(NS))
      H2(1) = H(1) + (H(2) - H(1)) * RS
      H2(2) = H(3) + (H(4) - H(3)) * RS
      H2(3) = H(5) + (H(6) - H(5)) * RS
      H2(4) = H(7) + (H(8) - H(7)) * RS
      H2(5) = H(9) + (H(10) - H(9)) * RS
      H2(6) = H(11) + (H(12) - H(11)) * RS
      H2(7) = H(13) + (H(14) - H(13)) * RS
      H2(8) = H(15) + (H(16) - H(15)) * RS
C-----------------------------------------------------------------------
C INTERPOLATE HEIGHT FOR CHANGE IN WIND SPEED
C-----------------------------------------------------------------------
      RW = (WINDSPEED - VWIND(NW))/(VWIND(NW+1) - VWIND(NW))
      H3(1) = H2(1) + (H2(2) - H2(1)) * RW
      H3(2) = H2(3) + (H2(4) - H2(3)) * RW
      H3(3) = H2(5) + (H2(6) - H2(5)) * RW
      H3(4) = H2(7) + (H2(8) - H2(7)) * RW
C-----------------------------------------------------------------------
C INTERPOLATE HEIGHT FOR CHANGE IN DEWPOINT DEPRESSION
C-----------------------------------------------------------------------
      RD = (DPD - VDPD(ND))/(VDPD(ND+1) - VDPD(ND))
      H4(1) = H3(1) + (H3(2) - H3(1)) * RD
      H4(2) = H3(3) + (H3(4) - H3(3)) * RD
C-----------------------------------------------------------------------
C INTERPOLATE HEIGHT FOR CHANGE IN SEA TEMPERATURE DEPRESSION
C-----------------------------------------------------------------------
      RT = (STD - VSTD(NT))/(VSTD(NT+1) - VSTD(NT))
      DUCT_HEIGHT = H4(1) + (H4(2) - H4(1)) * RT
C-----------------------------------------------------------------------
CL 10. Four dimensional interpolation of maximum wavelength
C-----------------------------------------------------------------------
      W(1) = IW(NT,ND,NW,NS)
      W(2) = IW(NT,ND,NW,NS+1)
C
      W(3) = IW(NT,ND,NW+1,NS)
      W(4) = IW(NT,ND,NW+1,NS+1)
C
      W(5) = IW(NT,ND+1,NW,NS)
      W(6) = IW(NT,ND+1,NW,NS+1)
      W(7) = IW(NT,ND+1,NW+1,NS)
      W(8) = IW(NT,ND+1,NW+1,NS+1)
C
      W(9) = IW(NT+1,ND,NW,NS)
      W(10) = IW(NT+1,ND,NW,NS+1)
      W(11) = IW(NT+1,ND,NW+1,NS)
      W(12) = IW(NT+1,ND,NW+1,NS+1)
      W(13) = IW(NT+1,ND+1,NW,NS)
      W(14) = IW(NT+1,ND+1,NW,NS+1)
      W(15) = IW(NT+1,ND+1,NW+1,NS)
      W(16) = IW(NT+1,ND+1,NW+1,NS+1)
C-----------------------------------------------------------------------
C INTERPOLATE WAVELENGTH CHANGE IN SEA TEMPERATURE
C-----------------------------------------------------------------------
      W2(1) = W(1) + (W(2) - W(1)) * RS
      W2(2) = W(3) + (W(4) - W(3)) * RS
      W2(3) = W(5) + (W(6) - W(5)) * RS
      W2(4) = W(7) + (W(8) - W(7)) * RS
      W2(5) = W(9) + (W(10) - W(9)) * RS
      W2(6) = W(11) + (W(12) - W(11)) * RS
      W2(7) = W(13) + (W(14) - W(13)) * RS
      W2(8) = W(15) + (W(16) - W(15)) * RS
C-----------------------------------------------------------------------
C INTERPOLATE WAVELENGTH FOR CHANGE IN WIND SPEED
C-----------------------------------------------------------------------
      W3(1) = W2(1) + (W2(2) - W2(1)) * RW
      W3(2) = W2(3) + (W2(4) - W2(3)) * RW
      W3(3) = W2(5) + (W2(6) - W2(5)) * RW
      W3(4) = W2(7) + (W2(8) - W2(7)) * RW
C-----------------------------------------------------------------------
C INTERPOLATE WAVELENGTH FOR CHANGE IN DEWPOINT DEPRESSION
C-----------------------------------------------------------------------
      W4(1) = W3(1) + (W3(2) - W3(1)) * RD
      W4(2) = W3(3) + (W3(4) - W3(3)) * RD
C-----------------------------------------------------------------------
C INTERPOLATE WAVELENGTH FOR CHANGE IN SEA TEMPERATURE DEPRESSION
C-----------------------------------------------------------------------
      MAX_WAVELENGTH = (W4(1) + (W4(2) - W4(1)) * RT)*0.01
      RETURN
      END
C-----------------------------------------------------------------------
C
