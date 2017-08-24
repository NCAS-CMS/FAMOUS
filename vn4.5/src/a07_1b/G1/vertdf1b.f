C ******************************COPYRIGHT******************************
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
C
CLL  SUBROUTINE VDIF_CTL and VERT_DIF --------------------------
CLL
CLL  PURPOSE:   CONTROL SECTION FOR VERTICAL DIFFUSION ROUTINE WHICH
CLL APPLIES VERTICAL DIFFUSION TO WIND COMPONENTS WITHIN A LATITUDE BAND
CLL SYMMETRIC ABOUT THE EQUATOR.  THE DIFFUSION COEFFICIENT TAILS OFF
CLL AWAY FROM THE EQUATOR AND IS ZERO OUTSIDE THE BAND.
CLL
CLL  FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY
CLL  SUITABLE FOR SINGLE COLUMN USE, CALL TO P_TO_UV BY-PASSED
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL  ORIGINAL VERSION FOR CRAY Y-MP
CLL
CLL  WRITTEN BY C. WILSON
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with VERTDF1A.
CLL   4.4    14/07/97 Gather to diffusion points removed, since it
CLL                   it becomes an overhead on T3E for equatorial
CLL                   PEs.
CLL                   A. Dickinson
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  Programming standard:
CLL
CLL  Logical components covered: P21
CLL
CLL  Project task:
CLL
CLL  Documentation:        The equation used is (2)
CLL                        in Unified Model documentation paper no.p21
CLL                        C. Wilson, version 2,dated 30/10/89
CLLEND-------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE VDIF_CTL
     *  (PSTAR,U,V,
     *   P_FIELD,U_FIELD,ROWS,FIRST_ROW,ROW_LENGTH,
     *   LEVEL_START,LEVEL_END,LEVELS_VD,P_LEVELS,
     *   AK,BK,DELTA_AK,DELTA_BK,COS_LAT, LATITUDE_BAND,
     *   VERTICAL_DIFFUSION, TIMESTEP,
     *   STASH_U_FLUX,FLUX_UD_ON,U_LIST,
     *   STASH_V_FLUX,FLUX_VD_ON,V_LIST,
     *   LEN_STASH_U_FLUX,LEN_STASH_V_FLUX,
     *   POINTS_FLUX_U,POINTS_FLUX_V,LEVELS_FLUX,
     *   IRET)

      IMPLICIT NONE

      INTEGER
     *  P_FIELD            !IN    1ST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD            !IN    1ST DIMENSION OF FIELD OF U,V
     *, ROWS               !IN    NUMBER OF ROWS TO BE UPDATED.
     *, FIRST_ROW          !IN    NUMBER OF FIRST ROW IN DATA ARRAY
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, LEVEL_START        !IN    BOTTOM LEVEL TO BE UPDATED.
     *, LEVEL_END          !IN    TOP    LEVEL TO BE UPDATED.
     *, LEVELS_VD          !IN    NO OF VERTICAL DIFFUSION LEVELS
     *, P_LEVELS           !IN    NUMBER OF MODEL LEVELS
     *, LEN_STASH_U_FLUX   !IN    DIMENSION OF STASH_U_FLUX
     *, LEN_STASH_V_FLUX   !IN    DIMENSION OF STASH_V_FLUX
     *, POINTS_FLUX_U      !IN    NO OF POINTS IN U FLUX FIELD
     *, POINTS_FLUX_V      !IN    NO OF POINTS IN V FLUX FIELD
     *, LEVELS_FLUX        !IN    NO OF FLUX LEVELS
     *, IRET               ! RETURN CODE      :    IRET=0   NORMAL EXIT

      REAL
     * PSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR PSTAR FIELD
     *,U(U_FIELD,P_LEVELS)    !INOUT PRIMARY MODEL ARRAY FOR U FIELD
     *,V(U_FIELD,P_LEVELS)    !INOUT PRIMARY MODEL ARRAY FOR V FIELD
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,
     *,DELTA_AK(P_LEVELS)     !IN    LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)     !IN    LAYER THICKNESS
     *,AK (P_LEVELS)          !IN    VALUE AT LAYER CENTRE
     *,BK (P_LEVELS)          !IN    VALUE AT LAYER CENTRE
     *,COS_LAT(U_FIELD)       !IN    COS(LAT) AT U POINTS
     *,LATITUDE_BAND          !IN    LATITUDE(RADIANS)
     *                        !      WHERE DIFFUSION PRESCRIBED ZERO
     *,VERTICAL_DIFFUSION     !IN    VALUE OF DIFFUSION COEFFICIENT
     *,TIMESTEP               !IN    TIMESTEP
     *,STASH_U_FLUX(LEN_STASH_U_FLUX,*) !U MOMENTUM FLUX - Diagnostic
     *,STASH_V_FLUX(LEN_STASH_V_FLUX,*) !V MOMENTUM FLUX - Diagnostic

C WARNING : Storage is only assigned by the controling routine
C           for the number of levels requested.

      LOGICAL
     * FLUX_UD_ON                !U momentum diagnostic switch
     *,FLUX_VD_ON                !V momentum diagnostic switch
     *,U_LIST(P_LEVELS)          ! List of levels required
     *,V_LIST(P_LEVELS)          ! List of levels required

C*---------------------------------------------------------------------

! Parameters and Common blocks
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
C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS: 4+(LEVEL_END-LEVEL_START+1)*2

C    +(LEVEL_END-LEVEL_START)*2

C   REAL ARRAYS REQUIRED AT FULL FIELD LENGTH
C   1 INTEGER INDEX ARRAY

      REAL
     * PSTAR_UV( ROWS*ROW_LENGTH )   ! INTERPOLATED PSTAR ON UV-GRID
     *,VERT_DIF_LAT(ROWS*ROW_LENGTH) ! LAT. DEPENDENT DIFFUSION*TIMESTEP
     *,FLUX_U_DG(POINTS_FLUX_U,LEVELS_FLUX) !FLUX
     *,FLUX_V_DG(POINTS_FLUX_V,LEVELS_FLUX) !FLUX



C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL P_TO_UV ,VERT_DIF,TIMER
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS+1) * ROWLENGTH
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS      !     NUMBER OF P POINTS NEEDED
     *, ROWS_P        !     NUMBER OF P ROWS   NEEDED
     *, U_POINTS      !     NUMBER OF U POINTS UPDATED
     *, START_P       !     START POSITION FOR P POINTS NEEDED
     *, START_U       !     START POSITION FOR U POINTS UPDATED
     *, POINTS_VD     !     NUMBER OF POINTS NON-ZERO DIFFUSION COEFFS

      REAL
     *  COS_LAT_BAND    ! COS LAT AT WHICH DIFFUSION SET TO ZERO
     *, COEFF           ! LATITUDE-DEPENDENT DIFFUSION * TIMESTEP
C
      INTEGER    K,I,II,IK,! LOOP COUNTERS IN ROUTINE
     *           KOUT_U,KOUT_V
C

C-------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL    1.     INITIALISATION
C--------------------------

        START_P       = (FIRST_ROW-1)*ROW_LENGTH
        START_U       = START_P
      IF (atbase) THEN
        ROWS_P=ROWS
      ELSE
        ROWS_P=ROWS+1
      ENDIF
        P_POINTS      = (ROWS_P)*ROW_LENGTH
      IF (atbase) THEN
        U_POINTS=(ROWS-1)*ROW_LENGTH
      ELSE
        U_POINTS=ROWS*ROW_LENGTH
      ENDIF

C------------------------------------------------------------------
CL    1.1  CALCULATE LATITUDE-DEPENDENT DIFFUSION COEFFICIENTS*TIMESTEP
CL         AND SET TO ZERO WHERE NO DIFFUSION REQUIRED
C------------------------------------------------------------------

        COS_LAT_BAND = COS(LATITUDE_BAND)
        COEFF = VERTICAL_DIFFUSION * TIMESTEP / (1.- COS_LAT_BAND)
C
        DO I=1,U_POINTS
         VERT_DIF_LAT(I) = COEFF * ( COS_LAT(START_U+I) - COS_LAT_BAND)
         IF(VERT_DIF_LAT(I).LT.0.0)VERT_DIF_LAT(I)=0.0
        END DO

C------------------------------------------------------------------
CL    1.2 INTERPOLATE PSTAR TO UV GRID
C------------------------------------------------------------------


      CALL P_TO_UV(PSTAR(START_P+1),PSTAR_UV,P_POINTS,U_POINTS,
     * ROW_LENGTH,ROWS_P)



      LEVELS_VD=LEVEL_END-LEVEL_START+1

C------------------------------------------------------------------
CL  2.  CALL VERT_DIF AND UPDATE WINDS
C------------------------------------------------------------------

      CALL VERT_DIF(PSTAR_UV,
     &        U(START_U+1,LEVEL_START),V(START_U+1,LEVEL_START),
     &              LEVELS_VD,U_POINTS,U_FIELD,
     &              AK(LEVEL_START),BK(LEVEL_START),
     &              DELTA_AK(LEVEL_START),DELTA_BK(LEVEL_START),
     &              VERT_DIF_LAT,FLUX_U_DG,FLUX_V_DG,
     &              POINTS_FLUX_U,POINTS_FLUX_V,LEVELS_FLUX,
     &              FLUX_UD_ON,FLUX_VD_ON)


      IF (FLUX_UD_ON .OR. FLUX_VD_ON) THEN

        KOUT_U=0
        KOUT_V=0

        DO K = LEVEL_START,LEVEL_END-1

          IF (U_LIST(K)) THEN
            KOUT_U=KOUT_U+1
          END IF
          IF (V_LIST(K)) THEN
            KOUT_V=KOUT_V+1
          END IF

          IK = K-LEVEL_START+1

          DO I=1,U_POINTS

          IF (FLUX_UD_ON .AND. U_LIST(K)) THEN
           STASH_U_FLUX(START_U+I,KOUT_U) = FLUX_U_DG(I,IK)
          ENDIF
          IF (FLUX_VD_ON .AND. V_LIST(K)) THEN
           STASH_V_FLUX(START_U+I,KOUT_V) = FLUX_V_DG(I,IK)
          ENDIF

          ENDDO

        ENDDO

      ENDIF

      IRET=0

1000  CONTINUE
      RETURN
      END
CLL  SUBROUTINE VERT_DIF--------------------------------------------
CLL
CLL  PURPOSE:   TO APPLY VERTICAL DIFFUSION TO WIND COMPONENTS
CLL             WITHIN A LATITUDE BAND SYMMETRIC ABOUT THE EQUATOR.
CLL             THE DIFFUSION COEFFICIENT TAILS OFF AWAY FROM THE
CLL             EQUATOR AND IS ZERO OUTSIDE THE BAND.
CLL             THIS ROUTINE APPLIES A PRECALCULATED
CLL             DIFFUSION COEFFICIENT TO ALL POINTS PASSED TO IT
CLL  SUITABLE FOR SINGLE COLUMN USE
CLL  SUITABLE FOR ROTATED GRIDS
CLL  FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY
CLL  ORIGINAL VERSION FOR CRAY Y-MP
CLL
CLL  WRITTEN BY C. WILSON
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARD:
CLL
CLL  LOGICAL COMPONENTS COVERED: P21
CLL
CLL  PROJECT TASK:
CLL
CLL  DOCUMENTATION:        THE EQUATIONS USED ARE (1) TO (4)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO.P21
CLL                        C. WILSON, VERSION 2,DATED 30/10/89
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE VERT_DIF
     *  (PSTAR,U,V,LEVELS_VD,POINTS_VD,UFIELD,AK,BK,DELTA_AK,DELTA_BK,
     *   DIFFUSION_K,FLUX_U_DG,FLUX_V_DG,POINTS_FLUX_U,POINTS_FLUX_V,
     *   LEVELS_FLUX,FLUX_UD_ON,FLUX_VD_ON)

      IMPLICIT NONE

      INTEGER
     *  POINTS_VD          !IN    NUMBER OF POINTS TO BE UPDATED
     *, LEVELS_VD          !IN    NUMBER OF LEVELS TO BE UPDATED
     *, UFIELD             !IN    DIMENSION OF HORIZ FIELD
     *, POINTS_FLUX_U      !IN    NUMBER OF LEVELS TO BE UPDATED
     *, POINTS_FLUX_V      !IN    NUMBER OF LEVELS TO BE UPDATED
     *, LEVELS_FLUX        !IN    NUMBER OF LEVELS TO BE UPDATED

      REAL
     * PSTAR(UFIELD)       !IN    PSTAR FIELD
     *,U(UFIELD,LEVELS_VD) !INOUT ARRAY FOR U FIELD
     *,V(UFIELD,LEVELS_VD) !INOUT ARRAY FOR V FIELD
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,
     *,DELTA_AK(LEVELS_VD)     !IN    LAYER THICKNESS
     *,DELTA_BK(LEVELS_VD)     !IN    LAYER THICKNESS
     *,AK (LEVELS_VD)          !IN    VALUE AT LAYER CENTRE
     *,BK (LEVELS_VD)          !IN    VALUE AT LAYER CENTRE
     *,DIFFUSION_K(UFIELD)     ! LAT. DEPENDENT DIFFUSION*TIMESTEP
     *,FLUX_U_DG(POINTS_FLUX_U,LEVELS_FLUX) ! U MOMENTUM FLUX
     *                                      ! DIAGNOSTIC
     *,FLUX_V_DG(POINTS_FLUX_V,LEVELS_FLUX) ! V MOMENTUM FLUX
     *                                      ! DIAGNOSTIC
      LOGICAL
     * FLUX_UD_ON                !U momentum diagnostic switch
     *,FLUX_VD_ON                !V momentum diagnostic switch

C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS: 4 REAL ARRAYS REQUIRED
C   AT FULL FIELD LENGTH (=POINTS)
C
      REAL
     * FLUX_U(POINTS_VD,2)           ! DOWNWARD FLUXES U-MOMENTUM
     *,FLUX_V(POINTS_VD,2)           ! DOWNWARD FLUXES V-MOMENTUM


C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
C     NONE
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED =POINTS
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      REAL
     *  DEL_AK          ! DIFFERENCE OF AK ACROSS FULL-LEVELS
     *, DEL_BK          ! DIFFERENCE OF BK ACROSS FULL-LEVELS
     *,DELTA_P          ! P(K+1/2) - P(K-1/2)
     *,DELTA_PL         ! P(K+1)   - P(K)
C
      INTEGER    K,I      ! LOOP COUNTERS IN ROUTINE
      INTEGER    KL,KU,KK ! LEVEL COUNTERS IN ROUTINE
C

C-------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
C------------------------------------------------------------------
CL    1. CALCULATE VERTICAL FLUX OF MOMENTUM , EQN(1) DOCUMENTATION
CL       AND UPDATE U,V
C------------------------------------------------------------------

      KL = 1
      KU = 2
      DO I=1,POINTS_VD
       FLUX_U(I,KL) = 0.0
       FLUX_V(I,KL) = 0.0
      END DO

CL    LOOP OVER LEVELS

      DO K = 1,LEVELS_VD-1

CL      1.1  CALCULATE DELTA_P(K) AND DELTA_PL(K)
        DEL_AK=AK(K+1) - AK(K)
        DEL_BK=BK(K+1) - BK(K)

        DO I=1,POINTS_VD
          DELTA_P=DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
          DELTA_PL=DEL_AK+DEL_BK*PSTAR(I)

CL      1.2  COMPUTE FLUX (+VE UP) AND INCREMENT

          FLUX_U(I,KU)=(U(I,K+1) - U(I,K))*DIFFUSION_K(I)/DELTA_PL
          FLUX_V(I,KU)=(V(I,K+1) - V(I,K))*DIFFUSION_K(I)/DELTA_PL

          U(I,K) = U(I,K) + (FLUX_U(I,KU) - FLUX_U(I,KL))/DELTA_P
          V(I,K) = V(I,K) + (FLUX_V(I,KU) - FLUX_V(I,KL))/DELTA_P

        END DO

        IF (FLUX_UD_ON) THEN   !  SF(201,7)
          DO I=1,POINTS_VD
            FLUX_U_DG(I,K)= FLUX_U(I,KU)
          ENDDO
        ENDIF
        IF (FLUX_VD_ON) THEN   !  SF(202,7)
          DO I=1,POINTS_VD
            FLUX_V_DG(I,K)= FLUX_V(I,KU)
          ENDDO
        ENDIF

C       SWAP STORAGE LOCATIONS FOR LOWER AND UPPER FLUXES
        KK = KL
        KL = KU
        KU = KK

      END DO
CL  END LOOP OVER LEVELS

CL    LAST LEVEL
      K=LEVELS_VD
      DO I=1,POINTS_VD
        DELTA_P=DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
        U(I,K) = U(I,K) - FLUX_U(I,KL)/DELTA_P
        V(I,K) = V(I,K) - FLUX_V(I,KL)/DELTA_P
      END DO

      RETURN
      END
