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
CLL   SUBROUTINE TH_ADV -------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES MASS-WEIGHTED INCREMENTS TO THETAL
CLL DUE TO ADVECTION  BY USING EQUATION (35) TO CALCULATE PROVISIONAL
CLL VALUES OF THETAL AT THE NEW TIME-LEVEL, AND THEN RECALCULATING THE
CLL ADVECTION TERMS ON THE RIGHT-HAND SIDE OF (35) USING THESE
CLL PROVISIONAL VALUES. THE FINAL INCREMENTS ARE CALCULATED AS IN
CLL EQUATION (40). THOSE REQUIRING FILTERING ARE FILTERED AND ALL THE
CLL INCREMENTS ARE ADDED ONTO THE FIELDS USING (40).  IF RUNNING A
CLL GLOBAL MODEL POLAR IS CALLED TO UPDATE POLAR VALUES.
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL MM, DR      <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL MPP CODE ADDED BY P.BURTON
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 4.1:
CLL VERSION  DATE
CLL
CLL 4.1      07/12/95 New version of routine specifically for MPP
CLL                   P.Burton
!LL   4.2    16/08/96  Add TYPFLDPT arguments to FILTER subroutine
!LL                    and make the FILTER_WAVE_NUMBER arrays
!LL                    globally sized                       P.Burton
!LL   4.2    10/01/97  Initialise unprocessed points in THETAL_PROV.
!LL                    D. Robinson.
!LL 4.3      24/04/97 Fixes to 4th order calculations   P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P121
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (35) AND (40)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE TH_ADV
     1              (THETAL,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,
     2              SEC_P_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,
     3              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     4              LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5              SOUTHERN_FILTERED_P_ROW,P_LEVELS,
     6              U_FIELD,P_FIELD,ROW_LENGTH,
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
     7              TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,SEC_U_LATITUDE,
     8              AKH,BKH,QCL,QCF,P_EXNER,OMEGA,
     9              Q_LEVELS,AK,BK,L_SECOND,
     &              extended_address,
     &              LWHITBROM)

      IMPLICIT NONE

! All TYPFLDPT arguments are intent IN
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

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN.
     &, FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &       ! LAST WAVE NUMBER NOT TO BE CHOPPED
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.

C LOGICAL VARIABLE
      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     & ,LWHITBROM
      INTEGER extended_address(P_FIELD)

      REAL
     * THETAL(P_FIELD,P_LEVELS)  !INOUT THETAL FIELD
     *                           ! MASS-WEIGHTED ON OUTPUT.

      REAL
     * U_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                          !   FROM ADJUSTMENT STEP
     *,V_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                          !   * COS(LAT) FROM ADJUSTMENT STEP
     *,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *,PSTAR(P_FIELD)           !IN PSTAR FIELD AT NEW TIME-LEVEL
     *,PSTAR_OLD(P_FIELD)       !IN PSTAR AT PREVIOUS TIME-LEVEL
     *,RS(P_FIELD,P_LEVELS)     !IN RS FIELD
     *,QCL(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCL
     *,QCF(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCF
     *,OMEGA(U_FIELD,P_LEVELS)  !IN. TRUE VERTICAL VELOCITY DP/DT
     *,P_EXNER(P_FIELD,P_LEVELS+1) !IN. PRIMARY ARRAY FOR EXNER FUNCTION

      REAL
     * DELTA_AK(P_LEVELS)      !IN LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)      !IN LAYER THICKNESS
     *,AK(P_LEVELS)            !IN HYBRID CO-ORDINATE AT FULL LEVELS
     *,BK(P_LEVELS)            !IN HYBRID CO-ORDINATE AT FULL LEVELS
     *,AKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,BKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,SEC_P_LATITUDE(P_FIELD) !IN 1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN 1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,LONGITUDE_STEP_INVERSE  !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE   !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP      !IN
     *,NU_BASIC                !IN STANDARD NU TERM FOR MODEL RUN.
     *,TRIGS(ROW_LENGTH)       !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                         ! IN FILTERING.

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 24 ARE REQUIRED.
      REAL
     &    OMEGA_P(P_FIELD,P_LEVELS)    ! HOLDS OMEGA AT P POINTS.

      REAL
     * THETAL_FIRST_INC(P_FIELD,P_LEVELS) ! HOLDS THETAL INCREMENT
     *                           ! RETURNED BY FIRST CALL TO ADV_P_GD
     *,THETAL_SECOND_INC(P_FIELD)! HOLDS THETAL INCREMENT
     *                           ! RETURNED BY SECOND CALL TO ADV_P_GD
     *,THETAL_PROV(P_FIELD,P_LEVELS) ! HOLDS PROVISIONAL VALUE OF
     *                           ! THETAL

      REAL
     * THETAL_INCREMENT(P_FIELD,P_LEVELS) !HOLDS INCREMENT TO THETAL
     *,ZERO(P_FIELD)             !A FIELD OF ZEROES USED WHERE VERTICAL
     *                           !VELOCITY IS ZERO.

      REAL
     * NUX(P_FIELD,P_LEVELS) !COURANT NBR DEPENDENT NU AT P PTS USED
     *                    ! IN EAST-WEST ADVECTION.
     *,NUY(P_FIELD,P_LEVELS) !COURANT NBR DEPENDENT NU AT P PTS USED
     *                    ! IN NORTH-SOUTH ADVECTION.

      REAL NUX_MIN(upd_P_ROWS), ! minimum value of NUX along a row
     &     NUY_MIN(ROW_LENGTH)  ! min of NUY along a column

      REAL
     * BRSP(P_FIELD,P_LEVELS) !MASS TERM AT LEVEL K

! Work space required to allow the use of Fourth Order Advection
! U/V_MEAN_COPY and T_COPY arrays are defined with an extra halo
! this is required for the bigger stencil of the 4th order operator.
      REAL U_MEAN_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                 (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of U_MEAN with extra halo space for 4th order
     &      V_MEAN_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                  (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of V_MEAN with extra halo space for 4th order
     &      T_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &             (tot_P_ROWS+2*extra_NS_Halo),P_LEVELS)
     &  !    Copy of THETAL with extra halo space for 4th order

      INTEGER  extended_P_FIELD,
     &         extended_U_FIELD
!  These are the sizes of the arrays with the extra halos
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS_UPDATE    ! NUMBER OF P POINTS TO BE UPDATED.
     *                     !  = ROWS*ROWLENGTH
     *, U_POINTS_UPDATE    ! NUMBER OF U POINTS TO BE UPDATED.
     *                     !  = (ROWS-1)*ROWLENGTH
     *, P_POINTS_REQUIRED  ! NUMBER OF P POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT P_POINTS_UPDATE
     *, U_POINTS_REQUIRED  ! NUMBER OF U POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT U_POINTS_UPDATE
     *, START_U_REQUIRED   ! FIRST U POINT OF VALUES REQUIRED TO UPDATE
     *                     ! AT P POINTS UPDATE.
     *, END_U_REQUIRED     ! LAST U POINT OF REQUIRED VALUES.

      INTEGER I_start,I_end  ! loop bounds
      INTEGER info  ! return code from comms

C REAL SCALARS
      REAL
     & SCALAR1,SCALAR2,CONST1,LC_LF,TIMESTEP
     &,PK, PK1         ! Pressure at half levels k and k1 (k1=k-1)
     &,P_EXNER_FULL    ! Exner pressure at full model level

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J,K1,IK,K
     *, FILTER_SPACE ! HORIZONTAL DIMENSION OF SPACE NEEDED IN FILTERING
     *               ! ROUTINE.

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL ADV_P_GD,POLAR,UV_TO_P,FILTER
C*---------------------------------------------------------------------
CL    CALL COMDECK TO GET PHYSICAL CONSTANTS USED.

CL    COMDECK C_THADV HOLDS PHYSICAL CONSTANTS REQUIRED BY ROUTINE
CL    TH_ADV.
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
CL    END OF COMDECK C_THADV.

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL

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


CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      LC_LF = LC + LF
      P_POINTS_UPDATE   = upd_P_ROWS*ROW_LENGTH
      U_POINTS_UPDATE   = upd_U_ROWS*ROW_LENGTH
      P_POINTS_REQUIRED = (upd_P_ROWS+2)*ROW_LENGTH
      U_POINTS_REQUIRED = (upd_U_ROWS+2)*ROW_LENGTH
      START_U_REQUIRED  = START_POINT_NO_HALO-ROW_LENGTH
      END_U_REQUIRED    = END_U_POINT_NO_HALO+ROW_LENGTH


C *IF -DEF,NOWHBR replaced by LWHITBROM logical
      IF (LWHITBROM) THEN
CL    CALCULATE BRSP TERM AT LEVEL K

      K=1
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=(3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))
     *                *BKH(K+1)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO
      K=P_LEVELS
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=-(3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))
     *                *BKH(K)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO

      DO K=2,P_LEVELS -1
! Loop over entire field
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          BRSP(I,K)=((3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))*BKH(K+1)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
     *              -((3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))*BKH(K)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
        ENDDO

      ENDDO
      END IF
C *ENDIF

      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        ZERO(I) = 0.
      ENDDO

! In order to use the same call to adv_p_gd for both the second and
! fourth order advection, U/V_MEAN are copied into _COPY arrays.
! In the case of second order advection some of the work space is
! wasted as there is more halo than we need.

! Calculate the size of the extended arrays which contain an
! extra halo:
      extended_U_FIELD=(ROW_LENGTH+2*extra_EW_Halo)*
     &                 (tot_U_ROWS+2*extra_NS_Halo)
      extended_P_FIELD=(ROW_LENGTH+2*extra_EW_Halo)*
     &                 (tot_P_ROWS+2*extra_NS_Halo)

      IF (L_SECOND) THEN

! Copy U/V_MEAN to U/V_MEAN_COPY with the same sized halos
        CALL COPY_FIELD(U_MEAN,U_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)
        CALL COPY_FIELD(V_MEAN,V_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)

      ELSE  ! if its fourth order:

        CALL COPY_FIELD(U_MEAN,U_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(V_MEAN,V_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(THETAL,T_COPY,
     &                  P_FIELD,extended_P_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

       ENDIF ! IF (L_SECOND)

CL LOOP OVER P_LEVELS+1.
CL    ON 1 TO P_LEVELS PROVISIONAL VALUES OF THE FIELD ARE CALCULATED.
CL    ON 2 TO P_LEVELS+1 THE FINAL INCREMENTS ARE CALCULATED AND ADDED
CL    ON. THE REASON FOR THIS LOGIC IS THAT THE PROVISIONAL VALUE AT
CL    LEVEL K+1 IS NEEDED BEFORE THE FINAL INCREMENT AT LEVEL K CAN BE
CL    CALCULATED.

      DO K=1,P_LEVELS+1
CL SET TIMESTEP APPROPRIATE TO LEVEL

        TIMESTEP = ADVECTION_TIMESTEP
CL IF NOT AT P_LEVELS+1 THEN
        IF(K.LE.P_LEVELS) THEN

CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE COURANT NUMBER DEPENDENT NU IF IN
CL                   FORECAST MODE.
CL                   CALCULATE
CL                   PROVISIONAL VALUES OF THETAL AT NEW TIME-LEVEL.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 2.1    SET NU TO NU_BASIC DEPENDENT ON MAX COURANT
CL                   NUMBER.
C ---------------------------------------------------------------------
CL    IF NU_BASIC IS ZERO THEN DO NOT BOTHER TO CALCULATE NU
          IF(.NOT.L_SECOND) THEN
CL    CALCULATE COURANT NUMBER
C NOTE: RS AND TRIG TERMS WILL BE INCLUDED AFTER INTERPOLATION TO P
C       GRID.
CL    CALL UV_TO_P TO MOVE MEAN VELOCITIES ONTO P GRID

          CALL UV_TO_P(U_MEAN(START_U_REQUIRED,K),
     *                 NUX(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

          CALL UV_TO_P(V_MEAN(START_U_REQUIRED,K),
     *                 NUY(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

CL    CALCULATE NU FROM COURANT NUMBER INCLUDING TRIG AND RS TERMS.
          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            NUX(I,K) = NUX(I,K)*LONGITUDE_STEP_INVERSE
            NUY(I,K) = NUY(I,K)*LATITUDE_STEP_INVERSE
            SCALAR1 = TIMESTEP/(RS(I,K)*
     *                RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            SCALAR2 = SEC_P_LATITUDE(I)*SCALAR1
            SCALAR1 = SCALAR1*SCALAR1
            SCALAR2 = SCALAR2*SCALAR2
            NUX(I,K) = (1. - NUX(I,K)*NUX(I,K)*SCALAR2)*NU_BASIC
            NUY(I,K) = (1. - NUY(I,K)*NUY(I,K)*SCALAR1)*NU_BASIC
          ENDDO

! Set NUX equal to minimum value along each row

          DO J=FIRST_ROW,FIRST_ROW+upd_P_ROWS-1
            I_start=(J-1)*ROW_LENGTH+FIRST_ROW_PT ! start and end of rpw
            I_end=(J-1)*ROW_LENGTH+LAST_ROW_PT    ! missing out halos
! Calculate minimum along this row
            SCALAR1=NUX(I_start,K)
            DO I=I_start+1,I_end
              IF (NUX(I,K) .LT. SCALAR1) SCALAR1=NUX(I,K)
            ENDDO
            NUX_MIN(J-FIRST_ROW+1)=SCALAR1
! The indexing of NUX_MIN goes from 1..ROWS
          ENDDO ! J : loop over rows

! So far we have only calculated the minimum along our local
! part of the row. Now we must find the minimum of all the
! local minimums along the row
          CALL GCG_RMIN(upd_P_ROWS,GC_ROW_GROUP,info,NUX_MIN)

! and now set all values of NUX to the minimum along the row
          DO J=FIRST_ROW,FIRST_ROW+upd_P_ROWS-1
            IF (NUX_MIN(J-FIRST_ROW+1) .LT. 0.0)
     &        NUX_MIN(J-FIRST_ROW+1)=0.0

            I_start=(J-1)*ROW_LENGTH+1  ! beginning and
            I_end=J*ROW_LENGTH          ! end of row

            DO I=I_start,I_end
              NUX(I,K)=NUX_MIN(J-FIRST_ROW+1)
            ENDDO

          ENDDO ! J : loop over rows

! Set NUY equal to minimum value along each column

          DO J=FIRST_ROW_PT,LAST_ROW_PT
            I_start=(FIRST_ROW-1)*ROW_LENGTH+J
! I_start points to the beginning of column J

! Calculate the minimum along this column
            I_end=I_start+(upd_P_ROWS-1)*ROW_LENGTH
! I_end points to the end of column J
            SCALAR1=NUY(I_start,K)
            DO I=I_start+ROW_LENGTH,I_end,ROW_LENGTH
              IF (NUY(I,K) .LT. SCALAR1) SCALAR1=NUY(I,K)
            ENDDO
            NUY_MIN(J)=SCALAR1

          ENDDO ! J : loop over columns
! Once again, this is only the minimum along our local part
! of each column. We must now find the miniumum of all the local
! minimums along the column
          CALL GCG_RMIN(ROW_LENGTH-2*EW_Halo,GC_COL_GROUP,info,
     &                  NUY_MIN(EW_Halo+1))

! and now set all values of NUY to the minimum along the column
          DO J=FIRST_ROW_PT,LAST_ROW_PT
            IF (NUY_MIN(J) .LT. 0.0) NUY_MIN(J)=0.0

            I_start=(FIRST_ROW-1)*ROW_LENGTH+J
            I_end=I_start+(upd_P_ROWS-1)*ROW_LENGTH

            DO I=I_start,I_end,ROW_LENGTH
              NUY(I,K)=NUY_MIN(J)
            ENDDO

          ENDDO ! J : loop over columns

        ENDIF  ! IF its fourth order advection


C ---------------------------------------------------------------------
CL    SECTION 2.2    CALL ADV_P_GD TO OBTAIN FIRST INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

CL    CALL ADV_P_GD FOR THETAL.
          K1=K+1

          IF(K.EQ.P_LEVELS) THEN
C PASS ANY THETAL VALUES AS THOSE APPARENTLY AT LEVEL K+1 AS ETADOT
C IS SET TO ZERO BY USING ARRAY ZERO.
            K1 = K-1

          CALL ADV_P_GD(THETAL(1,K1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),
     &                  ETADOT_MEAN(1,K),ZERO,SEC_P_LATITUDE,
     *                  THETAL_FIRST_INC(1,K),NUX(1,K),NUY(1,K),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,SEC_U_LATITUDE,
     *                  BRSP(1,K),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  
          ELSE IF(K.EQ.1)THEN

C PASS ANY THETAL VALUES FOR LEVEL K-1 AS ETADOT AT LEVEL 1
C IS SET TO ZERO BY USING ARRAY ZERO.
          CALL ADV_P_GD(THETAL(1,K1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),ZERO,
     *                  ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,THETAL_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  
          ELSE
          CALL ADV_P_GD(THETAL(1,K-1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),
     &                  ETADOT_MEAN(1,K),ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,THETAL_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  

          END IF


C ---------------------------------------------------------------------
CL    SECTION 2.3    REMOVE MASS-WEIGHTING FROM INCREMENT AND ADD ONTO
CL                   FIELD TO OBTAIN PROVISIONAL VALUE.
C ---------------------------------------------------------------------

          DO I=1,START_POINT_NO_HALO-1
            THETAL_PROV(I,K) = 0.0
          ENDDO

          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            SCALAR1 = RS(I,K)*RS(I,K)
     *                      *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I))
            THETAL_FIRST_INC(I,K)=THETAL_FIRST_INC(I,K)/SCALAR1
            THETAL_PROV(I,K) = THETAL(I,K)- THETAL_FIRST_INC(I,K)
          ENDDO

          DO I=END_P_POINT_NO_HALO+1,P_FIELD
            THETAL_PROV(I,K) = 0.0
          ENDDO

CL   GLOBAL MODEL CALCULATE PROVISIONAL POLAR VALUE.
          IF (at_top_of_LPG) THEN
! North Pole
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              THETAL_PROV(I,K) = THETAL(I,K)
              THETAL_FIRST_INC(I,K) = -THETAL_FIRST_INC(I,K)/
     &                         (RS(I,K)*RS(I,K)*
     &                         (DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            ENDDO
          ENDIF

          IF (at_base_of_LPG) THEN
! South Pole
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              THETAL_PROV(I,K) = THETAL(I,K)
              THETAL_FIRST_INC(I,K) = -THETAL_FIRST_INC(I,K)/
     &                         (RS(I,K)*RS(I,K)*
     &                         (DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            ENDDO
          ENDIF

        END IF
CL END CONDITIONAL ON LEVEL BEING LESS THAN P_LEVELS+1
      ENDDO


      CALL POLAR(THETAL_PROV,THETAL_FIRST_INC,THETAL_FIRST_INC,
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
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,P_LEVELS)

      IF (L_SECOND) THEN

! Do a halo update on the THETAL_PROV array
! that has just been calculated
        CALL SWAPBOUNDS(THETAL_PROV,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
!        CALL SET_SIDES(THETAL_PROV,P_FIELD,ROW_LENGTH,P_LEVELS,
!     &                 fld_type_p)

      ELSE  ! fourth order advection

! Copy THETAL_PROV into T_COPY which has double halos for fourth
! order advection, and do swap to fill these halos
        CALL COPY_FIELD(THETAL_PROV,T_COPY,
     &                  P_FIELD,extended_P_FIELD,
     &                  ROW_LENGTH,tot_P_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

      ENDIF

! Set up OMEGA_P array
! (Was SECTION 3.2):
!    SECTION 3.2    INTERPOLATE OMEGA TO P GRID AND CALCULATE
!                   REMAINING TERM IN ADVECTION EQUATION.
!                   CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD.

            DO K1=1,P_LEVELS
              CALL UV_TO_P(OMEGA(START_U_REQUIRED,K1),
     &                     OMEGA_P(START_POINT_NO_HALO,K1),
     &                     U_POINTS_REQUIRED,
     &                     P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)
              IF (at_top_of_LPG) THEN
                DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
                  OMEGA_P(I,K1)=0.
                ENDDO
              ENDIF

              IF (at_base_of_LPG) THEN
                DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
                  OMEGA_P(I,K1)=0.
                ENDDO
              ENDIF
            ENDDO

              CALL POLAR(OMEGA_P,OMEGA_P,OMEGA_P,
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
     &                   P_FIELD,P_FIELD,P_FIELD,
     &                   START_POINT_NO_HALO,
     &                   END_P_POINT_NO_HALO-ROW_LENGTH+1,
     &                   ROW_LENGTH,P_LEVELS)

CL BEGIN CONDITIONAL ON LEVEL BEING GREATER THAN 1

      DO K=1,P_LEVELS+1
        IF(K.GT.1) THEN
CL---------------------------------------------------------------------
CL    SECTION 3.     ALL WORK IN THIS SECTION PERFORMED AT LEVEL-1.
CL                   CALCULATE SECOND INCREMENT DUE TO ADVECTION.
CL                   CALCULATE TOTAL INCREMENT TO FIELD AND FILTER
CL                   WHERE NECESSARY THEN UPDATE FIELD.
CL                   THE POLAR INCREMENTS ARE THEN CALCULATED AND ADDED
CL                   ON BY CALLING POLAR.
CL---------------------------------------------------------------------

         TIMESTEP = ADVECTION_TIMESTEP

         CONST1 = R/(CP*CP)*TIMESTEP
C ---------------------------------------------------------------------
CL    SECTION 3.1    CALL ADV_P_GD TO OBTAIN SECOND INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

CL    CALL ADV_P_GD FOR THETAL.
          K1=K-1
C K1 HOLDS K-1.
          IF(K.GT.P_LEVELS) THEN
C THE ZERO VERTICAL FLUX AT THE TOP IS ENSURED BY PASSING ETADOT AS
C ZERO.

          CALL ADV_P_GD(THETAL_PROV(1,K-2),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K-2),
     *                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),
     &                  ETADOT_MEAN(1,K-1),ZERO,SEC_P_LATITUDE,
     *                  THETAL_SECOND_INC,NUX(1,K-1),NUY(1,K-1),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)

          ELSE IF(K.EQ.2) THEN
C THE ZERO VERTICAL FLUX AT THE BOTTOM IS ENSURED BY PASSING ETADOT AS
C ZERO.
          CALL ADV_P_GD(THETAL_PROV(1,K),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K),
     *                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),ZERO,
     *                  ETADOT_MEAN(1,K),
     *                  SEC_P_LATITUDE,THETAL_SECOND_INC,
     *                  NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)
          ELSE

          CALL ADV_P_GD(THETAL_PROV(1,K-2),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K),
     *                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),
     &                  ETADOT_MEAN(1,K-1),ETADOT_MEAN(1,K),
     *                  SEC_P_LATITUDE,THETAL_SECOND_INC,
     *                  NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  T_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)
          END IF


C TOTAL MASS-WEIGHTED HORIZONTAL AND VERTICAL INCREMENTS ARE CALCULATED
C SEPARATELY.

          IF(K.LT.Q_LEVELS+2) THEN
            DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO

              PK  = AKH(K)  + BKH(K) *PSTAR(I)
              PK1 = AKH(K1) + BKH(K1)*PSTAR(I)  !  K1 = K-1
              P_EXNER_FULL = P_EXNER_C
     *        (P_EXNER(I,K),P_EXNER(I,K1),PK,PK1,KAPPA)

              THETAL_INCREMENT(I,K1) = .5*(THETAL_SECOND_INC(I) +
     *                       THETAL_FIRST_INC(I,K-1)*RS(I,K1)*RS(I,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I)))
     *                      -(LC*QCL(I,K1)+LC_LF*QCF(I,K1))*CONST1*
     &                       OMEGA_P(I,K1)/((AK(K1)+BK(K1)*PSTAR(I))
     *                       *(P_EXNER_FULL))

            ENDDO
          ELSE
            DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
              THETAL_INCREMENT(I,K1) = .5*(THETAL_SECOND_INC(I) +
     *                       THETAL_FIRST_INC(I,K-1)*RS(I,K1)*RS(I,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I)))
            ENDDO
          END IF

C ---------------------------------------------------------------------
CL    SECTION 3.3    IF GLOBAL MODEL CALCULATE POLAR INCREMENTS.
CL                   IF LIMITED AREA MASS-WEIGHT BOUNDARIES.
C ---------------------------------------------------------------------

CL    GLOBAL MODEL CALCULATE POLAR INCREMENT.
CL    CALCULATE MERIDIONAL FLUX AROUND POLES BY ADDING THE TWO
CL    INCREMENTS AND ALSO MASS-WEIGHTING POLAR FIELDS.
C NEGATIVE SIGN BEFORE FIRST INCS IS DUE TO THEIR SIGN HAVING BEEN
C CHANGED PRIOR TO THE CALCULATION OF THE INTERMEDIATE VALUE.
          IF (at_top_of_LPG) THEN
! Northern boundary/pole
            IF (K.LT.Q_LEVELS+2) THEN
              DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
                SCALAR1 = RS(I,K1)*RS(I,K1)*
     &            (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
                PK  = AKH(K)  + BKH(K) *PSTAR(I)
                PK1 = AKH(K1) + BKH(K1)*PSTAR(I)  !  K1 = K-1
                P_EXNER_FULL = P_EXNER_C(P_EXNER(I,K),
     &                                   P_EXNER(I,K1),PK,PK1,KAPPA)

                THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     &                           - THETAL_FIRST_INC(I,K-1)*SCALAR1)
     &                           +(LC*QCL(I,K1)+LC_LF*QCF(I,K1))*CONST1*
     &                           OMEGA_P(I,K1)/((AK(K1)+BK(K1)*PSTAR(I))
     &                           *P_EXNER_FULL)
                THETAL(I,K1) = THETAL(I,K1)*SCALAR1
              ENDDO
            ELSE  ! (IF K.GE.Q_LEVELS+2)
              DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
                SCALAR1 = RS(I,K1)*RS(I,K1)*
     &                    (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
                THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     &                              - THETAL_FIRST_INC(I,K-1)*SCALAR1)
                THETAL(I,K1) = THETAL(I,K1)*SCALAR1
              ENDDO
            ENDIF ! (K.LT.Q_LEVELS+2)
          ENDIF ! (attop)

          IF (at_base_of_LPG) THEN
! Southern boundary/pole
            IF (K.LT.Q_LEVELS+2) THEN
              DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
                SCALAR2 = RS(I,K1)*RS(I,K1)*
     &            (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
                PK  = AKH(K)  + BKH(K) *PSTAR(I)
                PK1 = AKH(K1) + BKH(K1)*PSTAR(I)  !  K1 = K-1
                P_EXNER_FULL = P_EXNER_C(P_EXNER(I,K),
     &                                   P_EXNER(I,K1),PK,PK1,KAPPA)

                THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     &                           - THETAL_FIRST_INC(I,K-1)*SCALAR2)
     &                           +(LC*QCL(I,K1)+LC_LF*QCF(I,K1))*CONST1*
     &                           OMEGA_P(I,K1)/((AK(K1)+BK(K1)*PSTAR(I))
     &                           *P_EXNER_FULL)
                THETAL(I,K1) = THETAL(I,K1)*SCALAR2
              ENDDO
            ELSE  ! (IF K.GE.Q_LEVELS+2)
              DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
                SCALAR2 = RS(I,K1)*RS(I,K1)*
     &                    (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
                THETAL(I,K1) = THETAL(I,K1)*SCALAR2
                THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     &                              - THETAL_FIRST_INC(I,K-1)*SCALAR2)
              ENDDO
            ENDIF ! (K.LT.Q_LEVELS+2)
          ENDIF ! (atbase)

CL END CONDITIONAL LEVEL GREATER THAN ONE
        END IF

CL END LOOP OVER P_LEVELS+1
      enddo

CL---------------------------------------------------------------------
CL    SECTION 4      IF GLOBAL MODEL THEN FILTER INCREMENTS AND
CL                   UPDATE POLAR VALUES BY CALLING POLAR.
CL                   UPDATE ALL OTHER VALUES.
CL---------------------------------------------------------------------


C ---------------------------------------------------------------------
CL    SECTION 4.1    CALL FILTER TO DO FILTERING.
C ---------------------------------------------------------------------

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE = (ROW_LENGTH+2)*(NORTHERN_FILTERED_P_ROW-1+
     *                tot_P_ROWS-SOUTHERN_FILTERED_P_ROW)
CL    CALL FILTER FOR THETAL INCREMENTS

      CALL FILTER(THETAL_INCREMENT,P_FIELD,P_LEVELS,
     &            FILTER_SPACE,ROW_LENGTH,
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
     &            FILTER_WAVE_NUMBER_P_ROWS,TRIGS,IFAX,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW)

C ---------------------------------------------------------------------
CL    SECTION 4.2    CALL POLAR TO UPDATE POLAR VALUES
C ---------------------------------------------------------------------

      CALL POLAR(THETAL,THETAL_INCREMENT,THETAL_INCREMENT,
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
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,P_LEVELS)

C ---------------------------------------------------------------------
CL    SECTION 4.3    UPDATE ALL OTHER POINTS.
C   OUTPUT IS MASS-WEIGHTED.
C   INCREMENTS ARE ALREADY MASS-WEIGHTED
C ---------------------------------------------------------------------

      DO K=1,P_LEVELS
C UPDATE THETAL.
CFPP$ SELECT(CONCUR)
        DO I= START_POINT_NO_HALO,END_P_POINT_NO_HALO
          THETAL(I,K)=THETAL(I,K)*RS(I,K)*RS(I,K)*
     &        (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))-THETAL_INCREMENT(I,K)
        ENDDO
      ENDDO

CL    END OF ROUTINE TH_ADV

      RETURN
      END
