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
CLL   SUBROUTINE DIF_CTL ------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES AND ADDS DIFFUSION INCREMENTS TO U,V, QT
CLL              AND THETAL USING EQUATIONS (47) AND (48). ONE MORE
CLL              PRESSURE THAN VELOCITY ROW IS UPDATED.
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WRITTEN  BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with DIFCTL1A.
CLL   4.4    25/07/97 Calling sequence changed for UV_DIF, TH_Q_DIF,
CLL                   COEFF_UV, COEF_TH_Q from once per diffusion
CLL                   sweep per level to once per dynamics
CLL                   sweep, in order to improve MPP scalability.
CLL                   A. Dickinson
CLL
CLL
CLL   PROGRAMMING STANDARD:
CLL
CLL   SYSTEM COMPONENTS COVERED: P13
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (47) AND (48)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLLEND-------------------------------------------------------------
C
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE DIF_CTL
     1                  (PSTAR,U,V,THETAL,QT,RS_SQUARED_DELTAP,K1,K2,
     2                   KEXP_K1,KEXP_K2,
     &                   DELTA_AK,DELTA_BK,AK,BK,ADVECTION_TIMESTEP,
     3                   COS_U_LATITUDE,COS_P_LATITUDE,SEC_U_LATITUDE,
     4                   SEC_P_LATITUDE,LONGITUDE_STEP_INVERSE,P_FIELD,
     5                   LATITUDE_STEP_INVERSE,U_FIELD,ROW_LENGTH,
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
     6                   P_LEVELS,Q_LEVELS,
     7                   COS_U_LONGITUDE,SIN_U_LONGITUDE,
     8                   PRESSURE_ALTITUDE,L_TRACER_THETAL_QT)

      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, P_LEVELS           !IN NUMBER OF MODEL LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST MODEL LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     &, KEXP_K1(P_LEVELS)  !IN. EXPONENT OF DIFFUSION SCHEME FOR U,V
     &                     !    AND THETAL FIELDS.
     &, KEXP_K2(Q_LEVELS)  !IN. EXPONENT OF DIFFUSION SCHEME FOR
     &                     !    QT FIELD.

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

      REAL
     * U(U_FIELD,P_LEVELS)       !INOUT U VELOCITY FIELD.
     *,V(U_FIELD,P_LEVELS)       !INOUT V VELOCITY FIELD.
     *,THETAL(P_FIELD,P_LEVELS)  !INOUT THETAL FIELD.
     *,QT(P_FIELD,Q_LEVELS)      !INOUT QT FIELD.

      REAL
     * PSTAR(P_FIELD)            !IN PSTAR FIELD.
     *,RS_SQUARED_DELTAP(P_FIELD,P_LEVELS) !IN RS*RS*DELTA P
     *,COS_U_LATITUDE(U_FIELD)   !IN COS(LAT) AT U POINTS.
     *,COS_P_LATITUDE(P_FIELD)   !IN COS(LAT) AT P POINTS.
     *,SEC_U_LATITUDE(U_FIELD)   !IN 1/COS(LAT) AT U POINTS.
     *,SEC_P_LATITUDE(P_FIELD)   !IN 1/COS(LAT) AT P POINTS.
     *,COS_U_LONGITUDE(ROW_LENGTH) !IN COS(LONGITUDE) AT U POINTS.
     *,SIN_U_LONGITUDE(ROW_LENGTH) !IN SIN(LONGITUDE) AT U POINTS.

      REAL
     * DELTA_AK(P_LEVELS)       !IN LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)       !IN LAYER THICKNESS
     *,AK(P_LEVELS)             !LAYER AK'S
     *,BK(P_LEVELS)             !LAYER BK'S
     *,K1(P_LEVELS)             !IN DIFFUSION COEFF SEE EQ. (45)
     *,K2(P_LEVELS)             !IN DIFFUSION COEFF SEE EQ. (45)
     *,LONGITUDE_STEP_INVERSE   !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE    !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP       !IN
     *, PRESSURE_ALTITUDE       ! ALTITUDE FOR HIGHEST SLOPE TEST

      LOGICAL
     * L_TRACER_THETAL_QT       ! T if tracer advn. for thetal,qt
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 13 ARE REQUIRED
      REAL

     * RECIP_RS_SQUARED_DELTAP(P_FIELD,P_LEVELS)  ! 1./RS*RS*DELTA P
     *,PSTAR_UV(U_FIELD)                 ! HOLDS PRESSURE AT U POINTS.
     *, PRESSURE(P_FIELD,P_LEVELS)    !3-D PRESSURE ARRAY FOR TESTING
     *     ! SLOPE. LEVEL_P=1 IS SURFACE, LEVEL_P=K IS LEVEL K-1
     *     ! FOR UV POINTS PRESSURE RE-CALCULATED TO UV POINTS
     *, DIFFUSION_EW(P_FIELD,P_LEVELS)
!HOLDS EFFECTIVE EAST-WEST DIFFUSION
     *                                   ! COEFFICIENT
     *, DIFFUSION_NS(P_FIELD,P_LEVELS)
!HOLDS EFFECTIVE NORTH-SOUTH DIFFUSION
     *                                   ! COEFFICIENT
     *,COS_FUNCTION_U(U_FIELD)
     *,COS_FUNCTION_P(P_FIELD)

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     &  START_U_UPDATE   ! First U point to be updated
     &, END_U_UPDATE     ! Last U point to be updated

      REAL
     * SCALAR,PRESSURE_TEST
C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,J,K,LEVEL_BASE

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL
     *  TH_Q_DIF,UV_DIF,P_TO_UV,COEFF_TH_Q,COEFF_UV
C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISE LOCAL VARIABLES AND INTERPOLATE PSTAR
CL                   ONTO U-GRID.
CL---------------------------------------------------------------------

C****************************************************************
C     SET PRESSURE_TEST TO PRESSURE_ALTITUDE ABOVE WHICH HEIGHT
C     NO SLOPE TESTING FOR EFFECTIVE DIFFUSION
C***************************************************************
      PRESSURE_TEST=PRESSURE_ALTITUDE

! Diffusion is a bit different from the other dynamics routines.
! START_U_UPDATE and END_U_UPDATE are different for global and LAM
! models - for the global model they include the polar rows,
! but for the LAM they miss the Northern and Southern rows. So
! for this section of code only, we will keep the START_U_UPDATE
! and END_U_UPDATE, rather than using TYPFLDPT equivalents.

! Update U field over entire field, including poles
      START_U_UPDATE=FIRST_FLD_PT
      END_U_UPDATE=LAST_U_FLD_PT
      SCALAR = LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE/
     1           (LONGITUDE_STEP_INVERSE*LONGITUDE_STEP_INVERSE)
      DO I=FIRST_VALID_PT,LAST_U_VALID_PT
        COS_FUNCTION_U(I) = COS_U_LATITUDE(I)*COS_U_LATITUDE(I)*SCALAR
        COS_FUNCTION_P(I) = COS_P_LATITUDE(I)*COS_P_LATITUDE(I)*SCALAR
      END DO

      DO I=LAST_U_VALID_PT+1,LAST_P_VALID_PT
        COS_FUNCTION_P(I) = COS_P_LATITUDE(I)*COS_P_LATITUDE(I)*SCALAR
      END DO

CL    CALL P_TO_UV
C STORE PSTAR ON U GRID IN PSTAR_UV.

      CALL P_TO_UV(PSTAR,PSTAR_UV,P_FIELD,U_FIELD,ROW_LENGTH,tot_P_ROWS)
! Get correct values in halos
      CALL SWAPBOUNDS(PSTAR_UV,ROW_LENGTH,tot_U_ROWS,EW_Halo,NS_Halo,1)

CL    MAKE 3-D PRESSURE ARRAY AT P POINTS
CL    LEVEL_P=1 SURFACE, LEVEL_P=K IS LEVEL K-1
CL    ONLY NEED P_LEVELS AS SURFACES SHOULD BE PRESSURE SURFACES
CL    NEAR TOP OF MODEL SO TESTING UNNECESSARY
C****************************************************************
C    IF USING TRACER ADVECTION OF THETAL AND QT THEN DIFFUSION IS
C    CALLED FOR TOP LEVEL THETAL AND FOR ALL U'S AND V'S ONLY
C     NOTE STEEP SLOPE TEST SHOULD BE DISABLED BY APPROPRIATE
C     SETTING OF PRESSURE ALTITUDE WHICH MEANS THAT PRESSURES
C    DO NOT NEED CALCULATING FOR PRESSURE ARRAY
C***************************************************************

      IF(.NOT.L_TRACER_THETAL_QT)THEN

CL    FIRST LEVEL
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        RECIP_RS_SQUARED_DELTAP(I,1)=1./RS_SQUARED_DELTAP(I,1)
        PRESSURE(I,1)=PSTAR(I)
       END DO
CL OTHER LEVELS
      DO K=2,P_LEVELS
       DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        RECIP_RS_SQUARED_DELTAP(I,K)=1./RS_SQUARED_DELTAP(I,K)
        PRESSURE(I,K)=AK(K-1)+BK(K-1)*PSTAR(I)
       END DO
      END DO

C  POINTER FOR DIFFUSION LEVEL START
       LEVEL_BASE=1

      ELSE
CLL
C    IF USING TRACER ADVECTION OF THETAL AND QT THEN DIFFUSION IS
C    CALLED FOR TOP LEVEL THETAL AND FOR ALL U'S AND V'S ONLY
C     LEVEL_BASE IS THEN SET TO P_LEVELS OTHERWISE SET TO 1
C    IF NECESSARY THE TEST COULD BE MADE ON THE VALUE OF THE
C    DIFFUSION COEFFICIENT K1 FOR EACH LEVEL
C     NOTE STEEP SLOPE TEST SHOULD BE DISABLED BY APPROPRIATE
C     SETTING OF PRESSURE ALTITUDE WHICH MEANS THAT PRESSURES
C    DO NOT NEED CALCULATING FOR PRESSURE ARRAY
CLL
       LEVEL_BASE=P_LEVELS

      END IF


CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE DIFFUSION OF THETAL.
CL                   ADD ON INCREMENT TO ALL POINTS EXCEPT POLES
CL                   WHICH WOULD HAVE BEEN DONE INSIDE TH_Q_DIF.
CL---------------------------------------------------------------------


C CALL COEFF_TH_Q FOR EFFECTIVE DIFFUSION COEFFICIENT FOR THETAL
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE IN DIFFUSION_EW AND DIFFUSION_NS

       CALL COEFF_TH_Q
     1                  (DIFFUSION_EW,DIFFUSION_NS,
     2                   PRESSURE,LEVEL_BASE,PRESSURE_TEST,AK,BK,
     3                   COS_U_LATITUDE,ROW_LENGTH,
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
     5                   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                   P_FIELD,U_FIELD,P_LEVELS,
     7                   K1,DELTA_AK,DELTA_BK,PSTAR_UV,COS_FUNCTION_U)


CL---------------------------------------------------------------------
CL      CALL TH_Q_DIF
CL
CL---------------------------------------------------------------------
CL    NEW VERSION INCLUDES PRESSURE TEST ON SLOPES

          CALL TH_Q_DIF(THETAL,RECIP_RS_SQUARED_DELTAP,
     &                  SEC_P_LATITUDE,ROW_LENGTH,
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
     &                  LEVEL_BASE,P_LEVELS,
     &                  KEXP_K1,ADVECTION_TIMESTEP,
     &                  P_FIELD,U_FIELD,
     &                  DIFFUSION_EW,DIFFUSION_NS)


CL
CL---------------------------------------------------------------------
CL    SECTION 4.     CALCULATE DIFFUSION OF QT AND
CL                   ADD ON INCREMENT TO ALL POINTS EXCEPT POLES
CL                   WHICH WOULD HAVE BEEN DONE INSIDE TH_Q_DIF.
CL---------------------------------------------------------------------

      IF(.NOT.L_TRACER_THETAL_QT)THEN

C CALL COEFF_TH_Q FOR EFFECTIVE DIFFUSION COEFFICIENT FOR QT
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE IN DIFFUSION_EW AND DIFFUSION_NS
       CALL COEFF_TH_Q
     1                  (DIFFUSION_EW,DIFFUSION_NS,
     2                   PRESSURE,1,PRESSURE_TEST,AK,BK,
     3                   COS_U_LATITUDE,ROW_LENGTH,
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
     5                   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                   P_FIELD,U_FIELD,Q_LEVELS,
     7                   K2,DELTA_AK,DELTA_BK,PSTAR_UV,COS_FUNCTION_U)

CL---------------------------------------------------------------------
CL      CALL TH_Q_DIF AT A MOIST LEVEL.
CL
CL---------------------------------------------------------------------
CL    NEW VERSION INCLUDES PRESSURE TEST ON SLOPES

            CALL TH_Q_DIF(QT,RECIP_RS_SQUARED_DELTAP,
     &                    SEC_P_LATITUDE,ROW_LENGTH,
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
     &                    1,Q_LEVELS,
     &                    KEXP_K2,ADVECTION_TIMESTEP,
     &                    P_FIELD,U_FIELD,
     &                    DIFFUSION_EW,DIFFUSION_NS)
C
CL END IF TEST FOR NO DIFFUSION WITH TRACER ADVECTION
      END IF

CL   MAKE 3-D PRESSURE ARRAY AT UV POINTS
CL   LEVEL_P=1 SURFACE, LEVEL_P=K IS LEVEL K-1
CL    ONLY NEED P_LEVELS AS SURFACES SHOULD BE PRESSURE SURFACES
CL    NEAR TOP OF MODEL SO TESTING UNNECESSARY
CL   FIRST LEVEL
      DO I=FIRST_VALID_PT,LAST_U_VALID_PT
        PRESSURE(I,1)=PSTAR_UV(I)
       END DO
CL OTHER LEVELS
      DO K=2,P_LEVELS
       DO I=FIRST_VALID_PT,LAST_U_VALID_PT
        PRESSURE(I,K)=AK(K-1)+BK(K-1)*PSTAR_UV(I)
       END DO
      END DO


CL
CL---------------------------------------------------------------------
CL    SECTION 5.     SET DIFFUSION_COEFFICIENTS ON P GRID.
CL                   THEN CALCULATE DIFFUSION OF U AND V.
CL---------------------------------------------------------------------

C CALL COEFF_UV FOR EFFECTIVE DIFFUSION COEFFICIENT FOR U AND V
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE RETURNED IN DIFFUSION_EW AND DIFFUSION_NS
      CALL COEFF_UV
     1                 (DIFFUSION_EW,DIFFUSION_NS,
     2                 PRESSURE,PRESSURE_TEST,AK,BK,
     3                 COS_P_LATITUDE,START_U_UPDATE,END_U_UPDATE,
     &                 ROW_LENGTH,
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
     4                 LATITUDE_STEP_INVERSE,
     5                 LONGITUDE_STEP_INVERSE,P_FIELD,U_FIELD,P_LEVELS,
     6                 K1,DELTA_AK,DELTA_BK,PSTAR,COS_FUNCTION_P)


CL    CALL UV_DIF FOR U &V

        CALL UV_DIF(U,V,RS_SQUARED_DELTAP,
     *              SEC_U_LATITUDE,START_U_UPDATE,END_U_UPDATE,
     &              ROW_LENGTH,
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
     &              P_LEVELS,KEXP_K1,ADVECTION_TIMESTEP,
     *              P_FIELD,U_FIELD,
     *              DIFFUSION_EW,DIFFUSION_NS)

CL    END OF ROUTINE DIF_CTL

      RETURN
      END
