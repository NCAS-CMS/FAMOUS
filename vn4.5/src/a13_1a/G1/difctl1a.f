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
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   3.4    07/08/94 Directives inserted to improve parallel
CLL                   efficiency on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
!     3.5    28/03/95 MPP code: Change updateable area,
!                     add halo updates.   P.Burton
CLL

CLL  4.0  02/02/95: SET EFFECTIVE DIFFUSION TO ZERO WHENEVER
CLL             SLOPE IS CONSIDERED TOO STEEP. THIS REDUCES EXCESSIVE
CLL             PRECIPITATION OVER STEEP OROGRAPHY PROVIDED
CLL             PRESSURE_TEST IS SET TO AN APPROPRIATE ALTITUDE
CLL             E.G. 20000Pa (200hPa)
CLL             Author:  T. DAVIES FR.   Reviewer: M. MAWSON
!     4.1    07/05/96 Added MPP code and TYPFLDPT arguments  P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
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
     * DIFFUSION_COEFFICIENT(P_FIELD)    !HOLDS EAST-WEST DIFFUSION
     *                                   ! COEFFICIENT
     *,DIFFUSION_COEFFICIENT2(P_FIELD)   !HOLDS NORTH-SOUTH DIFFUSION
     *                                   ! COEFFICIENT
     *,QT_INC(P_FIELD)                   ! HOLDS QT INCREMENT
     *,THETAL_INC(P_FIELD)               ! HOLDS THETAL INCREMENT
     *,RS_SQUARED_DELTAP_U_GRID(U_FIELD) ! RS*RS*DELTA P AT U POINTS.
     *,RECIP_RS_SQUARED_DELTAP(P_FIELD)  ! 1./RS*RS*DELTA P
     *,PSTAR_UV(U_FIELD)                 ! HOLDS PRESSURE AT U POINTS.
     *,FIELD1(P_FIELD)                   ! GENERAL WORK-SPACE
     *,FIELD2(P_FIELD)                   ! GENERAL WORK-SPACE
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

! Update U field, missing top and bottom rows
      START_U_UPDATE=START_POINT_NO_HALO
      END_U_UPDATE=END_U_POINT_NO_HALO
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
        PRESSURE(I,1)=PSTAR(I)
       END DO
CL OTHER LEVELS
      DO K=2,P_LEVELS
       DO I=FIRST_VALID_PT,LAST_P_VALID_PT
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

cmic$ do all shared (advection_timestep, cos_function_u)
cmic$*       shared (cos_u_latitude, cos_u_longitude)
cmic$*       shared (delta_ak, delta_bk, end_u_update, k1)
cmic$*       shared (ak, bk, kexp_k1, latitude_step_inverse)
cmic$*       shared (longitude_step_inverse, p_field, u_field, p_levels)
cmic$*       shared (pstar_uv, thetal, row_length)
cmic$*       shared (rs_squared_deltap, sec_p_latitude)
cmic$*       shared (start_u_update)
CMIC@a   SHARED(FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
CMIC@b   P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
CMIC@c   FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
CMIC@d   FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
CMIC@e   VALID_P_ROWS, VALID_U_ROWS,
CMIC@f   START_POINT_NO_HALO, START_POINT_INC_HALO,
CMIC@g   END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
CMIC@h   END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
CMIC@i   FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
CMIC@j   GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD)
cmic$*       shared ( pressure, pressure_test, level_base)
cmic$*       private (diffusion_coefficient, diffusion_coefficient2)
cmic$*       private ( field1, i, j, k, scalar, thetal_inc)
cmic$*       private ( diffusion_ew, diffusion_ns)
cmic$*       private (recip_rs_squared_deltap)

      DO K=LEVEL_BASE,P_LEVELS

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE DIFFUSION OF THETAL.
CL                   ADD ON INCREMENT TO ALL POINTS EXCEPT POLES
CL                   WHICH WOULD HAVE BEEN DONE INSIDE TH_Q_DIF.
CL---------------------------------------------------------------------

C SET DIFFUSION COEFFICIENT AND COPY THETAL INTO FIELD1.
        DO  I=FIRST_VALID_PT,LAST_U_VALID_PT
          DIFFUSION_COEFFICIENT2(I) = K1(K)*
     1                           (DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I))
          DIFFUSION_COEFFICIENT(I) = COS_FUNCTION_U(I)*
     2                               DIFFUSION_COEFFICIENT2(I)
        END DO

C CALL COEFF_TH_Q FOR EFFECTIVE DIFFUSION COEFFICIENT FOR THETAL
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE IN DIFFUSION_EW AND DIFFUSION_NS
       CALL COEFF_TH_Q
     1                  (DIFFUSION_EW(1,K),DIFFUSION_NS(1,K),          
     2                   PRESSURE,K,PRESSURE_TEST,AK,BK,
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
     7                   DIFFUSION_COEFFICIENT,DIFFUSION_COEFFICIENT2)
      ENDDO

      CALL SWAPBOUNDS(DIFFUSION_EW,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)                    
      CALL SWAPBOUNDS(DIFFUSION_NS,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)                     


      DO K=LEVEL_BASE,P_LEVELS                                         
        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT                        
          FIELD1(I) = THETAL(I,K)                                 
          RECIP_RS_SQUARED_DELTAP(I) = 1./RS_SQUARED_DELTAP(I,K)       
        END DO                                                         
C LOOP THROUGH CODE KEXP_K1 TIMES. THE ORDER OF THE DIFFUSION SCHEME IS
C DEL TO THE POWER 2*KEXP_K1.
        DO  J=1,KEXP_K1(K)

CL   ZERO INCREMENTS FOR FIRST AND LAST ROW
CL   OVERWRITTEN BY POLAR IN GLOBAL MODELS
          IF (at_top_of_LPG) THEN
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              THETAL_INC(I)=0.0
            ENDDO
          ENDIF
          IF (at_base_of_LPG) THEN
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              THETAL_INC(I)=0.0
            ENDDO
          ENDIF

CL      CALL TH_Q_DIF
CL
CL---------------------------------------------------------------------
CL    NEW VERSION INCLUDES PRESSURE TEST ON SLOPES
          CALL TH_Q_DIF(FIELD1,THETAL_INC,
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
     &                  P_FIELD,U_FIELD,
     &                  DIFFUSION_EW(1,K),DIFFUSION_NS(1,K))       

C DE-MASS-WEIGHT INCREMENT AND COPY INTO FIELD1 SO THAT IT CAN BE FED
C BACK INTO TH_Q_DIF.
         DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            FIELD1(I) = THETAL_INC(I)*RECIP_RS_SQUARED_DELTAP(I)
         END DO

         if(j.ne.KEXP_K1(K))then
         CALL SWAPBOUNDS(FIELD1,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,1)
         endif

C  END OF DIFFUSION SWEEPS
      END DO

CL ADD FINAL INCREMENT ONTO THETAL FIELD.
        SCALAR = (-1)**KEXP_K1(K)
        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT
          THETAL(I,K) = THETAL(I,K) - FIELD1(I) * ADVECTION_TIMESTEP
     &                   *SCALAR
        END DO

CL END LOOP OVER P_LEVELS FOR THETAL
      END DO
         CALL SWAPBOUNDS
     1  (THETAL,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)

CL
CL---------------------------------------------------------------------
CL    SECTION 4.     CALCULATE DIFFUSION OF QT AND
CL                   ADD ON INCREMENT TO ALL POINTS EXCEPT POLES
CL                   WHICH WOULD HAVE BEEN DONE INSIDE TH_Q_DIF.
CL---------------------------------------------------------------------

      IF(.NOT.L_TRACER_THETAL_QT)THEN

cmic$ do all shared (advection_timestep, cos_function_u)
cmic$*       shared (cos_u_latitude, cos_u_longitude)
cmic$*       shared (delta_ak, delta_bk, end_u_update)
cmic$*       shared (ak, bk, k2, kexp_k2, latitude_step_inverse)
cmic$*       shared (longitude_step_inverse, p_field, u_field, p_levels)
cmic$*       shared (pstar_uv, q_levels, qt, row_length)
cmic$*       shared (rs_squared_deltap, sec_p_latitude)
cmic$*       shared (start_u_update)
CMIC@a   SHARED(FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
CMIC@b   P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
CMIC@c   FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
CMIC@d   FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
CMIC@e   VALID_P_ROWS, VALID_U_ROWS,
CMIC@f   START_POINT_NO_HALO, START_POINT_INC_HALO,
CMIC@g   END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
CMIC@h   END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
CMIC@i   FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
CMIC@j   GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD)
cmic$*       shared (pressure, pressure_test)
cmic$*       private ( field1, i, j, k, scalar, qt_inc)
cmic$*       private ( diffusion_ew, diffusion_ns)
cmic$*       private (diffusion_coefficient, diffusion_coefficient2)
cmic$*       private (recip_rs_squared_deltap)

      DO K=1,Q_LEVELS

C SET DIFFUSION COEFFICIENT AND COPY QT INTO FIELD1.
          DO  I=FIRST_VALID_PT,LAST_U_VALID_PT
            DIFFUSION_COEFFICIENT2(I) = K2(K)*
     1                            (DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I))
            DIFFUSION_COEFFICIENT(I) = COS_FUNCTION_U(I)*
     2                                 DIFFUSION_COEFFICIENT2(I)
          END DO


C CALL COEFF_TH_Q FOR EFFECTIVE DIFFUSION COEFFICIENT FOR QT
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE IN DIFFUSION_EW AND DIFFUSION_NS
       CALL COEFF_TH_Q
     1                  (DIFFUSION_EW(1,K),DIFFUSION_NS(1,K),         
     2                   PRESSURE,K,PRESSURE_TEST,AK,BK,
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
     7                   DIFFUSION_COEFFICIENT,DIFFUSION_COEFFICIENT2)
      ENDDO

      CALL SWAPBOUNDS(DIFFUSION_EW,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,Q_LEVELS)                     
      CALL SWAPBOUNDS(DIFFUSION_NS,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,Q_LEVELS)                     

C LOOP THROUGH CODE KEXP_K2 TIMES. THE ORDER OF THE DIFFUSION SCHEME IS
C DEL TO THE POWER 2*KEXP_K2.
      DO K=1,Q_LEVELS                                                  
          DO  I=FIRST_VALID_PT,LAST_P_VALID_PT                     
          FIELD1(I) = QT(I,K)                                         
          RECIP_RS_SQUARED_DELTAP(I) = 1./RS_SQUARED_DELTAP(I,K)       
          END DO                                                       
          DO J=1,KEXP_K2(K)

CL   ZERO INCREMENTS FOR FIRST AND LAST ROW
CL   OVERWRITTEN BY POLAR IN GLOBAL MODELS
          IF (at_top_of_LPG) THEN
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              QT_INC(I)=0.0
            ENDDO
          ENDIF
          IF (at_base_of_LPG) THEN
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              QT_INC(I)=0.0
            ENDDO
          ENDIF

CL      CALL TH_Q_DIF AT A MOIST LEVEL.

CL
CL---------------------------------------------------------------------
CL    NEW VERSION INCLUDES PRESSURE TEST ON SLOPES
            CALL TH_Q_DIF(FIELD1,QT_INC,
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
     &                    P_FIELD,U_FIELD,
     &                    DIFFUSION_EW(1,K),DIFFUSION_NS(1,K))         
C
CL---------------------------------------------------------------------
C DE-MASS-WEIGHT INCREMENT AND COPY INTO FIELD1 SO THAT IT CAN BE FED
C BACK INTO TH_Q_DIF.
            DO I=FIRST_FLD_PT,LAST_P_FLD_PT
              FIELD1(I) = QT_INC(I)*RECIP_RS_SQUARED_DELTAP(I)
            END DO
         if(J.ne.KEXP_K2(K))then
         CALL SWAPBOUNDS(FIELD1,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,1)
         endif

C  END OF DIFFUSION SWEEPS
         END DO

CL ADD FINAL INCREMENT ONTO QT FIELD.
          SCALAR = (-1)**KEXP_K2(K)
          DO I=FIRST_VALID_PT,LAST_P_VALID_PT
            QT(I,K) = QT(I,K) - FIELD1(I) * ADVECTION_TIMESTEP
     &                     *SCALAR
          END DO

CL END LOOP OVER P_LEVELS FOR  QT
       END DO
         CALL SWAPBOUNDS
     1  (QT,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,Q_LEVELS)
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

CL LOOP OVER P_LEVELS FOR U AND V
cmic$ do all shared (advection_timestep, cos_function_p)
cmic$*       shared (cos_p_latitude, cos_u_longitude, sin_u_longitude)
cmic$*       shared (delta_ak, delta_bk, end_u_update)
cmic$*       shared (ak, bk, k1, kexp_k1, latitude_step_inverse)
cmic$*       shared (longitude_step_inverse, p_field, u_field, p_levels)
cmic$*       shared (pstar, u, v, row_length)
cmic$*       shared (rs_squared_deltap, sec_u_latitude)
cmic$*       shared (start_u_update)
CMIC@a   SHARED(FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
CMIC@b   P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
CMIC@c   FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
CMIC@d   FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
CMIC@e   VALID_P_ROWS, VALID_U_ROWS,
CMIC@f   START_POINT_NO_HALO, START_POINT_INC_HALO,
CMIC@g   END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
CMIC@h   END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
CMIC@i   FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
CMIC@j   GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD)
cmic$*       shared (pressure, pressure_test)
cmic$*       private ( field1, field2, i, j, k, scalar)
cmic$*       private ( diffusion_ew, diffusion_ns)
cmic$*       private (diffusion_coefficient, diffusion_coefficient2)
cmic$*       private (rs_squared_deltap_u_grid)

      DO K=1,P_LEVELS
CL
CL---------------------------------------------------------------------
CL    SECTION 5.     SET DIFFUSION_COEFFICIENTS ON P GRID.
CL                   THEN CALCULATE DIFFUSION OF U AND V.
CL---------------------------------------------------------------------

C SET DIFFUSION COEFFICIENT
        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT
          DIFFUSION_COEFFICIENT2(I) = K1(K)*
     1                           (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
          DIFFUSION_COEFFICIENT(I) = COS_FUNCTION_P(I)*
     2                         DIFFUSION_COEFFICIENT2(I)
        END DO
C CALL COEFF_UV FOR EFFECTIVE DIFFUSION COEFFICIENT FOR U AND V
C AVERAGING IS DONE AS REQUIRED IN EQUATION(48).
C COEFFICIENTS ARE SET TO ZERO FOR STEEP SLOPES
C VALUES ARE RETURNED IN DIFFUSION_EW AND DIFFUSION_NS
      CALL COEFF_UV
     1                 (DIFFUSION_EW(1,K),DIFFUSION_NS(1,K),          
     2                 PRESSURE,K,PRESSURE_TEST,AK,BK,
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
     6                 DIFFUSION_COEFFICIENT,DIFFUSION_COEFFICIENT2)

      ENDDO

      CALL SWAPBOUNDS(DIFFUSION_EW,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)        
      CALL SWAPBOUNDS(DIFFUSION_NS,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)                    

      DO K=1,P_LEVELS                                                  
CL
CL---------------------------------------------------------------------
CL    SECTION 4.     INTERPOLATE RS_SQUARED_DELTAP TO U GRID.
CL---------------------------------------------------------------------

C INTERPOLATE RS_SQUARED_DELTAP TO U GRID.

        CALL P_TO_UV(RS_SQUARED_DELTAP(1,K),RS_SQUARED_DELTAP_U_GRID,
     *                P_FIELD,U_FIELD,ROW_LENGTH,tot_P_ROWS)
        DO I=FIRST_VALID_PT,LAST_U_VALID_PT                            
          FIELD1(I) = U(I,K)                                          
          FIELD2(I) = V(I,K)                                          
        END DO                                                         

                                                                      
C LOOP THROUGH CODE KEXP_K1 TIMES. THE ORDER OF THE DIFFUSION SCHEME IS
C DEL TO THE POWER 2*KEXP_K1.

        DO J=1,KEXP_K1(K)
CL    CALL UV_DIF FOR U &V

        CALL UV_DIF(FIELD1,FIELD2,RS_SQUARED_DELTAP_U_GRID,
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
     *              P_FIELD,U_FIELD,
     *              DIFFUSION_EW(1,K),DIFFUSION_NS(1,K))              
      if(j.ne.KEXP_K1(K))then
      CALL SWAPBOUNDS(FIELD1,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(FIELD2,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,1)
      endif

C    FIELD1 AND FIELD2 NOW CONTAIN DIFFUSED QUANTITIES WHICH CAN
C     BE USED IN FURTHER DIFFUSION SWEEPS

CL   END OF DIFFUSION SWEEPS
        END DO
CL ADD FINAL INCREMENT ONTO WIND FIELDS.
        SCALAR = (-1)**KEXP_K1(K)
! Loop over field, missing top and bottom rows and halos
        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          U(I,K) = U(I,K) - FIELD1(I) * ADVECTION_TIMESTEP
     &                     *SCALAR
          V(I,K) = V(I,K) - FIELD2(I) * ADVECTION_TIMESTEP
     &                     *SCALAR
        END DO
CL END LOOP OVER P_LEVELS

       END DO

      CALL SWAPBOUNDS
     1  (U,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)
      CALL SWAPBOUNDS
     1  (V,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)

CL    END OF ROUTINE DIF_CTL

      RETURN
      END
