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
CLL   SUBROUTINE P_TH_ADJ -------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES ADDS SURFACE PRESSURE INCREMENTS USING
CLL             EQUATION (27). CALCULATES AND ADDS POTENTIAL TEMPERATURE
CLL             INCREMENTS USING EQUATION (28).
CLL   NOT SUITABLE FOR I.B.M USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL M.MAWSON    <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL
CLL   3.4    07/08/94 Directives inserted to improve parallel
CLL                   efficiency on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
!     4.1    02/04/96 Added TYPFLDPT arguments to dynamics routines
!                     which allows many of the differences between
!                     MPP and "normal" code to be at top level
!                     P.Burton
CLL
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL
CLL   SYSTEM COMPONENTS COVERED: P113
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (27) AND (28)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 10 DATED 10/09/90.
CLLEND-------------------------------------------------------------

C
C*L   ARGUMENTS:---------------------------------------------------

      SUBROUTINE P_TH_ADJ
     1                   (PSTAR,PSTAR_OLD,THETA,THETA_REF,
     2                    ETADOT,RS,DELTA_AK,DELTA_BK,
     3                    P_FIELD,P_LEVELS,
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
     4                    CALL_NUMBER,ADJUSTMENT_TIMESTEP,
     5                    ERROR_CODE,ERROR_MESSAGE,
     *                    RECIP_RS_SQUARED_SURFACE,L_NEG_PSTAR)

      IMPLICIT NONE
      LOGICAL
     *  L_NEG_PSTAR    !IN SWITCH, IF TRUE THEN NEGATIVE PSTAR VALUES
     *                 ! WILL BE DETECTED AND OUTPUT.


      INTEGER
     *  P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS OF DATA
     *, P_FIELD            !IN    NUMBER OF POINTS IN PRESSURE FIELD.
     *, CALL_NUMBER        !IN    ADJUSTMENT STEP NUMBER ON WHICH CALL
     *                     !      TO ROUTINE IS BEING MADE.
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
     *  ERROR_CODE         !INOUT. 0 ON ENTRY. 1 ON EXIT IF NEGATIVE
     *                     ! PRESSURE DETECTED.

      CHARACTER*(80)
     *  ERROR_MESSAGE      !OUT. HOLDS ERROR MESSAGE IF ERROR_CODE
     *                     ! NON-ZERO.

      REAL
     * ETADOT(P_FIELD,P_LEVELS)  !IN. HOLDS MASS-WEIGHTED VERTICAL
     *                           ! VELOCITY. AT LEVEL ONE HOLDS SUM
     *                           ! OF DIVERGENCES IN THE COLUMN.
     *,RS(P_FIELD,P_LEVELS)      !IN. RADIUS OF EARTH AT EACH LEVEL.
     *,DELTA_AK(P_LEVELS)        !IN. DIFFERENCE BETWEEN AK'S AT HALF
     *                           !    LEVELS.
     *,DELTA_BK(P_LEVELS)        !IN. DIFFERENCE BETWEEN BK'S AT HALF
     *                           !    LEVELS.
     *,THETA_REF(P_LEVELS)       !IN. REFERENCE THETA PROFILE.
     *,ADJUSTMENT_TIMESTEP       !IN.
     *,RECIP_RS_SQUARED_SURFACE(P_FIELD) !IN. 1/(RS*RS) AT MODEL
     *                                   ! SURFACE.

      REAL
     * PSTAR(P_FIELD)            !INOUT. PRIMARY ARRAY FOR PSTAR
     *,THETA(P_FIELD,P_LEVELS)   !INOUT. PRIMARY ARRAY FOR THETA.

      REAL
     * PSTAR_OLD(P_FIELD)        !OUT. PSTAR AT OLD TIME-LEVEL.

C*---------------------------------------------------------------------

C*L   NO LOCAL ARRAYS NEEDED ------------------------------------------
C*---------------------------------------------------------------------

C DEFINE COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,K
     &, info ! return code from GCOM

C*L   NO EXTERNAL SUBROUTINE CALLS:---------------------------------
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1. IF CALL NUMBER one STORE VALUE OF PSTAR AT OLD
cl               TIME-LEVEL.
CL---------------------------------------------------------------------

      IF (CALL_NUMBER.EQ.1) THEN
C STORE PSTAR AT OLD TIME-LEVEL.
! loop over all points, including valid halos
        DO 100 I=FIRST_VALID_PT,LAST_P_VALID_PT
          PSTAR_OLD(I) = PSTAR(I)
 100    CONTINUE
      END IF


CL
CL---------------------------------------------------------------------
CL    SECTION 2. ADJUST THETA USING EQUATION (28).
CL               THETA ADJUSTMENT IS DONE BEFORE PSTAR AS THETA
CL               ADJUSTMENT REQUIRES PSTAR AT LAST TIME-LEVEL.
CL---------------------------------------------------------------------

C LOOP OVER ALL LEVELS
cmic$ do all shared (adjustment_timestep, delta_ak, delta_bk)
cmic$*       shared (etadot, p_field, p_levels, pstar_old, rs)
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
cmic$*       shared (theta, theta_ref)
cmic$*       private (i, k)
      DO 200 K=1,P_LEVELS

C FOR LIMITED AREA MODEL ADJUST ALL VALUES NOT ON POLEWARDS BOUNDARIES.
C AS ETADOT AT LEVEL 1 AND LEVEL P_LEVELS+1 ARE ZERO AND ALSO ARE
C NOT STORED SLIGHTLY DIFFERENT CODE IS REQUIRED.

        IF(K.EQ.1) THEN

C ADJUST ALL THETA VALUES.

CFPP$ SELECT(CONCUR)
! loop over all points, missing poleward bounds but including halos
          DO 210 I=START_POINT_INC_HALO,END_P_POINT_INC_HALO
            THETA(I,K) = THETA(I,K) - ADJUSTMENT_TIMESTEP * .5 *
     *                            (ETADOT(I,K+1)*(THETA_REF(K+1)-
     *                             THETA_REF(K)))/
     *                            (RS(I,K)*RS(I,K)*(DELTA_AK(K)+
     *                             DELTA_BK(K)*PSTAR_OLD(I)))
 210      CONTINUE

        ELSE IF (K.EQ.P_LEVELS) THEN

C ADJUST ALL THETA VALUES.

CFPP$ SELECT(CONCUR)
! loop over all points, missing poleward bounds but including halos
          DO 220 I=START_POINT_INC_HALO,END_P_POINT_INC_HALO
            THETA(I,K) = THETA(I,K) - ADJUSTMENT_TIMESTEP * .5 *
     *                            (ETADOT(I,K)*
     *                            (THETA_REF(K)-THETA_REF(K-1)))/
     *                            (RS(I,K)*RS(I,K)*(DELTA_AK(K)+
     *                             DELTA_BK(K)*PSTAR_OLD(I)))
 220      CONTINUE

        ELSE

C ADJUST ALL THETA VALUES.

CFPP$ SELECT(CONCUR)
! loop over all points, missing poleward bounds but including halos
          DO 230 I=START_POINT_INC_HALO,END_P_POINT_INC_HALO
            THETA(I,K) = THETA(I,K) - ADJUSTMENT_TIMESTEP * .5 *
     *                            (ETADOT(I,K+1)*(THETA_REF(K+1)-
     *                             THETA_REF(K)) + ETADOT(I,K)*
     *                            (THETA_REF(K)-THETA_REF(K-1)))/
     *                            (RS(I,K)*RS(I,K)*(DELTA_AK(K)+
     *                             DELTA_BK(K)*PSTAR_OLD(I)))
 230      CONTINUE

        END IF


C END LOOP OVER LEVELS
 200  CONTINUE

CL
CL---------------------------------------------------------------------
CL    SECTION 3. ADJUST PSTAR USING EQUATION (27).
CL---------------------------------------------------------------------


C IF LIMITED AREA MODEL ADJUST ALL PRESSURE VALUES NOT ON POLEWARDS
C BOUNDARIES.

! loop over all points, missing poleward bounds but including halos
          DO 300 I=START_POINT_INC_HALO,END_P_POINT_INC_HALO
        PSTAR(I) = PSTAR(I) + ADJUSTMENT_TIMESTEP * ETADOT(I,1)
     *                        *RECIP_RS_SQUARED_SURFACE(I)
 300  CONTINUE


      IF(L_NEG_PSTAR) THEN

CL    TEST FOR NEGATIVE PRESSURE VALUES.

! loop over all points, including valid halos
        DO 310 I=FIRST_VALID_PT,LAST_P_VALID_PT
          IF(PSTAR(I).LT.0.) THEN
            ERROR_CODE = 1
            WRITE(6,*)' NEGATIVE PRESSURE AT POINT ',I
            WRITE(6,*)' ON PROCESSOR ',MY_PROC_ID
          END IF
 310    CONTINUE
      CALL GC_IMAX(1,N_PROCS,info,ERROR_CODE)

        IF(ERROR_CODE.EQ.1)
     *    ERROR_MESSAGE='P_TH_ADJ : NEGATIVE PRESSURE VALUE CREATED.'

      ENDIF
CL    END OF ROUTINE P_TH_ADJ

      RETURN
      END
