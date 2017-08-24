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
CLL   SUBROUTINE MASS_UWT -----------------------------------------
CLL
CLL   PURPOSE:   CALCULATES RS AND REMOVES MASS WEIGHTING FROM THETAL,
CLL              QT U AND V FIELDS.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.4     22/06/94  Argument LLINTS added and passed to CALC_RS
CLL                     DEF NOWHBR replaced by LWHITBROM
CLL                                                  S.J.Swarbrick
!     3.5    28/03/95 MPP code: Take account of duff row
!                     of U_FIELD data at bottom.
!                                               P.Burton
!     4.1    02/04/96 Tidied up MPP code  P.Burton
CLL
CLL   4.4    14/07/97 RS output via arg list for use in FILTUV.
CLL                   A. Dickinson
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                       STANDARD B. VERSION 2, DATED 18/01/90
CLL
CLL   SYSTEM COMPONENTS COVERED: P12 (part)
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:        NIL.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE MASS_UWT
     1                   (RS_SQUARED_DELTAP,RS,THETAL,QT,U,V,PSTAR,AK,
     2                    BK,DELTA_AK,DELTA_BK,P_FIELD,U_FIELD,P_LEVELS,
     3                    Q_LEVELS,ROW_LENGTH,
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
     4                    LLINTS,LWHITBROM)

      IMPLICIT NONE
      LOGICAL  LLINTS, LWHITBROM

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS ON A ROW.
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
     * QT(P_FIELD,Q_LEVELS)           !INOUT. QT FIELD.
     *,THETAL(P_FIELD,P_LEVELS)       !INOUT THETAL FIELD.
     *,U(U_FIELD,P_LEVELS)            !INOUT U FIELD.
     *,V(U_FIELD,P_LEVELS)            !INOUT U FIELD.

      REAL
     * RS_SQUARED_DELTAP(P_FIELD,P_LEVELS) !OUT HOLDS RS*RS*DELTA P
     *,RS(P_FIELD,P_LEVELS)                !OUT HOLDS RS

      REAL
     * PSTAR(P_FIELD)                 !IN SURFACE PRESSURE
     *,AK(P_LEVELS)                   !IN FIRST TERM IN HYBRID CO-ORDS
     *,BK(P_LEVELS)                   !IN SECOND TERM IN HYBRID CO-ORDS
     *,DELTA_AK(P_LEVELS)             !IN LAYER THICKNESS TERM
     *,DELTA_BK(P_LEVELS)             !IN LAYER THICKNESS TERM

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 1 IS REQUIRED

      REAL
     *  WORK1(P_FIELD)  !GENERAL WORKSPACE
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,K,ROWS,POINTS

      REAL
     * SCALAR
C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------

      EXTERNAL
     * P_TO_UV
     * ,CALC_RS

C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------


C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
      POINTS=LAST_P_VALID_PT-FIRST_VALID_PT+1
! Number of points to be processed by CALC_RS. For non-MPP runs this
! is simply P_FIELD, for MPP, it is all the points, minus any
! unused halo areas (ie. the halo above North pole row, and beneath
! South pole row)

      IF (LWHITBROM) THEN
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALL CALC_RS TO CALCULATE RS.
CL---------------------------------------------------------------------

CL    CALL CALC_RS TO GET RS FOR LEVEL 1.
C TS IS RETURNED IN WORK1, RS AT LEVEL K-1 IS INPUT IN
C RS( ,2) AS AT K-1= 0 THE INPUT IS NOT USED BY CALC_RS.

      K=1
      CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,WORK1(FIRST_VALID_PT),
     &             RS(FIRST_VALID_PT,2),
     *             RS(FIRST_VALID_PT,K),
     &             POINTS,K,P_LEVELS,LLINTS)

CL LOOP FROM 2 TO P_LEVELS
      DO 100 K= 2,P_LEVELS

CL    CALL CALC_RS TO GET RS FOR LEVEL K.
C TS IS RETURNED IN WORK1, RS AT LEVEL K-1 IS INPUT AS
C RS(K-1).

        CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,WORK1(FIRST_VALID_PT),
     &               RS(FIRST_VALID_PT,K-1),
     &               RS(FIRST_VALID_PT,K),
     &               POINTS,K,P_LEVELS,LLINTS)
 100  CONTINUE

CL END LOOP FROM 2 TO P_LEVELS.

      END IF      !     LWHITBROM
CL
CL---------------------------------------------------------------------
C  IF (.NOT.LWHITBROM) THEN
CL    SECTION 1      CALCULATE A*A*DELTA P AND REMOVE MASS-WEIGHTING
C  ELSE
CL    SECTION 2      CALCULATE RS*RS*DELTA P AND REMOVE MASS-WEIGHTING
C  END IF
CL                   FROM THETAL AND QT.
CL---------------------------------------------------------------------

CL LOOP OVER MOIST LEVELS, IE: Q_LEVELS.
      DO 200 K=1,Q_LEVELS
CFPP$ SELECT(CONCUR)

! loop over all points, including valid halos
        DO 210 I=FIRST_VALID_PT,LAST_P_VALID_PT

      IF (.NOT.LWHITBROM) THEN

CL    CALCULATE A*A*DELTAP
                         RS(I,K) = A
          RS_SQUARED_DELTAP(I,K) = A*A*
     *                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
          SCALAR = 1./RS_SQUARED_DELTAP(I,K)

      ELSE

CL    CALCULATE RS*RS*DELTAP
          RS_SQUARED_DELTAP(I,K) = RS(I,K) *
     *                             RS(I,K) *
     *                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
          SCALAR = 1./RS_SQUARED_DELTAP(I,K)

      END IF

CL    REMOVE MASS-WEIGHTING FROM THETAL AND QT.
          THETAL(I,K) = THETAL(I,K)*SCALAR
          QT(I,K) = QT(I,K)*SCALAR
 210    CONTINUE
CL END LOOP OVER MOIST LEVELS.
 200  CONTINUE

CL LOOP OVER ANY REMAINING DRY LEVELS, IE: Q_LEVELS+1 TO P_LEVELS

      DO 220 K= Q_LEVELS+1, P_LEVELS
CFPP$ SELECT(CONCUR)
! loop over all points, including valid halos
        DO 230 I=FIRST_VALID_PT,LAST_P_VALID_PT

      IF (.NOT.LWHITBROM) THEN

CL    CALCULATE A*A*DELTAP
                         RS(I,K) = A
          RS_SQUARED_DELTAP(I,K) = A*A*
     *                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))

      ELSE

CL    CALCULATE RS*RS*DELTAP
          RS_SQUARED_DELTAP(I,K) = RS(I,K) *
     *                             RS(I,K) *
     *                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))

      END IF

CL    REMOVE MASS-WEIGHTING FROM THETAL.
          THETAL(I,K) = THETAL(I,K)/RS_SQUARED_DELTAP(I,K)
 230    CONTINUE

CL END LOOP OVER REMAINING DRY LEVELS.
 220  CONTINUE

CL
CL---------------------------------------------------------------------
C IF (.NOT.LWHITBROM) THEN
CL    SECTION 2      INTERPOLATE A*A*DELTA P ONTO U GRID AND REMOVE
C ELSE
CL    SECTION 3      INTERPOLATE RS*RS*DELTA P ONTO U GRID AND REMOVE
C END IF
CL                   MASS-WEIGHTING FROM U AND V.
CL---------------------------------------------------------------------

C SET ROWS
      ROWS = P_FIELD/ROW_LENGTH

CL LOOP OVER P_LEVELS
      DO 300 K=1,P_LEVELS

C IF (.NOT.LWHITBROM) THEN
CL    CALL P_TO_UV TO OBTAIN A*A*DELTA P ON U GRID.
C ELSE
CL    CALL P_TO_UV TO OBTAIN RS*RS*DELTA P ON U GRID.
C END IF

        CALL P_TO_UV(RS_SQUARED_DELTAP(1,K),WORK1,P_FIELD,U_FIELD,
     *               ROW_LENGTH,ROWS)

! loop over "local" points - not including top and bottom halos
        DO 310 I=FIRST_FLD_PT,LAST_U_FLD_PT
CL    REMOVE MASS-WEIGHTING FROM U AND V.

          SCALAR = 1./WORK1(I)
          U(I,K) = U(I,K)*SCALAR
          V(I,K) = V(I,K)*SCALAR
 310    CONTINUE

CL END LOOP OVER P_LEVELS.
 300  CONTINUE

CL    END OF ROUTINE MASS_UWT

      RETURN
      END
