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
CLL   SUBROUTINE FILT_UV ------------------------------------------
CLL
CLL   PURPOSE:   PERFORMS MASS-WEIGHTED FILTERING AND POLAR AVERAGING OF
CLL              U AND V FIELDS.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with FILTUV1A.
CLL   4.4    14/07/97  Simplify calculation of RS*DELTAP for efficiency.
CLL                    A.Dickinson
CLL
CLL   PROGRAMMING STANDARD:
CLL
CLL   SYSTEM COMPONENTS COVERED: P142
CLL   SYSTEM TASK: P1
CLL   DOCUMENTATION:       SECTION 3.5
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 8, DATED 10/09/90.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE FILT_UV
     1                  (PSTAR,U,V,RS_FUNCTIONS,DELTA_AK,DELTA_BK,
     2                   P_FIELD,U_FIELD,NORTHERN_FILTERED_P_ROW,
     3                   SOUTHERN_FILTERED_P_ROW,P_LEVELS,
     4                   ROW_LENGTH,
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
     &                   TRIGS,IFAX,
     4                   COS_LONGITUDE,SIN_LONGITUDE,
     5                   FILTER_WAVE_NUMBER_U_ROWS)

      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, P_LEVELS           !IN NUMBER OF MODEL LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *                          ! MOVING TOWARDS EQUATOR.
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN
     *                          ! MOVING TOWARDS SOUTH POLE.
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.

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
     &  FILTER_WAVE_NUMBER_U_ROWS(GLOBAL_U_FIELD/GLOBAL_ROW_LENGTH)
!       LAST WAVE NUMBER NOT TO BE CHOPPED
      REAL
     * U(U_FIELD,P_LEVELS) !INOUT U VELOCITY FIELD.
     *,V(U_FIELD,P_LEVELS) !INOUT V VELOCITY FIELD.

      REAL
     * PSTAR(P_FIELD)                 !IN PSTAR FIELD.
     *,RS_FUNCTIONS(P_FIELD,P_LEVELS) !IN RS
     *,DELTA_AK(P_LEVELS)             !IN LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)             !IN LAYER THICKNESS
     *,TRIGS(ROW_LENGTH)              !IN HOLDS TRIGONOMETRIC FUNCTIONS
     *                                ! USED IN FILTERING.
     *,COS_LONGITUDE(ROW_LENGTH)      !IN COSINE LONGITUDE AT U POINTS
     *,SIN_LONGITUDE(ROW_LENGTH)      !IN SINE   LONGITUDE AT U POINTS
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 1 IS REQUIRED
      REAL
     * WORK(P_FIELD)      ! GENERAL WORKSPACE.

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C REAL SCALARS
      REAL
     *  SCALAR3

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,K
     *, NORTHERN_FILTERED_U_ROW ! U ROW ON WHICH FILTERING STOPS
     *, SOUTHERN_FILTERED_U_ROW ! U ROW ON WHICH FILTERING STARTS AGAIN
     *, FILTER_SPACE ! HORIZONTAL DIMENSION OF SPACE NEEDED IN FILTERING
     *               ! ROUTINE.

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL
     * P_TO_UV,FILTER,POLAR_UV
C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     MASS-WEIGHT U AND V FIELDS.
CL---------------------------------------------------------------------

! QAN fix : blank out WORK array so halos don't contain junk
      DO I=1,P_FIELD
        WORK(I)=0.0
      ENDDO

CL LOOP OVER P_LEVELS.

      DO 100 K=1,P_LEVELS

CL    CALCULATE RS*DELTA P
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          WORK(I)   =   RS_FUNCTIONS(I,K)
     &              * (DELTA_AK(K) + DELTA_BK(K)*PSTAR(I))
        END DO


CL    CALL P_TO_UV TO TRANSFER RS*DELTA P TO U GRID.

        CALL P_TO_UV(WORK,RS_FUNCTIONS(1,K),P_FIELD,U_FIELD,ROW_LENGTH,
     &                tot_P_ROWS)

CL    MASS WEIGHT U AND V FIELDS.

        DO 120 I=FIRST_VALID_PT,LAST_U_VALID_PT
          U(I,K) = U(I,K) * RS_FUNCTIONS(I,K)
          V(I,K) = V(I,K) * RS_FUNCTIONS(I,K)
 120    CONTINUE

CL END LOOP OVER P_LEVELS.
 100  CONTINUE

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALL FILTER TO FILTER FIELDS.
CL---------------------------------------------------------------------

C SET NORTHERN AND SOUTHERN ROWS FOR U FILTERING.

      NORTHERN_FILTERED_U_ROW = NORTHERN_FILTERED_P_ROW
      SOUTHERN_FILTERED_U_ROW = SOUTHERN_FILTERED_P_ROW - 1

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE = (ROW_LENGTH+2)*(NORTHERN_FILTERED_U_ROW-1+
     *                U_FIELD/ROW_LENGTH-SOUTHERN_FILTERED_U_ROW)

CL    CALL FILTER FOR U

      CALL FILTER(U,U_FIELD,P_LEVELS,FILTER_SPACE,ROW_LENGTH,
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
     *            FILTER_WAVE_NUMBER_U_ROWS,TRIGS,IFAX,
     *            NORTHERN_FILTERED_U_ROW,SOUTHERN_FILTERED_U_ROW)

CL    CALL FILTER FOR V

      CALL FILTER(V,U_FIELD,P_LEVELS,FILTER_SPACE,ROW_LENGTH,
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
     *            FILTER_WAVE_NUMBER_U_ROWS,TRIGS,IFAX,
     *            NORTHERN_FILTERED_U_ROW,SOUTHERN_FILTERED_U_ROW)

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     REMOVE MASS-WEIGHTING FROM U AND V.
CL---------------------------------------------------------------------


! CALL POLAR_UV TO UPDATE THE POLAR VALUES OF MASS-WEIGHTED U AND V

      CALL POLAR_UV(U,V,ROW_LENGTH,U_FIELD,P_LEVELS,
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
     &              COS_LONGITUDE,SIN_LONGITUDE)
CL LOOP OVER P_LEVELS.

      DO 300 K=1,P_LEVELS


CL    REMOVE MASS-WEIGHTING.

CFPP$ SELECT(CONCUR)
        DO 310 I=FIRST_VALID_PT,LAST_U_VALID_PT
          SCALAR3 = 1./RS_FUNCTIONS(I,K)
          U(I,K) = U(I,K) * SCALAR3
          V(I,K) = V(I,K) * SCALAR3
 310    CONTINUE

CL END LOOP OVER P_LEVELS

 300  CONTINUE

CL    END OF ROUTINE FILT_UV

      RETURN
      END
