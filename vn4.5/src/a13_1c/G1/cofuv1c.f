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
CLL   SUBROUTINE COEFF_UV -----------------------------------------
CLL
CLL   PURPOSE:  CALCULATES EFFECTIVE DIFFUSIVE COEFFICIENTS FOR U AND V
CLL             IN NS AND EW DIRECTIONS
CLL              IF STEEP SLOPE THEN EFFECTIVE DIFFUSION IS ZERO.
CLL
CLL              NOTE PRESSURE ARRAY NEEDS TO BE GLOBAL (SHARED)
CLL              FOR MULTI-TASKING AT 3.4 UPWARDS.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   WAS VERSION FOR CRAY Y-MP
CLL
CLL  MODEL            MODIFICATION HISTORY
CLL VERSION  DATE
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with COFUV1A.
CLL   4.4    25/07/97 Calling sequence changed from once per level
CLL                   to once per dynamics sweep, in
CLL                   order to improve MPP scalability.
CLL                   A. Dickinson
CLL
CLL
CLL   PROGRAMMING STANDARD:
CLL
CLL   SYSTEM COMPONENTS COVERED: P132
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (47)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 16 DATED 09/01/91.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE COEFF_UV
     1                 (DIFFUSION_EW,DIFFUSION_NS,
     2                 PRESSURE,PRESSURE_TEST,AK,BK,
     3                 COS_P_LATITUDE,START_U_UPDATE,
     4                 END_U_UPDATE,ROW_LENGTH,
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
     &                 LATITUDE_STEP_INVERSE,
     5                 LONGITUDE_STEP_INVERSE,P_FIELD,U_FIELD,P_LEVELS,
     6                 KD,DELTA_AK,DELTA_BK,PSTAR,COS_FUNCTION_P)

      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, P_LEVELS           !IN NUMBER OF MODEL LEVELS
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, START_U_UPDATE     !IN FIRST POINT TO BE UPDATED.
     *, END_U_UPDATE       !IN LAST POINT TO BE UPDATED.

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
     * PRESSURE(P_FIELD,P_LEVELS)      !IN.3-D PRESSURE FIELD U POINTS
     *          ! LEVEL_P=1 SURFACE THEN LEVEL_P=K IS LEVEL K-1
     *,DIFFUSION_EW(P_FIELD,P_LEVELS)
            !OUT EFFECTIVE EW DIFFUSION COEFF
     *,DIFFUSION_NS(P_FIELD,P_LEVELS)
           !OUT EFFECTIVE NS DIFFUSION COEFF


      REAL
     * AK(P_LEVELS)                    !IN LAYER AK'S
     *,BK(P_LEVELS)                    !IN LAYER BK'S
     *,DELTA_AK(P_LEVELS)              !IN LAYER DELTA_AK'S
     *,DELTA_BK(P_LEVELS)              !IN LAYER DELTA_BK'S
     *,KD(P_LEVELS)                    !IN DIFFUSION COEFF SEE EQ. (45)
     *,PSTAR(P_FIELD)                  !IN PSTAR
     *,COS_P_LATITUDE(P_FIELD)         !IN COS(LAT) AT P POINTS
     *,COS_FUNCTION_P(P_FIELD)         !IN
     *,LATITUDE_STEP_INVERSE           !IN 1/(DELTA LAMDA)
     *,LONGITUDE_STEP_INVERSE          !IN 1/(DELTA PHI)
     *, PRESSURE_TEST      !IN PRESSURE ALTITUDE LIMIT FOR SLOPE TEST


C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
! Define local arrays
      LOGICAL MASK(P_FIELD) ! Indicates of EW_DIFFUSION to be set to
!                           ! zero at a point
      REAL
     * DIFFUSION_COEFFICIENT(P_FIELD)  !IN HOLD ON P GRID. FIRST POINT
     *                                 ! OF ARRAY IS FIRST P POINT ON
     *                                 ! SECOND P ROW. EAST-WEST
     *                                 ! DIFFUSION COEFFICIENT.
     *,DIFFUSION_COEFFICIENT2(P_FIELD) !IN HOLD ON P GRID. FIRST POINT
     *                                 ! OF ARRAY IS FIRST P POINT ON
     *                                 ! SECOND P ROW. NORTH-SOUTH
     *                                 ! DIFFUSION COEFFICIENT.

C DEFINE LOCAL VARIABLES

C LOCAL REALS.
      REAL
     *  PRESSURE_LEVEL

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,IJ,LEVEL,LEVEL_P
C   LEVEL_P=LEVEL+1 IS FOR PRESSURE TEST
C*L   EXTERNAL SUBROUTINE CALLS: NONE---------------------------------

C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS END_U_UPDATE-START_U_UPDATE+1+
CL                                   ROW_LENGTH
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL

      DO LEVEL=1,P_LEVELS

C SET DIFFUSION COEFFICIENT
        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT
          DIFFUSION_COEFFICIENT2(I) = KD(LEVEL)*
     1        (DELTA_AK(LEVEL)+DELTA_BK(LEVEL)*PSTAR(I))
          DIFFUSION_COEFFICIENT(I) = COS_FUNCTION_P(I)*
     2                         DIFFUSION_COEFFICIENT2(I)
        END DO



CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE FIRST TERM IN EQUATION (47)
CL---------------------------------------------------------------------

C   LEVEL_P=LEVEL+1 IS FOR PRESSURE TEST
      LEVEL_P=LEVEL+1
C----------------------------------------------------------------------
CL    TOP LEVEL LEVEL_P = P_LEVELS SINCE SLOPE TEST NEED NOT BE
CL     DONE FOR TOP MOST (PRESSURE) LEVELS
C----------------------------------------------------------------------
      IF(LEVEL_P.GT.P_LEVELS)LEVEL_P=P_LEVELS
C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE DELTALAMBDA TERMS
C                  DELTAPHIKLAMBDA*1/(DELTALAMBDA)SQUARED
C----------------------------------------------------------------------

      DO I= START_U_UPDATE,END_U_UPDATE
       DIFFUSION_EW(I,LEVEL) = 0.5*(DIFFUSION_COEFFICIENT(I+ROW_LENGTH)
     &            + DIFFUSION_COEFFICIENT(I))*LONGITUDE_STEP_INVERSE
     &             *LONGITUDE_STEP_INVERSE
       END DO



C----------------------------------------------------------------------
CL    SECTION 1.2    SET EFFECTIVE DIFFUSION COEFFICIENT TO ZERO
C                    IF STEEP SLOPE BELOW PRESSURE ALTITUDE LIMIT
C                    APPLY GENERAL TEST AT FIRST POINT ONLY
C----------------------------------------------------------------------

C      APPLY GENERAL TEST FOR REFERENCE SURFACE PRESSURE OF 1000HPA
       PRESSURE_LEVEL=AK(LEVEL)+100000.0*BK(LEVEL)
       IF(PRESSURE_LEVEL.GT.PRESSURE_TEST)THEN

      DO I= START_U_UPDATE+1,END_U_UPDATE
        MASK(I)=((PRESSURE(I-1,LEVEL_P).GT.PRESSURE(I,LEVEL_P-1)).OR.
     &           (PRESSURE(I-1,LEVEL_P).LT.PRESSURE(I,LEVEL_P+1)))
      ENDDO


! And zero appropriate points of EW_DIFFUSION
      DO I= START_U_UPDATE,END_U_UPDATE
        IF (MASK(I)) DIFFUSION_EW(I,LEVEL)=0.0
      ENDDO

       ENDIF


CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE SECOND TERM IN EQUATION (47)
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE DELTAPHI TERMS
CL        CALCULATE DELTALAMBDAK*COSLAT/(DELTAPHI)SQUARED
C----------------------------------------------------------------------

! Loop over field missing Northern row
      DO I=START_POINT_NO_HALO,LAST_U_FLD_PT-1
      DIFFUSION_NS(I,LEVEL)=0.5*(DIFFUSION_COEFFICIENT2(I)
     &           *COS_P_LATITUDE(I)
     &           +DIFFUSION_COEFFICIENT2(I+1)*COS_P_LATITUDE(I+1))*
     &            LATITUDE_STEP_INVERSE*LATITUDE_STEP_INVERSE
      END DO

C  RECALCULATE END POINTS.

      DIFFUSION_NS(LAST_U_FLD_PT,LEVEL)
     &                  =DIFFUSION_NS(LAST_U_FLD_PT-1,LEVEL)




C----------------------------------------------------------------------
CL    SECTION 2.2    SET EFFECTIVE DIFFUSION COEFFICIENT TO ZERO
C                    IF STEEP SLOPE BELOW PRESSURE ALTITUDE LIMIT
C                    APPLY GENERAL TEST AT FIRST POINT ONLY
C----------------------------------------------------------------------

C      APPLY GENERAL TEST FOR REFERENCE SURFACE PRESSURE OF 1000HPA
       IF(PRESSURE_LEVEL.GT.PRESSURE_TEST)THEN

! Loop over field, missing Northern row
      DO I=START_POINT_NO_HALO,LAST_U_FLD_PT
      IF((PRESSURE(I,LEVEL_P).GT.PRESSURE(I-ROW_LENGTH,LEVEL_P-1)).OR.
     &   (PRESSURE(I,LEVEL_P).LT.
     &        PRESSURE(I-ROW_LENGTH,LEVEL_P+1)))THEN
         DIFFUSION_NS(I,LEVEL)=0.0
       ENDIF

      END DO

      ENDIF
      ENDDO

      CALL SWAPBOUNDS(DIFFUSION_EW,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)
      CALL SWAPBOUNDS(DIFFUSION_NS,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,P_LEVELS)

CL    END OF ROUTINE COEFF_UV

      RETURN
      END
