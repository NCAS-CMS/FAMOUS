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
CLL  SUBROUTINE TESTDIAG------------------------------------------------
CLL
CLL  PURPOSE: CALCULATE SIMPLE TEST DIAGNOSTICS BASED ON A SIMPLE
CLL           ANALYTIC FORMULA:
CLL  VALUE=A*(LATITUDE+90.)+B*LONGITUDE+C*LEVEL+D*FORECAST_HRS
CLL  WHERE A=1.0, B=1.0E2, C=1.0E3, D=1.0E4
CLL  AND (LAT,LONG) ARE IN DEGREES, ACTUAL POSITION (ROTATED FOR LAM),
CLL  LEVEL IS EITHER MODEL LEVEL OR
CLL  PRESSURE LEVEL IN MB., AND FORECAST_HRS IN HOURS.
CLL  THESE DIAGNOSTICS ARE TO BE USED FOR CHECKING OUTPUT PROCEDURES
CLL  AFTER VARIOUS POST-PROCESSING ROUTES.
CLL  FOUR DIAGNOSTICS CAN BE CALCULATED:
CLL  1. SINGLE LEVEL FIELD (LEVEL=0.) AT U POINTS.
CLL  2. SINGLE LEVEL FIELD (LEVEL=0.) AT P POINTS.
CLL  3. MULTI- LEVEL FIELD (LEVEL=PRESS LEVEL) AT P POINTS
CLL  4. MULTI- LEVEL FIELD (LEVEL=MODEL LEVEL) AT P POINTS
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN USE
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   3.1   25/01/93  NEW DECK                             R.RAWLINS
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 2, DATED 18/01/90
CLL
CLL  SYSTEM TASK: D4
CLL
CLL  LOGICAL COMPONENT:
CLL
CLL  DOCUMENTATION:  UMDP NO. D7
CLL
CLLEND------------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE TESTDIAG
     1 (P_FIELD,U_FIELD,P_ROWS,U_ROWS,ROW_LENGTH,EW_SPACE,NS_SPACE,
     2  FIRST_LAT,FIRST_LONG,ELF,PHI_POLE,LAMBDA_POLE,
     3  PRESS_LEVELS_LIST,NO_PRESS_LEVELS,
     4  MODEL_LEVELS_LIST,NO_MODEL_LEVELS,FORECAST_HRS,
     5  DIAG1,DIAG2,DIAG3,DIAG4,
     6  QDIA1,QDIA2,QDIA3,QDIA4)
C
      IMPLICIT NONE
C
      INTEGER
     *  P_FIELD      !IN   FIRST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD      !IN   FIRST DIMENSION OF (U,V) FIELD
     *, P_ROWS       !IN   NO. OF ROWS FOR P FIELD
     *, U_ROWS       !IN   NO. OF ROWS FOR U FIELD
     *, ROW_LENGTH   !IN   NO. OF POINTS PER ROW
     *, NO_MODEL_LEVELS                    !IN MODEL LEVELS FOR OUTPUT
     *, NO_PRESS_LEVELS                    !IN PRESS LEVELS FOR OUTPUT
     *, FORECAST_HRS                       !IN FORECAST HOURS T+0, ETC
C
      LOGICAL
     *  ELF              !IN  TRUE IF MODEL IS LAM WITH ROTATED GRID
     * ,QDIA1            !IN  STASHFLAG FOR DIAG1
     * ,QDIA2            !IN  STASHFLAG FOR DIAG2
     * ,QDIA3            !IN  STASHFLAG FOR DIAG3
     * ,QDIA4            !IN  STASHFLAG FOR DIAG4
C
      REAL
     *  EW_SPACE         !IN  DELTA LONGITUDE (DEGREES)
     *, NS_SPACE         !IN  DELTA  LATITUDE (DEGREES)
     *, FIRST_LAT        !IN  LATITUDE  OF FIRST P ROW IN DEGREES
     *, FIRST_LONG       !IN  LONGITUDE OF FIRST P COL IN DEGREES
     *, PHI_POLE         !IN  LATITUDE OF THE PSEUDO POLE
     *, LAMBDA_POLE      !IN  LONGITUDE OF THE PSEUDO POLE
     *, MODEL_LEVELS_LIST(NO_MODEL_LEVELS) !IN LEVELS LIST (FOR DIAG3)
     *, PRESS_LEVELS_LIST(NO_PRESS_LEVELS) !IN LEVELS LIST (FOR DIAG4)
     *, DIAG1(U_FIELD)   !OUT DIAGNOSTIC 1
     *, DIAG2(P_FIELD)   !OUT DIAGNOSTIC 2
     *, DIAG3(P_FIELD,NO_PRESS_LEVELS)  !OUT DIAGNOSTIC 3
     *, DIAG4(P_FIELD,NO_MODEL_LEVELS)  !OUT DIAGNOSTIC 4
C
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
C*----------------------------------------------------------------------
      REAL
     *  LATITUDE(P_FIELD)    ! LATITUDE IN DEGREES
     *, LONGITUDE(P_FIELD)   ! LONGITUDE IN DEGREES
     *, LAT(P_FIELD)         ! LATITUDE  IN DEGREES ON EQUATORIAL GRID
     *, LONG(P_FIELD)        ! LONGITUDE IN DEGREES ON EQUATORIAL GRID
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL LLTOEQ
C*----------------------------------------------------------------------
C
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

C
C-----------------------------------------------------------------------
C DEFINE LOCAL CONSTANTS
C-----------------------------------------------------------------------
      REAL
     *  A,B,C,D  ! COEFFICIENTS FOR CALCULATING VALUES OF FIELD
C
      PARAMETER(A=1.0,B=1.0E2,C=1.0E3,D=1.0E4)
C
C-----------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  I,J,K        !  LOOP COUNTERS
     * ,L            !  LOOP INDEX
C-----------------------------------------------------------------------
CL  1. CALCULATE FIRST DIAGNOSTIC  (U GRID SINGLE LEVEL)
C-----------------------------------------------------------------------
      IF(QDIA1) THEN
CL
CL  1a. FIND EQUATORIAL LATITUDES,LONGITUDES
CL
         DO J=1,U_ROWS
           DO I=1,ROW_LENGTH
             L= I + (J-1)*ROW_LENGTH
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J-0.5)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I-0.5)
           ENDDO
         ENDDO
CL
CL  1b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
CL
         IF(ELF) THEN
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,
     *                 U_FIELD)
         ELSE
           DO I=1,U_FIELD
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
CL
CL  1c. CALCULATE VALUE FROM ANALYTIC FUNCTION
CL

         DO I=1,U_FIELD
           DIAG1(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +
     *              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF QDIA1 TEST

C-----------------------------------------------------------------------
CL  2. CALCULATE ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
C-----------------------------------------------------------------------
      IF(QDIA2.OR.QDIA3.OR.QDIA4) THEN
CL
CL  2a. FIND EQUATORIAL LATITUDES,LONGITUDES
CL
         DO J=1,P_ROWS
           DO I=1,ROW_LENGTH
             L= I + (J-1)*ROW_LENGTH
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J-1)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I-1)
           ENDDO
         ENDDO
CL
CL  2b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
CL
         IF(ELF) THEN
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,
     *                 P_FIELD)
         ELSE
           DO I=1,P_FIELD
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
      ENDIF                      ! END OF QDIA2-4 TEST
C-----------------------------------------------------------------------
CL  3. CALCULATE SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
C-----------------------------------------------------------------------
      IF(QDIA2) THEN

         DO I=1,P_FIELD
           DIAG2(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +
     *              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF QDIA2 TEST
C-----------------------------------------------------------------------
CL  4. CALCULATE THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
C-----------------------------------------------------------------------
      IF(QDIA3) THEN

         DO K=1,NO_PRESS_LEVELS
           DO I=1,P_FIELD
             DIAG3(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +
     *                  C*PRESS_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF QDIA3 TEST
C-----------------------------------------------------------------------
CL  5. CALCULATE FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
C-----------------------------------------------------------------------
      IF(QDIA4) THEN

         DO K=1,NO_MODEL_LEVELS
           DO I=1,P_FIELD
             DIAG4(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +
     *                  C*MODEL_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF QDIA4 TEST
C=======================================================================
C END OF TESTDIAG
C=======================================================================
      RETURN
      END
C=======================================================================
