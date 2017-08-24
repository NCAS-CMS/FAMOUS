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
CLL   SUBROUTINE CALC_RS ---------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES RS AS A FUNCTION OF PRESSURE USING
CLL              EQUATION (17) AND THE U.M. STANDARD ATMOSPHERE.
CLL              ALSO RETURNS U.M. STANDARD TEMPERATURE AT THE
CLL              INPUT PRESSURE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL   NOT SUITABLE FOR I.B.M. USE.
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.4     26/05/94  Argument LLINTS added and passed to CALC_TS
CLL                                                     S.J.Swarbrick
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD A. VERSION 2, DATED 18/01/90
CLL
CLL   LOGICAL COMPONENTS COVERED: P194
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (17)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN, T.DAVIES AND M.H.MAWSON,
CLL                        VERSION 9, DATED 27/06/90.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE CALC_RS
     1  (PSTAR,AK,BK,TS,RS_LOWER,RS,POINTS,LEVEL_REQUESTED,LEVELS,
     2   LLINTS)

      IMPLICIT NONE
      LOGICAL  LLINTS  ! Logical switch for linear TS in CALC_TS

      INTEGER
     *  POINTS            !IN. NUMBER OF POINTS OVER WHICH CALCULATION
     *                    !IS TO BE PERFORMED.
     *, LEVEL_REQUESTED   !IN. MODEL LEVEL AT WHICH ROUTINE IS BEING
     *                    !PERFORMED.
     *, LEVELS            !IN. NUMBER OF MODEL LEVELS

      REAL
     *  PSTAR(POINTS) !IN.   SURFACE PRESSURE VALUES.
     *, RS_LOWER(POINTS) !IN. HOLDS RS VALUES ONE LEVEL LOWER THAN
     *                ! REQUESTED IN CALL. IF CALLED FROM LEVEL 1
     *                ! THEN HOLDS DUMMY VALUES AND IS NOT USED.
     *, AK(LEVELS)    !IN. HOLDS AK VALUES AT PRESSURE LEVELS.
     *, BK(LEVELS)    !IN. HOLDS BK VALUES AT PRESSURE LEVELS.

      REAL
     *  RS(POINTS)    !OUT.  RS VALUES. NOTE THIS ARRAY HAS NO LEVELS
     *                !UNLIKE RS IN OTHER ROUTINES. THE RETURNED VALUES
     *                !ARE STORED IN THE LEVEL DETERMINED BY THE CALL.

      REAL
     *  TS(POINTS)    !INOUT.  U.M. STANDARD TEMPERATURE AT PRESSURE P.
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
CL    3 LOCAL ARRAYS NEEDED.
      REAL
     *  P_LEVEL(POINTS)        ! PRESSURE AT LEVEL_REQUESTED.
     *, P_LEVEL_MINUS1(POINTS) ! PRESSURE AT INPUT_REQUESTED MINUS 1
     *, TS_LEVEL(POINTS)       ! TS AT LEVEL_REQUESTED.
C*---------------------------------------------------------------------
C REAL SCALARS
      REAL TS0_BY_P0
C COUNT VARIABLES FOR DO LOOPS
      INTEGER
     *  I
C LOGICAL VARIABLE
      LOGICAL
     *  CONSTANT_PRESSURE      ! SET TO TRUE IF LEVEL_REQUIRED IS A
     *                         ! CONSTANT PRESSURE LEVEL.

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL CALC_TS
C*---------------------------------------------------------------------
CL    CALL COMDECK TO OBTAIN CONSTANTS USED.

CLL   COMDECK C_CALCRS HOLDS CONSTANTS FOR ROUTINE CALC_RS.
CLL   COMDECK CONTAINS SOME LOCAL CONSTANTS.
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

      REAL P0,TS0,HALF_R_OVER_G
      PARAMETER(P0=101325.,TS0=288.15,HALF_R_OVER_G=.5*R/G)
CL    END OF COMDECK C_CALCRS

CL    MAXIMUM VECTOR LENGTH IS DETERMINED BY POINTS.
CL
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
CL    ON A CALL TO CALC_TS EITHER SECTION 1 OR SECTION 2 IS USED
CL    DEPENDING ON WHETHER LEVEL_REQUESTED IS 1 OR NOT.
CL
CL---------------------------------------------------------------------

CL    CHECK TO SEE IF LEVEL_REQUESTED IS A CONSTANT PRESSURE LEVEL.

      IF(BK(LEVEL_REQUESTED).EQ.0.) THEN
        CONSTANT_PRESSURE= .TRUE.
      ELSE
        CONSTANT_PRESSURE= .FALSE.
      END IF
CL---------------------------------------------------------------------
CL    SECTION 1.  IF LEVEL_REQUESTED IS 1 CALCULATE RS USING
CL                EQUATION 17 WITH THE VALUES AT P0 BEING TAKEN FROM
CL                THE COMDECK.
CL    EITHER A) IF CONSTANT_PRESSURE USE SECTIONS 1.1 TO 1.3
CL    OR     B) IF NOT CONSTANT_PRESSURE THEN USE SECTIONS 1.4 - 1.6
CL
CL    NOTE THIS SECTION IS NOT VERY ACCURATE IF RUNNING
CL    STRATOSPHERIC MODEL.
CL---------------------------------------------------------------------

CL A)
      IF(LEVEL_REQUESTED.EQ.1) THEN

        IF(CONSTANT_PRESSURE) THEN
C----------------------------------------------------------------------
CL    SECTION 1.1 STORE PRESSURE AT LEVEL 1 IN P_LEVEL(1) ONLY.
CL             THIS IS BECAUSE NO OTHER ADDRESSES IN P_LEVEL WILL BE
CL             ACCESSED.
C----------------------------------------------------------------------

          P_LEVEL(1) = AK(1)

C----------------------------------------------------------------------
CL    SECTION 1.2 CALL CALC_TS TO OBTAIN TS AT LEVEL 1.
C----------------------------------------------------------------------

          CALL CALC_TS(P_LEVEL,TS,POINTS,CONSTANT_PRESSURE,LLINTS)

C----------------------------------------------------------------------
CL    SECTION 1.3 CALCULATE RS AT PRESSURE P USING EQUATION 17.
C----------------------------------------------------------------------

CL                 RS IS CALCULATED FOR POINT 1 AND THEN RS(I) FOR
CL                 I=2 TO POINTS IS SET EQUAL TO RS(1)
          RS(1)=A + HALF_R_OVER_G*(TS0 / P0 + TS(1)/P_LEVEL(1))
     *           *(P0-P_LEVEL(1))
          DO 130 I=2,POINTS
            RS(I) = RS(1)
 130      CONTINUE

CL B)
        ELSE
C----------------------------------------------------------------------
CL    SECTION 1.4 CALCULATE PRESSURE AT LEVEL 1.
C----------------------------------------------------------------------
          DO 140 I=1,POINTS
            P_LEVEL(I) = AK(1)+PSTAR(I)*BK(1)
140       CONTINUE

C----------------------------------------------------------------------
CL    SECTION 1.5 CALL CALC_TS TO OBTAIN TS AT LEVEL 1.
C----------------------------------------------------------------------

          CALL CALC_TS(P_LEVEL,TS,POINTS,CONSTANT_PRESSURE,LLINTS)

C----------------------------------------------------------------------
CL    SECTION 1.6 CALCULATE RS AT PRESSURE P USING EQUATION 17.
C----------------------------------------------------------------------

          TS0_BY_P0 = TS0 / P0
          DO 160 I=1,POINTS
            RS(I)=A + HALF_R_OVER_G*( TS0_BY_P0+ TS(I)/P_LEVEL(I))
     *           *(P0-P_LEVEL(I))
 160      CONTINUE

        ENDIF

      ELSE
CL---------------------------------------------------------------------
CL    SECTION 2.  IF LEVEL_REQUESTED IS NOT 1 THEN CALCULATE
CL                INCREMENT TO RS BETWEEN PRESSURE AT LEVEL_REQUESTED
CL                MINUS 1 AND LEVEL_REQUESTED USING EQUATION 17 AND
CL                ADD ON TO RS AT LEVEL_REQUESTED MINUS 1.
CL
CL    EITHER A) IF NOT CONSTANT_PRESSURE USE SECTIONS 2.1 TO 2.3.
CL    OR     B) IF CONSTANT_PRESSURE BUT LEVEL_REQUESTED - 1 IS
CL              NOT CONSTANT PRESSURE THEN USE SECTIONS 2.4 TO 2.6
CL    OR     C) BOTH CONSTANT PRESSURE USE SECTIONS 2.7 TO 2.9.
CL---------------------------------------------------------------------

CL A)
        IF(.NOT.CONSTANT_PRESSURE) THEN
C----------------------------------------------------------------------
CL    SECTION 2.1 CALCULATE PRESSURE AT LEVEL_REQUESTED AND
CL                LEVEL_REQUESTED MINUS 1.
C----------------------------------------------------------------------

          DO 210 I=1,POINTS
            P_LEVEL(I) = AK(LEVEL_REQUESTED) + BK(LEVEL_REQUESTED)*
     *                   PSTAR(I)
            P_LEVEL_MINUS1(I) = AK(LEVEL_REQUESTED-1) +
     *                   BK(LEVEL_REQUESTED-1)*PSTAR(I)
 210      CONTINUE

C----------------------------------------------------------------------
CL    SECTION 2.2 CALL CALC_TS TO GET TS AT LEVEL_REQUESTED AND STORE
CL                IN TS_LEVEL.
C----------------------------------------------------------------------

          CALL CALC_TS(P_LEVEL,TS_LEVEL,POINTS,CONSTANT_PRESSURE,
     &                 LLINTS)

C----------------------------------------------------------------------
CL    SECTION 2.3 CALCULATE INTEGRAL IN EQUATION 17 BETWEEN THE
CL                PRESSURES CALCULATED IN 2.1 AND ADD ONTO RS AT
CL                LEVEL_REQUESTED MINUS 1 TO GET VALUE AT
CL                LEVEL_REQUESTED. THEN OVER-WRITE OLD VALUE OF TS
CL                WITH VALUE CALCULATED IN CALL TO TS_CALC IN 2.2.
C----------------------------------------------------------------------

          DO 230 I=1,POINTS
            RS(I) = RS_LOWER(I) + (P_LEVEL_MINUS1(I)-P_LEVEL(I))
     *              *HALF_R_OVER_G*(TS(I)/P_LEVEL_MINUS1(I) +
     *                 TS_LEVEL(I)/P_LEVEL(I))
            TS(I) = TS_LEVEL(I)
 230      CONTINUE

        ELSE IF(BK(LEVEL_REQUESTED-1).NE.0.) THEN

CL B)
C----------------------------------------------------------------------
CL    SECTION 2.4 CALCULATE PRESSURE AT LEVEL_REQUESTED-1 AT ALL POINTS
CL                PRESSURE AT LEVEL_REQUESTED SET TO AK AT POINT 1 ONLY
CL                AS NO OTHER ADDRESSES ARE USED.
C----------------------------------------------------------------------

          P_LEVEL(1) = AK(LEVEL_REQUESTED)
          DO 240 I=1,POINTS
            P_LEVEL_MINUS1(I) = AK(LEVEL_REQUESTED-1) +
     *                   BK(LEVEL_REQUESTED-1)*PSTAR(I)
 240      CONTINUE

C----------------------------------------------------------------------
CL    SECTION 2.5 CALL CALC_TS TO GET TS AT LEVEL_REQUESTED AND STORE
CL                IN TS_LEVEL.
C----------------------------------------------------------------------

          CALL CALC_TS(P_LEVEL,TS_LEVEL,POINTS,CONSTANT_PRESSURE,
     &                 LLINTS)

C----------------------------------------------------------------------
CL    SECTION 2.6 CALCULATE INTEGRAL IN EQUATION 17 BETWEEN THE
CL                PRESSURES CALCULATED IN 2.4 AND ADD ONTO RS AT
CL                LEVEL_REQUESTED MINUS 1 TO GET VALUE AT
CL                LEVEL_REQUESTED. THEN OVER-WRITE OLD VALUE OF TS
CL                WITH VALUE CALCULATED IN CALL TO TS_CALC IN 2.5.
C----------------------------------------------------------------------

          TS0_BY_P0 = TS_LEVEL(1) / P_LEVEL(1)
          DO 260 I=1,POINTS
            RS(I) = RS_LOWER(I) + (P_LEVEL_MINUS1(I)-P_LEVEL(1))
     *              *HALF_R_OVER_G*(TS(I)/P_LEVEL_MINUS1(I) +
     *                TS0_BY_P0)
            TS(I) = TS_LEVEL(I)
 260      CONTINUE

        ELSE

CL C)
C----------------------------------------------------------------------
CL    SECTION 2.7 SET PRESSURE AT LEVEL_REQUESTED-1 TO
CL                AK(LEVEL_REQUESTED-1) AND PRESSURE AT LEVEL_REQUESTED
CL                TO AK(LEVEL_REQUESTED) AT POINT 1 ONLY
CL                AS NO OTHER ADDRESSES ARE USED.
C----------------------------------------------------------------------

          P_LEVEL(1) = AK(LEVEL_REQUESTED)
          P_LEVEL_MINUS1(1) = AK(LEVEL_REQUESTED-1)

C----------------------------------------------------------------------
CL    SECTION 2.8 CALL CALC_TS TO GET TS AT LEVEL_REQUESTED AND STORE
CL                IN TS_LEVEL.
C----------------------------------------------------------------------

          CALL CALC_TS(P_LEVEL,TS_LEVEL,POINTS,CONSTANT_PRESSURE,
     &                 LLINTS)

C----------------------------------------------------------------------
CL    SECTION 2.9 CALCULATE INTEGRAL IN EQUATION 17 BETWEEN THE
CL                PRESSURES CALCULATED IN 2.7 AND ADD ONTO RS AT
CL                LEVEL_REQUESTED MINUS 1 TO GET VALUE AT
CL                LEVEL_REQUESTED. THEN OVER-WRITE OLD VALUE OF TS
CL                WITH VALUE CALCULATED IN CALL TO TS_CALC IN 2.5.
C----------------------------------------------------------------------

          TS0_BY_P0 = HALF_R_OVER_G*(P_LEVEL_MINUS1(1)-P_LEVEL(1))
     *                   *(TS(1)/P_LEVEL_MINUS1(1) +
     *                     TS_LEVEL(1) / P_LEVEL(1))
          DO 290 I=1,POINTS
            RS(I) = RS_LOWER(I) + TS0_BY_P0
            TS(I) = TS_LEVEL(I)
 290      CONTINUE

        END IF

      END IF

CL    END OF ROUTINE CALC_RS

      RETURN
      END

