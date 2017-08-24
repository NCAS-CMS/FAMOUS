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
CLL   SUBROUTINE CALC_TS ---------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES TEMPERATURE AS A FUNCTION OF PRESSURE USING
CLL              THE UNIFIED MODEL STANDARD ATMOSPHERE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL   NOT SUITABLE FOR I.B.M. USE.
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   3.4  13/04/94   DEF LINEARTS REPLACED BY LOGICAL LLINTS
CLL                                                   S.J.SWARBRICK
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD A. VERSION 2, DATED 18/01/90
CLL
CLL   LOGICAL COMPONENTS COVERED: P195
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION: APPENDIX 2 OF DOCUMENTATION PAPER 10.
CLL                  BY M.J.P.CULLEN, T.DAVIES AND M.H.MAWSON,
CLL                  VERSION 9 ,DATED 27/06/90
CLLEND-------------------------------------------------------------
C
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE CALC_TS
     1  (P,TS,POINTS,CONSTANT_PRESSURE,LLINTS)

      IMPLICIT NONE
      LOGICAL  LLINTS  ! LOGICAL SWITCH FOR LINEAR TS CALC

      INTEGER
     *  POINTS           !IN NUMBER OF POINTS OVER WHICH CALCULATION
     *                   ! IS TO BE PERFORMED.

      REAL
     * P(POINTS)         !IN    PRESSURE VALUES.

      REAL
     * TS(POINTS)        !OUT U.M. STANDARD TEMPERATURE AT PRESSURE P.

      LOGICAL
     * CONSTANT_PRESSURE !IN. IF TRUE THEN P CONTAINS THE SAME VALUES
     *                   ! FOR ALL POINTS.
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
CL    NO LOCAL ARRAYS NEEDED.
C*---------------------------------------------------------------------
C REAL SCALARS
      REAL
     *  TERM1,TERM2,TERM3,TERM4,TERM5,TERM6,CONST
C COUNT VARIABLES FOR DO LOOPS
      INTEGER
     *  I

C*L   NO EXTERNAL SUBROUTINE CALLS.------------------------------------
C*---------------------------------------------------------------------
CL    CALL COMDECK TO OBTAIN CONSTANTS USED.

CLL  COMDECK C_CALCTS HOLDS CONSTANTS FOR ROUTINE CALC_TS.
CLL  COMDECK CONTAINS LOCAL CONSTANTS.
C
CC *IF DEF,LINEARTS  (DEF LINEARTS REPLACED BY LOGICAL LLINTS)
C
      REAL P0,TS0,AS0,BS0
C
C TS0 IS TEMPERATURE AT REFERENCE SURFACE PRESSURE P0.
C AS0 AND BS0 ARE THE CONSTANTS FOR THE LINEAR FIT IN THE LAYER ABOVE.
C
      PARAMETER(P0=101325,TS0=288.15,AS0=9.0859E-4,BS0=196.09)
C
CC *ELSE
C
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

C
C     REAL P0,TS0,L0
      REAL L0
C
C L0 IS LAPSE RATE IN TROPOSPHERE
      REAL P_ISO,T_ISO
C P_ISO IS THE PRESSURE AT WHICH THE ATMOSPHERE BECOMES ISOTHERMAL.
C T_ISO IS THE TEMPERATURE IN THE ISOTHERMAL LAYER.
      REAL P_LOW_STRAT,L_LOW_STRAT
C P_LOW_STRAT IS THE PRESSURE AT WHICH THE ATMOSPHERIC TEMPERATURE
C STARTS TO INCREASE WITH DECREASING PRESSURE.
C L_LOW_STRAT IS THE LAPSE RATE AT PRESSURES BELOW P_LOW_STRAT.
      REAL P_MID_STRAT,T_MID_STRAT,L_MID_STRAT
C P_MID_STRAT IS THE PRESSURE AT WHICH THE ATMOSPHERIC TEMPERATURE
C STARTS TO INCREASE MORE RAPIDLY WITH DECREASING PRESSURE.
C L_MID_STRAT IS THE LAPSE RATE AT PRESSURES BELOW P_MID_STRAT.
C T_MID_STRAT IS THE TEMPERATURE AT P_MID_STRAT.
      REAL P_UPPER_STRAT,T_UPPER_STRAT
C P_UPPER_STRAT IS THE PRESSURE AT WHICH THE ATMOSPHERIC IS ISOTHERMAL
C AGAIN.
C T_UPPER_STRAT IS THE TEMPERATURE AT P_UPPER_STRAT.
      REAL P_MESO,L_MESO
C P_MESO IS THE PRESSURE AT WHICH THE MESOSPHERE BEGINS.
C L_MESO IS THE LAPSE RATE IN THE MESOSPHERE.
      REAL P_MIN,T_MIN
C P_MIN IS THE PRESSURE BELOW WHICH IT IS ASSUMED NO MODEL SHOULD GET.
C T_MIN IS THE TEMPERATURE AT P_MIN.
C
C     PARAMETER(P0=101325,TS0=288.15,L0=-0.0065)
      PARAMETER(                     L0=-0.0065)
C
      PARAMETER(P_ISO=22632,T_ISO = 216.65)
      PARAMETER(P_LOW_STRAT=5475,L_LOW_STRAT=0.001)
      PARAMETER(P_MID_STRAT=868,L_MID_STRAT=0.0028,T_MID_STRAT=228.65)
      PARAMETER(P_UPPER_STRAT=110.9,T_UPPER_STRAT=270.65)
      PARAMETER(P_MESO=74.7,L_MESO=-0.002768)
      PARAMETER(P_MIN =0.0001,T_MIN = 89.309)
C TERMS USED IN TAYLOR SERIES EXPANSIONS. THEIR MEANINGS EXPLAINED
C IN THE SECTION USING THEM.
      REAL A1,A2,A3,A4,RECIP_A1,RECIP_A2,RECIP_A3,RECIP_A4
      REAL RECIP_P0,RECIP_P_LOW_STRAT,RECIP_P_MID_STRAT,RECIP_P_MESO
      REAL MINUS_RLG1,MINUS_RLG2,MINUS_RLG3,MINUS_RLG4
      REAL A1_TO_MINUS_RLG1,A2_TO_MINUS_RLG2,A3_TO_MINUS_RLG3
      REAL A4_TO_MINUS_RLG4
      REAL ONE_THIRD,ONE_SIXTH
      PARAMETER(A1=0.6117,A2=0.5793,A3=0.5638,A4=0.5000)
      PARAMETER(RECIP_A1=1./A1,RECIP_A2=1./A2,RECIP_A3=1./A3)
      PARAMETER(RECIP_A4=1./A4,RECIP_P0=1./P0,RECIP_P_MESO=1./P_MESO)
      PARAMETER(RECIP_P_LOW_STRAT=1./P_LOW_STRAT)
      PARAMETER(RECIP_P_MID_STRAT=1./P_MID_STRAT)
      PARAMETER(MINUS_RLG1=0.19,MINUS_RLG2=-0.029,MINUS_RLG3=-0.082)
      PARAMETER(MINUS_RLG4=0.082,A1_TO_MINUS_RLG1=0.9108)
      PARAMETER(A2_TO_MINUS_RLG2=1.016,A3_TO_MINUS_RLG3=1.048)
      PARAMETER(A4_TO_MINUS_RLG4=0.9448)
      PARAMETER(ONE_THIRD=1./3.,ONE_SIXTH=1./6.)
C
CC *ENDIF
C
CL    END OF COMDECK C_CALCTS

CL    MAXIMUM VECTOR LENGTH IS DETERMINED BY POINTS.
CL
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
C
      IF (LLINTS) THEN
CL
CL    CODE USES SIMPLE LINEAR APPROXIMATION TO TS FOR ALL PRESSURE
CL    VALUES. THIS IS TO MININISE COST.
CL
CL    SECTION 1. CALCULATE TS AT PRESSURE P.
CL---------------------------------------------------------------------

CL    IF CONSTANT_PRESSURE THEN CALCULATE TS FOR ONE POINT ONLY
CL    AND THEN SET ALL OTHER POINTS TO THIS VALUE.
CL    TS = AS0*P + BS0
CL

      IF(CONSTANT_PRESSURE) THEN
          TS(1) = AS0 * P(1) + BS0

        DO 100 I=2,POINTS
          TS(I) = TS(1)
 100    CONTINUE
      ELSE

CL    ELSE CALCULATE TS FOR ALL POINTS.

        DO 200 I=1,POINTS
          TS(I) = AS0 * P(I) + BS0
 200    CONTINUE
      END IF
CL
      ELSE      ! LLINTS
CL
CL    CALCULATE TS CLOSE TO U.M.S.A PROFILE.
CL    EXPONENTIATION HAS BEEN REPLACED BY A TAYLOR SERIES EXPANSION
CL    WHICH IS ACCURATE TO AT WORST ONE DEGREE KELVIN.
CL
CL    SECTION 1. CALCULATE TS AT PRESSURE P.
CL---------------------------------------------------------------------

CL    IF CONSTANT_PRESSURE THEN CALCULATE TS FOR FIRST POINT ONLY
CL    AND THEN SET ALL OTHER POINTS TO THIS VALUE.
CL    OTHERWISE PERFORM CALCULATION AT ALL POINTS.
CL

      IF(CONSTANT_PRESSURE) THEN
C----------------------------------------------------------------------
CL    SECTION 1.1 IF PRESSURE ABOVE P_ISO THEN
CL                TS = TS0*(P0/P)**(R*L0/G)
CL                IS REWRITTEN AS TS=TS0*(P/P0)**(-R*L0/G).
CL                THIS IS REPLACED BY THE TAYLOR EXPANSION STOPPED
CL                AFTER 5 TERMS EXPANDED ABOUT A WHERE A=(PB+PT)/2*PB
CL                THIS SECTION USES PARAMETERS SUFFIXED BY 1.
C----------------------------------------------------------------------

        IF(P(1).GT.P_ISO) THEN
          CONST = (P(1)*RECIP_P0-A1)*RECIP_A1
          TERM1 = 1.+ CONST*(MINUS_RLG1-5.)*ONE_SIXTH
          TERM2 = 1. +TERM1*CONST*(MINUS_RLG1-4.)*.2
          TERM3 = 1. +TERM2*CONST*(MINUS_RLG1-3.)*.25
          TERM4 = 1. +TERM3*CONST*(MINUS_RLG1-2.)*ONE_THIRD
          TERM5 = 1. +TERM4*CONST*(MINUS_RLG1-1.)*.5
          TERM6 = A1_TO_MINUS_RLG1*(1. +TERM5*CONST*MINUS_RLG1)
          TS(1) = TS0* TERM6

C----------------------------------------------------------------------
CL    SECTION 1.2 IF PRESSURE ABOVE P_ISO BUT LESS THAN P_LOW_STRAT THEN
CL                TS = T_ISO WHICH IS A CONSTANT.
C----------------------------------------------------------------------

        ELSE IF(P(1).LE.P_ISO.AND.P(1).GE.P_LOW_STRAT) THEN
          TS(1) = T_ISO

C----------------------------------------------------------------------
CL    SECTION 1.3 IF PRESSURE BELOW P_LOW_STRAT BUT GREATER THAN
CL                P_MID_STRAT
CL                TS = T_ISO*(P_LOW_STRAT/P)**(R*L_LOW_STRAT/G)
CL                TS IS EVALUATED AS DESCRIBED IN SECTION 1.1
CL                PARAMETERS IN THIS SECTION ARE SUFFIXED 2.
C----------------------------------------------------------------------

        ELSE IF(P(1).LT.P_LOW_STRAT.AND.P(1).GT.P_MID_STRAT) THEN
          CONST = (P(1)*RECIP_P_LOW_STRAT-A2)*RECIP_A2
          TERM1 = 1.+ CONST*(MINUS_RLG2-5.)*ONE_SIXTH
          TERM2 = 1. +TERM1*CONST*(MINUS_RLG2-4.)*.2
          TERM3 = 1. +TERM2*CONST*(MINUS_RLG2-3.)*.25
          TERM4 = 1. +TERM3*CONST*(MINUS_RLG2-2.)*ONE_THIRD
          TERM5 = 1. +TERM4*CONST*(MINUS_RLG2-1.)*.5
          TERM6 = A2_TO_MINUS_RLG2*(1. +TERM5*CONST*MINUS_RLG2)
          TS(1) = T_ISO* TERM6

C----------------------------------------------------------------------
CL    SECTION 1.4 IF PRESSURE BELOW P_MID_STRAT BUT GREATER THAN
CL                P_UPPER_STRAT
CL                TS = T_MID_STRAT*(P_MID_STRAT/P)**(R*L_MID_STRAT/G)
CL                TS IS EVALUATED AS DESCRIBED IN SECTION 1.1
CL                PARAMETERS IN THIS SECTION ARE SUFFIXED 3.
C----------------------------------------------------------------------

        ELSE IF(P(1).LT.P_MID_STRAT.AND.P(1).GT.P_UPPER_STRAT) THEN
          CONST = (P(1)*RECIP_P_MID_STRAT-A3)*RECIP_A3
          TERM1 = 1.+ CONST*(MINUS_RLG3-5.)*ONE_SIXTH
          TERM2 = 1. +TERM1*CONST*(MINUS_RLG3-4.)*.2
          TERM3 = 1. +TERM2*CONST*(MINUS_RLG3-3.)*.25
          TERM4 = 1. +TERM3*CONST*(MINUS_RLG3-2.)*ONE_THIRD
          TERM5 = 1. +TERM4*CONST*(MINUS_RLG3-1.)*.5
          TERM6 = A3_TO_MINUS_RLG3*(1. +TERM5*CONST*MINUS_RLG3)
          TS(1) = T_MID_STRAT* TERM6

C----------------------------------------------------------------------
CL    SECTION 1.5 IF PRESSURE BELOW P_UPPER_STRAT BUT GREATER THAN
CL                P_MESO
CL                TS = T_UPPER_STRAT
C----------------------------------------------------------------------

        ELSE IF(P(1).LT.P_UPPER_STRAT.AND.P(1).GT.P_MESO) THEN
          TS(1) = T_UPPER_STRAT

C----------------------------------------------------------------------
CL    SECTION 1.6 IF PRESSURE BELOW P_MESO BUT GREATER THAN P_MIN.
CL                USE STANDARD MESOSPHERE.
CL                TS IS EVALUATED AS DESCRIBED IN SECTION 1.1
CL                PARAMETERS IN THIS SECTION ARE SUFFIXED 4.
C----------------------------------------------------------------------

        ELSE IF(P(1).LT.P_MESO.AND.P(1).GT.P_MIN) THEN
          CONST = (P(1)*RECIP_P_MESO-A4)*RECIP_A4
          TERM1 = 1.+ CONST*(MINUS_RLG4-5.)*ONE_SIXTH
          TERM2 = 1. +TERM1*CONST*(MINUS_RLG4-4.)*.2
          TERM3 = 1. +TERM2*CONST*(MINUS_RLG4-3.)*.25
          TERM4 = 1. +TERM3*CONST*(MINUS_RLG4-2.)*ONE_THIRD
          TERM5 = 1. +TERM4*CONST*(MINUS_RLG4-1.)*.5
          TERM6 = A4_TO_MINUS_RLG4*(1. +TERM5*CONST*MINUS_RLG4)
          TS(1) = T_UPPER_STRAT* TERM6

C----------------------------------------------------------------------
CL    SECTION 1.7 IF PRESSURE BELOW P_MIN.
CL                SET TO T_MIN.
C----------------------------------------------------------------------

        ELSE
          TS(1) = T_MIN

        ENDIF

C----------------------------------------------------------------------
CL    SECTION 1.8  IF CONSTANT_PRESSURE SET TS VALUES FROM POINT
CL                 NUMBER 2 TO POINTS EQUAL TO TS VALUE CALCULATED
CL                 FOR POINT NUMBER 1.
C----------------------------------------------------------------------

        DO 180 I=2,POINTS
          TS(I) = TS(1)
 180    CONTINUE
      ELSE

C NOT CONSTANT PRESSURE SO LOOP OVER ALL POINTS.
C CODE SECTIONS ARE AS ABOVE.

CMIC$ DO ALL VECTOR AUTOSCOPE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 300 I=1,POINTS
C----------------------------------------------------------------------
C     SECTION 2.1 IF PRESSURE ABOVE P_ISO THEN
C                 TS = TS0*(P0/P)**(R*L0/G)
C----------------------------------------------------------------------

          IF(P(I).GT.P_ISO) THEN
            CONST = (P(I)*RECIP_P0-A1)*RECIP_A1
            TERM1 = 1.+ CONST*(MINUS_RLG1-5.)*ONE_SIXTH
            TERM2 = 1. +TERM1*CONST*(MINUS_RLG1-4.)*.2
            TERM3 = 1. +TERM2*CONST*(MINUS_RLG1-3.)*.25
            TERM4 = 1. +TERM3*CONST*(MINUS_RLG1-2.)*ONE_THIRD
            TERM5 = 1. +TERM4*CONST*(MINUS_RLG1-1.)*.5
            TERM6 = A1_TO_MINUS_RLG1*(1. +TERM5*CONST*MINUS_RLG1)
            TS(I) = TS0* TERM6

C----------------------------------------------------------------------
C     SECTION 2.2 IF PRESSURE ABOVE P_ISO BUT LESS THAN P_LOW_STRAT THEN
C                 TS = T_ISO WHICH IS A CONSTANT.
C----------------------------------------------------------------------

          ELSE IF(P(I).LE.P_ISO.AND.P(I).GE.P_LOW_STRAT) THEN
            TS(I) = T_ISO

C----------------------------------------------------------------------
C     SECTION 2.3 IF PRESSURE BELOW P_LOW_STRAT BUT GREATER THAN
C                 P_MID_STRAT
C                 TS = T_ISO*(P_LOW_STRAT/P)**(R*L_LOW_STRAT/G)
C----------------------------------------------------------------------

          ELSE IF(P(I).LT.P_LOW_STRAT.AND.P(I).GT.P_MID_STRAT) THEN
            CONST = (P(I)*RECIP_P_LOW_STRAT-A2)*RECIP_A2
            TERM1 = 1.+ CONST*(MINUS_RLG2-5.)*ONE_SIXTH
            TERM2 = 1. +TERM1*CONST*(MINUS_RLG2-4.)*.2
            TERM3 = 1. +TERM2*CONST*(MINUS_RLG2-3.)*.25
            TERM4 = 1. +TERM3*CONST*(MINUS_RLG2-2.)*ONE_THIRD
            TERM5 = 1. +TERM4*CONST*(MINUS_RLG2-1.)*.5
            TERM6 = A2_TO_MINUS_RLG2*(1. +TERM5*CONST*MINUS_RLG2)
            TS(I) = T_ISO* TERM6

C----------------------------------------------------------------------
C     SECTION 2.4 IF PRESSURE BELOW P_MID_STRAT BUT GREATER THAN
C                 P_UPPER_STRAT
C                 TS = T_MID_STRAT*(P_MID_STRAT/P)**(R*L_MID_STRAT/G)
C----------------------------------------------------------------------

          ELSE IF(P(I).LT.P_MID_STRAT.AND.P(I).GT.P_UPPER_STRAT) THEN
            CONST = (P(I)*RECIP_P_MID_STRAT-A3)*RECIP_A3
            TERM1 = 1.+ CONST*(MINUS_RLG3-5.)*ONE_SIXTH
            TERM2 = 1. +TERM1*CONST*(MINUS_RLG3-4.)*.2
            TERM3 = 1. +TERM2*CONST*(MINUS_RLG3-3.)*.25
            TERM4 = 1. +TERM3*CONST*(MINUS_RLG3-2.)*ONE_THIRD
            TERM5 = 1. +TERM4*CONST*(MINUS_RLG3-1.)*.5
            TERM6 = A3_TO_MINUS_RLG3*(1. +TERM5*CONST*MINUS_RLG3)
            TS(I) = T_MID_STRAT* TERM6

C----------------------------------------------------------------------
C     SECTION 2.5 IF PRESSURE BELOW P_UPPER_STRAT BUT GREATER THAN
C                 P_MESO
C                 TS = T_UPPER_STRAT
C----------------------------------------------------------------------

          ELSE IF(P(I).LT.P_UPPER_STRAT.AND.P(I).GT.P_MESO) THEN
            TS(I) = T_UPPER_STRAT

C----------------------------------------------------------------------
C     SECTION 2.6 IF PRESSURE BELOW P_MESO BUT GREATER THAN P_MIN.
C                 USE STANDARD MESOSPHERE.
C----------------------------------------------------------------------

          ELSE IF(P(I).LT.P_MESO.AND.P(I).GT.P_MIN) THEN
            CONST = (P(I)*RECIP_P_MESO-A4)*RECIP_A4
            TERM1 = 1.+ CONST*(MINUS_RLG4-5.)*ONE_SIXTH
            TERM2 = 1. +TERM1*CONST*(MINUS_RLG4-4.)*.2
            TERM3 = 1. +TERM2*CONST*(MINUS_RLG4-3.)*.25
            TERM4 = 1. +TERM3*CONST*(MINUS_RLG4-2.)*ONE_THIRD
            TERM5 = 1. +TERM4*CONST*(MINUS_RLG4-1.)*.5
            TERM6 = A4_TO_MINUS_RLG4*(1. +TERM5*CONST*MINUS_RLG4)
            TS(I) = T_UPPER_STRAT* TERM6


C----------------------------------------------------------------------
C     SECTION 2.7 IF PRESSURE BELOW P_MIN.
C                 SET TO T_MIN.
C----------------------------------------------------------------------

          ELSE
            TS(I) = T_MIN

          ENDIF

C  END LOOP OVER POINTS
 300    CONTINUE
      END IF

      END IF    ! LLINTS

CL    END OF ROUTINE CALC_TS
      RETURN
      END

