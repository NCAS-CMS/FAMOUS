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
CLL  SUBROUTINE DDRAUGHT-----------------------------------------------
CLL
CLL  PURPOSE : DOWNDRAUGHT ROUTINE
CLL
CLL            CONVECTIVE DOWNDRAUGHT BASED ON PARCEL THEORY
CLL
CLL            CARRY OUT DRY DESCENT
CLL
CLL            CALCULATE SUBSATURATION
CLL
CLL            CALCULATE EFFECT ON THE ENVIRONMENT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DDRAUGHT (NPNTS,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1,QE_K,
     &                     QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &                     DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,
     &                     DELPKM1,EXK,EXKM1,DELTD,DELQD,AMDETK,EKM14,
     &                     EKM34,RAIN,SNOW,BDD_START,BDDWT_K,
     &                     BDDWT_KM1,BDD_ON,B_DD_END,CCA)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL DET_LYR  ! THICKNESS LEVEL USED IN CALCULATION OF MIXING
                    ! DETRAINMENT FOR DOWNDRAUGHT  (PA)
      PARAMETER (DET_LYR = 10000.0)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                     ! LOOP COUNTER
C
      INTEGER NPNTS                 ! IN NUMBER OF POINTS
C
      INTEGER K                     ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                   ! IN CONVECTIVE CLOUD TOP
C
      REAL THE_KM1(NPNTS)           ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_KM1(NPNTS)            ! IN MIXING RATIO OF ENVIRONMENT IN
                                    !    LAYER K-1 (KG/KG)
C
      REAL P_KM1(NPNTS)             ! IN PRESSURE OF LAYER K-1 (PA)
C
      REAL DELPK(NPNTS)             ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
C
      REAL DELPKM1(NPNTS)           ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
C
      REAL EXK(NPNTS)               ! IN EXNER RATIO IN LAYER K
C
      REAL EXKM1(NPNTS)             ! IN EXNER RATIO IN LAYER K-1
C
      REAL AMDETK(NPNTS)            ! IN MIXING DETRAINMENT RATE
C
      REAL EKM14(NPNTS)             ! IN EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(NPNTS)             ! IN EXNER RATIO AT LAYER K-3/4
C
      REAL DELTD(NPNTS)             ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
C
      REAL DELQD(NPNTS)             ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
C
      LOGICAL BDDWT_K(NPNTS)        ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K
C
      LOGICAL BDDWT_KM1(NPNTS)      ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K-1
C
      REAL CCA(NPNTS)               ! IN CONVECTIVE CLOUD AMOUNT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(NPNTS)            ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     DOWNDRAUGHT IN LAYER K (K)
                                    ! OUT POTENTIAL TEMPERATURE RESET
                                    !     FOR NEXT LAYER (K)
C
      REAL QDD_K(NPNTS)             ! INOUT
                                    ! IN  DOWNDRAUGHT MIXING RATIO OF
                                    !     LAYER K (KG/KG)
                                    ! OUT MIXING RATIO RESET FOR NEXT
                                    !     LAYER (KG/KG)
C
      REAL THE_K(NPNTS)             ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     ENVIRONMENT IN LAYER K (K)
                                    ! OUT ENVIRONMENT POTENTIAL
                                    !     TEMPERATURE RESET FOR NEXT
                                    !     LAYER (K)
C
      REAL QE_K(NPNTS)              ! INOUT
                                    ! IN  MIXING RATIO OF ENVIRONMENT
                                    !     LAYER K (KG/KG)
                                    ! OUT ENVIRONMENT MIXING RATIO
                                    !     RESET FOR NEXT LAYER (KG/KG)
C
      REAL FLX_DD_K(NPNTS)          ! INOUT
                                    ! IN  DOWNDRAUGHT MASS FLUX OF
                                    !     LAYER K (PA/S)
                                    ! OUT DOWNDRAUGHT MASS FLUX RESET
                                    !      FOR NEXT LAYER (PA/S)
C
      REAL RAIN(NPNTS)              ! INOUT
                                    ! IN  AMOUNT OF RAIN (KG/M**2/S)
                                    ! OUT UPDATED RAINFALL (KG/M**2/S)
C
      REAL SNOW(NPNTS)              ! INOUT
                                    ! IN  AMOUNT OF SNOW(KG/M**2/S)
                                    ! OUT UPDATED SNOWFALL (KG/M**2/S)
C
      REAL DTHBYDT_K(NPNTS)         ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE OF LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K (K/S)
C
      REAL DTHBYDT_KM1(NPNTS)       ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE IN LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K-1 (K/S)
C
      REAL DQBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K
                                    !     (KG/KG/S)
C
      REAL DQBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K-1
                                    !     (KG/KG/S)
C
      LOGICAL BDD_ON(NPNTS)         ! INOUT
                                    ! IN  MASK FOR THOSE POINTS WHERE DD
                                    !     HAS CONTINUED FROM LAYER K+1
                                    ! OUT MASK FOR THOSE POINTS WHERE DD
                                    !     CONTINUES TO LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(NPNTS)      ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT MAY START IN
                                    !     LAYER K-1
C
      LOGICAL B_DD_END(NPNTS)       ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT IS ENDING IN
                                    !     LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL THDD_KM1(NPNTS)          ! POTENTIAL TEMPERATURE OF
                                    ! DOWNDRAUGHT IN LAYER K-1 (K)
C
      REAL QDD_KM1(NPNTS)           ! DOWNDRAUGHT MIXING RATIO OF
                                    ! LAYER K-1 (KG/KG)
C
      REAL QSATDD(NPNTS)            ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
C
      REAL TDD_KM1(NPNTS)           ! TEMPERATURE OF DOWNDRAUGHT
                                    ! IN LAYER K-1 (K)
C
      REAL THDDS(NPNTS)             ! POTENTIAL TEMPERATURE OF
                                    ! SATURATED DOWNDRAUGHT (K)
C
      REAL QDDS(NPNTS)              ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
C
      REAL FLX_DD_KM1(NPNTS)        ! DOWNDRAUGHT MASS FLUX IN
                                    ! LAYER K-1 (PA/S)
C
      REAL RAIN_TMP(NPNTS)          ! LIQUID PRECIPITATION STORE
C
      REAL SNOW_TMP(NPNTS)          ! SNOW STORE
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL SATCAL, CRS_FRZL, QSAT, DEVAP, TERMDD,
     *         DD_ENV, EVP
C
C-----------------------------------------------------------------------
C CALCULATE MASK FOR THOSE POINTS IN DOWNDRAUGHT WHERE PRECIPITATION
C IS LIQUID
C
C STORE PRECIPITATION IN LAYER K IN TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (K .EQ. KCT+1 .OR. BDD_START(I)) THEN
          BDDWT_K(I) = THDD_K(I) .GT. TM/EXK(I)
        ELSE
          BDDWT_K(I) = BDDWT_KM1(I)
        END IF
          RAIN_TMP(I) = RAIN(I)
          SNOW_TMP(I) = SNOW(I)
C
C-----------------------------------------------------------------------
C DRY DESCENT FROM LAYER K TO K-1
C
C ENTRAINMENT CALCULATION
C-----------------------------------------------------------------------
C
          THDD_KM1(I) = (THDD_K(I)+(EKM14(I)*THE_K(I)) +
     *                  (1.0+EKM14(I))*EKM34(I)*THE_KM1(I)) /
     *                  ((1.0+EKM14(I))*(1.0+EKM34(I)))
          QDD_KM1(I) = (QDD_K(I)+(EKM14(I)*QE_K(I)) +
     *                 (1.0+EKM14(I))*EKM34(I)*QE_KM1(I))/
     *                 ((1.0+EKM14(I))*(1.0+EKM34(I)))
C
C-----------------------------------------------------------------------
C UPDATE MASS FLUX  AND CALCULATE TEMPERATURE OF LAYER K-1
C-----------------------------------------------------------------------
C
          FLX_DD_KM1(I) = FLX_DD_K(I)*(1.0+EKM34(I))*(1.0+EKM14(I))*
     *                (1.0-AMDETK(I))
C
          TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE SUBSATURATION
C CALCULATE TEMPERATURE IF BROUGHT TO SATURATION
C-----------------------------------------------------------------------
C
       CALL SATCAL(NPNTS,TDD_KM1,THDD_KM1,P_KM1,QDDS,THDDS,
     &             K,EXKM1,QDD_KM1,THE_KM1)
C
      DO I=1,NPNTS
        BDDWT_KM1(I) = THDDS(I) .GT. TM/EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE CHANGE OF PHASE DUE TO DOWNDRAUGHT SATURATION TEMPERATURE
C-----------------------------------------------------------------------
C
       CALL CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,FLX_DD_KM1,
     &                BDDWT_KM1)
C
      DO I=1,NPNTS
        TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C RECALCULATE SUBSATURATION TEMPERATURE
C-----------------------------------------------------------------------
C
       CALL SATCAL(NPNTS,TDD_KM1,THDD_KM1,P_KM1,QDDS,THDDS,
     &             K,EXKM1,QDD_KM1,THE_KM1)
C
C-----------------------------------------------------------------------
C CALCULATE MOISTURE SUBSATURATION
C-----------------------------------------------------------------------
C
       CALL QSAT(QSATDD,TDD_KM1,P_KM1,NPNTS)
C
C-----------------------------------------------------------------------
C EVAPORATION CALCULATION AND ADJUSTMENT OF DOWNDRAUGHT TEMPERATURE
C AND MOISTURE
C-----------------------------------------------------------------------
C
       CALL DEVAP (NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,
     &             FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,
     &             DELPKM1,BDDWT_KM1,CCA,P_KM1)
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT CAN
C CONTINUE TO K-1
C-----------------------------------------------------------------------
C
       CALL TERMDD (NPNTS,BDD_START,THDD_KM1,QDD_KM1,THE_KM1,
     &              QE_KM1,K,B_DD_END,BDD_ON)
C
C-----------------------------------------------------------------------
C CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
C-----------------------------------------------------------------------
C
       CALL DD_ENV (NPNTS,THDD_K,THDD_KM1,QDD_K,QDD_KM1,THE_K,THE_KM1,
     &              QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &              DQBYDT_KM1,FLX_DD_K,FLX_DD_KM1,DELPK,DELPKM1,
     &              DELTD,DELQD,AMDETK,EKM14,B_DD_END,
     &              BDD_START,BDD_ON)
C
C-----------------------------------------------------------------------
C RESET DOWNDRAUGHT BIT VECTORS
C
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       BDD_START(I) = .FALSE.
       IF (.NOT. BDD_ON(I)) THEN
         RAIN(I) = RAIN_TMP(I)
         SNOW(I) = SNOW_TMP(I)
       END IF
       IF (B_DD_END(I)) BDD_ON(I) = .FALSE.
      END DO
C
C-----------------------------------------------------------------------
C SWITCH POTENTIAL TEMPERATURE, MIXING RATIO AND MASS FLUX FOR
C CALCULATION AT NEXT MODEL LAYER
C-----------------------------------------------------------------------
C
      IF (K.GT.2) THEN
        DO I=1,NPNTS
         IF (BDD_ON(I)) THEN
          THDD_K(I) = THDD_KM1(I)
          QDD_K(I) = QDD_KM1(I)
          FLX_DD_K(I) = FLX_DD_KM1(I)
         END IF
        END DO
      END IF

      RETURN
      END
C
