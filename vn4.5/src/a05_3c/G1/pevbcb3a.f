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
CLL  SUBROUTINE PEVP_BCB-----------------------------------------------
CLL
CLL  PURPOSE : EVAPORATE RAIN BELOW CLOUD BASE IF NO DOWNDRAUGHT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.0   5/05/95   New deck at version 4.0 to include a pressure
CLL                   dependency in the calculation of evaporation
CLL                   of convective precipitation.
CLL                   Pete Inness.
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4 DATED 5/2/92
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
      SUBROUTINE PEVP_BCB (NPNTS,K,ICCB,TH,PK,Q,DELP,RAIN,SNOW,
     *                     DTHBYDT,DQBYDT,EXK,TIMESTEP,CCA)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C  CONSTANTS
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
      REAL CLDAREA                ! FRACTIONAL CLOUD AREA WHEN NO
                                  ! DD, USED IN EVAPORATION CALC
      PARAMETER (CLDAREA=1.0)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! IN LOOP COUNTER
C
      INTEGER NPNTS              ! VECTOR LENGTH
C
      INTEGER K                  ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER ICCB(NPNTS)        ! IN CONVECTIVE CLOUD BASE LAYER
C
      REAL PK(NPNTS)             ! IN PRESSURE (PA)
C
      REAL Q(NPNTS)              ! IN MIXING RATIO (KG/KG)
C
      REAL TH(NPNTS)             ! IN POTENTIAL TEMPERATURE (K)
C
      REAL DELP(NPNTS)           ! IN CHANGE IN PRESSURE ACROSS
                                 !    LAYER K-1 (PA)
C
      REAL EXK(NPNTS)            ! IN EXNER RATIO OF LAYER K
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL CCA(NPNTS)            ! IN CONVECTIVE CLOUD AMOUNT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHBYDT(NPNTS)        ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE (K/S)
C
      REAL DQBYDT(NPNTS)         ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING RATIO
                                 !     (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MIXING RATIO
                                 !     AFTER EVAPORATION BELOW CLOUD
                                 !     BASE (KG/KG/S)
C
      REAL RAIN(NPNTS)           ! INOUT
                                 ! IN  AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
C
      REAL SNOW(NPNTS)           ! INOUT
                                 ! IN  AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL T(NPNTS)              ! MODEL TEMPERATURE (K)
C
      REAL EVAP_RAIN(NPNTS)      ! AMOUNT OF EVAPORATION OF RAIN
C
      REAL SUB_SNOW(NPNTS)       ! AMOUNT OF SNOW SUBLIMATION
C
      REAL QSATE(NPNTS)          ! SATURATED MIXING RATIO IN
                                 ! ENVIRONMENT (KG/KG)
C
      REAL DELQ(NPNTS)           ! CHANGE IN MIXING RATIO ACROSS LAYER K
                                 ! (KG/KG)
C
      REAL THS(NPNTS)            ! SATURATED PARCEL POTENTIAL
                                 ! TEMPERATURE (K)
C
      REAL QS(NPNTS)             ! SATURATED PARCEL MIXING RATIO
C
      LOGICAL BEVAP(NPNTS)       ! MASK FOR THOSE POINTS WHERE
                                 ! EVAPORATION OCCURS
C
      REAL DTHBYDT_EVP(NPNTS)    ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO EVAPORATION (K)
C
      REAL DQBYDT_EVP(NPNTS)     ! INCREMENT TO MIXING RATIO DUE TO
                                 ! EVAPORATION (KG/KG)
C
      REAL DTHBYDT_SAT(NPNTS)    ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO SATURATION (K)
C
      REAL FACTOR(NPNTS)         ! DTHBYDT_SAT / DTHBYDT_EVP
C
      REAL RHO(NPNTS)            ! DENSITY OF AIR IN PARCEL
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL QSAT, EVP, SATCAL
C
C-----------------------------------------------------------------------
C EVAPORATE RAIN IN LAYER K IF LAYER K IS BELOW CLOUD BASE
C CALCULATE MOISTURE SUB-SATURATION
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        T(I) = TH(I)*EXK(I)
        BEVAP(I) = .FALSE.
      END DO
C
      CALL QSAT(QSATE,T,PK,NPNTS)
C
      DO I=1,NPNTS
       IF (K .LT. ICCB(I)) THEN
         DELQ(I) = QSATE(I)-Q(I)
C
C-----------------------------------------------------------------------
C CHECK IF EVAPORATION POSSIBLE
C-----------------------------------------------------------------------
C
         IF ((RAIN(I).GT.0.0 .OR. SNOW(I).GT.0.0) .AND.
     &        DELQ(I) .GT. 0.0) THEN
C
            BEVAP(I) = .TRUE.
            RHO(I) = PK(I) / (R*T(I))
         END IF
       END IF
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE EVAPORATION
C-----------------------------------------------------------------------
C
        CALL EVP (NPNTS,RAIN,T,CCA,RHO,DELQ,DELP,EVAP_RAIN,
     &            BEVAP,1,CLDAREA,PK)
C
        CALL EVP (NPNTS,SNOW,T,CCA,RHO,DELQ,DELP,SUB_SNOW,
     &            BEVAP,2,CLDAREA,PK)
C
C-----------------------------------------------------------------------
C CALCULATE TEMPERATURE AND MIXING RATIO IF LAYER BROUGHT TO
C SATURATION BY EVAPORATION AND SUBLIMATION
C-----------------------------------------------------------------------
C
      CALL SATCAL(NPNTS,T,TH,PK,QS,THS,K,EXK,Q,TH)
C
C
      DO I=1,NPNTS
        IF (BEVAP(I)) THEN
          DTHBYDT_EVP(I) = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/
     &                   (CP*EXK(I)*DELP(I))
          DQBYDT_EVP(I) = (EVAP_RAIN(I)+SUB_SNOW(I))*G/DELP(I)
C
          DTHBYDT_SAT(I) = (THS(I)-TH(I))/TIMESTEP
C
          IF (DTHBYDT_EVP(I).LT.DTHBYDT_SAT(I)) THEN
C
C---------------------------------------------------------------------
C  ADJUST EVAPORATION AND SUBLIMATION RATES TO GIVE SATURATION
C---------------------------------------------------------------------
C
            FACTOR(I) = DTHBYDT_SAT(I)/DTHBYDT_EVP(I)
            DTHBYDT_EVP(I) = DTHBYDT_SAT(I)
            DQBYDT_EVP(I) = DQBYDT_EVP(I)*FACTOR(I)
            EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR(I)
            SUB_SNOW(I) = SUB_SNOW(I)*FACTOR(I)
          END IF
C
C---------------------------------------------------------------------
C  UPDATE INCREMENTS AND RAINFALL AND ADJUST BACK TO GRIDBOX MEANS
C---------------------------------------------------------------------
C
          DTHBYDT(I) = DTHBYDT(I)+DTHBYDT_EVP(I)*CCA(I)*CLDAREA
          DQBYDT(I) = DQBYDT(I)+DQBYDT_EVP(I)*CCA(I)*CLDAREA
          RAIN(I) = RAIN(I)-EVAP_RAIN(I)*CCA(I)*CLDAREA
          SNOW(I) = SNOW(I)-SUB_SNOW(I)*CCA(I)*CLDAREA
        END IF
      END DO
C
      RETURN
      END
C
