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
CLL  SUBROUTINE DEVAP--------------------------------------------------
CLL
CLL  PURPOSE : EVAPORATION ROUTINE
CLL
CLL            CARRIES OUT EVAPORATION AND UPDATES PRECIPITATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.0    5/05/95   New deck at version 4.0 to include pressure
CLL                   dependency into calculation of evaporation of
CLL                   convective precipitation, and to introduce
CLL                   traps for negative precipitation.
CLL                   Pete Inness.
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DEVAP(NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,
     *                 FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,
     *                 DELPKM1,BDDWT_KM1,CCA,PKM1)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
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

      REAL DDCLDFRA               ! FRACTIONAL CLOUD AREA OF DD
      PARAMETER (DDCLDFRA=0.5)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I               ! LOOP COUNTER
C
      INTEGER NPNTS           ! IN VECTOR LENGTH
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDDS(NPNTS)       ! IN SATURATED POTENTIAL
                              !    TEMPERATURE OF DOWNDRAUGHT
                              !    (K)
C
      REAL QDDS(NPNTS)        ! IN MIXING RATIO OF SATURATED
                              !    DOWNDRAUGHT (KG/KG)
C
      REAL FLX_DD_KM1(NPNTS)  ! IN DOWNDRAUGHT MASS FLUX IN
                              !    LAYER K-1 (PA/S)
C
      REAL THDD_K(NPNTS)      ! IN POTENTIAL TEMPERATURE OF
                              !    DOWNDRAUGHT IN LAYER K (K)
C
      REAL EXK(NPNTS)         ! IN EXNER RATIO OF LAYER K
C
      REAL EXKM1(NPNTS)       ! IN EXNER RATIO OF LAYER K-1
C
      REAL QSATDD(NPNTS)      ! IN SATURATED DOWNDRAUGHT
                              !    MIXING RATIO (KG/KG)
C
      REAL DELPKM1(NPNTS)     ! IN CHANGE IN PRESSURE ACROSS
                              !    LAYER K-1 (PA)
C
      LOGICAL BDDWT_KM1(NPNTS)! IN MASK WHERE PRECIPITATION IN
                              !    DOWNDRAUGHT IS LIQUID
C
      REAL CCA(NPNTS)         ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL PKM1(NPNTS)        ! IN PRESSURE OF LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_KM1(NPNTS)    ! INOUT
                              ! IN  POTENTIAL TEMPERATURE OF
                              !     DOWNDRAUGHT IN LAYER K-1 (K)
                              ! OUT UPDATED POTENTIAL TEMPERATURE
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (K)
C
      REAL QDD_KM1(NPNTS)     ! INOUT
                              ! IN  MODEL MIXING RATIO OF
                              !     DOWNDRAUGHT IN LAYER K-1
                              !     (KG/KG)
                              ! OUT UPDATED MODEL MIXING RATIO
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (KG/KG)
C
      REAL RAIN(NPNTS)        ! INOUT
                              ! IN  AMOUNT OF RAIN (KG/M**2/S)
                              ! OUT UPDATED RAINFALL (KG/M**2/S)
C
      REAL SNOW(NPNTS)        ! INOUT
                              ! IN  AMOUNT OF SNOW (KG/M**2/S)
                              ! OUT UPDATED SNOWFALL (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
C
      REAL TEVP(NPNTS)        ! TEMPERATURE USED IN EVAPORATION
                              ! CALCULATION (K)
C
      LOGICAL BEVAP(NPNTS)    ! MASK FOR THOSE POINTS AT WHICH
                              ! EVAPORATION CALCULATION IS TO
                              ! BE CARRIED OUT
C
      LOGICAL BSAT(NPNTS)     ! MASK FOR THOSE POINTS WHICH
                              ! ARE SUBSATURATED
C
      REAL EVAP_RAIN(NPNTS)   ! AMOUNT OF EVAPORATION OF RAIN
C
      REAL SUB_SNOW(NPNTS)    ! AMOUNT OF SNOW SUBLIMATION
C
      REAL DELQ(NPNTS)        ! DIFFERENCE IN MIXING RATIOS
                              ! (KG/KG)
C
      REAL DELTH(NPNTS)       ! INCREMENT TO DOWNDRAUGHT POTENTIAL
                              ! TEMPERATURE IN LAYER K-1 DUE TO
                              ! EVAPORATION
C
      REAL DELQE(NPNTS)       ! INCREMENT TO DOWNDRAUGHT MIXING RATIO
                              ! IN LAYER K-1 DUE TO EVAPORATION
C
      REAL DELTHS(NPNTS)      ! SATURATED POTENTIAL TEMPERATURE MINUS
                              ! POTENTIAL TEMPERATURE OF DOWNDRAUGHT
C
      REAL FACTOR(NPNTS)      ! DELTHS / DELTH
C
      REAL PINCR(NPNTS)       ! INCREASE IN PRECIPITATION IF PARCEL
                              ! SUPERSATURATES
C
      REAL RHO(NPNTS)         ! DENSITY OF AIR IN PARCEL
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL EVP
C
C-----------------------------------------------------------------------
C CHECK IF EVAPORATION POSSIBLE
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       DELQ(I) = QSATDD(I)-QDD_KM1(I)
C
       BEVAP(I) =((RAIN(I).GT.0.0) .OR. (SNOW(I).GT.0.0))
     *             .AND. (DELQ(I).GT.0.0)
       BSAT(I) = DELQ(I) .LT. 0.0
C
C-----------------------------------------------------------------------
C CALCULATE TEMPERATURE USED IN CALCULATION OF EVAPORATION CONSTANTS
C BASED ON TEMPERATURE OF PARCEL AFTER UNSATURATED DESCENT
C-----------------------------------------------------------------------
C
        IF (BEVAP(I)) THEN
          TEVP(I) = ((THDD_K(I)*EXK(I))+(THDD_KM1(I)*EXKM1(I)))*0.5
          RHO(I) = PKM1(I) / (R*TEVP(I))
        END IF
      END DO
C
C-----------------------------------------------------------------------
C EVAPORATION CALCULATION - CALCULATE RATES FOR RAIN AND SNOW
C-----------------------------------------------------------------------
C
      CALL EVP(NPNTS,RAIN,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP_RAIN,
     &         BEVAP,1,DDCLDFRA,PKM1)
C
      CALL EVP(NPNTS,SNOW,TEVP,CCA,RHO,DELQ,DELPKM1,SUB_SNOW,
     &         BEVAP,2,DDCLDFRA,PKM1)
C
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
C
C-----------------------------------------------------------------------
C ADJUST EVAPORATION AND SUBLIMATION RATES BACK TO GRID BOX MEANS
C-----------------------------------------------------------------------
C
       EVAP_RAIN(I) = EVAP_RAIN(I) * CCA(I) * DDCLDFRA
       SUB_SNOW(I) = SUB_SNOW(I) * CCA(I) * DDCLDFRA
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL SUPERSATURATED
C-----------------------------------------------------------------------
C
        DELTH(I) = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/
     &           (CP*EXKM1(I)*FLX_DD_KM1(I))
        DELQE(I) = (EVAP_RAIN(I)+SUB_SNOW(I))*G/FLX_DD_KM1(I)
C
        DELTHS(I) = THDDS(I)-THDD_KM1(I)
        IF (DELTH(I).LT.DELTHS(I)) THEN
C
C-----------------------------------------------------------------------
C ADJUST EVAP AND SUBLIMATION RATES TO GIVE SATURATION
C-----------------------------------------------------------------------
C
          FACTOR(I) = DELTHS(I)/DELTH(I)
          DELTH(I) = DELTHS(I)
          DELQE(I) = DELQE(I)*FACTOR(I)
          EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR(I)
          SUB_SNOW(I) = SUB_SNOW(I)*FACTOR(I)
        END IF
C
C-----------------------------------------------------------------------
C UPDATE T,Q AND PRECIPITATION
C-----------------------------------------------------------------------
C
        RAIN(I) = RAIN(I)-EVAP_RAIN(I)
        IF (RAIN(I).LT.0.0) RAIN(I)=0.0
        SNOW(I) = SNOW(I)-SUB_SNOW(I)
        IF (SNOW(I).LT.0.0) SNOW(I)=0.0
        THDD_KM1(I) = THDD_KM1(I)+DELTH(I)
        QDD_KM1(I) = QDD_KM1(I)+DELQE(I)
C
C-----------------------------------------------------------------------
C PARCEL IS SUPERSATURATED BEFORE EVAPORATION OCCURS
C BRING PARCEL TO SATURATION AND PRECIPITATE WATER
C-----------------------------------------------------------------------
C
      ELSE IF (BSAT(I)) THEN
         PINCR(I) = (QDD_KM1(I)-QDDS(I))*FLX_DD_KM1(I)/G
         QDD_KM1(I) = QDDS(I)
         THDD_KM1(I) = THDDS(I)
         IF (BDDWT_KM1(I)) THEN
           RAIN(I) = RAIN(I)+PINCR(I)
         ELSE
           SNOW(I) = SNOW(I)+PINCR(I)
         END IF
      END IF
      END DO
C
      RETURN
      END
C
