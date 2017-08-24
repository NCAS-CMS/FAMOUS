C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE DOWND--------------------------------------------------
CLL
CLL  PURPOSE : CALL DOWNDRAUGHT CALCULATION
CLL
CLL            CHANGE OF PHASE CALCULATION WHERE NO DOWNDRAUGHT OCCURS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL:
CLL VERSION  DATE
CLL   4.2   1/11/96   New deck version based on DOWND2A with HADCM2
CLL                   specific modifications: R Jones
!LL  4.5   23/02/98  Call comdecks. D. Robinson
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DOWND (NPNTS,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1,QE_K,
     &                  QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &                  DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,DELPKM1,EXK,
     &                  EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K,
     &                  RAIN,SNOW,ICCB,BWATER_K,BDD_START,
     &                  BDDWT_K,BDDWT_KM1,BDD_ON,RAIN_ENV,SNOW_ENV,
     &                  RAIN_DD,SNOW_DD,FLX_UD_K,TIMESTEP,CCA,NDDON_A)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! LOOP COUNTER
C
      INTEGER K                  ! IN PRESENT MODEL LAYER
C
      INTEGER NPNTS              ! IN NUMBER OF POINTS
C
      INTEGER NDDON,NDDON_A      ! NUMBER OF POINTS AT WHICH
                                 ! DOWNDRAUGHT DOES OCCUR
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
C
      REAL THDD_K(NPNTS)         ! IN MODEL POTENTIAL TEMPERATURE
                                 !    OF DOWNDRAUGHT IN LAYER K (K)
C
      REAL QDD_K(NPNTS)          ! IN MIXING RATIO OF DOWNDRAUGHT IN
                                 !    LAYER K (KG/KG)
C
      REAL THE_K(NPNTS)          ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(NPNTS)        ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_K(NPNTS)           ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K (KG/KG)
C
      REAL QE_KM1(NPNTS)         ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K-1 (KG/KG)
C
      REAL FLX_DD_K(NPNTS)       ! IN DOWNDRAUGHT MASS FLUX OF LAYER K
                                 !    (PA/S)
C
      REAL P_KM1(NPNTS)          ! IN PRESSURE OF LAYER K-1 (PA)
C
      REAL DELPK(NPNTS)          ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K (PA)
C
      REAL DELPKM1(NPNTS)        ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K-1 (PA)
C
      REAL EXK(NPNTS)            ! IN EXNER RATIO FOR LAYER K
C
      REAL EXKM1(NPNTS)          ! IN EXNER RATIO FOR LAYER K-1
C
      REAL PRECIP_K(NPNTS)       ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
C
      REAL AMDETK(NPNTS)         ! IN MIXING DETRAINMENT AT LEVEL K
                                 !    MULTIPLIED BY APPROPRIATE LAYER
                                 !    THICKNESS
C
      REAL EKM14(NPNTS)          ! IN EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(NPNTS)          ! IN EXNER RATIO AT LAYER K-3/4
C
      REAL DELTD(NPNTS)          ! IN COOLING NECESSARY TO
                                 !    ACHIEVE SATURATION (K)
C
      REAL DELQD(NPNTS)          ! IN MOISTENING NECESSARY TO
                                 !    ACHIEVE SATURATION (KG/KG)
C
      REAL ICCB(NPNTS)           ! IN CLOUD BASE LEVEL
C
      LOGICAL BWATER_K(NPNTS)    ! IN MASK FOR THOSE POINTS AT WHICH
                                 !    CONDENSATE IS WATER IN LAYER K
C
      LOGICAL BDDWT_K(NPNTS)     ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K
C
      LOGICAL BDDWT_KM1(NPNTS)   ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K-1
C
      REAL RAIN_ENV(NPNTS)       ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL SNOW_ENV(NPNTS)       ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL RAIN_DD(NPNTS)        ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL SNOW_DD(NPNTS)        ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL FLX_UD_K(NPNTS)       ! IN UPDRAUGHT MASSFLUX AT LAYER K
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL CCA(NPNTS)            ! IN CONVECTIVE CLOUD AMOUNT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(NPNTS)   ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER K
                                 ! OUT MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER
                                 !     K-1
C
      REAL DTHBYDT_K(NPNTS)      ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF LAYER K
                                 !     (K/S)
C
      REAL DTHBYDT_KM1(NPNTS)    ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K-1 (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (K/S)
C
      REAL DQBYDT_K(NPNTS)       ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1(NPNTS)     ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (KG/KG/S)
C
      REAL RAIN (NPNTS)          ! INOUT
                                 ! IN  INITIALISED RAINFALL (KG/M**2/S)
                                 ! OUT SURFACE RAINFALL (KG/M**2/S)
C
      REAL SNOW(NPNTS)           ! INOUT
                                 ! IN  INITIALISED SNOWFALL (KG/M**2/S)
                                 ! OUT SURFACE SNOWFALL (KG/M**2/S)
C
      LOGICAL BDD_ON(NPNTS)      ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE DD
                                 !     HAS CONTINUED FROM PREVIOUS LAYER
                                 ! OUT MASK FOR THOSE POINTS WHERE DD
                                 !     CONTINUES TO LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL WORK(NDDON_A,24)      !  WORK SPACE
C
      LOGICAL BWORK(NDDON_A,5)   !  WORK SPACE FOR 'BIT' MASKS
C
      INTEGER INDEX1(NDDON_A)    !  INDEX FOR COMPRESS AND
C
      LOGICAL B_DD_END(NPNTS)    !  MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! HAS ENDED
C
C
      REAL FACTOR                !  PROPORTION OF RAINFALL GOING
                                 !  THROUGH DOWNDRAUGHT
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL CHG_PHSE, PEVP_BCB, DDRAUGHT
C
C-----------------------------------------------------------------------
C START OF MAIN LOOP
C   UPDATE PRECIPITATION AND CALCULATE MASK FOR WHERE PRECIPITATION
C   IS LIQUID
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        B_DD_END(I) = .FALSE.
      END DO
C
      IF (K.EQ.KCT+1) THEN
        DO I=1,NPNTS
         RAIN_DD(I) = 0.0
         RAIN_ENV(I) = 0.0
         SNOW_DD(I) = 0.0
         SNOW_ENV(I) = 0.0
        END DO
      END IF
C
      DO I=1,NPNTS
C
         IF (BDD_ON(I)) THEN
           FACTOR = 1.0
         ELSE
           FACTOR = 0.0
         END IF
C
         IF (BWATER_K(I)) THEN
           RAIN_DD(I) = RAIN_DD(I)+(PRECIP_K(I)+RAIN_ENV(I))*FACTOR
           RAIN_ENV(I) = (RAIN_ENV(I)+PRECIP_K(I))*(1.0-FACTOR)
           SNOW_DD(I) = SNOW_DD(I) + SNOW_ENV(I)*FACTOR
           SNOW_ENV(I) = SNOW_ENV(I)*(1.0-FACTOR)
         ELSE
           SNOW_DD(I) = SNOW_DD(I)+(PRECIP_K(I)+SNOW_ENV(I))*FACTOR
           SNOW_ENV(I) = (SNOW_ENV(I)+PRECIP_K(I))*(1.0-FACTOR)
           RAIN_DD(I) = RAIN_DD(I) + RAIN_ENV(I)*FACTOR
           RAIN_ENV(I) = RAIN_ENV(I)*(1.0-FACTOR)
         END IF
         PRECIP_K(I) = 0.0
C
      END DO
C
C
C-----------------------------------------------------------------------
C COMPRESS OUT ON BASIS OF BIT VECTOR BDDON - THOSE POINTS WITH A
C DOWNDRAUGHT
C-----------------------------------------------------------------------
C
      NDDON=0
C
      DO I=1,NPNTS
        IF (BDD_ON(I)) THEN
           NDDON = NDDON+1
           INDEX1(NDDON) = I
        END IF
      END DO
C
      IF (NDDON .NE. 0) THEN
         DO I=1,NDDON
          WORK(I,1) = THDD_K(INDEX1(I))
          WORK(I,2) = QDD_K(INDEX1(I))
          WORK(I,3) = THE_K(INDEX1(I))
          WORK(I,4) = THE_KM1(INDEX1(I))
          WORK(I,5) = QE_K(INDEX1(I))
          WORK(I,6) = QE_KM1(INDEX1(I))
          WORK(I,7) = DTHBYDT_K(INDEX1(I))
          WORK(I,8) = DTHBYDT_KM1(INDEX1(I))
          WORK(I,9) = DQBYDT_K(INDEX1(I))
          WORK(I,10) = DQBYDT_KM1(INDEX1(I))
          WORK(I,11) = FLX_DD_K(INDEX1(I))
          WORK(I,12) = P_KM1(INDEX1(I))
          WORK(I,13) = DELPK(INDEX1(I))
          WORK(I,14) = DELPKM1(INDEX1(I))
          WORK(I,15) = EXK(INDEX1(I))
          WORK(I,16) = EXKM1(INDEX1(I))
          WORK(I,17) = DELTD(INDEX1(I))
          WORK(I,18) = DELQD(INDEX1(I))
          WORK(I,19) = AMDETK(INDEX1(I))
          WORK(I,20) = EKM14(INDEX1(I))
          WORK(I,21) = EKM34(INDEX1(I))
          WORK(I,22) = RAIN_DD(INDEX1(I))
          WORK(I,23) = SNOW_DD(INDEX1(I))
          WORK(I,24) = CCA(INDEX1(I))
          BWORK(I,1) = BDD_START(INDEX1(I))
          BWORK(I,2) = BDDWT_K(INDEX1(I))
          BWORK(I,3) = BDDWT_KM1(INDEX1(I))
          BWORK(I,4) = BDD_ON(INDEX1(I))
          BWORK(I,5) = B_DD_END(INDEX1(I))
      END DO
C
C-----------------------------------------------------------------------
C START DOWNDRAUGHT CALCULATION
C-----------------------------------------------------------------------
C
C
         CALL DDRAUGHT (NDDON,K,KCT,WORK(1,1),WORK(1,2),WORK(1,3),
     &                  WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7),
     &                  WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11),
     &                  WORK(1,12),WORK(1,13),WORK(1,14),
     &                  WORK(1,15),WORK(1,16),WORK(1,17),WORK(1,18),
     &                  WORK(1,19),WORK(1,20),WORK(1,21),WORK(1,22),
     &                  WORK(1,23),BWORK(1,1),BWORK(1,2),BWORK(1,3),
     &                  BWORK(1,4),BWORK(1,5),WORK(1,24))
C
C-----------------------------------------------------------------------
C EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
C-----------------------------------------------------------------------
C
      DO I=1,NDDON
       THDD_K(INDEX1(I)) = WORK(I,1)
       QDD_K(INDEX1(I)) = WORK(I,2)
       DTHBYDT_K(INDEX1(I)) = WORK(I,7)
       DTHBYDT_KM1(INDEX1(I)) = WORK(I,8)
       DQBYDT_K(INDEX1(I)) = WORK(I,9)
       DQBYDT_KM1(INDEX1(I)) = WORK(I,10)
       FLX_DD_K(INDEX1(I)) = WORK(I,11)
       RAIN_DD(INDEX1(I)) = WORK(I,22)
       SNOW_DD(INDEX1(I)) = WORK(I,23)
       BDD_START(INDEX1(I)) = BWORK(I,1)
       BDDWT_K(INDEX1(I)) = BWORK(I,2)
       BDDWT_KM1(INDEX1(I)) = BWORK(I,3)
       BDD_ON(INDEX1(I)) = BWORK(I,4)
       B_DD_END(INDEX1(I)) = BWORK(I,5)
      END DO
      END IF
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C DID NOT FORM
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (.NOT.BDD_ON(I).AND..NOT.B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C CARRY OUT CHANGE OF PHASE CALCULATION FOR PRECIPITATION FALLING
C THROUGH ENVIRONMENT
C-----------------------------------------------------------------------
C
         CALL CHG_PHSE (NPNTS,K,RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,
     &                  EXK,EXKM1,DELPKM1,THE_K,THE_KM1)
C
C-----------------------------------------------------------------------
C EVAPORATE RAIN FALLING THROUGH ENVIRONMENT IF LAYER K BELOW
C CLOUD BASE
C-----------------------------------------------------------------------
C
         CALL PEVP_BCB (NPNTS,K-1,ICCB,THE_KM1,P_KM1,QE_KM1,DELPKM1,
     &                  RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,DQBYDT_KM1,
     &                  EXKM1,TIMESTEP,CCA)
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C TERMINATES
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C UPDATE RAIN AND SNOW
C-----------------------------------------------------------------------
C
       IF (K.EQ.2) THEN
         DO I=1,NPNTS
           RAIN(I) = RAIN(I)+RAIN_DD(I)+RAIN_ENV(I)
           SNOW(I) = SNOW(I)+SNOW_DD(I)+SNOW_ENV(I)
         END DO
       END IF
C
      RETURN
      END
C
