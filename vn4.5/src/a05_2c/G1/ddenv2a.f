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
CLL  SUBROUTINE DD_ENV-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE THE EFFECT OF THE DOWNDRAUGHT
CLL            ON THE LARGE_SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
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
      SUBROUTINE DD_ENV (NPNTS,THDD_K,THDD_KM1,QDD_K,QDD_KM1,THE_K,
     *                   THE_KM1,QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,
     *                   DQBYDT_K,DQBYDT_KM1,FLX_DD_K,FLX_DD_KM1,DELPK,
     *                   DELPKM1,DELTD,DELQD,AMDETK,EKM14,
     *                   B_DD_END,BDD_START,BDD_ON)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS                 ! IN VECTOR LENGTH
C
      INTEGER I                     ! LOOP COUNTER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(NPNTS)            ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K (K)
C
      REAL THDD_KM1(NPNTS)          ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K-1 (K)
C
      REAL QDD_K(NPNTS)             ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K (KG/KG)
C
      REAL QDD_KM1(NPNTS)           ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K-1 (KG/KG)
C
      REAL THE_K(NPNTS)             ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(NPNTS)           ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_K(NPNTS)              ! IN MIXING RATIO AT LAYER K (KG/KG)
C
      REAL QE_KM1(NPNTS)            ! IN MIXING RATIO AT LAYER K-1
                                    !    (KG/KG)
C
      REAL FLX_DD_K(NPNTS)          ! IN MASS FLUX IN LAYER K (PA/S)
C
      REAL FLX_DD_KM1(NPNTS)        ! IN MASS FLUX IN LAYER K-1 (PA/S)
C
      REAL DELPK(NPNTS)             ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
C
      REAL DELPKM1(NPNTS)           ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
C
      REAL DELTD(NPNTS)             ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
C
      REAL DELQD(NPNTS)             ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
C
      REAL AMDETK(NPNTS)            ! IN MIXING DETRAINMENT AT LEVEL K
                                    !    MULTIPLIED BY APPROPRIATE LAYER
                                    !    THICKNESS
C
      REAL EKM14(NPNTS)             ! IN EXNER RATIO AT LAYER K-1/4
C
      LOGICAL B_DD_END(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS TERMINATING
C
      LOGICAL BDD_START(NPNTS)      ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS STARTING
C
      LOGICAL BDD_ON(NPNTS)         ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS ON
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHBYDT_K(NPNTS)         ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
C
      REAL DTHBYDT_KM1(NPNTS)       ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
C
      REAL DQBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO AT
                                    !     LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K-1 (KG/KG/S)
C
      REAL DQBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO
                                    !     AT LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K (KG/KG/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL TEMPRY                   ! USED IN CALCULATIONS OF THE
                                    ! EFFECT ON THE ENVIRONMENT
C
C-----------------------------------------------------------------------
C CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (BDD_ON(I)) THEN
C
C-----------------------------------------------------------------------
C SUBTRACT ENERGY USED TO FORM DOWNDRAUGHT
C-----------------------------------------------------------------------
C
       TEMPRY = FLX_DD_K(I)/DELPK(I)
       IF (BDD_START(I)) THEN
         DTHBYDT_K(I) = DTHBYDT_K(I)-TEMPRY*DELTD(I)
         DQBYDT_K(I) = DQBYDT_K(I)-TEMPRY*DELQD(I)
       END IF
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON POTENTIAL TEMPERATURE OF
C LAYER K
C-----------------------------------------------------------------------
C
       DTHBYDT_K(I) = DTHBYDT_K(I) + TEMPRY * (
     *
     *          (1.0+EKM14(I)) * (1.0-AMDETK(I)) *      ! COMPENSATING
     *           (THE_KM1(I)-THE_K(I))                  ! SUBSIDENCE
     *        +
     *          AMDETK(I)* (THDD_K(I)-THE_K(I))         ! MIXING
     *        )                                         ! DETRAINMENT
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON MIXING RATIO OF
C LAYER K
C-----------------------------------------------------------------------
C
       DQBYDT_K(I) = DQBYDT_K(I) + TEMPRY * (
     *
     *      (1.0+EKM14(I)) * (1.0-AMDETK(I)) *       ! COMPENSATING
     *      (QE_KM1(I)-QE_K(I))                      ! SUBSIDENCE
     *    +
     *      AMDETK(I)* (QDD_K(I)-QE_K(I))            ! MIXING
     *    )                                          ! DETRAINMENT
C
C-----------------------------------------------------------------------
C TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
C-----------------------------------------------------------------------
C
         IF (B_DD_END(I)) THEN
           TEMPRY = FLX_DD_KM1(I)/DELPKM1(I)
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)+TEMPRY*
     *                      (THDD_KM1(I)-THE_KM1(I))
           DQBYDT_KM1(I) = DQBYDT_KM1(I)+TEMPRY*(QDD_KM1(I)-QE_KM1(I))
         END IF
C
       END IF
      END DO
C
      RETURN
      END
C
