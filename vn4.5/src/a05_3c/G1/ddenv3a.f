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
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.0   5/5/95    new deck added for vesrion 3A of convection
CLL                   scheme. Includes tracers and momentum in the
CLL                   downdraught.
CLL                   Pete Inness.
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
      SUBROUTINE DD_ENV (NPNTS,NP_FULL,THDD_K,THDD_KM1,QDD_K,QDD_KM1,
     *                   THE_K,THE_KM1,QE_K,QE_KM1,DTHBYDT_K,
     *                   DTHBYDT_KM1,DQBYDT_K,DQBYDT_KM1,FLX_DD_K,
     *                   FLX_DD_KM1,DELPK,DELPKM1,DELTD,DELQD,AMDETK,
     *                   EKM14,B_DD_END,BDD_START,BDD_ON,L_MOM,UDD_K,
     *                   VDD_K,UDD_KM1,VDD_KM1,UE_K,VE_K,UE_KM1,VE_KM1,
     *                   DUBYDT_K,DUBYDT_KM1,DVBYDT_K,DVBYDT_KM1,DELUD,
     *                   DELVD,EFLUX_U_DD,EFLUX_V_DD,
     *                   L_TRACER,NTRA,TRADD_K,TRADD_KM1,TRAE_K,
     *                   TRAE_KM1,DTRABYDT_K,DTRABYDT_KM1,DELTRAD)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS                 ! IN VECTOR LENGTH
C
      INTEGER NP_FULL               ! IN FULL VECTOR LENGTH
C
      INTEGER I,KTRA                ! LOOP COUNTERS
C
      INTEGER NTRA                  ! IN NUMBER OF TRACERS
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
      REAL UDD_K(NPNTS)             ! IN DOWNDRAUGHT U IN LAYER
                                    !    K (M/S)
C
      REAL UDD_KM1(NPNTS)           ! IN DOWNDRAUGHT U IN LAYER
                                    !     K-1 (M/S)
C
      REAL VDD_K(NPNTS)             ! IN DOWNDRAUGHT V IN LAYER
                                    !    K (M/S)
C
      REAL VDD_KM1(NPNTS)           ! IN DOWNDRAUGHT V IN LAYER
                                    !     K-1 (M/S)
C
      REAL TRADD_K(NP_FULL,NTRA)    ! IN DOWNDRAUGHT TRACER CONTENT
                                    !    IN LAYER K (KG/KG)
C
      REAL TRADD_KM1(NPNTS,NTRA)    ! IN DOWNDRAUGHT TRACER CONTENT
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
      REAL UE_K(NPNTS)              ! IN U AT LAYER K (M/S)
C
      REAL UE_KM1(NPNTS)            ! IN U AT LAYER K-1 (M/S)
C
      REAL VE_K(NPNTS)              ! IN V AT LAYER K (M/S)
C
      REAL VE_KM1(NPNTS)            ! IN V AT LAYER K-1 (M/S)
C
      REAL TRAE_K(NP_FULL,NTRA)     ! IN TRACER CONTENT OF
                                    !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL TRAE_KM1(NP_FULL,NTRA)   ! IN TRACER CONTENT OF ENVIRONMENT
                                    !    IN LAYER K-1 (KG/KG)
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
      REAL DELUD(NPNTS)             ! IN CHANGE TO ENVIRONMENT U DUE TO
                                    !    DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELVD(NPNTS)             ! IN CHANGE TO ENVIRONMENT V DUE TO
                                    !    DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELTRAD(NP_FULL,NTRA)    ! IN DEPLETION OF ENVIRONMENT TRACER
                                    !    DUE TO DOWNDRAUGHT FORMATION
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
      LOGICAL L_TRACER              ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM                 ! IN SWITCH FOR INCLUSION OF
                                    !    MOMENTUM TRANSPORTS
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
      REAL DQBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO
                                    !     AT LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO AT
                                    !     LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K-1 (KG/KG/S)
C
      REAL DUBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO U AT LAYER K
                                    !     (M/S**2)
                                    ! OUT UPDATED INCREMENT TO U
                                    !     AT LAYER K (M/S**2)
C
      REAL DUBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO U AT LAYER K-1
                                    !     (M/S**2)
                                    ! OUT UPDATED INCREMENT TO U
                                    !     AT LAYER K-1 (M/S**2)
C
      REAL DVBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO V AT LAYER K
                                    !     (M/S**2)
                                    ! OUT UPDATED INCREMENT TO V
                                    !     AT LAYER K (M/S**2)
C
      REAL DVBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO V AT LAYER K-1
                                    !     (M/S**2)
                                    ! OUT UPDATED INCREMENT TO V
                                    !     AT LAYER K-1 (M/S**2)
C
      REAL DTRABYDT_K(NP_FULL,NTRA) ! INOUT
                                    ! IN  INCREMENT TO TRACER CONTENT
                                    !     IN LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO TRACER
                                    !     CONTENT IN LAYER K (KG/KG/S)
C
      REAL DTRABYDT_KM1(NP_FULL,    ! INOUT
     *                  NTRA)       ! IN  INCREMENT TO TRACER CONTENT
                                    !     AT LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO TRACER
                                    !     CONTENT AT LAYER K-1
                                    !     (KG/KG/S)
C
      REAL EFLUX_U_DD(NPNTS),       ! INOUT
     *     EFLUX_V_DD(NPNTS)        ! IN  EDDY FLUX OF MOMENTUM DUE TO
                                    !     DD AT TOP OF LAYER
                                    ! OUT EDDY FLUX OF MOMENTUM DUE TO
                                    !     DD AT BOTTOM OF LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL TEMPRY                   ! USED IN CALCULATIONS OF THE
                                    ! EFFECT ON THE ENVIRONMENT
C
      REAL FLX_U_KM0P5,             ! EDDY FLUX OF U AND V AT
     *     FLX_V_KM0P5              ! BOTTOM OF LAYER
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
C
         END IF
C
       END IF
      END DO
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON MOMENTUM OF
C LAYER K
C-----------------------------------------------------------------------
C
      IF(L_MOM)THEN
C
      DO I=1,NPNTS
      IF(BDD_ON(I))THEN
       TEMPRY = FLX_DD_K(I)/DELPK(I)
       IF(BDD_START(I))THEN
         DUBYDT_K(I) = DUBYDT_K(I)-TEMPRY*DELUD(I)
         DVBYDT_K(I) = DVBYDT_K(I)-TEMPRY*DELVD(I)
       END IF
C
C----------------------------------------------------------------------
C CALCULATE MOMENTUM FLUX AT BOTTOM OF LAYER
C----------------------------------------------------------------------
C
       FLX_U_KM0P5 = FLX_DD_K(I) * (1.0-AMDETK(I)) * (1.0+EKM14(I)) *
     *               (UDD_K(I) - UE_KM1(I))
       FLX_V_KM0P5 = FLX_DD_K(I) * (1.0-AMDETK(I)) * (1.0+EKM14(I)) *
     *               (VDD_K(I) - VE_KM1(I))
C
C---------------------------------------------------------------------
C CALCULATE INCREMENT TO LAYER K
C---------------------------------------------------------------------
C
       IF (BDD_START(I)) THEN
C
        DUBYDT_K(I) = DUBYDT_K(I) - FLX_U_KM0P5/DELPK(I)
        DVBYDT_K(I) = DVBYDT_K(I) - FLX_V_KM0P5/DELPK(I)
C
C STORE FLUX READY FOR UPDATING NEXT LAYER
C
        EFLUX_U_DD(I) = FLX_U_KM0P5
        EFLUX_V_DD(I) = FLX_V_KM0P5
C
       ELSE
C
        DUBYDT_K(I) = DUBYDT_K(I) + ((EFLUX_U_DD(I) - FLX_U_KM0P5)/
     *                                                    DELPK(I))
        DVBYDT_K(I) = DVBYDT_K(I) + ((EFLUX_V_DD(I) - FLX_V_KM0P5)/
     *                                                    DELPK(I))
C
C STORE FLUX READY FOR UPDATING NEXT LAYER
C
        EFLUX_U_DD(I) = FLX_U_KM0P5
        EFLUX_V_DD(I) = FLX_V_KM0P5
C
       END IF
C
CL--------------------------------------------------------------------
CL TERMINAL DETRAINMENT OF MOMENTUM
CL-------------------------------------------------------------------
C
        IF(B_DD_END(I))THEN
           DUBYDT_KM1(I) = FLX_U_KM0P5/DELPKM1(I)
           DVBYDT_KM1(I) = FLX_V_KM0P5/DELPKM1(I)
C
C ZERO FLUX AT BOTTOM OF LAYER
C
           EFLUX_U_DD(I) = 0.0
           EFLUX_V_DD(I) = 0.0
        END IF
      END IF
      END DO
C
      END IF
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON TRACER CONTENT OF
C LAYER K
C-----------------------------------------------------------------------
C
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
       DO I=1,NPNTS
        IF(BDD_ON(I))THEN
C
        TEMPRY = FLX_DD_K(I)/DELPK(I)
        IF(BDD_START(I))THEN
        DTRABYDT_K(I,KTRA) = DTRABYDT_K(I,KTRA)-TEMPRY*DELTRAD(I,KTRA)
        END IF
        DTRABYDT_K(I,KTRA) = DTRABYDT_K(I,KTRA) + TEMPRY * (
     *
     *   (1.0+EKM14(I)) * (1.0-AMDETK(I)) *             ! COMPENSATING
     *   (TRAE_KM1(I,KTRA)-TRAE_K(I,KTRA))              ! SUBSIDENCE
     *        +
     *   AMDETK(I)* (TRADD_K(I,KTRA)-TRAE_K(I,KTRA))    ! MIXING
     *        )                                         ! DETRAINMENT
C
C-----------------------------------------------------------------------
C TERMINAL DETRAINMENT OF TRACER
C-----------------------------------------------------------------------
C
           IF(B_DD_END(I))THEN
             TEMPRY = FLX_DD_KM1(I)/DELPKM1(I)
             DTRABYDT_KM1(I,KTRA)=DTRABYDT_KM1(I,KTRA)+TEMPRY*
     *                            (TRADD_KM1(I,KTRA)-TRAE_KM1(I,KTRA))
           END IF
C
        END IF
       END DO
      END DO
C
      END IF
C
      RETURN
      END
C
