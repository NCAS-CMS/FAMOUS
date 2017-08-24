C ******************************COPYRIGHT******************************
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
CLL  SUBROUTINE LAYER_DD--------------------------------------------
CLL
CLL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
CLL            -PRESSURE
CLL            -LAYER THICKNESS
CLL            -ENTRAINMENT COEFFICIENTS
CLL            -DETRAINMENT COEFFICIENTS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  MODEL            MODIFICATION HISTORY::
CLL VERSION  DATE
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with LAYERD2A.
!LL   4.5   18/02/98  Correct code to match LAYERD2A and call
!LL                   comdecks. D. Robinson.
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
      SUBROUTINE LAYER_DD(NPNTS,K,KCT,THE_K,THE_KM1,FLX_STRT,AK,
     *                    BK,AKM12,BKM12,DELAK,DELBK,EXNER_KM12,
     *                    EXNER_KP12,EXNER_KM32,PSTAR,PK,PKM1,DELPK,
     *                    DELPKM1,EXK,EXKM1,AMDETK,EKM14,EKM34,KMIN,
     *                    BDDI,recip_pstar)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
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

      REAL AE1,AE2,       ! COEFFICIENTS USED IN CALCULATION
     *     ENTCOEF,       ! OF ENTRAINMENT RATE
     *     SH_FAC
C
      PARAMETER(AE1=1.0,AE2=1.5)
      PARAMETER (ENTCOEF = 3.0,SH_FAC=1.0)
C
      REAL DDCOEF1, ! COEFFICIENTS USED IN CALCULATION OF DOWNDRAUGHT
     &     DDCOEF2  ! ENTRAINMENT RATES
      PARAMETER (DDCOEF1 = 1.8E6, DDCOEF2 = 3.0)
C
      REAL DET_LYR  ! THICKNESS LEVEL USED IN CALCULATION OF MIXING
                    ! DETRAINMENT FOR DOWNDRAUGHT  (PA)
      PARAMETER (DET_LYR = 10000.0)
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTER
C----------------------------------------------------------------------
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
      INTEGER I                 ! COUNTER FOR DO LOOPS
C
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP LAYER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL AK(K)                ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(K)                ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(K+1)           ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(K+1)           ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(K)             ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(K)             ! IN ) FOR THICKNESS OF LAYER K
C
      REAL PSTAR(NPNTS)         ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER_KM12(NPNTS)    ! IN EXNER FUNCTION AT LAYER K-1/2
C
      REAL EXNER_KP12(NPNTS)    ! IN EXNER FUNCTION AT LAYER K+1/2
C
      REAL EXNER_KM32(NPNTS)    ! IN EXNER FUNCTION AT LAYER K-3/2
C
      REAL FLX_STRT(NPNTS)      ! IN UPDRAUGHT MASSFLUX AT LEVEL WHERE
                                !    DOWNDRAUGHT STARTS (PA/S)
C
      REAL THE_K(NPNTS)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(NPNTS)       ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K-1 (K)
C
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR POINTS WHERE DOWNDRAUGHT
                                !    MAY INITIATE
      REAL recip_PSTAR(NPNTS)! Reciprocal of pstar array
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      INTEGER KMIN(NPNTS)       ! INOUT
                                ! FREEZING LEVEL
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PK(NPNTS)            ! OUT PRESSURE AT LAYER K (PA)
C
      REAL PKM1(NPNTS)          ! OUT PRESSURE AT LAYER K-1 (PA)
C
      REAL DELPK(NPNTS)         ! OUT THICKNESS OF LAYER K (PA)
C
      REAL DELPKM1(NPNTS)       ! OUT THICHNESS OF LAYER K-1 (PA)
C
      REAL EKM14(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EKM34(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(NPNTS)        ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EXK(NPNTS)           ! OUT EXNER FUNCTION AT LEVEL K
C
      REAL EXKM1(NPNTS)         ! OUT EXNER FUNCTION AT LEVEL K-1
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL TTK                  ! TEMPERATURE STORE AT LAYER K
C
      REAL TTKM1                ! TEMPERATURE STORE AT LAYER K-1
C
      REAL THKM12               ! POTENTIAL TEMPERATURE STORE AT
                                ! LAYER K-1/2
C
      REAL TTKM12               ! TEMPERATURE STORE AT LAYER K-1/2
C
      REAL INCR_FAC             ! INCREMENT FACTOR FOR ENTRAINMENT
                                ! RATES AT FREEZING LEVEL
C
      REAL
     &    PU,PL
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------



C----------------------------------------------------------------------
C SET KMIN TO INITIAL VALUE
CL CALCULATE PK, DELPK AND EXNER FUNCTION - IF K = KCT THEN
CL VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
C----------------------------------------------------------------------
C
      IF (K.EQ.KCT+1) THEN
       DO I=1,NPNTS
        KMIN(I) = KCT+2
        PK(I) = AK(K) + BK(K)*PSTAR(I)
        DELPK(I) = - DELAK(K) - DELBK(K)*PSTAR(I)
        PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        PL=PSTAR(I)*BKM12(K) + AKM12(K)
        EXK(I) = P_EXNER_C(EXNER_KP12(I),EXNER_KM12(I),PU,PL,KAPPA)
       END DO
      ELSE
       DO I=1,NPNTS
        PK(I) = PKM1(I)
        DELPK(I) = DELPKM1(I)
        EXK(I) = EXKM1(I)
       END DO
      END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE PKM1, DELPKM1
CL CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K-1, AND
CL DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
CL---------------------------------------------------------------------
CL
      DO I=1,NPNTS
        PKM1(I) = AK(K-1) + BK(K-1)*PSTAR(I)
        DELPKM1(I) = - DELAK(K-1) - DELBK(K-1)*PSTAR(I)
        PU=PSTAR(I)*BKM12(K) + AKM12(K)
        PL=PSTAR(I)*BKM12(K-1) + AKM12(K-1)
        EXKM1(I) = P_EXNER_C(EXNER_KM12(I),EXNER_KM32(I),PU,PL,KAPPA)
C
CL
CL---------------------------------------------------------------------
CL CALCULATE FREEZING LEVEL : CHECK IF FREEZING LEVEL IN THIS LAYER
CL---------------------------------------------------------------------
CL
       IF (KMIN(I).EQ.KCT+2) THEN
        TTK = THE_K(I)*EXK(I)
        TTKM1 = THE_KM1(I)*EXKM1(I)
        THKM12 = (THE_KM1(I)+THE_K(I))*0.5
        TTKM12 = THKM12*EXNER_KM12(I)
        IF (TTKM12 .GE. TM .AND. TTK .LT. TM) THEN
           KMIN(I) = K
        ELSE IF (TTKM1 .GE. TM .AND. TTKM12 .LT. TM) THEN
           KMIN(I) = K-1
        END IF
       END IF
C
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER 27
CL SECTION (2C), EQUATION(14)
CL---------------------------------------------------------------------
CL
      IF (PK(I).LT.PSTAR(I)-DET_LYR) THEN
        EKM14(I) = AE2 * (AKM12(K)+BKM12(K)*PSTAR(I)-PK(I)) *
     &                                                 recip_PSTAR(I)
        EKM34(I) = AE2 * (PKM1(I)-AKM12(K)-BKM12(K)*PSTAR(I)) *
     &                                                 recip_PSTAR(I)
        AMDETK(I) = (EKM14(I)+EKM34(I)) * (1.0-1.0/AE2)
      ELSE
        EKM14(I) = 0.0
        EKM34(I) = 0.0
        AMDETK(I) = DELPK(I) / (PSTAR(I)*(1.0-BKM12(K+1))-AKM12(K+1))
      END IF
C
      IF (BDDI(I)) THEN
C
      IF (K.EQ.KMIN(I) .AND. PK(I).LT.PSTAR(I)-DET_LYR) THEN
        INCR_FAC = FLX_STRT(I)*DDCOEF1*recip_pstar(I)
        IF (INCR_FAC.GT.6.0) INCR_FAC=6.0
        EKM14(I) = EKM14(I)*INCR_FAC
        EKM34(I) = EKM34(I)*INCR_FAC
      ELSE
        EKM14(I) = EKM14(I)*DDCOEF2
        EKM34(I) = EKM34(I)*DDCOEF2
        IF (KMIN(I).NE.KCT+2 .AND. K.LT.KMIN(I) .AND. PK(I).LT.
     * PSTAR(I)-DET_LYR)  AMDETK(I) = AMDETK(I)*DDCOEF2
      END IF
C
      END IF
      END DO
C
      RETURN
      END
C
