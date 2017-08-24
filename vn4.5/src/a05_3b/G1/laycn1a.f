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
CLL  SUBROUTINE LAYER_CN-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
CLL            -PRESSURE
CLL            -LAYER THICKNESS
CLL            -ENTRAINMENT COEFFICIENTS
CLL            -DETRAINMENT COEFFICIENTS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL C.W. , D.G. <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.3  Feb. 97   T3E optimisation: introduce recip_pstar to 
CLL                  eliminate divisions by pstar.      S.J.Swarbrick
!LL  4.5   20/02/98  Remove redundant code. A. Dickinson
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED:P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LAYER_CN(K,NP_FIELD,NPNTS,NLEV,EXNER,AK,BK,AKM12,BKM12,
     *                    DELAK,DELBK,PSTAR,PK,PKP1,DELPK,DELPKP1,
     *                    DELPKP12,EKP14,EKP34,AMDETK,EXK,EXKP1,
     *                    DELEXKP1,recip_PSTAR)    
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      REAL AE1,AE2,       ! COEFFICIENTS USED IN CALCULATION
     *     ENTCOEF,       ! OF ENTRAINMENT RATE
     *     SH_FAC
C
      PARAMETER(AE1=1.0,AE2=1.5)
      PARAMETER (ENTCOEF = 3.0,SH_FAC=1.0)
C
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
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTER
C----------------------------------------------------------------------
C
      INTEGER NP_FIELD          ! IN FULL LENGTH OF DATA
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
      INTEGER NLEV              ! IN NUMBER OF MODEL LEVELS
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
      INTEGER I                 ! COUNTER FOR DO LOOPS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL AK(NLEV)             ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(NLEV)             ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(NLEV+1)        ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(NLEV+1)        ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(NLEV)          ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(NLEV)          ! IN ) FOR THICKNESS OF LAYER K
C
      REAL PSTAR(NP_FIELD)      ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER FUNCTION AT LAYER
                                  ! BOUNDARIES STARTING AT LEVEL K-1/2
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PK(NPNTS)            ! OUT PRESSURE AT LAYER K (PA)
C
      REAL PKP1(NPNTS)          ! OUT PRESSURE AT LAYER K+1 (PA)
C
      REAL DELPK(NPNTS)         ! OUT THICKNESS OF LAYER K (PA)
C
      REAL DELPKP1(NPNTS)       ! OUT THICHNESS OF LAYER K+1 (PA)
C
      REAL DELPKP12(NPNTS)      ! OUT THICKNESS BETWEEN LAYER K AND K+1
                                !     (PA)
C
      REAL EKP14(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EKP34(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(NPNTS)        ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EXK(NPNTS)           ! EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1(NPNTS)         ! EXNER FUNCTION AT LEVEL K+1
C
      REAL DELEXKP1(NPNTS)      ! DIFFERENCE IN EXNER FUNCTION
                                ! BETWEEN K+3/2 AND K+1/2
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL AEKP14,AEKP34        ! USED IN CALCULATION OF ENTRAINMENT
                                ! RATE
C

      REAL
     &    PU,PL
      REAL recip_PSTAR(NP_FIELD)  ! Reciprocal of pstar array 
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


C*---------------------------------------------------------------------
C
C----------------------------------------------------------------------
C SET CONSTANT AE USED IN CALCULATION OF ENTARINMENT AND
C DETRAINMENT RATES DEPENDING UPON LEVEL
C----------------------------------------------------------------------
C
      IF(K.EQ.1)THEN
        AEKP14 = AE1
        AEKP34 = AE2
      ELSE
        AEKP14 = AE2
        AEKP34 = AE2
      END IF
C
      DO 10 I=1,NPNTS
CL
CL---------------------------------------------------------------------
CL CALCULATE PK AND DELPK - IF K = 1 (LOWEST MODEL LAYER) THEN
CL VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          PK(I) = AK(K) + BK(K)*PSTAR(I)
          DELPK(I) = -DELAK(K) - DELBK(K)*PSTAR(I)
        ELSE
          PK(I) = PKP1(I)
          DELPK(I) = DELPKP1(I)
        END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE PKP1, DELPKP1 AND DELPK+1/2
CL---------------------------------------------------------------------
CL
        PKP1(I) = AK(K+1) + BK(K+1)*PSTAR(I)
        DELPKP1(I) = -DELAK(K+1) - DELBK(K+1)*PSTAR(I)
        DELPKP12(I) = PK(I) - PKP1(I)
CL
CL---------------------------------------------------------------------
CL CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K+1, AND
CL DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
          PL=PSTAR(I)*BKM12(K) + AKM12(K)
          EXK(I) = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        ELSE
          EXK(I) = EXKP1(I)
        END IF
        PU=PSTAR(I)*BKM12(K+2) + AKM12(K+2)
        PL=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        EXKP1(I) = P_EXNER_C(EXNER(I,K+2),EXNER(I,K+1),PU,PL,KAPPA)
        DELEXKP1(I) = EXNER(I,K+1)-EXNER(I,K+2)
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT AND MIXING DETRAINMENT COEFFICIENTS
CL---------------------------------------------------------------------
CL
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2C), EQUATION(14)
CL---------------------------------------------------------------------
CL
        EKP14(I) = ENTCOEF * AEKP14 * PK(I) *                           
     *             (PK(I) - AKM12(K+1) - BKM12(K+1)*PSTAR(I)) *         
     *              recip_PSTAR(I)*recip_PSTAR(I)   
        EKP34(I) = ENTCOEF * AEKP34 * (AKM12(K+1)+BKM12(K+1)*PSTAR(I)) *
     *             (AKM12(K+1) + BKM12(K+1)*PSTAR(I) - PKP1(I)) *       
     *              recip_PSTAR(I)*recip_PSTAR(I)    
CL
CL---------------------------------------------------------------------
CL CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2C), EQUATION(15)
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          AMDETK(I) = 0.0
        ELSE
          AMDETK(I) = (EKP14(I) + EKP34(I)) * (1.0-1.0/AEKP34)
        END IF
 10   CONTINUE
C
      RETURN
      END
