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
CLL  SUBROUTINE DQSDTH-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES GARDIENT OF SATURATION MIXING RATIO
CLL            WITH POTENTIAL TEMPERATURE FORM THE
CLL            CLAUSIUS-CLAPEYRON EQUATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Removed divide
!LL                                                   D.Salmond
CLL
CLL  PROGRAMMING STANDARDS :
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL                  SECTION(4), EQUATION (20)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DQS_DTH (DQS,K,THEK,QSEK,EXK,NPNTS)
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      REAL RV  !  GAS CONSTANT FOR WATER VAPOUR (J/KG/K)
      PARAMETER ( RV = 461.1 )
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS        ! IN VECTOR LENGTH
C
      INTEGER K            ! IN PRESENT MODEL LAYER
C
      INTEGER I            ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(NPNTS)     ! IN POTENTIAL TEMPERATURE FOR LAYER K (K)
C
      REAL QSEK(NPNTS)     ! IN SATURATION MIXING RATIO FOR LAYER K (K)
C
      REAL EXK(NPNTS)      ! IN EXNER RATIO FOR LAYER K
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL DQS(NPNTS)      ! OUT GRADIENT OF SATURATION MIXING RATIO
                           !     WITH POTENTIAL TEMPERATURE FOR LAYER K
                           !     (KG/KG/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EL              ! LATENT HEATING OF CONDENSATION OR
                           ! (CONDENSATION PLUS HEATING) (J/KG)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C  CREATE A VECTOR OF LATENT HEATS ACCORDING TO WHETHER QSAT IS WRT
C  ICE OR WATER
C-----------------------------------------------------------------------
C
      IF (THEK(I)*EXK(I) .GT. 273.) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
C
C-----------------------------------------------------------------------
C CALCULATE D(QSAT)/D(THETA)
C-----------------------------------------------------------------------
C
       DQS(I) = EL*QSEK(I)/(EXK(I)*RV*THEK(I)*THEK(I))
   10  CONTINUE
C
      RETURN
      END
