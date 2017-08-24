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
CLL  SUBROUTINE THETAR-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE POTENTIAL TEMPERATURE OF THE DETRAINING
CLL            AIR IN LAYER K AND ALSO THE DIFFERENCE IN THE
CLL            WATER VAPOUR CONTENT OF THE DETRAINING AIR FROM THAT
CLL            OF THE MEAN PARCEL IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  4.1  6/6/96     Extra check added to ensure that
CLL                  negative values of parcel water
CLL                  content are not generated.
CLL                        Pete Inness.
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,
     *                   DQSK,BWKP1,EXK,PK)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(NPNTS)       ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                             !    IN LAYER K (K)
C
      REAL QEK(NPNTS)        ! IN ENVIRONMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
C
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K
                             !    (KG/KG)
C
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF THE
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    ENVIRONMENT OF LAYER K (KG/KG/K)
C
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS SATURATED IN LAYER K
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR POINTS AT WHICH CONDENSATE
                             !    IS LIQUID IN LAYER K+1
C
      REAL EXK(NPNTS)        ! IN EXNER RATIO FOR LEVEL K
C
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL XSQR(NPNTS)       ! OUT EXCESS WATER VAPOUR OF
                             !     DETRAINING AIR (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C
      REAL TT(NPNTS)         ! TEMPORARY TEMPERATURE FOR CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSAT
C
C*----------------------------------------------------------------------
C
      DO 20 I=1,NPNTS
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE POTENTIAL TEMPERATURE OF DETRAINING AIR
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (26)
CL----------------------------------------------------------------------
CL
       IF (.NOT.BGMK(I)) THEN
          THRK(I)=THEK(I) * (1. + C_VIRTUAL*QEK(I)) /
     *                 (1. + C_VIRTUAL*QPK(I))
       ELSE
          THRK(I) = THEK(I)*(1.0 + C_VIRTUAL*(QEK(I)-QSEK(I))/
     *                 (1.0 + C_VIRTUAL*THEK(I)*DQSK(I)))
       ENDIF
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE MIXING RATIO OF THE DETRAINING AIR AIR THE
CL  DIFFERENCE BETWEEN THIS AND THE MIXING RATIO OF THE MEAN
CL  PARCEL IN LAYER K
CL
CL  THE MOISTURE DIFFERENCE IS USED TO CALCULATE THE
CL  COND_DET_K TERM OF EQUATION (30), SECTION (6),
CL  UM DOCUMENTATIONM PAPER P27
CL----------------------------------------------------------------------
CL
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
       TT(I) = THRK(I)*EXK(I)
   20  CONTINUE
      CALL QSAT (XSQR,TT,PK,NPNTS)
C
      DO 30 I=1,NPNTS
CL----------------------------------------------------------------------
CL  SMALL NUMERICAL APPROXIMATIONS IN THE ABOVE CALCULATIONS CAN MEAN
CL  THAT THE DETRAINING PARCEL IS NO LONGER SATURATED AT THRK. ADD A
CL  CHECK TO SEE IF THE PARCEL IS STILL SATURATED, AND RESET BGMK TO
CL  FALSE IF IT IS NOT.
CL---------------------------------------------------------------------
         IF(XSQR(I).GT.QPK(I))BGMK(I)=.FALSE.


       IF (BGMK(I)) THEN
          QRK(I)  = XSQR(I)
          XSQR(I) = QPK(I) - XSQR(I)
       ELSE
          QRK(I)  = QPK(I)
          XSQR(I) = 0.
       ENDIF
   30  CONTINUE
C
      RETURN
      END
