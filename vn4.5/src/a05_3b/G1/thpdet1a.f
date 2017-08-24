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
CLL  SUBROUTINE THP_DET------------------------------------------------
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED : P27
CLL
CLL  PURPOSE : CALCULATES POTENTIAL TEMPERATURE OF THE
CLL            PARCEL IN LAYER K+1 AFTER FORCED DETRAINMENT
CLL            IN LAYER K
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (6), EQUATION (28)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,
     *                    DQSKP1,BGMKP1,BCALC)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      REAL XSBMIN !  MINIMUM EXCESS BUOYANCY TO CONTINUE PARCEL ASCENT
                  !  (K)
      PARAMETER (XSBMIN = 0.2)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS           ! IN VECTOR LENGTH
C
      INTEGER I               ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARAIBLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THEKP1(NPNTS)      ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                              !    IN LAYER K+1 (K)
C
      REAL QPKP1(NPNTS)       ! IN PARCEL MIXING RATIO IN LAYER K+1
                              !    (KG/KG)
C
      REAL QSEKP1(NPNTS)      ! IN ENVIRONMENT SATURATED MIXING RATIO
                              !    IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(NPNTS)      ! IN GRADIENT OF SATURATION MIXING RATIO
                              !    POTENTIAL TEMPERATURE FOR THE
                              !    ENVIRONMENT IN LAYER K+1 (KG/KG/K)
C
      REAL QEKP1(NPNTS)       ! IN ENVIRONMENT MIXING RATIO IN
                              !    LAYER K+1 (KG/KG)
C
      LOGICAL BGMKP1(NPNTS)   ! IN MASK FOR PARCELS WHICH ARE SATURATED
                              !    IN LAYER K+1
C
      LOGICAL BCALC(NPNTS)    ! IN MASK FOR PARCELS AT WHICH
                              !    CALCULATIONS OF THIS SUBROUTINE ARE
                              !    TO BE CARRIED OUT
C
C
C-----------------------------------------------------------------------
C VARAIBLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THPKP1(NPNTS)      ! OUT PARCEL POTENTIAL TEMPERATURE
                              !     IN LAYER K+1 AFTER FORCED
                              !     DETRAINMENT (K)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
C
      DO 10 I=1,NPNTS
        IF (BCALC(I))THEN
         IF (BGMKP1(I)) THEN
           THPKP1(I) = THEKP1(I) +
     *                (C_VIRTUAL*THEKP1(I)*
     *                           (QEKP1(I)-QSEKP1(I)) + XSBMIN)
     *               /( 1. + C_VIRTUAL*THEKP1(I)*DQSKP1(I) )
C
         ELSE
           THPKP1(I) = (THEKP1(I)*(1. + C_VIRTUAL*QEKP1(I))
     *                                                        + XSBMIN)
     *                    /(1. + C_VIRTUAL*QPKP1(I))
         END IF
        END IF
  10  CONTINUE
C
      RETURN
      END
