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
CLL  SUBROUTINE TERMDD-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE WHETHER DOWNDRAUGHT IS ABLE TO CONTINUE
CLL
CLL            CALCULATE BUOYANCY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
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
      SUBROUTINE TERMDD (NPNTS,BDD_START,THDD_K,QDD_K,THE_K,QE_K,K,
     *                   B_DD_END,BDD_ON)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS USED IN THIS ROUTINE
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      REAL DET_LYR  ! THICKNESS LEVEL USED IN CALCULATION OF MIXING
                    ! DETRAINMENT FOR DOWNDRAUGHT  (PA)
      PARAMETER (DET_LYR = 10000.0)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
      INTEGER K                    ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(NPNTS)           ! IN MODEL POTENTIAL TEMPERATURE
                                   !    OF DOWNDRAUGHT AT LAYER K (K)
C
      REAL QDD_K(NPNTS)            ! IN MODEL MIXING RATIO OF
                                   !    DOWNDRAUGHT AT LAYER K
C
      REAL THE_K(NPNTS)            ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENTAL AIR IN LAYER K
C
      REAL QE_K(NPNTS)             ! IN MODEL MIXING RATIO AT LAYER K
C
      LOGICAL BDD_START(NPNTS)     ! IN MASK FOR THOSE POINTS WHERE
                                   !    DOWNDRAUGHT MAY OCCUR IN
                                   !    LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL B_DD_END(NPNTS)      ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT IS TERMINATING
C
      LOGICAL BDD_ON(NPNTS)        ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT CONTINUES TO LAYER
                                   !     K-1 (AS BDD_START HERE)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL BUOY1                   ! BUOYANCY OF PARCEL
C
      REAL THDD_V                  ! USED IN CALCULATION OF BUOYANCY
C
      REAL THE_V                   ! USED IN CALCULATION OF BUOYANCY
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT
C CAN CONTINUE TO NEXT LAYER
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
         THDD_V = THDD_K(I)*(1.0+C_VIRTUAL*QDD_K(I))
         THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
         BUOY1 = THDD_V - THE_V
C
C-----------------------------------------------------------------------
C CALCULATE STATE OF DOWNDRAUGHT
C-----------------------------------------------------------------------
C
         IF (BDD_START(I) .AND. BUOY1.GT.0.5) THEN
            BDD_ON(I) = .FALSE.
         ELSE IF (BUOY1.GT.0.5 .OR. K.EQ.2) THEN
            B_DD_END(I) = .TRUE.
         END IF
      END DO
C
      RETURN
      END
C
