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
CLL  SUBROUTINE FLX_INIT-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY SUMMER 1992
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG020893 : DUE TO CHANGE IN WAY CLOUD TOP IS
CLL                               ESTIMATED BECAUSE OF CHANGES TO THE
CLL                               CALCULATION OF FORCED DETRAINMENT
CLL
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
      SUBROUTINE FLX_INIT (NPNTS,KCT,ICCB,ICCT,FLX,FLX_DD_K,BDDI,
     *                     FLX_STRT)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! IN NUMBER OF POINTS
C
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      INTEGER ICCB(NPNTS)       ! IN CONVECTIVE CLOUD BASE
C
      INTEGER ICCT(NPNTS)       ! IN CONVECTIVE CLOUD TOP
C
      REAL FLX(NPNTS,KCT+1)     ! IN CONVECTIVE MASSFLUX (PA/S)
C
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL FLX_DD_K(NPNTS)      ! OUT DOWNDRAUGHT MASSFLUX OF LAYER K
                                !     (PA/S)
C
      REAL FLX_STRT(NPNTS)      ! OUT UPDRAUGHT MASSFLUX AT LEVEL
                                !     DOWNDRAUGHT STARTS (PA/S)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      INTEGER KDDREF            ! REFERENCE LEVEL FOR DOWNDRAUGHT
                                ! MASSFLUX
C
C----------------------------------------------------------------------
C CALCULATE DOWNDRAUGHT MASSFLUX BASED ON A REFERENCE LEVEL WHICH IS
C 3/4 CLOUD DEPTH
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (BDDI(I)) THEN
          KDDREF = ICCB(I) + 0.75*(ICCT(I) - ICCB(I))
          IF (KDDREF .GE. ICCT(I)-1) KDDREF=ICCT(I)-1
          FLX_STRT(I) = FLX(I,KDDREF)
          FLX_DD_K(I) = FLX_STRT(I) * 0.05
       END IF
      END DO
C
      RETURN
      END
C
