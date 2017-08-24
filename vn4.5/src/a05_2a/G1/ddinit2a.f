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
CLL  SUBROUTINE DD_INIT------------------------------------------------
CLL
CLL  PURPOSE : ROUTINE TO INITIALISE THE DOWNDRAUGHT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
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
      SUBROUTINE DD_INIT(NPNTS,TH_UD_K,Q_UD_K,THE_K,QE_K,PK,EXK,THDD_K,
     &                   QDD_K,DELTD,DELQD,BDD_START,K,BDDI,BDD_ON)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! VECTOR LENGTH
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THE_K(NPNTS)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
C
      REAL TH_UD_K(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE OF
                                !    UPDRAUGHT, LAYER K (K)
C
      REAL QE_K(NPNTS)          ! IN MIXING RATIO OF ENVIRONMENT IN
                                !    LAYER K (KG/KG)
C
      REAL Q_UD_K(NPNTS)        ! IN PARCEL MIXING RATIO OF UPDRAUGHT,
                                !    LAYER K (KG/KG)
C
      REAL EXK(NPNTS)           ! IN EXNER RATIO OF LAYER K
C
      REAL PK(NPNTS)            ! IN PRESSURE OF LAYER K (PA)
C
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
C
      LOGICAL BDD_ON(NPNTS)     ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT IS ON
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(NPNTS)  ! INOUT
                                ! IN  MASK FOR THOSE POINT WHERE
                                !     DOWNDRAUGHT MAY START
                                ! OUT MASK FOR THOSE POINTS WHERE
                                !
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(NPNTS)        ! OUT DOWNDRAUGHT POTENTIAL TEMPERATURE
                                !     OF LAYER K
C
      REAL QDD_K(NPNTS)         ! OUT DOWNDRAUGHT MIXING RATIO OF
                                !     LAYER K
C
      REAL DELTD(NPNTS)         ! OUT COOLING NECESSARY TO ACHIEVE
                                !     SATURATION
C
      REAL DELQD(NPNTS)         ! OUT MOISTENING NECESSARY TO ACHIEVE
                                !     SATURATION
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL TH_MEAN(NPNTS)       ! MEAN POTENTIAL TEMPERATURE USED IN
                                ! CALCULATION OF SATURATED DOWNDRAUGHT
                                ! POTENTIAL TEMPERATURE IN LAYER K
C
      REAL Q_MEAN(NPNTS)        ! MEAN MIXING RATIO USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT
                                ! MIXING RATIO FOR LAYER K
C
      REAL T_MEAN(NPNTS)        ! MEAN TEMPERATURE USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE OF LAYER K (K)
C
      REAL THDDS(NPNTS)         ! SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE IN LAYER K (K)
C
      REAL QDDS(NPNTS)          ! SATURATED DOWNDRAUGHT MIXING RATIO
                                ! IN LAYER K (KG/KG)
C
C
      REAL BUOY                 ! BUOYANCY OF PARCEL IN LAYER K
C
      REAL THDD_V               ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! PARCEL IN LAYER K
C
      REAL THE_V                ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! ENVIRONMENT IN LAYER K
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL SATCAL
C
C-----------------------------------------------------------------------
C CALCULATE MEAN TEMPERATURE AND MIXING RATIO
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TH_MEAN(I) = (THE_K(I)+TH_UD_K(I))*0.5
       Q_MEAN(I) = (QE_K(I)+Q_UD_K(I))*0.5
       T_MEAN(I) = TH_MEAN(I)*EXK(I)
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE SATURATED DOWNDRAUGHT POTENTIAL TEMPERATURE FOR LAYER K
C-----------------------------------------------------------------------
C
      CALL SATCAL(NPNTS,T_MEAN,TH_MEAN,PK,QDDS,THDDS,K,EXK,Q_MEAN,
     *            THE_K)
C
C-----------------------------------------------------------------------
C IS SATURATED PARCEL NEGATIVELY BUOYANT COMPARED TO ENVIRONMENT
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (.NOT. BDD_ON(I) .AND. BDDI(I) .AND. K.GE.4) THEN
          THDD_V = THDDS(I)*(1.0+0.61*QDDS(I))
          THE_V = THE_K(I)*(1.0+0.61*QE_K(I))
          BUOY = THDD_V - THE_V
C
          IF (BUOY .LT. 0.5 ) THEN
C
C-----------------------------------------------------------------------
C INITIATE DOWNDRAUGHT
C-----------------------------------------------------------------------
C
             THDD_K(I) = THDDS(I)
             QDD_K(I) = QDDS(I)
             BDD_START(I) = .TRUE.
C
C-----------------------------------------------------------------------
C CALCULATE COOLING AND MOISTENING TO ACHIEVE SATURATION
C-----------------------------------------------------------------------
C
             DELTD(I) = THDDS(I)-THE_K(I)
             DELQD(I) = QDDS(I)-QE_K(I)
          END IF
       END IF
      END DO
C
      RETURN
      END
C
