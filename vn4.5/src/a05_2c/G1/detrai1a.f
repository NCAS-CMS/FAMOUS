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
CLL  SUBROUTINE DETRAIN------------------------------------------------
CLL
CLL  PURPOSE : FORCED DETRAINMENT CALCULATION
CLL
CLL            SUBROUTINE THP_DET CALCULATES THE POTENTIAL
CLL            TEMPERATURE OF THE PARCEL IN LAYER K+1
CLL            AFTER FORCED DETRAINMENT
CLL
CLL            SUBROUTINE THETAR CALCULATES THE POTENTIAL TEMPERATURE
CLL            OF THE AIR IN LAYER K UNDERGOING FORCED DETRAINMENT
CLL
CLL            SUBROUTINE DET_RATE CALCULATES THE FORCED DETRAINMENT
CLL            RATE OF THE ENSEMBLE IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93  : DG020893 : TO MAKE CALCULATIONS OF FORCED
CLL                                DETRAINMENT RATE LESS PRONE TO
CLL                                FAILURE
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DETRAIN (NPNTS,THEK,QEK,THPK,QPK,QSEK,DQSK,BGMK,
     *                     THEKP1,QEKP1,THPKP1,QPKP1,QSEKP1,DQSKP1,
     *                     BGMKP1,BWKP1,XSQKP1,
     *                     DELTAK,THRK,QRK,EKP14,EKP34,PK,PKP1,
     *                     EXK,EXKP1)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
      INTEGER NREDO          ! NUMBER OF POINTS FOR WHICH FORCED
                             ! DETRAINMENT CALCULATION MUST BE
                             ! AS THE PROCESSES EITHER CAUSES THE
                             ! PARCEL TO BECOME SATURATED OR
                             ! SUB-SATURATED
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      REAL THEK(NPNTS)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(NPNTS)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
C
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT LEVEL K+1
C
      REAL EXK(NPNTS)        ! IN EXNER RATIO AT LEVEL K
C
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THPKP1(NPNTS)     ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (K)
C
      REAL QPKP1(NPNTS)      ! INOUT
                             ! IN  PARCEL MIXING RATIO IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (KG/KG)
C
      REAL XSQKP1(NPNTS)     ! INOUT
                             ! IN  EXCESS WATER IN PARCEL AFTER
                             !     LIFTING FROM LAYER K TO K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
                             ! OUT EXCESS WATER IN PARCEL IN LAYER
                             !     K+1 AFTER FORCED DETRAINMENT
                             !     (KG/KG)
C
      LOGICAL BGMKP1(NPNTS)  ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     FORCED DETRAINMENT
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL DELTAK(NPNTS)     ! OUT PARCEL FORCED DETRAINMENT RATE
                             !     IN LAYER K
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
      LOGICAL BDETK(NPNTS)   ! MASK FOR PARCELS WHICH ARE
                             ! UNDERGOING FORCED DETRAINMENT
                             ! IN THEIR ASCENT FROM LAYER K
                             ! TO K+1
C
      REAL XSQR(NPNTS)       ! EXCESS PARCEL WATER VAPOUR
                             ! DURING DETRAINMENT (KG/KG)
C
      REAL THPKP1W(NPNTS) ,  ! TEMPORARY STOREAGE FOR PARCEL
     *     QPKP1W(NPNTS) ,   ! POTENTIAL TEMPERATURE (K), MIXING
     *     XSQK1W(NPNTS)     ! RATIO (KG/KG), EXCESS WATER VAPOUR
      LOGICAL BGKP1W(NPNTS)  ! (KG/KG) AND MASK FOR SATURATION
                             ! IN LAYER K+1
C
      LOGICAL BRECAL(NPNTS)  ! MASK FOR THOSE POINTS AT WHICH THE
                             ! THE DETRAINMENT CALCULATION NEEDS
                             ! REPEATING
C
      REAL TT(NPNTS)         ! TEMPORARY STORE FOR TEMPERATURE
                             ! FOR THE CALCULATION OF SATURATED
                             ! MIXING RATIO (K)
C
      REAL EPSS              ! (1+EKP14)*(1+EKP34)
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL THP_DET,QSAT,THETAR,DET_RATE
C
C*---------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
C
C----------------------------------------------------------------------
C AT START OF ROUTINE FORCED DETARINMENT DONE AT ALL POINTS SO
C SET ARRAY BDETK EQUAL TO .TRUE.
C SET FORCED DETRAINMENT RATE EQUAL TO ZERO
C----------------------------------------------------------------------
C
       BDETK(I) = .TRUE.
       DELTAK(I) = 0.0
C
C-----------------------------------------------------------------------
C   SAVE THE CURRENT VALUES OF QPKP1, XSQKP1 AND BGMKP1
C-----------------------------------------------------------------------
C
       THPKP1W(I) = THPKP1(I)
       QPKP1W(I) = QPKP1(I)
       XSQK1W(I) = XSQKP1(I)
       BGKP1W(I) = BGMKP1(I)
C
C-----------------------------------------------------------------------
C   ADD THE EXCESS WATER VAPOUR BACK INTO THE DETRAINING PARCELS
C-----------------------------------------------------------------------
C
       QPKP1(I) = QPKP1(I) + XSQKP1(I)
   10 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE IN LAYER K+1
CL  AT THE POINTS WHERE DETRAINMENT IS TAKING PLACE
CL
CL  SUBROUTINE THP_DET
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (28)
CL----------------------------------------------------------------------
CL
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     *              BGMKP1,BDETK)
CL
CL---------------------------------------------------------------------
CL  CHECK TO SEE IF SUFFICIENT EXCESS WATER VAPOUR IN THE
CL  INITIAL DRY ASCENT TO ALLOW PARCEL TO BE SATURATED
CL  IN LAYER K+1 AFTER FORCED DETRAINMENT
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (29)
CL
CL  NOTE : ONLY ALLOW PARCEL TO BE SATURATED IN LAYER K+1 IF
CL         SATURATED INITIALLY.  IT IS POSSIBLE FOR SMALL
CL         SUPERSATURATIONS TO IF SUBROUTINE LATENT_H CAUSES
CL         PARCEL TO BE COME UNSATURATED.  IN THIS CASE TREAT
CL         THE PARCEL AS UNSATURATED IN LAYER K+1
CL---------------------------------------------------------------------
CL
C
C-----------------------------------------------------------------------
C   CALCULATE THE EXCESS WATER VAPOUR IN LAYER K+1 AND RECALCULATE
C   BGMKP1 AND QPKP1.
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
      DO 25 I = 1,NPNTS
       TT(I) = THPKP1(I)*EXKP1(I)
   25 CONTINUE
      CALL QSAT (XSQKP1,TT,PKP1,NPNTS)
C
      DO 30 I=1,NPNTS
       XSQKP1(I) = QPKP1(I) - XSQKP1(I)
C
       BRECAL(I) = BGMKP1(I)
C
C----------------------------------------------------------------------
C ONLY ALLOW PARCEL TO BE SATURATED IN INITIAL BGMKP1 = .TRUE.
C (STORED IN BRECAL AT THIS POINT)
C----------------------------------------------------------------------
C
       IF ( BGMK(I) .OR.( (XSQKP1(I) .GT. 0.) .AND. BRECAL(I) ) ) THEN
         BGMKP1(I) = .TRUE.
       ELSE
         BGMKP1(I) = .FALSE.
         XSQKP1(I) = 0.0
       END IF
C
       QPKP1(I) = QPKP1(I) - XSQKP1(I)
CL
CL----------------------------------------------------------------------
CL  RECALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE AT POINTS
CL  WHERE THE ENSEMBLE HAS BECOME UNSATURATED.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (28)
CL----------------------------------------------------------------------
CL
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
   30 CONTINUE
C
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     *             BGMKP1,BRECAL)
CL
CL----------------------------------------------------------------------
CL  BECAUSE OF THE REMOVAL OF LATENT HEATING, THE NEW PARCEL POTENTIAL
CL  TEMPERATURE MAY BE LOWER THAN ITS VALUE BEFORE THE DETRAINMENT
CL  CALCULATION. IN THIS CASE ABANDON THE DETRAINMENT CALCULATION.
CL----------------------------------------------------------------------
CL
      DO 90 I=1,NPNTS
       BDETK(I) = THPKP1(I) .GT. THPKP1W(I)
   90 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE POTENTIAL TEMPERATURE AND MIXING RATIO  OF DETRAINING
CL  AIR AND THE EXCESS WATER VAPOUR CONDESED FROM DETRAINING AIR
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (26)
CL----------------------------------------------------------------------
CL
      CALL THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,DQSK,
     *             BWKP1,EXK,PK)
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE DETRAINMENT RATE, DELTAK.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (31)
CL----------------------------------------------------------------------
CL
      CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     *             XSQKP1,THPKP1,BWKP1,BDETK,EKP14,EKP34,EXK,EXKP1)
C
      NREDO = 0
CL
CL----------------------------------------------------------------------
CL  ADD WATER VAPOUR WHICH WAS REMOVED FROM DETRAINING AIR INTO XSQKP1
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION 86), EQUATION (11C)
CL----------------------------------------------------------------------
CL
      DO 120 I=1,NPNTS
C
       EPSS = (1.+EKP14(I))*(1.+EKP34(I))
C
       IF (BDETK(I))
     * XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/
     *               (EPSS*(1.-DELTAK(I))))
CL
CL----------------------------------------------------------------------
CL  IF THE EXCESS WATER VAPOUR IN LAYER K+1 IS LESS THAN ZERO
CL  I.E. THE PARCEL HAS BECOME UNSATURATED THROUGH THE FORCED
CL  DETRAINMENT PROCESS THEN ABANDON THE CALCULATION
CL----------------------------------------------------------------------
CL
       BRECAL(I) = BGMKP1(I)
C
       BGMKP1(I) = XSQKP1(I) .GT. 0.
C
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
C
       IF (BRECAL(I)) THEN
          QPKP1(I)  = QPKP1(I) + XSQKP1(I)
     *               - (DELTAK(I)*XSQR(I)/(EPSS*(1.-DELTAK(I))))
          XSQKP1(I) = 0.
       ENDIF
C
C----------------------------------------------------------------------
C COUNT POINTS AT WHICH DETRAINMENT CALCULATION NEEDS REPEATING
C----------------------------------------------------------------------
C
       IF (BRECAL(I)) NREDO = NREDO + 1
  120 CONTINUE
CL
CL---------------------------------------------------------------------
CL  REPEAT CALCULATION OF PARCEL POTENTIAL TEMPERATURE, DETRAINMENT
CL  RATE AND EXCESS PARCEL WATER IF THE PARCEL BECOMES UNSATURATED
CL  IN LAYER K+1 AFTER FORCED DETARINMENT
CL---------------------------------------------------------------------
CL
      IF (NREDO .GT. 0) THEN
C
C----------------------------------------------------------------------
C  CALCULATE NEW PARCEL POTENTIAL TEMPERATURE IN LAYER K+1
C  AFTER FORCED DETRAINMENT
C----------------------------------------------------------------------
C
        CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     *               BGMKP1,BRECAL)
C
C----------------------------------------------------------------------
C  CHECK IF FORCED DETRAINMENT STILL POSSIBLE AND RESET RECALCUATION
C  MASK TO FALSE IF IT IS NOT
C----------------------------------------------------------------------
C
        DO 130 I=1,NPNTS
          IF (BRECAL(I)) THEN
            BDETK(I) = THPKP1(I) .GT. THPKP1W(I)
            BRECAL(I) = BDETK(I)
          END IF
  130   CONTINUE
C
C----------------------------------------------------------------------
C  RCALCULATE FORCED DETRAINEMNT RATE
C----------------------------------------------------------------------
C
        CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     *             XSQKP1,THPKP1,BWKP1,BRECAL,EKP14,EKP34,EXK,EXKP1)
C
C----------------------------------------------------------------------
C  RECALCULATE EXCESS WATER VAPOUR IN LAYER K+1
C  AFTER FORCED DETRAINMENT
C----------------------------------------------------------------------
C
        DO 140 I=1,NPNTS
         IF (BRECAL(I)) THEN
            EPSS = (1.+EKP14(I))*(1.+EKP34(I))
            XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/
     *                       (EPSS*(1.-DELTAK(I))))
         END IF
  140   CONTINUE
C
      END IF
CL
CL----------------------------------------------------------------------
CL  MAKE SURE THAT THE DETRAINMENT RATE IS BETWEEN 0 AND 1
CL
CL  IF <0 THEN NO DETRAINMENT OCCURS AND ORIGINAL VALUES ARE
CL  RESTORED
CL
CL  IF >1 THEN SET TO 1 AND THRK = THPK, QRK = QPK AND VALUES
CL  IN LAYER K+1 ARE RESTORED.  ALTHOUGH THESE ARE NOT USED
CL  IN ANY THERMODYNAMIC CALCULATION THEY ARE USED TO SPECIFY
CL CLOUD TOP IN SUBROUTIBE CONRAD
CL----------------------------------------------------------------------
CL
      DO 180 I=1,NPNTS
C
       IF (BDETK(I)) THEN
C
        IF (DELTAK(I).LE.0.0) THEN
           BDETK(I) = .FALSE.
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
           DELTAK(I) = 0.0
        ELSE IF (DELTAK(I).GT.1.0) THEN
           DELTAK(I) = 1.0
           THRK(I) = THPK(I)
           QRK(I) = QPK(I)
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
        END IF
C
       ENDIF
  180  CONTINUE
C
      RETURN
      END
