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
CLL  SUBROUTINE SATCAL-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES SATURATED TEMPERATURE
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
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE SATCAL (NPNTS,T,TH,PK,QS,THDDS,K,EXK,Q_K,THE_K)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER IC                ! LOOP COUNTER
C
      INTEGER NPNTS             ! VECTOR LENGTH
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL TH(NPNTS)            ! IN POTENTIAL TEMPERATURE (K)
C
      REAL T(NPNTS)             ! IN TEMPERATURE (K)
C
      REAL PK(NPNTS)            ! IN PRESSURE OF LAYER K (PA)
C
      REAL Q_K(NPNTS)           ! IN MIXING RATIO OF LAYER K (KG/KG)
C
      REAL EXK(NPNTS)           ! IN EXNER RATIO OF LAYER K
C
      REAL THE_K(NPNTS)         ! IN ENVIRONMENTAL POTENTIAL TEMPERATURE
                                !    IN LAYER K
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL QS(NPNTS)            ! OUT SATURATED SPECIFIC HUMIDITY
                                !     (KG/KG)
C
      REAL THDDS(NPNTS)         ! OUT SATURATED ENVIRONMENTAL
                                !     POTENTIAL TEMPERATURE (K)
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
      REAL TS(NPNTS)            ! SATURATED TEMPERATURE (K)
C
      REAL T_FG(NPNTS)          ! TEMPERATURE FIRST GUESS (K)
C
      REAL TH_FG(NPNTS)         ! POTENTIAL TEMPERATURE FIRST GUESS (K)
C
      REAL DQBYDT(NPNTS)        ! FIRST GUESS AT MIXING RATIO INCREMENT
                                ! (KG/KG/K)
C
      REAL L                    ! LATENT HEAT
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL QSAT, DQS_DTH
C
C-----------------------------------------------------------------------
C SET INITIAL FIRST GUESS TEMPERATURE AND THETA - BASED UPON
C ENVIRONMENTAL TEMPERATURE IN LAYER K
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TH_FG(I) = THE_K(I)
       T_FG(I) = TH_FG(I)*EXK(I)
      END DO
C
C----------------------------------------------------------------------
C CALCULATE QSAT FOR INITIAL FIRST GUESS TEMPERATURE
C----------------------------------------------------------------------
C
      CALL QSAT(QS,T_FG,PK,NPNTS)
C
C----------------------------------------------------------------------
C DO TWO ITERATIONS TO FIND SATURATION POINT DUE TO EVAPORATION
C----------------------------------------------------------------------
C
      DO IC=1,2
C
C----------------------------------------------------------------------
C CALCULATE DQSAT/DT FOR FIRST GUESS TEMPERATURE
C----------------------------------------------------------------------
C
       CALL DQS_DTH(DQBYDT,K,TH_FG,QS,EXK,NPNTS)
C
C----------------------------------------------------------------------
C CALCULATE UPDATED TEMPERATURE AT SATURATION
C----------------------------------------------------------------------
C
       DO I=1,NPNTS
C
        IF (T_FG(I).GT.TM) THEN
         L=LC
        ELSE
         L=LC+LF
        END IF
C
        THDDS(I) = (TH(I) - (L/(CP*EXK(I)))*(QS(I)-Q_K(I)-
     *                  TH_FG(I)*DQBYDT(I))) /
     *                  (1+(L/(CP*EXK(I)))*DQBYDT(I))
C
C----------------------------------------------------------------------
C CALCULATE TEMPERATURE AT SATURATION AND UPDATE FIRST GUESS
C----------------------------------------------------------------------
C
        TH_FG(I) = THDDS(I)
        T_FG(I) = TH_FG(I)*EXK(I)
C
       END DO
C
C----------------------------------------------------------------------
C CALCULATE REVISED SATURATION MIXING RATIO AT SATURATION
C---------------------------------------------------------------------
C
       CALL QSAT(QS,T_FG,PK,NPNTS)
C
      END DO
C
      RETURN
      END
C
