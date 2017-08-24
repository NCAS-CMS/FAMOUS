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
CLL  SUBROUTINE ADD_ENG_CORR--------------------------------------
CLL
CLL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
CLL            - TO ADD IN TEMPERATURE CORRECTION TO
CLL              GLOBAL TEMPERATURE FIELD SO TO
CLL              CONSERVE TOTAL ENERGY GLOBALLY
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL LOGICAL COMPONENTS COVERED:
CLL
CLL  PROJECT TASK :
CLL
CLL  DOCUMENTATION :
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE ADD_ENG_CORR (ENERGY_CORR,T,P_FIELD,NPNTS,P_LEVELS,
     &                         TSTEP,TOT_MASS,TOT_FLUXES)
C
      IMPLICIT NONE
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
C VECTOR LENGTHS
C----------------------------------------------------------------------
C
      INTEGER P_FIELD          ! IN VECTOR LENGTH OF P GRID
C
      INTEGER NPNTS            ! IN VECTOR LENGTH FOR CALCULATION
C
      INTEGER P_LEVELS         ! IN NUMBER OF LEVELS IN VERTICAL
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL ENERGY_CORR         ! IN ENERGY CORRECTION
C
      REAL TSTEP               ! IN TIMESTEP
C
      REAL TOT_MASS            ! IN MASS OF ATMOSPHERE
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE IN AND OUT
C----------------------------------------------------------------------
C
      REAL T(P_FIELD,P_LEVELS) ! INOUT TEMPERATURE
C
      REAL TOT_FLUXES          !INOUT SUM OF DIABATIC FLUXES
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY  -  NONE
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C INTERNAL LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER I                ! LOOP COUNTER
C
      INTEGER K                ! LOOP COUNTER
C
C----------------------------------------------------------------------
C EXTERNAL SUBROUTINE CALLS  -  NONE
C----------------------------------------------------------------------
C
C*---------------------------------------------------------------------
C
C----------------------------------------------------------------------
C CORRECT TEMPERATURE FOR ERROR IN ENERGY BUDGET OF THE
C PREVIOUS DAY
C----------------------------------------------------------------------
C
      DO K=1,P_LEVELS
       DO I=1,NPNTS
        T(I,K) = T(I,K) + ENERGY_CORR*TSTEP
       END DO
      END DO
C
C
C----------------------------------------------------------------------
C ADD ENERGY CORRECTION INTO SUM OF DIABATIC FLUXES
C----------------------------------------------------------------------
C
      TOT_FLUXES = TOT_FLUXES + CP*ENERGY_CORR*TOT_MASS*TSTEP
C
      RETURN
      END
