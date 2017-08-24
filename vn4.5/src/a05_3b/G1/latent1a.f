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
CLL  SUBROUTINE LATENT_H-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES A NEW PARCEL TEMPERATURE AFTER
CLL            CONDENSATION HAS OCCURRED
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
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LATENT_H (NPNTS,THPKP1,QPKP1,THEKP1,
     *                     QSEKP1,DQSKP1,BGMKP1,BWKP1,EXKP1)
C
C-----------------------------------------------------------------------
C   ADJUSTS PARCEL POTENTIAL TEMPERATURES TO ACCOUNT FOR THE LATENT
C   HEAT RELEASE FROM SUPERSATURATED PARCELS CONDENSING/DEPOSITING OUT
C----------------------------------------------- AUTHOR: M FISHER 1987 -
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
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
      INTEGER NPNTS            ! VECTOR LENGTHS
C
      INTEGER I                ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL QPKP1(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K+1
                               !    (KG/KG/S)
C
      REAL THEKP1(NPNTS)       ! IN ENVIRONMENT POTENTIAL FOR LAYER K+1
                               !    (K/S)
C
      REAL QSEKP1(NPNTS)       ! IN SATURATION MIXING RATIO OF THE
                               !    ENVIRONMENT IN LAYER K+1 (KG/KG/S)
C
      REAL DQSKP1(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                               !    WITH POTENTIAL TEMPERATURE FOR THE
                               !    ENVIRONMENT IN LAYER K+1 (KG/KG/K)
C
      LOGICAL BGMKP1(NPNTS)    ! IN MASK FOR PARCELS WHICH ARE
                               !    SATURATED LAYER K+1
C
      LOGICAL BWKP1(NPNTS)     ! IN MASK FOR POINTS IN WHICH
                               ! CONDENSATE IS LIQUID IN LAYER K+1
C
      REAL EXKP1(NPNTS)        ! IN EXNER RATIO AT MID-POINT OF
                               !    LAYER K+1
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL THPKP1(NPNTS)       ! INOUT
                               ! IN  INITIAL ESTIMATE OF PARCEL
                               !     POTENTIAL TEMPERATURE IN
                               !     LAYER K+1 (K/S)
                               ! OUT PARCEL POTENTIAL TEMPERATURE
                               !     IN LAYER K+1 AFTER LATENT
                               !     HEATING (K/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EL                  ! LATENT HEAT OF CONDENSATION OR
                               ! (CONDENSATION + FUSION) (J/KG)
C
C*----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  ADJUST PARCEL POTENTIAL TEMPERATURES TO ACCOUNT FOR LATENT HEATING
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (4), EQUATION(21)
CL----------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
C
       IF (BGMKP1(I)) THEN
         THPKP1(I) = ( THPKP1(I) +
     *   (EL/(EXKP1(I)*CP))*(QPKP1(I)-QSEKP1(I)+THEKP1(I)*DQSKP1(I))
     *              ) / (1.+(EL/(EXKP1(I)*CP))*DQSKP1(I))
       ENDIF
   10  CONTINUE
C
      RETURN
      END
