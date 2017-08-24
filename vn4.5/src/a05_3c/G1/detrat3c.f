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
CLL  SUBROUTINE DET_RATE-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE FORCED DETRAINMENT RATE IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
CLL
CLL  PROGRAMMING STANDARDS :
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL                  SECTION (6), EQUATION (31)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     *                   XSQKP1,THPKP1,BWKP1,BCALC,EKP14,EKP34,
     *                   EXK,EXKP1)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS            ! VECTOR LENGTH
C
      INTEGER I                ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THRK(NPNTS)         ! IN PARCEL DETRAINMENT POTENTIAL
                               !    TEMPERATURE IN LAYER K (K)
C
      REAL XSQR(NPNTS)         ! IN EXCESS WATER VAPOUR OF THE
                               !    DETRAINING AIR IN LAYER K (KG/KG)
C
      REAL THPK(NPNTS)         ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
C
      REAL THEK(NPNTS)         ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
C
      REAL THEKP1(NPNTS)       ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
C
      REAL XSQKP1(NPNTS)       ! IN EXCESS WATER VAPOUR OF THE PARCEL
                               !    IN LAYER K+1 (KG/KG)
C
      REAL THPKP1(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
C
      LOGICAL BCALC(NPNTS)     ! IN MASK FOR POINTS AT WHICH
                               !    CALCULATIONS OF THIS ROUTINE
                               !    ARE NEEDED
C
      LOGICAL BWKP1(NPNTS)     ! IN MASK FOR THOSE POINTS AT WHICH
                               !    CONDENSATE IS LIQUID IN LAYER K+1
C
      REAL EKP14(NPNTS)        ! IN ENTRAINEMNT RATE FOR LEVEL K+1/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
C
      REAL EKP34(NPNTS)        ! IN ENTRAINEMNT RATE FOR LEVEL K+3/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
C
      REAL EXK(NPNTS)          ! IN EXNER RATIO FOR LEVEL K
C
      REAL EXKP1(NPNTS)        ! IN EXNER RATIO FOR LEVEL K+1
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL DELTAK(NPNTS)       ! OUT PARCEL FORCED DETRAINMENT RATE
                               !     IN LAYER K MULTIPLIED BY
                               !     APPROPRIATE LAYER THICKNESS
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                  ! LATENT HEAT OF CONDENSATION OR
                               ! (CONDENSATION + FUSION) (J/KG)
C
      REAL EPSS                ! (1+EKP14)*(1+EKP34)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
C
      DO 10 I=1,NPNTS
       IF (BCALC(I)) THEN
C
C-----------------------------------------------------------------------
C   CREATE A VECTOR OF LATENT HEATS
C-----------------------------------------------------------------------
C
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
C
C-----------------------------------------------------------------------
C   CALCULATE DETRAINMENT RATES
C-----------------------------------------------------------------------
C
       EPSS = (1. + EKP14(I)) * (1. + EKP34(I))
          DELTAK(I) = EKP14(I)*THEK(I)
     *        + EKP34(I)*(1.+EKP14(I))*THEKP1(I)
     *        - EPSS*(THPKP1(I) - EL/(EXKP1(I)*CP) * XSQKP1(I))
C
          DELTAK(I) =   (DELTAK(I) + THPK(I))
     *  *(EXK(I)*CP)/((DELTAK(I) + THRK(I))*(EXK(I)*CP) - EL*XSQR(I))
C
C----------------------------------------------------------------------
C  FROM A THEORETICAL VIEW POINT DELTAK CANNOT = 1 . HOWEVER
C  BECAUSE OF APPROXIMATION USED IN THE CALCULATION NUMERICALLY IT
C  MAY BE POSSIBLE.  HENCE IF DELTAK = 1 SET IT TO SLIGHTLY SMALLER
C  THAN 1
C----------------------------------------------------------------------
C
          DELTAK(I) = MIN(0.99999,DELTAK(I))
C
       ENDIF
   10 CONTINUE
C
      RETURN
      END
