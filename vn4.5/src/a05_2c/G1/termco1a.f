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
CLL  SUBROUTINE TERM_CON-----------------------------------------------
CLL
CLL  PURPOSE : RETURENS A MASK FOR POINTS AT WHICH CONVECTION
CLL            IS TERMINATING
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
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE TERM_CON(NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1,
     *                    THPI,QPI,QSEKP1,DELTAK,EXKP1,EKP14,EKP34,
     *                    PSTAR)
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

      REAL MPARFL  !  MINIMUM PARCEL MASS FLUX
                   !  = 1E-3 * MINIMUM PARCEL BUOYANCY *
                   !              MASS FLUX PARAMETER C
      PARAMETER (MPARFL = 1.0E-3 * 1.0 * 3.33E-4)
C
      REAL XSBMIN !  MINIMUM EXCESS BUOYANCY TO CONTINUE PARCEL ASCENT
                  !  (K)
      PARAMETER (XSBMIN = 0.2)
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      REAL QSTICE   !  APPROXIMATION TO SATURATION MIXING RATIO
                    !  AT TEMPERATURE AT WHICH LIQUID WATER TURNS TO
                    !  ICE (SEE COMDECK TICE) (KG/KG)
      PARAMETER (QSTICE = 3.5E-3)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYER
C
      INTEGER K              ! IN PRESENT MODEL LAYER
C
      INTEGER I              ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL THPI(NPNTS)       ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
C
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
C
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      REAL DELTAK(NPNTS)     ! IN FORCED DETRAINMENT IN LAYER K
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO FOR LEVEL K+1
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BTERM(NPNTS)   ! OUT MASK OF THOSE POINTS AT WHICH
                             !     CONVECTION IS ENDING
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
C
      REAL FLXMIN            ! MINIMUM CONVECTIVE MASSFLUX BELOW
                             ! WHICH TERMINAL DETRAINMENT OCCURS
                             ! (PA/S)
C
      REAL THVUNDI           ! POTENTIAL TEMPERATURE OF AN
                             ! UNDILUTE PARCEL IN LAYER K+1
                             ! FROM THE STARTING LAYER OF
                             ! CONVECTION (K)
C
      REAL THVEKP1           ! VIRTUAL POTENTIAL TEMPERATURE
                             ! OF ENVIRONMENT IN LAYER K+1 (K)
C
C*---------------------------------------------------------------------
C
C----------------------------------------------------------------------
C  CALCULATE MINIMUM MASS FLUX BELOW WHICH CONVECTION IS TERMINATED
C----------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
        FLXMIN = MPARFL*(1.+EKP14(I))*(1.+EKP34(I))*PSTAR(I)
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
CL
CL----------------------------------------------------------------------
CL  PARCELS ARE ONLY CONSIDERED FOR TERMINATION IF THEY ARE DETRAINING
CL  EXCEPT AT THE TOP MODEL LAYER, WHERE ALL CONVECTION TERMINATES
CL
CL  IF THE PARCEL HAS A POTENTIAL TEMPETURE GREATER THAN THE
CL  POTENTIAL TEMPERATURE OF AN UNDILUTE PARCEL FORM THE STARTING
CL  LAYER OF CONVECION IN LAYER K+1 THEN CONVECTION IS TERMINATED
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (7), EQUATION (32)
CL
CL  CONVECTION IS ALSO TERMINATED IF MASS FLUX IN LAYER K+1 IS LESS
CL  IS LESS THAN A MINIMUM VALUE
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (7), EQUATION (33)
CL----------------------------------------------------------------------
CL
       THVUNDI=( THPI(I) + (EL/(EXKP1(I)*CP)) *(QPI(I) - QSEKP1(I))
     *         +((LC-EL)/(EXKP1(I)*CP))*MAX(0.0,(QPI(I)-QSTICE))
     *         )*(1.+C_VIRTUAL*QSEKP1(I))
C
       THVEKP1 = (THEKP1(I)*(1.+C_VIRTUAL*QEKP1(I)) + XSBMIN)
C
       BTERM(I) = (((FLXKP1(I) .LT. FLXMIN) .OR. (THVUNDI .LT. THVEKP1))
     *               .AND. (DELTAK(I).GT.0.0)) .OR. (K+1) .EQ. NLEV
C
  10  CONTINUE
C
      RETURN
      END
