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
CLL  SUBROUTINE V_INT_ZH----------------------------------------------
CLL
CLL  Purpose:  Calculates height of each layer boundary (half level)
CLL            using the hydrostatic approximation.
CLL
CLL A.Dickinson <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   4.2  16/07/96  THETAV made into scalar temporary. Saves
CLL                   dynamic memory allocation and memory accesses.
CLL
CLL                  Author: A. Dickinson     Reviewer: R. Rawlins
CLL
CLL Programming standard :
CLL
CLL Logical components covered : D473
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE V_INT_ZH
     *(P_EXNER_HALF,THETA,Q,PHI_STAR,ZH,POINTS,P_LEVELS,Q_LEVELS)
C     ,START,END)  ! These arguments not yet implemented

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points to be processed
     *,P_LEVELS  !IN Number of model levels
     *,Q_LEVELS  !IN Number of wet levels
     *,START     !IN First point to be processed in POINTS dimension
     *,END       !IN Last point to be processed in POINTS dimension


      REAL
     * P_EXNER_HALF(POINTS,P_LEVELS+1)!IN Exner pressure at model
     *                                !    half levels
     *,THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,Q(POINTS,Q_LEVELS)     !IN Specific humidity at full levels
     *,PHI_STAR(POINTS)       !IN Geopotential height of topography
     *,ZH(POINTS,P_LEVELS+1)  !OUT Height of model half levels

C Workspace usage:-----------------------------------------------------
C None
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,K
      REAL DEL_EXNER,THETAV
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

C----------------------------------------------------------------------

CL 1. Define local constants and initialise height of first layer
CL    boundary.

      REAL ONE_OVER_G,CP_OVER_G
      PARAMETER(CP_OVER_G=CP/G)
      PARAMETER(ONE_OVER_G=1/G)

C Height of orography in metres
C  Interim fix for START,END arguments
      START=1
      END  =POINTS

      DO 100 I=1,POINTS
      ZH(I,1)=ONE_OVER_G*PHI_STAR(I)
100   CONTINUE

CL 2. Compute height of remaining layer boundaries

C Loop over wet levels
      DO 200 K=1,Q_LEVELS
       DO 210 I=START,END

C Calculate virtual potential temperature
        THETAV=THETA(I,K)*(1.0+C_VIRTUAL*Q(I,K))

C Accumulate heights using equation (3.5)
        DEL_EXNER=P_EXNER_HALF(I,K)-P_EXNER_HALF(I,K+1)
        ZH(I,K+1)=ZH(I,K)+CP_OVER_G*THETAV*DEL_EXNER

210    CONTINUE
200   CONTINUE
C Loop over dry levels
      DO 230 K=Q_LEVELS+1,P_LEVELS
       DO 220 I=START,END

C Accumulate heights using equation (3.5)
        DEL_EXNER=P_EXNER_HALF(I,K)-P_EXNER_HALF(I,K+1)
        ZH(I,K+1)=ZH(I,K)+CP_OVER_G*THETA(I,K)*DEL_EXNER

220    CONTINUE
230   CONTINUE


      RETURN
      END


