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
CLL  SUBROUTINE PMSL--------------------------------------------------
CLL
CLL  Purpose:  Calculates mean sea level pressure
CLL
CLL AD, DR      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  4.2    23/07/96  Revised for CRAY T3E. Option to use fast vector 
CLL                   like function for raising to a power introduced
CLL                   under DEF T3E.
CLL                   New arguments START and END introduced to
CLL                   facilitate the removal of duplicate calculations
CLL                   when using domain decomposition.
CLL                   Author A. Dickinson     Reviewer: F. Rawlins
!LL  4.5    20/4/98   Implement the START,END described above.
!LL                   S.D.Mullerworth
CLL
CLL  4.5    09/01/98  CRAY T3E optimisation: replace rtor_v by powr_v
CLL                                                    Deborah Salmond
CLL Programming standard :
CLL
CLL Logical components covered : D441
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE PMSL
     &  (P_MSL,PL,PSTAR,P_EXNER_HALF,THETA,Q,PHI_STAR
     &  ,POINTS,P_LEVELS,Q_LEVELS,L,AKH,BKH
     &  ,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points per level
     *,P_LEVELS  !IN Number of model levels
     *,Q_LEVELS  !IN Number of wet levels
     *,L         !IN Reference level for below surface T extrapolation.
     *           ! = No of B.L. levels plus one
     *,START     !IN First point to be processed
     *,END       !IN Last point to be processed

      REAL
     * P_MSL(POINTS)    !OUT Mean sea level pressure
     *,PHI_STAR(POINTS) !IN Geopotential height of topography
     *,PL(POINTS)       !IN Reference pressure at level L
     *,PSTAR(POINTS)    !IN Surface pressure
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,Q(POINTS,Q_LEVELS)     !IN Specific humidity at full levels
     *,AKH(P_LEVELS+1)        !IN Hybrid Coords. A and B values
     *,BKH(P_LEVELS+1)        !IN at half levels.

C Workspace usage:-----------------------------------------------------
      REAL TEMP(POINTS),POWER
C ---------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,K
      REAL PTOP    ! Pressure at top of layer
      REAL PBOT    ! Pressure at bottom of layer
      REAL P_EXNER_FULL  ! Exner Pressure at full model level
      REAL TS            ! Surface Temperature
      REAL ALOGHF,EXPHF
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

C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------
C----------------------------------------------------------------------

CL 1. Define local constants

      REAL LAPSE_R_OVER_G,G_OVER_LAPSE_R,ONE_OVER_G
      PARAMETER(LAPSE_R_OVER_G=LAPSE*R/G)
      PARAMETER(G_OVER_LAPSE_R=1./LAPSE_R_OVER_G)
      PARAMETER(ONE_OVER_G=1./G)

C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------


CL 2. Calculate mean sea level pressure: equations (3.8) & (3.11)


      DO I=START,END
      TEMP(I)=(PSTAR(I)/PL(I))**LAPSE_R_OVER_G

C Estimate surface temperature: equation (3.8)
      PTOP = AKH(L+1) + BKH(L+1)*PSTAR(I)
      PBOT = AKH(L)   + BKH(L)  *PSTAR(I)
      P_EXNER_FULL = P_EXNER_C
     +(P_EXNER_HALF(I,L+1),P_EXNER_HALF(I,L),PTOP,PBOT,KAPPA)
      TS=THETA(I,L)*P_EXNER_FULL*TEMP(I)
      TS=TS*(1.0+C_VIRTUAL*Q(I,1))
      TEMP(I)=(TS+LAPSE*ONE_OVER_G*PHI_STAR(I))/TS
      ENDDO

C Calculate PMSL using equation (3.11)


      DO I=START,END
      TEMP(I)=TEMP(I)**G_OVER_LAPSE_R

      P_MSL(I)=PSTAR(I) * TEMP(I)

      ENDDO

      RETURN
      END

