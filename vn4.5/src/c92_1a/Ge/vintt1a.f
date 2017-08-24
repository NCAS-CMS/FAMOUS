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
CLL  SUBROUTINE V_INT_T-----------------------------------------------
CLL
CLL  Purpose:  To calculate the temperature along an
CLL            arbitrary pressure level. Assumes input data is
CLL            on model levels.
CLL
CLL A.Dickinson <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 4.2     01/07/96 
CLL              : Revised for CRAY T3E. IC introduced to force exit
CLL              : from loop over levels when all points have processed.
CLL              : New arguments START and END introduced to
CLL              : facilitate the removal of duplicate calculations
CLL              : when using domain decomposition in MPP mode.     
CLL              : Author: A. Dickinson    Reviewer: F. Rawlins     
!LL 4.5     15/07/98
!LL                Use assumption that neighbouring points are
!LL                likely to be on or near same level. Jump out
!LL                of loop-over-levels once level found. Results
!LL                in a 40 percent speedup on 19 levels for
!LL                non-vector machines. S.D.Mullerworth
CLL  4.5    09/01/98  CRAY T3E optimisation: replace rtor_v by powr_v
CLL                                                    Deborah Salmond
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE V_INT_T
     &  (T,P,PL,PSTAR,P_EXNER_HALF,THETA,POINTS,P_LEVELS,L,AKH,BKH
     &  ,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points per level
     *,P_LEVELS  !IN Number of model levels
     *,L         !IN Reference level for below-surface T extrapolation
                 ! Use L=2
     *,START     !IN First point to be processed in POINTS dimension
     *,END       !IN Last point to be processed in POINTS dimension

      REAL
     * T(POINTS)     !OUT Temperature along input pressure surface
     *,P(POINTS)     !IN Pressure surface on which results required
     *,PL(POINTS)    !IN Pressure at reference level L
     *,PSTAR(POINTS) !IN Surface pressure
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,AKH(P_LEVELS+1) !IN Hybrid coords. A values at half levels.
     *,BKH(P_LEVELS+1) !IN Hybrid coords. B values at half levels.

C Workspace usage:-----------------------------------------------------
C
      REAL P_EXNER(POINTS)
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      REAL PTOP,PBOT,PKP2,PKP1,PK,PKM1,PPP1,PP,PPM1,P1,P2,P3
      INTEGER I,K
     &,LAST      ! Used to store level number of preceding point

      REAL P_EXNER_FULL_1,P_EXNER_FULL_2,TK,TERM1,TERM2,
     * P_EXNER_FULL_K,P_EXNER_FULL_KP1,P_EXNER_FULL_KM1
     *,P_EXNER_FULL_L,P_EXNER_FULL_LM1
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

      REAL CP_OVER_G,LAPSE_R_OVER_G
      PARAMETER(CP_OVER_G=CP/G)
      PARAMETER(LAPSE_R_OVER_G=LAPSE*R/G)

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


CL 2. Special cases   (i) Below ground
CL                   (ii) Bottom of bottom layer
CL                  (iii) Top of top layer and above



      LAST=2 ! Arbitrary initialisation
      DO I=START,END

C Convert target pressure to Exner pressure

        P_EXNER(I)=(P(I)/PREF)**KAPPA

! Start from same level as last point. First check whether this point
! is above or below, then continue search in appropriate direction
        IF(P_EXNER(I).GT.P_EXNER_HALF(I,LAST))THEN

! These next two loops exit immediately once level found.
! GOTO cuts out needless looping once level is found, reducing the
! cost of the routine by about 40 percent for 19 level runs.
          DO K=LAST,3,-1
            IF(P_EXNER(I).LE.P_EXNER_HALF(I,K-1))THEN
              GOTO 240
            ENDIF
          ENDDO
        ELSE
          DO K=LAST+1,P_LEVELS
            IF(P_EXNER(I).GT.P_EXNER_HALF(I,K))THEN
              GOTO 240
            ENDIF
          ENDDO
        ENDIF
 240    CONTINUE

! At this point, K is:
!    2           for below bottom level.
!    P_LEVELS+1  for above top level
!    Otherwise K is the level just above the point


        IF(K.EQ.2)THEN
CL (i) Below ground: equation (3.15)
CL A lapse rate of 6.5 deg/km is assumed. L is a
CL reference level - usually the first level above the model
CL boundary layer.

          IF(P(I).GT.PSTAR(I))THEN

            PTOP = AKH(L+1) + BKH(L+1)*PSTAR(I)
            PBOT = AKH(L)   + BKH(L)  *PSTAR(I)
            P_EXNER_FULL_L = P_EXNER_C
     +        (P_EXNER_HALF(I,L+1),P_EXNER_HALF(I,L),PTOP,PBOT,KAPPA)
            T(I)=THETA(I,L)*P_EXNER_FULL_L
     *        *(P(I)/PL(I))**LAPSE_R_OVER_G

CL (ii) Bottom layer: equation (3.14)

          ELSE

            P1 = AKH(1) + BKH(1)*PSTAR(I)
            P2 = AKH(2) + BKH(2)*PSTAR(I)
            P3 = AKH(3) + BKH(3)*PSTAR(I)
            P_EXNER_FULL_1 = P_EXNER_C
     +        (P_EXNER_HALF(I,2),P_EXNER_HALF(I,1),P2,P1,KAPPA)
            P_EXNER_FULL_2 = P_EXNER_C
     +        (P_EXNER_HALF(I,3),P_EXNER_HALF(I,2),P3,P2,KAPPA)

            TK=THETA(I,1)*P_EXNER_FULL_1

            TERM1=(TK-THETA(I,2)*P_EXNER_FULL_2)
     *        *THETA(I,1)*(P_EXNER(I)-P_EXNER_FULL_1)

            TERM2=THETA(I,2)*(P_EXNER_HALF(I,2)-P_EXNER_FULL_2)
     *        +THETA(I,1)*(P_EXNER_FULL_1-P_EXNER_HALF(I,2))

            T(I)=TK+TERM1/TERM2

          ENDIF

          LAST=2                ! Next point, start from level 2
CL (iii) Top layer and above: equation (3.13)
        ELSEIF(K.EQ.P_LEVELS+1)THEN
          PPP1 = AKH(P_LEVELS+1) + BKH(P_LEVELS+1)*PSTAR(I)
          PP   = AKH(P_LEVELS  ) + BKH(P_LEVELS  )*PSTAR(I)
          PPM1 = AKH(P_LEVELS-1) + BKH(P_LEVELS-1)*PSTAR(I)

          P_EXNER_FULL_L   = P_EXNER_C (P_EXNER_HALF(I,P_LEVELS+1),
     +      P_EXNER_HALF(I,P_LEVELS),PPP1,PP,KAPPA)
          P_EXNER_FULL_LM1 = P_EXNER_C (P_EXNER_HALF(I,P_LEVELS),
     +      P_EXNER_HALF(I,P_LEVELS-1),PP,PPM1,KAPPA)

          TK=THETA(I,P_LEVELS)*P_EXNER_FULL_L

          TERM1=(TK-THETA(I,P_LEVELS-1)*P_EXNER_FULL_LM1)
     *      *THETA(I,P_LEVELS)*(P_EXNER_FULL_L-P_EXNER(I))

          TERM2=THETA(I,P_LEVELS)*(P_EXNER_HALF(I,P_LEVELS)
     &      -P_EXNER_FULL_L)+THETA(I,P_LEVELS-1)
     &      *(P_EXNER_FULL_LM1-P_EXNER_HALF(I,P_LEVELS))

          T(I)=TK+TERM1/TERM2
          LAST=P_LEVELS         ! Next point, start from top
        ELSE
CL 3. Middle levels: equation (3.12)
CL Two alternatives are used depending on whether P_EXNER(I) falls in
CL the top or bottom half of layer k.

          PKP1 = AKH(K) + BKH(K)*PSTAR(I)
          PK   = AKH(K-1)   + BKH(K-1)*PSTAR(I)

          P_EXNER_FULL_K   = P_EXNER_C (P_EXNER_HALF(I,K),
     +      P_EXNER_HALF(I,K-1),PKP1,PK,KAPPA)

C Top half of layer k.

          IF(P_EXNER(I).LE.P_EXNER_FULL_K)THEN

            PKP2 = AKH(K+1) + BKH(K+1)*PSTAR(I)
            P_EXNER_FULL_KP1 = P_EXNER_C (P_EXNER_HALF(I,K+1),
     +        P_EXNER_HALF(I,K),PKP2,PKP1,KAPPA)

            TK=THETA(I,K-1)*P_EXNER_FULL_K

            TERM1=(THETA(I,K)*P_EXNER_FULL_KP1-TK)
     *        *THETA(I,K-1)*(P_EXNER_FULL_K-P_EXNER(I))

            TERM2=THETA(I,K-1)*(P_EXNER_FULL_K-P_EXNER_HALF(I,K))
     *        +THETA(I,K)*(P_EXNER_HALF(I,K)-P_EXNER_FULL_KP1)

            T(I)=TK+TERM1/TERM2

          ENDIF

C Bottom half of layer k.

          IF(P_EXNER(I).GT.P_EXNER_FULL_K)THEN

            PKM1 = AKH(K-2) + BKH(K-2)*PSTAR(I)
            P_EXNER_FULL_KM1 = P_EXNER_C (P_EXNER_HALF(I,K-1),
     +        P_EXNER_HALF(I,K-2),PK,PKM1,KAPPA)

            TK=THETA(I,K-1)*P_EXNER_FULL_K

            TERM1=(THETA(I,K-2)*P_EXNER_FULL_KM1-TK)
     *        *THETA(I,K-1)*(P_EXNER(I)-P_EXNER_FULL_K)

            TERM2=THETA(I,K-1)*(P_EXNER_HALF(I,K-1)-P_EXNER_FULL_K)
     *        +THETA(I,K-2)*(P_EXNER_FULL_KM1-P_EXNER_HALF(I,K-1))

            T(I)=TK+TERM1/TERM2

          ENDIF
          LAST=K                ! Next point, start from level K

        ENDIF                   ! IF(K.EQ.2)...ELSEIF...ELSE

      ENDDO                     ! DO I=START,END

      RETURN
      END

