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
CLL  SUBROUTINE V_INT_Z-----------------------------------------------
CLL
CLL  Purpose:  Calculates the height of an arbitrary pressure
CLL            surface. Since version 2, the top model level is
CLL            ignored in making the calculation.
CLL
CLL  28/04/92 calculation of P_EXNER_FULL consistent with
CLL     the geopotential eqn. New arguments AKH and BKH.
CLL
CLL RR, DR, AD  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  4.2    01/08/96  
CLL              : *DEF CRAY removed. alogh and exph 32-bit functions
CLL              : no longer required.
CLL              : New arguments START and END introduced to
CLL              : facilitate the removal of duplicate calculations
CLL              : when using domain decomposition in MPP mode.
CLL              : New variable IC introduced to remove superfluous  
CLL              : tests over levels.
CLL
CLL                   Author: A. Dickinson   Reviewer: R. Rawlins
!LL 4.5     15/07/98
!LL                Use assumption that neighbouring points are
!LL                likely to be on or near same level. Jump out
!LL                of loop-over-levels once level found. Results
!LL                in a 30 percent speedup on 19 levels for
!LL                non-vector machines. S.D.Mullerworth
CLL  4.5    09/01/98  CRAY T3E optimisation: replace rtor_v by powr_v
CLL                                                    Deborah Salmond
CLL
CLL Programming standard :
CLL
CLL Logical components covered : D471
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE V_INT_Z(P,PL,PSTAR,P_EXNER_HALF,THETA,Q,ZH,Z,POINTS
     *  ,P_LEVELS_MODEL,Q_LEVELS_MODEL,L,AKH,BKH
     &  ,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS          !IN Number of points to be processed
     *,P_LEVELS_MODEL  !IN Number of model levels
     *,Q_LEVELS_MODEL  !IN Number of model wet levels
     *,L         !IN Reference level for below surface T extrapolation.
     *           ! = No of B.L. levels plus one
     *,START     !IN First point to be processed in POINTS dimension
     *,END       !IN Last point to be processed in POINTS dimension


      REAL
     * P(POINTS)     !IN Pressure surface on which results required
     *,PL(POINTS)    !IN Pressure at reference level L
     *,PSTAR(POINTS) !IN Surface pressure
     *,P_EXNER_HALF(POINTS,P_LEVELS_MODEL+1) !IN Exner pressure at
     *                                       !   model half levels
     *,Z(POINTS)                   !OUT Height of pressure surface P
     *,ZH(POINTS,P_LEVELS_MODEL)   !IN Height of model half levels
     *,THETA(POINTS,P_LEVELS_MODEL)!IN Potential temp at full levels
     *,Q(POINTS,Q_LEVELS_MODEL)    !IN Specific humidity at full levels
     *,AKH(P_LEVELS_MODEL+1)       !IN Hybrid coord. A at half levels
     *,BKH(P_LEVELS_MODEL+1)       !IN Hybrid coord. B at half levels

C Local arrays:--------------------------------------------------------
      REAL P_EXNER(POINTS) ! Exner pressure at required pressure level
C ---------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,K    !DO loop indices
     *,K_BELOW       !K-1 or bottom level
     *,P_LEVELS      !No of model levels minus one
     *,Q_LEVELS      !No of wet levels (minus one if same as P_LEVELS)
     *,LAST          !Stores level of preceding point


      REAL
     * P_EXNER_FULL      ! Exner pressure on full level nearest P
     *,P_EXNER_FULL_UP   ! Exner pressure on full level above
     *,P_EXNER_FULL_LOW  ! Exner pressure on full level below
     *,DEL_EXNER         ! Vertical difference of exner pressure
     *,THETAV            ! Virtual potential temperature
     *,LOCAL_GRADIENT    ! Local potential temperature gradient
     *,T_GRADIENT        ! Temperature gradient
     *,SECOND_ORDER_TERM ! 2nd order correction to hydrostatic integral
     *,P_EXNER_FULL_L    ! Full level exner on level L
     *,TS                ! Surface temperature
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

      REAL LAPSE_R_OVER_G,CP_OVER_G
      PARAMETER(LAPSE_R_OVER_G=LAPSE*R/G)
      PARAMETER(CP_OVER_G=CP/G)

      REAL
     &    PUP1,PUP,PLOW,PLOW11,PLOW1

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


C SET TOP LEVEL FOR INTERPOLATION EQUAL TO TOP MODEL LAYER
C (CAUTION! THIS HAS BEEN OF DOUBTFUL QUALITY IN THE PAST FOR 20-LEV M)
      P_LEVELS=P_LEVELS_MODEL
      Q_LEVELS=MIN(P_LEVELS,Q_LEVELS_MODEL)

!L 2. Special cases   (i) Below ground
!L                   (ii) Top layer

      DO I=START,END

! Convert target pressure to Exner pressure
! ie  P_EXNER(I)=(P(I)/PREF)**KAPPA

        P_EXNER(I)=P(I)/PREF
        P_EXNER(I)=P_EXNER(I)**KAPPA

      ENDDO


      LAST=2 ! Arbitrary initialisation
      DO I=START,END
! Start from level of last point. First check whether this point
! is above or below, then continue search in appropriate direction
        IF(P_EXNER(I).LT.P_EXNER_HALF(I,LAST))THEN
          DO K=LAST,P_LEVELS-1
            IF(P_EXNER(I).GE.P_EXNER_HALF(I,K+1))THEN
              GOTO 210
            ENDIF
          ENDDO
        ELSE
          DO K=LAST-1,1,-1
            IF(P_EXNER(I).LT.P_EXNER_HALF(I,K))THEN
              GOTO 210
            ENDIF
          ENDDO
        ENDIF
 210    CONTINUE
! Here, K=0 for below level 1. K=P_LEVELS for above top level
! Otherwise K is set to level in range 1 to P_LEVELS-1 just below point.

        IF (K.EQ.0)THEN
!L (i) Below ground: equation (3.11).
!L A lapse rate of 6.5 deg/km is assumed. L is a
!L reference level - usually the first level above the model
!L boundary layer.
! Test for P>P* using Exner, for consistency with other sections
          IF(P_EXNER(I).GE.P_EXNER_HALF(I,1))THEN
            PUP=PSTAR(I)*BKH(L+1) + AKH(L+1)
            PLOW=PSTAR(I)*BKH(L) + AKH(L)
            P_EXNER_FULL_L= P_EXNER_C(
     &        P_EXNER_HALF(I,L+1),P_EXNER_HALF(I,L),PUP,PLOW,KAPPA)
            TS=THETA(I,L)*P_EXNER_FULL_L
     &        *(PSTAR(I)/PL(I))**LAPSE_R_OVER_G
            TS=TS*(1.0+C_VIRTUAL*Q(I,1))
            Z(I)=ZH(I,1)+(TS/LAPSE)
     &        *(1.-(P(I)/PSTAR(I))**LAPSE_R_OVER_G)

          ENDIF


          LAST=1 ! Start from bottom level next time

        ELSEIF(K.EQ.P_LEVELS)THEN
!L (iii) Top layer and above: equation (3.6) with local gradient of
!L theta given by equation (3.7) with k=top.
          IF(P_LEVELS.GT.Q_LEVELS)THEN

            PUP1=PSTAR(I)*BKH(P_LEVELS+1) + AKH(P_LEVELS+1)
            PUP =PSTAR(I)*BKH(P_LEVELS)   + AKH(P_LEVELS)
            PLOW=PSTAR(I)*BKH(P_LEVELS-1) + AKH(P_LEVELS-1)
            P_EXNER_FULL_UP= P_EXNER_C(
     &        P_EXNER_HALF(I,P_LEVELS+1),P_EXNER_HALF(I,P_LEVELS),
     &        PUP1,PUP,KAPPA)
            P_EXNER_FULL= P_EXNER_C(
     &        P_EXNER_HALF(I,P_LEVELS),P_EXNER_HALF(I,P_LEVELS-1),
     &        PUP,PLOW,KAPPA)

            T_GRADIENT=( THETA(I,P_LEVELS)*P_EXNER_FULL_UP
     &        -THETA(I,P_LEVELS-1)*P_EXNER_FULL )/
     &        ( P_EXNER_FULL_UP-P_EXNER_FULL )

            LOCAL_GRADIENT=(T_GRADIENT-THETA(I,P_LEVELS))/
     &        P_EXNER_FULL_UP

            SECOND_ORDER_TERM=.5*LOCAL_GRADIENT*( P_EXNER(I)*(P_EXNER(I)
     &        -2.*P_EXNER_FULL_UP)
     &        -P_EXNER_HALF(I,P_LEVELS)*(P_EXNER_HALF(I,P_LEVELS)
     &        -2.*P_EXNER_FULL_UP) )

            DEL_EXNER=P_EXNER_HALF(I,P_LEVELS)-P_EXNER(I)
            THETAV=THETA(I,P_LEVELS)
            Z(I)=ZH(I,P_LEVELS)+
     &        CP_OVER_G*(THETAV*DEL_EXNER-SECOND_ORDER_TERM)

          ELSE

            PUP1=PSTAR(I)*BKH(P_LEVELS+1) + AKH(P_LEVELS+1)
            PUP =PSTAR(I)*BKH(P_LEVELS) + AKH(P_LEVELS)
            PLOW=PSTAR(I)*BKH(P_LEVELS-1) + AKH(P_LEVELS-1)
            P_EXNER_FULL_UP= P_EXNER_C(
     &        P_EXNER_HALF(I,P_LEVELS+1),P_EXNER_HALF(I,P_LEVELS),
     &        PUP1,PUP,KAPPA)
            P_EXNER_FULL= P_EXNER_C(
     &        P_EXNER_HALF(I,P_LEVELS),P_EXNER_HALF(I,P_LEVELS-1),
     &        PUP,PLOW,KAPPA)

            T_GRADIENT=( THETA(I,P_LEVELS)*P_EXNER_FULL_UP
     &        -THETA(I,P_LEVELS-1)*P_EXNER_FULL )/
     &        ( P_EXNER_FULL_UP-P_EXNER_FULL )

            LOCAL_GRADIENT=(T_GRADIENT-THETA(I,P_LEVELS))/
     &        P_EXNER_FULL_UP

            SECOND_ORDER_TERM=.5*LOCAL_GRADIENT*(P_EXNER(I)*(P_EXNER(I)
     &        -2.*P_EXNER_FULL_UP)
     &        -P_EXNER_HALF(I,P_LEVELS)*(P_EXNER_HALF(I,P_LEVELS)
     &        -2.*P_EXNER_FULL_UP) )

            DEL_EXNER=P_EXNER_HALF(I,P_LEVELS)-P_EXNER(I)
            THETAV=THETA(I,P_LEVELS)*(1.0+C_VIRTUAL*Q(I,P_LEVELS))
            Z(I)=ZH(I,P_LEVELS)+
     &        CP_OVER_G*(THETAV*DEL_EXNER-SECOND_ORDER_TERM)

          ENDIF

          LAST=P_LEVELS         ! Start from top level next point
        ELSE
!L 3. Middle layers: equation (3.6) with local theta gradient given by
!L    equation (3.7) with 1 <= k < top.
          K_BELOW=MAX(1,K-1)

          IF(K.GT.Q_LEVELS)THEN

            PUP1   = PSTAR(I)*BKH(K+2) + AKH(K+2)
            PUP    = PSTAR(I)*BKH(K+1) + AKH(K+1)
            PLOW   = PSTAR(I)*BKH(K)   + AKH(K)
            PLOW11 = PSTAR(I)*BKH(K_BELOW+1) + AKH(K_BELOW+1)
            PLOW1  = PSTAR(I)*BKH(K_BELOW)   + AKH(K_BELOW)

            P_EXNER_FULL_UP= P_EXNER_C(
     &        P_EXNER_HALF(I,K+2),P_EXNER_HALF(I,K+1),
     &        PUP1,PUP,KAPPA)
            P_EXNER_FULL= P_EXNER_C(
     &        P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &        PUP,PLOW,KAPPA)
            P_EXNER_FULL_LOW= P_EXNER_C(
     &        P_EXNER_HALF(I,K_BELOW+1),P_EXNER_HALF(I,K_BELOW),
     &        PLOW11,PLOW1,KAPPA)

            T_GRADIENT=( THETA(I,K+1)*P_EXNER_FULL_UP
     &        -THETA(I,K_BELOW)*P_EXNER_FULL_LOW )/
     &        ( P_EXNER_FULL_UP-P_EXNER_FULL_LOW )

            LOCAL_GRADIENT=(T_GRADIENT-THETA(I,K))/P_EXNER_FULL

            SECOND_ORDER_TERM=0.5*LOCAL_GRADIENT
     &        *( P_EXNER(I)*(P_EXNER(I)-2.*P_EXNER_FULL)
     &        -P_EXNER_HALF(I,K)*(P_EXNER_HALF(I,K)-2.*P_EXNER_FULL))

            DEL_EXNER=P_EXNER_HALF(I,K)-P_EXNER(I)
            THETAV=THETA(I,K)
            Z(I)=ZH(I,K)+CP_OVER_G*(THETAV*DEL_EXNER-SECOND_ORDER_TERM)

          ELSE

! Calculation using virtual potential temperature

            PUP1   = PSTAR(I)*BKH(K+2) + AKH(K+2)
            PUP    = PSTAR(I)*BKH(K+1) + AKH(K+1)
            PLOW   = PSTAR(I)*BKH(K)   + AKH(K)
            PLOW11 = PSTAR(I)*BKH(K_BELOW+1) + AKH(K_BELOW+1)
            PLOW1  = PSTAR(I)*BKH(K_BELOW)   + AKH(K_BELOW)

            P_EXNER_FULL_UP= P_EXNER_C(
     &        P_EXNER_HALF(I,K+2),P_EXNER_HALF(I,K+1),
     &        PUP1,PUP,KAPPA)
            P_EXNER_FULL= P_EXNER_C(
     &        P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &        PUP,PLOW,KAPPA)
            P_EXNER_FULL_LOW= P_EXNER_C(
     &        P_EXNER_HALF(I,K_BELOW+1),P_EXNER_HALF(I,K_BELOW),
     &        PLOW11,PLOW1,KAPPA)

            T_GRADIENT=( THETA(I,K+1)*P_EXNER_FULL_UP
     &        -THETA(I,K_BELOW)*P_EXNER_FULL_LOW )/
     &        ( P_EXNER_FULL_UP-P_EXNER_FULL_LOW )

            LOCAL_GRADIENT=(T_GRADIENT-THETA(I,K))/P_EXNER_FULL

            SECOND_ORDER_TERM=0.5*LOCAL_GRADIENT
     &        *(P_EXNER(I)*(P_EXNER(I)-2.*P_EXNER_FULL)
     &        -P_EXNER_HALF(I,K)*(P_EXNER_HALF(I,K)-2.*P_EXNER_FULL) )

            DEL_EXNER=P_EXNER_HALF(I,K)-P_EXNER(I)
            THETAV=THETA(I,K)*(1.0+C_VIRTUAL*Q(I,K))
            Z(I)=ZH(I,K)+CP_OVER_G*(THETAV*DEL_EXNER-SECOND_ORDER_TERM)

          ENDIF
          LAST=K                ! Start from level K next point
        ENDIF
      ENDDO                     ! DO I=START,END

      RETURN
      END
