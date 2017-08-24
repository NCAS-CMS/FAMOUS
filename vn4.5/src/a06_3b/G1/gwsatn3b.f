C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! SUBROUTINE GW_SATN: SATURATION HYPOTHESIS VERT. STRESS DISTRIBUTION
!
      SUBROUTINE GW_SATN
     1  (PSTAR,PEXNER,THETA,U,V,S_X_STRESS,S_Y_STRESS,START_L,LEVELS
     2  ,POINTS,AKH,BKH,DELTA_AK,DELTA_BK,KAY,SD_OROG,H_O_LEV,H_JUMP
     3  ,H_CRIT,S_X_OROG,S_Y_OROG,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_SATN,POINTS_DU_DT_SATN,DU_DT_SATN_ON
     7  ,DV_DT_SATN,POINTS_DV_DT_SATN,DV_DT_SATN_ON )

      IMPLICIT NONE
! Description:
!             TO CALCULATE STRESS PROFILE DUE TO SUBGRID-SCALE
!             OROGRAPHIC LONG HYDROSTATIC WAVES.
!             THE WAVES PROPOGATE VERTICALLY WITH STRESS INDEPENDENT
!             OF HEIGHT UNLESS A CRITICAL LEVEL OR WAVE BREAKING IS
!             DIAGNOSED. THE CRITICAL STRESS IS CALCULATED
!             FROM WIND COMPONENT PARALLEL TO THE ORIGINAL SURFACE
!             STRESS , NOT NECESSARILY PARALLEL TO SURFACE WIND.
!             THE X AND Y COMPONENTS OF STRESS ARE TREATED
!             INDEPENDANTLY BUT THE VECTOR CAN NOT TURN.
!             IF HYDROLIC JUMP HAS BEEN DIAGNOSED THEN THIS
!             ROUTINE STARTS FROM H_O_LEV AND EQUIVALENT STARTING
!             STRESS OF A THIRD 'SURFACE' STRESS ,UNLESS A
!             CRITICAL LAYER HAS ALREADY BEEN DIAGNOSED.
!             DRAG ON MEAN FLOW IS CALCULATED FROM STRESS PROFILE.
!
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE (1),(2),(3),(4),(6)
!
! Current code owner: S.Webster
!
! History:
! Version  Date      Comment
!  4.5   03/06/98   Original Code. Copy of 4.4 GWSATN3A with operational
!                   changes.
!                   Equal acceleration in bottom 3 layers. Gamma factor
!                   introduced into critical stress formula.
!                   D. Robinson
!
! Code Description:
! Language: Fortran 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
! System component covered: ORIGINAL VERSION FOR CRAY Y-MP
! System task covered: PART OF P22
! SUITABLE FOR SINGLE COLUMN USE,ROTATED GRIDS
! FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY

! Global Variables
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

! Local constants
!
!  Description: This comdeck defines the constants for the 3A version
!               of the Gravity Wave Drag Code. These are tuneable 
!               parameters but are unlikely to be changed.
!
!  History:
!  Version    Date     Comment
!  -------    ----     -------
!    3.4     18/10/94  Original Version    J.R. Mitchell
!    4.3      7/03/97  Remove KAY_LEE (now set in RUNCNST) S.Webster
!    4.5     03/08/98  Add GAMMA_SATN (Used in 06_3B). D. Robinson
!
      REAL ALPHA,BETA,LEE_PHASE
      PARAMETER (
     & ALPHA = 4.0E-1   ! Values for tunable constants of eqn(55) of
     &,BETA  = 1.0      ! GWAVE_93 required in GW_SURF
     &,LEE_PHASE = 0.6  ! Phase across lee height
     & )
      REAL GAMMA_SATN     !  Critical Stress Factor
      PARAMETER (GAMMA_SATN = 0.5)
! Subroutine arguements;

      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,START_L             !IN    START LEVEL FOR WAVE-BREAKING TEST
     *,POINTS              !IN    NUMBER OF POINTS
     *,POINTS_STRESS_UD    !IN    ) No of land points in diagnostic
     *,POINTS_STRESS_VD    !IN    ) arrays for GW stress - u and v
     *,POINTS_DU_DT_SATN   !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_SATN   !IN    ) arrays for GW satn - du and dv
     *,H_O_LEV(POINTS)     !IN    LEVEL OF CRITICAL/JUMP HEIGHT

      LOGICAL
     * H_JUMP(POINTS)      !IN    TRUE IF POINT IS TO BE LINEARIZED
     *,H_CRIT(POINTS)      !IN    TRUE IF CRITICAL HEIGHT BEFORE JUMP
     *,STRESS_UD_ON        !IN U stress diagnostic switch
     *,STRESS_VD_ON        !IN V stress diagnostic switch
     *,DU_DT_SATN_ON       !IN U accel (saturation) diagnostic switch
     *,DV_DT_SATN_ON       !IN V accel (saturation) diagnostic switch

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,THETA(POINTS,LEVELS)             !IN    THETA FIELD
     *,U(POINTS,LEVELS)                 !IN    U FIELD
     *,V(POINTS,LEVELS)                 !IN    V FIELD
     *,S_X_STRESS(POINTS)               !IN    'SURFACE' X_STRESS
     *,S_Y_STRESS(POINTS)               !IN    'SURFACE' Y_STRESS
     *,S_X_OROG(POINTS)                 !IN    'SURFACE' X_OROG
     *,S_Y_OROG(POINTS)                 !IN    'SURFACE' Y_OROG
     *,SD_OROG(POINTS)   !IN  STANDARD DEVIATION OF OROGRAPHY
!      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
!      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY
     *,DELTA_AK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,DELTA_BK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,KAY                    !IN    stress constant (m-1)
     *,DU_DT(POINTS,LEVELS)   !OUT   U-ACCELERATION
     *,DV_DT(POINTS,LEVELS)   !OUT   V-ACCELERATION
! Diagnostics
      REAL
     * DU_DT_SATN(POINTS_DU_DT_SATN,LEVELS)  !U-ACCELN DIAGNOSTIC
     *,DV_DT_SATN(POINTS_DV_DT_SATN,LEVELS)  !V-ACCELN DIAGNOSTIC
     *,STRESS_UD(POINTS_STRESS_UD,LEVELS+1)  !U-STRESS DIAGNOSTIC
     *,STRESS_VD(POINTS_STRESS_VD,LEVELS+1)  !V-STRESS DIAGNOSTIC

! Local parameters
      REAL CPBYG
      PARAMETER(CPBYG=CP/G)

! Local scalers
      REAL
     * RHO                  ! DENSITY AT LAYER BOUNDARY
     *,TB                   ! TEMPERATURE AT LAYER BOUNDARY
     *,DZB                  ! HEIGHT DIFFERENCE ACROSS LAYER BOUNDARY
     *,UB                   ! U-WIND AT LAYER BOUNDARY
     *,VB                   ! V-WIND AT LAYER BOUNDARY
     *,N                    ! BRUNT_VAISALA FREQUENCY
     *,N_SQ                 ! SQUARE OF BRUNT_VAISALA FREQUENCY
     *,C_X_STRESS           ! CRITICAL X_STRESS (EQN 56)
     *,C_Y_STRESS           ! CRITICAL Y_STRESS (EQN 56)
     *,S_STRESS_SQ          ! SQUARE OF SURFACE STRESS
     *,S_STRESS             ! MAGNITUDE OF SURFACE STRESS
     *,SPEEDCALC            ! DOT PRODUCT CALCULATION FOR SPEED/STRESS
     *,DELTA_P              ! DIFFERENCE IN PRESSURE ACROSS LAYER
     *,ALPHA1               ! ALLOWS SWAP OF ALPHA AND BETA
     *,BETA1                !             "
     *,GAMMA_SQ             ! Parameter for scaling critical stress
     *,DELTA_AK_SUM ! DELTA_AK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,DELTA_BK_SUM ! DELTA_BK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,PU,PL,P_EXNER_CENTRE
      INTEGER   I,K    ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU,KT ! LEVEL COUNTERS IN ROUTINE
      INTEGER   H_O_L ! DUMMY FOR H_O_LEV(I)

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: 11  ARRAYS OF FULL FIELD LENGTH
!
      REAL
     * DZ(POINTS,3)         ! HEIGHT DIFFERENCES IN EACH HALF LAYER
     *,T(POINTS,2)          ! TEMPERATURES (LEVELS)
     *,X_STRESS(POINTS,2)   ! X_STRESSES (LAYER BOUNDARIES)
     *,Y_STRESS(POINTS,2)   ! Y_STRESSES (LAYER BOUNDARIES)
     *,X_S_CONST(POINTS)    ! LEVEL INDEPEDANT CONSTS FOR CALCULATION
     *,Y_S_CONST(POINTS)    ! OF CRITICAL STRESSES

! Function and subroutine calls
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


!-------------------------------------------------------------------
!    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
!   1. START LEVEL  PRELIMINARIES
!------------------------------------------

CFPP$ NOCONCUR L
!      TREAT LAYERS BELOW AND INCLUDING START_L AS ONE LAYER
        DELTA_AK_SUM = 0.0
        DELTA_BK_SUM = 0.0
      DO K=1,START_L
        DELTA_AK_SUM = DELTA_AK_SUM + DELTA_AK(K)
        DELTA_BK_SUM = DELTA_BK_SUM + DELTA_BK(K)
      END DO
CFPP$ CONCUR

!-----------------------------------------------------------------
!     CODE ASSUMES ALPHA < BETA . SWAP IS POSSIBLE BECAUSE OF
!     SYMMETRY OF CALCULATION ( SEE EQN(55), DOC )
!----------------------------------------------------------------
      IF( ALPHA.GT.BETA ) THEN
         ALPHA1 = BETA
         BETA1  = ALPHA
      ELSE
         ALPHA1 = ALPHA
         BETA1  = BETA
      ENDIF

      GAMMA_SQ = GAMMA_SATN * GAMMA_SATN

      KL=1
      KU=2
      KT=3

      DO I=1,POINTS

        IF( H_JUMP(I) .AND. H_CRIT(I) ) THEN
          X_STRESS(I,KL) = 0.0
          Y_STRESS(I,KL) = 0.0
          H_O_L=H_O_LEV(I)

        ELSE IF( H_JUMP(I) ) THEN
          X_STRESS(I,KL) = S_X_STRESS(I)/6.0
          Y_STRESS(I,KL) = S_Y_STRESS(I)/6.0
          H_O_L=H_O_LEV(I)
          PU=PSTAR(I)*BKH(H_O_L+1) + AKH(H_O_L+1)
          PL=PSTAR(I)*BKH(H_O_L) + AKH(H_O_L)
          P_EXNER_CENTRE=
     &    P_EXNER_C( PEXNER(I,H_O_L+1),PEXNER(I,H_O_L),PU,PL,KAPPA)

          DZ(I,KL)     = (P_EXNER_CENTRE - PEXNER(I,H_O_L+1))
     *                   *THETA(I,H_O_L)*CPBYG
          T(I,KL)      = P_EXNER_CENTRE*THETA(I,H_O_L)
          DZ(I,KT)     = DZ(I,KL)
          T(I,KU)      = T(I,KL)

        ELSE
          X_STRESS(I,KL) = S_X_STRESS(I)
          Y_STRESS(I,KL) = S_Y_STRESS(I)
          PU=PSTAR(I)*BKH(START_L+1) + AKH(START_L+1)
          PL=PSTAR(I)*BKH(START_L) + AKH(START_L)
          P_EXNER_CENTRE=
     &    P_EXNER_C( PEXNER(I,START_L+1),PEXNER(I,START_L),PU,PL,KAPPA)

          DZ(I,KL)     = (P_EXNER_CENTRE - PEXNER(I,START_L+1))
     *                   *THETA(I,START_L)*CPBYG
          T(I,KL)      = P_EXNER_CENTRE*THETA(I,START_L)

        ENDIF

!------------------------------------------------------------------
! 1.1 CALCULATE LEVEL INDEPENDANT STRESS CONSTANTS FOR  SECTION 2.2
!------------------------------------------------------------------
        S_STRESS_SQ = S_X_STRESS(I)**2 + S_Y_STRESS(I)**2
        S_STRESS=SQRT(S_STRESS_SQ)
        IF((BETA*SD_OROG(I)*SD_OROG(I)*S_STRESS_SQ*S_STRESS).LE.1.0E-30
     &     .OR.SD_OROG(I).LE.0.0 .OR. S_STRESS_SQ.LE.0.0 )THEN
          X_S_CONST(I) = 0.0
          Y_S_CONST(I) = 0.0
        ELSE
          Y_S_CONST(I) = KAY*ALPHA*GAMMA_SQ/
     *              (BETA*SD_OROG(I)*SD_OROG(I)*S_STRESS_SQ*S_STRESS)
          X_S_CONST(I) = Y_S_CONST(I)*S_X_OROG(I)
          Y_S_CONST(I) = Y_S_CONST(I)*S_Y_OROG(I)


        ENDIF

      END DO

      IF( STRESS_UD_ON ) THEN
        DO I=1,POINTS
          IF( H_JUMP(I) ) THEN
            STRESS_UD(I,H_O_LEV(I)) = X_STRESS(I,KL)
          ELSE
            STRESS_UD(I,START_L) = X_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

      IF( STRESS_VD_ON ) THEN
        DO I=1,POINTS
          IF( H_JUMP(I) ) THEN
            STRESS_VD(I,H_O_LEV(I)) = Y_STRESS(I,KL)
          ELSE
            STRESS_VD(I,START_L) = Y_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

!------------------------------------------------------------------
!   2  LOOP LEVELS
!------------------------------------------------------------------

      DO K=START_L+1,LEVELS


        DO I=1,POINTS

          X_STRESS(I,KU) = X_STRESS(I,KL)
          Y_STRESS(I,KU) = Y_STRESS(I,KL)

          IF( K .EQ. START_L+1 ) THEN
            DELTA_P = DELTA_AK_SUM+DELTA_BK_SUM*PSTAR(I)
          ELSE
            DELTA_P = DELTA_AK(K-1)+DELTA_BK(K-1)*PSTAR(I)
          END IF

          IF( (.NOT.H_JUMP(I)) .OR. K.GT.H_O_LEV(I) ) THEN

            IF( (X_STRESS(I,KL) .NE. 0.0)
     *        .OR.(Y_STRESS(I,KL) .NE. 0.0) )   THEN

              PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
              PL=PSTAR(I)*BKH(K) + AKH(K)
              P_EXNER_CENTRE=
     &                P_EXNER_C( PEXNER(I,K+1),PEXNER(I,K),PU,PL,KAPPA)

! lower half height of upper layer
              DZ(I,KU)    = (PEXNER(I,K) - P_EXNER_CENTRE)*THETA(I,K)
     *                       *CPBYG
! upper half height of upper layer
              DZ(I,KT)    = (P_EXNER_CENTRE - PEXNER(I,K+1))*THETA(I,K)
     *                       *CPBYG
! model level height difference
              DZB          = DZ(I,KU) + DZ(I,KL)
              UB       = (DZ(I,KU)*U(I,K-1)+DZ(I,KL)*U(I,K)) / DZB
              VB       = (DZ(I,KU)*V(I,K-1)+DZ(I,KL)*V(I,K)) / DZB
              T(I,KU)  = P_EXNER_CENTRE*THETA(I,K)
              TB       = (DZ(I,KU)*T(I,KL) + DZ(I,KL)*T(I,KU))/DZB
              RHO      = ( AKH(K) + BKH(K)*PSTAR(I) )/(R*TB)

!------------------------------------------------------------------
!            2.2 CALCULATE BRUNT-VAISALA FREQUENCY
!------------------------------------------------------------------

              N_SQ = G*( THETA(I,K) - THETA(I,K-1) )*PEXNER(I,K)/
     *                   ( TB*DZB )

              IF( N_SQ .LE. 0.0 ) THEN
!            SET STRESS TO ZERO IF UNSTABLE
                N_SQ = 0.0
                X_STRESS(I,KU) = 0.0
                Y_STRESS(I,KU) = 0.0
              ELSE
                N   = SQRT( N_SQ )
                SPEEDCALC = UB*S_X_STRESS(I) + VB*S_Y_STRESS(I)
                C_Y_STRESS = (SPEEDCALC**3)*RHO/N
                C_X_STRESS = X_S_CONST(I)*C_Y_STRESS
                C_Y_STRESS = Y_S_CONST(I)*C_Y_STRESS

!------------------------------------------------------------------
!           2.3    CALCULATE CRITICAL STRESS FOR
!                EACH COMPONENT    (EQN 6)
!                  TEST FOR WAVE-BREAKING
!                AND MODIFY STRESS AT UPPER LAYER BOUNDARY
!------------------------------------------------------------------

                IF( X_STRESS(I,KL) .GT. 0.0 ) THEN
                  IF( C_X_STRESS .LT. 0.0 ) THEN
                    C_X_STRESS     = 0.0
                  ENDIF

                  IF( C_X_STRESS .LT. X_STRESS(I,KU) ) THEN
                    X_STRESS(I,KU) = C_X_STRESS
                  ENDIF
                ENDIF
                IF( X_STRESS(I,KL) .LT. 0.0 ) THEN
                  IF( C_X_STRESS .GT. 0.0 ) THEN
                    C_X_STRESS     = 0.0
                  ENDIF

                  IF( C_X_STRESS .GT. X_STRESS(I,KU) ) THEN
                    X_STRESS(I,KU) = C_X_STRESS
                  ENDIF
                ENDIF

                IF( Y_STRESS(I,KL) .GT. 0.0 ) THEN
                  IF( C_Y_STRESS .LT. 0.0 ) THEN
                    C_Y_STRESS     = 0.0
                  ENDIF

                  IF( C_Y_STRESS .LT. Y_STRESS(I,KU) ) THEN
                    Y_STRESS(I,KU) = C_Y_STRESS
                  ENDIF
                ENDIF
                IF( Y_STRESS(I,KL) .LT. 0.0 ) THEN
                  IF( C_Y_STRESS .GT. 0.0 ) THEN
                    C_Y_STRESS     = 0.0
                  ENDIF

                  IF( C_Y_STRESS .GT. Y_STRESS(I,KU) ) THEN
                    Y_STRESS(I,KU) = C_Y_STRESS
                  ENDIF
                ENDIF

              END IF     ! (N_SQ < 0) ELSE N_SQ > 0

            END IF     ! STRESS X OR Y NE 0

          END IF     ! no jump or above jump height

!------------------------------------------------------------------
!              2.4 CALCULATE DRAG FROM VERTICAL STRESS CONVERGENCE
!                AND ACCELERATIONS FOR WIND COMPONENTS
!------------------------------------------------------------------

          DU_DT(I,K-1) = G*(X_STRESS(I,KL) - X_STRESS(I,KU))/DELTA_P
          DV_DT(I,K-1) = G*(Y_STRESS(I,KL) - Y_STRESS(I,KU))/DELTA_P

        END DO

! Diagnostics
        IF( STRESS_UD_ON ) THEN
          DO I=1,POINTS
            STRESS_UD(I,K) = X_STRESS(I,KU)
          END DO
        ENDIF

        IF( STRESS_VD_ON ) THEN
          DO I=1,POINTS
            STRESS_VD(I,K) = Y_STRESS(I,KU)
          END DO
        ENDIF

        IF( DU_DT_SATN_ON ) THEN
          DO I=1,POINTS
            DU_DT_SATN(I,K-1) = DU_DT(I,K-1)
          END DO
        ENDIF

        IF( DV_DT_SATN_ON ) THEN
          DO I=1,POINTS
            DV_DT_SATN(I,K-1) = DV_DT(I,K-1)
          END DO
        ENDIF

! Swap storage for lower and upper layers
        KK=KL
        KL=KU
        KU=KK

! Replace top half height of lower layer ready for next pass
        DO I=1,POINTS
          DZ(I,KL)=DZ(I,KT)
        END DO

      END DO
!   END LOOP LEVELS

!
!------------------------------------------------------------------
!  3.0 TOP OF MODEL. SET ACCELERATION SAME AS PENULTIMATE LAYER
!      WITH PROVISO THAT STRESS COMPONENTS DO NOT PASS THROUGH 0
!------------------------------------------------------------------

      DO I=1,POINTS
        DELTA_P   = DELTA_AK(LEVELS) + DELTA_BK(LEVELS)*PSTAR(I)

        X_STRESS(I,KU) = X_STRESS(I,KL) - DU_DT(I,LEVELS-1)*DELTA_P/G
        IF( (X_STRESS(I,KU).LT.0.0) .AND. (X_STRESS(I,KL).GT.0.0) )
     &    X_STRESS(I,KU) = 0.0
        IF( (X_STRESS(I,KU).GT.0.0) .AND. (X_STRESS(I,KL).LT.0.0) )
     &    X_STRESS(I,KU) = 0.0

        Y_STRESS(I,KU) = Y_STRESS(I,KL) - DV_DT(I,LEVELS-1)*DELTA_P/G
        IF( (Y_STRESS(I,KU).LT.0.0) .AND. (Y_STRESS(I,KL).GT.0.0) )
     &    Y_STRESS(I,KU) = 0.0
        IF( (Y_STRESS(I,KU).GT.0.0) .AND. (Y_STRESS(I,KL).LT.0.0) )
     &    Y_STRESS(I,KU) = 0.0

        DU_DT(I,LEVELS) = G*(X_STRESS(I,KL) - X_STRESS(I,KU))/DELTA_P
        DV_DT(I,LEVELS) = G*(Y_STRESS(I,KL) - Y_STRESS(I,KU))/DELTA_P

      END DO

! Diagnostics
      IF( STRESS_UD_ON ) THEN
        DO I=1,POINTS
          STRESS_UD(I,LEVELS+1) = X_STRESS(I,KU)
        END DO
      ENDIF

      IF( STRESS_VD_ON ) THEN
        DO I=1,POINTS
          STRESS_VD(I,LEVELS+1) = Y_STRESS(I,KU)
        END DO
      ENDIF

      IF( DU_DT_SATN_ON ) THEN
        DO I=1,POINTS
          DU_DT_SATN(I,LEVELS) = DU_DT(I,LEVELS)
        END DO
      ENDIF

      IF( DV_DT_SATN_ON ) THEN
        DO I=1,POINTS
          DV_DT_SATN(I,LEVELS) = DV_DT(I,LEVELS)
        END DO
      ENDIF

      RETURN
      END

