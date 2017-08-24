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
! SUBROUTINE GW_JUMP TO CALCULATE VERT DISTRIBUTION BELOW A JUMP HEIGHT
!
      SUBROUTINE GW_JUMP
     1  (PSTAR,PEXNER,S_X_STRESS,S_Y_STRESS,START_L,LEVELS
     2   ,POINTS,AKH,BKH,DELTA_AK,DELTA_BK,H_O_LEV,H_JUMP
     3   ,H_CRIT,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_JUMP,POINTS_DU_DT_JUMP,DU_DT_JUMP_ON
     7  ,DV_DT_JUMP,POINTS_DV_DT_JUMP,DV_DT_JUMP_ON )

      IMPLICIT NONE
! Description:
!             TO CALCULATE STRESS PROFILE DUE TO SUBGRID-SCALE
!             OROGRAPHIC GRAVITY WAVES FOR HIGH DRAG HYDRAULIC JUMP
!             STATES OR IF CRITICAL LAYER FOUND WITHIN JUMP.
!             STRESS AND DRAG LINEARISED WITH PRESSURE UPTO JUMP/
!             CRITICAL HEIGHT FROM STARTING LEVEL.
!
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE (???) TO (???)
!
! Current code owner: S.Webster
!
! History:
! Version  Date      Comment
!  4.5   03/06/98   Original Code. Copy of 4.4 GWJUMP3A with operational
!                   changes.
!                   Equal acceleration in bottom 3 layers. Ratio for
!                   Hydraulic Jump Stress increased to 0.833.
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
     *,POINTS_DU_DT_JUMP   !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_JUMP   !IN    ) arrays for GW satn - du and dv
     *,H_O_LEV(POINTS)     !IN    LEVEL OF CRITICAL/JUMP HEIGHT

      LOGICAL
     * H_JUMP(POINTS)      !IN    TRUE IF POINT IS TO BE LINEARIZED
     *,H_CRIT(POINTS)      !IN    TRUE IF CRITICAL HEIGHT BEFORE JUMP
     *,STRESS_UD_ON           !U stress diagnostic switch
     *,STRESS_VD_ON           !V stress diagnostic switch
     *,DU_DT_JUMP_ON          !U accel (hydr jump) diagnostic switch
     *,DV_DT_JUMP_ON          !V accel (hydr jump) diagnostic switch

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,S_X_STRESS(POINTS)               !IN    'SURFACE' X_STRESS
     *,S_Y_STRESS(POINTS)               !IN    'SURFACE' Y_STRESS
!      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
!      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY
     *,DELTA_AK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,DELTA_BK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,DU_DT(POINTS,LEVELS)   !OUT   U-ACCELERATION
     *,DV_DT(POINTS,LEVELS)   !OUT   V-ACCELERATION
! Diagnostics
      REAL
     * DU_DT_JUMP(POINTS_DU_DT_JUMP,LEVELS) !U-ACCELN  DIAGNOSTIC
     *,DV_DT_JUMP(POINTS_DV_DT_JUMP,LEVELS) !V-ACCELN  DIAGNOSTIC
     *,STRESS_UD(POINTS_STRESS_UD,LEVELS+1) !U STRESS  DIAGNOSTIC
     *,STRESS_VD(POINTS_STRESS_VD,LEVELS+1) !V STRESS  DIAGNOSTIC

! Local parameters
      REAL CPBYG
      PARAMETER(CPBYG=CP/G)

! Local scalers
      REAL
     * DELTA_P              ! DIFFERENCE IN PRESSURE ACROSS LAYER(S)
     *,DELTA_AK_SUM ! DELTA_AK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,DELTA_BK_SUM ! DELTA_BK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,PU,PL
      INTEGER   I,K    ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU,KT ! LEVEL COUNTERS IN ROUTINE

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: 6  ARRAYS OF FULL FIELD LENGTH
!

      REAL
     * X_STRESS(POINTS,2)   ! X_STRESSES (LAYER BOUNDARIES)
     *,Y_STRESS(POINTS,2)   ! Y_STRESSES (LAYER BOUNDARIES)
     *,DP_X_STRESS(POINTS)  ! STRESS GRADIENT
     *                      ! SURFACE X_STRESS
     *,DP_Y_STRESS(POINTS)  ! STRESS GRADIENT
     *                      ! SURFACE Y_STRESS

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
!   1. START LEVEL  PRELIMINARIES
!-------------------------------------------------------------------

CFPP$ NOCONCUR L
!      TREAT LAYERS BELOW AND INCLUDING START_L AS ONE LAYER
        DELTA_AK_SUM = 0.0
        DELTA_BK_SUM = 0.0
      DO K=1,START_L
        DELTA_AK_SUM = DELTA_AK_SUM + DELTA_AK(K)
        DELTA_BK_SUM = DELTA_BK_SUM + DELTA_BK(K)
      END DO
CFPP$ CONCUR

      KL=1
      KU=2
      KT=3

      DO I=1,POINTS
        IF( H_JUMP(I) ) THEN

          PU=PSTAR(I)*BKH(H_O_LEV(I)) + AKH(H_O_LEV(I))
          PL=PSTAR(I)*BKH(START_L) + AKH(START_L)
          DELTA_P= PL - PU

          IF( H_CRIT(I) ) THEN
            DP_X_STRESS(I) = S_X_STRESS(I)/ ( DELTA_P )
            DP_Y_STRESS(I) = S_Y_STRESS(I)/ ( DELTA_P )
          ELSE
            DP_X_STRESS(I) = 0.83333333*S_X_STRESS(I)/
     &                                   ( DELTA_P )
            DP_Y_STRESS(I) = 0.83333333*S_Y_STRESS(I)/
     &                                   ( DELTA_P )
          ENDIF

        ENDIF ! if H_JUMP

      END DO  ! Points

      IF( STRESS_UD_ON ) THEN
        DO I=1,POINTS
          IF( H_JUMP(I) ) THEN
            X_STRESS(I,KL) = S_X_STRESS(I)
            STRESS_UD(I,START_L) = X_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

      IF( STRESS_VD_ON ) THEN
        DO I=1,POINTS
          IF( H_JUMP(I) ) THEN
            Y_STRESS(I,KL) = S_Y_STRESS(I)
            STRESS_VD(I,START_L) = Y_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

!------------------------------------------------------------------
!    2 LOOP LEVELS
!------------------------------------------------------------------

      DO K=START_L+1,LEVELS


        DO I=1,POINTS
          IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN

            IF( K .EQ. START_L+1 ) THEN
              DELTA_P =(DELTA_AK(START_L)+DELTA_BK(START_L)*PSTAR(I))
     &               /( DELTA_AK_SUM     +DELTA_BK_SUM*PSTAR(I) )
              DU_DT(I,START_L) = - G*DP_X_STRESS(I)*DELTA_P
              DV_DT(I,START_L) = - G*DP_Y_STRESS(I)*DELTA_P
            ELSE
              DU_DT(I,K-1) = - G*DP_X_STRESS(I)
              DV_DT(I,K-1) = - G*DP_Y_STRESS(I)
            ENDIF

          ENDIF   ! if H_JUMP(I) .and. K<=H_O_LEV(I)

        END DO ! points

! Diagnostics
        IF( DU_DT_JUMP_ON ) THEN
          DO I=1,POINTS
            IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN
              DU_DT_JUMP(I,K-1) = DU_DT(I,K-1)
            ENDIF
          END DO
        ENDIF

        IF( DV_DT_JUMP_ON ) THEN
          DO I=1,POINTS
            IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN
              DV_DT_JUMP(I,K-1) = DV_DT(I,K-1)
            ENDIF
          END DO
        ENDIF

        IF( STRESS_UD_ON ) THEN
          DO I=1,POINTS
            IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN
              DELTA_P = DELTA_AK(K-1)+DELTA_BK(K-1)*PSTAR(I)
              X_STRESS(I,KU) = X_STRESS(I,KL)
! Stress at upper boundary  N.B. DELTA_P is -VE
              X_STRESS(I,KU) = X_STRESS(I,KL)+DP_X_STRESS(I)*DELTA_P
              STRESS_UD(I,K) = X_STRESS(I,KU)
! Top of model. Set acceln same as penultimate layer if stress >0
!     FOR COMPLETION- TOP LEVEL NEVER REACHED IN H_JUMP
              IF( K .EQ. LEVELS ) THEN
                DELTA_P   = DELTA_AK(LEVELS) + DELTA_BK(LEVELS)*PSTAR(I)
                X_STRESS(I,KL) = X_STRESS(I,KU) - DU_DT(I,LEVELS-1)
     *                                                     *DELTA_P/G
                IF( X_STRESS(I,KL) .LT. 0.0) X_STRESS(I,KL) = 0.0
                DU_DT(I,LEVELS) = G*(X_STRESS(I,KU) - X_STRESS(I,KL))
     *                                                      /DELTA_P
                STRESS_UD(I,LEVELS+1) = X_STRESS(I,KL)
              ENDIF ! top of model
            ENDIF   ! H_Jump & K < H_O_Lev
          END DO
        ENDIF       ! Stress_ud on

        IF( STRESS_VD_ON ) THEN
          DO I=1,POINTS
            IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN
              DELTA_P = DELTA_AK(K-1)+DELTA_BK(K-1)*PSTAR(I)
              Y_STRESS(I,KU) = Y_STRESS(I,KL)
! Stress at upper boundary  N.B. DELTA_P is -VE
              Y_STRESS(I,KU) = Y_STRESS(I,KL)+DP_Y_STRESS(I)*DELTA_P
              STRESS_VD(I,K) = Y_STRESS(I,KU)
! Top of model. Set acceln same as penultimate layer if stress >0
!     FOR COMPLETION- TOP LEVEL NEVER REACHED IN H_JUMP
              IF( K .EQ. LEVELS ) THEN
                DELTA_P   = DELTA_AK(LEVELS) + DELTA_BK(LEVELS)*PSTAR(I)
                Y_STRESS(I,KL) = Y_STRESS(I,KU) - DV_DT(I,LEVELS-1)
     *                                                     *DELTA_P/G
                IF( Y_STRESS(I,KL) .LT. 0.0) Y_STRESS(I,KL) = 0.0
                DV_DT(I,LEVELS) = G*(Y_STRESS(I,KU) - Y_STRESS(I,KL))
     *                                                      /DELTA_P
                STRESS_VD(I,LEVELS+1) = Y_STRESS(I,KL)
              ENDIF ! top of model
            ENDIF   ! H_Jump & K < H_O_Lev
          END DO
        ENDIF       ! Stress_vd on

! Swap storage for lower and upper layers
        KK=KL
        KL=KU
        KU=KK

      END DO
!   END LOOP LEVELS

      RETURN
      END
