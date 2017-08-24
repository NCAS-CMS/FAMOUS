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
! SUBROUTINE GW_VERT TO CALCULATE VERTICAL DISTRIBUTION OF GW DRAG VECTR
!
      SUBROUTINE GW_VERT
     1 (PSTAR,PEXNER,THETA,Q,U,V,S_X_STRESS,S_Y_STRESS,START_L,LEVELS
     2 ,Q_LEVELS,POINTS,AKH,BKH,DELTA_AK,DELTA_BK,KAY,KAY_LEE,SD_OROG
     3 ,S_X_OROG,S_Y_OROG,SIGMA_XX,SIGMA_XY,SIGMA_YY,TEST,DU_DT,DV_DT
     4 ,K_LIFT,U_S,V_S,RHO_S
! Diagnostics
     5  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     6  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     7  ,DU_DT_SATN,POINTS_DU_DT_SATN,DU_DT_SATN_ON
     8  ,DV_DT_SATN,POINTS_DV_DT_SATN,DV_DT_SATN_ON
     9  ,DU_DT_JUMP,POINTS_DU_DT_JUMP,DU_DT_JUMP_ON
     &  ,DV_DT_JUMP,POINTS_DV_DT_JUMP,DV_DT_JUMP_ON
     &  ,DU_DT_LEE ,POINTS_DU_DT_LEE ,DU_DT_LEE_ON
     &  ,DV_DT_LEE ,POINTS_DV_DT_LEE ,DV_DT_LEE_ON
     &  ,TRANS_D   ,POINTS_TRANS_D   ,TRANS_D_ON   )

      IMPLICIT NONE
! Description: TO CALCULATE VERTICAL STRESS PROFILE DUE TO SUBGRID-SCALE
!        ANISOTOPIC GRAVITY WAVES AND HENCE DRAG ON MEAN FLOW.
!        HYDRAULIC JUMP IS DIAGNOSED WITH TEST CONTAINING ALPHA.
!        THE HEIGHT OF THE UPSTREAM DIVIDING STREAMLINE IS
!        CALCULATED FOR JUMP POINTS, AND STRESS LINEARISED TO A
!        THIRD OF SURFACE STRESS AT THIS HEIGHT. THE REMAINING
!        WAVES AND NON_JUMP POINTS PROPOGATE VERTICALLY WITH
!        STRESS INDEPENDENT OF HEIGHT UNLESS A CRITICAL LEVEL OR
!        WAVE BREAKING IS DIAGNOSED. THE CRITICAL STRESS IS CALCULATED
!        BY A LAYER SATURATION HYPOTHESIS USING WIND COMPONENT PARALLEL
!             TO THE ORIGINAL SURFACE STRESS INSTEAD OF SURFACE WIND.
!
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE (4),(5),(7),(8),(9)
!
! Current code owner: S.Webster
!
! History:
! Version  Date      Comment
!  4.5   03/06/98   Original Code. Copy of 4.4 GWVERT3A with
!                   operational changes.
!                   Equal acceleration in bottom 3 layers. Correction
!                   to hydaulic jump height. D. Robinson.
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

! Subroutine arguements:

      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,Q_LEVELS            !IN    NUMBER OF WET LEVELS
     *,START_L             !IN    START LEVEL FOR WAVE-BREAKING TEST
     *,POINTS              !IN    NUMBER OF POINTS
     *,K_LIFT(POINTS)      !IN    MODEL LEVEL AT TOP OF BLOCKED LAYER
     *,POINTS_STRESS_UD    !IN    ) No of land points in diagnostic
     *,POINTS_STRESS_VD    !IN    ) arrays for GW stress - u and v
     *,POINTS_DU_DT_SATN   !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_SATN   !IN    ) arrays for GW satn - du and dv
     *,POINTS_DU_DT_JUMP   !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_JUMP   !IN    ) arrays for GW satn - du and dv
     *,POINTS_DU_DT_LEE    !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_LEE    !IN    ) arrays for GW lee - du and dv
     *,POINTS_TRANS_D      !IN    ) No of land points for trans diag

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,THETA(POINTS,LEVELS)             !IN    THETA FIELD
     *,Q(POINTS,Q_LEVELS)               !IN    SATURATION FIELD
     *,U(POINTS,LEVELS)                 !IN    U FIELD
     *,V(POINTS,LEVELS)                 !IN    V FIELD
     *,U_S(POINTS)                      !IN    'SURFACE' U FIELD
     *,V_S(POINTS)                      !IN    'SURFACE' V FIELD
     *,RHO_S(POINTS)                    !IN    'SURFACE' DENSITY
     *,S_X_STRESS(POINTS)               !IN    'SURFACE' X_STRESS
     *,S_Y_STRESS(POINTS)               !IN    'SURFACE' Y_STRESS
     *,S_X_OROG(POINTS)                 !IN    'SURFACE' X_STRESS
     *,S_Y_OROG(POINTS)                 !IN    'SURFACE' Y_STRESS
     *,SIGMA_XX(POINTS)  !IN    DH/DX SQUARED GRADIENT OROGRAPHY
     *,SIGMA_XY(POINTS)  !IN   (DH/DX)(DH/DY) GRADIENT OROGRAPHY
     *,SIGMA_YY(POINTS)  !IN    DH/DY SQUARED GRADIENT OROGRAPHY
     *,TEST(POINTS)      !IN  TEST HYDROLOIC JUMP (SIMILAR TO FROUDE)
     *,SD_OROG(POINTS)   !IN  STANDARD DEVIATION OF OROGRAPHY
!      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
!      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY
     *,DELTA_AK (LEVELS)      !IN    DIFFERENCE ACROSS LAYER
     *,DELTA_BK (LEVELS)      !IN    DIFFERENCE ACROSS LAYER
     *,KAY                    !IN    stress constant (m-1)
     *,KAY_LEE                !IN    TRAPPED LEE WAVE CONSTANT
     *,DU_DT(POINTS,LEVELS)   !OUT   U-ACCELERATION
     *,DV_DT(POINTS,LEVELS)   !OUT   V-ACCELERATION

! Diagnostics
      REAL
     * STRESS_UD(POINTS_STRESS_UD,LEVELS+1) !U STRESS DIAG
     *,STRESS_VD(POINTS_STRESS_VD,LEVELS+1) !V STRESS DIAG
     *,DU_DT_SATN(POINTS_DU_DT_SATN,LEVELS) !U ACCELN DIAG (SATURATION)
     *,DV_DT_SATN(POINTS_DV_DT_SATN,LEVELS) !V ACCELN DIAG (SATURATION)
     *,DU_DT_JUMP(POINTS_DU_DT_JUMP,LEVELS) !U ACCELN DIAG (HYDR JUMP)
     *,DV_DT_JUMP(POINTS_DV_DT_JUMP,LEVELS) !V ACCELN DIAG (HYDR JUMP)
     *,DU_DT_LEE(POINTS_DU_DT_LEE,LEVELS)   !U ACCELN DIAG (LEE WAVE)
     *,DV_DT_LEE(POINTS_DV_DT_LEE,LEVELS)   !V ACCELN DIAG (LEE WAVE)
     *,TRANS_D(POINTS_TRANS_D)   ! TRANSMITTION COEFFICIENT DIAGNOSTIC

      LOGICAL
     * STRESS_UD_ON           !U stress diagnostic switch
     *,STRESS_VD_ON           !V stress diagnostic switch
     *,DU_DT_SATN_ON          !U accel (saturation) diagnostic switch
     *,DV_DT_SATN_ON          !V accel (saturation) diagnostic switch
     *,DU_DT_JUMP_ON          !U accel (hydr jump) diagnostic switch
     *,DV_DT_JUMP_ON          !V accel (hydr jump) diagnostic switch
     *,DU_DT_LEE_ON           !U accel (lee wave) diagnostic switch
     *,DV_DT_LEE_ON           !V accel (lee wave) diagnostic switch
     *,TRANS_D_ON             !Transmittion coefficient diag switch

! Local parameters
      REAL CPBYG
      PARAMETER(CPBYG=CP/G)
! Local scalers
      REAL
     * UCPTSPD              ! |U|COS(.) COMPONENT SPEEED DIRN STRESS
     *,S_STRESS_SQ          ! SURFACE STRESS SQUARE MAGNITUDE
     *,S_STRESS             ! SURFACE STRESS MAGNITUDE
     *,ALPHA1               ! ALLOWS SWAP OF ALPHA AND BETA
     *,BETA1                !             "
     *,SPEED                ! WIND SPEED IN DIR OF STRESS AT LEVEL
     *,N_SQAV               ! AVERAGE OF BRUNT VAISALLA FREQ SQ
     *,NOVERU               ! NBYU FOR ONE LAYER
     *,DEL_EXNER            ! EXNER DIFFERENCE ACROSS LAYER
     *,TEST_CALC            ! CALCULATION FOR JUMP HEIGHT TEST
     *,PU,PL,PB             ! PRESSURES

      LOGICAL   FLAG

      INTEGER   I,K       ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU  ! LEVEL COUNTERS IN ROUTINE
      INTEGER   K_TROP    ! LIMIT OF LEVELS FOR H_JUMP

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: 21  ARRAYS OF FULL FIELD LENGTH
!
      LOGICAL
     * H_JUMP(POINTS)       ! TRUE IF HYDROLIC JUMP REGIME
     *,H_CRIT(POINTS)       ! TRUE IF CRITICAL LEVEL WITHIN JUMP
     *,L_CONT(POINTS)       ! LEVEL CONTINUE
     *,L_LEE(POINTS)        ! TRUE IF TRAPPED LEE WAVE DIAGNOSED

      INTEGER
     * H_O_LEV(POINTS)           ! MODEL LEVEL HEIGHT OF H_JUMP/H_CRIT
     *,K_LEE(POINTS,2)           ! MODEL LEVEL OF TRAPPED LEE WAVE
     *                           ! 'HEIGHT' AND TOP OF WAVE

      REAL
     * NBYU_P(POINTS)         ! U/N FOR CALCULATION OF H_O; AVERAGED
     *,UNIT_X(POINTS)         ! X_COMPNT OF UNIT STRESS VECTOR
     *,UNIT_Y(POINTS)         ! Y_COMPNT OF UNIT STRESS VECTOR
     *,H_O(POINTS)            ! GEOPOTENTIAL HEIGHT ABOVE SURFACE OF
     *                        ! HYDROLIC JUMP
     *,P_EXNER_CENTRE(POINTS,2) ! EXNER PRESSURE AT LAYER CENTRES
     *,N_SQ(POINTS,2)         ! SQUARE OF BRUNT_VAISALA FREQUENCY
     *,ZH(POINTS)             ! TOTAL HEIGHT OF JUMP CALCUALTION
     *,P0(POINTS)             ! PSTAR OR PRESS AT TOP OF K_LIFT
     *,TRANS(POINTS)          ! COEFFICIENT FOR TRANSMITTION OF
     *                        ! SURFACE STRESS
     *,H_LEE(POINTS)          ! TRAPPED LEE WAVE 'HEIGHT' (SEE DOC)
     *,LSQ_LEE(POINTS,2)      ! SCORER PARAMETER AVERAGED BELOW
     *                        ! AND ABOVE TRAPPED LEE WAVE HEIGHT

! Function and subroutine calls:
      EXTERNAL GW_SCOR,GW_SATN,GW_JUMP,GW_LEE
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
!   1.0 START  PRELIMINARIES
! Initialise increment and increment diagnostics
!------------------------------------------------------------
      DO K=1,LEVELS

        DO I=1,POINTS
          DU_DT(I,K)=0.0
          DV_DT(I,K)=0.0
        END DO

        IF( DU_DT_SATN_ON ) THEN
          DO I=1,POINTS
            DU_DT_SATN(I,K)=0.0
          END DO
        ENDIF

        IF( DV_DT_SATN_ON ) THEN
          DO I=1,POINTS
            DV_DT_SATN(I,K)=0.0
          END DO
        ENDIF

        IF( DU_DT_JUMP_ON ) THEN
          DO I=1,POINTS
            DU_DT_JUMP(I,K)=0.0
          END DO
        ENDIF

        IF( DV_DT_JUMP_ON ) THEN
          DO I=1,POINTS
            DV_DT_JUMP(I,K)=0.0
          END DO
        ENDIF

        IF( DU_DT_LEE_ON ) THEN
          DO I=1,POINTS
            DU_DT_LEE(I,K)=0.0
          END DO
        ENDIF

        IF( DV_DT_LEE_ON ) THEN
          DO I=1,POINTS
            DV_DT_LEE(I,K)=0.0
          END DO
        ENDIF

      ENDDO ! Levels
!-----------------------------------------------------------------
!     Code assumes ALPHA < BETA . Swap is possible because of
!     symmetry of calculation( SEE EQN(55), DOC )
!----------------------------------------------------------------
      IF( ALPHA.GT.BETA ) THEN
         ALPHA1 = BETA
         BETA1  = ALPHA
      ELSE
         ALPHA1 = ALPHA
         BETA1  = BETA
      ENDIF

      IF( START_L.LE.2 ) THEN
        WRITE(6,*) 'ERROR G_WAVE: ** START_L MUST BE GREATER THAN 2 ** '
        START_L=3
      ENDIF

      KL=1
      KU=2

      DO I=1,POINTS

!------------------------------------------------------------------
! Calculate logical array for hydraulic jump regime.
!------------------------------------------------------------------
        IF( TEST(I).GE.ALPHA1 ) THEN
          H_JUMP(I)=.TRUE.
        ELSE
          H_JUMP(I)=.FALSE.
        ENDIF
!-------------------------------------------------------------------
! Initialisation. UNIT_X is x_compnt of unit surface stress vector
!-------------------------------------------------------------------
        L_CONT(I) = .TRUE.
        NBYU_P(I) = 0.0
        S_STRESS_SQ = S_X_STRESS(I)**2 + S_Y_STRESS(I)**2
        IF ( S_STRESS_SQ .LE. 0.0 ) THEN
          UNIT_X(I) = 0.0
          UNIT_Y(I) = 0.0
        ELSE
          S_STRESS = SQRT( S_STRESS_SQ )
          UNIT_X(I) = S_X_STRESS(I) / S_STRESS
          UNIT_Y(I) = S_Y_STRESS(I) / S_STRESS
        ENDIF

      ENDDO   ! Points

!--------------------------------------------------------------------
!  2.0 Assess the vertical structure by calculating Scorer parameter
!      for each level. Determine transmittion factor allowing
!      reduction of surface stress from reflection of wave energy
!      off contrast in averaged Scoror profile. Determine trapped
!      lee wave height ( if exists ) and associated paramters
!----------------------------------------------------------------
      CALL GW_SCOR
     1 (PSTAR,PEXNER,THETA,U,V,LEVELS,START_L,H_JUMP,POINTS,AKH,BKH
     2 ,UNIT_X,UNIT_Y,TRANS,K_LEE,H_LEE,LSQ_LEE,L_LEE)

      DO I=1,POINTS
        S_X_STRESS(I)=S_X_STRESS(I)*TRANS(I)
        S_Y_STRESS(I)=S_Y_STRESS(I)*TRANS(I)
      ENDDO

      IF( TRANS_D_ON ) THEN
        DO I=1,POINTS
          TRANS_D(I)=TRANS(I)
        END DO
      ENDIF

!---------------------------------------------------------------------
! 3.0 Find approximate height of tropopause for maximum jump height
!     limit and level limit of orography
!---------------------------------------------------------------------
      FLAG = .TRUE.
      K_TROP = LEVELS-2
      DO K= 3,LEVELS-2
         IF (FLAG) THEN
           PU=100000.*BKH(K+1) + AKH(K+1)
           IF ( PU .LT. 25000. ) THEN
             K_TROP = K
             FLAG = .FALSE.
           ENDIF
         ENDIF
      ENDDO

!---------------------------------------------------------------------
! 3.2 Calculate N by U averaged over levels K_LIFT to a max of K_TROP
!     to test if N/UdeltaZ is greater than 3PI/2. Where this occurs
!     is the jump height, H_O_LEVEL (eqn 8,9)
!     N_SQAV is linearised from N_SQ at layer boundaries
!---------------------------------------------------------------------
      DO K=2,K_TROP
        DO I=1,POINTS
          IF( H_JUMP(I) .AND. L_CONT(I)
     &        .AND. K.GT.K_LIFT(I) ) THEN

            IF( K.EQ.K_LIFT(I)+1 .OR. K_LIFT(I).EQ.0) THEN
              ZH(I)=0.0
              P0(I)=PSTAR(I)*BKH(K_LIFT(I)+1) +AKH(K_LIFT(I)+1)
              PU=PSTAR(I)*BKH(K) + AKH(K)
              PL=PSTAR(I)*BKH(K-1) + AKH(K-1)
! lower layer labelled KU
              P_EXNER_CENTRE(I,KU)=
     &        P_EXNER_C( PEXNER(I,K),PEXNER(I,K-1),PU,PL,KAPPA )
              PL=PU
              PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
! upper layer labelled KL ready for next level stage
              P_EXNER_CENTRE(I,KL)= P_EXNER_C(
     &         PEXNER(I,K+1),PEXNER(I,K),PU,PL,KAPPA)
              N_SQ(I,KL) = G*(THETA(I,K)-THETA(I,K-1))/(THETA(I,K)*
     &         THETA(I,K-1)*(P_EXNER_CENTRE(I,KU)-P_EXNER_CENTRE(I,KL))*
     &         CPBYG)
              IF( N_SQ(I,KL).LE. 0.0 ) THEN
                H_JUMP(I)=.FALSE.
              ENDIF
            ENDIF

! next level stage
            PU=PSTAR(I)*BKH(K+2) + AKH(K+2)
            PL=PSTAR(I)*BKH(K+1) + AKH(K+1)
            P_EXNER_CENTRE(I,KU)=
     &             P_EXNER_C( PEXNER(I,K+2),PEXNER(I,K+1),PU,PL,KAPPA)
            N_SQ(I,KU) = G*(THETA(I,K+1)-THETA(I,K))/(THETA(I,K+1)*
     &       THETA(I,K)*(P_EXNER_CENTRE(I,KL)-P_EXNER_CENTRE(I,KU))*
     &       CPBYG)
            N_SQAV = ( (PEXNER(I,K)-P_EXNER_CENTRE(I,KL))*N_SQ(I,KU) +
     &             (P_EXNER_CENTRE(I,KL) - PEXNER(I,K+1))*N_SQ(I,KL) )
     &                   / ( PEXNER(I,K) - PEXNER(I,K+1) )
            IF( N_SQAV .LE. 0.0 ) THEN
              H_JUMP(I)=.FALSE.
              TEST_CALC = 0.0
            ELSE
!--------------------------------------------------------------------
!   Note U is component parallel to stress vector
!--------------------------------------------------------------------
              UCPTSPD = U(I,K)*UNIT_X(I) + V(I,K)*UNIT_Y(I)
              IF ( UCPTSPD .LE. 0.0 ) THEN
                NOVERU =  0.0
              ELSE
                NOVERU =  SQRT( N_SQAV ) / UCPTSPD
              ENDIF
              IF ( K_LIFT(I).EQ.0 ) THEN
                PB=PSTAR(I)
                DEL_EXNER = PEXNER(I,1) - PEXNER(I,2)
                ZH(I) = CPBYG*THETA(I,1)*DEL_EXNER
                K_LIFT(I)=1
              ELSE
                PB=PSTAR(I)*BKH(K) + AKH(K)
              ENDIF
              NBYU_P(I) = NBYU_P(I) + NOVERU*(PB-PL)
              DEL_EXNER = PEXNER(I,K) - PEXNER(I,K+1)
              ZH(I) = ZH(I) + CPBYG*THETA(I,K)*DEL_EXNER
              TEST_CALC = ZH(I)*NBYU_P(I)/ ( P0(I)-PL )
            ENDIF
!------------------------------------------------------------------
!  Test to see if jump height is reached
!  Note:  (3*PI) / 2 = 4.712389
!  Jump height is defined above LIFT (height of blocked layer)
!------------------------------------------------------------------
            IF( TEST_CALC .GT. 4.712389 ) THEN
              H_O_LEV(I) = K
              L_CONT(I) = .FALSE.

              IF (  H_O_LEV(I) .LE. START_L ) THEN
                H_JUMP(I) = .FALSE.
              ENDIF

            ENDIF   ! Test > 4.712

            IF ( K .EQ. K_TROP .AND. L_CONT(I) ) THEN
              H_O_LEV(I) = K
              L_CONT(I)  = .FALSE.
            ENDIF

          ENDIF ! H_Jump and L_Cont
        ENDDO   ! Points
!  Rename lower centre array as upper centre ready for next level
        KK=KU
        KU=KL
        KL=KK
      ENDDO     ! Levels 2 to K_Trop

!------------------------------------------------------------------
! 3.3 Find if critical layer occurs before H_O_LEV(I)
!------------------------------------------------------------------
       DO I=1,POINTS
          H_CRIT(I)=.FALSE.
       ENDDO

      DO K=START_L+1,LEVELS
        DO I=1,POINTS
          IF( H_JUMP(I) .AND. K.LE.H_O_LEV(I) ) THEN
            SPEED=S_X_STRESS(I)*U(I,K)+S_Y_STRESS(I)*V(I,K)
            IF(SPEED .LE. 0.0) THEN
              H_CRIT(I)=.TRUE.
              H_O_LEV(I)=K
            ENDIF
          ENDIF
        ENDDO   ! Points
      ENDDO     ! Levels 1 to 5

!---------------------------------------------------------------------
!  4.0 If no hydraulic jump the saturation hypothesis is applied from
!      START_L with S_STRESS.
!      Else for jump points, saturation is applied from H_O_LEV with
!      S_STRESS/6. If a critical level is found GW_SATN skipped.
!---------------------------------------------------------------------
      CALL GW_SATN
     1  (PSTAR,PEXNER,THETA,U,V,S_X_STRESS,S_Y_STRESS,START_L,LEVELS
     2  ,POINTS,AKH,BKH,DELTA_AK,DELTA_BK,KAY,SD_OROG,H_O_LEV,H_JUMP
     3  ,H_CRIT,S_X_OROG,S_Y_OROG,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_SATN,POINTS_DU_DT_SATN,DU_DT_SATN_ON
     7  ,DV_DT_SATN,POINTS_DV_DT_SATN,DV_DT_SATN_ON )


!
!------------------------------------------------------------------
! 5.0 Linearize stress profile with pressure up to H_O_LEV and
!     S_STRESS/3 if H_JUMP true. If H_CRIT true then linearise
!     upto zero stress. Skip for non-jump, non-critical points
!------------------------------------------------------------------
      CALL GW_JUMP
     1  (PSTAR,PEXNER,S_X_STRESS,S_Y_STRESS,START_L,LEVELS
     2   ,POINTS,AKH,BKH,DELTA_AK,DELTA_BK,H_O_LEV,H_JUMP
     3   ,H_CRIT,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_JUMP,POINTS_DU_DT_JUMP,DU_DT_JUMP_ON
     7  ,DV_DT_JUMP,POINTS_DV_DT_JUMP,DV_DT_JUMP_ON )


!---------------------------------------------------------------------
! 6.0 Calculate linearized stress profile for trapped lee wave points
!     Lee surface stress is calculated independantly of S_X_STRESS.
!     Lee Stress is distributed vertically upto K_LEE(I,1) where its
!     value at K_LEE(I,1) is reduced by a ratio also calculated
!     within GW_LEE. The remaining stress is deposited by a second
!     gradient, upto K_LEE(I,2). Drags calculated are ADDITIONAL.
!---------------------------------------------------------------------
      CALL GW_LEE
     1  (PSTAR,START_L,LEVELS,POINTS,AKH,BKH,DELTA_AK,DELTA_BK
     2  ,U_S,V_S,RHO_S,L_LEE,LSQ_LEE,H_LEE,K_LEE,KAY_LEE
     3  ,SIGMA_XX,SIGMA_XY,SIGMA_YY,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_LEE,POINTS_DU_DT_LEE,DU_DT_LEE_ON
     7  ,DV_DT_LEE,POINTS_DV_DT_LEE,DV_DT_LEE_ON )

!------------------------------------------------------------------
!  7.0 SET ACCELERATION SAME IN ALL LAYERS 1 UP TO START_L
!------------------------------------------------------------------
      DO KK=1,START_L-1
        DO I=1,POINTS
          DU_DT(I,KK) = DU_DT(I,START_L)
          DV_DT(I,KK) = DV_DT(I,START_L)
        END DO
      END DO

      IF( DU_DT_SATN_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DU_DT_SATN(I,KK) = DU_DT_SATN(I,START_L)
          END DO
        END DO
      ENDIF

      IF( DV_DT_SATN_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DV_DT_SATN(I,KK) = DV_DT_SATN(I,START_L)
          END DO
        END DO
      ENDIF

      IF( DU_DT_JUMP_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DU_DT_JUMP(I,KK) = DU_DT_JUMP(I,START_L)
          END DO
        END DO
      ENDIF

      IF( DV_DT_JUMP_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DV_DT_JUMP(I,KK) = DV_DT_JUMP(I,START_L)
          END DO
        END DO
      ENDIF

      IF( DU_DT_LEE_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DU_DT_LEE(I,KK) = DU_DT_LEE(I,START_L)
          END DO
        END DO
      ENDIF

      IF( DV_DT_LEE_ON ) THEN
        DO KK=1,START_L-1
          DO I=1,POINTS
            DV_DT_LEE(I,KK) = DV_DT_LEE(I,START_L)
          END DO
        END DO
      ENDIF

      RETURN
      END

