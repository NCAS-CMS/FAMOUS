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
!  SUBROUTINE GW_SURF TO CALCULATE SURFACE STRESS VECTOR FOR GWD
!
      SUBROUTINE GW_SURF
     1  (PSTAR,PEXNER,THETA,U,V,SD_OROG,SIGMA_XX,SIGMA_XY,SIGMA_YY,
     2   S_X_STRESS,S_Y_STRESS,S_X_OROG,S_Y_OROG,LEVELS,
     3   POINTS,AK,BK,AKH,BKH,KAY,TEST,K_LIFT,U_S,V_S,RHO_S)

      IMPLICIT NONE
! Description:
! Finds lift height for a blocked layer. Averages wind & stability from
! the largest of lift height or level 2, upto a height of 1.4*sig(h) or
! lev 3. Calculates anisotropic 'surface' stress due to sub-grid scale
! orography via squared gradients and variance. The wave amplitude is
! dependant on variance and may be proportional or follow a square or
! cubic law w.r.t wind speed. These regimes are controlled by
! parameters alpha and beta.
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE (1),(2)
!
! Current code owner: S.Webster
!
! History:
! Version  Date      Comment
!  3.4   18/10/94   Original Code. J.R.Mitchell
!  4.4   19/09/97   Remove *IF -DEF,CRAY compile options. S.Webster 
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

! Subroutine arguements
      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,POINTS              !IN    NUMBER OF POINTS
     *,K_LIFT(POINTS)      !OUT   MODEL LEVEL HEIGHT OF BLOCKED LAYER

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,THETA(POINTS,LEVELS)             !IN    THETA FIELD
     *,U(POINTS,LEVELS)                 !IN    U FIELD
     *,V(POINTS,LEVELS)                 !IN    V FIELD
!            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
      REAL
     * AK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,BK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER boundary
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER boundary
     *,KAY                    !IN    surface stress constant (m-1)
     *,SD_OROG(POINTS)        !IN    STANDARD DEVIATION OF OROGRAPHY
     *,SIGMA_XX(POINTS)       !IN    DH/DX SQUARED GRADIENT OROG
     *,SIGMA_XY(POINTS)       !IN    (DH/DX)(DH/DY) GRADIENT OROG
     *,SIGMA_YY(POINTS)       !IN    DH/DY SQUARED GRADIENT OROG
     *,S_X_STRESS(POINTS)     !OUT   'SURFACE' STRESS IN X-DIRN
     *,S_Y_STRESS(POINTS)     !OUT   'SURFACE' STRESS IN Y-DIRN
     *,S_X_OROG(POINTS)       !OUT   'SURFACE' STRESS/OROG IN X-DIRN
     *,S_Y_OROG(POINTS)       !OUT   'SURFACE' STRESS/OROG IN Y-DIRN
     *,TEST(POINTS)       !OUT  TEST HYDROLOIC JUMP (SIMILAR TO FROUDE)
     *,U_S(POINTS)        !OUT  U-WINDS AVERAGED OVER 'SURFACE'
     *,V_S(POINTS)        !OUT  V-WINDS AVERAGED OVER 'SURFACE'
     *,RHO_S(POINTS)      !OUT DENSITY AVERAGED OVER 'SURFACE'

! Local parameters
      REAL CPBYG
      PARAMETER(CPBYG=CP/G)
! Local scalers
      LOGICAL   FLAG

      INTEGER   I          ! LOOP COUNTER IN ROUTINE
      INTEGER   K,KK,KL,KU ! LEVEL COUNTERS IN ROUTINE
      INTEGER   K_MOUNT    ! LIMIT OF MOUNTAIN LEVELS

      REAL
     * RHO          ! DENSITY
     *,SPEED        ! WIND SPEED / WIND SPEED IN DIRN STRESS
     *,SPEEDCALC    ! NUMERATOR OF CALCUATION FOR SPEED
     *,S_STRESS_SQ  ! DENOMINATER OF CALCUATION FOR SPEED
     *,N            ! BRUNT_VAISALA FREQUENCY
     *,N_SQAV       ! AVERAGE OF BRUNT VAISALLA FREQ SQ
     *,H_SQ         ! H1*H2 OF NEW PARMS FORMULA (55)
     *,ALPHA1       ! DUMMY FOR ALPHA (C_GWAWE)
     *,BETA1        ! DUMMY FOR BETA (C_GWAVE)
     *,R_ALPHA      ! RECIPRICOL OF ALPHA1
     *,R_BETA       ! RECIPRICOL OF BETA1
     *,CALC         ! CALCULATION FOR SURFACE STRESS MAGNITUDE
     *,UOVERN       ! UBYN FOR ONE LAYER
     *,DEL_EXNER    ! EXNER DIFFERENCE ACROSS LAYER
     *,DELTA_P      ! PRESSURE DIFFERENCE
     *,PU,PL        ! PRESSURES AT LAYER BOUNDARIES

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: 13  ARRAYS OF FULL FIELD LENGTH
!
      LOGICAL
     * L_CONT(POINTS)       ! LEVEL CONTINUE
     *,L_CONT2(POINTS)      ! LEVEL CONTINUE

      INTEGER
     * K_BOT(POINTS)          ! MODEL LEVEL HEIGHT OF BLOCKED LAYER
     *,K_TOP(POINTS)          ! MODEL LEVEL AT MOUNTAIN TOPS

      REAL
     * LIFT(POINTS)           ! OROGRAPHIC CONTRIBUTION TO H_JUMP HEIGHT
     *,P_EXNER_CENTRE(POINTS,2) ! EXNER PRESSURE AT LAYER CENTRES
     *,N_SQ(POINTS,2)         ! SQUARE OF BRUNT_VAISALA FREQUENCY
     *,ZH(POINTS)             ! TOTAL HEIGHT ABOVE GROUND (M)
     *,DELTA_P_SUM(POINTS)    ! DELTA PRESSURE CUMULATIVE SUM
     *,NSQ_S(POINTS)          ! N SQUARED AVERAGED OVER 'SURFACE'

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


!-----------------------------------------------------------------
!     CODE ASSUMES ALPHA < BETA . SWAP IS POSSIBLE BECAUSE OF
!     SYMMETRY OF CALCULATION ( SEE EQN(55), DOC )
!----------------------------------------------------------------
      IF( ALPHA.GT.BETA ) THEN
         ALPHA1=BETA
         BETA1=ALPHA
      ELSE
         ALPHA1=ALPHA
         BETA1=BETA
      ENDIF


      R_ALPHA=1.0/ALPHA1
      R_BETA=1.0/BETA1

      KL=1
      KU=2

!---------------------------------------------------------------------
! 1.0 Find approximate level limit of orography. Initialation
!---------------------------------------------------------------------
      FLAG = .TRUE.
      K_MOUNT = LEVELS-2
      DO K= 3,LEVELS-2
         IF (FLAG) THEN
           PU=100000.*BKH(K+1) + AKH(K+1)
           IF ( PU .LT. 50000. ) THEN
             K_MOUNT = K
             FLAG = .FALSE.
           ENDIF
         ENDIF
      ENDDO
      DO I=1,POINTS
        DELTA_P_SUM(I) = 0.0
        NSQ_S(I)       = 0.0
        RHO_S(I)       = 0.0
        U_S(I)         = 0.0
        V_S(I)         = 0.0
        L_CONT(I)     =.TRUE.
        L_CONT2(I)    =.TRUE.
      ENDDO
!---------------------------------------------------------------------
! 1.1 Calculate N by U averaged over boundary values of level 2
!     Use this to calculate level depth of any blocked layer -
!     K_LIFT (eqn 7)
!---------------------------------------------------------------------
      DO K=2,K_MOUNT
        DO I=1,POINTS
          IF( K.EQ.2 ) THEN
! level 1.5
            PU=PSTAR(I)*BKH(K) + AKH(K)
            PL=PSTAR(I)*BKH(K-1) + AKH(K-1)
            P_EXNER_CENTRE(I,KL)=
     &      P_EXNER_C( PEXNER(I,K),PEXNER(I,K-1),PU,PL,KAPPA)
            PL=PU
            PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
            P_EXNER_CENTRE(I,KU)=
     &             P_EXNER_C( PEXNER(I,K+1),PEXNER(I,K),PU,PL,KAPPA)
            N_SQ(I,KL) = G*(THETA(I,K)-THETA(I,K-1))/(THETA(I,K)*
     &       THETA(I,K-1)*(P_EXNER_CENTRE(I,KL)-P_EXNER_CENTRE(I,KU))*
     &       CPBYG)
! level 2.5
            P_EXNER_CENTRE(I,KL)=P_EXNER_CENTRE(I,KU)
            PL=PU
            PU=PSTAR(I)*BKH(K+2) + AKH(K+2)
            P_EXNER_CENTRE(I,KU)=
     &             P_EXNER_C( PEXNER(I,K+2),PEXNER(I,K+1),PU,PL,KAPPA)
            N_SQ(I,KU) = G*(THETA(I,K+1)-THETA(I,K))/(THETA(I,K+1)*
     &       THETA(I,K)*(P_EXNER_CENTRE(I,KL)-P_EXNER_CENTRE(I,KU))*
     &       CPBYG)
! n_squared averaged interpolating boundary values of level 2
            N_SQAV = ( (PEXNER(I,K)-P_EXNER_CENTRE(I,KL))*N_SQ(I,KU)+
     &            (P_EXNER_CENTRE(I,KL) - PEXNER(I,K+1))*N_SQ(I,KL) )
     &                 / ( PEXNER(I,K) - PEXNER(I,K+1) )

            SPEED = U(I,2)*U(I,2) + V(I,2)*V(I,2)
            IF ( N_SQAV .LE. 0.0 .OR. SPEED .LE. 0.0 ) THEN
              UOVERN =  0.0
              LIFT(I)=0.0
            ELSE
              UOVERN=  SPEED / N_SQAV
              UOVERN =  SQRT( UOVERN )
              LIFT(I) = 1.4*SD_OROG(I) - 0.985*UOVERN
            ENDIF
            IF ( LIFT(I) .LE. 0.0 ) THEN
              K_LIFT(I) = 0
              L_CONT(I) = .FALSE.
            ENDIF

            K_BOT(I) = 1
            DEL_EXNER = PEXNER(I,K-1) - PEXNER(I,K)
            ZH(I) = CPBYG*THETA(I,K-1)*DEL_EXNER

          ENDIF     ! k=2
!------------------------------------------------------------------
!  Substitute in Eqn 7. Orographic height Hm = 1.4*SD_OROG(I)
!------------------------------------------------------------------
          DEL_EXNER = PEXNER(I,K) - PEXNER(I,K+1)
          ZH(I) = ZH(I) + CPBYG*THETA(I,K)*DEL_EXNER
          IF ( L_CONT(I) .AND.
     *      (ZH(I).GT.LIFT(I) .OR. ZH(I).GT.750.) ) THEN
            K_LIFT(I) = K
            K_BOT(I)  = K
            L_CONT(I) = .FALSE.
          ENDIF
!------------------------------------------------------------------
!  Find top of mountains where flow is still impacted upon.
!------------------------------------------------------------------
          IF (  L_CONT2(I) .AND.
     *      (ZH(I).GT.1.4*SD_OROG(I) .OR. ZH(I).GT.750.) ) THEN
            K_TOP(I)   = K
            IF( K_TOP(I).LT.3 )         K_TOP(I)=3
            IF( K_BOT(I).GE.K_TOP(I) )  K_BOT(I)=K_TOP(I)-1
            L_CONT2(I) = .FALSE.
          ENDIF

        ENDDO       ! Points
      ENDDO         ! Levels 2 to K_Mount
!---------------------------------------------------------------------
! 1.2 Calculate average N, U ,V ,RHO over 'surface' integral between
!     K_BOT and K_TOP as diagnosed above
!     N_SQAV is linearised from N_SQ at layer boundaries
!---------------------------------------------------------------------
      DO K=2,K_MOUNT
        DO I=1,POINTS
          IF(  K.GT.K_BOT(I) .AND. K.LE.K_TOP(I) ) THEN

            IF( K.EQ.K_BOT(I)+1 ) THEN
              ZH(I)=0.0
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
            PU=PSTAR(I)*BKH(K) + AKH(K)
            RHO=( AK(K) + BK(K)*PSTAR(I) )/( R*THETA(I,K)*
     &         P_EXNER_C( PEXNER(I,K+1),PEXNER(I,K),PL,PU,KAPPA))

!--------------------------------------------------------------------
!   Average 'surface' quantities with pressure
!--------------------------------------------------------------------
            DELTA_P    =  PSTAR(I)*BKH(K) + AKH(K) - PL
            DELTA_P_SUM(I) = DELTA_P_SUM(I) + DELTA_P
            NSQ_S(I)=NSQ_S(I)+( N_SQAV-NSQ_S(I))*DELTA_P/DELTA_P_SUM(I)
            RHO_S(I)=RHO_S(I)+( RHO - RHO_S(I) )*DELTA_P/DELTA_P_SUM(I)
            U_S(I)  =U_S(I)  +( U(I,K)-U_S(I)  )*DELTA_P/DELTA_P_SUM(I)
            V_S(I)  =V_S(I)  +( V(I,K)-V_S(I)  )*DELTA_P/DELTA_P_SUM(I)

          ENDIF ! H_Jump and L_Cont
        ENDDO   ! Points
!  Rename lower centre array as upper centre ready for next level
        KK=KU
        KU=KL
        KL=KK
      ENDDO     ! Levels 2 to K_Trop

!------------------------------------------------------------------
! 2.0 Begin calulation of surface stress
!-----------------------------------------------------------------
      DO I=1,POINTS

        SPEED = U_S(I)*U_S(I) + V_S(I)*V_S(I)

        IF ( SPEED .LE. 0.0 ) THEN
          S_X_OROG(I) = 0.0
          S_Y_OROG(I) = 0.0
        ELSE
          SPEED = SQRT(SPEED)
! Surf X/Y orog define orography as seen by the surface gravity wave
          S_X_OROG(I)= (U_S(I)*SIGMA_XX(I)+V_S(I)*SIGMA_XY(I)) /SPEED
          S_Y_OROG(I)= (U_S(I)*SIGMA_XY(I)+V_S(I)*SIGMA_YY(I)) /SPEED

          S_STRESS_SQ= S_X_OROG(I)*S_X_OROG(I)+S_Y_OROG(I)*S_Y_OROG(I)
          SPEEDCALC = U_S(I)*S_X_OROG(I) + V_S(I)*S_Y_OROG(I)
          IF ( S_STRESS_SQ .LE. 0.0 ) THEN
            SPEED    = 0.0
          ELSE
            SPEED    = SPEEDCALC / SQRT( S_STRESS_SQ )
! Speed component in dirn. of surface gravity wave
          ENDIF
        ENDIF

!------------------------------------------------------------------
!  2.1  Calculate BRUNT-VAISALA frequency  Eqn 1,2
!------------------------------------------------------------------
        IF(  NSQ_S(I).LE.0.0  .OR.  SPEED.LE.0.0
     &       .OR.  SD_OROG(I).LE.0.0  ) THEN
          S_X_STRESS(I) = 0.0
          S_Y_STRESS(I) = 0.0
          TEST(I)=1.0
!       NO WAVES IF UNSTABLE OR IF NO OROGRAPHY
        ELSE
          N   = SQRT( NSQ_S(I) )

!------------------------------------------------------------------
!  2.2  Impose semi-emperical formula. LINEAR/SQUARE/CUBIC DRAG
!------------------------------------------------------------------
          H_SQ = 1.0
          TEST(I) = N*SD_OROG(I) / SPEED
          IF( BETA1.GT.TEST(I) .AND. TEST(I).GE.ALPHA1 ) THEN
            H_SQ = TEST(I)*R_BETA
          ELSE IF( TEST(I).LT.ALPHA1 ) THEN
            H_SQ = TEST(I)*R_BETA*TEST(I)*R_ALPHA
          ENDIF

!------------------------------------------------------------------
!  2.3  Calculate 'SURFACE' STRESS
!------------------------------------------------------------------
          CALC= KAY*RHO_S(I)*(SPEED**3)*H_SQ/(N*SD_OROG(I)*SD_OROG(I))
          S_X_STRESS(I) = S_X_OROG(I) * CALC
          S_Y_STRESS(I) = S_Y_OROG(I) * CALC

        ENDIF     ! SPEED OR N OR OROG .LE. 0.0 ELSE
      END DO      ! I=POINTS
!   END LOOP POINTS


      RETURN
      END
