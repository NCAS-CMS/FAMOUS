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
! SUBROUTINE GW_SCOR TO DELINEATE VERT SCORER PROFILE FOR GW_LEE/GW_VERT
!
      SUBROUTINE GW_SCOR
     1 (PSTAR,PEXNER,THETA,U,V,LEVELS,START_L,H_JUMP,POINTS,AKH,BKH
     2 ,UNIT_X,UNIT_Y,TRANS,K_LEE,H_LEE,LSQ_LEE,L_LEE)

      IMPLICIT NONE
! Description:   TO CALCULATE THE PROFILE OF THE SCORER PARAMETER
!                TO CALCULATE LEE HEIGHT AND TRANSMITTION COEFFICIENT
!
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE eq 11,30,49.
!
! Current code owner: S.Webster

!
! History:
! Version  Date      Comment
!  3.4   18/10/94   Original Code. J.R.Mitchell
!  4.4   19/09/97   Remove *IF -DEF,CRAY compile options. S.Webster 
!LL  4.5  14/05/98  Add Fujitsu directive to stop potential zero
!LL                 divide caused by optimizer.  RBarnes@ecmwf.int
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
     * LEVELS            !IN  NUMBER OF MODEL LEVELS
     *,START_L           !IN  START LEVEL FOR GWD CALCULATIONS
     *,POINTS            !IN  NUMBER OF POINTS
     *,K_LEE(POINTS,2)   !OUT MODEL LEVEL OF LEE HEIGHT AND TOP

      REAL
     * PSTAR(POINTS)              !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)    !IN    PEXNER
     *,THETA(POINTS,LEVELS)       !IN    THETA FIELD
     *,U(POINTS,LEVELS)           !IN    U FIELD
     *,V(POINTS,LEVELS)           !IN    V FIELD
     *,UNIT_X(POINTS)             !IN X_COMPNT OF UNIT STRESS VECTOR
     *,UNIT_Y(POINTS)             !IN Y_COMPNT OF UNIT STRESS VECTOR
     *,TRANS(POINTS)              !OUT TRANSMITION COEFFICIENT
     *,LSQ_LEE(POINTS,2)          !OUT L1 SQUARED AND L2 SQUARED FOR
     *                            !    LEE WAVE CALCULATIONS
     *,H_LEE(POINTS)              !OUT LEE HEIGHT TO TOP OF LEVEL
!      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
!      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY

      LOGICAL
     * H_JUMP(POINTS)         !IN  TRUE IF HYDROLIC JUMP REGIME
     *,L_LEE(POINTS)          !OUT TRUE IF LEE WAVE/TRANS POINT

! Local parameters
      REAL CPBYG
      PARAMETER(CPBYG=CP/G)
! Local scalers
      REAL
     * DEL_EXNER,DEL_EXNER_KL,DEL_EXNER_KU,  ! EXNER DIFFERENCE
     * PU,PL

      LOGICAL   FLAG

      INTEGER   I,K,M    ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU ! LEVEL COUNTERS IN ROUTINE

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: LEVELS+8  ARRAYS OF FULL FIELD LENGTH
!
      INTEGER
     * K_LEE_LIM(LEVELS-3) ! TOP MODEL LEVEL HEIGHT FOR CALCULATION
     *                     ! OF L_SQ(I,2), FOR EACH POSSIBLE LEE LEVEL
     *,K_LEE_MAX           ! MAXIMUM LEVEL FOR FINDING LEE WAVE HEIGHT
     *,K_TRANS             ! LEVEL OF L_SQ SPLIT FOR CALCULATION OF
     *                     ! TRANSMITION COEFFICIENT

      REAL
     * ALPHA_L,PHASE            ! LINKED LEE WAVE PARAMETER
     *,DELTA_Z                  ! THICKNESS OF LAYER BETWEEN HALF LEVELS
     *,P_EXNER_CENTRE(POINTS,2) ! EXNER PRESSURE AT LAYER CENTRES
     *,N_SQ(POINTS,2)           ! SQUARE OF BRUNT_VAISALA FREQUENCY AT
     *                          ! LAYER BOUNDARIES
     *,N_SQAV                   ! CENTRE LAYER AV. OF BOUYANCY FREQ SQ
     *,ZH(POINTS)               ! TOTAL HEIGHT OF JUMP CALCUALTION
     *,P(LEVELS-3)              ! PRESSURE AT LEVEL TOP BOUNDARIES
     *                          ! FOR STANDARD POINT
     *,P_LIM(LEVELS-3)          ! PRESSURE TOP BOUNDARIES FOR
     *                          ! K_LEE_LIM
     *,UKP1                     ! COMPONENT OF WIND AT LAYER K+1 IN
     *                          ! DIRN. OF SURFACE STRESS
     *,UK                       ! COMPONENT OF WIND AT LAYER K
     *,UKM1                     ! COMPONENT OF WIND AT LAYER K-1
     *,DU_DZ(POINTS,2)          ! DELTAN OVER DELTAU - LAYER BOUNDARY
     *,LSQ_SUM(POINTS,2)        ! L1 SQUARED AND L2 SQUARED SUMMED
     *,LSQ(POINTS,LEVELS-3)     ! L SQUARED AT A LEVEL.
     *,LSQ1, LSQ2, L1, L2       ! L CALCULATIONS FOR TRANSMITTION COEFF
     *,CONST,CALC               ! CALCULATIONS

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
! 1. Parameter calculation. Tan(0.75*PI) = -1.0
!    Note ATAN produces answer in range -PI/2 to PI/2
!--------------------------------------------------------------------
      PHASE=LEE_PHASE*3.14159
      ALPHA_L= TAN(PHASE)
      CONST= ALPHA_L/( PHASE*SQRT( 1+ALPHA_L**2 ) )
      KL=1
      KU=2
!---------------------------------------------------------------------
! Find approximate level of 0.6 tropopause height - maximum lee limit
! Find approx half tropopause height - transmittion coefficient top
!---------------------------------------------------------------------
      FLAG = .TRUE.
      K_LEE_MAX = LEVELS-4
      K_TRANS = K_LEE_MAX
      DO K= 3,LEVELS-4
         IF (FLAG) THEN
           PU=100000.*BKH(K+1) + AKH(K+1)
           IF ( PU .GT. 50000. ) THEN
             K_TRANS = K
           ELSE IF ( PU .LT. 45000. ) THEN
             K_LEE_MAX = K
             FLAG = .FALSE.
           ENDIF
         ENDIF
      ENDDO
!---------------------------------------------------------------------
! Find array of top level scorer calculation limits K_LEE_LIM, for each
! possible lee level by using the same standard point as above.
!---------------------------------------------------------------------
      DO K= 2,LEVELS-3
        P(K)=100000.*BKH(K+1) + AKH(K+1)
        IF ( K .LE. K_LEE_MAX ) THEN
          P_LIM(K) = 2*P(K) - 100000.
          K_LEE_LIM(K) = LEVELS-3
        ENDIF
      ENDDO

      DO M= 2,K_LEE_MAX
        DO K= M+1,LEVELS-3
          IF ( P(K) .LT. P_LIM(M) .AND.
     *       K_LEE_LIM(M).EQ.LEVELS-3  ) THEN
             K_LEE_LIM(M) = K
          ENDIF
        ENDDO
      ENDDO
      K_LEE_LIM(1)=2
!---------------------------------------------------------------------
!  Initialisation
!---------------------------------------------------------------------
      DO I=1,POINTS
        K_LEE(I,1)=0
        K_LEE(I,2)=0
        IF( H_JUMP(I) ) THEN
          L_LEE(I)=.FALSE.
        ELSE
          L_LEE(I)=.TRUE.
        ENDIF
        TRANS(I)=1.0
        ZH(I)=0.0
      ENDDO

!---------------------------------------------------------------------
!  Bottom Level calculation - Calculates L_SQ at level 2
!---------------------------------------------------------------------
      DO I=1,POINTS
        IF( L_LEE(I) ) THEN
          PU=PSTAR(I)*BKH(2) + AKH(2)
          PL=PSTAR(I)*BKH(1) + AKH(1)
!-------------- KU / KL reversed ready for next stage------------
          P_EXNER_CENTRE(I,KU)=
     &       P_EXNER_C( PEXNER(I,2),PEXNER(I,1),PU,PL,KAPPA)
          PL=PU
          PU=PSTAR(I)*BKH(3) + AKH(3)
          P_EXNER_CENTRE(I,KL)=
     &       P_EXNER_C( PEXNER(I,3),PEXNER(I,2),PU,PL,KAPPA)
!-------------------------------------------------------------
          N_SQ(I,KL) = G*(THETA(I,2)-THETA(I,1))/(THETA(I,2)*
     &      THETA(I,1)*(P_EXNER_CENTRE(I,KU)-P_EXNER_CENTRE(I,KL))
     &      *CPBYG)
          PL=PU
          PU=PSTAR(I)*BKH(4) + AKH(4)
          P_EXNER_CENTRE(I,KU)=
     &        P_EXNER_C( PEXNER(I,4),PEXNER(I,3),PU,PL,KAPPA)
          N_SQ(I,KU) = G*(THETA(I,3)-THETA(I,2))/(THETA(I,3)*
     &      THETA(I,2)*(P_EXNER_CENTRE(I,KL)-P_EXNER_CENTRE(I,KU))
     &      *CPBYG)
!--------------------------------------------------------------------
!  N_SQAV at layer centre, interpolated from N_SQ(KL/KU) at boundaries
!--------------------------------------------------------------------
          N_SQAV = ( (PEXNER(I,2)-P_EXNER_CENTRE(I,KL))*N_SQ(I,KU)
     &         + (P_EXNER_CENTRE(I,KL) - PEXNER(I,3))*N_SQ(I,KL) )
     &         / ( PEXNER(I,2) - PEXNER(I,3) )
!--------------------------------------------------------------------
!   Note U is component parallel to stress vector
!--------------------------------------------------------------------
          UKP1   = U(I,3)*UNIT_X(I) + V(I,3)*UNIT_Y(I)
          UK     = U(I,2)*UNIT_X(I) + V(I,2)*UNIT_Y(I)
          UKM1   = U(I,1)*UNIT_X(I) + V(I,1)*UNIT_Y(I)
          IF( N_SQAV .LE. 0.0 .OR. UKP1.LE.0.3
     *         .OR. UK.LE.0.2 .OR. UKM1.LE.0.1 ) THEN
            L_LEE(I)=.FALSE.
          ELSE
            DEL_EXNER_KU = PEXNER(I,3) - P_EXNER_CENTRE(I,KU)
            DEL_EXNER_KL = P_EXNER_CENTRE(I,KL) - PEXNER(I,3)
            DELTA_Z= CPBYG*(THETA(I,2)*DEL_EXNER_KL +
     *                        THETA(I,3)*DEL_EXNER_KU)
            DU_DZ(I,KU)= ( UKP1-UK )/ DELTA_Z
            PU=PSTAR(I)*BKH(2) + AKH(2)
            PL=PSTAR(I)*BKH(1) + AKH(1)
            CALC=
     &         P_EXNER_C( PEXNER(I,2),PEXNER(I,1),PU,PL,KAPPA)
            DEL_EXNER_KU = PEXNER(I,2) - P_EXNER_CENTRE(I,KL)
            DEL_EXNER_KL = CALC - PEXNER(I,2)
            DELTA_Z= CPBYG*(THETA(I,1)*DEL_EXNER_KL +
     *                      THETA(I,2)*DEL_EXNER_KU )

            DU_DZ(I,KL)= ( UK-UKM1 )/ DELTA_Z
            DEL_EXNER_KU = PEXNER(I,2) - PEXNER(I,3)
            DEL_EXNER_KL = PEXNER(I,1) - PEXNER(I,2)
            ZH(I) = CPBYG*(THETA(I,1)*DEL_EXNER_KL +
     *                     THETA(I,2)*DEL_EXNER_KU )
            LSQ(I,2)=N_SQAV /UK**2 - ( DU_DZ(I,KU)-DU_DZ(I,KL) )
     *                          / (UK*CPBYG*THETA(I,2)*DEL_EXNER_KU)
            LSQ_SUM(I,1)= 0.0
            LSQ_SUM(I,2)= LSQ(I,2)
          ENDIF
        ENDIF   ! lee point
      ENDDO     ! Points

!---------------------------------------------------------------------
!  Loop Levels
!---------------------------------------------------------------------
      DO M= 2,K_LEE_MAX
        DO K= K_LEE_LIM(M-1), K_LEE_LIM(M)

          IF( K.NE.K_LEE_LIM(M-1) ) THEN
            KK=KL
            KL=KU
            KU=KK
          ENDIF

! Optimizer precomputes values which can cause zero divide, so:-
!OCL NOPREEX,NOEVAL
          DO I=1,POINTS
            IF(   L_LEE(I) .AND. ( K.NE.K_LEE_LIM(M-1) )
     &         .AND. ( K_LEE(I,1).EQ.0 .OR. M.LE.K_TRANS ) ) THEN
! next level stage
              PU=PSTAR(I)*BKH(K+2) + AKH(K+2)
              PL=PSTAR(I)*BKH(K+1) + AKH(K+1)
              P_EXNER_CENTRE(I,KU)=
     &            P_EXNER_C( PEXNER(I,K+2),PEXNER(I,K+1),PU,PL,KAPPA)
              N_SQ(I,KU) = G*(THETA(I,K+1)-THETA(I,K))/(THETA(I,K+1)*
     &        THETA(I,K)*(P_EXNER_CENTRE(I,KL)-P_EXNER_CENTRE(I,KU))*
     &          CPBYG)
              N_SQAV = ( (PEXNER(I,K)-P_EXNER_CENTRE(I,KL))*N_SQ(I,KU)
     &          + (P_EXNER_CENTRE(I,KL) - PEXNER(I,K+1))*N_SQ(I,KL) )
     &                  / ( PEXNER(I,K) - PEXNER(I,K+1) )
!--------------------------------------------------------------------
!   Note U is component parallel to stress vector
!--------------------------------------------------------------------
              UKP1 = U(I,K+1)*UNIT_X(I) + V(I,K+1)*UNIT_Y(I)
              UK   = U(I,K)  *UNIT_X(I) + V(I,K)  *UNIT_Y(I)
              UKM1 = U(I,K-1)*UNIT_X(I) + V(I,K-1)*UNIT_Y(I)
              IF( N_SQAV .LE. 0.0 .OR. UKP1 .LE. 0.3) THEN
                L_LEE(I)=.FALSE.
              ELSE
                DEL_EXNER_KU = PEXNER(I,K+1) - P_EXNER_CENTRE(I,KU)
                DEL_EXNER_KL = P_EXNER_CENTRE(I,KL) - PEXNER(I,K+1)
                DELTA_Z= CPBYG*(THETA(I,K)*DEL_EXNER_KL +
     *                          THETA(I,K+1)*DEL_EXNER_KU)
                DU_DZ(I,KU)= ( UKP1-UK )/ DELTA_Z
                DEL_EXNER = PEXNER(I,K) - PEXNER(I,K+1)
                LSQ(I,K)=N_SQAV/UK**2 - ( DU_DZ(I,KU)-DU_DZ(I,KL) )
     *                         / (UK*CPBYG*THETA(I,K)*DEL_EXNER)
                LSQ_SUM(I,2)=LSQ_SUM(I,2)+LSQ(I,K)
              ENDIF
            ENDIF  ! ...K.NE.K_Lee_Lim(M-1)..

            IF(   L_LEE(I)      .AND.
     &          ( K_LEE(I,1).EQ.0 .OR. M.LE.K_TRANS )  ) THEN

              IF( K .EQ. K_LEE_LIM(M) ) THEN

                DEL_EXNER = PEXNER(I,M) - PEXNER(I,M+1)
                ZH(I) = ZH(I) + CPBYG*THETA(I,M)*DEL_EXNER

                LSQ_SUM(I,1)=LSQ_SUM(I,1)+LSQ(I,M)
                LSQ_SUM(I,2)=LSQ_SUM(I,2)-LSQ(I,M)
                IF( LSQ_SUM(I,2).LT.0.0 ) LSQ_SUM(I,2)=0.0
                CALC= LSQ_SUM(I,1)/(M-1) - LSQ_SUM(I,2)/(K-M)
                IF( CALC .GT. 0.0 .AND. K_LEE(I,1).EQ.0
     &              .AND. M.GE.START_L ) THEN
                  CALC=CONST*ZH(I)*SQRT(CALC)
                  IF( CALC.LT.-1.0 ) THEN
                    K_LEE(I,1)=M
                    K_LEE(I,2)=K
                    H_LEE(I)=ZH(I)
                    LSQ_LEE(I,1)=LSQ_SUM(I,1)/(M-1)
                    LSQ_LEE(I,2)=LSQ_SUM(I,2)/(K-M)
                  ENDIF
                ENDIF
              ENDIF

              IF( M .EQ. K_TRANS ) THEN
                IF( LSQ_SUM(I,1).GT.0.0 .AND. LSQ_SUM(I,2).GT.0.0 )THEN
                  LSQ1=LSQ_SUM(I,1)/(M-1)
                  LSQ2=LSQ_SUM(I,2)/(K-M)
                  L1=SQRT(LSQ1)
                  L2=SQRT(LSQ2)
                  CALC=  LSQ1+LSQ2+(LSQ1-LSQ2)*COS(2.*L1*ZH(I))
                  IF( CALC.EQ.0.0) THEN
                    TRANS(I)=1.0
                  ELSE
                    TRANS(I)=2*L1*L2/ CALC
                  ENDIF
                  IF( TRANS(I).GT.1.0 ) THEN
                    TRANS(I)=1.0
                  ENDIF
                ENDIF
              ENDIF

              IF( K_LEE(I,1).EQ.0 .AND. M.EQ.K_LEE_MAX ) THEN
                L_LEE(I)= .FALSE.
              ENDIF

            ENDIF    ! K_Lee or K_trans = 0
          ENDDO      ! Points
        ENDDO        ! K-Loop
      ENDDO          ! M-Loop

! FIND CRITICAL LEVELS?
      RETURN
      END

