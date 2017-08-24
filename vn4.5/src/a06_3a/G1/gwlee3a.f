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
!  SUBROUTINE GW_LEE TO CALCULATE ADDITIONAL TRAPPED LEE WAVE DRAG
!
      SUBROUTINE GW_LEE
     1  (PSTAR,START_L,LEVELS,POINTS,AKH,BKH,DELTA_AK,DELTA_BK
     2  ,U_S,V_S,RHO_S,L_LEE,LSQ,H_LEE,K_LEE,KAY_LEE                    
     3  ,SIGMA_XX,SIGMA_XY,SIGMA_YY,DU_DT,DV_DT
! Diagnostics
     4  ,STRESS_UD,POINTS_STRESS_UD,STRESS_UD_ON
     5  ,STRESS_VD,POINTS_STRESS_VD,STRESS_VD_ON
     6  ,DU_DT_LEE,POINTS_DU_DT_LEE,DU_DT_LEE_ON
     7  ,DV_DT_LEE,POINTS_DV_DT_LEE,DV_DT_LEE_ON )

      IMPLICIT NONE
! Description:
!             TO CALCULATE DRAG/STRESS PROFILE DUE TO SUBGRID-SCALE
!             OROGRAPHIC GRAVITY WAVES FOR TRAPPED LEE WAVE STATES
!
! Method: UNIFIED MODEL DOCUMENTATION PAPER NO. ?
!         THE EQUATIONS USED ARE (???) TO (???)
!
! Current code owner: S.Webster
!
! History:
! Version  Date      Comment
!  3.4   18/10/94   Original Code. J.R.Mitchell
!  4.3    7/03/97   KAY_LEE passed in from namelist. S.Webster
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
! Subroutine arguements;
      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,START_L             !IN    START LEVEL FOR WAVE-BREAKING TEST
     *,POINTS              !IN    NUMBER OF POINTS
     *,POINTS_STRESS_UD    !IN    ) No of land points in diagnostic
     *,POINTS_STRESS_VD    !IN    ) arrays for GW stress - u and v
     *,POINTS_DU_DT_LEE    !IN    ) No of land points in diagnostic
     *,POINTS_DV_DT_LEE    !IN    ) arrays for GW lee - du and dv
     *,K_LEE(POINTS,2)     !IN    1. LEVEL OF LEE HEIGHT INTERFACE
     *                     !      2. LEVEL OF LEE HEIGHT TOP

      LOGICAL
     * L_LEE(POINTS)       !IN    TRUE IF POINT IS TO BE LINEARIZED
     *,STRESS_UD_ON        !IN    U STRESS DIAGNOSTIC SWITCH
     *,STRESS_VD_ON        !IN    V STRESS DIAGNOSTIC SWITCH
     *,DU_DT_LEE_ON        !IN    U-ACCELN LEE WAVE COMPT SWITCH
     *,DV_DT_LEE_ON        !IN    V-ACCELN LEE WAVE COMPT SWITCH

      REAL
     * PSTAR(POINTS)               !IN    PSTAR FIELD
     *,U_S(POINTS)                 !IN    'SURFACE' U FIELD
     *,V_S(POINTS)                 !IN    'SURFACE' V FIELD
     *,RHO_S(POINTS)               !IN    'SURFACE' LAYER DENSITY
!      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
!      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY
     *,DELTA_AK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,DELTA_BK(LEVELS)       !IN    DIFFERENCE ACROSS LAYER
     *,DU_DT(POINTS,LEVELS)   !IN/OUT   U-ACCELERATION
     *,DV_DT(POINTS,LEVELS)   !IN/OUT   V-ACCELERATION
     *,H_LEE(POINTS)          !IN    LEE HEIGHT
     *,LSQ(POINTS,2)          !IN    SCORER PARAMETER ABOVE AND
     *                        !IN    BELOW LEE HEIGHT
     *,KAY_LEE                !IN    TRAPPED LEE WAVE CONSTANT
     *,SIGMA_XX(POINTS)       !IN    OROGRAPHIC GRADIENT IN X-DIRN.
     *,SIGMA_XY(POINTS)       !IN    OROGRAPHIC GRADIENT IN X/Y-DIRN
     *,SIGMA_YY(POINTS)       !IN    OROGRAPHIC GRADIENT IN Y-DIRN.
     *,STRESS_UD(POINTS_STRESS_UD,LEVELS+1) !U STRESS  DIAGNOSTIC
     *,STRESS_VD(POINTS_STRESS_VD,LEVELS+1) !V STRESS  DIAGNOSTIC
     *,DU_DT_LEE(POINTS_DU_DT_LEE,LEVELS)   !U-ACCELN LEE WAVE
     *,DV_DT_LEE(POINTS_DV_DT_LEE,LEVELS)   !V-ACCELN LEE WAVE

! Local parameters
! Local scalers
      REAL
     * DELTA_P,DELTA_P1,DELTA_P2 ! DIFFERENCES IN PRESSURE
     *,PHASE                 ! PHASE ACROSS LEE INTERFACE (TUNABLE)
     *,ALPHA_L               ! TANGENT OF PHASE
     *,AL_SQ                 ! ALPHA_L SQUARED
     *,AL_SQ_PLUS1           ! ALPHA_L SQUARED PLUS 1
     *,GAMMA                 ! PHASE OVER ALPHA
     *,CONST                 ! CALCULATION CONSTANT OVER POINTS
     *,M2_SQ                 ! M_2 (UPPER VERTICAL WAVE NO) SQUARED
     *,K_SQ                  ! HORIZONTAL WAVE NO. SQUARED
     *,CALC1,CALC2           ! CALCULATIONS
     *,SPEED_SQ,SPEED        ! 'SURFACE' WIND SPEED CALCULATIONS
     *,COS,SIN               ! DIRECTION OF SURFACE WIND
     *,DELTA_PHI             ! ANGULAR WIDTH COEFFICIENT
     *,DELTA_AK_SUM ! DELTA_AK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,DELTA_BK_SUM ! DELTA_BK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,PU,PL
     *,KayLee_x_DeltaPhi     
      PARAMETER ( DELTA_PHI = 1.0 )      

      INTEGER   I,K    ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU,KT ! LEVEL COUNTERS IN ROUTINE

! Local dynamic arrays
! LOCAL WORKSPACE ARRAYS: 12  ARRAYS OF FULL FIELD LENGTH
!

      REAL
     * X_STRESS(POINTS,2)    ! X_STRESSES (LAYER BOUNDARIES)
     *,Y_STRESS(POINTS,2)    ! Y_STRESSES (LAYER BOUNDARIES)
     *,DP_X_STRESS(POINTS,2) ! STRESS X_GRADIENTS
     *                       ! 1.BELOW  2.ABOVE LEE HEIGHT INTERFACE
     *,DP_Y_STRESS(POINTS,2) ! STRESS Y_GRADIENTS
     *,X_LEE_STRESS(POINTS)  ! LEE WAVE X-STRESS AT SURFACE
     *,Y_LEE_STRESS(POINTS)  ! LEE WAVE Y-STRESS AT SURFACE
     *,RATIO(POINTS)         ! RATIO OF STRESS AT LEE HEIGHT TO
     *                       ! TRAPPED LEE WAVE STRESS AT SURF

! Function and subroutine calls: NONE

!-----------------------------------------------------------------
!   1. PRELIMINARIES
!-----------------------------------------------------------------

CFPP$ NOCONCUR L
!      TREAT LAYERS BELOW AND INCLUDING START_L AS ONE LAYER
        DELTA_AK_SUM = 0.0
        DELTA_BK_SUM = 0.0
      DO K=2,START_L
        DELTA_AK_SUM = DELTA_AK_SUM + DELTA_AK(K)
        DELTA_BK_SUM = DELTA_BK_SUM + DELTA_BK(K)
      END DO
CFPP$ CONCUR

      KayLee_x_DeltaPhi = KAY_LEE * DELTA_PHI 
      PHASE= LEE_PHASE*3.14159

      ALPHA_L    = TAN(PHASE)
      AL_SQ      = ALPHA_L*ALPHA_L
      AL_SQ_PLUS1= AL_SQ+1.0
      GAMMA      = PHASE/ALPHA_L
      CONST      = (GAMMA-1.) /(2*AL_SQ*GAMMA**3)

      KL=1
      KU=2
      KT=3
!---------------------------------------------------------------------
!   2. CALCULATE TRAPPED LEE WAVE SURFACE STRESS, UW(SURF)
!      CALCULATE INTERFACE STRESS RATIO  UW(Z=H)/UW(SURF)
!      CALCLATED STRESS GRADIENTS USING LINEAR INTERPOLATION
!      BETWEEN UW(SURF)/UW(Z=H) AND UW(Z=H)/0.0(Z=2H)
!---------------------------------------------------------------------

      DO I=1,POINTS
        IF( L_LEE(I) ) THEN

! Calculate Ratio
          M2_SQ= (LSQ(I,1)-LSQ(I,2)) / AL_SQ_PLUS1
          RATIO(I)= AL_SQ*LSQ(I,2) /
     *      ( (LSQ(I,1)+AL_SQ*LSQ(I,2))*(H_LEE(I)*SQRT(M2_SQ)+1.) )
! Calculate surface lee wave stress
          K_SQ = ( LSQ(I,1)+AL_SQ*LSQ(I,2) ) / AL_SQ_PLUS1
          CALC1= CONST*(H_LEE(I)**3)*(K_SQ**0.75)
          SPEED_SQ = U_S(I)*U_S(I)+V_S(I)*V_S(I)
          SPEED    = SQRT(SPEED_SQ)
          IF( SPEED.LE.0.0 ) THEN
             L_LEE(I)=.FALSE.
          ELSE
            COS      = U_S(I)/SPEED
            SIN      = V_S(I)/SPEED
            CALC2= 0.75*(  SIGMA_XX(I)*(4.*COS*COS - 1.)
     *                  +  SIGMA_XY(I)*(8.*COS*SIN)
     *                  +  SIGMA_YY(I)*(4.*SIN*SIN - 1.)  )
! If CALC2 is negavive then stress opposses surface wind
! This can occur in error if surface wind is at acute angle to a ridge.
            IF ( CALC2.LT.0.0 ) THEN
              CALC2 = 0.0
            ENDIF
            Y_LEE_STRESS(I)=
     *          SPEED*RHO_S(I)*KayLee_x_DeltaPhi*CALC2/CALC1          
            X_LEE_STRESS(I)= Y_LEE_STRESS(I)*U_S(I)
            Y_LEE_STRESS(I)= Y_LEE_STRESS(I)*V_S(I)
! Calculate stress gradients
            PU=PSTAR(I)*BKH(K_LEE(I,1)+1) + AKH(K_LEE(I,1)+1)
            PL=PSTAR(I)*BKH(START_L) + AKH(START_L)
            DELTA_P1= PL - PU
            PL=PSTAR(I)*BKH(K_LEE(I,2)+1) + AKH(K_LEE(I,2)+1)
            DELTA_P2= PU - PL
            DP_X_STRESS(I,1) = (1.-RATIO(I))*X_LEE_STRESS(I)/
     &                                   ( DELTA_P1 )
            DP_Y_STRESS(I,1) = (1.-RATIO(I))*Y_LEE_STRESS(I)/
     &                                   ( DELTA_P1 )
            DP_X_STRESS(I,2) = (RATIO(I))*X_LEE_STRESS(I)/
     &                                   ( DELTA_P2 )
            DP_Y_STRESS(I,2) = (RATIO(I))*Y_LEE_STRESS(I)/
     &                                   ( DELTA_P2 )

          ENDIF ! Speed<=0

        ENDIF ! if L_LEE(I)

      END DO  ! Points

      IF( STRESS_UD_ON ) THEN
        DO I=1,POINTS
          IF( L_LEE(I) ) THEN
            X_STRESS(I,KL) = X_LEE_STRESS(I)
            STRESS_UD(I,START_L) = STRESS_UD(I,START_L)+X_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

      IF( STRESS_VD_ON ) THEN
        DO I=1,POINTS
          IF( L_LEE(I) ) THEN
            Y_STRESS(I,KL) = Y_LEE_STRESS(I)
            STRESS_VD(I,START_L) = STRESS_VD(I,START_L)+Y_STRESS(I,KL)
          ENDIF
        END DO
      ENDIF

!------------------------------------------------------------------
!    3 LOOP LEVELS
!------------------------------------------------------------------
      DO K=START_L,LEVELS


        DO I=1,POINTS
          IF( L_LEE(I) .AND. K.LE.K_LEE(I,1) ) THEN

            IF( K .EQ. START_L ) THEN
              DELTA_P =(DELTA_AK(START_L)+DELTA_BK(START_L)*PSTAR(I))
     &               /( DELTA_AK_SUM     +DELTA_BK_SUM*PSTAR(I) )
              DU_DT(I,START_L) = DU_DT(I,START_L) -
     &                           G*DP_X_STRESS(I,1)*DELTA_P
              DV_DT(I,START_L) = DV_DT(I,START_L) -
     &                           G*DP_Y_STRESS(I,1)*DELTA_P
            ELSE
              DU_DT(I,K) = DU_DT(I,K) - G*DP_X_STRESS(I,1)
              DV_DT(I,K) = DV_DT(I,K) - G*DP_Y_STRESS(I,1)
            ENDIF

          ELSE IF ( L_LEE(I) .AND. K.LE.K_LEE(I,2) ) THEN

            DU_DT(I,K) = DU_DT(I,K) - G*DP_X_STRESS(I,2)
            DV_DT(I,K) = DV_DT(I,K) - G*DP_Y_STRESS(I,2)

          ENDIF  ! if L_LEE .and. K<=K_LEE(I,1) .else. K<=K_LEE(I,2)
        END DO ! points

! Diagnostics
        IF( DU_DT_LEE_ON ) THEN
          DO I=1,POINTS
            IF( L_LEE(I) .AND. K.LE.K_LEE(I,1) ) THEN
              DU_DT_LEE(I,K) = -G*DP_X_STRESS(I,1)
              IF( K.EQ.START_L ) THEN
               DELTA_P =(DELTA_AK(START_L)+DELTA_BK(START_L)*PSTAR(I))
     &                /( DELTA_AK_SUM     +DELTA_BK_SUM*PSTAR(I) )
                DU_DT_LEE(I,K)=DU_DT_LEE(I,K)*DELTA_P
              ENDIF
            ELSE IF ( L_LEE(I) .AND. K.LE.K_LEE(I,2) ) THEN
              DU_DT_LEE(I,K) = -G*DP_X_STRESS(I,2)
            ENDIF
          END DO
        ENDIF

        IF( DV_DT_LEE_ON ) THEN
          DO I=1,POINTS
            IF( L_LEE(I) .AND. K.LE.K_LEE(I,1) ) THEN
              DV_DT_LEE(I,K) = -G*DP_Y_STRESS(I,1)
              IF( K.EQ.START_L ) THEN
               DELTA_P =(DELTA_AK(START_L)+DELTA_BK(START_L)*PSTAR(I))
     &                /( DELTA_AK_SUM     +DELTA_BK_SUM*PSTAR(I) )
                DV_DT_LEE(I,K)=DV_DT_LEE(I,K)*DELTA_P
              ENDIF
            ELSE IF ( L_LEE(I) .AND. K.LE.K_LEE(I,2) ) THEN
              DV_DT_LEE(I,K) = -G*DP_Y_STRESS(I,2)
            ENDIF
          END DO
        ENDIF

! Calculate stress at upper boundary. NB DELTA_P is -ve
        IF( STRESS_UD_ON ) THEN
          DO I=1,POINTS
            IF( L_LEE(I) .AND. K.LE.K_LEE(I,1) ) THEN
              DELTA_P = DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
              X_STRESS(I,KU)=X_STRESS(I,KL)+DP_X_STRESS(I,1)*DELTA_P
              STRESS_UD(I,K+1) = STRESS_UD(I,K+1) + X_STRESS(I,KU)
            ELSE IF ( L_LEE(I) .AND. K.LE.K_LEE(I,2) ) THEN
              DELTA_P = DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
              X_STRESS(I,KU)=X_STRESS(I,KL)+DP_X_STRESS(I,2)*DELTA_P
              STRESS_UD(I,K+1) = STRESS_UD(I,K+1) + X_STRESS(I,KU)
            ENDIF
          END DO
        ENDIF

        IF( STRESS_VD_ON ) THEN
          DO I=1,POINTS
            IF( L_LEE(I) .AND. K.LE.K_LEE(I,1) ) THEN
              DELTA_P = DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
              Y_STRESS(I,KU)=Y_STRESS(I,KL)+DP_Y_STRESS(I,1)*DELTA_P
              STRESS_VD(I,K+1) = STRESS_VD(I,K+1) + Y_STRESS(I,KU)
            ELSE IF ( L_LEE(I) .AND. K.LE.K_LEE(I,2) ) THEN
              DELTA_P = DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)
              Y_STRESS(I,KU)=Y_STRESS(I,KL)+DP_Y_STRESS(I,2)*DELTA_P
              STRESS_VD(I,K+1) = STRESS_VD(I,K+1) + Y_STRESS(I,KU)
            ENDIF
          END DO
        ENDIF

! Swap storage for lower and upper layers
        KK=KL
        KL=KU
        KU=KK

      END DO  ! Levels

      RETURN
      END

