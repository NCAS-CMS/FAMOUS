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
CLL  SUBROUTINE GW_RICH------------------------------------------
CLL
CLL  PURPOSE:   TO CALCULATE STRESS PROFILE DUE TO SUBGRID-SCALE
CLL             OROGRAPHIC GRAVITY WAVES AND HENCE DRAG ON MEAN FLOW.
CLL             THE WAVES PROPOGATE VERTICALLY WITH STRESS INDEPENDENT
CLL             HEIGHT UNLESS A CRITICAL LEVEL OR WAVE BREAKING IS
CLL             DIAGNOSED. A MINIMUM RICHARDSON NUMBER CRITERION
CLL             DETERMINES WHERE THIS OCCURS AND THE WAVE AMPLITUDE
CLL             AND STRESS REDUCED SO THAT THE WAVES ARE MAINTAINED
CLL             AT MARGINAL STABILITY.
CLL  SUITABLE FOR SINGLE COLUMN USE
CLL
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL  ORIGINAL VERSION FOR CRAY Y-MP
CLL  WRITTEN BY C. WILSON
CLL  FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   3.3   25/10/93  Removal of DIAG06 directive. New arguments to
CLL                   dimension diagnostic arrays. D. Robinson.
CLL   3.4   24/09/94  Test correct diagnostic switch for GW Stress
CLL                   v component. D. Robinson.
CLL
CLL   4.4   19/09/97  Remove *IF -DEF,CRAY compile options. S.Webster 
CLL                                                                     
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 1, DATED 12/09/89
CLL
CLL  logical components covered: P22
CLL
CLL  SYSTEM TASK: PART OF P22
CLL
CLL  DOCUMENTATION:        THE EQUATIONS USED ARE (2.6) TO (2.10)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO  22
CLL                        C. A. WILSON AND R. SWINBANK
CLL                        VERSION 1,DATED 15/12/89
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE GW_RICH
     1  (PSTAR,PEXNER,THETA,U,V,S_STRESS,START_L,LEVELS,POINTS,
     2   AKH,BKH,DELTA_AK,DELTA_BK,KAY,SIN_A,COS_A,
     3   DU_DT,DV_DT,
     4   STRESS_UD_LAND,LAND_POINTS_UD,STRESS_UD_ON,
     5   STRESS_VD_LAND,LAND_POINTS_VD,STRESS_VD_ON)

      IMPLICIT NONE

      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,START_L             !IN    START LEVEL FOR WAVE-BREAKING TEST
     *,POINTS              !IN    NUMBER OF POINTS
     *,LAND_POINTS_UD      !IN    ) No of land points in diagnostic
     *,LAND_POINTS_VD      !IN    ) arrays for GW stress - u and v.

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,THETA(POINTS,LEVELS)             !IN    THETA FIELD
     *,U(POINTS,LEVELS)                 !IN    U FIELD
     *,V(POINTS,LEVELS)                 !IN    V FIELD
     *,S_STRESS(POINTS)                 !IN    'SURFACE' STRESS
C      AKH,BKH  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*-LAYER EDGES,
C      DELTA_AK,DELTA_BK  DEFINE PRESSURE DIFFERENCES ACROSS LAYERS
      REAL
     * AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUMDARY
     *,DELTA_AK (LEVELS)      !IN    DIFFERENCE ACROSS LAYER
     *,DELTA_BK (LEVELS)      !IN    DIFFERENCE ACROSS LAYER
     *,KAY                    !IN    stress constant (m-1)
     *,SIN_A(POINTS)          !IN    SIN (STRESS DIRECTION FROM NORTH)
     *,COS_A(POINTS)          !IN    COS (STRESS DIRECTION FROM NORTH)
     *,DU_DT(POINTS,LEVELS)   !OUT   U-ACCELERATION
     *,DV_DT(POINTS,LEVELS)   !OUT   V-ACCELERATION
     *,STRESS_UD_LAND(LAND_POINTS_UD,LEVELS+1) !U STRESS DIAGNOSTIC
     *,STRESS_VD_LAND(LAND_POINTS_VD,LEVELS+1) !V STRESS DIAGNOSTIC

      LOGICAL
     +   STRESS_UD_ON  ! )  Diagnostic switches for GW stress -
     +  ,STRESS_VD_ON  ! )  u and v  (Item Nos 201 and 202)

C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS:
C  10 REAL ARRAYS OF FULL FIELD LENGTH REQUIRED
C

      REAL
     * DZ(POINTS,3)         ! HEIGHT DIFFERENCES IN EACH HALF LAYER
     *,T(POINTS,2)          ! TEMPERATURES (LEVELS)
     *,SPEED(POINTS,2)      ! WIND SPEEDS (LEVELS)
     *,STRESS(POINTS,2)     ! STRESSES (LAYER BOUNDARIES)
     *,DRAG(POINTS)         ! DRAG EXERTED ON LAYER(IN DIRECTION OF
     *                      ! SURFACE STRESS

C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
C     NONE
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED = POINTS
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C
C   DEFINE LOCAL VARIABLES
C   LOCAL VARIABLES:
C
      REAL
     * RHO                  ! DENSITY AT LAYER BOUNDARY
     *,TB                   ! TEMPERATURE AT LAYER BOUNDARY
     *,SPEEDB               ! WIND SPEEDS AT LAYER BOUNDARY
     *,DZB                  ! HEIGHT DIFFERENCE ACROSS LAYER BOUNDARY
     *,N                    ! BRUNT_VAISALA FREQUENCY
     *,N_SQ                 ! SQUARE OF BRUNT_VAISALA FREQUENCY
     *,AMPL                 ! WAVE-AMPLITUDE
     *,DV_DZ                ! MAGNITUDE OF WIND SHEAR
     *,RI                   ! MINIMUM RICHARDSON NUMBER
     *,EPSILON              ! MAX WAVE-AMPLITUDE*N**2/WIND SPEED
     *,DELTA_P              ! DIFFERENCE IN PRESSURE ACROSS LAYER
      REAL
     * DELTA_AK_SUM ! DELTA_AK SUMMED OVER LOWEST LAYERS UP TO START_L
     *,DELTA_BK_SUM ! DELTA_BK SUMMED OVER LOWEST LAYERS UP TO START_L
      INTEGER   I,K    ! LOOP COUNTER IN ROUTINE
      INTEGER   KK,KL,KU,KT ! LEVEL COUNTERS IN ROUTINE
C
C INCLUDE PHYSICAL CONSTANTS
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


C LOCAL CONSTANTS
      REAL VAR_MAX,RIC
      PARAMETER(
     * VAR_MAX = 160000.    ! Maximum variance of orography (m**2)
     *)
      PARAMETER(
     * RIC = 2.5E-1     ! Critical Richardson number for wave-breaking
     *)

      REAL CPBYG
      PARAMETER(CPBYG=CP/G)

      REAL
     &    PU,PL,P_EXNER_CENTRE
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


C-------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL   1. START LEVEL  PRELIMINARIES
C------------------------------------------

CFPP$ NOCONCUR L
C      TREAT LAYERS BELOW AND INCLUDING START_L AS ONE LAYER
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

        SPEED(I,KL)  = U(I,START_L)*SIN_A(I) + V(I,START_L)*COS_A(I)
        STRESS(I,KL) = S_STRESS(I)

        PU=PSTAR(I)*BKH(START_L+1) + AKH(START_L+1)
        PL=PSTAR(I)*BKH(START_L) + AKH(START_L)
        P_EXNER_CENTRE=
     &    P_EXNER_C( PEXNER(I,START_L+1),PEXNER(I,START_L),PU,PL,KAPPA)

        DZ(I,KL)     = (P_EXNER_CENTRE - PEXNER(I,START_L+1))
     *                 *THETA(I,START_L)*CPBYG
        T(I,KL)      = P_EXNER_CENTRE*THETA(I,START_L)

      END DO

      IF (STRESS_UD_ON) THEN
        DO I=1,POINTS
          STRESS_UD_LAND(I,START_L) = STRESS(I,KL)*SIN_A(I)
        ENDDO
      ENDIF
      IF (STRESS_VD_ON) THEN
        DO I=1,POINTS
          STRESS_VD_LAND(I,START_L) = STRESS(I,KL)*COS_A(I)
        ENDDO
      ENDIF

C------------------------------------------------------------------
CL    2 LOOP LEVELS
C------------------------------------------------------------------

      DO K=START_L+1,LEVELS


        DO I=1,POINTS
          STRESS(I,KU) = STRESS(I,KL)
          IF(STRESS(I,KL) .GT. 0.0) THEN
            SPEED(I,KU)  = U(I,K)*SIN_A(I) + V(I,K)*COS_A(I)

            PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
            PL=PSTAR(I)*BKH(K) + AKH(K)
            P_EXNER_CENTRE=
     &                P_EXNER_C( PEXNER(I,K+1),PEXNER(I,K),PU,PL,KAPPA)

C lower half height of upper layer
            DZ(I,KU)     = (PEXNER(I,K) - P_EXNER_CENTRE)*THETA(I,K)
     *                     *CPBYG
C upper half height of upper layer
            DZ(I,KT)     = (P_EXNER_CENTRE - PEXNER(I,K+1))*THETA(I,K)
     *                     *CPBYG
C model level height difference
            DZB          = DZ(I,KU) + DZ(I,KL)
            SPEEDB       = (DZ(I,KU)*SPEED(I,KL)+DZ(I,KL)*SPEED(I,KU))/
     *                     DZB

C------------------------------------------------------------------
CL          2.1 TEST FOR CRITICAL LEVEL   V < or = 0
C------------------------------------------------------------------

            IF( SPEEDB .LE. 0.0) THEN
              STRESS(I,KU) = 0.0
            ELSE
              T(I,KU)  = P_EXNER_CENTRE*THETA(I,K)
              TB       = (DZ(I,KU)*T(I,KL) + DZ(I,KL)*T(I,KU))/DZB
              RHO      = ( AKH(K) + BKH(K)*PSTAR(I) )/(R*TB)

C------------------------------------------------------------------
CL            2.2 CALCULATE BRUNT-VAISALA FREQUENCY
C------------------------------------------------------------------

              N_SQ = G*( THETA(I,K) - THETA(I,K-1) )*PEXNER(I,K)/
     *                   ( TB*DZB )

              IF( N_SQ .LE.0.0 ) THEN
CL            SET STRESS TO ZERO IF UNSTABLE
                N_SQ = 0.0
                STRESS(I,KU) = 0.0
              ELSE
                N   = SQRT( N_SQ )

C------------------------------------------------------------------
CL              2.2 CALCULATE MINIMUM RICHARDSON NO. DUE TO GRAVITY
CL                WAVES     EQNS 2.6 AND 2.8
C------------------------------------------------------------------

                AMPL =
     *            SQRT(STRESS(I,KU)/( KAY*RHO*SPEEDB*N ) )
                DV_DZ = SQRT( (U(I,K)-U(I,K-1))*(U(I,K)-U(I,K-1)) +
     *                   (V(I,K)-V(I,K-1))*(V(I,K)-V(I,K-1)) )/DZB
                RI   = N_SQ*(1. - N*AMPL/SPEEDB)/
     *           ((DV_DZ + N_SQ*AMPL/SPEEDB)*(DV_DZ + N_SQ*AMPL/SPEEDB))

C------------------------------------------------------------------
CL              2.3 TEST FOR WAVE-BREAKING (RI < RIC )
CL                AND MODIFY STRESS AT UPPER LAYER BOUNDARY
CL                EQNS 2.9 , 2.10 AND 2.6
C------------------------------------------------------------------

                IF( RI.LT.RIC ) THEN
                  EPSILON =
     *              ( SQRT(4.*RIC*DV_DZ*N + N_SQ*(1.+4.*RIC))
     *             -(2.*RIC*DV_DZ + N ) )/(2.*RIC)
                  AMPL = EPSILON*SPEEDB/ N_SQ
                  IF( AMPL.LT.0.0 ) THEN
                    STRESS(I,KU) = 0.0
                  ELSE
                    STRESS(I,KU) = KAY*RHO*N*SPEEDB*AMPL*AMPL

                  END IF     ! AMPL < 0 ELSE AMPL > 0

                END IF     ! RI < RIC

              END IF     ! N_SQ < 0 ELSE N_SQ > 0

            END IF     ! SPEED < 0 ELSE SPEED  >0

          END IF     ! STRESS >0

            IF( K .EQ. START_L+1 ) THEN
              DELTA_P = DELTA_AK_SUM+DELTA_BK_SUM*PSTAR(I)
            ELSE
              DELTA_P = DELTA_AK(K-1)+DELTA_BK(K-1)*PSTAR(I)
            END IF

C------------------------------------------------------------------
CL              2.4 CALCULATE DRAG FROM VERTICAL STRESS CONVERGENCE
CL                AND ACCELERATIONS FOR WIND COMPONENTS
CL                EQN 2.1
C------------------------------------------------------------------

            DRAG(I) = G*(STRESS(I,KU) - STRESS(I,KL) )/DELTA_P
            DU_DT(I,K-1) = -DRAG(I)*SIN_A(I)
            DV_DT(I,K-1) = -DRAG(I)*COS_A(I)

        END DO

        IF (STRESS_UD_ON) THEN
          DO I=1,POINTS
            STRESS_UD_LAND(I,K) = STRESS(I,KU)*SIN_A(I)
          ENDDO
        ENDIF
        IF (STRESS_VD_ON) THEN
          DO I=1,POINTS
            STRESS_VD_LAND(I,K) = STRESS(I,KU)*COS_A(I)
          ENDDO
        ENDIF

        IF( K .EQ. START_L+1 .AND. START_L .GT.1 ) THEN
C         SET ACCELERATION SAME IN ALL LAYERS UP TO START_L
          DO KK=1,START_L-1
            DO I=1,POINTS
              DU_DT(I,KK) = DU_DT(I,START_L)
              DV_DT(I,KK) = DV_DT(I,START_L)
            END DO
          END DO
        END IF

C Swap storage for lower and upper layers
        KK=KL
        KL=KU
        KU=KK

C Replace top half height of lower layer ready for next pass
        DO I=1,POINTS
          DZ(I,KL)=DZ(I,KT)
        ENDDO

      END DO
CL   END LOOP LEVELS

C------------------------------------------------------------------
CL    3 TOP OF MODEL. SET ACCELERATION SAME AS PENULTIMATE LAYER
CL      WITH PROVISO  THAT STRESS >= 0
C------------------------------------------------------------------

      DO I=1,POINTS
        DELTA_P   = DELTA_AK(LEVELS) + DELTA_BK(LEVELS)*PSTAR(I)
        STRESS(I,KU) = STRESS(I,KL) + DRAG(I)*DELTA_P/G
        IF( STRESS(I,KU) .LT. 0.0) STRESS(I,KU) = 0.0
        DRAG(I) = G*(STRESS(I,KU) - STRESS(I,KL) )/DELTA_P
        DU_DT(I,LEVELS) = -DRAG(I)*SIN_A(I)
        DV_DT(I,LEVELS) = -DRAG(I)*COS_A(I)
      END DO

      IF (STRESS_UD_ON) THEN
        DO I=1,POINTS
          STRESS_UD_LAND(I,LEVELS+1) = STRESS(I,KU)*SIN_A(I)
        END DO
      ENDIF
      IF (STRESS_VD_ON) THEN
        DO I=1,POINTS
          STRESS_VD_LAND(I,LEVELS+1) = STRESS(I,KU)*COS_A(I)
        END DO
      ENDIF

      RETURN
      END

