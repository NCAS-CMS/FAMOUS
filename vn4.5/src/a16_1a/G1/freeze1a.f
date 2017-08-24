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
CLL  SUBROUTINE FREEZE------------------------------------------
CLL
CLL  PURPOSE:   Calculates the true height and pressure of a temperature
CLL             surface with temperature T0.
CLL             T0 set in PHY_DIAG. T0=273.16K for freezing level
CLL                                 T0=253.16K for -20 degree C level
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL J.Heming    <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL P.Smith     <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.4     23/8/94  Correct search for T0 to find highest occurence
CLL                   avioding low level inversions               PJS
!LL   4.5    15/04/98 Start-end args added to enable dupicate halo
!LL                   calculations to be avoided. S.D.Mullerworth
CLL
CLL  Logical components covered D422,D421,D423
CLL  Project TASK: D4
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation
CLL
CLLEND-------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE FREEZE(
C data in
     & T0,P,THETA,T,P_EXNER_HALF,PSTAR,Q,MODEL_HALF_HEIGHT,OROG,
C data out
     & Z_AT_T0,P_AT_T0,
C constants in
     & POINTS,P_LEVELS,Q_LEVELS,L,AKH,BKH,START,END)
C*
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      EXTERNAL  V_INT_Z
C-----------------------------------------------------------------------
      INTEGER
     * POINTS         ! IN  NO OF POINTS
     *,P_LEVELS       ! IN  NO OF MODEL LEVELS
     *,Q_LEVELS       ! IN  NO OF HUMIDITY LEVELS
     *,L              ! IN  LEVEL OF MODEL USED TO CALC HEIGHT
     &,START,END      ! IN  Range of points to calculate
C-----------------------------------------------------------------------
      REAL
     * P(POINTS,P_LEVELS)              ! IN  PRESSURE  AT FULL LEVELS
     *,THETA(POINTS,P_LEVELS)          ! IN  POTENTIAL TEMPERATURE " " "
     *,T(POINTS,P_LEVELS)              ! IN  TEMPERATURE AT FULL LEVELS
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) ! IN  EXNER PRESSURE AT MODEL
     *                                 !     HALF LEVELS
     *,PSTAR(POINTS)                   ! IN  SURFACE PRESSURE
     *,Q(POINTS,Q_LEVELS)              ! SPECIFIC HUM AT FULL LEVELS
     *,MODEL_HALF_HEIGHT(POINTS,P_LEVELS+1)!IN  HEIGHT OF HALF LVLS
     *,OROG(POINTS)                    ! IN  MODEL OROGRAPHY
     *,T0                              ! IN  TEMP OF SURFACE IN K
     *,AKH(P_LEVELS+1)                 ! IN  A values at half levels.
     *,BKH(P_LEVELS+1)                 ! IN  B values at half levels.
     *,P_AT_T0(POINTS)             ! OUT PRESSURE AT LEVEL WITH TEMP T0
     *,Z_AT_T0(POINTS)             ! OUT HEIGHT AT LEVEL WITH TEMP T0
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,K               !  LOOP COUNTERS
     *,J                 !  LOWEST ETA LEVEL WITH TEMP<T0
      REAL
     * PJP1, PJ, PJM1    ! Pressure at half levels J+1,J,J-1
     *,P_EXNER_FULL_JM1  ! EXNER PRESSURE AT FULL LEVEL J-1
     *,P_EXNER_FULL_J    !   "      "     "   "     "   J
     *,DEL_EXNER_JM1     ! EXNER PRESSURE DIFF BET LEVELS
     *                   !                  J-3/2 AND J-1/2
     *,DEL_EXNER_J       !   "     "  " " "J-1/2 AND J+1/2
     *,TERM1,TERM2       !
     *,LAPSE_RATE        ! LAPSE RATE BETWEEN LEVELS J-1 AND J
     *,Z(POINTS,P_LEVELS)! HEIGHT OF POINTS WITH PRESSURE P
C-----------------------------------------------------------------------
C   Note: these variables are temporary
C-----------------------------------------------------------------------
      INTEGER I_F
      REAL T1,MULT
C-----------------------------------------------------------------------
C Constants
C-----------------------------------------------------------------------
C*
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C-----------------------------------------------------------------------
      REAL CP_OVER_G
      PARAMETER(CP_OVER_G=CP/G)

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

C-----------------------------------------------------------------------
CL Calculate heights of full levels
C-----------------------------------------------------------------------
      DO K=1,P_LEVELS
        CALL V_INT_Z(P(1,K),P(1,L),PSTAR,P_EXNER_HALF,THETA,Q,
     &  MODEL_HALF_HEIGHT,Z(1,K),POINTS,P_LEVELS,Q_LEVELS,L,AKH,BKH
     &  ,START,END)
      ENDDO
C-----------------------------------------------------------------------
      T1=T0
      I_F=0
C-----------------------------------------------------------------------
CL Search upwards for levels with temperature less than T0 where level
CL  below is greater than T0. Highest of these is used for calulations.
C-----------------------------------------------------------------------
      DO 111 I=START,END
        J=0
        DO K=1,P_LEVELS
          IF (T(I,K).LT.T0) THEN
            IF (K.EQ.1) THEN
              J = K
            ELSEIF (T(I,K-1).GT.T0) THEN
              J = K
            ENDIF
          ENDIF
        ENDDO
C-----------------------------------------------------------------------
        IF (J.GE.2) THEN
          I_F=I_F+1
C-----------------------------------------------------------------------
CL Exner pressure at full levels
C-----------------------------------------------------------------------
          PJP1 = AKH(J+1) + BKH(J+1)*PSTAR(I)
          PJ   = AKH(J)   + BKH(J)  *PSTAR(I)
          PJM1 = AKH(J-1) + BKH(J-1)*PSTAR(I)
          P_EXNER_FULL_J = P_EXNER_C
     &    (P_EXNER_HALF(I,J+1),P_EXNER_HALF(I,J),PJP1,PJ,KAPPA)
          P_EXNER_FULL_JM1 = P_EXNER_C
     &    (P_EXNER_HALF(I,J),P_EXNER_HALF(I,J-1),PJ,PJM1,KAPPA)
C-----------------------------------------------------------------------
CL Exner pressure difference across half layers
C-----------------------------------------------------------------------
          DEL_EXNER_J=P_EXNER_HALF(I,J)-P_EXNER_FULL_J
          DEL_EXNER_JM1=P_EXNER_FULL_JM1-P_EXNER_HALF(I,J)
C-----------------------------------------------------------------------
CL Denominator
C-----------------------------------------------------------------------
          TERM2=CP_OVER_G*(THETA(I,J-1)*DEL_EXNER_JM1
     *       +THETA(I,J)*DEL_EXNER_J)
C-----------------------------------------------------------------------
CL Numerator
C-----------------------------------------------------------------------
          TERM1=THETA(I,J-1)*P_EXNER_FULL_JM1-THETA(I,J)*P_EXNER_FULL_J
C-----------------------------------------------------------------------
CL Lapse rate between level j-1 and j
C-----------------------------------------------------------------------
          LAPSE_RATE=TERM1/TERM2
          MULT=(T(I,J-1)-T0)/LAPSE_RATE
C-----------------------------------------------------------------------
CL Calculate the height and pressure at level with temperature T0
C-----------------------------------------------------------------------
          Z_AT_T0(I)=Z(I,J-1)+(T(I,J-1)-T0)/LAPSE_RATE
          P_AT_T0(I)=P(I,J-1)*(T0/T(I,J-1))**(G/(R*LAPSE_RATE))
        ELSE IF(J.EQ.1) THEN
          Z_AT_T0(I)=OROG(I)
          P_AT_T0(I)=PSTAR(I)
        ELSE ! J=0
          Z_AT_T0(I)=-1
          P_AT_T0(I)=-1
        ENDIF
 111  CONTINUE
C-----------------------------------------------------------------------
      RETURN
      END
