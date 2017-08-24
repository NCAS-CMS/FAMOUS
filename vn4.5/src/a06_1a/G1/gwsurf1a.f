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
CLL  SUBROUTINE GW_SURF------------------------------------------
CLL
CLL  PURPOSE:   TO CALCULATE 'SURFACE' STRESS DUE TO SUBGRID-SCALE
CLL             OROGRAPHIC GRAVITY WAVES.
CLL             THE SURFACE STRESS IS CALCULATED FROM THE VARIANCE OF
CLL             THE OROGRAPHY.THE MAXIMUM WAVE AMPLITUDE IS EITHER
CLL             LIMITED BY A LIMIT TO THE S.D. OF THE OROGRAPHY OR
CLL             BY THE FROUDE NUMBER RESTRICTED TO < OR = 1
CLL  SUITABLE FOR SINGLE COLUMN USE
CLL
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL  FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 1, DATED 12/09/89
CLL
CLL  SYSTEM TASK: PART OF P22
CLL
CLL C.WILSON    <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL D.GREGORY   <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.4   12/04/94  DEF FROUDE replaced by LOGICAL LFROUDE
CLL                                                 S.J.Swarbrick
CLL
CLL  DOCUMENTATION:        THE EQUATIONS USED ARE (2.2) TO (2.5)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO  22
CLL                        C. A. WILSON AND R. SWINBANK
CLL                        VERSION 1,DATED 15/12/89
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE GW_SURF
     1  (PSTAR,PEXNER,THETA,U,V,SD_OROG,S_STRESS,LEVELS,POINTS,
     2   AK,BK,AKH,BKH,KAY,SIN_A,COS_A,LFROUDE)

      IMPLICIT NONE

      INTEGER
     * LEVELS              !IN    NUMBER OF MODEL LEVELS
     *,POINTS              !IN    NUMBER OF POINTS
C

      REAL
     * PSTAR(POINTS)                    !IN    PSTAR FIELD
     *,PEXNER(POINTS,LEVELS+1)          !IN    PEXNER
     *,THETA(POINTS,LEVELS)             !IN    THETA FIELD
     *,U(POINTS,LEVELS)                 !IN    U FIELD
     *,V(POINTS,LEVELS)                 !IN    V FIELD
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
      REAL
     * AK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,BK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER boundary
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER boundary
     *,KAY                    !IN    surface stress constant (m-1)
     *,SD_OROG(POINTS)        !IN    STANDARD DEVIATION OF OROGRAPHY
     *,S_STRESS(POINTS)       !OUT   'SURFACE' STRESS
     *,SIN_A(POINTS)          !OUT   SIN (WIND DIRECTION FROM NORTH)
     *,COS_A(POINTS)          !OUT   COS (WIND DIRECTION FROM NORTH)
C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS:
C   0 REAL ARRAYS OF FULL FIELD LENGTH REQUIRED
C
      REAL
     * Z            ! HEIGHT DIFFERENCE
     *,RHO          ! DENSITY
     *,SPEED        ! WIND SPEED
     *,SPEED_SQ     ! SQUARE OF WIND SPEED
     *,N            ! BRUNT_VAISALA FREQUENCY
     *,N_SQ         ! SQUARE OF BRUNT_VAISALA FREQUENCY
     *,AMPL_SQ      ! SQUARE OF WAVE-AMPLITUDE
     *,FROUDEC      ! MAX SQUARE OF WAVE-AMPLITUDE (FROUDE <=1)
C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
C     NONE
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED = POINTS
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C
      INTEGER   I      ! LOOP COUNTER IN ROUTINE
C
      LOGICAL LFROUDE  ! Logical switch
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


      REAL
     &    PU,PL,PU3,PL3
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
CL    START LOOP POINTS
C------------------------------------------

      DO I=1,POINTS

C-------------------------------------------------------------------
CL    1.  CALCULATE HEIGHT DIFFERENCE Z3-Z1  EQN 2.4
C------------------------------------------ -------

        PU=PSTAR(I)*BKH(2) + AKH(2)
        PL=PSTAR(I)*BKH(1) + AKH(1)
        PU3=PSTAR(I)*BKH(4) + AKH(4)
        PL3=PSTAR(I)*BKH(3) + AKH(3)

        Z=((P_EXNER_C( PEXNER(I,2),PEXNER(I,1),PU,PL,KAPPA)
     &      -PEXNER(I,2))*THETA(I,1)
     *          +(PEXNER(I,2)-PEXNER(I,3))*THETA(I,2)
     *      +(PEXNER(I,3)-
     &         P_EXNER_C( PEXNER(I,4),PEXNER(I,3),PU3,PL3,KAPPA))
     *       *THETA(I,3))*CP/G

C------------------------------------------------------------------
CL    2 CALCULATE DENSITY AT LEVEL 2
C------------------------------------------------------------------

        RHO=( AK(2) + BK(2)*PSTAR(I) )/( R*THETA(I,2)*
     &         P_EXNER_C( PEXNER(I,3),PEXNER(I,2),PL3,PU,KAPPA))

C------------------------------------------------------------------
CL    3 CALCULATE WIND SPEED AND DIRECTION AT LEVEL 2
C------------------------------------------------------------------

        SPEED_SQ = U(I,2)*U(I,2) + V(I,2)*V(I,2)
        SPEED    = SQRT( SPEED_SQ )
        IF(SPEED.GT.0.0) THEN
          SIN_A(I)    = U(I,2)/SPEED
          COS_A(I)    = V(I,2)/SPEED
        ENDIF

C------------------------------------------------------------------
CL    4 CALCULATE BRUNT-VAISALA FREQUENCY EQN 2.3
C------------------------------------------------------------------

        N_SQ = G*( THETA(I,3) - THETA(I,1) )/
     &             ( 0.5*( THETA(I,3) + THETA(I,1) )*Z )
        IF( N_SQ .LT.0.0 ) N_SQ = 0.0
C       NO WAVES IF UNSTABLE
        N   = SQRT( N_SQ )

C------------------------------------------------------------------
CL    5 LIMIT SQUARE OF WAVE AMPLITUDE BY FROUDE NO.CRITERION
CL      OR BY IMPOSED VARIANCE MAXIMUM
C------------------------------------------------------------------

       AMPL_SQ = SD_OROG(I) * SD_OROG(I)

      IF (LFROUDE) THEN

       IF( N_SQ .GT.0.0 ) THEN
         FROUDEC = SPEED_SQ / N_SQ
       ELSE
         FROUDEC = 0.0
       ENDIF
       IF( AMPL_SQ . GT. FROUDEC ) AMPL_SQ = FROUDEC

      ELSE

       IF( AMPL_SQ . GT. VAR_MAX ) AMPL_SQ = VAR_MAX

      END IF

C------------------------------------------------------------------
CL    6 CALCULATE 'SURFACE' STRESS EQN 2.2a
C------------------------------------------------------------------

        S_STRESS(I) = KAY*RHO*N*SPEED*AMPL_SQ
      END DO
CL   END LOOP POINTS

      RETURN
      END
