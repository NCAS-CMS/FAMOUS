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
CLL  SUBROUTINE V_INT_TH----------------------------------------------
CLL
CLL  Purpose:  Calculates potential temperature of model layers using
CLL            temperatures input on an arbitray set of pressure levels
CLL
CLL R.Swinbank  <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  4.2     01/08/96  *DEF CRAY removed. Use of 32-bit log and exp
CLL                     for CRAY shared memory computers no longer 
CLL                     supported.
CLL
CLL                    Author: A. Dickinson    Reviewer: R. Rawlins
CLL
CLL Programming standard :
CLL
CLL Logical components covered : S113
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE V_INT_TH
     *(P_HALF,P_EXNER_HALF,THETA,T_SRCE,POINTS,P_SRCE,P_LEVELS,
     * SRCE_LEVELS,PSTAR,AKH,BKH)

      IMPLICIT NONE

      INTEGER
     * POINTS      !IN Number of points to be processed
     *,P_LEVELS    !IN Number of model levels
     *,SRCE_LEVELS !IN Number of input levels

      REAL
     * P_HALF(POINTS,P_LEVELS+1) !IN Pressure of model half levels
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,PSTAR(POINTS)              !IN surface pressure
     *,AKH(P_LEVELS+1)            !IN Hybrid coord. A at half levels
     *,BKH(P_LEVELS+1)            !IN Hybrid coord. B at half levels
     *,THETA(POINTS,P_LEVELS) !OUT Potential temperature at full levels
     *,T_SRCE(POINTS,SRCE_LEVELS) !IN Input temperature fields
     *,P_SRCE(POINTS,SRCE_LEVELS) !IN Pressure of input temperature
     *                            !   fields

C Workspace usage:-----------------------------------------------------
      REAL TL(POINTS),PL(POINTS)
      INTEGER KI(POINTS)
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,J,K,JK,JI
      REAL PUB           !  Pressure at top of layer
      REAL PLB           !  Pressure at bottom of layer
      REAL P_EXNER_FULL  !  Exner pressure at full model levels
      REAL PKPH,TKPH,THICK,PJ,PJM1,TJ,TJM1
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
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
      REAL LAPSE_R_OVER_G
      PARAMETER(LAPSE_R_OVER_G=LAPSE*R/G)
C----------------------------------------------------------------------

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


CL 1. Initialise arrays PL,TL,KI and THETA

      DO I=1,POINTS
        KI(I)=1
        TL(I)=0.0
        PL(I)=P_HALF(I,1)
      ENDDO

C THETA is used to accumulate model layer thicknesses (x2) in section 2
C The values are converted to potential temperatures in section 3.
      DO J=1,P_LEVELS
        DO I=1,POINTS
          THETA(I,J)=0.0
        ENDDO
      ENDDO

CL 2. Loop over JK. Each time around we either increment one input
CL    level (J) or one model level (K), as appropriate; JK=J+K.
CL    For each combination j,k we calculate the thickness of the
CL    overlap between the interval P_SRCE(j-1 to j) and model layer
CL    k (i.e. P_HALF(k to k+1)).

      DO 100 JK=2,P_LEVELS+SRCE_LEVELS+1

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO I=1,POINTS

C Get J (input level)
C Max. value of J is SRCE_LEVELS+1 (above top input level).
C JI is J value for interpolation - limited to SRCE_LEVELS, since
C we extrapolate above top level.

      K=KI(I)
      J=MIN(JK-K,SRCE_LEVELS+1)
      JI=MIN(J,SRCE_LEVELS)

C Gather P,T values for interpolation
      PJ=P_SRCE(I,JI)
      TJ=T_SRCE(I,JI)
      PJM1=P_SRCE(I,JI-1)
      TJM1=T_SRCE(I,J-1)

C If K>P_LEVELS integration is complete;
C if input pressure is greater than range of model levels integration is
C not yet started.
C In both cases no action is required.
      IF (K.LE.P_LEVELS .AND. PJ.LT.P_HALF(I,1)) THEN

C If TL is not set, integration has just started, so we require T at
C model level 1/2 (P_HALF(1)), which must be in current input pressure
C interval.

        IF(TL(I).EQ.0.0) THEN
          IF(J.EQ.1) THEN
C If J=1, extrapolate below input level 1
            TL(I)=TJ*(P_HALF(I,1)/PJ)**LAPSE_R_OVER_G
          ELSE
C Otherwise interpolate between levels J-1 and J
            TL(I)=TJM1+ALOG(P_HALF(I,1)/PJM1)*(TJ-TJM1)
     *    /ALOG(PJ/PJM1)
          END IF
        END IF

        PKPH=P_HALF(I,K+1)

C Check whether next level up (from PL) is a model layer boundary, or
C an input level; then integrate from PL to this level.

        IF(PJ.GE.PKPH.AND.J.LE.SRCE_LEVELS) THEN
C Integrate from PL to input level J

          THICK=(TL(I)+TJ)*ALOG(PL(I)/PJ)

C Save P & T values at current level
          PL(I)=PJ
          TL(I)=TJ

        ELSE
C Get T at model level K+1/2
          IF(J.EQ.1) THEN
C If J=1, extrapolate below input level 1

            TKPH=TJ*(PKPH/PJ)**LAPSE_R_OVER_G

          ELSE
C Otherwise interpolate between levels J-1 and J

            TKPH=TJM1+ALOG(PKPH/PJM1)*(TJ-TJM1)/ALOG(PJ/PJM1)

          END IF

C Integrate from PL to model level k+1/2
 
         THICK=(TL(I)+TKPH)*ALOG(PL(I)/PKPH)
C Save P & T values and increment K, so next model level is used
C on next iteration of JK
          PL(I)=PKPH
          TL(I)=TKPH
          KI(I)=K+1
        END IF

C Integration results are accumulated in THETA; they are converted
C to theta values in section 3.
        THETA(I,K)=THETA(I,K)+THICK

      END IF

      ENDDO
 100  CONTINUE


CL 3. Convert thickness values to potential temperature values.

      DO J=1,P_LEVELS
        DO I=1,POINTS

C Model layer above source data coverage: equation (3.28)
          IF(P_HALF(I,J).LE.P_SRCE(I,SRCE_LEVELS))THEN
            PUB=PSTAR(I)*BKH(J+1) + AKH(J+1)
            PLB=PSTAR(I)*BKH(J) + AKH(J)
            THETA(I,J)=T_SRCE(I,SRCE_LEVELS)/
     &      P_EXNER_C( P_EXNER_HALF(I,J+1),P_EXNER_HALF(I,J),
     &      PUB,PLB,KAPPA )

          ELSE
C Convert layer thickness to theta: finalise equation (3.27)
            THETA(I,J)=KAPPA*THETA(I,J)
     *          /(2.*(P_EXNER_HALF(I,J)-P_EXNER_HALF(I,J+1)))
          ENDIF

C If model layer straddles top srce level and interpolated temperature
C falls outside the range of the top two source level temperatures,
C then use equ 3.28

          IF((P_HALF(I,J).GT.P_SRCE(I,SRCE_LEVELS)).AND.
     *          (P_HALF(I,J+1).LE.P_SRCE(I,SRCE_LEVELS)))THEN

            PUB = PSTAR(I)*BKH(J+1) + AKH(J+1)
            PLB = PSTAR(I)*BKH(J)   + AKH(J)
            P_EXNER_FULL = P_EXNER_C
     *      ( P_EXNER_HALF(I,J+1),P_EXNER_HALF(I,J),PUB,PLB,KAPPA )
            TJ=THETA(I,J)*P_EXNER_FULL
            IF((TJ.GT.T_SRCE(I,SRCE_LEVELS).AND.
     *      TJ.GT.T_SRCE(I,SRCE_LEVELS-1)).OR.
     *      (TJ.LT.T_SRCE(I,SRCE_LEVELS).AND.
     *      TJ.LT.T_SRCE(I,SRCE_LEVELS-1)))THEN
              THETA(I,J)=T_SRCE(I,SRCE_LEVELS)/P_EXNER_FULL
            ENDIF

          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
