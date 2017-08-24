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
CLL   SUBROUTINE THETL_QT -----------------------------------------
CLL
CLL   PURPOSE: CALCULATES THETAL AND QT FROM THETA,Q, AND QC USING
CLL            EQUATIONS (10) AND (11),SUBTRACTS THETA_REF FROM THETAL.
CLL   NOT SUITABLE FOR I.B.M USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN  BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD A. VERSION 2, DATED 18/01/90
CLL
CLL   SYSTEM COMPONENTS COVERED: P192
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       EQUATIONS (10) AND (11)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE THETL_QT
     1                   (PSTAR,THETA,Q,QCL,QCF,P_EXNER,AKH,BKH,
     2                    P_FIELD,P_LEVELS,Q_LEVELS)

      IMPLICIT NONE

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, P_LEVELS           !IN NUMBER OF MODEL LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST MODEL LEVELS.

      REAL
     * PSTAR(P_FIELD)            !IN    PSTAR FIELD
     *,THETA(P_FIELD,P_LEVELS)   !INOUT THETA FIELD IN, THETAL OUT.
     *,Q(P_FIELD,Q_LEVELS)       !INOUT Q FIELD IN, QT FIELD OUT.

      REAL
     * QCL(P_FIELD,Q_LEVELS)       !IN QCL FIELD.
     *,QCF(P_FIELD,Q_LEVELS)       !IN QCF FIELD.
     *,P_EXNER(P_FIELD,P_LEVELS+1) !IN Exner Pressure on half levels.
     *,AKH(P_LEVELS+1)             !IN Hybrid Coords. A and B values
     *,BKH(P_LEVELS+1)             !IN for half levels.
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: NONE ARE REQUIRED

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

      REAL
     *  P_EXNER_FULL !HOLDS EXNER PRESSURE AT FULL LEVEL.
     *  ,PKP1,PK     !Pressures at half levels k+1 and k.

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,K

C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------
CL    CALL COMDECK TO GET CONSTANTS USED.

CLL   COMDECK C_THETLQ HOLDS CONSTANTS FOR ROUTINE THETL_QT
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

CL    END OF COMDECK C_THETLQT

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
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


CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE THETAL AND QT.
CL---------------------------------------------------------------------

CL LOOP OVER MOIST LEVELS IE: Q_LEVELS.

      DO 100 K=1,Q_LEVELS

CFPP$ SELECT(CONCUR)
        DO 110 I= 1,P_FIELD

          PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I)
          PK   = AKH(K)   + BKH(K)  *PSTAR(I)
          P_EXNER_FULL = P_EXNER_C
     *    (P_EXNER(I,K+1),P_EXNER(I,K),PKP1,PK,KAPPA)
          THETA(I,K) = THETA(I,K) -  (LC*QCL(I,K)+(LC+LF)*QCF(I,K))
     *                 /(CP*P_EXNER_FULL)
          Q(I,K) = Q(I,K) + QCL(I,K) + QCF(I,K)
 110    CONTINUE

CL END LOOP OVER MOIST LEVELS.
 100  CONTINUE

CL    END OF ROUTINE THETL_QT

      RETURN
      END
