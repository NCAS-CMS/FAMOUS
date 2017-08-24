C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE FLAG_WET-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES A MASK FOR WHEN CONDENSATION IS LIQUID
CLL
CLL            IF 0.5 * (TK + TK+1) > TICE THEN ANY CONDENSATION
CLL                                        IN LAYER K+1 IS LIQUID
CLL
CLL            IF 0.5 * (TK + TK+1) < TICE THEN ANY CONDENSATION
CLL                                        IN LAYER K+1 IS ICE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
!LL   4.5   03/03/98  Correct level swapping in loop. R. Rawlins
CLL
CLL  PROGRAMMING STANDARDS :
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL                  SECTION (2B)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE FLAG_WET (BWATER,TH,EXNER,PSTAR,AKH,BKH,
     &                     NP_FIELD,NPNTS,NLEV)
C
C-----------------------------------------------------------------------
C   RETURNS 'BWATER' - A BIT VECTOR OF POINTS WHERE CONDENSATE IS WATER
C   RATHER THAN ICE.
C----------------------------------------------- AUTHOR: M FISHER 1987 -
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      REAL TICE  ! TEMPERATURE IN KELVIN AT WHICH LIQUID WATER TURNS
                 ! TO ICE IN THE CONVECTIVE DOWNDRAUGHT SCHEME (K)
      PARAMETER (TICE = 273.15)
C
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

C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NP_FIELD           ! IN FULL VECTOR LENGTH
C
      INTEGER NPNTS              ! IN VECTOR LENGTH
C
      INTEGER NLEV               ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I,K                ! LOOP COUNTERS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL TH(NP_FIELD,NLEV)        ! IN POTENTIAL TEMPERATURE (K)
C
      REAL EXNER(NP_FIELD,NLEV+1)   ! IN EXNER RATIO AT LAYER
                                    ! BOUNDARIES (STARTING WITH THE
                                    ! SURFACE)
C
      REAL PSTAR(NPNTS)             ! IN Surface pressure
C
      REAL AKH(NLEV+1)              ! IN Hybrid coordinate A at
                                    !    layer boundary
      REAL BKH(NLEV+1)              ! IN Hybrid coordinate B at
                                    !    layer boundary
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      LOGICAL BWATER(NPNTS,2:NLEV)  ! OUT MASK FOR THOSE POINTS AT
                                    !     WHICH CONDENSATE IS LIQUID
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EXK                      ! EXNER RATIO FOR LEVEL K
      REAL EXKP1                    ! EXNER RATIO FOR LEVEL K+1
      REAL EX(NPNTS,2)  ! EXNER pressure at levels K & K+1
      INTEGER KL,K1,K2
C

      REAL
     &    PU,PL,PU2
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


C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
      K=1
      K1=1
      K2=2
      DO I=1,NPNTS
        PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
        PL=PSTAR(I)*BKH(K) + AKH(K)
        EX(I,K1) = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
      ENDDO
      DO K=1,NLEV-1
       DO I=1,NPNTS

        PU2=PSTAR(I)*BKH(K+2) + AKH(K+2)
        PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
        EX(I,K2) = P_EXNER_C(EXNER(I,K+2),EXNER(I,K+1),PU2,PU,KAPPA)

        BWATER(I,K+1) = 0.5*(TH(I,K)*EX(I,K1) +
     *           TH(I,K+1)*EX(I,K2)) .GT. TICE
       ENDDO
      KL=K2
      K2=K1
      K1=KL
      ENDDO
C
      RETURN
      END
