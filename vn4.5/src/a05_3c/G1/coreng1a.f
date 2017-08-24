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
CLL  SUBROUTINE COR_ENGY-----------------------------------------------
CLL
CLL  PURPOSE : TO ADJUST THE POTENTIAL TEMPERATURE INCREMENTS
CLL            TO ENSURE THE CONSERVATION OF MOIST STATIC ENERGY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.5  22/7/98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (12)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE COR_ENGY
     *          (NP_FIELD,NPNTS,NCORE,NLEV,DTHBYDT,DQBYDT,SNOW,         
     *                   EXNER,PSTAR,DELAK,DELBK,AKH,BKH,INDEX4)        
C
      IMPLICIT NONE
C
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO     
                                  ! SPECIFY STARTING POINT OF        
                                  ! DATA PASSED IN)             
C
      INTEGER NCORE               ! IN VECTOR LENGTHS
C
      INTEGER NPNTS               ! IN FULL VECTOR LENGTH
C
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I,K                 ! LOOP COUNTERS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      INTEGER INDEX4(NPNTS) 
      REAL DQBYDT(NP_FIELD,NLEV)  ! IN INCREMENT TO MODEL MIXING       
                                  !    RATIO DUE TO CONVECTION
                                  !    (KG/KG/S)
C
      REAL SNOW(NP_FIELD)         ! IN SNOW AT SURFACE (KG/M**2/S)      
C
      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER RATIO                      
C
      REAL PSTAR(NP_FIELD)        ! IN SURFACE PRESSURE (PA)            
C
      REAL DELAK(NLEV),           ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     *     DELBK(NLEV)            !    COEFFICIENTS A AND B
                                  !    ACROSS LAYER K
C
      REAL AKH(NLEV+1)              ! IN Hybrid coordinate A at
                                    !    layer boundary
      REAL BKH(NLEV+1)              ! IN Hybrid coordinate B at
                                    !    layer boundary
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL DTHBYDT(NP_FIELD,NLEV) ! INOUT              
                                  ! IN  INCREMENT TO MODEL POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
                                  ! OUT CORRECTED INCREMENT TO MODEL
                                  !     POTENTIAL TEMPERATURE DUE TO
                                  !     CONVECTION (K/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C
      REAL QSUM(NCORE)            ! SUMMATION OF INCREMENTS TO MODEL
                                  ! MIXING RATIO DUE TO CONVECTION
                                  ! IN THE VERTICAL, WEIGHTED
                                  ! ACCORDING TO THE MASS OF THE
                                  ! LAYER (KG/M**2/S)
C
      REAL TSPOS(NCORE)           ! SUMMATION OF POSITIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      REAL TSNEG(NCORE)           ! SUMMATION OF NEGATIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      REAL TERR(NCORE)            ! SUMMATION OF ALL INCREMENTS TO
                                  ! MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      LOGICAL BPOSER(NCORE)       ! MASK FOR POINTS IN LAYER K AT WHICH
                                  ! INCREMENTS TO MODEL POTENTIAL
                                  ! TEMPERATURE DUE TO CONVECTION ARE
                                  ! POSITIVE
C
      LOGICAL BCORR(NCORE)        ! MASK FOR POINTS AT WHICH ENTHALPY
                                  ! CORRECTION IS NECESSARY
C
      REAL DELPK                  ! DIFFERENCE IN PRESSURE ACROSS A
                                  ! LAYER (PA)
C
      REAL EXTEMPK                ! EXNER RATIO AT THE MID-POINT OF
                                  ! LAYER K
C

      REAL
     &    PU,PL
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


C*----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  SUM UP MIXING RATIO AND +VE AND -VE TEMPERATURE INCREMENTS
CL----------------------------------------------------------------------
CL
      DO 20 I=1,NCORE
       QSUM (I) = 0.0
       TSPOS(I) = 0.0
       TSNEG(I) = 0.0
   20  CONTINUE
C
      DO 40 K=1,NLEV
       DO 30 I=1,NCORE
C
        DELPK = -DELAK(K) - DELBK(K)*PSTAR(INDEX4(I))   
C
        PU=PSTAR(INDEX4(I))*BKH(K+1) + AKH(K+1)                         
        PL=PSTAR(INDEX4(I))*BKH(K) + AKH(K)                             
        EXTEMPK  = 
     &    P_EXNER_C(EXNER(INDEX4(I),K+1),
     &              EXNER(INDEX4(I),K),PU,PL,KAPPA)    
C
        QSUM(I) = QSUM(I) + DQBYDT(INDEX4(I),K)*DELPK      
C
        IF (DTHBYDT(INDEX4(I),K) .GT. 0.0) THEN                         
           TSPOS(I) = TSPOS(I) + 
     &                DTHBYDT(INDEX4(I),K)*(CP*DELPK*EXTEMPK)        
        ELSE
           TSNEG(I) = TSNEG(I) + 
     &                DTHBYDT(INDEX4(I),K)*(CP*DELPK*EXTEMPK)      
        ENDIF
   30  CONTINUE
   40 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE ERROR AND APPLY THE NECESSARY CORRECTION
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (12), EQUATION (48), (49)
CL----------------------------------------------------------------------
CL
      DO 50 I=1,NCORE
C
       TERR(I) = LC*QSUM(I) - LF*G*SNOW(INDEX4(I)) + 
     &                                   TSPOS(I) + TSNEG(I)   
C
       BPOSER(I) = TERR(I) .GT. 0.0
C
       IF (BPOSER(I) .AND. (TSPOS(I) .EQ. 0.0)) THEN
          BPOSER(I) = .FALSE.
       ELSE IF (.NOT.BPOSER(I) .AND. (TSNEG(I) .EQ. 0.0)) THEN
          BPOSER(I) = .TRUE.
       ENDIF
C
       BCORR(I) = (TSPOS(I) .NE. 0.0) .OR. (TSNEG(I) .NE. 0.0)
C
       IF (BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSPOS(I)
       ELSE IF (.NOT.BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSNEG(I)
       ENDIF
C
  50  CONTINUE
C
      DO 100 K=1,NLEV
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
       DO 100 I=1,NCORE
        IF (BCORR(I) .AND. (( BPOSER(I) .AND.
     *   (DTHBYDT(INDEX4(I),K) .GT. 0.0)) .OR. ( .NOT.BPOSER(I)         
     *   .AND. (DTHBYDT(INDEX4(I),K) .LT. 0.0))))                       
     *       DTHBYDT(INDEX4(I),K) = DTHBYDT(INDEX4(I),K)*TERR(I)        
  100      CONTINUE
C
      RETURN
      END
