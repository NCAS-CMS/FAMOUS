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
C*LL  SUBROUTINE LSP_EVAP-----------------------------------------------
CLL
CLL  Purpose: Calculate the amount of evaporation from precipitation
CLL           falling through one model layer, and the effects of this
CLL           evaporation on Q and T in the layer.
CLL           This version uses the revised constants.
CLL
CLL  Rewritten by S BETT
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed, HF functions
CLL                    replaced.
CLL                                    S.J.Swarbrick
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 4, dated 5/2/92.
CLL
CLL  System component covered: Part of P26.
CLL
CLL  Documentation: Unified Model Documentation Paper No 26.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_EVAP
     +(P,RHODZ,TIMESTEP,POINTS,Q,RAIN,SNOW,T)
      IMPLICIT NONE
      INTEGER          ! Input integer scalar :-
     + POINTS          ! IN Number of gridpoints being processed.
      REAL             ! Input real arrays :-
     + P(POINTS)       ! IN pressure N /sq m
     +,RHODZ(POINTS)   ! IN Air mass p.u.a. in layer (kg per sq m).
      REAL             ! Updated real arrays :-
     + Q(POINTS)       ! INOUT Specific humidity (kg water per kg air).
     +,RAIN(POINTS)    ! INOUT Rainfall rate (kg per sq m per s).
     +,SNOW(POINTS)    ! INOUT Snowfall rate (kg per sq m per s).
     +,T(POINTS)       ! INOUT Temperature (K).
      REAL             ! Input real scalar :-
     + TIMESTEP        ! IN Timestep (s).
C*
C*L  Workspace usage: 1 real array--------------------------------------
      REAL
     + QS(POINTS)      !  Saturated sp humidity for (T,p) in layer
C*L  external subprograms are called ---------------------------------
      EXTERNAL QSAT
C*
C*
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

C  Local (derived) physical constants ----------------------------------
C*LL--------------------------------------------------------------------
CLL  Constants used by LSP_EVAP.
CLL
CLL  These are tunables in the calculations of bulk coefficients of
CLL  evaporation from precipitation.
CLL
CLL  Naming convention: Cxxn where xx is EV for evaporation from
CLL   rain, SB for evaporation from snow (i.e. sublimation); n  is power
CLL   of T to which the coefficient is applied, or A for the factor by
CLL   which the whole expression is multiplied.
C
      REAL CEVA,CEV0,CEV1,CEV2,CSBA,CSB0,CSB1,CSB2
      PARAMETER (           !
     + CEVA=567.            ! For
     +,CEV0=2.424E-4        ! naming
     +,CEV1=-1.385E-6       ! convention
     +,CEV2=2.008E-9        ! see
     +,CSBA=681.            ! COMDECK
     +,CSB0=-2.9111E-4      ! header
     +,CSB1=2.5332E-6       ! documentation.
     +,CSB2=-5.2E-9         !
     +)
C*----------------------------------------------------------------------
C
      REAL LSRN_A, LSRN_B       ! CONSTANTS USED IN EVAPORATION
                                ! OF LARGE-SCALE RAIN
      PARAMETER (LSRN_A=75.67, LSRN_B=605.17)
C
      REAL LSRN_P1, LSRN_P2, LSRN_P3, LSRN_P4  ! EXPONENTS USED IN
                                               ! CALCULATION OF EVAP
      PARAMETER (LSRN_P1=0.21, LSRN_P2=0.42, LSRN_P3=0.55,
     &           LSRN_P4=0.6)
C
      REAL LSSW_A, LSSW_B       ! CONSTANTS USED IN EVAPORATION
                                ! OF LARGE-SCALE RAIN
      PARAMETER (LSSW_A=1765.55, LSSW_B=34784.06)
C
      REAL LSSW_P1, LSSW_P2, LSSW_P3, LSSW_P4  ! EXPONENTS USED IN
                                               ! CALCULATION OF EVAP
      PARAMETER (LSSW_P1=0.28, LSSW_P2=0.55, LSSW_P3=0.63,
     &           LSSW_P4=0.76)
C
!
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      REAL LCRCP,LFRCP,LSRCP
      PARAMETER(
     + LCRCP=LC/CP           ! Latent heat of condensation / Cp (K).
     +,LFRCP=LF/CP           ! Latent heat of fusion / Cp (K).
     +,LSRCP=LCRCP+LFRCP     ! Sum of the above (S for Sublimation).
     +)
      REAL ALPHF,ALPHL                  ! Derived parameters.
      PARAMETER (
     + ALPHF=EPSILON*(LF+LC)/R          ! For frozen AlphaL calculation.
     +,ALPHL=EPSILON*LC/R               ! For liquid AlphaL calculation.
     +)
C
C  Define local variables-----------------------------------------------
      INTEGER I        ! Loop counter (horizontal field index).
C   6 local variables which will effectively be expanded to workspace
C   by the Cray (using vector registers) are required :-
      REAL             ! Real "workspace".  Contents at end of loop :-
     + CEV             ! Bulk evporation coefficient (rain).
     +,CSB             ! Bulk evporation coefficient (snow).
     +,QEV             ! Evap rate from rain (kg wat per kg air per s).
     +,QSB             ! Evap rate from snow (kg wat per kg air per s).
     +,ALPHAL          ! factor from Clausius-Claperyon eqn
     +,BL              ! factor due to implicit treatment
     +,C1              ! temporary store
     +,C2              ! temporary store
     +,RHO             ! density of air in layer
     +,SR_RHO          ! square root of density
C-----------------------------------------------------------------------
CL  Internal structure.
CL  0 Call qsat
C-----------------------------------------------------------------------
      CALL QSAT(QS,T,P,POINTS)

C-----------------------------------------------------------------------
CL  Loop over gridpoints :-
C-----------------------------------------------------------------------
      DO 1 I=1,POINTS
C
C-----------------------------------------------------------------------
C Calculate density of air in layer
C-----------------------------------------------------------------------
C
       RHO = P(I) / (R*T(I))
       SR_RHO = SQRT(RHO)
C
C-----------------------------------------------------------------------
CL  1. Perform calculations for rain.
C
C-----------------------------------------------------------------------
C
        IF(RAIN(I).GT.0.0)THEN
C-----------------------------------------------------------------------
CL  1.1 Calculate evaporation coefficient for rain (CEV).
CL      See eqs P26.9, P26.11.
C-----------------------------------------------------------------------
          CEV = ((CEV2*T(I)+CEV1)*T(I)+CEV0)
          C1 = LSRN_A*(RAIN(I)*SR_RHO)**LSRN_P2
          C2 = LSRN_B*(RHO**LSRN_P3)*(RAIN(I)**LSRN_P4)
          CEV = CEV*(C1+C2)
C-----------------------------------------------------------------------
CL  1.2 Calculate implicit treatment factors alphal and bl
CL      See eqs P26.??, P26.??.
C-----------------------------------------------------------------------
          ALPHAL=ALPHL*QS(I)/(T(I)*T(I))
          BL=1. + CEV*TIMESTEP*(1. + LCRCP*ALPHAL)
C-----------------------------------------------------------------------
CL  1.3 Calculate evaporation rate, adjusted to be .LE. precip rate,
CL      in kg water per sq m  per sec.  See eq P26.1.
C       Store result in QEV.  NB this is QEV*RHODZ in documentation.
C-----------------------------------------------------------------------
          QEV=MIN ( RAIN(I) , RHODZ(I)*CEV*MAX(0.0,QS(I)-Q(I))/BL )
C-----------------------------------------------------------------------
CL  1.4 Calculate effects of evaporation on precip rates and on Q and T.
CL     See eqs P26.7,  P26.5, P26.6 respectively.
C
C-----------------------------------------------------------------------
C  Increment precipitation rates (kg per sq m per sec), Q (kg per kg)
C  and T (K).  For the last 2 the change is integrated over the
C  timestep.
C
          RAIN(I)=RAIN(I)-QEV
          Q(I)=Q(I)+QEV*TIMESTEP/RHODZ(I)
          T(I)=T(I)-LCRCP*QEV*TIMESTEP/RHODZ(I)
         ENDIF ! rain>0
    1 CONTINUE
C-----------------------------------------------------------------------
CL  2 Call qsat again since evap of rain may have altered T
C-----------------------------------------------------------------------
      CALL QSAT(QS,T,P,POINTS)

      DO 3 I=1,POINTS
C
C-----------------------------------------------------------------------
C Recalculate density of air in layer
C-----------------------------------------------------------------------
C
        RHO = P(I) / (R*T(I))
        SR_RHO = SQRT(RHO)
c
C-----------------------------------------------------------------------
CL  3. Perform calculations for snow.
C
C-----------------------------------------------------------------------
C
        IF(SNOW(I).GT.0.0)THEN
C-----------------------------------------------------------------------
CL  3.1 Calculate evaporation coefficient for snow (CSB).
CL      See eqs P26.10, P26.12
C       MAX value of ASB(T) is at 243.58. This value is used for
C       T below this. The quadratic fit has a root at about 301K, when
C       there shouldn't be any snow.(The other root is at 185K)
C-----------------------------------------------------------------------
          IF(T(I).LE.243.58) THEN
            CSB=1.7405E-5
          ELSE
            CSB = ((CSB2*T(I)+CSB1)*T(I)+CSB0)
          ENDIF
          C1 = LSSW_A*(SNOW(I)*SR_RHO)**LSSW_P2
          C2 = LSSW_B*(RHO**LSSW_P3)*(SNOW(I)**LSSW_P4)
          CSB = CSB*(C1+C2)
C-----------------------------------------------------------------------
CL  3.2 Calculate implicit treatment factors alphal and bl
CL      See eqs P26.??, P26.??.
C-----------------------------------------------------------------------
          ALPHAL=ALPHF*QS(I)/(T(I)*T(I))
          BL=1. + CSB*TIMESTEP*(1. + LSRCP*ALPHAL)
C-----------------------------------------------------------------------
CL  3.3 Calculate evaporation rate, adjusted to be .LE. precip rate,
CL      in kg water per sq m per sec.  See eq P26.2.
C       Store result in QSB.  NB this is QSB*RHODZ in documentation.
C-----------------------------------------------------------------------
          QSB=MIN ( SNOW(I) , RHODZ(I)*CSB*MAX(0.0,QS(I)-Q(I))/BL )
C
C-----------------------------------------------------------------------
CL  3.4 Calculate effects of evaporation on precip rates and on Q and T.
CL     See eqs  P26.8, P26.5, P26.6 respectively.
C
C-----------------------------------------------------------------------
C
C  Increment precipitation rates (kg per sq m per sec), Q (kg per kg)
C  and T (K).  For the last 2 the change is integrated over the
C  timestep.
C
          SNOW(I)=SNOW(I)-QSB
          Q(I)=Q(I)+QSB*TIMESTEP/RHODZ(I)
          T(I)=T(I)-LSRCP*QSB*TIMESTEP/RHODZ(I)
        ENDIF ! snow>0
    3 CONTINUE
      RETURN
      END
