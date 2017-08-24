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
CLL  Subroutine THETAW--------------------------------------------------
CLL
CLL Purpose: To calculate from temperatures,humidities the wet bulb
CLL potential temperature.
CLL
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
!LL   4.5    15/04/98 Added START,END arguments to enable duplicate
!LL                   halo calculations to be avoided. S.D.Mullerworth
CLL
CLL Programming standard : UM Doc Paper no 3
CLL
CLL Logical components covered : D482
CLL
CLL Project task :  D482
CLL
CLL External documentation : UMDP no
CLL
CLLEND -----------------------------------------------------------------
C
C*L ARGUMENTS:-----------------------------------------------------
      SUBROUTINE THETAW
     1 (PRESS_REQD,THETA,Q,P,PSTAR,P_EXNER_HALF,TW,
     2  AK,BK,AKH,BKH,P_FIELD,P_LEVELS,Q_LEVELS,T_REF,
     &  START,END,
     3  CMESSAGE,ICODE)
C
      IMPLICIT NONE
C
      INTEGER
     * P_FIELD        ! IN  No of points on a field
     *,P_LEVELS       ! IN  NO OF MODEL LEVELS
     *,Q_LEVELS       ! IN  NO OF HUMIDITY LEVELS
     *,T_REF          ! IN  LEVEL OF MODEL USED TO CALC HEIGHT
     &,START,END      ! IN  Range of points to calculate
     *,ICODE          ! OUT RETURN CODE
C
      CHARACTER*(80) CMESSAGE
C
      REAL
     * P(P_FIELD,P_LEVELS)    ! IN  Pressure at each level and point
     *,THETA(P_FIELD,P_LEVELS)! IN  Potential temperature on all pts
     *,P_EXNER_HALF(P_FIELD,P_LEVELS+1) ! IN EXNER Pressure at model
     *                        !     half levels
     *,PSTAR(P_FIELD)         ! IN  Surface pressure
     *,Q(P_FIELD,Q_LEVELS)    ! IN  Specific humidity at full levels
     *,MODEL_HALF_HEIGHT(P_FIELD,P_LEVELS+1)! IN Height of 1/2 levels
     *,AKH(P_LEVELS+1)        ! IN  Value of "A" at mid layer
     *,BKH(P_LEVELS+1)        ! IN  Value of "B" at mid layer
     *,AK (P_LEVELS)          ! IN  Value of "A" at model level
     *,BK (P_LEVELS)          ! IN  Value of "B" at model level
     *,PRESS_REQD             ! IN  Pressure level required for THETAW
C
      REAL
     * TW(P_FIELD)            ! OUT The WET BULB temperature. On
C ! output the TW at 1000mb ie thetaw
C
C*---------------------------------------------------------------
C
C*L WORKSPACE USAGE----------------------------------------------
C     10*P_FIELD
C DEFINE LOCAL WORKSPACE ARRAYS
      REAL
     * T(P_FIELD)             ! Intial temperature at Press_reqd (Ta)
     *,TW_TEMP(P_FIELD)       ! A temporary store for temperature.
     *,QAD(P_FIELD)           ! SPECIFIC HUMIDITY AT PRESS_REQD
     *,Q_SAT(P_FIELD)         ! Saturation specific humidityEQD
     *,PRESS_REQD_HORIZ(P_FIELD) ! PRESS_REQD BROADCAST OVER DOMAIN
     *,LATENT_HEAT(P_FIELD)    ! Latent heat of evaporation  (fn(T))
     *,GG(P_FIELD)             ! The 'G' used in equation (1.3)
     *,GS                      ! The 'G' used in equation (1.3)
     *,DIFF(P_FIELD)           ! The difference between G(Tw) & G(Ti)
     *,DGBYDT                  ! The function DG/DT  Eqn (1.6)
     *,FF(P_FIELD)             ! The function DT/DP  Eqn (1.5)
     *,FF_NEXT(P_FIELD)        ! The function DT/DP
C*---------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED----------------------------------
      EXTERNAL  V_INT,V_INT_T,QSAT
C*---------------------------------------------------------------
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
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C
C----------------------------------------------------------------
C   DEFINE LOCAL VARIABLES
C----------------------------------------------------------------
      INTEGER
     * I,K            !  LOOP COUNTERS
     *,J              !  LOWEST ETA LEVEL WITH TEMP<T0
     *,LOOP
     *,NLOOP
     *,LOOP_FIELD(P_FIELD)
C
      REAL
     * DUMMY1
     *,DUMMY2
     *,COEFF          !  COEFF USED IN THE CALCULATION OF LATENT HEAT
     *,CPV            !  SPECIFIC HEAT FOR WATER VAPOUR
     *,MV             !  Mol wt of water vapour
     *,MD             !  Mol wt of dry air
     *,RSTAR          !  Universal gas constant
     *,TEMP1          !
     *,TEMP2          !
     *,PP             !
     *,DELTAP         !  Pressure increment
     *,P_TARGET       !  The target pressure (1000mb)
     *,REAL_NLOOP     !  Real NLOOP
     *,REM            !  Temp local variable

      DATA COEFF/2.34E3/      ! Used in the calculation of LATENT_H
      DATA CPV/1850.0/        ! Specific heat of water vapour
      DATA MV/0.01801/        ! Mol wt of water vapour KG/MOL
      DATA MD/0.02896/        ! Mol wt of dry air      KG/MOL
      DATA RSTAR/8.314/       ! Universal gas constant
C
C
C------------------------------------------------------------------
CL  1. Calculate the temperature and specific hum at the PRESS_REQ
C------------------------------------------------------------------

      DO 1 I=START,END
      PRESS_REQD_HORIZ(I)=PRESS_REQD
    1 CONTINUE

        CALL V_INT_T(T,PRESS_REQD_HORIZ,P(1,T_REF),PSTAR,
     &  P_EXNER_HALF,THETA,P_FIELD,P_LEVELS,T_REF,AKH,BKH
     &  ,START,END)


        CALL V_INT(P,PRESS_REQD_HORIZ,Q,QAD,P_FIELD,Q_LEVELS,
     &  DUMMY1,DUMMY2,.FALSE.,START,END)


C---------------------------------------------------------------------
CL---------------------SECTION 2 -------------------------------------
CL      Calculate the function G for TA and QA (see eqn 3) doc no
C       G(Tw)=Qa(Ta)*L(Ta)+Ta(Cp+QaCpv)
C       Subscript a indicates initial values
C---------------------------------------------------------------------


      DO 20 I=START,END
      LATENT_HEAT(I)=LC-COEFF*(T(I)-ZERODEGC)
      GG(I)=QAD(I)*LATENT_HEAT(I)+T(I)*(CP+QAD(I)*CPV)  ! Eqn 1.2
  20  CONTINUE

C---------------------------------------------------------------------
CL---------------------SECTION 3 -------------------------------------
CL      Calculate the function DG/DT  (see eqn 5) doc no
CL      G'(T)=(Mv*L(Ta)**2*Qs(T)/R*T**2) +Cp+Qa*Cpv
CL      Iterate to find the Tw using 6
CL      T(i+1)=T(i)+(G(Tw)-G(Ti))/DG/DT(i)
C---------------------------------------------------------------------

      DO 30 I=START,END
      TW(I)=T(I)
 30   CONTINUE

      LOOP=1
10000 CONTINUE

      CALL QSAT(Q_SAT(START),TW(START),PRESS_REQD_HORIZ(START)
     &  ,END-START+1)           !Qs @ required Ta & P
C
      DO 31 I=START,END
      GS=Q_SAT(I)*LATENT_HEAT(I)+TW(I)*(CP+QAD(I)*CPV)   ! Eqn 1.3

CL---------------------SECTION 3.1------------------------------------
CL  Find the difference between G(Tw)-G(Ti)
C---------------------------------------------------------------------

      DIFF(I)=ABS(GG(I)-GS)      ! Eqn 1.8
      IF(DIFF(I).GT.1.0) THEN

CL---------------------SECTION 3.1------------------------------------
CL  Calculate the function DG/DT
C---------------------------------------------------------------------

        TEMP1=RSTAR*TW(I)**2
        TEMP2=Q_SAT(I)*MV*LATENT_HEAT(I)**2
        DGBYDT=TEMP2/TEMP1 + CP + CPV*QAD(I)  ! Eqn 1.6

CL---------------------SECTION 3.2------------------------------------
CL  Using the function DG/DT calculate an updated Temperature Tw
C---------------------------------------------------------------------

        TW(I) = TW(I) - (GS-GG(I)) / DGBYDT  ! Eqn 1.7

CL---------------------SECTION 3.3------------------------------------
CL  Using the new temperature Tw re-calculate GS First update LATENT_H
C---------------------------------------------------------------------

        LATENT_HEAT(I)=LC-COEFF*(TW(I)-ZERODEGC)

      ENDIF
  31  CONTINUE

      LOOP=LOOP+1

      IF(LOOP.GT.10) THEN
        ICODE=1
        CMESSAGE='THETAW  Convergence failure in THETAW'
        GOTO 9999
      ENDIF

      DO 32 I=START,END
      IF(DIFF(I).GT.1.0) GOTO 10000
   32 CONTINUE

C---------------------------------------------------------------------
CL---------------------SECTION 4 -------------------------------------
CL      Having calculated now have to descend the sat adiabat down
CL      to 1000mb solving:-
CL      DT=(L*EPS/P*Es +R*T/Md)/(Cp*P+(EPS*Mv*L**2*Es/(R*T**2)) *DP
C---------------------------------------------------------------------

      PP=PRESS_REQD
      DELTAP=5000.0     ! Pressure increment in PASCALS
      P_TARGET=100000.0  ! 1000 mb
      REAL_NLOOP=(P_TARGET-PP)/DELTAP + 0.000001
      NLOOP=REAL_NLOOP
      REM=REAL_NLOOP-NLOOP
      IF(REM.GT.0.00001) THEN
        NLOOP=NLOOP+1
        DELTAP=(P_TARGET-PP)/NLOOP
      ENDIF

CL Calculate DT/DP = F(T,P)

      DO 41 K=1,NLOOP
      DO 40 I=START,END
      TEMP1=(Q_SAT(I)*PP*MV*LATENT_HEAT(I)**2/(RSTAR*TW(I)**2)) + CP*PP
      TEMP2=Q_SAT(I)*LATENT_HEAT(I)+RSTAR*TW(I)/MD
      FF(I)=TEMP2/TEMP1         !    Eqn 1.5
      TW_TEMP(I)=FF(I)*DELTAP + TW(I)  ! Eqn 1.10
  40  CONTINUE

      PP=PP+DELTAP
      DO 42 I=START,END
      PRESS_REQD_HORIZ(I)=PP
      LATENT_HEAT(I)=LC-COEFF*(TW_TEMP(I)-ZERODEGC)
   42 CONTINUE

      CALL QSAT(Q_SAT(START),TW_TEMP(START),PRESS_REQD_HORIZ(START)
     &  ,END-START+1)           !Qs @  Tw & P

      DO 43 I=START,END
      TEMP1=(Q_SAT(I)*PP*MV*LATENT_HEAT(I)**2/(RSTAR*TW_TEMP(I)**2))
     * + CP*PP
      TEMP2=Q_SAT(I)*LATENT_HEAT(I)+RSTAR*TW_TEMP(I)/MD
      FF_NEXT(I)=TEMP2/TEMP1
      TW(I)=(FF(I)+FF_NEXT(I))*0.5*DELTAP + TW(I)  ! Eqn 1.11
      LATENT_HEAT(I)=LC-COEFF*(TW(I)-ZERODEGC)
  43  CONTINUE
  41  CONTINUE

 9999 CONTINUE
      RETURN
      END
C
C
