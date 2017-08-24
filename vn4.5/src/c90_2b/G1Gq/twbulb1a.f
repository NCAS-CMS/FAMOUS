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
CLL  Subroutine TWBULB--------------------------------------------------
CLL
CLL Purpose: To calculate from temperatures,humidities the wet bulb
CLL potential temperature on model levels and output the wet bulb
CCL freezing level height.
CLL
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    10/06/94 Written by S.A.Woltering
CLL  4.0  08/09/95  Speed up by vectorising.  RTHBarnes.
!LL  4.3 26/02/97  Add first & last points to arg.list. RTHBarnes.
CLL
CLL Programming standard : UM Doc Paper no 3
CLL
CLL External documentation : UMDP no
CLL
CLLEND -----------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE TWBULB
     1 (Q,PSTAR,T,                           ! IN Model field array.
     2  AK,BK,                               ! IN AK/BK array.
     3  P_FIELD,P_LEVELS,Q_LEVELS,           ! IN Field scalars.
     4  TW                                   ! OUT WBT array.
     5  ,FIRST_POINT,LAST_POINT)             ! IN loop start & end
C
      IMPLICIT NONE
C
      INTEGER
     * P_FIELD                ! IN  No of points on a field.
     *,P_LEVELS               ! IN  No of model levels.
     *,Q_LEVELS               ! IN  No of humidity levels.
     *,FIRST_POINT,LAST_POINT ! IN 1st & last pts for calc
C
      REAL
     * T(P_FIELD,P_LEVELS)    ! IN  Intial temperature at all points.
     *,PSTAR(P_FIELD)         ! IN  Surface pressure.
     *,Q(P_FIELD,Q_LEVELS)    ! IN  Specific humidity at full levels.
     *,AK(P_LEVELS)           ! IN  Value of "A" at model level.
     *,BK(P_LEVELS)           ! IN  Value of "B" at model level.
C
      REAL
     * TW(P_FIELD,Q_LEVELS)   ! OUT The WET BULB temperature at all
C                               levels and points.
C DEFINE LOCAL WORKSPACE ARRAYS-----------------------------------------
      REAL
     * P(P_FIELD)             ! Pressure at each point.
     *,LATENT_HEAT(P_FIELD)   ! Latent heat of evaporation  (fn(T)).
     *,GG(P_FIELD)            ! The 'G' used in equation (1.3).
     *,Q_SAT(P_FIELD)         ! Saturation specific humidity
     *,DIFF(P_FIELD)          ! The difference between G(Tw) & G(Ti).
C*----------------------------------------------------------------------
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL  QSAT
C*----------------------------------------------------------------------
C   CALL COMDECKS.
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
C   DEFINE LOCAL SCALAR VARIABLES.
C-----------------------------------------------------------------------
      INTEGER
     * I,K,L          ! Loop counters.
     *,LOOP
C
      REAL
     * COEFF          ! Coeff used in latent heat calculation.
     *,CPV            ! Specific heat for water vapour.
     *,MV             ! Mol wt of water vapour.
     *,RSTAR          ! Universal gas constant.
     *,GS             ! The 'G' used in equation (1.3).
     *,DGBYDT         ! The function DG/DT  Eqn (1.6).
     *,TEMP1          !
     *,TEMP2          !
C-----------------------------------------------------------------------
C   DEFINE PARAMETER STATEMENTS
C      COEFF -  Used in the calculation of LATENT_H.
C      CPV   -  Specific heat of water vapour.
C      MV    -  Mol wt of water vapour KG/MOL.
C      RSTAR -  Universal gas constant.
C-----------------------------------------------------------------------
      PARAMETER (COEFF=2.34E3, CPV=1850.0, MV=0.01801, RSTAR=8.314)
C-----------------------------------------------------------------------
C       Begin by looping over the wet-levels.
C-----------------------------------------------------------------------
      DO L=1,Q_LEVELS         ! Loop over all wet-levels.
C-----------------------------------------------------------------------
C       Calculate pressure for all points on that wet-level.
C-----------------------------------------------------------------------
        DO I=1,P_FIELD        ! Loop over points.
          P(I) = AK(L) + BK(L)*PSTAR(I)
CL-------------------- SECTION 1 ---------------------------------------
CL      Calculate the function G for TA and QA (see eqn 3) doc no
C       G(Tw)=Qa(Ta)*L(Ta)+Ta(Cp+QaCpv)
C       Subscript a indicates initial values.
C-----------------------------------------------------------------------
          LATENT_HEAT(I)=LC-COEFF*(T(I,L)-ZERODEGC)
          GG(I)=Q(I,L)*LATENT_HEAT(I)+T(I,L)*(CP+Q(I,L)*CPV) ! Eqn 1.2
C-----------------------------------------------------------------------
          TW(I,L)=T(I,L)        ! Initialise TW.
        ENDDO                   ! End of points loop.
CL-------------------- SECTION 2 ---------------------------------------
CL      Calculate the function DG/DT
CL      G'(T)=(Mv*L(Ta)**2*Qs(T)/R*T**2) +Cp+Qa*Cpv
CL      Iterate to find the Tw
CL      T(i+1)=T(i)+(G(Tw)-G(Ti))/DG/DT(i)
C-----------------------------------------------------------------------
CL---- Set loop counter.
        LOOP=1
 1000   CONTINUE
C
        CALL QSAT(Q_SAT,TW(1,L),P(1),P_FIELD)
C
        DO I=FIRST_POINT,LAST_POINT
          GS=Q_SAT(I)*LATENT_HEAT(I)+TW(I,L)*(CP+Q(I,L)*CPV) ! Eqn 1.3

CL--------------------- SECTION 2.1 ------------------------------------
CL  Find the difference between G(Tw)-G(Ti)
C-----------------------------------------------------------------------
          DIFF(I) = ABS(GG(I)-GS)      ! Eqn 1.8
          IF (DIFF(I).GT.1.0) THEN

CL--------------------- SECTION 2.2 ------------------------------------
CL  Calculate the function DG/DT
C-----------------------------------------------------------------------
            TEMP1 = RSTAR*TW(I,L)*TW(I,L)
            TEMP2 = Q_SAT(I)*MV*LATENT_HEAT(I)*LATENT_HEAT(I)
            DGBYDT = TEMP2/TEMP1 + CP + CPV*Q(I,L)  ! Eqn 1.6

CL--------------------- SECTION 2.3 ------------------------------------
CL  Using the function DG/DT calculate an updated Temperature Tw
C-----------------------------------------------------------------------
            TW(I,L) = TW(I,L) - (GS-GG(I)) / DGBYDT  ! Eqn 1.7

CL--------------------- SECTION 2.4 ------------------------------------
CL  Using the new temperature Tw re-calculate GS First update LATENT_H
C-----------------------------------------------------------------------
            LATENT_HEAT(I)=LC-COEFF*(TW(I,L)-ZERODEGC)
          ENDIF
        END DO ! I

CL----Increment iteration loop counter.
        LOOP=LOOP+1
CL----Test for convergence.
        IF(LOOP.GT.10) THEN
          WRITE(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          WRITE(6,*)'>>> TWBULB - Convergence failure, level ',L
          WRITE(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          GOTO 9999
        ENDIF
CL----Difference test.
        DO I=FIRST_POINT,LAST_POINT
          IF (DIFF(I).GT.1.0) GOTO 1000
        ENDDO ! I
 9999   CONTINUE
!  Set points outside calculation range to sensible values.
        DO I=1,FIRST_POINT-1
          TW(I,L) = TW(FIRST_POINT,L)
        END DO
        DO I=LAST_POINT+1,P_FIELD
          TW(I,L) = TW(LAST_POINT,L)
        END DO
      ENDDO                     ! End of loop over levels.
      RETURN
      END
