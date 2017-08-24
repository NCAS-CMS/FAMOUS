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
C*LL  SUBROUTINE DEWPNT-------------------------------------------------
CLL
CLL  Purpose: Calculates the 1.5 metre dewpoint from 1.5 metre specific
CLL           humidity, 1.5 metre temperature and 1.5 metre pressure.
CLL
CLL  Suitable for single column usage.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL    3.3  28/04/94 Created by Steve Woltering
CLL    4.4  Sept 97  Avoid crash if negative Q input. Damian Wilson.
CLL
CLL  Programming standard:  Unified Model Documentation Paper No 3,
CLL                         Version 5, dated 08/12/92
CLL Documentation:  To be added to UM Doc Paper ?
CLL
CLLEND-----------------------------------------------------------------
C
C*L
C*LArguments:----------------------------------------------------------
      SUBROUTINE DEWPNT(
     + Q, P, T,      ! IN
     + P_FIELD,      ! IN
     + TD            ! OUT
     +)
      IMPLICIT NONE
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

      INTEGER P_FIELD         ! IN Size of field arrays.
      REAL P(P_FIELD),        ! IN Pressure.
     +     Q(P_FIELD),        ! IN Specific humidity.
     +     T(P_FIELD)         ! IN Temperature.
      REAL RV,                ! LOCAL Gas constant for water vapour.
     +     RL1,               ! LOCAL Latent heat of evaporation.
     +     RT,                ! LOCAL.
     +     P1(P_FIELD),       ! LOCAL Pressure.
C                               j/Kg at 0 deg C.
     +     RL(P_FIELD),       ! LOCAL.
     +     Q0(P_FIELD),       ! LOCAL local SH.
     +     ES0,               ! LOCAL Saturated vapour pressure.
     +     V_PRES(P_FIELD)    ! LOCAL Vapour pressure.
      INTEGER I               ! LOCAL loop variable.
      REAL TD(P_FIELD)        ! OUT Dew point.
      PARAMETER ( RV = R / EPSILON )
      PARAMETER ( RL1 = -2.73E3 )
C*----------------------------------------------------------------------
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL  QSAT_WAT
C----------------------------------------------------------------------
C  Calculate P in HPa.
C
      DO I=1,P_FIELD
        P1(I) = P(I) / 100.0
C----------------------------------------------------------------------
C  Calculate RL - The latent heat of evaporation.
        RL(I) = LC + RL1 * ( T(I) - TM )
C----------------------------------------------------------------------
C  Calculate Vapour pressure, and from that the dewpoint in Kelvins.
        V_PRES(I) = Q(I) * P1(I) / ( EPSILON + Q(I))
      ENDDO
      CALL QSAT_WAT(Q0,T,P,P_FIELD)
      DO I=1,P_FIELD
        IF (V_PRES(I) .GT. 0.0) THEN
          ES0=(Q0(I) * P1(I)) / (EPSILON + Q0(I))
          RT = (1 / T(I)) - ( RV * ALOG(V_PRES(I)/ES0) )/RL(I)
          TD(I)=1.0/RT
          IF (TD(I) .GT. T(I)) TD(I) = T(I)
        ELSE
          TD(I)=0.0
!         print*,'WARNING. Neg or zero Q in dewpoint calc.'
        ENDIF
      ENDDO
      RETURN
      END
