C *****************************COPYRIGHT******************************
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
!
! SUBROUTINE BOUY_TQ

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels
!
! METHOD:
!
! ORIGINAL PROGRAMMER: J. James
! CURRENT CODE OWNER: R.N.B. Smith
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! new deck
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.
!



      SUBROUTINE BOUY_TQ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,T,Q,QCF,QCL
     &,BT,BQ,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,LTIMER
     & )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field.
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be processed.
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed  <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA

      REAL
     & P(P_FIELD,BL_LEVELS)   ! IN Pressure at pressure points.
     &,T(P_FIELD,BL_LEVELS)   ! IN Temperature (K). At P points
     &,Q(P_FIELD,BL_LEVELS)   ! IN Sp humidity (kg water per kg air).
     &,QCL(P_FIELD,BL_LEVELS) ! IN Cloud liq water (kg per kg air).
     &,QCF(P_FIELD,BL_LEVELS) ! IN Cloud liq water (kg per kg air).


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL
     & BQ(P_FIELD,BL_LEVELS)  ! OUT A buoyancy parameter for clear air
     &,BT(P_FIELD,BL_LEVELS)  ! OUT A buoyancy parameter for clear air
     &,BQ_CLD(P_FIELD,BL_LEVELS)
!                             ! OUT A buoyancy parameter for cloudy air
     &,BT_CLD(P_FIELD,BL_LEVELS)
!                             ! OUT A buoyancy parameter for cloudy air
     &,A_QS(P_FIELD,BL_LEVELS)
!                             ! OUT Saturated lapse rate factor
     &,A_DQSDT(P_FIELD,BL_LEVELS)
!                             ! OUT Saturated lapse rate factor
     &,DQSDT(P_FIELD,BL_LEVELS)
!                             ! OUT Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.

      REAL
     & QS(P_FIELD)            ! WORK Saturated mixing ratio.

      INTEGER
     &  I
     &, K

      REAL
     &  BC

      EXTERNAL
     &  QSAT,TIMER

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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------


      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     & ETAR=1.0/(1.0-EPSILON)   ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                ! Used in DZTL, FTL calculations.
     &,LCRCP=LC/CP              ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP              ! Latent heat of fusion / CP.
     &,LS=LC+LF                 ! Latent heat of sublimation.
     &,LSRCP=LS/CP              ! Latent heat of sublimation / CP.
     &)

      IF (LTIMER) THEN
        CALL TIMER('BOUY_TQ ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round levels.
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
!-----------------------------------------------------------------------
!! 1.1 Calculate saturated specific humidity at pressure and
!!     temperature of current level.
!-----------------------------------------------------------------------
        CALL QSAT(QS(P1),T(P1,K),P(P1,K),P_POINTS)
!
        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 1.2 Calculate buoyancy parameters BT and BQ, required for the
!!     calculation of stability.
!-----------------------------------------------------------------------

          BT(I,K) = 1.0/T(I,K)
          BQ(I,K) = C_VIRTUAL/(1.0+C_VIRTUAL*Q(I,K)-QCL(I,K)-QCF(I,K))
!
          IF (T(I,K) .GT. TM) THEN
            DQSDT(I,K) = (EPSILON * LC * QS(I))
     &                   / ( R * T(I,K) * T(I,K) )
!                      ...  (Clausius-Clapeyron) for T above freezing
!
            A_QS(I,K) = 1.0 / (1.0 + LCRCP*DQSDT(I,K))
!
            A_DQSDT(I,K) = A_QS(I,K) * DQSDT(I,K)
!
            BC = LCRCP*BT(I,K) - ETAR*BQ(I,K)
!
          ELSE
            DQSDT(I,K) = (EPSILON * LS * QS(I))
     &                   / ( R * T(I,K) * T(I,K) )
!                      ...  (Clausius-Clapeyron) for T below freezing
!
            A_QS(I,K) = 1.0 / (1.0 + LSRCP*DQSDT(I,K))
!
            A_DQSDT(I,K) = A_QS(I,K) * DQSDT(I,K)
!
            BC = LSRCP*BT(I,K) - ETAR*BQ(I,K)
!
          ENDIF
!
!-----------------------------------------------------------------------
!! 1.3 Calculate in-cloud buoyancy parameters.
!-----------------------------------------------------------------------
!
          BT_CLD(I,K) = BT(I,K) - A_DQSDT(I,K) * BC
          BQ_CLD(I,K) = BQ(I,K) + A_QS(I,K) * BC
!
        ENDDO ! p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('BOUY_TQ ',4)
      ENDIF
      RETURN
      END
