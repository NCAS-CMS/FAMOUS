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

! PURPOSE: To calculate buoyancy parameters BT and BQ
!
! METHOD:
!
! ORIGINAL PROGRAMMER: J James
! CURRENT CODE OWNER: RNB Smith
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! new deck
! 4.4   29/10/97   K loops extended to give surface buoyancy parameters
!                                                             R. Essery
! 8/9/97   4.4   L_BL_LSPICE specifies mixed phase precipitation
!                scheme                    D.Wilson
!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.
!
! SYSTEM COMPONENT COVERED: ??
! SYSTEM TASK:              ??



      SUBROUTINE BOUY_TQ (
     & P_FIELD,P1
     &,P_POINTS,BL_LEVELS,P
     &,CF,Q,QCF,QCL,T
     &,TL,BT,BQ,BF
     &,L_BL_LSPICE,LTIMER
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics
     &,       L_BL_LSPICE             !IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations

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
     &,TL(P_FIELD,BL_LEVELS)  ! IN Liquid/frozen water temperature
     &,Q(P_FIELD,BL_LEVELS)   ! IN Sp humidity (kg water per kg air).
     &,QCL(P_FIELD,BL_LEVELS) ! IN Cloud liq water (kg per kg air).
     &,QCF(P_FIELD,BL_LEVELS) ! IN Cloud liq water (kg per kg air).
     &,CF(P_FIELD,BL_LEVELS)  ! IN Cloud fractions for boundary levs.


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL
     & BQ(P_FIELD,BL_LEVELS)  ! OUT A buoyancy parameter (beta q tilde)
     &,BT(P_FIELD,BL_LEVELS)  ! OUT A buoyancy parameter (beta T tilde)
     &,BF(P_FIELD,BL_LEVELS)  ! OUT A buoyancy parameter (beta F tilde)

! LOCAL VARIABLES.

      REAL
     & QSL(P_FIELD,BL_LEVELS) ! WORK Cloud liq water (kg per kg air).

      INTEGER
     &  I
     &, K

      REAL
     &  BETAT
     &, BETAQ
     &, BETAC
     &, ALPHAL
     &, AL

      EXTERNAL
     &  QSAT,QSAT_WAT,TIMER

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
!! 1.  First loop round "boundary" levels.
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS

!-----------------------------------------------------------------------
!! 1.1 Calculate saturated specific humidity at pressure and ice/liquid
!!     temperature of current level (TL).
!!     Store pressure temporarily in BQ(*,K).
!-----------------------------------------------------------------------

        IF (L_BL_LSPICE) THEN
          CALL QSAT_WAT(QSL(P1,K),TL(P1,K),P(P1,K),P_POINTS)
        ELSE
          CALL QSAT(QSL(P1,K),TL(P1,K),P(P1,K),P_POINTS)
        ENDIF
      ENDDO ! bl_levels

!-----------------------------------------------------------------------
!! 1.2 Calculate buoyancy parameters BT and BQ, required for the
!!     calculation of Richardson numbers above the surface.
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1


          BETAT = 1.0/T(I,K)                         ! P243.19 (1st eqn)
          BETAQ = C_VIRTUAL/(1.0+C_VIRTUAL*Q(I,K)-QCL(I,K)-QCF(I,K))
!                                                  ... P243.19 (2nd eqn)

          IF (TL(I,K).GT.TM.OR.L_BL_LSPICE) THEN
            ALPHAL=(EPSILON * LC * QSL(I,K))/(R * TL(I,K) * TL(I,K) )
!                 ... P243.20 (Clausius-Clapeyron) for TL above freezing

            AL = 1.0 / (1.0 + LCRCP*ALPHAL)
!                                      ... P243.21 for TL above freezing

            BETAC = CF(I,K) * AL * ( LCRCP*BETAT - ETAR*BETAQ )
!                            ... P243.19 (3rd eqn) for TL above freezing

          ELSE
            ALPHAL = (EPSILON * LS * QSL(I,K))/(R*TL(I,K) * TL(I,K))
!                 ... P243.20 (Clausius-Clapeyron) for TL below freezing

            AL = 1.0 / (1.0 + LSRCP*ALPHAL)
!                                      ... P243.21 for TL below freezing

            BETAC = CF(I,K) * AL * ( LSRCP*BETAT - ETAR*BETAQ )
!                            ... P243.19 (3rd eqn) for TL below freezing

          ENDIF
          BT(I,K) = BETAT - ALPHAL*BETAC             ! P243.18 (1st eqn)
          BQ(I,K) = BETAQ + BETAC                    ! P243.18 (2nd eqn)
          BF(I,K) = BETAQ*EPSILON*ETAR
        ENDDO !p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('BOUY_TQ ',4)
      ENDIF
      RETURN
      END
