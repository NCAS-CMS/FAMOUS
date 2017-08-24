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

!-----------------------------------------------------------------------
!
! Subroutines SF_FLUX_LAND and SF_FLUX_SEA to calculate explicit surface
! fluxes of heat and moisture
!
!-----------------------------------------------------------------------

!     SUBROUTINE SF_FLUX_LAND-------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over
!     land tiles
!
!     ------------------------------------------------------------------
      SUBROUTINE SF_FLUX_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,
     & ASHTF,LH,QS1,QSTAR,QW_1,RADNET,RESFT,RHOKH_1,TILE_FRAC,
     & TL_1,TS1,TSTAR,Z0H,Z0M_EFF,Z1_TQ,
     & FQW_1_GB,FTL_1_GB,
     & ALPHA1,FQW_1,FTL_1,RHOKPM,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD             ! IN Total number of P-grid points.
     &,LAND_FIELD          ! IN Total number of land points.
     &,TILE_PTS            ! IN Number of tile points.
     &,LAND_INDEX(P_FIELD) ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD)! IN Index of tile points.

      LOGICAL
     & LTIMER              ! IN Logical for TIMER

      REAL
     & ASHTF(P_FIELD)      ! IN Coefficient to calculate surface
!                          !    heat flux into soil (W/m2/K).
     &,LH                  ! IN Latent heat (J/K/kg).
     &,QS1(P_FIELD)        ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSTAR(LAND_FIELD)   ! IN Surface qsat.
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
     &,RADNET(P_FIELD)     ! IN Net surface radiation (W/m2) positive
!                          !    downwards
     &,RESFT(LAND_FIELD)   ! IN Total resistance factor.
     &,RHOKH_1(LAND_FIELD) ! IN Surface exchange coefficient.
     &,TILE_FRAC(LAND_FIELD)
!                          ! IN Tile fraction.
     &,TL_1(P_FIELD)       ! IN Liquid/frozen water temperature for
!                          !    lowest atmospheric layer (K).
     &,TS1(LAND_FIELD)     ! IN Temperature of surface layer (K).
     &,TSTAR(LAND_FIELD)   ! IN Surface temperature (K).
     &,Z0H(LAND_FIELD)     ! IN Roughness length for heat and moisture
     &,Z0M_EFF(LAND_FIELD) ! IN Effective roughness length for momentum
     &,Z1_TQ(P_FIELD)      ! IN Height of lowest atmospheric level (m).

      REAL
     & FQW_1_GB(P_FIELD)   ! INOUT GBM surface flux of QW (kg/m2/s).
     &,FTL_1_GB(P_FIELD)   ! INOUT GBM surface flux of TL.

      REAL
     & ALPHA1(LAND_FIELD)  ! OUT Gradient of saturated specific humidity
!                          !     with respect to temperature between the
!                          !     bottom model layer and the surface.
     &,FQW_1(LAND_FIELD)   ! OUT Local surface flux of QW (kg/m2/s).
     &,FTL_1(LAND_FIELD)   ! OUT Local surface flux of TL.
     &,RHOKPM(LAND_FIELD)  ! OUT Modified surface exchange coefficient.

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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


! Derived local parameters
! Derived local parameters
      REAL GRCP,LS
      PARAMETER (
     & GRCP=G/CP
     &,LS=LF+LC            ! Latent heat of sublimation.
     & )

! Scalars
      INTEGER
     & I                   ! Horizontal field index.
     &,J                   ! Tile field index.
     &,L                   ! Land point field index.

      REAL
     & DQ1                 ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
     &,D_T                 ! Temporary in calculation of alpha1.
     &,RAD_REDUC           ! Radiation term required for surface flux

      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        D_T = TSTAR(L) - TL_1(I)
        IF (D_T .GT. 0.05 .OR. D_T .LT. -0.05) THEN
          ALPHA1(L) = (QSTAR(L) - QS1(I)) / D_T
        ELSEIF (TL_1(I) .GT. TM) THEN
          ALPHA1(L) = EPSILON*LC*QS1(I)*(1. + C_VIRTUAL*QS1(I))
     &                                            / ( R*TL_1(I)*TL_1(I))
        ELSE
          ALPHA1(L) = EPSILON*LS*QS1(I)*(1. + C_VIRTUAL*QS1(I))
     &                                            / ( R*TL_1(I)*TL_1(I))
        ENDIF
      ENDDO

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)

        RHOKPM(L) = RHOKH_1(L) / ( ASHTF(I)  +
     &                         RHOKH_1(L)*(LH*ALPHA1(L)*RESFT(L) + CP) )
        RAD_REDUC = RADNET(I) - ASHTF(I) * ( TL_1(I) - TS1(L)
     &                         + GRCP*(Z1_TQ(I) + Z0M_EFF(L) - Z0H(L)) )
        DQ1 = QS1(I) - QW_1(I) +
     &                   GRCP*ALPHA1(L)*(Z1_TQ(I) + Z0M_EFF(L) - Z0H(L))
        FQW_1(L) = RESFT(L)*RHOKPM(L)*( ALPHA1(L)*RAD_REDUC
     &                                + (CP*RHOKH_1(L) + ASHTF(I))*DQ1 )
        FTL_1(L) = RHOKPM(L)*(RAD_REDUC - LH*RESFT(L)*RHOKH_1(L)*DQ1)

        FTL_1_GB(I) = FTL_1_GB(I) + TILE_FRAC(L)*FTL_1(L)
        FQW_1_GB(I) = FQW_1_GB(I) + TILE_FRAC(L)*FQW_1(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',4)
      ENDIF

      RETURN
      END

!     SUBROUTINE SF_FLUX_SEA--------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over sea
!     and sea-ice
!
!     ------------------------------------------------------------------
      SUBROUTINE SF_FLUX_SEA (
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,LAND_MASK,
     & ASHTF,ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1,TI,
     & TL_1,TSTAR_ICE,TSTAR_SEA,Z0H_ICE,Z0M_ICE,Z0H_SEA,Z0M_SEA,Z1_TQ,
     & ALPHA1,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,RHOKPM,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,P_FIELD             ! IN Total number of P-grid points.
     &,P1                  ! IN First P-point to be processed.
     &,NSICE               ! IN Number of sea-ice points.
     &,SICE_INDEX(P_FIELD) ! IN Index of sea-ice points

      LOGICAL
     & LTIMER              ! IN  Logical for TIMER
     &,LAND_MASK(P_FIELD)  ! IN .TRUE. for land, .FALSE. elsewhere.

      REAL
     & ASHTF(P_FIELD)      ! IN Coefficient to calculate surface
!                          !    heat flux into sea-ice (W/m2/K).
     &,ICE_FRACT(P_FIELD)  ! IN Fraction of gridbox which is sea-ice.
     &,QS1(P_FIELD)        ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSTAR_ICE(P_FIELD)  ! IN Surface qsat for sea-ice.
     &,QSTAR_SEA(P_FIELD)  ! IN Surface qsat for sea or sea-ice leads.
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
     &,RADNET(P_FIELD)     ! IN Net surface radiation (W/m2) positive
!                          !    downwards
     &,RHOKH_1(P_FIELD)    ! IN Surface exchange coefficient.
     &,TI(P_FIELD)         ! IN Temperature of sea-ice surface layer (K)
     &,TL_1(P_FIELD)       ! IN Liquid/frozen water temperature for
!                          !    lowest atmospheric layer (K).
     &,TSTAR_ICE(P_FIELD)  ! IN Sea-ice surface temperature (K).
     &,TSTAR_SEA(P_FIELD)  ! IN Sea surface temperature (K).
     &,Z0H_ICE(P_FIELD)    ! IN Sea-ice heat and moisture roughness
!                          !    length (m).
     &,Z0M_ICE(P_FIELD)    ! IN Sea-ice momentum roughness length (m).
     &,Z0H_SEA(P_FIELD)    ! IN Sea and lead heat and moisture roughness
!                          !    length (m).
     &,Z0M_SEA(P_FIELD)    ! IN Sea and lead momentum roughness length.
     &,Z1_TQ(P_FIELD)      ! IN Height of lowest atmospheric level (m).

      REAL
     & ALPHA1(P_FIELD)     ! OUT Gradient of saturated specific humidity
!                          !     with respect to temperature between the
!                          !     bottom model layer and the surface.
     &,E_SEA(P_FIELD)      ! OUT Evaporation from sea times leads
!                          !     fraction (kg/m2/s).
     &,FQW_ICE(P_FIELD)    ! OUT Surface flux of QW for sea-ice.
     &,FQW_1(P_FIELD)      ! OUT GBM surface flux of QW (kg/m2/s).
     &,FTL_ICE(P_FIELD)    ! OUT Surface flux of TL for sea-ice.
     &,FTL_1(P_FIELD)      ! OUT GBM surface flux of TL.
     &,H_SEA(P_FIELD)      ! OUT Surface sensible heat flux over sea
!                          !     times leads fraction (W/m2).
     &,RHOKPM(P_FIELD)     ! OUT Modified surface exchange coefficient.

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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


! Derived local parameters.
      REAL GRCP,LS
      PARAMETER (
     & GRCP=G/CP
     &,LS=LF+LC            ! Latent heat of sublimation.
     & )

! Scalars
      INTEGER
     & I                   ! Horizontal field index.
     &,J                   ! Sea-ice field index.
      REAL
     & DQ1                 ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
     &,D_T                 ! Temporary in calculation of alpha1.
     &,RAD_REDUC           ! Radiation term required for surface flux
!                          ! calcs.

      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',3)
      ENDIF

      DO I=P1,P1+P_POINTS-1
        ALPHA1(I) = 0.
        E_SEA(I) = 0.
        H_SEA(I) = 0.
        FQW_ICE(I) = 0.
        FTL_ICE(I) = 0.
        RHOKPM(I) = 0.
      ENDDO

!----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes - only required for sea-ice points
!----------------------------------------------------------------------
      DO J=1,NSICE
        I = SICE_INDEX(J)
        D_T = TSTAR_ICE(I) - TL_1(I)
        IF (D_T .GT. 0.05 .OR. D_T .LT. -0.05) THEN
          ALPHA1(I) = (QSTAR_ICE(I) - QS1(I)) / D_T
        ELSEIF (TL_1(I) .GT. TM) THEN
          ALPHA1(I) = EPSILON*LC*QS1(I)*( 1.0+C_VIRTUAL*QS1(I) )
     &                                             /(R*TL_1(I)*TL_1(I))
        ELSE
          ALPHA1(I) = EPSILON*LS*QS1(I)*( 1.0+C_VIRTUAL*QS1(I) )
     &                                             /(R*TL_1(I)*TL_1(I))
        ENDIF
      ENDDO

      DO I=P1,P1+P_POINTS-1
        IF ( .NOT. LAND_MASK(I) ) THEN

          E_SEA(I) = - (1. - ICE_FRACT(I)) *
     &                               RHOKH_1(I)*(QW_1(I) - QSTAR_SEA(I))
          H_SEA(I) = - (1. - ICE_FRACT(I))*CP*RHOKH_1(I) *
     &                 ( TL_1(I) - TSTAR_SEA(I)
     &                     + GRCP*(Z1_TQ(I) + Z0M_SEA(I) - Z0H_SEA(I)) )

          IF ( ICE_FRACT(I) .GT. 0. ) THEN
! Sea-ice
            RHOKPM(I) = RHOKH_1(I) / ( ASHTF(I) +
     &                                  RHOKH_1(I)*(LS*ALPHA1(I) + CP) )
            RAD_REDUC = RADNET(I) - ICE_FRACT(I) * ASHTF(I) *
     &                  ( TL_1(I) - TI(I) +
     &                       GRCP*(Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I)) )
            DQ1 = QS1(I) - QW_1(I) +
     &               GRCP*ALPHA1(I)*(Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I))
            FQW_ICE(I) = RHOKPM(I) * ( ALPHA1(I)*RAD_REDUC +
     &                     (CP*RHOKH_1(I) + ASHTF(I))*DQ1*ICE_FRACT(I) )
            FTL_ICE(I) = RHOKPM(I) * ( RAD_REDUC -
     &                                  ICE_FRACT(I)*LS*RHOKH_1(I)*DQ1 )

          ENDIF

          FTL_1(I) = FTL_ICE(I) + H_SEA(I) / CP
          FQW_1(I) = FQW_ICE(I) + E_SEA(I)

        ENDIF
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',4)
      ENDIF

      RETURN
      END

