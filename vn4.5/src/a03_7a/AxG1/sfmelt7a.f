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
C SUBROUTINE SF_MELT----------------------------------------------------
C Purpose : Calculates surface melting (snow and sea-ice) and increments
C           surface fluxes to satisfy energy balance.
C           Sub-surface snowmelt is calculated and snowdepth incremented
C           by melt and sublimation in P251.
C           R.Essery 19/1/95
C-----------------------------------------------------------------------
      SUBROUTINE SF_MELT (
     & POINTS,P_FIELD,P1,LAND_FIELD,LAND_INDEX
     &,SNOW_INDEX,NSNOW,LAND_MASK,LTIMER,SIMLT,SMLT
     &,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_SNOW,DTRDZ_1,ICE_FRACT
     &,LYING_SNOW,RHOKH_1,RHOKH_1_SICE,SNOW_FRAC,TIMESTEP
     &,FQW_1,FQW_ICE,FQW_SNOW,FTL_1,FTL_SNOW,QW_1
     &,TL_1,TSTAR,TSTAR_SNOW,TI
     &,EI,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT
     & )

      IMPLICIT NONE

      INTEGER
     & POINTS               ! IN Number of P-grid points to be
!                           !    processed.
     &,P_FIELD              ! IN Total number of P-grid points.
     &,P1                   ! IN First P-point to be processed.
     &,LAND_FIELD           ! IN Total number of land points..
     &,LAND_INDEX(P_FIELD)  !IN Index of land points.
     &,SNOW_INDEX(LAND_FIELD)!IN Index of snow points.
     &,NSNOW                !IN Number of snow points.

      LOGICAL
     & LAND_MASK(P_FIELD)   ! IN T for land points, F otherwise.
     &,LTIMER               ! IN Logical for TIMER.
     &,SIMLT                ! IN STASH flag for sea-ice melting ht flux.
     &,SMLT                 ! IN STASH flag for snow melting ht flux.

       REAL
     & ALPHA1(LAND_FIELD)   ! IN Gradient of saturated specific
!                           !    humidity with respect to temp.
!                           !    between the bottom model layer
!                           !    and the snow surface.
     &,ALPHA1_SICE(P_FIELD) ! IN ALPHA1 for sea-ice.
     &,ASHTF(P_FIELD)       ! IN Coefficient to calculate surface
!                           !    heat flux into sea-ice (W/m2/K).
     &,ASHTF_SNOW(P_FIELD)  ! IN Coefficient to calculate surface
!                           !    heat flux into snow (W/m2/K).
     &,DTRDZ_1(P_FIELD)     ! IN -g.dt/dp for surface layer
     &,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is covered
!                           !    by sea-ice.
     &,LYING_SNOW(P_FIELD)  ! IN Lying snow (kg/m2).
     &,RHOKH_1(LAND_FIELD)  ! IN Surface exchange coefficient for snow.
     &,RHOKH_1_SICE(P_FIELD)! IN Surface exchange coefficient for
!                           !    sea-ice.
     &,SNOW_FRAC(LAND_FIELD)! IN Fraction of gridbox which is covered
!                           !    by snow.
     &,TIMESTEP             ! IN Timestep (sec).

      REAL
     & FQW_1(P_FIELD)       ! INOUT GBM surface moisture flux (kg/m2/s).
     &,FQW_ICE(P_FIELD)     ! INOUT FQW for sea-ice.
     &,FQW_SNOW(LAND_FIELD) ! INOUT FQW for snow.
     &,FTL_1(P_FIELD)       ! INOUT GBM surface sens. heat flux (W/m2).
     &,FTL_SNOW(LAND_FIELD) ! INOUT FTL for snow.
     &,QW_1(P_FIELD)        ! INOUT Total water content of lowest
!                           !       atmospheric layer (kg per kg air).
     &,TL_1(P_FIELD)        ! INOUT Liquid/frozen water temperature for
!                           !       lowest atmospheric layer (K).
     &,TSTAR(P_FIELD)       ! INOUT GBM surface temperature (K).
     &,TSTAR_SNOW(LAND_FIELD)!INOUT Snow surface temperature (K).
     &,TI(P_FIELD)          ! INOUT Sea-ice surface layer temp. (K).

      REAL
     & EI(P_FIELD)          ! OUT Sublimation from lying snow or
!                           !     sea-ice (kg/m2/s).
     &,SICE_MLT_HTF(P_FIELD)! OUT Heat flux due to melting of sea-ice
!                           !     (W/m2).
     &,SNOMLT_SURF_HTF(P_FIELD)
!                           ! OUT Heat flux due to surface melting
!                           !     of snow (W/m2).
     &,SNOWMELT(P_FIELD)    ! OUT Surface snowmelt (kg/m2/s).

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
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

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------

      REAL
     & DFQW                 ! Moisture flux increment.
     &,DFTL                 ! Sensible heat flux increment.
     &,DTSTAR               ! Surface temperature increment.
     &,LCMELT               ! Temporary in melt calculations.
     &,LSMELT               ! Temporary in melt calculations.
     &,RHOKH1_PRIME         ! Modified forward time-weighted
!                           ! transfer coefficient.
     &,SNOW_MAX             ! Snow available for melting.
     &,TSTARMAX             ! Maximum gridbox mean surface temperature
!                           ! at sea points with ice.
      INTEGER
     & I                    ! Loop counter - full horizontal field.
     &,J                    !
     &,L                    ! Loop counter - land field.
C
      IF (LTIMER) THEN
      CALL TIMER('SFMELT  ',3)
      ENDIF

      DO I=P1,P1+POINTS-1
        IF (SIMLT) SICE_MLT_HTF(I) = 0.0
        IF (SMLT) SNOMLT_SURF_HTF(I) = 0.0
        SNOWMELT(I) = 0.0
        EI(I) = 0.0
      ENDDO

      DO J=1,NSNOW
        L = SNOW_INDEX(J)
        I = LAND_INDEX(L)
!-----------------------------------------------------------------------
!  Melt snow if TSTAR_SNOW is greater than TM.
!-----------------------------------------------------------------------
        EI(I) = SNOW_FRAC(L)*FQW_SNOW(L)
        SNOW_MAX = MAX( 0.0,
     &               LYING_SNOW(I)/SNOW_FRAC(L) - FQW_SNOW(L)*TIMESTEP )
        IF ( SNOW_MAX.GT.0.0 .AND. TSTAR_SNOW(L).GT.TM ) THEN
          RHOKH1_PRIME = 1. / ( 1. / RHOKH_1(L) +
     &                                GAMMA(1)*SNOW_FRAC(L)*DTRDZ_1(I) )
          LCMELT = (CP + LC*ALPHA1(L))*RHOKH1_PRIME + ASHTF_SNOW(I)
          LSMELT = (CP + (LC+LF)*ALPHA1(L))*RHOKH1_PRIME + ASHTF_SNOW(I)
          SNOWMELT(I) = LSMELT * MIN( (TSTAR_SNOW(L) - TM)/LF ,
     &                                      SNOW_MAX/(LCMELT*TIMESTEP) )
          DFTL = - CP*RHOKH1_PRIME*LF*SNOWMELT(I) / LSMELT
          DFQW = - ALPHA1(L)*RHOKH1_PRIME*LF*SNOWMELT(I) / LSMELT
          FTL_SNOW(L) = FTL_SNOW(L) + DFTL
          FQW_SNOW(L) = FQW_SNOW(L) + DFQW
          TSTAR_SNOW(L) = TM
!-----------------------------------------------------------------------
!  Update gridbox-mean quantities
!-----------------------------------------------------------------------
          SNOWMELT(I) = SNOW_FRAC(L)*SNOWMELT(I)
          IF (SMLT) SNOMLT_SURF_HTF(I) = LF*SNOWMELT(I)
          DFTL = SNOW_FRAC(L)*DFTL
          DFQW = SNOW_FRAC(L)*DFQW
          TL_1(I) = TL_1(I) + DTRDZ_1(I) * DFTL / CP
          QW_1(I) = QW_1(I) + DTRDZ_1(I) * DFQW
          FTL_1(I) = FTL_1(I) + DFTL
          FQW_1(I) = FQW_1(I) + DFQW
          EI(I) = EI(I) + DFQW
        ENDIF
      ENDDO

      DO I=P1,P1+POINTS-1
        IF ( .NOT.LAND_MASK(I) .AND. ICE_FRACT(I).GT.0.0 ) THEN
!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR > TSTARMAX or TI > TM.
!-----------------------------------------------------------------------
          EI(I) = FQW_ICE(I)
          TSTARMAX = ICE_FRACT(I)*TM + (1.0 - ICE_FRACT(I))*TFS
          IF ( TSTAR(I) .GT. TSTARMAX ) THEN
            RHOKH1_PRIME = 1. / ( 1. / RHOKH_1_SICE(I)
     &                              + ICE_FRACT(I)*GAMMA(1)*DTRDZ_1(I) )
            DTSTAR = TSTARMAX - TSTAR(I)
            LSMELT = (CP + (LC + LF)*ALPHA1_SICE(I))*RHOKH1_PRIME
     &                                                        + ASHTF(I)
            DFTL = CP * RHOKH1_PRIME * DTSTAR
            DFQW = ALPHA1_SICE(I) * RHOKH1_PRIME * DTSTAR
            TI(I) =TI(I) + AI*ASHTF(I)*DTSTAR*TIMESTEP / ICE_FRACT(I)
            TSTAR(I) = TSTARMAX
            IF (SIMLT) SICE_MLT_HTF(I) = - LSMELT * DTSTAR
            TL_1(I) = TL_1(I) + DTRDZ_1(I) * DFTL / CP
            QW_1(I) = QW_1(I) + DTRDZ_1(I) * DFQW
            FTL_1(I) = FTL_1(I) + DFTL
            FQW_1(I) = FQW_1(I) + DFQW
            EI(I) = EI(I) + DFQW
          ENDIF
          IF ( TI(I) .GT. TM ) THEN
            IF (SIMLT) SICE_MLT_HTF(I) = SICE_MLT_HTF(I) +
     &                           ICE_FRACT(I)*(TI(I) - TM)/(AI*TIMESTEP)
            TI(I) = TM
          ENDIF
        ENDIF
      ENDDO

      IF (LTIMER) THEN
      CALL TIMER('SFMELT  ',4)
      ENDIF

      RETURN
      END
