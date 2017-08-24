C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C           in surface fluxes to satisfy energy balance.
C           Sub-surface snowmelt is calculated and snowdepth incremented
C           by melt and sublimation in P251.
C           R.Essery 19/1/95
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
C-----------------------------------------------------------------------
      SUBROUTINE SF_MELT(
     + P_FIELD,P1,LAND_FIELD,LAND1
     +,POINTS,LAND_MASK,LAND_PTS,LAND_INDEX
     +,ALPHA1,ASHTF,ASURF,ICE_FRACT
     +,RHOKH1_PRIME,TIMESTEP,SIMLT,SMLT,DFQW,DIFF_SENS_HTF
     +,EI,LYING_SNOW,SURF_HT_FLUX,TSTAR,TI
     +,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + P_FIELD              ! IN No. of gridpoints in the whole grid.
     +,P1                   ! IN 1st P-pt in full field to be processed.
     +,LAND_FIELD           ! IN No. of land points in the whole grid.
     +,LAND1                ! IN 1st L-pt in full field to be processed.
     +,POINTS               ! IN No. of gridpoints to be processed.
     +,LAND_PTS             ! IN No. of land points to be processed.
      LOGICAL
     + LAND_MASK(P_FIELD)   ! IN T for land points, F otherwise.
      INTEGER
     + LAND_INDEX(P_FIELD)  ! IN Index of land points on the P-grid.
C                           !    The ith element contains the position
C                           !    in whole grid of the ith land point.
       REAL
     + ALPHA1(P_FIELD)      ! IN Gradient of saturated specific
C                           !    humidity with respect to temp.
C                           !    between the bottom model layer
C                           !    and the surface.
     +,ASHTF(P_FIELD)       ! IN Forward time weighted coeff.
C                           !    to calculate the soil heat flux
C                           !    between the surface and top soil
C                           !    layer (W/m2/K).
     +,ASURF(P_FIELD)       ! IN Reciprocal areal heat capacity of
C                           !    top soil layer or sea-ice surface
C                           !    layer (m2 K / J).
     +,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is covered
C                           !    by sea-ice.
     +,RHOKH1_PRIME(P_FIELD)! IN Modified forward time-weighted
C                           !    transfer coefficient.
     +,TIMESTEP             ! IN Timestep (sec).
      LOGICAL
     + SIMLT                ! IN STASH flag for sea-ice melting ht flux.
     +,SMLT                 ! IN STASH flag for snow melting ht flux.
      REAL
     + DFQW(P_FIELD)        ! INOUT Increment to the flux of total
C                           !       water.
     +,DIFF_SENS_HTF(P_FIELD)! INOUT Increment to the sensible heat
C                           !        flux (W/m2).
     +,EI(P_FIELD)          ! INOUT Sublimation from lying snow or
C                           !       sea-ice (Kg/m2/s).
     +,LYING_SNOW(P_FIELD)  ! INOUT Lying snow (kg/m2).
     +,SURF_HT_FLUX(P_FIELD)! INOUT Net downward heat flux at surface
C                           !       over land or sea-ice fraction of
C                           !       gridbox (W/m2).
     +,TSTAR(P_FIELD)       ! INOUT Surface temperature (K).
     +,TI(P_FIELD)          ! INOUT Sea-ice surface layer temp. (K).
     +,SICE_MLT_HTF(P_FIELD)! OUT Heat flux due to melting of sea-ice
C                           !     (W/m2).
     +,SNOMLT_SURF_HTF(P_FIELD)! OUT Heat flux due to surface melting
C                              !     of snow (W/m2).
     +,SNOWMELT(P_FIELD)       ! OUT Surface snowmelt (kg/m2/s).
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

      REAL
     + DMELT                ! Temporary in calculations of melting
C                           ! heat fluxes.
     +,DIFF_EI              ! Increment to sublimation.
     +,DTSTAR               ! Increment to surface temperature.
     +,DIFF_SURF_HTF        ! Increment to surface heat flux.
     +,SNOW_MAX             ! Snow available for melting at land
C                           ! points.
     +,TSTARMAX             ! Maximum gridbox mean surface temperature
C                           ! at sea points with ice.
      INTEGER
     + I                    ! Loop counter - full horizontal field.
     +,L                    ! Loop counter - land field.
C
      IF (LTIMER) THEN
      CALL TIMER('SFMELT  ',3)
      ENDIF
      DO 1 I=P1,P1+POINTS-1
        IF (SIMLT) SICE_MLT_HTF(I) = 0.0
        IF (SMLT) SNOMLT_SURF_HTF(I) = 0.0
        SNOWMELT(I) = 0.0
    1 CONTINUE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 10 L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
C
C-----------------------------------------------------------------------
C  Melt snow if TSTAR is greater than TM.
C-----------------------------------------------------------------------
        SNOW_MAX = MAX(0.0, LYING_SNOW(I) - EI(I)*TIMESTEP )
        IF ( SNOW_MAX.GT.0.0 .AND. TSTAR(I).GT.TM ) THEN
          DMELT = ( CP + LC * ALPHA1(I) ) * RHOKH1_PRIME(I) + ASHTF(I)
          DTSTAR = - MIN ( TSTAR(I) - TM ,
     &                       LF * SNOW_MAX / ( TIMESTEP * DMELT ) )
          DMELT = DMELT + LF * ALPHA1(I) * RHOKH1_PRIME(I)
          SNOWMELT(I) = - DMELT * DTSTAR / LF
          DIFF_SENS_HTF(I) = DIFF_SENS_HTF(I) +
     &                        CP * RHOKH1_PRIME(I) * DTSTAR
          DIFF_EI = ALPHA1(I) * RHOKH1_PRIME(I) * DTSTAR
          EI(I) = EI(I) + DIFF_EI
          DFQW(I) = DFQW(I) + DIFF_EI
          DIFF_SURF_HTF = ASHTF(I) * DTSTAR
          SURF_HT_FLUX(I) = SURF_HT_FLUX(I) + DIFF_SURF_HTF
          TSTAR(I) = TSTAR(I) + DTSTAR
          IF (SMLT) SNOMLT_SURF_HTF(I) = LF*SNOWMELT(I)
        ENDIF
10    CONTINUE ! End of loop over land points
      DO 20 I=P1,P1+POINTS-1
        IF ( .NOT. LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I) .GT. 0.0 ) THEN
C-----------------------------------------------------------------------
C   Melt sea-ice if TSTAR > TSTARMAX or TI > TM.
C-----------------------------------------------------------------------
            TSTARMAX = ICE_FRACT(I)*TM + (1.0 - ICE_FRACT(I))*TFS
            IF ( TSTAR(I) .GT. TSTARMAX ) THEN
              DTSTAR = TSTARMAX - TSTAR(I)
              DMELT = (CP + (LC + LF)*ALPHA1(I))*RHOKH1_PRIME(I)
     &                    + ASHTF(I)
              DIFF_SENS_HTF(I) = CP * RHOKH1_PRIME(I) * DTSTAR
              DIFF_EI = ALPHA1(I) * RHOKH1_PRIME(I) * DTSTAR
              EI(I) = EI(I) + DIFF_EI
              DFQW(I) = DFQW(I) + DIFF_EI
              DIFF_SURF_HTF = ASHTF(I) * DTSTAR
              TI(I) =TI(I) + ASURF(I) * TIMESTEP * DIFF_SURF_HTF
              TSTAR(I) = TSTARMAX
              IF (SIMLT) SICE_MLT_HTF(I) = - DMELT * DTSTAR
            ENDIF
            IF ( TI(I) .GT. TM ) THEN
              IF (SIMLT) SICE_MLT_HTF(I) = SICE_MLT_HTF(I) +
     &                                  (TI(I) - TM)/(ASURF(I)*TIMESTEP)
              TI(I) = TM
            ENDIF

          ENDIF              ! Sea-ice
        ENDIF                ! Sea
20    CONTINUE               ! End of loop over sea points
      IF (LTIMER) THEN
      CALL TIMER('SFMELT  ',4)
      ENDIF
      RETURN
      END
