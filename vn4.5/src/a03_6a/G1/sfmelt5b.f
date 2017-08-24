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

! SUBROUTINE SF_MELT----------------------------------------------------
! Purpose : Calculates surface melting (snow and sea-ice) and increments
!           in surface fluxes to satisfy energy balance.
!           Sub-surface snowmelt is calculated and snowdepth incremented
!           by melt and sublimation in P251.
!           R.Essery 19/1/95
!
!
!
!  Model            Modification history:
! version  Date
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!-----------------------------------------------------------------------
      SUBROUTINE SF_MELT(
     & P_FIELD,P1,N_TYPES,LAND_FIELD,LAND1
     &,POINTS,LAND_MASK,LAND_PTS,LAND_INDEX
     &,ALPHA1,ASHTF,ASURF,TILE_FRAC,ICE_FRACT
     &,RHOKH1_PRIME,TIMESTEP,SIMLT,SMLT,DFQW,DIFF_SENS_HTF
     &,EI,LYING_SNOW,SURF_HT_FLUX,TSTAR_TILE,TI
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,LTIMER
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & P_FIELD              ! IN No. of gridpoints in the whole grid.
     &,P1                   ! IN 1st P-pt in full field to be processed.
     &,N_TYPES              ! IN max number of tiles per grid-box
     &,LAND_FIELD           ! IN No. of land points in the whole grid.
     &,LAND1                ! IN 1st L-pt in full field to be processed.
     &,POINTS               ! IN No. of gridpoints to be processed.
     &,LAND_PTS             ! IN No. of land points to be processed.

      LOGICAL
     & LAND_MASK(P_FIELD)   ! IN T for land points, F otherwise.
      INTEGER
     & LAND_INDEX(P_FIELD)  ! IN Index of land points on the P-grid.
!                                The ith element contains the position
!                                in whole grid of the ith land point.

       REAL
     & ALPHA1(P_FIELD,N_TYPES)
!                             IN Gradient of saturated specific
!                                humidity with respect to temp.
!                                between the bottom model layer
!                                and the surface.
     &,ASHTF(P_FIELD)       ! IN Forward time weighted coeff.
!                                to calculate the soil heat flux
!                                between the surface and top soil
!                                layer (W/m2/K).
     &,ASURF(P_FIELD)       ! IN Reciprocal areal heat capacity of
!                                top soil layer or sea-ice surface
!                                layer (m2 K / J).
     &,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is covered
!                                by sea-ice.
     &,RHOKH1_PRIME(P_FIELD,N_TYPES)
!                             IN Modified forward time-weighted
!                                transfer coefficient.
     &,TILE_FRAC(P_FIELD,N_TYPES)
!                             IN Fraction of gridbox which is covered
!                                by a tile.
     &,TIMESTEP             ! IN Timestep (sec).

      LOGICAL
     & SIMLT                ! IN STASH flag for sea-ice melting ht flux.
     &,SMLT                 ! IN STASH flag for snow melting ht flux.

      REAL
     & DFQW(P_FIELD,N_TYPES)! INOUT Increment to the flux of total
!                                   water.
     &,DIFF_SENS_HTF(P_FIELD,N_TYPES)
!                             INOUT Increment to the sensible heat
!                                    flux (W/m2).
     &,EI(P_FIELD,N_TYPES)  ! INOUT Sublimation from lying snow or
!                                   sea-ice (Kg/m2/s).
     &,LYING_SNOW(P_FIELD)  ! INOUT Lying snow (kg/m2).
     &,SURF_HT_FLUX(P_FIELD,N_TYPES)
!                             INOUT Net downward heat flux at surface
!                                   over land or sea-ice fraction of
!                                   gridbox (W/m2).
     &,TSTAR_TILE(P_FIELD,N_TYPES)
!                             INOUT Surface temperature (K).
     &,TI(P_FIELD)          ! INOUT Sea-ice surface layer temp. (K).
     &,SICE_MLT_HTF(P_FIELD)! OUT Heat flux due to melting of sea-ice
!                                 (W/m2).
     &,SNOMLT_SURF_HTF(P_FIELD)
!                             OUT Heat flux due to surface melting
!                                 of snow (W/m2).
     &,SNOWMELT(P_FIELD,N_TYPES)
!                            OUT Surface snowmelt (kg/m2/s).

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
     & DMELT          ! Temporary in calculations of melting heat fluxes
     &,DIFF_EI        ! Increment to sublimation.
     &,DTSTAR         ! Increment to surface temperature.
     &,DIFF_SURF_HTF  ! Increment to surface heat flux.
     &,SNOW_MAX       ! Snow available for melting at land points.
     &,TSTARMAX       ! Maximum gridbox mean surface temperature at sea
!                       points with ice.

      INTEGER
     & I                    ! Loop counter - full horizontal field.
     &,L                    ! Loop counter - land field.
     &,ITILE                ! Loop counter - land tiles.


      IF (LTIMER) THEN
        CALL TIMER('SFMELT  ',3)
      ENDIF

      DO I=P1,P1+POINTS-1
        IF (SIMLT) SICE_MLT_HTF(I) = 0.0
        IF (SMLT) SNOMLT_SURF_HTF(I) = 0.0
      ENDDO


!-----------------------------------------------------------------------
!  Melt land snow if TSTAR_TILE is greater than TM.
!-----------------------------------------------------------------------
      DO ITILE=1,N_TYPES

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)


          SNOW_MAX = MAX(0.0, LYING_SNOW(I) - EI(I,ITILE)*TIMESTEP )

          IF ( SNOW_MAX.GT.0.0 .AND. TSTAR_TILE(I,ITILE).GT.TM ) THEN

            DMELT = ( CP + LC * ALPHA1(I,ITILE) )
     &                    * RHOKH1_PRIME(I,ITILE) + ASHTF(I)

            DTSTAR = - MIN ( TSTAR_TILE(I,ITILE) - TM ,
     &                       LF * SNOW_MAX / ( TIMESTEP * DMELT ) )

            DMELT = DMELT + LF * ALPHA1(I,ITILE) * RHOKH1_PRIME(I,ITILE)

            SNOWMELT(I,ITILE) = - DMELT * DTSTAR / LF

            DIFF_SENS_HTF(I,ITILE) = DIFF_SENS_HTF(I,ITILE) +
     &                        CP * RHOKH1_PRIME(I,ITILE) * DTSTAR

            TSTAR_TILE(I,ITILE) = TSTAR_TILE(I,ITILE) + DTSTAR

            DIFF_SURF_HTF = ASHTF(I) * DTSTAR

            SURF_HT_FLUX(I,ITILE) = SURF_HT_FLUX(I,ITILE) +
     &                                 DIFF_SURF_HTF

            DIFF_EI = ALPHA1(I,ITILE) * RHOKH1_PRIME(I,ITILE) * DTSTAR

            EI(I,ITILE) = EI(I,ITILE) + DIFF_EI

            DFQW(I,ITILE) = DFQW(I,ITILE) + DIFF_EI

            IF (SMLT)
     &        SNOMLT_SURF_HTF(I) = SNOMLT_SURF_HTF(I) +
     &                        LF*SNOWMELT(I,ITILE) * TILE_FRAC(I,ITILE)


          ENDIF
        ENDDO ! End of loop over land points
      ENDDO ! loop over land tiles



!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR_TILE > TSTARMAX or TI > TM.
!-----------------------------------------------------------------------
      DO I=P1,P1+POINTS-1
        IF ( .NOT. LAND_MASK(I) .AND. ICE_FRACT(I) .GT. 0.0 ) THEN

          TSTARMAX = ICE_FRACT(I)*TM + (1.0 - ICE_FRACT(I))*TFS

          IF ( TSTAR_TILE(I,1) .GT. TSTARMAX ) THEN
            DTSTAR = TSTARMAX - TSTAR_TILE(I,1)
            DMELT = (CP + (LC + LF)*ALPHA1(I,1))*RHOKH1_PRIME(I,1)
     &                  + ASHTF(I)
            DIFF_SENS_HTF(I,1) = CP * RHOKH1_PRIME(I,1) * DTSTAR
            DIFF_EI = ALPHA1(I,1) * RHOKH1_PRIME(I,1) * DTSTAR
            EI(I,1) = EI(I,1) + DIFF_EI
            DFQW(I,1) = DFQW(I,1) + DIFF_EI
            DIFF_SURF_HTF = ASHTF(I) * DTSTAR
            TI(I) =TI(I) + ASURF(I) * TIMESTEP * DIFF_SURF_HTF
            TSTAR_TILE(I,1) = TSTARMAX
            IF (SIMLT) SICE_MLT_HTF(I) = - DMELT * DTSTAR

          ENDIF !end of TSTAR_TILE > TSTARMAX block

          IF ( TI(I) .GT. TM ) THEN
            IF (SIMLT) SICE_MLT_HTF(I) = SICE_MLT_HTF(I) +
     &                                (TI(I) - TM)/(ASURF(I)*TIMESTEP)
            TI(I) = TM
          ENDIF ! end of TI > TM block

        ENDIF  ! Sea-ice points

      ENDDO ! End of loop over p_points


      IF (LTIMER) THEN
        CALL TIMER('SFMELT  ',4)
      ENDIF

      RETURN
      END
