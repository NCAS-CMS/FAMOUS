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
C*LL  SUBROUTINE SF_EVAP------------------------------------------------
CLL
CLL  Purpose: Calculate surface evaporation and sublimation amounts
CLL           (without applying them to the surface stores).
CLL
CLL
CLL  Suitable for single column usage.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   4.1             New deck.
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL
CLL   4.4   Jul. 97   MOSES II         Richard Essery
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        version 2, dated 18/1/90.
CLL
CLL  Logical component covered: P245.
CLL
CLL  System task:
CLL
CLL  Documentation: UMDP 24
CLL
CLL---------------------------------------------------------------------
C*
C*L Arguments :---------------------------------------------------------
      SUBROUTINE SF_EVAP (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,
     & LAND_INDEX,TILE_INDEX,TILE_PTS,NSHYD,LTIMER,
     & ASHTF,ASHTF_SNOW,CANOPY,DTRDZ_1,FRACA,LYING_SNOW,RESFS,
     & RESFT,RHOKH_1,TILE_FRAC,SMC,WT_EXT,TIMESTEP,
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,QW_1,TL_1,TSTAR_TILE,
     & ECAN,ECAN_TILE,ESOIL,ESOIL_TILE,EXT
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS              ! IN Number of P-grid points to be
!                            !    processed.
     &,P_FIELD               ! IN Total number of P-grid points.
     &,P1                    ! IN First P-point to be processed.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.
     &,LAND_FIELD            ! IN Total number of land points.
     &,NTYPE                 ! IN Number of tiles per land point.
     &,LAND_INDEX(P_FIELD)   ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD,NTYPE)
!                            ! IN Index of tile points.
     &,TILE_PTS(NTYPE)       ! IN Number of tile points.
     &,NSHYD                 ! IN Number of soil moisture levels.

      LOGICAL
     & LTIMER                ! IN Logical for TIMER.

      REAL
     & ASHTF(P_FIELD)        ! IN Coefficient to calculate surface
!                            !    heat flux into soil or sea-ice.
     &,ASHTF_SNOW(P_FIELD)   ! IN ASHTF for snow
     &,CANOPY(LAND_FIELD,NTYPE-1)
!                            ! IN Surface/canopy water on snow-free
!                            !    land tiles (kg/m2).
     &,DTRDZ_1(P_FIELD)      ! IN -g.dt/dp for surface layer
     &,FRACA(LAND_FIELD,NTYPE-1)
!                            ! IN Fraction of surface moisture flux
!                            !    with only aerodynamic resistance
!                            !    for snow-free land tiles.
     &,LYING_SNOW(P_FIELD)   ! IN Lying snow amount (kg per sq metre).
     &,RESFS(LAND_FIELD,NTYPE-1)
!                            ! IN Combined soil, stomatal and aerodynam.
!                            !    resistance factor for fraction 1-FRACA
!                            !    of snow-free land tiles.
     &,RESFT(LAND_FIELD,NTYPE)
!                            ! IN Total resistance factor
!                            !    FRACA+(1-FRACA)*RESFS.
     &,RHOKH_1(LAND_FIELD,NTYPE)
!                            ! IN Surface exchange coefficients.
     &,TILE_FRAC(LAND_FIELD,NTYPE)
!                            ! IN Tile fractions.
     &,SMC(LAND_FIELD)       ! IN Available soil moisture (kg/m2).
     &,WT_EXT(LAND_FIELD,NSHYD)
!                            ! IN Fraction of transpiration
!                            !    extracted from each soil layer.
     &,TIMESTEP              ! IN Timestep in seconds.

      REAL
     & FQW_1(P_FIELD)        ! INOUT Surface moisture flux (kg/m2/s).
     &,FQW_TILE(LAND_FIELD,NTYPE)
!                            ! INOUT Local FQW_1 for tiles.
     &,FTL_1(P_FIELD)        ! INOUT Surface sensible heat flux (W/m2).
     &,FTL_TILE(LAND_FIELD,NTYPE)
!                            ! INOUT Local FTL_1 for tiles.
     &,QW_1(P_FIELD)         ! INOUT Total water content of lowest
!                            !       atmospheric layer (kg per kg air).
     &,TL_1(P_FIELD)         ! INOUT Liquid/frozen water temperature for
!                            !       lowest atmospheric layer (K).
     &,TSTAR_TILE(LAND_FIELD,NTYPE)
!                            ! INOUT Tile surface temperatures (K).

      REAL
     & ECAN(P_FIELD)         ! OUT Gridbox mean evaporation from canopy/
!                            !     surface store (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
     &,ECAN_TILE(LAND_FIELD,NTYPE-1)
!                            ! OUT ECAN for snow-free land tiles.
     &,ESOIL(P_FIELD)        ! OUT Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
     &,ESOIL_TILE(LAND_FIELD,NTYPE-1)
!                            ! OUT ESOIL for snow-free land tiles.
     &,EXT(LAND_FIELD,NSHYD) ! OUT Extraction of water from each
!                            !     soil layer (kg/m2/s).

!  Local and other symbolic constants :-
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
     & DFQW(LAND_FIELD)      ! Increment in GBM moisture flux.
     &,DFTL(LAND_FIELD)      ! Increment in GBM sensible heat flux.
     &,FQW_TILE_OLD(LAND_FIELD,NTYPE)
!                            ! FQW_TILE before adjustment.
     &,SNOW_FRAC(LAND_FIELD) ! Fractional snow coverage.

      REAL
     & DIFF_LAT_HTF          ! Increment in local latent heat flux.
     &,DIFF_SENS_HTF         ! Increment in local sensible heat flux.
     &,DTSTAR                ! Increment in local surface temperature.
     &,EDT                   ! Moisture flux x timestep
     &,FQW_ADJ               ! Adjustment of moisture fluxes.

      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Loop counter (land, snow or land-ice field index).
     &,K           ! Loop counter (soil level index).
     &,L           ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).

      IF (LTIMER) THEN
        CALL TIMER('SFEVAP  ',3)
      ENDIF

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          FQW_TILE_OLD(L,N) = FQW_TILE(L,N)
        ENDDO
      ENDDO

      DO N=1,NTYPE-1
        DO L=1,LAND_FIELD
          ECAN_TILE(L,N) = 0.
          ESOIL_TILE(L,N) = 0.
        ENDDO
      ENDDO

      DO L=LAND1,LAND1+LAND_PTS-1
        SNOW_FRAC(L) = TILE_FRAC(L,NTYPE)
      ENDDO

!-----------------------------------------------------------------------
! Sublimation from snow (tile NTYPE)
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS(NTYPE)
        L = TILE_INDEX(J,NTYPE)
        I = LAND_INDEX(L)
        EDT = SNOW_FRAC(L)*FQW_TILE(L,NTYPE)*TIMESTEP
        IF ( EDT .GT. LYING_SNOW(I) ) THEN
          FQW_ADJ = ( 1. - LYING_SNOW(I)/(FQW_TILE(L,NTYPE)*TIMESTEP) )
     &                                             / (1. - SNOW_FRAC(L))
          FQW_TILE(L,NTYPE) = LYING_SNOW(I) / (SNOW_FRAC(L)*TIMESTEP)
          DO N=1,NTYPE-1
            FQW_TILE(L,N) = FQW_ADJ*FQW_TILE(L,N)
          ENDDO
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Surface evaporation from and condensation onto snow-free land
! (tiles 1 to NTYPE-1)
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        ECAN(I) = 0.
        ESOIL(I) = 0.
      ENDDO

      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          IF ( FQW_TILE(L,N) .GT. 0.0 ) THEN
            ECAN_TILE(L,N) = FRACA(L,N) * FQW_TILE(L,N) / RESFT(L,N)
            ESOIL_TILE(L,N) = (1. - FRACA(L,N))*RESFS(L,N)*FQW_TILE(L,N)
     &                                                      / RESFT(L,N)
            EDT = ECAN_TILE(L,N)*TIMESTEP
            IF ( EDT .GT. CANOPY(L,N) ) THEN
              ESOIL_TILE(L,N) = ( 1. - FRACA(L,N)*CANOPY(L,N)/EDT ) *
     &                               RESFS(L,N)*FQW_TILE(L,N)/RESFT(L,N)
              ECAN_TILE(L,N) = CANOPY(L,N) / TIMESTEP
            ENDIF
          ELSE
            ECAN_TILE(L,N) = FQW_TILE(L,N)
            ESOIL_TILE(L,N) = 0.
          ENDIF
          ECAN(I) = ECAN(I) + TILE_FRAC(L,N)*ECAN_TILE(L,N)
          ESOIL(I) = ESOIL(I) + TILE_FRAC(L,N)*ESOIL_TILE(L,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Soil evapotranspiration
!-----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        EDT = ESOIL(I)*TIMESTEP
        IF ( EDT .GT. SMC(L) ) THEN
          DO N=1,NTYPE-1
            ESOIL_TILE(L,N) = SMC(L)*ESOIL_TILE(L,N) / EDT
          ENDDO
          ESOIL(I) = SMC(L) / TIMESTEP
        ENDIF
      ENDDO

      DO K=1,NSHYD
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          EXT(L,K) = WT_EXT(L,K)*ESOIL(I)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate increments to surface heat fluxes, moisture fluxes and
! temperatures
!-----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        DFTL(L) = 0.
        DFQW(L) = 0.
      ENDDO

! Snow-free land tiles
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          DIFF_LAT_HTF = LC * ( FQW_TILE(L,N) - FQW_TILE_OLD(L,N) )
          DIFF_SENS_HTF = - DIFF_LAT_HTF /
     &                               ( 1. + ASHTF(I)/(CP*RHOKH_1(L,N)) )
          FTL_TILE(L,N) = FTL_TILE(L,N) + DIFF_SENS_HTF
          DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF(I)
          TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR
          DFTL(L) = DFTL(L) + TILE_FRAC(L,N)*DIFF_SENS_HTF
          DFQW(L) = DFQW(L) + TILE_FRAC(L,N)*DIFF_LAT_HTF / LC
        ENDDO
      ENDDO

! Snow tile
      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        DIFF_LAT_HTF = (LC + LF) * ( FQW_TILE(L,N) - FQW_TILE_OLD(L,N) )
        DIFF_SENS_HTF = - DIFF_LAT_HTF /
     &                          ( 1. + ASHTF_SNOW(I)/(CP*RHOKH_1(L,N)) )
        FTL_TILE(L,N) = FTL_TILE(L,N) + DIFF_SENS_HTF
        DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF_SNOW(I)
        TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR
        DFTL(L) = DFTL(L) + SNOW_FRAC(L)*DIFF_SENS_HTF
        DFQW(L) = DFQW(L) + SNOW_FRAC(L)*DIFF_LAT_HTF / (LC + LF)
      ENDDO

!-----------------------------------------------------------------------
! Update level 1 temperature and humidity and GBM heat and moisture
! fluxes due to limited moisture availability
!-----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        TL_1(I) = TL_1(I) + DTRDZ_1(I) * DFTL(L) / CP
        QW_1(I) = QW_1(I) + DTRDZ_1(I) * DFQW(L)
        FTL_1(I) = FTL_1(I) + DFTL(L)
        FQW_1(I) = FQW_1(I) + DFQW(L)
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SFEVAP  ',4)
      ENDIF

      RETURN
      END
