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
!!!   SUBROUTINE SF_EVAP------------------------------------------------
!!!
!!!  Purpose: Calculate surface evaporation and sublimation amounts
!!!           (without applying them to the surface stores).
!!!           Also calculate heat flux due to sea-ice melting.
!!!           Also calculate 1.5 metre T and Q.
!!!
!!!
!!!  Suitable for single column usage.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        version 2, dated 18/1/90.
!!!
!!!   4.3             New deck.
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  Logical component covered: P245.
!!!
!!!  System task:
!!!
!!!  Documentation: UMDP 24
!!!
!!!---------------------------------------------------------------------

!  Arguments :---------------------------------------------------------
      SUBROUTINE SF_EVAP (
     & P_FIELD,P1,N_TYPES,LAND_FIELD,LAND1,GAMMA
     &,POINTS,BL_LEVELS,LAND_MASK,LAND_PTS,LAND_INDEX
     &,TILE_FRAC,ALPHA1,ASURF,ASHTF,CANOPY,CATCH
     &,DTRDZ,DTRDZ_RML,E_SEA,FRACA
     &,ICE_FRACT,NRML,RHOKH_1,SMC,TIMESTEP,CER1P5M,CHR1P5M
     &,PSTAR,RESFS,RESFT,Z0M,Z0H,SQ1P5,ST1P5,SIMLT,SMLT
     &,FTL,FTL_TILE,FQW,FQW_TILE,LYING_SNOW,QW,SURF_HT_FLUX
     &,TL,TSTAR_TILE,TSTAR_GB,TI,ECAN_GB,ES,EI_GB
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT_GB
     &,H_BLEND,HEAT_BLEND_FACTOR,QCL_1,QCF_1,Z1_TQ    
     &,Q1P5M,T1P5M,LTIMER
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & P_FIELD                ! IN No. of gridpoints in the whole grid.
     &,P1                     ! IN 1st P-pt in full field to be
!                                  processed.
     &,N_TYPES                ! IN No. of land tiles
     &,LAND_FIELD             ! IN No. of landpoints in the whole grid.
     &,LAND1                  ! IN 1st L-pt in full field to be
!                                  processed.
     &,POINTS                 ! IN No. of gridpoints to be processed.
     &,BL_LEVELS              ! IN No. of levels treated by b.l. scheme.
     &,LAND_PTS               ! IN No. of land points to be processed.

      LOGICAL
     & LAND_MASK(P_FIELD)     ! IN T for land points, F otherwise.

      INTEGER
     & LAND_INDEX(P_FIELD)    ! IN Index of land points on the P-grid.
!                                  The ith element contains the position
!                                  in whole grid of the ith land point.

      REAL

     & ALPHA1(P_FIELD,N_TYPES)! IN Gradient of saturated specific
!                                  humidity with respect to temp.
!                                  between the bottom model layer
!                                  and the surface.
     &,ASURF(P_FIELD)         ! IN Soil coefficient from P242 (m2 K per
!                                  per Joule * timestep).
     &,ASHTF(P_FIELD)         ! IN Coefficient to calculate the soil
!                                  heat flux between the surface and
!                                  top soil layer (W/m2/K)
     &,CANOPY(LAND_FIELD)     ! IN Gridbox mean canopy / surface water
!                                  store (kg/m2).
     &,CATCH(LAND_FIELD,N_TYPES)
!                               IN Canopy / surface water store capacity
!                                  (kg per sq m).
     &,CER1P5M(P_FIELD)       ! IN Interpolation coefficient, from P243
     &,CHR1P5M(P_FIELD)       ! IN Interpolation coefficient, from P243.
     &,DTRDZ(P_FIELD,BL_LEVELS)!IN -g.dt/dp for each model layer on
!                                  p-grid From P244 ((kg/m2/s)**-1).
     &,DTRDZ_RML(P_FIELD)     ! IN -g.dt/dp for the rapidly mixing layer
!                                  (if it exists) on the p-grid from
!                                  P244
     &,E_SEA(P_FIELD)         ! IN Evaporation from sea (weighted with
!                                  leads fraction at sea-ice points).
     &,FRACA(P_FIELD,N_TYPES) ! IN Fraction of surface moisture flux
!                                  with only aerodynamic resistance.
!                                  Diagnostics defined on land and sea.
     &,GAMMA(BL_LEVELS)       ! IN Weights for implicit BL scheme.
     &,H_BLEND(P_FIELD)       ! IN Blending height
     &,HEAT_BLEND_FACTOR(P_FIELD)
!                               IN Used for tile adjustment
     &,Z1_TQ(P_FIELD)         ! IN Height of lowest tq level (m).  
     &,ICE_FRACT(P_FIELD)     ! IN Fraction of gridbox which is covered
!                                  by sea-ice (decimal fraction, but
!                                  mostly this sub-component assumes it
!                                  to be either 1.0 or 0.0 precisely).
!                                  NB Dimension is PFIELD not LAND_FIELD
!                                  for snow on sea-ice in coupled model
!                                  runs.
     &,PSTAR(P_FIELD)         ! IN Surface pressure (Pa).
     &,QCL_1(P_FIELD)         ! IN Liquid water at level 1
     &,QCF_1(P_FIELD)         ! IN frozen water at level 1
     &,RESFS(P_FIELD,N_TYPES) ! IN Combined soil, stomatal and
!                                  aerodynamic resistance factor
     &,RESFT(P_FIELD,N_TYPES) ! IN Total resistance factor
!                                  FRACA+(1-FRACA)*RESFS.
     &,RHOKH_1(P_FIELD,N_TYPES)!IN Turbulent surface exchange
!                                   coefficient for sensible heat.
     &,SMC(LAND_FIELD,N_TYPES)! IN Soil moisture content (kg per sq m).
     &,TILE_FRAC(P_FIELD,N_TYPES)
!                               IN fractional coverage for each tile
     &,TIMESTEP               ! IN Timestep (sec).
     &,Z0M(P_FIELD,N_TYPES)   ! IN Roughness length for momentum (m)
     &,Z0H(P_FIELD,N_TYPES)   ! IN Roughness length for heat and
!                                  moisture

      INTEGER
     & NRML(P_FIELD)          ! IN  The Number of model layers in the
!                                   Rapidly Mixing Layer.

      LOGICAL
     & SQ1P5                  ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                  ! IN STASH flag for 1.5-metre temperature.
     &,SIMLT                  ! IN STASH flag for sea-ice melting ht
!                                  flux.
     &,SMLT                   ! IN STASH flag for snow melting ht flux.

      REAL
     & FTL(P_FIELD,BL_LEVELS) ! INOUT Sensible heat flux from layer k-1
!                                     to layer k (W/sq m).  From P243
!                                     and P244, units changed in P24
!                                     top level.
     &,FTL_TILE(P_FIELD,N_TYPES)
!                               INOUT Sensible surf heat flux for tile
     &,FQW(P_FIELD,BL_LEVELS) ! INOUT Turbulent moisture flux from level
!                                     k-1 to k (kg/sq m/s). From P243/4.
!                                     Diagnostics defined on land and
!                                     sea
     &,FQW_TILE(P_FIELD,N_TYPES)
!                               INOUT Moisture flux for tile
     &,LYING_SNOW(P_FIELD)    ! INOUT Lying snow (kg per sq m).
     &,QW(P_FIELD,BL_LEVELS)  ! INOUT Total water content (kg(water)/
!                                     kg(air)).  From P243/4.
     &,SURF_HT_FLUX(P_FIELD,N_TYPES)
!                               INOUT Net downward heat flux at surface
!                                     over land or sea-ice fraction of
!                                     gridbox (W/m2).
     &,TL(P_FIELD,BL_LEVELS)  ! INOUT Liquid/frozen water temperature K.
     &,TSTAR_GB(P_FIELD)      ! INOUT mean land Surface temperature (K)
     &,TSTAR_TILE(P_FIELD,N_TYPES)
!                               INOUT Tile surface temperature (K).
     &,TI(P_FIELD)            ! INOUT Sea-ice surface layer temp. (K).

! OUTPUT

      REAL
     & ECAN_GB(P_FIELD)       ! OUT Gridbox mean evap. from canopy/
!                                   surface store (kg/m2/s).
!                                   Zero over sea and sea-ice.
     &,ES_GB(P_FIELD)         ! OUT Surface evapotranspiration (through
!                                   a resistance which is not entirely
!                                   aerodynamic).  Always non-negative.
!                                   Kg per sq m per sec.
!                                   Diagnostics defined on land and sea.
     &,EI_GB(P_FIELD)         ! OUT Sublimation from lying snow or sea-
!                                   ice (kg per sq m per s).
      REAL
     & SICE_MLT_HTF(P_FIELD)  ! OUT Heat flux due to melting of sea-ice
!                                   (Watts per square metre).
     &,SNOMLT_SURF_HTF(P_FIELD)!OUT Heat flux due to surface melting
!                                   of snow (W/m2).
     &,SNOWMELT_GB(P_FIELD)   ! OUT Surface snowmelt (kg/m2/s).
     &,Q1P5M(P_FIELD)         ! OUT Specific humidity at screen height
!                                    of 1.5 metres (kg water / kg air).
     &,T1P5M(P_FIELD)         ! OUT Temperature at 1.5 metres above the
!                                   surface (K).


!  External subprogram(s) required :-
      EXTERNAL QSAT,SF_MELT
      EXTERNAL TIMER


C*L  Local and other symbolic constants used :-
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
C*L-----------COMDECK C_HT_M FOR SUBROUTINE SF_EXCH----------
C Z10M  = height of 10m level for diagnostic calculations (m).
C Z1P5M = height of 1.5m level for diagnostic calculations (m).
      REAL Z10M,Z1P5M

      PARAMETER(Z10M  = 10.0,
     &          Z1P5M = 1.5)
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

      REAL KAPPAI
      PARAMETER (
     + KAPPAI=2.09          ! Thermal conductivity of sea-ice (W per
C                           ! m per K).
     +)
      REAL DE
      PARAMETER (
     + DE = 0.1             ! Effective thickness of sea-ice surface
C                           ! layer (m).
     +)
C-----------------------------------------------------------------------

      REAL GRCP,LS,LCRCP,LSRCP
      PARAMETER (
     & GRCP=G/CP              ! Accn due to gravity / standard heat
!                               capacity of air at const pressure.
     &,LS=LF+LC               ! Latent heat of sublimation.
     &,LCRCP=LC/CP            ! Evaporation-to-dT conversion factor.
     &,LSRCP=LS/CP            ! Sublimation-to-dT conversion factor.
     &)


!! Workspace

      REAL
     & DFQW(P_FIELD,N_TYPES)  ! Adjustment increment to the flux of
!                               total water for tile
     &,DFQW_GB(P_FIELD)       ! Adjustment increment to the flux of
!                               total water
     &,DIFF_SENS_HTF(P_FIELD,N_TYPES)
!                               Adjustment increment to the sensible
!                               heat flux
     &,DQW(P_FIELD)           ! Increment to specific humidity for
!                               current tile
     &,DTL(P_FIELD)           ! Increment to temperature for current
!                               tile
     &,DQW_GB(P_FIELD)        ! Increment to specific humidity
     &,DTL_GB(P_FIELD)        ! Increment to temperature
     &,D_S_H_GB(P_FIELD)      ! Change in sens. heat flux over gridbox
     &,DTRDZ_1(P_FIELD)       ! -g.dt/dp for surface layer or rml if it
!                               exists from P244 ((kg/sq m/s)**-1).
     &,ECAN(P_FIELD,N_TYPES)  ! Tile evaporation from canopy/
!                               surface store (kg per sq m per s).
!                               Zero over sea and sea-ice.
     &,EOLD(P_FIELD,N_TYPES)  ! Used to store initial value of evap.
!                               for current tile from P244
     &,EOLD_GB(P_FIELD)       ! Used to store initial mean value of
!                               evap.for gridbox from P244
     &,EI(P_FIELD,N_TYPES)    ! Sublimation from lying snow or sea-
!                               ice (kg per sq m per s).
     &,ES(P_FIELD,N_TYPES)    ! Surface evapotranspiration (through
!                               a resistance which is not entirely
!                               aerodynamic).  Always non-negative.
!                               Kg per sq m per sec
     &,EW(P_FIELD)            ! Total surface flux of water, excluding
!                               sublimation/frost deposition, over land.
     &,LEOLD(P_FIELD)         ! Used to store initial value of latent
!                               heat flux from P244
     &,QS(P_FIELD)            ! Used for saturated specific humidity
!                               at surface, in Q1P5M calculation.
     &,QSTAR_GB(P_FIELD)      ! Qstar in Q1P5M calculation.
     &,RHOKH1_PRIME(P_FIELD,N_TYPES)
!                               Modified forward time-weighted transfer
!                               coefficient
     &,SNOWMELT(P_FIELD,N_TYPES)
!                               Surface snowmelt (kg/m2/s).

!  Local scalars
      REAL
     & DIFF_LAT_HTF        ! Increment to the latent heat flux.
     &,DIFF_SURF_HTF       ! Increment to the surface heat flux.
     &,DTSTAR              ! Increment for surface temperature.
     &,EA                  ! Surface evaporation with only aero-
!                            dynamic resistance (+ve), or condens-
!                            ation (-ve), averaged over gridbox
!                            (kg/m2/s).
     &,EADT                ! EA (q.v.) integrated over timestep.
     &,ECANDT              ! ECAN (q.v.) integrated over timestep.
     &,EDT                 ! E=FQW(,1) (q.v.) integrated over timestep.
     &,EIDT                ! EI (q.v.) integrated over timestep.
     &,ESDT                ! ES (q.v.) integrated over timestep.
     &,ESL                 ! ES (q.v.) without fractional weighting
!                            factor FRACS ('L' is for 'local')
!                            (kg/m2/s).
     &,ESLDT               ! ESL (q.v.) integrated over timestep.
     &,FRACS               ! Fraction of gridbox at which moisture flux
!                            is additionally impeded by a surface and/or
!                            stomatal resistance.
     &,QW_BLEND            ! QW at blending height
     &,TL_BLEND            ! TL at blending height

      INTEGER
     & I                   ! Loop counter - full horizontal field index.
     &,ITILE               ! Loop counter - land tile index.
     &,L                   ! Loop counter - land field index.
     &,K                   ! Loop counter in the vertical.
     &,KM1                 ! K - 1

      IF (LTIMER) THEN
        CALL TIMER('SFEVAP  ',3)
      ENDIF

!-----------------------------------------------------------------------
!! 1. Initialise some output variables and flux increments to zero.
!-----------------------------------------------------------------------


      DO I=P1,P1+POINTS-1
        ECAN_GB(I) = 0.0
        ES_GB(I) = 0.0
        EI_GB(I) = 0.0
        D_S_H_GB(I) = 0.0
        DFQW_GB(I) = 0.0
        SNOWMELT_GB(I) = 0.0
      ENDDO

      DO ITILE=1,N_TYPES
        DO I=P1,P1+POINTS-1
          DIFF_SENS_HTF(I,ITILE) = 0.0
          DFQW(I,ITILE) = 0.0
          EI(I,ITILE) = 0.0
          SNOWMELT(I,ITILE) = 0.0
          ECAN(I,ITILE) = 0.0
        ENDDO
      ENDDO

!---------------------------------------------------------------------
!! 2. Do calculations for land points.
!---------------------------------------------------------------------

CMIC$ DO ALL VECTOR SHARED(P_FIELD, LAND_FIELD, BL_LEVELS, LAND1,
CMIC$1   LAND_PTS, LAND_INDEX, ESL, TIMESTEP, ES, LYING_SNOW, ECAN,
CMIC$2   EA, CATCH, CANOPY, SMC, EI, TSTAR_TILE, FQW_TILE, EOLD,
CMIC$3   LEOLD, P1,POINTS,LC,LF,TM,LAND_MASK,EW,FRACA,RESFT,RESFS)
CMIC$4   PRIVATE(I, L, ESLDT,
CMIC$5   ESDT, EADT, EDT, ECANDT, FRACS, EIDT)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC

      DO ITILE=1,N_TYPES

        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)

          IF (FQW_TILE(I,ITILE).EQ.0.0) THEN
            EA = 0.0
            ESL = 0.0
          ELSE
            EA = FQW_TILE(I,ITILE) / RESFT(I,ITILE) * FRACA(I,ITILE)
            ESL = FQW_TILE(I,ITILE) / RESFT(I,ITILE) * RESFS(I,ITILE)
          END IF
          ES(I,ITILE) = ESL * (1. - FRACA(I,ITILE))

!-----------------------------------------------------------------------
!! 2.1 Calculate fluxes integrated over timestep.
!-----------------------------------------------------------------------

          ESLDT = ESL * TIMESTEP
          EADT = EA * TIMESTEP
          ESDT = ES(I,ITILE) * TIMESTEP
          EDT = EADT + ESDT

!-----------------------------------------------------------------------
!! 2.2 Do calculations for snow-free land.  Canopy processes operate.
!!     LYING_SNOW is defined on sea and land points for snow on sea-ice
!!     in coupled model runs.
!-----------------------------------------------------------------------

          IF (LYING_SNOW(I).LE.0.0) THEN

!**********************************************************************
! Store initial value of evaporation and latent heat flux
!**********************************************************************

            EOLD(I,ITILE) = FQW_TILE(I,ITILE)
            EOLD_GB(I) = FQW(I,1)
            LEOLD(I) = FQW_TILE(I,ITILE) * LC
            IF (EDT.GE.0.0) THEN

!-----------------------------------------------------------------------
!! 2.2.1 Non-negative moisture flux over snow-free land.
!-----------------------------------------------------------------------

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!   (a) Water in canopy and soil is assumed to be liquid, so all
!!       positive moisture flux over snow-free land is evaporation
!!       rather than sublimation, even if TSTAR_TILE is less than or
!!       equal to TM.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

              ECAN(I,ITILE) = EA
              ECANDT = EADT

!  If EDT is non-negative, then ECANDT must be non-negative.

              FRACA(I,ITILE) = 0.0
              IF (CATCH(L,ITILE).GT.0.0)
     &          FRACA(I,ITILE) = CANOPY(L) / CATCH(L,ITILE)
              IF (CANOPY(L).LT.ECANDT) THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!   (b) It is assumed that any 'canopy' moisture flux in excess of the
!!       current canopy water amount is in fact soil evaporation.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!        This situation is highly improbable - it will occur at, at
!        most, a few gridpoints in any given timestep.

                FRACS = 1.0 - FRACA(I,ITILE)*( CANOPY(L) / ECANDT )
                ESDT = ESLDT * FRACS
                ECANDT = CANOPY(L)
                ECAN(I,ITILE) = ECANDT / TIMESTEP
                ES(I,ITILE) = ESDT / TIMESTEP
              ENDIF

!  (The canopy store is depleted by evaporation in P252, and not here,
!   according to the formula: CANOPY=CANOPY-ECANDT)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!   (c) Adjustments to evaporation from soil as calculated so far :-
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

              IF (SMC(L,ITILE).LE.0.0) THEN

!!   (i) If there is currently no soil moisture, there must be no
!!       evaporation of soil moisture, so this flux is set to zero.

                ESDT = 0.0
                ES(I,ITILE) = 0.0
              ELSEIF (SMC(L,ITILE).LT.ESDT) THEN


!!  (ii) Ensure that the soil evaporation is not greater than the
!!       current soil moisture store.
!        This situation is extremely unlikely at any given gridpoint
!        at any given timestep.

                ESDT = SMC(L,ITILE)
                ES(I,ITILE) = ESDT / TIMESTEP
              ENDIF

!  (The soil moisture store is depleted by evaporation in P253, and not
!   here, using the formula:  SMC=SMC-ESDT)

              EW(I) = ECAN(I,ITILE) + ES(I,ITILE)
              EI(I,ITILE) = 0.0

!-----------------------------------------------------------------------
!! 2.2.2 Negative moisture flux onto snow-free land above freezing
!-----------------------------------------------------------------------
!!       (i.e. condensation onto snow-free land).  The whole flux is
!!       into the surface/canopy store.

            ELSEIF (TSTAR_TILE(I,ITILE).GT.TM) THEN ! ELSE of
!                                     !  evaporation/condensation block.

!  Condensation implies ES=0, so ECAN=EA=EW=E (=FQW(,1))

              ECAN(I,ITILE) = FQW_TILE(I,ITILE)
              ES(I,ITILE) = 0.0
              EW(I) = ECAN(I,ITILE)
              EI(I,ITILE) = 0.0

!  (The canopy store is augmented by interception of condensation at
!   P252, and not here.)

!-----------------------------------------------------------------------
!! 2.2.3 Negative moisture flux onto snow-free land below freezing
!!       (i.e. deposition of frost).
!-----------------------------------------------------------------------

            ELSE      ! ELSE of condensation / frost deposition block.
              EI(I,ITILE) = FQW_TILE(I,ITILE)
              ES(I,ITILE) = 0.0
              EW(I) = 0.0

!  (Negative EI is used to increment the snowdepth store - there is
!   no separate "frost" store.  This incrementing is done in P251,
!   according to:  LYING_SNOW = LYING_SNOW - EI*TIMESTEP)

            ENDIF  ! End of evaporation/condensation/deposition block.

!-----------------------------------------------------------------------
!! 2.3 Do calculations for snow-covered land.
!-----------------------------------------------------------------------

          ELSEIF (LYING_SNOW(I).LE.EDT) THEN     ! ELSEIF of no-snow.

!**********************************************************************
! Store initial value of evaporation and latent heat flux
!**********************************************************************

            EOLD(I,ITILE) = FQW_TILE(I,ITILE)
            EOLD_GB(I) = FQW(I,1)
            LEOLD(I) = FQW(I,1) * ( LC + LF )

!-----------------------------------------------------------------------
!! 2.3.1 Shallow snow (lying snow or frost which is being exhausted
!!       by evaporation).  All the snow is sublimated, the remaining
!!       moisture flux being taken from the canopy and soil, with all
!!       the palaver of section 1.2.1 above.
!-----------------------------------------------------------------------

!        This is extremely unlikely at more than one or two gridpoints
!        at any given timestep, yet the complicated logic probably
!        slows down the routine considerably - this section is a
!        suitable candidate for further consideration as regards
!        making the model optimally efficient.

            EI(I,ITILE) = LYING_SNOW(I) / TIMESTEP
            EIDT = LYING_SNOW(I)

!  Set EDT = ( E - SNOSUB ) * TIMESTEP.  This is the moisture in kg per
!  square metre left over to be evaporated from the canopy and soil.
!  N.B.  E=FQW(,1)

            EDT = EDT - EIDT

!  (Snowdepth is decreased using EI at P251, and not here.  The formula
!   used is simply:  LYING_SNOW = LYING_SNOW - EI*TIMESTEP.)

!  Now that all the snow has sublimed, canopy processes come into
!  operation (FRACA no longer necessarily equal to 1).

            FRACA(I,ITILE) = 0.0
            IF (CATCH(L,ITILE).GT.0.0)
     &        FRACA(I,ITILE) = CANOPY(L) / CATCH(L,ITILE)
            ECANDT = EDT * FRACA(I,ITILE)
            IF (CANOPY(L).LT.ECANDT) THEN

!  Dry out the canopy completely and assume the remaining moisture flux
!  is soil evaporation.

              FRACS = 1.0 - FRACA(I,ITILE)*( CANOPY(L) / ECANDT )
              ESDT = EDT * FRACS
              ECANDT = CANOPY(L)
            ELSE

!  Calculate soil evaporation.

              FRACS = 1.0 - FRACA(I,ITILE)
              ESDT = EDT * FRACS
            ENDIF
            ECAN(I,ITILE) = ECANDT / TIMESTEP
            ES(I,ITILE) = ESDT / TIMESTEP

!  (ECAN is used to deplete the canopy store at P252, and not here.  The
!   formula used is simply:  CANOPY = CANOPY - ECAN*TIMESTEP.)

!  Evaporation from soil.

            IF (SMC(L,ITILE).LE.0.0) THEN

!  No evaporation from soil possible when there is no soil moisture.

              ESDT = 0.0
              ES(I,ITILE) = 0.0
            ELSEIF (SMC(L,ITILE).LT.ESDT) THEN

!  Limit evaporation of soil moisture in the extremely unlikely event
!  that soil moisture is exhausted by the evaporation left over from
!  sublimation which exhausted the snow store.

              ESDT = SMC(L,ITILE)
              ES(I,ITILE) = ESDT / TIMESTEP
            ENDIF

!  (ES is used to deplete the soil moisture store at P253, and not here,
!   according to the formula:  SMC = SMC - ES*TIMESTEP.)

            EW(I) = ECAN(I,ITILE) + ES(I,ITILE)

!-----------------------------------------------------------------------
!! 2.3.2 Deep snow (i.e. not being exhausted by evaporation).  This
!!       covers two cases: (a) sublimation from deep snow (if total
!!       moisture flux over the timestep is non-negative but less than
!!       the lying snow amount), and (b) deposition onto an already
!!       snowy surface (if the total moisture flux is negative and
!!       the lying snow amount is positive).
!-----------------------------------------------------------------------

          ELSE          ! ELSE of shallow snow / deep snow block.
            EI(I,ITILE) = FQW_TILE(I,ITILE)
            EW(I) = 0.0

!**********************************************************************
! Store initial value of evaporation and latent heat flux
!**********************************************************************

            EOLD(I,ITILE) = FQW_TILE(I,ITILE)
            EOLD_GB(I) = FQW(I,1)
            LEOLD(I) = FQW_TILE(I,ITILE) * ( LC + LF )

!  (EI is used to increase or decrease the snowdepth at P251, and not
!   here, according to the formula:
!   LYING_SNOW = LYING_SNOW - EI*TIMESTEP . )

          ENDIF         ! End of no snow/shallow snow/deep snow block.
          FQW_TILE(I,ITILE) = EW(I) + EI(I,ITILE)

        ENDDO ! end of loop over land points

!  Split loop 2 here so that it will vectorise.

CMIC$ DO ALL VECTOR SHARED(DTRDZ_1,DTRDZ,RHOKH_1,GAMMA,
CMIC$1 NRML,DTRDZ_RML,EI,EW,LEOLD,DIFF_LAT_HTF,FQW,
CMIC$2 EOLD,DFQW,ASHTF,DIFF_SENS_HTF,DIFF_SURF_HTF,
CMIC$3 ASURF,TIMESTEP,TSTAR_TILE,LAND_INDEX,RHOKH1_PRIME,SURF_HT_FLUX)
CMIC$4 PRIVATE(DTSTAR,I)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC

        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)

!***********************************************************************
!  2.4 Calculate increments to surface and subsurface temperatures,
!      surface heat and moisture fluxes and soil heat flux. Apply
!      increments to TSTAR_TILE to give interim values before any
!      snowmelt.
!***********************************************************************
          IF (NRML(I).GE.2) THEN
            DTRDZ_1(I) = DTRDZ_RML(I)
          ELSE
            DTRDZ_1(I) = DTRDZ(I,1)
          ENDIF
          RHOKH1_PRIME(I,ITILE) = 1.0 /
     &    ( 1.0 / RHOKH_1(I,ITILE) + GAMMA(1) * DTRDZ_1(I))

          DIFF_LAT_HTF = (LC + LF) * EI(I,ITILE) +
     &                    LC * EW(I) - LEOLD(I)
          DFQW(I,ITILE) = FQW_TILE(I,ITILE) - EOLD(I,ITILE)

          DIFF_SENS_HTF(I,ITILE) = - DIFF_LAT_HTF /
     &               ( 1. + ASHTF(I) /(RHOKH1_PRIME(I,ITILE) * CP) )

          DIFF_SURF_HTF = - DIFF_LAT_HTF / ( 1.0 +
     &                      RHOKH1_PRIME(I,ITILE) * CP / ASHTF(I) )

          SURF_HT_FLUX(I,ITILE) = SURF_HT_FLUX(I,ITILE) +
     &                                    DIFF_SURF_HTF
          DTSTAR = DIFF_SURF_HTF / ASHTF(I)
          TSTAR_TILE(I,ITILE) = TSTAR_TILE(I,ITILE) + DTSTAR

        ENDDO !End of loop over land points
      ENDDO  !End of tile loop

!-----------------------------------------------------------------------
!! 2.5 Do calculations for sea points.
!-----------------------------------------------------------------------

CMIC$ DO ALL VECTOR SHARED(P_FIELD, BL_LEVELS, P1, POINTS,NRML,
CMIC$1  LAND_MASK, ES, EI, EOLD,
CMIC$2  ICE_FRACT, FQW, E_SEA,DTRDZ_RML,
CMIC$3  TSTAR_TILE, TSTAR_GB, SMLT, SICE_MLT_HTF, KAPPAI,
CMIC$4  DTRDZ_1,DTRDZ,RHOKH_1,GAMMA,RHOKH1_PRIME,
CMIC$5  TIMESTEP,TM,TFS) PRIVATE(I, TSTARMAX)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC


      DO I=P1,P1+POINTS-1
        IF (.NOT.LAND_MASK(I)) THEN

!-----------------------------------------------------------------------
!! 2.5.1 Set soil and canopy evaporation amounts to zero, and set
!!       sublimation to zero for liquid sea points.
!-----------------------------------------------------------------------

          ES(I,1) = 0.0
          EI(I,1) = 0.0
!-----------------------------------------------------------------------
!! 2.5.3 For sea-ice points :-
!-----------------------------------------------------------------------

          IF (ICE_FRACT(I).GT.0.0) THEN
            EOLD_GB(I) = FQW(I,1)
            EI(I,1) = FQW(I,1) - E_SEA(I)
            IF (NRML(I).GE.2) THEN
              DTRDZ_1(I) = DTRDZ_RML(I)
            ELSE
              DTRDZ_1(I) = DTRDZ(I,1)
            ENDIF
            RHOKH1_PRIME(I,1) = 1.0 / ( 1.0 / RHOKH_1(I,1)
     &                          + ICE_FRACT(I)*GAMMA(1)*DTRDZ_1(I) )
          ENDIF     ! End of liquid sea/sea-ice block.

        ENDIF       ! End of sea point calculations.
      ENDDO  !End of loop over points

!-----------------------------------------------------------------------
!  Calculate fluxes and increments associated with melting of snow
!  or sea-ice.
!-----------------------------------------------------------------------


      CALL SF_MELT(P_FIELD,P1,N_TYPES,LAND_FIELD,LAND1
     &,POINTS,LAND_MASK,LAND_PTS,LAND_INDEX
     &,ALPHA1,ASHTF,ASURF,TILE_FRAC,ICE_FRACT
     &,RHOKH1_PRIME,TIMESTEP,SIMLT,SMLT,DFQW,DIFF_SENS_HTF
     &,EI,LYING_SNOW,SURF_HT_FLUX,TSTAR_TILE,TI
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,LTIMER)

!-----------------------------------------------------------------------
! 3. Update heat and moisture fluxes due to limited evaporation and snow
!    or sea-ice melting.
!-----------------------------------------------------------------------

      DO I = P1,P1+POINTS-1

        TSTAR_GB(I) = TSTAR_TILE(I,1)
        EI_GB(I) = EI(I,1)
        SNOWMELT_GB(I)=SNOWMELT(I,1)


        IF ( ICE_FRACT(I).GT.0.0 ) THEN
          DQW_GB(I) = DTRDZ_1(I) * DFQW(I,1)
          DTL_GB(I) = DTRDZ_1(I) * DIFF_SENS_HTF(I,1) / CP
          TL(I,1) = TL(I,1) + DTL_GB(I)
          QW(I,1) = QW(I,1) + DQW_GB(I)
          FTL(I,1) = FTL(I,1) + DIFF_SENS_HTF(I,1)
          FQW(I,1) = EOLD_GB(I) + DFQW(I,1)

          do itile=1,n_types
             ftl_tile(i,itile)=ftl_tile(i,itile) + DIFF_SENS_HTF(I,1)
             fqw_tile(i,itile)=fqw_tile(i,itile) + DFQW(I,1)
          enddo

          D_S_H_GB(I) = DIFF_SENS_HTF(I,1)
          DFQW_GB(I) = DFQW(I,1)

        ENDIF ! ice_fract .gt. 0

        IF ( LAND_MASK(I) ) THEN

          DQW_GB(I) = 0.0
          EI_GB(I) = 0.0
          DTL_GB(I) = 0.0
          D_S_H_GB(I) = 0.0
          DFQW_GB(I) = 0.0
          SNOWMELT_GB(I) = 0.0
          TSTAR_GB(I) = 0.0
        ENDIF  ! land
      ENDDO ! POINTS


      DO ITILE=1,N_TYPES
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
       DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)

          DQW(I) = DTRDZ_1(I) * DFQW(I,ITILE)
          DTL(I) = DTRDZ_1(I) * DIFF_SENS_HTF(I,ITILE) / CP
          EI_GB(I) = EI_GB(I) + EI(I,ITILE) * TILE_FRAC(I,ITILE)
          DQW_GB(I) = DQW_GB(I) + DQW(I) * TILE_FRAC(I,ITILE)
          DTL_GB(I) = DTL_GB(I) + DTL(I) * TILE_FRAC(I,ITILE)
          D_S_H_GB(I) = D_S_H_GB(I) + DIFF_SENS_HTF(I,ITILE) *
     &                                TILE_FRAC(I,ITILE)
          DFQW_GB(I) = DFQW_GB(I) + DFQW(I,ITILE) *
     &                              TILE_FRAC(I,ITILE)

          TSTAR_GB(I) = TSTAR_GB(I) + TSTAR_TILE(I,ITILE) *
     &                                TILE_FRAC(I,ITILE)

          ECAN_GB(I) = ECAN_GB(I) + ECAN(I,ITILE) *
     &                              TILE_FRAC(I,ITILE)

          SNOWMELT_GB(I)=SNOWMELT_GB(I) + SNOWMELT(I,ITILE) *
     &                                    TILE_FRAC(I,ITILE)

          TL(I,1) = TL(I,1) + DTL(I) * TILE_FRAC(I,ITILE)
          QW(I,1) = QW(I,1) + DQW(I) * TILE_FRAC(I,ITILE)

          FTL_TILE(I,ITILE) = FTL_TILE(I,ITILE) +
     &                        DIFF_SENS_HTF(I,ITILE)
          FQW_TILE(I,ITILE) = EOLD(I,ITILE) + DFQW(I,ITILE)

        ENDDO ! land points
      ENDDO ! Tile loop


      DO I=P1,P1+POINTS-1
        IF ( LAND_MASK(I)) THEN
            FTL(I,1) = FTL(I,1) + D_S_H_GB(I)
            FQW(I,1) = EOLD_GB(I) + DFQW_GB(I)
        ENDIF ! Land block
      ENDDO


!-----------------------------------------------------------------------
!!  Apply increments to rapidly mixing layer.
!-----------------------------------------------------------------------

      DO K = 2,BL_LEVELS-1
        KM1 = K - 1
        DO I=P1,P1+POINTS-1

          IF ( LAND_MASK(I) .OR. ICE_FRACT(I).GT.0.0 ) THEN
            IF ( K .LE. NRML(I) ) THEN
              TL(I,K) = TL(I,K) + DTL_GB(I)
              QW(I,K) = QW(I,K) + DQW_GB(I)
              D_S_H_GB(I) = D_S_H_GB(I)
     &                           - CP * DTL_GB(I) / DTRDZ(I,KM1)
              DFQW_GB(I) = DFQW_GB(I) - DQW_GB(I) / DTRDZ(I,KM1)
              FTL(I,K) = FTL(I,K) + D_S_H_GB(I)
              FQW(I,K) = FQW(I,K) + DFQW_GB(I)
            ENDIF  ! Rapidly mixing layer
          ENDIF      ! Land or sea-ice
        ENDDO     ! Loop over points
      ENDDO ! Loop over levels

!-----------------------------------------------------------------------
!! 4. Diagnose temperature and/or specific humidity at screen height
!!    (1.5 metres), as requested via the STASH flags.
!-----------------------------------------------------------------------

      IF (SQ1P5 .OR. ST1P5) THEN
          ITILE=1  ! when using more than 1 tile, use short grass
        IF (SQ1P5) THEN                                                 
          CALL QSAT(QS(P1),TSTAR_TILE(P1,ITILE),PSTAR(P1),POINTS)
          CALL QSAT(QSTAR_GB(P1),TSTAR_GB(P1),PSTAR(P1),POINTS)
        ENDIF
        DO I=P1,P1+POINTS-1

          IF (ST1P5) THEN

            TL_BLEND = TSTAR_GB(I) - G/CP * (H_BLEND(I) - Z0H(I,ITILE))
     &                 + (TL(I,1)
     &                    + G/CP * (Z1_TQ(I)+Z0M(I,ITILE)-Z0H(I,ITILE))
     &                    - TSTAR_GB(I) ) * HEAT_BLEND_FACTOR(I)
     &                 + ( HEAT_BLEND_FACTOR(I) - 1.0 )
     &                 * ( LCRCP*QCL_1(I) + LSRCP*QCF_1(I) )

            T1P5M(I) = TSTAR_TILE(I,ITILE) - GRCP*Z1P5M + CHR1P5M(I) *
     &                   ( TL_BLEND - TSTAR_TILE(I,ITILE)     
     &                     + GRCP * (H_BLEND(I) - Z0H(I,ITILE)) )

!            T1P5M(I) = TSTAR_TILE(I,1) - GRCP*Z1P5M + CHR1P5M(I) *
!     &             ( TL_BLEND - TSTAR_TILE(I,1) +
!     &                    GRCP*(H_BLEND(I)+Z0M(I,1)-Z0H(I,1)) )

          ENDIF ! st1p5
          IF (SQ1P5) THEN
            QW_BLEND = HEAT_BLEND_FACTOR(I) * (QW(I,1) - QSTAR_GB(I)) +
     &                 QSTAR_GB(I) - ( HEAT_BLEND_FACTOR(I) - 1.0 ) *
     &                ( QCL_1(I) + QCF_1(I) )

            Q1P5M(I) = QW_BLEND + CER1P5M(I)*( QW_BLEND - QS(I) )
          ENDIF !sq1p5
        ENDDO ! POINTS
      ENDIF ! sq1p5 or qt1p5

      IF (LTIMER) THEN
        CALL TIMER('SFEVAP  ',4)
      ENDIF

      RETURN
      END
