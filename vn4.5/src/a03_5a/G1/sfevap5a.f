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
C*LL  SUBROUTINE SF_EVAP------------------------------------------------
CLL
CLL  Purpose: Calculate surface evaporation and sublimation amounts
CLL           (without applying them to the surface stores).
CLL           Also calculate heat flux due to sea-ice melting.
CLL           Also calculate 1.5 metre T and Q.
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
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
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
     + P_FIELD,P1,LAND_FIELD,LAND1
     +,POINTS,BL_LEVELS,LAND_MASK,LAND_PTS,LAND_INDEX
     +,ALPHA1,ASURF,ASHTF,CANOPY,CATCH
     +,DTRDZ,DTRDZ_RML,E_SEA,FRACA
     +,ICE_FRACT,NRML,RHOKH_1,SMC,TIMESTEP,CER1P5M,CHR1P5M
     +,PSTAR,RESFS,RESFT,Z1,Z0M,Z0H,SQ1P5,ST1P5,SIMLT,SMLT
     +,FTL,FQW,LYING_SNOW,QW,SURF_HT_FLUX
     +,TL,TSTAR,TI,ECAN,ES,EI
     +,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT
     +,Q1P5M,T1P5M,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + P_FIELD              ! IN No. of gridpoints in the whole grid.
     +,P1                   ! IN 1st P-pt in full field to be processed.
     +,LAND_FIELD           ! IN No. of landpoints in the whole grid.
     +,LAND1                ! IN 1st L-pt in full field to be processed.
     +,POINTS               ! IN No. of gridpoints to be processed.
     +,BL_LEVELS            ! IN No. of levels treated by b.l. scheme.
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
     +,ASURF(P_FIELD)       ! IN Soil coefficient from P242 (sq m K per
C                           !    per Joule * timestep).
     +,ASHTF(P_FIELD)       ! IN Coefficient to calculate
C                           !    the soil heat flux
C                           !    between the surface and top soil
C                           !    layer (W/m2/K)
     +,CANOPY(LAND_FIELD)   ! IN Gridbox mean canopy / surface water
C                           !    store (kg per sq m).
     +,CATCH(LAND_FIELD)    ! IN Canopy / surface water store capacity
C                           !    (kg per sq m).
     +,DTRDZ(P_FIELD,       ! IN -g.dt/dp for each model layer on p-grid
     +       BL_LEVELS)     !    From P244 ((kg/sq m/s)**-1).
     +,DTRDZ_RML(P_FIELD)   ! IN -g.dt/dp for the rapidly mixing layer
C                           !    (if it exists) on the p-grid from P244.
     +,E_SEA(P_FIELD)       ! IN Evaporation from sea (weighted with
C                           !    leads fraction at sea-ice points).
     +,FRACA(P_FIELD)       ! IN Fraction of surface moisture flux
C                           !    with only aerodynamic resistance.
C                           !       Diagnostics defined on land and sea.
     +,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is covered by
C                           !    sea-ice (decimal fraction, but most of
C                           !    this sub-component assumes it to be
C                           !    either 1.0 or 0.0 precisely).
C                           !    NB Dimension is PFIELD not LAND_FIELD f
C                           !    snow on sea-ice in coupled model runs.
     +,SMC(LAND_FIELD)      ! IN Soil moisture content (kg per sq m).
     +,TIMESTEP             ! IN Timestep (sec).
      LOGICAL
     + SQ1P5                ! IN STASH flag for 1.5-metre sp humidity.
     +,ST1P5                ! IN STASH flag for 1.5-metre temperature.
     +,SIMLT                ! IN STASH flag for sea-ice melting ht flux.
     +,SMLT                 ! IN STASH flag for snow melting ht flux.
      REAL
     + CER1P5M(P_FIELD)     ! IN Transfer coefficient ratio, from P243.
     +,CHR1P5M(P_FIELD)     ! IN Transfer coefficient ratio, from P243.
     +,PSTAR(P_FIELD)       ! IN Surface pressure (Pa).
     +,RESFS(P_FIELD)       ! IN Combined soil, stomatal and
C                           !    aerodynamic resistance factor
     +,RESFT(P_FIELD)       ! IN Total resistance factor
C                           !     FRACA+(1-FRACA)*RESFS.
     +,Z1(P_FIELD)          ! IN Height of lowest atmospheric level
C                           !    (i.e. middle of lowest layer).  Metres.
     +,Z0M(P_FIELD)         ! IN Roughness length for momentum (m)
     +,Z0H(P_FIELD)         ! IN Roughness length for heat and moisture
      INTEGER
     & NRML(P_FIELD)        ! IN  The Number of model layers in the
C                           !     Rapidly Mixing Layer.
      REAL
     + RHOKH_1(P_FIELD)     ! IN    Turbulent surface exchange
C                           !       coefficient for sensible heat.
     +,FTL(P_FIELD,         ! INOUT Sensible heat flux from layer k-1 to
     +     BL_LEVELS)       !       layer k (W/sq m).  From P243 and
C                           !      P244, units changed in P24 top level.
     +,FQW(P_FIELD,         ! INOUT Turbulent moisture flux from level
     +     BL_LEVELS)       !       k-1 to k (kg/sq m/s). From P243/4.
C                           !       Diagnostics defined on land and sea.
     +,LYING_SNOW(P_FIELD)  ! INOUT Lying snow (kg per sq m).
     +,QW(P_FIELD,BL_LEVELS)! INOUT Total water content (kg(water)/
C                           !       kg(air)).  From P243/4.
C
     +,SURF_HT_FLUX(P_FIELD)! INOUT Net downward heat flux at surface
C                           !       over land or sea-ice fraction of
C                           !       gridbox (W/m2).
     +,TSTAR(P_FIELD)       ! INOUT Surface temperature (K).
     +,TI(P_FIELD)          ! INOUT Sea-ice surface layer temp. (K).
     +,TL(P_FIELD,BL_LEVELS)! INOUT Liquid/frozen water temperature (K).
      REAL
     + ECAN(P_FIELD)        ! OUT Gridbox mean evaporation from canopy/
C                           !     surface store (kg per sq m per s).
C                           !     Zero over sea and sea-ice.
     +,ES(P_FIELD)          ! OUT Surface evapotranspiration (through
C                           !       a resistance which is not entirely
C                           !       aerodynamic).  Always non-negative.
C                           !       Kg per sq m per sec.
C                           !     Diagnostics defined on land and sea.
     +,EI(P_FIELD)          ! OUT Sublimation from lying snow or sea-
C                           !     ice (kg per sq m per s).
      REAL
     + SICE_MLT_HTF(P_FIELD)! OUT Heat flux due to melting of sea-ice
C                           !     (Watts per square metre).
     +,SNOMLT_SURF_HTF(P_FIELD)! OUT Heat flux due to surface melting
C                              !     of snow (W/m2).
     +,SNOWMELT(P_FIELD)    ! OUT Surface snowmelt (kg/m2/s).
     +,Q1P5M(P_FIELD)       ! OUT Specific humidity at screen height of
C                           !     1.5 metres (kg water per kg air).
     +,T1P5M(P_FIELD)       ! OUT Temperature at 1.5 metres above the
C                           !     surface (K).
C*
C*L  External subprogram(s) required :-
      EXTERNAL QSAT,SF_MELT
      EXTERNAL TIMER
C*
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

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
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
      REAL GRCP
      PARAMETER (
     + GRCP=G/CP   ! Accn due to gravity / standard heat capacity of
C                  ! air at const pressure.  Used in diagnosis of 1.5
C                  ! metre temperature.
     +)
C*
      REAL
     + DFQW(P_FIELD)         ! Adjustment increment to the flux of
C                            ! total water
     +,DIFF_SENS_HTF(P_FIELD)! Adjustment increment to the sensible
C                            ! heat flux
     +,DQW(P_FIELD)          ! Increment to specific humidity
     +,DTL(P_FIELD)          ! Increment to temperature
     +,DTRDZ_1(P_FIELD)      ! -g.dt/dp for surface layer or rml if it
C                            ! exists from P244 ((kg/sq m/s)**-1).
     +,EOLD(P_FIELD)         ! Used to store initial value of evap.
C                            ! from P244
     +,EW(P_FIELD)           ! Total surface flux of water, excluding
C                            ! sublimation/frost deposition, over land.
     +,LEOLD(P_FIELD)        ! Used to store initial value of latent
C                            ! heat flux from P244
     +,QS(P_FIELD)           ! Used for saturated specific humidity
C                            ! at surface, in Q1P5M calculation.
     +,RHOKH1_PRIME(P_FIELD) ! Modified forward time-weighted transfer
C                            ! coefficient
C  Local scalars
      REAL
     + DIFF_LAT_HTF        ! Increment to the latent heat flux.
     +,DIFF_SURF_HTF       ! Increment to the surface heat flux.
     +,DTSTAR              ! Increment for surface temperature.
     +,EA                  ! Surface evaporation with only aero-
C                          ! dynamic resistance (+ve), or condens-
C                          ! ation (-ve), averaged over gridbox
C                          ! (kg/m2/s).
     +,EADT                ! EA (q.v.) integrated over timestep.
     +,ECANDT              ! ECAN (q.v.) integrated over timestep.
     +,EDT                 ! E=FQW(,1) (q.v.) integrated over timestep.
     +,EIDT                ! EI (q.v.) integrated over timestep.
     +,ESDT                ! ES (q.v.) integrated over timestep.
     +,ESL                 ! ES (q.v.) without fractional weighting
C                          ! factor FRACS ('L' is for 'local')
C                          ! (kg/m2/s).
     +,ESLDT               ! ESL (q.v.) integrated over timestep.
     +,FRACS               ! Fraction of gridbox at which moisture flux
C                          ! is additionally impeded by a surface and/or
C                          ! stomatal resistance.
      INTEGER
     + I                   ! Loop counter - full horizontal field index.
     +,L                   ! Loop counter - land field index.
     +,K                   ! Loop counter in the vertical.
     +,KM1                 ! K - 1
C
C-----------------------------------------------------------------------
CL 1. Initialise some output variables and flux increments to zero.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
      CALL TIMER('SFEVAP  ',3)
      ENDIF
      DO 1 I=P1,P1+POINTS-1
        ECAN(I) = 0.0
        EI(I) = 0.0
        DIFF_SENS_HTF(I) = 0.0
        DFQW(I) = 0.0
    1 CONTINUE
C
C---------------------------------------------------------------------
CL 2. Do calculations for land points.
C---------------------------------------------------------------------
C
CMIC$ DO ALL VECTOR SHARED(P_FIELD, LAND_FIELD, BL_LEVELS, LAND1,
CMIC$1   LAND_PTS, LAND_INDEX, ESL, TIMESTEP, ES, LYING_SNOW, ECAN,
CMIC$2   EA, CATCH, CANOPY, SMC, EI, TSTAR, FQW, EOLD, LEOLD,
CMIC$3   P1,POINTS,LC,LF,TM,LAND_MASK,EW,FRACA,RESFT,RESFS)
CMIC$4   PRIVATE(I, L, ESLDT,
CMIC$5   ESDT, EADT, EDT, ECANDT, FRACS, EIDT)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 2 L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          IF (FQW(I,1).EQ.0.0) THEN
            EA = 0.0
            ESL = 0.0
          ELSE
            EA = FQW(I,1) / RESFT(I) * FRACA(I)
            ESL = FQW(I,1) / RESFT(I) * RESFS(I)
          END IF
          ES(I) = ESL * (1. - FRACA(I))
C
C-----------------------------------------------------------------------
CL 2.1 Calculate fluxes integrated over timestep.
C-----------------------------------------------------------------------
C
            ESLDT = ESL * TIMESTEP
            EADT = EA * TIMESTEP
            ESDT = ES(I) * TIMESTEP
            EDT = EADT + ESDT
C
C-----------------------------------------------------------------------
CL 2.2 Do calculations for snow-free land.  Canopy processes operate.
CL     LYING_SNOW is defined on sea and land points for snow on sea-ice
CL     in coupled model runs.
C-----------------------------------------------------------------------
C
            IF (LYING_SNOW(I).LE.0.0) THEN
C
C**********************************************************************
C Store initial value of evaporation and latent heat flux
C**********************************************************************
C
              EOLD(I) = FQW(I,1)
              LEOLD(I) = FQW(I,1) * LC
              IF (EDT.GE.0.0) THEN
C
C-----------------------------------------------------------------------
CL 2.2.1 Non-negative moisture flux over snow-free land.
C-----------------------------------------------------------------------
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL   (a) Water in canopy and soil is assumed to be liquid, so all
CL       positive moisture flux over snow-free land is evaporation
CL       rather than sublimation, even if TSTAR is less than or equal
CL       to TM.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
                ECAN(I) = EA
                ECANDT = EADT
C
C  If EDT is non-negative, then ECANDT must be non-negative.
C
                FRACA(I) = 0.0
                IF (CATCH(L).GT.0.0)
     +            FRACA(I) = CANOPY(L) / CATCH(L)
                IF (CANOPY(L).LT.ECANDT) THEN
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL   (b) It is assumed that any 'canopy' moisture flux in excess of the
CL       current canopy water amount is in fact soil evaporation.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C        This situation is highly improbable - it will occur at, at
C        most, a few gridpoints in any given timestep.
C
                  FRACS = 1.0 - FRACA(I)*( CANOPY(L) / ECANDT )
                  ESDT = ESLDT * FRACS
                  ECANDT = CANOPY(L)
                  ECAN(I) = ECANDT / TIMESTEP
                  ES(I) = ESDT / TIMESTEP
                ENDIF
C
C  (The canopy store is depleted by evaporation in P252, and not here,
C   according to the formula: CANOPY=CANOPY-ECANDT)
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL   (c) Adjustments to evaporation from soil as calculated so far :-
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
              IF (SMC(L).LE.0.0) THEN
C
CL   (i) If there is currently no soil moisture, there must be no
CL       evaporation of soil moisture, so this flux is set to zero.
C
                  ESDT = 0.0
                  ES(I) = 0.0
                ELSEIF (SMC(L).LT.ESDT) THEN
C
CL  (ii) Ensure that the soil evaporation is not greater than the
CL       current soil moisture store.
C        This situation is extremely unlikely at any given gridpoint
C        at any given timestep.
C
                  ESDT = SMC(L)
                  ES(I) = ESDT / TIMESTEP
                ENDIF
C
C  (The soil moisture store is depleted by evaporation in P253, and not
C   here, using the formula:  SMC=SMC-ESDT)
C
                EW(I) = ECAN(I) + ES(I)
                EI(I) = 0.0
C
C-----------------------------------------------------------------------
CL 2.2.2 Negative moisture flux onto snow-free land above freezing
C-----------------------------------------------------------------------
CL       (i.e. condensation onto snow-free land).  The whole flux is
CL       into the surface/canopy store.
C
              ELSEIF (TSTAR(I).GT.TM) THEN     ! ELSE of evaporation /
C                                              ! condensation block.
C
C  Condensation implies ES=0, so ECAN=EA=EW=E (=FQW(,1))
C
                ECAN(I) = FQW(I,1)
                ES(I) = 0.0
                EW(I) = ECAN(I)
                EI(I) = 0.0
C
C  (The canopy store is augmented by interception of condensation at
C   P252, and not here.)
C
C-----------------------------------------------------------------------
CL 2.2.3 Negative moisture flux onto snow-free land below freezing
CL       (i.e. deposition of frost).
C-----------------------------------------------------------------------
C
              ELSE      ! ELSE of condensation / frost deposition block.
                EI(I) = FQW(I,1)
                ES(I) = 0.0
                EW(I) = 0.0
C
C  (Negative EI is used to increment the snowdepth store - there is
C   no separate "frost" store.  This incrementing is done in P251,
C   according to:  LYING_SNOW = LYING_SNOW - EI*TIMESTEP)
C
              ENDIF  ! End of evaporation/condensation/deposition block.
C
C-----------------------------------------------------------------------
CL 2.3 Do calculations for snow-covered land.
C-----------------------------------------------------------------------
C
            ELSEIF (LYING_SNOW(I).LE.EDT) THEN     ! ELSEIF of no-snow.
C
C**********************************************************************
C Store initial value of evaporation and latent heat flux
C**********************************************************************
C
              EOLD(I) = FQW(I,1)
              LEOLD(I) = FQW(I,1) * ( LC + LF )
C
C-----------------------------------------------------------------------
CL 2.3.1 Shallow snow (lying snow or frost which is being exhausted
CL       by evaporation).  All the snow is sublimated, the remaining
CL       moisture flux being taken from the canopy and soil, with all
CL       the palaver of section 1.2.1 above.
C-----------------------------------------------------------------------
C
C        This is extremely unlikely at more than one or two gridpoints
C        at any given timestep, yet the complicated logic probably
C        slows down the routine considerably - this section is a
C        suitable candidate for further consideration as regards
C        making the model optimally efficient.
C
              EI(I) = LYING_SNOW(I) / TIMESTEP
              EIDT = LYING_SNOW(I)
C
C  Set EDT = ( E - SNOSUB ) * TIMESTEP.  This is the moisture in kg per
C  square metre left over to be evaporated from the canopy and soil.
C  N.B.  E=FQW(,1)
C
              EDT = EDT - EIDT
C
C  (Snowdepth is decreased using EI at P251, and not here.  The formula
C   used is simply:  LYING_SNOW = LYING_SNOW - EI*TIMESTEP.)
C
C  Now that all the snow has sublimed, canopy processes come into
C  operation (FRACA no longer necessarily equal to 1).
C
              FRACA(I) = 0.0
              IF (CATCH(L).GT.0.0)
     +          FRACA(I) = CANOPY(L) / CATCH(L)
              ECANDT = EDT * FRACA(I)
              IF (CANOPY(L).LT.ECANDT) THEN
C
C  Dry out the canopy completely and assume the remaining moisture flux
C  is soil evaporation.
C
                FRACS = 1.0 - FRACA(I)*( CANOPY(L) / ECANDT )
                ESDT = EDT * FRACS
                ECANDT = CANOPY(L)
              ELSE
C
C  Calculate soil evaporation.
C
                FRACS = 1.0 - FRACA(I)
                ESDT = EDT * FRACS
              ENDIF
              ECAN(I) = ECANDT / TIMESTEP
              ES(I) = ESDT / TIMESTEP
C
C  (ECAN is used to deplete the canopy store at P252, and not here.  The
C   formula used is simply:  CANOPY = CANOPY - ECAN*TIMESTEP.)
C
C  Evaporation from soil.
C
              IF (SMC(L).LE.0.0) THEN
C
C  No evaporation from soil possible when there is no soil moisture.
C
                ESDT = 0.0
                ES(I) = 0.0
              ELSEIF (SMC(L).LT.ESDT) THEN
C
C  Limit evaporation of soil moisture in the extremely unlikely event
C  that soil moisture is exhausted by the evaporation left over from
C  sublimation which exhausted the snow store.
C
                ESDT = SMC(L)
                ES(I) = ESDT / TIMESTEP
              ENDIF
C
C  (ES is used to deplete the soil moisture store at P253, and not here,
C   according to the formula:  SMC = SMC - ES*TIMESTEP.)
C
              EW(I) = ECAN(I) + ES(I)
C
C-----------------------------------------------------------------------
CL 2.3.2 Deep snow (i.e. not being exhausted by evaporation).  This
CL       covers two cases: (a) sublimation from deep snow (if total
CL       moisture flux over the timestep is non-negative but less than
CL       the lying snow amount), and (b) deposition onto an already
CL       snowy surface (if the total moisture flux is negative and
CL       the lying snow amount is positive).
C-----------------------------------------------------------------------
C
            ELSE          ! ELSE of shallow snow / deep snow block.
              EI(I) = FQW(I,1)
              EW(I) = 0.0
C
C**********************************************************************
C Store initial value of evaporation and latent heat flux
C**********************************************************************
C
              EOLD(I) = FQW(I,1)
              LEOLD(I) = FQW(I,1) * ( LC + LF )
C
C  (EI is used to increase or decrease the snowdepth at P251, and not
C   here, according to the formula:
C            LYING_SNOW = LYING_SNOW - EI*TIMESTEP . )
C
            ENDIF         ! End of no snow/shallow snow/deep snow block.
            FQW(I,1) = EW(I) + EI(I)
    2   CONTINUE
C
C  Split loop 2 here so that it will vectorise.
C
CMIC$ DO ALL VECTOR SHARED(DTRDZ_1,DTRDZ,RHOKH_1,GAMMA,
CMIC$1 NRML,DTRDZ_RML,EI,EW,LEOLD,DIFF_LAT_HTF,FQW,
CMIC$2 EOLD,DFQW,ASHTF,DIFF_SENS_HTF,DIFF_SURF_HTF,
CMIC$3 ASURF,TIMESTEP,TSTAR,LAND_INDEX,RHOKH1_PRIME,
CMIC$4 SURF_HT_FLUX) PRIVATE(DTSTAR,I,L)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 24 L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
C
C***********************************************************************
C  2.4 Calculate increments to surface and subsurface temperatures,
C      surface heat and moisture fluxes and soil heat flux. Apply
C      increments to TSTAR to give interim values before any
C      snowmelt.
C***********************************************************************
          IF (NRML(I).GE.2) THEN
            DTRDZ_1(I) = DTRDZ_RML(I)
          ELSE
            DTRDZ_1(I) = DTRDZ(I,1)
          ENDIF
          RHOKH1_PRIME(I) = 1.0 / ( 1.0 / RHOKH_1(I)
     &                                  + GAMMA(1) * DTRDZ_1(I) )
          DIFF_LAT_HTF = (LC + LF) * EI(I) + LC * EW(I) - LEOLD(I)
          DFQW(I) = FQW(I,1) - EOLD(I)
          DIFF_SENS_HTF(I) = - DIFF_LAT_HTF /
     &               ( 1. + ASHTF(I) /(RHOKH1_PRIME(I) * CP) )
          DIFF_SURF_HTF = - DIFF_LAT_HTF / ( 1.0 +
     &                      RHOKH1_PRIME(I) * CP / ASHTF(I) )
          SURF_HT_FLUX(I) = SURF_HT_FLUX(I) + DIFF_SURF_HTF
          DTSTAR = DIFF_SURF_HTF / ASHTF(I)
          TSTAR(I) = TSTAR(I) + DTSTAR
   24   CONTINUE
C
C-----------------------------------------------------------------------
CL 2.5 Do calculations for sea points.
C-----------------------------------------------------------------------
C
CMIC$ DO ALL VECTOR SHARED(P_FIELD, BL_LEVELS, P1, POINTS,NRML,
CMIC$1  LAND_MASK, ES, ECAN, EI, EOLD,
CMIC$2  ICE_FRACT, FQW, E_SEA,DTRDZ_RML,
CMIC$3  TSTAR, SMLT, SICE_MLT_HTF, KAPPAI,
CMIC$4  DTRDZ_1,DTRDZ,RHOKH_1,GAMMA,RHOKH1_PRIME,
CMIC$5  TIMESTEP,TM,TFS) PRIVATE(I, TSTARMAX)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 25 I=P1,P1+POINTS-1
          IF (.NOT.LAND_MASK(I)) THEN
C
C-----------------------------------------------------------------------
CL 2.5.1 Set soil and canopy evaporation amounts to zero, and set
CL       sublimation to zero for liquid sea points.
C-----------------------------------------------------------------------
C
            ES(I) = 0.0
            ECAN(I) = 0.0
            EI(I) = 0.0
C-----------------------------------------------------------------------
CL 2.5.3 For sea-ice points :-
C-----------------------------------------------------------------------
C
            IF (ICE_FRACT(I).GT.0.0) THEN
              EOLD(I) = FQW(I,1)
              EI(I) = FQW(I,1) - E_SEA(I)
              IF (NRML(I).GE.2) THEN
                DTRDZ_1(I) = DTRDZ_RML(I)
              ELSE
                DTRDZ_1(I) = DTRDZ(I,1)
              ENDIF
              RHOKH1_PRIME(I) = 1.0 / ( 1.0 / RHOKH_1(I)
     &                          + ICE_FRACT(I)*GAMMA(1)*DTRDZ_1(I) )
            ENDIF     ! End of liquid sea/sea-ice block.
          ENDIF       ! End of sea point calculations.
   25   CONTINUE
C
C-----------------------------------------------------------------------
C  Calculate fluxes and increments associated with melting of snow
C  or sea-ice.
C-----------------------------------------------------------------------
      CALL SF_MELT(P_FIELD,P1,LAND_FIELD,LAND1
     +,POINTS,LAND_MASK,LAND_PTS,LAND_INDEX
     +,ALPHA1,ASHTF,ASURF,ICE_FRACT
     +,RHOKH1_PRIME,TIMESTEP,SIMLT,SMLT,DFQW,DIFF_SENS_HTF
     +,EI,LYING_SNOW,SURF_HT_FLUX,TSTAR,TI
     +,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,LTIMER)
C
C-----------------------------------------------------------------------
C 3. Update heat and moisture fluxes due to limited evaporation and snow
C    or sea-ice melting.
C-----------------------------------------------------------------------
C
      DO I = P1,P1+POINTS-1
        IF ( LAND_MASK(I) .OR. ICE_FRACT(I).GT.0.0 ) THEN
          DQW(I) = DTRDZ_1(I) * DFQW(I)
          DTL(I) = DTRDZ_1(I) * DIFF_SENS_HTF(I) / CP
          TL(I,1) = TL(I,1) + DTL(I)
          QW(I,1) = QW(I,1) + DQW(I)
          FTL(I,1) = FTL(I,1) + DIFF_SENS_HTF(I)
          FQW(I,1) = EOLD(I) + DFQW(I)
        ENDIF                    ! LAND_MASK etc.
      ENDDO                      ! P1+POINTS-1
C-----------------------------------------------------------------------
C     Apply increments to rapidly mixing layer.
C-----------------------------------------------------------------------
      DO K = 2,BL_LEVELS-1
        KM1 = K - 1
        DO I = P1,P1+POINTS-1
          IF ( LAND_MASK(I) .OR. ICE_FRACT(I).GT.0.0 ) THEN
            IF ( K .LE. NRML(I) ) THEN
              TL(I,K) = TL(I,K) + DTL(I)
              QW(I,K) = QW(I,K) + DQW(I)
              DIFF_SENS_HTF(I) = DIFF_SENS_HTF(I)
     &                           - CP * DTL(I) / DTRDZ(I,KM1)
              DFQW(I) = DFQW(I) - DQW(I) / DTRDZ(I,KM1)
              FTL(I,K) = FTL(I,K) + DIFF_SENS_HTF(I)
              FQW(I,K) = FQW(I,K) + DFQW(I)
            ENDIF  ! Rapidly mixing layer
          ENDIF                 ! Land or sea-ice
        ENDDO                   ! Loop over points
      ENDDO                     ! Loop over levels
C
C-----------------------------------------------------------------------
CL 4. Diagnose temperature and/or specific humidity at screen height
CL    (1.5 metres), as requested via the STASH flags.
C-----------------------------------------------------------------------
C
      IF (SQ1P5 .OR. ST1P5) THEN
        IF (SQ1P5) CALL QSAT(QS(P1),TSTAR(P1),PSTAR(P1),POINTS)
        DO 4 I=P1,P1+POINTS-1
          IF (ST1P5) T1P5M(I) = TSTAR(I) - GRCP*Z1P5M + CHR1P5M(I) *
     +               ( TL(I,1) - TSTAR(I) + GRCP*(Z1(I)+Z0M(I)-Z0H(I)) )
          IF (SQ1P5) Q1P5M(I) = QW(I,1) + CER1P5M(I)*( QW(I,1) - QS(I) )
    4   CONTINUE
      ENDIF
      IF (LTIMER) THEN
      CALL TIMER('SFEVAP  ',4)
      ENDIF
      RETURN
      END
