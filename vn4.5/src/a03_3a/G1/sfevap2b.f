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
C*LL  SUBROUTINE SF_EVAP------------------------------------------------
CLL
CLL  Purpose: Calculate surface evaporation and sublimation amounts
CLL           (without applying them to the surface stores).
CLL           Also calculate heat flux due to sea-ice melting.
CLL           Also calculate 1.5 metre T and Q.
CLL
CLL         Change to the calculation of soil evaporation when
CLL    lying snow disappears during a timestep. Also two new
CLL    common decks C_HT_M and C_SICEHC
CLL
CLL  Suitable for single column usage.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   3.1   11/1/93   Co-ordinate transformation so that screen
CLL                   temperature is calculated at 1.5m above
CLL                   surface + 1.5m            S.Jackson
CLL   3.4   06/06/94  DEF TIMER replaced by LOGICAL LTIMER
CLL                   Argument LTIMER added
CLL                                                 S.J.Swarbrick
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
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
     + P_FIELD,P1,LAND_FIELD,LAND1,
     + POINTS,BL_LEVELS,LAND_MASK,LAND_PTS,LAND_INDEX
     +,ASOIL_1,CANOPY,CATCH,DTRDZ,DTRDZ_RML,EA,ESL,SOIL_HT_FLUX
     +,E_SEA,ICE_FRACT,LYING_SNOW,RADNET,SMC,TIMESTEP
     +,CER1P5M,CHR1P5M,PSTAR,Z1,Z0M,Z0H,SQ1P5,ST1P5,SMLT
     +,NRML,RHOKH_1,DQW_1,DQW_RML,DTSTAR,FTL,FQW,ES,QW,TL,TSTAR
     +,ECAN,EI
     +,SICE_MLT_HTF,Q1P5M,T1P5M,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + P_FIELD              ! IN No. of gridpoints in the whole grid.
     +,P1                   ! IN 1st P-pt in full field to be processed.
     +,LAND_FIELD           ! IN No. of LANDpoints in the whole grid.
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
     + ASOIL_1(P_FIELD)     ! IN Soil coefficient from P242 (sq m K per
C                           !    per Joule * timestep).
     +,CANOPY(LAND_FIELD)   ! IN Gridbox mean canopy / surface water
C                           !    store (kg per sq m).
     +,CATCH(LAND_FIELD)    ! IN Canopy / surface water store capacity
C                           !    (kg per sq m).
     +,DTRDZ(P_FIELD,       ! IN -g.dt/dp for each model layer on p-grid
     +       BL_LEVELS)     !    From P244 ((kg/sq m/s)**-1).
     +,DTRDZ_RML(P_FIELD)   ! IN -g.dt/dp for the rapidly mixing layer
C                           !    (if it exists) on the p-grid from P244.
     +,EA(P_FIELD)          ! IN Surface evaporation with only aero-
C                           !    dynamic resistance (+ve), or condens-
C                           !    ation (-ve), averaged over gridbox.
C                           !    From P243 and P244.  Kg per sq m per s.
     +,ESL(P_FIELD)         ! IN ES (q.v.) without fractional weighting
C                           !    factor FRACS ('L' is for 'local').
C                           !    From P243 and P244.  Kg per sq m per s.
     +,SOIL_HT_FLUX(P_FIELD) ! IN Heat flux from surface to deep soil
C                           !    layer 1 (Watts per sq m).  From P242.
     +,E_SEA(P_FIELD)       ! IN Evaporation from sea (weighted with
C                           !    leads fraction at sea-ice points).
C                           !       Diagnostics defined on land and sea.
     +,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is covered by
C                           !    sea-ice (decimal fraction, but most of
C                           !    this sub-component assumes it to be
C                           !    either 1.0 or 0.0 precisely).
     +,LYING_SNOW(P_FIELD)  ! IN Lying snow (kg per sq m).
C                           !    NB Dimension is PFIELD not LAND_FIELD f
C                           !    snow on sea-ice in coupled model runs.
     +,RADNET(P_FIELD)      ! IN Surface net radiation (Watts per sq m).
     +,SMC(LAND_FIELD)      ! IN Soil moisture content (kg per sq m).
     +,TIMESTEP             ! IN Timestep (sec).
      LOGICAL
     + SQ1P5                ! IN STASH flag for 1.5-metre sp humidity.
     +,ST1P5                ! IN STASH flag for 1.5-metre temperature.
     +,SMLT                 ! IN STASH flag for sea-ice melting ht flux.
      REAL
     + CER1P5M(P_FIELD)     ! IN Transfer coefficient ratio, from P243.
     +,CHR1P5M(P_FIELD)     ! IN Transfer coefficient ratio, from P243.
     +,PSTAR(P_FIELD)       ! IN Surface pressure (Pa).
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
     +,DQW_1(P_FIELD)       ! IN    Increment for lowest-level total
C                           !       water content.  From P244.
     +,DQW_RML(P_FIELD)     ! IN    Increment to QW due to rapid mixing.
     +,DTSTAR(P_FIELD)      ! IN    Increment for surface temperature.
C                           !       From P244.
     +,FTL(P_FIELD,         ! INOUT Sensible heat flux from layer k-1 to
     +     BL_LEVELS)       !       layer k (W/sq m).  From P243 and
C                           !      P244, units changed in P24 top level.
     +,FQW(P_FIELD,         ! INOUT Turbulent moisture flux from level
     +     BL_LEVELS)       !       k-1 to k (kg/sq m/s). From P243/4.
     +,ES(P_FIELD)          ! INOUT Surface evapotranspiration (through
C                           !       a resistance which is not entirely
C                           !       aerodynamic).  Always non-negative.
C                           !       Kg per sq m per sec.  From P243/4.
C                           !       Diagnostics defined on land and sea.
     +,QW(P_FIELD,BL_LEVELS)! INOUT Total water content (kg(water)/
C                           !       kg(air)).  From P243/4.
     +,TL(P_FIELD,BL_LEVELS)! INOUT Liquid/frozen water temperature (K).
     +,TSTAR(P_FIELD)       ! INOUT Surface temperature (K).
      REAL
     + ECAN(P_FIELD)        ! OUT Gridbox mean evaporation from canopy/
C                           !     surface store (kg per sq m per s).
C                           !     Zero over sea and sea-ice.
C                           !     Diagnostics defined on land and sea.
     +,EI(P_FIELD)          ! OUT Sublimation from lying snow or sea-
C                           !     ice (kg per sq m per s).
      REAL
     + SICE_MLT_HTF(P_FIELD)! OUT Heat flux due to melting of sea-ice
C                           !     (Watts per square metre).
     +,Q1P5M(P_FIELD)       ! OUT Specific humidity at screen height of
C                           !     1.5 metres (kg water per kg air).
     +,T1P5M(P_FIELD)       ! OUT Temperature at 1.5 metres above the
C                           !     surface (K).
C*
C*L  External subprogram(s) required :-
      EXTERNAL QSAT
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

C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------
      REAL GRCP
      PARAMETER (
     + GRCP=G/CP   ! Accn due to gravity / standard heat capacity of
C                  ! air at const pressure.  Used in diagnosis of 1.5
C                  ! metre temperature.
     +)
C*
      REAL
     + QS(P_FIELD)           ! Used for saturated specific humidity
C                            ! at surface, in Q1P5M calculation.
      REAL
     + DTL_RML(P_FIELD)      ! Adjustment increment to TL in the
C                            ! rapidly mixing layer.
     +,DFTL(P_FIELD,         ! Adjustment increment to the sensible
     +      BL_LEVELS)       ! heat flux.
     +,DFQW(P_FIELD,         ! Adjustment increment to the flux of
     +      BL_LEVELS)       ! total water.
C  Local scalars
      REAL
     + EADT                ! EA (q.v.) integrated over timestep.
     +,ECANDT              ! ECAN (q.v.) integrated over timestep.
     +,EDT                 ! E=FQW(,1) (q.v.) integrated over timestep.
     +,EIDT                ! EI (q.v.) integrated over timestep.
     +,ESDT                ! ES (q.v.) integrated over timestep.
     +,ESLDT               ! ESL (q.v.) integrated over timestep.
     +,EW                  ! Total surface flux of water, excluding
C                          ! sublimation/frost deposition, over land.
     +,FRACA               ! Fraction of gridbox at which moisture flux
C                          ! is going at the potential rate, i.e. with
C                          ! only aerodynamic resistance.
     +,FRACS               ! Fraction of gridbox at which moisture flux
C                          ! is additionally impeded by a surface and/or
C                          ! stomatal resistance.
     +,CT_1                ! Parameter in the implicit adjustment of
C                          ! TSTAR and TL(,1).
     +,CT_RML              ! Parameter in the implicit adjustment of
C                          ! TSTAR and TL in the rapidly mixing layer.
     +,TSTARMAX            ! Maximum gridbox mean surface temperature in
C                          ! sea-ice melting calculation.  See comments
C                          ! to section 1.6.3(b) or 2.5.3(b) below.
      INTEGER
     + I                   ! Loop counter - full horizontal field index.
     +,L                   ! Loop counter - land field index.
     +,K                   ! Loop counter in the vertical.
     +,KM1                 ! K - 1
     +,NRMLP1              ! NRML + 1
C
C-----------------------------------------------------------------------
CL 1. Initialise some output variables to zero.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
      CALL TIMER('SFEVAP  ',3)
      ENDIF
      DO 1 I=P1,P1+POINTS-1
        ECAN(I) = 0.0
        EI(I) = 0.0
        IF (SMLT)
     +    SICE_MLT_HTF(I) = 0.0
    1 CONTINUE
C
C---------------------------------------------------------------------
CL 2. Do calculations for land points.
C---------------------------------------------------------------------
C
CMIC$ DO ALL VECTOR SHARED(P_FIELD, LAND_FIELD, BL_LEVELS, LAND1,
CMIC$1   LAND_PTS, LAND_INDEX, ESL, TIMESTEP, ES, LYING_SNOW, ECAN,
CMIC$2   EA, CATCH, CANOPY, SMC, EI, TSTAR, FQW) PRIVATE(I, L, ESLDT,
CMIC$3   ESDT, EADT, EDT, ECANDT, FRACA, FRACS, EW, EIDT)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 2 L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
C
C-----------------------------------------------------------------------
CL 2.1 Calculate fluxes integrated over timestep.
C-----------------------------------------------------------------------
C
            ESLDT = ESL(I) * TIMESTEP
            EADT = EA(I) * TIMESTEP
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
                ECAN(I) = EA(I)
                ECANDT = EADT
C
C  If EDT is non-negative, then ECANDT must be non-negative.
C
                FRACA = 0.0
                IF (CATCH(L).GT.0.0)
     +            FRACA = CANOPY(L) / CATCH(L)
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
                  FRACS = 1.0 - FRACA*( CANOPY(L) / ECANDT )
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
                EW = ECAN(I) + ES(I)
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
                EW = ECAN(I)
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
                EW = 0.0
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
              FRACA = 0.0
              IF (CATCH(L).GT.0.0)
     +          FRACA = CANOPY(L) / CATCH(L)
              ECANDT = EDT * FRACA
              IF (CANOPY(L).LT.ECANDT) THEN
C
C  Dry out the canopy completely and assume the remaining moisture flux
C  is soil evaporation.
C
                FRACS = 1.0 - FRACA*( CANOPY(L) / ECANDT )
                ESDT = EDT * FRACS
                ECANDT = CANOPY(L)
              ELSE
C
C  Calculate soil evaporation.
C
                FRACS = 1.0 - FRACA
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
              EW = ECAN(I) + ES(I)
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
              EW = 0.0
C
C  (EI is used to increase or decrease the snowdepth at P251, and not
C   here, according to the formula:
C            LYING_SNOW = LYING_SNOW - EI*TIMESTEP . )
C
            ENDIF         ! End of no snow/shallow snow/deep snow block.
            FQW(I,1) = EW + EI(I)
    2   CONTINUE
C
C  Split loop 2 here so that it will vectorise.
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 24 L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
C
C-----------------------------------------------------------------------
CL 2.4 Do calculations which are the same at all land points.  Calculate
CL     and apply increments to (a) total atmospheric water content QW,
CL     (b) surface temperature TSTAR, (c) TL
C-----------------------------------------------------------------------
C
            DTSTAR(I) = ASOIL_1(I) * ( RADNET(I) - FTL(I,1)
     &                  - LC*FQW(I,1) - LF*EI(I) - SOIL_HT_FLUX(I) )
     &                  - DTSTAR(I)
            IF ( NRML(I) .GE. 2 ) THEN
              NRMLP1 = NRML(I) + 1
              DQW_RML(I) = -DTRDZ_RML(I) * ( FQW(I,NRMLP1) - FQW(I,1) )
     &                       -DQW_RML(I)
              DFQW(I,1) = DQW_RML(I) / DTRDZ_RML(I)
              CT_RML = -DTRDZ_RML(I) * RHOKH_1(I)
              DTSTAR(I) = DTSTAR(I) /
     &                     ( 1.0 - CT_RML + ASOIL_1(I)*CP*RHOKH_1(I) )
              DFTL(I,1) = CP * RHOKH_1(I) * DTSTAR(I)
              FTL(I,1) = FTL(I,1) + DFTL(I,1)
              DTL_RML(I) = -CT_RML * DTSTAR(I)
              DTSTAR(I) = DTSTAR(I) + DTL_RML(I)
            ELSE
              DQW_1(I) = -DTRDZ(I,1) * ( FQW(I,2) - FQW(I,1) )
     &                       - DQW_1(I)
              QW(I,1) = QW(I,1) + DQW_1(I)
              CT_1 = -DTRDZ(I,1) * RHOKH_1(I)
              DTSTAR(I) = DTSTAR(I) /
     &                      ( 1.0 - CT_1 + ASOIL_1(I)*CP*RHOKH_1(I) )
              FTL(I,1) = FTL(I,1) + CP * RHOKH_1(I) * DTSTAR(I)
              TL(I,1) = TL(I,1) - CT_1 * DTSTAR(I)
              DTSTAR(I) = (1.0 - CT_1) * DTSTAR(I)
            ENDIF
            TSTAR(I) = TSTAR(I) + DTSTAR(I)
   24   CONTINUE
C
C-----------------------------------------------------------------------
CL 2.5 Do calculations for sea points.
C-----------------------------------------------------------------------
C
CMIC$ DO ALL VECTOR SHARED(P_FIELD, BL_LEVELS, P1, POINTS,
CMIC$1  LAND_MASK, ES, ECAN, EI, ICE_FRACT, FQW, E_SEA,
CMIC$2  TSTAR, SMLT, SICE_MLT_HTF, TIMESTEP) PRIVATE(I, TSTARMAX)
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
C
CL   (a) Set sublimation to total moisture flux less evaporation
CL       from sea.
C
              EI(I) = FQW(I,1) - E_SEA(I)
C
CL   (b) If the surface temperature is above the melting point of
CL       sea-ice, it is reset to be equal to the melting point.
C        (Strictly speaking, the gridbox mean surface temperature is
C        adjusted so that the temperature of the icy fraction is not
C        above the melting point; unlike the rest of this brick, this
C        sub-section is valid for any ice fraction between 1 and 0
C        inclusive, and not just for fractions of precisely 1 or 0.)
CL       This corresponds to a melting of some of the sea-ice.  If
CL       requested, the heat flux associated with this melting is
CL       calculated.  This code is indifferent to whether the melting
CL       reduces the thickness or areal extent of the sea-ice - neither
CL       is altered here.
C        (This too works correctly for a general ice fraction, by
C        weighting the heat flux with the ice fraction, to give the
C        gridbox mean heat flux.)
C
              TSTARMAX = ICE_FRACT(I)*TM + (1.0-ICE_FRACT(I))*TFS
              IF (TSTAR(I).GT.TSTARMAX) THEN
                IF (SMLT)
     +            SICE_MLT_HTF(I) = (TSTAR(I)-TSTARMAX) / (AI*TIMESTEP)
                TSTAR(I) = TSTARMAX
              ENDIF   ! End of TSTAR adjustment block.
            ENDIF     ! End of liquid sea/sea-ice block.
          ENDIF       ! End of sea point calculations.
   25   CONTINUE
C
C-----------------------------------------------------------------------
CL 3. Update heat and moisture fluxes over land due to a rapidly
CL    layer.
C-----------------------------------------------------------------------
C
      DO 3 K = 1,BL_LEVELS-1
        KM1 = K - 1
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          IF ( (K .LE. NRML(I)) .AND. (NRML(I) .GE. 2) ) THEN
            TL(I,K) = TL(I,K) + DTL_RML(I)
            QW(I,K) = QW(I,K) + DQW_RML(I)
            IF ( K .GE. 2 ) THEN
              DFTL(I,K) = DFTL(I,KM1) - CP * DTL_RML(I) / DTRDZ(I,KM1)
              DFQW(I,K) = DFQW(I,KM1) - DQW_RML(I) / DTRDZ(I,KM1)
              FTL(I,K) = FTL(I,K) + DFTL(I,K)
              FQW(I,K) = FQW(I,K) + DFQW(I,K)
            ENDIF ! K .GE. 2
          ENDIF  ! Rapidly mixing layer exists and is more than one
C                ! model layer deep and current layer is in it.
        ENDDO ! Loop over points
    3 CONTINUE! Loop over levels
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
