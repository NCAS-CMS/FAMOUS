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
C*LL  SUBROUTINE SF_RIB------------------------------------------------
CLL
CLL  Purpose: Calculate bulk Richardson number for surface layer
CLL
CLL
CLLEND-----------------------------------------------------------------
C*    Vn                                                                
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL   4.4   08/09/97   L_BL_LSPICE specifies mixed phase precipitation
CLL                    scheme.                   D.Wilson
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
C*L  Arguaments -------------------------------------------------------
      SUBROUTINE SF_RIB (
     & P_POINTS,LAND_PTS,LAND_MASK,INT_STOM,
     & GATHER,P1,LAND_INDEX,
     & NSICE,SICE_INDEX,ICE_FRACT,
     & PSTAR,AK_1,BK_1,Q_1,QW_1,QCL_1,QCF_1,
     & CF_1,T_1,TL_1,QSL,QSTAR,QSTAR_LEAD,
     & QS1,TSTAR_NL,Z1,Z0M_EFF,Z0M,Z0H,Z0HS,Z0MSEA,
     & WIND_PROFILE_FACTOR,U_1_P,U_0_P,V_1_P,V_0_P,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,
     & LYING_SNOW,GC,RESIST,RIB,RIB_LEAD,PSIS,VSHR,ALPHA1,
     & BT_1,BQ_1,BF_1,FRACA,RESFS,DQ,DQ_LEAD,DTEMP,
     & DTEMP_LEAD,L_BL_LSPICE,
     & LTIMER)
      IMPLICIT NONE


      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed
     &,LAND_PTS            ! IN Number of land points to be processed.
     &,NSICE               ! IN Number of sea-ice points.
     &,SICE_INDEX(P_POINTS) ! IN Index vector for gather to sea-ice
C                          !     points

     &,LAND_INDEX(LAND_PTS)! IN Index for compressed land point array;
C                          !    i'th element holds position in the FULL
C                          !    field of the ith land pt to be
C                          !    processed
     &,P1                  ! IN First P-point to be processed.
      LOGICAL
     & GATHER              ! IN If true then leads variables are comp-
C                          !    ressed for sea-ice calculations. This
C                          !    saves duplicating calculations if there
C                          !    are a relatively few of sea-ice points.
C                          !    Set to false for a limited area run
C                          !    with a high proportion of sea-ice.
      LOGICAL
     & INT_STOM            ! IN T for interactive stomatal resistance.
     &,L_BL_LSPICE             ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
C
      REAL
     & AK_1                ! IN Hybrid "A" for lowest model layer.
     &,BK_1                ! IN Hybrid "B" for lowest model layer.
     &,CANOPY(LAND_PTS)    ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_PTS)     ! IN Surface capacity (max. surface water)
C                          !    (kg per sq metre).  F6416.
     &,CF_1(P_POINTS)      ! IN Cloud fraction for lowest atmospheric
C                          !    layer (decimal fraction).
     &,ICE_FRACT(P_POINTS) ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_POINTS)! IN Lying snow amount (kg per sq metre).
     &,PSTAR(P_POINTS)     ! IN Surface pressure (Pascals).
     &,Q_1(P_POINTS)       ! IN Specific humidity for lowest
C                          !    atmospheric layer (kg water per kg air)
     &,QCF_1(P_POINTS)     ! IN Cloud ice for lowest atmospheric layer
C                          !    (kg water per kg air).
     &,QCL_1(P_POINTS)     ! IN Cloud liquid water for lowest atm layer
C                          !    (kg water per kg air).
     &,QS1(P_POINTS)       ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSL(P_POINTS)       ! IN Saturated sp humidity at liquid/ice
C                          !    temperature and pressure of lowest
C                          !    atmospheric level.
     &,QSTAR(P_POINTS)     ! IN Surface saturated sp humidity. Holds
C                          !    value over sea-ice where ICE_FRACT > 0.
C                          !    i.e. Leads contribution not included.
     &,QSTAR_LEAD(P_POINTS) ! IN QSTAR for sea-ice leads.
C                          !    Missing data indicator over non sea-ice.
     &,GC(LAND_PTS)        ! IN Interactive canopy conductance
C                          !    to evaporation (m/s)
     &,RESIST(LAND_PTS)    ! IN Fixed "stomatal" resistance
C                          !    to evaporation (s/m)
     &,ROOTD(LAND_PTS)     ! IN "Root depth" (metres).  F6412.
     &,SMC(LAND_PTS)       ! IN Soil moisture content (kg per sq m).
C                          !    F621.
     &,SMVCCL(LAND_PTS)    ! IN Critical volumetric SMC (cubic metres
C                          !    per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_PTS)    ! IN Volumetric wilting point (cubic m of
C                          !    water per cubic m of soil).  F6231.
C
C    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
C          moisture content.                            ~~~~~~~~~
C
     &,T_1(P_POINTS)      ! IN Temperature for lowest atmospheric layer
C                         !    (Kelvin).
     &,TL_1(P_POINTS)     ! IN Liquid/frozen water temperature for
C                         !    lowest atmospheric layer (K).
     &,TSTAR_NL(P_POINTS) ! IN TSTAR No Leads: surface temperature
C                         !    over sea-ice fraction of gridsquare.
C                         !    =TSTAR over non sea-ice points.
     &,U_1_P(P_POINTS)    ! IN West-to-east wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
     &,V_1_P(P_POINTS)    ! IN South-to-north wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
     &,U_0_P(P_POINTS)    ! IN West-to-east component of ocean surface
C                         !    current (m/s; zero over land if U_0 OK).
C                         !    P grid.  F615.
     &,V_0_P(P_POINTS)    ! IN South-to-north component of ocean surface
C                         !    current (m/s; zero over land if V_0 OK).
C                         !    P grid.  F616.
     &,V_SOIL(LAND_PTS)   ! IN Volumetric soil moisture concentration
C                         !    in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_PTS)    ! IN Vegetated fraction.
     &,WIND_PROFILE_FACTOR(P_POINTS)
C                         ! IN For transforming effective surface transf
C                         !    coefficients to those excluding form drag
     +,Z0H(P_POINTS)      ! IN Roughness length for heat and moisture m
     +,Z0HS(P_POINTS)     ! IN Roughness length for heat and moisture
C                         !    transport over sea.
     &,Z0M(P_POINTS)      ! IN Roughness length for heat and moisture m
     &,Z0MSEA(P_POINTS)   ! IN Sea-surface roughness length for
C                         !    momentum (m).  F617.
     &,Z0M_EFF(P_POINTS)  ! IN Effective roughness length for momentum
     &,Z1(P_POINTS)       ! IN Height of lowest atmospheric level (m).
      LOGICAL
     & LAND_MASK(P_POINTS) ! IN .TRUE. for land; .FALSE. elsewhere. F60.
C
C  Output variables.
C
      REAL
     & ALPHA1(P_POINTS)   ! OUT Gradient of saturated specific humidity
C                         !     with respect to temperature between the
C                         !     bottom model layer and the surface
     &,BQ_1(P_POINTS)     ! OUT A buoyancy parameter for lowest atm
C                         !     level. ("beta-q twiddle").
     &,BF_1(P_POINTS)     
!        OUT A buoyancy parameter for lowest atm level. 
!            ("beta-q twiddle").
     &,BT_1(P_POINTS)     ! OUT A buoyancy parameter for lowest atm
C                         !     level. ("beta-T twiddle").
     &,DQ(P_POINTS)       ! OUT Sp humidity difference between surface
C                         !     and lowest atmospheric level (Q1 - Q*).
C                         !     Holds value over sea-ice where
C                         !     ICE_FRACT>0 i.e. Leads contribution not
C                         !     included.
     &,DQ_LEAD(P_POINTS)  ! OUT DQ for leads fraction of gridsquare.
     &,DTEMP(P_POINTS)    ! OUT Liquid/ice static energy difference
C                         !     between surface and lowest atmospheric
C                         !     level, divided by CP (a modified
C                         !     temperature difference).
C                         !     Holds value over sea-ice where
C                         !     ICE_FRACT>0 i.e. Leads contribution not
C                         !     included.
     &,DTEMP_LEAD(P_POINTS) ! OUT DTEMP for leads fraction of
C                         !       gridsquare.
     &,FRACA(P_POINTS)    ! OUT Fraction of surface moisture flux with
C                         !     only aerodynamic resistance.
     &,PSIS(P_POINTS)     !     Soil moisture availability factor.
     &,QW_1(P_POINTS)     ! OUT Total water content of lowest
C                         !     atmospheric layer (kg per kg air).
     &,RESFS(P_POINTS)    ! OUT Combined soil, stomatal and aerodynamic
C                         !     resistance factor = PSIS/(1+RS/RA) for
C                         !     fraction (1-FRACA)
     &,RESFT(P_POINTS)    ! OUT Total resistance factor
C                         !     FRACA+(1-FRACA)*RESFS.
     &,RIB(P_POINTS)      ! OUT Bulk Richardson number for lowest layer
     &,RIB_LEAD(P_POINTS) ! OUT Bulk Richardson no. for sea-ice leads
C                         !     at lowest layer. At non sea-ice points
C                         !     holds RIB for FCDCH calculation, then
C                         !      set to missing data indicator.
     &,VSHR(P_POINTS)     ! OUT Magnitude of surface-to-lowest-lev. wind
C
      LOGICAL LTIMER      ! Logical switch for TIMER diags
C*
C*L  Symbolic constants -----------------------------------------------
C

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

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C
C   (3) Derived local parameters.
C
      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     & ETAR=1./(1.-EPSILON) ! Used in calc of buoyancy parameter BETAC
     &,GRCP=G/CP           ! Used in calc of dT across surface layer.
     &,LCRCP=LC/CP         ! Evaporation-to-dT conversion factor.
     &,LFRCP=LF/CP         ! Freezing-to-dT conversion factor.
     &,LS=LF+LC            ! Latent heat of sublimation.
     &,LSRCP=LS/CP         ! Sublimation-to-dT conversion factor.
     & )

C   Define local storage.
C
C   (a) Workspace.
C
C*L  Workspace --------------------------------------------------------
      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,L           ! Loop counter (land field index).
     &,SI          ! Loop counter (sea-ice field index).
      REAL
     & AL          ! Temporary in calculation of buoyancy parameters.
     &,ALPHAL      ! Temporary in calculation of buoyancy parameters.
     &,BETAC       ! Temporary in calculation of buoyancy parameters.
     &,BETAQ       ! Temporary in calculation of buoyancy parameters.
     &,BETAT       ! Temporary in calculation of buoyancy parameters.
     &,D_T         ! Temporary in calculation of alpha1
     &,USHEAR      ! U-component of surface-to-lowest-level wind shear.
     &,VSHEAR      ! V-component of surface-to-lowest-level wind shear.
     &,VSHR2       ! Square of magnitude of surface-to-lowest-level
C                  ! wind shear.
     &,ZETAM       ! Temporary in calculation of CHN.
     &,ZETAH       ! Temporary in calculation of CHN.

      REAL
     & CHN(P_POINTS)     ! Neutral-stability value of CH, used as a firs
C                        ! approximation to the "true" CH.
     &,EPDT(P_POINTS)    ! "Potential" Evaporation * Timestep - dummy
C                        ! variale = 0.
     &,F_SE(P_POINTS)    ! Dummy variable - actual value set in
C                        ! 2nd Call to SF_RESIST in SF_EXCH


      EXTERNAL SF_RESIST
      IF (LTIMER) THEN
        CALL TIMER('SFRIB   ',3)
      ENDIF

C----------------------------------------------------------------------
CL  1 Calculate buoyancy parameters for the lowest model level.
C----------------------------------------------------------------------
      DO I=1,P_POINTS
        IF (L_BL_LSPICE) THEN
          QW_1(I) = Q_1(I) + QCL_1(I)                         ! P243.10
        ELSE
          QW_1(I) = Q_1(I) + QCL_1(I) + QCF_1(I)              ! P243.10
        ENDIF
        BETAT = 1.0 / T_1(I)                         ! P243.19 (1st eqn)
        BETAQ = C_VIRTUAL /
     &     ( 1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I) )
C                                                  ... P243.19 (2nd eqn)
C
       IF (TL_1(I).GT.TM.OR.L_BL_LSPICE) THEN
          ALPHAL = (EPSILON * LC * QSL(I)) / (R * TL_1(I) * TL_1(I))
C                                       ... P243.20 (Clausius-Clapeyron)
C
          AL = 1.0 / (1.0 + LCRCP*ALPHAL)                      ! P243.21
          BETAC = CF_1(I) * AL * (LCRCP*BETAT - ETAR*BETAQ)
C                                                  ... P243.19 (3rd eqn)
C
        ELSE
          ALPHAL = (EPSILON * LS * QSL(I)) / (R * TL_1(I) * TL_1(I))
C                                       ... P243.20 (Clausius-Clapeyron)
C
          AL = 1.0 / (1.0 + LSRCP*ALPHAL)                      ! P243.21
          BETAC = CF_1(I) * AL * (LSRCP*BETAT - ETAR*BETAQ)
C                                                  ... P243.19 (3rd eqn)
C
        ENDIF
        BT_1(I) = BETAT - ALPHAL*BETAC               ! P243.18 (1st eqn)
        BQ_1(I) = BETAQ + BETAC                      ! P243.18 (2nd eqn)
        BF_1(I) = BETAQ *EPSILON*ETAR               ! P243.18 (2nd eqn)
C
C
C***********************************************************************
C   2. Calculate ALPHA1 {qsat(T*,P*) -qsat(TL1,P*)}/
C      {T*-TL1} for land and sea-ice points only -
C      set to zero over sea points
C***********************************************************************
C
        D_T = TSTAR_NL(I)-TL_1(I)
        IF (D_T .GT. 0.05 .OR. D_T .LT. -0.05) THEN
          ALPHA1(I) = (QSTAR(I) - QS1(I)) / D_T
        ELSE IF (TL_1(I).GT.TM) THEN
          ALPHA1(I) = (EPSILON*LC*QS1(I)*(1.0+C_VIRTUAL*QS1(I)))
     &              / (R*TL_1(I)*TL_1(I))
        ELSE
          ALPHA1(I) = (EPSILON*LS*QS1(I)*(1.0+C_VIRTUAL*QS1(I)))
     &              / (R*TL_1(I)*TL_1(I))
        ENDIF
        IF (.NOT.LAND_MASK(I).AND.ICE_FRACT(I).LE.0.0)
     &  ALPHA1(I) = 0.0   ! Sea points
      ENDDO
C-----------------------------------------------------------------------
CL  3 Calculate temperature (strictly, liquid/ice static energy) and
CL    humidity jumps, and wind shear, across the surface layer.
CL    Separate values are required for the leads and ice fractions
CL    of sea-ice grid-squares.
C-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=1,P_POINTS
          DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1(I) + Z0M_EFF(I) - Z0H(I) )
C                                                             ! P243.118
          DQ(I) = QW_1(I) - QSTAR(I)                          ! P243.119
C
          DTEMP_LEAD(I) = 1.E30                  ! Missing data indicato
          DQ_LEAD(I) = 1.E30                     ! Missing data indicato
          USHEAR = U_1_P(I) - U_0_P(I)
          VSHEAR = V_1_P(I) - V_0_P(I)
          VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
          VSHR(I) = SQRT(VSHR2)
C                                ... P243.120 (previous 4 lines of code)
        ENDDO
C-----------------------------------------------------------------------
CL  4 Calculate leads values by looping round sea-ice points only.
C     Avoids an if test in the above loop, so code can run faster.
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI=1,NSICE
          I = SICE_INDEX(SI)
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*(Z1(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(SI)
        ENDDO
      ELSE
      DO I=1,P_POINTS
        USHEAR = U_1_P(I) - U_0_P(I)
        VSHEAR = V_1_P(I) - V_0_P(I)

        VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
        VSHR(I) = SQRT(VSHR2)
C                                ... P243.120 (previous 4 lines of code)
        DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1(I) + Z0M_EFF(I) - Z0H(I) )
C                                                             ! P243.118
        DQ(I) = QW_1(I) - QSTAR(I)                            ! P243.119
C
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*(Z1(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(I)
        ELSE
          DTEMP_LEAD(I) = RMDI            ! Missing data indicator
          DQ_LEAD(I) = RMDI               ! Missing data indicator
        ENDIF
C
      ENDDO
      ENDIF                   ! End of IF (GATHER) THEN... ELSE...


C-----------------------------------------------------------------------
CL  6 Evaporation over land surfaces without snow is limited by
CL    soil moisture availability and stomatal resistance.
C     Set FRACA (= fA in the documentation) according to P243.68,
C     PSIS according to P243.65, and RESFS (= fS) according to P243.75
C     and P243.61, using neutral-stability value of CH (as explained
C     in section (v) of the P243 documentation).
C-----------------------------------------------------------------------

      DO I=1,P_POINTS

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  6.2.1 Calculate neutral stability value of CH (CHN), as an
CL        approximation to CH.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ZETAM = LOG ( (Z1(I) + Z0M_EFF(I))/Z0M_EFF(I) )
        ZETAH = LOG ( (Z1(I) + Z0M_EFF(I))/Z0H(I) )
        CHN(I) = (VKMAN/ZETAH) * (VKMAN/ZETAM) * WIND_PROFILE_FACTOR(I)

        EPDT(I) = 0.0  !Dummy variable for SF_RESIST

      ENDDO         ! Calc of CHN


      CALL SF_RESIST (
     & P_POINTS,LAND_PTS,LAND_MASK,INT_STOM,
     & P1,LAND_INDEX,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,DQ,
     & EPDT,LYING_SNOW,GC,RESIST,VSHR,CHN,PSIS,FRACA,
     & RESFS,F_SE,RESFT,LTIMER
     & )


C-----------------------------------------------------------------------
CL  7 Calculate bulk Richardson numbers for the surface layer.
CL    At sea-ice points RIB contains value for ice only (not leads).
CL      Initialise RIB_LEAD to RIB so that it contains sensible
CL      values at non sea ice points for the FCDCH calculation below.
C-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=1,P_POINTS
          RIB(I) = G * Z1(I) *
     &                 ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) ) /
     &                 ( VSHR(I)*VSHR(I) )
          RIB_LEAD(I) = RIB(I)
        ENDDO
C-----------------------------------------------------------------------
CL  7.1  Calculate bulk Richardson no. for leads at sea-ice points
CL       only.
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI = 1,NSICE
          I = SICE_INDEX(SI)
          RIB_LEAD(I) = G * Z1(I) *
     &                      ( BT_1(I) * DTEMP_LEAD(I) +
     &                        BQ_1(I) * RESFT(I) * DQ_LEAD(I) ) /
     &                      ( VSHR(I) * VSHR(I) )
C                            ... P2430.2, for sea-ice leads.
        ENDDO
      ELSE
      DO I=1,P_POINTS
        RIB(I) = G * Z1(I) *
     &               ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) ) /
     &               ( VSHR(I)*VSHR(I) )
C                            ... P243.43 (G times middle line is surface
C                                layer buoyancy difference, P243.25)
C
        RIB_LEAD(I) = RIB(I)
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          RIB_LEAD(I) = G * Z1(I) *
     &                      ( BT_1(I) * DTEMP_LEAD(I) +
     &                        BQ_1(I) * RESFT(I) * DQ_LEAD(I) ) /
     &                      ( VSHR(I) * VSHR(I) )
C                            ... P2430.2, for sea-ice leads.
        ENDIF
      ENDDO
      ENDIF           ! End of IF (GATHER) THEN... ELSE...

C
      IF (LTIMER) THEN
        CALL TIMER('SFRIB   ',4)
      ENDIF
      RETURN
      END

