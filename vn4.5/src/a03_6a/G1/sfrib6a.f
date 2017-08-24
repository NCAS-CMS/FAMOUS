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
!!!  SUBROUTINE SF_RIB------------------------------------------------
!!!
!!!  Purpose: Calculate bulk Richardson number for surface layer
!!!
!!!          <- programmer of some or all of previous code changes
!!!  Simon Jackson, Roderick Smith
!!!
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!END-----------------------------------------------------------------

!!  Arguaments -------------------------------------------------------
      SUBROUTINE SF_RIB (
     & P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,LAND_MASK,L_LAND,INT_STOM,
     & P1,LAND1,
     & GATHER,LAND_INDEX,
     & NSICE,SICE_INDEX,ICE_FRACT,Q_1,QW_1,QCL_1,QCF_1,
     & T_1,TL_1,QSL,QSTAR,QSTAR_LEAD,
     & QS1,TSTAR_NL,Z1_TQ,Z1_UV,Z0M_EFF,Z0M,Z0H,Z0HS,Z0MSEA,
     & WIND_PROFILE_FACTOR,U_1_P,U_0_P,V_1_P,V_0_P,
     & ROOTD,SMVCCL,SMVCWT,SMC,VFRAC,V_SOIL,CANOPY,CATCH,
     & LYING_SNOW,GC,RESIST,DB,DB_LEAD,RIB,RIB_LEAD,PSIS,VSHR,ALPHA1,
     & BT_1,BQ_1,FRACA,RESFS,DQ,DQ_LEAD,DTEMP,DTEMP_LEAD,LTIMER
     & )

      IMPLICIT NONE


      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed
     &,LAND_PTS            ! IN Number of land points to be processed.
     &,P_FIELD             ! IN Total number of P-grid points
     &,LAND_FIELD          ! IN Total number of land points.
     &,NSICE               ! IN Number of sea-ice points.
     &,SICE_INDEX(P_FIELD) ! IN Index vector for gather to sea-ice
!                          !     points
     &,P1                  ! IN First P-point to be processed.
     &,LAND1               ! IN First land-point to be processed.
     &,LAND_INDEX(LAND_FIELD)
!                            IN Index for compressed land point array;
!                               i'th element holds position in the FULL
!                               field of the ith land pt to be
!                               processed
      LOGICAL
     & GATHER              ! IN If true then leads variables are comp-
!                               ressed for sea-ice calculations. This
!                               saves duplicating calculations if there
!                               are a relatively few of sea-ice points.
!                               Set to false for a limited area run
!                               with a high proportion of sea-ice.

      LOGICAL
     & INT_STOM            ! IN T for interactive stomatal resistance.
     &,LTIMER              ! IN logical for TIMER
     &,L_LAND              ! IN calculate only over land if .true.
     &,LAND_MASK(P_FIELD)  ! IN .TRUE. for land; .FALSE. elsewhere. F60.

      REAL
     & CANOPY(LAND_FIELD)  ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_FIELD)   ! IN Surface capacity (max. surface water)
!                               (kg per sq metre).  F6416.
     &,GC(LAND_FIELD)      ! IN Interactive canopy conductance
!                               to evaporation (m/s)
     &,ICE_FRACT(P_FIELD)  ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_FIELD) ! IN Lying snow amount (kg per sq metre).
     &,Q_1(P_FIELD)        ! IN Specific humidity for lowest
!                               atmospheric layer (kg water per kg air)
     &,QCF_1(P_FIELD)      ! IN Cloud ice for lowest atmospheric layer
!                               (kg water per kg air).
     &,QCL_1(P_FIELD)      ! IN Cloud liquid water for lowest atm layer
!                               (kg water per kg air).
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                               atmospheric layer (kg per kg air).
     &,QS1(P_FIELD)        ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSL(P_FIELD)        ! IN Saturated sp humidity at liquid/ice
!                               temperature and pressure of lowest
!                               atmospheric level.
     &,QSTAR(P_FIELD)      ! IN Surface saturated sp humidity. Holds
!                               value over sea-ice where ICE_FRACT > 0.
!                               i.e. Leads contribution not included.
     &,QSTAR_LEAD(P_FIELD) ! IN QSTAR for sea-ice leads.
!                               Missing data indicator over non sea-ice.
     &,RESIST(LAND_FIELD)  ! IN "Stomatal" resistance to evaporation
!                               (seconds per metre).  F6415.
     &,ROOTD(LAND_FIELD)   ! IN "Root depth" (metres).  F6412.
     &,SMC(LAND_FIELD)     ! IN Soil moisture content (kg per sq m).
!                               F621.
     &,SMVCCL(LAND_FIELD)  ! IN Critical volumetric SMC (cubic metres
!                               per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_FIELD)  ! IN Volumetric wilting point (cubic m of
!                               water per cubic m of soil).  F6231.

!    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
!          moisture content.                            ~~~~~~~~~

     &,T_1(P_FIELD)       ! IN Temperature for lowest atmospheric layer
!                              (Kelvin).
     &,TL_1(P_FIELD)      ! IN Liquid/frozen water temperature for
!                              lowest atmospheric layer (K).
     &,TSTAR_NL(P_FIELD)  ! IN TSTAR No Leads: surface temperature
!                              over sea-ice fraction of gridsquare.
!                              =TSTAR over non sea-ice points.
     &,U_1_P(P_FIELD)     ! IN West-to-east wind component for lowest
!                              atmospheric layer (m/s).  On P grid.
     &,V_1_P(P_FIELD)     ! IN South-to-north wind component for lowest
!                              atmospheric layer (m/s).  On P grid.
     &,U_0_P(P_FIELD)     ! IN West-to-east component of ocean surface
!                              current (m/s; zero over land if U_0 OK).
!                              P grid.  F615.
     &,V_0_P(P_FIELD)     ! IN South-to-north component of ocean surface
!                              current (m/s; zero over land if V_0 OK).
!                              P grid.  F616.
     &,V_SOIL(LAND_FIELD) ! IN Volumetric soil moisture concentration
!                              in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_FIELD)  ! IN Vegetated fraction.
     &,WIND_PROFILE_FACTOR(P_FIELD)
!                         ! IN For transforming effective surface
!                              transfer coefficients to those excluding
!                              form drag.
     &,Z0H(P_FIELD)       ! IN Roughness length for heat and moisture m
     &,Z0HS(P_FIELD)      ! IN Roughness length for heat and moisture
!                              transport over sea.
     &,Z0M(P_FIELD)       ! IN Roughness length for heat and moisture m
     &,Z0MSEA(P_FIELD)    ! IN Sea-surface roughness length for
!                              momentum (m).  F617.
     &,Z0M_EFF(P_FIELD)   ! IN Effective roughness length for momentum
     &,Z1_TQ(P_FIELD)     ! IN Height of lowest TQ level (m).
     &,Z1_UV(P_FIELD)     ! IN Height of lowest UV level (m).


!  Output variables.

      REAL
     & ALPHA1(P_FIELD)    ! OUT Gradient of saturated specific humidity
!                               with respect to temperature between the
!                               bottom model layer and the surface
     &,BQ_1(P_FIELD)      ! OUT A buoyancy parameter for lowest atm
!                               level. ("beta-q twiddle").
     &,BT_1(P_FIELD)      ! OUT A buoyancy parameter for lowest atm
!                               level. ("beta-T twiddle").
     &,DQ(P_FIELD)        ! OUT Sp humidity difference between surface
!                               and lowest atmospheric level (Q1 - Q*).
!                               Holds value over sea-ice where
!                               ICE_FRACT>0 i.e. Leads contribution not
!                               included.
     &,DQ_LEAD(P_FIELD)   ! OUT DQ for leads fraction of gridsquare.
     &,DTEMP(P_FIELD)     ! OUT Liquid/ice static energy difference
!                               between surface and lowest atmospheric
!                               level, divided by CP (a modified
!                               temperature difference).
!                               Holds value over sea-ice where
!                               ICE_FRACT>0 i.e. Leads contribution not
!                               included.
     &,DTEMP_LEAD(P_FIELD)! OUT DTEMP for leads fraction of gridsquare.
     &,FRACA(P_FIELD)     ! OUT Fraction of surface moisture flux with
!                               only aerodynamic resistance.
     &,PSIS(P_FIELD)      ! OUT Soil moisture availability factor.
     &,RESFS(P_FIELD)     ! OUT Combined soil, stomatal and aerodynamic
!                               resistance factor = PSIS/(1+RS/RA) for
!                               fraction (1-FRACA)
     &,RESFT(P_FIELD)     ! OUT Total resistance factor
!                               FRACA+(1-FRACA)*RESFS.
     &,DB(P_POINTS)       ! OUT Buoyancy difference between surface
!                         !     and lowest atmospheric level.
!                         !     Holds value over sea-ice where ICE_FRACT
!                         !     >0 i.e. Leads contribution not included.
     &,DB_LEAD(P_POINTS)  ! OUT DB for leads fraction of gridsquare.
!                         !     Missing data indicator over non sea-ice.
     &,RIB(P_FIELD)       ! OUT Bulk Richardson number for lowest layer
     &,RIB_LEAD(P_FIELD)  ! OUT Bulk Richardson no. for sea-ice leads
!                               at lowest layer. At non sea-ice points
!                               holds RIB for FCDCH calculation, then
!                               set to missing data indicator.
     &,VSHR(P_FIELD)      ! OUT Magnitude of surface-to-lowest-lev. wind


!  Symbolic constants -----------------------------------------------


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

!   (3) Derived local parameters.

      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     & ETAR=1./(1.-EPSILON)! Used in calc of buoyancy parameter BETAC
     &,GRCP=G/CP           ! Used in calc of dT across surface layer.
     &,LCRCP=LC/CP         ! Evaporation-to-dT conversion factor.
     &,LFRCP=LF/CP         ! Freezing-to-dT conversion factor.
     &,LS=LF+LC            ! Latent heat of sublimation.
     &,LSRCP=LS/CP         ! Sublimation-to-dT conversion factor.
     & )

!   Define local storage.

!   (a) Workspace.

!  Workspace --------------------------------------------------------
      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,L           ! Loop counter (land field index).
     &,SI          ! Loop counter (sea-ice field index).
      REAL
     & VIRT_FACTOR ! Factor for converting temperature to virtual temp.
     &,D_T         ! Temporary in calculation of alpha1.
     &,USHEAR      ! U-component of surface-to-lowest-level wind shear.
     &,VSHEAR      ! V-component of surface-to-lowest-level wind shear.
     &,VSHR2       ! Square of magnitude of surface-to-lowest-level
!                  ! wind shear.
     &,ZETAM       ! Temporary in calculation of CHN.
     &,ZETAH       ! Temporary in calculation of CHN.

      REAL
     & CHN(P_FIELD)      ! Neutral-stability value of CH, used as a
!                          first approximation to the "true" CH.
     &,EPDT(P_FIELD)     ! "Potential" Evaporation * Timestep - dummy
!                          variable = 0.
     &,F_SE(P_FIELD)     ! Dummy variable - actual value set in
!                          2nd Call to SF_RESIST in SF_EXCH


      EXTERNAL SF_RESIST

!----------------------------------------------------------------------
!!  1 Calculate buoyancy parameters for the lowest model level.
!----------------------------------------------------------------------
      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',3)
      ENDIF

      DO I=P1,P1+P_POINTS-1

          VIRT_FACTOR =  1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I)
          BT_1(I) = 1.0 / T_1(I)                     ! P243.19 (1st eqn)
          BQ_1(I) = C_VIRTUAL / VIRT_FACTOR
!                                                    ! P243.19 (2nd eqn)
!***********************************************************************
!!   2. Calculate gradient of saturated specific
!!      humidity for use in calculation of surface fluxes
!***********************************************************************
!
          D_T = TSTAR_NL(I)-TL_1(I)

          IF (D_T .GT. 0.05 .OR. D_T .LT. -0.05) THEN
            ALPHA1(I) = (QSTAR(I) - QS1(I)) / D_T
           ELSEIF (TL_1(I) .GT. TM) THEN

             ALPHA1(I)=(EPSILON * LC * QS1(I) *
     &                 (1.0+C_VIRTUAL*QS1(I)) )/(R*TL_1(I)*TL_1(I))

           ELSE

             ALPHA1(I)=(EPSILON * LS * QS1(I) *
     &                 (1.0+C_VIRTUAL*QS1(I)) )/(R*TL_1(I)*TL_1(I))
          ENDIF

! ALPHA1 set to zero over open sea, so that P-M only applies over land

          IF(.NOT.LAND_MASK(I).AND.ICE_FRACT(I).LE.0) ALPHA1(I)=0.0

      ENDDO

!-----------------------------------------------------------------------
!!  3 Calculate temperature (strictly, liquid/ice static energy) and
!!    humidity jumps, and wind shear, across the surface layer.
!!    Separate values are required for the leads and ice fractions
!!    of sea-ice grid-squares.
!-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=P1,P1+P_POINTS-1
          DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1_TQ(I) + Z0M_EFF(I) - Z0H(I) )
!                                                    ! P243.118
          DQ(I) = QW_1(I) - QSTAR(I)                 ! P243.119

          DTEMP_LEAD(I) = 1.E30                 ! Missing data indicator
          DQ_LEAD(I) = 1.E30                    ! Missing data indicator
          USHEAR = U_1_P(I) - U_0_P(I)
          VSHEAR = V_1_P(I) - V_0_P(I)
          VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
          VSHR(I) = SQRT(VSHR2)
!                                 ! P243.120 (previous 4 lines of code)
        ENDDO

!-----------------------------------------------------------------------
!!  4 Calculate leads values by looping round sea-ice points only.
!!    Avoids an if test in the above loop, so code can run faster.
!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI=1,NSICE
          I = SICE_INDEX(SI)
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*
     &                                 (Z1_TQ(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(SI)
        ENDDO
      ELSE
      DO I=P1,P1+P_POINTS-1
        USHEAR = U_1_P(I) - U_0_P(I)
        VSHEAR = V_1_P(I) - V_0_P(I)

        VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
        VSHR(I) = SQRT(VSHR2)
!                                ... P243.120 (previous 4 lines of code)
        DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1_TQ(I) + Z0M_EFF(I) - Z0H(I) )
!                                                             ! P243.118
        DQ(I) = QW_1(I) - QSTAR(I)                            ! P243.119
!
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*
     &                                  (Z1_TQ(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(I)
        ELSE
          DTEMP_LEAD(I) = RMDI            ! Missing data indicator
          DQ_LEAD(I) = RMDI               ! Missing data indicator
        ENDIF

      ENDDO
      ENDIF                   ! End of IF (GATHER) THEN... ELSE...

!-----------------------------------------------------------------------
!!  5 Evaporation over land surfaces without snow is limited by
!!    soil moisture availability and stomatal resistance.
!!    Set FRACA (= fA in the documentation) according to P243.68,
!!    PSIS according to P243.65, and RESFS (= fS) according to P243.75
!!    and P243.61, using neutral-stability value of CH (as explained
!!    in section (v) of the P243 documentation).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        ZETAM = LOG ( (Z1_UV(I) + Z0M_EFF(I))/Z0M_EFF(I) )
        ZETAH = LOG ( (Z1_TQ(I) + Z0M_EFF(I))/Z0H(I) )
        CHN(I) = (VKMAN/ZETAH) * (VKMAN/ZETAM) * WIND_PROFILE_FACTOR(I)
        CHN(I) = CHN(I)*VSHR(I)

        EPDT(I) = 0.0  !Dummy variable for SF_RESIST

      ENDDO         ! Calc of CHN


      CALL SF_RESIST (
     & P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,LAND_MASK,INT_STOM,P1,LAND1,
     & LAND_INDEX,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,DQ,
     & EPDT,LYING_SNOW,GC,RESIST,CHN,PSIS,FRACA,
     & RESFS,F_SE,RESFT,LTIMER)


!-----------------------------------------------------------------------
!!  6 Calculate bulk Richardson numbers for the surface layer.
!!    At sea-ice points RIB contains value for ice only (not leads).
!!      Initialise RIB_LEAD to RIB so that it contains sensible
!!      values at non sea ice points for the FCDCH calculation below.
!-----------------------------------------------------------------------
      IF (GATHER) THEN
       DO I=P1,P1+P_POINTS-1

            DB(I) = G * ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) )
            RIB(I) = Z1_UV(I) * DB(I) / ( VSHR(I)*VSHR(I) )
            DB_LEAD(I) = DB(I)
            RIB_LEAD(I) = RIB(I)

        ENDDO
!-----------------------------------------------------------------------
!!  6.1  Calculate bulk Richardson no. for leads at sea-ice points
!!       only.
!-----------------------------------------------------------------------

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI = 1,NSICE
          I = SICE_INDEX(SI)
          DB_LEAD(I) = G * ( BT_1(I) * DTEMP_LEAD(I) +
     &                       BQ_1(I) * RESFT(I) * DQ_LEAD(I) )
          RIB_LEAD(I) = Z1_UV(I) * DB_LEAD(I) / ( VSHR(I) * VSHR(I) )
!                            ... P2430.2, for sea-ice leads.
        ENDDO
      ELSE
      DO I=P1,P1+P_POINTS-1

        DB(I) = G * ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) )
        RIB(I) = Z1_UV(I) * DB(I) / ( VSHR(I)*VSHR(I) )
!
        DB_LEAD(I) = DB(I)
        RIB_LEAD(I) = RIB(I)
!
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          DB_LEAD(I) = G * ( BT_1(I) * DTEMP_LEAD(I) +
     &                       BQ_1(I) * RESFT(I) * DQ_LEAD(I) )
          RIB_LEAD(I) = Z1_UV(I) * DB_LEAD(I) / ( VSHR(I) * VSHR(I) )
!                            ... P2430.2, for sea-ice leads.

        ENDIF ! Ice_fract
      ENDDO
      ENDIF           ! End of IF (GATHER) THEN... ELSE...

      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',4)
      ENDIF
      RETURN
      END

