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
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)

      SUBROUTINE SF_FLUX (
     & P_POINTS,P_FIELD,LAND_PTS,LAND_FIELD,LAND_MASK,L_LAND,P1,LAND1,
     & LAND_INDEX,
     & ALPHA1,DQ,DQ_LEAD,DTEMP,DTEMP_LEAD,DZ1,HCONS,ICE_FRACT,
     & LYING_SNOW,QS1,QW_1,RADNET_C,RESFT,RHOKE,RHOKH_1,TI,TL_1,TS1,
     & Z0H,Z0M_EFF,Z1_TQ,Z1_UV,
     & ASHTF,E_SEA,EPOT,FQW_1,FTL_1,H_SEA,RHOKPM,RHOKPM_POT,LTIMER 
     &,TSTAR,VFRAC,TIMESTEP,CANCAP 
     & )

      IMPLICIT NONE

      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,P_FIELD             ! IN Total number of P-grid points.
     &,P1                  ! IN First P-point to be processed.
     &,LAND1               ! IN First land point to be processed.
     &,LAND_PTS            ! IN Number of land points to be processed.
     &,LAND_FIELD          ! IN Total number of land points.
     &,LAND_INDEX(LAND_FIELD) ! IN Index for compressed land point array
!                               ith element holds position in the FULL
!                               field of the ith land pt to be processed

      LOGICAL
     & L_LAND              ! IN switch to only calculate for land points
     &,LTIMER
      REAL
     & ALPHA1(P_FIELD)     ! IN Gradient of saturated specific humidity
!                                with respect to temperature between the
!                                bottom model layer and the surface
     &,DQ(P_FIELD)         ! IN Sp humidity difference between surface
!                               and lowest atmospheric level (Q1 - Q*).
!                               Holds value over sea-ice where ICE_FRACT
!                               >0 i.e. Leads contribution not included.
     &,DQ_LEAD(P_FIELD)    ! IN DQ for leads fraction of gridsquare.
!                               Missing data indicator over non sea-ice.
     &,DTEMP(P_FIELD)      ! IN Liquid/ice static energy difference
!                               between surface and lowest atmospheric
!                               level, divided by CP (a modified
!                               temperature difference).
!                               Holds value over sea-ice where ICE_FRACT
!                               >0 i.e. Leads contribution not included.
     &,DTEMP_LEAD(P_FIELD) ! IN DTEMP for leads fraction of gridsquare.
!                               Missing data indicator over non sea-ice.
     &,DZ1                 ! IN  Thickness of the top soil layer (m).
     &,HCONS(LAND_FIELD)   ! IN Soil thermal conductivity including
!                               the effects of water and ice(W/m/K).
     &,ICE_FRACT(P_FIELD)  ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_FIELD) ! IN Lying snow amount (kg per sq metre).
     &,QS1(P_FIELD)        ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                               atmospheric layer (kg per kg air).
     &,RESFT(P_FIELD)      ! IN Total resistance factor
!                                FRACA+(1-FRACA)*RESFS.
     &,RHOKE(P_FIELD)      ! IN For FQW, then *GAMMA(1) for implicit
!                               calcs
     &,RHOKH_1(P_FIELD)    ! IN For FTL, then *GAMMA(1) for implicit
!                          !    calcs
     &,TI(P_FIELD)         ! IN Temperature of sea-ice surface layer (K)
     &,TL_1(P_FIELD)       ! IN Liquid/frozen water temperature for
!                               lowest atmospheric layer (K).
     &,TS1(LAND_FIELD)     ! IN Temperature of top soil layer (K)
     &,Z0H(P_FIELD)        ! IN Roughness length for heat and moisture m
     &,Z0M_EFF(P_FIELD)    ! IN Effective roughness length for momentum
     &,Z1_TQ(P_FIELD)      ! IN Height of lowest atmospheric level (m).
     &,Z1_UV(P_FIELD)      ! IN Height of lowest atmospheric level (m).

      LOGICAL
     & LAND_MASK(P_FIELD)  ! IN .TRUE. for land; .FALSE. elsewhere. F60.


! output
      REAL
     & ASHTF(P_FIELD)      ! OUT Coefficient to calculate surface
!                             heat flux into soil or sea-ice (W/m2/K).
     &,E_SEA(P_FIELD)      ! OUT Evaporation from sea times leads
!                             fraction (kg/m2/s). Zero over land.
     &,EPOT(P_FIELD)       ! OUT potential evaporation on P-grid
!                             (kg/m2/s).
     &,FQW_1(P_FIELD)      ! OUT "Explicit" surface flux of QW (i.e.
!                             evaporation), on P-grid (kg/m2/s).
     &,FTL_1(P_FIELD)      ! OUT "Explicit" surface flux of TL = H/CP.
!                            (sensible heat / CP).
     &,H_SEA(P_FIELD)      ! OUT Surface sensible heat flux over sea
!                            times leads fraction (W/m2).
!                            Zero over land.
     &,RHOKPM(P_FIELD)     ! OUT NB NOT * GAMMA for implicit calcs.
     &,RHOKPM_POT(P_FIELD) ! OUT Surface exchange coeff. for
!                            potential evaporation.

!-----------------------------------------------------------------------
! Extra variables required for the thermal canopy options.
!-----------------------------------------------------------------------
      REAL
     & TSTAR(P_FIELD)   ! IN Mean gridsquare surface temperature (K).
     &,VFRAC(LAND_FIELD)! IN Vegetation fraction.
     &,TIMESTEP         ! IN Timestep (s).
     &,CANCAP(P_FIELD)  ! INOUT Volumetric heat capacity of
C                       !       vegetation canopy (J/Kg/m3).
     &,RADNET_C(P_FIELD)! INOUT Adjusted net radiation for vegetation
C                       !       over land (W/m2).

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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
       REAL
     & SBCON                       ! Stefan-Boltzmann constant
!                                  ! (W/m**2/K**4).
      PARAMETER ( SBCON=5.67E-8 )
      INTEGER
     & CAN_MODEL            ! 1 for no thermal canopy (pre 4.5 UM).
C                           ! 2 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature.
C                           ! 3 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature, plus canopy thermal
C                           !   capacity.
     &,REX_MODEL            ! 1 for uniform root density profile
C                           !   (pre 4.5 UM) with MOSES I
C                           !   rootdepths.
C                           ! 2 for exponential root density
C                           !   profile with "old" rootdepths.
     &,TF_MODEL             ! 1 for uniformily distributed canopy
C                           !   water (pre 4.5 UM).
C                           ! 2 for bimodally distributed canopy
C                           !   water with random overlap.
C                           ! 3 for bimodally distributed canopy
C                           !   water with maximum overlap.
C
C For pre 4.5 UM MOSES I choose:
C     PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)
C
      PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)


!   (3) Derived local parameters.
      REAL GRCP,LS
      PARAMETER (
     & GRCP=G/CP           ! Used in calc of dT across surface layer.
     &,LS=LF+LC            ! Latent heat of sublimation.
     & )

!   (b) Scalars.

      INTEGER
     & I                   ! Loop counter (horizontal field index).
     &,L                   ! Loop counter (land point field index).

      REAL
     & DS_RATIO            ! 2 * snowdepth / depth of top soil layer.
     &,FQW_ICE             ! "Explicit" surface flux of QW for sea-ice
!                            fraction of gridsquare.
     &,FTL_ICE             ! "Explicit" surface flux of TL for sea-ice
!                            fraction of gridsquare.
     &,LAT_HEAT
     &,RAD_REDUC           ! Radiation term required for surface flux
!                            calcs.

!! Workspace
      REAL
     & DQ1(P_FIELD)        ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1


      EXTERNAL TIMER

!***********************************************************************
!  Calculate sensible and latent heat fluxes for land points.
!***********************************************************************
      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',3)
      ENDIF

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)

        ASHTF(I) = 2.0 * HCONS(L) / DZ1
        IF (LYING_SNOW(I) .GT. 0.0) THEN
          LAT_HEAT = LS
          DS_RATIO = 2.0 * LYING_SNOW(I) / (RHO_SNOW * DZ1)
          IF (DS_RATIO.LE.1.0) THEN
            ASHTF(I) =  ASHTF(I) /
     &                  (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
          ELSE IF (LYING_SNOW(I) .LT. 5.0E3) THEN
            ASHTF(I) =  ASHTF(I)*SNOW_HCON / HCONS(L)
          ENDIF
        ELSE
          LAT_HEAT = LC
        ENDIF
        E_SEA(I) = 0.0
        H_SEA(I) = 0.0

!-----------------------------------------------------------------------
! Options for treatment of vegetation thermal canopy
!-----------------------------------------------------------------------
        IF (CAN_MODEL .EQ. 1) THEN
          ASHTF(I) = ASHTF(I)
          CANCAP(I) = 0.0

        ELSEIF (CAN_MODEL .EQ. 2) THEN
          ASHTF(I) = (1.0-VFRAC(L))*ASHTF(I) +
     &                VFRAC(L)*4.0*SBCON*(TS1(L)**3)
          CANCAP(I) = 0.0

        ELSEIF (CAN_MODEL .EQ. 3) THEN
          ASHTF(I) = (1.0-VFRAC(L))*ASHTF(I) +
     &                VFRAC(L)*4.0*SBCON*(TS1(L)**3)
          CANCAP(I) = (1.0-VFRAC(L))*0.0 + VFRAC(L)*10.0*1.0E4

        ENDIF
        ASHTF(I) = ASHTF(I) + CANCAP(I)/TIMESTEP

        RHOKPM(I) = RHOKH_1(I) / ( RHOKH_1(I) *
     &              (LAT_HEAT*ALPHA1(I)*RESFT(I) + CP) + ASHTF(I) )
        RADNET_C(I)=RADNET_C(I) + CANCAP(I)*(TSTAR(I)-TS1(L))/TIMESTEP
        RAD_REDUC = RADNET_C(I) - ASHTF(I) *
     &        ( TL_1(I) - TS1(L) + GRCP * (Z1_TQ(I)
     &                                   + Z0M_EFF(I) - Z0H(I)) )
        DQ1(I) = (QS1(I)-QW_1(I)) +
     &            GRCP * ALPHA1(I) * (Z1_TQ(I) + Z0M_EFF(I) - Z0H(I))
        FQW_1(I) = RESFT(I)*RHOKPM(I)*( ALPHA1(I) *
     &             RAD_REDUC + (CP*RHOKH_1(I) + ASHTF(I))* DQ1(I) )
        RHOKPM_POT(I)=RHOKH_1(I) / ( RHOKH_1(I) *
     &              (LAT_HEAT*ALPHA1(I) + CP) + ASHTF(I) )
        EPOT(I) = RHOKPM_POT(I)*( ALPHA1(I) *
     &             RAD_REDUC + (CP*RHOKH_1(I) + ASHTF(I))* DQ1(I) )
        FTL_1(I) = RHOKPM(I) *
     &          ( RAD_REDUC - LAT_HEAT*RESFT(I)*RHOKH_1(I)*DQ1(I) )

!***********************************************************************
!  Calculate sensible and latent heat fluxes for sea and sea-ice points
!***********************************************************************

      ENDDO ! loop over land points

      IF (.NOT.L_LAND) THEN
        DO I=P1,P1+P_POINTS-1
          IF ( .NOT. LAND_MASK(I).AND.ICE_FRACT(I).GT.0.0) THEN !sea-ice
            ASHTF(I) = 2 * KAPPAI / DE
            E_SEA(I) = - (1. - ICE_FRACT(I))*RHOKH_1(I)*DQ_LEAD(I)

            H_SEA(I) = - (1. - ICE_FRACT(I))*RHOKH_1(I)*CP*DTEMP_LEAD(I)

!***********************************************************************
! Calculate the sensible and latent heat fluxes from sea-ice portion
! of gridbox. Weight RHOKPM by ICE_FRACT for use in IMPL_CAL.
!***********************************************************************

            RHOKPM(I) = RHOKH_1(I) / ( RHOKH_1(I) *
     &                               (LS * ALPHA1(I) + CP) + ASHTF(I) )
            RAD_REDUC = RADNET_C(I) - ICE_FRACT(I) * ASHTF(I) *
     &          ( TL_1(I) - TI(I) + GRCP * (Z1_TQ(I)
     &                                     + Z0M_EFF(I) - Z0H(I)) )
            DQ1(I) = (QS1(I)-QW_1(I)) +
     &              GRCP * ALPHA1(I) * (Z1_TQ(I) + Z0M_EFF(I) - Z0H(I))
            FQW_ICE = RHOKPM(I) * ( ALPHA1(I) * RAD_REDUC +
     &           (CP * RHOKH_1(I) + ASHTF(I)) * DQ1(I) * ICE_FRACT(I) )
            FTL_ICE = RHOKPM(I) * ( RAD_REDUC -
     &                     ICE_FRACT(I) * LS * RHOKH_1(I) * DQ1(I) )
            RHOKPM(I) = ICE_FRACT(I) * RHOKPM(I)
            RHOKPM_POT(I)=RHOKPM(I)

!***********************************************************************
! Calculate the total flux over the gridbox
!***********************************************************************

            FTL_1(I) = H_SEA(I)/CP + FTL_ICE
            FQW_1(I) = E_SEA(I) + FQW_ICE
            EPOT(I) = E_SEA(I) + FQW_ICE
!       Sea points
          ELSE IF( .NOT.LAND_MASK(I) .AND. .NOT.ICE_FRACT(I).GT.0.0 )
     &    THEN
            E_SEA(I) = - RHOKH_1(I) * DQ(I)
            H_SEA(I) = - RHOKH_1(I) * CP * DTEMP(I)
            FQW_1(I) = E_SEA(I)
            EPOT(I) = E_SEA(I)
            FTL_1(I) =  H_SEA(I) / CP
            RHOKPM(I) = 0.0
            RHOKPM_POT(I)=RHOKPM(I)
            ASHTF(I) = 1.0

        ENDIF        ! sea/sea-ice block
        ENDDO

      ENDIF ! L_LAND
      IF (LTIMER) THEN
        CALL TIMER('SF_FLUX ',4)
      ENDIF

      RETURN
      END
