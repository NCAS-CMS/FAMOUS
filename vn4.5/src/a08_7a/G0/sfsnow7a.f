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
C
CLL  SUBROUTINE SFSNOW ------------------------------------------------
CLL
CLL  Purpose:  Calculates the decrease/increase in snowdepth due to the
CLL            sublimation/deposition of lying snow; adds the large-
CLL            scale and convective snowfall to the snowdepth;
CLL            updates the snow layer temperature;
CLL            melts snow when the snow layer temperature is above
CLL            the melting point of ice.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  Programming standard: Unified Model Documentation Paper No.4
CLL                        version no. 2, dated 18/1/90.
CLL
CLL  Logical component covered: P251.
CLL
CLL  System task:
CLL
CLL  Documentation: um documentation paper no 25
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE SFSNOW(
     & NPNTS,SOIL_PTS,SOIL_INDEX,
     & CONV_SNOW,LS_SNOW,DZ_1,HCONS,SNOW_FRAC,SNOW_SUB,SNOW_SURF_HTF,
     & TSOIL_1,TSTAR,TIMESTEP,
     & LYING_SNOW,RGRAIN,L_SNOW_ALBEDO,SNOWMELT,TSNOW,
     & SNOMLT_SUB_HTF,SNOW_SOIL_HTF,STF_HF_SNOW_MELT,LTIMER)

      IMPLICIT NONE

      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.
     &,SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & CONV_SNOW(NPNTS)     ! IN Convective snowfall (kg/m2/s).
     &,LS_SNOW(NPNTS)       ! IN Large-scale snowfall (kg/m2/s).
     &,DZ_1                 ! IN Soil surface layer depth (m).
     &,HCONS(NPNTS)         ! IN Thermal conductivity of surface soil
!                           !    layer (W/m/K).
     &,SNOW_FRAC(NPNTS)     ! IN Snow-cover fraction.
     &,SNOW_SUB(NPNTS)      ! IN Sublimation of lying snow (kg/m2/s).
     &,SNOW_SURF_HTF(NPNTS) ! IN Snow surface heat flux (W/m2).
     &,TSOIL_1(NPNTS)       ! IN Soil surface layer temperature (K).
     &,TIMESTEP             ! IN Timestep (s).
     &,TSTAR(NPNTS)         ! IN Snow surface temperature (K).

      REAL
     & LYING_SNOW(NPNTS)    ! INOUT Snow on the ground (kg/m2).
     &,RGRAIN(NPNTS)        ! INOUT Snow grain size (microns).
     &,SNOWMELT(NPNTS)      ! IN    Surface snowmelt (kg/m2/s).
!                           ! OUT   Total snowmelt (kg/m2/s).
     &,TSNOW(NPNTS)         ! INOUT Snow surface layer temperature (K).

      REAL
     & SNOMLT_SUB_HTF(NPNTS)! OUT Sub-surface snowmelt heat flux (W/m2).
     &,SNOW_SOIL_HTF(NPNTS) ! OUT Heat flux from snow to soil (W/m2).

      LOGICAL
     & STF_HF_SNOW_MELT     ! IN Stash flag for snow melt heat flux.
     &,L_SNOW_ALBEDO        ! IN Flag for prognostic snow albedo.
     &,LTIMER               ! IN Logical for TIMER.

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

! Local variables
      REAL
     & ASNOW               ! Reciprocal areal heat capacity of surface
!                          ! snow layer (K m2 / J).
     &,R0                  ! Grain size for fresh snow (microns).
     &,RMAX                ! Maximum snow grain size (microns).
     &,RATE                ! Grain area growth rate (microns**2 / s).
     &,SNOWFALL            ! Snowfall in timestep (kg/m2).
     &,SNOWDEPTH           ! Local snowdepth (m).
     &,SNOMLT_SUB          ! Sub-surface snow melt (kg/m2).
      PARAMETER (R0 = 50., RMAX = 2000.)
      INTEGER I,J          ! Loop counters.

! NO EXTERNAL SUBROUTINES CALLED

      IF (LTIMER) THEN
        CALL TIMER('SFSNOW  ',103)
      ENDIF
      ASNOW = 1./(SNOW_HCAP*DEFF_SNOW)

!-----------------------------------------------------------------------
! Update TSNOW for land points without permanent ice cover.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        IF (SNOW_FRAC(I) .GT. 0.0) THEN
          SNOWDEPTH = LYING_SNOW(I)/(RHO_SNOW*SNOW_FRAC(I))
          SNOWDEPTH = MAX (SNOWDEPTH, DEFF_SNOW/2)
          SNOW_SOIL_HTF(I) = SNOW_FRAC(I)*(TSNOW(I) - TSOIL_1(I)) /
     &                       ( (2*SNOWDEPTH - DEFF_SNOW)/(2*SNOW_HCON)
     &                                           + DZ_1 / (2*HCONS(I)) )
          TSNOW(I) = TSNOW(I) + TIMESTEP*(ASNOW/SNOW_FRAC(I))* 
     &                             (SNOW_SURF_HTF(I) - SNOW_SOIL_HTF(I))
        ELSE
          SNOW_SOIL_HTF(I) = 0.
          TSNOW(I) = TSOIL_1(I)
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Increment snowdepth by sublimation and surface melt.
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        LYING_SNOW(I) = LYING_SNOW(I) - TIMESTEP *
     &                                     ( SNOW_SUB(I) + SNOWMELT(I) )
        LYING_SNOW(I) = MAX( LYING_SNOW(I), 0. )
      ENDDO

!-----------------------------------------------------------------------
! Melt snow over land if TSNOW is above freezing.
! Increment snowdepth by subsurface melt.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        SNOMLT_SUB = 0.0
        IF (TSNOW(I).GT.TM .AND. LYING_SNOW(I).GT.0.0) THEN
          SNOMLT_SUB = MIN( LYING_SNOW(I),
     &                         SNOW_FRAC(I)*(TSNOW(I) - TM)/(LF*ASNOW) )
          TSNOW(I) = TM
          LYING_SNOW(I) = LYING_SNOW(I) - SNOMLT_SUB
          SNOWMELT(I) = SNOWMELT(I) + SNOMLT_SUB/TIMESTEP
        ENDIF
        IF (STF_HF_SNOW_MELT) SNOMLT_SUB_HTF(I) = LF*SNOMLT_SUB
      ENDDO

!-----------------------------------------------------------------------
! Increment snowdepth by snowfall.
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        LYING_SNOW(I) = LYING_SNOW(I) + TIMESTEP *
     &                                     ( LS_SNOW(I) + CONV_SNOW(I) )
      ENDDO

!-----------------------------------------------------------------------
! Increment snow grain size used in albedo calculations
!-----------------------------------------------------------------------
      IF ( L_SNOW_ALBEDO ) THEN
        DO I=1,NPNTS
          IF ( LYING_SNOW(I) .GT. 0.) THEN
            SNOWFALL = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I))
            RATE = 0.6
            IF (TSTAR(I) .LT. TM) THEN
              IF (RGRAIN(I) .LT. 150.) THEN
                RATE = 0.06
              ELSE
                RATE = 0.23E6*EXP(-3.7E4/(8.13451*TSTAR(I)))
              ENDIF
            ENDIF
            RGRAIN(I) = SQRT( RGRAIN(I)**2 + (RATE/3.14159)*TIMESTEP )
     &                                   - (RGRAIN(I) - R0)*SNOWFALL/2.5
            RGRAIN(I) = MIN( RMAX, RGRAIN(I) )
            RGRAIN(I) = MAX( R0, RGRAIN(I) )
          ELSE
            RGRAIN(I) = R0
          ENDIF
        ENDDO
      ENDIF

      IF (LTIMER) THEN
        CALL TIMER('SFSNOW  ',104)
      ENDIF

      RETURN
      END
