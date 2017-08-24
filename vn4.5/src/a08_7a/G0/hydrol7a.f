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
!    SUBROUTINE HYDROL-------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYDROL (
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,NTYPE,
     &                   SOIL_INDEX,TILE_PTS,TILE_INDEX,
     &                   NPNTS,NSHYD,B,CAN_CPY,CON_RAIN,CON_SNOW,
     &                   E_CANOPY,EXT,HCAP,HCON,LS_RAIN,
     &                   LS_SNOW,SATCON,SATHH,SNOWSUB,SOIL_SURF_HTF,
     &                   SNOW_SURF_HTF,TSTAR_SNOW,FRAC,SNOW_FRAC,
     &                   TIMESTEP,V_SAT,V_WILT,
     &                   CAN_WCNT,HF_SNOW_MELT,STF_HF_SNOW_MELT,
     &                   RGRAIN,L_SNOW_ALBEDO,SMCL,STHF,STHU,
     &                   SNOW_DEP,TSOIL,TSNOW,
     &                   CAN_WCNT_GB,INFIL,SMC,SNOW_MELT,
     &                   SNOMLT_SUB_HTF,STF_SUB_SURF_ROFF,
     &                   SUB_SURF_ROFF,SURF_ROFF,TOT_TFALL,LTIMER
     & )

      IMPLICIT NONE
!
! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme. Calls the following:
!
!     HEAT_CON - to calculate the thermal conductivity of the top
!                soil layer                           (Cox, 6/95)
!
!     SFSNOW - to calculate the sub-surface snowmelt
!              and update the lying snow amount    (Essery, 1/95)
!
!     INFILT - to calculate the maximum surface infiltration rate
!                                                     (Cox, 6/95)
!
!     SURF_HYD - to calculate canopy interception and
!                surface runoff         (Allen-Bett, Gregory, 90)
!
!     SOIL_HYD - to update the layer soil moisture contents
!                and calculate the drainage            (Cox 6/95)
!
!     SOIL_HTC - to update the soil layer temperatures and the
!                layer ice contents                    (Cox 6/95)
!
!     ICE_HTC - to update the layer temperatures for land ice
!                                                      (Cox 10/95)
!
!     SOIL_MC - to diagnose the soil moisture in the top metre
!                                                    (Essery 7/97)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1      6/96     New deck.  Peter Cox
!  4.4      7/97     MOSES II.  Richard Essery
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
      REAL
     + DZSOIL(4)               ! Soil layer thicknesses (m).
      DATA DZSOIL /0.100, 0.250, 0.650, 2.000 /
C-----------------------------------------------------------------------

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & LICE_PTS            ! IN Number of land ice points.
     &,NPNTS               ! IN Number of gridpoints.
     &,NSHYD               ! IN Number of soil moisture levels.
     &,SOIL_PTS            ! IN Number of soil points.
     &,NTYPE               ! IN Number of tiles.

      REAL
     & TIMESTEP            ! IN Model timestep (s).

      LOGICAL LTIMER       ! Logical switch for TIMER diags

      LOGICAL
     & STF_HF_SNOW_MELT    ! IN Stash flag for snowmelt heat flux.
     &,STF_SUB_SURF_ROFF   ! IN Stash flag for sub-surface runoff.
     &,L_SNOW_ALBEDO       ! IN Flag for prognostic snow albedo


!   Array arguments with intent(IN) :
      INTEGER
     & LICE_INDEX(NPNTS)   ! IN Array of land ice points.
     &,SOIL_INDEX(NPNTS)   ! IN Array of soil points.
     &,TILE_PTS(NTYPE)     ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTYPE)
!                          ! IN Index of tile points.

      REAL
     & B(NPNTS)            ! IN Clapp-Hornberger exponent.
     &,CAN_CPY(NPNTS,NTYPE-1)
!                          ! IN Canopy/surface capacity of
!                          !    snow-free land tiles (kg/m2).
     &,CON_RAIN(NPNTS)     ! IN Convective rain (kg/m2/s).
     &,CON_SNOW(NPNTS)     ! IN Convective snowfall (kg/m2/s).
     &,E_CANOPY(NPNTS,NTYPE-1)
!                          ! IN Canopy evaporation from snow-free
!                          !    land tiles (kg/m2/s).
     &,EXT(NPNTS,NSHYD)    ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
     &,HCAP(NPNTS)         ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS)         ! IN Soil thermal conductivity (W/m/K).
     &,LS_RAIN(NPNTS)      ! IN Large-scale rain (kg/m2/s).
     &,LS_SNOW(NPNTS)      ! IN Large-scale snowfall (kg/m2/s).
     &,SATCON(NPNTS)       ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
     &,SATHH(NPNTS)        ! IN Saturated soil water pressure (m).
     &,SNOWSUB(NPNTS)      ! IN Sublimation of lying snow (kg/m2/s).
     &,SOIL_SURF_HTF(NPNTS)! IN Net downward surface heat flux (W/m2)
!                          !    - snow-free land.
     &,SNOW_SURF_HTF(NPNTS)! IN Net downward surface heat flux (W/m2)
!                          !    - snow.
     &,TSTAR_SNOW(NPNTS)   ! IN Surface temperature (K).
     &,FRAC(NPNTS,NTYPE)   ! IN Tile fractions.
     &,SNOW_FRAC(NPNTS)    ! IN Fraction of gridbox with snow cover.
     &,V_SAT(NPNTS)        ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)       ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
!
!   Array arguments with intent(INOUT) :
!
      REAL
     & CAN_WCNT(NPNTS,NTYPE-1)
!                          ! INOUT Canopy water content for snow-free
!                          !       land tiles (kg/m2).
     &,HF_SNOW_MELT(NPNTS) ! INOUT Total snowmelt heat flux (W/m2).
     &,RGRAIN(NPNTS)       ! INOUT Snow grain size (microns).
     &,SMCL(NPNTS,NSHYD)   ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
     &,SNOW_DEP(NPNTS)     ! INOUT Snowmass (kg/m2).
     &,STHF(NPNTS,NSHYD)   ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)   ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,TSOIL(NPNTS,NSHYD)  ! INOUT Sub-surface temperatures (K).
     &,TSNOW(NPNTS)        ! INOUT Snow layer temperature (K).


!   Array arguments with intent(OUT) :
      REAL
     & CAN_WCNT_GB(NPNTS)   ! OUT Gridbox canopy water content (kg/m2).
     &,INFIL(NPNTS)         ! OUT Maximum surface infiltration
!                           !     rate (kg/m2/s).
     &,SMC(NPNTS)           ! OUT Soil moisture in the top metre (kg/m2)
     &,SNOW_MELT(NPNTS)     ! OUT Snowmelt (kg/m2/s).
     &,SNOMLT_SUB_HTF(NPNTS)! OUT Sub-surface snowmelt heat
!                           !     flux (W/m2).
     &,SUB_SURF_ROFF(NPNTS) ! OUT Sub-surface runoff (kg/m2/s).
     &,SURF_ROFF(NPNTS)     ! OUT Surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)     ! OUT Total throughfall (kg/m2/s).

! Local scalars:
      INTEGER
     & I,J                  ! WORK Loop counters.
     &,N                    ! WORK Tile loop counter.

! Local arrays:

      REAL
     & DSMC_DT(NPNTS)       ! WORK Rate of change of soil moisture
!                           !      due to water falling onto the
!                           !      surface after surface runoff
!                           !      (kg/m2/s).
     &,HCONS(NPNTS)         ! WORK Soil surface layer thermal
!                           !      conductivity including the effects
!                           !      of water and ice (W/m/K).
     &,INFIL_TILE(NPNTS,NTYPE)
!                           ! WORK Maximum surface infiltration
!                           !      rate for tiles (kg/m2/s).
     &,SNOW_SOIL_HTF(NPNTS) ! WORK Heat flux from snow to soil (W/m2).
     &,W_FLUX(NPNTS,0:NSHYD)! WORK Fluxes of water between layers
!                           !      (kg/m2/s).

! Tile parameters :
      REAL
     + INFIL_FAC(9)               ! Infiltration enhancement factor.
C----------------------------------------------------------------------
C                         BT   NT   C3G  C4G  Shr  Urb  Wat  Soil Ice
C----------------------------------------------------------------------
      DATA INFIL_FAC   /  4.0, 4.0, 2.0, 2.0, 2.0, 0.1, 1.0, 0.5, 0.0 /



! Function & Subroutine calls:
      EXTERNAL
     & HEAT_CON,SFSNOW,INFILT,SURF_HYD,SOIL_HYD,SOIL_HTC,ICE_HTC,SOILMC

! End of header--------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('HYDROL ',103)
      ENDIF

!----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer
!----------------------------------------------------------------------
      CALL HEAT_CON(NPNTS,HCON,STHU,STHF,V_SAT,HCONS,LTIMER)

!----------------------------------------------------------------------
! Calculate the subsurface snowmelt and update the snow mass
!----------------------------------------------------------------------
      CALL SFSNOW(NPNTS,SOIL_PTS,SOIL_INDEX,
     &            CON_SNOW,LS_SNOW,DZSOIL(1),HCONS,SNOW_FRAC,
     &            SNOWSUB,SNOW_SURF_HTF,TSOIL,TSTAR_SNOW,TIMESTEP,
     &            SNOW_DEP,RGRAIN,L_SNOW_ALBEDO,SNOW_MELT,TSNOW,
     &            SNOMLT_SUB_HTF,SNOW_SOIL_HTF,STF_HF_SNOW_MELT,LTIMER)

!----------------------------------------------------------------------
! Update the total snowmelt heat flux
!----------------------------------------------------------------------
      IF (STF_HF_SNOW_MELT) THEN
        DO I=1,NPNTS
          HF_SNOW_MELT(I)=LF*SNOW_MELT(I)
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Calculate the maximum surface infiltration rate
!----------------------------------------------------------------------
      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          INFIL_TILE(I,N) = INFIL_FAC(N)*SATCON(I)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------
      CALL SURF_HYD (NPNTS,NTYPE,TILE_PTS,TILE_INDEX,
     &               LICE_PTS,LICE_INDEX,
     &               CAN_CPY,E_CANOPY,FRAC,INFIL_TILE,CON_RAIN,LS_RAIN,
     &               SNOW_FRAC,SNOW_MELT,TIMESTEP,
     &               CAN_WCNT,CAN_WCNT_GB,DSMC_DT,SURF_ROFF,TOT_TFALL)

!-----------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------
      CALL SOIL_HYD (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,B,DZSOIL,
     &               EXT,DSMC_DT,SATCON,SATHH,TIMESTEP,V_SAT,
     &               SUB_SURF_ROFF,SMCL,STHU,W_FLUX,
     &               STF_SUB_SURF_ROFF,LTIMER)

!-----------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------
      IF (SOIL_PTS.NE.0) THEN
        CALL SOIL_HTC (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,B,
     &                 DZSOIL,HCAP,HCON,SATHH,SNOW_SOIL_HTF, 
     &                 SOIL_SURF_HTF,TIMESTEP,V_SAT,
     &                 W_FLUX,SMCL,STHU,STHF,TSOIL,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
! Update the sub-surface temperatures for land ice
!-----------------------------------------------------------------------
      IF (LICE_PTS.NE.0) THEN
        CALL ICE_HTC (NPNTS,NSHYD,LICE_PTS,LICE_INDEX,DZSOIL,
     &                SNOW_SURF_HTF,TIMESTEP,
     &                TSOIL,LTIMER)
! Copy surface layer temperature to TSNOW at ice points
        DO J=1,LICE_PTS
          I=LICE_INDEX(J)
          TSNOW(I)=TSOIL(I,1)
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Diagnose the soil moisture in the top metre.
!-----------------------------------------------------------------------
      CALL SOILMC ( NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,
     &              DZSOIL,STHU,V_SAT,V_WILT,SMC )

      IF (LTIMER) THEN
        CALL TIMER('HYDROL ',104)
      ENDIF

      RETURN
      END
