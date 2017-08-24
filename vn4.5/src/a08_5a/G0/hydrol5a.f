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
!    SUBROUTINE HYDROL-------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYDROL (
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,
     &                   NPNTS,NSHYD,B,CAN_CPY,CON_RAIN,CON_SNOW,
     &                   E_CANOPY,EXT,HCAP,HCON,INFIL_FAC,LS_RAIN,
     &                   LS_SNOW,ROOTD,SATCON,SATHH,SNOWSUB,
     &                   SURF_HT_FLUX,TIMESTEP,VFRAC,V_SAT,V_WILT,
     &                   TSTAR,RGRAIN,L_SNOW_ALBEDO,
     &                   CAN_WCNT,HF_SNOW_MELT,STF_HF_SNOW_MELT,
     &                   SMCL,STHF,STHU,SNOW_DEP,TSOIL,
     &                   INFIL,SMC,SNOW_MELT,SNOMLT_SUB_HTF,
     &                   STF_SUB_SURF_ROFF,SUB_SURF_ROFF,SURF_ROFF,
     &                   TOT_TFALL,LTIMER
     &)

      IMPLICIT NONE
!
! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with the Penman-Monteith
!     surface flux scheme. Calls the following:
!
!     HEAT_CAP - to calculate the heat capacity of the top
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
!     SMC_ROOT - to diagnose the available soil moisture in the
!                rootzone - Temporarary                (Cox 2/96)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1      6/96     New deck.  Peter Cox
!  4.4      24/11/97 Don't call SOIL_HTC/ICE_HTC if there are no
!                    soil/ice points. S.D.Mullerworth
!  4.4  29/10/97     Modified for prognostic snow albedo scheme
!                                                     R. Essery
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

      REAL
     & B(NPNTS)            ! IN Clapp-Hornberger exponent.
     &,CAN_CPY(NPNTS)      ! IN Canopy capacity (kg/m2).
     &,CON_RAIN(NPNTS)     ! IN Convective rain (kg/m2/s).
     &,CON_SNOW(NPNTS)     ! IN Convective snowfall (kg/m2/s).
     &,E_CANOPY(NPNTS)     ! IN Canopy evaporation (kg/m2/s).
     &,EXT(NPNTS,NSHYD)    ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
     &,HCAP(NPNTS)         ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS)         ! IN Soil thermal conductivity (W/m/K).
     &,INFIL_FAC(NPNTS)    ! IN Infiltration enhancement factor.
     &,LS_RAIN(NPNTS)      ! IN Large-scale rain (kg/m2/s).
     &,LS_SNOW(NPNTS)      ! IN Large-scale snowfall (kg/m2/s).
     &,ROOTD(NPNTS)        ! IN Rootdepth (m).
     &,SATCON(NPNTS)       ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
     &,SATHH(NPNTS)        ! IN Saturated soil water pressure (m).
     &,SNOWSUB(NPNTS)      ! IN Sublimation of lying snow (kg/m2/s).
     &,SURF_HT_FLUX(NPNTS) ! IN Net downward surface heat flux (W/m2).
     &,TSTAR(NPNTS)        ! IN Surface temperature (K).   
     &,VFRAC(NPNTS)        ! IN Vegetated fraction.
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
     & CAN_WCNT(NPNTS)     ! INOUT Canopy water content (kg/m2).
     &,HF_SNOW_MELT(NPNTS) ! INOUT Total snowmelt heat flux (W/m2).
     &,RGRAIN(NPNTS)       ! INOUT Snow grain size (microns).
     &,SMCL(NPNTS,NSHYD)   ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
     &,SNOW_DEP(NPNTS)     ! INOUT Snowmass (kg/m2).
     &,SNOW_MELT(NPNTS)    ! INOUT Snowmelt (kg/m2/s). 
     &,STHF(NPNTS,NSHYD)   ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)   ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,TSOIL(NPNTS,NSHYD)  ! INOUT Sub-surface temperatures (K).


!   Array arguments with intent(OUT) :
      REAL
     & INFIL(NPNTS)         ! OUT Maximum surface infiltration
!                                 rate (kg/m2/s).
     &,SMC(NPNTS)           ! OUT Available soil moisture in the
!                                 rootzone (kg/m2).



     &,SNOMLT_SUB_HTF(NPNTS)! OUT Sub-surface snowmelt heat
!                                 flux (W/m2).
     &,SUB_SURF_ROFF(NPNTS) ! OUT Sub-surface runoff (kg/m2/s).
     &,SURF_ROFF(NPNTS)     ! OUT Surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)     ! OUT Total throughfall (kg/m2/s).

! Local scalars:
      INTEGER
     & I                    ! WORK Loop counter.


! Local arrays:
      INTEGER
     & F_TYPE(NPNTS)        ! WORK Plant functional type:
!                           !       1 - Broadleaf Tree
!                           !       2 - Needleleaf Tree
!                           !       3 - C3 Grass
!                           !       4 - C4 Grass

      REAL
     & ASOIL(NPNTS)         ! WORK Reciprocal areal heat capacity
!                                  of the top soil layer (K m2/J).
     &,DSMC_DT(NPNTS)       ! WORK Rate of change of soil moisture
!                                  due to water falling onto the
!                                  surface after surface runoff
!                                  (kg/m2/s).
     &,HCAPS(NPNTS)         ! WORK Soil thermal capacity
!                                  including the effects of water
!                                  and ice (W/m/K).
     &,SOIL_HT_FLUX12(NPNTS)! WORK Heat flux between soil layers
!                                  1 and 2 (W/m2).
     &,V_ROOT(NPNTS)        ! WORK Volumetric soil moisture
!                           !      concentration in the rootzone
!                           !      (m3 H2O/m3 soil).
     &,V_SOIL(NPNTS)        ! WORK Volumetric soil moisture
!                           !      concentration in the top
!                           !      soil layer (m3 H2O/m3 soil).
     &,W_FLUX(NPNTS,0:NSHYD)! WORK Fluxes of water between layers
!                                  (kg/m2/s).
     &,WT_EXT(NPNTS,NSHYD)  ! WORK Fraction of transpiration which
!                                  is extracted from each soil layer.
! Function & Subroutine calls:
      EXTERNAL
     & HEAT_CAP,SFSNOW,INFILT,SURF_HYD,SOIL_HYD,SOIL_HTC,ICE_HTC

! End of header--------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('HYDROL ',103)
      ENDIF
!----------------------------------------------------------------------
! Calculate the heat capacity of the top soil layer
!----------------------------------------------------------------------
      CALL HEAT_CAP (NPNTS,SOIL_PTS,SOIL_INDEX,B,DZSOIL,
     &               HCAP,SATHH,SMCL,STHF,TSOIL,V_SAT,HCAPS,LTIMER)

!----------------------------------------------------------------------
! Calculate the reciprocal areal heat capacity of the top soil layer
!----------------------------------------------------------------------
      DO I=1,NPNTS
        ASOIL(I)=1./(DZSOIL(1)*HCAPS(I))
      ENDDO

!----------------------------------------------------------------------
! Calculate the subsurface snowmelt and update the snow mass
!----------------------------------------------------------------------
      CALL SFSNOW (ASOIL,CON_SNOW,LS_SNOW,SNOWSUB,TSTAR,
     &             TIMESTEP,NPNTS,SNOW_DEP,RGRAIN,L_SNOW_ALBEDO,
     &             TSOIL,SNOW_MELT,SNOMLT_SUB_HTF,STF_HF_SNOW_MELT,
     &             LTIMER)

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
      CALL INFILT (NPNTS,SOIL_PTS,SOIL_INDEX,B,SATCON,INFIL_FAC,STHF,
     &             INFIL,LTIMER)

!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------
      CALL SURF_HYD (NPNTS,E_CANOPY,SNOW_MELT,LS_RAIN,
     &               CON_RAIN,DSMC_DT,SURF_ROFF,CAN_WCNT,
     &               CAN_CPY,INFIL,TOT_TFALL,TIMESTEP)

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
     &               DZSOIL,HCAP,HCON,SNOW_DEP,SATHH,
     &               SURF_HT_FLUX,TIMESTEP,V_SAT,
     &               W_FLUX,SMCL,STHU,STHF,TSOIL,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
! Update the sub-surface temperatures for land ice
!-----------------------------------------------------------------------
      IF (LICE_PTS.NE.0)THEN
        CALL ICE_HTC (NPNTS,NSHYD,LICE_PTS,LICE_INDEX,DZSOIL,
     &              SURF_HT_FLUX,TIMESTEP,
     &              TSOIL,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
! Diagnose the plant functional types at each location (temporary).
! Assume : Broadleaf Trees if rootdepth > 0.8m
!          C3 Grass        if rootdepth < 0.8m
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        IF (ROOTD(I).GT.0.8) THEN
          F_TYPE(I)=1
        ELSE
          F_TYPE(I)=3
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Diagnose the soil moisture in the root zone.
!-----------------------------------------------------------------------
      CALL SMC_ROOT (NPNTS,NSHYD,F_TYPE,DZSOIL,ROOTD,STHU,VFRAC,        
     &               V_SAT,V_WILT,SMC,V_ROOT,V_SOIL,WT_EXT,LTIMER)


      IF (LTIMER) THEN
        CALL TIMER('HYDROL ',104)
      ENDIF

      RETURN
      END
