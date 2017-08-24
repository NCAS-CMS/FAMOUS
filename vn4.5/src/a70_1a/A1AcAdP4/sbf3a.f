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
!+ Subroutine to calculate the fluxes within the band using FESFT.
!
! Method:
!       The fluxes in the band including the grey processes are
!       calculated. Effective transmissions are found for each of the
!       gases by calculating the flux including the gas dividing
!       these by the grey fluxes. These effective transmissions
!       are then used to scale the grey fluxes. This is as described
!       by Geleyn and Ritter (1992).
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Initialization of
!                                               partial fluxes.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code for vertically
!                                               coherent convective
!                                               cloud added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_FESFT(IERR
!                       Atmospheric Column
     &   , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &   , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Propeties
     &   , I_BAND, N_GAS
     &   , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &   , K_ESFT, W_ESFT, SCALE_VECTOR
     &   , P_REFERENCE, T_REFERENCE
     &   , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &   , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &   , ISOLIR
!                       Solar Properties
     &   , SEC_0, SOLAR_FLUX
!                       Infra-red Properties
     &   , PLANCK_SOURCE_TOP, PLANCK_SOURCE_BOTTOM
     &   , DIFF_PLANCK_BAND
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &   , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &   , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &   , N_CLOUD_TOP
     &   , N_CLOUD_TYPE, FRAC_CLOUD
     &   , I_REGION_CLOUD, FRAC_REGION
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &   , CLOUD_OVERLAP
     &   , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &   , FLUX_DIRECT_BAND, FLUX_DIFFUSE_BAND
!                       Flags for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_DIFFUSE_CLEAR_BAND
!                       Planckian Function
     &   , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , NPD_BAND, NPD_SPECIES
     &   , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_BAND
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_SPECIES
!             MAXIMUM NUMBER OF SPECIES
     &   , NPD_ESFT_TERM
!             MAXIMUM NUMBER OF ESFT TERMS
     &   , NPD_SCALE_VARIABLE
!             MAXIMUM NUMBER OF SCALE VARIABLES
     &   , NPD_SCALE_FNC
!             MAXIMUM NUMBER OF SCALING FUNCTIONS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     INCLUDE COMDECKS.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INTERNAL DIMENSIONS TIED TO ALGORITHMS,
!     MOSTLY FOR CLOUDS.
!
      INTEGER
     &     NPD_CLOUD_COMPONENT
!             NUMBER OF COMPONENTS OF CLOUDS
     &   , NPD_CLOUD_TYPE
!             NUMBER OF PERMITTED TYPES OF CLOUDS.
     &   , NPD_CLOUD_REPRESENTATION
!             NUMBER OF PERMITTED REPRESENTATIONS OF CLOUDS.
     &   , NPD_OVERLAP_COEFF
!             NUMBER OF OVERLAP COEFFICIENTS FOR CLOUDS
     &   , NPD_SOURCE_COEFF
!             NUMBER OF COEFFICIENTS FOR TWO-STREAM SOURCES
     &   , NPD_REGION
!             NUMBER OF REGIONS IN A LAYER
!
      PARAMETER(
     &     NPD_CLOUD_COMPONENT=4
     &   , NPD_CLOUD_TYPE=4
     &   , NPD_CLOUD_REPRESENTATION=4
     &   , NPD_OVERLAP_COEFF=18
     &   , NPD_SOURCE_COEFF=2
     &   , NPD_REGION=3
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE SETTING TYPES OF ESFT SCALING
!
      INTEGER
     &     IP_SCALE_NULL
!             NO SCALING AT ALL
     &   , IP_SCALE_BAND
!             SAME SCALING THROUGHOUT BAND
     &   , IP_SCALE_TERM
!             DIFFERENT SCALING FOR EACH ESFT
!
      PARAMETER(
     &     IP_SCALE_NULL=0
     &   , IP_SCALE_BAND=1
     &   , IP_SCALE_TERM=2
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET FLAGS FOR DIFFERENT PORTIONS
!     OF THE SPECTRUM.
!
      INTEGER
     &     IP_SOLAR
!             SOLAR REGION
     &   , IP_INFRA_RED
!             INFRA-RED REGION
!
      PARAMETER(
     &     IP_SOLAR=1
     &   , IP_INFRA_RED=2
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP OF VERTICAL GRID
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             PROPERTIES GIVEN IN LAYERS
      REAL  !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
!
!                       Angular Integration
      INTEGER   !, INTENT(IN)
     &     I_ANGULAR_INTEGRATION
!             ANGULAR INTEGRATION SCHEME
     &   , I_2STREAM
!             TWO-STREAM SCHEME
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
      LOGICAL   !, INTENT(IN)
     &     L_2_STREAM_CORRECT
!             USE AN EDGE CORRECTION
     &   , L_RESCALE
!             RESCALE OPTICAL PROPERTIES
!
!                       Treatment of Scattering
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             METHOD OF TREATING SCATTERING
!
!                       Options for solver
      INTEGER   !, INTENT(IN)
     &     I_SOLVER
!             SOLVER USED
     &   , N_AUGMENT
!             LENGTH OF LONG FLUX VECTOR
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES
!
!                       Gaseous Properties
      INTEGER   !, INTENT(IN)
     &     I_BAND
!             BAND BEING CONSIDERED
     &   , N_GAS
!             NUMBER OF GASES IN BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN BANDS
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)
!             NUMBER OF TERMS IN BAND
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
      LOGICAL   !, INTENT(IN)
     &     L_DOPPLER(NPD_SPECIES)
!             DOPPLER BROADENING INCLUDED
      REAL  !, INTENT(IN)
     &     K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             EXPONENTIAL ESFT TERMS
     &   , W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             WEIGHTS FOR ESFT
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND
     &        , NPD_SPECIES)
!             ABSORBER SCALING PARAMETERS
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE SCALING PRESSURE
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE SCALING TEMPERATURE
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             GAS MASS MIXING RATIOS
     &   , GAS_FRAC_RESCALED(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             RESCALED GAS MASS FRACTIONS
     &   , DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER BROADENING TERMS
!
!                       Spectral Region
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             VISIBLE OR IR
!
!                       Solar Properties
      REAL  !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF SOLAR ZENITH ANGLE
     &   , SOLAR_FLUX(NPD_PROFILE)
!             INCIDENT SOLAR FLUX IN BAND
!
!                       Infra-red Properties
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE FUNCTION
      REAL  !, INTENT(IN)
     &     PLANCK_SOURCE_TOP(NPD_PROFILE)
!             PLANCKIAN SOURCE AT TOP
     &   , PLANCK_SOURCE_BOTTOM(NPD_PROFILE)
!             PLANCKIAN SOURCE AT BOTTOM
     &   , DIFF_PLANCK_BAND(NPD_PROFILE, NPD_LAYER)
!             THERMAL SOURCE FUNCTION
     &   , DIFF_PLANCK_BAND_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCE OF PLANCKIAN IN BAND
!
!                       Surface Properties
      REAL  !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
     &   , THERMAL_GROUND_BAND(NPD_PROFILE)
!             THERMAL SOURCE FUNCTION AT GROUND
!
!                       Clear-sky Optical Properties
      REAL  !, INTENT(IN)
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE SCATTERING EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY ASYMMETRY
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE FORWARD SCATTERING
!
!                       Cloudy Properties
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD
!             CLOUDS REQUIRED
      INTEGER   !, INTENT(IN)
     &     I_CLOUD
!             CLOUD SCHEME USED
!
!                       Cloud Geometry
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP
!             TOP CLOUDY LAYER
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
     &   , N_FREE_PROFILE(NPD_LAYER)
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)
!             INDICES OF FREE PROFILES
     &   , N_CLOUD_PROFILE(NPD_LAYER)
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             INDICES OF CLOUDY PROFILES
     &   , N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS REQUIRED
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
      LOGICAL   !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             COLUMN FLAGS FOR COLUMNS
      REAL  !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FRACTION
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF TYPES OF CLOUD
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY FRACTION
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER FOR ENERGY AT INTERFACES
     &   , AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREAS OF COLUMNS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!                       Cloudy Optical Propeties
      REAL  !, INTENT(IN)
     &     K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SCATTERING EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ASYMMETRY
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
!                       Fluxes Calculated
      REAL  !, INTENT(OUT)
     &     FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX IN BAND
     &   , FLUX_DIFFUSE_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUX IN BAND
!
!                       Flags for Clear-sky Calculations
      LOGICAL   !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE CLEAR-SKY FLUXES
      INTEGER   !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Clear-sky Fluxes Calculated
      REAL  !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DIRECT FLUX IN BAND
     &   , FLUX_DIFFUSE_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR-SKY DIFFUSE FLUX IN BAND
!
!                       Planckian Source Function
      REAL  !, INTENT(IN)
     &     PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE IN BAND
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      INTEGER
     &     I_GAS_BAND
!             INDEX OF ACTIVE GAS
     &   , I_GAS_POINTER(NPD_SPECIES)
!             POINTER ARRAY FOR MONOCHROMATIC ESFTs
     &   , IEX
!             INDEX OF ESFT TERM
      REAL
     &     K_ESFT_MONO(NPD_SPECIES)
!             ESFT MONOCHROMATIC EXPONENTS
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS ABSORPTION
     &   , SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD FLUX
     &   , ESFT_WEIGHT
!             ESFT WEIGHT FOR CURRENT CALCULATION
     &   , DUMMY_KE(NPD_PROFILE, NPD_LAYER)
!             DUMMY ARRAY (NOT USED)
      REAL
     &     FLUX_DIRECT_PART(NPD_PROFILE, 0: NPD_LAYER)
!             PARTIAL DIRECT FLUX
     &   , FLUX_DIFFUSE_PART(NPD_PROFILE, 2*NPD_LAYER+2)
!             PARTIAL DIFFUSE FLUX
     &   , FLUX_DIRECT_CLEAR_PART(NPD_PROFILE, 0: NPD_LAYER)
!             PARTIAL CLEAR-SKY DIRECT FLUX
     &   , FLUX_DIFFUSE_CLEAR_PART(NPD_PROFILE, 2*NPD_LAYER+2)
!             PARTIAL CLEAR-SKY DIFFUSE FLUX
     &   , FLUX_DIRECT_GREY(NPD_PROFILE, 0: NPD_LAYER)
!             GREY DIRECT FLUX
     &   , FLUX_DIFFUSE_GREY(NPD_PROFILE, 2*NPD_LAYER+2)
!             GREY DIFFUSE FLUX
     &   , FLUX_DIRECT_CLEAR_GREY(NPD_PROFILE, 0: NPD_LAYER)
!             GREY CLEAR-SKY DIRECT FLUX
     &   , FLUX_DIFFUSE_CLEAR_GREY(NPD_PROFILE, 2*NPD_LAYER+2)
!             GREY CLEAR-SKY DIFFUSE FLUX
     &   , FLUX_RATIO_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             RATIO OF DIRECT FLUXES
     &   , FLUX_RATIO_DIFFUSE(NPD_PROFILE, 2*NPD_LAYER+2)
!             RATIO OF DIFFUSE FLUXES
     &   , FLUX_RATIO_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             RATIO OF CLEAR-SKY DIRECT FLUXES
     &   , FLUX_RATIO_DIFFUSE_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             RATIO OF CLEAR-SKY DIFFUSE FLUXES
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     INITIALIZE_FLUX, SCALE_ABSORB, GAS_OPTICAL_PROPERTIES
     &   , MONOCHROMATIC_FLUX, AUGMENT_FLUX
!
!
!
!     FAST EXPONENTIAL OVERLAP, SUPERPOSING ONE GAS AT A TIME.
!
!     INITIAL SOLUTION FOR GREY FLUXES.
!     SET THE APPROPRIATE BOUNDARY CONDITIONS FOR THE TOTAL
!     UPWARD AND DOWNWARD FLUXES.
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
!        VISIBLE REGION.
         DO L=1, N_PROFILE
            SOURCE_GROUND(L)=0.0E+00
            FLUX_INC_DIRECT(L)=SOLAR_FLUX(L)
            FLUX_INC_DOWN(L)=SOLAR_FLUX(L)
         ENDDO
      ELSEIF (ISOLIR.EQ.IP_INFRA_RED) THEN
!        INFRA-RED REGION.
         DO L=1, N_PROFILE
            FLUX_INC_DIRECT(L)=0.0E+00
            FLUX_DIRECT_GREY(L, N_LAYER)=0.0E+00
            FLUX_INC_DOWN(L)=-PLANCK_SOURCE_TOP(L)
            SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)
     &         -(1.0E+00-ALBEDO_SURFACE_DIFF(L))
     &         *PLANCK_SOURCE_BOTTOM(L)
         ENDDO
         IF (L_CLEAR) THEN
            DO L=1, N_PROFILE
               FLUX_DIRECT_CLEAR_GREY(L, N_LAYER)=0.0E+00
            ENDDO
         ENDIF
      ENDIF
!
!     THIS CALL CONTAINS NO GASEOUS ABSORPTION SO K_GAS_ABS MUST
!     BE ZEROED.
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            K_GAS_ABS(L, I)=0.0E+00
         ENDDO
      ENDDO
!
!
      CALL MONOCHROMATIC_FLUX(IERR
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER, D_MASS
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &   , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &   , K_GAS_ABS
!                       Options for Equivalent Extinction
     &   , .FALSE., DUMMY_KE
!                       Spectral Region
     &   , ISOLIR
!                       Infra-red Properties
     &   , DIFF_PLANCK_BAND
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Conditions at TOA
     &   , SEC_0, SOLAR_FLUX, FLUX_INC_DOWN
!                       Surface Properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , SOURCE_GROUND, THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &   , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &   , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &   , N_CLOUD_TOP
     &   , N_CLOUD_TYPE, FRAC_CLOUD
     &   , I_REGION_CLOUD, FRAC_REGION
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &   , CLOUD_OVERLAP
     &   , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &   , FLUX_DIRECT_GREY, FLUX_DIFFUSE_GREY
!                       Flags for Clear-sky Calculations
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR_GREY, FLUX_DIFFUSE_CLEAR_GREY
!                       Planckian Function
     &   , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!     FESFT IS BASED ON THE DIFFUSE FLUXES, WHEREAS TOTAL FLUXES
!     ARE USED HERE IN THE SOLAR: CONVERSION BACK TO DIFFUSE
!     FLUXES IS REQUIRED.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         IF (L_NET) THEN
            DO I=1, N_AUGMENT
               DO L=1, N_PROFILE
                  FLUX_DIFFUSE_GREY(L, I)=FLUX_DIFFUSE_GREY(L, I)
     &               -FLUX_DIRECT_GREY(L, I-1)
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_AUGMENT
               DO L=1, N_PROFILE
                  FLUX_DIFFUSE_GREY(L, 2*I)=FLUX_DIFFUSE_GREY(L, 2*I)
     &               -FLUX_DIRECT_GREY(L, I-1)
               ENDDO
            ENDDO
         ENDIF
         IF (L_CLEAR) THEN
            IF (L_NET) THEN
               DO I=1, N_AUGMENT
                  DO L=1, N_PROFILE
                     FLUX_DIFFUSE_CLEAR_GREY(L, I)
     &                  =FLUX_DIFFUSE_CLEAR_GREY(L, I)
     &                  -FLUX_DIRECT_CLEAR_GREY(L, I-1)
                  ENDDO
               ENDDO
            ELSE
               DO I=1, N_AUGMENT
                  DO L=1, N_PROFILE
                     FLUX_DIFFUSE_CLEAR_GREY(L, 2*I)
     &                  =FLUX_DIFFUSE_CLEAR_GREY(L, 2*I)
     &                  -FLUX_DIRECT_CLEAR_GREY(L, I-1)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!
!
!     THE FLUX RATIOS ARE THE RATIOS OF THE GASEOUS FLUXES TO
!     THE GREY FLUXES. THE PRODUCT OVER ALL GASES IS USED TO
!     CALCULATE THE OVERALL FLUX.
      CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &   , ISOLIR
     &   , FLUX_RATIO_DIRECT, FLUX_RATIO_DIFFUSE
     &   , L_CLEAR
     &   , FLUX_RATIO_DIRECT_CLEAR, FLUX_RATIO_DIFFUSE_CLEAR
     &   , 1.0E+00
     &   , NPD_PROFILE, NPD_LAYER
     &   , L_NET
     &   )
!
      DO J=1, N_GAS
!
!        INITIALIZE THE FLUX IN THE BAND TO ZERO. IN THIS
!        LOOP FLUX_..._BAND IS USED AS A TEMPORARY VARIABLE
!        TO HOLD THE FLUXES FOR ONE GAS.
         CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &      , ISOLIR
     &      , FLUX_DIRECT_BAND, FLUX_DIFFUSE_BAND
     &      , L_CLEAR
     &      , FLUX_DIRECT_CLEAR_BAND, FLUX_DIFFUSE_CLEAR_BAND
     &      , 0.0E+00
     &      , NPD_PROFILE, NPD_LAYER
     &      , L_NET
     &      )
!
         I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
         I_GAS_POINTER(1)=I_GAS_BAND
         DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS_BAND)
!
!           STORE THE ESFT WEIGHT FOR FUTURE USE.
            ESFT_WEIGHT=W_ESFT(IEX, I_BAND,  I_GAS_BAND)
!
!           RESCALE THE AMOUNT OF GAS FOR THIS ABSORBER IF REQUIRED.
            IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND).EQ.IP_SCALE_TERM) THEN
               CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &            , GAS_MIX_RATIO(1, 0, I_GAS_BAND), P, T
     &            , L_LAYER, I_TOP
     &            , GAS_FRAC_RESCALED(1, 0, I_GAS_BAND)
     &            , I_SCALE_FNC(I_BAND, I_GAS_BAND)
     &            , P_REFERENCE(I_GAS_BAND, I_BAND)
     &            , T_REFERENCE(I_GAS_BAND, I_BAND)
     &            , SCALE_VECTOR(1, IEX, I_BAND, I_GAS_BAND)
     &            , L_DOPPLER(I_GAS_BAND)
     &            , DOPPLER_CORRECTION(I_GAS_BAND)
     &            , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &            , NPD_SCALE_VARIABLE
     &            )
               IF (IERR.NE.I_NORMAL) RETURN
            ENDIF
!
!           SET THE APPROPRIATE SOURCE TERMS FOR THE TWO-STREAM
!           EQUATIONS.
!           SET THE GROUND SOURCE FUNCTION
            IF (ISOLIR.EQ.IP_SOLAR) THEN
!              VISIBLE REGION.
               DO L=1, N_PROFILE
                  SOURCE_GROUND(L)=0.0E+00
                  FLUX_INC_DOWN(L)=SOLAR_FLUX(L)
                  FLUX_INC_DIRECT(L)=SOLAR_FLUX(L)
               ENDDO
            ELSEIF (ISOLIR.EQ.IP_INFRA_RED) THEN
!              INFRA-RED REGION.
               DO L=1, N_PROFILE
                  FLUX_INC_DIRECT(L)=0.0E+00
                  FLUX_DIRECT_PART(L, N_LAYER)=0.0E+00
                  FLUX_INC_DOWN(L)=-PLANCK_SOURCE_TOP(L)
                  SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)
     &               -(1.0E+00-ALBEDO_SURFACE_DIFF(L))
     &               *PLANCK_SOURCE_BOTTOM(L)
               ENDDO
               IF (L_CLEAR) THEN
                  DO L=1, N_PROFILE
                     FLUX_DIRECT_CLEAR_PART(L, N_LAYER)=0.0E+00
                  ENDDO
               ENDIF
            ENDIF
!
!           ASSIGN THE MONOCHROMATIC ABSORPTION COEFFICIENT.
            K_ESFT_MONO(I_GAS_BAND)=K_ESFT(IEX, I_BAND, I_GAS_BAND)
!
            CALL GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER
     &         , 1, I_GAS_POINTER, K_ESFT_MONO
     &         , GAS_FRAC_RESCALED
     &         , K_GAS_ABS
     &         , NPD_PROFILE, NPD_LAYER, NPD_SPECIES
     &         )
!
            CALL MONOCHROMATIC_FLUX(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER, D_MASS
!                       Angular Integration
     &         , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &         , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &         , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &         , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Propreties
     &         , K_GAS_ABS
!                       Options for Equivalent Extinction
     &         , .FALSE., DUMMY_KE
!                       Spectral Region
     &         , ISOLIR
!                       Infra-red Properties
     &         , DIFF_PLANCK_BAND
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Conditions at TOA
     &         , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN
!                       Surface Properties
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
     &         , THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &         , K_GREY_TOT_FREE, K_EXT_SCAT_FREE
     &         , ASYMMETRY_FREE, FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &         , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &         , N_CLOUD_TOP
     &         , N_CLOUD_TYPE, FRAC_CLOUD
     &        , I_REGION_CLOUD, FRAC_REGION
     &         , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &         , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &         , CLOUD_OVERLAP
     &         , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &         , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &         , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Flxues Calculated
     &         , FLUX_DIRECT_PART, FLUX_DIFFUSE_PART
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR_PART, FLUX_DIFFUSE_CLEAR_PART
!                       Planckian Function
     &         , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
!
!           CONVERT THE TOTAL FLUXES BACK TO DIFFUSE FLUXES.
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               IF (L_NET) THEN
                  DO I=1, N_AUGMENT
                     DO L=1, N_PROFILE
                        FLUX_DIFFUSE_PART(L, I)
     &                     =FLUX_DIFFUSE_PART(L, I)
     &                     -FLUX_DIRECT_PART(L, I-1)
                     ENDDO
                  ENDDO
               ELSE
                  DO I=1, N_AUGMENT
                     DO L=1, N_PROFILE
                        FLUX_DIFFUSE_PART(L, 2*I)
     &                     =FLUX_DIFFUSE_PART(L, 2*I)
     &                     -FLUX_DIRECT_PART(L, I-1)
                     ENDDO
                  ENDDO
               ENDIF
               IF (L_CLEAR) THEN
                  IF (L_NET) THEN
                     DO I=1, N_AUGMENT
                        DO L=1, N_PROFILE
                           FLUX_DIFFUSE_CLEAR_PART(L, I)
     &                        =FLUX_DIFFUSE_CLEAR_PART(L, I)
     &                        -FLUX_DIRECT_CLEAR_PART(L, I-1)
                        ENDDO
                     ENDDO
                  ELSE
                     DO I=1, N_AUGMENT
                        DO L=1, N_PROFILE
                           FLUX_DIFFUSE_CLEAR_PART(L, 2*I)
     &                        =FLUX_DIFFUSE_CLEAR_PART(L, 2*I)
     &                        -FLUX_DIRECT_CLEAR_PART(L, I-1)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
!
!           INCREMENT THE FLUXES WITHIN THE BAND.
            CALL AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &         , ISOLIR, L_CLEAR
     &         , ESFT_WEIGHT
     &         , FLUX_DIRECT_BAND, FLUX_DIFFUSE_BAND
     &         , FLUX_DIRECT_PART, FLUX_DIFFUSE_PART
     &         , FLUX_DIRECT_CLEAR_BAND, FLUX_DIFFUSE_CLEAR_BAND
     &         , FLUX_DIRECT_CLEAR_PART, FLUX_DIFFUSE_CLEAR_PART
     &         , NPD_PROFILE, NPD_LAYER
     &         )
!
          ENDDO
!
!         CALCULATE THE FLUX RATIO.
          IF (ISOLIR.EQ.IP_SOLAR) THEN
             DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  IF (ABS(FLUX_DIRECT_GREY(L, I)).GT.0.0E+00) THEN
                     FLUX_RATIO_DIRECT(L, I)=FLUX_RATIO_DIRECT(L, I)
     &                  *FLUX_DIRECT_BAND(L, I)/FLUX_DIRECT_GREY(L, I)
                  ELSE
                     FLUX_RATIO_DIRECT(L, I)=1.0E+00
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         DO I=1, N_AUGMENT
            DO L=1, N_PROFILE
               IF (ABS(FLUX_DIFFUSE_GREY(L, I)).GT.0.0E+00) THEN
                  FLUX_RATIO_DIFFUSE(L, I)=FLUX_RATIO_DIFFUSE(L, I)
     &               *FLUX_DIFFUSE_BAND(L, I)/FLUX_DIFFUSE_GREY(L, I)
               ELSE
                  FLUX_RATIO_DIFFUSE(L, I)=1.0E+00
               ENDIF
            ENDDO
         ENDDO
!
         IF (L_CLEAR) THEN
!           CALCULATE THE FLUX RATIOS FOR THE CLEAR FLUXES.
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO I=0, N_LAYER
                  DO L=1, N_PROFILE
                     IF (ABS(FLUX_DIRECT_CLEAR_GREY(L, I))
     &                  .GT.0.0E+00) THEN
                        FLUX_RATIO_DIRECT_CLEAR(L, I)
     &                     =FLUX_RATIO_DIRECT_CLEAR(L, I)
     &                     *FLUX_DIRECT_CLEAR_BAND(L, I)
     &                     /FLUX_DIRECT_CLEAR_GREY(L, I)
                     ELSE
                        FLUX_RATIO_DIRECT_CLEAR(L, I)=1.0E+00
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
!
            DO I=1, N_AUGMENT
               DO L=1, N_PROFILE
                  IF (ABS(FLUX_DIFFUSE_CLEAR_GREY(L, I))
     &               .GT.0.0E+00) THEN
                     FLUX_RATIO_DIFFUSE_CLEAR(L, I)
     &                  =FLUX_RATIO_DIFFUSE_CLEAR(L, I)
     &                  *FLUX_DIFFUSE_CLEAR_BAND(L, I)
     &                  /FLUX_DIFFUSE_CLEAR_GREY(L, I)
                  ELSE
                     FLUX_RATIO_DIFFUSE_CLEAR(L, I)=1.0E+00
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
      ENDDO
!
!     THE OVERALL FLUX IN THE BAND IS CALCULATED FROM THE
!     GREY FLUX AND THE FLUX RATIOS.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT_BAND(L, I)=FLUX_RATIO_DIRECT(L, I)
     &            *FLUX_DIRECT_GREY(L, I)
            ENDDO
         ENDDO
      ENDIF
      DO I=1, N_AUGMENT
         DO L=1, N_PROFILE
            FLUX_DIFFUSE_BAND(L, I)=FLUX_RATIO_DIFFUSE(L, I)
     &         *FLUX_DIFFUSE_GREY(L, I)
         ENDDO
      ENDDO
      IF (L_CLEAR) THEN
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR_BAND(L, I)
     &               =FLUX_RATIO_DIRECT_CLEAR(L, I)
     &               *FLUX_DIRECT_CLEAR_GREY(L, I)
               ENDDO
            ENDDO
         ENDIF
         DO I=1, N_AUGMENT
            DO L=1, N_PROFILE
               FLUX_DIFFUSE_CLEAR_BAND(L, I)
     &            =FLUX_RATIO_DIFFUSE_CLEAR(L, I)
     &            *FLUX_DIFFUSE_CLEAR_GREY(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
!     CONVERT BACK TO TOTAL FLUXES.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         IF (L_NET) THEN
            DO I=1, N_AUGMENT
               DO L=1, N_PROFILE
                  FLUX_DIFFUSE_BAND(L, I)=FLUX_DIFFUSE_BAND(L, I)
     &               +FLUX_DIRECT_BAND(L, I-1)
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_AUGMENT
               DO L=1, N_PROFILE
                  FLUX_DIFFUSE_BAND(L, 2*I)=FLUX_DIFFUSE_BAND(L, 2*I)
     &               +FLUX_DIRECT_BAND(L, I-1)
               ENDDO
            ENDDO
         ENDIF
         IF (L_CLEAR) THEN
            IF (L_NET) THEN
               DO I=1, N_AUGMENT
                  DO L=1, N_PROFILE
                     FLUX_DIFFUSE_CLEAR_BAND(L, I)
     &                  =FLUX_DIFFUSE_CLEAR_BAND(L, I)
     &                  +FLUX_DIRECT_CLEAR_BAND(L, I-1)
                  ENDDO
               ENDDO
            ELSE
               DO I=1, N_AUGMENT
                  DO L=1, N_PROFILE
                     FLUX_DIFFUSE_CLEAR_BAND(L, 2*I)
     &                  =FLUX_DIFFUSE_CLEAR_BAND(L, 2*I)
     &                  +FLUX_DIRECT_CLEAR_BAND(L, I-1)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!
      RETURN
      END
