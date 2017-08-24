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
!+ Subroutine to solve for the monochromatic fluxes.
!
! Method:
!       The final single scattering properties are calculated
!       and rescaled. An appropriate subroutine is called to
!       calculate the fluxes depending on the treatment of
!       cloudiness.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code for vertically
!                                               coherent cloud added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed.
!                                               Unused variables
!                                               removed from call
!                                               to TRPILE_COLUMN.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_FLUX(IERR
!                       Atmospheric Propetries
     &   , N_PROFILE, N_LAYER, D_MASS
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &   , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Propeties
     &   , K_GAS_ABS
!                       Options for Equivalent Extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &   , ISOLIR
!                       Infra-red Properties
     &   , DIFF_PLANCK
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &   , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN
!                       Surface Propeties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
     &   , GROUND_EMISSION
!                       Clear-sky Optical Propeties
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
!                       Cloudy Optical Propeties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculation
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Planckian Source
     &   , PLANCK_SOURCE
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   )
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
!     MODULE TO DEFINE REFERENCE NUMBERS FOR CLOUD SCHEMES.
!
      INTEGER
     &     IP_CLOUD_MIX_MAX
!             MAXIMUM/RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_MIX_RANDOM
!             RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_COLUMN_MAX
!             MAXIMUM OVERLAP IN A COLUMN MODEL
     &   , IP_CLOUD_CLEAR
!             CLEAR COLUMN
     &   , IP_CLOUD_TRIPLE
!             MIXED COLUMN WITH SPLIT BETWEEN 
!             CONVECTIVE AND LAYER CLOUD.
!
      PARAMETER(
     &     IP_CLOUD_MIX_MAX=2
     &   , IP_CLOUD_MIX_RANDOM=4
     &   , IP_CLOUD_COLUMN_MAX=3
     &   , IP_CLOUD_CLEAR=5
     &   , IP_CLOUD_TRIPLE=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET THE TYPES OF ANGULAR INTEGRATION.
!
      INTEGER
     &     IP_TWO_STREAM
!             TWO STREAM SCHEME
     &   , IP_IR_GAUSS
!             GAUSSIAN INTEGRATION IN THE IR
!
      PARAMETER(
     &     IP_TWO_STREAM=1
     &   , IP_IR_GAUSS=2
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
!
!                       Atmospheric Properties
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
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
!             CORRECTION TO TWO-STREAM SCHEME
     &   , L_RESCALE
!             RESCALE OPTICAL PROPERTIES
!
!                       Treatment of Scattering
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!
!                       Options for Solver
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
      REAL      !, INTENT(IN)
     &     K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS ABSORPTIVE EXTINCTIONS
!
!                       Variables for Equivalent Extinction
      LOGICAL   !, INTENT(IN)
     &     L_SCALE_SOLAR
!             APPLY SCALING TO SOLAR FLUX
      REAL      !, INTENT(IN)
     &     ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
!
!                       Spectral Region
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             VISIBLE OR IR
!
!                       Infra-red Properties
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             FLAG FOR QUADRATIC IR-SOURCE
      REAL      !, INTENT(IN)
     &     PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)
!             MONOCHROMATIC PLANCKIAN SOURCE
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             THERMAL SOURCE FUNCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2ND DIFF. OF THERMAL SOURCE FUNCTION
!
!                       Conditions at TOA
      REAL      !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF SOLAR ZENITH ANGLE
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD FLUX
!
!                       Surface Propeties
      REAL      !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
     &   , SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
      REAL      !, INTENT(IN)
     &     GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
!                       Optical Properties
      REAL      !, INTENT(IN)
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
!             TOPMOST CLOUDY LAYER
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
!             FLAGS FOR CONTENTS OF COLUMNS
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FRACTION
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY FRACTION
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR ENERGY TRANSFER AT INTERFACES
     &   , AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREAS OF COLUMNS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!                       Cloudy Optical Properties
      REAL      !, INTENT(IN)
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
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!                       Flags for Clear-sky Calculations
      LOGICAL   !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE CLEAR-SKY PROPERTIES
      INTEGER   !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Clear-sky Fluxes Calculated
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR-SKY TOTAL FLUX
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     K
!             LOOP VARIABLE
      REAL
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ALBEDO OF S. S.
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SINGLE SCATTERING ALBEDO
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     SINGLE_SCATTERING_ALL, RESCALE_TAU_OMEGA
     &   , TWO_STREAM
     &   , MIX_COLUMN, CLOUD_COLUMN
     &   , GAUSS_ANGLE
!
!
!
!     CALCULATE SINGLE SCATTERING PROPERTIES FOR ALL ATMOSPHERIC
!     CONSTITUENTS.
!
      CALL SINGLE_SCATTERING_ALL(I_SCATTER_METHOD_BAND
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER, D_MASS
!                       Cloudy Properties
     &   , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE
!                       Optical Properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , K_GAS_ABS
!                       Single Scattering Properties
     &   , TAU_FREE, OMEGA_FREE
     &   , TAU_CLOUD, OMEGA_CLOUD
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IF (I_ANGULAR_INTEGRATION.EQ.IP_TWO_STREAM) THEN
!
!        RESCALE TAU AND OMEGA. THE ASYMMETRY HAS ALREADY BEEN RESCALED.
!
         IF (L_RESCALE) THEN
!
            CALL RESCALE_TAU_OMEGA(N_PROFILE, 1, N_LAYER
     &         , TAU_FREE, OMEGA_FREE, FORWARD_SCATTER_FREE
     &         , NPD_PROFILE, NPD_LAYER
     &         )
!
            IF (L_CLOUD) THEN
!
               DO K=1, N_CLOUD_TYPE
                  CALL RESCALE_TAU_OMEGA(N_PROFILE, N_CLOUD_TOP
     &               , N_LAYER
     &               , TAU_CLOUD(1, 1, K), OMEGA_CLOUD(1, 1, K)
     &               , FORWARD_SCATTER_CLOUD(1, 1, K)
     &               , NPD_PROFILE, NPD_LAYER
     &               )
               ENDDO
!
            ENDIF
!
         ENDIF
!
!
!        SOLVE THE EQUATIONS USING THE SCHEME INDICATED BY THE VALUES
!        OF I_CLOUD AND I_SOLVER.
         IF (I_CLOUD.EQ.IP_CLOUD_CLEAR) THEN
!
!           A TWO-STREAM SCHEME WITH NO CLOUDS.
            CALL TWO_STREAM(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &         , I_2STREAM
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &         , L_NET, I_SOLVER
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &         , ISOLIR
!                       Infra-red Properties
     &         , DIFF_PLANCK
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Surface Conditions
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Single Scattering Propeties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL
!                       Flag for Clear-sky Fluxes
     &         , L_CLEAR
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Sizes of Arrays
     &         , NPD_PROFILE, NPD_LAYER
     &         )
               IF (IERR.NE.I_NORMAL) RETURN
!
         ELSEIF ( (I_CLOUD.EQ.IP_CLOUD_MIX_MAX)
     &          .OR. (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) ) THEN
!
!           CLOUDS ARE TREATED USING ZDUNKOWSKI'S MIXED-COLUMN SCHEME.
!           THE GEOMETRY HAS BEEN SET BEFORE.
!
            CALL MIX_COLUMN(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &         , I_2STREAM
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &         , I_SOLVER, L_NET
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &         , ISOLIR
!                       Infra-red Properties
     &         , DIFF_PLANCK
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE
!                       Cloud Geometry
     &         , N_CLOUD_TOP
     &         , N_CLOUD_TYPE, FRAC_CLOUD
     &         , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &         , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &         , CLOUD_OVERLAP
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
         ELSEIF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
!           CLOUDS ARE TREATED USING A DECOMPOSITION OF THE COLUMN
!           INTO CLEAR-SKY, STRATIFORM AND CONVECTIVE REGIONS, ALL
!           MAXIMALLY OVERLAPPED.
!
            CALL TRIPLE_COLUMN(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &         , I_2STREAM
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &         , I_SOLVER, L_NET
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &         , ISOLIR
!                       Infra-red Properties
     &         , DIFF_PLANCK
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE
!                       Cloud Geometry
     &         , N_CLOUD_TOP
     &         , N_CLOUD_TYPE, FRAC_CLOUD
     &         , I_REGION_CLOUD, FRAC_REGION
     &         , W_FREE, W_CLOUD
     &         , CLOUD_OVERLAP
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
!
         ELSEIF (I_CLOUD.EQ.IP_CLOUD_COLUMN_MAX) THEN
!           CLOUDS ARE TREATED ON THE ASSUMPTION OF MAXIMUM OVERLAP
!           IN A COLUMN MODEL.
            CALL CLOUD_COLUMN(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &         , I_2STREAM
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &         , I_SOLVER, N_AUGMENT
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &         , ISOLIR
!                       Infra-red Properties
     &         , DIFF_PLANCK
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE
!                       Cloud Geometry
     &         , N_CLOUD_TOP
     &         , N_CLOUD_TYPE, FRAC_CLOUD
     &         , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
!
         ENDIF
!
      ELSE IF (I_ANGULAR_INTEGRATION.EQ.IP_IR_GAUSS) THEN
!
!        FULL ANGULAR RESOLUTION USING GASUSSIAN INTEGRATION.
         CALL GAUSS_ANGLE(N_PROFILE, N_LAYER, L_NET, N_AUGMENT
     &      , N_ORDER_GAUSS
     &      , TAU_FREE
     &      , FLUX_INC_DOWN
     &      , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF
     &      , FLUX_TOTAL
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
!
      ENDIF
!
!
      RETURN
      END
