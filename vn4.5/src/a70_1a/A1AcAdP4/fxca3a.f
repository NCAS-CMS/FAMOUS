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
!+ Subroutine to calculate radiative fluxes.
!
! Method:
!       Properties independent of the spectral bands are set.
!       A loop over bands is then entered. Grey optical properties
!       are set and an appropriate subroutine is called to treat
!       the gaseous overlaps. The final fluxes are assigned.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New solvers added.
!                                               Revised formulation for
!                                               sea-ice. Renaming of
!                                               diagnostic and logical
!                                               switch for flux
!                                               below 690 nm. Pointer
!                                               to water vapour added.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Generalization for
!                                               vertically coherent
!                                               convective cloud.
!       4.4             30-09-96                Effective Radius
!                                               relabelled as charact-
!                                               eristic dimension for
!                                               generality to cover
!                                               parametrizations of
!                                               non-spherical ice.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Removal of variable for
!                                               obsolete solver.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE FLUX_CALC(IERR
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM
     &   , L_CLOUD, L_DROP, L_ICE
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &   , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD, L_SWITCH_SCATTER
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL
!                       Options for Solver
     &   , I_SOLVER
!                       General Spectral Properties
     &   , N_BAND, I_FIRST_BAND, I_LAST_BAND, WEIGHT_BAND
!                       General Atmospheric Properties
     &   , N_PROFILE, N_LAYER
     &   , L_LAYER, L_CLOUD_LAYER
     &   , P, T, T_GROUND, T_LEVEL, D_MASS
!                       Spectral Region
     &   , ISOLIR
!                       Solar Fields
     &   , SEC_0, SOLAR_TOA, SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT
!                       Infra-red Fields
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK
     &   , L_IR_SOURCE_QUAD
!                       Gaseous Absorption
     &   , N_ABSORB, I_GAS_OVERLAP, I_GAS
     &   , GAS_MIX_RATIO, N_BAND_ABSORB, INDEX_ABSORB
     &   , I_BAND_ESFT, W_ESFT, K_ESFT, I_SCALE_ESFT
     &   , I_SCALE_FNC, SCALE_VECTOR
     &   , P_REFERENCE, T_REFERENCE
!                       Doppler Broadening
     &   , L_DOPPLER, DOPPLER_CORRECTION
!                       Surface Fields
     &   , L_SURFACE, I_SURFACE, I_SPEC_SURFACE, SURFACE_ALBEDO
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM
     &   , EMISSIVITY_GROUND, EMISSIVITY_FIELD
!                       Continuum Absorption
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM
!                       Properties of Aerosols
     &   , N_AEROSOL, AEROSOL_MIX_RATIO
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , I_AEROSOL_PARAMETRIZATION, NHUMIDITY, HUMIDITIES
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, W_CLOUD, FRAC_CLOUD
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_DOWN, FLUX_UP
!                       Options for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR
!                       Arrays specific to the UM
!                       Arrays for Coupling
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR
     &   , SEA_FLUX
!                       Arrays for diagnostics specific to the UM
     &   , L_FLUX_BELOW_690NM_SURF, WEIGHT_690NM
     &   , FLUX_BELOW_690NM_SURF
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &   , L_SURF_UP_CLR, SURF_UP_CLR
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , NPD_BAND
     &   , NPD_SPECIES
     &   , NPD_ESFT_TERM, NPD_SCALE_FNC, NPD_SCALE_VARIABLE
     &   , NPD_CONTINUUM
     &   , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES
     &   , NPD_CLOUD_PARAMETER
     &   , NPD_THERMAL_COEFF
     &   , NPD_SURFACE, NPD_ALBEDO_PARM
     &   )
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_BAND
!             NUMBER OF BANDS
     &   , NPD_SPECIES
!             NUMBER OF SPECIES
     &   , NPD_CONTINUUM
!             NUMBER OF CONTINUA
     &   , NPD_AEROSOL_SPECIES
!             NUMBER OF AEROSOL SPECIES
     &   , NPD_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_ESFT_TERM
!             MAXIMUM NUMBER OF ESFT TERMS
     &   , NPD_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             NUMBER OF SCALING VARIABLES
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
     &   , NPD_THERMAL_COEFF
!             MAXIMUM NUMBER OF THERMAL COEFFICIENTS
     &   , NPD_SURFACE
!             NUMBER OF SURFACE TYPES
     &   , NPD_ALBEDO_PARM
!             NUMBER OF PARAMETERS FOR DIRECT ALB.
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
!     MODULE TO SET TREATMENTS OF OVERLAPPING GASAEOUS ABSORPTION.
!
      INTEGER
     &     IP_OVERLAP_SINGLE
!             ONE SPECIES ONLY
     &   , IP_OVERLAP_RANDOM
!             RANDOM OVERLAP
     &   , IP_OVERLAP_FESFT
!             FAST ESFT
     &   , IP_OVERLAP_CLR_FESFT
!             CLEAR-SKY FAST ESFT
     &   , IP_OVERLAP_K_EQV
!             EQUIVALENT EXTINCTION
     &   , IP_OVERLAP_SINGLE_INT
!             INTERPOLATED TREATMENT FOR PRINCIPAL SPECIES
!
      PARAMETER(
     &     IP_OVERLAP_SINGLE=1
     &   , IP_OVERLAP_RANDOM=2
     &   , IP_OVERLAP_FESFT=3
     &   , IP_OVERLAP_CLR_FESFT=4
     &   , IP_OVERLAP_K_EQV=5
     &   , IP_OVERLAP_SINGLE_INT=6
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
!     MODULE TO SET REPRESENTATIONS OF CLOUDS.
!
!     REPRESENTATIONS
      INTEGER
     &     IP_CLOUD_HOMOGEN
!             ALL COMPONENTS ARE MIXED HOMOGENEOUSLY
     &   , IP_CLOUD_ICE_WATER
!             ICE AND WATER CLOUDS ARE TREATED SEPARATELY
     &   , IP_CLOUD_CONV_STRAT
!             CLOUDS ARE DIVIDED INTO HOMOGENEOUSLY MIXED
!             STRATIFORM AND CONVECTIVE PARTS
     &   , IP_CLOUD_CSIW
!             CLOUDS DIVIDED INTO ICE AND WATER PHASES AND
!             INTO STRATIFORM AND CONVECTIVE COMPONENTS.
!
      PARAMETER(
     &     IP_CLOUD_HOMOGEN=1
     &   , IP_CLOUD_ICE_WATER=2
     &   , IP_CLOUD_CONV_STRAT=3
     &   , IP_CLOUD_CSIW=4
     &   )
!
!     TYPES OF CLOUDS:
      INTEGER
     &     NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)
!             NUMBER OF TYPE OF CLOUDS IN REPRESENTATION
!
      INTEGER
     &     IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT
     &       , NPD_CLOUD_REPRESENTATION)
!            MAP OF COMPONENTS CONTRIBUTING TO TYPES OF CLOUDS
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET COMPONENTS OF CLOUDS.
!
      INTEGER
     &     IP_CLCMP_ST_WATER
!             STRATIFORM WATER DROPLETS
     &   , IP_CLCMP_ST_ICE
!             STRATIFORM ICE CRYSTALS
     &   , IP_CLCMP_CNV_WATER
!             CONVECTIVE WATER DROPLETS
     &   , IP_CLCMP_CNV_ICE
!             CONVECTIVE ICE CRYSTALS
!
      PARAMETER(
     &     IP_CLCMP_ST_WATER=1
     &   , IP_CLCMP_ST_ICE=2
     &   , IP_CLCMP_CNV_WATER=3
     &   , IP_CLCMP_CNV_ICE=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET TYPES OF CLOUDS.
!
      INTEGER
     &     IP_CLOUD_TYPE_HOMOGEN
!             CLOUD COMPOSED OF MIXED WATER AND ICE
     &   , IP_CLOUD_TYPE_WATER
!             CLOUD COMPOSED ONLY OF WATER
     &   , IP_CLOUD_TYPE_ICE
!             CLOUD COMPOSED ONLY OF ICE
     &   , IP_CLOUD_TYPE_STRAT
!             MIXED-PHASE STRATIFORM CLOUD
     &   , IP_CLOUD_TYPE_CONV
!             MIXED-PHASE CONVECTIVE CLOUD
     &   , IP_CLOUD_TYPE_SW
!             STRATIFORM WATER CLOUD
     &   , IP_CLOUD_TYPE_SI
!             STRATIFORM ICE CLOUD
     &   , IP_CLOUD_TYPE_CW
!             CONVECTIVE WATER CLOUD
     &   , IP_CLOUD_TYPE_CI
!             CONVECTIVE ICE CLOUD
!
      PARAMETER(
     &     IP_CLOUD_TYPE_HOMOGEN=1
     &   , IP_CLOUD_TYPE_WATER=1
     &   , IP_CLOUD_TYPE_ICE=2
     &   , IP_CLOUD_TYPE_STRAT=1
     &   , IP_CLOUD_TYPE_CONV=2
     &   , IP_CLOUD_TYPE_SW=1
     &   , IP_CLOUD_TYPE_SI=2
     &   , IP_CLOUD_TYPE_CW=3
     &   , IP_CLOUD_TYPE_CI=4
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
!     MODULE TO DEFINE REFERENCE NUMBERS FOR SOLVERS.
!
      INTEGER
     &     IP_SOLVER_PENTADIAGONAL
!             PENTADIAGONAL SCHEME
     &   , IP_SOLVER_MIX_11
!             MIXED COLUMN SCHEME USING FULL ENDEKADIAGONAL MATRIX
     &   , IP_SOLVER_MIX_APP_SCAT
!             MIXED COLUMN SCHEME WITH APPROXIMATE SCATTERING
     &   , IP_SOLVER_MIX_DIRECT
!             DIRECT MIXED COLUMN SCHEME FOR FULL FLUXES 
     &   , IP_SOLVER_HOMOGEN_DIRECT
!             DIRECT SOLVER FOR A HOMOGENEOUS COLUMN
     &   , IP_SOLVER_TRIPLE
!             DIRECT SOLVER FOR TRIPLE COLUMN
     &   , IP_SOLVER_TRIPLE_APP_SCAT
!             DIRECT SOLVER FOR TRIPLE COLUMN APPROXIMATING SCATTERING
!
      PARAMETER(
     &     IP_SOLVER_PENTADIAGONAL=1
     &   , IP_SOLVER_MIX_11=6
     &   , IP_SOLVER_MIX_APP_SCAT=9
     &   , IP_SOLVER_MIX_DIRECT=11
     &   , IP_SOLVER_HOMOGEN_DIRECT=13
     &   , IP_SOLVER_TRIPLE=14
     &   , IP_SOLVER_TRIPLE_APP_SCAT=15
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET THE DEFINING NUMBERS FOR THE TWO-STREAM
!     SCHEMES.
!
      INTEGER
     &     IP_EDDINGTON
!             EDDINGTON APPROXIMATION
     &   , IP_DISCRETE_ORD
!             DISCRETE ORDINATE METHOD
     &   , IP_IFM
!             IMPROVED FLUX METHOD
     &   , IP_PIFM85
!             PRACTICAL IMPROVED FLUX METHOD
!             (VERSION OF ZDUNKOWSKI ET AL. 1985)
     &   , IP_ZDK_FLUX
!             ZDUNKOWSKI'S FLUX METHOD
     &   , IP_KRSCHG_FLUX
!             KERSCHGEN'S FLUX METHOD
     &   , IP_COAKLEY_CHYLEK_1
!             COAKLEY & CHYLEK'S 1ST METHOD
     &   , IP_COAKLEY_CHYLEK_2
!             COAKLEY & CHYLEK'S 2ND METHOD
     &   , IP_MEADOR_WEAVER
!             MEADOR & WEAVER'S METHOD
     &   , IP_ELSASSER
!             ELSASSER'S DIFFUSIVITY SCHEME
     &   , IP_2S_TEST
!             USER'S DEFINED TEST APPROXIMATION.
     &   , IP_HEMI_MEAN
!             HEMISPHERIC MEAN APPROXIMATION.
     &   , IP_PIFM80
!             PRACTICAL IMPROVED FLUX METHOD
!             (VERSION OF ZDUNKOWSKI ET AL. 1980)
!
      PARAMETER(
     &     IP_EDDINGTON=2
     &   , IP_DISCRETE_ORD=4
     &   , IP_IFM=5
     &   , IP_PIFM85=6
     &   , IP_ZDK_FLUX=7
     &   , IP_KRSCHG_FLUX=8
     &   , IP_COAKLEY_CHYLEK_1=9
     &   , IP_COAKLEY_CHYLEK_2=10
     &   , IP_MEADOR_WEAVER=11
     &   , IP_ELSASSER=12
     &   , IP_2S_TEST=14
     &   , IP_HEMI_MEAN=15
     &   , IP_PIFM80=16
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
!     MODULE TO SET THE AEROSOL PARAMETRIZATIONS.
!
      INTEGER
     &     IP_AEROSOL_PARAM_DRY
!             PARAMETRIZATION FOR DRY AEROSOLS
     &   , IP_AEROSOL_PARAM_MOIST
!             PARAMETRIZATION FOR MOIST AEROSOLS
     &   , IP_AEROSOL_UNPARAMETRIZED
!             OBSERVATIONAL AEROSOL DATA
!
      PARAMETER(
     &     IP_AEROSOL_PARAM_DRY=1
     &   , IP_AEROSOL_PARAM_MOIST=2
     &   , IP_AEROSOL_UNPARAMETRIZED=3
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
!STR  GENERAL LOGICAL SWITCHES:
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             VALUES GIVEN IN LAYERS
     &   , L_CLOUD_LAYER
!             CLOUD VALUES GIVEN IN LAYERS
     &   , L_CLEAR
!             CALCULATE CLEAR-SKY FLUXES
     &   , L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE FUNCTION
     &   , L_RESCALE
!             FLAG FOR DELTA-RESCALING
     &   , L_2_STREAM_CORRECT
!             CORRECTION TO 2-STREAM SCHEME
!
!STR  PARAMETERS CONTROLLING ALGORITHMS:
!     REPRESENTATION OF CLOUDS:
      INTEGER   !, INTENT(IN)
     &     I_CLOUD
!             CLOUD SCHEME USED
      LOGICAL   !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP_GLOBAL
!             GLOBAL TOPMOST CLOUDY LAYER
!
!     NUMERICAL ALGORITHMS:
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             VISIBLE OR IR
     &   , I_SOLVER
!             SOLVER USED
     &   , I_SOLVER_CLEAR
!             CLEAR SOLVER USED
     &   , I_2STREAM
!             TWO-STREAM SCHEME
     &   , I_ANGULAR_INTEGRATION
!             ANGULAR INTEGRATION SCHEME
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
!     RANGE OF SPECTRAL BANDS:
      INTEGER   !, INTENT(IN)
     &     I_FIRST_BAND
!             FIRST BAND
     &   , I_LAST_BAND
!             LAST BAND
!
!     GENERAL PROPERTIES OF SPECTRUM:
      INTEGER   !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , N_ABSORB
!             NUMBER OF ABSORBERS
     &   , N_AEROSOL
!             NUMBER OF AEROSOL SPECIES
!
!STR  SOLAR FIELDS:
      REAL      !, INTENT(IN)
     &     SOLAR_TOA(NPD_PROFILE)
!             INCIDENT SOLAR RADIATION
     &   , SOLAR_FLUX_BAND(NPD_BAND)
!             NORMALIZED FLUX IN BAND
     &   , SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!STR  ATMOSPHERIC PROFILES:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
     &   , T_GROUND(NPD_PROFILE)
!             TEMPERATURE OF GROUND
     &   , T_LEVEL(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE ON LEVELS
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             GASEOUS MASS MIXING RATIOS
!
!STR  SURFACE PROPERTIES:
      LOGICAL   !, INTENT(IN)
     &     L_SURFACE(NPD_SURFACE)
!             TYPES OF SURFACES FOR WHICH DATA ARE PRESENT
      INTEGER   !, INTENT(IN)
     &     I_SURFACE(NPD_PROFILE)
!             TYPE OF SURFACE AT THE FOOT OF EACH PROFILE
     &   , I_SPEC_SURFACE(NPD_SURFACE)
!             METHOD OF SPECIFYING ALBEDO
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS IN FIT TO DIRECT ALBEDO
      REAL      !, INTENT(IN)
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDO
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND)
!             SPECIFIED DIFFUSE ALBEDOS
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND)
!             SPECIFIED DIRECT ALBEDOS
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR DIRECT ALBEDOS
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
     &   , EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND)
!             SPECIFIED EMISSIVITIES
!
!STR  RAYLEIGH SCATTERING:
      LOGICAL   !, INTENT(IN)
     &     L_RAYLEIGH
!             INCLUDE RAYLEIGH SCATTERING IN THE CALCULATION.
      REAL      !, INTENT(IN)
     &     RAYLEIGH_COEFFICIENT(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!STR  FIELDS FOR GASEOUS ABSORPTION:
      LOGICAL   !, INTENT(IN)
     &     L_GAS
!             INCLUDE GAS ABSORPTION IN THE CALCULATION
!     GASEOUS OVERLAPS:
      INTEGER   !, INTENT(IN)
     &     I_GAS_OVERLAP(NPD_BAND)
!             GAS OVERLAP ASSUMPTION
     &   , I_GAS
!             GAS TO BE CONSIDERED (ONE GAS ONLY)
!     ESFTS:
      INTEGER   !, INTENT(IN)
     &     N_BAND_ABSORB(NPD_BAND)
!             NUMBER OF ABSORBERS IN BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN BANDS
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)
!             NUMBER OF TERMS IN BAND
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
      REAL      !, INTENT(IN)
     &     W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             WEIGHTS FOR ESFT
     &   , K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             EXPONENTIAL ESFT TERMS
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND
     &        , NPD_SPECIES)
!             ABSORBER SCALING PARAMETERS
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE SCALING PRESSURE
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE SCALING TEMPERATURE
!
!STR  SPECTRAL DATA FOR THE CONTINUUM:
      LOGICAL   !, INTENT(IN)
     &     L_CONTINUUM
!             INCLUDE CONTINUUM ABSORPTION IN THE CALCULATION
      INTEGER   !, INTENT(IN)
     &     N_BAND_CONTINUUM(NPD_BAND)
!             NUMBER OF CONTINUA IN BANDS
     &   , INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             INDICES OF CONTINUA
     &   , INDEX_WATER
!             INDEX OF WATER
     &   , I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
      REAL      !, INTENT(IN)
     &     K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             CONTINUUM EXTINCTION COEFFICIENTS
     &   , SCALE_CONTINUUM(NPD_SCALE_VARIABLE, NPD_BAND, NPD_CONTINUUM)
!             CONTINUUM SCALING PARAMETERS
     &   , P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             CONTINUUM REFERENCE PRESSURE
     &   , T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             CONTINUUM REFERENCE TEMPERATURE
!
!
!STR  GENERAL CLOUD FIELDS:
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD
!             CLOUDS ARE REQUIRED IN THE CALCULATION
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             AMOUNT OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
!
!STR  FIELDS FOR MICROPHYSICAL QUANTITIES:
!
      LOGICAL   !, INTENT(IN)
     &     L_DROP
!             INCLUDE DROPLETS IN THE CALCULATION
     &   , L_ICE
!             INCLUDE ICE IN THE CALCULATION
      INTEGER   !, INTENT(IN)
     &     N_CONDENSED
!             NUMBER OF CONDENSED COMPONENTS IN CLOUDS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF CONDENSED COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATION SCHEMES FOR COMPONENTS
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF MIXING RULE CHOSEN
!
      REAL      !, INTENT(IN)
     &     CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             MIXING RATIOS OF CONDENSED COMPONENTS
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CONDENSED COMPONENTS
     &   , CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_CLOUD_COMPONENT, NPD_BAND)
!             COEFFICIENTS IN PARAMETRIZATIONS OF CONDENSED PHASES
!
!
!
!STR  FIELDS FOR AEROSOLS:
      LOGICAL   !, INTENT(IN)
     &     L_AEROSOL
!             INCLUDE AEROSOLS IN THE CALCULATION
      INTEGER   !, INTENT(IN)
     &     I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION FLAGS FOR AEROSOL
      INTEGER   !, INTENT(IN)
     &     NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBER OF HUMIDITIES
      REAL      !, INTENT(IN)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_AEROSOL_SPECIES)
!             NUMBER DENSITY OF AEROSOLS
      REAL      !, INTENT(IN)
     &     AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ASYMMETRY BY AEROSOLS
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR SPECIES
!
!
!STR  FITTING OF THE PLANCKIAN FUNCTION:
      INTEGER   !, INTENT(IN)
     &     N_DEG_FIT
!             DEGREE OF THERMAL FITTING FNC.
      REAL      !, INTENT(IN)
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)
!             COEFFICIENTS OF SOURCE TERMS
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
!
!STR  DOPPLER BROADENING
      LOGICAL   !, INTENT(IN)
     &     L_DOPPLER(NPD_SPECIES)
!             FLAGS TO ACTIVATE DOPPLER CORRECTIONS
      REAL      !, INTENT(IN)
     &     DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER BROADENING TERM
      REAL      !, INTENT(OUT)
     &     WEIGHT_BAND(NPD_BAND)
!             WEIGHTING FUNCTION FOR BANDS
!
!STR  CONTROL OF SCATTERING:
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD
!             METHOD OF TREATING SCATTERING
      LOGICAL   !, INTENT(IN)
     &     L_SWITCH_SCATTER(NPD_BAND)
!             SWITCHES FOR SCATTERING IN BANDS
!
!
!STR  FLUXES CALCULATED:
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_DOWN_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DOWNWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR UPWARD FLUX
!
!STR  ARRAYS SPECIFIC TO THE UNIFIED MODEL
!
!     SWITCHES FOR DIAGNOSTICS:
      LOGICAL   !, INTENT(IN)
     &     L_FLUX_BELOW_690NM_SURF
!             FLUX BELOW 690NM AT SURFACE TO BE CALCULATED
     &   , L_SURFACE_DOWN_FLUX
!             DOWNWARD SURFACE FLUX REQUIRED
     &   , L_SURF_DOWN_CLR
!             CALCULATE DOWNWARD CLEAR FLUX
     &   , L_SURF_UP_CLR
!             CALCULATE UPWARD CLEAR FLUX
!
!     ARRAYS FOR USE WITH COUPLING:
      INTEGER   !, INTENT(IN)
     &     N_FRAC_ICE_POINT
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER
      REAL      !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!             ICE FRACTION
      REAL      !, INTENT(IN)
     &     ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND)
!             DIFFUSE ALBEDO FOR OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND)
!             DIRECT ALBEDO FOR OPEN SEA
!
!     ARRAYS FOR USE WITH DIAGNOSTICS:
      REAL      !, INTENT(IN)
     &     WEIGHT_690NM(NPD_BAND)
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL      !, INTENT(OUT)
     &     SEA_FLUX(NPD_PROFILE)
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX(NPD_PROFILE)
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR(NPD_PROFILE)
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR(NPD_PROFILE)
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF(NPD_PROFILE)
!             SURFACE FLUX BELOW 690NM
!
!
!
!     LOCAL ARGUMENTS.
!     GENERAL POINTERS:
      INTEGER
     &     I_TOP
!             TOP LEVEL OF PROFILES
     &   , I_BAND
!             SPECTRAL BAND
     &   , N_AUGMENT
!             LENGTH OF LONG FLUX VECTOR
     &   , N_GAS
!             NUMBER OF ACTIVE GASES
     &   , I_GAS_BAND
!             SINGLE VARIABLE FOR GAS IN BAND
     &   , N_CONTINUUM
!             NUMBER OF CONTINUA IN BAND
     &   , I_CONTINUUM
!             CONTINUUM NUMBER
     &   , I_CONTINUUM_POINTER(NPD_CONTINUUM)
!             POINTERS TO CONTINUA
!
!     SCATTERING IN INDIVIDUAL BANDS
      INTEGER
     &     I_SCATTER_METHOD_BAND
!             METHOD OF TREATING SCATTERING IN GIVEN BAND
!
!     VARIABLES FOR SURFACE PROPERTIES:
      INTEGER
     &     N_POINT_TYPE(NPD_SURFACE)
!             NUMBER OF POINTS OF EACH TYPE
     &   , INDEX_SURFACE(NPD_PROFILE, NPD_SURFACE)
!             INDICES OF EACH SURFACE TYPE
!
!     POINTERS TO THE CONTENTS OF LAYERS:
      INTEGER
     &     N_FREE_PROFILE(NPD_LAYER)
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)
!             COLUMNS CONTAINING FREE PROFILES
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
     &   , N_CLOUD_PROFILE(NPD_LAYER)
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
!
!     POINTERS TO TYPES OF CLOUDS:
      LOGICAL
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             LOGICAL SWITCHES TO INCLUDE COMPONENTS
      INTEGER
     &     I_PHASE_CMP(NPD_CLOUD_COMPONENT)
!             PHASES OF COMPONENTS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH PARTICULAR TYPE OF CLOUD FALL
!
!     FRACTIONAL COVERAGE OF DIFFERENT REGIONS:
      REAL
     &     FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTION OF TOTAL CLOUD OCCUPIED BY SPECIFIC REGIONS
!
!     POINTER TO TABLE OF HUMIDITY:
      INTEGER
     &     I_HUMIDITY_POINTER(NPD_PROFILE, NPD_LAYER)
!             POINTER TO LOOK-UP TABLE FOR AEROSOLS
!
!     CONTROLLING VARIABLES:
      INTEGER
     &     I
!             LOOP INDEX
     &   , J
!             LOOP INDEX
     &   , K
!             LOOP INDEX
     &   , L
!             LOOP INDEX
!
!     LOGICAL SWITCHES:
      LOGICAL
     &     L_GAS_BAND
!             FLAG TO INCLUDE GASEOUS ABSORPTION IN A PARTICULAR BAND
     &   , L_MOIST_AEROSOL
!             FLAG FOR MOIST AEROSOL
     &   , L_AEROSOL_DENSITY
!             FLAG FOR CALCULATION OF ATMOSPHERIC DENSITY FOR AEROSOLS
     &   , L_NET
!             FLAG FOR NET FLUXES
!
!     SURFACE PROPERTIES:
      REAL
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
!
      REAL
     &     INC_SOLAR_FLUX_BAND(NPD_PROFILE)
!             INCIDENT SOLAR FLUX IN BAND
      REAL
     &     GAS_FRAC_RESCALED(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             RESCALED GAS MIXING RATIOS
     &   , AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER, NPD_CONTINUUM)
!             AMOUNTS OF CONTINUA
     &   , K_CONTINUUM_MONO(NPD_CONTINUUM)
!             MONOCHROMATIC CONTINUUM COMPONENTS
!
!     THERMAL ARRAYS:
      REAL
     &     PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCK FUNCTION IN BAND AT LEVELS
     &   , DIFF_PLANCK_BAND(NPD_PROFILE, NPD_LAYER)
!             DIFFERENTIAL THERMAL SOURCE IN BAND
     &   , DIFF_PLANCK_BAND_2(NPD_PROFILE, NPD_LAYER)
!             2 x 2ND DIFF. THERMAL SOURCE IN BAND
     &   , THERMAL_GROUND_BAND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION IN BAND
!
!     ATMOSPHERIC DENSITIES:
      REAL
     &     DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             OVERALL DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF WATER
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF FOREIGN SPECIES
!
!     FIELDS FOR MOIST AEROSOLS:
      REAL
     &     DELTA_HUMIDITY
!             INCREMENT IN LOOK-UP TABLE FOR HUM.
     &   , MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)
!             MEAN RELATIVE HUMIDITY OF LAYERS
!
!     FUNDAMENTAL OPTICAL PROPERTIES OF LAYERS:
      REAL
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE TOTAL GREY EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE SCATTERING EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ASYMMETRIES
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE FORWARD SCATTERING
     &   , K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             TOTAL CLOUDY GREY EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SCATTERING EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
!     LOCAL RADIATIVE FLUXES:
      REAL
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
     &   , FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX IN BAND
     &   , FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX IN BAND
     &   , PLANCK_FLUX(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN FLUX IN BAND
!
!     COEFFICIENTS FOR THE TRANSFER OF ENERGY BETWEEN
!     PARTIALLY CLOUDY LAYERS:
      REAL
     &     CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS DEFINING OVERLAPPING OPTIONS FOR CLOUDS:
!             THESE ALSO DEPEND ON THE SOLVER SELECTED.
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY FRACTION
      INTEGER
     &     N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS REQUIRED
      LOGICAL
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             COLUMN FLAGS FOR COLUMNS
      REAL
     &     AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREAS OF COLUMNS
!
!     LOCAL VARIABLES SPECIFIC TO THE UNIFIED MODEL.
      REAL
     &     PLANCK_FREEZE_SEA
!             PLANCK FUNCTION OVER FREEZING SEA
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     SET_CLOUD_POINTER, SET_CLOUD_GEOMETRY, SET_SCATTERING
     &   , OVERLAP_MIX_MAXIMUM, OVERLAP_MIX_RANDOM
     &   , SPLIT_MAXIMUM, COLLECT_SURFACE, INITIALIZE_FLUX
     &   , CALCULATE_DENSITY, SET_MOIST_AEROSOL_PROPERTIES
     &   , SCALE_ABSORB, RESCALE_CONTINUUM, GREY_EXTINCTION
     &   , RESCALE_ASYMMETRY
     &   , DIFF_PLANCK_SOURCE, SET_SURFACE_PROPERTIES
     &   , SOLVE_BAND_WITHOUT_GAS, SOLVE_BAND_ONE_GAS
     &   , SOLVE_BAND_RANDOM_OVERLAP, SOLVE_BAND_FESFT
     &   , SOLVE_BAND_CLR_FESFT, SOLVE_BAND_K_EQV
     &   , AUGMENT_TOTAL_FLUX, ASSIGN_FLUX
     &   , R2_INIT_COUPLE_DIAG, R2_COUPLE_DIAG
     &   , OVERLAP_TRIPLE
!     FUNCTIONS CALLED:
      LOGICAL
     &     L_CLOUD_DENSITY
!             FLAG FOR CALCULATION OF ATMOSPHERIC DENSITIES FOR CLOUDS
      EXTERNAL
     &     L_CLOUD_DENSITY
!
!
!     SETTING OF PROPERTIES OF ARRAYS:
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET PROPERTIES OF REPRESENTATIONS OF CLOUDS.
!
      DATA
     &    NP_CLOUD_TYPE(IP_CLOUD_HOMOGEN)/1/
     &  , NP_CLOUD_TYPE(IP_CLOUD_ICE_WATER)/2/
     &  , NP_CLOUD_TYPE(IP_CLOUD_CONV_STRAT)/2/
     &  , NP_CLOUD_TYPE(IP_CLOUD_CSIW)/4/
!
!
!     THE ARRAY IP_CLOUD_TYPE_MAP INDICATES TO WHICH TYPE OF CLOUD
!     EACH COMPONENT BELONGS IN A PARTICULAR REPRESENTATION. AN
!     ENTRY OF 0 INDICATES THAT THAT COMPONENT SHOULD NOT BE
!     PRESENT IN THE REPRESENTATION.
!
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_HOMOGEN)
     &       /IP_CLOUD_TYPE_HOMOGEN/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_HOMOGEN)
     &       /IP_CLOUD_TYPE_HOMOGEN/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_HOMOGEN)
     &       /0/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_HOMOGEN)
     &       /0/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_ICE_WATER)
     &       /IP_CLOUD_TYPE_WATER/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_ICE_WATER)
     &       /IP_CLOUD_TYPE_ICE/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_ICE_WATER)
     &       /0/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_ICE_WATER)
     &       /0/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_STRAT/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_CONV/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_STRAT/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_CONV/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_SW/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_SI/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_CW/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_CI/
!
!     ------------------------------------------------------------------
!
!
!
!
!     INITIAL DETERMINATION OF FLAGS AND SWITCHES:
!
      IF (I_ANGULAR_INTEGRATION.EQ.IP_TWO_STREAM) THEN
!
!        CURRENTLY, WE DO NOT ALLOW THE USE OF SOLVERS FOR THE
!        NET FLUX ALONE AS THIS MAKES IT DIFFICULT TO INCORPORATE
!        SOME DIAGNOSTICS.
!
         L_NET=.FALSE.
!
      ELSE IF (I_ANGULAR_INTEGRATION.EQ.IP_IR_GAUSS) THEN
!
         L_NET=.FALSE.
!
      ENDIF
!
!     THE LENGTH OF THE LONG FLUX VECTOR DEPENDS ON WHETHER WE SOLVE
!     FOR THE FULL OR THE NET FLUX.
      IF (L_NET) THEN
         N_AUGMENT=N_LAYER+1
      ELSE
         N_AUGMENT=2*(N_LAYER+1)
      ENDIF
!
!     SET THE TOP LEVEL OF THE PROFILES.
      IF (L_LAYER) THEN
         I_TOP=1
      ELSE
         I_TOP=0
      ENDIF
!
!
!
!     INITIAL CALCULATIONS FOR SURFACE PROPERTIES:
!
!     COLLECT POINTS WITH THE SAME SURFACE SPECIFICATION.
      CALL COLLECT_SURFACE(N_PROFILE
     &   , I_SURFACE
     &   , N_POINT_TYPE, INDEX_SURFACE
     &   , NPD_PROFILE, NPD_SURFACE
     &   )
!
!
!
!     INITIAL CALCULATIONS FOR AEROSOLS:
!
!     SET THE SPECTRALLY INDEPENDENT PROPERTIES OF MOIST AEROSOLS.
      L_MOIST_AEROSOL=.FALSE.
      DO J=1, N_AEROSOL
         L_MOIST_AEROSOL=L_MOIST_AEROSOL.OR.
     &      (I_AEROSOL_PARAMETRIZATION(J)
     &      .EQ.IP_AEROSOL_PARAM_MOIST)
      ENDDO
!
      IF (L_MOIST_AEROSOL) THEN
         CALL SET_MOIST_AEROSOL_PROPERTIES(IERR
     &      , N_PROFILE, N_LAYER, L_LAYER
     &      , N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY
     &      , GAS_MIX_RATIO(1, 0, INDEX_WATER), T, P, DELTA_HUMIDITY
     &      , MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER
     &      , NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
      ENDIF
!
!
!     CHECK WHETHER THE DENSITIES WILL BE NEEDED FOR
!     UNPARAMETRIZED AEROSOLS.
      L_AEROSOL_DENSITY=.FALSE.
      IF (L_AEROSOL) THEN
         DO J=1, N_AEROSOL
            L_AEROSOL_DENSITY=L_AEROSOL_DENSITY.OR.
     &         (I_AEROSOL_PARAMETRIZATION(J).EQ.
     &            IP_AEROSOL_PARAM_MOIST)
     &         .OR.(I_AEROSOL_PARAMETRIZATION(J).EQ.
     &            IP_AEROSOL_UNPARAMETRIZED)
         ENDDO
      ENDIF
!
!
!
!     INITIAL CALCULATIONS FOR CLOUDS:
!
!     SET POINTERS TO THE TYPES OF CLOUD.
      CALL SET_CLOUD_POINTER(IERR
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION
     &   , L_DROP, L_ICE
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!     SET THE GEOMETRY OF THE CLOUDS.
      CALL SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL
     &   , W_CLOUD
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &   , N_CLOUD_TOP
     &   , N_FREE_PROFILE, I_FREE_PROFILE
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
      IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!        AGGREGATE CLOUDS INTO REGIONS FOR SOLVING.
         CALL AGGREGATE_CLOUD(IERR
     &      , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , I_CLOUD, I_CLOUD_REPRESENTATION
     &      , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)
     &      , FRAC_CLOUD
     &      , I_REGION_CLOUD, FRAC_REGION
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!     CALCULATE ENERGY TRANSFER COEFFICIENTS IN A MIXED COLUMN,
!     OR SPLIT THE ATMOSPHERE INTO COLUMNS WITH A COLUMN MODEL:
!
      IF (I_CLOUD.EQ.IP_CLOUD_MIX_MAX) THEN
!
         CALL OVERLAP_MIX_MAXIMUM(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , ISOLIR, I_SOLVER
     &      , W_CLOUD, W_FREE
     &      , CLOUD_OVERLAP
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) THEN
         CALL OVERLAP_MIX_RANDOM(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , ISOLIR, I_SOLVER
     &      , W_CLOUD, W_FREE
     &      , CLOUD_OVERLAP
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
!        CALCULATE OVERLAPS FOR THE TRIPLE DECOMPOSITION
!        OF THE COLUMN INTO STRATIFORM AND CONVECTIVE PARTS.
         CALL OVERLAP_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , W_CLOUD, W_FREE, FRAC_REGION
     &      , CLOUD_OVERLAP
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
!
      ELSE IF (I_CLOUD.EQ.IP_CLOUD_COLUMN_MAX) THEN
!
         CALL SPLIT_MAXIMUM(N_PROFILE, N_LAYER
     &      , W_CLOUD
     &      , N_COLUMN, AREA_COLUMN, L_COLUMN
     &      , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &      )
      ENDIF
!
!
!     CALCULATE THE ATMOSPHERIC DENSITIES:
!
      IF ( L_CONTINUUM
     &      .OR.L_AEROSOL_DENSITY
     &      .OR.(L_CLOUD
     &      .AND.L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP
     &                        , I_CONDENSED_PARAM
     &                        ) ) ) THEN
         CALL CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM
     &      , GAS_MIX_RATIO(1, 0, INDEX_WATER)
     &      , P, T, I_TOP
     &      , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!
!     INITIALIZE THE TOTAL FLUXES.
!
      CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &   , ISOLIR
     &   , FLUX_DIRECT, FLUX_TOTAL
     &   , L_CLEAR
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &   , 0.0E+00
     &   , NPD_PROFILE, NPD_LAYER
     &   , L_NET
     &   )
!
!     INITIALIZE THE PLANCKIAN FLUXES IF REQUIRED.
      IF ( (ISOLIR.EQ.IP_INFRA_RED).AND.(.NOT.L_NET) ) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               PLANCK_FLUX(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
!
!     INITIALIZATION OF DIAGNOSTICS AND COUPLING ARRAYS FOR THE
!     UNIFIED MODEL.
      CALL R2_INIT_COUPLE_DIAG(N_PROFILE
     &   , SEA_FLUX
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &   , L_SURF_UP_CLR, SURF_UP_CLR
     &   , L_FLUX_BELOW_690NM_SURF, FLUX_BELOW_690NM_SURF
     &   , NPD_PROFILE
     &   )
!
!
!
!
!
!     SOLVE THE EQUATION OF TRANSFER IN EACH BAND AND
!     INCREMENT THE FLUXES.
!
      DO I_BAND=I_FIRST_BAND, I_LAST_BAND
!
!
!        SET THE FLAG FOR THE TREATMENT OF SCATTERING IN THIS BAND
         CALL SET_SCATTERING(I_SCATTER_METHOD
     &      , L_SWITCH_SCATTER(I_BAND)
     &      , I_SCATTER_METHOD_BAND
     &      )
!
!        DETERMINE WHETHER GASEOUS ABSORPTION IS INCLUDED IN THIS BAND.
         IF ( (L_GAS).AND.(N_BAND_ABSORB(I_BAND).GT.0) ) THEN
!
!           NOTE: I_GAS_BAND IS USED EXTENSIVELY BELOW SINCE NESTED
!           ARRAY ELEMENTS IN A SUBROUTINE CALL (SEE LATER) CAN
!           CONFUSE SOME COMPILERS.
!
!           NORMALLY THE NUMBER OF GASES IN THE CALCULATION WILL BE
!           AS IN THE SPECTRAL FILE, BUT PARTICULAR OPTIONS MAY RESULT
!           IN THE OMISSION OF SOME GASES.
!
            N_GAS=N_BAND_ABSORB(I_BAND)
!
            IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_SINGLE) THEN
!
!              THERE WILL BE NO GASEOUS ABSORPTION IN THIS BAND
!              UNLESS THE SELECTED GAS APPEARS.
               N_GAS=0
!
               DO I=1, N_BAND_ABSORB(I_BAND)
                  IF (INDEX_ABSORB(I, I_BAND).EQ.I_GAS) N_GAS=1
               ENDDO
!
            ENDIF
!
!
            IF (N_GAS.GT.0) THEN
!
!              SET THE FLAG FOR GASEOUS ABSORPTION IN THE BAND.
               L_GAS_BAND=.TRUE.
!
               DO J=1, N_GAS
!
                  I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
!
!                 RESET THE POINTER IF THERE IS JUST ONE GAS.
!
                  IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_SINGLE)
     &               THEN
!                    ONLY THE SELECTED GAS IS ACTIVE IN THE BAND.
                     I_GAS_BAND=I_GAS
!
                  ENDIF
!
                  IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)
     &               .EQ.IP_SCALE_BAND) THEN
!                    RESCALE THE AMOUNT OF GAS FOR THIS BAND NOW.
                     CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &                  , GAS_MIX_RATIO(1, 0, I_GAS_BAND), P, T
     &                  , L_LAYER, I_TOP
     &                  , GAS_FRAC_RESCALED(1, 0, I_GAS_BAND)
     &                  , I_SCALE_FNC(I_BAND, I_GAS_BAND)
     &                  , P_REFERENCE(I_GAS_BAND, I_BAND)
     &                  , T_REFERENCE(I_GAS_BAND, I_BAND)
     &                  , SCALE_VECTOR(1, 1, I_BAND, I_GAS_BAND)
     &                  , L_DOPPLER(I_GAS_BAND)
     &                  , DOPPLER_CORRECTION(I_GAS_BAND)
     &                  , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &                  , NPD_SCALE_VARIABLE
     &                  )
                     IF (IERR.NE.I_NORMAL) RETURN
!
                  ELSE IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)
     &               .EQ.IP_SCALE_NULL) THEN
!                    COPY ACROSS THE UNSCALED ARRAY.
                        DO I=1, N_LAYER
                           DO L=1, N_PROFILE
                              GAS_FRAC_RESCALED(L, I, I_GAS_BAND)
     &                           =GAS_MIX_RATIO(L, I, I_GAS_BAND)
                           ENDDO
                        ENDDO
                  ENDIF
               ENDDO
            ELSE
               L_GAS_BAND=.FALSE.
            ENDIF
!
         ELSE
            L_GAS_BAND=.FALSE.
         ENDIF
!
!
!
!        RESCALE AMOUNTS OF CONTINUA.
!
         IF (L_CONTINUUM) THEN
            N_CONTINUUM=N_BAND_CONTINUUM(I_BAND)
            DO I=1, N_CONTINUUM
               I_CONTINUUM_POINTER(I)=INDEX_CONTINUUM(I_BAND, I)
               I_CONTINUUM=I_CONTINUUM_POINTER(I)
               K_CONTINUUM_MONO(I_CONTINUUM)
     &            =K_CONTINUUM(I_BAND, I_CONTINUUM)
               CALL RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM
     &            , P, T, L_LAYER, I_TOP
     &            , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN
     &            , GAS_MIX_RATIO(1, 0, INDEX_WATER)
     &            , AMOUNT_CONTINUUM(1, 0, I_CONTINUUM)
     &            , I_SCALE_FNC_CONT(I_BAND, I_CONTINUUM)
     &            , P_REF_CONTINUUM(I_CONTINUUM, I_BAND)
     &            , T_REF_CONTINUUM(I_CONTINUUM, I_BAND)
     &            , SCALE_CONTINUUM(1, I_BAND, I_CONTINUUM)
     &            , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &            , NPD_SCALE_VARIABLE
     &            )
            ENDDO
         ENDIF
!
!
!
!        CALCULATE THE GREY EXTINCTION WITHIN THE BAND.
!
         CALL GREY_EXTINCTION(IERR
     &      , N_PROFILE, N_LAYER, L_LAYER, P, T, DENSITY
     &      , L_RESCALE
     &      , L_RAYLEIGH, RAYLEIGH_COEFFICIENT(I_BAND)
     &      , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER
     &      , K_CONTINUUM_MONO, AMOUNT_CONTINUUM
     &      , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO
     &      , I_AEROSOL_PARAMETRIZATION
     &      , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY
     &      , MEAN_REL_HUMIDITY
     &      , AEROSOL_ABSORPTION(1, 1, I_BAND)
     &      , AEROSOL_SCATTERING(1, 1, I_BAND)
     &      , AEROSOL_ASYMMETRY(1, 1, I_BAND)
     &      , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP
     &      , L_CLOUD_LAYER, I_CLOUD
     &      , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP
     &      , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST(1, 1, I_BAND)
     &      , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &      , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)
     &      , I_CLOUD_TYPE
     &      , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &      , FORWARD_SCATTER_FREE
     &      , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &      , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
     &      , NPD_PROFILE, NPD_LAYER, NPD_CONTINUUM
     &      , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES
     &      , NPD_CLOUD_PARAMETER
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
!
!
!
         IF (I_ANGULAR_INTEGRATION.EQ.IP_TWO_STREAM) THEN
!
!           RESCALE THE ASYMMETRY AND CALCULATE THE SCATTERING
!           FRACTIONS. (THESE ARE GREY AND MAY BE CALCULATED OUTSIDE
!           A LOOP OVER GASES).
!
            IF (L_RESCALE) THEN
!
!              RESCALE FREE ASYMMETRY:
!
               CALL RESCALE_ASYMMETRY(N_PROFILE, 1, N_LAYER
     &            , ASYMMETRY_FREE, FORWARD_SCATTER_FREE
     &            , NPD_PROFILE, NPD_LAYER
     &            )
!
!
               IF (L_CLOUD) THEN
!
!                 RESCALE CLOUDY ASYMMETRY:
!
                  DO K=1, NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)
                     CALL RESCALE_ASYMMETRY(N_PROFILE, N_CLOUD_TOP
     &                  , N_LAYER
     &                  , ASYMMETRY_CLOUD(1, 1, K)
     &                  , FORWARD_SCATTER_CLOUD(1, 1, K)
     &                  , NPD_PROFILE, NPD_LAYER
     &                  )
                  ENDDO
!
               ENDIF
!
            ENDIF
!
         ENDIF
!
!
!
!
!        PRELIMINARY CALCULATIONS FOR SOURCE TERMS:
!
         IF (ISOLIR.EQ.IP_SOLAR) THEN
!           CONVERT NORMALIZED BAND FLUXES TO ACTUAL ENERGY FLUXES.
            DO L=1, N_PROFILE
               INC_SOLAR_FLUX_BAND(L)=SOLAR_TOA(L)
     &            *SOLAR_FLUX_BAND(I_BAND)/SEC_0(L)
            ENDDO
!
         ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
!
!           CALCULATE THE CHANGE IN THE THERMAL SOURCE FUNCTION
!           ACROSS EACH LAYER FOR THE INFRA-RED PART OF THE SPECTRUM.
!
            CALL DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER
     &         , N_DEG_FIT, THERMAL_COEFFICIENT(0, I_BAND)
     &         , T_REF_PLANCK, T_LEVEL, T_GROUND
     &         , PLANCK_SOURCE_BAND, DIFF_PLANCK_BAND
     &         , THERMAL_GROUND_BAND
     &         , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_BAND_2
     &         , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &         , PLANCK_FREEZE_SEA
     &         , NPD_PROFILE, NPD_LAYER, NPD_THERMAL_COEFF
     &         )
         ENDIF
!
!
!
!
!        SET THE SURFACE PROPERTIES:
!
         CALL SET_SURFACE_PROPERTIES(N_POINT_TYPE, INDEX_SURFACE
     &      , I_SPEC_SURFACE
     &      , ISOLIR, I_BAND
     &      , SURFACE_ALBEDO
     &      , ALBEDO_FIELD_DIFF(1, I_BAND), ALBEDO_FIELD_DIR(1, I_BAND)
     &      , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, SEC_0
     &      , EMISSIVITY_GROUND, EMISSIVITY_FIELD(1, I_BAND)
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , THERMAL_GROUND_BAND
     &      , NPD_PROFILE, NPD_BAND, NPD_SURFACE, NPD_ALBEDO_PARM
     &      )
!
!
!
!
!
!        CALL A SOLVER APPROPRIATE TO THE PRESENCE OF GASES AND
!        THE OVERLAP ASSUMED:
!
         IF (.NOT.L_GAS_BAND) THEN
!
!           THERE IS NO GASEOUS ABSORPTION. SOLVE FOR THE
!           FLUXES DIRECTLY.
!
            CALL SOLVE_BAND_WITHOUT_GAS(IERR
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER, D_MASS
!                       Angular integration
     &         , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &         , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of scattering
     &         , I_SCATTER_METHOD_BAND
!                       Options for solver
     &         , I_SOLVER, L_NET, N_AUGMENT
!                       Spectral region
     &         , ISOLIR
!                       Solar Properties
     &         , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &         , PLANCK_SOURCE_BAND(1, 0)
     &         , PLANCK_SOURCE_BAND(1, N_LAYER)
     &         , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &         , THERMAL_GROUND_BAND
!                       Clear-sky optical properties
     &         , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &         , FORWARD_SCATTER_FREE
!                       Cloudy properties
     &         , L_CLOUD, I_CLOUD
!                       Cloudy Geometry
     &         , N_CLOUD_TOP
     &         , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &         , I_REGION_CLOUD, FRAC_REGION
     &         , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &         , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &         , CLOUD_OVERLAP
     &         , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy optical properties
     &         , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &         , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Calculated Fluxes
     &         , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear-sky Fluxes
     &         , L_CLEAR, I_SOLVER_CLEAR
!                       Calculated Clear-sky Fluxes
     &         , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
!                       Planckian Function
     &         , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
!
!
         ELSE
!
!           GASES ARE INCLUDED.
!
!           INITIALIZE THE FLUX IN THE BAND TO ZERO.
            CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &         , ISOLIR
     &         , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
     &         , L_CLEAR
     &         , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
     &         , 0.0E+00
     &         , NPD_PROFILE, NPD_LAYER
     &         , L_NET
     &         )
!
!           TREAT THE GASEOUS OVERLAPS AS DIRECTED BY
!           THE OVERLAP SWITCH.
!
            IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_SINGLE) THEN
!
               CALL SOLVE_BAND_ONE_GAS(IERR
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &            , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND
!                       Options for solver
     &            , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &            , I_BAND, I_GAS
     &            , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &            , K_ESFT, W_ESFT, SCALE_VECTOR
     &            , P_REFERENCE, T_REFERENCE
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &            , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &            , ISOLIR
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)
     &            , DIFF_PLANCK_BAND
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , THERMAL_GROUND_BAND
!                       Clear-sky optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &            , FORWARD_SCATTER_FREE
!                       Cloudy properties
     &            , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &            , N_CLOUD_TOP
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &            , I_REGION_CLOUD, FRAC_REGION
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , CLOUD_OVERLAP
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &            , ASYMMETRY_CLOUD
     &            , FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for clear-sky calculations
     &            , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
!                       Planckian Function
     &            , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &            , NPD_BAND, NPD_SPECIES
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_RANDOM) THEN
!
               CALL SOLVE_BAND_RANDOM_OVERLAP(IERR
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &            , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND
!                       Options for solver
     &            , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &            , I_BAND, N_GAS
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &            , K_ESFT, W_ESFT, SCALE_VECTOR
     &            , P_REFERENCE, T_REFERENCE
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &            , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &            , ISOLIR
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)
     &            , DIFF_PLANCK_BAND
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , THERMAL_GROUND_BAND
!                       Clear-sky optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &            , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &            , N_CLOUD_TOP
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &            , I_REGION_CLOUD, FRAC_REGION
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , CLOUD_OVERLAP
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &            , ASYMMETRY_CLOUD
     &            , FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear-sky Calculations
     &            , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
!                       Planckian Function
     &            , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &            , NPD_BAND, NPD_SPECIES
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_FESFT) THEN
!
               CALL SOLVE_BAND_FESFT(IERR
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &            , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND
!                       Options for solver
     &            , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &            , I_BAND, N_GAS
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &            , K_ESFT, W_ESFT, SCALE_VECTOR
     &            , P_REFERENCE, T_REFERENCE
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &            , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &            , ISOLIR
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)
     &            , DIFF_PLANCK_BAND
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &            , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &            , N_CLOUD_TOP
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &            , I_REGION_CLOUD, FRAC_REGION
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , CLOUD_OVERLAP
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &            , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear Fluxes
     &            , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
!                       Planckian Source Function
     &            , PLANCK_SOURCE_BAND
!                       Sizes of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &            , NPD_BAND, NPD_SPECIES
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_CLR_FESFT) THEN
!
               CALL SOLVE_BAND_CLR_FESFT(IERR
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &            , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &            , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &            , I_BAND, N_GAS
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &            , K_ESFT, W_ESFT, SCALE_VECTOR
     &            , P_REFERENCE, T_REFERENCE
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &            , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &            , ISOLIR
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND
     &            , DIFF_PLANCK_BAND
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &            , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &            , N_CLOUD_TOP
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &            , I_REGION_CLOUD, FRAC_REGION
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , CLOUD_OVERLAP
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &            , ASYMMETRY_CLOUD
     &            , FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear-sky Fluxes
     &            , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &            , NPD_BAND, NPD_SPECIES
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND).EQ.IP_OVERLAP_K_EQV) THEN
!
               CALL SOLVE_BAND_K_EQV(IERR
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &            , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &            , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
     &            , I_BAND, N_GAS
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &            , K_ESFT, W_ESFT, SCALE_VECTOR
     &            , P_REFERENCE, T_REFERENCE
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED
     &            , L_DOPPLER, DOPPLER_CORRECTION
!                       Spectral Region
     &            , ISOLIR
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND
     &            , DIFF_PLANCK_BAND
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , THERMAL_GROUND_BAND
!                       Clear-sky optical properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &            , FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &            , N_CLOUD_TOP
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD
     &            , I_REGION_CLOUD, FRAC_REGION
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , CLOUD_OVERLAP
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &            , ASYMMETRY_CLOUD
     &            , FORWARD_SCATTER_CLOUD
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear-sky Calculations
     &            , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &            , NPD_BAND, NPD_SPECIES
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC
     &            )
!
            ENDIF
         ENDIF
!
!
!
!        INCREMENT THE TOTAL FLUXES.
!
         CALL AUGMENT_TOTAL_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &      , ISOLIR, L_CLEAR, L_NET
     &      , WEIGHT_BAND(I_BAND), PLANCK_SOURCE_BAND
     &      , FLUX_DIRECT, FLUX_TOTAL
     &      , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
     &      , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &      , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
     &      , PLANCK_FLUX
     &      , NPD_PROFILE, NPD_LAYER
     &      , ALBEDO_SURFACE_DIFF
     &      , ALBEDO_SURFACE_DIR
     &      )
!
!
!
!        INCREMENT THE BAND-DEPENDENT DIAGNOSTICS FOR THE
!        UNIFIED MODEL.
         CALL R2_COUPLE_DIAG(N_PROFILE, L_NET, ISOLIR
     &      , ALBEDO_FIELD_DIFF(1, I_BAND), ALBEDO_FIELD_DIR(1, I_BAND)
     &      , ALBEDO_SEA_DIFF(1, I_BAND), ALBEDO_SEA_DIR(1, I_BAND)
     &      , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &      , PLANCK_FREEZE_SEA
     &      , PLANCK_SOURCE_BAND(1, N_LAYER), THERMAL_GROUND_BAND
     &      , FLUX_TOTAL_BAND(1, N_AUGMENT)
     &      , FLUX_TOTAL_BAND(1, N_AUGMENT-1)
     &      , FLUX_DIRECT_BAND(1, N_LAYER)
     &      , FLUX_TOTAL_CLEAR_BAND(1, N_AUGMENT)
     &      , FLUX_TOTAL_CLEAR_BAND(1, N_AUGMENT-1)
     &      , FLUX_DIRECT_CLEAR_BAND(1, N_LAYER)
     &      , WEIGHT_690NM(I_BAND)
     &      , SEA_FLUX
     &      , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &      , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &      , L_SURF_UP_CLR, SURF_UP_CLR
     &      , L_FLUX_BELOW_690NM_SURF, FLUX_BELOW_690NM_SURF
     &      , NPD_PROFILE
     &      )
!
      ENDDO
!
!
!
!
!     PASS THE CALCULATED FLUXES INTO THE OUTPUT ARRAYS.
!
      CALL ASSIGN_FLUX(N_PROFILE, N_LAYER
     &   , FLUX_TOTAL, FLUX_TOTAL_CLEAR
     &   , ISOLIR
     &   , PLANCK_FLUX
     &   , L_CLEAR, L_NET
     &   , FLUX_DOWN, FLUX_UP, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      RETURN
      END
