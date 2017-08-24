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
!+ Shortwave Interface to Edwards-Slingo radiation scheme.
!
! Purpose:
!   This subroutine interface the Edwards-Slingo radiation scheme
!   in the shortwave.
!
! Method:
!   Principally, arrays are transferred to the appropriate formats.
!   Separate subroutines are called for each physical process.
!
! Current Owner of Code: J. M. Edwards
!
! History:
! Version   Date                    Comment
!  4.0    27-07-95                Original Code
!                                (J. M. Edwards)
!       4.1             10-06-96                Checking code for data
!                                               in the spectral files
!                                               added. L_AEROSOL_CCN
!                                               introduced. Coupling
!                                               flux scaled by open-sea
!                                               fraction at sea-ice
!                                               points. Correction of
!                                               heating rates by
!                                               fraction of time for
!                                               which a point is illum-
!                                               inated.
!                                               (J. M. Edwards)
!
!  4.1    22.5.96   Use surface flux at wavelengths below 690nm to
!                   provide photosynthetically active radiation for use
!                   in MOSES boundary layer scheme.  This is added
!                   to the SWOUT array as an 'extra level', without
!                   Zenith Angle adjustement, to enable use in all
!                   physics timesteps.            R.A.Betts
!
!
!       4.2             10-10-96                Climatological aerosol
!                                               introduced.
!                                               (J. M. Edwards)
!  4.4    08-04-97  Changes for new precip scheme (qCF prognostic)
!                                               (A. C. Bushell)
!
!  4.4    29/10/97  Optional prognostic snow albedo scheme introduced
!                                                           R. Essery
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!
!       4.4             04-09-97                Changes to the passing
!                                               of arguments introduced.
!                                               Dissolved sulphate is
!                                               now included in the
!                                               indirect effect.
!                                               Fluxes at the tropopause
!                                               can be diagnosed.
!                                               (J. M. Edwards)
!       4.5     April 1998    Pass soot variables to FILL3A routines
!                                                      Luke Robinson.
!
!       4.5             18-05-98                Obsolete solvers
!                                               removed. New partitioni-
!                                               ing in convective cloud
!                                               introduced.
!                                               (J. M. Edwards)
!
!  4.5    13/05/98  Various changes to argument list to pass an extended
!                   'area' cloud fraction into R2_SET_CLOUD.
!                                                  S.Cusack
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SWRAD(IERR
!                       Mixing Ratios
     &   , H2O, CO2, O3, O2_MIX_RATIO
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D
!                       Pressure Fields
     &   , PSTAR, AB, BB, AC, BC
!                       Temperatures
     &   , TAC
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, GLOBAL_CLOUD_TOP
!                       Stratiform Cloud Fields
     &   , L_CLOUD_WATER_PARTITION
     &   , LCA_AREA, LCA_BULK, LCCWC1, LCCWC2
!                       Convective Cloud Fields
     &   , CCA, CCCWP, CCB, CCT, L_3D_CCA
!                       Surface Fields
     &   , SAL_VIS, SAL_NIR
     &   , LAND_ICE_ALBEDO, OPEN_SEA_ALBEDO, ICE_FRACTION, LAND
     &   , LYING_SNOW
!                       Prognostic snow albedo flag
     &   , L_SNOW_ALBEDO, SAL_DIM
!                       Solar Fields
     &   , COSZIN, LIT, LIST, SCS
!                       Aerosol Fields
     &   , L_CLIMAT_AEROSOL, N_LEVELS_BL
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT
     &   , SULP_DIM1, SULP_DIM2
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE, DISS_SULPHATE
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT
!                       Level of tropopause
     &   , TRINDX
!                       Spectrum
!     ------------------------------------------------------------------
!     ARGUMENT LIST FOR THE REDUCED SW SPECTRAL FILE.
!     (NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW
     &   , NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW
     &   , NPD_SCALE_VARIABLE_SW, NPD_THERMAL_COEFF_SW
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW, NPD_CONTINUUM_SW
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW
     &   , NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW
     &   , L_PRESENT_SW, N_BAND_SW, WAVE_LENGTH_SHORT_SW
     &   , WAVE_LENGTH_LONG_SW, N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW
     &   , SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW
     &   , N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW
     &   , TYPE_ABSORB_SW, I_BAND_ESFT_SW, I_SCALE_ESFT_SW
     &   , I_SCALE_FNC_SW, K_ESFT_SW, W_ESFT_SW
     &   , SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW
     &   , N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW
     &   , I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW
     &   , SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW
     &   , EMISSIVITY_GROUND_SW, N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW
     &   , INDEX_WATER_SW, I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW
     &   , SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW
     &   , I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW
     &   , DROP_PARAMETER_LIST_SW
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW
     &   , N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW
     &   , NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW
     &   , AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW
     &   , AEROSOL_ASYMMETRY_SW
     &   , I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW
     &   , ICE_PARAMETER_LIST_SW
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW
     &   , L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW
!
!     ------------------------------------------------------------------
!                       Algorithmic options
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF CONTROLLING OPTIONS FOR SHORTWAVE RADIATION.
!
     &,
!     ------------------------------------------------------------------
!     VARIABLES FOR CONTROLLING OPTIONS FOR SHORTWAVE RADIATION.
!
!                       Algorithmic Options
     &     I_2STREAM_SW, I_GAS_OVERLAP_SW, I_CLOUD_SW
     &   , I_CLOUD_REPRESENTATION_SW, I_SOLVER_SW
     &   , L_O2_SW
     &   , I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW
     &   , L_LOCAL_CNV_PARTITION_SW
!
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
     &   , PTS
!                       General Diagnostics
     &   , SOLAR_OUT_TOA, L_SOLAR_OUT_TOA
     &   , SOLAR_OUT_CLEAR, L_SOLAR_OUT_CLEAR
     &   , FLUX_BELOW_690NM_SURF, L_FLUX_BELOW_690NM_SURF
     &   , SURFACE_DOWN_FLUX, L_SURFACE_DOWN_FLUX
     &   , SURF_DOWN_CLR, L_SURF_DOWN_CLR
     &   , SURF_UP_CLR, L_SURF_UP_CLR
     &   , LAYER_CLOUD_LIT, L_LAYER_CLOUD_LIT
     &   , CONV_CLOUD_LIT, L_CONV_CLOUD_LIT
     &   , TOTAL_CLOUD_COVER, L_TOTAL_CLOUD_COVER
     &   , CLEAR_HR, L_CLEAR_HR
     &   , NET_FLUX_TROP, L_NET_FLUX_TROP
     &   , UP_FLUX_TROP, L_UP_FLUX_TROP
!                       Microphysical Flag
     &   , L_MICROPHYSICS
!                       Microphysical Diagnostics
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG
     &   , LWP_STRAT, LWP_STRAT_FLAG
     &   , WEIGHTED_RE, WEIGHTED_RE_FLAG
     &   , SUM_WEIGHT_RE, SUM_WEIGHT_RE_FLAG
     &   , NTOT_DIAG, NTOT_DIAG_FLAG
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG
!                       Physical Dimensions
     &   , NLIT
     &   , N_PROFILE, NLEVS, NCLDS
     &   , NWET, NOZONE
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , N_CCA_LEV
!                       Working Dimensions for Diagnostics
     &   , NPDWD_CL_PROFILE
!                       Output
     &   , NETSW, SWSEA, SWOUT
     &   )
!
!
!
      IMPLICIT NONE
!
!
!
!     COMDECKS INCLUDED
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
      REAL SC                           !  Solar constant
      PARAMETER ( SC = 1365. )
!     INTERNAL DIMENSIONS OF THE CODE
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
!     SPECTRAL REGIONS
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
!     ANGULAR INTEGRATION
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
!     TREATMENT OF SCATTERING
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET THE METHODS OF TREATING SCATTERING.
!
      INTEGER
     &     IP_SCATTER_FULL
!             FULL TREATMENT OF SCATTERING
     &   , IP_NO_SCATTER_ABS
!             SCATTERING IGNORED COMPLETELY.
     &   , IP_NO_SCATTER_EXT
!             SCATTERING TREATED AS ABSORPTION
!
      PARAMETER(
     &     IP_SCATTER_FULL=1
     &   , IP_NO_SCATTER_ABS=2
     &   , IP_NO_SCATTER_EXT=3
     &   )
!
!     ------------------------------------------------------------------
!     OPTIONS TO THE CODE ALTERABLE IN THE UM.
!     ------------------------------------------------------------------
!     OPTIONS FOR 3A-RADIATION CODE: VERSION FOR SHORTWAVE
!     CALCULATIONS.
!     (NOTE: SWOPT3A AND SWCAVR3A MUST BE CONSISTENT)
!
      INTEGER
     &     I_2STREAM_SW
!             TWO_STREAM SCHEME SELECTED
     &   , I_GAS_OVERLAP_SW
!             TREATMENT OF GASEOUS OVERLAPS
     &   , I_CLOUD_SW
!             TREATMENT OF CLOUD OVERLAPS
     &   , I_CLOUD_REPRESENTATION_SW
!             REPRESENTATION OF CLOUDS
     &   , I_SOLVER_SW
!             SOLVER SELECTED
!
      LOGICAL
     &     L_O2_SW
!             FLAG FOR OXYGEN
!
!     TYPES OF DROPLETS OR ICE CRYSTALS USED FOR PARAMETRIZATIONS
      INTEGER
     &     I_ST_WATER_SW
!             TYPE FOR STRATIFORM WATER
     &   , I_CNV_WATER_SW
!             TYPE FOR CONVECTIVE WATER
     &   , I_ST_ICE_SW
!             TYPE FOR STRATIFORM ICE
     &   , I_CNV_ICE_SW
!             TYPE FOR CONVECTIVE ICE
      LOGICAL
     &     L_LOCAL_CNV_PARTITION_SW
!             FLAG TO PARTITION CONVECTIVE CLOUDS USING THE LOCAL
!             TEMPERATURE
!
!     ------------------------------------------------------------------
!     OPTIONS TO THE CODE FIXED IN THE UM.
!     ------------------------------------------------------------------
!     MODULE DEFINING OPTIONS TO THE EDWARDS-SLINGO RADIATION CODE
!     FIXED IN THE UNIFIED MODEL. OPTIONS FOR SHORTWAVE CALCULATIONS.
!
!     ALGORITHMIC OPTIONS:
      INTEGER
     &     ISOLIR_SW
!             SPECTRAL REGION
     &   , I_ANGULAR_INTEGRATION_SW
!             METHOD OF ANGULAR INTEGRATION
     &   , I_SCATTER_METHOD_SW
!             TREATMENT OF SCATTERING
!
      LOGICAL
     &     L_LAYER_SW
!             FLAG FOR PROPERTIES IN LAYERS
     &   , L_CLOUD_LAYER_SW
!             FLAG FOR CLOUDY PROPERTIES IN LAYERS
     &   , L_2_STREAM_CORRECT_SW
!             FLAG FOR CORRECTIONS TO 2-STREAM SCHEME
     &   , L_RESCALE_SW
!             FLAG FOR RESCALING OF OPTICAL PROPERTIES
!
!
      PARAMETER(
     &     ISOLIR_SW=IP_SOLAR
     &   , I_ANGULAR_INTEGRATION_SW=IP_TWO_STREAM
     &   , I_SCATTER_METHOD_SW=IP_SCATTER_FULL
     &   , L_LAYER_SW=.TRUE.
     &   , L_CLOUD_LAYER_SW=.TRUE.
     &   , L_2_STREAM_CORRECT_SW=.FALSE.
     &   , L_RESCALE_SW=.TRUE.
     &   )
!
!
!
!     OPTIONS INVOKING PROCESSES:
!
      LOGICAL
     &     L_GAS_SW
!             FLAG FOR GASEOUS ABSORPTION
     &   , L_RAYLEIGH_SW
!             FLAG FOR RAYLEIGH SCATTERING
     &   , L_CONTINUUM_SW
!             FLAG FOR CONTINUUM ABSORPTION
     &   , L_CLOUD_SW
!             FLAG FOR CLOUDS
     &   , L_DROP_SW
!             FLAG FOR DROPLETS
     &   , L_ICE_SW
!             FLAG FOR ICE CRYSTALS
     &   , L_AEROSOL_SW
!             FLAG FOR AEROSOLS
     &   , L_AEROSOL_CCN_SW
!             FLAG AEROSOLS TO DETERMINE CCN
!
      PARAMETER(
     &     L_GAS_SW=.TRUE.
     &   , L_RAYLEIGH_SW=.TRUE.
     &   , L_CONTINUUM_SW=.TRUE.
     &   , L_CLOUD_SW=.TRUE.
     &   , L_DROP_SW=.TRUE.
     &   , L_ICE_SW=.TRUE.
     &   , L_AEROSOL_SW=.TRUE.
     &   , L_AEROSOL_CCN_SW=.TRUE.
     &   )
!
!     ------------------------------------------------------------------
!     NUMERICAL PRECISION
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
     &   )
!
!     ------------------------------------------------------------------
!     SOLVERS
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
!     ERROR FLAGS
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
!     UNIT NUMBERS FOR PRINTED OUTPUT
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             ARRAY SIZES FOR LAYERS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     DIMENSIONS FOR DIAGNOSTIC WORKSPACE
      INTEGER   !, INTENT(IN)
     &     NPDWD_CL_PROFILE
!             NUMBER OF PROFILES ALLOWED IN WORKSPACE FOR
!             CLOUD DIAGNOSTICS
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NWET
!             NUMBER OF WET LEVELS
     &   , NOZONE
!             NUMBER OF LEVELS WITH OZONE
     &   , NLEVS
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
     &   , N_LEVELS_BL
!             NUMBER OF LEVELS IN THE BOUNDARY LAYER
     &   , N_CCA_LEV
!             NUMBER OF CONVECTIVE CLOUD LEVELS
!
!     SPECTRAL DATA:
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED SW-SPECTRAL FILE.
!     (NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT)
!     ------------------------------------------------------------------
!
!
!     DIMENSIONS FOR THE SPECTRUM
!
      INTEGER
     &     NPD_TYPE_SW
!             NUMBER OF TYPES OF DATA IN SW SPECTRUM
     &   , NPD_BAND_SW
!             NUMBER OF SPECTRAL BANDS IN SW SPECTRUM
     &   , NPD_EXCLUDE_SW
!             NUMBER OF EXCLUDED BANDS IN SW SPECTRUM
     &   , NPD_SPECIES_SW
!             NUMBER OF GASEOUS SPECIES IN SW SPECTRUM
     &   , NPD_ESFT_TERM_SW
!             NUMBER OF ESFT TERMS IN SW SPECTRUM
     &   , NPD_SCALE_FNC_SW
!             NUMBER OF SCALING FUNCTIONS IN SW SPECTRUM
     &   , NPD_SCALE_VARIABLE_SW
!             NUMBER OF SCALING VARIABLES IN SW SPECTRUM
     &   , NPD_SURFACE_SW
!             NUMBER OF SURFACE TYPES IN SW SPECTRUM
     &   , NPD_ALBEDO_PARM_SW
!             NUMBER OF ALBEDO PARAMETERS IN SW SPECTRUM
     &   , NPD_CONTINUUM_SW
!             NUMBER OF CONTINUA IN SW SPECTRUM
     &   , NPD_DROP_TYPE_SW
!             NUMBER OF DROP TYPES IN SW SPECTRUM
     &   , NPD_ICE_TYPE_SW
!             NUMBER OF ICE CRYSTAL TYPES IN SW SPECTRUM
     &   , NPD_AEROSOL_SPECIES_SW
!             NUMBER OF AEROSOL SPECIES IN SW SPECTRUM
     &   , NPD_CLOUD_PARAMETER_SW
!             MAX NUMBER OF CLOUD PARAMETERS IN SW SPECTRUM
     &   , NPD_HUMIDITIES_SW
!             MAXIMUM NUMBER OF HUMIDITIES IN SW SPECTRUM
     &   , NPD_THERMAL_COEFF_SW
!             NUMBER OF THERMAL COEFFICIENTS IN SW SPECTRUM
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT_SW(0: NPD_TYPE_SW)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND_SW
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT_SW(NPD_BAND_SW)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_SW(NPD_BAND_SW)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE_SW(NPD_BAND_SW)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_SW(NPD_EXCLUDE_SW, NPD_BAND_SW)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND_SW(NPD_BAND_SW)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT_SW(NPD_BAND_SW)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB_SW
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_SW(NPD_BAND_SW)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_SW(NPD_SPECIES_SW, NPD_BAND_SW)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_SW(NPD_SPECIES_SW)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_SW(NPD_BAND_SW, NPD_SPECIES_SW)
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_SW(NPD_BAND_SW, NPD_SPECIES_SW)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_SW(NPD_BAND_SW, NPD_SPECIES_SW)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT_SW(NPD_ESFT_TERM_SW, NPD_BAND_SW, NPD_SPECIES_SW)
!             ESFT EXPONENTS
     &   , W_ESFT_SW(NPD_ESFT_TERM_SW, NPD_BAND_SW, NPD_SPECIES_SW)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_SW(NPD_SCALE_VARIABLE_SW, NPD_ESFT_TERM_SW
     &        , NPD_BAND_SW, NPD_SPECIES_SW)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_SW(NPD_SPECIES_SW, NPD_BAND_SW)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_SW(NPD_SPECIES_SW, NPD_BAND_SW)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT_SW
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT_SW(0: NPD_THERMAL_COEFF_SW-1
     &   , NPD_BAND_SW)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_SW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE_SW(NPD_SURFACE_SW)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_SW(NPD_SURFACE_SW)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE_SW(NPD_SURFACE_SW)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO_SW(NPD_BAND_SW, NPD_SURFACE_SW)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_SW(0: NPD_ALBEDO_PARM_SW
     &      , NPD_BAND_SW, NPD_SURFACE_SW)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_SW(NPD_BAND_SW, NPD_SURFACE_SW)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM_SW(NPD_BAND_SW)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_SW
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_SW(NPD_SCALE_VARIABLE_SW
     &      , NPD_BAND_SW, NPD_CONTINUUM_SW)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_SW(NPD_CONTINUUM_SW, NPD_BAND_SW)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_SW(NPD_CONTINUUM_SW, NPD_BAND_SW)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION_SW(NPD_DROP_TYPE_SW)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE_SW(NPD_DROP_TYPE_SW)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER_SW
     &        , NPD_BAND_SW, NPD_DROP_TYPE_SW)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_SW(NPD_DROP_TYPE_SW)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_SW(NPD_DROP_TYPE_SW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL_SW
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_SW(NPD_AEROSOL_SPECIES_SW)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_SW(NPD_AEROSOL_SPECIES_SW)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_SW(NPD_AEROSOL_SPECIES_SW)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES_SW(NPD_AEROSOL_SPECIES_SW)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION_SW(NPD_HUMIDITIES_SW
     &        , NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_SW(NPD_HUMIDITIES_SW
     &        , NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_SW(NPD_HUMIDITIES_SW
     &        , NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_SW(NPD_HUMIDITIES_SW, NPD_AEROSOL_SPECIES_SW)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION_SW(NPD_ICE_TYPE_SW)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE_SW(NPD_ICE_TYPE_SW)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER_SW
     &        , NPD_BAND_SW, NPD_ICE_TYPE_SW)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_SW(NPD_ICE_TYPE_SW)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_SW(NPD_ICE_TYPE_SW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT_SW(NPD_SPECIES_SW)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION_SW(NPD_SPECIES_SW)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!    ------------------------------------------------------------------
!
!
!
!     GASEOUS MIXING RATIOS
      REAL      !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)
!             MASS MIXING RATIO OF WATER
     &   , CO2
!             MASS MIXING RATIO OF CO2
     &   , O3(NPD_FIELD, NOZONE)
!             MASS MIXING RATIOS OF OZONE
     &   , O2_MIX_RATIO
!             MASS MIXING RATIO OF OXYGEN
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL      !, INTENT(IN)
     &     PSTAR(NPD_FIELD)
!             SURFACE PRESSURES
     &   , AB(NLEVS+1)
!             A AT BOUNDARIES OF LAYERS
     &   , BB(NLEVS+1)
!             B AT BOUNDARIES OF LAYERS
     &   , AC(NLEVS)
!             A AT CENTRES OF LAYERS
     &   , BC(NLEVS)
!             B AT CENTRES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
!
!     INCIDENT SOLAR RADIATION:
      INTEGER   !, INTENT(IN)
     &     NLIT
!             NUMBER OF LIT POINTS
     &   , LIST(NPD_FIELD)
!             LIST OF LIT POINTS
      REAL      !, INTENT(IN)
     &     COSZIN(NPD_FIELD)
!             COSINES OF ZENITH ANGLE
     &   , SCS
!             SCALING OF SOLAR INCIDENT FIELD
     &   , LIT(NPD_FIELD)
!             FRACTION OF TIME POINT IS LIT
!
!     MICROPHYSICAL FLAG:
      LOGICAL   !, INTENT(IN)
     &     L_MICROPHYSICS
!             FLAG FOR PARAMETRIZED MICROPHYSICS
!
!     OPTIONS FOR TREATING CLOUDS
      LOGICAL   !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
!             TO ENSURE REPRODUCIBLE RESULTS
      INTEGER   !, INTENT(IN)
     &     GLOBAL_CLOUD_TOP
!             GLOBAL TOPMOST CLOUDY LAYER
!
!     PROPERTIES OF STRATIFORM CLOUDS:
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD_WATER_PARTITION
!             FLAG TO USE PROGNOSTIC CLOUD ICE CONTENTS
      REAL      !, INTENT(IN)
     &     LCCWC1(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             NOMINAL LIQUID WATER CONTENTS
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             NOMINAL ICE WATER CONTENTS
     &   , LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             AREA FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
     &   , LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
!
!     PROPERTIES OF CONVECTIVE CLOUDS:
      INTEGER   !, INTENT(IN)
     &     CCB(NPD_FIELD)
!             BASE OF CONVECTIVE CLOUD
     &   , CCT(NPD_FIELD)
!             TOP OF CONVECTIVE CLOUD
      REAL      !, INTENT(IN)
     &     CCCWP(NPD_FIELD)
!             WATER PATH OF CONVECTIVE CLOUD
     &   , CCA(NPD_FIELD,N_CCA_LEV)
!             FRACTION OF CONVECTIVE CLOUD
      LOGICAL   !, INTENT(IN)
     &     L_3D_CCA
!             FLAG FOR 3D convective cloud amount
!
!     AEROSOLS:
      LOGICAL   !, INTENT(IN)
     &     L_CLIMAT_AEROSOL
!             FLAG FOR CLIMATOLOGICAL AEROSOL
      LOGICAL   !, INTENT(IN)
     &     L_USE_SULPC_DIRECT
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
     &   , L_USE_SULPC_INDIRECT
!             FLAG TO USE SULPHUR CYCLE FOR INDIRECT EFFECT
     &   , L_USE_SOOT_DIRECT ! USE DIRECT RAD. EFFECT OF SOOT AEROSOL
      INTEGER   !, INTENT(IN)
     &     SULP_DIM1,SULP_DIM2
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
     &   , SOOT_DIM1, SOOT_DIM2
!          DIMENSIONS FOR SOOT ARRAYS (P_FIELD,P_LEVELS or 1,1)
      REAL      !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIO OF ACCUMULATION MODE AEROSOL
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIO OF AITKEN MODE AEROSOL
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIO OF DISSOLVED SULPHATE
     &,FRESH_SOOT(SOOT_DIM1,SOOT_DIM2),AGED_SOOT(SOOT_DIM1,SOOT_DIM2)
!             SOOT MIXING RATIOS
!
!     CARBON CYCLE:
      LOGICAL   L_CO2_3D    !  controls use of 3D co2 field
      INTEGER   !, INTENT(IN)
     &     CO2_DIM1, CO2_DIM2
!             DIMENSIONS FOR CO2 ARRAY, (P_FIELD,P_LEVELS or 1,1)
      REAL      !, INTENT(IN)
     &     CO2_3D(CO2_DIM1, CO2_DIM2)
!             MASS MIXING RATIO OF CARBON DIOXIDE
!     PROPERTIES OF THE SURFACE:
      LOGICAL   !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND SEA MASK
     &   , L_SNOW_ALBEDO
!             FLAG FOR PROGNOSTIC SNOW ALBEDO
      INTEGER   !, INTENT(IN)
     &     SAL_DIM
!             DIMENSION FOR SAL_VIS AND SAL_NIR 
      REAL      !, INTENT(IN)
     &     ICE_FRACTION(NPD_FIELD)
!             SEA ICE FRACTION
     &   , SAL_VIS(SAL_DIM,2)
!             SURFACE VISIBLE ALBEDO FIELD
     &   , SAL_NIR(SAL_DIM,2)
!             SURFACE NEAR-IR ALBEDO FIELD
     &   , LAND_ICE_ALBEDO(NPD_FIELD)
!             SURFACE ALBEDO OF LAND OR SEA-ICE
     &   , OPEN_SEA_ALBEDO(NPD_FIELD, 2)
!             SURFACE ALBEDO FIELD OF OPEN SEA
!             (DIRECT AND DIFFUSE COMPONENTS)
     &   , LYING_SNOW(NPD_FIELD)
!             MASS LOADING OF LYING SNOW
!
!                       Level of tropopause
      INTEGER
     &     TRINDX(NPD_FIELD)
!             THE LAYER BOUNDARY OF THE TROPOPAUSE
!
!     INCREMENT OF TIME:
      REAL      !, INTENT(IN)
     &     PTS
!             TIME INCREMENT
!
!
!     CALCULATED FLUXES:
      REAL      !, INTENT(OUT)
     &     SWOUT(NPD_FIELD, NLEVS+2)
!             NET DOWNWARD FLUXES
     &   , SWSEA(NPD_FIELD)
!             SEA-SURFACE COMPONENTS OF FLUX
     &   , NETSW(NPD_FIELD)
!             NET ABSORBED SHORTWAVE RADIATION
!
!
!
!     DIAGNOSTICS:
!
!     INPUT SWITCHES:
      LOGICAL   !, INTENT(IN)
     &     L_SOLAR_OUT_TOA
!             REFLECTED SOLAR TOA REQUIRED
     &   , L_SOLAR_OUT_CLEAR
!             CLEAR REFLECTED SOLAR REQUIRED
     &   , L_FLUX_BELOW_690NM_SURF
!             FLUX BELOW 690NM AT SURFACE TO BE DIAGNOSED
     &   , L_SURFACE_DOWN_FLUX
!             DOWNWARD SURFACE FLUX REQUIRED
     &   , L_SURF_DOWN_CLR
!             CALCULATE DOWNWARD CLEAR FLUX
     &   , L_SURF_UP_CLR
!             CALCULATE UPWARD CLEAR FLUX
     &   , L_TOTAL_CLOUD_COVER
!             CALCULATE CLOUD COVER
     &   , L_CLEAR_HR
!             CALCULATE CLEAR-SKY HEATING RATES
     &   , L_NET_FLUX_TROP
!             CALCULATE NET DOWNWARD FLUX AT THE TROPOPAUSE
     &   , L_UP_FLUX_TROP
!             CALCULATE UPWARD FLUX AT THE TROPOPAUSE
!
!     CALCULATED DIAGNOSTICS:
      REAL      !, INTENT(OUT)
     &     SOLAR_OUT_TOA(NPD_FIELD)
!             REFLECTED SOLAR TOA
     &   , SOLAR_OUT_CLEAR(NPD_FIELD)
!             CLEAR REFLECTED SOLAR
     &   , FLUX_BELOW_690NM_SURF(NPD_FIELD)
!             NET SURFACE FLUX BELOW 690NM (AT POINTS WHERE THERE
!             IS SEA-ICE THIS IS WEIGHTED BY THE FRACTION OF OPEN SEA.)
     &   , SURFACE_DOWN_FLUX(NPD_FIELD)
!             DOWNWARD SURFACE FLUX
     &   , SURF_DOWN_CLR(NPD_FIELD)
!             DOWNWARD CLEAR SURFACE FLUX
     &   , SURF_UP_CLR(NPD_FIELD)
!             UPWARD CLEAR SURFACE FLUX
     &   , TOTAL_CLOUD_COVER(NPD_FIELD)
!             TOTAL CLOUD AMOUNT
     &   , CLEAR_HR(NPD_FIELD, NLEVS)
!             CLEAR-SKY HEATING RATES
     &   , NET_FLUX_TROP(NPD_FIELD)
!             NET DOWNWARD FLUX AT THE TROPOPAUSE
     &   , UP_FLUX_TROP(NPD_FIELD)
!             UPWARD FLUX AT THE TROPOPAUSE
!
      LOGICAL   !, INTENT(IN)
     &     L_LAYER_CLOUD_LIT
!             LAYER CLOUD AT LIT POINTS WANTED
     &   , L_CONV_CLOUD_LIT
!             CONVECTIVE CLOUD AT LIT POINTS WANTED
      REAL      !, INTENT(IN)
     &     LAYER_CLOUD_LIT(NPD_FIELD, NCLDS)
!             FRACTION OF LAYER CLOUD LIT
     &   , CONV_CLOUD_LIT(NPD_FIELD)
!             FRACTION OF CONVECTIVE CLOUD LIT
!
!     DIAGNOSTICS FOR THE MRF/UMIST PARAMETRIZATION
!
      LOGICAL   !, INTENT(IN)
     &     RE_CONV_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV_FLAG
!             DIAGNOSE WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT_FLAG
!             DIAGNOSE WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT_FLAG
!             DIAGNOSE LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , WEIGHTED_RE_FLAG
!             CALCULATE OBSERVED EFFECTIVE RADIUS
     &   , SUM_WEIGHT_RE_FLAG
!             CALCULATE SUM OF WEIGHTS FOR EFFECTIVE RADIUS
     &   , NTOT_DIAG_FLAG
!             DIAGNOSE DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG_FLAG
!             DIAGNOSE STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG_FLAG
!             DIAGNOSE SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT_FLAG
!             DIAGNOSE CONDITIONAL SAMPLING WEIGHT
!
      REAL      !, INTENT(OUT)
     &     RE_CONV(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV(NPD_FIELD, NCLDS)
!             WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT(NPD_FIELD, NCLDS)
!             WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT(NPD_FIELD, NCLDS)
!             LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , WEIGHTED_RE(NPD_FIELD)
!             WEIGHTED SUM OF EFFECTIVE RADII
     &   , SUM_WEIGHT_RE(NPD_FIELD)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS
     &   , NTOT_DIAG(NPD_FIELD, NCLDS)
!             DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG(NPD_FIELD, NCLDS)
!             STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG(NPD_FIELD, NCLDS)
!             SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT(NPD_FIELD, NCLDS)
!             CONDITIONAL SAMPLING WEIGHT
!
!
!
!
!
!     LOCAL VARIABLES.
!
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      LOGICAL
     &     L_CLEAR
!             CALCULATE CLEAR-SKY FIELDS
!     FLAGS FOR PROCESSES ACTUALLY ENABLED.
      LOGICAL
     &     L_RAYLEIGH
!             LOCAL FLAG FOR RAYLEIGH SCATTERING
     &   , L_GAS
!             LOCAL FLAG FOR GASEOUS ABSORPTION
     &   , L_CONTINUUM
!             LOCAL FLAG FOR CONTINUUM ABSORPTION
     &   , L_DROP
!             LOCAL FLAG FOR SCATTERING BY DROPLETS
     &   , L_AEROSOL
!             LOCAL FLAG FOR SCATTERING BY AEROSOLS
     &   , L_AEROSOL_CCN
!             LOCAL FLAG TO USE AEROSOLS TO DETERMINE CCN
     &   , L_ICE
!             LOCAL FLAG FOR SCATTERING BY ICE CRYSTALS
      INTEGER
     &     I_SOLVER_CLEAR
!             SOLVER FOR CLEAR-SKY FLUXES
     &   , I_GAS_OVERLAP(NPD_BAND_SW)
!             OVERLAPS IN EACH BAND
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE FIELD
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES_SW)
!             MASS FRACTIONS OF GASES
     &   , NULLMMR
!             NULL MASS MIXING RATIO
      PARAMETER(
     &     NULLMMR=0.0E+00
     &   )
!
!     CLOUDY PROPERTIES:
      INTEGER
     &     N_CONDENSED
!             NUMBER OF CONDENSED PHASES
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF CONDENSED COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATION SCHEMES FOR COMPONENTS
     &   , N_CLOUD_TOP_GLOBAL
!             INVERTED GLOBAL TOPMOST CLOUDY LAYER
      REAL
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER_SW
     &        , NPD_CLOUD_COMPONENT, NPD_BAND_SW)
!             PARAMETERS FOR CONDENSED PHASES
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             CHARACTERISTIC DIMENSIONS OF CONDENSED SPECIES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             MASS FRACTIONS OF CONDENSED SPECIES
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)
!             MINIMUM DIMENSIONS OF CONDENSED COMPONENTS
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSIONS OF CONDENSED COMPONENTS
!
!     PROPERTIES OF AEROSOLS:
      REAL
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_AEROSOL_SPECIES_SW)
!             MIXING RATIOS OF AEROSOLS
!
!     SOLAR FIELDS:
      REAL
     &     SEC_0(NPD_PROFILE)
!             SECANTS OF ZENITH ANGLE
     &   , SOLAR_INCIDENT_NORM(NPD_PROFILE)
!             NORMALLY INCIDENT SOLAR IRRADIANCE
!
!     SURFACE PROPERTIES:
      LOGICAL
     &     LAND_G(NPD_PROFILE)
!             GATHERED SURFACE MASK
      INTEGER
     &     I_SURFACE(NPD_PROFILE)
!             TYPES OF SURFACE AT EACH POINT
      REAL
     &     ALBEDO_FIELD_DIFF_GREY(NPD_PROFILE)
!             DIFFUSE ALBEDO FIELD
     &   , ALBEDO_FIELD_DIR_GREY(NPD_PROFILE)
!             DIRECT ALBEDO FIELD
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_SW)
!             DIFFUSE ALBEDO FIELD
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_SW)
!             DIRECT ALBEDO FIELD
     &   , EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_SW)
!             EMISSIVITY FIELD
     &   , ALBEDO_SEA_DIFF_G(NPD_PROFILE, NPD_BAND_SW)
!             GATHERED DIFFUSE ALBEDO FOR OPEN SEA
     &   , ALBEDO_SEA_DIR_G(NPD_PROFILE, NPD_BAND_SW)
!             GATHERED DIRECT ALBEDO FOR OPEN SEA
!
!     FLUXES:
      REAL
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_NET(NPD_PROFILE, 0: NPD_LAYER)
!             NET/DOWNWARD FLUX
     &   , FLUX_NET_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY NET/DOWNWARD TOTAL FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY UPWARD FLUX
!
!     ARRAYS FOR USE WITH DIAGNOSTICS:
      REAL
     &     WEIGHT_690NM(NPD_BAND_SW)
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM
     &   , W_CLOUD_DIAG(NPDWD_CL_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS FOR DIAGNOSTIC USE
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL
     &     SEA_FLUX_G(NPD_PROFILE)
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX_G(NPD_PROFILE)
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR_G(NPD_PROFILE)
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR_G(NPD_PROFILE)
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF_G(NPD_PROFILE)
!             GATHERED SURFACE FLUX BELOW 690NM
!
!     FIELDS REQUIRED FOR CALL TO RADIATION CODE BUT NOT USED
      INTEGER
     &     N_ORDER_GAUSS
     &   , I_GAS
      LOGICAL
     &     L_SWITCH_SCATTER(NPD_BAND_SW)
!
!     AUXILIARY VARIABLES:
      REAL
     &     CPBYG
!             SPECIFIC HEAT BY GRAVITY
     &   , DACON
!             DIFFERENCE IN A's
     &   , DBCON
!             DIFFERENCE IN B's
     &   , WEIGHT_BAND(NPD_BAND_SW)
!             WEIGHTING FACTORS FOR BANDS
      PARAMETER(CPBYG=CP/G)
!
!     VARIABLES REQUIRED FOR COMPATIBILITY WITH SUBROUTINES:
      INTEGER
     &     N_FRAC_ICE_POINT
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
      REAL
     &     DUMMY
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     R2_SET_GAS_MIX_RATIO, R2_SET_THERMODYNAMIC
     &   , R2_SET_AEROSOL_FIELD, R2_SET_CLOUD_FIELD
     &   , R2_SET_CLOUD_PARAMETRIZATION
     &   , R2_SET_SURFACE_FIELD_SW, R2_ZERO_1D
     &   , R2_INIT_MRF_UMIST_DIAG
     &   , R2_COMPARE_PROC
!
!
!
!
!
!
!     INITIALIZE THE ERROR FLAG FOR THE RADIATION CODE.
      IERR=I_NORMAL
!
!     INITIALIZATIONS FOR DIAGNOSTICS DEPENDING ON BANDS
!
      IF ( L_FLUX_BELOW_690NM_SURF .OR. L_SNOW_ALBEDO ) THEN
         CALL R2_SET_690NM_WEIGHT(N_BAND_SW
     &      , L_PRESENT_SW
     &      , N_BAND_EXCLUDE_SW
     &      , INDEX_EXCLUDE_SW
     &      , WAVE_LENGTH_SHORT_SW
     &      , WAVE_LENGTH_LONG_SW
     &      , WEIGHT_690NM
     &      , NPD_BAND_SW, NPD_EXCLUDE_SW, NPD_TYPE_SW
     &      )
      ENDIF
!                                                                       
!     COMPARE PROCESSES IN THE SPECTRAL FILE WITH THOSE ENABLED IN
!     THE CODE.
      CALL R2_COMPARE_PROC(IERR, L_PRESENT_SW
     &   , L_RAYLEIGH_SW, L_GAS_SW, L_CONTINUUM_SW
     &   , L_DROP_SW, L_AEROSOL_SW, L_AEROSOL_CCN_SW, L_ICE_SW
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT
     &   , L_USE_SOOT_DIRECT
     &   , L_CLIMAT_AEROSOL
     &   , L_RAYLEIGH, L_GAS, L_CONTINUUM
     &   , L_DROP, L_AEROSOL, L_AEROSOL_CCN, L_ICE
     &   , NPD_TYPE_SW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!
!     SET THE PROPERTIES OF THE SURFACE
      CALL R2_SET_SURFACE_FIELD_SW(
     &     N_BAND_SW
     &   , NLIT, LIST
     &   , I_SURFACE, I_SPEC_SURFACE_SW
     &   , L_SURFACE_SW
     &   , L_MICROPHYSICS, L_SNOW_ALBEDO, SAL_DIM
     &   , LAND, OPEN_SEA_ALBEDO, LAND_ICE_ALBEDO, ICE_FRACTION
     &   , SAL_VIS, SAL_NIR, WEIGHT_690NM   
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF
     &   , LAND_G, ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G
     &   , NPD_FIELD, NPD_PROFILE, NPD_BAND_SW, NPD_SURFACE_SW
     &   )
!
!     SET THE MIXING RATIOS OF GASES.
      CALL R2_SET_GAS_MIX_RATIO(IERR
     &   , NLIT, NLEVS, NWET, NOZONE
     &   , LIST
     &   , N_ABSORB_SW, TYPE_ABSORB_SW
     &   , .FALSE., .FALSE., .FALSE., .FALSE., L_O2_SW
     &   , .FALSE., .FALSE., .FALSE., .FALSE.
     &   , H2O, CO2, O3, NULLMMR, NULLMMR, NULLMMR, NULLMMR
     &   , O2_MIX_RATIO
     &   , NULLMMR, NULLMMR, NULLMMR, NULLMMR
     &   , GAS_MIX_RATIO
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_SPECIES_SW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!     SET THE THERMODYNAMIC PROPERTIES OF THE ATMOSPHERE.
      CALL R2_SET_THERMODYNAMIC(NLIT, NLEVS, LIST, .FALSE.
     &   , PSTAR, DUMMY, AB, BB, AC, BC
     &   , DUMMY, TAC
     &   , P, T, DUMMY, DUMMY, D_MASS
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &   )
!
!
!     SET THE MIXING RATIOS OF AEROSOLS.
      IF (L_AEROSOL.OR.L_AEROSOL_CCN) THEN
         CALL R2_SET_AEROSOL_FIELD(IERR
     &      , NLIT, NLEVS, N_AEROSOL_SW, TYPE_AEROSOL_SW
     &      , LIST
     &      , L_CLIMAT_AEROSOL, N_LEVELS_BL
     &      , L_USE_SULPC_DIRECT
     &      , SULP_DIM1, SULP_DIM2
     &      , ACCUM_SULPHATE, AITKEN_SULPHATE
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT
     &      , LAND, LYING_SNOW, PSTAR, AB, BB, TRINDX
     &      , AEROSOL_MIX_RATIO
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_SW
     &      )
      ENDIF
!
!
!     ASSIGN THE PROPERTIES OF CLOUDS.
!
      CALL R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND_SW
     &   , I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW
     &   , L_DROP_TYPE_SW
     &   , I_DROP_PARAMETRIZATION_SW
     &   , DROP_PARAMETER_LIST_SW
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW
     &   , L_ICE_TYPE_SW
     &   , I_ICE_PARAMETRIZATION_SW
     &   , ICE_PARAMETER_LIST_SW
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , NPD_BAND_SW
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
      CALL R2_INIT_MRF_UMIST_DIAG(IERR
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG
     &   , LWP_STRAT, LWP_STRAT_FLAG
     &   , NTOT_DIAG, NTOT_DIAG_FLAG
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG
     &   , NPD_FIELD, NPD_PROFILE, NCLDS
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
      CALL R2_SET_CLOUD_FIELD(NLIT, NLEVS, NCLDS
     &   , LIST
     &   , P, T, D_MASS
     &   , CCB, CCT, CCA, CCCWP
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK
     &   , L_MICROPHYSICS, L_AEROSOL_CCN
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE
     &   , L_CLOUD_WATER_PARTITION,  LAND_G
     &   , I_CLOUD_REPRESENTATION_SW, I_CONDENSED_PARAM
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , W_CLOUD, FRAC_CLOUD, L_LOCAL_CNV_PARTITION_SW
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG
     &   , LWP_STRAT, LWP_STRAT_FLAG
     &   , NTOT_DIAG, NTOT_DIAG_FLAG
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_SW
     &   , N_CCA_LEV, L_3D_CCA
     &   )
!
      IF (WEIGHTED_RE_FLAG.AND.SUM_WEIGHT_RE_FLAG) THEN
         CALL R2_CLOUD_LEVEL_DIAG(IERR, NLIT, NLEVS, NCLDS
     &      , LIST
     &      , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW
     &      , W_CLOUD, FRAC_CLOUD
     &      , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &      , WEIGHTED_RE_FLAG, WEIGHTED_RE, SUM_WEIGHT_RE
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
      ENDIF
!
!
!
!     SET THE INCIDENT SOLAR FLUX.
      DO L=1, NLIT
         SOLAR_INCIDENT_NORM(L)=SCS*SC*LIT(LIST(L))
         SEC_0(L)=1.0E+00/COSZIN(LIST(L))
      ENDDO
!
!
!     CHECK THAT A VALID NUMBER HAS BEEN SUPPLIED FOR THE SOLVER.
      IF ( (I_SOLVER_SW.NE.IP_SOLVER_PENTADIAGONAL).AND.
     &     (I_SOLVER_SW.NE.IP_SOLVER_MIX_11).AND.
     &     (I_SOLVER_SW.NE.IP_SOLVER_MIX_DIRECT).AND.
     &     (I_SOLVER_SW.NE.IP_SOLVER_HOMOGEN_DIRECT).AND.
     &     (I_SOLVER_SW.NE.IP_SOLVER_TRIPLE)
     &   ) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: AN INVALID SOLVER HAS BEEN SELECTED '
     &      , 'IN THE SHORTWAVE REGION.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET CLEAR-SKY CALCULATIONS.
      L_CLEAR=L_SOLAR_OUT_CLEAR.OR.
     &        L_SURF_DOWN_CLR.OR.
     &        L_CLEAR_HR
!
      IF (L_CLEAR) THEN
!
!        SELECT A CLEAR-SKY SOLVER TO MATCH THE MAIN SOLVER.
         IF (I_SOLVER_SW.EQ.IP_SOLVER_PENTADIAGONAL) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_SW.EQ.IP_SOLVER_MIX_11) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_SW.EQ.IP_SOLVER_MIX_DIRECT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW.EQ.IP_SOLVER_HOMOGEN_DIRECT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW.EQ.IP_SOLVER_TRIPLE) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ENDIF
!
      ENDIF
!
!
!     SET PROPERTIES FOR INDIVIDUAL BANDS.
      DO I=1, N_BAND_SW
         WEIGHT_BAND(I)=1.0E+00
         I_GAS_OVERLAP(I)=I_GAS_OVERLAP_SW
      ENDDO
!
!
!     INVERT THE TOPMOST CLOUDY LAYER IF USING A GLOBAL VALUE.
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP_GLOBAL=NLEVS+1-GLOBAL_CLOUD_TOP
      ENDIF
!
!
!
!
      CALL FLUX_CALC(IERR
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM
     &   , L_CLOUD_SW, L_DROP, L_ICE
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION_SW, I_2STREAM_SW, L_2_STREAM_CORRECT_SW
     &   , L_RESCALE_SW, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_SW, L_SWITCH_SCATTER
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL
!                       Options for Solver
     &   , I_SOLVER_SW
!                       General Spectral Properties
     &   , N_BAND_SW, 1, N_BAND_SW
     &   , WEIGHT_BAND
!                       General Atmospheric Properties
     &   , NLIT, NLEVS
     &   , L_LAYER_SW, L_CLOUD_LAYER_SW
     &   , P, T, DUMMY, DUMMY, D_MASS
!                       Spectral Region
     &   , ISOLIR_SW
!                       Solar Fields
     &   , SEC_0, SOLAR_INCIDENT_NORM, SOLAR_FLUX_BAND_SW
     &   , RAYLEIGH_COEFFICIENT_SW
!                       Infra-red Fields
     &   , N_DEG_FIT_SW
     &   , THERMAL_COEFFICIENT_SW
     &   , T_REF_PLANCK_SW, .FALSE.
!                       Gaseous Absorption
     &   , N_ABSORB_SW, I_GAS_OVERLAP, I_GAS
     &   , GAS_MIX_RATIO
     &   , N_BAND_ABSORB_SW, INDEX_ABSORB_SW
     &   , I_BAND_ESFT_SW
     &   , W_ESFT_SW, K_ESFT_SW
     &   , I_SCALE_ESFT_SW, I_SCALE_FNC_SW
     &   , SCALE_VECTOR_SW
     &   , P_REFERENCE_SW, T_REFERENCE_SW
!                       Doppler Broadening
     &   , L_DOPPLER_PRESENT_SW
     &   , DOPPLER_CORRECTION_SW
!                       Surface Fields
     &   , L_SURFACE_SW, I_SURFACE
     &   , I_SPEC_SURFACE_SW
     &   , SURFACE_ALBEDO_SW
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR
     &   , N_DIR_ALBEDO_FIT_SW
     &   , DIRECT_ALBEDO_PARM_SW
     &   , EMISSIVITY_GROUND_SW
     &   , EMISSIVITY_FIELD
!                       Continuum Absorption
     &   , N_BAND_CONTINUUM_SW
     &   , INDEX_CONTINUUM_SW, INDEX_WATER_SW
     &   , K_CONTINUUM_SW, I_SCALE_FNC_CONT_SW
     &   , SCALE_CONTINUUM_SW
     &   , P_REF_CONTINUUM_SW
     &   , T_REF_CONTINUUM_SW
!                       Properties of Aerosols
     &   , N_AEROSOL_SW
     &   , AEROSOL_MIX_RATIO
     &   , AEROSOL_ABSORPTION_SW
     &   , AEROSOL_SCATTERING_SW
     &   , AEROSOL_ASYMMETRY_SW
     &   , I_AEROSOL_PARAMETRIZATION_SW
     &   , NHUMIDITY_SW
     &   , HUMIDITIES_SW
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW, W_CLOUD, FRAC_CLOUD
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_NET, FLUX_UP
!                       Options for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_NET_CLEAR, FLUX_UP_CLEAR
!                       Arrays specific to the UM
!                       Arrays for Coupling
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G
     &   , SEA_FLUX_G
!                       Arrays for diagnostics specific to the UM
     &   , L_FLUX_BELOW_690NM_SURF, WEIGHT_690NM
     &   , FLUX_BELOW_690NM_SURF_G
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX_G
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR_G
     &   , L_SURF_UP_CLR, SURF_UP_CLR_G
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , NPD_BAND_SW
     &   , NPD_SPECIES_SW
     &   , NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW
     &   , NPD_SCALE_VARIABLE_SW
     &   , NPD_CONTINUUM_SW
     &   , NPD_AEROSOL_SPECIES_SW
     &   , NPD_HUMIDITIES_SW
     &   , NPD_CLOUD_PARAMETER_SW
     &   , NPD_THERMAL_COEFF_SW
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!     PREPARE THE OUTPUT ARRAYS:
!
!     ZERO SWOUT SO THAT POINTS LYING IN THE NIGHT WILL CONTAIN VALID
!     FLUXES AFTER SCATTERING.
      DO I=1, NLEVS+1
         CALL R2_ZERO_1D(N_PROFILE, SWOUT(1, I))
      ENDDO
      IF (L_CLEAR_HR) THEN
         DO I=1, NLEVS
            CALL R2_ZERO_1D(N_PROFILE, CLEAR_HR(1, I))
         ENDDO
      ENDIF
!
!     SCATTER THE NET DOWNWARD FLUX AT EACH LEVEL INTO SWOUT.
      DO I=1, NLEVS+1
         DO L=1, NLIT
            SWOUT(LIST(L), I)=FLUX_NET(L, NLEVS+1-I)
         ENDDO
      ENDDO
!
!
!     NET SHORTWAVE RADIATION ABSORBED BY THE PLANET
!     (I. E. EARTH AND ATMOSPHERE TOGETHER):
!
      CALL R2_ZERO_1D(N_PROFILE, NETSW)
      DO L=1, NLIT
         NETSW(LIST(L))=SWOUT(LIST(L), NLEVS+1)
      ENDDO
!
!
!
!
!     ASSIGNMENT OF DIAGNOSTICS:
!
!     TOTAL CLOUD COVER:
!
      IF (L_TOTAL_CLOUD_COVER) THEN
!
!        THE CLOUD AMOUNTS MUST BE RECALCULATED SINCE W_CLOUD
!        AS DEFINED ABOVE HOLDS VALUES ONLY AT LIT POINTS.
!        A DIFFERENTLY DEFINED DIAGNOSTIC ARRAY IS USED TO PREVENT
!        OUR HAVING TO DECLARE A LOT OF SPACE FOR W_CLOUD.
         IF (L_3D_CCA) THEN
         DO I=NLEVS+1-NCLDS, NLEVS
            DO L=1, N_PROFILE
               W_CLOUD_DIAG(L,I) = CCA(L,NLEVS+1-I)
     &                +(1.0E+00-CCA(L,NLEVS+1-I))*LCA_AREA(L,NLEVS+1-I)
              ENDDO
           ENDDO
         ELSE
           DO I=NLEVS+1-NCLDS, NLEVS
             DO L=1, N_PROFILE
               IF ( (CCT(L).GE.NLEVS+2-I).AND.(CCB(L).LE.NLEVS+1-I) )
     &            THEN
                  W_CLOUD_DIAG(L, I)
     &               =CCA(L,1)+(1.0E+00-CCA(L,1))*LCA_AREA(L, NLEVS+1-I)
               ELSE
                  W_CLOUD_DIAG(L, I)=LCA_AREA(L, NLEVS+1-I)
               ENDIF
            ENDDO
         ENDDO
         ENDIF
!
         CALL R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, NLEVS, NCLDS
     &      , I_CLOUD_SW, W_CLOUD_DIAG, TOTAL_CLOUD_COVER
     &      , NPDWD_CL_PROFILE, NPD_LAYER
     &      )
!
      ENDIF
!
!
!     AMOUNT OF CONVECTIVE CLOUD AT DAYLIT POINTS.
      IF (L_CONV_CLOUD_LIT) THEN
!        ZERO THE ARRAY EVERYWHERE AND FILL ONLY AT LIT POINTS.
         CALL R2_ZERO_1D(N_PROFILE, CONV_CLOUD_LIT)
        IF (L_3D_CCA) THEN
         DO L=1, NLIT
            CONV_CLOUD_LIT(LIST(L))=CCA( LIST(L),CCT(LIST(L)) )
         ENDDO
        ELSE
          DO L=1, NLIT
            CONV_CLOUD_LIT(LIST(L))=CCA(LIST(L),1)
          ENDDO
        ENDIF
      ENDIF
!
!
!     AMOUNT OF STRATIFORM CLOUD AT DAYLIT POINTS.
      IF (L_LAYER_CLOUD_LIT) THEN
         DO I=1, NCLDS
            CALL R2_ZERO_1D(N_PROFILE, LAYER_CLOUD_LIT(1, I))
            DO L=1, NLIT
               LAYER_CLOUD_LIT(LIST(L), I)=LCA_AREA(LIST(L), I)
            ENDDO
         ENDDO
      ENDIF
!
!
!     OUTGOING SOLAR RADIATION AT TOA:
!
      IF (L_SOLAR_OUT_TOA) THEN
         CALL R2_ZERO_1D(N_PROFILE, SOLAR_OUT_TOA)
         DO L=1, NLIT
            SOLAR_OUT_TOA(LIST(L))=SOLAR_INCIDENT_NORM(L)/SEC_0(L)
     &         -FLUX_NET(L, 0)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY OUTGOING SOLAR RADIATION AT TOA:
!
      IF (L_SOLAR_OUT_CLEAR) THEN
         CALL R2_ZERO_1D(N_PROFILE, SOLAR_OUT_CLEAR)
         DO L=1, NLIT
            SOLAR_OUT_CLEAR(LIST(L))=SOLAR_INCIDENT_NORM(L)/SEC_0(L)
     &         -FLUX_NET_CLEAR(L, 0)
         ENDDO
      ENDIF
!
!
!     SURFACE FLUX BELOW 690NM.
!
      IF (L_FLUX_BELOW_690NM_SURF) THEN
         CALL R2_ZERO_1D(N_PROFILE, FLUX_BELOW_690NM_SURF)
         DO L=1, NLIT
            IF (LAND(LIST(L))) THEN
               FLUX_BELOW_690NM_SURF(LIST(L))
     &            =FLUX_BELOW_690NM_SURF_G(L)
            ELSE
               FLUX_BELOW_690NM_SURF(LIST(L))
     &            =FLUX_BELOW_690NM_SURF_G(L)
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
            ENDIF
         ENDDO
      ENDIF
!
!
!     DOWNWARD FLUX AT THE SURFACE:
!
      IF (L_SURFACE_DOWN_FLUX) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURFACE_DOWN_FLUX)
         DO L=1, NLIT
            SURFACE_DOWN_FLUX(LIST(L))=SURFACE_DOWN_FLUX_G(L)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY DOWNWARD FLUX AT THE SURFACE:
!
      IF (L_SURF_DOWN_CLR) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_CLR)
         DO L=1, NLIT
            SURF_DOWN_CLR(LIST(L))=SURF_DOWN_CLR_G(L)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY UPWARD FLUX AT THE SURFACE:
!
      IF (L_SURF_UP_CLR) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURF_UP_CLR)
         DO L=1, NLIT
            SURF_UP_CLR(LIST(L))=SURF_UP_CLR_G(L)
         ENDDO
      ENDIF
!
!
!     NET FLUX AT THE TROPOPAUSE:
!
      IF (L_NET_FLUX_TROP) THEN
         CALL R2_ZERO_1D(N_PROFILE, NET_FLUX_TROP)
         DO L=1, NLIT
            NET_FLUX_TROP(LIST(L))
     &         =FLUX_NET(L, NLEVS+1-TRINDX(LIST(L)))
         ENDDO
      ENDIF
!
!
!     UPWARD FLUX AT THE TROPOPAUSE:
!
      IF (L_UP_FLUX_TROP) THEN
         CALL R2_ZERO_1D(N_PROFILE, UP_FLUX_TROP)
         DO L=1, NLIT
            UP_FLUX_TROP(LIST(L))
     &         =FLUX_UP(L, NLEVS+1-TRINDX(LIST(L)))
         ENDDO
      ENDIF
!
!
!
!
!
!     FINAL PROCESSING OF OUTPUT FIELDS
!
!     CONVERT THE FLUXES TO INCREMENTS.
      DO I=NLEVS, 1, -1
!
         DACON=(AB(I)-AB(I+1))*CPBYG/PTS
         DBCON=(BB(I)-BB(I+1))*CPBYG/PTS
         DO L=1, N_PROFILE
            SWOUT(L, I+1)=(SWOUT(L, I+1)-SWOUT(L, I))
     &         /(DACON+PSTAR(L)*DBCON)
         ENDDO
!
         IF (L_CLEAR_HR) THEN
            DO L=1, NLIT
               CLEAR_HR(LIST(L), I)=(FLUX_NET_CLEAR(L, NLEVS-I)
     &            -FLUX_NET_CLEAR(L, NLEVS+1-I))
     &            /(PTS*(DACON+PSTAR(LIST(L))*DBCON))
            ENDDO
         ENDIF
!
      ENDDO
!
!
!
!     SEPARATE CONTRIBUTIONS OVER OPEN SEA.
!     SEA_FLUX_G IS NOT WEIGHTED BY THE FRACTION OF ICE.
      CALL R2_ZERO_1D(N_PROFILE, SWSEA)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L=1, NLIT
         IF (.NOT.LAND(LIST(L))) THEN
            SWSEA(LIST(L))=(1.0E+00-ICE_FRACTION(LIST(L)))
     &         *SEA_FLUX_G(L)
            SWOUT(LIST(L), 1)=SWOUT(LIST(L), 1)-SWSEA(LIST(L))
         ENDIF
      ENDDO
!
!
!     DIVIDE FLUX_BELOW_690NM_SURF BY LAND ALBEDO TO GIVE TOTAL
!     DOWNWARD FLUX OF PHOTOSYTHETICALLY ACTIVE RADIATION.  ADD THIS
!     TO THE SWOUT ARRAY AS AN EXTRA 'LEVEL' TO ENABLE USE IN NON-
!     RADIATION TIMESTEPS.
      IF (L_FLUX_BELOW_690NM_SURF) THEN
        DO L=1, N_PROFILE
           SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /
     &        (1 - LAND_ICE_ALBEDO(L))
        ENDDO
      ELSE
        DO L=1, N_PROFILE
           SWOUT(L, NLEVS+2)=0.0
        ENDDO
      ENDIF
!
!
!     DIVIDE BY COSINE OF SOLAR ZENITH ANGLE TO PROVIDE VALUES FOR
!     UPPER ROUTINES. THIS APPLIES ONLY TO SWOUT. THE MACHINE TOLERANCE
!     IS ADDED TO MAINTAIN CONDITIONING.
      DO I=1, NLEVS+2
         DO L=1, N_PROFILE
            SWOUT(L, I)=SWOUT(L, I)/(COSZIN(L)*LIT(L)+TOL_MACHINE)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to set surface fields.
!
! Purpose:
!   The albedos and emissivity of the surface are set.
!
! Method:
!   Straightforward. Though the arrays passed to the code may depend
!   on the spectral band, the input arrays have no spectral dependence.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_SURFACE_FIELD_SW(
     &     N_BAND
     &   , NLIT, LIST
     &   , I_SURFACE, I_SPEC_SURFACE, L_SURFACE
     &   , L_MICROPHYSICS, L_SNOW_ALBEDO, SAL_DIM
     &   , LAND, OPEN_SEA_ALBEDO, LAND_ICE_ALBEDO, ICE_FRACTION
     &   , SAL_VIS, SAL_NIR, WEIGHT_690NM  
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF
     &   , LAND_G, ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR
     &   , NPD_FIELD, NPD_PROFILE, NPD_BAND_SW, NPD_SURFACE_SW
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET PERMITTED METHODS OF SPECIFYING THE
!     SURFACE ALBEDO AND EMISSIVITY.
!
      INTEGER
     &     IP_SURFACE_SPECIFIED
!             PROPERTIES SPECIFIED BY SURFACE TYPE
     &   , IP_SURFACE_INTERNAL
!             PROPERTIES PASSED INTO CODE
     &   , IP_SURFACE_POLYNOMIAL
!             DIRECT ALBEDO FITTED AS POLYNOMIAL
     &   , IP_SURFACE_PAYNE
!             FIT IN THE FUNCTIONAL FORM USED BY PAYNE
!
      PARAMETER(
     &     IP_SURFACE_SPECIFIED=1
     &   , IP_SURFACE_INTERNAL=2
     &   , IP_SURFACE_POLYNOMIAL=3
     &   , IP_SURFACE_PAYNE=4
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES:
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF INPUT FIELDS
     &   , NPD_PROFILE
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
     &   , NPD_BAND_SW
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_SURFACE_SW
!             MAXIMUM NUMBER OF SURFACES
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , SAL_DIM
!             DIMENSION OF SAL_VIS AND SAL_NIR
!
!     LIT POINTS:
      INTEGER   !, INTENT(IN)
     &     NLIT
!             NUMBER OF LIT POINTS
     &   , LIST(NPD_FIELD)
!             LIST OF SUNLIT POINTS
!
!     PROPERTIES OF SURFACES
      INTEGER   !, INTENT(OUT)
     &     I_SURFACE(NPD_PROFILE)
!             TYPES OF SURFACES
     &   , I_SPEC_SURFACE(NPD_SURFACE_SW)
      LOGICAL   !, INTENT(OUT)
     &     L_SURFACE(NPD_SURFACE_SW)
!             FLAGS FOR TYPES OF SURFACES
!
!     PHYSICAL PROPERTIES OF SURFACES:
      LOGICAL   !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND MASK
      REAL      !, INTENT(IN)
     &     OPEN_SEA_ALBEDO(NPD_FIELD, 2)
!             DIFFUSE ALBEDO FIELD
     &   , LAND_ICE_ALBEDO(NPD_FIELD)
!             DIRECT ALBEDO FIELD
     &   , SAL_VIS(SAL_DIM,2)
!             VISIBLE ALBEDO FIELD
     &   , SAL_NIR(SAL_DIM,2)
!             NEAR-IR ALBEDO FIELD
     &   , WEIGHT_690NM(NPD_BAND_SW)
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM
     &   , ICE_FRACTION(NPD_FIELD)
!             FRACTION OF SEA ICE
!
!     MISCELLANEOUS INPUTS
      LOGICAL   !, INTENT(IN)
     &     L_MICROPHYSICS
!             FLAG TO CALCULATE MICROPHYSICS
     &   , L_SNOW_ALBEDO
!             FLAG FOR PROGNOSTIC SNOW ALBEDO
!
!
!     SURFACE PROPERTIES SET.
      REAL      !, INTENT(OUT)
     &     EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_SW)
!             EMISSIVITIES OF SURFACES
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_SW)
!             DIFFUSE ALBEDO OF SURFACE
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_SW)
!             DIRECT ALBEDO OF SURFACE
!
!     GATHERED SURFACE FIELDS
      LOGICAL   !, INTENT(OUT)
     &     LAND_G(NPD_PROFILE)
!             GATHERED LAND FLAGS
      REAL      !, INTENT(OUT)
     &     ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND_SW)
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND_SW)
!             DIRECT ALBEDO OF OPEN SEA
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
!     OVERRIDE ANY SURFACE PROERTIES READ IN FROM THE SPECTRAL FILE.
      DO L=1, NLIT
         I_SURFACE(L)=1
      ENDDO
      L_SURFACE(1)=.TRUE.
      I_SPEC_SURFACE(1)=IP_SURFACE_INTERNAL
!
!
      IF (L_MICROPHYSICS) THEN
!        GATHER THE ARRAY OF SURFACE FLAGS IF THE MICROPHYSICS
!        IS PARAMETRIZED.
         DO L=1, NLIT
            LAND_G(L)=LAND(LIST(L))
         ENDDO
      ENDIF
!
!
!     SET THE ALBEDO FIELDS: AN AVERAGE ALBEDO IS REQUIRED OVER WHERE
!     THERE IS SEA-ICE. SEPARATE ALBEDOS ARE PROVIDED FOR LAND/ICE
!     OR FOR OPEN SEA. BAND-DEPENDENT COPIES OF THE ALBEDOS MUST BE
!     MADE FOR CALCULATING COUPLING FLUXES.
!
      DO I=1, N_BAND
         DO L=1, NLIT
!
            EMISSIVITY_FIELD(L, I)=0.0E+00
!
            IF (.NOT.LAND(LIST(L))) THEN
               ALBEDO_FIELD_DIFF(L, I)
     &            =LAND_ICE_ALBEDO(LIST(L))*ICE_FRACTION(LIST(L))
     &            +OPEN_SEA_ALBEDO(LIST(L), 2)
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
               ALBEDO_FIELD_DIR(L, I)
     &            =LAND_ICE_ALBEDO(LIST(L))*ICE_FRACTION(LIST(L))
     &            +OPEN_SEA_ALBEDO(LIST(L), 1)
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
               ALBEDO_SEA_DIR(L, I)=OPEN_SEA_ALBEDO(LIST(L), 1)
               ALBEDO_SEA_DIFF(L, I)=OPEN_SEA_ALBEDO(LIST(L), 2)
            ELSE
               IF ( L_SNOW_ALBEDO ) THEN
                 ALBEDO_FIELD_DIFF(L,I) = 
     &                                WEIGHT_690NM(I)*SAL_VIS(LIST(L),2)
     &                       + (1. - WEIGHT_690NM(I))*SAL_NIR(LIST(L),2)
                 ALBEDO_FIELD_DIR(L,I) = 
     &                                WEIGHT_690NM(I)*SAL_VIS(LIST(L),1)
     &                       + (1. - WEIGHT_690NM(I))*SAL_NIR(LIST(L),1)
               ELSE    
               ALBEDO_FIELD_DIFF(L, I)=LAND_ICE_ALBEDO(LIST(L))
               ALBEDO_FIELD_DIR(L, I)=LAND_ICE_ALBEDO(LIST(L))
               ENDIF
               ALBEDO_SEA_DIR(L, I)=0.0E+00
               ALBEDO_SEA_DIFF(L, I)=0.0E+00
            ENDIF
!
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to calculate weights for the flux below 690 nm.
!
! Purpose:
!   Weights to calculate the flux below 690 nm are set.
!
! Method:
!   Straightforward. The flux is assumed to be linearly distributed
!   across bands.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_690NM_WEIGHT(N_BAND
     &   , L_PRESENT
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
     &   , WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG
     &   , WEIGHT_690NM
     &   , NPD_BAND_SW, NPD_EXCLUDE_SW, NPD_TYPE_SW
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY VARIABLES:
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_BAND_SW
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_EXCLUDE_SW
!             MAXIMUM NUMBER OF EXCLUDED REGIONS
     &   , NPD_TYPE_SW
!             MAXIMUM NUMBER OF TYPES OF SPECTRAL DATA
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , N_BAND_EXCLUDE(NPD_BAND_SW)
!             NUMBER OF EXCLUDED REGIONS IN BANDS
     &   , INDEX_EXCLUDE(NPD_EXCLUDE_SW, NPD_BAND_SW)
!             INDICES OF EXCLUDED REGIONS IN BANDS
!
      LOGICAL   !, INTENT(IN)
     &     L_PRESENT(0: NPD_TYPE_SW)
!             FLAG FOR TYPES OF SPECTRAL DATA PRESENT
!
      REAL      !, INTENT(IN)
     &     WAVE_LENGTH_SHORT(NPD_BAND_SW)
!             SHORT WAVELENGTH LIMITS OF BANDS
     &   , WAVE_LENGTH_LONG(NPD_BAND_SW)
!             LONG WAVELENGTH LIMITS OF BANDS
!
!
!     WEIGHTS SET.
      REAL      !, INTENT(OUT)
     &     WEIGHT_690NM(NPD_BAND_SW)
!             WEIGHTS APPLYING TO EACH BAND
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
      REAL
     &     TOTAL_ENERGY_RANGE
!             TOTAL RANGE OF ENERGIES COVERED BY BAND
     &   , ENERGY_RANGE_BELOW_690NM
!             RANGE OF ENERGIES IN BAND BELOW 690 NM
!
!
!
      DO I=1, N_BAND
         IF (WAVE_LENGTH_LONG(I).LT.6.9E-07) THEN
            WEIGHT_690NM(I)=1.0E+00
         ELSE IF (WAVE_LENGTH_SHORT(I).GT.6.9E-07) THEN
            WEIGHT_690NM(I)=0.0E+00
         ELSE
!
            ENERGY_RANGE_BELOW_690NM=1.0E+00/WAVE_LENGTH_SHORT(I)
     &         -1.0E+00/6.9E-07
            TOTAL_ENERGY_RANGE=1.0E+00/WAVE_LENGTH_SHORT(I)
     &         -1.0E+00/WAVE_LENGTH_LONG(I)
            IF (L_PRESENT(14)) THEN
!              REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
               DO J=1, N_BAND_EXCLUDE(I)
                  IF (WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I)).LT.
     &               6.9E-07) THEN
                     ENERGY_RANGE_BELOW_690NM=ENERGY_RANGE_BELOW_690NM
     &                  -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))
     &                  +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                  ELSE IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)).LT.
     &               6.9E-07) THEN
                     ENERGY_RANGE_BELOW_690NM=ENERGY_RANGE_BELOW_690NM
     &                  -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))
     &                  +1.0E+00/6.9E-07
                  ENDIF
                  TOTAL_ENERGY_RANGE=TOTAL_ENERGY_RANGE
     &               -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))
     &               +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
               ENDDO
            ENDIF
!
            WEIGHT_690NM(I)=ENERGY_RANGE_BELOW_690NM/TOTAL_ENERGY_RANGE
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to initialize diagnostics for MRF/UMIST parametrization.
!
! Purpose:
!   Checks are made for consistency of the diagnostic requests and the
!   arrays are filled with zeros at all points.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_INIT_MRF_UMIST_DIAG(IERR
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG
     &   , LWP_STRAT, LWP_STRAT_FLAG
     &   , NTOT_DIAG, NTOT_DIAG_FLAG
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG
     &   , NPD_FIELD, N_PROFILE, NCLDS
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
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
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES:
!
!     ERROR FLAG
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             ACTUAL SIZE OF INPUT ARRAY
!
!     SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS

!     DIAGNOSTICS FOR THE MRF/UMIST PARAMETRIZATION
!
      LOGICAL
     &     RE_CONV_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV_FLAG
!             DIAGNOSE WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT_FLAG
!             DIAGNOSE WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT_FLAG
!             DIAGNOSE LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG_FLAG
!             DIAGNOSE DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG_FLAG
!             DIAGNOSE STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG_FLAG
!             DIAGNOSE SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT_FLAG
!             DIAGNOSE CONDITIONAL SAMPLING WEIGHT
!
      REAL
     &     RE_CONV(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT(NPD_FIELD, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV(NPD_FIELD, NCLDS)
!             WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT(NPD_FIELD, NCLDS)
!             WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT(NPD_FIELD, NCLDS)
!             LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG(NPD_FIELD, NCLDS)
!             DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG(NPD_FIELD, NCLDS)
!             STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG(NPD_FIELD, NCLDS)
!             SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT(NPD_FIELD, NCLDS)
!             CONDITIONAL SAMPLING WEIGHT
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
!
!
!
      IF (RE_CONV_FLAG) THEN
         IF (.NOT.WGT_CONV_FLAG) THEN
            WRITE(IU_ERR, '(/A, /A)')
     &         '*** ERROR: MICROPHYSICAL DIAGNOSTICS FOR CONVECTIVE'
     &         , 'CLOUD MUST INCLUDE THE CLOUD WEIGHTING.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF ( (RE_STRAT_FLAG).OR.(LWP_STRAT_FLAG) ) THEN
         IF (.NOT.WGT_STRAT_FLAG) THEN
            WRITE(IU_ERR, '(/A, /A)')
     &         '*** ERROR: MICROPHYSICAL DIAGNOSTICS FOR STRATIFORM'
     &         , 'CLOUD MUST INCLUDE THE CLOUD WEIGHTING.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
      DO I=1, NCLDS
         IF (WGT_CONV_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, WGT_CONV(1, I))
         IF (RE_CONV_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, RE_CONV(1, I))
         IF (WGT_STRAT_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, WGT_STRAT(1, I))
         IF (RE_STRAT_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, RE_STRAT(1, I))
         IF (LWP_STRAT_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, LWP_STRAT(1, I))
         IF (NTOT_DIAG_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, NTOT_DIAG(1, I))
         IF (STRAT_LWC_DIAG_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, STRAT_LWC_DIAG(1, I))
         IF (SO4_CCN_DIAG_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, SO4_CCN_DIAG(1, I))
         IF (COND_SAMP_WGT_FLAG)
     &      CALL R2_ZERO_1D(N_PROFILE, COND_SAMP_WGT(1, I))
      ENDDO
!
!
!
      RETURN
      END
