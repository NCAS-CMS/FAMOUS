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
!+ Longwave Interface to the Edwards-Slingo Radiation Scheme.
!
! Purpose:
!   This routine prepares the call to the Edwards-Slingo radiation
!   scheme in the longwave.
!
! Method:
!   Principally, this routine transfers arrays into the correct formats.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Revised formulation
!                                               over sea-ice. Testing
!                                               of spectral options
!                                               introduced. New solvers
!                                               added.
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.2             08-08-96                Climatological aerosols
!                                               introduced.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!
!       4.4             04-09-96                Changes to the passing
!                                               of arguments into the
!                                               routine. Dissolved
!                                               sulphate aerosol is
!                                               now included in the
!                                               indirect effect.
!                                               Diagnostics of fluxes
!                                               at the tropopause
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for new (H)(C)FCs
!                                               added. New option
!                                               for treating convective
!                                               partitioning added.
!                                               Code for obsolete
!                                               solvers removed.

!                                               (J. M. Edwards)
!
!       4.5     April 1998    Pass soot variables to FILL3A routines
!                                                      Luke Robinson.
!       4.5     June  1998    Various changes to argument list to pass
!                             an extended 'area' cloud fraction into
!                             R2_SET_CLOUD.              S. Cusack
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_LWRAD(IERR
!                       Gaseous Mixing Ratios
     &   , H2O, CO2, O3
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D
     &   , N2O_MIX_RATIO, CH4_MIX_RATIO
     &   , CFC11_MIX_RATIO, CFC12_MIX_RATIO, CFC113_MIX_RATIO
     &   , HCFC22_MIX_RATIO, HFC125_MIX_RATIO, HFC134A_MIX_RATIO
!                       Thermodynamic Variables
     &   , TAC, PEXNER, TSTAR, PSTAR, AB, BB, AC, BC
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, GLOBAL_CLOUD_TOP
!                       Stratiform Cloud Fields
     &   , L_CLOUD_WATER_PARTITION
     &   , LCA_AREA, LCA_BULK, LCCWC1, LCCWC2
!                       Convective Cloud Fields
     &   , CCA, CCCWP, CCB, CCT, L_3D_CCA
!                       Surface Fields
     &   , LAND, ICE_FRACTION
     &   , LYING_SNOW
!                       Aerosol Fields
     &   , L_CLIMAT_AEROSOL, N_LEVELS_BL
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT
     &   , SULP_DIM1,SULP_DIM2
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE, DISS_SULPHATE
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT
!                       Level of tropopause
     &   , TRINDX
!                       Spectrum
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF LW SPECTRAL DATA.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW
     &   , NPD_SCALE_VARIABLE_LW
     &   , NPD_THERMAL_COEFF_LW
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW
     &   , NPD_CONTINUUM_LW
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW
     &   , L_PRESENT_LW
     &   , N_BAND_LW, WAVE_LENGTH_SHORT_LW, WAVE_LENGTH_LONG_LW
     &   , N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW
     &   , N_ABSORB_LW, N_BAND_ABSORB_LW, INDEX_ABSORB_LW
     &   , TYPE_ABSORB_LW
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW
     &   , K_ESFT_LW, W_ESFT_LW
     &   , SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW
     &   , EMISSIVITY_GROUND_LW
     &   , N_BAND_CONTINUUM_LW, INDEX_CONTINUUM_LW, INDEX_WATER_LW
     &   , I_SCALE_FNC_CONT_LW, K_CONTINUUM_LW
     &   , SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW, T_REF_CONTINUUM_LW
     &   , I_DROP_PARAMETRIZATION_LW, L_DROP_TYPE_LW
     &   , DROP_PARAMETER_LIST_LW
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW
     &   , I_AEROSOL_PARAMETRIZATION_LW
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW
     &   , AEROSOL_ASYMMETRY_LW
     &   , I_ICE_PARAMETRIZATION_LW, L_ICE_TYPE_LW
     &   , ICE_PARAMETER_LIST_LW
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW
!
!     ------------------------------------------------------------------
!                       Algorithmic Options
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF CONTROLLING OPTIONS FOR THE LONGWAVE RADIATION.
!
     &   ,
!     ------------------------------------------------------------------
!     VARIABLES FOR CONTROLLING OPTIONS FOR THE LONGWAVE RADIATION.
!     (NOTE: LWOPT3A AND LWCAVR3A MUST BE CONSISTENT)
!
     &     I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW
     &   , I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW, I_SOLVER_LW
     &   , L_N2O_LW, L_CH4_LW, L_CFC11_LW , L_CFC12_LW, L_CFC113_LW
     &   , L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW
     &   , I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW
     &   , I_CNV_ICE_LW, L_MICROPHYSICS_LW, L_LOCAL_CNV_PARTITION_LW
!
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
     &   , PTS
!                       General Diagnostics
     &   , TOTAL_CLOUD_COVER, L_TOTAL_CLOUD_COVER
     &   , CLEAR_OLR, L_CLEAR_OLR
     &   , SURFACE_DOWN_FLUX, L_SURFACE_DOWN_FLUX
     &   , SURF_DOWN_CLR, L_SURF_DOWN_CLR
     &   , CLEAR_HR, L_CLEAR_HR
     &   , NET_FLUX_TROP, L_NET_FLUX_TROP
     &   , DOWN_FLUX_TROP, L_DOWN_FLUX_TROP
!                       Physical Dimensions
     &   , N_PROFILE, NLEVS, NCLDS
     &   , NWET, NOZONE, NPD_FIELD
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , N_CCA_LEV
!                       Output Fields
     &   , OLR, LWSEA, LWOUT
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
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
!     METHODS OF INTEGRATION
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
!     METHODS OF SCATTERING
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
!     OPTIONS TO THE CODE ALTERABLE IN THE UM
!     ------------------------------------------------------------------
!     OPTIONS FOR 3A-RADIATION CODE: VERSION FOR LONGWAVE
!     CALCULATIONS.
!     (NOTE: LWOPT3A AND LWCAVR3A MUST BE CONSISTENT)
!
      INTEGER
     &     I_2STREAM_LW
!             TWO-STREAM SCHEME
     &   , I_GAS_OVERLAP_LW
!             TREATMENT OF GASEOUS OVERLAPS
     &   , I_CLOUD_LW
!             TREATMENT OF CLOUDY OVERLAPS
     &   , I_CLOUD_REPRESENTATION_LW
!             REPRESENTATION OF CLOUDS
     &   , I_SOLVER_LW
!             SOLVER SELECTED
      LOGICAL
     &     L_IR_SOURCE_QUAD_LW
!             REPRESENTATION OF THE IR-SOURCE TERM
     &   , L_MICROPHYSICS_LW
!             FLAG FOR MICROPHYSICS IN LONG WAVE
     &   , L_LOCAL_CNV_PARTITION_LW
!             FLAG TO PARTITION CONVECTIVE CLOUD USING THE
!             LOCAL TEMPERATURE.
!     OPTIONS FOR TRACE GASES:
      LOGICAL
     &     L_N2O_LW
!             FLAG FOR NITROUS OXIDE
     &   , L_CH4_LW
!             FLAG FOR METHANE
     &   , L_CFC11_LW
!             FLAG FOR CFC11
     &   , L_CFC12_LW
!             FLAG FOR CFC12
     &   , L_CFC113_LW
!             FLAG FOR CFC113
     &   , L_HCFC22_LW
!             FLAG FOR HCFC22
     &   , L_HFC125_LW
!             FLAG FOR HFC125
     &   , L_HFC134A_LW
!             FLAG FOR HFC134A
!
!     TYPES OF DROPLETS OR ICE CRYSTALS USED FOR PARAMETRIZATIONS
      INTEGER
     &     I_ST_WATER_LW
!             TYPE FOR STRATIFORM WATER
     &   , I_CNV_WATER_LW
!             TYPE FOR CONVECTIVE WATER
     &   , I_ST_ICE_LW
!             TYPE FOR STRATIFORM ICE
     &   , I_CNV_ICE_LW
!             TYPE FOR CONVECTIVE ICE
!
!     ------------------------------------------------------------------
!     OPTIONS TO THE CODE FIXED IN THE UM
!     ------------------------------------------------------------------
!     MODULE DEFINING OPTIONS TO THE EDWARDS-SLINGO RADIATION CODE
!     FIXED IN THE UNIFIED MODEL. OPTIONS FOR LONGWAVE CALCULATIONS.
!
!     ALGORITHMIC OPTIONS:
      INTEGER
     &     ISOLIR_LW
!             SPECTRAL REGION
     &   , I_ANGULAR_INTEGRATION_LW
!             METHOD OF ANGULAR INTEGRATION
     &   , I_SCATTER_METHOD_LW
!             TREATMENT OF SCATTERING
!
      LOGICAL
     &     L_LAYER_LW
!             FLAG FOR PROPERTIES IN LAYERS
     &   , L_CLOUD_LAYER_LW
!             FLAG FOR CLOUDY PROPERTIES IN LAYERS
     &   , L_2_STREAM_CORRECT_LW
!             FLAG FOR CORRECTIONS TO 2-STREAM SCHEME
     &   , L_RESCALE_LW
!             FLAG FOR RESCALING OF OPTICAL PROPERTIES
!
!
      PARAMETER(
     &     ISOLIR_LW=IP_INFRA_RED
     &   , I_ANGULAR_INTEGRATION_LW=IP_TWO_STREAM
     &   , I_SCATTER_METHOD_LW=IP_SCATTER_FULL
     &   , L_LAYER_LW=.TRUE.
     &   , L_CLOUD_LAYER_LW=.TRUE.
     &   , L_2_STREAM_CORRECT_LW=.FALSE.
     &   , L_RESCALE_LW=.TRUE.
     &   )
!
!
!
!     OPTIONS INVOKING PROCESSES:
!
      LOGICAL
     &     L_GAS_LW
!             FLAG FOR GASEOUS ABSORPTION
     &   , L_RAYLEIGH_LW
!             FLAG FOR RAYLEIGH SCATTERING
     &   , L_CONTINUUM_LW
!             FLAG FOR CONTINUUM ABSORPTION
     &   , L_CLOUD_LW
!             FLAG FOR CLOUDS
     &   , L_DROP_LW
!             FLAG FOR DROPLETS
     &   , L_ICE_LW
!             FLAG FOR ICE CRYSTALS
     &   , L_AEROSOL_LW
!             FLAG FOR AEROSOLS
     &   , L_AEROSOL_CCN_LW
!             FLAG TO USE AEROSOLS TO DETERMINE CCN
!
      PARAMETER(
     &     L_GAS_LW=.TRUE.
     &   , L_RAYLEIGH_LW=.FALSE.
     &   , L_CONTINUUM_LW=.TRUE.
     &   , L_CLOUD_LW=.TRUE.
     &   , L_DROP_LW=.TRUE.
     &   , L_ICE_LW=.TRUE.
     &   , L_AEROSOL_LW=.TRUE.
     &   , L_AEROSOL_CCN_LW=.TRUE.
     &   )
!
!     ------------------------------------------------------------------
!     NUMERICAL PRECISIONS
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
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
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
!     PHYSICAL CONSTANTS
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE SETTING PHYSICAL CONSTANTS.
!
      REAL
     &     MOL_WEIGHT_AIR
!             MOLAR WEIGHT OF DRY AIR
     &   , N2_MASS_FRAC
!             MASS FRACTION OF NITROGEN
!
      PARAMETER(
     &     MOL_WEIGHT_AIR=28.966E-3
     &   , N2_MASS_FRAC=0.781E+00
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
!
!     SPECTRAL DATA:
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED LW-SPECTRAL FILE.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!     ------------------------------------------------------------------
!
!
!     DIMENSIONS FOR THE REDUCED LW SPECTRAL FILE
!
      INTEGER
     &     NPD_TYPE_LW
!             NUMBER OF TYPES OF DATA IN LW SPECTRUM
     &   , NPD_BAND_LW
!             NUMBER OF SPECTRAL BANDS IN LW SPECTRUM
     &   , NPD_EXCLUDE_LW
!             NUMBER OF EXCLUDED BANDS IN LW SPECTRUM
     &   , NPD_SPECIES_LW
!             NUMBER OF GASEOUS SPECIES IN LW SPECTRUM
     &   , NPD_ESFT_TERM_LW
!             NUMBER OF ESFT TERMS IN LW SPECTRUM
     &   , NPD_SCALE_FNC_LW
!             NUMBER OF SCALING FUNCTIONS IN LW SPECTRUM
     &   , NPD_SCALE_VARIABLE_LW
!             NUMBER OF SCALING VARIABLES IN LW SPECTRUM
     &   , NPD_SURFACE_LW
!             NUMBER OF SURFACE TYPES IN LW SPECTRUM
     &   , NPD_ALBEDO_PARM_LW
!             NUMBER OF ALBEDO PARAMETERS IN LW SPECTRUM
     &   , NPD_CONTINUUM_LW
!             NUMBER OF CONTINUA IN LW SPECTRUM
     &   , NPD_DROP_TYPE_LW
!             NUMBER OF DROP TYPES IN LW SPECTRUM
     &   , NPD_ICE_TYPE_LW
!             NUMBER OF ICE CRYSTAL TYPES IN LW SPECTRUM
     &   , NPD_AEROSOL_SPECIES_LW
!             NUMBER OF AEROSOL SPECIES IN LW SPECTRUM
     &   , NPD_CLOUD_PARAMETER_LW
!             MAX NUMBER OF CLOUD PARAMETERS IN LW SPECTRUM
     &   , NPD_HUMIDITIES_LW
!             MAXIMUM NUMBER OF HUMIDITIES IN LW SPECTRUM
     &   , NPD_THERMAL_COEFF_LW
!             NUMBER OF THERMAL COEFFICIENTS IN LW SPECTRUM
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT_LW(0: NPD_TYPE_LW)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND_LW
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT_LW(NPD_BAND_LW)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_LW(NPD_BAND_LW)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE_LW(NPD_BAND_LW)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_LW(NPD_EXCLUDE_LW, NPD_BAND_LW)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND_LW(NPD_BAND_LW)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT_LW(NPD_BAND_LW)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB_LW
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_LW(NPD_BAND_LW)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_LW(NPD_SPECIES_LW, NPD_BAND_LW)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_LW(NPD_SPECIES_LW)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_LW(NPD_BAND_LW, NPD_SPECIES_LW)
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_LW(NPD_BAND_LW, NPD_SPECIES_LW)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_LW(NPD_BAND_LW, NPD_SPECIES_LW)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT_LW(NPD_ESFT_TERM_LW, NPD_BAND_LW, NPD_SPECIES_LW)
!             ESFT EXPONENTS
     &   , W_ESFT_LW(NPD_ESFT_TERM_LW, NPD_BAND_LW, NPD_SPECIES_LW)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_LW(NPD_SCALE_VARIABLE_LW, NPD_ESFT_TERM_LW
     &        , NPD_BAND_LW, NPD_SPECIES_LW)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_LW(NPD_SPECIES_LW, NPD_BAND_LW)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_LW(NPD_SPECIES_LW, NPD_BAND_LW)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT_LW
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT_LW(0: NPD_THERMAL_COEFF_LW-1
     &   , NPD_BAND_LW)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_LW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE_LW(NPD_SURFACE_LW)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_LW(NPD_SURFACE_LW)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE_LW(NPD_SURFACE_LW)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO_LW(NPD_BAND_LW, NPD_SURFACE_LW)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_LW(0: NPD_ALBEDO_PARM_LW
     &      , NPD_BAND_LW, NPD_SURFACE_LW)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_LW(NPD_BAND_LW, NPD_SURFACE_LW)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM_LW(NPD_BAND_LW)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_LW
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_LW(NPD_SCALE_VARIABLE_LW
     &      , NPD_BAND_LW, NPD_CONTINUUM_LW)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_LW(NPD_CONTINUUM_LW, NPD_BAND_LW)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_LW(NPD_CONTINUUM_LW, NPD_BAND_LW)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION_LW(NPD_DROP_TYPE_LW)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE_LW(NPD_DROP_TYPE_LW)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER_LW
     &        , NPD_BAND_LW, NPD_DROP_TYPE_LW)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_LW(NPD_DROP_TYPE_LW)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_LW(NPD_DROP_TYPE_LW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL_LW
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_LW(NPD_AEROSOL_SPECIES_LW)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_LW(NPD_AEROSOL_SPECIES_LW)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_LW(NPD_AEROSOL_SPECIES_LW)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES_LW(NPD_AEROSOL_SPECIES_LW)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION_LW(NPD_HUMIDITIES_LW
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_LW(NPD_HUMIDITIES_LW
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_LW(NPD_HUMIDITIES_LW
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_LW(NPD_HUMIDITIES_LW, NPD_AEROSOL_SPECIES_LW)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION_LW(NPD_ICE_TYPE_LW)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE_LW(NPD_ICE_TYPE_LW)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER_LW
     &        , NPD_BAND_LW, NPD_ICE_TYPE_LW)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_LW(NPD_ICE_TYPE_LW)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_LW(NPD_ICE_TYPE_LW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT_LW(NPD_SPECIES_LW)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION_LW(NPD_SPECIES_LW)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!    ------------------------------------------------------------------
!
!     GASEOUS MIXING RATIOS:
      REAL      !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)
!             MASS MIXING RATIO OF WATER
     &   , CO2
!             MASS MIXING RATIO OF CO2
     &   , O3(NPD_FIELD, NOZONE)
!             MASS MIXING RATIOS OF OZONE
     &   , N2O_MIX_RATIO
!             MASS MIXING RATIO OF NITROUS OXIDE
     &   , CH4_MIX_RATIO
!             MASS MIXING RATIO OF METHANE
     &   , CFC11_MIX_RATIO
!             MASS MIXING RATIO OF CFC11
     &   , CFC12_MIX_RATIO
!             MASS MIXING RATIO OF CFC12
     &   , CFC113_MIX_RATIO
!             MASS MIXING RATIO OF CFC113
     &   , HCFC22_MIX_RATIO
!             MASS MIXING RATIO OF HCFC22
     &   , HFC125_MIX_RATIO
!             MASS MIXING RATIO OF HFC125
     &   , HFC134A_MIX_RATIO
!             MASS MIXING RATIO OF HFC134A
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL      !, INTENT(IN)
     &     AB(NLEVS+1)
!             A AT BOUNDARIES OF LAYERS
     &   , BB(NLEVS+1)
!             B AT BOUNDARIES OF LAYERS
     &   , AC(NLEVS)
!             A AT CENTRES OF LAYERS
     &   , BC(NLEVS)
!             B AT CENTRES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
     &   , PEXNER(NPD_FIELD, NLEVS+1)
!             Exner FUNCTION AT BOUNDARIES
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
!             LIQUID WATER CONTENTS (THESE ARE NOT USED DIRECTLY IN
!             THE RADIATION: THE TOTAL CONDENSED WATER CONTENT IS
!             REPARTITIONED USING FOCWWIL).
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             ICE WATER CONTENTS (THESE ARE NOT USED DIRECTLY IN
!             THE RADIATION: THE TOTAL CONDENSED WATER CONTENT IS
!             REPARTITIONED USING FOCWWIL).
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
!             FRACTION OF GRID-BOX COVERED BY CONVECTIVE CLOUD
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
      INTEGER   !,INTENT (IN)
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
!     SURFACE FIELDS:
      LOGICAL   !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND SEA MASK
      REAL      !, INTENT(IN)
     &     PSTAR(NPD_FIELD)
!             SURFACE PRESSURES
     &   , TSTAR(NPD_FIELD)
!             SURFACE TEMPERATURES
     &   , ICE_FRACTION(NPD_FIELD)
!             SEA ICE FRACTION
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
     &     OLR(NPD_FIELD)
!             NET OUTGOING RADIATION
     &   , LWOUT(NPD_FIELD, NLEVS+1)
!             NET DOWNWARD FLUXES OR HEATING RATES
     &   , LWSEA(NPD_FIELD)
!             SEA-SURFACE COMPONENTS OF FLUX
!
!
!
!     DIAGNOSTICS:
!
!     INPUT SWITCHES:
      LOGICAL   !, INTENT(IN)
     &     L_TOTAL_CLOUD_COVER
!             TOTAL CLOUD AMOUNT DIAGNOSED
     &   , L_CLEAR_OLR
!             CLEAR OLR DIAGNOSED
     &   , L_SURFACE_DOWN_FLUX
!             SURFACE DOWNWARD FLUX DIAGNOSED
     &   , L_SURF_DOWN_CLR
!             SURFACE DOWNWARD CLEAR FLUX DIAG.
     &   , L_CLEAR_HR
!             CALCULATE CLEAR-SKY HEATING RATES
     &   , L_NET_FLUX_TROP
!             CALCULATE NET DOWNWARD FLUX AT THE TROPOPAUSE
     &   , L_DOWN_FLUX_TROP
!             CALCULATE DOWNWARD FLUX AT THE TROPOPAUSE
!
!     CALCULATED DIAGNOSTICS:
      REAL      !, INTENT(OUT)
     &     TOTAL_CLOUD_COVER(NPD_FIELD)
!             TOTAL CLOUD COVER
     &   , CLEAR_OLR(NPD_FIELD)
!             CLEAR-SKY OLR
     &   , SURFACE_DOWN_FLUX(NPD_FIELD)
!             DOWNWARD SURFACE FLUX
     &   , SURF_DOWN_CLR(NPD_FIELD)
!             DOWNWARD SURFACE CLEARFLUX
     &   , CLEAR_HR(NPD_FIELD, NLEVS)
!             CLEAR-SKY HEATING RATES
     &   , NET_FLUX_TROP(NPD_FIELD)
!             NET DOWNWARD FLUX AT THE TROPOPAUSE
     &   , DOWN_FLUX_TROP(NPD_FIELD)
!             DOWNWARD FLUX AT THE TROPOPAUSE
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
      INTEGER
     &     I_GATHER(NPD_FIELD)
!             GATHERING ARRAY
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
!             LOCAL FLAG FOR SCATTERING BY AEROSOLS
     &   , L_ICE
!             LOCAL FLAG FOR SCATTERING BY ICE CRYSTALS
      INTEGER
     &     I_SOLVER_CLEAR
!             SOLVER FOR CLEAR-SKY FLUXES
     &   , I_GAS_OVERLAP(NPD_BAND_LW)
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
     &   , T_BDY(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE FIELD AT BOUNDARIES
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES_LW)
!             MASS FRACTIONS OF GASES
!
!     SURFACE FIELDS:
      INTEGER
     &     I_SURFACE(NPD_PROFILE)
!             TYPE OF SURFACE AT THE FOOT OF EACH PROFILE
      REAL      !, INTENT(IN)
     &     ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_LW)
!             DIFFUSE ALBEDOS
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_LW)
!             DIRECT ALBEDOS
     &   , EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_LW)
!             EMISSIVITIES
     &   , ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND_LW)
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND_LW)
!             DIRECT ALBEDO OF OPEN SEA
     &   , T_SURFACE(NPD_PROFILE)
!             GATHERED TEMPERATURE OF SURFACE
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
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER_LW
     &        , NPD_CLOUD_COMPONENT, NPD_BAND_LW)
!             PARAMETERS FOR CONDENSED PHASES
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             CHARACTERISTIC DIMENSIONS OF CONDENSED SPECIES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             MASS FRACTIONS OF LIQUID WATER
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUD AMOUNTS
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)
!             MINIMUM DIMENSIONS OF CONDENSED COMPONENTS
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSIONS OF CONDENSED COMPONENTS
!
!     PROPERTIES OF AEROSOLS:
      REAL
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_AEROSOL_SPECIES_LW)
!             MIXING RATIOS OF AEROSOLS
!
!     COUPLING FIELDS:
      INTEGER
     &     N_FRAC_ICE_POINT
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER
!
!     FLUXES:
      REAL
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_NET(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD/NET FLUX
     &   , FLUX_NET_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DOWNWARD/NET FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY UPWARD FLUX
!
!     FIELDS REQUIRED FOR CALL TO RADIATION CODE BUT NOT USED
      INTEGER
     &     N_ORDER_GAUSS
     &   , I_GAS
      LOGICAL
     &     L_SWITCH_SCATTER(NPD_BAND_LW)
      REAL
     &     SEC_0(NPD_PROFILE)
     &   , SOLAR_CONSTANT(NPD_PROFILE)
!
!
!     AUXILIARY VARIABLES:
      REAL
     &     CPBYG
!             SPECIFIC HEAT BY GRAVITY
     &   , DACON
!             DIFFERENCE IN A's
     &   , DBCON
!             DIFFERENCE IN B's
     &   , WEIGHT_BAND(NPD_BAND_LW)
!             WEIGHTING FACTORS FOR BANDS
     &   , NULLMMR
!             NULL MASS MIXING RATIO
      PARAMETER(CPBYG=CP/G)
      PARAMETER(NULLMMR=0.0E+00)
!
!     DUMMY FIELDS FOR RADIATION CODE
      LOGICAL
     &     L_DUMMY
      REAL
     &     DUMMY
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     R2_SET_GAS_MIX_RATIO, R2_SET_THERMODYNAMIC
     &   , R2_SET_AEROSOL_FIELD, R2_SET_CLOUD_FIELD
     &   , R2_SET_CLOUD_PARAMETRIZATION
     &   , R2_SET_SURFACE_FIELD_LW, R2_ZERO_1D
     &   , R2_COMPARE_PROC
!
!
!
!
!     INITIALIZE THE ERROR FLAG FOR THE RADIATION CODE.
      IERR=I_NORMAL
!     SET THE LOGICAL FLAG FOR DUMMY DIAGNOSTICS NOT AVAILABLE FROM
!     THE LOWER CODE IN THE LONG-WAVE TO .FALSE..
      L_DUMMY=.FALSE.
!
!
!     COMPARE PROCESSES IN THE SPECTRAL FILE WITH THOSE ENABLED IN
!     THE CODE.
      CALL R2_COMPARE_PROC(IERR, L_PRESENT_LW
     &   , L_RAYLEIGH_LW, L_GAS_LW, L_CONTINUUM_LW
     &   , L_DROP_LW, L_AEROSOL_LW, L_AEROSOL_CCN_LW, L_ICE_LW
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT
     &   , L_USE_SOOT_DIRECT
     &   , L_CLIMAT_AEROSOL
     &   , L_RAYLEIGH, L_GAS, L_CONTINUUM
     &   , L_DROP, L_AEROSOL, L_AEROSOL_CCN, L_ICE
     &   , NPD_TYPE_LW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!     CHECK THAT A VALID NUMBER HAS BEEN SUPPLIED FOR THE SOLVER.
      IF ( (I_SOLVER_LW.NE.IP_SOLVER_PENTADIAGONAL).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_MIX_11).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_MIX_APP_SCAT).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_MIX_DIRECT).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_HOMOGEN_DIRECT).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_TRIPLE).AND.
     &     (I_SOLVER_LW.NE.IP_SOLVER_TRIPLE_APP_SCAT)
     &   ) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: AN INVALID SOLVER HAS BEEN SELECTED '
     &      , 'IN THE LONGWAVE REGION.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!
!     THE GATHERING ARRAY IS REQUIRED BY THE SETTING SUBROUTINES (FOR
!     COMPATIBILITY WITH THE SHORTWAVE), BUT IS FILLED WITH INTEGERS
!     FROM 1 TO N_PROFILE SINCE ALL POINTS WILL BE CONSIDERED.
      DO L=1, N_PROFILE
         I_GATHER(L)=L
      ENDDO
!
!
!     SET THE MIXING RATIOS OF GASES.
      CALL R2_SET_GAS_MIX_RATIO(IERR
     &   , N_PROFILE, NLEVS, NWET, NOZONE
     &   , I_GATHER
     &   , N_ABSORB_LW, TYPE_ABSORB_LW
     &   , L_N2O_LW, L_CH4_LW, L_CFC11_LW, L_CFC12_LW,. FALSE.
     &   , L_CFC113_LW, L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW
     &   , H2O, CO2, O3, N2O_MIX_RATIO, CH4_MIX_RATIO
     &   , CFC11_MIX_RATIO, CFC12_MIX_RATIO, NULLMMR
     &   , CFC113_MIX_RATIO, HCFC22_MIX_RATIO, HFC125_MIX_RATIO
     &   , HFC134A_MIX_RATIO
     &   , GAS_MIX_RATIO
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_SPECIES_LW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!     CALCULATE PRESSURES AND TEMPERATURES.
      CALL R2_SET_THERMODYNAMIC(N_PROFILE, NLEVS, I_GATHER, .TRUE.
     &   , PSTAR, TSTAR, AB, BB, AC, BC, PEXNER, TAC
     &   , P, T, T_BDY, T_SURFACE, D_MASS
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &   )
!
!
!     SET THE MIXING RATIOS OF AEROSOLS.
      IF (L_AEROSOL.OR.L_AEROSOL_CCN) THEN
         CALL R2_SET_AEROSOL_FIELD(IERR
     &      , N_PROFILE, NLEVS, N_AEROSOL_LW, TYPE_AEROSOL_LW
     &      , I_GATHER
     &      , L_CLIMAT_AEROSOL, N_LEVELS_BL
     &      , L_USE_SULPC_DIRECT
     &      , SULP_DIM1, SULP_DIM2
     &      , ACCUM_SULPHATE, AITKEN_SULPHATE
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT
     &      , LAND, LYING_SNOW, PSTAR, AB, BB, TRINDX
     &      , AEROSOL_MIX_RATIO
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_LW
     &      )
      ENDIF
!
!
!     ASSIGN THE PROPERTIES OF CLOUDS. A DUMMY ARRAY MUST BE PASSED
!     FOR THE MICROPHYSICAL DIAGNOSTICS SINCE THEY ARE NOT AVAILABLE
!     THROUGH STASH IN THE LONG-WAVE.
!
      CALL R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND_LW
     &   , I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW, I_CNV_ICE_LW
     &   , L_DROP_TYPE_LW
     &   , I_DROP_PARAMETRIZATION_LW
     &   , DROP_PARAMETER_LIST_LW
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW
     &   , L_ICE_TYPE_LW
     &   , I_ICE_PARAMETRIZATION_LW
     &   , ICE_PARAMETER_LIST_LW
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , NPD_BAND_LW, NPD_DROP_TYPE_LW
     &   , NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
      CALL R2_SET_CLOUD_FIELD(N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , P, T, D_MASS
     &   , CCB, CCT, CCA, CCCWP
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK
     &   , L_MICROPHYSICS_LW, L_AEROSOL_CCN
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE
     &   , L_CLOUD_WATER_PARTITION,  LAND
     &   , I_CLOUD_REPRESENTATION_LW, I_CONDENSED_PARAM
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , W_CLOUD, FRAC_CLOUD, L_LOCAL_CNV_PARTITION_LW
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
!                       Microphysical Diagnostics are not available
!                       in this spectral region.
     &   , DUMMY, .FALSE., DUMMY, .FALSE.
     &   , DUMMY, .FALSE., DUMMY, .FALSE.
     &   , DUMMY, .FALSE.
     &   , DUMMY, .FALSE., DUMMY, .FALSE.
     &   , DUMMY, .FALSE., DUMMY, .FALSE.
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_LW
     &   , N_CCA_LEV, L_3D_CCA
     &   )
!
!
      CALL R2_SET_SURFACE_FIELD_LW(
     &     N_PROFILE, N_BAND_LW
     &   , I_SURFACE, I_SPEC_SURFACE_LW
     &   , L_SURFACE_LW
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF
     &   , ALBEDO_SEA_DIR, ALBEDO_SEA_DIFF
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , NPD_PROFILE, NPD_BAND_LW, NPD_SURFACE_LW
     &   )
!
!     SET CLEAR-SKY CALCULATIONS.
      L_CLEAR=L_CLEAR_OLR.OR.
     &        L_SURF_DOWN_CLR.OR.
     &        L_CLEAR_HR
!
      IF (L_CLEAR) THEN
!
!        SELECT A CLEAR-SKY SOLVER TO MATCH THE MAIN SOLVER.
         IF (I_SOLVER_LW.EQ.IP_SOLVER_PENTADIAGONAL) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_LW.EQ.IP_SOLVER_MIX_11) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_LW.EQ.IP_SOLVER_MIX_APP_SCAT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_LW.EQ.IP_SOLVER_MIX_DIRECT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_LW.EQ.IP_SOLVER_TRIPLE) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_LW.EQ.IP_SOLVER_TRIPLE_APP_SCAT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ENDIF
!
      ENDIF
!
!
!     SET PROPERTIES OF INDIVIDUAL BANDS.
      DO I=1, N_BAND_LW
         WEIGHT_BAND(I)=1.0E+00
         I_GAS_OVERLAP(I)=I_GAS_OVERLAP_LW
      ENDDO
!
!     INVERT THE TOPMOST CLOUDY LAYER IF USING A GLOBAL VALUE.
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP_GLOBAL=NLEVS+1-GLOBAL_CLOUD_TOP
      ENDIF
!
!
!
      CALL FLUX_CALC(IERR
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM
     &   , L_CLOUD_LW, L_DROP, L_ICE
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION_LW, I_2STREAM_LW, L_2_STREAM_CORRECT_LW
     &   , L_RESCALE_LW, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_LW, L_SWITCH_SCATTER
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL
!                       Options for Solver
     &   , I_SOLVER_LW
!                       General Spectral Properties
     &   , N_BAND_LW, 1, N_BAND_LW
     &   , WEIGHT_BAND
!                       General Atmospheric Properties
     &   , N_PROFILE, NLEVS
     &   , L_LAYER_LW, L_CLOUD_LAYER_LW
     &   , P, T, T_SURFACE, T_BDY, D_MASS
!                       Spectral Region
     &   , ISOLIR_LW
!                       Solar Fields
     &   , SEC_0, SOLAR_CONSTANT, SOLAR_FLUX_BAND_LW
     &   , RAYLEIGH_COEFFICIENT_LW
!                       Infra-red Fields
     &   , N_DEG_FIT_LW
     &   , THERMAL_COEFFICIENT_LW
     &   , T_REF_PLANCK_LW, L_IR_SOURCE_QUAD_LW
!                       Gaseous Absorption
     &   , N_ABSORB_LW, I_GAS_OVERLAP, I_GAS
     &   , GAS_MIX_RATIO
     &   , N_BAND_ABSORB_LW, INDEX_ABSORB_LW
     &   , I_BAND_ESFT_LW
     &   , W_ESFT_LW, K_ESFT_LW
     &   , I_SCALE_ESFT_LW, I_SCALE_FNC_LW
     &   , SCALE_VECTOR_LW
     &   , P_REFERENCE_LW, T_REFERENCE_LW
!                       Doppler Broadening
     &   , L_DOPPLER_PRESENT_LW
     &   , DOPPLER_CORRECTION_LW
!                       Surface Fields
     &   , L_SURFACE_LW, I_SURFACE
     &   , I_SPEC_SURFACE_LW
     &   , SURFACE_ALBEDO_LW
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR
     &   , N_DIR_ALBEDO_FIT_LW
     &   , DIRECT_ALBEDO_PARM_LW
     &   , EMISSIVITY_GROUND_LW
     &   , EMISSIVITY_FIELD
!                       Continuum Absorption
     &   , N_BAND_CONTINUUM_LW
     &   , INDEX_CONTINUUM_LW, INDEX_WATER_LW
     &   , K_CONTINUUM_LW
     &   , I_SCALE_FNC_CONT_LW
     &   , SCALE_CONTINUUM_LW
     &   , P_REF_CONTINUUM_LW
     &   , T_REF_CONTINUUM_LW
!                       Properties of Aerosols
     &   , N_AEROSOL_LW
     &   , AEROSOL_MIX_RATIO
     &   , AEROSOL_ABSORPTION_LW
     &   , AEROSOL_SCATTERING_LW
     &   , AEROSOL_ASYMMETRY_LW
     &   , I_AEROSOL_PARAMETRIZATION_LW
     &   , NHUMIDITY_LW
     &   , HUMIDITIES_LW
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED
     &   , I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW, W_CLOUD, FRAC_CLOUD
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
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, LWSEA
!                       Arrays for diagnostics specific to the UM
     &   , L_DUMMY, DUMMY, DUMMY
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &   , L_DUMMY, DUMMY
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &   , NPD_BAND_LW
     &   , NPD_SPECIES_LW
     &   , NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW, NPD_SCALE_VARIABLE_LW
     &   , NPD_CONTINUUM_LW
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW
     &   , NPD_CLOUD_PARAMETER_LW
     &   , NPD_THERMAL_COEFF_LW
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!
!     ASSIGNMENT OF DIAGNOSTICS:
!
!
!     OLR:
!
      DO L=1, N_PROFILE
         OLR(L)=-FLUX_NET(L, 0)
      ENDDO
      IF (L_CLEAR_OLR) THEN
         DO L=1, N_PROFILE
            CLEAR_OLR(L)=-FLUX_NET_CLEAR(L, 0)
         ENDDO
      ENDIF
!
!
!     TOTAL CLOUD COVER:
!
      IF (L_TOTAL_CLOUD_COVER) THEN
         CALL R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, NLEVS, NCLDS
     &      , I_CLOUD_LW, W_CLOUD, TOTAL_CLOUD_COVER
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!     NET FLUX AT THE TROPOPAUSE:
!
      IF (L_NET_FLUX_TROP) THEN
         DO L=1, N_PROFILE
            NET_FLUX_TROP(L)
     &         =FLUX_NET(L, NLEVS+1-TRINDX(L))
         ENDDO
      ENDIF
!
!
!     DOWNWARD FLUX AT THE TROPOPAUSE:
!
      IF (L_DOWN_FLUX_TROP) THEN
         DO L=1, N_PROFILE
            DOWN_FLUX_TROP(L)
     &         =FLUX_NET(L, NLEVS+1-TRINDX(L))
     &         +FLUX_UP(L, NLEVS+1-TRINDX(L))
         ENDDO
      ENDIF
!
!
!
!
!
!     OUTPUT ARRAYS:
!
!     CONVERT THE FLUXES TO INCREMENTS IN THE HEATING RATE EXCEPT AT
!     THE SURFACE: THERE, THE NET DOWNWARD FLUX IS ASSIGNED TO LWOUT.
      DO I=NLEVS, 1, -1
         DACON=(AB(I)-AB(I+1))*CPBYG/PTS
         DBCON=(BB(I)-BB(I+1))*CPBYG/PTS
         DO L=1, N_PROFILE
            LWOUT(L, I+1)=(FLUX_NET(L, NLEVS-I)
     &         -FLUX_NET(L, NLEVS+1-I))/(DACON+PSTAR(L)*DBCON)
         ENDDO
         IF (L_CLEAR_HR) THEN
!           THE FACTOR OF PTS IS INCLUDED HERE TO YIELD A RATE FROM AN
!           INCREMENT.
            DO L=1, N_PROFILE
               CLEAR_HR(L, I)=(FLUX_NET_CLEAR(L, NLEVS-I)
     &            -FLUX_NET_CLEAR(L, NLEVS+1-I))/(PTS
     &            *(DACON+PSTAR(L)*DBCON))
            ENDDO
         ENDIF
      ENDDO

      DO L=1, N_PROFILE
         LWOUT(L, 1)=FLUX_NET(L, NLEVS)
      ENDDO
!
!     SEPARATE THE CONTRIBUTIONS OVER OPEN SEA AND SEA-ICE.
!     LWSEA MUST BE WEIGHTED WITH THE FRACTION OF OPEN SEA.
      DO L=1, N_PROFILE
         IF (LAND(L)) THEN
            LWSEA(L)=0.0
         ELSE IF (ICE_FRACTION(L).LT.TOL_TEST) THEN
            LWSEA(L)=LWOUT(L, 1)
            LWOUT(L, 1)=0.0
         ELSE
!           LWSEA MUST BE SCALED BY THE FRACTION OF OPEN SEA FOR
!           CONSISTENCY WITH UPPER LEVELS IN THE MODEL.
            LWSEA(L)=(1.0E+00-ICE_FRACTION(L))*LWSEA(L)
            LWOUT(L, 1)=LWOUT(L, 1)-LWSEA(L)
         ENDIF
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
      SUBROUTINE R2_SET_SURFACE_FIELD_LW(
     &     N_PROFILE, N_BAND
     &   , I_SURFACE, I_SPEC_SURFACE, L_SURFACE
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , NPD_PROFILE, NPD_BAND_LW, NPD_SURFACE_LW
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
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES:
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
     &   , NPD_BAND_LW
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_SURFACE_LW
!             MAXIMUM NUMBER OF SURFACES
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , N_BAND
!             NUMBER OF SPECTRAL BANDS
!
!     PROPERTIES OF SURFACES
      INTEGER   !, INTENT(OUT)
     &     I_SURFACE(NPD_PROFILE)
!             TYPES OF SURFACES
     &   , I_SPEC_SURFACE(NPD_SURFACE_LW)
      LOGICAL   !, INTENT(OUT)
     &     L_SURFACE(NPD_SURFACE_LW)
!             FLAGS FOR TYPES OF SURFACES
!
!     SURFACE PROPERTIES SET.
      REAL      !, INTENT(OUT)
     &     EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_LW)
!             EMISSIVITIES OF SURFACES
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_LW)
!             DIFFUSE ALBEDO OF SURFACE
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_LW)
!             DIRECT ALBEDO OF SURFACE
     &   , ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND_LW)
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND_LW)
!             DIRECT ALBEDO OF OPEN SEA
!
!     VARIABLES CONCERNED WITH FRACTIONAL SEA ICE
      REAL      !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!
      INTEGER   !, INTENT(OUT)
     &     N_FRAC_ICE_POINT
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     SEARCH_ARRAY(NPD_PROFILE)
!             ARRAY FOR SEARCHING
     &   , TARGET
!             TARGET TO SEARCH FOR
!
!
!
!     OVERRIDE ANY SURFACE PROERTIES READ IN FROM THE SPECTRAL FILE.
      DO L=1, N_PROFILE
         I_SURFACE(L)=1
      ENDDO
      L_SURFACE(1)=.TRUE.
      I_SPEC_SURFACE(1)=IP_SURFACE_INTERNAL
!
!     SET THE EMISSIVITY FIELD.
      DO I=1, N_BAND
         DO L=1, N_PROFILE
            EMISSIVITY_FIELD(L, I)=1.0E+00
         ENDDO
      ENDDO
!
!     ZERO THE SURFACE ALBEDOS.
      DO I=1, N_BAND
         DO L=1, N_PROFILE
            ALBEDO_FIELD_DIFF(L, I)=0.0E+00
            ALBEDO_FIELD_DIR(L, I)=0.0E+00
            ALBEDO_SEA_DIFF(L, I)=0.0E+00
            ALBEDO_SEA_DIR(L, I)=0.0E+00
         ENDDO
      ENDDO
!
!     SET THE FRACTIONAL ICE COVERAGE. POINTS ARE REQUIRED WHERE
!     THE ICE FRACTION IS NEITHER 0 NOR 1.
      DO L=1, N_PROFILE
         SEARCH_ARRAY(L)=ICE_FRACTION(L)*(1.0E+00-ICE_FRACTION(L))
      ENDDO
      TARGET=TOL_TEST
!                                                                       
      N_FRAC_ICE_POINT=0
      DO L   =1,N_PROFILE
        IF (SEARCH_ARRAY(L).GT.TARGET) THEN
          N_FRAC_ICE_POINT                  =N_FRAC_ICE_POINT+1
          I_FRAC_ICE_POINT(N_FRAC_ICE_POINT)=L
        END IF
      END DO
!
!
!
      RETURN
      END
