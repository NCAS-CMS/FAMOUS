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
!+ Subroutine to read a shortwave spectral namelist.
!
! Purpose:
!   To read a shortwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-05-96                Set lower limits
!                                               for reduced dimensions
!                                               to ensure that they
!                                               may never be 0.
!                                               (J. M. Edwards)
!       4.4             02-09-97                Aerosol flags passed
!                                               in to the code to
!                                               enable only those
!                                               required to be
!                                               selected. Spectral
!                                               data are now longer
!                                               compressed into a
!                                               single array.
!                                               Actual IOS code put
!                                               into CMESSAGE.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!
!       4.5        April 1998   Allow soot spectral data to be read.
!                                                     Luke Robinson.
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SW_SPECIN(IERR, CMESSAGE
     &   , L_O2
     &   , L_CLIMAT_AEROSOL, L_USE_SULPC_DIRECT
     &  , L_USE_SOOT_DIRECT
     &   )
!
!
      IMPLICIT NONE
!
!
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
      INTEGER
     &     NPD_TYPE
!             NUMBER OF TYPES OF DATA
     &   , NPD_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , NPD_EXCLUDE
!             NUMER OF EXCLUDED BANDS
     &   , NPD_SPECIES
!             NUMBER OF GASEOUS SPECIES
     &   , NPD_ESFT_TERM
!             NUMBER OF ESFT TERMS
     &   , NPD_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             NUMBER OF SCALING VARIABLES
     &   , NPD_SURFACE
!             NUMBER OF SURFACE TYPES
     &   , NPD_ALBEDO_PARM
!             NUMBER OF ALBEDO PARAMETERS
     &   , NPD_CONTINUUM
!             NUMBER OF CONTINUA
     &   , NPD_DROP_TYPE
!             NUMBER OF DROP TYPES
     &   , NPD_ICE_TYPE
!             NUMBER OF ICE CRYSTAL TYPES
     &   , NPD_AEROSOL_SPECIES
!             NUMBER OF AEROSOL SPECIES
     &   , NPD_CLOUD_PARAMETER
!             MAX NUMBER OF CLOUD PARAMETERS
     &   , NPD_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
!
      PARAMETER(
     &     NPD_TYPE=15
     &   , NPD_BAND=20
     &   , NPD_EXCLUDE=2
     &   , NPD_SPECIES=11
     &   , NPD_ESFT_TERM=16
     &   , NPD_SCALE_FNC=3
     &   , NPD_SCALE_VARIABLE=4
     &   , NPD_THERMAL_COEFF=9
     &   , NPD_SURFACE=1
     &   , NPD_ALBEDO_PARM=4
     &   , NPD_CONTINUUM=2
     &   , NPD_DROP_TYPE=6
     &   , NPD_ICE_TYPE=7
     &   , NPD_CLOUD_PARAMETER=30
     &   , NPD_AEROSOL_SPECIES=9
     &   , NPD_HUMIDITIES=21
     &   )
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
      LOGICAL   !, INTENT(IN)
     &     L_O2
!             ABSORPTION BY OXYGEN IS TO BE INCLUDED.
     &   , L_CLIMAT_AEROSOL
!             CLIMATOLOGICAL AEROSOLS ARE TO BE INCLUDED
     &   , L_USE_SULPC_DIRECT
!             THE DIRECT EFFECTS OF SULPHATE AEROSOLS ARE
!             TO BE INCLUDED
     &   , L_USE_SOOT_DIRECT
!             USE THE DIRECT RAD EFFECTS OF SOOT IN THE SW
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      CHARACTER*80      !, INTENT(OUT)
     &     CMESSAGE
!
!
!
!     LOCAL VARIABLES.
!
!
!     RADIATIVE VARIABLES FOR REDUCING THE SPECTRUM
!
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
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
      INTEGER
     &     NPD_AEROSOL_COMPONENT
!             MAXIMUM NUMBER OF AEROSOL COMPONENTS
      INTEGER
     &     IP_WATER_SOLUBLE
!             WATER SOLUBLE AEROSOL
     &   , IP_DUST_LIKE
!             DUST-LIKE AEROSOL
     &   , IP_OCEANIC
!             OCEANIC AEROSOL
     &   , IP_SOOT
!             SOOT AEROSOL
     &   , IP_ASH
!             VOLCANIC ASH
     &   , IP_SULPHURIC
!             SULPHURIC ACID
     &   , IP_ACCUM_SULPHATE
!             ACCUMULATION MODE SULPHATE
     &   , IP_AITKEN_SULPHATE
!             AITKEN MODE SULPHATE
     &   , IP_FRESH_SOOT
     &   , IP_AGED_SOOT
!
!
      PARAMETER(
     &     NPD_AEROSOL_COMPONENT=13
     &   )
      PARAMETER(
     &     IP_WATER_SOLUBLE=1
     &   , IP_DUST_LIKE=2
     &   , IP_OCEANIC=3
     &   , IP_SOOT=4
     &   , IP_ASH=5
     &   , IP_SULPHURIC=6
     &   , IP_ACCUM_SULPHATE=10
     &   , IP_AITKEN_SULPHATE=11
     &   , IP_FRESH_SOOT=12
     &   , IP_AGED_SOOT=13

     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INDEXING NUMBERS OF GASEOUS ABSORBING SPECIES.
!     THE NUMBERING 1-12 CORRESPONDS TO LOWTRAN 7.
!
      INTEGER
     &     NPD_GASES
!             NUMBER OF INDEXED GASES
      PARAMETER(NPD_GASES=19)
!
      INTEGER
     &     IP_H2O
!             INDEX NUMBER OF WATER VAPOUR
     &   , IP_CO2
!             INDEX NUMBER OF CARBON DIOXIDE
     &   , IP_O3
!             INDEX NUMBER OF OZONE
     &   , IP_N2O
!             INDEX NUMBER OF DINITROGEN OXIDE
     &   , IP_CO
!             INDEX NUMBER OF CARBON MONOXIDE
     &   , IP_CH4
!             INDEX NUMBER OF METHANE
     &   , IP_O2
!             INDEX NUMBER OF OXYGEN
     &   , IP_NO
!             INDEX NUMBER OF NITROGEN MONOXIDE
     &   , IP_SO2
!             INDEX NUMBER OF SULPHUR DIOXIDE
     &   , IP_NO2
!             INDEX NUMBER OF NITROGEN DIOXIDE
     &   , IP_NH3
!             INDEX NUMBER OF AMMONIA
     &   , IP_HNO3
!             INDEX NUMBER OF NITRIC ACID
     &   , IP_N2
!             INDEX NUMBER OF NITROGEN
     &   , IP_CFC11
!             INDEX NUMBER OF CFC11
     &   , IP_CFC12
!             INDEX NUMBER OF CFC12
     &   , IP_CFC113
!             INDEX NUMBER OF CFC113
     &   , IP_HCFC22
!             INDEX NUMBER OF HCFC22
     &   , IP_HFC125
!             INDEX NUMBER OF HFC125
     &   , IP_HFC134A
!             INDEX NUMBER OF HCF134A
!
      PARAMETER(
     &     IP_H2O=1
     &   , IP_CO2=2
     &   , IP_O3=3
     &   , IP_N2O=4
     &   , IP_CO=5
     &   , IP_CH4=6
     &   , IP_O2=7
     &   , IP_NO=8
     &   , IP_SO2=9
     &   , IP_NO2=10
     &   , IP_NH3=11
     &   , IP_HNO3=12
     &   , IP_N2=13
     &   , IP_CFC11=14
     &   , IP_CFC12=15
     &   , IP_CFC113=16
     &   , IP_HCFC22=17
     &   , IP_HFC125=18
     &   , IP_HFC134A=19
     &   )
!
!     ------------------------------------------------------------------
!
      CHARACTER*80
     &     SW_SPECTRAL_FILE
!             NAME OF FILE CONTAINING THE SPECTRAL DATA
      INTEGER
     &     IERR_GET_FILE
!             ERROR FLAG RETURNED BY GET_FILE (NOT NECESSARILY
!             CONSISTENT WITH THE FLAGS IN ERROR3A).
     &   , IOS
!             STATUS OF I/O
!
      LOGICAL
     &     L_RETAIN_ABSORB(NPD_SPECIES)
!             FLAG SET TO .TRUE. IF THE ABSORBER IS TO BE RETAINED
     &   , L_GAS_INCLUDED(NPD_GASES)
!             LOGICAL TO TEST FOR ACTUAL GASES INCLUDED
      INTEGER
     &     N_ABSORB_RETAIN
!             NUMBER OF ABSORBERS TO RETAIN
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , COMPRESSED_INDEX(NPD_SPECIES)
!             MAPPING FROM ORIGINAL TO COMPRESSED INDICES OF ABSORBERS
     &   , N_AEROSOL_RETAIN
!             NUMBER OF AEROSOLS IN THE SPECTRAL FILE TO BE RETAINED
!             FOR THE RADIATIVE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)
!             INDEXING NUMBERS OF THE RETAINED AEROSOLS
     &   , N_AEROSOL_FOUND
!             NUMBER OF AEROSOLS FOR THE CURRENT GROUP OF PROCESSES
!             FOUND IN THE SPECTRAL FILE
!
!
!
!     DECLARE THE ELEMENTS OF THE INITIAL SPECTRUM FOR DYNAMIC
!     ALLOCATION AND SET UP AN APPROPRIATE NAMELIST.
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR SPECTRAL FILE.
!     ------------------------------------------------------------------
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT(0: NPD_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT(NPD_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE(NPD_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB(NPD_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB(NPD_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)
!             NUMBER OF ESFT TERMS IN BAND FOR EACH GAS
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND
     &        , NPD_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE(NPD_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM(NPD_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             LIST OF CONTINUA CONTINUUA IN EACH BAND
     &   , INDEX_WATER
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM(NPD_SCALE_VARIABLE, NPD_BAND, NPD_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER CORRECTION TERMS
!
!
!
!    ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE TO DECLARE ELEMENTS OF A SPECTRAL FILE AS A NAMELIST.
!
      NAMELIST/R2SWSP/
!                       Blocks Present
     &     L_PRESENT
!                       Block 0
     &   , N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL
!                       Block 1
     &   , WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG
!                       Block 2
     &   , SOLAR_FLUX_BAND
!                       Block 3
     &   , RAYLEIGH_COEFFICIENT
!                       Block 4
     &   , N_BAND_ABSORB, INDEX_ABSORB
!                       Block 5
     &   , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &   , P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR
!                       Block 6
     &   , N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT
!                       Block 7
     &   , I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE
     &   , SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND
!                       Block 8
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
!                       Block 9
     &   , I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM
     &   , K_CONTINUUM, SCALE_CONTINUUM
!                       Block 10
     &   , I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
!                       Block 11
     &   , I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , HUMIDITIES
!                       Block 12
     &   , I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
!                       Block 13
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION
!                       Block 14
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
!
!     ------------------------------------------------------------------
!
!
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED SW-SPECTRAL FILE.
!     NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT
!     NOTE: SINCE THE ARRAYS HERE WILL BE PASSED IN A COMMON BLOCK
!     THEIR SIZES MUST BE FIXED, EVEN THOUGH VARIABLE SIZES ARE USED
!     LOWER IN THE CODE. THEY ARE ACCORDINGLY DEFINED AS 1-DIMENSIONAL
!     ARRAYS WITH FIXED MAXIMUM SIZES AT THIS LEVEL.
!
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
     &     L_PRESENT_SW(0: NPD_TYPE)
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
     &     WAVE_LENGTH_SHORT_SW(NPD_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_SW(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE_SW(NPD_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_SW(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND_SW(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT_SW(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB_SW
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_SW(NPD_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_SW(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_SW(NPD_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_SW(NPD_BAND, NPD_SPECIES)
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_SW(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_SW(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_SW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM
     &        , NPD_BAND, NPD_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_SW(NPD_SPECIES, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_SW(NPD_SPECIES, NPD_BAND)
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
     &     THERMAL_COEFFICIENT_SW(0: NPD_THERMAL_COEFF-1
     &        , NPD_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_SW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE_SW(NPD_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_SW(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE_SW(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO_SW(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_SW(0: NPD_ALBEDO_PARM
     &      , NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_SW(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM_SW(NPD_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_SW
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_SW(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_SW(NPD_SCALE_VARIABLE
     &      , NPD_BAND, NPD_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION_SW(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE_SW(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_SW(NPD_DROP_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_SW(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL_SW
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_SW(NPD_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_SW(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_SW(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES_SW(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION_SW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_SW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_SW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION_SW(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE_SW(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_SW(NPD_ICE_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_SW(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT_SW(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION_SW(NPD_SPECIES)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!    ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE DECLARING COMMON BLOCK CONTAINING THE REDUCED SW SPECTRAL
!     FILE.
!     (NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT)
!
!
      COMMON/R2SWSPCM/
!
!     DIMENSIONS OF ARRAYS:
     &     NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW
     &   , NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW
     &   , NPD_SCALE_VARIABLE_SW
     &   , NPD_THERMAL_COEFF_SW
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW
     &   , NPD_CONTINUUM_SW
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW
     &   , NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW
!
!     GENERAL ARRAYS:
     &   , L_PRESENT_SW
!
!     PROPERTIES OF BANDS:
     &   , N_BAND_SW, WAVE_LENGTH_SHORT_SW, WAVE_LENGTH_LONG_SW
!
!     EXCLUSIONS FROM BANDS:
     &   , N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW
!
!     SOLAR FIELDS:
     &   , SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW
!
!     GASEOUS ABSORPTION:
     &   , N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW
     &   , TYPE_ABSORB_SW
     &   , I_BAND_ESFT_SW, I_SCALE_ESFT_SW, I_SCALE_FNC_SW
     &   , K_ESFT_SW, W_ESFT_SW
     &   , SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW
!
!     THERMAL SOURCE FUNCTION:
     &   , N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW
!
!     SURFACE PROPERTIES:
     &   , I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW
     &   , SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW
     &   , EMISSIVITY_GROUND_SW
!
!     CONTINUA:
     &   , N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW, INDEX_WATER_SW
     &   , I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW
     &   , SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW
!
!     WATER DROPLETS:
     &   , I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW
     &   , DROP_PARAMETER_LIST_SW
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW
!
!     AEROSOLS:
     &   , N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW
     &   , NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW
     &   , AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW
     &   , AEROSOL_ASYMMETRY_SW
!
!     ICE CRYSTALS:
     &   , I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW
     &   , ICE_PARAMETER_LIST_SW
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW
!
!     DOPPLER BROADENING:
     &   , L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW
!
!     ------------------------------------------------------------------
!
!
!
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
!
      CHARACTER
     &     CH_IOS*5
!             CHARACTER STRING FOR IOS ERROR
!
!
!     SUBROUTINES CALLED
      EXTERNAL
     &     R2_COMPRESS_SPECTRUM
!
!
!     EACH BLOCK IS INITIALIZED AS MISSING:
      DATA L_PRESENT/.FALSE., NPD_TYPE*.FALSE./
!
!     INITIALIZE THE RANGE OF VALIDITY OF THE PARAMETRIZATIONS OF
!     DROPLETS AND ICE CRYSTALS. OLD SPECTRAL FILES WILL NOT CONTAIN
!     SUCH DATA, SO THE LIMITS FOR DROPLETS ARE INITIALIZED TO THOSE
!     FORMERLY SET IN THE MICROPHYSICAL SCHEME (MRF/UMIST
!     PARAMETRIZATION) TO ENSURE THAT THE RESULTS ARE BIT-REPRODUCIBLE.
!     VALUES FOR ICE COVER THE RANGE OF EFFECTIVE RADII USED IN
!     GENERATING THE DATA FOR THE ORIGINAL PARAMETRIZATION OF ICE
!     CRYSTALS.
!     AT SOME FUTURE RELEASE IT MAY BE DESIRABLE TO REMOVE DEFAULT
!     SETTINGS.
      DATA DROP_PARM_MIN_DIM/NPD_DROP_TYPE*3.5E-07/
      DATA DROP_PARM_MAX_DIM/NPD_DROP_TYPE*3.7E-05/
      DATA ICE_PARM_MIN_DIM/NPD_ICE_TYPE*3.75E-07/
      DATA ICE_PARM_MAX_DIM/NPD_ICE_TYPE*8.0E-05/
!
!
!
!     READ THE SHORTWAVE SPECTRUM AS A NAMELIST.
      CALL GET_FILE(57, SW_SPECTRAL_FILE, 80, IERR_GET_FILE)
      IF (IERR_GET_FILE.NE.0) THEN
!        CONVERT THE ERROR FLAG FROM GET_FILE TO A FLAG RECOGNISED
!        BY THE RADIATION CODE.
         IERR=I_ERR_IO
         CMESSAGE='Error reading name of shortwave spectral file.'
         RETURN
      ENDIF
      OPEN(UNIT=57, FILE=SW_SPECTRAL_FILE, IOSTAT=IOS,
     & DELIM='APOSTROPHE')
      IF (IOS.NE.0) THEN
         IERR=I_ERR_IO
      WRITE(CH_IOS, '(I5)') IOS
         CMESSAGE='Error opening shortwave spectral file.'
     &      //' IOSTAT='//CH_IOS
         RETURN
      ENDIF
      READ(57, R2SWSP)
      CLOSE(57)
!
!     TEST FOR MINIMAL REQUISITE INFORMATION.
      IF ( .NOT.(L_PRESENT(0).AND.
     &           L_PRESENT(2) ) ) THEN
         CMESSAGE='Shortwave spectrum is deficient.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET REDUCED DIMENSIONS, EITHER FROM THE SIZES OF THE FIXED ARRAYS
!     OR FROM THE ARRAYS READ IN.
!
      NPD_TYPE_SW=NPD_TYPE
      NPD_BAND_SW=MAX(N_BAND, 1)
      NPD_SPECIES_SW=MAX(N_ABSORB, 1)
      NPD_ALBEDO_PARM_SW=NPD_ALBEDO_PARM
      NPD_SCALE_FNC_SW=NPD_SCALE_FNC
      NPD_SCALE_VARIABLE_SW=NPD_SCALE_VARIABLE
      NPD_SURFACE_SW=NPD_SURFACE
      NPD_CONTINUUM_SW=NPD_CONTINUUM
      NPD_CLOUD_PARAMETER_SW=NPD_CLOUD_PARAMETER
      NPD_THERMAL_COEFF_SW=1
!
!
!     SEARCH THE SPECTRUM TO FIND MAXIMUM DIMENSIONS.
!
      NPD_EXCLUDE_SW=1
      IF (L_PRESENT(14)) THEN
         DO I=1, N_BAND
            NPD_EXCLUDE_SW=MAX(NPD_EXCLUDE_SW, N_BAND_EXCLUDE(I))
         ENDDO
      ENDIF
!
!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      DO I=1, NPD_GASES
         L_GAS_INCLUDED(I)=.FALSE.
      ENDDO
      N_ABSORB_RETAIN=0
!
      DO I=1, N_ABSORB
!
         L_RETAIN_ABSORB(I)=.FALSE.
         COMPRESSED_INDEX(I)=0
!
         IF ( (TYPE_ABSORB(I).EQ.IP_H2O).OR.
     &        (TYPE_ABSORB(I).EQ.IP_CO2).OR.
     &        (TYPE_ABSORB(I).EQ.IP_O3).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_O2).AND.L_O2 ) ) THEN
            N_ABSORB_RETAIN=N_ABSORB_RETAIN+1
            INDEX_ABSORB_RETAIN(N_ABSORB_RETAIN)=I
            COMPRESSED_INDEX(I)=N_ABSORB_RETAIN
            L_RETAIN_ABSORB(I)=.TRUE.
            L_GAS_INCLUDED(TYPE_ABSORB(I))=.TRUE.
         ENDIF
!
      ENDDO
!
!
!     Print warning messages if those gases normally expected
!     are not present.
      IF (.NOT.L_GAS_INCLUDED(IP_H2O)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Water vapour is not included in the '
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_CO2)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Carbon dioxide is not included in the '
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_O3)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Ozone is not included in the '
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_O2)).AND.L_O2) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: Oxygen is not included in the shortwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_SPECIES_SW=MAX(N_ABSORB_RETAIN, 1)
!
!
      NPD_ESFT_TERM_SW=1
      IF (L_PRESENT(5)) THEN
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I)))
     &            NPD_ESFT_TERM_SW=MAX(NPD_ESFT_TERM_SW
     &            , I_BAND_ESFT(I, INDEX_ABSORB(J, I)))
            ENDDO
         ENDDO
      ENDIF
!
      NPD_DROP_TYPE_SW=1
      IF (L_PRESENT(10)) THEN
         DO I=1, NPD_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               NPD_DROP_TYPE_SW=MAX(NPD_DROP_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
      NPD_ICE_TYPE_SW=1
      IF (L_PRESENT(12)) THEN
         DO I=1, NPD_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               NPD_ICE_TYPE_SW=MAX(NPD_ICE_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      NPD_HUMIDITIES_SW=1
      N_AEROSOL_RETAIN=0
!
!     Check the spectral file for climatological aerosols
      IF (L_CLIMAT_AEROSOL) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I).EQ.IP_WATER_SOLUBLE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_DUST_LIKE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_OCEANIC).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_SOOT).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_SULPHURIC) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.5) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '
     &            //'climatological aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF

!
!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).
!
      IF (L_USE_SULPC_DIRECT) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I).EQ.IP_ACCUM_SULPHATE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_AITKEN_SULPHATE) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '
     &            //'sulphate aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for soot aerosol modes. (Also only
!     required for the direct effect).
!
      IF (L_USE_SOOT_DIRECT) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I).EQ.IP_FRESH_SOOT).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_AGED_SOOT) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '
     &            //'soot aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_AEROSOL_SPECIES_SW=MAX(N_AEROSOL_RETAIN, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      IF (L_PRESENT(11)) THEN
         DO I=1, N_AEROSOL_RETAIN
            IF (I_AEROSOL_PARAMETRIZATION(INDEX_AEROSOL_RETAIN(I)).EQ.
     &         IP_AEROSOL_PARAM_MOIST) THEN
               NPD_HUMIDITIES_SW=MAX(NPD_HUMIDITIES_SW
     &            , NHUMIDITY(INDEX_AEROSOL_RETAIN(I)))
            ENDIF
         ENDDO
      ENDIF
!
!
!
!
!     TRANSFER THE LARGE NAMELIST TO THE REDUCED SPECTRUM.
!
!
      CALL R2_COMPRESS_SPECTRUM(
!                       Spectral Array in Namelist
     &     L_PRESENT
     &   , N_BAND, WAVE_LENGTH_SHORT , WAVE_LENGTH_LONG
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
     &   , SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT
     &   , N_ABSORB, N_BAND_ABSORB, INDEX_ABSORB, TYPE_ABSORB
     &   , L_RETAIN_ABSORB, N_ABSORB_RETAIN, INDEX_ABSORB_RETAIN
     &   , COMPRESSED_INDEX, I_BAND_ESFT, K_ESFT, W_ESFT, I_SCALE_ESFT
     &   , I_SCALE_FNC, SCALE_VECTOR, P_REFERENCE, T_REFERENCE
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK
     &   , I_SPEC_SURFACE, L_SURFACE, SURFACE_ALBEDO
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
     &   , N_AEROSOL, TYPE_AEROSOL
     &   , N_AEROSOL_RETAIN, INDEX_AEROSOL_RETAIN
     &   , L_AEROSOL_SPECIES, AEROSOL_ABSORPTION
     &   , AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , NHUMIDITY, HUMIDITIES, I_AEROSOL_PARAMETRIZATION
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION
!                       Reduced Spectral Array
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
     &   )
!
!
!
      RETURN
      END
!+ Subroutine to read a longwave spectral namelist.
!
! Purpose:
!   To read a longwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
!       4.4             02-09-97                Aerosol flags passed
!                                               in to the code to
!                                               enable only those
!                                               required to be
!                                               selected. Spectral
!                                               data are no longer
!                                               compressed into a
!                                               single array.
!                                               IOSTAT error code
!                                               returned as part of
!                                               CMESSAGE.
!                                               (J. M. Edwards)
!       4.5        April 1998   Allow soot spectral data to be read.
!                                                     Luke Robinson.
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_LW_SPECIN(IERR, CMESSAGE
     &   , L_CH4, L_N2O, L_CFC11, L_CFC12
     &   , L_CFC113, L_HCFC22, L_HFC125, L_HFC134A
     &   , L_CLIMAT_AEROSOL, L_USE_SULPC_DIRECT
     &  , L_USE_SOOT_DIRECT
     &   )
!
!
      IMPLICIT NONE
!
!
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
      INTEGER
     &     NPD_TYPE
!             NUMBER OF TYPES OF DATA
     &   , NPD_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , NPD_EXCLUDE
!             NUMER OF EXCLUDED BANDS
     &   , NPD_SPECIES
!             NUMBER OF GASEOUS SPECIES
     &   , NPD_ESFT_TERM
!             NUMBER OF ESFT TERMS
     &   , NPD_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             NUMBER OF SCALING VARIABLES
     &   , NPD_SURFACE
!             NUMBER OF SURFACE TYPES
     &   , NPD_ALBEDO_PARM
!             NUMBER OF ALBEDO PARAMETERS
     &   , NPD_CONTINUUM
!             NUMBER OF CONTINUA
     &   , NPD_DROP_TYPE
!             NUMBER OF DROP TYPES
     &   , NPD_ICE_TYPE
!             NUMBER OF ICE CRYSTAL TYPES
     &   , NPD_AEROSOL_SPECIES
!             NUMBER OF AEROSOL SPECIES
     &   , NPD_CLOUD_PARAMETER
!             MAX NUMBER OF CLOUD PARAMETERS
     &   , NPD_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
!
      PARAMETER(
     &     NPD_TYPE=15
     &   , NPD_BAND=20
     &   , NPD_EXCLUDE=2
     &   , NPD_SPECIES=11
     &   , NPD_ESFT_TERM=16
     &   , NPD_SCALE_FNC=3
     &   , NPD_SCALE_VARIABLE=4
     &   , NPD_THERMAL_COEFF=9
     &   , NPD_SURFACE=1
     &   , NPD_ALBEDO_PARM=4
     &   , NPD_CONTINUUM=2
     &   , NPD_DROP_TYPE=6
     &   , NPD_ICE_TYPE=7
     &   , NPD_CLOUD_PARAMETER=30
     &   , NPD_AEROSOL_SPECIES=9
     &   , NPD_HUMIDITIES=21
     &   )
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
      LOGICAL   !, INTENT(IN)
     &     L_CH4
!             ABSORPTION BY METHANE IS INCLUDED
     &   , L_N2O
!             ABSORPTION BY NITROUS OXIDE IS INCLUDED
     &   , L_CFC11
!             ABSORPTION BY CFC11 IS INCLUDED
     &   , L_CFC12
!             ABSORPTION BY CFC12 IS INCLUDED
     &   , L_CFC113
!             ABSORPTION BY CFC113 IS INCLUDED
     &   , L_HCFC22
!             ABSORPTION BY HCFC22 IS INCLUDED
     &   , L_HFC125
!             ABSORPTION BY HFC125 IS INCLUDED
     &   , L_HFC134A
!             ABSORPTION BY HFC134A IS INCLUDED
     &   , L_CLIMAT_AEROSOL
!             CLIMATOLOGICAL AEROSOLS ARE TO BE INCLUDED
     &   , L_USE_SULPC_DIRECT
!             THE DIRECT EFFECTS OF SULPHATE AEROSOLS ARE
!             TO BE INCLUDED
     &   , L_USE_SOOT_DIRECT
!             USE THE DIRECT RAD EFFECTS OF SOOT IN THE LW
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      CHARACTER*80      !, INTENT(OUT)
     &     CMESSAGE
!
!
!
!     LOCAL VARIABLES.
!
!
!     RADIATIVE VARIABLES FOR REDUCING THE SPECTRUM
!
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
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
      INTEGER
     &     NPD_AEROSOL_COMPONENT
!             MAXIMUM NUMBER OF AEROSOL COMPONENTS
      INTEGER
     &     IP_WATER_SOLUBLE
!             WATER SOLUBLE AEROSOL
     &   , IP_DUST_LIKE
!             DUST-LIKE AEROSOL
     &   , IP_OCEANIC
!             OCEANIC AEROSOL
     &   , IP_SOOT
!             SOOT AEROSOL
     &   , IP_ASH
!             VOLCANIC ASH
     &   , IP_SULPHURIC
!             SULPHURIC ACID
     &   , IP_ACCUM_SULPHATE
!             ACCUMULATION MODE SULPHATE
     &   , IP_AITKEN_SULPHATE
!             AITKEN MODE SULPHATE
     &   , IP_FRESH_SOOT
     &   , IP_AGED_SOOT
!
!
      PARAMETER(
     &     NPD_AEROSOL_COMPONENT=13
     &   )
      PARAMETER(
     &     IP_WATER_SOLUBLE=1
     &   , IP_DUST_LIKE=2
     &   , IP_OCEANIC=3
     &   , IP_SOOT=4
     &   , IP_ASH=5
     &   , IP_SULPHURIC=6
     &   , IP_ACCUM_SULPHATE=10
     &   , IP_AITKEN_SULPHATE=11
     &   , IP_FRESH_SOOT=12
     &   , IP_AGED_SOOT=13

     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INDEXING NUMBERS OF GASEOUS ABSORBING SPECIES.
!     THE NUMBERING 1-12 CORRESPONDS TO LOWTRAN 7.
!
      INTEGER
     &     NPD_GASES
!             NUMBER OF INDEXED GASES
      PARAMETER(NPD_GASES=19)
!
      INTEGER
     &     IP_H2O
!             INDEX NUMBER OF WATER VAPOUR
     &   , IP_CO2
!             INDEX NUMBER OF CARBON DIOXIDE
     &   , IP_O3
!             INDEX NUMBER OF OZONE
     &   , IP_N2O
!             INDEX NUMBER OF DINITROGEN OXIDE
     &   , IP_CO
!             INDEX NUMBER OF CARBON MONOXIDE
     &   , IP_CH4
!             INDEX NUMBER OF METHANE
     &   , IP_O2
!             INDEX NUMBER OF OXYGEN
     &   , IP_NO
!             INDEX NUMBER OF NITROGEN MONOXIDE
     &   , IP_SO2
!             INDEX NUMBER OF SULPHUR DIOXIDE
     &   , IP_NO2
!             INDEX NUMBER OF NITROGEN DIOXIDE
     &   , IP_NH3
!             INDEX NUMBER OF AMMONIA
     &   , IP_HNO3
!             INDEX NUMBER OF NITRIC ACID
     &   , IP_N2
!             INDEX NUMBER OF NITROGEN
     &   , IP_CFC11
!             INDEX NUMBER OF CFC11
     &   , IP_CFC12
!             INDEX NUMBER OF CFC12
     &   , IP_CFC113
!             INDEX NUMBER OF CFC113
     &   , IP_HCFC22
!             INDEX NUMBER OF HCFC22
     &   , IP_HFC125
!             INDEX NUMBER OF HFC125
     &   , IP_HFC134A
!             INDEX NUMBER OF HCF134A
!
      PARAMETER(
     &     IP_H2O=1
     &   , IP_CO2=2
     &   , IP_O3=3
     &   , IP_N2O=4
     &   , IP_CO=5
     &   , IP_CH4=6
     &   , IP_O2=7
     &   , IP_NO=8
     &   , IP_SO2=9
     &   , IP_NO2=10
     &   , IP_NH3=11
     &   , IP_HNO3=12
     &   , IP_N2=13
     &   , IP_CFC11=14
     &   , IP_CFC12=15
     &   , IP_CFC113=16
     &   , IP_HCFC22=17
     &   , IP_HFC125=18
     &   , IP_HFC134A=19
     &   )
!
!     ------------------------------------------------------------------
!
      CHARACTER*80
     &     LW_SPECTRAL_FILE
!             NAME OF FILE CONTAINING THE SPECTRAL DATA
      INTEGER
     &     IERR_GET_FILE
!             ERROR FLAG RETURNED BY GET_FILE (NOT NECESSARILY
!             CONSISTENT WITH THE FLAGS IN ERROR3A).
     &   , IOS
!             STATUS OF I/O
!
      LOGICAL
     &     L_RETAIN_ABSORB(NPD_SPECIES)
!             FLAG SET TO .TRUE. IF THE ABSORBER IS TO BE RETAINED
     &   , L_GAS_INCLUDED(NPD_GASES)
!             LOGICAL TO TEST FOR ACTUAL GASES INCLUDED
      INTEGER
     &     N_ABSORB_RETAIN
!             NUMBER OF ABSORBERS TO RETAIN
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , N_AEROSOL_RETAIN
!             NUMBER OF AEROSOLS IN THE SPECTRAL FILE TO BE RETAINED
!             FOR THE RADIATIVE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)
!             INDEXING NUMBERS OF THE RETAINED AEROSOLS
     &   , COMPRESSED_INDEX(NPD_SPECIES)
!             MAPPING FROM OLD TO NEW INDICES OF ABSORBERS
     &   , N_AEROSOL_FOUND
!             NUMBER OF AEROSOLS FOR THE CURRENT GROUP OF PROCESSES
!             FOUND IN THE SPECTRAL FILE
!
!
!     DECLARE THE ELEMENTS OF THE INITIAL SPECTRUM FOR DYNAMIC
!     ALLOCATION AND SET UP AN APPROPRIATE NAMELIST.
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR SPECTRAL FILE.
!     ------------------------------------------------------------------
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT(0: NPD_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT(NPD_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE(NPD_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB(NPD_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB(NPD_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)
!             NUMBER OF ESFT TERMS IN BAND FOR EACH GAS
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND
     &        , NPD_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE(NPD_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM(NPD_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             LIST OF CONTINUA CONTINUUA IN EACH BAND
     &   , INDEX_WATER
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM(NPD_SCALE_VARIABLE, NPD_BAND, NPD_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER CORRECTION TERMS
!
!
!
!    ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE TO DECLARE ELEMENTS OF A SPECTRAL FILE AS A NAMELIST.
!
      NAMELIST/R2LWSP/
!                       Blocks Present
     &     L_PRESENT
!                       Block 0
     &   , N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL
!                       Block 1
     &   , WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG
!                       Block 2
     &   , SOLAR_FLUX_BAND
!                       Block 3
     &   , RAYLEIGH_COEFFICIENT
!                       Block 4
     &   , N_BAND_ABSORB, INDEX_ABSORB
!                       Block 5
     &   , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC
     &   , P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR
!                       Block 6
     &   , N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT
!                       Block 7
     &   , I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE
     &   , SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND
!                       Block 8
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
!                       Block 9
     &   , I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM
     &   , K_CONTINUUM, SCALE_CONTINUUM
!                       Block 10
     &   , I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
!                       Block 11
     &   , I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , HUMIDITIES
!                       Block 12
     &   , I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
!                       Block 13
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION
!                       Block 14
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
!
!     ------------------------------------------------------------------
!
!
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED LW-SPECTRAL FILE.
!     NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT.
!     NOTE: SINCE THE ARRAYS HERE WILL BE PASSED IN A COMMON BLOCK
!     THEIR SIZES MUST BE FIXED, EVEN THOUGH VARIABLE SIZES ARE USED
!     LOWER IN THE CODE. THEY ARE ACCORDINGLY DEFINED AS 1-DIMENSIONAL
!     ARRAYS WITH FIXED MAXIMUM SIZES AT THIS LEVEL.
!
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
     &     L_PRESENT_LW(0: NPD_TYPE)
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
     &     WAVE_LENGTH_SHORT_LW(NPD_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_LW(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE_LW(NPD_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_LW(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND_LW(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT_LW(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB_LW
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_LW(NPD_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_LW(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_LW(NPD_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_LW(NPD_BAND, NPD_SPECIES)
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_LW(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_LW(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_LW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM
     &        , NPD_BAND, NPD_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_LW(NPD_SPECIES, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_LW(NPD_SPECIES, NPD_BAND)
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
     &     THERMAL_COEFFICIENT_LW(0: NPD_THERMAL_COEFF-1
     &        , NPD_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_LW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE_LW(NPD_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_LW(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE_LW(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO_LW(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_LW(0: NPD_ALBEDO_PARM
     &      , NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_LW(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM_LW(NPD_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_LW
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_LW(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_LW(NPD_SCALE_VARIABLE
     &      , NPD_BAND, NPD_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION_LW(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE_LW(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_LW(NPD_DROP_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_LW(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL_LW
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_LW(NPD_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_LW(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_LW(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES_LW(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION_LW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_LW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_LW(NPD_HUMIDITIES
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION_LW(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE_LW(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_LW(NPD_ICE_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_LW(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT_LW(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION_LW(NPD_SPECIES)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!    ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE DECLARING COMMON BLOCK CONTAINING THE REDUCED LW SPECTRAL
!     FILE.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
!
      COMMON/R2LWSPCM/
!
!     DIMENSIONS OF ARRAYS:
     &     NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW
     &   , NPD_SCALE_VARIABLE_LW
     &   , NPD_THERMAL_COEFF_LW
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW
     &   , NPD_CONTINUUM_LW
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW
!
!     GENERAL ARRAYS:
     &   , L_PRESENT_LW
!
!     PROPERTIES OF BANDS:
     &   , N_BAND_LW, WAVE_LENGTH_SHORT_LW, WAVE_LENGTH_LONG_LW
!
!     EXCLUSIONS FROM BANDS:
     &   , N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW
!
!     SOLAR FIELDS:
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW
!
!     GASEOUS ABSORPTION:
     &   , N_ABSORB_LW, N_BAND_ABSORB_LW, INDEX_ABSORB_LW
     &   , TYPE_ABSORB_LW
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW
     &   , K_ESFT_LW, W_ESFT_LW
     &   , SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW
!
!     THERMAL SOURCE FUNCTION:
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW
!
!     SURFACE PROPERTIES:
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW
     &   , EMISSIVITY_GROUND_LW
!
!     CONTINUA:
     &   , N_BAND_CONTINUUM_LW, INDEX_CONTINUUM_LW, INDEX_WATER_LW
     &   , I_SCALE_FNC_CONT_LW, K_CONTINUUM_LW
     &   , SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW, T_REF_CONTINUUM_LW
!
!     WATER DROPLETS:
     &   , I_DROP_PARAMETRIZATION_LW, L_DROP_TYPE_LW
     &   , DROP_PARAMETER_LIST_LW
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW
!
!     AEROSOLS:
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW, I_AEROSOL_PARAMETRIZATION_LW
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW
     &   , AEROSOL_ASYMMETRY_LW
!
!     ICE CRYSTALS:
     &   , I_ICE_PARAMETRIZATION_LW, L_ICE_TYPE_LW
     &   , ICE_PARAMETER_LIST_LW
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW
!
!     DOPPLER BROADENING:
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW
!
!     ------------------------------------------------------------------
!
!
!
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
!
      CHARACTER
     &     CH_IOS*5
!             CHARACTER STRING FOR IOSTAT ERROR
!
!     SUBROUTINES CALLED
      EXTERNAL
     &     R2_COMPRESS_SPECTRUM
!
!
!     EACH BLOCK IS INITIALIZED AS MISSING:
      DATA L_PRESENT/.FALSE., NPD_TYPE*.FALSE./
!
!     INITIALIZE THE RANGE OF VALIDITY OF THE PARAMETRIZATIONS OF
!     DROPLETS AND ICE CRYSTALS. OLD SPECTRAL FILES WILL NOT CONTAIN
!     SUCH DATA, SO THE LIMITS FOR DROPLETS ARE INITIALIZED TO THOSE
!     FORMERLY SET IN THE MICROPHYSICAL SCHEME (MRF/UMIST
!     PARAMETRIZATION) TO ENSURE THAT THE RESULTS ARE BIT-REPRODUCIBLE.
!     VALUES FOR ICE COVER THE RANGE OF EFFECTIVE RADII USED IN
!     GENERATING THE DATA FOR THE ORIGINAL PARAMETRIZATION OF ICE
!     CRYSTALS.
!     AT SOME FUTURE RELEASE IT MAY BE DESIRABLE TO REMOVE DEFAULT
!     SETTINGS.
      DATA DROP_PARM_MIN_DIM/NPD_DROP_TYPE*3.5E-07/
      DATA DROP_PARM_MAX_DIM/NPD_DROP_TYPE*3.7E-05/
      DATA ICE_PARM_MIN_DIM/NPD_ICE_TYPE*3.75E-07/
      DATA ICE_PARM_MAX_DIM/NPD_ICE_TYPE*8.0E-05/
!
!
!
!     READ THE LONGWAVE SPECTRUM AS A NAMELIST.
      CALL GET_FILE(80, LW_SPECTRAL_FILE, 80, IERR_GET_FILE)
      IF (IERR_GET_FILE.NE.0) THEN
!        CONVERT THE ERROR FLAG FROM GET_FILE TO A FLAG RECOGNISED
!        BY THE RADIATION CODE.
         IERR=I_ERR_IO
         CMESSAGE='Error reading name of longwave spectral file.'
         RETURN
      ENDIF
      OPEN(UNIT=80, FILE=LW_SPECTRAL_FILE, IOSTAT=IOS,
     & DELIM='APOSTROPHE')
      IF (IOS.NE.0) THEN
         IERR=I_ERR_IO
      WRITE(CH_IOS, '(I5)') IOS
         CMESSAGE='Error opening longwave spectral file.'
     &      //' IOSTAT='//CH_IOS
         RETURN
      ENDIF
      READ(80, R2LWSP)
      CLOSE(80)
!
!     TEST FOR MINIMAL REQUISITE INFORMATION.
      IF ( .NOT.(L_PRESENT(0).AND.
     &           L_PRESENT(6) ) ) THEN
         CMESSAGE='Longwave spectrum is deficient.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET REDUCED DIMENSIONS, EITHER FROM THE SIZES OF THE FIXED ARRAYS
!     OR FROM THE ARRAYS READ IN.
!
      NPD_TYPE_LW=NPD_TYPE
      NPD_BAND_LW=MAX(N_BAND, 1)
      NPD_SPECIES_LW=MAX(N_ABSORB, 1)
      NPD_ALBEDO_PARM_LW=NPD_ALBEDO_PARM
      NPD_SCALE_FNC_LW=NPD_SCALE_FNC
      NPD_SCALE_VARIABLE_LW=NPD_SCALE_VARIABLE
      NPD_SURFACE_LW=NPD_SURFACE
      NPD_CONTINUUM_LW=NPD_CONTINUUM
      NPD_THERMAL_COEFF_LW=N_DEG_FIT+1
      NPD_CLOUD_PARAMETER_LW=NPD_CLOUD_PARAMETER
!
!
!     SEARCH THE SPECTRUM TO FIND MAXIMUM DIMENSIONS.
!
      NPD_EXCLUDE_LW=1
      IF (L_PRESENT(14)) THEN
         DO I=1, N_BAND
            NPD_EXCLUDE_LW=MAX(NPD_EXCLUDE_LW, N_BAND_EXCLUDE(I))
         ENDDO
      ENDIF
!
!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      DO I=1, NPD_GASES
         L_GAS_INCLUDED(I)=.FALSE.
      ENDDO
      N_ABSORB_RETAIN=0
!
      DO I=1, N_ABSORB
!
         L_RETAIN_ABSORB(I)=.FALSE.
         COMPRESSED_INDEX(I)=0
!
         IF ( (TYPE_ABSORB(I).EQ.IP_H2O).OR.
     &        (TYPE_ABSORB(I).EQ.IP_CO2).OR.
     &        (TYPE_ABSORB(I).EQ.IP_O3).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_CH4).AND.L_CH4 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_N2O).AND.L_N2O ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_CFC11).AND.L_CFC11 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_CFC12).AND.L_CFC12 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_CFC113).AND.L_CFC113 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_HCFC22).AND.L_HCFC22 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_HFC125).AND.L_HFC125 ).OR.
     &        ( (TYPE_ABSORB(I).EQ.IP_HFC134A).AND.L_HFC134A ) ) THEN
            N_ABSORB_RETAIN=N_ABSORB_RETAIN+1
            INDEX_ABSORB_RETAIN(N_ABSORB_RETAIN)=I
            COMPRESSED_INDEX(I)=N_ABSORB_RETAIN
            L_RETAIN_ABSORB(I)=.TRUE.
            L_GAS_INCLUDED(TYPE_ABSORB(I))=.TRUE.
         ENDIF
!
      ENDDO
!
!
!     Print warning messages if those gases normally expected
!     are not present.
      IF (.NOT.L_GAS_INCLUDED(IP_H2O)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Water vapour is not included in the '
     &      , 'longwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_CO2)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Carbon dioxide is not included in the '
     &      , 'longwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_O3)) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** WARNING: Ozone is not included in the '
     &      , 'longwave spectral file.'
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CH4)).AND.L_CH4) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: Methane is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_N2O)).AND.L_N2O) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: Nitrous oxide is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC11)).AND.L_CFC11) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: CFC11 is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC12)).AND.L_CFC12) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: CFC12 is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC113)).AND.L_CFC113) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: CFC113 is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HCFC22)).AND.L_HCFC22) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: HCFC22 is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HFC125)).AND.L_HFC125) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: HFC125 is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HFC134A)).AND.L_HFC134A) THEN
         WRITE(IU_ERR, '(/A, /A)')
     &      '*** ERROR: HFC134A is not included in the longwave '
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_SPECIES_LW=MAX(N_ABSORB_RETAIN, 1)
!
      NPD_ESFT_TERM_LW=1
      IF (L_PRESENT(5)) THEN
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I)))
     &            NPD_ESFT_TERM_LW=MAX(NPD_ESFT_TERM_LW
     &            , I_BAND_ESFT(I, INDEX_ABSORB(J, I)))
            ENDDO
         ENDDO
      ENDIF
!
      NPD_DROP_TYPE_LW=1
      IF (L_PRESENT(10)) THEN
         DO I=1, NPD_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               NPD_DROP_TYPE_LW=MAX(NPD_DROP_TYPE_LW, I)
            ENDIF
         ENDDO
      ENDIF
!
      NPD_ICE_TYPE_LW=1
      IF (L_PRESENT(12)) THEN
         DO I=1, NPD_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               NPD_ICE_TYPE_LW=MAX(NPD_ICE_TYPE_LW, I)
            ENDIF
         ENDDO
      ENDIF
!
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      NPD_HUMIDITIES_LW=1
      N_AEROSOL_RETAIN=0
!
!     Check the spectral file for climatological aerosols
      IF (L_CLIMAT_AEROSOL) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I).EQ.IP_WATER_SOLUBLE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_DUST_LIKE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_OCEANIC).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_SOOT).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_SULPHURIC) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.5) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '
     &            //'climatological aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!     Check the spectral file for soot aerosols.
      IF (L_USE_SOOT_DIRECT) THEN
         IF (L_PRESENT(11)) THEN
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
               IF ((TYPE_AEROSOL(I).EQ.IP_FRESH_SOOT) .OR.
     &             (TYPE_AEROSOL(I).EQ.IP_AGED_SOOT)) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '
     &            //'soot aerosol data.'
               RETURN
!
            ENDIF
!
         ELSE
!
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no soot data.'
            RETURN
!
         ENDIF
!
      ENDIF
!

!
!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).
!
      IF (L_USE_SULPC_DIRECT) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I).EQ.IP_ACCUM_SULPHATE).OR.
     &              (TYPE_AEROSOL(I).EQ.IP_AITKEN_SULPHATE) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND.NE.2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '
     &            //'sulphate aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_AEROSOL_SPECIES_LW=MAX(N_AEROSOL_RETAIN, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      IF (L_PRESENT(11)) THEN
         DO I=1, N_AEROSOL_RETAIN
            IF (I_AEROSOL_PARAMETRIZATION(INDEX_AEROSOL_RETAIN(I)).EQ.
     &         IP_AEROSOL_PARAM_MOIST) THEN
               NPD_HUMIDITIES_LW=MAX(NPD_HUMIDITIES_LW
     &            , NHUMIDITY(INDEX_AEROSOL_RETAIN(I)))
            ENDIF
         ENDDO
      ENDIF
!
!
!
!
!     TRANSFER THE LARGE NAMELIST TO THE REDUCED SPECTRUM.
!
!
      CALL R2_COMPRESS_SPECTRUM(
!                       Spectral Array in Namelist
     &     L_PRESENT
     &   , N_BAND, WAVE_LENGTH_SHORT , WAVE_LENGTH_LONG
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
     &   , SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT
     &   , N_ABSORB, N_BAND_ABSORB, INDEX_ABSORB, TYPE_ABSORB
     &   , L_RETAIN_ABSORB, N_ABSORB_RETAIN, INDEX_ABSORB_RETAIN
     &   , COMPRESSED_INDEX, I_BAND_ESFT, K_ESFT, W_ESFT, I_SCALE_ESFT
     &   , I_SCALE_FNC, SCALE_VECTOR, P_REFERENCE, T_REFERENCE
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK
     &   , I_SPEC_SURFACE, L_SURFACE, SURFACE_ALBEDO
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
     &   , N_AEROSOL, TYPE_AEROSOL
     &   , N_AEROSOL_RETAIN, INDEX_AEROSOL_RETAIN
     &   , L_AEROSOL_SPECIES, AEROSOL_ABSORPTION
     &   , AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , NHUMIDITY, HUMIDITIES, I_AEROSOL_PARAMETRIZATION
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION
!                       Reduced Spectral Array
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
     &   )
!
!
!
      RETURN
      END
!+ Subroutine to transfer spectrum to reduced array.
!
! Purpose:
!       Spectral data from the large dynamically allocated array
!       are transferred to the reduced array.
!
! Method:
!       Elements are copied across.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             03-09-97                Coding changes
!                                               associated with the
!                                               removal of pointers
!                                               into the spectral data.
!                                               Capability to select
!                                               aerosols added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COMPRESS_SPECTRUM(
!                       Original Spectrum
     &     L_PRESENT
     &   , N_BAND, WAVE_LENGTH_SHORT , WAVE_LENGTH_LONG
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE
     &   , SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT
     &   , N_ABSORB, N_BAND_ABSORB, INDEX_ABSORB, TYPE_ABSORB
     &   , L_RETAIN_ABSORB, N_ABSORB_RETAIN, INDEX_ABSORB_RETAIN
     &   , COMPRESSED_INDEX, I_BAND_ESFT, K_ESFT, W_ESFT, I_SCALE_ESFT
     &   , I_SCALE_FNC, SCALE_VECTOR, P_REFERENCE, T_REFERENCE
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK
     &   , I_SPEC_SURFACE, L_SURFACE, SURFACE_ALBEDO
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM
     &   , N_AEROSOL, TYPE_AEROSOL
     &   , N_AEROSOL_RETAIN, INDEX_AEROSOL_RETAIN
     &   , L_AEROSOL_SPECIES, AEROSOL_ABSORPTION
     &   , AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , NHUMIDITY, HUMIDITIES, I_AEROSOL_PARAMETRIZATION
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION
!                       Reduced Spectral Array
     &   , NPDR_TYPE, NPDR_BAND, NPDR_EXCLUDE
     &   , NPDR_SPECIES, NPDR_ESFT_TERM, NPDR_SCALE_FNC
     &   , NPDR_SCALE_VARIABLE, NPDR_THERMAL_COEFF
     &   , NPDR_SURFACE, NPDR_ALBEDO_PARM
     &   , NPDR_CONTINUUM, NPDR_DROP_TYPE, NPDR_ICE_TYPE
     &   , NPDR_CLOUD_PARAMETER, NPDR_AEROSOL_SPECIES
     &   , NPDR_HUMIDITIES
     &   , L_PRESENT_RD
     &   , N_BAND_RD, WAVE_LENGTH_SHORT_RD , WAVE_LENGTH_LONG_RD
     &   , N_BAND_EXCLUDE_RD, INDEX_EXCLUDE_RD
     &   , SOLAR_FLUX_BAND_RD, RAYLEIGH_COEFFICIENT_RD
     &   , N_ABSORB_RD, N_BAND_ABSORB_RD, INDEX_ABSORB_RD
     &   , TYPE_ABSORB_RD
     &   , I_BAND_ESFT_RD, I_SCALE_ESFT_RD, I_SCALE_FNC_RD
     &   , K_ESFT_RD, W_ESFT_RD, SCALE_VECTOR_RD
     &   , P_REFERENCE_RD, T_REFERENCE_RD
     &   , N_DEG_FIT_RD, THERMAL_COEFFICIENT_RD, T_REF_PLANCK_RD
     &   , I_SPEC_SURFACE_RD, N_DIR_ALBEDO_FIT_RD
     &   , L_SURFACE_RD, SURFACE_ALBEDO_RD, DIRECT_ALBEDO_PARM_RD
     &   , EMISSIVITY_GROUND_RD
     &   , N_BAND_CONTINUUM_RD, INDEX_CONTINUUM_RD, INDEX_WATER_RD
     &   , I_SCALE_FNC_CONT_RD, K_CONTINUUM_RD, SCALE_CONTINUUM_RD
     &   , P_REF_CONTINUUM_RD, T_REF_CONTINUUM_RD
     &   , I_DROP_PARAMETRIZATION_RD, L_DROP_TYPE_RD
     &   , DROP_PARAMETER_LIST_RD
     &   , DROP_PARM_MIN_DIM_RD, DROP_PARM_MAX_DIM_RD
     &   , N_AEROSOL_RD, TYPE_AEROSOL_RD, I_AEROSOL_PARAMETRIZATION_RD
     &   , NHUMIDITY_RD, HUMIDITIES_RD
     &   , L_AEROSOL_SPECIES_RD, AEROSOL_ABSORPTION_RD
     &   , AEROSOL_SCATTERING_RD, AEROSOL_ASYMMETRY_RD
     &   , I_ICE_PARAMETRIZATION_RD, L_ICE_TYPE_RD
     &   , ICE_PARAMETER_LIST_RD
     &   , ICE_PARM_MIN_DIM_RD, ICE_PARM_MAX_DIM_RD
     &   , L_DOPPLER_PRESENT_RD, DOPPLER_CORRECTION_RD
     &   )
!
!
      IMPLICIT NONE
!
!
!
!     ------------------------------------------------------------------
!     DECLARATION OF INITIAL SPECTRUM.
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
      INTEGER
     &     NPD_TYPE
!             NUMBER OF TYPES OF DATA
     &   , NPD_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , NPD_EXCLUDE
!             NUMER OF EXCLUDED BANDS
     &   , NPD_SPECIES
!             NUMBER OF GASEOUS SPECIES
     &   , NPD_ESFT_TERM
!             NUMBER OF ESFT TERMS
     &   , NPD_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             NUMBER OF SCALING VARIABLES
     &   , NPD_SURFACE
!             NUMBER OF SURFACE TYPES
     &   , NPD_ALBEDO_PARM
!             NUMBER OF ALBEDO PARAMETERS
     &   , NPD_CONTINUUM
!             NUMBER OF CONTINUA
     &   , NPD_DROP_TYPE
!             NUMBER OF DROP TYPES
     &   , NPD_ICE_TYPE
!             NUMBER OF ICE CRYSTAL TYPES
     &   , NPD_AEROSOL_SPECIES
!             NUMBER OF AEROSOL SPECIES
     &   , NPD_CLOUD_PARAMETER
!             MAX NUMBER OF CLOUD PARAMETERS
     &   , NPD_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
!
      PARAMETER(
     &     NPD_TYPE=15
     &   , NPD_BAND=20
     &   , NPD_EXCLUDE=2
     &   , NPD_SPECIES=11
     &   , NPD_ESFT_TERM=16
     &   , NPD_SCALE_FNC=3
     &   , NPD_SCALE_VARIABLE=4
     &   , NPD_THERMAL_COEFF=9
     &   , NPD_SURFACE=1
     &   , NPD_ALBEDO_PARM=4
     &   , NPD_CONTINUUM=2
     &   , NPD_DROP_TYPE=6
     &   , NPD_ICE_TYPE=7
     &   , NPD_CLOUD_PARAMETER=30
     &   , NPD_AEROSOL_SPECIES=9
     &   , NPD_HUMIDITIES=21
     &   )
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR SPECTRAL FILE.
!     ------------------------------------------------------------------
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT(0: NPD_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT(NPD_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER
     &     N_BAND_EXCLUDE(NPD_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB(NPD_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB(NPD_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)
!             NUMBER OF ESFT TERMS IN BAND FOR EACH GAS
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND
     &        , NPD_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE(NPD_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM(NPD_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             LIST OF CONTINUA CONTINUUA IN EACH BAND
     &   , INDEX_WATER
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM(NPD_SCALE_VARIABLE, NPD_BAND, NPD_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES
     &        , NPD_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_BAND, NPD_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER CORRECTION TERMS
!
!
!
!    ------------------------------------------------------------------
!
!     AUXILIARY VARIABLES USED TO SELECT PARTS OF THE INITIAL SPECTRUM
      LOGICAL   !, INTENT(IN)
     &     L_RETAIN_ABSORB(NPD_SPECIES)
!             FLAGS FOR THE RETENTION OF GASES IN THE SPECTRAL FILE
      INTEGER   !, INTENT(IN)
     &     N_ABSORB_RETAIN
!             NUMBER OF ABSORBERS TO BE RETAINED
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , COMPRESSED_INDEX(NPD_SPECIES)
!             MAPPING FROM OLD TO NEW INDICES OF ABSORBERS
     &   , N_AEROSOL_RETAIN
!             NUMBER OF AEROSOLS IN THE INITIAL SPECTRUM TO BE USED
!             IN THE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)
!             INDICES OF THE RETAINED AEROSOLS
!
!
!     ------------------------------------------------------------------
!     DECLARATION OF REDUCED SPECTRUM.
!     ------------------------------------------------------------------
!
!     DIMENSIONS OF REDUCED ARRAY:
!
      INTEGER
     &     NPDR_BAND
!             NUMBER OF SPECTRAL BANDS
     &   , NPDR_EXCLUDE
!             NUMER OF EXCLUDED BANDS
     &   , NPDR_ESFT_TERM
!             NUMBER OF ESFT TERMS
     &   , NPDR_TYPE
!             NUMBER OF DATA TYPES
     &   , NPDR_SPECIES
!             NUMBER OF GASEOUS SPECIES
     &   , NPDR_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPDR_SCALE_VARIABLE
!             NUMBER OF SCALING VARIABLES
     &   , NPDR_SURFACE
!             NUMBER OF SURFACE TYPES
     &   , NPDR_ALBEDO_PARM
!             NUMBER OF ALBEDO PARAMETERS
     &   , NPDR_CONTINUUM
!             NUMBER OF CONTINUA
     &   , NPDR_DROP_TYPE
!             NUMBER OF DROP TYPES
     &   , NPDR_ICE_TYPE
!             NUMBER OF ICE CRYSTAL TYPES
     &   , NPDR_AEROSOL_SPECIES
!             NUMBER OF AEROSOL SPECIES
     &   , NPDR_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
     &   , NPDR_CLOUD_PARAMETER
!             MAX NUMBER OF CLOUD PARAMETERS
     &   , NPDR_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
!
!
!
!     GENERAL FIELDS:
!
      LOGICAL
     &     L_PRESENT_RD(0: NPDR_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER
     &     N_BAND_RD
!             NUMBER OF SPECTRAL BANDS
!
      REAL
     &     WAVE_LENGTH_SHORT_RD(NPDR_BAND)
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_RD(NPDR_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
      INTEGER
     &     N_BAND_EXCLUDE_RD(NPDR_BAND)
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_RD(NPDR_EXCLUDE, NPDR_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL
     &     SOLAR_FLUX_BAND_RD(NPDR_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL
     &     RAYLEIGH_COEFFICIENT_RD(NPDR_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER
     &     N_ABSORB_RD
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_RD(NPDR_BAND)
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_RD(NPDR_SPECIES, NPDR_BAND)
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_RD(NPDR_SPECIES)
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_RD(NPDR_BAND, NPDR_SPECIES)
!             NUMBER OF ESFT TERMS IN BAND FOR EACH GAS
     &   , I_SCALE_ESFT_RD(NPDR_BAND, NPDR_SPECIES)
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_RD(NPDR_BAND, NPDR_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL
     &     K_ESFT_RD(NPDR_ESFT_TERM, NPDR_BAND, NPDR_SPECIES)
!             ESFT EXPONENTS
     &   , W_ESFT_RD(NPDR_ESFT_TERM, NPDR_BAND, NPDR_SPECIES)
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_RD(NPDR_SCALE_VARIABLE, NPDR_ESFT_TERM
     &        , NPDR_BAND, NPDR_SPECIES)
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_RD(NPDR_SPECIES, NPDR_BAND)
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_RD(NPDR_SPECIES, NPDR_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER
     &     N_DEG_FIT_RD
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL
     &     THERMAL_COEFFICIENT_RD(0: NPDR_THERMAL_COEFF-1, NPDR_BAND)
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_RD
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER
     &     I_SPEC_SURFACE_RD(NPDR_SURFACE)
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_RD(NPDR_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL
     &     L_SURFACE_RD(NPDR_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL
     &     SURFACE_ALBEDO_RD(NPDR_BAND, NPDR_SURFACE)
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_RD(0: NPDR_ALBEDO_PARM
     &        , NPD_BAND, NPD_SURFACE)
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_RD(NPDR_BAND, NPDR_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER
     &     N_BAND_CONTINUUM_RD(NPDR_BAND)
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_RD(NPDR_BAND, NPDR_CONTINUUM)
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_RD
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_RD(NPDR_BAND, NPDR_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL
     &     K_CONTINUUM_RD(NPDR_BAND, NPDR_CONTINUUM)
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_RD(NPDR_SCALE_VARIABLE
     &        , NPDR_BAND, NPDR_CONTINUUM)
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_RD(NPDR_CONTINUUM, NPDR_BAND)
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_RD(NPDR_CONTINUUM, NPDR_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER
     &     I_DROP_PARAMETRIZATION_RD(NPDR_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL
     &     L_DROP_TYPE_RD(NPDR_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL
     &     DROP_PARAMETER_LIST_RD(NPDR_CLOUD_PARAMETER
     &        , NPDR_BAND, NPDR_DROP_TYPE)
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_RD(NPDR_DROP_TYPE)
!             MINIMUM SIZE OF DROPLET PERMITTED IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_RD(NPDR_DROP_TYPE)
!             MAXIMUM SIZE OF DROPLET PERMITTED IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER
     &     N_AEROSOL_RD
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_RD(NPDR_AEROSOL_SPECIES)
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_RD(NPDR_AEROSOL_SPECIES)
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_RD(NPDR_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL
     &     L_AEROSOL_SPECIES_RD(NPDR_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL
     &     AEROSOL_ABSORPTION_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES
     &        , NPDR_BAND)
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES
     &        , NPDR_BAND)
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES
     &        , NPDR_BAND)
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER
     &     I_ICE_PARAMETRIZATION_RD(NPDR_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL
     &     L_ICE_TYPE_RD(NPDR_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL
     &     ICE_PARAMETER_LIST_RD(NPDR_CLOUD_PARAMETER
     &        , NPDR_BAND, NPDR_ICE_TYPE)
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_RD(NPDR_ICE_TYPE)
!             MINIMUM SIZE OF ICE CRYSTAL PERMITTED
!             IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_RD(NPDR_ICE_TYPE)
!             MAXIMUM SIZE OF ICE CRYSTAL PERMITTED
!             IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL
     &     L_DOPPLER_PRESENT_RD(NPDR_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL
     &     DOPPLER_CORRECTION_RD(NPDR_SPECIES)
!             DOPPLER CORRECTION TERMS
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , N_PARAMETER
!             NUMBER OF PARAMETERS IN SCHEME.
     &   , I_SPECIES
!             SPECIES OF GAS
     &   , I_CONTINUUM
!             TYPE OF CONTINUUM
     &   , I_INITIAL
!             INDEXING NUMBER IN INITIAL SPECTRAL FILE
!
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET TYPES OF SCALING FOR ABSORBER AMOUNTS
!
      INTEGER
     &     IP_SCALE_FNC_NULL
!             NULL SCALING FUNCTION
     &   , IP_SCALE_POWER_LAW
!             POWER LAW SCALING FUNCTION
     &   , IP_SCALE_POWER_QUAD
!             POWER LAW FOR P; QUADRATIC FOR T
     &   , IP_SCALE_DOPPLER_QUAD
!             POWER LAW FOR P; QUADRATIC FOR T
!              WITH IMPLICIT DOPPLER CORRECTION
     &   , N_SCALE_VARIABLE(0: NPD_SCALE_FNC)
!             NUMBER OF SCALING VARIABLES
!
      PARAMETER(
     &     IP_SCALE_FNC_NULL=0
     &   , IP_SCALE_POWER_LAW=1
     &   , IP_SCALE_POWER_QUAD=2
     &   , IP_SCALE_DOPPLER_QUAD=3
     &   )
!
!     -----------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET NUMBERS FOR WATER CLOUD SCHEMES.
!
      INTEGER
     &     NPD_CLOUD_FIT
!             NUMBER OF CLOUD FITTING SCHEMES
     &   , IP_SLINGO_SCHRECKER
!             PARAMETRIZATION OF SLINGO-SCHRECKER
     &   , IP_ACKERMAN_STEPHENS
!             PARAMETRIZATION OF ACKERMAN & STEPHENS
     &   , IP_DROP_UNPARAMETRIZED
!             UNPARAMETRIZED DROPLET DATA
     &   , IP_DROP_PADE_2
!             PADE APPROXIMATION OF THE SECOND ORDER
!             (THIRD ORDER FOR THE EXTINCTION)
!
!
      PARAMETER(
     &     NPD_CLOUD_FIT=3
     &   , IP_SLINGO_SCHRECKER=1
     &   , IP_ACKERMAN_STEPHENS=2
     &   , IP_DROP_UNPARAMETRIZED=3
     &   , IP_DROP_PADE_2=5
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET NUMBERS FOR ICE CLOUD SCHEMES.
!
      INTEGER
     &     NPD_ICE_CLOUD_FIT
!             NUMBER OF CLOUD FITTING SCHEMES
     &   , IP_SLINGO_SCHRECKER_ICE
!             PARAMETRIZATION OF SLINGO AND SCHRECKER.
     &   , IP_ICE_UNPARAMETRIZED
!             UNPARAMETRIZED ICE CRYSTAL DATA
     &   , IP_SUN_SHINE_VN2_VIS
!             SUN AND SHINE'S PARAMETRIZATION IN THE VISIBLE (VERSION 2)
     &   , IP_SUN_SHINE_VN2_IR
!             SUN AND SHINE'S PARAMETRIZATION IN THE IR (VERSION 2)
     &   , IP_ICE_ADT
!             SCHEME BASED ON ANOMALOUS DIFFRACTION THEORY
!             FOR ICE CRYSTALS

!
      PARAMETER(
     &     NPD_ICE_CLOUD_FIT=6
     &   , IP_SLINGO_SCHRECKER_ICE=1
     &   , IP_ICE_UNPARAMETRIZED=3
     &   , IP_SUN_SHINE_VN2_VIS=4
     &   , IP_SUN_SHINE_VN2_IR=5
     &   , IP_ICE_ADT=6
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
!     MODULE TO SET DATA FOR SCALING FUNCTIONS.
!
      DATA N_SCALE_VARIABLE(IP_SCALE_POWER_LAW)/2/
     &  ,  N_SCALE_VARIABLE(IP_SCALE_FNC_NULL)/0/
     &  ,  N_SCALE_VARIABLE(IP_SCALE_POWER_QUAD)/3/
     &  ,  N_SCALE_VARIABLE(IP_SCALE_DOPPLER_QUAD)/4/
!
!     ------------------------------------------------------------------
!
!
!
!
!
!     INITAILIZE ALL BLOCKS OF THE COMPRESSED SPECTRUM TO .FALSE.
      DO I=1, NPDR_TYPE
         L_PRESENT_RD(I)=.FALSE.
      ENDDO
!
!
!     PROCEED THROUGH EACH BLOCK OF THE SPECTRAL FILE TRANSFERRING
!     THE DATA FROM THE INPUT ARRAY TO THE REDUCED ARRAY.
!
!
!     BLOCK 0:
!
      IF (L_PRESENT(0)) THEN
         L_PRESENT_RD(0)=.TRUE.
         N_BAND_RD=N_BAND
         N_ABSORB_RD=N_ABSORB_RETAIN
         N_AEROSOL_RD=N_AEROSOL_RETAIN
         DO I=1, N_ABSORB_RETAIN
            TYPE_ABSORB_RD(I)=TYPE_ABSORB(INDEX_ABSORB_RETAIN(I))
         ENDDO
         DO I=1, N_AEROSOL_RETAIN
            TYPE_AEROSOL_RD(I)=TYPE_AEROSOL(INDEX_AEROSOL_RETAIN(I))
         ENDDO
      ENDIF
!
!     BLOCK 1:
      IF (L_PRESENT(1)) THEN
         L_PRESENT_RD(1)=.TRUE.
         DO I=1, N_BAND
            WAVE_LENGTH_SHORT_RD(I)=WAVE_LENGTH_SHORT(I)
            WAVE_LENGTH_LONG_RD(I)=WAVE_LENGTH_LONG(I)
         ENDDO
      ENDIF
!
!     BLOCK 2:
      IF (L_PRESENT(2)) THEN
         L_PRESENT_RD(2)=.TRUE.
         DO I=1, N_BAND
            SOLAR_FLUX_BAND_RD(I)=SOLAR_FLUX_BAND(I)
         ENDDO
      ENDIF
!
!     BLOCK 3:
      IF (L_PRESENT(3)) THEN
         L_PRESENT_RD(3)=.TRUE.
         DO I=1, N_BAND
            RAYLEIGH_COEFFICIENT_RD(I)=RAYLEIGH_COEFFICIENT(I)
         ENDDO
      ENDIF
!
!     BLOCK 4:
      IF (L_PRESENT(4)) THEN
         L_PRESENT_RD(4)=.TRUE.
         DO I=1, N_BAND
            N_BAND_ABSORB_RD(I)=0
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I))) THEN
                  N_BAND_ABSORB_RD(I)=N_BAND_ABSORB_RD(I)+1
                  INDEX_ABSORB_RD(N_BAND_ABSORB_RD(I), I)
     &               =COMPRESSED_INDEX(INDEX_ABSORB(J, I))
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 5:
      IF (L_PRESENT(5)) THEN
         L_PRESENT_RD(5)=.TRUE.
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB_RD(I)
               I_SPECIES=INDEX_ABSORB_RD(J, I)
               I_INITIAL=INDEX_ABSORB_RETAIN(I_SPECIES)
               I_BAND_ESFT_RD(I, I_SPECIES)=I_BAND_ESFT(I, I_INITIAL)
               I_SCALE_ESFT_RD(I, I_SPECIES)=I_SCALE_ESFT(I, I_INITIAL)
               I_SCALE_FNC_RD(I, I_SPECIES)=I_SCALE_FNC(I, I_INITIAL)
               P_REFERENCE_RD(I_SPECIES, I)=P_REFERENCE(I_INITIAL, I)
               T_REFERENCE_RD(I_SPECIES, I)=T_REFERENCE(I_INITIAL, I)
               DO K=1, I_BAND_ESFT(I, I_INITIAL)
                  K_ESFT_RD(K, I, I_SPECIES)=K_ESFT(K, I, I_INITIAL)
                  W_ESFT_RD(K, I, I_SPECIES)=W_ESFT(K, I, I_INITIAL)
                  DO L=1, N_SCALE_VARIABLE(I_SCALE_FNC(I, I_INITIAL))
                     SCALE_VECTOR_RD(L, K, I, I_SPECIES)
     &                  =SCALE_VECTOR(L, K, I, I_INITIAL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 6:
      IF (L_PRESENT(6)) THEN
         L_PRESENT_RD(6)=.TRUE.
         N_DEG_FIT_RD=N_DEG_FIT
         T_REF_PLANCK_RD=T_REF_PLANCK
         DO I=1, N_BAND
            DO J=0, N_DEG_FIT
               THERMAL_COEFFICIENT_RD(J, I)=THERMAL_COEFFICIENT(J, I)
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 7:
!
!     OMITTED SINCE SURFACE ALBEDOS ARE PROVIDED BY THE MODEL.
!
!     BLOCK 8:
      IF (L_PRESENT(8)) THEN
         L_PRESENT_RD(8)=.TRUE.
         DO I=1, N_BAND
            N_BAND_CONTINUUM_RD(I)=N_BAND_CONTINUUM(I)
            DO J=1, N_BAND_CONTINUUM(I)
               INDEX_CONTINUUM_RD(I, J)=INDEX_CONTINUUM(I, J)
            ENDDO
         ENDDO
!
         INDEX_WATER_RD=0
         DO I=1, N_ABSORB_RETAIN
            IF (INDEX_ABSORB_RETAIN(I).EQ.INDEX_WATER) THEN
               INDEX_WATER_RD=I
            ENDIF
         ENDDO
!
      ENDIF
!
!     BLOCK 9:
      IF (L_PRESENT(9)) THEN
         L_PRESENT_RD(9)=.TRUE.
         DO I=1, N_BAND
            DO J=1, N_BAND_CONTINUUM(I)
               I_CONTINUUM=INDEX_CONTINUUM(I, J)
               I_SCALE_FNC_CONT_RD(I, I_CONTINUUM)
     &            =I_SCALE_FNC_CONT(I, I_CONTINUUM)
               P_REF_CONTINUUM_RD(I_CONTINUUM, I)
     &            =P_REF_CONTINUUM(I_CONTINUUM, I)
               T_REF_CONTINUUM_RD(I_CONTINUUM, I)
     &            =T_REF_CONTINUUM(I_CONTINUUM, I)
               K_CONTINUUM_RD(I, I_CONTINUUM)
     &            =K_CONTINUUM(I, I_CONTINUUM)
               DO L=1, N_SCALE_VARIABLE(I_SCALE_FNC_CONT
     &               (I, I_CONTINUUM))
                  SCALE_CONTINUUM_RD(L, I, I_CONTINUUM)
     &               =SCALE_CONTINUUM(L, I, I_CONTINUUM)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 10:
      IF (L_PRESENT(10)) THEN
         L_PRESENT_RD(10)=.TRUE.
         DO I=1, NPDR_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               L_DROP_TYPE_RD(I)=.TRUE.
               I_DROP_PARAMETRIZATION_RD(I)=I_DROP_PARAMETRIZATION(I)
               DROP_PARM_MIN_DIM_RD(I)=DROP_PARM_MIN_DIM(I)
               DROP_PARM_MAX_DIM_RD(I)=DROP_PARM_MAX_DIM(I)
               IF (I_DROP_PARAMETRIZATION(I)
     &            .EQ.IP_SLINGO_SCHRECKER) THEN
                  N_PARAMETER=6
               ELSE IF (I_DROP_PARAMETRIZATION(I)
     &            .EQ.IP_ACKERMAN_STEPHENS) THEN
                  N_PARAMETER=9
               ELSE IF (I_DROP_PARAMETRIZATION(I)
     &            .EQ.IP_DROP_PADE_2) THEN
                  N_PARAMETER=16
               ENDIF
!
               DO J=1, N_PARAMETER
                  DO K=1, N_BAND
                     DROP_PARAMETER_LIST_RD(J, K, I)
     &                  =DROP_PARAMETER_LIST(J, K, I)
                  ENDDO
               ENDDO
            ELSE
               L_DROP_TYPE_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 11:
      IF (L_PRESENT(11)) THEN
         L_PRESENT_RD(11)=.TRUE.
         DO I=1, N_AEROSOL_RETAIN
            I_INITIAL=INDEX_AEROSOL_RETAIN(I)
            IF (L_AEROSOL_SPECIES(I_INITIAL)) THEN
               L_AEROSOL_SPECIES_RD(I)=.TRUE.
               I_AEROSOL_PARAMETRIZATION_RD(I)
     &            =I_AEROSOL_PARAMETRIZATION(I_INITIAL)
               IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)
     &            .EQ.IP_AEROSOL_PARAM_DRY) THEN
                  NHUMIDITY_RD(I)=0
                  DO K=1, N_BAND
                     AEROSOL_ABSORPTION_RD(1, I, K)
     &                  =AEROSOL_ABSORPTION(1, I_INITIAL, K)
                     AEROSOL_SCATTERING_RD(1, I, K)
     &                  =AEROSOL_SCATTERING(1, I_INITIAL, K)
                     AEROSOL_ASYMMETRY_RD(1, I, K)
     &                  =AEROSOL_ASYMMETRY(1, I_INITIAL, K)
                  ENDDO
               ELSE IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)
     &            .EQ.IP_AEROSOL_PARAM_MOIST) THEN
                  INDEX_WATER_RD=INDEX_WATER
                  NHUMIDITY_RD(I)=NHUMIDITY(I_INITIAL)
                  DO J=1, NHUMIDITY(I_INITIAL)
                     HUMIDITIES_RD(J, I)=HUMIDITIES(J, I_INITIAL)
                     DO K=1, N_BAND
                        AEROSOL_ABSORPTION_RD(J, I, K)
     &                     =AEROSOL_ABSORPTION(J, I_INITIAL, K)
                        AEROSOL_SCATTERING_RD(J, I, K)
     &                     =AEROSOL_SCATTERING(J, I_INITIAL, K)
                        AEROSOL_ASYMMETRY_RD(J, I, K)
     &                     =AEROSOL_ASYMMETRY(J, I_INITIAL, K)
                     ENDDO
                  ENDDO
               ENDIF
!
            ELSE
               L_AEROSOL_SPECIES_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 12:
      IF (L_PRESENT(12)) THEN
         L_PRESENT_RD(12)=.TRUE.
         DO I=1, NPDR_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               L_ICE_TYPE_RD(I)=.TRUE.
               ICE_PARM_MIN_DIM_RD(I)=ICE_PARM_MIN_DIM(I)
               ICE_PARM_MAX_DIM_RD(I)=ICE_PARM_MAX_DIM(I)
!
               I_ICE_PARAMETRIZATION_RD(I)=I_ICE_PARAMETRIZATION(I)
               IF (I_ICE_PARAMETRIZATION(I)
     &            .EQ.IP_SLINGO_SCHRECKER_ICE) THEN
                  N_PARAMETER=6
               ELSE IF (I_ICE_PARAMETRIZATION(I)
     &            .EQ.IP_SUN_SHINE_VN2_VIS) THEN
                  N_PARAMETER=6
               ELSE IF (I_ICE_PARAMETRIZATION(I)
     &            .EQ.IP_SUN_SHINE_VN2_IR) THEN
                  N_PARAMETER=0
               ELSE IF (I_ICE_PARAMETRIZATION(I)
     &            .EQ.IP_ICE_ADT) THEN
                  N_PARAMETER=30
               ENDIF
!
               DO J=1, N_PARAMETER
                  DO K=1, N_BAND
                     ICE_PARAMETER_LIST_RD(J, K, I)
     &                  =ICE_PARAMETER_LIST(J, K, I)
                  ENDDO
               ENDDO
            ELSE
               L_ICE_TYPE_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 13:
      IF (L_PRESENT(13)) THEN
         L_PRESENT_RD(13)=.TRUE.
         DO I=1, N_ABSORB
            IF (L_RETAIN_ABSORB(I)) THEN
               L_DOPPLER_PRESENT_RD(COMPRESSED_INDEX(I))
     &            =L_DOPPLER_PRESENT(I)
               IF (L_DOPPLER_PRESENT(I))
     &            DOPPLER_CORRECTION_RD(COMPRESSED_INDEX(I))
     &               =DOPPLER_CORRECTION(I)
            ENDIF
         ENDDO
      ENDIF
!
!
!     BLOCK 14:
      IF (L_PRESENT(14)) THEN
         L_PRESENT_RD(14)=.TRUE.
         DO I=1, N_BAND
             N_BAND_EXCLUDE_RD(I)=N_BAND_EXCLUDE(I)
             DO J=1, N_BAND_EXCLUDE(I)
                INDEX_EXCLUDE_RD(J, I)=INDEX_EXCLUDE(J, I)
             ENDDO
         ENDDO
      ENDIF
!
!
!
      RETURN
      END
