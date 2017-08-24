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
!+ Subroutine to calculate fluxes using equivalent extinction.
!
! Method:
!       For each minor gas an equivalent extinction is calculated
!       from a clear-sky calculation. These equivalent extinctions
!       are then used in a full calculation involving the major gas.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             04-03-96                Calculation of equival-
!                                               ent extinction removed
!                                               from within the loop
!                                               over the major gas since
!                                               no dependence on
!                                               the ESFT
!                                               coefficient of the major
!                                               gas is currently used.
!                                               Half-precision exponent-
!                                               ial introduced.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code for vertically
!                                               coherent convective
!                                               cloud added.
!                                               (J. M. Edwards)
!       4.2             Oct. 96     T3E migration: EXPHF function
!                                   replaced.   (S.J.Swarbrick)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_K_EQV(IERR
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &   , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Properties
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
     &   , PLANCK_SOURCE_BAND
     &   , DIFF_PLANCK_BAND
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Surface Properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, THERMAL_GROUND_BAND
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
     &   , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
!                       Flags for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
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
!     MODULE TO SET THE DIFFUSIVITY FACTOR FOR USE WITH
!     EQUIVALENT EXTINCTION
!
      REAL
     &     DIFFUSIVITY_FACTOR_MINOR
!             MINOR DIFFUSIVITY FACTOR
!
      PARAMETER(
     &     DIFFUSIVITY_FACTOR_MINOR=1.66E+00
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
     &     PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE IN BAND
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
!             COLUMN FLAGS FOR COLUMNS
      REAL  !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FRACTION
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY FRACTION
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER FOR ENERGY AT INTERFACES
     &   , AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREAS OF COLUMNS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!                       Cloudy Optical Properties
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
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX IN BAND
!
!                       Flags for Clear-sky Fluxes
      LOGICAL   !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE CLEAR-SKY PROPERTIES
      INTEGER   !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Clear-sky Fluxes Calculated
      REAL  !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR-SKY TOTAL FLUX IN BAND
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
      INTEGER
     &     I_GAS
!             INDEX OF MAIN GAS
     &   , I_GAS_BAND
!             INDEX OF ACTIVE GAS
     &   , I_GAS_POINTER(NPD_SPECIES)
!             POINTER ARRAY FOR MONOCHROMATIC ESFTs
     &   , IEX
!             INDEX OF ESFT TERM
      REAL
     &     SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD FLUX
     &   , ESFT_WEIGHT
!             ESFT WEIGHT FOR CURRENT CALCULATION
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR TRANSMISSION TO INCLUDE EFFECTS
!             OF MINOR GASES AND TAKE OUT EQUIVALENT EXTINCTION
     &   , K_EQV(NPD_PROFILE, NPD_LAYER)
!             EQUIVALENT EXTINCTION
     &   , TAU_GAS(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH OF GAS
     &   , K_ESFT_MONO(NPD_SPECIES)
!             MONOCHROMATIC EXPONENTS
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS EXTINCTION
      REAL
     &     FLUX_DIRECT_PART(NPD_PROFILE, 0: NPD_LAYER)
!             PARTIAL DIRECT FLUX
     &   , FLUX_TOTAL_PART(NPD_PROFILE, 2*NPD_LAYER+2)
!             PARTIAL TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR_PART(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR PARTIAL DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR_PART(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR PARTIAL TOTAL FLUX
     &   , SUM_FLUX(NPD_PROFILE, 2*NPD_LAYER+2, NPD_SPECIES)
!             SUM OF FLUXES FOR WEIGHTING
     &   , SUM_K_FLUX(NPD_PROFILE, 2*NPD_LAYER+2, NPD_SPECIES)
!             SUM OF K*FLUXES FOR WEIGHTING
     &   , FLUX_TERM(NPD_PROFILE, 0: NPD_LAYER)
!             FLUX WITH ONE TERM
     &   , FLUX_GAS(NPD_PROFILE, 0: NPD_LAYER)
!             FLUX WITH ONE GAS
      REAL
     &     MEAN_NET_FLUX
!             MEAN NET FLUX
     &   , MEAN_K_NET_FLUX
!             MEAN K-WEIGHTED NET FLUX
     &   , K_WEAK
!             WEAK ABSORPTION FOR MINOR GAS
      REAL
     &     KE_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)
!             EQUIVALENT FREE ABSORPTIVE EXTINCTION
     &   , KE_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             EQUIVALENT FREE ABSORPTIVE EXTINCTION
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     SCALE_ABSORB, GAS_OPTICAL_PROPERTIES
     &   , MONOCHROMATIC_GAS_FLUX, MONOCHROMATIC_FLUX
     &   , AUGMENT_FLUX
!
!
!
      I_GAS=INDEX_ABSORB(1, I_BAND)
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
!
!        AN APPROPRIATE SCALING FACTOR IS CALCULATED FOR THE DIRECT
!        BEAM, WHILST THE EQUIVALENT EXTINCTION FOR THE DIFFUSE BEAM
!        IS WEIGHTED WITH THE SOLAR SCALING FACTOR AS EVALUATED
!        AT THE SURFACE.
!
!        INITIALIZE THE SCALING FACTORS:
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               ADJUST_SOLAR_KE(L, I)=1.0E+00
               K_EQV(L, I)=0.0E+00
            ENDDO
         ENDDO
!
         DO J=2, N_GAS
!
!           INITIALIZE THE NORMALIZED FLUX FOR THE GAS.
            DO L=1, N_PROFILE
               FLUX_GAS(L, 0)=1.0E+00
            ENDDO
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_GAS(L, I)=0.0E+00
                  SUM_K_FLUX(L, N_LAYER, J)=0.0E+00
                  SUM_FLUX(L, N_LAYER, J)=0.0E+00
               ENDDO
            ENDDO
!
            I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
            DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS_BAND)
!
!              STORE THE ESFT WEIGHT FOR FUTURE USE.
               ESFT_WEIGHT=W_ESFT(IEX, I_BAND,  I_GAS_BAND)
!
!              RESCALE THE AMOUNT OF GAS FOR THIS ABSORBER IF REQUIRED.
               IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND).EQ.IP_SCALE_TERM)
     &            THEN
                  CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &               , GAS_MIX_RATIO(1, 0, I_GAS_BAND), P, T
     &               , L_LAYER, I_TOP
     &               , GAS_FRAC_RESCALED(1, 0, I_GAS_BAND)
     &               , I_SCALE_FNC(I_BAND, I_GAS_BAND)
     &               , P_REFERENCE(I_GAS_BAND, I_BAND)
     &               , T_REFERENCE(I_GAS_BAND, I_BAND)
     &               , SCALE_VECTOR(1, IEX, I_BAND, I_GAS_BAND)
     &               , L_DOPPLER(I_GAS_BAND)
     &               , DOPPLER_CORRECTION(I_GAS_BAND)
     &               , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &               , NPD_SCALE_VARIABLE
     &               )
                  IF (IERR.NE.I_NORMAL) RETURN
               ENDIF
!
               DO L=1, N_PROFILE
                  FLUX_TERM(L, 0)=ESFT_WEIGHT
               ENDDO
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     FLUX_TERM(L, I)=FLUX_TERM(L, I-1)
     &                  *EXP(-K_ESFT(IEX, I_BAND, I_GAS_BAND)          
     &                  *GAS_FRAC_RESCALED(L, I, I_GAS_BAND)
     &                  *D_MASS(L, I)*SEC_0(L))
                     FLUX_GAS(L, I)=FLUX_GAS(L, I)+FLUX_TERM(L, I)
                  ENDDO
               ENDDO
!
!              CALCULATE THE INCREMENT IN THE ABSORPTIVE EXTINCTION
               DO L=1, N_PROFILE
                  SUM_K_FLUX(L, N_LAYER, J)
     &               =SUM_K_FLUX(L, N_LAYER, J)
     &               +K_ESFT(IEX, I_BAND, I_GAS_BAND)
     &               *FLUX_TERM(L, N_LAYER)
                  SUM_FLUX(L, N_LAYER, J)
     &               =SUM_FLUX(L, N_LAYER, J)+FLUX_TERM(L, N_LAYER)
               ENDDO
!
            ENDDO
!
!           SET THE EQUIVALENT EXTINCTION FOR THE DIFFUSE BEAM,
!           WEIGHTING WITH THE DIRECT SURFACE FLUX.
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (SUM_FLUX(L, N_LAYER, J).GT.0.0E+00) THEN
                     K_EQV(L, I)=K_EQV(L, I)
     &                  +GAS_FRAC_RESCALED(L, I, I_GAS_BAND)
     &                  *SUM_K_FLUX(L, N_LAYER, J)
     &                  /SUM_FLUX(L, N_LAYER, J)
                  ELSE
!                    THIS CASE CAN ARISE ONLY WHEN THE SUN IS CLOSE
!                    TO THE HORIZON WHEN THE EXPONENTIAL MAY UNDERFLOW
!                    TO 0. WE USE THE WEAKEST ESFT-TERM.
                     K_EQV(L, I)=K_EQV(L, I)
     &                  *K_ESFT(1, I_BAND, I_GAS_BAND)
     &                  *GAS_FRAC_RESCALED(L, I, I_GAS_BAND)
                  ENDIF
                  IF (FLUX_GAS(L, I-1).GT.0.0E+00) THEN
!                    IF THE FLUX HAS BEEN REDUCED TO 0 AT THE UPPER
!                    LEVEL THE ADJUSTING FACTOR IS NOT OF IMPORTANCE
!                    AND NEED NOT BE ADJUSTED. THIS WILL PREVENT
!                    POSSIBLE FAILURES.
                     ADJUST_SOLAR_KE(L, I)
     &                  =ADJUST_SOLAR_KE(L, I)*FLUX_GAS(L, I)
     &                  /FLUX_GAS(L, I-1)
                  ENDIF

               ENDDO
            ENDDO
!
         ENDDO
!
!        SINCE THE GREY EXTINCTION WILL LATER BE MODIFIED WE MUST
!        INCREASE THE TRANSMISSION OF THE SOLAR BEAM TO COMPENSATE.
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               ADJUST_SOLAR_KE(L, I)=ADJUST_SOLAR_KE(L, I)
     &            *EXP(K_EQV(L, I)*D_MASS(L, I)*SEC_0(L)) 
            ENDDO
         ENDDO
!
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
!
!        EQUIVALENT ABSORPTION IS USED FOR THE MINOR GASES.
!
         DO J=2, N_GAS
!
!
!           INITIALIZE THE SUMS TO FORM THE RATIO TO 0.
            DO I=1, 2*N_LAYER+2
               DO L=1, N_PROFILE
                  SUM_FLUX(L, I, J)=0.0E+00
                  SUM_K_FLUX(L, I, J)=0.0E+00
               ENDDO
            ENDDO
!
            I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
            DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS_BAND)
!
!              STORE THE ESFT WEIGHT FOR FUTURE USE.
               ESFT_WEIGHT=W_ESFT(IEX, I_BAND,  I_GAS_BAND)
!
!
!              RESCALE THE AMOUNT OF GAS FOR THIS ABSORBER IF REQUIRED.
               IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND).EQ.IP_SCALE_TERM)
     &            THEN
                  CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &               , GAS_MIX_RATIO(1, 0, I_GAS_BAND), P, T
     &               , L_LAYER, I_TOP
     &               , GAS_FRAC_RESCALED(1, 0, I_GAS_BAND)
     &               , I_SCALE_FNC(I_BAND, I_GAS_BAND)
     &               , P_REFERENCE(I_GAS_BAND, I_BAND)
     &               , T_REFERENCE(I_GAS_BAND, I_BAND)
     &               , SCALE_VECTOR(1, IEX, I_BAND, I_GAS_BAND)
     &               , L_DOPPLER(I_GAS_BAND)
     &               , DOPPLER_CORRECTION(I_GAS_BAND)
     &               , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &               , NPD_SCALE_VARIABLE
     &               )
                  IF (IERR.NE.I_NORMAL) RETURN
               ENDIF
!
!              SET THE APPROPRIATE BOUNDARY TERMS FOR THE
!              TOTAL UPWARD AND DOWNWARD FLUXES AT THE BOUNDARIES.
!
               DO L=1, N_PROFILE
                  FLUX_INC_DIRECT(L)=0.0E+00
                  FLUX_DIRECT_PART(L, N_LAYER)=0.0E+00
                  FLUX_INC_DOWN(L)=-PLANCK_SOURCE_BAND(L, 0)
                  SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)
     &               -(1.0E+00-ALBEDO_SURFACE_DIFF(L))
     &               *PLANCK_SOURCE_BAND(L, N_LAYER)
               ENDDO
!
!              SET THE OPTICAL DEPTHS OF EACH LAYER.
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     TAU_GAS(L, I)=K_ESFT(IEX, I_BAND, I_GAS_BAND)
     &                  *GAS_FRAC_RESCALED(L, I, I_GAS_BAND)
     &                  *D_MASS(L, I)
                  ENDDO
               ENDDO
!
!              CALCULATE THE FLUXES WITH JUST THIS GAS.
               CALL MONOCHROMATIC_GAS_FLUX(N_PROFILE, N_LAYER
     &            , .FALSE.
     &            , TAU_GAS
     &            , ISOLIR, SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN
     &            , DIFF_PLANCK_BAND, SOURCE_GROUND
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &            , DIFFUSIVITY_FACTOR_MINOR
     &            , FLUX_DIRECT_PART, FLUX_TOTAL_PART
     &            , NPD_PROFILE, NPD_LAYER
     &            )
!
               DO I=1, 2*N_LAYER+2
                  DO L=1, N_PROFILE
                     SUM_K_FLUX(L, I, J)=SUM_K_FLUX(L, I, J)
     &                  +K_ESFT(IEX, I_BAND, I_GAS_BAND)
     &                  *ESFT_WEIGHT*FLUX_TOTAL_PART(L, I)
                     SUM_FLUX(L, I, J)=SUM_FLUX(L, I, J)
     &                  +ESFT_WEIGHT*FLUX_TOTAL_PART(L, I)
                  ENDDO
               ENDDO
!
            ENDDO
!
         ENDDO
!
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               K_EQV(L, I)=0.0E+00
            ENDDO
         ENDDO
!
         DO J=2, N_GAS
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  MEAN_K_NET_FLUX=0.5E+00*(SUM_K_FLUX(L, 2*I, J)
     &               +SUM_K_FLUX(L, 2*I+2, J)
     &               -SUM_K_FLUX(L, 2*I-1, J)
     &               -SUM_K_FLUX(L, 2*I+1, J))
                  MEAN_NET_FLUX=0.5E+00*(SUM_FLUX(L, 2*I, J)
     &               +SUM_FLUX(L, 2*I+2, J)
     &               -SUM_FLUX(L, 2*I-1, J)
     &               -SUM_FLUX(L, 2*I+1, J))
!                 NEGATIVE EFFECTIVE EXTINCTIONS MUST BE REMOVED.
                  K_WEAK=MAX(0.0E+00, MEAN_K_NET_FLUX/MEAN_NET_FLUX)
                  K_EQV(L, I)=K_EQV(L, I)
     &               +K_WEAK*GAS_FRAC_RESCALED(L, I
     &               , INDEX_ABSORB(J, I_BAND))
               ENDDO
            ENDDO
         ENDDO

      ENDIF
!
!
!     THE ESFT TERMS FOR THE MAJOR GAS IN THE BAND ARE USED WITH
!     APPROPRIATE WEIGHTED TERMS FOR THE MINOR GASES.
      I_GAS_POINTER(1)=I_GAS
      DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS)
!
!        STORE THE ESFT WEIGHT FOR FUTURE USE.
         ESFT_WEIGHT=W_ESFT(IEX, I_BAND,  I_GAS)
!
!        RESCALE FOR EACH ESFT TERM IF THAT IS REQUIRED.
         IF (I_SCALE_ESFT(I_BAND, I_GAS).EQ.IP_SCALE_TERM) THEN
            CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &         , GAS_MIX_RATIO(1, 0, I_GAS), P, T
     &         , L_LAYER, I_TOP
     &         , GAS_FRAC_RESCALED(1, 0, I_GAS)
     &         , I_SCALE_FNC(I_BAND, I_GAS)
     &         , P_REFERENCE(I_GAS, I_BAND)
     &         , T_REFERENCE(I_GAS, I_BAND)
     &         , SCALE_VECTOR(1, IEX, I_BAND, I_GAS)
     &         , L_DOPPLER(I_GAS), DOPPLER_CORRECTION(I_GAS)
     &         , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &         , NPD_SCALE_VARIABLE
     &         )
            IF (IERR.NE.I_NORMAL) RETURN
         ENDIF
!
!        SET THE APPROPRIATE BOUNDARY TERMS FOR THE TOTAL
!        UPWARD AND DOWNWARD FLUXES.
!
         IF (ISOLIR.EQ.IP_SOLAR) THEN
!           VISIBLE REGION.
            DO L=1, N_PROFILE
               SOURCE_GROUND(L)=0.0E+00
               FLUX_INC_DOWN(L)=SOLAR_FLUX(L)
               FLUX_INC_DIRECT(L)=SOLAR_FLUX(L)
            ENDDO
         ELSEIF (ISOLIR.EQ.IP_INFRA_RED) THEN
!           INFRA-RED REGION.
            DO L=1, N_PROFILE
               FLUX_INC_DIRECT(L)=0.0E+00
               FLUX_DIRECT_PART(L, N_LAYER)=0.0E+00
               FLUX_INC_DOWN(L)=-PLANCK_SOURCE_BAND(L, 0)
               SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)
     &            -(1.0E+00-ALBEDO_SURFACE_DIFF(L))
     &            *PLANCK_SOURCE_BAND(L, N_LAYER)
            ENDDO
            IF (L_CLEAR) THEN
               DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR_PART(L, N_LAYER)=0.0E+00
               ENDDO
            ENDIF
         ENDIF
!
!
!        AUGMENT THE GREY EXTINCTION WITH AN EFFECTIVE VALUE
!        FOR EACH GAS.
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               KE_GREY_TOT_FREE(L, I)=K_GREY_TOT_FREE(L, I)
            ENDDO
         ENDDO
         IF (L_CLOUD) THEN
            DO K=1, N_CLOUD_TYPE
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     KE_GREY_TOT_CLOUD(L, I, K)
     &                  =K_GREY_TOT_CLOUD(L, I, K)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               KE_GREY_TOT_FREE(L, I)=KE_GREY_TOT_FREE(L, I)
     &            +K_EQV(L, I)
            ENDDO
         ENDDO
         IF (L_CLOUD) THEN
            DO K=1, N_CLOUD_TYPE
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     KE_GREY_TOT_CLOUD(L, I, K)
     &                  =KE_GREY_TOT_CLOUD(L, I, K)+K_EQV(L, I)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
!        ASSIGN THE MONOCHROMATIC ABSORPTION COEFFICIENT.
         K_ESFT_MONO(I_GAS)=K_ESFT(IEX, I_BAND, I_GAS)
!
         CALL GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER
     &      , 1, I_GAS_POINTER, K_ESFT_MONO
     &      , GAS_FRAC_RESCALED
     &      , K_GAS_ABS
     &      , NPD_PROFILE, NPD_LAYER, NPD_SPECIES
     &      )
!
!
         CALL MONOCHROMATIC_FLUX(IERR
!                       Atmospheric Properties
     &      , N_PROFILE, N_LAYER, D_MASS
!                       Angular Integration
     &      , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT
     &      , L_RESCALE, N_ORDER_GAUSS
!                       Treatment of Scattering
     &      , I_SCATTER_METHOD_BAND
!                       Options for Solver
     &      , I_SOLVER, L_NET, N_AUGMENT
!                       Gaseous Propreties
     &      , K_GAS_ABS
!                       Options for Equivalent Extinction
     &      , .TRUE., ADJUST_SOLAR_KE
!                       Spectral Region
     &      , ISOLIR
!                       Infra-red Properties
     &      , DIFF_PLANCK_BAND
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2
!                       Conditions at TOA
     &      , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN
!                       Surface Properties
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
     &      , THERMAL_GROUND_BAND
!                       Clear-sky Optical Properties
     &      , KE_GREY_TOT_FREE, K_EXT_SCAT_FREE
     &      , ASYMMETRY_FREE, FORWARD_SCATTER_FREE
!                       Cloudy Properties
     &      , L_CLOUD, I_CLOUD
!                       Cloud Geometry
     &      , N_CLOUD_TOP
     &      , N_CLOUD_TYPE, FRAC_CLOUD
     &      , I_REGION_CLOUD, FRAC_REGION
     &      , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &      , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &      , CLOUD_OVERLAP
     &      , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &      , KE_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &      , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
!                       Flxues Calculated
     &      , FLUX_DIRECT_PART, FLUX_TOTAL_PART
!                       Flags for Clear-sky Calculations
     &      , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &      , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART
!                       Planckian Function
     &      , PLANCK_SOURCE_BAND
!                       Dimensions of Arrays
     &      , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
!
!        INCREMENT THE FLUXES WITHIN THE BAND.
         CALL AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &      , ISOLIR, L_CLEAR
     &      , ESFT_WEIGHT
     &      , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND
     &      , FLUX_DIRECT_PART, FLUX_TOTAL_PART
     &      , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND
     &      , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDDO
!
!
      RETURN
      END
