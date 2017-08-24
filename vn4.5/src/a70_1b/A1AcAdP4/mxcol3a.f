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
!+ Subroutine to solve the two-stream equations in a mixed column.
!
! Method:
!       The two-stream coefficients are calculated in clear regions
!       and in stratiform and convective clouds. From these
!       coefficients transmission and reflection coefficients are
!       determined. The coefficients for convective and stratiform
!       clouds are appropriately mixed to form single cloudy values
!       and an appropriate solver is called.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-04-96                New solver added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for obsolete
!                                               solver removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIX_COLUMN(IERR
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &   , I_2STREAM
!                       Corrections to Two-stream Equations
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &   , I_SOLVER, L_NET
!                       Options for Equivalent Extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &   , ISOLIR
!                       Infra-red Properties
     &   , DIFF_PLANCK
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &   , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Conditions at Surface
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Clear-sky Single Scattering Properties
     &   , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE
!                       Cloud Geometry
     &   , N_CLOUD_TOP
     &   , N_CLOUD_TYPE, FRAC_CLOUD
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &   , CLOUD_OVERLAP
!                       Cloudy Optical Properties
     &   , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculations
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER
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
!
!     INCLUDE COMDECKS.
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
!     MODULE TO SET POINTERS IN CLOUD_OVERLAP.
!
!     NOTE THAT SEVERAL POINTERS ARE IDENTICAL SINCE ONLY CERTAIN
!     GROUPS OF COEFFICIENTS ARE RELEVANT TO A PARTICULAR SCHEME.
!
      INTEGER
     &     IP_CLOVLP_GFF
!             POINTER TO GAMMA-FREE-FREE
     &   , IP_CLOVLP_GFC
!             POINTER TO GAMMA-FREE-CLOUD
     &   , IP_CLOVLP_GCF
!             POINTER TO GAMMA-CLOUD-FREE
     &   , IP_CLOVLP_GCC
!             POINTER TO GAMMA-CLOUD-CLOUD
     &   , IP_CLOVLP_BFF
!             POINTER TO BETA-FREE-FREE
     &   , IP_CLOVLP_BFC
!             POINTER TO BETA-FREE-CLOUD
     &   , IP_CLOVLP_BCF
!             POINTER TO BETA-CLOUD-FREE
     &   , IP_CLOVLP_BCC
!             POINTER TO BETA-CLOUD-CLOUD
     &   , IP_CLOVLP_GFM
!             POINTER TO GAMMA_F-
     &   , IP_CLOVLP_GFP
!             POINTER TO GAMMA_F+
     &   , IP_CLOVLP_BFM
!             POINTER TO BETA_F-
     &   , IP_CLOVLP_BFP
!             POINTER TO BETA_F+
     &   , IP_CLOVLP_GM
!             POINTER TO GAMMA_-
     &   , IP_CLOVLP_GP
!             POINTER TO GAMMA_+
     &   , IP_CLOVLP_BM
!             POINTER TO BETA_-
     &   , IP_CLOVLP_BP
!             POINTER TO BETA_+
!
!     POINTERS FOR TRIPLE OVERLAPS:
      INTEGER
     &     IP_CLOVLP_V11
     &   , IP_CLOVLP_V12
     &   , IP_CLOVLP_V13
     &   , IP_CLOVLP_V21
     &   , IP_CLOVLP_V22
     &   , IP_CLOVLP_V23
     &   , IP_CLOVLP_V31
     &   , IP_CLOVLP_V32
     &   , IP_CLOVLP_V33
     &   , IP_CLOVLP_U11
     &   , IP_CLOVLP_U12
     &   , IP_CLOVLP_U13
     &   , IP_CLOVLP_U21
     &   , IP_CLOVLP_U22
     &   , IP_CLOVLP_U23
     &   , IP_CLOVLP_U31
     &   , IP_CLOVLP_U32
     &   , IP_CLOVLP_U33
!
      PARAMETER(
     &     IP_CLOVLP_GFF=1
     &   , IP_CLOVLP_GFC=2
     &   , IP_CLOVLP_GCF=3
     &   , IP_CLOVLP_GCC=4
     &   , IP_CLOVLP_BFF=5
     &   , IP_CLOVLP_BFC=6
     &   , IP_CLOVLP_BCF=7
     &   , IP_CLOVLP_BCC=8
     &   , IP_CLOVLP_GFM=5
     &   , IP_CLOVLP_GFP=6
     &   , IP_CLOVLP_BFM=7
     &   , IP_CLOVLP_BFP=8
     &   , IP_CLOVLP_GM=5
     &   , IP_CLOVLP_GP=6
     &   , IP_CLOVLP_BM=7
     &   , IP_CLOVLP_BP=8
     &   )
!
      PARAMETER(
     &     IP_CLOVLP_V11=1
     &   , IP_CLOVLP_V12=2
     &   , IP_CLOVLP_V13=3
     &   , IP_CLOVLP_V21=4
     &   , IP_CLOVLP_V22=5
     &   , IP_CLOVLP_V23=6
     &   , IP_CLOVLP_V31=7
     &   , IP_CLOVLP_V32=8
     &   , IP_CLOVLP_V33=9
     &   , IP_CLOVLP_U11=10
     &   , IP_CLOVLP_U12=11
     &   , IP_CLOVLP_U13=12
     &   , IP_CLOVLP_U21=13
     &   , IP_CLOVLP_U22=14
     &   , IP_CLOVLP_U23=15
     &   , IP_CLOVLP_U31=16
     &   , IP_CLOVLP_U32=17
     &   , IP_CLOVLP_U33=18
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
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
     &   , ISOLIR
!             SPECTRAL REGION
     &   , I_2STREAM
!             TWO-STREAM SCHEME
     &   , I_SOLVER
!             SOLVER USED
     &   , I_SOLVER_CLEAR
!             SOLVER FOR CLEAR-SKY FLUXES
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES
     &   , L_CLEAR
!             CALCULATE CLEAR-SKY FLUXES
     &   , L_SCALE_SOLAR
!             FLAG TO SCALE SOLAR
     &   , L_IR_SOURCE_QUAD
!             USE QUADRATIC SOURCE TERM
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION TO 2-STREAM
!
!     OPTICAL PROPERTIES:
      REAL      !, INTENT(IN)
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY ASYMMETRY
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ASYMMETRY
!
!     CLOUD GEOMETRY:
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FRACTIONS IN EACH LAYER
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR SKY FRACTIONS IN EACH LAYER
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             ENERGY TRANSFER COEFFICIENTS
      REAL      !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF SOLAR ZENITH ANGLE
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             CHANGE IN PLANCK FUNCTION
     &   , SOURCE_GROUND(NPD_PROFILE)
!             FLUX FROM SURFACE
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCE OF PLANCKIAN
     &   , PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE FUNCTION
     &   , GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
!     FLUXES CALCULATED
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             LONG FLUX VECTOR
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
!
!
!
!     LOCAL VARIABALES.
      INTEGER
     &     N_SOURCE_COEFF
!             NUMBER OF SOURCE COEFFICIENTS
     &   , N_EQUATION
!             NUMBER OF EQUATIONS
     &   , I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!     CLEAR-SKY COEFFICIENTS:
      REAL
     &     TRANS_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE TRANSMISSION OF LAYER
     &   , REFLECT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE REFLECTANCE OF LAYER
     &   , TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DIRECT TRANSMISSION OF LAYER
     &   , SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             FREE SOURCE COEFFICIENTS
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DOWNWARD SOURCE
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE UPWARD SOURCE
     &   , S_DOWN_CLEAR(NPD_PROFILE, NPD_LAYER)
!             CLEAR DOWNWARD SOURCE
     &   , S_UP_CLEAR(NPD_PROFILE, NPD_LAYER)
!             CLEAR UPWARD SOURCE
!
!     CLOUDY COEFFICIENTS:
      REAL
     &     TRANS_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY TRANSMISSION OF LAYER
     &   , REFLECT_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY REFLECTANCE OF LAYER
     &   , TRANS_0_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DIRECT TRANSMISSION OF LAYER
     &   , SOURCE_COEFF_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             CLOUDY SOURCE COEFFICIENTS
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DOWNWARD SOURCE
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY UPWARD SOURCE
!
!     SOURCE FUNCTIONS AT THE CROUND
      REAL
     &     SOURCE_GROUND_FREE(NPD_PROFILE)
!             SOURCE FROM GROUND UNDER CLEAR SKIES
     &   , SOURCE_GROUND_CLOUD(NPD_PROFILE)
!             SOURCE FROM GROUND UNDER CLOUDY SKIES
     &   , FLUX_DIRECT_GROUND_CLOUD(NPD_PROFILE)
!             DIRECT FLUX AT GROUND UNDER CLOUDY SKIES
!
!
!     NUMERICAL ARRAYS:
      REAL
     &     A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             RHS OF MATRIX EQUATION
     &   , WORK(NPD_PROFILE)
!             WORKING ARRAY FOR SOLVER
!
!     FUNCTIONS CALLED:
      INTEGER
     &     SET_N_SOURCE_COEFF
!             FUNCTION TO SET NUMBER OF SOURCE COEFFICIENTS
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     TWO_COEFF, TWO_COEFF_CLOUD, IR_SOURCE, MIXED_SOLAR_SOURCE
     &   , BAND_SOLVER, MIX_COLUMN_FULL, MIX_APP_SCAT
     &   , CLEAR_SUPPLEMENT
!
!
!
!     CALCULATE THE TRANSMISSION AND REFLECTION COEFFICIENTS AND
!     SOURCE TERMS FOR THE CLEAR AND CLOUDY PARTS OF THE COLUMN
!
!     SET THE NUMBER OF SOURCE COEFFICIENTS FOR THE APPROXIMATION
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
      CALL TWO_COEFF(IERR
     &   , N_PROFILE, 1, N_LAYER
     &   , I_2STREAM, L_IR_SOURCE_QUAD
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE
     &   , ISOLIR, SEC_0
     &   , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE
     &   , SOURCE_COEFF_FREE
     &   , NPD_PROFILE, NPD_LAYER
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!     INFRA-RED SOURCE TERMS DEPEND ONLY ON THE LAYER AND MAY BE
!     CALCULATED NOW. SOLAR TERMS DEPEND ON CONDITIONS IN CLOUD
!     IN OVERLYING LAYERS AND MUST BE CALCULATED LATER.
!
      IF (ISOLIR.EQ.IP_INFRA_RED) THEN
!
         CALL IR_SOURCE(N_PROFILE, 1, N_LAYER
     &      , SOURCE_COEFF_FREE, DIFF_PLANCK
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &      , GROUND_EMISSION, N_LAYER
     &      , TAU_FREE, TRANS_FREE
     &      , S_DOWN_FREE, S_UP_FREE
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
!        IF A CLEAR-SKY CALCULATION IS REQUIRED THESE SOURCE TERMS MUST
!        BE STORED.
         IF (L_CLEAR) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  S_DOWN_CLEAR(L, I)=S_DOWN_FREE(L, I)
                  S_UP_CLEAR(L, I)=S_UP_FREE(L, I)
               ENDDO
            ENDDO
         ENDIF
!
!        SCALE THE SOURCES BY THE CLEAR-SKY FRACTIONS IN THE CLOUDY
!        LAYERS. IN HIGHER LAYERS THE CLEAR-SKY FRACTION IS 1.
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               S_DOWN_FREE(L, I)=W_FREE(L, I)*S_DOWN_FREE(L, I)
               S_UP_FREE(L, I)=W_FREE(L, I)*S_UP_FREE(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     REPEAT THE CALCULATION FOR CLOUDY REGIONS.
!
!
      CALL TWO_COEFF_CLOUD(IERR
     &   , N_PROFILE, N_CLOUD_TOP, N_LAYER
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF
     &   , N_CLOUD_TYPE, FRAC_CLOUD
     &   , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD
     &   , ISOLIR, SEC_0
     &   , TRANS_CLOUD, REFLECT_CLOUD, TRANS_0_CLOUD
     &   , SOURCE_COEFF_CLOUD
     &   , NPD_PROFILE, NPD_LAYER
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
      IF (ISOLIR.EQ.IP_INFRA_RED) THEN
!
!        EDGE CORRECTIONS FOR THE TWO-STREAM EQUATIONS DO NOT
!        REALLY FIT WITH THIS METHOD OF TREATING CLOUDS. OPTICAL
!        DEPTHS AND TRANSMISSIONS MUST BE PASSED TO THE SUBROUTINE
!        TO FILL THE ARGUMENT LIST, BUT IT IS NOT INTENDED THAT
!        THESE ARRAYS WILL BE USED.
!
         CALL IR_SOURCE(N_PROFILE, N_CLOUD_TOP, N_LAYER
     &      , SOURCE_COEFF_CLOUD, DIFF_PLANCK
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &      , GROUND_EMISSION, N_LAYER
     &      , TAU_CLOUD, TRANS_CLOUD
     &      , S_DOWN_CLOUD, S_UP_CLOUD
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
!
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               S_DOWN_CLOUD(L, I)=W_CLOUD(L, I)*S_DOWN_CLOUD(L, I)
               S_UP_CLOUD(L, I)=W_CLOUD(L, I)*S_UP_CLOUD(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     CALCULATE THE APPROPRIATE SOURCE TERMS FOR THE SOLAR: CLOUDY
!     AND CLEAR PROPERTIES ARE BOTH NEEDED HERE.
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
!
         CALL MIXED_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , FLUX_INC_DIRECT
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &      , TRANS_0_FREE, SOURCE_COEFF_FREE
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)
     &      , TRANS_0_CLOUD, SOURCE_COEFF_CLOUD
     &      , FLUX_DIRECT
     &      , FLUX_DIRECT_GROUND_CLOUD
     &      , S_UP_FREE, S_DOWN_FREE
     &      , S_UP_CLOUD, S_DOWN_CLOUD
     &      , NPD_PROFILE, NPD_LAYER
     &   )
      ENDIF
!
!
!
!     FORMULATE THE MATRIX EQUATIONS FOR THE FLUXES.
!
      IF (I_SOLVER.EQ.IP_SOLVER_MIX_APP_SCAT) THEN
!
         CALL MIX_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE
     &      , TRANS_CLOUD, REFLECT_CLOUD
     &      , S_DOWN_CLOUD, S_UP_CLOUD
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCC)
     &      , L_NET
     &      , FLUX_INC_DOWN
     &      , SOURCE_GROUND, ALBEDO_SURFACE_DIFF
     &      , FLUX_TOTAL
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE IF (I_SOLVER.EQ.IP_SOLVER_MIX_DIRECT) THEN
!
!        SET THE PARTITIONED SOURCE FUNCTIONS AT THE GROUND.
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            DO L=1, N_PROFILE
               SOURCE_GROUND_FREE(L)=(ALBEDO_SURFACE_DIR(L)
     &            -ALBEDO_SURFACE_DIFF(L))
     &            *(FLUX_DIRECT(L, N_LAYER)
     &            -FLUX_DIRECT_GROUND_CLOUD(L))
               SOURCE_GROUND_CLOUD(L)=(ALBEDO_SURFACE_DIR(L)
     &            -ALBEDO_SURFACE_DIFF(L))
     &            *FLUX_DIRECT_GROUND_CLOUD(L)
            ENDDO
         ELSE
            DO L=1, N_PROFILE
               SOURCE_GROUND_FREE(L)
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BFF)
     &            *SOURCE_GROUND(L)
               SOURCE_GROUND_CLOUD(L)
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BCF)
     &            *SOURCE_GROUND(L)
            ENDDO
         ENDIF
!
         CALL SOLVER_MIX_DIRECT(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE
     &      , TRANS_CLOUD, REFLECT_CLOUD
     &      , S_DOWN_CLOUD, S_UP_CLOUD
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFC)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCF)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCC)
     &      , L_NET
     &      , FLUX_INC_DOWN
     &      , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD
     &      , ALBEDO_SURFACE_DIFF
     &      , FLUX_TOTAL
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE IF (I_SOLVER.EQ.IP_SOLVER_MIX_11) THEN
!
         CALL MIX_COLUMN_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE
     &      , TRANS_CLOUD, REFLECT_CLOUD
     &      , S_DOWN_CLOUD, S_UP_CLOUD
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GM)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GP)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BM)
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BP)
     &      , FLUX_INC_DOWN
     &      , SOURCE_GROUND, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , FLUX_DIRECT(1, N_LAYER)
     &      , FLUX_TOTAL
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE
!
         WRITE(IU_ERR, '(/A)')
     &      '***ERROR: THE SOLVER SPECIFIED IS NOT VALID HERE.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      IF (L_CLEAR) THEN
!
         CALL CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER, I_SOLVER_CLEAR
     &      , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE
     &      , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN
     &      , S_DOWN_CLEAR, S_UP_CLEAR
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , SOURCE_GROUND
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &      , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!
      RETURN
      END
