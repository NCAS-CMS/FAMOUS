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
!+ Subroutine to calculate grey extinctions.
!
! Method:
!       For each activated optical process, excluding gaseous
!       absorption, increments are calculated for the total and
!       scattering extinctions, and the products of the asymmetry
!       factor and the forward scattering factor in clear and
!       cloudy regions. These increments are summed, and the grey
!       total and scattering extinctions and the asymmetry and forward
!       scattering factors are thus calculated.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             06-06-96                Indentation made
!                                               consistent
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.4             30-09-96                Effective radius
!                                               relabelled as
!                                               characteristic
!                                               dimension for
!                                               generality to cover
!                                               parametrizations
!                                               of non-spherical
!                                               ice.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Removal of test
!                                               for deleted cloud
!                                               scheme.
!                                               (J. M. Edwards)
!LL  4.5  27/04/98  Add Fujitsu vectorization directive. 
!LL                                           RBarnes@ecmwf.int
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
! Fujitsu directive to encourage vectorization for whole routine
!OCL NOVREC
      SUBROUTINE GREY_EXTINCTION(IERR
     &   , N_PROFILE, N_LAYER, L_LAYER, P, T, DENSITY
     &   , L_RESCALE
     &   , L_RAYLEIGH, RAYLEIGH_COEFF
     &   , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER, K_CONTINUUM
     &   , AMOUNT_CONTINUUM
     &   , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO
     &   , I_AEROSOL_PARAMETRIZATION
     &   , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY
     &   , MEAN_REL_HUMIDITY
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY
     &   , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP
     &   , L_CLOUD_LAYER, I_CLOUD
     &   , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR
     &   , N_CLOUD_TYPE, I_CLOUD_TYPE
     &   , K_EXT_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE
     &   , FORWARD_SCATTER_FREE
     &   , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD
     &   , NPD_PROFILE, NPD_LAYER, NPD_CONTINUUM
     &   , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES
     &   , NPD_CLOUD_PARAMETER
     &   )
!
!
!
      IMPLICIT NONE
!
!
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOLS
     &   , NPD_HUMIDITIES
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_CONTINUUM
!             MAXIMUM NUMBER OF CONTINUA
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
!
!     INCLUDE COMDECKS
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
!     MODULE TO SET INDICES FOR PHASES.
!
      INTEGER
     &     IP_PHASE_WATER
!             LIQUID PHASE
     &   , IP_PHASE_ICE
!             ICE PHASE
!
      PARAMETER(
     &     IP_PHASE_WATER=1
     &   , IP_PHASE_ICE=2
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
!
!     BASIC ATMOSPHERIC PROPERTIES:
!
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             VARIABLES GIVEN IN LAYERS
!
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      REAL      !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
     &   , DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             DENSITY AT LEVELS
!
!
!     OPTICAL SWITCHES:
      LOGICAL   !, INTENT(IN)
     &     L_RESCALE
!             DELTA-RESCALING REQUIRED
!
!
!     RAYLEIGH SCATTERING:
!
      LOGICAL   !, INTENT(IN)
     &     L_RAYLEIGH
!             RAYLEIGH SCATTERING ACTIVATED
!
      REAL      !, INTENT(IN)
     &     RAYLEIGH_COEFF
!             RAYLEIGH COEFFICIENT
!
!
!     CONTINUUM PROCESSES:
      LOGICAL   !, INTENT(IN)
     &     L_CONTINUUM
!             CONTINUUM ABSORPTION ACTIVATED
!
      INTEGER   !, INTENT(IN)
     &     N_CONTINUUM
!             NUMBER OF CONTINUA
     &   , I_CONTINUUM_POINTER(NPD_CONTINUUM)
!             POINTERS TO ACTIVE CONTINUA
!
      REAL      !, INTENT(IN)
     &     K_CONTINUUM(NPD_CONTINUUM)
!             CONTINUUM EXTINCTION
     &   , AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER, NPD_CONTINUUM)
!             AMOUNTS FOR CONTINUA
!
!
!     PROPERTIES OF AEROSOLS:
!
      LOGICAL   !, INTENT(IN)
     &     L_AEROSOL
!             AEROSOLS ACTIVATED
!
      INTEGER   !, INTENT(IN)
     &     N_AEROSOL
!             NUMBER OF AEROSOL SPECIES
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATIONS OF AEROSOLS
     &   , I_HUMIDITY_POINTER(NPD_PROFILE,  NPD_LAYER)
!             POINTER TO AEROSOL LOOK-UP TABLE
!
      REAL      !, INTENT(IN)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_AEROSOL_SPECIES)
!             NUMBER DENSTY OF AEROSOLS
     &   , AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             AEROSOL ABSORPTION IN BAND/MIX RT.
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             AEROSOL SCATTERING IN BAND/MIX RT.
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             AEROSOL ASYMMETRY IN BAND
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             ARRAY OF HUMIDITIES
     &   , DELTA_HUMIDITY
!             INCREMENT IN HUMIDITY
     &   , MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)
!             MIXING RATIO OF WATER VAPOUR
!
!
!
!     PROPERTIES OF CLOUDS:
!
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD
!             CLOUDS ACTIVATED
!
!     GEOMETRY OF CLOUDS:
!
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD_LAYER
!             CLOUD VARIABLES GIVEN IN LAYERS
!
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
     &   , I_CLOUD
!             CLOUD SCHEME TO BE USED
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
     &   , N_CLOUD_PROFILE(NPD_LAYER)
!             NUMBER OF CLOUDY PROFILES IN EACH LAYER
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
!
!     MICROPHYSICAL QUANTITIES:
      INTEGER   !, INTENT(IN)
     &     N_CONDENSED
!             NUMBER OF CONDENSED COMPONENTS
     &   , I_PHASE_CMP(NPD_CLOUD_COMPONENT)
!             PHASES OF CLOUDY COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATION SCHEMES FOR CLOUDY COMPONENTS
!
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             FLAGS TO ACTIVATE CLOUDY COMPONENTS
!
      REAL      !, INTENT(IN)
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER
     &        , NPD_CLOUD_COMPONENT)
!             COEFFICIENTS IN PARAMETRIZATION SCHEMES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &        , NPD_CLOUD_COMPONENT)
!             MIXING RATIOS OF CLOUDY COMPONENTS
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER
     &      , NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CLOUDY COMPONENTS
!
!
!
!     CALCULATED OPTICAL PROPETIES:
!
      REAL      !, INTENT(OUT)
     &     K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE SCATTERING EXTINCTION
     &   , K_EXT_TOT_FREE(NPD_PROFILE, NPD_LAYER)
!             TOTAL FREE EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ASYMMETRIES
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE FORWARD SCATTERING
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SCATTERING EXTINCTION
     &   , K_EXT_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             TOTAL CLOUDY EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I_CONTINUUM
!             TEMPORARY CONTINUUM INDEX
     &   , I_POINTER
!             TEMPORARY POINTER
     &   , L
!             LOOP VARIABLE
     &   , LL
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , N_INDEX
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
!
!     TEMPORARY OPTICAL PROPERTIES:
!
      REAL
     &     K_EXT_SCAT_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)
!             SCATTERING EXTINCTION OF CLOUDY COMPONENT
     &   , K_EXT_TOT_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)
!             TOTAL EXTINCTION OF CLOUDY COMPONENT
     &   , ASYMMETRY_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRIES OF CLOUDY COMPONENT
     &   , FORWARD_SCATTER_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)
!             FORWARD SCATTERING OF CLOUDY COMPONENT
     &   , K_SCATTER(NPD_PROFILE)
!             SCATTERING VARIABLE
     &   , ASYMMETRY_PROCESS(NPD_PROFILE)
!             ASYMMETRY FACTOR FOR CURRENT PROC.
!
!
      REAL
     &     WEIGHT_UPPER
!             UPPER WEIGHT FOR INTERPOLATION
     &   , WEIGHT_LOWER
!             LOWER WEIGHT FOR INTERPOLATION
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     OPT_PROP_WATER_CLOUD, OPT_PROP_ICE_CLOUD
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
Cfpp$ NODEPCHK R
!
!
!
!     INITIALIZE THE EXTINCTION COEFFICIENTS AND THE ASYMMETRY PRODUCT.
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            K_EXT_TOT_FREE(L, I)=0.0E+00
            K_EXT_SCAT_FREE(L, I)=0.0E+00
            ASYMMETRY_FREE(L, I)=0.0E+00
         ENDDO
      ENDDO
!     FORWARD SCATTERING IS REQUIRED ONLY IN THE VISIBLE WHERE
!     DELTA-RESCALING IS PERFORMED.
      IF (L_RESCALE) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FORWARD_SCATTER_FREE(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_RAYLEIGH) THEN
!        INCLUDE RAYLEIGH SCATTERING.
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               K_EXT_SCAT_FREE(L, I)
     &             =K_EXT_SCAT_FREE(L, I)+RAYLEIGH_COEFF
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_AEROSOL) THEN
!        INCLUDE THE EFFECTS OF AEROSOL.
         DO J=1, N_AEROSOL
            IF (I_AEROSOL_PARAMETRIZATION(J)
     &         .EQ.IP_AEROSOL_PARAM_DRY) THEN
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)
     &                  +AEROSOL_MIX_RATIO(L, I, J)
     &                  *AEROSOL_ABSORPTION(1, J)
                     K_SCATTER(L)=AEROSOL_MIX_RATIO(L, I, J)
     &                  *AEROSOL_SCATTERING(1, J)
                     K_EXT_SCAT_FREE(L, I)=K_EXT_SCAT_FREE(L, I)
     &                  +K_SCATTER(L)
                     ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)
     &                  +K_SCATTER(L)*AEROSOL_ASYMMETRY(1, J)
                  ENDDO
                  IF (L_RESCALE) THEN
!                    THIS BLOCK IS PLACED WITHIN THE LOOP OVER I TO SAVE
!                    STORAGE. THE COST OF RE-EXECUTING THE TEST IS QUITE
!                    SMALL.
                     DO L=1, N_PROFILE
                        FORWARD_SCATTER_FREE(L, I)
     &                     =FORWARD_SCATTER_FREE(L, I)+K_SCATTER(L)
     &                     *(AEROSOL_ASYMMETRY(1, J))**2
                     ENDDO
                  ENDIF
               ENDDO
            ELSE IF (I_AEROSOL_PARAMETRIZATION(J)
     &         .EQ.IP_AEROSOL_PARAM_MOIST) THEN
               DO I=1, N_LAYER
! Optimizer precomputes 1/DELTA_HUMIDITY and causes zero divide, so:-
!OCL NOPREEX,NOEVAL
                  DO L=1, N_PROFILE
                     I_POINTER=I_HUMIDITY_POINTER(L, I)
                     WEIGHT_UPPER=(MEAN_REL_HUMIDITY(L, I)
     &                 -HUMIDITIES(I_POINTER, J))
     &                 /DELTA_HUMIDITY
                     WEIGHT_LOWER=1.0E+00-WEIGHT_UPPER
                     K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)
     &                  +AEROSOL_MIX_RATIO(L, I, J)
     &                  *(AEROSOL_ABSORPTION(I_POINTER, J)
     &                  *WEIGHT_LOWER+WEIGHT_UPPER
     &                  *AEROSOL_ABSORPTION(I_POINTER+1, J))
                     K_SCATTER(L)=
     &                  AEROSOL_MIX_RATIO(L, I, J)
     &                  *(AEROSOL_SCATTERING(I_POINTER, J)
     &                  *WEIGHT_LOWER+WEIGHT_UPPER
     &                  *AEROSOL_SCATTERING(I_POINTER+1, J))
                     K_EXT_SCAT_FREE(L, I)=K_EXT_SCAT_FREE(L, I)
     &                  +K_SCATTER(L)
                     ASYMMETRY_PROCESS(L)=
     &                  AEROSOL_ASYMMETRY(I_POINTER, J)
     &                  *WEIGHT_LOWER+WEIGHT_UPPER
     &                  *AEROSOL_ASYMMETRY(I_POINTER+1, J)
                     ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)
     &                  +K_SCATTER(L)*ASYMMETRY_PROCESS(L)
                  ENDDO
                  IF (L_RESCALE) THEN
                     DO L=1, N_PROFILE
                        FORWARD_SCATTER_FREE(L, I)
     &                     =FORWARD_SCATTER_FREE(L, I)+K_SCATTER(L)
     &                     *(ASYMMETRY_PROCESS(L))**2
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               WRITE(IU_ERR, '(/A, I3, A)')
     &            '*** ERROR : I_AEROSOL_PARAMETRIZATION FOR SPECIES '
     &            , J, ' HAS BEEN SET TO AN ILLEGAL VALUE.'
               IERR=I_ERR_FATAL
               RETURN

            ENDIF
         ENDDO
      ENDIF
!
      IF (L_CONTINUUM) THEN
!        INCLUDE CONTINUUM ABSORPTION.
         DO J=1, N_CONTINUUM
            I_CONTINUUM=I_CONTINUUM_POINTER(J)
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)
     &               +K_CONTINUUM(I_CONTINUUM)
     &               *AMOUNT_CONTINUUM(L, I, I_CONTINUUM)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!
!     ADD THE SCATTERING ON TO THE TOTAL EXTINCTION. THE FINAL FREE
!     ASYMMETRY IS NOT CALCULATED HERE SINCE THE PRODUCT OF ASYMMETRY
!     AND SCATTERING IS ALSO NEEDED TO CALCULATE THE CLOUDY ASYMMETRY.
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)
     &         +K_EXT_SCAT_FREE(L, I)
         ENDDO
      ENDDO
!
!
!     IF THERE ARE NO CLOUDS CALCULATE THE FINAL OPTICAL PROPERTIES
!     AND RETURN TO THE CALLING ROUTINE.
!
      IF (.NOT.L_CLOUD) THEN
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               IF (K_EXT_SCAT_FREE(L, I).GT.TOL_DIV) THEN
                  ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)
     &               /K_EXT_SCAT_FREE(L, I)
               ENDIF
            ENDDO
         ENDDO
         IF (L_RESCALE) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (K_EXT_SCAT_FREE(L, I).GT.TOL_DIV) THEN
                     FORWARD_SCATTER_FREE(L, I)
     &                  =FORWARD_SCATTER_FREE(L, I)
     &                  /K_EXT_SCAT_FREE(L, I)
                  ENDIF
               ENDDO
            ENDDO
!
         ENDIF
!
         RETURN
!
      ENDIF
!
!
!
!
!     ADDITION OF CLOUDY PROPERTIES:
!
!
!     ADD IN BACKGROUND CONTIBUTIONS:
!
!
!     ALL THE PROCESSES OCCURRING OUTSIDE CLOUDS ALSO OCCUR
!     WITHIN THEM.
      DO K=1, N_CLOUD_TYPE
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               K_EXT_TOT_CLOUD(L, I, K)=K_EXT_TOT_FREE(L, I)
               K_EXT_SCAT_CLOUD(L, I, K)=K_EXT_SCAT_FREE(L, I)
               ASYMMETRY_CLOUD(L, I, K)=ASYMMETRY_FREE(L, I)
               FORWARD_SCATTER_CLOUD(L, I, K)
     &            =FORWARD_SCATTER_FREE(L, I)
            ENDDO
         ENDDO
      ENDDO
!
!
!
!     ADD ON THE TERMS REPRESENTING PROCESSES WITHIN CLOUDS.
!
!     LOOP OVER THE CONDENSED COMPONENTS, CALCULATING THEIR OPTICAL
!     PROPERTIES AND THEN ASSIGN THEM TO THE ARRAYS FOR THE TYPES OF
!     CLOUD.
!
      DO K=1, N_CONDENSED
!
!        FLAGS FOR DEALING WITH COMPONENTS WERE SET IN THE SUBROUTINE
!        SET_CLOUD_POINTER. WE NOW DETERMINE WHETHER THE COMPONENT IS
!        TO BE INCLUDED AND CALCULATE ITS OPTICAL PROPERTIES ACCORDING
!        TO THE PHASE OF THE COMPONENT. THESE CONTRIBUTIONS ARE ADDED
!        TO THE ARRAYS FOR THE SELECTED TYPE OF CLOUD.
!
         IF (L_CLOUD_CMP(K)) THEN
!
            IF (I_PHASE_CMP(K).EQ.IP_PHASE_WATER) THEN
!
!              INCLUDE SCATTERING BY WATER DROPLETS.
!
               CALL OPT_PROP_WATER_CLOUD(IERR
     &            , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &            , N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , L_RESCALE, L_LAYER, L_CLOUD_LAYER
     &            , I_CONDENSED_PARAM(K), CONDENSED_PARAM_LIST(1, K)
     &            , CONDENSED_MIX_RATIO(1, 0, K)
     &            , CONDENSED_DIM_CHAR(1, 0, K)
     &            , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP
     &            , ASYMMETRY_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP
     &            , NPD_PROFILE, NPD_LAYER
     &            , NPD_CLOUD_PARAMETER
     &            )
!
            ELSE IF (I_PHASE_CMP(K).EQ.IP_PHASE_ICE) THEN
!
!              INCLUDE SCATTERING BY ICE CRYSTALS.
!
               CALL OPT_PROP_ICE_CLOUD(IERR
     &            , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &            , N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &            , L_RESCALE, L_LAYER, L_CLOUD_LAYER
     &            , I_CONDENSED_PARAM(K), CONDENSED_PARAM_LIST(1, K)
     &            , CONDENSED_MIX_RATIO(1, 0, K)
     &            , CONDENSED_DIM_CHAR(1, 0, K)
     &            , T, DENSITY
     &            , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP
     &            , ASYMMETRY_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP
     &            , NPD_PROFILE, NPD_LAYER
     &            , NPD_CLOUD_PARAMETER
     &            )
!
            ENDIF
!
!
!
!           INCREMENT THE ARRAYS OF OPTICAL PROPERTIES.
!
!
            DO I=N_CLOUD_TOP, N_LAYER
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  K_EXT_TOT_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               =K_EXT_TOT_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               +K_EXT_TOT_CLOUD_COMP(L, I)
                  K_EXT_SCAT_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               =K_EXT_SCAT_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               +K_EXT_SCAT_CLOUD_COMP(L, I)
                  ASYMMETRY_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               =ASYMMETRY_CLOUD(L, I, I_CLOUD_TYPE(K))
     &               +ASYMMETRY_CLOUD_COMP(L, I)
               ENDDO
            ENDDO
            IF (L_RESCALE) THEN
               DO I=N_CLOUD_TOP, N_LAYER
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     FORWARD_SCATTER_CLOUD(L, I, I_CLOUD_TYPE(K))
     &                  =FORWARD_SCATTER_CLOUD(L, I, I_CLOUD_TYPE(K))
     &                  +FORWARD_SCATTER_CLOUD_COMP(L, I)
                  ENDDO
               ENDDO
            ENDIF
!
         ENDIF
!
      ENDDO
!
!
!
!
!     CALCULATE THE FINAL OPTICAL PROPERTIES.
!     THE SCATTERING WAS INCLUDED IN THE FREE TOTAL EXTINCTION EARLIER,
!     BUT WE HAVE YET TO DIVIDE THE PRODUCT OF THE ASYMMETRY AND THE
!     SCATTERING BY THE MEAN SCATTERING.
!
      DO I=1, N_LAYER
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (K_EXT_SCAT_FREE(L,I).GT.TOL_DIV) THEN
             N_INDEX       =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
               ASYMMETRY_FREE(INDEX(K), I)=ASYMMETRY_FREE(INDEX(K), I)
     &            /K_EXT_SCAT_FREE(INDEX(K), I)
         ENDDO
!
         IF (L_RESCALE) THEN
            DO K=1, N_INDEX
               FORWARD_SCATTER_FREE(INDEX(K), I)
     &            =FORWARD_SCATTER_FREE(INDEX(K), I)
     &            /K_EXT_SCAT_FREE(INDEX(K), I)
            ENDDO
         ENDIF
      ENDDO
!
!
!     REPEAT FOR CLOUDS.
      DO K=1, N_CLOUD_TYPE
         DO I=N_CLOUD_TOP, N_LAYER
!
            J      =1
            N_INDEX=0
            DO L   =1,N_PROFILE
              IF (K_EXT_SCAT_CLOUD(L,I,K).GT.TOL_DIV) THEN
                INDEX(J)=L
                J       =J+1
                N_INDEX =N_INDEX+1
              END IF
            END DO

            DO J=1, N_INDEX
               ASYMMETRY_CLOUD(INDEX(J), I, K)
     &            =ASYMMETRY_CLOUD(INDEX(J), I, K)
     &            /K_EXT_SCAT_CLOUD(INDEX(J), I, K)
            ENDDO
            IF (L_RESCALE) THEN
               DO J=1, N_INDEX
                  FORWARD_SCATTER_CLOUD(INDEX(J), I, K)
     &               =FORWARD_SCATTER_CLOUD(INDEX(J), I, K)
     &               /K_EXT_SCAT_CLOUD(INDEX(J), I, K)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!
!
!
      RETURN
      END
