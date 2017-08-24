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
!+ Subroutine to calculate cloudy fluxes by division into columns.
!
! Method:
!       A number of atmospheric profiles are taken and split into
!       columns in which each layer is homogeneous. The areal
!       coverages of these columns has been calculated before. The
!       sub-columns are passed into a long vector to be passed to a
!       multicolumn solver.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New solvers added.
!                                               (J. M. Edwards)
!       4.4             19-09-97                Addressing for long
!                                               rows corrected and
!                                               missing initialization
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete solvers
!                                               removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CLOUD_COLUMN(IERR
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &   , I_2STREAM
!                       Corrections to Two-stream Equations
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &   , I_SOLVER, N_AUGMENT
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
     &   , N_COLUMN, L_COLUMN, AREA_COLUMN
!                       Cloudy Optical Properties
     &   , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL
!                       Flags for Clear-sky Calculations
     &   , L_CLEAR, I_SOLVER_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
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
!     INCLUDE COMDECKS
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             SPECTRAL REGION
      INTEGER   !, INTENT(IN)
     &     I_2STREAM
!             TWO STREAM SCHEME
     &   , I_SOLVER
!             SOLVER FOR TWO-STREAM EQUATIONS
     &   , N_AUGMENT
!             LENGTH OF LONG VECTOR
     &   , I_SOLVER_CLEAR
!             SOLVER FOR CLEAR FLUXES
      LOGICAL   !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE CLEAR-SKY FLUXES
     &   , L_SCALE_SOLAR
!             SCALE SOLAR BEAM
     &   , L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE TERM
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION TO 2-STREAM
!
!     FIELDS FOR EQUIVALENT EXTINCTION
      REAL  !, INTENT(IN)
     &     ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
!
!     CLEAR-SKY OPTICAL PROPETIES
      REAL      !, INTENT(IN)
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE FRACTIONAL FORWARD SCATTER
!
!     CLOUDY OPTICAL PROPETIES
      REAL      !, INTENT(IN)
     &     TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             OPTICAL DEPTH IN CLOUD
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FRACTIONAL FORWARD SCATTER
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             ALBEDO OF SINGLE SCATTERING IN CLOUD
!
!     PLANCKIAN TERMS:
      REAL      !, INTENT(IN)
     &     DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             CHANGE IN PLANCK FUNCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             TWICE 2ND DIFFERENCES IN PLANCKIAN
     &   , PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE FUNCTION
!
!     CONDITIONS AT TOA
      REAL      !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
!
!     CONDITIONS AT SURFACE
      REAL      !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO OF GROUND
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO OF GROUND
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FUNCTION OF GROUND
!
!     CLOUD GEOMETRY
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP
!             TOP CLOUDY LAYER
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
     &   , N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS
      LOGICAL   !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             TYPE FLAG FOR EACH LAYER/COLUMN
      REAL      !, INTENT(IN)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_COLUMN, NPD_CLOUD_TYPE)
!             AREA OF EACH COLUMN
     &   , AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREA OF EACH COLUMN
!
      REAL      !, INTENT(IN)
     &     GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
!     FLUXES CALCULATED:
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     N_SOURCE_COEFF
!             NUMBER OF SOURCE COEFFICIENTS
     &   , N_EQUATION
!             NUMBER OF EQUATIONS SOLVED
     &   , I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , JS
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
      INTEGER
     &     N_LONG
!             LENGTH OF LONG VECTOR
     &   , N_PROFILE_SOLVED
!             NUMBER OF PROFILES SOLVED AT ONCE
     &   , N_COLUMN_DONE
!             NUMBER OF COLUMNS ALREADY ASSIGNED
     &   , OFFSET
!             OFFSET IN LIST OF PROFILES
     &   , I_PROFILE
!             PROFILE BEING CONSIDERED
!
      REAL
     &     SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             FREE SOURCE COEFFICIENTS
     &   , TRANS_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DIFFUSE TRANSMISSION
     &   , TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DIRECT TRANSMISSION
     &   , REFLECT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE REFLECTANCE
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DOWNWARD SOURCE
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE UPWARD SOURCE
!
      REAL
     &     SOURCE_COEFF_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             CLOUDY TWO-STREAM SOURCE COEFFICIENTS
     &   , TRANS_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DIFFUSE TRANSMISSION
     &   , REFLECT_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY REFLECTANCE
     &   , TRANS_0_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DIRECT TRANSMISSION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DOWNWARD SOURCE
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY UPWARD SOURCE
!
      REAL
     &     FLUX_LONG(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUXES IN COLUMNS
     &   , FLUX_DIRECT_LONG(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUXES IN COLUMNS
!
      REAL
     &     SOURCE_GROUND_LONG(NPD_PROFILE)
!             SOURCE FUNCTION OF GROUND IN COL.
     &   , FLUX_INC_DIRECT_LONG(NPD_PROFILE)
!             DIRECT FLUX INCIDENT ON COLUMN
     &   , FLUX_INC_DOWN_LONG(NPD_PROFILE)
!             DIRECT FLUX INCIDENT ON COLUMN
     &   , ALBEDO_SURFACE_DIFF_LONG(NPD_PROFILE)
!             DIFFUSE ALBEDO OF GROUND
     &   , ALBEDO_SURFACE_DIR_LONG(NPD_PROFILE)
!             DIRECT ALBEDO OF GROUND
     &   , TRANS_LONG(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION IN LONG ARRAY
     &   , REFLECT_LONG(NPD_PROFILE, NPD_LAYER)
!             REFLECTION IN LONG ARRAY
     &   , TRANS_0_LONG(NPD_PROFILE, NPD_LAYER)
!             SOLAR COEFFICIENT IN LONG ARRAY
     &   , SOURCE_COEFF_LONG(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS IN LONG ARRAY
     &   , S_UP_LONG(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
     &   , S_DOWN_LONG(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE FUNCTION
     &   , SCALE_SOLAR_LONG(NPD_PROFILE, NPD_LAYER)
!             LONG VECTOR OF SOLAR SCALINGS
     &   , WORK_1(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORK ARRAY
     &   , WORK_2(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORK ARRAY
!
      REAL
     &     A3(NPD_PROFILE, 3, 2*NPD_LAYER+2)
!             TRIDIAGONAL MATRIX
     &   , A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             RIGHT-HAND SIDES OF EQUATIONS
!
!
!     FUNCTIONS CALLED:
      INTEGER
     &     SET_N_SOURCE_COEFF
!             FUNCTION TO SET NUMBER OF SOURCE COEFFICIENTS
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     TWO_COEFF, TWO_COEFF_CLOUD, IR_SOURCE, SOLAR_SOURCE
     &  , SET_MATRIX_NET, TRIDIAG_SOLVER_UP, SET_MATRIX_FULL      
     &   , SET_MATRIX_PENTADIAGONAL, BAND_SOLVER
     &   , SOLVER_HOMOGEN_DIRECT, CLEAR_SUPPLEMENT
!
!
!
!
!     ENTER A SUMMING LOOP TO CALCULATE THE TOTAL FLUX BY ADDING UP
!     THE FLOW OF ENERGY IN EACH COLUMN.
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
      DO I=1, N_AUGMENT
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=0.0E+00
         ENDDO
      ENDDO
!
!     SET THE NUMBER OF SOURCE COEFFICIENTS FOR THE APPROXIMATION
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
!     THE FUNDAMENTAL PARAMETERS OF THE TWO-STREAM EQUATIONS CAN BE
!     PRECALCULATED.
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
!     THE INFRA-RED SOURCE FUNCTIONS DEPEND ONLY ON THE LAYER IN WHICH
!     THEY ARE EVALUATED AND, UNLIKE THE VISIBLE SOURCE FUNCTIONS, THEY
!     MAY BE PRECALCULATED.
      IF (ISOLIR.EQ.IP_INFRA_RED) THEN
         CALL IR_SOURCE(N_PROFILE, 1, N_LAYER
     &      , SOURCE_COEFF_FREE, DIFF_PLANCK
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &      , GROUND_EMISSION, N_LAYER
     &      , TAU_FREE, TRANS_FREE
     &      , S_DOWN_FREE, S_UP_FREE
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         CALL IR_SOURCE(N_PROFILE, N_CLOUD_TOP, N_LAYER
     &      , SOURCE_COEFF_CLOUD, DIFF_PLANCK
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &      , GROUND_EMISSION, N_LAYER
     &      , TAU_CLOUD, TRANS_CLOUD
     &      , S_DOWN_CLOUD, S_UP_CLOUD
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!     THE MAIN LOOP: PROFILES ARE ADDED ON TO THE LONG ARRAY UNTIL
!     IT IS NO LONGER POSSIBLE TO SOLVE FOR THEM ALL IN ONE GO.
      OFFSET=0
      N_LONG=0
      I_PROFILE=1
!
10       IF (N_LONG+N_COLUMN(I_PROFILE).LE.NPD_PROFILE) THEN
!           CONTINUE FEEDING PROFILES INTO THE LONG ARRAY:
!
            DO J=1, N_COLUMN(I_PROFILE)
!
!              ASSIGN THE OPTICAL PROPERTIES TO EACH LAYER WITHIN THE
!              COLUMN. J IS THE INDEX OF COLUMNS AT A GRID-POINT:
!              K INDEXES POINTS IN THE LONG VECTOR.
               K=N_LONG+J
!
               DO I=1, N_LAYER
                  IF (L_COLUMN(I_PROFILE, I, J)) THEN
!                    THE LAYER IS CLOUDY.
                     TRANS_LONG(K, I)=TRANS_CLOUD(I_PROFILE, I)
                     REFLECT_LONG(K, I)=REFLECT_CLOUD(I_PROFILE, I)
                  ELSE
!                    THE LAYER IS FREE OF CLOUD.
                     TRANS_LONG(K, I)=TRANS_FREE(I_PROFILE, I)
                     REFLECT_LONG(K, I)=REFLECT_FREE(I_PROFILE, I)
                  ENDIF
               ENDDO
               IF (ISOLIR.EQ.IP_SOLAR) THEN
                  FLUX_INC_DIRECT_LONG(K)=FLUX_INC_DIRECT(I_PROFILE)
                  SOURCE_GROUND_LONG(K)=0.0E+00
                  DO I=1, N_LAYER
                     IF (L_COLUMN(I_PROFILE, I, J)) THEN
                        TRANS_0_LONG(K, I)
     &                     =TRANS_0_CLOUD(I_PROFILE, I)
                        DO JS=1, N_SOURCE_COEFF
                           SOURCE_COEFF_LONG(K, I, JS)
     &                        =SOURCE_COEFF_CLOUD(I_PROFILE, I, JS)
                        ENDDO
                     ELSE
                        TRANS_0_LONG(K, I)=TRANS_0_FREE(I_PROFILE, I)
                        DO JS=1, N_SOURCE_COEFF
                           SOURCE_COEFF_LONG(K, I, JS)
     &                        =SOURCE_COEFF_FREE(I_PROFILE, I, JS)
                        ENDDO
                     ENDIF
                  ENDDO
                  IF (L_SCALE_SOLAR) THEN
                     DO I=1, N_LAYER
                        SCALE_SOLAR_LONG(K, I)
     &                     =ADJUST_SOLAR_KE(I_PROFILE, I)
                     ENDDO
                  ENDIF
               ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
                  DO I=1, N_LAYER
                     IF (L_COLUMN(I_PROFILE, I, J)) THEN
                        S_UP_LONG(K, I)=S_UP_CLOUD(I_PROFILE, I)
                        S_DOWN_LONG(K, I)=S_DOWN_CLOUD(I_PROFILE, I)
                     ELSE
                        S_DOWN_LONG(K, I)=S_DOWN_FREE(I_PROFILE, I)
                        S_UP_LONG(K, I)=S_UP_FREE(I_PROFILE, I)
                     ENDIF
                  ENDDO
                  SOURCE_GROUND_LONG(K)=SOURCE_GROUND(I_PROFILE)
               ENDIF
               FLUX_INC_DOWN_LONG(K)=FLUX_INC_DOWN(I_PROFILE)
               ALBEDO_SURFACE_DIFF_LONG(K)
     &            =ALBEDO_SURFACE_DIFF(I_PROFILE)
               ALBEDO_SURFACE_DIR_LONG(K)
     &            =ALBEDO_SURFACE_DIR(I_PROFILE)
            ENDDO
!
            N_LONG=N_LONG+N_COLUMN(I_PROFILE)
!
!           INCREMENT I_PROFILE AND RETURN TO SEE IF ANOTHER PROFILE
!           MAY BE ADDED.
            IF (I_PROFILE.LT.N_PROFILE) THEN
               I_PROFILE=I_PROFILE+1
               GOTO 10
            ENDIF
!
         ELSE IF (N_COLUMN(I_PROFILE).GT.NPD_PROFILE) THEN
            WRITE(IU_ERR, '(/A, I5, A, /A)')
     &         '*** ERROR: PROFILE ', I_PROFILE
     &         , ' CONTAINS TOO MANY COLUMNS.'
     &         , 'INCREASE NPD_PROFILE AND RECOMPILE.'
         ELSE
!           IF NO MORE PROFILES CAN BE ADDED THE EQUATIONS ARE SOLVED.
!           RESET THE POINTER TO POINT TO THE LAST PROFILE CONSIDERED.
            I_PROFILE=I_PROFILE-1
         ENDIF
!
         N_PROFILE_SOLVED=I_PROFILE-OFFSET
!
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            CALL SOLAR_SOURCE(N_LONG, N_LAYER
     &         , FLUX_INC_DIRECT_LONG
     &         , TRANS_0_LONG, SOURCE_COEFF_LONG
     &         , L_SCALE_SOLAR, SCALE_SOLAR_LONG
     &         , FLUX_DIRECT_LONG
     &         , S_DOWN_LONG, S_UP_LONG
     &         , NPD_PROFILE, NPD_LAYER
     &         )
         ELSE
!           SET THE DIRECT FLUX AT THE GROUND FOR USE IN THE SOLVER.
            DO K=1, N_LONG
               FLUX_DIRECT_LONG(K, N_LAYER)=0.0E+00
            ENDDO
         ENDIF
!
!
!        SELECT AN APPROPRIATE SOLVER FOR THE EQUATIONS.
!
         IF (I_SOLVER.EQ.IP_SOLVER_PENTADIAGONAL) THEN
!
            CALL SET_MATRIX_PENTADIAGONAL(N_LONG, N_LAYER
     &         , TRANS_LONG, REFLECT_LONG
     &         , S_DOWN_LONG, S_UP_LONG
     &         , ALBEDO_SURFACE_DIFF_LONG, ALBEDO_SURFACE_DIR_LONG
     &         , FLUX_DIRECT_LONG(1, N_LAYER), FLUX_INC_DOWN_LONG
     &         , SOURCE_GROUND_LONG
     &         , A5, B
     &         , NPD_PROFILE, NPD_LAYER
     &         )
            N_EQUATION=2*N_LAYER+2
!
            CALL BAND_SOLVER(N_LONG, N_EQUATION
     &         , 2, 2
     &         , A5, B
     &         , FLUX_LONG
     &         , NPD_PROFILE, 2*NPD_LAYER+2
     &         , WORK_1
     &         )
!
         ELSE IF (I_SOLVER.EQ.IP_SOLVER_HOMOGEN_DIRECT) THEN
!
            CALL SOLVER_HOMOGEN_DIRECT(N_LONG, N_LAYER
     &         , TRANS_LONG, REFLECT_LONG
     &         , S_DOWN_LONG, S_UP_LONG
     &         , ALBEDO_SURFACE_DIFF_LONG, ALBEDO_SURFACE_DIR_LONG
     &         , FLUX_DIRECT_LONG(1, N_LAYER), FLUX_INC_DOWN_LONG
     &         , SOURCE_GROUND_LONG
     &         , FLUX_LONG
     &         , NPD_PROFILE, NPD_LAYER
     &         )
!
         ELSE
!
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: THE SOLVER AND CLOUD SCHEME '
     &         //'ARE NOT COMPATIBLE.'
            IERR=I_ERR_FATAL
            RETURN
!
         ENDIF
!
!
!        ADD THE PARTIAL FLUX FOR THE COLUMN ONTO THE CUMULATIVE TOTAL.
!
         N_COLUMN_DONE=0
         DO L=OFFSET+1, N_PROFILE_SOLVED+OFFSET
            DO J=1, N_COLUMN(L)
               K=J+N_COLUMN_DONE
               IF (ISOLIR.EQ.IP_SOLAR) THEN
                  DO I=0, N_LAYER
                     FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)
     &                  +FLUX_DIRECT_LONG(K, I)*AREA_COLUMN(L, J)
                  ENDDO
               ENDIF
               DO I=1, N_AUGMENT
                  FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)
     &               +FLUX_LONG(K, I)*AREA_COLUMN(L, J)
               ENDDO
            ENDDO
            N_COLUMN_DONE=N_COLUMN_DONE+N_COLUMN(L)
         ENDDO
!
!        THE OFFSET IS NOW ADVANCED FOR THE NEXT LOOP AND A NEW
!        LONG VECTOR IS FORMED UNLESS WE HAVE ALREADY SOLVED FOR ALL
!        THE PROFILES.
         OFFSET=OFFSET+N_PROFILE_SOLVED
         IF (OFFSET.LT.N_PROFILE) THEN
            N_LONG=0
            I_PROFILE=I_PROFILE+1
            GOTO 10
         ENDIF
!
!
!
      IF (L_CLEAR) THEN
         CALL CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER, I_SOLVER_CLEAR
     &      , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE
     &      , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN
     &      , S_DOWN_FREE, S_UP_FREE
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , SOURCE_GROUND
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &      , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
      RETURN
      END
