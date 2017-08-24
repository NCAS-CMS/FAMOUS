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
!+ Subroutine to solve the two-stream equations in a column.
!
! Method:
!       The coefficients of the two-stream equations are calculated.
!       From these we obtain the transmission and reflection
!       coefficients and the source terms. Depending on the solver
!       selected, an appropriate set of matrix equations is formulated
!       and solved to give the fluxes.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-04-96                New solver added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete solvers
!                                               removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_STREAM(IERR
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER
!                       Two-stream Scheme
     &   , I_2STREAM
!                       Corrections to Two-stream scheme
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION
!                       Options for Solver
     &   , L_NET, I_SOLVER
!                       Options for Equivalent Extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
!                       Spectral Region
     &   , ISOLIR
!                       Infra-red Properties
     &   , DIFF_PLANCK
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
!                       Conditions at TOA
     &   , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0
!                       Surface Conditions
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND
!                       Single Scattering Properties
     &   , TAU, OMEGA, ASYMMETRY
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL
!                       Flag for Clear-sky Fluxes
     &   , L_CLEAR
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
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
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , ISOLIR
!             SPECTRAL REGION
     &   , I_SOLVER
!             SOLVER EMPLOYED
     &   , I_2STREAM
!             TWO-STREAM SCHEME
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES
     &   , L_CLEAR
!             CALCULATE CLEAR FLUXES
     &   , L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR FLUX
     &   , L_IR_SOURCE_QUAD
!             USE QUADRATIC SOURCE TERM
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION FOR 2-STREAM
      REAL      !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRY
     &   , SEC_0(NPD_PROFILE)
!             SECANTS OF SOLAR ZENITH ANGLES
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE IN PI*PLANCK FUNCTION
     &   , SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCES OF PLANCKIAN
     &   , PLANCK_SOURCE(NPD_PROFILE, NPD_LAYER)
!             PLANCKIAN SOURCE FUNCTION
     &   , GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUXES
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUXES
!
!
!     LOCAL VARIABALES.
      INTEGER
     &     N_EQUATION
!             NUMBER OF EQUATIONS
     &   , I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION OF LAYER
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)
!             REFLECTANCE OF LAYER
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)
!             DIRECT TRANSMITTANCE
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE
      REAL
     &     A3(NPD_PROFILE, 3, 2*NPD_LAYER+2)
!             TRIDIAGONAL MATRIX
     &   , A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)
!             PENTADIGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             RHS OF MATRIX EQUATION
     &   , WORK_1(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORKING ARRAY FOR SOLVER
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     TWO_COEFF, SOLAR_SOURCE, IR_SOURCE
     &   , SET_MATRIX_PENTADIAGONAL
     &   , BAND_SOLVER, SOLVER_HOMOGEN_DIRECT
!
!
!
!     CALCULATE THE TWO-STREAM COEFFICIENTS.
      CALL TWO_COEFF(IERR
     &   , N_PROFILE, 1, N_LAYER
     &   , I_2STREAM, L_IR_SOURCE_QUAD
     &   , ASYMMETRY, OMEGA, TAU
     &   , ISOLIR, SEC_0
     &   , TRANS, REFLECT, TRANS_0
     &   , SOURCE_COEFF
     &   , NPD_PROFILE, NPD_LAYER
     &   )
      IF (IERR.NE.I_NORMAL) THEN
         RETURN
      ENDIF
!
!     CALCULATE THE APPROPRIATE SOURCE TERMS.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         CALL SOLAR_SOURCE(N_PROFILE, N_LAYER
     &      , FLUX_INC_DIRECT
     &      , TRANS_0, SOURCE_COEFF
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &      , FLUX_DIRECT
     &      , S_DOWN, S_UP
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
         CALL IR_SOURCE(N_PROFILE, 1, N_LAYER
     &      , SOURCE_COEFF, DIFF_PLANCK
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &      , GROUND_EMISSION, N_LAYER
     &      , TAU, TRANS
     &      , S_DOWN, S_UP
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!     SELECT AN APPROPRIATE SOLVER FOR THE EQUATIONS OF TRANSFER.
!
      IF (I_SOLVER.EQ.IP_SOLVER_PENTADIAGONAL) THEN
         CALL SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER
     &      , TRANS, REFLECT
     &      , S_DOWN, S_UP
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , FLUX_DIRECT(1, N_LAYER), FLUX_INC_DOWN
     &      , SOURCE_GROUND
     &      , A5, B
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         N_EQUATION=2*N_LAYER+2
!
         CALL BAND_SOLVER(N_PROFILE, N_EQUATION
     &      , 2, 2
     &      , A5, B
     &      , FLUX_TOTAL
     &      , NPD_PROFILE, 2*NPD_LAYER+2
     &      , WORK_1
     &      )
!
      ELSE IF (I_SOLVER.EQ.IP_SOLVER_HOMOGEN_DIRECT) THEN
!
         CALL SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER
     &      , TRANS, REFLECT
     &      , S_DOWN, S_UP
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , FLUX_DIRECT(1, N_LAYER), FLUX_INC_DOWN
     &      , SOURCE_GROUND
     &      , FLUX_TOTAL
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE
!
         WRITE(IU_ERR, '(/A)')
     &      '***ERROR: THE SOLVER AND CLOUD SCHEME ARE NOT COMPATIBLE.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      IF (L_CLEAR) THEN
!        THE CLEAR FLUXES HERE CAN BE COPIED DIRECTLY WITHOUT
!        ANY FURTHER CALCULATION.
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR(L, I)=FLUX_DIRECT(L, I)
               ENDDO
            ENDDO
         ENDIF
         IF (L_NET) THEN
            DO I=1, N_LAYER+1
               DO L=1, N_PROFILE
                  FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL(L, I)
               ENDDO
            ENDDO
         ELSE
            DO I=1, 2*N_LAYER+2
               DO L=1, N_PROFILE
                  FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL(L, I)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!
!
      RETURN
      END
