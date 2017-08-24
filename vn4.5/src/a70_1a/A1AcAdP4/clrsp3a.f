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
!+ Subroutine to calculate clear-sky fluxes.
!
! Method:
!       This subroutine is called after fluxes including clouds have
!       been calculated to find the corresponding clear-sky fluxes.
!       The optical properties of the column are already known.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-04-96                New solver added
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete solvers
!                                               removed.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER
     &   , I_SOLVER_CLEAR
     &   , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE
     &   , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN
     &   , S_DOWN_FREE, S_UP_FREE
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , SOURCE_GROUND
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &   , NPD_PROFILE, NPD_LAYER
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
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , ISOLIR
!             SPECTRAL REGION
     &   , I_SOLVER_CLEAR
!             SOLVER FOR CLEAR FLUXES
      LOGICAL   !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR BEAM
      REAL  !, INTENT(IN)
     &     TRANS_FREE(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION COEFFICIENTS
     &   , REFLECT_FREE(NPD_PROFILE, NPD_LAYER)
!             REFLECTION COEFFICIENTS
     &   , TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)
!             DIRECT TRANSMISSION COEFFICIENTS
     &   , SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER)
!             COEFFICIENTS IN SOURCE TERMS
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             SCALING OF SOLAR BEAM
!
!
      REAL  !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUXES
!
!
!     DUMMY VARIABALES.
      INTEGER
     &     N_EQUATION
!             NUMBER OF EQUATIONS
      REAL
     &     A3(NPD_PROFILE, 3, 2*NPD_LAYER+2)
!             TRIDIAGONAL MATRIX
     &   , A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             RHS OF MATRIX EQUATION
     &   , WORK_1(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORKING ARRAY FOR SOLVER
     &   , WORK_2(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORKING ARRAY FOR SOLVER
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &  SOLAR_SOURCE, SET_MATRIX_NET, TRIDIAG_SOLVER_UP           
     &   , SET_MATRIX_FULL, SET_MATRIX_PENTADIAGONAL
     &   , BAND_SOLVER, SOLVER_HOMOGEN_DIRECT
!
!
!     THE SOURCE FUNCTIONS ONLY NEED TO BE RECALCULATED IN THE VISIBLE.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         CALL SOLAR_SOURCE(N_PROFILE, N_LAYER
     &      , FLUX_INC_DIRECT
     &      , TRANS_0_FREE, SOURCE_COEFF_FREE
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &      , FLUX_DIRECT_CLEAR
     &      , S_DOWN_FREE, S_UP_FREE
     &      , NPD_PROFILE, NPD_LAYER
     &      )
      ENDIF
!
!
!     SELECT AN APPROPRIATE SOLVER FOR THE EQUATIONS OF TRANSFER.
!
      IF (I_SOLVER_CLEAR.EQ.IP_SOLVER_PENTADIAGONAL) THEN
!
!        CALCULATE THE ELEMENTS OF THE MATRIX EQUATIONS.
         CALL SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER
     &      , TRANS_FREE, REFLECT_FREE
     &      , S_DOWN_FREE, S_UP_FREE
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN
     &      , SOURCE_GROUND
     &      , A5, B
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         N_EQUATION=2*N_LAYER+2
!
         CALL BAND_SOLVER(N_PROFILE, N_EQUATION
     &      , 2, 2
     &      , A5, B
     &      , FLUX_TOTAL_CLEAR
     &      , NPD_PROFILE, 2*NPD_LAYER+2
     &      , WORK_1
     &      )
!
      ELSE IF (I_SOLVER_CLEAR.EQ.IP_SOLVER_HOMOGEN_DIRECT) THEN
!
!        SOLVE FOR THE FLUXES IN THE COLUMN DIRECTLY.
         CALL SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER
     &      , TRANS_FREE, REFLECT_FREE
     &      , S_DOWN_FREE, S_UP_FREE
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &      , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN
     &      , SOURCE_GROUND
     &      , FLUX_TOTAL_CLEAR
     &      , NPD_PROFILE, NPD_LAYER
     &      )
!
      ELSE
!
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: THE SOLVER SPECIFIED IS NOT VALID '
     &      //'FOR CLEAR FLUXES.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      RETURN
      END
