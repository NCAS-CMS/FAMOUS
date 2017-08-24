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
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!       The basic two-stream coefficients in the differential equations
!       are calculated. These are then used to determine the
!       transmission and reflection coefficients. Coefficients for
!       determining the solar or infra-red source terms are calculated.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!LL  4.5  27/04/98  Add Fujitsu vectorization directive.
!LL                                           RBarnes@ecmwf.int
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
! Fujitsu directive to encourage vectorization for whole routine
!OCL NOVREC
      SUBROUTINE TWO_COEFF(IERR
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , I_2STREAM, L_IR_SOURCE_QUAD
     &   , ASYMMETRY, OMEGA, TAU
     &   , ISOLIR, SEC_0
     &   , TRANS, REFLECT, TRANS_0
     &   , SOURCE_COEFF
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
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST
!             LAST LAYER TO CONSIDER
     &   , ISOLIR
!             SPECTRAL REGION
     &   , I_2STREAM
!             TWO STREAM SCHEME
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!               USE A QUADRATIC SOURCE FUNCTION
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL      !, INTENT(IN)
     &     ASYMMETRY(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRY FACTOR
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
     &   , TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH
!
!     SOLAR BEAM
      REAL      !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL      !, INTENT(OUT)
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , N_INDEX
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL
     &     LAMBDA(NPD_PROFILE, NPD_LAYER)
!             COEFFICIENTS IN TWO-STREAM EQUATIONS
     &   , SUM(NPD_PROFILE, NPD_LAYER)
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
     &   , GAMMA_UP(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR DOWNWARD RADIATION
!
      REAL
     &     TARGET
!             TARGET TO SEARCH FOR
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     TWO_COEFF_BASIC, SOLAR_COEFFICIENT_BASIC 
     &   , TRANS_SOURCE_COEFF
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
Cfpp$ NODEPCHK R
!
!
!
!     PERTURB THE SINGLE SCATTERING ALBEDO AWAY FROM 1 TO AVOID
!     LATER DIVISION BY 0.
      TARGET=1.0E+00-TOL_DIV
      DO I=I_LAYER_FIRST, I_LAYER_LAST
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (OMEGA(L,I).GT.TARGET) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            OMEGA(INDEX(K), I)=TARGET
         ENDDO
      ENDDO
!
!     CALCULATE THE BASIC TWO-STREAM COEFFICIENTS.
      CALL TWO_COEFF_BASIC(IERR
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , I_2STREAM
     &   , ASYMMETRY, OMEGA
     &   , SUM, DIFF
     &   , NPD_PROFILE, NPD_LAYER
     &   )
      IF (IERR.NE.I_NORMAL) THEN
         RETURN
      ENDIF
!
!     LAMBDA IS NOW CALCULATED.
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            LAMBDA(L, I)=SQRT(SUM(L, I)*DIFF(L, I))
         ENDDO
      ENDDO
!
!
!     CALCULATE THE BASIC COEFFICIENTS FOR THE SOLAR SOURCE TERMS.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
!        LAMBDA MAY BE PERTURBED BY THIS ROUTINE TO AVOID
!        ILL-CONDITIONING FOR THE SINGULAR ZENITH ANGLE.
         CALL SOLAR_COEFFICIENT_BASIC(IERR
     &      , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &      , OMEGA, ASYMMETRY, SEC_0
     &      , I_2STREAM
     &      , SUM, DIFF, LAMBDA
     &      , GAMMA_UP, GAMMA_DOWN
     &      , NPD_PROFILE, NPD_LAYER
     &      )
         IF (IERR.NE.I_NORMAL) RETURN
      ENDIF
!
!
!     DETERMINE THE TRANSMISSION AND REFLECTION COEFFICIENTS.
      CALL TRANS_SOURCE_COEFF(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , ISOLIR, L_IR_SOURCE_QUAD
     &   , TAU, SUM, DIFF, LAMBDA, SEC_0
     &   , GAMMA_UP, GAMMA_DOWN
     &   , TRANS, REFLECT, TRANS_0, SOURCE_COEFF
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      RETURN
      END
