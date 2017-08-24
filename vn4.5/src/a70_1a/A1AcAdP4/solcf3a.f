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
!+ Subroutine to calculate the basic coefficients for the solar beam.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFLT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.2             08-08-96                Extra two-stream option
!                                               added.
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
      SUBROUTINE SOLAR_COEFFICIENT_BASIC(IERR
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , OMEGA, ASYMMETRY, SEC_0
     &   , I_2STREAM
     &   , SUM, DIFF, LAMBDA
     &   , GAMMA_UP, GAMMA_DOWN
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
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
!     MODULE TO SET THE DEFINING NUMBERS FOR THE TWO-STREAM
!     SCHEMES.
!
      INTEGER
     &     IP_EDDINGTON
!             EDDINGTON APPROXIMATION
     &   , IP_DISCRETE_ORD
!             DISCRETE ORDINATE METHOD
     &   , IP_IFM
!             IMPROVED FLUX METHOD
     &   , IP_PIFM85
!             PRACTICAL IMPROVED FLUX METHOD
!             (VERSION OF ZDUNKOWSKI ET AL. 1985)
     &   , IP_ZDK_FLUX
!             ZDUNKOWSKI'S FLUX METHOD
     &   , IP_KRSCHG_FLUX
!             KERSCHGEN'S FLUX METHOD
     &   , IP_COAKLEY_CHYLEK_1
!             COAKLEY & CHYLEK'S 1ST METHOD
     &   , IP_COAKLEY_CHYLEK_2
!             COAKLEY & CHYLEK'S 2ND METHOD
     &   , IP_MEADOR_WEAVER
!             MEADOR & WEAVER'S METHOD
     &   , IP_ELSASSER
!             ELSASSER'S DIFFUSIVITY SCHEME
     &   , IP_2S_TEST
!             USER'S DEFINED TEST APPROXIMATION.
     &   , IP_HEMI_MEAN
!             HEMISPHERIC MEAN APPROXIMATION.
     &   , IP_PIFM80
!             PRACTICAL IMPROVED FLUX METHOD
!             (VERSION OF ZDUNKOWSKI ET AL. 1980)
!
      PARAMETER(
     &     IP_EDDINGTON=2
     &   , IP_DISCRETE_ORD=4
     &   , IP_IFM=5
     &   , IP_PIFM85=6
     &   , IP_ZDK_FLUX=7
     &   , IP_KRSCHG_FLUX=8
     &   , IP_COAKLEY_CHYLEK_1=9
     &   , IP_COAKLEY_CHYLEK_2=10
     &   , IP_MEADOR_WEAVER=11
     &   , IP_ELSASSER=12
     &   , IP_2S_TEST=14
     &   , IP_HEMI_MEAN=15
     &   , IP_PIFM80=16
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
!
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST
!             FIRST LAYER TO CONSIDER
     &   , I_2STREAM
!             TWO-STREAM SCHEME
!
      REAL      !, INTENT(IN)
     &     OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRY
     &   , SEC_0(NPD_PROFILE)
!             SECANT OF SOLAR ZENITH ANGLE
     &   , SUM(NPD_PROFILE, NPD_LAYER)
!             SUM OF TWO-STREAM COEFFICIENTS
     &   , DIFF(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE OF TWO-STREAM COEFFICIENTS
     &   , LAMBDA(NPD_PROFILE, NPD_LAYER)
!             LAMBDA
!
!     BASIC TWO-STREAM COEFFICIENTS:
      REAL      !, INTENT(OUT)
     &     GAMMA_UP(NPD_PROFILE, NPD_LAYER)
!             COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             COEFFICIENT FOR DOWNWAD RADIATION
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
     &   , KK
!             TEMPORARY VARIABLE
     &   , N_INDEX
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
      REAL
     &     KSI_0(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE IN SOLAR SCATTERING FRACTIONS
     &   , TEST_ARRAY(NPD_PROFILE)
!             ARRAY TO TEST
     &   , FACTOR
!             TEMPORARY VARIABLE
      REAL
     &     ROOT_3
!             SQUARE ROOT OF 3
      PARAMETER(
     &     ROOT_3=1.7320508075688772E+00
     &   )
!
!     SUBROUTINES CALLED:
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
Cfpp$ NODEPCHK R
!
!
!
!     IF LAMBDA IS TOO CLOSE TO SEC_0 IT MUST BE PERTURBED
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            TEST_ARRAY(L)=ABS(LAMBDA(L, I)-SEC_0(L))
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (TEST_ARRAY(L).LT.TOL_DIV) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
CDIR$ IVDEP
         DO K=1, N_INDEX
            KK=INDEX(K)
            SUM(KK, I)=(1.0E+00+TOL_DIV)*SUM(KK, I)
            DIFF(KK, I)=(1.0E+00+TOL_DIV)*DIFF(KK, I)
            LAMBDA(KK, I)=(1.0E+00+TOL_DIV)*LAMBDA(KK, I)
         ENDDO
      ENDDO
!
      IF ( (I_2STREAM.EQ.IP_EDDINGTON).OR.
     &     (I_2STREAM.EQ.IP_ELSASSER).OR.
     &     (I_2STREAM.EQ.IP_PIFM85).OR.
     &     (I_2STREAM.EQ.IP_2S_TEST).OR.
     &     (I_2STREAM.EQ.IP_HEMI_MEAN).OR.
     &     (I_2STREAM.EQ.IP_PIFM80) ) THEN
!
          DO I=I_LAYER_FIRST, I_LAYER_LAST
             DO L=1, N_PROFILE
                KSI_0(L, I)=1.5E+00*ASYMMETRY(L, I)/SEC_0(L)
             ENDDO
          ENDDO
!
       ELSE IF (I_2STREAM.EQ.IP_DISCRETE_ORD) THEN
!
          DO I=I_LAYER_FIRST, I_LAYER_LAST
             DO L=1, N_PROFILE
                KSI_0(L, I)=ROOT_3*ASYMMETRY(L, I)/SEC_0(L)
             ENDDO
          ENDDO
!
       ELSE
!
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: AN ILLEGAL SOLAR TWO-STREAM SCHEME HAS '
     &      //'BEEN SELECTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!     DETERMINE THE BASIC SOLAR COEFFICIENTS FOR THE
!     TWO-STREAM EQUATIONS.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            FACTOR=0.5E+00*OMEGA(L, I)*SEC_0(L)
     &         /((LAMBDA(L, I)-SEC_0(L))*(LAMBDA(L, I)+SEC_0(L)))
            GAMMA_UP(L, I)=FACTOR*(SUM(L, I)-SEC_0(L)
     &         -KSI_0(L, I)*(DIFF(L, I)-SEC_0(L)))
            GAMMA_DOWN(L, I)=FACTOR*(SUM(L, I)+SEC_0(L)
     &         +KSI_0(L, I)*(DIFF(L, I)+SEC_0(L)))
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
