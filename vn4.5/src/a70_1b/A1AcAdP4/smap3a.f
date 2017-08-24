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
!+ Subroutine to set moist aerosol properties independent of bands.
!
! Method:
!       The mean relative humidities are calculated and pointers to
!       the lookup tables are set.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.3             17-12-96                Code extended to permit
!                                               use with both moist
!                                               and dry aerosols.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_MOIST_AEROSOL_PROPERTIES(IERR
     &   , N_PROFILE, N_LAYER
     &   , L_LAYER, N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY
     &   , WATER_MIX_RATIO, T, P, DELTA_HUMIDITY
     &   , MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER
     &   , NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES
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
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOLS
!
!     INCLUDE COMDECKS.
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_AEROSOL
!             NUMBER OF AEROSOL SPECIES
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATIONS OF AEROSOL
!             SPECIES
     &   , NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBER OF HUMIDITY VALUES
      INTEGER   !, INTENT(OUT)
     &     I_HUMIDITY_POINTER(NPD_PROFILE, NPD_LAYER)
!             POINTERS TO LOOK-UP TABLES
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             LAYER FLAG
      REAL      !, INTENT(IN)
     &     WATER_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER)
!             MIXING RATIO OF WATER VAPOUR
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURES
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURES
      REAL      !, INTENT(OUT)
     &     MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)
!             MEAN HUMIDITIES OF LAYERS
     &   , DELTA_HUMIDITY
!             INCREMENT IN HUMIDITY
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , NHUMIDITY_COMMON
!             COMMON NUMBER OF HUMIDITIES FOR MOIST AEROSOLS
      REAL
     &     MIX_RATIO_SAT(NPD_PROFILE, 0: NPD_LAYER)
!             SATURATED HUMIDITY MIXING RATIO
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     QSAT_WAT
!
!
!
!     SET UP ARRAY OF POINTERS TO INCLUDE THE EFFECTS OF HUMIDITY.
!     CALCULATE THE SATURATED MIXING RATIO.
      DO I=1, N_LAYER
         CALL QSAT_WAT(MIX_RATIO_SAT(1, I), T(1, I), P(1, I)
     &      , N_PROFILE)
      ENDDO
!
!     DETERMINE THE NUMBER OF HUMIDITIES TO BE USED FOR MOIST
!     AEROSOLS. THIS MUST BE THE SAME FOR ALL MOIST AEROSOLS
!     IN THE CURRENT VERSION OF THE CODE.
      NHUMIDITY_COMMON=0
      DO J=1, N_AEROSOL
         IF (I_AEROSOL_PARAMETRIZATION(J).EQ.IP_AEROSOL_PARAM_MOIST)
     &         THEN
            IF (NHUMIDITY(J).GT.0) THEN
!              SET THE ACTUAL COMMON VALUE.
               IF (NHUMIDITY_COMMON.EQ.0) THEN
                  NHUMIDITY_COMMON=NHUMIDITY(J)
               ELSE IF (NHUMIDITY(J).NE.NHUMIDITY_COMMON) THEN
!                 THERE IS AN INCONSISTENCY.
                  WRITE(IU_ERR, '(/A)')
     &               '***ERROR: THE LOOK-UP TABLES FOR MOIST AEROSOLS '
     &               , 'ARE OF DIFFERENT SIZES. THIS IS NOT PERMITTED.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!     THE LOOK-UP TABLE IS ASSUMED TO BE UNIFORM IN HUMIDITY.
      DELTA_HUMIDITY=1.0E+00/(REAL(NHUMIDITY_COMMON)-1.0E+00)
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            MEAN_REL_HUMIDITY(L, I)
     &         =WATER_MIX_RATIO(L, I)*(1.0E+00-MIX_RATIO_SAT(L, I))
     &         /((1.0E+00-WATER_MIX_RATIO(L, I))*MIX_RATIO_SAT(L, I))
!           CHECK THAT THE MEAN RELATIVE HUMIDITY
!           DOES NOT EXCEED OR EQUAL 1.0.
            MEAN_REL_HUMIDITY(L, I)=MIN(MEAN_REL_HUMIDITY(L, I)
     &        , 0.99999)
            I_HUMIDITY_POINTER(L, I)=1
     &         +INT(MEAN_REL_HUMIDITY(L, I)*(NHUMIDITY_COMMON-1))
         ENDDO
      ENDDO
!
!
!

      RETURN
      END
