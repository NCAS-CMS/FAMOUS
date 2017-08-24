C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Subroutine to scale amounts of absorbers.
!
! Method:
!       The mixing ratio is multiplied by a factor determined
!       by the type of scaling selected.
!
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.5             11-06-98                Optimised Code
!                                               (P. Burton)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SCALE_ABSORB(IERR, N_PROFILE, N_LAYER
     &   , GAS_MIX_RATIO, P, T, L_LAYER, I_TOP
     &   , GAS_FRAC_RESCALED
     &   , I_FNC, P_REFERENCE, T_REFERENCE, SCALE_PARAMETER
     &   , L_DOPPLER, DOPPLER_CORRECTION
     &   , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC
     &   , NPD_SCALE_VARIABLE
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
     &   , NPD_SCALE_FNC
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             MAX. NUMBER OF SCALING VARIABLES
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , I_FNC
!             TYPE OF SCALING FUNCTION
     &   , I_TOP
!             UPPERMOST INDEX FOR SCALING (THIS WILL BE 1 FOR FIELDS
!             GIVEN IN LAYERS, AS IN THE UNIFIED MODEL, OR 0 FOR
!             FIELDS GIVEN AT THE BOUNDARIES OF LAYERS)
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             DATA SPECIFIED IN LAYERS
     &   , L_DOPPLER
!             FLAG FOR DOPPLER TERM
      REAL      !, INTENT(IN)
     &     GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER)
!             MASS MIXING RATIO OF GAS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
     &   , P_REFERENCE
!             REFERENCE PRESSURE
     &   , T_REFERENCE
!             REFERENCE TEMPERATURE
     &   , SCALE_PARAMETER(NPD_SCALE_VARIABLE)
!             SCALING PARAMTERS
     &   , DOPPLER_CORRECTION
!             DOPPLER-BROADENING CORRECTION
      REAL      !, INTENT(OUT)
     &     GAS_FRAC_RESCALED(NPD_PROFILE, 0: NPD_LAYER)
!             MASS FRACTION OF GAS
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
      REAL
     &     PRESSURE_OFFSET
!             OFFSET TO PRESSURE
      REAL    PWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      REAL    TWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      REAL TMP, T_INV,P_REF_OFF_INV

!
!     SET THE OFFSET TO THE PRESSURE FOR THE DOPPLER CORRECTION.
      IF (L_DOPPLER) THEN
         PRESSURE_OFFSET=DOPPLER_CORRECTION
      ELSE
         PRESSURE_OFFSET=0.0E+00
      ENDIF

      T_INV = 1.0/T_REFERENCE
      P_REF_OFF_INV = 1.0/(P_REFERENCE+PRESSURE_OFFSET)

!     THE ARRAY GAS_FRAC_RESCALED IS USED INITIALLY TO HOLD ONLY THE
!     SCALING FUNCTIONS, AND ONLY LATER IS IT MULTIPLIED BY THE
!     MIXING RATIOS
!
      IF (I_FNC.EQ.IP_SCALE_POWER_LAW) THEN
!
         IF(L_DOPPLER) THEN
           DO I=   1, N_LAYER-I_TOP+1
              DO L=1, N_PROFILE
                PWK(L,I)=(P(L,I_TOP+I-1)+PRESSURE_OFFSET)
     &                 *P_REF_OFF_INV
                TWK(L,I)=T(L,I_TOP+I-1)*T_INV
              END DO
           END DO
         ELSE
            DO I=   1, N_LAYER-I_TOP+1
               DO L=1, N_PROFILE
                  PWK(L,I)=P(L,I_TOP+I-1)
     &                 *P_REF_OFF_INV
                  TWK(L,I)=T(L,I_TOP+I-1)*T_INV
               END DO
            END DO
         END IF
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=PWK(L,I)**SCALE_PARAMETER(1)
              TWK(L,I)=TWK(L,I)**SCALE_PARAMETER(2)
            ENDDO
         ENDDO
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               GAS_FRAC_RESCALED(L, I)
     &                       =PWK(L,I-I_TOP+1)*TWK(L,I-I_TOP+1)
            ENDDO
         ENDDO
      ELSE IF (I_FNC.EQ.IP_SCALE_FNC_NULL) THEN
         RETURN
      ELSE IF (I_FNC.EQ.IP_SCALE_POWER_QUAD) THEN
!
         IF(L_DOPPLER) THEN
           DO I=   1, N_LAYER-I_TOP+1
              DO L=1, N_PROFILE
                PWK(L,I)=(P(L,I_TOP+I-1)+PRESSURE_OFFSET)
     &                 *P_REF_OFF_INV
              END DO
           END DO
         ELSE
            DO I=   1, N_LAYER-I_TOP+1
               DO L=1, N_PROFILE
                  PWK(L,I)=P(L,I_TOP+I-1)
     &                 *P_REF_OFF_INV
               END DO
            END DO
         END IF
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=PWK(L,I)**SCALE_PARAMETER(1)
            ENDDO
         ENDDO
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TMP = T(L,I)*T_INV - 1.0
               GAS_FRAC_RESCALED(L, I)=PWK(L,I-I_TOP+1)*
     &    (1.0E+00+TMP*SCALE_PARAMETER(2)
     &            +SCALE_PARAMETER(3)*TMP*TMP)
            ENDDO
         ENDDO
      ELSE IF (I_FNC.EQ.IP_SCALE_DOPPLER_QUAD) THEN
!        THERE IS NO DOPPLER TERM HERE SINCE IT IS IMPLICITLY INCLUDED
!        IN THE SCALING.
!
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=(P(L,I_TOP+I-1)+SCALE_PARAMETER(2))
     &                       /(P_REFERENCE+SCALE_PARAMETER(2))
            END DO
         END DO
         DO I=1,N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=PWK(L,I)**SCALE_PARAMETER(1)
            ENDDO
         ENDDO
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TMP = T(L,I)*T_INV - 1.0
               GAS_FRAC_RESCALED(L, I)=PWK(L,I-I_TOP+1)
     &            *(1.0E+00
     &            +TMP*SCALE_PARAMETER(3)
     &            +SCALE_PARAMETER(4)*TMP*TMP)
            ENDDO
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: AN ILLEGAL TYPE OF SCALING HAS BEEN GIVEN.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     MULTIPLY BY THE MIXING RATIO AND LIMIT NEGATIVE SCALINGS.
      IF (L_LAYER) THEN
         DO I=N_LAYER, 1, -1
            DO L=1, N_PROFILE
               GAS_FRAC_RESCALED(L, I)=MAX(REAL(0.0E+00)
     &            , GAS_FRAC_RESCALED(L, I)*GAS_MIX_RATIO(L, I))
            ENDDO
         ENDDO
      ELSE
!        CONVERT TO VALUES IN LAYERS.
         DO I=N_LAYER, 1, -1
            DO L=1, N_PROFILE
               GAS_FRAC_RESCALED(L, I)
     &            =0.5E+00*(GAS_FRAC_RESCALED(L, I-1)
     &            *GAS_MIX_RATIO(L, I-1)
     &            +GAS_FRAC_RESCALED(L, I)*GAS_MIX_RATIO(L, I))
               GAS_FRAC_RESCALED(L, I)
     &            =MAX(REAL(0.0E+00), GAS_FRAC_RESCALED(L, I))
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
