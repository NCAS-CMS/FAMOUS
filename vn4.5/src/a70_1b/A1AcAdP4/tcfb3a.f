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
!+ Subroutine to calculate basic coefficients in two-stream equations.
!
! Method:
!       Depending on the two-stream equations employed, the
!       appropriate coefficients for the fluxes are calculated.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                PIFM85 restored to
!                                               original form. PIFM80
!                                               introduced.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_BASIC(IERR
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , I_2STREAM
     &   , ASYMMETRY, OMEGA
     &   , SUM, DIFF
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
!     MODULE TO SET DIFFUSIVITY FOR ELSASSER'S SCHEME.
!
      REAL
     &     ELSASSER_FACTOR
!             DIFFUSIVITY FACTOR FOR ELSASSER'S SCHEME
!
      PARAMETER(
     &     ELSASSER_FACTOR=1.66E+00
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
     &   , I_2STREAM
!             TWO STREAM SCHEME
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL      !, INTENT(IN)
     &     ASYMMETRY(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRY FACTOR
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL      !, INTENT(OUT)
     &     SUM(NPD_PROFILE, NPD_LAYER)
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     ROOT_3
!             SQUARE ROOT OF 3
!
!
      PARAMETER(
     &     ROOT_3=1.7320508075688772E+00
     &   )
!
!
!
!
      IF (I_2STREAM.EQ.IP_EDDINGTON) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=1.5E+00*(1.0E+00
     &            -OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_ELSASSER) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=ELSASSER_FACTOR
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=ELSASSER_FACTOR*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_DISCRETE_ORD) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=ROOT_3*(1.0E+00
     &            -OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=ROOT_3*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_PIFM85) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_2S_TEST) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=1.5E+00
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=1.5E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_HEMI_MEAN) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00
     &            *(1.0E+00-OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_PIFM80) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
     &            -0.5E+00*OMEGA(L, I)
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM.EQ.IP_IFM) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: THE IMPROVED FLUX METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM.EQ.IP_ZDK_FLUX) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: ZDUNKOWSKI''S FLUX METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM.EQ.IP_KRSCHG_FLUX) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: KERSCHGEN''S FLUX METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM.EQ.IP_COAKLEY_CHYLEK_1) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: COAKLEY-CHYLEK''S FIRST METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM.EQ.IP_COAKLEY_CHYLEK_2) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: COAKLEY-CHYLEK''S SECOND METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM.EQ.IP_MEADOR_WEAVER) THEN
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: MEADOR & WEAVER''S METHOD HAS '
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE
         WRITE(IU_ERR, '(/A)')
     &      '*** ERROR: AN ILLEGAL PARAMETER HAS BEEN SUPPLIED '
     &      //'TO DEFINE THE 2-STREAM SCHEME.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      RETURN
      END
