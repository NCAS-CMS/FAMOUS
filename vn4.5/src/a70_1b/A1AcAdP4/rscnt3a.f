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
!+ Subroutine to apply a path-length scaling to the continuum.
!
! Method:
!       The scaling function is calculated. This is multpiled by a
!       suitable "amount" of continuum incorporating a broadening
!       density.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Oct. 96     T3E migration: HF functions
!                                   replaced by T3E vec_lib function
!                                   rtor_v      (S.J.Swarbrick)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM
     &   , P, T, L_LAYER, I_TOP
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN
     &   , WATER_FRAC
     &   , AMOUNT_CONTINUUM
     &   , I_FNC
     &   , P_REFERENCE, T_REFERENCE, SCALE_PARAMETER
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
!     INCLUDE COMDECKS
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE SETTING PHYSICAL CONSTANTS.
!
      REAL
     &     MOL_WEIGHT_AIR
!             MOLAR WEIGHT OF DRY AIR
     &   , N2_MASS_FRAC
!             MASS FRACTION OF NITROGEN
!
      PARAMETER(
     &     MOL_WEIGHT_AIR=28.966E-3
     &   , N2_MASS_FRAC=0.781E+00
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE SETTING PARAMETERS FOR CONTINUUM DATA.
!
      INTEGER
     &     IP_SELF_CONTINUUM
!             SELF-BROADENED CONTINUUM
     &   , IP_FRN_CONTINUUM
!             FOREIGN-BROADENED CONTINUUM
     &   , IP_N2_CONTINUUM
!             NITROGEN CONTINUUM
!
      PARAMETER(
     &     IP_SELF_CONTINUUM=1
     &   , IP_FRN_CONTINUUM=2
     &   , IP_N2_CONTINUUM=3
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , I_CONTINUUM
!             CONTINUUM TYPE
     &   , I_FNC
!             SCALING FUNCTION
     &   , I_TOP
!             TOP INDEX OF ARRAYS
      LOGICAL   !, INTENT(IN)
     &     L_LAYER
!             DATA ARE SUPPLIED IN LAYERS
      REAL      !, INTENT(IN)
     &     WATER_FRAC(NPD_PROFILE, 0: NPD_LAYER)
!             MASS FRACTION OF WATER
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
     &   , DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             OVERALL DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF WATER VAPOUR
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF FOREIGN SPECIES
     &   , P_REFERENCE
!             REFERENCE PRESSURE
     &   , T_REFERENCE
!             REFERENCE PRESSURE
     &   , SCALE_PARAMETER(NPD_SCALE_VARIABLE)
!             SCALING PARAMTERS
      REAL      !, INTENT(OUT)
     &     AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER)
!             AMOUNT OF CONTINUUM
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
      REAL    PWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      REAL    TWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
!
!
         DO I=   1, N_LAYER-I_TOP+1    
            DO L=1, N_PROFILE                                           
              PWK(L,I)=P(L, I_TOP+I-1)/P_REFERENCE
            END DO
         END DO
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE   
              PWK(L,I)=PWK(L,I)**SCALE_PARAMETER(1)
            ENDDO                                                       
         ENDDO                
!
      IF (I_FNC.EQ.IP_SCALE_POWER_LAW) THEN
!
         DO I=   1, N_LAYER-I_TOP+1  
            DO L=1, N_PROFILE                                           
              TWK(L,I)=T(L, I_TOP+I-1)/T_REFERENCE
            END DO
         END DO
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE   
              TWK(L,I)=TWK(L,I)**SCALE_PARAMETER(2)
            ENDDO                                                       
         ENDDO                
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               AMOUNT_CONTINUUM(L, I)
     &                       =PWK(L,I-I_TOP+1)*TWK(L,I-I_TOP+1) 
            ENDDO
         ENDDO
      ELSE IF(I_FNC.EQ.IP_SCALE_POWER_QUAD) THEN
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               AMOUNT_CONTINUUM(L, I)
     &            =PWK(L,I-I_TOP+1)
     &            *(1.0E+00+SCALE_PARAMETER(2)*(T(L, I)
     &            /T_REFERENCE-1.0E+00)
     &            +SCALE_PARAMETER(3)*(T(L, I)
     &            /T_REFERENCE-1.0E+00)**2)
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_LAYER) THEN
         IF (I_CONTINUUM.EQ.IP_SELF_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)
     &               *MOLAR_DENSITY_WATER(L, I)*WATER_FRAC(L, I)
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM.EQ.IP_FRN_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)
     &               *MOLAR_DENSITY_FRN(L, I)*WATER_FRAC(L, I)
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM.EQ.IP_N2_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)
     &               *N2_MASS_FRAC*DENSITY(L, I)
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        IF VALUES ARE GIVEN ON LEVELS WE NOW INTERPOLATE TO AVERAGES
!        ACROSS THE LAYER.
         IF (I_CONTINUUM.EQ.IP_SELF_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00
     &               *(AMOUNT_CONTINUUM(L, I)
     &               *MOLAR_DENSITY_WATER(L, I)*WATER_FRAC(L, I)
     &               +AMOUNT_CONTINUUM(L, I-1)
     &               *MOLAR_DENSITY_WATER(L, I-1)*WATER_FRAC(L, I-1))
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM.EQ.IP_FRN_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00
     &               *(AMOUNT_CONTINUUM(L, I)
     &               *MOLAR_DENSITY_FRN(L, I)*WATER_FRAC(L, I)
     &               +AMOUNT_CONTINUUM(L, I-1)
     &               *MOLAR_DENSITY_FRN(L, I-1)*WATER_FRAC(L, I-1))
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM.EQ.IP_N2_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00
     &               *(AMOUNT_CONTINUUM(L, I)
     &               *N2_MASS_FRAC*DENSITY(L, I)
     &               +AMOUNT_CONTINUUM(L, I-1)
     &               *N2_MASS_FRAC*DENSITY(L, I-1))
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!
      RETURN
      END
