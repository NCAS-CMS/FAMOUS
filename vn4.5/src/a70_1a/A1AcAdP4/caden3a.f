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
!+ Subroutine to calculate densities.
!
! Method:
!       This routine calculates the density of air and the molar
!       densities of the broadening species for the self and foreign-
!       broadened continua using the gas law including the effect
!       of water vapour.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM
     &   , WATER_FRAC, P, T, I_TOP
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN
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
!     INCLUDE COMDECKS
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP VERTICAL INDEX
      LOGICAL
     &     L_CONTINUUM
!             CONTINUUM FLAG
      REAL      !, INTENT(IN)
     &     WATER_FRAC(NPD_PROFILE, 0: NPD_LAYER)
!             MASS FRACTION OF WATER
     &   , P(NPD_PROFILE, 0: NPD_LAYER)
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
      REAL      !, INTENT(OUT)
     &     DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             AIR DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF WATER
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF FOREIGN SPECIES
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
!
!
!     FIND THE AIR DENSITY FIRST.
      DO I=I_TOP, N_LAYER
         DO L=1, N_PROFILE
            DENSITY(L, I)=P(L, I)/(R*T(L, I)
     &         *(1.0E+00+C_VIRTUAL*WATER_FRAC(L, I)))
         ENDDO
      ENDDO
!
      IF (L_CONTINUUM) THEN
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               MOLAR_DENSITY_FRN(L, I)=DENSITY(L, I)
     &            *(1.0E+00-WATER_FRAC(L, I))/MOL_WEIGHT_AIR
               MOLAR_DENSITY_WATER(L, I)=DENSITY(L, I)
     &            *WATER_FRAC(L, I)/(EPSILON*MOL_WEIGHT_AIR)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
