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
!+ Subroutine to calculate fluxes in a homogeneous column directly.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             09-04-96                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER
     &   , TRANS, REFLECT
     &   , S_DOWN, S_UP
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND, FLUX_INC_DOWN
     &   , SOURCE_GROUND
     &   , FLUX_TOTAL
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)
!             REFLECTION COEFFICIENT
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD DIFFUSE SOURCE
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD DIFFUSE SOURCE
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FUNCTION OF GROUND
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE)
!             DIRECT FLUX AT
!                     GROUND LEVEL
!
      REAL      !, INTENT(OUT)
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!     DECLARATION OF LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     ALPHA(NPD_PROFILE, NPD_LAYER+1)
!             COMBINED ALBEDO OF LOWER LAYERS
     &   , BETA(NPD_PROFILE, NPD_LAYER)
!             WORKING ARRAY
     &   , GAMMA(NPD_PROFILE, NPD_LAYER)
!             WORKING ARRAY
     &   , H(NPD_PROFILE, NPD_LAYER)
!             WORKING ARRAY
     &   , S_UP_PRIME(NPD_PROFILE, NPD_LAYER+1)
!             MODIFIED UPWARD SOURCE FUNCTION
!
!
!
!     INITIALIZATION AT THE BOTTOM FOR UPWARD ELIMINATION:
      DO L=1, N_PROFILE
         ALPHA(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         S_UP_PRIME(L, N_LAYER+1)=SOURCE_GROUND(L)
     &      +(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))
     &      *FLUX_DIRECT_GROUND(L)
      ENDDO
!
!     ELIMINATING LOOP:
      DO I=N_LAYER, 1, -1
         DO L=1, N_PROFILE
            BETA(L, I)=1.0E+00/(1.0E+00-ALPHA(L, I+1)*REFLECT(L, I))
            GAMMA(L, I)=ALPHA(L, I+1)*TRANS(L, I)
            H(L, I)=S_UP_PRIME(L, I+1)+ALPHA(L, I+1)*S_DOWN(L, I)
            ALPHA(L, I)=REFLECT(L, I)
     &         +BETA(L, I)*GAMMA(L, I)*TRANS(L, I)
            S_UP_PRIME(L, I)=S_UP(L, I)+BETA(L, I)*TRANS(L, I)*H(L, I)
         ENDDO
      ENDDO
!
!     INITIALIZE FOR BACKWARD SUBSTITUTION.
      DO L=1, N_PROFILE
         FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
         FLUX_TOTAL(L, 1)=ALPHA(L, 1)*FLUX_TOTAL(L, 2)+S_UP_PRIME(L, 1)
      ENDDO
!
!     BACKWARD SUBSTITUTION:
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
!           UPWARD FLUX
            FLUX_TOTAL(L, 2*I+1)
     &         =BETA(L, I)*(H(L, I)+GAMMA(L, I)*FLUX_TOTAL(L, 2*I))
!           DOWNWARD FLUX
            FLUX_TOTAL(L, 2*I+2)=S_DOWN(L, I)
     &         +TRANS(L, I)*FLUX_TOTAL(L, 2*I)
     &         +REFLECT(L, I)*FLUX_TOTAL(L, 2*I+1)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
