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
!+ Subroutine to set the matrix equations for the fluxes.
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
!       4.1             05-03-96                Surface albedo and
!                                               reflection coefficients
!                                               perturbed to avoid
!                                               failure above a
!                                               non-reflecting surface
!                                               or in a strongly
!                                               absorbing atmosphere.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_MATRIX_FULL(N_PROFILE, N_LAYER
     &   , TRANS, REFLECT
     &   , S_DOWN, S_UP
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND, FLUX_INC_DOWN
     &   , SOURCE_GROUND
     &   , A3, B
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
      IMPLICIT NONE
!
!
!      COMDECKS INCLUDED.
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
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION COEFFICIENTS
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)
!             REFLECT COEFFICIENTS
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE FUNCTION
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE REFLECTION COEFFICIENTS
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT REFLECTION COEFFICIENTS
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE)
!             DIRECT FLUX AT SURFACE
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FROM GROUND
      REAL      !, INTENT(OUT)
     &     A3(NPD_PROFILE, 3, 2*NPD_LAYER+2)
!             TRIDIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             RHS OF EQUATIONS
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     PERTURBATION
!             PERTURBATION TO RESTORE CONDITIONING
!
!
!     CODE THE EQUATIONS INTO THE MATRIX:
!     UPPER SURFACE:
      DO L=1, N_PROFILE
         A3(L, 2, 1)=0.0E+00
         A3(L, 3, 1)=1.0E+00
         B(L, 1)=FLUX_INC_DOWN(L)
      ENDDO
!     THE INTERIOR EQUATIONS:
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            PERTURBATION=(1.0E+00-TRANS(L, I)-REFLECT(L, I))
     &         *TOL_MACHINE/(SQRT_TOL_MACHINE+REFLECT(L, I))
            A3(L, 1, 2*I)=1.0E+00
            A3(L, 2, 2*I)=-REFLECT(L, I)-PERTURBATION
            A3(L, 3, 2*I)=-TRANS(L, I)
            B(L, 2*I)=S_UP(L, I)
            A3(L, 1, 2*I+1)=-TRANS(L, I)
            A3(L, 2, 2*I+1)=-REFLECT(L, I)-PERTURBATION
            A3(L, 3, 2*I+1)=1.0E+00
            B(L, 2*I+1)=S_DOWN(L, I)
         ENDDO
      ENDDO
!     LOWER BOUNDARY
      DO L=1, N_PROFILE
         A3(L, 1, 2*N_LAYER+2)=1.0E+00
         A3(L, 2, 2*N_LAYER+2)=-ALBEDO_SURFACE_DIFF(L)
     &      -TOL_MACHINE/(ALBEDO_SURFACE_DIFF(L)+SQRT_TOL_MACHINE)
         B(L, 2*N_LAYER+2)
     &      =(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))
     &      *FLUX_DIRECT_GROUND(L)+SOURCE_GROUND(L)
      ENDDO
!
!
      RETURN
      END
