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
!+ Subroutine to calculate the infra-red radiance ignoring scattering.
!
! Method:
!       Using the secant of the ray transmission coefficients for
!       each layer may be defined and source terms may be calculated.
!       The upward and downward radiances are integrated along
!       their paths.
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
      SUBROUTINE MONOCHROMATIC_IR_RADIANCE(N_PROFILE, N_LAYER
     &   , L_NET
     &   , TAU
     &   , RAD_INC_DOWN
     &   , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF
     &   , SECANT_RAY
     &   , RADIANCE
     &   , NPD_PROFILE, NPD_LAYER
     &         )
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES.
      REAL      !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTHS OF LAYERS
     &   , RAD_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD RADIANCE
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FUNCTION OF GROUND
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE IN PLANCK FUNCTION
     &   , SECANT_RAY
!             SECANT OF ANGLE WITH VERTICAL
      REAL      !, INTENT(OUT)
     &     RADIANCE(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE RADIANCE
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSIVITIES
     &   , SOURCE_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
     &   , SOURCE_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE FUNCTION
!
!
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            TRANS(L, I)=EXP(-SECANT_RAY*TAU(L, I))
         ENDDO
      ENDDO
!
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            SOURCE_UP(L, I)=(1.0E+00-TRANS(L, I)+SQRT_TOL_MACHINE)
     &         *DIFF_PLANCK(L, I)
     &         /(SECANT_RAY*TAU(L, I)+SQRT_TOL_MACHINE)
            SOURCE_DOWN(L, I)=-SOURCE_UP(L, I)
         ENDDO
      ENDDO
!
!     DOWNWARD RADIANCE.
      DO L=1, N_PROFILE
         RADIANCE(L, 2)=RAD_INC_DOWN(L)
      ENDDO
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            RADIANCE(L, 2*I+2)=TRANS(L, I)*RADIANCE(L, 2*I)
     &         +SOURCE_DOWN(L, I)
         ENDDO
      ENDDO
!
!     UPWARD RADIANCE.
      DO L=1, N_PROFILE
         RADIANCE(L, 2*N_LAYER+1)=SOURCE_GROUND(L)
     &      +ALBEDO_SURFACE_DIFF(L)*RADIANCE(L, 2*N_LAYER+2)
      ENDDO
      DO I=N_LAYER, 1, -1
         DO L=1, N_PROFILE
            RADIANCE(L, 2*I-1)=TRANS(L, I)*RADIANCE(L, 2*I+1)
     &         +SOURCE_UP(L, I)
         ENDDO
      ENDDO
!
!     REDUCE TO A NET RADIANCE IF THIS IS REQUIRED.
      IF (L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               RADIANCE(L, I+1)=RADIANCE(L, 2*I+2)
     &            -RADIANCE(L, 2*I+1)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
