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
!+ Subroutine implementing an endekadiagonal solution of the equations.
!
! Method:
!       This subroutine calls subroutines to set the endekadiagonal
!       matrix for the clear and cloudy upward and downward fluxes
!       and to solve those equations. It is kept separate to avoid
!       allocating space for so large a matrix in the upper routine
!       when that is used in a GCM and this routine is very
!       unlikely to be used.
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
      SUBROUTINE MIX_COLUMN_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T_FREE, R_FREE, S_DOWN_FREE, S_UP_FREE
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD
     &   , G_M, G_P, B_M, B_P
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND
     &   , FLUX_TOTAL
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      REAL      !, INTENT(IN)
     &     T_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE TRANSMISSION
     &   , R_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE REFLECTION
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DOWNWARD SOURCE FUNCTION
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE UPWARD SOURCE FUNCTION
     &   , T_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY TRANSMISSION
     &   , R_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY REFLECTION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD CLOUDY SOURCE FUNCTION
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             UPWARD CLOUDY SOURCE FUNCTION
      REAL      !, INTENT(IN)
     &     B_M(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , B_P(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_M(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_P(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
      REAL      !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT TOTAL FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FROM GROUND
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE)
!             DIRECT FLUX AT GROUND
      REAL      !, INTENT(OUT)
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , N_EQUATION
!             NUMBER OF EQUATIONS
      REAL
     &     A11(NPD_PROFILE, 11, 4*NPD_LAYER+4)
!             MATRIX TO BE SOLVED
     &   , B(NPD_PROFILE, 4*NPD_LAYER+4)
!             RIGHT-HAND SIDE OF EQUATION
     &   , X(NPD_PROFILE, 4*NPD_LAYER+4)
!             SOLUTION TO MATRIX EQUATION
     &   , WORK(NPD_PROFILE)
!             WORKING ARRAY
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     MIX_MATRIX_ELEM, BAND_SOLVER
!
!
!
!     ASSIGN THE ELEMENTS OF THE MATRIX A11 AND THE VECTOR B.
      CALL MIX_MATRIX_ELEM(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T_FREE, R_FREE, S_DOWN_FREE, S_UP_FREE
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD
     &   , G_M, G_P, B_M, B_P
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND
     &   , A11, B
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
      N_EQUATION=4*N_LAYER+4
      CALL BAND_SOLVER(N_PROFILE, N_EQUATION
     &   , 5, 5
     &   , A11, B
     &   , X
     &   , NPD_PROFILE, 4*NPD_LAYER+4
     &   , WORK
     &   )
!
!     PICK OUT THE SUMMED FLUXES FROM THE LONG VECTOR
      DO I=1, 2*N_LAYER+2
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=X(L, 2*I)
         ENDDO
      ENDDO
!
!
      RETURN
      END
