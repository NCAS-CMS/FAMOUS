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
!+ Subroutine to set the pentadiagonal matrix for the fluxes.
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
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER
     &   , TRANS, REFLECT
     &   , S_DOWN, S_UP
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND, FLUX_INC_DOWN
     &   , SOURCE_GROUND
     &   , A5, B
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
      REAL      !, INTENT(OUT)
     &     A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)
!             SOURCE TERMS FOR EQUATIONS
!
!     DECLARATION OF LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
!     THE TOP BOUNDARY CONDITION:
      DO L=1, N_PROFILE
         A5(L, 4, 2)=0.0E+00
         A5(L, 3, 2)=1.0E+00
         A5(L, 2, 2)=0.0E+00
         A5(L, 1, 2)=0.0E+00
         B(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
!
!     INTERIOR ROWS: ODD AND EVEN ROWS CORRESPOND TO DIFFERENT BOUNDARY
!     CONDITIONS.
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            A5(L, 5, 2*I-1)=0.0E+00
            A5(L, 4, 2*I-1)=0.0E+00
            A5(L, 3, 2*I-1)=-1.0E+00
            A5(L, 2, 2*I-1)=REFLECT(L, I)
            A5(L, 1, 2*I-1)=TRANS(L, I)
            B(L, 2*I-1)=-S_UP(L, I)
            A5(L, 5, 2*I+2)=TRANS(L, I)
            A5(L, 4, 2*I+2)=REFLECT(L, I)
            A5(L, 3, 2*I+2)=-1.0E+00
            A5(L, 2, 2*I+2)=0.0E+00
            A5(L, 1, 2*I+2)=0.0E+00
            B(L, 2*I+2)=-S_DOWN(L, I)
         ENDDO
      ENDDO
!
!     THE SURFACE BOUNDARY CONDITION:
      DO L=1, N_PROFILE
         A5(L, 5, 2*N_LAYER+1)=0.0E+00
         A5(L, 4, 2*N_LAYER+1)=0.0E+00
         A5(L, 3, 2*N_LAYER+1)=1.0E+00
         A5(L, 2, 2*N_LAYER+1)=-ALBEDO_SURFACE_DIFF(L)
         B(L, 2*N_LAYER+1)=SOURCE_GROUND(L)
     &      +(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))
     &      *FLUX_DIRECT_GROUND(L)
      ENDDO
!
!
      RETURN
      END
