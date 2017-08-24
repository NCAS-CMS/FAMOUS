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
!+ Subroutine to set the coefficients of the endekadiagonal matrix.
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
      SUBROUTINE MIX_MATRIX_ELEM(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T_FREE, R_FREE, S_DOWN_FREE, S_UP_FREE
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD
     &   , G_M, G_P, B_M, B_P
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , FLUX_DIRECT_GROUND
     &   , A, B
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
     &   , N_LAYER
     &   , N_CLOUD_TOP
      REAL      !, INTENT(IN)
     &     T_FREE(NPD_PROFILE, NPD_LAYER)
     &   , R_FREE(NPD_PROFILE, NPD_LAYER)
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
     &   , T_CLOUD(NPD_PROFILE, NPD_LAYER)
     &   , R_CLOUD(NPD_PROFILE, NPD_LAYER)
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
      REAL      !, INTENT(IN)
     &     B_M(NPD_PROFILE, 0: NPD_LAYER)
     &   , B_P(NPD_PROFILE, 0: NPD_LAYER)
     &   , G_M(NPD_PROFILE, 0: NPD_LAYER)
     &   , G_P(NPD_PROFILE, 0: NPD_LAYER)
      REAL      !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)
     &   , SOURCE_GROUND(NPD_PROFILE)
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE)
      REAL      !, INTENT(OUT)
     &     A(NPD_PROFILE, 11, 4*NPD_LAYER+4)
     &   , B(NPD_PROFILE, 4*NPD_LAYER+4)
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
     &   , L
      REAL
     &     TP(NPD_PROFILE, NPD_LAYER)
     &   , TM(NPD_PROFILE, NPD_LAYER)
     &   , RP(NPD_PROFILE, NPD_LAYER)
     &   , RM(NPD_PROFILE, NPD_LAYER)
!
!
!
!     PRECALCULATE THE MEAN SUMMED AND DIFFERENCED TRANSMISSIVITIES
!     AND REFLECTION COEFFICIENTS FOR EFFICIENCY.
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            TP(L, I)=0.5E+00*(T_FREE(L, I)+T_CLOUD(L, I))
            TM(L, I)=0.5E+00*(T_FREE(L, I)-T_CLOUD(L, I))
            RP(L, I)=0.5E+00*(R_FREE(L, I)+R_CLOUD(L, I))
            RM(L, I)=0.5E+00*(R_FREE(L, I)-R_CLOUD(L, I))
         ENDDO
      ENDDO
!
!     CONDITIONS ON THE DOWNWARD FLUX ENTERING THE TOP LAYER:
      DO L=1, N_PROFILE
!        INCIDENT DIFFERENCE FLUX:
         A(L, 11, 3)=0.0E+00
         A(L, 10, 3)=0.0E+00
         A(L, 9, 3)=0.0E+00
         A(L, 8, 3)=0.0E+00
         A(L, 7, 3)=0.0E+00
         A(L, 6, 3)=1.0E+00
         A(L, 5, 3)=0.0E+00
         A(L, 4, 3)=0.0E+00
         A(L, 3, 3)=0.0E+00
         A(L, 2, 3)=0.0E+00
         A(L, 1, 3)=0.0E+00
         B(L, 3)=FLUX_INC_DOWN(L)
!        INCIDENT SUMMED FLUX:
         A(L, 11, 4)=0.0E+00
         A(L, 10, 4)=0.0E+00
         A(L, 9, 4)=0.0E+00
         A(L, 8, 4)=0.0E+00
         A(L, 7, 4)=0.0E+00
         A(L, 6, 4)=1.0E+00
         A(L, 5, 4)=0.0E+00
         A(L, 4, 4)=0.0E+00
         A(L, 3, 4)=0.0E+00
         A(L, 2, 4)=0.0E+00
         A(L, 1, 4)=0.0E+00
         B(L, 4)=FLUX_INC_DOWN(L)
      ENDDO
!
!     ASSIGN SETS OF EQUATIONS FOR EACH LAYER:
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
!
!           DIFFERENCE IN FLUXES JUST BELOW THE TOP OF THE LAYER:
            A(L, 11, 4*I-3)=0.0E+00
            A(L, 10, 4*I-3)=0.0E+00
            A(L, 9, 4*I-3)=0.0E+00
            A(L, 8, 4*I-3)=0.0E+00
            A(L, 7, 4*I-3)=0.0E+00
            A(L, 6, 4*I-3)=-1.0E+00
            A(L, 5, 4*I-3)=1.0E+00
            A(L, 4, 4*I-3)=0.0E+00
            A(L, 3, 4*I-3)=0.0E+00
            A(L, 2, 4*I-3)=0.0E+00
            A(L, 1, 4*I-3)=0.0E+00
            B(L, 4*I-3)=0.0E+00
!
!           SUMMED FLUX AT THE TOP OF THE LAYER:
            A(L, 11, 4*I-2)=0.0E+00
            A(L, 10, 4*I-2)=0.0E+00
            A(L, 9, 4*I-2)=0.0E+00
            A(L, 8, 4*I-2)=0.0E+00
            A(L, 7, 4*I-2)=0.0E+00
            A(L, 6, 4*I-2)=-1.0E+00
            A(L, 5, 4*I-2)=0.0E+00
            A(L, 4, 4*I-2)=R_FREE(L, I)
            A(L, 3, 4*I-2)=0.0E+00
            A(L, 2, 4*I-2)=T_FREE(L, I)
            A(L, 1, 4*I-2)=0.0E+00
            B(L, 4*I-2)=-S_UP_FREE(L, I)
!
!           DIFFERENCE FLUX JUST ABOVE THE BASE OF THE LAYER:
            A(L, 11, 4*I+3)=0.0E+00
            A(L, 10, 4*I+3)=0.0E+00
            A(L, 9, 4*I+3)=0.0E+00
            A(L, 8, 4*I+3)=0.0E+00
            A(L, 7, 4*I+3)=0.0E+00
            A(L, 6, 4*I+3)=-1.0E+00
            A(L, 5, 4*I+3)=1.0E+00
            A(L, 4, 4*I+3)=0.0E+00
            A(L, 3, 4*I+3)=0.0E+00
            A(L, 2, 4*I+3)=0.0E+00
            A(L, 1, 4*I+3)=0.0E+00
            B(L, 4*I+3)=0.0E+00
!
!           SUMMED FLUXES AT THE BASE OF THE LAYER:
            A(L, 11, 4*I+4)=0.0E+00
            A(L, 10, 4*I+4)=T_FREE(L, I)
            A(L, 9, 4*I+4)=0.0E+00
            A(L, 8, 4*I+4)=R_FREE(L, I)
            A(L, 7, 4*I+4)=0.0E+00
            A(L, 6, 4*I+4)=-1.0E+00
            A(L, 5, 4*I+4)=0.0E+00
            A(L, 4, 4*I+4)=0.0E+00
            A(L, 3, 4*I+4)=0.0E+00
            A(L, 2, 4*I+4)=0.0E+00
            A(L, 1, 4*I+4)=0.0E+00
            B(L, 4*I+4)=-S_DOWN_FREE(L, I)
         ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
!
!           DIFFERENCE IN FLUXES JUST BELOW THE TOP OF THE LAYER:
            A(L, 11, 4*I-3)=0.0E+00
            A(L, 10, 4*I-3)=0.0E+00
            A(L, 9, 4*I-3)=0.0E+00
            A(L, 8, 4*I-3)=0.0E+00
            A(L, 7, 4*I-3)=0.0E+00
            A(L, 6, 4*I-3)=-1.0E+00
            A(L, 5, 4*I-3)=0.0E+00
            A(L, 4, 4*I-3)=0.5E+00*G_M(L, I-1)*RP(L, I)
            A(L, 3, 4*I-3)=RM(L, I)+0.5E+00*G_P(L, I-1)*RP(L, I)
            A(L, 2, 4*I-3)=0.5E+00*B_M(L, I)*TP(L, I)
            A(L, 1, 4*I-3)=TM(L, I)+0.5E+00*B_P(L, I)*TP(L, I)
            B(L, 4*I-3)=-S_UP_FREE(L, I)+S_UP_CLOUD(L, I)
!
!           SUMMED FLUX AT THE TOP OF THE LAYER:
            A(L, 11, 4*I-2)=0.0E+00
            A(L, 10, 4*I-2)=0.0E+00
            A(L, 9, 4*I-2)=0.0E+00
            A(L, 8, 4*I-2)=0.0E+00
            A(L, 7, 4*I-2)=0.0E+00
            A(L, 6, 4*I-2)=-1.0E+00
            A(L, 5, 4*I-2)=0.5E+00*G_M(L, I-1)*RM(L, I)
            A(L, 4, 4*I-2)=RP(L, I)+0.5E+00*G_P(L, I-1)*RM(L, I)
            A(L, 3, 4*I-2)=0.5E+00*B_M(L, I)*TM(L, I)
            A(L, 2, 4*I-2)=TP(L, I)+0.5E+00*B_P(L, I)*TM(L, I)
            A(L, 1, 4*I-2)=0.0E+00
            B(L, 4*I-2)=-S_UP_FREE(L, I)-S_UP_CLOUD(L, I)
!
!           DIFFERENCE FLUX JUST ABOVE THE BASE OF THE LAYER:
            A(L, 11, 4*I+3)=0.0E+00
            A(L, 10, 4*I+3)=0.5E+00*G_M(L, I-1)*TP(L, I)
            A(L, 9, 4*I+3)=TM(L, I)+0.5E+00*G_P(L, I-1)*TP(L, I)
            A(L, 8, 4*I+3)=0.5E+00*B_M(L, I)*RP(L, I)
            A(L, 7, 4*I+3)=RM(L, I)+0.5E+00*B_P(L, I)*RP(L, I)
            A(L, 6, 4*I+3)=-1.0E+00
            A(L, 5, 4*I+3)=0.0E+00
            A(L, 4, 4*I+3)=0.0E+00
            A(L, 3, 4*I+3)=0.0E+00
            A(L, 2, 4*I+3)=0.0E+00
            A(L, 1, 4*I+3)=0.0E+00
            B(L, 4*I+3)=-S_DOWN_FREE(L, I)+S_DOWN_CLOUD(L, I)
!
!           SUMMED FLUXES AT THE BASE OF THE LAYER:
            A(L, 11, 4*I+4)=0.5E+00*G_M(L, I-1)*TM(L, I)
            A(L, 10, 4*I+4)=TP(L, I)+0.5E+00*G_P(L, I-1)*TM(L, I)
            A(L, 9, 4*I+4)=0.5E+00*B_M(L, I)*RM(L, I)
            A(L, 8, 4*I+4)=RP(L, I)+0.5E+00*B_P(L, I)*RM(L, I)
            A(L, 7, 4*I+4)=0.0E+00
            A(L, 6, 4*I+4)=-1.0E+00
            A(L, 5, 4*I+4)=0.0E+00
            A(L, 4, 4*I+4)=0.0E+00
            A(L, 3, 4*I+4)=0.0E+00
            A(L, 2, 4*I+4)=0.0E+00
            A(L, 1, 4*I+4)=0.0E+00
            B(L, 4*I+4)=-S_DOWN_FREE(L, I)-S_DOWN_CLOUD(L, I)
!
         ENDDO
      ENDDO
!
!     CONDITIONS AT THE SURFACE:
      DO L=1, N_PROFILE
!
         A(L, 11, 4*N_LAYER+1)=0.0E+00
         A(L, 10, 4*N_LAYER+1)=0.0E+00
         A(L, 9, 4*N_LAYER+1)=0.0E+00
         A(L, 8, 4*N_LAYER+1)=0.0E+00
         A(L, 7, 4*N_LAYER+1)=0.0E+00
         A(L, 6, 4*N_LAYER+1)=1.0E+00
         A(L, 5, 4*N_LAYER+1)=-1.0E+00
         A(L, 4, 4*N_LAYER+1)=0.0E+00
         A(L, 3, 4*N_LAYER+1)=0.0E+00
         A(L, 2, 4*N_LAYER+1)=0.0E+00
         A(L, 1, 4*N_LAYER+1)=0.0E+00
         B(L, 4*N_LAYER+1)=0.0E+00
!
         A(L, 11, 4*N_LAYER+2)=0.0E+00
         A(L, 10, 4*N_LAYER+2)=0.0E+00
         A(L, 9, 4*N_LAYER+2)=0.0E+00
         A(L, 8, 4*N_LAYER+2)=0.0E+00
         A(L, 7, 4*N_LAYER+2)=0.0E+00
         A(L, 6, 4*N_LAYER+2)=1.0E+00
         A(L, 5, 4*N_LAYER+2)=0.0E+00
         A(L, 4, 4*N_LAYER+2)=-ALBEDO_SURFACE_DIFF(L)
         A(L, 3, 4*N_LAYER+2)=0.0E+00
         A(L, 2, 4*N_LAYER+2)=0.0E+00
         A(L, 1, 4*N_LAYER+2)=0.0E+00
         B(L, 4*N_LAYER+2)=SOURCE_GROUND(L)
     &      +(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))
     &      *FLUX_DIRECT_GROUND(L)
      ENDDO
!
!
!
      RETURN
      END
