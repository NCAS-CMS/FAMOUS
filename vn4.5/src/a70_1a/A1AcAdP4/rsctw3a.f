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
!+ Subroutine to rescale optical depth and albedo.
!
! Method:
!       The standard rescaling formulae are applied.
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
      SUBROUTINE RESCALE_TAU_OMEGA(N_PROFILE
     &   , I_LAYER_FIRST, I_LAYER_LAST
     &   , TAU, OMEGA, FORWARD_SCATTER
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST
!             FIRST LAYER TO RESCALE
     &   , I_LAYER_LAST
!             FIRST LAYER TO RESCALE
      REAL      !, INTENT(IN)
     &     FORWARD_SCATTER(NPD_PROFILE, NPD_LAYER)
!             FORWARD SCATTERING
      REAL      !, INTENT(INOUT)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            TAU(L, I)=TAU(L, I)*(1.0E+00
     &         -OMEGA(L, I)*FORWARD_SCATTER(L, I))
            OMEGA(L, I)=OMEGA(L, I)*(1.0E+00-FORWARD_SCATTER(L, I))
     &         /(1.0E+00-OMEGA(L, I)*FORWARD_SCATTER(L, I))
         ENDDO
      ENDDO
!
!
      RETURN
      END
