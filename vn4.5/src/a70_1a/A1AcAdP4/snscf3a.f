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
!+ Function to set number of source coefficients.
!
! Method:
!       The two-stream approximation is examined and the number
!       of coefficients is set accordingly.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET FLAGS FOR DIFFERENT PORTIONS
!     OF THE SPECTRUM.
!
      INTEGER
     &     IP_SOLAR
!             SOLAR REGION
     &   , IP_INFRA_RED
!             INFRA-RED REGION
!
      PARAMETER(
     &     IP_SOLAR=1
     &   , IP_INFRA_RED=2
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             SPECTRAL REGION
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             FLAG FOR QUADRATIC INFRA-RED SOURCE
!
      INTEGER   !, INTENT(OUT)
     &     SET_N_SOURCE_COEFF
!             RETURNED NUMBER OF SOURCE COEFFICIENTS
!
!
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         SET_N_SOURCE_COEFF=2
      ELSE
         IF (L_IR_SOURCE_QUAD) THEN
            SET_N_SOURCE_COEFF=2
         ELSE
            SET_N_SOURCE_COEFF=1
         ENDIF
      ENDIF
!
!
!
      RETURN
      END
