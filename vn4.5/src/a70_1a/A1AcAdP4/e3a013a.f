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
!+ Function to calculate third expoential integral.
!
! Method:
!       For small arguments a power series is used. For larger
!       arguments a Pade approximant derived from the asymptotic
!       expansion is used.
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
      FUNCTION E3_ACC01(X)
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
      REAL
     &     X
!             POINT OF EVALUATION
     &   , E3_ACC01
!             NAME OF FUNCTION
!
!     LOCAL VARIABLES
      REAL
     &     EULER
!             EULER'S CONSTANT
!
      PARAMETER(EULER=0.5772156E+00)
!
!
      IF (X.LT.1.0E-06) THEN
         E3_ACC01=0.5E+00
      ELSE IF (X.LT.2.0E+00) THEN
         E3_ACC01=-0.5E+00*X*X*LOG(X)+0.5E+00
     &      +X*(-1.0E+00+X*(0.75E+00-0.5E+00*EULER
     &      +X*(1.0E+00/6.0E+00-X*(1.0E+00/48.0E+00
     &      -X*(1.0E+00/3.60E+02-X*(1.0E+00/2.880E+03
     &      -X*(1.0E+00/2.5200E+04)))))))
      ELSE
!        WE USE A DOUBLY CUBIC PADE APPROXIMANT DERIVED FROM THE
!        ASYMPTOTIC EXPRESSION.
         E3_ACC01=(EXP(-X)/X)
     &      *(6.0E+00+X*(48.0E+00+X*(15.0E+00+X)))
     &      /(120.0E+00+X*(90.0E+00+X*(18.0E+00+X)))
      ENDIF
!
!
      RETURN
      END
