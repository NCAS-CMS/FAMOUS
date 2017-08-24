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
!+ Subroutine to solve a set of banded matrix equations.
!
! Method:
!       A set of bands matrix equations is solved using the
!       standard method of Gaussian elimination.
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
      SUBROUTINE BAND_SOLVER(M, N
     &   , IU, IL
     &   , A, B
     &   , X
     &   , PM, PN
     &   , RHO
     &   )
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
      INTEGER   !, INTENT(IN)
     &     M
!             NUMBER OF PROFILES
     &   , N
!             NUMBER OF EQUATIONS
     &   , IU
!             NUMBER OF SUPERDIAGONALS
     &   , IL
!             NUMBER OF SUBDIAGONALS
     &   , PM
!             FIRST ARRAY SIZE
     &   , PN
!             SECOND ARRAY SIZE
      REAL      !, INTENT(IN)
     &     A(PM, (1+IU+IL), PN)
!             MATRIX OF COEFFICIENTS
     &   , B(PM, PN)
!             RIGHTHAND SIDES
      REAL      !, INTENT(OUT)
     &     X(PM, PN)
!             SOLUTION VECTOR
      REAL      !, INTENT(WORK)
     &     RHO(PM)
!             TEMPORARY ARRAY
!
!     LOCAL VARIABLES
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , IU1
!             LOCAL SCALAR
!
!
      IU1=IU+1
!     ELIMINATIVE PHASE.
      DO I=N, 2, -1
         DO J=1, MIN(IU, I-1)
            DO L=1, M
               RHO(L)=A(L, IU1-J, I-J)/A(L, IU1, I)
               B(L, I-J)=B(L, I-J)-RHO(L)*B(L, I)
            ENDDO
            DO K=1, MIN(IL, I-1)
               DO L=1, M
                  A(L, IU1+K-J, I-J)=A(L, IU1+K-J, I-J)
     &               -RHO(L)*A(L, IU1+K, I)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     SOLUTION AND BACK-SUBSTITUTION:
!
      IF ( (IU.EQ.2).AND.(IL.EQ.2) ) THEN
!        A SPECIAL VERSION IS USED FOR THE PENTADIAGONAL CASE TO ALLOW
!        US TO CHAIN OPERATIONS TOGETHER FOR EFFICIENCY ON THE CRAY.
!        FIRST EQUATION:
         DO L=1, M
            X(L, 1)=B(L, 1)/A(L, 3, 1)
         ENDDO
!        SECOND EQUATION:
         DO L=1, M
            X(L, 2)=(B(L, 2)-A(L, 4, 2)*X(L, 1))/A(L, 3, 2)
         ENDDO
!        REMAINING EQUATIONS:
         DO I=3, N
            DO L=1, M
               X(L, I)=(B(L, I)-A(L, 4, I)*X(L, I-1)
     &            -A(L, 5, I)*X(L, I-2))/A(L, 3, I)
            ENDDO
         ENDDO
      ELSE
!        GENERAL CASE:
         DO I=1, N
            DO L=1, M
               X(L, I)=B(L, I)
            ENDDO
            DO K=1, MIN(IL, I-1)
               DO L=1, M
                  X(L, I)=X(L, I)-A(L, IU1+K, I)*X(L, I-K)
               ENDDO
            ENDDO
            DO L=1, M
               X(L, I)=X(L, I)/A(L, IU1, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
