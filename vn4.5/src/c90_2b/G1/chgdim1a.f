C ******************************COPYRIGHT******************************
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
!+ Change 1st dimension of 2 dimensional array.
!
! Subroutine Interface:
      SUBROUTINE CHANGE_DIMENS(X,INSIZE,OUTSIZE,LEVELS,ICODE)

      IMPLICIT NONE
!
! Description and Method:
! Convert array a(insize,levels) with elements a((i=1,outsize),levels)
!  defined to an array of contiguous elements such that it is
!  equivalent to an array of dimension (outsize,levels).
!  Note that outsize must be le insize.
!
! Current Code Owner: Rick Rawlins (FR)
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  15/07/94  Original code. RR; implemented by RTHBarnes.
!  4.5  27/04/98  Add Fujitsu vectorization directive.
!                                    RBarnes@ecmwf.int
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER
     &       INSIZE           ! IN Input first dimension
     &      ,OUTSIZE          ! IN Output first dimension
     &      ,LEVELS           ! IN Input second dimension
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
      INTEGER  ICODE          ! INOUT Return code
!   Array  arguments with intent(InOut):
      REAL   X(INSIZE*LEVELS) ! INOUT Array for redimensioning
!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
      INTEGER
     &       I,LEVEL,I1,I2  ! Local loops and counters
! Local dynamic arrays:

! Function & Subroutine calls:
!     External - NONE

      IF (OUTSIZE.GT.INSIZE) THEN
        write(6,*) 'CHANGE_DIMENS: ERROR, OUTSIZE GT INSIZE'
        ICODE = 1
        GO TO 9999
      END IF
      DO LEVEL = 1,LEVELS
        I1 = (LEVEL-1)*INSIZE
        I2 = (LEVEL-1)*OUTSIZE
! Fujitsu vectorization directive
!OCL NOVREC
        DO  I = 1,OUTSIZE
          X(I+I2) = X(I+I1)
        END DO
      END DO
 9999 CONTINUE
      RETURN
      END
