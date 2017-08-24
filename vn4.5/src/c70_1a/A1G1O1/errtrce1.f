!+ Provides calling tree for errors.
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
!
! Subroutine Interface:
      SUBROUTINE ErrorTrace (NameOfRoutine, message)

      IMPLICIT NONE
!
! Description:
!   Outputs NameOfRoutine and message. Should be called if ErrorStatus
!   .ne. ErrorStatus_OK just before exiting a routine.
!
! Method:
!   Straight forward.
!
! Current Code Owner: Phil Andrews
!
! History:
! Version  Date      Comment
! -------  ----      -------
! 3.4      11/3/94   Original code. Phil Andrews
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!
! Global variables (*CALLed common blocks etc.):

! Subroutine arguments
!   Scalar arguments with intent(in):
      CHARACTER*80   NameOfRoutine      ! that called this one.
      CHARACTER*80   message            ! to output

!- End of header

      WRITE(6,*) ' '
      WRITE(6,*) 'ErrorTrace called by: '//NameOfRoutine
      WRITE(6,*) message

      Return
      End
!+
