!+ Outputs error messages
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
! Subroutine Interface:
      SUBROUTINE ErrorReport (NoOfLines, NameOfRoutine, Message,
     &  ErrorStatus)

      IMPLICIT NONE
!
! Description:
!   Displays error messages. Severity of error is determined from the
!  value of ErrorStatus. ErrorStatus is reset to ErrorStatus_OK after
!  reporting a warning message.
!
!  ErrorStatus = ErrorStatus_OK   No error
!  ErrorStatus < ErrorStatus_OK   Warning message
!  ErrorStatus > ErrorStatus_OK   Fatal Error message
!
!  ErrorStatus_OK is defined in the comdec C_ErrVal
!
! Method:
!
! Current Code Owner: Phil Andrews
!
! History:
! Version  Date    Comment
! -------  ----    -------
! 3.4      11/3/94 Original code. Phil Andrews
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
! 1.0 Global variables (*CALLed common blocks etc.):
!
! Description:
!   Sets parameters with possible ErrorStatus values
!
! Current Code Owner: Phil Andrews
!
! History:
! Date    Version     Comment
! ----    -------     -------
!24/12/93   3.4       Original code. Phil Andrews
!
! Declarations:
! Global parameters:
      INTEGER      ErrorStatus_OK         ! Good value of ErrorStatus
        PARAMETER (ErrorStatus_OK = 0)

!- End of COMDECK declaration

! 2.0 Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      NoOfLines            ! No of output lines

      CHARACTER*80 NameOfRoutine        ! calling this one

!   Array  arguments with intent(in):
      CHARACTER*80 message(NoOfLines)   ! Message to output

!   ErrorStatus
      INTEGER      ErrorStatus          ! Error flag (0 = no error)

! 3.0 Local scalars:
      INTEGER      line                 ! loop counter

!- End of header

      WRITE(6,*) ' '
      WRITE(6,*) 'ErrorReport:'

      If (ErrorStatus .ne. ErrorStatus_OK) then
        If (ErrorStatus .lt. ErrorStatus_OK) then
          WRITE(6,*) '!! WARNING from '//NameOfRoutine

        Else ! ErrorStatus .gt. ErrorStatus_OK
          WRITE(6,*) '!! FATAL ERROR in '//NameOfRoutine

        End if

        Do line = 1, NoOfLines
          WRITE(6,*) message(line)

        End do
      End if

      Return
      End
