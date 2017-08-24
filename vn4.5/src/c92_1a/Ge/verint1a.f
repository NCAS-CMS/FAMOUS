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
        Subroutine VerticalInterp
     &          ( Scheme,
     &        len1In,
     &          len2In,
     &          len3In,
     &          len3Out,
     &          RIn,
     &          ROut,
     &          DataIn,
     &          DataOut,
     &          ErrorStatus,ErrorMessage)

!   Description
!   Performs linear, cubic Lagrange or quintic lagrange vertical interp
!   on 3D fields. The horizontal location of points must be the same for
!   both the input and the output data arrays.
!   N.B. This routine should NOT be used for vertically interpolating
!        density (or pressure?).
!
! Method:
!   This code is a simplifed version of the full 3D interpolation
!   routine written by Mark Mawson
!                The proposed semi-Lagrangian advection scheme for the
!                   semi-Implicit Unified Model integration scheme.
!                         F.R. Division working paper No 162.
!
! Current Code Owner: Stuart Bell
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      11/7/95   Equiv. to VAR code as at time of build:Stuart Bell
! 4.1      1/2/96    Duff Spelling of subroutine name. Stuart Bell
!
! Code Description:
!   Language:           Fortran 77 plus
!   Software Standards: "UM and Met O standards".
!
! Declarations:
!

        Implicit None

! Subroutine arguments
! Scalar arguments with : intent(in)
        Integer   Scheme     ! a code saying which order
!                            ! of scheme to use:
!                            ! 1 = linear;3 = cubic;5 = quintic
        Integer   len1In    ! Extent of DataIn /out in i direction.
        Integer   len2In    ! Extent of DataIn /out in j direction.
        Integer   len3In    ! Extent of DataIn in k direction.
        Integer   len3Out   ! Extent of DataOut in k direction.

! Array  arguments with :intent(in)
        Real  RIn  (len1In,len2In,len3In)  ! Vertical coordinate input

        Real  ROut (len1In,len2In,len3Out)  ! Vertical coordinate output

        Real  DataIn (len1In,len2In,len3In) ! Data on the original
!                                           ! vertical levels.

! Array  arguments with : intent(out)
        Real     DataOut (len1In,len2In,len3Out) ! data interpolated to
!                                            ! new vertical levels.
!  ErrorStatus
      INTEGER       ErrorStatus
      CHARACTER*256 ErrorMessage

! Local scalars:
        Integer              i, j, k, index  ! Loop indices
        Integer              VIorder         ! order of interpolation
        Integer              ErrorStat       ! Error status

! Local arrays:
        Integer   KOut (len1In,len2In,len3In) !level in RIn below ROut

!- End of header

! ----------------------------------------------------------------------
!  Section 1.   Initialize
! ----------------------------------------------------------------------
      ErrorStatus = 0
      ErrorMessage= "  "

! Check that a valid interpolation scheme has been specified.
        If (( Scheme .ne. 1 .and.
     &      Scheme .ne. 3 .and.
     &      Scheme .ne. 5       ) ) then

       VIorder = 1

         ErrorMessage="Invalid value of 'Scheme':linear interp. done  "
         ErrorStat     = -1

        Else
       VIorder = Scheme

        End If

! ----------------------------------------------------------------------
! Section 2.   For each output point find k so that the point on the
!              output grid lies between k and k+1
! ----------------------------------------------------------------------

! Find k point.
! Set minimum value to level one.
        Do k = 1, len3Out
         Do j = 1, len2In
          Do i = 1, len1In
         KOut(i,j,k) = 1
        End Do
       End Do
      End Do

! Find level which is just below ROut value
        Do index = 2, len3In - 1
         Do k = 1, len3Out
          Do j = 1, len2In
            Do i = 1, len1In
           If ( ROut(i,j,k) .gt. RIn(i,j,index) ) KOut(i,j,k) = index

          End Do
         End Do
        End Do
       End Do

! ----------------------------------------------------------------------
! Section 3.   Perform required Interpolations.
! ----------------------------------------------------------------------

! Call the specified interpolation scheme:
         If(VIorder.eq.1)Then
             Call VerticalInterpLinear
     &      ( len1In,
     &        len2In,
     &        len3In,
     &        len3Out,
     &        KOut,
     &        RIn,
     &        ROut,
     &        DataIn,
     &        DataOut )

         ElseIf(VIorder.eq.3)Then
             Call VerticalInterpCubic
     &      ( len1In,
     &        len2In,
     &        len3In,
     &        len3Out,
     &        KOut,
     &        RIn,
     &        ROut,
     &        DataIn,
     &        DataOut )

         ElseIf(VIorder.eq.5)Then
             Call VerticalInterpQuintic
     &      ( len1In,
     &        len2In,
     &        len3In,
     &        len3Out,
     &        KOut,
     &        RIn,
     &        ROut,
     &        DataIn,
     &        DataOut )

         EndIf


         Return
         End
